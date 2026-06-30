#include <iostream>
#include <chrono>
#include <random>
#include <string>
#include <vector>
#include "SeedGenerator/Container.h"
#include "SeedGenerator/SeedGenerator.h"
#include "SeedGenerator/RadiusCalculator.h"
#include "ScaffoldGenerator/GeneratorLewiner.h"
#include "Logger/Logger.h"

// Standalone CLI harness (no GUI/OpenGL) to measure where time goes in the
// scaffold generation pipeline using a real mesh container, matching the
// actual GUI workflow settings as closely as possible.
//
// Usage: scaffoldProfile [meshPath|--box] [scale] [resX resY resZ] [flags]
//   meshPath        path to an STL container, or "--box" for a trivial O(1)
//                    BoxContainer (isolates the seed-kdtree cost)
//   scale            uniform scale applied to the mesh about its bbox center
//                    (default 100, e.g. to bring a cm-ish STL into mm)
//   resX resY resZ   explicit grid resolution; overrides --voxel-size
//   --no-thickness   disable the "distance from container" varied-thickness
//                    function, to isolate its added cost
//   --repro-bug      use Poisson-uniform seeding + varied thickness settings
//                    for reproducing reported generation defects
//   --voxel-size v   target voxel size in mm (default 0.05); derives
//                    resolution from the (scaled) mesh bounds
//   --start-thickness v / --end-thickness v / --transition v / --stretch-y v
//                    override the --repro-bug thickness/anisotropy defaults
//   --poisson-radius v
//                    Poisson-uniform seed spacing (default 0.583)
int main(int argc, char** argv) {

	std::vector<std::string> args(argv + 1, argv + argc);

	bool noThickness = false;
	bool reproBug = false;
	float voxelSizeArg = 0.05f;
	float startThicknessArg = 5.0f;
	float endThicknessArg = 0.3f;
	float transitionArg = 5.0f;
	float stretchYArg = 2.0f;
	float poissonRadius = 0.583f;
	std::vector<std::string> positional;

	for (size_t i = 0; i < args.size(); i++) {
		const std::string& a = args[i];
		if (a == "--no-thickness") noThickness = true;
		else if (a == "--repro-bug") reproBug = true;
		else if (a == "--voxel-size" && i + 1 < args.size()) voxelSizeArg = std::stof(args[++i]);
		else if (a == "--start-thickness" && i + 1 < args.size()) startThicknessArg = std::stof(args[++i]);
		else if (a == "--end-thickness" && i + 1 < args.size()) endThicknessArg = std::stof(args[++i]);
		else if (a == "--transition" && i + 1 < args.size()) transitionArg = std::stof(args[++i]);
		else if (a == "--stretch-y" && i + 1 < args.size()) stretchYArg = std::stof(args[++i]);
		else if (a == "--poisson-radius" && i + 1 < args.size()) poissonRadius = std::stof(args[++i]);
		else positional.push_back(a);
	}

	std::string meshPath = positional.size() > 0 ? positional[0] : "data/femoral_head_container.stl";
	bool useBox = (meshPath == "--box");
	float scale = positional.size() > 1 ? std::stof(positional[1]) : 100.0f;

	Logger& logger = Logger::get_instance();

	IContainer* container = nullptr;
	if (useBox) {
		std::cout << "Using trivial BoxContainer (isolates seed-kdtree cost)" << std::endl;
		// same bounds as data/femoral_head_container.stl scaled by 100, hardcoded
		// so this works without loading any mesh
		container = new BoxContainer(Vec3(45.14f, 49.11f, 50.92f), Vec3(0.0f, -2.67f, 3.40f), false);
	}
	else {
		std::cout << "Loading mesh: " << meshPath << std::endl;
		AbstractContainer* mesh = new AbstractContainer(meshPath, false);
		IContainer* meshAsContainer = mesh;
		Bounds rawBounds = meshAsContainer->compute_bounds();
		std::cout << "Raw mesh bounds: [" << rawBounds.xMin << ", " << rawBounds.xMax << "]" << std::endl;
		std::cout << "Applying scale x" << scale << std::endl;
		mesh->set_scale(scale);
		container = mesh;
	}

	Bounds b = container->compute_bounds();
	std::cout << "Scaled bounds: ["
		<< b.xMin << ", " << b.xMax << "] x ["
		<< b.yMin << ", " << b.yMax << "] x ["
		<< b.zMin << ", " << b.zMax << "]" << std::endl;

	std::array<float, 6> bounds = {
		(float)b.xMin, (float)b.xMax, (float)b.yMin, (float)b.yMax, (float)b.zMin, (float)b.zMax
	};

	std::array<int, 3> res;
	if (positional.size() > 4) {
		res = { std::stoi(positional[2]), std::stoi(positional[3]), std::stoi(positional[4]) };
	}
	else {
		res = {
			std::max(2, (int)((b.xMax - b.xMin) / voxelSizeArg)),
			std::max(2, (int)((b.yMax - b.yMin) / voxelSizeArg)),
			std::max(2, (int)((b.zMax - b.zMin) / voxelSizeArg))
		};
	}
	size_t scaledVoxels = (size_t)res[0] * res[1] * res[2];
	std::cout << "Profiling resolution (voxel size " << voxelSizeArg << "): "
		<< res[0] << "x" << res[1] << "x" << res[2]
		<< " = " << scaledVoxels << " voxels" << std::endl;

	std::vector<Vec3> seeds;
	if (reproBug) {
		// real Poisson-disk uniform sampling (rMin==rMax), matching the
		// GUI's "Poisson Uniform" seed generator - NOT a naive random
		// scatter, which has no minimum-spacing guarantee and can produce
		// very different local cell topology
		Poisson3D poisson(poissonRadius, poissonRadius, 30, false);
		poisson.run(*container);
		seeds = poisson.get_seeds();
		std::cout << "Poisson-uniform seed count (radius=" << poissonRadius << "): " << seeds.size() << std::endl;
	}
	else {
		// representative seed scatter for non-repro profiling runs
		std::mt19937 rng(42);
		std::uniform_real_distribution<float> dx(b.xMin, b.xMax), dy(b.yMin, b.yMax), dz(b.zMin, b.zMax);
		int seedNr = 2000;
		seeds.reserve(seedNr);
		for (int i = 0; i < seedNr; i++) {
			seeds.push_back(Vec3(dx(rng), dy(rng), dz(rng)));
		}
		std::cout << "Seed count: " << seedNr << std::endl;
	}

	float isoLevel = 0.255f;
	GeneratorLewiner gen(seeds, bounds, res, &logger, 0.5f, isoLevel, 0, false);

	if (reproBug) {
		gen.set_stretch(1.0f, stretchYArg, 1.0f);
		std::cout << "Varied thickness (repro): start=" << startThicknessArg
			<< " end=" << endThicknessArg << " transitionDist=" << transitionArg
			<< " stretchY=" << stretchYArg << std::endl;
		gen.set_thickness_functions(
			container->sdf,
			std::make_shared<LinearFunction>(transitionArg),
			startThicknessArg, endThicknessArg, transitionArg
		);
	}
	else if (!noThickness && !useBox) {
		// matches the GUI's "varied thickness, distance from container" mode:
		// thicknessSDF is literally the container's own SDF (con.sdf)
		float startThickness = 10.0f;     // at the container surface (solid)
		float endThickness = isoLevel;    // deep interior (porous)
		float transitionDistance = 2.0f;  // mm from the surface
		std::cout << "Varied thickness ENABLED: start=" << startThickness
			<< " end=" << endThickness << " transitionDist=" << transitionDistance << std::endl;
		gen.set_thickness_functions(
			container->sdf,
			std::make_shared<LinearFunction>(transitionDistance),
			startThickness, endThickness, transitionDistance
		);
	}
	else {
		std::cout << "Varied thickness DISABLED" << std::endl;
	}

	std::cout << "\n--- compute_scalar_field ---" << std::endl;
	auto t0 = std::chrono::steady_clock::now();
	gen.compute_scalar_field(*container);
	auto t1 = std::chrono::steady_clock::now();

	std::cout << "\n--- marching_cubes ---" << std::endl;
	gen.marching_cubes();
	auto t2 = std::chrono::steady_clock::now();

	std::cout << "\n=== SUMMARY ===" << std::endl;
	std::cout << "compute_scalar_field TOTAL: " << std::chrono::duration<double>(t1 - t0).count() << " s" << std::endl;
	std::cout << "marching_cubes TOTAL:       " << std::chrono::duration<double>(t2 - t1).count() << " s" << std::endl;

	delete container;
	return 0;
}