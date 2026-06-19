#include <iostream>
#include <chrono>
#include <random>
#include <string>
#include <vector>
#include "SeedGenerator/Container.h"
#include "SeedGenerator/RadiusCalculator.h"
#include "ScaffoldGenerator/GeneratorLewiner.h"
#include "Logger/Logger.h"

// Standalone CLI harness (no GUI/OpenGL) to measure where time goes in the
// scaffold generation pipeline using a real mesh container, matching the
// actual GUI workflow settings as closely as possible.
//
// Usage: scaffoldProfile [meshPath|--box] [scale] [resX resY resZ] [--no-thickness]
//   meshPath        path to an STL container, or "--box" for a trivial O(1)
//                    BoxContainer (isolates the seed-kdtree cost)
//   scale            uniform scale applied to the mesh about its bbox center
//                    (default 100, e.g. to bring a cm-ish STL into mm)
//   resX resY resZ   explicit grid resolution; if omitted, derived from the
//                    (scaled) mesh bounds at a 0.05 (mm) voxel size
//   --no-thickness   disable the "distance from container" varied-thickness
//                    function, to isolate its added cost
int main(int argc, char** argv) {

	std::vector<std::string> args(argv + 1, argv + argc);

	bool noThickness = false;
	bool reproBug = false;
	std::vector<std::string> positional;
	for (auto& a : args) {
		if (a == "--no-thickness") noThickness = true;
		else if (a == "--repro-bug") reproBug = true;
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
		float voxelTarget = 0.05f;
		res = {
			std::max(2, (int)((b.xMax - b.xMin) / voxelTarget)),
			std::max(2, (int)((b.yMax - b.yMin) / voxelTarget)),
			std::max(2, (int)((b.zMax - b.zMin) / voxelTarget))
		};
	}
	size_t scaledVoxels = (size_t)res[0] * res[1] * res[2];
	std::cout << "Profiling resolution: "
		<< res[0] << "x" << res[1] << "x" << res[2]
		<< " = " << scaledVoxels << " voxels" << std::endl;

	// generate representative seed scatter (porous lattice seed count)
	std::mt19937 rng(42);
	std::uniform_real_distribution<float> dx(b.xMin, b.xMax), dy(b.yMin, b.yMax), dz(b.zMin, b.zMax);
	std::vector<Vec3> seeds;
	int seedNr = 2000;
	seeds.reserve(seedNr);
	for (int i = 0; i < seedNr; i++) {
		seeds.push_back(Vec3(dx(rng), dy(rng), dz(rng)));
	}
	std::cout << "Seed count: " << seedNr << std::endl;

	float isoLevel = 0.255f;
	GeneratorLewiner gen(seeds, bounds, res, &logger, 0.5f, isoLevel, 0, false);

	if (reproBug) {
		// reproduces the reported "millions of floating pieces" scenario:
		// startThickness=5, endThickness=0.3, transitionDistance=5 (distance
		// from container), stretchY=2.0
		gen.set_stretch(1.0f, 2.0f, 1.0f);
		std::cout << "Varied thickness (repro): start=5 end=0.3 transitionDist=5 stretchY=2.0" << std::endl;
		gen.set_thickness_functions(
			container->sdf,
			std::make_shared<LinearFunction>(5.0f),
			5.0f, 0.3f, 5.0f
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