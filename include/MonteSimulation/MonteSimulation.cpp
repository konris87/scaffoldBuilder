//#define EIGEN_DONT_PARALLELIZE      // Prevents thread explosion with OpenMP
//#define EIGEN_MAX_ALIGN_BYTES 0

#include "SeedGenerator/Container.h"
#include "SeedGenerator/SeedGenerator.h"
//#include "ScaffoldGenerator/ScaffoldGenerator.h"
//#include "ScaffoldGenerator/ScaffoldGeneratorBox.h"
#include "MonteSimulation.h"
#include "voro++.hh"
#include <omp.h>
#include <limits>
#include "Eigen/Dense"

LhsSampler::LhsSampler(int nr, int dims) : sampleNr{nr}, dimNr(dims) {};

std::vector<std::vector<double>> LhsSampler::generate_parameters() {
	
	std::vector<std::vector<double >> table(sampleNr, std::vector<double>(dimNr));

	// this is the vector containing the index of simulation set i
	std::vector<int> index(sampleNr);

	// create a random engine
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0.0, 1.0);

	// create columns for each dimension
	for (int dim{ 0 }; dim < dimNr; dim++) {

		// create indices
		std::iota(index.begin(), index.end(), 0);

		// randomly shuffle them to randomize rows
		std::shuffle(index.begin(), index.end(), gen);

        // generate samples
        for (int j = 0; j < sampleNr; ++j) {
            // Formula: (BinIndex + RandomOffset) / NumBins
            // This ensures we pick exactly one value from the i-th slice of the histogram
            double binVal = index[j];
            double noise = dis(gen);

            table[j][dim] = (binVal + noise) / static_cast<double>(sampleNr);
        }
	}
	return table;
};

//-----------------------------------------------------------
// generate lhs sampling for the scaffold
std::vector<SimInput> generate_lhs_sampling(int nr, int dims, SamplingType t) {

	// ranges
	const double vol = 1000; // we keep creating everything inside a 10x10x10 cube
	const int seedMin = 50;
	const int seedMax = 2000;
	const double rMin = 0.6; // pore separation mm
	const double rMax = 3.0;  // pore separation mm
	//const double tMin = 0.1; // minimum thickness mm
	//const double tMax = 3.0; // maximum thickness mm

	const double ratioMin = 0.1;
	const double ratioMax = 0.30;

	const double pMin = 0.5;
	const double pMax = 1.0;

	std::vector<SimInput> inputs(nr);

	dims = (t == SamplingType::varied) ? 4 : 3;

	// first generate the table using Lhs 

	LhsSampler sampler(nr, dims);

	std::vector<std::vector<double>> table = sampler.generate_parameters();

	// apply the limits
	for (int i = 0; i < nr; ++i) {
		inputs[i].id = i;
		inputs[i].type = t;

		double characteristicScale = 0.0;

		switch (t) {
		case random:{
			double count = seedMin + (seedMax - seedMin) * table[i][0];
			inputs[i].p1 = (int)count; // Density
			inputs[i].p2 = 0.0;         // Unused

			characteristicScale = std::pow(vol / count, 1.0 / 3.0);
			break;
		}
		case uniform: {
			double radius = rMin + (rMax - rMin) * table[i][0];
			inputs[i].p1 = radius;
			inputs[i].p2 = 0.0; // Unused

			characteristicScale = radius;
			break;
		}
		case varied: {
				double valA = rMin + (rMax - rMin) * table[i][0];
				double valB = rMin + (rMax - rMin) * table[i][1];
				// always rmin should be lower than rmax
				inputs[i].p1 = std::min(valA, valB); // Min Radius
				inputs[i].p2 = std::max(valA, valB); // Max Radius

				// set characteristic scale to the average
				characteristicScale = (inputs[i].p1 + inputs[i].p2) * 0.5;
				break;
			}
		}

		// thickness idx
		int tIdx = (t == SamplingType::varied) ? 2 : 1;
		// pulling back idx
		int pIdx = (t == SamplingType::varied) ? 3 : 2;

		// we use a thickness as a scale of the radius
		double tau = ratioMin + (ratioMax - ratioMin) * table[i][tIdx];

		double thickness = characteristicScale * tau;
		inputs[i].thickness = std::clamp(thickness, 0.1, 5.0);
		inputs[i].pullbackRatio = pMin + (pMax - pMin) * table[i][pIdx];
	}

	return inputs;
};


// monte carlo simulation
void run_monte_carlo_simulations(int nr, int dims, SamplingType t) {

	//// get the inputs
	//std::vector<SimInput> inputs = generate_lhs_sampling(nr, dims, t);

	//// create a vector for the outputs
	//std::vector<SimOutput> outputs(nr);

	//std::cout << "Starting " << nr << " simulations..." << std::endl;

	//#pragma omp parallel for schedule(dynamic)
	//for (int i = 0; i < nr; i++) {
	//	try {

	//		BoxContainer boxContainer;
	//		IContainer* activeContainer = &boxContainer;
	//		// Re-initialize adapter with the fresh container dimensions (0-10)
	//		ContainerAdapter adapter = { *activeContainer, 10.0, 10.0, 10.0 };
	//		vtkSmartPointer<vtkPolyData> scaffoldPoly = vtkSmartPointer<vtkPolyData>::New();

	//		std::vector<std::array<double, 3>> seeds;
	//		std::unique_ptr<Graph> graph;
	//		// ---------------------------------------------

	//		const auto& params = inputs[i];
	//		
	//		//std:cout << params.id << " " << params.p1 << " " <<
	//		//params.p2 << " " << params.pullbackRatio << " " << params.thickness << std::endl;
	//		
	//		SimOutput out;
	//		out.id = params.id;

	//		// Generate points
	//		switch (t) {
	//		case random: {
	//			int nrSeeds = static_cast<int>(params.p1);
	//			Random rnd(nrSeeds);
	//			rnd.run(adapter, seeds);
	//			break;
	//		}
	//		case uniform: {
	//			double r = params.p1;
	//			Poisson3D gnr(r, r, { 5.0f, 5.0f, 5.0f }, 30);
	//			gnr.run(adapter, seeds);
	//			break;
	//		}
	//		case varied: {
	//			double rMin = params.p1;
	//			double rMax = params.p2;
	//			Poisson3D gnr(rMin, rMax, { 5.0f, 5.0f, 5.0f }, 30);

	//			// Gradient along X
	//			std::array<double, 3> center = { 0.0, 5.0, 5.0 };
	//			std::array<double, 3> normal = { 1.0, 0.0, 0.0 };

	//			// Need a local config object
	//			RunConfig cfg;
	//			cfg.dist = std::make_shared<PlaneDistEstimator>(center, normal);

	//			gnr.run(adapter, seeds, cfg);
	//			break;
	//			}
	//		}

	//		Bounds bounds = activeContainer->compute_bounds();

	//		std::array<double, 6> bnds{
	//				static_cast<double>(bounds.xMin),
	//				static_cast<double>(bounds.xMax),
	//				static_cast<double>(bounds.yMin),
	//				static_cast<double>(bounds.yMax),
	//				static_cast<double>(bounds.zMin),
	//				static_cast<double>(bounds.zMax)
	//		};

	//		//std::unique_ptr<GeneratorInterface> generator;
	//		//
	//		//std::array<int, 3> modelDims = { 50, 50, 50 };

	//		//auto* box = dynamic_cast<const BoxContainer*>(activeContainer);
	//		////std::cout << box->xDim << " " << box->yDim << " " << box->zDim << std::endl;

	//		//auto sgb = std::make_unique<ScaffoldGeneratorBox>(
	//		//	seeds, bnds, modelDims,
	//		//	params.thickness,
	//		//	params.pullbackRatio, 
	//		//	box->is_inside()
	//		//);

	//		//generator = std::move(sgb);

	//		//// Run Voro Analysis
	//		//// Safety check: Don't run Voro++ on empty seeds (causes crash)
	//		//if (!seeds.empty()) {
	//		//	generator->run_process_faces();
	//		//	generator->run_generate_mesh(
	//		//		params.thickness, 
	//		//		scaffoldPoly,
	//		//		modelDims
	//		//	);

	//			vtkNew<vtkMassProperties> massProperties;
	//			massProperties->SetInputData(scaffoldPoly);
	//			massProperties->Update();

	//			double volume = massProperties->GetVolume();
	//			double domainVolume = 1000.0;

	//			out.porosity = 1 - (volume / domainVolume);
	//			//std::cout << "porosity: " << out.porosity << std::endl;

	//			out.connectivity = generator->run_estimate_connectivity();
	//			//std::cout << "connectivity: " << out.connectivity << std::endl;
	//		}
	//		else {
	//			out.porosity = 1.0;
	//			out.connectivity = 0.0;
	//		}

	//		outputs[i] = out;

	//		if (i % 100 == 0) std::cout << i << "\n" << std::flush;
	//	}
	//	catch (const std::exception& e) {
	//		std::cerr << "CRASH in Sim " << i << ": " << e.what() << std::endl;
	//	}
	//	catch (...) {
	//		std::cerr << "UNKNOWN CRASH in Sim " << i << std::endl;
	//	}
	//}
	//std::cout << std::endl;

	//std::string fileName;
	//// save to csv
	//switch (t) {
	//case random: {
	//	fileName = "random_sampling_results.csv";
	//	break;
	//	}
	//case uniform: {
	//	fileName = "uniform_sampling_results.csv";
	//	break;
	//}
	//case varied: {
	//	fileName = "varied_sampling_results.csv";
	//	break;
	//}
	//}
	//save_csv(fileName, inputs, outputs, t);
	//std::cout << "Done. Saved to " << fileName << std::endl;
};

void run_voro(
	std::vector<std::array<double, 3>>& seeds,
	std::unique_ptr<Graph>& graph,
	const SimInput& params,
	SimOutput& out
) {
	if (seeds.empty()) return;
	double solidVolume = 0.0;
	
	size_t seedNr = seeds.size();
	int nGrid = std::max(5, static_cast<int>(std::pow(seedNr / 8.0, 1.0 / 3.0)));

	voro::container con(
		0.0, 10.0, 0.0, 10.0, 0.0, 10.0,
		nGrid, nGrid, nGrid,
		false, false, false, 16
	);

	for (int i{ 0 }; i < seeds.size(); i++) {
		con.put(i, seeds[i][0], seeds[i][1], seeds[i][2]);
	}

	// populate the graph
	voro::c_loop_all cla(con);
	voro::voronoicell_neighbor cell;

	graph = std::make_unique<Graph>();

	// populate the graph nodes
	if (cla.start()) do {
		graph->add_vertex(cla.pid());
	} while (cla.inc());

	// loop through the cells
	if (cla.start()) do if (con.compute_cell(cell, cla)) {

		// get the id of the seed
		int seedId = cla.pid();

		// get the neighbors and the corresponding face area
		std::vector<int> cellNeighs;
		cell.neighbors(cellNeighs);

		std::vector<int> faceIndices;
		cell.face_vertices(faceIndices);

		double px{ 0.0 }, py{ 0.0 }, pz{ 0.0 };
		std::vector<double> cellVertices;
		cla.pos(px, py, pz);
		cell.vertices(px, py, pz, cellVertices);

		// get each face areas
		std::vector<double> faceAreas;
		cell.face_areas(faceAreas);

		// get each face perimeter
		std::vector<double> faceNormalsRaw;
		cell.normals(faceNormalsRaw);

		const int numFaces = static_cast<int>(cellNeighs.size());
		std::vector<double> faceMinWidths(numFaces, 0.0);
		//std::vector<double> faceMaxWidths(numFaces, 0.0);
		std::vector<Eigen::Vector3d> faceNormals(numFaces, Eigen::Vector3d::Zero());
		std::vector<double> planeCoeffs(numFaces, 0.0);
		std::vector<Eigen::Vector3d> faceCentroids(numFaces, Eigen::Vector3d::Zero());

		int idxOffset = 0;

		// --- Pass 1: Compute minWidth and normals for each face ---
		for (int i = 0; i < numFaces; ++i) {

			int neighborId = cellNeighs[i];
			int order = faceIndices[idxOffset];

			// avoid process of cells on the boundaries
			//if (neighborId >= 0) {

			std::vector<double> faceVertices3D;

			for (int j = 1; j <= order; ++j) {
				int localVertIdx = faceIndices[idxOffset + j];
				faceVertices3D.push_back(cellVertices[3 * localVertIdx + 0]);
				faceVertices3D.push_back(cellVertices[3 * localVertIdx + 1]);
				faceVertices3D.push_back(cellVertices[3 * localVertIdx + 2]);
			}

			if (3 * i + 2 < static_cast<int>(faceNormalsRaw.size())) {
				faceNormals[i] = Eigen::Vector3d(
					faceNormalsRaw[3 * i + 0],
					faceNormalsRaw[3 * i + 1],
					faceNormalsRaw[3 * i + 2]
				).normalized();
			}
			else {
				faceNormals[i] = Eigen::Vector3d::Zero();
			}

			// add the to connected vertices also the mean (centroid) of the face only once
			Eigen::Vector3d centroid = vector_to_matrix3D(faceVertices3D).colwise().mean().transpose();
			faceCentroids[i] = centroid;

			// estimate the plane coefficient d that sits on the centroid
			planeCoeffs[i] = faceNormals[i].dot(centroid);
			idxOffset += order + 1;
		}

		// --- Pass 2: For each face check if the eroded cell with block it or not
		idxOffset = 0;
		for (int i = 0; i < numFaces; ++i) {

			int neighborId = cellNeighs[i];
			int order = faceIndices[idxOffset];

			//if (neighborId >= 0) {

			Eigen::Vector3d normal = faceNormals[i];

			std::vector<double> faceVertices3D;

			for (int j = 1; j <= order; ++j) {
				int localVertIdx = faceIndices[idxOffset + j];
				faceVertices3D.push_back(cellVertices[3 * localVertIdx + 0]);
				faceVertices3D.push_back(cellVertices[3 * localVertIdx + 1]);
				faceVertices3D.push_back(cellVertices[3 * localVertIdx + 2]);
			}

			Eigen::MatrixXd verts2D; Eigen::Vector3d u, v, origin;
			project_vertices_on_plane(faceVertices3D, origin, u, v, verts2D);

			Eigen::Matrix<double, 3, 2> basis;
			basis.col(0) = u;
			basis.col(1) = v;

			double pullbackRatio = params.pullbackRatio;

			// estimate the delta as the inlet radius
			float inradius = polygon_inradius(verts2D);
			float delta = pullbackRatio * inradius;
			float dstar = delta;

			/*double dstar = erosion_margin_for_face(
				i,
				faceNormals,
				planeCoeffs,
				basis,
				origin,
				verts2D,
				delta
			);*/

			double thickness = params.thickness;

			bool isOpen = (dstar >= thickness * 0.5);
			
			if (isOpen && neighborId >= 0) {

				graph->add_edge(seedId, neighborId,
					{ dstar, delta });

				graph->add_edge(neighborId, seedId, { dstar, delta });
			}

			// avoid process of duplicates
			if (neighborId < 0 || seedId < neighborId) {
				double surfaceArea = faceAreas[i];
				double holeArea = 0.0;

				if (isOpen) {

					double ratio = params.pullbackRatio;
					holeArea = surfaceArea * (ratio * ratio);

					//Eigen::MatrixXd holeVertices;
					//auto rot90 = [](const Eigen::Vector2d& v) { return Eigen::Vector2d(-v.y(), v.x()); };
					//// get hole vertices, this is the inner polygon and we will subtract
					//// its face from the total surface area

					//// get the halfspaces for the face
					//auto hs = get_poly_half_spaces(verts2D);

					//holeVertices.resize(verts2D.rows(), 2);

					//Eigen::Vector2d centroid = verts2D.colwise().mean();

					//for (int vIdx = 0; vIdx < verts2D.rows(); ++vIdx) {
					//	int idx1 = (vIdx - 1 + verts2D.rows()) % verts2D.rows();
					//	int idx2 = vIdx;
					//	int idx3 = (vIdx + 1) % verts2D.rows();

					//	Eigen::Vector2d n1 = verts2D.row(idx2) - verts2D.row(idx1);
					//	Eigen::Vector2d n2 = verts2D.row(idx3) - verts2D.row(idx2);

					//	if (n1.norm() == 0 || n2.norm() == 0) {
					//		holeVertices.row(vIdx) = verts2D.row(vIdx);
					//		continue;
					//	}

					//	n1.normalize();
					//	n2.normalize();

					//	Eigen::Vector2d t1 = rot90(n1);
					//	Eigen::Vector2d t2 = rot90(n2);

					//	// current vertex
					//	Eigen::Vector2d xi = verts2D.row(idx2);

					//	// ensure they are inward the origin should lie above the half space
					//	if (t1.dot(centroid) > t1.dot(xi)) t1 = -t1;
					//	if (t2.dot(centroid) > t2.dot(xi)) t2 = -t2;

					//	Eigen::Vector2d normal = t1 + t2;
					//	if (normal.norm() > 1e-12) normal = t2;
					//	else normal.normalize();

					//	double smax = std::numeric_limits<double>::infinity();
					//	for (const auto& L : hs) {
					//		double num = L.d - L.n.dot(xi);  // distance to line along its normal (>=0 inside)
					//		double den = L.n.dot(normal);    // how much moving along b approaches that line
					//		//if (den > 1e-12) smax = std::min(smax, num / (- den));
					//		if (den < -1e-12) {               // moving TOWARD violating the line
					//			smax = std::min(smax, num / (-den));
					//		}
					//	}
					//	if (!std::isfinite(smax)) smax = 0.0;

					//	holeVertices.row(vIdx) = xi - smax * normal;
					//}

					//// estimate the area of the hole
					//holeArea = polygon_area(holeVertices);
				}
				double solidFaceVol = std::max(0.0, surfaceArea - holeArea) * params.thickness;
				solidVolume += solidFaceVol;

			}
			idxOffset += order + 1;
		}

	} while (cla.inc());
		
	// now estimate porosity, fixed for a cube 10x10x10
	double porosity = 1 - (solidVolume / 1000.0f);
	out.porosity = porosity;

	// estimate connectivity from the graph
	double connectivity = graph->find_longest_network_new();
	//out.connectivity = connectivity / seeds.size();
};

void save_csv(
	const std::string fileName,
	const std::vector<SimOutput>& outputs,
	const SamplingType& t) {

	std::ofstream fout;
	fout.open(fileName); 

	// Header
	switch (t) {
	case random:
		fout << "SeedDensity,Thickness,PullbackRatio,Porosity,Connectivity\n";
		break;
	case uniform:
		fout << "Id,Thickness,MinRadius,MaxRadius,Openess,Porosity,Volume,TotalSurface,SurfaceToVolume,Connectivity Density, Local Thickness, Local Thickness Std, Local Separation, Local Separation Std, trabecular Nr, Anisotropy\n";
		break;
	case varied:
		fout << "MinRadius,MaxRadius,Thickness,PullbackRatio,Porosity,Connectivity\n";
		break;
	}

	// Values
	for (size_t i = 0; i < outputs.size(); i++) {
		switch (t) {
		case random:
			//fout << inputs[i].p1 << "," << inputs[i].thickness << "," << inputs[i].pullbackRatio << ","
			//	<< outputs[i].porosity << "," << outputs[i].connectivity << "\n";
			break;
		case uniform:
			fout << outputs[i].id << "," << outputs[i].thickness << "," << outputs[i].minRadius << ","
				<< outputs[i].maxRadius << "," << outputs[i].openess << "," << outputs[i].porosity << "," << outputs[i].volume << "," << outputs[i].totalSurface << "," << outputs[i].surfaceToVolume << "," << outputs[i].connectivityDenisty << "," << outputs[i].localThickness << "," << outputs[i].localThicknessStd << "," << outputs[i].localSeparation << "," << outputs[i].localSeparationStd << "," << outputs[i].trabecularNr << "," << outputs[i].anisotropy << "\n";
			break;
		case varied:
			//fout << inputs[i].p1 << "," << inputs[i].p2 << "," << inputs[i].thickness << "," << inputs[i].pullbackRatio << ","
			//	<< outputs[i].porosity << "," << outputs[i].connectivity << "\n";
			break;
		}
	}
	fout.close();
};

std::vector<std::vector<float>> read_csv(const std::string& fileName) {

	std::vector<std::vector<float>> data;

	std::ifstream file(fileName);

	if (!file.is_open()) {
		std::cerr << "Failed to open file: " << fileName << std::endl;
		return data;
	}

	std::string line;

	// skip first line
	std::getline(file, line);

	while (std::getline(file, line)) {

		// row contains the values of a single line
		std::vector<float> row;
		std::stringstream ss(line);
		std::string cell;

		 //skip first cell which is the id
		//std::getline(ss, cell, ';');

		while (std::getline(ss, cell, ',')) {
			row.push_back(std::stof(cell));
		}

		data.push_back(row);
	}
	file.close();
	return data;
};
