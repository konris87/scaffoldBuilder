//#define EIGEN_DONT_PARALLELIZE      // Prevents thread explosion with OpenMP
//#define EIGEN_MAX_ALIGN_BYTES 0

#include "SeedGenerator/Container.h"
#include "SeedGenerator/SeedGenerator.h"
//#include "ScaffoldGenerator/ScaffoldGenerator.h"
//#include "ScaffoldGenerator/ScaffoldGeneratorBox.h"
#include "MonteSimulation.h"
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
		fout << "Id,Thickness,Radius,MaxRadius,Openess,Porosity,Volume,TotalSurface,SurfaceToVolume,Connectivity Density, Local Thickness, Local Thickness Std, Local Separation, Local Separation Std, trabecular Nr, Anisotropy\n";
		break;
	case varied:
		fout << "Id,MinThickness,MaxThickness,MinRadius,MaxRadius,Openess,Porosity,Volume,TotalSurface,SurfaceToVolume,Connectivity Density, Local Thickness, Local Thickness Std, Local Separation, Local Separation Std, trabecular Nr, Anisotropy\n";
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
			fout << outputs[i].id << "," << outputs[i].minThickness << "," << outputs[i].minRadius << "," << outputs[i].openess << "," << outputs[i].porosity << "," << outputs[i].volume << "," << outputs[i].totalSurface << "," << outputs[i].surfaceToVolume << "," << outputs[i].connectivityDenisty << "," << outputs[i].localThickness << "," << outputs[i].localThicknessStd << "," << outputs[i].localSeparation << "," << outputs[i].localSeparationStd << "," << outputs[i].trabecularNr << "," << outputs[i].anisotropy << "\n";
			break;
		case varied:
			fout << outputs[i].id << "," << outputs[i].minThickness << "," << outputs[i].maxThickness << "," << outputs[i].minRadius << "," << outputs[i].maxRadius << "," << outputs[i].openess << "," << outputs[i].porosity << "," << outputs[i].volume << "," << outputs[i].totalSurface << "," << outputs[i].surfaceToVolume << "," << outputs[i].connectivityDenisty << "," << outputs[i].localThickness << "," << outputs[i].localThicknessStd << "," << outputs[i].localSeparation << "," << outputs[i].localSeparationStd << "," << outputs[i].trabecularNr << "," << outputs[i].anisotropy << "\n";
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
