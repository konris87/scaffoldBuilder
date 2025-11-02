#include <json.hpp>
#include <random>
#include <cmath>
#include <array>
#include <sstream>
#include <iostream>
#include <fstream>
#include <chrono>
#include "buildScaffold.h"
#include <Optimization/bfgs.h>
#include <Optimization/Objective.h>
#include <voro++.hh>
#include <Misc/cell_process.h>
#include <Visualize/visualize.h>

buildScaffold::buildScaffold(
	const Eigen::VectorXd targets,
	const Eigen::VectorXd initV,
	const std::string file,
	std::function<void(const std::string&)> logCallback) :
	settingsFile(file),
	wInit(initV),
	targetVols(targets),
	log_callback(logCallback) {

	load_settings();

	if (genOption == "generate") {

		std::cout << "Generating Random Points!" << std::endl;

		// Will be used to obtain a seed for the random number engine
		std::random_device rd;
		std::mt19937 gen(rd());

		// Standard mersenne_twister_engine seeded with rd()
		std::uniform_real_distribution<> disx(xMin + 0.1, xMax - 0.1);
		std::uniform_real_distribution<> disy(yMin + 0.1, yMax - 0.1);
		std::uniform_real_distribution<> disz(zMin + 0.1, zMax - 0.1);

		for (int i = 0; i < seedNr; i++) {
			std::array<double, 3> row{ disx(gen), disy(gen), disz(gen) };
			std::cout << disx(gen) << " " << disy(gen) << " " << disz(gen) << std::endl;
			currSeeds.emplace_back(row);
		}
	}
	else {
		// vector of cell vertices
		std::fstream myfile("../vertices500.txt");
		int nr_cells{ 0 };
		if (myfile.is_open()) {

			//std::cout << "reading line" << std::endl;

			std::string line;
			while (std::getline(myfile, line)) {
				std::istringstream iss(line);

				double a, b, c;
				iss >> a >> b >> c;
				std::array<double, 3> lineVec{ a, b ,c };

				currSeeds.emplace_back(lineVec);

				nr_cells++;

			}
			myfile.close();
		}

	}

};

void buildScaffold::loop() {

	std::cout << wInit << std::endl;
	std::cout << targetVols << std::endl;

	for (int step = 0; step < regSteps; step++) {

		if (log_callback) {
			log_callback("Regularization step: " + std::to_string(step + 1));
		}
		else {
			std::cout << "Regularization step : " << step + 1 << std::endl;
		}	

		// create an instance of the objective function
		myFunc func(
			bounds,
			currSeeds,
			wInit,
			targetVols
		);

		double f0;
		Eigen::VectorXd g0 = Eigen::VectorXd::Zero(wInit.size());
		f0 = func(wInit, g0);

		// estimate the tolerance based on the initial step gradient
		double tol = std::min(0.1 * 1e-2 * targetVols.minCoeff(), 0.1 * 1e-2 * targetVols.minCoeff() / g0.array().abs().maxCoeff());

		if (log_callback) {
			log_callback("Optimality Tolerance: " + std::to_string(tol));
		}
		else {
			std::cout << "Optimality Tolerance: " << tol << std::endl;
		}		

		// create a solver object
		BFGS<myFunc> solver(func);

		solver.maxIter = 1000;
		//solver.epsilon = 0.00056;
		solver.epsilon = tol;
		solver.minimize(wInit);

		//std::cout << func.w << std::endl;
		log_callback("" + std::to_string(func.volError));
		std::cout << "Volume Error: " << func.volError << std::endl;
		//std::cout << func.w.transpose() << std::endl;

		// get the final cells
		// Convert weights to radii
		Eigen::VectorXd radii = convert_radii(func.w);

		// create container
		voro::container_poly* conp = new voro::container_poly(
			bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5],
			6, 6, 6,
			false, false, false, 16
		);

		for (int i = 0; i < currSeeds.size(); i++) {
			conp->put(i, func.centroids[i][0], func.centroids[i][1], func.centroids[i][2], radii[i]);
		}

		Eigen::VectorXd currVols = Eigen::VectorXd::Zero(currSeeds.size());
		Eigen::VectorXd currCosts = Eigen::VectorXd::Zero(currSeeds.size());

		// loop in cells
		voro::c_loop_all cla(*conp);
		voro::voronoicell_neighbor cell;

		if (cla.start()) do if (conp->compute_cell(cell, cla)) {

			// vector to store the cell vertexes in global coordinates, cell neighbors and face vertices
			std::vector<double> cellVertices;
			std::vector<int> faceVertices;
			std::vector<int> cellNeighs;

			int seedId = cla.pid();

			// get position of seed and store it to an array
			double x = 0.0, y = 0.0, z = 0.0;
			cla.pos(x, y, z);

			// cell vertices in global system
			cell.vertices(x, y, z, cellVertices);

			// get cell faces and neighbors
			cell.face_vertices(faceVertices);
			cell.neighbors(cellNeighs);

			vtkSmartPointer<vtkPolyData> poly = cell_2_vtk(cellNeighs, cellVertices, faceVertices);
			polys.push_back(poly);

		} while (cla.inc());

		delete conp;

		// check if the centroids and the currSeeds are close enough
		double maxDiff{ 0 };
		for (int i{ 0 }; i < wInit.size(); i++) {
			double diff = std::sqrt(
				std::pow(currSeeds[i][0] - func.centroids[i][0], 2) +
				std::pow(currSeeds[i][1] - func.centroids[i][1], 2) +
				std::pow(currSeeds[i][2] - func.centroids[i][2], 2));

			if (maxDiff < diff) { maxDiff = diff; }
		}
		//std::cout << "max diff between centroids and previous seeds: " << maxDiff << std::endl;
		
		log_callback("max diff between centroids and previous seeds: " + std::to_string(maxDiff));

		if (maxDiff < 0.001) {
			log_callback("Regularization ending!");
			break;
		}
		else if (step < regSteps - 1) {
			polys.clear();
		}
		else {
			break;
		}

		// update seeds

		currSeeds = func.centroids;

		log_callback("\n ----------- \n");
	}

	std::cout << "polys size: " << polys.size() << std::endl;

	const std::string fileName = "../data/" + version + ".stl";
	create_mesh(polys, thickness, bounds, fileName);

};

void buildScaffold::load_settings() {

	//log_callback("Loading Settings.. ");
	std::cout << "loading" << std::endl;

	std::ifstream input_file(settingsFile);

    nlohmann::json config;

	input_file >> config;

	xMin = config["Scaffold"]["Domain"]["xMin"];
	xMax = config["Scaffold"]["Domain"]["xMax"];
	yMin = config["Scaffold"]["Domain"]["yMin"];
	yMax = config["Scaffold"]["Domain"]["yMax"];
	zMin = config["Scaffold"]["Domain"]["zMin"];
	zMax = config["Scaffold"]["Domain"]["zMax"];
	
	genOption = config["Scaffold"]["Pores"]["genOption"];
	seedNr = config["Scaffold"]["Pores"]["poreNr"];
	thickness = config["Scaffold"]["thickness"];
	regSteps = config["Scaffold"]["regSteps"];

	version = config["Scaffold"]["version"];

	bounds = std::array<double, 6>{xMin, xMax, yMin, yMax, zMin, zMax};

};
