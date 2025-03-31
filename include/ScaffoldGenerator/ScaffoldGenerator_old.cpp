#include "ScaffoldGenerator/ScaffoldGenerator.h"
#include "Utils/Utils.h"
#include "Visualize/visualize.h"
#include <vtkNamedColors.h>
#include <vtkAppendFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkImplicitModeller.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkContourFilter.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkKdTreePointLocator.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <voro++.hh>
#include <random>
#include <Eigen/Dense>
#include <Optimization/bfgs.h>
#include <Optimization/objective.h>

ScaffoldGenerator::ScaffoldGenerator(
	const double& xMin,
	const double& xMax,
	const double& yMin,
	const double& yMax,
	const double& zMin,
	const double& zMax,
	const int& nX, const int& nY, const int& nZ
) : xMin(xMin), xMax(xMax), yMin(yMin), yMax(yMax), zMin(zMin), zMax(zMax), nX(nX), nY(nY), nZ(nZ){};

ScaffoldGenerator::ScaffoldGenerator(
	const std::string& meshFile,
	const int& nX, const int& nY, const int& nZ
	) : nX(nX), nY(nY), nZ(nZ), meshFile(meshFile) {};

ScaffoldGenerator::ScaffoldGenerator(
	const std::vector<std::array<double, 3>>& seeds,
	const std::array<float, 6>& bounds,
	const std::array<int, 3>& blockDim
) : seeds(seeds) {

	xMin = bounds[0];
	xMax = bounds[1];
	yMin = bounds[2];
	yMax = bounds[3];
	zMin = bounds[4];
	zMax = bounds[5];

	nX = blockDim[0];
	nY = blockDim[1];
	nZ = blockDim[2];
};



void ScaffoldGenerator::generate_mesh(
	const double& thickness,
	const std::string& fileName) {

	vtkNew<vtkNamedColors> colors;

	// decide mesh resolution
	int dim[3] = {};

	if (thickness < 0.3) {
		dim[0] = 300;
		dim[1] = 300;
		dim[2] = 300;
	}
	else if (0.3 <= thickness && thickness < 0.5) {
		dim[0] = 200;
		dim[1] = 200;
		dim[2] = 200;
	}
	else {
		dim[0] = 100;
		dim[1] = 100;
		dim[2] = 100;
	}

	// build along line
	vtkNew<vtkImplicitModeller> implictModeller;
	implictModeller->AddInputData(scaffoldMesh);
	implictModeller->SetMaximumDistance(thickness / 2 * 1.01);
	implictModeller->SetSampleDimensions(dim[0], dim[1], dim[2]);
	implictModeller->SetModelBounds(xMin, xMax, yMin, yMax, zMin, zMax);

	// extract isosurface
	vtkNew<vtkContourFilter> isoFilter;
	isoFilter->SetInputConnection(implictModeller->GetOutputPort());
	isoFilter->SetValue(0, thickness / 2);
	isoFilter->Update();

	// update normals
	vtkNew<vtkPolyDataNormals> norms;
	//vtkPolyDataNormals* norms = vtkPolyDataNormals::New();

	norms->SetInputConnection(isoFilter->GetOutputPort());
	norms->ComputePointNormalsOff();
	norms->ComputeCellNormalsOn();
	norms->ConsistencyOn();
	norms->AutoOrientNormalsOn();
	norms->Update();
	
	// use window sinc filter
	//vtkNew<vtkWindowedSincPolyDataFilter> sincSmoother;
	//sincSmoother->SetInputConnection(norms->GetOutputPort());
	//sincSmoother->SetNumberOfIterations(30);
	//sincSmoother->BoundarySmoothingOn();
	//sincSmoother->FeatureEdgeSmoothingOff();
	//sincSmoother->Update();

	vtkNew<vtkSTLWriter> stlWriter;
	stlWriter->SetFileName(fileName.c_str());
	//stlWriter->SetInputConnection(sincSmoother->GetOutputPort());
	stlWriter->SetInputConnection(norms->GetOutputPort());
	stlWriter->Write();
};

void ScaffoldGenerator::get_seeds(std::vector<std::array<double,3>>& outSeeds) {
	outSeeds = seeds;
};

void ScaffoldGenerator::process_triangles(
	std::vector<std::array<double, 3 >>& seeds,
	std::vector<std::array<double, 3>>& bCenters,
	std::vector<std::array<double, 3>>& normals) {

	for (vtkIdType i = 0; i < meshContainer->GetNumberOfCells(); ++i) {
		vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(meshContainer->GetCell(i));
		if (!triangle) continue;

		double p1[3], p2[3], p3[3];
		triangle->GetPoints()->GetPoint(0, p1);
		triangle->GetPoints()->GetPoint(1, p2);
		triangle->GetPoints()->GetPoint(2, p3);

		//trianglePt.push_back({p1[0], p1[1], p1[2]});

		// convert pts to Eigen
		Eigen::Vector3d pe1{ p1[0], p1[1], p1[2] };
		Eigen::Vector3d pe2{ p2[0], p2[1], p2[2] };
		Eigen::Vector3d pe3{ p3[0], p3[1], p3[2] };

		// estimate the barycenter
		double bc1 = (p1[0] + p2[0] + p3[0]) / 3;
		double bc2 = (p1[1] + p2[1] + p3[1]) / 3;
		double bc3 = (p1[2] + p2[2] + p3[2]) / 3;
		std::array<double, 3> bc = { bc1, bc2, bc3 };
		//Eigen::Vector3d bc{ bc1, bc2, bc3 };
		bCenters.push_back(bc);
		seeds.push_back(bc);

		Eigen::Vector3d normal = (pe2 - pe1).cross(pe3 - pe1);
		normal.normalize();
		normals.push_back({ normal[0], normal[1], normal[2] });
	}
};

void ScaffoldGenerator::get_bounds(
	float& xMn, float& xMx, float& yMn, float& yMx, float& zMn, float& zMx) {
	xMn = this->xMin;
	xMx = this->xMax;
	yMn = this->yMin;
	yMx = this->yMax;
	zMn = this->zMin;
	zMx = this->zMax;
};

// ---------------------------------------------------------------------------------------------------
// simple random inside a box
Random::Random(
	const double& xMin,
	const double& xMax,
	const double& yMin,
	const double& yMax,
	const double& zMin,
	const double& zMax,
	const int& nX, const int& nY, const int& nZ,
	const int& seedNr 
	) : ScaffoldGenerator(xMin, xMax, yMin, yMax, zMin, zMax, nX, nY, nZ), seedNr(seedNr) {
	seeds.reserve(seedNr);
};

void Random::generate_seeds() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> disX(xMin + 0.1, xMax - 0.1);
	std::uniform_real_distribution<> disY(yMin + 0.1, yMax - 0.1);
	std::uniform_real_distribution<> disZ(zMin + 0.1, zMax - 0.1);

	for (int i{ 0 }; i < seedNr; i++) {

		double x = disX(gen);
		double y = disY(gen);
		double z = disZ(gen);

		seeds.push_back({ x, y, z });
	}
};

void Random::generate_voro(const int& regSteps) {

	generate_seeds();

	con = new voro::container(
		xMin, xMax, yMin, yMax, zMin, zMax, nX, nY, nZ, false, false, false, 16);

	std::cout << "Looping for steps: " << regSteps << std::endl;
	std::vector<vtkSmartPointer<vtkPolyData>> polys;

	for (int step{ 0 }; step < regSteps; ++step) {

		std::cout << "Regularization Step: " << step + 1 << std::endl;
		
		std::vector<std::array<double, 3>> newSeeds;

		// clear previous container
		con->clear();

		for (int i{ 0 }; i < seeds.size(); i++) {
			con->put(i, seeds[i][0], seeds[i][1], seeds[i][2]);
		}

		voro::c_loop_all cla(*con);
		voro::voronoicell_neighbor cell;

		if (cla.start()) do if (con->compute_cell(cell, cla)) {

			//std::cout << "Entered" << std::endl;
			// vector to store the cell vertexes in global coordinates, cell neighbors and face vertices
			std::vector<double> cellVertices;
			std::vector<int> faceVertices;
			std::vector<int> cellNeighs;

			int seedId = cla.pid();
			// get position of seed and store it to an array
			double px = 0.0, py = 0.0, pz = 0.0;

			cla.pos(px, py, pz);

			// --------------------

			// cell vertices in global system
			cell.vertices(px, py, pz, cellVertices);

			// get cell faces and neighbors
			cell.face_vertices(faceVertices);
			cell.neighbors(cellNeighs);

			std::array<double, 3> centroid;
			cell.centroid(centroid[0], centroid[1], centroid[2]);

			centroid[0] += px;
			centroid[1] += py;
			centroid[2] += pz;

			newSeeds.push_back(centroid);

			//seeds[seedId] = centroid;

			vtkSmartPointer<vtkPolyData> poly = cell_2_vtk(cellNeighs, cellVertices, faceVertices);
			polys.push_back(poly);

		} while (cla.inc());

		//seeds.clear();

		if (step < regSteps - 1) {
			polys.clear();
		}
		else {
			seeds = newSeeds;
		}
	}

	vtkNew<vtkAppendPolyData> appendFilter;
	for (int i = 0; i < polys.size(); i++) {
		appendFilter->AddInputData(polys[i]);
	}
	appendFilter->Update();
	//scaffoldMesh = vtkSmartPointer<vtkPolyData>::New();
	scaffoldMesh = appendFilter->GetOutput();
};

void Random::get_seeds(std::vector<std::array<double, 3>>& outSeeds) {
	outSeeds = seeds;
};

// ---------------------------------------------------------------------------------------------------
// Random inside a mesh

RandomWall::RandomWall(
	const std::string fileName, const int& seedNr,
	const int& nX, const int& nY, const int& nZ, const int& neihgbors
) : ScaffoldGenerator(fileName, nX, nY, nZ), seedNr(seedNr), neighbors(neihgbors) {
	
	
	std::string stlPath = fileName;
	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(stlPath.c_str());
	reader->Update();

	meshContainer = vtkSmartPointer<vtkPolyData>::New();
	meshContainer = reader->GetOutput();

	double stlBounds[6];
	meshContainer->GetBounds(stlBounds);

	double scaleFactor = 1.1;
	xMin = stlBounds[0] * scaleFactor;
	xMax = stlBounds[1] * scaleFactor;
	yMin = stlBounds[2] * scaleFactor;
	yMax = stlBounds[3] * scaleFactor;
	zMin = stlBounds[4] * scaleFactor;
	zMax = stlBounds[5] * scaleFactor;

	seeds.resize(seedNr);
}

void RandomWall::generate_seeds() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> disX(xMin + 0.1, xMax - 0.1);
	std::uniform_real_distribution<> disY(yMin + 0.1, yMax - 0.1);
	std::uniform_real_distribution<> disZ(zMin + 0.1, zMax - 0.1);

	int idx{ 0 };
	while ( idx < seedNr) {

		double x = disX(gen);
		double y = disY(gen);
		double z = disZ(gen);
		Eigen::Vector3d pt{ x, y, z };
		Eigen::Vector3d dir = { 1.0, 1.0, 1.0 };

		if (is_inside_mesh(meshContainer, pt)) {
			idx++;
			std::array<double, 3> t{ x, y, z };
			seeds.push_back(t);
		}
	}
};

void RandomWall::generate_voro(const int& regSteps) {

	std::vector<std::array<double, 3>> bCenters;
	std::vector<std::array<double, 3>> normals;

	bool kdsurface{ false };

	process_triangles(seeds, bCenters, normals);

	generate_seeds();

	std::cout << "Seeds are generated: " << seeds.size() << std::endl;

	// create kdtree to locate points 
	// Create and build the k-d tree
	vtkSmartPointer<vtkKdTreePointLocator> kdTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
	if (kdsurface) {
		kdTree->SetDataSet(meshContainer);
	}
	else {
		vtkSmartPointer<vtkPoints> bcvtk = vtkSmartPointer<vtkPoints>::New();
		for (int i{ 0 }; i < bCenters.size(); i++) {
			bcvtk->InsertNextPoint(bCenters[i][0], bCenters[i][1], bCenters[i][2]);
		}

		vtkNew<vtkPolyData> bcvtkdata;
		bcvtkdata->SetPoints(bcvtk);
		kdTree->SetDataSet(bcvtkdata);
	}
	kdTree->BuildLocator();

	// create container
	con = new voro::container(
		xMin, xMax, yMin, yMax, zMin, zMax, nX, nY, nZ, false, false, false, 16);

	for (int i{ 0 }; i < seeds.size(); i++) {
		con->put(i, seeds[i][0], seeds[i][1], seeds[i][2]);
	}

	MeshWall wall(meshContainer, normals, kdTree, neighbors);
	con->add_wall(wall);

	voro::c_loop_all cla(*con);
	voro::voronoicell_neighbor cell;
	std::vector<vtkSmartPointer<vtkPolyData>> polys;

	std::cout << "Looping in cells" << std::endl;
	if (cla.start()) do if (con->compute_cell(cell, cla)) {

		//std::cout << "Entered" << std::endl;
		// vector to store the cell vertexes in global coordinates, cell neighbors and face vertices
		std::vector<double> cellVertices;
		std::vector<int> faceVertices;
		std::vector<int> cellNeighs;

		int seedId = cla.pid();
		// get position of seed and store it to an array
		double px = 0.0, py = 0.0, pz = 0.0;

		cla.pos(px, py, pz);

		//
		// --------------------

		// cell vertices in global system
		cell.vertices(px, py, pz, cellVertices);

		// get cell faces and neighbors
		cell.face_vertices(faceVertices);
		cell.neighbors(cellNeighs);

		vtkSmartPointer<vtkPolyData> poly = cell_2_vtk(cellNeighs, cellVertices, faceVertices);
		polys.push_back(poly);

	} while (cla.inc());

	delete con;

	vtkNew<vtkAppendPolyData> appendFilter;
	for (int i = 0; i < polys.size(); i++) {
		appendFilter->AddInputData(polys[i]);
	}
	appendFilter->Update();
	//scaffoldMesh = vtkSmartPointer<vtkPolyData>::New();
	scaffoldMesh = appendFilter->GetOutput();
}

void RandomWall ::get_seeds(std::vector<std::array<double, 3>>& outSeeds) {
	outSeeds = seeds;
};


// ---------------------------------------------------------------------------------
// Volume optimization

VolOpt::VolOpt(
	const Eigen::VectorXd targets,
	const Eigen::VectorXd initV,
	const std::array<double, 6>& bounds,
	std::function<void(const std::string&)> logCallback) :
		wInit(initV), targetVols(targets), bounds(bounds), log_callback(logCallback) {

	xMin = bounds[0];
	xMax = bounds[1];
	yMin = bounds[2];
	yMax = bounds[3];
	zMin = bounds[4];
	zMax = bounds[5];
	seedNr = targets.size();
};

void VolOpt::generate_random_seeds() {

	std::cout << "Generate Random Seeds" << std::endl;
	std::random_device rd;
	std::mt19937 gen(rd());

	// Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> disx(xMin + 0.1, xMax - 0.1);
	std::uniform_real_distribution<> disy(yMin + 0.1, yMax - 0.1);
	std::uniform_real_distribution<> disz(zMin + 0.1, zMax - 0.1);

	for (int i = 0; i < seedNr; i++) {
		std::array<double, 3> row{ disx(gen), disy(gen), disz(gen) };
		std::cout << disx(gen) << " " << disy(gen) << " " << disz(gen) << std::endl;
		currSeeds.push_back(row);
	}
};

void VolOpt::generate_random_container_seeds(vtkSmartPointer<vtkPolyData>& container) {
	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> disX(xMin + 0.1, xMax - 0.1);
	std::uniform_real_distribution<> disY(yMin + 0.1, yMax - 0.1);
	std::uniform_real_distribution<> disZ(zMin + 0.1, zMax - 0.1);

	int idx{ 0 };
	while (idx < seedNr) {
		double x = disX(gen);
		double y = disY(gen);
		double z = disZ(gen);
		Eigen::Vector3d pt{ x, y, z };
		Eigen::Vector3d dir = { 1.0, 1.0, 1.0 };
		if (is_inside_mesh(container, pt)) {
			idx++;
			std::array<double, 3> t{ x, y, z };
			currSeeds.push_back(t);
		}
	}
};


void VolOpt::loop(const int regSteps) {

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
		//log_callback("" + std::to_string(func.volError));
		std::cout << "Volume Error: " << func.volError << std::endl;
		//std::cout << func.w.transpose() << std::endl;

		// get the final cells
		// Convert weights to radii
		Eigen::VectorXd radii = convert_radii(func.w);

		std::cout << "Creating Container" << std::endl;
		// create container
		voro::container_poly* conp = new voro::container_poly(
			bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5],
			10, 10, 10,
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

		//log_callback("max diff between centroids and previous seeds: " + std::to_string(maxDiff));

		if (maxDiff < 0.001) {
			//log_callback("Regularization ending!");
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

		//log_callback("\n ----------- \n");
	}
	vtkNew<vtkAppendPolyData> appendFilter;
	for (int i = 0; i < polys.size(); i++) {
		appendFilter->AddInputData(polys[i]);
	}
	appendFilter->Update();
	//scaffoldMesh = vtkSmartPointer<vtkPolyData>::New();
	scaffoldMesh = appendFilter->GetOutput();
};

void VolOpt::generate_mesh(
	const double& thickness,
	const std::string& fileName) {

	std::cout << "Generating Mesh..." << std::endl;

	vtkNew<vtkNamedColors> colors;

	// decide mesh resolution
	int dim[3] = {};

	if (thickness < 0.3) {
		dim[0] = 300;
		dim[1] = 300;
		dim[2] = 300;
	}
	else if (0.3 <= thickness && thickness < 0.5) {
		dim[0] = 200;
		dim[1] = 200;
		dim[2] = 200;
	}
	else {
		dim[0] = 100;
		dim[1] = 100;
		dim[2] = 100;
	}

	// build along line
	vtkNew<vtkImplicitModeller> implictModeller;
	implictModeller->AddInputData(scaffoldMesh);
	implictModeller->SetMaximumDistance(thickness / 2 * 1.01);
	implictModeller->SetSampleDimensions(dim[0], dim[1], dim[2]);
	implictModeller->SetModelBounds(xMin, xMax, yMin, yMax, zMin, zMax);

	// extract isosurface
	vtkNew<vtkContourFilter> isoFilter;
	isoFilter->SetInputConnection(implictModeller->GetOutputPort());
	isoFilter->SetValue(0, thickness / 2);
	isoFilter->Update();

	// update normals
	vtkNew<vtkPolyDataNormals> norms;
	//vtkPolyDataNormals* norms = vtkPolyDataNormals::New();

	norms->SetInputConnection(isoFilter->GetOutputPort());
	norms->ComputePointNormalsOff();
	norms->ComputeCellNormalsOn();
	norms->ConsistencyOn();
	norms->AutoOrientNormalsOn();
	norms->Update();

	// use window sinc filter
	//vtkNew<vtkWindowedSincPolyDataFilter> sincSmoother;
	//sincSmoother->SetInputConnection(norms->GetOutputPort());
	//sincSmoother->SetNumberOfIterations(30);
	//sincSmoother->BoundarySmoothingOn();
	//sincSmoother->FeatureEdgeSmoothingOff();
	//sincSmoother->Update();

	vtkNew<vtkSTLWriter> stlWriter;
	stlWriter->SetFileName(fileName.c_str());
	//stlWriter->SetInputConnection(sincSmoother->GetOutputPort());
	stlWriter->SetInputConnection(norms->GetOutputPort());
	stlWriter->Write();
};

void VolOpt::get_seeds(std::vector<std::array<double, 3>>& outSeeds) {
	outSeeds = currSeeds;
};