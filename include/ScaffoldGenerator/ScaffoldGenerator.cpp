#include "ScaffoldGenerator.h"
#include "Visualize/visualize.h"
#include "Wall/Wall.h"
#include "Optimization/objective.h"
#include "Optimization/bfgs.h"
#include <vtkAppendPolyData.h>
#include <vtkNamedColors.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkImplicitModeller.h>
#include <vtkCleanPolyData.h>
#include <vtkSTLWriter.h>
#include <vtkKdTreePointLocator.h>
#include <vtkTriangle.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <voro++.hh>
#include <Eigen/Dense>
#include <Utils/Utils.h>

void ScaffoldGenerator::generate_mesh(
	const double& thickness, vtkSmartPointer<vtkPolyData>& finalPolyData, const std::vector<int>& res) {

	std::cout << "generating mesh" << std::endl;

	// first get the center of the scaffold mesh
	double center[3] = { 0.0, 0.0, 0.0 };

	vtkIdType numPoints = scaffoldMesh->GetNumberOfPoints();

	vtkSmartPointer<vtkPoints> points = scaffoldMesh->GetPoints();
	for (vtkIdType i = 0; i < numPoints; i++) {
		double p[3];
		points->GetPoint(i, p);
		center[0] += p[0];
		center[1] += p[1];
		center[2] += p[2];
	}
	center[0] /= numPoints;
	center[1] /= numPoints;
	center[2] /= numPoints;

	// Transform the polydata
	vtkNew<vtkTransform> transform;
	transform->Translate(center[0], center[1], center[2]);
	float sf = 1.0f -  (thickness + 0.1)  * 0.1f;
	transform->Scale(sf, sf, sf);
	transform->Translate(-center[0], -center[1], -center[2]);

	vtkNew<vtkTransformPolyDataFilter> transformPD;
	transformPD->SetTransform(transform);
	transformPD->SetInputData(scaffoldMesh);
	transformPD->Update();
	scaffoldMesh->DeepCopy(transformPD->GetOutput());

	vtkNew<vtkNamedColors> colors;

	// decide mesh resolution
	int dim[3] = {};

	if (res.empty()) {
		if (thickness < 0.3) {
			dim[0] = 300;
			dim[1] = 300;
			dim[2] = 300;
		}
		else if (0.3 <= thickness && thickness < 0.5) {
			dim[0] = 100;
			dim[1] = 100;
			dim[2] = 100;
		}
		else {
			dim[0] = 100;
			dim[1] = 100;
			dim[2] = 100;
		}
	}
	else {
		// If resolution was provided, use it
		dim[0] = res[0];
		dim[1] = res[1];
		dim[2] = res[2];
	}

	std::cout << "Dims: " << dim[0] << " " << dim[1] << " " << dim[2] << std::endl;
	
	// build along line
	vtkNew<vtkImplicitModeller> implictModeller;
	implictModeller->AddInputData(scaffoldMesh);
	implictModeller->SetMaximumDistance(thickness / 2 * 1.01);
	implictModeller->SetSampleDimensions(dim[0], dim[1], dim[2]);
	implictModeller->SetModelBounds(
		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

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

	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputData(norms->GetOutput());
	cleaner->Update();
	finalPolyData = cleaner->GetOutput();

	// use window sinc filter
	//vtkNew<vtkWindowedSincPolyDataFilter> sincSmoother;
	//sincSmoother->SetInputConnection(norms->GetOutputPort());
	//sincSmoother->SetNumberOfIterations(30);
	//sincSmoother->BoundarySmoothingOn();
	//sincSmoother->FeatureEdgeSmoothingOff();
	//sincSmoother->Update();
	//finalPolyData = norms->GetOutput();

	//vtkNew<vtkSTLWriter> stlWriter;
	//stlWriter->SetFileName(fileName.c_str());
	////stlWriter->SetInputConnection(sincSmoother->GetOutputPort());
	//stlWriter->SetInputConnection(norms->GetOutputPort());
	//stlWriter->Write();

};

// -------------------------------------------------------------------------------------
// Scaffold Generator inside a rectangle
ScaffoldGeneratorBox::ScaffoldGeneratorBox(
	std::vector<std::array<double, 3>>& seeds,
	const std::array<float, 6>& bounds,
	const std::array<int, 3>& blockDim) : ScaffoldGenerator(seeds, bounds, blockDim) {

	//std::cout << blockDim[0] << blockDim[1] << blockDim[2] << std::endl;
	//std::cout << bounds[0] << bounds[1] << bounds[2] << std::endl;

};

void ScaffoldGeneratorBox::generate_voro(const int regSteps) {

	con = new voro::container(
		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5],
		blockDim[0], blockDim[1], blockDim[2], false, false, false, 16);

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

			if (step == regSteps - 1) {
				vtkSmartPointer<vtkPolyData> poly = cell_2_vtk(cellNeighs, cellVertices, faceVertices);
				polys.push_back(poly);
			}
		} while (cla.inc());

		seeds = newSeeds;
	}

	std::cout << "Polys Size: " << polys.size() << std::endl;
	vtkNew<vtkAppendPolyData> appendFilter;
	for (int i = 0; i < polys.size(); i++) {
		appendFilter->AddInputData(polys[i]);
	}
	appendFilter->Update();
	//scaffoldMesh = vtkSmartPointer<vtkPolyData>::New();
	scaffoldMesh = appendFilter->GetOutput();
};

// ------------------------------------------
// generate mesh inside wall

ScaffoldGeneratorWall::ScaffoldGeneratorWall(
	std::vector<std::array<double, 3>>& seeds,
	vtkSmartPointer<vtkPolyData>& containerPoly,
	const std::array<int, 3>& blockDim, const int neighbors,
	const float minDist) : 
	ScaffoldGenerator(seeds, { 0,0,0,0,0,0 }, blockDim), containerMesh(containerPoly), neighbors(neighbors), minDist(minDist) {

	double bds[6];

	containerMesh->GetBounds(bds);
	bounds[0] = bds[0] - 1.0;
	bounds[1] = bds[1] + 1.0;
	bounds[2] = bds[2] - 1.0;
	bounds[3] = bds[3] + 1.0;
	bounds[4] = bds[4] - 1.0;
	bounds[5] = bds[5] + 1.0;

	_process_triangles();
};

void ScaffoldGeneratorWall::generate_voro(const int regSteps) {

	bool kdsurface{ false };

	std::cout << "Seed Nr: " << seeds.size() << std::endl;
	std::cout << "Barycenter Nr: " << bCenters.size() << std::endl;
	std::cout << "Normals Nr: " << normals.size() << std::endl;

	// create kdtree to locate points 
	// Create and build the k-d tree
	vtkSmartPointer<vtkKdTreePointLocator> kdTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
	if (kdsurface) {
		kdTree->SetDataSet(containerMesh);
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
		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5],
		blockDim[0], blockDim[1], blockDim[2], false, false, false, 16);

	MeshWall wall(containerMesh, normals, kdTree, neighbors);
	con->add_wall(wall);

	std::vector<vtkSmartPointer<vtkPolyData>> polys;
	std::cout << "Regularization Steps: " << regSteps << std::endl;
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

			if (is_inside_mesh(containerMesh, {centroid[0], centroid[1], centroid[2]})) {
				newSeeds.push_back(centroid);
			}
			
			//seeds[seedId] = centroid;

			//if (step == regSteps - 1) {
			//	vtkSmartPointer<vtkPolyData> poly = cell_2_vtk(cellNeighs, cellVertices, faceVertices);
			//	polys.push_back(poly);
			//}
		} while (cla.inc());

		seeds = newSeeds;

	}
	
	std::cout << "------------------" << std::endl;

	std::cout << "after reg seeds: " << seeds.size() << std::endl;

	// clear previous container
	con->clear();

	for (int i{ 0 }; i < seeds.size(); i++) {
		con->put(i, seeds[i][0], seeds[i][1], seeds[i][2]);
	}

	for (int i{ 0 }; i < bCenters.size(); i++) {
		seeds.push_back({ bCenters[i][0], bCenters[i][1], bCenters[i][2] });
		con->put(i, bCenters[i][0], bCenters[i][1], bCenters[i][2]);
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

};

void ScaffoldGeneratorWall::_process_triangles() {

	std::cout << "processing triangles" << std::endl;
	
	std::vector<std::array<double, 3>> centers;
	std::vector<std::array<double, 3>> tempNormals;
	for (vtkIdType i = 0; i < containerMesh->GetNumberOfCells(); ++i) {

		vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(containerMesh->GetCell(i));
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

		centers.push_back(bc);
		//seeds.push_back(bc);

		Eigen::Vector3d normal = (pe2 - pe1).cross(pe3 - pe1);
		normal.normalize();
		tempNormals.push_back({ normal[0], normal[1], normal[2] });
	}

	// filter centers based on the minimum distance

	// Filtering step
	for (size_t i = 0; i < centers.size(); i++) {
		bool keep = true;
		for (size_t j = 0; j < bCenters.size(); j++) {
			if (distance(centers[i], bCenters[j]) < minDist) {
				keep = false;
				break;
			}
		}
		if (keep) {
			bCenters.push_back(centers[i]);
			normals.push_back(tempNormals[i]);
		}
	}
	std::cout << "Before filtering: " << centers.size() << std::endl;
	std::cout << "After filtering: " << bCenters.size() << std::endl;
};


// ------------------------------------------------------------------
// Volume Optimization

VolOpt::VolOpt(
	std::vector<std::array<double, 3>>& seeds,
	const Eigen::VectorXd targets,
	const Eigen::VectorXd initV,
	const std::array<double, 6>& bounds,
	std::function<void(const std::string&)> logCallback) :
		currSeeds(seeds), targetVols(targets), wInit(initV), bounds(bounds), log_callback(logCallback) {};

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
	vtkSmartPointer<vtkPolyData>& finalPolyData,
	const std::vector<int>& res) {

	vtkNew<vtkNamedColors> colors;

	// decide mesh resolution
	int dim[3] = {};

	if (res.empty()) {
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
	}
	else {
		dim[0] = res[0];
		dim[1] = res[1];
		dim[2] = res[2];
	}
	// build along line
	vtkNew<vtkImplicitModeller> implictModeller;
	implictModeller->AddInputData(scaffoldMesh);
	implictModeller->SetMaximumDistance(thickness / 2 * 1.01);
	implictModeller->SetSampleDimensions(dim[0], dim[1], dim[2]);
	implictModeller->SetModelBounds(
		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

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
	
	finalPolyData = norms->GetOutput();
	//vtkNew<vtkSTLWriter> stlWriter;
	//stlWriter->SetFileName(fileName.c_str());
	////stlWriter->SetInputConnection(sincSmoother->GetOutputPort());
	//stlWriter->SetInputConnection(norms->GetOutputPort());
	//stlWriter->Write();
};

void VolOpt::get_seeds(std::vector<std::array<double, 3>>& outSeeds) {
	outSeeds = currSeeds;
};

// ----------------------------------------------------------------
VolOptWall::VolOptWall(
	std::vector<std::array<double, 3>>& seeds,
	const Eigen::VectorXd targets,
	const Eigen::VectorXd initV,
	vtkSmartPointer<vtkPolyData>& containerPoly,
	int neighbors,
	std::function<void(const std::string&)> logCallback) :
	currSeeds(seeds), targetVols(targets), wInit(initV), containerMesh(containerPoly),
	neighbors(neighbors), log_callback(logCallback) {

	double bds[6];
	containerMesh->GetBounds(bds);

	bounds[0] = bds[0];
	bounds[1] = bds[1];
	bounds[2] = bds[2];
	bounds[3] = bds[3];
	bounds[4] = bds[4];
	bounds[5] = bds[5];

	_process_triangles();
};

void VolOptWall::loop(const int regSteps) {

	// create container
	bool kdsurface{ false };
	vtkSmartPointer<vtkKdTreePointLocator> kdTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
	if (kdsurface) {
		kdTree->SetDataSet(containerMesh);
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

	MeshWall abstractWall(containerMesh, normals, kdTree, neighbors);

	for (int step = 0; step < regSteps; step++) {

		if (log_callback) {
			log_callback("Regularization step: " + std::to_string(step + 1));
		}
		else {
			std::cout << "Regularization step : " << step + 1 << std::endl;
		}

		// create an instance of the objective function
		myFuncWall func(
			abstractWall,
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
		BFGS<myFuncWall> solver(func);

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

		voro::container_poly* conp = new voro::container_poly(
			bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5],
			10, 10, 10,
			false, false, false, 16
		);

		for (int i = 0; i < currSeeds.size(); i++) {
			conp->put(i, func.centroids[i][0], func.centroids[i][1], func.centroids[i][2], radii[i]);
		}

		conp->add_wall(abstractWall);

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

void VolOptWall::_process_triangles() {

	std::cout << "processing triangles" << std::endl;

	for (vtkIdType i = 0; i < containerMesh->GetNumberOfCells(); ++i) {

		vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(containerMesh->GetCell(i));
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

		bCenters.push_back(bc);
		currSeeds.push_back(bc);

		Eigen::Vector3d normal = (pe2 - pe1).cross(pe3 - pe1);
		normal.normalize();
		normals.push_back({ normal[0], normal[1], normal[2] });
	}
};

void VolOptWall::generate_mesh(
	const double& thickness,
	const std::string& fileName, 
	const std::vector<int>& res) {

	vtkNew<vtkNamedColors> colors;

	// decide mesh resolution
	int dim[3] = {};

	if (res.empty()) {
		if (thickness < 0.3) {
			dim[0] = 300;
			dim[1] = 300;
			dim[2] = 300;
		}
		else if (0.3 <= thickness && thickness < 0.5) {
			dim[0] = 100;
			dim[1] = 100;
			dim[2] = 100;
		}
		else {
			dim[0] = 100;
			dim[1] = 100;
			dim[2] = 100;
		}
	}
	else {
		dim[0] = res[0];
		dim[1] = res[1];
		dim[2] = res[2];
	}

	// build along line
	vtkNew<vtkImplicitModeller> implictModeller;
	implictModeller->AddInputData(scaffoldMesh);
	implictModeller->SetMaximumDistance(thickness / 2 * 1.01);
	implictModeller->SetSampleDimensions(dim[0], dim[1], dim[2]);
	implictModeller->SetModelBounds(
		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

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

// -----------------------------------------------------------------
// Scaffold Face Builder

// constructor
ScaffoldGeneratorFaceBox::ScaffoldGeneratorFaceBox(
	std::vector<std::array<double, 3>>& seeds,
	const std::array<float, 6>& bounds,
	const std::array<int, 3>& blockDim,
	double minRad, double maxRad, double edgeSize
) : minHoleRadius(minRad), maxHoleRadius(maxRad), edgeSize(edgeSize), ScaffoldGenerator(seeds, bounds, blockDim) {};

// generate voronoi 

void ScaffoldGeneratorFaceBox::generate_voro() {

	con = new voro::container(
		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5],
		blockDim[0], blockDim[1], blockDim[2], false, false, false, 16
	);

	std::vector<vtkSmartPointer<vtkPolyData>> polys;

	// clear previous container
	con->clear();

	for (int i{ 0 }; i < seeds.size(); i++) {
		con->put(i, seeds[i][0], seeds[i][1], seeds[i][2]);
	}

	voro::c_loop_all cla(*con);
	voro::voronoicell_neighbor cell;

	int cellNr{ 1 };

	if (cla.start()) do if (con->compute_cell(cell, cla)) {

		std::cout << " --------------- " << std::endl;
		
		std::cout << " Processing cell: " << cellNr << std::endl;
		cellNr += 1;

		// vector to store the cell vertexes in global coordinates, cell neighbors and face vertices
		std::vector<double> cellVertices;
		std::vector<int> faceIndices;
		std::vector<int> cellNeighs;

		int seedId = cla.pid();

		// get position of seed and store it to an array
		double px = 0.0, py = 0.0, pz = 0.0;

		cla.pos(px, py, pz);

		//-------------------------------

		// cell vertices in global system
		cell.vertices(px, py, pz, cellVertices);

		// nr of cell vertices
		int cellVertexNr = cellVertices.size() / 3;

		// create also a mapping from local to global
		std::unordered_map <int, int> localToGlobal;

		// get cell faces and neighbors
		cell.face_vertices(faceIndices);
		cell.neighbors(cellNeighs);
		
		// start processing faces 
		// loop in faces 
		int idx{ 0 };

		for (int i = 0; i < cellNeighs.size(); i++) {

			std::cout << " ------------------- " << std::endl;
			std::cout << " Processing Face: " << i << std::endl;
			std::cout << " Neighbour Id: " << cellNeighs[i] << std::endl;

			bool faceProcessed = false;

			// create a vector to store face vertices
			std::vector<double> faceVertices;

			// each face is a vector of values where the first value denotes the number of vertices in the faces
			int order = faceIndices[idx];

			if (cellNeighs[i] < 0) {
				idx += order + 1;

				continue;
			}

			std::cout << "nr of verts: " << order << std::endl;

			std::vector<int> localFace;
			std::vector<int> globalFace;

			for (int j = idx + 1; j < order + idx + 1; j++) {

				int localVertIdx = faceIndices[j];
				localFace.push_back(localVertIdx);

			};
			idx += order + 1;

			// get the face vertices
			for (const auto vIdx : localFace){
				faceVertices.push_back(cellVertices[3 * vIdx]);
				faceVertices.push_back(cellVertices[3 * vIdx + 1]);
				faceVertices.push_back(cellVertices[3 * vIdx + 2]);
				//std::cout << "local v " << vIdx << " : " << cellVertices[3 * vIdx] << " " << cellVertices[3 * vIdx + 1] << " " << cellVertices[3 * vIdx + 2] << std::endl;

			}

			Eigen::MatrixXd verts2D;
			Eigen::Vector3d u;
			Eigen::Vector3d v;
			Eigen::Vector3d origin;

			// first project the vertices
			project_vertices_on_plane(faceVertices, origin, u, v, verts2D);

			// ensure clockwise order
			ensure_ccw(localFace, verts2D);

			// now that we have clockwise order we can check if the vertices are already added
			Eigen::MatrixXd tempVerts(verts2D.rows(), 3);
			back_to_3d(tempVerts, verts2D, origin, u, v);

			//// for each vertex check if is inside the global hash map
			for (int rowIdx{ 0 }; rowIdx < tempVerts.rows(); rowIdx++) {

				//std::cout << tempVerts(rowIdx, 0) << " " << tempVerts(rowIdx, 1) << " " << tempVerts(rowIdx, 2) << std::endl;

				std::array<double, 3> coords = { tempVerts(rowIdx, 0), tempVerts(rowIdx, 1), tempVerts(rowIdx, 2) };

				GlobalVertex key(coords);

				auto it = globalVertexMap.find(key);

				if (it == globalVertexMap.end()) {
					//std::cout << "vertex is not inside " << std::endl;
					//std::cout << "Adding vertex " << globalIndex << " with global vertex coords: " << key.coords[0] << " " << key.coords[1] << " " << key.coords[2] << std::endl;

					// add to the map
					globalVertexMap[key] = globalIndex;

					localToGlobal[rowIdx] = globalIndex;
					
					globalFace.push_back(globalIndex);

					globalVertices.push_back({ coords[0], coords[1], coords[2] });

					globalIndex++;

				}
				else {
					// the vertex is already 
					//std::cout << "Found vertex with global vertex coords: " << key.coords[0] << " " << key.coords[1] << " " << key.coords[2] << " at index: " << it->second << std::endl;
					localToGlobal[rowIdx] = it->second;
					globalFace.push_back(it->second);
				}
			}

			//// print the local and global faces
			//std::cout << " local face: " << std::endl;
			//for (int vIdx{ 0 }; vIdx < localFace.size(); vIdx++) {

			//		std::cout << localFace[vIdx] << " ";

			//}
			//std::cout << "\n global face: " << std::endl;
			//for (int vIdx{ 0 }; vIdx < globalFace.size(); vIdx++) {

			//	std::cout << globalFace[vIdx] << " ";

			//}		

			// using the face with global idxs check if the face has already been checked
			std::vector<int> globalFaceOrdered = globalFace;
			std::sort(globalFaceOrdered.begin(), globalFaceOrdered.end());

			//------------------------------------------------------------------------------------------------
			// if not inside continue with the face
			if (globalFaceMap.find(globalFaceOrdered) == globalFaceMap.end()) {
				
				globalFaceMap.insert(globalFaceOrdered);
				//std::cout << "Added Face" << std::endl;

				Eigen::Vector2d holeCenter;
				double radius1{ 0.0 };
				int resolution{ 50 };
				int holePtNr{ 12 };

				//render_vtk_points(verts2D, "Before Sampling");

				sample_face_polygon(verts2D, holeCenter, radius1, 50, 0.5);

				// 4. We have the center and the radius, we can interpolate each edge to 
				// assist ear clipping triangulation
				Eigen::MatrixXd interpolatedVerts;
				std::vector<int> newLocalFace;

				interpolate_edges(verts2D, interpolatedVerts, newLocalFace, edgeSize);
			
				//render_vtk_points(interpolatedVerts, "After Sampling");
				//for (auto trIdx : newLocalFace) {
				//	std::cout << trIdx << " " << std::endl;
				//}

				//render_vtk_points(interpolatedVerts, "Interpolated");

				// now we can create the hole
				Eigen::MatrixXd holeVertices;
				hole_points(radius1, holeCenter, holePtNr, holeVertices);

				Eigen::Vector2d pair = {0, interpolatedVerts.rows()};

				double minDist = (interpolatedVerts.row(0) - holeVertices.row(0)).norm();

				// 6. Find the connector pair this is the pair between 
				for (int polyIdx{ 0 }; polyIdx < interpolatedVerts.rows(); polyIdx++) {
					for (int holeIdx{0}; holeIdx < holeVertices.rows();	holeIdx++) {
						double dist = (holeVertices.row(holeIdx) - interpolatedVerts.row(polyIdx)).norm();
						if (dist < minDist) {
							pair = {polyIdx, holeIdx + interpolatedVerts.rows()};
							minDist = dist;
						}
					}
				}

				Eigen::MatrixXd finalVerts(interpolatedVerts.rows() + holeVertices.rows(), 2);

				// Copy verts2D
				for (int cIdx = 0; cIdx < interpolatedVerts.rows(); ++cIdx) {
					finalVerts.row(cIdx) = interpolatedVerts.row(cIdx);
				}

				// Copy holeVertices
				for (int cIdx = 0; cIdx < holeVertices.rows(); ++cIdx) {
					finalVerts.row(interpolatedVerts.rows() + cIdx) = holeVertices.row(cIdx);
				}
				
				//std::cout << "Connector: " << pair[0] << " " << pair[1] << std::endl;

				//render_vtk_points(finalVerts, "final vertices");

				// 6. Create a dll
				Cdll holeDll;

				// hold the nodes
				std::vector<Node*> nodes;

				for (int vIdx{ 0 }; vIdx < holeVertices.rows(); vIdx++) {
					Node* ni = holeDll.append(vIdx + interpolatedVerts.rows());
					nodes.push_back(ni);
				}

				// rotate it to start from connector pair[1]
				holeDll.rotate(pair[1] - interpolatedVerts.rows());

				// retrieve the data
				std::vector<int> holeIdxs;
				holeDll.get_data(holeIdxs);

				// print the length
				int length{ 0 };
				holeDll.get_length(length);

				// now create the correct row of indices to pass to the ear clipping algorithm
				std::vector<int> triIdxs;

				//for (auto hIdx : newLocalFace) {
				//		std::cout << hIdx << std::endl;
				//}

				std::cout << " - " << std::endl;

				// first add the indices up to the first point of the pair  
				auto it = std::find(newLocalFace.begin(), newLocalFace.end(), pair[0]);

				int pIdx1 = it - newLocalFace.begin();
				
				for (int cIdx{ 0 }; cIdx < pIdx1 + 1; cIdx++) {
					triIdxs.push_back(newLocalFace[cIdx]);
				}

				// then add the indices of the rotated hole dll
				for (auto hIdx : holeIdxs) {
					triIdxs.push_back(hIdx);
				}

				// then add the next idxs of the face
				for (int cIdx{pIdx1}; cIdx < interpolatedVerts.rows(); cIdx++) {
					//std::cout << cIdx << std::endl;
					triIdxs.push_back(newLocalFace[cIdx]);
				}

				//for (auto trIdx : triIdxs) {
				//	std::cout << trIdx << " " << std::endl;
				//}

				// ear clipping triangulation
				std::vector<std::vector<int>> faceCells;
				ear_clipping(finalVerts, triIdxs, faceCells);

				//render_vtk_face(finalVerts, faceCells, "Face");

				// bring the vertices of the face back to 3d
				Eigen::MatrixXd finalVerts3d(finalVerts.rows(), 3);
				back_to_3d(finalVerts3d, finalVerts, origin, u, v);

				std::unordered_map<int, int> newLocalToGlobal;

				// add the finalVerts to the global list

				for (int rowIdx{ 0 }; rowIdx < finalVerts3d.rows(); rowIdx++) {
					
					std::array<double, 3> coords = {
						finalVerts3d(rowIdx, 0),
						finalVerts3d(rowIdx, 1),
						finalVerts3d(rowIdx, 2),
					};

					GlobalVertex key(coords);

					auto it = globalVertexMap.find(key);

					//std::cout << "Looking for local vertex " << rowIdx << " " << key.coords[0] << " " << key.coords[1] << " " << key.coords[2] << std::endl;

					if (it == globalVertexMap.end()) {
						//std::cout << "vertex is not inside " << std::endl;
						//std::cout << "Adding vertex " << globalIndex << " with global vertex coords: " << key.coords[0] << " " << key.coords[1] << " " << key.coords[2] << std::endl;
						// add to the map
						globalVertexMap[key] = globalIndex;
						localToGlobal[rowIdx] = globalIndex;

						globalVertices.push_back({ coords[0], coords[1], coords[2] });
						globalIndex++;
					}
					else {
						// the vertex is already 
						//std::cout << "Found vertex with global vertex coords: " << key.coords[0] << " " << key.coords[1] << " " << key.coords[2] << " at index: " << it->second << std::endl;
						localToGlobal[rowIdx] = it->second;
					}
				}
				// add the global faces
				for (auto face : faceCells) {
					std::vector<int> gIdxs;
					for (auto vIdx : face) {
						gIdxs.push_back(localToGlobal[vIdx]);
					}
					globalFaces.push_back(gIdxs);
				}

			}
			else {
				//std::cout << "Skipped Face " << std::endl;
				continue;
			}			
		}

		//render_vtk_triangular_cell(facePolys);

		//std::cout << "Polys Size: " << facePolys.size() << std::endl;
		//vtkNew<vtkAppendPolyData> appendFilter;
		//for (int pp = 0; pp < facePolys.size(); pp++) {
		//	appendFilter->AddInputData(facePolys[pp]);
		//}
		//appendFilter->Update();

		//polys.push_back(appendFilter->GetOutput());
		
	} while (cla.inc());

	//std::cout << 

	//render_vtk_face(globalVertices, globalFaces, "FinalMesh");
	std::cout << " Creating Mesh " << std::endl;

	Eigen::MatrixXd meshVertices(globalVertices.size(), 3);
	for (int vrow{ 0 }; vrow < globalVertices.size(); vrow++) {
		Eigen::Vector3d row = { globalVertices[vrow][0], globalVertices[vrow][1], globalVertices[vrow][2] };
		meshVertices.row(vrow) = row;
	}

	vtkNew<vtkNamedColors> colors;

	// Create points
	vtkNew<vtkPoints> points;

	if (meshVertices.cols() == 3) {
		for (int i = 0; i < meshVertices.rows(); ++i) {
			points->InsertNextPoint(meshVertices(i, 0), meshVertices(i, 1), meshVertices(i, 2));
		}
	}
	else if (meshVertices.cols() == 2) {
		for (int i = 0; i < meshVertices.rows(); ++i) {
			points->InsertNextPoint(meshVertices(i, 0), meshVertices(i, 1), 0.0);
		}
	}

	vtkNew<vtkCellArray> cells;
	for (int idx{ 0 }; idx < globalFaces.size(); idx++) {
		cells->InsertNextCell(3);
		cells->InsertCellPoint(globalFaces[idx][0]);
		cells->InsertCellPoint(globalFaces[idx][1]);
		cells->InsertCellPoint(globalFaces[idx][2]);
	}

	// Create polydata and set points
	vtkNew<vtkPolyData> polyData;
	polyData->SetPoints(points);
	polyData->SetPolys(cells);

	// decide mesh resolution
	int dim[3] = {};
	double thickness{ 0.3 };

	if (thickness < 0.3) {
		dim[0] = 100;
		dim[1] = 100;
		dim[2] = 100;
	}
	else if (0.3 <= thickness && thickness < 0.5) {
		dim[0] = 50;
		dim[1] = 50;
		dim[2] = 50;
	}
	else {
		dim[0] = 100;
		dim[1] = 100;
		dim[2] = 100;
	}

	// build along line
	vtkNew<vtkImplicitModeller> implictModeller;
	implictModeller->AddInputData(polyData);
	implictModeller->SetMaximumDistance(thickness / 2 * 1.01);
	implictModeller->SetSampleDimensions(dim[0], dim[1], dim[2]);
	implictModeller->SetModelBounds(
		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

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
	std::string fileName = "test.stl";
	stlWriter->SetFileName(fileName.c_str());
	//stlWriter->SetInputConnection(sincSmoother->GetOutputPort());
	stlWriter->SetInputConnection(norms->GetOutputPort());
	stlWriter->Write();

	vtkSmartPointer<vtkPolyData> pd;
	pd = norms->GetOutput();
	render_vtk_polydata(pd);

	//for (auto face : globalFaces) {
	//	std::cout << face[0] << " " << face[1] << " " << face[2] << std::endl;
	//}

	//render_vtk_points(meshVertices, "FinalPoints");
	//render_vtk_face(meshVertices, globalFaces, "FinalPoints");

	// create scaffold

	//render_vtk_triangular_cell(polys);
};

