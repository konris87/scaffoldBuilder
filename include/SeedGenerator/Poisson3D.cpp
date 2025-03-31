#include <Eigen/Dense>
#include "Poisson3D.h"
#include "Utils/Utils.h"
#include <vtkPolyDataNormals.h>
#include <vtkTriangle.h>
#include <vtkSTLReader.h>
#include <vtkKdTreePointLocator.h>
#include "Wall/Wall.h"

// base class implementation

void Poisson3DSampler::generate_mesh(
	const double& thickness,
	const std::string & fileName) {

	vtkNew<vtkNamedColors> colors;

	// decide mesh resolution
	int dim[3] = {};

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

	std::cout << "Dim: " << dim[0] << std::endl;

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

void Poisson3DSampler::generate_voro() {

		std::cout << "Seed Number: " << seeds.size() << std::endl;

		con = new voro::container(
			xMin, xMax, yMin, yMax, zMin, zMax, nX, nY, nZ, false, false, false, 16);

		std::vector<vtkSmartPointer<vtkPolyData>> polys;

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

			//seeds[seedId] = centroid;

			vtkSmartPointer<vtkPolyData> poly = cell_2_vtk(cellNeighs, cellVertices, faceVertices);
			polys.push_back(poly);

		} while (cla.inc());

		vtkNew<vtkAppendPolyData> appendFilter;
		for (int i = 0; i < polys.size(); i++) {
			appendFilter->AddInputData(polys[i]);
		}
		appendFilter->Update();
		//scaffoldMesh = vtkSmartPointer<vtkPolyData>::New();
		scaffoldMesh = appendFilter->GetOutput();
};


// --------------------------------------------------------------------------
// Poisson 3D inside box

Poisson3D::Poisson3D(
	const float& rMinVal,
	const float& rMaxVal,
	const std::array<double, 3>& rootVal,
	const std::array<float, 6>& bounds) : Poisson3DSampler(rMinVal, rMaxVal, rootVal) {

	xMin = bounds[0];
	xMax = bounds[1];
	yMin = bounds[2];
	yMax = bounds[3];
	zMin = bounds[4];
	zMax = bounds[5];
	
	scale = sqrt(pow(xMax - xMin, 2) + pow(yMax - yMin, 2) + pow(zMax - zMin, 2));

	cellSize = rMin / sqrt(3.0f);
	std::cout << "cell size: " << cellSize << std::endl;
};

void Poisson3D::generate_seeds() {

	try {

		float radius = rMin;

		int currSize = { 0 };
		if (!seeds.empty()) {
			currSize = seeds.size();
		}
		// this is the active list
		std::vector<int> active;

		// initialize with root and update cell
		seeds.push_back(root);

		int N{ 0 };
		active.push_back(N);

		// the centroid is in the central cell
		std::array<int, 3> centralIdx{ 0, 0, 0 };
		grid[centralIdx].seedIdx = N;

		// finally push in the neighbors the index of the seed! The number of neighbors is determined
		// based on the corresponding radius
		double rxi = radius;
		int n = ceil(rxi / radius) + 1;

		pushIdxs(n, centralIdx, N);

		while (!active.empty()) {

			// step 1: chose randomly an index from the active list
			int randomIndex = rand() % active.size();
			int idx = active[randomIndex];

			// step 2: identify the poisson disc parameter for the point,
			std::array<double, 3> xi = seeds[idx + currSize];

			// identify the cell of the grid that the point is located
			std::array<int, 3> cellIdx = getGridIndex(xi, cellSize);

			// step 3: generate k uniformly random points inside the spherical annulus
			std::random_device rd;
			std::mt19937 gen(rd());

			std::uniform_real_distribution<> disx(0.0, 1.0);
			std::uniform_real_distribution<> disphi(0.0, 2 * PI);
			std::uniform_real_distribution<> distheta(-1.0, 1.0);

			std::vector<std::array<double, 3>> kPts;

			//std::cout << "step 3: creating k neighbors" << std::endl;
			for (int i{ 0 }; i < neighNr; i++) {

				double r = cbrt((disx(gen) * pow(2 * rxi, 3))) + ((1 - disx(gen) * pow(rxi, 3)));
				double phi = disphi(gen);
				double theta = acos(distheta(gen));

				double x = r * sin(theta) * cos(phi);
				double y = r * sin(theta) * sin(phi);
				double z = r * cos(theta);

				// add a point with center the xi point
				kPts.push_back({ x + xi[0], y + xi[1], z + xi[2] });
			}

			// a boolean flag to check if all candidates are false
			bool valid = false;

			// for each candidate point get the index and check if it is inside a distance r of existing samples
			for (const auto& pt : kPts) {
				Eigen::Vector3d ept{ pt[0], pt[1], pt[2] };

				// get the index of the point under consideration
				std::array<int, 3> ptIdx = getGridIndex(pt, cellSize);

				// check if inside
				if (pt[0] < xMin || pt[0] > xMax ||
					pt[1] < yMin || pt[1] > yMax ||
					pt[2] < zMin || pt[2] > zMax) {
					continue;
				}

				// also check if a cell contains already a pt
				if (grid[ptIdx].seedIdx > 0) {
					//std::cout << "Already contains a point" << std::endl;
					continue;
					//std::cout << "Point inside Mesh!" << std::endl;
				}

				// loop through the cellIdxs of the cell
				std::vector<double> dists;

				if (grid[ptIdx].cellIdxs.empty()) {
					dists.push_back(distance(pt, xi));
				}
				else {
					for (const auto& nIdx : grid[ptIdx].cellIdxs) {
						double dist = distance(pt, seeds[nIdx + currSize]);
						// get the neighbors and find the distance
						dists.push_back(dist);
					}
				}

				double minDist = *min_element(dists.begin(), dists.end());

				// if all neighbors are valid push the candidate to the seeds
				// this means that the distance is larger than the radius of the 
				// point under consideration
				if (minDist > radius) {
					seeds.push_back(pt);
					active.push_back(N + 1);
					int newN = ceil(radius / rMin) + 1;
					pushIdxs(newN, ptIdx, N + 1);
					grid[ptIdx].seedIdx = N + 1;
					++N;
					valid = true;
				}
			}

			// if none is valid then erase i from the active list
			if (!valid) {
				//std::cout << "None of the points is valid!" << std::endl;
				//std::cout << "Removing idx " << idx << " from active list. " << std::endl;
				active.erase(std::remove(active.begin(), active.end(), idx), active.end());
			}
			//std::cout << "Active list length: " << active.size() << std::endl;
			//std::cout << "------------------------" << std::endl;
		}
	}
	catch (const std::exception& e) {
		std::cout << e.what() << std::endl;
	}

};

void Poisson3D::generate_seeds(DistanceEstimator& distEstimator, RadiusFunction& radFunc) {

	try {

		//std::cout << root[0] << " " << root[1] << " " << root[2] << std::endl;

		int currSize = { 0 };

		// this is the active list
		std::vector<int> active;

		seeds.push_back(root);

		int N{ 0 };

		active.push_back(N);

		// the centroid is in the central cell
		std::array<int, 3> centralIdx{ 0, 0, 0 };
		grid[centralIdx].seedIdx = N;

		// finally push in the neighbors the index of the seed! The number of neighbors is determined
		// based on the corresponding radius
		double rxi = rMin;
		int n = ceil(rxi / rMin) + 1;

		pushIdxs(n, centralIdx, N);

		radii.push_back(rxi);

		while (!active.empty()) {

			// step 1: chose randomly an index from the active list
			int randomIndex = rand() % active.size();
			int idx = active[randomIndex];

			// step 2: identify the poisson disc parameter for the point,
			std::array<double, 3> xi = seeds[idx + currSize];

			// identify the cell of the grid that the point is located
			std::array<int, 3> cellIdx = getGridIndex(xi, cellSize);

			// estimate also the radius based on the distance from the wall
			double d = distEstimator.compute_distance({ xi[0], xi[1], xi[2] });
			//std::cout << "distance from surface: " << d << std::endl;

			double rxi = radFunc.estimate_radius(d, rMin, rMax);
			//std::cout << "distance: " << abs(d) << " radius: " << rxi << std::endl;

			// step 3: generate k uniformly random points inside the spherical annulus
			std::random_device rd;
			std::mt19937 gen(rd());

			std::uniform_real_distribution<> disx(0.0, 1.0);
			std::uniform_real_distribution<> disphi(0.0, 2 * PI);
			std::uniform_real_distribution<> distheta(-1.0, 1.0);

			std::vector<std::array<double, 3>> kPts;

			for (int i{ 0 }; i < neighNr; i++) {

				double r = cbrt((disx(gen) * pow(2 * rxi, 3))) + ((1 - disx(gen) * pow(rxi, 3)));

				double phi = disphi(gen);
				double theta = acos(distheta(gen));

				double x = r * sin(theta) * cos(phi);
				double y = r * sin(theta) * sin(phi);
				double z = r * cos(theta);

				// add a point with center the xi point
				kPts.push_back({ x + xi[0], y + xi[1], z + xi[2] });
			}

			// a boolean flag to check if all candidates are false
			bool valid = false;

			// for each candidate point get the index and check if it is inside a distance r of existing samples
			for (const auto& pt : kPts) {
				// get the index of the point under consideration
				std::array<int, 3> ptIdx = getGridIndex(pt, cellSize);
				//std::cout << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
				if (pt[0] < xMin || pt[0] > xMax ||
					pt[1] < yMin || pt[1] > yMax ||
					pt[2] < zMin || pt[2] > zMax) {
					continue;
				}

				// also check if a cell contains already a pt
				if (grid[ptIdx].seedIdx > 0) {
					continue;
				}

				//std::cout << "check distance" << std::endl;
				double dFromSurf = distEstimator.compute_distance({ pt[0], pt[1], pt[2] });

				//std::cout << "get radius" << std::endl;
				double ryi = radFunc.estimate_radius(dFromSurf, rMin, rMax);

				// loop through the cellIdxs of the cell
				std::vector<double> dists;

				if (grid[ptIdx].cellIdxs.empty()) {
					dists.push_back(distance(pt, xi));
				}
				else {
					for (const auto& nIdx : grid[ptIdx].cellIdxs) {
						//std::cout << nIdx << " " << std::endl;
						double dist = distance(pt, seeds[nIdx + currSize]);
						// get the neighbors and find the distance
						dists.push_back(dist);
					}
				}
				//std::cout << "find min dist" << std::endl;
				//if (dists.empty()) {
				//	std::cout << "empty" << std::endl;
				//}

				double minDist = *min_element(dists.begin(), dists.end());

				// if all neighbors are valid push the candidate to the seeds
				// this means that the distance is larger than the radius of the 
				// point under consideration
				//std::cout << "check distance tol" << std::endl;
				if (minDist > ryi) {
					N++;
					seeds.push_back(pt);
					radii.push_back(ryi);
					active.push_back(N);
					int newN = ceil(ryi / rMin) + 1;
					pushIdxs(newN, ptIdx, N);
					grid[ptIdx].seedIdx = N;
					valid = true;
				}
			}
			//std::cout << "Active List: ";
			//for (const auto& t : active) {
			//    std::cout << t << " ";
			//}
			// if none is valid then erase i from the active list
			if (!valid) {
				//std::cout << "None of the points is valid!" << std::endl;
				//std::cout << "Removing idx " << idx << " from active list. " << std::endl;
				active.erase(std::remove(active.begin(), active.end(), idx), active.end());
				//std::cout << "Active List After Removing: ";
				//for (const auto& t : active) {
				//    std::cout << t << " ";
				//}

			}
			//std::cout << "Active list length: " << active.size() << std::endl;
			//std::cout << "------------------------" << std::endl;
		}
	}
	catch (const std::exception& e) {
		std::cout << e.what() << std::endl;
	}
};

// ---------------------------------------------------------------------------------
// Poisson 3D wall

Poisson3DWall::Poisson3DWall(
	vtkSmartPointer<vtkPolyData>& containerMesh,
	const float& rMinVal,
	const float& rMaxVal,
	const int& neighbors,
	const int& nXval, const int& nYval, const int& nZval
	) : neighbors(neighbors){

	container = vtkSmartPointer<vtkPolyData > ::New();
	container = containerMesh;

	// get bounds
	double bounds[6];
	container->GetBounds(bounds);
	xMin = bounds[0];
	xMax = bounds[1];
	yMin = bounds[2];
	yMax = bounds[3];
	zMin = bounds[4];
	zMax = bounds[5];

	nX = nXval;
	nY = nYval;
	nZ = nZval;

	rMin = rMinVal;
	rMax = rMaxVal;

	double center[3];
	container->GetCenter(center);

	root[0] = center[0];
	root[1] = center[1];
	root[2] = center[2];

	scale = sqrt(pow(xMax - xMin, 2) + pow(yMax - yMin, 2) + pow(zMax - zMin, 2));

	cellSize = rMin / sqrt(3);

	//_populate_with_barycenters();
	
	std::cout << "seed Nr: " << seeds.size();
};

void Poisson3DWall::generate_seeds() {

	try {

		float radius = rMin;

		int currSize = { 0 };
		if (!seeds.empty()) {
			currSize = seeds.size();
		}
		// this is the active list
		std::vector<int> active;

		// initialize with root and update cell
		seeds.push_back(root);

		int N{ 0 };
		active.push_back(N);

		// the centroid is in the central cell
		std::array<int, 3> centralIdx{ 0, 0, 0 };
		grid[centralIdx].seedIdx = N;

		// finally push in the neighbors the index of the seed! The number of neighbors is determined
		// based on the corresponding radius
		double rxi = radius;
		int n = ceil(rxi / radius) + 1;

		pushIdxs(n, centralIdx, N);

		while (!active.empty()) {

			// step 1: chose randomly an index from the active list
			int randomIndex = rand() % active.size();
			int idx = active[randomIndex];

			// step 2: identify the poisson disc parameter for the point,
			std::array<double, 3> xi = seeds[idx + currSize];

			// identify the cell of the grid that the point is located
			std::array<int, 3> cellIdx = getGridIndex(xi, cellSize);

			// step 3: generate k uniformly random points inside the spherical annulus
			std::random_device rd;
			std::mt19937 gen(rd());

			std::uniform_real_distribution<> disx(0.0, 1.0);
			std::uniform_real_distribution<> disphi(0.0, 2 * PI);
			std::uniform_real_distribution<> distheta(-1.0, 1.0);

			std::vector<std::array<double, 3>> kPts;

			//std::cout << "step 3: creating k neighbors" << std::endl;
			for (int i{ 0 }; i < neighNr; i++) {

				double r = cbrt((disx(gen) * pow(2 * rxi, 3))) + ((1 - disx(gen) * pow(rxi, 3)));
				double phi = disphi(gen);
				double theta = acos(distheta(gen));

				double x = r * sin(theta) * cos(phi);
				double y = r * sin(theta) * sin(phi);
				double z = r * cos(theta);

				// add a point with center the xi point
				kPts.push_back({ x + xi[0], y + xi[1], z + xi[2] });
			}

			// a boolean flag to check if all candidates are false
			bool valid = false;

			// for each candidate point get the index and check if it is inside a distance r of existing samples
			for (const auto& pt : kPts) {
				Eigen::Vector3d ept{ pt[0], pt[1], pt[2] };

				// get the index of the point under consideration
				std::array<int, 3> ptIdx = getGridIndex(pt, cellSize);

				// first check if is inside
				if (!is_inside_mesh(container, ept)) {
					//std::cout << "Point outside of Mesh!" << std::endl;
					continue;
				}

				// also check if a cell contains already a pt
				if (grid[ptIdx].seedIdx > 0) {
					//std::cout << "Already contains a point" << std::endl;
					continue;
					//std::cout << "Point inside Mesh!" << std::endl;
				}

				// loop through the cellIdxs of the cell
				std::vector<double> dists;

				if (grid[ptIdx].cellIdxs.empty()) {
					dists.push_back(distance(pt, xi));
				}
				else {
					for (const auto& nIdx : grid[ptIdx].cellIdxs) {
						double dist = distance(pt, seeds[nIdx + currSize]);
						// get the neighbors and find the distance
						dists.push_back(dist);
					}
				}

				double minDist = *min_element(dists.begin(), dists.end());

				// if all neighbors are valid push the candidate to the seeds
				// this means that the distance is larger than the radius of the 
				// point under consideration
				if (minDist > radius) {
					seeds.push_back(pt);
					active.push_back(N + 1);
					int newN = ceil(radius / rMin) + 1;
					pushIdxs(newN, ptIdx, N + 1);
					grid[ptIdx].seedIdx = N + 1;
					++N;
					valid = true;
				}
			}

			// if none is valid then erase i from the active list
			if (!valid) {
				//std::cout << "None of the points is valid!" << std::endl;
				//std::cout << "Removing idx " << idx << " from active list. " << std::endl;
				active.erase(std::remove(active.begin(), active.end(), idx), active.end());
			}
			//std::cout << "Active list length: " << active.size() << std::endl;
			//std::cout << "------------------------" << std::endl;
		}
	}
	catch (const std::exception& e) {
		std::cout << e.what() << std::endl;
	}

};

void Poisson3DWall::generate_seeds(DistanceEstimator& distEstimator, RadiusFunction& radFunc) {

	try {

		int currSize = { 0 };
		if (!seeds.empty()) {
			currSize = seeds.size();
		}
		// this is the active list
		std::vector<int> active;

		// initialize with root and update cell
		seeds.push_back(root);

		int N{ 0 };
		active.push_back(N);

		// the centroid is in the central cell
		std::array<int, 3> centralIdx{ 0, 0, 0 };
		grid[centralIdx].seedIdx = N;

		// finally push in the neighbors the index of the seed! The number of neighbors is determined
		// based on the corresponding radius
		double rxi = rMin;
		int n = ceil(rxi / rMin) + 1;

		pushIdxs(n, centralIdx, N);

		radii.push_back(rxi);

		while (!active.empty()) {

			// step 1: chose randomly an index from the active list
			int randomIndex = rand() % active.size();
			int idx = active[randomIndex];

			// step 2: identify the poisson disc parameter for the point,
			std::array<double, 3> xi = seeds[idx + currSize];
			
			// identify the cell of the grid that the point is located
			std::array<int, 3> cellIdx = getGridIndex(xi, cellSize);

			//std::cout << "distance from mesh face" << std::endl;
			double d = distEstimator.compute_distance({ xi[0], xi[1], xi[2] });
			double rxi = radFunc.estimate_radius(abs(d), rMin, rMax);
			//double rxi = quad_radius_func(abs(d), rMax, rMin, 5);
			//std::cout << "distance: " << abs(d) << " radius: " << rxi << std::endl;
		
			// step 3: generate k uniformly random points inside the spherical annulus
			std::random_device rd;
			std::mt19937 gen(rd());

			std::uniform_real_distribution<> disx(0.0, 1.0);
			std::uniform_real_distribution<> disphi(0.0, 2 * PI);
			std::uniform_real_distribution<> distheta(-1.0, 1.0);

			std::vector<std::array<double, 3>> kPts;

			//std::cout << "step 3: creating k neighbors" << std::endl;
			for (int i{ 0 }; i < neighNr; i++) {

				double r = cbrt((disx(gen) * pow(2 * rxi, 3))) + ((1 - disx(gen) * pow(rxi, 3)));
				double phi = disphi(gen);
				double theta = acos(distheta(gen));

				double x = r * sin(theta) * cos(phi);
				double y = r * sin(theta) * sin(phi);
				double z = r * cos(theta);

				// add a point with center the xi point
				kPts.push_back({ x + xi[0], y + xi[1], z + xi[2] });
			}

			// a boolean flag to check if all candidates are false
			bool valid = false;

			// for each candidate point get the index and check if it is inside a distance r of existing samples
			for (const auto& pt : kPts) {
				Eigen::Vector3d ept{ pt[0], pt[1], pt[2] };

				// get the index of the point under consideration
				std::array<int, 3> ptIdx = getGridIndex(pt, cellSize);

				// first check if is inside
				if (!is_inside_mesh(container, ept)) {
					//std::cout << "Point outside of Mesh!" << std::endl;
					continue;
				}

				// also check if a cell contains already a pt
				if (grid[ptIdx].seedIdx > 0) {
					//std::cout << "Already contains a point" << std::endl;
					continue;
					//std::cout << "Point inside Mesh!" << std::endl;
				}

				double dFromSurf = distEstimator.compute_distance({ pt[0], pt[1], pt[2] });
				double ryi = radFunc.estimate_radius(dFromSurf, rMin, rMax);

				// loop through the cellIdxs of the cell
				std::vector<double> dists;

				if (grid[ptIdx].cellIdxs.empty()) {
					dists.push_back(distance(pt, xi));
				}
				else {
					for (const auto& nIdx : grid[ptIdx].cellIdxs) {
						double dist = distance(pt, seeds[nIdx + currSize]);
						// get the neighbors and find the distance
						dists.push_back(dist);
					}
				}

				double minDist = *min_element(dists.begin(), dists.end());

				// if all neighbors are valid push the candidate to the seeds
				// this means that the distance is larger than the radius of the 
				// point under consideration
				if (minDist > ryi) {
					seeds.push_back(pt);
					radii.push_back(ryi);
					active.push_back(N + 1);
					int newN = ceil(ryi / rMin) + 1;
					pushIdxs(newN, ptIdx, N + 1);
					grid[ptIdx].seedIdx = N + 1;
					++N;
					valid = true;
				}
			}
			
			// if none is valid then erase i from the active list
			if (!valid) {
				//std::cout << "None of the points is valid!" << std::endl;
				//std::cout << "Removing idx " << idx << " from active list. " << std::endl;
				active.erase(std::remove(active.begin(), active.end(), idx), active.end());
			}
			//std::cout << "Active list length: " << active.size() << std::endl;
			//std::cout << "------------------------" << std::endl;
		}
	}
	catch (const std::exception& e) {
		std::cout << e.what() << std::endl;
	}

};

void Poisson3DWall::_populate_with_barycenters() {

	for (vtkIdType i = 0; i < container->GetNumberOfCells(); ++i) {
		vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(container->GetCell(i));
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

void Poisson3DWall::generate_voro() {

	// create a kdtree to quickly find close points
	vtkSmartPointer<vtkKdTreePointLocator> kdTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
	vtkSmartPointer<vtkPoints> bcvtk = vtkSmartPointer<vtkPoints>::New();
	for (int i{ 0 }; i < bCenters.size(); i++) {
		bcvtk->InsertNextPoint(bCenters[i][0], bCenters[i][1], bCenters[i][2]);
	}

	vtkNew<vtkPolyData> bcvtkdata;
	bcvtkdata->SetPoints(bcvtk);
	kdTree->SetDataSet(bcvtkdata);
	kdTree->BuildLocator();

	// create container
	con = new voro::container(
		xMin, xMax, yMin, yMax, zMin, zMax, nX, nY, nZ, false, false, false, 16);

	for (int i{ 0 }; i < seeds.size(); i++) {
		con->put(i, seeds[i][0], seeds[i][1], seeds[i][2]);
	}

	MeshWall wall(container, normals, kdTree, neighbors);
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