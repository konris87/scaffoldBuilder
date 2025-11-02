#include "SeedGenerator.h"
#include <random>

void Random::run(const ContainerAdapter& adapter, std::vector<std::array<double, 3>>& seeds, const RunConfig&) {

	seeds.clear();
	seeds.reserve(seedNr);

	auto bounds = adapter.get_bounds();

	std::cout <<
		"xMax: " << bounds.xMin << " xMax: " << bounds.xMax <<
		"yMin: " << bounds.yMin << " yMax: " << bounds.yMax <<
		"zMin: " << bounds.zMin << " zMax: " << bounds.zMax << std::endl;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> disX(bounds.xMin + 0.1, bounds.xMax - 0.1);
	std::uniform_real_distribution<> disY(bounds.yMin + 0.1, bounds.yMax - 0.1);
	std::uniform_real_distribution<> disZ(bounds.zMin + 0.1, bounds.zMax - 0.1);

	while (seeds.size() < seedNr) {
		Pt pt{ disX(gen), disY(gen), disZ(gen) };
		if (adapter.is_inside(pt)) {
			seeds.push_back(pt);
		}
	}

};

void Poisson3D::run(const ContainerAdapter& adapter, std::vector<std::array<double, 3>>& seeds, const RunConfig& cfg) {

	seeds.clear();
	radii.clear();
	grid.clear();

	Bounds bounds = adapter.get_bounds();

	double scale = std::sqrt(
		std::pow(bounds.xMax - bounds.xMin, 2) +
		std::pow(bounds.yMax - bounds.yMin, 2) +
		std::pow(bounds.zMax - bounds.zMin, 2));

	double cellSize = rMin / sqrt(3.0f);

	// check if there are values in the RunConfig struct, if yes we have varied poisson
	const bool varied = (bool)cfg.dist && (bool)cfg.rad;

	try {

		int currSize = { 0 };

		// this is the active list
		std::vector<int> active;

		seeds.push_back(root);

		int N{ 0 };

		active.push_back(N);

		// the centroid is in the central cell
		std::array<int, 3> centralIdx{ 0, 0, 0 };
		grid[centralIdx].seedIdx = N;

		// finally push in the neighbors the index of the seed! 
		// The number of neighbors is determined
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
			const double d = varied ? cfg.dist->compute_distance({ xi[0], xi[1], xi[2] }) : 0.0;
			
			// if varied get the value from the function otherwise rMin = rMax
			const double rxi = varied ? cfg.rad->estimate_radius(d, rMin, rMax) : rMin;

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
				//if (pt[0] < xMin || pt[0] > xMax ||
				//	pt[1] < yMin || pt[1] > yMax ||
				//	pt[2] < zMin || pt[2] > zMax) {
				//	continue;
				//}

				if (!adapter.is_inside(pt)) {
					continue;
				}

				// also check if a cell contains already a pt
				if (grid[ptIdx].seedIdx > 0) {
					continue;
				}

				//std::cout << "check distance" << std::endl;
				double dFromSurf = varied ? cfg.dist->compute_distance({ pt[0], pt[1], pt[2] }) : 0.0f;

				//std::cout << "get radius" << std::endl;
				double ryi = varied ? cfg.rad->estimate_radius(dFromSurf, rMin, rMax) : rMin;

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



