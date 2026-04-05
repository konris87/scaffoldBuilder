#include "SeedGenerator.h"
#include <random>

void Random::run(const IContainer& adapter) {

	seeds.clear();
	seeds.reserve(seedNr);

	auto bounds = adapter.compute_bounds();

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> disX(bounds.xMin + 0.1, bounds.xMax - 0.1);
	std::uniform_real_distribution<> disY(bounds.yMin + 0.1, bounds.yMax - 0.1);
	std::uniform_real_distribution<> disZ(bounds.zMin + 0.1, bounds.zMax - 0.1);

	while (seeds.size() < seedNr) {
		Vec3 pt{ (float)disX(gen), (float)disY(gen), (float)disZ(gen) };
		if (adapter.is_inside(pt)) {
			seeds.push_back(pt);
		}
	}

	if (!seeds.empty()) {
		model = std::make_unique<VisualizeSeeds>(seeds);
	}

	version++;
};

void Random::render_gui() {
	
	std::string nr = std::to_string(seeds.size());
	ImGui::Text(("Seeds: " + nr).c_str());

	ImGui::Separator();
	ImGui::SetNextItemWidth(200);
	ImGui::InputInt("Seed Number", &seedNr);
};

void Poisson3D::run(const IContainer& adapter, const RunConfig& cfg) {

	config = cfg;

	seeds.clear();
	radii.clear();
	grid.clear();

	// define root

	Bounds bounds = adapter.compute_bounds();

	root = bounds.center;

	// if the root is not inside randomly create it 
	if (!adapter.is_inside(root)) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> disX(bounds.xMin + 0.1, bounds.xMax - 0.1);
		std::uniform_real_distribution<> disY(bounds.yMin + 0.1, bounds.yMax - 0.1);
		std::uniform_real_distribution<> disZ(bounds.zMin + 0.1, bounds.zMax - 0.1);

		while (!adapter.is_inside(root)) {
			root.x = (float)disX(gen);
			root.y = (float)disY(gen);
			root.z = (float)disZ(gen);
		}
	}

	double scale = std::sqrt(
		std::pow(bounds.xMax - bounds.xMin, 2) +
		std::pow(bounds.yMax - bounds.yMin, 2) +
		std::pow(bounds.zMax - bounds.zMin, 2));

	double cellSize = rMin / sqrt(3.0f);

	// check if there are values in the RunConfig struct, if yes we have varied poisson
	//const bool varied = (bool)cfg.dist && (bool)cfg.rad;

	// lambda function to get squared distance
	auto distance = [&](const Vec3& p1, const Vec3& p2) {
		return (p2 - p1).norm();
	};

	// generate k uniformly random points inside the spherical annulus
	std::random_device rd;
	std::mt19937 gen(rd());

	std::uniform_real_distribution<> disx(0.0, 1.0);
	std::uniform_real_distribution<> disphi(0.0, 2 * PI);
	std::uniform_real_distribution<> distheta(-1.0, 1.0);

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
			Vec3 xi = seeds[idx + currSize];

			// identify the cell of the grid that the point is located
			std::array<int, 3> cellIdx = getGridIndex(xi, cellSize);

			// estimate also the radius based on the distance from the metric,
			// since it will be inside the distance will be negative, thus,
			// we need the absolute value
			double rxi = 0.0f;

			if (cfg.dist) {
				const float d = std::abs(cfg.dist->compute_distance(xi));
				rxi = cfg.rad->estimate_radius(d, rMin, rMax);
			}
			else {
				// this is for the uniform random between rMin rMax
				rxi = cfg.rad->estimate_radius(0, rMin, rMax);
			}
			
			// if varied get the value from the function otherwise rMin = rMax

			std::vector<Vec3> kPts;

			for (int i{ 0 }; i < neighNr; i++) {

				float r = cbrt((disx(gen) * pow(2 * rxi, 3))) + ((1 - disx(gen) * pow(rxi, 3)));

				float phi = disphi(gen);
				float theta = acos(distheta(gen));

				float x = r * sin(theta) * cos(phi);
				float y = r * sin(theta) * sin(phi);
				float z = r * cos(theta);

				// add a point with center the xi point
				kPts.push_back(Vec3(x + xi.x, y + xi.y, z + xi.z ));
			}

			// a boolean flag to check if all candidates are false
			bool valid = false;

			// for each candidate point get the index and check if it is inside a distance r of existing samples
			for (const auto& pt : kPts) {
				// get the index of the point under consideration
				std::array<int, 3> ptIdx = getGridIndex(pt, cellSize);

				if (!adapter.is_inside(pt)) {
					continue;
				}

				// also check if a cell contains already a pt
				if (grid[ptIdx].seedIdx > 0) {
					continue;
				}


				double ryi = 0.0f;

				if (cfg.dist) {
					//std::cout << "check distance" << std::endl;
					double dFromSurf = std::abs(cfg.dist->compute_distance(pt));
					//std::cout << "get radius" << std::endl;
					ryi = cfg.rad->estimate_radius(dFromSurf, rMin, rMax);
				}
				else {
					ryi = cfg.rad->estimate_radius(0.0, rMin, rMax);

				}

				////std::cout << "check distance" << std::endl;
				//double dFromSurf = std::abs(cfg.dist->compute_distance(pt));

				////std::cout << "get radius" << std::endl;
				//double ryi = cfg.rad->estimate_radius(dFromSurf, rMin, rMax);

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

	if (!seeds.empty()) {
		model = std::make_unique<VisualizeSeeds>(seeds);
	}
	version++;
};

void Poisson3D::run(const IContainer& adapter) {

	seeds.clear();
	radii.clear();
	grid.clear();

	// define root

	Bounds bounds = adapter.compute_bounds();

	root = bounds.center;

	// if the root is not inside randomly create it 
	if (!adapter.is_inside(root)) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> disX(bounds.xMin + 0.1, bounds.xMax - 0.1);
		std::uniform_real_distribution<> disY(bounds.yMin + 0.1, bounds.yMax - 0.1);
		std::uniform_real_distribution<> disZ(bounds.zMin + 0.1, bounds.zMax - 0.1);

		while (!adapter.is_inside(root)) {
			root.x = (float)disX(gen);
			root.y = (float)disY(gen);
			root.z = (float)disZ(gen);
		}
	}

	double scale = std::sqrt(
		std::pow(bounds.xMax - bounds.xMin, 2) +
		std::pow(bounds.yMax - bounds.yMin, 2) +
		std::pow(bounds.zMax - bounds.zMin, 2));

	double cellSize = rMin / sqrt(3.0f);

	// lambda function to get squared distance
	auto distance = [&](const Vec3& p1, const Vec3& p2) {
		return (p2 - p1).norm();
	};

	// generate k uniformly random points inside the spherical annulus
	std::random_device rd;
	std::mt19937 gen(rd());

	std::uniform_real_distribution<> disx(0.0, 1.0);
	std::uniform_real_distribution<> disphi(0.0, 2 * PI);
	std::uniform_real_distribution<> distheta(-1.0, 1.0);

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
			Vec3 xi = seeds[idx + currSize];

			// identify the cell of the grid that the point is located
			std::array<int, 3> cellIdx = getGridIndex(xi, cellSize);

			// estimate also the radius based on the distance from the wall
			const double d = 0.0;
			const double rxi = rMin;

			std::vector<Vec3> kPts;

			for (int i{ 0 }; i < neighNr; i++) {

				float r = cbrt((disx(gen) * pow(2 * rxi, 3))) + ((1 - disx(gen) * pow(rxi, 3)));

				float phi = disphi(gen);
				float theta = acos(distheta(gen));

				float x = r * sin(theta) * cos(phi);
				float y = r * sin(theta) * sin(phi);
				float z = r * cos(theta);

				Vec3 newP = { x + xi.x, y + xi.y, z + xi.z };

				// add a point with center the xi point
				kPts.push_back(newP);
			}

			// a boolean flag to check if all candidates are false
			bool valid = false;

			// for each candidate point get the index and check if it is inside a distance r of existing samples
			for (const auto& pt : kPts) {
				// get the index of the point under consideration
				std::array<int, 3> ptIdx = getGridIndex(pt, cellSize);

				if (!adapter.is_inside(pt)) {
					continue;
				}

				// also check if a cell contains already a pt
				if (grid[ptIdx].seedIdx > 0) {
					continue;
				}

				//std::cout << "check distance" << std::endl;
				double dFromSurf = 0.0f;

				//std::cout << "get radius" << std::endl;
				double ryi = rMin;

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

	if (!seeds.empty()) {
		model = std::make_unique<VisualizeSeeds>(seeds);
	}
	version++;
};

void Poisson3D::render_gui() {

	if (type == ObjectType::UniformGeneratorType) {

		std::string seedNr = std::to_string(seeds.size());
		ImGui::Text(("Seeds: " + seedNr).c_str());

		ImGui::Separator();
		ImGui::SetNextItemWidth(200);
		ImGui::InputDouble("Radius", &rMin, 0.001, 1.0, "%.3f");
	}

	else if (type == ObjectType::VariedGeneratorType) {

		std::string seedNr = std::to_string(seeds.size());
		ImGui::Text(("Seeds: " + seedNr).c_str());

		ImGui::Separator();

		ImGui::SetNextItemWidth(200);
		ImGui::InputDouble("Minimum Radius", &rMin, 0.001, 1.0, "%.3f");
		ImGui::SetNextItemWidth(200);
		ImGui::InputDouble("Maximum Radius", &rMax, 0.001, 1.0, "%.3f");

		ImGui::SeparatorText("");

		ImGui::SeparatorText("");
	}
};

void Poisson3D::set_min_radius(const double newRadius) { rMin = newRadius; };

void Poisson3D::set_max_radius(const double newRadius) { rMax = newRadius; };

void Poisson3D::set_center(const Vec3& newCenter) { root = newCenter; };
