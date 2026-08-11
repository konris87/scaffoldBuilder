#include "SeedGenerator.h"
#include "Misc/Imgui_Stdlib.h"
#include <random>

// =============================================================================
// Random Generator Class Implementation
// =============================================================================
void Random::run(const IContainer& adapter) {

	containerType = adapter.get_type();
	seeds.clear();
	seeds.reserve(seedNr);

	bounds = adapter.compute_bounds();

	// rngSeed == 0 -> a fresh random_device draw (a new realization every run);
	// rngSeed != 0 -> reproducible. See InterfaceSeedGenerator::rngSeed.
	std::mt19937 gen(resolve_seed());
	std::uniform_real_distribution<> disX(bounds.xMin + 0.1, bounds.xMax - 0.1);
	std::uniform_real_distribution<> disY(bounds.yMin + 0.1, bounds.yMax - 0.1);
	std::uniform_real_distribution<> disZ(bounds.zMin + 0.1, bounds.zMax - 0.1);

	while (seeds.size() < seedNr) {
		Vec3 pt{ (float)disX(gen), (float)disY(gen), (float)disZ(gen) };
		if (adapter.is_inside(pt)) {
			seeds.push_back(pt);
		}
	}

	update_model();

	version++;
};

void Random::render_gui() {
	
	std::string nr = std::to_string(seeds.size());
	ImGui::Text(("Seeds: " + nr).c_str());

	ImGui::Separator();
	ImGui::SetNextItemWidth(200);
	ImGui::InputInt("Seed Number", &seedNr);
};

void Random::update_model() {

	if (!seeds.empty() && renderMode) {
		std::array<float, 6> tempBounds = {
			bounds.xMin, bounds.xMax, bounds.yMin, bounds.yMax, bounds.zMin, bounds.zMax };
		model = std::make_unique<VisualizeSeeds>(seeds, tempBounds);
		modelSeedSize = model->initialCalculatedSize;
	}
};


// =============================================================================
// Poisson 3D Generator Class Implementation
// =============================================================================
void Poisson3D::run(const IContainer& adapter) {

	if (type == ObjectType::VariedGeneratorType){
		rebuild_config(adapter);
		varied_run(adapter, config);
	}
	else{
		uniform_run(adapter);
	}
};

void Poisson3D::varied_run(const IContainer& adapter, const RunConfig& cfg) {

	containerType = adapter.get_type();

	config = cfg;
	isUniform = false;

	seeds.clear();
	radii.clear();
	grid.clear();

	// define root. Assign the MEMBER 'bounds' (not a local) so update_model()
	// below reads the real container extents; a local here would shadow the
	// member, leaving it uninitialized and giving the seed spheres a garbage
	// scale (huge or invisibly tiny).
	bounds = adapter.compute_bounds();
	root = bounds.center;

	// in the case it is the center of a sphere it will destroy the min distance from sdf
	// we can a add a small noise
	root += Vec3(1e-5f, 1e-5f, 1e-5f);

	// One RNG stream for the whole run, so the root fallback and the candidate
	// sampling below draw from the same sequence and a given rngSeed reproduces
	// the entire seed cloud exactly. rngSeed == 0 -> fresh random_device draw.
	std::mt19937 gen(resolve_seed());

	// if the root is not inside randomly create it
	if (!adapter.is_inside(root)) {
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
	// (reuses the single 'gen' stream seeded at the top of this run)
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
		std::array<int, 3> centralIdx = getGridIndex(root, cellSize);
		grid[centralIdx].seedIdx = N;

		// finally push in the neighbors the index of the seed! 
		// The number of neighbors is determined
		// based on the corresponding radius
		// Stochastic-radius draw: truncated normal N(rMean, radiusStd) on
		// [rMin, rMax]. Rejection (bounded tries, then clamp) keeps the Gaussian
		// shape instead of piling mass at the bounds. Draws from the same seeded
		// 'gen' stream, so the whole cloud (positions AND radii) reproduces for a
		// given rngSeed. rMean <= 0 defaults to the range midpoint.
		const double rMean = (radiusMean > 0.0) ? radiusMean : 0.5 * (rMin + rMax);
		auto draw_radius = [&]() -> double {
			std::normal_distribution<double> nd(rMean, std::max(0.0, radiusStd));
			for (int t = 0; t < 16; ++t) {
				double r = nd(gen);
				if (r >= rMin && r <= rMax) return r;
			}
			return std::clamp(nd(gen), rMin, rMax);
		};

		double rxi;

		// stochastic radius ignores any distance field; else use the graded function
		if (stochasticRadius) {
			rxi = draw_radius();
		}
		else if (cfg.dist) {
			const float d = std::abs(cfg.dist->compute_distance(root));
			rxi = cfg.rad->estimate_radius(d, rMin, rMax);
		}
		else {
			rxi = cfg.rad->estimate_radius(0.0, rMin, rMax);
		}

		int n = ceil(rxi / cellSize) + 1;

		pushIdxs(n, centralIdx, N);

		radii.push_back(rxi);

		while (!active.empty()) {

			// step 1: chose randomly an index from the active list
			// Draw from the seeded 'gen' stream (NOT the unseeded global rand(),
			// which made the seed cloud non-reproducible even at a fixed rngSeed).
			int randomIndex = std::uniform_int_distribution<int>(0, static_cast<int>(active.size()) - 1)(gen);
			int idx = active[randomIndex];

			// step 2: identify the poisson disc parameter for the point,
			Vec3 xi = seeds[idx + currSize];

			// identify the cell of the grid that the point is located
			std::array<int, 3> cellIdx = getGridIndex(xi, cellSize);

			// estimate also the radius based on the distance from the metric,
			// since it will be inside the distance will be negative, thus,
			// we need the absolute value
			double rxi = 0.0f;

			if (stochasticRadius) {
				// reuse this seed's already-assigned radius so its territory (the
				// candidate annulus) is stable across revisits, not re-drawn.
				rxi = radii[idx + currSize];
			}
			else if (cfg.dist) {
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

				float u = disx(gen); // store it to avoid calling the generator twice!
				float r = std::cbrt(u * std::pow(2 * rxi, 3) + (1.0f - u) * std::pow(rxi, 3));

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

				if (stochasticRadius) {
					// each candidate draws its own constraint; if accepted it is
					// stored below and becomes this seed's fixed radius.
					ryi = draw_radius();
				}
				else if (cfg.dist) {
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

	update_model();

	version++;
};

void Poisson3D::uniform_run(const IContainer& adapter) {

	containerType = adapter.get_type();

	seeds.clear();
	radii.clear();
	grid.clear();

	// set max
	rMax = rMin;
	isUniform = true;

	// define root. Assign the MEMBER 'bounds' (not a local) so update_model()
	// below reads the real container extents; a local here would shadow the
	// member, leaving it uninitialized and giving the seed spheres a garbage
	// scale (huge or invisibly tiny).
	bounds = adapter.compute_bounds();

	root = bounds.center;
	// in the case it is the center of a sphere it will destroy the min distance from sdf
	// we can a add a small noise
	//root += Vec3(100.f, 100.f, 100.f);

	// One RNG stream for the whole run, so the root fallback and the candidate
	// sampling below draw from the same sequence and a given rngSeed reproduces
	// the entire seed cloud exactly. rngSeed == 0 -> fresh random_device draw.
	std::mt19937 gen(resolve_seed());

	// if the root is not inside randomly create it
	if (!adapter.is_inside(root)) {
		std::uniform_real_distribution<> disX(bounds.xMin + 0.1, bounds.xMax - 0.1);
		std::uniform_real_distribution<> disY(bounds.yMin + 0.1, bounds.yMax - 0.1);
		std::uniform_real_distribution<> disZ(bounds.zMin + 0.1, bounds.zMax - 0.1);

		while (!adapter.is_inside(root)) {
			root.x = (float)disX(gen);
			root.y = (float)disY(gen);
			root.z = (float)disZ(gen);
		}
	}

	std::cout << "final root: " << root << std::endl;

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
	// (reuses the single 'gen' stream seeded at the top of this run)
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
		std::array<int, 3> centralIdx = getGridIndex(root, cellSize);
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
			// Draw from the seeded 'gen' stream (NOT the unseeded global rand(),
			// which made the seed cloud non-reproducible even at a fixed rngSeed).
			int randomIndex = std::uniform_int_distribution<int>(0, static_cast<int>(active.size()) - 1)(gen);
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

				float u = disx(gen); // store it to avoid calling the generator twice
				float r = std::cbrt(u * std::pow(2 * rxi, 3) + (1.0f - u) * std::pow(rxi, 3));

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

	update_model();

	version++;
};

void Poisson3D::render_gui() {

	if (type == ObjectType::UniformGeneratorType) {

		std::string seedNr = std::to_string(seeds.size());
		ImGui::Text(("Seeds: " + seedNr).c_str());

		ImGui::Separator();
		ImGui::SetNextItemWidth(200);
		ImGui::InputDouble("Radius", &rMin, 0.001, 1.0, "%.3f");

		// ImGui::Checkbox("Fit", &fit);
		// if (fit){
		// 	ImGui::InputDouble("Tb.Sp.", &targetSp);
		// 	ImGui::InputDouble("Tb.Th.", &targetTh);
		// 	ImGui::InputDouble("Offset", &offset);

		// 	// check the container
		// 	if (containerType == ObjectType::BoxContainerType){
		// 		rMin = (targetSp + targetTh + offset) / 1.61;
		// 	}
		// 	else if (containerType == ObjectType::CylinderContainerType){
		// 		rMin = (targetSp + targetTh + offset) / 1.084;
		// 	}
		// }
	}

	else if (type == ObjectType::VariedGeneratorType) {

		std::string seedNr = std::to_string(seeds.size());
		ImGui::Text(("Seeds: " + seedNr).c_str());

		ImGui::Separator();

		ImGui::SeparatorText("Select Distance Function");
		ImGui::RadioButton("Distance From Plane", &distIdx, 0);
		if (distIdx== 0) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat3("Normal", planeNormal);
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat3("Center", planeCenter);
		};
		ImGui::RadioButton("Distance From Point", &distIdx, 1);
		if (distIdx == 1) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat3("Point", point);
		}
		ImGui::RadioButton("Distance From Container", &distIdx, 2);

		ImGui::SeparatorText("Select Radius Function");
		ImGui::RadioButton("Linear", &radiusIdx, 0);
		ImGui::RadioButton("Quadratic", &radiusIdx, 1);
		ImGui::RadioButton("Sigmoid", &radiusIdx, 2);
		ImGui::RadioButton("Random (Normal)", &radiusIdx, 3);
		if (radiusIdx == 3) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputDouble("Mean Radius (0 = midpoint)", &radiusMean);
			ImGui::SetNextItemWidth(200);
			ImGui::InputDouble("Radius Std", &radiusStd);
			ImGui::TextDisabled("Each seed draws r ~ N(mean, std) clamped to [rMin, rMax].");
		}

		ImGui::SeparatorText("Parameters");

		ImGui::SetNextItemWidth(200);
		ImGui::InputDouble("Start Radius", &rMin);
		ImGui::SetNextItemWidth(200);
		ImGui::InputDouble("End Radius", &rMax);
		ImGui::SetNextItemWidth(200);
		ImGui::InputDouble("Transition Distance", &transitionDist);
	}
};

void Poisson3D::set_min_radius(const double newRadius) { rMin = newRadius; };

void Poisson3D::set_max_radius(const double newRadius) { rMax = newRadius; };

void Poisson3D::set_center(const Vec3& newCenter) { root = newCenter; };

void Poisson3D::update_model() {
	if (!seeds.empty() && renderMode) {
		std::array<float, 6> tempBounds = {
			bounds.xMin, bounds.xMax, bounds.yMin, bounds.yMax, bounds.zMin, bounds.zMax };
		model = std::make_unique<VisualizeSeeds>(seeds, tempBounds);
		modelSeedSize = model->initialCalculatedSize;
	}
};

void Poisson3D::rebuild_config(const IContainer& con){
	switch (radiusIdx) {
		// linear radius function
		case 0: {
			config.rad = std::make_shared<LinearFunction>();
			break;
		}
		case 1: {
			config.rad = std::make_shared<QuadraticFunction>(transitionDist);
			break;
		}
		case 2: {
			config.rad = std::make_shared<SmoothStep>(transitionDist);
			break;
		}
		case 3: {
			// stochastic: no radius function, handled inside run()
			break;
		}
	}

	stochasticRadius = (radiusIdx == 3);

	if (!stochasticRadius)
	switch (distIdx) {
		// distance from plane
		case 0: {
			config.dist = std::make_shared<PlaneSDF>(planeCenter, planeNormal);
			break;
		}
		// distance from point
		case 1: {
			config.dist = std::make_shared<PointSDF>(point);
			break;
		}
		// distance from container surface
		case 2: {
			config.dist = con.get_distance_estimator();
			break;
		}
	}

	//ensure distance func is nullptr for random
	if (radiusIdx == 3) {
		config.dist = nullptr;
	}	
};

// =============================================================================
// Seed Generator Factory
// =============================================================================

void SeedGeneratorFactory::launch(){

	selectedCon = nullptr;
	name = "";
	warningFlashTimer = 0.0f;
};

IContainer* SeedGeneratorFactory::gui_draw(
	Logger* logger,
    const char* popupName, bool& showPopup,
	SelectedObject* selectedPanelObj,
	std::vector<std::shared_ptr<InterfaceSeedGenerator>>& generators,
	const std::vector<std::shared_ptr<IContainer>>& containers
){

	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}
	else return nullptr;

	// showPopup is only ever set true alongside a set_type() call that builds
	// pendingGenerator; guard anyway so a stray flag can never dereference null.
	if (!pendingGenerator) return nullptr;

	// always centered
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
	ImGui::SetNextWindowSize(ImVec2(600, 500), ImGuiCond_Appearing);

	if (ImGui::BeginPopupModal(popupName, NULL))
	{

		ImGui::InputText("Name", &name);

		ImGui::SeparatorText("Select Container");

		if (warningFlashTimer > 0.0f) {
			warningFlashTimer -= ImGui::GetIO().DeltaTime;
		}

		bool isFlashing = (warningFlashTimer > 0.0f);
		if (isFlashing) {
			float pulseAlpha = (float)(std::sin(ImGui::GetTime() * 15.0f) * 0.5f + 0.5f);
			ImVec4 flashColor = ImVec4(1.0f, 0.0f, 0.0f, pulseAlpha);
			ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 5.0f);
			ImGui::PushStyleColor(ImGuiCol_Border, flashColor);
		}

		ImGui::BeginChild("Containers", ImVec2(ImGui::GetContentRegionAvail().x, 80), ImGuiChildFlags_Borders);

		for (const auto& md : containers) {

			IContainer* cptr = dynamic_cast<IContainer*>(md.get());

			if (cptr) {
				bool isSelected = (selectedCon && selectedCon == cptr);

				if (ImGui::Selectable(md->name.c_str(), isSelected)) {
					selectedCon = md.get();
				};
			}
		}
		if (isFlashing) {
			ImGui::PopStyleColor();
			ImGui::PopStyleVar();
		}
		
		ImGui::EndChild();

		// The pending generator renders its own type-specific controls, editing
		// its members in place (same render_gui() the properties panel uses).
		pendingGenerator->render_gui();

		ImGui::Separator();
		ImGui::NewLine();
		if (ImGui::Button("Create")){

			if(selectedCon){

				std::string msg =
					(genType == ObjectType::RandomGeneratorType)  ? " random "  :
					(genType == ObjectType::UniformGeneratorType) ? " uniform " :
					(genType == ObjectType::VariedGeneratorType)  ? " varied "  : " ";

				pendingGenerator->container = selectedCon;
				pendingGenerator->run(*selectedCon);

				if (!name.empty()) {
					pendingGenerator->name = name;
				}
				else {
					pendingGenerator->name = "Generator" + std::to_string(generators.size() + 1);
				}

				size_t nr = pendingGenerator->get_seeds().size();
				generators.push_back(std::move(pendingGenerator));
				selectedPanelObj->ptr = generators.back().get();
				selectedPanelObj->type = genType;
				
				logger->log(
					LogPriority::SUCCESS,
					std::to_string(nr) + msg + " seeds created inside " + selectedCon->name + "!" 
				);

				IContainer* con = selectedCon;

				ImGui::CloseCurrentPopup();
				showPopup = false;
				ImGui::EndPopup();
				return con;
			}
			else{
				warningFlashTimer = 1.5f;
			}
		}

		ImGui::SameLine();

		if (ImGui::Button("Cancel")){
			showPopup = false;
			ImGui::CloseCurrentPopup();
			ImGui::EndPopup();
			return nullptr;
		}

		ImGui::EndPopup();
	}
	return nullptr;
};

void SeedGeneratorFactory::set_type(const ObjectType& newType){
	genType = newType;

	// Build the pending generator up front, with its type set, so gui_draw can
	// render its controls (render_gui) and run it (run) polymorphically. The
	// generator's own members hold the parameters - the factory keeps no copy.
	switch (genType)
	{
	case ObjectType::RandomGeneratorType:
		pendingGenerator = std::make_shared<Random>(100);
		break;
	case ObjectType::UniformGeneratorType:
		pendingGenerator = std::make_shared<Poisson3D>(1.0, 1.0, 30, true);
		break;
	case ObjectType::VariedGeneratorType:
		pendingGenerator = std::make_shared<Poisson3D>(1.0, 1.2, 30, true);
		break;
	default:
		pendingGenerator = nullptr;
		return;
	}

	pendingGenerator->type = genType;
};
