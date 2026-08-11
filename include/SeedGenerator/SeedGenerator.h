#ifndef SEED_GENERATOR_H
#define SEED_GENERATOR_H

#include <array>
#include <vector>
#include <functional>
#include <optional>
#include <random>
#include "Utils/Utils.h"
#include <Visualize/VisualizeSeeds.h>
#include "Container.h"
#include "Math/Vec.h"
#include "Logger/Logger.h"

using DistFunc = std::function<double(const Pt&)>;
using RadFunc = std::function<double(double, double, double)>;

struct RunConfig {
	std::shared_ptr<const SDF> dist; 
	std::shared_ptr<const RadiusFunction> rad;
};

class InterfaceSeedGenerator {
public:
	virtual ~InterfaceSeedGenerator() = default;

	virtual void run(const IContainer& adapter) = 0;

	virtual void render_gui() = 0;

	virtual ObjectType get_type() = 0;

	virtual void update_model() = 0;

	void draw() { 
		if (model) model->draw();
		else { std::cout << "no model found" << std::endl; }
	};

	void set_renderMode(bool render) { renderMode = render; };

	// RNG seed for the stochastic seed placement.
	//   non-zero (DEFAULT 1) -> deterministic: the same value reproduces the exact
	//                  same seed cloud, so regenerating at identical parameters
	//                  gives an identical scaffold (and identical metrics, incl.
	//                  DA). Distinct values give independent realizations
	//                  (use 1..N for N replicates).
	//   0 -> opt-in non-deterministic: a fresh std::random_device draw per run, so
	//                  every run is a different realization.
	// The Poisson radius fixes the MINIMUM SPACING; this fixes WHICH of the many
	// valid arrangements at that spacing you get.
	//
	// The default is a FIXED seed (not 0) so the GUI - which does not expose this -
	// is reproducible by default: a deterministic modelling tool should map
	// identical inputs to identical output. For a distribution of realizations use the profiler (--replicates / --seed), which sets this explicitly.
	uint32_t rngSeed = 1;

	void set_rng_seed(uint32_t s) { rngSeed = s; };

	// Resolves the seed actually handed to the generator for one run.
	uint32_t resolve_seed() const {
		if (rngSeed != 0) return rngSeed;
		std::random_device rd;
		return static_cast<uint32_t>(rd());
	};

	std::vector<Vec3> get_seeds() { return seeds; };

	virtual void set_seeds(const std::vector<Vec3>& newSeeds) { 
		seeds = newSeeds;
	};

	std::string name = "";
	bool hidden = false;
	ObjectType type = ObjectType::NoneType;
	ObjectType containerType = ObjectType::NoneType;

	IContainer* container = nullptr;
	uint32_t version = 1;
	float modelSeedSize = 1.0f;

	bool renderMode = true;

protected:
	Bounds bounds;
	std::unique_ptr<VisualizeSeeds> model;
	std::vector<Vec3> seeds;
};

class Random final : public InterfaceSeedGenerator {
public:

	Random() {};

	explicit Random(int n, const bool render = true) : seedNr(n) { renderMode = render; };
	
	void run(const IContainer& adapter) override;

	void render_gui() override;

	void update_model() override;

	void set_seeds(const std::vector<Vec3>& newSeeds) override {
		seeds = newSeeds;
		seedNr = (int)seeds.size();
	};

	ObjectType get_type() override { return ObjectType::RandomGeneratorType; };
	
	int seedNr{ 0 };

};

// -----------------------------------------------------------------------------------------
// cell structure to hold information of background cell
struct Cell {
	std::vector<int> cellIdxs;
	float radius;
	int seedIdx = -1;
};

//custom hash to work with keys of type std::array<int, 3> these are the
//indices of our cell based on the position of a point
struct GridHashFunc {
	std::size_t operator()(const std::array<int, 3>& arr) const {
		std::size_t h1 = std::hash<int>{}(arr[0]);
		std::size_t h2 = std::hash<int>{}(arr[1]);
		std::size_t h3 = std::hash<int>{}(arr[2]);
		return h1 ^ (h2 << 1) ^ (h3 << 2);
	}
};

class Poisson3D final : public InterfaceSeedGenerator {
public:

	Poisson3D() {};

	Poisson3D(const double& rMinVal, const double& rMaxVal, int neighs, const bool render=true)
		: rMin(rMinVal), rMax(rMaxVal), neighNr(neighs) {
		
		renderMode = render;
	};

	void varied_run(const IContainer& adapter, const RunConfig& cfg);

	void uniform_run(const IContainer& adapter);

	void run(const IContainer& adapter) override;

	void render_gui() override;

	void update_model() override;

	void set_min_radius(const double newRadius);

	void set_max_radius(const double newRadius);

	void set_center(const Vec3& newCenter);

	float get_min_radius() const { return (float)rMin; };

	float get_max_radius() const { return  (float)rMax; };

	RunConfig get_config() const { return config; };

	ObjectType get_type() override { return ObjectType::PoissonGeneratorType; };

	bool is_uniform() const { return isUniform; };

	void rebuild_config(const IContainer& con);

	// keep here the indices for distance and radius function
	int radiusIdx = 0;
	int distIdx = 0;
	double rMin{ 0.0 }, rMax{ 0.0 };

	// distance references
	double transitionDist = 1.0;
	Vec3 planeNormal{ 0.0f, 0.0f, 1.0f };
	Vec3 planeCenter{ 0.0f, 0.0f, 0.0f };
	Vec3 point{ 0.0f, 0.0f, 0.0f };

	// Stochastic-radius mode (a third varied-Poisson option, no distance field):
	// each seed draws its own minimum-distance constraint from a truncated normal
	// N(radiusMean, radiusStd) clamped to [rMin, rMax]. This injects INTRA-sample
	// spacing heterogeneity (a spread of cell sizes) for structural realism, e.g.
	// to approach a literature Tb.Sp distribution. The draws come from the run's
	// single seeded RNG stream, so a given rngSeed still reproduces the whole cloud.
	// NOTE: input radiusStd is NOT the output Tb.Sp std (the packing regularises the
	// variance and radius->Tb.Sp is a monotone transfer, not identity) - use the
	// profiler radius-std sweep to map it. radiusMean <= 0 defaults to (rMin+rMax)/2.
	bool stochasticRadius = false;
	double radiusMean = 0.0;
	double radiusStd = 0.0;

private:
	Vec3 root{ 0.0, 0.0, 0.0 };
	std::vector<float> radii;
	std::unordered_map < std::array<int, 3>, Cell, GridHashFunc> grid;
	int neighNr{ 30 };
	float PI = 3.14159265358979323846f;
	RunConfig config;
	bool isUniform = true;
	bool fit = false;
	double targetSp = 1.0;
	double targetTh = 0.4;
	double offset = 0.35;

	// protected functions
	void pushIdxs(const int& n, const std::array<int, 3> cellIdx, const int& N) {
		for (int i = cellIdx[0] - n; i < cellIdx[0] + n + 1; i++) {
			for (int j = cellIdx[1] - n; j < cellIdx[1] + n + 1; j++) {
				for (int k = cellIdx[2] - n; k < cellIdx[2] + n + 1; k++) {
					std::array<int, 3> idx = { i, j, k };
					grid[idx].cellIdxs.push_back(N);
				}
			}
		}
	};

	std::array<int, 3> getGridIndex(const Vec3& p, float cellSize) {
		return {
			static_cast<int>(std::floor(p.x / cellSize)),
			static_cast<int>(std::floor(p.y / cellSize)),
			static_cast<int>(std::floor(p.z / cellSize)) };
	};
};

// =============================================================================
// Seed Generator Factory
// =============================================================================

class SeedGeneratorFactory{
public:
	SeedGeneratorFactory(){};
	~SeedGeneratorFactory(){};

	void launch();

	IContainer* gui_draw(
		Logger* logger,
        const char* popupName, bool& showPopup,
        SelectedObject* selectedPanelObj,
		std::vector<std::shared_ptr<InterfaceSeedGenerator>>& generators,
		const std::vector<std::shared_ptr<IContainer>>& containers
	);

	void set_type(const ObjectType& newType);
private:

	// The generator being built. Created (with its type) in set_type(), edited
	// in-place through pendingGenerator->render_gui(), and moved into the scene
	// list on "Create". This is the single source of truth for the new
	// generator's parameters - there is no parallel copy of them on the factory.
	std::shared_ptr<InterfaceSeedGenerator> pendingGenerator = nullptr;

	ObjectType genType = ObjectType::NoneType;
	IContainer* selectedCon = nullptr;

	// rendering gui
	std::string name = "";
	float warningFlashTimer = 0.0f;
};

#endif