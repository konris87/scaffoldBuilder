#ifndef SEED_GENERATOR_H
#define SEED_GENERATOR_H

#include <array>
#include <functional>
#include <optional>
#include "Utils/Utils.h"
#include "Container.h"

using DistFunc = std::function<double(const Pt&)>;
using RadFunc = std::function<double(double, double, double)>;

struct RunConfig {
	std::shared_ptr<const DistanceEstimator> dist; 
	std::shared_ptr<const RadiusFunction> rad;
};

class InterfaceSeedGenerator {
public:
	virtual ~InterfaceSeedGenerator() = default;
	virtual void run(
		const ContainerAdapter& adapter,
		std::vector<std::array<double, 3>>& out,
		const RunConfig& cfg = {}) = 0;
};

class Random final : InterfaceSeedGenerator {
public:
	explicit Random(int n) : seedNr(n) {};
	
	void run(
		const ContainerAdapter& adapter,
		std::vector<std::array<double, 3>>& out,
		const RunConfig& cfg = {}) override;

private:
	int seedNr{ 0 };
};

// -----------------------------------------------------------------------------------------
// cell structure to hold information of background cell
struct Cell {
	std::vector<int> cellIdxs;
	double radius;
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
	Poisson3D(const float& rMinVal, const float& rMaxVal, const std::array<double, 3>& rootVal, int neighs)
		: rMin(rMinVal), rMax(rMaxVal), root(rootVal), neighNr(neighs) {};

	void run(
		const ContainerAdapter& adapter,
		std::vector<std::array<double, 3>>& out,
		const RunConfig& cfg = {}) override;

private:
	std::array<double, 3> root{ 0.0, 0.0, 0.0 };
	float rMin{ 0.0 }, rMax{ 0.0 };
	std::vector<float> radii;
	std::unordered_map < std::array<int, 3>, Cell, GridHashFunc> grid;
	int neighNr{ 30 };
	double PI = 3.14159265358979323846;

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

	std::array<int, 3> getGridIndex(const std::array<double, 3>& p, double cellSize) {
		return {
			static_cast<int>(std::floor(p[0] / cellSize)),
			static_cast<int>(std::floor(p[1] / cellSize)),
			static_cast<int>(std::floor(p[2] / cellSize)) };
	};
};


#endif