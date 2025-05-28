#ifndef POISSON3D_H // include guard
#define POISSON3D_H

#include <vector>
#include <array>
#include <random>
#include <unordered_set>
#include <algorithm>
#include <functional>
#include <vtkPolyData.h>
#include <vtkAppendPolyData.h>
#include <vtkImplicitModeller.h>
#include <vtkNamedColors.h>
#include <vtkContourFilter.h>
#include <vtkSTLWriter.h>
#include "DistanceCalculator.h"
#include "Visualize/visualize.h"
#include "voro++.hh"

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

//struct GridHashFunc {
//	std::size_t operator()(const std::array<int, 3>& arr) const {
//		std::size_t seed = 0;
//		for (int i = 0; i < 3; i++) {
//			seed ^= std::hash<int>{}(arr[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//		}
//		return seed;
//	}
//};

class RadiusFunction {
public:
	virtual double estimate_radius(double distance, double rMin, double rMax) const = 0;
	virtual ~RadiusFunction() = default;
};

class LinearFunction : public RadiusFunction {

public:
	LinearFunction() {};
	~LinearFunction() {};
	LinearFunction(double distMax) : distMax(distMax) {};

	double estimate_radius(double distance, double rMin, double rMax) const override {
		if (distance < 0) {
			return rMin;
		}
		else if (distance > distMax) {
			return rMax;
		}
		else {
			return ((rMax - rMin) / distMax) * distance + rMin;
		}
	}
	//double estimate_radius(double distance, double rMin, double rMax) const override {
	//	return	rMax - (rMax - rMin) * std::min(distance / 5.0, 1.0);
	//};

private:
	double distMax{ 10.0 };
};

// -----------------------------------------------------------------------------------------------
// base class definition

class Poisson3DSampler {

public:
	Poisson3DSampler() {};
	~Poisson3DSampler() {};
	Poisson3DSampler(
		const float rMin,
		const float rMax,
		const std::array<double, 3>& root) : rMin(rMin), rMax(rMax), root(root) {};

	virtual void generate_seeds(DistanceEstimator& distEstimator, RadiusFunction& radFunc) = 0;
	void generate_voro();
	void generate_mesh(const double& thickness, const std::string& fileName);
	void get_seeds(std::vector<std::array<double, 3>>& outSeeds) {
		outSeeds = seeds;
	};
	void get_radii(std::vector<float>& outRadii) {
		outRadii = radii;
	};

protected:
	std::array<double, 3> root;
	int neighNr{ 30 };
	double PI = 3.14159265358979323846;
	int seedNr{ 0 };
	std::vector<float> radii;
	std::unordered_map < std::array<int, 3>, Cell, GridHashFunc> grid;
	float cellSize{ 0.0f };
	float rMin{ 0.0f }, rMax{ 0.0f };
	double xMin, xMax, yMin, yMax, zMin, zMax;
	double scale{ 1.0 };
	std::vector<std::array<double, 3>> seeds;
	int nX{ 1 }, nY{ 1 }, nZ{ 1 };
	voro::container* con;
	vtkSmartPointer<vtkPolyData> scaffoldMesh;

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


class Poisson3D : public Poisson3DSampler {

public:
	Poisson3D() {};
	~Poisson3D() {};
	Poisson3D(
		const float& rMinVal,
		const float& rMaxVal,
		const std::array<double, 3>& rootVal,
		const std::array<float, 6>& bounds,
		std::function<bool(const std::array<double, 3>&)> is_inside
	);

	void generate_seeds(DistanceEstimator& distEstimator, RadiusFunction& radFunc) override;
	void generate_seeds();

	// create a function pointer
	std::function<bool(const std::array<double, 3>&)> is_inside;
};


class Poisson3DWall : public Poisson3DSampler {

public:
	Poisson3DWall() {};
	~Poisson3DWall() {};
	Poisson3DWall(
		vtkSmartPointer<vtkPolyData>& containerMesh,
		const float& rMinVal,
		const float& rMaxval,
		const int& neighbors,
		const std::array<int, 3>& blockDim, 
		const double wallResolution);
	void generate_seeds(DistanceEstimator& distEstimator, RadiusFunction& radFunc) override;
	void generate_seeds();
	void generate_voro();

private:
	vtkSmartPointer<vtkPolyData> container;
	vtkSmartPointer<vtkImplicitPolyDataDistance> distanceCalculator;
	int neighbors{ 1 };
	double wallResolution{ 1.0 };
	std::vector<std::array<double, 3>> bCenters;
	std::vector<std::array<double, 3>> normals;
	void _populate_with_barycenters();

};


#endif
