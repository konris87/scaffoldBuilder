#ifndef SCAFFOLDGENERATOR_H // include guard
#define SCAFFOLDGENERATOR_H

#include <vector>
#include <array>
#include <string>
#include <voro++.hh>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <Eigen/Dense>
#include <unordered_map>
#include <unordered_set>

// hash function to compute hash from key
struct GlobalVertexKey {
	
	std::size_t operator()(const std::array<double, 3>& coords) const {
		const double scale = 1e6;

		// scale the coordinates to convert to integers
		// for example 0.3456778 -> 345677
		std::size_t h1 = std::hash<int>{}(static_cast<int>(coords[0] * scale));
		std::size_t h2 = std::hash<int>{}(static_cast<int>(coords[1] * scale));
		std::size_t h3 = std::hash<int>{}(static_cast<int>(coords[2] * scale));

		// create the key by combining hashes
		std::size_t seed = h1;
		
		seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);

		return seed;
	};
};

// function to compare keys for the vertex map
struct GlobalVertexInside {

	const double precision = 1e-6;

	bool operator()(const std::array<double, 3>& a1, const std::array<double, 3>& a2) const {

		return
			std::fabs(a1[0] - a2[0]) < precision &&
			std::fabs(a1[1] - a2[1]) < precision &&
			std::fabs(a1[2] - a2[2]) < precision;
	};
};

// Map for global faces

// function to create the key
struct GlobalFaceKey {
	std::size_t operator()(const std::vector<int>& face) const {
		std::size_t seed = face.size();
		for (int i : face) {
			seed ^= std::hash<int>{}(i)+0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

struct GlobalFaceInside {

	bool operator()(const std::vector<int>& face1, const std::vector<int>& face2) const {
		return face1 == face2;
	};

};


class ScaffoldGenerator{
public:
	
	virtual ~ScaffoldGenerator() {};
	ScaffoldGenerator(
		std::vector<std::array<double, 3>>& seeds,
		const std::array<float, 6>& bounds,
		const std::array<int, 3>& blockDim
	) : seeds(seeds), bounds(bounds), blockDim(blockDim) {};
	
	void generate_mesh(const double& thickness, vtkSmartPointer<vtkPolyData>& finalPolyData, const std::vector<int>& res = {});
	virtual void generate_voro(const int regSteps) = 0;
	void get_seeds(std::vector<std::array<double, 3>>& outseeds) { seeds = outseeds; };

protected:
	std::vector<std::array<double, 3>>& seeds;
	std::array<float, 6> bounds;
	std::array<int, 3> blockDim;
	
	// voro container
	voro::container* con;

	// vtk scaffold mesh
	vtkSmartPointer<vtkPolyData> scaffoldMesh;
};

// -----------------------------------------------------
// scaffold generator inside box
class ScaffoldGeneratorBox : public ScaffoldGenerator {

public:
	ScaffoldGeneratorBox(
		std::vector<std::array<double, 3>>& seeds,
		const std::array<float, 6>& bounds,
		const std::array<int, 3>& blockDim
	);
	~ScaffoldGeneratorBox() override {}

	void generate_voro(const int regSteps) override;
};

// -----------------------------------------------------
// scaffold generator with custom Wall
class ScaffoldGeneratorWall : public ScaffoldGenerator {
public:
	~ScaffoldGeneratorWall() override {}
	ScaffoldGeneratorWall(
		std::vector<std::array<double, 3>>& seeds,
		vtkSmartPointer<vtkPolyData>& containerPoly,
		const std::array<int, 3>& blockDim,
		const int neighbors,
		const float minDist);

	void generate_voro(const int regSteps) override;

private:

	std::vector<std::array<double, 3>> bCenters;
	std::vector<std::array<double, 3>> normals;
	vtkSmartPointer<vtkPolyData> containerMesh;
	float minDist{ 0.0 };
	int neighbors{ 1 };
	void _process_triangles();
};


// -------------------------------------------------------------
// Volume Optimization

class VolOpt {

public:
	VolOpt() {};
	~VolOpt() {};

	// constructor overloading on for the generate points option one for the load points from file
	VolOpt(
		std::vector<std::array<double, 3>>& seeds,
		const Eigen::VectorXd targets,
		const Eigen::VectorXd initV,
		const std::array<double, 6>& bounds,
		std::function<void(const std::string&)> logCallback = nullptr
	);

	void loop(const int regSteps);
	void generate_mesh(
		const double& thickness,
		vtkSmartPointer<vtkPolyData>& finalPolyData,
		const std::vector<int>& res);
		//const std::string& fileName);
	void get_seeds(std::vector<std::array<double, 3>>& outseeds);

private:
	std::vector<std::array<double, 3>> currSeeds;
	std::array<double, 6> bounds{ 0.0, 10.0, 0.0, 10.0, 0.0, 10.0 };
	std::function<void(const std::string&)> log_callback;
	Eigen::VectorXd targetVols;
	Eigen::VectorXd wInit;
	std::vector<vtkSmartPointer<vtkPolyData>> polys;
	vtkSmartPointer<vtkPolyData> scaffoldMesh;
};

// ----------------------------------------------------------------------
// Volume Optimization inside a Mesh
class VolOptWall {

public:
	VolOptWall() {};
	~VolOptWall() {};

	// constructor overloading on for the generate points option one for the load points from file
	VolOptWall(
		std::vector<std::array<double, 3>>& seeds,
		const Eigen::VectorXd targets,
		const Eigen::VectorXd initV,
		vtkSmartPointer<vtkPolyData>& containerPoly,
		int neighbors = 1,
		std::function<void(const std::string&)> logCallback = nullptr
	);

	void loop(const int regSteps);
	void generate_mesh(
		const double& thickness,
		const std::string& fileName,
		const std::vector<int>& res);
	void get_seeds(std::vector<std::array<double, 3>>& outseeds) {
		outseeds = currSeeds;
	};

private:
	std::vector<std::array<double, 3>> currSeeds;
	std::array<double, 6> bounds{ 0.0, 10.0, 0.0, 10.0, 0.0, 10.0 };
	std::function<void(const std::string&)> log_callback;
	Eigen::VectorXd targetVols;
	Eigen::VectorXd wInit;
	std::vector<vtkSmartPointer<vtkPolyData>> polys;
	vtkSmartPointer<vtkPolyData> scaffoldMesh;
	vtkSmartPointer<vtkPolyData> containerMesh;
	std::vector<std::array<double, 3>> bCenters;
	std::vector<std::array<double, 3>> normals;
	int neighbors{ 1 };
	void _process_triangles();
};

//---------------------------------------------------------------------
// Scaffold Generator using face editing

class ScaffoldGeneratorFaceBox : public ScaffoldGenerator {

public:

	// constructor using range of hole radius and interpolation for edge size
	ScaffoldGeneratorFaceBox(
		std::vector<std::array<double, 3>>& seeds,
		const std::array<float, 6>& bounds,
		const std::array<int, 3>& blockDim,
		double minRad, double maxRad, double edgeSize);

	~ScaffoldGeneratorFaceBox() override {};
	void generate_voro(const int regSteps) override {};
	void generate_voro();

protected:

	double minHoleRadius{ 0.0 };
	double maxHoleRadius{ 0.0 };
	double edgeSize{ 0.0 };

private: 
	std::vector<int> globalIndices;
	std::vector<std::array<double, 3>> globalVertices;
	std::vector<std::vector<int>> globalFaces;
	int globalIndex{ 0 };
	std::unordered_map<std::array<double, 3>, int, GlobalVertexKey, GlobalVertexInside> globalVertexMap;
	std::unordered_set<std::vector<int>, GlobalFaceKey, GlobalFaceInside> globalFaceMap;
};

#endif