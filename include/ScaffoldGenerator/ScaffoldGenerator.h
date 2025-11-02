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
#include <functional>
#include <Utils/Utils.h>


//Implementing the strategy pattern
//@brief class to act as the interface
class GeneratorInterface {
public:
	virtual ~GeneratorInterface() = default;
	virtual void run_generate_voro(int regSteps) = 0;
	virtual float run_estimate_connectivity() = 0;
	virtual void run_generate_mesh(
		const double& thickness, 
		vtkSmartPointer<vtkPolyData>& finalPolyData,
		const std::array<int, 3>& res = {}) = 0;
	virtual void run_get_seeds(std::vector<std::array<double, 3>>& outseeds) = 0;
	virtual void run_process_faces() = 0;
	virtual std::vector<float> run_get_connected_vertices() = 0;
	virtual std::vector<unsigned int> run_get_connected_edges() = 0;
};

// hash function to compute hash from key
struct GlobalVertex {

	std::array<double, 3> coords;
	
	GlobalVertex(std::array<double, 3>& key) : coords(key) {};

	bool operator==(const GlobalVertex& other) const {
		const double precision = 1e-6;

		return std::fabs(coords[0] - other.coords[0]) < precision &&
			std::fabs(coords[1] - other.coords[1]) < precision &&
			std::fabs(coords[2] - other.coords[2]) < precision;
	};
};

// function to compare keys for the vertex map
struct GlobalVertexHash {

	std::size_t operator()(const GlobalVertex& v) const {
		const double inv_tol = 1e6;
		int x = static_cast<int>(std::llround(v.coords[0] * inv_tol));
		int y = static_cast<int>(std::llround(v.coords[1] * inv_tol));
		int z = static_cast<int>(std::llround(v.coords[2] * inv_tol));

		return static_cast<std::size_t>(
			(x * 73856093) ^ (y * 19349663) ^ (z * 83492791)
			);
	}
};

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

struct CellFace {
	std::vector<double> vertices;
	std::vector<int> localFace;
	bool createHole = true;
	float delta{ 0.0f };
};

struct NewFace {
	Eigen::MatrixXd vertices;
	std::vector<std::vector<int>> faces;
};

class ScaffoldGenerator : public GeneratorInterface{
public:
	virtual ~ScaffoldGenerator() {};
	ScaffoldGenerator(
		std::vector<std::array<double, 3>>& seeds,
		const std::array<float, 6>& bounds,
		const std::array<int, 3>& blockDim,
		const double thickness,
		const double pullbackRatio,
		std::function<bool(const std::array<double, 3>&)> is_inside
	) : seeds(seeds), bounds(bounds), blockDim(blockDim), thickness(thickness), pullbackRatio(pullbackRatio), is_inside(is_inside) {};
	
	void generate_mesh(
		const double& thickness,
		vtkSmartPointer<vtkPolyData>& finalPolyData,
		const std::array<int, 3>& res = {});

	void run_generate_mesh(
		const double& thickness,
		vtkSmartPointer<vtkPolyData>& finalPolyData,
		const std::array<int, 3>& res = {}) override {
		this->generate_mesh(thickness, finalPolyData, res);
	}

	float estimate_connectivity();

	float run_estimate_connectivity() override {
		return this->estimate_connectivity();
	}

	void get_seeds(std::vector<std::array<double, 3>>& outseeds) { seeds = outseeds; };
	void run_get_seeds(std::vector<std::array<double, 3>>& outSeeds) override {
		this->get_seeds(outSeeds);
	};

	void process_faces();

	void run_process_faces() override {
		this->process_faces();
	};

	std::vector<float> get_connected_vertices();
	std::vector<float> run_get_connected_vertices() {
		return this->get_connected_vertices();
	}

	std::vector<unsigned int> get_connected_edges();
	std::vector<unsigned int> run_get_connected_edges() {
		return this->get_connected_edges();
	};


protected:
	std::vector<std::array<double, 3>>& seeds;
	std::array<float, 6> bounds;
	std::array<int, 3> blockDim;

	// connectivity network
	std::vector<float> connectedVertices;
	std::vector<unsigned int> connectedIndices;
	
	// voro container
	voro::container* con = nullptr;

	// vtk scaffold mesh
	vtkSmartPointer<vtkPolyData> scaffoldMesh;

	double thickness{ 0.0 };
	double pullbackRatio{ 0.0 };

	std::vector<int> globalIndices;
	std::vector<std::array<double, 3>> globalVertices;
	std::vector<std::vector<int>> globalFaces;
	std::vector<CellFace> allCellFaces;

	std::unordered_map<GlobalVertex, int, GlobalVertexHash> globalVertexMap;

	// this map is only used for the faces that consist of the vertices of voro cells
	std::unordered_set<std::vector<int>, GlobalFaceKey, GlobalFaceInside> globalFaceMap;

	std::unordered_map <int, int> localToGlobal;

	int globalPolyIndex{ 0 };
	int globalIndex{ 0 };

	std::function<bool(const std::array<double, 3>&)> is_inside;

	std::unique_ptr<Graph> graph;

	std::unordered_map<int, std::unordered_map<int, Eigen::Vector3d>> centroids;

	//std::unordered_map<int, std::vector<Eigen::Vector3d

	void regularize_voro(int regSteps);
};

// -----------------------------------------------------
// scaffold generator inside box
class ScaffoldGeneratorBox : public ScaffoldGenerator {

public:

	// constructor using range of hole radius and interpolation for edge size
	ScaffoldGeneratorBox(
		std::vector<std::array<double, 3>>& seeds,
		const std::array<float, 6>& bounds,
		const std::array<int, 3>& blockDim,
		const double thickness,
		const double pullbackRatio,
		std::function<bool(const std::array<double, 3>&)> is_inside
	);

	~ScaffoldGeneratorBox() override {};

	void populate_voro(const int regSteps);

	// update the interface 
	void run_generate_voro(const int regSteps) override {
		this->populate_voro(regSteps);
	}

	void add_cylindrical_wall(
		const double pt0, const double pt1, const double pt2,
		const double axis0,
		const double axis1,
		const double axis2,
		const double radius);

protected:

	double minHoleRadius{ 0.0 };
	double maxHoleRadius{ 0.0 };

private:

	std::unique_ptr<voro::wall> wall;

	std::unordered_map<GlobalVertex, int, GlobalVertexHash> globalVertexMap;

	// this map is only used for the faces that consist of the vertices of voro cells
	std::unordered_set<std::vector<int>, GlobalFaceKey, GlobalFaceInside> globalFaceMap;
	
	std::unordered_map <int, int> localToGlobal;

	void check_global_vertices(const Eigen::MatrixXd& vertices, std::vector<int>* face = nullptr);
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
		const float minDist,
		const double thickness,
		const double pullbackRatio,
		std::function<bool(const std::array<double, 3>&)> is_inside
	);

	void generate_voro(const int regSteps);

	void run_generate_voro(const int regSteps) override {
		this->generate_voro(regSteps);
	}

private:

	std::vector<std::array<double, 3>> bCenters;
	std::vector<std::array<double, 3>> normals;
	vtkSmartPointer<vtkPolyData> containerMesh;
	float minDist{ 0.0 };
	int neighbors{ 1 };

	void _process_triangles();
	void check_global_vertices(const Eigen::MatrixXd& vertices, std::vector<int>* face = nullptr);
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
		const std::array<int, 3>& res);

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
		const std::array<int, 3>& res);
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


#endif