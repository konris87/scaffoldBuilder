#ifndef SCAFFOLDGENERATOR_H // include guard
#define SCAFFOLDGENERATOR_H

#include <vector>
#include <array>
#include <string>
#include <voro++.hh>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <Eigen/Dense>

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



#endif