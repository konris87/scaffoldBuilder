#ifndef SCAFFOLDGENERATOR_H // include guard
#define SCAFFOLDGENERATOR_H

#include <array>
#include <vector>
#include <string>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include "Wall/Wall.h"

class ScaffoldGenerator{
	
public:
	ScaffoldGenerator() {};
	ScaffoldGenerator(
		const std::string& meshFile, const int& nX, const int& nY, const int& nZ);
	ScaffoldGenerator(
		const double& xMin,
		const double& xMax,
		const double& yMin,
		const double& yMax,
		const double& zMin,
		const double& zMax,
		const int& nX, const int& nY, const int& nZ
	);
	ScaffoldGenerator(
		const std::vector<std::array<double, 3>>& seeds,
		const std::array<float, 6>& bounds,
		const std::array<int, 3>& blockDim
	);

	~ScaffoldGenerator() {};

	//void get_voro_cells(vtkSmartPointer<vtkPolyData> polys);
	virtual void generate_seeds() = 0;
	virtual void generate_voro(const int& regSteps) = 0;
	void generate_mesh(const double& thickness, const std::string& fileName);
	void get_seeds(std::vector<std::array<double, 3>>& outseeds);
	void process_triangles(
		std::vector<std::array<double, 3 >>& seeds,
		std::vector<std::array<double, 3>>& bCenters,
		std::vector<std::array<double, 3>>& normals);
	void get_bounds(float& xMn, float& xMx, float& yMn, float& yMx, float& zMn, float& zMx);

protected:
	// container settings for Voro++
	double xMin{ 0.0 };
	double xMax{ 0.0 };
    double yMin{ 0.0 };
	double yMax{ 0.0 };
	double zMin{ 0.0 };
    double zMax{ 0.0 };

	// block size
	int nX{ 1 }, nY{ 1 }, nZ{ 1 };
	
	// voro container
	voro::container* con;

	// mesh file name for container
	const std::string meshFile;

	// vtk scaffold mesh
	vtkSmartPointer<vtkPolyData> scaffoldMesh;

	vtkSmartPointer<vtkPolyData> meshContainer;

	// voro seeds
	std::vector<std::array<double, 3>> seeds;
};

// random without wall
class Random : public ScaffoldGenerator {
public:
	Random() {};
	~Random() {};
	Random(
		const double& xMin,
		const double& xMax,
		const double& yMin,
		const double& yMax,
		const double& zMin,
		const double& zMax,
		const int& nX, const int& nY, const int& nZ,
		const int& seedNr
	);

	void generate_seeds() override;
	void generate_voro(const int& regSteps) override;
	void get_seeds(std::vector<std::array<double, 3>>& outseeds);

private:
	int seedNr{ 100 };
	std::vector<std::array<double, 3>> seeds;
};

// random with wall
class RandomWall : public ScaffoldGenerator {
public:
	RandomWall() {};
	~RandomWall() {};
	RandomWall(const std::string fileName, const int& seedNr,
		const int& nX, const int& nY, const int& nZ, const int& neighbors);

	void RandomWall::generate_seeds() override;
	void generate_voro(const int& regSteps) override;
	void get_seeds(std::vector<std::array<double, 3>>& outseeds);

private:
	int seedNr{ 100 };
	int neighbors{ 1 };
	std::vector<std::array<double, 3>> seeds;
	/* Function that processes the triangles of the surface mesh container
	* 
	*/
};

// --------------------------------------------------------------------------------
// Volume Optimization

class VolOpt {

public:
	VolOpt() {};
	~VolOpt() {};

	// constructor overloading on for the generate points option one for the load points from file
	VolOpt(
		const Eigen::VectorXd targets,
		const Eigen::VectorXd initV, 
		const std::array<double, 6>& bounds,
		std::function<void(const std::string&)> logCallback = nullptr
	);
	void generate_random_seeds();
	void generate_random_container_seeds(vtkSmartPointer<vtkPolyData>& container);
	void loop(const int regSteps);
	void generate_mesh(
		const double& thickness,
		const std::string& fileName);
	void get_seeds(std::vector<std::array<double, 3>>& outseeds);

private:
	double xMin, xMax, yMin, yMax, zMin, zMax;
	std::vector<std::array<double, 3>> currSeeds;
	int seedNr{ 1 };
	std::array<double, 6> bounds{0.0, 10.0, 0.0, 10.0, 0.0, 10.0};
	std::function<void(const std::string&)> log_callback;
	Eigen::VectorXd targetVols;
	Eigen::VectorXd wInit;
	std::vector<vtkSmartPointer<vtkPolyData>> polys;
	vtkSmartPointer<vtkPolyData> scaffoldMesh;
};

#endif