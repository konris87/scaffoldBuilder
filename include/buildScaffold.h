#ifndef SCAFFOLD_H
#define SCAFFOLD_H

#include <string>
#include <vector>
#include <Eigen/Dense>
//#include <vtkPolyData.h>
//#include <vtkSmartPointer.h>
#include <json.hpp>


/*class buildScaffold {
public:	

	std::vector<vtkSmartPointer<vtkPolyData>> polys;

	Eigen::VectorXd targetVols;
	Eigen::VectorXd wInit;

	// constructor overloading on for the generate points option one for the load points from file
	buildScaffold(const Eigen::VectorXd targets,
		const Eigen::VectorXd initV, const std::string file, std::function<void(const std::string&)> logCallback = nullptr
	);

	void loop();

private:
	// json settings
	std::string settingsFile;

	float xMin{ 0 };
	float xMax{ 0 };
	float yMin{ 0 };
	float yMax{ 0 };
	float zMin{ 0 };
	float zMax{ 0 };
	std::array<double, 6> bounds;

	std::vector<std::array <double, 3>> currSeeds;
	int seedNr{ 0 };
	float thickness{ 0.0f };
	int regSteps{ 0 };

	// seed generation option
	std::string genOption;

	// voro stuff
	const int nX{6}, nY{ 6 }, nZ{ 6 };

	// scaffold name
	std::string version{ "" };

	void check_volume();

	void generate_seeds();

	void load_vertices();

	void load_settings();

	std::function<void(const std::string&)> log_callback;
};

*/
#endif