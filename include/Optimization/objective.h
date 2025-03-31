#ifndef OBJECTIVE_H // include guard
#define OBJECTIVE_H

#include <Eigen/Dense>
#include <Misc/cell_process.h>
#include <Wall/Wall.h>
#include <voro++.hh>

class ObjectiveFunction {

public:
	// Overload call operator and create a virtual function
	virtual double operator () (const Eigen::VectorXd& x, Eigen::VectorXd& gradient) = 0;
};

class myFunc : public ObjectiveFunction {
public:
	// Constructor
	myFunc(
		const std::array<double, 6>& bounds,
		const std::vector<std::array <double, 3>>& points,
		Eigen::VectorXd& weights,
		const Eigen::VectorXd& targetVols
	) : domBounds(bounds),
		seeds(points),
		w(weights),
		m(targetVols),
		centroids(points),
		L(Eigen::VectorXd::Zero(weights.size())),
		costs(Eigen::VectorXd::Zero(weights.size()))
		{};

	double operator()(const Eigen::VectorXd& weights, Eigen::VectorXd& gradient) override {

		// Convert weights to radii
		Eigen::VectorXd radii = convert_radii(weights);

		// create container
		voro::container_poly* conp = new voro::container_poly(
			domBounds[0], domBounds[1], domBounds[2], domBounds[3], domBounds[4], domBounds[5],
			6, 6, 6,
			false, false, false, 16
		);

		for (int i = 0; i < seeds.size(); i++) {
			conp->put(i, seeds[i][0], seeds[i][1], seeds[i][2], radii[i]);
		}

		Eigen::VectorXd currVols = Eigen::VectorXd::Zero(seeds.size());
		Eigen::VectorXd currCosts = Eigen::VectorXd::Zero(seeds.size());

		// loop in cells
		voro::c_loop_all cla(*conp);
		voro::voronoicell_neighbor cell;

		//std::cout << "Start Loop" << std::endl;
		if (cla.start()) do if (conp->compute_cell(cell, cla)) {

			// vector to store the cell vertexes in global coordinates, cell neighbors and face vertices
			std::vector<double> cellVertices;

			int seedId = cla.pid();

			// get position of seed and store it to an array
			double x = 0.0, y = 0.0, z = 0.0;
			cla.pos(x, y, z);

			// cell vertices in global system
			cell.vertices(x, y, z, cellVertices);

			// estimate volume and add it to currVolumes
			double volume{ cell.volume() };
			//std::cout << "Cell Volume: " << volume << std::endl;
			currVols(seedId) = volume;

			// estimate costs
			double cost = transport(cellVertices, std::vector<double> {x, y, z});
			currCosts(seedId) = cost;

			double cx, cy, cz;
			cell.centroid(cx, cy, cz);

			//std::cout << cx + x << " " << cy + y << " " << cz + z << std::endl;
			//std::cout << "------------------" << std::endl;

			std::array<double, 3> centroid{ cx + x, cy + y, cz + z };

			centroids[seedId] = centroid;

		} while (cla.inc());

		delete conp;

		//std::cout << weights << std::endl;

		//std::cout << "minimum volume: " << currVols.minCoeff() << std::endl;
		//std::cout << "sum of costs: " << currCosts.sum() << std::endl;

		//// dot product of volume difference with weights
		gradient = currVols - m;

		// function value
		double f = (currVols-m).dot(weights) - currCosts.sum();

		// update volume Error
		volError = ((currVols - m).array() / m.array()).abs().maxCoeff() * 100;

		//std::cout << "Vol Error " << volError << std::endl;
		w = weights;

		return f;
	}

	// volume error
	double volError{ 0 };
	// weights
	Eigen::VectorXd w;
	// array to store centroids for each cell
	std::vector<std::array<double, 3>> centroids;


private:
	// points
	const std::vector < std::array<double, 3>> seeds;
	// target volumes
	const Eigen::VectorXd m;
	// current volumes
	Eigen::VectorXd L;
	// transport costs
	Eigen::VectorXd costs;
	// domain bounds for pyvoro
	const std::array<double, 6> domBounds;

};

class myFuncWall : public ObjectiveFunction {
public:
	// Constructor
	myFuncWall(
		MeshWall& wallObj,
		const std::array<double, 6>& bounds,
		const std::vector<std::array <double, 3>>& points,
		Eigen::VectorXd& weights,
		const Eigen::VectorXd& targetVols) : 
		meshWall(wallObj),
		domBounds(bounds),
		seeds(points),
		w(weights),
		m(targetVols),
		centroids(points),
		L(Eigen::VectorXd::Zero(weights.size())),
		costs(Eigen::VectorXd::Zero(weights.size()))
	{};

	double operator()(const Eigen::VectorXd& weights, Eigen::VectorXd& gradient) override {

		// Convert weights to radii
		Eigen::VectorXd radii = convert_radii(weights);

		// create container
		voro::container_poly* conp = new voro::container_poly(
			domBounds[0], domBounds[1], domBounds[2], domBounds[3], domBounds[4], domBounds[5],
			10, 10, 10,
			false, false, false, 16
		);
		conp->add_wall(meshWall);

		for (int i = 0; i < seeds.size(); i++) {
			conp->put(i, seeds[i][0], seeds[i][1], seeds[i][2], radii[i]);
		}

		Eigen::VectorXd currVols = Eigen::VectorXd::Zero(seeds.size());
		Eigen::VectorXd currCosts = Eigen::VectorXd::Zero(seeds.size());

		// loop in cells
		voro::c_loop_all cla(*conp);
		voro::voronoicell_neighbor cell;

		//std::cout << "Start Loop" << std::endl;
		if (cla.start()) do if (conp->compute_cell(cell, cla)) {

			// vector to store the cell vertexes in global coordinates, cell neighbors and face vertices
			std::vector<double> cellVertices;

			int seedId = cla.pid();

			// get position of seed and store it to an array
			double x = 0.0, y = 0.0, z = 0.0;
			cla.pos(x, y, z);

			// cell vertices in global system
			cell.vertices(x, y, z, cellVertices);

			// estimate volume and add it to currVolumes
			double volume{ cell.volume() };
			//std::cout << "Cell Volume: " << volume << std::endl;
			currVols(seedId) = volume;

			// estimate costs
			double cost = transport(cellVertices, std::vector<double> {x, y, z});
			currCosts(seedId) = cost;

			double cx, cy, cz;
			cell.centroid(cx, cy, cz);

			//std::cout << cx + x << " " << cy + y << " " << cz + z << std::endl;
			//std::cout << "------------------" << std::endl;

			std::array<double, 3> centroid{ cx + x, cy + y, cz + z };

			centroids[seedId] = centroid;

		} while (cla.inc());

		delete conp;

		//std::cout << weights << std::endl;

		//std::cout << "minimum volume: " << currVols.minCoeff() << std::endl;
		//std::cout << "sum of costs: " << currCosts.sum() << std::endl;

		//// dot product of volume difference with weights
		gradient = currVols - m;

		// function value
		double f = (currVols - m).dot(weights) - currCosts.sum();

		// update volume Error
		volError = ((currVols - m).array() / m.array()).abs().maxCoeff() * 100;

		//std::cout << "Vol Error " << volError << std::endl;
		w = weights;

		return f;
	}

	// volume error
	double volError{ 0 };
	// weights
	Eigen::VectorXd w;
	// array to store centroids for each cell
	std::vector<std::array<double, 3>> centroids;


private:

	MeshWall& meshWall;
	// points
	const std::vector < std::array<double, 3>> seeds;
	// target volumes
	const Eigen::VectorXd m;
	// current volumes
	Eigen::VectorXd L;
	// transport costs
	Eigen::VectorXd costs;
	// domain bounds for pyvoro
	std::array<double, 6> domBounds;

};

#endif