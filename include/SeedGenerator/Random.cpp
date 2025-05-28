#include "Random.h"
#include "Utils/Utils.h"
#include <random>


// Simple Random Generator
RandomGenerator::RandomGenerator(
	const std::array<float, 6>& bounds,
	const int seedNr) : bounds(bounds), seedNr(seedNr) {};

void RandomGenerator::generate_seeds() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> disX(bounds[0] + 0.1, bounds[1] - 0.1);
	std::uniform_real_distribution<> disY(bounds[2] + 0.1, bounds[3] - 0.1);
	std::uniform_real_distribution<> disZ(bounds[4] + 0.1, bounds[5] - 0.1);

	for (int i{ 0 }; i < seedNr; i++) {

		double x = disX(gen);
		double y = disY(gen);
		double z = disZ(gen);

		seeds.push_back({ x, y, z });
	}
};

void RandomGenerator::generate_seeds(
	const std::array<float, 3>& cylinderPt, const std::array<float, 3>& cylinderAxis, double cylinderRadius) {

	std::vector<Eigen::Vector3d> points;
	points.reserve(seedNr);

	double height = std::abs(bounds[5] - bounds[4]);

	// ensure the axis is normalized
	Eigen::Vector3d axis = { cylinderAxis[0], cylinderAxis[1], cylinderAxis[2] };
	axis = axis.normalized();

	// Create orthonormal basis
	Eigen::Vector3d x;
	if (std::abs(axis.x()) < 0.99)
		x = axis.cross(Eigen::Vector3d::UnitX()).normalized();
	else
		x = axis.cross(Eigen::Vector3d::UnitY()).normalized();

	Eigen::Vector3d y = axis.cross(x);

		std::mt19937 gen(std::random_device{}());
	std::uniform_real_distribution<double> distAngle(0, 2 * M_PI);
	std::uniform_real_distribution<double> distHeight(0, height);
	std::uniform_real_distribution<double> distRadius(0, 1);

	Eigen::Vector3d startPoint = { cylinderPt[0], cylinderPt[1], cylinderPt[2] };

	for (int i = 0; i < seedNr; ++i) {
		double theta = distAngle(gen);
		double u = distRadius(gen);
		double r_prime = std::sqrt(u) * cylinderRadius;
		double h = distHeight(gen);

		Eigen::Vector3d localPoint =
			r_prime * std::cos(theta) * x +
			r_prime * std::sin(theta) * y +
			h * axis;

		Eigen::Vector3d finalPt = startPoint + localPoint;
		seeds.push_back({ finalPt.x(), finalPt.y(), finalPt.z()});
	}

};


// --------------------------------------------------------------------------------------------------
// Simple Random Generator Inside Mesh
RandomGeneratorWall::RandomGeneratorWall(
	const int seedNr, vtkSmartPointer<vtkPolyData>& container) : seedNr(seedNr), container(container)
{
	double bs[6];

	container->GetBounds(bs);

	bounds[0] = bs[0];
	bounds[1] = bs[1];
	bounds[2] = bs[2];
	bounds[3] = bs[3];
	bounds[4] = bs[4];
	bounds[5] = bs[5];
	
};

void RandomGeneratorWall::generate_seeds() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> disX(bounds[0] + 0.1, bounds[1] - 0.1);
	std::uniform_real_distribution<> disY(bounds[2] + 0.1, bounds[3] - 0.1);
	std::uniform_real_distribution<> disZ(bounds[4] + 0.1, bounds[5] - 0.1);

	int idx{ 0 };
	while (idx < seedNr) {
		double x = disX(gen);
		double y = disY(gen);
		double z = disZ(gen);
		Eigen::Vector3d pt{ x, y, z };
		Eigen::Vector3d dir = { 1.0, 1.0, 1.0 };

		if (is_inside_mesh(container, pt)) {
			idx++;
			std::array<double, 3> t{ x, y, z };
			seeds.push_back(t);
		}
	}
};
