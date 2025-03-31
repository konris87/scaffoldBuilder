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
