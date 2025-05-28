#ifndef SEEDGENERATOR_H // include guard
#define SEEDGENERATOR_H

#include <vector>
#include <array>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

class SeedGenerator{
public:
	SeedGenerator(){};
	~SeedGenerator(){};
	
	virtual void generate_seeds() = 0;
	void get_seeds(std::vector<std::array<double, 3>>& outSeeds) {
		outSeeds = seeds;
	};
	
protected:
	// voro seeds
	std::vector<std::array<double, 3>> seeds;	
};
	
class RandomGenerator : public SeedGenerator{
	
public:
	RandomGenerator(){};
	~RandomGenerator(){};
	RandomGenerator(const std::array<float, 6>& bounds, const int seedNr);
	void generate_seeds() override;
	void generate_seeds(const std::array<float, 3>& cylinderPt, const std::array<float, 3>& cylinderAxis, double cylinderRadius);

private:
	std::array<float, 6> bounds{ 0.0f, 10.0f, 0.0f, 10.0f, 0.0f, 10.0f };
	int seedNr{ 100 };
	
};

// -------------------------------------------------------------------------------
// Random inside a wall
class RandomGeneratorWall : public SeedGenerator {

public:
	RandomGeneratorWall() {};
	~RandomGeneratorWall() {};
	RandomGeneratorWall(const int seedNr, vtkSmartPointer<vtkPolyData>& container);
	void generate_seeds() override;
private:
	std::array<float, 6> bounds{ 0.0f, 10.0f, 0.0f, 10.0f, 0.0f, 10.0f };
	vtkSmartPointer<vtkPolyData> container;
	int seedNr{ 100 };
};

#endif