#ifndef DISTANCECALCULATOR_H // include guard
#define DISTANCECALCULATOR_H

#include <array>
#include <vector>
#include <Eigen/Dense>
#include <vtkSmartPointer.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkPolyData.h>
#include <vtkBox.h>
#include <vtkCylinder.h>

class DistanceEstimator {
public:
    virtual double compute_distance(const std::array<double, 3>& point) const = 0;
    virtual ~DistanceEstimator() = default;
};

class PlaneDistEstimator : public DistanceEstimator {
public:
	PlaneDistEstimator() {};
	~PlaneDistEstimator() {};
	PlaneDistEstimator(const std::array<double, 3>& originVec, const std::array<double, 3>& normalVec) : normal(normalVec.data()), origin(originVec.data()) {
	
		normalNormSquared = normal.squaredNorm();
		
	};
	
	double compute_distance(const std::array<double, 3>& point) const override {
		Eigen::Vector3d pt(point.data());
		return std::abs((pt- origin).dot(normal)) / std::sqrt(normalNormSquared);
		
	};

private:
    Eigen::Vector3d normal;
    Eigen::Vector3d origin;
	double normalNormSquared;
};

class MeshDistEstimator : public DistanceEstimator {

public:
	MeshDistEstimator() {};
	~MeshDistEstimator() {};
	MeshDistEstimator(
		vtkSmartPointer<vtkPolyData>& containerMesh) {
		
		container = vtkSmartPointer<vtkPolyData>::New();
		container = containerMesh;

		distanceCalculator = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
		distanceCalculator->SetInput(containerMesh);
	};

	double compute_distance(const std::array<double, 3>& point) const override {
		return distanceCalculator->EvaluateFunction(point[0], point[1], point[2]);
	};

private:
	vtkSmartPointer<vtkPolyData> container;
	vtkSmartPointer<vtkImplicitPolyDataDistance> distanceCalculator;
};

//class BoxDistEstimator : public DistanceEstimator {
//
//public:
//	BoxDistEstimator() {};
//	~BoxDistEstimator() {};
//	BoxDistEstimator(
//		vtkSmartPointer<vtkBox>& box) {
//
//		container = vtkSmartPointer<vtkBox>::New();
//		container = box;
//	};
//
//	double compute_distance(const std::array<double, 3>& point) const override {
//		return std::abs(container->EvaluateFunction(point[0], point[1], point[2]));
//	};
//
//private:
//	vtkSmartPointer<vtkBox> container;
//};
//
//// Cylinder box
//class CylinderDistEstimator : public DistanceEstimator {
//
//public:
//	CylinderDistEstimator() {};
//	~CylinderDistEstimator() {};
//	CylinderDistEstimator(
//		vtkSmartPointer<vtkCylinder>& cylinder) {
//
//		container = vtkSmartPointer<vtkCylinder>::New();
//		container = cylinder;
//	};
//
//	double compute_distance(const std::array<double, 3>& point) const override {
//		return std::abs(container->EvaluateFunction(point[0], point[1], point[2]));
//	};
//
//private:
//	vtkSmartPointer<vtkCylinder> container;
//};

class ImplicitFunctionDistEstimator : public DistanceEstimator {
public:
	ImplicitFunctionDistEstimator() = default;
	~ImplicitFunctionDistEstimator() override = default;

	ImplicitFunctionDistEstimator(vtkSmartPointer<vtkImplicitFunction> implicitFunc)
		: container(implicitFunc) {
	}

	double compute_distance(const std::array<double, 3>& point) const override {
		return std::abs(container->EvaluateFunction(point[0], point[1], point[2]));
	}

private:
	vtkSmartPointer<vtkImplicitFunction> container;
};

#endif