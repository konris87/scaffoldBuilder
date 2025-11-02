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
	PlaneDistEstimator() = default;
	~PlaneDistEstimator() {};
	explicit PlaneDistEstimator(
		const std::array<double, 3>& originVec,
		const std::array<double, 3>& normalVec) : origin(originVec.data()){

		Eigen::Vector3d n(normalVec.data());
		const double norm = n.norm();
		if (norm == 0.0) {
			normal.setZero();
			dot = 0.0;
		}
		else {
			normal = n / norm;
			dot = normal.dot(Eigen::Vector3d(origin.data()));
		}		
	};
	
	double compute_distance(const std::array<double, 3>& point) const override {
		return std::abs(normal.dot(Eigen::Vector3d(point.data())) - dot);
	};

private:
	Eigen::Vector3d normal{0.0, 0.0, 1.0};
    Eigen::Vector3d origin;
	double dot = { 0.0 };
};

class MeshDistEstimator : public DistanceEstimator {

public:
	MeshDistEstimator() {};
	~MeshDistEstimator() {};
	explicit MeshDistEstimator(
		vtkSmartPointer<vtkPolyData> containerMesh) : container(std::move(containerMesh)) {
		distanceCalculator = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
		distanceCalculator->SetInput(container);
	};

	double compute_distance(const std::array<double, 3>& point) const override {
		return distanceCalculator->EvaluateFunction(point[0], point[1], point[2]);
	};

private:
	vtkSmartPointer<vtkPolyData> container;
	vtkSmartPointer<vtkImplicitPolyDataDistance> distanceCalculator;
};

class ImplicitFunctionDistEstimator : public DistanceEstimator {
public:
	ImplicitFunctionDistEstimator() = default;
	~ImplicitFunctionDistEstimator() override = default;

	ImplicitFunctionDistEstimator(vtkSmartPointer<vtkImplicitFunction> implicitFunc)
		: func(std::move(implicitFunc)) {
	}

	double compute_distance(const std::array<double, 3>& point) const override {
		return std::abs(func->EvaluateFunction(point[0], point[1], point[2]));
	}

private:
	vtkSmartPointer<vtkImplicitFunction> func;
};

class PointDistEstimator : public DistanceEstimator {
public:
	explicit PointDistEstimator(const std::array<double, 3>& q) : Q(q.data()) {}
	double compute_distance(const std::array<double, 3>& p) const override {
		Eigen::Vector3d P(p.data());
		return (P - Q).norm();
	}
private:
	Eigen::Vector3d Q{ 0,0,0 };
};

#endif