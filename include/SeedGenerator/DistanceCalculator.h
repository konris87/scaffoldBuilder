#ifndef DISTANCECALCULATOR_H // include guard
#define DISTANCECALCULATOR_H

#include <array>
#include <vector>
#include <Math/Vec.h>

class SDF {
public:
    virtual float compute_distance(const Vec3& point) const = 0;
    virtual ~SDF() = default;
};

class CylinderSDF : public SDF {
public:
	explicit CylinderSDF(const float radius, const float height, const Vec3& center = { 0,0,0 }) : radius(radius), height(height) {
	
		halfHeight = height * 0.5f;
	};

	~CylinderSDF() override = default;

	float compute_distance(const Vec3& pt) const override {

		// offset the point to the center of the cylinder
		Vec3 ptOffset;
		ptOffset.x = pt.x - center.x;
		ptOffset.y = pt.y - (center.y + halfHeight);
		ptOffset.z = pt.z - center.z;

		float dr = std::sqrt(ptOffset.x * ptOffset.x + ptOffset.z * ptOffset.z) - radius;
		float da = std::abs(ptOffset.y) - halfHeight;

		float dout = std::sqrt(std::max(dr, 0.0f) * std::max(dr, 0.0f) + std::max(da, 0.0f) * std::max(da, 0.0f));
		float din = std::min(std::max(dr, da), 0.0f);

		return dout + din;
	};

private:
	float radius = 0.0f;
	float height = 0.0f;
	float halfHeight = 0.0f;
	Vec3 center;
};

class PlaneSDF : public SDF {
public:
	PlaneSDF() = default;
	~PlaneSDF() override = default;

	explicit PlaneSDF(const Vec3& pointVec, const Vec3& normalVec) : normal(normalVec), point(pointVec) {
		// ensure normal
		normal = normal.normalized();
	};

	float compute_distance(const Vec3& pt) const override {
		// find diff = pt - point;
		Vec3 diff = pt - point;

		return normal.dot(diff);
	};

	Vec3 get_point() const { return point; }
	Vec3 get_normal() const { return normal; }

private:
	Vec3 normal;
	Vec3 point;
};

class BoxSDF : public SDF {
public:
	BoxSDF() = default;
	~BoxSDF() override = default;
	explicit BoxSDF(const Bounds& bounds) {
	
		// find center
		center = bounds.center;
	
		// estimate half space distances
		half.x = (bounds.xMax - bounds.xMin) * 0.5f;
		half.y = (bounds.yMax - bounds.yMin) * 0.5f;
		half.z = (bounds.zMax - bounds.zMin) * 0.5f;
	};

	float compute_distance(const Vec3& pt) const override {
		float dx = std::abs(pt.x - center.x) - half.x;
		float dy = std::abs(pt.y - center.y) - half.y;
		float dz = std::abs(pt.z - center.z) - half.z;

		// pythagorean
		float outX = std::max(dx, 0.0f);
		float outY = std::max(dy, 0.0f);
		float outZ = std::max(dz, 0.0f);
		float dout = std::sqrt(outX * outX + outY * outY + outZ * outZ);

		float din = std::min(std::max(dx, std::max(dy, dz)), 0.0f);

		return dout + din;
	}
private:
	Vec3 center{};
	Vec3 half{};
};

// just for convenience
class PointSDF : public SDF {
public:
	PointSDF() = default;
	~PointSDF() override = default;
	explicit PointSDF(const Vec3& q) : p(q) {};

	float compute_distance(const Vec3& pt) const override {
		return (p - pt).norm();
	}

	Vec3 get_point() const { return p; };
	
private:
	Vec3 p{ 0.0f, 0.0f, 0.0f };
};

//
//class PlaneDistEstimator : public DistanceEstimator {
//public:
	//PlaneDistEstimator() = default;
	//~PlaneDistEstimator() {};
//	explicit PlaneDistEstimator(
//		const std::array<double, 3>& originVec,
//		const std::array<double, 3>& normalVec) : origin(originVec.data()){
//
//		Eigen::Vector3d n(normalVec.data());
//		const float norm = n.norm();
//		if (norm == 0.0) {
//			normal.setZero();
//			dot = 0.0;
//		}
//		else {
//			normal = n / norm;
//			dot = normal.dot(Eigen::Vector3d(origin.data()));
//		}		
//	};
//	
//	float compute_distance(const std::array<double, 3>& point) const override {
//		return std::abs(normal.dot(Eigen::Vector3d(point.data())) - dot);
//	};
//
//	Vec3 get_normal() const { 
//		Vec3 vec = {
//			(float)normal.x(),
//			(float)normal.y(),
//			(float)normal.z(),
//		};
//		return vec; 
//	};
//
//	Vec3 get_origin() const { 
//		Vec3 vec = {
//			(float)origin.x(),
//			(float)origin.y(),
//			(float)origin.z(),
//		};
//		return vec;
//	};
//
//private:
//	Eigen::Vector3d normal{0.0, 0.0, 1.0};
//    Eigen::Vector3d origin;
//	float dot = { 0.0 };
//};
//
//class MeshDistEstimator : public DistanceEstimator {
//
//public:
//	MeshDistEstimator() {};
//	~MeshDistEstimator() {};
//	explicit MeshDistEstimator(
//		vtkSmartPointer<vtkPolyData> containerMesh) : container(std::move(containerMesh)) {
//		distanceCalculator = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
//		distanceCalculator->SetInput(container);
//	};
//
//	float compute_distance(const std::array<double, 3>& point) const override {
//		return distanceCalculator->EvaluateFunction((double)point[0], (double)point[1], (double)point[2]);
//	};
//
//private:
//	vtkSmartPointer<vtkPolyData> container;
//	vtkSmartPointer<vtkImplicitPolyDataDistance> distanceCalculator;
//};
//
//class ImplicitFunctionDistEstimator : public DistanceEstimator {
//public:
//	ImplicitFunctionDistEstimator() = default;
//	~ImplicitFunctionDistEstimator() override = default;
//
//	ImplicitFunctionDistEstimator(vtkSmartPointer<vtkImplicitFunction> implicitFunc)
//		: func(std::move(implicitFunc)) {
//	}
//
//	float compute_distance(const std::array<double, 3>& point) const override {
//		return std::abs(func->EvaluateFunction(point[0], point[1], point[2]));
//	}
//
//private:
//	vtkSmartPointer<vtkImplicitFunction> func;
//};
//
//class PointDistEstimator : public DistanceEstimator {
//public:
//	explicit PointDistEstimator(const std::array<double, 3>& q) : Q(q.data()) {}
//	float compute_distance(const std::array<double, 3>& p) const override {
//		Eigen::Vector3d P(p.data());
//		return (P - Q).norm();
//	}
//
//	Vec3 get_point() const {
//		Vec3 vec = {
//			(float)Q.x(),
//			(float)Q.y(),
//			(float)Q.z()
//		};
//
//		return vec;
//	};
//
//private:
//	Eigen::Vector3d Q{ 0,0,0 };
//};

#endif