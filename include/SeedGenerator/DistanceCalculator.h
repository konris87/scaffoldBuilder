#ifndef DISTANCECALCULATOR_H // include guard
#define DISTANCECALCULATOR_H

#include <array>
#include <vector>
#include <Math/Vec.h>
#include <Math/AabbTree.h>
#include "Utils/Utils.h"

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

//------------------------------------------
class MeshSDF : public SDF {
public: 
	MeshSDF() = default;
	~MeshSDF() override = default;

	MeshSDF(const std::vector<Triangle>& tris, const std::vector<Vec3>& psNormals);

	float compute_distance(const Vec3& pt) const override;

	bool is_inside_raycast(const Vec3& pt) const;

private:

	std::unique_ptr<AabbTree> tree;
};


#endif