#include "DistanceCalculator.h"

MeshSDF::MeshSDF(const std::vector<Triangle>& tris, const std::vector<Vec3>& psNormals) {

	// use the triangles to build the Aabb tree
	tree = std::make_unique<AabbTree>(tris, psNormals);
};

float MeshSDF::compute_distance(const Vec3& pt) const {
	if (!tree) return 9999.f;
	return tree->get_closest_distance(pt);
};