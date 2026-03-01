#ifndef AABBTREE_H
#define AABBTREE_H

#include <vector>
#include <array>
#include "Utils/Utils.h"
#include "Vec.h"

struct AabbNode {
	std::array<float, 6> bounds = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
	int leftChild = -1;
	int rightChild = -1;
	std::vector<unsigned int> triangleIds;
};

class AabbTree {
public:

	AabbTree(const std::vector<Triangle>& tris, const std::vector<Vec3>& psNormals);
	
	float get_closest_distance(const Vec3& pt) const;

private:
	const std::vector<Triangle>& triangles;
	const std::vector<Vec3>& normals;

	// the actual tree as a vector of packed nodes
	std::vector<AabbNode> tree;

	int rootIdx = -1;
	const int maxTrianglesPerLeaf = 8;

	int build(std::vector<int>& indices);

	float distance_to_node(const Vec3& pt, const std::array<float, 6>& bounds) const;
};


#endif