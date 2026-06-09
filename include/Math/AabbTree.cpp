#include "AabbTree.h"
#include <limits>
#include <algorithm>

AabbTree::AabbTree(
	const std::vector<Triangle>& tris, 
	const std::vector<Vec3>& psNormals) : triangles(tris), normals(psNormals) {

	// let's build the tree
	std::vector<int> initialIndices(triangles.size());
	std::iota(initialIndices.begin(), initialIndices.end(), 0);

	if (!initialIndices.empty()) {
		rootIdx = build(initialIndices);
	}
};

int AabbTree::build(std::vector<int>& indices) {

	// create the root 
	AabbNode node;
	node.bounds = { 1e9f, -1e9f, 1e9f, -1e9f, 1e9f, -1e9f };
	std::array<float, 6> centroidBounds = { 1e9f, -1e9f, 1e9f, -1e9f, 1e9f, -1e9f };

	// create the nodes for each triangle that is in the branch
	for (int idx : indices) {

		// get the triangle
		const auto& tri = triangles[idx];

		// this is the aabb of the node
		node.bounds[0] = std::min({ node.bounds[0], tri.v1.x, tri.v2.x, tri.v3.x });
		node.bounds[1] = std::max({ node.bounds[1], tri.v1.x, tri.v2.x, tri.v3.x });
		node.bounds[2] = std::min({ node.bounds[2], tri.v1.y, tri.v2.y, tri.v3.y });
		node.bounds[3] = std::max({ node.bounds[3], tri.v1.y, tri.v2.y, tri.v3.y });
		node.bounds[4] = std::min({ node.bounds[4], tri.v1.z, tri.v2.z, tri.v3.z });
		node.bounds[5] = std::max({ node.bounds[5], tri.v1.z, tri.v2.z, tri.v3.z });

		// get also the centroid
		float cx = (tri.v1.x + tri.v2.x + tri.v3.x) / 3.0f;
		float cy = (tri.v1.y + tri.v2.y + tri.v3.y) / 3.0f;
		float cz = (tri.v1.z + tri.v2.z + tri.v3.z) / 3.0f;

		// update the centroid bounds
		centroidBounds[0] = std::min(centroidBounds[0], cx);
		centroidBounds[1] = std::max(centroidBounds[1], cx);
		centroidBounds[2] = std::min(centroidBounds[2], cy);
		centroidBounds[3] = std::max(centroidBounds[3], cy);
		centroidBounds[4] = std::min(centroidBounds[4], cz);
		centroidBounds[5] = std::max(centroidBounds[5], cz);
	}

	if (indices.size() <= maxTrianglesPerLeaf) {
		// Copy the indices into the node and mark it as a leaf (children remain -1)
		node.triangleIds.assign(indices.begin(), indices.end());

		tree.push_back(node);
		return tree.size() - 1; // Return where this node lives in the array
	}

	// split the box along the longest axis
	float dx = centroidBounds[1] - centroidBounds[0];
	float dy = centroidBounds[3] - centroidBounds[2];
	float dz = centroidBounds[5] - centroidBounds[4];

	int splitAxis = 0; // default is axis X
	if (dy > dx && dy > dz) splitAxis = 1; // Y is longest
	else if (dz > dx && dz > dy) splitAxis = 2; // Z is longest

	// find the median, use nth element where the median is in the middle and everything on the left is smaller
	int midIndex = indices.size() / 2;

	std::nth_element(indices.begin(), indices.begin() + midIndex, indices.end(),
		[&](int a, int b) {
		const Triangle& tA = triangles[a];
		const Triangle& tB = triangles[b];

		float centroidA = 0.0f, centroidB = 0.0f;
		// axis X
		if (splitAxis == 0) {
			centroidA = (tA.v1.x + tA.v2.x + tA.v3.x);
			centroidB = (tB.v1.x + tB.v2.x + tB.v3.x);
		}
		// axis Y
		else if (splitAxis == 1) {
			centroidA = (tA.v1.y + tA.v2.y + tA.v3.y);
			centroidB = (tB.v1.y + tB.v2.y + tB.v3.y);
		}
		// axis Z
		else {
			centroidA = (tA.v1.z + tA.v2.z + tA.v3.z);
			centroidB = (tB.v1.z + tB.v2.z + tB.v3.z);
		}
		return centroidA < centroidB;
	}
	);
	// apply the split
	std::vector<int> leftIndices(indices.begin(), indices.begin() + midIndex);
	std::vector<int> rightIndices(indices.begin() + midIndex, indices.end());

	// assign the childs and recursively build the tree
	int leftChildIdx = build(leftIndices);
	int rightChildIdx = build(rightIndices);

	node.leftChild = leftChildIdx;
	node.rightChild = rightChildIdx;

	tree.push_back(node);
	return tree.size() - 1;
};

//@brief distance squared from a querry point to the bounds of a node
float AabbTree::distance_to_node(const Vec3& pt, const std::array<float, 6>& bounds) const {

	float dist = 0.0f;

	if (pt.x < bounds[0]) dist += (bounds[0] - pt.x) * (bounds[0] - pt.x); // pt on the left
	if (pt.x > bounds[1]) dist += (pt.x - bounds[1]) * (pt.x - bounds[1]); // on the right 

	if (pt.y < bounds[2]) dist += (bounds[2] - pt.y) * (bounds[2] - pt.y);
	if (pt.y > bounds[3]) dist += (pt.y - bounds[3]) * (pt.y - bounds[3]);
	
	if (pt.z < bounds[4]) dist += (bounds[4] - pt.z) * (bounds[4] - pt.z); // on the right 
	if (pt.z > bounds[5]) dist += (pt.z - bounds[5]) * (pt.z - bounds[5]); // on the right 

	return dist;
};

//@brief find the closest distance from the querry point to the entire tree
float AabbTree::get_closest_distance(const Vec3& pt) const {
	
	if (rootIdx == -1) return 9999.9f; // Empty tree

	float minDistSq = std::numeric_limits<float>::max();
	Vec3 bestClosestPoint;
	int bestTriangleIdx = -1;

	std::vector<int> stack;
	stack.reserve(64);
	stack.push_back(rootIdx);

	while (!stack.empty()) {

		int nodeIdx = stack.back();

		// pop out
		stack.pop_back();

		// grab the node
		AabbNode node = tree[nodeIdx];

		// check the distance from the querry point to the node box
		float dist = distance_to_node(pt, node.bounds);

		// if the distance is larger than the current best skip
		if (dist >= minDistSq) {
			continue;
		}

		// check if it is a leaf node loop for each triangle (the max contained triangles are 8)
		if (node.leftChild == -1 && node.rightChild == -1) {
			for (const auto triIdx : node.triangleIds) {
				const Triangle& t = triangles[triIdx];

				Vec3 closest = closest_triangle_point(pt, t);

				// get the difference between this point and the querry pt
				Vec3 diff = closest - pt;
				float distSq = diff.dot(diff);

				// check if the found distance is smaller than the best found so far
				if (distSq < minDistSq) {
					minDistSq = distSq;
					bestClosestPoint = closest; // best point
					bestTriangleIdx = triIdx; // best triangle
				}
			}
		}
		// if not a leaf add the childer to the stack
		else {
			float distLeft = distance_to_node(pt, tree[node.leftChild].bounds);
			float distRight = distance_to_node(pt, tree[node.rightChild].bounds);

			// push first the child with the larger distance, this will allow the closer child
			// to be popped first 
			if (distLeft < distRight) {
				stack.push_back(node.rightChild);
				stack.push_back(node.leftChild);
			}
			else {
				stack.push_back(node.leftChild);
				stack.push_back(node.rightChild);
			}
		}
	}

	// if no candidate return outside
	if (bestTriangleIdx == -1) return 9999.9f;
	return std::sqrt(minDistSq);

	//// grab the best triangle
	//const Triangle& bestTri = triangles[bestTriangleIdx];

	//Vec3 bc = get_barycentric_point(bestClosestPoint, bestTri);
	//bc.x = std::clamp(bc.x, 0.0f, 1.0f);
	//bc.y = std::clamp(bc.y, 0.0f, 1.0f);
	//bc.z = std::clamp(bc.z, 0.0f, 1.0f);
	//float bcSum = bc.x + bc.y + bc.z;
	//bc = bc / bcSum;

	//// get the pseudonormal coordinates
	//Vec3 n0 = normals[bestTri.i1];
	//Vec3 n1 = normals[bestTri.i2];
	//Vec3 n2 = normals[bestTri.i3];

	//Vec3 interpolatedNormal = (n0 * bc.x) + (n1 * bc.y) + (n2 * bc.z);	
	//if (interpolatedNormal.norm() > 1e-6f) {
	//	interpolatedNormal.normalize();
	//}
	//else {
	//	interpolatedNormal = n0; // Fallback
	//}

	//Vec3 dir = pt - bestClosestPoint;

	//// use the weighted angle criterion
	//float sign = (dir.dot(interpolatedNormal) >= 0.0f) ? 1.0f : -1.0f;

	//// return the signed distance
	//return sign * std::sqrt(minDistSq);
};
