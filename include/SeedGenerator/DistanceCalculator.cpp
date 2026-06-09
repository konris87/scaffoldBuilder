#include "DistanceCalculator.h"

MeshSDF::MeshSDF(const std::vector<Triangle>& tris, const std::vector<Vec3>& psNormals) {

	// use the triangles to build the Aabb tree
	tree = std::make_unique<AabbTree>(tris, psNormals);
};

//float MeshSDF::compute_distance(const Vec3& pt) const {
//	if (!tree) return 9999.f;
//	return tree->get_closest_distance(pt);
//};
float MeshSDF::compute_distance(const Vec3& pt) const {
    if (!tree) return 9999.f;
    
    // 1. Get the raw distance magnitude (always positive)
    float unsigned_dist = tree->get_closest_distance(pt);
    
    // 2. Use the robust raycast to determine if we are inside
    bool inside = is_inside_raycast(pt);
    
    // 3. Apply the sign (Negative = Inside, Positive = Outside)
    return inside ? -unsigned_dist : unsigned_dist;
};

//bool MeshSDF::is_inside_raycast(const Vec3& pt) const {
//    if (!tree || tree->get_root_idx() == -1) return false;
//
//    // Use a slightly irrational direction so the ray doesn't perfectly slide along flat edges
//    Vec3 dir(0.7813f, 0.5209f, 0.3421f);
//    dir.normalize();
//    Vec3 invDir(1.0f / dir.x, 1.0f / dir.y, 1.0f / dir.z);
//
//    int hitCount = 0;
//    std::vector<int> stack;
//    stack.reserve(64);
//    stack.push_back(tree->get_root_idx());
//
//    while (!stack.empty()) {
//        int nodeIdx = stack.back();
//        stack.pop_back();
//
//        const AabbNode& node = tree->get_node(nodeIdx);
//
//        // Fast AABB rejection (assume you have a standard Ray-AABB intersect function)
//        if (!ray_aabb_intersection(pt, invDir, node.bounds)) continue;
//
//        if (node.leftChild == -1 && node.rightChild == -1) {
//            // Leaf node: check the actual triangles
//            for (int triIdx : node.triangleIds) {
//                const Triangle& tri = tree->get_triangle(triIdx);
//                if (ray_triangle_intersection(pt, dir, tri.v1, tri.v2, tri.v3)) {
//                    hitCount++;
//                }
//            }
//        }
//        else {
//            stack.push_back(node.leftChild);
//            stack.push_back(node.rightChild);
//        }
//    }
//
//    // Even-Odd Rule: Odd = Inside (true), Even = Outside (false)
//    return (hitCount % 2) != 0;
//}

bool MeshSDF::is_inside_raycast(const Vec3& pt) const {
    if (!tree || tree->get_root_idx() == -1) return false;

    // Cast 3 rays in wildly different prime directions
    Vec3 dirs[3] = {
        Vec3(0.7813f, 0.5209f, 0.3421f).normalized(),
        Vec3(-0.6123f, 0.7345f, -0.2987f).normalized(),
        Vec3(0.3542f, -0.8123f, 0.4635f).normalized()
    };

    int insideVotes = 0;

    for (int i = 0; i < 3; i++) {
        Vec3 dir = dirs[i];
        Vec3 invDir(1.0f / dir.x, 1.0f / dir.y, 1.0f / dir.z);

        int hitCount = 0;
        std::vector<int> stack;
        stack.reserve(64);
        stack.push_back(tree->get_root_idx());

        while (!stack.empty()) {
            int nodeIdx = stack.back();
            stack.pop_back();

            const AabbNode& node = tree->get_node(nodeIdx);

            if (!ray_aabb_intersection(pt, invDir, node.bounds)) continue;

            if (node.leftChild == -1 && node.rightChild == -1) {
                for (int triIdx : node.triangleIds) {
                    const Triangle& tri = tree->get_triangle(triIdx);
                    if (ray_triangle_intersection(pt, dir, tri.v1, tri.v2, tri.v3)) {
                        hitCount++;
                    }
                }
            }
            else {
                stack.push_back(node.leftChild);
                stack.push_back(node.rightChild);
            }
        }

        // If this specific ray says "Inside" (Odd), add a vote
        if (hitCount % 2 != 0) {
            insideVotes++;
        }
    }

    // Majority Rules: If 2 or more rays say Inside, it is Inside.
    return insideVotes >= 2;
}