#ifndef QUADRIC_SIMPLIFIER_H
#define QUADRIC_SIMPLIFIER_H

#include <limits>
#include <vector>
#include <unordered_map>
#include "Vec.h"

namespace mesh_simplify{

    // collapsing options
    struct Options {
        size_t targetTriangles = 0; // target nr of triangles. if > 0 it is the default option
        float targetRatio = 0.25f; // percentage of maintained triangles

        // early exit if the cost of collapsing a face is larger than this
        float maxError = std::numeric_limits<float>::max();

        // prevent collapsing a face that inverts an incident triangle
        bool  preventFlips = true;

        // recompute area-weighted vertex normals after simplification.
        bool  recomputeNormals = true;

        // hard lock boundary edges near container
        bool lockBoundary = true;
    };

    // collapsing info
    struct Info {
        std::size_t inputVertices   = 0;
        std::size_t inputTriangles  = 0;
        std::size_t outputVertices  = 0;
        std::size_t outputTriangles = 0;
        std::size_t collapses       = 0;
        double      maxCollapseError = 0.0;
    };

    Info simplify(
        std::vector<float>& vertices,
        std::vector<unsigned int>& indices,
        std::vector<float>& normals,
        const Options& opts);

    // A hash function used to hash a pair of any kind
    struct hash_pair {
        size_t operator()(const std::pair<unsigned,unsigned>& p) const {
            return (size_t(p.first) << 32) ^ size_t(p.second);
        }
    };

    // create a struct that for each edge holds the information we want
    struct Edge {
        double cost;
        unsigned int v0; // position of v1
        unsigned int v1; // position of v2
        Vec3 targetPoint; // collapsing point
        
        uint32_t version0; // vertex version at insertion
        uint32_t version1;
    };

    // create also a comparator
    struct EdgeComparator {
        bool operator()(const Edge& edge1, const Edge& edge2) const {
            return edge1.cost > edge2.cost;
        };
    };
}

#endif // QUADRIC_SIMPLIFIER_H
