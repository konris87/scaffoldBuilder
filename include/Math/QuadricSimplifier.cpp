#include "QuadricSimplifier.h"
#include "Vec.h"
#include "Eigen/Dense"
#include <queue>
#include <array>
#include <map>
#include <unordered_set>

using namespace mesh_simplify;

Info mesh_simplify::simplify(
    std::vector<float>& vertices,
    std::vector<unsigned int>& indices,
    std::vector<float>& normals,
    const Options& opts){

        // create positions
        std::vector<Vec3> positions(vertices.size() / 3);
        for (size_t i{0}; i < positions.size(); i++){
            positions[i] = {
                vertices[3 * i],
                vertices[3 * i + 1],
                vertices[3 * i + 2]
            };
        }

        // create triangles as unsigned triads from indices
        std::vector<std::array<unsigned,3>> tris(indices.size() / 3);
        for (size_t i{0}; i < tris.size(); i++){
            tris[i] = {
                indices[3 * i],
                indices[3 * i + 1],
                indices[3 * i + 2]
            };
        }

        // now proceed with the Q metrics for each vertex
        // each vertex has an Σkp
        std::vector<Eigen::Matrix4d> Q(positions.size(), Eigen::Matrix4d::Zero());

        // loop for each triangle. It defines a plane that we are going to use to build the q metric
        for(const auto& tri: tris){

            Vec3 p1 = positions[tri[0]];
            Vec3 p2 = positions[tri[1]];
            Vec3 p3 = positions[tri[2]];

            // weight by area
            Vec3 cr   = (p2 - p1).cross(p3 - p1);
            double area = 0.5 * cr.norm();
            Vec3 normal = cr.normalized();

            // find constant d
            double d = -normal.dot(p1);

            Eigen::Vector4d p{normal.x, normal.y, normal.z, d};

            // find kp weighted by area and add to each vertex
            Eigen::Matrix4d kp = area * (p * p.transpose());

            Q[tri[0]] += kp;
            Q[tri[1]] += kp;
            Q[tri[2]] += kp;
        }

        // valid pairs = mesh edges (undirected, deduplicated). The same map
        // also gives us boundary detection: an edge used by a single triangle
        // is on an open boundary.
        auto makeEdge = [](unsigned int a, unsigned int b){
            return a < b ? std::make_pair(a,b) : std::make_pair(b,a);
        };

        std::unordered_map<std::pair<unsigned int, unsigned int>, int, hash_pair> edgeSet;

        for (const auto& tri : tris){
            for(int e{0}; e < 3; e++){
                edgeSet[makeEdge(tri[e], tri[(e + 1) % 3])]++;
            }
        }

        // hard-lock boundary vertices: any edge with a single incident triangle
        std::vector<char> locked(positions.size(), 0);
        if (opts.lockBoundary){
            for (auto& [edge, c] : edgeSet){
                if (c == 1) {
                    locked[edge.first]  = 1;
                    locked[edge.second] = 1;
                }
            }
        }

        // min-heap of collapse candidates ordered by cost
        std::priority_queue<Edge, std::vector<Edge>, EdgeComparator> minq;

        // per-vertex version stamps for lazy invalidation of stale heap entries
        std::vector<uint32_t> version(positions.size(), 0);

        // Score edge (v0,v1): optimal contraction point + its quadric cost.
        // Shared by the initial heap build and the post-collapse re-push.
        auto scoreEdge = [&](unsigned v0, unsigned v1) -> Edge {
            Eigen::Matrix4d Qbar = Q[v0] + Q[v1];
            Eigen::Matrix3d A = Qbar.topLeftCorner<3,3>();
            Eigen::Vector3d b = Qbar.block<3,1>(0,3);

            Eigen::Vector4d vbar;
            if (std::abs(A.determinant()) > 1e-9) {
                vbar << A.ldlt().solve(-b), 1.0;          // optimal point
            } else {
                Eigen::Vector3d p0(positions[v0].x, positions[v0].y, positions[v0].z);
                Eigen::Vector3d p1(positions[v1].x, positions[v1].y, positions[v1].z);
                vbar << 0.5*(p0 + p1), 1.0;               // fallback: midpoint
            }

            double cost = (vbar.transpose() * Qbar * vbar).value();
            if (cost < 0) cost = 0;

            Vec3 target{ (float)vbar.x(), (float)vbar.y(), (float)vbar.z() };
            return Edge{ cost, v0, v1, target, version[v0], version[v1] };
        };

        // build the initial heap: one candidate per non-locked edge
        for (auto& [edge, c] : edgeSet) {
            unsigned v0 = edge.first, v1 = edge.second;
            if (locked[v0] || locked[v1]) continue;
            minq.push(scoreEdge(v0, v1));
        }

        // per-vertex / per-triangle liveness
        std::vector<char> vertexAlive(positions.size(), 1);
        std::vector<char> triAlive(tris.size(), 1);

        // vertex -> incident triangle ids
        std::vector<std::vector<unsigned>> vtris(positions.size());
        for (unsigned t = 0; t < tris.size(); ++t)
            for (unsigned k = 0; k < 3; ++k)
                vtris[tris[t][k]].push_back(t);

        // Would collapsing (vKeep,vDrop) onto `tgt` invert any incident triangle?
        auto wouldFlip = [&](unsigned vKeep, unsigned vDrop, const Vec3& tgt) -> bool {
            for (unsigned side : { vKeep, vDrop }) {
                for (unsigned t : vtris[side]) {
                    if (!triAlive[t]) continue;
                    const auto& tri = tris[t];

                    // triangles containing BOTH endpoints vanish on collapse -> skip
                    bool hasKeep = false, hasDrop = false;
                    for (unsigned k = 0; k < 3; ++k) {
                        if (tri[k] == vKeep) hasKeep = true;
                        if (tri[k] == vDrop) hasDrop = true;
                    }
                    if (hasKeep && hasDrop) continue;

                    // corner positions before vs. after the move
                    Vec3 A0[3], A1[3];
                    for (unsigned k = 0; k < 3; ++k) {
                        A0[k] = positions[tri[k]];
                        A1[k] = (tri[k] == vKeep || tri[k] == vDrop)
                            ? tgt : positions[tri[k]];
                    }

                    Vec3 nOld = (A0[1]-A0[0]).cross(A0[2]-A0[0]);
                    Vec3 nNew = (A1[1]-A1[0]).cross(A1[2]-A1[0]);

                    if (nOld.dot(nNew) < 0.0f) return true;   // orientation inverted
                }
            }
            return false;
        };

        // stop when we drop below the target triangle count
        size_t target = opts.targetTriangles > 0
              ? opts.targetTriangles
              : size_t(opts.targetRatio * tris.size());

        size_t liveTris       = triAlive.size();
        size_t collapses      = 0;
        double maxCollapseErr = 0.0;

        while (target < liveTris && !minq.empty()){

            // get the edge with the minimum cost and pop it
            Edge e = minq.top();
            minq.pop();

            // stale entry? versions moved on since it was pushed
            if (e.version0 != version[e.v0] || e.version1 != version[e.v1]) continue;

            // either endpoint already dead
            if (!vertexAlive[e.v0] || !vertexAlive[e.v1]) continue;

            // cheapest remaining collapse exceeds the budget -> nothing better left
            if (e.cost > opts.maxError) break;

            unsigned vKeep = e.v0, vDrop = e.v1;

            // reject collapses that would invert an incident triangle
            if (opts.preventFlips && wouldFlip(vKeep, vDrop, e.targetPoint)) continue;

            // record statistics for this (accepted) collapse
            ++collapses;
            if (e.cost > maxCollapseErr) maxCollapseErr = e.cost;

            // contract vDrop into vKeep
            positions[vKeep] = e.targetPoint;   // move survivor to optimal point
            Q[vKeep] += Q[vDrop];               // survivor inherits both quadrics
            vertexAlive[vDrop] = 0;

            // rewire every triangle that referenced vDrop
            for (unsigned t : vtris[vDrop]){
                if (!triAlive[t]) continue;

                for (unsigned k = 0; k < 3; ++k)
                    if (tris[t][k] == vDrop) tris[t][k] = vKeep;

                bool degenerate =
                    tris[t][0]==tris[t][1] ||
                    tris[t][1]==tris[t][2] ||
                    tris[t][0]==tris[t][2];

                if (degenerate) { triAlive[t] = 0; --liveTris; }
                else            { vtris[vKeep].push_back(t); }
            }
            vtris[vDrop].clear();

            version[vKeep]++;   // invalidates old heap entries touching vKeep
            version[vDrop]++;

            // re-push updated edges around vKeep
            std::unordered_set<unsigned> nbrs;
            for (unsigned t : vtris[vKeep]) {
                if (!triAlive[t]) continue;
                for (unsigned k = 0; k < 3; ++k)
                    if (tris[t][k] != vKeep) nbrs.insert(tris[t][k]);
            }
            for (unsigned n : nbrs) {
                if (!vertexAlive[n] || locked[n]) continue;
                minq.push(scoreEdge(vKeep, n));
            }
        }

        // ---------------------------------------------------------------
        // Compaction: keep only vertices still referenced by a live
        // triangle, remap indices into a contiguous range.
        // ---------------------------------------------------------------
        std::vector<char> used(positions.size(), 0);
        for (unsigned t = 0; t < tris.size(); ++t) {
            if (!triAlive[t]) continue;
            used[tris[t][0]] = used[tris[t][1]] = used[tris[t][2]] = 1;
        }

        std::vector<unsigned> remap(positions.size(), 0xFFFFFFFFu);
        std::vector<float> outVertices;
        outVertices.reserve(positions.size() * 3);

        unsigned newIndex = 0;
        for (unsigned v = 0; v < positions.size(); ++v) {
            if (!used[v]) continue;
            remap[v] = newIndex++;
            outVertices.push_back(positions[v].x);
            outVertices.push_back(positions[v].y);
            outVertices.push_back(positions[v].z);
        }

        std::vector<unsigned int> outIndices;
        outIndices.reserve(liveTris * 3);
        for (unsigned t = 0; t < tris.size(); ++t) {
            if (!triAlive[t]) continue;
            outIndices.push_back(remap[tris[t][0]]);
            outIndices.push_back(remap[tris[t][1]]);
            outIndices.push_back(remap[tris[t][2]]);
        }

        // ---------------------------------------------------------------
        // Recompute area-weighted vertex normals on the simplified mesh.
        // ---------------------------------------------------------------
        std::vector<float> outNormals;
        if (opts.recomputeNormals) {
            std::vector<Vec3> vn(newIndex, Vec3(0.0f, 0.0f, 0.0f));
            for (size_t i = 0; i < outIndices.size(); i += 3) {
                unsigned i0 = outIndices[i], i1 = outIndices[i+1], i2 = outIndices[i+2];
                Vec3 a{ outVertices[3*i0], outVertices[3*i0+1], outVertices[3*i0+2] };
                Vec3 b{ outVertices[3*i1], outVertices[3*i1+1], outVertices[3*i1+2] };
                Vec3 c{ outVertices[3*i2], outVertices[3*i2+1], outVertices[3*i2+2] };
                Vec3 fn = (b - a).cross(c - a);   // unnormalized => area weighted
                vn[i0] += fn; vn[i1] += fn; vn[i2] += fn;
            }
            outNormals.resize(vn.size() * 3);
            for (size_t i = 0; i < vn.size(); ++i) {
                Vec3 n = vn[i].norm() > 1e-12f ? vn[i].normalized()
                                               : Vec3(0.0f, 0.0f, 1.0f);
                outNormals[3*i]   = n.x;
                outNormals[3*i+1] = n.y;
                outNormals[3*i+2] = n.z;
            }
        }

        // ---------------------------------------------------------------
        // Fill statistics, then hand the new buffers back to the caller.
        // ---------------------------------------------------------------
        Info information;
        information.inputVertices    = positions.size();
        information.inputTriangles   = tris.size();
        information.outputVertices   = newIndex;
        information.outputTriangles  = liveTris;
        information.collapses        = collapses;
        information.maxCollapseError = maxCollapseErr;

        vertices.swap(outVertices);
        indices.swap(outIndices);
        if (opts.recomputeNormals) normals.swap(outNormals);

        return information;
    }
