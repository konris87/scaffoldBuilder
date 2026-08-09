#include "GeneratorLewiner.h"
#include "LookUpTable.h"
#include "Math/Vec.h"
#include "Math/Kdtree.h"
#include "Misc/Imgui_Stdlib.h"
#include <chrono>
#include <memory>
#include <random>
#include <iostream>
#include <map>
#include <algorithm>
#include <limits>
#include <omp.h>
#include <Eigen/Dense>

using namespace std::chrono_literals;

GeneratorLewiner::GeneratorLewiner(
	const std::vector<Vec3>& seeds,
	const std::array<float, 6>& bounds,
	const std::array<int, 3>& dims,
	Logger* uiLogger,
	const float threshold,
	const float isoLevel,
	const bool renderMode
	) : seeds(seeds), bounds(bounds), blockDims(dims), logger(uiLogger), threshold(threshold), isoLevel(isoLevel), renderMode(renderMode) {

	// setup opengl
	if (renderMode) {
		_setup_mesh();

		_setup_edges();
	}

	isLoadedFromFile = false;
};

GeneratorLewiner::GeneratorLewiner(
	const std::string fileName,
	Logger* uiLogger,
	const bool renderMode
) : logger(uiLogger) {

	std::ifstream file(fileName, std::ios::in | std::ios::binary);
	std::vector<openstl::Triangle> meshTris = openstl::deserializeStl(file);
	file.close();

	// Deduplicate vertices & get faces
	auto [meshVs, meshFs] = convertToVerticesAndFaces(meshTris);

	// populate vertices, indices and faces
	meshVertices.clear();
	meshVertices.reserve(meshVs.size());
	meshTriangles.clear();
	meshTriangles.reserve(meshFs.size());

	// pass the vertices and triangles
	for (const auto& v : meshVs) {
		LVertex vertex;
		vertex.x = v.x;
		vertex.y = v.y;
		vertex.z = v.z;
		vertex.nx = 0.0f;
		vertex.ny = 0.0f;
		vertex.nz = 0.0f;

		meshVertices.push_back(vertex);
	}

	for (const auto& f : meshFs) {
		LTriangle triangle;
		triangle.v1 = f[0];
		triangle.v2 = f[1];
		triangle.v3 = f[2];

		// add also the normal
		Vec3 a = { meshVertices[triangle.v1].x, meshVertices[triangle.v1].y, meshVertices[triangle.v1].z };
		Vec3 b = { meshVertices[triangle.v2].x, meshVertices[triangle.v2].y, meshVertices[triangle.v2].z };
		Vec3 c = { meshVertices[triangle.v3].x, meshVertices[triangle.v3].y, meshVertices[triangle.v3].z };

		Vec3 normal = (b - a).cross(c - a).normalized();

		triangle.normal = normal;

		// add to all three vertices the face normal to average
		meshVertices[triangle.v1].nx += normal.x;
		meshVertices[triangle.v1].ny += normal.y;
		meshVertices[triangle.v1].nz += normal.z;
		meshVertices[triangle.v2].nx += normal.x;
		meshVertices[triangle.v2].ny += normal.y;
		meshVertices[triangle.v2].nz += normal.z;
		meshVertices[triangle.v3].nx += normal.x;
		meshVertices[triangle.v3].ny += normal.y;
		meshVertices[triangle.v3].nz += normal.z;

		meshTriangles.push_back(triangle);
	}

	// normalize the vertex normals
	for (auto& v : meshVertices) {
		Vec3 normal = { v.nx, v.ny, v.nz };
		float len = normal.norm();

		if (len > 1e-6f) {
			normal = normal.normalized();
			v.nx = normal.x;
			v.ny = normal.y;
			v.nz = normal.z;
		}
		else {
			// Fallback for isolated vertices
			v.nx = 0.0f;
			v.ny = 1.0f;
			v.nz = 0.0f;
		}
	}

	if (renderMode) {
		_setup_mesh();

		_setup_edges();

		// update opengl objects
		update_render();
	}
	// update axis aligned bounding box
	_update_bounding_box();

	// for a scaffold loaded from STL the mesh AABB is the only geometry we
	// have, so it intentionally also becomes the grid/domain definition
	bounds[0] = aabb.pMin.x;
	bounds[1] = aabb.pMax.x;
	bounds[2] = aabb.pMin.y;
	bounds[3] = aabb.pMax.y;
	bounds[4] = aabb.pMin.z;
	bounds[5] = aabb.pMax.z;

	isLoadedFromFile = true;

};

/**
 * @brief Classifies every voxel of the generation grid as definitely OUTSIDE
 * the container, definitely INSIDE it, or NeedsExact (close enough to the
 * surface, within 'margin', that it must be evaluated exactly).
 *
 * Why this is exact and not a heuristic:
 * A signed distance function is 1-Lipschitz: |SDF(p) - SDF(q)| <= |p - q|
 * for any two points p, q. Consequently, if the EXACT distance is known at
 * every corner of an axis-aligned box, and the smallest of those corner
 * distances exceeds (margin + boxDiagonal), then no point inside that box -
 * however close to a corner - can have a distance smaller than 'margin'
 * (it would require a smaller-than-boxDiagonal step to close that much of a
 * gap). The mirror argument holds for proving a box lies entirely more than
 * 'margin' inside.
 *
 * This lets us evaluate the real (and, for a mesh container, expensive:
 * BVH nearest-point query + raycast sign test) SDF only on a coarse
 * sub-sampling of the grid (every 'stride' voxels), and use that coarse
 * sampling to safely resolve the classification of every coarse CELL's
 * interior in O(1), without ever calling the real SDF on those interior
 * (fine) voxels. Only voxels in cells that remain ambiguous - i.e. within
 * 'margin' of the surface - fall back to per-voxel exact evaluation.
 *
 * @param con     the container whose distance field is being classified
 * @param margin  the same distance threshold used by compute_scalar_field's
 *                "more than margin outside, don't bother" early-out; must
 *                match it for the classification to be consistent with how
 *                the caller actually uses the result
 */
std::vector<GeneratorLewiner::NarrowBandClass> GeneratorLewiner::classify_container_narrow_band(
	const IContainer& con, float margin) {

	size_t totalVoxels = static_cast<size_t>(blockDims[0]) *
		static_cast<size_t>(blockDims[1]) *
		static_cast<size_t>(blockDims[2]);

	// default to the safe fallback everywhere; only cells we can PROVE are
	// far enough from the surface get reclassified below
	std::vector<NarrowBandClass> classification(totalVoxels, NarrowBandClass::NeedsExact);

	// coarse sampling stride, in fine voxels. Larger = cheaper coarse pass
	// but a thicker band of "ambiguous" cells around the surface (since the
	// safety margin below grows with the coarse cell's diagonal).
	const int stride = 8;

	auto ceil_div = [](int a, int b) { return (a + b - 1) / b; };

	std::array<int, 3> coarseDims = {
		ceil_div(blockDims[0] - 1, stride) + 1,
		ceil_div(blockDims[1] - 1, stride) + 1,
		ceil_div(blockDims[2] - 1, stride) + 1
	};

	auto coarse_idx = [&](int ci, int cj, int ck) -> size_t {
		return static_cast<size_t>(ci) +
			static_cast<size_t>(cj) * static_cast<size_t>(coarseDims[0]) +
			static_cast<size_t>(ck) * static_cast<size_t>(coarseDims[0]) * static_cast<size_t>(coarseDims[1]);
	};

	// maps a coarse-grid index back to the fine voxel index it samples,
	// clamped so the last coarse layer always lands exactly on the grid edge
	auto coarse_to_fine = [&](int c, int dim) {
		return std::min(c * stride, dim - 1);
	};

	std::vector<float> coarseDist(
		static_cast<size_t>(coarseDims[0]) * coarseDims[1] * coarseDims[2]);

	// 1. Evaluate the EXACT container distance, but only on the coarse grid.
	// This is the only place we still pay the real (possibly expensive) SDF
	// cost up front - on roughly 1/stride^3 of the points the fine grid has.
	#pragma omp parallel for collapse(3)
	for (int ci = 0; ci < coarseDims[0]; ci++) {
		for (int cj = 0; cj < coarseDims[1]; cj++) {
			for (int ck = 0; ck < coarseDims[2]; ck++) {

				int i = coarse_to_fine(ci, blockDims[0]);
				int j = coarse_to_fine(cj, blockDims[1]);
				int k = coarse_to_fine(ck, blockDims[2]);

				Vec3 point(
					bounds[0] + i * stepX,
					bounds[2] + j * stepY,
					bounds[4] + k * stepZ
				);

				coarseDist[coarse_idx(ci, cj, ck)] = con.sdf->compute_distance(point);
			}
		}
	}

	// 2. Classify every coarse CELL (the 8 coarse samples at its corners)
	// using the Lipschitz bound described above, then stamp that
	// classification onto every fine voxel inside the cell.
	const float cellDiagonal = std::sqrt(
		(stride * stepX) * (stride * stepX) +
		(stride * stepY) * (stride * stepY) +
		(stride * stepZ) * (stride * stepZ)
	);
	const float safeOutsideThreshold = margin + cellDiagonal;
	const float safeInsideThreshold = -(margin + cellDiagonal);

	#pragma omp parallel for collapse(3)
	for (int ci = 0; ci < coarseDims[0] - 1; ci++) {
		for (int cj = 0; cj < coarseDims[1] - 1; cj++) {
			for (int ck = 0; ck < coarseDims[2] - 1; ck++) {

				float cellMin = std::numeric_limits<float>::max();
				float cellMax = std::numeric_limits<float>::lowest();

				for (int dc = 0; dc < 2; dc++) {
					for (int dj = 0; dj < 2; dj++) {
						for (int dk = 0; dk < 2; dk++) {
							float d = coarseDist[coarse_idx(ci + dc, cj + dj, ck + dk)];
							cellMin = std::min(cellMin, d);
							cellMax = std::max(cellMax, d);
						}
					}
				}

				NarrowBandClass cls = NarrowBandClass::NeedsExact;
				if (cellMin > safeOutsideThreshold) {
					cls = NarrowBandClass::Outside;
				}
				else if (cellMax < safeInsideThreshold) {
					cls = NarrowBandClass::Inside;
				}

				// leaving ambiguous cells at their NeedsExact default avoids
				// writing anything for the common (near-surface) case
				if (cls == NarrowBandClass::NeedsExact) {
					continue;
				}

				int iStart = coarse_to_fine(ci, blockDims[0]);
				int iEnd = coarse_to_fine(ci + 1, blockDims[0]);
				int jStart = coarse_to_fine(cj, blockDims[1]);
				int jEnd = coarse_to_fine(cj + 1, blockDims[1]);
				int kStart = coarse_to_fine(ck, blockDims[2]);
				int kEnd = coarse_to_fine(ck + 1, blockDims[2]);

				for (int i = iStart; i <= iEnd; i++) {
					for (int j = jStart; j <= jEnd; j++) {
						for (int k = kStart; k <= kEnd; k++) {
							// benign race on shared cell-boundary voxels: every
							// write here is independently a mathematically
							// valid classification for that exact voxel (see
							// the Lipschitz argument above), so whichever
							// neighbouring cell's write wins is still correct
							classification[find_vertex_index(i, j, k)] = cls;
						}
					}
				}
			}
		}
	}

	return classification;
}

/**
 * @brief Fills scalarField with the signed lattice density value at every
 * voxel of the generation grid, ready for marching_cubes() to extract an
 * isosurface from. This is the dominant cost of scaffold generation, so its
 * structure (and the optimizations below) directly determine how long
 * generation takes.
 *
 * Per voxel, two distance queries normally compete for cost:
 *  1. distance to the container surface (con.sdf), used to know whether the
 *     voxel is inside the container at all, and (for "distance from
 *     container" varied thickness) how close it is to the wall;
 *  2. distance to the nearest seed points (the kdtree query), which defines
 *     the actual porous lattice value and is unavoidable for every voxel
 *     that lies inside the container, since the lattice fills the interior.
 *
 * For a mesh container, (1) is by far the more expensive of the two - it is
 * a BVH nearest-point query plus a multi-ray raycast just to get the sign
 * (inside/outside). Three optimizations target it specifically, all
 * preserving the output exactly (verified by identical isolated-island
 * counts before/after at matched resolutions):
 *  - containerDistField caches the result of (1) so the boundary-clamp pass
 *    later in this function doesn't repeat the same query a second time;
 *  - when the thickness function's distance source is the container itself,
 *    its query reuses the same cached value instead of repeating (1) again;
 *  - classify_container_narrow_band() proves, for the large majority of
 *    voxels, that they are clearly inside or outside the container using a
 *    cheap coarse pre-pass, skipping (1) entirely for those voxels and
 *    falling back to the real query only in a thin band around the surface.
 */
bool GeneratorLewiner::compute_cached_field_values(const IContainer& con) {

	// ensure cached voxels is free
	cachedVoxels.clear();

	if (seeds.size() < 3) {
		if (logger) logger->log(LogPriority::ERROR, "We need at least three seeds!");
		return false;
	}

	// The iso-level is the wall/strut half-thickness anchor; a non-positive
	// value produces a field that never crosses the surface, i.e. an empty
	// mesh. Reject it here rather than silently generating nothing.
	if (isoLevel <= 0.0f) {
		if (logger) logger->log(LogPriority::ERROR,
			"Thickness (iso-level) must be positive; got " + std::to_string(isoLevel) + ".");
		return false;
	}

    // Control the number of anisotropy sources
    const bool variedAnisotropy = anisotropySources.size() >= 1;

    // Use the general anisotropy tensor as a background
    AnisotropySource background;
    background.direction = anisotropyVec;
	background.angle = anisotropyAngle;
    background.stretch = Vec3(stretchX, stretchY, stretchZ);
    create_metric(background); 

    // Ensure every source's Metric is up-to-date 
    for (auto& src : anisotropySources) {
        src->update_metric();
    }

    domainVolume = con.get_volume();
    update_steps();

    Eigen::Matrix3f rot = rotation_from_direction_roll(anisotropyVec, anisotropyAngle);
    Vec3 center = con.compute_bounds().center;

    std::vector<Vec3> warpedSeeds = this->seeds;

    if (!variedAnisotropy) {
        for (auto& seed : warpedSeeds) {
            Vec3 local = seed - center;
            Vec3 rotated = Vec3(rot * Eigen::Vector3f{ local.x, local.y, local.z });

            seed.x = rotated.x / stretchX;
            seed.y = rotated.y / stretchY;
            seed.z = rotated.z / stretchZ;
        }
    }

    // Populate the kdtree with the seeds (warped if uniform, unwarped if varied)
    std::unique_ptr<Kdtree> kdtree = std::make_unique<Kdtree>(warpedSeeds);

    // Decide the number 'k' of nearest neighbors for the kdtree
    std::vector<float> stretches = { stretchX, stretchY, stretchZ };
    for (const auto& src : anisotropySources) {
        stretches.push_back(src->stretch.x);
        stretches.push_back(src->stretch.y);
        stretches.push_back(src->stretch.z);
    }

    auto minmax = std::minmax_element(stretches.begin(), stretches.end());
    float sMin = *minmax.first;
    float sMax = *minmax.second;

    size_t Kn = variedAnisotropy ? choose_candidate_number(sMin, sMax, seeds.size()) : 3;

    // Update scalar field capacity
    scalarField.clear();
    size_t totalVoxels = static_cast<size_t>(blockDims[0]) *
                         static_cast<size_t>(blockDims[1]) *
                         static_cast<size_t>(blockDims[2]);

    std::cout << "Voxel Nr: " << totalVoxels << std::endl;
    scalarField.assign(totalVoxels, 9999.9f);

    // member (not local): assemble_field reuses it for the clamp + porosity mask
    containerDistField.assign(totalVoxels, 0.0f);

    // Numerical bands are expressed in grid spacings (voxels), not absolute
    // millimetres, so they track the model scale and resolution automatically.
    // surfaceMargin is the thin shell of valid lattice field kept just outside
    // the container wall so the boundary clamp, Taubin smoothing and marching
    // cubes all see a smooth gradient across the surface. Three voxels keeps it
    // strictly wider than the air-skin band (air_skip_level, two voxels).
    const float hMax = std::max(stepX, std::max(stepY, stepZ));
    const float surfaceMargin = 3.0f * hMax;

    const bool thicknessReusesContainerSDF = thicknessFunction && thicknessSDF && (thicknessSDF.get() == con.sdf.get());
    const float classificationMargin = thicknessReusesContainerSDF 
                                       ? std::max(surfaceMargin, transitionDistance) 
                                       : surfaceMargin;

    // Narrow-Band Classification Phase
    const auto fillStart = std::chrono::steady_clock::now();
    std::vector<NarrowBandClass> narrowBand = classify_container_narrow_band(con, classificationMargin);
    const auto classifyEnd = std::chrono::steady_clock::now();
    std::cout << "  [compute_scalar_field] narrow-band classification: "
              << std::chrono::duration<double>(classifyEnd - fillStart).count() << " s" << std::endl;

	// RAM guard for the per-voxel cache. compute_scalar_field always caches now
	// (calibration re-assembles off it), so every caller - the GUI, the CLI
	// profiler, the .scaf loader - reaches this allocation. The GUI has its own
	// pre-flight check, but the non-GUI callers do not, so guard here too: a
	// too-large grid should fail cleanly with an actionable message instead of
	// throwing std::bad_alloc (or OOM-killing) deep inside resize(). Convention
	// mirrors the MC sizing guard: a 0 from the query means "unknown, proceed"
	// and the bad_alloc backstop still applies.
	{
		const uint64_t bytesPerVoxel =
			sizeof(VoxelValue) + sizeof(float) /*scalarField*/ + sizeof(float) /*containerDist*/;
		const uint64_t requiredBytes = static_cast<uint64_t>(totalVoxels) * bytesPerVoxel;
		const uint64_t freeBytes = available_physical_memory_bytes();
		const double   budgetFrac = 0.80; // leave headroom for the OS / other buffers
		const bool memoryOk =
			(freeBytes == 0) ||
			(requiredBytes <= static_cast<uint64_t>(budgetFrac * static_cast<double>(freeBytes)));
		if (!memoryOk) {
			if (logger) logger->log(LogPriority::ERROR,
				"Scalar-field cache aborted: grid is " + std::to_string(totalVoxels) +
				" voxels (~" + std::to_string(requiredBytes >> 30) + " GB needed, " +
				std::to_string(freeBytes >> 30) + " GB free). Increase the generation "
				"voxel size (coarser grid) or restrict the domain.");
			else
				std::cerr << "Scalar-field cache aborted: grid too large ("
				          << totalVoxels << " voxels, ~" << (requiredBytes >> 30)
				          << " GB needed, " << (freeBytes >> 30) << " GB free)." << std::endl;
			cachedVoxels.clear();
			return false;
		}
	}

	// resize the vector store the cached voxels
	cachedVoxels.resize(totalVoxels);

	// Lambda for face hash
	auto stable_face_hash = [](size_t a, size_t b){
		size_t u = std::min(a, b);
		size_t v = std::max(a, b);
		
		// Simple pairing and bit-mixing (e.g., Murmur3 style integer hash)
		uint64_t hash = u;
		hash ^= v + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
		hash ^= hash >> 33;
		hash *= 0xff51afd7ed558ccdULL;
		hash ^= hash >> 33;
		hash *= 0xc4ceb9fe1a85ec53ULL;
		hash ^= hash >> 33;

		// this returns a value [0, 1.0], the probability of the face (a, b)
		return static_cast<float>(hash >> 40) * (1.0f / 16777216.0f);  
	};

    // Primary Voxel Evaluation Loop
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < blockDims[0]; i++) {
        for (int j = 0; j < blockDims[1]; j++) {
            for (int k = 0; k < blockDims[2]; k++) {

                size_t idx = find_vertex_index(i, j, k);
                float x = bounds[0] + i * stepX;
                float y = bounds[2] + j * stepY;
                float z = bounds[4] + k * stepZ;
                Vec3 point(x, y, z);

                // --- 1. Evaluate Container Boundary ---
                float containerDist;
                switch (narrowBand[idx]) {
                    case NarrowBandClass::Outside:
                        // any value safely past the classification band; sign is
                        // what matters (positive => outside, excluded from the
                        // domain and killed by the early-out below)
                        containerDist = classificationMargin + hMax;
                        break;
                    case NarrowBandClass::Inside:
                        // safely negative and beyond the thickness-grading
                        // transition, so graded thickness saturates correctly
                        containerDist = -(classificationMargin + hMax);
                        break;
                    default:
                        containerDist = con.sdf->compute_distance(point);
                        break;
                }

                containerDistField[idx] = containerDist;
                
                if (containerDist > surfaceMargin) {
                    scalarField[idx] = air_skip_level();
                    continue;
                }

                // --- 2. Thickness SDF distance (invariant to the thickness knob) ---
                // Only the DISTANCE is cached; the per-voxel iso level is re-derived
                // from it in assemble_field, so a thickness-scale change during
                // calibration does not need a fresh (expensive) cache pass.
                float rawDist = 0.0f;
                if (thicknessFunction && thicknessSDF) {
                    rawDist = std::abs(static_cast<float>((thicknessSDF.get() == con.sdf.get())
                        ? containerDist
                        : thicknessSDF->compute_distance(point)));
                }

                float d1 = 0.0f, d2 = 0.0f, d3 = 0.0f, d4 = 0.0f;
                Vec3 grad1, grad2, grad3, grad4;
				size_t id1 = 0, id2 = 0;

                // --- 3. Gather Anisotropic Nearest Neighbors ---
                if (variedAnisotropy) {
                    // Estimate the blended metric
                    Eigen::Matrix3f M = blend_metric(point, anisotropySources, background, backgroundWeight);

                    auto cand = kdtree->knn(point, Kn, [](const Vec3& a, const Vec3& b) {
                        Vec3 v = b - a; 
                        return double(v.x * v.x + v.y * v.y + v.z * v.z);
                    });

                    // Re-rank using local anisotropic metric
                    for (auto& c : cand) {
                        c.second = aniso_distance_sq(M, point, seeds[c.first]);
                    }
                    
                    std::partial_sort(
						cand.begin(), cand.begin() + std::min<size_t>(4, cand.size()), cand.end(),
						[](const auto& a, const auto& b) { return a.second < b.second; });

					id1 = cand[0].first;
					id2 = cand[1].first;

                    d1 = std::sqrt((float)cand[0].second);
                    d2 = std::sqrt((float)cand[1].second);
                    d3 = std::sqrt((float)cand[2].second);
                    
                    Vec3 p1 = seeds[cand[0].first];  
                    Vec3 p2 = seeds[cand[1].first];
                    Vec3 p3 = seeds[cand[2].first];
                    
                    grad1 = aniso_distance_grad(M, point, p1, d1);
                    grad2 = aniso_distance_grad(M, point, p2, d2);
                    grad3 = aniso_distance_grad(M, point, p3, d3);

                    // Fourth neighbour (for the optional edge-rounding smin).
                    // Falls back to the third if the seed set is too small, so
                    // the edge-vs-strut smin degenerates to a no-op there.
                    if (cand.size() >= 4) {
                        d4 = std::sqrt((float)cand[3].second);
                        grad4 = aniso_distance_grad(M, point, seeds[cand[3].first], d4);
                    } else {
                        d4 = d3;
                        grad4 = grad3;
                    }
                }
                else {
                    Vec3 localPt = point - center;
                    Vec3 rotatedPt = Vec3(rot * Eigen::Vector3f{ localPt.x, localPt.y, localPt.z });
                    Vec3 wrapped(rotatedPt.x / stretchX, rotatedPt.y / stretchY, rotatedPt.z / stretchZ);

                    auto cand = kdtree->knn(wrapped, 4, [this](const Vec3& p1, const Vec3& p2) {
                        Vec3 v = p2 - p1;
                        return (v.x * v.x) + (v.y * v.y) + (v.z * v.z);
                    });

					id1 = cand[0].first;
					id2 = cand[1].first;

                    d1 = (float)std::sqrt(cand[0].second);
                    d2 = (float)std::sqrt(cand[1].second);
                    d3 = (float)std::sqrt(cand[2].second);

                    Vec3 p1 = warpedSeeds[cand[0].first];
                    Vec3 p2 = warpedSeeds[cand[1].first];
                    Vec3 p3 = warpedSeeds[cand[2].first];

                    auto calc_grad = [&](const Vec3& p, float d) -> Vec3 {
                        if (d < 1e-6f) return Vec3(0.0f, 0.0f, 0.0f); 
                        Vec3 local(
                            (wrapped.x - p.x) / (d * stretchX),
                            (wrapped.y - p.y) / (d * stretchY),
                            (wrapped.z - p.z) / (d * stretchZ)
                        );
                        return Vec3(rot.transpose() * Eigen::Vector3f(local.x, local.y, local.z));
                    };

                    grad1 = calc_grad(p1, d1);
                    grad2 = calc_grad(p2, d2);
                    grad3 = calc_grad(p3, d3);

                    // Fourth neighbour (for the optional edge-rounding smin).
                    // Falls back to the third if the seed set is too small.
                    if (cand.size() >= 4) {
                        d4 = (float)std::sqrt(cand[3].second);
                        grad4 = calc_grad(warpedSeeds[cand[3].first], d4);
                    } else {
                        d4 = d3;
                        grad4 = grad3;
                    }
                }
			
			// 2. The two geometric extremes for this voxel.
			// wallVal vanishes on the Voronoi FACES (bisector of the nearest two
			// seeds) -> thickening it gives plates.
			// strutVal vanishes on the Voronoi EDGES (nearest three equidistant)
			// -> thickening it gives rods. Note wallVal <= strutVal always.
			float wallVal = d2 - d1;
			Vec3 wallGrad = grad2 - grad1;

			float strutVal = d3 - d1;
			Vec3 strutGrad = grad3 - grad1;

			// edgeVal vanishes on the Voronoi VERTICES (nearest four equidistant).
			// Only the optional edge-rounding smin reads it; strutVal <= edgeVal.
			float edgeVal = d4 - d1;
			Vec3 edgeGrad = grad4 - grad1;

			// 3. Per-face fenestration.
			// A Voronoi face is identified by its nearest seed PAIR, so hashing
			// that pair gives every face a fixed pseudo-random value in [0,1)
			// that all voxels near the face agree on. That makes each plate a
			// coherent membrane rather than speckle, needs no storage, and is
			// thread-safe and reproducible.
			size_t a = std::min(id1, id2);
			size_t b = std::max(id1, id2);

			float faceP = stable_face_hash(a, b);

			// cached voxel (all invariant to the calibration knobs iso/threshold/spread)
			VoxelValue cachedVoxel;
			cachedVoxel.wallVal = wallVal;
			cachedVoxel.wallGrad = wallGrad;
			cachedVoxel.strutVal = strutVal;
			cachedVoxel.strutGrad = strutGrad;
			cachedVoxel.edgeVal = edgeVal;
			cachedVoxel.edgeGrad = edgeGrad;
			cachedVoxel.faceP = faceP;
			cachedVoxel.rawDist = rawDist;

			cachedVoxels[idx] = cachedVoxel;
			// Each face gets its OWN openness, centred on the global threshold:
			//
			//   threshold ("openness") = MEAN openness of the structure.
			//        0 -> foam (solid plates on every face)
			//        1 -> lattice (bare rods on the Voronoi edges)
			//        It drives how much material there is, i.e. BV/TV.
			//
			//   spread = per-face DIVERSITY of that openness (rod/plate mixture).
			//        0 -> every face identical: the plain uniform morph, where
			//             all plates are fenestrated to the same degree.
			//        1 -> faces range from near-full plate to bare rod, i.e.
			//             genuine rod-and-plate coexistence like real trabecular bone, where some trabeculae are plate-like and others
			//             rod-like.
			//
			// An intermediate localTau eats a plate back toward its bounding
			// edges: that is a FENESTRATED plate. Different faces therefore end
			// up perforated to different degrees.
			//
			// E[localTau] = threshold exactly (faceP is uniform on [0,1)), so
			// 'spread' reshapes the rod/plate mixture without shifting the mean
			// openness. It is, however, NOT neutral for BV/TV: material is
			// nonlinear in tau, so spreading tau about a fixed mean nets more
			// solid (Jensen) - expect porosity to fall as spread rises, and
			// re-tune 'threshold' after changing it.
			//
			// Measured behaviour (4mm box, thickness 0.13, threshold 0.5):
			//   spread 0.0 -> 0.8 leaves Tb.Th within 0.2% (0.1357 -> 0.1360),
			//   while porosity drops 0.878 -> 0.855 and Conn.D 3.95 -> 3.02.
			//   So spread is orthogonal to thickness but not to BV/TV.
			// At MATCHED porosity (openness re-calibrated to compensate the leak),
			// spread is a genuine rod/plate control: SMI +49% (0.50 -> 0.75) over
			// spread 0 -> 1 with Tb.Th flat. openness + spread jointly span the
			// (porosity, SMI) plane - openness alone traces the perforation curve
			// where SMI and BV/TV are locked; spread moves off it. See §8.4.
            }
        }
    }

    const auto fillEnd = std::chrono::steady_clock::now();
    std::cout << "  [compute_cached_field_values] fill loop (container SDF + seed kdtree): "
              << std::chrono::duration<double>(fillEnd - fillStart).count() << " s" << std::endl;

    // Post-processing (smoothing, clamp, island removal, seal, porosity) is NOT
    // done here. It depends on the ASSEMBLED field and on the calibration knobs
    // (isoLevel / threshold / spread), so it lives in assemble_field(), which runs
    // once per iteration off this cache.
    return true;
}

// Cheap per-iteration pass: reconstruct the scalar field from the cached,
// knob-invariant per-voxel data (no kdtree, no container SDF), then run the
// knob-dependent post-processing. Called repeatedly by the calibrators after a
// single compute_cached_field_values(). Mirrors the tail of the old
// compute_scalar_field exactly, so the output is identical.
bool GeneratorLewiner::assemble_field(){

	if (cachedVoxels.empty()) {
		if (logger) logger->log(LogPriority::ERROR,
			"assemble_field: no cached field values; call compute_cached_field_values first.");
		return false;
	}

	const auto fillStart = std::chrono::steady_clock::now();

	// exterior band, identical to the cache pass's early-out
	const float hMax = std::max(stepX, std::max(stepY, stepZ));
	const float surfaceMargin = 3.0f * hMax;
	const float airLevel = air_skip_level();
	const bool variedThickness = thicknessFunction && thicknessSDF;

	const long long totalVoxels = static_cast<long long>(cachedVoxels.size());

	#pragma omp parallel for
	for (long long idx = 0; idx < totalVoxels; idx++) {

		// Exterior voxels are air: keep the stamp instead of writing a bogus
		// (all-zero-cache) value that would flip them to solid.
		if (containerDistField[idx] > surfaceMargin) {
			scalarField[idx] = airLevel;
			continue;
		}

		const VoxelValue& cachedVoxel = cachedVoxels[idx];

		// Re-derive the per-voxel iso level from the cached (invariant) distance
		// and the CURRENT thickness. Uniform: localIsoLevel = isoLevel, so the
		// shift below cancels to localRaw. Varied: tracks the calibrated scale.
		float localIsoLevel = isoLevel;
		if (variedThickness) {
			localIsoLevel = static_cast<float>(
				thicknessFunction->estimate_radius(cachedVoxel.rawDist, startThickness, endThickness));
		}

		// local tau (openness + per-face spread)
		float localTau = std::clamp(
			threshold + spread * (cachedVoxel.faceP - 0.5f), 0.0f, 1.0f);

		// Boundary frame: within frameDepth of the container wall, ramp the openness
		// toward zero so the outermost cells close into full Voronoi walls. This
		// leaves a connected honeycomb rim (the "frame") that ties the cut boundary
		// struts together. containerDistField is the container SDF (<=0 inside), so
		// depth = -SDF grows inward from the wall.
		// The band depth is expressed as a MULTIPLE of the wall thickness (isoLevel),
		// so it stays proportionate at any scale. This matters because tau -> 0 is the
		// plate field: a large ABSOLUTE depth turns a thick boundary shell plate-like
		// ("plates appear"), whereas tying it to the thickness keeps the frame a thin
		// rim. Default frameDepth = 1.0 => one wall layer (a band of isoLevel), which
		// is the shallowest depth that still reads as a solid honeycomb rim.
		if (frameBoundary && frameDepth > 1e-6f) {
			const float frameBand = frameDepth * isoLevel;
			const float depth = -containerDistField[idx];   // 0 at wall, >0 inward
			if (depth < frameBand) {
				float w = 1.0f - depth / frameBand;          // 1 at wall -> 0 at frameBand
				w = std::clamp(w, 0.0f, 1.0f);               // guard the thin outside band
				w = w * w * (3.0f - 2.0f * w);               // smoothstep for a gentle blend
				localTau *= (1.0f - w);                      // drive tau -> 0 at the wall
			}
		}

		// Base geometry per face. Two paths that agree exactly when edgeK -> 0:
		//
		//   value = (1-tau)(sd2 - sd1) + tau(sd3 - sd1)
		//
		// Linear blend (edge rounding off): sd_k = d_k, i.e. sd2-sd1 = wallVal,
		// sd3-sd1 = strutVal -- the validated crease-free field.
		//
		// Edge rounding on: sd2 = smin(d2,d3;k), sd3 = smin(d3,d4;k). By the
		// translation-equivariance of smin, sd2-sd1 = smin(wallVal, strutVal; k)
		// and sd3-sd1 = smin(strutVal, edgeVal; k) -- computed purely from the
		// cached contrasts. This rounds the Voronoi EDGES (d2=d3) and VERTICES
		// (d3=d4), independent of tau (so it also rounds tau=1 lattices). Since
		// wallVal <= strutVal <= edgeVal, at k=0 each smin picks the min branch
		// and the field reduces EXACTLY to the linear blend above.
		float value;
		Vec3 gradValue;
		if (roundEdges && edgeK > 1e-6f) {
			SmoothDist softWall = smin_gradient(
				cachedVoxel.wallVal, cachedVoxel.strutVal,
				cachedVoxel.wallGrad, cachedVoxel.strutGrad, edgeK);   // sd2 - sd1
			SmoothDist softStrut = smin_gradient(
				cachedVoxel.strutVal, cachedVoxel.edgeVal,
				cachedVoxel.strutGrad, cachedVoxel.edgeGrad, edgeK);   // sd3 - sd1
			value     = (1.0f - localTau) * softWall.val  + localTau * softStrut.val;
			gradValue = (1.0f - localTau) * softWall.grad + localTau * softStrut.grad;
		} else {
			// Exact geometric interpolation for this face (linear blend, crease-free).
			value     = (1.0f - localTau) * cachedVoxel.wallVal  + localTau * cachedVoxel.strutVal;
			gradValue = (1.0f - localTau) * cachedVoxel.wallGrad + localTau * cachedVoxel.strutGrad;
		}

		// Optional junction smoothing (smin fillet of the rod/plate fields).
		// Independent of edge rounding above; the two can be stacked.
		if (smoothJunctions && smoothK > 1e-6f) {
			float rodField = cachedVoxel.strutVal + (1.0f - localTau) * (cachedVoxel.strutVal - cachedVoxel.wallVal);
			SmoothDist sd = smin_gradient(value, rodField, gradValue, cachedVoxel.strutGrad, smoothK);
			value = sd.val;
			gradValue = sd.grad;
		}

		// Safe normalization: protect against gradient cancellation at nodes.
		float gradMag = gradValue.norm();
		float localRaw = value;
		if (gradMag > 1e-4f) {
			localRaw = (value / gradMag) * 2.0f;
		}

		// Shift the scalar field
		scalarField[idx] = localRaw - localIsoLevel + isoLevel;
	}

	const auto fillEnd = std::chrono::steady_clock::now();
	std::cout << "  [assemble_field] blend + shift: "
	          << std::chrono::duration<double>(fillEnd - fillStart).count() << " s" << std::endl;

	// --- Post-processing (knob-dependent; runs every assemble) ---
	smooth_scalar_field_taubin(iter, lambda, mu);

	// Boundary clamp: contiguous flat walk (see the note formerly on this loop).
	#pragma omp parallel for
	for (long long idx = 0; idx < totalVoxels; idx++) {
		scalarField[idx] = std::max(scalarField[idx], containerDistField[idx] + isoLevel);
	}

	// Container edge frame: union solid beams along the container's own edges
	// (a box's 12 edges, a cylinder's two circular rims) so the specimen gains a
	// rigid outer cage for gripping/loading. edgeDist is a true Euclidean distance
	// to the edge set, so (edgeDist - frameBeam + isoLevel) is an SDF-shifted field
	// whose iso-surface is a beam of radius frameBeam. Gated to the interior
	// (containerDist <= 0) so the beams stay flush with the outer wall. Off by
	// default; only box/cylinder are supported (a mesh has no simple edge set).
	if (frameContainerEdges && frameBeam > 1e-6f) {
		auto lockedCon = container.lock();
		if (lockedCon) {
			const ObjectType ctype = lockedCon->get_type();
			const long long nx = blockDims[0];
			const long long ny = blockDims[1];
			if (ctype == ObjectType::BoxContainerType) {
				const float x0 = bounds[0], x1 = bounds[1];
				const float y0 = bounds[2], y1 = bounds[3];
				const float z0 = bounds[4], z1 = bounds[5];
				#pragma omp parallel for
				for (long long idx = 0; idx < totalVoxels; idx++) {
					if (containerDistField[idx] > 0.0f) continue;   // interior only
					const long long i = idx % nx;
					const long long j = (idx / nx) % ny;
					const long long k = idx / (nx * ny);
					const float x = bounds[0] + i * stepX;
					const float y = bounds[2] + j * stepY;
					const float z = bounds[4] + k * stepZ;
					// distance to the nearest box edge = hypot of the two smallest
					// face distances (an edge is where those two faces meet).
					float a = std::min(x - x0, x1 - x);
					float b = std::min(y - y0, y1 - y);
					float c = std::min(z - z0, z1 - z);
					if (a > b) std::swap(a, b);
					if (b > c) std::swap(b, c);
					if (a > b) std::swap(a, b);              // a <= b <= c
					const float edgeDist = std::sqrt(a * a + b * b);
					scalarField[idx] = std::min(scalarField[idx], edgeDist - frameBeam + isoLevel);
				}
			}
			else if (ctype == ObjectType::CylinderContainerType) {
				const float cx = 0.5f * (bounds[0] + bounds[1]);
				const float cy = 0.5f * (bounds[2] + bounds[3]);
				const float zB = bounds[4], zT = bounds[5];
				const float R  = 0.5f * (bounds[1] - bounds[0]);   // radius from the AABB
				#pragma omp parallel for
				for (long long idx = 0; idx < totalVoxels; idx++) {
					if (containerDistField[idx] > 0.0f) continue;   // interior only
					const long long i = idx % nx;
					const long long j = (idx / nx) % ny;
					const long long k = idx / (nx * ny);
					const float x = bounds[0] + i * stepX;
					const float y = bounds[2] + j * stepY;
					const float z = bounds[4] + k * stepZ;
					const float dr = std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy)) - R;
					const float dB = std::sqrt(dr * dr + (z - zB) * (z - zB));
					const float dT = std::sqrt(dr * dr + (z - zT) * (z - zT));
					const float edgeDist = std::min(dB, dT);
					scalarField[idx] = std::min(scalarField[idx], edgeDist - frameBeam + isoLevel);
				}
			}
			else if (logger) {
				logger->log(LogPriority::WARNING,
					"Container edge frame is only supported for box and cylinder containers; skipped.");
			}
		}
	}

	const auto clampEnd = std::chrono::steady_clock::now();
	std::cout << "  [assemble_field] boundary clamp: "
	          << std::chrono::duration<double>(clampEnd - fillEnd).count() << " s" << std::endl;

	if (!isROI) {
		remove_isolated_islands();
	}

	seal_grid_boundaries();

	// Voxel-based porosity: numerator and denominator counted on the same grid so
	// the discretization bias cancels (see estimate_metrics for the mesh estimate).
	long long solidVoxels = 0;
	long long domainVoxels = 0;
	#pragma omp parallel for reduction(+:solidVoxels, domainVoxels)
	for (long long v = 0; v < totalVoxels; v++) {
		if (containerDistField[v] <= 0.0f) {
			domainVoxels++;
			if (scalarField[v] < isoLevel) solidVoxels++;
		}
	}
	if (domainVoxels > 0) {
		porosity = 1.0f - static_cast<float>(solidVoxels) / static_cast<float>(domainVoxels);
	}

	const auto cleanupEnd = std::chrono::steady_clock::now();
	std::cout << "  [assemble_field] island removal + seal: "
	          << std::chrono::duration<double>(cleanupEnd - clampEnd).count() << " s" << std::endl;

	return true;
}

// Full regeneration: cache pass + one assemble. Every non-calibration caller uses
// this; calibration caches once and re-assembles per iteration.
bool GeneratorLewiner::compute_scalar_field(const IContainer& con) {
	return compute_cached_field_values(con) && assemble_field();
}

// void GeneratorLewiner::smooth_scalar_field() {
// 	std::vector<float> smoothed = scalarField; // Copy

// 	// Simple 3x3x3 Box Blur
// 	#pragma omp parallel for collapse(3)
// 	for (int z = 1; z < blockDims[2] - 1; z++) {
// 		for (int y = 1; y < blockDims[1] - 1; y++) {
// 			for (int x = 1; x < blockDims[0] - 1; x++) {

// 				float sum = 0.0;
// 				int count = 0;

// 				size_t centerIdx = find_vertex_index(x, y, z);
// 				if (scalarField[centerIdx] > isoLevel + 0.5f) {
// 					// It's safely outside the container and the boundary margin, skip smoothing
// 					continue;
// 				}

// 				// Average neighbors
// 				for (int kz = -1; kz <= 1; kz++) {
// 					for (int ky = -1; ky <= 1; ky++) {
// 						for (int kx = -1; kx <= 1; kx++) {
// 							size_t idx = find_vertex_index(x + kx, y + ky, z + kz);
// 							sum += scalarField[idx];
// 							count++;
// 						}
// 					}
// 				}

// 				smoothed[find_vertex_index(x, y, z)] = sum / count;

// 			}
// 		}
// 	}

// 	scalarField = smoothed; // Swap back
// }

void GeneratorLewiner::smooth_scalar_field() {
    std::vector<float> smoothed = scalarField; // Copy

    // 3D Binomial/Gaussian Filter (1-2-1 kernel)
    // The sum of all weights in this 3x3x3 kernel is exactly 64.
    const float invWeight = 1.0f / 64.0f;

    // scale-relative air-skip band, identical to the value stamped into
    // exterior voxels by compute_scalar_field
    const float skipLevel = air_skip_level();

    #pragma omp parallel for collapse(3)
    for (int z = 1; z < blockDims[2] - 1; z++) {
        for (int y = 1; y < blockDims[1] - 1; y++) {
            for (int x = 1; x < blockDims[0] - 1; x++) {

                size_t centerIdx = find_vertex_index(x, y, z);

                // If the voxel is safely outside the surface margin, skip it to save massive CPU time
                // (>= because exterior voxels are set to exactly air_skip_level())
                if (scalarField[centerIdx] >= skipLevel) {
                    continue;
                }

                float sum = 0.0f;

                // Convolve with the 1-2-1 weighted kernel
                for (int kz = -1; kz <= 1; kz++) {
                    int wz = (kz == 0) ? 2 : 1;
                    
                    for (int ky = -1; ky <= 1; ky++) {
                        int wy = (ky == 0) ? 2 : 1;
                        
                        for (int kx = -1; kx <= 1; kx++) {
                            int wx = (kx == 0) ? 2 : 1;
                            
                            int weight = wx * wy * wz;
                            size_t idx = find_vertex_index(x + kx, y + ky, z + kz);
                            
                            sum += scalarField[idx] * weight;
                        }
                    }
                }

                smoothed[centerIdx] = sum * invWeight;
            }
        }
    }

    scalarField = smoothed; // Swap back
}

void GeneratorLewiner::remove_isolated_islands() {
	size_t nx = static_cast<size_t>(blockDims[0]);
	size_t ny = static_cast<size_t>(blockDims[1]);
	size_t nz = static_cast<size_t>(blockDims[2]);
	size_t totalVoxels = nx * ny * nz;

	std::vector<bool> visited(totalVoxels, false);
	std::vector<std::vector<size_t>> islands;

	// 1. Find all connected components
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {

				size_t idx = find_vertex_index(i, j, k);

				// If voxel is solid and hasn't been grouped into an island yet
				if (scalarField[idx] < isoLevel && !visited[idx]) {

					std::vector<size_t> currentIsland;
					std::queue<Vec3i> q;

					q.push({ i, j, k });
					visited[idx] = true;

					// Standard Breadth-First Search (Flood Fill)
					while (!q.empty()) {
						Vec3i curr = q.front();
						q.pop();

						currentIsland.push_back(find_vertex_index(curr.x, curr.y, curr.z));

						// The 6 adjacent neighbors (up, down, left, right, front, back)
						Vec3i neighbors[6] = {
							{curr.x + 1, curr.y, curr.z}, {curr.x - 1, curr.y, curr.z},
							{curr.x, curr.y + 1, curr.z}, {curr.x, curr.y - 1, curr.z},
							{curr.x, curr.y, curr.z + 1}, {curr.x, curr.y, curr.z - 1}
						};

						for (const auto& n : neighbors) {
							// Ensure neighbor is within grid bounds
							if (n.x >= 0 && n.x < nx && n.y >= 0 && n.y < ny && n.z >= 0 && n.z < nz) {
								size_t nIdx = find_vertex_index(n.x, n.y, n.z);
								if (scalarField[nIdx] < isoLevel && !visited[nIdx]) {
									visited[nIdx] = true;
									q.push(n);
								}
							}
						}
					}
					islands.push_back(currentIsland);
				}
			}
		}
	}

	if (islands.empty()) return;

	// 2. Find the largest island (The main structure)
	size_t maxIslandIdx = 0;
	size_t maxSize = 0;
	for (size_t i = 0; i < islands.size(); i++) {
		if (islands[i].size() > maxSize) {
			maxSize = islands[i].size();
			maxIslandIdx = i;
		}
	}

	// 3. Erase all other islands by filling them with "Air"
	int removedCount = 0;
	for (size_t i = 0; i < islands.size(); i++) {
		if (i != maxIslandIdx) {
			for (size_t idx : islands[i]) {
				scalarField[idx] = 9999.9f;
				removedCount++;
			}
		}
	}

	std::cout << "Cleaned up: Removed " << islands.size() - 1
		<< " floating pieces (" << removedCount << " voxels)." << std::endl;
}

//-----------------------------------------------------------------------------

bool GeneratorLewiner::marching_cubes(bool supress) {

	size_t expectedVoxels = static_cast<size_t>(blockDims[0]) *
		static_cast<size_t>(blockDims[1]) *
		static_cast<size_t>(blockDims[2]);
	if (scalarField.size() != expectedVoxels) {
		if (logger) logger->log(LogPriority::ERROR,
			"Cannot run marching cubes: scalar field is empty or does not match the grid!");
		return false;
	}

	meshVertices.clear();
	meshTriangles.clear();
	x_verts.clear();
	y_verts.clear();
	z_verts.clear();

	// get current time
	const auto start{ std::chrono::steady_clock::now() };

	size_t totalVoxels = static_cast<size_t>(blockDims[0]) *
		static_cast<size_t>(blockDims[1]) *
		static_cast<size_t>(blockDims[2]);
	x_verts.assign(totalVoxels, SIZE_MAX);
	y_verts.assign(totalVoxels, SIZE_MAX);
	z_verts.assign(totalVoxels, SIZE_MAX);

	// Pre-allocating meshVertices/meshTriangles for the theoretical worst
	// case (every grid edge crossing the surface, every cell maximally
	// subdivided) requires buffers proportional to totalVoxels - at
	// billion-voxel grids that is hundreds of GB and reliably exhausts
	// available RAM, even though the real (sparse) surface only ever
	// touches a small fraction of the grid.
	//
	// Instead, two cheap counting passes mirror - exactly, not
	// heuristically - the crossing tests compute_intersection_points() and
	// process_cube() use below, just without allocating or writing
	// anything, to get an EXACT lower bound on how many vertices/triangles
	// the real passes will produce. The buffers are then sized from that
	// real count instead of the full grid size. compute_intersection_points()
	// and process_cube()/add_triangle() are otherwise unchanged and still
	// append via the same thread-safe atomic fetch_add into these buffers.

	// pass 1: count edge vertices (mirrors compute_intersection_points())
	std::atomic<size_t> edgeVertexCount{ 0 };

	#pragma omp parallel for collapse(3)
	for (int i = 0; i < blockDims[0]; i++) {
		for (int j = 0; j < blockDims[1]; j++) {
			for (int k = 0; k < blockDims[2]; k++) {

				float val0 = get_data(i, j, k) - isoLevel;
				float val1 = (i < blockDims[0] - 1) ? get_data(i + 1, j, k) - isoLevel : val0;
				float val2 = (j < blockDims[1] - 1) ? get_data(i, j + 1, k) - isoLevel : val0;
				float val3 = (k < blockDims[2] - 1) ? get_data(i, j, k + 1) - isoLevel : val0;

				if (fabs(val0) < FLT_EPSILON) val0 = (val0 < 0) ? -FLT_EPSILON : FLT_EPSILON;
				if (fabs(val1) < FLT_EPSILON) val1 = (val1 < 0) ? -FLT_EPSILON : FLT_EPSILON;
				if (fabs(val2) < FLT_EPSILON) val2 = (val2 < 0) ? -FLT_EPSILON : FLT_EPSILON;
				if (fabs(val3) < FLT_EPSILON) val3 = (val3 < 0) ? -FLT_EPSILON : FLT_EPSILON;

				size_t crossings = 0;
				if ((val0 < 0 && val1 >= 0) || (val0 >= 0 && val1 < 0)) crossings++;
				if ((val0 < 0 && val2 >= 0) || (val0 >= 0 && val2 < 0)) crossings++;
				if ((val0 < 0 && val3 >= 0) || (val0 >= 0 && val3 < 0)) crossings++;

				if (crossings > 0) {
					edgeVertexCount.fetch_add(crossings, std::memory_order_relaxed);
				}
			}
		}
	}

	// pass 2: count surface-crossing cells (mirrors the lut_entry
	// classification in the process_cube loop below)
	std::atomic<size_t> crossingCellCount{ 0 };

	#pragma omp parallel for collapse(3)
	for (int i = 0; i < blockDims[0] - 1; i++) {
		for (int j = 0; j < blockDims[1] - 1; j++) {
			for (int k = 0; k < blockDims[2] - 1; k++) {

				int lut_entry = 0;
				for (int p = 0; p < 8; ++p) {
					float v = get_data(i + ((p ^ (p >> 1)) & 1), j + ((p >> 1) & 1), k + ((p >> 2) & 1)) - isoLevel;
					if (v > 0) lut_entry |= (1 << p);
				}
				if (lut_entry != 0 && lut_entry != 255) {
					crossingCellCount.fetch_add(1, std::memory_order_relaxed);
				}
			}
		}
	}

	const auto countEnd = std::chrono::steady_clock::now();
	std::cout << "  [marching_cubes] surface-crossing pre-count (" << edgeVertexCount.load() << " edges, "
		<< crossingCellCount.load() << " of " << totalVoxels << " cells): "
		<< std::chrono::duration<double>(countEnd - start).count() << " s" << std::endl;

	// exact bound: edgeVertexCount vertices are created by
	// compute_intersection_points() below, plus at most one additional
	// "c_vertex" per crossing cell created later by process_cube(); triangle
	// count is bounded by the LUT's documented true worst case of 12
	// triangles for a single cell (case 13.4 - see add_triangle() call
	// sites). "+64" only matters for degenerate near-empty grids.
	size_t vertexCapacity = edgeVertexCount.load() + crossingCellCount.load() + 64;
	size_t triangleCapacity = crossingCellCount.load() * 12 + 64;
	meshVertices.resize(vertexCapacity);
	meshTriangles.resize(triangleCapacity);
	vertexCount = 0;
	triangleCount = 0;

	compute_intersection_points();

	const auto intersectEnd = std::chrono::steady_clock::now();
	std::cout << "  [marching_cubes] compute_intersection_points: "
		<< std::chrono::duration<double>(intersectEnd - countEnd).count() << " s" << std::endl;

	for (int i = 0; i < blockDims[0] - 1; i++) {
		for (int j = 0; j < blockDims[1] - 1; j++) {
			for (int k = 0; k < blockDims[2] - 1; k++) {
				
				int lut_entry = 0;

				float _cube[8];

				for (int p = 0; p < 8; ++p) {

					// ^ XOR this line

					_cube[p] = get_data(i + ((p ^ (p >> 1)) & 1), j + ((p >> 1) & 1), k + ((p >> 2) & 1)) - isoLevel;
					if (fabs(_cube[p]) < FLT_EPSILON) {
						_cube[p] = (_cube[p] < 0) ? -FLT_EPSILON : FLT_EPSILON;
					}
					
					if (_cube[p] > 0) lut_entry |= (1 << p);
				}
				if (lut_entry == 0 || lut_entry == 255) {
					continue;
				}

				process_cube(i, j, k, _cube, lut_entry);
			}
		}
	}

	const auto cubesEnd = std::chrono::steady_clock::now();
	std::cout << "  [marching_cubes] process_cube loop: "
		<< std::chrono::duration<double>(cubesEnd - intersectEnd).count() << " s" << std::endl;

	// update to exact size
	meshVertices.resize(vertexCount);
	meshTriangles.resize(triangleCount);

	validate_topology(true);

	// build adjacency
	build_topology();

	// update axis aligned bounding box
	_update_bounding_box();

	// update mesh version
	meshVersion++;

	const auto topoEnd = std::chrono::steady_clock::now();
	std::cout << "  [marching_cubes] topology + bbox: "
		<< std::chrono::duration<double>(topoEnd - cubesEnd).count() << " s" << std::endl;
	std::cout << "  [marching_cubes] TOTAL: "
		<< std::chrono::duration<double>(topoEnd - start).count() << " s" << std::endl;

	return true;
};

void GeneratorLewiner::_update_bounding_box() {

	aabb.pMin.x = std::numeric_limits<float>::max();
	aabb.pMin.y = std::numeric_limits<float>::max();
	aabb.pMin.z = std::numeric_limits<float>::max();
	aabb.pMax.x = std::numeric_limits<float>::lowest();
	aabb.pMax.y = std::numeric_limits<float>::lowest();
	aabb.pMax.z = std::numeric_limits<float>::lowest();

	for (const auto& v : meshVertices) {
		aabb.pMin.x = std::min(aabb.pMin.x, v.x);
		aabb.pMax.x = std::max(aabb.pMax.x, v.x);

		aabb.pMin.y = std::min(aabb.pMin.y, v.y);
		aabb.pMax.y = std::max(aabb.pMax.y, v.y);

		aabb.pMin.z = std::min(aabb.pMin.z, v.z);
		aabb.pMax.z = std::max(aabb.pMax.z, v.z);
	}

	// NOTE: 'bounds' (the generation-grid definition) is deliberately NOT
	// touched here. The mesh AABB is always strictly inside the container
	// (sealed shell + sub-voxel marching-cubes inset), so copying it into
	// 'bounds' shrank the sampling domain on every regeneration.
	// External callers that want the tight mesh box use get_bounds()/get_aabb().
};

void GeneratorLewiner::compute_intersection_points() {
	
	#pragma omp parallel for collapse(3)
	for (int i = 0; i < blockDims[0]; i++) {
		for (int j = 0; j < blockDims[1]; j++) {
			for (int k = 0; k < blockDims[2]; k++) {

				// get current voxel index
				size_t currentIdx = find_vertex_index(i, j, k);

				float val0 = get_data(i, j, k) - isoLevel;

				// get physical position in space
				Vec3 pos = get_position(i, j, k);

				// x edge
				float val1 = 0.0f;
				if (i < blockDims[0] - 1) {

					val1 = get_data(i + 1, j, k) - isoLevel;
				}
				else {
					val1 = val0;
				}

				// y edge
				float val2 = 0.0f;
				if (j < blockDims[1] - 1) {

					val2 = get_data(i, j + 1, k) - isoLevel;
				}
				else {
					val2 = val0;
				}

				float val3 = 0.0f;
				if (k < blockDims[2] - 1) {

					val3 = get_data(i, j, k + 1) - isoLevel;
				}
				else {
					val3 = val0;
				}

				if (fabs(val0) < FLT_EPSILON) val0 = (val0 < 0) ? -FLT_EPSILON : FLT_EPSILON;
				if (fabs(val1) < FLT_EPSILON) val1 = (val1 < 0) ? -FLT_EPSILON : FLT_EPSILON;
				if (fabs(val2) < FLT_EPSILON) val2 = (val2 < 0) ? -FLT_EPSILON : FLT_EPSILON;
				if (fabs(val3) < FLT_EPSILON) val3 = (val3 < 0) ? -FLT_EPSILON : FLT_EPSILON;


				// check if the surface crosses the edge, they have different signs
				if ((val0 < 0 && val1 >= 0) || (val0 >= 0 && val1 < 0)) {
					size_t vIdx = vertexCount.fetch_add(1, std::memory_order_relaxed);
					
					// add to vector of Vertices
					meshVertices[vIdx] = add_x_vertex(i, j, k, val0, val1);
					
					// cash current to x vertices
					x_verts[currentIdx] = vIdx;
				}

				if ((val0 < 0 && val2 >= 0) || (val0>=0 && val2 < 0)) {
					size_t vIdx = vertexCount.fetch_add(1, std::memory_order_relaxed);
					meshVertices[vIdx] = add_y_vertex(i, j, k, val0, val2);
					y_verts[currentIdx] = vIdx;
				}

				if ((val0 < 0 && val3 >= 0) || (val0 >= 0 && val3 < 0)) {
					size_t vIdx = vertexCount.fetch_add(1, std::memory_order_relaxed);
					meshVertices[vIdx] = add_z_vertex(i, j, k, val0, val3);
					z_verts[currentIdx] = vIdx;
				};
			}
		}
	}
};

// get data from scalar field
float GeneratorLewiner::get_data(const int i, const int j, const int k) const {

	size_t idx = static_cast<size_t>(i) +
		static_cast<size_t>(j) * static_cast<size_t>(blockDims[0]) +
		static_cast<size_t>(k) * static_cast<size_t>(blockDims[0]) * static_cast<size_t>(blockDims[1]);

	return scalarField[idx];
};

//@brief function to find the index of the vertex in the scalar field vector based on its position in the grid
size_t GeneratorLewiner::find_vertex_index(int x, int y, int z) {
	return static_cast<size_t>(x) +
		static_cast<size_t>(y) * static_cast<size_t>(blockDims[0]) +
		static_cast<size_t>(z) * static_cast<size_t>(blockDims[0]) * static_cast<size_t>(blockDims[1]);
}

Vec3 GeneratorLewiner::get_position(int x, int y, int z) {

	return Vec3(
		bounds[0] + x * stepX,
		bounds[2] + y * stepY,
		bounds[4] + z * stepZ);
};

void GeneratorLewiner::process_cube(int i, int j, int k, const float cube[8], int lut_entry) {

	size_t v12 = SIZE_MAX;

	// get te case from the LUT table
	int baseCase = cases[lut_entry][0];

	// get the configuration
	int baseConfig = cases[lut_entry][1];

	int subConfig = 0;

	int tunnelOrientation = 0;

	switch (baseCase) {

	case 0:
		break;
	case 1:
		add_triangle(tiling1[baseConfig], i, j, k, 1);
		break;

	case 2:
		add_triangle(tiling2[baseConfig], i, j, k, 2);
		break;

	case 3:
		if (test_face(test3[baseConfig], cube))
			add_triangle(tiling3_2[baseConfig], i, j, k, 4); // 3.2
		else
			add_triangle(tiling3_1[baseConfig], i, j, k, 2); // 3.1
		break;

	case 4:
		if (modified_test_interior(test4[baseConfig], baseCase, baseConfig, cube))
			add_triangle(tiling4_1[baseConfig], i, j, k, 2); // 4.1.1
		else
			add_triangle(tiling4_2[baseConfig], i, j, k, 6); // 4.1.2
		break;

	case 5:
		add_triangle(tiling5[baseConfig], i, j, k, 3);
		break;

	case 6:
		if (test_face(test6[baseConfig][0], cube))
			add_triangle(tiling6_2[baseConfig], i, j, k, 5); // 6.2
		else {
			if (modified_test_interior(test6[baseConfig][1], baseCase, baseConfig, cube))
				add_triangle(tiling6_1_1[baseConfig], i, j, k, 3); // 6.1.1
			else {
				v12 = add_c_vertex(i, j, k);
				add_triangle(tiling6_1_2[baseConfig], i, j, k, 9, v12); // 6.1.2
			}
		}
		break;

	case 7:
		if (test_face(test7[baseConfig][0], cube))
			subConfig += 1;
		if (test_face(test7[baseConfig][1], cube))
			subConfig += 2;
		if (test_face(test7[baseConfig][2], cube))
			subConfig += 4;
		switch (subConfig) {
		case 0:
			add_triangle(tiling7_1[baseConfig], i, j, k, 3);
			break;
		case 1:
			add_triangle(tiling7_2[baseConfig][0], i, j, k, 5);
			break;
		case 2:
			add_triangle(tiling7_2[baseConfig][1], i, j, k, 5);
			break;
		case 3:
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling7_3[baseConfig][0], i, j, k, 9, v12);
			break;
		case 4:
			add_triangle(tiling7_2[baseConfig][2], i, j, k, 5);
			break;
		case 5:
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling7_3[baseConfig][1], i, j, k, 9, v12);
			break;
		case 6:
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling7_3[baseConfig][2], i, j, k, 9, v12);
			break;
		case 7:
			if (modified_test_interior(test7[baseConfig][3], baseCase, baseConfig, cube))
				add_triangle(tiling7_4_1[baseConfig], i, j, k, 5);
			else
				add_triangle(tiling7_4_2[baseConfig], i, j, k, 9);
			break;
		}
		break;

	case 8:
		add_triangle(tiling8[baseConfig], i, j, k, 2);
		break;

	case 9:
		add_triangle(tiling9[baseConfig], i, j, k, 4);
		break;

	case 10:
		if (test_face(test10[baseConfig][0], cube)) {
			if (test_face(test10[baseConfig][1], cube)) {
				if (modified_test_interior(-(signed char)test10[baseConfig][2], baseCase, baseConfig, cube))
					add_triangle(tiling10_1_1_[baseConfig], i, j, k, 4); // 10.1.1
				else
					add_triangle(tiling10_1_2[5 - baseConfig], i, j, k, 8); // 10.1.2

			}
			else {
				v12 = add_c_vertex(i, j, k);
				add_triangle(tiling10_2[baseConfig], i, j, k, 8, v12); // 10.2
			}
		}
		else {
			if (test_face(test10[baseConfig][1], cube)) {
				v12 = add_c_vertex(i, j, k);
				add_triangle(tiling10_2_[baseConfig], i, j, k, 8, v12); // 10.2
			}
			else {
				if (modified_test_interior(test10[baseConfig][2], baseCase, baseConfig, cube))
					add_triangle(tiling10_1_1[baseConfig], i, j, k, 4); // 10.1.1
				else
					add_triangle(tiling10_1_2[baseConfig], i, j, k, 8); // 10.1.2
			}
		}
		break;

	case 11:
		add_triangle(tiling11[baseConfig], i, j, k, 4);
		break;

	case 12:
		if (test_face(test12[baseConfig][0], cube)) {
			if (test_face(test12[baseConfig][1], cube)) {
				if (modified_test_interior(-(signed char)test12[baseConfig][2], baseCase, baseConfig, cube))
					add_triangle(tiling12_1_1_[baseConfig], i, j, k, 4); // 12.1.1
				else
					add_triangle(tiling12_1_2[23 - baseConfig], i, j, k, 8); // 12.1.2
			}
			else {
				v12 = add_c_vertex(i, j, k);
				add_triangle(tiling12_2[baseConfig], i, j, k, 8, v12); // 12.2
			}
		}
		else {
			if (test_face(test12[baseConfig][1], cube)) {
				v12 = add_c_vertex(i, j, k);
				add_triangle(tiling12_2_[baseConfig], i, j, k, 8, v12); // 12.2
			}
			else {
				if (modified_test_interior(test12[baseConfig][2], baseCase, baseConfig, cube))
					add_triangle(tiling12_1_1[baseConfig], i, j, k, 4); // 12.1.1
				else
					add_triangle(tiling12_1_2[baseConfig], i, j, k, 8); // 12.1.2
			}
		}
		break;

	case 13:
		if (test_face(test13[baseConfig][0], cube))
			subConfig += 1;
		if (test_face(test13[baseConfig][1], cube))
			subConfig += 2;
		if (test_face(test13[baseConfig][2], cube))
			subConfig += 4;
		if (test_face(test13[baseConfig][3], cube))
			subConfig += 8;
		if (test_face(test13[baseConfig][4], cube))
			subConfig += 16;
		if (test_face(test13[baseConfig][5], cube))
			subConfig += 32;
		//interior_test_case13_2(0.0f, cube, tunnelOrientation);
		switch (subconfig13[subConfig]) {
		case 0:/* 13.1 */
			add_triangle(tiling13_1[baseConfig], i, j, k, 4);
			break;

		case 1:/* 13.2 */
			add_triangle(tiling13_2[baseConfig][0], i, j, k, 6);
			break;
		case 2:/* 13.2 */
			add_triangle(tiling13_2[baseConfig][1], i, j, k, 6);
			break;
		case 3:/* 13.2 */
			add_triangle(tiling13_2[baseConfig][2], i, j, k, 6);
			break;
		case 4:/* 13.2 */
			add_triangle(tiling13_2[baseConfig][3], i, j, k, 6);
			break;
		case 5:/* 13.2 */
			add_triangle(tiling13_2[baseConfig][4], i, j, k, 6);
			break;
		case 6:/* 13.2 */
			add_triangle(tiling13_2[baseConfig][5], i, j, k, 6);
			break;

		case 7:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3[baseConfig][0], i, j, k, 10, v12);
			break;
		case 8:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3[baseConfig][1], i, j, k, 10, v12);
			break;
		case 9:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3[baseConfig][2], i, j, k, 10, v12);
			break;
		case 10:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3[baseConfig][3], i, j, k, 10, v12);
			break;
		case 11:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3[baseConfig][4], i, j, k, 10, v12);
			break;
		case 12:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3[baseConfig][5], i, j, k, 10, v12);
			break;
		case 13:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3[baseConfig][6], i, j, k, 10, v12);
			break;
		case 14:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3[baseConfig][7], i, j, k, 10, v12);
			break;
		case 15:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3[baseConfig][8], i, j, k, 10, v12);
			break;
		case 16:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3[baseConfig][9], i, j, k, 10, v12);
			break;
		case 17:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3[baseConfig][10], i, j, k, 10, v12);
			break;
		case 18:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3[baseConfig][11], i, j, k, 10, v12);
			break;

		case 19:/* 13.4 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_4[baseConfig][0], i, j, k, 12, v12);
			break;
		case 20:/* 13.4 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_4[baseConfig][1], i, j, k, 12, v12);
			break;
		case 21:/* 13.4 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_4[baseConfig][2], i, j, k, 12, v12);
			break;
		case 22:/* 13.4 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_4[baseConfig][3], i, j, k, 12, v12);
			break;

		case 23:/* 13.5 */
			subConfig = 0;
			if (baseConfig == 0) {
				if (interior_test_case13(cube))
					add_triangle(tiling13_5_1[0][0], i, j, k, 6);
				else {
					if (tunnelOrientation == 1)
						add_triangle(tiling13_5_2[0][0], i, j, k, 10);
					else
						add_triangle(tiling13_5_2[1][2], i, j, k, 10);
				}
			}
			else {
				if (interior_test_case13(cube))
					add_triangle(tiling13_5_1[1][0], i, j, k, 6);
				else {
					if (tunnelOrientation == 1)
						add_triangle(tiling13_5_2[1][0], i, j, k, 10);
					else
						add_triangle(tiling13_5_2[0][2], i, j, k, 10);
				}
			}
			break;

		case 24:/* 13.5 */
			subConfig = 1;
			if (baseConfig == 0) {
				if (interior_test_case13(cube))
					add_triangle(tiling13_5_1[0][1], i, j, k, 6);
				else {
					if (tunnelOrientation == 1)
						add_triangle(tiling13_5_2[0][1], i, j, k, 10);
					else
						add_triangle(tiling13_5_2[1][0], i, j, k, 10);
				}
			}
			else
			{
				if (interior_test_case13(cube))
					add_triangle(tiling13_5_1[1][1], i, j, k, 6);
				else {
					if (tunnelOrientation == 1)
						add_triangle(tiling13_5_2[1][1], i, j, k, 10);
					else
						add_triangle(tiling13_5_2[0][3], i, j, k, 10);
				}
			}

			break;

		case 25:/* 13.5 */
			subConfig = 2;
			if (baseConfig == 0)
			{
				if (interior_test_case13(cube))
					add_triangle(tiling13_5_1[0][2], i, j, k, 6);
				else {
					if (tunnelOrientation == 1)
						add_triangle(tiling13_5_2[0][2], i, j, k, 10);
					else
						add_triangle(tiling13_5_2[1][3], i, j, k, 10);
				}
			}
			else
			{
				if (interior_test_case13(cube))
					add_triangle(tiling13_5_1[1][2], i, j, k, 6);
				else {
					if (tunnelOrientation == 1)
						add_triangle(tiling13_5_2[1][2], i, j, k, 10);
					else
						add_triangle(tiling13_5_2[0][0], i, j, k, 10);
				}

			}
			break;

		case 26: /* 13.5 */
			subConfig = 3;
			if (baseConfig == 0)
			{
				if (interior_test_case13(cube))
					add_triangle(tiling13_5_1[0][3], i, j, k, 6);
				else {
					if (tunnelOrientation == 1)
						add_triangle(tiling13_5_2[0][3], i, j, k, 10);
					else
						add_triangle(tiling13_5_2[1][1], i, j, k, 10);
				}
			}
			else
			{
				if (interior_test_case13(cube))
					add_triangle(tiling13_5_1[1][3], i, j, k, 6);
				else {
					if (tunnelOrientation == 1)
						add_triangle(tiling13_5_2[1][3], i, j, k, 10);
					else
						add_triangle(tiling13_5_2[0][2], i, j, k, 10);
				}
			}
			/* 13.4  common node is negative*/
			// v12 = add_c_vertex() ;
			// add_triangle( tiling13_4[baseConfig][3], 12, v12 ) ;
			break;

		case 27:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3_[baseConfig][0], i, j, k, 10, v12);
			break;
		case 28:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3_[baseConfig][1], i, j, k, 10, v12);
			break;
		case 29:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3_[baseConfig][2], i, j, k, 10, v12);
			break;
		case 30:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3_[baseConfig][3], i, j, k, 10, v12);
			break;
		case 31:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3_[baseConfig][4], i, j, k, 10, v12);
			break;
		case 32:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3_[baseConfig][5], i, j, k, 10, v12);
			break;
		case 33:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3_[baseConfig][6], i, j, k, 10, v12);
			break;
		case 34:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3_[baseConfig][7], i, j, k, 10, v12);
			break;
		case 35:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3_[baseConfig][8], i, j, k, 10, v12);
			break;
		case 36:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3_[baseConfig][9], i, j, k, 10, v12);
			break;
		case 37:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3_[baseConfig][10], i, j, k, 10, v12);
			break;
		case 38:/* 13.3 */
			v12 = add_c_vertex(i, j, k);
			add_triangle(tiling13_3_[baseConfig][11], i, j, k, 10, v12);
			break;

		case 39:/* 13.2 */
			add_triangle(tiling13_2_[baseConfig][0], i, j, k, 6);
			break;
		case 40:/* 13.2 */
			add_triangle(tiling13_2_[baseConfig][1], i, j, k, 6);
			break;
		case 41:/* 13.2 */
			add_triangle(tiling13_2_[baseConfig][2], i, j, k, 6);
			break;
		case 42:/* 13.2 */
			add_triangle(tiling13_2_[baseConfig][3], i, j, k, 6);
			break;
		case 43:/* 13.2 */
			add_triangle(tiling13_2_[baseConfig][4], i, j, k, 6);
			break;
		case 44:/* 13.2 */
			add_triangle(tiling13_2_[baseConfig][5], i, j, k, 6);
			break;

		case 45:/* 13.1 */
			add_triangle(tiling13_1_[baseConfig], i, j, k, 4);
			break;

		default:
			printf("Marching Cubes: Impossible case 13?\n");
			//print_cube();
		}
		break;

	case 14:
		add_triangle(tiling14[baseConfig], i, j, k, 4);
		break;
	};
};

// const uint8_t* points to the first element of the array
void GeneratorLewiner::add_triangle(const uint8_t* trig, int i, int j, int k, int n, int v12) {

	size_t tv[3] = {};

	size_t startTriIdx = triangleCount.fetch_add(n, std::memory_order_relaxed);

	for (int t = 0; t < 3 * n; t++) {

		int edge = trig[t];

		switch (edge) {
			// x edges 
		case 0:
			tv[t % 3] = get_x_vert(i, j, k);
			break;
		case 1:
			tv[t % 3] = get_y_vert(i + 1, j, k);
			break;
		case 2:
			tv[t % 3] = get_x_vert(i, j + 1, k);
			break;
		case 3:
			tv[t % 3] = get_y_vert(i, j, k);
			break;
		case 4:
			tv[t % 3] = get_x_vert(i, j, k + 1);
			break;
		case 5:
			tv[t % 3] = get_y_vert(i + 1, j, k + 1);
			break;
		case 6:
			tv[t % 3] = get_x_vert(i, j + 1, k + 1);
			break;
		case 7:
			tv[t % 3] = get_y_vert(i, j, k + 1);
			break;
		case 8:
			tv[t % 3] = get_z_vert(i, j, k);
			break;
		case 9:
			tv[t % 3] = get_z_vert(i + 1, j, k);
			break;
		case 10:
			tv[t % 3] = get_z_vert(i + 1, j + 1, k);
			break;
		case 11:
			tv[t % 3] = get_z_vert(i, j + 1, k);
			break;
		case 12:
			tv[t % 3] = v12;
			break;
		default:
			break;
		}
		
		// create and add the Triangle to the list
		if (t % 3 == 2) {
			LTriangle triangle;
			triangle.v1 = tv[0];
			triangle.v2 = tv[1];
			triangle.v3 = tv[2];

			// add to the list
			meshTriangles[startTriIdx + (t / 3)] = triangle;
		}
	}
};

LVertex GeneratorLewiner::add_x_vertex(
	const int& i, const int& j, const int& k,
	const float& val0, const float& val1) {

	LVertex vertex;

	float u = val0 / (val0 - val1);

	vertex.x = bounds[0] + ((float)i + u) * stepX;
	vertex.y = bounds[2] + (float)j * stepY;
	vertex.z = bounds[4] + (float)k * stepZ;

	vertex.nx = (1 - u) * get_x_grad(i, j, k) + u * get_x_grad(i + 1, j, k);
	vertex.ny = (1 - u) * get_y_grad(i, j, k) + u * get_y_grad(i + 1, j, k);
	vertex.nz = (1 - u) * get_z_grad(i, j, k) + u * get_z_grad(i + 1, j, k);

	u = (float)std::sqrt(vertex.nx * vertex.nx + vertex.ny * vertex.ny + vertex.nz * vertex.nz);
	if (u > 0) {
		vertex.nx /= u;
		vertex.ny /= u;
		vertex.nz /= u;
	}

	return vertex;
};

LVertex GeneratorLewiner::add_y_vertex(
	const int& i, const int& j, const int& k,
	const float& val0, const float& val2) {

	LVertex vertex;

	float u = val0 / (val0 - val2);

	vertex.x = bounds[0] + (float)i * stepX;
	vertex.y = bounds[2] + ((float)j + u) * stepY;
	vertex.z = bounds[4] + (float)k * stepZ;

	vertex.nx = (1 - u) * get_x_grad(i, j, k) + u * get_x_grad(i, j + 1, k);
	vertex.ny = (1 - u) * get_y_grad(i, j, k) + u * get_y_grad(i, j + 1, k);
	vertex.nz = (1 - u) * get_z_grad(i, j, k) + u * get_z_grad(i, j + 1, k);

	u = (float)std::sqrt(vertex.nx * vertex.nx + vertex.ny * vertex.ny + vertex.nz * vertex.nz);
	if (u > 0) {
		vertex.nx /= u;
		vertex.ny /= u;
		vertex.nz /= u;
	}

	return vertex;
};

LVertex GeneratorLewiner::add_z_vertex(
	const int& i, const int& j, const int& k,
	const float& val0, const float& val3) {

	LVertex vertex;

	float u = val0 / (val0 - val3);

	vertex.x = bounds[0] + (float)i * stepX;
	vertex.y = bounds[2] + (float)j * stepY;
	vertex.z = bounds[4] + ((float)k + u) * stepZ;

	vertex.nx = (1 - u) * get_x_grad(i, j, k) + u * get_x_grad(i, j, k + 1);
	vertex.ny = (1 - u) * get_y_grad(i, j, k) + u * get_y_grad(i, j, k + 1);
	vertex.nz = (1 - u) * get_z_grad(i, j, k) + u * get_z_grad(i, j, k + 1);

	u = (float)std::sqrt(vertex.nx * vertex.nx + vertex.ny * vertex.ny + vertex.nz * vertex.nz);
	if (u > 0) {
		vertex.nx /= u;
		vertex.ny /= u;
		vertex.nz /= u;
	}

	return vertex;
};

float GeneratorLewiner::get_x_grad(const int i, const int j, const int k) {

	if (i > 0) {
		if (i < blockDims[0] - 1) {
			return(get_data(i + 1, j, k) - get_data(i - 1, j, k)) * 0.5f;
		}
		else {
			return get_data(i, j, k) - get_data(i - 1, j, k);
		}
	}
	else {
		return get_data(i + 1, j, k) - get_data(i, j, k);
	}

};

float GeneratorLewiner::get_y_grad(const int i, const int j, const int k) {
	if (j > 0) {
		if (j < blockDims[1] - 1)
			return (get_data(i, j + 1, k) - get_data(i, j - 1, k)) * 0.5f;
		else
			return get_data(i, j, k) - get_data(i, j - 1, k);
	}
	else
		return get_data(i, j + 1, k) - get_data(i, j, k);
};

float GeneratorLewiner::get_z_grad(const int i, const int j, const int k) {
	if (k > 0) {
		if (k < blockDims[2] - 1)
			return (get_data(i, j, k + 1) - get_data(i, j, k - 1)) * 0.5f;
		else
			return get_data(i, j, k) - get_data(i, j, k - 1);
	}
	else
		return get_data(i, j, k + 1) - get_data(i, j, k);
};

size_t GeneratorLewiner::get_x_vert(int i, int j, int k) {
	return x_verts[find_vertex_index(i, j, k)];
};

size_t GeneratorLewiner::get_y_vert(int i, int j, int k) {
	return y_verts[find_vertex_index(i, j, k)];
};

size_t GeneratorLewiner::get_z_vert(int i, int j, int k) {
	return z_verts[find_vertex_index(i, j, k)];
};

size_t GeneratorLewiner::add_c_vertex(const int i, const int j, const int k) {

	// add test_vertex_addition();

	LVertex vertex;

	// average
	float u = 0.0f;

	// set all to zero
	vertex.x = vertex.y = vertex.z = vertex.nx = vertex.ny = vertex.nz = 0;

	// estimate the average of the intersection points of the cube
	size_t vid = get_x_vert(i, j, k);
	if (vid != SIZE_MAX) {
		u++;
		const LVertex& v = meshVertices[vid];
		vertex.x += v.x;
		vertex.y += v.y;
		vertex.z += v.z;
		vertex.nx += v.nx;
		vertex.ny += v.ny;
		vertex.nz += v.nz;
	}
	vid = get_y_vert(i + 1, j, k);
	if (vid != SIZE_MAX) {
		u++;
		const LVertex& v = meshVertices[vid];
		vertex.x += v.x;
		vertex.y += v.y;
		vertex.z += v.z;
		vertex.nx += v.nx;
		vertex.ny += v.ny;
		vertex.nz += v.nz;
	}
	vid = get_x_vert(i, j + 1, k);
	if (vid != SIZE_MAX) {
		u++;
		const LVertex& v = meshVertices[vid];
		vertex.x += v.x;
		vertex.y += v.y;
		vertex.z += v.z;
		vertex.nx += v.nx;
		vertex.ny += v.ny;
		vertex.nz += v.nz;
	}
	vid = get_y_vert(i, j, k);
	if (vid != SIZE_MAX) {
		++u;
		const LVertex& v = meshVertices[vid];
		vertex.x += v.x;
		vertex.y += v.y;
		vertex.z += v.z;
		vertex.nx += v.nx;
		vertex.ny += v.ny;
		vertex.nz += v.nz;
	}
	vid = get_x_vert(i, j, k + 1);
	if (vid != SIZE_MAX) {
		++u;
		const LVertex& v = meshVertices[vid];
		vertex.x += v.x;
		vertex.y += v.y;
		vertex.z += v.z;
		vertex.nx += v.nx;
		vertex.ny += v.ny;
		vertex.nz += v.nz;
	}
	vid = get_y_vert(i + 1, j, k + 1);
	if (vid != SIZE_MAX) {
		++u;
		const LVertex& v = meshVertices[vid];
		vertex.x += v.x;
		vertex.y += v.y;
		vertex.z += v.z;
		vertex.nx += v.nx;
		vertex.ny += v.ny;
		vertex.nz += v.nz;
	}
	vid = get_x_vert(i, j + 1, k + 1);
	if (vid != SIZE_MAX) {
		++u;
		const LVertex& v = meshVertices[vid];
		vertex.x += v.x;
		vertex.y += v.y;
		vertex.z += v.z;
		vertex.nx += v.nx;
		vertex.ny += v.ny;
		vertex.nz += v.nz;
	}
	vid = get_y_vert(i, j, k + 1);
	if (vid != SIZE_MAX) {
		++u;
		const LVertex& v = meshVertices[vid];
		vertex.x += v.x;
		vertex.y += v.y;
		vertex.z += v.z;
		vertex.nx += v.nx;
		vertex.ny += v.ny;
		vertex.nz += v.nz;
	}
	vid = get_z_vert(i, j, k);
	if (vid != SIZE_MAX) {
		++u;
		const LVertex& v = meshVertices[vid];
		vertex.x += v.x;
		vertex.y += v.y;
		vertex.z += v.z;
		vertex.nx += v.nx;
		vertex.ny += v.ny;
		vertex.nz += v.nz;
	}
	vid = get_z_vert(i + 1, j, k);
	if (vid != SIZE_MAX) {
		++u;
		const LVertex& v = meshVertices[vid];
		vertex.x += v.x;
		vertex.y += v.y;
		vertex.z += v.z;
		vertex.nx += v.nx;
		vertex.ny += v.ny;
		vertex.nz += v.nz;
	}
	vid = get_z_vert(i + 1, j + 1, k);
	if (vid != SIZE_MAX) {
		++u;
		const LVertex& v = meshVertices[vid];
		vertex.x += v.x;
		vertex.y += v.y;
		vertex.z += v.z;
		vertex.nx += v.nx;
		vertex.ny += v.ny;
		vertex.nz += v.nz;
	}
	vid = get_z_vert(i, j + 1, k);
	if (vid != SIZE_MAX) {
		++u;
		const LVertex& v = meshVertices[vid];
		vertex.x += v.x;
		vertex.y += v.y;
		vertex.z += v.z;
		vertex.nx += v.nx;
		vertex.ny += v.ny;
		vertex.nz += v.nz;
	}

	if (u > 0) {
		vertex.x /= u;
		vertex.y /= u;
		vertex.z /= u;

		u = (float)std::sqrt(
			vertex.nx * vertex.nx + vertex.ny * vertex.ny + vertex.nz * vertex.nz);

		vertex.nx /= u;
		vertex.ny /= u;
		vertex.nz /= u;
	}

	// use atom to add the vertex into the correct slot
	size_t vIdx = vertexCount.fetch_add(1, std::memory_order_relaxed);
	meshVertices[vIdx] = vertex;

	return vIdx;
};

//-----------------------------------------------------------------------------
bool GeneratorLewiner::test_face(signed char face, const float _cube[8])
{
	float A, B, C, D;

	switch (face) {
	case -1:
	case 1:
		A = _cube[0];
		B = _cube[4];
		C = _cube[5];
		D = _cube[1];
		break;
	case -2:
	case 2:
		A = _cube[1];
		B = _cube[5];
		C = _cube[6];
		D = _cube[2];
		break;
	case -3:
	case 3:
		A = _cube[2];
		B = _cube[6];
		C = _cube[7];
		D = _cube[3];
		break;
	case -4:
	case 4:
		A = _cube[3];
		B = _cube[7];
		C = _cube[4];
		D = _cube[0];
		break;
	case -5:
	case 5:
		A = _cube[0];
		B = _cube[3];
		C = _cube[2];
		D = _cube[1];
		break;
	case -6:
	case 6:
		A = _cube[4];
		B = _cube[7];
		C = _cube[6];
		D = _cube[5];
		break;

	default:
		printf("Invalid face code %d\n", face);
		//print_cube();
		A = B = C = D = 0;
		break;
	};

	if (fabs(A * C - B * D) < FLT_EPSILON)
		return face >= 0;

	return face * A * (A * C - B * D) >= 0; // face and A invert signs
}

//-----------------------------------------------------------------------------
bool GeneratorLewiner::modified_test_interior(signed char s, int _case, int _config, const float _cube[8])
{

	char edge = -1;
	int amb_face;

	int inter_amb = 0;

	switch (_case) {
	case 4:

		amb_face = 1;
		edge = interior_ambiguity(amb_face, s, _cube);
		inter_amb += interior_ambiguity_verification(edge, _cube);

		amb_face = 2;
		edge = interior_ambiguity(amb_face, s, _cube);
		inter_amb += interior_ambiguity_verification(edge, _cube);

		amb_face = 5;
		edge = interior_ambiguity(amb_face, s, _cube);
		inter_amb += interior_ambiguity_verification(edge, _cube);

		if (inter_amb == 0) return false;
		else                return true;
		break;

	case 6:

		amb_face = abs((signed char)test6[_config][0]);

		edge = interior_ambiguity(amb_face, s, _cube);
		inter_amb = interior_ambiguity_verification(edge, _cube);

		if (inter_amb == 0) return false;
		else				return true;

		break;

	case 7:
		s = s * -1;

		amb_face = 1;
		edge = interior_ambiguity(amb_face, s, _cube);
		inter_amb += interior_ambiguity_verification(edge, _cube);

		amb_face = 2;
		edge = interior_ambiguity(amb_face, s, _cube);
		inter_amb += interior_ambiguity_verification(edge, _cube);

		amb_face = 5;
		edge = interior_ambiguity(amb_face, s, _cube);
		inter_amb += interior_ambiguity_verification(edge, _cube);

		if (inter_amb == 0) return false;
		else                return true;
		break;

	case 10:

		amb_face = abs((signed char)test10[_config][0]);

		edge = interior_ambiguity(amb_face, s, _cube);
		inter_amb = interior_ambiguity_verification(edge, _cube);

		if (inter_amb == 0) return false;
		else                return true;
		break;

	case 12:
		amb_face = abs((signed char)test12[_config][0]);
		edge = interior_ambiguity(amb_face, s, _cube);
		inter_amb += interior_ambiguity_verification(edge, _cube);


		amb_face = abs((signed char)test12[_config][1]);
		edge = interior_ambiguity(amb_face, s, _cube);
		inter_amb += interior_ambiguity_verification(edge, _cube);

		if (inter_amb == 0) return false;
		else				return true;
		break;
	}

}

//________________________________________________________________________________________________________
int GeneratorLewiner::interior_ambiguity(int amb_face, int s, const float _cube[8]) {
	int edge;

	switch (amb_face) {
	case 1:
	case 3:
		if (((_cube[1] * s) > 0) && ((_cube[7] * s) > 0))
			edge = 4;
		if (((_cube[0] * s) > 0) && ((_cube[6] * s) > 0))
			edge = 5;
		if (((_cube[3] * s) > 0) && ((_cube[5] * s) > 0))
			edge = 6;
		if (((_cube[2] * s) > 0) && ((_cube[4] * s) > 0))
			edge = 7;

		break;

	case 2:
	case 4:
		if (((_cube[1] * s) > 0) && ((_cube[7] * s) > 0))
			edge = 0;
		if (((_cube[2] * s) > 0) && ((_cube[4] * s) > 0))
			edge = 1;
		if (((_cube[3] * s) > 0) && ((_cube[5] * s) > 0))
			edge = 2;
		if (((_cube[0] * s) > 0) && ((_cube[6] * s) > 0))
			edge = 3;
		break;

	case 5:
	case 6:
	case 0:
		if (((_cube[0] * s) > 0) && ((_cube[6] * s) > 0))
			edge = 8;
		if (((_cube[1] * s) > 0) && ((_cube[7] * s) > 0))
			edge = 9;
		if (((_cube[2] * s) > 0) && ((_cube[4] * s) > 0))
			edge = 10;
		if (((_cube[3] * s) > 0) && ((_cube[5] * s) > 0))
			edge = 11;

		break;
	}

	return edge;
}
//-----------------------------------------------------------------------------
int GeneratorLewiner::interior_ambiguity_verification(int edge, const float _cube[8]) {
	float t, At = 0, Bt = 0, Ct = 0, Dt = 0, a = 0, b = 0;
	float verify;

	switch (edge) {

	case 0:
		a = (_cube[0] - _cube[1]) * (_cube[7] - _cube[6])
			- (_cube[4] - _cube[5]) * (_cube[3] - _cube[2]);
		b = _cube[6] * (_cube[0] - _cube[1]) + _cube[1] * (_cube[7] - _cube[6])
			- _cube[2] * (_cube[4] - _cube[5])
			- _cube[5] * (_cube[3] - _cube[2]);

		if (a > 0)
			return 1;

		t = -b / (2 * a);
		if (t < 0 || t > 1)
			return 1;

		At = _cube[1] + (_cube[0] - _cube[1]) * t;
		Bt = _cube[5] + (_cube[4] - _cube[5]) * t;
		Ct = _cube[6] + (_cube[7] - _cube[6]) * t;
		Dt = _cube[2] + (_cube[3] - _cube[2]) * t;

		verify = At * Ct - Bt * Dt;

		if (verify > 0)
			return 0;
		if (verify < 0)
			return 1;

		break;

	case 1:
		a = (_cube[3] - _cube[2]) * (_cube[4] - _cube[5])
			- (_cube[0] - _cube[1]) * (_cube[7] - _cube[6]);
		b = _cube[5] * (_cube[3] - _cube[2]) + _cube[2] * (_cube[4] - _cube[5])
			- _cube[6] * (_cube[0] - _cube[1])
			- _cube[1] * (_cube[7] - _cube[6]);

		if (a > 0)
			return 1;

		t = -b / (2 * a);
		if (t < 0 || t > 1)
			return 1;

		At = _cube[2] + (_cube[3] - _cube[2]) * t;
		Bt = _cube[1] + (_cube[0] - _cube[1]) * t;
		Ct = _cube[5] + (_cube[4] - _cube[5]) * t;
		Dt = _cube[6] + (_cube[7] - _cube[6]) * t;

		verify = At * Ct - Bt * Dt;

		if (verify > 0)
			return 0;
		if (verify < 0)
			return 1;
		break;

	case 2:
		a = (_cube[2] - _cube[3]) * (_cube[5] - _cube[4])
			- (_cube[6] - _cube[7]) * (_cube[1] - _cube[0]);
		b = _cube[4] * (_cube[2] - _cube[3]) + _cube[3] * (_cube[5] - _cube[4])
			- _cube[0] * (_cube[6] - _cube[7])
			- _cube[7] * (_cube[1] - _cube[0]);
		if (a > 0)
			return 1;

		t = -b / (2 * a);
		if (t < 0 || t > 1)
			return 1;

		At = _cube[3] + (_cube[2] - _cube[3]) * t;
		Bt = _cube[7] + (_cube[6] - _cube[7]) * t;
		Ct = _cube[4] + (_cube[5] - _cube[4]) * t;
		Dt = _cube[0] + (_cube[1] - _cube[0]) * t;

		verify = At * Ct - Bt * Dt;

		if (verify > 0)
			return 0;
		if (verify < 0)
			return 1;
		break;

	case 3:
		a = (_cube[1] - _cube[0]) * (_cube[6] - _cube[7])
			- (_cube[2] - _cube[3]) * (_cube[5] - _cube[4]);
		b = _cube[7] * (_cube[1] - _cube[0]) + _cube[0] * (_cube[6] - _cube[7])
			- _cube[4] * (_cube[2] - _cube[3])
			- _cube[3] * (_cube[5] - _cube[4]);
		if (a > 0)
			return 1;

		t = -b / (2 * a);
		if (t < 0 || t > 1)
			return 1;

		At = _cube[0] + (_cube[1] - _cube[0]) * t;
		Bt = _cube[3] + (_cube[2] - _cube[3]) * t;
		Ct = _cube[7] + (_cube[6] - _cube[7]) * t;
		Dt = _cube[4] + (_cube[5] - _cube[4]) * t;

		verify = At * Ct - Bt * Dt;

		if (verify > 0)
			return 0;
		if (verify < 0)
			return 1;
		break;

	case 4:

		a = (_cube[2] - _cube[1]) * (_cube[7] - _cube[4])
			- (_cube[3] - _cube[0]) * (_cube[6] - _cube[5]);
		b = _cube[4] * (_cube[2] - _cube[1]) + _cube[1] * (_cube[7] - _cube[4])
			- _cube[5] * (_cube[3] - _cube[0])
			- _cube[0] * (_cube[6] - _cube[5]);

		if (a > 0)
			return 1;

		t = -b / (2 * a);
		if (t < 0 || t > 1)
			return 1;

		At = _cube[1] + (_cube[2] - _cube[1]) * t;
		Bt = _cube[0] + (_cube[3] - _cube[0]) * t;
		Ct = _cube[4] + (_cube[7] - _cube[4]) * t;
		Dt = _cube[5] + (_cube[6] - _cube[5]) * t;

		verify = At * Ct - Bt * Dt;

		if (verify > 0)
			return 0;
		if (verify < 0)
			return 1;
		break;

	case 5:

		a = (_cube[3] - _cube[0]) * (_cube[6] - _cube[5])
			- (_cube[2] - _cube[1]) * (_cube[7] - _cube[4]);
		b = _cube[5] * (_cube[3] - _cube[0]) + _cube[0] * (_cube[6] - _cube[5])
			- _cube[4] * (_cube[2] - _cube[1])
			- _cube[1] * (_cube[7] - _cube[4]);
		if (a > 0)
			return 1;

		t = -b / (2 * a);
		if (t < 0 || t > 1)
			return 1;

		At = _cube[0] + (_cube[3] - _cube[0]) * t;
		Bt = _cube[1] + (_cube[2] - _cube[1]) * t;
		Ct = _cube[5] + (_cube[6] - _cube[5]) * t;
		Dt = _cube[4] + (_cube[7] - _cube[4]) * t;

		verify = At * Ct - Bt * Dt;

		if (verify > 0)
			return 0;
		if (verify < 0)
			return 1;
		break;

	case 6:
		a = (_cube[0] - _cube[3]) * (_cube[5] - _cube[6])
			- (_cube[4] - _cube[7]) * (_cube[1] - _cube[2]);
		b = _cube[6] * (_cube[0] - _cube[3]) + _cube[3] * (_cube[5] - _cube[6])
			- _cube[2] * (_cube[4] - _cube[7])
			- _cube[7] * (_cube[1] - _cube[2]);
		if (a > 0)
			return 1;

		t = -b / (2 * a);
		if (t < 0 || t > 1)
			return 1;

		At = _cube[3] + (_cube[0] - _cube[3]) * t;
		Bt = _cube[7] + (_cube[4] - _cube[7]) * t;
		Ct = _cube[6] + (_cube[5] - _cube[6]) * t;
		Dt = _cube[2] + (_cube[1] - _cube[2]) * t;

		verify = At * Ct - Bt * Dt;

		if (verify > 0)
			return 0;
		if (verify < 0)
			return 1;
		break;

	case 7:
		a = (_cube[1] - _cube[2]) * (_cube[4] - _cube[7])
			- (_cube[0] - _cube[3]) * (_cube[5] - _cube[6]);
		b = _cube[7] * (_cube[1] - _cube[2]) + _cube[2] * (_cube[4] - _cube[7])
			- _cube[6] * (_cube[0] - _cube[3])
			- _cube[3] * (_cube[5] - _cube[6]);
		if (a > 0)
			return 1;

		t = -b / (2 * a);
		if (t < 0 || t > 1)
			return 1;

		At = _cube[2] + (_cube[1] - _cube[2]) * t;
		Bt = _cube[3] + (_cube[0] - _cube[3]) * t;
		Ct = _cube[7] + (_cube[4] - _cube[7]) * t;
		Dt = _cube[6] + (_cube[5] - _cube[6]) * t;

		verify = At * Ct - Bt * Dt;

		if (verify > 0)
			return 0;
		if (verify < 0)
			return 1;
		break;

	case 8:
		a = (_cube[4] - _cube[0]) * (_cube[6] - _cube[2])
			- (_cube[7] - _cube[3]) * (_cube[5] - _cube[1]);
		b = _cube[2] * (_cube[4] - _cube[0]) + _cube[0] * (_cube[6] - _cube[2])
			- _cube[1] * (_cube[7] - _cube[3])
			- _cube[3] * (_cube[5] - _cube[1]);
		if (a > 0)
			return 1;

		t = -b / (2 * a);
		if (t < 0 || t > 1)
			return 1;

		At = _cube[0] + (_cube[4] - _cube[0]) * t;
		Bt = _cube[3] + (_cube[7] - _cube[3]) * t;
		Ct = _cube[2] + (_cube[6] - _cube[2]) * t;
		Dt = _cube[1] + (_cube[5] - _cube[1]) * t;

		verify = At * Ct - Bt * Dt;

		if (verify > 0)
			return 0;
		if (verify < 0)
			return 1;
		break;

	case 9:
		a = (_cube[5] - _cube[1]) * (_cube[7] - _cube[3])
			- (_cube[4] - _cube[0]) * (_cube[6] - _cube[2]);
		b = _cube[3] * (_cube[5] - _cube[1]) + _cube[1] * (_cube[7] - _cube[3])
			- _cube[2] * (_cube[4] - _cube[0])
			- _cube[0] * (_cube[6] - _cube[2]);
		if (a > 0)
			return 1;

		t = -b / (2 * a);
		if (t < 0 || t > 1)
			return 1;

		At = _cube[1] + (_cube[5] - _cube[1]) * t;
		Bt = _cube[0] + (_cube[4] - _cube[0]) * t;
		Ct = _cube[3] + (_cube[7] - _cube[3]) * t;
		Dt = _cube[2] + (_cube[6] - _cube[2]) * t;

		verify = At * Ct - Bt * Dt;

		if (verify > 0)
			return 0;
		if (verify < 0)
			return 1;
		break;

	case 10:
		a = (_cube[6] - _cube[2]) * (_cube[4] - _cube[0])
			- (_cube[5] - _cube[1]) * (_cube[7] - _cube[3]);
		b = _cube[0] * (_cube[6] - _cube[2]) + _cube[2] * (_cube[4] - _cube[0])
			- _cube[3] * (_cube[5] - _cube[1])
			- _cube[1] * (_cube[7] - _cube[3]);
		if (a > 0)
			return 1;

		t = -b / (2 * a);
		if (t < 0 || t > 1)
			return 1;

		At = _cube[2] + (_cube[6] - _cube[2]) * t;
		Bt = _cube[1] + (_cube[5] - _cube[1]) * t;
		Ct = _cube[0] + (_cube[4] - _cube[0]) * t;
		Dt = _cube[3] + (_cube[7] - _cube[3]) * t;

		verify = At * Ct - Bt * Dt;

		if (verify > 0)
			return 0;
		if (verify < 0)
			return 1;
		break;

	case 11:
		a = (_cube[7] - _cube[3]) * (_cube[5] - _cube[1])
			- (_cube[6] - _cube[2]) * (_cube[4] - _cube[0]);
		b = _cube[1] * (_cube[7] - _cube[3]) + _cube[3] * (_cube[5] - _cube[1])
			- _cube[0] * (_cube[6] - _cube[2])
			- _cube[2] * (_cube[4] - _cube[0]);
		if (a > 0)
			return 1;

		t = -b / (2 * a);
		if (t < 0 || t > 1)
			return 1;

		At = _cube[3] + (_cube[7] - _cube[3]) * t;
		Bt = _cube[2] + (_cube[6] - _cube[2]) * t;
		Ct = _cube[1] + (_cube[5] - _cube[1]) * t;
		Dt = _cube[0] + (_cube[4] - _cube[0]) * t;

		verify = At * Ct - Bt * Dt;

		if (verify > 0)
			return 0;
		if (verify < 0)
			return 1;
		break;
	}

	return 0;
};


//_____________________________________________________________________________
// NEWS INTERIOR TEST FOR CASE 13
// Return true if the interior is empty(two faces)

bool GeneratorLewiner::interior_test_case13(const float _cube[8])
{
	float t1, t2, At1 = 0, Bt1 = 0, Ct1 = 0, Dt1 = 0, At2 = 0, Bt2 = 0, Ct2 = 0, Dt2 = 0, a = 0, b = 0, c = 0;

	a = (_cube[0] - _cube[1]) * (_cube[7] - _cube[6])
		- (_cube[4] - _cube[5]) * (_cube[3] - _cube[2]);
	b = _cube[6] * (_cube[0] - _cube[1]) + _cube[1] * (_cube[7] - _cube[6])
		- _cube[2] * (_cube[4] - _cube[5])
		- _cube[5] * (_cube[3] - _cube[2]);

	c = _cube[1] * _cube[6] - _cube[5] * _cube[2];

	double delta = b * b - 4 * a * c;

	t1 = (-b + sqrt(delta)) / (2 * a);
	t2 = (-b - sqrt(delta)) / (2 * a);

	printf("t1 = %f, t2 = %f\n", t1, t2);

	if ((t1 < 1) && (t1 > 0) && (t2 < 1) && (t2 > 0))
	{

		At1 = _cube[1] + (_cube[0] - _cube[1]) * t1;
		Bt1 = _cube[5] + (_cube[4] - _cube[5]) * t1;
		Ct1 = _cube[6] + (_cube[7] - _cube[6]) * t1;
		Dt1 = _cube[2] + (_cube[3] - _cube[2]) * t1;

		float x1 = (At1 - Dt1) / (At1 + Ct1 - Bt1 - Dt1);
		float y1 = (At1 - Bt1) / (At1 + Ct1 - Bt1 - Dt1);

		At2 = _cube[1] + (_cube[0] - _cube[1]) * t2;
		Bt2 = _cube[5] + (_cube[4] - _cube[5]) * t2;
		Ct2 = _cube[6] + (_cube[7] - _cube[6]) * t2;
		Dt2 = _cube[2] + (_cube[3] - _cube[2]) * t2;

		float x2 = (At2 - Dt2) / (At2 + Ct2 - Bt2 - Dt2);
		float y2 = (At2 - Bt2) / (At2 + Ct2 - Bt2 - Dt2);

		if ((x1 < 1) && (x1 > 0) && (x2 < 1) && (x2 > 0) && (y1 < 1) && (y1 > 0) && (y2 < 1) && (y2 > 0))
			return false;

		return true;
	}

	else return true;

}
//--------------------------------------------------------------------------------------------------------------------
// control the tunnel orientation triangulation
bool GeneratorLewiner::interior_test_case13_2(float isovalue, const float _cube[8], int& tunnelOrientation) {

	double critival_point_value1, critival_point_value2;

	double a = -_cube[0] + _cube[1] + _cube[3] - _cube[2] + _cube[4] - _cube[5] - _cube[7] + _cube[6],
		b = _cube[0] - _cube[1] - _cube[3] + _cube[2],
		c = _cube[0] - _cube[1] - _cube[4] + _cube[5],
		d = _cube[0] - _cube[3] - _cube[4] + _cube[7],
		e = -_cube[0] + _cube[1],
		f = -_cube[0] + _cube[3],
		g = -_cube[0] + _cube[4],
		h = _cube[0];

	double x1, y1, z1, x2, y2, z2;
	int numbercritivalpoints = 0;

	double dx = b * c - a * e, dy = b * d - a * f, dz = c * d - a * g;

	if (dx != 0.0f && dy != 0.0f && dz != 0.0f) {
		if (dx * dy * dz < 0)
			return true;

		double disc = sqrt(dx * dy * dz);

		x1 = (-d * dx - disc) / (a * dx);
		y1 = (-c * dy - disc) / (a * dy);
		z1 = (-b * dz - disc) / (a * dz);

		if ((x1 > 0) && (x1 < 1) && (y1 > 0) && (y1 < 1)
			&& (z1 > 0) && (z1 < 1)) {
			numbercritivalpoints++;

			critival_point_value1 = a * x1 * y1 * z1 + b * x1 * y1 + c * x1 * z1
				+ d * y1 * z1 + e * x1 + f * y1 + g * z1 + h - isovalue;
		}

		x2 = (-d * dx + disc) / (a * dx);
		y2 = (-c * dy + disc) / (a * dy);
		z2 = (-b * dz + disc) / (a * dz);

		if ((x2 > 0) && (x2 < 1) && (y2 > 0) && (y2 < 1)
			&& (z2 > 0) && (z2 < 1)) {
			numbercritivalpoints++;

			critival_point_value2 = a * x2 * y2 * z2 + b * x2 * y2 + c * x2 * z2
				+ d * y2 * z2 + e * x2 + f * y2 + g * z2 + h - isovalue;

		}

		if (numbercritivalpoints < 2)
			return true;
		else
		{
			if ((critival_point_value1 * critival_point_value2 > 0))
			{
				if (critival_point_value1 > 0)
					tunnelOrientation = 1;
				else
					tunnelOrientation = -1;
			}

			return critival_point_value1 * critival_point_value2 < 0;
		}

	}
	else
		return true;
}

void GeneratorLewiner::seal_grid_boundaries() {
	// Force the 1-voxel thick outer shell of the grid to be empty space (air)
	// This forces Marching Cubes to cap off any open holes and make the mesh watertight.

	for (int i = 0; i < blockDims[0]; i++) {
		for (int j = 0; j < blockDims[1]; j++) {
			for (int k = 0; k < blockDims[2]; k++) {
				// If we are on the outermost boundary...
				if (i == 0 || i == blockDims[0] - 1 ||
					j == 0 || j == blockDims[1] - 1 ||
					k == 0 || k == blockDims[2] - 1)
				{
					scalarField[find_vertex_index(i, j, k)] = 9999.9f;
				}
			}
		}
	}
}

//@brief function to compute a smooth minimum between two values, 
// this is used to create a smooth transition between the distance fields of the seeds,
// and to create a smoother mesh. [Inigo Quilez]
//@param a: the first value
//@param b: the second value
//@param k: the fillet radius (in field units). Blending is active only where
// |a-b| < k, so LARGER k gives a WIDER, smoother transition and k -> 0 reduces
// exactly to the hard min (a sharp crease).
//@returns the smooth minimum between a and b
	SmoothDist GeneratorLewiner::smin_gradient(	
		float a, float b, 
		const Vec3& gradA, const Vec3& gradB, float k) {

	float diff = std::abs(a-b);

	if (diff >= k){
		if (a < b) return {a , gradA};
		else return {b, gradB};

	}

	float h = (k - diff) / k;
	float m = 0.5f * h * h;

	SmoothDist result;
	result.val = std::min(a, b) - h * h * h * k * (1.0f / 6.0f);

	// estimate gradient
	if (a < b) 	result.grad = gradA * (1.0f - m) + gradB * m;
	else result.grad = m * gradA + (1.0f - m) * gradB;

	return result;
}

void GeneratorLewiner::update_steps() {
	stepX = (bounds[1] - bounds[0]) / (blockDims[0] - 1);
	stepY = (bounds[3] - bounds[2]) / (blockDims[1] - 1);
	stepZ = (bounds[5] - bounds[4]) / (blockDims[2] - 1);
}

float GeneratorLewiner::air_skip_level() const {
	// Two voxels above the iso-level: wide enough to cover the reach of the
	// 6-/26-neighbour smoothing stencils, narrow enough not to erode the
	// solid. Scales with the grid spacing so it behaves the same whether the
	// model is 0.5 mm or 500 mm across.
	const float hMax = std::max(stepX, std::max(stepY, stepZ));
	return isoLevel + 2.0f * hMax;
}

void GeneratorLewiner::_setup_mesh() {

	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);
	glGenBuffers(1, &normalsVBO);

	// We just bind the VAO configuration here, but we don't upload data yet
	glBindVertexArray(VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	// Vertex attributes (pointers) are state of the VAO, so we define them once.
	// Note: It's safe to define pointers even if buffer is empty, as long as VBO is bound.
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// normals
	glBindBuffer(GL_ARRAY_BUFFER, normalsVBO);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(1);

	// Bind EBO to VAO state
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);

	glBindVertexArray(0);

};

void GeneratorLewiner::_setup_edges() {

	// generate Vertex Array Object
	glGenVertexArrays(1, &edgeVAO);
	// generate Vertex Buffer Object to store vertex attributes
	glGenBuffers(1, &edgeVBO);
	// generate element buffer Object
	glGenBuffers(1, &edgeEBO);

	// bind to vertex array
	glBindVertexArray(edgeVAO);

	// bind array buffer and send data
	glBindBuffer(GL_ARRAY_BUFFER, edgeVBO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

	// bind element array buffer and send data
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edgeEBO);

	glEnableVertexAttribArray(0);
};

void GeneratorLewiner::draw() {

	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
};

void GeneratorLewiner::draw_edges() {
	glBindVertexArray(edgeVAO);
	glDrawElements(GL_LINES, edgeIndices.size(), GL_UNSIGNED_INT, 0);
};

void GeneratorLewiner::update_render() {

	// update vertices
	vertices.clear();
	normals.clear();
	indices.clear();
	edgeSet.clear();
	edgeIndices.clear();

	for (const auto& v : meshVertices) {
		vertices.push_back(v.x);
		vertices.push_back(v.y);
		vertices.push_back(v.z);

		// add normals
		normals.push_back(v.nx);
		normals.push_back(v.ny);
		normals.push_back(v.nz);

	}

	// use the triangles to add edges
	for (const auto& tri : meshTriangles) {

		indices.push_back(tri.v1);
		indices.push_back(tri.v2);
		indices.push_back(tri.v3);

		add_edge(tri.v1, tri.v2, tri.v3);
	}

	update_buffers();

};

void GeneratorLewiner::update_buffers(){

glBindVertexArray(VAO); 

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, normalsVBO);
	glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(float), normals.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

	glBindVertexArray(0);

	// update edges
	glBindBuffer(GL_ARRAY_BUFFER, edgeVBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, edgeEBO);
	glBufferData(GL_ARRAY_BUFFER, edgeIndices.size() * sizeof(int), edgeIndices.data(), GL_STATIC_DRAW);

	glBindVertexArray(0);

};

void GeneratorLewiner::add_edge(int idx1, int idx2, int idx3) {

	std::pair<unsigned int, unsigned int> edges[3] = {
		{ std::min(idx1, idx2), std::max(idx1, idx2) },
		{ std::min(idx2, idx3), std::max(idx2, idx3) },
		{ std::min(idx3, idx1), std::max(idx3, idx1) }
	};

	for (auto& e : edges) {
		// insert into edgeSet (guarantees uniqueness)
		if (edgeSet.insert(e).second) {
			edgeIndices.push_back(e.first);
			edgeIndices.push_back(e.second);
		}
	}
};

void GeneratorLewiner::export_nrrd(const std::string fileName, float voxelSize, std::array<float, 6> blockSize) {

	// get a new field
	std::vector<uint8_t> field = get_image_field(voxelSize, blockSize, false, 255);

	// Calculate required Grid Dimensions in the same way as the get_image_field function
	int nx = static_cast<int>(std::ceil((blockSize[1] - blockSize[0]) / voxelSize));
	int ny = static_cast<int>(std::ceil((blockSize[3] - blockSize[2]) / voxelSize));
	int nz = static_cast<int>(std::ceil((blockSize[5] - blockSize[4]) / voxelSize));

	//std::cout << "Starting High-Res Export..." << std::endl;
	//std::cout << "Physical Box: " << (blockSize[1] - blockSize[0]) << " mm" << std::endl;
	//std::cout << "Target Voxel: " << voxelSize << " mm" << std::endl;
	//std::cout << "Grid Size: " << nx << " x " << ny << " x " << nz << std::endl;

	// 4. Write NRRD
	std::ofstream file(fileName, std::ios::binary);
	file.imbue(std::locale::classic()); // Force "." decimals

	file << "NRRD0004\n";
	file << "# High-Res Scaffold Export\n";
	file << "type: unsigned char\n";
	file << "dimension: 3\n";
	file << "sizes: " << nx << " " << ny << " " << nz << "\n";

	// Explicitly write spacing in mm
	file << "spacings: " << std::fixed << std::setprecision(5)
		<< voxelSize << " " << voxelSize << " " << voxelSize << "\n";

	file << "units: \"mm\" \"mm\" \"mm\"\n"; // Explicit Units
	file << "encoding: raw\n";
	file << "\n"; // End of header

	file.write(reinterpret_cast<char*>(field.data()), field.size());
	file.close();

	logger->log(LogPriority::SUCCESS, "Export to Nrrd Completed: " + fileName);
};

void GeneratorLewiner::export_mhd(std::filesystem::path& path, float voxelSize, std::array<float, 6> blockBounds) {

	std::filesystem::path basePath = std::filesystem::absolute(path);
	basePath.replace_extension(""); // Strip extension to be safe

	std::string rawFileName = basePath.string() + ".raw";
	std::string mhdFileName = basePath.string() + ".mhd";

	// We need just the filename for the header (no full path)
	std::string rawBaseName = std::filesystem::path(rawFileName).filename().string();

	std::vector<uint8_t> field = get_image_field(voxelSize, blockBounds, false, 255);

	// Calculate required Grid Dimensions in the same way as the get_image_field function
	int nx = static_cast<int>(std::ceil((blockBounds[1] - blockBounds[0]) / voxelSize));
	int ny = static_cast<int>(std::ceil((blockBounds[3] - blockBounds[2]) / voxelSize));
	int nz = static_cast<int>(std::ceil((blockBounds[5] - blockBounds[4]) / voxelSize));

	std::ofstream rawFile(rawFileName, std::ios::binary);
	if (!rawFile) {
		logger->log(LogPriority::ERROR, "Error: Could not open " + rawFileName);
		return;
	}
	rawFile.write(reinterpret_cast<char*>(field.data()), field.size());
	rawFile.close();

	std::ofstream mhdFile(mhdFileName);
	if (!mhdFile) {
		logger->log(LogPriority::ERROR, "Error: Could not open " + mhdFileName);
		return;
	}
	mhdFile.imbue(std::locale::classic());

	mhdFile << "ObjectType = Image\n";
	mhdFile << "NDims = 3\n";
	mhdFile << "BinaryData = True\n";
	mhdFile << "BinaryDataByteOrderMSB = False\n";
	mhdFile << "CompressedData = False\n";
	mhdFile << "TransformMatrix = 1 0 0 0 1 0 0 0 1\n";
	mhdFile << "Offset = " << blockBounds[0] << " " << blockBounds[2] << " " << blockBounds[4] << "\n";
	mhdFile << "CenterOfRotation = 0 0 0\n";
	mhdFile << "ElementSpacing = " << std::fixed << std::setprecision(6)
		<< voxelSize << " " << voxelSize << " " << voxelSize << "\n";
	mhdFile << "DimSize = " << nx << " " << ny << " " << nz << "\n";
	mhdFile << "ElementType = MET_UCHAR\n";
	mhdFile << "ElementDataFile = " << rawBaseName << "\n";
	mhdFile.close();

	logger->log(LogPriority::SUCCESS, "Export to Mhd Completed: " + mhdFileName);

};

void GeneratorLewiner::export_stl(std::string fileName) {


	std::ofstream out(fileName, std::ios::binary);

	// header (80 bytes of zeroes/comments)
	char header[80] = { 0 };
	out.write(header, 80);

	// number of triangles (4 bytes)
	uint32_t numTriangles = static_cast<uint32_t>(meshTriangles.size());
	out.write(reinterpret_cast<const char*>(&numTriangles), sizeof(uint32_t));

	// pass triangle data
	for (const auto& tri : meshTriangles) {

		// Fetch the actual vertices using the integer indices
		const LVertex& p1 = meshVertices[tri.v1];
		const LVertex& p2 = meshVertices[tri.v2];
		const LVertex& p3 = meshVertices[tri.v3];

		// Calculate the Face Normal (Cross Product of Edge 1 and Edge 2)
		float e1x = p2.x - p1.x;
		float e1y = p2.y - p1.y;
		float e1z = p2.z - p1.z;

		float e2x = p3.x - p1.x;
		float e2y = p3.y - p1.y;
		float e2z = p3.z - p1.z;

		float nx = (e1y * e2z) - (e1z * e2y);
		float ny = (e1z * e2x) - (e1x * e2z);
		float nz = (e1x * e2y) - (e1y * e2x);

		// Normalize the face normal
		float len = std::sqrt(nx * nx + ny * ny + nz * nz);
		if (len > 0) {
			nx /= len;
			ny /= len;
			nz /= len;
		}

		// 3. Populate the 12-float STL data block
		float data[12];

		// Face Normal
		data[0] = nx;   data[1] = ny;   data[2] = nz;

		// Vertex 1
		data[3] = p1.x; data[4] = p1.y; data[5] = p1.z;
		// Vertex 2
		data[6] = p2.x; data[7] = p2.y; data[8] = p2.z;
		// Vertex 3
		data[9] = p3.x; data[10] = p3.y; data[11] = p3.z;

		// Write the 48 bytes (12 floats)
		out.write(reinterpret_cast<const char*>(data), 12 * sizeof(float));

		// Write the 2-byte attribute catch (usually just 0 for standard STLs)
		uint16_t attr = 0;
		out.write(reinterpret_cast<const char*>(&attr), sizeof(uint16_t));
	}

	out.close();
	logger->log(LogPriority::SUCCESS, "Exported stl as " + fileName);
	
	// std::filesystem::path parent = std::filesystem::path(fileName).parent_path();
	// std::string stlName = std::filesystem::path(fileName).stem().string();
	// std::string parameterFileName = stlName + "_parameters.csv";
	// std::filesystem::path parameterPath = parent / parameterFileName;

	// export_parameters(parameterPath.string());
	// export_metrics(parameterPath)
};

void GeneratorLewiner::validate_topology(bool supress) {

	// Map to count how many faces share each edge
	std::map<std::pair<int, int>, int> edgeCounts;

	for (const auto& tri : meshTriangles) {
		// Create 3 edges for each triangle. 
		// We use std::min and std::max so edge (A,B) and edge (B,A) are treated as the same edge.
		std::pair<int, int> e1 = { std::min(tri.v1, tri.v2), std::max(tri.v1, tri.v2) };
		std::pair<int, int> e2 = { std::min(tri.v2, tri.v3), std::max(tri.v2, tri.v3) };
		std::pair<int, int> e3 = { std::min(tri.v3, tri.v1), std::max(tri.v3, tri.v1) };

		// Increment the usage count for each edge
		edgeCounts[e1]++;
		edgeCounts[e2]++;
		edgeCounts[e3]++;
	}

	int boundaryEdges = 0;
	int nonManifoldEdges = 0;
	int manifoldEdges = 0;

	// Evaluate the counts
	for (const auto& pair : edgeCounts) {
		if (pair.second == 1) {
			boundaryEdges++;
		}
		else if (pair.second == 2) {
			manifoldEdges++;
		}
		else if (pair.second > 2) {
			nonManifoldEdges++;
		}
	}

	//std::cout << "Vertices: " << meshVertices.size() << std::endl;
	//std::cout << "Faces: " << meshTriangles.size() << std::endl;
	//std::cout << "Total Unique Edges: " << edgeCounts.size() << std::endl;
	//std::cout << "Boundary Edges (Holes): " << boundaryEdges << std::endl;
	//std::cout << "Non-Manifold Edges: " << nonManifoldEdges << std::endl;

	if (supress) return;
	if (boundaryEdges == 0 && nonManifoldEdges == 0) {
		logger->log(LogPriority::SUCCESS, "Mesh is 100% Watertight and 2-Manifold!");
	}
	else {
		logger->log(LogPriority::ERROR, "Mesh has topological errors.");
	}
}

// void GeneratorLewiner::render_properties(bool& updateScaffold, GenerationTask* task) 
// {
// 	std::shared_ptr<IContainer> lockedCon = container.lock();
// 	std::shared_ptr<InterfaceSeedGenerator> lockedGen = generator.lock();

// 	ImGui::ColorEdit4("Appearance", (float*)&color);

// 	// first render the applied generator and container
// 	if (lockedCon) {
// 		ImGui::Text("Container: %s", lockedCon->name.c_str());
// 	}
// 	if (lockedGen) {
// 		ImGui::Text("Generator: %s", lockedGen->name.c_str());
// 	}

// 	// --------------------------------------------------------------------------------
// 	ImGui::SeparatorText("Thickness");

// 	ImGui::RadioButton("Apply Uniform Thickness", &selectedThicknessOption, 0);
// 	ImGui::RadioButton("Apply Varied Thickness", &selectedThicknessOption, 1);

// 	// this is the uniform case
// 	if (selectedThicknessOption == 0) {
// 		ImGui::SetNextItemWidth(200);
// 		ImGui::InputFloat("Thickness", &isoLevel);
// 	}
// 	// this is the varied case
// 	else {
// 		ImGui::SetNextItemWidth(200);
// 		ImGui::InputFloat("Start Thickness", &startThickness);
// 		ImGui::SetNextItemWidth(200);
// 		ImGui::InputFloat("End Thickness", &endThickness);
// 		ImGui::SetNextItemWidth(200);
// 		ImGui::InputFloat("Transition Distance", &transitionDistance);

// 		ImGui::SeparatorText("Select Distance Function");
// 		ImGui::RadioButton("Distance From Plane", &selectedDist, 0);
// 		if (selectedDist == 0) {
// 			ImGui::SetNextItemWidth(200);
// 			ImGui::InputFloat3("Normal", distancePlaneNormal);
// 			ImGui::SetNextItemWidth(200);
// 			ImGui::InputFloat3("Center", distancePlaneCenter);
// 		};
// 		ImGui::RadioButton("Distance From Point", &selectedDist, 1);
// 		if (selectedDist == 1) {
// 			ImGui::SetNextItemWidth(200);
// 			ImGui::InputFloat3("Point", distancePoint);
// 		}
// 		ImGui::RadioButton("Distance From Container", &selectedDist, 2);

// 		ImGui::SeparatorText("Select Radius Function");
// 		ImGui::RadioButton("Linear", &selectedFunc, 0);
// 		ImGui::RadioButton("Quadratic", &selectedFunc, 1);
// 		ImGui::RadioButton("Constant", &selectedFunc, 2);
// 		ImGui::RadioButton("Random", &selectedFunc, 3);
// 	}

// 	ImGui::InputFloat("Voxel Size", &voxelSize);

// 	// other parameters -----------------------------------------------------------
// 	ImGui::SeparatorText("Parameters");
	
// 	ImGuiTableFlags flags = ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_BordersInnerV;
// 	// create a table
// 	if (ImGui::BeginTable("", 2, flags = flags)) {
		
// 		if (lockedGen) {
// 			ImGui::TableNextRow();

// 			if (lockedGen->get_type() == ObjectType::RandomGeneratorType) {
// 				ImGui::TableNextColumn(); ImGui::Text("Random Seeds");
// 				ImGui::TableNextColumn();
// 				ImGui::Text("%d", lockedGen->get_seeds().size());
// 			}
// 			else if (lockedGen->get_type() == ObjectType::PoissonGeneratorType) {
// 				ImGui::TableNextColumn(); ImGui::Text("Poisson 3D");
// 				ImGui::TableNextColumn();

// 				Poisson3D* dummy = static_cast<Poisson3D*>(lockedGen.get());
// 				if (dummy->is_uniform()) {
// 					ImGui::Text("Radius (Uniform) %.4f", dummy->get_min_radius());
// 				}
// 				else {
// 					ImGui::BeginTable("##", 2);
// 					ImGui::TableNextRow();
// 					ImGui::TableNextColumn();
// 					ImGui::Text("Rmin %4.f", dummy->get_min_radius());
// 					ImGui::TableNextColumn();
// 					ImGui::Text("Rmax %4.f", dummy->get_max_radius());
// 					ImGui::EndTable();
// 				}
// 			}
// 		}
// 		else {
// 			ImGui::TableNextRow();
// 			ImGui::TableNextColumn(); ImGui::Text("Source");
// 			ImGui::TableNextColumn(); ImGui::TextColored(ImVec4(0.4f, 0.8f, 0.4f, 1.0f), "Loaded from CSV");
// 		}
		
// 		ImGui::TableNextRow();
// 		ImGui::TableNextColumn(); ImGui::Text("Openess");
// 		ImGui::TableNextColumn();
// 		ImGui::SliderFloat("##Openess", &threshold, 0.0f, 1.0f, "%.3f");

// 		ImGui::TableNextRow();
// 		ImGui::TableNextColumn(); ImGui::Text("Stretch X");
// 		ImGui::TableNextColumn();
// 		ImGui::InputFloat("##Stretch X", &stretchX, 0.01f, 100.0f, "%.3f");

// 		ImGui::TableNextRow();
// 		ImGui::TableNextColumn(); ImGui::Text("Stretch Y");
// 		ImGui::TableNextColumn();
// 		ImGui::InputFloat("##Stretch Y", &stretchY, 0.01f, 100.0f, "%.3f");

// 		ImGui::TableNextRow();
// 		ImGui::TableNextColumn(); ImGui::Text("Stretch Z");
// 		ImGui::TableNextColumn();
// 		ImGui::InputFloat("##Stretch Z", &stretchZ, 0.01f, 100.0f, "%.3f");

// 		ImGui::TableNextRow();
// 		ImGui::TableNextColumn(); ImGui::Text("Material Direction");
// 		ImGui::TableNextColumn();
// 		ImGui::SetNextItemWidth(200.0f);
// 		ImGui::InputFloat3("##Material Direction", anisotropyVec, "%.4f");

// 		ImGui::TableNextRow();
// 		ImGui::TableNextColumn(); ImGui::Text("Angle");
// 		ImGui::TableNextColumn();
// 		ImGui::InputFloat("##Angle", &anisotropyAngle, 0.01f, 10.0f, "%.4f");
	
// 		ImGui::EndTable();
// 	}
	
// 	ImGui::SeparatorText("Metrics");

// 	render_metrics();

// 	// if the user pressed the update button from the gui
// 	if (updateScaffold) {
// 		if (isLoadedFromFile) {
// 			logger->log(LogPriority::WARNING, "Cannot update a static mesh loaded from a file.");
// 			updateScaffold = false;
// 		}
// 		if (lockedCon && lockedGen) {

// 			auto startTime = std::chrono::steady_clock::now();

// 			std::vector<Vec3> seeds = lockedGen->get_seeds();
// 			Bounds bds = lockedCon->compute_bounds();

// 			std::array<float, 6> bounds = {
// 				bds.xMin,
// 				bds.xMax,
// 				bds.yMin,
// 				bds.yMax,
// 				bds.zMin,
// 				bds.zMax
// 			};

// 			// check the resolution
// 			blockDims = {
// 				static_cast<int>(std::ceil((bds.xMax - bds.xMin) / voxelSize)) + 1,
// 				static_cast<int>(std::ceil((bds.yMax - bds.yMin) / voxelSize)) + 1,
// 				static_cast<int>(std::ceil((bds.zMax - bds.zMin) / voxelSize)) + 1
// 			};

// 			if (selectedThicknessOption == 1) {

// 				switch (selectedFunc) {
// 					// linear radius function
// 					case 0: {
// 						thicknessFunction = std::make_shared<LinearFunction>(transitionDistance);
// 						break;
// 					}
// 					case 1: {
// 						thicknessFunction = std::make_shared<QuadraticFunction>(transitionDistance);
// 						break;
// 					}
// 					case 2: {
// 						thicknessFunction = std::make_shared<ConstantRadiusFunction>();
// 						break;
// 					}
// 					case 3: {
// 						thicknessFunction = std::make_shared<RandomRadiusFunction>();
// 					}
// 				}

// 				switch (selectedDist) {
// 					// distance from plane
// 					case 0: {
// 						thicknessSDF = std::make_shared<PlaneSDF>(distancePlaneCenter, distancePlaneNormal);
// 						break;
// 					}
// 						  // distance from point
// 					case 1: {
// 						thicknessSDF = std::make_shared<PointSDF>(distancePoint);
// 						break;
// 					}
// 						  // distance from container surface
// 					case 2: {
// 						thicknessSDF = lockedCon->get_distance_estimator();
// 						break;
// 					}
// 				}
// 			}
// 			else {
// 				thicknessSDF.reset();
// 				thicknessFunction.reset();
// 			}

// 			set_thickness_functions(thicknessSDF, thicknessFunction, startThickness, endThickness, transitionDistance);

// 			tortuosityPathModel.reset();
// 			set_bounds(bounds);
// 			set_seeds(seeds);
			
// 			//std::cout << "uniform flag" << selectedThicknessOption << " minT: " << startThickness << " maxT: " << maxThickness << " tDist: " << transitionDistance << " resolution: " << resolution[0] << " " << resolution[1] << " " << resolution[2] << std::endl;
			
// 			compute_scalar_field(*lockedCon);
// 			marching_cubes();
// 			estimate_metrics(*lockedCon);

// 			auto endTime = std::chrono::steady_clock::now();

// 			// Calculate the duration in milliseconds
// 			auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

// 			std::ostringstream oss;
// 			oss << std::fixed << std::setprecision(3) // Set precision to 3 decimal places
// 				<< duration_ms.count() / 1000.0   // Convert ms to seconds
// 				<< " seconds!";

// 			logger->log(LogPriority::SUCCESS, "Updated Scaffold Successfully in " + oss.str());

// 			// reset
// 			updateScaffold = false;
// 		}
// 	}
// };

void GeneratorLewiner::render_properties(
	bool& updateScaffold,
	 GenerationTask* task,
	 std::vector<std::shared_ptr<AnisotropySource>>& globalSources
)
{
    std::shared_ptr<IContainer> lockedCon = container.lock();
    std::shared_ptr<InterfaceSeedGenerator> lockedGen = generator.lock();

    const bool mineBusy = task && task->is_running_for(this);

	ImGui::SeparatorText("Container & Generator");
    if (lockedCon) ImGui::Text("Container: %s", lockedCon->name.c_str());
    if (lockedGen){
		ImGui::Text("Generator: %s", lockedGen->name.c_str());
		if (lockedGen->get_type() == ObjectType::RandomGeneratorType) {
			ImGui::SameLine();
			ImGui::Text("Random Seeds %d", lockedGen->get_seeds().size());
		}
		else if (lockedGen->get_type() == ObjectType::PoissonGeneratorType) {
			Poisson3D* dummy = static_cast<Poisson3D*>(lockedGen.get());
			if (dummy->is_uniform()) {
				ImGui::SameLine();
				ImGui::Text("Radius (Uniform) %.4f", dummy->get_min_radius());
            } 
			else {
				ImGui::Text(
					"Rmin %.4f Rmax %4.f",
					dummy->get_min_radius(), dummy->get_max_radius());
            }			
		}
	} 	

	ImGui::SeparatorText("Settings");

    // ---- editable parameters: locked while this scaffold is regenerating ----
    ImGui::BeginDisabled(mineBusy);

	ImGui::BeginTabBar("Items");

	main_properties();

	thickness_properties();

	anisotropy_properties(globalSources);

	smooth_properties();

	ImGui::EndTabBar();


    ImGui::EndDisabled();

    // ---- metrics / progress ----
    ImGui::SeparatorText("Metrics");
    if (mineBusy) {
        ImGui::ProgressBar(task->get_progress(), ImVec2(-1, 0));
        ImGui::TextDisabled("Regenerating...");
    } else {
        render_metrics();
    }

    // ---- launch (one-shot trigger from the caller's Update button) ----
    if (updateScaffold) {
        updateScaffold = false;  // consume the trigger immediately

		logger->log(LogPriority::INFO,
			"update: calT=" + std::to_string(calibrateThickness) +
			" calP=" + std::to_string(calibratePorosity) +
			" targetTh=" + std::to_string(targetThickness) +
			" targetPor%=" + std::to_string(targetPorosity) +
			" con=" + std::to_string((bool)lockedCon) +
			" gen=" + std::to_string((bool)lockedGen) +
			" voxel=" + std::to_string(voxelSize) +
			" measVoxel=" + std::to_string(measurementVoxelSize));

        if (mineBusy) {
            // already regenerating this scaffold; ignore
        }
        else if (isLoadedFromFile) {
            logger->log(LogPriority::WARNING, "Cannot update a static mesh loaded from a file.");
        }
        else if (lockedCon && lockedGen && task && !task->get_running()) {

            std::vector<Vec3> newSeeds = lockedGen->get_seeds();
            Bounds bds = lockedCon->compute_bounds();
            std::array<float, 6> newBounds = {
                bds.xMin, bds.xMax, bds.yMin, bds.yMax, bds.zMin, bds.zMax
            };

            blockDims = {
                static_cast<int>(std::ceil((bds.xMax - bds.xMin) / voxelSize)) + 1,
                static_cast<int>(std::ceil((bds.yMax - bds.yMin) / voxelSize)) + 1,
                static_cast<int>(std::ceil((bds.zMax - bds.zMin) / voxelSize)) + 1
            };

            if (selectedThicknessOption == 1) {
                switch (selectedFunc) {
                    case 0: thicknessFunction = std::make_shared<LinearFunction>(transitionDistance);    break;
                    case 1: thicknessFunction = std::make_shared<QuadraticFunction>(transitionDistance); break;
                    case 2: thicknessFunction = std::make_shared<SmoothStep>(); break;
                    case 3: thicknessFunction = std::make_shared<ConstantRadiusFunction>();              break;
                }
                switch (selectedDist) {
                    case 0: thicknessSDF = std::make_shared<PlaneSDF>(distancePlaneCenter, distancePlaneNormal); break;
                    case 1: thicknessSDF = std::make_shared<PointSDF>(distancePoint);                            break;
                    case 2: thicknessSDF = lockedCon->get_distance_estimator();                                  break;
                }
            } else {
                thicknessSDF.reset();
                thicknessFunction.reset();
            }

            set_thickness_functions(thicknessSDF, thicknessFunction,
                                    startThickness, endThickness, transitionDistance);
            tortuosityPathModel.reset();
            set_bounds(newBounds);
            set_seeds(newSeeds);

            // keep the container alive for the whole job; capture `this` for the compute
            auto conShared = lockedCon;
            GeneratorLewiner* self = this;
            // wire cooperative cancellation to the running task's cancel flag
            self->cancelRequested = [task]{ return task->is_cancel_requested(); };
            start_time = std::chrono::steady_clock::now();

			const bool anyCalibration =
				calibrateThickness || calibratePorosity /* || calibrateSmi */ || calibrateDA;

			if (anyCalibration) {
				auto cStep = measurementVoxelSize; // metrics at the measurement voxel
				task->start(
					[this, conShared, cStep, t = task]() {
					t->set_progress(0.00f);

					// Snapshot the knobs the calibration may move, so a cancel truly
					// changes nothing (the popup also skips the GPU upload on cancel).
					const float snapIso = isoLevel, snapThr = threshold, snapSpread = spread;
					const float snapS0 = startThickness, snapS1 = endThickness;
					const float snapSx = stretchX, snapSy = stretchY, snapSz = stretchZ;

					this->build_calibration_stages(cStep);
					bool ok = !t->is_cancel_requested()
						&& this->solve_calibration(this->stages, 0, *conShared);

					if (t->is_cancel_requested()) {
						// restore the pre-calibration state; nothing is applied
						isoLevel = snapIso; threshold = snapThr; spread = snapSpread;
						startThickness = snapS0; endThickness = snapS1;
						set_stretch(snapSx, snapSy, snapSz);
					}
					else if (ok && this->marching_cubes()) {
						t->set_progress(0.90f);
						// full-resolution DA (the loop used a reduced ray budget)
						if (this->calibrateDA)
							this->estimate_anisotropy(cStep, 2000, 10000, this->targetFormulaIdx);
						this->estimate_metrics(*conShared);
						// stamp calibrated metric versions to this final mesh so the
						// inner-loop targets (Tb.Th, SMI) no longer read stale.
						this->stamp_calibrated_versions();
					}
					t->set_progress(1.00f);
				}, this);
			}
			else {
				task->start(
					[this, conShared, t = task]() {
					t->set_progress(0.00f);
					if (!t->is_cancel_requested() &&
						this->compute_scalar_field(*conShared)) {
						t->set_progress(0.50f);
						if (!t->is_cancel_requested() && 
							this->marching_cubes()) {
							t->set_progress(0.90f);
							if (!t->is_cancel_requested()) 
							this->estimate_metrics(*conShared);
						}
					}
					t->set_progress(1.00f);
				}, this);
			}
                                 
			ImGui::OpenPopup("Updating Scaffold..."); // open the progress bar
		}
    }

	// -- progress bar --
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
    ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

   if (ImGui::BeginPopupModal("Updating Scaffold...", NULL, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove)) {
        
        ImGui::Text("Regenerating '%s' in the background.", name.empty() ? "Scaffold" : name.c_str());
        ImGui::Dummy(ImVec2(0.0f, 5.0f)); // Small spacer
        
        // Fixed width for the progress bar looks better in a floating window
        ImGui::ProgressBar(task->get_progress(), ImVec2(300.0f, 0.0f));

        // ---- cancel: request cooperative cancellation of the worker ----
        ImGui::Dummy(ImVec2(0.0f, 5.0f));
        const bool cancelling = task && task->is_cancel_requested();
        ImGui::BeginDisabled(cancelling);
        if (ImGui::Button(cancelling ? "Cancelling..." : "Cancel", ImVec2(300.0f, 0.0f))) {
            if (task) task->request_cancel();
        }
        ImGui::EndDisabled();

        // ---- completion (main thread: GL upload + log) ----
        if (task && task->poll(this)) {
            update_render();                               // GPU upload, GL thread

            if (task->is_cancel_requested()) {
                logger->log(LogPriority::WARNING,
                    "Scaffold update cancelled; showing the last valid state.");
            }
            else {
                auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::steady_clock::now() - start_time);
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(3) << dur.count() / 1000.0 << " seconds!";
                logger->log(LogPriority::SUCCESS, "Updated Scaffold Successfully in " + oss.str());
            }

            ImGui::CloseCurrentPopup();
        }

        ImGui::EndPopup();
    }
}

void GeneratorLewiner::main_properties(){

	if(ImGui::BeginTabItem("Main")){
		
		ImGui::ColorEdit4("Appearance", (float*)&color);  
		ImGui::InputFloat("Voxel Size", &voxelSize);
		ImGui::InputFloat("Measurement Voxel Size", &measurementVoxelSize);
		
		ImGui::SeparatorText("Porosity Calibration");

		ImGui::Checkbox("Calibrate Porosity", &calibratePorosity);
		ImGui::InputFloat("Target Porosity", &targetPorosity);

		ImGui::EndTabItem();
	}

};

void GeneratorLewiner::thickness_properties(){

	if (ImGui::BeginTabItem("Thickness")){

		ImGui::Checkbox("Calibrate Thickness", &calibrateThickness);

		ImGui::Separator();

		ImGui::RadioButton("Apply Uniform Thickness", &selectedThicknessOption, 0);
		ImGui::RadioButton("Apply Varied Thickness", &selectedThicknessOption, 1);

		if (selectedThicknessOption == 0) {
			ImGui::InputFloat("Iso Surface Level", &isoLevel, 0.001f, 1.0f);
			ImGui::InputFloat(
				"Target Thickness", &targetThickness, 0.001f, 1.0f);
			// uniform: the Thickness value IS the calibration target (Tb.Th)
			if (calibrateThickness)
				ImGui::Text("Calibration target (Tb.Th): %.3f mm", targetThickness);
		}
		else {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Start Thickness", &startThickness);
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("End Thickness", &endThickness);
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Transition Distance", &transitionDistance);
			// varied: the calibration target is the MEAN of the (scaled) range
			if (calibrateThickness)
				ImGui::Text("Calibration target (mean Tb.Th): %.3f mm",
					(startThickness + endThickness) * 0.5f);

			ImGui::SeparatorText("Select Distance Function");
			ImGui::RadioButton("Distance From Plane", &selectedDist, 0);
			if (selectedDist == 0) {
				ImGui::SetNextItemWidth(200);
				ImGui::InputFloat3("Normal", distancePlaneNormal);
				ImGui::SetNextItemWidth(200);
				ImGui::InputFloat3("Center", distancePlaneCenter);
			};
			ImGui::RadioButton("Distance From Point", &selectedDist, 1);
			if (selectedDist == 1) {
				ImGui::SetNextItemWidth(200);
				ImGui::InputFloat3("Point", distancePoint);
			}
			ImGui::RadioButton("Distance From Container", &selectedDist, 2);
			
			ImGui::SeparatorText("Select Radius Function");
			ImGui::RadioButton("Linear", &selectedFunc, 0);
			ImGui::RadioButton("Quadratic", &selectedFunc, 1);
			ImGui::RadioButton("Constant", &selectedFunc, 2);
			ImGui::RadioButton("Random", &selectedFunc, 3);
		}

		ImGui::SliderFloat("Openess", &threshold, 0.0f, 1.0f);
		ImGui::SliderFloat("Spread", &spread, 0.0f, 1.0f);
		ImGui::EndTabItem();
	}
}

void GeneratorLewiner::anisotropy_properties(
	std::vector<std::shared_ptr<AnisotropySource>>& globalSources
){

	static int selectedIdx = -1;
	if(ImGui::BeginTabItem("Anisotropy")){

		ImGui::Checkbox("Calibrate", &calibrateDA);
		ImGui::InputFloat("Calibration Target", &targetDa);

		ImGui::Separator();
		
		ImGui::SeparatorText("Global Background");
		ImGui::InputFloat("Stretch X", &stretchX, 0.01f, 5.0f, "%.3f");
		ImGui::InputFloat("Stretch Y", &stretchY, 0.01f, 5.0f, "%.3f");
		ImGui::InputFloat("Stretch Z", &stretchZ, 0.01f, 5.0f, "%.3f");
		ImGui::InputFloat("Angle", &anisotropyAngle, 0.01f, 10.0f, "%.4f");
		ImGui::InputFloat("Background Weight", &backgroundWeight, 0.01f, 0.5f, "%.3f");
		ImGui::InputFloat3("Direction (Local X)", anisotropyVec, "%.3f");
		
		ImGui::SameLine();
		if (ImGui::Button("Normalize")) {
			float len = std::sqrt(
				anisotropyVec.x*anisotropyVec.x + 
				anisotropyVec.y*anisotropyVec.y + 
				anisotropyVec.z*anisotropyVec.z);
			if (len > 1e-6f) {
				anisotropyVec.x /= len;
				anisotropyVec.y /= len;
				anisotropyVec.z /= len;
			}
		}

		for (int i{ 0 }; i < (int)anisotropySources.size(); i++) {
			std::string label = anisotropySources[i]->name.empty()
				? "Anisotropy Source " + std::to_string(i + 1)
				: anisotropySources[i]->name;
			if (ImGui::TreeNode(label.c_str())) {
				anisotropySources[i]->render_properties();
				selectedIdx = i;
				ImGui::TreePop();
			}
		}

		if (ImGui::Button("Add Anisotropy Source")){

			std::shared_ptr<AnisotropySource> source = std::make_shared<AnisotropySource>();
			source->stretch = Vec3(1.0f, 1.0f, 1.0f);
			source->sigma = 5.0f;
			source->update_metric();
			source->update_model();
			anisotropySources.push_back(source);
			globalSources.push_back(source);
		}

		ImGui::SameLine();

		if(ImGui::Button("Delete Selected")){
			anisotropySources.erase(
				anisotropySources.begin() + selectedIdx
			);
			globalSources.erase(
				globalSources.begin() + selectedIdx
			);
		}

		ImGui::EndTabItem();
	}
	else{
		selectedIdx = -1;
	}
}

void GeneratorLewiner::smooth_properties(){

	if (ImGui::BeginTabItem("Smoothness")){

		ImGui::SeparatorText("Mesh smoothing (Taubin)");
		ImGui::InputInt("Iterations", &iter, 1);

		ImGui::InputFloat("Lambda", &lambda, 0.01f, 10.0f);

		ImGui::InputFloat("Mu", &mu, 0.01f, 10.0f);

		ImGui::SeparatorText("Junction smoothing");
		ImGui::Checkbox("Smooth rod-plate junctions", &smoothJunctions);
		ImGui::SetItemTooltip("Fillet/fatten the trabecular nodes (more organic, "
			"bone-like). Tb.Th inflation is absorbed by thickness calibration.");
		if (smoothJunctions) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Junction k", &smoothK, 0.001f, 0.1f, "%.4f");
		}

		ImGui::SeparatorText("Edge rounding");
		ImGui::Checkbox("Round cell edges", &roundEdges);
		ImGui::SetItemTooltip("Soften the raw Voronoi distance field so the cell "
			"EDGES and VERTICES become round rather than polygonal (rounder "
			"pores). Unlike junction smoothing it acts before the openness blend, "
			"so it also rounds bare (open-lattice) struts. Tb.Th inflation is "
			"absorbed by thickness calibration. Requires re-generation (recache).");
		if (roundEdges) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Edge k", &edgeK, 0.001f, 0.1f, "%.4f");
		}

		ImGui::SeparatorText("Boundary frame");
		ImGui::Checkbox("Frame boundary cells", &frameBoundary);
		ImGui::SetItemTooltip("Close the outermost Voronoi cells into full walls, "
			"forming a connected honeycomb rim around the scaffold. Ties the cut "
			"boundary struts together for gripping/loading test specimens.");
		if (frameBoundary) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Frame depth (x thickness)", &frameDepth, 0.1f, 0.5f, "%.2f");
			ImGui::SetItemTooltip("Rim depth as a multiple of the wall thickness. "
				"Keep it small (~0.5) so the boundary stays a thin honeycomb rim; "
				"large values turn the boundary shell plate-like.");
		}
		ImGui::Checkbox("Frame container edges", &frameContainerEdges);
		ImGui::SetItemTooltip("Add solid beams along the container's own edges "
			"(a box's 12 edges, a cylinder's two rims) for a rigid outer cage. "
			"Box and cylinder containers only.");
		if (frameContainerEdges) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Edge beam (mm)", &frameBeam, 0.05f, 0.5f, "%.3f");
		}

		ImGui::EndTabItem();
	}

};

void GeneratorLewiner::estimate_metrics(const IContainer& container) {

	surfaceArea = 0.0f;
	volume = 0.0f;
	float signedVolume = 0.0f;

	for (const auto& tri : meshTriangles) {

		// grab the three vertices
		const LVertex& v1 = meshVertices[tri.v1];
		const LVertex& v2 = meshVertices[tri.v2];
		const LVertex& v3 = meshVertices[tri.v3];
	
		// estimate the surface area
		Vec3 edge1 = Vec3(v2.x, v2.y, v2.z) - Vec3(v1.x, v1.y, v1.z);
		Vec3 edge2 = Vec3(v3.x, v3.y, v3.z) - Vec3(v1.x, v1.y, v1.z);
	
		surfaceArea += 0.5 * edge1.cross(edge2).norm();

		// get the volume using the signed tetrahedron 
		float vol = Vec3(v1.x, v1.y, v1.z).dot(Vec3(v2.x, v2.y, v2.z).cross(Vec3(v3.x, v3.y, v3.z)));
		signedVolume += vol;
	}

	// we need to divide the volume by 1/6
	volume = std::abs(signedVolume) / 6.0f;

	// mesh-based (legacy) porosity: mesh volume vs analytic container volume.
	// Overestimates porosity systematically - the mesh can never reach the
	// container surface. The primary 'porosity' is the voxel-based estimate
	// set by compute_scalar_field().
	float domainVolume = container.get_volume();
	porosityMesh = (1 - volume / domainVolume);

	// scaffolds loaded from STL have no scalar field, so the mesh estimate
	// is the only porosity available for them
	if (scalarField.empty()) {
		porosity = porosityMesh;
	}

	// Two distinct normalizations, both standard (Parfitt/ASBMR). Do not mix
	// them up when comparing against uCT literature - they differ by BV/TV,
	// roughly an order of magnitude:
	//   BS/BV  surface per BONE volume. ~2/Tb.Th (plates) .. ~4/Tb.Th (rods).
	//   BS/TV  surface per TOTAL sample volume = (BS/BV) * (BV/TV).
	surfaceToVolume = surfaceArea / volume;
	surfaceToTotalVolume = (domainVolume > 0.0f) ? (surfaceArea / domainVolume) : 0.0f;
};

float GeneratorLewiner::estimate_smi(float dilation) {

	// Structure Model Index (Hildebrand & Ruegsegger 1997):
	//     SMI = 6 * (BS' * BV) / BS^2
	// BS' = dBS/dr is the rate at which the surface area grows when the surface
	// is dilated. The intuition: dilating a PLATE barely changes its area (the
	// two faces just move apart), so BS' ~ 0 -> SMI ~ 0. Dilating a ROD grows
	// its lateral area in proportion to the radius, so BS' > 0 -> SMI ~ 3.
	// A sphere gives ~4.
	//
	// BS' is estimated by a forward difference: push every vertex a small step
	// along its (already normalized) vertex normal and re-measure the area.
	if (meshTriangles.empty() || meshVertices.empty()) {
		if (logger) logger->log(LogPriority::ERROR, "Cannot estimate SMI: no mesh.");
		smi = 0.0f;
		return smi;
	}

	// A scale-relative default: a fraction of a voxel. Too large and the
	// dilation is no longer a derivative; too small and it drowns in float
	// noise. Must be well below the strut thickness.
	if (dilation <= 0.0f) {
		dilation = 0.05f * std::max(stepX, std::max(stepY, stepZ));
	}

	auto surface_area_of = [&](float offset) -> double {
		double area = 0.0;
		#pragma omp parallel for reduction(+:area)
		for (long long t = 0; t < static_cast<long long>(meshTriangles.size()); ++t) {
			const LTriangle& tri = meshTriangles[t];
			const LVertex& a = meshVertices[tri.v1];
			const LVertex& b = meshVertices[tri.v2];
			const LVertex& c = meshVertices[tri.v3];

			// displace each vertex along its own normal
			Vec3 p1(a.x + a.nx * offset, a.y + a.ny * offset, a.z + a.nz * offset);
			Vec3 p2(b.x + b.nx * offset, b.y + b.ny * offset, b.z + b.nz * offset);
			Vec3 p3(c.x + c.nx * offset, c.y + c.ny * offset, c.z + c.nz * offset);

			area += 0.5 * (p2 - p1).cross(p3 - p1).norm();
		}
		return area;
	};

	const double bs = surface_area_of(0.0f);
	const double bsDilated = surface_area_of(dilation);
	const double bsPrime = (bsDilated - bs) / static_cast<double>(dilation);

	// volume/surfaceArea are filled by estimate_metrics(); recompute the volume
	// here only if it has not been measured yet.
	double bv = static_cast<double>(volume);
	if (bv <= 0.0) {
		double signedVolume = 0.0;
		for (const auto& tri : meshTriangles) {
			const LVertex& v1 = meshVertices[tri.v1];
			const LVertex& v2 = meshVertices[tri.v2];
			const LVertex& v3 = meshVertices[tri.v3];
			signedVolume += Vec3(v1.x, v1.y, v1.z).dot(
				Vec3(v2.x, v2.y, v2.z).cross(Vec3(v3.x, v3.y, v3.z)));
		}
		bv = std::abs(signedVolume) / 6.0;
	}

	smi = (bs > 1e-12) ? static_cast<float>(6.0 * bsPrime * bv / (bs * bs)) : 0.0f;

	smiVersion = meshVersion;

	if (logger) {
		std::ostringstream oss;
		oss << "Estimated SMI = " << smi << " (0=plate, 3=rod)";
		logger->log(LogPriority::SUCCESS, oss.str());
	}

	return smi;
};

void GeneratorLewiner::render_metrics() {

	// define a lambda function for checking version
	auto draw_metric_row = [&](const char* label, uint32_t metricVersion, auto value_func) {
		bool needsUpdate = metricVersion != meshVersion;
		ImVec4 color = needsUpdate ? ImVec4(1.0f, 0.6f, 0.0f, 1.0f) : ImVec4(1.0f, 1.0f, 1.0f, 1.0f);
		ImGui::TableNextRow();
		ImGui::TableNextColumn();
		ImGui::TextColored(color, "%s", label);

		ImGui::TableNextColumn();
		ImGui::PushStyleColor(ImGuiCol_Text, color);
		value_func();
		ImGui::PopStyleColor();
	};

	ImGuiTableFlags flags = ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_BordersInnerV;
	// create a table
	if (ImGui::BeginTable("", 2, flags = flags)) {
		ImGui::TableNextRow();
		ImGui::TableNextColumn(); ImGui::Text("Volume (mm^3)"); // Note: Volume is mm^3
		ImGui::TableNextColumn(); ImGui::Text("%.4f", volume);

		ImGui::TableNextRow();
		ImGui::TableNextColumn(); ImGui::Text("Total Surface (mm^2)");
		ImGui::TableNextColumn(); ImGui::Text("%.4f", surfaceArea);

		ImGui::TableNextRow();
		ImGui::TableNextColumn(); ImGui::Text("BS/BV - surface / bone vol (1/mm)");
		ImGui::TableNextColumn(); ImGui::Text("%.4f", surfaceToVolume);
		ImGui::SetItemTooltip("Bone surface per BONE volume. ~2/Tb.Th for plates, ~4/Tb.Th for rods.");

		ImGui::TableNextRow();
		ImGui::TableNextColumn(); ImGui::Text("BS/TV - surface / total vol (1/mm)");
		ImGui::TableNextColumn(); ImGui::Text("%.4f", surfaceToTotalVolume);
		ImGui::SetItemTooltip("Bone surface per TOTAL sample volume = (BS/BV) * (BV/TV).\nThis is the 'BS/TV' reported in the uCT literature.");

		ImGui::TableNextRow();
		ImGui::TableNextColumn(); ImGui::Text("Porosity (%%)");
		ImGui::TableNextColumn(); ImGui::Text("%.4f", porosity * 100.0f);

		ImGui::TableNextRow();
		ImGui::TableNextColumn(); ImGui::Text("Porosity - mesh, legacy (%%)");
		ImGui::TableNextColumn(); ImGui::Text("%.4f", porosityMesh * 100.0f);

		draw_metric_row("Local Thickness (mm)", thicknessVersion, [&]() {
			ImGui::Text("%.4f std: %.4f", localThickness, localThicknessStd);
		});

		draw_metric_row("Local Separation (mm)", separationVersion, [&]() {
			ImGui::Text("%.4f std: %.4f", localSeparation, localSeparationStd);
		});

		draw_metric_row("Trabecular Number (1/mm)", trabecularNrVersion, [&]() {
			ImGui::Text("%.4f", trabecularNr);
		});

		draw_metric_row("Connectivity Density (1/mm^3)", connectivityVersion, [&]() {
			ImGui::Text("%.4f", connectivityDensity);
		});

		draw_metric_row("Tortuosity", tortuosityVersion, [&]() {
			ImGui::Text("%.4f", tortuosity);
		});

		draw_metric_row("Degree of Anisotropy", anisotropyVersion, [&]() {
			ImGui::Text("%.4f", anisotropyDegree);
		});

		draw_metric_row("Structure Model Index", smiVersion, [&]() {
			ImGui::Text("%.4f", smi);
		});

		ImGui::EndTable();
	};
};

void GeneratorLewiner::set_thickness_functions(
	std::shared_ptr<const SDF> sdf,
	std::shared_ptr<const RadiusFunction> radFunc,
	float tMin, float tMax, float distance
) {
	thicknessFunction = radFunc;
	thicknessSDF = sdf;
	startThickness = tMin;
	endThickness = tMax;
	transitionDistance = distance;
};

void GeneratorLewiner::set_options_from_factory(
	int distOption, 
	int distFunc, int thicknessOption, 
	float voxSize, 
	float measureVoxel) {
	selectedDist = distOption;
	selectedFunc = distFunc;
	selectedThicknessOption = thicknessOption;
	voxelSize = voxSize;
	measurementVoxelSize = measureVoxel;
};

void GeneratorLewiner::set_distance_plane_options(Vec3 center, Vec3 normal) {
	distancePlaneCenter = center;
	distancePlaneNormal = normal;
};

void GeneratorLewiner::set_distance_point_options(Vec3 point) {
	distancePoint = point;
};

// These change the per-voxel cache's inputs (grid / seeds / warp), so the cache
// is invalidated; the next calibrate_* / compute_scalar_field rebuilds it.
void GeneratorLewiner::set_resolution(const std::array<int, 3>& newResolution) {
	blockDims = newResolution;
	cachedVoxels.clear();
};

void GeneratorLewiner::set_bounds(const std::array<float, 6>& newBounds) {
	bounds = newBounds;
	stepX = (bounds[1] - bounds[0]) / (blockDims[0] - 1);
	stepY = (bounds[3] - bounds[2]) / (blockDims[1] - 1);
	stepZ = (bounds[5] - bounds[4]) / (blockDims[2] - 1);
	cachedVoxels.clear();
};

void GeneratorLewiner::set_seeds(const std::vector<Vec3>& newSeeds) {
	seeds = newSeeds;
	cachedVoxels.clear();
};

void GeneratorLewiner::set_stretch(float newStretchX, float newStretchY, float newStretchZ) {
	this->stretchX = newStretchX;
	this->stretchY = newStretchY;
	this->stretchZ = newStretchZ;
	cachedVoxels.clear();
};

void GeneratorLewiner::set_thickness(float newThickness) {
	isoLevel = newThickness;
};

std::array<float, 6> GeneratorLewiner::get_bounds() const {
	// Tight bounds around the generated mesh - what every external caller
	// (camera framing, cut plane, thickness analysis) expects. The internal
	// 'bounds' member keeps defining the generation grid and is returned
	// only while no mesh exists yet.
	if (!meshVertices.empty()) {
		return { aabb.pMin.x, aabb.pMax.x,
				 aabb.pMin.y, aabb.pMax.y,
				 aabb.pMin.z, aabb.pMax.z };
	}
	return bounds;
};

Aabb GeneratorLewiner::get_aabb() const { return aabb; };

void GeneratorLewiner::estimate_local_thickness(
	float voxelSize, std::array<float, 6>& blockBounds, bool separation, bool supress) {

	// Clip the requested window to the GENERATION grid - the domain the scalar
	// field is defined on - NOT the mesh AABB. This makes the measurement extent
	// identical whether or not a mesh exists (calibration runs without one), and
	// immune to a stale mesh AABB from a previous build. Passing the generation
	// bounds therefore measures the whole container interior (out-of-container
	// voxels are still excluded by the is_inside gate below), while an explicit
	// ROI (a sub-box) is preserved. get_image_field already treats anything past
	// the grid as air, so the clip only guards against a degenerate window.
	blockBounds[0] = std::max(blockBounds[0], bounds[0]);
	blockBounds[1] = std::min(blockBounds[1], bounds[1]);
	blockBounds[2] = std::max(blockBounds[2], bounds[2]);
	blockBounds[3] = std::min(blockBounds[3], bounds[3]);
	blockBounds[4] = std::max(blockBounds[4], bounds[4]);
	blockBounds[5] = std::min(blockBounds[5], bounds[5]);

	if (blockBounds[0] >= blockBounds[1] ||
		blockBounds[2] >= blockBounds[3] ||
		blockBounds[4] >= blockBounds[5]) {
		std::cerr << "Error: The requested bounds do not overlap the generation grid." << std::endl;
		return; // Abort early to prevent voxel creation crashes
	}

	std::vector<uint8_t> field = get_image_field(voxelSize, blockBounds, separation);

	// estimate new block dimensions
	int nx = static_cast<int>(std::ceil((blockBounds[1] - blockBounds[0]) / voxelSize));
	int ny = static_cast<int>(std::ceil((blockBounds[3] - blockBounds[2]) / voxelSize));
	int nz = static_cast<int>(std::ceil((blockBounds[5] - blockBounds[4]) / voxelSize));

	size_t totalVoxels = field.size();

	// initialize all values to infinity
	std::vector<float> squaredDistanceField;
	if (separation) {
		float maxLinearDist = (float)(nx + ny + nz);
		squaredDistanceField.resize(totalVoxels, maxLinearDist);
	}
	else {
		squaredDistanceField.resize(totalVoxels, std::numeric_limits<float>::max());
	}

	auto get_idx = [&](int x, int y, int z) -> size_t {
		return static_cast<size_t>(x) +
			static_cast<size_t>(y) * static_cast<size_t>(nx) +
			static_cast<size_t>(z) * static_cast<size_t>(nx) * static_cast<size_t>(ny);
	};

	// initialize to zero the empty voxels
	auto con = container.lock();

	// initialize to zero the empty voxels AND voxels outside the container
	#pragma omp parallel for collapse(3)
	for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx; x++) {
				size_t idx = get_idx(x, y, z);

				bool isInsideROI = true;
				if (con) {
					// Calculate the physical center of the voxel
					float px = blockBounds[0] + (x + 0.5f) * voxelSize;
					float py = blockBounds[2] + (y + 0.5f) * voxelSize;
					float pz = blockBounds[4] + (z + 0.5f) * voxelSize;
					isInsideROI = con->is_inside(Vec3(px, py, pz));
				}

				// If it's outside the container, it acts as a hard boundary (0.0f).
				// Otherwise, check the image field as usual.
				if (!isInsideROI || field[idx] == 0) {
					squaredDistanceField[idx] = 0.0f;
				}
			}
		}
	}

	// since we have 3d data we need three phases to estimate the distance field one along each dimension, we can use a simple 1D distance transform along each dimension
	// 1st pass: forward scan along x
	#pragma omp parallel for
	for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
			// forward scan
			for (int x = 1; x < nx; x++) {
				size_t idx = get_idx(x, y, z);
				squaredDistanceField[idx] = std::min(
					squaredDistanceField[idx], squaredDistanceField[get_idx(x - 1, y, z)] + 1.0f);
			}
			// backward scan
			for (int x = nx - 2; x >= 0; x--) {
				size_t idx = get_idx(x, y, z);
				squaredDistanceField[idx] = std::min(
					squaredDistanceField[idx], squaredDistanceField[get_idx(x + 1, y, z)] + 1.0f);
			}

			// use the squared distance for next pass 
			for (int x = 0; x < nx; x++) {
				//int idx = find_vertex_index(x, y, z);
				size_t idx = get_idx(x, y, z);
				squaredDistanceField[idx] = squaredDistanceField[idx] * squaredDistanceField[idx];
			}
		}
	}
	// pass along y
	#pragma omp parallel for
	for (int z{ 0 }; z < nz; z++) {
		for (int x{ 0 }; x < nx; x++) {
			std::vector<float> g(ny);
			for (int y = 0; y < ny; y++) {
				g[y] = squaredDistanceField[get_idx(x, y, z)];
			}
			// we need to find for each column the minimum distance for the current x, z coordinate
			// create the list of segments, in these segmentes we will check the intersection of the parabolas
			// since the euclidean distance is a parabola
			std::vector<int> s(ny);
			std::vector<int> t(ny);
			int q = 0;
			s[0] = 0;
			t[0] = -1e10f;
			for (int u = 1; u < ny; u++) {
				while (q >= 0) {
					int w = s[q];
					float f_u_w = ((g[u] + u * u) - (g[w] + w * w)) / (2.0f * (u - w));
					if (f_u_w <= t[q]) {
						q--;
					}
					else {
						q++;
						s[q] = u;
						t[q] = f_u_w;
						break;
					}
				}

				if (q < 0) {
					q = 0;
					s[0] = u;
					t[0] = -1e10f;
				}
			}
			// second pass 
			for (int u = ny - 1; u >= 0; u--) {
				while (u < t[q]) q--;
				float dy = (float)(u - s[q]);
				squaredDistanceField[get_idx(x, u, z)] = g[s[q]] + dy * dy;
			}
		}
	}

	// finally pass along z it is the same as the y pass but we need to iterate along z and keep x, y fixed
	#pragma omp parallel for
	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++) {

			std::vector<float> g(nz);

			for (int z = 0; z < nz; z++) g[z] = squaredDistanceField[get_idx(x, y, z)];
			std::vector<int> s(nz);
			std::vector<float> t(nz);

			int q = 0; s[0] = 0; t[0] = -1e10f;
			for (int u = 1; u < nz; u++) {
				while (q >= 0) {
					int w = s[q];
					float f_u_w = ((g[u] + u * u) - (g[w] + w * w)) / (2.0f * (u - w));
					if (f_u_w <= t[q]) q--;
					else {
						q++;
						s[q] = u;
						t[q] = f_u_w;
						break;
					}
				}
				if (q < 0) { q = 0; s[0] = u; t[0] = -1e10f; }
			}
			for (int u = nz - 1; u >= 0; u--) {
				while (u < t[q]) q--;
				float dz = (float)(u - s[q]);
				squaredDistanceField[get_idx(x, y, u)] = g[s[q]] + dz * dz;
			}
		}
	}
	// now that we have our distance field we can find redundant voxels, these are not local maxima in their neighborhood, we can remove them by setting their distance to zero, this will give us a skeleton of the solid part that represents the local thickness
	std::vector<size_t> redundantVoxels;

	#pragma omp parallel
	{
		std::vector<size_t> localRedundant;
		#pragma omp for
		for (int z = 1; z < nz - 1; z++) {
			for (int y = 1; y < ny - 1; y++) {
				for (int x = 1; x < nx - 1; x++) {
					size_t idx = get_idx(x, y, z);

					float r = std::sqrt(squaredDistanceField[idx]);
					//if (scalarField[idx] >= isoLevel) continue; // only consider solid voxels
					if (field[idx] == 0) continue; // only consider solid voxels

					bool isLocalMax = true;
					// check 26 neighbors we already have the indices in the const std::array<Neighbor, 26> neighbors> , we can just iterate through them and check if any of the neighbors has a higher distance value than the current voxel, if it does, then the current voxel is not a local maximum and we can mark it as redundant
					for (const auto& nb : neighbors) {
						size_t nidx = get_idx(x + nb.dx, y + nb.dy, z + nb.dz);
						//if (scalarField[nidx] >= isoLevel) continue; // only consider solid voxels
						if (field[nidx] == 0) continue; // only consider solid voxels

						float r_nb = std::sqrt(squaredDistanceField[nidx]);
						float dist_centers = std::sqrt(nb.dx * nb.dx + nb.dy * nb.dy + nb.dz * nb.dz);
						// Paper Inclusion Test: If this sphere is inside neighbor's sphere
						if (r + dist_centers <= r_nb + 0.001f) { // epsilon for float stability
							isLocalMax = false;
							break;
						}
					}
					// if it is a local maximum, we keep it, otherwise we mark it as redundant
					if (isLocalMax) localRedundant.push_back(idx);
				}
			}
		}

	#pragma omp critical
		redundantVoxels.insert(redundantVoxels.end(), localRedundant.begin(), localRedundant.end());

	}

	// now estimate the local thickness at any voxel as the maximum diameter of all the spheres that can fit inside the solid part at that point, which is just twice the distance value at that voxel in the distance field
	std::vector<omp_lock_t> z_locks(nz);
	for (int z = 0; z < nz; z++) {
		omp_init_lock(&z_locks[z]);
	}

	std::vector<float> thicknessMap(totalVoxels, 0.0f);

	#pragma omp parallel for schedule(dynamic)
	for (long long i = 0; i < static_cast<long long>(redundantVoxels.size()); i++) {

		size_t ridgeIdx = redundantVoxels[i];

		float radius = std::sqrt(squaredDistanceField[ridgeIdx]);
		float diameter = 2.0f * radius;
		float rSq = radius * radius;

		// decode ridgeIdx to grid coords, matching get_idx()'s layout
		// idx = x + y*nx + z*nx*ny  (x fastest)
		int rx = static_cast<int>(ridgeIdx % static_cast<size_t>(nx));
		int ry = static_cast<int>((ridgeIdx / static_cast<size_t>(nx)) % static_cast<size_t>(ny));
		int rz = static_cast<int>(ridgeIdx / (static_cast<size_t>(nx) * static_cast<size_t>(ny)));

		// Define a bounding box for the sphere to limit the search
		int R = (int)std::ceil(radius);

		int zMin = std::max(0, rz - R);
		int zMax = std::min(nz - 1, rz + R);

		for (int z = zMin; z <= zMax; z++) {

			float dz = (float)(z - rz);
			float dzSq = dz * dz;

			int yMin = std::max(0, ry - R);
			int yMax = std::min(ny - 1, ry + R);

			for (int y = yMin; y <= yMax; y++) {

				float dy = (float)(y - ry);
				float dySq = dy * dy;

				// Calculate how much squared radius is left for the X dimension
				float remainingSq = rSq - dzSq - dySq;

				// If < 0, this Y,Z coordinate is completely outside the sphere! Skip it.
				if (remainingSq < 0) continue;

				// Calculate the exact start and end of the sphere on this X-line
				int max_dx = (int)std::sqrt(remainingSq);

				int x_start = std::max(0, rx - max_dx);
				int x_end = std::min(nx - 1, rx + max_dx);

				// Threads processing different Z-planes can work simultaneously, lock the thread
				omp_set_lock(&z_locks[z]);

				for (int x = x_start; x <= x_end; x++) {
					size_t targetIdx = get_idx(x, y, z);

					// No distance check needed! Every 'x' in this loop is guaranteed inside.
					if (diameter > thicknessMap[targetIdx]) {
						thicknessMap[targetIdx] = diameter;
					}
				}

				omp_unset_lock(&z_locks[z]);
			}
		}
	}

	// clean up locks
	for (int z = 0; z < nz; z++) {
		omp_destroy_lock(&z_locks[z]);
	}

	float totalThicknessSum = 0.0f;
	size_t solidVoxelCount = 0;
	// Average the local thickness over the foreground phase only (field != 0),
	// matching BoneJ. Gating on thicknessMap > 0 instead would leak the thin
	// shell of background voxels that ridge spheres over-paint at their
	// boundary (the sphere of radius r reaches the nearest background voxel at
	// distance r) into both the mean and the sample count.

	#pragma omp parallel for reduction(+:totalThicknessSum, solidVoxelCount) collapse(3)
	for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx; x++) {
				size_t idx = get_idx(x, y, z);

				bool isInsideROI = true;
				if (con) {
					float px = blockBounds[0] + (x + 0.5f) * voxelSize;
					float py = blockBounds[2] + (y + 0.5f) * voxelSize;
					float pz = blockBounds[4] + (z + 0.5f) * voxelSize;
					isInsideROI = con->is_inside(Vec3(px, py, pz));
				}

				// Only accumulate foreground voxels physically inside the container ROI
				if (isInsideROI && field[idx] != 0) {
					totalThicknessSum += thicknessMap[idx];
					solidVoxelCount++;
				}
			}
		}
	}

	// Model-Independent Mean Thickness (Tb.Th for bone)
	float meanThicknessVoxels = (solidVoxelCount > 0) ? (totalThicknessSum / solidVoxelCount) : 0.0f;

	// estimate standard deviation
	float deviationSum = 0.0f;

	#pragma omp parallel for reduction(+:deviationSum) collapse(3)
	for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx; x++) {
				size_t idx = get_idx(x, y, z);

				bool isInsideROI = true;
				if (con) {
					float px = blockBounds[0] + (x + 0.5f) * voxelSize;
					float py = blockBounds[2] + (y + 0.5f) * voxelSize;
					float pz = blockBounds[4] + (z + 0.5f) * voxelSize;
					isInsideROI = con->is_inside(Vec3(px, py, pz));
				}

				if (isInsideROI && field[idx] != 0) {
					float diff = thicknessMap[idx] - meanThicknessVoxels;
					deviationSum += (diff * diff);
				}
			}
		}
	}

	float stdDevVoxels = (solidVoxelCount > 0) ? std::sqrt(deviationSum / solidVoxelCount) : 0.0f;

	if (separation) {
		localSeparation = meanThicknessVoxels * voxelSize;

		localSeparationStd = stdDevVoxels * voxelSize; // Final conversion to mm
	}
	else {
		localThickness = meanThicknessVoxels * voxelSize;

		localThicknessStd = stdDevVoxels * voxelSize; // Final conversion to mm
	}

	if (separation) {
		separationVersion = meshVersion;
		if (logger && !supress) logger->log(LogPriority::SUCCESS, "Estimated Local Separation!");
	}
	else {
		thicknessVersion = meshVersion;
		if (logger && !supress) logger->log(LogPriority::SUCCESS, "Estimated Local Thickness!");
	}
};

void GeneratorLewiner::build_phantom_field(
	int shape, float feature, float voxelSize, float domainSize) {

	// All phantoms are bounded in ALL three axes: an infinite slab is degenerate
	// for the separable distance transform (two axes have no background, so the
	// squared distance overflows). Each shape's feature size is set so its local
	// thickness equals 'feature' everywhere it is measured.
	if (domainSize <= 0.0f) {
		domainSize = std::max(1.0f, 6.0f * feature);
	}

	bounds = { 0.0f, domainSize, 0.0f, domainSize, 0.0f, domainSize };
	int n = std::max(4, static_cast<int>(std::ceil(domainSize / voxelSize)) + 1);
	blockDims = { n, n, n };
	update_steps();

	const float c = 0.5f * domainSize;   // domain centre
	const float half = 0.5f * feature;   // feature half-thickness / radius

	// A generous but inset lateral half-extent so ends never dominate the
	// measured thickness and a background margin exists for the transform.
	const float latHalf = std::min(0.4f * domainSize, std::max(4.0f * half, 0.3f * domainSize));
	// Torus: tube radius = half, ring radius a few tube-radii, kept inside.
	const float majorR = std::min(0.30f * domainSize, std::max(3.0f * half, 0.20f * domainSize));

	std::cout << "lathalf: " << latHalf << " , half: " << half << std::endl; 

	size_t total = static_cast<size_t>(n) * static_cast<size_t>(n) * static_cast<size_t>(n);
	scalarField.assign(total, 0.0f);

	// SDF-like field: solid where the shape distance is negative. The isoLevel
	// offset cancels in the get_image_field / marching_cubes threshold.
	#pragma omp parallel for collapse(3)
	for (int k = 0; k < n; ++k) {
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < n; ++i) {
				float x = bounds[0] + i * stepX;
				float y = bounds[2] + j * stepY;
				float z = bounds[4] + k * stepZ;
				float dx = x - c, dy = y - c, dz = z - c;

				float d;
				switch (shape) {
					case PHANTOM_CYLINDER: {
						// axis along z, diameter 'feature'
						float radial = std::sqrt(dx * dx + dy * dy) - half;
						float axial = std::fabs(dz) - latHalf;
						d = std::max(radial, axial);
						break;
					}
					case PHANTOM_SPHERE: {
						d = std::sqrt(dx * dx + dy * dy + dz * dz) - half;
						break;
					}
					case PHANTOM_TORUS: {
						// ring in the xy-plane; tube diameter 'feature'
						float q = std::sqrt(dx * dx + dy * dy) - majorR;
						d = std::sqrt(q * q + dz * dz) - half;
						break;
					}
					case PHANTOM_SLAB:
					default: {
						// finite plate, thin along z
						float ax = std::fabs(dx) - latHalf;
						float ay = std::fabs(dy) - latHalf;
						float az = std::fabs(dz) - half;
						d = std::max(ax, std::max(ay, az));
						break;
					}
				}
				scalarField[find_vertex_index(i, j, k)] = isoLevel + d;
			}
		}
	}

	// AABB spans the whole domain so estimators do not clip it, and no container
	// is attached so every voxel counts as inside the ROI.
	aabb.pMin = Vec3(0.0f, 0.0f, 0.0f);
	aabb.pMax = Vec3(domainSize, domainSize, domainSize);
	container.reset();

	// domainVolume so estimate_connectivity_density_voxel can normalize; volume
	// reset so estimate_smi recomputes BV from the fresh phantom mesh.
	domainVolume = domainSize * domainSize * domainSize;
	volume = 0.0f;
};

std::array<float, 2> GeneratorLewiner::run_slab_phantom(
	float slabThickness, float voxelSize, float domainSize) {

	build_phantom_field(PHANTOM_SLAB, slabThickness, voxelSize, domainSize);
	std::array<float, 6> blockBounds = bounds;
	estimate_local_thickness(voxelSize, blockBounds, false);
	return { localThickness, localThicknessStd };
};

bool GeneratorLewiner::secant_1d(
	const std::function<float(float)>& eval,
	float guess, float target, float tol, float lo, float hi, int maxIter,
	bool ascending, float& xOut,
	const std::string& label, int depth
){

	// default the output to the starting guess so xOut is always valid, even if
	// the very first evaluation already meets the tolerance.
	xOut = guess;

	// we need the following varibales
	float c = guess; // the current iteration value
	float cNext; // the next iteration value
	float cPrev = c; // the previous iteration value

	float m = eval(c); // this is the current metric, use the guess
	// A NaN means a sub-stage (assemble / marching cubes / a nested stage) failed
	// or the run was cancelled; abort cleanly instead of iterating on garbage.
	if (std::isnan(m)) { xOut = c; return false; }
	float mPrev = m; // previous iteration value

	for (int it{0}; it < maxIter; it++){

		// Cancellation: stop the ENTIRE calibration the instant the user cancels,
		// rather than finishing this (and every outer) secant first.
		if (cancelRequested && cancelRequested()) {
			if (logger) logger->log(LogPriority::WARNING,
				label + " calibration cancelled by user.");
			xOut = c;
			return false;
		}

		// per-iteration progress (mirrors the old calibrate_thickness/openness logs)
		if (logger) logger->log(LogPriority::INFO,
			std::string(2 * depth, ' ') + label + " [" + std::to_string(it) +
			"] knob=" + std::to_string(c) + " -> " + std::to_string(m) +
			" (target " + std::to_string(target) + ")");

		if (std::fabs(m - target) <= tol) {
			xOut = c;
			if (logger) logger->log(
				LogPriority::SUCCESS, label + " calibration converged.");
			return true;
		}

		if (it == 0) {
			// proportional first correction (handles the ~linear-through-origin part)
			cNext = c * target / std::max(m, 1e-6f);
		}
		else {
			float denom = m - mPrev; // f(c) - f(cPrev)

			// check if denom is zero
			if (std::fabs(denom) < 1e-9){
				float ratio = target / std::max(m, 1e-6f);
                ratio = std::clamp(ratio, 0.5f, 2.0f); 
                cNext = c * ratio;

				if (std::fabs(cNext - c) < 1e-6f) {
                    cNext = c * ((m < target) ? 1.1f : 0.9f);
                }
			}
			// else use the secant formula on the RESIDUAL f(c) = measure(c) - target
			// (using m alone would drive the metric to 0, not to the target).
			else{
				cNext = c - (m - target) * (c - cPrev) / denom;
			}
		}

		// clamp the NEXT candidate (not the already-evaluated current one)
		cNext = std::clamp(cNext, lo, hi);

		// Pinned-at-bound early-out: if we are already sitting on a bound and the
		// step cannot move us off it, the target lies outside the reachable range
		// (e.g. openness at 1 but the porosity target still higher). Continuing
		// only burns iterations and, for a nested solve, would return a spurious
		// failure that NaNs the outer stage. Stop and accept the closest feasible
		// value as best effort (converged in the reachable sense).
		if (cNext == c && (c == lo || c == hi)) {
			if (logger) logger->log(LogPriority::WARNING,
				label + " target unreachable within its range; using the closest "
				"feasible knob value (" + std::to_string(c) + ").");
			xOut = c;
			return true;
		}

		// update for next iteration
		cPrev = c; mPrev = m;
		c = cNext;
		m = eval(c);
		// sub-stage failed / cancelled at this candidate: keep the last good knob
		if (std::isnan(m)) { xOut = cPrev; return false; }
	}

	// max iterations reached: keep the last value, report whether it met tol
	xOut = c;
	return std::fabs(m - target) <= tol;
};

bool GeneratorLewiner::solve_calibration(
		std::vector<CalibrationStage>& stageSet,
		int lvl, 
		const IContainer& con
){

	// we have reached the level
	if (lvl == (int)stageSet.size()) return true;

	// Stop the whole nested solve immediately on cancel (before the expensive
	// cache build or the next secant).
	if (cancelRequested && cancelRequested()) {
		if (logger) logger->log(LogPriority::WARNING, "Calibration cancelled by user.");
		return false;
	}

	// Ensure the per-voxel cache exists before any Reassemble stage tries to
	// assemble_field(). A DA (Recache) outer stage rebuilds it itself, so only
	// pre-build when the outermost stage won't. Built once; inner stages reuse it.
	if (lvl == 0 && cachedVoxels.empty() &&
		(stageSet.empty() || stageSet[0].cost != Recompute::Recache)) {
		if (!compute_cached_field_values(con)) {
			if (logger) logger->log(LogPriority::ERROR,
				"Calibration aborted: scalar-field cache failed.");
			return false;
		}
	}

	// get the current stage
	CalibrationStage& k = stageSet[lvl];

	// add the logger info
	if (logger) logger->log(LogPriority::INFO,
        "[" + k.name + "] target " + std::to_string(k.target) +
        " (start " + std::to_string(k.get_knob()) + ")");

	// create a lambda that evaluates the stage. Explicit -> float so the NAN
	// (double) failure returns and the float metric returns don't conflict.
	auto eval = [&](float x) -> float {
		k.set_knob(x);

		// this is for DA if the cached 
		if(k.cost == Recompute::Recache && !compute_cached_field_values(con)) return NAN;

		// recursive call of solve if it did not converge return nan
		if(!solve_calibration(stageSet, lvl + 1, con)) return NAN;

		// chechk those the need assemble field
		if (lvl + 1 == (int)stageSet.size() && !assemble_field()) return NAN;
		
		// else return the measure
		return k.measure();
	};

	float xOut;
	bool ok = secant_1d(
		eval,
		k.get_knob(),
		k.target,
		k.tol, k.lo, k.hi, k.maxIter, k.ascending, xOut,
		k.name, lvl);

	// record convergence so the version can be stamped once the final mesh exists
	k.converged = ok;

	// get last value
	k.set_knob(xOut);

	// update logger
    if (ok && logger) logger->log(LogPriority::SUCCESS,
        "[" + k.name + "] converged: knob=" + std::to_string(xOut) +
        " metric=" + std::to_string(k.measure()));
    else if (!ok && logger) logger->log(LogPriority::ERROR,
        "[" + k.name + "] failed to converge");

	return ok;
};

// =============================================================================
// Calibration stage factories + assembly
// =============================================================================
// Each stage binds a scalar knob to its metric. solve_calibration drives them:
// an outer stage's every evaluation re-solves the inner stages, so coupling is
// absorbed automatically and warm starts (get_knob) keep each inner re-solve to
// 1-2 iterations. cost = Recache only for DA (stretch changes the warped seeds,
// so the per-voxel cache must be rebuilt); every other knob just re-assembles.

CalibrationStage GeneratorLewiner::make_thickness_stage(
	float target, float voxelSize, float tol, int maxIter) {

	CalibrationStage s;
	s.name = "Thickness";
	s.cost = Recompute::Recompute;          // iso-level / range scale -> re-assemble only
	s.target = target; s.tol = tol; s.maxIter = maxIter;
	s.lo = 1e-4f; s.hi = 100.0f;            // iso-level must stay positive
	s.ascending = true;                      // Tb.Th increases with iso-level

	// UNIFORM: knob is the iso-level. VARIED: knob is a scale on [start,end] that
	// preserves the grading ratio (mirrors calibrate_thickness).
	const bool  varied = (thicknessFunction && thicknessSDF);
	const float start0 = startThickness;
	const float end0   = endThickness;

	s.set_knob = [this, varied, start0, end0](float c) {
		if (varied) set_thickness_functions(thicknessSDF, thicknessFunction,
			c * start0, c * end0, transitionDistance);
		else        set_thickness(c);        // isoLevel = c
	};
	s.get_knob = [this, varied, start0, end0]() -> float {
		if (varied) {
			const float mean0 = 0.5f * (start0 + end0);
			return (std::fabs(mean0) > 1e-6f)
				? 0.5f * (startThickness + endThickness) / mean0 : 1.0f;
		}
		return isoLevel;
	};
	s.measure = [this, voxelSize]() -> float {
		std::array<float, 6> bb = bounds;    // measure on the generation domain...
		estimate_local_thickness(voxelSize, bb, false, true); // ...at the measurement voxel
		return localThickness;
	};
	s.stampVersion = [this]() { thicknessVersion = meshVersion; };
	return s;
}

CalibrationStage GeneratorLewiner::make_porosity_stage(
	float target, float tol, int maxIter) {

	CalibrationStage s;
	s.name = "Porosity";
	s.cost = Recompute::Recompute;
	s.target = target; s.tol = tol; s.maxIter = maxIter;
	s.lo = 0.0f; s.hi = 1.0f;                // openness in [0,1]
	s.ascending = true;                      // more openness -> more porosity
	s.set_knob = [this](float x) { threshold = std::clamp(x, 0.0f, 1.0f); };
	s.get_knob = [this]() -> float { return threshold; };
	// porosity (voxel-based) is set by assemble_field, which solve_calibration
	// runs for the innermost stage before measure() - no mesh needed.
	s.measure  = [this]() -> float { return porosity; };
	return s;
}

// DISABLED: SMI is not a reliable calibration target -- it is governed by openness
// (which also fixes porosity) and cannot be steered independently at fixed
// radius/thickness; spread is too weak/confounded to move it. Kept (disabled) in
// case a coupled SMI solve is revisited. SMI is a reported outcome, tuned by hand
// via openness / radius / thickness. See the parameter-control discussion.
#if 0
CalibrationStage GeneratorLewiner::make_smi_stage(
	float target, float tol, int maxIter) {

	CalibrationStage s;
	s.name = "SMI";
	s.cost = Recompute::Recompute;
	s.target = target; s.tol = tol; s.maxIter = maxIter;
	s.lo = 0.0f; s.hi = 1.0f;                // spread in [0,1]
	s.ascending = true;                      // SMI rises with spread at matched porosity
	s.set_knob = [this](float x) { spread = std::clamp(x, 0.0f, 1.0f); };
	s.get_knob = [this]() -> float { return spread; };
	s.measure  = [this]() -> float {
		if (!marching_cubes(true)) return NAN;
		return estimate_smi();
	};
	s.stampVersion = [this]() { smiVersion = meshVersion; };
	return s;
}
#endif

CalibrationStage GeneratorLewiner::make_da_stage(
	float target, float voxelSize, float tol, int maxIter, int mode) {

	CalibrationStage s;
	s.name = "DA";
	s.cost = Recompute::Recache;             // stretch changes warped seeds -> rebuild cache
	s.target = target; s.tol = tol; s.maxIter = maxIter;
	s.lo = 1.0f; s.hi = 8.0f;                // along-fibre stretch (>=1 elongates)
	s.ascending = true;                      // DA increases with along-fibre stretch
	// The anisotropy axis is rotated onto local X, so stretchX is the along-fibre
	// stretch; keep Y/Z as-is so DA is a single scalar knob.
	s.set_knob = [this](float x) { set_stretch(x, stretchY, stretchZ); };
	s.get_knob = [this]() -> float { return stretchX; };
	s.measure  = [this, voxelSize, mode]() -> float {
		if (!marching_cubes(true)) return NAN;
		// Reduced ray budget inside the loop (MIL ray-casting dominates cost);
		// a full-resolution DA is measured once after solve_calibration returns.
		estimate_anisotropy(voxelSize, /*dirs*/ 500, /*lines*/ 2000, mode, nullptr);
		return anisotropyDegree;
	};
	s.stampVersion = [this]() { anisotropyVersion = meshVersion; };
	return s;
}

void GeneratorLewiner::build_calibration_stages(float voxelSize) {

	stages.clear();

	const float tol = calibrationTol;
	const int   it  = calibrationIter;

	// Priority order (outer -> inner). DA/SMI use a looser tolerance because their
	// metrics are coarser/noisier than the resampled Tb.Th and voxel porosity.
	if (calibrateDA)
		stages.push_back(make_da_stage(targetDa, voxelSize, std::max(0.02f, 10.0f * tol), it, targetFormulaIdx));
	// SMI is intentionally NOT a calibration target: it is set by openness (which
	// also fixes porosity) and cannot be steered independently at fixed
	// radius/thickness. It is reported, and tuned by hand via openness / radius /
	// thickness (see the parameter-control discussion in the paper). Disabled:
	// if (calibrateSmi)
	//     stages.push_back(make_smi_stage(targetSmi, std::max(0.05f, 10.0f * tol), it));
	if (calibrateThickness)
		stages.push_back(make_thickness_stage(targetThickness, voxelSize, tol, it));
	if (calibratePorosity)
		stages.push_back(make_porosity_stage(targetPorosity * 0.01f, tol, it)); // percent -> fraction
}

void GeneratorLewiner::stamp_calibrated_versions() {
	for (auto& s : stages)
		if (s.converged && s.stampVersion) s.stampVersion();
}


bool GeneratorLewiner::calibrate_thickness(
	const IContainer& con, float targetThickness, float voxelSize,
	float tol, int maxIter, const std::function<void(float)>& onProgress) {

	if (targetThickness <= 0.0f) {
		if (logger) logger->log(LogPriority::ERROR, "Calibration target thickness must be positive.");
		return false;
	}

	// UNIFORM thickness: the search variable c IS the iso-level. VARIED thickness
	// (a graded iso-level from thicknessFunction over [startThickness,endThickness])
	// ignores isoLevel entirely, so instead c is a multiplicative SCALE on the
	// whole range - preserving the grading ratio while shifting its magnitude.
	const bool varied = (thicknessFunction && thicknessSDF);
	const float start0 = startThickness;
	const float end0 = endThickness;

	// The knob-invariant per-voxel cache (the expensive kdtree pass) is built once
	// and reused every iteration through assemble_field. If a caller (two_knob, the
	// GUI) already populated it we reuse it; otherwise self-cache so standalone use
	// and the profiler work. set_seeds/set_bounds/set_resolution clear it, so a
	// stale cache from a different config is never silently reused.
	if (cachedVoxels.empty() && !compute_cached_field_values(con)) {
		if (logger) logger->log(LogPriority::ERROR, "Scalar field cache failed.");
		return false;
	}

	// Generate for search value c and return the resulting MEASURED mean Tb.Th
	// (mm), or -1 on generation failure. Leaves the field assembled at c (no mesh;
	// standalone callers run marching_cubes once after convergence).
	auto measure = [&](float c) -> float {
		if (varied) {
			set_thickness_functions(thicknessSDF, thicknessFunction,
				c * start0, c * end0, transitionDistance);
		}
		else {
			set_thickness(c); // isoLevel = c
		}
		if (!assemble_field()) return -1.0f;

		// Measure over the generation grid (not the mesh AABB / get_bounds(), which
		// flips to the aabb once a mesh exists) so the calibration substrate is
		// stable and matches the coupling re-check / reported Tb.Th.
		std::array<float, 6> blockBounds = bounds;
		estimate_local_thickness(voxelSize, blockBounds, false, true);
		return localThickness;
	};

	auto report = [&](float p) { if (onProgress) onProgress(std::clamp(p, 0.0f, 1.0f)); };

	// First guess. Uniform: the target iso-level itself. Varied: the scale that
	// puts the mean of the range on the target (c * mean(range) = target).
	float c;
	if (varied) {
		float mean0 = 0.5f * (start0 + end0);
		c = (std::fabs(mean0) > 1e-6f) ? (targetThickness / mean0) : 1.0f;
	}
	else {
		c = targetThickness;
	}

	// The c -> Tb.Th map is monotonic increasing, so a secant search over the
	// actual (re-measured) output converges without a fitted model.
	float m = measure(c);
	if (m < 0.0f) {
		if (logger) logger->log(LogPriority::ERROR, "Calibration aborted: generation failed.");
		return false;
	}
	report(1.0f / (maxIter + 1.0f));

	float cPrev = c, mPrev = m;
	for (int it = 0; it < maxIter; ++it) {
		if (cancelRequested && cancelRequested()) {
			if (logger) logger->log(LogPriority::WARNING, "Thickness calibration cancelled by user.");
			report(1.0f);
			return true; // leave the last valid built scaffold
		}
		std::ostringstream oss;
		oss << "Calibrate [" << it << "] " << (varied ? "scale=" : "c_glob=") << c
			<< " -> Tb.Th=" << m << " (target " << targetThickness << ")";
		if (logger) logger->log(LogPriority::INFO, oss.str());

		if (std::fabs(m - targetThickness) <= tol) {
			if (logger) logger->log(LogPriority::SUCCESS, "Thickness calibration converged.");
			report(1.0f);
			return true;
		}

		float cNext;
		if (it == 0) {
			// proportional first correction (handles the ~linear-through-origin part)
			cNext = c * targetThickness / std::max(m, 1e-6f);
		}
		else {
			float denom = m - mPrev;

			if (std::fabs(denom) < 1e-9){
				float ratio = targetThickness / std::max(m, 1e-6f);
                ratio = std::clamp(ratio, 0.5f, 2.0f); 
                cNext = c * ratio;

				if (std::fabs(cNext - c) < 1e-6f) {
                    cNext = c * ((m < targetThickness) ? 1.1f : 0.9f);
                }
			}
			else{
				cNext = c - (m - targetThickness) * (c - cPrev) / denom;
			}
		}

		// keep positive and inside a generous bracket. For uniform c is an
		// iso-level (~target); for varied c is a dimensionless scale (~1).
		cNext = varied ? std::clamp(cNext, 0.01f, 100.0f)
			: std::clamp(cNext, 0.1f * targetThickness, 10.0f * targetThickness);

		cPrev = c; mPrev = m;
		c = cNext;
		m = measure(c);
		if (m < 0.0f) {
			// restore the last good build so the generator holds a valid scaffold
			if (logger) logger->log(LogPriority::WARNING,
				"Calibration step failed; restoring the previous valid scaffold.");
			bool restored = (measure(cPrev) >= 0.0f);
			report(1.0f);
			return restored;
		}
		report((it + 2.0f) / (maxIter + 1.0f));
	}

	if (logger) logger->log(LogPriority::WARNING,
		"Thickness calibration hit the iteration limit; returning best estimate.");
	report(1.0f);
	return true; // generator is left built at this c (last measured)
};

bool GeneratorLewiner::calibrate_openness(
	const IContainer& con,
	float targetPorosity,
	float voxelSize,
	float tol, int maxIter,
	const std::function<void(float)>& onProgress) {

	if (targetPorosity <= 0.0f || targetPorosity >= 1.0f) {
		if (logger) logger->log(
			LogPriority::ERROR, "Porosity target must be in (0,1).");
		return false;
	}

	// Reuse the caller's cache (two_knob, GUI) or self-cache for standalone use.
	if (cachedVoxels.empty() && !compute_cached_field_values(con)) {
		if (logger) logger->log(LogPriority::ERROR, "Scalar field cache failed.");
		return false;
	}

	// Vary the openness (threshold) at the CURRENT, fixed thickness (iso-level or
	// graded range). Porosity is the voxel-based estimate that assemble_field
	// sets, so no marching cubes is needed to measure it. Porosity is monotone
	// increasing in openness (foam -> lattice), so a secant converges.
	auto measure = [&](float op) -> float {
		threshold = std::clamp(op, 0.0f, 1.0f);
		if (!assemble_field()) return -1.0f;
		return porosity;
	};

	auto report = [&](float p) { if (onProgress) onProgress(std::clamp(p, 0.0f, 1.0f)); };

	// Start the search from the NEUTRAL midpoint, not the current openness. A
	// previous calibration may have left threshold pinned at a bound (0 or 1),
	// from which the secant + "unreachable" guard would fail immediately even
	// though the target is reachable. 0.5 is a robust, well-conditioned start.
	float c = 0.5f;
	float m = measure(c);
	if (m < 0.0f) {
		if (logger) logger->log(LogPriority::ERROR, "Openness calibration aborted: generation failed.");
		return false;
	}
	report(1.0f / (maxIter + 1.0f));

	float cPrev = c, mPrev = m;
	for (int it = 0; it < maxIter; ++it) {
		if (cancelRequested && cancelRequested()) {
			if (logger) logger->log(LogPriority::WARNING, "Openness calibration cancelled by user.");
			report(1.0f);
			return true; // leave the last valid field
		}
		std::ostringstream oss;
		oss << "Calibrate openness [" << it << "] tau=" << c
			<< " -> porosity=" << m << " (target " << targetPorosity << ")";
		if (logger) logger->log(LogPriority::INFO, oss.str());

		if (std::fabs(m - targetPorosity) <= tol) {
			if (logger) logger->log(LogPriority::SUCCESS, "Openness calibration converged.");
			report(1.0f);
			return true;
		}

		float cNext;
		if (it == 0) {
			// nudge in the correct monotone direction (porosity increases with tau)
			cNext = c + ((targetPorosity > m) ? 0.1f : -0.1f);
		}
		else {
			float denom = m - mPrev;
			if (std::fabs(denom) < 1e-9f) {
				float nudge = (targetPorosity > m) ? 0.05f : -0.05f;
                cNext = c + nudge;
			}
			else{
				cNext = c - (m - targetPorosity) * (c - cPrev) / denom;
			}
		}
		cNext = std::clamp(cNext, 0.0f, 1.0f);   // openness is a fraction in [0,1]

		// If we are pinned at a bound and still off-target, the porosity is
		// unreachable at this thickness/spacing - stop and report.
		if (cNext == cPrev && (c == 0.0f || c == 1.0f)) {
			if (logger) logger->log(LogPriority::WARNING,
				"Target porosity unreachable within openness [0,1] at this thickness/spacing.");
			report(1.0f);
			return true;
		}

		cPrev = c; mPrev = m;
		c = cNext;
		m = measure(c);
		if (m < 0.0f) {
			bool restored = (measure(cPrev) >= 0.0f);
			report(1.0f);
			return restored;
		}
		report((it + 2.0f) / (maxIter + 1.0f));
	}

	if (logger) logger->log(LogPriority::WARNING,
		"Openness calibration hit the iteration limit; returning best estimate.");
	report(1.0f);
	return true;
};

namespace {
// RAII restore of the junction-smoothing knobs (smoothJunctions/smoothK). The
// joint calibrations may back smoothK off during the thickness solve; on any
// failure/cancel path the object must return to its pre-call state, but on a
// SUCCESSFUL calibration the backed-off value is part of the result and must
// persist (so a later regeneration reproduces the calibrated Tb.Th). Restores on
// scope exit unless commit() was called. Replaces the hand-written cache/restore
// that had to be repeated at every early return.
struct SmoothKnobGuard {
	bool& jRef; float& kRef; const bool origJ; const float origK; bool committed = false;
	SmoothKnobGuard(bool& j, float& k) : jRef(j), kRef(k), origJ(j), origK(k) {}
	SmoothKnobGuard(const SmoothKnobGuard&) = delete;
	SmoothKnobGuard& operator=(const SmoothKnobGuard&) = delete;
	void commit() { committed = true; }
	~SmoothKnobGuard() { if (!committed) { jRef = origJ; kRef = origK; } }
};
} // namespace

// Two-knob calibration (Tb.Th, porosity). The (thickness, openness) -> (Tb.Th,
// porosity) Jacobian is near-triangular (dTb.Th/dopenness ~ 0, measured), so the
// solve is a back-substitution: calibrate thickness, then openness. The outer
// loop absorbs the tiny residual coupling (an openness change nudges Tb.Th by
// ~the noise floor) and normally breaks after one pass.
bool GeneratorLewiner::two_knob_calibration(
	const IContainer& con,
	float targetThickness,
	float targetPorosity,
	float voxelSize,
	float tol, int maxIter,
	const std::function<void(float)>& onProgress) {

	auto report = [&](float p) { if (onProgress) onProgress(std::clamp(p, 0.0f, 1.0f)); };

	// Restore smoothJunctions/smoothK on any early return; commit() (just before
	// the final success return) keeps the possibly backed-off values.
	SmoothKnobGuard knobGuard(smoothJunctions, smoothK);

	// Reset the knobs to neutral so the solve does not inherit a bad state (e.g. openness pinned at a bound) from a previous calibration. Thickness is reset inside calibrate_thickness (first guess = target); reset openness here so the thickness step also runs at a well-conditioned, mid-range openness.
	threshold = 0.5f;

	// The system is near-triangular (dTb.Th/dopenness ~ 0), so ONE outer pass
	// suffices; cap at 3 so a tight tolerance can never cause a maxIter-pass
	// runaway (each pass is ~20 generations).
	const int outerMax = std::min(3, maxIter);
	for (int it = 0; it < outerMax; ++it) {
		if (cancelRequested && cancelRequested()) break; // sub-solves already bailed; stop looping

		if (logger) {
            std::ostringstream oss;
            oss << "=== Two-knob search -> Iter " << (it + 1) << " of " << outerMax << " ===";
            logger->log(LogPriority::INFO, oss.str());
            logger->log(LogPriority::INFO, "Step 1: Calibrating thickness (iso-level/scale).");
        }

		// 1. thickness -> Tb.Th (independent of openness); leaves the mesh built.
		if (!calibrate_thickness(con, targetThickness, voxelSize, tol, 8,
            [&](float p) { report(0.05f + 0.40f * p); }))
            return false;   // knobGuard restores

		// Junction smoothing (the smin fillet of radius k = smoothK) fattens the
		// rod/plate nodes, which raises the MINIMUM Tb.Th the geometry can reach.
		// If k is too large that floor sits ABOVE the target, so the secant in
		// calibrate_thickness drives the iso-level down yet Tb.Th stays pinned
		// high - it "fails" (returns built, but off-target). The cure is to shrink
		// the fillet: back k off geometrically and recalibrate until Tb.Th is
		// reachable, or k is effectively zero (the plain linear blend, whose floor
		// is the lowest achievable). Only runs when smoothing is actually on.
		const int kBackoffMax = 6;
		for (int r = 0; smoothJunctions && smoothK > 1e-6f
			&& std::fabs(localThickness - targetThickness) > tol
			&& r < kBackoffMax; ++r) {

			const float kPrev = smoothK;
			smoothK = (0.5f * smoothK > 1e-6f) ? 0.5f * smoothK : 0.0f;
			if (smoothK == 0.0f) smoothJunctions = false; // fell through to linear blend

			if (logger) {
                std::ostringstream oss;
                oss << "Tb.Th target " << targetThickness
                    << " unreachable at smoothK=" << kPrev << " (measured "
                    << localThickness << "); reducing k to " << smoothK
                    << " and recalibrating thickness.";
                logger->log(LogPriority::WARNING, oss.str());
            }

			if (!calibrate_thickness(con, targetThickness, voxelSize, tol, 8,
                [&](float p) { report(0.05f + 0.40f * p); }))
                return false;   // knobGuard restores
		}

		if (logger) {
            logger->log(LogPriority::INFO, "Step 2: Calibrating openness (porosity).");
        }

		// 2. openness -> porosity at the fixed thickness (field only, no mesh).
		if (!calibrate_openness(con, targetPorosity, voxelSize, tol, maxIter,
            [&](float p) { report(0.45f + 0.40f * p); }))
            return false;   // knobGuard restores

		// 3. re-measure Tb.Th after the openness change to check the coupling. This
		// is field-based (no mesh needed) and uses the SAME generation-grid
		// substrate as calibrate_thickness, so the coupling check and the reported
		// Tb.Th are consistent. The final mesh is built once after the loop.
		std::array<float, 6> bnds = bounds;
		estimate_local_thickness(voxelSize, bnds, false, true);

		// Break as soon as Tb.Th is on target. calibrate_openness already drove
		// openness to its optimum for this thickness (converged, or pinned at a
		// bound if the porosity target is infeasible), so once the thickness drift
		// is gone another pass would only REPEAT the identical openness/porosity.
		bool okT = std::fabs(localThickness - targetThickness) <= tol;
		if (okT) {
			if (std::fabs(porosity - targetPorosity) > tol && logger)
				logger->log(LogPriority::WARNING,
					"two_knob: porosity target not reachable at the calibrated Tb.Th; "
					"returning best-effort openness.");
			break;
		}
		else if (it + 1 < outerMax) {
            // LOG: Explain the restart reason
            if (logger) {
                std::ostringstream oss;
                oss << "Coupling drift detected: changing openness shifted Tb.Th to "
                    << localThickness << " (target " << targetThickness 
                    << "). Restarting outer loop to correct.";
                logger->log(LogPriority::INFO, oss.str());
            }
        }
	}

	// If cancelled, skip the heavy final pass; the last valid scaffold is built. knobGuard restores the original knobs (cancel discards the backed-off k).
	if (cancelRequested && cancelRequested()) {
        if (logger) logger->log(LogPriority::WARNING, "Two-knob calibration cancelled by user.");
        report(1.0f);
        return true;
    }

	// final: full metric set on the calibrated scaffold.
	if (!marching_cubes()) return false;   // knobGuard restores
	estimate_metrics(con);
	report(1.0f);

	thicknessVersion = meshVersion;
	// success: keep the calibrated (possibly backed-off) smoothK
	knobGuard.commit();   
	return true;
};

// Three-knob calibration (Tb.Th, porosity, SMI). Thickness is calibrated and
// frozen (it decouples from openness/spread). The (openness, spread)->(porosity,
// SMI) block is ILL-CONDITIONED (the knobs move both outputs nearly parallel), so
// instead of inverting a near-singular 2x2 we solve it NESTED: secant on spread
// to hit SMI, with openness re-calibrated inside each step to hold porosity. Two
// robust monotone 1-D solves; porosity is held exactly, SMI is best-effort within
// its reachable band (see §8.4 - SMI is only weakly controllable at fixed density).
bool GeneratorLewiner::three_knob_calibration(
	const IContainer& con,
	float targetThickness,
	float targetPorosity,
	float targetSMI,
	float voxelSize,
	float tol, int maxIter,
	const std::function<void(float)>& onProgress) {

	auto report = [&](float p) { if (onProgress) onProgress(std::clamp(p, 0.0f, 1.0f)); };

	// Restore smoothJunctions/smoothK on any early return; commit() (just before
	// the final success return) keeps the possibly backed-off values.
	SmoothKnobGuard knobGuard(smoothJunctions, smoothK);

	// Neutral starting knobs so the solve does not inherit a bad accumulated state.
	threshold = 0.5f;
	spread = 0.5f;

	if (logger) {
		logger->log(LogPriority::INFO, "=== Three-knob search ===");
		logger->log(LogPriority::INFO, "Step 1: Calibrating thickness (iso-level/scale).");
	}

	// 1. thickness -> Tb.Th, frozen.
	if (!calibrate_thickness(con, targetThickness, voxelSize, tol, 8,
		[&](float p) { report(0.05f + 0.25f * p); }))
		return false;   // knobGuard restores

	// Junction-fillet back-off (identical to two_knob_calibration): the smin fillet
	// of radius k = smoothK fattens the nodes and raises the MINIMUM reachable
	// Tb.Th; if that floor sits above the target, shrink k and recalibrate until the
	// target is reachable, or k underflows to the plain linear blend.
	const int kBackoffMax = 6;
	for (int r = 0; smoothJunctions && smoothK > 1e-6f
		&& std::fabs(localThickness - targetThickness) > tol
		&& r < kBackoffMax; ++r) {

		const float kPrev = smoothK;
		smoothK = (0.5f * smoothK > 1e-6f) ? 0.5f * smoothK : 0.0f;
		if (smoothK == 0.0f) smoothJunctions = false; // fell through to linear blend

		if (logger) {
			std::ostringstream oss;
			oss << "Tb.Th target " << targetThickness
				<< " unreachable at smoothK=" << kPrev << " (measured "
				<< localThickness << "); reducing k to " << smoothK
				<< " and recalibrating thickness.";
			logger->log(LogPriority::WARNING, oss.str());
		}

		if (!calibrate_thickness(con, targetThickness, voxelSize, tol, 8,
			[&](float p) { report(0.05f + 0.25f * p); }))
			return false;   // knobGuard restores
	}

	if (logger) logger->log(LogPriority::INFO,
		"Step 2: Calibrating openness + spread (porosity, SMI).");

	// Step 2: nested 1-D solve, robust to the (openness, spread) ILL-CONDITIONING.
	// openness and spread move (porosity, SMI) along nearly parallel directions,
	// so the 2x2 Jacobian is near-singular and a direct Newton/LM overshoots.
	// Instead we exploit the structure the sensitivity study found: SMI(spread) AT
	// MATCHED POROSITY is monotone. So we secant on SPREAD to hit SMI, and each
	// evaluation re-calibrates OPENNESS (1-D, well-conditioned) to hold porosity.
	// Two nested robust 1-D solves, no matrix inversion; porosity is held EXACTLY
	// (a firm literature target) and SMI is best-effort within its reachable band.
	const float tolS = std::max(0.05f, 10.0f * tol);

	// spread -> SMI, with openness re-solved to keep porosity on target. Leaves
	// the generator built (field + mesh) at the calibrated (openness, spread).
	auto smi_at_target_porosity = [&](float sp) -> float {
		spread = std::clamp(sp, 0.0f, 1.0f);
		if (!calibrate_openness(con, targetPorosity, voxelSize, tol, maxIter, nullptr))
			return -1e9f;
		if (!marching_cubes(true)) return -1e9f;
		estimate_metrics(con);
		estimate_smi();
		return smi;
	};

	float sp = spread;
	float S = smi_at_target_porosity(sp);
	if (S < -1e8f) {
		if (logger) logger->log(LogPriority::ERROR, "three_knob: openness sub-calibration failed.");
		return false;
	}
	report(0.4f);

	float spPrev = sp, SPrev = S;
	for (int it = 0; it < maxIter; ++it) {
		if (cancelRequested && cancelRequested()) {
			if (logger) logger->log(LogPriority::WARNING, "Three-knob calibration cancelled by user.");
			break;
		}
		std::ostringstream oss;
		oss << "3-knob (nested) [" << it << "] spread=" << sp << " -> SMI=" << S
			<< " (target " << targetSMI << ", porosity held at " << porosity << ")";
		if (logger) logger->log(LogPriority::INFO, oss.str());

		if (std::fabs(S - targetSMI) <= tolS) {
			if (logger) logger->log(LogPriority::SUCCESS, "Three-knob calibration converged.");
			break;
		}

		float spNext;
		if (it == 0) {
			// SMI increases with spread at matched porosity (measured, see 8.4)
			spNext = sp + ((targetSMI > S) ? 0.1f : -0.1f);
		}
		else {
			float denom = S - SPrev;
			spNext = (std::fabs(denom) < 1e-9f)
				? 0.5f * (sp + spPrev)
				: sp - (S - targetSMI) * (sp - spPrev) / denom;
		}
		spNext = std::clamp(spNext, 0.0f, 1.0f);

		if (spNext == spPrev && (sp == 0.0f || sp == 1.0f)) {
			if (logger) logger->log(LogPriority::WARNING,
				"three_knob: target SMI outside the reachable band at this porosity; "
				"returning closest feasible spread.");
			break;
		}

		spPrev = sp; SPrev = S;
		sp = spNext;
		S = smi_at_target_porosity(sp);
		if (S < -1e8f) { smi_at_target_porosity(spPrev); break; }  // restore last good
		report(0.4f + 0.5f * (it + 1.0f) / maxIter);
	}

	// generator already holds the calibrated scaffold; final full metric pass so
	// all reported metrics are consistent (skipped if the user cancelled).
	// knobGuard restores on cancel/failure; commit() keeps the calibrated knobs.
	if (cancelRequested && cancelRequested()) { report(1.0f); return true; }
	if (!marching_cubes()) return false;
	estimate_metrics(con);
	report(1.0f);
	knobGuard.commit();   // success: keep the calibrated (possibly backed-off) smoothK
	return true;
};

//@Function to get a subregion of the created mesh to compute the image metrics, we also should add
// an origin (e.g. the centroid).
std::vector<uint8_t> GeneratorLewiner::get_image_field(
	float voxelSize, std::array<float, 6>& blockBounds, bool inverse, uint8_t solidValue) {

	float sizeX = blockBounds[1] - blockBounds[0];
	float sizeY = blockBounds[3] - blockBounds[2];
	float sizeZ = blockBounds[5] - blockBounds[4];

	// get number of voxels
	size_t nx = static_cast<size_t>(std::ceil(sizeX / voxelSize));
	size_t ny = static_cast<size_t>(std::ceil(sizeY / voxelSize));
	size_t nz = static_cast<size_t>(std::ceil(sizeZ / voxelSize));

	float step = voxelSize;
	size_t totalVoxels = nx * ny * nz;

	std::vector<uint8_t> field(totalVoxels, 0);

	auto clamp_idx = [](int val, int maxVal) {
		if (val < 0) return 0;
		if (val >= maxVal) return maxVal - 1;
		return val;
	};

	#pragma omp parallel for
	for (long long idx = 0; idx < static_cast<long long>(totalVoxels); idx++) {

		// 2. Use size_t for the extracted grid coordinates
		size_t z = idx / (nx * ny);
		size_t y = (idx % (nx * ny)) / nx;
		size_t x = idx % nx;

		// find physical position
		float pX = blockBounds[0] + x * voxelSize;
		float pY = blockBounds[2] + y * voxelSize;
		float pZ = blockBounds[4] + z * voxelSize;

		// find index of voxel in the scalar field and then interpolate
		float oldX = (pX - bounds[0]) / stepX;
		float oldY = (pY - bounds[2]) / stepY;
		float oldZ = (pZ - bounds[4]) / stepZ;

		if (oldX < 0.0f || oldX >= (blockDims[0] - 1.0f) ||
			oldY < 0.0f || oldY >= (blockDims[1] - 1.0f) ||
			oldZ < 0.0f || oldZ >= (blockDims[2] - 1.0f)) {
			field[idx] = 0;
			continue;
		}

		// this is like 5.1, 7.8, 9.1, we need to find the interpolated value of the scalar field to assing 0 or 255
		// get the 8 corners
		int x0 = static_cast<int>(std::floor(oldX));
		int y0 = static_cast<int>(std::floor(oldY));
		int z0 = static_cast<int>(std::floor(oldZ));

		int x1 = clamp_idx(x0 + 1, blockDims[0]);
		int y1 = clamp_idx(y0 + 1, blockDims[1]);
		int z1 = clamp_idx(z0 + 1, blockDims[2]);

		// this measures how far we are from x0,y0,z0
		float tx = oldX - x0;
		float ty = oldY - y0;
		float tz = oldZ - z0;

		double c000 = scalarField[find_vertex_index(x0, y0, z0)];
		double c100 = scalarField[find_vertex_index(x1, y0, z0)];
		double c010 = scalarField[find_vertex_index(x0, y1, z0)];
		double c110 = scalarField[find_vertex_index(x1, y1, z0)];
		double c001 = scalarField[find_vertex_index(x0, y0, z1)];
		double c101 = scalarField[find_vertex_index(x1, y0, z1)];
		double c011 = scalarField[find_vertex_index(x0, y1, z1)];
		double c111 = scalarField[find_vertex_index(x1, y1, z1)];

		double c00 = c000 * (1 - tx) + c100 * tx;
		double c10 = c010 * (1 - tx) + c110 * tx;
		double c01 = c001 * (1 - tx) + c101 * tx;
		double c11 = c011 * (1 - tx) + c111 * tx;

		// 2. Interpolate along Y (Reduce 4 points to 2)
		double c0 = c00 * (1 - ty) + c10 * ty;
		double c1 = c01 * (1 - ty) + c11 * ty;

		// 3. Interpolate along Z (Reduce 2 points to 1)
		double val = c0 * (1 - tz) + c1 * tz;

		bool isSolid = (val < isoLevel);

		if (inverse) {
			isSolid = !isSolid;
		}

		field[idx] = isSolid ? solidValue : 0;
	}

	return field;
};

void GeneratorLewiner::export_scaf(const std::string& fileName) {

	std::ofstream out(fileName, std::ios::out | std::ios::binary);
	if (!out.is_open()) {
		logger->log(LogPriority::ERROR, "Failed to open file for binary serialization: " + fileName);
		return;
	}

	// 1. Header Block
	char header[4] = { 'S', 'C', 'A', 'F' };
	out.write(header, 4);
	uint32_t fileVersion = 3;   // v2 adds the extension block (spread, junction/frame, calibration); v3 appends edge rounding (roundEdges, edgeK)
	out.write(reinterpret_cast<const char*>(&fileVersion), sizeof(fileVersion));

	// 2. Core Parameters Block
	out.write(reinterpret_cast<const char*>(&selectedThicknessOption), sizeof(selectedThicknessOption));
	out.write(reinterpret_cast<const char*>(&selectedDist), sizeof(selectedDist));
	out.write(reinterpret_cast<const char*>(&selectedFunc), sizeof(selectedFunc));
	out.write(reinterpret_cast<const char*>(&startThickness), sizeof(startThickness));
	out.write(reinterpret_cast<const char*>(&endThickness), sizeof(endThickness));
	out.write(reinterpret_cast<const char*>(&transitionDistance), sizeof(transitionDistance));
	out.write(reinterpret_cast<const char*>(&threshold), sizeof(threshold));
	out.write(reinterpret_cast<const char*>(&isoLevel), sizeof(isoLevel));
	out.write(reinterpret_cast<const char*>(&voxelSize), sizeof(voxelSize));
	
	// smoothness features
	out.write(reinterpret_cast<const char*>(&iter), sizeof(iter));
	out.write(reinterpret_cast<const char*>(&lambda), sizeof(lambda));
	out.write(reinterpret_cast<const char*>(&mu), sizeof(mu));
	
	out.write(reinterpret_cast<const char*>(&stretchX), sizeof(stretchX));
	out.write(reinterpret_cast<const char*>(&stretchY), sizeof(stretchY));
	out.write(reinterpret_cast<const char*>(&stretchZ), sizeof(stretchZ));
	out.write(reinterpret_cast<const char*>(&anisotropyAngle), sizeof(anisotropyAngle));

	out.write(reinterpret_cast<const char*>(&anisotropyVec.x), sizeof(float));
	out.write(reinterpret_cast<const char*>(&anisotropyVec.y), sizeof(float));
	out.write(reinterpret_cast<const char*>(&anisotropyVec.z), sizeof(float));

	out.write(reinterpret_cast<const char*>(&renderMode), sizeof(renderMode));
	out.write(reinterpret_cast<const char*>(&backgroundWeight), sizeof(backgroundWeight));

	// Anisotropy Sources Block
	uint32_t numSources = static_cast<uint32_t>(anisotropySources.size());
    out.write(reinterpret_cast<const char*>(&numSources), sizeof(numSources));
    
    for (const auto& src : anisotropySources) {
        // Write the properties required to reconstruct the source
        out.write(reinterpret_cast<const char*>(&src->origin.x), sizeof(float) * 3);
        out.write(reinterpret_cast<const char*>(&src->direction.x), sizeof(float) * 3);
        out.write(reinterpret_cast<const char*>(&src->stretch.x), sizeof(float) * 3);
        out.write(reinterpret_cast<const char*>(&src->sigma), sizeof(float));
        out.write(reinterpret_cast<const char*>(&src->angle), sizeof(float));         
    }

	// 3. Container Context Block
	std::shared_ptr<IContainer> lockedCon = container.lock();
	int32_t containerTypeID = 0;

	if (lockedCon) {
		if (lockedCon->get_type() == ObjectType::BoxContainerType) {
			containerTypeID = 1;
			out.write(reinterpret_cast<const char*>(&containerTypeID), sizeof(containerTypeID));
			auto* box = static_cast<BoxContainer*>(lockedCon.get());

			// Fixed: Component-wise Vec3 decomposition
			out.write(reinterpret_cast<const char*>(&box->size.x), sizeof(float));
			out.write(reinterpret_cast<const char*>(&box->size.y), sizeof(float));
			out.write(reinterpret_cast<const char*>(&box->size.z), sizeof(float));
			out.write(reinterpret_cast<const char*>(&box->origin.x), sizeof(float));
			out.write(reinterpret_cast<const char*>(&box->origin.y), sizeof(float));
			out.write(reinterpret_cast<const char*>(&box->origin.z), sizeof(float));
		}
		else if (lockedCon->get_type() == ObjectType::CylinderContainerType) {
			containerTypeID = 2;
			out.write(reinterpret_cast<const char*>(&containerTypeID), sizeof(containerTypeID));
			auto* cylinder = static_cast<CylinderContainer*>(lockedCon.get());
			out.write(reinterpret_cast<const char*>(&cylinder->cylinderHeight), sizeof(cylinder->cylinderHeight));
			out.write(reinterpret_cast<const char*>(&cylinder->cylinderRadius), sizeof(cylinder->cylinderRadius));
		}
		else if (lockedCon->get_type() == ObjectType::AbstractContainerType) {
			containerTypeID = 3;
			out.write(reinterpret_cast<const char*>(&containerTypeID), sizeof(containerTypeID));
			auto* meshContainer = static_cast<AbstractContainer*>(lockedCon.get());
			std::string pathStr = meshContainer->fileName;

			uint32_t pathLength = static_cast<uint32_t>(pathStr.size());
			out.write(reinterpret_cast<const char*>(&pathLength), sizeof(pathLength));
			if (pathLength > 0) {
				out.write(pathStr.data(), pathLength);
			}
			float sf = meshContainer->get_scale();
			out.write(reinterpret_cast<const char*>(&sf), sizeof(sf));
		}

		std::string conName = lockedCon->name;
		uint32_t nameLength = static_cast<uint32_t>(conName.size());
		out.write(reinterpret_cast<const char*>(&nameLength), sizeof(nameLength));
		if (nameLength > 0) {
			out.write(conName.data(), nameLength);
		}
	}
	else {
		out.write(reinterpret_cast<const char*>(&containerTypeID), sizeof(containerTypeID));
	}

	// 4. Generator Context Block
	std::shared_ptr<InterfaceSeedGenerator> lockedGen = generator.lock();
	int32_t genTypeID = -1;

	if (lockedGen) {
		if (lockedGen->get_type() == ObjectType::RandomGeneratorType) {
			genTypeID = 0;
			out.write(reinterpret_cast<const char*>(&genTypeID), sizeof(genTypeID));
			auto* randGen = static_cast<Random*>(lockedGen.get());
			int seedSize = randGen->seedNr;
			out.write(reinterpret_cast<const char*>(&seedSize), sizeof(seedSize));
		}
		else if (lockedGen->get_type() == ObjectType::PoissonGeneratorType ||
			lockedGen->get_type() == ObjectType::UniformGeneratorType ||
			lockedGen->get_type() == ObjectType::VariedGeneratorType) {
			genTypeID = 1;
			out.write(reinterpret_cast<const char*>(&genTypeID), sizeof(genTypeID));
			auto* poissonGen = static_cast<Poisson3D*>(lockedGen.get());
			float minRad = static_cast<float>(poissonGen->get_min_radius());
			float maxRad = static_cast<float>(poissonGen->get_max_radius());
			out.write(reinterpret_cast<const char*>(&minRad), sizeof(minRad));
			out.write(reinterpret_cast<const char*>(&maxRad), sizeof(maxRad));
		}

		std::string genName = lockedGen->name;
		uint32_t nameLength = static_cast<uint32_t>(genName.size());
		out.write(reinterpret_cast<const char*>(&nameLength), sizeof(nameLength));
		if (nameLength > 0) {
			out.write(genName.data(), nameLength);
		}
	}
	else {
		out.write(reinterpret_cast<const char*>(&genTypeID), sizeof(genTypeID));
	}

	// 5. Metrics Block
	out.write(reinterpret_cast<const char*>(&porosity), sizeof(porosity));
	out.write(reinterpret_cast<const char*>(&volume), sizeof(volume));
	out.write(reinterpret_cast<const char*>(&surfaceArea), sizeof(surfaceArea));
	out.write(reinterpret_cast<const char*>(&surfaceToVolume), sizeof(surfaceToVolume));
	out.write(reinterpret_cast<const char*>(&connectivityDensity), sizeof(connectivityDensity));
	out.write(reinterpret_cast<const char*>(&localThickness), sizeof(localThickness));
	out.write(reinterpret_cast<const char*>(&localThicknessStd), sizeof(localThicknessStd));
	out.write(reinterpret_cast<const char*>(&localSeparation), sizeof(localSeparation));
	out.write(reinterpret_cast<const char*>(&localSeparationStd), sizeof(localSeparationStd));
	out.write(reinterpret_cast<const char*>(&trabecularNr), sizeof(trabecularNr));
	out.write(reinterpret_cast<const char*>(&anisotropyDegree), sizeof(anisotropyDegree));
	out.write(reinterpret_cast<const char*>(&tortuosity), sizeof(tortuosity));

	// 6. Versions Block
	out.write(reinterpret_cast<const char*>(&thicknessVersion), sizeof(thicknessVersion));
	out.write(reinterpret_cast<const char*>(&separationVersion), sizeof(separationVersion));
	out.write(reinterpret_cast<const char*>(&trabecularNrVersion), sizeof(trabecularNrVersion));
	out.write(reinterpret_cast<const char*>(&connectivityVersion), sizeof(connectivityVersion));
	out.write(reinterpret_cast<const char*>(&tortuosityVersion), sizeof(tortuosityVersion));
	out.write(reinterpret_cast<const char*>(&anisotropyVersion), sizeof(anisotropyVersion));
	out.write(reinterpret_cast<const char*>(&meshVersion), sizeof(meshVersion));

	// 7. Spatial Structure Configuration Block
	out.write(reinterpret_cast<const char*>(blockDims.data()), sizeof(int) * 3);
	out.write(reinterpret_cast<const char*>(bounds.data()), sizeof(float) * 6);

	// 8. Seeds Array Block
	uint64_t seedCount = static_cast<uint64_t>(seeds.size());
	out.write(reinterpret_cast<const char*>(&seedCount), sizeof(seedCount));
	if (seedCount > 0) {
		out.write(reinterpret_cast<const char*>(seeds.data()), sizeof(Vec3) * seedCount);
	}

	// 9. Raw Grid Density Array Block
	uint64_t voxelCount = static_cast<uint64_t>(scalarField.size());
	out.write(reinterpret_cast<const char*>(&voxelCount), sizeof(voxelCount));
	if (voxelCount > 0) {
		out.write(reinterpret_cast<const char*>(scalarField.data()), sizeof(float) * voxelCount);
	}

	// 10. v2 Extension Block: spread, junction smoothing, boundary/edge frame, and
	// the calibration state (checkboxes + targets). Appended at the end and gated on
	// fileVersion so older v1 files still load.
	out.write(reinterpret_cast<const char*>(&spread), sizeof(spread));
	out.write(reinterpret_cast<const char*>(&smoothJunctions), sizeof(smoothJunctions));
	out.write(reinterpret_cast<const char*>(&smoothK), sizeof(smoothK));
	out.write(reinterpret_cast<const char*>(&frameBoundary), sizeof(frameBoundary));
	out.write(reinterpret_cast<const char*>(&frameDepth), sizeof(frameDepth));
	out.write(reinterpret_cast<const char*>(&frameContainerEdges), sizeof(frameContainerEdges));
	out.write(reinterpret_cast<const char*>(&frameBeam), sizeof(frameBeam));
	out.write(reinterpret_cast<const char*>(&calibrateThickness), sizeof(calibrateThickness));
	out.write(reinterpret_cast<const char*>(&targetThickness), sizeof(targetThickness));
	out.write(reinterpret_cast<const char*>(&calibratePorosity), sizeof(calibratePorosity));
	out.write(reinterpret_cast<const char*>(&targetPorosity), sizeof(targetPorosity));
	out.write(reinterpret_cast<const char*>(&calibrateDA), sizeof(calibrateDA));
	out.write(reinterpret_cast<const char*>(&targetDa), sizeof(targetDa));
	out.write(reinterpret_cast<const char*>(&targetFormulaIdx), sizeof(targetFormulaIdx));
	out.write(reinterpret_cast<const char*>(&calibrationStep), sizeof(calibrationStep));
	out.write(reinterpret_cast<const char*>(&calibrationTol), sizeof(calibrationTol));
	out.write(reinterpret_cast<const char*>(&calibrationIter), sizeof(calibrationIter));

	// 11. v3 Extension Block: optional edge rounding (roundEdges, edgeK). Appended
	// after the v2 block and gated on fileVersion >= 3 so v1/v2 files still load.
	out.write(reinterpret_cast<const char*>(&roundEdges), sizeof(roundEdges));
	out.write(reinterpret_cast<const char*>(&edgeK), sizeof(edgeK));

	out.close();
	logger->log(LogPriority::SUCCESS, "Procedural scaffold serialized smoothly to " + fileName);
}

bool GeneratorLewiner::load_scaf(const std::string& fileName,
	std::vector<std::shared_ptr<IContainer>>& containerList,
	std::vector<std::shared_ptr<InterfaceSeedGenerator>>& generatorList,
	std::vector<std::shared_ptr<AnisotropySource>>& globalSources,
	std::atomic<int>* stage
) {
	if (stage) stage->store(0, std::memory_order_relaxed);

	std::ifstream in(fileName, std::ios::in | std::ios::binary);
	if (!in.is_open()) {
		logger->log(LogPriority::ERROR, "Failed to open file for binary deserialization: " + fileName);
		return false;
	}

	char magic[4];
	in.read(magic, 4);
	if (magic[0] != 'S' || magic[1] != 'C' || magic[2] != 'A' || magic[3] != 'F') {
		logger->log(LogPriority::ERROR, "Invalid file format! Not a true .scaf file.");
		return false;
	}

	uint32_t fileVersion = 0;
	in.read(reinterpret_cast<char*>(&fileVersion), sizeof(fileVersion));
	std::cout << fileVersion << std::endl;
	if (fileVersion != 1 && fileVersion != 2 && fileVersion != 3) {
		logger->log(LogPriority::ERROR, "Unsupported file version.");
		return false;
	}

	if (stage) stage->store(1, std::memory_order_relaxed);

	// Read Core Parameters Block
	in.read(reinterpret_cast<char*>(&selectedThicknessOption), sizeof(selectedThicknessOption));
	in.read(reinterpret_cast<char*>(&selectedDist), sizeof(selectedDist));
	in.read(reinterpret_cast<char*>(&selectedFunc), sizeof(selectedFunc));
	in.read(reinterpret_cast<char*>(&startThickness), sizeof(startThickness));
	in.read(reinterpret_cast<char*>(&endThickness), sizeof(endThickness));
	in.read(reinterpret_cast<char*>(&transitionDistance), sizeof(transitionDistance));
	in.read(reinterpret_cast<char*>(&threshold), sizeof(threshold));
	in.read(reinterpret_cast<char*>(&isoLevel), sizeof(isoLevel));
	in.read(reinterpret_cast<char*>(&voxelSize), sizeof(voxelSize));

	// read smoothness
	in.read(reinterpret_cast<char*>(&iter), sizeof(iter));
	in.read(reinterpret_cast<char*>(&lambda), sizeof(lambda));
	in.read(reinterpret_cast<char*>(&mu), sizeof(mu));
	
	in.read(reinterpret_cast<char*>(&stretchX), sizeof(stretchX));
	in.read(reinterpret_cast<char*>(&stretchY), sizeof(stretchY));
	in.read(reinterpret_cast<char*>(&stretchZ), sizeof(stretchZ));
	in.read(reinterpret_cast<char*>(&anisotropyAngle), sizeof(anisotropyAngle));

	// Explicit element parsing for Vec3 elements
	in.read(reinterpret_cast<char*>(&anisotropyVec.x), sizeof(float));
	in.read(reinterpret_cast<char*>(&anisotropyVec.y), sizeof(float));
	in.read(reinterpret_cast<char*>(&anisotropyVec.z), sizeof(float));

	in.read(reinterpret_cast<char*>(&renderMode), sizeof(renderMode));
	in.read(reinterpret_cast<char*>(&backgroundWeight), sizeof(backgroundWeight));

	if (stage) stage->store(2, std::memory_order_relaxed);

	// Load Anisotropy Sources
	anisotropySources.clear();

    uint32_t numSources = 0;
	in.read(reinterpret_cast<char*>(&numSources), sizeof(numSources));

	for (uint32_t i = 0; i < numSources; ++i) {
		auto src = std::make_shared<AnisotropySource>();

		in.read(reinterpret_cast<char*>(&src->origin.x), sizeof(float) * 3);
		in.read(reinterpret_cast<char*>(&src->direction.x), sizeof(float) * 3);
		in.read(reinterpret_cast<char*>(&src->stretch.x), sizeof(float) * 3);
		in.read(reinterpret_cast<char*>(&src->sigma), sizeof(float));

		in.read(reinterpret_cast<char*>(&src->angle), sizeof(float));

		src->name = "Anisotropy Source " + std::to_string(i + 1);
		src->update_metric();
		// update_model() (OpenGL) is deferred to the main thread after load completes

		anisotropySources.push_back(src);
		globalSources.push_back(src);
	}

	if (stage) stage->store(3, std::memory_order_relaxed);

	// Load Container Context Block
	int32_t containerTypeID = 0;
	in.read(reinterpret_cast<char*>(&containerTypeID), sizeof(containerTypeID));
	std::shared_ptr<IContainer> rebuiltContainer = nullptr;

	if (containerTypeID == 1) {
		Vec3 boxSize, boxOrigin;
		in.read(reinterpret_cast<char*>(&boxSize.x), sizeof(float));
		in.read(reinterpret_cast<char*>(&boxSize.y), sizeof(float));
		in.read(reinterpret_cast<char*>(&boxSize.z), sizeof(float));
		in.read(reinterpret_cast<char*>(&boxOrigin.x), sizeof(float));
		in.read(reinterpret_cast<char*>(&boxOrigin.y), sizeof(float));
		in.read(reinterpret_cast<char*>(&boxOrigin.z), sizeof(float));
		rebuiltContainer = std::make_shared<BoxContainer>(boxSize, boxOrigin);
	}
	else if (containerTypeID == 2) {
		float height = 0.0f, radius = 0.0f;
		in.read(reinterpret_cast<char*>(&height), sizeof(height));
		in.read(reinterpret_cast<char*>(&radius), sizeof(radius));
		rebuiltContainer = std::make_shared<CylinderContainer>(height, radius);
	}
	else if (containerTypeID == 3) {
		uint32_t pathLength = 0;
		in.read(reinterpret_cast<char*>(&pathLength), sizeof(pathLength));
		std::string pathStr;
		if (pathLength > 0) {
			pathStr.resize(pathLength);
			in.read(&pathStr[0], pathLength);
		}
		float scaleFactor = 1.0f;
		in.read(reinterpret_cast<char*>(&scaleFactor), sizeof(scaleFactor));
		if (!pathStr.empty() && std::filesystem::exists(pathStr)) {
			rebuiltContainer = std::make_shared<AbstractContainer>(pathStr);
			rebuiltContainer->set_scale(scaleFactor);
		}
	}

	if (containerTypeID != 0) {
		uint32_t nameLength = 0;
		in.read(reinterpret_cast<char*>(&nameLength), sizeof(nameLength));
		std::string conName;
		if (nameLength > 0) {
			conName.resize(nameLength);
			in.read(&conName[0], nameLength);
		}
		if (rebuiltContainer) {
			rebuiltContainer->name = conName.empty() ? "Restored_Container" : conName;
		}
	}

	if (rebuiltContainer) {
		containerList.push_back(rebuiltContainer);
		this->container = rebuiltContainer;
	}

	if (stage) stage->store(4, std::memory_order_relaxed);

	// Read Generator Context Block
	int32_t genTypeID = -1;
	in.read(reinterpret_cast<char*>(&genTypeID), sizeof(genTypeID));

	int tempSeedNr = 0;
	float tempRmin = 0.0f;
	float tempRmax = 0.0f;

	if (genTypeID == 0) {
		in.read(reinterpret_cast<char*>(&tempSeedNr), sizeof(tempSeedNr));
	}
	else if (genTypeID == 1) {
		in.read(reinterpret_cast<char*>(&tempRmin), sizeof(tempRmin));
		in.read(reinterpret_cast<char*>(&tempRmax), sizeof(tempRmax));
	}

	std::string genName = "Restored_Generator";
	if (genTypeID != -1) {
		uint32_t nameLength = 0;
		in.read(reinterpret_cast<char*>(&nameLength), sizeof(nameLength));
		if (nameLength > 0) {
			genName.resize(nameLength);
			in.read(&genName[0], nameLength);
		}
	}

	// --- STREAM HEALTH VERIFICATION ---
	// If our stream collapsed due to layout offsets, trap it before loading garbage metrics!
	if (in.fail()) {
		logger->log(LogPriority::ERROR, "SCAF Stream Deserialization corrupted prior to Metrics parsing block!");
		in.close();
		return false;
	}

	// Read Metrics Block
	in.read(reinterpret_cast<char*>(&porosity), sizeof(porosity));
	in.read(reinterpret_cast<char*>(&volume), sizeof(volume));
	in.read(reinterpret_cast<char*>(&surfaceArea), sizeof(surfaceArea));
	in.read(reinterpret_cast<char*>(&surfaceToVolume), sizeof(surfaceToVolume));
	in.read(reinterpret_cast<char*>(&connectivityDensity), sizeof(connectivityDensity));
	in.read(reinterpret_cast<char*>(&localThickness), sizeof(localThickness));
	in.read(reinterpret_cast<char*>(&localThicknessStd), sizeof(localThicknessStd));
	in.read(reinterpret_cast<char*>(&localSeparation), sizeof(localSeparation));
	in.read(reinterpret_cast<char*>(&localSeparationStd), sizeof(localSeparationStd));
	in.read(reinterpret_cast<char*>(&trabecularNr), sizeof(trabecularNr));
	in.read(reinterpret_cast<char*>(&anisotropyDegree), sizeof(anisotropyDegree));
	in.read(reinterpret_cast<char*>(&tortuosity), sizeof(tortuosity));

	// Read Versions Block
	in.read(reinterpret_cast<char*>(&thicknessVersion), sizeof(thicknessVersion));
	in.read(reinterpret_cast<char*>(&separationVersion), sizeof(separationVersion));
	in.read(reinterpret_cast<char*>(&trabecularNrVersion), sizeof(trabecularNrVersion));
	in.read(reinterpret_cast<char*>(&connectivityVersion), sizeof(connectivityVersion));
	in.read(reinterpret_cast<char*>(&tortuosityVersion), sizeof(tortuosityVersion));
	in.read(reinterpret_cast<char*>(&anisotropyVersion), sizeof(anisotropyVersion));
	in.read(reinterpret_cast<char*>(&meshVersion), sizeof(meshVersion));

	// Map verification states
	bool thicknessValid = (thicknessVersion == meshVersion);
	bool separationValid = (separationVersion == meshVersion);
	bool trabecularValid = (trabecularNrVersion == meshVersion);
	bool connectivityValid = (connectivityVersion == meshVersion);
	bool tortuosityValid = (tortuosityVersion == meshVersion);
	bool anisotropyValid = (anisotropyVersion == meshVersion);

	// Read Grid Bounds Setup Block
	in.read(reinterpret_cast<char*>(blockDims.data()), sizeof(int) * 3);
	in.read(reinterpret_cast<char*>(bounds.data()), sizeof(float) * 6);
	update_steps();

	// Read Seeds Array Block
	uint64_t seedCount = 0;
	in.read(reinterpret_cast<char*>(&seedCount), sizeof(seedCount));
	seeds.resize(seedCount);
	if (seedCount > 0) {
		in.read(reinterpret_cast<char*>(seeds.data()), sizeof(Vec3) * seedCount);
	}

	// Reconstruct Generator Instance Graph
	std::shared_ptr<InterfaceSeedGenerator> rebuiltGen = nullptr;
	if (genTypeID == 0) {
		auto randGen = std::make_shared<Random>();
		randGen->seedNr = tempSeedNr;
		randGen->name = genName;
		rebuiltGen = randGen;
	}
	else if (genTypeID == 1) {
		auto poissonGen = std::make_shared<Poisson3D>();
		poissonGen->name = genName;
		poissonGen->set_min_radius(static_cast<double>(tempRmin));
		poissonGen->set_max_radius(static_cast<double>(tempRmax));
		poissonGen->type = (tempRmin == tempRmax) ?
			ObjectType::UniformGeneratorType : ObjectType::VariedGeneratorType;
		rebuiltGen = poissonGen;
	}

	if (rebuiltGen) {
		rebuiltGen->set_seeds(seeds);
		rebuiltGen->set_renderMode(renderMode);

		// update_model() is an OpenGL call — only safe on the main thread
		if (!stage && this->renderMode && !this->seeds.empty()) {
			rebuiltGen->update_model();
		}

		generatorList.push_back(rebuiltGen);
		this->generator = rebuiltGen;
	}
	// Robustness: a scaffold exported without an attached seed generator (e.g. a
	// batch export built from a raw seed vector) stores genTypeID == -1 and no
	// generator data, which would leave the loaded model uneditable ("generator
	// not found"). Wrap the saved seeds in a fallback uniform generator so the
	// model stays editable — the seeds are preserved exactly; only the spacing
	// radius is unknown (left at its default, so a *re-seed* would use it).
	else if (!seeds.empty()) {
		auto fallbackGen = std::make_shared<Poisson3D>();
		fallbackGen->name = genName;
		fallbackGen->set_seeds(seeds);
		fallbackGen->type = ObjectType::UniformGeneratorType;
		fallbackGen->set_renderMode(renderMode);
		if (!stage && this->renderMode && !this->seeds.empty()) {
			fallbackGen->update_model();
		}
		generatorList.push_back(fallbackGen);
		this->generator = fallbackGen;
		if (logger) logger->log(LogPriority::WARNING,
			"SCAF had no stored seed generator; reconstructed a uniform generator "
			"from the saved seeds (spacing radius unknown).");
	}

	if (stage) stage->store(5, std::memory_order_relaxed);

	// Read Raw Grid Array Block
	uint64_t voxelCount = 0;
	in.read(reinterpret_cast<char*>(&voxelCount), sizeof(voxelCount));
	scalarField.resize(voxelCount);
	if (voxelCount > 0) {
		in.read(reinterpret_cast<char*>(scalarField.data()), sizeof(float) * voxelCount);
	}

	// v2 Extension Block (see export_scaf). Present only in fileVersion >= 2; on v1
	// files the members keep their defaults.
	if (fileVersion >= 2) {
		in.read(reinterpret_cast<char*>(&spread), sizeof(spread));
		in.read(reinterpret_cast<char*>(&smoothJunctions), sizeof(smoothJunctions));
		in.read(reinterpret_cast<char*>(&smoothK), sizeof(smoothK));
		in.read(reinterpret_cast<char*>(&frameBoundary), sizeof(frameBoundary));
		in.read(reinterpret_cast<char*>(&frameDepth), sizeof(frameDepth));
		in.read(reinterpret_cast<char*>(&frameContainerEdges), sizeof(frameContainerEdges));
		in.read(reinterpret_cast<char*>(&frameBeam), sizeof(frameBeam));
		in.read(reinterpret_cast<char*>(&calibrateThickness), sizeof(calibrateThickness));
		in.read(reinterpret_cast<char*>(&targetThickness), sizeof(targetThickness));
		in.read(reinterpret_cast<char*>(&calibratePorosity), sizeof(calibratePorosity));
		in.read(reinterpret_cast<char*>(&targetPorosity), sizeof(targetPorosity));
		in.read(reinterpret_cast<char*>(&calibrateDA), sizeof(calibrateDA));
		in.read(reinterpret_cast<char*>(&targetDa), sizeof(targetDa));
		in.read(reinterpret_cast<char*>(&targetFormulaIdx), sizeof(targetFormulaIdx));
		in.read(reinterpret_cast<char*>(&calibrationStep), sizeof(calibrationStep));
		in.read(reinterpret_cast<char*>(&calibrationTol), sizeof(calibrationTol));
		in.read(reinterpret_cast<char*>(&calibrationIter), sizeof(calibrationIter));
	}

	// v3 Extension Block (see export_scaf). Present only in fileVersion >= 3; on
	// v1/v2 files roundEdges/edgeK keep their defaults (off).
	if (fileVersion >= 3) {
		in.read(reinterpret_cast<char*>(&roundEdges), sizeof(roundEdges));
		in.read(reinterpret_cast<char*>(&edgeK), sizeof(edgeK));
	}

	in.close();

	// Re-verify stream status right before triggering reconstruction
	if (in.fail()) {
		logger->log(LogPriority::ERROR, "SCAF Stream parsing collapsed during binary data array extraction phases!");
		return false;
	}

	if (stage) stage->store(6, std::memory_order_relaxed);

	// Run final level-set boundary tracking reconstruction
	if (!marching_cubes()) {
		if (logger) logger->log(LogPriority::ERROR, "SCAF reconstruction failed: invalid scalar field!");
		return false;
	}

	if (stage) stage->store(7, std::memory_order_relaxed);

	// Re-snap accurate visual state matches post-marching cubes execution
	if (thicknessValid)    thicknessVersion = meshVersion;
	if (separationValid)   separationVersion = meshVersion;
	if (trabecularValid)   trabecularNrVersion = meshVersion;
	if (connectivityValid) connectivityVersion = meshVersion;
	if (tortuosityValid)   tortuosityVersion = meshVersion;
	if (anisotropyValid)   anisotropyVersion = meshVersion;

	isLoadedFromFile = false;

	return true;
}

//@brief function to estimate tortuosity of the porous structure, using the A* algorithm on the grid, we can estimate the shortest path between two points in the porous structure, and compare it to the straight line distance between those points to get an estimate of the tortuosity
bool GeneratorLewiner::estimate_tortuosity(float voxelSize) {

	if (isLoadedFromFile) return false;

    tortuosityPathModel.reset();
    tortuosityPathVertices.clear();
    tortuosityPathEdges.clear();

    // use the bounds of the aabb
    std::array<float, 6> aabbBounds = {
        aabb.pMin.x, aabb.pMax.x,
        aabb.pMin.y, aabb.pMax.y,
        aabb.pMin.z, aabb.pMax.z
    };

    // FIX 1: Pass 'false' (0) instead of 1. 
    // Now Bone = 255 and Air/Pores = 0, which perfectly matches A* logic.
    std::vector<uint8_t> field = get_image_field(voxelSize, aabbBounds, false);

    // estimate new block dimensions
    int nx = static_cast<int>(std::ceil((aabbBounds[1] - aabbBounds[0]) / voxelSize));
    int ny = static_cast<int>(std::ceil((aabbBounds[3] - aabbBounds[2]) / voxelSize));
    int nz = static_cast<int>(std::ceil((aabbBounds[5] - aabbBounds[4]) / voxelSize));

    // Safety check: ensure the field size matches our expected dimensions
    if (field.size() != (size_t)nx * ny * nz) {
        std::cerr << "Dimension mismatch in tortuosity estimation!" << std::endl;
        return false;
    }

    auto get_idx = [&](int x, int y, int z) -> size_t {
        return static_cast<size_t>(x) +
            static_cast<size_t>(y) * static_cast<size_t>(nx) +
            static_cast<size_t>(z) * static_cast<size_t>(nx) * static_cast<size_t>(ny);
    };

    size_t totalVoxels = field.size();
    std::vector<size_t> parentMap(totalVoxels, SIZE_MAX);
    std::vector<float> gScore(totalVoxels, std::numeric_limits<float>::max());
    std::vector<bool> visited(totalVoxels, false);

    // FIX 3: Define a 26-connected Moore Neighborhood for 3D diagonal traversal
    struct Neighbor { int dx, dy, dz; float cost; };
    std::vector<Neighbor> neighbors26;
    for (int dz = -1; dz <= 1; dz++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                if (dx == 0 && dy == 0 && dz == 0) continue;
                // Cost is geometric distance: 1.0 for straight, 1.41 for 2D diag, 1.73 for 3D diag
                float cost = std::sqrt(static_cast<float>(dx * dx + dy * dy + dz * dz));
                neighbors26.push_back({ dx, dy, dz, cost });
            }
        }
    }

    std::priority_queue<AStarNode, std::vector<AStarNode>, std::greater<AStarNode>> openSet;

    int startZ = 1;
    int targetZ = nz - 2;

    int xSize = static_cast<int>(nx * 0.15f);
    int ySize = static_cast<int>(ny * 0.15f);
    if (xSize == 0) xSize = 1;
    if (ySize == 0) ySize = 1;

    bool foundStart = false;

    // Helper lambda to search for inlets
    auto add_inlets = [&](int marginX, int marginY) {
        for (int x = marginX; x < nx - marginX - 1; x++) {
            for (int y = marginY; y < ny - marginY - 1; y++) {
                size_t idx = get_idx(x, y, startZ);
                float h = (targetZ - startZ) * voxelSize;

                // field == 0 is now correctly mapped to Air
                if (field[idx] == 0 && gScore[idx] > 0.0f) { 
                    gScore[idx] = 0.0f;
                    openSet.push({ idx, h });
                    foundStart = true;
                }
            }
        }
    };

    // Try central region first
    add_inlets(xSize, ySize);

    // FIX 2: Fallback to the entire Z slice if the center is blocked by bone
    if (!foundStart) {
        add_inlets(1, 1); 
    }

    if (!foundStart) {
        std::cerr << "Inlet is completely blocked! No starting points found." << std::endl;
        return false;
    }

    float minPathLength = std::numeric_limits<float>::infinity();
    size_t goalIndex = SIZE_MAX;

    while (!openSet.empty()) {
        AStarNode current = openSet.top();
        openSet.pop();

        int idx = current.idx;

        if (visited[idx]) continue;
        visited[idx] = true;

        int x = idx % nx;
        int y = (idx / nx) % ny;
        int z = idx / (nx * ny);

        if (z >= targetZ) {
            minPathLength = gScore[idx];
            goalIndex = idx;
            break;
        }

        // Use the newly defined 26-connected array
        for (const auto& nb : neighbors26) { 
            int nx_ = x + nb.dx;
            int ny_ = y + nb.dy;
            int nz_ = z + nb.dz;

            if (nx_ >= 0 && nx_ < nx && ny_ >= 0 && ny_ < ny && nz_ >= 0 && nz_ < nz) {
                size_t nbIdx = get_idx(nx_, ny_, nz_);

                if (!visited[nbIdx] && field[nbIdx] == 0) {
                    
                    float g = gScore[idx] + (nb.cost * voxelSize);

                    if (g < gScore[nbIdx]) {
                        gScore[nbIdx] = g;

                        // Target exact distance to end plane
                        float h = (targetZ - nz_) * voxelSize; 
                        openSet.push({ nbIdx, g + h});

                        parentMap[nbIdx] = idx;
                    }
                }
            }
        }
    }

    if (goalIndex == SIZE_MAX) {
        std::cerr << "No connected path found from inlet to outlet! Porosity is closed." << std::endl;
        tortuosity = -1;
        return false;
    }

    // --- Path Reconstruction ---
    size_t currIdx = goalIndex;
    int vertexCount = 0;

    while (currIdx != SIZE_MAX) {
        int cx = static_cast<int>(currIdx % nx);
        int cy = static_cast<int>((currIdx / nx) % ny);
        int cz = static_cast<int>(currIdx / (nx * ny));

        Vec3 pos(
            aabbBounds[0] + cx * voxelSize,
            aabbBounds[2] + cy * voxelSize,
            aabbBounds[4] + cz * voxelSize);

        tortuosityPathVertices.push_back(pos.x);
        tortuosityPathVertices.push_back(pos.y);
        tortuosityPathVertices.push_back(pos.z);

        if (vertexCount > 0) { 
            tortuosityPathEdges.push_back(vertexCount - 1);
            tortuosityPathEdges.push_back(vertexCount);
        }

        currIdx = parentMap[currIdx];
        vertexCount++;
    }

    // create a model
    tortuosityPathModel = std::make_unique<PoreNetwork>(tortuosityPathVertices, tortuosityPathEdges);
    
    // tortuosity is the ratio of the actual path length to the straight line distance
    float straightLineDist = (targetZ - startZ) * voxelSize;
    tortuosity = minPathLength / straightLineDist; 
    
    if (tortuosity < 1.0f) {
        tortuosity = 1.0f;
    }

    tortuosityVersion = meshVersion;
	
	logger->log(LogPriority::SUCCESS, "Estimated Tortuosity!");

    return true;
};

void GeneratorLewiner::draw_tortuosity_path() {

	if (tortuosityPathModel) {
		tortuosityPathModel->draw();
	}
};

void GeneratorLewiner::apply_scale() {
	
	int vertNr = vertices.size() / 3;

	for (int i{ 0 }; i < vertNr; i++) {
		vertices[i] *= scaleVec.x;
		vertices[i + 1] *= scaleVec.y;
		vertices[i + 2] *= scaleVec.z;
	}

	glBindVertexArray(VAO); // Optional but good practice to ensure we target correct state

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

}

void GeneratorLewiner::estimate_anisotropy(float voxelSize, int daDirectionNr, int linesPerDirection, int mode, ROI* roi) {

	std::vector<float> milValues(daDirectionNr, 0.0f);

	// Create Random Uniform Directions (Fibonacci Sphere)
	std::vector<Vec3> dirs(daDirectionNr);
	const float PI = 3.14159265359f;
	const float goldenRatio = (1.0f + std::sqrt(5.0f)) * 0.5f;

	for (int i = 0; i < daDirectionNr; i++) {
		float theta = 2.0f * PI * i / goldenRatio;
		float phi = std::acos(1.0f - 2.0f * (i + 0.5f) / daDirectionNr);
		dirs[i] = Vec3(std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta), std::cos(phi));
	}

	// Resample the structure to a binary image at the requested voxel size, the
	// same substrate used by local thickness/separation. The MIL rays then walk
	// this image, so every voxel-based metric shares one sampling resolution.
	std::array<float, 6> blockBounds = { aabb.pMin.x, aabb.pMax.x, aabb.pMin.y, aabb.pMax.y, aabb.pMin.z, aabb.pMax.z };
	if (roi) {
		std::array<float, 6> reqBounds = roi->get_bounds();
		blockBounds[0] = std::max(reqBounds[0], aabb.pMin.x);
		blockBounds[1] = std::min(reqBounds[1], aabb.pMax.x);
		blockBounds[2] = std::max(reqBounds[2], aabb.pMin.y);
		blockBounds[3] = std::min(reqBounds[3], aabb.pMax.y);
		blockBounds[4] = std::max(reqBounds[4], aabb.pMin.z);
		blockBounds[5] = std::min(reqBounds[5], aabb.pMax.z);
	}

	std::vector<uint8_t> field = get_image_field(voxelSize, blockBounds, false);
	const long nx = static_cast<long>(std::ceil((blockBounds[1] - blockBounds[0]) / voxelSize));
	const long ny = static_cast<long>(std::ceil((blockBounds[3] - blockBounds[2]) / voxelSize));
	const long nz = static_cast<long>(std::ceil((blockBounds[5] - blockBounds[4]) / voxelSize));

	auto is_solid = [&](long vx, long vy, long vz) -> bool {
		return field[static_cast<size_t>(vx) + static_cast<size_t>(vy) * nx +
			static_cast<size_t>(vz) * nx * ny] != 0;
	};

	// Container mask on the resampled grid. For a non-box container (e.g. a
	// cylinder) the AABB corners lie OUTSIDE the wall; their structure-free ray
	// length would bias each direction's MIL DIFFERENTLY (axial rays skirt the
	// corners, radial rays plough through them), distorting the fitted ellipsoid
	// and hence DA. Restrict the intercept length and crossing count to voxels
	// inside the container - matching the Tb.N MIL path and local thickness.
	// Because DA is a RATIO of per-direction MIL (normalized by maxMIL below), a
	// box container - whose mask is all-inside - is left essentially unchanged.
	auto conMask = container.lock();
	std::vector<uint8_t> insideField;
	if (conMask) {
		insideField.assign(field.size(), 1);
		#pragma omp parallel for collapse(3)
		for (long z = 0; z < nz; z++) {
			for (long y = 0; y < ny; y++) {
				for (long x = 0; x < nx; x++) {
					float px = blockBounds[0] + (x + 0.5f) * voxelSize;
					float py = blockBounds[2] + (y + 0.5f) * voxelSize;
					float pz = blockBounds[4] + (z + 0.5f) * voxelSize;
					insideField[static_cast<size_t>(x) + static_cast<size_t>(y) * nx +
						static_cast<size_t>(z) * nx * ny] =
						conMask->is_inside(Vec3(px, py, pz)) ? 1 : 0;
				}
			}
		}
	}
	auto is_inside_vox = [&](long vx, long vy, long vz) -> bool {
		if (insideField.empty()) return true; // no container bound: measure the whole box
		return insideField[static_cast<size_t>(vx) + static_cast<size_t>(vy) * nx +
			static_cast<size_t>(vz) * nx * ny] != 0;
	};

	// ray box spans the whole resampled image (image-voxel coordinates)
	float boundsVoxel[6] = {
		0.0f, static_cast<float>(nx),
		0.0f, static_cast<float>(ny),
		0.0f, static_cast<float>(nz)
	};

	Vec3 boxCenter((boundsVoxel[1] + boundsVoxel[0]) * 0.5f,
		(boundsVoxel[3] + boundsVoxel[2]) * 0.5f,
		(boundsVoxel[5] + boundsVoxel[4]) * 0.5f);

	float d_plane = std::sqrt(std::pow(boundsVoxel[1] - boundsVoxel[0], 2) + std::pow(boundsVoxel[3] - boundsVoxel[2], 2) + std::pow(boundsVoxel[5] - boundsVoxel[4], 2));

	float R = d_plane * 0.5f;

	int gridN = static_cast<int>(std::ceil(std::sqrt(linesPerDirection)));
	float gridStep = d_plane / std::max(1, gridN);

	//BoneJ's exact sampling increment
	float rayStepSize = std::sqrt(3.0f);

#pragma omp parallel 
	{
		std::mt19937 rng(1337 + omp_get_thread_num());
		std::uniform_real_distribution<float> dist(0.0f, 1.0f);

#pragma omp for
		for (int i = 0; i < daDirectionNr; i++) {
			Vec3 d = dirs[i];
			Vec3 w = (std::abs(d.x) > 0.9f) ? Vec3(0, 1, 0) : Vec3(1, 0, 0);
			Vec3 u = d.cross(w).normalized();
			Vec3 v = d.cross(u).normalized();

			long long localTransitions = 0;
			float localBoxLen = 0.0;

			for (int uIdx = 0; uIdx < gridN; uIdx++) {
				for (int vIdx = 0; vIdx < gridN; vIdx++) {
					float uPos = -R + (uIdx + dist(rng)) * gridStep;
					float vPos = -R + (vIdx + dist(rng)) * gridStep;
					Vec3 rayOrigin = boxCenter + (u * uPos) + (v * vPos) - (d * R);

					float tMin = 0.0f;
					float tMax = 1e9f;
					bool hit = true;

					// The loop runs identically, but uses the dynamically updated boundsVoxel array!
					for (int axis = 0; axis < 3; ++axis) {
						float invD = 1.0f / (axis == 0 ? d.x : (axis == 1 ? d.y : d.z));
						float t0 = (boundsVoxel[axis * 2] - (axis == 0 ? rayOrigin.x : (axis == 1 ? rayOrigin.y : rayOrigin.z))) * invD;
						float t1 = (boundsVoxel[axis * 2 + 1] - (axis == 0 ? rayOrigin.x : (axis == 1 ? rayOrigin.y : rayOrigin.z))) * invD;
						if (invD < 0.0f) std::swap(t0, t1);
						tMin = std::max(tMin, t0);
						tMax = std::min(tMax, t1);
						if (tMax <= tMin) { hit = false; break; }
					}

					if (hit && tMax > 0.0f) {
						tMin = std::max(0.0f, tMin);
						float startT = tMin + dist(rng) * rayStepSize;
						long samples = static_cast<long>(std::ceil((tMax - startT) / rayStepSize));
						bool previousPhase = false;
						bool previousInside = false;

						for (long s = 0; s < samples; s++) {
							Vec3 pt = rayOrigin + d * (startT + s * rayStepSize);
							long vx = std::clamp<long>(static_cast<long>(pt.x), static_cast<long>(boundsVoxel[0]), static_cast<long>(boundsVoxel[1]) - 1);
							long vy = std::clamp<long>(static_cast<long>(pt.y), static_cast<long>(boundsVoxel[2]), static_cast<long>(boundsVoxel[3]) - 1);
							long vz = std::clamp<long>(static_cast<long>(pt.z), static_cast<long>(boundsVoxel[4]), static_cast<long>(boundsVoxel[5]) - 1);

							bool currentInside = is_inside_vox(vx, vy, vz);
							bool currentPhase = is_solid(vx, vy, vz);

							// Count length and crossings only inside the container, so the
							// empty AABB corners (and the wall itself) never enter the MIL.
							if (currentInside) {
								localBoxLen += rayStepSize;
								if (previousInside && currentPhase != previousPhase) { localTransitions++; }
							}
							previousPhase = currentPhase;
							previousInside = currentInside;
						}
					}
				}
			}
			milValues[i] = (localBoxLen > 0.0) ? ((localTransitions > 0) ? static_cast<float>(localBoxLen / localTransitions) : static_cast<float>(localBoxLen)) : d_plane;
		}
	}

	// get the max mean interception value we found so far
	float maxMIL = 1e-9f;
	for (float m : milValues) {
		if (m > maxMIL) maxMIL = m;
	}

	// quadratic fit
	Eigen::MatrixXf A_mat(daDirectionNr, 9);
	Eigen::VectorXf b_vec(daDirectionNr);

	for (int v = 0; v < daDirectionNr; v++) {
		float mil = milValues[v] / maxMIL;

		float px = dirs[v].x * mil;
		float py = dirs[v].y * mil;
		float pz = dirs[v].z * mil;

		A_mat(v, 0) = px * px;
		A_mat(v, 1) = py * py;
		A_mat(v, 2) = pz * pz;
		A_mat(v, 3) = px * py;
		A_mat(v, 4) = px * pz;
		A_mat(v, 5) = py * pz;
		A_mat(v, 6) = px;
		A_mat(v, 7) = py;
		A_mat(v, 8) = pz;

		b_vec(v) = 1.0;
	}

	Eigen::VectorXf beta = A_mat.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_vec);

	float a = beta(0), b = beta(1), c = beta(2);
	float d_val = beta(3) / 2.0;
	float e_val = beta(4) / 2.0;
	float f_val = beta(5) / 2.0;
	float g = beta(6) / 2.0;
	float h = beta(7) / 2.0;
	float i_val = beta(8) / 2.0;

	Eigen::Matrix4f quadric;
	quadric << a, d_val, e_val, g,
		d_val, b, f_val, h,
		e_val, f_val, c, i_val,
		g, h, i_val, -1.0;

	// Find Center
	Eigen::Matrix3f sub;
	sub << a, d_val, e_val,
		d_val, b, f_val,
		e_val, f_val, c;
	sub = -1.0 * sub;
	Eigen::Vector3f translationVec(g, h, i_val);
	Eigen::Vector3f center = sub.inverse() * translationVec;

	// Translate to origin
	Eigen::Matrix4f tMat = Eigen::Matrix4f::Identity();
	tMat(0, 3) = center.x();
	tMat(1, 3) = center.y();
	tMat(2, 3) = center.z();

	Eigen::Matrix4f translated = tMat * quadric * tMat.transpose();

	// Eigendecomposition 
	Eigen::Matrix3f input;
	float scale = -1.0 / translated(3, 3);
	for (int row = 0; row < 3; ++row) {
		for (int col = 0; col < 3; ++col) {
			input(row, col) = translated(row, col) * scale;
		}
	}

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> solver(input);
	Eigen::Vector3f evals = solver.eigenvalues();

	if (evals(0) <= 0.0 || evals(1) <= 0.0 || evals(2) <= 0.0) {
		anisotropyDegree = 0.0f;
	}
	else {
		float lambdaMin = evals(0);
		float lambdaMax = evals(2);

		switch (mode) {
		case 0: {
			// ratio of ellipsoid radii (Lmax/Lmin) = sqrt(lambdaMax/lambdaMin), >= 1
			float lenMax = (lambdaMin > 1e-9f) ? (1.0f / std::sqrt(lambdaMin)) : 0.0f;
			float lenMin = (lambdaMax > 1e-9f) ? (1.0f / std::sqrt(lambdaMax)) : 0.0f;
			anisotropyDegree = (lenMax > 1e-9f) ? (lenMax / lenMin) : 0.0f;
			break;
		}
		// mode 1 (1 - Lmax/Lmin) removed: Lmax/Lmin >= 1 so it is always <= 0,
		// which is not a meaningful degree of anisotropy.
		case 2: {
			// eigenvalue ratio lambdaMin/lambdaMax in (0,1], 1 = isotropic
			anisotropyDegree = (lambdaMax > 1e-9f) ? static_cast<float>(lambdaMin / lambdaMax) : 0.0f;
			break;
		}
		case 3:
		default: {
			// 1 - lambdaMin/lambdaMax in [0,1), 0 = isotropic (BoneJ convention)
			anisotropyDegree = (lambdaMax > 1e-9f) ? static_cast<float>(1.0 - (lambdaMin / lambdaMax)) : 0.0f;
			break;
		}
		}
	}

	// create also the ellipsoid pointcloud
	float r1 = (evals(0) > 1e-9) ? maxMIL / std::sqrt(evals(0)) : 0.0f;
	float r2 = (evals(1) > 1e-9) ? maxMIL / std::sqrt(evals(1)) : 0.0f;
	float r3 = (evals(2) > 1e-9) ? maxMIL / std::sqrt(evals(2)) : 0.0f;
	Eigen::Matrix3f v = solver.eigenvectors();

	Eigen::Vector3f origin(
		(0.5f * bounds[1] - bounds[0]),
		(0.5f * bounds[3] - bounds[2]),
		(0.5f * bounds[5] - bounds[4]));

	//ellipsoidModel = std::make_unique<Ellipsoid>(origin, v, r1, r2, r3);
	anisotropyVersion = meshVersion;

	logger->log(LogPriority::SUCCESS, "Estimated Anisotropy Degree!");
};

// Requires: #include <limits>  (for std::numeric_limits<float>::infinity())
// Also assumes <random>, <cmath>, <algorithm>, <vector>, and OpenMP are already included.
//
// Key changes vs. the original:
//   1) Tb.N = P_L / 2: the DDA counts every interface crossing, and each
//      trabecula gives two (entry + exit), so the trabecular count per length
//      is half the crossing density. (An earlier version dropped this /2 and
//      read ~2x high vs. reference values.)
//   2) Replaced the sqrt(3) point-sampling with an exact 3D DDA (Amanatides & Woo)
//      voxel walk, so no thin strut is ever skipped and the count is step-size-free.
//   3) Baseline phase is taken from the FIRST voxel, removing the phantom entry transition.
//   4) Rays walk a binary image resampled at the caller's voxelSize, so the
//      1/voxel -> 1/mm conversion uses that voxelSize (mm/voxel) directly.

void GeneratorLewiner::estimate_trabecular_number(
	float voxelSize, 
	int formula, 
	int daDirectionNr, int linesPerDirection, ROI* roi) {

	// this is if the formula is set to use the Tb.N = (BV/TV) / Tb.Th	
	if (
		(formula == 0 || formula == 1) && trabecularNrVersion < thicknessVersion){
		logger->log(LogPriority::WARNING, "The Estimated Trabecular Number is measured with the previously estimated local thickness. If this was estimated inside a ROI, the estimation is wrong. Try creating the scaffold inside the ROI first.");
	}

	if (formula == 2 && trabecularNrVersion < separationVersion){
		logger->log(LogPriority::WARNING, "The Estimated Trabecular Number is measured with the previously estimated local separation. If this was estimated inside a ROI, the estimation is wrong. Try creating the scaffold inside the ROI first.");
	} 

	if (formula == 1) {
		float bvtv = 1.0f - porosity;
		trabecularNr = bvtv / localThickness;
		logger->log(LogPriority::SUCCESS, "Estimated Trabecular Number with Tb.N = (BV/TV) / Tb.Th");
		trabecularNrVersion = meshVersion;
		return;
	}

	if (formula == 2) {
		trabecularNr = 1.0f / localSeparation;
		logger->log(LogPriority::SUCCESS, "Estimated Trabecular Number with Tb.N = 1/ Tb.Sp");
		trabecularNrVersion = meshVersion;
		return;
	}

	// else use the MIL method
	std::vector<Vec3> dirs(daDirectionNr);
	const float PI = 3.14159265359f;
	const float goldenRatio = (1.0f + std::sqrt(5.0f)) * 0.5f;
	for (int i = 0; i < daDirectionNr; i++) {
		float theta = 2.0f * PI * i / goldenRatio;
		float phi = std::acos(1.0f - 2.0f * (i + 0.5f) / daDirectionNr);
		dirs[i] = Vec3(std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta), std::cos(phi));
	}

	// Resample the structure to a binary image at the requested voxel size (the
	// same substrate as local thickness/separation) and cast the MIL rays over
	// it, so the trabecular number is measured at the target resolution and the
	// physical conversion below uses that same voxel size.
	std::array<float, 6> blockBounds = { aabb.pMin.x, aabb.pMax.x, aabb.pMin.y, aabb.pMax.y, aabb.pMin.z, aabb.pMax.z };
	if (roi) {
		std::array<float, 6> reqBounds = roi->get_bounds();
		blockBounds[0] = std::max(reqBounds[0], aabb.pMin.x);
		blockBounds[1] = std::min(reqBounds[1], aabb.pMax.x);
		blockBounds[2] = std::max(reqBounds[2], aabb.pMin.y);
		blockBounds[3] = std::min(reqBounds[3], aabb.pMax.y);
		blockBounds[4] = std::max(reqBounds[4], aabb.pMin.z);
		blockBounds[5] = std::min(reqBounds[5], aabb.pMax.z);
	}

	std::vector<uint8_t> field = get_image_field(voxelSize, blockBounds, false);
	const long nx = static_cast<long>(std::ceil((blockBounds[1] - blockBounds[0]) / voxelSize));
	const long ny = static_cast<long>(std::ceil((blockBounds[3] - blockBounds[2]) / voxelSize));
	const long nz = static_cast<long>(std::ceil((blockBounds[5] - blockBounds[4]) / voxelSize));

	auto is_solid = [&](long vx, long vy, long vz) -> bool {
		return field[static_cast<size_t>(vx) + static_cast<size_t>(vy) * nx +
			static_cast<size_t>(vz) * nx * ny] != 0;
	};

	// Container mask on the SAME resampled grid. For a non-box container (e.g. a
	// cylinder) the resampled AABB has empty corners OUTSIDE the container wall.
	// Those corners are structure-free, so counting their ray length in the MIL
	// denominator (P_L = crossings / length) dilutes P_L and biases Tb.N LOW - for
	// a cylinder inscribed in its box the box is ~4/pi larger, a ~-21% error.
	// Restrict both the intercept length and the crossing count to voxels inside
	// the container, exactly as local thickness / separation / porosity already do.
	// A box container fills its AABB, so the mask is all-inside and Tb.N is
	// unchanged for the box case.
	auto conMask = container.lock();
	std::vector<uint8_t> insideField;
	if (conMask) {
		insideField.assign(field.size(), 1);
		#pragma omp parallel for collapse(3)
		for (long z = 0; z < nz; z++) {
			for (long y = 0; y < ny; y++) {
				for (long x = 0; x < nx; x++) {
					float px = blockBounds[0] + (x + 0.5f) * voxelSize;
					float py = blockBounds[2] + (y + 0.5f) * voxelSize;
					float pz = blockBounds[4] + (z + 0.5f) * voxelSize;
					insideField[static_cast<size_t>(x) + static_cast<size_t>(y) * nx +
						static_cast<size_t>(z) * nx * ny] =
						conMask->is_inside(Vec3(px, py, pz)) ? 1 : 0;
				}
			}
		}
	}
	auto is_inside_vox = [&](long vx, long vy, long vz) -> bool {
		if (insideField.empty()) return true; // no container bound: measure the whole box
		return insideField[static_cast<size_t>(vx) + static_cast<size_t>(vy) * nx +
			static_cast<size_t>(vz) * nx * ny] != 0;
	};

	float boundsVoxel[6] = {
		0.0f, static_cast<float>(nx),
		0.0f, static_cast<float>(ny),
		0.0f, static_cast<float>(nz)
	};

	// Center plane tracks the specific targeted box sub-window
	Vec3 boxCenter((boundsVoxel[1] + boundsVoxel[0]) * 0.5f,
		(boundsVoxel[3] + boundsVoxel[2]) * 0.5f,
		(boundsVoxel[5] + boundsVoxel[4]) * 0.5f);

	// create the plane that passes through the center
	float d_plane = std::sqrt(
		std::pow(boundsVoxel[1] - boundsVoxel[0], 2) +
		std::pow(boundsVoxel[3] - boundsVoxel[2], 2) +
		std::pow(boundsVoxel[5] - boundsVoxel[4], 2)
	);

	float R = d_plane * 0.5f;
	int gridN = static_cast<int>(std::ceil(std::sqrt(static_cast<float>(linesPerDirection))));
	float gridStep = d_plane / std::max(1, gridN);

	const float INF = std::numeric_limits<const float>::infinity();

	// Global accumulators
	long long globalTransitions = 0;
	double globalBoxLen = 0.0;

	// 3. Parallelize Ray Marching
#pragma omp parallel reduction(+:globalTransitions, globalBoxLen)
	{
		std::mt19937 rng(1337 + omp_get_thread_num());
		std::uniform_real_distribution<float> dist(0.0f, 1.0f);
#pragma omp for
		for (int i = 0; i < daDirectionNr; i++) {
			Vec3 d = dirs[i];
			Vec3 w = (std::abs(d.x) > 0.9f) ? Vec3(0, 1, 0) : Vec3(1, 0, 0);
			Vec3 u = d.cross(w).normalized();
			Vec3 v = d.cross(u).normalized();

			// Local tally for THIS specific direction 'i'
			long long localTransitions = 0;
			double localBoxLen = 0.0;

			for (int uIdx = 0; uIdx < gridN; uIdx++) {
				for (int vIdx = 0; vIdx < gridN; vIdx++) {
					float uPos = -R + (uIdx + dist(rng)) * gridStep;
					float vPos = -R + (vIdx + dist(rng)) * gridStep;
					Vec3 rayOrigin = boxCenter + (u * uPos) + (v * vPos) - (d * R);

					// --- Slab test: clip ray to the voxel box (same as before) ---
					float tMin = 0.0f;
					float tMax = 1e9f;
					bool hit = true;
					for (int axis = 0; axis < 3; ++axis) {
						float dd = (axis == 0 ? d.x : (axis == 1 ? d.y : d.z));
						float oo = (axis == 0 ? rayOrigin.x : (axis == 1 ? rayOrigin.y : rayOrigin.z));
						float invD = 1.0f / dd;
						float t0 = (boundsVoxel[axis * 2] - oo) * invD;
						float t1 = (boundsVoxel[axis * 2 + 1] - oo) * invD;
						if (invD < 0.0f) std::swap(t0, t1);
						tMin = std::max(tMin, t0);
						tMax = std::min(tMax, t1);
						if (tMax <= tMin) { hit = false; break; }
					}
					if (!hit || tMax <= 0.0f) continue;

					float tEntry = std::max(0.0f, tMin);
					float tExit = tMax;

					// --- Exact 3D DDA Voxel Walk within the Subgrid ---
					Vec3 p = rayOrigin + d * tEntry;
					int ix = std::clamp<int>(static_cast<int>(std::floor(p.x)), static_cast<int>(boundsVoxel[0]), static_cast<int>(boundsVoxel[1]) - 1);
					int iy = std::clamp<int>(static_cast<int>(std::floor(p.y)), static_cast<int>(boundsVoxel[2]), static_cast<int>(boundsVoxel[3]) - 1);
					int iz = std::clamp<int>(static_cast<int>(std::floor(p.z)), static_cast<int>(boundsVoxel[4]), static_cast<int>(boundsVoxel[5]) - 1);

					int sx = (d.x > 0.0f) ? 1 : ((d.x < 0.0f) ? -1 : 0);
					int sy = (d.y > 0.0f) ? 1 : ((d.y < 0.0f) ? -1 : 0);
					int sz = (d.z > 0.0f) ? 1 : ((d.z < 0.0f) ? -1 : 0);

					float tDeltaX = (d.x != 0.0f) ? std::abs(1.0f / d.x) : INF;
					float tDeltaY = (d.y != 0.0f) ? std::abs(1.0f / d.y) : INF;
					float tDeltaZ = (d.z != 0.0f) ? std::abs(1.0f / d.z) : INF;

					float nextX = (sx > 0) ? static_cast<float>(ix + 1) : static_cast<float>(ix);
					float nextY = (sy > 0) ? static_cast<float>(iy + 1) : static_cast<float>(iy);
					float nextZ = (sz > 0) ? static_cast<float>(iz + 1) : static_cast<float>(iz);
					float tMaxX = (d.x != 0.0f) ? (tEntry + (nextX - p.x) / d.x) : INF;
					float tMaxY = (d.y != 0.0f) ? (tEntry + (nextY - p.y) / d.y) : INF;
					float tMaxZ = (d.z != 0.0f) ? (tEntry + (nextZ - p.z) / d.z) : INF;

					float tCur = tEntry;
					bool previousPhase = is_solid(ix, iy, iz);
					bool previousInside = is_inside_vox(ix, iy, iz);

					while (true) {
						float tNext = std::min(tMaxX, std::min(tMaxY, tMaxZ));
						float segEnd = std::min(tNext, tExit);

						// Length of the ray inside the CURRENT voxel; only count it
						// toward the MIL denominator when that voxel lies inside the
						// container (empty AABB corners are excluded).
						if (previousInside && segEnd > tCur) localBoxLen += (segEnd - tCur);

						if (tNext > tExit) break;
						tCur = tNext;

						if (tMaxX <= tMaxY && tMaxX <= tMaxZ) { ix += sx; tMaxX += tDeltaX; }
						else if (tMaxY <= tMaxZ) { iy += sy; tMaxY += tDeltaY; }
						else { iz += sz; tMaxZ += tDeltaZ; }

						// Break if out of the scoped boundary limits
						if (ix < static_cast<int>(boundsVoxel[0]) || ix >= static_cast<int>(boundsVoxel[1]) ||
							iy < static_cast<int>(boundsVoxel[2]) || iy >= static_cast<int>(boundsVoxel[3]) ||
							iz < static_cast<int>(boundsVoxel[4]) || iz >= static_cast<int>(boundsVoxel[5])) break;

						bool currentPhase = is_solid(ix, iy, iz);
						bool currentInside = is_inside_vox(ix, iy, iz);
						// Count a solid<->pore interface only when BOTH the entered and
						// the left voxel are inside the container, so the container wall
						// (interior solid -> exterior air) is never miscounted as a
						// trabecular crossing.
						if (currentInside && previousInside && currentPhase != previousPhase) localTransitions++;
						previousPhase = currentPhase;
						previousInside = currentInside;
					}
				}
			}

			globalTransitions += localTransitions;
			globalBoxLen += localBoxLen;
		}
	}

	// final values
	if (globalBoxLen > 0.0) {
		double PL_voxels = static_cast<double>(globalTransitions) / globalBoxLen; // interface crossings / voxel-length
		// Each trabecula crossed by a test line produces TWO interface crossings
		// (entry + exit), so the trabecular count per length is P_L / 2.
		double TbN_voxels = 0.5 * PL_voxels;
		trabecularNr = static_cast<float>(TbN_voxels / voxelSize); // 1/voxel -> 1/mm  (voxelSize = mm/voxel)
	}
	else {
		trabecularNr = 0.0f;
	}
	trabecularNrVersion = meshVersion;

	logger->log(LogPriority::SUCCESS, "Estimated Trabecular Number using MIL method!");
}

void GeneratorLewiner::estimate_connectivity_density() {

	// ensure that edgeset is actually added
	if (edgeSet.empty()) {
		// use the triangles to add edges
		for (const auto& tri : meshTriangles) {
			add_edge(tri.v1, tri.v2, tri.v3);
		}
	}

	// we can also estimate the connectivity density since we have also the mesh
	long long V = meshVertices.size();
	long long F = meshTriangles.size();
	long long E = edgeSet.size();
	long long eulerCharacteristic = V - E + F;

	// estimate the genus
	float genus = 1.0f - (static_cast<float>(eulerCharacteristic) / 2.0f);

	// connectivity density is genus / domain volume
	connectivityDensity = genus / domainVolume;
	
	connectivityVersion = meshVersion;

	logger->log(LogPriority::SUCCESS, "Estimated Connectivity Density!");
};

float GeneratorLewiner::estimate_connectivity_density_voxel(float voxelSize, int connectivity) {

	if (connectivity != 6 && logger) {
		logger->log(LogPriority::WARNING,
			"Voxel connectivity: only 6-connectivity is implemented; using 6.");
	}

	// Resample the solid to a binary image at the target voxel size, the same
	// substrate as thickness / separation / Tb.N / DA.
	std::array<float, 6> blockBounds = { aabb.pMin.x, aabb.pMax.x, aabb.pMin.y, aabb.pMax.y, aabb.pMin.z, aabb.pMax.z };
	std::vector<uint8_t> field = get_image_field(voxelSize, blockBounds, false);
	const long nx = static_cast<long>(std::ceil((blockBounds[1] - blockBounds[0]) / voxelSize));
	const long ny = static_cast<long>(std::ceil((blockBounds[3] - blockBounds[2]) / voxelSize));
	const long nz = static_cast<long>(std::ceil((blockBounds[5] - blockBounds[4]) / voxelSize));

	// foreground accessor; out-of-range = background, which closes the object at
	// the image boundary (implicit one-voxel padding).
	auto fg = [&](long i, long j, long k) -> int {
		if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) return 0;
		return field[static_cast<size_t>(i) + static_cast<size_t>(j) * nx +
			static_cast<size_t>(k) * nx * ny] != 0 ? 1 : 0;
	};

	// Precompute the per-window Euler contribution for all 256 configurations of
	// a 2x2x2 voxel block (bit index = dx + 2*dy + 4*dz). Each cubical cell is
	// attributed to the window whose max corner is the lattice point:
	//   dchi = vertex - (3 edges) + (3 faces) - cube   [chi = V - E + F - C]
	int lut[256];
	for (int code = 0; code < 256; ++code) {
		auto bit = [&](int dx, int dy, int dz) { return (code >> (dx + 2 * dy + 4 * dz)) & 1; };
		int v000 = bit(0, 0, 0), v100 = bit(1, 0, 0), v010 = bit(0, 1, 0), v110 = bit(1, 1, 0);
		int v001 = bit(0, 0, 1), v101 = bit(1, 0, 1), v011 = bit(0, 1, 1), v111 = bit(1, 1, 1);

		int vertex = (code != 0) ? 1 : 0;                       // any of the 8 voxels
		int edgeX = (v000 | v010 | v001 | v011) ? 1 : 0;        // edge +x: voxels v(0,dy,dz)
		int edgeY = (v000 | v100 | v001 | v101) ? 1 : 0;        // edge +y: voxels v(dx,0,dz)
		int edgeZ = (v000 | v100 | v010 | v110) ? 1 : 0;        // edge +z: voxels v(dx,dy,0)
		int faceX = (v000 | v100) ? 1 : 0;                      // face perp x
		int faceY = (v000 | v010) ? 1 : 0;                      // face perp y
		int faceZ = (v000 | v001) ? 1 : 0;                      // face perp z
		int cube = v000;

		lut[code] = vertex - (edgeX + edgeY + edgeZ) + (faceX + faceY + faceZ) - cube;
	}

	// Sum the contribution over every lattice window (px,py,pz) in [0,n]. The
	// window at (px,py,pz) holds the 8 voxels (px-1..px, py-1..py, pz-1..pz).
	long long chi = 0;
	#pragma omp parallel for reduction(+:chi)
	for (long pz = 0; pz <= nz; ++pz) {
		for (long py = 0; py <= ny; ++py) {
			for (long px = 0; px <= nx; ++px) {
				int code = 0;
				for (int dz = 0; dz < 2; ++dz)
					for (int dy = 0; dy < 2; ++dy)
						for (int dx = 0; dx < 2; ++dx)
							if (fg(px - 1 + dx, py - 1 + dy, pz - 1 + dz))
								code |= 1 << (dx + 2 * dy + 4 * dz);
				chi += lut[code];
			}
		}
	}

	// Connectivity = 1 - chi (= beta_1, the number of independent loops, for a
	// single connected cavity-free object). Density normalizes by the container
	// volume, matching the reported units of 1/mm^3.
	double connectivity_beta1 = 1.0 - static_cast<double>(chi);
	connectivityDensity = (domainVolume > 0.0f)
		? static_cast<float>(connectivity_beta1 / domainVolume)
		: 0.0f;

	connectivityVersion = meshVersion;

	if (logger) {
		std::ostringstream oss;
		oss << "Estimated voxel Connectivity Density (Euler=" << chi
			<< ", Conn=" << connectivity_beta1 << "): " << connectivityDensity << " 1/mm^3";
		logger->log(LogPriority::SUCCESS, oss.str());
	}

	return connectivityDensity;
};

void GeneratorLewiner::estimate_connectivity_network() {

	logger->log(LogPriority::ERROR, "Not implemented yet!");
	return;

};

void GeneratorLewiner::apply_taubin_smooth(int iter, float lambda, float mu) {

	size_t vertNr = meshVertices.size();

	std::vector<Vec3> currentVerts(vertNr);
	std::vector<Vec3> tempVerts(vertNr);

	for (int i = 0; i < vertNr; ++i) {
		currentVerts[i] = Vec3(meshVertices[i].x, meshVertices[i].y, meshVertices[i].z);
	}

	// 3. Taubin Smoothing Loop
	for (int k = 0; k < iter; ++k) {

		// Pass 1: Shrink (using lambda > 0)
		#pragma omp parallel for
		for (int i = 0; i < vertNr; ++i) {
			const auto& nbrs = adjacency[i];
			if (nbrs.empty()) {
				tempVerts[i] = currentVerts[i];
				continue;
			}

			Vec3 avg(0.0f, 0.0f, 0.0f);
			for (int n : nbrs) {
				avg = avg + currentVerts[n];
			}
			float inv = 1.0f / static_cast<float>(nbrs.size());
			avg = avg * inv;

			tempVerts[i] = currentVerts[i] + (avg - currentVerts[i]) * lambda;
		}

		// Pass 2: Inflate (using mu < 0)
		#pragma omp parallel for
		for (int i = 0; i < vertNr; ++i) {
			const auto& nbrs = adjacency[i];
			if (nbrs.empty()) {
				currentVerts[i] = tempVerts[i];
				continue;
			}

			Vec3 avg(0.0f, 0.0f, 0.0f);
			for (int n : nbrs) {
				avg = avg + tempVerts[n];
			}
			float inv = 1.0f / static_cast<float>(nbrs.size());
			avg = avg * inv;

			currentVerts[i] = tempVerts[i] + (avg - tempVerts[i]) * mu;
		}
	}

	// 4. Write back to the mesh
	for (int i = 0; i < vertNr; ++i) {
		meshVertices[i].x = currentVerts[i].x;
		meshVertices[i].y = currentVerts[i].y;
		meshVertices[i].z = currentVerts[i].z;
	}

	// build the model again
	_update_bounding_box();

	// update opengl objects
	update_render();

	logger->log(LogPriority::SUCCESS, "Applied Taubin Smoothing!");

};

void GeneratorLewiner::apply_mesh_simplification(
	const mesh_simplify::Options& options){

	mesh_simplify::Info information = mesh_simplify::simplify(
		vertices,
		indices,
		normals,
		options
	);

	edgeSet.clear();
	edgeIndices.clear();
	// update edges from indices
	
	for(size_t i{0}; i < indices.size() / 3; i++){
		add_edge(
			indices[3 * i],
			indices[3 * i + 1],
			indices[3 * i + 2]
		);
	}

	_setup_edges();

	// rebuild the mesh 
	update_buffers();
};

void GeneratorLewiner::smooth_scalar_field_taubin(int iterations, float lambda, float mu) {
    
    std::vector<float> tempField = scalarField;

    // Standard Taubin parameters: lambda > 0 (shrink), mu < 0 (inflate).
    // The requirement for stability is 0 < lambda < -mu < 1
    // Example defaults: lambda = 0.5f, mu = -0.53f

    // Per-axis Laplacian weights: with anisotropic voxels
    // (stepX != stepY != stepZ) a uniform 6-neighbour average smooths more,
    // in physical units, along the coarser axes. Weighting each axis by
    // 1/step^2 makes the filter isotropic in world space; for cubic voxels
    // this reduces exactly to the previous uniform (sum - 6f)/6 average.
    float wx = 1.0f, wy = 1.0f, wz = 1.0f;
    if (stepX > 0.0f && stepY > 0.0f && stepZ > 0.0f) {
        wx = 1.0f / (stepX * stepX);
        wy = 1.0f / (stepY * stepY);
        wz = 1.0f / (stepZ * stepZ);
    }
    const float wSum = wx + wy + wz;
    const float invNorm = 1.0f / (2.0f * wSum);

    // Exterior voxels are set to exactly air_skip_level() in
    // compute_scalar_field, so the skip test below must be >= (a strict >
    // never fires and the whole exterior gets convolved for nothing). The
    // band is scale-relative (a few voxels above the iso-level).
    const float skipLevel = air_skip_level();

    for (int iter = 0; iter < iterations; ++iter) {

        // --- PASS 1: SHRINK (Lambda) ---
        #pragma omp parallel for collapse(3)
        for (int z = 1; z < blockDims[2] - 1; z++) {
            for (int y = 1; y < blockDims[1] - 1; y++) {
                for (int x = 1; x < blockDims[0] - 1; x++) {

                    size_t idx = find_vertex_index(x, y, z);

                    // Optimization: Skip voxels safely outside the interaction margin
                    if (scalarField[idx] >= skipLevel) {
                        tempField[idx] = scalarField[idx];
                        continue;
                    }

                    // 3D discrete Laplacian using 6 face-neighbors, weighted per axis
                    float laplacian =
                        wx * (scalarField[find_vertex_index(x - 1, y, z)] +
                              scalarField[find_vertex_index(x + 1, y, z)]) +
                        wy * (scalarField[find_vertex_index(x, y - 1, z)] +
                              scalarField[find_vertex_index(x, y + 1, z)]) +
                        wz * (scalarField[find_vertex_index(x, y, z - 1)] +
                              scalarField[find_vertex_index(x, y, z + 1)]) -
                        (2.0f * wSum * scalarField[idx]);

                    // Normalize and apply lambda
                    laplacian *= invNorm;
                    tempField[idx] = scalarField[idx] + lambda * laplacian;
                }
            }
        }

        #pragma omp parallel for collapse(3)
        for (int z = 1; z < blockDims[2] - 1; z++) {
            for (int y = 1; y < blockDims[1] - 1; y++) {
                for (int x = 1; x < blockDims[0] - 1; x++) {

                    size_t idx = find_vertex_index(x, y, z);

                    if (tempField[idx] >= skipLevel) {
                        scalarField[idx] = tempField[idx];
                        continue;
                    }

                    // Calculate Laplacian using the tempField from Pass 1
                    float laplacian =
                        wx * (tempField[find_vertex_index(x - 1, y, z)] +
                              tempField[find_vertex_index(x + 1, y, z)]) +
                        wy * (tempField[find_vertex_index(x, y - 1, z)] +
                              tempField[find_vertex_index(x, y + 1, z)]) +
                        wz * (tempField[find_vertex_index(x, y, z - 1)] +
                              tempField[find_vertex_index(x, y, z + 1)]) -
                        (2.0f * wSum * tempField[idx]);

                    // Normalize and apply mu
                    laplacian *= invNorm;
                    scalarField[idx] = tempField[idx] + mu * laplacian;
                }
            }
        }
    }
}

void GeneratorLewiner::build_topology() {

	adjacency.clear();
	adjacency.resize(meshVertices.size());

	//int currentIdx = 0;
	for (const auto& f : meshTriangles) {

		adjacency[f.v1].push_back(f.v2);
		adjacency[f.v1].push_back(f.v3);
		adjacency[f.v2].push_back(f.v1);
		adjacency[f.v2].push_back(f.v3);
		adjacency[f.v3].push_back(f.v1);
		adjacency[f.v3].push_back(f.v2);

		//currentIdx++;
	}

	// Remove duplicates in adjacency
	for (auto& neighbors : adjacency) {
		std::sort(neighbors.begin(), neighbors.end());
		neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
	}
};

//@brief export the metrics in csv
void GeneratorLewiner::export_metrics(std::string fileName) {

	std::ofstream fout;
	fout.open(fileName);

	// add header
	fout << "Porosity, Mesh Porosity, Volume, TotalSurface, SurfaceToVolume, SurfaceToTotalVolume, Connectivity Density, SMI, Local Thickness, Local Thickness Std, Local Separation, Local Separation Std, trabecular Nr, Anisotropy, Tortuosity\n";

	// pass values
	fout << porosity << "," << porosityMesh << "," << volume << "," << surfaceArea << "," << surfaceToVolume << "," << surfaceToTotalVolume << "," << connectivityDensity << "," << smi << "," << localThickness << "," << localThicknessStd << "," << localSeparation << "," << localSeparationStd << "," << trabecularNr << "," << anisotropyDegree << "," << tortuosity << "\n";

	fout.close();

	logger->log(LogPriority::SUCCESS, "Exported metrics to " + fileName + "!");

};

void GeneratorLewiner::read_metrics(const std::string fileName) {
	std::ifstream fin(fileName);

	if (!fin.is_open()) {
		if (logger) logger->log(LogPriority::ERROR, "File " + fileName + " not successfully read. Metrics Not Loaded!");
		return;
	}

	// Lambda to trim whitespace
	auto trim = [](std::string& s) {
		s.erase(0, s.find_first_not_of(" \t\r\n"));
		s.erase(s.find_last_not_of(" \t\r\n") + 1);
	};

	std::string line;

	// 1. Read Header
	if (!std::getline(fin, line)) {
		if (logger) logger->log(LogPriority::ERROR, "Metrics file is empty!");
		return;
	}
	std::vector<std::string> header;
	std::stringstream headerStream(line);
	std::string cell;
	while (std::getline(headerStream, cell, ',')) {
		trim(cell);
		header.push_back(cell);
	}

	// 2. Read Data
	if (!std::getline(fin, line)) {
		if (logger) logger->log(LogPriority::ERROR, "Metrics data row is missing!");
		return;
	}
	std::vector<std::string> data;
	std::stringstream dataStream(line);
	while (std::getline(dataStream, cell, ',')) {
		trim(cell);
		data.push_back(cell);
	}

	if (header.size() != data.size()) {
		if (logger) logger->log(LogPriority::ERROR, "CSV Column mismatch in metrics file!");
		return;
	}

	// Instantiate your struct
	ScaffoldMetrics loadedMetrics;

	// 3. Map Data dynamically based on the exact strings in your export_metrics
	for (size_t i = 0; i < header.size(); ++i) {
		std::string key = header[i];
		float val = 0.0f;

		try {
			val = std::stof(data[i]);
		}
		catch (...) {
			continue; // Skip corrupted float data
		}

		if (key == "Porosity") { loadedMetrics.porosity = val; }
		else if (key == "Mesh Porosity") { loadedMetrics.porosityMesh = val; }
		else if (key == "Volume") { loadedMetrics.volume = val; }
		else if (key == "TotalSurface") { loadedMetrics.totalSurface = val; }
		else if (key == "SurfaceToVolume") { loadedMetrics.surfToVol = val; }
		else if (key == "SurfaceToTotalVolume") { loadedMetrics.surfToTotalVol = val; }
		else if (key == "SMI") { loadedMetrics.smi = val; }
		else if (key == "Connectivity Density") { loadedMetrics.connectivityDensity = val; }
		else if (key == "Local Thickness") { loadedMetrics.thickness = val; }
		else if (key == "Local Thickness Std") { loadedMetrics.thicknessStd = val; }
		else if (key == "Local Separation") { loadedMetrics.separation = val; }
		else if (key == "Local Separation Std") { loadedMetrics.separationStd = val; }
		else if (key == "trabecular Nr") { loadedMetrics.trNumber = val; }
		else if (key == "Anisotropy") { loadedMetrics.anisotropyDeg = val; }
		else if (key == "Tortuosity") { loadedMetrics.tortuosity = val; }
	}

	// 4. Apply the struct values back to your class members
	porosity = loadedMetrics.porosity;
	// old metric files have no "Mesh Porosity" column (-1 sentinel stays)
	if (loadedMetrics.porosityMesh >= 0.0f) {
		porosityMesh = loadedMetrics.porosityMesh;
	}
	volume = loadedMetrics.volume;
	surfaceArea = loadedMetrics.totalSurface;
	surfaceToVolume = loadedMetrics.surfToVol;
	// older metric files predate these columns (-1 sentinel stays)
	if (loadedMetrics.surfToTotalVol >= 0.0f) surfaceToTotalVolume = loadedMetrics.surfToTotalVol;
	if (loadedMetrics.smi >= 0.0f) smi = loadedMetrics.smi;
	connectivityDensity = loadedMetrics.connectivityDensity;
	localThickness = loadedMetrics.thickness;
	localThicknessStd = loadedMetrics.thicknessStd;
	localSeparation = loadedMetrics.separation;
	localSeparationStd = loadedMetrics.separationStd;
	trabecularNr = loadedMetrics.trNumber;
	anisotropyDegree = loadedMetrics.anisotropyDeg;
	tortuosity = loadedMetrics.tortuosity;

	if (logger) {
		logger->log(LogPriority::SUCCESS, "Successfully loaded metrics from " + fileName);
	}
}

//@brief export the applied parameters for scaffold in csv
void GeneratorLewiner::export_parameters(std::string fileName) {
	std::shared_ptr<InterfaceSeedGenerator> lockedGen = generator.lock();
	if (!lockedGen) return;

	std::ofstream fout(fileName);
	if (!fout.is_open()) return;

	ScaffoldParameters cfg;

	// Populate Core Parameters
	cfg.thicknessOption = selectedThicknessOption;
	cfg.openess = threshold;
	cfg.spread = spread;
	cfg.stretchX = stretchX; cfg.stretchY = stretchY; cfg.stretchZ = stretchZ;
	cfg.anisotropyAngle = anisotropyAngle;
	cfg.dirX = anisotropyVec.x; cfg.dirY = anisotropyVec.y; cfg.dirZ = anisotropyVec.z;
	cfg.backgroundWeight = backgroundWeight;
	cfg.smoothJunctions = smoothJunctions ? 1 : 0;
	cfg.smoothK = smoothK;
	cfg.roundEdges = roundEdges ? 1 : 0;
	cfg.edgeK = edgeK;
	cfg.frameBoundary = frameBoundary ? 1 : 0;
	cfg.frameDepth = frameDepth;
	cfg.frameContainerEdges = frameContainerEdges ? 1 : 0;
	cfg.frameBeam = frameBeam;
	cfg.calibrateThickness = calibrateThickness ? 1 : 0;
	cfg.targetThickness = targetThickness;
	cfg.calibratePorosity = calibratePorosity ? 1 : 0;
	cfg.targetPorosity = targetPorosity;
	cfg.calibrateDA = calibrateDA ? 1 : 0;
	cfg.targetDa = targetDa;
	cfg.targetFormulaIdx = targetFormulaIdx;
	cfg.calibrationStep = calibrationStep;
	cfg.calibrationTol = calibrationTol;
	cfg.calibrationIter = calibrationIter;
	cfg.transitionDistance = transitionDistance;

	// Populate Thickness & Distance Parameters
	if (selectedThicknessOption == 0) {
		cfg.uniformThickness = isoLevel;
	}
	else {
		cfg.startThickness = startThickness;
		cfg.endThickness = endThickness;
		cfg.distFunction = selectedDist;
		cfg.radFunction = selectedFunc;

		// Capture Vector Components
		cfg.planeOriginX = distancePlaneCenter.x; cfg.planeOriginY = distancePlaneCenter.y; cfg.planeOriginZ = distancePlaneCenter.z;
		cfg.planeNormalX = distancePlaneNormal.x; cfg.planeNormalY = distancePlaneNormal.y; cfg.planeNormalZ = distancePlaneNormal.z;
		cfg.pointX = distancePoint.x; cfg.pointY = distancePoint.y; cfg.pointZ = distancePoint.z;
	}

	// Populate Generator Parameters
	if (lockedGen->get_type() == ObjectType::RandomGeneratorType) {
		cfg.generatorType = 0;
		auto rgn = std::dynamic_pointer_cast<Random>(lockedGen);
		if (rgn) cfg.seedNr = static_cast<int>(rgn->get_seeds().size());
	}
	else if (lockedGen->get_type() == ObjectType::PoissonGeneratorType) {
		cfg.generatorType = 1;
		auto pgn = std::dynamic_pointer_cast<Poisson3D>(lockedGen);
		if (pgn) {
			cfg.minRadius = pgn->get_min_radius();
			if (!pgn->is_uniform()) {
				cfg.maxRadius = pgn->get_max_radius();
			}
		}
	}

	// Write Fixed Header (Expanded). New columns are appended at the END so files
	// written by older builds still load (the reader maps by key, not position).
	fout << "ThicknessOption,UniformThickness,StartThickness,EndThickness,DistFunction,RadFunction,"
		<< "PlaneOriginX,PlaneOriginY,PlaneOriginZ,PlaneNormalX,PlaneNormalY,PlaneNormalZ,PointX,PointY,PointZ,"
		<< "GeneratorType,SeedNr,MinRadius,MaxRadius,Openess,StretchX,StretchY,StretchZ,"
		<< "AnisotropyAngle,DirX,DirY,DirZ,"
		<< "Spread,BackgroundWeight,SmoothJunctions,SmoothK,TransitionDistance,FrameBoundary,FrameDepth,FrameContainerEdges,FrameBeam,"
		<< "CalibrateThickness,TargetThickness,CalibratePorosity,TargetPorosity,CalibrateDA,TargetDa,TargetFormulaIdx,CalibrationStep,CalibrationTol,CalibrationIter,"
		<< "RoundEdges,EdgeK\n";

	// Write Fixed Data Row
	fout << cfg.thicknessOption << "," << cfg.uniformThickness << "," << cfg.startThickness << "," << cfg.endThickness << "," << cfg.distFunction << "," << cfg.radFunction << ","
		<< cfg.planeOriginX << "," << cfg.planeOriginY << "," << cfg.planeOriginZ << ","
		<< cfg.planeNormalX << "," << cfg.planeNormalY << "," << cfg.planeNormalZ << ","
		<< cfg.pointX << "," << cfg.pointY << "," << cfg.pointZ << ","
		<< cfg.generatorType << "," << cfg.seedNr << "," << cfg.minRadius << "," << cfg.maxRadius << ","
		<< cfg.openess << "," << cfg.stretchX << "," << cfg.stretchY << "," << cfg.stretchZ << ","
		<< cfg.anisotropyAngle << "," << cfg.dirX << "," << cfg.dirY << "," << cfg.dirZ << ","
		<< cfg.spread << "," << cfg.backgroundWeight << "," << cfg.smoothJunctions << "," << cfg.smoothK << "," << cfg.transitionDistance << "," << cfg.frameBoundary << "," << cfg.frameDepth << "," << cfg.frameContainerEdges << "," << cfg.frameBeam << ","
		<< cfg.calibrateThickness << "," << cfg.targetThickness << "," << cfg.calibratePorosity << "," << cfg.targetPorosity << "," << cfg.calibrateDA << "," << cfg.targetDa << "," << cfg.targetFormulaIdx << "," << cfg.calibrationStep << "," << cfg.calibrationTol << "," << cfg.calibrationIter << ","
		<< cfg.roundEdges << "," << cfg.edgeK << "\n";

	fout.close();

	logger->log(LogPriority::SUCCESS, "Exported parameters to " + fileName + "!");
};

void GeneratorLewiner::read_parameters(const std::string fileName) {
	std::ifstream fin(fileName);

	if (!fin.is_open()) {
		logger->log(LogPriority::ERROR, "File " + fileName + " not successfully read");
		return;
	}

	auto trim = [](std::string& s) {
		s.erase(0, s.find_first_not_of(" \t\r\n"));
		s.erase(s.find_last_not_of(" \t\r\n") + 1);
	};

	std::string line;

	// Read Header
	if (!std::getline(fin, line)) return;
	std::vector<std::string> header;
	std::stringstream headerStream(line);
	std::string cell;
	while (std::getline(headerStream, cell, ',')) { trim(cell); header.push_back(cell); }

	// Read Data
	if (!std::getline(fin, line)) return;
	std::vector<std::string> data;
	std::stringstream dataStream(line);
	while (std::getline(dataStream, cell, ',')) { trim(cell); data.push_back(cell); }

	if (header.size() != data.size()) {
		logger->log(LogPriority::ERROR, "CSV Column mismatch!");
		return;
	}

	// Unified mapping loop
	for (size_t i = 0; i < header.size(); ++i) {
		std::string key = header[i];
		float val = 0.0f;
		try { val = std::stof(data[i]); }
		catch (...) { continue; }

		if (val == -1.0f) continue; // Ignore unused placeholders

		if (key == "ThicknessOption") { selectedThicknessOption = static_cast<int>(val); }
		else if (key == "UniformThickness") { isoLevel = val; }
		else if (key == "StartThickness") { startThickness = val; }
		else if (key == "EndThickness") { endThickness = val; }
		else if (key == "DistFunction") { selectedDist = static_cast<int>(val); }
		else if (key == "RadFunction") { selectedFunc = static_cast<int>(val); }

		// Map Vector Components Back to Vec3
		else if (key == "PlaneOriginX") { distancePlaneCenter.x = val; }
		else if (key == "PlaneOriginY") { distancePlaneCenter.y = val; }
		else if (key == "PlaneOriginZ") { distancePlaneCenter.z = val; }
		else if (key == "PlaneNormalX") { distancePlaneNormal.x = val; }
		else if (key == "PlaneNormalY") { distancePlaneNormal.y = val; }
		else if (key == "PlaneNormalZ") { distancePlaneNormal.z = val; }
		else if (key == "PointX") { distancePoint.x = val; }
		else if (key == "PointY") { distancePoint.y = val; }
		else if (key == "PointZ") { distancePoint.z = val; }

		else if (key == "Openess") { threshold = val; }
		else if (key == "StretchX") { stretchX = val; }
		else if (key == "StretchY") { stretchY = val; }
		else if (key == "StretchZ") { stretchZ = val; }
		else if (key == "AnisotropyAngle") { anisotropyAngle = val; }
		else if (key == "DirX") { anisotropyVec.x = val; }
		else if (key == "DirY") { anisotropyVec.y = val; }
		else if (key == "DirZ") { anisotropyVec.z = val; }

		else if (key == "Spread") { spread = val; }
		else if (key == "BackgroundWeight") { backgroundWeight = val; }
		else if (key == "SmoothJunctions") { smoothJunctions = (val != 0.0f); }
		else if (key == "SmoothK") { smoothK = val; }
		else if (key == "RoundEdges") { roundEdges = (val != 0.0f); }
		else if (key == "EdgeK") { edgeK = val; }
		else if (key == "FrameBoundary") { frameBoundary = (val != 0.0f); }
		else if (key == "FrameDepth") { frameDepth = val; }
		else if (key == "FrameContainerEdges") { frameContainerEdges = (val != 0.0f); }
		else if (key == "FrameBeam") { frameBeam = val; }
		else if (key == "CalibrateThickness") { calibrateThickness = (val != 0.0f); }
		else if (key == "TargetThickness") { targetThickness = val; }
		else if (key == "CalibratePorosity") { calibratePorosity = (val != 0.0f); }
		else if (key == "TargetPorosity") { targetPorosity = val; }
		else if (key == "CalibrateDA") { calibrateDA = (val != 0.0f); }
		else if (key == "TargetDa") { targetDa = val; }
		else if (key == "TargetFormulaIdx") { targetFormulaIdx = static_cast<int>(val); }
		else if (key == "CalibrationStep") { calibrationStep = val; }
		else if (key == "CalibrationTol") { calibrationTol = val; }
		else if (key == "CalibrationIter") { calibrationIter = static_cast<int>(val); }
		else if (key == "TransitionDistance") { transitionDistance = val; }
	}

	logger->log(LogPriority::SUCCESS, "Successfully loaded parameters from " + fileName);
}

std::unique_ptr<GeneratorLewiner> GeneratorLewiner::extract_from_ROI(ROI* roi) {

	std::array<float, 6> reqBounds = roi->get_bounds();

	// Calculate the exact integer grid indices in the parent mesh that enclose the ROI
	int min_i = static_cast<int>(std::floor((reqBounds[0] - this->bounds[0]) / this->stepX));
	int max_i = static_cast<int>(std::ceil((reqBounds[1] - this->bounds[0]) / this->stepX));

	int min_j = static_cast<int>(std::floor((reqBounds[2] - this->bounds[2]) / this->stepY));
	int max_j = static_cast<int>(std::ceil((reqBounds[3] - this->bounds[2]) / this->stepY));

	int min_k = static_cast<int>(std::floor((reqBounds[4] - this->bounds[4]) / this->stepZ));
	int max_k = static_cast<int>(std::ceil((reqBounds[5] - this->bounds[4]) / this->stepZ));

	// Safety Clamp: Ensure we don't accidentally expand outside the parent's actual dimensions
	min_i = std::max(0, min_i);  max_i = std::min(this->blockDims[0] - 1, max_i);
	min_j = std::max(0, min_j);  max_j = std::min(this->blockDims[1] - 1, max_j);
	min_k = std::max(0, min_k);  max_k = std::min(this->blockDims[2] - 1, max_k);

	// Reconstruct the perfectly snapped physical bounds
	std::array<float, 6> snappedBounds = {
		this->bounds[0] + min_i * this->stepX,
		this->bounds[0] + max_i * this->stepX,
		this->bounds[2] + min_j * this->stepY,
		this->bounds[2] + max_j * this->stepY,
		this->bounds[4] + min_k * this->stepZ,
		this->bounds[4] + max_k * this->stepZ
	};

	// Calculate exact integer resolution (No floats, no std::ceil needed here)
	std::array<int, 3> localResolution = {
		max_i - min_i + 1,
		max_j - min_j + 1,
		max_k - min_k + 1
	};

	auto roiScaffold = std::make_unique<GeneratorLewiner>(
		this->seeds,
		snappedBounds,
		localResolution,
		this->logger,
		this->threshold,
		this->isoLevel,
		this->renderMode
	);

	roiScaffold->set_stretch(this->stretchX, this->stretchY, this->stretchZ);
	roiScaffold->anisotropyAngle = this->anisotropyAngle;
	roiScaffold->anisotropyVec = this->anisotropyVec;
	roiScaffold->anisotropySources = this->anisotropySources;
	roiScaffold->set_thickness_functions(this->thicknessSDF, this->thicknessFunction, this->startThickness, this->endThickness, this->transitionDistance);
	roiScaffold->isROI = true;
	roiScaffold->name = this->name + " (ROI)";

	// Generate the mesh using the parent's container
	std::shared_ptr<IContainer> parentCon = this->container.lock();
	if (parentCon) {
		if (roiScaffold->compute_scalar_field(*parentCon) &&
			roiScaffold->marching_cubes()) {
			roiScaffold->estimate_metrics(*parentCon);
		}
	}

	roiScaffold->update_render();
	return roiScaffold;
}

// =========================================================================
// Scaffold Factory Class Implementation
// =========================================================================
void ScaffoldFactory::launch() {

	selectedCon.reset();
	selectedGen.reset();
	lockedCon.reset();
	lockedGen.reset();
	name="";
	genContainerName = "";
	genGeneratorName = "";
	thickness = { 0.3f };
	openess = { 0.5f };
	stretchX = { 1.0f };
	stretchY = { 1.0f };
	stretchZ = { 1.0f };
	anisotropyVec = { 1.0f, 0.0f, 0.0f };
	anisotropyAngle = { 0.0f };
	anisotropySources.clear();
	voxelSize = 0.05f;

	// for thickness function
	selectedThicknessOption = 0;
	distancePlaneNormal = { 0.0f, 0.0f, 1.0f };
	distancePlaneCenter = { 0.0f, 0.0f, 0.0f };
	distancePoint = { 0.0f, 0.0f, 0.0f };
	transitionDistance = 10.0f;

	iter = 15;
	lambda = 0.5f;
	mu = -0.53f;

	warningFlashTimer1 = 0.0f;
	warningFlashTimer2 = 0.0f;
	thicknessRadFunc.reset();
	thicknessSDF.reset();
};

void ScaffoldFactory::gui_draw(
	GenerationTask* task,
	Logger* logger,
	const char* popupName, bool& showPopup,
	SelectedObject* selectedPanelObj, void*& selectedSceneObj,
	std::vector<std::unique_ptr<GeneratorLewiner>>& scaffoldList,
	std::vector<std::shared_ptr<IContainer>>& containers,
	std::vector<std::shared_ptr<InterfaceSeedGenerator>>& generators,
	std::vector<std::shared_ptr<AnisotropySource>>& anisoSources
) {
	
	// set to already created sources
	anisotropySources = anisoSources;

	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}

	// always centered
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
	ImGui::SetNextWindowSize(ImVec2(600, 500), ImGuiCond_Appearing);

	if (ImGui::BeginPopupModal("Scaffold Creator", NULL))
	{
		
        const bool anyBusy = task->get_running();
		const bool mineBusy = task->is_running_for(this);

        // ---------- Generate ----------
        ImGui::BeginDisabled(anyBusy);
        
		ImGuiTabBarFlags tab_bar_flags = ImGuiTabBarFlags_None;
		if (ImGui::BeginTabBar("MyTabBar", tab_bar_flags)){

			main_options(containers, generators);

			// calibration_options();

			thickness_options();

			anisotropy_options(anisoSources);

			smooth_options();

			ImGui::EndTabBar();
		}

		ImGui::Separator();
		
		if (ImGui::Button("Generate")) {
            if (!lockedCon) warningFlashTimer1 = 1.5f;
            if (!lockedGen) warningFlashTimer2 = 1.5f;

            if (lockedCon && lockedGen) {
                auto conShared = lockedCon;                 
                auto seeds     = lockedGen->get_seeds();    
                auto bds       = conShared->compute_bounds();
				auto targetThickness = -1.0f;
				auto targetP = targetPorosity * 0.01f;

				if(selectedThicknessOption == 0){
					targetThickness = thickness;
				}
				else{
					// varied thickness: calibrate the MEAN of the graded range 
					targetThickness = (startThickness + endThickness) * 0.5f;
				}

                std::array<float, 6> bounds = {
                    bds.xMin, bds.xMax, bds.yMin, bds.yMax, bds.zMin, bds.zMax
                };
                resolution = {
                    static_cast<int>(std::ceil((bds.xMax - bds.xMin) / voxelSize)) + 1,
                    static_cast<int>(std::ceil((bds.yMax - bds.yMin) / voxelSize)) + 1,
                    static_cast<int>(std::ceil((bds.zMax - bds.zMin) / voxelSize)) + 1
                };

                // --- Memory guard --------------------------------------------------
                // The scalar-field allocations scale with the voxel count; a fine grid
                // over a large container can exceed RAM and hard-crash with bad_alloc.
                // Estimate the peak footprint and refuse up front with an actionable
                // message instead of crashing on the worker thread.
                const uint64_t projectedVoxels =
                    (uint64_t)resolution[0] * (uint64_t)resolution[1] * (uint64_t)resolution[2];
                const bool     willCalibrate = calibrateThickness || calibratePorosity;
                // peak bytes/voxel: cache(40)+field(4)+containerDist(4)+narrowBand(4)+MC verts(24)
                // when calibrating; without the cache the plain path drops the 40.
                const uint64_t bytesPerVoxel = willCalibrate ? 76ull : 36ull;
                const uint64_t requiredBytes = projectedVoxels * bytesPerVoxel;
                const uint64_t freeBytes     = available_physical_memory_bytes();
                const double   budgetFrac    = 0.80; // leave headroom for the OS/other buffers

                const bool memoryOk =
                    (freeBytes == 0) || // query unavailable: don't block (bad_alloc backstop applies)
                    (requiredBytes <= (uint64_t)(budgetFrac * (double)freeBytes));

                if (!memoryOk) {
                    const double dimVol = (double)(bds.xMax - bds.xMin)
                                        * (double)(bds.yMax - bds.yMin)
                                        * (double)(bds.zMax - bds.zMin);
                    const double minVoxel = std::cbrt(
                        dimVol * (double)bytesPerVoxel / (budgetFrac * (double)freeBytes));
                    const double minVoxelR = std::ceil(minVoxel * 1000.0) / 1000.0;
                    std::ostringstream oss;
                    oss << "Generation aborted: grid is " << projectedVoxels << " voxels (~"
                        << (requiredBytes >> 30) << " GB needed, " << (freeBytes >> 30)
                        << " GB free). Increase the generation voxel size to >= "
                        << minVoxelR << " mm"
                        << (willCalibrate
                            ? "; metrics can still be measured at the measurement voxel inside an ROI."
                            : ".");
                    if (logger) logger->log(LogPriority::ERROR, oss.str());
                }
                else {
                auto scaffold = std::make_unique<GeneratorLewiner>(
                    seeds, bounds, resolution, logger, openess, thickness);

				scaffold->spread = spread;
				scaffold->smoothJunctions = smoothJunctions;
				scaffold->smoothK = smoothK;
				scaffold->roundEdges = roundEdges;
				scaffold->edgeK = edgeK;
				scaffold->frameBoundary = frameBoundary;
				scaffold->frameDepth = frameDepth;
				scaffold->frameContainerEdges = frameContainerEdges;
				scaffold->frameBeam = frameBeam;

                if (selectedThicknessOption == 1) {
                    switch (selectedFunc) {
                        case 0: thicknessRadFunc = std::make_shared<LinearFunction>(transitionDistance);    break;
                        case 1: thicknessRadFunc = std::make_shared<QuadraticFunction>(transitionDistance); break;
						case 2: thicknessRadFunc = std::make_shared<SmoothStep>(); break;
                        case 3: thicknessRadFunc = std::make_shared<ConstantRadiusFunction>();              break;
                    }
                    switch (selectedDist) {
                        case 0:
                            thicknessSDF = std::make_shared<PlaneSDF>(distancePlaneCenter, distancePlaneNormal);
                            scaffold->set_distance_plane_options(distancePlaneCenter, distancePlaneNormal);
                            break;
                        case 1:
                            thicknessSDF = std::make_shared<PointSDF>(distancePoint);
                            scaffold->set_distance_point_options(distancePoint);
                            break;
                        case 2:
                            thicknessSDF = conShared->get_distance_estimator();
                            break;
                    }
                } else {
                    thicknessSDF.reset();
                    thicknessRadFunc.reset();
                }

                scaffold->set_options_from_factory(selectedDist, selectedFunc, selectedThicknessOption, voxelSize, measurementVoxelSize);
                scaffold->set_thickness_functions(
					thicknessSDF, thicknessRadFunc,
                    startThickness, endThickness, transitionDistance);
                scaffold->set_stretch(stretchX, stretchY, stretchZ);
				
				// update calibration settings (flags + targets consumed by
				// build_calibration_stages, which runs on the worker thread)
				scaffold->calibrateThickness = calibrateThickness;
				scaffold->calibratePorosity  = calibratePorosity;
				// scaffold->calibrateSmi       = calibrateSmi;   // SMI calibration disabled
				scaffold->calibrateDA        = calibrateDA;
				scaffold->targetThickness    = targetThickness;
				scaffold->targetPorosity     = targetPorosity;
				// scaffold->targetSmi          = targetSmi;      // SMI calibration disabled
				scaffold->targetDa           = targetDa;
				scaffold->targetFormulaIdx   = targetFormulaIdx;
				scaffold->calibrationIter    = calibrationIter;
				scaffold->calibrationTol     = calibrationTol;

                scaffold->anisotropyAngle = anisotropyAngle;
                scaffold->anisotropyVec   = anisotropyVec;
				scaffold->anisotropySources = anisoSources;
				scaffold->backgroundWeight = backgroundWeight;
                scaffold->container = lockedCon;
                scaffold->generator = lockedGen;
				scaffold->iter = iter;
				scaffold->mu = mu;
				scaffold->lambda = lambda;

                genContainerName = lockedCon->name;   // capture identity now, for the log later
                genGeneratorName = lockedGen->name;

                GeneratorLewiner* raw = scaffold.get();
                pendingScaffold = std::move(scaffold);

                // wire cooperative cancellation to the running task's cancel flag
                raw->cancelRequested = [task]{ return task->is_cancel_requested(); };

                start_time = std::chrono::steady_clock::now();

				// NOTE: the expensive cache pass runs on the WORKER thread, not here
				// (it would freeze the UI). Calibrators lazy-cache; the plain path
				// calls compute_scalar_field (cache + assemble).

					// Calibration: build the enabled knob stack (DA -> SMI ->
					// Thickness -> Porosity) and solve it with one recursive driver;
					// any subset reproduces the old one/two/three-knob behaviour.
				const bool anyCalibration =
					calibrateThickness || calibratePorosity /* || calibrateSmi */ || calibrateDA;

				if (anyCalibration) {
					auto cStep = measurementVoxelSize; // metrics at the measurement voxel
					task->start(
						[raw, conShared, cStep, t = task]() {
						t->set_progress(0.00f);

						raw->build_calibration_stages(cStep);
						bool ok = !t->is_cancel_requested()
							&& raw->solve_calibration(raw->stages, 0, *conShared);

						if (ok && !t->is_cancel_requested() && raw->marching_cubes()) {
							t->set_progress(0.90f);
							// full-resolution DA (the loop used a reduced ray budget)
							if (raw->calibrateDA)
								raw->estimate_anisotropy(cStep, 2000, 10000, raw->targetFormulaIdx);
							raw->estimate_metrics(*conShared);
							// stamp calibrated metric versions to this final mesh so the
							// inner-loop targets (Tb.Th, SMI) no longer read stale.
							raw->stamp_calibrated_versions();
						}
						t->set_progress(1.00f);
					}, this);
				}
				else {
					task->start(
						[raw, conShared, t = task]() {
						t->set_progress(0.00f);
						if (!t->is_cancel_requested() &&
							raw->compute_scalar_field(*conShared)) {
							t->set_progress(0.50f);
							if (!t->is_cancel_requested() && raw->marching_cubes()) {
								t->set_progress(0.90f);
								if (!t->is_cancel_requested()) raw->estimate_metrics(*conShared);
							}
						}
						t->set_progress(1.00f);
					}, this);
				}
                } // end memory guard (memoryOk)
            }
        }
        ImGui::EndDisabled();

		ImGui::SameLine();

		// While a job is running, Cancel ends the parallel task (the popup stays
		// open until poll() joins the thread, then discards the partial scaffold).
		// With nothing running, Cancel just closes the creator dialog.
		const bool cancelling = mineBusy && task->is_cancel_requested();
		ImGui::BeginDisabled(cancelling);
		if (ImGui::Button(cancelling ? "Cancelling..." : "Cancel")) {
			if (mineBusy) task->request_cancel();
			else          showPopup = false;
		};
		ImGui::EndDisabled();

        // ---------- Progress ----------
        if (mineBusy) {
            ImGui::ProgressBar(task->get_progress(),
			 ImVec2(-FLT_MIN, 0.0f));
        }
		else{
			ImGui::Dummy(ImVec2(0.0f, ImGui::GetFrameHeight()));
		}

        // ---------- Completion (main thread: GL upload + register) ----------
        if (task->poll(this)) {
          if (task->is_cancel_requested()) {
            // user cancelled: discard the partial scaffold, do not register it
            pendingScaffold.reset();
            logger->log(LogPriority::WARNING, "Scaffold creation cancelled.");
            showPopup = false;
            ImGui::CloseCurrentPopup();
          }
          else {
            pendingScaffold->update_render();    // GPU upload happens here, on the GL thread

            pendingScaffold->name = name.empty()
                ? "Scaffold" + std::to_string(scaffoldList.size() + 1)
                : name;

            scaffoldList.push_back(std::move(pendingScaffold));
            selectedSceneObj       = scaffoldList.back().get();
            selectedPanelObj->ptr  = scaffoldList.back().get();
            selectedPanelObj->type = ObjectType::ScaffoldType;

            auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - start_time);
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(3) << dur.count() / 1000.0 << " seconds!";
            logger->log(LogPriority::SUCCESS,
                "Created scaffold successfully using container " + genContainerName +
                " and generator " + genGeneratorName + " in " + oss.str());

            showPopup = false;
            ImGui::CloseCurrentPopup();
          }
        }
		ImGui::EndPopup();
	}
};

void ScaffoldFactory::thickness_options(){
	if (ImGui::BeginTabItem("Thickness")){

		ImGui::Checkbox("Calibrate Thickness", &calibrateThickness);

		ImGui::Separator();

		ImGui::RadioButton("Apply Uniform Thickness", &selectedThicknessOption, 0);
		ImGui::RadioButton("Apply Varied Thickness", &selectedThicknessOption, 1);

		if (selectedThicknessOption == 0) {
			ImGui::InputFloat("Thickness", &thickness, 0.001f, 1.0f);
			// uniform: the Thickness value IS the calibration target (Tb.Th)
			if (calibrateThickness)
				ImGui::Text("Calibration target (Tb.Th): %.3f mm", thickness);
		}
		else {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Start Thickness", &startThickness);
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("End Thickness", &endThickness);
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Transition Distance", &transitionDistance);
			// varied: the calibration target is the MEAN of the (scaled) range
			if (calibrateThickness)
				ImGui::Text("Calibration target (mean Tb.Th): %.3f mm",
					(startThickness + endThickness) * 0.5f);

			ImGui::SeparatorText("Select Distance Function");
			ImGui::RadioButton("Distance From Plane", &selectedDist, 0);
			if (selectedDist == 0) {
				ImGui::SetNextItemWidth(200);
				ImGui::InputFloat3("Normal", distancePlaneNormal);
				ImGui::SetNextItemWidth(200);
				ImGui::InputFloat3("Center", distancePlaneCenter);
			};
			ImGui::RadioButton("Distance From Point", &selectedDist, 1);
			if (selectedDist == 1) {
				ImGui::SetNextItemWidth(200);
				ImGui::InputFloat3("Point", distancePoint);
			}
			ImGui::RadioButton("Distance From Container", &selectedDist, 2);
			
			ImGui::SeparatorText("Select Radius Function");
			ImGui::RadioButton("Linear", &selectedFunc, 0);
			ImGui::RadioButton("Quadratic", &selectedFunc, 1);
			ImGui::RadioButton("Constant", &selectedFunc, 2);
			ImGui::RadioButton("Random", &selectedFunc, 3);
		}

		ImGui::SliderFloat("Openess", &openess, 0.0f, 1.0f);
		ImGui::SliderFloat("Spread", &spread, 0.0f, 1.0f);
		ImGui::EndTabItem();
	}
};

void ScaffoldFactory::anisotropy_options(
	std::vector<std::shared_ptr<AnisotropySource>>& globalSources
){

	if(ImGui::BeginTabItem("Anisotropy")){

		ImGui::Checkbox("Calibrate Anisotropy", &calibrateDA);
		ImGui::InputFloat("Calibration Target", &targetDa);		

		ImGui::SeparatorText("Global Background");
		ImGui::InputFloat("Stretch X", &stretchX, 0.01f, 5.0f, "%.3f");
		ImGui::InputFloat("Stretch Y", &stretchY, 0.01f, 5.0f, "%.3f");
		ImGui::InputFloat("Stretch Z", &stretchZ, 0.01f, 5.0f, "%.3f");
		ImGui::InputFloat("Angle", &anisotropyAngle, 0.01f, 5.0f, "%.3f");
		ImGui::InputFloat("Background Weight", &backgroundWeight, 0.01f, 0.5f, "%.3f");
		ImGui::InputFloat3("Stretch Direction", anisotropyVec, "%.3f");
		
		ImGui::SameLine();
		if (ImGui::Button("Normalize")) {
			float len = std::sqrt(anisotropyVec.x*anisotropyVec.x + 
								anisotropyVec.y*anisotropyVec.y + 
								anisotropyVec.z*anisotropyVec.z);
			if (len > 1e-6f) {
				anisotropyVec.x /= len;
				anisotropyVec.y /= len;
				anisotropyVec.z /= len;
			}
		}
		
		for (int i{ 0 }; i < (int)anisotropySources.size(); i++) {
			std::string label = anisotropySources[i]->name.empty()
				? "Anisotropy Source " + std::to_string(i + 1)
				: anisotropySources[i]->name;
			if (ImGui::TreeNode(label.c_str())) {
				anisotropySources[i]->render_properties();
				ImGui::TreePop();
			}
		}

		if (ImGui::Button("Add Anisotropy Source")){
			std::shared_ptr<AnisotropySource> source = std::make_shared<AnisotropySource>();

			if(source->name.empty()){
				source->name = "Anisotropy Source" + std::to_string(
					globalSources.size() + 1
				);
			}

			source->update_metric();
			source->update_model();
			anisotropySources.push_back(source);
			globalSources.push_back(source);
		}
		ImGui::EndTabItem();
	};
};

void ScaffoldFactory::main_options(
	std::vector<std::shared_ptr<IContainer>>& containers,
	std::vector<std::shared_ptr<InterfaceSeedGenerator>>& generators
){

	if(ImGui::BeginTabItem("Main")){

		ImGui::InputText("Name", &name);
		ImGui::InputFloat("Voxel Size", &voxelSize);
		ImGui::InputFloat("Measurement Voxel Size", &measurementVoxelSize);
		
		lockedCon = selectedCon.lock();
		lockedGen = selectedGen.lock();

		// here we should get the seeds from the corresponding creator
		ImGui::SeparatorText("Select Container and generator");

		// timer for flashing
		warningFlashTimer1 = 0.0f;

		if (warningFlashTimer1 > 0.0f) {
			warningFlashTimer1 -= ImGui::GetIO().DeltaTime;
		}

		bool isFlashing1 = (warningFlashTimer1 > 0.0f);
		if (isFlashing1) {
			float pulseAlpha = (float)(std::sin(ImGui::GetTime() * 15.0f) * 0.5f + 0.5f);
			ImVec4 flashColor = ImVec4(1.0f, 0.0f, 0.0f, pulseAlpha);
			ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 5.0f);
			ImGui::PushStyleColor(ImGuiCol_Border, flashColor);
		}

		ImGui::BeginChild("Containers", ImVec2(ImGui::GetContentRegionAvail().x * 0.5f, 80), ImGuiChildFlags_Borders);

		for (const auto& md : containers) {

			bool isSelected = (lockedCon && lockedCon == md);
			if (ImGui::Selectable(md->name.c_str(), isSelected)) {
				selectedCon = md;
				lockedCon = md;
			};
		}

		if (isFlashing1) {
			ImGui::PopStyleColor(); // Pop the red border color
			ImGui::PopStyleVar();   // Pop the thick border size
		}

		ImGui::EndChild();

		ImGui::SameLine();

		warningFlashTimer2 = 0.0f;

		if (warningFlashTimer2 > 0.0f) {
			warningFlashTimer2 -= ImGui::GetIO().DeltaTime;
		}
		bool isFlashing2 = (warningFlashTimer2 > 0.0f);
		if (isFlashing2) {
			float pulseAlpha = (float)(std::sin(ImGui::GetTime() * 15.0f) * 0.5f + 0.5f);
			ImVec4 flashColor = ImVec4(1.0f, 0.0f, 0.0f, pulseAlpha);
			ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 5.0f);
			ImGui::PushStyleColor(ImGuiCol_Border, flashColor);
		}

		ImGui::BeginChild("Generators", ImVec2(0.0, 80), ImGuiChildFlags_Borders);

		for (const auto& md : generators) {

			bool isSelected = (lockedGen && lockedGen == md);

			if (ImGui::Selectable(md->name.c_str(), isSelected)) {
				selectedGen = md;
				lockedGen = md;
			};
		}

		if (isFlashing2) {
			ImGui::PopStyleColor();
			ImGui::PopStyleVar();  
		}

		ImGui::EndChild();	

		ImGui::SeparatorText("Porosity Calibration");

		ImGui::Checkbox("Calibrate Porosity", &calibratePorosity);
		ImGui::InputFloat("Target Porosity", &targetPorosity);

		// SMI calibration disabled -- SMI is a reported outcome, not a target:
		// ImGui::Checkbox("Calibrate SMI", &calibrateSmi);
		// ImGui::InputFloat("Target SMI", &targetSmi);

		ImGui::EndTabItem();
	};
};

void ScaffoldFactory::smooth_options(){

	if(ImGui::BeginTabItem("Smoothness")){

		ImGui::SeparatorText("Mesh smoothing (Taubin)");
		ImGui::InputInt("Iterations", &iter, 1);

		ImGui::InputFloat("Lambda", &lambda, 0.01f, 10.0f);

		ImGui::InputFloat("Mu", &mu, 0.01f, 10.0f);

		ImGui::SeparatorText("Junction smoothing");
		ImGui::Checkbox("Smooth rod-plate junctions", &smoothJunctions);
		ImGui::SetItemTooltip("Fillet/fatten the trabecular nodes (more organic, "
			"bone-like). Tb.Th inflation is absorbed by thickness calibration.");
		if (smoothJunctions) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Junction k", &smoothK, 0.001f, 0.1f, "%.4f");
		}

		ImGui::SeparatorText("Edge rounding");
		ImGui::Checkbox("Round cell edges", &roundEdges);
		ImGui::SetItemTooltip("Soften the raw Voronoi distance field so the cell "
			"EDGES and VERTICES become round rather than polygonal (rounder "
			"pores). Unlike junction smoothing it acts before the openness blend, "
			"so it also rounds bare (open-lattice) struts. Tb.Th inflation is "
			"absorbed by thickness calibration. Requires re-generation (recache).");
		if (roundEdges) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Edge k", &edgeK, 0.001f, 0.1f, "%.4f");
		}

		ImGui::SeparatorText("Boundary frame");
		ImGui::Checkbox("Frame boundary cells", &frameBoundary);
		ImGui::SetItemTooltip("Close the outermost Voronoi cells into full walls, "
			"forming a connected honeycomb rim around the scaffold. Ties the cut "
			"boundary struts together for gripping/loading test specimens.");
		if (frameBoundary) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Frame depth (x thickness)", &frameDepth, 0.1f, 0.5f, "%.2f");
			ImGui::SetItemTooltip("Rim depth as a multiple of the wall thickness. "
				"Keep it small (~0.5) so the boundary stays a thin honeycomb rim; "
				"large values turn the boundary shell plate-like.");
		}
		ImGui::Checkbox("Frame container edges", &frameContainerEdges);
		ImGui::SetItemTooltip("Add solid beams along the container's own edges "
			"(a box's 12 edges, a cylinder's two rims) for a rigid outer cage. "
			"Box and cylinder containers only.");
		if (frameContainerEdges) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Edge beam (mm)", &frameBeam, 0.05f, 0.5f, "%.3f");
		}

		ImGui::EndTabItem();
	}
};