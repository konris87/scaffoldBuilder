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
	const int foam,
	const bool renderMode
	) : seeds(seeds), bounds(bounds), blockDims(dims), logger(uiLogger), threshold(threshold), isoLevel(isoLevel), foam(foam), renderMode(renderMode) {

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
void GeneratorLewiner::compute_scalar_field(const IContainer& con) {

	// we have to control the number of anisotropy sources. 
	const bool variedAnisotropy = anisotropySources.size() >= 1;

	// we use the general anisotropy tensor as a background
	AnisotropySource background;
	background.direction = anisotropyVec;
	background.stretch = Vec3(stretchX, stretchY, stretchZ);
	create_metric(background);

	// ensure every source's M is up-to-date (GUI edits don't call create_metric)
	for (auto& src : anisotropySources) {
		src->update_metric();
	}

	// use also the domain volume
	domainVolume = con.get_volume();

	update_steps();

	// Eigen::Matrix3f rot = rotation_axis_angle(anisotropyVec, anisotropyAngle);
	Eigen::Matrix3f rot = rotation_from_direction_roll(
		anisotropyVec, anisotropyAngle);

	Vec3 center = con.compute_bounds().center;

	std::vector<Vec3> warpedSeeds = this->seeds;

	if (!variedAnisotropy){
		for (auto& seed : warpedSeeds) {
			Vec3 local = seed - center;
			Vec3 rotated = Vec3(rot * Eigen::Vector3f{ local.x, local.y, local.z });

			seed.x = rotated.x / stretchX;
			seed.y = rotated.y / stretchY;
			seed.z = rotated.z / stretchZ;
		}
	}

	// first populate the kdtree with the seeds (warped in case of uniform anisotropy, unwarped in case of varied anisotropy)
	std::unique_ptr<Kdtree>kdtree = std::make_unique<Kdtree>(warpedSeeds);

	// decide the number k of nn for the kdtree
	std::vector<float> stretches = {stretchX, stretchY, stretchZ};
	for (const auto& src: anisotropySources){
		stretches.push_back(src->stretch.x);
		stretches.push_back(src->stretch.y);
		stretches.push_back(src->stretch.z);
	}

	auto minmax = std::minmax_element(stretches.begin(), stretches.end());
	float sMin = *minmax.first;
	float sMax = *minmax.second;

	size_t Kn = variedAnisotropy ? choose_candidate_number(
		sMin, sMax, seeds.size()) : 3;

	// update scalar field
	scalarField.clear();

	// it will be more efficient to reserve the size of the scalar field vector before filling it, 
	// since we know the size of the grid in advance
	size_t totalVoxels = static_cast<size_t>(blockDims[0]) *
		static_cast<size_t>(blockDims[1]) *
		static_cast<size_t>(blockDims[2]);

	std::cout << "Voxel Nr: " << totalVoxels << std::endl;

	scalarField.resize(totalVoxels, 9999.9f);

	// cache the container distance computed in the first pass so the
	// boundary-clamp pass below doesn't have to query the SDF a second time
	std::vector<float> containerDistField(totalVoxels);

	// "more than this far outside, don't bother" margin for the early-out /
	// boundary-clamp use of the container distance
	const float surfaceMargin = 2.0f;

	// When "distance from container" varied thickness is active, the
	// thickness lookup below reuses this same container distance and needs
	// it to be meaningful out to transitionDistance, not just surfaceMargin.
	// If the narrow-band classification's placeholder for "Inside" voxels
	// only guaranteed >surfaceMargin (e.g. 2mm) while transitionDistance is
	// larger (e.g. 5mm), every voxel between those two depths would be fed
	// the same fixed placeholder instead of its real distance, collapsing
	// the intended smooth thickness taper into a wrong constant with a
	// discontinuity at the surfaceMargin shell. Classifying against the
	// larger of the two margins keeps the placeholder valid for both uses
	// (narrower/cheaper when thickness doesn't depend on it, conservatively
	// wider - more voxels fall back to exact evaluation - when it does).
	const bool thicknessReusesContainerSDF =
		thicknessFunction && thicknessSDF && (thicknessSDF.get() == con.sdf.get());
	const float classificationMargin = thicknessReusesContainerSDF
		? std::max(surfaceMargin, transitionDistance)
		: surfaceMargin;

	const auto fillStart = std::chrono::steady_clock::now();

	std::vector<NarrowBandClass> narrowBand = classify_container_narrow_band(con, classificationMargin);

	const auto classifyEnd = std::chrono::steady_clock::now();
	std::cout << "  [compute_scalar_field] narrow-band classification: "
		<< std::chrono::duration<double>(classifyEnd - fillStart).count() << " s" << std::endl;

	#pragma omp parallel for collapse(3)
	// create the grid points based on the bounds and block dimensions
	for (int i{ 0 }; i < blockDims[0]; i++) {
		for (int j{ 0 }; j < blockDims[1]; j++) {
			for (int k{ 0 }; k < blockDims[2]; k++) {

				size_t idx = find_vertex_index(i, j, k);
				float x = bounds[0] + i * stepX;
				float y = bounds[2] + j * stepY;
				float z = bounds[4] + k * stepZ;
				Vec3 point(x, y, z);

				// for voxels the narrow-band classification has already
				// proven are well clear of the surface, skip the real
				// (possibly expensive) SDF call entirely and use a placeholder
				// that behaves identically to a real value of that magnitude
				// for every downstream use (early-out, thickness lookup,
				// boundary clamp)
				float containerDist;
				switch (narrowBand[idx]) {
				case NarrowBandClass::Outside:
					containerDist = classificationMargin + 1.0f;
					break;
				case NarrowBandClass::Inside:
					containerDist = -(classificationMargin + 1.0f);
					break;
				default:
					containerDist = con.sdf->compute_distance(point);
					break;
				}

				containerDistField[idx] = containerDist;
				if (containerDist > surfaceMargin) { // If it's more than surfaceMargin outside, don't bother
					scalarField[idx] = isoLevel + 1.0f;
					continue; // scalarField[idx] remains 9999.9f
				}

				// evaluate distance and radius function
				float localIsoLevel = isoLevel;

				if (thicknessFunction && thicknessSDF) {

					// "distance from container" thickness mode reuses the
					// container's own SDF; reuse the value just computed above
					// instead of paying for a second identical BVH+raycast query
					double rawDist = (thicknessSDF.get() == con.sdf.get())
						? containerDist
						: thicknessSDF->compute_distance(point);

					// convert to unsigned
					rawDist = std::abs(rawDist);

					localIsoLevel = static_cast<float>(
						thicknessFunction->estimate_radius(rawDist, startThickness, endThickness));
				}
				float d1 = 0.0f, d2 = 0.0f, d3 = 0.0f;
				Vec3 grad1, grad2, grad3;

				if(variedAnisotropy){
					
					// estimate the blended metric
					Eigen::Matrix3f M = blend_metric(
						point, anisotropySources, background, backgroundWeight);

					// use the estimated candidate number to find knn
					auto cand = kdtree->knn(
						point, Kn,
						[](const Vec3& a, const Vec3& b){
       						Vec3 v = b - a; return double(v.x*v.x + v.y*v.y + v.z*v.z);   // Euclidean
   						}
					);

					// re-rank them
					// first change the distance of each candidate to anisotropic
					for (auto& c : cand){
						c.second = aniso_distance_sq(M, point, seeds[c.first]);
					}
					std::partial_sort(
						cand.begin(), cand.begin() + 3, cand.end(),
						[](const auto& a, const auto& b){ return a.second < b.second; });
					
					d1 = std::sqrt((float)cand[0].second);
					d2 = std::sqrt((float)cand[1].second);
					d3 = std::sqrt((float)cand[2].second);
					Vec3 p1 = seeds[cand[0].first];  
					Vec3 p2 = seeds[cand[1].first];
					Vec3 p3 = seeds[cand[2].first];
					grad1 = aniso_distance_grad(M, point, p1, d1);
					grad2 = aniso_distance_grad(M, point, p2, d2);
					grad3 = aniso_distance_grad(M, point, p3, d3);
				}
				else{
					Vec3 localPt = point - center;

					Vec3 rotatedPt = Vec3(rot * Eigen::Vector3f{ localPt.x, localPt.y, localPt.z });

					// wrap the point
					//Vec3 wrapped(point.x / stretchX, point.y / stretchY, point.z / stretchZ);
					Vec3 wrapped(rotatedPt.x / stretchX, rotatedPt.y / stretchY, rotatedPt.z / stretchZ);

					// we need to find the two nearest seeds to the point, and compute the distance to the nearest seed, and the distance to the second nearest seed
					auto neighbors = kdtree->knn(wrapped, 3, [this](const Vec3& p1, const Vec3& p2) {
						// apply the formula (p - q)^T M (p - q), where M = diag(stretchX,stretchY,stretchZ)
						Vec3 v = p2 - p1;
						return (v.x * v.x) + (v.y * v.y) + (v.z * v.z);
					});

					// get_distances of three closest seeds
					d1 = (float)std::sqrt(neighbors[0].second);
					d2 = (float)std::sqrt(neighbors[1].second);
					d3 = (float)std::sqrt(neighbors[2].second);

					//float kf = 0.5f; // Smoothing radius (tunable)
					//float smoothd1 = smin(d1, d2, kf);
					//float smoothd2 = smin(d2, d3, kf);

					//float value = 0.0;
					//if (foam) {
					//	//value = (smoothd2 - smoothd1) + threshold * (d3 - d2);
					//	value = d2 - d1;
					//}
					//else {
					//	//value = (d3 - d1) + threshold * (smoothd2 - smoothd1);
					//	value = (d3 - d1) + threshold * (d2 - d1);
					//}

					//scalarField[idx] = value;
					
					// Retrieve the actual warped seed coordinates using the indices
					//Vec3 p1 = seeds[neighbors[0].first];
					//Vec3 p2 = seeds[neighbors[1].first];
					//Vec3 p3 = seeds[neighbors[2].first];
					Vec3 p1 = warpedSeeds[neighbors[0].first];
					Vec3 p2 = warpedSeeds[neighbors[1].first];
					Vec3 p3 = warpedSeeds[neighbors[2].first];

					// Lambda to calculate the exact analytical gradient of a single distance
					auto calc_grad = [&](const Vec3& p, float d) -> Vec3 {
						if (d < 1e-6f) return Vec3(0.0f, 0.0f, 0.0f); // Prevent division by zero at the exact seed center
						
						Vec3 local(
							(wrapped.x - p.x) / (d * stretchX),
							(wrapped.y - p.y) / (d * stretchY),
							(wrapped.z - p.z) / (d * stretchZ)
						);

						Vec3 res = Vec3(rot.transpose() * Eigen::Vector3f(local.x, local.y, local.z));

						return res;
					};

					// Calculate the individual distance gradients
					grad1 = calc_grad(p1, d1);
					grad2 = calc_grad(p2, d2);
					grad3 = calc_grad(p3, d3);
				}
				float value = 0.0f;
				Vec3 gradValue;

				// Combine the values and gradients analytically
				if (foam) {
					value = d2 - d1;
					gradValue = grad2 - grad1;
				}
				else {
					value = (d3 - d1) + (1.0 - threshold) * (d2 - d1);
					gradValue = (grad3 - grad1) + (grad2 - grad1) * (1.0 - threshold);
					//value = (d2 - d1) + threshold * (d3 - d1);
					//gradValue = (grad2 - grad1) + (grad3 - grad1) * threshold;
				}

				// normalize the gradient
				float gradMag = gradValue.norm();

				// this is the local raw value
				float localRaw = (gradMag > 1e-5f) ? (value / gradMag * 2.0f) : value;

				// shift the scalar field based on the local iso level, if no distance functions
				// it will keep the global iso level
				scalarField[idx] = localRaw - localIsoLevel + isoLevel;
			}
		}
	}

	const auto fillEnd = std::chrono::steady_clock::now();
	std::cout << "  [compute_scalar_field] fill loop (container SDF + seed kdtree): "
		<< std::chrono::duration<double>(fillEnd - fillStart).count() << " s" << std::endl;

	smooth_scalar_field();

	const auto smoothEnd = std::chrono::steady_clock::now();
	std::cout << "  [compute_scalar_field] smoothing: "
		<< std::chrono::duration<double>(smoothEnd - fillEnd).count() << " s" << std::endl;

	#pragma omp parallel for collapse(3)
	for (int i = 0; i < blockDims[0]; i++) {
		for (int j = 0; j < blockDims[1]; j++) {
			for (int k = 0; k < blockDims[2]; k++) {

				int idx = find_vertex_index(i, j, k);

				float mappedContainer = containerDistField[idx] + isoLevel;

				// Intersection: Take the maximum (most "Air-like") value
				scalarField[idx] = std::max(scalarField[idx], mappedContainer);
			}
		}
	}

	const auto clampEnd = std::chrono::steady_clock::now();
	std::cout << "  [compute_scalar_field] boundary clamp (cached SDF lookup): "
		<< std::chrono::duration<double>(clampEnd - smoothEnd).count() << " s" << std::endl;

	if (!isROI) {
		remove_isolated_islands();
	}

	seal_grid_boundaries();

	const auto cleanupEnd = std::chrono::steady_clock::now();
	std::cout << "  [compute_scalar_field] island removal + seal: "
		<< std::chrono::duration<double>(cleanupEnd - clampEnd).count() << " s" << std::endl;

	size_t solidVoxels = 0;
	for (float val : scalarField) {
		if (val < isoLevel) {
			solidVoxels++;
		}
	}
}

void GeneratorLewiner::smooth_scalar_field() {
	std::vector<float> smoothed = scalarField; // Copy

	// Simple 3x3x3 Box Blur
	#pragma omp parallel for collapse(3)
	for (int z = 1; z < blockDims[2] - 1; z++) {
		for (int y = 1; y < blockDims[1] - 1; y++) {
			for (int x = 1; x < blockDims[0] - 1; x++) {

				float sum = 0.0;
				int count = 0;

				size_t centerIdx = find_vertex_index(x, y, z);
				if (scalarField[centerIdx] > isoLevel + 0.5f) {
					// It's safely outside the container and the boundary margin, skip smoothing
					continue;
				}

				// Average neighbors
				for (int kz = -1; kz <= 1; kz++) {
					for (int ky = -1; ky <= 1; ky++) {
						for (int kx = -1; kx <= 1; kx++) {
							size_t idx = find_vertex_index(x + kx, y + ky, z + kz);
							sum += scalarField[idx];
							count++;
						}
					}
				}

				smoothed[find_vertex_index(x, y, z)] = sum / count;

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

void GeneratorLewiner::marching_cubes() {

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

	validate_topology();

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
	
	bounds[0] = aabb.pMin.x;
	bounds[1] = aabb.pMax.x;
	bounds[2] = aabb.pMin.y;
	bounds[3] = aabb.pMax.y;
	bounds[4] = aabb.pMin.z;
	bounds[5] = aabb.pMax.z;

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
// and to create a smoother mesh
//@param a: the first value
//@param b: the second value
//@param k: the smoothing factor, higher values create a sharper transition,
// lower values create a smoother transition
//@returns the smooth minimum between a and b
float GeneratorLewiner::smin(float a, float b, float k) {
	float h = std::max(k - std::abs(a - b), 0.0f) / k;
	return std::min(a, b) - h * h * h * k * (1.0f / 6.0f);
}

void GeneratorLewiner::update_steps() {
	stepX = (bounds[1] - bounds[0]) / (blockDims[0] - 1);
	stepY = (bounds[3] - bounds[2]) / (blockDims[1] - 1);
	stepZ = (bounds[5] - bounds[4]) / (blockDims[2] - 1);
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
	
	std::filesystem::path parent = std::filesystem::path(fileName).parent_path();
	std::string stlName = std::filesystem::path(fileName).stem().string();
	std::string parameterFileName = stlName + "_parameters.csv";
	std::filesystem::path parameterPath = parent / parameterFileName;

	export_parameters(parameterPath.string());
};

void GeneratorLewiner::validate_topology() {

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

    ImGui::ColorEdit4("Appearance", (float*)&color);  

    if (lockedCon) ImGui::Text("Container: %s", lockedCon->name.c_str());
    if (lockedGen) ImGui::Text("Generator: %s", lockedGen->name.c_str());

    // ---- editable parameters: locked while this scaffold is regenerating ----
    ImGui::BeginDisabled(mineBusy);

	ImGui::BeginTabBar("Items");

	thickness_properties();

	anisotropy_properties(globalSources);

	ImGui::EndTabBar();

    ImGui::SeparatorText("Parameters");
    ImGuiTableFlags flags = ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_BordersInnerV;
    if (ImGui::BeginTable("", 2, flags)) {
        if (lockedGen) {
            ImGui::TableNextRow();
            if (lockedGen->get_type() == ObjectType::RandomGeneratorType) {
                ImGui::TableNextColumn(); ImGui::Text("Random Seeds");
                ImGui::TableNextColumn(); ImGui::Text("%d", lockedGen->get_seeds().size());
            }
            else if (lockedGen->get_type() == ObjectType::PoissonGeneratorType) {
                ImGui::TableNextColumn(); ImGui::Text("Poisson 3D");
                ImGui::TableNextColumn();
                Poisson3D* dummy = static_cast<Poisson3D*>(lockedGen.get());
                if (dummy->is_uniform()) {
                    ImGui::Text("Radius (Uniform) %.4f", dummy->get_min_radius());
                } else {
                    ImGui::BeginTable("##", 2);
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn(); ImGui::Text("Rmin %4.f", dummy->get_min_radius());
                    ImGui::TableNextColumn(); ImGui::Text("Rmax %4.f", dummy->get_max_radius());
                    ImGui::EndTable();
                }
            }
        } else {
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::Text("Source");
            ImGui::TableNextColumn(); ImGui::TextColored(ImVec4(0.4f, 0.8f, 0.4f, 1.0f), "Loaded from CSV");
        }

        ImGui::EndTable();
    }

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
                    case 2: thicknessFunction = std::make_shared<ConstantRadiusFunction>();              break;
                    case 3: thicknessFunction = std::make_shared<RandomRadiusFunction>();                break;
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
            start_time = std::chrono::steady_clock::now();

            task->start([self, conShared, t = task]() {
                t->set_progress(0.00f);  self->compute_scalar_field(*conShared);
                t->set_progress(0.50f);  self->marching_cubes();           // CPU only — no GL
                t->set_progress(0.90f);  self->estimate_metrics(*conShared);
                t->set_progress(1.00f);
            }, this);                      
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

        // ---- completion (main thread: GL upload + log) ----
        if (task && task->poll(this)) {
            update_render();                               // GPU upload, GL thread

            auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - start_time);
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(3) << dur.count() / 1000.0 << " seconds!";
            logger->log(LogPriority::SUCCESS, "Updated Scaffold Successfully in " + oss.str());

            ImGui::CloseCurrentPopup(); 
        }

        ImGui::EndPopup();
    }
}

void GeneratorLewiner::thickness_properties(){

	if (ImGui::BeginTabItem("Thickness")){
	 	ImGui::RadioButton(
			"Apply Uniform Thickness", &selectedThicknessOption, 0);
    	ImGui::RadioButton(
			"Apply Varied Thickness", &selectedThicknessOption, 1);

		if (selectedThicknessOption == 0) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Thickness", &isoLevel);
		}
		else {
			ImGui::SetNextItemWidth(200); ImGui::InputFloat("Start Thickness", &startThickness);
			ImGui::SetNextItemWidth(200); ImGui::InputFloat("End Thickness", &endThickness);
			ImGui::SetNextItemWidth(200); ImGui::InputFloat("Transition Distance", &transitionDistance);

			ImGui::SeparatorText("Select Distance Function");
			ImGui::RadioButton("Distance From Plane", &selectedDist, 0);
			if (selectedDist == 0) {
				ImGui::SetNextItemWidth(200); ImGui::InputFloat3("Normal", distancePlaneNormal);
				ImGui::SetNextItemWidth(200); ImGui::InputFloat3("Center", distancePlaneCenter);
			}
			ImGui::RadioButton("Distance From Point", &selectedDist, 1);
			if (selectedDist == 1) {
				ImGui::SetNextItemWidth(200); ImGui::InputFloat3("Point", distancePoint);
			}
			ImGui::RadioButton("Distance From Container", &selectedDist, 2);

			ImGui::SeparatorText("Select Radius Function");
			ImGui::RadioButton("Linear", &selectedFunc, 0);
			ImGui::RadioButton("Quadratic", &selectedFunc, 1);
			ImGui::RadioButton("Constant", &selectedFunc, 2);
			ImGui::RadioButton("Random", &selectedFunc, 3);
		}

		ImGui::InputFloat("Voxel Size", &voxelSize);
		ImGui::SliderFloat("##Openess", &threshold, 0.0f, 1.0f, "%.3f");
		ImGui::EndTabItem();
	}
}

void GeneratorLewiner::anisotropy_properties(
	std::vector<std::shared_ptr<AnisotropySource>>& globalSources
){

	static int selectedIdx = -1;
	if(ImGui::BeginTabItem("Anisotropy")){
		
		ImGui::SeparatorText("Global Background");
		ImGui::InputFloat("Stretch X", &stretchX, 0.01f, 5.0f, "%.3f");
		ImGui::InputFloat("Stretch Y", &stretchY, 0.01f, 5.0f, "%.3f");
		ImGui::InputFloat("Stretch Z", &stretchZ, 0.01f, 5.0f, "%.3f");
		ImGui::InputFloat("Angle", &anisotropyAngle, 0.01f, 10.0f, "%.4f");
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
		ImGui::SliderFloat("Background Weight", &backgroundWeight, 0.01f, 1.0f, "%.3f");

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

	// get the domain volume to update the porosity
	float domainVolume = container.get_volume();
	porosity = (1 - volume / domainVolume);

	surfaceToVolume = surfaceArea / volume;
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
		ImGui::TableNextColumn(); ImGui::Text("Surface to Volume Ratio (1/mm)");
		ImGui::TableNextColumn(); ImGui::Text("%.4f", surfaceToVolume);

		ImGui::TableNextRow();
		ImGui::TableNextColumn(); ImGui::Text("Porosity (%%)");
		ImGui::TableNextColumn(); ImGui::Text("%.4f", porosity * 100.0f);

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

void GeneratorLewiner::set_options_from_factory(int distOption, int distFunc, int thicknessOption, float voxSize) {
	selectedDist = distOption;
	selectedFunc = distFunc;
	selectedThicknessOption = thicknessOption;
	voxelSize = voxSize;
};

void GeneratorLewiner::set_distance_plane_options(Vec3 center, Vec3 normal) {
	distancePlaneCenter = center;
	distancePlaneNormal = normal;
};

void GeneratorLewiner::set_distance_point_options(Vec3 point) {
	distancePoint = point;
};

void GeneratorLewiner::set_resolution(const std::array<int, 3>& newResolution) {
	blockDims = newResolution;
};

void GeneratorLewiner::set_bounds(const std::array<float, 6>& newBounds) {
	bounds = newBounds;
	stepX = (bounds[1] - bounds[0]) / (blockDims[0] - 1);
	stepY = (bounds[3] - bounds[2]) / (blockDims[1] - 1);
	stepZ = (bounds[5] - bounds[4]) / (blockDims[2] - 1);
};

void GeneratorLewiner::set_seeds(const std::vector<Vec3>& newSeeds) {
	seeds = newSeeds;
};

void GeneratorLewiner::set_stretch(float newStretchX, float newStretchY, float newStretchZ) {
	this->stretchX = newStretchX;
	this->stretchY = newStretchY;
	this->stretchZ = newStretchZ;
};

void GeneratorLewiner::set_thickness(float newThickness) {
	isoLevel = newThickness;
};

std::array<float, 6> GeneratorLewiner::get_bounds() const {
	return bounds;
};

Aabb GeneratorLewiner::get_aabb() const { return aabb; };

void GeneratorLewiner::estimate_local_thickness(
	float voxelSize, std::array<float, 6>& blockBounds, bool separation) {

	// we should ensure that the inserted blockbounds are clipped inside the aligned bounding box
	blockBounds[0] = std::max(blockBounds[0], aabb.pMin.x); // Min X
	blockBounds[1] = std::min(blockBounds[1], aabb.pMax.x); // Max X
	blockBounds[2] = std::max(blockBounds[2], aabb.pMin.y); // Min Y
	blockBounds[3] = std::min(blockBounds[3], aabb.pMax.y); // Max Y
	blockBounds[4] = std::max(blockBounds[4], aabb.pMin.z); // Min Z
	blockBounds[5] = std::min(blockBounds[5], aabb.pMax.z); // Max Z

	if (blockBounds[0] >= blockBounds[1] ||
		blockBounds[2] >= blockBounds[3] ||
		blockBounds[4] >= blockBounds[5]) {
		std::cerr << "Error: The requested bounds do not overlap with the mesh AABB." << std::endl;
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

		// get the position of the voxel in the grid
		int rz = static_cast<int>(ridgeIdx % nz);
		int ry = static_cast<int>((ridgeIdx / nz) % ny);
		int rx = static_cast<int>(ridgeIdx / (ny * nz));

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
	float T_min = 1.5f; // Minimum thickness threshold to avoid surface noise

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

				// Only accumulate if the voxel is physically inside the container ROI
				if (isInsideROI) {
					float t = thicknessMap[idx];
					if (t >= T_min) {
						totalThicknessSum += t;
						solidVoxelCount++;
					}
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

				if (isInsideROI) {
					float t = thicknessMap[idx];
					if (t >= T_min) {
						float diff = t - meanThicknessVoxels;
						deviationSum += (diff * diff);
					}
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
		logger->log(LogPriority::SUCCESS, "Estimated Local Separation!");
	}
	else {
		thicknessVersion = meshVersion;
		logger->log(LogPriority::SUCCESS, "Estimated Local Thickness!");
	}
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
	uint32_t fileVersion = 1;
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
	out.write(reinterpret_cast<const char*>(&foam), sizeof(foam));
	out.write(reinterpret_cast<const char*>(&voxelSize), sizeof(voxelSize));
	out.write(reinterpret_cast<const char*>(&stretchX), sizeof(stretchX));
	out.write(reinterpret_cast<const char*>(&stretchY), sizeof(stretchY));
	out.write(reinterpret_cast<const char*>(&stretchZ), sizeof(stretchZ));
	out.write(reinterpret_cast<const char*>(&anisotropyAngle), sizeof(anisotropyAngle));

	// FIXED: Explicitly write components of Vec3 to avoid struct padding offsets
	out.write(reinterpret_cast<const char*>(&anisotropyVec.x), sizeof(float));
	out.write(reinterpret_cast<const char*>(&anisotropyVec.y), sizeof(float));
	out.write(reinterpret_cast<const char*>(&anisotropyVec.z), sizeof(float));

	out.write(reinterpret_cast<const char*>(&renderMode), sizeof(renderMode));

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

	out.close();
	logger->log(LogPriority::SUCCESS, "Procedural scaffold serialized smoothly to " + fileName);
}

bool GeneratorLewiner::load_scaf(const std::string& fileName,
	std::vector<std::shared_ptr<IContainer>>& containerList,
	std::vector<std::shared_ptr<InterfaceSeedGenerator>>& generatorList,
	std::vector<std::shared_ptr<AnisotropySource>>& globalSources
) {
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
	if (fileVersion != 1) {
		logger->log(LogPriority::ERROR, "Unsupported file version.");
		return false;
	}

	// Read Core Parameters Block
	in.read(reinterpret_cast<char*>(&selectedThicknessOption), sizeof(selectedThicknessOption));
	in.read(reinterpret_cast<char*>(&selectedDist), sizeof(selectedDist));
	in.read(reinterpret_cast<char*>(&selectedFunc), sizeof(selectedFunc));
	in.read(reinterpret_cast<char*>(&startThickness), sizeof(startThickness));
	in.read(reinterpret_cast<char*>(&endThickness), sizeof(endThickness));
	in.read(reinterpret_cast<char*>(&transitionDistance), sizeof(transitionDistance));
	in.read(reinterpret_cast<char*>(&threshold), sizeof(threshold));
	in.read(reinterpret_cast<char*>(&isoLevel), sizeof(isoLevel));
	in.read(reinterpret_cast<char*>(&foam), sizeof(foam));
	in.read(reinterpret_cast<char*>(&voxelSize), sizeof(voxelSize));
	in.read(reinterpret_cast<char*>(&stretchX), sizeof(stretchX));
	in.read(reinterpret_cast<char*>(&stretchY), sizeof(stretchY));
	in.read(reinterpret_cast<char*>(&stretchZ), sizeof(stretchZ));
	in.read(reinterpret_cast<char*>(&anisotropyAngle), sizeof(anisotropyAngle));

	// Explicit element parsing for Vec3 elements
	in.read(reinterpret_cast<char*>(&anisotropyVec.x), sizeof(float));
	in.read(reinterpret_cast<char*>(&anisotropyVec.y), sizeof(float));
	in.read(reinterpret_cast<char*>(&anisotropyVec.z), sizeof(float));

	in.read(reinterpret_cast<char*>(&renderMode), sizeof(renderMode));

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
		
		// Read angle if you kept it in the export
		in.read(reinterpret_cast<char*>(&src->angle), sizeof(float));
		
		src->name = "Anisotropy Source " + std::to_string(i + 1);
		
		// Rebuild the math tensor and the visual geometry
		src->update_metric();
		src->update_model();
		
		anisotropySources.push_back(src);
		globalSources.push_back(src);
	}

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
		if (std::filesystem::exists(pathStr)) {
			rebuiltContainer = std::make_shared<AbstractContainer>(pathStr);
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

		if (this->renderMode && !this->seeds.empty()) {
			rebuiltGen->update_model();
		}

		generatorList.push_back(rebuiltGen);
		this->generator = rebuiltGen;
	}

	// Read Raw Grid Array Block
	uint64_t voxelCount = 0;
	in.read(reinterpret_cast<char*>(&voxelCount), sizeof(voxelCount));
	scalarField.resize(voxelCount);
	if (voxelCount > 0) {
		in.read(reinterpret_cast<char*>(scalarField.data()), sizeof(float) * voxelCount);
	}

	in.close();

	// Re-verify stream status right before triggering reconstruction
	if (in.fail()) {
		logger->log(LogPriority::ERROR, "SCAF Stream parsing collapsed during binary data array extraction phases!");
		return false;
	}

	// 12. Run final level-set boundary tracking reconstruction
	marching_cubes();

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

void GeneratorLewiner::estimate_anisotropy(int daDirectionNr, int linesPerDirection, int mode, ROI* roi) {

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

	float boundsVoxel[6];
	if (roi) {
		std::array<float, 6> reqBounds = roi->get_bounds();
		int min_i = std::max(0, static_cast<int>(std::floor((reqBounds[0] - bounds[0]) / stepX)));
		int max_i = std::min(blockDims[0] - 1, static_cast<int>(std::ceil((reqBounds[1] - bounds[0]) / stepX)));
		int min_j = std::max(0, static_cast<int>(std::floor((reqBounds[2] - bounds[2]) / stepY)));
		int max_j = std::min(blockDims[1] - 1, static_cast<int>(std::ceil((reqBounds[3] - bounds[2]) / stepY)));
		int min_k = std::max(0, static_cast<int>(std::floor((reqBounds[4] - bounds[4]) / stepZ)));
		int max_k = std::min(blockDims[2] - 1, static_cast<int>(std::ceil((reqBounds[5] - bounds[4]) / stepZ)));

		boundsVoxel[0] = static_cast<float>(min_i); boundsVoxel[1] = static_cast<float>(max_i);
		boundsVoxel[2] = static_cast<float>(min_j); boundsVoxel[3] = static_cast<float>(max_j);
		boundsVoxel[4] = static_cast<float>(min_k); boundsVoxel[5] = static_cast<float>(max_k);
	}
	else {
		// Fallback to full grid voxel dimensions
		boundsVoxel[0] = 0.0f; boundsVoxel[1] = static_cast<float>(blockDims[0]);
		boundsVoxel[2] = 0.0f; boundsVoxel[3] = static_cast<float>(blockDims[1]);
		boundsVoxel[4] = 0.0f; boundsVoxel[5] = static_cast<float>(blockDims[2]);
	}

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
						localBoxLen += (tMax - tMin);
						float startT = tMin + dist(rng) * rayStepSize;
						long samples = static_cast<long>(std::ceil((tMax - startT) / rayStepSize));
						bool previousPhase = false;

						for (long s = 0; s < samples; s++) {
							Vec3 pt = rayOrigin + d * (startT + s * rayStepSize);
							long vx = std::clamp<long>(static_cast<long>(pt.x), static_cast<long>(boundsVoxel[0]), static_cast<long>(boundsVoxel[1]) - 1);
							long vy = std::clamp<long>(static_cast<long>(pt.y), static_cast<long>(boundsVoxel[2]), static_cast<long>(boundsVoxel[3]) - 1);
							long vz = std::clamp<long>(static_cast<long>(pt.z), static_cast<long>(boundsVoxel[4]), static_cast<long>(boundsVoxel[5]) - 1);

							bool currentPhase = (get_data(vx, vy, vz) < isoLevel);
							if (currentPhase != previousPhase) { localTransitions++; }
							previousPhase = currentPhase;
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
			float lenMax = (lambdaMin > 1e-9f) ? (1.0f / std::sqrt(lambdaMin)) : 0.0f;
			float lenMin = (lambdaMax > 1e-9f) ? (1.0f / std::sqrt(lambdaMax)) : 0.0f;
			anisotropyDegree = (lenMax > 1e-9f) ? (lenMax / lenMin) : 0.0f;
			break;
		}
		case 1: {
			float lenMax = (lambdaMin > 1e-9f) ? (1.0f / std::sqrt(lambdaMin)) : 0.0f;
			float lenMin = (lambdaMax > 1e-9f) ? (1.0f / std::sqrt(lambdaMax)) : 0.0f;
			anisotropyDegree = (lenMax > 1e-9f) ? (1.0f - (lenMax / lenMin)) : 0.0f;
			break;
			}
		case 2: {
			anisotropyDegree = (lambdaMax > 1e-9f) ? static_cast<float>(lambdaMin / lambdaMax) : 0.0f;
			break;
		}
		case 3: {
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
//   1) Removed the spurious "/ 2.0".  Tb.N = S_V/2 and S_V = 2*P_L  ==>  Tb.N = P_L.
//      The two factors of 2 cancel, so isotropically-averaged transitions/length
//      already IS Tb.N. No extra division.
//   2) Replaced the sqrt(3) point-sampling with an exact 3D DDA (Amanatides & Woo)
//      voxel walk, so no thin strut is ever skipped and the count is step-size-free.
//   3) Baseline phase is taken from the FIRST voxel, removing the phantom entry transition.
//   4) Unit conversion unchanged: valid IFF stepX is the physical voxel size in mm/voxel.

void GeneratorLewiner::estimate_trabecular_number(int formula, int daDirectionNr, int linesPerDirection, ROI* roi) {

	// this is if the formula is set to use the Tb.N = (BV/TV) / Tb.Th
	if (formula == 1 && trabecularNrVersion < thicknessVersion) {
		float bvtv = 1.0f - porosity;
		trabecularNr = bvtv / localThickness;
		logger->log(LogPriority::SUCCESS, "Estimated Trabecular Number with Tb.N = (BV/TV) / Tb.Th");
		logger->log(LogPriority::WARNING, "The Estimated Trabecular Number is measured with the previously estimated local thickness. If this was estimated inside a ROI, the estimation is wrong. Try creating the scaffold inside the ROI first.");
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

	float boundsVoxel[6];
	if (roi) {
		std::array<float, 6> reqBounds = roi->get_bounds();
		int min_i = std::max(0, static_cast<int>(std::floor((reqBounds[0] - bounds[0]) / stepX)));
		int max_i = std::min(blockDims[0] - 1, static_cast<int>(std::ceil((reqBounds[1] - bounds[0]) / stepX)));
		int min_j = std::max(0, static_cast<int>(std::floor((reqBounds[2] - bounds[2]) / stepY)));
		int max_j = std::min(blockDims[1] - 1, static_cast<int>(std::ceil((reqBounds[3] - bounds[2]) / stepY)));
		int min_k = std::max(0, static_cast<int>(std::floor((reqBounds[4] - bounds[4]) / stepZ)));
		int max_k = std::min(blockDims[2] - 1, static_cast<int>(std::ceil((reqBounds[5] - bounds[4]) / stepZ)));

		boundsVoxel[0] = static_cast<float>(min_i); boundsVoxel[1] = static_cast<float>(max_i);
		boundsVoxel[2] = static_cast<float>(min_j); boundsVoxel[3] = static_cast<float>(max_j);
		boundsVoxel[4] = static_cast<float>(min_k); boundsVoxel[5] = static_cast<float>(max_k);
	}
	else {
		// Fallback to full grid voxel dimensions
		boundsVoxel[0] = 0.0f; boundsVoxel[1] = static_cast<float>(blockDims[0]);
		boundsVoxel[2] = 0.0f; boundsVoxel[3] = static_cast<float>(blockDims[1]);
		boundsVoxel[4] = 0.0f; boundsVoxel[5] = static_cast<float>(blockDims[2]);
	}

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
					localBoxLen += (tExit - tEntry);

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

					bool previousPhase = (get_data(ix, iy, iz) < isoLevel);

					while (true) {
						float tNext = std::min(tMaxX, std::min(tMaxY, tMaxZ));
						if (tNext > tExit) break;

						if (tMaxX <= tMaxY && tMaxX <= tMaxZ) { ix += sx; tMaxX += tDeltaX; }
						else if (tMaxY <= tMaxZ) { iy += sy; tMaxY += tDeltaY; }
						else { iz += sz; tMaxZ += tDeltaZ; }

						// Break if out of the scoped boundary limits
						if (ix < static_cast<int>(boundsVoxel[0]) || ix >= static_cast<int>(boundsVoxel[1]) ||
							iy < static_cast<int>(boundsVoxel[2]) || iy >= static_cast<int>(boundsVoxel[3]) ||
							iz < static_cast<int>(boundsVoxel[4]) || iz >= static_cast<int>(boundsVoxel[5])) break;

						bool currentPhase = (get_data(ix, iy, iz) < isoLevel);
						if (currentPhase != previousPhase) localTransitions++;
						previousPhase = currentPhase;
					}
				}
			}

			globalTransitions += localTransitions;
			globalBoxLen += localBoxLen;
		}
	}

	// final values
	if (globalBoxLen > 0.0) {
		double PL_voxels = static_cast<double>(globalTransitions) / globalBoxLen; // intersections / voxel-length
		double TbN_voxels = PL_voxels;            // Tb.N = P_L  (NO extra /2)
		trabecularNr = static_cast<float>(TbN_voxels / stepX); // 1/voxel -> 1/mm  (stepX must be mm/voxel)
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

	//std::cout << eulerCharacteristic << std::endl;
	//std::cout << domainVolume << std::endl;

	// estimate the genus
	float genus = 1.0f - (static_cast<float>(eulerCharacteristic) / 2.0f);

	// connectivity density is genus / domain volume
	connectivityDensity = genus / domainVolume;
	
	connectivityVersion = meshVersion;

	logger->log(LogPriority::SUCCESS, "Estimated Connectivity Density!");
};

void GeneratorLewiner::estimate_connectivity_network() {

	logger->log(LogPriority::ERROR, "Not implemented yet!");
	return;

};

void GeneratorLewiner::apply_taubin_smooth(int iter, float lambda, float mu) {

	int vertNr = meshVertices.size();

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
	fout << "Porosity, Volume, TotalSurface, SurfaceToVolume, Connectivity Density, Local Thickness, Local Thickness Std, Local Separation, Local Separation Std, trabecular Nr, Anisotropy, Tortuosity\n";

	// pass values
	fout << porosity << "," << volume << "," << surfaceArea << "," << surfaceToVolume << "," << connectivityDensity << "," << localThickness << "," << localThicknessStd << "," << localSeparation << "," << localSeparationStd << "," << trabecularNr << "," << anisotropyDegree << "," << tortuosity << "\n";

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
		else if (key == "Volume") { loadedMetrics.volume = val; }
		else if (key == "TotalSurface") { loadedMetrics.totalSurface = val; }
		else if (key == "SurfaceToVolume") { loadedMetrics.surfToVol = val; }
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
	volume = loadedMetrics.volume;
	surfaceArea = loadedMetrics.totalSurface;
	surfaceToVolume = loadedMetrics.surfToVol;
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
	cfg.stretchX = stretchX; cfg.stretchY = stretchY; cfg.stretchZ = stretchZ;
	cfg.anisotropyAngle = anisotropyAngle;
	cfg.dirX = anisotropyVec.x; cfg.dirY = anisotropyVec.y; cfg.dirZ = anisotropyVec.z;

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

	// Write Fixed Header (Expanded)
	fout << "ThicknessOption,UniformThickness,StartThickness,EndThickness,DistFunction,RadFunction,"
		<< "PlaneOriginX,PlaneOriginY,PlaneOriginZ,PlaneNormalX,PlaneNormalY,PlaneNormalZ,PointX,PointY,PointZ,"
		<< "GeneratorType,SeedNr,MinRadius,MaxRadius,Openess,StretchX,StretchY,StretchZ,"
		<< "AnisotropyAngle,DirX,DirY,DirZ\n";

	// Write Fixed Data Row
	fout << cfg.thicknessOption << "," << cfg.uniformThickness << "," << cfg.startThickness << "," << cfg.endThickness << "," << cfg.distFunction << "," << cfg.radFunction << ","
		<< cfg.planeOriginX << "," << cfg.planeOriginY << "," << cfg.planeOriginZ << ","
		<< cfg.planeNormalX << "," << cfg.planeNormalY << "," << cfg.planeNormalZ << ","
		<< cfg.pointX << "," << cfg.pointY << "," << cfg.pointZ << ","
		<< cfg.generatorType << "," << cfg.seedNr << "," << cfg.minRadius << "," << cfg.maxRadius << ","
		<< cfg.openess << "," << cfg.stretchX << "," << cfg.stretchY << "," << cfg.stretchZ << ","
		<< cfg.anisotropyAngle << "," << cfg.dirX << "," << cfg.dirY << "," << cfg.dirZ << "\n";

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
		this->foam,
		this->renderMode
	);

	roiScaffold->set_stretch(this->stretchX, this->stretchY, this->stretchZ);
	roiScaffold->anisotropyAngle = this->anisotropyAngle;
	roiScaffold->anisotropyVec = this->anisotropyVec;
	roiScaffold->anisotropySources = this->anisotropySources;
	roiScaffold->set_thickness_functions(this->thicknessSDF, this->thicknessFunction, this->startThickness, this->endThickness, this->transitionDistance);
	roiScaffold->isROI = true;
	roiScaffold->name = this->name + " (ROI)";

	// 5. Generate the mesh using the parent's container
	std::shared_ptr<IContainer> parentCon = this->container.lock();
	if (parentCon) {
		roiScaffold->compute_scalar_field(*parentCon);
		roiScaffold->marching_cubes();
		roiScaffold->estimate_metrics(*parentCon);
	}

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
	foam = 0;
	voxelSize = 0.05f;

	// for thickness function
	selectedThicknessOption = 0;
	distancePlaneNormal = { 0.0f, 0.0f, 1.0f };
	distancePlaneCenter = { 0.0f, 0.0f, 0.0f };
	distancePoint = { 0.0f, 0.0f, 0.0f };
	transitionDistance = 10.0f;

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
		
		ImGuiTabBarFlags tab_bar_flags = ImGuiTabBarFlags_None;
		if (ImGui::BeginTabBar("MyTabBar", tab_bar_flags)){

			main_options(containers, generators);

			thickness_options();

			anisotropy_options(anisoSources);

			ImGui::EndTabBar();
		}

		ImGui::Separator();

        const bool anyBusy = task->get_running();
		const bool mineBusy = task->is_running_for(this);

        // ---------- Generate ----------
        ImGui::BeginDisabled(anyBusy);
        if (ImGui::Button("Generate")) {
            if (!lockedCon) warningFlashTimer1 = 1.5f;
            if (!lockedGen) warningFlashTimer2 = 1.5f;

            if (lockedCon && lockedGen) {
                auto conShared = lockedCon;                 
                auto seeds     = lockedGen->get_seeds();    
                auto bds       = conShared->compute_bounds();

                std::array<float, 6> bounds = {
                    bds.xMin, bds.xMax, bds.yMin, bds.yMax, bds.zMin, bds.zMax
                };
                resolution = {
                    static_cast<int>(std::ceil((bds.xMax - bds.xMin) / voxelSize)) + 1,
                    static_cast<int>(std::ceil((bds.yMax - bds.yMin) / voxelSize)) + 1,
                    static_cast<int>(std::ceil((bds.zMax - bds.zMin) / voxelSize)) + 1
                };

                auto scaffold = std::make_unique<GeneratorLewiner>(
                    seeds, bounds, resolution, logger, openess, thickness, foam);

                if (selectedThicknessOption == 1) {
                    switch (selectedFunc) {
                        case 0: thicknessRadFunc = std::make_shared<LinearFunction>(transitionDistance);    break;
                        case 1: thicknessRadFunc = std::make_shared<QuadraticFunction>(transitionDistance); break;
                        case 2: thicknessRadFunc = std::make_shared<ConstantRadiusFunction>();              break;
                        case 3: thicknessRadFunc = std::make_shared<RandomRadiusFunction>();                break;
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

                scaffold->set_options_from_factory(selectedDist, selectedFunc, selectedThicknessOption, voxelSize);
                scaffold->set_thickness_functions(
					thicknessSDF, thicknessRadFunc,
                    startThickness, endThickness, transitionDistance);
                scaffold->set_stretch(stretchX, stretchY, stretchZ);
                scaffold->anisotropyAngle = anisotropyAngle;
                scaffold->anisotropyVec   = anisotropyVec;
				scaffold->backgroundWeight  = backgroundWeight;
				scaffold->anisotropySources = anisoSources;
                scaffold->container = lockedCon;
                scaffold->generator = lockedGen;
                if (foam == 1) scaffold->foam = true;

                genContainerName = lockedCon->name;   // capture identity now, for the log later
                genGeneratorName = lockedGen->name;

                GeneratorLewiner* raw = scaffold.get();
                pendingScaffold = std::move(scaffold);

                start_time = std::chrono::steady_clock::now();
                task->start([raw, conShared, t = task]() {
                    t->set_progress(0.00f);  
					raw->compute_scalar_field(*conShared);
                    t->set_progress(0.50f);  
					raw->marching_cubes();
                    t->set_progress(0.90f);
					raw->estimate_metrics(*conShared);
                    t->set_progress(1.00f);
                }, this);
            }
        }
        ImGui::EndDisabled();

		ImGui::SameLine();

		if (ImGui::Button("Cancel")) {
			showPopup = false;
		};

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
		ImGui::EndPopup();
	}
};

void ScaffoldFactory::thickness_options(){
	if (ImGui::BeginTabItem("Thickness")){
		ImGui::RadioButton("Apply Uniform Thickness", &selectedThicknessOption, 0);
		ImGui::RadioButton("Apply Varied Thickness", &selectedThicknessOption, 1);

		if (selectedThicknessOption == 0) {
			ImGui::InputFloat("Thickness", &thickness, 0.001f, 1.0f);
		}
		else {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Start Thickness", &startThickness);
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("End Thickness", &endThickness);
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat("Transition Distance", &transitionDistance);

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

		ImGui::SliderFloat("Openess", &openess, 0.0f, 1.0f, "%.3f");
		ImGui::EndTabItem();
	}
};

void ScaffoldFactory::anisotropy_options(
	std::vector<std::shared_ptr<AnisotropySource>>& globalSources
){

	if(ImGui::BeginTabItem("Anisotropy")){
		ImGui::SeparatorText("Global Background");
		ImGui::InputFloat("Stretch X", &stretchX, 0.01f, 5.0f, "%.3f");
		ImGui::InputFloat("Stretch Y", &stretchY, 0.01f, 5.0f, "%.3f");
		ImGui::InputFloat("Stretch Z", &stretchZ, 0.01f, 5.0f, "%.3f");
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
		
		ImGui::SliderFloat("Background Weight", &backgroundWeight, 0.01f, 1.0f, "%.3f");

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

		ImGui::SeparatorText("Mode");
		ImGui::RadioButton("Porous", &foam, 0);
		ImGui::RadioButton("Foam", &foam, 1);
		ImGui::InputFloat("Voxel Size", &voxelSize);
		
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
		
		ImGui::EndTabItem();
	};
};
