#include "GeneratorLewiner.h"
#include "LookUpTable.h"
#include "Math/Vec.h"
#include "Math/Kdtree.h"
#include <chrono>
#include <memory>
#include <random>
#include <iostream>
#include <map>
#include <algorithm>
#include <omp.h>
#include <Eigen/Dense>

using namespace std::chrono_literals;

GeneratorLewiner::GeneratorLewiner(
	const std::vector<Vec3>& seeds,
	const std::array<float, 6>& bounds,
	const std::array<int, 3>& dims,
	const float threshold,
	const float isoLevel,
	const int foam) : seeds(seeds), bounds(bounds), blockDims(dims), threshold(threshold), isoLevel(isoLevel), foam(foam) {

	// setup opengl
	_setup_mesh();

	_setup_edges();

	aabb.pMin = Vec3(bounds[0], bounds[2], bounds[4]);
	aabb.pMax = Vec3(bounds[1], bounds[3], bounds[5]);
};


void GeneratorLewiner::compute_scalar_field(const IContainer& con) {

	// use also the domain volume
	domainVolume = con.get_volume();

	update_steps();

	// first populate the kdtree with the seeds
	std::unique_ptr<Kdtree>kdtree = std::make_unique<Kdtree>(seeds);

	// update scalar field
	scalarField.clear();

	// it will be more efficient to reserve the size of the scalar field vector before filling it, 
	// since we know the size of the grid in advance

	int totalVoxels = blockDims[0] * blockDims[1] * blockDims[2];

	//scalarField.reserve(totalVoxels);
	scalarField.resize(totalVoxels, 9999.9);

	#pragma omp parallel for collapse(3)
	// create the grid points based on the bounds and block dimensions
	for (int i{ 0 }; i < blockDims[0]; i++) {
		for (int j{ 0 }; j < blockDims[1]; j++) {
			for (int k{ 0 }; k < blockDims[2]; k++) {

				int idx = find_vertex_index(i, j, k);
				float x = bounds[0] + i * stepX;
				float y = bounds[2] + j * stepY;
				float z = bounds[4] + k * stepZ;
				Vec3 point(x, y, z);

				float containerDist = con.sdf->compute_distance(point);
				if (containerDist > 2.0f) { // If it's more than 2mm outside, don't bother
					continue; // scalarField[idx] remains 9999.9f
				}

				// we need to find the two nearest seeds to the point, and compute the distance to the nearest seed, and the distance to the second nearest seed
				auto neighbors = kdtree->knn(point, 3, [](const Vec3& p1, const Vec3& p2) {
					double dx = p1[0] - p2[0];
					double dy = p1[1] - p2[1];
					double dz = p1[2] - p2[2];
					return dx * dx + dy * dy + dz * dz;
				});

				// get_distances of three closest seeds
				float d1 = std::sqrt(neighbors[0].second);
				float d2 = std::sqrt(neighbors[1].second);
				float d3 = std::sqrt(neighbors[2].second);

				float kf = 0.5f; // Smoothing radius (tunable)
				float smoothd1 = smin(d1, d2, kf);
				float smoothd2 = smin(d2, d3, kf);

				float value = 0.0;
				if (foam) {
					//value = (smoothd2 - smoothd1) + threshold * (d3 - d2);
					value = d2 - d1;
				}
				else {
					value = (d3 - d1) + threshold * (smoothd2 - smoothd1);
				}

				scalarField[idx] = value;
			}
		}
	}

	smooth_scalar_field();

	#pragma omp parallel for collapse(3)
	for (int i = 0; i < blockDims[0]; i++) {
		for (int j = 0; j < blockDims[1]; j++) {
			for (int k = 0; k < blockDims[2]; k++) {

				int idx = find_vertex_index(i, j, k);

				float x = bounds[0] + i * stepX;
				float y = bounds[2] + j * stepY;
				float z = bounds[4] + k * stepZ;
				Vec3 point(x, y, z);

				float containerDist = con.sdf->compute_distance(point);
				float mappedContainer = containerDist + isoLevel;

				// Intersection: Take the maximum (most "Air-like") value
				scalarField[idx] = std::max(scalarField[idx], mappedContainer);
			}
		}
	}

	remove_isolated_islands();

	seal_grid_boundaries();

	int solidVoxels = 0;
	for (float val : scalarField) {
		if (val < isoLevel) {
			solidVoxels++;
		}
	}
	std::cout << "finished grids" << std::endl;

	std::cout << "scalar field size: " << scalarField.size() << std::endl;
	//volumeFraction = (float)solidVoxels / (float)totalVoxels;
	//porosity = 1.0 - volumeFraction;

	//volume = volumeFraction * (paddedBounds[1] - paddedBounds[0]) * (paddedBounds[3] - paddedBounds[2]) * (paddedBounds[5] - paddedBounds[4]);

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

				// Average neighbors
				for (int kz = -1; kz <= 1; kz++) {
					for (int ky = -1; ky <= 1; ky++) {
						for (int kx = -1; kx <= 1; kx++) {
							int idx = find_vertex_index(x + kx, y + ky, z + kz);
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
	int nx = blockDims[0];
	int ny = blockDims[1];
	int nz = blockDims[2];
	int totalVoxels = nx * ny * nz;

	std::vector<bool> visited(totalVoxels, false);
	std::vector<std::vector<int>> islands;

	// 1. Find all connected components
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {

				int idx = find_vertex_index(i, j, k);

				// If voxel is solid and hasn't been grouped into an island yet
				if (scalarField[idx] < isoLevel && !visited[idx]) {

					std::vector<int> currentIsland;
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
								int nIdx = find_vertex_index(n.x, n.y, n.z);
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
			for (int idx : islands[i]) {
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

	int totalVoxels = blockDims[0] * blockDims[1] * blockDims[2];
	x_verts.assign(totalVoxels, -1);
	y_verts.assign(totalVoxels, -1);
	z_verts.assign(totalVoxels, -1);
	meshVertices.resize(totalVoxels * 3);
	meshTriangles.resize(totalVoxels * 5); 
	vertexCount = 0;
	triangleCount = 0;

	compute_intersection_points();
	for (int k = 0; k < blockDims[2] - 1; k++) {
		for (int j = 0; j < blockDims[1] - 1; j++) {
			for (int i = 0; i < blockDims[0] - 1; i++) {
				
				int lut_entry = 0;

				float _cube[8];

				for (int p = 0; p < 8; ++p) {

					// ^ XOR this line is the 

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

	// update to exact size
	meshVertices.resize(vertexCount);
	meshTriangles.resize(triangleCount);

	// update opengl objects
	_update_render();

	validate_topology();
};

void GeneratorLewiner::compute_intersection_points() {
	
	#pragma omp parallel for collapse(3)
	for (int i = 0; i < blockDims[0]; i++) {
		for (int j = 0; j < blockDims[1]; j++) {
			for (int k = 0; k < blockDims[2]; k++) {

				// get current voxel index
				int currentIdx = find_vertex_index(i, j, k);

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
					int vIdx = vertexCount.fetch_add(1, std::memory_order_relaxed);
					
					// add to vector of Vertices
					meshVertices[vIdx] = add_x_vertex(i, j, k, val0, val1);
					
					// cash current to x vertices
					x_verts[currentIdx] = vIdx;
				}

				if ((val0 < 0 && val2 >= 0) || (val0>=0 && val2 < 0)) {
					int vIdx = vertexCount.fetch_add(1, std::memory_order_relaxed);
					meshVertices[vIdx] = add_y_vertex(i, j, k, val0, val2);
					y_verts[currentIdx] = vIdx;
				}

				if ((val0 < 0 && val3 >= 0) || (val0 >= 0 && val3 < 0)) {
					int vIdx = vertexCount.fetch_add(1, std::memory_order_relaxed);
					meshVertices[vIdx] = add_z_vertex(i, j, k, val0, val3);
					z_verts[currentIdx] = vIdx;
				};
			}
		}
	}
};

// get data from scalar field
float GeneratorLewiner::get_data(const int i, const int j, const int k) const {

	return scalarField[i + j * blockDims[0] + k * blockDims[0] * blockDims[1]];
};

//@brief function to find the index of the vertex in the scalar field vector based on its position in the grid
int GeneratorLewiner::find_vertex_index(int x, int y, int z) {

	return x * blockDims[1] * blockDims[2] + y * blockDims[2] + z;
};

Vec3 GeneratorLewiner::get_position(int x, int y, int z) {

	return Vec3(
		bounds[0] + x * stepX,
		bounds[2] + y * stepY,
		bounds[4] + z * stepZ);
};

void GeneratorLewiner::process_cube(int i, int j, int k, const float cube[8], int lut_entry) {

	int v12 = -1;

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

	int tv[3] = {};

	int startTriIdx = triangleCount.fetch_add(n, std::memory_order_relaxed);

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

int GeneratorLewiner::get_x_vert(int i, int j, int k) {
	return x_verts[find_vertex_index(i, j, k)];
};

int GeneratorLewiner::get_y_vert(int i, int j, int k) {
	return y_verts[find_vertex_index(i, j, k)];
};

int GeneratorLewiner::get_z_vert(int i, int j, int k) {
	return z_verts[find_vertex_index(i, j, k)];
};

int GeneratorLewiner::add_c_vertex(const int i, const int j, const int k) {

	// add test_vertex_addition();

	LVertex vertex;

	// average
	float u = 0.0f;

	// set all to zero
	vertex.x = vertex.y = vertex.z = vertex.nx = vertex.ny = vertex.nz = 0;

	// estimate the average of the intersection points of the cube
	int vid = get_x_vert(i, j, k);
	if (vid != -1) {
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
	if (vid != -1) {
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
	if (vid != -1) {
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
	if (vid != -1) {
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
	if (vid != -1) {
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
	if (vid != -1) {
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
	if (vid != -1) {
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
	if (vid != -1) {
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
	if (vid != -1) {
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
	if (vid != -1) {
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
	if (vid != -1) {
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
	if (vid != -1) {
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
	int vIdx = vertexCount.fetch_add(1, std::memory_order_relaxed);
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

void GeneratorLewiner::_update_render() {

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

	glBindVertexArray(VAO); // Optional but good practice to ensure we target correct state

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
	std::vector<uint8_t> field = get_image_field(voxelSize, blockSize);

	// 1. Calculate required Grid Dimensions
	// We want voxels to be 'desiredVoxelSize' (e.g. 0.03 mm)
	int nx = static_cast<int>((blockSize[1] - blockSize[0]) / voxelSize);
	int ny = static_cast<int>((blockSize[3] - blockSize[2]) / voxelSize);
	int nz = static_cast<int>((blockSize[5] - blockSize[4]) / voxelSize);

	std::cout << "Starting High-Res Export..." << std::endl;
	std::cout << "Physical Box: " << (blockSize[1] - blockSize[0]) << " mm" << std::endl;
	std::cout << "Target Voxel: " << voxelSize << " mm" << std::endl;
	std::cout << "Grid Size: " << nx << " x " << ny << " x " << nz << std::endl;

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

	std::cout << "Export Complete: " << fileName << std::endl;
};

void GeneratorLewiner::export_mhd(std::filesystem::path& path, float voxelSize, std::array<float, 6> blockBounds) {

	std::filesystem::path basePath = std::filesystem::absolute(path);
	basePath.replace_extension(""); // Strip extension to be safe

	std::string rawFileName = basePath.string() + ".raw";
	std::string mhdFileName = basePath.string() + ".mhd";

	// We need just the filename for the header (no full path)
	std::string rawBaseName = std::filesystem::path(rawFileName).filename().string();

	// 2. Get Data (Reuse the helper!)
	std::vector<uint8_t> field = get_image_field(voxelSize, blockBounds);

	float sizeX = blockBounds[1] - blockBounds[0];
	float sizeY = blockBounds[3] - blockBounds[2];
	float sizeZ = blockBounds[5] - blockBounds[4];

	int nx = (int)std::ceil(sizeX / voxelSize);
	int ny = (int)std::ceil(sizeY / voxelSize);
	int nz = (int)std::ceil(sizeZ / voxelSize);

	std::ofstream rawFile(rawFileName, std::ios::binary);
	if (!rawFile) {
		std::cerr << "Error: Could not open " << rawFileName << std::endl;
		return;
	}
	rawFile.write(reinterpret_cast<char*>(field.data()), field.size());
	rawFile.close();

	std::ofstream mhdFile(mhdFileName);
	if (!mhdFile) {
		std::cerr << "Error: Could not open " << mhdFileName << std::endl;
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

	std::cout << "MHD Export Complete: " << mhdFileName << std::endl;
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

		// 1. Fetch the actual vertices using the integer indices
		const LVertex& p1 = meshVertices[tri.v1];
		const LVertex& p2 = meshVertices[tri.v2];
		const LVertex& p3 = meshVertices[tri.v3];

		// 2. Calculate the Face Normal (Cross Product of Edge 1 and Edge 2)
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
	std::cout << "Successfully exported STL to: " << fileName << std::endl;
};

void GeneratorLewiner::validate_topology() {
	std::cout << "\n--- Running Topology Validation ---" << std::endl;

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

	std::cout << "Vertices: " << meshVertices.size() << std::endl;
	std::cout << "Faces: " << meshTriangles.size() << std::endl;
	std::cout << "Total Unique Edges: " << edgeCounts.size() << std::endl;
	std::cout << "Boundary Edges (Holes): " << boundaryEdges << std::endl;
	std::cout << "Non-Manifold Edges: " << nonManifoldEdges << std::endl;

	if (boundaryEdges == 0 && nonManifoldEdges == 0) {
		std::cout << "Status: SUCCESS - Mesh is 100% Watertight and 2-Manifold!" << std::endl;
	}
	else {
		std::cout << "Status: FAILED - Mesh has topological errors." << std::endl;
	}
	std::cout << "-----------------------------------\n" << std::endl;
}

void GeneratorLewiner::render_properties() {

	ImGui::SeparatorText("Parameters");
	//ImGui::InputText("Name", buffer, sizeof(buffer));
	ImGui::InputFloat("Thickness", &isoLevel, 0.01f, 1.0f, "%.3f");
	ImGui::SliderFloat("Openess", &threshold, 0.0f, 1.0f, "%.3f");

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

	// get the domain volume to update the porosity
	float domainVolume = container.get_volume();
	porosity = (1 - volume / domainVolume);

	surfaceToVolume = surfaceArea / volume;
};

void GeneratorLewiner::render_metrics() {

	ImGui::Text(("Volume: " + std::to_string(volume) + ("mm^3")).c_str(), "%.4f");
	ImGui::Text(("Total Surface: " + std::to_string(surfaceArea) + (" mm^2")).c_str(), "&.4f");
	ImGui::Text(("Surface to Volume Ratio: " + std::to_string(surfaceToVolume) + (" 1/mm")).c_str(), "&.4f");
	ImGui::Text(("Porosity: " + std::to_string(porosity * 100.0f) + " %").c_str(), "%.4f");

	// convert local thickness and std to string
	std::string ts = std::to_string(localThickness) + " std: " + std::to_string(localThicknessStd);
	std::string ps = std::to_string(localSeparation) + " std: " + std::to_string(localSeparationStd);
	ImGui::Text(("Local Thickness: " + ts + " mm").c_str(), "%.4f");
	ImGui::Text(("Pore Separation: " + ps + " mm").c_str(), "%.4f");

	ImGui::Text(("Trabecular Number: " + std::to_string(trabecularNr) + (" 1/mm")).c_str(), "&.4f");

	ImGui::Text(("Connectivity Density: " + std::to_string(connectivityDensity) + (" 1/mm^3")).c_str(), "&.4f");
	
	ImGui::Text(("Tortuosity: " + std::to_string(tortuosity)).c_str(), "%.4f");

	ImGui::Text(("Degree of Anisotropy: " + std::to_string(anisotropyDegree)).c_str(), "%.4f");

};

void GeneratorLewiner::set_resolution(const std::array<int, 3>& newResolution) {
	blockDims = newResolution;
};

void GeneratorLewiner::set_bounds(const std::array<float, 6>& newBounds) {
	bounds = newBounds;
};

void GeneratorLewiner::set_seeds(const std::vector<Vec3>& newSeeds) {
	seeds = newSeeds;
};

std::array<float, 6> GeneratorLewiner::get_bounds() const {
	return bounds;
};

Aabb GeneratorLewiner::get_aabb() const { return aabb; };

void GeneratorLewiner::estimate_local_thickness(
	float voxelSize, std::array<float, 6>& blockBounds, bool separation) {

	std::vector<uint8_t> field = get_image_field(voxelSize, blockBounds, separation);

	int nx = (int)std::ceil((blockBounds[1] - blockBounds[0]) / voxelSize);
	int ny = (int)std::ceil((blockBounds[3] - blockBounds[2]) / voxelSize);
	int nz = (int)std::ceil((blockBounds[5] - blockBounds[4]) / voxelSize);

	//int totalVoxels = nx * ny * nz;
	int totalVoxels = (int)field.size();

	// initialize all values to infinity
	std::vector<float> squaredDistanceField;
	if (separation) {
		float maxLinearDist = (float)(nx + ny + nz);
		squaredDistanceField.resize(totalVoxels, maxLinearDist);
	}
	else {
		squaredDistanceField.resize(totalVoxels, std::numeric_limits<float>::max());
	}

	// initialize to zero the empty voxels
	#pragma omp parallel for
	for (int i{ 0 }; i < totalVoxels; i++) {
		// for measuring thickness set to zero the empty space voxels
		//if (scalarField[i] >= isoLevel) {
		//	squaredDistanceField[i] = 0.0f;
		//}
		if (field[i] == 0) squaredDistanceField[i] = 0.0f;
	}

	auto find_idx = [&](int x, int y, int z) { return z * nx * ny + y * nx + x; };

	// since we have 3d data we need three phases to estimate the distance field one along each dimension, we can use a simple 1D distance transform along each dimension
	// 1st pass: forward scan along x
	#pragma omp parallel for
	for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
			// forward scan
			for (int x = 1; x < nx; x++) {
				int idx = find_idx(x, y, z);
				squaredDistanceField[idx] = std::min(
					squaredDistanceField[idx], squaredDistanceField[find_idx(x - 1, y, z)] + 1.0f);
			}
			// backward scan
			for (int x = nx - 2; x >= 0; x--) {
				int idx = find_idx(x, y, z);
				squaredDistanceField[idx] = std::min(
					squaredDistanceField[idx], squaredDistanceField[find_idx(x + 1, y, z)] + 1.0f);
			}

			// use the squared distance for next pass 
			for (int x = 0; x < nx; x++) {
				//int idx = find_vertex_index(x, y, z);
				int idx = find_idx(x, y, z);
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
				g[y] = squaredDistanceField[find_idx(x, y, z)];
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
				squaredDistanceField[find_idx(x, u, z)] = g[s[q]] + dy * dy;
			}
		}
	}

	// finally pass along z it is the same as the y pass but we need to iterate along z and keep x, y fixed
	#pragma omp parallel for
	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++) {

			std::vector<float> g(nz);

			for (int z = 0; z < nz; z++) g[z] = squaredDistanceField[find_idx(x, y, z)];
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
				squaredDistanceField[find_idx(x, y, u)] = g[s[q]] + dz * dz;
			}
		}
	}
	// now that we have our distance field we can find redundant voxels, these are not local maxima in their neighborhood, we can remove them by setting their distance to zero, this will give us a skeleton of the solid part that represents the local thickness
	std::vector<int> redundantVoxels;

	#pragma omp parallel
	{
		std::vector<int> localRedundant;
	#pragma omp for
		for (int z = 1; z < nz - 1; z++) {
			for (int y = 1; y < ny - 1; y++) {
				for (int x = 1; x < nx - 1; x++) {
					int idx = find_idx(x, y, z);

					float r = std::sqrt(squaredDistanceField[idx]);
					//if (scalarField[idx] >= isoLevel) continue; // only consider solid voxels
					if (field[idx] == 0) continue; // only consider solid voxels

					bool isLocalMax = true;
					// check 26 neighbors we already have the indices in the const std::array<Neighbor, 26> neighbors> , we can just iterate through them and check if any of the neighbors has a higher distance value than the current voxel, if it does, then the current voxel is not a local maximum and we can mark it as redundant
					for (const auto& nb : neighbors) {
						int nidx = find_idx(x + nb.dx, y + nb.dy, z + nb.dz);
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
	for (int i = 0; i < redundantVoxels.size(); i++) {

		int ridgeIdx = redundantVoxels[i];

		float radius = std::sqrt(squaredDistanceField[ridgeIdx]);
		float diameter = 2.0f * radius;
		float rSq = radius * radius;

		// get the position of the voxel in the grid
		int rz = ridgeIdx % nz;
		int ry = (ridgeIdx / nz) % ny;
		int rx = ridgeIdx / (ny * nz);

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

				// --- CRITICAL SECTION: Lock only this Z-plane ---
				// Threads processing different Z-planes can work simultaneously!

				omp_set_lock(&z_locks[z]);

				for (int x = x_start; x <= x_end; x++) {
					int targetIdx = find_idx(x, y, z);

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
	int solidVoxelCount = 0;
	float T_min = 1.5f; // Minimum thickness threshold to avoid surface noise
	// Physical voxel size

	//float voxelSize = (bounds[5] - bounds[4]) / (float)blockDim[2];

	#pragma omp parallel for reduction(+:totalThicknessSum, solidVoxelCount)
	for (int i = 0; i < totalVoxels; i++) {

		//if (scalarField[i] <= isoLevel) {
		float t = thicknessMap[i];
		// Filter out surface noise (Equation in paper refers to filtering small spheres)
		if (t >= T_min) {
			totalThicknessSum += t;
			solidVoxelCount++;
		}
	}

	// Model-Independent Mean Thickness (Tb.Th for bone)
	float meanThicknessVoxels = (solidVoxelCount > 0) ? (totalThicknessSum / solidVoxelCount) : 0.0f;

	// estimate standard deviation
	float deviationSum = 0.0f;

	#pragma omp parallel for reduction(+:deviationSum)
	for (int i = 0; i < totalVoxels; i++) {

		// Only consider voxels that are part of the SOLID structure
		float t = thicknessMap[i];
		if (t >= T_min) {
			float diff = t - meanThicknessVoxels;
			deviationSum += (diff * diff);
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

};

//@Function to get a subregion of the created mesh to compute the image metrics, we also should add
// an origin (e.g. the centroid).
std::vector<uint8_t> GeneratorLewiner::get_image_field(
	float voxelSize, std::array<float, 6>& blockBounds, bool inverse) {

	float sizeX = blockBounds[1] - blockBounds[0];
	float sizeY = blockBounds[3] - blockBounds[2];
	float sizeZ = blockBounds[5] - blockBounds[4];

	// get number of voxels
	int nx = (int)std::ceil(sizeX / voxelSize);
	int ny = (int)std::ceil(sizeY / voxelSize);
	int nz = (int)std::ceil(sizeZ / voxelSize);

	// step is the voxel size
	float step = voxelSize;

	int totalVoxels = nx * ny * nz;

	std::vector<uint8_t> field(totalVoxels);

	// Helper to clamp indices to avoid segfaults
	auto clamp_idx = [](int val, int maxVal) {
		if (val < 0) return 0;
		if (val >= maxVal) return maxVal - 1;
		return val;
	};

	#pragma omp parallel for
	for (int idx{ 0 }; idx < totalVoxels; idx++) {

		// find the indices in the new grid
		int z = idx / (nx * ny);
		int y = (idx % (nx * ny)) / nx;
		int x = idx % nx;

		// find physical position
		float pX = blockBounds[0] + x * voxelSize;
		float pY = blockBounds[2] + y * voxelSize;
		float pZ = blockBounds[4] + z * voxelSize;

		// find index of voxel in the scalar field and then interpolate
		float oldX = (pX - bounds[0]) / stepX;
		float oldY = (pY - bounds[2]) / stepY;
		float oldZ = (pZ - bounds[4]) / stepZ;

		// this is like 5.1, 7.8, 9.1, we need to find the interpolated value of the scalar field to assing 0 or 255
		// get the 8 corners
		int x0 = (int)std::floor(oldX);
		int y0 = (int)std::floor(oldY);
		int z0 = (int)std::floor(oldZ);

		x0 = clamp_idx(x0, blockDims[0]);
		y0 = clamp_idx(y0, blockDims[1]);
		z0 = clamp_idx(z0, blockDims[2]);

		int x1 = clamp_idx(x0 + 1, blockDims[0]);
		int y1 = clamp_idx(y0 + 1, blockDims[1]);
		int z1 = clamp_idx(z0 + 1, blockDims[2]);

		// this measures how far we are from x0,y0,z0
		float tx = oldX - x0;
		float ty = oldY - y0;
		float tz = oldZ - z0;

		// clamp weights
		if (tx < 0) tx = 0; if (tx > 1) tx = 1;
		if (ty < 0) ty = 0; if (ty > 1) ty = 1;
		if (tz < 0) tz = 0; if (tz > 1) tz = 1;

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

		// H. Threshold the Interpolated Value
		if (inverse) {
			if (val >= isoLevel) {
				field[idx] = 255;
			}
			else {
				field[idx] = 0;
			}
		}
		else {
			if (val < isoLevel) {
				field[idx] = 255;
			}
			else {
				field[idx] = 0;
			}
		}


	}

	return field;
};

//@brief function to estimate tortuosity of the porous structure, using the A* algorithm on the grid, we can estimate the shortest path between two points in the porous structure, and compare it to the straight line distance between those points to get an estimate of the tortuosity
bool GeneratorLewiner::estimate_tortuosity() {

	int nx = blockDims[0];
	int ny = blockDims[1];
	int nz = blockDims[2];

	int totalVoxels = scalarField.size();
	tortuosityPathEdges.clear();

	std::vector<int> parentMap(totalVoxels, -1);

	// vectors to store the gScore and fScore for each voxel, initialized to infinity
	// gscore is the cost from the start voxel to the current voxel, 
	// fscore is the estimated total cost from the start voxel to the goal voxel through the current voxel,
	// which is gscore + heuristic cost to goal
	// f(x,y,z) = g(x,y,z) + h(x,y,z). h is the heuristic function
	std::vector<float> gScore(totalVoxels, std::numeric_limits<float>::max());

	// define inlet as z = 0 plane, outlet as z = max plane
	float height = bounds[5] - bounds[4];
	float voxelSize = height / (float)nz;

	// we will use a priority queue to store the open set of voxels to explore, ordered by their fScore
	std::priority_queue<AStarNode, std::vector<AStarNode>, std::greater<AStarNode>> openSet;

	// initialize the open set with the inlet voxels (z = 0 plane), so here we add all voxels 
	// in the z=0 plane that are below the isoLevel (solid voxels) as starting points for the A* search,
	// since we want to find paths through the porous structure
	for (int x{ 10 }; x < nx - 10; x++) {
		for (int y{ 10 }; y < ny - 10; y++) {
			int idx = find_vertex_index(x, y, 0);

			// we only consider solid voxels as starting points for the A* search, 
			// since we want to find paths through the porous structure, 
			// if the voxel is above the isoLevel, it is considered solid

			// if it is empty we can start from it, since we want to find paths through the porous structure
			if (scalarField[idx] > isoLevel) {
				gScore[idx] = 0.0f;
				// straight line distance from inlet to outlet, since we are starting at z=0 and want to reach z=max, the heuristic is just the height of the box
				float heuristic = height;
				openSet.push({ idx, heuristic });

				parentMap[idx] = -1; // Roots have no parent
			}
		}
	}

	// now we can perform the A* search to find the shortest path from the inlet to the outlet, 
	// we will keep track of the best path length found to reach the outlet

	float minPathLength = std::numeric_limits<float>::infinity();
	float pathFound = false;
	int goalIndex = -1;

	while (openSet.size() > 0) {
		// get the current voxel with the lowest fScore from the open set
		AStarNode current = openSet.top();
		openSet.pop();

		// get its index
		int idx = current.idx;

		// first check if we have already found a shorter path to this voxel, if so we can skip it
		if (current.fScore > gScore[idx] + height) {
			continue;
		}

		// convert 1d to 3d indices
		int z = idx % nz;
		int y = (idx / nz) % ny;
		int x = idx / (ny * nz);

		// check if we have reached the outlet (z = max plane), if so we can update the minimum path length found
		if (z == nz - 1) {
			if (gScore[idx] < minPathLength) {
				minPathLength = (float)gScore[idx];
				pathFound = true;
				goalIndex = idx;
				break;
			}
		}

		// if we haven't reached the outlet, we can explore the neighbors of the current voxel
		for (const auto& nb : neighbors) {
			int nx_ = x + nb.dx;
			int ny_ = y + nb.dy;
			int nz_ = z + nb.dz;

			// check if the neighbor is within the grid bounds
			if (nx_ >= 1 && nx_ < nx - 1 &&
				ny_ >= 1 && ny_ < ny - 1 &&
				nz_ >= 0 && nz_ < nz) {
				int nbIdx = find_vertex_index(nx_, ny_, nz_);
				// we only consider empty voxels as valid neighbors to explore, 
				// since we want to find paths through the porous structure
				if (scalarField[nbIdx] > isoLevel) {
					// the cost to move from the current voxel to the neighbor is just the voxel size, since we are moving through a regular grid
					double tentative_gScore = gScore[idx] + (nb.cost * voxelSize);
					if (tentative_gScore < gScore[nbIdx]) {
						gScore[nbIdx] = tentative_gScore;
						// heuristic is the straight line distance from the neighbor to the outlet, which is just the remaining height in the z direction
						float heuristic = (nz - 1 - nz_) * voxelSize;
						openSet.push({ nbIdx, (float)gScore[nbIdx] + heuristic });

						parentMap[nbIdx] = idx; // Update parent map for path reconstruction
						//std::cout << pos << std::endl;
					}
				}
			}
		}
	}

	if (goalIndex == -1) {
		std::cerr << "No path found from inlet to outlet!" << std::endl;
		return false;
	}

	// update the model
	int currIdx = goalIndex;
	int vertexCount = 0;

	while (currIdx != -1) {
		int z = currIdx % nz;
		int y = (currIdx / nz) % ny;
		int x = currIdx / (ny * nz);

		Vec3 pos = get_position(x, y, z);

		tortuosityPathVertices.push_back(pos.x);
		tortuosityPathVertices.push_back(pos.y);
		tortuosityPathVertices.push_back(pos.z);

		if (vertexCount > 0) { // Add edge from previous vertex to current vertex
			tortuosityPathEdges.push_back(vertexCount - 1);
			tortuosityPathEdges.push_back(vertexCount);
		}
		currIdx = parentMap[currIdx];
		vertexCount++;
	}

	tortuosityPathModel = std::make_unique<PoreNetwork>(tortuosityPathVertices, tortuosityPathEdges);

	tortuosity = minPathLength / height; // tortuosity is the ratio of the actual path length to the straight line distance (height)

	return true;
};

void GeneratorLewiner::draw_tortuosity_path() {

	if (tortuosityPathModel) {
		tortuosityPathModel->draw();
	}
};

void GeneratorLewiner::estimate_anisotropy(int daDirectionNr, int daMinsteps, int daMaxsteps, float vcLimit) {
	
	// Global tally for MIL
	std::vector<float> totalHits(daDirectionNr, 0.0f);
	std::vector<float> totalLengths(daDirectionNr, 0.0f);

	// Create Random Directions (Fibonacci Sphere)
	std::vector<Vec3> dirs(daDirectionNr);
	const float PI = 3.14159265359f;
	const float goldenRatio = (1.0f + std::sqrt(5.0f)) * 0.5f;

	for (int i = 0; i < daDirectionNr; i++) {
		float theta = 2.0f * PI * i / goldenRatio;
		float phi = std::acos(1.0f - 2.0f * (i + 0.5f) / daDirectionNr);
		dirs[i] = Vec3(std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta), std::cos(phi));
	}

	// create a generator to sample points inside the domain
	std::mt19937 rng(1337);
	std::uniform_real_distribution<float> distX(bounds[0], bounds[1]);
	std::uniform_real_distribution<float> distY(bounds[2], bounds[3]);
	std::uniform_real_distribution<float> distZ(bounds[4], bounds[5]);
	std::uniform_real_distribution<float> distRand(0.0f, 1.0f);

	
	float maxRayLength = std::min({(bounds[1] - bounds[0]), (bounds[3] - bounds[2]) , (bounds[5] - bounds[4]) }) * 0.45;
	float rayStepSize = std::min({ stepX, stepY, stepZ }) * 0.5f;
	int maxRaySteps = static_cast<int>(maxRayLength / rayStepSize);

	float vf = 999.0f;
	int iteration = 0;

	// a vector to keep da history
	std::vector<float> daHistory;

	while ((iteration < daMinsteps || vf > vcLimit) && iteration < daMaxsteps) {

		// pick a random point
		Vec3 center(distX(rng), distY(rng), distZ(rng));

		// ensure that the point is inside
		int sI = std::clamp((int)((center.x - bounds[0]) / stepX), 0, blockDims[0] - 1);
		int sJ = std::clamp((int)((center.y - bounds[2]) / stepY), 0, blockDims[1] - 1);
		int sK = std::clamp((int)((center.z - bounds[4]) / stepZ), 0, blockDims[2] - 1);

		if (get_data(sI, sJ, sK) >= isoLevel) continue; // skip

		// increase the iteration since we continue
		iteration++;

		// parallelize the for loop to create the vectors
		#pragma omp parallel for
		for (int i = 0; i < daDirectionNr; i++) {
			int hits = 0;
			bool currentlyIn = true; // already checked that

			// perform the maximum ray steps
			for (int s = 1; s <= maxRaySteps; s++) {

				Vec3 pt = center + (dirs[i] * (rayStepSize * s));

				// get index if voxel and check its value
				int gX = std::clamp((int)((pt.x - bounds[0]) / stepX), 0, blockDims[0] - 1);
				int gY = std::clamp((int)((pt.y - bounds[2]) / stepY), 0, blockDims[1] - 1);
				int gZ = std::clamp((int)((pt.z - bounds[4]) / stepZ), 0, blockDims[2] - 1);
				bool isInside = (get_data(gX, gY, gZ) < isoLevel);

				// if we were inside and we are still inside we dont have a hit
				if (isInside != currentlyIn) {
					hits++;
					currentlyIn = isInside;
				}
			}
			// update global hits
			totalHits[i] += hits;
			totalLengths[i] += maxRayLength;
		}

		// estimate the point cloud and the covariance matrix
		Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
		for (int v = 0; v < daDirectionNr; v++) {
			// estimate mean interception length so far
			double mil = (totalHits[v] > 0) ? (totalLengths[v] / (double)totalHits[v]) : (double)totalLengths[v];
			Vec3 p = dirs[v] * (float)mil;
			cov(0, 0) += p.x * p.x; cov(0, 1) += p.x * p.y; cov(0, 2) += p.x * p.z;
			cov(1, 0) += p.y * p.x; cov(1, 1) += p.y * p.y; cov(1, 2) += p.y * p.z;
			cov(2, 0) += p.z * p.x; cov(2, 1) += p.z * p.y; cov(2, 2) += p.z * p.z;
		}
		cov /= daDirectionNr;

		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);

		// get the eigenvalues, the first is the smallest
		Eigen::Vector3d evals = solver.eigenvalues();
		float Lmin = std::sqrt(std::abs(evals(0)));
		float Lmax = std::sqrt(std::abs(evals(2)));

		// estimate degree of anisotropy
		float currentDa = (Lmax > 1e-6) ? (1.0f - (Lmin / Lmax)) : 0.0f;

		// push back to the da estimated so far
		daHistory.push_back(currentDa);

		// Step 5: Update Coefficient of Variation
		if (daHistory.size() >= 5) {
			float sum = 0, sq_sum = 0;
			for (float val : daHistory) { sum += val; sq_sum += val * val; }
			float mean = sum / daHistory.size();
			float variance = (sq_sum / daHistory.size()) - (mean * mean);
			vf = (mean > 0) ? (std::sqrt(std::abs(variance)) / mean) : 999.0f;
		}

		std::cout << "Points Sampled: " << iteration << " | DA: " << currentDa << " | CV: " << vf << "\r" << std::flush;
	};
	anisotropyDegree = daHistory.back();
}

void GeneratorLewiner::estimate_trabecular_number() {

	// a vector holding the values of mil values
	// to iteratively estimate anisotropy until the coefficient of variation falls down a limit
	std::vector<float> milVector;

	// vector counts
	int vecNr = 2000;

	// we should keep vectors for directions, hits and vector lengths
	std::vector<Vec3> dirs(vecNr);
	std::vector<float> totalHits(vecNr, 0.0f);
	std::vector<float> totalLengths(vecNr, 0.0f);

	// populate directions
	const float PI = 3.14159265359f;
	const float goldenRatio = (1.0f + std::sqrt(5.0f)) * 0.5f;
	for (int i = 0; i < vecNr; i++) {

		float theta = 2.0f * PI * i / goldenRatio;
		float phi = std::acos(1.0f - 2.0f * (i + 0.5f) / vecNr);

		Vec3 dir(std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta), std::cos(phi));
		dirs[i] = dir;
	}

	std::mt19937 rng(1337);
	std::uniform_real_distribution<float> distX(bounds[0], bounds[1]);
	std::uniform_real_distribution<float> distY(bounds[2], bounds[3]);
	std::uniform_real_distribution<float> distZ(bounds[4], bounds[5]);

	// this is the maximum ray length (almost half the minimum along dimensions)
	float maxRayLength = std::min({
		(bounds[1] - bounds[0]), (bounds[3] - bounds[2]), (bounds[5] - bounds[4])
		}) * 0.45f;
	// this the step size along the minimum dimension
	float rayStepSize = std::min({ stepX, stepY, stepZ }) * 0.5f;
	int maxSteps = static_cast<int>(maxRayLength / rayStepSize);
	
	std::vector<float> milHistory;
	float vf = 999.0f;
	int iteration = 0;
	int minIters = 10; // minimum iteration
	int maxIters = 1000; // maximum iteration

	while ((iteration < minIters || vf > 1e-2) && iteration < maxIters) {

		// 1. Pick a random point in the bounding box
		Vec3 center(distX(rng), distY(rng), distZ(rng));

		// Get the grid index of this random point
		int startI = std::clamp((int)((center.x - bounds[0]) / stepX), 0, blockDims[0] - 1);
		int startJ = std::clamp((int)((center.y - bounds[2]) / stepY), 0, blockDims[1] - 1);
		int startK = std::clamp((int)((center.z - bounds[4]) / stepZ), 0, blockDims[2] - 1);

		// If the random point is in AIR, skip it and try again. We only want to measure FOAM.
		if (get_data(startI, startJ, startK) >= isoLevel) {
			continue;
		}

		iteration++; // We found a valid foam point!

		// parallelize the vector casting
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < vecNr; i++) {

			int hits = 0;
			bool currentlyIn = true; // We already checked that the center is inside foam
			float lengthMarched = maxRayLength;

			// Traverse the grid for this specific ray
			for (int step = 1; step <= maxSteps; step++) {
				Vec3 pt = center + (dirs[i] * (rayStepSize * step));

				int gridX = std::clamp((int)((pt.x - bounds[0]) / stepX), 0, blockDims[0] - 1);
				int gridY = std::clamp((int)((pt.y - bounds[2]) / stepY), 0, blockDims[1] - 1);
				int gridZ = std::clamp((int)((pt.z - bounds[4]) / stepZ), 0, blockDims[2] - 1);

				bool isInside = (get_data(gridX, gridY, gridZ) < isoLevel);

				if (isInside != currentlyIn) {
					hits++;
					currentlyIn = isInside;
				}
			}

			totalHits[i] += hits;
			totalLengths[i] += lengthMarched;
		}

		// Build mean interception length
		float currentIterationMeanMil = 0.0f;
		for (int v = 0; v < vecNr; v++) {
			// Global vector length divided by global vector hits
			float mil = (totalHits[v] > 0) ? (totalLengths[v] / totalHits[v]) : totalLengths[v];
			currentIterationMeanMil += mil;
		}
		currentIterationMeanMil /= vecNr;
		milHistory.push_back(currentIterationMeanMil);

		if (milHistory.size() >= 3) {
			float mean = 0.0f;
			for (float val : milHistory) mean += val;
			mean /= milHistory.size();

			float variance = 0.0f;
			for (float val : milHistory) variance += (val - mean) * (val - mean);
			variance /= (milHistory.size() - 1); // Sample variance

			vf = (mean > 0.0f) ? (std::sqrt(variance) / mean) : 999.0f;

			// Print progress cleanly on one line using \r
			std::cout << "Sample: " << iteration << " | Current DA: " << mean << " | Variation Coeff: " << vf << "\r" << std::flush;
		}
	}

	float aggregateTotalLength = 0.0f;
	float aggregateTotalHits = 0.0f;

	for (int v = 0; v < vecNr; v++) {
		aggregateTotalLength += totalLengths[v];
		aggregateTotalHits += totalHits[v];
	}

	if (aggregateTotalHits > 0) {
		// MIL = Total Path Length / Total Intercepts
		float finalMIL = aggregateTotalLength / aggregateTotalHits;
		trabecularNr = 1.0f / finalMIL;
	}
	else {
		trabecularNr = 0.0f; // Truly no boundaries found
	}
};

void GeneratorLewiner::estimate_connectivity_density() {

	// we can also estimate the connectivity density since we have also the mesh
	long long V = meshVertices.size();
	long long F = meshTriangles.size();
	long long E = edgeSet.size();
	long long eulerCharacteristic = V - E + F;

	// estimate the genus
	float genus = 1.0f - (static_cast<float>(eulerCharacteristic) / 2.0f);

	// connectivity density is genus / domain volume
	connectivityDensity = genus / domainVolume;
};