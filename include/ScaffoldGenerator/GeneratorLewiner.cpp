#include "GeneratorLewiner.h"
#include "LookUpTable.h"
#include "Math/Vec.h"
#include "Math/Kdtree.h"
#include <chrono>
#include <memory>
#include <iostream>
#include <map>
#include <algorithm>

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

};


void GeneratorLewiner::compute_scalar_field(const IContainer& con) {

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

	std::cout << "mesh vertices size: " << meshVertices.size() << std::endl;
	std::cout << "mesh triangles size: " << meshTriangles.size() << std::endl;

	// update opengl objects
	_update_render();

	validate_topology();

	const auto finish{ std::chrono::steady_clock::now() };
	const std::chrono::duration<double> elapsed_seconds{ finish - start };

	std::cout << "Elapsed time (s): " << 
		std::chrono::duration_cast<std::chrono::seconds>(elapsed_seconds).count() << std::endl;

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

