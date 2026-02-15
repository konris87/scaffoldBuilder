#include "Generator.h"
#include <fstream>
#include <queue>
#include <limits>
#include <iomanip>
#include <omp.h>

Generator::Generator(
	std::vector<Vec3>& seeds,
	const std::array<float, 6>& bounds,
	const std::array<int, 3>& blockDim,
	const float threshold, const float isoLevel, const int foam
) : seeds(seeds), bounds(bounds), blockDim(blockDim), threshold(threshold), isoLevel(isoLevel), foam(foam) {

	paddedBounds = bounds;
	workingBlockDim = blockDim;

	_setup_mesh();

	_setup_edges();
}

void Generator::populate_grids(const IContainer& con, bool pad) {
	if (pad) {
		// Pad bounds by 2 steps on each side
		float padX = stepX * 2.0f;
		float padY = stepY * 2.0f;
		float padZ = stepZ * 2.0f;

		paddedBounds[0] = bounds[0] - padX;
		paddedBounds[1] = bounds[1] + padX;
		paddedBounds[2] = bounds[2] - padY;
		paddedBounds[3] = bounds[3] + padY;
		paddedBounds[4] = bounds[4] - padZ;
		paddedBounds[5] = bounds[5] + padZ;

		// CRITICAL FIX: Increase grid resolution to match the physical expansion
		// We added 2 steps on left and 2 on right = 4 extra steps
		workingBlockDim[0] = blockDim[0] + 4;
		workingBlockDim[1] = blockDim[1] + 4;
		workingBlockDim[2] = blockDim[2] + 4;
	}
	else {
		paddedBounds = bounds;
		workingBlockDim = blockDim;
	}
	
	if (seeds.size() < 3) {
		return;
	}

	stepX = (paddedBounds[1] - paddedBounds[0]) / (workingBlockDim[0] - 1);
	stepY = (paddedBounds[3] - paddedBounds[2]) / (workingBlockDim[1] - 1);
	stepZ = (paddedBounds[5] - paddedBounds[4]) / (workingBlockDim[2] - 1);

	// first populate the kdtree with the seeds
	kdtree = std::make_unique<Kdtree>(seeds);

	// update scalar field
	scalarField.clear();

	// it will be more efficient to reserve the size of the scalar field vector before filling it, 
	// since we know the size of the grid in advance
	
	int totalVoxels = workingBlockDim[0] * workingBlockDim[1] * workingBlockDim[2];
	//scalarField.reserve(totalVoxels);
	scalarField.resize(totalVoxels, 9999.9);

	//#pragma omp parallel for collapse(3)
	// create the grid points based on the bounds and block dimensions
	for (int i{ 0 }; i < workingBlockDim[0]; i++) {
		for (int j{ 0 }; j < workingBlockDim[1]; j++) {
			for (int k{ 0 }; k < workingBlockDim[2]; k++) {

				int idx = find_vertex_index(i, j, k);
				float x = paddedBounds[0] + i * stepX;
				float y = paddedBounds[2] + j * stepY;
				float z = paddedBounds[4] + k * stepZ;
				Vec3 point(x, y, z);

				//float containerDist = con.sdf->compute_distance(point);
				//if (containerDist > 2.0f) { // If it's more than 2mm outside, don't bother
				//	continue; // scalarField[idx] remains 9999.9f
				//}

				// we need to find the two nearest seeds to the point, and compute the distance to the nearest seed, and the distance to the second nearest seed
				auto neighbors = kdtree->knn(point, 3, [](const Vec3 &p1, const Vec3& p2) {
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
				if(foam){
					//value = (smoothd2 - smoothd1) + threshold * (d3 - d2);
					value = d2 - d1;
				}
				else {
					value = (d3 - d1) + threshold * (smoothd2 - smoothd1);
				}
				
				//double value = (d3 - d1) + threshold * (smoothd2 - smoothd1); // this works
				//float containerDist = con.sdf->compute_distance(point);
				//float mappedContainer = containerDist + isoLevel;

				//scalarField.push_back(std::max(value, mappedContainer));
				scalarField[idx] = value;
			}
		}
	}

	std::cout << "scalar field size: " << scalarField.size() << std::endl;

	smooth_scalar_field();

	// apply clamping
	int solidVoxels = 0;
	//#pragma omp parallel for collapse(3) reduction(+:solidVoxels)
	for (int i = 0; i < workingBlockDim[0]; i++) {
		for (int j = 0; j < workingBlockDim[1]; j++) {
			for (int k = 0; k < workingBlockDim[2]; k++) {

				int idx = find_vertex_index(i, j, k);

				float x = paddedBounds[0] + i * stepX;
				float y = paddedBounds[2] + j * stepY;
				float z = paddedBounds[4] + k * stepZ;
				Vec3 point(x, y, z);

				float containerDist = con.sdf->compute_distance(point);
				float mappedContainer = containerDist + isoLevel;

				// Mathematical Boolean Intersection
				scalarField[idx] = std::max(scalarField[idx], mappedContainer);

				//// Optional Padding Check (Uncomment if needed)
				//if (pad && (i == 0 || i == workingBlockDim[0] - 1 ||
				//	j == 0 || j == workingBlockDim[1] - 1 ||
				//	k == 0 || k == workingBlockDim[2] - 1)) {
				//	scalarField[idx] = 9999.9f;
				//}

				// Count metrics now (after final modification)
				if (scalarField[idx] < isoLevel) {
					solidVoxels++;
				}
			}
		}
	}

	volumeFraction = (float)solidVoxels / (float)totalVoxels;
	porosity = 1.0 - volumeFraction;

	volume = volumeFraction * (paddedBounds[1] - paddedBounds[0]) * (paddedBounds[3] - paddedBounds[2]) * (paddedBounds[5] - paddedBounds[4]);
};


void Generator::smooth_scalar_field() {
	std::vector<float> smoothed = scalarField; // Copy

	// Simple 3x3x3 Box Blur
	for (int z = 1; z < workingBlockDim[2] - 1; z++) {
		for (int y = 1; y < workingBlockDim[1] - 1; y++) {
			for (int x = 1; x < workingBlockDim[0] - 1; x++) {

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

void Generator::marching_cubes() {

	// reset
	globalIndex = 0;
	globalVertexMap.clear();
	vertices.clear();
	indices.clear();
	edgeIndices.clear();
	normals.clear();
	tempNormals.clear();
	surfaceArea = 0.0;
	triangles.clear();

	int solidVoxels = 0;

	graph.reset(new Graph());

	// iterate over the grid cells
	for (int i{ 0 }; i < workingBlockDim[0] - 1; i++) {
		for (int j{ 0 }; j < workingBlockDim[1] - 1; j++) {
			for (int k{ 0 }; k < workingBlockDim[2] - 1; k++) {

				// estimate positions of the cube vertices
				Vec3 p[8] = {};

				p[0] = get_position(i, j, k);
				p[1] = get_position(i + 1, j, k);
				p[2] = get_position(i + 1, j, k + 1);
				p[3] = get_position(i, j, k + 1);
				p[4] = get_position(i, j + 1, k);
				p[5] = get_position(i + 1, j + 1, k);
				p[6] = get_position(i + 1, j + 1, k + 1);
				p[7] = get_position(i, j + 1, k + 1);

				// estimate values at each vertex
				float values[8];
				values[0] = scalarField[find_vertex_index(i, j, k)];
				values[1] = scalarField[find_vertex_index(i + 1, j, k)];
				values[2] = scalarField[find_vertex_index(i + 1, j, k + 1)];
				values[3] = scalarField[find_vertex_index(i, j, k + 1)];
				values[4] = scalarField[find_vertex_index(i, j + 1, k)];
				values[5] = scalarField[find_vertex_index(i + 1, j + 1, k)];
				values[6] = scalarField[find_vertex_index(i + 1, j + 1, k + 1)];
				values[7] = scalarField[find_vertex_index(i, j + 1, k + 1)];


				// determine the cube index based on the values at the vertices and the isoLevel
				int idx = find_cube_index(isoLevel, values);

				// check if zero vertices are inside the surface, if they are, we skip this cube
				if (edgeTable[idx] == 0) continue;

				// for each edge of the cube, if it is intersected by the surface, 
				// we estimate the vertex position on the edge based on the scalar values at the endpoints and the isoLevel
				// using the interpolate_vertex function
				// the cube has 12 edges, and we check each edge based on the edge table
				Vec3 vertList[12] = {};
				
				if (edgeTable[idx] & 1) vertList[0] = interpolate_vertex(isoLevel, p[0], p[1], values[0], values[1]);
				if (edgeTable[idx] & 2) vertList[1] = interpolate_vertex(isoLevel, p[1], p[2], values[1], values[2]);
				if (edgeTable[idx] & 4) vertList[2] = interpolate_vertex(isoLevel, p[2], p[3], values[2], values[3]);
				if (edgeTable[idx] & 8) vertList[3] = interpolate_vertex(isoLevel, p[3], p[0], values[3], values[0]);
				if (edgeTable[idx] & 16) vertList[4] = interpolate_vertex(isoLevel, p[4], p[5], values[4], values[5]);
				if (edgeTable[idx] & 32) vertList[5] = interpolate_vertex(isoLevel, p[5], p[6], values[5], values[6]);
				if (edgeTable[idx] & 64) vertList[6] = interpolate_vertex(isoLevel, p[6], p[7], values[6], values[7]);
				if (edgeTable[idx] & 128) vertList[7] = interpolate_vertex(isoLevel, p[7], p[4], values[7], values[4]);
				if (edgeTable[idx] & 256) vertList[8] = interpolate_vertex(isoLevel, p[0], p[4], values[0], values[4]);
				if (edgeTable[idx] & 512) vertList[9] = interpolate_vertex(isoLevel, p[1], p[5], values[1], values[5]);
				if (edgeTable[idx] & 1024) vertList[10] = interpolate_vertex(isoLevel, p[2], p[6], values[2], values[6]);
				if (edgeTable[idx] & 2048) vertList[11] = interpolate_vertex(isoLevel, p[3], p[7], values[3], values[7]);

				// now we have the vertex positions for the edges that are intersected by the surface,
				// we can create the triangles based on the triangle table

				for (int t=0; triangleTable[idx][t] != -1; t += 3) {
					
					Triangle tri;
					tri.v1 = vertList[triangleTable[idx][t]];
					tri.v2 = vertList[triangleTable[idx][t + 1]];
					tri.v3 = vertList[triangleTable[idx][t + 2]];

					// estimate face normal
					Vec3 normal = (tri.v2 - tri.v1).cross(tri.v3 - tri.v1).normalized();

					if (std::isnan(normal.x)) continue;

					tri.normal = normal;

					// push also the vertices for rendering, 
					// we should also check if the vertex is already in the list of vertices, to avoid duplicates
					triangles.push_back(tri);

					// estimate triangle area and add it to the total area
					surfaceArea += 0.5 * (tri.v2 - tri.v1).cross(tri.v3 - tri.v1).norm();

					int idx1 = check_vertex(tri.v1);
					int idx2 = check_vertex(tri.v2);
					int idx3 = check_vertex(tri.v3);

					//if (tempNormals.size() <= std::max({ idx1, idx2, idx3 })) {
					//	tempNormals.resize(std::max({ idx1, idx2, idx3 }) + 1, { 0,0,0 });
					//}

					// push to the temporary normals vector, we will average them later for each vertex
					tempNormals[idx1] = tempNormals[idx1] + normal;
					tempNormals[idx2] = tempNormals[idx2] + normal;
					tempNormals[idx3] = tempNormals[idx3] + normal;

					indices.push_back(idx1);
					indices.push_back(idx2);
					indices.push_back(idx3);

					// add also edges for rendering but we need to check the indices
					add_edge(idx1, idx2, idx3);
				}	
			}
		}
	}

	// one normal per vertex, we average the face normals for each vertex
	normals.reserve(vertices.size());

	for (int i{0}; i < tempNormals.size(); i++) {
		Vec3 n = tempNormals[i].normalized();
		normals.push_back(n.x);
		normals.push_back(n.y);
		normals.push_back(n.z);
	}

	//std::cout << "Triangles: " << triangles.size() << std::endl;

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

void Generator::render_properties() {

	ImGui::SeparatorText("Parameters");
	//ImGui::InputText("Name", buffer, sizeof(buffer));
	ImGui::InputFloat("Thickness", &isoLevel, 0.01f, 1.0f, "%.3f");
	ImGui::SliderFloat("Openess", &threshold, 0.0f, 1.0f, "%.3f");

};

int Generator::check_vertex(Vec3& vec) {

	int idx = -1;

	GlobalVertexVec3 key(vec);

	// check if it is inside
	auto it = globalVertexMap.find(key);
	if (it == globalVertexMap.end()) {

		idx = globalIndex;

		// This is a new vertex
		globalVertexMap[key] = globalIndex;
		
		globalIndex++;

		// push also the rendering vertices
		vertices.push_back(vec.x);
		vertices.push_back(vec.y);
		vertices.push_back(vec.z);

		tempNormals.emplace_back(0.0f, 0.0f, 0.0f);
	}
	else {
		idx = it->second;
	}

	return idx;
};

void Generator::add_edge(int idx1, int idx2, int idx3) {

	// Add also edges
	std::pair<unsigned int, unsigned int> edges[3] = {
		{ std::min(idx1, idx3), std::max(idx1, idx2) },
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

//@brief function to find the index of the vertex in the scalar field vector based on its position in the grid
int Generator::find_vertex_index(int x, int y, int z) {

	return x * workingBlockDim[1] * workingBlockDim[2] + y * workingBlockDim[2] + z;
};

//@brief function to determine the cube index based on the values at the vertices and the isoLevel,
// this is used to look up the edge table for marching cubes
int Generator::find_cube_index(float isoLevel, float values[8]) {
	int idx = 0;

	// loop for every vertex and check if its value is below the isoLevel, if it is, we set the corresponding bit in the index
	for (int i{ 0 }; i < 8; i++) {
		if (values[i] < isoLevel) {

			// shift the bit to the left by i positions and set it in the index 8->256
			idx |= (1 << i);
		}
	}
	return idx;
};

Vec3 Generator::interpolate_vertex(
	float isoLevel,
	const Vec3& p1,
	const Vec3& p2,
	float val1, float val2) {

	if (std::abs(val1 - val2) < 0.00001) return p1;
	if (std::abs(val1 - isoLevel) < 0.00001) return p1;
	if (std::abs(val2 - isoLevel) < 0.00001) return p2;

	// estimate the interpolation factor, this is linear interpolation between
	// the two points based on the scalar values at the endpoints and the isoLevel
	float mu = (isoLevel - val1) / (val2 - val1);

	// interpolate the vertex position
	return Vec3(
		p1.x + mu * (p2.x - p1.x),
		p1.y + mu * (p2.y - p1.y),
		p1.z + mu * (p2.z - p1.z)
	);
};

//@function to get position of vertex in the 3D grid
Vec3 Generator::get_position(int x, int y, int z) {

	//float stepX = (paddedBounds[1] - paddedBounds[0]) / (workingBlockDim[0] - 1);
	//float stepY = (paddedBounds[3] - paddedBounds[2]) / (workingBlockDim[1] - 1);
	//float stepZ = (paddedBounds[5] - paddedBounds[4]) / (workingBlockDim[2] - 1);

	return Vec3(
		paddedBounds[0] + x * stepX,
		paddedBounds[2] + y * stepY,
		paddedBounds[4] + z * stepZ);
};

void Generator::_setup_mesh() {

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

void Generator::_setup_edges() {

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

void Generator::draw() {

	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
};

void Generator::draw_edges() {
	glBindVertexArray(edgeVAO);
	glDrawElements(GL_LINES, edgeIndices.size(), GL_UNSIGNED_INT, 0);
};

void Generator::draw_tortuosity_path() {

	if (tortuosityPathModel) {
		tortuosityPathModel->draw();
	}
};

void Generator::draw_pore_network() {
	if (connectivityPathModel) {
		connectivityPathModel->draw();
	}
};

//@brief function to compute a smooth minimum between two values, 
// this is used to create a smooth transition between the distance fields of the seeds,
// and to create a smoother mesh
//@param a: the first value
//@param b: the second value
//@param k: the smoothing factor, higher values create a sharper transition,
// lower values create a smoother transition
//@returns the smooth minimum between a and b
float smin(float a, float b, float k) {
	float h = std::max(k - std::abs(a - b), 0.0f) / k;
	return std::min(a, b) - h * h * h * k * (1.0f / 6.0f);
}

void Generator::render_metrics() {

	ImGui::Text(("Volume: " + std::to_string(volume)).c_str(), "%.4f");
	ImGui::Text(("Porosity: " + std::to_string(porosity)).c_str(), "%.4f");

	// convert local thickness and std to string
	std::string ts = std::to_string(localThickness) + " std: " + std::to_string(localThicknessStd);
	std::string ps = std::to_string(localSeparation) + " std: " + std::to_string(localSeparationStd);
	ImGui::Text(("Local Thickness: " + ts).c_str(), "%.4f");
	ImGui::Text(("Pore Separation: " + ps).c_str(), "%.4f");
	ImGui::Text(("Tortuosity: " + std::to_string(tortuosity)).c_str(), "%.4f");

};

Metrics Generator::get_metrics() {

	return Metrics{
		surfaceArea,
		porosity,
		volumeFraction,
		volume,
		tortuosity,
		localThickness,
		localThicknessStd,
		localSeparation,
		localSeparationStd
	};

};

void Generator::export_stl(std::string fileName) {

	std::ofstream out(fileName, std::ios::binary);

	// header
	char header[80] = { 0 };
	out.write(header, 80);

	// number of triangles 4 bytes
	uint32_t numTriangles = static_cast<uint32_t>(triangles.size());
	out.write(reinterpret_cast<char*>(&numTriangles), sizeof(uint32_t));

	// pass triangle data
	for (const auto& tri : triangles) {

		float data[12];

		data[0] = (float)tri.normal.x; data[1] = (float)tri.normal.y; data[2] = (float)tri.normal.z;
		data[3] = (float)tri.v1.x; data[4] = (float)tri.v1.y; data[5] = (float)tri.v1.z;
		data[6] = (float)tri.v2.x; data[7] = (float)tri.v2.y; data[8] = (float)tri.v2.z;
		data[9] = (float)tri.v3.x; data[10] = (float)tri.v3.y; data[11] = (float)tri.v3.z;

		out.write(reinterpret_cast<char*>(data), 12 * sizeof(float));

		uint16_t attr = 0;
		out.write(reinterpret_cast<char*>(&attr), sizeof(uint16_t));
	}

	out.close();
};

//@brief function to estimate tortuosity of the porous structure, using the A* algorithm on the grid, we can estimate the shortest path between two points in the porous structure, and compare it to the straight line distance between those points to get an estimate of the tortuosity
bool Generator::estimate_tortuosity() {

	int nx = blockDim[0];
	int ny = blockDim[1];
	int nz = blockDim[2];

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
	for(int x{ 10 }; x < nx - 10; x++) {
		for(int y{ 10 }; y < ny - 10; y++) {
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

	while(openSet.size() > 0) {
		// get the current voxel with the lowest fScore from the open set
		AStarNode current = openSet.top();
		openSet.pop();

		// get its index
		int idx = current.idx;
		
		// first check if we have already found a shorter path to this voxel, if so we can skip it
		if(current.fScore > gScore[idx] + height) {
			continue;
		}

		// convert 1d to 3d indices
		int z = idx % nz;
		int y = (idx / nz) % ny;
		int x = idx / (ny * nz);

		// check if we have reached the outlet (z = max plane), if so we can update the minimum path length found
		if(z == nz - 1) {
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
			if (nx_ >= 1 && nx_ < nx - 1&&
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

	while(currIdx!= -1) {
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

	// edges should be half the number of vertices since each edge connects two vertices
	std::cout << "tortuosity vertices: " << tortuosityPathVertices.size() / 3 << std::endl;	
	std::cout << "tortuosity edges: " << tortuosityPathEdges.size() / 2 << std::endl;

	tortuosity = minPathLength / height; // tortuosity is the ratio of the actual path length to the straight line distance (height)

	return true;
};

bool Generator::estimate_connectivity_network() {

	connectivityPathVertices.clear();
	connectivityPathEdges.clear();	

	for (const auto& s : seeds) {
		connectivityPathVertices.push_back(s[0]);
		connectivityPathVertices.push_back(s[1]);
		connectivityPathVertices.push_back(s[2]);
	}

	// create a voro container
	voro::container con = voro::container(
		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5],
		blockDim[0], blockDim[1], blockDim[2], false, false, false, 16
	);

	for (int i{ 0 }; i < seeds.size(); i++) {
		con.put(i, seeds[i][0], seeds[i][1], seeds[i][2]);
	}

	//create a kdtree for the seed points to find nearest neighbors
	Kdtree tree(seeds);

	voro::c_loop_all cla(con);
	voro::voronoicell_neighbor cell;

	if (cla.start()) do if (con.compute_cell(cell, cla)) {
	
		int seedId = cla.pid();

		// get position of seed and store it to an array
		double px = 0.0, py = 0.0, pz = 0.0;

		// seed position
		cla.pos(px, py, pz);

		// cell vertices
		std::vector<double> cellVertices;
		cell.vertices(px, py, pz, cellVertices);

		// get the neigbor ids
		std::vector<int> neighbors;
		cell.neighbors(neighbors);

		// cell face indices
		std::vector<int> faceIndices;
		cell.face_vertices(faceIndices);

		int idxOffset = 0;

		const int numFaces = static_cast<int>(neighbors.size());

		for (int i = 0; i < numFaces; ++i) {

			int neighborId = neighbors[i];
			int order = faceIndices[idxOffset];

			if (neighborId > seedId) {
			
				double cx = 0.0, cy = 0.0, cz = 0.0;

				for (int j = 1; j <= order; ++j) {
					int vertIdx = faceIndices[idxOffset + j];

					// Look up coordinates (Stride is 3)
					cx += cellVertices[3 * vertIdx + 0];
					cy += cellVertices[3 * vertIdx + 1];
					cz += cellVertices[3 * vertIdx + 2];
				}

				// centroid of the face
				cx /= order;
				cy /= order;
				cz /= order;

				// sample scalar field, first find grid cell that contains the centroid
				int gx = static_cast<int>((cx - bounds[0]) / (bounds[1] - bounds[0]) * (blockDim[0] - 1));
				int gy = static_cast<int>((cy - bounds[2]) / (bounds[3] - bounds[2]) * (blockDim[1] - 1));
				int gz = static_cast<int>((cz - bounds[4]) / (bounds[5] - bounds[4]) * (blockDim[2] - 1));


				// check if it is within bounds
				if(gx >=0 && gx < blockDim[0] &&
					gy >= 0 && gy < blockDim[1] &&
					gz >= 0 && gz < blockDim[2]) {
					int idx = find_vertex_index(gx, gy, gz);
					double value = scalarField[idx];
					// if the value is above the isoLevel, we consider it a valid connection between the two seeds
					if (value > isoLevel) {
						connectivityPathEdges.push_back(seedId);
						connectivityPathEdges.push_back(neighborId);
					}
				}
			}			

			// update offset
			idxOffset += 1 + order; // Move past the count and the vertices of this face

		}
	
	} while (cla.inc());

	// prepare the network
	connectivityPathModel = std::make_unique<PoreNetwork>(connectivityPathVertices, connectivityPathEdges);

	return true;
};

void Generator::set_bounds(std::array<float, 6>& newBounds) {
	bounds = newBounds;
};

void Generator::set_thickness(const float newThickness) { isoLevel = newThickness; };

void Generator::set_openess(const float newOpeness) { threshold = newOpeness; };

void Generator::set_seeds(const std::vector<Vec3>& newSeeds) {
	seeds.clear();
	seeds = newSeeds;
};

void Generator::set_resolution(const std::array<int, 3>& newResolution) { blockDim = newResolution; };

void Generator::export_mhd(std::filesystem::path& path, float voxelSize, std::array<float, 6> blockBounds) {

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

void Generator::export_nrrd(const std::string fileName, float voxelSize, std::array<float, 6> blockSize) {

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

void Generator::estimate_local_separation() {

	int nx = workingBlockDim[0];
	int ny = workingBlockDim[1];
	int nz = workingBlockDim[2];

	int totalVoxels = nx * ny * nz;
	float maxLinearDist = (float)(nx + ny + nz);
	// initialize all values to infinity
	std::vector<float> squaredDistanceField(totalVoxels, maxLinearDist);

	// initialize to zero the empty voxels
	for (int i{ 0 }; i < totalVoxels; i++) {
		// for measuring thickness set to zero the empty space voxels
		if (scalarField[i] < isoLevel) {
			squaredDistanceField[i] = 0.0f;
		}
	}

	// since we have 3d data we need three phases to estimate the distance field one along each dimension, we can use a simple 1D distance transform along each dimension
	// 1st pass: forward scan along x
	#pragma omp parallel for
	for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
			// forward scan
			for (int x = 1; x < nx; x++) {
				int idx = find_vertex_index(x, y, z);
				squaredDistanceField[idx] = std::min(
					squaredDistanceField[idx], squaredDistanceField[find_vertex_index(x - 1, y, z)] + 1.0f);
			}
			// backward scan
			for (int x = nx - 2; x >= 0; x--) {
				int idx = find_vertex_index(x, y, z);
				squaredDistanceField[idx] = std::min(
					squaredDistanceField[idx], squaredDistanceField[find_vertex_index(x + 1, y, z)] + 1.0f);
			}

			// use the squared distance for next pass 
			for (int x = 0; x < nx; x++) {
				int idx = find_vertex_index(x, y, z);
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
				g[y] = squaredDistanceField[find_vertex_index(x, y, z)];
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
				squaredDistanceField[find_vertex_index(x, u, z)] = g[s[q]] + dy * dy;
			}
		}
	}

	// finally pass along z it is the same as the y pass but we need to iterate along z and keep x, y fixed
	#pragma omp parallel for
	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++) {

			std::vector<float> g(nz);

			for (int z = 0; z < nz; z++) g[z] = squaredDistanceField[find_vertex_index(x, y, z)];
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
				squaredDistanceField[find_vertex_index(x, y, u)] = g[s[q]] + dz * dz;
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
					int idx = find_vertex_index(x, y, z);

					float r = std::sqrt(squaredDistanceField[idx]);
					if (scalarField[idx] < isoLevel) continue; // skip solid voxels, only find ridges in empty space
					
					bool isLocalMax = true;
					// check 26 neighbors we already have the indices in the const std::array<Neighbor, 26> neighbors> , we can just iterate through them and check if any of the neighbors has a higher distance value than the current voxel, if it does, then the current voxel is not a local maximum and we can mark it as redundant
					for (const auto& nb : neighbors) {
						int nidx = find_vertex_index(x + nb.dx, y + nb.dy, z + nb.dz);
						if (scalarField[nidx] < isoLevel) continue; // skip solid voxels, only find ridges in empty space

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

	std::cout << redundantVoxels.size() << std::endl;

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
					int targetIdx = find_vertex_index(x, y, z);

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

	float voxelSize = (bounds[5] - bounds[4]) / (float)blockDim[2];

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

	localSeparation = meanThicknessVoxels * voxelSize;

	localSeparationStd = stdDevVoxels * voxelSize; // Final conversion to mm

	std::cout << "done" << std::endl;
};

void Generator::estimate_local_thickness(float voxelSize, std::array<float, 6>& blockBounds, bool separation) {

	std::vector<uint8_t> field = get_image_field(voxelSize, blockBounds, separation);

	int nx = (int)std::ceil((blockBounds[1] - blockBounds[0]) / voxelSize);
	int ny = (int)std::ceil((blockBounds[3] - blockBounds[2]) / voxelSize);
	int nz = (int)std::ceil((blockBounds[5] - blockBounds[4]) / voxelSize);

	//int totalVoxels = nx * ny * nz;
	int totalVoxels = (int)field.size();

	// initialize all values to infinity
	std::vector<float> squaredDistanceField;
	if (separation){
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

	std::cout << redundantVoxels.size() << std::endl;

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
std::vector<uint8_t> Generator::get_image_field(float voxelSize, std::array<float, 6>& blockBounds, bool inverse) {

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

		x0 = clamp_idx(x0, blockDim[0]);
		y0 = clamp_idx(y0, blockDim[1]);
		z0 = clamp_idx(z0, blockDim[2]);

		int x1 = clamp_idx(x0 + 1, blockDim[0]);
		int y1 = clamp_idx(y0 + 1, blockDim[1]);
		int z1 = clamp_idx(z0 + 1, blockDim[2]);

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