#include "Container.h"
#include <limits>

//@brief constructor, we use this to 
AbstractContainer::AbstractContainer(const std::string& fileName) {

	std::ifstream file(fileName, std::ios::binary);
	std::vector<openstl::Triangle> meshTriangles = openstl::deserializeStl(file);
	file.close();

	// Deduplicate vertices & get faces
	std::tuple<std::vector<openstl::Vec3>, std::vector<openstl::Face>> entries;

	//auto [meshVerts, meshFaces] = convertToVerticesAndFaces(meshTriangles);
	entries = convertToVerticesAndFaces(meshTriangles);

	meshVerts = std::get<0>(entries);
	meshFaces = std::get<1>(entries);

	// generate
	generate();

	// create the visualization model
	create();

	// get the pseudonormals
	estimate_pseudonormals();

	// create an sdf
	sdf = std::make_shared<MeshSDF>(this->triangles, this->pseudonormals);
};

void AbstractContainer::generate() {

	// Fill vertex buffer (flattened floats)
	vertices.clear();
	vertices.reserve(meshVerts.size() * 3);
	for (const auto& v : meshVerts) {
		vertices.push_back(v.x);
		vertices.push_back(v.y);
		vertices.push_back(v.z);
	}

	// Fill index buffer
	indices.clear();
	for (const auto& f : meshFaces) {
		indices.push_back(f[0]);
		indices.push_back(f[1]);
		indices.push_back(f[2]);

	}

	// Compute vertex normals (simple averaging of face normals)
	vertexNormals.assign(meshVerts.size() * 3, 0.0f);
	for (const auto& f : meshFaces) {
		Vec3 a = { meshVerts[f[0]].x, meshVerts[f[0]].y, meshVerts[f[0]].z };
		Vec3 b = { meshVerts[f[1]].x, meshVerts[f[1]].y, meshVerts[f[1]].z };
		Vec3 c = { meshVerts[f[2]].x, meshVerts[f[2]].y, meshVerts[f[2]].z };

		Vec3 n = (b - a).cross(c - a).normalize();

		//std::cout << f[0] << " " << f[1] << " " << f[2] << std::endl;

		for (int i = 0; i < 3; ++i) {
			vertexNormals[f[i] * 3 + 0] += n.x;
			vertexNormals[f[i] * 3 + 1] += n.y;
			vertexNormals[f[i] * 3 + 2] += n.z;
		}

		// Store triangle with Vec3
		triangles.push_back({
			Vec3(a.x, a.y, a.z),
			Vec3(b.x, b.y, b.z),
			Vec3(c.x, c.y, c.z),
			static_cast<unsigned int>(f[0]),
			static_cast<unsigned int>(f[1]),
			static_cast<unsigned int>(f[2]),
			});
	}
	// Normalize vertex normals
	for (size_t i = 0; i < meshVerts.size(); ++i) {
		Vec3 n(vertexNormals[i * 3 + 0],
			vertexNormals[i * 3 + 1],
			vertexNormals[i * 3 + 2]);
		n.normalize();
		vertexNormals[i * 3 + 0] = n.x;
		vertexNormals[i * 3 + 1] = n.y;
		vertexNormals[i * 3 + 2] = n.z;
	}

	// compute adjacency
	adjacency.clear();
	adjacency.resize(meshVerts.size());

	faceAdjacency.clear();
	faceAdjacency.resize(meshVerts.size());

	int currentIdx = 0;
	for (const auto& f : meshFaces) {

		adjacency[f[0]].push_back(f[1]);
		adjacency[f[0]].push_back(f[2]);
		adjacency[f[1]].push_back(f[0]);
		adjacency[f[1]].push_back(f[2]);
		adjacency[f[2]].push_back(f[0]);
		adjacency[f[2]].push_back(f[1]);

		// update also the face adjacency if needed (not required for current operations but useful for future extensions)
		faceAdjacency[f[0]].push_back(currentIdx);
		faceAdjacency[f[1]].push_back(currentIdx);
		faceAdjacency[f[2]].push_back(currentIdx);

		currentIdx++;
	}

	// Remove duplicates in adjacency
	for (auto& neighbors : adjacency) {
		std::sort(neighbors.begin(), neighbors.end());
		neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
	}

	// Compute bounds
	if (!meshVerts.empty()) {
		openstl::Vec3 v0 = meshVerts[0];
		bounds = { v0.x, v0.x, v0.y, v0.y,v0.z, v0.z };
		for (const auto& v : meshVerts) {
			bounds[0] = std::min(bounds[0], v.x);
			bounds[1] = std::max(bounds[1], v.x);
			bounds[2] = std::min(bounds[2], v.y);
			bounds[3] = std::max(bounds[3], v.y);
			bounds[4] = std::min(bounds[4], v.z);
			bounds[5] = std::max(bounds[5], v.z);
		}
	}

	// Fill edge indices if needed
	edgeIndices.clear();
	std::set<std::pair<unsigned int, unsigned int>> edgeSet;
	for (const auto& f : meshFaces) {
		std::pair<unsigned int, unsigned int> edges[3] = {
			{std::min(f[0], f[1]), std::max(f[0], f[1])},
			{std::min(f[1], f[2]), std::max(f[1], f[2])},
			{std::min(f[2], f[0]), std::max(f[2], f[0])},
		};
		for (auto& e : edges) {
			if (edgeSet.insert(e).second) {
				edgeIndices.push_back(e.first);
				edgeIndices.push_back(e.second);
			}
		}
	}
};

void AbstractContainer::gui_setup() {

	ImGui::ColorEdit4("Container Color", (float*)&color);

	ImGui::SeparatorText("Parameters");

	ImGui::Text("Vertices: %d ", vertices.size() / 3);

	ImGui::Text("Triangles: %d ", triangles.size() / 3);

	ImGui::Text("Bound Box Dimensions: x %.4f, y %.4f, z %.4f",
		bounds[1] - bounds[0],
		bounds[3] - bounds[2],
		bounds[5] - bounds[4]);

	ImGui::InputFloat("Scale", &scaleFactor, 0.01f, 100.0f);

	if (ImGui::Button("Update Scale")) {
		apply_scale();
	}
};

void AbstractContainer::estimate_pseudonormals() {

	size_t numVerts = vertices.size() / 3;

	pseudonormals.assign(numVerts, Vec3(0.0f, 0.0f, 0.0f));

	for (const auto& t : triangles) {

		// Edge vectors emanating from vertex A (for angle A)
		Vec3 e01 = t.v2 - t.v1;
		Vec3 e02 = t.v3 - t.v1;

		// Edge vectors emanating from vertex B (for angle B)
		Vec3 e10 = t.v1 - t.v2;
		Vec3 e12 = t.v3 - t.v2;

		// Edge vectors emanating from vertex C (for angle C)
		Vec3 e20 = t.v1 - t.v3;
		Vec3 e21 = t.v2 - t.v3;

		// Compute the face normal using the cross product
		Vec3 faceNormal = e01.cross(e02);
		float normalLen = faceNormal.norm();

		// Safety check: Skip mathematically degenerate (zero-area) triangles
		if (normalLen < 1e-6f) continue;

		// Normalize the face normal
		faceNormal = faceNormal / normalLen;

		// Helper lambda to safely compute the angle between two edge vectors
		auto compute_angle = [](const Vec3& vA, const Vec3& vB) {
			float lenA = vA.norm();
			float lenB = vB.norm();
			if (lenA < 1e-6f || lenB < 1e-6f) return 0.0f;

			float dotVal = vA.dot(vB) / (lenA * lenB);
			// Crucial: Clamp between -1.0 and 1.0 to prevent std::acos from returning NaN 
			// due to minor floating-point inaccuracies (e.g., dotVal = 1.0000001)
			dotVal = std::clamp(dotVal, -1.0f, 1.0f);
			return std::acos(dotVal);
		};

		// Calculate the interior angles in radians
		float angle0 = compute_angle(e01, e02);
		float angle1 = compute_angle(e10, e12);
		float angle2 = compute_angle(e20, e21);

		// Accumulate the weighted normals using the triangle's original vertex indices
		pseudonormals[t.i1] += faceNormal * angle0;
		pseudonormals[t.i2] += faceNormal * angle1;
		pseudonormals[t.i3] += faceNormal * angle2;
	}

	// 4. Normalize all the accumulated pseudonormals
	for (auto& pn : pseudonormals) {
		float len = pn.norm();
		if (len > 1e-6f) {
			pn = pn / len;
		}
		else {
			// Fallback for completely isolated/unconnected vertices in messy STLs
			pn = Vec3(0.0f, 1.0f, 0.0f);
		}
	}
};

void AbstractContainer::apply_scale() {

	for (auto& v : meshVerts) {
		v.x *= scaleFactor;
		v.y *= scaleFactor;
		v.z *= scaleFactor;
	}

	// generate
	generate();

	// create the visualization model
	create();

	// get the pseudonormals
	estimate_pseudonormals();

	// create an sdf
	sdf = std::make_shared<MeshSDF>(this->triangles, this->pseudonormals);

	updated = true;
};
