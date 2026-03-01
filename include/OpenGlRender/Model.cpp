#include "Model.h"

Model::Model(const std::vector<float>& vertices, const std::vector<unsigned int>& indices, const std::vector<float>& normals) : vertices(vertices), indices(indices), vertexNormals(normals) {

	setup();

};

void Model::draw() {
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
};

void Model::clean() {
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	glDeleteBuffers(1, &normalsVBO);
	glDeleteBuffers(1, &EBO);
};

void Model::setup() {
	
	// generate Vertex Array Object
	glGenVertexArrays(1, &VAO);
	// bind to vertex array
	glBindVertexArray(VAO);

	// generate Vertex Buffer Object to store vertex attributes
	// bind array buffer and send data
	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

	// vertex position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// vertex normal VBO
	glGenBuffers(1, &normalsVBO);
	glBindBuffer(GL_ARRAY_BUFFER, normalsVBO);
	glBufferData(GL_ARRAY_BUFFER, vertexNormals.size() * sizeof(float), vertexNormals.data(), GL_STATIC_DRAW);

	// vertex normal attribute
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(1);

	// generate element buffer Object
	// bind element array buffer and send data
	glGenBuffers(1, &EBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

	glBindVertexArray(0);
};

// --------------------------------------------------------------------------------------
CutPlane::CutPlane(float size) : size(size) {

	_setup();

};

CutPlane::~CutPlane() {
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	glDeleteBuffers(1, &EBO);
}

void CutPlane::draw(){
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
}

void CutPlane::update_model_matrix() {

	normal = normal.normalized();

	glm::vec3 normalVec = glm::vec3(normal.x, normal.y, normal.z);
	glm::vec3 centerVec = glm::vec3(center.x, center.y, center.z);

	glm::vec3 n = glm::normalize(normalVec);
	if (glm::length(n) < 0.01f) n = glm::vec3(0.0f, 0.0f, 1.0f);

	// 1. Start with an identity matrix
	glm::mat4 model = glm::mat4(1.0f);

	// 2. Translate to the center
	model = glm::translate(model, centerVec);

	glm::vec3 defaultDir(0.0f, 0.0f, 1.0f);

	if (glm::length(normalVec + defaultDir) < 0.001f) {
		// Rotate 180 degrees around Y axis
		model = glm::rotate(model, glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));
	}
	else {
		// Create a quaternion that rotates from defaultDir to n
		glm::quat rot = glm::quat(defaultDir, normalVec);
		// Convert quaternion to mat4 and apply it
		model = model * glm::mat4_cast(rot);
	}
	
	modelMatrix = model;

};

void CutPlane::_setup() {
	// generate Vertex Array Object
	glGenVertexArrays(1, &VAO);
	// bind to vertex array
	glBindVertexArray(VAO);

	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	glGenBuffers(1, &EBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	glBindVertexArray(0);
};

// -------------------------------------------------------------------
BoundingBox::BoundingBox() {
	_setup();
};

void BoundingBox::draw() {
	glBindVertexArray(VAO);
	glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
};

void BoundingBox::_setup() {

	// generate Vertex Array Object
	glGenVertexArrays(1, &VAO);
	// bind to vertex array
	glBindVertexArray(VAO);

	// generate Vertex Buffer Object to store vertex attributes
	// bind array buffer and send data
	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

	glGenBuffers(1, &EBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	// vertex position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	glBindVertexArray(0);

};

void BoundingBox::clean() {
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	VAO = 0;
	VBO = 0;
};


// -------------------------------------------------------------------------------------
BBox::BBox(float xMin, float xMax, float yMin, float yMax, float zMin, float zMax){
		
	// first point
	vertices.push_back(xMin);
	vertices.push_back(yMin);
	vertices.push_back(zMin);
		
	// second point
	vertices.push_back(xMax);
	vertices.push_back(yMin);
	vertices.push_back(zMin);

	// third point
	vertices.push_back(xMax);
	vertices.push_back(yMin);
	vertices.push_back(zMax);

	// fourth point
	vertices.push_back(xMin);
	vertices.push_back(yMin);
	vertices.push_back(zMax);

	// fifth point
	vertices.push_back(xMin);
	vertices.push_back(yMax);
	vertices.push_back(zMin);
	
	// sixth point
	vertices.push_back(xMax);
	vertices.push_back(yMax);
	vertices.push_back(zMin);

	// seventh point
	vertices.push_back(xMax);
	vertices.push_back(yMax);
	vertices.push_back(zMax);

	// eigth point
	vertices.push_back(xMin);
	vertices.push_back(yMax);
	vertices.push_back(zMax);

	_setup();
};

void BBox::draw() {
	glBindVertexArray(VAO);
	glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
};

void BBox::clean() {
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
};

void BBox::_setup() {
	// generate Vertex Array Object
	glGenVertexArrays(1, &VAO);
	// bind to vertex array
	glBindVertexArray(VAO);

	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

	glGenBuffers(1, &EBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
};

// -------------------------------------------------------------------------
void Arrow::draw() {
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, 198, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
}

void Arrow::clean() {
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	glDeleteBuffers(1, &EBO);
};

void Arrow::_setup() {

	// generate Vertex Array Object
	glGenVertexArrays(1, &VAO);
	// bind to vertex array
	glBindVertexArray(VAO);

	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	glGenBuffers(1, &EBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
};

// ------------------------------------------------------------------------------
PoreNetwork::PoreNetwork() {
	_setup();
};

PoreNetwork::PoreNetwork(
	const std::vector<float>& vertices, const std::vector<unsigned int>& indices)
	: vertices(vertices), indices(indices) {

	_setup();
};

PoreNetwork::~PoreNetwork() {
	clean();
};

void PoreNetwork::_setup() {

	// generate Vertex Array Object
	glGenVertexArrays(1, &VAO);
	// generate Vertex Buffer Object to store vertex attributes
	glGenBuffers(1, &VBO);
	// generate element buffer Object
	glGenBuffers(1, &EBO);

	// bind to vertex array
	glBindVertexArray(VAO);

	// bind array buffer and send data
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

	// bind elementa array buffer and send data
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

	// position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

};

void PoreNetwork::draw() {
	glBindVertexArray(VAO);
	glDrawElements(GL_LINES, indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
};

void PoreNetwork::clean() {
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	glDeleteBuffers(1, &EBO);
};

// ------------------------------------------------------------
Cylinder::Cylinder(const Vec3 base, const Vec3 direction, const double height, const double radius) {

	// create a basis u, v, w
	Vec3 w = direction.normalized();

	// temp vec
	Vec3 a;
	if (std::abs(w.z) < 0.999) {
		a = Vec3(0.0f, 0.0f, 1.0f);
	}
	else {
		a = Vec3(0.0f, 1.0f, 0.0f);
	}
	Vec3 u = a.cross(w).normalized();
	Vec3 v = w.cross(u);

	float thetaStep = 2.0f * static_cast<float>(PI / res);

	Vec3 top = base + w * height;

	std::cout << top << std::endl;

	std::cout << base << std::endl;

	// create the base circle
	for (int i{ 0 }; i < res; i++) {

		int idx1 = i;
		int idx2 = (i + 1) % res;

		float theta = i * thetaStep;
		float cosTheta = std::cos(theta);
		float sinTheta = std::sin(theta);

		Vec3 p = base + radius * (u * cosTheta + v * sinTheta);

		vertices.push_back(p.x);
		vertices.push_back(p.y);
		vertices.push_back(p.z);

		indices.push_back(idx1);
		indices.push_back(idx2);
	}

	// create the top circle
	for (int i{ 0 }; i < res; i++) {
		
		int idx1 = i;
		int idx2 = (i + 1) % res;

		indices.push_back(idx1 + res);
		indices.push_back(idx2 + res);

		float theta = i * thetaStep;
		float cosTheta = std::cos(theta);
		float sinTheta = std::sin(theta);

		Vec3 p = top + radius * (u * cosTheta + v * sinTheta);

		vertices.push_back(p.x);
		vertices.push_back(p.y);
		vertices.push_back(p.z);
	}

	// add 4 vertical lines connecting the two circles
	indices.push_back(0);
	indices.push_back(res);

	indices.push_back(8);
	indices.push_back(res + 8);

	indices.push_back(16);
	indices.push_back(res + 16);

	indices.push_back(24);
	indices.push_back(res + 24);

	_setup();
};

void Cylinder::_setup() {
	// generate Vertex Array Object
	glGenVertexArrays(1, &VAO);
	// bind to vertex array
	glBindVertexArray(VAO);

	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

	glGenBuffers(1, &EBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
};

void Cylinder::draw() {
	glBindVertexArray(VAO);
	glDrawElements(GL_LINES, indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
};

void Cylinder::clean() {
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
};

// --------------------------------------------------------------------------------------------
void LineModel::_setup() {

	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);

	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_DYNAMIC_DRAW);

	// position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
};

void LineModel::draw() {
	glBindVertexArray(VAO);
	glDrawArrays(GL_LINES, 0, 2);
	glBindVertexArray(0);
};

void LineModel::clean() {
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	VAO = 0;
	VBO = 0;
};

void LineModel::set_vertices(const glm::vec3& newpt1, const glm::vec3& newpt2) {

	pt1 = newpt1;
	pt2 = newpt2;

	vertices[0] = pt1.x;
	vertices[1] = pt1.y;
	vertices[2] = pt1.z;
	vertices[3] = pt2.x;
	vertices[4] = pt2.y;
	vertices[5] = pt2.z;

	//std::cout << " updating " << std::endl;
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * vertices.size(), vertices.data());
	glBindBuffer(GL_ARRAY_BUFFER, 0);
};

void LineModel::set_vertices(const Vec3& newpt1, const Vec3& newpt2) {

	pt1 = glm::vec3(newpt1.x, newpt1.y, newpt1.z);
	pt2 = glm::vec3(newpt2.x, newpt2.y, newpt2.z);

	vertices[0] = pt1.x;
	vertices[1] = pt1.y;
	vertices[2] = pt1.z;
	vertices[3] = pt2.x;
	vertices[4] = pt2.y;
	vertices[5] = pt2.z;

	//std::cout << " updating " << std::endl;
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * vertices.size(), vertices.data());
	glBindBuffer(GL_ARRAY_BUFFER, 0);
};