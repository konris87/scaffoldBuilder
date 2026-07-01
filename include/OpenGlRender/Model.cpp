# define _USE_MATH_DEFINES
#include "Model.h"
#include <cmath>

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
BBox::BBox(const Vec3 size, const Vec3 center) {

	float xMin = center.x - size.x / 2.0f;
	float xMax = center.x + size.x / 2.0f;
	float yMin = center.y - size.y / 2.0f;
	float yMax = center.y + size.y / 2.0f;
	float zMin = center.z - size.z / 2.0f;
	float zMax = center.z + size.z / 2.0f;
		
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

// ----------------------------------------------------------------------------
Ellipsoid::Ellipsoid(const Eigen::Vector3f& center, const Eigen::Matrix3f& Rot, float r1, float r2, float r3) {

	vertices.clear();
	Eigen::Vector3f radii(r1, r2, r3);

	const float PI = 3.14159265359f;
	const float goldenRatio = (1.0f + std::sqrt(5.0f)) * 0.5f;

	for (int i = 0; i < 1000; i++) {
		float theta = 2.0f * PI * i / goldenRatio;
		float phi = std::acos(1.0f - 2.0f * (i + 0.5f) / 300);
		// Unit sphere coordinates
		Eigen::Vector3f p_unit(
			std::sin(phi) * std::cos(theta),
			std::sin(phi) * std::sin(theta),
			std::cos(phi)
		);

		// Scale by radii, Rotate by eigenvectors, Translate to center
		Eigen::Vector3f p_scaled = p_unit.cwiseProduct(radii);
		Eigen::Vector3f p_ellipsoid = center + Rot * p_scaled;

		// push to the vector list
		vertices.push_back(static_cast<float>(p_ellipsoid.x()));
		vertices.push_back(static_cast<float>(p_ellipsoid.y()));
		vertices.push_back(static_cast<float>(p_ellipsoid.z()));
	}

	_setup();

	// create the lines
	Vec3 c(center);
	Vec3 xDir(Rot.col(0).x(), Rot.col(0).y(), Rot.col(0).z()) ;
	Vec3 yDir(Rot.col(1).x(), Rot.col(1).y(), Rot.col(1).z());
	Vec3 zDir(Rot.col(2).x(), Rot.col(2).y(), Rot.col(2).z());

	// check the order
	Vec3 diff = xDir.cross(yDir);

	// check the first component. if they have the same sign the leave as it is
	if ((diff.x * zDir.x) < 0){
		Vec3 temp = yDir;
		yDir = zDir;
		zDir = temp;
	}

	Vec3 p1(c + xDir * r1);
	Vec3 p2(c + yDir * r2);
	Vec3 p3(c + zDir * r3);

	xAxis = std::make_unique<LineModel>(c, p1);
	yAxis = std::make_unique<LineModel>(c, p2);
	zAxis = std::make_unique<LineModel>(c, p3);

};

void Ellipsoid::_setup() {
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);

	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

	// position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
};

void Ellipsoid::draw() {
	glPointSize(2.0f);
	glBindVertexArray(VAO);
	glDrawArrays(GL_POINTS, 0, vertices.size() / 3);
	glBindVertexArray(0);
};

void Ellipsoid::_clean() {
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	VAO = 0;
	VBO = 0;
};

// =============================================================
// Sphere Class Implementation
// =============================================================
Sphere::Sphere(float radius, int sectorCount, int stackCount)
    : radius(radius), sectorCount(sectorCount), stackCount(stackCount) 
{
    create();
    setup();
}

void Sphere::create() {
    float x, y, z, xy;                              // vertex position
    float nx, ny, nz, lengthInv = 1.0f / radius;    // vertex normal

    float sectorStep = 2 * M_PI / sectorCount;
    float stackStep = M_PI / stackCount;
    float sectorAngle, stackAngle;

    // Generate Vertices, Normals, and TexCoords
    for(int i = 0; i <= stackCount; ++i) {
        stackAngle = M_PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
        xy = radius * cosf(stackAngle);             // r * cos(u)
        z = radius * sinf(stackAngle);              // r * sin(u)

        // add (sectorCount+1) vertices per stack
        // the first and last vertices have same position and normal, but different tex coords
        for(int j = 0; j <= sectorCount; ++j) {
            sectorAngle = j * sectorStep;           // starting from 0 to 2pi

            // vertex position (x, y, z)
            x = xy * cosf(sectorAngle);             // r * cos(u) * cos(v)
            y = xy * sinf(sectorAngle);             // r * cos(u) * sin(v)
            
            // Note: Assuming Z-up or Y-up here. This is Z-up. 
            // If your engine is Y-up, swap the y and z assignments!
            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);

            // normalized vertex normal (nx, ny, nz)
            nx = x * lengthInv;
            ny = y * lengthInv;
            nz = z * lengthInv;
            normals.push_back(nx);
            normals.push_back(ny);
            normals.push_back(nz);
        }
    }

    // Generate Indices for Triangle Strips/Lists
    int k1, k2;
    for(int i = 0; i < stackCount; ++i) {
        k1 = i * (sectorCount + 1);     // beginning of current stack
        k2 = k1 + sectorCount + 1;      // beginning of next stack

        for(int j = 0; j < sectorCount; ++j, ++k1, ++k2) {
            // 2 triangles per sector excluding first and last stacks
            // k1 => k2 => k1+1
            if(i != 0) {
                indices.push_back(k1);
                indices.push_back(k2);
                indices.push_back(k1 + 1);
            }

            // k1+1 => k2 => k2+1
            if(i != (stackCount - 1)) {
                indices.push_back(k1 + 1);
                indices.push_back(k2);
                indices.push_back(k2 + 1);
            }
        }
    }
}

void Sphere::update_radius(float newRadius) {
    if (this->radius == newRadius) return;
    
    this->radius = newRadius;
    
    // Clear old CPU data
    vertices.clear();
    normals.clear();
    indices.clear();
    
    // Generate new CPU data
    create(); 
    
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, normalsVBO);
    glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(float), normals.data(), GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_DYNAMIC_DRAW);

    glBindVertexArray(0);
}

void Sphere::setup(){

	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);

	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	glGenBuffers(1, &normalsVBO);
	glBindBuffer(GL_ARRAY_BUFFER, normalsVBO);
	glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(float), normals.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(1);

	glGenBuffers(1, &EBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

	glBindVertexArray(0);
}

void Sphere::draw() {
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
}

void Sphere::clean() {
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	glDeleteBuffers(1, &normalsVBO);
	glDeleteBuffers(1, &EBO);
	VAO = 0;
	VBO = 0;
	normalsVBO = 0;
	EBO = 0;
}