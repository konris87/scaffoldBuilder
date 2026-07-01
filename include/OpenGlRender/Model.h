#ifndef MODEL_H
#define MODEL_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <vector>
#include <set>
#include <iostream>
#include <filesystem>
#include "Math/Quaternion.h"
#include "Math/Vec.h"

class Model {

public:

	unsigned int VAO{ 0 }, edgeVAO{ 0 };

	Model() { setup(); };
	~Model() { clean(); };

	Model(const std::vector<float>& vertices, const std::vector<unsigned int>& indices, const std::vector<float>& normals);

	void draw();

private:

	std::vector<float> vertices;
	std::vector<unsigned int> indices;
	std::vector<float> vertexNormals;
	unsigned int VBO{ 0 }, EBO{ 0 }, normalsVBO{ 0 };

	void clean();

	void setup();
};
//private:
//	vtkSmartPointer<vtkPolyData> polyData;
	//std::vector<float> vertices;
	//std::vector<unsigned int> indices;
//	std::vector<unsigned int> edgeIndices;
//	std::vector<float> vertexNormals;
//	unsigned int VBO{ 0 }, EBO{ 0 }, normalsVBO{ 0 };
//	unsigned int edgeVBO{ 0 }, edgeEBO{ 0 };
//
//	int vertexSize{ 0 };
//	int cellSize{ 0 };
//	int edgeSize{ 0 };
//	double volume{ 0.0 };
//	double porosity{ 0.0 };
//
//	void setup_data();
//
//	void setup_mesh(); 
//
//	void setup_edges();
//};

// ------------------------------------------------------------------------------------------

class CutPlane {

public:
	CutPlane(float size = 1000.0f);
	~CutPlane();

	void draw();

	void clean();

	void update_model_matrix();

	glm::mat4 modelMatrix = glm::mat4(1.0f);

	Vec3 center{ 0.0f, 0.0f, 0.0f };
	Vec3 normal{ 0.0f, 0.0f, 1.0f };

private:

	unsigned int VAO{ 0 }, VBO{ 0 }, EBO{ 0 };
	
	float size{ 15.0f };

	float vertices[12] = {
		-size, -size, 0.0f,
		size, -size, 0.0f,
		size, size, 0.0f,
		-size, size, 0.0f
	};

	unsigned int indices[6] = {
		0, 1, 2, 2, 3, 0
	};

	void _setup();

};

// 
class BoundingBox {

public:
	BoundingBox();
	~BoundingBox() {};

	void draw();

	void clean();

private:
	std::array<float, 6> bounds;

	std::array<float, 24> vertices{

		-0.5f, -0.5f, -0.5f,
		0.5f, -0.5f, -0.5f,
		0.5f, 0.5f, -0.5f,
		-0.5f, 0.5f, -0.5f,
		-0.5f, -0.5f, 0.5f,
		0.5f, -0.5f, 0.5f,
		0.5f, 0.5f, 0.5f,
		-0.5f, 0.5f, 0.5f,
	};

	unsigned int indices[24] = {
	0, 1, 1, 2, 2, 3, 3, 0,
	4, 5, 5, 6, 6, 7, 7, 4,
	0, 4, 3, 7, 2, 6, 1, 5
	};

	unsigned int VAO{ 0 }, VBO{ 0 }, EBO{ 0 };

	void _setup();
};

// 
class BBox {
public:
	BBox() { _setup(); };
	~BBox() { clean(); };

	BBox(const Vec3 size, const Vec3 center = {0.0f, 0.0f, 0.0f});

	void draw();

	void clean();

private:
	Vec3 center { 0.0f, 0.0f, 0.0f };
	std::vector<float> vertices;
	unsigned int indices[24] = {
		0, 1, 1, 2, 2, 3, 3, 0,
		4, 5, 5, 6, 6, 7, 7, 4,
		0, 4, 3, 7, 2, 6, 1, 5
	};
	unsigned int VAO{ 0 }, VBO{ 0 }, EBO{ 0 };

	void _setup();
};

class Arrow {
public:

	Arrow() {
		_setup();
	};
	~Arrow() {
		clean();
	};

	void draw();

	void clean();

private:
	float vertices[183] = {
	0.0f, -0.1f, 0.65f, 0.0f, -0.1f, 0.0f, 0.05f, -0.0866f, 0.65f, 
	0.05f, -0.0866f, 0.0f, 0.0866f, -0.05f, 0.65f, 0.0866f, -0.05f, 0.0f, 0.1f, -0.0f, 0.65f, 0.1f, -0.0f, 0.0f, 0.0866f, 0.05f, 0.65f, 0.0866f, 0.05f, 0.0f, 0.05f, 0.0866f, 0.65f, 0.05f, 0.0866f, 0.0f, 0.0f, 0.1f, 0.65f, 0.0f, 0.1f, 0.0f, -0.05f, 0.0866f, 0.65f, -0.05f, 0.0866f, 0.0f, -0.0866f, 0.05f, 0.65f, -0.0866f, 0.05f, 0.0f, -0.1f, 0.0f, 0.65, -0.1f, 0.0f, 0.0f, -0.0866f, -0.05f, 0.65f, -0.0866f, -0.05f, 0.0f, -0.05f, -0.0866f, 0.65f, -0.05f, -0.0866f, 0.0f, 0.0f, -0.1f, 0.65f, 0.05f, -0.0866f, 0.65f, 0.0866f, -0.05f, 0.65f, 0.1f, -0.0f, 0.65f, 0.0866f, 0.05f, 0.65f, 0.05f, 0.0866f, 0.65f, 0.0f, 0.1f, 0.65f, -0.05f, 0.0866f, 0.65f, -0.0866f, 0.05f, 0.65, -0.1, 0.0, 0.65, -0.0866, -0.05, 0.65f, -0.05f, -0.0866f, 0.65f, -0.05f, -0.0866f, 0.0f, -0.0866f, -0.05f, 0.0f, -0.1f, 0.0f, 0.0f, -0.0866f, 0.05f, 0.0f, -0.05f, 0.0866f, 0.0f, 0.0f, 0.1f, 0.0, 0.05f, 0.0866, 0.0f, 0.0866, 0.05, 0.0f, 0.1f, -0.0f, 0.0f, 0.0866f, -0.05f, 0.0f, 0.05f, -0.0866f, 0.0f, 0.0f, -0.1f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.2f, 0.65, -0.1, 0.1732, 0.65, -0.1732, 0.1, 0.65f, -0.2f, 0.0f, 0.65f, -0.1732f, -0.1f, 0.65f, -0.1, -0.1732f, 0.65f, -0.0f, -0.2f, 0.65f, 0.1f, -0.1732, 0.65f, 0.1732f, -0.1f, 0.65f, 0.2f, -0.0f, 0.65f, 0.1732f, 0.1f, 0.65f, 0.1f, 0.1732f, 0.65f
	};
	//std::vector<unsigned int> indices;
	unsigned int indices[198] = {
		0, 1, 2, 1, 3, 2, 2, 3, 4, 3, 5, 4, 4, 5, 6, 5, 7, 6, 6, 7, 8, 7, 9, 8, 8, 9, 10, 9, 11, 10,
		10, 11, 12, 11, 13, 12, 12, 13, 14, 13, 15, 14, 14, 15, 16, 15, 17, 16, 16, 17, 18, 17, 19, 18, 18,
		19, 20, 19, 21, 20, 20, 21, 22, 21, 23, 22, 22, 23, 0, 23, 1, 0, 24, 25, 35, 35, 25, 34, 34, 25, 33,
		33, 25, 32, 32, 25, 31, 31, 25, 30, 30, 25, 29, 29, 25, 28, 28, 25, 27, 27, 25, 26, 38, 39, 37, 37,
		39, 36, 36, 39, 47, 47, 39, 46, 46, 39, 45, 45, 39, 44, 44, 39, 43, 43, 39, 42, 42, 39, 41, 41, 39,
		40, 58, 57, 59, 59, 57, 60, 60, 57, 49, 49, 57, 50, 50, 57, 51, 51, 57, 52, 52, 57, 53, 53, 57, 54,
		54, 57, 55, 55, 57, 56, 48, 49, 50, 48, 50, 51, 48, 51, 52, 48, 52, 53, 48, 53, 54, 48, 54, 55, 48,
		55, 56, 48, 56, 57, 48, 57, 58, 48, 58, 59, 48, 59, 60, 48, 60, 49
	};
	unsigned int VAO{ 0 }, VBO{ 0 }, EBO{ 0 };

	void _setup();
};

class PoreNetwork {
public:
	PoreNetwork();
	~PoreNetwork();

	PoreNetwork(const std::vector<float>& vertices, const std::vector<unsigned int>& indices);

	void draw();

	void clean();

private:
	unsigned int VAO{ 0 }, VBO{ 0 }, EBO{ 0 };
	std::vector<float> vertices;
	std::vector<unsigned int> indices;

	void _setup();
};

class Cylinder {
public:
	Cylinder() {};
	~Cylinder() {};

	Cylinder(const Vec3 base, const Vec3 direction, const double height, const double radius);

	void draw();

	void clean();

private:
	unsigned int VAO{ 0 }, VBO{ 0 }, EBO{ 0 };
	std::vector<float> vertices;
	std::vector<unsigned int> indices;
	double PI = 3.14159265358979323846;
	int res = 32;
	void _setup();
};

// -----------------------------------------------------------------------------------
// @brief class that creates a simple line for rendering between two points
class LineModel {
public:

	~LineModel() {};

	LineModel(const glm::vec3& pt1, const glm::vec3& pt2) : pt1(pt1), pt2(pt2), vertices({
		pt1.x, pt1.y, pt1.z,
		pt2.x, pt2.y, pt2.z
		}) {
		_setup();
	};

	LineModel(const Eigen::Vector3f& p1, const Eigen::Vector3f& p2) {

		vertices[0] = p1.x();
		vertices[1] = p1.y();
		vertices[2] = p1.z();
		vertices[3] = p2.x();
		vertices[4] = p2.y();
		vertices[5] = p2.z();

		_setup();
	};

	LineModel(const Vec3& p1, const Vec3& p2) {

		vertices[0] = p1.x;
		vertices[1] = p1.y;
		vertices[2] = p1.z;
		vertices[3] = p2.x;
		vertices[4] = p2.y;
		vertices[5] = p2.z;

		_setup();

	};

	void draw();

	void clean();

	void set_vertices(const glm::vec3& pt1, const glm::vec3& pt2);

	void set_vertices(const Vec3& newpt1, const Vec3& newpt2);

private:

	glm::vec3 pt1;
	glm::vec3 pt2;

	std::array<float, 6> vertices;

	unsigned int VAO{ 0 }, VBO{ 0 };

	void _setup();
};

// -----------------------------------------------
class Ellipsoid {
public:
	Ellipsoid(const Eigen::Vector3f& center, const Eigen::Matrix3f& Rot, float r1, float r2, float r3);
	~Ellipsoid() { _clean(); };

	void draw();

	std::unique_ptr<LineModel> xAxis;
	std::unique_ptr<LineModel> yAxis;
	std::unique_ptr<LineModel> zAxis;

private:
	std::vector<float> vertices;

	unsigned int VAO{ 0 }, VBO{ 0 };

	int stacks = 20;
	int slices = 20;

	void _setup();
	
	void _clean();
};

// ==========================================================
// Sphere Class Definition
// ==========================================================

class Sphere {
public:
    // Constructor matches the signature we used in AnisotropySource
    Sphere(float radius = 1.0f, int sectorCount = 36, int stackCount = 18);
	~Sphere(){clean();};

    void draw(); 
    void set_position(const Vec3& pos);
	void clean();
	void update_radius(float newRadius);

private:
    float radius;
    int sectorCount;
    int stackCount;

    // Standard raw data containers for OpenGL
    std::vector<float> vertices;
    std::vector<float> normals;
    std::vector<unsigned int> indices;

	// opengl buffers
	unsigned int VAO{ 0 }, VBO{ 0 }, normalsVBO{ 0 }, EBO{ 0 };

    // Helper function to generate the geometry
    void create();
    
    // Helper to push data to the GPU (adapt to your engine's specific VAO/VBO setup)
    void setup(); 
};

#endif MODEL_H