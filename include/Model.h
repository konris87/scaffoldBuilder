#ifndef MODEL_H
#define MODEL_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <vtkSTLReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkCellArrayIterator.h>
#include <vtkPolyDataNormals.h>
#include <vtkMassProperties.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
//#include <vtkFillHolesFilter.h>
#include <vector>
#include <set>
#include <iostream>
#include <filesystem>
#include "Math/Quaternion.h"
#include "Math/Vec.h"

class Model {
	
public:

	unsigned int VAO{ 0 }, edgeVAO{ 0 };

	Model() {};
	~Model() {};

	Model(std::string& modelFile) {
		
		// check model file extension
		std::filesystem::path filePath = modelFile;

		std::cout << filePath << std::endl;

		if (filePath.extension() == ".stl") {
			vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
			reader->SetFileName(modelFile.c_str());
			reader->Update();

			// get header
			std::string header{ "" };
			header = reader->GetHeader();

			// get polydata, points and cells
			polyData = vtkSmartPointer<vtkPolyData>::New();
			polyData = reader->GetOutput();
		}
		else if (filePath.extension() == ".vtk") {
			vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
			reader->SetFileName(modelFile.c_str());
			reader->Update();

			// get polydata, points and cells
			polyData = vtkSmartPointer<vtkPolyData>::New();
			polyData = reader->GetOutput();
		}
		else {
			std::cerr << "Invalid format. Provide a proper .stl or .vtk mesh" << std::endl;
		}
	
		setup_data();
		setup_mesh();
		setup_edges();
	};

	Model(vtkSmartPointer<vtkPolyData>& inputData) {

		vtkNew<vtkCleanPolyData> cleaner;
		cleaner->SetInputData(inputData);
		cleaner->Update();

		vtkNew<vtkTriangleFilter> triangleFilter;
		triangleFilter->SetInputConnection(cleaner->GetOutputPort());
		triangleFilter->Update();

		//vtkNew<vtkFillHolesFilter> fillHoles;
		//fillHoles->SetInputConnection(triangleFilter->GetOutputPort());
		//fillHoles->SetHoleSize(1e6);
		//fillHoles->Update();

		vtkNew<vtkCleanPolyData> finalClean;
		finalClean->SetInputConnection(triangleFilter->GetOutputPort());
		finalClean->Update();

		polyData = finalClean->GetOutput();

		vtkNew<vtkMassProperties> massProperties;
		massProperties->SetInputData(polyData);
		massProperties->Update();

		volume = massProperties->GetVolume();
		std::cout << "Volume: " << volume << std::endl;

		setup_data();
		setup_mesh();
		setup_edges();
	};

	void draw() {
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	};

	void draw_edges() {
		glBindVertexArray(edgeVAO);
		glDrawElements(GL_LINES, edgeIndices.size(), GL_UNSIGNED_INT, 0);
		//glBindVertexArray(1);
	};

	void clean() {
		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		glDeleteVertexArrays(1, &VAO);
		glDeleteBuffers(1, &VBO);
		glDeleteBuffers(1, &normalsVBO);
		glDeleteBuffers(1, &EBO);
		glDeleteVertexArrays(1, &edgeVAO);
		glDeleteBuffers(1, &edgeVBO);
		glDeleteBuffers(1, &edgeEBO);
		VAO, VBO, normalsVBO, EBO, edgeVAO, edgeVBO, edgeEBO = 0;
	};

	double xMin{ 0.0 }, xMax{ 0.0 }, yMin{ 0.0 }, yMax{ 0.0 }, zMin{ 0.0 }, zMax{ 0.0 };

	//void get_bounds(double& xMin, double& xMax, double& yMin, double& yMax, double& zMin, double& zMax) {};

	void get_details(int& faceNr, int& vertexNr, int& edgeNr, double& vol, double& por) {

		//std::cout << cellSize << vertexSize << edgeSize << std::endl;
		faceNr = cellSize;
		vertexNr = vertexSize;
		edgeNr = edgeSize;
		vol = volume;
		por = porosity;
	}

private:
	vtkSmartPointer<vtkPolyData> polyData;
	std::vector<float> vertices;
	std::vector<unsigned int> indices;
	std::vector<unsigned int> edgeIndices;
	std::vector<float> vertexNormals;
	unsigned int VBO{ 0 }, EBO{ 0 }, normalsVBO{ 0 };
	unsigned int edgeVBO{ 0 }, edgeEBO{ 0 };

	int vertexSize{ 0 };
	int cellSize{ 0 };
	int edgeSize{ 0 };
	double volume{ 0.0 };
	double porosity{ 0.0 };

	void setup_data() {

		//std::cout << "setup_data" << std::endl;

		vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
		vtkSmartPointer<vtkCellArray> cells = polyData->GetPolys();

		// update also bounds
		double bounds[6];
		polyData->GetBounds(bounds);
		xMin = bounds[0];
		xMax = bounds[1];
		yMin = bounds[2];
		yMax = bounds[3];
		zMin = bounds[4];
		zMax = bounds[5];

		//std::cout << "volume" << std::endl;
		// get volume

		vtkSmartPointer<vtkFieldData> fieldData = polyData->GetFieldData();
		if (fieldData) {

			// Read Porosity
			vtkDoubleArray* porosityArray = vtkDoubleArray::SafeDownCast(fieldData->GetArray("Porosity"));
			if (porosityArray && porosityArray->GetNumberOfTuples() > 0)
			{
				porosity = porosityArray->GetValue(0);
				std::cout << "Scaffold Porosity: " << porosity << std::endl;
			}
			else
			{
				std::cerr << "Porosity data not found!" << std::endl;
			}
		}

		double center[3];
		polyData->GetCenter(center);
		//std::cout << "Actual center: " << center[0] << " " << center[1] << center[2] << std::endl;

		// get vertex data and pass them to the vector of vertices
		for (vtkIdType i{ 0 }; i < points->GetNumberOfPoints(); i++) {
			vertexSize++;
			double pt[3];
			points->GetPoint(i, pt);
			vertices.push_back(pt[0]);
			vertices.push_back(pt[1]);
			vertices.push_back(pt[2]);
		}

		// get face indices and edge indices
		// use a set to avoid duplicate edges
		std::set<std::pair<unsigned int, unsigned int>> edgeSet;

		auto iter = vtk::TakeSmartPointer(cells->NewIterator());
		for (iter->GoToFirstCell(); !iter->IsDoneWithTraversal(); iter->GoToNextCell())
		{
			cellSize += 1;
			vtkSmartPointer<vtkIdList> cell = vtkSmartPointer<vtkIdList>::New();
			iter->GetCurrentCell(cell);

			unsigned int v1 = static_cast<unsigned int>(cell->GetId(0));
			unsigned int v2 = static_cast<unsigned int>(cell->GetId(1));
			unsigned int v3 = static_cast<unsigned int>(cell->GetId(2));

			indices.push_back(v1);
			indices.push_back(v2);
			indices.push_back(v3);
			//std::cout << cell->GetId(0) << " " << cell->GetId(1) << " " << cell->GetId(2) << std::endl;

			// push also edges, use a set to not save duplicates
			std::pair<unsigned int, unsigned int> e1 = { std::min(v1, v2), std::max(v1, v2) };
			std::pair<unsigned int, unsigned int> e2 = { std::min(v2, v3), std::max(v2, v3) };
			std::pair<unsigned int, unsigned int> e3 = { std::min(v3, v1), std::max(v3, v1) };

			if (edgeSet.insert(e1).second) {
				edgeIndices.push_back(e1.first);
				edgeIndices.push_back(e1.second);
			}

			if (edgeSet.insert(e2).second) {
				edgeIndices.push_back(e2.first);
				edgeIndices.push_back(e2.second);
			}

			if (edgeSet.insert(e3).second) {
				edgeIndices.push_back(e3.first);
				edgeIndices.push_back(e3.second);
			}
		}

		// get also the mesh normals for now we use the face normals
		//std::cout << "Estimating norms" << std::endl;
		vtkNew<vtkPolyDataNormals> norms;

		norms->SetInputData(polyData);
		norms->ComputePointNormalsOn();
		//norms->ComputeCellNormalsOn();
		//norms->ConsistencyOn();
		//norms->AutoOrientNormalsOn();
		norms->Update();

		// Get normals from the updated polydata
		vtkSmartPointer<vtkDataArray> normalData = norms->GetOutput()->GetPointData()->GetNormals();
		if (!normalData) {
			std::cerr << "Error: No normals found!" << std::endl;
			return;
		}
		else {
			//std::cout << normalData->GetNumberOfTuples() << std::endl;
			for (vtkIdType i{ 0 }; i < normalData->GetNumberOfTuples(); i++) {
				double normal[3];

				normalData->GetTuple(i, normal);

				vertexNormals.push_back(normal[0]);
				vertexNormals.push_back(normal[1]);
				vertexNormals.push_back(normal[2]);
			}
		}

		//std::cout << "Loaded vertices: " << vertices.size() / 3 << std::endl;
		//std::cout << "Loaded vertex normals: " << vertexNormals.size() << std::endl;
		//std::cout << "Loaded faces (triangles): " << indices.size() / 3 << std::endl;
		//std::cout << "Edge count: " << edgeIndices.size() / 2 << std::endl;

		vertexSize = vertices.size() / 3;
		cellSize = indices.size() / 3;
		edgeSize = edgeIndices.size() / 2;

	}

	void setup_mesh() {
	
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

	void setup_edges() {

		//std::cout << "setup edges" << std::endl;

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
		glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

		// bind elementa array buffer and send data
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edgeEBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, edgeIndices.size() * sizeof(unsigned int), edgeIndices.data(), GL_STATIC_DRAW);

		// position attribute
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);

		//glBindVertexArray(0);

	};
};

class CutPlane {

public:
	CutPlane(float size) : size(size) {
		_setup();
	};
	~CutPlane() {
		clean();
	};

	CutPlane(float normal[3], float& offset) : offset(offset) {
		_setup();
	};

	CutPlane(Vec3 origin){
		
		center = glm::vec3(origin.x, origin.y, origin.z);

		//std::cout << "Passed center: " << origin.x << " " << origin.y << " " << origin.z << std::endl;

		// change coordinates of the plane
		for (int i{ 0 }; i < 4; i++) {
			vertices[i] += origin.x;
			vertices[i+1] += origin.y;
			vertices[i+2] += origin.z;
		}

		_setup();
	}

	void draw(){
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}

	void clean() {
		glDeleteVertexArrays(1, &VAO);
		glDeleteBuffers(1, &VBO);
		glDeleteBuffers(1, &EBO);
	};

	void updateModelMatrix() {

		glm::vec3 newNormal = glm::normalize(glm::vec3(normal[0], normal[1], normal[2]));
		prevNormal = glm::normalize(prevNormal);

		std::cout << "Previous: " << prevNormal[0] << " " << prevNormal[1] << " " << prevNormal[2] << std::endl;
		std::cout << "Current: " << newNormal[0] << " " << newNormal[1] << " " << newNormal[2] << std::endl;

		if (glm::length(newNormal - prevNormal) < 1e-6f) {
			return;
		}

		glm::vec3 axis = glm::cross(prevNormal, newNormal);
		float angle = acos(glm::clamp(glm::dot(prevNormal, newNormal), -1.0f, 1.0f));

		if (glm::length(axis) < 1e-6f) {
			modelMatrix = glm::mat4(1.0f);  // No rotation needed
		}
		

		rotMatrix = glm::rotate(glm::mat4(1.0f), angle, glm::normalize(axis)) * rotMatrix;

	};

	void updateTranslation() {
		tMatrix = glm::translate(glm::mat4(1.0f), normal * offset);
	}

	glm::vec3 center = { 0.0f, 0.0f, 0.0f };
	glm::vec3 normal = { 0.0f, 0.0f, 1.0f };
	glm::vec3 prevNormal = { 0.0f, 0.0f, 1.0f };
	glm::mat4 modelMatrix = glm::mat4(1.0);
	glm::mat4 rotMatrix = glm::mat4(1.0);
	glm::mat4 tMatrix = glm::mat4(1.0);
	glm::mat4 initMatrix = glm::mat4(1.0);
	float offset{ 0.0f };
	float prevOffset{ 0.0f };

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
		0, 1, 2, 0, 2, 3
	};

	void _setup() {
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
};

class BBox {
public:
	BBox() { _setup(); };
	~BBox() { clean(); };

	BBox(float xMin, float xMax, float yMin, float yMax, float zMin, float zMax){
		
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

	void draw() {
		glBindVertexArray(VAO);
		glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	};

	void clean() {
		glDeleteVertexArrays(1, &VAO);
		glDeleteBuffers(1, &VBO);
	};

private:
	std::vector<float> vertices;
	unsigned int indices[24] = {
		0, 1, 1, 2, 2, 3, 3, 0,
		4, 5, 5, 6, 6, 7, 7, 4,
		0, 4, 3, 7, 2, 6, 1, 5
	};
	unsigned int VAO, VBO, EBO;

	void _setup() {
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
};

class Arrow {
public:

	Arrow() {
		_setup();
	};
	~Arrow() {
		clean();
	};

	void draw() {
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, 198, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}

	void clean() {
		glDeleteVertexArrays(1, &VAO);
		glDeleteBuffers(1, &VBO);
		glDeleteBuffers(1, &EBO);
	};

private:
	float vertices[183] = {
	0.0, -0.1, 0.65, 0.0, -0.1, 0.0, 0.05, -0.0866, 0.65, 0.05, -0.0866, 0.0, 0.0866, -0.05, 0.65, 0.0866, -0.05, 0.0, 0.1, -0.0, 0.65, 0.1, -0.0, 0.0, 0.0866, 0.05, 0.65, 0.0866, 0.05, 0.0, 0.05, 0.0866, 0.65, 0.05, 0.0866, 0.0, 0.0, 0.1, 0.65, 0.0, 0.1, 0.0, -0.05, 0.0866, 0.65, -0.05, 0.0866, 0.0, -0.0866, 0.05, 0.65, -0.0866, 0.05, 0.0, -0.1, 0.0, 0.65, -0.1, 0.0, 0.0, -0.0866, -0.05, 0.65, -0.0866, -0.05, 0.0, -0.05, -0.0866, 0.65, -0.05, -0.0866, 0.0, 0.0, -0.1, 0.65, 0.05, -0.0866, 0.65, 0.0866, -0.05, 0.65, 0.1, -0.0, 0.65, 0.0866, 0.05, 0.65, 0.05, 0.0866, 0.65, 0.0, 0.1, 0.65, -0.05, 0.0866, 0.65, -0.0866, 0.05, 0.65, -0.1, 0.0, 0.65, -0.0866, -0.05, 0.65, -0.05, -0.0866, 0.65, -0.05, -0.0866, 0.0, -0.0866, -0.05, 0.0, -0.1, 0.0, 0.0, -0.0866, 0.05, 0.0, -0.05, 0.0866, 0.0, 0.0, 0.1, 0.0, 0.05, 0.0866, 0.0, 0.0866, 0.05, 0.0, 0.1, -0.0, 0.0, 0.0866, -0.05, 0.0, 0.05, -0.0866, 0.0, 0.0, -0.1, 0.0, 0.0, 0.0, 1.0, 0.0, 0.2, 0.65, -0.1, 0.1732, 0.65, -0.1732, 0.1, 0.65, -0.2, 0.0, 0.65, -0.1732, -0.1, 0.65, -0.1, -0.1732, 0.65, -0.0, -0.2, 0.65, 0.1, -0.1732, 0.65, 0.1732, -0.1, 0.65, 0.2, -0.0, 0.65, 0.1732, 0.1, 0.65, 0.1, 0.1732, 0.65
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

	void _setup() {
		
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
	}
};

#endif MODEL_H