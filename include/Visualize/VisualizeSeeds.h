#ifndef VISUALIZESEEDS_H
#define VISUALIZESEEDS_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
#include <array>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCellArrayIterator.h>

class VisualizeSeeds {

public:
	VisualizeSeeds() {};
	~VisualizeSeeds() {};
	VisualizeSeeds(std::vector<std::array<double, 3>>& coords) {

		auto sphere = vtkSmartPointer<vtkSphereSource>::New();
		sphere->SetRadius(1.0);
		sphere->SetPhiResolution(10);
		sphere->SetThetaResolution(10);
		sphere->Update();

		vtkSmartPointer<vtkPolyData> spherePoly = sphere->GetOutput();
		vtkSmartPointer<vtkPoints> points = spherePoly->GetPoints();
		vtkSmartPointer<vtkCellArray> polys = spherePoly->GetPolys();

		for (int i{ 0 }; i < points->GetNumberOfPoints(); i++) {
			double p[3];
			points->GetPoint(i, p);
			sphereVertices.push_back(p[0]);
			sphereVertices.push_back(p[1]);
			sphereVertices.push_back(p[2]);
		}

		// get indices
		auto cellIter = vtk::TakeSmartPointer(polys->NewIterator());

		//std::cout << cellNr << std::endl;
		for (cellIter->GoToFirstCell(); !cellIter->IsDoneWithTraversal();
			cellIter->GoToNextCell())
		{
			vtkSmartPointer<vtkIdList> cell = cellIter->GetCurrentCell();
			sphereIndices.emplace_back(cell->GetId(0));
			sphereIndices.emplace_back(cell->GetId(1));
			sphereIndices.emplace_back(cell->GetId(2));
		}

		for (unsigned int i{ 0 }; i < coords.size(); i++) {
			glm::mat4 model = glm::mat4(1.0f);
			model = glm::translate(model, glm::vec3(coords[i][0], coords[i][1], coords[i][2]));
			modelMatrices.push_back(model);
		}

		_create();
	};

	void updateSize(float newSize) {
		sphereSize = newSize;

		for (unsigned int i = 0; i < modelMatrices.size(); i++) {
			modelMatrices[i] = glm::scale(glm::mat4(1.0f), glm::vec3(sphereSize));
		}

		glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
		glBufferData(GL_ARRAY_BUFFER, modelMatrices.size() * sizeof(glm::mat4), modelMatrices.data(), GL_DYNAMIC_DRAW);
	}

	void draw() {
		// Update instance positions if they change
		glBindVertexArray(VAO);
		glDrawElementsInstanced(
			GL_TRIANGLES, sphereIndices.size(), GL_UNSIGNED_INT, 0, modelMatrices.size());
		glBindVertexArray(0);
	};

	void deleteObj() {
		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		glDeleteVertexArrays(1, &VAO);
		glDeleteVertexArrays(1, &instanceVBO);
		glDeleteBuffers(1, &VBO);
	};

private:
	std::vector<float> seeds;
	std::vector<float> sphereVertices;
	std::vector<unsigned int> sphereIndices;
	unsigned int VAO, VBO, instanceVBO, EBO;
	std::vector<glm::mat4> modelMatrices;
	float sphereSize = 0.1f;

	void _create() {

		// this is for the sphere object vertices
		glGenVertexArrays(1, &VAO);
		// bind vertex array
		glBindVertexArray(VAO);

		// vertex data
		glGenBuffers(1, &VBO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(
			GL_ARRAY_BUFFER,
			sphereVertices.size() * sizeof(float),
			sphereVertices.data(),
			GL_STATIC_DRAW);

		// indices
		glGenBuffers(1, &EBO);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(
			GL_ELEMENT_ARRAY_BUFFER,
			sphereIndices.size() * sizeof(unsigned int),
			sphereIndices.data(),
			GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

		// instancing buffer
		glGenBuffers(1, &instanceVBO);
		glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
		glBufferData(GL_ARRAY_BUFFER, modelMatrices.size() * sizeof(glm::mat4), modelMatrices.data(), GL_STATIC_DRAW);

		std::size_t vec4Size = sizeof(glm::vec4);
		for (unsigned int i{ 0 }; i < 4; i++) {
			
			glEnableVertexAttribArray(i + 1);
			glVertexAttribPointer(i + 1, 4, GL_FLOAT, GL_FALSE, 4 * vec4Size, (void*)(i * vec4Size));
			glVertexAttribDivisor(i + 1, 1);
		}
		glBindVertexArray(0);
	};
};	

#endif