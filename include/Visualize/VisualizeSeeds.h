#ifndef VISUALIZESEEDS_H
#define VISUALIZESEEDS_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
#include <array>
#include "Math/Vec.h"
#include <algorithm>

class VisualizeSeeds {

public:
	VisualizeSeeds() {};
	~VisualizeSeeds() {};
	VisualizeSeeds(const std::vector<Vec3>& coords, const std::array<float, 6>& bounds) : bounds(bounds), seedCoords(coords) {

		// 1. Calculate automatic size
		float dx = bounds[1] - bounds[0];
		float dy = bounds[3] - bounds[2];
		float dz = bounds[5] - bounds[4];
		float maxDim = std::max({ dx, std::max(dy, dz) });

		initialCalculatedSize = maxDim * 0.01f;
		if (initialCalculatedSize <= 0.0f) initialCalculatedSize = 0.1f;

		_generate_points();

		// 2. Generate model matrices (TRANSLATION ONLY)
		for (unsigned int i{ 0 }; i < seedCoords.size(); i++) {
			glm::mat4 model = glm::mat4(1.0f);
			model = glm::translate(model, glm::vec3(seedCoords[i].x, seedCoords[i].y, seedCoords[i].z));
			// NO SCALING HERE! Leave it to the shader.
			modelMatrices.push_back(model);
		}

		_create();
	};

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

	float initialCalculatedSize = 0.1f;

private:
	std::array<float, 6> bounds = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
	std::vector<float> sphereVertices;
	std::vector<Vec3> seedCoords;
	std::vector<unsigned int> sphereIndices;
	unsigned int VAO{ 0 }, VBO{ 0 }, instanceVBO{ 0 }, EBO{ 0 };
	std::vector<glm::mat4> modelMatrices;
	float sphereSize = 0.1f;
	float radius = 1.0f;
	int sectorCount = 5; // Theta resolution (longitude)
	int stackCount = 5;

	void _generate_points() {

		// reset
		sphereVertices.clear();
		sphereIndices.clear();

		const float PI = 3.14159265359f;
		float sectorStep = 2 * PI / sectorCount;
		float stackStep = PI / stackCount;
		float sectorAngle, stackAngle;

		for (int i = 0; i <= stackCount; ++i) {
			stackAngle = PI / 2.0f - i * stackStep;     // from pi/2 to -pi/2
			float xy = radius * cosf(stackAngle);       // r * cos(u)
			float z = radius * sinf(stackAngle);        // r * sin(u)

			for (int j = 0; j <= sectorCount; ++j) {
				sectorAngle = j * sectorStep;           // from 0 to 2pi

				// vertex position (x, y, z)
				float x = xy * cosf(sectorAngle);
				float y = xy * sinf(sectorAngle);

				sphereVertices.push_back(x);
				sphereVertices.push_back(y);
				sphereVertices.push_back(z);
			}
		}

		// Generate indices connecting the vertices to form triangles
		int k1, k2;
		for (int i = 0; i < stackCount; ++i) {
			k1 = i * (sectorCount + 1);     // beginning of current stack
			k2 = k1 + sectorCount + 1;      // beginning of next stack

			for (int j = 0; j < sectorCount; ++j, ++k1, ++k2) {
				// 2 triangles per sector excluding the first and last stacks (poles)

				// Triangle 1 (k1 -> k2 -> k1+1)
				if (i != 0) {
					sphereIndices.push_back(k1);
					sphereIndices.push_back(k2);
					sphereIndices.push_back(k1 + 1);
				}

				// Triangle 2 (k1+1 -> k2 -> k2+1)
				if (i != (stackCount - 1)) {
					sphereIndices.push_back(k1 + 1);
					sphereIndices.push_back(k2);
					sphereIndices.push_back(k2 + 1);
				}
			}
		}
	}

	void update_size(float newSize) {
		sphereSize = newSize;

		for (unsigned int i = 0; i < seedCoords.size(); i++) {
			glm::mat4 model = glm::mat4(1.0f);
			// Apply translation first, THEN scale
			model = glm::translate(model, glm::vec3(seedCoords[i].x, seedCoords[i].y, seedCoords[i].z));
			model = glm::scale(model, glm::vec3(sphereSize));
			modelMatrices[i] = model;
		}

		// Send the updated matrices to the GPU
		glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
		glBufferData(GL_ARRAY_BUFFER, modelMatrices.size() * sizeof(glm::mat4), modelMatrices.data(), GL_DYNAMIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	void _create() {
		glGenVertexArrays(1, &VAO);
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
		glBufferData(GL_ARRAY_BUFFER, modelMatrices.size() * sizeof(glm::mat4), modelMatrices.data(), GL_DYNAMIC_DRAW);

		std::size_t vec4Size = sizeof(glm::vec4);
		for (unsigned int i{ 0 }; i < 4; i++) {
			glEnableVertexAttribArray(i + 1);
			glVertexAttribPointer(i + 1, 4, GL_FLOAT, GL_FALSE, 4 * (GLsizei)vec4Size, (void*)(i * vec4Size));
			glVertexAttribDivisor(i + 1, 1);
		}
		glBindVertexArray(0);
	};
};
#endif