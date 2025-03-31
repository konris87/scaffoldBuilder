#ifndef GRID_H
#define GRID_H
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <vector>


class Grid{

public:
	Grid(){};
	~Grid(){};

	enum gridMode { XY };

	Grid(gridMode mode) : mode(mode) {
		
		if (mode == XY) {

			// top left corner
			vertices.push_back(0.0f);
			vertices.push_back(0.0f);
			vertices.push_back(0.0f);
			
			//// 
			//vertices.push_back(1.0f);
			//vertices.push_back(1.0f);
			//vertices.push_back(0.0f);
			//// third
			//vertices.push_back(1.0f);
			//vertices.push_back(-1.0f);
			//vertices.push_back(0.0f);
			//// fourth
			//vertices.push_back(-1.0f);
			//vertices.push_back(-1.0f);
			//vertices.push_back(0.0f);

			_setup();
		}	
	} ;

	void draw() {
		glBindVertexArray(VAO);
		glDrawArrays(GL_POINTS, 0, 1);
		glBindVertexArray(0);
	};

	void clean() {
		glDeleteVertexArrays(1, &VAO);
		glDeleteBuffers(1, &VBO);
	};

private:
	gridMode mode;
	std::vector<float> vertices;
	unsigned int VAO{ 0 }, VBO{ 0 };

	void _setup() {
		
		// create VAO
		glGenVertexArrays(1, &VAO);
		glBindVertexArray(VAO);

		glGenBuffers(1, &VBO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

		// position attribute for grid shader program
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);
		glBindVertexArray(0);
	};
};


#endif