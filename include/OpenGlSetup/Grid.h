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

		_setup();
	} ;

	void draw() {
		glBindVertexArray(VAO);
		glDrawArrays(GL_TRIANGLES, 0, 6);
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
		// Attribute-less draw: the vertex shader generates the quad positions
		// from gl_VertexID (its own Pos[]/Indices[] arrays), so we only need a
		// valid, non-zero VAO bound at draw time - no VBO, no attribute pointers.
		// (Enabling an attribute bound to an empty VBO would fetch out of bounds
		// the moment the VS declares a matching input, so we deliberately omit it.)
		glGenVertexArrays(1, &VAO);
	};
};


#endif