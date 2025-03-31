#ifndef MESH_H
#define MESH_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <Visualize/loadMesh.h>

class Mesh {
	
public:

	std::vector<float> vertices;
	std::vector<unsigned int> indices;
	std::vector<float> normals;
	
	// this VAO is going to be used to the draw loop so its public
	unsigned int VAO{ 0 };
	
	// default constructor
	Mesh() {};
	
	// Constructor
	Mesh(std::vector<float> verts, std::vector<unsigned int> ids, std::vector<float> norms) 
		: vertices(verts), indices(ids), normals(norms) {
		
		std::cout << "Creating Mesh" << std::endl;
		_create_mesh();
		
	}

	// constructor overloading for loading from stl
	Mesh(const std::string fileName){
		

		std::cout << "loading Mesh" << std::endl;
		loadStlMesh(fileName, vertices, indices, normals);

		std::cout << "Creating Mesh" << std::endl;
		_create_mesh();

	}
	
	void draw(){
		// draw mesh
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);
	}

	void draw_edges() {
		// 1. Draw the filled mesh
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glEnable(GL_DEPTH_TEST);

		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);

		// 2. Enable polygon offset to avoid Z-fighting
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(-1.0f, -1.0f);

		// 3. Draw the edges as wireframe
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);

		// 4. Restore settings
		glDisable(GL_POLYGON_OFFSET_FILL);
		//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);  // Reset to default
	}

	void deleteMesh() {
		// Delete the vertex buffer, element buffer, and vertex array object
		glDeleteBuffers(1, &VBO);
		glDeleteBuffers(1, &NBO);
		glDeleteBuffers(1, &EBO);
		glDeleteVertexArrays(1, &VAO);
	}

private:
	unsigned int VBO{ 0 }, EBO{ 0 }, NBO{ 0 };
	
	void _create_mesh(){

		// create buffers
		glGenVertexArrays(1, &VAO);
		// bind vertex array
		glBindVertexArray(VAO);

		// vertex data
		glGenBuffers(1, &VBO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

		// normal data
		glGenBuffers(1, &NBO);
		glBindBuffer(GL_ARRAY_BUFFER, NBO);
		glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(float), normals.data(), GL_STATIC_DRAW);

		glGenBuffers(1, &EBO); 
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

		// set vertex attribute pointers
		// vertex locations
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		
		// normals
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, NBO);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)(3 * sizeof(float)));

		// unbind VAO
		glBindVertexArray(0);
	}
	
};

#endif