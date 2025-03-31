#ifndef UNIFORMMANAGER_H // include guard
#define UNIFORMMANAGER_H

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include "shader.h"

class UniformManager {

public:
	UniformManager() {};
	~UniformManager() {};

	// add a uniform to the manager
	void add_uniform(Shader& shader, const std::string& uniformTag) {

		// get location of uniform
		GLuint loc = glGetUniformLocation(shader.programID, uniformTag.c_str());

		if (loc != -1) {
			uniforms[shader.programID][uniformTag] = loc;
		}
		else {
			std::cerr << "Uniform: " << uniformTag << " not found in shader: " << shader.programID << std::endl;
		}

	};

	// function to set a uniform, overloading for different uniform types
	// pass a float
	void setUniform(Shader& shader, const std::string& uniformTag, const float value) {
		GLuint loc = uniforms[shader.programID][uniformTag];
		if (loc != -1) {
			glUniform1f(loc, value);
		}
	}
	// pass an int
	void setUniform(Shader& shader, const std::string& uniformTag, const int value) {
		GLuint loc = uniforms[shader.programID][uniformTag];
		if (loc != -1) {
			glUniform1i(loc, value);
		}
	}
	// pass a vec3
	void setUniform(Shader& shader, const std::string& uniformTag, const glm::vec3& value) {
		GLuint loc = uniforms[shader.programID][uniformTag];
		if (loc != -1) {
			glUniform3fv(loc, 1, glm::value_ptr(value));
		}
	};
	// three floats
	void setUniform(Shader& shader,
		const std::string& uniformTag,
		const float& value1,
		const float& value2,
		const float& value3
	) {

		GLuint loc = uniforms[shader.programID][uniformTag];
		if (loc != -1) {
			glUniform3f(loc, value1, value2, value3);
		}
	};

	void setUniform(Shader& shader,
		const std::string& uniformTag,
		const float& value1,
		const float& value2,
		const float& value3,
		const float& value4
	) {

		GLuint loc = uniforms[shader.programID][uniformTag];
		if (loc != -1) {
			glUniform4f(loc, value1, value2, value3, value4);
		}
	};

	// pass a flatten array
	void setUniform(Shader& shader, const std::string& uniformTag, const std::vector<float>& value) {
		GLuint loc = uniforms[shader.programID][uniformTag];
		if (loc != -1) {
			glUniform3fv(loc, 1, glm::value_ptr(value));
		}
	}
	// pass a vec4
	void setUniform(Shader& shader, const std::string& uniformTag, const glm::vec4& value) {
		GLuint loc = uniforms[shader.programID][uniformTag];
		if (loc != -1) {
			glUniform4fv(loc, 1, glm::value_ptr(value));
		}
	}
	void setUniform(Shader& shader, const std::string& uniformTag, const glm::vec2& value) {
		GLuint loc = uniforms[shader.programID][uniformTag];
		if (loc != -1) {
			glUniform2fv(loc, 1, glm::value_ptr(value));
		}
	}
	// pass a mat4
	void setUniform(Shader& shader, const std::string& uniformTag, const glm::mat4& value) {
		GLuint loc = uniforms[shader.programID][uniformTag];
		if (loc != -1) {
			glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(value));
		}
	}

private:

	// create an unorderded map of structure [shader id nr] [uniform name] [uniform loc]
	std::unordered_map<GLuint, std::unordered_map<std::string, GLuint>> uniforms;

};

#endif