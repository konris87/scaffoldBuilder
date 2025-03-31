#ifndef SHADER_H
#define SHADER_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class Shader {

public:	

	GLuint programID{ 0 };

	Shader() {};
	~Shader() {};

	Shader(
	const char* vShaderFile,
	const char* fShaderFile,
	const char* gShaderFile
	)   {
		
		std::string vShaderCodeStr = _read_source(vShaderFile);
		std::string fShaderCodeStr = _read_source(fShaderFile);
		
		const char* vShaderCode = vShaderCodeStr.c_str();
		const char* fShaderCode = fShaderCodeStr.c_str();
		
		vertexShader = glCreateShader(GL_VERTEX_SHADER);
		fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
		
		_compile_shader(vertexShader, vShaderCode, "VERTEX");
		_compile_shader(fragmentShader, fShaderCode, "FRAGMENT");
		
		if (gShaderFile){
			std::string gShaderCodeStr = _read_source(gShaderFile);
			const char* gShaderCode = gShaderCodeStr.c_str();
			geometryShader = glCreateShader(GL_GEOMETRY_SHADER);
			_compile_shader(geometryShader, gShaderCode, "GEOMETRY");
		}		
		
	// create program
	programID = glCreateProgram();
	glAttachShader(programID, vertexShader);
	glAttachShader(programID, fragmentShader);
	if (gShaderFile){
		glAttachShader(programID, geometryShader);
	}
	
	glLinkProgram(programID);
	_check_progError();
	
	// delete shaders
	_delete_shaders();
	
	};

	void use(){
		glUseProgram(programID);
	}

	void deleteProg(){
		glDeleteProgram(programID);
	}

	// unlink and delete
	void _delete_shaders() {
		glDeleteShader(vertexShader);
		glDeleteShader(fragmentShader);
		if (geometryShader) {
			glDeleteShader(geometryShader);
		}
	}
	
private:
	GLuint vertexShader{ 0 };
	GLuint fragmentShader{ 0 };
	GLuint geometryShader{ 0 };
	
	std::string _read_source(const char* filePath) {
    
		std::string code;
		std::ifstream shaderFile;

		// throw exceptions
		shaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

		try
		{
			shaderFile.open(filePath);
			std::stringstream shaderStream;
			shaderStream << shaderFile.rdbuf();

			// close files
			shaderFile.close();

			// convert to string
			code = shaderStream.str();
			return code;
		}
		catch (std::ifstream::failure e)
		{
			std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl; 
		}	
	}
	
	// Compile Shader
	void _compile_shader(const unsigned int shader, const char* src, std::string type){
		
		glShaderSource(shader, 1, &src, NULL);
		glCompileShader(shader);
		
		_check_shaderError(shader, type);
	}
	
    // Check for compile errors
	void _check_shaderError(unsigned int shader, std::string type){
		int success;
		char infoLog[512];
		
		glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
		if (!success) {
			glGetShaderInfoLog(shader, 512, NULL, infoLog);
			if (type=="VERTEX"){
				std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
			}
			if (type == "FRAGMENT") {
				std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
			}
			if (type == "GEOMETRY") {
				std::cout << "ERROR::SHADER::GEOMETRY::COMPILATION_FAILED\n" << infoLog << std::endl;
			}
		}
		
	}
	
	void _check_progError(){
		
		int success;
		char infoLog[512];
				
		glGetProgramiv(programID, GL_LINK_STATUS, &success);
		if (!success) {
			glGetProgramInfoLog(programID, 512, NULL, infoLog);
			std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
		}
	}
	

};

#endif