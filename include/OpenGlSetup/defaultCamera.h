#ifndef DEAFULTCAMERA_H
#define DEFAULTCAMERA_H
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <Math/Quaternion.h>
#include <Math/Vec.h>

class defaultCamera {

public: 

	defaultCamera() {};
	~defaultCamera() {};

	// constructor
	defaultCamera(GLFWwindow* window, double split,
		glm::vec3 cameraPos, glm::vec3 cameraTarget, float speed);

	void update();
	void set_window_size(int newWidth, int newHeight);

	int width{ 0 }, height{ 0 };
	GLFWwindow* window = nullptr;
	
	glm::mat4 viewMatrix = glm::mat4(1.0f);
	glm::mat4 projectionMatrix = glm::mat4(1.0f);

	float horizontalAngle;
	float verticalAngle;
	float mouseSpeed;
	float fovSpeed;

	double xCurr{ 0 }, yCurr{ 0 };
	double xPrev{ 0 }, yPrev{ 0 };

	glm::vec3 position = glm::vec3(0.0f, 0.0f, -10.0f);

	glm::vec3 target = glm::vec3(0.0f, 0.0f, 0.0f);

	void scroll_setup(GLFWwindow* window);

private:
	
	glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);
	//glm::vec3 up = Vec3(0.0f, 1.0f, 0.0f);
	
	//Vec3 direction = Vec3(0.0f, 0.0f, -1.0f);
	glm::vec3 direction = glm::vec3(0.0f, 0.0f, -1.0f);
	
	glm::vec3 right = glm::vec3(1.0f, 0.0f, 0.0f);
	//Vec3 right = Vec3(1.0f, 0.0f, 0.0f);
	
	float speed{ 0.0f };

	float FoV{ 45.0f };

	double split{ 1.0f };

	void to_viewport(double& xpos, double& ypos);

	void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

	bool in_viewport();
};

glm::mat4 quaternion_2_matrix(const Quaternion q);

#endif