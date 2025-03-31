#ifndef CAMERA_H
#define CAMERA_H
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <Misc/Quaternion.h>

class Camera {

public: 

	// constructor
	Camera(GLFWwindow* glWindow, double split, const Vec3 targetCamera);
	void update();
	void update_scroll(float yoffset);
	Vec3 sphere_projection(double xpos, double ypos, Vec3 target);
	Vec3 arc(double xpos, double ypos);
	void set_window_size(int newWidth, int newHeight);

	int width{ 0 }, height{ 0 };
	GLFWwindow* window;
	double split;

	glm::mat4 viewMatrix;
	glm::mat4 projectionMatrix;

	glm::vec3 position;

	// Initial horizontal angle : toward -Z
	float horizontalAngle;
	// Initial vertical angle : none
	float verticalAngle;

	float FoV;
	float speed; // units / second
	float mouseSpeed;
	float fovSpeed;

	// trackball
	float sphereRadius = 300.0f;
	Vec3 targetCamera;

	Quaternion q;
	Quaternion qPrev;

	double xCurr{ 0 }, yCurr{ 0 };
	double xPrev{ 0 }, yPrev{ 0 };

	bool rightPressed{ false };
	bool leftPressed{ false };
	int arcball_on = false;

private:
	glm::vec3 up;
	glm::vec3 direction;
	glm::vec3 right;

	void to_viewport(double& xpos, double& ypos);
};

glm::mat4 quaternion_2_matrix(const Quaternion q);

#endif