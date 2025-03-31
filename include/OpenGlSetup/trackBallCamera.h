#ifndef TRACKBALLCAMERA_H
#define TRACKBALLCAMERA_H
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <imgui.h>
#include <Misc/Quaternion.h>

void ImGui_ImplGlfw_ScrollCallback(GLFWwindow* window, double xoffset, double yoffset);

class trackBallCamera {

public:

	// constructor
	trackBallCamera() {};
	~trackBallCamera() {};
	trackBallCamera(GLFWwindow* glWindow, double split, glm::vec3 targetCamera, glm::vec3 posCamera);
	void update();
	Vec3 sphere_projection(double xpos, double ypos, Vec3 target);
	Vec3 arc(double xpos, double ypos);
	void set_window_size(int newWidth, int newHeight);

	int width{ 0 }, height{ 0 };
	GLFWwindow* window;
	double split;

	glm::mat4 viewMatrix= glm::mat4(1.0f);
	glm::mat4 projectionMatrix = glm::mat4(1.0f);
	glm::vec3 position{ 0.0, 0., 20.0 };

    // Initial horizontal angle : toward -Z
	float horizontalAngle;
	// Initial vertical angle : none
	float verticalAngle;

	float FoV { 45.0f };
	float speed; // units / second
	float mouseSpeed;
	float fovSpeed;

	// trackball
	float sphereRadius = 300.0f;

	Quaternion q;
	Quaternion qPrev;

	double xCurr{ 0 }, yCurr{ 0 };
	double xPrev{ 0 }, yPrev{ 0 };

	bool rightPressed{ false };
	bool leftPressed{ false };
	int arcball_on = false;

	glm::vec3 target{0.0f, 0.0f, 0.f};

	void scroll_setup(GLFWwindow* window);

private:

	glm::vec3 up;
	glm::vec3 direction;
	glm::vec3 right;

	void to_viewport(double& xpos, double& ypos);

	void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

	bool trackBallCamera::in_viewport();
};

glm::mat4 quaternion_2_matrix(const Quaternion q);

#endif