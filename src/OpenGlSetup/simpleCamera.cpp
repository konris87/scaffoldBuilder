#include "OpenGlSetup/simpleCamera.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>

SimpleCamera::SimpleCamera() {};

SimpleCamera::SimpleCamera(int vW, int vH) :

	viewportWidth(vW), viewportHeight(vH) {

	direction = glm::normalize(position - target);

	right = glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), direction));

	up = glm::normalize(glm::cross(direction, right));

	projectionMatrix = glm::perspective(glm::radians(FoV), (float)viewportWidth / viewportHeight, near, far);

	viewMatrix = glm::lookAt(position, target, up);

};

glm::mat4 SimpleCamera::get_view_matrix() const {

	glm::mat4 view;

	view = viewMatrix;

	return view;

};

glm::mat4 SimpleCamera::get_projection_matrix() const {

	glm::mat4 proj;

	proj = projectionMatrix;

	return proj;

};

void SimpleCamera::set_dimensions(const int vW, const int vH) {
	viewportWidth = vW;
	viewportHeight = vH;

	// update projection matrix
	projectionMatrix = glm::perspective(glm::radians(FoV), (float)viewportWidth / viewportHeight, near, far);
};