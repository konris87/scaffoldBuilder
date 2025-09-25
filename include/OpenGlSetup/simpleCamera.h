#ifndef CAMERA_H
#define CAMERA_H
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <Math/Quaternion.h>
#include <Math/Vec.h>

class SimpleCamera {
public:
	SimpleCamera();
	~SimpleCamera() {};
	SimpleCamera(int vW, int vH);

	glm::mat4 get_view_matrix() const;

	glm::mat4 get_projection_matrix() const;

	void set_dimensions(const int vW, const int vH);

private:

	int viewportWidth{ 800 };
	int viewportHeight{ 600 };

	float FoV{ 45.0f };
	float near{ 0.01f };
	float far{ 10.0f };

	glm::mat4 viewMatrix = glm::mat4(1.0f);
	glm::mat4 projectionMatrix = glm::mat4(1.0f);

	glm::vec3 direction = glm::vec3(0.0f, 0.0f, -1.0f);
	glm::vec3 position{ 0.0f, 0.0f, 3.0f };
	glm::vec3 target{ 0.0f, 0.0f, 0.0f };
	glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);
	glm::vec3 right = glm::vec3(1.0f, 0.0f, 0.0f);

};

#endif