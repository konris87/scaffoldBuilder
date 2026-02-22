#ifndef TRACKBALLCAMERA_H
#define TRACKBALLCAMERA_H
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <imgui.h>

// custom headers
#include <Math/Quaternion.h>
#include "Math/Vec.h"

enum ProjectionMode {
	Ortho, Perspective
};

class TrackBall {

public:
	TrackBall();
	TrackBall(GLFWwindow* window, float radius, int viewportWidth, int viewportHeight, int viewportX, int viewportY);
	~TrackBall();

	// Functions
	void update();
	void reset();

	// Setters
	void set(float r, int w, int h, int originX);
	void set_radius(float r);
	void set_viewport(int viewX, int viewY, int viewW, int viewH);
	void set_screen_size(int w, int h, int originX);
	void set_position(float x, float y, float z);
	void set_target(float x, float y, float z);
	void set_projection_plane(float nearPlane, float farPlane);
	void set_view(glm::vec3 newDir, glm::vec3 newUp);
	void set_projection_ortho();
	void set_projection_perspective();
	void set_projection_mode(ProjectionMode newMode);

	// Getters
	int get_screen_width() const;
	int get_screen_height() const;
	float get_radius() const;
	glm::mat4 get_view_matrix() const;
	glm::mat4 get_projection_matrix() const;
	glm::mat4 TrackBall::get_rotation_matrix() const;
	glm::mat4 TrackBall::get_translation_matrix() const;
	Vec3 get_point_on_sphere();
	Vec3 get_position() const;
	Vec3 get_target() const;

private:

	GLFWwindow* window;

	// viewport info
	int viewportX = 0;
	int viewportY = 0;
	int viewportW = 0;
	int viewportH = 0;

	float radius;
	float FoV{ 60.0f };
	float near{ 0.01f };
	float far{ 100.0f };

	glm::mat4 viewMatrix = glm::mat4(1.0f);
	glm::mat4 projectionMatrix = glm::mat4(1.0f);

	glm::mat4 rotMatrix = glm::mat4(1.0f);
	glm::mat4 translationMatrix = glm::mat4(1.0f);

	glm::vec3 direction = glm::vec3(0.0f, 0.0f, -1.0f);
	glm::vec3 position{ 0.0f, 0.0f, 10.0f };
	glm::vec3 target{ 0.0f, 0.0f, 0.0f };
	glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);
	glm::vec3 right = glm::vec3(1.0f, 0.0f, 0.0f);

	float speed{ 100.0f };

	bool isMiddleDragging = false;
	bool isRightDragging = false;
	bool isLeftDragging = false;
	double lastMouseX = 0.0, lastMouseY = 0.0;

	double lastZoomY = 0.0;

	float mouseRotX = 0.0, mouseRotY = 0.0;
	double lastRotX = 0.0, lastRotY = 0.0;

	// these are rotation parameters
	float halfViewWidth{ 0.0f };
	float halfViewHeight{ 0.0f };

	Quaternion quat;
	Quaternion prevQuat;

	ProjectionMode mode = ProjectionMode::Ortho;

	// ----------------------------------------------
	// Functions
	Vec3 TrackBall::get_vector(float x, float y);

	Vec3 TrackBall::get_unit_vector(float x, float y);

	Vec3 sphereVector{ 0.0f, 0.0f, 0.0f };

	glm::mat4 mT{ 1.0f };

	// window split
	float split = 600.0f;

	glm::vec3 translationVec{ 0.0f, 0.0f, 0.0f };

	glm::vec3 initialOffset;

	glm::vec3 initialDir;

	glm::vec3 initialUp = glm::vec3(0.f, 1.0f, 0.f);

	float distance;

	bool is_inside_viewport(double mouseX, double mouseY) {
		
		return mouseX >= viewportX && mouseX < (viewportX + viewportW) &&
			mouseY >= viewportY && mouseY < (viewportY + viewportH);
	};

	void map_to_viewport(double xPos, double yPos, float& localX, float& localY) const;
};

#endif