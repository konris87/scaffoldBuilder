#include "OpenGlSetup/defaultCamera.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <imgui_impl_glfw.h>

defaultCamera::defaultCamera(
    GLFWwindow* window, double split,
    glm::vec3 cameraPos,
    glm::vec3 cameraTarget,
    float speed) : split(split), window(window), position(cameraPos), target(cameraTarget), speed(speed) {

    FoV = 45.0f;

    direction = glm::normalize(position - target);

    //direction = (position - target).normalize();

    //right = Vec3(0.0f, 1.0f, 0.0f).cross(direction).normalize();

    right = glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), direction));

    up = glm::normalize(glm::cross(direction, right));
    //up = direction.cross(right).normalize();

    horizontalAngle = 3.14f;
    verticalAngle = 0.0f;
    mouseSpeed = 0.001f;
    fovSpeed = 2.0f;

    glfwGetWindowSize(window, &width, &height);

    // Automatically set up the scroll callback
    //scroll_setup(window);
};

void defaultCamera::update() {

    // Get mouse position

    glm::vec3 currPos = position;

    if (in_viewport()) {

        // glfwGetTime is called only once, the first time this function is called
        static double lastTime = glfwGetTime();

        // Compute time difference between current and last frame
        double currentTime = glfwGetTime();
        float deltaTime = float(currentTime - lastTime);
        direction = glm::normalize(position - target);

        //direction = (position - target).normalize();

        //right = Vec3(0.0f, 1.0f, 0.0f).cross(direction).normalize();

        right = glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), direction));

        up = glm::normalize(glm::cross(direction, right));

        // Move forward: If 'W' is pressed, move the camera forward along the direction vector
        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
            position = position - (speed * deltaTime * direction);
        }

        // Move backward: If 'S' is pressed, move the camera backward along the direction vector
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
            position = position + (speed * deltaTime * direction);
        }

        // Strafe right: If 'D' is pressed, move the camera to the right along the right vector
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
            position = position + (speed * deltaTime * right);
        }

        // Strafe left: If 'A' is pressed, move the camera to the left along the right vector
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
            position = position - (speed * deltaTime * right);
            //position = position - (speed * deltaTime * aglm::vec3(1.0, 0.0, 0.0));
        }

        //// Zoom in if the up arrow key is pressed
        //if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) {
        //    FoV -= fovSpeed * deltaTime; // Narrow the FoV for zoom in effect
        //}

        //// Zoom out if down arrow key is pressed
        //if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) {
        //    FoV += fovSpeed * deltaTime; // Widen the FoV for zoom out effect
        //}

        // For the next frame, the "last time" will be "now"
        lastTime = currentTime;
    }
    else {
        position = currPos;
    }

    projectionMatrix = glm::perspective(glm::radians(FoV), (float)width / (float)height, 0.1f, 100.0f);

    //Vec3 r = position 
    // Update the view matrix to reflect the camera's current position and orientation
    viewMatrix = glm::lookAt(
        glm::vec3(position.x, position.y, position.z),
        glm::vec3(target.x, target.y, target.z),
        glm::vec3(up.x, up.y, up.z)
    );
};

void defaultCamera::to_viewport(double& xpos, double& ypos) {

    xpos += split * width;
    
};

void defaultCamera::set_window_size(int newWidth, int newHeight) {
    
    width = newWidth;
    height = newHeight;

};

void defaultCamera::scroll_setup(GLFWwindow* window) {
    glfwSetWindowUserPointer(window, this);
    glfwSetScrollCallback(window, [](GLFWwindow* window, double xoffset, double yoffset) {
        defaultCamera* camera = (defaultCamera*)glfwGetWindowUserPointer(window);
        camera->scroll_callback(window, xoffset, yoffset);
    });
};

void defaultCamera::scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    ImGui_ImplGlfw_ScrollCallback(window, xoffset, yoffset);

    if (in_viewport()) {
        this->FoV -= (float)yoffset;
        if (this->FoV < 1.0f)
            this->FoV = 1.0f;
        if (this->FoV > 45.0f)
            this->FoV = 45.0f;
    }
}

bool defaultCamera::in_viewport() {

    double xPos, yPos;
    glfwGetCursorPos(window, &xPos, &yPos);

    int width, height;
    glfwGetWindowSize(window, &width, &height);

    return (xPos >= split * width && xPos <= width &&
        yPos >= 0 && yPos <= height);
};