#include "OpenGlSetup/Camera.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <Misc/Quaternion.h>

Camera::Camera(
    GLFWwindow* glWindow, const double split, const Vec3 targetCamera) : window(glWindow), split(split), targetCamera(targetCamera) {

    position = glm::vec3(0.0, 0.0, 20.0);
    horizontalAngle = 3.14f;
    verticalAngle = 0.0f;
    FoV = 45.0f;
    speed = 3.0f;
    mouseSpeed = 0.001f;
    fovSpeed = 2.0f;

    glfwGetWindowSize(window, &width, &height);

    //targetCamera.x = width * split + (1 - split) * width / 2;
    //targetCamera.x = width / 2.0f;
    //targetCamera.y = height / 2.0f;
    //targetCamera = Vec3(0.0, 0.0, 0.0);

    //xPrev = width * split + (1 - split) * width / 2;

    Quaternion q = Quaternion(0.0, 0.0, 0.0, 1.0);
    Quaternion qPrev;

    up = glm::vec3(0.0, 1.0, 0.0);
};

void Camera::update() {

    static double lastTime = glfwGetTime();

    // Compute time difference between current and last frame
    double currentTime = glfwGetTime();
    float deltaTime = float(currentTime - lastTime);

    // Get mouse position
    //glfwSetMouseButtonCallback(window, mouse_button_callback);

    double xPos, yPos;
    glfwGetCursorPos(window, &xPos, &yPos);

    xCurr = xPos;
    yCurr = yPos;

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {

        if (!leftPressed) {

            leftPressed = true;

            // get the current position for the next frame
            xPrev = xCurr;
            yPrev = yCurr;
            qPrev = q;
        }

        /*std::cout << "x Prev: " << xPrev << "y Prev: " << yPrev << std::endl;
        std::cout << "x Curr: " << xCurr << "y Curr: " << yCurr << std::endl;*/

        // Sphere projection for dragging (or just use dx, dy directly for simple rotation)
        /*Vec3 v1 = sphere_projection(xPrev, yPrev, targetCamera);
        Vec3 v2 = sphere_projection(xCurr, yCurr, targetCamera);*/

        Vec3 v1 = arc(xPrev, yPrev);
        Vec3 v2 = arc(xCurr, yCurr);

        // Compute quaternion rotation between v1 and v2
        Quaternion qd = get_from_v1_to_v2(v1, v2);

        // Update the current quaternion
        q = qd * qPrev;
        q.normalize();
    }
    else {
        leftPressed = false;
        qPrev = q;
        xPrev = xCurr;
        yPrev = yCurr;
    }
    
    // reset if press R
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
        q = Quaternion(0.0, 0.0, 0.0, 1.0);
    }

    // Task 5.7: construct projection and view matrices
    //projectionMatrix = glm::perspective(glm::radians(FoV), 4.0f / 3.0f, 0.1f, 100.0f);
    projectionMatrix = glm::perspective(glm::radians(FoV), (float)width/height, 0.1f, 1000.0f);
    //viewMatrix = lookAt(
    //    position,
    //    direction,
    //    up
    //);

    glm::mat4 mt(1);
    glm::mat4 mr(1);

    Mat4 mrot = q.to_matrix();
    
    mr[0] = glm::vec4(mrot.r00, mrot.r10, mrot.r20, 0.0f);
    mr[1] = glm::vec4(mrot.r01, mrot.r11, mrot.r21, 0.0f);
    mr[2] = glm::vec4(mrot.r02, mrot.r12, mrot.r22, 0.0f);

    Quaternion pos = Quaternion(position.x, position.y, position.z, 0.0);

    Quaternion rotPos = q.inverse() * pos * q;

    mt[3] = glm::vec4(targetCamera.x - rotPos.x, targetCamera.y - rotPos.y, targetCamera.z - rotPos.z, 1.0);

    //// new direction
    //direction = glm::normalize(- glm::vec3(rotPos.x, rotPos.y, rotPos.z));

    //// new right
    //right = glm::normalize(glm::cross(direction, glm::vec3(0.0f, 1.0f, 0.0f)));

    //// new up
    //up = glm::cross(right, direction);

    viewMatrix = mr * mt;

    lastTime = currentTime;
};

void Camera::update_scroll(float yoffset) {
    FoV -= (float)yoffset;
    if (FoV < 1.0f)
        FoV = 1.0f;
    if (FoV > 80.0f)
        FoV = 80.0f;
};

Vec3 Camera::sphere_projection(double xpos, double ypos, Vec3 target) {

    // https://www.khronos.org/opengl/wiki/Object_Mouse_Trackball

    float mx = xpos - width / 2.0f;
    float my = height / 2.0f - ypos;
    //float my = ypos - height/2.0f;

    // this initial vector in screen coordinates
    Vec3 pj = Vec3(mx, my, 0.0f);

    float lenSquared = pj.x * pj.x + pj.y * pj.y;
    float radSquared = sphereRadius * sphereRadius;

    if (lenSquared <= (0.5f * radSquared)) {
        pj.z = sqrtf(radSquared - lenSquared);
    }
    else {
        pj.z = 0.5f * radSquared / sqrtf(lenSquared);

        // now scale x, y so they fit to the sphere. Scaling is performed
        // based on the aspect ratio between these two aspectRatio = y / x
        float newX, newY, aspectRatio;
        if (mx == 0.0) {

            // use the sphere equation to get the new y
            newX = 0.0;
            newY = sqrtf(radSquared - pj.z * pj.z);

            if (newY < 0) { newY *= -1; }
        }
        else {
            aspectRatio = my / mx;
            newX = sqrtf((radSquared - pj.z * pj.z) / (1 + aspectRatio * aspectRatio));
            if (newX < 0.0) {
                newX *= -1;
            }
            newY = aspectRatio * newX;
        }
        pj.x = newX;
        pj.y = newY;
    }

    // normalize p
    pj.normalize();
       
    return pj;
    
};

Vec3 Camera::arc(double xpos, double ypos) {

    xpos -= width * split;

    //float mx = xpos - width / 2.0f;
    float mx = xpos - (1 - split) * width / 2.0f ;
    float my = height / 2.0f - ypos;

    float arc = sqrtf(mx * mx + my * my);
    float a = arc / sphereRadius;
    float b = atan2f(my, mx);
    float newX = sphereRadius * sinf(a);

    Vec3 v;
    v.x = newX * cosf(b);
    v.y = newX * sinf(b);
    v.z = sphereRadius * cosf(a);

    v.normalize();

    return v;
}

void Camera::to_viewport(double& xpos, double& ypos) {

    xpos += split * width;
    
};

void Camera::set_window_size(int newWidth, int newHeight) {
    
    width = newWidth;
    height = newHeight;

};