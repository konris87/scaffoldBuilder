#include "OpenGLSetup/trackBallCamera.h"
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>

const float EPSILON = 0.00001f;

TrackBall::TrackBall() : radius(2.0f), viewportX(0), viewportY(0), viewportW(800), viewportH(600) {

    halfViewWidth = viewportW / 2.0f;
    halfViewHeight = viewportH / 2.0f;

    initialOffset = position - target;

    initialDir = glm::normalize(position - target);

    set_projection_ortho();
};

TrackBall::~TrackBall() {};

TrackBall::TrackBall(GLFWwindow* window, float radius, int viewportWidth, int viewportHeight, int viewportX, int viewportY) :
    window(window),
    viewportW(viewportWidth),
    viewportH(viewportHeight),
    viewportX(viewportX),
    viewportY(viewportY),
    radius(radius) {

    halfViewWidth = viewportWidth / 2.0f;
    halfViewHeight = viewportHeight / 2.0f;

    initialOffset = position - target;
    initialDir = glm::normalize(position - target);

    direction = glm::normalize(position - target);

    distance = glm::length(position - target);

    right = glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), direction));

    up = glm::normalize(glm::cross(direction, right));

    rotMatrix[0][0] = right.x; // First column, first row
    rotMatrix[1][0] = right.y; // First column, first row
    rotMatrix[2][0] = right.z;
    rotMatrix[0][1] = up.x; // First column, second row
    rotMatrix[1][1] = up.y;
    rotMatrix[2][1] = up.z;
    rotMatrix[0][2] = direction.x; // First column, third row
    rotMatrix[1][2] = direction.y;
    rotMatrix[2][2] = direction.z;
};

void TrackBall::update() {

    static double lastTime = glfwGetTime();

    // Compute time difference between current and last frame
    double currentTime = glfwGetTime();
    float deltaTime = float(currentTime - lastTime);

    double xPos, yPos;
    glfwGetCursorPos(window, &xPos, &yPos);

    // update camera only if the cursor and the action is captured inside the viewport
    if (is_inside_viewport(xPos, yPos)) {

        // 1. Translation Handling
        // check if middle button is dragged
        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS) {
            if (!isMiddleDragging) {
                isMiddleDragging = true;
                lastMouseX = xPos;
                lastMouseY = yPos;
            }
            else {
                float viewportScale = distance * 0.002f;
                float dx = float(xPos - lastMouseX);
                float dy = float(yPos - lastMouseY);
                glm::vec3 panRight = -dx * viewportScale * right;
                glm::vec3 panUp = dy * viewportScale * up;

                position += panRight + panUp;
                target += panRight + panUp;

                lastMouseX = xPos;
                lastMouseY = yPos;
            }
        }
        else {
            isMiddleDragging = false;
        }

        // 2. Rotation Handling is controlled by left mouse button dragging 
        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {

            // if the left mouse was not dragged previously
            if (!isLeftDragging) {
                isLeftDragging = true;
                lastRotX = xPos;
                lastRotY = yPos;
                prevQuat = quat;
            }
            // else it was dragged
            else {

                // create the quaternion from current to next position
                Vec3 v1 = get_unit_vector(lastRotX, lastRotY);

                Vec3 v2 = get_unit_vector(xPos, yPos);

                Quaternion delta = get_from_v1_to_v2(v1, v2);

                quat = delta * prevQuat;
            }
        }
        // if releashed change mode
        else {
            isLeftDragging = false;
            lastRotX = xPos;
            lastRotY = yPos;
        }

        // 3. Zoom handling right mouse
        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
            if (!isRightDragging) {
                isRightDragging = true;
                lastZoomY = yPos;
            }
            else {
                float dy = float(yPos - lastZoomY);
                FoV += dy * speed * deltaTime;

                // clamp
                if (FoV < 0.1f)
                    FoV = 0.1f;
                if (FoV > 120.0f)
                    FoV = 120.0f;
                lastZoomY = yPos;
            }
        }
        else {
            isRightDragging = false;
        }

    }

    rotMatrix = quat.to_matrix().to_glm();

    direction = target - position;
    //double distance = glm::length(direction) * 1.5f;
    distance = glm::length(direction) * 1.5f;
    direction = glm::normalize(direction);

    right = glm::vec3(rotMatrix[0][0], rotMatrix[1][0], rotMatrix[2][0]);

    up = glm::vec3(rotMatrix[0][1], rotMatrix[1][1], rotMatrix[2][1]);

    if (mode == ProjectionMode::Ortho) {
        set_projection_ortho();
    }
    else {
        set_projection_perspective();
    }
    viewMatrix = glm::translate(glm::mat4(1.0), glm::vec3(0.0, 0.0, -distance)) * rotMatrix * glm::translate(glm::mat4(1.0), -target);

    lastTime = currentTime;

};

// Setters
void TrackBall::set(float r, int w, int h, int originX) {
    radius = r;
    viewportW = w;
    viewportH = h;
    halfViewWidth = w * 0.5f;
    halfViewHeight = h * 0.5f;
    viewportX = originX;
};

void TrackBall::set_radius(float r) {
    radius = r;
};

void TrackBall::set_screen_size(int w, int h, int originX) {
    viewportW = w;
    viewportH = h;
    halfViewWidth = w * 0.5f;
    halfViewHeight = h * 0.5f;
    viewportX = originX;
};

void TrackBall::set_position(float x, float y, float z) {
    position.x = x;
    position.y = y;
    position.z = z;
};

void TrackBall::set_target(float x, float y, float z) {
    target.x = x;
    target.y = y;
    target.z = z;
};

void TrackBall::set_viewport(int viewX, int viewY, int viewW, int viewH) {
    viewportX = viewX;
    viewportY = viewY;
    viewportW = viewW;
    viewportH = viewH;
};


void TrackBall::set_view(glm::vec3 newDir, glm::vec3 newUp) {

    direction = glm::normalize(newDir);

    right = glm::normalize(glm::cross(newUp, direction));

    up = glm::normalize(glm::cross(direction, right));

    right = glm::normalize(glm::cross(up, direction));

    glm::vec3 test = glm::cross(up, right);

    float distance = glm::length(target - position) * 1.1f;

    // Update rotation matrix from the new basis
    rotMatrix = glm::mat4(1.0f);
    rotMatrix[0][0] = right.x;
    rotMatrix[1][0] = right.y;
    rotMatrix[2][0] = right.z;
    rotMatrix[0][1] = up.x;
    rotMatrix[1][1] = up.y;
    rotMatrix[2][1] = up.z;
    rotMatrix[0][2] = direction.x;
    rotMatrix[1][2] = direction.y;
    rotMatrix[2][2] = direction.z;

    // Reset quaternion (optional, prevents conflict with arcball rotation)
    quat = Quaternion(rotMatrix);
    prevQuat = quat;
};

void TrackBall::set_projection_plane(float nearPlane, float farPlane) {
    near = nearPlane;
    far = farPlane;
};

void TrackBall::set_projection_ortho() {

    float orthoHeight = tan(glm::radians(FoV / 2.0f)) * distance;

    float aspect = (float)viewportW / (float)viewportH;

    projectionMatrix = glm::ortho(
        -orthoHeight * aspect, orthoHeight * aspect,
        -orthoHeight, orthoHeight,
        near, 1000.0f
    );
};

void TrackBall::set_projection_mode(ProjectionMode newMode) {
    mode = newMode;
};


void TrackBall::set_projection_perspective() {
    projectionMatrix = glm::perspective(glm::radians(FoV), (float)viewportW / viewportH, near, far);
};

// Getters
int TrackBall::get_screen_width() const { return viewportW; };

int TrackBall::get_screen_height() const { return viewportH; };

float TrackBall::get_radius() const { return radius; };

glm::mat4 TrackBall::get_rotation_matrix() const { return rotMatrix; };

glm::mat4 TrackBall::get_translation_matrix() const { return translationMatrix; };

glm::mat4 TrackBall::get_view_matrix() const {

    glm::mat4 view;

    view = viewMatrix;

    return view;

};

glm::mat4 TrackBall::get_projection_matrix() const {

    glm::mat4 proj;

    proj = projectionMatrix;

    return proj;

};

Vec3 TrackBall::get_unit_vector(float x, float y) {

    float localX, localY;
    map_to_viewport(x, y, localX, localY);

    Vec3 vec = get_vector(localX, localY);
    //Vec3 vec = get_vector(x, y);

    vec = vec.normalized();

    //std::cout << "vec normal: vx " << vec.x << " " << vec.y << " " << vec.z << std::endl;

    return vec;
};

Vec3 TrackBall::get_vector(float x, float y) {

    float nx = (2.0f * x / viewportW) - 1.0f;
    float ny = (2.0f * y / viewportH) - 1.0f;  // now y=-1 bottom, y=1 top
    
    float length2 = nx * nx + ny * ny;
    Vec3 vec;
    if (length2 <= 1.0f) {
        vec = Vec3(nx, ny, sqrtf(1.0f - length2));
    }
    else {
        float norm = 1.0f / sqrtf(length2);
        vec = Vec3(nx * norm, ny * norm, 0.0f);
    }

    //// center origin at viewport middle
    //float mx = x - viewportW * 0.5f;
    //float my = y - viewportH * 0.5f;

    //float arc = sqrtf(mx * mx + my * my);
    //float a = arc / radius;
    //float b = atan2f(my, mx);
    //float x2 = radius * sinf(a);

    //Vec3 vec;
    //vec.x = x2 * cosf(b);
    //vec.y = x2 * sinf(b);
    //vec.z = radius * cosf(a);

    return vec;
};

void TrackBall::map_to_viewport(double xPos, double yPos, float& localX, float& localY) const {
    // Convert from window coords to viewport-local coords
    // GLFW: origin is top-left, OpenGL viewport: origin is bottom-left
    localX = static_cast<float>(xPos - viewportX);
    localY = static_cast<float>((viewportH + viewportY) - yPos);

}

Vec3 TrackBall::get_point_on_sphere() {

    double xPos, yPos;
    glfwGetCursorPos(window, &xPos, &yPos);

    // Ensure the cursor is inside the active viewport
    if (!is_inside_viewport(xPos, yPos)){
        return Vec3(0, 0, 0);
    }

    float localX, localY;
    map_to_viewport(xPos, yPos, localX, localY);

    Vec3 vec = get_vector(xPos, yPos);

    return vec;

};

Vec3 TrackBall::get_position() const {

    glm::vec3 camOffset = glm::vec3(rotMatrix * glm::vec4(0, 0, distance, 1.0));
    glm::vec3 camPos = target + camOffset;

    Vec3 pos;

    //pos.x = camPos.x;
    //pos.y = camPos.y;
    //pos.z = camPos.z;

    pos.x = position.x;
    pos.y = position.y;
    pos.z = position.z;

    return pos;
};

Vec3 TrackBall::get_target() const {
    Vec3 vec(target.x, target.y, target.z);

    return vec;
};
