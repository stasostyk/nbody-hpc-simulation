#include "Camera.hpp"

Camera::Camera(float radius, Shader& shader, int windowWidth, int windowHeight): 
    radius(radius), 
    shader(shader), 
    windowWidth(windowWidth),
    windowHeight(windowHeight),
    pitch(0.0f),
    yaw(0.0f),
    zoom(1.0f)
    {
}

void Camera::rotateCamera(float pitchDiff, float yawDiff) {
    pitch += pitchDiff;
    pitch = glm::clamp(pitch, -glm::radians(89.0f), glm::radians(89.0f));
    yaw += yawDiff;
    updateCamera();
}

void Camera::zoomCamera(float newZoom) {
    zoom = newZoom;
    updateCamera();
}

void Camera::updateCamera() {
    glm::vec3 camPos = glm::vec3(
        radius * std::cos(pitch) * std::cos(yaw), 
        radius * std::sin(pitch), 
        radius * std::cos(pitch) * std::sin(yaw)) / zoom;
    
    glm::vec3 camTarget = glm::vec3(0.0f, 0.0f, 0.0f);
    glm::vec3 camUp = glm::vec3(0.0f, 1.0f, 0.0f);

    glm::mat4 view = glm::lookAt(camPos, camTarget, camUp);
    float aspect = (windowHeight == 0) ? 1.0f : (float)windowWidth / (float)windowHeight;
    glm::mat4 projection = glm::perspective(glm::radians(45.0f), aspect, 0.1f, 2000.0f);
    shader.setMat4("view", view);
    shader.setMat4("projection", projection);
}