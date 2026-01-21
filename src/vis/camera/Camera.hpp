#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "../shaders/shader.hpp"

class Camera {
    public:
        Camera(float radius, Shader& shader, int windowWidth, int windowHeight);

        void rotateCamera(float pitch, float yaw);
        void zoomCamera(float offset);
        
        glm::mat4 view;
        glm::mat4 projection;
    private:
        void updateCamera();

        Shader shader;
        float radius;
        float pitch;
        float yaw;
        float zoom;


        int windowWidth;
        int windowHeight;

        float pos;
};