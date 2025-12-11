#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <vector>
#include <cmath>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "shaders/shader.hpp"
#include "particle/Particle.hpp"

int WINDOW_HEIGHT  = 700;
int WINDOW_WIDTH = 700;

int main() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return -1;
    }

    // set version to 3.3
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(WINDOW_HEIGHT, WINDOW_WIDTH, "N-Body Simulation", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Failed to initialize GLAD\n";
        glfwTerminate();
        return -1;
    }

    Shader shader("./shaders/basic.vert", "./shaders/basic.frag");

    std::vector<Particle> particles;
    for (int i = 0; i < 5; i++) {
        Particle p(glm::vec3(i*0.3f, 0.0f, 0.0f), glm::vec3(i*1.0f/5, 0.0f, 1.0f), 0.1);
        p.init();
        particles.push_back(p);
    }

    glm::vec3 camPos    = glm::vec3(3.0f, 3.0f, 3.0f);
    glm::vec3 camTarget = glm::vec3(0.0f, 0.5f, 0.0f);
    glm::vec3 camUp     = glm::vec3(0.0f, 1.0f, 0.0f);

    while (!glfwWindowShouldClose(window)) {
        float time = static_cast<float>(glfwGetTime());

        glClearColor(0.05f, 0.05f, 0.08f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        shader.use();

        glm::mat4 view = glm::lookAt(camPos, camTarget, camUp);
        float aspect = (WINDOW_HEIGHT == 0) ? 1.0f : (float)WINDOW_WIDTH / (float)WINDOW_HEIGHT;
        glm::mat4 projection = glm::perspective(glm::radians(45.0f), aspect, 0.1f, 100.0f);
        shader.setMat4("view", view);
        shader.setMat4("projection", projection);

        for (Particle& p : particles) {
            p.position = p.position + time*glm::vec3(0.0f,0.0f,0.0001f);
            p.drawTrace(shader);
            p.draw(shader);
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
