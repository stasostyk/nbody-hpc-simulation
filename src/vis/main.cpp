#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "shaders/shader.hpp"
#include "particle/Particle.hpp"

#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include "../utils.hpp"

constexpr int DIM = 2;
const int NUM_FILES = 1000;
int WINDOW_HEIGHT  = 700;
int WINDOW_WIDTH = 700;

void loadParticles(int index, std::vector<Particle>& particles) {
    particles.clear();
    int steps;
    double dt;
    bodies<DIM> bodies;
    utils::readFromFile("../build/test-MPI-reduced-"+std::to_string(index)+".out", steps, dt, bodies);

    for (int i = 0; i < bodies.globalSize(); i++) {
        Particle p(glm::vec3(bodies.position[i][0], bodies.position[i][1], 0.0f), glm::vec3(i*1.0f/bodies.globalSize(), 0.0f, 1.0f), bodies.mass[i]/20.0f);
        p.init();
        particles.push_back(p);
    }
}

int main() {
    // GLFW setup
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return -1;
    }
    // set version to 3.3 (needed for WSL support)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    
    
    // window setup
    GLFWwindow* window = glfwCreateWindow(WINDOW_HEIGHT, WINDOW_WIDTH, "N-Body Simulation", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);


    // GLAD setup
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Failed to initialize GLAD\n";
        glfwTerminate();
        return -1;
    }


    // IMGUI setup
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
    

    // scene setup
    Shader shader("./shaders/basic.vert", "./shaders/basic.frag");
    glm::vec3 camPos    = glm::vec3(0.0f, 0.0f, 400.0f);
    glm::vec3 camTarget = glm::vec3(0.0f, 0.0f, 0.0f);
    glm::vec3 camUp     = glm::vec3(0.0f, 1.0f, 0.0f);
    

    // particle setup
    int currentParticleSet = 0;
    std::vector<Particle> particles;
    loadParticles(currentParticleSet, particles);
    
    
    // main render loop
    while (!glfwWindowShouldClose(window)) {
        float time = static_cast<float>(glfwGetTime());

        glClearColor(0.05f, 0.05f, 0.08f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        shader.use();

        // camera setup
        glm::mat4 view = glm::lookAt(camPos, camTarget, camUp);
        float aspect = (WINDOW_HEIGHT == 0) ? 1.0f : (float)WINDOW_WIDTH / (float)WINDOW_HEIGHT;
        glm::mat4 projection = glm::perspective(glm::radians(45.0f), aspect, 0.1f, 1000.0f);
        shader.setMat4("view", view);
        shader.setMat4("projection", projection);

        // draw particles and trace
        for (Particle& p : particles) {
            p.drawTrace(shader);
            p.draw(shader);
        }

        // IMGUI scene
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin("Particle Loader");
        if (ImGui::Button("Next Timestep")) {
            currentParticleSet = (currentParticleSet + 20) % NUM_FILES;
            loadParticles(currentParticleSet, particles);
        }
        ImGui::End();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
