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

constexpr int DIM = 3;
const int NUM_FILES = 100;
int WINDOW_HEIGHT  = 700;
int WINDOW_WIDTH = 700;

void updateParticles(int index, std::vector<Particle>& particles) {
    int n;
    int steps;
    double dt;
    std::vector<double> masses;
    std::vector<Vec<DIM>> positions;
    std::vector<Vec<DIM>> velocities;
    std::vector<Vec<DIM>> forces;
    utils::readFromFile("../build/test1-MPI."+std::to_string(index)+".out", n, steps, dt, masses, positions, velocities);
    

    for (int i = 0; i < n; i++) {;
        particles[i].position = glm::vec3(positions[i][0], positions[i][1], DIM==2 ? 0.0f : positions[i][2]);
        if (particles[i].trace.size() < NUM_FILES-1)
            particles[i].trace.push_back(particles[i].position);
    }
}

void drawAxes3D(Shader& shader, float length = 1000.0f) {
    glm::mat4 model = glm::mat4(1.0f);
    shader.setMat4("model", model);

    // lines in X, Y, Z directions
    float verts[] = {
        -length, 0.0f, 0.0f,
         length, 0.0f, 0.0f,

         0.0f, -length, 0.0f,
         0.0f,  length, 0.0f,

         0.0f, 0.0f, -length,
         0.0f, 0.0f,  length
    };

    GLuint VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(verts), verts, GL_STATIC_DRAW);

    glVertexAttribPointer(
        0, 3, GL_FLOAT, GL_FALSE,
        3 * sizeof(float), (void*)0
    );
    glEnableVertexAttribArray(0);

    shader.setVec3("color", glm::vec3(1, 0, 0));
    glDrawArrays(GL_LINES, 0, 2);


    shader.setVec3("color", glm::vec3(0, 1, 0));
    glDrawArrays(GL_LINES, 2, 2);

    shader.setVec3("color", glm::vec3(0, 0, 1));
    glDrawArrays(GL_LINES, 4, 2);

    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
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
    float camRadius = 400.0f;
    float camSpeed  = 0.1f;
    

    // particle setup
    int currentParticleSet = 1; // first one already loaded
    int n;
    int steps;
    double dt;
    std::vector<Particle> particles;
    std::vector<double> masses;
    std::vector<Vec<DIM>> positions;
    std::vector<Vec<DIM>> velocities;
    std::vector<Vec<DIM>> forces;
    
    utils::readFromFile("../build/test1-MPI.0.out", n, steps, dt, masses, positions, velocities);
    particles.reserve(n);

    // initialize all particles
    for (int i = 0; i < n; i++) {
        particles.emplace_back(Particle(
            glm::vec3(positions[i][0], positions[i][1], DIM==2 ? 0.0f : positions[i][2]), 
            glm::vec3(i*1.0f/n, 0.0f, 1.0f), 
            masses[i]/20.0f));
        particles[i].trace.push_back(particles[i].position);
        particles[i].init();
    }
    
    // control simulation speed when playing it
    double prevFrameTime = glfwGetTime();
    bool playSim = false;
    // main render loop
    while (!glfwWindowShouldClose(window)) {
        double time = glfwGetTime();
        if (playSim && time - prevFrameTime > dt*10) { // slowdown factor of 10
            prevFrameTime = time;
            currentParticleSet = (currentParticleSet + 1) % NUM_FILES;
            updateParticles(currentParticleSet, particles);
        }

        glClearColor(0.05f, 0.05f, 0.08f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        shader.use();

        // rotating camera setup
        camPos = glm::vec3(
            camRadius * std::cos(camSpeed * time), 
            camRadius * 0.2f,
            camRadius * std::sin(camSpeed * time));
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

        // draw the 3d axis for reference
        drawAxes3D(shader);

        // IMGUI scene
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin("Particle Loader");
        if (ImGui::Button("Next Timestep")) {
            currentParticleSet = (currentParticleSet + 1) % NUM_FILES;
            updateParticles(currentParticleSet, particles);
        }

        ImGui::Checkbox("Play Simulation", &playSim);
        ImGui::End();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
