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
#include "Camera.hpp"

constexpr int DIM = 3;
const int NUM_FILES = 100;
constexpr int WINDOW_HEIGHT = 600;
constexpr int WINDOW_WIDTH = 900;

void processKeyboardInputs(GLFWwindow* window, float deltaTime, Camera& camera) {
    float yawDiff = 0, pitchDiff = 0;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        yawDiff += deltaTime;
    }

    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        yawDiff -= deltaTime;
    }

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        pitchDiff += deltaTime;
    }

    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        pitchDiff -= deltaTime;
    }

    camera.rotateCamera(pitchDiff, yawDiff);
}

void updateParticles(int index, std::vector<Particle>& particles, std::vector<glm::vec3>& combinedTraces) {
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
        if (combinedTraces.size() < n*(NUM_FILES-1)) {
            combinedTraces.push_back(particles[i].position);
            if (combinedTraces.size() > 10*n) {
                combinedTraces.erase(combinedTraces.begin(), combinedTraces.begin()+n);
            }
        }
    }
}


void createSphere(std::vector<float>& vertices,
                  std::vector<unsigned int>& indices,
                  unsigned int X_SEGMENTS,
                  unsigned int Y_SEGMENTS)
{
    vertices.clear();
    indices.clear();

    for (unsigned int y = 0; y <= Y_SEGMENTS; ++y) {
        for (unsigned int x = 0; x <= X_SEGMENTS; ++x) {
            float xSegment = (float)x / (float)X_SEGMENTS;
            float ySegment = (float)y / (float)Y_SEGMENTS;
            float xPos = std::cos(xSegment * 2.0f * M_PI) * std::sin(ySegment * M_PI);
            float yPos = std::cos(ySegment * M_PI);
            float zPos = std::sin(xSegment * 2.0f * M_PI) * std::sin(ySegment * M_PI);

            vertices.push_back(xPos);
            vertices.push_back(yPos);
            vertices.push_back(zPos);
        }
    }

    for (unsigned int y = 0; y < Y_SEGMENTS; ++y) {
        for (unsigned int x = 0; x < X_SEGMENTS; ++x) {
            unsigned int i0 = y * (X_SEGMENTS + 1) + x;
            unsigned int i1 = (y + 1) * (X_SEGMENTS + 1) + x;

            indices.push_back(i0);
            indices.push_back(i1);
            indices.push_back(i0 + 1);

            indices.push_back(i0 + 1);
            indices.push_back(i1);
            indices.push_back(i1 + 1);
        }
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
    GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "N-Body Simulation", nullptr, nullptr);
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
    Camera camera(400.0f, shader, WINDOW_WIDTH, WINDOW_HEIGHT);
    

    // particle setup
    int currentParticleSet = 1; // first one already loaded
    int n;
    int steps;
    double dt;
    std::vector<Particle> particles;
    std::vector<glm::vec3> combinedTraces; // all combined to improve performance
    std::vector<double> masses;
    std::vector<Vec<DIM>> positions;
    std::vector<Vec<DIM>> velocities;
    std::vector<Vec<DIM>> forces;
    
    utils::readFromFile("../build/test1-MPI.0.out", n, steps, dt, masses, positions, velocities);
    particles.reserve(n);

    // Set up render objects for particles
    GLuint VAO, VBO, EBO;
    GLsizei indexCount;
    std::vector<float> sphereVerts;
    std::vector<unsigned int> sphereInds;
    createSphere(sphereVerts, sphereInds, 8, 8);
    indexCount = (GLsizei)sphereInds.size();

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER,
                 sphereVerts.size() * sizeof(float),
                 sphereVerts.data(),
                 GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 sphereInds.size() * sizeof(unsigned int),
                 sphereInds.data(),
                 GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);

    // Set up render objects for traces
    GLuint traceVAO, traceVBO;
    glGenVertexArrays(1, &traceVAO);
    glGenBuffers(1, &traceVBO);

    glBindVertexArray(traceVAO);
    glBindBuffer(GL_ARRAY_BUFFER, traceVBO);

    glBufferData(
        GL_ARRAY_BUFFER,
        NUM_FILES * n * sizeof(glm::vec3),
        nullptr,
        GL_DYNAMIC_DRAW
    );

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), 0);
    glEnableVertexAttribArray(0);

    // initialize all particles
    for (int i = 0; i < n; i++) {
        particles.emplace_back(Particle(
            glm::vec3(positions[i][0], positions[i][1], DIM==2 ? 0.0f : positions[i][2]), 
            glm::vec3(i*1.0f/n, 0.0f, 1.0f), 
            masses[i]/20.0f));
        particles[i].VAO = VAO;
        particles[i].VBO = VBO;
        particles[i].EBO = EBO;
        particles[i].indexCount = indexCount;
    }
    
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    
    // control simulation speed when playing it
    double prevSimFrameTime = glfwGetTime(); // last time we update simulation timestep
    double prevLoopTime = glfwGetTime(); // last time opengl drew onto screen (mainloop)
    bool playSim = false;
    bool disableTrace = false;
    float scale = 1.0;
    float zoom = 1.0;
    // main render loop
    while (!glfwWindowShouldClose(window)) {
        double time = glfwGetTime();
        processKeyboardInputs(window,time-prevLoopTime, camera);

        if (playSim && time - prevSimFrameTime > dt) { // slowdown factor of 10
            prevSimFrameTime = time;
            currentParticleSet = (currentParticleSet + 1) % NUM_FILES;
            updateParticles(currentParticleSet, particles, combinedTraces);
        }
        
        glClearColor(0.05f, 0.05f, 0.08f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        shader.use();
        
        
        // draw particles and trace
        for (Particle& p : particles) {
            p.draw(shader, scale);
        }
        
        if (!disableTrace) {
            glm::mat4 model = glm::mat4(1.0f);
            shader.setMat4("model", model);
            shader.setVec3("color", glm::vec3(1.0f, 1.0f, 1.0f));
            int traceCount = combinedTraces.size();
            glBindBuffer(GL_ARRAY_BUFFER, traceVBO);
            glBufferSubData(
                GL_ARRAY_BUFFER,
                0,
                traceCount * sizeof(glm::vec3),
                combinedTraces.data()
            );
            glBindVertexArray(traceVAO);
            glDrawArrays(GL_POINTS, 0, combinedTraces.size());
        }
        
        
        
        // IMGUI scene
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        
        ImGui::Begin("Control Panel");
        ImGui::Spacing();
        ImGui::Text("Time Settings");
        if (ImGui::Button("Next Timestep")) {
            currentParticleSet = (currentParticleSet + 1) % NUM_FILES;
            updateParticles(currentParticleSet, particles, combinedTraces);
        }
        
        ImGui::Checkbox("Loop Simulation", &playSim);
        
        ImGui::Spacing();
        ImGui::Text("Visual Settings");
        ImGui::Checkbox("Disable Trace", &disableTrace);
        ImGui::SliderFloat("Particle Scale", &scale, 0.1f, 2.0f, "ratio = %.1f");
        
        ImGui::Spacing();
        ImGui::Text("Movement");
        if (ImGui::SliderFloat("Zoom", &zoom, 0.5f, 10.0f, "ratio = %.1f")) {
            camera.zoomCamera(zoom);
        }
        
        
        ImGui::End();
        
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        prevLoopTime = time;
        
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    glfwTerminate();
    return 0;
}
