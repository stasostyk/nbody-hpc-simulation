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

#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include "../utils.hpp"
#include "renderBufferObjects/ParticleBufferObject.hpp"
#include "renderBufferObjects/TraceBufferObject.hpp"
#include "camera/Camera.hpp"

constexpr int DIM = 3;
const int NUM_FILES = 49;
constexpr int WINDOW_HEIGHT = 600;
constexpr int WINDOW_WIDTH = 900;
std::string FILE_PATH = "../build/test-MPI-reduced-";
// std::string FILE_PATH = "../build/test1-MPI.";


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

glm::vec3 particlePositionVec(Vec<DIM>& position) {
    return glm::vec3(position[0], position[1], DIM==3?position[2]:0.0f);
}

void updateParticles(int index, Bodies<DIM>& bodies, std::vector<glm::vec3>& combinedTraces) {
    double dt;
    int steps;
    utils::readFromFile(FILE_PATH+std::to_string(index)+".out", steps, dt, bodies, false);
    int n = bodies.position.size();

    for (int i = 0; i < n; i++) {;
        if (combinedTraces.size() < n*(NUM_FILES-1)) {
            combinedTraces.push_back(particlePositionVec(bodies.position[i]));
            if (combinedTraces.size() > 10*n) {
                combinedTraces.erase(combinedTraces.begin(), combinedTraces.begin()+n);
            }
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
    Shader basicShader("./shaders/basic.vert", "./shaders/basic.frag");
    Shader particleShader("./shaders/particle.vert", "./shaders/particle.frag");
    Camera camera(400.0f, basicShader, WINDOW_WIDTH, WINDOW_HEIGHT);
    

    // particle setup
    int currentParticleSet = 1; // first one already loaded
    double dt;
    int steps;
    std::vector<glm::vec3> combinedTraces; // all combined to improve performance
    Bodies<DIM> bodies;
    
    utils::readFromFile(FILE_PATH + "0.out", steps, dt, bodies, false);
    int n = bodies.position.size();
    // utils::readFromFile("../build/test-MPI-reduced-0.out", n, steps, dt, masses, positions, velocities);

    // Set up render objects for particles
    ParticleBufferObject particleBufferObject;
    particleBufferObject.initialize(n);

    // Set up render objects for traces
    TraceBufferObject traceBufferObject;
    traceBufferObject.initialize(n);
    
    // force rendering to take depth into account
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    
    // mainloop variables
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

        // Load new simulation timestep if dt time has passed
        if (playSim && time - prevSimFrameTime > dt) {
            prevSimFrameTime = time;
            currentParticleSet = (currentParticleSet + 1) % NUM_FILES;
            updateParticles(currentParticleSet, bodies, combinedTraces);
        }
        
        // Reset each frame
        glClearColor(0.05f, 0.05f, 0.08f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // Camera set up
        basicShader.use();
        basicShader.setMat4("view", camera.view);
        basicShader.setMat4("projection", camera.projection);
        particleShader.use();
        particleShader.setMat4("view", camera.view);
        particleShader.setMat4("projection", camera.projection);
        
        
        // Draw particles using parallel instancing
        std::vector<ParticleBufferObject::ParticleInstanceData> instances;
        instances.reserve(n);

        for (int i = 0; i < n; i++)
        {
            glm::mat4 model(1.0f);
            model = glm::translate(model, particlePositionVec(bodies.position[i]));
            float radius = scale * bodies.mass[i] / 20.0f;
            model = glm::scale(model, glm::vec3(radius));
            glm::vec3 color(i*1.0f/n, 0.0f, 1.0f); 

            instances.push_back({ model, color });
        }

        particleBufferObject.loadNewData(instances);
        particleBufferObject.draw(instances.size());


        // Draw trace using parallel instancing
        if (!disableTrace) {
            basicShader.use();
            glm::mat4 model = glm::mat4(1.0f);
            basicShader.setMat4("model", model);
            basicShader.setVec3("color", glm::vec3(1.0f, 1.0f, 1.0f));
            traceBufferObject.loadNewData(combinedTraces);
            traceBufferObject.draw(combinedTraces.size());
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
            updateParticles(currentParticleSet, bodies, combinedTraces);
        }
        
        ImGui::Checkbox("Loop Simulation", &playSim);
        
        ImGui::Spacing();
        ImGui::Text("Visual Settings");
        ImGui::Checkbox("Disable Trace", &disableTrace);
        ImGui::SliderFloat("Particle Scale", &scale, 0.1f, 20.0f, "ratio = %.1f");
        
        ImGui::Spacing();
        ImGui::Text("Movement");
        if (ImGui::SliderFloat("Zoom", &zoom, 0.5f, 10.0f, "zoom = %.1f")) {
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
