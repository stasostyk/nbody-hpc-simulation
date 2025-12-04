#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <vector>
#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "../shaders/shader.hpp"

void createSphere(std::vector<float>& vertices,
                  std::vector<unsigned int>& indices,
                  unsigned int X_SEGMENTS = 32,
                  unsigned int Y_SEGMENTS = 32);

class Particle {
public:
    glm::vec3 position;
    glm::vec3 color;
    float radius;

    GLuint VAO, VBO, EBO;
    GLsizei indexCount;

    Particle(glm::vec3 pos, glm::vec3 col, float r)
        : position(pos), color(col), VAO(0), VBO(0), EBO(0), indexCount(0), radius(r) {}

    void init();
    void draw(Shader& shader);
};

#endif
