#pragma once

#include "RenderBufferObject.hpp"
#include <glm/glm.hpp>
#include <vector>
#include <glad/glad.h>

class ParticleBufferObject : public RenderBufferObject
{
public:
    struct ParticleInstanceData {
        glm::mat4 model;
        glm::vec3 color;
    };
    

    void initialize(int n);
    void loadNewData(std::vector<ParticleInstanceData>& instances);
    void draw(int n);

private:
    void createSphere(std::vector<float>& vertices,
                  std::vector<unsigned int>& indices,
                  unsigned int X_SEGMENTS,
                  unsigned int Y_SEGMENTS);
    GLuint VAO;
    GLuint meshVBO;
    GLuint EBO;
    GLuint instanceVBO;
    GLsizei indexCount;
};
