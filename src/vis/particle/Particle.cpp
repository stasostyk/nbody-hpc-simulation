#include "Particle.hpp"
#include <cmath>

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

void Particle::init()
{
    std::vector<float> sphereVerts;
    std::vector<unsigned int> sphereInds;
    createSphere(sphereVerts, sphereInds, 32, 32);
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
}

void Particle::draw(Shader& shader)
{
    glm::mat4 model = glm::mat4(1.0f);
    model = glm::translate(model, position);
    model = glm::scale(model, glm::vec3(radius));

    shader.setMat4("model", model);
    shader.setVec3("color", color);

    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, indexCount, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}


void Particle::drawTrace(Shader& shader)
{
    if (trace.size() < 2)
        return;

    glm::mat4 model = glm::mat4(1.0f);
    shader.setMat4("model", model);

    shader.setVec3("color", glm::vec3(1.0f, 1.0f, 1.0f));

    std::vector<float> verts;
    verts.reserve(trace.size() * 3);

    for (auto& p : trace) {
        verts.push_back(p.x);
        verts.push_back(p.y);
        verts.push_back(p.z);
    }

    GLuint traceVAO, traceVBO;
    glGenVertexArrays(1, &traceVAO);
    glGenBuffers(1, &traceVBO);

    glBindVertexArray(traceVAO);
    glBindBuffer(GL_ARRAY_BUFFER, traceVBO);
    glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float),
                 verts.data(), GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glDrawArrays(GL_LINE_STRIP, 0, trace.size());

    glDeleteBuffers(1, &traceVBO);
    glDeleteVertexArrays(1, &traceVAO);
}