#include "ParticleBufferObject.hpp"

void ParticleBufferObject::initialize(int n) {
    // Set up particle sphere geometry
    std::vector<float> sphereVerts;
    std::vector<unsigned int> sphereInds;
    createSphere(sphereVerts, sphereInds, 8, 8);
    indexCount = (GLsizei)sphereInds.size();

    // Set up buffer objects
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &meshVBO);
    glGenBuffers(1, &EBO);
    glGenBuffers(1, &instanceVBO);

    glBindVertexArray(VAO);

    // Sphere mesh positions
    glBindBuffer(GL_ARRAY_BUFFER, meshVBO);
    glBufferData(GL_ARRAY_BUFFER,
                sphereVerts.size() * sizeof(float),
                sphereVerts.data(),
                GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

    // Sphere indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                sphereInds.size() * sizeof(unsigned int),
                sphereInds.data(),
                GL_STATIC_DRAW);

    // Instance data
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferData(GL_ARRAY_BUFFER,
                n * sizeof(ParticleInstanceData),
                nullptr,
                GL_DYNAMIC_DRAW);

    // mat4 occupies locations 3,4,5,6
    std::size_t offset = 0;
    for (int i = 0; i < 4; ++i) {
        glEnableVertexAttribArray(3 + i);
        glVertexAttribPointer(3 + i, 4, GL_FLOAT, GL_FALSE,
            sizeof(ParticleInstanceData), (void*)offset);
        glVertexAttribDivisor(3 + i, 1);
        offset += sizeof(glm::vec4);
    }

    // color at location 7
    glEnableVertexAttribArray(7);
    glVertexAttribPointer(7, 3, GL_FLOAT, GL_FALSE,
        sizeof(ParticleInstanceData), (void*)offsetof(ParticleInstanceData, color));
    glVertexAttribDivisor(7, 1);

    glBindVertexArray(0);
}

void ParticleBufferObject::loadNewData(std::vector<ParticleInstanceData>& instances) {
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0,
        instances.size() * sizeof(ParticleInstanceData),
        instances.data());
}

void ParticleBufferObject::draw(int n) {
    glBindVertexArray(VAO);
    glDrawElementsInstanced(
        GL_TRIANGLES,
        indexCount,
        GL_UNSIGNED_INT,
        0,
        n
    );
}

void ParticleBufferObject::createSphere(std::vector<float>& vertices,
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