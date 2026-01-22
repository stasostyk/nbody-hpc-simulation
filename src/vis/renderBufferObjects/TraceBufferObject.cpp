#include "TraceBufferObject.hpp"

const int NUM_FILES = 99;

void TraceBufferObject::initialize(int n) {
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
    
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void TraceBufferObject::loadNewData(std::vector<glm::vec3>& combinedTraces) {
    int traceCount = combinedTraces.size();
    glBindBuffer(GL_ARRAY_BUFFER, traceVBO);
    glBufferSubData(
        GL_ARRAY_BUFFER,
        0,
        traceCount * sizeof(glm::vec3),
        combinedTraces.data()
    );
}

void TraceBufferObject::draw(int n) {
    glBindVertexArray(traceVAO);
    glDrawArrays(GL_POINTS, 0, n);
}