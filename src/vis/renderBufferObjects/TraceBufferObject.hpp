#pragma once

#include "RenderBufferObject.hpp"
#include <glm/glm.hpp>
#include <vector>
#include <glad/glad.h>

class TraceBufferObject : public RenderBufferObject
{
public:

    void initialize(int n);
    void loadNewData(std::vector<glm::vec3>& combinedTraces);
    void draw(int n);

private:
    GLuint traceVAO;
    GLuint traceVBO;
};
