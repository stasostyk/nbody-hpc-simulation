#include "Particle.hpp"

void Particle::draw(Shader& shader, float scale)
{
    glm::mat4 model = glm::mat4(1.0f);
    model = glm::translate(model, position);
    model = glm::scale(model, glm::vec3(radius * scale));

    shader.setMat4("model", model);
    shader.setVec3("color", color);

    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, indexCount, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}