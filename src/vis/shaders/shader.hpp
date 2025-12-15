#pragma once

#include <string>
#include <glm/glm.hpp>
#include <glad/glad.h>

class Shader {
public:
    unsigned int id;

    Shader(const std::string& vertexPath, const std::string& fragmentPath);
    ~Shader() = default;

    void use() const;

    void setMat4(const std::string& name, const glm::mat4& mat) const;
    void setVec3(const std::string& name, const glm::vec3& v) const;
    void setFloat(const std::string& name, float value) const;
    void setInt(const std::string& name, int value) const;

private:
    static std::string loadFileToString(const std::string& path);
    static unsigned int compile(GLenum type, const std::string& src);
};
