#include "shader.hpp"

#include <fstream>
#include <sstream>
#include <iostream>

std::string Shader::loadFileToString(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        std::cerr << "Failed to open shader file: " << path << std::endl;
        return "";
    }
    std::stringstream ss;
    ss << file.rdbuf();
    return ss.str();
}

unsigned int Shader::compile(GLenum type, const std::string& src) {
    unsigned int shader = glCreateShader(type);
    const char* csrc = src.c_str();
    glShaderSource(shader, 1, &csrc, nullptr);
    glCompileShader(shader);

    int success = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        int length = 0;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
        std::string log(length, '\0');
        glGetShaderInfoLog(shader, length, nullptr, log.data());
        std::cerr << "Shader compile error:\n" << log << std::endl;
    }

    return shader;
}

Shader::Shader(const std::string& vertexPath, const std::string& fragmentPath) {
    std::string vertSrc = loadFileToString(vertexPath);
    std::string fragSrc = loadFileToString(fragmentPath);

    unsigned int v = compile(GL_VERTEX_SHADER, vertSrc);
    unsigned int f = compile(GL_FRAGMENT_SHADER, fragSrc);

    id = glCreateProgram();
    glAttachShader(id, v);
    glAttachShader(id, f);
    glLinkProgram(id);

    int success = 0;
    glGetProgramiv(id, GL_LINK_STATUS, &success);
    if (!success) {
        int length = 0;
        glGetProgramiv(id, GL_INFO_LOG_LENGTH, &length);
        std::string log(length, '\0');
        glGetProgramInfoLog(id, length, nullptr, log.data());
        std::cerr << "Program link error:\n" << log << std::endl;
    }

    glDeleteShader(v);
    glDeleteShader(f);
}

void Shader::use() const {
    glUseProgram(id);
}

void Shader::setMat4(const std::string& name, const glm::mat4& mat) const {
    int loc = glGetUniformLocation(id, name.c_str());
    glUniformMatrix4fv(loc, 1, GL_FALSE, &mat[0][0]);
}

void Shader::setVec3(const std::string& name, const glm::vec3& v) const {
    int loc = glGetUniformLocation(id, name.c_str());
    glUniform3fv(loc, 1, &v[0]);
}

void Shader::setFloat(const std::string& name, float value) const {
    int loc = glGetUniformLocation(id, name.c_str());
    glUniform1f(loc, value);
}

void Shader::setInt(const std::string& name, int value) const {
    int loc = glGetUniformLocation(id, name.c_str());
    glUniform1i(loc, value);
}
