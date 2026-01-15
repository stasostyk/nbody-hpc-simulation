#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 3) in mat4 instanceModel;
layout(location = 7) in vec3 instanceColor;

out vec3 fragColor;

uniform mat4 view;
uniform mat4 projection;

void main() {
    fragColor = instanceColor;
    gl_Position = projection * view * instanceModel * vec4(aPos, 1.0);
}