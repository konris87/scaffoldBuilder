#version 330 core
layout (location = 0) in vec3 aPos;

out VS_OUT {
    vec4 posWorld;
} vs_out;

uniform mat4 view;
uniform mat4 projection;
uniform mat4 model;

void main()
{
	vs_out.posWorld = model * vec4(aPos.x, aPos.y, aPos.z, 1.0); 
    gl_Position = projection * view * model * vec4(aPos.x, aPos.y, aPos.z, 1.0); 
}