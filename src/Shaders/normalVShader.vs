#version 330 core
layout (location = 0) in vec3 vertexPosModel;

// view and model matrix uniforms
uniform mat4 view;
uniform mat4 model;

void main(){
	gl_Position = view * model * vec4(vertexPosModel, 1.0);
}
