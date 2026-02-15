#version 330 core
layout (location = 0) in vec3 vertexPosModel;
layout (location = 1) in vec3 vertexNormModel;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

// the position of light in world space
uniform vec3 lightPosWorld;

// outputs for fragment shader
// 1. position of vertices in camera space
out vec4 vertexPosCamera;

// 2. vertex normals in camera space
out vec4 vertexNormCamera;

// 3. position of light in camera space
out vec4 lightPosCamera;

// 4. position of vertices in world
out vec4 vertexPosWorld;

void main()
{
	// 1. position of vertices in camera space
	vertexPosCamera = view * model * vec4(vertexPosModel, 1.0);
	
	// 2. vertex normals in camera space
	mat3 normalMatrix = transpose(inverse(mat3(view * model)));
	vertexNormCamera = vec4(normalMatrix * vertexNormModel, 0.0);
	
	// 3. position of light in camera space
	lightPosCamera = vec4(lightPosWorld, 1.0);
	
	// 4.
	vertexPosWorld = model * vec4(vertexPosModel, 1.0);
	
	// position of fragments in clip space
	gl_Position = projection * vertexPosCamera;
}