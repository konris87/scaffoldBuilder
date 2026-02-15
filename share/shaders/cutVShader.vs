#version 330 core
layout (location = 0) in vec3 planePos;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

out vec4 vertexPosWorld;
uniform vec3 minBounds;
uniform vec3 maxBounds;

void main()
{
	// position of fragments in clip space
	// gl_Position = projection * view * model * vec4(planePos, 1.0);
	vertexPosWorld = model * vec4(planePos, 1.0);
	
    // Clamp the vertex position to the bounding box
    vertexPosWorld.x = clamp(vertexPosWorld.x, minBounds.x, maxBounds.x);
    vertexPosWorld.y = clamp(vertexPosWorld.y, minBounds.y, maxBounds.y);
    vertexPosWorld.z = clamp(vertexPosWorld.z, minBounds.z, maxBounds.z);
	
	// gl_Position = projection * view * worldPos;
	gl_Position = projection * view * vertexPosWorld;
}