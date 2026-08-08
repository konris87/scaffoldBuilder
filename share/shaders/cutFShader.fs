#version 330 core

// input from vertexshader
in vec4 vertexPosWorld;
uniform vec4 planeColor;
uniform vec3 minBounds;
uniform vec3 maxBounds;

// output: fragment color
out vec4 FragColor;

void main(){
		
	if (vertexPosWorld.x < minBounds.x || vertexPosWorld.x > maxBounds.x ||
		vertexPosWorld.y < minBounds.y || vertexPosWorld.y > maxBounds.y ||
		vertexPosWorld.z < minBounds.z || vertexPosWorld.z > maxBounds.z)
		{
			discard;
		} 
		
	FragColor = planeColor;

}