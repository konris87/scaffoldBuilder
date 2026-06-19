#version 330 core

// input from vertexshader
in vec4 vertexPosWorld;

uniform vec3 minBounds;
uniform vec3 maxBounds;
uniform vec4 planeColor;

// output: fragment color
out vec4 FragColor;

void main(){
		
//	if (vertexPosWorld.x < minBounds.x || vertexPosWorld.x > maxBounds.x ||
//		vertexPosWorld.y < minBounds.y || vertexPosWorld.y > maxBounds.y ||
//		vertexPosWorld.z < minBounds.z || vertexPosWorld.z > maxBounds.)
//		{
//			discard;
//		} 
		
	FragColor = planeColor;

}