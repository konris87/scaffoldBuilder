#version 330 core

// uniforms for lighting options
uniform float Ka; 
uniform vec3 lightColor;
uniform vec4 objectColor;

// inputs from vertex shader
in vec4 vertexPosWorld;
in vec4 vertexPosCamera;
in vec4 vertexNormCamera;
in vec4 lightPosCamera;

// output: fragment color
out vec4 FragColor;

void main()
{
	
	// 1. ambient light
	vec4 ambientColor = vec4(1.0, 1.0, 1.0, 1.0); 
	vec4 diffuseColor = vec4(1.0, 1.0, 1.0, 1.0); 
	vec4 specularColor = vec4(1.0, 1.0, 1.0, 1.0); 
	
	ambientColor = Ka * vec4(lightColor, 1.0);
	
	// 2. diffuse color
	// 2.1 ensure that the vertex normals from the vertex shader is normal :)
	vec4 normals = normalize(vec4(vertexNormCamera.xyz, 0.0));
	
	// 2.2 normalized light direction
	vec4 direction = normalize(lightPosCamera - vertexPosCamera);
	
	// 2.3 compute the cos using the dot product
	float cosTheta = max(dot(normals, direction), 0.0);
	// cosTheta = clamp(cosTheta, 0, 1);
	
	// calculate diffuse
	diffuseColor = vec4(lightColor, 1.0) * cosTheta;	
	
	// final fragment color
	vec4 outColor = (ambientColor + diffuseColor) * vec4(objectColor.xyz, 1.0);
	FragColor = vec4(outColor.xyz, objectColor.w);

}