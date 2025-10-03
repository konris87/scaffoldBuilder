#version 330 core
out vec4 FragColor;

uniform vec4 outColor;

void main()
{
	FragColor = vec4(outColor.x, outColor.y, outColor.z, 1.0f);
} 