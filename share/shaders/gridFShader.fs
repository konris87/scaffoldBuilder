#version 330 core
out vec4 FragColor;

uniform vec3 gridColor;

void main()
{
	FragColor = vec4(gridColor, 1.0);
    // FragColor = vec4(1.0, 1.0, 0.0, 1.0);   
}  