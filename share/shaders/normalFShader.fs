#version 330 core

uniform vec3 normalColor;

out vec4 FragColor;

void main()
{
    // FragColor = vec4(1.0, 0.5, 0.2, 1.0);
    FragColor = vec4(normalColor, 1.0f);
}  