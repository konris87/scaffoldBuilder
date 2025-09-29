#version 330 core
out vec4 FragColor;

//uniform vec4 vertexColor;
in vec3 outColor;

void main()
{
    FragColor = vec4(0.0f, 1.0f, 0.0f, 1.0f);
    // FragColor = vec4(outColor, 1.0f);
} 