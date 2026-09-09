#version 330 core
layout (location = 0) in vec3 aPos;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

void main() {
    // 1. Get the position in View Space
    vec4 viewPos = view * model * vec4(aPos, 1.0);
    
    // 2. Project it
    gl_Position = projection * viewPos;

    // 3. The Trick: Force the Z to be exactly at the far clipping plane 
    // but keep the X and Y projection logic. 
    // This makes the line look like it touches the horizon.
    gl_Position.z = gl_Position.w; 
}