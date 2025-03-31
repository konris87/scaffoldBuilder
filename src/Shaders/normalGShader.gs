#version 330 core
layout (triangles) in;
layout (line_strip, max_vertices = 2) out;

// we need the projection matrix
uniform mat4 projection;

// set up some boolean flags for mesh features we wish to render
// uniform boolean showFaceNormals;
// uniform float normalLength;

// write a function that estimates the normal of each triangle
vec3 get_normal(){
	
	vec3 v1 = vec3(gl_in[0].gl_Position) - vec3(gl_in[1].gl_Position);
	vec3 v2 = vec3(gl_in[0].gl_Position) - vec3(gl_in[2].gl_Position);
	
	return normalize(cross(v1, v2));
}


// function for estimating the centroid to draw face normals

vec3 get_centroid() {

    vec3 c = (vec3(gl_in[0].gl_Position) + vec3(gl_in[1].gl_Position) + vec3(gl_in[2].gl_Position)) / 3.0;

    return c;
}

void main(){
	
	vec3 center = get_centroid();
	vec3 normal = get_normal();

    // Emit the centroid
    gl_Position = projection * vec4(center, 1.0);
    EmitVertex();
    
    // Emit a point along the normal direction for the face normal line
    gl_Position = projection * (vec4(center, 1.0) + vec4(normal, 0.0) * 0.05);  // Scale the normal as needed
    EmitVertex();
    
    EndPrimitive();
};
