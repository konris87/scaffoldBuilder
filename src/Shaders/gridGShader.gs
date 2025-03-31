#version 330 core
layout(points) in;
layout(line_strip, max_vertices = 200) out;

float resolutionX = 1;
float resolutionY = 1;
int gridSizeX = 20;
int gridSizeY = 20;

in VS_OUT {
    vec4 posWorld;
} gs_in[];

uniform mat4 projection;
uniform mat4 view;

void main(){
	
	int stepsX = int(ceil(2 * gridSizeX / resolutionX)); 	
	int stepsY = int(ceil(2 * gridSizeY / resolutionY)); 	
	
	vec4 center = gs_in[0].posWorld;
	
	for (int i = 0 ; i < stepsX + 1; i++){
	
		vec4 p1 = center +  vec4(-gridSizeX, gridSizeY - (i * resolutionY), 0.0, 0.0);
		//gl_Position = gl_in[0].gl_Position + vec4(-gridSizeX, gridSizeY - (i * resolutionY), 0.0, 0.0);
		gl_Position = projection * view * p1;
		
		EmitVertex();
	
		vec4 p2 = center + vec4(gridSizeX, gridSizeY - (i * resolutionY), 0.0, 0.0);
		// gl_Position = gl_in[0].gl_Position + vec4(gridSizeX, gridSizeY - (i * resolutionY), 0.0, 0.0);
		gl_Position = projection * view * p2;
		EmitVertex();
	
		EndPrimitive();
		
	}
	
	for (int i = 0 ; i < stepsY + 1; i++){
		
		vec4 p1 = center + vec4(gridSizeX - (i * resolutionX), -gridSizeY, 0.0, 0.0);
		//gl_Position = gl_in[0].gl_Position + vec4(gridSizeX - (i * resolutionX), -gridSizeY, 0.0, 0.0);
		
		gl_Position = projection * view * p1;
		EmitVertex();
	
		vec4 p2 = center + vec4(gridSizeX - (i * resolutionX), gridSizeY, 0.0, 0.0);
		//gl_Position = gl_in[0].gl_Position + vec4(gridSizeX - (i * resolutionX), gridSizeY, 0.0, 0.0);
		gl_Position = projection * view * p2;
		EmitVertex();
	
		EndPrimitive();
		
	}
	
	// horizontal lines
	//for (int i = -gridSizeX ; i < stepsX + 1; i+=resolutionX) {
		
		//gl_Position = gl_in[0].gl_Position + vec4(0.0, i * resolutionX, 0.0, 0.0);
		//EmitVertex();
			
		//gl_Position = gl_in[0].gl_Position + vec4(200.0, i * resolutionX, 0.0, 0.0);
		//EmitVertex();
			
		//EndPrimitive();
		
	//}
	
	// vertical lines
	//for (int i = 0 ; i < stepsY + 1; ++i) {
		
		//gl_Position = gl_in[0].gl_Position + vec4(i * resolutionY, 0.0, 0.0, 0.0);
		//EmitVertex();
			
		//gl_Position = gl_in[0].gl_Position + vec4(i * resolutionY, 200.0, 0.0, 0.0);
		//EmitVertex();
			
		//EndPrimitive();
		
	//}	
}