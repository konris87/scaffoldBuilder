#version 330 core
out vec4 FragColor;

uniform vec4 cutPlaneCoeffs;
uniform bool cutPlane;
uniform vec4 edgeColor;

//uniform vec4 vertexColor;
in vec4 vertexPosWorld;

void main()
{
    // FragColor = vec4(outColor, 1.0f);
	
	if (cutPlane){
		float d = dot(cutPlaneCoeffs, vertexPosWorld);
	
		if (d < 0){
			
			discard;
		
		}
	}
	
	FragColor = edgeColor;
} 