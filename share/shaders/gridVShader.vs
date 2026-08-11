#version 330 core

// -----------------------------------------------------------------------------
// Infinite ground-plane grid - vertex stage.
//
// Instead of a world-space quad (which clips against the near plane once it gets
// large), we draw a screen-covering quad in clip space and, for each corner,
// unproject one point on the near plane and one on the far plane into world
// space. The fragment shader turns those two points into a per-pixel view ray
// and intersects it with the ground plane, so the plane is mathematically
// infinite and there is no geometry to clip.
//
// The quad is generated from gl_VertexID (Pos[]/Indices[]), so the Grid object
// keeps issuing 6 attribute-less vertices via glDrawArrays(GL_TRIANGLES, 0, 6).
// -----------------------------------------------------------------------------

uniform mat4 view;
uniform mat4 projection;

out vec3 nearPoint;   // world-space point on the near plane for this pixel
out vec3 farPoint;    // world-space point on the far plane  for this pixel

const vec2 Pos[4] = vec2[4](
	vec2(-1.0, -1.0),   // bottom left
	vec2( 1.0, -1.0),   // bottom right
	vec2( 1.0,  1.0),   // top right
	vec2(-1.0,  1.0)    // top left
);

// CCW winding so the full-screen quad is front-facing (backface culling is on
// when the grid draws). Pos is laid out CCW: BL, BR, TR, TL.
const int Indices[6] = int[6](0, 1, 2, 0, 2, 3);

vec3 unproject(vec2 ndc, float ndcZ, mat4 invViewProj)
{
	vec4 p = invViewProj * vec4(ndc, ndcZ, 1.0);
	return p.xyz / p.w;
}

void main()
{
	vec2 ndc = Pos[Indices[gl_VertexID]];

	mat4 invViewProj = inverse(projection * view);

	// Default OpenGL NDC depth: near plane = -1, far plane = +1.
	nearPoint = unproject(ndc, -1.0, invViewProj);
	farPoint  = unproject(ndc,  1.0, invViewProj);

	// The quad is already in NDC and covers the whole screen; pass it through.
	// (Depth here is irrelevant - the fragment shader writes gl_FragDepth.)
	gl_Position = vec4(ndc, 0.0, 1.0);
}
