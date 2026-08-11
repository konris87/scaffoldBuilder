#version 330 core
layout(location = 0) out vec4 FragColor;

// Per-pixel view ray endpoints (world space), from the VS.
in vec3 nearPoint;
in vec3 farPoint;

// Needed to reproject the ground-plane hit for correct depth.
uniform mat4 view;
uniform mat4 projection;

// Live camera world position (used only for the horizon fade).
uniform vec3 camWorldPos;

// Radius (world units) at which the grid fades toward the horizon. This no
// longer bounds any geometry - the plane is infinite - it only controls where
// the grid dissolves so the far moire (cells below a pixel) is hidden. Lower it
// if you ever see line wobble far from the origin (float precision in mod()).
uniform float gGridSize = 1000.0;   // fade radius < far (100), not 10000

// Size of the finest grid cell in world units (1.0 = one line per unit).
uniform float gGridCellSize = 1.0;

// Minimum pixels between visible cells before the LOD fades the finer grid out.
uniform float gGridMinPixelsBetweenCells = 2.0;

// Line thickness, roughly in pixels (dudv is world-units-per-pixel, so this
// scales the anti-aliased line width). The in-plane axes reuse it, so they
// scale together. ~1.0-2.0 = hairline, 4.0 = default, higher = bold.
uniform float gGridLineWidth = 4.0;

uniform vec4 gGridColorThin  = vec4(0.5, 0.5, 0.5, 1.0);
uniform vec4 gGridColorThick = vec4(0.2, 0.2, 0.2, 1.0);

// --- helpers -----------------------------------------------------------------
float satf(float x) { return clamp(x, 0.0, 1.0); }
vec2  satv(vec2  x) { return clamp(x, vec2(0.0), vec2(1.0)); }
float max2(vec2  v) { return max(v.x, v.y); }
float log10(float x) { return log(x) / log(10.0); }

// -----------------------------------------------------------------------------
// PLANE CHOICE: the ground plane is y = 0 (the XZ plane). To use z = 0 (XY):
//   - intersect on .z:  float t = -nearPoint.z / (farPoint.z - nearPoint.z);
//   - set             vec2 P = WorldPos.xy;   at the line marked [PLANE]
//   - the in-plane axes then become X = red / Y = green (see the axis block).
// The out-of-plane (vertical) axis cannot be produced here - draw it separately.
// -----------------------------------------------------------------------------

void main()
{
	// --- reconstruct the infinite ground plane for this pixel ----------------
	// The view ray runs nearPoint -> farPoint. Find where it crosses y = 0.
	// NOTE: we compute WorldPos and the grid unconditionally (WorldPos is a
	// continuous function of the ray), then discard below - this keeps the
	// screen-space derivatives well-defined for pixels straddling the horizon.
	// Component of the ray along the plane normal. When the view is parallel to
	// the ground (edge-on), this is ~0 and there is no intersection - bail out
	// before dividing by it, otherwise the NaN propagates and satf(NaN) lights
	// the whole screen red via the axis test. NOTE: the DEFAULT camera looks
	// straight down -z with the grid on y=0, i.e. exactly edge-on, so at startup
	// every pixel hits this guard and the grid is invisible until you orbit to
	// look down at it (see the message accompanying this change).
	float denom = farPoint.y - nearPoint.y;
	if (abs(denom) < 1e-6) discard;
	float t = -nearPoint.y / denom;
	vec3 WorldPos = nearPoint + t * (farPoint - nearPoint);

	vec2 P = WorldPos.xz;                        // [PLANE] in-plane coordinates

	// Screen-space derivative of the plane coords: how many world units one
	// pixel spans here. Makes line width and AA resolution/distance-independent.
	vec2 dudv = vec2(
		length(vec2(dFdx(P.x), dFdy(P.x))),
		length(vec2(dFdx(P.y), dFdy(P.y)))
	);
	float pixelUnit = length(dudv);

	// Continuous level of detail: as cells shrink below a few pixels we step up
	// to the next 10x-coarser grid. LOD is fractional so we can cross-fade.
	float LOD = max(0.0,
		log10(pixelUnit * gGridMinPixelsBetweenCells / gGridCellSize) + 1.0);
	float LODfade = fract(LOD);

	float cellLOD0 = gGridCellSize * pow(10.0, floor(LOD));
	float cellLOD1 = cellLOD0 * 10.0;
	float cellLOD2 = cellLOD1 * 10.0;

	// Scale the derivative to set line thickness (~pixels). See gGridLineWidth.
	dudv *= gGridLineWidth;

	// Coverage of a gridline for each LOD: 1 on a line, 0 between.
	vec2 m0 = mod(P, cellLOD0) / dudv;
	float lod0a = max2(vec2(1.0) - abs(satv(m0) * 2.0 - vec2(1.0)));
	vec2 m1 = mod(P, cellLOD1) / dudv;
	float lod1a = max2(vec2(1.0) - abs(satv(m1) * 2.0 - vec2(1.0)));
	vec2 m2 = mod(P, cellLOD2) / dudv;
	float lod2a = max2(vec2(1.0) - abs(satv(m2) * 2.0 - vec2(1.0)));

	// Compose the visible LODs: coarse lines stay solid ("thick"), fine lines
	// fade in/out with the fractional LOD ("thin").
	vec4 Color;
	if (lod2a > 0.0) {
		Color = gGridColorThick;
		Color.a *= lod2a;
	} else if (lod1a > 0.0) {
		Color = mix(gGridColorThick, gGridColorThin, LODfade);
		Color.a *= lod1a;
	} else {
		Color = gGridColorThin;
		Color.a *= lod0a * (1.0 - LODfade);
	}

	// --- in-plane axes (X = red, Z = blue for the XZ plane) ------------------
	// An axis is just "the line where the other coordinate is ~0". Reuse the
	// same derivative width so axes match the gridlines and are anti-aliased.
	float xAxis = 1.0 - satf(abs(P.y) / dudv.y);   // P.y == world z == 0 -> X axis
	float zAxis = 1.0 - satf(abs(P.x) / dudv.x);   // P.x == world x == 0 -> Z axis
	if (zAxis > 0.0) {                             // Z axis (blue)
		Color.rgb = vec3(0.0, 0.0, 1.0);
		Color.a   = max(Color.a, zAxis);
	}
	if (xAxis > 0.0) {                             // X axis (red) on top
		Color.rgb = vec3(1.0, 0.0, 0.0);
		Color.a   = max(Color.a, xAxis);
	}

	// --- horizon fade --------------------------------------------------------
	float dist = length(P - camWorldPos.xz);
	float falloff = 1.0 - satf(dist / gGridSize);
	Color.a *= falloff;

	if (Color.a <= 0.0) discard;

	// --- projection-agnostic visibility + depth ------------------------------
	// The old `if (t <= 0.0) discard;` was PERSPECTIVE-only: it assumed every ray
	// starts at the eye, so t <= 0 meant "behind the camera". This camera is
	// ORTHOGRAPHIC by default (TrackBall mode = Ortho), where rays are parallel
	// and the near-plane is a large rectangle whose lower half sits BELOW the
	// y = 0 ground - so foreground pixels legitimately give t < 0 and were being
	// wrongly clipped (the bottom-of-screen band, growing with zoom because the
	// ortho height scales with the camera distance/radius).
	//
	// Reproject the hit and discard only when it is genuinely behind a
	// perspective eye (clip.w <= 0). In ortho clip.w is always 1, so the
	// foreground is kept; clamp the depth so the near/far planes can't clip it.
	vec4 clip = projection * view * vec4(WorldPos, 1.0);
	if (clip.w <= 0.0) discard;
	gl_FragDepth = clamp(0.5 + 0.5 * (clip.z / clip.w), 0.0, 1.0);

	FragColor = Color;
}
