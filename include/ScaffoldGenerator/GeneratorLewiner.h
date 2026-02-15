/**
 * @file    MarchingCubes.h
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @author  Math Dept, PUC-Rio
 * @version 0.2
 * @date    12/08/2002
 * 
 * @brief   MarchingCubes Algorithm
 */
//________________________________________________

#ifndef GENERATOR_LEWINER_H
#define GENERATOR_LEWINER_H

#include <vector>
#include <array>
#include <atomic>
#include <string>
#include <queue>

#include "SeedGenerator/Container.h"
#include "Math/Vec.h"

typedef struct
{
  float  x,  y,  z ;  /**< Vertex coordinates */
  float nx, ny, nz ;  /**< Vertex normal */
} LVertex ;

typedef struct
{
  int v1,v2,v3 ;
  Vec3 normal;
  /**< Triangle vertices */
} LTriangle ;
//____________

struct Vec3i { int x, y, z; };

class GeneratorLewiner{
public:

	GeneratorLewiner() {};
	
	GeneratorLewiner(
		const std::vector<Vec3>& seeds,
		const std::array<float, 6>& bounds,
		const std::array<int, 3>& dims,
		const float threshold = 0.0f,
		const float isoLevel = 0.1f,
		const int foam = 0
	);
	
	~GeneratorLewiner() {};

// functions
public:
	void compute_scalar_field(const IContainer& con);

	void marching_cubes();

	void export_stl(std::string fileName);

	//--------------------------------
	// rendering
	void draw();
	void draw_edges();

private:
	void smooth_scalar_field();

	void seal_grid_boundaries();

	void remove_isolated_islands();

	void compute_intersection_points();

	void process_cube(int i, int j, int k, const float cube[8], int lut_entry);

	void update_steps();

	int find_vertex_index(int x, int y, int z);
	
	Vec3 get_position(int x, int y, int z);

	float smin(float a, float b, float k);

	void add_triangle(const uint8_t* trig, int i, int j, int k, int n, int v12 = -1);

	LVertex add_x_vertex(
		const int& i, const int& j, const int& k,
		const float& val0, const float& val1);
	LVertex add_y_vertex(
		const int& i, const int& j, const int& k,
		const float& val0, const float& val2);
	LVertex add_z_vertex(
		const int& i, const int& j, const int& k,
		const float& val0, const float& val3);

	float get_data(const int i, const int j, const int k) const;

	float get_x_grad(const int i, const int j, const int k);

	float get_y_grad(const int i, const int j, const int k);
	
	float get_z_grad(const int i, const int j, const int k);

	int get_x_vert(int i, int j, int k);

	int get_y_vert(int i, int j, int k);
	
	int get_z_vert(int i, int j, int k);

	int add_c_vertex(const int i, const int j, const int k);

	bool test_face(signed char face, const float _cube[8]);

	bool modified_test_interior(signed char s, int _case, int _config, const float _cube[8]);

	int interior_ambiguity(int amb_face, int s, const float _cube[8]);

	int interior_ambiguity_verification(int edge, const float _cube[8]);

	bool interior_test_case13(const float _cube[8]);

	bool interior_test_case13_2(float isovalue, const float _cube[8], int& tunnelOrientation);

	void validate_topology();
	//-------------------------------
	// rendering
	void _setup_mesh();
	void _setup_edges();
	void _update_render();
	void add_edge(int idx1, int idx2, int idx3);

// members
public:
	void* container = nullptr;
	void* generator = nullptr;
	bool foam = false;
	std::string name = "";
	std::array<float, 4> color = { 0.5f, 0.5f, 0.5f, 1.0f };
	bool hidden = false;

private:	
	const std::vector<Vec3> seeds;
	float stepX, stepY, stepZ;
	std::vector<float> scalarField;
	std::array<int, 3> blockDims;
	std::array<float, 6> bounds;
	float isoLevel{ 0.1f };
	std::vector<int> x_verts;
	std::vector<int> y_verts;
	std::vector<int> z_verts;

	std::vector<LTriangle> meshTriangles;
	std::vector<LVertex> meshVertices;

	std::atomic<int> vertexCount{ 0 };
	std::atomic<int> triangleCount{ 0 };

	// scaffold properties
	float threshold = { 0.5f };

	// ----------------------------------------------------
	// these are for opengl
	std::vector<float> vertices;
	std::vector<float> normals;
	std::set<std::pair<unsigned int, unsigned int>> edgeSet;
	// these are for rendering mesh edges
	std::vector<unsigned int> edgeIndices;
	// this is for rendering the EBO
	std::vector<unsigned int> indices;
	unsigned int VAO{ 0 }, VBO{ 0 }, EBO{ 0 };
	unsigned int edgeVAO{ 0 }, edgeVBO{ 0 }, edgeEBO{ 0 }, normalsVBO{ 0 };
	unsigned int tortuosityPathVAO{ 0 }, tortuosityPathVBO{ 0 };
};

#endif // ! "GENERATOR_LEWINER_H"
