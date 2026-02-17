#ifndef GENERATOR_H // include guard
#define GENERATOR_H

// opengl
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>

#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <filesystem>
#include <set>
#include <string>
#include <cmath>
#include <iostream>
#include <voro++.hh>
#include <Math/Kdtree.h>
#include <Math/Vec.h>
#include <Utils/Utils.h>
#include "OpenGlRender/Model.h"
#include "SeedGenerator/Container.h"

struct Metrics {
    float surfaceArea{ 0.0f };
    float porosity{ 0.0f };
    float volumeFraction{ 0.0f };
    float volume{ 0.0f };
    float tortuosity{ 0.0f };
    float localThickness{ 0.0f };
    float localThicknessStd{ 0.0f };
    float localSeperation{ 0.0f };
    float localSeperationStd{ 0.0f };
};

// hash function to compute hash from key
struct GlobalVertexVec3 {
    Vec3 coords; // Assuming Vec3 has x, y, z floats

    // Constructor
    GlobalVertexVec3(const Vec3& v) : coords(v) {}

    // EQUALITY OPERATOR (Required for unordered_map collision handling)
    // Returns true if two vertices snap to the same integer grid point.
    bool operator==(const GlobalVertexVec3& other) const {
        const float inv_tol = 100000.0f; // Precision: 0.00001

        // Quantize inputs to long integers
        long long x1 = static_cast<long long>(std::llround(coords.x * inv_tol));
        long long y1 = static_cast<long long>(std::llround(coords.y * inv_tol));
        long long z1 = static_cast<long long>(std::llround(coords.z * inv_tol));

        long long x2 = static_cast<long long>(std::llround(other.coords.x * inv_tol));
        long long y2 = static_cast<long long>(std::llround(other.coords.y * inv_tol));
        long long z2 = static_cast<long long>(std::llround(other.coords.z * inv_tol));

        return x1 == x2 && y1 == y2 && z1 == z2;
    }
};

// 2. The Hash Functor
struct GlobalVertexVec3Hash {
    // HASH FUNCTION (Required for unordered_map bucket finding)
    // Must use the EXACT SAME quantization logic as operator==
    std::size_t operator()(const GlobalVertexVec3& v) const {
        const float inv_tol = 100000.0f; // MUST match operator==

        long long x = static_cast<long long>(std::llround(v.coords.x * inv_tol));
        long long y = static_cast<long long>(std::llround(v.coords.y * inv_tol));
        long long z = static_cast<long long>(std::llround(v.coords.z * inv_tol));

        // Spatial Hashing using large primes
        return static_cast<std::size_t>(
            (x * 73856093) ^ (y * 19349663) ^ (z * 83492791)
            );
    }
};

class Generator {

public:

    ~Generator() {};

    Generator() : bounds({ 0.0f, 10.0f, 0.0f, 10.0f, 0.0f, 10.0f }), blockDim({ 64, 64, 64 }), threshold(0.0f), isoLevel(1.0f) {};

    Generator(
        std::vector<Vec3>& seeds,
        const std::array<float, 6>& bounds,
        const std::array<int, 3>& blockDim,
        const float threshold = 0.0f,
        const float isoLevel = 0.1f,
        const int foam = 0
    );

    void populate_grids(const IContainer& con, bool pad = false);

    //void _pad_bounds();

    void smooth_scalar_field();
    //void estimate_scalar_field();

    void marching_cubes();

    void draw();

    void draw_edges();

    void draw_tortuosity_path();

    void draw_pore_network();

    void export_stl(std::string fileName);

    void export_nrrd(const std::string fileName, float voxelSize, std::array<float, 6> blockSize);

    void export_mhd(std::filesystem::path& path, float voxelSize, std::array<float, 6> blockBounds);

    void render_properties();

    Metrics get_metrics();

    bool estimate_tortuosity();

    bool estimate_connectivity_network();

    void set_bounds(std::array<float, 6>& newBounds);

    void set_thickness(const float newThickness);

    void set_seeds(const std::vector<Vec3>& newSeeds);

    void set_resolution(const std::array<int, 3>& newResolution);

    void set_openess(const float newOpeness);

    void estimate_local_thickness(float voxelSize, std::array<float, 6>& blockBounds, bool separation = false);

    std::vector<uint8_t> get_image_field(float voxelSize, std::array<float, 6>& blockBounds, bool inverse = false);

    void estimate_local_separation();

    void render_metrics();

    std::array<float, 6> get_bounds() const {
        return bounds;
    }

    std::array<float, 4> color = { 0.5f, 0.5f, 0.5f, 1.0f };

    bool hidden = false;

    std::string name = "";

    void* container = nullptr;

    void* generator = nullptr;

    int foam = false;

private:

    //@brief function to interpolate the vertex position on the edge based on the scalar values at the endpoints and the isoLevel
    Vec3 interpolate_vertex(
        float isoLevel,
        const Vec3& p1,
        const Vec3& p2,
        float val1, float val2);

    int find_vertex_index(int x, int y, int z);

    int find_cube_index(float isoLevel, float values[8]);

    Vec3 get_position(int x, int y, int z);

    int check_vertex(Vec3& vec);

    void add_edge(int idx1, int idx2, int idx3);

    void _setup_mesh();

    void _setup_edges();

    // parameters
    std::vector<Vec3> seeds;
    std::array<float, 6> bounds;
    std::array<float, 6> paddedBounds{};
    std::array<int, 3> blockDim;
    std::array<int, 3> workingBlockDim;

    float threshold{ 0.6f };
    float isoLevel{ 0.01f };
    float stepX{ 0.0 }, stepY{ 0.0 }, stepZ{ 0.0 };

    // graph for connectivity between seeds
    std::unique_ptr<Graph> graph;

    // kdtree for nearest neighbor search
    std::unique_ptr<Kdtree> kdtree;

    // scalar field
    std::vector<float> scalarField;

    // hold triangles
    std::vector<Triangle> triangles;

    // these are for rendering
    std::vector<float> vertices;
    std::vector<float> normals;

    // temporary storage for normals, since we need to compute the normal for each triangle and then average them for each vertex, we need to store them in a temporary vector before we push them to the normals vector for rendering
    std::vector<Vec3> tempNormals;

    // these are to keep track of the global vertex indices and positions, 
    // since the marching cubes will create vertices that are shared between cells, we need to keep track of them to avoid duplicates
    std::unordered_map<GlobalVertexVec3, int, GlobalVertexVec3Hash> globalVertexMap;
    int globalIndex{ 0 };

    // we need to also have a map for edge indices, since the marching cubes will create edges that are shared between cells, we need to keep track of them to avoid duplicates
    std::set<std::pair<unsigned int, unsigned int>> edgeSet;

    // these are for rendering mesh edges
    std::vector<unsigned int> edgeIndices;

    // this is for rendering the EBO
    std::vector<unsigned int> indices;

    // store also a list of edges to draw the tortuosity path
    std::vector<unsigned int> tortuosityPathEdges;
    std::vector<float> tortuosityPathVertices;
    std::unique_ptr<PoreNetwork> tortuosityPathModel;

    // connectivity network model for rendering
    std::unique_ptr<PoreNetwork> connectivityPathModel;
    std::vector<float> connectivityPathVertices;
    std::vector<unsigned int> connectivityPathEdges;

    // metrics
    float surfaceArea{ 0.0f };
    float porosity{ 0.0f };
    float volumeFraction{ 0.0f };
    float volume{ 0.0f };
    float tortuosity{ 0.0f };
    float localThickness{ 0.0f };
    float localThicknessStd{ 0.0f };
    float localSeparation{ 0.0f };
    float localSeparationStd{ 0.0f };

    // opengl
    unsigned int VAO{ 0 }, VBO{ 0 }, EBO{ 0 };
    unsigned int edgeVAO{ 0 }, edgeVBO{ 0 }, edgeEBO{ 0 }, normalsVBO{ 0 };
    unsigned int tortuosityPathVAO{ 0 }, tortuosityPathVBO{ 0 };
};

float smin(float a, float b, float k);

#endif