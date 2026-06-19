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
#include <string>
#include <atomic>
#include <queue>

#include "SeedGenerator/Container.h"
#include <SeedGenerator/SeedGenerator.h>
#include <SeedGenerator/DistanceCalculator.h>
#include <SeedGenerator/RadiusCalculator.h>
#include "Math/Vec.h"
#include "Utils/Utils.h"
#include "Logger/Logger.h"

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

struct ScaffoldParameters {
	// 0: Uniform, 1: Varied
	int thicknessOption = -1;

	// Thickness
	float uniformThickness = -1.0f;
	float startThickness = -1.0f;
	float endThickness = -1.0f;
	int distFunction = -1;
	int radFunction = -1;

	// Distance Related Vectors (Flattened)
	float planeOriginX = 0.0f, planeOriginY = 0.0f, planeOriginZ = 0.0f;
	float planeNormalX = 0.0f, planeNormalY = 0.0f, planeNormalZ = 1.0f;
	float pointX = 0.0f, pointY = 0.0f, pointZ = 0.0f;

	// Generator (0: Random, 1: Poisson)
	int generatorType = -1;
	int seedNr = -1;
	float minRadius = -1.0f;
	float maxRadius = -1.0f;

	// Global
	float openess = -1.0f;
	float stretchX = 1.0f, stretchY = 1.0f, stretchZ = 1.0f;
	float anisotropyAngle = 0.0f;
	float dirX = 1.0f, dirY = 0.0f, dirZ = 0.0f;
};

struct ScaffoldMetrics {
	float porosity = -1.0f;
	float volume = -1.0f;
	float totalSurface = -1.0f;
	float surfToVol = -1.0f;
	float connectivityDensity = -1.0f;
	float thickness = -1.0f;
	float thicknessStd = -1.0f;
	float separation = -1.0f;
	float separationStd = -1.0f;
	float trNumber = -1.0f;
	float anisotropyDeg = -1.0f;
	float tortuosity = -1.0f;
};

class GeneratorLewiner {
public:

	GeneratorLewiner(bool render = true) : renderMode(render) {
		
		if (renderMode) {
			_setup_mesh();

			_setup_edges();
		}
	};

	GeneratorLewiner(
		const std::vector<Vec3>& seeds,
		const std::array<float, 6>& bounds,
		const std::array<int, 3>& dims,
		Logger* uiLogger,
		const float threshold = 0.0f,
		const float isoLevel = 0.1f,
		const int foam = 0,
		const bool renderMode = true
	);

	GeneratorLewiner(
		const std::string fileName,
		Logger* uiLogger,
		const bool renderMode = true
	);

	~GeneratorLewiner() {};

	// functions
public:
	void compute_scalar_field(const IContainer& con);

	void marching_cubes();

	void export_stl(std::string fileName);

	void estimate_metrics(const IContainer& container);

	void set_resolution(const std::array<int, 3>& newResolution);

	void set_bounds(const std::array<float, 6>& newBounds);

	void set_seeds(const std::vector<Vec3>& newSeeds);

	void set_stretch(float newStretchX, float newStretchY, float newStretchZ);

	void set_thickness(float newThickness);

	std::array<float, 6> get_bounds() const;

	Aabb get_aabb() const;

	void estimate_local_thickness(float voxelSize, std::array<float, 6>& blockBounds, bool separation = false);

	bool estimate_tortuosity(float voxelSize);

	//void estimate_anisotropy(int daDirectionNr, int daMinsteps, int daMaxsteps, float vcLimit, int mode = 0);

	void estimate_anisotropy(int daDirectionNr = 2000, int linesPerDirection = 10000, int mode = 3, ROI* roi = nullptr);

	void estimate_trabecular_number(
		int formula = 0, int daDirectionNr = 2000, int linesPerDirection = 10000, ROI* roi = nullptr);

	void estimate_connectivity_density();

	void estimate_connectivity_network();

	void export_nrrd(const std::string fileName, float voxelSize, std::array<float, 6> blockSize);

	void export_mhd(std::filesystem::path& path, float voxelSize, std::array<float, 6> blockBounds);

	void apply_taubin_smooth(int iter, float lambda, float mu);

	void export_metrics(std::string fileName);
	
	void read_metrics(const std::string fileName);

	void export_parameters(std::string fileName);

	void read_parameters(const std::string fileName);

	void set_thickness_functions(
		std::shared_ptr<const SDF> sdf,
		std::shared_ptr<const RadiusFunction> radFunc,
		float tMin, float tMax, float distance
	);

	void set_options_from_factory(int distOption, int distFunc, int thicknessOption, float voxSize);

	void set_distance_plane_options(Vec3 center, Vec3 normal);

	void set_distance_point_options(Vec3 point);

	// create sub scaffold from ROI
	std::unique_ptr<GeneratorLewiner> extract_from_ROI(ROI* roi);

	void set_logger(Logger* newLogger) { logger = newLogger; };

	// saving / loading
	void export_scaf(const std::string& fileName);
	bool load_scaf(
		const std::string& fileName,
		std::vector<std::shared_ptr<IContainer>>& containerList,
		std::vector<std::shared_ptr<InterfaceSeedGenerator>>& generatorList
	);

	//--------------------------------
	// rendering
	void draw();
	void draw_edges();
	void render_properties(bool& updateScaffold);
	void render_metrics();
	void draw_tortuosity_path();
	void apply_scale();

private:
	// Per-voxel classification produced by classify_container_narrow_band(),
	// used to skip the expensive exact container-distance evaluation
	// (BVH nearest-point query + raycast inside/outside test) for voxels
	// that a cheap coarse pre-pass already proves are safely far from the
	// container surface. See classify_container_narrow_band() for the
	// derivation of why this is exact, not an approximation.
	enum class NarrowBandClass : uint8_t {
		NeedsExact = 0,	// within 'margin' of the surface (or unproven) - must call the real SDF
		Outside    = 1,	// proven > margin outside the container
		Inside     = 2	// proven > margin inside the container
	};

	std::vector<NarrowBandClass> classify_container_narrow_band(const IContainer& con, float margin);

	void smooth_scalar_field();

	void seal_grid_boundaries();

	void remove_isolated_islands();

	void compute_intersection_points();

	void process_cube(int i, int j, int k, const float cube[8], int lut_entry);

	void update_steps();

	size_t find_vertex_index(int x, int y, int z);

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

	size_t get_x_vert(int i, int j, int k);

	size_t get_y_vert(int i, int j, int k);

	size_t get_z_vert(int i, int j, int k);

	size_t add_c_vertex(const int i, const int j, const int k);

	bool test_face(signed char face, const float _cube[8]);

	bool modified_test_interior(signed char s, int _case, int _config, const float _cube[8]);

	int interior_ambiguity(int amb_face, int s, const float _cube[8]);

	int interior_ambiguity_verification(int edge, const float _cube[8]);

	bool interior_test_case13(const float _cube[8]);

	bool interior_test_case13_2(float isovalue, const float _cube[8], int& tunnelOrientation);

	void validate_topology();

	void build_topology();

	void _update_bounding_box();

	std::vector<uint8_t> get_image_field(
		float voxelSize, std::array<float, 6>& blockBounds, bool inverse = false, uint8_t solidValue = 1);

	//-------------------------------
	// rendering
	void _setup_mesh();
	void _setup_edges();
	void _update_render();
	void add_edge(int idx1, int idx2, int idx3);

	// members
public:
	std::weak_ptr<IContainer> container;
	std::weak_ptr<InterfaceSeedGenerator> generator;
	std::shared_ptr<const SDF> thicknessSDF = nullptr;
	std::shared_ptr<const RadiusFunction> thicknessFunction = nullptr;
	bool foam = false;
	std::string name = "";
	std::array<float, 4> color = { 0.5f, 0.5f, 0.5f, 1.0f };
	bool hidden = false;
	bool hiddenTortuosityPath = false;
	bool hiddenEllipsoid = false;
	bool hiddenNetworkPath = false;
	bool isROI = false;

	float volume{ 0.0f };
	float domainVolume{ 0.0f };
	float porosity{ 0.0f };
	float surfaceArea{ 0.0f };
	float surfaceToVolume{ 0.0f };
	float localThickness{ 0.0f };
	float localThicknessStd{ 0.0f };
	float localSeparation{ 0.0f };
	float localSeparationStd{ 0.0f };
	float tortuosity{ 0.0f };
	float anisotropyDegree{ 0.0f };
	float trabecularNr{ 0.0f };
	float connectivityDensity{ 0.0f };

	Vec3 translateVec{ 0.0f, 0.0f, 0.0f };
	Vec3 scaleVec{ 0.0f, 0.0f, 0.0f };
	std::unique_ptr<PoreNetwork> tortuosityPathModel;

	float stretchX{ 1.0f };
	float stretchY{ 1.0f };
	float stretchZ{ 1.0f };

	Vec3 anisotropyVec{ 1.0f, 0.0f, 0.0f };
	float anisotropyAngle{ 0.0f };
	std::unique_ptr<Ellipsoid> ellipsoidModel;

	// loaded from csv
	bool isLoadedFromFile = false;

private:
	std::vector<Vec3> seeds;
	float stepX{ 0.0f }, stepY{ 0.0f }, stepZ{ 0.0f };
	std::vector<float> scalarField;
	std::array<int, 3> blockDims = { 100, 100, 100 };
	float voxelSize{ 0.05f };
	std::array<float, 6> bounds = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
	
	int selectedThicknessOption = 0;
	int selectedDist = 0;
	int selectedFunc = 0;
	float startThickness = 0.3f;
	float endThickness = 1.0f;
	float transitionDistance = 10.0f;
	float isoLevel{ 0.1f };
	Vec3 distancePlaneNormal;
	Vec3 distancePlaneCenter;
	Vec3 distancePoint;

	std::vector<size_t> x_verts;
	std::vector<size_t> y_verts;
	std::vector<size_t> z_verts;

	std::vector<LTriangle> meshTriangles;
	std::vector<LVertex> meshVertices;

	std::atomic<size_t> vertexCount{ 0 };
	std::atomic<size_t> triangleCount{ 0 };

	// scaffold properties
	float threshold = { 0.5f };

	Aabb aabb;

	Logger* logger = nullptr;

	// ----------------------------------------------------
	// these are for opengl
	std::vector<float> vertices;
	std::vector<float> normals;
	std::set<std::pair<unsigned int, unsigned int>> edgeSet;
	// these are for rendering mesh edges
	std::vector<unsigned int> edgeIndices;
	// this is for rendering the EBO
	std::vector<unsigned int> indices;
	std::vector<std::vector<unsigned int>> adjacency;

	// store also a list of edges to draw the tortuosity path
	std::vector<unsigned int> tortuosityPathEdges;
	std::vector<float> tortuosityPathVertices;

	unsigned int VAO{ 0 }, VBO{ 0 }, EBO{ 0 };
	unsigned int edgeVAO{ 0 }, edgeVBO{ 0 }, edgeEBO{ 0 }, normalsVBO{ 0 };
	unsigned int tortuosityPathVAO{ 0 }, tortuosityPathVBO{ 0 };

	bool renderMode = true;
	// ---------------------------------------------------
	// version handling
	uint32_t meshVersion = 0;
	uint32_t thicknessVersion = 0;
	uint32_t separationVersion = 0;
	uint32_t trabecularNrVersion = 0;
	uint32_t connectivityVersion = 0;
	uint32_t tortuosityVersion = 0;
	uint32_t anisotropyVersion = 0;
};

class ScaffoldFactory {

public:
	void launch();

	void gui_draw(
		Logger* logger,
		const char* popupName, bool& showPopup,
		SelectedObject* selectedPanelObj, void*& selectedSceneObj,
		std::vector<std::unique_ptr<GeneratorLewiner>>& scaffoldList,
		std::vector<std::shared_ptr<IContainer>>& containers,
		std::vector<std::shared_ptr<InterfaceSeedGenerator>>& generators);

private:
	std::weak_ptr<IContainer> selectedCon;
	std::weak_ptr<InterfaceSeedGenerator> selectedGen;
	char buffer[256] = "";
	float thickness = { 0.3f };
	float openess = { 0.5f };
	float stretchX = { 1.0f };
	float stretchY = { 1.0f };
	float stretchZ = { 1.0f };
	Vec3 anisotropyVec = { 1.0f, 0.0f, 0.0f };
	float anisotropyAngle = { 0.0f };
	int foam = 0;
	std::array<int, 3> resolution = { 100, 100, 100 };
	float voxelSize{ 0.05f };
	uint32_t lastUsedContainerVersion = 0;
	uint32_t lastUsedGeneratorVersion = 0;

	int selectedThicknessOption = false;
	Vec3 distancePlaneNormal;
	Vec3 distancePlaneCenter;
	Vec3 distancePoint;

	int selectedDist = 0;
	int selectedFunc = 0;
	std::shared_ptr<const SDF> thicknessSDF;
	std::shared_ptr<const RadiusFunction> thicknessRadFunc;

	float startThickness = 0.3f;
	float endThickness = 1.0f;
	float transitionDistance = 20.0f;
};

#endif // ! "GENERATOR_LEWINER_H"
