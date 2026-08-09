/**/
//________________________________________________

#ifndef GENERATOR_LEWINER_H
#define GENERATOR_LEWINER_H

#include <vector>
#include <array>
#include <string>
#include <atomic>
#include <queue>
#include <functional>

#include "SeedGenerator/Container.h"
#include <SeedGenerator/SeedGenerator.h>
#include <SeedGenerator/DistanceCalculator.h>
#include <SeedGenerator/RadiusCalculator.h>
#include <ScaffoldGenerator/GenerationTask.h>
#include "Math/Vec.h"
#include "Math/QuadricSimplifier.h"
#include "Utils/Utils.h"
#include "Anisotropy.h"
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
	float spread = 0.5f;
	float stretchX = 1.0f, stretchY = 1.0f, stretchZ = 1.0f;
	float anisotropyAngle = 0.0f;
	float dirX = 1.0f, dirY = 0.0f, dirZ = 0.0f;
	float backgroundWeight = 0.1f;

	// Junction smoothing (smin fillet of the rod/plate fields)
	int smoothJunctions = 0;      // 0/1 flag
	float smoothK = 0.01f;

	// Edge rounding (smin softening of the raw distance order statistics)
	int roundEdges = 0;           // 0/1 flag
	float edgeK = 0.01f;

	// Boundary frame: close the Voronoi cells near the container wall into full
	// walls (a honeycomb rim that ties the cut struts together).
	int frameBoundary = 0;        // 0/1 flag
	float frameDepth = 1.0f;      // band depth as a MULTIPLE of wall thickness (isoLevel)

	// Container edge frame: solid beams along the container's own edges (box: 12
	// edges, cylinder: 2 rims), forming a rigid outer cage for test specimens.
	int frameContainerEdges = 0;  // 0/1 flag
	float frameBeam = 0.3f;       // mm; radius of the edge beams

	// Calibration state (checkboxes + targets), so a saved setup restores which
	// metrics were being driven and to what values.
	int calibrateThickness = 1;   // 0/1 flag
	float targetThickness = 0.3f;
	int calibratePorosity = 1;    // 0/1 flag
	float targetPorosity = 80.0f;
	int calibrateDA = 1;          // 0/1 flag
	float targetDa = 1.2f;
	int targetFormulaIdx = 0;
	float calibrationStep = 0.02f;
	float calibrationTol = 0.001f;
	int calibrationIter = 10;

	// Varied-thickness transition distance (only meaningful when thicknessOption == 1)
	float transitionDistance = 10.0f;
};

struct ScaffoldMetrics {
	float porosity = -1.0f;
	float porosityMesh = -1.0f;
	float volume = -1.0f;
	float totalSurface = -1.0f;
	float surfToVol = -1.0f;
	float surfToTotalVol = -1.0f;
	float connectivityDensity = -1.0f;
	float smi = -1.0f;
	float thickness = -1.0f;
	float thicknessStd = -1.0f;
	float separation = -1.0f;
	float separationStd = -1.0f;
	float trNumber = -1.0f;
	float anisotropyDeg = -1.0f;
	float tortuosity = -1.0f;
};

struct VoxelValue {
	float strutVal = 0.0f;
	float wallVal = 0.0f;
	Vec3 strutGrad = {0.0f, 0.0f, 0.0f};
	Vec3 wallGrad = {0.0f, 0.0f, 0.0f};
	// Fourth-nearest contrast (d4 - d1) and its gradient. Only used by the
	// optional edge-rounding (k_edge) softening, which smins strutVal against
	// edgeVal to round the Voronoi VERTICES (where d3 = d4) exactly the way
	// wallVal-vs-strutVal rounds the EDGES (where d2 = d3). When edge rounding is
	// off this is inert. edgeVal >= strutVal >= wallVal always (d2<=d3<=d4).
	float edgeVal = 0.0f;
	Vec3 edgeGrad = {0.0f, 0.0f, 0.0f};
	float faceP = 0.0f;
	// |distance| from the thickness SDF, needed to RE-derive the per-voxel iso
	// level in assemble_field for VARIED thickness (localIsoLevel itself depends
	// on the thickness knob being calibrated, so it must not be cached). For
	// uniform thickness this is unused (localIsoLevel = isoLevel).
	float rawDist = 0.0f;
};

struct SmoothDist {
	float val;
	Vec3 grad;
};

// enum for the knob recache for DA because it needs stretch computations, recompute for others
enum class Recompute {Recache, Recompute};

struct CalibrationStage {

	std::string name = "";
	Recompute cost;
	// set stage
	std::function<void(float)> set_knob;
	// current value
	std::function<float()> get_knob;
	// return metric re-estimate mesh if needed 
	std::function<float()> measure;
	float target, tol, lo, hi;
	int maxIter=5;

	bool ascending;

	// Stamp this metric's version to meshVersion. The inner stages measure on an
	// earlier meshVersion than the final mesh, so their versions read stale even
	// though the VALUE is correct for the converged config (re-mesh is
	// deterministic). Applied once, after the final mesh is built, for every
	// stage that converged. Null for metrics with no version (voxel porosity).
	std::function<void()> stampVersion;
	bool converged = false;
};

// =============================================================================
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
	// Both return false when their preconditions are not met (too few seeds /
	// empty scalar field), so callers can skip the rest of the pipeline
	// instead of running marching cubes or metrics on garbage.
	bool compute_cached_field_values(const IContainer& con);

	bool assemble_field();

	// Full regeneration: the cache pass (expensive kdtree/wall-strut/gradient,
	// invariant to the calibration knobs) followed by one assemble pass. Kept as
	// the single entry point every non-calibration caller uses. Calibration calls
	// compute_cached_field_values once and then assemble_field() per iteration.
	bool compute_scalar_field(const IContainer& con);

	bool marching_cubes(bool supress = false);

	void export_stl(std::string fileName);

	void estimate_metrics(const IContainer& container);

	void set_resolution(const std::array<int, 3>& newResolution);

	void set_bounds(const std::array<float, 6>& newBounds);

	void set_seeds(const std::vector<Vec3>& newSeeds);

	void set_stretch(float newStretchX, float newStretchY, float newStretchZ);

	void set_thickness(float newThickness);

	std::array<float, 6> get_bounds() const;

	Aabb get_aabb() const;

	// Calibrated "background" values: the target Tb.Th / porosity / SMI in the GUI
	// stay fixed (they are goals); calibration moves these internal knobs to hit
	// them. Expose them so the GUI can show the achieved iso-level / openness.
	float get_iso_level() const { return isoLevel; }
	float get_openness() const { return threshold; }
	float get_spread() const { return spread; }

	void estimate_local_thickness(float voxelSize, std::array<float, 6>& blockBounds, bool separation = false, bool supress = false);

	// Validate the local-thickness estimator on a junction-free geometry: a
	// solid slab of known physical thickness. Returns {mean, std} in mm. A slab
	// has no junctions, so the estimator should return the true thickness (up to
	// the ~1-voxel distance-transform bias), isolating measurement accuracy from
	// the junction-bleeding that elevates Tb.Th on foam/lattice structures.
	std::array<float, 2> run_slab_phantom(float slabThickness, float voxelSize, float domainSize = 0.0f);

	// Analytic phantoms with known ground truth, used to verify the metric
	// estimators independently of the generator. 'feature' is the characteristic
	// size the estimators should recover:
	//   SLAB     plate of thickness   'feature'  -> Tb.Th = feature, SMI ~ 0
	//   CYLINDER rod of diameter      'feature'  -> Tb.Th = feature, SMI ~ 3
	//   SPHERE   ball of diameter     'feature'  -> Tb.Th = feature, SMI ~ 4
	//   TORUS    tube of diameter     'feature'  -> Tb.Th = feature, Conn = 1
	// build_phantom_field only lays down the scalar field (so the caller may run
	// any estimator or marching_cubes); the generator is left ready to measure.
	enum PhantomShape { PHANTOM_SLAB = 0, PHANTOM_CYLINDER = 1, PHANTOM_SPHERE = 2, PHANTOM_TORUS = 3 };
	void build_phantom_field(int shape, float feature, float voxelSize, float domainSize = 0.0f);

	bool solve_calibration(
		std::vector<CalibrationStage>& stageSet,
		int lvl,
		const IContainer& con
	);

	// Assemble the ENABLED calibration knobs into `stages`, ordered outer->inner:
	// DA -> SMI -> Thickness -> Porosity. Outermost = most expensive metric and
	// most coupling (fewest evaluations); innermost = cheapest (evaluated most).
	// voxelSize is the MEASUREMENT voxel (metrics are measured there, decoupled
	// from the generation grid). Any subset reproduces the old two/three-knob.
	void build_calibration_stages(float voxelSize);

	// Stamp every converged calibration stage's metric version to the CURRENT
	// meshVersion. Call once AFTER the final mesh is built: the inner stages were
	// measured on an earlier meshVersion, so their versions read stale even though
	// their values already match the converged config (re-mesh is deterministic).
	void stamp_calibrated_versions();

	// Stage factories: each binds a knob (set/get) and its metric (measure) to
	// THIS generator. DA re-caches (stretch changes the warped seeds); the rest
	// only re-assemble.
	CalibrationStage make_da_stage(float target, float voxelSize, float tol, int maxIter, int mode);
	// CalibrationStage make_smi_stage(float target, float tol, int maxIter);   // SMI calibration disabled
	CalibrationStage make_thickness_stage(float target, float voxelSize, float tol, int maxIter);
	CalibrationStage make_porosity_stage(float target, float tol, int maxIter);

	bool secant_1d(
		const std::function<float(float)>& eval,
	 	float guess, float target, float tol, float lo, float hi, int maxIter,
		bool ascending, float& xOut, 
		const std::string& label = "", int depth = 0
	);

	// Auto-calibrate the iso-level (c_glob) so the MEASURED local thickness at
	// the given voxel size matches targetThickness (mm). A secant search over the
	// monotonic map, regenerating and re-measuring each step, so it needs no
	// fitted model and is agnostic to openness, seed spacing, voxel size and
	// container shape. UNIFORM thickness: tunes isoLevel. VARIED thickness
	// (thicknessFunction set): scales the whole [startThickness, endThickness]
	// range by one factor, preserving the grading, so the measured MEAN Tb.Th
	// hits the target. Leaves the FIELD assembled at the calibrated value but does
	// NOT build the mesh (the secant skips marching cubes for speed); standalone
	// callers run marching_cubes() once after it returns. two_knob/three_knob
	// build the mesh themselves. Builds/reuses the per-voxel cache automatically.
	// onProgress, if given, is called with a fraction in [0,1] each iteration.
	// Returns true if the generator holds a valid calibrated field.
	bool calibrate_thickness(const IContainer& con, float targetThickness,
		float voxelSize, float tol = 0.002f, int maxIter = 8,
		const std::function<void(float)>& onProgress = nullptr);

	// calibrate openness to achieve target porosity with fixed thickness
	bool calibrate_openness(const IContainer& con, float targetPorosity, 
		float voxelSize, float tol=0.001, int maxIter = 10,
		const std::function<void(float)>& onProgress = nullptr);

		
	// two knob calibration of thickness and porosity. This adjusts the isoLevel to reach targeted Tb.Th, then freezes it and adjusts openness to achieve porosity
	bool two_knob_calibration(
		const IContainer& con, 
		float targetThickness, 
		float targetPorosity,
		float voxelSize,
		float tol=0.001, int maxIter = 10,
		const std::function<void(float)>& onProgress = nullptr
	);

	// three knob calibration of thickness, porosity and SMI. This adjusts the isoLevel to reach targeted Tb.Th, then freezes it and creates a 2x2 Newton scheme to adjust openness and spread to achieve a targeted porosity & SMI
	bool three_knob_calibration(
		const IContainer& con, 
		float targetThickness, 
		float targetPorosity,
		float targetSMI,
		float voxelSize,
		float tol=0.001, int maxIter = 10,
		const std::function<void(float)>& onProgress = nullptr
	);

	bool estimate_tortuosity(float voxelSize);

	//void estimate_anisotropy(int daDirectionNr, int daMinsteps, int daMaxsteps, float vcLimit, int mode = 0);

	// voxelSize: the target (e.g. microCT) voxel size the structure is resampled
	// to before the MIL rays are cast, so this metric shares the same sampling
	// resolution as local thickness/separation rather than the generation grid.
	void estimate_anisotropy(float voxelSize, int daDirectionNr = 2000, int linesPerDirection = 10000, int mode = 3, ROI* roi = nullptr);

	void estimate_trabecular_number(
		float voxelSize, int formula = 0, int daDirectionNr = 2000, int linesPerDirection = 10000, ROI* roi = nullptr);

	void estimate_connectivity_density();

	// Connectivity density from the Euler characteristic of the resampled solid
	// image (BoneJ-style), robust to multiple components and internal cavities
	// unlike the surface-genus method. Samples the same image substrate as the
	// other voxel metrics. connectivity: 6 (implemented) or 26 (future).
	// Returns Conn.D in 1/mm^3.
	float estimate_connectivity_density_voxel(float voxelSize, int connectivity = 6);

	// Structure Model Index (Hildebrand & Ruegsegger 1997): the rod-vs-plate
	// descriptor. SMI = 6*(BS'*BV)/BS^2, where BS' = dBS/dr is obtained by
	// dilating the surface a small step along the vertex normals.
	//   ~0 -> ideal plate, ~3 -> ideal rod, ~4 -> sphere.
	// This is the only metric that isolates the rod/plate axis that 'spread'
	// controls; the others (BV/TV, Conn.D, Tb.N) confound it with openness and
	// seed density. Mesh-based, so it needs marching_cubes() to have run.
	// dilation: offset in mm; <= 0 picks a scale-relative default.
	// Returns the SMI.
	float estimate_smi(float dilation = 0.0f);

	void estimate_connectivity_network();

	void export_nrrd(const std::string fileName, float voxelSize, std::array<float, 6> blockSize);

	void export_mhd(std::filesystem::path& path, float voxelSize, std::array<float, 6> blockBounds);

	void smooth_scalar_field_taubin(int iterations, float lambda, float mu);

	void apply_taubin_smooth(int iter, float lambda, float mu);

	void apply_mesh_simplification(const mesh_simplify::Options& options);

	void export_metrics(std::string fileName);
	
	void read_metrics(const std::string fileName);

	void export_parameters(std::string fileName);

	void read_parameters(const std::string fileName);

	void set_thickness_functions(
		std::shared_ptr<const SDF> sdf,
		std::shared_ptr<const RadiusFunction> radFunc,
		float tMin, float tMax, float distance
	);

	void set_options_from_factory(int distOption, int distFunc, int thicknessOption, float voxSize, float measureVoxel);

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
		std::vector<std::shared_ptr<InterfaceSeedGenerator>>& generatorList,
		std::vector<std::shared_ptr<AnisotropySource>>& globalSources,
		std::atomic<int>* stage = nullptr
	);

	//--------------------------------
	// rendering
	void draw();
	void draw_edges();
	void render_properties(
		bool& updateScaffold, GenerationTask* task,
		std::vector<std::shared_ptr<AnisotropySource>>& globalSources
	);
	void render_metrics();
	void draw_tortuosity_path();
	void apply_scale();
	void update_buffers();
	void update_render();

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

	// Field value that marks a voxel as "safely air" near the surface: a few
	// voxels above the iso-level. Expressed in grid spacings (not a fixed mm
	// offset) so the smoothing/early-out skip band tracks the model scale and
	// resolution. Shared by compute_scalar_field and both smoothers so their
	// "set to air" and "skip if air" tests stay identical.
	float air_skip_level() const;

	size_t find_vertex_index(int x, int y, int z);

	Vec3 get_position(int x, int y, int z);

	SmoothDist smin_gradient(	
		float a, float b, 
		const Vec3& gradA, const Vec3& gradB, float k = 0.5f);

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

	void validate_topology(bool supress = false);

	void build_topology();

	void _update_bounding_box();

	std::vector<uint8_t> get_image_field(
		float voxelSize, std::array<float, 6>& blockBounds, bool inverse = false, uint8_t solidValue = 1);

	void main_properties();
	void smooth_properties();
	void thickness_properties();
	void anisotropy_properties(
		std::vector<std::shared_ptr<AnisotropySource>>& globalSources);

	//-------------------------------
	// rendering
	void _setup_mesh();
	void _setup_edges();
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
	std::vector<CalibrationStage> stages;

	bool calibrateDA = true;
	float targetDa = 1.2f;
	int targetFormulaIdx = 0;
	bool calibrateThickness = true;
	float targetThickness = 0.3f;
	bool calibratePorosity = true;
	float targetPorosity = 80.0f;
	// bool calibrateSmi = false;   // SMI calibration disabled (SMI is reported, not targeted)
	// float targetSmi = 0.5f;
	float calibrationStep = 0.02f;
	float calibrationTol = 0.001f;
	int calibrationIter = 10;

	float measurementVoxelSize = 0.02f;
	float volume{ 0.0f };
	float domainVolume{ 0.0f };
	// primary porosity: voxel-based, computed in compute_scalar_field on the
	// generation grid (numerator and denominator share the discretization)
	float porosity{ 0.0f };
	// legacy porosity: mesh volume vs analytic container volume; systematically
	// overestimates because the mesh never reaches the container surface
	float porosityMesh{ 0.0f };
	float surfaceArea{ 0.0f };
	// BS/BV: bone surface per BONE volume. Tied to thickness: ~2/Tb.Th for
	// ideal plates, ~4/Tb.Th for ideal rods, so it is a useful cross-check.
	float surfaceToVolume{ 0.0f };
	// BS/TV: bone surface per TOTAL (sample) volume. This is the quantity
	// reported as "BS/TV" in the uCT literature (Parfitt/ASBMR nomenclature).
	// BS/TV = (BS/BV) * (BV/TV), so it is ~an order of magnitude smaller.
	float surfaceToTotalVolume{ 0.0f };
	float localThickness{ 0.0f };
	float localThicknessStd{ 0.0f };
	float localSeparation{ 0.0f };
	float localSeparationStd{ 0.0f };
	float tortuosity{ 0.0f };
	float anisotropyDegree{ 0.0f };
	float trabecularNr{ 0.0f };
	float connectivityDensity{ 0.0f };      
	float smi{ 0.0f };                // Structure Model Index: 0=plate, 3=rod

	Vec3 translateVec{ 0.0f, 0.0f, 0.0f };
	Vec3 scaleVec{ 0.0f, 0.0f, 0.0f };
	std::unique_ptr<PoreNetwork> tortuosityPathModel;

	float spread = 0.5f;

	// Optional junction smoothing (smooth-min fillet of the rod/plate fields).
	// Off by default (the validated linear blend); when on, smoothK sets the
	// fillet radius (~0.01 works with typical thicknesses). Tb.Th inflation from
	// the fattened nodes is absorbed by calibrate_thickness.
	bool smoothJunctions = false;
	float smoothK = 0.01f;

	// Optional edge rounding (smooth-min softening of the raw distance order
	// statistics, applied BEFORE the tau blend -- distinct from the junction
	// fillet above, which acts AFTER it). Rounds the intrinsic Voronoi cell
	// edges/vertices, so it works at every openness (including tau=1 bare
	// lattices, where the junction fillet degenerates to a uniform offset).
	// Off by default; edgeK=0 reduces EXACTLY to the linear blend. Tb.Th
	// inflation from the fattened edges is absorbed by calibrate_thickness.
	bool roundEdges = false;
	float edgeK = 0.01f;

	// Boundary frame. When on, the openness is ramped to zero within frameDepth of
	// the container wall, so the outermost cells close into full Voronoi walls and
	// form a connected honeycomb rim around the scaffold (a "frame" that ties the
	// otherwise cut boundary struts together, useful for gripping test specimens).
	bool frameBoundary = false;
	float frameDepth = 1.0f;      // band depth as a multiple of wall thickness (isoLevel)

	// Container edge frame. When on, solid beams of radius frameBeam are unioned
	// along the container's own edges (a box's 12 edges, a cylinder's two rims),
	// giving the specimen a rigid outer cage. Box and cylinder only.
	bool frameContainerEdges = false;
	float frameBeam = 0.3f;       // mm; edge beam radius

	// Cooperative cancellation hook. The GUI wires this to the running
	// GenerationTask's cancel flag before launching a job; the long-running
	// solves (calibration secants, generation pipeline) poll it at their
	// checkpoints and bail early, leaving the last valid built scaffold. Empty
	// (default) means "never cancel".
	std::function<bool()> cancelRequested;

	float stretchX{ 1.0f };
	float stretchY{ 1.0f };
	float stretchZ{ 1.0f };

	Vec3 anisotropyVec{ 1.0f, 0.0f, 0.0f };
	float anisotropyAngle{ 0.0f };
	std::unique_ptr<Ellipsoid> ellipsoidModel;
	std::vector<std::shared_ptr<AnisotropySource>> anisotropySources;
	// Constant weight of the background (global) covariance in the
	// partition-of-unity blend (see blend_metric in Anisotropy.h).
	// Far from all sources the metric is exactly the background; at a source
	// centre the background retains a residual wb/(1+wb) influence.
	// Raise to bleed the background into source regions more aggressively;
	// values are clamped away from 0. Typical range 0.05 - 0.5.
	float backgroundWeight{ 0.1f };

	int iter = 15;
	float lambda = 0.5f;
	float mu = -0.53f;

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

	std::vector<VoxelValue> cachedVoxels;

	// Per-voxel signed container distance, filled by compute_cached_field_values
	// and reused by assemble_field for the boundary clamp and the porosity domain
	// mask. A member (not a local) so the assemble pass can run standalone.
	std::vector<float> containerDistField;

	// scaffold properties
	float threshold = { 0.5f };
	float kSmooth = {0.5f};

	Aabb aabb;

	Logger* logger = nullptr;

	std::chrono::steady_clock::time_point start_time;
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
	uint32_t smiVersion = 0;

};

class ScaffoldFactory {

public:
	void launch();

	void gui_draw(
		GenerationTask* task,
		Logger* logger,
		const char* popupName, bool& showPopup,
		SelectedObject* selectedPanelObj, void*& selectedSceneObj,
		std::vector<std::unique_ptr<GeneratorLewiner>>& scaffoldList,
		std::vector<std::shared_ptr<IContainer>>& containers,
		std::vector<std::shared_ptr<InterfaceSeedGenerator>>& generators,
		std::vector<std::shared_ptr<AnisotropySource>>& anisoSources
	);

private:
	std::unique_ptr<GeneratorLewiner> pendingScaffold = nullptr;
	std::weak_ptr<IContainer> selectedCon;
	std::weak_ptr<InterfaceSeedGenerator> selectedGen;
	std::shared_ptr<IContainer> lockedCon;
	std::shared_ptr<InterfaceSeedGenerator> lockedGen;
	std::string name = "";
	std::string genContainerName = "";
	std::string genGeneratorName = "";
	std::chrono::steady_clock::time_point start_time;
	float thickness = { 0.3f };
	float openess = { 0.5f };
	float spread = { 0.5f };
	bool smoothJunctions = false;
	float smoothK = 0.01f;
	bool roundEdges = false;
	float edgeK = 0.01f;
	bool frameBoundary = false;
	float frameDepth = 1.0f;
	bool frameContainerEdges = false;
	float frameBeam = 0.3f;
	float stretchX = { 1.0f };
	float stretchY = { 1.0f };
	float stretchZ = { 1.0f };
	Vec3 anisotropyVec = { 1.0f, 0.0f, 0.0f };
	float anisotropyAngle = { 0.0f };
	std::vector<std::shared_ptr<AnisotropySource>> anisotropySources;
	float backgroundWeight = { 0.1f };

	std::array<int, 3> resolution = { 100, 100, 100 };
	float voxelSize{ 0.05f };
	float measurementVoxelSize{ 0.03f };
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

	std::vector<CalibrationStage> stages;

	bool calibrateDA = true;
	float targetDa = 1.2f;
	int targetFormulaIdx = 0;
	bool calibrateThickness = true;
	float targetThickness = 0.3f;
	bool calibratePorosity = true;
	float targetPorosity = 80.0f;
	// bool calibrateSmi = false;   // SMI calibration disabled (SMI is reported, not targeted)
	// float targetSmi = 0.5f;
	float calibrationStep = 0.02f;
	float calibrationTol = 0.001f;
	int calibrationIter = 10;
	float startThickness = 0.3f;
	float endThickness = 1.0f;
	float transitionDistance = 20.0f;
	
	int iter = 15;
	float lambda = 0.5f;
	float mu = -0.53f;

	float warningFlashTimer1 = 0.0f;
	float warningFlashTimer2 = 0.0f;

	void thickness_options();
	void smooth_options();
	void anisotropy_options(
		std::vector<std::shared_ptr<AnisotropySource>>& globalSources);
	void main_options(
		std::vector<std::shared_ptr<IContainer>>& containers,
		std::vector<std::shared_ptr<InterfaceSeedGenerator>>& generators
	);
};

#endif // ! "GENERATOR_LEWINER_H"
