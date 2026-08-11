#include <iostream>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <random>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib")
// wingdi.h defines ERROR as a macro, which collides with LogPriority::ERROR in
// Logger.h. Drop it (we don't use GDI here) before any project header.
#undef ERROR
#endif
#include "SeedGenerator/Container.h"
#include "SeedGenerator/SeedGenerator.h"
#include "SeedGenerator/RadiusCalculator.h"
#include "ScaffoldGenerator/GeneratorLewiner.h"
#include "ScaffoldGenerator/Anisotropy.h"
#include "Logger/Logger.h"

// Standalone CLI harness (no GUI/OpenGL window) to measure where time goes in
// the scaffold generation pipeline using a real mesh container.
//
// Usage: scaffoldProfile [meshPath|--box] [scale] [resX resY resZ] [flags]
//   meshPath              path to an STL container, or "--box" for a trivial
//                          BoxContainer (isolates seed-kdtree cost)
//   scale                  uniform scale applied to the mesh about its bbox
//                          center (default 100, e.g. cm STL → mm)
//   resX resY resZ         explicit grid resolution; overrides --voxel-size
//   --no-thickness         disable "distance from container" varied-thickness
//   --repro-bug            Poisson seeding + varied thickness for defect repro
//   --voxel-size v         target voxel size in mm (default 0.05)
//   --start-thickness v / --end-thickness v / --transition v / --stretch-y v
//   --poisson-radius v     Poisson seed spacing (default 0.583)
//   --aniso-source ox oy oz sx sy sz sigma
//                          add one anisotropy source (repeatable); origin
//                          (ox,oy,oz), stretch (sx,sy,sz), influence sigma
//   --background-weight w  background metric weight in blend (default 0.1)
//   --openness t           openness parameter tau in [0,1] (default 0.5)
//   --phantom T v          local-thickness self-test on an analytic slab of
//                          thickness T at voxel size v; prints measured vs true
//                          and exits (no container/seeds needed)
//   --calibrate t0 t1 n v  sweep iso-level (thickness) over n steps in [t0,t1],
//                          measuring Tb.Th at voxel size v; prints a CSV curve
//                          of iso-level -> Tb.Th and exits (use with --box)
//   --calibrate-target T v auto-tune iso-level so measured Tb.Th = T at voxel
//                          size v (secant loop); prints calibrated c_glob and
//                          final metrics, then exits
//   --box-size s           cubic box side in mm for --box (default: legacy box)
//   --poisson              uniform Poisson seeding at --poisson-radius
//   --stretch sx sy sz     anisotropy stretch factors
//   --spread s             per-face fenestration spread (rod/plate diversity)
//   --thickness t          iso-level / nominal strut thickness
//   --seed N               seed-placement RNG seed (0 = non-deterministic).
//                          Fixing it makes a run exactly reproducible.
//   --bonej-suite N v      emit a cohort of N scaffolds for the BoneJ cross-tool
//                          comparison, measured at voxel size v: per scaffold an
//                          STL, a binary NRRD (exported at the SAME v, so BoneJ
//                          sees the identical image), plus generation_parameters
//                          .csv and scaffoldbuilder_metrics.csv. Thickness and
//                          spacing are sampled from physiological ranges with
//                          spacing >= 2.5*thickness. Seeds are 1..N (reproducible).
//   --bonej-out dir        output directory (default data/bonej_comparison)
//   --bonej-thickness-range a b / --bonej-spacing-range a b   sampling ranges
//   --bonej-stretch-range a b   per-scaffold dominant-y anisotropy stretch range
//                          (default 1.0..2.5) so the cohort spans a DA range;
//                          a fixed --stretch overrides it. DA is stored in
//                          BoneJ's convention (1 - lambda_min/lambda_max).
//   --runtime-suite        12.3c performance: time generation (field + marching
//                          cubes) vs grid resolution (fixed seeds) and vs seed
//                          count (fixed grid). Use with --box --poisson. Writes
//                          runtime_vs_resolution.csv and runtime_vs_seeds.csv.
//   --runtime-repeats n    timed repeats per config, median reported (default 3).
//   --runtime-out dir      output dir (default doc/paper/experiments/runtime).
//   --sweep-input <name> a b n   sensitivity OAT: sweep one input over [a,b] in n
//                          steps (others at the CLI baseline), measure the full
//                          metric set. name = thickness|openness|stretch|spread|
//                          radius. Requires --poisson and a fixed --seed.
//   --sweep-grid n1 a b s n2 a b s   2D grid over two inputs (e.g. thickness x
//                          openness) for the coupling surface / Jacobian.
//   --sweep-radius-std rMin rMax a b n   stochastic-radius Poisson: sweep the
//                          per-seed radius std over [a,b] in n steps and report
//                          the resulting Tb.Sp mean and intra-sample std (the
//                          input->output transfer). Uses --replicates (default 3).
//   --sweep-out dir        output dir (default doc/paper/experiments/sensitivity).
//   --sweep-extra v        add an explicit voxel size to --sweep-voxel (repeatable).
//                          Use it for the uCT resolution so the comparison point
//                          is a real measurement, not an interpolation.
//   --sweep-voxel a b s    12.2b resolution convergence: generate ONCE, then
//                          re-measure every metric at image voxel sizes a..b in
//                          steps of s. Prints a CSV. SMI/BS/porosity are mesh or
//                          generation-grid quantities and stay flat (controls).
//                          Useful range ~0.01..0.05; past ~0.06 a 0.127mm wall is
//                          under 2 voxels and the metrics are meaningless.
//   --phantom-suite f v    analytic phantoms (slab/cylinder/sphere/torus) with
//                          known Tb.Th / SMI / Conn ground truth
//   --mil-directions N     MIL directions for Tb.N/DA (default 10000, = GUI)
//   --mil-lines N          MIL rays per direction (default 2000, = GUI)
//   --replicates N v       generate N independent realizations (seeds 1..N),
//                          measure all metrics at voxel size v, and report a
//                          per-run CSV plus mean +- SD. Use this rather than a
//                          single run: the seeding is stochastic, so one run is
//                          a single draw from a distribution of scaffolds.
int main(int argc, char** argv) {

	std::vector<std::string> args(argv + 1, argv + argc);

	bool noThickness = false;
	bool reproBug = false;
	float voxelSizeArg = 0.05f;
	float startThicknessArg = 5.0f;
	float endThicknessArg = 0.3f;
	float transitionArg = 5.0f;
	float stretchYArg = 2.0f;
	float poissonRadius = 0.583f;
	float backgroundWeightArg = 0.1f;
	float opennessArg = 0.5f;

	// Optional edge rounding (smin softening of the raw distance order
	// statistics). Off unless --round-edges is passed; edgeKArg=0 is a no-op.
	bool roundEdgesArg = false;
	float edgeKArg = 0.0f;

	bool phantomMode = false;
	float phantomThickness = 0.0f, phantomVoxel = 0.0f;

	bool phantomSuite = false;
	float phantomSuiteFeature = 0.3f, phantomSuiteVoxel = 0.0f;


	bool calibrateMode = false;
	float calStart = 0.0f, calEnd = 0.0f, calVoxel = 0.05f;
	int calSteps = 0;

	bool calibrateTargetMode = false;
	float calTarget = 0.0f, calTargetVoxel = 0.05f;

	// joint calibration test/scripting hooks
	bool twoKnobMode = false, threeKnobMode = false;
	float jkTh = 0.0f, jkPor = 0.0f, jkSMI = 0.0f, jkVoxel = 0.05f;
	float jkTol = 0.001f; int jkIter = 10;

	bool metricsMode = false;
	float metricsVoxel = 0.05f;

	// 12.2b resolution convergence: generate ONCE, re-measure at each voxel size
	bool sweepVoxel = false;
	float sweepMin = 0.01f, sweepMax = 0.05f, sweepStep = 0.01f;
	// extra voxel sizes injected into the sweep (e.g. the exact uCT resolution,
	// which must be a real measurement rather than an interpolation)
	std::vector<float> sweepExtra;

	// Input/output sensitivity study: OAT sweep of one input, or a 2D grid over
	// two inputs, measuring the full metric set each config. Names:
	// thickness | openness | stretch | spread | radius.
	bool sweepInput = false, sweepGrid = false;
	std::string sweepName, gridName1, gridName2;
	float sweepA = 0.0f, sweepB = 0.0f; int sweepN = 0;

	// --sweep-radius-std rMin rMax a b n : map the stochastic-radius input std to
	// the resulting Tb.Sp (mean and intra-sample std), since input std != output
	// Tb.Sp std (the packing regularises the variance).
	bool sweepRadiusStd = false;
	float rsRMin = 0.0f, rsRMax = 0.0f, rsA = 0.0f, rsB = 0.0f; int rsN = 0;
	float gridA1 = 0.0f, gridB1 = 0.0f, gridA2 = 0.0f, gridB2 = 0.0f;
	int gridN1 = 0, gridN2 = 0;
	std::string sweepOut = "doc/paper/experiments/sensitivity";

	// MIL sampling for Tb.N / DA (GUI defaults). Exposed because the sweep at
	// fine voxel sizes is dominated by ray marching.
	int milDirections = 10000;
	int milLines = 2000;

	// 12.2c BoneJ cross-tool comparison: emit a cohort of scaffolds as meshes +
	// binary images for independent analysis in BoneJ.
	bool bonejSuite = false;
	int bonejCount = 20;
	float bonejVoxel = 0.025f;
	std::string bonejOut = "data/bonej_comparison";

	// Literature replication suite (11 samples): builds each sample's container +
	// seeds and runs BOTH a manual (single-shot: DA + porosity iterated, iso fixed
	// at reported Tb.Th, no thickness calibration) and a fully-calibrated trial.
	// Edge rounding is on with edgeK = kEdgeFrac * reported Tb.Th per sample.
	bool literatureSuite = false;
	std::string literatureCsv = "doc/paper/experiments/literature/literature_samples.csv";
	std::string literatureOut = "doc/paper/experiments/literature";
	float literatureGenVoxel = 0.0f;    // generation grid spacing (mm); 0 = auto (Tb.Th/4)
	float kEdgeFrac = 1.0f;             // edgeK = kEdgeFrac * reported Tb.Th
	// physiologically plausible sampling ranges (Ulrich 1999 min..max)
	float bonejThMin = 0.10f, bonejThMax = 0.26f;
	float bonejSpMin = 0.45f, bonejSpMax = 1.00f;
	// per-scaffold anisotropy stretch (dominant y-axis) so the cohort spans a
	// real DA range instead of sitting at ~isotropic. 1.0 = isotropic.
	float bonejStretchMin = 1.0f, bonejStretchMax = 2.5f;
	// 12.3c runtime/performance benchmark: generation time vs grid resolution
	// and vs seed count.
	bool runtimeSuite = false;
	int runtimeRepeats = 3;
	std::string runtimeOut = "doc/paper/experiments/runtime";
	float thicknessArg = -1.0f; // override iso-level when > 0
	float spreadArg = -1.0f;    // override per-face fenestration spread when >= 0
	uint32_t seedArg = 0;       // seed-placement RNG seed; 0 = non-deterministic
	int replicates = 0;         // >0 -> run N realizations and report mean +- SD

	bool poissonFlag = false;
	float boxSize = 0.0f; // cubic box side for --box (0 = use hardcoded default)
	bool useCylinder = false; float cylR = 2.5f, cylH = 10.0f; // --cylinder r h
	float stretchX = 1.0f, stretchY = 1.0f, stretchZ = 1.0f;
	bool stretchSet = false;

	// each entry: {ox, oy, oz, sx, sy, sz, sigma}
	std::vector<std::array<float, 7>> anisoSourceArgs;

	std::vector<std::string> positional;

	for (size_t i = 0; i < args.size(); i++) {
		const std::string& a = args[i];
		if (a == "--no-thickness") noThickness = true;
		else if (a == "--repro-bug") reproBug = true;
		else if (a == "--voxel-size" && i + 1 < args.size()) voxelSizeArg = std::stof(args[++i]);
		else if (a == "--start-thickness" && i + 1 < args.size()) startThicknessArg = std::stof(args[++i]);
		else if (a == "--end-thickness" && i + 1 < args.size()) endThicknessArg = std::stof(args[++i]);
		else if (a == "--transition" && i + 1 < args.size()) transitionArg = std::stof(args[++i]);
		else if (a == "--stretch-y" && i + 1 < args.size()) stretchYArg = std::stof(args[++i]);
		else if (a == "--poisson-radius" && i + 1 < args.size()) poissonRadius = std::stof(args[++i]);
		else if (a == "--background-weight" && i + 1 < args.size()) backgroundWeightArg = std::stof(args[++i]);
		else if (a == "--openness" && i + 1 < args.size()) opennessArg = std::stof(args[++i]);
		else if (a == "--round-edges" && i + 1 < args.size()) { edgeKArg = std::stof(args[++i]); roundEdgesArg = edgeKArg > 0.0f; }
		else if (a == "--phantom" && i + 2 < args.size()) {
			phantomMode = true;
			phantomThickness = std::stof(args[++i]);
			phantomVoxel = std::stof(args[++i]);
		}
		else if (a == "--phantom-suite" && i + 2 < args.size()) {
			phantomSuite = true;
			phantomSuiteFeature = std::stof(args[++i]);
			phantomSuiteVoxel = std::stof(args[++i]);
		}
		else if (a == "--calibrate" && i + 4 < args.size()) {
			calibrateMode = true;
			calStart = std::stof(args[++i]);
			calEnd = std::stof(args[++i]);
			calSteps = std::stoi(args[++i]);
			calVoxel = std::stof(args[++i]);
		}
		else if (a == "--calibrate-target" && i + 2 < args.size()) {
			calibrateTargetMode = true;
			calTarget = std::stof(args[++i]);
			calTargetVoxel = std::stof(args[++i]);
		}
		else if (a == "--two-knob" && i + 3 < args.size()) {
			twoKnobMode = true;
			jkTh = std::stof(args[++i]);
			jkPor = std::stof(args[++i]);
			jkVoxel = std::stof(args[++i]);
			// optional: tol iter (to replicate GUI calibrationTol/Iter)
			if (i + 2 < args.size() && args[i + 1][0] != '-') {
				jkTol = std::stof(args[++i]);
				jkIter = std::stoi(args[++i]);
			}
		}
		else if (a == "--three-knob" && i + 4 < args.size()) {
			threeKnobMode = true;
			jkTh = std::stof(args[++i]);
			jkPor = std::stof(args[++i]);
			jkSMI = std::stof(args[++i]);
			jkVoxel = std::stof(args[++i]);
		}
		else if (a == "--metrics" && i + 1 < args.size()) {
			metricsMode = true;
			metricsVoxel = std::stof(args[++i]);
		}
		else if (a == "--sweep-voxel" && i + 3 < args.size()) {
			sweepVoxel = true;
			sweepMin = std::stof(args[++i]);
			sweepMax = std::stof(args[++i]);
			sweepStep = std::stof(args[++i]);
		}
		else if (a == "--bonej-suite" && i + 2 < args.size()) {
			bonejSuite = true;
			bonejCount = std::stoi(args[++i]);
			bonejVoxel = std::stof(args[++i]);
		}
		else if (a == "--bonej-out" && i + 1 < args.size()) bonejOut = args[++i];
		else if (a == "--literature-suite") literatureSuite = true;
		else if (a == "--literature-csv" && i + 1 < args.size()) literatureCsv = args[++i];
		else if (a == "--literature-out" && i + 1 < args.size()) literatureOut = args[++i];
		else if (a == "--literature-genvoxel" && i + 1 < args.size()) literatureGenVoxel = std::stof(args[++i]);
		else if (a == "--kedge-frac" && i + 1 < args.size()) kEdgeFrac = std::stof(args[++i]);
		else if (a == "--bonej-thickness-range" && i + 2 < args.size()) {
			bonejThMin = std::stof(args[++i]);
			bonejThMax = std::stof(args[++i]);
		}
		else if (a == "--bonej-spacing-range" && i + 2 < args.size()) {
			bonejSpMin = std::stof(args[++i]);
			bonejSpMax = std::stof(args[++i]);
		}
		else if (a == "--bonej-stretch-range" && i + 2 < args.size()) {
			bonejStretchMin = std::stof(args[++i]);
			bonejStretchMax = std::stof(args[++i]);
		}
		else if (a == "--runtime-suite") runtimeSuite = true;
		else if (a == "--runtime-repeats" && i + 1 < args.size()) runtimeRepeats = std::stoi(args[++i]);
		else if (a == "--runtime-out" && i + 1 < args.size()) runtimeOut = args[++i];
		else if (a == "--sweep-input" && i + 4 < args.size()) {
			sweepInput = true;
			sweepName = args[++i];
			sweepA = std::stof(args[++i]);
			sweepB = std::stof(args[++i]);
			sweepN = std::stoi(args[++i]);
		}
		else if (a == "--sweep-grid" && i + 8 < args.size()) {
			sweepGrid = true;
			gridName1 = args[++i]; gridA1 = std::stof(args[++i]); gridB1 = std::stof(args[++i]); gridN1 = std::stoi(args[++i]);
			gridName2 = args[++i]; gridA2 = std::stof(args[++i]); gridB2 = std::stof(args[++i]); gridN2 = std::stoi(args[++i]);
		}
		else if (a == "--sweep-radius-std" && i + 5 < args.size()) {
			sweepRadiusStd = true;
			rsRMin = std::stof(args[++i]);
			rsRMax = std::stof(args[++i]);
			rsA = std::stof(args[++i]);
			rsB = std::stof(args[++i]);
			rsN = std::stoi(args[++i]);
		}
		else if (a == "--sweep-out" && i + 1 < args.size()) sweepOut = args[++i];
		else if (a == "--sweep-extra" && i + 1 < args.size()) sweepExtra.push_back(std::stof(args[++i]));
		else if (a == "--mil-directions" && i + 1 < args.size()) milDirections = std::stoi(args[++i]);
		else if (a == "--mil-lines" && i + 1 < args.size()) milLines = std::stoi(args[++i]);
		else if (a == "--thickness" && i + 1 < args.size()) thicknessArg = std::stof(args[++i]);
		else if (a == "--spread" && i + 1 < args.size()) spreadArg = std::stof(args[++i]);
		else if (a == "--seed" && i + 1 < args.size()) seedArg = static_cast<uint32_t>(std::stoul(args[++i]));
		else if (a == "--replicates" && i + 2 < args.size()) {
			replicates = std::stoi(args[++i]);
			metricsVoxel = std::stof(args[++i]);
		}
		else if (a == "--poisson") poissonFlag = true;
		else if (a == "--box-size" && i + 1 < args.size()) boxSize = std::stof(args[++i]);
		else if (a == "--cylinder" && i + 2 < args.size()) {
			useCylinder = true;
			cylR = std::stof(args[++i]);
			cylH = std::stof(args[++i]);
		}
		else if (a == "--stretch" && i + 3 < args.size()) {
			stretchX = std::stof(args[++i]);
			stretchY = std::stof(args[++i]);
			stretchZ = std::stof(args[++i]);
			stretchSet = true;
		}
		else if (a == "--aniso-source" && i + 7 < args.size()) {
			std::array<float, 7> s;
			for (int k = 0; k < 7; k++) s[k] = std::stof(args[++i]);
			anisoSourceArgs.push_back(s);
		}
		else positional.push_back(a);
	}

	std::string meshPath = positional.size() > 0 ? positional[0] : "data/femoral_head_container.stl";
	bool useBox = (meshPath == "--box");
	float scale = positional.size() > 1 ? std::stof(positional[1]) : 100.0f;

	Logger& logger = Logger::get_instance();

	// --- Literature replication suite (11 samples) -------------------------
	// Per CSV row: build the sample's container + uniform-Poisson seeds (radius =
	// reported Tb.Sp), turn edge rounding on (edgeK = kEdgeFrac * Tb.Th), and run
	// TWO trials:
	//   manual     : DA + porosity solved (a user iterating stretch/openness), but
	//                iso FIXED at the reported Tb.Th (no thickness calibration) ->
	//                exposes the +junction-bleeding Tb.Th bias.
	//   calibrated : full closed loop (DA + thickness + porosity).
	// DA uses each sample's own formula (Lmax/Lmin idx 0, or lambda_min/max idx 2),
	// and Tb.N is reported with the source-matched formula for a fair comparison.
	if (literatureSuite) {
		std::filesystem::create_directories(literatureOut);
		std::ifstream fin(literatureCsv);
		if (!fin.is_open()) { std::cerr << "cannot open " << literatureCsv << "\n"; return 1; }

		std::vector<std::string> header;
		std::vector<std::vector<std::string>> rows;
		std::string line;
		while (std::getline(fin, line)) {
			if (line.empty() || line[0] == '#') continue;
			std::vector<std::string> cells; std::stringstream ss(line); std::string c;
			while (std::getline(ss, c, ',')) cells.push_back(c);
			if (header.empty()) header = cells; else rows.push_back(cells);
		}
		auto col = [&](const std::vector<std::string>& r, const std::string& name) -> std::string {
			for (size_t i = 0; i < header.size(); ++i) if (header[i] == name && i < r.size()) return r[i];
			return "";
		};
		auto colf = [&](const std::vector<std::string>& r, const std::string& name, float def) -> float {
			std::string v = col(r, name); return v.empty() ? def : std::stof(v);
		};

		std::ofstream fout(literatureOut + "/literature_results.csv");
		fout << "sample,trial,container,gen_voxel,meas_voxel,edgeK,converged,"
			<< "iso,stretch,openness,TbTh,TbTh_SD,TbTh_ref,TbSp,TbSp_SD,TbSp_ref,porosity,porosity_ref,"
			<< "DA,DA_ref,TbN_mil,TbN_src,TbN_ref,ConnD,SMI,BS_BV,BS_TV\n";
		std::cout << "\n=== LITERATURE SUITE (" << rows.size() << " samples, kEdge="
			<< kEdgeFrac << "*Tb.Th) ===" << std::endl;

		for (auto& r : rows) {
			const std::string sid = col(r, "sample_id");
			const std::string ctype = col(r, "container_type");
			const float tbth = colf(r, "tb_th", 0.1f);
			const float tbsp = colf(r, "tb_sp", 0.7f);
			const float porRef = colf(r, "porosity_pct", 85.0f);   // percent
			const float daRef = colf(r, "da", 1.5f);
			const float voxel = colf(r, "microct_voxel_mm", 0.02f);
			const std::string daForm = col(r, "da_formula");
			const std::string tbnForm = col(r, "tbn_formula");
			const int formIdx = (daForm == "lmin_lmax") ? 2 : 0;
			const float edgeK = kEdgeFrac * tbth;

			IContainer* con = nullptr; std::string cdesc;
			if (ctype == "cylinder") {
				float d = colf(r, "cyl_diameter_mm", 5.0f), h = colf(r, "cyl_height_mm", 10.0f);
				con = new CylinderContainer(d * 0.5f, h, false);
				cdesc = "cyl_d" + col(r, "cyl_diameter_mm") + "_h" + col(r, "cyl_height_mm");
			} else {
				float s = colf(r, "box_side_mm", 4.0f);
				con = new BoxContainer(Vec3(s, s, s), Vec3(0.5f * s, 0.5f * s, 0.5f * s), false);
				cdesc = "box" + col(r, "box_side_mm");
			}
			std::shared_ptr<IContainer> conSh(con, [](IContainer*) {});
			Bounds b = con->compute_bounds();
			std::array<float, 6> bnds = { (float)b.xMin,(float)b.xMax,(float)b.yMin,(float)b.yMax,(float)b.zMin,(float)b.zMax };
			// Generation grid = reported Tb.Th / 4: four voxels across the thinnest
			// feature (Nyquist-comfortable), and since edgeK = Tb.Th = 4h the edge
			// rounding sits well above its aliasing floor (2h). Self-scales with the
			// sample's thickness, so big containers still stay tractable. A positive
			// --literature-genvoxel overrides it.
			float gVox = (literatureGenVoxel > 0.0f) ? literatureGenVoxel : (tbth / 4.0f);
			std::array<int, 3> res = {
				std::max(2,(int)((b.xMax - b.xMin) / gVox)),
				std::max(2,(int)((b.yMax - b.yMin) / gVox)),
				std::max(2,(int)((b.zMax - b.zMin) / gVox)) };

			// Build the seed generator as a shared object and ATTACH it to each
			// scaffold (g.generator) so export_scaf stores a valid generator block
			// (genTypeID = uniform Poisson + radius + name). Without this the .scaf
			// would carry genTypeID = -1 and load as "generator not found",
			// blocking edits in the GUI.
			auto poisson = std::make_shared<Poisson3D>(tbsp, tbsp, 30, false);
			poisson->set_rng_seed(1);
			poisson->run(*con);
			std::vector<Vec3> seeds = poisson->get_seeds();
			if (seeds.size() < 4) { std::cerr << sid << ": too few seeds\n"; delete con; continue; }
			poisson->set_seeds(seeds);
			poisson->name = sid + "_poisson";
			poisson->type = ObjectType::UniformGeneratorType;

			std::cout << "\n--- " << sid << " (" << cdesc << ", genVox=" << gVox
				<< ", seeds=" << seeds.size() << ", edgeK=" << edgeK << ") ---" << std::endl;

			auto run_trial = [&](const std::string& trial, bool calThick, const std::string& exportPrefix) {
				GeneratorLewiner g(seeds, bnds, res, &logger, 0.5f, tbth, false);
				g.container = conSh;
				g.generator = poisson;        // attach seed generator so .scaf is editable on load
				g.iter = 100;                 // heavy Taubin smoothing (release-quality surfaces)
				g.roundEdges = true; g.edgeK = edgeK;
				g.set_stretch(1.0f, 1.0f, 1.0f);
				g.targetFormulaIdx = formIdx;
				g.calibrateDA = true;        g.targetDa = daRef;
				g.calibratePorosity = true;  g.targetPorosity = porRef;  // percent
				g.calibrateThickness = calThick; g.targetThickness = tbth;
				g.build_calibration_stages(voxel);
				bool ok = g.solve_calibration(g.stages, 0, *con);
				g.marching_cubes(true);
				g.estimate_metrics(*con);                 // mesh metrics (BS/BV, BS/TV, SMI)
				Aabb bb = g.get_aabb();
				std::array<float, 6> bA = { bb.pMin.x, bb.pMax.x, bb.pMin.y, bb.pMax.y, bb.pMin.z, bb.pMax.z };
				std::array<float, 6> bB = bA;
				g.estimate_local_thickness(voxel, bA, false);
				g.estimate_local_thickness(voxel, bB, true);
				g.estimate_anisotropy(voxel, milDirections, milLines, formIdx);
				g.estimate_trabecular_number(voxel, 0, milDirections, milLines);
				g.estimate_connectivity_density_voxel(voxel, 6);
				g.estimate_smi();
				const float bvtv = 1.0f - g.porosity;
				const float tbnSrc = (tbnForm == "inv_tbsp") ? 1.0f / g.localSeparation
					: (tbnForm == "bvtv_tbth") ? bvtv / g.localThickness
					: g.trabecularNr;
				fout << sid << "," << trial << "," << cdesc << "," << gVox << "," << voxel << ","
					<< edgeK << "," << (ok ? 1 : 0) << "," << g.get_iso_level() << "," << g.stretchX << ","
					<< g.get_openness() << "," << g.localThickness << "," << g.localThicknessStd << "," << tbth << ","
					<< g.localSeparation << "," << g.localSeparationStd << "," << tbsp << "," << (g.porosity * 100.0f) << "," << porRef << ","
					<< g.anisotropyDegree << "," << daRef << "," << g.trabecularNr << "," << tbnSrc << ","
					<< col(r, "tb_n") << "," << g.connectivityDensity << "," << g.smi << ","
					<< g.surfaceToVolume << "," << g.surfaceToTotalVolume << "\n";
				fout.flush();
				std::cout << "  [" << trial << "] conv=" << ok << " iso=" << g.get_iso_level()
					<< " stretch=" << g.stretchX << " -> TbTh=" << g.localThickness
					<< " (ref " << tbth << ") por=" << (g.porosity * 100.0f) << " (ref " << porRef
					<< ") DA=" << g.anisotropyDegree << " (ref " << daRef << ")" << std::endl;

				// Release artifacts (from the calibrated trial): scaffold + metrics +
				// a self-contained generation-recipe CSV (export_parameters needs a
				// seed-generator handle we do not have here, so it is written inline).
				if (!exportPrefix.empty()) {
					g.export_scaf(exportPrefix + ".scaf");
					g.export_metrics(exportPrefix + "_metrics.csv");
					std::ofstream fp(exportPrefix + "_parameters.csv");
					if (fp.is_open()) {
						fp << "key,value\n"
							<< "sample," << sid << "\n"
							<< "container," << cdesc << "\n"
							<< "generation_voxel_mm," << gVox << "\n"
							<< "measurement_voxel_mm," << voxel << "\n"
							<< "poisson_radius_mm," << tbsp << "\n"
							<< "seed_count," << seeds.size() << "\n"
							<< "iso_level_mm," << g.get_iso_level() << "\n"
							<< "stretch_x," << g.stretchX << "\n"
							<< "openness," << g.get_openness() << "\n"
							<< "spread," << g.get_spread() << "\n"
							<< "round_edges,1\n"
							<< "edgeK_mm," << edgeK << "\n"
							<< "taubin_iter," << g.iter << "\n"
							<< "taubin_lambda," << g.lambda << "\n"
							<< "taubin_mu," << g.mu << "\n"
							<< "target_TbTh_mm," << tbth << "\n"
							<< "target_porosity_pct," << porRef << "\n"
							<< "target_DA," << daRef << "\n"
							<< "da_formula_idx," << formIdx << "\n";
						fp.close();
					}
				}
			};

			run_trial("manual", false, "");
			run_trial("calibrated", true, literatureOut + "/" + sid);
			delete con;
		}
		std::cout << "\nWrote " << literatureOut << "/literature_results.csv" << std::endl;
		return 0;
	}

	// --- Phantom self-test: measure an analytic slab of known thickness ---
	// A slab has no junctions, so the estimator should return the true
	// thickness (up to the ~1-voxel distance-transform bias). This isolates
	// measurement accuracy from the junction-bleeding that elevates Tb.Th on
	// real foam/lattice structures.
	if (phantomMode) {
		GeneratorLewiner gen(false);
		gen.set_logger(&logger);
		std::cout << "\n=== SLAB PHANTOM SELF-TEST ===" << std::endl;
		std::cout << "true thickness = " << phantomThickness << " mm, voxel size = "
			<< phantomVoxel << " mm (" << (phantomThickness / phantomVoxel)
			<< " voxels across the slab)" << std::endl;
		std::array<float, 2> r = gen.run_slab_phantom(phantomThickness, phantomVoxel);
		float bias = r[0] - phantomThickness;
		std::cout << "measured Tb.Th   = " << r[0] << " +- " << r[1] << " mm" << std::endl;
		std::cout << "measurement bias = " << bias << " mm ("
			<< (100.0f * bias / phantomThickness) << " %, "
			<< (bias / phantomVoxel) << " voxels)" << std::endl;
		return 0;
	}

	// --- Phantom SUITE: verify every estimator against analytic ground truth --
	// Each shape has a known local thickness (= feature) and a known SMI/Conn.
	// This isolates measurement accuracy from the generator entirely.
	if (phantomSuite) {
		GeneratorLewiner gen(false);
		gen.set_logger(&logger);
		const float f = phantomSuiteFeature;
		const float v = phantomSuiteVoxel;
		const float dom = std::max(1.0f, 6.0f * f);   // matches build_phantom_field default
		const float domVol = dom * dom * dom;
		std::cout << "\n=== PHANTOM SUITE (feature " << f << " mm, voxel " << v
			<< " mm = " << (f / v) << " voxels) ===" << std::endl;
		std::cout << "shape,Tb.Th,Tb.Th_true,SMI,SMI_true,Conn,Conn_true" << std::endl;

		struct Case { int shape; const char* name; float smiTrue; float connTrue; };
		Case cases[] = {
			{ GeneratorLewiner::PHANTOM_SLAB,     "slab",     0.0f, 0.0f },
			{ GeneratorLewiner::PHANTOM_CYLINDER, "cylinder", 3.0f, 0.0f },
			{ GeneratorLewiner::PHANTOM_SPHERE,   "sphere",   4.0f, 0.0f },
			{ GeneratorLewiner::PHANTOM_TORUS,    "torus",    3.0f, 1.0f },  // tube~rod, 1 loop
		};

		for (const Case& c : cases) {
			gen.build_phantom_field(c.shape, f, v);
			Aabb bb = gen.get_aabb();
			std::array<float, 6> bnds = { bb.pMin.x, bb.pMax.x, bb.pMin.y, bb.pMax.y, bb.pMin.z, bb.pMax.z };
			gen.estimate_local_thickness(v, bnds, false);   // voxel metric
			gen.marching_cubes();                           // mesh for SMI
			gen.estimate_smi();
			gen.estimate_connectivity_density_voxel(v, 6);
			// raw loop count = Conn.D * domainVolume (ground truth is an integer)
			float rawConn = gen.connectivityDensity * domVol;
			std::cout << c.name << "," << gen.localThickness << "," << f << ","
				<< gen.smi << "," << c.smiTrue << ","
				<< rawConn << "," << c.connTrue << std::endl;
		}
		return 0;
	}

	IContainer* container = nullptr;
	if (useCylinder) {
		std::cout << "Using CylinderContainer, radius " << cylR << " mm, height " << cylH << " mm" << std::endl;
		container = new CylinderContainer(cylR, cylH, false);
	}
	else if (useBox) {
		if (boxSize > 0.0f) {
			std::cout << "Using cubic BoxContainer, side " << boxSize << " mm" << std::endl;
			container = new BoxContainer(Vec3(boxSize, boxSize, boxSize),
				Vec3(0.5f * boxSize, 0.5f * boxSize, 0.5f * boxSize), false);
		}
		else {
			std::cout << "Using trivial BoxContainer (isolates seed-kdtree cost)" << std::endl;
			container = new BoxContainer(Vec3(45.14f, 49.11f, 50.92f), Vec3(0.0f, -2.67f, 3.40f), false);
		}
	}
	else {
		std::cout << "Loading mesh: " << meshPath << std::endl;
		AbstractContainer* mesh = new AbstractContainer(meshPath, false);
		IContainer* meshAsContainer = mesh;
		Bounds rawBounds = meshAsContainer->compute_bounds();
		std::cout << "Raw mesh bounds: [" << rawBounds.xMin << ", " << rawBounds.xMax << "]" << std::endl;
		std::cout << "Applying scale x" << scale << std::endl;
		mesh->set_scale(scale);
		container = mesh;
	}

	// Non-owning shared_ptr so we can populate each generator's `container`
	// weak_ptr. Without it, container.lock() is null in the metric estimators and
	// the outside-container mask is skipped - which inflates Tb.Sp/Tb.N for
	// non-box containers (their AABB corners get counted as marrow). The GUI sets
	// this via scaffold->container; the profiler must too. Empty deleter: the raw
	// `delete container` at each mode's exit still owns the lifetime.
	std::shared_ptr<IContainer> containerShared(container, [](IContainer*) {});

	Bounds b = container->compute_bounds();
	std::cout << "Scaled bounds: ["
		<< b.xMin << ", " << b.xMax << "] x ["
		<< b.yMin << ", " << b.yMax << "] x ["
		<< b.zMin << ", " << b.zMax << "]" << std::endl;

	std::array<float, 6> bounds = {
		(float)b.xMin, (float)b.xMax, (float)b.yMin, (float)b.yMax, (float)b.zMin, (float)b.zMax
	};

	std::array<int, 3> res;
	if (positional.size() > 4) {
		res = { std::stoi(positional[2]), std::stoi(positional[3]), std::stoi(positional[4]) };
	}
	else {
		res = {
			std::max(2, (int)((b.xMax - b.xMin) / voxelSizeArg)),
			std::max(2, (int)((b.yMax - b.yMin) / voxelSizeArg)),
			std::max(2, (int)((b.zMax - b.zMin) / voxelSizeArg))
		};
	}
	size_t scaledVoxels = (size_t)res[0] * res[1] * res[2];
	std::cout << "Profiling resolution (voxel size " << voxelSizeArg << "): "
		<< res[0] << "x" << res[1] << "x" << res[2]
		<< " = " << scaledVoxels << " voxels" << std::endl;

	// Build one realization of the seed cloud. 's' is the placement RNG seed:
	// 0 -> non-deterministic (a different arrangement every call), non-zero ->
	// reproducible. The Poisson radius fixes the minimum spacing; 's' picks
	// which of the many valid arrangements at that spacing you get.
	auto make_seeds = [&](uint32_t s, bool verbose) -> std::vector<Vec3> {
		std::vector<Vec3> out;
		if (reproBug || poissonFlag) {
			Poisson3D poisson(poissonRadius, poissonRadius, 30, false);
			poisson.set_rng_seed(s);
			poisson.run(*container);
			out = poisson.get_seeds();
			if (verbose) std::cout << "Poisson-uniform seed count (radius=" << poissonRadius
				<< ", seed=" << s << "): " << out.size() << std::endl;
		}
		else {
			std::mt19937 rng(s != 0 ? s : 42u);
			std::uniform_real_distribution<float> dx(b.xMin, b.xMax), dy(b.yMin, b.yMax), dz(b.zMin, b.zMax);
			int seedNr = 2000;
			out.reserve(seedNr);
			for (int i = 0; i < seedNr; i++) {
				out.push_back(Vec3(dx(rng), dy(rng), dz(rng)));
			}
			if (verbose) std::cout << "Seed count: " << seedNr << " (seed=" << s << ")" << std::endl;
		}
		return out;
	};

	float isoLevel = (thicknessArg > 0.0f) ? thicknessArg : 0.255f;

	// --- Radius-std sweep: stochastic-radius input std -> resulting Tb.Sp ----
	// The stochastic varied-Poisson mode draws each seed's spacing from
	// N(midpoint(rMin,rMax), radiusStd). Because radius->Tb.Sp is a monotone
	// transfer and the packing regularises variance, the INPUT std is not the
	// OUTPUT Tb.Sp std - this sweep measures that transfer so a target
	// intra-sample Tb.Sp spread can be dialled in. TbSp_intraStd is the metric's
	// own spatial std within one sample (the heterogeneity this feature adds).
	if (sweepRadiusStd) {
		std::cout << "\n=== RADIUS-STD SWEEP (rMin=" << rsRMin << ", rMax=" << rsRMax
			<< ", measure voxel=" << metricsVoxel << ") ===" << std::endl;
		std::cout << "radiusStd,seedCount,TbSp_mean,TbSp_intraStd" << std::endl;

		const int reps = (replicates > 0) ? replicates : 3;
		for (int is = 0; is < rsN; ++is) {
			float rstd = (rsN <= 1) ? rsA : (rsA + (rsB - rsA) * is / (rsN - 1));

			double spSum = 0.0, stdSum = 0.0; long long cntSum = 0; int ok = 0;
			for (int r = 1; r <= reps; ++r) {
				Poisson3D poisson(rsRMin, rsRMax, 30, false);
				poisson.stochasticRadius = true;
				poisson.radiusMean = 0.0;              // 0 -> range midpoint
				poisson.radiusStd = rstd;
				poisson.set_rng_seed(static_cast<uint32_t>(r));
				RunConfig cfg;                          // null dist/rad: stochastic draw in varied_run()
				// Explicit-config entry: honours the stochasticRadius member set
				// above. The single-arg run() would dispatch on type (NoneType ->
				// uniform_run) and ignore stochastic spacing entirely.
				poisson.varied_run(*container, cfg);
				std::vector<Vec3> s = poisson.get_seeds();
				if (s.size() < 3) continue;

				GeneratorLewiner g(s, bounds, res, &logger, opennessArg, isoLevel, false);
				g.container = containerShared;          // container mask for Tb.Sp
				if (stretchSet) g.set_stretch(stretchX, stretchY, stretchZ);
				if (spreadArg >= 0.0f) g.spread = spreadArg;
				if (!g.compute_scalar_field(*container) || !g.marching_cubes()) continue;

				Aabb bb = g.get_aabb();
				std::array<float, 6> bB = { bb.pMin.x, bb.pMax.x, bb.pMin.y, bb.pMax.y, bb.pMin.z, bb.pMax.z };
				g.estimate_local_thickness(metricsVoxel, bB, true);   // separation = Tb.Sp

				spSum += g.localSeparation;
				stdSum += g.localSeparationStd;
				cntSum += static_cast<long long>(s.size());
				ok++;
			}
			if (ok > 0) {
				std::cout << rstd << "," << (cntSum / ok) << ","
					<< (spSum / ok) << "," << (stdSum / ok) << std::endl;
			}
		}

		delete container;
		return 0;
	}

	// --- Replicates: N independent realizations of the same recipe ---------
	// The seed placement is stochastic, so one run is a single draw from a
	// distribution of statistically-equivalent scaffolds. Report mean +- SD
	// across realizations, the analogue of the literature's mean +- SD across
	// specimens.
	if (replicates > 0) {
		std::cout << "\n=== REPLICATES (n=" << replicates
			<< ", measure voxel=" << metricsVoxel << ") ===" << std::endl;
		std::cout << "seed,TbTh,TbSp,TbN,DA,SMI,BS/BV,BS/TV,ConnD,porosity" << std::endl;

		std::vector<float> tbth, tbsp, tbn, da, smiV, bsbv, bstv, connd, por;

		for (int r = 1; r <= replicates; ++r) {
			std::vector<Vec3> s = make_seeds(static_cast<uint32_t>(r), false);
			GeneratorLewiner g(s, bounds, res, &logger, opennessArg, isoLevel, false);
			g.container = containerShared;   // enable the outside-container mask
			if (stretchSet) g.set_stretch(stretchX, stretchY, stretchZ);
			if (spreadArg >= 0.0f) g.spread = spreadArg;

			if (!g.compute_scalar_field(*container) || !g.marching_cubes()) {
				std::cerr << "replicate " << r << " failed, skipping." << std::endl;
				continue;
			}
			g.estimate_metrics(*container);
			Aabb bb = g.get_aabb();
			std::array<float, 6> bA = { bb.pMin.x, bb.pMax.x, bb.pMin.y, bb.pMax.y, bb.pMin.z, bb.pMax.z };
			std::array<float, 6> bB = bA;
			g.estimate_local_thickness(metricsVoxel, bA, false);
			g.estimate_local_thickness(metricsVoxel, bB, true);
			g.estimate_anisotropy(metricsVoxel, milDirections, milLines, 0);
			g.estimate_trabecular_number(metricsVoxel, 0, milDirections, milLines);
			g.estimate_connectivity_density_voxel(metricsVoxel, 6);
			g.estimate_smi();

			tbth.push_back(g.localThickness);   tbsp.push_back(g.localSeparation);
			tbn.push_back(g.trabecularNr);      da.push_back(g.anisotropyDegree);
			smiV.push_back(g.smi);              bsbv.push_back(g.surfaceToVolume);
			bstv.push_back(g.surfaceToTotalVolume);
			connd.push_back(g.connectivityDensity); por.push_back(g.porosity);

			std::cout << r << "," << g.localThickness << "," << g.localSeparation << ","
				<< g.trabecularNr << "," << g.anisotropyDegree << "," << g.smi << ","
				<< g.surfaceToVolume << "," << g.surfaceToTotalVolume << ","
				<< g.connectivityDensity << "," << g.porosity << std::endl;
		}

		// sample SD (n-1), matching how the literature reports spread across specimens
		auto mean_sd = [](const std::vector<float>& v) {
			double m = 0.0;
			for (float x : v) m += x;
			m /= (v.empty() ? 1 : v.size());
			double s = 0.0;
			for (float x : v) s += (x - m) * (x - m);
			s = (v.size() > 1) ? std::sqrt(s / (v.size() - 1)) : 0.0;
			return std::make_pair(m, s);
		};

		auto report = [&](const char* label, const std::vector<float>& v) {
			auto ms = mean_sd(v);
			std::cout << "  " << label << " = " << ms.first << " +- " << ms.second << std::endl;
		};

		std::cout << "\n--- mean +- SD over " << tbth.size() << " realizations ---" << std::endl;
		report("Tb.Th   ", tbth);
		report("Tb.Sp   ", tbsp);
		report("Tb.N    ", tbn);
		report("DA      ", da);
		report("SMI     ", smiV);
		report("BS/BV   ", bsbv);
		report("BS/TV   ", bstv);
		report("Conn.D  ", connd);
		report("porosity", por);

		delete container;
		return 0;
	}

	// --- Input/output sensitivity study ------------------------------------
	// Characterize how each generation input moves each output, to (1) validate
	// the workflow decoupling (structure/anisotropy set first, then a 2D solve
	// over thickness x openness), (2) extract the local Jacobian for the joint
	// calibrator, and (3) produce the parameter-sensitivity figure. Seeds are
	// held fixed (--seed) so a change in output comes from the input, not seed
	// noise. Measurement voxel = generation voxel (--metrics-voxel) for speed;
	// relative sensitivities are resolution-robust. Requires --poisson.
	if (sweepInput || sweepGrid) {
		std::filesystem::create_directories(sweepOut);
		if (!poissonFlag) {
			std::cerr << "sensitivity study requires --poisson (uniform Poisson seeding at Tb.Sp)\n";
			delete container; return 1;
		}

		// Baseline operating point from the CLI args (default: Ulrich sample 1).
		const float th0 = isoLevel;
		const float op0 = opennessArg;
		const float sy0 = stretchSet ? stretchY : 1.0f;
		const float sp0 = (spreadArg >= 0.0f) ? spreadArg : 0.5f;
		const float rad0 = poissonRadius;

		// Generate one config and append a CSV row of all inputs + all outputs.
		auto measure_and_write = [&](std::ofstream& fout, const std::string& swept,
			float th, float openness, float sy, float spread, float radius) {
			poissonRadius = radius;                  // make_seeds reads this (by ref)
			std::vector<Vec3> s = make_seeds(seedArg, false);
			if (s.size() < 3) { std::cerr << "  too few seeds, skipping\n"; return; }
			GeneratorLewiner g(s, bounds, res, &logger, openness, th, false);
			g.container = containerShared;           // enable the outside-container mask
			g.set_stretch(1.0f, sy, 1.0f);           // dominant-y anisotropy
			g.spread = spread;
			if (roundEdgesArg) { g.roundEdges = true; g.edgeK = edgeKArg; }
			if (!g.compute_scalar_field(*container) || !g.marching_cubes()) {
				std::cerr << "  generation failed, skipping\n"; return;
			}
			g.estimate_metrics(*container);
			Aabb bb = g.get_aabb();
			std::array<float, 6> bA = { bb.pMin.x, bb.pMax.x, bb.pMin.y, bb.pMax.y, bb.pMin.z, bb.pMax.z };
			std::array<float, 6> bB = bA;
			// measurement voxel = generation voxel (--voxel-size) per the study
			// design: relative sensitivities are resolution-robust, so we avoid
			// paying for a finer µCT resolution here.
			const float mVox = voxelSizeArg;
			g.estimate_local_thickness(mVox, bA, false);
			g.estimate_local_thickness(mVox, bB, true);
			g.estimate_anisotropy(mVox, milDirections, milLines, 0);   // ratio (Ulrich)
			g.estimate_trabecular_number(mVox, 0, milDirections, milLines);
			g.estimate_connectivity_density_voxel(mVox, 6);
			g.estimate_smi();
			fout << swept << "," << th << "," << openness << "," << sy << "," << spread << ","
				<< radius << "," << s.size() << ","
				<< g.porosity << "," << (1.0f - g.porosity) << ","
				<< g.localThickness << "," << g.localThicknessStd << ","
				<< g.localSeparation << "," << g.localSeparationStd << ","
				<< g.trabecularNr << "," << g.anisotropyDegree << ","
				<< g.connectivityDensity << "," << g.smi << ","
				<< g.surfaceToVolume << "," << g.surfaceToTotalVolume << "\n";
			fout.flush();
			std::cout << "  [" << swept << "] th=" << th << " op=" << openness << " sy=" << sy
				<< " sp=" << spread << " r=" << radius << " -> por=" << g.porosity
				<< " TbTh=" << g.localThickness << " DA=" << g.anisotropyDegree << std::endl;
		};

		const char* header =
			"swept,thickness,openness,stretch_y,spread,radius,seeds,"
			"porosity,BVTV,TbTh,TbTh_SD,TbSp,TbSp_SD,TbN,DA,ConnD,SMI,BS_BV,BS_TV\n";

		// apply the named input to the baseline, leaving the rest at baseline
		auto apply = [&](const std::string& name, float v,
			float& th, float& op, float& sy, float& sp, float& rad) {
			if (name == "thickness") th = v;
			else if (name == "openness") op = v;
			else if (name == "stretch") sy = v;
			else if (name == "spread") sp = v;
			else if (name == "radius") rad = v;
			else { std::cerr << "unknown sweep input '" << name << "'\n"; }
		};

		if (sweepInput) {
			std::ofstream fout(sweepOut + "/sweep_" + sweepName + ".csv");
			fout << header;
			std::cout << "\n=== SENSITIVITY OAT: " << sweepName << " in ["
				<< sweepA << "," << sweepB << "], " << sweepN << " steps ===" << std::endl;
			for (int i = 0; i < sweepN; ++i) {
				float v = (sweepN <= 1) ? sweepA
					: sweepA + (sweepB - sweepA) * i / static_cast<float>(sweepN - 1);
				float th = th0, op = op0, sy = sy0, sp = sp0, rad = rad0;
				apply(sweepName, v, th, op, sy, sp, rad);
				measure_and_write(fout, sweepName, th, op, sy, sp, rad);
			}
			std::cout << "Wrote " << sweepOut << "/sweep_" << sweepName << ".csv" << std::endl;
		}

		if (sweepGrid) {
			std::ofstream fout(sweepOut + "/sweep_grid_" + gridName1 + "_" + gridName2 + ".csv");
			fout << header;
			std::cout << "\n=== SENSITIVITY GRID: " << gridName1 << " x " << gridName2
				<< " (" << gridN1 << "x" << gridN2 << ") ===" << std::endl;
			for (int i = 0; i < gridN1; ++i) {
				float v1 = (gridN1 <= 1) ? gridA1
					: gridA1 + (gridB1 - gridA1) * i / static_cast<float>(gridN1 - 1);
				for (int j = 0; j < gridN2; ++j) {
					float v2 = (gridN2 <= 1) ? gridA2
						: gridA2 + (gridB2 - gridA2) * j / static_cast<float>(gridN2 - 1);
					float th = th0, op = op0, sy = sy0, sp = sp0, rad = rad0;
					apply(gridName1, v1, th, op, sy, sp, rad);
					apply(gridName2, v2, th, op, sy, sp, rad);
					measure_and_write(fout, gridName1 + "x" + gridName2, th, op, sy, sp, rad);
				}
			}
			std::cout << "Wrote " << sweepOut << "/sweep_grid_" << gridName1 << "_" << gridName2 << ".csv" << std::endl;
		}

		delete container;
		return 0;
	}

	// --- 12.3c Runtime / performance ---------------------------------------
	// Two sweeps that establish practicality: generation time vs (A) grid
	// resolution at a fixed seed cloud, and (B) seed count at a fixed grid.
	// Each configuration is timed runtimeRepeats times and the MEDIAN reported
	// (robust to OS scheduling noise). Run with a box container, e.g.:
	//   scaffoldProfile --box --box-size 5 --poisson --runtime-suite
	if (runtimeSuite) {
		namespace ch = std::chrono;
		std::filesystem::create_directories(runtimeOut);

		int threads = 1;
#ifdef _OPENMP
		threads = omp_get_max_threads();
#endif
		const float boxSide = (boxSize > 0.0f) ? boxSize : 5.0f;
		std::cout << "\n=== RUNTIME SUITE (threads=" << threads << ", box " << boxSide
			<< " mm, repeats=" << runtimeRepeats << ") ===" << std::endl;

		auto median = [](std::vector<double> v) {
			std::sort(v.begin(), v.end());
			return v.empty() ? 0.0 : v[v.size() / 2];
		};
		// Peak resident memory (working set) in MB. Cumulative process peak, so
		// with the resolution sweep run small->large grid each row reflects that
		// grid's footprint; the finest-grid row is the headline number.
		auto peak_ram_mb = []() -> double {
#ifdef _WIN32
			PROCESS_MEMORY_COUNTERS pmc;
			if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
				return pmc.PeakWorkingSetSize / (1024.0 * 1024.0);
#endif
			return 0.0;
		};
		auto res_for = [&](float v) {
			int n = std::max(2, (int)(boxSide / v));
			return std::array<int, 3>{ n, n, n };
		};
		// Time field + marching cubes for one grid/seed pair.
		auto time_run = [&](const std::array<int, 3>& r, const std::vector<Vec3>& s,
			double& fieldS, double& mcS) -> bool {
			GeneratorLewiner g(s, bounds, r, &logger, opennessArg, isoLevel, false);
			if (stretchSet) g.set_stretch(stretchX, stretchY, stretchZ);
			auto ta = ch::steady_clock::now();
			if (!g.compute_scalar_field(*container)) return false;
			auto tb = ch::steady_clock::now();
			if (!g.marching_cubes()) return false;
			auto tc = ch::steady_clock::now();
			fieldS = ch::duration<double>(tb - ta).count();
			mcS = ch::duration<double>(tc - tb).count();
			return true;
		};

		// ---- Sweep A: runtime vs grid resolution (fixed seed cloud) ----
		std::vector<Vec3> fixedSeeds = make_seeds(1u, false);
		std::vector<float> voxelList = { 0.05f, 0.04f, 0.03f, 0.025f, 0.02f, 0.015f, 0.01f };
		std::ofstream fa(runtimeOut + "/runtime_vs_resolution.csv");
		fa << "voxel_mm,nx,ny,nz,voxels,seeds,field_s,mc_s,total_s,peak_ram_mb,threads\n";
		std::cout << "\n-- Sweep A: runtime vs grid resolution (seeds=" << fixedSeeds.size()
			<< ") --\nvoxel,voxels,field_s,mc_s,total_s,peak_ram_mb\n";
		for (float v : voxelList) {
			std::array<int, 3> r = res_for(v);
			std::vector<double> fs, ms;
			for (int rep = 0; rep < runtimeRepeats; ++rep) {
				double f = 0, m = 0;
				if (time_run(r, fixedSeeds, f, m)) { fs.push_back(f); ms.push_back(m); }
			}
			if (fs.empty()) continue;
			double f = median(fs), m = median(ms);
			double ram = peak_ram_mb();
			size_t nv = (size_t)r[0] * r[1] * r[2];
			fa << v << "," << r[0] << "," << r[1] << "," << r[2] << "," << nv << ","
				<< fixedSeeds.size() << "," << f << "," << m << "," << (f + m) << ","
				<< ram << "," << threads << "\n";
			fa.flush();
			std::cout << v << "," << nv << "," << f << "," << m << "," << (f + m) << ","
				<< ram << std::endl;
		}

		// ---- Sweep B: runtime vs seed count (fixed grid), timing seed gen too ----
		const float vFixed = 0.025f;
		std::array<int, 3> rFixed = res_for(vFixed);
		// Keep the Poisson radius >= 2x the wall thickness (isoLevel) so every
		// structure stays an OPEN lattice. Below ~2x, adjacent walls merge into a
		// near-solid: the iso-surface area collapses and marching-cubes time
		// becomes non-monotonic in seed count (a structural artifact, not a
		// scaling law). Expressing radii as multiples of thickness enforces this
		// for any --thickness. Smaller multiples -> denser sampling -> more seeds.
		std::vector<float> radFactors = { 2.0f, 2.5f, 3.0f, 4.0f, 5.0f, 7.0f };
		std::vector<float> radii;
		for (float k : radFactors) radii.push_back(k * isoLevel);
		std::ofstream fb(runtimeOut + "/runtime_vs_seeds.csv");
		fb << "poisson_radius_mm,seeds,voxels,seed_gen_s,field_s,mc_s,total_s,"
			<< "surface_mm2,porosity,threads\n";
		std::cout << "\n-- Sweep B: runtime vs seed count (voxel " << vFixed << ", grid "
			<< rFixed[0] << "^3) --\nradius,seeds,seed_gen_s,field_s,mc_s,total_s,surface_mm2,porosity\n";
		for (float rad : radii) {
			std::vector<double> gs, fs, ms;
			std::vector<Vec3> s;
			for (int rep = 0; rep < runtimeRepeats; ++rep) {
				Poisson3D poisson(rad, rad, 30, false);
				poisson.set_rng_seed(1u);
				auto tg0 = ch::steady_clock::now();
				poisson.run(*container);
				auto tg1 = ch::steady_clock::now();
				s = poisson.get_seeds();
				gs.push_back(ch::duration<double>(tg1 - tg0).count());
				double f = 0, m = 0;
				if (time_run(rFixed, s, f, m)) { fs.push_back(f); ms.push_back(m); }
			}
			if (fs.empty()) continue;
			double g = median(gs), f = median(fs), m = median(ms);
			size_t nv = (size_t)rFixed[0] * rFixed[1] * rFixed[2];
			// One untimed generation to record the structure (surface area,
			// porosity) so the marching-cubes cost can be read against the
			// iso-surface area it tracks.
			float surf = 0.0f, por = 0.0f;
			{
				GeneratorLewiner gm(s, bounds, rFixed, &logger, opennessArg, isoLevel, false);
				if (stretchSet) gm.set_stretch(stretchX, stretchY, stretchZ);
				if (gm.compute_scalar_field(*container) && gm.marching_cubes()) {
					gm.estimate_metrics(*container);
					surf = gm.surfaceArea;
					por = gm.porosity;
				}
			}
			fb << rad << "," << s.size() << "," << nv << "," << g << "," << f << "," << m << ","
				<< (g + f + m) << "," << surf << "," << por << "," << threads << "\n";
			fb.flush();
			std::cout << rad << "," << s.size() << "," << g << "," << f << "," << m << ","
				<< (g + f + m) << "," << surf << "," << por << std::endl;
		}

		std::cout << "\nWrote " << runtimeOut << "/runtime_vs_resolution.csv and runtime_vs_seeds.csv" << std::endl;
		delete container;
		return 0;
	}

	std::vector<Vec3> seeds = make_seeds(seedArg, true);
	GeneratorLewiner gen(seeds, bounds, res, &logger, opennessArg, isoLevel, false);

	if (stretchSet) {
		gen.set_stretch(stretchX, stretchY, stretchZ);
		std::cout << "Stretch: (" << stretchX << ", " << stretchY << ", " << stretchZ << ")" << std::endl;
	}

	if (spreadArg >= 0.0f) {
		gen.spread = spreadArg;
		std::cout << "Fenestration spread: " << spreadArg << std::endl;
	}

	if (roundEdgesArg) {
		gen.roundEdges = true;
		gen.edgeK = edgeKArg;
		std::cout << "Edge rounding: ON (edgeK = " << edgeKArg << ")" << std::endl;
	}

	// --- Auto-calibrate iso-level so measured Tb.Th hits a target ---
	if (calibrateTargetMode) {
		std::cout << "\n=== THICKNESS CALIBRATION (target Tb.Th=" << calTarget
			<< " mm, measure voxel=" << calTargetVoxel
			<< ", openness tau=" << opennessArg << ") ===" << std::endl;
		bool ok = gen.calibrate_thickness(*container, calTarget, calTargetVoxel);
		std::cout << "\nCalibration " << (ok ? "converged" : "FAILED") << std::endl;
		if (ok) gen.marching_cubes(true);   // calibrate_thickness leaves only the field
		gen.estimate_metrics(*container);
		std::cout << "final Tb.Th = " << gen.localThickness << " +- " << gen.localThicknessStd
			<< " mm, porosity = " << gen.porosity << std::endl;
		delete container;
		return 0;
	}

	// --- Joint calibration (2- and 3-knob) ---------------------------------
	if (twoKnobMode || threeKnobMode) {
		bool ok;
		if (threeKnobMode) {
			std::cout << "\n=== THREE-KNOB (Tb.Th=" << jkTh << ", porosity=" << jkPor
				<< ", SMI=" << jkSMI << ", voxel=" << jkVoxel << ") ===" << std::endl;
			ok = gen.three_knob_calibration(*container, jkTh, jkPor, jkSMI, jkVoxel);
		}
		else {
			std::cout << "\n=== TWO-KNOB (Tb.Th=" << jkTh << ", porosity=" << jkPor
				<< ", voxel=" << jkVoxel << ") ===" << std::endl;
			ok = gen.two_knob_calibration(*container, jkTh, jkPor, jkVoxel, jkTol, jkIter);
		}
		std::cout << "\nCalibration " << (ok ? "completed" : "FAILED")
			<< " (compare achieved vs target; a knob pinned at a bound = target "
			<< "outside its reachable band)" << std::endl;
		std::cout << "knobs:  iso-level=" << gen.get_iso_level()
			<< "  openness=" << gen.get_openness()
			<< "  spread=" << gen.get_spread() << std::endl;
		std::cout << "achieved:  Tb.Th=" << gen.localThickness
			<< "  porosity=" << gen.porosity
			<< "  SMI=" << gen.smi << std::endl;

		// Full outcome-metric report at the calibration voxel, so the calibrated
		// scaffold's Tb.N / Conn.D / BS-TV / DA (not tuned) can be read off. Used
		// by the k_edge investigation, which calibrates Tb.Th+porosity at each
		// edgeK and compares these outcomes against the reference.
		if (ok) {
			gen.marching_cubes(true);
			Aabb bb = gen.get_aabb();
			std::array<float, 6> bA = { bb.pMin.x, bb.pMax.x, bb.pMin.y, bb.pMax.y, bb.pMin.z, bb.pMax.z };
			std::array<float, 6> bB = bA;
			gen.estimate_local_thickness(jkVoxel, bA, false);
			gen.estimate_local_thickness(jkVoxel, bB, true);
			gen.estimate_anisotropy(jkVoxel, milDirections, milLines, 0);
			gen.estimate_trabecular_number(jkVoxel, 0, milDirections, milLines);
			gen.estimate_connectivity_density_voxel(jkVoxel, 6);
			gen.estimate_smi();
			std::cout << "KEDGE_ROW," << edgeKArg << "," << (ok ? 1 : 0) << ","
				<< gen.get_iso_level() << "," << gen.localThickness << "," << gen.porosity << ","
				<< gen.localSeparation << "," << gen.trabecularNr << "," << gen.anisotropyDegree << ","
				<< gen.connectivityDensity << "," << gen.smi << ","
				<< gen.surfaceToVolume << "," << gen.surfaceToTotalVolume << std::endl;
		}
		delete container;
		return 0;
	}

	// --- BoneJ cross-tool comparison cohort --------------------------------
	// Emits, per scaffold: an STL mesh, a binary NRRD image, one row of
	// generation parameters and one row of our metrics.
	//
	// The NRRD is produced by the SAME get_image_field() call at the SAME voxel
	// size that our metrics are measured on, so BoneJ analyses bit-for-bit the
	// image we analysed. Any disagreement is therefore algorithmic, not a
	// resolution or resampling difference - which is the entire point of the
	// experiment. Do not change the voxel size between the export and the
	// measurement.
	if (bonejSuite) {
		std::filesystem::create_directories(bonejOut);

		std::cout << "\n=== BoneJ COMPARISON COHORT (n=" << bonejCount
			<< ", measurement voxel " << bonejVoxel << " mm) ===" << std::endl;
		std::cout << "output: " << bonejOut << std::endl;

		std::ofstream fparams(bonejOut + "/generation_parameters.csv");
		std::ofstream fmetrics(bonejOut + "/scaffoldbuilder_metrics.csv");
		if (!fparams.is_open() || !fmetrics.is_open()) {
			std::cerr << "cannot open output CSVs in " << bonejOut << std::endl;
			delete container; return 1;
		}

		fparams << "id,seed,box_mm,thickness_mm,poisson_radius_mm,openness,spread,"
			<< "stretch_x,stretch_y,stretch_z,generation_voxel_mm,measurement_voxel_mm,seed_count\n";
		fmetrics << "id,TbTh,TbTh_SD,TbSp,TbSp_SD,TbN,DA,ConnD,SMI,BS_BV,BS_TV,porosity,"
			<< "volume_mm3,surface_mm2\n";

		// Parameter sampling is itself seeded, so the whole cohort is reproducible.
		std::mt19937 prng(20260715u);
		std::uniform_real_distribution<float> dTh(bonejThMin, bonejThMax);
		std::uniform_real_distribution<float> dSp(bonejSpMin, bonejSpMax);
		std::uniform_real_distribution<float> dStretch(bonejStretchMin, bonejStretchMax);

		const float boxSide = (boxSize > 0.0f) ? boxSize : 4.0f;

		for (int id = 1; id <= bonejCount; ++id) {
			// sample thickness and spacing; spacing must comfortably exceed
			// thickness or the struts merge and the structure is not a lattice
			float th = dTh(prng);
			float sp = dSp(prng);
			int guard = 0;
			while (sp < 2.5f * th && guard++ < 100) sp = dSp(prng);
			if (sp < 2.5f * th) { std::cerr << "id " << id << ": no valid spacing, skipping\n"; continue; }

			Poisson3D poisson(sp, sp, 30, false);
			poisson.set_rng_seed(static_cast<uint32_t>(id));
			poisson.run(*container);
			std::vector<Vec3> s = poisson.get_seeds();
			if (s.size() < 3) { std::cerr << "id " << id << ": too few seeds, skipping\n"; continue; }

			GeneratorLewiner g(s, bounds, res, &logger, opennessArg, th, false);
			g.container = containerShared;   // enable the outside-container mask
			// Anisotropy: a fixed --stretch (if given) applies to all scaffolds;
			// otherwise sample a per-scaffold dominant-y stretch so the cohort
			// spans a real DA range for the anisotropy comparison.
			if (stretchSet) g.set_stretch(stretchX, stretchY, stretchZ);
			else g.set_stretch(1.0f, dStretch(prng), 1.0f);
			if (spreadArg >= 0.0f) g.spread = spreadArg;

			if (!g.compute_scalar_field(*container) || !g.marching_cubes()) {
				std::cerr << "id " << id << ": generation failed, skipping\n";
				continue;
			}

			g.estimate_metrics(*container);
			Aabb bb = g.get_aabb();
			std::array<float, 6> bnds = { bb.pMin.x, bb.pMax.x, bb.pMin.y, bb.pMax.y, bb.pMin.z, bb.pMax.z };
			std::array<float, 6> bA = bnds, bB = bnds;
			g.estimate_local_thickness(bonejVoxel, bA, false);
			g.estimate_local_thickness(bonejVoxel, bB, true);
			// mode 3 = 1 - lambda_min/lambda_max, matching BoneJ's DA definition
			// directly so no post-hoc formula conversion is needed.
			g.estimate_anisotropy(bonejVoxel, milDirections, milLines, 3);
			g.estimate_trabecular_number(bonejVoxel, 0, milDirections, milLines);
			g.estimate_connectivity_density_voxel(bonejVoxel, 6);
			g.estimate_smi();

			char stem[64];
			std::snprintf(stem, sizeof(stem), "/scaffold_%02d", id);
			// g.export_stl(bonejOut + stem + ".stl");
			// Export the .nrrd over the FULL container box, not the mesh AABB:
			// TV (total sample volume) must equal the container so that BoneJ's
			// TV-normalized indices (porosity, BS/TV, BV/TV) match ours, which
			// normalize by domainVolume. Intensive metrics (Tb.Th, DA, ...) are
			// unaffected since the extra air border holds no solid voxels.
			g.export_nrrd(bonejOut + stem + ".nrrd", bonejVoxel, bounds);

			fparams << id << "," << id << "," << boxSide << "," << th << "," << sp << ","
				<< opennessArg << "," << g.spread << ","
				<< g.stretchX << "," << g.stretchY << "," << g.stretchZ << ","
				<< voxelSizeArg << "," << bonejVoxel << "," << s.size() << "\n";
			fparams.flush();

			fmetrics << id << "," << g.localThickness << "," << g.localThicknessStd << ","
				<< g.localSeparation << "," << g.localSeparationStd << ","
				<< g.trabecularNr << "," << g.anisotropyDegree << ","
				<< g.connectivityDensity << "," << g.smi << ","
				<< g.surfaceToVolume << "," << g.surfaceToTotalVolume << ","
				<< g.porosity << "," << g.volume << "," << g.surfaceArea << "\n";
			fmetrics.flush();

			std::cout << "  [" << id << "/" << bonejCount << "] th=" << th << " sp=" << sp
				<< " -> Tb.Th=" << g.localThickness << " Tb.Sp=" << g.localSeparation
				<< " por=" << g.porosity << std::endl;
		}

		std::cout << "\nWrote generation_parameters.csv, scaffoldbuilder_metrics.csv,"
			<< " scaffold_NN.stl and scaffold_NN.nrrd to " << bonejOut << std::endl;
		delete container;
		return 0;
	}

	// --- 12.2b Resolution convergence -------------------------------------
	// Generate the scaffold ONCE, then re-measure it at a range of image voxel
	// sizes. Regenerating per step would confound generation randomness with
	// measurement resolution; reusing one geometry isolates the measurement.
	//
	// Expect: Tb.Th / Tb.Sp / Tb.N / Conn.D vary with voxel size and converge as
	// v -> 0. SMI, BS/BV, BS/TV and porosity are computed on the mesh or the
	// generation grid, so they are INDEPENDENT of the image voxel and should stay
	// flat across the sweep - they act as controls.
	if (sweepVoxel) {
		std::cout << "\n=== RESOLUTION CONVERGENCE (generate once, re-measure) ===" << std::endl;
		std::cout << "generation voxel = " << voxelSizeArg
			<< " mm, MIL " << milDirections << " dirs x " << milLines << " lines" << std::endl;

		if (!gen.compute_scalar_field(*container) || !gen.marching_cubes()) {
			std::cerr << "generation failed, aborting." << std::endl;
			delete container; return 1;
		}
		gen.estimate_metrics(*container);
		gen.estimate_smi();
		const float smiFixed = gen.smi;
		const float bsbvFixed = gen.surfaceToVolume;
		const float bstvFixed = gen.surfaceToTotalVolume;
		const float porFixed = gen.porosity;

		// assemble the voxel sizes: the range, plus any explicit extras (the uCT
		// resolution must be measured, not interpolated), sorted and de-duplicated
		std::vector<float> vlist;
		for (float v = sweepMin; v <= sweepMax + 1e-6f; v += sweepStep) vlist.push_back(v);
		for (float v : sweepExtra) vlist.push_back(v);
		std::sort(vlist.begin(), vlist.end());
		vlist.erase(std::unique(vlist.begin(), vlist.end(),
			[](float a, float b) { return std::fabs(a - b) < 1e-6f; }), vlist.end());

		std::cout << "voxel,TbTh,TbTh_SD,TbSp,TbSp_SD,TbN,DA,ConnD,SMI,BS/BV,BS/TV,porosity" << std::endl;

		Aabb bb = gen.get_aabb();
		for (float v : vlist) {
			std::array<float, 6> bA = { bb.pMin.x, bb.pMax.x, bb.pMin.y, bb.pMax.y, bb.pMin.z, bb.pMax.z };
			std::array<float, 6> bB = bA;
			gen.estimate_local_thickness(v, bA, false);
			gen.estimate_local_thickness(v, bB, true);
			gen.estimate_anisotropy(v, milDirections, milLines, 0);
			gen.estimate_trabecular_number(v, 0, milDirections, milLines);
			gen.estimate_connectivity_density_voxel(v, 6);

			std::cout << v << "," << gen.localThickness << "," << gen.localThicknessStd << ","
				<< gen.localSeparation << "," << gen.localSeparationStd << ","
				<< gen.trabecularNr << "," << gen.anisotropyDegree << ","
				<< gen.connectivityDensity << ","
				<< smiFixed << "," << bsbvFixed << "," << bstvFixed << "," << porFixed
				<< std::endl;
		}
		delete container;
		return 0;
	}

	// --- One-shot: generate and report all voxel-based metrics at one resolution ---
	if (metricsMode) {
		std::cout << "\n=== METRICS (voxel size " << metricsVoxel
			<< ", openness tau=" << opennessArg << ") ===" << std::endl;
		if (!gen.compute_scalar_field(*container) || !gen.marching_cubes()) {
			std::cerr << "generation failed, aborting." << std::endl;
			delete container; return 1;
		}
		gen.estimate_metrics(*container);
		Aabb bb = gen.get_aabb();
		std::array<float, 6> bA = { bb.pMin.x, bb.pMax.x, bb.pMin.y, bb.pMax.y, bb.pMin.z, bb.pMax.z };
		std::array<float, 6> bB = bA;
		gen.estimate_local_thickness(metricsVoxel, bA, false);
		gen.estimate_local_thickness(metricsVoxel, bB, true);
		gen.estimate_anisotropy(metricsVoxel, milDirections, milLines, 0);      // mode 0: Lmax/Lmin
		gen.estimate_trabecular_number(metricsVoxel, 0, milDirections, milLines); // formula 0: MIL
		// both connectivity methods write the same member, so only call the
		// preferred (voxel Euler) one here
		gen.estimate_connectivity_density_voxel(metricsVoxel, 6);
		gen.estimate_smi();                                         // rod(3) vs plate(0)
		std::cout << "Tb.Th = " << gen.localThickness << " +- " << gen.localThicknessStd << " mm" << std::endl;
		std::cout << "Tb.Sp = " << gen.localSeparation << " +- " << gen.localSeparationStd << " mm" << std::endl;
		std::cout << "Tb.N  = " << gen.trabecularNr << " 1/mm" << std::endl;
		std::cout << "DA    = " << gen.anisotropyDegree << std::endl;
		std::cout << "SMI   = " << gen.smi << " (0=plate, 3=rod)" << std::endl;
		std::cout << "BS/BV = " << gen.surfaceToVolume << " 1/mm" << std::endl;
		std::cout << "BS/TV = " << gen.surfaceToTotalVolume << " 1/mm" << std::endl;
		std::cout << "Conn.D  = " << gen.connectivityDensity << " 1/mm^3" << std::endl;
		std::cout << "porosity = " << gen.porosity << " (BV/TV = " << (1.0f - gen.porosity) << ")" << std::endl;
		delete container; return 0;
	}

	// --- Calibration sweep: iso-level (thickness) -> measured Tb.Th ---
	// Uniform thickness only; run with --box so the whole domain is fillable.
	if (calibrateMode) {
		std::cout << "\n=== CALIBRATION SWEEP (openness tau=" << opennessArg
			<< ", voxel size=" << calVoxel << ") ===" << std::endl;
		std::cout << "isoLevel_mm,TbTh_mm,TbTh_std_mm,porosity" << std::endl;
		for (int s = 0; s < calSteps; ++s) {
			float t = (calSteps > 1)
				? calStart + (calEnd - calStart) * s / (calSteps - 1)
				: calStart;
			gen.set_thickness(t);
			if (!gen.compute_scalar_field(*container)) continue;
			if (!gen.marching_cubes()) continue;
			Aabb bb = gen.get_aabb();
			std::array<float, 6> blockBounds = {
				bb.pMin.x, bb.pMax.x, bb.pMin.y, bb.pMax.y, bb.pMin.z, bb.pMax.z
			};
			gen.estimate_local_thickness(calVoxel, blockBounds, false);
			std::cout << t << "," << gen.localThickness << ","
				<< gen.localThicknessStd << "," << gen.porosity << std::endl;
		}
		delete container;
		return 0;
	}

	if (reproBug) {
		gen.set_stretch(1.0f, stretchYArg, 1.0f);
		std::cout << "Varied thickness (repro): start=" << startThicknessArg
			<< " end=" << endThicknessArg << " transitionDist=" << transitionArg
			<< " stretchY=" << stretchYArg << std::endl;
		gen.set_thickness_functions(
			container->sdf,
			std::make_shared<LinearFunction>(transitionArg),
			startThicknessArg, endThicknessArg, transitionArg
		);
	}
	else if (!noThickness && !useBox) {
		float startThickness = 10.0f;
		float endThickness = isoLevel;
		float transitionDistance = 2.0f;
		std::cout << "Varied thickness ENABLED: start=" << startThickness
			<< " end=" << endThickness << " transitionDist=" << transitionDistance << std::endl;
		gen.set_thickness_functions(
			container->sdf,
			std::make_shared<LinearFunction>(transitionDistance),
			startThickness, endThickness, transitionDistance
		);
	}
	else {
		std::cout << "Varied thickness DISABLED" << std::endl;
	}

	// Build anisotropy sources. update_model() is intentionally NOT called here
	// because there is no OpenGL context in the profiler. update_metric() is
	// pure math and is safe; compute_scalar_field also calls it at startup.
	if (!anisoSourceArgs.empty()) {
		std::cout << "Anisotropy sources (" << anisoSourceArgs.size()
			<< "), background weight=" << backgroundWeightArg << ":" << std::endl;
		for (size_t i = 0; i < anisoSourceArgs.size(); i++) {
			const auto& s = anisoSourceArgs[i];
			auto src = std::make_shared<AnisotropySource>();
			src->origin  = Vec3(s[0], s[1], s[2]);
			src->stretch = Vec3(s[3], s[4], s[5]);
			src->sigma   = s[6];
			src->update_metric();
			gen.anisotropySources.push_back(src);
			std::cout << "  [" << i << "] origin=(" << s[0] << "," << s[1] << "," << s[2]
				<< ") stretch=(" << s[3] << "," << s[4] << "," << s[5]
				<< ") sigma=" << s[6] << std::endl;
		}
		gen.backgroundWeight = backgroundWeightArg;
	}
	else {
		std::cout << "No anisotropy sources (isotropic run)" << std::endl;
	}

	std::cout << "\n--- compute_scalar_field ---" << std::endl;
	auto t0 = std::chrono::steady_clock::now();
	if (!gen.compute_scalar_field(*container)) {
		std::cerr << "compute_scalar_field failed (too few seeds?), aborting." << std::endl;
		return 1;
	}
	auto t1 = std::chrono::steady_clock::now();

	std::cout << "\n--- marching_cubes ---" << std::endl;
	if (!gen.marching_cubes()) {
		std::cerr << "marching_cubes failed (invalid scalar field), aborting." << std::endl;
		return 1;
	}
	auto t2 = std::chrono::steady_clock::now();

	std::cout << "\n=== SUMMARY ===" << std::endl;
	std::cout << "compute_scalar_field TOTAL: " << std::chrono::duration<double>(t1 - t0).count() << " s" << std::endl;
	std::cout << "marching_cubes TOTAL:       " << std::chrono::duration<double>(t2 - t1).count() << " s" << std::endl;

	delete container;
	return 0;
}
