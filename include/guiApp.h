#ifndef GUIAPP_H
#define GUIAPP_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <memory>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <Model.h>
#include <Visualize/VisualizeSeeds.h>
#include "ScaffoldGenerator/ScaffoldGenerator.h"
#include "SeedGenerator/DistanceCalculator.h"
#include "OpenGlSetup/trackBallCamera.h"
#include "OpenGlSetup/defaultCamera.h"
#include "OpenGlSetup/simpleCamera.h"
#include "OpenGlSetup/UniformManager.h"
#include "OpenGlSetup/Grid.h"
#include "Shader.h"
#include "SeedGenerator/Poisson3D.h"
#include "Logger/Logger.h"
#include "OpenGlSetup/Misc.h"

class myGUI {

public:

	ImGuiIO io;

	// render settings
	float scaffoldColor[4] = { 1.0f, 0.5f, 0.5f, 1.0f };
	float containerColor[4] = { 1.0f, 0.5f, 0.5f, 0.5f };
	float seedColor[4] = { 1.0f, 0.0f, 0.0f, 1.0f };
	float gridColor[4] = { 0.5f, 0.5f, 0.5f, 1.0f };
	float meshColor[3] = { 1.0f, 1.0f, 1.0f }; // RGB
	float seedSize{ 0.05f };
	float normalColor[4] = { 0.0f, 1.0f, 0.0f, 1.0f };

	// lighting
	float lightColor[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
	float lightPosCamera[3] = { 0.0f, 5.0f, 10.0f };
	float Ka{ 0.2f };
	float Ks{ 2.0f };

	// mesh objects
	vtkSmartPointer<vtkPolyData> scaffoldPolyData;
	std::unique_ptr<Model> scaffoldModel;
	//Model* containerModel{ nullptr };
	std::unique_ptr<Model> containerModel;
	std::vector<std::array<double, 3>> seeds;
	//VisualizeSeeds* seedObj{ nullptr };
	std::unique_ptr<VisualizeSeeds> seedObj;
	//CutPlane* cutPlane{ nullptr };
	std::unique_ptr<CutPlane> cutPlane;
	//BBox* box;
	std::unique_ptr<BBox> box;
	//BBox* container;
	std::unique_ptr<BBox> container;
	//Grid* grid;
	std::unique_ptr<Grid> grid;
	//CutPlane* distPlane;
	std::unique_ptr<CutPlane> distPlane;
	std::unique_ptr<Arrow> xArrow;
	std::unique_ptr<Arrow> yArrow;
	std::unique_ptr<Arrow> zArrow;

	PlaneDistEstimator planeDistance;
	MeshDistEstimator meshDistance;
	//BoxDistEstimator boxDistance;
	//CylinderDistEstimator cylinderDistance;
	ImplicitFunctionDistEstimator containerDistance;

	// gui split viewports
	//float split{ 0.0f };

	// flags for scaffold mesh
	bool scaffoldReady{ false };
	bool scaffoldCreated{ false };

	// filenames
	std::string scaffoldFilePath{ "" };
	std::string scaffoldFileName{ "" };
	std::string boneFileName{ "" };

	// flag for bone mesh
	bool boneReady{ false };

	// flags for view panel
	bool showSeeds{ false };
	bool showNormals{ false };
	bool showEdges{ false };
	bool showGrid{ false };
	bool showScaffold{ true };
	bool showContainer{ false };
	bool showDistancePlane{ false };

	// flags for tool panel
	bool showCutPlane{ false };
	bool changeCutPlane{ false };
	float planeOrigin[3] = { 0.0f, 0.0f, 0.0f };
	float planeNormal[3] = { 0.0f, 0.0f, 1.0f };
	float tempNormal[3] = { 0.0f, 0.0f, 1.0f };
	float planeOffset = { 0.0f };
	glm::vec4 planeCoeffs{ 1.0f, 0.0f, 0.0f, 0.0f };

	// flag for seeds obj
	bool seedsReady{ false };

	// window flags
	bool showDisplaySettingsWin{ false };
	bool showDisplayMeshSettingsWin{ false };
	bool showPlaneCutSettings{ false };

	// for console logging
	//std::vector<std::string> log;
	Logger& logger = Logger::get_instance();

	// pointers for seed generator
	//Random* sgr;

	// constructor
	myGUI() {};

	myGUI(int width, int height);

	void init();

	void clean();

	void update_render();

	void framebuffer_size_callback_imp(int width, int height);

	void run();

	// destructor
	~myGUI() {};

private:

	// -----------------------------------------------------
	// 1. Opengl window settings
	int height{ 1200 };
	int width{ 800 };

		// background color
	glm::vec4 fontColor = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);

	// framebuffer
	FrameBuffer framebuffer;

	float lineWidth = 0.1f;
	float pointSize = 1.0f;
	//float split{ 0.3 };

	// interpolation size
	double edgeSize{ 2.0 };
	double scaleFactor{ 0.5f };

	// static pointer
	//static myGUI* instance;

	// opengl stuff
	GLFWwindow* window;
	std::string glsl_version;

	// 1a Camera Managment
	// camera options
	bool defCameraFlag{ false };
	bool trackCameraFlag{ true };
	// camera type
	enum cameraType {
		defaultOption, trackOption
	};
	cameraType cameraOption = trackOption;

	glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 5.0f);
	glm::vec3 cameraCenter = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 cameraTarget = glm::vec3(0.0f, 0.0f, 0.0f);

	std::unique_ptr<TrackBall> trackCamera;
	defaultCamera* defCamera;
	std::unique_ptr<SimpleCamera> sCamera;
	ProjectionMode cameraProjectionMode = ProjectionMode::Ortho;

	glm::mat4 projection = glm::mat4(1.0f);
	glm::mat4 view = glm::mat4(1.0f);
	glm::mat4 model = glm::mat4(1.0f);
	glm::mat4 modelPlane = glm::mat4(1.0f);

	// 1b Shader Management
	UniformManager uniManager;

	Shader scaffoldShader;
	Shader normalShader;
	Shader edgeShader;
	Shader gridShader;
	Shader seedShader;
	Shader cutShader;
	Shader bboxShader;
	Shader containerShader;
	Shader frameShader;

	// -----------------------------------------------------------
	// 2. Seed Generator
	//RandomGenerator* randomGenerator;

	// -----------------------------------------------------------
	// 3. Scaffold settings
	char version[256]{ "model" };
	float xDim{ 10.0 };
	float yDim{ 10.0 };
	float zDim{ 10.0 };
	int seedNr{ 100 };
	float xMin{ 0 }, xMax{ 0 }, yMin{ 0 }, yMax{ 0 }, zMin{ 0 }, zMax{ 0 };
	float thickness{ 0.3 };
	int conOption{ 0 };
	int generateOption{ 0 };
	int distFunc{ 0 };
	int radiusFunc{ 0 };
	int radiusOpt{ 0 };
	int regSteps{ 0 };
	int nX{ 10 }, nY{ 10 }, nZ{ 10 };
	std::array<int, 3> resolution = {};
	bool customResolutionFlag{ false };
	float wallResolution{ 1.0f };

	// model details
	double scaffoldPorosity{ 50.0 };
	double scaffoldVolume{ 0.0 };
	int faceNr{ 0 }, vertexNr{ 0 }, edgeNr{ 0 };

	// double domain volume
	double domainVolume{ 0.0 };

	// container file
	std::string containerFile{ "" };
	int neighbors{ 1 };
	vtkSmartPointer<vtkPolyData> containerMesh;
	double containerCenter[3]{ 0.0, 0.0, 0.0 };

	// Cylinder Container
	int cylinderDir{ 0 };
	float cylinderPt[3] = { 0.0f, 0.0f, 0.0f };
	float cylinderAxis[3] = {0.0f, 0.0f, 1.0f};
	double cylinderRadius{ 2.0f };
	double cylinderHeight{ 1.0f };

	// Poisson3d
	float rMin{ 0.8f };
	float rMax{ 2.0f };
	LinearFunction linearFunc;
	double maxDist{ 5.0 };
	float distPlaneCenter[3] = { 0.0f, 0.0f, 0.0f };
	float distPlaneNormal[3] = { 0.0f, 0.0f, 1.0f };
	std::vector<float> radii;

	// Volume Optimization
	Eigen::VectorXd wInit;
	Eigen::VectorXd targetVols;
	int volOption{ 0 };
	bool runVolumeOptimization = false;

	// loading flags
	enum meshType {
		scaffold, bone
	};
	meshType loadedMesh;

	// ----------------------------------------------------------------------
	// 4. Functions

	void _init_opengl();

	void _init_imgui();

	void _create_shaders();

	std::string _get_glsl_version();

	void _render_settings_panel();

	void write_settings();

	static void help_marker(const char* descr) {
		ImGui::TextDisabled("(?)");
		if (ImGui::BeginItemTooltip())
		{
			ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
			ImGui::TextUnformatted(descr);
			ImGui::PopTextWrapPos();
			ImGui::EndTooltip();
		}
	}

	// console stuff
	void _render_console();

	void add_log(LogPriority priority, const std::string& message);

	void _generate_scaffold();

	void _render_scaffold();

	void _render_container_options(int& conOption);

	void _render_menu_bar();

	void _update_cameras();

	void _update_bounds_center(int& conOption);
	
	void _export_mesh();

	void _render_mesh_settings();

	void _render_axes_viewport();

	void _render_main_menu_bar();

	void _render_display_settings();

	void _render_cut_tool();

	void _render_seed_generator();

	void _render_volume_optimization();
	
	void _render_scaffold_settings();

	void _action_generate_seeds();

	void _action_generate_scaffold();

	void _create_dockspace();

	void _render_visualizer();
};

#endif