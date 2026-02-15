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
#include <atomic>
#include <thread>
#include <functional>

// custom
#include <OpenGlRender/Model.h>
#include <SeedGenerator/Container.h>
#include <SeedGenerator/SeedGenerator.h>
#include <Visualize/VisualizeSeeds.h>
//#include "ScaffoldGenerator/ScaffoldGenerator.h"
#include "SeedGenerator/DistanceCalculator.h"
#include "OpenGlSetup/trackBallCamera.h"
#include "OpenGlSetup/defaultCamera.h"
#include "OpenGlSetup/simpleCamera.h"
#include "OpenGlSetup/UniformManager.h"
#include "OpenGlSetup/Grid.h"
#include "ScaffoldGenerator/Generator.h"
#include "Shader.h"
#include "Utils/Utils.h"
//#include "SeedGenerator/Poisson3D.h"
#include "Logger/Logger.h"
#include "OpenGlSetup/Misc.h"

struct poreNetworkObject {
	std::unique_ptr<PoreNetwork> model;
	std::array<float, 4> lineColor{ 0.0f, 0.0f, 0.0f, 1.0f };
	std::string name = "PoreNetwork_";
	float lineSize = 1.0f;
};

//struct scaffoldObject {
//	std::unique_ptr<Model> model;
//	std::array<float, 4> objectColor{1.0f, 0.5f, 0.5f, 1.0f};
//	std::string name = "Scaffold_";
//};

struct RenderSettings {
	std::array<float, 4> poreNetworkColor = { 0.0f, 0.0f, 0.0f, 1.0f };
	float poreNetworkLineSize = 1.0f;
};

struct ContainerModel {
	void* obj = nullptr;
	bool show = true;
	std::array<float, 3> center = { 0.0, 0.0, 0.0 };
};

//// seed generator properties
//enum SeedGeneratorType {
//	Box, Cylinder, Mesh, None
//};
//
//struct SeedGeneratorObject {
//	std::string name = "";
//	SeedGeneratorType type = SeedGeneratorType::None;
//	std::unique_ptr<IContainer> container;
//	std::unique_ptr<InterfaceSeedGenerator> generator;
//};

//struct SelectedObject {
//	ObjectType type = ObjectType::None;
//	void* ptr = nullptr;
//};

struct SelectedObject {
	void* ptr = nullptr;
	ObjectType type = ObjectType::NoneType;
};

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
	float lightColor[3] = { 1.0f, 1.0f, 1.0f };
	float lightPosCamera[3] = { 0.0f, 5.0f, 10.0f };
	float Ka{ 0.2f };
	float Ks{ 2.0f };

	// selection objects
	SelectedObject selectedPanelObj;
	void* selectedSceneObj = nullptr;

	// mesh objects
	std::vector<std::unique_ptr<Generator>> scaffolds;

	std::vector<std::unique_ptr<InterfaceSeedGenerator>> seedGenerators;

	std::vector<std::unique_ptr<IContainer>> containers;

	//std::unique_ptr<Generator> scaffoldGenerator;

	std::vector<poreNetworkObject> poreNetworkList;
	int activePoreNetworkIndex{ -1 };

	std::unique_ptr<VisualizeSeeds> seedObj;
	std::unique_ptr<CutPlane> cutPlane;
	std::unique_ptr<BBox> box;
	std::unique_ptr<Grid> grid;
	std::unique_ptr<CutPlane> distPlane;
	std::unique_ptr<Arrow> xArrow;
	std::unique_ptr<Arrow> yArrow;
	std::unique_ptr<Arrow> zArrow;
	std::unique_ptr<PoreNetwork> poreNetwork;

	// containers
	CylinderContainer cylContainer;
	BoxContainer boxContainer;
	//MeshContainer meshContainer;
	ContainerModel container;

	// create the active container and initialize it to a box container
	IContainer* activeContainer = nullptr;
	std::array<float, 6> bounds = {};

	//PlaneDistEstimator planeDistance;
	//MeshDistEstimator meshDistance;
	//BoxDistEstimator boxDistance;
	//CylinderDistEstimator cylinderDistance;
	//ImplicitFunctionDistEstimator containerDistance;

	// gui split viewports
	//float split{ 0.0f };

	// flags for scaffold mesh
	bool scaffoldReady{ false };
	bool scaffoldCreated{ false };

	// filenames
	std::string scaffoldFilePath{ "" };
	std::string scaffoldFileName{ "" };
	std::string boneFileName{ "" };

	std::array<float, 6> voxelBounds = {};
	float voxelSize = 0.05f;

	// flag for bone mesh
	bool boneReady{ false };

	// flags for view panel
	bool showSeeds{ true };
	bool showNormals{ false };
	bool showEdges{ false };
	bool showGrid{ false };
	bool showScaffold{ true };
	bool showContainer{ false };
	bool showPoreNetwork{ true };
	bool showDistancePlane{ false };
	bool showBinaryImageWindow{ false };
	bool showGeometryExportWindow{ false };
	bool showVisualizer{ true };
	bool showTortuosityPath{ false };
	bool showScaffoldList{ false };
	bool showConsole{ false };	
	bool showExitConfirmation{ false };
	bool showToolbar{ true };
	bool showBoxSeedCreator{ false };
	bool showBoxContainerCreator{ false };
	bool showCylinderContainerCreator{ false };
	bool showAbstractContainerCreator{ false };
	bool showRandomSeedCreator{ false };
	bool showUniformSeedCreator{ false };
	bool showVariedSeedCreator{ false };
	bool showScaffoldCreator{ false };
	bool showProperties{ false };
	bool measureThickness{ false };
	bool measureSeparation{ false };
	
	bool estimatePoreNetwork = { false };
	bool estimateTortuosity = { false };

	bool updateScaffold = { false };

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
	bool scrollToBottom = false;

	// a variable for adding some message to the log
	std::string suffix{ "" };

	// rendering settings
	RenderSettings renderSettings;

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
	bool cameraUpdate{ false };
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
	std::unique_ptr<SeedGeneratorInterface> seedGenerator;

	// -----------------------------------------------------------
	// 3. Scaffold settings
	//std::unique_ptr<GeneratorInterface> generator;
	char version[256]{ "model" };
	float xDim{ 10.0 };
	float yDim{ 10.0 };
	float zDim{ 10.0 };
	int seedNr{ 100 };
	float xMin{ 0 }, xMax{ 0 }, yMin{ 0 }, yMax{ 0 }, zMin{ 0 }, zMax{ 0 };
	float thickness{ 0.3 };
	float scaffoldConnectivity{ 0.0 };
	float connectivityThreshold{ 0.5f };
	int conOption{ 0 };
	int prevConOption{ -1 };
	int generateOption{ 0 };
	int distFunc{ 0 };
	int radiusFunc{ 0 };
	int radiusOpt{ 0 };
	int regSteps{ 0 };
	int nX{ 10 }, nY{ 10 }, nZ{ 10 };
	std::array<int, 3> resolution = { 100, 100, 100 };
	std::array<int, 3> imageResolution = {200, 200, 200};
	float wallResolution{ 1.0f };

	// determine function to check if a point is inside a container
	std::function<bool(const std::array<double, 3>&)> inside_check;

	// model details
	double scaffoldPorosity{ 50.0 };
	double scaffoldVolume{ 0.0 };
	int faceNr{ 0 }, vertexNr{ 0 }, edgeNr{ 0 };

	// double domain volume
	double domainVolume{ 0.0 };

	// container file


	// Cylinder Container


	// Poisson3d
	float rMin{ 0.8f };
	float rMax{ 2.0f };
	LinearFunction linearFunc;
	double maxDist{ 5.0 };
	float distPlaneCenter[3] = { 0.0f, 0.0f, 0.0f };
	float distPlaneNormal[3] = { 0.0f, 0.0f, 1.0f };
	float distPoint[3] = {0.0f, 0.0f, 0.0f};
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

	// thread handling
	//std::atomic<bool> scaffoldGenerating{ false };
	std::atomic<float> scaffoldProgress{ 0.0f };
	std::thread scaffoldGenerationThread;

	//
	std::array<float, 4> logColor = { 1.0f, 1.0f, 1.0f, 1.0f };

	// textures
	GLuint boxContainerTexture = 0;
	GLuint cylinderContainerTexture = 0;
	GLuint abstractContainerTexture = 0;

	GLuint randomSeedTexture = 0;
	GLuint uniformSeedTexture = 0;
	GLuint variedSeedTexture = 0;

	GLuint createScaffoldTexture = 0;
	GLuint localThicknessTexture = 0;
	GLuint separationTexture = 0;
	GLuint tortuosityTexture = 0;
	GLuint poreNetworkTexture = 0;
	GLuint updateScaffoldTexture = 0;

	// ----------------------------------------------------------------------
	// 4. Functions

	void _init_opengl();

	void _init_imgui();

	void _create_shaders();

	std::string _get_glsl_version();

	void _render_toolbar();

	void _render_settings_panel();

	void _render_seed_generator_list();

	void _render_container_list();

	void _render_object_list();

	void _render_properties_panel();

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

	void _render_menu_bar();

	void _update_cameras();

	void _update_cameras(IContainer& container);

	void _update_bounds_center(int& conOption);
	
	void _render_mesh_settings();

	void _render_axes_viewport();

	void _render_main_menu_bar();

	void _render_display_settings();

	void _render_seed_generator();

	void _render_volume_optimization();
	
	void _render_scaffold_settings();

	void _render_box_container_creator(const char* popupName, bool& showPopup);

	void _render_cylinder_container_creator(const char* popupName, bool& showPopup);

	void _render_abstract_container_creator(const char* popupName, bool& showPopup);

	void _render_random_seed_creator(const char* popupName, bool& showPopup);

	void _render_random_seed_generator_properties();

	void _render_uniform_seed_creator(const char* popupName, bool& showPopup);
	
	void _render_uniform_seed_generator_properties();

	void _render_varied_seed_creator(const char* popupName, bool& showPopup);

	void _render_varied_seed_generator_properties();

	void _render_scaffold_creator(const char* popupName, bool& showPopup);

	void _render_scaffold_properties();

	void _render_binary_image_window(const char* popupName, bool& showPopup);

	void _action_update_voronoi();

	void _action_estimate_connectivity();

	void _local_thickness_measure(const char* popupName, bool& showPopup, bool flag);

	void _action_estimate_tortuosity();
	
	void _action_estimate_pore_network();

	void _action_generate_scaffold();

	void _render_cutting_plane_settings(const char* popupName, bool& showPopup);

	void _create_dockspace();

	void _on_close_request();

	//@brief static function for window close callback
	static void window_close_callback(GLFWwindow* window) {
		myGUI* app = static_cast<myGUI*>(glfwGetWindowUserPointer(window));

		if (app) {
			app->_on_close_request();
		}
	};

	void _render_visualizer();
};

bool create_single_button_textured(const char* name, bool& flag, const std::string& tooltip, const GLuint textureId, bool enabled=true);

// @brief Function to load texture from file, acquired by https://github.com/ocornut/imgui/wiki/Image-Loading-and-Displaying-Examples
bool load_texture_from_file(const char* fileName, GLuint* outTexture, int* outWidth, int* outHeight);

bool load_texture_from_memory(const void* data, size_t dataSize, GLuint* outTexture, int* outWidth, int* outHeight);

//template a function for selectables
template <typename L>
void render_selectable_list(std::vector<std::unique_ptr<L>>& list, L*& selected, const std::string name) {

	for (const auto& md : list) {

		auto obj = md.get();

		if (obj) {

			bool isSelected = (selected && selected == obj);

			if (ImGui::Selectable(md->name.c_str(), isSelected)) {
				selected = md.get();
			};

			if (isSelected) ImGui::SetItemDefaultFocus();
		}
	}
};
#endif