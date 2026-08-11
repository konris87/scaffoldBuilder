#define _CRT_SECURE_NO_WARNINGS
#define STB_IMAGE_IMPLEMENTATION
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <iostream>
#include <fstream> 
#include <json.hpp>
#include <cmath>
#include <filesystem>
#include <chrono>
#include <sstream>
#include <iomanip> 

// custom
#include "guiApp.h"
#include <ImGuiFileDialog/ImGuiFileDialog.h>
#include "buildScaffold.h"
//#include "ScaffoldGenerator/ScaffoldGenerator.h"
#include "SeedGenerator/DistanceCalculator.h"
#include "OpenGlRender/Model.h"
#include "OpenGlSetup/stb_image.h"
#include "Misc/Imgui_Stdlib.h"
#include "Math/QuadricSimplifier.h"

//#include "SeedGenerator/Random.h"
//#include "SeedGenerator/Poisson3D.h"
#include "Logger/Logger.h"


myGUI::myGUI(int width, int height) : width(width), height(height) {
	
	_init_opengl();

	glsl_version = _get_glsl_version().c_str();

	_init_imgui();
};

void myGUI::_init_opengl() {

	//std::cout << "starting opengl" << std::endl;
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Get the primary monitor
	GLFWmonitor* primary = glfwGetPrimaryMonitor();

	// Get the video mode (resolution, refresh rate, etc.)
	const GLFWvidmode* mode = glfwGetVideoMode(primary);

	glfwWindowHint(GLFW_MAXIMIZED, GLFW_TRUE);

	// Create a window the size of the monitor
	window = glfwCreateWindow(mode->width, mode->height, "Scaffold Builder", NULL, NULL);
	if (!window) {
		glfwTerminate();
		throw std::runtime_error(std::string(std::string("Failed to open GLFW window.") +
			" If you have an Intel GPU, they are not 3.3 compatible." +
			"Try the 2.1 version.\n"));
	}

	width = mode->width;
	height = mode->height;

	// store the pointer
	glfwSetWindowUserPointer(window, this);

	// add the callback
	glfwSetWindowCloseCallback(window, window_close_callback);

	// create context
	glfwMakeContextCurrent(window);

	// glad: load all OpenGL function pointers, these are OS-specific
	// ----------------------------------------------------------------
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		logger.log(LogPriority::ERROR, "Failed to initialize GLAD!");
		//std::cerr << "Failed to initialize GLAD" << std::endl;
		return;
	}

	// Now you can check the OpenGL version
	std::string version = std::string((const char*)glGetString(GL_VERSION));
	std::string glsVersion = std::string((const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));

	// Print the OpenGL version
	logger.log(LogPriority::INFO, "OpenGL Version: " + version);
	logger.log(LogPriority::INFO, "GLSL Version: " + glsVersion);

	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

	// Set the mouse at the center of the screen
	glfwPollEvents();
	glfwSetCursorPos(window, width / 2, height / 2);

	// Gray background color
	glClearColor(fontColor[0], fontColor[1], fontColor[2], fontColor[3]);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);

	// enable blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_CULL_FACE);
	
	// load shaders
	scaffoldShader = Shader(
		"./share/shaders/mainVShader.vertexshader",
		"./share/shaders/mainFShader.fragmentshader",
		NULL
	);

	// add uniforms to the shaders
	uniManager.add_uniform(scaffoldShader, "projection");
	uniManager.add_uniform(scaffoldShader, "view");
	uniManager.add_uniform(scaffoldShader, "model");
	uniManager.add_uniform(scaffoldShader, "cutPlaneCoeffs");
	uniManager.add_uniform(scaffoldShader, "cutPlane");
	uniManager.add_uniform(scaffoldShader, "lightColor");
	uniManager.add_uniform(scaffoldShader, "objectColor");
	uniManager.add_uniform(scaffoldShader, "lightPosWorld");
	uniManager.add_uniform(scaffoldShader, "Ka");
	uniManager.add_uniform(scaffoldShader, "Ks");

	scaffoldShader.use();
	// pass light
	uniManager.setUniform(scaffoldShader, "lightColor", lightColor[0], lightColor[1], lightColor[2]);
	uniManager.setUniform(scaffoldShader,
		"lightPosWorld", lightPosCamera[0], lightPosCamera[1], lightPosCamera[2]);
	uniManager.setUniform(scaffoldShader,
		"objectColor", scaffoldColor[0], scaffoldColor[1], scaffoldColor[2], scaffoldColor[3]);
	uniManager.setUniform(scaffoldShader, "Ka", Ka);
	uniManager.setUniform(scaffoldShader, "Ks", Ks);

	normalShader = Shader(
		"./share/shaders/normalVShader.vs",
		"./share/shaders/normalFShader.fs",
		"./share/shaders/normalGShader.gs"
	);

	uniManager.add_uniform(normalShader, "view");
	uniManager.add_uniform(normalShader, "model");
	uniManager.add_uniform(normalShader, "projection");
	uniManager.add_uniform(normalShader, "normalColor");

	edgeShader = Shader(
		"./share/shaders/edgeVShader.vs",
		"./share/shaders/edgeFShader.fs",
		NULL
	);
	uniManager.add_uniform(edgeShader, "view");
	uniManager.add_uniform(edgeShader, "model");
	uniManager.add_uniform(edgeShader, "projection");
	uniManager.add_uniform(edgeShader, "cutPlaneCoeffs");
	uniManager.add_uniform(edgeShader, "cutPlane");
	uniManager.add_uniform(edgeShader, "edgeColor");

	gridShader = Shader(
		"./share/shaders/gridVShader.vs",
		"./share/shaders/gridFShader.fs",
		"./share/shaders/gridGShader.gs"
	);
	uniManager.add_uniform(gridShader, "view");
	uniManager.add_uniform(gridShader, "model");
	uniManager.add_uniform(gridShader, "projection");
	uniManager.add_uniform(gridShader, "gridColor");

	seedShader = Shader(
		"./share/shaders/seedVShader.vertexshader",
		"./share/shaders/seedFShader.fragmentshader",
		NULL
	);
	uniManager.add_uniform(seedShader, "projection");
	uniManager.add_uniform(seedShader, "view");
	uniManager.add_uniform(seedShader, "seedColor");
	uniManager.add_uniform(seedShader, "seedSize");

	cutShader = Shader(
		"share/shaders/cutVShader.vs",
		"share/shaders/cutFShader.fs",
		NULL
	);
	uniManager.add_uniform(cutShader, "projection");
	uniManager.add_uniform(cutShader, "view");
	uniManager.add_uniform(cutShader, "model");
	uniManager.add_uniform(cutShader, "minBounds");
	uniManager.add_uniform(cutShader, "maxBounds");
	uniManager.add_uniform(cutShader, "planeColor");

	// bounding box shader
	boxShader = Shader(
		"./share/shaders/bboxVShader.vs",
		"./share/shaders/bboxFShader.fs",
		NULL
	);
	uniManager.add_uniform(boxShader, "view");
	uniManager.add_uniform(boxShader, "model");
	uniManager.add_uniform(boxShader, "projection");
	uniManager.add_uniform(boxShader, "boxColor");

	containerShader = Shader(
		"./share/shaders/containerVShader.vs",
		"./share/shaders/containerFShader.fs",
		NULL
	);
	uniManager.add_uniform(containerShader, "projection");
	uniManager.add_uniform(containerShader, "view");
	uniManager.add_uniform(containerShader, "model");
	uniManager.add_uniform(containerShader, "lightColor");
	uniManager.add_uniform(containerShader, "objectColor");
	uniManager.add_uniform(containerShader, "lightPosWorld");
	uniManager.add_uniform(containerShader, "Ka");

	frameShader = Shader(
		std::filesystem::relative(std::filesystem::path("./share/shaders/frameVShader.vs")).string().c_str(),
		std::filesystem::relative(std::filesystem::path("./share/shaders/frameFShader.fs")).string().c_str(),
		NULL
	);
	frameShader.use();
	uniManager.add_uniform(frameShader, "projection");
	uniManager.add_uniform(frameShader, "view");
	uniManager.add_uniform(frameShader, "model");
	uniManager.add_uniform(frameShader, "outColor");

	xArrow = std::make_unique<Arrow>();
	yArrow = std::make_unique<Arrow>();
	zArrow = std::make_unique<Arrow>();

	lineShader = Shader(
		"./share/shaders/lineShader.vs",
		"./share/shaders/lineShader.fs",
		NULL
	);
	lineShader.use();
	uniManager.add_uniform(lineShader, "projection");
	uniManager.add_uniform(lineShader, "view");
	uniManager.add_uniform(lineShader, "model");
	uniManager.add_uniform(lineShader, "lineColor");
	uniManager.setUniform(lineShader, "view", view);
	uniManager.setUniform(lineShader, "projection", projection);

	// axe lines
	lineX = std::make_unique<LineModel>(Vec3(-1.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f));
	lineY = std::make_unique<LineModel>(Vec3(0.0f, -1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f));
	lineZ = std::make_unique<LineModel>(Vec3(0.0f, 0.0f, -1.0f), Vec3(0.0f, 0.0f, 1.0f));

	//_draw_axes_lines();

	frameShader.use();
	uniManager.add_uniform(frameShader, "projection");
	uniManager.add_uniform(frameShader, "view");
	uniManager.add_uniform(frameShader, "model");
	uniManager.add_uniform(frameShader, "outColor");

	float sphereRadius{ 0.0f };
	if (width > height)
		sphereRadius = height * 0.5f;
	else
		sphereRadius = width * 0.5f;

	defCamera = new defaultCamera(window, 0.3f, glm::vec3(0.0f, 0.0f, 10.0f), cameraTarget, 2.0f);
	trackCamera = std::make_unique<TrackBall>(window, sphereRadius, framebuffer.width, framebuffer.height, 0, 0);
	
	scaffoldShader.use();
	Vec3 cameraPos = trackCamera->get_position();
	uniManager.setUniform(scaffoldShader, "lightPosWorld", cameraPos.x, cameraPos.y, cameraPos.z);

	// create frame buffer
	create_frame_buffer(framebuffer);

	// this is for the bounding box
	box = std::make_unique<BoundingBox>();

	// load textures
	int dummyW = 0, dummyH = 0;
	load_texture_from_file("./share/textures/boxContainer.png", &boxContainerTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/cylinderContainer.png", &cylinderContainerTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/abstractContainer.png", &abstractContainerTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/randomSeeds.png", &randomSeedTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/uniformSeeds.png", &uniformSeedTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/variedSeeds.png", &variedSeedTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/scaffoldCreator.png", &createScaffoldTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/thickness.png", &localThicknessTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/separation.png", &separationTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/tortuosity.png", &tortuosityTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/anisotropy.png", &anisotropyTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/network.png", &poreNetworkTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/trabecularNumber.png", &trabecularNumberTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/connectivityDensity.png", &connectivityDensityTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/update.png", &updateScaffoldTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/translate.png", &translateTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/scale.png", &scaleTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/tsmooth.png", &taubinSmoothTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/simplification.png", &simplifyMeshTexture, &dummyW, &dummyH);
};

void myGUI::_init_imgui() {

	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();

	io = ImGui::GetIO(); (void)io;
	io.Fonts->AddFontFromFileTTF("./share/fonts/DroidSans.ttf", fontSize);
	ImGui::GetIO().ConfigFlags |= ImGuiConfigFlags_DockingEnable;
	io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;         // Enable Docking
	io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;

	ImGui::StyleColorsDark();

	// When viewports are enabled we tweak WindowRounding/WindowBg so platform windows can look identical to regular ones.
	ImGuiStyle& style = ImGui::GetStyle();
	if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
	{
		style.WindowRounding = 0.0f;
		style.Colors[ImGuiCol_WindowBg].w = 1.0f;
	}

	ImGui::LoadIniSettingsFromDisk("./share/imgui.ini");

	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init(glsl_version.c_str());
}

void myGUI::run() {

	grid = std::make_unique<Grid>(Grid::XY);

	io = ImGui::GetIO();
	
	newFontSize = fontSize;

	while (!glfwWindowShouldClose(window))
	{

		glfwPollEvents();

		if (rebuildFont) {
			// Destroy the old OpenGL font texture to prevent memory leaks
			ImGui_ImplOpenGL3_DestroyFontsTexture();

			// Clear the CPU-side font data
			io.Fonts->Clear();

			// Load the newly sized font
			io.Fonts->AddFontFromFileTTF("./share/fonts/DroidSans.ttf", newFontSize);

			// Rebuild the OpenGL texture with the backend
			ImGui_ImplOpenGL3_CreateFontsTexture();

			// Reset the flag and sync the current size
			rebuildFont = false;
			fontSize = newFontSize;
		}

		// Start a new ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		_create_dockspace();

		if (showToolbar) {
			_render_toolbar();
		}

		// check here if the user pressed x to close the window
		if (showExitConfirmation) {
			ImGui::OpenPopup("Exit");
			showExitConfirmation = false;
		}

		// Draw the Modal
		if (ImGui::BeginPopupModal("Exit", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {

			ImGui::Text("Are you sure you want to quit?");

			if (ImGui::Button("Yes")) {
				// save imgui
				ImGui::SaveIniSettingsToDisk("./share/imgui.ini");
				_reset_scene();
				glfwSetWindowShouldClose(window, GLFW_TRUE);
				ImGui::CloseCurrentPopup();
			}
			ImGui::SameLine();
			if (ImGui::Button("No")) {
				ImGui::CloseCurrentPopup();
			}
			ImGui::EndPopup();
		}

		// render the console
		if (showConsole) {
			_render_console();
		}

		_render_settings_panel();

		_render_object_list();

		_render_properties_panel();

		_render_tool_panel();

		_render_console();

		if (showBinaryImageWindow) {
			_render_binary_image_window("Export Binary Image", showBinaryImageWindow);
		}

		glBindFramebuffer(GL_FRAMEBUFFER, framebuffer.id);
		glViewport(0, 0, framebuffer.width, framebuffer.height);

		glClearColor(fontColor[0], fontColor[1], fontColor[2], fontColor[3]);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if (showVisualizer) {
			_render_visualizer();
		}

		double xPos, yPos;
		glfwGetCursorPos(window, &xPos, &yPos);

		// ==============================================
		// First pass -> objects without alpha

		glEnable(GL_DEPTH_TEST);
		glDepthMask(GL_TRUE); 
		glDisable(GL_BLEND);
		glEnable(GL_CULL_FACE);

		if (cameraUpdate) {
			float sphereRadius;
			sphereRadius = std::min(static_cast<float>(framebuffer.width), static_cast<float>(framebuffer.height)) * 0.5f;
			trackCamera->set_radius(sphereRadius);
			//trackCamera->set_screen_size(0.3f * framebuffer.width, framebuffer.height, 0);
			trackCamera->update();
			Vec3 cameraPos = trackCamera->get_position();
			projection = trackCamera->get_projection_matrix();
			view = trackCamera->get_view_matrix();
			model = model;
		}

		// pass the standard values to the mesh shader
		scaffoldShader.use();
		// pass light
		uniManager.setUniform(scaffoldShader, "lightColor", lightColor[0], lightColor[1], lightColor[2]);
		uniManager.setUniform(scaffoldShader, "lightPosWorld", lightPosCamera[0], lightPosCamera[1], lightPosCamera[2]);
		//uniManager.setUniform(scaffoldShader, "lightPosWorld", cameraPos.x, cameraPos.y, cameraPos.z);
		uniManager.setUniform(scaffoldShader, "Ka", Ka);
		uniManager.setUniform(scaffoldShader, "Ks", Ks);

		// loop inside the list of scaffold generetors
		for (const auto& gen : scaffolds) {
			if (!gen->hidden && gen->color[3] == 1.0f) {
				
				uniManager.setUniform(scaffoldShader, "projection", projection);
				uniManager.setUniform(scaffoldShader, "view", view);
				uniManager.setUniform(scaffoldShader, "model", glm::translate(
					glm::mat4(1.0), 
					glm::vec3(gen->translateVec.x, gen->translateVec.y, gen->translateVec.z)));
				uniManager.setUniform(
					scaffoldShader, "objectColor",
					gen->color[0], gen->color[1], gen->color[2], gen->color[3]);
				glEnable(GL_POLYGON_OFFSET_FILL);
				glPolygonOffset(1.0f, 1.0f);

				if (showCutPlane && gen.get() == selectedSceneObj) {

					glDisable(GL_CULL_FACE);

					float d = cutPlane->normal.dot(cutPlane->center);

					planeCoeffs = glm::vec4(
						cutPlane->normal.x,
						cutPlane->normal.y,
						cutPlane->normal.z, d);

					uniManager.setUniform(scaffoldShader, "cutPlane", 1);
					uniManager.setUniform(
						scaffoldShader, "cutPlaneCoeffs", planeCoeffs);
					gen->draw();
					glEnable(GL_CULL_FACE);
				}
				else {
					uniManager.setUniform(scaffoldShader, "cutPlane", 0);
					gen->draw();
				}
				glDisable(GL_POLYGON_OFFSET_FILL);
			}
		}

		for (const auto& gen : scaffolds) {
			// if the object is the selected one draw specific
			
			if (showEdges && gen.get() == selectedSceneObj) {
				edgeShader.use();

				if (showCutPlane) {
					uniManager.setUniform(edgeShader, "cutPlane", 1);
					uniManager.setUniform(edgeShader, "cutPlaneCoeffs", planeCoeffs);
				}
				else {
					uniManager.setUniform(edgeShader, "cutPlane", 0);
				}
				glDepthFunc(GL_LEQUAL);
				glLineWidth(0.05);
				uniManager.setUniform(edgeShader, "projection", projection);
				uniManager.setUniform(edgeShader, "view", view);
				uniManager.setUniform(edgeShader, "model", glm::translate(
					glm::mat4(1.0),
					glm::vec3(gen->translateVec.x, gen->translateVec.y, gen->translateVec.z)));
				uniManager.setUniform(edgeShader, "edgeColor", 0.0f, 0.0f, 0.0f, 1.0f);
				if (showScaffold) {
					gen->draw_edges();
				}
			}

			if (showNormals && gen.get() == selectedSceneObj) {
				normalShader.use();
				uniManager.setUniform(normalShader, "projection", projection);
				uniManager.setUniform(normalShader, "view", view);
				uniManager.setUniform(normalShader, "model", glm::translate(
					glm::mat4(1.0),
					glm::vec3(gen->translateVec.x, gen->translateVec.y, gen->translateVec.z)));
				uniManager.setUniform(normalShader, "normalColor", normalColor[0], normalColor[1], normalColor[2]);
				if (showScaffold) {
					gen->draw();
				}
			}

			if (!gen->hiddenTortuosityPath && gen->tortuosityPathModel && gen.get() == selectedSceneObj) {
				glDepthFunc(GL_LEQUAL);
				glLineWidth(renderSettings.tortuosityPathSize);
				edgeShader.use();
				uniManager.setUniform(edgeShader, "projection", projection);
				uniManager.setUniform(edgeShader, "view", view);
				uniManager.setUniform(edgeShader, "model", glm::translate(
					glm::mat4(1.0),
					glm::vec3(gen->translateVec.x, gen->translateVec.y, gen->translateVec.z)));
				uniManager.setUniform(edgeShader, "cutPlane", 0);
				uniManager.setUniform(
					edgeShader,
					"edgeColor",
					renderSettings.tortuosityPathColor[0],
					renderSettings.tortuosityPathColor[1],
					renderSettings.tortuosityPathColor[2],
					renderSettings.tortuosityPathColor[3]
				);
				gen->draw_tortuosity_path();
			}

			if (!gen->hiddenEllipsoid && gen->ellipsoidModel && gen.get() == selectedSceneObj) {
				boxShader.use();
				uniManager.setUniform(boxShader, "projection", projection);
				uniManager.setUniform(boxShader, "view", view);
				uniManager.setUniform(boxShader, "model", model);
				uniManager.setUniform(boxShader, "boxColor", 0.1f, 1.0f, 0.1f, 1.0f);
				gen->ellipsoidModel->draw();

				glDepthFunc(GL_LEQUAL);
				glLineWidth(2.0f);
				lineShader.use();
				uniManager.setUniform(lineShader, "projection", projection);
				uniManager.setUniform(lineShader, "view", view);
				uniManager.setUniform(
					lineShader,
					"lineColor",
					0.0f, 0.0f, 0.0f, 1.0f
				);
				
				uniManager.setUniform(lineShader, "model", glm::mat4(1.0f));
				gen->ellipsoidModel->xAxis->draw();

				uniManager.setUniform(
					lineShader,
					"lineColor",
					0.0f, 0.0f, 0.0f, 1.0f
				);

				uniManager.setUniform(lineShader, "model", glm::mat4(1.0f));
				gen->ellipsoidModel->yAxis->draw();

				uniManager.setUniform(
					lineShader,
					"lineColor",
					0.0f, 0.0f, 0.0f, 1.0f
				);

				uniManager.setUniform(lineShader, "model", glm::mat4(1.0f));
				gen->ellipsoidModel->zAxis->draw();
			};

			if (showPoreNetwork && gen.get() == selectedSceneObj) {
				glDepthFunc(GL_LEQUAL);
				glLineWidth(renderSettings.poreNetworkLineSize);
				edgeShader.use();
				uniManager.setUniform(edgeShader, "projection", projection);
				uniManager.setUniform(edgeShader, "view", view);
				uniManager.setUniform(edgeShader, "model", glm::translate(
					glm::mat4(1.0),
					glm::vec3(gen->translateVec.x, gen->translateVec.y, gen->translateVec.z)));
				uniManager.setUniform(edgeShader, "cutPlane", 0);
				uniManager.setUniform(
					edgeShader,
					"edgeColor",
					renderSettings.poreNetworkColor[0],
					renderSettings.poreNetworkColor[1],
					renderSettings.poreNetworkColor[2],
					renderSettings.poreNetworkColor[3]
				);
				//gen->draw_pore_network();
			}
		}

		// containers
		for (const auto& con : containers) {
			if (con->get_type() == ObjectType::AbstractContainerType) {
				if (con && !con->hidden && con->color[3] == 1.0f) {
					scaffoldShader.use();
					uniManager.setUniform(scaffoldShader, "projection", projection);
					uniManager.setUniform(scaffoldShader, "view", view);
					uniManager.setUniform(scaffoldShader, "model", model);
					uniManager.setUniform(scaffoldShader, "cutPlane", 0);
					uniManager.setUniform(
						scaffoldShader,
						"objectColor",
						con->color[0], con->color[1], con->color[2], con->color[3]);
					con->render();
				}
			}
			else {
				if (con && !con->hidden) {
					boxShader.use();
					uniManager.setUniform(boxShader, "projection", projection);
					uniManager.setUniform(boxShader, "view", view);
					uniManager.setUniform(boxShader, "model", model);
					uniManager.setUniform(
						boxShader,
						"boxColor",
						con->color[0],
						con->color[1],
						con->color[2],
						con->color[3]
					);
					con->render();
				}
			}
		}

		// ROIs
		for (const auto& roi : rois) {
			if (roi && !roi->hidden) {
				boxShader.use();
				uniManager.setUniform(boxShader, "projection", projection);
				uniManager.setUniform(boxShader, "view", view);
				uniManager.setUniform(boxShader, "model", model);
				uniManager.setUniform(boxShader, "boxColor",
					roi->color[0],
					roi->color[1],
					roi->color[2],
					roi->color[3]);
				roi->render_model();
			}
		}

		if (!io.WantCaptureKeyboard) {
			if (ImGui::IsKeyPressed(ImGuiKey_M)) {
				showEdges = !showEdges;
			}

			if (ImGui::IsKeyPressed(ImGuiKey_G)) {
				showGrid = !showGrid;
			}

			if (ImGui::IsKeyPressed(ImGuiKey_N)) {
				showNormals = !showNormals;
			}
		}

		if (showGrid) {
			gridShader.use();
			uniManager.setUniform(gridShader, "projection", projection);
			uniManager.setUniform(gridShader, "view", view);
			uniManager.setUniform(gridShader, "model", model);
			uniManager.setUniform(gridShader, "gridColor", gridColor[0], gridColor[1], gridColor[2]);
			grid->draw();
		}

		for (const auto& seedGen : seedGenerators) {

			if (showSeeds && !seedGen->hidden && (seedGen.get() == selectedPanelObj.ptr)) {

				seedSize = seedGen->modelSeedSize;

				seedShader.use();
				uniManager.setUniform(seedShader, "projection", projection);
				uniManager.setUniform(seedShader, "view", view);
				uniManager.setUniform(seedShader, "seedSize", seedSize);
				uniManager.setUniform(seedShader, "seedColor", seedColor[0], seedColor[1], seedColor[2]);
				seedGen->draw();				
			}
		}

		if (showCutPlane && renderSettings.cutPlaneColor[3] == 1.0f) {
			cutShader.use();
			uniManager.setUniform(cutShader, "projection", projection);
			uniManager.setUniform(cutShader, "view", view);
			uniManager.setUniform(cutShader, "model", cutPlane->modelMatrix);
			uniManager.setUniform(cutShader, "minBounds", bounds[0], bounds[2], bounds[4]);
			uniManager.setUniform(cutShader, "maxBounds", bounds[1], bounds[3], bounds[5]);
			uniManager.setUniform(cutShader, "planeColor",
				renderSettings.cutPlaneColor[0],
				renderSettings.cutPlaneColor[1],
				renderSettings.cutPlaneColor[2],
				renderSettings.cutPlaneColor[3]
			);
			glDisable(GL_CULL_FACE);
			cutPlane->draw();
			glEnable(GL_CULL_FACE);
		}

		if (showAxesLines) {
			_draw_axes_lines();
		}

		if (box && showBbox) {
			_draw_selected_box();
		}

		// ---------------------------------------------------------------------
		// Pass 2, draw transparent objects

		glEnable(GL_BLEND);       // Enable blending
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glDepthMask(GL_FALSE);

		// draw all models
		for (auto& scaffold : scaffolds) {
			
			if (!scaffold->hidden && scaffold->color[3] < 1.0f) {
				scaffoldShader.use();
				uniManager.setUniform(scaffoldShader, "projection", projection);
				uniManager.setUniform(scaffoldShader, "view", view);
				uniManager.setUniform(scaffoldShader, "model", glm::translate(
					glm::mat4(1.0),
					glm::vec3(scaffold->translateVec.x, scaffold->translateVec.y, scaffold->translateVec.z)));
				uniManager.setUniform(
					scaffoldShader,
					"objectColor",
					scaffold->color[0], scaffold->color[1], scaffold->color[2], scaffold->color[3]);
				glEnable(GL_POLYGON_OFFSET_FILL);
				glPolygonOffset(1.0f, 1.0f);
				scaffold->draw();
				glDisable(GL_POLYGON_OFFSET_FILL);
			}
		}

		if (showCutPlane && renderSettings.cutPlaneColor[3] < 1.0f) {
			cutShader.use();
			uniManager.setUniform(cutShader, "projection", projection);
			uniManager.setUniform(cutShader, "view", view);
			uniManager.setUniform(cutShader, "model", cutPlane->modelMatrix);
			uniManager.setUniform(cutShader, "minBounds", bounds[0], bounds[2], bounds[4]);
			uniManager.setUniform(cutShader, "maxBounds", bounds[1], bounds[3], bounds[5]);
			uniManager.setUniform(cutShader, "planeColor",
				renderSettings.cutPlaneColor[0],
				renderSettings.cutPlaneColor[1],
				renderSettings.cutPlaneColor[2],
				renderSettings.cutPlaneColor[3]
			);
			glDisable(GL_CULL_FACE);
			cutPlane->draw();
			glEnable(GL_CULL_FACE);
		}

		for (const auto& con : containers) {
			if (con->get_type() == ObjectType::AbstractContainerType) {
				if (con && !con->hidden && con->color[3] < 1.0f) {
					scaffoldShader.use();
					uniManager.setUniform(scaffoldShader, "projection", projection);
					uniManager.setUniform(scaffoldShader, "view", view);
					uniManager.setUniform(scaffoldShader, "model", model);
					uniManager.setUniform(scaffoldShader, "cutPlane", 0);
					uniManager.setUniform(
						scaffoldShader,
						"objectColor",
						con->color[0], con->color[1], con->color[2], con->color[3]);
					con->render();
				}
			}
		}

		// Anisotropy Sources
		for (const auto& source : anisoSources) {
			if (source && !source->hidden && source->color[3] < 1.0f) {
				boxShader.use();
				uniManager.setUniform(boxShader, "projection", projection);
				uniManager.setUniform(boxShader, "view", view);
				// change the translation
				uniManager.setUniform(
					boxShader, "model", 
					glm::translate(
					glm::mat4(1.0),
					glm::vec3(
						source->origin.x, 
						source->origin.y, 
						source->origin.z 
					)));
				uniManager.setUniform(boxShader,
					"boxColor",
					source->color[0],
					source->color[1],
					source->color[2],
					source->color[3]
				);
				source->render_sphere_model();
				uniManager.setUniform(boxShader, "view", view);
				uniManager.setUniform(boxShader,
					"boxColor",
					source->lineColor[0],
					source->lineColor[1],
					source->lineColor[2],
					source->lineColor[3]
				);
				source->render_line_model();
			}
		}

		io = ImGui::GetIO();

		if (ImGui::IsKeyPressed(ImGuiKey::ImGuiKey_Delete) && !io.WantCaptureKeyboard) {
			
			logger.log(LogPriority::INFO, "Deleting everything!");

			scaffolds.clear();
			containers.clear();
			seedGenerators.clear();
			selectedPanelObj.ptr = nullptr;
			selectedPanelObj.type = ObjectType::NoneType;
			selectedSceneObj = nullptr;
			box.reset();
			_reset_camera();
		}

		// update for next frame
		glDepthMask(GL_TRUE);    
		glDisable(GL_BLEND);
		glDepthFunc(GL_LESS);

		_render_axes_viewport();

		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		// Render ImGui into the OpenGL context
 		ImGui::Render();

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		update_render();

		// Swap buffers
		glfwSwapBuffers(window);

		glfwPollEvents();
	}

}

void myGUI::clean() {
	// clean the gui
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
};

void myGUI::update_render() {
	if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
	{
		GLFWwindow* backup_current_context = glfwGetCurrentContext();
		ImGui::UpdatePlatformWindows();
		ImGui::RenderPlatformWindowsDefault();
		glfwMakeContextCurrent(backup_current_context);
	}
}

void myGUI::_render_toolbar() {
	ImGuiWindowFlags toolbarFlags =
		ImGuiWindowFlags_NoTitleBar |
		//ImGuiWindowFlags_NoResize |
		//ImGuiWindowFlags_NoMove |
		ImGuiWindowFlags_NoScrollbar |
		ImGuiWindowFlags_NoScrollWithMouse |
		//	ImGuiWindowFlags_NoSavedSettings |
		//	ImGuiWindowFlags_NoCollapse |
		ImGuiWindowFlags_AlwaysAutoResize;

	ImGuiTableFlags tableFlags = ImGuiTableFlags_BordersInnerV |
		ImGuiTableFlags_NoHostExtendX |
		ImGuiTableFlags_SizingFixedFit;

	if (ImGui::Begin("Toolbar", &showToolbar, toolbarFlags)) {

		ImGui::PushStyleVar(ImGuiStyleVar_CellPadding, ImVec2(20.0f, 4.0f));

		ImGui::BeginTable("", 5, tableFlags, ImVec2(0.0f, ImGui::GetFrameHeight()));

		ImGui::TableNextColumn();

		create_single_button_textured(
			"Create Box Container", showBoxContainerCreator, "Create a box container", boxContainerTexture);
		ImGui::SameLine();

		create_single_button_textured(
			"Create Cylinder Container", showCylinderContainerCreator, "Create a cylinder container", cylinderContainerTexture);

		ImGui::SameLine();

		create_single_button_textured(
			"Create Abstract Container", showAbstractContainerCreator, "Create an abstract container loading an .stl file.", abstractContainerTexture);

		ImGui::SameLine();
		
		ImGui::TableNextColumn();
		if(create_single_button_textured(
			"Create Random Seed Generator", showRandomSeedCreator, "create random seeds inside a container", randomSeedTexture)){
				seedFactory->set_type(ObjectType::RandomGeneratorType);
				seedFactory->launch();
			};
		ImGui::SameLine();

		if(create_single_button_textured(
			"Create Uniform Seed Generator", showUniformSeedCreator, "Create seeds inside a container using uniform Poisson sampling", uniformSeedTexture)){
				seedFactory->set_type(ObjectType::UniformGeneratorType);
				seedFactory->launch();
			};
		ImGui::SameLine();

		if(create_single_button_textured(
			"Create Varied Seed Generator", showVariedSeedCreator, "Create seeds inside a container using varied Poisson sampling", variedSeedTexture)){
				seedFactory->set_type(ObjectType::VariedGeneratorType);
				seedFactory->launch();
			};
		ImGui::SameLine();

		ImGui::TableNextColumn();


		create_single_button_textured(
			"Measure Local Thickness", measureThickness,
			"Measure the local thickness of the scaffold", localThicknessTexture);

		ImGui::SameLine();
		create_single_button_textured(
			"Measure Pore Separation", measureSeparation,
			"Measure pore separation of the scaffold", separationTexture);

		ImGui::SameLine();
		create_single_button_textured(
			"Measure Tortuoisity", estimateTortuosity,
			"Measure the tortuosity of the scaffold and visualize the path", tortuosityTexture);

		ImGui::SameLine();
		create_single_button_textured(
			"Measure Anisotropy", estimateAnisotropy,
			"Measure the anisotropy of the scaffold using the Mean Intercept Length method.", anisotropyTexture);

		ImGui::SameLine();
		
		// create_single_button_textured(
		// 	"Measure Porosity Network", estimatePoreNetwork,
		// 	"Create a graph that visualizes the connectivity network and estimates the percentage of interconneted seeds.", poreNetworkTexture
		// );
		create_single_button_textured(
			"Estimate Structural Model Index", estimateSmi,
			"Estimate Structure Model Index (SMI)", poreNetworkTexture
		);
		ImGui::SameLine();

		create_single_button_textured(
			"Estimate Connectivity Density", estimateConnectivityDensity,
			"Estimate connectivity density of the active scaffold.", connectivityDensityTexture
		);
		ImGui::SameLine();

		create_single_button_textured(
			"Estimate Trabecular Number", estimateTrabecularNr,
			"Estimate trabecular number of the active scaffold.", trabecularNumberTexture
		);
		ImGui::SameLine();

		ImGui::TableNextColumn();

		if (create_single_button_textured(
			"Create Scaffold", showScaffoldCreator,
			"Create the scaffold selecting a container and a generator", createScaffoldTexture
		)) {
			factory->launch();
		};
		
		ImGui::SameLine();

		create_single_button_textured(
			"Update Current Scaffold", updateScaffold,
			"Update Current Scaffold.", updateScaffoldTexture
		);

		if (updateScaffold && selectedPanelObj.type != ObjectType::ScaffoldType){
			updateScaffold = false;
		}

		ImGui::SameLine();
		ImGui::TableNextColumn();

		create_button_textured("Translate Object", translateScaffold, "Translate selected object.", translateTexture);
		ImGui::SameLine();

		create_button_textured("Scale Object", scaleScaffold, "Scale selected object.", scaleTexture);
		ImGui::SameLine();

		create_single_button_textured("Taubin Mesh", taubinSmooth, "Smooth With Taubin", taubinSmoothTexture);

		ImGui::SameLine();
		create_single_button_textured(
			"Simplify Mesh", simplify,
			 "Simplify mesh with QEM", simplifyMeshTexture);

		ImGui::EndTable();

		ImGui::PopStyleVar();
		//if (ImGui::Button("Estimate Connectivity")) {
		//	_action_estimate_connectivity();
		//}

	}

	ImGui::End();
};

void myGUI::_render_settings_panel() {

	// display settings window
	if (showDisplaySettingsWin) {
		_render_display_settings();
	}

	if (showDisplayMeshSettingsWin) {
		_render_mesh_settings();
	}

	if (showAlgorithmSettings) {
		_render_algorithm_settings();
	}

	ImVec2 maxSize = ImVec2((float)width, (float)height);
	ImVec2 minSize = ImVec2(500.0f, 500.0f);
	if (ImGuiFileDialog::Instance()->Display("Export Scaffold", ImGuiWindowFlags_NoCollapse, minSize, maxSize)) {
		if (ImGuiFileDialog::Instance()->IsOk()) { 
			scaffoldFilePath = ImGuiFileDialog::Instance()->GetFilePathName();
			scaffoldFileName = ImGuiFileDialog::Instance()->GetCurrentFileName();
			
			// get the active scaffold
			GeneratorLewiner* gen = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);

			std::filesystem::path filePath = scaffoldFilePath;

			if (gen && filePath.extension() == ".scaf") {

				gen->export_scaf(scaffoldFilePath);

				logger.log(LogPriority::SUCCESS, "Exported scaffold as " + scaffoldFilePath);

			}
		}
		ImGuiFileDialog::Instance()->Close();
	}

	if (ImGuiFileDialog::Instance()->Display("Save Project", ImGuiWindowFlags_NoCollapse, minSize, maxSize)) {
		if (ImGuiFileDialog::Instance()->IsOk()) {
			std::string filePath = ImGuiFileDialog::Instance()->GetFilePathName();
			std::filesystem::path p = filePath;
			if (p.extension() != ".sbproj") p.replace_extension(".sbproj");
			save_project(p.string());
		}
		ImGuiFileDialog::Instance()->Close();
	}

	if (ImGuiFileDialog::Instance()->Display("Open Project", ImGuiWindowFlags_NoCollapse, minSize, maxSize)) {
		if (ImGuiFileDialog::Instance()->IsOk()) {
			std::string filePath = ImGuiFileDialog::Instance()->GetFilePathName();
			load_project(filePath);
		}
		ImGuiFileDialog::Instance()->Close();
	}

	if (ImGuiFileDialog::Instance()->Display("Export Mesh Scaffold", ImGuiWindowFlags_NoCollapse, minSize, maxSize)) {
		if (ImGuiFileDialog::Instance()->IsOk()) { // action if OK
			scaffoldFilePath = ImGuiFileDialog::Instance()->GetFilePathName();
			scaffoldFileName = ImGuiFileDialog::Instance()->GetCurrentFileName();

			// get the active scaffold
			GeneratorLewiner* gen = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);

			std::filesystem::path filePath = scaffoldFilePath;

			if (filePath.extension() == ".stl") {

				gen->export_stl(scaffoldFilePath);

				logger.log(LogPriority::SUCCESS, "Exported scaffold mesh to " + scaffoldFilePath);
				showGeometryExportWindow = false;

			}
		}
		ImGuiFileDialog::Instance()->Close();
		showGeometryExportWindow = false;
	}

	if (ImGuiFileDialog::Instance()->Display("Export Binary Scaffold", ImGuiWindowFlags_NoCollapse, minSize, maxSize)) {
		if (ImGuiFileDialog::Instance()->IsOk()) { // action if OK
			scaffoldFilePath = ImGuiFileDialog::Instance()->GetFilePathName();
			scaffoldFileName = ImGuiFileDialog::Instance()->GetCurrentFileName();
			std::string folder = ImGuiFileDialog::Instance()->GetCurrentPath();

			GeneratorLewiner* gen = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);

			std::filesystem::path filePath = scaffoldFilePath;

			if (filePath.extension() == ".mhd") {

				gen->export_mhd(filePath, voxelSize, voxelBounds);

				logger.log(LogPriority::SUCCESS, "Exported scaffold as a binary image to " + scaffoldFilePath);

			}

			else if (filePath.extension() == ".nrrd") {

				gen->export_nrrd(scaffoldFilePath, voxelSize, voxelBounds);

				logger.log(LogPriority::SUCCESS, "Exported scaffold as a binary image to " + scaffoldFilePath +
					" with voxelSize: " + std::to_string(voxelSize) + ".");

			}
		}
		ImGuiFileDialog::Instance()->Close();
	}

	maxSize = ImVec2(width, height);
	minSize = ImVec2(500.0f, 500.0f);
	if (ImGuiFileDialog::Instance()->Display("Load Container", ImGuiWindowFlags_NoCollapse, minSize, maxSize)) {
		if (ImGuiFileDialog::Instance()->IsOk()) { 

			std::string filePath = ImGuiFileDialog::Instance()->GetFilePathName();
			std::string fileName = ImGuiFileDialog::Instance()->GetCurrentFileName();
			
			// create a container object
			std::unique_ptr<AbstractContainer> con = std::make_unique<AbstractContainer>(filePath);

			con->name = fileName;
			con->hidden = false;
			con->color[3] = 0.1f;
			_update_cameras(*con);

			// push it
			containers.push_back(std::move(con));

			selectedPanelObj.ptr = containers.back().get();
			selectedPanelObj.type = ObjectType::AbstractContainerType;

			showAbstractContainerCreator = false;
		}
		ImGuiFileDialog::Instance()->Close();
		showAbstractContainerCreator = false;
	}


	if (ImGuiFileDialog::Instance()->Display("Export Metrics", ImGuiWindowFlags_NoCollapse, minSize, maxSize)) {
		if (ImGuiFileDialog::Instance()->IsOk()) { // action if OK
			std::string filePath = ImGuiFileDialog::Instance()->GetFilePathName();
			std::string fileName = ImGuiFileDialog::Instance()->GetCurrentFileName();
			std::string folder = ImGuiFileDialog::Instance()->GetCurrentPath();

			GeneratorLewiner* gen = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);

			std::filesystem::path path = scaffoldFilePath;

			gen->export_metrics(filePath);

			logger.log(LogPriority::SUCCESS, "Metrics exported to: " + filePath + "!");
		}
		ImGuiFileDialog::Instance()->Close();
	}

		if (ImGuiFileDialog::Instance()->Display("Export Parameters", ImGuiWindowFlags_NoCollapse, minSize, maxSize)) {
		if (ImGuiFileDialog::Instance()->IsOk()) { // action if OK
			std::string filePath = ImGuiFileDialog::Instance()->GetFilePathName();
			std::string fileName = ImGuiFileDialog::Instance()->GetCurrentFileName();
			std::string folder = ImGuiFileDialog::Instance()->GetCurrentPath();

			GeneratorLewiner* gen = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);

			std::filesystem::path path = scaffoldFilePath;

			gen->export_parameters(filePath);

			logger.log(
				LogPriority::SUCCESS, "Parameters exported to: " + filePath + "!");
		}
		ImGuiFileDialog::Instance()->Close();
	}

	if (ImGuiFileDialog::Instance()->Display("Load Scaffold", ImGuiWindowFlags_NoCollapse, minSize, maxSize)) {
		if (ImGuiFileDialog::Instance()->IsOk()) {
			std::string filePath = ImGuiFileDialog::Instance()->GetFilePathName();
			std::filesystem::path p = filePath;

			loadTask = std::make_unique<ScaffoldLoadTask>();
			loadTask->scaffold = std::make_unique<GeneratorLewiner>();
			loadTask->scaffold->set_logger(&logger);
			loadTask->scaffold->name = p.stem().string();

			GeneratorLewiner* raw = loadTask->scaffold.get();
			ScaffoldLoadTask*  lts = loadTask.get();

			loadTask->task.start([raw, lts, filePath]() {
				raw->load_scaf(
					filePath,
					lts->newContainers,
					lts->newGenerators,
					lts->newAnisoSources,
					&lts->stage
				);
			});

			ImGui::OpenPopup("##ScaffoldLoadProgress");
		}
		ImGuiFileDialog::Instance()->Close();
	}

	// Progress modal — shown while a scaffold is loading in the background
	if (ImGui::BeginPopupModal("##ScaffoldLoadProgress", nullptr,
		ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar)) {

		static const char* const kStageMessages[] = {
			"Initializing...",
			"Reading parameters...",
			"Loading anisotropy sources...",
			"Loading container...",
			"Loading seed generator...",
			"Reading scalar field...",
			"Running marching cubes...",
			"Complete!"
		};
		static const float kStageProgress[] = {
			0.02f, 0.10f, 0.20f, 0.35f, 0.45f, 0.60f, 0.65f, 1.0f
		};

		if (loadTask) {
			int s = std::clamp(loadTask->stage.load(std::memory_order_relaxed), 0, 7);

			ImGui::Text("Loading Scaffold");
			ImGui::Separator();
			ImGui::Spacing();
			ImGui::TextUnformatted(kStageMessages[s]);
			ImGui::ProgressBar(kStageProgress[s], ImVec2(360.0f, 0.0f));
			ImGui::Spacing();

			if (loadTask->task.poll()) {
				_finalize_scaffold_load();
				ImGui::CloseCurrentPopup();
			}
		} else {
			ImGui::CloseCurrentPopup();
		}

		ImGui::EndPopup();
	}

	if (ImGuiFileDialog::Instance()->Display("Load Mesh Scaffold", ImGuiWindowFlags_NoCollapse, minSize, maxSize)) {
		if (ImGuiFileDialog::Instance()->IsOk()) { // action if OK
			std::string filePath = ImGuiFileDialog::Instance()->GetFilePathName();
			std::string fileName = ImGuiFileDialog::Instance()->GetCurrentFileName();
			
			std::unique_ptr<GeneratorLewiner> scaffold = std::make_unique<GeneratorLewiner>(filePath, &logger);
			scaffold->name = fileName;

			std::filesystem::path p = filePath;
			std::string paramFileName = p.stem().string() + "_parameters.csv";
			std::string parametersFile = (p.parent_path() / paramFileName).string();

			scaffold->read_parameters(parametersFile);

			std::string metricsFileName = p.stem().string() + "_metrics.csv";
			std::string metricsFile = (p.parent_path() / metricsFileName).string();
			scaffold->read_metrics(metricsFile);

			_update_cameras(*scaffold);

			// push to the scaffold list
			scaffolds.push_back(std::move(scaffold));

			// set it as the selected object
			selectedSceneObj = scaffolds.back().get();
			selectedPanelObj.ptr = scaffolds.back().get();
			selectedPanelObj.type = ObjectType::ScaffoldType;

		}
		ImGuiFileDialog::Instance()->Close();
	}

};

void myGUI::_finalize_scaffold_load() {
	if (!loadTask) return;

	// OpenGL calls must run on the main thread
	for (auto& src : loadTask->newAnisoSources)
		src->update_model();
	for (auto& gen : loadTask->newGenerators)
		gen->update_model();

	for (auto& con : loadTask->newContainers)
		containers.push_back(con);
	for (auto& gen : loadTask->newGenerators)
		seedGenerators.push_back(gen);
	for (auto& src : loadTask->newAnisoSources)
		anisoSources.push_back(src);

	// Upload mesh data to GPU and sync cameras
	loadTask->scaffold->update_render();
	_update_cameras(*loadTask->scaffold);

	scaffolds.push_back(std::move(loadTask->scaffold));

	selectedSceneObj = scaffolds.back().get();
	selectedPanelObj.ptr = scaffolds.back().get();
	selectedPanelObj.type = ObjectType::ScaffoldType;

	loadTask.reset();
}

// ===========================================================================
// Project (.sbproj) scene serialization
// ===========================================================================
// The .sbproj file is a manifest that stores the shared scene entities once
// (containers, generators, anisotropy sources, ROIs). Each scaffold is written
// as its own full-fidelity .scaf next to the project (guaranteeing a byte
// identical scalar field on reload) and referenced by filename + indices into
// the shared entity lists. Mesh containers embed their geometry so the project
// is self-contained and portable.

void myGUI::save_project(const std::string& path) {
	namespace fs = std::filesystem;

	std::ofstream out(path, std::ios::out | std::ios::binary);
	if (!out.is_open()) {
		logger.log(LogPriority::ERROR, "Failed to open project file for writing: " + path);
		return;
	}

	auto writeU32 = [&](uint32_t v) { out.write(reinterpret_cast<const char*>(&v), sizeof(v)); };
	auto writeI32 = [&](int32_t v) { out.write(reinterpret_cast<const char*>(&v), sizeof(v)); };
	auto writeU64 = [&](uint64_t v) { out.write(reinterpret_cast<const char*>(&v), sizeof(v)); };
	auto writeF   = [&](float v) { out.write(reinterpret_cast<const char*>(&v), sizeof(v)); };
	auto writeStr = [&](const std::string& s) { writeU32((uint32_t)s.size()); if (!s.empty()) out.write(s.data(), s.size()); };
	auto writeVec3 = [&](const Vec3& v) { writeF(v.x); writeF(v.y); writeF(v.z); };

	// Header
	const char magic[4] = { 'S', 'B', 'P', 'J' };
	out.write(magic, 4);
	writeU32(1); // format version

	// Block 1: containers
	writeU32((uint32_t)containers.size());
	for (auto& c : containers) {
		ObjectType t = c->get_type();
		if (t == ObjectType::BoxContainerType) {
			writeI32(1);
			auto* box = static_cast<BoxContainer*>(c.get());
			writeVec3(box->size);
			writeVec3(box->origin);
		}
		else if (t == ObjectType::CylinderContainerType) {
			writeI32(2);
			auto* cyl = static_cast<CylinderContainer*>(c.get());
			writeF(cyl->cylinderRadius);
			writeF(cyl->cylinderHeight);
		}
		else if (t == ObjectType::AbstractContainerType) {
			writeI32(3);
			auto* mesh = static_cast<AbstractContainer*>(c.get());
			writeF(mesh->get_scale());
			const auto& verts = mesh->get_mesh_verts();
			const auto& faces = mesh->get_mesh_faces();
			writeU64((uint64_t)verts.size());
			if (!verts.empty()) out.write(reinterpret_cast<const char*>(verts.data()), sizeof(openstl::Vec3) * verts.size());
			writeU64((uint64_t)faces.size());
			if (!faces.empty()) out.write(reinterpret_cast<const char*>(faces.data()), sizeof(openstl::Face) * faces.size());
		}
		else {
			writeI32(0); // unknown type placeholder (keeps indices aligned)
		}
		writeStr(c->name);
	}

	// Block 2: generators
	writeU32((uint32_t)seedGenerators.size());
	for (auto& g : seedGenerators) {
		if (g->get_type() == ObjectType::RandomGeneratorType) {
			writeI32(0);
			auto* r = static_cast<Random*>(g.get());
			writeI32((int32_t)r->seedNr);
		}
		else { // Poisson / Uniform / Varied all reconstruct as Poisson3D
			writeI32(1);
			auto* p = static_cast<Poisson3D*>(g.get());
			writeF(p->get_min_radius());
			writeF(p->get_max_radius());
		}
		writeStr(g->name);
	}

	// Block 3: anisotropy sources
	writeU32((uint32_t)anisoSources.size());
	for (auto& s : anisoSources) {
		writeVec3(s->origin);
		writeVec3(s->direction);
		writeVec3(s->stretch);
		writeF(s->sigma);
		writeF(s->angle);
		writeStr(s->name);
	}

	// Block 4: ROIs
	writeU32((uint32_t)rois.size());
	for (auto& r : rois) {
		writeVec3(r->get_size());
		writeVec3(r->get_center());
		for (int i = 0; i < 4; i++) writeF(r->color[i]);
		writeStr(r->name);
		writeStr(r->relatedMeshName);
	}

	// Block 5: scaffolds (each written as its own .scaf next to the project)
	fs::path folder = fs::path(path).parent_path();

	auto indexOf = [](const auto& list, const auto& item) -> int32_t {
		for (size_t i = 0; i < list.size(); ++i) if (list[i] == item) return (int32_t)i;
		return -1;
	};

	writeU32((uint32_t)scaffolds.size());
	std::vector<std::string> usedNames;
	for (auto& gen : scaffolds) {
		// derive a unique, filesystem-safe .scaf filename from the scaffold name
		std::string base = gen->name.empty() ? "scaffold" : gen->name;
		for (auto& ch : base)
			if (ch == '/' || ch == '\\' || ch == ':' || ch == '*' || ch == '?' ||
				ch == '"' || ch == '<' || ch == '>' || ch == '|') ch = '_';
		std::string fname = base + ".scaf";
		int dup = 1;
		while (std::find(usedNames.begin(), usedNames.end(), fname) != usedNames.end())
			fname = base + "_" + std::to_string(dup++) + ".scaf";
		usedNames.push_back(fname);

		writeStr(fname);
		writeI32(indexOf(containers, gen->container.lock()));
		writeI32(indexOf(seedGenerators, gen->generator.lock()));
		writeU32((uint32_t)gen->anisotropySources.size());
		for (auto& s : gen->anisotropySources) writeI32(indexOf(anisoSources, s));

		gen->set_logger(&logger);
		gen->export_scaf((folder / fname).string());
	}

	out.close();
	logger.log(LogPriority::SUCCESS, "Project saved to " + path);
}

bool myGUI::load_project(const std::string& path) {
	namespace fs = std::filesystem;

	std::ifstream in(path, std::ios::in | std::ios::binary);
	if (!in.is_open()) {
		logger.log(LogPriority::ERROR, "Failed to open project file: " + path);
		return false;
	}

	auto readU32 = [&]() { uint32_t v = 0; in.read(reinterpret_cast<char*>(&v), sizeof(v)); return v; };
	auto readI32 = [&]() { int32_t v = 0; in.read(reinterpret_cast<char*>(&v), sizeof(v)); return v; };
	auto readU64 = [&]() { uint64_t v = 0; in.read(reinterpret_cast<char*>(&v), sizeof(v)); return v; };
	auto readF   = [&]() { float v = 0.0f; in.read(reinterpret_cast<char*>(&v), sizeof(v)); return v; };
	auto readStr = [&]() { uint32_t n = readU32(); std::string s; if (n) { s.resize(n); in.read(&s[0], n); } return s; };
	auto readVec3 = [&]() { Vec3 v; v.x = readF(); v.y = readF(); v.z = readF(); return v; };

	char magic[4] = { 0 };
	in.read(magic, 4);
	if (magic[0] != 'S' || magic[1] != 'B' || magic[2] != 'P' || magic[3] != 'J') {
		logger.log(LogPriority::ERROR, "Invalid file format! Not a .sbproj project.");
		return false;
	}
	uint32_t version = readU32();
	if (version != 1) {
		logger.log(LogPriority::ERROR, "Unsupported .sbproj version.");
		return false;
	}

	// Block 1: containers  (null entries keep index alignment for scaffold refs)
	std::vector<std::shared_ptr<IContainer>> loadedContainers;
	uint32_t nc = readU32();
	for (uint32_t i = 0; i < nc; i++) {
		int32_t tid = readI32();
		std::shared_ptr<IContainer> con;
		if (tid == 1) {
			Vec3 size = readVec3(); Vec3 origin = readVec3();
			con = std::make_shared<BoxContainer>(size, origin);
		}
		else if (tid == 2) {
			float radius = readF(); float height = readF();
			con = std::make_shared<CylinderContainer>(radius, height);
		}
		else if (tid == 3) {
			float scale = readF();
			uint64_t vc = readU64();
			std::vector<openstl::Vec3> verts(vc);
			if (vc) in.read(reinterpret_cast<char*>(verts.data()), sizeof(openstl::Vec3) * vc);
			uint64_t fc = readU64();
			std::vector<openstl::Face> faces(fc);
			if (fc) in.read(reinterpret_cast<char*>(faces.data()), sizeof(openstl::Face) * fc);
			con = std::make_shared<AbstractContainer>(std::move(verts), std::move(faces), scale);
		}
		std::string name = readStr();
		if (con) {
			con->name = name;
			con->hidden = false;
			containers.push_back(con);
		}
		loadedContainers.push_back(con);
	}

	// Block 2: generators
	std::vector<std::shared_ptr<InterfaceSeedGenerator>> loadedGenerators;
	uint32_t ng = readU32();
	for (uint32_t i = 0; i < ng; i++) {
		int32_t tid = readI32();
		std::shared_ptr<InterfaceSeedGenerator> gen;
		if (tid == 0) {
			int32_t seedNr = readI32();
			gen = std::make_shared<Random>(seedNr);
		}
		else {
			float rmin = readF(); float rmax = readF();
			gen = std::make_shared<Poisson3D>((double)rmin, (double)rmax, 30);
		}
		std::string name = readStr();
		if (gen) {
			gen->name = name;
			gen->update_model();
			seedGenerators.push_back(gen);
		}
		loadedGenerators.push_back(gen);
	}

	// Block 3: anisotropy sources
	std::vector<std::shared_ptr<AnisotropySource>> loadedSources;
	uint32_t ns = readU32();
	for (uint32_t i = 0; i < ns; i++) {
		auto s = std::make_shared<AnisotropySource>();
		s->origin = readVec3();
		s->direction = readVec3();
		s->stretch = readVec3();
		s->sigma = readF();
		s->angle = readF();
		s->name = readStr();
		s->update_metric();
		s->update_model();
		loadedSources.push_back(s);
		anisoSources.push_back(s);
	}

	// Block 4: ROIs
	uint32_t nr = readU32();
	for (uint32_t i = 0; i < nr; i++) {
		Vec3 size = readVec3(); Vec3 origin = readVec3();
		auto roi = std::make_shared<ROI>(size, origin, true);
		for (int k = 0; k < 4; k++) roi->color[k] = readF();
		roi->name = readStr();
		roi->relatedMeshName = readStr();
		rois.push_back(roi);
	}

	// Block 5: scaffolds
	fs::path folder = fs::path(path).parent_path();
	uint32_t nsc = readU32();
	for (uint32_t i = 0; i < nsc; i++) {
		std::string fname = readStr();
		int32_t ci = readI32();
		int32_t gi = readI32();
		uint32_t anc = readU32();
		std::vector<int32_t> anis(anc);
		for (uint32_t a = 0; a < anc; a++) anis[a] = readI32();

		auto scaffold = std::make_unique<GeneratorLewiner>();
		scaffold->set_logger(&logger);
		scaffold->name = fs::path(fname).stem().string();

		// load_scaf rebuilds its own container/generator/sources into throwaway
		// lists; we discard those and relink to the shared project entities.
		std::vector<std::shared_ptr<IContainer>> tmpCon;
		std::vector<std::shared_ptr<InterfaceSeedGenerator>> tmpGen;
		std::vector<std::shared_ptr<AnisotropySource>> tmpSrc;
		if (!scaffold->load_scaf((folder / fname).string(), tmpCon, tmpGen, tmpSrc)) {
			logger.log(LogPriority::ERROR, "Failed to load scaffold file: " + fname);
			continue;
		}

		if (ci >= 0 && ci < (int)loadedContainers.size() && loadedContainers[ci])
			scaffold->container = loadedContainers[ci];
		else
			scaffold->container.reset();

		if (gi >= 0 && gi < (int)loadedGenerators.size() && loadedGenerators[gi])
			scaffold->generator = loadedGenerators[gi];
		else
			scaffold->generator.reset();

		scaffold->anisotropySources.clear();
		for (int32_t a : anis)
			if (a >= 0 && a < (int)loadedSources.size() && loadedSources[a])
				scaffold->anisotropySources.push_back(loadedSources[a]);

		scaffold->update_render();
		scaffolds.push_back(std::move(scaffold));
	}

	in.close();

	if (!scaffolds.empty()) {
		_update_cameras(*scaffolds.back());
		selectedSceneObj = scaffolds.back().get();
		selectedPanelObj.ptr = scaffolds.back().get();
		selectedPanelObj.type = ObjectType::ScaffoldType;
	}

	logger.log(LogPriority::SUCCESS, "Project loaded from " + path);
	return true;
}

void myGUI::_render_object_list() {

	if (ImGui::Begin("Objects", &showScaffoldList)) {

		bool open = ImGui::TreeNodeEx("Scaffolds", ImGuiDockNodeFlags_None | ImGuiTreeNodeFlags_DefaultOpen);
		if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort | ImGuiHoveredFlags_NoSharedDelay)) {
			ImGui::SetTooltip("Scaffold objects");
		}
		if (open) {

			static RenameState state;
			static bool doDelete = false;

			for (int i = 0; i < scaffolds.size(); ++i) {

				ImGui::PushID(i);

				auto& gen = scaffolds[i];

				bool isSelected = (selectedPanelObj.ptr == gen.get());

				ImGui::Selectable(gen->name.c_str(), isSelected, ImGuiSelectableFlags_AllowDoubleClick);
			
				if (ImGui::BeginPopupContextItem()) {
					
					if (ImGui::MenuItem("Select")) {
						selectedSceneObj = gen.get();
					};
					if (ImGui::MenuItem(!gen->hidden ? "Hide" : "Show")) {
						gen->hidden = !gen->hidden;
					}
					if (ImGui::MenuItem(!gen->hiddenTortuosityPath ? "Hide Tortuosity Model" : "Show Tortuosity Model")) {
						gen->hiddenTortuosityPath = !gen->hiddenTortuosityPath;
					}
					if (ImGui::MenuItem(!gen->hiddenEllipsoid ? "Hide Ellipsoid Model" : "Show Ellipsoid Model")) {
						gen->hiddenEllipsoid = !gen->hiddenEllipsoid;
					}
					if (ImGui::MenuItem("Delete")) {
						doDelete = true;
					}
					if (ImGui::MenuItem("Export as Mesh")) {
						showGeometryExportWindow = true;
						IGFD::FileDialogConfig config;
						config.path = "../data";
						ImGuiFileDialog::Instance()->OpenDialog(
							"Export Mesh Scaffold", "Export Scaffold Geometry", ".stl, .vtk", config);
					}
					if (ImGui::MenuItem("Export as Image")) {
						showBinaryImageWindow = true;
					}
					if (ImGui::MenuItem("Export Metrics")) {
						showMetricsExportWindow = true;
						IGFD::FileDialogConfig config;
						config.path = "../data";
						ImGuiFileDialog::Instance()->OpenDialog("Export Metrics", "Export Scaffold Metrics", ".csv", config);
					}
					if (ImGui::MenuItem("Export Parameters")) {
						showMetricsExportWindow = true;
						IGFD::FileDialogConfig config;
						config.path = "../data";
						ImGuiFileDialog::Instance()->OpenDialog(
							"Export Parameters", 
							"Export Scaffold Parameters", ".csv", config);
					}

					if (ImGui::MenuItem("Edit Name")) {
						// index to change the name
						state.targetIdx = i;
						strncpy_s(state.buffer, gen->name.c_str(), sizeof(state.buffer));
						state.buffer[sizeof(state.buffer) - 1] = 0;
						state.showPopup = true;
					}

					ImGui::EndPopup();
				}

				if (doDelete) {
					auto it = scaffolds.erase(scaffolds.begin() + i);
					selectedSceneObj = nullptr;
					selectedPanelObj.ptr = nullptr;
					selectedPanelObj.type = ObjectType::NoneType;
					doDelete = false;
				}

				if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && ImGui::IsItemHovered()) {
					selectedPanelObj.ptr = gen.get();
					selectedPanelObj.type = ObjectType::ScaffoldType;
				}

				// double click selects the active scene object
				if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {

					if (ImGui::IsItemFocused() || ImGui::IsItemHovered()) {
						selectedSceneObj = gen.get();
					}
				}

				ImGui::PopID();

			}

			render_change_name_popup(scaffolds, state);

			ImGui::TreePop();

		}
		
		_render_container_list();

		_render_seed_generator_list();

		_render_roi_list();

		_render_anisource_list();
	}

	ImGui::End();
};

void myGUI::_render_tool_panel() {

	if (showPlaneCutSettings) {
		_render_cutting_plane_settings(
			"Cutting With Plane" , showPlaneCutSettings);
	}

	if (translateScaffold) {
		_render_translate_panel("Translate Object", translateScaffold);
	}

	if (scaleScaffold) {
		_render_scale_panel("Scale Object", scaleScaffold);
	}

	if (taubinSmooth) {
		_render_taubin_smooth_panel("Taubin Smoothing", taubinSmooth);
	}

	if (simplify){
		_render_simplify_mesh_panel("Simplify Mesh", simplify);
	}

};

void myGUI::_render_properties_panel() {

	if (ImGui::Begin("Properties", &showProperties)) {

		switch (selectedPanelObj.type) {
			case ObjectType::BoxContainerType: {
				BoxContainer* box = static_cast<BoxContainer*>(selectedPanelObj.ptr);
				box->gui_setup();
				break;
			}
			case ObjectType::CylinderContainerType: {
				CylinderContainer* cyl = static_cast<CylinderContainer*>(selectedPanelObj.ptr);
				cyl->gui_setup();
				break;
			}
			case ObjectType::AbstractContainerType: {
				AbstractContainer* con = static_cast<AbstractContainer*>(selectedPanelObj.ptr);
				con->gui_setup();
				if (con->updated) {
					_update_cameras(*con);
					con->updated = false;
				}
				break;
			}
			case ObjectType::RandomGeneratorType:
			case ObjectType::UniformGeneratorType:
			case ObjectType::VariedGeneratorType:{
				_render_seed_generator_properties();
				break;
			}
			// case ObjectType::UniformGeneratorType: {
			// 	_render_uniform_seed_generator_properties();
			// 	break;
			// }
			// case ObjectType::VariedGeneratorType: {
			// 	_render_varied_seed_generator_properties();
			// 	break;
			// }
			case ObjectType::ScaffoldType: {
				GeneratorLewiner* scaffold = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);
				scaffold->render_properties(
					updateScaffold, gTask.get(), anisoSources);
				break;
			}
			case ObjectType::Roi: {
				ROI* roi = static_cast<ROI*>(selectedPanelObj.ptr);
				roi->render_properties();
				break;
			}
			case ObjectType::AnisoSource: {
				AnisotropySource* roi = static_cast<AnisotropySource*>(selectedPanelObj.ptr);
				roi->render_properties();
				break;
			}
			case ObjectType::NoneType: {
				break;
			}
		}

	}
	ImGui::End();
};

void myGUI::_render_seed_generator_list() {
	bool open = ImGui::TreeNodeEx("Seed Generators", ImGuiDockNodeFlags_None | ImGuiTreeNodeFlags_DefaultOpen);
	if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort | ImGuiHoveredFlags_NoSharedDelay)) {
		ImGui::SetTooltip("Created Seed Generators");
	}
	if (open) {

		static RenameState state;

		for (int i = 0; i < seedGenerators.size(); ++i) {

			auto& gen = seedGenerators[i];

			bool isSelected = (selectedPanelObj.ptr == gen.get());

			ImGui::Selectable(gen->name.c_str(), isSelected, ImGuiSelectableFlags_AllowDoubleClick);

			if (ImGui::BeginPopupContextItem()) {

				if (ImGui::MenuItem(!gen->hidden ? "Hide" : "Show")) {
					gen->hidden = !gen->hidden;
				}
				if (ImGui::MenuItem("Delete")) {
					auto it = seedGenerators.erase(seedGenerators.begin() + i);
					ImGui::EndPopup();
					selectedPanelObj.ptr = nullptr;
					selectedPanelObj.type = ObjectType::NoneType;
					break;
				}

				if (ImGui::MenuItem("Edit Name")) {
					// index to change the name
					state.targetIdx = i;
					strncpy_s(state.buffer, gen->name.c_str(), sizeof(state.buffer));
					state.buffer[sizeof(state.buffer) - 1] = 0;
					state.showPopup = true;
				}
				ImGui::EndPopup();
			}
			if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && ImGui::IsItemHovered()) {
				selectedPanelObj.ptr = gen.get();
				selectedPanelObj.type = gen->type;
			}
		}
		render_change_name_popup(seedGenerators, state);

		ImGui::TreePop();
	}
};

void myGUI::_render_roi_list() {
	bool open = ImGui::TreeNodeEx("ROIs", ImGuiDockNodeFlags_None | ImGuiTreeNodeFlags_DefaultOpen);
	if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort | ImGuiHoveredFlags_NoSharedDelay)) {
		ImGui::SetTooltip("Created ROIs");
	}
	if (open) {

		static RenameState state;

		for (int i = 0; i < rois.size(); ++i) {

			auto& roi = rois[i];

			bool isSelected = (selectedPanelObj.ptr == roi.get());

			ImGui::Selectable(roi->name.c_str(), isSelected, ImGuiSelectableFlags_AllowDoubleClick);

			if (ImGui::BeginPopupContextItem()) {

				if (ImGui::MenuItem(!roi->hidden ? "Hide" : "Show")) {
					roi->hidden = !roi->hidden;
				}
				if (ImGui::MenuItem("Delete")) {
					auto it = rois.erase(rois.begin() + i);
					ImGui::EndPopup();
					selectedPanelObj.ptr = nullptr;
					selectedPanelObj.type = ObjectType::NoneType;
					break;
				}

				ImGui::EndPopup();
			}
			if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && ImGui::IsItemHovered()) {
				selectedPanelObj.ptr = roi.get();
				selectedPanelObj.type = roi->type;
			}
		}

		ImGui::TreePop();
	}
};

void myGUI::_render_anisource_list() {
	bool open = ImGui::TreeNodeEx("Anisotropy Sources", ImGuiDockNodeFlags_None | ImGuiTreeNodeFlags_DefaultOpen);
	if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort | ImGuiHoveredFlags_NoSharedDelay)) {
		ImGui::SetTooltip("Created Anisotropy Sources");
	}
	if (open) {

		static RenameState state;

		for (int i = 0; i < anisoSources.size(); ++i) {

			ImGui::PushID(i);

			auto& source = anisoSources[i];

			bool isSelected = (selectedPanelObj.ptr == source.get());

			ImGui::Selectable(source->name.c_str(), isSelected, ImGuiSelectableFlags_AllowDoubleClick);

			if (ImGui::BeginPopupContextItem()) {

				if (ImGui::MenuItem(!source->hidden ? "Hide" : "Show")) {
					source->hidden = !source->hidden;
				}
				if (ImGui::MenuItem("Delete")) {
					auto it = anisoSources.erase(
						anisoSources.begin() + i);

					ImGui::EndPopup();
					selectedPanelObj.ptr = nullptr;
					selectedPanelObj.type = ObjectType::NoneType;
					break;
				}

				ImGui::EndPopup();
			}
			if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && ImGui::IsItemHovered()) {
				selectedPanelObj.ptr = source.get();
				selectedPanelObj.type = source->get_type();
			}

			ImGui::PopID();
		}

		ImGui::TreePop();
	}
};

void myGUI::_render_container_list() {

	bool open = ImGui::TreeNodeEx("Containers", ImGuiDockNodeFlags_None | ImGuiTreeNodeFlags_DefaultOpen);
	if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort | ImGuiHoveredFlags_NoSharedDelay)) {
		ImGui::SetTooltip("Created container objects");
	}
	if (open) {

		static RenameState state;

		for (int i = 0; i < containers.size(); ++i) {

			auto& con = containers[i];

			ImGui::PushID(i);

			bool isSelected = (selectedPanelObj.ptr == con.get() && selectedPanelObj.type == con.get()->get_type());

			ImGui::Selectable(con->name.c_str(), isSelected, ImGuiSelectableFlags_AllowDoubleClick);

			if (ImGui::BeginPopupContextItem()) {

				if (ImGui::MenuItem(!con->hidden ? "Hide" : "Show")) {
					con->hidden = !con->hidden;
				}
				if (ImGui::MenuItem("Delete")) {
					auto it = containers.erase(containers.begin() + i);
					ImGui::EndPopup();
					selectedPanelObj.ptr = nullptr;
					selectedPanelObj.type = ObjectType::NoneType;
					break;
				}

				if (ImGui::MenuItem("Edit Name")) {
					// index to change the name
					state.targetIdx = i;
					strncpy_s(state.buffer, con->name.c_str(), sizeof(state.buffer));
					state.buffer[sizeof(state.buffer) - 1] = 0;
					state.showPopup = true;
				}

				ImGui::EndPopup();
			}

			if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && ImGui::IsItemHovered()) {
				selectedPanelObj.ptr = con.get();
				selectedPanelObj.type = con->get_type();
			}

			ImGui::PopID();
		}
		render_change_name_popup(containers, state);

		ImGui::TreePop();
	}

};

void myGUI::write_settings() {

	nlohmann::json settings;

	//settings["Scaffold"]["version"] = version;
	//settings["Scaffold"]["Domain"]["xMin"] = 0.0;
	//settings["Scaffold"]["Domain"]["xMax"] = xDim;
	//settings["Scaffold"]["Domain"]["yMin"] = 0.0;
	//settings["Scaffold"]["Domain"]["yMax"] = yDim;
	//settings["Scaffold"]["Domain"]["zMin"] = 0.0;
	//settings["Scaffold"]["Domain"]["zMax"] = zDim;
	//settings["Scaffold"]["Pores"]["genOption"] = generateOption;
	//settings["Scaffold"]["Pores"]["poreNr"] = seedNr;
	//settings["Scaffold"]["thickness"] = thickness;
	//settings["Scaffold"]["regSteps"] = regSteps;
	//settings["Voro"]["nX"] = 6;
	//settings["Voro"]["nY"] = 6;
	//settings["Voro"]["nZ"] = 6;

	//std::ofstream file("settings.json");
	//file << settings.dump(4);
	//file.close();

};

std::string myGUI::_get_glsl_version() {
	// Get OpenGL version
	const char* glVersion = reinterpret_cast<const char*>(glGetString(GL_VERSION));

	if (!glVersion) {
		std::cerr << "Error: Could not retrieve OpenGL version." << std::endl;
		return "#version 130";  // Default fallback
	}

	//std::cout << "OpenGL version: " << glVersion << std::endl;
	//logger.log(LogPriority::INFO, "OpenGL version: " + std::string(glVersion));

	// Extract the major and minor OpenGL version numbers
	int major, minor;
	sscanf_s(glVersion, "%d.%d", &major, &minor);

	// Map OpenGL version to corresponding GLSL version
	if (major == 3 && minor >= 3) {
		return "#version 330";
	}
	else if (major == 4) {
		if (minor >= 5) return "#version 450";
		if (minor >= 4) return "#version 440";
		if (minor >= 3) return "#version 430";
		if (minor >= 2) return "#version 420";
		if (minor >= 1) return "#version 410";
		return "#version 400";
	}
	else if (major >= 4) {
		return "#version 460";
	}

	// Default to a safe GLSL version if nothing else matches
	return "#version 130";  // Safe fallback for OpenGL 3.0-3.2
}

// console

void myGUI::_render_console() {

	if (ImGui::Begin("Console", NULL)) {
		// Start a scrolling region inside the console
		ImGui::BeginChild("ScrollingConsole", ImVec2(0, 120), true, ImGuiWindowFlags_AlwaysVerticalScrollbar | ImGuiWindowFlags_AlwaysHorizontalScrollbar);

		for (int i{ 0 }; i < logger.get_logs().size(); i ++){

			ImVec4 color = {
				logger.get_colors().at(i)[0],
				logger.get_colors().at(i)[1],
				logger.get_colors().at(i)[2],
				logger.get_colors().at(i)[3],
			};

			ImGui::TextColored(color, logger.get_logs().at(i).c_str());
		}
		
		if (logger.check_and_clear_new_messages() || scrollToBottom) {
			ImGui::SetScrollHereY(1.0f);
			scrollToBottom = false;
		}
		
		ImGui::EndChild();

		if (ImGui::Button("Clear")) {
			logger.clear();
			logger.log(LogPriority::INFO, "Cleared.");
		}
	}

	ImGui::End();

}

//void myGUI::logger.log(LogPriority priority, const std::string& message) {
//
//	if (priority == LogPriority::SUCCESS) {
//		logColor = { 0.0f, 1.0f, 0.0f, 1.0f };
//	}
//
//	else if (priority == LogPriority::ERROR) {
//		logColor = { 1.0f, 0.0f, 0.0f, 1.0f };
//	}
//
//	else {
//		logColor = { 1.0f, 1.0f, 1.0f, 1.0f };
//	}
//
//	logger.log(priority, message, logColor);
//
//	scrollToBottom = true;
//}

void myGUI::_update_cameras(const GeneratorLewiner& gen) {

	auto bnds = gen.get_bounds();

	//float xc = static_cast<float>((xMax + xMin) * 0.5f);
	//float yc = static_cast<float>((yMax + yMin) * 0.5f);
	//float zc = static_cast<float>((zMax + zMin) * 0.5f);

	float dx = bnds[1] - bnds[0];
	float dy = bnds[3] - bnds[2];
	float dz = bnds[5] - bnds[4];
	float diagonal = std::sqrt(dx * dx + dy * dy + dz * dz);
	float distance = 1.1f * diagonal;

	float cx = (bnds[0] + bnds[1]) * 0.5f;
	float cy = (bnds[2] + bnds[3]) * 0.5f;
	float cz = (bnds[4] + bnds[5]) * 0.5f;

	cameraUpdate = true;
	cameraTarget = glm::vec3(cx, cy, cz);

	std::cout << cx << " " << cy << " " << cz << std::endl;

	defCamera->position = cameraPos;
	defCamera->target = cameraTarget;
	//trackCamera->position = cameraPos;
	trackCamera->set_target(cameraTarget.x, cameraTarget.y, cameraTarget.z);
	trackCamera->set_position(cx, cy, cz + distance);
	trackCamera->update();
	//lightPosCamera[0] = cx;
	//lightPosCamera[1] = cy;
	//lightPosCamera[2] = cz + distance;
	projection = trackCamera->get_projection_matrix();
	view = trackCamera->get_view_matrix();
	cameraUpdate = false;
}

void myGUI::_update_cameras(IContainer& con) {

	Bounds bnds = con.compute_bounds();

	//float xc = static_cast<float>((xMax + xMin) * 0.5f);
	//float yc = static_cast<float>((yMax + yMin) * 0.5f);
	//float zc = static_cast<float>((zMax + zMin) * 0.5f);

	float dx = bnds.xMax - bnds.xMin;
	float dy = bnds.yMax - bnds.yMin;
	float dz = bnds.zMax - bnds.zMin;
	float diagonal = std::sqrt(dx * dx + dy * dy + dz * dz);
	float distance = 1.1f * diagonal;

	cameraUpdate = true;
	cameraTarget = glm::vec3(bnds.center.x, bnds.center.y, bnds.center.z);
	defCamera->position = cameraPos;
	defCamera->target = cameraTarget;
	//trackCamera->position = cameraPos;
	trackCamera->set_target(cameraTarget.x, cameraTarget.y, cameraTarget.z);
	trackCamera->set_position(bnds.center.x, bnds.center.y, bnds.center.z + distance);
	trackCamera->update();
	projection = trackCamera->get_projection_matrix();
	view = trackCamera->get_view_matrix();
	cameraUpdate = false;
}

void myGUI::_reset_camera() {
	cameraUpdate = true;
	trackCamera->reset();
	trackCamera->update();
	projection = trackCamera->get_projection_matrix();
	view = trackCamera->get_view_matrix();
	cameraUpdate = false;
};

void myGUI::_reset_scene(){
	seedGenerators.clear();
	containers.clear();
	scaffolds.clear();
	anisoSources.clear();
	rois.clear();
	_reset_camera();
	selectedSceneObj = nullptr;
	selectedPanelObj.ptr = nullptr;
	selectedPanelObj.type = ObjectType::NoneType;
}

void myGUI::_render_mesh_settings() {

	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
	ImGui::SetNextWindowSize(ImVec2(400, 300), ImGuiCond_Appearing);

	if (ImGui::Begin("Mesh Settings", &showDisplayMeshSettingsWin, ImGuiWindowFlags_AlwaysAutoResize)) {

		static int picked = { -99 };
		// Left side - Selectable items
		{
			ImGui::BeginChild("LeftPanel", ImVec2(ImGui::GetContentRegionAvail().x * 0.5f, 300), ImGuiChildFlags_Borders);
			if (ImGui::Selectable("Resolution Settings", picked == 0)) picked = 0;
			ImGui::EndChild();
		}
		ImGui::SameLine();

		// Right side - Color and size controls
		{
			ImGui::BeginGroup();
			{
				ImGui::BeginChild("settings item", ImVec2(0, 300), ImGuiChildFlags_Borders);

				if (picked == 0) {
					ImGui::Text("Resolution Settings");
					ImGui::InputInt("x Res", &resolution[0]);
					ImGui::InputInt("y Res", &resolution[1]);
					ImGui::InputInt("z Res", &resolution[2]);
				}

				ImGui::EndChild();

				ImGui::SameLine();

				ImGui::EndGroup();
			}
		}
		ImGui::End();

	}
};

void myGUI::_render_algorithm_settings() {
	
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
	ImGui::SetNextWindowSize(ImVec2(600, 300), ImGuiCond_Appearing);

	if (ImGui::Begin("Algorithm Settings", &showAlgorithmSettings, ImGuiWindowFlags_AlwaysAutoResize)) {

		static int picked = { -99 };
		// Left side - Selectable items
		{
			ImGui::BeginChild("LeftPanel", ImVec2(ImGui::GetContentRegionAvail().x * 0.3f, 300), ImGuiChildFlags_Borders);
			if (ImGui::Selectable("Anisotropy Settings", picked == 0)) picked = 0;
			if (ImGui::Selectable("Tortuosity Settings", picked == 1)) picked = 1;
			if (ImGui::Selectable("Trabecular Number Settings", picked == 2)) picked = 2;
			ImGui::EndChild();
		}
		ImGui::SameLine();

		// Right side - Color and size controls
		{
			ImGui::BeginGroup();
			{
				ImGui::BeginChild("settings item", ImVec2(0, 300), ImGuiChildFlags_Borders);

				if (picked == 0) {
					ImGui::Text("Anisotropy Algorithm Settings");

					ImGui::InputInt("daMinVectors", &daMinLines);
					ImGui::SetItemTooltip("Set the number of parallel lines along each dimension.");
					ImGui::InputInt("Directions", &daMinDirections);
					ImGui::SetItemTooltip("Set the number of tested directions.");
					// Ratio form (>=1, isotropic = 1): matches ratio-convention
					// literature (e.g. Ulrich 1999). 1 - Lmax/Lmin dropped: always <= 0.
					ImGui::RadioButton("MaxRadius/MinRadius", &daFormulaIdx, 0);
					// [0,1] forms: eigenvalue ratio (isotropic = 1) and its complement
					// (isotropic = 0, the BoneJ convention).
					ImGui::RadioButton("MinEigValue/MaxEigValue", &daFormulaIdx, 2);
					ImGui::RadioButton("1 - MinEigValue/MaxEigValue", &daFormulaIdx, 3);
				}


				if (picked == 1) {
					ImGui::Text("Tortuosity Algorithm Settings");
					ImGui::InputFloat("Voxel Size (mm)", &tortuosityVoxelSize);
				}

				if (picked == 1) {
					ImGui::Text("Trabecular Number Algorithm Settings");

					// Index 0 matches formula == 0 (MIL/DDA), Index 1 matches formula == 1 (Derived Proxy)
					const char* formulaOptions[] = { "Mean Intercept Length (DDA)", "Derived Proxy (BV/TV / Tb.Th)" };

					ImGui::Combo("Algorithm", &trabecularNrFormula, formulaOptions, IM_ARRAYSIZE(formulaOptions));
				}

				ImGui::EndChild();

				ImGui::SameLine();

				ImGui::EndGroup();
			}
		}
		ImGui::End();
	}
};

void myGUI::_render_axes_viewport() {

	int windowWidth{ 0 }, windowHeight{ 0 };
	glfwGetWindowSize(window, &windowWidth, &windowHeight);

	// Save depth state
	glEnable(GL_SCISSOR_TEST);

	// Define the overlay area (e.g., 200x200 in the top right)
	int overlayW = 100;
	int overlayH = 100;
	int overlayX = 0;
	int overlayY = 0; 

	// create a camera
	sCamera = std::make_unique<SimpleCamera>(overlayW, overlayH);

	glm::mat4 rotMatrix = trackCamera->get_rotation_matrix();
	glm::mat4 frameProjection = sCamera->get_projection_matrix();

	glm::mat4 frameView = glm::translate(glm::mat4(1.0), glm::vec3(0.0, 0.0, -5.0)) * rotMatrix;

	glScissor(overlayX, overlayY, overlayW, overlayH);

	// Render the inset scene
	glViewport(overlayX, overlayY, overlayW, overlayH);

	frameShader.use();
	glm::mat4 frameModel = glm::scale(glm::mat4(1.0), glm::vec3(0.5f));
	frameShader.use();
	uniManager.setUniform(frameShader, "projection", frameProjection);
	uniManager.setUniform(frameShader, "view", frameView);
	uniManager.setUniform(frameShader, "model", glm::mat4(1.0f));
	uniManager.setUniform(frameShader, "height", height);
	uniManager.setUniform(frameShader, "outColor", glm::vec4(0.0f, 0.0f, 1.0f, 1.0f));
	zArrow->draw();
	uniManager.setUniform(frameShader, "model", glm::rotate(glm::mat4(1.0f), glm::radians(-90.0f), glm::vec3(1.0f, 0.0f, 0.0f)));
	uniManager.setUniform(frameShader, "outColor", glm::vec4(0.0f, 1.0f, 0.0f, 1.0f));
	yArrow->draw();
	uniManager.setUniform(frameShader, "model", glm::rotate(glm::mat4(1.0f), glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f)));
	uniManager.setUniform(frameShader, "outColor", glm::vec4(1.0f, 0.0f, 0.0f, 1.0f));
	xArrow->draw();

	glDisable(GL_SCISSOR_TEST);

};

void myGUI::_render_main_menu_bar() {

	// Add Menu Bar
	if (ImGui::BeginMainMenuBar()) {

		//ImGui::PushItemWidth(150);

		if (ImGui::BeginMenu("File")) {

			if (ImGui::IsItemHovered()) {
				ImGui::SetMouseCursor(ImGuiMouseCursor_Arrow);
			}

			if (ImGui::MenuItem("Save Settings", "CTRL+S")) {
				// Call your function to save model settings here
				write_settings();
				logger.log(LogPriority::INFO, "Settings saved.");
			}

			ImGui::Separator();

			if (ImGui::MenuItem("Save Project", "Save the whole scene as .sbproj")) {
				IGFD::FileDialogConfig config;
				config.path = "../data";
				config.flags = ImGuiFileDialogFlags_ConfirmOverwrite;
				ImGuiFileDialog::Instance()->OpenDialog("Save Project", "Save Project as .sbproj", ".sbproj", config);
			}

			if (ImGui::MenuItem("Open Project", "Load a scene from .sbproj")) {
				IGFD::FileDialogConfig config;
				config.path = "../data";
				ImGuiFileDialog::Instance()->OpenDialog("Open Project", "Open Project", ".sbproj", config);
			}

			ImGui::Separator();

			if (ImGui::MenuItem("Load Scaffold", "Load Scaffold from .scaf file")) {
				// Call your function to load a model mesh here
				IGFD::FileDialogConfig config;
				config.path = "..//data";
				ImGuiFileDialog::Instance()->OpenDialog("Load Scaffold", "Load Scaffold", ".scaf", config);
			}

			if (ImGui::MenuItem("Load Mesh Scaffold", "Load Scaffold from .stl file")) {
				// Call your function to load a model mesh here
				IGFD::FileDialogConfig config;
				config.path = "..//data";
				ImGuiFileDialog::Instance()->OpenDialog("Load Mesh Scaffold", "Load Scaffold Geometry", ".vtk, .stl", config);
			}

			if (ImGui::MenuItem("Export Scaffold", "Save Scaffold as .scaf")) {
				
				if(selectedPanelObj.type == ObjectType::ScaffoldType){
					IGFD::FileDialogConfig config;
					config.path = "../data";
					ImGuiFileDialog::Instance()->OpenDialog("Export Scaffold", "Save Scaffold as .scaf", ".scaf", config);
				}
			}

			if (ImGui::MenuItem("Save Scaffold as geometry", "Export scaffold geometry.")) {
				showGeometryExportWindow = true;
				IGFD::FileDialogConfig config;
				config.path = "../data";
				ImGuiFileDialog::Instance()->OpenDialog("Export Mesh Scaffold", "Export Scaffold Geometry", ".stl, .vtk", config);
			}

			if (ImGui::MenuItem("Save Scaffold as binary image", "Export scaffold as binary image.")) {
				showBinaryImageWindow = true;
			}

			ImGui::EndMenu();
		}

		// menu for view settings like normals, edges etc
		if (ImGui::BeginMenu("View")) {

			if (ImGui::IsItemHovered()) {
				ImGui::SetMouseCursor(ImGuiMouseCursor_Arrow);
			}

			if (ImGui::MenuItem("Show Scaffold", NULL, showScaffold)) {
				showScaffold = !showScaffold;
			}

			if (ImGui::MenuItem("Show Seeds", NULL, showSeeds)) {
				showSeeds = !showSeeds;
			}

			if (ImGui::MenuItem("Show Scaffold Face Normals", NULL, showNormals)) {
				showNormals = !showNormals;
			}

			if (ImGui::MenuItem("Show Scaffold Edges", NULL, showEdges)) {
				showEdges = !showEdges;
			}

			if (ImGui::MenuItem("Show Grid", NULL, showGrid)) {
				showGrid = !showGrid;
			}

			if (ImGui::MenuItem("Show Container", NULL, showContainer)) {
				showContainer = !showContainer;
			}

			if (ImGui::MenuItem("Show Pore Network", NULL, showPoreNetwork)) {
				showPoreNetwork = !showPoreNetwork;
			}

			if (ImGui::MenuItem("Show ToolBar", NULL, showToolbar)) {
				showToolbar = !showToolbar;
			}

			if (ImGui::MenuItem("Show Bounding Box", NULL, showBbox)) {
				showBbox = !showBbox;
			}

			if (ImGui::MenuItem("Show Axes Lines", NULL, showAxesLines)) {
				showAxesLines = !showAxesLines;
			}

			if (ImGui::MenuItem("Reset Scene")){
				_reset_scene();
			}

			ImGui::EndMenu();
		}

		if (ImGui::BeginMenu("Window")) {
			if (ImGui::MenuItem("Visualizer Window", NULL, showVisualizer)) {
				showVisualizer = !showVisualizer;
			}
			ImGui::EndMenu();
		}

		if (ImGui::BeginMenu("Settings")) {
			if (ImGui::IsItemHovered()) {
				ImGui::SetMouseCursor(ImGuiMouseCursor_Arrow);
			}

			if (ImGui::MenuItem("Display")) {
				showDisplaySettingsWin = true;
			}

			if (ImGui::MenuItem("Mesh")) {
				showDisplayMeshSettingsWin = true;
			}

			if (ImGui::MenuItem("Utilities")) {
				showAlgorithmSettings = true;
			}

			if (ImGui::BeginMenu("Units")) {
				if (ImGui::MenuItem("mm")) {

				}
				ImGui::EndMenu();
			}

			if (ImGui::BeginMenu("Select Camera")) {
				if (ImGui::MenuItem("Default", NULL, &defCameraFlag))
				{
					cameraOption = defaultOption;
					trackCameraFlag = false;
					defCamera = new defaultCamera(window, 0.3f, glm::vec3(0.0f, 0.0f, 10.0f), cameraTarget, 2.0f);
				}
				if (ImGui::MenuItem("TrackBall", NULL, &trackCameraFlag))
				{
					cameraOption = trackOption;
					defCameraFlag = false;

					float sphereRadius{ 0.0f };
					if (width > height) {
						sphereRadius = height * 0.5f;
					}
					else {
						sphereRadius = width * 0.5f;
					}

					trackCamera = std::make_unique<TrackBall>(
						window, sphereRadius, framebuffer.width, framebuffer.height, framebuffer.posx, framebuffer.posy);
				}

				ImGui::EndMenu();
			}

			ImGui::EndMenu();
		}

		if (ImGui::BeginMenu("Add")) {
			if (ImGui::MenuItem("Box ROI")) {
				if (selectedPanelObj.type != ObjectType::ScaffoldType) {
					logger.log(LogPriority::ERROR, "Selected a scaffold from the panel!");
				}
				else {
					showROICreator = true;
				}
			}

			if (ImGui::MenuItem("Add Anisotropy Source")) {
				showAnisotropySourceCreator = true;
			}

			ImGui::EndMenu();
			
		}
		if (ImGui::BeginMenu("Tools")) {

			if (ImGui::IsItemHovered()) {
				ImGui::SetMouseCursor(ImGuiMouseCursor_Arrow);
			}
			if (ImGui::MenuItem("Cut With Plane")) {

				if (selectedPanelObj.type != ObjectType::ScaffoldType) {
					logger.log(LogPriority::ERROR, "Selected a scaffold from the panel!");
				}
				else {
					showPlaneCutSettings = true;
					showCutPlane = true;
					cutPlane = std::make_unique<CutPlane>();
					GeneratorLewiner* gen = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);
					bounds = gen->get_bounds();
				}
			}			

			if (ImGui::MenuItem("Cut With ROI")) { 
				if (selectedPanelObj.type != ObjectType::ScaffoldType) {
					logger.log(LogPriority::ERROR, "Selected a scaffold from the panel!");
				}
				else {
					showROICutter = true;
				}
			}

			ImGui::EndMenu();
		}


		//selectFileButton("Export Scaffold", "../data/", "Export File", ".stl, .vtk");

		// Close the menu bar
		ImGui::EndMainMenuBar();
	}

	if (showBoxContainerCreator) {
		_render_box_container_creator("Box Container Creator", showBoxContainerCreator);
	}

	if (showCylinderContainerCreator) {
		_render_cylinder_container_creator("Cylindrical Container Creator", showCylinderContainerCreator);
	}

	if (showAbstractContainerCreator) {
		IGFD::FileDialogConfig config;
		config.path = "../data";
		ImGuiFileDialog::Instance()->OpenDialog("Load Container", "Load Container Geometry", ".stl", config);
	}

	// if (showRandomSeedCreator) {
	// 	_render_random_seed_creator("Random Seed Creator", showRandomSeedCreator);
	// }

	// if (showUniformSeedCreator) {
	// 	_render_uniform_seed_creator("Uniform Seed Creator", showUniformSeedCreator);
	// }

	// if (showVariedSeedCreator) {
	// 	_render_varied_seed_creator("Varied Seed Creator", showVariedSeedCreator);
	// }

	if (showRandomSeedCreator){
		IContainer* con = nullptr;
		con = seedFactory->gui_draw(
			&logger,
			"Random Seed Generator", showRandomSeedCreator,
			&selectedPanelObj, seedGenerators, containers
		);
		if(con) _update_cameras(*con);
	}
	
	if(showUniformSeedCreator){
		IContainer* con = nullptr;
		con = seedFactory->gui_draw(
			&logger,
			"Uniform Seed Generator", showUniformSeedCreator,
			&selectedPanelObj, seedGenerators, containers
		);
		if(con) _update_cameras(*con);
	}
	
	if(showVariedSeedCreator){
		IContainer* con = nullptr;
		con = seedFactory->gui_draw(
			&logger,
			"Varied Seed Generator", showVariedSeedCreator,
			&selectedPanelObj, seedGenerators, containers
		);
		if(con) _update_cameras(*con);
	} 

	if (showScaffoldCreator) {
		//_render_scaffold_creator("Scaffold Creator", showScaffoldCreator);
		factory->gui_draw(
			gTask.get(),
			&logger,
			"Scaffold Creator", showScaffoldCreator,
			&selectedPanelObj, selectedSceneObj, 
			scaffolds, containers, seedGenerators, anisoSources
		);
	}

	if (measureThickness) {
		_action_estimate_local_thickness("Local Thickness Measure", measureThickness, false);
	}

	if (measureSeparation) {
		_action_estimate_local_thickness("Local Separation Measure", measureSeparation, true);
	}

	if (estimateTortuosity) {
		_action_estimate_tortuosity();
	}

	if (estimateAnisotropy) {
		_action_estimate_anisotropy();
	}

	if (estimateConnectivityDensity) {
		_action_estimate_connectivity_density();
	}
	
	if (estimateTrabecularNr) {
		_action_estimate_trabecular_number();
	}

	// if (estimatePoreNetwork) {
	// 	_action_estimate_pore_network();
	// }
	if (estimateSmi){
		_action_estimate_smi();
	}

	if (showROICreator) {
		_render_roi_creator("Create a ROI", showROICreator);
	}

	if (showAnisotropySourceCreator) {
		_render_aniso_source_creator("Create Anisotropy Source", showAnisotropySourceCreator);
	}

	if (showROICutter) {
		_render_roi_cutter("Get a scaffold region inside a ROI", showROICutter);
	}
};

void myGUI::_render_display_settings() {

	static int pickedItem{ -99 };
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
	ImGui::SetNextWindowSize(ImVec2(700, 400), ImGuiCond_Appearing);
	if (ImGui::Begin("Display Settings", &showDisplaySettingsWin, ImGuiChildFlags_AlwaysAutoResize)) {

		// Left side - Selectable items
		{
			ImGui::BeginChild("LeftPanel", 
				ImVec2(ImGui::GetContentRegionAvail().x * 0.3f, 300), ImGuiChildFlags_Borders);
			if (ImGui::Selectable("Mesh", pickedItem == 0)) pickedItem = 0;
			if (ImGui::Selectable("Seeds", pickedItem == 1)) pickedItem = 1;
			if (ImGui::Selectable("Grid", pickedItem == 2)) pickedItem = 2;
			if (ImGui::Selectable("Lighting", pickedItem == 3)) pickedItem = 3;
			if (ImGui::Selectable("Container", pickedItem == 5)) pickedItem = 5;
			if (ImGui::Selectable("Pore Network", pickedItem == 6)) pickedItem = 6;
			if (ImGui::Selectable("Tortuosity Rendering", pickedItem == 7)) pickedItem = 7;
			if (ImGui::Selectable("Cutting Tool Settings", pickedItem == 8)) pickedItem = 8;
			if (ImGui::Selectable("Gui Settings", pickedItem == 9)) pickedItem = 9;
			ImGui::EndChild();
		}
		ImGui::SameLine();

		// Right side - Color and size controls
		{
			ImGui::BeginGroup();
			{
				ImGui::BeginChild("settings item", ImVec2(0, 300));

				if (pickedItem == 0) {
					ImGui::Text("Mesh Settings");
					ImGui::ColorEdit4("Mesh Color", (float*)&scaffoldColor);
				}
				else if (pickedItem == 1) {
					ImGui::Text("Seed Settings");
					ImGui::ColorEdit3("Seed Color", (float*)&seedColor);

					// The point size the renderer actually uses is the selected
					// generator's per-object modelSeedSize (the render loop copies
					// it into seedSize every frame). Bind the control to that so the
					// edit takes effect instead of being overwritten.
					InterfaceSeedGenerator* selGen = nullptr;
					for (auto& g : seedGenerators) {
						if (g.get() == selectedPanelObj.ptr) { selGen = g.get(); break; }
					}
					if (selGen) {
						ImGui::InputFloat("Point Size", &selGen->modelSeedSize, 0.001f, 1.0f, "%.3f");
					}
					else {
						ImGui::BeginDisabled();
						ImGui::InputFloat("Point Size", &seedSize, 0.001f, 1.0f, "%.3f");
						ImGui::EndDisabled();
						ImGui::TextDisabled("Select a seed generator to adjust its point size.");
					}
				}
				else if (pickedItem == 2) {
					ImGui::Text("Grid Settings");
					ImGui::ColorEdit3("Grid Color", (float*)&gridColor);
					ImGui::Checkbox("Use Grid", &showGrid);
				}
				else if (pickedItem == 3) {
					ImGui::Text("Lighting Settings");
					ImGui::ColorEdit3("Light Color", (float*)&lightColor);
					ImGui::ColorEdit3("Font Color", (float*)&fontColor);
					ImGui::ColorEdit3("Normal Color", (float*)&normalColor);
					ImGui::InputFloat("Ambient Strength", &Ka, 0.01f, 100.0f, "%.3f");
					ImGui::InputFloat("Specular Strength", &Ks, 0.01f, 100.0f, "%.3f");
				}
				if (pickedItem == 5) {
					ImGui::Text("Container Settings");
					ImGui::ColorEdit4("Container Color", (float*)&containerColor);
				}
				if (pickedItem == 6) {
					ImGui::Text("Pore Network Display Settings");
					ImGui::ColorEdit4("Line Color", renderSettings.poreNetworkColor.data());
					ImGui::InputFloat("Line Size", &renderSettings.poreNetworkLineSize, 0.1f, 1.0f, "%.3f");
				}

				if (pickedItem == 7) {
					ImGui::Text("Tortuosity Display Settings");
					ImGui::ColorEdit4("Line Color", renderSettings.tortuosityPathColor.data());
					ImGui::InputFloat("Line Size", &renderSettings.tortuosityPathSize, 0.1f, 1.0f, "%.3f");
				}
				if (pickedItem == 8) {
					ImGui::Text("Cutting Tool Settings");
					ImGui::ColorEdit4("Line Color", renderSettings.cutPlaneColor.data());
				}

				if (pickedItem == 9) {
					ImGui::Text("Gui Panel Settings");
					if (ImGui::InputFloat("Font Size", &newFontSize, 0.1f, 1.0f, "%.3f")) {
						if (newFontSize > 4.0f && newFontSize != fontSize) { // Prevent crashes from 0 or negative sizes
							rebuildFont = true;
						}
					};
				}

				ImGui::EndChild();

				ImGui::SameLine();

				ImGui::EndGroup();
			}
		}
		ImGui::End();

	}

};

void myGUI::_render_cutting_plane_settings(const char* popupName, bool& showPopup) {

	// always centered
	// ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	// ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

	if (ImGui::Begin(popupName, &showPopup, ImGuiWindowFlags_AlwaysAutoResize)) {
		ImGui::ColorEdit4("Color", renderSettings.cutPlaneColor.data());
		ImGui::Text("Origin");
		ImGui::SetNextItemWidth(100);
		ImGui::InputFloat("##x", &cutPlane->center.x, 0.01f, 1000.0f);
		ImGui::SameLine();
		ImGui::SetNextItemWidth(100);
		ImGui::InputFloat("##y", &cutPlane->center.y, 0.01f, 1000.0f); 
		ImGui::SameLine();
		ImGui::SetNextItemWidth(100);
		ImGui::InputFloat("##z", &cutPlane->center.z, 0.01f, 1000.0f);

		ImGui::Text("Normal");
		ImGui::SetNextItemWidth(100);
		ImGui::InputFloat("##nx", &cutPlane->normal.x, 0.01f, 1000.0f); 
		ImGui::SameLine();
		ImGui::SetNextItemWidth(100);
		ImGui::InputFloat("##ny", &cutPlane->normal.y, 0.01f, 1000.0f);
		ImGui::SameLine();
		ImGui::SetNextItemWidth(100);
		ImGui::InputFloat("##nz", &cutPlane->normal.z, 0.01f, 1000.0f);

		cutPlane->normal.normalize();

		cutPlane->update_model_matrix();

		if (ImGui::Button("Close")) {
			showPopup = false;
			showCutPlane = false;
			cutPlane.reset();
			ImGui::End();
		}
	}

	if (!showPopup) showCutPlane = false;

	ImGui::End();

};

void myGUI::_render_box_container_creator(const char* popupName, bool& showPopup) {

	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}

	// always centered
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

	if (ImGui::BeginPopupModal(popupName, NULL, ImGuiWindowFlags_AlwaysAutoResize))
	{
		static char buffer[256] = "";
		ImGui::InputText("Name", buffer, sizeof(buffer));

		// Display here the settings for a box using Size and Center
		// Default size of 10 and center at 5 matches your old [0, 10] bounds
		static float tempSize[3] = { 10.0f, 10.0f, 10.0f };
		static float tempOrigin[3] = { 5.0f, 5.0f, 5.0f };

		ImGui::SeparatorText("Dimensions");

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Width (X) mm", &tempSize[0], 0.1f, 100.0f, "%.3f");

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Height (Y) mm", &tempSize[1], 0.1f, 100.0f, "%.3f");

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Depth (Z) mm", &tempSize[2], 0.1f, 100.0f, "%.3f");

		ImGui::SeparatorText("Position");

		ImGui::InputFloat3("Center", tempOrigin);
		ImGui::SetItemTooltip("Center of the Box Container");

		ImGui::Separator();
		ImGui::Spacing();

		if (ImGui::Button("Create", ImVec2(120, 0))) {

			Vec3 sizeVec{ tempSize[0], tempSize[1], tempSize[2] };
			Vec3 originVec{ tempOrigin[0], tempOrigin[1], tempOrigin[2] };

			std::unique_ptr<BoxContainer> container = std::make_unique<BoxContainer>(sizeVec, originVec);

			std::string name = std::string(buffer);

			if (!name.empty()) {
				container->name = name;
			}
			else {
				container->name = "Container" + std::to_string(containers.size());
			}

			selectedPanelObj.type = ObjectType::BoxContainerType;

			_update_cameras(*container);

			// push to the container list
			containers.push_back(std::move(container));
			buffer[0] = '\0';
			selectedPanelObj.ptr = containers.back().get();

			logger.log(LogPriority::SUCCESS, "Created box container!");

			showPopup = false;
		}

		ImGui::SameLine();

		if (ImGui::Button("Cancel", ImVec2(120, 0))) {
			buffer[0] = '\0';
			showPopup = false;
		};

		ImGui::EndPopup();
	}
};

void myGUI::_render_cylinder_container_creator(const char* popupName, bool& showPopup) {

	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}
	// always centered
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

	if (ImGui::BeginPopupModal(popupName, NULL, ImGuiWindowFlags_AlwaysAutoResize))
	{
		static float tempRadius = { 5.0f };
		static float tempHeight = { 10.0f };
		
		static char buffer[256] = "";
		ImGui::InputText("Name", buffer, sizeof(buffer));

		ImGui::InputFloat("Cylinder Radius", &tempRadius, 0.01f, 1.0f);
		ImGui::InputFloat("Cylinder Height", &tempHeight, 0.01f, 1.0f);

		if (ImGui::Button("Create")) {

			std::unique_ptr<CylinderContainer> cyl = std::make_unique<CylinderContainer>(tempRadius, tempHeight);

			std::string name = buffer;

			if (name.empty()) {
				cyl->name = "Container" + std::to_string(containers.size());
			}
			else {
				cyl->name = name;
			}

			_update_cameras(*cyl);

			containers.push_back(std::move(cyl));

			selectedPanelObj.ptr = containers.back().get();
			selectedPanelObj.type = ObjectType::CylinderContainerType;
			buffer[0] = '\0';

			logger.log(LogPriority::SUCCESS, "Created cylindrical container!");

			showPopup = false;
		}

		ImGui::SameLine();
		if (ImGui::Button("Cancel")) {
			buffer[0] = '\0';
			showPopup = false;
		};

		ImGui::EndPopup();
	}
};

void myGUI::_render_abstract_container_creator(const char* popupName, bool& showPopup) {

	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}

	showPopup = false;
};

void myGUI::_render_container_picker(IContainer*& container){

	ImGui::BeginChild("Containers", ImVec2(ImGui::GetContentRegionAvail().x, 80), ImGuiChildFlags_Borders);

	for (const auto& md : containers) {

		IContainer* con = md.get();

		if(!con) continue;

		if(container == nullptr) container = con;

		bool isSelected = container == con;

		if (ImGui::Selectable(md->name.c_str(), isSelected)) {
			container = con;
		};

		if (isSelected) ImGui::SetItemDefaultFocus();
	}

	ImGui::EndChild();
};


void myGUI::_render_seed_generator_properties(){

	// get the object
	auto* gen = static_cast<InterfaceSeedGenerator*>(selectedPanelObj.ptr);

	// render the container selection panel
	
	IContainer* con = gen->container;
	_render_container_picker(con);

	// render the properties
	gen->render_gui();

	// update
	if (ImGui::Button("Update")){
		if(con){
			gen->container = con;
			gen->run(*con);
			logger.log(LogPriority::SUCCESS, "Seeds Updated!");
		}
	}
};


void myGUI::_create_dockspace() {

	ImGui::DockSpaceOverViewport(0, ImGui::GetMainViewport());

	_render_main_menu_bar();

};

void myGUI::_render_visualizer() {

	ImGui::SetNextWindowSize(ImVec2((float)framebuffer.width, (float)framebuffer.height), ImGuiCond_Once);
	ImGui::Begin("Visualizer", &showVisualizer);

	cameraUpdate = (ImGui::IsWindowHovered() && ImGui::IsWindowFocused());

	ImVec2 availableSize = ImGui::GetContentRegionAvail();
	ImVec2 pos = ImGui::GetCursorScreenPos();

	if ((int)availableSize.x != framebuffer.width || (int)availableSize.y != framebuffer.height) {
		framebuffer.width = (int)availableSize.x;
		framebuffer.height = (int)availableSize.y;
		framebuffer.posx = pos.x;
		framebuffer.posy = pos.y;
		create_frame_buffer(framebuffer);
	}

	// Update Camera Viewport
	trackCamera->set_viewport((int)pos.x, (int)pos.y, (int)availableSize.x, (int)availableSize.y);

	// Render the main scene texture
	ImGui::Image((ImTextureID)(intptr_t)framebuffer.textureId, availableSize, ImVec2(0, 1), ImVec2(1, 0));

	ImGui::End();

};

void myGUI::_render_binary_image_window(const char* popupName, bool& showPopup) {

	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
	ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_Always);

	if (ImGui::BeginPopup(popupName, ImGuiWindowFlags_AlwaysAutoResize)) {

		static float tempVoxelSize{ 0.05f };
		static std::array<float, 2> xDims = { 0.0f, 5.0f };
		static std::array<float, 2> yDims = { 0.0f, 5.0f };
		static std::array<float, 2> zDims = { 0.0f, 5.0f };

		ImGui::InputFloat("Voxel Size (mm)", &tempVoxelSize, 0.01f, 10.0f);

		ImGui::InputFloat2("X Dimensions (mm)", xDims);
		ImGui::InputFloat2("Y Dimensions (mm)", yDims);
		ImGui::InputFloat2("Z Dimensions (mm)", zDims);

		if (ImGui::Button("Save")) {

			voxelBounds = {xDims[0], xDims[1], yDims[0], yDims[1], zDims[0], zDims[1]};
			voxelSize = tempVoxelSize;

			IGFD::FileDialogConfig config;
			config.path = "../data";
			ImGuiFileDialog::Instance()->OpenDialog("Export Binary Scaffold", "Export binary scaffold", ".mhd, .nrrd", config);
			showPopup = false;
			ImGui::CloseCurrentPopup();

		}
		ImGui::SameLine();
		if (ImGui::Button("Close")) {
			showPopup = false;
			ImGui::CloseCurrentPopup();
		}
		ImGui::EndPopup();
	}
};

void myGUI::_on_close_request() {

	glfwSetWindowShouldClose(window, GLFW_FALSE);
	showExitConfirmation = true;
};

void myGUI::_action_estimate_local_thickness(const char* popupName, bool& showPopup, bool flag) {
	
	// grab the active generatior
	if (selectedPanelObj.type != ObjectType::ScaffoldType) {

		logger.log(LogPriority::ERROR, "Select a scaffold in the panel.");
		showPopup = false;
		return;
	}
	GeneratorLewiner* gen = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);

	if (!gen) {
		showPopup = false;
		return;
	}

	if (gen->isLoadedFromFile) {
		logger.log(LogPriority::ERROR, "Scaffold is loaded from file, cannot recalculate metrics!");
		showPopup = false;
		return;
	}

	static int analysisScope = 0;
	static float tempVoxelSize = gen->measurementVoxelSize;
	static std::array<float, 2> xDims = { 0.0f, 5.0f };
	static std::array<float, 2> yDims = { 0.0f, 5.0f };
	static std::array<float, 2> zDims = { 0.0f, 5.0f };
	static ROI* selectedROI = nullptr;

	if (showPopup) {
		if (!ImGui::IsPopupOpen(popupName)) {
			ImGui::OpenPopup(popupName);
			std::array<float, 6> bds = gen->get_bounds();
			xDims[0] = bds[0];
			xDims[1] = bds[1];
			yDims[0] = bds[2];
			yDims[1] = bds[3];
			zDims[0] = bds[4];
			zDims[1] = bds[5];
		}
	}

	if (ImGui::BeginPopupModal(popupName, nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
		
		ImGui::SeparatorText("Parameters");
		
		ImGui::InputFloat("Voxel Size (mm)", &tempVoxelSize, 0.001f, 1.0f);
		ImGui::RadioButton("Analyze Entire Scaffold", &analysisScope, 0);
		ImGui::RadioButton("Analyze Inside ROI", &analysisScope, 1);
		
		if (analysisScope == 1){
		
			std::string previewValue = selectedROI ? selectedROI->name : "Choose...";

			if (ImGui::BeginCombo("##ROIChooser", previewValue.c_str()))
			{
				for (size_t i = 0; i < rois.size(); i++)
				{
					const bool isSelected = (selectedROI == rois[i].get());

					if (ImGui::Selectable(rois[i]->name.c_str(), isSelected))
					{
						selectedROI = rois[i].get();
					}

					if (isSelected)
					{
						ImGui::SetItemDefaultFocus();
					}
				}
				ImGui::EndCombo();
			}
		
		}
	
		ImGui::NewLine();
		if (ImGui::Button("Estimate")) {

			if (analysisScope == 0) {
				std::array<float, 6> bds = gen->get_bounds();
				gen->estimate_local_thickness(tempVoxelSize, bds, flag);
			}
			else {
				ROI* roi = static_cast<ROI*>(selectedROI);
				std::array<float, 6> roiBounds = roi->get_bounds();
				gen->estimate_local_thickness(tempVoxelSize, roiBounds, flag);
			}
			showPopup = false;
			selectedROI = nullptr;
			voxelSize = 0.05f;
		}
		
		ImGui::SameLine();
		if (ImGui::Button("Cancel")) {
			showPopup = false;
			voxelSize = 0.05f;
			selectedROI = nullptr;
			ImGui::CloseCurrentPopup();
		}

		ImGui::EndPopup();
	}
};

void myGUI::_action_estimate_tortuosity() {

	// get the selected scaffold
	if (selectedPanelObj.type != ObjectType::ScaffoldType) {
		logger.log(LogPriority::ERROR, "Select a scaffold from the left panel.");
		estimateTortuosity = false;
		return;
	}

	// get the scaffold
	GeneratorLewiner* scaffold = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);

	bool flag = false;

	if (scaffold->isLoadedFromFile) {
		logger.log(LogPriority::ERROR, "Scaffold is loaded from file, cannot recalculate metrics!");
		estimateTortuosity = false;
		return;
	}

	if (scaffold) {
		auto start = std::chrono::high_resolution_clock::now();

		flag = scaffold->estimate_tortuosity(tortuosityVoxelSize);
		
		auto end = std::chrono::high_resolution_clock::now();

		auto duration_ms = std::chrono::duration_cast<std::chrono::seconds>(end - start);

		std::ostringstream oss;
		oss << std::fixed << std::setprecision(3) // Set precision to 3 decimal places
			<< duration_ms.count()   // Convert ms to seconds
			<< " seconds!";
		if (flag) {
			logger.log(
				LogPriority::SUCCESS, "Estimated tortuosity in " + oss.str());
		}
		else {
			logger.log(
				LogPriority::ERROR, "Error in estimating tortuosity!");
		}
		
	}
	
	estimateTortuosity = false;
};

void myGUI::_action_estimate_anisotropy() {

	// get the selected scaffold
	if (selectedPanelObj.type != ObjectType::ScaffoldType) {
		logger.log(LogPriority::ERROR, "Select a scaffold from the left panel.");
		estimateAnisotropy = false;
		return;
	}

	// get the scaffold
	GeneratorLewiner* scaffold = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);

	// pass settings
	if(scaffold->isLoadedFromFile){
		logger.log(LogPriority::ERROR, "Scaffold is loaded from file, cannot recalculate metrics!");
		estimateAnisotropy = false;
		return;
	}

	if (estimateAnisotropy) {
		if (!ImGui::IsPopupOpen("Anisotropy Settings")) {
			ImGui::OpenPopup("Anisotropy Settings");
		}
	}
	
	if (ImGui::BeginPopupModal("Anisotropy Settings", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {

		static float tempVoxelSize = scaffold->measurementVoxelSize;
		static int analysisScope = 0;
		static ROI* selectedROI = nullptr;

		ImGui::SeparatorText("Parameters");

		ImGui::InputFloat("Voxel Size (mm)", &tempVoxelSize, 0.001f, 1.0f);
		ImGui::RadioButton("Analyze Entire Scaffold", &analysisScope, 0);
		ImGui::RadioButton("Analyze Inside ROI", &analysisScope, 1);

		if (analysisScope == 1) {

			std::string previewValue = selectedROI ? selectedROI->name : "Choose...";

			if (ImGui::BeginCombo("##ROIChooser", previewValue.c_str()))
			{
				for (size_t i = 0; i < rois.size(); i++)
				{
					const bool isSelected = (selectedROI == rois[i].get());

					if (ImGui::Selectable(rois[i]->name.c_str(), isSelected))
					{
						selectedROI = rois[i].get();
					}

					if (isSelected)
					{
						ImGui::SetItemDefaultFocus();
					}
				}
				ImGui::EndCombo();
			}
		}

		ImGui::NewLine();
		if (ImGui::Button("Estimate")) {

			scaffold->estimate_anisotropy(tempVoxelSize, daMinDirections, daMinLines, daFormulaIdx, selectedROI);
			estimateAnisotropy = false;
			selectedROI = nullptr;
		}

		ImGui::SameLine();
		if (ImGui::Button("Cancel")) {
			estimateAnisotropy = false;
			selectedROI = nullptr;
			ImGui::CloseCurrentPopup();
		}

		ImGui::EndPopup();
	}

};

void myGUI::_action_estimate_connectivity_density() {

	// get the selected scaffold
	if (selectedPanelObj.type != ObjectType::ScaffoldType) {
		logger.log(LogPriority::ERROR, "Select a scaffold from the left panel.");
		estimateConnectivityDensity = false;
		return;
	}

	// get the scaffold
	GeneratorLewiner* scaffold = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);

	// pass settings
	if(scaffold->isLoadedFromFile){
		logger.log(LogPriority::ERROR, "Scaffold is loaded from file, cannot recalculate metrics!");
		estimateConnectivityDensity = false;
		return;
	}

	if (estimateConnectivityDensity) {
		if (!ImGui::IsPopupOpen("Connectivity Density Settings")) {
			ImGui::OpenPopup("Connectivity Density Settings");
		}
	}
	
	if (ImGui::BeginPopupModal("Connectivity Density Settings", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {

		static int selectedMethod = 1;
		static float tempVoxelSize = 0.05f;

		const char* formulaOptions[] = { "Mesh - based", "Voxel - based" };

		ImGui::Combo("Algorithm", &selectedMethod, formulaOptions, IM_ARRAYSIZE(formulaOptions));

		if (selectedMethod == 1){
			ImGui::InputFloat("Voxel Size (mm)", &tempVoxelSize, 0.001f, 1.0f);
		}

		ImGui::NewLine();
		if (ImGui::Button("Estimate")) {

			if (selectedMethod == 0){
				scaffold->estimate_connectivity_density();
			}
			else if (selectedMethod == 1){
				scaffold->estimate_connectivity_density_voxel(tempVoxelSize);
			}
			estimateConnectivityDensity = false;
			selectedMethod = 1;
		}

		ImGui::SameLine();
		if (ImGui::Button("Cancel")) {
			estimateConnectivityDensity = false;
			selectedMethod = 1;
			ImGui::CloseCurrentPopup();
		}

		ImGui::EndPopup();

	}
};

void myGUI::_action_estimate_trabecular_number() {

	// get the selected scaffold
	if (selectedPanelObj.type != ObjectType::ScaffoldType) {
		logger.log(LogPriority::ERROR, "Select a scaffold from the left panel.");
		estimateTrabecularNr = false;
		return;
	}

	// get the scaffold
	GeneratorLewiner* scaffold = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);

	// pass settings
	if (scaffold->isLoadedFromFile) {
		logger.log(LogPriority::ERROR, "Scaffold is loaded from file, cannot recalculate metrics!");
		estimateTrabecularNr = false;
		return;
	}

	if (estimateTrabecularNr) {
		if (!ImGui::IsPopupOpen("Trabecular Number Settings")) {
			ImGui::OpenPopup("Trabecular Number Settings");
		}
	}

	if (ImGui::BeginPopupModal("Trabecular Number Settings", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {

		static float tempVoxelSize = 0.05f;
		static int analysisScope = 0;
		static ROI* selectedROI = nullptr;

		ImGui::SeparatorText("Parameters");

		ImGui::InputFloat("Voxel Size (mm)", &tempVoxelSize, 0.001f, 1.0f);
		
		const char* formulaOptions[] = { "Mean Intercept Length (DDA)", "Derived Proxy (BV/TV / Tb.Th)", "Derived Proxy (1 / Tb.Sp.)" };
		ImGui::Combo("Algorithm", &trabecularNrFormula, formulaOptions, IM_ARRAYSIZE(formulaOptions));
		
		ImGui::RadioButton("Analyze Entire Scaffold", &analysisScope, 0);
		ImGui::RadioButton("Analyze Inside ROI", &analysisScope, 1);

		if (analysisScope == 1) {

			std::string previewValue = selectedROI ? selectedROI->name : "Choose...";

			if (ImGui::BeginCombo("##ROIChooser", previewValue.c_str()))
			{
				for (size_t i = 0; i < rois.size(); i++)
				{
					const bool isSelected = (selectedROI == rois[i].get());

					if (ImGui::Selectable(rois[i]->name.c_str(), isSelected))
					{
						selectedROI = rois[i].get();
					}

					if (isSelected)
					{
						ImGui::SetItemDefaultFocus();
					}
				}
				ImGui::EndCombo();
			}
		}

		ImGui::NewLine();
		if (ImGui::Button("Estimate")) {
			scaffold->estimate_trabecular_number(tempVoxelSize, trabecularNrFormula, daMinDirections, daMinLines, selectedROI);
			estimateTrabecularNr = false;
			selectedROI = nullptr;
			trabecularNrFormula = 0;
		}

		ImGui::SameLine();
		if (ImGui::Button("Cancel")) {
			estimateTrabecularNr = false;
			trabecularNrFormula = 0;
			selectedROI = nullptr;
			ImGui::CloseCurrentPopup();
		}

		ImGui::EndPopup();
	}
};



void myGUI::_action_estimate_pore_network() {

	// get the selected scaffold
	if (selectedPanelObj.type != ObjectType::ScaffoldType) {
		logger.log(LogPriority::ERROR, "Select a scaffold from the left panel.");
		estimatePoreNetwork = false;
		return;
	}

	// get the scaffold
	GeneratorLewiner* scaffold = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);
	if (scaffold->isLoadedFromFile) {
		logger.log(LogPriority::ERROR, "Scaffold is loaded from file, cannot recalculate metrics!");
		estimatePoreNetwork = false;
		return;
	}

	if (scaffold) {
		scaffold->estimate_connectivity_network();
	}
};

void myGUI::_action_estimate_smi() {

	// get the selected scaffold
	if (selectedPanelObj.type != ObjectType::ScaffoldType) {
		logger.log(LogPriority::ERROR, "Select a scaffold from the left panel.");
		estimateSmi = false;
		return;
	}

	// get the scaffold
	GeneratorLewiner* scaffold = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);

	// pass settings
	if (scaffold->isLoadedFromFile) {
		logger.log(LogPriority::ERROR, "Scaffold is loaded from file, cannot recalculate metrics!");
		estimateSmi = false;
		return;
	}

	if (estimateSmi) {
		if (!ImGui::IsPopupOpen("Estimate Structural Model Index")) {
			ImGui::OpenPopup("Estimate Structural Model Index");
		}
	}
	if (ImGui::BeginPopupModal("Estimate Structural Model Index", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {

		static float dilation = 0.0f;

		// SMI needs dBS/dr, so the offset must stay far below the strut
		// thickness or the forward difference stops being a derivative.
		// 0 lets estimate_smi pick a scale-relative default (0.05 * voxel).
		ImGui::InputFloat("Dilation (0 = auto)", &dilation, 0.001f, 0.01f, "%.4f");
		ImGui::SetItemTooltip(
			"Surface offset used for dBS/dr. Leave at 0 for an automatic, "
			"scale-relative value. Must be much smaller than the strut thickness.");

		if (dilation < 0.0f) dilation = 0.0f;

		ImGui::NewLine();
		if (ImGui::Button("Estimate")) {
			scaffold->estimate_smi(dilation);
			estimateSmi = false;
			dilation = 0.0f;
		}

		ImGui::SameLine();
		if (ImGui::Button("Cancel")) {
			estimateSmi = false;
			dilation = 0.0f;
			ImGui::CloseCurrentPopup();
		}

		ImGui::EndPopup();

	}

};

void myGUI::_render_translate_panel(const char* popupName, bool& showPopup) {

	if (!selectedSceneObj) {
		showPopup = false;
		return;
	}

	if (ImGui::Begin(popupName, NULL)) {
	
		static float tempX = { 0.0f };
		static float tempY = { 0.0f };
		static float tempZ = { 0.0f };

		ImGui::InputFloat("Translate X", &tempX, 0.01f, 1000.0f);
		ImGui::InputFloat("Translate Y", &tempY, 0.01f, 1000.0f);
		ImGui::InputFloat("Translate Z", &tempZ, 0.01f, 1000.0f);

		// grab the selected panel
		GeneratorLewiner* gen = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);
		gen->translateVec.x = tempX;
		gen->translateVec.y = tempY;
		gen->translateVec.z = tempZ;

		if (ImGui::Button("Close")) {
			showPopup = false;
		}

		ImGui::End();
	
	};

};

void myGUI::_render_scale_panel(const char* popupName, bool& showPopup) {
	
	if (!selectedSceneObj) {
		showPopup = false;
		return;
	}

	if (ImGui::Begin(popupName, NULL)) {

		static float tempX = { 1.0f };
		static float tempY = { 1.0f };
		static float tempZ = { 1.0f };

		ImGui::InputFloat("Scale X", &tempX, 0.001f, 1000.0f);
		ImGui::InputFloat("Scale Y", &tempY, 0.001f, 1000.0f);
		ImGui::InputFloat("Scale Z", &tempZ, 0.001f, 1000.0f);

		// grab the selected panel
		GeneratorLewiner* gen = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);
		gen->scaleVec.x = tempX;
		gen->scaleVec.y = tempY;
		gen->scaleVec.z = tempZ;

		if (ImGui::Button("Scale")) {
			gen->apply_scale();
		}
		ImGui::SameLine();
		
		if (ImGui::Button("Close")) {
			showPopup = false;
		}

		ImGui::End();

	};
};

void myGUI::_render_roi_creator(const char* popupName, bool& showPopup) {
	
	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}

	// always centered
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

	if (ImGui::BeginPopupModal(popupName, NULL, ImGuiWindowFlags_AlwaysAutoResize))
	{
		static std::string buffer = "";
		ImGui::InputText("Name", &buffer);

		// display here the settings for a box set
		static Vec3 size{ 10.0f, 10.0f, 10.0f };;
		static Vec3 origin;

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Width (X)", &size.x, 0.1f, 100.0f, "%.3f");

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Height (Y)", &size.y, 0.1f, 100.0f, "%.3f");

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Depth (Z)", &size.z, 0.1f, 100.0f, "%.3f");
		
		ImGui::InputFloat3("Center", origin);
		ImGui::SetItemTooltip("Center of ROI");

		ImGui::Separator();
		ImGui::SameLine();

		if (ImGui::Button("Create")) {
			std::shared_ptr<ROI> roi = std::make_shared<ROI>(size, origin);
			if (buffer.empty()) {
				roi->name = "roi" + std::to_string(rois.size() + 1);
			}
			else {
				roi->name = buffer;
			};

			rois.push_back(roi);

			size = { 10.0f, 10.0f, 10.0f };
			buffer[0] = '\0';
			origin = Vec3{ 0.0f, 0.0f, 0.0f };
			showPopup = false;
			showPopup = false;
		}
	
		ImGui::SameLine();
		if (ImGui::Button("Cancel")) {
			buffer[0] = '\0';
			size = { 10.0f, 10.0f, 10.0f };
			origin = Vec3{ 0.0f, 0.0f, 0.0f };
			showPopup = false;
		};
		ImGui::EndPopup();
	}
};

void myGUI::_render_aniso_source_creator(
	const char* popupName, bool& showPopup) {
	
	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}

	// always centered
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

	if (ImGui::BeginPopupModal(popupName, NULL, ImGuiWindowFlags_AlwaysAutoResize))
	{
		static int selected = -1;
		static std::vector<std::shared_ptr<AnisotropySource>> tempSources;
		for (int i{ 0 }; i < tempSources.size(); i++) {
			std::string name = "Anisotropy Source" + std::to_string(i + 1);
			if (ImGui::TreeNode(name.c_str())){
				selected = i;
				tempSources[i]->render_properties();
				ImGui::TreePop();
			}
			else{
				selected = -1;
			}	
		}
		
		ImGui::Separator();
		
		if (ImGui::Button("Add")){
			std::shared_ptr<AnisotropySource> source = 
			std::make_shared<AnisotropySource>();
			source->update_metric();
			source->update_model();
			if (source->name.empty()){
				source->name = "Anisotropy Source " + std::to_string(tempSources.size() + 1);
			}
			tempSources.push_back(std::move(source));
			selected = -1;
		}		

		ImGui::SameLine();
		if (ImGui::Button("Delete Selected")){
			if(selected != -1){
				auto it = tempSources.erase(tempSources.begin() + selected);
			}
		}

		ImGui::SameLine();
		if (ImGui::Button("Cancel")) {
			tempSources.clear();
			showPopup = false;
			selected = -1;
		};

		// if ok push to the anisosources
		ImGui::SameLine();
		if (ImGui::Button("Ok")) {
			anisoSources.insert(
				anisoSources.end(),
				tempSources.begin(),
				tempSources.end());

			showPopup = false;
			selected = -1;
		};

		ImGui::EndPopup();
	}
};

void myGUI::_render_roi_cutter(const char* popupName, bool& showPopup) {

	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}

	// always centered
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

	if (ImGui::BeginPopupModal(popupName, NULL, ImGuiWindowFlags_AlwaysAutoResize))
	{
		ImGui::Text("Select the Region of Interest to extract:");
		ImGui::Separator();
		ImGui::Spacing();

		static ROI* selectedROI = nullptr;
		std::string previewValue = selectedROI ? selectedROI->name : "Select ROI...";

		if (ImGui::BeginCombo("##ROIChooser", previewValue.c_str()))
		{
			for (size_t i = 0; i < rois.size(); i++)
			{
				const bool isSelected = (selectedROI == rois[i].get());

				if (ImGui::Selectable(rois[i]->name.c_str(), isSelected))
				{
					selectedROI = rois[i].get();
				}

				if (isSelected)
				{
					ImGui::SetItemDefaultFocus();
				}
			}
			ImGui::EndCombo();
		}

		ImGui::Spacing();
		ImGui::Separator();
		ImGui::Spacing();

		bool hasSelection = (selectedROI != nullptr);
		if (!hasSelection) ImGui::BeginDisabled();

		if (ImGui::Button("Cut")) {

			auto* gen = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);
			
			if (gen->isLoadedFromFile) {
				logger.log(LogPriority::ERROR, "The Scaffold is loaded from a scalar field that is no more available!");
			}
			std::unique_ptr<GeneratorLewiner> newMesh;
			newMesh = gen->extract_from_ROI(selectedROI);
			if (selectedROI) {
				newMesh->name = gen->name + "_" + selectedROI->name;
			}
			else {
				newMesh->name = "Scaffold" + scaffolds.size() + 1;
			}
			scaffolds.push_back(std::move(newMesh));
			selectedPanelObj.ptr = scaffolds.back().get();
			selectedPanelObj.type = ObjectType::ScaffoldType;
			selectedSceneObj = scaffolds.back().get();
			selectedROI = nullptr;
			showPopup = false;
		}
		if(!hasSelection) ImGui::EndDisabled();
		ImGui::SameLine();
		if (ImGui::Button("Cancel")) {
			selectedROI = nullptr;
			showPopup = false;
		};
		ImGui::EndPopup();
	}
};

void myGUI::_render_taubin_smooth_panel(const char* popupName, bool& showPopup) {

	if (!selectedSceneObj) {
		showPopup = false;
		return;
	}

	// apply to the vertices of the selected object
	auto model = static_cast<GeneratorLewiner*>(selectedSceneObj);

	if (ImGui::Begin(popupName, &showPopup, ImGuiWindowFlags_AlwaysAutoResize)) {

		static int iter = 1;
		static float lambda = 0.5f;
		static float mu = -0.3f;

		ImGui::InputInt("Iterations", &iter, 1);

		ImGui::InputFloat("Lambda", &lambda, 0.01f, 10.0f);

		ImGui::InputFloat("Mu", &mu, 0.01f, 10.0f);
	
		if (ImGui::Button("Apply")) {
			model->apply_taubin_smooth(iter, lambda, mu);
			// model->smooth_scalar_field_taubin(iter, lambda, mu);
			// std::shared_ptr<IContainer> parentCon = model->container.lock();
			// model->marching_cubes();
			// model->estimate_metrics(*parentCon);
			// model->update_render();
		}

		if (ImGui::Button("Cancel")) {
			showPopup = false;
		}

		ImGui::End();
	}
};

void myGUI::_render_simplify_mesh_panel(const char* popupName, bool& showPopup){
	
	if (!selectedSceneObj) {
		showPopup = false;
		return;
	}

	// apply to the vertices of the selected object
	auto model = static_cast<GeneratorLewiner*>(selectedSceneObj);

	if (ImGui::Begin(popupName, &showPopup, ImGuiWindowFlags_AlwaysAutoResize)) {
		static float targetRatio = 0.25;
		static float maxError = 10.0;
		static int targetTris = 5000;
		static bool preventFlip = true;
		static bool  recomputeNormals = true;
		static bool lockBoundary = true;

		ImGui::InputInt("Target Triangles", &targetTris, 1);
		ImGui::InputFloat("Percentage", &targetRatio);
		ImGui::Checkbox("Prevent Flipping Triangles", &preventFlip);
		ImGui::Checkbox("Lock Boundary Vertices", &lockBoundary);
		ImGui::Checkbox("Recompute Normals", &recomputeNormals);

		if (ImGui::Button("Apply")) {

			mesh_simplify::Options options;
			options.targetTriangles  = targetTris > 0 ? (size_t)targetTris : 0;
			options.targetRatio      = targetRatio;
			options.maxError         = maxError;
			options.preventFlips     = preventFlip;
			options.recomputeNormals = recomputeNormals;
			options.lockBoundary     = lockBoundary;
			model->apply_mesh_simplification(options);

			showPopup = false;
		}
		ImGui::SameLine();
		if (ImGui::Button("Cancel")) {
			showPopup = false;
		}

		ImGui::End();

	}
	
};


void myGUI::_draw_selected_box() {

	if (!selectedSceneObj) return;

	auto md = static_cast<GeneratorLewiner*>(selectedSceneObj);
	if (md->hidden) return;

	boxShader.use();
	glDepthFunc(GL_LEQUAL);
	glLineWidth(1.0f);

	Aabb aabb = md->get_aabb();

	Vec3 size = aabb.pMax - aabb.pMin;
	Vec3 center = (aabb.pMin + aabb.pMax) * 0.5f;

	glm::mat4 mt = glm::translate(glm::mat4(1.0f), glm::vec3(
		center.x + md->translateVec.x,
		center.y + md->translateVec.y,
		center.z + md->translateVec.z
	));
	mt = glm::scale(mt, glm::vec3(size.x, size.y, size.z));
	uniManager.setUniform(boxShader, "projection", projection);
	uniManager.setUniform(boxShader, "view", view);
	uniManager.setUniform(boxShader, "model", mt);
	uniManager.setUniform(
		boxShader,
		"boxColor",
		0.0f, 1.0f, 0.0f, 1.0f
	);

	box->draw();
};

void myGUI::_draw_axes_lines() {

	float infinity = 10000.0f;

	glDisable(GL_DEPTH_CLAMP);
	glDepthFunc(GL_LEQUAL);
	glLineWidth(0.5f);
	lineShader.use();
	uniManager.setUniform(lineShader, "projection", projection);
	uniManager.setUniform(lineShader, "view", view);
	uniManager.setUniform(lineShader, "model", glm::mat4(1.0f));
	uniManager.setUniform(
		lineShader,
		"lineColor",
		1.0f, 0.0f, 0.0f, 1.0f
	);
	glm::mat4 modelX = glm::scale(glm::mat4(1.0f), glm::vec3(infinity, 1.0f, 1.0f));
	uniManager.setUniform(lineShader, "model", modelX);
	lineX->draw();

	// ---------------------
	glm::mat4 modelY = glm::scale(glm::mat4(1.0f), glm::vec3(1.0f, infinity, 1.0f));
	uniManager.setUniform(lineShader, "model", modelY);
	uniManager.setUniform(
		lineShader,
		"lineColor",
		0.0f, 1.0f, 0.0f, 1.0f
	);
	lineY->draw();

	// ------------------
	glm::mat4 modelZ = glm::scale(glm::mat4(1.0f), glm::vec3(1.0f, 1.0f, infinity));
	uniManager.setUniform(
		lineShader,
		"lineColor",
		0.0f, 0.0f, 1.0f, 1.0f
	);
	uniManager.setUniform(lineShader, "model", modelZ);
	lineZ->draw();
	glEnable(GL_DEPTH_CLAMP);

};

// ---------------------------------------------------------------------------------------------------------------

bool create_single_button_textured(const char* name, bool& flag, const std::string& tooltip, const GLuint textureId, bool enabled) {

	bool pressed = false;

	// Disable if not enabled
	if (!enabled) {
		ImGui::BeginDisabled();
	}

	if (ImGui::ImageButton(
		name, (ImTextureID)(intptr_t)textureId, ImVec2(40.0, 40.0)
	)) {
		flag = !flag;
		pressed = true;
	}

	if (ImGui::IsItemHovered() && !tooltip.empty()) {
		ImGui::SetTooltip(tooltip.c_str());
	}

	if (!enabled) {
		ImGui::EndDisabled();
	}

	return pressed;
};

// @brief Function to load texture from file, acquired by https://github.com/ocornut/imgui/wiki/Image-Loading-and-Displaying-Examples
bool load_texture_from_file(const char* fileName, GLuint* outTexture, int* outWidth, int* outHeight)
{
	FILE* f = fopen(fileName, "rb");
	if (f == NULL)
		return false;
	fseek(f, 0, SEEK_END);
	size_t fileSize = (size_t)ftell(f);
	if (fileSize == -1)
		return false;
	fseek(f, 0, SEEK_SET);
	void* fileData = IM_ALLOC(fileSize);
	fread(fileData, 1, fileSize, f);
	fclose(f);
	bool ret = load_texture_from_memory(fileData, fileSize, outTexture, outWidth, outHeight);
	IM_FREE(fileData);
	return ret;
};

// @brief function to create OpenGL texture, acquired by 
// https://github.com/ocornut/imgui/wiki/Image-Loading-and-Displaying-Examples
bool load_texture_from_memory(const void* data, size_t dataSize, GLuint* outTexture, int* outWidth, int* outHeight) {

	int width, height;

	// use stbi_load
	unsigned char* imgData = stbi_load_from_memory((const unsigned char*)data, int(dataSize), &width, &height, NULL, 4);

	// check if there are data
	if (!imgData) {
		std::cerr << "No data found" << std::endl;
		return false;
	}

	// opengl stuff -> from opengl tutorial
	unsigned int texture;
	glGenTextures(1, &texture);

	// bind texture
	glBindTexture(GL_TEXTURE_2D, texture);

	// texture wrapping filtering options
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	// create image
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, imgData);
	glGenerateMipmap(GL_TEXTURE_2D);
	stbi_image_free(imgData);

	*outTexture = texture;
	*outWidth = width;
	*outHeight = height;

	return true;
};

bool create_button_textured(const char* name, bool& flag, const std::string& tooltip, const GLuint textureId, bool enabled) {

	bool clicked = false;

	// Disable if not enabled
	if (!enabled) {
		ImGui::BeginDisabled();
	}

	// this will keep the state 
	const bool wasOn = flag;

	// Push style if it was active
	if (wasOn) {
		ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.4f, 0.8f, 0.4f, 1.0f));
		ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.5f, 0.9f, 0.5f, 1.0f));
		ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.3f, 0.7f, 0.3f, 1.0f));
	}

	if (ImGui::ImageButton(
		name, (ImTextureID)(intptr_t)textureId, ImVec2(40.0, 40.0)
	)) {
		flag = !flag;
		clicked = true;
	}

	if (ImGui::IsItemHovered() && !tooltip.empty()) {
		ImGui::SetTooltip(tooltip.c_str());
	}

	// Pop only if pushed
	if (wasOn) {
		ImGui::PopStyleColor(3);
	}

	if (!enabled) {
		ImGui::EndDisabled();
	}

	return clicked && !wasOn;
};
