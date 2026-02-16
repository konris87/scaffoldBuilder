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
		add_log(LogPriority::ERROR, "Failed to initialize GLAD!");
		//std::cerr << "Failed to initialize GLAD" << std::endl;
		return;
	}

	// Now you can check the OpenGL version
	std::string version = std::string((const char*)glGetString(GL_VERSION));
	std::string glsVersion = std::string((const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));

	// Print the OpenGL version
	add_log(LogPriority::INFO, "OpenGL Version: " + version);
	add_log(LogPriority::INFO, "GLSL Version: " + glsVersion);

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
		std::filesystem::absolute(std::filesystem::path("./share/shaders/mainVShader.vertexshader")).string().c_str(),
		std::filesystem::absolute(std::filesystem::path("./share/shaders/mainFShader.fragmentshader")).string().c_str(),
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
		//"shaders/frameVShader.vs",
		//"shaders/frameFShader.fs",
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
	//load_texture_from_file("./share/textures/abstractContainer.png", &abstractContainerTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/randomSeeds.png", &randomSeedTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/uniformSeeds.png", &uniformSeedTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/variedSeeds.png", &variedSeedTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/scaffoldCreator.png", &createScaffoldTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/thickness.png", &localThicknessTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/separation.png", &separationTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/tortuosity.png", &tortuosityTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/network.png", &poreNetworkTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/update.png", &updateScaffoldTexture, &dummyW, &dummyH);
	load_texture_from_file("./share/textures/translate.png", &translateTexture, &dummyW, &dummyH);
};

void myGUI::_init_imgui() {

	//std::cout << "Starting imgui" << std::endl;

	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();

	io = ImGui::GetIO(); (void)io;
	io.Fonts->AddFontFromFileTTF("./share/fonts/DroidSans.ttf", 18.0f);
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

	while (!glfwWindowShouldClose(window))
	{

		//glfwGetFramebufferSize(window, &width, &height);

		glfwPollEvents();

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

		// -------------------------------------------------------------------------------------------
		// First pass -> objects without alpha

		glEnable(GL_DEPTH_TEST);
		glDepthMask(GL_TRUE);     // <--- IMPORTANT: Opaque objects WRITE to depth buffer
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
					uniManager.setUniform(scaffoldShader, "cutPlaneCoeffs", planeCoeffs);
					gen->draw();
					glEnable(GL_CULL_FACE);
				}
				else {
					uniManager.setUniform(scaffoldShader, "cutPlane", 0);
					gen->draw();
				}
				glDisable(GL_POLYGON_OFFSET_FILL);

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
					uniManager.setUniform(edgeShader, "model", model);
					uniManager.setUniform(edgeShader, "edgeColor", 0.0f, 0.0f, 0.0f, 1.0f);
					if (showScaffold) {
						gen->draw_edges();
					}
				}

				if (showNormals && gen.get() == selectedSceneObj) {
					normalShader.use();
					uniManager.setUniform(normalShader, "projection", projection);
					uniManager.setUniform(normalShader, "view", view);
					uniManager.setUniform(normalShader, "model", model);
					uniManager.setUniform(normalShader, "normalColor", normalColor[0], normalColor[1], normalColor[2]);
					if (showScaffold) {
						gen->draw();
					}
				}

				if (showTortuosityPath && gen.get() == selectedSceneObj) {
					glDepthFunc(GL_LEQUAL);
					glLineWidth(renderSettings.poreNetworkLineSize);
					edgeShader.use();
					uniManager.setUniform(edgeShader, "projection", projection);
					uniManager.setUniform(edgeShader, "view", view);
					uniManager.setUniform(edgeShader, "model", model);
					uniManager.setUniform(edgeShader, "cutPlane", 0);
					uniManager.setUniform(
						edgeShader,
						"edgeColor",
						renderSettings.poreNetworkColor[0],
						renderSettings.poreNetworkColor[1],
						renderSettings.poreNetworkColor[2],
						renderSettings.poreNetworkColor[3]
					);
					//gen->draw_tortuosity_path();				
				}

				if (showPoreNetwork && gen.get() == selectedSceneObj) {
					glDepthFunc(GL_LEQUAL);
					glLineWidth(renderSettings.poreNetworkLineSize);
					edgeShader.use();
					uniManager.setUniform(edgeShader, "projection", projection);
					uniManager.setUniform(edgeShader, "view", view);
					uniManager.setUniform(edgeShader, "model", model);
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
		}

		// containers
		for (const auto& con : containers) {
			if (con && !con->hidden) {
				boxShader.use();
				uniManager.setUniform(boxShader, "projection", projection);
				uniManager.setUniform(boxShader, "view", view);
				uniManager.setUniform(boxShader, "model", model);
				con->render();
			}
		}

		io = ImGui::GetIO();
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

				seedShader.use();
				uniManager.setUniform(seedShader, "projection", projection);
				uniManager.setUniform(seedShader, "view", view);
				uniManager.setUniform(seedShader, "seedSize", seedSize);
				uniManager.setUniform(seedShader, "seedColor", seedColor[0], seedColor[1], seedColor[2]);
				seedGen->draw();				
			}
		}

		_draw_selected_box();

		// ----------------------------------------------------------------------------
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
				uniManager.setUniform(scaffoldShader, "model", glm::mat4(1.0));
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

		if (showCutPlane) {
			cutShader.use();
			uniManager.setUniform(cutShader, "projection", projection);
			uniManager.setUniform(cutShader, "view", view);
			uniManager.setUniform(cutShader, "model", cutPlane->modelMatrix);
			uniManager.setUniform(cutShader, "minBounds", bounds[0], bounds[2], bounds[4]);
			uniManager.setUniform(cutShader, "maxBounds", bounds[1], bounds[3], bounds[5]);
			glDisable(GL_CULL_FACE);
			cutPlane->draw();
			glEnable(GL_CULL_FACE);
		}

		io = ImGui::GetIO();

		if (ImGui::IsKeyPressed(ImGuiKey::ImGuiKey_Delete) && !io.WantCaptureKeyboard) {
			
			add_log(LogPriority::INFO, "Deleting everything!");

			scaffolds.clear();
			containers.clear();
			seedGenerators.clear();
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

		ImGui::TableNextColumn();
		create_single_button_textured(
			"Create Random Seed Generator", showRandomSeedCreator, "create random seeds inside a container", randomSeedTexture);
		ImGui::SameLine();

		create_single_button_textured(
			"Create Uniform Seed Generator", showUniformSeedCreator, "Create seeds inside a container using uniform Poisson sampling", uniformSeedTexture);
		ImGui::SameLine();

		create_single_button_textured(
			"Create Varied Seed Generator", showVariedSeedCreator, "Create seeds inside a container using varied Poisson sampling", variedSeedTexture);
		ImGui::SameLine();

		ImGui::TableNextColumn();

		create_single_button_textured(
			"Create Scaffold", showScaffoldCreator,
			"Create the scaffold selecting a container and a generator", createScaffoldTexture
		);
		ImGui::SameLine();

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
			"Measure Porosity Network", estimatePoreNetwork,
			"Create a graph that visualizes the connectivity network and estimates the percentage of interconneted seeds.", poreNetworkTexture
		);

		ImGui::SameLine();

		ImGui::TableNextColumn();

		create_single_button_textured(
			"Update Current Scaffold", updateScaffold,
			"Update Current Scaffold.", updateScaffoldTexture
		);

		ImGui::SameLine();
		ImGui::TableNextColumn();

		create_button_textured("Translate Object", translateScaffold, "Translate selected object.", translateTexture);

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

	ImVec2 maxSize = ImVec2(width, height);
	ImVec2 minSize = ImVec2(500.0f, 500.0f);
	if (ImGuiFileDialog::Instance()->Display("Export Scaffold", ImGuiWindowFlags_NoCollapse, minSize, maxSize)) {
		if (ImGuiFileDialog::Instance()->IsOk()) { // action if OK
			scaffoldFilePath = ImGuiFileDialog::Instance()->GetFilePathName();
			scaffoldFileName = ImGuiFileDialog::Instance()->GetCurrentFileName();
			//_export_mesh();

			// get the active scaffold
			//Generator* gen = static_cast<Generator*>(selectedPanelObj.ptr);
			GeneratorLewiner* gen = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);

			std::filesystem::path filePath = scaffoldFilePath;

			if (filePath.extension() == ".stl") {

				gen->export_stl(scaffoldFilePath);

				add_log(LogPriority::SUCCESS, "Exported scaffold mesh to " + scaffoldFilePath);

			}
		}
		ImGuiFileDialog::Instance()->Close();
	}

	maxSize = ImVec2(width, height);
	minSize = ImVec2(500.0f, 500.0f);
	if (ImGuiFileDialog::Instance()->Display("Export Binary Scaffold", ImGuiWindowFlags_NoCollapse, minSize, maxSize)) {
		if (ImGuiFileDialog::Instance()->IsOk()) { // action if OK
			scaffoldFilePath = ImGuiFileDialog::Instance()->GetFilePathName();
			scaffoldFileName = ImGuiFileDialog::Instance()->GetCurrentFileName();
			std::string folder = ImGuiFileDialog::Instance()->GetCurrentPath();
			/*_export_binary_image(scaffoldFilePath,
					voxelSize, backVal, forVal);*/
					// get the active scaffold
			Generator* gen = static_cast<Generator*>(selectedPanelObj.ptr);

			std::filesystem::path filePath = scaffoldFilePath;

			if (filePath.extension() == ".mhd") {

				//std::string baseName = filePath.stem()
				//std::filesystem::path basePath = std::filesystem::absolute(filePath);
				//basePath.replace_extension("");

				gen->export_mhd(filePath, voxelSize, voxelBounds);

				add_log(LogPriority::SUCCESS, "Exported scaffold as a binary image to " + scaffoldFilePath);

			}

			else if (filePath.extension() == ".nrrd") {

				gen->export_nrrd(scaffoldFilePath, voxelSize, voxelBounds);

				add_log(LogPriority::SUCCESS, "Exported scaffold as a binary image to " + scaffoldFilePath);


			}
		}
		ImGuiFileDialog::Instance()->Close();
	}

	//_render_console();

	//ImGui::End();
};

void myGUI::_render_object_list() {

	if (ImGui::Begin("Objects", &showScaffoldList)) {

		bool open = ImGui::TreeNodeEx("Scaffolds", ImGuiDockNodeFlags_None | ImGuiTreeNodeFlags_DefaultOpen);
		if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort | ImGuiHoveredFlags_NoSharedDelay)) {
			ImGui::SetTooltip("Scaffold objects");
		}
		if (open) {
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
					if (ImGui::MenuItem("Delete")) {
						auto it = scaffolds.erase(scaffolds.begin() + i);
						selectedSceneObj = nullptr;
						selectedPanelObj.ptr = nullptr;
						selectedPanelObj.type = ObjectType::NoneType;
						ImGui::EndPopup();
						break;
					}
					if (ImGui::MenuItem("Export as Mesh")) {
						showGeometryExportWindow = true;
						IGFD::FileDialogConfig config;
						config.path = "../data";
						ImGuiFileDialog::Instance()->OpenDialog("Export Scaffold", "Export Scaffold Geometry", ".stl, .vtk", config);
					}
					if (ImGui::MenuItem("Export as Image")) {
						showBinaryImageWindow = true;
					}

					ImGui::EndPopup();
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
			ImGui::TreePop();

		}
		_render_container_list();

		_render_seed_generator_list();
	}

	ImGui::End();
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
			case ObjectType::RandomGeneratorType:{
				_render_random_seed_generator_properties();
				break;
			}
			case ObjectType::UniformGeneratorType: {
				_render_uniform_seed_generator_properties();
				break;
			}
			case ObjectType::VariedGeneratorType: {
				_render_varied_seed_generator_properties();
				break;
			}
			case ObjectType::ScaffoldType: {
				_render_scaffold_properties();
				break;
			}
			case ObjectType::NoneType: {
				//std::cout << "None" << std::endl;
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
				ImGui::EndPopup();
			}
			if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && ImGui::IsItemHovered()) {
				selectedPanelObj.ptr = gen.get();
				selectedPanelObj.type = gen->type;
			}
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
		for (int i = 0; i < containers.size(); ++i) {

			auto& con = containers[i];

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
				ImGui::EndPopup();
			}

			if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && ImGui::IsItemHovered()) {
				selectedPanelObj.ptr = con.get();
				selectedPanelObj.type = con->get_type();
			}
		}
		ImGui::TreePop();
	}

};

void myGUI::write_settings() {

	nlohmann::json settings;

	settings["Scaffold"]["version"] = version;
	settings["Scaffold"]["Domain"]["xMin"] = 0.0;
	settings["Scaffold"]["Domain"]["xMax"] = xDim;
	settings["Scaffold"]["Domain"]["yMin"] = 0.0;
	settings["Scaffold"]["Domain"]["yMax"] = yDim;
	settings["Scaffold"]["Domain"]["zMin"] = 0.0;
	settings["Scaffold"]["Domain"]["zMax"] = zDim;
	settings["Scaffold"]["Pores"]["genOption"] = generateOption;
	settings["Scaffold"]["Pores"]["poreNr"] = seedNr;
	settings["Scaffold"]["thickness"] = thickness;
	settings["Scaffold"]["regSteps"] = regSteps;
	settings["Voro"]["nX"] = 6;
	settings["Voro"]["nY"] = 6;
	settings["Voro"]["nZ"] = 6;

	std::ofstream file("settings.json");
	file << settings.dump(4);
	file.close();

};

std::string myGUI::_get_glsl_version() {
	// Get OpenGL version
	const char* glVersion = reinterpret_cast<const char*>(glGetString(GL_VERSION));

	if (!glVersion) {
		std::cerr << "Error: Could not retrieve OpenGL version." << std::endl;
		return "#version 130";  // Default fallback
	}

	//std::cout << "OpenGL version: " << glVersion << std::endl;
	//add_log(LogPriority::INFO, "OpenGL version: " + std::string(glVersion));

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
		
		if (scrollToBottom) {
			ImGui::SetScrollHereY(1.0f);  // scroll to bottom
			scrollToBottom = false;
		}
		
		ImGui::EndChild();

		if (ImGui::Button("Clear")) {
			logger.clear();
			add_log(LogPriority::INFO, "Cleared.");
		}
	}

	ImGui::End();

}

void myGUI::add_log(LogPriority priority, const std::string& message) {

	if (priority == LogPriority::SUCCESS) {
		logColor = { 0.0f, 1.0f, 0.0f, 1.0f };
	}

	else if (priority == LogPriority::ERROR) {
		logColor = { 1.0f, 0.0f, 0.0f, 1.0f };
	}

	else {
		logColor = { 1.0f, 1.0f, 1.0f, 1.0f };
	}

	logger.log(priority, message, logColor);

	scrollToBottom = true;
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

void myGUI::_render_mesh_settings() {

	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
	ImGui::SetNextWindowSize(ImVec2(400, 300), ImGuiCond_Appearing);

	if (ImGui::Begin("Display Settings", &showDisplayMeshSettingsWin, ImGuiWindowFlags_AlwaysAutoResize)) {

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
	uniManager.setUniform(frameShader, "outColor", glm::vec4(0.0f, 1.0f, 0.0f, 1.0f));
	//uniManager.setUniform(frameShader, "model", glm::scale(glm::mat4(1.0), glm::vec3(2, 2, 2)));
	zArrow->draw();
	uniManager.setUniform(frameShader, "model", glm::rotate(glm::mat4(1.0f), glm::radians(-90.0f), glm::vec3(1.0f, 0.0f, 0.0f)));
	uniManager.setUniform(frameShader, "outColor", glm::vec4(0.0f, 0.0f, 1.0f, 1.0f));
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
				add_log(LogPriority::INFO, "Settings saved.");
			}
			if (ImGui::MenuItem("Load Scaffold", "load mesh file by specifying the path")) {
				// Call your function to load a model mesh here
				IGFD::FileDialogConfig config;
				config.path = "..//data";
				ImGuiFileDialog::Instance()->OpenDialog("ChooseFileDlgKey", "Choose File", ".stl, .vtk", config);

				loadedMesh = scaffold;
				//add_log("Model mesh loaded.");
			}

			if (ImGui::MenuItem("Load Bone", "load bone mesh file by specifying the path")) {
				// Call your function to load a model mesh here
				IGFD::FileDialogConfig config;
				config.path = "..//data";
				ImGuiFileDialog::Instance()->OpenDialog("ChooseFileDlgKey", "Choose File", ".stl", config);
				//add_log("Model mesh loaded.");
				loadedMesh = bone;
			}

			ImGui::SameLine(); help_marker("So far only .stl files are supported");

			if (ImGui::MenuItem("Save Scaffold as geometry", "Export scaffold geometry.")) {
				showGeometryExportWindow = true;
				IGFD::FileDialogConfig config;
				config.path = "../data";
				ImGuiFileDialog::Instance()->OpenDialog("Export Scaffold", "Export Scaffold Geometry", ".stl, .vtk", config);
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

		if (ImGui::BeginMenu("Tools")) {

			if (ImGui::IsItemHovered()) {
				ImGui::SetMouseCursor(ImGuiMouseCursor_Arrow);
			}
			if (ImGui::MenuItem("Cut With Plane")) {

				if (selectedPanelObj.type != ObjectType::ScaffoldType) {
					add_log(LogPriority::ERROR, "Selected a scaffold from the panel!");
				}
				else {
					showPlaneCutSettings = true;
					showCutPlane = true;
					cutPlane = std::make_unique<CutPlane>();
					GeneratorLewiner* gen = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);

					IContainer* con = static_cast<IContainer*>(gen->container);
					con->hidden = false;

					bounds = gen->get_bounds();
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

	if (showRandomSeedCreator) {
		_render_random_seed_creator("Random Seed Creator", showRandomSeedCreator);
	}

	if (showUniformSeedCreator) {
		_render_uniform_seed_creator("Uniform Seed Creator", showUniformSeedCreator);
	}

	if (showVariedSeedCreator) {
		_render_varied_seed_creator("Varied Seed Creator", showVariedSeedCreator);
	}

	if (showScaffoldCreator) {
		_render_scaffold_creator("Scaffold Creator", showScaffoldCreator);
	}

	if (measureThickness) {
		_local_thickness_measure("Local Thickness Measure", measureThickness, false);
	}

	if (measureSeparation) {
		_local_thickness_measure("Local Separation Measure", measureSeparation, true);
	}

	if (estimateTortuosity) {
		_action_estimate_tortuosity();
	}

	if (estimatePoreNetwork) {
		_action_estimate_pore_network();
	}

	if (showPlaneCutSettings) {
		_render_cutting_plane_settings("Cutting With Plane" , showPlaneCutSettings);
	}

	if (translateScaffold) {
		_render_translate_panel("Translate Object", translateScaffold);
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
			if (ImGui::Selectable("Linear Function", pickedItem == 4)) pickedItem = 4;
			if (ImGui::Selectable("Container", pickedItem == 5)) pickedItem = 5;
			if (ImGui::Selectable("Pore Network", pickedItem == 6)) pickedItem = 6;
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
					ImGui::InputFloat("Point Size", &seedSize, 0.001f, 1.0f, "%.3f");
				}
				else if (pickedItem == 2) {
					ImGui::Text("Grid Settings");
					ImGui::ColorEdit3("Grid Color", (float*)&gridColor);
					ImGui::Checkbox("Use Grid", &showGrid);
					//ImGui::SliderFloat("Point Size", &seedSize, 0.001f, 1.0f);
				}
				else if (pickedItem == 3) {
					ImGui::Text("Lighting Settings");
					ImGui::ColorEdit3("Light Color", (float*)&lightColor);
					ImGui::ColorEdit3("Font Color", (float*)&fontColor);
					ImGui::ColorEdit3("Normal Color", (float*)&normalColor);
					ImGui::InputFloat("Ambient Strength", &Ka, 0.01f, 100.0f, "%.3f");
					ImGui::InputFloat("Specular Strength", &Ks, 0.01f, 100.0f, "%.3f");
				}
				else if (pickedItem == 4) {
					ImGui::Text("Linear Function Settings");
					ImGui::InputDouble("Max Distance", &maxDist, 0.01f, 1.0f, "%.3f");
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
	//ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	//ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

	if (ImGui::Begin("Plane Cut Tool", NULL)) {
		
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
	ImGui::End();

};

void myGUI::_render_seed_generator() {

	if (ImGui::Begin("Seed Generator", NULL)) {

		ImGui::SeparatorText("Generator");
		ImGui::RadioButton("Random Seed Generator", &generateOption, 0);
		ImGui::SameLine(); help_marker("Generate random points inside the container");

		ImGui::RadioButton("Poisson3D", &generateOption, 1);
		ImGui::SameLine(); help_marker("Generate seeds with distance constraints");

		ImGui::SeparatorText("Settings");

		if (generateOption == 0) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputInt("Seeds", &seedNr, 1, 1000);
		}

		if (generateOption == 1) {
			//ImGui::InputFloat("Min Radius", &rMin);
			//ImGui::InputFloat("Max Radius", &rMax);

			//ImGui::SeparatorText("Varying Radius");
			//ImGui::RadioButton("Uniform Radius", &distFunc, 0);

			ImGui::RadioButton("Uniform Radius", &radiusOpt, 0);
			if (radiusOpt == 0) {
				ImGui::SetNextItemWidth(200);
				ImGui::InputFloat("Radius", &rMin, 0.01f, 1.0f, "%.3f");
			}

			ImGui::RadioButton("Varied Radius", &radiusOpt, 1);

			if (radiusOpt == 1) {
				ImGui::SetNextItemWidth(200);
				ImGui::InputFloat("Min Radius", &rMin, 0.01f, 1.0f, "%.3f");
				ImGui::SetNextItemWidth(200);
				ImGui::InputFloat("Max Radius", &rMax, 0.01f, 1.0f, "%.3f");
				if (ImGui::TreeNode("Distance Metric")) {
					ImGui::RadioButton("Distance From Plane", &distFunc, 0);
					ImGui::SameLine();
					help_marker("Estimate radius function as distance from a plane!");
					if (distFunc == 0) {
						ImGui::SetNextItemWidth(200);
						ImGui::InputFloat3("Normal", distPlaneNormal);
						ImGui::SetNextItemWidth(200);
						ImGui::InputFloat3("Center", distPlaneCenter);
					};
					ImGui::RadioButton("Distance From Mesh Face", &distFunc, 1);
					ImGui::SameLine();
					help_marker("Estimate radius function as distance from the container face!");

					ImGui::RadioButton("Distance From Point", &distFunc, 2);
					if (distFunc == 2) {
						ImGui::SetNextItemWidth(200);
						ImGui::InputFloat3("Point", distPoint);
					}
					ImGui::SameLine();
					help_marker("Estimate radius function as distance from a single point!");
					ImGui::TreePop();
				}
				if (ImGui::TreeNode("Distance - Radius Function")) {
					ImGui::RadioButton("Linear", &radiusFunc, 0);
					ImGui::RadioButton("Quadratic", &radiusFunc, 1);
					ImGui::TreePop();
				}
			}
		}
	}
	ImGui::End();
};

void myGUI::_render_volume_optimization() {

	if (ImGui::Begin("Volume Optimization", NULL)) {

		ImGui::Checkbox("Volume Optimization", &runVolumeOptimization);
		ImGui::SameLine(); help_marker("An optimization approach to enforce pore volumes");
		if (runVolumeOptimization) {
			//ImGui::SeparatorText("Volume Optimization Settings"); 
			ImGui::Separator();
			static int initOpt = 0;

			//Eigen::VectorXd 

			if (ImGui::TreeNode("Weight Initialization")) {

				ImGui::RadioButton("All zero", &initOpt, 0);
				ImGui::RadioButton("Load", &initOpt, 1);
				ImGui::SameLine(); help_marker("Not implemented yet");

				if (initOpt == 0) {
					wInit = Eigen::VectorXd::Zero(seedNr);
				}
				else if (initOpt == 1) {
				};

				ImGui::TreePop();
			}

			if (ImGui::TreeNode("Target Volumes")) {

				ImGui::RadioButton("All equal", &volOption, 0);
				ImGui::RadioButton("Use Radii From Poisson 3D", &volOption, 1);
				ImGui::RadioButton("Load", &volOption, 2);
				ImGui::SameLine(); help_marker("Not implemented yet");

				if (volOption == 0) {

				}

				if (volOption == 2) {
					add_log(LogPriority::ERROR, "Not implemented yet!");
				};

				ImGui::TreePop();
			}
			ImGui::Separator();
		}
	}
	ImGui::End();

};

void myGUI::_render_scaffold_settings() {

	if (ImGui::Begin("Scaffold Settings", NULL)) {

		ImGui::InputText("Model Name", version, IM_ARRAYSIZE(version));

		ImGui::SeparatorText("Scaffold Mesh Generator");
		ImGui::SetNextItemWidth(200);
		ImGui::InputInt("Regularization Steps", &regSteps, 1, 1000);
		ImGui::SameLine(); help_marker("More regularization steps lead to a more regular voronoi grid");

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Thickness", &thickness, 0.1f, 1.0f, "%.3f");
		ImGui::SameLine(); help_marker("Thickness of Scaffold");

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Connectivity Threshold", &connectivityThreshold, 0.01f, 1.0f, "%.3f");
		ImGui::SameLine(); help_marker("0 full face open - 1 maximum pullback");

		connectivityThreshold = std::clamp(connectivityThreshold, 0.001f, 0.999f);

		//ImGui::InputDouble("Hole Scale Factor", &scaleFactor, 0.1, 0.99, "%.3f");
		//ImGui::SameLine(); help_marker("To create holes to each face we estimate the maximum inscribed circle, this factor scales its radius. Default value is 0.5");
	}
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

		// display here the settings for a box set
		static float tempMinX{ 0.0f };
		static float tempMaxX{ 10.0f };
		static float tempMinY{ 0.0f };
		static float tempMaxY{ 10.0f };
		static float tempMinZ{ 0.0f };
		static float tempMaxZ{ 10.0f };

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("x Min (mm)", &tempMinX, 0.1f, 100.0f);
		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("x Max (mm)", &tempMaxX, 0.1f, 100.0f);
		ImGui::SetItemTooltip("Dimension of Scaffold Along X (mm)");

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("y Min (mm)", &tempMinY, 0.1f, 100.0f);
		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("y Max (mm)", &tempMaxY, 0.1f, 100.0f);
		ImGui::SetItemTooltip("Dimension of Scaffold Along Y (mm)");

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("z Min (mm)", &tempMinZ, 0.1f, 100.0f);
		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("z Max (mm)", &tempMaxZ, 0.1f, 100.0f);
		ImGui::SetItemTooltip("Dimension of Scaffold Along Z (mm)");

		ImGui::Separator();
		ImGui::SameLine();

		if (ImGui::Button("Create")) {

			std::unique_ptr<BoxContainer> container = std::make_unique<BoxContainer>(
				tempMinX, tempMaxX, tempMinY, tempMaxY, tempMinZ, tempMaxZ
			);
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

			add_log(LogPriority::SUCCESS, "Created box container!");

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
		// we need the cylinder center, axis, radius

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

			selectedPanelObj.ptr = &containers.back();
			selectedPanelObj.type = ObjectType::CylinderContainerType;
			buffer[0] = '\0';

			add_log(LogPriority::SUCCESS, "Created cylindrical container!");

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

void myGUI::_render_random_seed_creator(const char* popupName, bool& showPopup) {

	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}

	// always centered
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));


	static IContainer* selectedCon = nullptr;
	static int tempSeedNr = { 100 };

	if (ImGui::BeginPopupModal("Random Seed Creator", NULL, ImGuiWindowFlags_AlwaysAutoResize))
	{

		ImGui::SeparatorText("Select Container");
		ImGui::BeginChild("Containers", ImVec2(ImGui::GetContentRegionAvail().x, 80), ImGuiChildFlags_Borders);

		for (const auto& md : containers) {

			IContainer* con = dynamic_cast<IContainer*>(md.get());

			if (con) {
				ImGui::PushID(con);
				bool isSelected = (selectedCon && selectedCon == con);

				if (ImGui::Selectable(md->name.c_str(), isSelected)) {
					selectedCon = md.get();
				};
				ImGui::PopID();
			}
		}
		ImGui::EndChild();

		ImGui::SeparatorText("Parameters");

		static char buffer[256] = "";
		ImGui::InputText("Name", buffer, sizeof(buffer));

		ImGui::InputInt("Number of seeds", &tempSeedNr);
		
		ImGui::Separator();
		ImGui::NewLine();

		if (ImGui::Button("Generate")) {

			if (selectedCon) {
				// create the generator
				std::unique_ptr<Random> rnd = std::make_unique<Random>(tempSeedNr);
				std::string name = std::string(buffer);

				if (!name.empty()) {
					rnd->name = name;
				}
				else {
					rnd->name = "Generator" + std::to_string(seedGenerators.size());
				}
				
				//ContainerAdapter adapter = { *selectedCon, xDim, yDim, zDim };
				//rnd->run(adapter, seeds);
				rnd->run(*selectedCon);
				rnd->container = selectedCon;
				rnd->type = ObjectType::RandomGeneratorType;

				size_t nr = rnd->get_seeds().size();

				// push to the list
				seedGenerators.push_back(std::move(rnd));
				selectedPanelObj.ptr = seedGenerators.back().get();
				selectedPanelObj.type = ObjectType::RandomGeneratorType;
				buffer[0] = '\0';

				add_log(
					LogPriority::SUCCESS,
					std::to_string(nr) + " seeds created randomly inside " + selectedCon->name + "!");

				_update_cameras(*selectedCon);
				
				selectedCon = nullptr;
				showPopup = false;
			}
		}

		ImGui::SameLine();

		if (ImGui::Button("Cancel")) {
			buffer[0] = '\0';
			selectedCon = nullptr;
			showPopup = false;
		};
		ImGui::EndPopup();
	}
};

void myGUI::_render_random_seed_generator_properties() {

	// get selected object
	Random* random = static_cast<Random*>(selectedPanelObj.ptr);
	
	// display the available containers and set the container to the internal pointer
	ImGui::BeginChild("Containers", ImVec2(ImGui::GetContentRegionAvail().x, 80), ImGuiChildFlags_Borders);

	static IContainer* selectedContainer = random->container;

	for (const auto& md : containers) {

		IContainer* con = static_cast<IContainer*>(md.get());

		if (con) {
			if (selectedContainer == nullptr) {
				selectedContainer = md.get();
			}

			bool isSelected = (selectedContainer && selectedContainer == con);

			if (ImGui::Selectable(md->name.c_str(), isSelected)) {
				selectedContainer = md.get();
			};

			if (isSelected) ImGui::SetItemDefaultFocus();
		}
	}

	ImGui::EndChild();
	
	random->render_gui();

	if (ImGui::Button("Update")) {
		if (selectedContainer) {
			//ContainerAdapter adapter = { *selectedContainer, xDim, yDim, zDim };
			random->run(*selectedContainer);
			random->container = selectedContainer;
			add_log(LogPriority::SUCCESS, "Seeds Updated");
		}
	}
};


void myGUI::_render_uniform_seed_creator(const char* popupName, bool& showPopup) {

	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}
	// always centered
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

	if (ImGui::BeginPopupModal("Uniform Seed Creator", NULL, ImGuiWindowFlags_AlwaysAutoResize)){
		ImGui::SameLine();

		static IContainer* selectedCon = nullptr;
		static float tempRadius = { 1.0f };

		ImGui::SeparatorText("Select Container");

		ImGui::BeginChild("Containers", ImVec2(ImGui::GetContentRegionAvail().x, 80), ImGuiChildFlags_Borders);

		for (const auto& md : containers) {

			IContainer* con = static_cast<IContainer*>(md.get());

			if (con) {
				bool isSelected = (selectedCon && selectedCon == con);

				if (ImGui::Selectable(md->name.c_str(), isSelected)) {
					selectedCon = md.get();
				};
			}
		}
		ImGui::EndChild();

		ImGui::SeparatorText("Parameters");

		static char buffer[256] = "";
		ImGui::InputText("Name", buffer, sizeof(buffer));

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Radius", &tempRadius);

		ImGui::Separator();
		ImGui::NewLine();
		if (ImGui::Button("Generate")) {

			if (selectedCon) {

				//ContainerAdapter adapter = { *selectedCon, xDim, yDim, zDim };

				// from the container get the center
				//Bounds bnds = selectedCon->compute_bounds();

				//std::array<double, 3> center = bnds.center;

				std::unique_ptr<Poisson3D> rnd = std::make_unique<Poisson3D>(
					tempRadius, tempRadius, 30);
				std::string name = std::string(buffer);


				if (!name.empty()) {
					rnd->name = name;
				}
				else {
					rnd->name = "Generator" + std::to_string(seedGenerators.size());
				}
				rnd->type = ObjectType::UniformGeneratorType;

				rnd->run(*selectedCon);

				size_t nr = rnd->get_seeds().size();

				// push to the list
				seedGenerators.push_back(std::move(rnd));
				selectedPanelObj.ptr = seedGenerators.back().get();
				selectedPanelObj.type = ObjectType::UniformGeneratorType;
				buffer[0] = '\0';

				add_log(LogPriority::SUCCESS, std::to_string(nr) + " uniform seeds created inside " + selectedCon->name + "!");

				_update_cameras(*selectedCon);
				selectedCon = nullptr;

				showPopup = false;
			}

			buffer[0] = '\0';
			tempRadius = 1.0f;
			showPopup = false;
		}
		ImGui::SameLine();
		if (ImGui::Button("Cancel")) {
			buffer[0] = '\0';
			tempRadius = 1.0f;
			selectedCon = nullptr;
			showPopup = false;
		};
		ImGui::EndPopup();
	}
};

void myGUI::_render_uniform_seed_generator_properties() {

	// get selected object
	Poisson3D* poisson = static_cast<Poisson3D*>(selectedPanelObj.ptr);

	// display the available containers and set the container to the internal pointer
	ImGui::BeginChild("Containers", ImVec2(ImGui::GetContentRegionAvail().x, 100), ImGuiChildFlags_Borders);

	static IContainer* selectedContainer = poisson->container;

	for (const auto& md : containers) {

		IContainer* con = static_cast<IContainer*>(md.get());

		if (con) {
			if (selectedContainer == nullptr) {
				selectedContainer = md.get();
			}

			bool isSelected = (selectedContainer && selectedContainer == con);

			if (ImGui::Selectable(md->name.c_str(), isSelected)) {
				selectedContainer = md.get();
			};

			if (isSelected) ImGui::SetItemDefaultFocus();
		}
	}

	ImGui::EndChild();

	poisson->render_gui();

	if (ImGui::Button("Update")) {
		if (selectedContainer) {
			//ContainerAdapter adapter = { *selectedContainer, xDim, yDim, zDim };
			poisson->run(*selectedContainer);
			poisson->container = selectedContainer;
			add_log(LogPriority::SUCCESS, "Seeds Updated");
		}
	}
};

void myGUI::_render_varied_seed_creator(const char* popupName, bool& showPopup) {

	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}
	// always centered
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

	static IContainer* selectedCon = nullptr;
	static float minRadius = { 0.8f };
	static float maxRadius = { 1.2f };
	static int distanceFunc = 0;
	static int radiusFunc = 0;
	static Vec3 normal{ 0.0f, 0.0f, 1.0f };
	static Vec3 planeCenter{ 0.0f, 0.0f, 1.0f };
	static Vec3 point{ 0.0f, 0.0f, 0.0f };

	if (ImGui::BeginPopupModal("Varied Seed Creator", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
		ImGui::SameLine();

		static char buffer[256] = "";
		ImGui::InputText("Name", buffer, sizeof(buffer));

		ImGui::SeparatorText("Select Container");

		ImGui::BeginChild("Containers", ImVec2(ImGui::GetContentRegionAvail().x, 80), ImGuiChildFlags_Borders);

		for (const auto& md : containers) {

			IContainer* con = dynamic_cast<IContainer*>(md.get());

			if (con) {
				bool isSelected = (selectedCon && selectedCon == con);

				if (ImGui::Selectable(md->name.c_str(), isSelected)) {
					selectedCon = md.get();
				};
			}
		}
		ImGui::EndChild();
		ImGui::SeparatorText("Select Distance Function");
		ImGui::RadioButton("Distance From Plane", &distanceFunc, 0);
		if (distanceFunc == 0) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat3("Normal", normal);
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat3("Center", planeCenter);
		};
		ImGui::RadioButton("Distance From Point", &distanceFunc, 1);
		if (distanceFunc == 1) {
			ImGui::SetNextItemWidth(200);
			ImGui::InputFloat3("Point", point);
		}
		ImGui::RadioButton("Distance From Container", &distanceFunc, 2);

		ImGui::SeparatorText("Select Radius Function");
		ImGui::RadioButton("Linear", &radiusFunc, 0);
		ImGui::RadioButton("Quadratic", &radiusFunc, 1);
		ImGui::RadioButton("Constant", &radiusFunc, 2);

		ImGui::SeparatorText("Parameters");

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Minimum Radius", &minRadius);
		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Maximum Radius", &maxRadius);

		ImGui::Separator();
		ImGui::NewLine();
		if (ImGui::Button("Generate")) {

			if (selectedCon) {

				RunConfig cfg;

				switch (radiusFunc) {
					// linear radius function
					case 0: {
						cfg.rad = std::make_shared<LinearFunction>();
						break;
					}
					case 1: {
						//cfg.rad = std::make_shared<QuadraticFunction>();
						break;
					}
					case 2: {
						cfg.rad = std::make_shared<ConstantRadiusFunction>();
						break;
					}
				}

				switch (distanceFunc) {
					// distance from plane
					case 0: {
						//std::array<double, 3> center = { planeCenter.x, planeCenter.y, planeCenter.z };
						//std::array<double, 3> norm = { normal.x, normal.y, normal.z};
						cfg.dist = std::make_shared<PlaneSDF>(planeCenter, normal);
						break;
					}
						  // distance from point
					case 1: {
						cfg.dist = std::make_shared<PointSDF>(point);
						break;
					}
					// distance from container surface
					case 2: {
						cfg.dist = selectedCon->get_distance_estimator();
						break;
					}
				}

				// from the container get the center
				Bounds bnds = selectedCon->compute_bounds();

				std::unique_ptr<Poisson3D> rnd = std::make_unique<Poisson3D>(
					minRadius, maxRadius, 30);
				std::string name = std::string(buffer);

				if (!name.empty()) {
					rnd->name = name;
				}
				else {
					rnd->name = "Generator" + std::to_string(seedGenerators.size() + 1);
				}
				rnd->type = ObjectType::VariedGeneratorType;
				rnd->distIdx = distanceFunc;
				rnd->radiusIdx = radiusFunc;
				rnd->run(*selectedCon, cfg);

				size_t nr = (int)rnd->get_seeds().size();

				// push to the list
				seedGenerators.push_back(std::move(rnd));
				selectedPanelObj.ptr = seedGenerators.back().get();
				selectedPanelObj.type = ObjectType::VariedGeneratorType;
				buffer[0] = '\0';

				add_log(
					LogPriority::SUCCESS,
					std::to_string(nr) + " varied seeds created inside " + selectedCon->name + "!");

				_update_cameras(*selectedCon);
				selectedCon = nullptr;

				showPopup = false;
			}
		
		}
		ImGui::SameLine();
		if (ImGui::Button("Cancel")) {
			buffer[0] = '\0';
			selectedCon = nullptr;
			showPopup = false;
		};
		ImGui::EndPopup();
	}
};

void myGUI::_render_varied_seed_generator_properties() {
	// get selected object
	Poisson3D* poisson = static_cast<Poisson3D*>(selectedPanelObj.ptr);

	// display the available containers and set the container to the internal pointer
	ImGui::BeginChild("Containers", ImVec2(ImGui::GetContentRegionAvail().x, 100), ImGuiChildFlags_Borders);

	static IContainer* selectedContainer = poisson->container;
	static float minRadius = poisson->get_min_radius();
	static float maxRadius = poisson->get_max_radius();
	static int distanceFunc = poisson->distIdx;
	static int radiusFunc = poisson->radiusIdx;

	static Vec3 normal{ 0.0f, 0.0f, 1.0f };
	static Vec3 planeCenter{ 0.0f, 0.0f, 1.0f };
	static Vec3 point{ 0.0f, 0.0f, 0.0f };

	RunConfig cfg = poisson->get_config();

	if (cfg.dist) {
		if (distanceFunc == 0){
			const PlaneSDF* pd = dynamic_cast<const PlaneSDF*>(cfg.dist.get());
			if (pd) {
				normal = pd->get_normal();
				planeCenter = pd->get_point();
			}
		}
		else if (distanceFunc == 1) {
			const PointSDF* pd = dynamic_cast<const PointSDF*>(cfg.dist.get());
			if (pd) {
				point = pd->get_point();
			}
		}
	}

	for (const auto& md : containers) {

		IContainer* con = static_cast<IContainer*>(md.get());

		if (con) {
			if (selectedContainer == nullptr) {
				selectedContainer = md.get();
			}

			bool isSelected = (selectedContainer && selectedContainer == con);

			if (ImGui::Selectable(md->name.c_str(), isSelected)) {
				selectedContainer = md.get();
			};

			if (isSelected) ImGui::SetItemDefaultFocus();
		}
	}

	ImGui::EndChild();

	ImGui::SeparatorText("Select Distance Function");
	ImGui::RadioButton("Distance From Plane", &distanceFunc, 0);
	if (distanceFunc == 0) {
		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat3("Normal", normal);
		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat3("Center", planeCenter);
	};
	ImGui::RadioButton("Distance From Point", &distanceFunc, 1);
	if (distanceFunc == 1) {
		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat3("Point", point);
	}
	ImGui::RadioButton("Distance From Container", &distanceFunc, 2);

	ImGui::SeparatorText("Select Radius Function");
	ImGui::RadioButton("Linear", &radiusFunc, 0);
	ImGui::RadioButton("Quadratic", &radiusFunc, 1);
	ImGui::RadioButton("Constant", &radiusFunc, 2);

	ImGui::SeparatorText("Parameters");
	ImGui::InputFloat("Minimum Distance", &minRadius);
	ImGui::InputFloat("Maximum Distance", &maxRadius);
	//poisson->render_gui();

	if (ImGui::Button("Update")) {
		if (selectedContainer) {

			RunConfig cfg;

			switch (radiusFunc) {
				// linear radius function
				case 0: {
					cfg.rad = std::make_shared<LinearFunction>();
					break;
				}
				case 1: {
					//cfg.rad = std::make_shared<QuadraticFunction>();
					break;
				}
				case 2: {
					cfg.rad = std::make_shared<ConstantRadiusFunction>();
					break;
				}
				}

				switch (distanceFunc) {
					// distance from plane
				case 0: {
					cfg.dist = std::make_shared<PlaneSDF>(planeCenter, normal);
					break;
				}
					  // distance from container surface
				case 2: {
					cfg.dist = selectedContainer->get_distance_estimator();
					break;
				}
					  // distance from point
				case 1: {

					cfg.dist = std::make_shared<PointSDF>(point);
					break;
				}
			}

			// from the container get the center
			Bounds bnds = selectedContainer->compute_bounds();

			poisson->set_min_radius(minRadius);
			poisson->set_max_radius(maxRadius);
			poisson->run(*selectedContainer, cfg);
			add_log(LogPriority::SUCCESS, "Seeds Updated");
		}
	}
};

void myGUI::_render_scaffold_creator(const char* popupName, bool& showPopup) {

	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}

	// always centered
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

	if (ImGui::BeginPopupModal("Scaffold Creator", NULL, ImGuiWindowFlags_AlwaysAutoResize))
	{
		static IContainer* selectedCon = nullptr;
		static InterfaceSeedGenerator* selectedGen = nullptr;
		static char buffer[256] = "";
		static float tempThickness = { 0.3f };
		static float tempOpeness = {0.5f};
		static int foam = 0;

		ImGui::InputText("Name", buffer, sizeof(buffer));

		ImGui::SeparatorText("Parameters");
		ImGui::InputFloat("Thickness", &tempThickness, 0.001f, 1.0f);
		ImGui::SliderFloat("Openess", &tempOpeness, 0.0f, 1.0f, "%.3f");
		
		ImGui::RadioButton("Porous", &foam, 0);
		ImGui::RadioButton("Foam", &foam, 1);

		// here we should get the seeds from the corresponding creator
		ImGui::SeparatorText("Select Container and generator");

		// timer for flashing
		static float warningFlashTimer1 = 0.0f;

		if (warningFlashTimer1 > 0.0f) {
			warningFlashTimer1 -= ImGui::GetIO().DeltaTime;
		}

		bool isFlashing1 = (warningFlashTimer1 > 0.0f);
		if (isFlashing1) {
			float pulseAlpha = (float)(std::sin(ImGui::GetTime() * 15.0f) * 0.5f + 0.5f);
			ImVec4 flashColor = ImVec4(1.0f, 0.0f, 0.0f, pulseAlpha);
			ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 5.0f);
			ImGui::PushStyleColor(ImGuiCol_Border, flashColor);
		}

		ImGui::BeginChild("Containers", ImVec2(ImGui::GetContentRegionAvail().x * 0.5f, 80), ImGuiChildFlags_Borders);

		for (const auto& md : containers) {

			IContainer* con = dynamic_cast<IContainer*>(md.get());

			if (con) {
				bool isSelected = (selectedCon && selectedCon == con);

				if (ImGui::Selectable(md->name.c_str(), isSelected)) {
					selectedCon = md.get();
				};
			}
		}

		if (isFlashing1) {
			ImGui::PopStyleColor(); // Pop the red border color
			ImGui::PopStyleVar();   // Pop the thick border size
		}

		ImGui::EndChild();

		ImGui::SameLine();

		static float warningFlashTimer2 = 0.0f;

		if (warningFlashTimer2 > 0.0f) {
			warningFlashTimer2 -= ImGui::GetIO().DeltaTime;
		}
		bool isFlashing2 = (warningFlashTimer2 > 0.0f);
		if (isFlashing2) {
			float pulseAlpha = (float)(std::sin(ImGui::GetTime() * 15.0f) * 0.5f + 0.5f);
			ImVec4 flashColor = ImVec4(1.0f, 0.0f, 0.0f, pulseAlpha);
			ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 5.0f);
			ImGui::PushStyleColor(ImGuiCol_Border, flashColor);
		}

		ImGui::BeginChild("Generators", ImVec2(0.0, 80), ImGuiChildFlags_Borders);

		for (const auto& md : seedGenerators) {

			InterfaceSeedGenerator* gen = dynamic_cast<InterfaceSeedGenerator*>(md.get());

			if (gen) {
				bool isSelected = (selectedGen && selectedGen == gen);

				if (ImGui::Selectable(md->name.c_str(), isSelected)) {
					selectedGen = md.get();
				};
			}
		}
		
		if (isFlashing2) {
			ImGui::PopStyleColor(); // Pop the red border color
			ImGui::PopStyleVar();   // Pop the thick border size
		}

		ImGui::EndChild();
		ImGui::Separator();
		ImGui::NewLine();

		if (ImGui::Button("Generate")) {

			if (!selectedCon) {
				warningFlashTimer1 = 1.5f;
			}

			if (!selectedGen) {
				warningFlashTimer2 = 1.5f;
			}

			else if (selectedCon && selectedGen) {

				auto start_time = std::chrono::steady_clock::now();
				
				std::string name = std::string(buffer);
				std::vector<Vec3> seeds = selectedGen->get_seeds();
				Bounds bds = selectedCon->compute_bounds();

				std::array<float, 6> bounds = {
					bds.xMin,
					bds.xMax,
					bds.yMin,
					bds.yMax,
					bds.zMin,
					bds.zMax
				};

				//std::unique_ptr<Generator> scaffold = std::make_unique<Generator>(
				//	seeds, bounds, resolution, tempOpeness, tempThickness, foam
				//);

				std::unique_ptr<GeneratorLewiner> scaffold = std::make_unique<GeneratorLewiner>(
					seeds, bounds, resolution, tempOpeness, tempThickness, foam
				);

				// estimate the scalar field
				scaffold->compute_scalar_field(*selectedCon);

				//scaffold->populate_grids(*selectedCon);

				scaffold->container = selectedCon;
				scaffold->generator = selectedGen;

				if (foam == 1) {
					std::cout << " create foam " << std::endl;
					scaffold->foam = true;
				}

				if (!name.empty()) {
					scaffold->name = name;
				}
				else {
					scaffold->name = "Scaffold" + std::to_string(scaffolds.size() + 1);
				}

				scaffold->marching_cubes();

				scaffold->estimate_metrics(*selectedCon);

				//// push to the scaffold list
				scaffolds.push_back(std::move(scaffold));

				// set it as the selected object
				selectedSceneObj = scaffolds.back().get();
				selectedPanelObj.ptr = scaffolds.back().get();
				selectedPanelObj.type = ObjectType::ScaffoldType;

				auto end_time = std::chrono::steady_clock::now();

				// Calculate the duration in milliseconds
				auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
				
				std::ostringstream oss;
				oss	<< std::fixed << std::setprecision(3) // Set precision to 3 decimal places
					<< duration_ms.count() / 1000.0   // Convert ms to seconds
					<< " seconds!";

				add_log(
					LogPriority::SUCCESS,
					"Created scaffold succesffully using container " + selectedCon->name +
					" and generator " + selectedGen->name + " in" + oss.str());

				// restore ptrs
				selectedCon = nullptr;
				selectedGen = nullptr;
				buffer[0] = '\0';
				showPopup = false;
			}
		}

		ImGui::SameLine();

		if (ImGui::Button("Cancel")) {
			buffer[0] = '\0';
			selectedCon = nullptr;
			selectedGen = nullptr;
			showPopup = false;
		};
		ImGui::EndPopup();
	}
};

void myGUI::_render_scaffold_properties() {

	GeneratorLewiner* gen = static_cast<GeneratorLewiner*>(selectedPanelObj.ptr);

	ImGui::ColorEdit4("Color", (float*)&gen->color);

	static IContainer* selectedContainer = static_cast<IContainer*>(gen->container);

	ImGui::BeginChild("Containers", ImVec2(ImGui::GetContentRegionAvail().x * 0.5f, 100), ImGuiChildFlags_Borders);

	render_selectable_list(containers, selectedContainer, "Containers");

	ImGui::EndChild();

	static InterfaceSeedGenerator* selectedGenerator = static_cast<InterfaceSeedGenerator*>(gen->generator);
	ImGui::SameLine();

	ImGui::BeginChild("Generator List", ImVec2(0, 100), ImGuiChildFlags_Borders);
	
	render_selectable_list(seedGenerators, selectedGenerator, "Generators");

	ImGui::EndChild();

	//ImGui::SeparatorText("Parameters");
	
	gen->render_properties();

	//ImGui::InputFloat("Thickness", &tempThickness, 0.001f, 1.0f);
	//ImGui::SliderFloat("Openess", &tempOpeness, 0.0f, 1.0f, "%.3f");
	
	ImGui::SeparatorText("Metrics");

	//static Metrics metrics = gen->get_metrics();

	gen->render_metrics();

	if (updateScaffold) {

		if (selectedContainer && selectedGenerator) {
			std::vector<Vec3> seeds = selectedGenerator->get_seeds();
			Bounds bds = selectedContainer->compute_bounds();

			std::array<float, 6> bounds = {
				bds.xMin,
				bds.xMax,
				bds.yMin,
				bds.yMax,
				bds.zMin,
				bds.zMax
			};

			gen->set_bounds(bounds);
			gen->set_seeds(seeds);
			gen->container = selectedContainer;
			gen->generator = selectedGenerator;

			gen->set_resolution(resolution);
			gen->compute_scalar_field(*selectedContainer);
			gen->marching_cubes();
			gen->estimate_metrics(*selectedContainer);

			add_log(LogPriority::SUCCESS, "Updated Scaffold Successfully");

			// reset
			updateScaffold = false;
		}
	}


	//if (ImGui::Button("Update")) {	}

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

	// --- RESIZE LOGIC ---
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

	//// --- SCALABLE BOTTOM-RIGHT OVERLAY ---
	//if (scaffoldModel && showScaffold) {

	//	// 1. Configuration
	//	float padding = 15.0f;           // Distance from the edge of the viewport
	//	float widthPct = 0.15f;          // Width is 25% of viewport
	//	float heightPct = 0.20f;         // Height is 28% of viewport

	//	// 2. Calculate Anchor Point (Bottom-Right of the visualizer area)
	//	ImVec2 bottomRightAnchor;
	//	bottomRightAnchor.x = pos.x + availableSize.x - padding;
	//	bottomRightAnchor.y = pos.y + availableSize.y - padding;

	//	// 3. Calculate Size
	//	ImVec2 overlaySize = ImVec2(availableSize.x * widthPct, availableSize.y * heightPct);

	//	// 4. Calculate Font Scale (Based on height, reference 1000px)
	//	float fontScale = availableSize.y / 1000.0f;
	//	if (fontScale < 0.6f) fontScale = 0.6f; // Minimum readable size
	//	if (fontScale > 2.0f) fontScale = 2.0f; // Max size cap

	//	// 5. Set Window Properties
	//	ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoDecoration |
	//		ImGuiWindowFlags_NoDocking |
	//		ImGuiWindowFlags_NoSavedSettings |
	//		ImGuiWindowFlags_NoFocusOnAppearing |
	//		ImGuiWindowFlags_NoNav |
	//		ImGuiWindowFlags_NoMove |
	//		ImGuiWindowFlags_NoInputs; // Optional: Click-through

	//	// PIVOT MAGIC: (1.0f, 1.0f) means the coordinates provided are the Bottom-Right corner
	//	ImGui::SetNextWindowPos(bottomRightAnchor, ImGuiCond_Always, ImVec2(1.0f, 1.0f));
	//	ImGui::SetNextWindowSize(overlaySize, ImGuiCond_Always);
	//	ImGui::SetNextWindowBgAlpha(0.45f); // Slightly darker for readability

	//	if (ImGui::Begin("Mesh Properties", nullptr, window_flags))
	//	{
	//		ImGui::SetWindowFontScale(fontScale);

	//		ImGui::Text("Mesh Details");
	//		ImGui::Separator();

	//		// Use WrapPos to prevent text spilling if window gets narrow
	//		ImGui::PushTextWrapPos(ImGui::GetContentRegionAvail().x);

	//		ImGui::Text("Name: %s", scaffoldFileName.c_str());
	//		ImGui::Text("Vertices: %d", vertexNr);
	//		ImGui::Text("Faces: %d", faceNr);
	//		ImGui::Text("Vol: %.3f", scaffoldVolume);
	//		ImGui::Text("Porosity: %.3f", scaffoldPorosity);
	//		ImGui::Text("Conn: %.3f", scaffoldConnectivity);

	//		ImGui::PopTextWrapPos();
	//		ImGui::SetWindowFontScale(1.0f); // Reset scale
	//		ImGui::End();
	//	}
	//}

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

		ImGui::InputFloat("Voxel Size (mm)", &voxelSize, 0.01f, 10.0f);

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

void myGUI::_local_thickness_measure(const char* popupName, bool& showPopup, bool flag) {
	
	// grab the active generatior
	if (selectedPanelObj.type != ObjectType::ScaffoldType) {

		add_log(LogPriority::ERROR, "Select a scaffold in the panel.");
		showPopup = false;
		return;
	}
	Generator* gen = static_cast<Generator*>(selectedPanelObj.ptr);

	if (!gen) {
		showPopup = false;
		return;
	}

	// fire up a modal
	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}
	static float tempVoxelSize = 0.05f;
	static std::array<float, 2> xDims = { 0.0f, 5.0f};
	static std::array<float, 2> yDims = { 0.0f, 5.0f };
	static std::array<float, 2> zDims = { 0.0f, 5.0f };
	if (ImGui::BeginPopupModal(popupName, nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
		
		ImGui::SeparatorText("Parameters");

		ImGui::InputFloat("Voxel Size (mm)", &tempVoxelSize, 0.001f, 1.0f);
		ImGui::InputFloat2("X Dimensions (mm)", xDims);
		ImGui::InputFloat2("Y Dimensions (mm)", yDims);
		ImGui::InputFloat2("Z Dimensions (mm)", zDims);

		ImGui::NewLine();
		if (ImGui::Button("Estimate")) {

			voxelBounds = {
				xDims[0], xDims[1],
				yDims[0], yDims[1],
				zDims[0], zDims[1]
			};

			gen->estimate_local_thickness(voxelSize, voxelBounds, flag);

			showPopup = false;
		}
		
		ImGui::SameLine();
		if (ImGui::Button("Cancel")) {
			showPopup = false;
			ImGui::CloseCurrentPopup();
		}

		ImGui::EndPopup();
	}
};

void myGUI::_action_estimate_tortuosity() {

	// get the selected scaffold
	if (selectedPanelObj.type != ObjectType::ScaffoldType) {
		add_log(LogPriority::ERROR, "Select a scaffold from the left panel.");
		return;
	}

	// get the scaffold
	Generator* scaffold = static_cast<Generator*>(selectedPanelObj.ptr);

	if (scaffold) {
		scaffold->estimate_tortuosity();
	}
};

void myGUI::_action_estimate_pore_network() {

	// get the selected scaffold
	if (selectedPanelObj.type != ObjectType::ScaffoldType) {
		add_log(LogPriority::ERROR, "Select a scaffold from the left panel.");
		return;
	}

	// get the scaffold
	Generator* scaffold = static_cast<Generator*>(selectedPanelObj.ptr);

	if (scaffold) {
		scaffold->estimate_connectivity_network();
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

		ImGui::InputFloat("Translate X", &tempX, 0.001f, 1000.0f);
		ImGui::InputFloat("Translate Y", &tempY, 0.001f, 1000.0f);
		ImGui::InputFloat("Translate Z", &tempZ, 0.001f, 1000.0f);

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

void myGUI::_draw_selected_box() {

	if (!selectedSceneObj) return;

	boxShader.use();
	glDepthFunc(GL_LEQUAL);
	glLineWidth(1.0f);

	auto md = static_cast<GeneratorLewiner*>(selectedSceneObj);

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
