#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <iostream>
#include <fstream> 
#include <json.hpp>
#include <cmath>
#include <functional>
#include <filesystem>
#include <vtkBox.h>
#include <vtkCylinder.h>
#include <vtkMassProperties.h>
#include <vtkSTLWriter.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkDoubleArray.h>
#include "guiApp.h"
#include <ImGuiFileDialog/ImGuiFileDialog.h>
#include "buildScaffold.h"
#include "Utils/Utils.h"
#include "ScaffoldGenerator/ScaffoldGenerator.h"
#include "SeedGenerator/DistanceCalculator.h"
#include "Model.h"
#include "SeedGenerator/Random.h"
#include "SeedGenerator/Poisson3D.h"
#include "Logger/Logger.h"


myGUI::myGUI(int width, int height) : width(width), height(height) {
	
	_init_opengl();

	glsl_version = _get_glsl_version().c_str();

	_init_imgui();
};

void myGUI::_init_opengl() {

	std::cout << "starting opengl" << std::endl;

	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// create opengl window
	window = glfwCreateWindow(width, height, "Scaffold Builder", NULL, NULL);
	if (window == NULL)
	{
		glfwTerminate();
		throw std::runtime_error(std::string(std::string("Failed to open GLFW window.") +
			" If you have an Intel GPU, they are not 3.3 compatible." +
			"Try the 2.1 version.\n"));
	}

	// create context
	glfwMakeContextCurrent(window);

	//glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	// glad: load all OpenGL function pointers, these are OS-specific
	// ----------------------------------------------------------------
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cerr << "Failed to initialize GLAD" << std::endl;
		return;
	}

	// Now you can check the OpenGL version
	const GLubyte* version = glGetString(GL_VERSION);
	const GLubyte* glslVersion = glGetString(GL_SHADING_LANGUAGE_VERSION);

	// Print the OpenGL version
	std::cout << "OpenGL Version: " << version << std::endl;
	std::cout << "GLSL Version: " << glslVersion << std::endl;

	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

	// Hide the mouse and enable unlimited movement
	//glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	// Set the mouse at the center of the screen
	glfwPollEvents();
	glfwSetCursorPos(window, width / 2, height / 2);

	// scroll camera
	//glfwSetScrollCallback(window, scroll_callback);

	// Gray background color
	glClearColor(fontColor[0], fontColor[1], fontColor[2], fontColor[3]);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);

	// enable blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	//glEnable(GL_CULL_FACE);
	
	// load shaders
	scaffoldShader = Shader(
		"../src/shaders/mainVShader.vertexshader",
		"../src/shaders/mainFShader.fragmentshader",
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

	//// pass light
	scaffoldShader.use();
	uniManager.setUniform(scaffoldShader, "lightColor", lightColor[0], lightColor[1], lightColor[2]);
	uniManager.setUniform(scaffoldShader, "lightPosWorld", lightPosCamera[0], lightPosCamera[1], lightPosCamera[2]);
	uniManager.setUniform(scaffoldShader, "Ka", Ka);

	normalShader = Shader(
		"../src/shaders/normalVShader.vs",
		"../src/shaders/normalFShader.fs",
		"../src/shaders/normalGShader.gs"
	);

	uniManager.add_uniform(normalShader, "view");
	uniManager.add_uniform(normalShader, "model");
	uniManager.add_uniform(normalShader, "projection");
	uniManager.add_uniform(normalShader, "normalColor");

	edgeShader = Shader(
		"../src/shaders/edgeVShader.vs",
		"../src/shaders/edgeFShader.fs",
		NULL
	);
	uniManager.add_uniform(edgeShader, "view");
	uniManager.add_uniform(edgeShader, "model");
	uniManager.add_uniform(edgeShader, "projection");
	uniManager.add_uniform(edgeShader, "cutPlaneCoeffs");
	uniManager.add_uniform(edgeShader, "cutPlane");

	gridShader = Shader(
		"../src/shaders/gridVShader.vs",
		"../src/shaders/gridFShader.fs",
		"../src/shaders/gridGShader.gs"
	);
	uniManager.add_uniform(gridShader, "view");
	uniManager.add_uniform(gridShader, "model");
	uniManager.add_uniform(gridShader, "projection");
	uniManager.add_uniform(gridShader, "gridColor");

	seedShader = Shader(
		"../src/shaders/seedVShader.vertexshader",
		"../src/shaders/seedFShader.fragmentshader",
		NULL
	);
	uniManager.add_uniform(seedShader, "projection");
	uniManager.add_uniform(seedShader, "view");
	uniManager.add_uniform(seedShader, "seedColor");
	uniManager.add_uniform(seedShader, "seedSize");

	cutShader = Shader(
		"../src/shaders/cutVShader.vs",
		"../src/shaders/cutFShader.fs",
		NULL
	);
	uniManager.add_uniform(cutShader, "projection");
	uniManager.add_uniform(cutShader, "view");
	uniManager.add_uniform(cutShader, "model");
	uniManager.add_uniform(cutShader, "minBounds");
	uniManager.add_uniform(cutShader, "maxBounds");

	bboxShader = Shader(
		"../src/shaders/bboxVShader.vs",
		"../src/shaders/bboxFShader.fs",
		NULL
	);
	uniManager.add_uniform(bboxShader, "view");
	uniManager.add_uniform(bboxShader, "model");
	uniManager.add_uniform(bboxShader, "projection");

	containerShader = Shader(
		"../src/shaders/containerVShader.vs",
		"../src/shaders/containerFShader.fs",
		NULL
	);
	uniManager.add_uniform(containerShader, "projection");
	uniManager.add_uniform(containerShader, "view");
	uniManager.add_uniform(containerShader, "model");
	uniManager.add_uniform(containerShader, "lightColor");
	uniManager.add_uniform(containerShader, "objectColor");
	uniManager.add_uniform(containerShader, "lightPosWorld");
	uniManager.add_uniform(containerShader, "Ka");

	containerShader.use();
	uniManager.setUniform(containerShader, "lightColor", lightColor[0], lightColor[1], lightColor[2]);
	uniManager.setUniform(containerShader, "lightPosWorld", lightPosCamera[0], lightPosCamera[1], lightPosCamera[2]);
	uniManager.setUniform(containerShader, "Ka", Ka);

	frameShader = Shader(
		"../src/shaders/frameVShader.vs",
		"../src/shaders/frameFShader.fs",
		NULL
	);
	frameShader.use();
	uniManager.add_uniform(frameShader, "projection");
	uniManager.add_uniform(frameShader, "view");
	uniManager.add_uniform(frameShader, "model");
	uniManager.add_uniform(frameShader, "outColor");
	uniManager.add_uniform(frameShader, "screenSize");

	std::cout << "Creating arrows: " << std::endl;
	xArrow = std::make_unique<Arrow>();
	yArrow = std::make_unique<Arrow>();
	zArrow = std::make_unique<Arrow>();
	std::cout << "Created arrows: " << std::endl;

	float sphereRadius{ 0.0f };
	if (width > height)
		sphereRadius = height * 0.5f;
	else
		sphereRadius = width * 0.5f;

	defCamera = new defaultCamera(window, 0.3f, glm::vec3(0.0f, 0.0f, 10.0f), cameraTarget, 2.0f);
	//trackCamera = new trackBallCamera(window, split, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0, 0.0, 20.0));
	trackCamera = std::make_unique<TrackBall>(window, sphereRadius, framebuffer.width, framebuffer.height, 0, 0);

	//defCamera.scroll_setup();
	//framebuffer.width = width;
	//framebuffer.height = height;

	// create frame buffer

	std::cout << "init opengl create buffer" << std::endl;
	create_frame_buffer(framebuffer);
};

void myGUI::_init_imgui() {

	std::cout << "Starting imgui" << std::endl;

	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();

	io = ImGui::GetIO(); (void)io;
	io.Fonts->AddFontFromFileTTF("../lib/imgui/misc/fonts/DroidSans.ttf", 20.0f);
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

	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init(glsl_version.c_str());
}

void myGUI::run() {

	grid = std::make_unique<Grid>(Grid::XY);

	while (!glfwWindowShouldClose(window))
	{

		glfwGetFramebufferSize(window, &width, &height);

		glfwPollEvents();

		// Start a new ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		_create_dockspace();

		_render_settings_panel();

		glBindFramebuffer(GL_FRAMEBUFFER, framebuffer.id);
		glViewport(0, 0, framebuffer.width, framebuffer.height);

		glClearColor(fontColor[0], fontColor[1], fontColor[2], fontColor[3]);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		_render_visualizer();

		double xPos, yPos;
		glfwGetCursorPos(window, &xPos, &yPos);

		//std::cout << xPos << " " << yPos << std::endl;

		bool cameraUpdateFlags = !ImGuiFileDialog::Instance()->IsOpened() && !ImGui::IsAnyItemActive() && !ImGui::IsAnyItemHovered();

		//if (cameraUpdateFlags) {
		//	if (cameraOption == defaultOption) {
		//		defCamera->scroll_setup(window);
		//		defCamera->update();
		//		projection = defCamera->projectionMatrix;
		//		view = defCamera->viewMatrix;
		//		model = model;
		//	}
		//	else if (cameraOption == trackOption) {

		//		float sphereRadius{ 0.0f };
		//		sphereRadius = std::min(0.3f * framebuffer.width, static_cast<float>(framebuffer.height)) * 0.5f;
		//		trackCamera->set_radius(sphereRadius);
		//		//trackCamera->set_screen_size(0.3f * framebuffer.width, framebuffer.height, 0);
		//		trackCamera->update();
		//		Vec3 cameraPos = trackCamera->get_position();
		//		projection = trackCamera->get_projection_matrix();
		//		view = trackCamera->get_view_matrix();
		//		model = model;
		//	}
		//}

		float sphereRadius;
		sphereRadius = std::min(static_cast<float>(framebuffer.width), static_cast<float>(framebuffer.height)) * 0.5f;
		trackCamera->set_radius(sphereRadius);
		//trackCamera->set_screen_size(0.3f * framebuffer.width, framebuffer.height, 0);
		trackCamera->update();
		Vec3 cameraPos = trackCamera->get_position();
		projection = trackCamera->get_projection_matrix();
		view = trackCamera->get_view_matrix();
		model = model;

		scaffoldShader.use();

		if (showCutPlane) {

			float offset = cutPlane->offset;
			glm::vec3 normal = cutPlane->normal;

			float d = -glm::dot(normal, cutPlane->center + normal * offset);

			planeCoeffs = glm::vec4(normal, d);

			uniManager.setUniform(scaffoldShader, "cutPlane", 1);
			uniManager.setUniform(scaffoldShader, "cutPlaneCoeffs", planeCoeffs);
		}
		else {
			uniManager.setUniform(scaffoldShader, "cutPlane", 0);
		}

		if (scaffoldModel){
			uniManager.setUniform(scaffoldShader, "projection", projection);
			uniManager.setUniform(scaffoldShader, "view", view);
			uniManager.setUniform(scaffoldShader, "model", glm::mat4(1.0));
			uniManager.setUniform(scaffoldShader, "objectColor", scaffoldColor[0], scaffoldColor[1], scaffoldColor[2]);
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(1.0f, 1.0f);
			if (showScaffold) {
				scaffoldModel->draw();
			}
			glDisable(GL_POLYGON_OFFSET_FILL);
		}

		if (showNormals && scaffoldModel) {
			normalShader.use();
			uniManager.setUniform(normalShader, "projection", projection);
			uniManager.setUniform(normalShader, "view", view);
			uniManager.setUniform(normalShader, "model", model);
			uniManager.setUniform(normalShader, "normalColor", normalColor[0], normalColor[1], normalColor[2]);
			if (showScaffold) {
				scaffoldModel->draw();
			}
		}

		if (showEdges && scaffoldModel) {
			edgeShader.use();

			if (showCutPlane) {
				uniManager.setUniform(edgeShader, "cutPlane", 1);
				uniManager.setUniform(edgeShader, "cutPlaneCoeffs", planeCoeffs);
			}
			else {
				uniManager.setUniform(edgeShader, "cutPlane", 0);
			}
			glDepthFunc(GL_LEQUAL);
			uniManager.setUniform(edgeShader, "projection", projection);
			uniManager.setUniform(edgeShader, "view", view);
			uniManager.setUniform(edgeShader, "model", model);
			if (showScaffold) {
				scaffoldModel->draw_edges();
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

		if (showCutPlane && scaffoldModel) {
			modelPlane = cutPlane->tMatrix * cutPlane->initMatrix * cutPlane->rotMatrix;
			float xc = (xMax + xMin) * 0.5f;
			float yc = (yMax + yMin) * 0.5f;
			float zc = (zMax + zMin) * 0.5f;
			glm::vec3 minBounds = { xMin, yMin, zMin };
			glm::vec3 maxBounds = { xMax, yMax, zMax };
			cutShader.use();
			uniManager.setUniform(cutShader, "projection", projection);
			uniManager.setUniform(cutShader, "view", view);
			uniManager.setUniform(cutShader, "model", modelPlane);
			uniManager.setUniform(cutShader, "minBounds", xMin, yMin, zMin);
			uniManager.setUniform(cutShader, "maxBounds", xMax, yMax, zMax);
			cutPlane->draw();

			bboxShader.use();
			uniManager.setUniform(bboxShader, "projection", projection);
			uniManager.setUniform(bboxShader, "view", view);
			uniManager.setUniform(bboxShader, "model", model);
			box->draw();
		}

		if (showContainer && containerModel) {
			if (conOption == 0) {
				bboxShader.use();
				uniManager.setUniform(bboxShader, "projection", projection);
				uniManager.setUniform(bboxShader, "view", view);
				uniManager.setUniform(bboxShader, "model", model);
				box->draw();
			}
			else if (conOption == 2) {
				containerShader.use();
				uniManager.setUniform(containerShader, "projection", projection);
				uniManager.setUniform(containerShader, "view", view);
				uniManager.setUniform(containerShader, "model", model);
				uniManager.setUniform(
					containerShader, "objectColor",
					containerColor[0], containerColor[1], containerColor[2], containerColor[3]);
				containerModel->draw();
			}
		}

		if (showSeeds && seedObj) {
			seedShader.use();
			uniManager.setUniform(seedShader, "projection", projection);
			uniManager.setUniform(seedShader, "view", view);
			uniManager.setUniform(seedShader, "seedSize", seedSize);
			uniManager.setUniform(seedShader, "seedColor", seedColor[0], seedColor[1], seedColor[2]);
			seedObj->draw();
		}

		if (showDistancePlane) {
			std::cout << "Draw distance plane " << std::endl;
			glm::mat4 modelDist = distPlane->tMatrix * distPlane->initMatrix * distPlane->rotMatrix;
			cutShader.use();
			uniManager.setUniform(cutShader, "projection", projection);
			uniManager.setUniform(cutShader, "view", view);
			uniManager.setUniform(cutShader, "model", modelDist);
			uniManager.setUniform(cutShader, "minBounds", 0, 0, 0);
			uniManager.setUniform(cutShader, "maxBounds", xDim, yDim, zDim);
			distPlane->draw();
		}

		if (glfwGetKey(window, GLFW_KEY_DELETE) == GLFW_PRESS) {
			
			containerCenter[0] = 0.0;
			containerCenter[1] = 0.0;
			containerCenter[2] = 0.0;
			_update_cameras();
			if (scaffoldModel) {
				scaffoldModel->clean();
				scaffoldModel.reset();
				//scaffoldReady = false;
				scaffoldCreated = false;
				showContainer = false;
				showSeeds = false;
			}
			if (seedObj) {
				seedObj->deleteObj();
				seedObj.reset();
				showSeeds = false;
			}
			//if (glfwGetKey(window, GLFW_KEY_DELETE) == GLFW_RELEASE) {
			//	add_log(LogPriority::INFO, "Deleting Everything");
			//}
		}

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

void myGUI::_render_settings_panel() {

	//// get window width and height
	//int windowWidth, windowHeight;

	//glfwGetFramebufferSize(window, &windowWidth, &windowHeight);

	//ImGuiIO& io = ImGui::GetIO();
	//ImGuiID dockspaceId = ImGui::GetID("MainDockSpace");

	//ImGuiWindowFlags flags = 0;

	//flags = ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoBringToFrontOnFocus;

	//// first render the main menu bar
	//_render_main_menu_bar();

	//// Set window position to the left and window dimensions based on window size
	//ImGui::SetNextWindowSize(ImVec2(300, height));
	//ImGui::SetNextWindowPos(ImVec2(0, 0));
	//ImGui::Begin("Scaffold Modeler", nullptr);

	////ImGui::DockSpace(dockspaceId, ImVec2(0.0f, 0.0f), ImGuiDockNodeFlags_PassthruCentralNode);

	//float menuWidth = ImGui::GetWindowWidth();
	//split = menuWidth / windowWidth;

	//static ImVec4 currentScaffoldColor = ImVec4(1.0f, 0.5f, 0.5f, 1.0f);
	//static ImVec4 seedColor = ImVec4(1.0f, 0.0f, 0.0f, 1.0f);
	//static float pointSize{ 5.0f };

	// display settings window
	if (showDisplaySettingsWin) {
		_render_display_settings();
	}

	if (showDisplayMeshSettingsWin) {
		_render_mesh_settings();
	}

	if (showCutPlane) {
		_render_cut_tool();
		//cutPlane = new CutPlane();
	}

	//if (showSeeds && !seeds.empty()) {
	//	std::cout << "creating seeds obj" << std::endl;
	//	seedObj = new VisualizeSeeds(seeds);
	//}

	// display file dialog
	// Display the file dialog
	if (ImGuiFileDialog::Instance()->Display("ChooseFileDlgKey")) {
		if (ImGuiFileDialog::Instance()->IsOk()) {
			std::string filePath = ImGuiFileDialog::Instance()->GetFilePathName();
			scaffoldFileName = ImGuiFileDialog::Instance()->GetCurrentFileName();

	//		// Process the file, open, load etc.
			if (loadedMesh == scaffold) {

				//scaffoldMesh = new Mesh(filePath);
				//scaffoldModel = new Model(filePath);
				scaffoldModel = std::make_unique<Model>(filePath);
				xMin = scaffoldModel->xMin;
				xMax = scaffoldModel->xMax;
				yMin = scaffoldModel->yMin;
				yMax = scaffoldModel->yMax;
				zMin = scaffoldModel->zMin;
				zMax = scaffoldModel->zMax;

				//box = new BBox(xMin, xMax, yMin, yMax, zMin, zMax);
				box = std::make_unique<BBox>(xMin, xMax, yMin, yMax, zMin, zMax);

				//scaffoldReady = true;
				scaffoldFilePath = filePath;
				float xc = static_cast<float>((xMax + xMin) * 0.5f);
				float yc = static_cast<float>((yMax + yMin) * 0.5f);
				float zc = static_cast<float>((zMax + zMin) * 0.5f);
				containerCenter[0] = xc;
				containerCenter[1] = yc;
				containerCenter[2] = zc;
				// change camera target and camera position
				_update_cameras();

				scaffoldModel->get_details(faceNr, vertexNr, edgeNr, scaffoldVolume, scaffoldPorosity);

				add_log(LogPriority::INFO, "Loading Scaffold Succesfull");
			}

			else if (loadedMesh == bone) {

				/*	boneMesh = new Mesh(filePath);
					boneReady = true;
					boneFileName = filePath;*/
				std::cout << boneFileName << std::endl;
			}

		}
		// Close the file dialog
		ImGuiFileDialog::Instance()->Close();
	}

	////ImGui::Spacing();

	//ImGui::Separator();

	_render_container_options(conOption);

	_render_seed_generator();

	_render_volume_optimization();

	_render_scaffold_settings();

	if (ImGuiFileDialog::Instance()->Display("Export File")) {
		if (ImGuiFileDialog::Instance()->IsOk()) { // action if OK
			scaffoldFilePath = ImGuiFileDialog::Instance()->GetFilePathName();
			scaffoldFileName = ImGuiFileDialog::Instance()->GetCurrentFileName();
			_export_mesh();
		}
		ImGuiFileDialog::Instance()->Close();
	}

	_render_console();

	//ImGui::End();
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

	std::cout << "OpenGL version: " << glVersion << std::endl;

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

	// Ensure the console box is always scrolled to the bottom when a new message is added
	//if (ImGui::Button("Clear")) {
	//	log.clear();
	//}

	if (ImGui::Begin("Console", NULL)) {
		// Start a scrolling region inside the console
		ImGui::BeginChild("ScrollingConsole", ImVec2(0, 200), true, ImGuiWindowFlags_AlwaysVerticalScrollbar | ImGuiWindowFlags_AlwaysHorizontalScrollbar);
		for (const auto& entry : logger.get_logs()) {
			ImGui::TextUnformatted(entry.c_str());
		}
		ImGui::EndChild();
	}
	ImGui::End();

}

void myGUI::add_log(LogPriority priority, const std::string& message) {
	logger.log(priority, message);
}

void myGUI::_generate_scaffold() {
	
	std::cout << wInit << std::endl;
	std::cout << " -------------------- " << std::endl;
	std::cout << targetVols << std::endl;

	// create scaffold builder object, use a lambda expression for 
	// passing the logging function
	buildScaffold bsc(
		targetVols,
		wInit,
		"settings.json"
		//[this](const std::string& message) {
		//	this->add_log(message);}
	);

	bsc.loop();
		
};

void myGUI::_render_scaffold(){
	
	std::string v{ version };
	//currentMesh = Mesh::Mesh(v + "_model.stl");

	//scaffold.draw();

};

void myGUI::_render_container_options(int& conOption) {

	if (ImGui::Begin("Container Options", NULL)) {

		ImGui::RadioButton("Box Container", &conOption, 0);
		ImGui::SameLine(); help_marker("An orthogonal box container");
		if (conOption == 0) {
			ImGui::SliderFloat("x Dimensions", &xDim, 0, 10);
			ImGui::SameLine(); help_marker("Dimension of Scaffold Along X");

			ImGui::SliderFloat("y Dimensions", &yDim, 0, 10);
			ImGui::SameLine(); help_marker("Dimension of Scaffold Along Y.");

			ImGui::SliderFloat("z Dimensions", &zDim, 0, 10);
			ImGui::SameLine(); help_marker("Dimension of Scaffold Along Z.");
		}

		ImGui::RadioButton("Cylindrical Container", &conOption, 1);
		ImGui::SameLine(); help_marker("A cylindrical container, starting from (0, 0, 0). Select Axis, Radius. The Height depends on the Bounds");

		if (conOption == 1) {
			ImGui::RadioButton("X", &cylinderDir, 0);
			ImGui::SameLine(); ImGui::RadioButton("Y", &cylinderDir, 1);
			ImGui::SameLine(); ImGui::RadioButton("Z", &cylinderDir, 2);

			if (cylinderDir == 0) {
				cylinderAxis[0] = 1.0;
				cylinderAxis[1] = 0.0;
				cylinderAxis[2] = 0.0;
				cylinderHeight = std::abs(xMax - xMin);
			}
			else if (cylinderDir == 1) {
				cylinderAxis[0] = 0.0;
				cylinderAxis[1] = 1.0;
				cylinderAxis[2] = 0.0;
				cylinderHeight = std::abs(yMax - yMin);
			}
			else if (cylinderDir == 2) {
				cylinderAxis[0] = 0.0;
				cylinderAxis[1] = 0.0;
				cylinderAxis[2] = 1.0;
				cylinderHeight = std::abs(zMax - zMin);
			}

			ImGui::InputDouble("Cylinder Radius", &cylinderRadius, 0.1, 0.0, "%.3f");
			ImGui::SameLine(); help_marker("The radius of the cylinder.");
			//ImGui::InputDouble("Cylinder Height", &cylinderHeight, 0.1, 1.0, "%.3f");
			//ImGui::SameLine(); help_marker("The height of the cylinder.");
		}

		ImGui::RadioButton("STL Mesh", &conOption, 2);
		ImGui::SameLine(); help_marker("A custom container geometry. In this case the bounding box of the mesh is used for the voro++ container and the scaffold is built inside the mesh domain");

		if (conOption == 2) {

			if (ImGuiFileDialog::Instance()->Display("Select STL Container")) {
				if (ImGuiFileDialog::Instance()->IsOk()) { // action if OK
					std::string filePathName = ImGuiFileDialog::Instance()->GetFilePathName();
					containerFile = filePathName;
					std::cout << containerFile << std::endl;

					containerMesh = vtkSmartPointer<vtkPolyData>::New();
					vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
					reader->SetFileName(containerFile.c_str());
					reader->Update();
					containerMesh = reader->GetOutput();

					double bounds[6];
					containerMesh->GetBounds(bounds);
					xMin = bounds[0] - 1.0;
					xMax = bounds[1] + 1.0;
					yMin = bounds[2] - 1.0;
					yMax = bounds[3] + 1.0;
					zMin = bounds[4] - 1.0;
					zMax = bounds[5] + 1.0;

					// set also the mesh center as the camera pos
					float xc = (xMax + xMin) * 0.5f;
					float yc = (yMax + yMin) * 0.5f;
					float zc = (zMax + zMin) * 0.5f;
					containerCenter[0] = xc;
					containerCenter[1] = yc;
					containerCenter[2] = zc;

				}
				ImGuiFileDialog::Instance()->Close();
			}
			selectFileButton("Select Geometry File", "../data/", "Select STL Container", ".stl");
			ImGui::InputInt("Neighbors", &neighbors, 1, 100);
			ImGui::InputFloat("Wall Resolution", &wallResolution, 0.1f, 5.0f);
		}
	}
	ImGui::End();
		
};

//void myGUI::framebuffer_size_callback_imp(int width, int height) {
//	defCamera.set_window_size(width, height);
//	trackCamera.set_window_size(width, height);
//};

void myGUI::_update_cameras() {
	cameraTarget = glm::vec3(containerCenter[0], containerCenter[1], containerCenter[2]);
	cameraPos = cameraTarget - 20.0f * glm::vec3(0.0f, 0.0f, 1.0f);
	defCamera->position = cameraPos;
	defCamera->target = cameraTarget;
	//trackCamera->position = cameraPos;
	trackCamera->set_target(cameraTarget.x, cameraTarget.y, cameraTarget.z);
}

void myGUI::_update_bounds_center(int& conOption) {

	if (conOption == 0 ) {
		xMin = 0.0;
		xMax = xDim;
		yMin = 0.0;
		yMax = yDim;
		zMin = 0.0;
		zMax = zDim;

		containerCenter[0] = xDim / 2.0f;
		containerCenter[1] = yDim / 2.0f;
		containerCenter[2] = zDim / 2.0f;
	}

	else if (conOption == 1) {

		// if the cylinder is along the x axis
		if (cylinderDir == 0) {

			xMin = 0.0;
			xMax = xDim;
			yMin = -cylinderRadius;
			yMax = cylinderRadius;
			zMin = -cylinderRadius;
			zMax = cylinderRadius;

			containerCenter[0] = xDim * 0.5;
			containerCenter[1] = 0.0;
			containerCenter[2] = 0.0;
		}
		else if (cylinderDir == 1) {

			xMin = -cylinderRadius;
			xMax = cylinderRadius;
			yMin = 0.0;
			yMax = yDim;
			zMin = -cylinderRadius;
			zMax = cylinderRadius;

			containerCenter[0] = 0.0;
			containerCenter[1] = yDim * 0.5;
			containerCenter[2] = 0.0;
		}
		else if (cylinderDir == 2) {

			xMin = -cylinderRadius;
			xMax = cylinderRadius;
			yMin = -cylinderRadius;
			yMax = cylinderRadius;
			zMin = 0.0;
			zMax = zDim;

			containerCenter[0] = 0.0;
			containerCenter[1] = 0.0;
			containerCenter[2] = zDim * 0.5;
		}
	}

	else if(conOption == 2) {

		double bounds[6];
		containerMesh->GetBounds(bounds);

		xMin = bounds[0];
		xMax = bounds[1];
		yMin = bounds[2];
		yMax = bounds[3];
		zMin = bounds[4];
		zMax = bounds[5];

		containerCenter[0] = static_cast<double>((xMin + xMax) * 0.5f);
		containerCenter[1] = static_cast<double>((yMin + yMax) * 0.5f);
		containerCenter[2] = static_cast<double>((zMin + zMax) * 0.5f);
	}
}

void myGUI::_render_mesh_settings() {

	ImGui::SetNextWindowSize(ImVec2(200, 250), ImGuiCond_FirstUseEver);
	if (ImGui::Begin("Display Settings", &showDisplayMeshSettingsWin)) {

		static int picked = { -99 };
		// Left side - Selectable items
		{
			ImGui::BeginChild("LeftPanel", ImVec2(150, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
			if (ImGui::Selectable("Resolution Settings", picked == 0)) picked = 0;
			ImGui::EndChild();
		}
		ImGui::SameLine();

		// Right side - Color and size controls
		{
			ImGui::BeginGroup();
			{
				ImGui::BeginChild("settings item");

				if (picked == 0) {
					ImGui::Text("Resolution Settings");
					ImGui::Checkbox("Custom", &customResolutionFlag);
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

void myGUI::_export_mesh() {

	add_log(LogPriority::INFO, "Exporting Mesh...");

	std::filesystem::path filePath = scaffoldFilePath;

	if (filePath.extension() == ".stl") {

		vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
		stlWriter->SetFileName(scaffoldFilePath.c_str());
		//stlWriter->SetInputConnection(sincSmoother->GetOutputPort());
		stlWriter->SetInputData(scaffoldPolyData);
		stlWriter->Write();

	}
	else if (filePath.extension() == ".vtk") {

		// create header to store volume porosity etc
		vtkSmartPointer<vtkDoubleArray> volumeArray = vtkSmartPointer<vtkDoubleArray>::New();
		volumeArray->SetName("Volume");
		volumeArray->InsertNextValue(scaffoldVolume);

		vtkSmartPointer<vtkDoubleArray> porosityArray = vtkSmartPointer<vtkDoubleArray>::New();
		porosityArray->SetName("Porosity");
		porosityArray->InsertNextValue(scaffoldPorosity);

		scaffoldPolyData->GetFieldData()->AddArray(volumeArray);
		scaffoldPolyData->GetFieldData()->AddArray(porosityArray);

		//vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
		writer->SetFileName(scaffoldFilePath.c_str());
		writer->SetInputData(scaffoldPolyData);
		writer->Update();
		writer->SetFileTypeToASCII();
		writer->SetFileVersion(2);
		writer->Write();
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
	uniManager.setUniform(frameShader, "screenSize", glm::vec2(width, height));
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

				if (scaffoldModel) {
					showCutPlane = true;
					//cutPlane = new CutPlane(5.0f);
					cutPlane = std::make_unique<CutPlane>(5.0f);
					cutPlane->center = glm::vec3(containerCenter[0], containerCenter[1], containerCenter[2]);
					cutPlane->initMatrix = glm::translate(
						glm::mat4(1.0),
						glm::vec3(containerCenter[0], containerCenter[1], containerCenter[2]));
				}
			}

			ImGui::EndMenu();
		}

		if (ImGui::Button("Generate Seeds")) {
			_action_generate_seeds();
		}

		if (ImGui::Button("Create Scaffold")) {
			_action_generate_scaffold();
		}

		selectFileButton("Export Scaffold", "../data/", "Export File", ".stl, .vtk");

		// Close the menu bar
		ImGui::EndMainMenuBar();
	}

};

void myGUI::_render_display_settings() {

	static int pickedItem{ -99 };

	ImGui::SetNextWindowSize(ImVec2(200, 250), ImGuiCond_FirstUseEver);
	if (ImGui::Begin("Display Settings", &showDisplaySettingsWin)) {

		// Left side - Selectable items
		{
			ImGui::BeginChild("LeftPanel", ImVec2(150, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
			if (ImGui::Selectable("Mesh", pickedItem == 0)) pickedItem = 0;
			if (ImGui::Selectable("Seeds", pickedItem == 1)) pickedItem = 1;
			if (ImGui::Selectable("Grid", pickedItem == 2)) pickedItem = 2;
			if (ImGui::Selectable("Lighting", pickedItem == 3)) pickedItem = 3;
			if (ImGui::Selectable("Linear Function", pickedItem == 4)) pickedItem = 4;
			if (ImGui::Selectable("Container", pickedItem == 5)) pickedItem = 5;
			ImGui::EndChild();
		}
		ImGui::SameLine();

		// Right side - Color and size controls
		{
			ImGui::BeginGroup();
			{
				ImGui::BeginChild("settings item");

				if (pickedItem == 0) {
					ImGui::Text("Mesh Settings");
					ImGui::ColorEdit3("Mesh Color", (float*)&scaffoldColor);
				}
				else if (pickedItem == 1) {
					ImGui::Text("Seed Settings");
					ImGui::ColorEdit3("Seed Color", (float*)&seedColor);
					ImGui::SliderFloat("Point Size", &seedSize, 0.001f, 1.0f);
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
					ImGui::InputFloat("Ambient Strength", &Ka);
					ImGui::InputFloat("Specular Strength", &Ks);
					//ImGui::SliderFloat("Point Size", &seedSize, 0.001f, 1.0f);
				}
				else if (pickedItem == 4) {
					ImGui::Text("Linear Function Settings");
					ImGui::InputDouble("Max Distance", &maxDist);
				}
				if (pickedItem == 5) {
					ImGui::Text("Container Settings");
					ImGui::ColorEdit4("Container Color", (float*)&containerColor);
				}

				ImGui::EndChild();

				ImGui::SameLine();

				ImGui::EndGroup();
			}
		}
		ImGui::End();

	}

};

void myGUI::_render_cut_tool() {

	if (ImGui::Begin("Plane Cut Tool", &showCutPlane)) {
		if (ImGui::Button("X")) {
			// update previous normal
			cutPlane->prevNormal = cutPlane->normal;
			// update current normal
			cutPlane->normal[0] = 1.0f;
			cutPlane->normal[1] = 0.0f;
			cutPlane->normal[2] = 0.0f;
			cutPlane->offset = 0.0f;
			planeOffset = 0.0f;
			cutPlane->updateModelMatrix();
		}
		ImGui::SameLine();
		if (ImGui::Button("Y")) {
			cutPlane->prevNormal = cutPlane->normal;
			cutPlane->normal[0] = 0.0f;
			cutPlane->normal[1] = 1.0f;
			cutPlane->normal[2] = 0.0f;
			cutPlane->offset = 0.0f;
			planeOffset = 0.0f;
			cutPlane->updateModelMatrix();
		}
		ImGui::SameLine();
		if (ImGui::Button("Z")) {
			cutPlane->prevNormal = cutPlane->normal;
			cutPlane->normal[0] = 0.0f;
			cutPlane->normal[1] = 0.0f;
			cutPlane->normal[2] = 1.0f;
			cutPlane->offset = 0.0f;
			planeOffset = 0.0f;
			cutPlane->updateModelMatrix();
		}
		ImGui::SliderFloat("Offset", &planeOffset, -10.0f, 10.0f);
		cutPlane->offset = planeOffset;
		cutPlane->updateTranslation();
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
				ImGui::InputFloat("Radius", &rMin);
			}

			ImGui::RadioButton("Varied Radius", &radiusOpt, 1);

			if (radiusOpt == 1) {
				ImGui::InputFloat("Min Radius", &rMin);
				ImGui::InputFloat("Max Radius", &rMax);
				if (ImGui::TreeNode("Distance Metric")) {
					ImGui::RadioButton("Distance From Plane", &distFunc, 0);
					if (distFunc == 0) {
						ImGui::InputFloat3("Normal", distPlaneNormal);
						ImGui::InputFloat3("Center", distPlaneCenter);
					};
					ImGui::RadioButton("Distance From Mesh Face", &distFunc, 1);
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

		ImGui::InputFloat("Thickness", &thickness, 0.1f, 1.0f, "%.3f");
		ImGui::SameLine(); help_marker("Thickness of Scaffold");

		ImGui::InputDouble("Hole Scale Factor", &scaleFactor, 0.1, 0.99, "%.3f");
		ImGui::SameLine(); help_marker("To create holes to each face we estimate the maximum inscribed circle, this factor scales its radius. Default value is 0.5");
	}
	ImGui::End();
};

void myGUI::_action_generate_seeds() {

	// ensure that seeds are empty
	seeds.clear();

	// decide boundaries and center of domain
	_update_bounds_center(conOption);

	if (radiusOpt == 0) {
		rMax = rMin;
	}

	// check if distance function is defined from a plane
	if (distFunc == 0) {

		planeDistance = PlaneDistEstimator(
			{ distPlaneCenter[0], distPlaneCenter[1], distPlaneCenter[2] },
			{ distPlaneNormal[0], distPlaneNormal[1], distPlaneNormal[2] });
	}

	// or from the container surface
	else if (distFunc == 1) {

		// if conOption is just a box
		if (conOption == 0) {
			vtkSmartPointer<vtkBox> box = vtkSmartPointer<vtkBox>::New();
			std::cout << xMin << " " << xMax << std::endl;
			box->SetBounds(xMin, xMax, yMin, yMax, zMin, zMax);
			containerDistance = ImplicitFunctionDistEstimator(box);
		}
		else if (conOption == 1) {
			vtkSmartPointer<vtkCylinder> cylinder = vtkSmartPointer<vtkCylinder>::New();
			cylinder->SetRadius(cylinderRadius);

			double center[3] = {
				cylinderHeight * cylinderAxis[0] * 0.5,
				cylinderHeight * cylinderAxis[1] * 0.5,
				cylinderHeight * cylinderAxis[2] * 0.5,
			};

			std::cout << "Center: " << center[0] << " " << center[1] << " " << center[2] << std::endl;

			cylinder->SetCenter(center);

			double axis[3] = {
				cylinderAxis[0], cylinderAxis[1], cylinderAxis[2]
			};
			cylinder->SetAxis(axis);
			containerDistance = ImplicitFunctionDistEstimator(cylinder);
		}
		else {
			meshDistance = MeshDistEstimator(containerMesh);
		}
	}

	// Random seed generation
	if (generateOption == 0) {

		if (conOption == 0) {

			add_log(LogPriority::INFO, "Generating Random Seeds Inside Box.");
			RandomGenerator rg(
				{ xMin, xMax, yMin, yMax, zMin, zMax }, seedNr
			);
			rg.generate_seeds();
			_update_cameras();
			rg.get_seeds(seeds);
		}
		else if (conOption == 1) {

			add_log(LogPriority::INFO, "Generating Random Seeds Inside Cylinder.");
			RandomGenerator rg(
				{ xMin, xMax, yMin, yMax, zMin, zMax }, seedNr
			);
			rg.generate_seeds(
				{ cylinderPt[0], cylinderPt[1], cylinderPt[2] },
				{ cylinderAxis[0], cylinderAxis[1], cylinderAxis[2] },
				cylinderRadius);
			_update_cameras();
			rg.get_seeds(seeds);
		}
		else {

			add_log(LogPriority::INFO, "Generating Random Seeds Inside Mesh Container.");

			RandomGeneratorWall sg(seedNr, containerMesh);
			sg.generate_seeds();
			_update_cameras();
			sg.get_seeds(seeds);
		}
	}

	// Poisson 3D seed generation
	else if (generateOption == 1) {

		if (conOption == 0 || conOption == 1) {
			// determine function to check if a point is inside a container
			std::function<bool(const std::array<double, 3>&)> inside_check;

			// a string to print messages
			std::string suffix{ "" };

			if (conOption == 0) {

				suffix += " inside a Box.";

				inside_check = [&](const std::array<double, 3>& pt) {
					return is_inside_box(pt, { xMin, xMax, yMin, yMax, zMin, zMax });
				};

			}
			else if (conOption == 1) {

				suffix += " inside a Cylinder.";

				inside_check = [&](const std::array<double, 3>& pt) {
					return is_inside_cylinder(
						pt,
						{ cylinderPt[0], cylinderPt[1], cylinderPt[2] },
						{ cylinderAxis[0], cylinderAxis[1], cylinderAxis[2] },
						cylinderRadius, zDim);
				};

			}

			Poisson3D sg(
				rMin, rMax,
				{ containerCenter[0], containerCenter[1], containerCenter[2] },
				{ 0, xDim, 0.0, yDim, 0.0, zDim },
				inside_check
			);

			if (radiusOpt == 0) {
				add_log(LogPriority::INFO, "Generating Seeds Using Uniform Poisson 3D " + suffix);
				sg.generate_seeds();
			}
			else if (radiusOpt == 1 && distFunc == 0) {
				//std::cout << "generating with Poisson 3D and dist from plane" << std::endl;
				std::string tag = "Generating Seeds Using Varied Poisson 3D " + suffix + " Distance is measured from plane.";
				add_log(LogPriority::INFO, tag);
				sg.generate_seeds(planeDistance, linearFunc);
			}
			else if (radiusOpt == 1 && distFunc == 1) {
				std::string tag = "Generating Seeds Using Varied Poisson 3D " + suffix + " Distance is measured from mesh outer face.";
				add_log(LogPriority::INFO, tag);
				sg.generate_seeds(containerDistance, linearFunc);
			}
			_update_cameras();
			sg.get_seeds(seeds);
		}

		// Poisson 3d inside a container
		if (conOption == 2) {
			Poisson3DWall sg(containerMesh, rMin, rMax, neighbors, { 10, 10, 10 }, wallResolution);
			if (radiusOpt == 0) {
				add_log(LogPriority::INFO, "Generating seeds using Uniform Poisson 3D inside Mesh Container");
				sg.generate_seeds();
			}
			else if (radiusOpt == 1 && distFunc == 0) {
				add_log(LogPriority::INFO, "Generating Seeds Using Varied Poisson 3D Inside Mesh Container. Distance is measured from plane.");
				sg.generate_seeds(planeDistance, linearFunc);
			}
			else if (radiusOpt == 1 && distFunc == 1) {
				add_log(LogPriority::INFO, "Generating Seeds Using Varied Poisson 3D Inside Mesh Container. Distance is measured from mesh outer face.");
				std::cout << "generating with Poisson 3D and dist from mesh" << std::endl;
				sg.generate_seeds(meshDistance, linearFunc);
			}
			_update_cameras();
			sg.get_radii(radii);
			sg.get_seeds(seeds);
		}
	}

	// create a bounding box
	//box = new BBox(xMin, xMax, yMin, yMax, zMin, zMax);
	box = std::make_unique<BBox>(xMin, xMax, yMin, yMax, zMin, zMax);

	// create the container if it is a mesh
	if (conOption == 2) {
		//containerModel = new Model(containerFile);
		containerModel = std::make_unique<Model>(containerFile);
	}

	//seedObj = new VisualizeSeeds(seeds);
	seedObj = std::make_unique<VisualizeSeeds>(seeds);
	showSeeds = true;
};

void myGUI::_action_generate_scaffold() {


	add_log(LogPriority::INFO, "Generating Scaffold");

	if (conOption == 0) {
		domainVolume = xDim * yDim * zDim;
	}

	else if (conOption == 1) {
		domainVolume = std::pow(cylinderRadius, 2) * 3.14f * cylinderHeight;
	}

	else {
		vtkNew<vtkMassProperties> massProperties;
		massProperties->SetInputData(containerMesh);
		massProperties->Update();
		domainVolume = massProperties->GetVolume();
	}

	// define the model file name
	scaffoldFilePath = "../data/" + std::string(version) + ".stl";
	scaffoldFileName = std::string(version) + ".stl";

	if (!runVolumeOptimization) {
		if (conOption == 0 || conOption == 1) {
			std::cout << "Generating Random Seeds Mesh inside rectangular" << std::endl;
			ScaffoldGeneratorBox sgb(seeds, { xMin, xMax, yMin, yMax, zMin, zMax }, { nX, nY, nZ }, edgeSize, scaleFactor);

			if (conOption == 1) {
				sgb.add_cylindrical_wall(
					static_cast<double>(cylinderPt[0]),
					static_cast<double>(cylinderPt[1]),
					static_cast<double>(cylinderPt[2]),
					static_cast<double>(cylinderAxis[0]),
					static_cast<double>(cylinderAxis[1]),
					static_cast<double>(cylinderAxis[2]),
					cylinderRadius
				);
			}

			sgb.generate_voro(regSteps);
			sgb.get_seeds(seeds);
			//sgb.generate_mesh(thickness, scaffoldFilePath);
			if (customResolutionFlag) {
				sgb.generate_mesh(thickness, scaffoldPolyData, { resolution[0], resolution[1], resolution[2] });
			}
			else {
				sgb.generate_mesh(thickness, scaffoldPolyData, {});
			}
			_update_cameras();
		}

		else if (conOption == 2) {
			std::cout << "Generating Random Seeds Mesh inside Container" << std::endl;
			//wallResolution = rMin / 2.0f;
			ScaffoldGeneratorWall sgb(seeds, containerMesh, { nX, nY, nZ }, neighbors, wallResolution, edgeSize, scaleFactor);
			sgb.generate_voro(regSteps);
			sgb.get_seeds(seeds);
			//sgb.generate_mesh(thickness, scaffoldFilePath);
			if (customResolutionFlag) {
				sgb.generate_mesh(thickness, scaffoldPolyData, { resolution[0], resolution[1], resolution[2] });
			}
			else {
				sgb.generate_mesh(thickness, scaffoldPolyData, {});
			}
			seedObj = std::make_unique<VisualizeSeeds>(seeds);
			_update_cameras();
		}
	}

	else if (runVolumeOptimization) {
		std::cout << "Generating with Volume Optimization" << std::endl;

		if (seeds.empty()) {
			std::cerr << "First Populate Seeds" << std::endl;
		}
		int seedSize = seeds.size();

		wInit.resize(seedSize);
		wInit.setZero();
		std::cout << " -------------------- " << std::endl;

		targetVols.resize(seedSize);

		// case 1: volumes all equal -> divide domain volume per seed nr
		if (volOption == 0) {

			for (int i = 0; i < seedSize; i++) {
				targetVols[i] = domainVolume / seedSize;
			}
		}

		// case 2: volumes are decided from Poisson 3d radii 

		else if (volOption == 1) {

			double volSum{ 0.0 };
			for (int i{ 0 }; i < radii.size(); i++) {
				targetVols[i] = std::pow(radii[i], 3);
				volSum += targetVols[i];
			}
			// Normalize
			for (int i{ 0 }; i < seeds.size(); i++) {
				targetVols[i] = (targetVols[i] / volSum) * domainVolume;
			}

			double targetSum = 0;
			for (int i{ 0 }; i < seeds.size(); i++) {
				targetSum += targetVols[i];
			}
			std::cout << "domain volume: " << domainVolume << std::endl;
			std::cout << "target volume sum: " << targetSum << std::endl;
			// std::cerr << "Not implemented yet!" << std::endl;
		}

		else if (volOption == 2) {
			std::cerr << "Not implemented yet!" << std::endl;
		}

		// just a simple box container
		if (conOption == 0) {

			std::cout << "Vol Opt inside rect Domain" << std::endl;
			VolOpt vo(
				seeds,
				targetVols,
				wInit,
				{ xMin, xMax, yMin, yMax, zMin, zMax }
			);
			vo.loop(regSteps);
			vo.get_seeds(seeds);
			//vo.generate_mesh(thickness, scaffoldFilePath);
			if (customResolutionFlag) {
				vo.generate_mesh(thickness, scaffoldPolyData, { resolution[0], resolution[1], resolution[2] });
			}
			else {
				vo.generate_mesh(thickness, scaffoldPolyData, {});
			}
		}

		// volume optimization inside a mesh wall
		else if (conOption == 2) {

			std::cout << "Vol Opt inside Mesh" << std::endl;

			VolOptWall vo(
				seeds,
				targetVols,
				wInit,
				containerMesh
			);
			vo.loop(regSteps);
			vo.get_seeds(seeds);
			if (customResolutionFlag) {
				vo.generate_mesh(thickness, scaffoldFilePath, { resolution[0], resolution[1], resolution[2] });
			}
			else {
				vo.generate_mesh(thickness, scaffoldFilePath, {});
			}
		}

	}

	// create the scaffold mesh

	scaffoldModel = std::make_unique<Model>(scaffoldPolyData);
	scaffoldModel->get_details(faceNr, vertexNr, edgeNr, scaffoldVolume, scaffoldPorosity);
	scaffoldPorosity = scaffoldVolume / domainVolume;

	add_log(LogPriority::INFO, "Scaffold Ready");

};

void myGUI::_create_dockspace() {

	ImGui::DockSpaceOverViewport(0, ImGui::GetMainViewport());

	_render_main_menu_bar();

};

void myGUI::_render_visualizer() {

	ImGui::SetNextWindowSize(ImVec2(800, 600), ImGuiCond_Once);
	ImGui::Begin("Visualizer");
	
	ImVec2 availableSize = ImGui::GetWindowSize();
	ImVec2 pos = ImGui::GetWindowPos();

	if ((int)availableSize.x != framebuffer.width || (int)availableSize.y != framebuffer.height) {

		//std::cout << "inside visualizer" << std::endl;
		framebuffer.width = (int)availableSize.x;
		framebuffer.height = (int)availableSize.y;
		framebuffer.posx = pos.x;
		framebuffer.posy = pos.y;
		create_frame_buffer(framebuffer);
	}

	int vw = static_cast<int>(availableSize.x);
	int vh = static_cast<int>(availableSize.y);
	int vx = static_cast<int>(pos.x);
	int vy = static_cast<int>(pos.y);

	trackCamera->set_viewport(vx, vy, vw, vh);

	ImGui::Image((ImTextureID)(intptr_t)framebuffer.textureId, availableSize, ImVec2(0, 1), ImVec2(1, 0));

	if (scaffoldModel && showScaffold) {
		ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav;
		const float PAD = 10.0f;
		//const ImGuiViewport* viewport = ImGui::GetMainViewport();
		//ImVec2 work_pos = viewport->WorkPos; // Use work area to avoid menu-bar/task-bar, if any!
		//ImVec2 work_size = viewport->WorkSize;
		ImVec2 window_pos, window_pos_pivot;
		//ImVec2 work_pos;
		window_pos.x = framebuffer.posx + framebuffer.width - PAD;
		window_pos.y = framebuffer.posy - framebuffer.height;
		window_pos_pivot.x = 1.0f;
		window_pos_pivot.y = 1.0f;
		ImGui::SetNextWindowPos(window_pos, ImGuiCond_Always, window_pos_pivot);
		//ImGui::SetNextWindowViewport(viewport->ID);
		window_flags |= ImGuiWindowFlags_NoMove;
		ImGui::SetNextWindowBgAlpha(0.35f); // Transparent background
		if (ImGui::Begin("Example: Simple overlay", nullptr, window_flags))
		{
			ImGui::Text("Mesh Details");
			ImGui::Separator();
			ImGui::Text("Mesh Name: %s", scaffoldFileName.c_str());
			ImGui::Text("Mesh Vertices: %d", vertexNr);
			ImGui::Text("Mesh Faces: %d", faceNr);
			ImGui::Text("Volume: %.3f", scaffoldVolume);
			ImGui::Text("Porosity: %.3f", scaffoldPorosity);
			ImGui::End();
		}
	}


	ImGui::End();

};
