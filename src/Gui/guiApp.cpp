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

// vtk
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
#include <vtkImageData.h>
#include <vtkImageData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkMetaImageWriter.h>

// custom
#include "guiApp.h"
#include <ImGuiFileDialog/ImGuiFileDialog.h>
#include "buildScaffold.h"
#include "Utils/Utils.h"
#include "ScaffoldGenerator/ScaffoldGenerator.h"
#include "SeedGenerator/DistanceCalculator.h"
#include "OpenGlRender/Model.h"
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
		std::filesystem::absolute(std::filesystem::path("shaders/mainVShader.vertexshader")).string().c_str(),
		std::filesystem::absolute(std::filesystem::path("shaders/mainFShader.fragmentshader")).string().c_str(),
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
		"shaders/normalVShader.vs",
		"shaders/normalFShader.fs",
		"shaders/normalGShader.gs"
	);

	uniManager.add_uniform(normalShader, "view");
	uniManager.add_uniform(normalShader, "model");
	uniManager.add_uniform(normalShader, "projection");
	uniManager.add_uniform(normalShader, "normalColor");

	edgeShader = Shader(
		"shaders/edgeVShader.vs",
		"shaders/edgeFShader.fs",
		NULL
	);
	uniManager.add_uniform(edgeShader, "view");
	uniManager.add_uniform(edgeShader, "model");
	uniManager.add_uniform(edgeShader, "projection");
	uniManager.add_uniform(edgeShader, "cutPlaneCoeffs");
	uniManager.add_uniform(edgeShader, "cutPlane");
	uniManager.add_uniform(edgeShader, "edgeColor");

	gridShader = Shader(
		"shaders/gridVShader.vs",
		"shaders/gridFShader.fs",
		"shaders/gridGShader.gs"
	);
	uniManager.add_uniform(gridShader, "view");
	uniManager.add_uniform(gridShader, "model");
	uniManager.add_uniform(gridShader, "projection");
	uniManager.add_uniform(gridShader, "gridColor");

	seedShader = Shader(
		"shaders/seedVShader.vertexshader",
		"shaders/seedFShader.fragmentshader",
		NULL
	);
	uniManager.add_uniform(seedShader, "projection");
	uniManager.add_uniform(seedShader, "view");
	uniManager.add_uniform(seedShader, "seedColor");
	uniManager.add_uniform(seedShader, "seedSize");

	cutShader = Shader(
		"shaders/cutVShader.vs",
		"shaders/cutFShader.fs",
		NULL
	);
	uniManager.add_uniform(cutShader, "projection");
	uniManager.add_uniform(cutShader, "view");
	uniManager.add_uniform(cutShader, "model");
	uniManager.add_uniform(cutShader, "minBounds");
	uniManager.add_uniform(cutShader, "maxBounds");

	bboxShader = Shader(
		"shaders/bboxVShader.vs",
		"shaders/bboxFShader.fs",
		NULL
	);
	uniManager.add_uniform(bboxShader, "view");
	uniManager.add_uniform(bboxShader, "model");
	uniManager.add_uniform(bboxShader, "projection");

	containerShader = Shader(
		"shaders/containerVShader.vs",
		"shaders/containerFShader.fs",
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
		std::filesystem::relative(std::filesystem::path("shaders/frameVShader.vs")).string().c_str(),
		std::filesystem::relative(std::filesystem::path("shaders/frameFShader.fs")).string().c_str(),
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
};

void myGUI::_init_imgui() {

	//std::cout << "Starting imgui" << std::endl;

	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();

	io = ImGui::GetIO(); (void)io;
	io.Fonts->AddFontFromFileTTF("share/fonts/DroidSans.ttf", 20.0f);
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

	ImGui::LoadIniSettingsFromDisk("share/imgui.ini");

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

		//std::cout << xPos << " " << yPos << std::endl;

		//bool cameraUpdateFlags = !ImGuiFileDialog::Instance()->IsOpened() && !ImGui::IsAnyItemActive() && !ImGui::IsAnyItemHovered();

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

		scaffoldShader.use();
		// pass light
		uniManager.setUniform(scaffoldShader, "lightColor", lightColor[0], lightColor[1], lightColor[2]);
		uniManager.setUniform(scaffoldShader, "lightPosWorld", lightPosCamera[0], lightPosCamera[1], lightPosCamera[2]);
		//uniManager.setUniform(scaffoldShader, "lightPosWorld", cameraPos.x, cameraPos.y, cameraPos.z);
		uniManager.setUniform(scaffoldShader, "Ka", Ka);
		uniManager.setUniform(scaffoldShader, "Ks", Ks);

		if (showCutPlane) {

			glDisable(GL_CULL_FACE);

			float offset = cutPlane->offset;
			glm::vec3 normal = cutPlane->normal;

			float d = -glm::dot(normal, cutPlane->center + normal * offset);

			planeCoeffs = glm::vec4(normal, d);

			uniManager.setUniform(scaffoldShader, "cutPlane", 1);
			uniManager.setUniform(scaffoldShader, "cutPlaneCoeffs", planeCoeffs);

			glEnable(GL_CULL_FACE);
		}
		else {
			uniManager.setUniform(scaffoldShader, "cutPlane", 0);
		}

		if (scaffoldReady) {

			scaffoldModel = std::make_unique<Model>(scaffoldPolyData);
			scaffoldModel->get_details(faceNr, vertexNr, edgeNr, scaffoldVolume, scaffoldPorosity);
			scaffoldPorosity = scaffoldVolume / domainVolume;

			add_log(LogPriority::SUCCESS, "Scaffold Ready");

			scaffoldReady = false;

			if (scaffoldGenerationThread.joinable()) {
				scaffoldGenerationThread.join(); 
			}
		}

		if (scaffoldModel){

			uniManager.setUniform(scaffoldShader, "projection", projection);
			uniManager.setUniform(scaffoldShader, "view", view);
			uniManager.setUniform(scaffoldShader, "model", glm::mat4(1.0));
			uniManager.setUniform(
				scaffoldShader, "objectColor", scaffoldColor[0], scaffoldColor[1], scaffoldColor[2], scaffoldColor[3]);
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

		if (ImGui::IsKeyPressed(ImGuiKey_M)) {
			showEdges = !showEdges;
		}

		if (ImGui::IsKeyPressed(ImGuiKey_G)) {
			showGrid = !showGrid;
		}

		if (ImGui::IsKeyPressed(ImGuiKey_N)) {
			showNormals = !showNormals;
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
			uniManager.setUniform(edgeShader, "edgeColor", 0.0f, 0.0f, 0.0f, 1.0f);
			if (showScaffold) {
				scaffoldModel->draw_edges();
			}
		}

		if (showPoreNetwork) {

			if (activePoreNetworkIndex != -1) {
				auto pnw = &poreNetworkList[activePoreNetworkIndex];

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
				pnw->model->draw();
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

		if (showContainer) {
			if (conOption == 0 || conOption == 1) {
				//auto cb = static_cast<BoxContainer*>(container.obj);
				bboxShader.use();
				uniManager.setUniform(bboxShader, "projection", projection);
				uniManager.setUniform(bboxShader, "view", view);
				uniManager.setUniform(bboxShader, "model", model);
				activeContainer->render();
			}
			else{
				containerShader.use();
				uniManager.setUniform(containerShader, "projection", projection);
				uniManager.setUniform(containerShader, "view", view);
				uniManager.setUniform(containerShader, "model", model);
				containerShader.use();
				uniManager.setUniform(containerShader, "lightColor", lightColor[0], lightColor[1], lightColor[2]);
				uniManager.setUniform(containerShader, "lightPosWorld", lightPosCamera[0], lightPosCamera[1], lightPosCamera[2]);
				uniManager.setUniform(containerShader, "Ka", Ka);
				uniManager.setUniform(
					containerShader, "objectColor",
					containerColor[0], containerColor[1], containerColor[2], containerColor[3]);
				activeContainer->render();
				//if (conOption == 1) {
				//	//auto cc = static_cast<CylinderContainer*>(container.obj);
				//	//cc->model->draw();
				//	activeContainer->render();
				//}
				//else {
				//	auto cc = static_cast<MeshContainer*>(container.obj);
				//	cc->model->draw();
				//}
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
			//std::cout << "Draw distance plane " << std::endl;
			glm::mat4 modelDist = distPlane->tMatrix * distPlane->initMatrix * distPlane->rotMatrix;
			cutShader.use();
			uniManager.setUniform(cutShader, "projection", projection);
			uniManager.setUniform(cutShader, "view", view);
			uniManager.setUniform(cutShader, "model", modelDist);
			uniManager.setUniform(cutShader, "minBounds", 0, 0, 0);
			uniManager.setUniform(cutShader, "maxBounds", xDim, yDim, zDim);
			distPlane->draw();
		}

		if (ImGui::IsKeyPressed(ImGuiKey::ImGuiKey_Delete)) {
			
			add_log(LogPriority::INFO, "Deleting everything!");

			container.obj = nullptr;
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

			activePoreNetworkIndex = -1;
			poreNetworkList.clear();
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

			// Process the file, open, load etc.
			if (loadedMesh == scaffold) {

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
				//float xc = static_cast<float>((xMax + xMin) * 0.5f);
				//float yc = static_cast<float>((yMax + yMin) * 0.5f);
				//float zc = static_cast<float>((zMax + zMin) * 0.5f);

				//float dx = xMax - xMin;
				//float dy = yMax - yMin;
				//float dz = zMax - zMin;
				//float diagonal = std::sqrt(dx * dx + dy * dy + dz * dz);
				//float distance = 1.1f * diagonal;

				//containerCenter[0] = xc;
				//containerCenter[1] = yc;
				//containerCenter[2] = zc;

				//trackCamera->set_position(xc, yc, zc + distance);

				// change camera target and camera position
				_update_cameras();

				scaffoldModel->get_details(faceNr, vertexNr, edgeNr, scaffoldVolume, scaffoldPorosity);

				add_log(LogPriority::SUCCESS, "Loading Scaffold Succesfull");
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

	if (ImGuiFileDialog::Instance()->Display("Export Scaffold")) {
		if (ImGuiFileDialog::Instance()->IsOk()) { // action if OK
			scaffoldFilePath = ImGuiFileDialog::Instance()->GetFilePathName();
			scaffoldFileName = ImGuiFileDialog::Instance()->GetCurrentFileName();
			_export_mesh();

		}
		ImGuiFileDialog::Instance()->Close();
	}

	//_render_console();

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

void myGUI::_generate_scaffold() {
	
	//std::cout << wInit << std::endl;
	//std::cout << " -------------------- " << std::endl;
	//std::cout << targetVols << std::endl;

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

		ImGui::RadioButton("Cylindrical Container", &conOption, 1);
		ImGui::SameLine(); help_marker("A cylindrical container, starting from (0, 0, 0). Select Axis, Radius. The Height depends on the Bounds");

		ImGui::RadioButton("STL Mesh", &conOption, 2);
		ImGui::SameLine(); help_marker("A custom container geometry. In this case the bounding box of the mesh is used for the voro++ container and the scaffold is built inside the mesh domain");

		if (conOption != prevConOption) {
			switch(conOption){
			case 0: {
				activeContainer = &boxContainer; 
				break;
			}
			case 1: activeContainer = &cylContainer; break;
			case 2: activeContainer = &meshContainer; break;
			}
			activeContainer->create();
			_update_cameras();
			prevConOption = conOption;
		}
		activeContainer->gui_setup();
	}
	ImGui::End();
};

//void myGUI::framebuffer_size_callback_imp(int width, int height) {
//	defCamera.set_window_size(width, height);
//	trackCamera.set_window_size(width, height);
//};

void myGUI::_update_cameras() {

	float xc = static_cast<float>((xMax + xMin) * 0.5f);
	float yc = static_cast<float>((yMax + yMin) * 0.5f);
	float zc = static_cast<float>((zMax + zMin) * 0.5f);

	float dx = xMax - xMin;
	float dy = yMax - yMin;
	float dz = zMax - zMin;
	float diagonal = std::sqrt(dx * dx + dy * dy + dz * dz);
	float distance = 1.1f * diagonal;

	container.center[0] = xc;
	container.center[1] = yc;
	container.center[2] = zc;

	cameraUpdate = true;
	cameraTarget = glm::vec3(container.center[0], container.center[1], container.center[2]);
	defCamera->position = cameraPos;
	defCamera->target = cameraTarget;
	//trackCamera->position = cameraPos;
	trackCamera->set_target(cameraTarget.x, cameraTarget.y, cameraTarget.z);
	trackCamera->set_position(xc, yc, zc + distance);
	trackCamera->update();
	projection = trackCamera->get_projection_matrix();
	view = trackCamera->get_view_matrix();
	cameraUpdate = false;
}

void myGUI::_update_bounds_center(int& conOption) {

	const Bounds bounds = activeContainer->compute_bounds();

	xMin = bounds.xMin; xMax = bounds.xMax;
	yMin = bounds.yMin; yMax = bounds.yMax;
	zMin = bounds.zMin; zMax = bounds.zMax;

	container.center = {
		static_cast<double>(bounds.center[0]),
		static_cast<double>(bounds.center[1]),
		static_cast<double>(bounds.center[2]),
	};
};

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

void myGUI::_export_binary_image(
	const std::string& filepath, const float& voxelSize, const int& backVal, const int& forVal) {

	vtkSmartPointer<vtkImageData> blackImage = vtkSmartPointer<vtkImageData>::New();
	double bounds[6];
	scaffoldPolyData->GetBounds(bounds);

	// set spacing
	double spacing[3];
	spacing[0] = voxelSize;
	spacing[1] = voxelSize;
	spacing[2] = voxelSize;
	blackImage->SetSpacing(spacing);

	// compute dimensions
	int dim[3];
	for (int i = 0; i < 3; i++)
	{
		dim[i] = static_cast<int>(
			std::ceil((bounds[i * 2 + 1] - bounds[i * 2]) / voxelSize));
	}
	blackImage->SetDimensions(dim);
	blackImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

	// set orgin
	double origin[3];
	origin[0] = bounds[0] + spacing[0] / 2;
	origin[1] = bounds[2] + spacing[1] / 2;
	origin[2] = bounds[4] + spacing[2] / 2;
	blackImage->SetOrigin(origin);
	blackImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

	// fill with black color
	vtkIdType count = blackImage->GetNumberOfPoints();
	for (vtkIdType i = 0; i < count; ++i)
	{
		blackImage->GetPointData()->GetScalars()->SetTuple1(i, forVal);
	}

	// polygonal data --> image stencil:
	vtkNew<vtkPolyDataToImageStencil> pol2stenc;
	pol2stenc->SetInputData(scaffoldPolyData);
	pol2stenc->SetOutputOrigin(blackImage->GetOrigin());
	pol2stenc->SetOutputSpacing(blackImage->GetSpacing());
	pol2stenc->SetOutputWholeExtent(blackImage->GetExtent());
	pol2stenc->Update();

	// Cut the corresponding white image and set the background:
	vtkNew<vtkImageStencil> imgstenc;
	imgstenc->SetInputData(blackImage);
	imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
	imgstenc->ReverseStencilOff();
	imgstenc->SetBackgroundValue(backVal);
	imgstenc->Update();

	vtkNew<vtkMetaImageWriter> writer;
	writer->SetFileName(filepath.c_str());
	writer->SetInputData(imgstenc->GetOutput());
	writer->Write();
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

				if (scaffoldModel) {
					showCutPlane = true;
					//cutPlane = new CutPlane(5.0f);
					cutPlane = std::make_unique<CutPlane>(5.0f);
					cutPlane->center = glm::vec3(container.center[0], container.center[1], container.center[2]);
					cutPlane->initMatrix = glm::translate(
						glm::mat4(1.0),
						glm::vec3(container.center[0], container.center[1], container.center[2]));
				}
			}

			ImGui::EndMenu();
		}

		if (ImGui::Button("Generate Seeds")) {

			_action_generate_seeds();
		}

		if (ImGui::Button("Update Voronoi")) {
			_action_update_voronoi();
		}

		if (ImGui::Button("Estimate Threshold")) {
			_action_estimate_connectivity();
		}

		if (ImGui::Button("Create Scaffold")) {

			if (seeds.empty()) {

				add_log(LogPriority::ERROR, "First Populate Seeds!");

			}
			else {

				if (scaffoldGenerationThread.joinable()) {
					scaffoldGenerationThread.join(); // join any previous thread first
				}

				scaffoldReady = false;

				scaffoldGenerationThread = std::thread([this]() {

					// keep time
					auto start_time = std::chrono::steady_clock::now();

					scaffoldProgress = 0.0f;
					_action_generate_scaffold();
					scaffoldReady = true;

					// --- TIMER END & LOGGING ---
					auto end_time = std::chrono::steady_clock::now();

					// Calculate the duration in milliseconds
					auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

					// Format the output string
					std::ostringstream oss;
					oss << "Scaffold creation finished in "
						<< std::fixed << std::setprecision(3) // Set precision to 3 decimal places
						<< duration_ms.count() / 1000.0   // Convert ms to seconds
						<< " seconds.";

					// Add the formatted string to your log
					add_log(LogPriority::INFO, oss.str());
				});
			}			
		}

		//selectFileButton("Export Scaffold", "../data/", "Export File", ".stl, .vtk");

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
			if (ImGui::Selectable("Pore Network", pickedItem == 6)) pickedItem = 6;
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
					ImGui::ColorEdit4("Mesh Color", (float*)&scaffoldColor);
				}
				else if (pickedItem == 1) {
					ImGui::Text("Seed Settings");
					ImGui::ColorEdit3("Seed Color", (float*)&seedColor);
					ImGui::InputFloat("Point Size", &seedSize, 0.001f, 1.0f);
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
					ImGui::InputFloat("Ambient Strength", &Ka, 0.01f, 100.0f);
					ImGui::InputFloat("Specular Strength", &Ks, 0.01f, 100.0f);
				}
				else if (pickedItem == 4) {
					ImGui::Text("Linear Function Settings");
					ImGui::InputDouble("Max Distance", &maxDist);
				}
				if (pickedItem == 5) {
					ImGui::Text("Container Settings");
					ImGui::ColorEdit4("Container Color", (float*)&containerColor);
				}

				if (pickedItem == 6) {
					ImGui::Text("Pore Network Display Settings");
					ImGui::ColorEdit4("Line Color", renderSettings.poreNetworkColor.data());
					ImGui::InputFloat("Line Size", &renderSettings.poreNetworkLineSize, 0.1f);
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
				ImGui::SetNextItemWidth(200);
				ImGui::InputFloat("Radius", &rMin);
			}

			ImGui::RadioButton("Varied Radius", &radiusOpt, 1);

			if (radiusOpt == 1) {
				ImGui::SetNextItemWidth(200);
				ImGui::InputFloat("Min Radius", &rMin);
				ImGui::SetNextItemWidth(200);
				ImGui::InputFloat("Max Radius", &rMax);
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
		ImGui::SliderFloat("Connectivity Threshold", &connectivityThreshold, 0.0f, 1.0f, "%.3f");
		ImGui::SameLine(); help_marker("0 full face open - 1 maximum pullback");

		connectivityThreshold = std::clamp(connectivityThreshold, 0.001f, 0.999f);

		//ImGui::InputDouble("Hole Scale Factor", &scaleFactor, 0.1, 0.99, "%.3f");
		//ImGui::SameLine(); help_marker("To create holes to each face we estimate the maximum inscribed circle, this factor scales its radius. Default value is 0.5");
	}
	ImGui::End();
};

void myGUI::_action_generate_seeds() {

	// ensure that seeds are empty
	seeds.clear();

	// decide boundaries and center of domain
	_update_bounds_center(conOption);

	std::array<float, 6> bounds = { xMin, xMax, yMin, yMax, zMin, zMax };

	ContainerAdapter adapter = { *activeContainer, xDim, yDim, zDim };

	// the configuration struct
	RunConfig cfg;

	// decide the message to the console
	switch (conOption) {
	case 0: {
		suffix = "Inside Box.";
		break;
	}
	case 1: {
		suffix = "Inside Cylinder.";
		break;
	}
	case 2: {
		suffix = "Inside Mesh Container.";
		break;
	}
	};

	// just random
	if (generateOption == 0) {

		add_log(LogPriority::INFO, "Generating Random Seeds " + suffix);

		// create the seeds using the adapter
		Random rnd(seedNr);
		rnd.run(adapter, seeds);
		_update_cameras();
	}
	else {
		// uniform sampling
		if (radiusOpt == 0) {

			add_log(LogPriority::INFO, "Generating Seeds Using Uniform Poisson 3D " + suffix);

			// create the generator
			Poisson3D gnr(rMin, rMax, { container.center[0], container.center[1], container.center[2] }, 30);
			gnr.run(adapter, seeds);
			_update_cameras();
		}
		// varied sampling
		else {

			// get the distance from container
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

			add_log(LogPriority::INFO, "Generating Seeds Using Varied Poisson 3D " + suffix);

			Poisson3D gnr(rMin, rMax, { container.center[0], container.center[1], container.center[2] }, 30);
			
			switch (distFunc) {
			// distance from plane
			case 0: {
				std::array<double, 3> center = { distPlaneCenter[0], distPlaneCenter[1], distPlaneCenter[2] };
				std::array<double, 3> normal = { distPlaneNormal[0], distPlaneNormal[1], distPlaneNormal[2] };
				cfg.dist = std::make_shared<PlaneDistEstimator>(center, normal);
				break;
				}
			// distance from container surface
			case 1: {
				cfg.dist = activeContainer->get_distance_estimator();
				break;
				}
			// distance from point
			case 2: {

				std::array<double, 3> point = {
					distPoint[0], distPoint[1], distPoint[2]
				};

				cfg.dist = std::make_shared<PointDistEstimator>(point);
				break; 
				}
			}
			
			gnr.run(adapter, seeds, cfg);
			_update_cameras();
		}
	}

	//if (radiusOpt == 0) {
	//	rMax = rMin;
	//}

	//// check if distance function is defined from a plane
	//if (distFunc == 0) {

	//	planeDistance = PlaneDistEstimator(
	//		{ distPlaneCenter[0], distPlaneCenter[1], distPlaneCenter[2] },
	//		{ distPlaneNormal[0], distPlaneNormal[1], distPlaneNormal[2] });
	//}

	//// or from the container surface
	//else if (distFunc == 1) {


	//	switch (conOption) {
	//	case 0: {
	//		vtkSmartPointer<vtkBox> box = vtkSmartPointer<vtkBox>::New();
	//		box->SetBounds(xMin, xMax, yMin, yMax, zMin, zMax);
	//		containerDistance = ImplicitFunctionDistEstimator(box);
	//		break;
	//	}
	//	case 1: {
	//		vtkSmartPointer<vtkCylinder> cylinder = vtkSmartPointer<vtkCylinder>::New();
	//		cylinder->SetRadius(cylinderRadius);
	//		
	//		break;
	//	}
	//	case 2: {

	//		break;
	//	}


	//	// if conOption is just a box
	//	if (conOption == 0) {
	//		vtkSmartPointer<vtkBox> box = vtkSmartPointer<vtkBox>::New();
	//		box->SetBounds(xMin, xMax, yMin, yMax, zMin, zMax);
	//		containerDistance = ImplicitFunctionDistEstimator(box);
	//	}
	//	else if (conOption == 1) {
	//		vtkSmartPointer<vtkCylinder> cylinder = vtkSmartPointer<vtkCylinder>::New();
	//		cylinder->SetRadius(cylinderRadius);

	//		double center[3] = {
	//			cylinderHeight * cylinderAxis[0] * 0.5,
	//			cylinderHeight * cylinderAxis[1] * 0.5,
	//			cylinderHeight * cylinderAxis[2] * 0.5,
	//		};

	//		std::cout << "Center: " << center[0] << " " << center[1] << " " << center[2] << std::endl;

	//		cylinder->SetCenter(center);

	//		double axis[3] = {
	//			cylinderAxis[0], cylinderAxis[1], cylinderAxis[2]
	//		};
	//		cylinder->SetAxis(axis);
	//		containerDistance = ImplicitFunctionDistEstimator(cylinder);
	//	}
	//	else {
	//		meshDistance = MeshDistEstimator(containerMesh);
	//	}
	//}

	//// Random seed generation
	//if (generateOption == 0) {

	//	if (conOption == 0) {

	//		inside_check = [&](const std::array<double, 3>& pt) {
	//			return is_inside_box(pt, bounds);
	//		};

	//		add_log(LogPriority::INFO, "Generating Random Seeds Inside Box.");
	//		auto rg = std::make_unique<RandomGenerator>(
	//			bounds, seedNr
	//		);
	//		//rg.generate_seeds();
	//		//_update_cameras();
	//		//rg.get_seeds(seeds);
	//		seedGenerator = std::move(rg);
	//	}

	//	else if (conOption == 1) {

	//		auto cc = static_cast<CylinderContainer*>(container.obj);

	//		inside_check = [&](const std::array<double, 3>& pt) {
	//			return is_inside_cylinder(
	//				pt,
	//				{ cc->cylinderPt[0], cc->cylinderPt[1], cc->cylinderPt[2] },
	//				{ cc->cylinderAxis[0], cc->cylinderAxis[1], cc->cylinderAxis[2] },
	//				cc->cylinderRadius, zDim);
	//		};

	//		add_log(LogPriority::INFO, "Generating Random Seeds Inside Cylinder.");
	//		auto rg = std::make_unique<RandomGeneratorCylinder>(
	//			bounds, seedNr
	//		);
	//		//rg.generate_seeds(
	//		//	{ cylinderPt[0], cylinderPt[1], cylinderPt[2] },
	//		//	{ cylinderAxis[0], cylinderAxis[1], cylinderAxis[2] },
	//		//	cylinderRadius);
	//		_update_cameras();
	//		//rg.get_seeds(seeds);
	//		seedGenerator = std::move(rg);
	//	}
	//	else {

	//		add_log(LogPriority::INFO, "Generating Random Seeds Inside Mesh Container.");
	//		auto con = static_cast<MeshContainer*>(container.obj);
	//		auto rg = std::make_unique<RandomGeneratorWall>(seedNr, con->containerMesh);
	//		_update_cameras();
	//	}
	//}

	//// Poisson 3D seed generation
	//else if (generateOption == 1) {

	//	if (conOption == 0 || conOption == 1) {

	//		// a string to print messages
	//		std::string suffix{ "" };

	//		if (conOption == 0) {

	//			suffix += " inside a Box.";

	//			inside_check = [&](const std::array<double, 3>& pt) {
	//				return is_inside_box(pt, { xMin, xMax, yMin, yMax, zMin, zMax });
	//			};

	//		}
	//		else if (conOption == 1) {

	//			suffix += " inside a Cylinder.";

	//			inside_check = [&](const std::array<double, 3>& pt) {
	//				return is_inside_cylinder(
	//					pt,
	//					{ cylinderPt[0], cylinderPt[1], cylinderPt[2] },
	//					{ cylinderAxis[0], cylinderAxis[1], cylinderAxis[2] },
	//					cylinderRadius, zDim);
	//			};
	//		}

	//		Poisson3D sg(
	//			rMin, rMax,
	//			{ containerCenter[0], containerCenter[1], containerCenter[2] },
	//			{ 0, xDim, 0.0, yDim, 0.0, zDim },
	//			inside_check
	//		);

	//		if (radiusOpt == 0) {
	//			add_log(LogPriority::INFO, "Generating Seeds Using Uniform Poisson 3D " + suffix);
	//			sg.generate_seeds();
	//		}
	//		else if (radiusOpt == 1 && distFunc == 0) {
	//			//std::cout << "generating with Poisson 3D and dist from plane" << std::endl;
	//			std::string tag = "Generating Seeds Using Varied Poisson 3D " + suffix + " Distance is measured from plane.";
	//			add_log(LogPriority::INFO, tag);
	//			sg.generate_seeds(planeDistance, linearFunc);
	//		}
	//		else if (radiusOpt == 1 && distFunc == 1) {
	//			std::string tag = "Generating Seeds Using Varied Poisson 3D " + suffix + " Distance is measured from mesh outer face.";
	//			add_log(LogPriority::INFO, tag);
	//			sg.generate_seeds(containerDistance, linearFunc);
	//		}
	//		_update_cameras();
	//		sg.get_seeds(seeds);
	//		if (!seeds.empty()) {
	//			seedNr = seeds.size();
	//		}
	//	}

	//	// Poisson 3d inside a container
	//	if (conOption == 2) {
	//		Poisson3DWall sg(containerMesh, rMin, rMax, neighbors, { 10, 10, 10 }, wallResolution);
	//		if (radiusOpt == 0) {
	//			add_log(LogPriority::INFO, "Generating seeds using Uniform Poisson 3D inside Mesh Container");
	//			sg.generate_seeds();
	//		}
	//		else if (radiusOpt == 1 && distFunc == 0) {
	//			add_log(LogPriority::INFO, "Generating Seeds Using Varied Poisson 3D Inside Mesh Container. Distance is measured from plane.");
	//			sg.generate_seeds(planeDistance, linearFunc);
	//		}
	//		else if (radiusOpt == 1 && distFunc == 1) {
	//			add_log(LogPriority::INFO, "Generating Seeds Using Varied Poisson 3D Inside Mesh Container. Distance is measured from mesh outer face.");
	//			std::cout << "generating with Poisson 3D and dist from mesh" << std::endl;
	//			sg.generate_seeds(meshDistance, linearFunc);
	//		}
	//		_update_cameras();
	//		sg.get_radii(radii);
	//		sg.get_seeds(seeds);
	//		if (!seeds.empty()) {
	//			seedNr = seeds.size();
	//		}
	//	} 
	//}
	//seedGenerator->run_generate_seeds();
	//seedGenerator->run_get_seeds(seeds);

	// create a bounding box

	// create the container if it is a mesh
	//if (conOption == 2) {
	//	//containerModel = new Model(containerFile);
	//	container.model = std::make_unique<Model>(containerFile);
	//}

	//seedObj = new VisualizeSeeds(seeds);
seedObj = std::make_unique<VisualizeSeeds>(seeds);
showSeeds = true;

add_log(LogPriority::SUCCESS, "Seeds are ready.");
};

void myGUI::_action_update_voronoi() {

	if (seeds.empty()) {
		add_log(LogPriority::ERROR, "First Populate Seeds!");
		return;
	}

	//switch (conOption) {
	//default:
	//	domainVolume = (xDim - 2.0f * thickness) * (yDim - 2.0f * thickness) * (zDim - 2.0f * thickness);
	//	break;
	//case 1:
	//	domainVolume = std::pow(cylinderRadius, 2) * 3.14f * cylinderHeight;
	//	break;
	//case 2:
	//	vtkNew<vtkMassProperties> massProperties;
	//	massProperties->SetInputData(containerMesh);
	//	massProperties->Update();
	//	domainVolume = massProperties->GetVolume();
	//	break;
	//}

	//std::cout << "Domain Volume: " << domainVolume << std::endl;

	if (!runVolumeOptimization) {

		Bounds bounds = activeContainer->compute_bounds();

		std::array<float, 6> bnds{
				static_cast<float>(bounds.xMin),
				static_cast<float>(bounds.xMax),
				static_cast<float>(bounds.yMin),
				static_cast<float>(bounds.yMax),
				static_cast<float>(bounds.zMin),
				static_cast<float>(bounds.zMax)
		};

		std::array<int, 3> blockDim{ nX, nY, nZ };

		switch (conOption){
			case 0: {
				auto* box = dynamic_cast<const BoxContainer*>(activeContainer);
				
				domainVolume = (
					box->xDim - 2.0f * thickness) * 
					(box->yDim - 2.0f * thickness) *
					(box->zDim - 2.0f * thickness);

				auto sgb = std::make_unique<ScaffoldGeneratorBox>(
					seeds, bnds, blockDim, thickness, connectivityThreshold, box->is_inside()
				);

				generator = std::move(sgb);
				_update_cameras();

				//add_log(LogPriority::INFO, "Generating scaffold inside rectangular container. Please wait...");
				break;
			}
			case 1: {

				auto* cylinder = dynamic_cast<const CylinderContainer*>(activeContainer);
				
				auto sgb = std::make_unique<ScaffoldGeneratorBox>(
					seeds, bnds, blockDim, thickness, connectivityThreshold, cylinder->is_inside()
				);

				sgb->add_cylindrical_wall(
					cylinder->cylinderPt[0], cylinder->cylinderPt[1], cylinder->cylinderPt[2],
					cylinder->cylinderAxis[0], cylinder->cylinderAxis[1], cylinder->cylinderAxis[2],
					cylinder->cylinderRadius
				);

				domainVolume = std::pow(cylinder->cylinderRadius, 2) * 3.14f * cylinder->cylinderHeight;
				
				generator = std::move(sgb);
				_update_cameras();
				break;
			}
			case 2: {

				auto* mesh = dynamic_cast<const MeshContainer*>(activeContainer);

				vtkSmartPointer<vtkPolyData> container = mesh->containerMesh;

				auto sgb = std::make_unique<ScaffoldGeneratorWall>(
					seeds,
					container,
					blockDim,
					mesh->neighbors,
					wallResolution,
					thickness,
					connectivityThreshold,
					mesh->is_inside()
				);

				vtkNew<vtkMassProperties> massProperties;
				massProperties->SetInputData(mesh->containerMesh);
				massProperties->Update();
				domainVolume = massProperties->GetVolume();

				generator = std::move(sgb);
				_update_cameras();

				break;
			}
		}



		//if (conOption == 0 || conOption == 1) {



		//	auto sgb = std::make_unique<ScaffoldGeneratorBox>(
		//		seeds, bnds, blockDim, thickness, connectivityThreshold, inside_check
		//	);

		//	if (conOption == 0) {

		//		domainVolume = (xDim - 2.0f * thickness) * (yDim - 2.0f * thickness) * (zDim - 2.0f * thickness);

		//		//add_log(LogPriority::INFO, "Generating scaffold inside rectangular container. Please wait...");
		//	}
		//	else if (conOption == 1) {
		//		auto* cyl = dynamic_cast<const CylinderContainer*>(activeContainer);

		//		add_log(LogPriority::INFO, "Generating scaffold inside cylindrical container. Please wait...");

		//		sgb->add_cylindrical_wall(
		//			cyl->cylinderPt[0], cyl->cylinderPt[1], cyl->cylinderPt[2],
		//			cyl->cylinderAxis[0], cyl->cylinderAxis[1], cyl->cylinderAxis[2],
		//			cyl->cylinderRadius
		//		);

		//		domainVolume = std::pow(cyl->cylinderRadius, 2) * 3.14f * cyl->cylinderHeight;
		//	}

		//	generator = std::move(sgb);
		//	_update_cameras();
		//}

	//	else if (conOption == 2) {

	//		add_log(LogPriority::INFO, "Generating scaffold inside custom mesh container. Please wait...");
	//		
	//		auto* mesh = dynamic_cast<const MeshContainer*>(activeContainer);
	//		
	//		vtkSmartPointer<vtkPolyData> container = mesh->containerMesh;

	//		auto sgb = std::make_unique<ScaffoldGeneratorWall>(
	//			seeds,
	//			container,
	//			blockDim,
	//			mesh->neighbors,
	//			wallResolution,
	//			thickness,
	//			connectivityThreshold,
	//			mesh->is_inside()
	//		);

	//		vtkNew<vtkMassProperties> massProperties;
	//		massProperties->SetInputData(mesh->containerMesh);
	//		massProperties->Update();
	//		domainVolume = massProperties->GetVolume();

	//		//ScaffoldGeneratorWall sgb(seeds, containerMesh, { nX, nY, nZ }, neighbors, wallResolution, edgeSize, scaleFactor);
	////		//sgb.generate_voro(regSteps);
	////		//sgb.get_seeds(seeds);
	////		//sgb.generate_mesh(thickness, scaffoldPolyData, { resolution[0], resolution[1], resolution[2] });
	////		//seedObj = std::make_unique<VisualizeSeeds>(seeds);

	//		generator = std::move(sgb);
	//		_update_cameras();
	//	}
	}

	generator->run_generate_voro(regSteps);
	add_log(LogPriority::INFO, "Generating Voronoi with " + std::to_string(regSteps) + " regularization steps");
};

void myGUI::_action_estimate_connectivity() {

	if (!generator) {
		add_log(LogPriority::ERROR, "First Create Seeds");
		return;
	}

	scaffoldConnectivity = generator->run_estimate_connectivity();
	std::vector<float> vertices = generator->run_get_connected_vertices();
	std::vector<unsigned int> indices = generator->run_get_connected_edges();

	poreNetworkObject poreNetwork;
	poreNetwork.model = std::make_unique<PoreNetwork>(vertices, indices);
	poreNetwork.name += std::to_string(poreNetworkList.size());
	poreNetworkList.push_back(std::move(poreNetwork));
	activePoreNetworkIndex = static_cast<int>(poreNetworkList.size() - 1);
	
	add_log(LogPriority::INFO, "Estimated Connectivity " + std::to_string(scaffoldConnectivity));
};

void myGUI::_action_generate_scaffold() {

	add_log(LogPriority::INFO, "Generating Scaffold " + suffix + " Please Wait...");

	std::array<int, 3> blockDim = { resolution[0], resolution[1], resolution[2] };

	if (!runVolumeOptimization) {

		generator->run_get_seeds(seeds);
		generator->run_process_faces();
		generator->run_generate_mesh(thickness, scaffoldPolyData, blockDim);
	}

//	else {
//		add_log(LogPriority::INFO, "Generating with Volume Optimization. Please wait...");
//
//		//std::cout << "Generating with Volume Optimization" << std::endl;
//
//		if (seeds.empty()) {
//
//			add_log(LogPriority::ERROR, "First Populate Seeds!");
//
//			return;
//		}
//		int seedSize = seeds.size();
//
//		wInit.resize(seedSize);
//		wInit.setZero();
//		//std::cout << " -------------------- " << std::endl;
//
//		targetVols.resize(seedSize);
//
//		// case 1: volumes all equal -> divide domain volume per seed nr
//		if (volOption == 0) {
//
//			for (int i = 0; i < seedSize; i++) {
//				targetVols[i] = domainVolume / seedSize;
//			}
//		}
//
//		// case 2: volumes are decided from Poisson 3d radii 
//
//		else if (volOption == 1) {
//
//			double volSum{ 0.0 };
//			for (int i{ 0 }; i < radii.size(); i++) {
//				targetVols[i] = std::pow(radii[i], 3);
//				volSum += targetVols[i];
//			}
//			// Normalize
//			for (int i{ 0 }; i < seeds.size(); i++) {
//				targetVols[i] = (targetVols[i] / volSum) * domainVolume;
//			}
//
//			double targetSum = 0;
//			for (int i{ 0 }; i < seeds.size(); i++) {
//				targetSum += targetVols[i];
//			}
//			//std::cout << "domain volume: " << domainVolume << std::endl;
//			//std::cout << "target volume sum: " << targetSum << std::endl;
//			// std::cerr << "Not implemented yet!" << std::endl;
//		}
//
//		else if (volOption == 2) {
//			add_log(LogPriority::ERROR, "Not implemented yet!");
//			return;
////;			std::cerr << "Not implemented yet!" << std::endl;
//		}
//
//		// just a simple box container
//		if (conOption == 0) {
//			add_log(LogPriority::INFO, "Vol Opt inside Box");
//			//std::cout <<  << std::endl;
//			VolOpt vo(
//				seeds,
//				targetVols,
//				wInit,
//				{ xMin, xMax, yMin, yMax, zMin, zMax }
//			);
//			vo.loop(regSteps);
//			vo.get_seeds(seeds);
//			//vo.generate_mesh(thickness, scaffoldFilePath);
//			vo.generate_mesh(thickness, scaffoldPolyData, { resolution[0], resolution[1], resolution[2] });
//		}
//
//		// volume optimization inside a mesh wall
//		else if (conOption == 2) {
//			add_log(LogPriority::INFO, "Vol Opt inside custom mesh container.");
//
//			//std::cout << "Vol Opt inside Mesh" << std::endl;
//
//			VolOptWall vo(
//				seeds,
//				targetVols,
//				wInit,
//				containerMesh
//			);
//			vo.loop(regSteps);
//			vo.get_seeds(seeds);
//			vo.generate_mesh(thickness, scaffoldFilePath, { resolution[0], resolution[1], resolution[2] });
//		}
//
//	}

	// update the thread status
	//scaffoldGenerating = false;

};

void myGUI::_create_dockspace() {

	ImGui::DockSpaceOverViewport(0, ImGui::GetMainViewport());

	_render_main_menu_bar();

};

void myGUI::_render_visualizer() {

	ImGui::SetNextWindowSize(ImVec2(framebuffer.width, framebuffer.height), ImGuiCond_Once);
	ImGui::Begin("Visualizer", &showVisualizer);
	
	cameraUpdate = (ImGui::IsWindowHovered() && ImGui::IsWindowFocused()) ? true : false;

	ImVec2 availableSize = ImGui::GetContentRegionAvail();
	ImVec2 pos = ImGui::GetCursorScreenPos();

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
		ImVec2 window_pos, window_pos_pivot;
		window_pos.x = framebuffer.posx + 0.82 * framebuffer.width;
		window_pos.y = framebuffer.posy + 0.76 * framebuffer.height;

		window_pos_pivot.x = 1.0f;
		window_pos_pivot.y = 1.0f;
		ImGui::SetNextWindowPos(window_pos, ImGuiCond_Always);
		ImGui::SetNextWindowSize(ImVec2(0.15f * framebuffer.width, 0.23f * framebuffer.height));

		//ImGui::SetNextWindowViewport(viewport->ID);
		window_flags |= ImGuiWindowFlags_NoMove;
		ImGui::SetNextWindowBgAlpha(0.35f); // Transparent background
		if (ImGui::Begin("Mesh Properties", nullptr, window_flags))
		{
			ImGui::Text("Mesh Details");
			ImGui::Separator();
			ImGui::SetNextItemWidth(ImGui::GetFontSize() * 0.8);
			ImGui::Text("Mesh Name: %s", scaffoldFileName.c_str());
			ImGui::SetNextItemWidth(ImGui::GetFontSize() * 0.8);
			ImGui::Text("Mesh Vertices: %d", vertexNr);
			ImGui::SetNextItemWidth(ImGui::GetFontSize() * 0.8);
			ImGui::Text("Mesh Faces: %d", faceNr);
			ImGui::SetNextItemWidth(ImGui::GetFontSize() * 0.8);
			ImGui::Text("Volume: %.3f", scaffoldVolume);
			ImGui::SetNextItemWidth(ImGui::GetFontSize() * 0.8);
			ImGui::Text("Porosity: %.3f", scaffoldPorosity);
			ImGui::SetNextItemWidth(ImGui::GetFontSize() * 0.8);
			ImGui::Text("Interconnectivity: %.3f", scaffoldConnectivity);
			ImGui::End();
		}
	}


	ImGui::End();

};

void myGUI::_render_binary_image_window(const char* popupName, bool& showPopup) {
	
	if (showPopup) {
		ImGui::OpenPopup(popupName);
	}
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
	ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_Appearing);
	if (ImGui::BeginPopup(popupName, ImGuiWindowFlags_AlwaysAutoResize)) {

		static float voxelSize{ 0.05f };
		static float backVal = { 0 };
		static float forVal = { 1 };
		static std::string name = "";

		ImGui::InputFloat("Voxel Size", &voxelSize, 0.01f, 10.0f);
		ImGui::InputFloat("Background Value", &backVal, 0.01f, 10.0f);
		ImGui::InputFloat("Foreground", &forVal, 0.01f, 10.0f);

		backVal = std::clamp(backVal, 0.0f, 1.0f);
		forVal = std::clamp(forVal, 0.0f, 1.0f);

		if (ImGui::Button("Save")) {
			IGFD::FileDialogConfig config;
			config.path = "..//data";
			config.flags = ImGuiFileDialogFlags_Modal;
			ImGuiFileDialog::Instance()->OpenDialog("Export Binary Scaffold", "Export Binary", ".mhd", config);
		}
		ImGui::SameLine();
		if (ImGui::Button("Close")) {
			showPopup = false;
			ImGui::CloseCurrentPopup();
		}

		ImVec2 maxSize = ImVec2(width, height);  
		ImVec2 minSize = ImVec2(500.0f, 500.0f); 
		if (ImGuiFileDialog::Instance()->Display("Export Binary Scaffold", ImGuiWindowFlags_NoCollapse, minSize, maxSize)) {
			if (ImGuiFileDialog::Instance()->IsOk()) { // action if OK
				scaffoldFilePath = ImGuiFileDialog::Instance()->GetFilePathName();
				scaffoldFileName = ImGuiFileDialog::Instance()->GetCurrentFileName();
				_export_binary_image(scaffoldFilePath,
					 voxelSize, backVal, forVal);
			}
			ImGuiFileDialog::Instance()->Close();
		}

		ImGui::EndPopup();
	}

};
