cmake_minimum_required(VERSION 3.10)

# Set the project name
project(scaffoldBuilder)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 17)

# Find OpenGL
find_package(OpenGL REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS} )

# Add GLFW library
# set (CMAKE_PREFIX_PATH "C:/Users/Konst/Desktop/Software/glfw_install/")
# find_package(glfw3 3.3 REQUIRED)
set(GLFW_DIR ./lib/glfw)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
add_subdirectory(${GLFW_DIR})
include_directories(${GLFW_DIR}/include)

# add glad
include_directories(${PROJECT_NAME} "./lib/glad/include/")
include_directories(${ProjectName} PRIVATE ${CMAKE_SOURCE_DIR}/include)
include_directories(
	${ProjectName} PRIVATE 
	${CMAKE_SOURCE_DIR}/lib
	${CMAKE_SOURCE_DIR}/lib/glm-1.0.1/
)

# Include imgui
set(IMGUI_DIR ./lib/imgui)
include_directories(${IMGUI_DIR} ${IMGUI_DIR}/backends ..)

# sources shared by the main app and the CLI profiling executable
set(CORE_SOURCES

	# src/Scaffold/buildScaffold.cpp

	# glad
	lib/glad/src/glad.c

	# other src
	# src/Misc/cell_process.cpp
	include/Math/Quaternion.cpp
	include/Math/Vec.cpp
	include/Math/AabbTree.cpp
	include/Math/Matrix.cpp
	include/Math/Bfgs.cpp
	#src/Visualize/visualize.cpp
	#src/Visualize/loadMesh.cpp
	src/OpenGlSetup/simpleCamera.cpp
	src/OpenGlSetup/trackBallCamera.cpp
	src/OpenGlSetup/defaultCamera.cpp
	src/OpenGlSetup/Misc.cpp
	include/Utils/Utils.cpp
	# include/Wall/Wall.cpp
	include/Misc/Imgui_Stdlib.cpp
	# include/SeedGenerator/Random.cpp
	# include/SeedGenerator/Poisson3D.cpp
	# include/ScaffoldGenerator/ScaffoldGenerator.cpp
	include/SeedGenerator/Container.cpp
	include/SeedGenerator/DistanceCalculator.cpp
	# include/ScaffoldGenerator/Generator.cpp
	include/ScaffoldGenerator/GeneratorLewiner.cpp
	include/OpenGlRender/Model.cpp
	include/SeedGenerator/SeedGenerator.cpp
	include/MonteSimulation/MonteSimulation.cpp

	# imgui
	${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
	${IMGUI_DIR}/backends/imgui_impl_opengl3
	${IMGUI_DIR}/imgui.cpp
	${IMGUI_DIR}/imgui_draw.cpp
 	${IMGUI_DIR}/imgui_widgets.cpp
	${IMGUI_DIR}/imgui_tables.cpp
	${IMGUI_DIR}/imgui_demo.cpp

	# imgui file dialog
	lib/ImGuiFileDialog/ImGuiFileDialog.cpp
)

set(SOURCES
	${CORE_SOURCES}
	src/Gui/guiApp.cpp
	src/main.cpp
)

# Add executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Standalone CLI profiling harness (no GUI/OpenGL window needed) used to
# measure where time goes in the scaffold generation pipeline on a real mesh.
add_executable(scaffoldProfile ${CORE_SOURCES} src/profile_main.cpp)

# settings shared by both the main app and the CLI profiling executable
foreach(TGT ${PROJECT_NAME} scaffoldProfile)

	if (MSVC)
		target_compile_options(${TGT} PRIVATE /bigobj)
	endif()

	# The codebase uses #pragma omp parallel for extensively (scalar field
	# generation, marching cubes, smoothing, etc.) but without this, MSVC/GCC/Clang
	# silently ignore those pragmas and the loops run single-threaded.
	find_package(OpenMP REQUIRED)
	target_link_libraries(${TGT} PRIVATE OpenMP::OpenMP_CXX)

	# Link libraries
	target_link_libraries(${TGT} PRIVATE
		${OPENGL_LIBRARIES}
		glfw
	)

	# Include directories
	target_include_directories(${TGT} PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
	)

	# Eigen
	target_include_directories(${TGT} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

endforeach()

# Copy the share/ folder (shaders, fonts, textures) next to the built executable
# so relative paths like "./share/shaders/..." resolve regardless of build type
# (Debug/Release/RelWithDebInfo) or generator (Ninja/Visual Studio).
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
		${CMAKE_SOURCE_DIR}/share
		$<TARGET_FILE_DIR:${PROJECT_NAME}>/share
)

# vtk_module_autoinit is needed
#vtk_module_autoinit(
#  TARGETS ${PROJECT_NAME}
#  MODULES ${VTK_LIBRARIES}
#)

# -----------------------------
# Installation rules
# -----------------------------

# Install the main executable
install(TARGETS scaffoldBuilder
    RUNTIME DESTINATION bin           # Windows/Linux executable
    BUNDLE DESTINATION .   
	LIBRARY DESTINATION	bin# macOS bundle
)

# Install the share/ folder (shaders, fonts, textures) alongside the executable
install(DIRECTORY ${CMAKE_SOURCE_DIR}/share/ DESTINATION bin/share)
# -----------------------------
# CPack setup (packaging)
# -----------------------------
include(InstallRequiredSystemLibraries)

set(CPACK_PACKAGE_NAME "scaffoldBuilder")
set(CPACK_PACKAGE_VENDOR "VVR")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Trabecular bone scaffold generator (ImGui + GLFW + OpenGL)")
set(CPACK_PACKAGE_VERSION "0.8")
set(CPACK_PACKAGE_CONTACT "krisvas@ece.upatras.gr")

# Choose package formats depending on platform
if(WIN32)
    set(CPACK_GENERATOR "ZIP")     # Windows users like ZIP
elseif(APPLE)
    set(CPACK_GENERATOR "DragNDrop") # macOS .dmg
else()
    set(CPACK_GENERATOR "TGZ")     # Linux tar.gz
endif()

include(CPack)

