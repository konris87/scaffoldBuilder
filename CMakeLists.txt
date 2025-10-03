cmake_minimum_required(VERSION 3.10)

# Set the project name
project(scaffoldBuilder)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 17)

# --------------------------
# VTK
find_package(VTK COMPONENTS 
  CommonColor
  CommonCore
  CommonDataModel
  CommonTransforms
  FiltersGeneral
  FiltersSources
  FiltersCore
  FiltersGeometry
  FiltersHybrid
  IOGeometry
  IOImage
  IOCore
  IOLegacy
  IOXML
  InteractionStyle
  RenderingAnnotation
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
  RenderingLabel
  ImagingStencil
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "${sName}: Unable to find the VTK build folder.")
endif()

# Find OpenGL
find_package(OpenGL REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS} )

# Add GLFW library
# set (CMAKE_PREFIX_PATH "C:/Users/Konst/Desktop/Software/glfw_install/")
# find_package(glfw3 3.3 REQUIRED)
set(GLFW_DIR ./lib/glfw)
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

# Include Voro ++ 
set (CMAKE_PREFIX_PATH "C:/Users/Konst/Desktop/Software/voro-install/bin/")
find_package(VORO CONFIG REQUIRED)

if (NOT VORO_FOUND)
  message(FATAL_ERROR "Unable to find VORO.")
endif()

set(SOURCES

	src/Gui/guiApp.cpp
    src/main.cpp
	src/Scaffold/buildScaffold.cpp
	
	# glad
	lib/glad/src/glad.c
	
	# other src
	src/Misc/cell_process.cpp
	include/Math/Quaternion.cpp
	include/Math/Vec.cpp
	src/Visualize/visualize.cpp
	src/Visualize/loadMesh.cpp
	src/OpenGlSetup/simpleCamera.cpp
	src/OpenGlSetup/trackBallCamera.cpp
	src/OpenGlSetup/defaultCamera.cpp
	src/OpenGlSetup/Misc.cpp
	include/Utils/Utils.cpp
	include/Wall/Wall.cpp
	include/SeedGenerator/Random.cpp
	include/SeedGenerator/Poisson3D.cpp
	include/ScaffoldGenerator/ScaffoldGenerator.cpp
	# include/ScaffoldGenerator/Poisson3D.cpp
	
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

# Add executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
    ${OPENGL_LIBRARIES}
	glfw
	${VTK_LIBRARIES}
)

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS ${PROJECT_NAME}
  MODULES ${VTK_LIBRARIES}
)

target_include_directories(${PROJECT_NAME} PUBLIC VORO::voro++)
target_link_libraries(${PROJECT_NAME} PUBLIC VORO::voro++)

# Eigen
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# SOURCE_GROUP(common REGULAR_EXPRESSION ".*/common/.*" )
# SOURCE_GROUP(shaders REGULAR_EXPRESSION ".*/.*shader$" )
# SOURCE_GROUP(obj REGULAR_EXPRESSION ".*/.*obj$" )

# Set the shader source directory
# set(SHADER_SRC_DIR "${CMAKE_SOURCE_DIR}/src/shaders")

# Set the directory where the executable will be located
# set(SHADER_DEST_DIR "${CMAKE_BINARY_DIR}/shaders")

# also copy the shaders and fonts to the build file during development
file(COPY ${CMAKE_SOURCE_DIR}/src/shaders DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/share DESTINATION ${CMAKE_BINARY_DIR})

# -----------------------------
# Installation rules
# -----------------------------

# Install the main executable
install(TARGETS scaffoldBuilder
    RUNTIME DESTINATION bin           # Windows/Linux executable
    BUNDLE DESTINATION .   
	LIBRARY DESTINATION	bin# macOS bundle
)

# Install shader files (uncomment if needed)
install(DIRECTORY src/shaders/ DESTINATION bin/shaders)
install(DIRECTORY share/ DESTINATION bin/share)
install(DIRECTORY lib/vtk/ DESTINATION bin)

# -----------------------------
# CPack setup (packaging)
# -----------------------------
include(InstallRequiredSystemLibraries)

set(CPACK_PACKAGE_NAME "scaffoldBuilder")
set(CPACK_PACKAGE_VENDOR "VVR")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Scaffold Builder with VTK, ImGui, GLFW, Voro++")
set(CPACK_PACKAGE_VERSION "0.5")
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

# -----------------------------
# Runtime library search paths
# -----------------------------
if(WIN32)
    # Windows: we’ll use SetDllDirectory in code, so no RPATH
elseif(APPLE)
    # macOS: search in ./vtk next to executable
    set(CMAKE_INSTALL_RPATH "@executable_path/vtk")
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
else()
    # Linux: search in ./vtk next to executable
    set(CMAKE_INSTALL_RPATH "$ORIGIN/vtk")
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif()