cmake_minimum_required(VERSION 3.10)

# Set the project name
project(scaffoldBuilder)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 17)

# --------------------------
# VTK
# find_package(VTK COMPONENTS 
  # AcceleratorsVTKmCore
  # AcceleratorsVTKmFilters
  # CommonColor
  # CommonCore
  # CommonDataModel
  # CommonTransforms
  # FiltersGeneral
  # FiltersSources
  # FiltersCore
  # FiltersSources
  # FiltersParallel
  # FiltersGeometry
  # FiltersHybrid
  # FiltersModeling
  # FiltersPoints
  # IOGeometry
  # IOImage
  # IOCore
  # IOLegacy
  # IOXML
  # InteractionStyle
  # RenderingAnnotation
  # RenderingContextOpenGL2
  # RenderingCore
  # RenderingFreeType
  # RenderingGL2PSOpenGL2
  # RenderingOpenGL2
  # RenderingLabel
  # ImagingStencil
# )

# if (NOT VTK_FOUND)
  # message(FATAL_ERROR "${sName}: Unable to find the VTK build folder.")
# endif()

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

# Include Voro ++ 
set (CMAKE_PREFIX_PATH "C:/Users/Konst/Desktop/Software/voro-install/bin/")
find_package(VORO CONFIG REQUIRED)

if (NOT VORO_FOUND)
  message(FATAL_ERROR "Unable to find VORO.")
endif()

set(SOURCES

	src/Gui/guiApp.cpp
    src/main.cpp
	# src/Scaffold/buildScaffold.cpp
	
	# glad
	lib/glad/src/glad.c
	
	# other src
	# src/Misc/cell_process.cpp
	include/Math/Quaternion.cpp
	include/Math/Vec.cpp
	include/Math/AabbTree.cpp
	#src/Visualize/visualize.cpp
	#src/Visualize/loadMesh.cpp
	src/OpenGlSetup/simpleCamera.cpp
	src/OpenGlSetup/trackBallCamera.cpp
	src/OpenGlSetup/defaultCamera.cpp
	src/OpenGlSetup/Misc.cpp
	include/Utils/Utils.cpp
	include/Wall/Wall.cpp
	include/Misc/Imgui_Stdlib.cpp
	# include/SeedGenerator/Random.cpp
	# include/SeedGenerator/Poisson3D.cpp
	# include/ScaffoldGenerator/ScaffoldGenerator.cpp
	include/SeedGenerator/Container.cpp
	include/SeedGenerator/DistanceCalculator.cpp
	include/ScaffoldGenerator/Generator.cpp
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

# Add executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
	${OPENGL_LIBRARIES}
	glfw
	${VTK_LIBRARIES}
)

# vtk_module_autoinit is needed
#vtk_module_autoinit(
#  TARGETS ${PROJECT_NAME}
#  MODULES ${VTK_LIBRARIES}
#)

target_include_directories(${PROJECT_NAME} PUBLIC VORO::voro++)
target_link_libraries(${PROJECT_NAME} PUBLIC VORO::voro++)

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Eigen
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# also copy the shaders and fonts to the build file during development
add_custom_command(TARGET scaffoldBuilder POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/share
        $<TARGET_FILE_DIR:scaffoldBuilder>/share
    COMMENT "Copying assets to output directory..."
)

# Repeat for shaders if they are separate
add_custom_command(TARGET scaffoldBuilder POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/src/shaders
        $<TARGET_FILE_DIR:scaffoldBuilder>/shaders
)

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
# install(DIRECTORY src/shaders/ DESTINATION bin/shaders)
install(DIRECTORY share/ DESTINATION bin/share)

# vtk dlls
# if(WIN32)
    # # Strategy 1: Try finding the Release DLL specifically
    # get_target_property(VTK_DLL_CORE VTK::CommonCore IMPORTED_LOCATION_RELEASE)

    # # Strategy 2: If that failed, try the generic location
    # if(NOT VTK_DLL_CORE)
        # get_target_property(VTK_DLL_CORE VTK::CommonCore IMPORTED_LOCATION)
    # endif()

    # # Strategy 3: Fallback - Guess based on where VTKConfig.cmake was found
    # if(NOT VTK_DLL_CORE)
        # # VTK_DIR usually points to .../lib/cmake/vtk-x.x
        # # We assume the bin folder is at .../bin
        # get_filename_component(VTK_ROOT_GUESS "${VTK_DIR}/../../.." ABSOLUTE)
        # set(VTK_RUNTIME_DIR "${VTK_ROOT_GUESS}/bin")
        # message(STATUS "Strategy 3 (Guess) - VTK Binaries might be at: ${VTK_RUNTIME_DIR}")
    # else()
        # # If Strategy 1 or 2 worked, get the directory from the file path
        # get_filename_component(VTK_RUNTIME_DIR "${VTK_DLL_CORE}" DIRECTORY)
    # endif()

    # # --------------------------------------------------
    # # DEBUG: Check the output
    # message(STATUS "--------------------------------------------------")
    # message(STATUS "FINAL DETECTED VTK BIN FOLDER: ${VTK_RUNTIME_DIR}")
    # message(STATUS "--------------------------------------------------")
    # # --------------------------------------------------

    # # SAFETY: Ensure we found a valid path and it is NOT your source folder
    # if(EXISTS "${VTK_RUNTIME_DIR}" AND NOT "${VTK_RUNTIME_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")
        
        # install(DIRECTORY "${VTK_RUNTIME_DIR}/" 
                # DESTINATION bin 
                # FILES_MATCHING PATTERN "*.dll")
    # else()
        # message(WARNING "Could not find valid VTK DLLs to install. Check the path above.")
    # endif()
# endif()

# -----------------------------
# CPack setup (packaging)
# -----------------------------
include(InstallRequiredSystemLibraries)

set(CPACK_PACKAGE_NAME "scaffoldBuilder")
set(CPACK_PACKAGE_VENDOR "VVR")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Scaffold Builder with ImGui, GLFW, Voro++")
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

