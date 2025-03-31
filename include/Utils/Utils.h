#ifndef UTILS_H // include guard
#define UTILS_H

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include "Eigen/Dense"
#include <string.h>
#include <vtkGlyph3D.h>
#include <array>
#include <vector>
#include <ImGuiFileDialog/ImGuiFileDialog.h>

bool ray_intersection(
	const Eigen::Vector3d& p,
	const Eigen::Vector3d& v1,
	const Eigen::Vector3d& v2,
	const Eigen::Vector3d& v3,
	const Eigen::Vector3d& dir
);

bool ray_intersection(
	const Eigen::Vector3d& p,
	const Eigen::Vector3d& v1,
	const Eigen::Vector3d& v2,
	const Eigen::Vector3d& v3,
	const Eigen::Vector3d& dir,
	Eigen::Vector3d& intersection
);

bool is_inside_mesh(
	const vtkSmartPointer<vtkPolyData>& mesh,
	const Eigen::Vector3d& point);

bool is_inside_mesh(
	const vtkSmartPointer<vtkPolyData>& mesh,
	const Eigen::Vector3d& point,
	const Eigen::Vector3d& rayDir,
	Eigen::Vector3d& intersection);

bool is_inside_box(
	const std::array<double, 6> bounds);

vtkSmartPointer<vtkActor> render_points(
	const int& particles,
	const std::vector<std::array<double, 3>>& currSeeds, 
	const std::string color);

void make_glyphs(vtkPolyData* src, double size, vtkGlyph3D* glyph);

vtkSmartPointer<vtkActor> render_cell(const std::vector<double>& cellVerts);

vtkSmartPointer<vtkActor> render_plane(
	const double& n1, const double& n2, const double& n3,
	const double& p1, const double& p2, const double& p3,
	const double& s1, const double& s2
);

double distance(const std::array<double, 3>& arr1, const std::array<double, 3>& arr2);

void selectFileButton(const char* title, const std::string path, const char* dialogName, const char* fileExt);
#endif