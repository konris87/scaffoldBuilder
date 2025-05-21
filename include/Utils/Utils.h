#ifndef UTILS_H // include guard
#define UTILS_H

# define _USE_MATH_DEFINES

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include "Eigen/Dense"
#include <string.h>
#include <vtkGlyph3D.h>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
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

void ensure_ccw(
	std::vector<int>& faceIdxs,
	Eigen::MatrixXd& vertices);

void project_vertices_on_plane(
	const std::vector<double>& vertices,
	Eigen::Vector3d& u,
	Eigen::Vector3d& v,
	Eigen::MatrixXd& vertices2D
);


// create an argsort function
// taken by https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes

template <typename T>
std::vector<size_t> sort_indices(const std::vector<T>& v) {

	// initialize original index locations
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	// using std::stable_sort instead of std::sort
	// to avoid unnecessary index re-orderings
	// when v contains elements of equal values 
	std::stable_sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx;
};

void sample_face_polygon(
	const Eigen::MatrixXd& verts2D,
	Eigen::Vector2d& center, double& radius1, int resolution);

bool is_inside_polygon(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices);

double distance_from_edges(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices);

void interpolate_edges(const Eigen::MatrixXd& vertices, Eigen::MatrixXd& interpolatedVertices, const double& edgeSize);

void hole_points(const double& radius, const Eigen::Vector2d& center, const int& ptNr, Eigen::MatrixXd& vertices);
