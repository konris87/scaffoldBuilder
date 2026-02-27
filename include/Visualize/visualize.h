#ifndef VISUALIZE_H // include guard
#define VISUALIZE_H

#include <vector>
#include <string>
#include <array>
#include <Eigen/Dense>
//#include <vtkPolyData.h>

//vtkSmartPointer<vtkPolyData> cell_2_vtk(std::vector<int>& cellNeighs, const std::vector<double>& cellVertices, std::vector<int>& faceVertices);
//
//void create_mesh(
//	const std::vector<vtkSmartPointer<vtkPolyData>>& polys,
//	const float& thickness, const std::array<double, 6>& bounds,
//	const std::string& fileName
//);

void render_vtk_points(const Eigen::MatrixXd& vertices, const std::string& name);

void render_vtk_face(const Eigen::MatrixXd& vertices, const std::vector<std::vector<int>>& indices, const std::string& name);

//void render_vtk_triangular_cell(const std::vector<vtkSmartPointer<vtkPolyData>>& polys);

//vtkSmartPointer<vtkPolyData> create_face_poly(const Eigen::MatrixXd& vertices, const std::vector<std::vector<int>>& indices);

//void render_vtk_polydata(vtkSmartPointer<vtkPolyData>& polyData);

//void render_vtk_mesh(const Eigen::MatrixXd& vertices, const std::vector<std::vector<int>>& indices, const std::string& name);

#endif