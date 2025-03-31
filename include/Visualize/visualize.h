#ifndef VISUALIZE_H // include guard
#define VISUALIZE_H

#include <vector>
#include <string>
#include <array>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

vtkSmartPointer<vtkPolyData> cell_2_vtk(std::vector<int>& cellNeighs, const std::vector<double>& cellVertices, std::vector<int>& faceVertices);

void create_mesh(
	const std::vector<vtkSmartPointer<vtkPolyData>>& polys,
	const float& thickness, const std::array<double, 6>& bounds,
	const std::string& fileName
);

#endif