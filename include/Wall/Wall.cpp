#include "Wall.h"
#include "Utils/Utils.h"
#include <vtkTriangle.h>

MeshWall::MeshWall(
	const vtkSmartPointer<vtkPolyData>& meshData,
	const std::vector<std::array<double, 3>>& barycenters,
	const std::vector<std::array<double, 3>>& faceNormals,
	const std::vector<std::array<double, 3>>& trianglePt,
	vtkSmartPointer<vtkKdTreePointLocator> kdTree,
	const int neighbors,
	int iw_id)
: barycenters(barycenters), normals(faceNormals), trianglePt(trianglePt),
	kdTree(kdTree), neighbors(neighbors), w_id(iw_id) {
	
	mesh = vtkSmartPointer<vtkPolyData>::New();
};

MeshWall::MeshWall(
	const vtkSmartPointer<vtkPolyData>& meshData,
	const std::vector<std::array<double, 3>>& vertexNormals,
	vtkSmartPointer<vtkKdTreePointLocator> kdTree,
	const int neighbors,
	int iw_id)
	: vertexNormals(vertexNormals), kdTree(kdTree), neighbors(neighbors), w_id(iw_id) {

	mesh = vtkSmartPointer<vtkPolyData>::New();
};

MeshWall::MeshWall(
	const vtkSmartPointer<vtkPolyData>& meshData,
	vtkSmartPointer<vtkKdTreePointLocator> kdTree,
	const int neighbors,
	int iw_id) : neighbors(neighbors), w_id(iw_id) {

	mesh = vtkSmartPointer<vtkPolyData>::New();

	for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); ++i) {
		vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(mesh->GetCell(i));
		if (!triangle) continue;

		double p1[3], p2[3], p3[3];
		triangle->GetPoints()->GetPoint(0, p1);
		triangle->GetPoints()->GetPoint(1, p2);
		triangle->GetPoints()->GetPoint(2, p3);

		std::array<double, 3> v1 = { p1[0], p1[1], p1[2] };
		trianglePt.push_back(v1);

		// convert pts to Eigen
		Eigen::Vector3d pe1{ p1[0], p1[1], p1[2] };
		Eigen::Vector3d pe2{ p2[0], p2[1], p2[2] };
		Eigen::Vector3d pe3{ p3[0], p3[1], p3[2] };

		// estimate the barycenter
		double bc1 = (p1[0] + p2[0] + p3[0]) / 3;
		double bc2 = (p1[1] + p2[1] + p3[1]) / 3;
		double bc3 = (p1[2] + p2[2] + p3[2]) / 3;
		std::array<double, 3> bc = { bc1, bc2, bc3 };
		//Eigen::Vector3d bc{ bc1, bc2, bc3 };
		barycenters.push_back(bc);

		// estimate the normal
		Eigen::Vector3d normal = (pe3 - pe1).cross(pe2 - pe1);
		normal.normalize();
		normals.push_back(bc);
	}
};

/**
* Check if a point is inside the Mesh. A random ray is casted and the number 
* of intersections with the triangles of the mesh is estimated
* If count == odd then the point is inside
* If count == even then the point is outside
* 
* It starts by checking if the point is inside the bounds of the mesh to avoid
* reduntant calculations
* 
* @param[point] EigenVector3d that contains the coords of the point
* @ returns true if point is inside the mesh
* 
* It uses the algorithm: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
* 
*/
bool MeshWall::isInside(Eigen::Vector3d& point){

	Eigen::Vector3d rayDir(1.0, 0.0, 0.0);
	int count{ 0 };

	// get the triangles
	for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); ++i) {
		vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(mesh->GetCell(i));
		if (!triangle) continue;

		double p1[3], p2[3], p3[3];
		triangle->GetPoints()->GetPoint(0, p1);
		triangle->GetPoints()->GetPoint(1, p2);
		triangle->GetPoints()->GetPoint(2, p3);

		// store the points in Eigen vectors
		Eigen::Vector3d v1(p1[0], p1[1], p1[2]);
		Eigen::Vector3d v2(p2[0], p2[1], p2[2]);
		Eigen::Vector3d v3(p3[0], p3[1], p3[2]);

		// use the algorithm
		if (ray_intersection(point, v1, v2, v3, rayDir)) {
			count++;
		}

	}

	if (count % 2 == 1) { return true; }
	else { return false; }
};

bool MeshWall::point_inside(double x, double y, double z) {
	Eigen::Vector3d p(x, y, z);

	return isInside(p);
}