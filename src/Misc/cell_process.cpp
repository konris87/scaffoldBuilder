#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Eigen>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkCellIterator.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDelaunay3D.h>
#include <Misc/cell_process.h>

double tet_integral(const std::vector<double>& tetVerts){
	/* tet_verts is a vector of 9 values, every three are the coords of a vertex, the apex is at zero
	 */
	
	double x1 = tetVerts[0], y1 = tetVerts[1], z1 = tetVerts[2]; 
	double x2 = tetVerts[3], y2 = tetVerts[4], z2 = tetVerts[5]; 
	double x3 = tetVerts[6], y3 = tetVerts[7], z3 = tetVerts[8]; 
	
	return pow(x1, 2) + pow(x2, 2) + pow(x3, 2) + x1*x2 + x1*x3 + x2*x3 +\
      pow(y1, 2) + pow(y2, 2) + pow(y3, 2) + y1*y2 + y1*y3 + y2*y3 +\
      pow(z1, 2) + pow(z2, 2) + pow(z3, 2) + z1*z2 + z1*z3 + z2*z3;
}

/*double transport(std::vector<double>& cell_verts, std::vector<double>& cell_center){

	double cell_cost{0};
	
	for(int idx=0; idx < cell_verts.size(); idx++){

		std::vector<double> tet_verts = {
			cell_verts[3*idx] - cell_center[0],
			cell_verts[3*idx + 1] - cell_center[1],
			cell_verts[3*idx + 2] - cell_center[2]} ; 
		
		cell_cost += tet_integral(tet_verts);
	}

	std::cout << "the cell cost is: " << cell_cost << std::endl;
	return cell_cost;
}*/

double transport(const std::vector<double>& cellVertices, const std::vector<double>& cellCenter) {
	

	double cellCost{0};

	int Nv = cellVertices.size() / 3;

	// create a vtk points obj
	vtkSmartPointer<vtkPoints> polyPts = vtkSmartPointer<vtkPoints>::New();

	// add each vertex
	for (int i = 0; i < Nv; i++) {
		polyPts->InsertNextPoint(
			cellVertices[3 * i],
			cellVertices[3 * i + 1],
			cellVertices[3 * i + 2]);
	}

	// create the poly
	vtkNew<vtkPolyData> pdata;
	pdata->SetPoints(polyPts);

	vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
	delaunay->SetInputData(pdata);
	delaunay->Update();

	// surface filter
	vtkNew<vtkDataSetSurfaceFilter> surfaceFilter;
	surfaceFilter->SetInputConnection(delaunay->GetOutputPort());
	surfaceFilter->Update();

	// get the polydata and cell data of triangulation
	vtkPolyData* polydata = surfaceFilter->GetOutput();


	vtkCellIterator* it = polydata->NewCellIterator();
	for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell())
	{
		vtkIdList* pointIds = it->GetPointIds();

		vtkPoints* points = it->GetPoints();

		std::vector<double> cellpts = {};

		for (vtkIdType j = 0; j < points->GetNumberOfPoints(); j++) {
			double pt[3];
			points->GetPoint(j, pt);
			cellpts.push_back(pt[0] - cellCenter[0]);
			cellpts.push_back(pt[1] - cellCenter[1]);
			cellpts.push_back(pt[2] - cellCenter[2]);
			/*std:cout << pt[0] << ", " << pt[1] << ", " << pt[2] << std::endl;*/
		}

		/*std::cout << cellpts.size() << std::endl;
		std::cout << "The volume is: " << volume << std::endl;*/

		cellCost += (tet_volume(cellpts) * tet_integral(cellpts) / 60);

	}
	it->Delete();


	return cellCost;
}


double tet_volume(const std::vector<double>& tetVerts){
	Eigen::Matrix<double, 3, 3> X;
    X(0, 0) = tetVerts[0];
    X(1, 0) = tetVerts[1];
    X(2, 0) = tetVerts[2];

    X(0, 1) = tetVerts[3];
    X(1, 1) = tetVerts[4];
	X(2, 1) = tetVerts[5];

	X(0, 2) = tetVerts[6];
    X(1, 2) = tetVerts[7];
    X(2, 2) = tetVerts[8];
    return X.determinant();

}

double euclidean_distance(std::array<double, 3>& pts1, std::array<double, 3>& pts2) {

	double distance{ 0 };

	distance = sqrt(
		pow(pts2[0] - pts1[0], 2) +
		pow(pts2[1] - pts1[1], 2) +
		pow(pts2[2] - pts1[2], 2) 
	); 

	return distance;
}


std::array<double, 3> compute_centroid(const std::vector<double>& vertices) {
    std::array<double, 3> centroid = { 0.0, 0.0, 0.0 };
    int Nv = vertices.size() / 3;  // Number of vertices
    for (int i = 0; i < Nv; ++i) {
        centroid[0] += vertices[3 * i];
        centroid[1] += vertices[3 * i + 1];
        centroid[2] += vertices[3 * i + 2];
    }
    centroid[0] /= Nv;
    centroid[1] /= Nv;
    centroid[2] /= Nv;
    return centroid;
}

Eigen::VectorXd convert_radii(const Eigen::VectorXd& weights) {

	// first find the minimum
	double wMin = weights.minCoeff();

	Eigen::VectorXd radii = Eigen::VectorXd::Zero(weights.size());

	for (int i = 0; i < weights.size(); i++) {
		radii(i) = sqrt(weights(i) - wMin);
	}
	return radii;

}