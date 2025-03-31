#include <Visualize/visualize.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkLine.h>
#include <vtkCellIterator.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkNamedColors.h>
#include <vtkAppendFilter.h>
#include <vtkImplicitModeller.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkStlWriter.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkPolyDataNormals.h>
#include <vtkReverseSense.h>


vtkSmartPointer<vtkPolyData> cell_2_vtk(std::vector<int>& cellNeighs, const std::vector<double>& cellVertices, std::vector<int>& faceVertices){

	vtkNew<vtkPolyData> cellPolyData;
	vtkNew<vtkCellArray> cellArray;
	vtkNew<vtkPoints> cellPoints;

	int Nv = cellVertices.size() / 3;

	for (int i = 0; i < Nv; i++) {
		cellPoints->InsertNextPoint(
			cellVertices[3 * i],
			cellVertices[3 * i + 1],
			cellVertices[3 * i + 2]);
	}

	int idx{ 0 };
	for (int i = 0; i < cellNeighs.size(); i++) {

		// check if the neihbor id is smaller than this cell id to avoid 
		// double checking

		int order = faceVertices[idx];

		std::vector<int> vertList;
		for (int j = idx + 1; j < order + idx + 1; j++) {
			vertList.push_back(faceVertices[j]);
		}
			
		for (int k = 0; k < vertList.size() - 1; k++) {

			vtkNew<vtkLine> line;
			line->GetPointIds()->SetId(0, vertList[k]);
			line->GetPointIds()->SetId(1, vertList[k + 1]);
			cellArray->InsertNextCell(line);
		}

		//vtkNew<vtkLine> line;
		//line->GetPointIds()->SetId(0, vertList.size() - 1);
		//line->GetPointIds()->SetId(1, vertList[0]);
		//cellArray->InsertNextCell(line);
		
		// update location of idx
		idx += order + 1;
	}
	cellPolyData->SetPoints(cellPoints);
	cellPolyData->SetLines(cellArray);

	return cellPolyData;	
}

// this function creates the mesh after the optimization problem is solved
void create_mesh(
	const std::vector<vtkSmartPointer<vtkPolyData>>& polys,
	const float& thickness, const std::array<double, 6>& bounds,
	const std::string& fileName
	) {

	vtkNew<vtkNamedColors> colors;

	// merge polys
	std::cout << "Merging Polys and Creating Implicit Meshing ... " << std::endl;

	vtkNew<vtkAppendFilter> appendFilter;
	for (int i = 0; i < polys.size(); i++) {
		appendFilter->AddInputData(polys[i]);
	}
	appendFilter->Update();

	int dim[3] = {};

	if (thickness < 0.3) {
		dim[0] = 300;
		dim[1] = 300;
		dim[2] = 300;
	}
	else if (0.3 <= thickness && thickness < 0.5) {
		dim[0] = 100;
		dim[1] = 100;
		dim[2] = 100;
	}
	else {
		dim[0] = 100;
		dim[1] = 100;
		dim[2] = 100;
	}
	std::cout << dim[0] << std::endl;
	// build along line
	vtkNew<vtkImplicitModeller> implictModeller;
	implictModeller->SetInputConnection(appendFilter->GetOutputPort());
	implictModeller->SetMaximumDistance(thickness / 2 * 1.01);
	implictModeller->SetSampleDimensions(dim[0], dim[1], dim[2]);
	implictModeller->SetModelBounds(
		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

	// extract isosurface
	vtkNew<vtkContourFilter> isoFilter;
	isoFilter->SetInputConnection(implictModeller->GetOutputPort());
	isoFilter->SetValue(0, thickness / 2);
	isoFilter->Update();

	// update normals
	vtkNew<vtkPolyDataNormals> norms;
	//vtkPolyDataNormals* norms = vtkPolyDataNormals::New();

	norms->SetInputConnection(isoFilter->GetOutputPort());
	norms->ComputePointNormalsOff();
	norms->ComputeCellNormalsOn();
	norms->ConsistencyOn();
	norms->AutoOrientNormalsOn();
	norms->Update();

	//// flip the normals
	//vtkNew<vtkReverseSense> reverseFilter;
	//reverseFilter->SetInputConnection(norms->GetOutputPort());
	//reverseFilter->ReverseCellsOn();
	//reverseFilter->ReverseNormalsOn();
	//reverseFilter->Update();

	vtkNew<vtkSTLWriter> stlWriter;
	stlWriter->SetFileName(fileName.c_str());
	//stlWriter->SetInputConnection(isoFilter->GetOutputPort());
	stlWriter->SetInputConnection(norms->GetOutputPort());
	//stlWriter->SetInputConnection(reverseFilter->GetOutputPort());
	stlWriter->Write();
};