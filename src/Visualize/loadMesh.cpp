#include "Visualize/loadMesh.h"

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkCellArrayIterator.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkNew.h>

void loadStlMesh(
	const std::string modelFile,
	std::vector<float>& vertices,
	std::vector<unsigned int>& indices,
	std::vector<float>& normals
) {
	vtkNew<vtkSTLReader> reader;
	reader->SetFileName(modelFile.c_str());
	reader->Update();

	std::cout << "Read Data" << std::endl;

	// get polydata
	vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();

	vtkSmartPointer<vtkPoints> points = polydata->GetPoints();
	vtkSmartPointer<vtkCellArray> cells = polydata->GetPolys();
	
	std::cout << "Estimating norms" << std::endl;
	vtkNew<vtkPolyDataNormals> norms;

	norms->SetInputData(polydata);
	norms->ComputePointNormalsOn();
	//norms->ComputeCellNormalsOn();
	norms->ConsistencyOn();
	//norms->AutoOrientNormalsOn();
	norms->Update();

	// Get normals from the updated polydata
	vtkSmartPointer<vtkDataArray> normalData = norms->GetOutput()->GetPointData()->GetNormals();
	if (!normalData) {
		std::cerr << "Error: No normals found!" << std::endl;
		return;
	}
	else {
		for (vtkIdType i{ 0 }; i < normalData->GetNumberOfComponents(); i++) {
			double normal[3];

			normalData->GetTuple(i, normal);

			normals.emplace_back(normal[0]);
			normals.emplace_back(normal[1]);
			normals.emplace_back(normal[2]);
		}
	}

	for (vtkIdType i{ 0 }; i < points->GetNumberOfPoints(); i++) {
		double point[3];

		points->GetPoint(i, point);

		vertices.emplace_back(point[0]);
		vertices.emplace_back(point[1]);
		vertices.emplace_back(point[2]);
	};

	std::cout << "Getting indices" << std::endl;

	// get indices
	auto cellIter = vtk::TakeSmartPointer(cells->NewIterator());

	//std::cout << cellNr << std::endl;
	for (cellIter->GoToFirstCell(); !cellIter->IsDoneWithTraversal();
		cellIter->GoToNextCell())
	{
		//std::cout << "Cell " << cellIter->GetCurrentCellId() << ":\n";

		vtkSmartPointer<vtkIdList> cell = cellIter->GetCurrentCell();

		indices.emplace_back(cell->GetId(0));
		indices.emplace_back(cell->GetId(1));
		indices.emplace_back(cell->GetId(2));
	}
	//std::cout << "exitiing loading" << std::endl;
};