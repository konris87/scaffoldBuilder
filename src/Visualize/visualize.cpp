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
#include <vtkAppendFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkImplicitModeller.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkStlWriter.h>
#include <vtkProperty.h>
#include <vtkPolyDataNormals.h>
#include <vtkReverseSense.h>
#include <vtkSphereSource.h>
#include <vtkVertex.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkLabeledDataMapper.h>
#include <vtkActor2D.h>
#include <vtkTextProperty.h>
#include <vtkNamedColors.h>
#include <vtkGlyph3DMapper.h>
#include <vtkAxesActor.h>


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

void render_vtk_points(const Eigen::MatrixXd& vertices, const std::string& name) {

	vtkNew<vtkNamedColors> colors;

	// Create points
	vtkNew<vtkPoints> points;

	if (vertices.cols() == 3) {
		for (int i = 0; i < vertices.rows(); ++i) {
			points->InsertNextPoint(vertices(i, 0), vertices(i, 1), vertices(i, 2));
		}
	}
	else if (vertices.cols() == 2) {
		for (int i = 0; i < vertices.rows(); ++i) {
			points->InsertNextPoint(vertices(i, 0), vertices(i, 1), 0.0);
		}
	}


	// Create polydata and set points
	vtkNew<vtkPolyData> polyData;
	polyData->SetPoints(points);

	// Create a mapper and actor.
	vtkNew<vtkPolyDataMapper> pointMapper;
	pointMapper->SetInputData(polyData);

	vtkNew<vtkActor> pointActor;
	pointActor->SetMapper(pointMapper);
	pointActor->GetProperty()->SetPointSize(1);
	pointActor->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());

	// create actor
	auto bounds = points->GetBounds();
	double maxLen = 0;
	for (int i = 1; i < 3; ++i)
	{
		maxLen = std::max(bounds[i + 1] - bounds[i], maxLen);
	}

	vtkNew<vtkSphereSource> sphereSource;
	sphereSource->SetRadius(0.01 * maxLen);

	vtkNew<vtkPolyData> pd;
	pd->SetPoints(points);

	vtkNew<vtkGlyph3DMapper> mapper;
	mapper->SetInputData(pd);
	mapper->SetSourceConnection(sphereSource->GetOutputPort());
	mapper->ScalarVisibilityOff();
	mapper->ScalingOff();

	vtkNew<vtkActor> actor;
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(colors->GetColor3d("Gold").GetData());

	vtkNew<vtkLabeledDataMapper> labelMapper;
	labelMapper->SetInputData(polyData);
	labelMapper->GetLabelTextProperty()->SetColor(
		colors->GetColor3d("Magenta").GetData());
	labelMapper->GetLabelTextProperty()->SetFontSize(20.0);

	vtkNew<vtkActor2D> labelActor;
	labelActor->SetMapper(labelMapper);

	vtkNew<vtkAxesActor> axes;
	axes->AxisLabelsOff();
	axes->SetTotalLength(5, 5, 5);

	// Create a renderer, render window, and interactor.
	vtkNew<vtkRenderer> renderer;
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->AddRenderer(renderer);
	
	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Add the actor to the scene.
	// renderer->AddActor(pointActor);
	renderer->AddActor(actor);
	renderer->AddActor(labelActor);
	renderer->AddActor(axes);

	renderer->SetBackground(colors->GetColor3d("DarkSlateGray").GetData());

	// Render and interact
	renderWindow->Render();
	renderWindow->SetWindowName(name.c_str());
	renderWindowInteractor->Start();

};

void render_vtk_face(const Eigen::MatrixXd& vertices, const std::vector<std::vector<int>>& indices, const std::string& name) {

	vtkNew<vtkNamedColors> colors;

	// Create points
	vtkNew<vtkPoints> points;

	if (vertices.cols() == 3) {
		for (int i = 0; i < vertices.rows(); ++i) {
			points->InsertNextPoint(vertices(i, 0), vertices(i, 1), vertices(i, 2));
		}
	}
	else if (vertices.cols() == 2) {
		for (int i = 0; i < vertices.rows(); ++i) {
			points->InsertNextPoint(vertices(i, 0), vertices(i, 1), 0.0);
		}
	}

	vtkNew<vtkCellArray> cells;
	for (int idx{ 0 }; idx < indices.size(); idx++) {
		cells->InsertNextCell(3);
		cells->InsertCellPoint(indices[idx][0]);
		cells->InsertCellPoint(indices[idx][1]);
		cells->InsertCellPoint(indices[idx][2]);
	}

	// Create polydata and set points
	vtkNew<vtkPolyData> polyData;
	polyData->SetPoints(points);
	polyData->SetPolys(cells);

	// Create a mapper and actor.
	vtkNew<vtkPolyDataMapper> mapper;
	mapper->SetInputData(polyData);

	vtkNew<vtkActor> actor;
	actor->SetMapper(mapper);
	actor->GetProperty()->SetLineWidth(2.0f);
	actor->GetProperty()->SetColor(
		colors->GetColor3d("Yellow").GetData()
	);

	// Wireframe actor (overlay)
	vtkNew<vtkActor> wireframeActor;
	wireframeActor->SetMapper(mapper);
	wireframeActor->GetProperty()->SetRepresentationToWireframe();
	wireframeActor->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
	wireframeActor->GetProperty()->SetLineWidth(2.0f);
	wireframeActor->GetProperty()->SetOpacity(1.0);

	// Create a renderer, render window, and interactor.
	vtkNew<vtkRenderer> renderer;
	renderer->AddActor(actor);
	renderer->AddActor(wireframeActor);
	renderer->SetBackground(colors->GetColor3d("MidnightBlue").GetData());
	
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->AddRenderer(renderer);

	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Render and interact
	renderWindow->Render();
	renderWindow->SetWindowName(name.c_str());
	renderWindowInteractor->Start();

};

void render_vtk_triangular_cell(const std::vector<vtkSmartPointer<vtkPolyData>>& polys) {

	vtkNew<vtkNamedColors> colors;

	// update the polydata filter and render

	std::cout << "Polys Size: " << polys.size() << std::endl;
	vtkNew<vtkAppendPolyData> appendFilter;
	for (int i = 0; i < polys.size(); i++) {
		appendFilter->AddInputData(polys[i]);
	}
	appendFilter->Update();

	// Create a mapper and actor.
	vtkNew<vtkPolyDataMapper> mapper;
	mapper->SetInputData(appendFilter->GetOutput());

	vtkNew<vtkActor> actor;
	actor->SetMapper(mapper);
	actor->GetProperty()->SetLineWidth(2.0f);
	actor->GetProperty()->SetColor(
		colors->GetColor3d("Yellow").GetData()
	);

	// Create a renderer, render window, and interactor.
	vtkNew<vtkRenderer> renderer;
	renderer->AddActor(actor);

	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->AddRenderer(renderer);

	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Render and interact
	renderWindow->Render();
	renderWindow->SetWindowName("Triangular Voronoi Cell");
	renderWindowInteractor->Start();

};

vtkSmartPointer<vtkPolyData> create_face_poly(const Eigen::MatrixXd& vertices, const std::vector<std::vector<int>>& indices) {

	// Create points
	vtkNew<vtkPoints> points;

	if (vertices.cols() == 3) {
		for (int i = 0; i < vertices.rows(); ++i) {
			points->InsertNextPoint(vertices(i, 0), vertices(i, 1), vertices(i, 2));
		}
	}
	else if (vertices.cols() == 2) {
		for (int i = 0; i < vertices.rows(); ++i) {
			points->InsertNextPoint(vertices(i, 0), vertices(i, 1), 0.0);
		}
	}

	vtkNew<vtkCellArray> cells;
	for (int idx{ 0 }; idx < indices.size(); idx++) {
		cells->InsertNextCell(3);
		cells->InsertCellPoint(indices[idx][0]);
		cells->InsertCellPoint(indices[idx][1]);
		cells->InsertCellPoint(indices[idx][2]);
	}

	// Create polydata and set points
	vtkNew<vtkPolyData> polyData;
	polyData->SetPoints(points);
	polyData->SetPolys(cells);

	return polyData;

};

void render_vtk_polydata(vtkSmartPointer<vtkPolyData>& polyData) {

	vtkNew<vtkNamedColors> colors;

	// Create a mapper and actor.
	vtkNew<vtkPolyDataMapper> mapper;
	mapper->SetInputData(polyData);

	vtkNew<vtkActor> actor;
	actor->SetMapper(mapper);
	actor->GetProperty()->SetLineWidth(2.0f);
	actor->GetProperty()->SetColor(
		colors->GetColor3d("Yellow").GetData()
	);

	// Create a renderer, render window, and interactor.
	vtkNew<vtkRenderer> renderer;
	renderer->AddActor(actor);

	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->AddRenderer(renderer);

	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Render and interact
	renderWindow->Render();
	renderWindow->SetWindowName("Final Mesh");
	renderWindowInteractor->Start();
};


