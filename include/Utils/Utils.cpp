#include "Utils.h"
#include "Visualize/visualize.h"
#include <Eigen/Dense>
#include <vtkTriangle.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkNamedColors.h>
#include <vtkPoints.h>
#include <vtkNew.h>
#include <vtkIdList.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkArrowSource.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include <cstdlib>
#include <vtkPlaneSource.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTriangleFilter.h>
#include <cmath>

bool is_inside_mesh(
	const vtkSmartPointer<vtkPolyData>& mesh,
	const Eigen::Vector3d& point,
	const Eigen::Vector3d& rayDir,
	Eigen::Vector3d& intersection) {
	int count{ 0 };

	Eigen::Vector3d intersec;
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

		if (ray_intersection(point, v1, v2, v3, rayDir, intersec)) {
			count++;
		}
	}

	if (count % 2 == 1) {
		intersection = intersec;
		return true;
	}
	else { return false; }
};

bool is_inside_mesh(
	const vtkSmartPointer<vtkPolyData>& mesh,
	const Eigen::Vector3d& point) {

	Eigen::Vector3d rayDir(1.0, 0.0, 0.0);
	int count{ 0 };

	Eigen::Vector3d intersec;
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

		if (ray_intersection(point, v1, v2, v3, rayDir, intersec)) {
			count++;
		}
	}

	if (count % 2 == 1) {
		return true;
	}
	else { return false; }
};

bool is_inside_box(const std::array<double, 3>& pt, const std::array<float, 6>& bounds) {

	if (pt[0] < bounds[0] || pt[0] > bounds[1] ||
		pt[1] < bounds[2] || pt[1] > bounds[3] ||
		pt[2] < bounds[4] || pt[2] > bounds[5]) {
		return false;
	}
	else {
		return true;
	}

}

bool is_inside_box(const std::array<double, 3>& p, const Bounds& b) {
	return (p[0] >= b.xMin && p[0] <= b.xMax) &&
		(p[1] >= b.yMin && p[1] <= b.yMax) &&
		(p[2] >= b.zMin && p[2] <= b.zMax);
};

Eigen::Vector3d unit_axis_from_dir(int dir) {
	Eigen::Vector3d a(0, 0, 0);
	if (dir == 0) a[0] = 1; else if (dir == 1) a[1] = 1; else a[2] = 1;
	return a;
}

bool is_inside_cylinder(
	const std::array<double, 3>& querryPt,
	const Eigen::Vector3d& basePt,
	const Eigen::Vector3d& axis,
	double radius, double height) {

	Eigen::Vector3d p(querryPt[0], querryPt[1], querryPt[2]);
	Eigen::Vector3d d = p - basePt;
	double projection = d.dot(axis.normalized());
	if (projection < 0 || projection > height)
		return false;

	Eigen::Vector3d radial = d - projection * axis.normalized();
	return radial.norm() <= radius;

};

bool ray_intersection(
	const Eigen::Vector3d& p,
	const Eigen::Vector3d& v1,
	const Eigen::Vector3d& v2,
	const Eigen::Vector3d& v3,
	const Eigen::Vector3d& dir
) {
	
	const float eps = 1e-8;

	// get the edges
	Eigen::Vector3d edge1 = v2 - v1;
	Eigen::Vector3d edge2 = v3 - v1;

	// find the normal
	Eigen::Vector3d rayCrossEdge2 = dir.cross(edge2);
	float det = edge1.dot(rayCrossEdge2);

	// if the dot product is close to zero then the direction is 
	// parallel to the plane that is defined by the triangle 
	// (and orthonormal to the triangle normal)
	if (det > eps && det < -eps) {
		return false;
	}

	// otherwise continue and check if the ray-plane intersection 
	// lies outside the triangle
	float invDet = 1.0 / det;
	Eigen::Vector3d s = p - v1;
	float u = invDet * s.dot(rayCrossEdge2);

	if ((u < 0 && abs(u) > eps) || (u > 1 && abs(u - 1) > eps)) {
		return false;
	}

	Eigen::Vector3d sCrossEdge1 = s.cross(edge1);
	float v = invDet * dir.dot(sCrossEdge1);

	if ((v < 0 && abs(v) > eps) || (u + v > 1 && abs(u + v - 1) > eps)) {
		return false;
	}

	// Now estimate t and check if it is not zero
	float t = invDet * edge2.dot(sCrossEdge1);

	if (t > eps) { return true; }
	else { return false; }
};

bool ray_intersection(
	const Eigen::Vector3d& p,
	const Eigen::Vector3d& v1,
	const Eigen::Vector3d& v2,
	const Eigen::Vector3d& v3,
	const Eigen::Vector3d& dir,
	Eigen::Vector3d& intersection
) {

	const float eps = 1e-8;
	//const float eps = 0.1;

	// get the edges
	Eigen::Vector3d edge1 = v2 - v1;
	Eigen::Vector3d edge2 = v3 - v1;

	// find the normal
	Eigen::Vector3d rayCrossEdge2 = dir.cross(edge2);
	float det = edge1.dot(rayCrossEdge2);

	// if the dot product is close to zero then the direction is 
	// parallel to the plane that is defined by the triangle 
	// (and orthonormal to the triangle normal)
	if (det > eps && det < -eps) {
		return false;
	}

	// otherwise continue and check if the ray-plane intersection 
	// lies outside the triangle
	float invDet = 1.0 / det;
	Eigen::Vector3d s = p - v1;
	float u = invDet * s.dot(rayCrossEdge2);

	if ((u < 0 && abs(u) > eps) || (u > 1 && abs(u - 1) > eps)) {
		return false;
	}

	Eigen::Vector3d sCrossEdge1 = s.cross(edge1);
	float v = invDet * dir.dot(sCrossEdge1);

	if ((v < 0 && abs(v) > eps) || (u + v > 1 && abs(u + v - 1) > eps)) {
		return false;
	}

	// Now estimate t and check if it is not zero
	float t = invDet * edge2.dot(sCrossEdge1);

	if (t > eps) { 
		intersection = p + t * dir;
		//std::cout << "Intersection at: " << intersection.transpose() << std::endl;
		return true; 
	}
	else { return false; }
};

vtkSmartPointer<vtkActor> render_points(
	const int& particles,
	const std::vector<std::array <double, 3>>& currSeeds,
	const std::string color) {

	// colors
	vtkNew<vtkNamedColors> colors;

	// Create a vtkPoints object and store the points in it.
	vtkNew<vtkPoints> seedPoints;

	// Create a polyData object
	vtkNew<vtkPolyData> seedData;

	// Create a cell array
	vtkSmartPointer<vtkCellArray> seedArray = vtkSmartPointer <vtkCellArray>::New();

	// Array of point ids
	//vtkIdType pid[particles];
	//vtkIdType pid[nr];
	std::vector<vtkIdType> pid(currSeeds.size());

	std::cout << "inserting points" << std::endl;
	for (int i{ 0 }; i < currSeeds.size(); i++) {
		
		pid[i] = seedPoints->InsertNextPoint(
			currSeeds[i][0],
			currSeeds[i][1],
			currSeeds[i][2]
		);
		//seedArray->InsertNextCell(i, pid);
		// Create a cell containing a single point
		vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
		idList->InsertNextId(pid[i]);

		// Insert this cell into the cell array
		seedArray->InsertNextCell(idList);
	}

	seedData->SetPoints(seedPoints);
	seedData->SetVerts(seedArray);

	// Create one sphere for all
	vtkNew<vtkSphereSource> sphere;
	sphere->SetPhiResolution(21);
	sphere->SetThetaResolution(21);
	sphere->SetRadius(.08);

	// create vtkglyph to plot points as spheres
	vtkNew<vtkGlyph3D> Glyph3D;
	Glyph3D->SetSourceConnection(sphere->GetOutputPort());
	Glyph3D->SetInputData(seedData);
	Glyph3D->Update();

	// Setup actor and mapper.
	vtkNew<vtkPolyDataMapper> pointMapper;
	pointMapper->SetInputData(seedData);
	pointMapper->SetInputConnection(Glyph3D->GetOutputPort());
	pointMapper->Update();

	vtkSmartPointer<vtkActor> pointActor = vtkSmartPointer<vtkActor>::New();
	pointActor->SetMapper(pointMapper);
	pointActor->GetProperty()->SetColor(colors->GetColor3d(color).GetData());
	pointActor->GetProperty()->SetPointSize(30);
	
	return pointActor;
};

void make_glyphs(vtkPolyData* src, double size, vtkGlyph3D* glyph)
{
	// Source for the glyph filter
	vtkNew<vtkArrowSource> arrow;
	arrow->SetTipResolution(16);
	arrow->SetTipLength(0.3);
	arrow->SetTipRadius(0.1);

	glyph->SetSourceConnection(arrow->GetOutputPort());
	glyph->SetInputData(src);
	glyph->SetVectorModeToUseNormal();
	glyph->SetScaleModeToScaleByVector();
	glyph->SetScaleFactor(size);
	glyph->OrientOn();
	glyph->Update();
}

vtkSmartPointer<vtkActor> render_cell(const std::vector<double>& cellVerts) {

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer< vtkPoints >::New();

	int Nv = cellVerts.size() / 3;

	for (int i = 0; i < Nv; i++) {
		points->InsertNextPoint(
			cellVerts[3 * i],
			cellVerts[3 * i + 1],
			cellVerts[3 * i + 2]);
	}

	vtkSmartPointer< vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);

	vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer< vtkDelaunay3D >::New();
	delaunay->SetInputData(polydata);
	delaunay->Update();

	vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
	surfaceFilter->SetInputConnection(delaunay->GetOutputPort());
	surfaceFilter->Update();

	vtkSmartPointer<vtkPolyDataMapper> surfaceMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	surfaceMapper->SetInputConnection(surfaceFilter->GetOutputPort());
	vtkSmartPointer<vtkActor> surfaceActor = vtkSmartPointer<vtkActor>::New();
	surfaceActor->SetMapper(surfaceMapper);

	double r = static_cast<double>(std::rand()) / RAND_MAX;
	double g = static_cast<double>(std::rand()) / RAND_MAX;
	double b = static_cast<double>(std::rand()) / RAND_MAX;

	surfaceActor->GetProperty()->SetColor(r, g, b);
	surfaceActor->GetProperty()->SetOpacity(0.99);
	surfaceActor->GetProperty()->SetLineWidth(1.0);

	return surfaceActor;
};

vtkSmartPointer<vtkActor> render_plane(
	const double& n1, const double& n2, const double& n3,
	const double& p1, const double& p2, const double& p3,
	const double& s1, const double& s2
) {
	double pi = 3.14159265;

	double normal[3] = { n1, n2, n3 }; // Normal to the plane
	// Ensure the normal is a unit vector
	vtkMath::Normalize(normal);

	vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();
	
	vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
	t->PostMultiply();
	t->Scale(s1, s2, 1);
	double zAxis[3] = { 0.0, 0.0, 1.0 }; // Default plane normal (along Z-axis)
	double rotationAxis[3];
	vtkMath::Cross(zAxis, normal, rotationAxis); // Axis of rotation
	double angle = vtkMath::DegreesFromRadians(acos(vtkMath::Dot(zAxis, normal))); // Rotation angle

	// rodriguez formula
	t->Scale(2, 2, 1);
	if (vtkMath::Norm(rotationAxis) > 1e-6) { // Avoid division by zero
		t->RotateWXYZ(angle, rotationAxis);
	}
	t->Translate(p1, p2, p3);

	vtkSmartPointer<vtkTransformPolyDataFilter> transformPD = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformPD->SetTransform(t);
	transformPD->SetInputConnection(plane->GetOutputPort());
	transformPD->Update();

	// create mapper and actor
	vtkSmartPointer<vtkPolyDataMapper> planeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	planeMapper->SetInputConnection(transformPD->GetOutputPort());
	vtkSmartPointer<vtkActor> planeActor = vtkSmartPointer<vtkActor>::New();
	planeActor->SetMapper(planeMapper);

	return planeActor;
};

double distance(const std::array<double, 3>& arr1, const std::array<double, 3>& arr2) {

	return sqrt(pow(arr2[0] - arr1[0], 2) + pow(arr2[1] - arr1[1], 2) + pow(arr2[2] - arr1[2], 2));
};

void select_file_button(const char* title, const std::string path, const char* dialogName, const char* fileExt) {

	if (ImGui::Button(title)) {
		IGFD::FileDialogConfig config;
		config.path = path;
		ImGuiFileDialog::Instance()->OpenDialog(dialogName, dialogName, fileExt, config);
	}
};

void ensure_ccw(
	std::vector<int>& faceIdxs,
	Eigen::MatrixXd& vertices) {

	// 1. Get the centroid of the face
	Eigen::Vector3d center = vertices.colwise().mean();

	// 2. Compute angles for sorting
	std::vector<double> angles;
	for (int i = 0; i < vertices.rows(); ++i) {
		double dx = vertices(i, 0) - center(0);
		double dy = vertices(i, 1) - center(1);
		angles.push_back(std::atan2(dy, dx)); // 2D angle
	}

	// 3. Sort indices by angle
	std::vector<size_t> sortedIdxs = sort_indices(angles);

	// 4. Reorder both faceIdxs and vertices consistently
	std::vector<int> sortedFace(faceIdxs.size());
	Eigen::MatrixXd orderedVertices(vertices.rows(), vertices.cols());

	for (int i = 0; i < sortedIdxs.size(); ++i) {
		sortedFace[i] = faceIdxs[sortedIdxs[i]];
		orderedVertices.row(i) = vertices.row(sortedIdxs[i]);
	}

	faceIdxs = sortedFace;
	vertices = orderedVertices;
};

Eigen::MatrixXd vector_to_matrix3D(const std::vector<double>& verts) {
	assert(verts.size() % 3 == 0 && "Vertex list must have multiple of 3 elements.");
	int n = static_cast<int>(verts.size() / 3);
	Eigen::MatrixXd mat(n, 3);
	for (int i = 0; i < n; ++i) {
		mat(i, 0) = verts[3 * i + 0];
		mat(i, 1) = verts[3 * i + 1];
		mat(i, 2) = verts[3 * i + 2];
	}
	return mat;
};

void project_vertices_on_plane(
	const std::vector<double>& vertices,
	Eigen::Vector3d& origin,
	Eigen::Vector3d& u,
	Eigen::Vector3d& v,
	Eigen::MatrixXd& vertices2D
){
	int vertNr = static_cast<int>(vertices.size() / 3);

	vertices2D.resize(vertNr, 2);

	// create an Eigen matrix for the vertices
	Eigen::MatrixXd verts(vertNr, 3);

	// populate it
	for (int i{ 0 }; i < vertNr; i++) {
		verts.row(i) = Eigen::Vector3d{
			vertices[3 * i],
			vertices[3 * i + 1],
			vertices[3 * i + 2]
		};
	}	

	//render_vtk_points(verts, "Initial");

	// create a local reference frame using two edges and the normal
	Eigen::Vector3d edge1 = verts.row(1) - verts.row(0);

	Eigen::Vector3d edge2 = verts.row(2) - verts.row(0);

	// get the cross product
	Eigen::Vector3d normal = edge1.cross(edge2).normalized();

	// u axis is 
	Eigen::Vector3d ue = edge1.normalized();

	// v axis
	Eigen::Vector3d ve = normal.cross(ue).normalized();

	// get the projected vertices
	u = ue;
	v = ve;

	// use the centroid as center
	Eigen::Vector3d center = verts.colwise().mean();
	origin = center;

	// use the center as origin
	Eigen::MatrixXd centeredVertices = verts.rowwise() - center.transpose();

	// create a matrix with columns the unit vector
	Eigen::MatrixXd unitVecs(2, 3);

	unitVecs.row(0) = ue.transpose();
	unitVecs.row(1) = ve.transpose();

	vertices2D = (unitVecs * centeredVertices.transpose()).transpose();
};

double erosion_margin_for_face(
	int f,
	const std::vector<Eigen::Vector3d>& normals,
	const std::vector<double>& planeB,
	const Eigen::Matrix<double, 3, 2>& U,
	const Eigen::Vector3d& x0,
	const Eigen::MatrixXd& verts2D, // N x 2
	double delta,                   // pullback
	double eps
) {
	double dstar = std::numeric_limits<double>::infinity();

	for (int j = 0; j < (int)normals.size(); ++j) {
		if (j == f) continue;

		Eigen::Vector2d m = U.transpose() * normals[j];      // in-plane direction
		double m2 = m.squaredNorm();
		if (m2 < 1e-16) continue;                            // plane parallel to face

		double cj = planeB[j] - normals[j].dot(x0);

		// support over polygon in direction m
		double maxDot = -std::numeric_limits<double>::infinity();
		for (int k = 0; k < verts2D.rows(); ++k) {
			double val = m.dot(verts2D.row(k));
			if (val > maxDot) maxDot = val;
		}

		// d*_j = (c_j - supp_P(m)) + delta * ||m||
		double dstar_j = (cj - maxDot) + delta * std::sqrt(m2);
		if (dstar_j < dstar) dstar = dstar_j;
	}
	return (dstar > eps) ? dstar : 0.0;
}


void project_vertices_on_plane(
	Eigen::MatrixXd& verts,
	Eigen::Vector3d& origin,
	Eigen::Vector3d& u,
	Eigen::Vector3d& v,
	Eigen::MatrixXd& vertices2D
) {
	//int vertNr = static_cast<int>(verts.size() / 3);

	//// create an Eigen matrix for the vertices
	//Eigen::MatrixXd verts(vertNr, 3);

	//// populate it
	//for (int i{ 0 }; i < vertNr; i++) {
	//	verts.row(i) = Eigen::Vector3d{
	//		vertices[3 * i],
	//		vertices[3 * i + 1],
	//		vertices[3 * i + 2]
	//	};
	//}

	//render_vtk_points(verts, "Initial");

	// create a local reference frame using two edges and the normal
	Eigen::Vector3d edge1 = verts.row(1) - verts.row(0);

	Eigen::Vector3d edge2 = verts.row(2) - verts.row(0);

	// get the cross product
	Eigen::Vector3d normal = edge1.cross(edge2).normalized();

	// u axis is 
	Eigen::Vector3d ue = edge1.normalized();

	// v axis
	Eigen::Vector3d ve = normal.cross(ue).normalized();

	// get the projected vertices
	u = ue;
	v = ve;

	// use the centroid as center
	Eigen::Vector3d center = verts.colwise().mean();
	origin = center;

	// use the center as origin
	Eigen::MatrixXd centeredVertices = verts.rowwise() - center.transpose();

	// create a matrix with columns the unit vector
	Eigen::MatrixXd unitVecs(2, 3);

	unitVecs.row(0) = ue.transpose();
	unitVecs.row(1) = ve.transpose();

	vertices2D = (unitVecs * centeredVertices.transpose()).transpose();
};

void back_to_3d(Eigen::MatrixXd& vertices3d, const Eigen::MatrixXd& vertices2d, const Eigen::Vector3d& center, const Eigen::Vector3d& u, const Eigen::Vector3d& v) {

	vertices3d.resize(vertices2d.rows(), 3);

	for (int i = 0; i < vertices2d.rows(); ++i) {
		Eigen::Vector3d row = center + u * vertices2d(i, 0) + v * vertices2d(i, 1);
		vertices3d.row(i) = row.transpose();
	}

};

double catmull_rom_get_t(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, double t0, double alpha) {

	double x = std::pow(p2.x() - p1.x(), 2);
	double y = std::pow(p2.y() - p1.y(), 2);

	double ti = std::pow(std::sqrt(x + y), alpha) + t0;

	return ti;
};

void catmull_rom_interpolation(const Eigen::MatrixXd& currentVerts, Eigen::MatrixXd& interpolatedVerts, double alpha) {

	std::vector<Eigen::Vector2d> samples;
	const int sampleNr = 3;

	for (int vIdx{ 0 }; vIdx < currentVerts.rows(); vIdx++) {

		int idx1 = (vIdx - 1 + currentVerts.rows()) % currentVerts.rows();
		int idx2 = vIdx;
		int idx3 = (vIdx + 1) % currentVerts.rows();
		int idx4 = (vIdx + 2) % currentVerts.rows();

		std::cout << idx1 << " " << idx2 << " " << idx3 << " " << idx4 << std::endl;

		Eigen::Vector2d p0 = currentVerts.row(idx1);
		Eigen::Vector2d p1 = currentVerts.row(idx2);
		Eigen::Vector2d p2 = currentVerts.row(idx3);
		Eigen::Vector2d p3 = currentVerts.row(idx4);

		double t0 = 0.0;
		double t1 = catmull_rom_get_t(p0, p1, t0, alpha);
		double t2 = catmull_rom_get_t(p1, p2, t1, alpha);
		double t3 = catmull_rom_get_t(p2, p3, t2, alpha);


		for (int i{ 0 }; i < sampleNr - 1; i ++) {

			double ti = t1 + (t2 - t1) * (double(i) / (sampleNr - 1));

			Eigen::VectorXd a1 = ((t1 - ti) / (t1 - t0)) * p0 + ((ti - t0) * (t1 - t0)) * p1;
			Eigen::VectorXd a2 = ((t2 - ti) / (t2 - t1)) * p1 + ((ti - t1) * (t2 - t1)) * p2;
			Eigen::VectorXd a3 = ((t3 - ti) / (t3 - t2)) * p2 + ((ti - t2) * (t3 - t2)) * p3;

			Eigen::VectorXd b1 = ((t2 - ti) / (t2 - t0)) * a1 + ((ti - t0) / (t2 - t0)) * a2;
			Eigen::VectorXd b2 = ((t3 - ti) / (t3 - t1)) * a2 + ((ti - t1) / (t3 - t1)) * a3;

			// this is the interpolation point
			Eigen::VectorXd c = ((t2 - ti) / (t2 - t1)) * b1 + ((ti - t1) / (t2 - t1)) * b2;

			//if (i < sampleNr - 1) {
			//	samples.push_back(c);
			//}
			samples.push_back(c);
		}
	}

	// Convert to Eigen::MatrixXd
	interpolatedVerts.resize(samples.size(), 2);
	for (int i = 0; i < samples.size(); ++i) {
		interpolatedVerts.row(i) = samples[i];
	}
};

void chaikin_subdivision(const Eigen::MatrixXd& currentVerts, Eigen::MatrixXd& interpolatedVerts, int num) {

	std::vector<Eigen::Vector2d> samples;

	for (int i{ 0 }; i < currentVerts.rows(); i++) {

		int idx1 = i;
		int idx2 = (i + 1) % currentVerts.rows();

		Eigen::Vector2d p0 = currentVerts.row(idx1);
		Eigen::Vector2d p1 = currentVerts.row(idx2);

		//double norm = (p0 - p1).norm();
		//if (norm < 0.1) {
		//	continue;
		//}

		samples.push_back(0.75 * p0 + 0.25 * p1);
		samples.push_back(0.25 * p0 + 0.75 * p1);
	}
	interpolatedVerts.resize(samples.size(), 2);
	for (int i{ 0 }; i < samples.size(); i++) {
		interpolatedVerts.row(i) = samples[i];
	}
}

PolygonWidth get_polygon_width(const Eigen::MatrixXd& verts2D) {

	PolygonWidth stats = { 0.0, 0.0 };

	const int vertNr = verts2D.rows();

	if (vertNr < 3) {
		return stats;
	}

	double minWidth = std::numeric_limits<double>::infinity();
	double maxWidth = 0.0;

	// loop for each edge
	for (int i{ 0 }; i < vertNr; i++) {
		int j = (i + 1) % vertNr;

		Eigen::Vector2d p1 = verts2D.row(i);
		Eigen::Vector2d p2 = verts2D.row(j);

		Eigen::Vector2d edge = p2 - p1;

		// check if the edge is very small
		double edgeNorm = edge.norm();

		if (edgeNorm < 1e-20) continue;

		// estimate the normal of the edge
		Eigen::Vector2d normal(-edge.y() / edgeNorm, edge.x() / edgeNorm);

		// now for all vertices estimate the scalar projections
		double smin = std::numeric_limits<double>::infinity();
		double smax = -std::numeric_limits<double>::infinity();

		for (int v{ 0 }; v < vertNr; v++) {

			Eigen::Vector2d p = verts2D.row(v);
			double sp = normal.dot(p);
			smin = std::min(smin, sp);
			smax = std::max(smax, sp);
		}

		// estimate the width
		double width = smax - smin;
		minWidth = std::min(minWidth, width);
		maxWidth = std::max(maxWidth, width);
	}

	// check if minWidth is not changed
	if (minWidth == std::numeric_limits<double>::infinity()) minWidth = 0.0;

	stats.minWidth = minWidth;
	stats.maxWidth = maxWidth;

	return stats;
};

void sample_face_polygon(
	const Eigen::MatrixXd& verts2D,
	Eigen::Vector2d& center,
	double& radius1, 
	int resolution,
	double scale) {

	// 1. find the min and max for the x and y coordinates of the vertices
	double xMin = verts2D.col(0).minCoeff();
	double xMax = verts2D.col(0).maxCoeff();

	double yMin = verts2D.col(1).minCoeff();
	double yMax = verts2D.col(1).maxCoeff();

	// 2. create x and y points
	Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(resolution, xMin, xMax);
	Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(resolution, yMin, yMax);

	// 3. now loop for each point and find the biggest radius estimating distance from edges
	double maxDist{ 0.0 };
	int ptIdx = 0;
	for (auto xi : x) {
		for (auto yi : y) {

			Eigen::Vector2d pt{ xi, yi };

			// 4. first check if the pt is inside the polygon defined by the vertices
			if (is_inside_polygon(pt, verts2D)) {

				// 5. find distance from edges, the minimum corresponds to the radius
				// we also use a scale factor to have some space between the hole and the face edges
				double distance = scale * min_distance_from_edges(pt, verts2D, ptIdx);

				// 6. if the distance is larger than maxDist then this point is the center
				// of the largest circle
				if (distance > maxDist) {
					maxDist = distance;
					center = pt;
				}
			}
		}
		ptIdx++;
	}
	radius1 = maxDist;
};


bool is_inside_polygon(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices) {

	bool isInside = false;
	const double eps = 1e-10;

	// 1. loop for each edge
	for (int i{ 0 }; i < vertices.rows(); i++) {

		// 2. get the two points of the edge
		Eigen::VectorXd p1 = vertices.row(i);
		Eigen::VectorXd p2 = vertices.row((i + 1) % vertices.rows());

		// 3. first check if the point is one of the edge vertices
		if ((pt - p1).norm() < eps || (pt - p2).norm() < eps) {
			return false;
		}

		// 4. check if a point is on the edge
		Eigen::Vector2d edge = p2 - p1;
		Eigen::Vector2d toPt = pt - p1;
		double cross = edge.x() * toPt.y() - edge.y() * toPt.x();

		if (std::abs(cross) < eps &&
			pt.x() >= std::min(p1.x(), p2.x()) &&
			pt.x() <= std::max(p1.x(), p2.x()) &&
			pt.y() >= std::min(p1.y(), p2.y()) &&
			pt.y() <= std::max(p1.y(), p2.y())) {
			return false;  // Consider point on edge as not inside
		}

		// 4. then check if the point is  between ymin and ymax of the edge
		if ((p1.y() > pt.y()) != (p2.y() > pt.y())) {

			// 5. estimate the x intersect. if the point has a smaller x then we have an intersection
			if ((p1.y() > pt.y()) != (p2.y() > pt.y())) {
				double xIntersect = p1.x() + (pt.y() - p1.y()) * (p2.x() - p1.x()) / (p2.y() - p1.y());
				if (pt.x() < xIntersect)
					isInside = !isInside;
			}

		}
	}
	return isInside;
};

bool is_inside_triangle(const Eigen::Vector2d& vPrev, const Eigen::Vector2d& vCurr, const Eigen::Vector2d& vNext, const Eigen::Vector2d& vTest) {

	Eigen::Vector2d v0 = vNext - vPrev;
	Eigen::Vector2d v1 = vCurr - vPrev;
	Eigen::Vector2d v2 = vTest - vPrev;

	// dot products
	double dot00 = v0.dot(v0);
	double dot01 = v0.dot(v1);
	double dot02 = v0.dot(v2);
	double dot11 = v1.dot(v1);
	double dot12 = v1.dot(v2);

	// barycentric coordinates
	double denom = (dot00 * dot11 - dot01 * dot01);

	if (std::abs(denom) - 1e-20) {
		return true;
	}
	
	double invDenom = 1.0 / denom;

	//u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	//v = (dot00 * dot12 - dot01 * dot02) * invDenom;

};

double min_distance_from_edges(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices, int excludeIdx) {

	double minDist = std::numeric_limits<double>::infinity();

	// 1. loop for each edge
	for (int i{ 0 }; i < vertices.rows(); i++) {

		int j = (i + 1) % vertices.rows();

		// dont check the edges that contain the vertex idx
		if (i == excludeIdx || j == excludeIdx) continue;

		// 2. get the two points of the edge
		Eigen::VectorXd p1 = vertices.row(i);
		Eigen::VectorXd p2 = vertices.row(j);

		// 3. estimate a vector for the edge
		Eigen::Vector2d edge = p2 - p1;

		// 4. estimate the vector from p1 to pt
		Eigen::Vector2d toPt = pt - p1;

		// 5. we get the projection of the pt to the edge
		double t = edge.dot(toPt) / edge.squaredNorm();

		// 6. clamp it to be on the segment
		t = std::clamp(t, 0.0, 1.0); 

		// 7. the vector projection is 
		Eigen::Vector2d projection = p1 + t * edge;
		
		// 8. the distance is the norm of the vertical vector
		double dist = (pt - projection).norm();

		// 9. check if the dist is smaller than minDists
		minDist = std::min(minDist, dist);
	}

	return minDist;

};

double max_distance_from_edges(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices, int excludeIdx) {

	double maxDist = std::numeric_limits<double>::infinity();

	// 1. loop for each edge
	for (int i{ 0 }; i < vertices.rows(); i++) {

		int j = (i + 1) % vertices.rows();

		// dont check the edges that contain the vertex idx
		if (i == excludeIdx || j == excludeIdx) continue;

		// 2. get the two points of the edge
		Eigen::VectorXd p1 = vertices.row(i);
		Eigen::VectorXd p2 = vertices.row(j);

		// 3. estimate a vector for the edge
		Eigen::Vector2d edge = p2 - p1;

		// 4. estimate the vector from p1 to pt
		Eigen::Vector2d toPt = pt - p1;

		// 5. we get the projection of the pt to the edge
		double t = edge.dot(toPt) / edge.squaredNorm();

		// 6. clamp it to be on the segment
		t = std::clamp(t, 0.0, 1.0);

		// 7. the vector projection is 
		Eigen::Vector2d projection = p1 + t * edge;

		// 8. the distance is the norm of the vertical vector
		double dist = (pt - projection).norm();

		// 9. check if the dist is smaller than minDists
		maxDist = std::max(maxDist, dist);
	}

	return maxDist;

};

double polygon_min_width(const Eigen::MatrixXd& vertices) {

	double minWidth = std::numeric_limits<double>::infinity();

	for (int i{ 0 }; i < vertices.rows(); i++) {

		double width = min_distance_from_edges(vertices.row(i), vertices, i);

		minWidth = std::min(minWidth, width);
	}

	return minWidth;
};

double polygon_max_width(const Eigen::MatrixXd& vertices) {

	double maxWidth = -1.0;

	for (int i{ 0 }; i < vertices.rows(); i++) {

		double width = max_distance_from_edges(vertices.row(i), vertices, i);

		maxWidth = std::max(maxWidth, width);
	}

	return maxWidth;
};

double polygon_average_edge_length(const Eigen::MatrixXd& vertices) {

	double avgLength{ 0.0 };

	if (vertices.rows() < 2) return 0.0;

	for (int i{ 0 }; i < vertices.rows(); i++) {

		int j = (i + 1) % vertices.rows();

		// 2. get the two points of the edge
		Eigen::VectorXd p1 = vertices.row(i);
		Eigen::VectorXd p2 = vertices.row(j);

		// 3. estimate edge length
		avgLength += (p2 - p1).norm();
	}

	return avgLength / vertices.rows();
};


void interpolate_edges(
	const Eigen::MatrixXd& vertices,
	Eigen::MatrixXd& interpolatedVertices,
	std::vector<int>& newLocalFace,
	const double& edgeSize) {

	// vector to store temp coords
	std::vector<Eigen::Vector2d> tempPoints;

	// add all the vertices to the tempPoints
	
	for (int i = 0; i < vertices.rows(); ++i) {
		tempPoints.push_back(vertices.row(i));
	}
	
	int vIdx = vertices.rows();

	//newLocalFace.push_back(vIdx);
	// 1. loop for each edge
	for (int i{ 0 }; i < vertices.rows(); i++) {
		
		// 2. get the two points of the edge
		int idx1 = i;
		int idx2 = (i + 1) % vertices.rows();

		Eigen::Vector2d p1 = vertices.row(idx1);
		Eigen::Vector2d p2 = vertices.row(idx2);

		// 3. find the interpolation point number
		int nr = static_cast<int>(std::ceil((p2 - p1).norm() / edgeSize));
		nr = std::max(nr, 2);

		// 4. get the alpha values
		Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(nr, 0.0, 1.0);

		//std::cout << idx1 << " " << std::endl;
		newLocalFace.push_back(idx1);
		for (int newIdx{ 0 }; newIdx < t.size() - 2; newIdx++) {
			//std::cout << vIdx << " " << std::endl;
			newLocalFace.push_back(vIdx);
			vIdx++;
		}

		for (int j{ 1 }; j < t.size() - 1; j++) {
			double alpha = t[j];
			// get the interpolated value
			Eigen::Vector2d interp = p1 * (1 - alpha) + p2 * alpha;
			tempPoints.push_back(interp);
		}
	}

	// 6. update the interpolated vertices
	interpolatedVertices.resize(tempPoints.size(), 2);

	for (int  i= 0; i < tempPoints.size(); ++i) {
		interpolatedVertices.row(i) = tempPoints[i];
	}
};

void hole_points(const double& radius, const Eigen::Vector2d& center, const int& ptNr, Eigen::MatrixXd& vertices) {

	vertices.resize(ptNr, 2);

	// 1. create the hole vertices
	for (int i{ 0 }; i < ptNr; i++) {

		double theta = i * 2 *  M_PI / ptNr;

		Eigen::Vector2d pt = { radius * std::cos(theta), -radius * std::sin(theta) };

		vertices(i, 0) = pt(0) + center(0);
		vertices(i, 1) = pt(1) + center(1);
	}
};

bool vertex_locally_convex(const Eigen::Vector2d& v1, const Eigen::Vector2d& v2, const Eigen::Vector2d& v3) {

	Eigen::Vector2d edge1 = v2 - v1;
	Eigen::Vector2d edge2 = v3 - v1;

	// estimate the angle
	double dot = edge1.dot(edge2);
	double det = edge1(0) * edge2(1) - edge1(1) * edge2(0);

	double angle = std::atan2(det, dot);

	if (angle < 0.0) {
		angle = 2.0 * M_PI + angle;
	}

	return angle < M_PI;
};

bool ear_clipping(const Eigen::MatrixXd& vertices, const std::vector<int>& idxs, std::vector<std::vector<int>>& cells) {

	// 1. create dll for the indices
	Cdll cdll;

	for (auto idx : idxs) {
		cdll.append(idx);
	}

	//cdll.display();

	int iterNr = 1;

	Node* node = cdll.head;

	while (cdll.length > 2) {

		iterNr += 1;
		if (iterNr > 200) {
			//render_vtk_points(vertices, "Face with issue");
			return false;
		}

		//std::cout << "..." << std::endl;
		// 2. get the triad indices
		int i{ node->prev->data };
		int j{ node->data };
		int k{ node->next->data };

		// 3. get the correspondin vertices
		Eigen::Vector2d vPrev = vertices.row(i);
		Eigen::Vector2d vCurr = vertices.row(j);
		Eigen::Vector2d vNext = vertices.row(k);

		// 4. check if the current vertex is convex
		bool isConvex = vertex_locally_convex(vPrev, vCurr, vNext);

		// 5. variable to check if cuurent vertex is an ear
		bool isEar{ true };
		
		// 6. if it is convex continue
		if (isConvex) {

			// now we have to check if the other nodes are inside the triangle
			Node* testNode = node->next->next;

			while (testNode != node->prev && isEar) {

				if ((testNode->data != i) && (testNode->data != j) && (testNode->data != k)) {

					// get the vertex
					Eigen::Vector2d vTest = vertices.row(testNode->data);

					Eigen::MatrixXd otherVerts (3, 2);

					otherVerts.row(0) = vPrev;
					otherVerts.row(1) = vCurr;
					otherVerts.row(2) = vNext;

					bool isInside = is_inside_polygon(vTest, otherVerts);

					if (isInside) {
						isEar = false;
					}
					//else {
					//	std::cout << "outside triangle" << std::endl;
					//}
				}
				testNode = testNode->next;

			}
		}
		else {
			isEar = false;
		}

		if (isEar) {
			//std::cout << "Adding triangle: " << i << " - " << j << " - " << k << std::endl;
			cells.push_back({ i, j, k });
			cdll.remove(node);
		}
		node = node->next;
	}

	return true;
	//std::cout << "..." << std::endl;
};

bool polygon_is_open(
	const std::vector<Eigen::Vector3d>& normals,
	const std::vector<double>& d,
	const Eigen::MatrixXd& faceVertices,
	int fIdx,
	double erosion,
	double eps
) {
	// loop for all face vertices
	//double dstar = std::numeric_limits<double>::infinity();

	for (int i{ 0 }; i < faceVertices.rows(); i++) {

		// get the vertex coords
		Eigen::Vector3d p = faceVertices.row(i);

		// now loop for all other faces
		for (int j{ 0 }; j < static_cast<int>(normals.size()); j++) {
			if (j == fIdx) continue;

			// now check if the vertex will be above or below the eroded plane
			// we will accept if normal_j * p > d_j - erosion + eps
			double diff = d[j] - normals[j].dot(p);
			if (normals[j].dot(p) < d[j] - erosion + eps) return false;
			//dstar = std::min(dstar, diff);
		}
	}
	//return (dstar > eps) ? dstar : 0.0;
	return true;
	//return (dstar > eps) ? dstar : 0.0;
};


// @function for half space clipping using the Sutherland - Hodgman Algorithm
Eigen::MatrixXd clip_polygon(
	const Eigen::MatrixXd& verts2D,
	const Eigen::Vector2d normal,
	const double d
) {

	// if it is empty return it
	if (verts2D.rows() == 0) {
		return verts2D;
	}

	// result
	std::vector<Eigen::Vector2d> temp;

	for (int i{ 0 }; i < verts2D.rows(); i++) {

		// test the edge between two consecutive points
		Eigen::Vector2d p = verts2D.row(i);
		Eigen::Vector2d q = verts2D.row((i+1)%verts2D.rows());
		
		// check if these are inside the halfspace
		bool pInside = (p.dot(normal) <= d + 1e-12);
		bool qInside = (q.dot(normal) <= d + 1e-12);

		// if both inside keep q
		if (pInside && qInside) { 
			temp.push_back(q); 
		}
		// if first inside and second out push the intersection
		else if (pInside && !qInside) {
			temp.push_back(intersection_of_lines(p, q, normal, d));
		}
		else if (!pInside && qInside) {
			temp.push_back(intersection_of_lines(p, q, normal, d));
			temp.push_back(q);
		}
	}
	
	// return a matrix
	Eigen::MatrixXd result(temp.size(), 2);
	for (int j{ 0 }; j < temp.size(); j++) {
		result.row(j) = temp[j];
	}
	return result;
};

// @function to find intersection of lines
Eigen::Vector2d intersection_of_lines(
	const Eigen::Vector2d& p,
	const Eigen::Vector2d& q,
	const Eigen::Vector2d& normal,
	double d
) {

	double denom = normal.dot(q - p);
	if (std::abs(denom) > 1e-15) return q;

	double t = (d - p.dot(normal)) / denom;

	return p + t * (q - p);
};

bool clipping_is_valid(
	const Eigen::MatrixXd& verts2D,
	const std::vector<HalfSpace> hspaces,
	double delta,
	Eigen::MatrixXd* outPoly
) {
	Eigen::MatrixXd poly = verts2D;
	for (const auto& hs : hspaces) {
		poly = clip_polygon(poly, hs.n, hs.d - delta);

		// if the resulted polygon has fewer that 3 vertices then
		// it is degenerated
		if (poly.rows() < 3) {
			if (outPoly) outPoly->resize(0, 2);
			return false;
		}
	}
	if (outPoly) {
		*outPoly = poly;
	}
	return true;
};

std::vector<HalfSpace> get_poly_half_spaces(const Eigen::MatrixXd& verts2D) {
	
	const int vertNr = verts2D.rows();
	std::vector<HalfSpace> hSpaces;
	hSpaces.reserve(vertNr);

	// get the centroid
	Eigen::Vector2d centroid = verts2D.colwise().mean();
	
	for (int i{ 0 }; i < vertNr; i++) {
		//get each edge
		Eigen::Vector2d p = verts2D.row(i);
		Eigen::Vector2d q = verts2D.row((i + 1) % vertNr);
		Eigen::Vector2d e = q - p;


		// get the inward normal it is the left rotated normalized edge
		Eigen::Vector2d normal{ -e.y(), e.x() };

		if (normal.dot(centroid) > normal.dot(p)) normal = -normal;

		double len = normal.norm();
		if (len < 1e-16) continue;
		normal /= len;

		// find the plane coeff
		double d = normal.dot(p);

		hSpaces.push_back({ normal, d });
	}

	return hSpaces;
};

double polygon_inradius(const Eigen::MatrixXd& verts2D) {

	// check if a face is degenerated
	if (verts2D.rows() < 3 || verts2D.cols() != 2) return 0.0;

	// first gather the halfspaces
	std::vector<HalfSpace> hSpaces = get_poly_half_spaces(verts2D);
	
	if (hSpaces.empty()) return 0.0;

	// bisection limits
	double low = 0.0;
	double upper = 1e-6;

	// find an infeasible upper bound
	Eigen::MatrixXd dummy;
	while (clipping_is_valid(verts2D, hSpaces, upper, &dummy)) {
		low = upper;
		upper *= 2.0;
		if (upper > 1e6) break; // safety
	}

	// clip to zero
	upper = std::max(0.0, upper);

	// keep also the best poly to provide the centroid as a center
	Eigen::MatrixXd best;

	// use a number of iterations to find the inradius
	for (int i{ 0 }; i < 25; i++) {
		// take midpoint
		double mid = 0.5 * (low + upper);

		Eigen::MatrixXd poly;
		if (clipping_is_valid(verts2D, hSpaces, mid, &poly)) {
			low = mid;
			best = poly;
		}
		else {
			upper = mid;
		}
	}

	return low;
};


// ----------------------------------------------------------------------------
bool Graph::add_vertex(const int idx) {

	// check if is inside the adjList
	if (adjList.count(idx) == 0) {
		adjList[idx];
		return true;
	}
	return false;
};

void Graph::print() {

	std::cout << " ----------- " << std::endl;
	for (auto [vertex, edges] : adjList) {

		std::string edgeNames{ "" };

		for (auto edge : edges) {
			edgeNames += " '" + std::to_string(edge.first) + " weight " + std::to_string(edge.second.dstar) + "' ";
		}

		std::cout << vertex << " : [" << edgeNames << " ]" << std::endl;
	};
};

//bool Graph::add_edge(const int vertex1, const int vertex2, const float weight) {
bool Graph::add_edge(const int vertex1, const int vertex2, const EdgeData& edge) {

	if (adjList.count(vertex1) != 0 && adjList.count(vertex2) != 0) {
	//if (adjList.count(vertex1) != 0) {
		//adjList.at(vertex1).insert({ vertex2, weight });
		//adjList.at(vertex2).insert({ vertex1, weight });
		adjList.at(vertex1).insert({ vertex2, edge });
		adjList.at(vertex2).insert({ vertex1, edge });
		//std::cout << "added edge: " << vertex1 << " - " << vertex2 << std::endl;
		return true;
	}
	return false;
};

bool Graph::remove_edge(const int vertex1, const int vertex2) {

	if (adjList.count(vertex1) != 0 && adjList.count(vertex2) != 0) {
		adjList.at(vertex1).erase(vertex2);
		adjList.at(vertex2).erase(vertex1);
		return true;
	}
	return false;
};


bool Graph::remove_vertex(const int vertex) {

	if (adjList.count(vertex) == 0) {
		return false;
	}
	// loop and delete edges
	for (auto otherVertex : adjList.at(vertex)) {
		adjList.at(otherVertex.first).erase(vertex);
	}
	adjList.erase(vertex);
	return true;

};

void Graph::remove_edges_below(const float weight) {

	std::vector<std::pair<int, int>> edges_to_remove;

	// 1. Find all edges below the weight
	for (auto const& [vertex1, edges] : adjList) {
		for (auto const& [vertex2, w] : edges) {
			if (w.dstar < weight && vertex1 < vertex2) {
				edges_to_remove.push_back({ vertex1, vertex2 });
			}
		}
	}

	// 2. Now, remove both directions for each edge
	for (auto const& edge_pair : edges_to_remove) {
		remove_edge(edge_pair.first, edge_pair.second); 
		remove_edge(edge_pair.second, edge_pair.first); 
	}
};

void Graph::remove_edges_above(const float weight) {
	std::vector<std::pair<int, int>> edges_to_remove;

	// 1. Find all edges below the weight
	for (auto const& [vertex1, edges] : adjList) {
		for (auto const& [vertex2, w] : edges) {
			if (w.dstar > weight) {
				if (vertex1 < vertex2) {
					edges_to_remove.push_back({ vertex1, vertex2 });
				}
			}
		}
	}

	// 2. Now, remove both directions for each edge
	for (auto const& edge_pair : edges_to_remove) {
		remove_edge(edge_pair.first, edge_pair.second);
		remove_edge(edge_pair.second, edge_pair.first);
	}
};

void Graph::dfs_traversal(int startNode, std::set<int>& visited) {

	// add the start node
	visited.insert(startNode);

	if (adjList.count(startNode)) {
		for (auto const& [idx, weight] : adjList.at(startNode)) {
			if (visited.count(idx) == 0) {
				dfs_traversal(idx, visited);
			}
		}
	}
};

int Graph::find_longest_network() {

	std::set<int> visited;

	int maxLength = 0;

	for (auto const& [vIdx, weight] : adjList) {

		// check inside the set if we already checked it or
		// it is a new component
		if (visited.count(vIdx) == 0) {
			std::set<int> components;

			dfs_traversal(vIdx, components);

			int currentSize = components.size();

			if (currentSize > maxLength) {
				maxLength = currentSize;
			}

			// add the members of the components to the visited
			// to avoid checking them again
			visited.insert(components.begin(), components.end());
		}

	}

	return maxLength;
};

int Graph::get_vertex_count() {
	return adjList.size();
};

//@brief function to get edge data. It is used inside process_faces
EdgeData Graph::get_edge_width(const int vertex1, const int vertex2) {

	// get the set of edges of the vertex1
	auto it1 = adjList.find(vertex1);

	// if not found return -1
	if (it1 == adjList.end()) {
		return { -1.0f, -1.0f };
	}

	auto it_v2 = it1->second.find(vertex2);
	if (it_v2 == it1->second.end()) {
		return {-1.0f, -1.0f}; // edge v1 -> v2 not in graph
	}

	// if ok return the weight
	return it_v2->second;
};

//std::unordered_map<int, std::unordered_map<int, float>> Graph::get_adj_list() {
std::unordered_map<int, std::unordered_map<int, EdgeData>> Graph::get_adj_list() {
	return adjList;
};

