#include "Utils.h"
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

bool is_inside_box(const Eigen::VectorXd& pt, const std::array<double, 6>& bounds) {

	if (pt[0] < bounds[0] || pt[0] > bounds[1] ||
		pt[1] < bounds[2] || pt[1] > bounds[3] ||
		pt[2] < bounds[4] || pt[2] > bounds[5]) {
		return false;
	}
	else {
		return true;
	}

}

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

void selectFileButton(const char* title, const std::string path, const char* dialogName, const char* fileExt) {

	if (ImGui::Button(title)) {
		IGFD::FileDialogConfig config;
		config.path = path;
		ImGuiFileDialog::Instance()->OpenDialog(dialogName, dialogName, fileExt, config);
	}
};

void ensure_ccw(
	std::vector<int>& faceIdxs,
	Eigen::MatrixXd& vertices) {

	// 1. get the centroid of the face
	Eigen::Vector3d center = vertices.colwise().mean();

	// 3. estimate the arctan2 between the columns of the projected and centered vertices

	std::vector<double> angles;
	for (int i = 0; i < vertices.rows(); ++i) {
		double dx = vertices(i, 0) - center(0);
		double dy = vertices(i, 1) - center(1);
		angles.push_back(std::atan2(dy, dx));
	}

	// 4. sort the based on the atan2 value
	auto sortedIdxs = sort_indices<double>(angles);

	std::vector<int> sortedFace(faceIdxs.size());
	for (int i{ 0 }; i < faceIdxs.size(); i++) {
		sortedFace[i] = faceIdxs[sortedIdxs[i]];
	}

	faceIdxs = sortedFace;

	// 5. Sort the 2D projected vertices accordingly
	Eigen::MatrixXd orderedVertices(vertices.rows(), vertices.cols());
	for (int i = 0; i < vertices.rows(); ++i) {
		orderedVertices.row(i) = vertices.row(sortedFace[i]);
	}
	vertices = orderedVertices;
};


void project_vertices_on_plane(
	const std::vector<double>& vertices,
	Eigen::Vector3d& u,
	Eigen::Vector3d& v,
	Eigen::MatrixXd& vertices2D
){
	int vertNr = static_cast<int>(vertices.size() / 3);

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

	// use the center as origin
	Eigen::MatrixXd centeredVertices = verts.rowwise() - center.transpose();

	// create a matrix with columns the unit vector
	Eigen::MatrixXd unitVecs(3, 2);

	unitVecs.col(0) = ue;
	unitVecs.col(1) = ve;

	vertices2D = centeredVertices * unitVecs;
};

void sample_face_polygon(
	const Eigen::MatrixXd& verts2D,
	Eigen::Vector2d& center,
	double& radius1, 
	int resolution) {

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
	for (auto xi : x) {
		for (auto yi : y) {

			Eigen::Vector2d pt{ xi, yi };

			// 4. first check if the pt is inside the polygon defined by the vertices
			if (is_inside_polygon(pt, verts2D)) {

				// 5. find distance from edges, the minimum corresponds to the radius
				// we also use a scale factor to have some space between the hole and the face edges
				double distance = 0.9 * distance_from_edges(pt, verts2D);

				// 6. if the distance is larger than maxDist then this point is the center
				// of the largest circle
				if (distance > maxDist) {
					maxDist = distance;
					center = pt;
				}
			}
		}
	}
	radius1 = maxDist;
};


bool is_inside_polygon(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices) {

	bool isInside = false;

	// 1. loop for each edge
	for (int i{ 0 }; i < vertices.rows(); i++) {

		// 2. get the two points of the edge
		Eigen::VectorXd p1 = vertices.row(i);
		Eigen::VectorXd p2 = vertices.row((i + 1) % vertices.rows());

		// 3. first check if the point is one of the edge vertices
		if (pt == p1 || pt == p2) {
			return false;
		}

		// 4. check if a point is on the edge
		Eigen::Vector2d edge = p2 - p1;
		Eigen::Vector2d toPt = pt - p1;
		double cross = edge.x() * toPt.y() - edge.y() * toPt.x();

		if (std::abs(cross) < 1e-10 &&
			pt.x() >= std::min(p1.x(), p2.x()) &&
			pt.x() <= std::max(p1.x(), p2.x()) &&
			pt.y() >= std::min(p1.y(), p2.y()) &&
			pt.y() <= std::max(p1.y(), p2.y())) {
			return false;  // Consider point on edge as not inside
		}

		// 4. then check if the point is  between ymin and ymax of the edge
		if ((p1.y() > pt.y()) != (p2.x() > pt.y())) {

			// 5. estimate the x intersect. if the point has a smaller x then we have an intersection
			double xIntersect = p1.x() + (pt.y() - p1.y()) * (p2.x() - p1.x()) / (p2.y() - p1.y());
			if (pt.x() < xIntersect) {
				isInside = !isInside;
			}
		}
	}
	return isInside;
};

double distance_from_edges(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices) {

	double minDist = 1e100;

	// 1. loop for each edge
	for (int i{ 0 }; i < vertices.rows(); i++) {

		// 2. get the two points of the edge
		Eigen::VectorXd p1 = vertices.row(i);
		Eigen::VectorXd p2 = vertices.row((i + 1) % vertices.rows());

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

void interpolate_edges(const Eigen::MatrixXd& vertices, Eigen::MatrixXd& interpolatedVertices, const double& edgeSize) {

	// vector to store temp coords
	std::vector<Eigen::Vector2d> tempPoints;

	// push the first point of the first edge
	tempPoints.push_back(vertices.row(0));

	// 1. loop for each edge
	for (int i{ 0 }; i < vertices.rows(); i++) {

		// 2. get the two points of the edge
		Eigen::VectorXd p1 = vertices.row(i);
		Eigen::VectorXd p2 = vertices.row((i + 1) % vertices.rows());

		// 3. find the interpolation point number
		int nr = static_cast<int>(std::ceil((p2 - p1).norm() / edgeSize));

		// 4. get the alpha values
		Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(nr, 0.0, 1.0);

		// 5. Estimate the interpolated values
		//for (int j{ 0 }; j < t.size(); j++) {

		//	double alpha = t[j];

		//	// get the interpolated value
		//	Eigen::Vector2d interp = p1 * (1 - alpha) + p2 * alpha;

		//	tempPoints.push_back(interp);
		//}

		for (int j{ 1 }; j < t.size(); j++) {
			double alpha = t[j];
			// get the interpolated value
			Eigen::Vector2d interp = p1 * (1 - alpha) + p2 * alpha;
			tempPoints.push_back(interp);
		}

	}
	// add also the last point
	tempPoints.push_back(vertices.row(-1));

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

		Eigen::Vector2d pt = { radius * std::cos(theta), radius * std::sin(theta) };

		vertices.row(i) = pt + center;
	}
};

// implementation of circular double linked list


