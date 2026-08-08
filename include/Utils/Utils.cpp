#define _USE_MATH_DEFINES
#include "Utils.h"
#include "Visualize/visualize.h"
#include <Eigen/Dense>
#include <cmath>
#include <math.h>

// Isolated here (rather than in the Eigen/OpenMP-heavy generator TU) so the
// Windows min/max macros never leak into std::min/max or Eigen.
#if defined(_WIN32)
	#ifndef NOMINMAX
		#define NOMINMAX
	#endif
	#include <windows.h>
#endif

uint64_t available_physical_memory_bytes() {
#if defined(_WIN32)
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	if (GlobalMemoryStatusEx(&status)) {
		return static_cast<uint64_t>(status.ullAvailPhys);
	}
	return 0;
#else
	return 0; // unsupported platform: caller treats 0 as "unknown, proceed"
#endif
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

bool ray_aabb_intersection(const Vec3& orig, const Vec3& invDir, const std::array<float, 6> bounds) {
	float t1 = (bounds[0] - orig.x) * invDir.x;
	float t2 = (bounds[1] - orig.x) * invDir.x;
	float tmin = std::min(t1, t2);
	float tmax = std::max(t1, t2);

	t1 = (bounds[2] - orig.y) * invDir.y;
	t2 = (bounds[3] - orig.y) * invDir.y;
	tmin = std::max(tmin, std::min(t1, t2));
	tmax = std::min(tmax, std::max(t1, t2));

	t1 = (bounds[4] - orig.z) * invDir.z;
	t2 = (bounds[5] - orig.z) * invDir.z;
	tmin = std::max(tmin, std::min(t1, t2));
	tmax = std::min(tmax, std::max(t1, t2));

	return tmax >= std::max(0.0f, tmin);
}

bool ray_triangle_intersection(
	const Vec3& orig, const Vec3& dir, const Vec3& v0, const Vec3& v1, const Vec3& v2
) {

	const float EPSILON = 1e-7f;
	Vec3 edge1 = v1 - v0;
	Vec3 edge2 = v2 - v0;
	Vec3 h = dir.cross(edge2);
	float a = edge1.dot(h);

	// If a is near 0, ray is parallel to triangle
	if (a > -EPSILON && a < EPSILON) return false;

	float f = 1.0f / a;
	Vec3 s = orig - v0;
	float u = f * s.dot(h);
	if (u < 0.0f || u > 1.0f) return false;

	Vec3 q = s.cross(edge1);
	float v = f * dir.dot(q);
	if (v < 0.0f || u + v > 1.0f) return false;

	float t = f * edge2.dot(q);
	return t > EPSILON; // valid intersection ahead of the ray

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

float polygon_area_2d(const Eigen::MatrixXd& verts2D) {
	float area = 0.0;
	int n = static_cast<int>(verts2D.rows());
	for (int i = 0; i < n; i++) {
		int next = (i + 1) % n;
		area += (verts2D(i, 0) * verts2D(next, 1)) - (verts2D(next, 0) * verts2D(i, 1));
	}
	return area * 0.5;
};

void fix_orientation(Eigen::MatrixXd& vertices, bool requestCCW) {
	float area = polygon_area_2d(vertices);
	bool isCCW = (area > 0);

	if (isCCW != requestCCW) {
		// Reverse the rows of the matrix
		Eigen::MatrixXd reversed = vertices.colwise().reverse().eval();
		// Note: colwise().reverse() flips rows in Eigen. 
		// If that's confusing, use a manual row swap loop.
		vertices = reversed;
	}
};

void ensure_ccw(
	std::vector<int>& faceIdxs,
	Eigen::MatrixXd& vertices) {

	// 1. Get the centroid of the face
	Eigen::Vector3d center = vertices.colwise().mean();

	// 2. Compute angles for sorting
	std::vector<float> angles;
	for (int i = 0; i < vertices.rows(); ++i) {
		float dx = vertices(i, 0) - center(0);
		float dy = vertices(i, 1) - center(1);
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

/*
* @brief The erosion margin $d^* $ is the largest offset distance such that eroding the cell by $d \le d^* $
* guarantees the target face $f$ remains a non - degenerate 2D polygon.This is critical for
* robustly generating shelled geometries without topological defects.
*
* The margin is determined by the minimum constraint imposed by every other non - face half - space $j$.
*
* @param f The index of the target face for which the erosion margin is being calculated.
* @param normals Vector of all cell face inward unit normals($\mathbf{ n }_j$).
* @param planeB Vector of all cell face plane coefficients($d_j$) such that $\mathbf{ n }_j \cdot \mathbf{ x } \leq d_j$.
* @param U The $3 \times 2$ matrix defining the orthonormal basis of the target face $f$ plane.
* @param x0 An origin point on the plane of face $f$(e.g., its centroid).
* @param verts2D The 2D coordinates of the target face $f$ vertices, projected onto its plane.
* @param delta A small 'pullback' distance used to stabilize the margin and account for the face's curvature/inradius.
* @param eps A tolerance for minimum margin(default $1e - 12$).
* @return float The computed maximum safe erosion margin, $\mathbf{ d^* }$.
*/

float erosion_margin_for_face(
	int f,
	const std::vector<Eigen::Vector3d>& normals,
	const std::vector<double>& planeB,
	const Eigen::Matrix<double, 3, 2>& U,
	const Eigen::Vector3d& x0,
	const Eigen::MatrixXd& verts2D, // N x 2
	float delta,                   // pullback
	float eps
) {
	float dstar = std::numeric_limits<float>::infinity();

	for (int j = 0; j < (int)normals.size(); ++j) {
		
		// skip if the face index f is the same as index j
		if (j == f) continue;

		Eigen::Vector2d m = U.transpose() * normals[j];      // in-plane direction
		float m2 = m.squaredNorm();
		if (m2 < 1e-16) continue;                            // plane parallel to face

		// estimate the inwards 
		float cj = planeB[j] - normals[j].dot(x0);

		// support over polygon in direction m
		float maxDot = -std::numeric_limits<float>::infinity();
		for (int k = 0; k < verts2D.rows(); ++k) {
			float val = m.dot(verts2D.row(k));
			if (val > maxDot) maxDot = val;
		}

		// d*_j = (c_j - supp_P(m)) + delta * ||m||
		float dstar_j = (cj - maxDot) + delta * std::sqrt(m2);
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
	//Eigen::matrixXd verts(vertNr, 3);

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

void back_to_3d(Eigen::MatrixXd& vertices3d, const Eigen::MatrixXd& vertices2d,
	const Eigen::Vector3d& center, const Eigen::Vector3d& u, const Eigen::Vector3d& v) {

	vertices3d.resize(vertices2d.rows(), 3);

	for (int i = 0; i < vertices2d.rows(); ++i) {
		Eigen::Vector3d row = center + u * vertices2d(i, 0) + v * vertices2d(i, 1);
		vertices3d.row(i) = row.transpose();
	}

};

float catmull_rom_get_t(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, float t0, float alpha) {

	float x = std::pow(p2.x() - p1.x(), 2);
	float y = std::pow(p2.y() - p1.y(), 2);

	float ti = std::pow(std::sqrt(x + y), alpha) + t0;

	return ti;
};

void catmull_rom_interpolation(const Eigen::MatrixXd& currentVerts, Eigen::MatrixXd& interpolatedVerts, float alpha) {

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

		float t0 = 0.0;
		float t1 = catmull_rom_get_t(p0, p1, t0, alpha);
		float t2 = catmull_rom_get_t(p1, p2, t1, alpha);
		float t3 = catmull_rom_get_t(p2, p3, t2, alpha);


		for (int i{ 0 }; i < sampleNr - 1; i ++) {

			float ti = t1 + (t2 - t1) * (float(i) / (sampleNr - 1));

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

	// Convert to Eigen::matrixXd
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

//void sample_face_polygon(
//	const Eigen::matrixXd& verts2D,
//	Eigen::Vector2f& center,
//	float& radius1, 
//	int resolution,
//	float scale) {
//
//	// 1. find the min and max for the x and y coordinates of the vertices
//	float xMin = verts2D.col(0).minCoeff();
//	float xMax = verts2D.col(0).maxCoeff();
//
//	float yMin = verts2D.col(1).minCoeff();
//	float yMax = verts2D.col(1).maxCoeff();
//
//	// 2. create x and y points
//	Eigen::VectorXf x = Eigen::VectorXf::LinSpaced(resolution, xMin, xMax);
//	Eigen::VectorXf y = Eigen::VectorXf::LinSpaced(resolution, yMin, yMax);
//
//	// 3. now loop for each point and find the biggest radius estimating distance from edges
//	float maxDist{ 0.0 };
//	int ptIdx = 0;
//	for (auto xi : x) {
//		for (auto yi : y) {
//
//			Eigen::Vector2f pt{ xi, yi };
//
//			// 4. first check if the pt is inside the polygon defined by the vertices
//			if (is_inside_polygon(pt, verts2D)) {
//
//				// 5. find distance from edges, the minimum corresponds to the radius
//				// we also use a scale factor to have some space between the hole and the face edges
//				float distance = scale * min_distance_from_edges(pt, verts2D, ptIdx);
//
//				// 6. if the distance is larger than maxDist then this point is the center
//				// of the largest circle
//				if (distance > maxDist) {
//					maxDist = distance;
//					center = pt;
//				}
//			}
//		}
//		ptIdx++;
//	}
//	radius1 = maxDist;
//};

void sample_face_polygon(
	const Eigen::MatrixXd& verts2D,
	Eigen::Vector2d& center, double& radius1, int resolution, double scale)
{
	// Initialize Bounds
	double xMin = verts2D.col(0).minCoeff();
	double xMax = verts2D.col(0).maxCoeff();
	double yMin = verts2D.col(1).minCoeff();
	double yMax = verts2D.col(1).maxCoeff();

	double bestDist = 0.0;
	Eigen::Vector2d bestCenter = verts2D.colwise().mean(); // Start guess at centroid

	// --- PASS 1: COARSE GRID (Global Search) ---
	// Scan the whole bounding box
	Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(resolution, xMin, xMax);
	Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(resolution, yMin, yMax);
	int ptIdx = 0;
	for (auto xi : x) {
		for (auto yi : y) {
			Eigen::Vector2d pt{ xi, yi };

			// Optimization: Bounding Box pre-check 
			// (If the point is too close to the box edge, it can't be better than what we found)
			if (xi - xMin < bestDist || xMax - xi < bestDist ||
				yi - yMin < bestDist || yMax - yi < bestDist) continue;

			if (is_inside_polygon(pt, verts2D)) {
				// Remove 'scale' and 'ptIdx' from the raw math
				double d = min_distance_from_edges(pt, verts2D);
				if (d > bestDist) {
					bestDist = d;
					bestCenter = pt;
				}
			}
		}
	}

	// --- PASS 2: FINE GRID (Local Refinement) ---
	// Zoom in on the best spot found in Pass 1.
	// The new search box size is roughly 2x the step size of the previous grid.
	double rangeX = (xMax - xMin) / (resolution - 1);
	double rangeY = (yMax - yMin) / (resolution - 1);

	double localXMin = std::max(xMin, bestCenter.x() - rangeX);
	double localXMax = std::min(xMax, bestCenter.x() + rangeX);
	double localYMin = std::max(yMin, bestCenter.y() - rangeY);
	double localYMax = std::min(yMax, bestCenter.y() + rangeY);

	Eigen::VectorXd xFine = Eigen::VectorXd::LinSpaced(resolution, localXMin, localXMax);
	Eigen::VectorXd yFine = Eigen::VectorXd::LinSpaced(resolution, localYMin, localYMax);
	ptIdx = 0;
	for (auto xi : xFine) {
		for (auto yi : yFine) {
			Eigen::Vector2d pt{ xi, yi };
			if (is_inside_polygon(pt, verts2D)) {
				float d = min_distance_from_edges(pt, verts2D, ptIdx);
				if (d > bestDist) {
					bestDist = d;
					bestCenter = pt;
				}
			}
		}
		ptIdx++;
	}

	// Final Output
	center = bestCenter;
	radius1 = bestDist * scale; // Apply scaling at the very end
};

bool is_point_in_triangle(const Eigen::Vector2d& p, const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c) {
	auto sign = [](const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& p3) {
		return (p1.x() - p3.x()) * (p2.y() - p3.y()) - (p2.x() - p3.x()) * (p1.y() - p3.y());
	};

	float d1 = sign(p, a, b);
	float d2 = sign(p, b, c);
	float d3 = sign(p, c, a);

	bool has_neg = (d1 < -1e-12) || (d2 < -1e-12) || (d3 < -1e-12);
	bool has_pos = (d1 > 1e-12) || (d2 > 1e-12) || (d3 > 1e-12);

	return !(has_neg && has_pos);
};

bool is_inside_polygon(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices) {

	bool isInside = false;
	const float eps = 1e-10;

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
		float cross = edge.x() * toPt.y() - edge.y() * toPt.x();

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
				float xIntersect = p1.x() + (pt.y() - p1.y()) * (p2.x() - p1.x()) / (p2.y() - p1.y());
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
	float dot00 = v0.dot(v0);
	float dot01 = v0.dot(v1);
	float dot02 = v0.dot(v2);
	float dot11 = v1.dot(v1);
	float dot12 = v1.dot(v2);

	// barycentric coordinates
	float denom = (dot00 * dot11 - dot01 * dot01);

	if (std::abs(denom) - 1e-20) {
		return true;
	}
	
	float invDenom = 1.0 / denom;

	//u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	//v = (dot00 * dot12 - dot01 * dot02) * invDenom;

};

float min_distance_from_edges(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices, int excludeIdx) {

	float minDist = std::numeric_limits<float>::infinity();

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
		float t = edge.dot(toPt) / edge.squaredNorm();

		// 6. clamp it to be on the segment
		t = std::clamp(t, 0.0f, 1.0f); 

		// 7. the vector projection is 
		Eigen::Vector2d projection = p1 + t * edge;
		
		// 8. the distance is the norm of the vertical vector
		float dist = (pt - projection).norm();

		// 9. check if the dist is smaller than minDists
		minDist = std::min(minDist, dist);
	}

	return minDist;

};

float min_distance_from_edges(const Eigen::Vector2d& pt, const Eigen::MatrixXd& verts) {
	float minD = std::numeric_limits<float>::max();

	int numVerts = verts.rows();
	for (int i = 0; i < numVerts; ++i) {
		// Get edge vertices (A -> B)
		Eigen::Vector2d A = verts.row(i);
		Eigen::Vector2d B = verts.row((i + 1) % numVerts); // Wrap around for last edge

		// Project point onto line segment AB
		Eigen::Vector2d AB = B - A;
		Eigen::Vector2d AP = pt - A;

		// Calculate factor 't' along the segment
		float t = AP.dot(AB) / AB.squaredNorm();

		// Clamp 't' to the segment [0, 1]
		// This handles the case where the closest point is a vertex, not the edge itself
		t = std::clamp(t, 0.0f, 1.0f);

		// Find closest point on segment
		Eigen::Vector2d closest = A + t * AB;

		// Distance
		float d = (pt - closest).norm();

		if (d < minD) {
			minD = d;
		}
	}
	return minD;
}

float max_distance_from_edges(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices, int excludeIdx) {

	float maxDist = std::numeric_limits<float>::infinity();

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
		float t = edge.dot(toPt) / edge.squaredNorm();

		// 6. clamp it to be on the segment
		t = std::clamp(t, 0.0f, 1.0f);

		// 7. the vector projection is 
		Eigen::Vector2d projection = p1 + t * edge;

		// 8. the distance is the norm of the vertical vector
		float dist = (pt - projection).norm();

		// 9. check if the dist is smaller than minDists
		maxDist = std::max(maxDist, dist);
	}

	return maxDist;

};

float polygon_min_width(const Eigen::MatrixXd& vertices) {

	float minWidth = std::numeric_limits<float>::infinity();

	for (int i{ 0 }; i < vertices.rows(); i++) {

		float width = min_distance_from_edges(vertices.row(i), vertices, i);

		minWidth = std::min(minWidth, width);
	}

	return minWidth;
};

float polygon_max_width(const Eigen::MatrixXd& vertices) {

	float maxWidth = -1.0;

	for (int i{ 0 }; i < vertices.rows(); i++) {

		float width = max_distance_from_edges(vertices.row(i), vertices, i);

		maxWidth = std::max(maxWidth, width);
	}

	return maxWidth;
};

float polygon_average_edge_length(const Eigen::MatrixXd& vertices) {

	float avgLength{ 0.0 };

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

float polygon_area(const Eigen::MatrixXd& verts) {
	float area = 0.0;
	int N = verts.rows();
	for (int i = 0; i < N; ++i) {
		int j = (i + 1) % N;
		area += verts(i, 0) * verts(j, 1);
		area -= verts(i, 1) * verts(j, 0);
	}
	return 0.5 * std::abs(area);
}

float polygon_perimeter(const Eigen::MatrixXd& verts) {
	float perim = 0.0;
	for (int i = 0; i < verts.rows(); ++i) {
		perim += (verts.row(i) - verts.row((i + 1) % verts.rows())).norm();
	}
	return perim;
};

void interpolate_edges(
	const Eigen::MatrixXd& vertices,
	Eigen::MatrixXd& interpolatedVertices,
	std::vector<int>& newLocalFace,
	const float& edgeSize) {

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
			float alpha = t[j];
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

void hole_points(const float& radius, const Eigen::Vector2d& center, const int& ptNr, Eigen::MatrixXd& vertices) {

	vertices.resize(ptNr, 2);

	// 1. create the hole vertices
	for (int i{ 0 }; i < ptNr; i++) {

		float theta = i * 2 *  M_PI / ptNr;

		Eigen::Vector2d pt = { radius * std::cos(theta), -radius * std::sin(theta) };

		vertices(i, 0) = pt(0) + center(0);
		vertices(i, 1) = pt(1) + center(1);
	}
};

bool vertex_locally_convex(const Eigen::Vector2d& v1, const Eigen::Vector2d& v2, const Eigen::Vector2d& v3) {

	Eigen::Vector2d edge1 = v2 - v1;
	Eigen::Vector2d edge2 = v3 - v1;

	// estimate the angle
	float dot = edge1.dot(edge2);
	float det = edge1(0) * edge2(1) - edge1(1) * edge2(0);

	float angle = std::atan2(det, dot);

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
		bool earFound = false;
		Node* startNode = node; // Track where we started this pass

		do {
			int i = node->prev->data, j = node->data, k = node->next->data;
			Eigen::Vector2d vPrev = vertices.row(i), vCurr = vertices.row(j), vNext = vertices.row(k);

			if (vertex_locally_convex(vPrev, vCurr, vNext)) {
				bool isEar = true;
				Node* testNode = node->next->next;

				while (testNode != node->prev) {
					// IMPORTANT: Skip checking if the test vertex is one of the triangle vertices
					int tid = testNode->data;
					if (tid != i && tid != j && tid != k) {
						if (is_point_in_triangle(vertices.row(tid), vPrev, vCurr, vNext)) {
							isEar = false;
							break;
						}
					}
					testNode = testNode->next;
				}

				if (isEar) {
					cells.push_back({ i, j, k });
					Node* toDelete = node;
					node = node->next; // Move to next before deleting
					cdll.remove(toDelete);
					earFound = true;
					break; // Found an ear, exit the 'do-while' and restart
				}
			}
			node = node->next;
		} while (node != startNode);

		if (!earFound) {
			// EMERGENCY FALLBACK: If no ear is found, the polygon is likely degenerate 
			// or self-intersecting. Clip the first convex vertex we find regardless of 
			// points inside, or exit to avoid infinite loop.
			return false;
		}
	}

	return true;
};

bool polygon_is_open(
	const std::vector<Eigen::Vector3d>& normals,
	const std::vector<float>& d,
	const Eigen::MatrixXd& faceVertices,
	int fIdx,
	float erosion,
	float eps
) {
	// loop for all face vertices
	//float dstar = std::numeric_limits<float>::infinity();

	for (int i{ 0 }; i < faceVertices.rows(); i++) {

		// get the vertex coords
		Eigen::Vector3d p = faceVertices.row(i);

		// now loop for all other faces
		for (int j{ 0 }; j < static_cast<int>(normals.size()); j++) {
			if (j == fIdx) continue;

			// now check if the vertex will be above or below the eroded plane
			// we will accept if normal_j * p > d_j - erosion + eps
			float diff = d[j] - normals[j].dot(p);
			if (normals[j].dot(p) < d[j] - erosion + eps) return false;
			//dstar = std::min(dstar, diff);
		}
	}
	//return (dstar > eps) ? dstar : 0.0;
	return true;
	//return (dstar > eps) ? dstar : 0.0;
};


/* @brief function for half space clipping using the Sutherland - Hodgman Algorithm
* 
* @param verts2D: A $N \times 2$ Eigen::matrixXd where each row represents a vertex $(x, y)$ of the polygon. It contains
* the 2D vertices	of the polygon (projected on the plane defined by the polygon)
*
* @param normal (Eigen::Vector2f): A 2d vector that represents the normal of the clipping plane
* 
* @param d (float): The coefficient of the clipping plane equation
* 
* @return 
*	result (Eigen::matrixXd): An matrix where each row corresponds to a vertex of the clipped polygon 
*
*/
Eigen::MatrixXd clip_polygon(
	const Eigen::MatrixXd& verts2D,
	const Eigen::Vector2d& normal,
	const float d
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

/* @brief Estimates the intersection of 2D lines
* 
* @param: p(Eigen::Vector2f) a point in lin
*
*/
Eigen::Vector2d intersection_of_lines(
	const Eigen::Vector2d& p,
	const Eigen::Vector2d& q,
	const Eigen::Vector2d& normal,
	float d
) {

	float denom = normal.dot(q - p);
	if (std::abs(denom) > 1e-15) return q;

	float t = (d - p.dot(normal)) / denom;

	return p + t * (q - p);
};

/* @brief Checks if the polygon clipping leads to a valid polygon
*	 a polygon is valid if it has more than 2 vertices (at least a triangle)
*	
* @param verts2D: A $N \times 2$ Eigen::matrixXd where each row represents a vertex $(x, y)$ of the polygon. It contains
* the 2D vertices	of the polygon (projected on the plane defined by the polygon)
* 
* @return true if the clipping leads to a valid polygon
*/ 
bool clipping_is_valid(
	const Eigen::MatrixXd& verts2D,
	const std::vector<HalfSpace>& hspaces,
	float delta,
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

/*@brief function to create half spaces for a convex 2D polygon that define its interior. 
* 
* @param verts2D: A $N \times 2$ Eigen::matrixXd where each row represents a vertex $(x, y)$ of the polygon. It contains 
* the 2D vertices	of the polygon (projected on the plane defined by the polygon)
* 
* @return A vector of HalfSpace structures. Each structure containes:
*	- `normal` (Eigen::Vector2f): The unit vector pointing inward from the edge.
*	- `d` (float): The plane coefficient such that the half-space is defined by $\mathbf{n} \cdot \mathbf{x} \leq d$.
* 
* @note The function assumes the polygon is **closed** (i.e., the last vertex connects to the first).
* @note The function uses the polygon's **centroid** to robustly determine the inward-pointing direction of the normal.
* @note Edges with length less than $1 \times 10^{-16}$ (degenerate edges) are skipped.
*/
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

		// if the normal points inwards flip it
		if (normal.dot(centroid) > normal.dot(p)) normal = -normal;

		float len = normal.norm();
		if (len < 1e-16) continue;
		normal /= len;

		// find the plane coeff
		float d = normal.dot(p);

		hSpaces.push_back({ normal, d });
	}

	return hSpaces;
};

/*@brief Estimates the inradius of a 2D convex polygon, using the halfspaces
* defined by the polygon's faces. 
* 
* @param: verts2D (Eigen::matrixXd&) an Eigen Matrix that contains the 2D vertices
*	of the polygon (projected on the plane defined by the polygon)
* 
* @returns 
*	-low (float) the polygon's inradius
*/
float polygon_inradius(const Eigen::MatrixXd& verts2D) {

	// check if a face is degenerated
	if (verts2D.rows() < 3 || verts2D.cols() != 2) return 0.0;

	// first gather the halfspaces
	std::vector<HalfSpace> hSpaces = get_poly_half_spaces(verts2D);
	
	if (hSpaces.empty()) return 0.0;

	// bisection limits
	float low = 0.0;
	float upper = 1e-6;

	// find an infeasible upper bound
	Eigen::MatrixXd dummy;
	while (clipping_is_valid(verts2D, hSpaces, upper, &dummy)) {
		low = upper;
		upper *= 2.0;
		if (upper > 1e6) break; // safety
	}

	// clip to zero
	upper = std::max(0.0f, upper);

	// keep also the best poly to provide the centroid as a center
	Eigen::MatrixXd best;

	// use a number of iterations to find the inradius
	for (int i{ 0 }; i < 25; i++) {
		// take midpoint
		float mid = 0.5 * (low + upper);

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

float get_ground_truth_margin(const Eigen::MatrixXd& verts2D) {

	float low = 0.0;
	float high = polygon_inradius(verts2D) * 1.5; // Search a bit beyond inradius

	std::vector<HalfSpace> hSpaces = get_poly_half_spaces(verts2D);

	// Binary Search
	for (int k = 0; k < 20; ++k) {
		float mid = 0.5 * (low + high);
		Eigen::MatrixXd dummyPoly;
		if (clipping_is_valid(verts2D, hSpaces, mid, &dummyPoly)) {
			low = mid;
		}
		else {
			high = mid;
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
			edgeNames += " '" + std::to_string(edge.first) + " weight " + std::to_string(edge.second.ratio) + "' ";
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
			if (w.ratio < weight && vertex1 < vertex2) {
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
			if (w.ratio > weight) {
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

int Graph::find_longest_network_new() {
	// Use a flat vector for visited status (faster than std::set)
	// Assuming node IDs are 0 to numVertices-1
	// If your IDs are sparse/random, keep using std::set or std::unordered_set

	// Find max node ID to size the vector
	int maxId = 0;
	for (const auto& pair : adjList) {
		if (pair.first > maxId) maxId = pair.first;
	}

	// Visited array (0 = unvisited, 1 = visited)
	std::vector<uint8_t> visited(maxId + 1, 0);

	int maxComponentSize = 0;

	// Iterate over all nodes in the adjacency list
	for (const auto& [startNode, neighbors] : adjList) {

		if (visited[startNode]) continue;

		// --- ITERATIVE DFS (Heap Stack) ---
		int currentComponentSize = 0;
		std::vector<int> stack;

		stack.push_back(startNode);
		visited[startNode] = 1;

		while (!stack.empty()) {
			int u = stack.back();
			stack.pop_back();
			currentComponentSize++;

			// Check neighbors
			if (adjList.count(u)) {
				for (const auto& edge : adjList.at(u)) {
					int v = edge.first;

					// Safety check if v is within bounds (if logic relies on external IDs)
					if (v < visited.size() && !visited[v]) {
						visited[v] = 1;
						stack.push_back(v);
					}
				}
			}
		}
		// ----------------------------------

		if (currentComponentSize > maxComponentSize) {
			maxComponentSize = currentComponentSize;
		}
	}

	return maxComponentSize;
}

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
		return { -1.0f, -1.0f }; // edge v1 -> v2 not in graph
	}

	// if ok return the weight
	return it_v2->second;
};

//std::unordered_map<int, std::unordered_map<int, float>> Graph::get_adj_list() {
std::unordered_map<int, std::unordered_map<int, EdgeData>> Graph::get_adj_list() {
	return adjList;
};

static float Median(std::vector<float>& v) {
	if (v.empty()) return 0.0;
	size_t mid = v.size() / 2;
	std::nth_element(v.begin(), v.begin() + mid, v.end());
	float m = v[mid];
	if (v.size() % 2 == 0) {
		auto it = std::max_element(v.begin(), v.begin() + mid);
		m = 0.5 * (m + *it);
	}
	return m;
};

Vec3 get_barycentric_point(const Vec3& pt, const Triangle& tri) {
	Vec3 v0 = tri.v2 - tri.v1;
	Vec3 v1 = tri.v3 - tri.v1;
	Vec3 v2 = pt - tri.v1;

	float d00 = v0.dot(v0);
	float d01 = v0.dot(v1);
	float d11 = v1.dot(v1);
	float d20 = v2.dot(v0);
	float d21 = v2.dot(v1);

	float denom = d00 * d11 - d01 * d01;

	// Fallback for degenerate triangles to avoid division by zero
	if (std::abs(denom) < 1e-8f) {
		return { 1.0f, 0.0f, 0.0f };
	}

	Vec3 barycenter;
	barycenter.y = (d11 * d20 - d01 * d21) / denom; // v
	barycenter.z = (d00 * d21 - d01 * d20) / denom; // w
	barycenter.x = 1.0f - barycenter.y - barycenter.z; // u

	return barycenter;
};

Vec3 closest_triangle_point(const Vec3& pt, const Triangle& tri) {

	Vec3 ab = tri.v2 - tri.v1;
	Vec3 ac = tri.v3 - tri.v1;
	Vec3 ap = pt - tri.v1;

	float d1 = ab.dot(ap);
	float d2 = ac.dot(ap);
	if (d1 <= 0.0f && d2 <= 0.0f) return tri.v1;

	Vec3 bp = pt - tri.v2;
	float d3 = ab.dot(bp);
	float d4 = ac.dot(bp);
	if (d3 >= 0.0f && d4 <= d3) return tri.v2;

	float vc = d1 * d4 - d3 * d2;
	if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
		float v = d1 / (d1 - d3);
		return tri.v1 + ab * v; // Edge AB
	}

	Vec3 cp = pt - tri.v3;
	float d5 = ab.dot(cp);
	float d6 = ac.dot(cp);
	if (d6 >= 0.0f && d5 <= d6) return tri.v3; // Vertex C

	float vb = d5 * d2 - d1 * d6;
	if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
		float w = d2 / (d2 - d6);
		return tri.v1 + ac * w; // Edge AC
	}

	float va = d3 * d6 - d5 * d4;
	if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
		float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		return tri.v2 + (tri.v3 - tri.v2) * w; // Edge BC
	}

	float denom = 1.0f / (va + vb + vc);
	float v = vb * denom;
	float w = vc * denom;
	return tri.v1 + ab * v + ac * w; // Inside face

};

Eigen::Matrix3f rotation_axis_angle(const Vec3& dir, float angle) {

	if (angle == 0.0f) return Eigen::Matrix3f::Identity();

	// enusre normalizednDir.x
	Vec3 nDir = dir.normalized();
	Eigen::Vector3f k = { nDir.x, nDir.y, nDir.z };

	// convet angle to radians
	float angleRad = angle * M_PI / 180.0f;

	return Eigen::AngleAxisf(angleRad, k).toRotationMatrix();
};

Eigen::Matrix3f rotation_from_direction(
	const Vec3& dir, float stretchX, float stretchY, float stretchZ) {

	// ensure direction is normalized
	Vec3 nDir = dir.normalized();
	Eigen::Vector3f targetDir = { nDir.x, nDir.y, nDir.z };
	Eigen::Vector3f defaultAxis;
	
	if (stretchZ >= stretchX && stretchZ >= stretchY) {
		defaultAxis = Eigen::Vector3f::UnitZ(); // Z is the main stretch
	}
	else if (stretchY >= stretchX && stretchY >= stretchZ) {
		defaultAxis = Eigen::Vector3f::UnitY(); // Y is the main stretch
	}
	else {
		defaultAxis = Eigen::Vector3f::UnitX(); // X is the main stretch
	}

	Eigen::Quaternionf quat = Eigen::Quaternionf::FromTwoVectors(
		targetDir, defaultAxis);

	return quat.toRotationMatrix();
};

Eigen::Matrix3f rotation_from_direction_roll(
	const Vec3& dir, float rollAngle
) {
    
	// convert to radians
	float rollRadians = rollAngle * (3.14159265359f / 180.0f);

    Eigen::Vector3f targetDir(dir.x, dir.y, dir.z);
    
    // Fallback to X-axis if the user inputs a zero-vector (0,0,0)
    if (targetDir.norm() < 1e-6f) {
        targetDir = Eigen::Vector3f::UnitX();
    } else {
        targetDir.normalize();
    }

    // 1. Point the Local X-axis toward the Target Direction
    Eigen::Quaternionf alignQuat = Eigen::Quaternionf::FromTwoVectors(targetDir, Eigen::Vector3f::UnitX());
    
    // 2. Roll around the Local X-axis
    Eigen::Quaternionf rollQuat(Eigen::AngleAxisf(
		rollRadians, Eigen::Vector3f::UnitX()));
    
    // 3. Combine (Roll in local space, then point to global space)
    Eigen::Quaternionf finalQuat = rollQuat * alignQuat;
    
    return finalQuat.toRotationMatrix();
}

