#ifndef UTILS_H // include guard
#define UTILS_H

# define _USE_MATH_DEFINES

#include "Eigen/Dense"
#include <string.h>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <ImGuiFileDialog/ImGuiFileDialog.h>
#include "Math/Vec.h"

enum ObjectType {
	BoxContainerType,
	CylinderContainerType,
	AbstractContainerType,
	RandomGeneratorType,
	UniformGeneratorType,
	VariedGeneratorType,
	ScaffoldType,
	NoneType
};

struct Bounds {
	double xMin, xMax, yMin, yMax, zMin, zMax;
	Vec3 center;
};

struct PolygonWidth {
	float minWidth, maxWidth;
};

struct Triangle {
	Vec3 v1, v2, v3;
	unsigned int i1, i2, i3;
	Vec3 normal;
};

struct Aabb {
	Vec3 pMin{};
	Vec3 pMax{};
};

struct AStarNode {
	int idx;
	float fScore;

	AStarNode(int idx, float fScore) : idx(idx), fScore(fScore) {};
	// we want the priority queue to order the nodes by their fScore, so we need to overload the < operator
	bool operator>(const AStarNode& other) const {
		return fScore > other.fScore; // we want the node with the lowest fScore to be at the top of the priority queue, so we use > instead of <
	}
};

// create a class to control seed generator using the strategy pattern
class SeedGeneratorInterface {
public:
	virtual ~SeedGeneratorInterface() = default;
	virtual void run_generate_seeds() = 0;
	virtual void run_get_seeds(std::vector<std::array<double, 3>>& outSeeds) const = 0;
	virtual void run_generate_seeds(
		const std::array<double, 3>&, const std::array<double, 3>&, double) = 0;
};

bool ray_intersection(
	const Eigen::Vector3d& p,
	const Eigen::Vector3d& v1,
	const Eigen::Vector3d& v2,
	const Eigen::Vector3d& v3,
	const Eigen::Vector3d& dir
);

bool ray_intersection(
	const Eigen::Vector3d& p,
	const Eigen::Vector3d& v1,
	const Eigen::Vector3d& v2,
	const Eigen::Vector3d& v3,
	const Eigen::Vector3d& dir,
	Eigen::Vector3d& intersection
);

template <typename T>
bool is_inside_box(const std::array<T, 3>& pt, const std::array<T, 6>& bounds) {

	if (pt[0] < bounds[0] || pt[0] > bounds[1] ||
		pt[1] < bounds[2] || pt[1] > bounds[3] ||
		pt[2] < bounds[4] || pt[2] > bounds[5]) {
		return false;
	}
	else {
		return true;
	}
}


Eigen::Vector3d unit_axis_from_dir(int dir);

bool is_inside_box(const std::array<double, 3>& p, const Bounds& b);

bool is_inside_cylinder(
	const std::array<double, 3>& querryPt,
	const Eigen::Vector3d& basePt,
	const Eigen::Vector3d& axis,
	double radius, double height
);

double distance(const std::array<double, 3>& arr1, const std::array<double, 3>& arr2);

void select_file_button(const char* title, const std::string path, const char* dialogName, const char* fileExt);

float polygon_area_2d(const Eigen::MatrixXd& verts2D);

void fix_orientation(Eigen::MatrixXd& vertices, bool requestCCW);

void ensure_ccw(
	std::vector<int>& faceIdxs,
	Eigen::MatrixXd& vertices);

Eigen::MatrixXd vector_to_matrix3D(const std::vector<double>& verts);

//void face_metrics(float& area, float& perimeter, float& minWidth);

void project_vertices_on_plane(
	const std::vector<double>& vertices,
	Eigen::Vector3d& origin,
	Eigen::Vector3d& u,
	Eigen::Vector3d& v,
	Eigen::MatrixXd& vertices2D
);

void project_vertices_on_plane(
	Eigen::MatrixXd& verts,
	Eigen::Vector3d& origin,
	Eigen::Vector3d& u,
	Eigen::Vector3d& v,
	Eigen::MatrixXd& vertices2D
);

// create an argsort function
// taken by https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes

template <typename T>
std::vector<size_t> sort_indices(const std::vector<T>& v) {

	// initialize original index locations
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	// using std::stable_sort instead of std::sort
	// to avoid unnecessary index re-orderings
	// when v contains elements of equal values 
	std::stable_sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx;
};

void sample_face_polygon(
	const Eigen::MatrixXd& verts2D,
	Eigen::Vector2d& center, double& radius1, int resolution, double scale);

bool is_point_in_triangle(const Eigen::Vector2d& p, const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c);

bool is_inside_polygon(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices);

bool is_inside_triangle(const Eigen::Vector2d& vPrev, const Eigen::Vector2d& vCurr, const Eigen::Vector2d& vNext, const Eigen::Vector2d& vTest);

float min_distance_from_edges(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices, int excludedIdx);

float min_distance_from_edges(const Eigen::Vector2d& pt, const Eigen::MatrixXd& verts);

float max_distance_from_edges(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices, int excludedIdx);

float polygon_min_width(const Eigen::MatrixXd& vertices);

float polygon_max_width(const Eigen::MatrixXd& vertices);

float polygon_average_edge_length(const Eigen::MatrixXd& vertices);

float polygon_area(const Eigen::MatrixXd& verts);

float polygon_perimeter(const Eigen::MatrixXd& verts);

void interpolate_edges(
	const Eigen::MatrixXd& vertices,
	Eigen::MatrixXd& interpolatedVertices,
	std::vector<int>& newLocalFace,
	const float& edgeSize);

void hole_points(const float& radius, const Eigen::Vector2d center, const int& ptNr, Eigen::MatrixXd& vertices);

bool ear_clipping(const Eigen::MatrixXd& vertices, const std::vector<int>& idxs, std::vector<std::vector<int>>& cells);

bool vertex_locally_convex(const Eigen::Vector2d& v1, const Eigen::Vector2d& v2, const Eigen::Vector2d& v3);

void back_to_3d(
	Eigen::MatrixXd& vertices3f,
	const Eigen::MatrixXd& vertices2d, const Eigen::Vector3d& center, const Eigen::Vector3d& u, const Eigen::Vector3d& v);
// implementation of circular float linked list

void catmull_rom_interpolation(const Eigen::MatrixXd& currentVerts, Eigen::MatrixXd& interpolatedVerts, float alpha=0.5);

float catmull_rom_get_t(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, float t0, float alpha = 0.5);

void chaikin_subdivision(const Eigen::MatrixXd& currentVerts, Eigen::MatrixXd& interpolatedVerts, int num = 5);

// @brief this functions returns a struct that holds the minimum and maximum width of a polygon face
// @brief for each edge we project all vertices to its normal. the width is projMax - projMin.  
PolygonWidth get_polygon_width(const Eigen::MatrixXd& verts2D);

// @brief function to decide if a polygon will be totally closed if we erode all edges by a factor
bool polygon_is_open(
	const std::vector<Eigen::Vector3d>& normals,
	const std::vector<float>& d,
	const Eigen::MatrixXd& faceVertices,
	int fIdx,
	float erosion,
	float eps = 1e-6
);

float erosion_margin_for_face(
	int f,
	const std::vector<Eigen::Vector3d>& normals,
	const std::vector<double>& planeB,
	const Eigen::Matrix<double, 3, 2>& U,
	const Eigen::Vector3d& x0,
	const Eigen::MatrixXd& verts2D, // N x 2
	float delta,                   // pullback
	float eps = 1e-9
);

//@brief struct for holding the half space properties normal, and offset
struct HalfSpace {
	Eigen::Vector2d n;
	float d;
};

// @function for half space clipping using the Sutherland - Hodgman Algorithm
Eigen::MatrixXd clip_polygon(
	const Eigen::MatrixXd& verts2D,
	const Eigen::Vector2d& normal,
	const float d
);

// @function to find intersection of lines
Eigen::Vector2d intersection_of_lines(
	const Eigen::Vector2d& p,
	const Eigen::Vector2d& q,
	const Eigen::Vector2d& normal,
	float d
);

// @function to test if a clipping by all shifted inward edges by a factor delta returns a polygon
bool clipping_is_valid(
	const Eigen::MatrixXd& verts2D,
	const std::vector<HalfSpace>& hspaces,
	float delta,
	Eigen::MatrixXd* outPoly = nullptr
);

// @function to get the half spaces for a polygon
std::vector<HalfSpace> get_poly_half_spaces(const Eigen::MatrixXd& verts2D);

// @brief function to find the max inset radius of a polygon using bisection
float polygon_inradius(const Eigen::MatrixXd& verts2D);

float get_ground_truth_margin(const Eigen::MatrixXd& verts2D);

Vec3 get_barycentric_point(const Vec3& pt, const Triangle& tri);

Vec3 closest_triangle_point(const Vec3& pt, const Triangle& tri);

// --------------------------------------------------------------------------------------------------------------
// create a node struct
class Node {

public:

	~Node() {};
	Node(int& data, Node* next, Node* prev) : data(data), next(next), prev(prev) {};

	// data of the node
	int data;

	// pointer to next and prev
	Node* next;
	Node* prev;
};

class Cdll {
public:

	Cdll() : head(nullptr) {};

	// append at the end of the list
	Node* append(int data) {

		Node* newNode = new Node(data, nullptr, nullptr);

		length++;

		// check if the dll is empty
		if (!head) {
			head = newNode;
			newNode->prev = newNode;
			newNode->next = newNode;
		}
		else {
			// get the final node, since we have a circular dll this the previous of the tail
			Node* temp = head->prev;
			temp->next = newNode;
			newNode->prev = temp;
			newNode->next = head;
			head->prev = newNode;
		}

		return newNode;
	};

	// a function to display the elements
	void display() {

		Node* node = head;

		while (node->next != head) {

			std::cout << node->data << " -- ";

			node = node->next;

		}

		// print also the last
		std::cout << node->data << " < -- > " << std::endl;
		std::cout << node->next->data << std::endl;
	};

	// a function to delete an element
	void remove(Node* node) {

		// first check if we have only a node
		if (node->next == node) {
			head = nullptr;
			head->prev = nullptr;
		}
		else {
			Node* prev = node->prev;
			Node* next = node->next;

			prev->next = next;
			next->prev = prev;

			if (node == head) {
				head = node->next;
			}
		}
		length -= 1;

	};

	// a function to rotate the dll by a position
	void rotate(int pos) {

		if (head == nullptr) {
			return;
		}

		else {
			Node* current = head;

			for (int i{ 0 }; i < pos; i++) {

				current = current->next;

				head = current;
				head->prev = current->prev;

			} 
		}

	}

	// a function to get the data from the nodes
	void get_data(std::vector<int>& data) {

		Node* node = head;

		while (true){

			data.push_back(node->data);
			node = node->next;
			if (node->next == head->next) {
				break;
			}
		}
		data.push_back(head->data);
	}

	void get_length(int& length) {

		int len{ 0 };

		Node* node = head;

		while (node->next != head) {

			len++;

			node = node->next;

		}

		// add one to get also the final node
		length = len + 1;
	}

	int length{ 0 };

	// an empty first node as a head
	Node* head;
};

struct EdgeData {
	//float dstar;   // margin (world units)
	//float delta;   // in-plane pullback actually used when computing d*
	float ratio;
	float targetArea;

	//float inradius;
	//std::array<float, 2> center{0.0, 0.0};
};

class Graph {
private:
	std::unordered_map<int, std::unordered_map<int, EdgeData>> adjList;

public:
	Graph() {};	
	~Graph() {};

	// copy constructor
	Graph(const Graph& otherGraph) = default;

	bool add_vertex(const int idx);

	//bool add_edge(const int vertex1, const int vertex2, const float weight);
	bool add_edge(const int vertex1, const int vertex2, const EdgeData& edge);

	bool remove_edge(const int vertex1, const int vertex2);

	bool remove_vertex(const int vertex1);

	void remove_edges_below(const float weight);

	void remove_edges_above(const float weight);

	void dfs_traversal(int startNode, std::set<int>& visited);

	int find_longest_network();

	int get_vertex_count();

	int find_longest_network_new();

	//float get_edge_width(const int vertex1, const int vertex2);
	EdgeData Graph::get_edge_width(const int vertex1, const int vertex2);

	std::unordered_map<int, std::unordered_map<int, EdgeData>> get_adj_list();

	void print();

};

// helper: median
static float Median(std::vector<float>& v);

#endif