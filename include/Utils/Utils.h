#ifndef UTILS_H // include guard
#define UTILS_H

# define _USE_MATH_DEFINES

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include "Eigen/Dense"
#include <string.h>
#include <vtkGlyph3D.h>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <ImGuiFileDialog/ImGuiFileDialog.h>

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

bool is_inside_mesh(
	const vtkSmartPointer<vtkPolyData>& mesh,
	const Eigen::Vector3d& point);

bool is_inside_mesh(
	const vtkSmartPointer<vtkPolyData>& mesh,
	const Eigen::Vector3d& point,
	const Eigen::Vector3d& rayDir,
	Eigen::Vector3d& intersection);

bool is_inside_box(
	const std::array<double, 3>& pt,
	const std::array<float, 6>& bounds);

bool is_inside_cylinder(

	const std::array<double, 3>& querryPt,
	const Eigen::Vector3d& basePt,
	const Eigen::Vector3d& axis,
	double radius, double height
);

vtkSmartPointer<vtkActor> render_points(
	const int& particles,
	const std::vector<std::array<double, 3>>& currSeeds, 
	const std::string color);

void make_glyphs(vtkPolyData* src, double size, vtkGlyph3D* glyph);

vtkSmartPointer<vtkActor> render_cell(const std::vector<double>& cellVerts);

vtkSmartPointer<vtkActor> render_plane(
	const double& n1, const double& n2, const double& n3,
	const double& p1, const double& p2, const double& p3,
	const double& s1, const double& s2
);

double distance(const std::array<double, 3>& arr1, const std::array<double, 3>& arr2);

void selectFileButton(const char* title, const std::string path, const char* dialogName, const char* fileExt);
#endif

void ensure_ccw(
	std::vector<int>& faceIdxs,
	Eigen::MatrixXd& vertices);

void project_vertices_on_plane(
	const std::vector<double>& vertices,
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

bool is_inside_polygon(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices);

bool is_inside_triangle(const Eigen::Vector2d& vPrev, const Eigen::Vector2d& vCurr, const Eigen::Vector2d& vNext, const Eigen::Vector2d& vTest);

double distance_from_edges(const Eigen::VectorXd& pt, const Eigen::MatrixXd& vertices);

void interpolate_edges(
	const Eigen::MatrixXd& vertices,
	Eigen::MatrixXd& interpolatedVertices,
	std::vector<int>& newLocalFace,
	const double& edgeSize);

void hole_points(const double& radius, const Eigen::Vector2d& center, const int& ptNr, Eigen::MatrixXd& vertices);

bool ear_clipping(const Eigen::MatrixXd& vertices, const std::vector<int>& idxs, std::vector<std::vector<int>>& cells);

bool vertex_locally_convex(const Eigen::Vector2d& v1, const Eigen::Vector2d& v2, const Eigen::Vector2d& v3);

void back_to_3d(Eigen::MatrixXd& vertices3d, const Eigen::MatrixXd& vertices2d, const Eigen::Vector3d& center, const Eigen::Vector3d& u, const Eigen::Vector3d& v);
// implementation of circular double linked list

void catmull_rom_interpolation(const Eigen::MatrixXd& currentVerts, Eigen::MatrixXd& interpolatedVerts, double alpha=0.5);

double catmull_rom_get_t(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, double t0, double alpha = 0.5);

void chaikin_subdivision(const Eigen::MatrixXd& currentVerts, Eigen::MatrixXd& interpolatedVerts, int num = 5);

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

