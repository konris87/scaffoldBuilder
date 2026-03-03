#ifndef KDTREE_H
#define KDTREE_H

#include <array>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <queue>
#include "Vec.h"

//template <typename pointType, size_t K>
class Kdtree {
public:
	Kdtree() { _root = nullptr; };
	~Kdtree() {};

	// constructor to build tree from points
	Kdtree(const std::vector<Vec3>& points) : _root(nullptr) {
		
		// change the size of the vector _nodes
		_nodes.reserve(points.size());

		// populate the vector _nodes
		for (size_t i{ 0 }; i < points.size(); i++) {
			_nodes.emplace_back(points[i], i);
		}

		// now build the tree
		_root = populateRecursively(_nodes.begin(), _nodes.end(), 0UL);

	};

	//@brief KNN Search
	//@param query: the query point for which we want to find the nearest neighbors
	//@param k: the number of nearest neighbors to find
	//@param distFunc: a function that takes two points and returns the distance between them, it should be defined as: double distFunc(const std::array<pointType, K>& p1, const std::array<pointType, K>& p2)
	//@returns a vector of pairs of index and distance, sorted by distance (closest first)
	template <typename DistFunction>
	std::vector<std::pair<size_t, double>> knn(
		const Vec3& query,
		size_t k,
		DistFunction distFunc) const
	{
		// Max-heap to store the k-nearest found so far
		std::priority_queue<Neighbor> pq;

		knnIterative(_root, query, k, pq, distFunc);

		// Convert priority queue to result vector
		std::vector<std::pair<size_t, double>> result;
		result.reserve(pq.size());

		// PQ pops usually from largest to smallest, but we usually want 
		// smallest to largest (d1, then d2). 
		while (!pq.empty()) {
			Neighbor n = pq.top();
			pq.pop();
			result.push_back({ n.node->index, n.distSq });
		}

		// Reverse so index 0 is the closest (d1), index 1 is second closest (d2)
		std::reverse(result.begin(), result.end());
		return result;
	}


	template <typename DistFunction>
	std::pair<Vec3, size_t>	nnp(const Vec3& query, DistFunction distFunc) const {

		// create an empty node that represents the best guess
		Node* best = nullptr;
		double bestDist = std::numeric_limits<double>::max();

		nearestRecurse(_root, query, 0, best, bestDist, distFunc);

		if (best) {
			return { best->point, best->index };
		}
		else{
			return { std::array<pointType, K> {}, size_t(-1) };
		};
	};

	template <typename DistFunction>
	std::vector<std::pair<Vec3, size_t>> radiusSearch(
		const Vec3& query, double radius, DistFunction distFunc) const {

		std::vector<std::pair<Vec3, size_t>> result;
		radiusRecurse(_root, query, 0, radius * radius, result, distFunc);
		return result;
	}

private:

	struct Node {

		// the corresponding point
		Vec3 point;

		// left child
		Node* left = nullptr;

		// right child
		Node* right = nullptr;

		// index if node
		size_t index{ 0 };;

		// Constructor
		Node(const Vec3& pt, size_t idx) : point(pt), index(idx), left(nullptr), right(nullptr) {};
	};

	//@brief structure to hold the neighbors during the search, we use a max-heap to keep track of the k nearest neighbors
	struct Neighbor {
		double distSq;
		Node* node;

		// Max-heap needs to order by largest distance first
		bool operator<(const Neighbor& other) const {
			return distSq < other.distSq;
		}
	};

	struct StackNode {
		Node* node = nullptr;
		int depth = 0;
	};

	using Iterator = typename std::vector<Node>::iterator;

	// root node
	Node* _root = nullptr;

	// a vector of Nodes
	std::vector<Node> _nodes;

	//@ we pass only the range of the nodes we want to process and the index
	Node* populateRecursively(
		Iterator begin,
		Iterator end,
		std::size_t index
	) {

		// first check if the range is empty
		if (begin >= end) {
			return nullptr;
		}

		// else keep going finding the median
		auto median = begin + std::distance(begin, end) / 2;

		// now partially sort so that:
		// all elements before median are less than or equal to the median for that dimension 
		// all elements before median are greater than or equal to the median for that dimension 
		// the nth element at position middle is the median based on dimension index
		// the comp function is defined with a lamda expression that captures index by reference
		std::nth_element(begin, median, end, 
			[&index](const Node &n1, const Node &n2) -> bool {
			return (n1.point[index] < n2.point[index]);
		});

		index = (index + 1) % 3;

		median->left = populateRecursively(begin, median, index);
		median->right = populateRecursively(median+1, end, index);

		return &(*median);

	};

	template <typename DistFunction>
	void nearestRecurse(
		Node* node, 
		const Vec3& query,
		int depth, 
		Node*& best, double& bestDist,
		DistFunction distFunc
		) const {
		
		if (!node) {
			return;
		}

		double dist = distFunc(node->point, query);

		if (dist < bestDist) {
			bestDist = dist;
			best = node;
		}

		// get next node to check
		int axis = depth % 3;

		double diff = query[axis] - node->point[axis];
        Node* nearChild = (diff < 0) ? node->left : node->right;
        Node* farChild  = (diff < 0) ? node->right : node->left;

		nearestRecurse(nearChild, query, depth + 1, best, bestDist, distFunc);

		double diff = query[axis] - node->point[axis];
		if (diff * diff < bestDist) {
			nearestRecurse(farChild, query, depth + 1, best, bestDist, distFunc);
		}
	};

	template <typename DistFunction>
	void radiusRecurse(
		Node* node,
		const Vec3& query,
		int depth,
		double radiusSq,
		std::vector<std::pair<Vec3, size_t>>& result,
		DistFunction distFunc) const {

		if (!node) return;

		double dist = distFunc(node->point, query);

		if (dist <= radiusSq) {
			result.push_back({ node->point, node->index });
		}

		int axis = depth % K;

		double diff = query[axis] - node->point[axis];

		Node* nearChild = (diff < 0) ? node->left : node->right;
		Node* farChild = (diff < 0) ? node->right : node->left;

		// Always check near child
		radiusRecurse(nearChild, query, depth + 1, radiusSq, result, distFunc);

		// Only check far child if the hypersphere crosses the splitting plane
		if (diff * diff <= radiusSq) {
			radiusRecurse(farChild, query, depth + 1, radiusSq, result, distFunc);
		}
	}

	template <typename DistFunction>
	void knnIterative(
		Node* root,
		const Vec3& query,
		size_t k,
		std::priority_queue<Neighbor>& pq,
		DistFunction distFunc) const
	{
		if (!root) return;

		// Use a manual stack to avoid recursion stack overflow
		// Max depth of 64 is enough for billions of points
		std::vector<StackNode> stack;
		stack.reserve(64);
		stack.push_back({ root, 0 });

		while (!stack.empty()) {
			StackNode current = stack.back();
			stack.pop_back();

			Node* node = current.node;
			if (!node) continue;

			int depth = current.depth;
			double distSq = distFunc(node->point, query);

			// Manage Priority Queue
			if (pq.size() < k) {
				pq.push({ distSq, node });
			}
			else if (distSq < pq.top().distSq) {
				pq.pop();
				pq.push({ distSq, node });
			}

			// Determine splitting axis
			int axis = depth % 3;
			double diff = query[axis] - node->point[axis];

			Node* nearChild = (diff < 0) ? node->left : node->right;
			Node* farChild = (diff < 0) ? node->right : node->left;

			// Logic for far child check
			bool checkFar = (pq.size() < k || (diff * diff) < pq.top().distSq);

			// Push children onto our manual stack
			// We push farChild first so that nearChild is processed first (LIFO)
			if (farChild && checkFar) {
				stack.push_back({ farChild, depth + 1 });
			}
			if (nearChild) {
				stack.push_back({ nearChild, depth + 1 });
			}
		}
	}

	// use squared distance to avoid sqrt
	double distanceSquared(const Vec3& point, const Vec3& query) {

		double sum{ 0.0 };

		for (size_t i = 0; i < 3; i++) {
			double diff = point[i] - query[i];
			sum += diff * diff;
		};

		return sum;
	};

};

#endif