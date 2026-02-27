#ifndef AABBTREE_H
#define AABBTREE_H

#include <vector>
#include <array>

struct AabbNode {
	std::array<float, 6> bounds = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
	int leftChild = -1;
	int rightChild = -1;
	std::vector<unsigned int> triangleIds;
};

class AabbTree {
public:

	AabbTree() {};
	~AabbTree() {};

	void populate(std::vector<int>& triIds);

private:

	std::vector<AabbNode*> tree;

};


#endif