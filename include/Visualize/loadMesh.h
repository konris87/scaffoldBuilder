#pragma once
#ifndef LOADMESH_H // include guard
#define LOADMESH_H

#include <vector>
#include <string>


void loadStlMesh(
	const std::string modelFile,
	std::vector<float>& vertices,
	std::vector<unsigned int>& indices,
	std::vector<float>& normals
);

#endif