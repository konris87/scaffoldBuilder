#ifndef TRANSPORT_H // include guard
#define TRANSPORT_H

#include <Eigen/Dense>
#include <vector>

double tet_integral(const std::vector<double>& tetVerts);
double transport(const std::vector<double>& cellVertices, const std::vector<double>& cellCenter);
double tet_volume(const std::vector<double>& tetVerts);
double euclidean_distance(std::array<double, 3>& pts1, std::array<double, 3>& pts2);
std::array<double, 3> compute_centroid(const std::vector<double>& vertices);
Eigen::VectorXd convert_radii(const Eigen::VectorXd& weights);

#endif