#ifndef WALL_H
#define WALL_H

#include <voro++.hh>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkKdTreePointLocator.h>
#include "Eigen/Dense"

using namespace voro;

class MeshWall : public voro::wall {
private:
    std::vector<std::array<double, 3>> barycenters;
    std::vector<std::array<double, 3>> normals;
    std::vector<std::array<double, 3>> trianglePt;
    std::vector<std::array<double, 3>> vertexNormals;
    vtkSmartPointer<vtkKdTreePointLocator> kdTree;
    
    bool isInside(Eigen::Vector3d& point);

    vtkSmartPointer<vtkPolyData> mesh;
    
    int w_id{ 0 };

    int neighbors{ 1 };

public:
    MeshWall() {};
    ~MeshWall() {};
    MeshWall(
        const vtkSmartPointer<vtkPolyData>& meshData,
        vtkSmartPointer<vtkKdTreePointLocator> kdTree,
        const int neighbors,
        int iw_id = 0);

    MeshWall(
        const vtkSmartPointer<vtkPolyData>& meshData,
        const std::vector<std::array<double, 3>>& vertexNormals,
        vtkSmartPointer<vtkKdTreePointLocator> kdTree,
        const int neighbors,
        int iw_id = 0);

    MeshWall(
        const vtkSmartPointer<vtkPolyData>& meshData,
        const std::vector<std::array<double, 3>>& barycenters,
        const std::vector<std::array<double, 3>>& faceNormals,
        const std::vector<std::array<double, 3>>& trianglePt,
        vtkSmartPointer<vtkKdTreePointLocator> kdTree,
        const int neighbors,
        int iw_id = 0);

    bool point_inside(double x, double y, double z);

    template < class vc_class>
    inline bool cut_cell_base(vc_class& c, double x, double y, double z) {

        double queryPoint[3] = { x, y, z };
        vtkSmartPointer<vtkIdList> closestIds = vtkSmartPointer<vtkIdList>::New();

        assert(!std::isnan(queryPoint[0]) && !std::isnan(queryPoint[1]) && !std::isnan(queryPoint[2]));

        // Find the nearest N points to the cell center
        kdTree->FindClosestNPoints(neighbors, queryPoint, closestIds);

        std::vector<std::array<double, 3>> closestPoints;

        for (vtkIdType i = 0; i < closestIds->GetNumberOfIds(); ++i) {
            
            double point[3];
            kdTree->GetDataSet()->GetPoint(closestIds->GetId(i), point);
            
            int closestId = closestIds->GetId(i);

            // get the neighbor point in local coordinates
            double px = point[0] - x;
            double py = point[1] - y;
            double pz = point[2] - z;
            
            // get its normal
            std::array<double, 3> norm = vertexNormals[closestId];
            double xc = norm[0];
            double yc = norm[1];
            double zc = norm[2];
 
            // find distance from plane and estimate the rsq
            double d = - (xc * px + yc * py + zc * pz);

            //std::cout << "distance: " << d << std::endl;
            double offset;
            offset = 0.5 * sqrt(pow(xc * d, 2) + pow(yc * d, 2) + pow(zc * d, 2));

            c.nplane(xc, yc, zc, offset, w_id);
            ++w_id;
        }
        return true;
    }

    // These virtual functions are called during the cell
    // computation in the container class. They call instances of
    // the template given above.
    bool cut_cell(voronoicell & c, double x, double y, double z) {
        return cut_cell_base(c, x, y, z);
    }
    bool cut_cell(voronoicell_neighbor & c, double x, double y, double z) {
        return cut_cell_base(c, x, y, z);
    }
};

#endif