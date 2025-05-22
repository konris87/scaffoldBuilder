#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>

// custom headers
//#include "Shader.h"
//#include "Mesh.h"
//#include "Model.h"
//#include "guiApp.h"
//#include "json.hpp"
//#include "OpenGlSetup/Grid.h"

// test global
#include "ScaffoldGenerator/ScaffoldGenerator.h"
#include "SeedGenerator/Poisson3D.h"

// ---------------------------------------------------------------------
// Globals
const unsigned int SCR_WIDTH = 1200;
const unsigned int SCR_HEIGHT = 800;

int main(){
    
    //myGUI myGui(SCR_WIDTH, SCR_HEIGHT);

    //myGui.run();

    // Map: Vec3 -> global index
    //std::unordered_map<std::array<double, 3>, int, GlobalVertexKey, GlobalVertexInside> vertex_to_index;

    // Vector: global index -> Vec3 (reverse lookup)
    //std::vector<std::array<double, 3>> global_vertices;

    //// Sample vertices (some duplicates)
    //std::vector<std::array<double, 3>> input_vertices = {
    //    {0.0000000e+00,  0.0000000e+00, 0.0} ,
    //    {1.8662844e+00,  4.3032864e-18, 0.0},
    //    {1.9789295e+00,  2.0352561e+00, 0.0},
    //    {- 2.7037561e-01,  2.1597481e+00, 0.0},
    //    {1.9226069e+00,  1.0176281e+00, 0.0 },
    //    {8.5427696e-01,  2.0975020e+00, 0.0},
    //    {-1.3518780e-01,  1.0798740e+00, 0.0},
    //    {1.3280561e+00,  1.0578358e+00, 0.0},
    //    {8.7722903e-01,  6.0700876e-01, 0.0},
    //    {4.2640197e-01,  1.0578358e+00, 0.0},
    //    {8.7722903e-01,  1.5086629e+00, 0.0}
    //};

    //for (const auto& v : input_vertices) {
    //    // Check if already in map
    //    auto it = vertex_to_index.find(v);
    //    if (it == vertex_to_index.end()) {
    //        int index = global_vertices.size();
    //        global_vertices.push_back(v);
    //        vertex_to_index[v] = index;
    //        std::cout << "Added vertex " << index << ": ("
    //            << v[0] << ", " << v[1] << ", " << v[2] << ")\n";
    //    }
    //    else {
    //        std::cout << "Duplicate vertex found, index = " << it->second << "\n";
    //    }
    //}

    //std::cout << "\n=== Global Vertex List ===\n";
    //for (std::size_t i = 0; i < global_vertices.size(); ++i) {
    //    const auto& v = global_vertices[i];
    //    std::cout << i << ": (" << v[0] << ", " << v[1] << ", " << v[2] << ")\n";
    //}

    // create some seeds
    double rMin = 2;
    
    std::array<double, 3> root{ 5.0, 5.0, 5.0 };

    std::array<float, 6> bounds{ 0.0, 10.0, 0.0, 10.0, 0.0, 10.0 };

    Poisson3D sg = Poisson3D(rMin, rMin, root, bounds);

    std::vector<std::array<double, 3>> test;
    sg.generate_seeds();
    sg.get_seeds(test);
    std::cout << test.size() << std::endl;


    std::vector<std::array<double, 3>> seeds = {
        
        {0.0, 0.1, 2.0},
        //{4.0, 1.1, 2.0},
        //{3.0, 2.0 , 1.0},
        //{1.23, 4.5, 3.1},
        //{3.0, 3.0, 3.8},
        //{2.1, 1.96, 7.8},
        //{4.1, 3.8, 8.1},
        //{ 8.1, 3.8, 7.1},
        //{9.1, 1.8, 2.1},
        //{6.1, 4.8, 5.1}
    
    };
    
    // get the seeds
    //std::vector<std::array<double, 3>> seeds;
    //sg.get_seeds(seeds);

    std::cout << seeds.size() << std::endl;
        
    // create voronoi
    std::array<int, 3> blockDim{ 100, 100, 100 };

    ScaffoldGeneratorFaceBox sgfb = ScaffoldGeneratorFaceBox(test, bounds, blockDim, 4, 4, 1.0);

    sgfb.generate_voro();

    return 0;
}