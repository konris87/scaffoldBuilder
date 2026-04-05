#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <memory>

// custom headers
#include "guiApp.h"
#include "json.hpp"
#include "MonteSimulation/MonteSimulation.h"
#include "ScaffoldGenerator/GeneratorLewiner.h"
#include <SeedGenerator/Container.h>
#include <SeedGenerator/SeedGenerator.h>
//#include <Visualize/VisualizeSeeds.h>
//#include "ScaffoldGenerator/ScaffoldGenerator.h"
#include "SeedGenerator/RadiusCalculator.h"
#include "Utils/Utils.h"

// test global
//#include "ScaffoldGenerator/ScaffoldGenerator.h"
//#include "SeedGenerator/Poisson3D.h"
//#include "SeedGenerator/Random.h"

// ---------------------------------------------------------------------
// Globals
const unsigned int SCR_WIDTH = 1200;
const unsigned int SCR_HEIGHT = 800;

int main(){
    
    myGUI myGui(SCR_WIDTH, SCR_HEIGHT);
    myGui.run();

     //run first the uniform
    //run_monte_carlo_simulations(2000, 3, SamplingType::random);
    
    //std::unique_ptr<BoxContainer> container = std::make_unique<BoxContainer>(
    //    0.0f, 5.0f, 0.0f, 5.0f, 0.0f, 5.0f
    //);
    //
    //// load data
    ////std::vector<std::vector<float>> csvData = read_csv("../data/uniform_random_scaffolds.csv");
    ////std::vector<std::vector<float>> csvData = read_csv("../data/bonej_comparison.csv");
    //std::vector<std::vector<float>> csvData = read_csv("../data/varied.csv");

    //std::array<int, 3> resolution = { 120, 120, 120 };
    //
    //float voxelSize = 0.025;
    //float openess = 0.5f;

    //std::vector<SimOutput> results(csvData.size());

    //// run as many simulations as the rows
    //for (int row = 0; row < csvData.size(); row++) {

    //    float thickness = csvData[row][0];
    //    float radiusMin = csvData[row][1];
    //    float radiusMax = csvData[row][2];
    //    float stretchY = csvData[row][3];

    //    Poisson3D rng = Poisson3D(radiusMin, radiusMax, 30);

    //    RunConfig cfg;
    //    cfg.dist = container->get_distance_estimator();
    //    cfg.rad = std::make_shared<QuadraticFunction>(radiusMin);

    //    // create the seeds
    //    rng.run(*container, cfg);

    //    std::vector<Vec3> seeds = rng.get_seeds();
    //    Bounds bds = container->compute_bounds();

    //    std::array<float, 6> bounds = {
    //        bds.xMin,
    //        bds.xMax,
    //        bds.yMin,
    //        bds.yMax,
    //        bds.zMin,
    //        bds.zMax
    //    };
    //    // create the model using the thickness
    //    std::unique_ptr<GeneratorLewiner> scaffold = std::make_unique<GeneratorLewiner>(
    //        seeds, bounds, resolution, openess, thickness
    //    );

    //    scaffold->stretchY = stretchY;

    //    // create the scaffold
    //    scaffold->compute_scalar_field(*container);

    //    scaffold->marching_cubes();

    //    scaffold->estimate_metrics(*container);

    //    scaffold->estimate_anisotropy(2000, 10000, 3);

    //    scaffold->estimate_local_thickness(voxelSize, bounds);

    //    // estimate separation
    //    scaffold->estimate_local_thickness(voxelSize, bounds, 1);

    ////    scaffold->estimate_trabecular_number();

    //    scaffold->estimate_connectivity_density();

    //    // write output to csv
    //    // create the output
    //    SimOutput res = {
    //        row + 1,
    //        thickness,
    //        radiusMin,
    //        radiusMax,
    //        openess,
    //        scaffold->porosity,
    //        scaffold->volume,
    //        scaffold->surfaceArea,
    //        scaffold->surfaceToVolume,
    //        scaffold->connectivityDensity,
    //        scaffold->localThickness,
    //        scaffold->localThicknessStd,
    //        scaffold->localSeparation,
    //        scaffold->localSeparationStd,
    //        scaffold->trabecularNr,
    //        scaffold->anisotropyDegree
    //    };

    //    results[row] = res;

    //    std::string stlFileName = "scaffold" + std::to_string(row) + ".stl";
    //    std::string nrrdFileName = "scaffold" + std::to_string(row) + ".nrrd";

    //    //std::string fName = std::format("");

    //    std::stringstream stlStream;
    //    stlStream << "../data/varied/scaffold" << row << ".stl";
    //    scaffold->export_stl(stlStream.str());

    //    stlStream.str(std::string());
    //    stlStream << "../data/varied/scaffold" << row << ".nrrd";
    //    scaffold->export_nrrd(stlStream.str(), voxelSize, bounds);

    //    std::cout << "scaffold " << row << " done!" << std::endl;
    //    std::cout << " ----------------------------------------- " << std::endl;
    //}

    //// write the results to a csv
    //std::stringstream folderName;
    //folderName << "../data/varied/op" << openess << "_vs" << voxelSize << "_varied_quadratic_results.csv";

    //std::cout << folderName.str() << std::endl;

    //save_csv(folderName.str(), results, SamplingType::uniform);
    

    return 0;
}