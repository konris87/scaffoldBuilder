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
#include "SeedGenerator/DistanceCalculator.h"
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

    /*std::unique_ptr<BoxContainer> container = std::make_unique<BoxContainer>(
        0.0f, 5.0f, 0.0f, 5.0f, 0.0f, 5.0f
    );
    
    // load data
    std::vector<std::vector<float>> csvData = read_csv("../data/uniform_random_scaffolds.csv");

    std::array<int, 3> resolution = { 120, 120, 120 };

    int daMinsteps = 50;
    int daMaxsteps = 200;
    int daDirectionNr = 1000;
    float daTolerance = 1e-2f;
    int daFormulaIdx = 0;

    std::vector<SimOutput> results(csvData.size());

    // run as many simulations as the rows
    for (int row = 0; row < csvData.size(); row++) {

        float radius = csvData[row][1];
        float thickness = csvData[row][0];
        float openess = 0.5f;

        Poisson3D rng = Poisson3D(radius, radius, 30);

        // create the seeds
        rng.run(*container);

        std::vector<Vec3> seeds = rng.get_seeds();
        Bounds bds = container->compute_bounds();

        std::array<float, 6> bounds = {
            bds.xMin,
            bds.xMax,
            bds.yMin,
            bds.yMax,
            bds.zMin,
            bds.zMax
        };

        std::cout << "thickness: " << thickness << " radius: " << radius << " openess: " << openess << std::endl;

        // create the model using the thickness
        std::unique_ptr<GeneratorLewiner> scaffold = std::make_unique<GeneratorLewiner>(
            seeds, bounds, resolution, openess, thickness
        );


        // create the scaffold
        scaffold->compute_scalar_field(*container);

        scaffold->marching_cubes();

        scaffold->estimate_metrics(*container);

        scaffold->estimate_anisotropy(daDirectionNr, daMinsteps, daMaxsteps, daTolerance);

        scaffold->estimate_local_thickness(0.05, bounds);

        // estimate separation
        scaffold->estimate_local_thickness(0.05, bounds, 1);

        scaffold->estimate_trabecular_number();

        scaffold->estimate_connectivity_density();

        // write output to csv
        // create the output
        SimOutput res = {
            row + 1,
            thickness,
            radius,
            radius,
            openess,
            scaffold->porosity,
            scaffold->volume,
            scaffold->surfaceArea,
            scaffold->surfaceToVolume,
            scaffold->connectivityDensity,
            scaffold->localThickness,
            scaffold->localThicknessStd,
            scaffold->localSeparation,
            scaffold->localSeparationStd,
            scaffold->trabecularNr,
            scaffold->anisotropyDegree
        };

        results[row] = res;

        std::string fileName = "scaffold" + std::to_string(row) + ".stl";

        scaffold->export_stl("../data/uniform/" + fileName);

        std::cout << "scaffold " << row << " done!" << std::endl;
        std::cout << " ----------------------------------------- " << std::endl;
    }

    // write the results to a csv
    save_csv("../data/uniform/uniform_results.csv", results, SamplingType::uniform);
    */
    return 0;
}