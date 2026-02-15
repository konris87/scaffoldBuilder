#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>

// custom headers
#include "guiApp.h"
#include "json.hpp"
#include "MonteSimulation/MonteSimulation.h"

#include <SeedGenerator/Container.h>
#include <SeedGenerator/SeedGenerator.h>
#include <Visualize/VisualizeSeeds.h>
//#include "ScaffoldGenerator/ScaffoldGenerator.h"
#include "SeedGenerator/DistanceCalculator.h"

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

    // run first the uniform
    //run_monte_carlo_simulations(2000, 3, SamplingType::random);

    //BoxContainer boxContainer;
    //IContainer* activeContainer = &boxContainer;
    //// Re-initialize adapter with the fresh container dimensions (0-10)
    //ContainerAdapter adapter = { *activeContainer, 10, 10, 10 };
    //std::vector<std::array<double, 3>> seeds;
    //int nrSeeds = 2;
    //Random rnd(nrSeeds);
    //rnd.run(adapter, seeds);
    //// loop inside cells
    //run_test_voro(seeds, 1.0);

    return 0;
}