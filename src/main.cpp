#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>

// custom headers
#include "Shader.h"
#include "Mesh.h"
#include "Model.h"
#include "guiApp.h"
#include "json.hpp"
#include "OpenGlSetup/Grid.h"

// test global
#include "ScaffoldGenerator/ScaffoldGenerator.h"
#include "SeedGenerator/Poisson3D.h"
#include "SeedGenerator/Random.h"

// ---------------------------------------------------------------------
// Globals
const unsigned int SCR_WIDTH = 1200;
const unsigned int SCR_HEIGHT = 800;

int main(){
    
    myGUI myGui(SCR_WIDTH, SCR_HEIGHT);

    myGui.run();

    //int seedNr{ 2 };
    //RandomGenerator rg({0.0, 10.0, 0.0, 10.0, 0.0, 10.0}, seedNr);
    //rg.generate_seeds();
    //std::vector<std::array<double, 3>> seeds;
    //rg.get_seeds(seeds);
    //ScaffoldGeneratorBox sgb(
    //    seeds, { 0.0, 10.0, 0.0, 10.0, 0.0, 10.0 }, { 10, 10, 10 }, 20, 0.5
    //);
    //sgb.generate_voro(0);

    return 0;
}