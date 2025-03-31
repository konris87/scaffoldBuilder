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

// ---------------------------------------------------------------------
// Globals
const unsigned int SCR_WIDTH = 1200;
const unsigned int SCR_HEIGHT = 800;

int main(){
    
    myGUI myGui(SCR_WIDTH, SCR_HEIGHT);

    myGui.run();

    return 0;
}