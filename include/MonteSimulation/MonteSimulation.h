#ifndef MONTESIMULATION_H
#define MONTESIMULATION_H

#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <array>
#include <memory>
#include "Utils/Utils.h"

enum SamplingType {
	random, uniform, varied
};

// @brief structure to hold simulation inputs
struct SimInput {

	SamplingType type;
    int id;
    
	double p1;
	double p2;
    double thickness;      
	double pullbackRatio;  
};

struct SimOutput {
	int id;
    double porosity;
    double svr;
    double connectivity;   // 0.0 to 1.0
};

// @class a class to create 
class LhsSampler {

public:
	LhsSampler() {};
	~LhsSampler() {};
	
	LhsSampler(int nr, int dims);
	
	std::vector<std::vector<double>> generate_parameters();

private:
	int sampleNr {100};
	int dimNr {3};
	
};

// ---------------------------------------------------------------------
// @brief generate lhs samples for scaffold
std::vector<SimInput> generate_lhs_sampling(int nr = 1000, int dims = 3, SamplingType t=uniform);

// @brief this function runs monte carlo simulations for the desired input parameters
void run_monte_carlo_simulations(int nr = 1000, int dims = 3, SamplingType t=uniform);

void run_voro(std::vector<std::array<double, 3>>& seeds, std::unique_ptr<Graph>& graph,
	const SimInput& params, SimOutput& out);

void save_csv(const std::string fileName,
	const std::vector<SimInput>& inputs,
	const std::vector<SimOutput>& outputs,
	const SamplingType& t);

#endif MONTESIMULATION_H