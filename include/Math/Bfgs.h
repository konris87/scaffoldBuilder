#ifndef BFGS_H // include guard
#define BFGS_H

#include <Eigen/Dense>
#include <cmath>

typedef std::function<double(const Eigen::VectorXd)> fName;
typedef std::function<Eigen::VectorXd(const Eigen::VectorXd&)> gName;

class BFGS {

public:	
	// constructor
	BFGS(const fName& f, const gName& df);

	// minimize function
	void minimize(Eigen::VectorXd& xInit);
			
	// LBFGS parameters
	double epsilon = 1e-5;
	int maxIter = 100;
		
	// line search parameters
	int searchMaxIter = 20;
	double minStep = 1e-20;
	double maxStep = 1e20;
	double c1 = 1e-4;
	double c2 = 0.9;
	
	// print information to the screen
	bool verboseFlag{true};	

private:
	// function and its gradient
	fName func;
	gName grad;

	// current solution
	Eigen::VectorXd x;

	// step solution for line search
	Eigen::VectorXd xStep;
	double fInit{0}, gInit{0};

	// current direction
	Eigen::VectorXd p;

	// current func and gradient evaluation
	double fEval;
	double gEval;

	// line search
	double line_search();

	// zoom phase
	double zoom(double aLow, double aHigh);

	// cubic interpolation
	double cubic_interpolation(
		double a0, double ai, double f0, double fi, double g0, double gi);

	// function that prints results in the screen
	void verbose(){};
};


#endif