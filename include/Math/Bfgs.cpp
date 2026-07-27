#include <iostream>
#include "Bfgs.h"

BFGS::BFGS(const fName& f, const gName& df) : func(f), grad(df) {};

double BFGS::line_search(){

	double aStar{0};

	xStep = x;
	Eigen::VectorXd xi;

	// function and gradient value at the start of line search
	fInit = func(xStep);
	gInit = grad(xStep).dot(p);

	// check that we have a descent direction
	if (gInit > 0) {
		std::logic_error("The moving direction increases the objective function");
	}

	double fi{ 0 }, gi{ 0 };
	double fStep{ 0 }, gStep{ 0 };

	// initial step interval a0 -> 0, a1->1 
	double a0{0};
	double ai{1};

	for (int i = 0; i < searchMaxIter; i++) {

		xi = xStep + ai * p;
		fi = func(xStep + ai * p);
		gi = grad(xStep + ai * p).dot(p);

		std::cout << "fi" << fi << std::endl;
		std::cout << "gi" << gi << std::endl;
		std::cout << "abs(gi)" << abs(gi) << std::endl;
		std::cout << "gInit" << -c2 * gInit << std::endl;

		// Condition 1, Armijo rule
		if (fi > fInit + c1 * ai * gInit || (fi >= fStep && i>0)) {
			std::cout << "condition 1" << std::endl;
			aStar = zoom(a0, ai);
			std::cout << aStar << std::endl;
			return aStar;
		}

		// Case 2, evaluate gi

		if (abs(gi) <= -c2 * gInit) {
			std::cout << "condition 2" << std::endl;
			aStar = ai;
			std::cout << aStar << std::endl;
			return aStar;
		}
		
		// case 3
		if (gi >= 0) {
			std::cout << "condition 3" << std::endl;
			aStar = zoom(ai, a0);
			return aStar;
		}	

		// if nothing of the above set ai to 2*ai and continue
		std::cout << "condition None" << std::endl;
		ai *= 2;
		xStep = xi;
		fStep = fi;
		gStep = gi;
		a0 = ai;
	}

	std::cout << "Max iterations reached for line search! " << std::endl;
	return ai;
};

double BFGS::zoom(double alo, double ahi){
	
	// final accepted step
	double aStar{0};

	for (int j = 0; j < searchMaxIter; j++){

		double flo = func(xStep + alo * p);
		double fhi = func(xStep + ahi * p);
		double glo = grad(xStep + alo * p).dot(p);
		double ghi = grad(xStep + ahi * p).dot(p);
	
		double aj = cubic_interpolation(alo, ahi, flo, fhi, glo, ghi);

		// test if nan use bisection
		const bool candid_nan = !(std::isfinite(aj));

		aj = candid_nan ? (alo + ahi) / 2 : aj;

		// evaluate fj and gj
		double fj = func(xStep + aj * p);
		double gj = grad(xStep + aj * p).dot(p);

		// condition 1
		if (fj > fInit + c1 * aj * gInit || fj >= flo) {
			
			// check numeric precision
			if (aj == ahi){
				std::runtime_error(
					"line search routine failed maybe due to insufficient numeric precision");
			}
			ahi = aj;
		}
			
		else {
			if (abs(gj) <= -c2 * gInit) {
				aStar = aj;
				return aStar;
			}
			if (gj * (ahi - alo) >= 0) {
				ahi = alo;
			}

			// check numeric precision
			if (aj == alo){
				std::runtime_error(
					"line search routine failed maybe due to insufficient numeric precision");
			}
			alo = aj;
		}
	}
	// if we have used all of our iterations return the lowest
	aStar = alo;
	return aStar;
};

double BFGS::cubic_interpolation(
	double a0, double ai, double f0, double fi, double g0, double gi){
	
	double alpha{0} ;
	double d1 = g0 + gi - 3 * ((f0 - fi) / (a0 - ai));
	double d2 = sqrt(pow(d1, 2) - g0 * gi);
	
	if (std::signbit(ai - a0)) {
		d2 *= -1;
	}

	alpha = ai - (ai - a0) * ((gi + d2 - d1) / (gi - g0 + 2 * d2));

	return alpha;
};

void BFGS::minimize(Eigen::VectorXd& xInit){

	x = xInit;

	// Hessian initialization
	Eigen::MatrixXd H = Eigen::MatrixXd::Identity(x.size(), x.size());
	const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(x.size(), x.size());

	Eigen::VectorXd g = grad(x);

	// search direction
	p = -H * g;

	// line search
	double alpha = 1.0;

	// find xk
	Eigen::VectorXd x0 = x + alpha * p;

	// find sk
	Eigen::VectorXd s0 = x0 - x;

	// find yk
	Eigen::VectorXd g0 = grad(x0);
	Eigen::VectorXd y0 = g0 - g;

	// initial Hessian
	H =  H * (y0.dot(s0) / y0.dot(y0));

	double ak{0};

	// start main loop
	for (int iter{0}; iter<maxIter; iter++){

		std::cout << "Iteration: " << iter + 1 << std::endl;

		p = - H * g;

		// line search to estimate alpha using Cubic Interpolation
		ak = line_search(); 

		std::cout << "step length: " << ak << std::endl;

		// after finding step update solution
		// find xk
		Eigen::VectorXd xk = x + ak * p;

		// find sk
		Eigen::VectorXd sk = xk - x;

		// find yk
		Eigen::VectorXd gk = grad(xk);
		Eigen::VectorXd yk = gk - g;

		// Update Hessian using BFGS
		double rho = 1.0 / yk.dot(sk);

		H = (I - rho * sk * yk.transpose()) * H * (I - rho * yk * sk.transpose()) + rho * sk * sk.transpose();

		x = xk;
		g = gk;

		if (g.norm() < epsilon) {
			std::cout << "Converged at iteration " << iter + 1 << std::endl;
			std::cout << "Function Value: " << func(xk) << std::endl;
			std::cout << "Grad Norm: " << g.norm() << std::endl;
			std::cout << "x = \n" << xk.transpose() << std::endl;
			
			break;
		}

		std::cout << "Function Value: " << func(xk) << std::endl;
		std::cout << "Grad Norm: " << g.norm() << std::endl;

		std::cout << "-----------------------------------" << std::endl;

	}

};