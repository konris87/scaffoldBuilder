#ifndef BFGS_H // include guard
#define BFGS_H

#include <cmath>
#include <Eigen/Dense>

template <typename Objective>

class BFGS {

public:	
	// constructor
	BFGS(Objective& foo) : func(foo) {}

	// minimize function
	int minimize(Eigen::VectorXd& xInit){
		
		x = xInit;

		// Hessian initialization
		Eigen::MatrixXd H = Eigen::MatrixXd::Identity(x.size(), x.size());
		const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(x.size(), x.size());

		f = func(x, df);
		
		// search direction
		p = -H * df;

		// line search
		double alpha = 1.0;

		// find xk
		Eigen::VectorXd x0 = x + alpha * p;

		// find sk
		Eigen::VectorXd s0 = x0 - x;

		// find yk
		Eigen::VectorXd df0 = Eigen::VectorXd::Zero(df.size());
		double f0 = func(x0, df0);
		/*std::cout << "g" << g(0) << std::endl;
		std::cout << "g0" << g0(0) << std::endl;*/
		Eigen::VectorXd y0 = df0 - df;

		// initial Hessian
		H =  H * (y0.dot(s0) / y0.dot(y0));

		std::cout << "-----------------------------------" << std::endl;

		// start main loop
		for (int iter{0}; iter<maxIter; iter++){

			std::cout << "Iteration: " << iter + 1 << std::endl;

			p = - H * df;

			// line search to estimate alpha using Cubic Interpolation
			double ak = line_search(); 
			if (ak == NULL) {
				break;
			}

			std::cout << "step length: " << ak << std::endl;

			// after finding step update solution
			// find xk
			Eigen::VectorXd xk = x + ak * p;

			// find sk
			Eigen::VectorXd sk = xk - x;

			// find yk
			double fk = func(xk, dfk);
			Eigen::VectorXd yk = dfk - df;

			// Update Hessian using BFGS
			double rho = 1.0 / (yk.transpose() * sk);

			Eigen::MatrixXd Hk = (I - rho * sk * yk.transpose()) * H * (I - rho * yk * sk.transpose()) + rho * sk * sk.transpose();
			
			//Eigen::MatrixXd t1 = yk * yk.transpose() / (yk.transpose() * sk);
			//Eigen::MatrixXd t2 = H * sk * sk.transpose() * H.transpose() / (sk.transpose() * H * sk);
			//Eigen::MatrixXd Hk = H + t1 - t2;
			
			df = dfk;
			x = xk;
			f = fk;
			H = Hk;

			if (dfk.lpNorm<Eigen::Infinity>() < epsilon) {
				std::cout << "Converged at iteration " << iter + 1 << std::endl;
				std::cout << "Function Value: " << fk << std::endl;
				std::cout << "Infinity Norm: " << dfk.lpNorm<Eigen::Infinity>() << std::endl;
				std::cout << "Grad Norm: " << dfk.norm() << std::endl;
				/*std::cout << "x = \n" << xk.transpose() << std::endl;
				std::cout << func.VolError() << std::endl;
				*/
				return 1;
			}
			if (verboseFlag) {
				std::cout << "Function Value: " << fk << std::endl;
				std::cout << "Infinity Norm: " << dfk.lpNorm<Eigen::Infinity>() << std::endl;
				std::cout << "Grad Norm: " << dfk.norm() << std::endl;
				std::cout << "-----------------------------------" << std::endl;
			}
		}
	};
			
	// LBFGS parameters
	double epsilon = 1e-5;
	int maxIter = 100;
		
	// line search parameters
	int searchMaxIter = 20;
	double c1 = 1e-4;
	double c2 = 0.9;
	
	// print information to the screen
	bool verboseFlag{true};	

private:
	// obj function object
	Objective& func;

	// current solution
	Eigen::VectorXd x;

	// gradient vec at current and next position
	Eigen::VectorXd df = Eigen::VectorXd::Zero(x.size());
	Eigen::VectorXd dfk = Eigen::VectorXd::Zero(x.size());

	// function value
	double f{0};

	// current direction
	Eigen::VectorXd p;

	// line search
	double line_search(){

		//std::cout << "Starting Line search" << std::endl;

		double aStar{ 0 };
		double aMax{1000};

		Eigen::VectorXd xi;
		Eigen::VectorXd xPrev = x;
		Eigen::VectorXd dfi;

		// function and gradient value at the start of line search
		const double fInit = f;
		const double gInit = df.dot(p);
		
		// check that we have a descent direction
		if (gInit > 0) {
			throw std::logic_error("The moving direction increases the objective function");
		}

		double fi{ 0 }, gi{ 0 };
		double fStep{ fInit }, gStep{ gInit };

		double fPrev{ 0 }, gPrev{ 0 };

		// initial step interval a0 -> 0, a1 -> 1 
		double a0{0};
		double ai{1};
		double aPrev{ 0 };

		for (int i = 0; i < searchMaxIter; i++) {

			xi = xPrev + ai * p;
						
			fi = func(xi, dfi);

			// Condition 1, Armijo rule
			if (fi > fInit + c1 * ai * gInit || (fi >= fPrev && i>0)) {
				std::cout << "Condition 1" << std::endl;
				aStar = zoom(a0, ai, xPrev, fInit, gInit);
				return aStar;
			}
			// evaluate gradient at ai
			gi = dfi.dot(p);
			//if (gi < c2 * gInit) {
			if (std::abs(gi) <= -c2 * gInit) {
				std::cout << "Condition 2" << std::endl;
				aStar = ai;
				return aStar;
			}
			
			// case 3
			if (gi >= 0) {
				std::cout << "Condition 3" << std::endl;
				aStar = zoom(ai, a0, xPrev, fInit, gInit);
				return aStar;
			}	
			
			if (2 * ai - aPrev >= aMax) {
				ai = aMax;
			}
			else {
				xPrev = xi;
				fPrev = fi;
				gPrev = gi;
				ai *= 2;
				aPrev = ai;
			}
			// if nothing of the above set ai to 2*ai and continue

			//Eigen::VectorXd dfMax;
			//double fMax = func(x, dfMax);
			//double gMax = dfMax.dot(p);
			//ai = cubic_interpolation(ai, aMax, fi, fMax, gi, gMax);
			//a0 = ai;

		}

		std::cout << "Max iterations reached for line search! " << std::endl;
		return NULL;
	};

	// zoom phase
	double zoom(
		double& alo, double& ahi,
		const Eigen::VectorXd& xStep,
		const double& fInit, const double& gInit)
		{
		Eigen::VectorXd dfj;

		// final accepted step
		double aStar{0};

		for (int j = 0; j < searchMaxIter; j++){

			double flo = func(xStep + alo * p , dfj);
			double glo = dfj.dot(p);
			double fhi = func(xStep + ahi * p, dfj);
			double ghi = dfj.dot(p);
			//double fhi = func(xStep + ahi * p);
			//double glo = func.gradient(xStep + alo * p).dot(p);
			//double ghi = func.gradient(xStep + ahi * p).dot(p);
		
			double aj = cubic_interpolation(alo, ahi, flo, fhi, glo, ghi);

			// test if nan 
			bool candid_nan = !(std::isfinite(aj));

			// test if too close in the limits
			bool nearEnd = std::min(std::abs(alo - aj), std::abs(ahi - aj)) < 0.1 * std::abs(alo - ahi);

			//if (candid_nan || nearEnd) {
			//	aj = quad_interpolation(alo, ahi, flo, fhi, glo);
			//}

			//// test if too close in the limits
			//candid_nan = !(std::isfinite(aj));
			//nearEnd = std::min(std::abs(alo - aj), std::abs(ahi - aj)) < 0.1 * std::abs(alo - ahi);
			if (candid_nan || nearEnd) {
				aj = (alo + ahi) / 2;
			}
			// evaluate fj and gj
			double fj = func(xStep + aj * p, dfj);
			double gj = dfj.dot(p);
			//double gj = func.gradient(xStep + aj * p).dot(p);

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
				if (std::abs(gj) <= -c2 * gInit) {
				//if (gj < c2 * gInit) {
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

	// cubic interpolation
	double cubic_interpolation(
		double& a0, double& ai,
		double& f0, double& fi,
		double& g0, double& gi){
	
		double alpha{0} ;
		double d1 = g0 + gi - 3 * ((f0 - fi) / (a0 - ai));
		double d2 = sqrt(pow(d1, 2) - g0 * gi);
		
		if (std::signbit(ai - a0)) {
			d2 *= -1;
		}

		alpha = ai - (ai - a0) * ((gi + d2 - d1) / (gi - g0 + 2 * d2));

		return alpha;
	};

	double quad_interpolation(
		double& a0, double& ai,
		double& f0, double& fi,
		double& g0
		) {

		double d1 = g0 * std::pow(a0, 2);
		double d2 = 2 * (fi - f0 - g0 * a0);

		return -d1 / d2;
	};

	// function that prints results in the screen
	void verbose(){};
};


#endif