#ifndef LINESEARCH_H // include guard
#define LINESEARCH_H

#include <iostream>
#include "Eigen/Dense"
#include <cmath>

template<typename obj>

class LineSearch {
public:

    LineSearch(
        obj& objF, 
        const double& param1=1e-4,
        const double& param2=0.9,
        const double& aInit=1.0)
         : func(objF), c1(param1), c2(param2), ai(aInit) {}

    // search step
    double search_step(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& p,
        double iNit){

            // values at iteration zero
            Eigen::VectorXd df0_vec = Eigen::VectorXd::Zero(x.size());
            const double f0 = func(x, df0_vec);
            const double df0 = df0_vec.dot(p);

            // std::cout << "f0 " << f0 << " df0 " << df0 << std::endl;

            // values at xi
            Eigen::VectorXd dfi_vec = Eigen::VectorXd::Zero(x.size());
            double fi;
            double dfi;

            // values at x i - 1
            double fPrev{f0}, dfPrev{df0}, aPrev{a0};

            // values for zoom phase
            double alo{0.0}, ahi{ai}, flo{0.0}, fhi{0.0}, dflo{0.0}, dfhi{0.0};
            Eigen::VectorXd dflo_vec = Eigen::VectorXd::Zero(x.size());
            Eigen::VectorXd dfhi_vec = Eigen::VectorXd::Zero(x.size());
            
            // std::cout << p << std::endl;
            ai = iNit;

            // bracketing face to find an accepted interval for steps
            for(int iter{0} ; iter<maxIter; iter++){
                // std::cout << "--------------------------------" << std::endl;
                // std::cout << "Start Bracketing" << std::endl;
                Eigen::VectorXd xi = x + ai * p;     
                fi = func(xi, dfi_vec);
                
                // check condition 1
                if ((fi > f0 + c1 * ai * df0) || (fi >= fPrev && iter>1)) {

                    alo = aPrev;
                    ahi = ai;
                    // break to go to the loop below
                    break;
                }
                dfi = dfi_vec.dot(p);
                if (std::abs(dfi) <= -c2 * df0){
                    return ai;
                }
                if (dfi >= 0){
                    ahi = aPrev;
                    alo = ai;
                    // break to start zoom phase
                    break;
                }
                // check if significant progress is made

                //if (std::abs(ai - aPrev) < searchTol){
                //    std::cerr << "Not significant progress in line search" << std::endl;
                //    return ai;
                //}
                // std::cout << "doubling" << std::endl;
                // if nothing of the above either multiply ai by a factor or use cubic interpolation
                aPrev = ai;
                ai *= 2;
                fPrev = fi;
                dfPrev = dfi;
            }

            // zoom phase
            for (int j{0}; j < maxIter; j++){
                
                // std::cout << "Zoom Phase" << std::endl;
                
                // double alot = alo;
                // double ahit = ahi;
                // alo = alot + 0.1 * (ahit - alot);
                // ahi = ahit - 0.5 * (ahit - alot);

                // std::cout << "alo " << alo << " ahi " << ahi << std::endl;

                flo = func(x + alo * p, dflo_vec);
                dflo = dflo_vec.dot(p);

                fhi = func(x + ahi * p, dfhi_vec);
                dfhi = dfhi_vec.dot(p);
                
                // first perform cubic intpd
                double aj = cubic_interpolation(
                    alo, ahi, flo, fhi, dflo, dfhi);

                bool candidNan =! (std::isfinite(aj));
                bool nearLimit = std::min(std::abs(alo - aj), std::abs(ahi - aj)) < 0.01 * std::abs(alo - ahi);
                bool outOfInterval = (aj < std::min(alo, ahi) || aj > std::max(alo, ahi));

                if (candidNan || nearLimit || outOfInterval){
                    aj = quad_interpolation(alo, ahi, flo, fhi, dflo);
                }

                // candidNan =! (std::isfinite(aj));
                // nearLimit = std::min(std::abs(alo - aj), std::abs(ahi - aj)) < 0.1 * std::abs(alo - ahi);

                // if (candidNan || nearLimit) {aj = (alo + ahi) / 2;}
                        // std::cout << "aj " << aj << std::endl;
                                
                Eigen::VectorXd xj = x + aj * p;     
                Eigen::VectorXd dfj_vec = Eigen::VectorXd::Zero(x.size());
                double fj = func(xj, dfj_vec);

                if(fj > f0 + c1 * aj * df0 || fj >= fPrev){
                    // std::cout << "zoom cond 1" << std::endl;
                    ahi = aj;
                }
                else{
                    double dfj = dfj_vec.dot(p);
                    if(std::abs(dfj) <= - c2 * df0){
                        // std::cout << "zoom cond 2" << std::endl;
                        return aj; 
                    }
                    if(dfj*(ahi - alo) >= 0){
                        // std::cout << "zoom cond 3" << std::endl;
                        ahi = alo;
                    }
                    alo = aj;
                }
                //if (std::abs(aj - aPrev) < searchTol){
                //std::cerr << "Not significant progress in line search" << std::endl;
                //return aj;
                //}
            }
            return ai;
        }

    obj& func;
    double searchTol = 1e-7;
    double c1{0};
    double c2{0};

private:

    double a0{0};
    double ai{1};

    int maxIter{20};

    // initial values
	double f0, df0;
	Eigen::VectorXd df0_vec;

    double quad_interpolation(
        double alo, double ahi,
        double flo, double fhi,
        double glo)
    {
        double alpha{0};

        double t1 = glo * (ahi - alo) * (ahi - alo);
        double t2 = 2 * (fhi - flo - glo * (ahi - alo));

        return alpha - (t1 / t2);
    }

    double cubic_interpolation(
        double alo, double ahi,
        double flo, double fhi,
        double glo, double ghi){

            double alpha{0};

            double d1 = glo + ghi - 3 * ((flo - fhi) / (alo - ahi));
		    double d2 = sqrt(pow(d1, 2) - glo * ghi);
		
            if (std::signbit(ahi - alo)) {
                d2 *= -1;
            }

            alpha = ahi - (ahi - alo) * ((ghi + d2 - d1) / (ghi - glo + 2 * d2));

            return alpha;
    }
};

#endif