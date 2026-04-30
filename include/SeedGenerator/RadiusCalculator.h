#ifndef RADIUSCALCULATOR_H
#define RADIUSCALCULATOR_H

#include <random>

class RadiusFunction {
public:
	virtual double estimate_radius(double distance, double tStart, double tEnd) const = 0;
	virtual ~RadiusFunction() = default;
};

class LinearFunction : public RadiusFunction {

public:
	LinearFunction() {};
	~LinearFunction() {};
	LinearFunction(double distMax) : distMax(distMax) {};

	double estimate_radius(double distance, double tStart, double tEnd) const override {
		if (distance < 0) {
			return tStart;
		}
		else if (distance > distMax) {
			return tEnd;
		}
		else {
			double t = distance / distMax;
			return tStart + t * (tEnd - tStart);
		}
	}

private:
	double distMax{ 10.0 };
};

class QuadraticFunction : public RadiusFunction {
public:
	QuadraticFunction() {};
	QuadraticFunction(double distMax) : distMax(distMax) {}
	double estimate_radius(double distance, double tStart, double tEnd) const override {
		const double t = std::clamp(distance / std::max(1e-9, distMax), 0.0, 1.0);
		return tStart + (tEnd - tStart) * (t * t);
	}
private:
	double distMax{ 10.0 };
};

class ConstantRadiusFunction : public RadiusFunction {
public:
	double estimate_radius(double, double t, double) const override { return t; }
};

class RandomRadiusFunction : public RadiusFunction {
public:

	// add a constructor and create a default 
	RandomRadiusFunction() : rng(std::random_device{}()) {};

	double estimate_radius(double, double rMin, double rMax) const override {
		std::uniform_real_distribution<double> dist(rMin, rMax);
		return dist(rng);
	}
private:
	mutable std::mt19937 rng;
};


#endif 
