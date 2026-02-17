#ifndef RADIUSCALCULATOR_H
#define RADIUSCALCULATOR_H

#include <random>

class RadiusFunction {
public:
	virtual double estimate_radius(double distance, double rMin, double rMax) const = 0;
	virtual ~RadiusFunction() = default;
};

class LinearFunction : public RadiusFunction {

public:
	LinearFunction() {};
	~LinearFunction() {};
	LinearFunction(double distMax) : distMax(distMax) {};

	double estimate_radius(double distance, double rMin, double rMax) const override {
		if (distance < 0) {
			return rMin;
		}
		else if (distance > distMax) {
			return rMax;
		}
		else {
			return ((rMax - rMin) / distMax) * distance + rMin;
		}
	}
	//double estimate_radius(double distance, double rMin, double rMax) const override {
	//	return	rMax - (rMax - rMin) * std::min(distance / 5.0, 1.0);
	//};

private:
	double distMax{ 10.0 };
};

class QuadraticFunction : public RadiusFunction {
public:
	explicit QuadraticFunction(double distMax) : distMax(distMax) {}
	double estimate_radius(double distance, double rMin, double rMax) const override {
		const double t = std::clamp(distance / std::max(1e-9, distMax), 0.0, 1.0);
		return rMin + (rMax - rMin) * (1.0 - t * t);
	}
private:
	double distMax{ 10.0 };
};

class ConstantRadiusFunction : public RadiusFunction {
public:
	double estimate_radius(double, double rMin, double) const override { return rMin; }
};

class RandomRadiusFunction : public RadiusFunction {
public:
	double estimate_radius(double, double rMin, double rMax) const override {

		std::uniform_real_distribution<double> dist(rMin, rMax);
		std::mt19937 rng;
		return dist(rng);
	}
};


#endif 
