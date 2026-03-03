#ifndef MATRIX_H
#define MATRIX_H

#include "Vec.h"
#include <vector>
#include <iostream>
#include <cmath>

class Matrix3 {
	
public:
	Matrix3();
	Matrix3(float val);
	~Matrix3() = default;
	
	// overload matrix multiplication
	Matrix3 operator*(Matrix3 const& mat2) const;

	Matrix3& operator+=(const Matrix3& mat2);

	Matrix3& operator-=(const Matrix3& mat2);
		
	// overload matrix addition / subtraction
	Matrix3 operator-(Matrix3 const& mat2) const;

	Matrix3 operator+(Matrix3 const& mat2) const;

	// overload matrix * scalar, scalar * matrix 
	Matrix3 operator*(float scalar) const;

	friend Matrix3 operator*(float scalar, const Matrix3& m);

	// determinant
	float det();

	// get row as vec3
	Vec3 row(int idx) const;
	
	// get col as vec3
	Vec3 col(int idx) const;
		
	// transpose
	Matrix3 transpose();
	
	// overload assign col
	
	// overload assign row
	
	// overload = operator (assign)	
	
	// overload []
	float& operator()(int, int);

	const float& operator()(int i, int j) const;

private:
	float mat[3][3];
};

// ----------------------------------------------------------------------
class Identity : public Matrix3 {
	
public:
	Identity() { 
		mat[0][0] = 1.0f;
		mat[1][1] = 1.0f;
		mat[2][2] = 1.0f;
	};
	~Identity() {};
private:
	unsigned int nRows = 3;
	unsigned int nCols = 3;
	std::vector<std::vector<float>> mat;
};

#endif