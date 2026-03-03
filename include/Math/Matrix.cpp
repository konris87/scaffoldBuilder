#include "Matrix.h"

Matrix3::Matrix3(){
	for (size_t nr{ 0 }; nr < 3; nr++) {
		for (size_t nc{ 0 }; nc < 3; nc++) {
			mat[nr][nc] = 0.0f;
		}
	}
};

Matrix3::Matrix3(float val) {
	for (size_t nr{ 0 }; nr < 3; nr++) {
		for (size_t nc{ 0 }; nc < 3; nc++) {
			mat[nr][nc] = val;
		}
	}
};

Matrix3 Matrix3::operator+(Matrix3 const& mat2) const {
	Matrix3 sum;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			sum(i, j) = mat[i][j] + mat2(i, j);
		}
	}
	return sum;
};

Matrix3 Matrix3::operator-(Matrix3 const& mat2) const {
	Matrix3 diff;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			diff(i, j) = mat[i][j] - mat2(i, j);
		}
	}
	return diff;
};

float& Matrix3::operator()(int i, int j) {
	return mat[i][j];
};

const float& Matrix3::operator()(int i, int j) const {
	return mat[i][j];
};

Matrix3 Matrix3::operator*(Matrix3 const& mat2) const {

	// we dont have to check dimensions for now both are 3x3
	Matrix3 res;
	for (int i{ 0 }; i < 3; i++) {
		for (int j{ 0 }; j < 3; j++) {
			res(i, j) = row(i).dot(col(j));
		}
	}

	return res;
};

Vec3 Matrix3::col(int idx) const {
	Vec3 column;
	column.x = mat[0][idx];
	column.y = mat[1][idx];
	column.z = mat[2][idx];

	return column;
};

Vec3 Matrix3::row(int idx) const {
	Vec3 row;
	row.x = mat[idx][0];
	row.y = mat[idx][1];
	row.z = mat[idx][2];

	return row;
};

Matrix3& Matrix3::operator+=(const Matrix3& mat2) {

	for (int i{ 0 }; i < 3; i++) {
		for (int j{ 0 }; j < 3; j++) {
			mat[i][j] += mat2(i, j);
		}
	}
	return *this;
};

Matrix3& Matrix3::operator-=(const Matrix3& mat2) {
	for (int i{ 0 }; i < 3; i++) {
		for (int j{ 0 }; j < 3; j++) {
			mat[i][j] -= mat2(i, j);
		}
	}
	return *this;
};

Matrix3 Matrix3::operator*(float scalar) const {
	Matrix3 res;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			res(i, j) = mat[i][j] * scalar;
	return res;
};

Matrix3 operator*(float scalar, const Matrix3& mat) {
	Matrix3 res;
	return res * scalar;
};

Matrix3 Matrix3::transpose() {
	Matrix3 res;
	for (int i{ 0 }; i < 3; i++) {
		for (int j{ 0 }; j < 3; j++) {
			res(i, j) = mat[j][i];
		}
	}
	return res;
}

float Matrix3::det() {

	float v1 = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]);
	float v2 = mat[0][1] * (mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2]);
	float v3 = mat[0][2] * (mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1]);

	return v1 - v2 + v3;
};

//Matrix3 Matrix3::operator*(const Matrix3& mat, float const& scalar) {};