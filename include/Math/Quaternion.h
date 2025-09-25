#ifndef QUATERNION_H
#define QUATERNION_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "Vec.h"

class Mat4 {

public:

	~Mat4() {};
	Mat4() {
		for (int c = 0; c < 4; ++c)
			for (int r = 0; r < 4; ++r)
				m[c][r] = (c == r) ? 1.0f : 0.0f;
	}

	glm::mat4 to_glm() const {
		return glm::mat4(
			m[0][0], m[1][0], m[2][0], m[3][0],  // column 0
			m[0][1], m[1][1], m[2][1], m[3][1],  // column 1
			m[0][2], m[1][2], m[2][2], m[3][2],  // column 2
			m[0][3], m[1][3], m[2][3], m[3][3]   // column 3
		);
	}

	// Optional: Add double operator[] support
	struct ColumnProxy {
		float* colData;
		ColumnProxy(float* data) : colData(data) {}
		float& operator[](int row) { return colData[row]; }
		const float& operator[](int row) const { return colData[row]; }
	};

	ColumnProxy operator[](int col) { return ColumnProxy(m[col]); }
	const ColumnProxy operator[](int col) const { return ColumnProxy((float*)m[col]); }

private:
	float& operator()(int col, int row) { return m[col][row]; }
	const float& operator()(int col, int row) const { return m[col][row]; }
	float m[4][4];
};



class Quaternion {

public:

	// default Constructor
	Quaternion();

	// constructor from values
	Quaternion(float x, float y, float z, float w);

	// constructor from axis angle
	Quaternion(Vec3 axis, float angle);

	// constructor from rotation matrix
	Quaternion(Mat4 matrix);

	Quaternion(const glm::mat4& matrix);

	// overload * operator
	Quaternion operator*(const Quaternion& q2) const;

	// overload - operator
	Quaternion Quaternion::operator-(const Quaternion& q2) const;

	// overload + operator
	Quaternion Quaternion::operator+(const Quaternion& q2) const;

	// invertion
	Quaternion inverse() const;

	Mat4 to_matrix() const;

	Mat4 to_matrix_test() const;

	void normalize();

	float norm();

	float x{ 0.0f };
	float y{ 0.0f };
	float z{ 0.0f };
	float w{ 1.0f };
};

Quaternion get_from_v1_to_v2(const Vec3& v1, const Vec3& v2);

Quaternion get_v1_v2_test(const Vec3& v1, const Vec3& v2);

#endif
