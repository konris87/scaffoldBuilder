#pragma once
#ifndef QUATERNION_H
#define QUATERNION_H

class Vec3 {

public:
	float x{ 0 };
	float y{ 0 };
	float z{ 0 };

	Vec3();

	Vec3(float x, float y, float z);

	Vec3 normalize();

	float norm();

	Vec3 operator-(Vec3 const& vec2) const;

	Vec3 operator+(Vec3 const& vec2) const;

	friend Vec3 operator*(const Vec3& v, float const& scalar);

	friend Vec3 operator*(float const& scalar, const Vec3& v);

	float dot(Vec3 v2);

	Vec3 cross(Vec3 v2);

};

struct Mat4 {

	float r00{ 1.0 };
	float r01{ 0.0 };
	float r02{ 0.0 };
	float r03{ 0.0 };
	float r10{ 0.0 };
	float r11{ 1.0 };
	float r12{ 0.0 };
	float r13{ 0.0 };
	float r20{ 0.0 };
	float r21{ 0.0 };
	float r22{ 1.0 };
	float r23{ 0.0 };
	float r30{ 0.0 };
	float r31{ 0.0 };
	float r32{ 0.0 };
	float r33{ 1.0 };
};

class Quaternion {

public:

	// default Constructor
	Quaternion();

	// constructor from values
	Quaternion(float x, float y, float z, float w);

	// constructor from axis angle
	Quaternion(Vec3 axis, float angle);

	// overload * operator
	Quaternion operator*(const Quaternion& q2) const;

	// overload - operator
	Quaternion Quaternion::operator-(const Quaternion& q2) const;

	// overload + operator
	Quaternion Quaternion::operator+(const Quaternion& q2) const;

	// invertion
	Quaternion inverse();

	Mat4 to_matrix();

	void normalize();

	float norm();

	float x{ 0.0f };
	float y{ 0.0f };
	float z{ 0.0f };
	float w{ 1.0f };
};

Quaternion get_from_v1_to_v2(const Vec3& v1, const Vec3& v2);

#endif
