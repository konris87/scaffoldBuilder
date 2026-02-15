#ifndef VEC_H
#define VEC_H

#include <iostream>
#include <array>

#include "Eigen/Dense"

class Vec3 {

public:
	union {
		// Option A: Named variables
		struct { float x, y, z; };
		// Option B: Array
		float data[3];
	};

	Vec3();

	Vec3(float x, float y, float z);

	Vec3(const std::array<float, 3>& arr);

	Vec3(const Eigen::Vector3f& arr);

	Vec3 normalize();

	Vec3 Vec3::normalized() const;

	float norm() const;

	Vec3 operator-(Vec3 const& vec2) const;

	Vec3 operator+(Vec3 const& vec2) const;

	Vec3 operator-() const;

	Vec3& operator*=(const float& scalar);

	Vec3& operator+=(const Vec3& v);

	friend Vec3 operator*(const Vec3& v, float const& scalar);

	friend Vec3 operator*(float const& scalar, const Vec3& v);

	Vec3 operator/(float scalar) const;

	// for reading
	float& operator[](int index);

	// for writing
	const float& operator[](int index) const;

	float dot(Vec3 v2) const;

	Vec3 cross(Vec3 v2) const;

	bool equal(const Vec3 rhs, float epsilon) const;

private:
	int* ptr;
};

std::ostream& operator << (std::ostream& out, const Vec3& v);


// ------------------------------------------------------------------
// 2d vector class 

class Vec2 {

public:
	float x{ 0 };
	float y{ 0 };

	Vec2();

	Vec2(float x, float y);

	Vec2(const std::array<float, 2>& arr);

	Vec2(const Eigen::Vector2f& arr);

	Vec2 normalize();

	Vec2 Vec2::normalized() const;

	float norm() const;

	Vec2 operator-(Vec2 const& vec2) const;

	Vec2 operator+(Vec2 const& vec2) const;

	Vec2 operator-() const;

	Vec2& operator*=(const float& scalar);

	Vec2& operator+=(const Vec2& v);

	friend Vec2 operator*(const Vec2& v, float const& scalar);

	friend Vec2 operator*(float const& scalar, const Vec2& v);

	Vec2 operator/(float scalar) const;

	int& operator[](int index) const;

	float dot(Vec2 v2) const;

	bool equal(const Vec2 rhs, float epsilon) const;

private:
	int* ptr;
};

std::ostream& operator << (std::ostream& out, const Vec2& v);

#endif