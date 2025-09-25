#ifndef VEC_H
#define VEC_H

#include <iostream>
#include <array>

class Vec3 {

public:
	float x{ 0 };
	float y{ 0 };
	float z{ 0 };

	Vec3();

	Vec3(float x, float y, float z);

	Vec3(const std::array<float, 3>& arr);

	Vec3 normalize();

	Vec3 Vec3::normalized() const;

	float norm() const;

	Vec3 operator-(Vec3 const& vec2) const;

	Vec3 operator+(Vec3 const& vec2) const;

	Vec3 operator-() const;

	friend Vec3 operator*(const Vec3& v, float const& scalar);

	friend Vec3 operator*(float const& scalar, const Vec3& v);

	Vec3 operator/(float scalar) const;

	int& operator[](int index) const;

	float dot(Vec3 v2) const;

	Vec3 cross(Vec3 v2) const;

	bool equal(const Vec3 rhs, float epsilon) const;

private:
	int* ptr;
};

std::ostream& operator << (std::ostream& out, const Vec3& v);

#endif