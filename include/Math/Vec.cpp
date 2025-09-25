#include "Vec.h"
#include <cmath>

Vec3::Vec3() {
	x = 0.0f;
	y = 0.0f;
	z = 0.0f;
};

Vec3::Vec3(const std::array<float, 3>& arr) {
	x = arr[0];
	y = arr[1];
	z = arr[2];
}

Vec3::Vec3(float x, float y, float z) : x(x), y(y), z(z) {};

Vec3 Vec3::operator-() const {
	return Vec3(-x, -y, -z);
};

Vec3 Vec3::operator-(Vec3 const& vec2) const {

	Vec3 summed;
	summed.x = x - vec2.x;
	summed.y = y - vec2.y;
	summed.z = z - vec2.z;
	return summed;
};

Vec3 Vec3::operator+(Vec3 const& vec2) const {

	Vec3 summed;
	summed.x = x + vec2.x;
	summed.y = y + vec2.y;
	summed.z = z + vec2.z;
	return summed;
};

Vec3 operator*(const Vec3& v, float const& scalar) {

	Vec3 scalarProd;
	scalarProd.x = scalar * v.x;
	scalarProd.y = scalar * v.y;
	scalarProd.z = scalar * v.z;

	return scalarProd;
};

Vec3 operator*(float const& scalar, const Vec3& v) {

	Vec3 scalarProd;
	scalarProd.x = scalar * v.x;
	scalarProd.y = scalar * v.y;
	scalarProd.z = scalar * v.z;

	return scalarProd;
};

Vec3 Vec3::operator/(float scalar) const {
	return Vec3(x / scalar, y / scalar, z / scalar);
}

int& Vec3::operator[](int index) const {
	
	if (index >= 3) {
		std::cout << "Array Index out ouf bound" << std::endl;
		return ptr[0];
	}

	return ptr[index];

};

std::ostream& operator << (std::ostream& out, const Vec3& v) {

	out << v.x << " " << v.y << " " << v.z;

	return out;
};

float Vec3::norm() const {

	return std::sqrt(x * x + y * y + z * z);

};

Vec3 Vec3::cross(Vec3 v2) const {

	Vec3 cross;

	cross.x = this->y * v2.z - this->z * v2.y;
	cross.y = this->z * v2.x - this->x * v2.z;
	cross.z = this->x * v2.y - v2.x * this->y;

	return cross;

};

Vec3 Vec3::normalize() {

	float norm = this->norm();

	if (norm > 0) {
		x /= norm;
		y /= norm;
		z /= norm;
	}
	return *this;
}

Vec3 Vec3::normalized() const {
	float n = this->norm();
	if (n > 0.0f) {
		return Vec3(x / n, y / n, z / n);
	}
	return Vec3(); // zero vector
};

float Vec3::dot(Vec3 v2) const {

	float dotp = this->x * v2.x + this->y * v2.y + this->z * v2.z;

	return dotp;
};

bool Vec3::equal(const Vec3 rhs, float epsilon) const {

	return fabs(x - rhs.x) < epsilon && fabs(y - rhs.y) < epsilon;

};
