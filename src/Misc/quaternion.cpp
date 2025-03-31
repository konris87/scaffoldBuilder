#include "Misc/Quaternion.h"
#include <cmath>
#include <iostream>

Vec3::Vec3() {
	x = 0.0f;
	y = 0.0f;
	z = 0.0f;
};

Vec3::Vec3(float x, float y, float z) : x(x), y(y), z(z) {};

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

float Vec3::norm() {

	return std::sqrt(x * x + y * y + z * z);

};

Vec3 Vec3::cross(Vec3 v2) {

	Vec3 cross;

	cross.x = this->y * v2.z - this-> z * v2.y;
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

float Vec3::dot(Vec3 v2) {

	float dotp = this->x * v2.x + this->y * v2.y + this->z * v2.z;

	return dotp;
};



// ------------------------------------------------------------------------

Quaternion::Quaternion() {
	x = 0.0;
	y = 0.0;
	z = 0.0;
	w = 1.0;
};

Quaternion::Quaternion(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {};

Quaternion::Quaternion(Vec3 axis, float angle) {

	// ensure the axis is normalized

	Vec3 nAxis = axis.normalize();

	float halfAngle{ angle / 2.0f };

	x = nAxis.x * sin(halfAngle);
	y = nAxis.y * sin(halfAngle);
	z = nAxis.z * sin(halfAngle);
	w = cos(halfAngle);

};

Mat4 Quaternion::to_matrix() {

	//this->normalize();

	Mat4 mat;

	mat.r00 = 1 - 2 * y * y - 2 * z * z;
	mat.r01 = 2 * x * y - 2 * w * z;
	mat.r02 = 2 * x * z + 2 * w * y;

	mat.r10 = 2 * x * y + 2 * w * z;
	mat.r11 = 1 - 2 * x * x - 2 * z * z;
	mat.r12 = 2 * y * z - 2 * w * x;

	mat.r20 = 2 * x * z - 2 * w * y;
	mat.r21 = 2 * y * z + 2 * w * x;
	mat.r22 = 1 - 2 * x * x - 2 * y * y;

	return mat;
};

void Quaternion::normalize() {
	float norm = sqrt(x * x + y * y + z * z + w * w);

	x /= norm;
	y /= norm;
	z /= norm;
	w /= norm;

}

float Quaternion::norm() {

	return sqrt(this->x * this->x + this->y * this->y + this->z * this->z + this->w * this->w);

}

Quaternion Quaternion::operator*(const Quaternion& q2) const {

	Quaternion result;

	result.x = this->w * q2.x + this->x * q2.w + this->y * q2.z - this->z * q2.y;
	result.y = this->w * q2.y - this->x * q2.z + this->y * q2.w + this->z * q2.x;
	result.z = this->w * q2.z + this->x * q2.y - this->y * q2.x + this->z * q2.w;

	result.w = this->w * q2.w - this->x * q2.x - this->y * q2.y - this->z * q2.z;

	//result.normalize();

	return result;

};

Quaternion Quaternion::operator-(const Quaternion& q2) const {

	Quaternion result;

	result.x = this->x - q2.x;
	result.y = this->y - q2.y;
	result.z = this->z - q2.z;
	result.w = this->w - q2.w;


	return result;
};

Quaternion Quaternion::operator+(const Quaternion& q2) const {
	
	Quaternion result;

	result.x = this->x + q2.x;
	result.y = this->y + q2.y;
	result.z = this->z + q2.z;
	result.w = this->w + q2.w;


	return result;
}

Quaternion Quaternion::inverse() {
	
	Quaternion q = Quaternion(this->x * -1, this->y * -1, this->z * -1, this->w);

	return q;
};

Quaternion get_from_v1_to_v2(const Vec3& v1, const Vec3& v2) {

	Vec3 u1 = v1;
	Vec3 u2 = v2;

	u1.normalize();
	u2.normalize();

	float dotp = u1.dot(u2);

	// first check if the two are equal
	if (dotp > 0.9999f) {
		return Quaternion(0.0, 0.0, 0.0, 1.0);
	}
	
	// then if opposite
	if (dotp < -0.9999f) {

		// estimate the perpendicular vector
		Vec3 ortho = Vec3(1.0, 0.0, 0.0).cross(u1);
		
		// if this has a norm close to zero then try
		// with y
		if (ortho.norm() < 0.000001f) {
			ortho = Vec3(0, 1, 0).cross(u1);
		}

		// normalize and return quaternion representing
		// 180.0 of rotation
		ortho.normalize();
		return Quaternion(ortho, 3.14159265f);
	}

	// otherwise proceed as usually

	float angle = acos(u1.dot(v2));

	//std::cout << angle << std::endl;

	// axis of rotation is the cross product between these two
	Vec3 axis = u1.cross(v2);
	axis.normalize();

	// convert it to quaternion from axis angle

	Quaternion q = Quaternion(axis, angle);
	q.normalize();

	return q;
};
