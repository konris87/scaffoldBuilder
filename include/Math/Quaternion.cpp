#include "Quaternion.h"
#include "Vec.h"
#include <cmath>
#include <iostream>

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

Quaternion::Quaternion(Mat4 mat) {

	float trace = mat[0][0] + mat[1][1] + mat[2][2]; // I removed + 1.0f; see discussion with Ethan
	if (trace > 0) {// I changed M_EPSILON to 0
		float s = 0.5f / sqrtf(trace + 1.0f);
		w = 0.25f / s;
		x = (mat[2][1] - mat[1][2]) * s;
		y = (mat[0][2] - mat[2][0]) * s;
		z = (mat[1][0] - mat[0][1]) * s;
	}
	else {
		if (mat[0][0] > mat[1][1] && mat[0][0] > mat[2][2]) {
			float s = 2.0f * sqrtf(1.0f + mat[0][0] - mat[1][1] - mat[2][2]);
			w = (mat[2][1] - mat[1][2]) / s;
			x = 0.25f * s;
			y = (mat[0][1] + mat[1][0]) / s;
			z = (mat[0][2] + mat[2][0]) / s;
		}
		else if (mat[1][1] > mat[2][2]) {
			float s = 2.0f * sqrtf(1.0f + mat[1][1] - mat[0][0] - mat[2][2]);
			w = (mat[0][2] - mat[2][0]) / s;
			x = (mat[0][1] + mat[1][0]) / s;
			y = 0.25f * s;
			z = (mat[1][2] + mat[2][1]) / s;
		}
		else {
			float s = 2.0f * sqrtf(1.0f + mat[2][2] - mat[0][0] - mat[1][1]);
			w = (mat[1][0] - mat[0][1]) / s;
			x = (mat[0][2] + mat[2][0]) / s;
			y = (mat[1][2] + mat[2][1]) / s;
			z = 0.25f * s;
		}
	}
};

Quaternion::Quaternion(const glm::mat4& mat) {
	// Only upper-left 3x3 is relevant for rotation
	float m00 = mat[0][0], m01 = mat[1][0], m02 = mat[2][0];
	float m10 = mat[0][1], m11 = mat[1][1], m12 = mat[2][1];
	float m20 = mat[0][2], m21 = mat[1][2], m22 = mat[2][2];

	float trace = m00 + m11 + m22;

	if (trace > 0.0f) {
		float s = 0.5f / sqrtf(trace + 1.0f);
		w = 0.25f / s;
		x = (m21 - m12) * s;
		y = (m02 - m20) * s;
		z = (m10 - m01) * s;
	}
	else {
		if (m00 > m11 && m00 > m22) {
			float s = 2.0f * sqrtf(1.0f + m00 - m11 - m22);
			w = (m21 - m12) / s;
			x = 0.25f * s;
			y = (m01 + m10) / s;
			z = (m02 + m20) / s;
		}
		else if (m11 > m22) {
			float s = 2.0f * sqrtf(1.0f + m11 - m00 - m22);
			w = (m02 - m20) / s;
			x = (m01 + m10) / s;
			y = 0.25f * s;
			z = (m12 + m21) / s;
		}
		else {
			float s = 2.0f * sqrtf(1.0f + m22 - m00 - m11);
			w = (m10 - m01) / s;
			x = (m02 + m20) / s;
			y = (m12 + m21) / s;
			z = 0.25f * s;
		}
	}

	// Normalize to avoid drift
	float norm = sqrtf(w * w + x * x + y * y + z * z);
	if (norm > 1e-8f) {
		w /= norm; x /= norm; y /= norm; z /= norm;
	}
}

Mat4 Quaternion::to_matrix() const {

	Mat4 mat;

	mat[0][0] = 1 - 2 * y * y - 2 * z * z;
	mat[0][1] = 2 * x * y - 2 * w * z;
	mat[0][2] = 2 * x * z + 2 * w * y;
	mat[0][3] = 0.0f;

	mat[1][0] = 2 * x * y + 2 * w * z;
	mat[1][1] = 1 - 2 * x * x - 2 * z * z;
	mat[1][2] = 2 * y * z - 2 * w * x;
	mat[1][3] = 0.0f;

	mat[2][0] = 2 * x * z - 2 * w * y;
	mat[2][1] = 2 * y * z + 2 * w * x;
	mat[2][2] = 1 - 2 * x * x - 2 * y * y;
	mat[2][3] = 0.0f;

	mat[3][0] = 0.0f;
	mat[3][1] = 0.0f;
	mat[3][2] = 0.0f;
	mat[3][3] = 1.0f;

	return mat;
};

Mat4 Quaternion::to_matrix_test() const {

	// NOTE: assume the quaternion is unit length
	// compute common values
	float x2 = x + x;
	float y2 = y + y;
	float z2 = z + z;
	float xx2 = x * x2;
	float xy2 = x * y2;
	float xz2 = x * z2;
	float yy2 = y * y2;
	float yz2 = y * z2;
	float zz2 = z * z2;
	float sx2 = w * x2;
	float sy2 = w * y2;
	float sz2 = w * z2;

	// build 4x4 matrix (column-major) and return
	Mat4 mat;
	mat[0][0] = 1 - (yy2 + zz2);
	mat[1][0] = xy2 + sz2;
	mat[2][0] = xz2 - sy2;
	mat[3][0] = 0.0;

	mat[0][1] = xy2 - sz2;
	mat[1][1] = 1 - (xx2 + zz2);
	mat[2][1] = yz2 + sx2;
	mat[3][1] = 0.0;

	mat[0][2] = xz2 + sy2;
	mat[1][2] = yz2 - sx2;
	mat[2][2] = 1 - (xx2 + yy2);
	mat[3][2] = 0.0;
	
	mat[0][3] = 0.0;
	mat[1][3] = 0.0;
	mat[2][3] = 0.0;
	mat[3][3] = 1.0;

	return mat;
};

void Quaternion::normalize() {
	float norm = sqrt(x * x + y * y + z * z + w * w);

	if (norm > 0.00001f) {
		x /= norm;
		y /= norm;
		z /= norm;
		w /= norm;
	}
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

Quaternion Quaternion::inverse() const {
	
	float n2 = x * x + y * y + z * z + w * w;
	if (n2 > 0.000001f) {
		return Quaternion(-x / n2, -y / n2, -z / n2, w / n2);
	}
	// else return indentity
	return Quaternion(); 
};

Quaternion get_v1_v2_test(const Vec3& v1, const Vec3& v2) {
	
	const float EPSILON = 0.001f;
	const float HALF_PI = acos(-1) * 0.5f;

	// if two vectors are equal return the vector with 0 rotation
	if (v1.equal(v2, EPSILON))
	{
		return Quaternion(v1, 0);
	}
	// if two vectors are opposite return a perpendicular vector with 180 angle
	else if (v1.equal(-v2, EPSILON))
	{
		Vec3 v;
		if (v1.x > -EPSILON && v1.x < EPSILON) {
			v.x = 1;
			v.y = 0;
			v.z = 0;
		}    // if x ~= 0
		else if (v1.y > -EPSILON && v1.y < EPSILON) {
			v.x = 0;
			v.y = 1;
			v.z = 0;
		}  // if y ~= 0
		else {
			v.x = 0;
			v.y = 0;
			v.z = 1;
		}                                  // if z ~= 0
		return Quaternion(v, HALF_PI);
	}

	Vec3 u1 = v1;                    // convert to normal vector
	Vec3 u2 = v2;
	u1.normalize();
	u2.normalize();

	Vec3 v = u1.cross(u2);           // compute rotation axis
	float angle = acosf(u1.dot(u2));    // rotation angle
	return Quaternion(v, angle * 0.5f); // half angle
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
