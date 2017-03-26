#ifndef __VECTOR_H__
#define __VECTOR_H__

#pragma once

class Vector3f
{
public:
	float		data_[3];
	//float		&x_,&y_,&z_;
public:
	Vector3f(void);
	~Vector3f(void);
	Vector3f(float x, float y, float z);
	Vector3f(float data[]);
	Vector3f(const Vector3f &v);
	void ones();
	void e(int i);

	inline float& x();
	inline float& y();
	inline float& z();
	/// Deprecated, not really reqd anymore
	//void setX(float x);
	//void setY(float y);
	//void setZ(float z);
	inline float& operator() (int i);

	float length();
	float lengthSquared();
	void operator=(const Vector3f& v);
	//Vector3f operator+(const Vector3f& v);
	//Vector3f operator-(const Vector3f& v);
	//float operator*(const Vector3f& v);
	//Vector3f operator*(float a);
	friend inline Vector3f operator+(const Vector3f &v1, const Vector3f &v2);
	friend inline Vector3f operator-(const Vector3f &v1, const Vector3f &v2);
	friend inline float operator*(const Vector3f &v1, const Vector3f &v2);
	friend inline Vector3f operator*(const Vector3f &v1, const float &a);
	friend inline Vector3f operator/(const Vector3f &v1, const float &a);
	friend inline Vector3f operator/(const Vector3f &v1, const Vector3f &v2);
	friend inline Vector3f operator*(const float &a, const Vector3f &v1);
	float dot(const Vector3f& v);
	inline Vector3f cross(const Vector3f& v);
	inline bool isNan();
};

inline Vector3f Vector3f::cross(const Vector3f& v) {
	return Vector3f(data_[1] * v.data_[2] - data_[2] * v.data_[1], data_[2] * v.data_[0] - data_[0] * v.data_[2], data_[0] * v.data_[1] - data_[1] * v.data_[0] );
}

inline void saxpby(Vector3f& dest, Vector3f& x, float a, Vector3f& y, float b) {
	dest=Vector3f(a * x.x() + b * y.x(), a * x.y() + b * y.y(), a * x.z() + b * y.z());
}

inline void isaxpby(Vector3f& x, float a, Vector3f& y, float b) {
	x = Vector3f(a * x.x() + b * y.x(), a * x.y() + b * y.y(), a * x.z() + b * y.z());
}

inline Vector3f saxpby(Vector3f& x, float a, Vector3f& y, float b) {
	return Vector3f(a * x.x() + b * y.x(), a * x.y() + b * y.y(), a * x.z() + b * y.z());
}

inline Vector3f operator+(const Vector3f &v1, const Vector3f &v2) {
	return Vector3f(v1.data_[0] + v2.data_[0], v1.data_[1] + v2.data_[1], v1.data_[2] + v2.data_[2]);
}

inline Vector3f operator-(const Vector3f &v1, const Vector3f &v2) {
	return Vector3f(v1.data_[0] - v2.data_[0], v1.data_[1] - v2.data_[1], v1.data_[2] - v2.data_[2]);
}

inline float operator*(const Vector3f &v1, const Vector3f &v2) {
	return (v1.data_[0] * v2.data_[0] + v1.data_[1] * v2.data_[1] + v1.data_[2] * v2.data_[2]);
}

inline Vector3f operator+(const Vector3f &v1, const float &b) {
	return Vector3f(v1.data_[0] * b, v1.data_[1] * b, v1.data_[2] * b);
}

inline Vector3f operator+(const float &b, const Vector3f &v1) {
	return Vector3f(v1.data_[0] * b, v1.data_[1] * b, v1.data_[2] * b);
}

inline float& Vector3f::x() {
	return data_[0];
}

inline float& Vector3f::y() {
	return data_[1];
}

inline float& Vector3f::z() {
	return data_[2];
}

inline float& Vector3f::operator () (int i) {
	return data_[i];
}

inline Vector3f operator*(const Vector3f &v1, const float &a) {
	return Vector3f( v1.data_[0] * a, v1.data_[1] * a, v1.data_[2] * a);
}

inline Vector3f operator*(const float &a, const Vector3f &v1) {
	return Vector3f( v1.data_[0] * a, v1.data_[1] * a, v1.data_[2] * a);
}

inline Vector3f operator/(const Vector3f &v1, const float &a) {
	return Vector3f( v1.data_[0] / a, v1.data_[1] / a, v1.data_[2] / a);
}

inline Vector3f operator/(const Vector3f &v1, const Vector3f &v2) {
	return Vector3f( v1.data_[0] / v2.data_[0], v1.data_[1] / v2.data_[1], v1.data_[2] / v2.data_[2]);
}

inline bool Vector3f::isNan() {
	return (data_[0] != data_[0]) || (data_[1] != data_[1]) || (data_[2] != data_[2]);
}
#endif