#pragma once


template<typename T> class Vec3
{
public:
	T x, y, z;
	Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
	Vec3(T a, T b, T c) : x(a), y(b), z(c) {}
	T length() {
		return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	}
	void Print_self();

	void ReConstruct(T a, T b, T c){
		x = a; y = b, z = c;
	};

	//标量*向量
	Vec3<T> operator * (const T& f) const { return Vec3<T>(x * f, y * f, z * f); }

	//向量/标量
	Vec3<T> operator / (const T& t) { return Vec3<T>(x / t, y / t, z / t); }

	//向量1*向量2
	Vec3<T> operator * (const Vec3<T>& v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }

	//向量1/向量2
	Vec3<T> operator / (const Vec3<T>& v) const { return Vec3<T>(x / v.x, y / v.y, z / v.z); }

	//向量1-向量2
	Vec3<T> operator - (const Vec3<T>& v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }

	//向量1+向量2
	Vec3<T> operator + (const Vec3<T>& v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }

};

void Vec3<double>::Print_self()
{
	printf("%f  %f  %f\n", x, y, z);
}
void Vec3<int>::Print_self()
{
	printf("%d  %d  %d\n", x, y, z);
}

template<typename T> T Vec3_Dot(Vec3<T> A, Vec3<T> B) {
	return A.x * B.x + A.y * B.y + A.z * B.z;
}

template<typename T> Vec3<T> Vec3_Cross(Vec3<T> A, Vec3<T> B) {
	Vec3<T> C;
	C.x = A.y * B.z - A.z * B.y;
	C.y = A.z * B.x - A.x * B.z;
	C.z = A.x * B.y - A.y * B.x;
	return C;
}

// 二维
template<typename T> class Vec2
{
public:
	T x, y;
	Vec2() : x(T(0)), y(T(0)) {}
	Vec2(T a, T b) : x(a), y(b) {}
	T length() {
		return sqrt(pow(x, 2) + pow(y, 2));
	}
	void Print_self();

	void ReConstruct(T a, T b) {
		x = a; y = b;
	};

	//标量*向量
	Vec2<T> operator * (const T& f) const { return Vec2<T>(x * f, y * f); }

	//向量/标量
	Vec2<T> operator / (const T& t) { return Vec2<T>(x / t, y / t); }

	//向量1*向量2
	Vec2<T> operator * (const Vec2<T>& v) const { return Vec2<T>(x * v.x, y * v.y); }

	//向量1/向量2
	Vec2<T> operator / (const Vec2<T>& v) const { return Vec2<T>(x / v.x, y / v.y); }

	//向量1-向量2
	Vec2<T> operator - (const Vec2<T>& v) const { return Vec2<T>(x - v.x, y - v.y); }

	//向量1+向量2
	Vec2<T> operator + (const Vec2<T>& v) const { return Vec2<T>(x + v.x, y + v.y); }

};
void Vec2<double>::Print_self()
{
	printf("%f  %f  %f\n", x, y);
}
void Vec2<int>::Print_self()
{
	printf("%d  %d  %d\n", x, y);
}
