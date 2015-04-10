#ifndef __VEC2D_H__
#define __VEC2D_H__

#include <math.h>

template<class T>
class Vec2D {
public:
	T x, y; // The coordinates of this vector
	
	Vec2D(T x, T y): x(x), y(y) {}
	Vec2D() {}
	
	T norm2() {
		return x*x + y*y;
	}
	
	T norm() {
		return sqrt(norm2());
	}
	
	Vec2D<T> &operator*=(T s) {
		x *= s;
		y *= s;
		return *this;
	}
	
	Vec2D<T> &operator/=(T s) {
		x /= s;
		y /= s;
		return *this;
	}
	
	Vec2D<T> operator*(T s) const {
		return Vec2D<T>(s*x, s*y);
	}
	
	friend Vec2D<T> operator*(T s, const Vec2D<T> &v) {
		return Vec2D<T>(s*v.x, s*v.y);
	}

	Vec2D<T> operator/(T s) const {
		return Vec2D<T>(x/s, y/s);
	}

	Vec2D<T> operator-() const {
		return Vec2D<T>(-x, -y);
	}

	Vec2D<T> &operator+=(const Vec2D<T> &v) {
		x += v.x;
		y += v.y;
		return *this;
	}
	
	Vec2D<T> &operator-=(const Vec2D<T> &v) {
		return (*this += -v);
	}
	
	Vec2D<T> operator+(const Vec2D<T> &v) const {
		return Vec2D<T>(x + v.x, y + v.y);
	}
	
	Vec2D<T> operator-(const Vec2D<T> &v) const {
		return Vec2D<T>(x - v.x, y - v.y);
	}

	T dot(const Vec2D<T> &v) const {
		return x*v.x + y*v.y;
	}
};

#endif // __VEC2D_H__
