#include "vec.h"

// Vec2 implementation
Vec2::Vec2() : x(0), y(0) {}
Vec2::Vec2(double x_, double y_) : x(x_), y(y_) {}

Vec2 Vec2::operator+(const Vec2& other) const {
    return {x + other.x, y + other.y};
}
Vec2 Vec2::operator-(const Vec2& other) const {
    return {x - other.x, y - other.y};
}
Vec2 Vec2::operator*(const double a) const {
    return {x * a, y * a};
}

double Vec2::norm2() const {
    return x*x + y*y;
}
double Vec2::norm() const {
    return std::sqrt(x*x + y*y);
}

// Vec2View implementation (basically same diff but with addresses)
// vector ops
double Vec2View::norm() const {
    return std::sqrt(x*x + y*y);
}

// vector operators on Vec2View
Vec2View& Vec2View::operator+=(const Vec2& v) {
    x += v.x;
    y += v.y;
    return *this;
}
Vec2View& Vec2View::operator-=(const Vec2& v) {
    x -= v.x;
    y -= v.y;
    return *this;
}
Vec2View& Vec2View::operator*=(double a) {
    x *= a;
    y *= a;
    return *this;
}

Vec2View::operator Vec2() const {
    return {x, y};
}


// Vec3 implementation
Vec3::Vec3() : x(0), y(0), z(0) {}
Vec3::Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

Vec3 Vec3::operator+(const Vec3& other) const {
    return {x + other.x, y + other.y, z + other.z};
}
Vec3 Vec3::operator-(const Vec3& other) const {
    return {x - other.x, y - other.y, z - other.z};
}
Vec3 Vec3::operator*(const double a) const {
    return {x * a, y * a, z * a};
}

double Vec3::norm2() const {
    return x*x + y*y + z*z;
}
double Vec3::norm() const {
    return std::sqrt(x*x + y*y + z*z);
}

// Vec3View implementation (basically same diff but with addresses)
// vector ops
double Vec3View::norm() const {
    return std::sqrt(x*x + y*y + z*z);
}

// vector operators on Vec2View
Vec3View& Vec3View::operator+=(const Vec3& v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}
Vec3View& Vec3View::operator-=(const Vec3& v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}
Vec3View& Vec3View::operator*=(double a) {
    x *= a;
    y *= a;
    z *= a;
    return *this;
}

Vec3View::operator Vec3() const {
    return {x, y, z};
}
