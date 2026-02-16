#include "vec2.h"

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
