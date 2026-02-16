// Header for Vec2 and Vec2View structs
#pragma once
#include <cmath>

struct Vec2 {
    double x, y;

    // constructors
    Vec2();
    Vec2(double x_, double y_);

    // operators on vectors
    Vec2 operator+(const Vec2& other) const;
    Vec2 operator-(const Vec2& other) const;
    Vec2 operator*(const double a) const;

    double norm2() const;
    double norm() const;
};

struct Vec2View {
    double& x;
    double& y;
    // vector operators on Vec2View
    Vec2View& operator+=(const Vec2& v);
    Vec2View& operator-=(const Vec2& v);
    Vec2View& operator*=(double a);

    double norm() const;

    operator Vec2() const;
};
