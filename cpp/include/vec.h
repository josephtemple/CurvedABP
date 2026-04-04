// Header for Vec2, Vec2View, Vec3, Vec3View structs
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


struct Vec3 {
    double x, y, z;

    // constructors
    Vec3();
    Vec3(double x_, double y_, double z_);

    // operators on vectors
    Vec3 operator+(const Vec3& other) const;
    Vec3 operator-(const Vec3& other) const;
    Vec3 operator*(const double a) const;

    double norm2() const;
    double norm() const;
};

struct Vec3View {
    double& x, &y, &z;
    // vector operators on Vec3View
    Vec3View& operator+=(const Vec3& v);
    Vec3View& operator-=(const Vec3& v);
    Vec3View& operator*=(double a);



    double norm() const;

    operator Vec3() const;
};
