// manifold.h
#pragma once
#include <cmath>
#include "vec.h"

constexpr double PI = 3.14159265358979323846;

// abstract base manifold struct
struct Manifold {
    virtual double vielbein1(double q1, double q2) const = 0;
    virtual double vielbein2(double q1, double q2) const = 0;
    virtual double g11(double q1, double q2) const = 0;
    virtual double g22(double q1, double q2) const = 0;
    virtual double connection(double q1, double q2, double dq1, double dq2) const = 0;
    virtual Vec2 gradV(double q1, double q2, double A) const = 0;
    virtual ~Manifold() = default;
};

// sphere manifold [q1 = phi (azimuthal angle), q2 = theta (polar angle)]
struct SphereManifold : Manifold {
    double R;
    SphereManifold(double R_);

    double vielbein1(double q1, double q2) const override;
    double vielbein2(double q1, double q2) const override;
    double g11(double q1, double q2) const override;
    double g22(double q1, double q2) const override;
    double connection(double q1, double q2, double dq1, double dq2) const override;
    virtual Vec2 gradV(double q1, double q2, double A) const override;
};

// torus manifold [q1 = phi (around the donut), q2 = psi (angle around cross-section)]
struct TorusManifold : Manifold {
    double R, r;
    TorusManifold(double R_, double r_);

    double vielbein1(double q1, double q2) const override;
    double vielbein2(double q1, double q2) const override;
    double g11(double q1, double q2) const override;
    double g22(double q1, double q2) const override;
    double connection(double q1, double q2, double dq1, double dq2) const override;
    virtual Vec2 gradV(double q1, double q2, double A) const override;
};

// euclidean manifold [q1 = x, q2 = y]
struct EuclideanManifold : Manifold {
    EuclideanManifold();

    double vielbein1(double q1, double q2) const override;
    double vielbein2(double q1, double q2) const override;
    double g11(double q1, double q2) const override;
    double g22(double q1, double q2) const override;
    double connection(double q1, double q2, double dq1, double dq2) const override;
    virtual Vec2 gradV(double q1, double q2, double A) const override;
};