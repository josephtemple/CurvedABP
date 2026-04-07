// manifold.cpp
#include "manifold.h"

// sphere manifold
SphereManifold::SphereManifold(double R_) : R(R_) {}

// q1 = phi (azimuthal angle), q2 = theta (polar angle)
double SphereManifold::g11(double q1, double q2) const { return R*R * std::sin(q2)*std::sin(q2); }
double SphereManifold::g22(double q1, double q2) const { return R*R; }

double SphereManifold::vielbein1(double q1, double q2) const { return 1.0 / (R * std::sin(q2)); }
double SphereManifold::vielbein2(double q1, double q2) const { return 1.0 / R; }

double SphereManifold::connection(double q1, double q2, double dq1, double dq2) const {
    return std::cos(q2) * dq1;
}


// torus manifold
TorusManifold::TorusManifold(double R_, double r_) : R(R_), r(r_) {}

double TorusManifold::g11(double q1, double q2) const { double f = R + r*std::cos(q2); return f*f; }
double TorusManifold::g22(double q1, double q2) const { return r*r; }

double TorusManifold::vielbein1(double q1, double q2) const { return 1.0 / (R + r*std::cos(q2)); }
double TorusManifold::vielbein2(double q1, double q2) const { return 1.0 / r; }

double TorusManifold::connection(double q1, double q2, double dq1, double dq2) const {
    return -std::sin(q2) * dq1;
}