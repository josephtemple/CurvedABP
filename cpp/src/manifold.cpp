// manifold.cpp
// note that the potential functions for each manifold were chosen to have two peaks and two wells
// with a maximum value of A
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

// V = Asin^2(\theta)cos(2\phi) -> four wells   
Vec2 SphereManifold::gradV(double q1, double q2, double A) const {
    double dVdq1 = -2*A * std::sin(q2)*std::sin(q2) * std::sin(2*q1);
    double dVdq2 =  2*A * std::sin(q2)*std::cos(q2) * std::cos(2*q1);
    return Vec2(dVdq1, dVdq2);
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

// V = Acos(theta)cos(phi)
Vec2 TorusManifold::gradV(double q1, double q2, double A) const {
    double dVdq1 = -A * std::sin(q1)*std::cos(q2);
    double dVdq2 = -A * std::sin(q2)*std::cos(q1);
    return Vec2(dVdq1, dVdq2);
}


// euclidean manifold
EuclideanManifold::EuclideanManifold() {}
double EuclideanManifold::g11(double q1, double q2) const { return 1; }
double EuclideanManifold::g22(double q1, double q2) const { return 1; }

double EuclideanManifold::vielbein1(double q1, double q2) const { return 1; }
double EuclideanManifold::vielbein2(double q1, double q2) const { return 1; }

double EuclideanManifold::connection(double q1, double q2, double dq1, double dq2) const { return 0; }

// V = A(xy/L^2) exp(1 - (x^2+y^2)/2L^2) -- four wells 
Vec2 EuclideanManifold::gradV(double q1, double q2, double A) const {
    double L = 2*sqrt(PI);

    double dVdq1 = A*q2/(L*L) * (1 + q1*q1/(L*L))*exp(1 - (q1*q1 + q2*q2)/(2*L*L));
    double dVdq2 = A*q1/(L*L) * (1 + q2*q2/(L*L))*exp(1 - (q1*q1 + q2*q2)/(2*L*L));
    return Vec2(dVdq1, dVdq2);
}