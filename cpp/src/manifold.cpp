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

// keep q1 = phi between (0,2pi) and q2 = theta between (0,pi)
// flipping phi if theta exceeds its bounds
void SphereManifold::wrap(double& q1, double& q2) const {
    // reflect over north pole
    if (q2 < 0) {
        q2 = -q2;
        q1 += M_PI;
    }
    // reflect over south pole
    if (q2 > M_PI) {
        q2 = 2.0 * M_PI - q2;
        q1 += M_PI;
    }
    
    // keep phi between (0,2pi)
    q1 = std::fmod(q1, 2.0 * PI);
    if (q1 < 0) {q1 += 2.0 * PI;}
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

// keep q1 = phi and q2 = theta between (0,2pi)
void TorusManifold::wrap(double& q1, double& q2) const {
    q1 = std::fmod(q1, 2.0 * PI);
    if (q1 < 0) {q1 += 2.0 * PI;}

    q2 = std::fmod(q2, 2.0 * PI);
    if (q2 < 0) {q2 += 2.0 * PI;}
}


// euclidean manifold
EuclideanManifold::EuclideanManifold() {}
double EuclideanManifold::g11(double q1, double q2) const { return 1; }
double EuclideanManifold::g22(double q1, double q2) const { return 1; }

double EuclideanManifold::vielbein1(double q1, double q2) const { return 1; }
double EuclideanManifold::vielbein2(double q1, double q2) const { return 1; }

double EuclideanManifold::connection(double q1, double q2, double dq1, double dq2) const { return 0; }

// keep particles in +/- 2 sqrt(pi)  (same surface area as R = 1 sphere), pacman style
void EuclideanManifold::wrap(double& q1, double& q2) const {
    double L = 2*std::sqrt(PI);

    if (q1 >  L/2) q1 -= L;
    if (q1 < -L/2) q1 += L;
    
    if (q2 >  L/2) q2 -= L;
    if (q2 < -L/2) q2 += L;
}