/*
# Particle class
class Particles:
    def __init__(self, x_i, y_i, v_i, theta_i):
        self.n = len(x_i)
        self.x = x_i
        self.y = y_i

        self.v = v_i
        self.theta = theta_i

        params.radius = 0.02
        self.diffusion_rate = 5
        self.mobility = .5
        self.dt = 1/500

    def update(self):
        # brownian-ly change direction
        self.theta += np.sqrt(2 * self.diffusion_rate * self.dt) * np.random.normal(size = self.n)

        # collision detection, checking each pair of particles
        for i in range(self.n):
            for j in range(i+1, self.n):
                x_separation = x(j) - x(i)
                y_separation = y(j) - y(i)

                separation = np.sqrt(x_separation**2 + y_separation**2)
                
                if separation < 2*params.radius:
                    # flip particle velocities
                    theta(i) = -theta(i)
                    self.theta[j] = -self.theta[j]

                    # move particles apart 
                    nx = x_separation / separation 
                    ny = y_separation / separation

                    dist_to_move = 2*params.radius - separation
                    dx = dist_to_move * nx
                    dy = dist_to_move * ny

                    x(i) = x(i) - dx / 2
                    y(i) = y(i) - dy / 2

                    x(j) = x(j) + dx / 2
                    y(j) = y(j) + dy / 2

            # boundary conditions: reverse along that direction
            if x(i) - params.radius < -params.box_length:
                x(i) = -params.box_length + params.radius
                theta(i) = M_PI - theta(i) # (vx, vy) -> (-vx, vy) 
            if x(i) + params.radius > params.box_length:
                x(i) = params.box_length - params.radius
                theta(i) = M_PI - theta(i) 
            if y(i) - params.radius < -params.box_length:
                y(i) = -params.box_length + params.radius
                theta(i) = -theta(i) # (vx, vy) -> (vx, -vy) 
            if y(i) + params.radius > params.box_length:
                y(i) = params.box_length - params.radius
                theta(i) = -theta(i) 

        self.theta %= 2*M_PI

    def step(self):
        gradV = del_V(self.x, self.y)
        Fx, Fy = -self.mobility * gradV

        self.x += (self.v * np.cos(self.theta) + Fx) * self.dt
        self.y += (self.v * np.sin(self.theta) + Fy) * self.dt

    def get_position(self):
        return self.x, self.y
*/

// core.cpp
//
// controls dynamics calculations for active brownian motion

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cmath>
#include <random>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace py = pybind11;

// Data structure for the state of all the particles
struct ParticleState {
    py::array_t<double> x;        // pointer to array of x positions
    py::array_t<double> y;        // same for y
    py::array_t<double> theta;    // same for orientations
};

// Data structure for the constant parameters of the simulation
struct SimParams {
    double v;                         // velocity of each particle            
    double radius;                    // radius of each particle in params.box_length
    double diffusion;                 // rotational diffusion constant in rad^2 / time
    double mobility;                  // 
    double dt;                        // time step size
    double box_length;                // size of box to which particles are constrained
    double potential_strength;        // change how strong the potential function acts
    unsigned int seed;          // random number generation seed
};

void step(ParticleState& state, const SimParams& params) {
    // construct C++ variables from numpy arrays stored in ParticleState
    auto x = state.x.mutable_unchecked<1>();
    auto y = state.y.mutable_unchecked<1>();
    auto theta = state.theta.mutable_unchecked<1>();
    int N = x.shape(0);

    // set random seed and define probability distribution for brownian motion
    static std::mt19937 rng(params.seed);
    std::normal_distribution<double> noise(0.0, 1.0);

    // update particle trajectories
    for (int i = 0; i < N; ++i) {
        // rotational diffusion
        theta(i) += std::sqrt(2.0 * params.diffusion * params.dt) * noise(rng);

        // collision detection
        for (int j = i + 1; j < N; ++j) {
            auto x_separation = x(j) - x(i);
            auto y_separation = y(j) - y(i);

            auto separation = std::sqrt(x_separation * x_separation + y_separation * y_separation);

            // if particle centers are closer than twice the radius, they are intersecting, so push them apart
            if (separation < 2*params.radius) {
                // flip particle's heading directions
                theta(i) = -theta(i);
                theta(j) = -theta(j);

                // move particles apart 
                auto nx = x_separation / separation;
                auto ny = y_separation / separation;

                auto dist_to_move = 2*params.radius - separation;
                auto dx = dist_to_move * nx;
                auto dy = dist_to_move * ny;

                x(i) -= dx / 2;
                y(i) -= dy / 2;

                x(j) += dx / 2;
                y(j) += dy / 2;
            }

            // if edge of particle tries to leave the box, reflect it back
            if (x(i) - params.radius < -params.box_length) {
                x(i) = -params.box_length + params.radius;
                theta(i) = M_PI - theta(i);       // (vx, vy) -> (-vx, vy) 
            }
            if (x(i) + params.radius > params.box_length) {
                x(i) = params.box_length - params.radius;
                theta(i) = M_PI - theta(i) ;
            }
            if (y(i) - params.radius < -params.box_length) {
                y(i) = -params.box_length + params.radius;
                theta(i) = -theta(i);             // (vx, vy) -> (vx, -vy) 
            }
            if (y(i) + params.radius > params.box_length) {
                y(i) = params.box_length - params.radius;
                theta(i) = -theta(i);
            }
        }

        // externial field with potential V(r) = 4r^3 - 2r for both r={x,y}
        double Fx = -params.mobility * params.potential_strength * (4.0*pow(x(i),3) - 2.0*x(i));
        double Fy = -params.mobility * params.potential_strength * (4.0*pow(y(i),3) - 2.0*y(i));
        x(i) += (params.v * std::cos(theta(i)) + Fx) * params.dt;
        y(i) += (params.v * std::sin(theta(i)) + Fy) * params.dt;
    }
};

// Expose module
PYBIND11_MODULE(abp_core, m) {
    py::class_<ParticleState>(m, "ParticleState")
        .def(py::init<>())
        .def_readwrite("x", &ParticleState::x)
        .def_readwrite("y", &ParticleState::y)
        .def_readwrite("theta", &ParticleState::theta);

    py::class_<SimParams>(m, "SimParams")
        .def(py::init<>())
        .def_readwrite("v", &SimParams::v)
        .def_readwrite("radius", &SimParams::radius)
        .def_readwrite("diffusion", &SimParams::diffusion)
        .def_readwrite("mobility", &SimParams::mobility)
        .def_readwrite("dt", &SimParams::dt)
        .def_readwrite("box_length", &SimParams::box_length)
        .def_readwrite("potential_strength", &SimParams::potential_strength)
        .def_readwrite("seed", &SimParams::seed);

    m.def("step", &step);
}