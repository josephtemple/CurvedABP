// core.cpp
//
// v1. 2D euclidean space
// 
// controls dynamics calculations for active brownian motion
// 1. stochastic rotational diffusion 

#include <cmath>
#include <random>
#include <vector>
//#include <vec2.h>

constexpr double PI = 3.14159265358979323846;

struct Vec2 {
    double x, y;

    // constructors
    Vec2() : x(0), y(0) {}
    Vec2(double x_, double y_) : x(x_), y(y_) {}

    // operators on vectors
    Vec2 operator+(const Vec2& other) const {
        return {x + other.x, y + other.y};
    }
    Vec2 operator-(const Vec2& other) const {
        return {x - other.x, y - other.y};
    }
    Vec2 operator*(const double a) const {
        return {x * a, y * a};
    }
    double norm2() const {
        return x*x + y*y;
    }
    double norm() const {
        return sqrt(x*x + y*y);
    }
};

struct Vec2View {
    double& x;
    double& y;

    // vector ops
    double norm() const {
        return std::sqrt(x*x + y*y);
    }

    // vector operators on Vec2View
    Vec2View& operator+=(const Vec2& v) {
        x += v.x;
        y += v.y;
        return *this;
    }
    Vec2View& operator-=(const Vec2& v) {
        x -= v.x;
        y -= v.y;
        return *this;
    }
    Vec2View& operator*=(double a) {
        x *= a;
        y *= a;
        return *this;
    }

    operator Vec2() const {
        return {x, y};
    }
};


// Heuristic gradient of potential function from V = A(x^4 - x^2 + y^4 - y^2)
Vec2 gradV(double x, double y, double A) {
    // grad(V) = A[4x^3 - 2x, 4y^3 - 2y]
    return Vec2(A*x * (4*x*x - 2), A*y * (4*y*y - 2));
}

// Data structure for the state of all the particles
struct ParticleState {
    int N;         
    std::vector<double> x; 
    std::vector<double> y;
    std::vector<double> theta; 

    // constructor, sets particule number N and allocates empty arrays of length N for x,y,theta
    ParticleState(int N_)
        : N(N_), x(N_), y(N_), theta(N_) {}

    // make available a way to do vector math on the coordinates, but keep SoA memory layout
    Vec2View pos(int i) {
        return Vec2View{ x[i], y[i] };
    }
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
    unsigned int seed;                // random number generation seed
};

// Data structure for packaging the whole simulation
struct Simulation {
    ParticleState state;
    SimParams params;
    std::mt19937 rng;
    std::size_t step_index = 0;

    Simulation(int N, const SimParams& p)
        : state(N), params(p), rng(p.seed) {}
};


// Step particle in time according to 
void step(Simulation& sim) {
    // set random seed and define probability distribution for brownian motion
    auto& state = sim.state;
    auto& params = sim.params;
    auto& rng = sim.rng;
    
    std::normal_distribution<double> noise(0.0, 1.0);

    // update particle trajectories
    for (int i = 0; i < state.N; ++i) {
        // rotational diffusion
        state.theta[i] += std::sqrt(2.0 * params.diffusion * params.dt) * noise(rng);
        
        Vec2View pi = state.pos(i);
        // collision detection
        for (int j = i + 1; j < state.N; ++j) {
            Vec2View pj = state.pos(j);

            Vec2 separation_vec = Vec2(pj) - Vec2(pi); 
            auto separation = separation_vec.norm();

            // if particle centers are closer than twice the radius, they are intersecting, so push them apart
            if (separation < 2*params.radius && separation > 1e-10) {
                // flip particle's heading directions
                state.theta[i] = -state.theta[i];
                state.theta[j] = -state.theta[j];

                // move particles apart 
                auto n = separation_vec * (1 / separation);  // unit vector along separation

                auto dist_to_move = 2*params.radius - separation; // particles move apart so their boundaries don't intersect
                Vec2 disp = n * dist_to_move * 0.5 ; // displacement vector for each is along n, half of total displacment

                pi -= disp;
                pj += disp;
            }

        // if edge of particle tries to leave the box, reflect it back
        if (state.x[i] - params.radius < -params.box_length) {
            state.x[i] = -params.box_length + params.radius;
            state.theta[i] = PI - state.theta[i];       // (vx, vy) -> (-vx, vy) 
        }
        if (state.x[i] + params.radius > params.box_length) {
            state.x[i] = params.box_length - params.radius;
            state.theta[i] = PI - state.theta[i] ;
        }
        if (state.y[i] - params.radius < -params.box_length) {
            state.y[i] = -params.box_length + params.radius;
            state.theta[i] = -state.theta[i];             // (vx, vy) -> (vx, -vy) 
        }
        if (state.y[i] + params.radius > params.box_length) {
            state.y[i] = params.box_length - params.radius;
            state.theta[i] = -state.theta[i];
        }
        }

        // update positions with external field along heading direction
        Vec2 grad_V = gradV(state.x[i], state.y[i], params.potential_strength);
        Vec2 unit_theta = Vec2(std::cos(state.theta[i]), std::sin(state.theta[i]));

        pi += (unit_theta * params.v - grad_V * params.mobility) * params.dt;
    }

    ++sim.step_index;
};
