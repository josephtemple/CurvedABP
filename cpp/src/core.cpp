// core.cpp
//
// v1. 2D euclidean space
// 
// controls dynamics calculations for active brownian motion
// 1. stochastic rotational diffusion 
// 2. external field, move along gradient (a la taxis)
// Those combined constitute an Euler-Maruyama discretization of the SDE
// 3. hard-sphere collision detection
// 4. application of boundary conditions
// 5. saving output to HDF5 file format


#include <cmath>
#include <random>
#include <vector>
#include <H5Cpp.h>

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
        return std::sqrt(x*x + y*y);
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

    SimParams(double v_, double r_, double diff_, double mobil_, double dt_, double boxlength_, double potenstr_, unsigned int seed_)
        : v(v_), radius(r_), diffusion(diff_), mobility(mobil_), dt(dt_), box_length(boxlength_), potential_strength(potenstr_), seed(seed_) {}
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

        // update positions with external field along heading direction
        Vec2View pi = state.pos(i);
        Vec2 grad_V = gradV(state.x[i], state.y[i], params.potential_strength);
        Vec2 unit_theta = Vec2(std::cos(state.theta[i]), std::sin(state.theta[i]));

        pi += (unit_theta * params.v - grad_V * params.mobility) * params.dt;
    }

    // collision detection
    for (int i = 0; i < state.N; ++i) {
        for (int j = i + 1; j < state.N; ++j) {
            Vec2View pi = state.pos(i);
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
        }
    }
    // boundary checking
    for (int i = 0; i < state.N; ++i) {
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

    ++sim.step_index;
};


// Template method for writing to HDF5
template<typename T>
struct H5TypeMap;

template<>
struct H5TypeMap<double> {
    static H5::PredType type() { return H5::PredType::NATIVE_DOUBLE; }
};

template<>
struct H5TypeMap<unsigned int> {
    static H5::PredType type() { return H5::PredType::NATIVE_UINT; }
};

template<>
struct H5TypeMap<int> {
    static H5::PredType type() { return H5::PredType::NATIVE_INT; }
};

template<typename T, typename H5Container>
void write_H5(const std::string& name, H5Container& write_to, const T& val) {
    hsize_t dims[1] = {1};
    H5::DataSpace ds(1, dims);

    auto dtype = H5TypeMap<T>::type();

    H5::DataSet dset = write_to.createDataSet(name, dtype, ds);
    dset.write(&val, dtype);
}


// main function to run timestepping and data saving
int main() {
    // define trial simulation
    int numParticles = 50;

    double propulsion_speed = 1;
    double particle_radius = 0.02;
    double diffusion = 5;
    double mobility = 0.05;
    double dt = 0.005;
    double box_length = 1;
    double potential_strength = 1;
    unsigned int seed = 10;
    SimParams simulationParameters(propulsion_speed, particle_radius, diffusion, mobility, dt, box_length, potential_strength, seed);

    Simulation simulation(numParticles, simulationParameters);

    // set up HDF5 file
    int numSteps = 5000;
    int saveEvery = 10;
    int T = numSteps / saveEvery;

    H5::H5File file("../../data/particles.h5", H5F_ACC_TRUNC);

    H5::Group g_params = file.createGroup("/params");
    H5::Group g_state  = file.createGroup("/state");
    H5::Group g_meta   = file.createGroup("/meta");

    // write system parameters
    write_H5("v", g_params, simulationParameters.v);
    write_H5("radius", g_params, simulationParameters.radius);
    write_H5("diffusion", g_params, simulationParameters.diffusion);
    write_H5("mobility", g_params, simulationParameters.mobility);
    write_H5("dt", g_params, simulationParameters.dt);
    write_H5("box_length", g_params, simulationParameters.box_length);
    write_H5("potential_strength", g_params, simulationParameters.potential_strength);
    write_H5("seed", g_params, simulationParameters.seed);


    write_H5("numParticles", g_meta, numParticles);
    write_H5("numSteps", g_meta, numSteps);
    write_H5("saveEvery", g_meta, saveEvery);

    // datasets for state parameters
    hsize_t dims[2] = {(hsize_t)T, (hsize_t)numParticles};
    H5::DataSpace traj_space(2, dims);

    H5::DataSet dset_x = g_state.createDataSet("x", H5::PredType::NATIVE_DOUBLE, traj_space);
    H5::DataSet dset_y = g_state.createDataSet("y", H5::PredType::NATIVE_DOUBLE, traj_space);
    H5::DataSet dset_theta = g_state.createDataSet("theta", H5::PredType::NATIVE_DOUBLE, traj_space);

    // memory dataspace for one frame
    hsize_t mem_dims[2] = {1, (hsize_t)numParticles};
    H5::DataSpace mem_space(2, mem_dims);


    // run for however many timesteps and save to HDF5
    int save_index = 0;

    for (int i = 0; i < numSteps; ++i) {
        if (i % saveEvery == 0){
            // hyperslab selection
            hsize_t start[2] = {(hsize_t)save_index, 0};
            hsize_t count[2] = {1, (hsize_t)numParticles};

            // define dataspaces
            H5::DataSpace fs_x = dset_x.getSpace();
            H5::DataSpace fs_y = dset_y.getSpace();
            H5::DataSpace fs_t = dset_theta.getSpace();

            fs_x.selectHyperslab(H5S_SELECT_SET, count, start);
            fs_y.selectHyperslab(H5S_SELECT_SET, count, start);
            fs_t.selectHyperslab(H5S_SELECT_SET, count, start);

            // write current positions to HDF5
            dset_x.write(simulation.state.x.data(), H5::PredType::NATIVE_DOUBLE, mem_space, fs_x);
            dset_y.write(simulation.state.y.data(), H5::PredType::NATIVE_DOUBLE, mem_space, fs_y);
            dset_theta.write(simulation.state.theta.data(), H5::PredType::NATIVE_DOUBLE, mem_space, fs_t);

            save_index++;
        }
        step(simulation);
    }

}