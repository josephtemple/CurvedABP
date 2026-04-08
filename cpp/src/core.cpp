// core.cpp
//
// v2. active brownian motion on a curved manifold
/*

varies from previous version by moving from a bounded 2D Euclidean
plane to a closed manifold (sphere or torus). 

some things stay the same
- stochastic rotational diffusion (now in the tangent space)
- 

*/

// my headers
#include "vec.h"       // vector math
#include "state.h"      // particle state, simulation parameters, and sim itself 
#include "dynamics.h"   // potential function and simulation timesteps
#include "h5sim.h"      // defining structure of and writing to hdf5 
#include "manifold.h"   // metric and related functions on manifolds

// std includes
#include <ctime>

std::string getCurrentDateTime() {
    std::time_t now = std::time(nullptr);
    std::tm* localTime = std::localtime(&now);

    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d_%H-%M-%S_", localTime);
    return std::string(buffer);
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
    double potential_strength = 1;
    std::string manifold_type = "torus";
    unsigned int seed = 10;
    SimParams simulationParameters(propulsion_speed, particle_radius, diffusion, mobility, dt, potential_strength, manifold_type, seed);

    // this doesn't work rn.
    auto manifold = std::make_unique<SphereManifold>(2.0);
    if (manifold_type == "torus") { auto manifold = std::make_unique<TorusManifold>(2.0, 0.5); } // R = 2, r =0.5
    else if (manifold_type == "sphere") { auto manifold = std::make_unique<SphereManifold>(2.0); } // R
    else if (manifold_type == "euclidean") {auto manifold = std::make_unique<EuclideanManifold>(); }
    else {return -1;} // failure, undefined manifold type

    Simulation simulation(numParticles, simulationParameters, std::move(manifold));

    // set up HDF5 file`
    int numSteps = 5000;
    int saveEvery = 10;

    std::string data_dir = "../../data/";
    std::string time = getCurrentDateTime();
    std::string filename = "particles.h5";
    auto h5 = init_hdf5(data_dir+time+filename, simulationParameters, numParticles, numSteps, saveEvery);

    // run for however many timesteps and save to HDF5
    int save_index = 0;

    for (int i = 0; i < numSteps; ++i) {
        if (i % saveEvery == 0){
            write_frame(h5, simulation, save_index);
            save_index++;
        }
        step(simulation);
    }
}