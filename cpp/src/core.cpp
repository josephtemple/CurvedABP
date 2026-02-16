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

// my headers
#include "vec2.h"       // vector math
#include "state.h"      // particle state, simulation parameters, and sim itself 
#include "dynamics.h"   // potential function and simulation timesteps
#include "h5sim.h"      // defining structure of and writing to hdf5 

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

    std::string data_dir = "../../data/";
    std::string filename = "particles.h5";
    auto h5 = init_hdf5(data_dir+filename, simulationParameters, numParticles, numSteps, saveEvery);

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