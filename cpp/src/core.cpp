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

// std includes
#include <cmath>
#include <vector>

// my headers
#include "vec2.h"
#include "state.h"
#include "dynamics.h"
#include "h5write.h"

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