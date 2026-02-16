#pragma once
#include "h5sim.h"
#include "h5write.h"

HDF5Handles init_hdf5(
    const std::string& path,
    const SimParams& params,
    int numParticles,
    int numSteps,
    int saveEvery
) {
    HDF5Handles h5 {
        H5::H5File(path, H5F_ACC_TRUNC),
        {}, {}, {},  // 3 empty groups (params, state, metadata)
        {}, {}, {},  // 3 empty datasets (x, y, theta)
        {}           // 1 empty dataspace (memory buffer)
    };

    // groups
    h5.g_params = h5.file.createGroup("/params");
    h5.g_state  = h5.file.createGroup("/state");
    h5.g_meta   = h5.file.createGroup("/meta");

    // parameters
    write_H5("v", h5.g_params, params.v);
    write_H5("radius", h5.g_params, params.radius);
    write_H5("diffusion", h5.g_params, params.diffusion);
    write_H5("mobility", h5.g_params, params.mobility);
    write_H5("dt", h5.g_params, params.dt);
    write_H5("box_length", h5.g_params, params.box_length);
    write_H5("potential_strength", h5.g_params, params.potential_strength);
    write_H5("seed", h5.g_params, params.seed);

    // metadata
    write_H5("numParticles", h5.g_meta, numParticles);
    write_H5("numSteps", h5.g_meta, numSteps);
    write_H5("saveEvery", h5.g_meta, saveEvery);

    int T = numSteps / saveEvery;

    // datasets
    hsize_t dims[2] = {(hsize_t)T, (hsize_t)numParticles};
    H5::DataSpace traj_space(2, dims);

    h5.dset_x = h5.g_state.createDataSet("x", H5::PredType::NATIVE_DOUBLE, traj_space);
    h5.dset_y = h5.g_state.createDataSet("y", H5::PredType::NATIVE_DOUBLE, traj_space);
    h5.dset_theta = h5.g_state.createDataSet("theta", H5::PredType::NATIVE_DOUBLE, traj_space);

    // memory dataspace
    hsize_t mem_dims[2] = {1, (hsize_t)numParticles};
    h5.mem_space = H5::DataSpace(2, mem_dims);

    return h5;
}

void write_frame(HDF5Handles& h5, const Simulation& sim, int save_index) {
    hsize_t start[2] = {(hsize_t)save_index, 0};
    hsize_t count[2] = {1, (hsize_t)sim.state.N};

    // get current dataset spaces
    H5::DataSpace fs_x = h5.dset_x.getSpace();
    H5::DataSpace fs_y = h5.dset_y.getSpace();
    H5::DataSpace fs_t = h5.dset_theta.getSpace();

    // select hyperslabs
    fs_x.selectHyperslab(H5S_SELECT_SET, count, start);
    fs_y.selectHyperslab(H5S_SELECT_SET, count, start);
    fs_t.selectHyperslab(H5S_SELECT_SET, count, start);

    // write current state
    h5.dset_x.write(sim.state.x.data(), H5::PredType::NATIVE_DOUBLE, h5.mem_space, fs_x);
    h5.dset_y.write(sim.state.y.data(), H5::PredType::NATIVE_DOUBLE, h5.mem_space, fs_y);
    h5.dset_theta.write(sim.state.theta.data(), H5::PredType::NATIVE_DOUBLE, h5.mem_space, fs_t);
}
