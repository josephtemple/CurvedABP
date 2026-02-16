#pragma once
#include <string>
#include <H5Cpp.h>
#include "state.h"

struct HDF5Handles {
    H5::H5File file;
    H5::Group g_params;
    H5::Group g_state;
    H5::Group g_meta;

    H5::DataSet dset_x;
    H5::DataSet dset_y;
    H5::DataSet dset_theta;

    H5::DataSpace mem_space;
};

HDF5Handles init_hdf5(
    const std::string& path,
    const SimParams& params,
    int numParticles,
    int numSteps,
    int saveEvery
);

void write_frame(HDF5Handles& h5, const Simulation& sim, int save_index);