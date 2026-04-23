// core.cpp
//
// v3. active brownian motion on a curved manifold, but better
/*

varies from previous version by moving from a bounded 2D Euclidean
plane to a closed manifold (sphere or torus). 

now includes git command line for experiment parameters and can be run by an external
bash script to perform parameter sweeps

*/

// my headers
#include "vec.h"       // vector math
#include "state.h"      // particle state, simulation parameters, and sim itself 
#include "dynamics.h"   // potential function and simulation timesteps
#include "h5sim.h"      // defining structure of and writing to hdf5 
#include "manifold.h"   // metric and related functions on manifolds

// std includes
#include <ctime>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <iostream>
#include <filesystem>

std::string getCurrentDateTime() {
    std::time_t now = std::time(nullptr);
    std::tm* localTime = std::localtime(&now);

    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d_%H-%M-%S", localTime);
    return std::string(buffer);
}

// Parse --key value pairs from argv
std::unordered_map<std::string, std::string> parseArgs(int argc, char* argv[]) {
    std::unordered_map<std::string, std::string> args;
    for (int i = 1; i < argc - 1; i += 2) {
        std::string key = argv[i];
        if (key.substr(0, 2) != "--") throw std::invalid_argument("Expected --key, got: " + key);
        args[key.substr(2)] = argv[i + 1];
    }
    return args;
}

int main(int argc, char* argv[]) {
    auto args = parseArgs(argc, argv);

    // helper lambdas for getting args with defaults
    auto getDouble = [&](const std::string& key, double def) {
        return args.count(key) ? std::stod(args[key]) : def;
    };
    auto getUInt = [&](const std::string& key, unsigned int def) {
        return args.count(key) ? (unsigned int)std::stoul(args[key]) : def;
    };
    auto getString = [&](const std::string& key, const std::string& def) {
        return args.count(key) ? args[key] : def;
    };

    // simulation parameters
    int numParticles         = args.count("N")        ? std::stoi(args["N"])       : 50;
    double propulsion_speed  = getDouble("v0",         1.0);
    double particle_radius   = getDouble("radius",     0.02);
    double diffusion         = getDouble("diffusion",  0.05);
    double mobility          = getDouble("mobility",   5.0);
    double dt                = getDouble("dt",         0.005);
    double potential_strength= getDouble("potential",  1.0);
    std::string manifold_type= getString("manifold",   "torus");
    unsigned int seed        = getUInt  ("seed",       10);
    int numSteps             = args.count("steps")    ? std::stoi(args["steps"])   : 10000;
    int saveEvery            = args.count("saveEvery") ? std::stoi(args["saveEvery"]): 1;
    std::string exp_dir      = getString("expdir",  "../../data/trash");
    std::string output       = getString("output", "particles.h5");

    SimParams simulationParameters(propulsion_speed, particle_radius, diffusion, mobility, dt, potential_strength, manifold_type, seed);

    std::unique_ptr<Manifold> manifold;
    if (manifold_type == "torus") { manifold = std::make_unique<TorusManifold>(1.0, 0.3183); } // R = 1, r = 1/pi
    else if (manifold_type == "sphere") { manifold = std::make_unique<SphereManifold>(1.0); } // R = 1
    else if (manifold_type == "euclidean") { manifold = std::make_unique<EuclideanManifold>(); }
    else { std::cerr << "Unknown manifold: " << manifold_type << "\n"; return -1; }

    Simulation simulation(numParticles, simulationParameters, std::move(manifold));

    // set up HDF5 file
    std::string filename = output;

    auto h5 = init_hdf5(exp_dir+"/"+manifold_type+"/"+filename, simulationParameters, numParticles, numSteps, saveEvery);

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