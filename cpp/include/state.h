#pragma once
#include <vector>
#include <random>
#include <memory>
#include "vec.h"
#include "manifold.h"

// Data structure for the state of all the particles
struct ParticleState {
    int N;         
    std::vector<double> q1; 
    std::vector<double> q2;
    std::vector<double> theta; 

    ParticleState(int N_);

    // make available a way to do vector math on the coordinates, but keep SoA memory layout
    Vec2View pos(int);
};

// Data structure for the constant parameters of the simulation
struct SimParams {
    double v;                         // velocity of each particle            
    double radius;                    // radius of each particle in params.box_length
    double diffusion;                 // rotational diffusion constant in rad^2 / time
    double coupling_str;              // 
    double interaction_rad;
    double dt;                        // time step size
    std::string manifold_type;        // type of manifold. allowed: {sphere, torus}
    unsigned int seed;                // random number generation seed

    SimParams(double v_, double r_, double diff_, double coup_, double inter_, double dt_, std::string manifold_type_, unsigned int seed_);
};

// Data structure for packaging the whole simulation
struct Simulation {
    ParticleState state;
    SimParams params;
    std::unique_ptr<Manifold> manifold;
    std::mt19937 rng;
    std::size_t step_index = 0;

    Simulation(int N, const SimParams& p, std::unique_ptr<Manifold> m);
};