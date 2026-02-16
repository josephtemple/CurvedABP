#pragma once
#include <vector>
#include <random>
#include "vec2.h"

// Data structure for the state of all the particles
struct ParticleState {
    int N;         
    std::vector<double> x; 
    std::vector<double> y;
    std::vector<double> theta; 

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