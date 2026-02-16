#pragma once
#include <state.h>

// particle state constructor
ParticleState::ParticleState(int N_): N(N_), x(N_), y(N_), theta(N_) {}

// memory efficient way to do vector math on particle coords
Vec2View ParticleState::pos(int i) { return Vec2View{ x[i], y[i] };}

// params constructor
SimParams::SimParams(double v_, double r_, double diff_, double mobil_, double dt_, double boxlength_, double potenstr_, unsigned int seed_)
        : v(v_), radius(r_), diffusion(diff_), mobility(mobil_), dt(dt_), box_length(boxlength_), potential_strength(potenstr_), seed(seed_) {}

// simulation constructor
Simulation::Simulation(int N, const SimParams& p)
        : state(N), params(p), rng(p.seed) {}