#include <state.h>

// particle state constructor
ParticleState::ParticleState(int N_): N(N_), q1(N_), q2(N_), theta(N_) {}

// memory efficient way to do vector math on particle coords
Vec2View ParticleState::pos(int i) { return Vec2View{ q1[i], q2[i] };}

// params constructor
SimParams::SimParams(double v_, double r_, double diff_, double mobil_, double dt_, double boq1length_, double potenstr_, unsigned int seed_)
        : v(v_), radius(r_), diffusion(diff_), mobility(mobil_), dt(dt_), box_length(boq1length_), potential_strength(potenstr_), seed(seed_) {}

// simulation constructor
Simulation::Simulation(int N, const SimParams& p)
        : state(N), params(p), rng(p.seed) {}