// state.cpp
#include <state.h>

// particle state constructor
ParticleState::ParticleState(int N_): N(N_), q1(N_), q2(N_), theta(N_) {}

// memory efficient way to do vector math on particle coords
Vec2View ParticleState::pos(int i) { return Vec2View{ q1[i], q2[i] };}

// params constructor
SimParams::SimParams(double v_, double r_, double diff_, double mobil_, double dt_, double potenstr_, std::string manifold_type_, unsigned int seed_)
        : v(v_), radius(r_), diffusion(diff_), mobility(mobil_), dt(dt_), potential_strength(potenstr_), manifold_type(manifold_type_), seed(seed_) {}

// simulation constructor
Simulation::Simulation(int N, const SimParams& p, std::unique_ptr<Manifold> m)
        : state(N), params(p), rng(p.seed), manifold(std::move(m)) {
    
    std::uniform_real_distribution<double> rand_angle(0.0, 2*PI);

    // randomly distribute
    if (p.manifold_type == "torus") {   
        for (int i = 0; i < N; ++i) {
                state.q1[i]    = rand_angle(rng);   // phi in (0,2pi)
                state.q2[i]    = rand_angle(rng);   // psi in (0,2pi)
                state.theta[i] = rand_angle(rng);   // tangent space theta in (0,2pi)
        }
    }
    else if (p.manifold_type == "sphere") {
        for (int i = 0; i < N; ++i) {
                state.q1[i]    = rand_angle(rng);      // phi in (0,2pi)
                state.q2[i]    = rand_angle(rng) / 2;  // theta in (0,pi)
                state.theta[i] = rand_angle(rng);      // tangent space theta in (0,2pi)
        }
    }
    else if (p.manifold_type == "euclidean") {
        double L = 2*sqrt(PI); 
        std::uniform_real_distribution<double> pos(-L/2, L/2);  // starting area of 4pi, side lengths of 2(pi)^1/2
        for (int i = 0; i < N; ++i) {
                state.q1[i]    = pos(rng);
                state.q2[i]    = pos(rng);
                state.theta[i] = rand_angle(rng);
        }
    }
}