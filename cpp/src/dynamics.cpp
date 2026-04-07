#include "dynamics.h"
#include "manifold.h"

// Heuristic gradient of potential function from V = A(x^4 - x^2 + y^4 - y^2)
Vec2 gradV(double x, double y, double A) {
    // grad(V) = A[4x^3 - 2x, 4y^3 - 2y]
    return Vec2(0,0);
}

// Step particle in time according to 
void step(Simulation& sim) {
    // set random seed and define probability distribution for brownian motion
    auto& state = sim.state;
    auto& params = sim.params;
    auto& rng = sim.rng;
    
    std::normal_distribution<double> noise(0.0, 1.0);

    // unpack manifold object
    auto& manifold = *sim.manifold;

    // update particle trajectories
    for (int i = 0; i < state.N; ++i) {
        // geometry things. vielbein translates between manifolds coord basis and orthonormal basis of tangent space
        Vec2View p_i = state.pos(i);
        double e1 = manifold.vielbein1(state.q1[i], state.q2[i]);
        double e2 = manifold.vielbein2(state.q1[i], state.q2[i]);

        // potential
        Vec2 dV = gradV(state.q1[i], state.q2[i], params.potential_strength);
        Vec2 grad_V = Vec2(e1*dV.x, e2*dV.y);

        // theta vector
        Vec2 unit_theta = Vec2(e1 * std::cos(state.theta[i]), e2 * std::sin(state.theta[i]));

        // displacement of coords
        Vec2 disp = (unit_theta * params.v - grad_V * params.mobility) * params.dt;

        // theta update with spin connection to account for tangent space rotating across manifold
        state.theta[i] += manifold.connection(state.q1[i], state.q2[i], disp.x, disp.y)
                          + std::sqrt(2.0 * params.diffusion * params.dt) * noise(rng);
        p_i += disp;
    }

    // collision detection
    for (int i = 0; i < state.N; ++i) {
        for (int j = i + 1; j < state.N; ++j) {
            Vec2View pi = state.pos(i);
            Vec2View pj = state.pos(j);

            Vec2 dq = Vec2(pj) - Vec2(pi);

            // metric-weighted geodesic distance (valid for small separations)
            double g1 = manifold.g11(state.q1[i], state.q2[i]);
            double g2 = manifold.g22(state.q1[i], state.q2[i]);
            double separation = std::sqrt(g1*dq.x*dq.x + g2*dq.y*dq.y);

            if (separation < 2*params.radius && separation > 1e-10) {
                // flip headings
                state.theta[i] = -state.theta[i];
                state.theta[j] = -state.theta[j];

                // push apart along geodesic
                Vec2 n = Vec2(std::sqrt(g1) * dq.x, std::sqrt(g2) * dq.y);
                n = n * (1.0 / n.norm());

                double dist_to_move = 2*params.radius - separation;
                Vec2 disp = n * dist_to_move * 0.5;

                pi -= disp;
                pj += disp;
            }
        }
    }
    
    ++sim.step_index;
};
