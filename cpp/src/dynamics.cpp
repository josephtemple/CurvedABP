#include "dynamics.h"
#include "manifold.h"

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

        // theta vector
        Vec2 unit_theta = Vec2(e1 * std::cos(state.theta[i]), e2 * std::sin(state.theta[i]));

        // displacement of coords
        Vec2 disp = (unit_theta * params.v) * params.dt;

        // theta update with spin connection to account for tangent space rotating across manifold
        state.theta[i] += manifold.connection(state.q1[i], state.q2[i], disp.x, disp.y)   // spin connection
                        + std::sqrt(2.0 * params.diffusion * params.dt) * noise(rng);   // rotational diffusion
        p_i += disp;
    }

    // pre-allocate vicsek accumulators, seeded with own orientation
    std::vector<double> sin_sum(state.N), cos_sum(state.N);
    for (int i = 0; i < state.N; ++i) {
        sin_sum[i] = std::sin(state.theta[i]);
        cos_sum[i] = std::cos(state.theta[i]);
    }

    // collision detection & vicsek orientation averaging
    for (int i = 0; i < state.N; ++i) {
        for (int j = i + 1; j < state.N; ++j) {
            Vec2View pi = state.pos(i);
            Vec2View pj = state.pos(j);

            Vec2 dq = Vec2(pj) - Vec2(pi);

            // metric-weighted geodesic distance (valid for small separations)
            double g1 = manifold.g11(state.q1[i], state.q2[i]);
            double g2 = manifold.g22(state.q1[i], state.q2[i]);
            double separation = std::sqrt(g1*dq.x*dq.x + g2*dq.y*dq.y);

            // move apart if too close
            if (separation < 2*params.radius && separation > 1e-10) {
                // flip heading vectors along collision normal
                Vec2 n = Vec2(std::sqrt(g1) * dq.x, std::sqrt(g2) * dq.y);
                n = n * (1.0 / n.norm()); // because i didn't define scalar division lol

                auto reflect = [&](double theta) {
                    Vec2 v = Vec2(std::cos(theta), std::sin(theta));
                    double vn = v.x * n.x + v.y * n.y;
                    Vec2 vr = Vec2(v.x - 2*vn*n.x, v.y - 2*vn*n.y);
                    return std::atan2(vr.y, vr.x);
                };
                state.theta[i] = reflect(state.theta[i]);
                state.theta[j] = reflect(state.theta[j]);

                // push apart along geodesic
                double dist_to_move = 2*params.radius - separation;
                Vec2 disp = dq * (1/separation) * dist_to_move * 0.5;

                pi -= disp;
                pj += disp;
            }

            // vicsek orientational alignment type beat i started this too late bruh
            if (separation < params.interaction_rad) {
                // parallel transport theta into shared tangent space via spin connection
                // (first order approx — good enough for small interaction_rad)
                double transported_theta_j = state.theta[j]
                    - manifold.connection(state.q1[j], state.q2[j], dq.x, dq.y);
                double transported_theta_i = state.theta[i]
                    + manifold.connection(state.q1[i], state.q2[i], dq.x, dq.y);

                // accumulate symmetrically so we only need the i<j pairs
                sin_sum[i] += std::sin(transported_theta_j);
                cos_sum[i] += std::cos(transported_theta_j);
                sin_sum[j] += std::sin(transported_theta_i);
                cos_sum[j] += std::cos(transported_theta_i);
            }
        }
    }

    // apply vicsek alignment
    for (int i = 0; i < state.N; ++i) {
        // compute mean theta from summed sin and cos components
        double theta_mean = std::atan2(sin_sum[i], cos_sum[i]);

        // shift theta towards mean of nearby particles sin(mean - i) means we add nothing if
        // they're close (sin 0 = 0) and max if theyre perpendicular (sin 90 )
        state.theta[i] += params.coupling_str * std::sin(theta_mean - state.theta[i]) * params.dt;
    }
    
    ++sim.step_index;
};