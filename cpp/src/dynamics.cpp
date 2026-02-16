#include "dynamics.h"

// Heuristic gradient of potential function from V = A(x^4 - x^2 + y^4 - y^2)
Vec2 gradV(double x, double y, double A) {
    // grad(V) = A[4x^3 - 2x, 4y^3 - 2y]
    return Vec2(A*x * (4*x*x - 2), A*y * (4*y*y - 2));
}

// Step particle in time according to 
void step(Simulation& sim) {
    // set random seed and define probability distribution for brownian motion
    auto& state = sim.state;
    auto& params = sim.params;
    auto& rng = sim.rng;
    
    std::normal_distribution<double> noise(0.0, 1.0);

    // update particle trajectories
    for (int i = 0; i < state.N; ++i) {
        // rotational diffusion
        state.theta[i] += std::sqrt(2.0 * params.diffusion * params.dt) * noise(rng);

        // update positions with external field along heading direction
        Vec2View pi = state.pos(i);
        Vec2 grad_V = gradV(state.x[i], state.y[i], params.potential_strength);
        Vec2 unit_theta = Vec2(std::cos(state.theta[i]), std::sin(state.theta[i]));

        pi += (unit_theta * params.v - grad_V * params.mobility) * params.dt;
    }

    // collision detection
    for (int i = 0; i < state.N; ++i) {
        for (int j = i + 1; j < state.N; ++j) {
            Vec2View pi = state.pos(i);
            Vec2View pj = state.pos(j);

            Vec2 separation_vec = Vec2(pj) - Vec2(pi); 
            auto separation = separation_vec.norm();

            // if particle centers are closer than twice the radius, they are intersecting, so push them apart
            if (separation < 2*params.radius && separation > 1e-10) {
                // flip particle's heading directions
                state.theta[i] = -state.theta[i];
                state.theta[j] = -state.theta[j];

                // move particles apart 
                auto n = separation_vec * (1 / separation);  // unit vector along separation

                auto dist_to_move = 2*params.radius - separation; // particles move apart so their boundaries don't intersect
                Vec2 disp = n * dist_to_move * 0.5 ; // displacement vector for each is along n, half of total displacment

                pi -= disp;
                pj += disp;
            }
        }
    }
    // boundary checking
    for (int i = 0; i < state.N; ++i) {
        if (state.x[i] - params.radius < -params.box_length) {
            state.x[i] = -params.box_length + params.radius;
            state.theta[i] = PI - state.theta[i];       // (vx, vy) -> (-vx, vy) 
        }
        if (state.x[i] + params.radius > params.box_length) {
            state.x[i] = params.box_length - params.radius;
            state.theta[i] = PI - state.theta[i] ;
        }
        if (state.y[i] - params.radius < -params.box_length) {
            state.y[i] = -params.box_length + params.radius;
            state.theta[i] = -state.theta[i];             // (vx, vy) -> (vx, -vy) 
        }
        if (state.y[i] + params.radius > params.box_length) {
            state.y[i] = params.box_length - params.radius;
            state.theta[i] = -state.theta[i];
        }
    }

    ++sim.step_index;
};
