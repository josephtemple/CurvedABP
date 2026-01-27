import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.rcParams['figure.dpi'] = 150

import abp_core

# ---------------------------
# Simulation parameters
# ---------------------------
L = 1
N = 1000
v0 = 2
particle_spawn_lim = 0.7 * L

# initialize positions and orientations
xs, ys = np.random.uniform(-particle_spawn_lim, particle_spawn_lim, (2, N))
thetas = np.random.uniform(0, 2*np.pi, N)

# C++ objects
state = abp_core.ParticleState()
state.x = xs.copy()
state.y = ys.copy()
state.theta = thetas.copy()

params = abp_core.SimParams()
params.v = v0
params.radius = 0.02
params.diffusion = 5.0
params.mobility = 0.5
params.dt = 1/500
params.box_length = L
params.potential_strength = 3.0
params.seed = 42

# ---------------------------
# Plot setup
# ---------------------------
fig, ax = plt.subplots()
ax.set_xlim(-L, L)
ax.set_ylim(-L, L)
ax.set_xticks([])
ax.set_yticks([])

# Potential field
n_pts = 500
x_range = np.linspace(-L, L, n_pts)
y_range = np.linspace(-L, L, n_pts)
X, Y = np.meshgrid(x_range, y_range)
potential = 3*(X**4 - X**2) + 3*(Y**4 - Y**2)  # scaled 4r^3-2r
im = ax.imshow(potential, extent=(-L, L, -L, L), cmap='PiYG')
cbar = fig.colorbar(im, ax=ax, label='Potential')

# Use a single scatter for all particles
fig_width_inch = fig.get_size_inches()[0]
ax_width_data = ax.get_xlim()[1] - ax.get_xlim()[0]
diam_points = params.radius * (fig_width_inch / ax_width_data) * 72
s = np.pi * (diam_points / 2)**2

scatter = ax.scatter(state.x, state.y, s=s, c='white', edgecolors='black')

# ---------------------------
# Animation functions
# ---------------------------
def init():
    scatter.set_offsets(np.c_[state.x, state.y])
    return (scatter,)

def animate(frame):
    abp_core.step(state, params)
    scatter.set_offsets(np.c_[state.x, state.y])
    return (scatter,)

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=300, interval=5, blit=True)
plt.show()
