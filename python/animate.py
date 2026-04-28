"""
animate.py

animate the motion of active brownian particles on curved manifolds under the effects of orientational alignment and rotational diffusion
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from dataclasses import dataclass

import os
from pathlib import Path

### LOAD THAT DATA TWIN ###
script_path = Path(__file__).resolve()
script_dir = script_path.parent
os.chdir(script_dir)  # set script directory as cwd for relative paths


@dataclass
class SimData:
    q1: np.ndarray      # (n_frames, N) - phi
    q2: np.ndarray      # (n_frames, N) - theta
    theta: np.ndarray   # (n_frames, N) - orientation

def load(path: Path, frames=None) -> SimData:
    with h5py.File(path, 'r') as f:
        sl = frames if frames is not None else slice(None)
        return SimData(
            q1    = np.asarray(f['state/q1'])[sl],
            q2    = np.asarray(f['state/q2'])[sl],
            theta = np.asarray(f['state/theta'])[sl]
        )
    
data_dir = script_path.parent.parent / 'data' / '2026-04-28_16-25-59'
manifolds = ['sphere', 'torus', 'euclidean']
ds, cs, ss = [0.0, 0.2, 1.0, 5.0], [0.0, 0.2, 1.0, 5.0], [10, 20, 30]

manifold = manifolds[1]
d = ds[3]
c = cs[1]
s = ss[2]

dataset = f'd{d}_c{c}_s{s}.h5'

file_path = data_dir / manifold / dataset

data = load(file_path)


# Convert to Cartesian on sphere and torus, and set some plotting params
phi, theta = data.q1, data.q2
R, r = 0,0
if (manifold == manifolds[0]) :
    R = 1
    x = R*np.sin(theta) * np.cos(phi)
    y = R*np.sin(theta) * np.sin(phi)
    z = R*np.cos(theta)
elif (manifold == manifolds[1]) :
    # ignore the bad coding practices (hard coded R = 1, r = 1/pi for torus)
    R, r = 1, 1/np.pi
    x = (R + r*np.cos(theta)) * np.cos(phi)
    y = (R + r*np.cos(theta)) * np.sin(phi)
    z = r*np.sin(theta)
else:
    # this is a not very clear way to do this (better would just be x=q1 etc but phi and theta are already defined)
    x = phi
    y = theta
    z = np.zeros_like(phi)


### MATPLOTLIB THAT JAWN ###
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')


# set limits based on geometry
limit = 2 * np.sqrt(np.pi) - 1
if (manifold == manifolds[0]):
    limit = R + 0.1
elif (manifold == manifolds[1]):
    limit = R + r + 0.1


limit *= 0.8
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.set_zlim(-limit, limit)

ax.set_axis_off()

# add wireframe surface
frame_density = 30
if (manifold == manifolds[0]):
    u = np.linspace(0, 2*np.pi, frame_density)
    v = np.linspace(0, np.pi, frame_density)
    U, V = np.meshgrid(u, v)

    X_surf = R * np.cos(U) * np.sin(V)
    Y_surf = R * np.sin(U) * np.sin(V)
    Z_surf = R * np.cos(V)
elif (manifold == manifolds[1]):
    u = np.linspace(0, 2*np.pi, frame_density)
    v = np.linspace(0, 2*np.pi, frame_density)
    U, V = np.meshgrid(u, v)

    X_surf = (R + r * np.cos(V)) * np.cos(U)
    Y_surf = (R + r * np.cos(V)) * np.sin(U)
    Z_surf = r * np.sin(V)
else :
    u = np.linspace(-limit, limit, frame_density)
    v = np.linspace(-limit, limit, frame_density)
    U, V = np.meshgrid(u, v)

    X_surf = U
    Y_surf = V
    Z_surf = np.zeros_like(U)

ax.plot_wireframe(X_surf, Y_surf, Z_surf, color='gray', alpha=0.12)


# helper function for naming
def anim_title(m: str, d: float, c: float, s: int) -> str:
    if (m == manifolds[0]):
        m_str = "Sphere"
    elif (m == manifolds[1]):
        m_str = "Torus"
    else:
        m_str = "2D Euclidean Space"


    dc_map = {
        0.0: "No",
        0.2: "Weak",
        1.0: "Intermediate",
        5.0: "Strong"
    }
    
    return f"Active Brownian Particles\n{dc_map[d]} diffusion (d = {d}) and {dc_map[c]} coupling (c = {c}) on a {m_str}. \nRandom seed = {s}"

ax.set_title(anim_title(manifold, d, c, s))


# different colors for different manifolds bc why not
if (manifold == 'sphere'):
    color = 'green'
    color_channel_idx = 1
elif (manifold == 'torus'):
    color = 'blue'
    color_channel_idx = 2
else:
    color = 'red'
    color_channel_idx = 0

# initial frame
scat = ax.scatter(x[0], y[0], z[0], c=color)
time_text = ax.text2D(0.95, 0.05, "Step: 0", transform=ax.transAxes, 
                      fontsize=12, color='black')

# update for all frames after that
def update(frame_idx):
    xi, yi, zi = x[frame_idx], y[frame_idx], z[frame_idx]
    
    scat._offsets3d = (xi, yi, zi)  # type: ignore

    # transform points to camera coordinates to get depth
    # M is the 4x4 projection matrix that updates when you rotate the view
    M = ax.get_proj()
    points = np.vstack([xi, yi, zi, np.ones_like(xi)])
    transformed = M @ points
    
    # use the transformed z-coordinate for depth
    depths = transformed[2]
    d_norm = (depths - depths.min()) / (depths.max() - depths.min() + 1e-10)

    # make particles on the back be lighter and smaller
    alphas = 0.8 + 0.2 * d_norm    # higher alpha = more opaque at front
    sizes = 20 + 40 * d_norm        # bigger at front

    scat.set_sizes(sizes)

    # adjust brightness wrt alphas
    colors = np.zeros((len(xi), 4))
    colors[:, color_channel_idx] = 0.5             # color channel
    colors[:, 3] = alphas          # alpha channel
    scat.set_facecolors(colors)    # type: ignore
    scat.set_edgecolors(colors)    # type: ignore


    # update the timestep counter
    time_text.set_text(f"Step: {frame_idx}")
    
    return scat, time_text

ani = FuncAnimation(
    fig, 
    update, 
    frames=range(0, len(x), 2), 
    interval=30, 
    blit=False
)

plt.show()