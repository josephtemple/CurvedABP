"""
animate.py

animate the motion of active brownian particles on curved manifolds under the effects of orientational alignment and rotational diffusion

ensure the cwd is the python dir or else the interpreter will get confused poor buddy

"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from dataclasses import dataclass

import os
from pathlib import Path

# 1. Get the path to the current file
script_path = Path(__file__).resolve()

# 2. Get the directory of that file
script_dir = script_path.parent

# 3. Change the working directory
os.chdir(script_dir)

### LOAD THAT DATA TWIN ###
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
    
data_dir = script_path.parent.parent / 'data' / '2026-04-27_23-57-25'
metrics = ['sphere', 'torus', 'euclidean']
ds, cs, ss = [0.0, 0.2, 1.0, 5.0], [0.0, 0.2, 1.0, 5.0], [10, 20, 30]

metric = metrics[2]
d = ds[0]
c = cs[3]
s = ss[2]

dataset = f'd{d}_c{c}_s{s}.h5'

file_path = data_dir / metric / dataset

data = load(file_path)


# Convert to Cartesian on sphere and torus, and set some plotting params
phi, theta = data.q1, data.q2
R, r = 0,0
if (metric == metrics[0]) :
    R = 1
    x = R*np.sin(theta) * np.cos(phi)
    y = R*np.sin(theta) * np.sin(phi)
    z = R*np.cos(theta)
elif (metric == metrics[1]) :
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

# initial frame
scat = ax.scatter(x[0], y[0], z[0], c='blue', s=10, alpha=0.8)

# set limits based on geometry
limit = 5
if (metric == metrics[0]):
    limit = R + 0.1
elif (metric == metrics[1]):
    limit = R + r + 0.1

ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.set_zlim(-limit, limit)

ax.set_axis_off()

# add wireframe surface
if (metric == metrics[0]):
    u = np.linspace(0, 2*np.pi, 30)
    v = np.linspace(0, np.pi, 30)
    U, V = np.meshgrid(u, v)

    X_surf = R * np.cos(U) * np.sin(V)
    Y_surf = R * np.sin(U) * np.sin(V)
    Z_surf = R * np.cos(V)

elif (metric == metrics[1]):
    u = np.linspace(0, 2*np.pi, 30)
    v = np.linspace(0, 2*np.pi, 30)
    U, V = np.meshgrid(u, v)

    X_surf = (R + r * np.cos(V)) * np.cos(U)
    Y_surf = (R + r * np.cos(V)) * np.sin(U)
    Z_surf = r * np.sin(V)

else :
    u = np.linspace(-limit, limit, 30)
    v = np.linspace(-limit, limit, 30)
    U, V = np.meshgrid(u, v)

    X_surf = U
    Y_surf = V
    Z_surf = np.zeros_like(U)

ax.plot_wireframe(X_surf, Y_surf, Z_surf, color='gray', alpha=0.1)


# helper function for naming
def anim_title(m: str, d: float, c: float, s: int) -> str:
    if (m == metrics[0]):
        m_str = "Sphere"
    elif (m == metrics[1]):
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

ax.set_title(anim_title(metric, d, c, s))


# build the animation
def update(frame_idx):
    scat._offsets3d = (x[frame_idx], y[frame_idx], z[frame_idx]) # type: ignore
    return scat,

ani = FuncAnimation(
    fig, 
    update, 
    frames=range(0, len(x), 5), 
    interval=30, 
    blit=False
)

plt.show()