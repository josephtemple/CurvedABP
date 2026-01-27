# abp_core.pyi
from typing import Any, Tuple

# -------------------------
# Particle state
# -------------------------
class ParticleState:
    x: Any      # numpy array of shape (N,)
    y: Any      # numpy array of shape (N,)
    theta: Any  # numpy array of shape (N,)

    def __init__(self) -> None: ...
    
# -------------------------
# Simulation parameters
# -------------------------
class SimParams:
    v: float
    radius: float
    diffusion: float
    mobility: float
    dt: float
    box_length: float
    potential_strength: float
    seed: int

    def __init__(self) -> None: ...

# -------------------------
# Step function
# -------------------------
def step(state: ParticleState, params: SimParams) -> None: ...
