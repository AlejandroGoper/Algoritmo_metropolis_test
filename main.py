# Import modules
from Energy_module import pair_energy, total_energy
from Monte_Carlo_module import displacement_move, volume_move, transfer_move
from Box_module import init_positions

# Import python modules
import random
import math

# Simulation control parameters
T = 1.1                 # Reduced temperature
beta = 1.0 / T           
Ntot = 500              # Total number of particles
rho_tot = 0.6           # Total reduced density
Vtot = Ntot / rho_tot   # Total reduced volume

# Box parameters (Split into two boxes)
N1 = Ntot // 2
N2 = Ntot - N1
V1 = Vtot / 2
V2 = Vtot - V1

# Move step sizes (to be tuned for ~50% acceptance)
dr_max = 0.2            # Max displacement
dlnV_max = 0.01         # Max log-volume change
swap_attempts = int(0.02 * Ntot)  # ~2% of particles per cycle

# Number of Monte Carlo cycles
n_cycles = 50000


# =============================================================================
# Gibbs Ensemble: Monte Carlo Loop
# =============================================================================

def run_gibbs_ensemble():
    # Initialize positions and box sizes
    # Work in progress...
    
    # Storage for averages
    densities = []

    for cycle in range(n_cycles):
        # 1) Displacement moves
        
        # 2) Volume exchange move
        
        # 3) Particle transfer moves
        
        # Record densities
        return densities


# ===============================================================================
# End of script: To execute `run_gibbs_ensemble()` 
# ===============================================================================

run_gibbs_ensemble()