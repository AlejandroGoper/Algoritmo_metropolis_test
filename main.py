# Import modules
from Energy_module import pair_energy, total_energy
from Monte_Carlo_module import displacement_move, volume_move, transfer_move
from Box_module import init_positions, plot_boxes

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
N1 = Ntot // 2 # Number of particles in box 1
N2 = Ntot - N1 # Number of particles in box 2
V1 = Vtot / 2 # Volume of box 1 
V2 = Vtot - V1 # Volume of box 2

# Number of Monte Carlo cycles
n_cycles = 1000

# =============================================================================
# Gibbs Ensemble: Monte Carlo Loop
# =============================================================================

def run_gibbs_ensemble():
    # Initialize positions and box sizes
    positions_box_1, box_1_length = init_positions(N1,V1)
    positions_box_2, box_2_length = init_positions(N2,V2)
    # For the moment ignore this line
    #plot_boxes(positions_box_1, box_1_length, positions_box_2, box_2_length) 

    # Calculating energy
    total_energy_box_1 = total_energy(positions_box_1, box_1_length) 
    total_energy_box_2 = total_energy(positions_box_2, box_2_length)

    count = 0

    # Storage for averages
    densities = []

    for cycle in range(n_cycles):
        result = displacement_move(positions_box_1,box_1_length)
        if result == True:
            count +=1
        # 1) Displacement moves
        
        # 2) Volume exchange move
        
        # 3) Particle transfer moves
        
        # Record densities
        #return densities
    
    print(count)


# ===============================================================================
# End of script: To execute `run_gibbs_ensemble()` 
# ===============================================================================

run_gibbs_ensemble()