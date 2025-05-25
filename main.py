# Import modules
from MonteCarloSimulation import MonteCarloSimulation
from Box_module import init_positions, plot_boxes, logic_adjustment

# Import python modules
import random
import math

# Simulation control parameters
T = 1.1                 # Reduced temperature
beta = 1.0 / T           
Ntot = 100              # Total number of particles
rho_tot = 0.3           # Total reduced density
Vtot = Ntot / rho_tot   # Total reduced volume
rcut = 2.5              # Cut-off distance for tailing correction
LJ_epsilon = 1.0        # Parameter of Lennard-Jones potential
LJ_sigma = 1.0          # Parameter of Lennard-Jones potential
dr_max = 0.084          # Maximum step in random displacemt 
dlnV_max = 0.01*Vtot        # Maximum step in the ratio of volumes using the parametrization v_new1/v_new2 = exp(dlV_max) 

# Box parameters (Split into two boxes)
N1 = Ntot // 2 # Number of particles in box 1
N2 = Ntot - N1 # Number of particles in box 2
V1 = Vtot / 2 # Volume of box 1 
V2 = Vtot - V1 # Volume of box 2

# Number of Monte Carlo cycles
n_cycles = 100

# =============================================================================
#  Class Initialization
# =============================================================================

simulation = MonteCarloSimulation(particles=Ntot, red_temperature=T,LJ_epsilon=LJ_epsilon, LJ_sigma=LJ_sigma,
                                  LJ_rcut=rcut, dr_max=dr_max, dlnV_max=dlnV_max)

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
    total_energy_box_1 = simulation.total_energy(positions_box_1, box_1_length) 
    total_energy_box_2 = simulation.total_energy(positions_box_2, box_2_length)

    # Storage for averages
    densities = []
    print(total_energy_box_1)
    count=0
    for cycle in range(n_cycles):
        
        # 1) Displacement moves
        #accepted = simulation.displacement_move(positions=positions_box_1, box_length=box_1_length)

        # 2) Volume exchange move
        x, cajas, npart, vmax = logic_adjustment(positions_box_1, positions_box_2, box_1_length, box_2_length)
        box1n, box2n, accepted = simulation.volume_move(x,cajas,npart, beta, vmax)
        if(not accepted):
            continue
        else:
            print("Accepted")
            count +=1
        # 3) Particle transfer moves
        

        #pos1, pos2, L1, L2, accepted = transfer_move(pos1=positions_box_1, pos2=positions_box_2, L1= box_1_length, L2=box_2_length)
        #print("Originals", len(positions_box_1), len(positions_box_2), "Accepted:", accepted)
        
        #if(accepted):
        #    print("New ones:", len(pos1), len(pos2))
        #    plot_boxes(pos1, L1, pos2, L2)

        # Record densities
        #return densities
    print("Acceptance:", count/n_cycles)
    print(simulation.total_energy(positions_box_1, box_1_length))



# ===============================================================================
# End of script: To execute `run_gibbs_ensemble()` 
# ===============================================================================

run_gibbs_ensemble()