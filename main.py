# Import modules
from MonteCarloSimulation import MonteCarloSimulation
from Box_module import init_positions, plot_boxes, logic_adjustment,reverse_logic_adjustment

# Import python modules
from numpy.random import rand
import math

# Simulation control parameters  
T = 1.1
beta = 1/T
Ntot = 100              # Total number of particles
rho_tot = 0.3           # Total reduced density
Vtot = Ntot / rho_tot   # Total reduced volume
rcut = 2.5              # Cut-off distance for tailing correction
LJ_epsilon = 1.0        # Parameter of Lennard-Jones potential
LJ_sigma = 1.0          # Parameter of Lennard-Jones potential
dr_max = 0.084          # Maximum step in random displacemt 
dV_max = 1*Vtot     # Maximum step in the ratio of volumes using the parametrization v_new1/v_new2 = exp(dlV_max) 
swap_attempts = int(0.02*Vtot) # Article reccomends this value
volumen_attemps_perswap = 0.1 #Number of volumen chance per swap attemps 
cycle_averaging = 50 #every "cycle_averaging" cylces we take the avarage and save the information
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

    # Melting each box by applying sucessive displacement moves
    for _ in range(10000):
        simulation.displacement_move(positions_box_1, box_1_length)
        simulation.displacement_move(positions_box_2,box_2_length)

    # Calculating energy
    total_energy_box_1 = simulation.total_energy(positions_box_1, box_1_length) 
    total_energy_box_2 = simulation.total_energy(positions_box_2, box_2_length)

    # Storage variables
    densities = []
    v_fraction = []
    n_fraction = [] 
    avg_densities = []
    avg_v_fraction = []
    avg_n_fraction = []
    # Main cycle
    for cycle in range(n_cycles):
        
        # 1) Displacement moves
        for _ in range(N1):
            accepted_dmb1 = simulation.displacement_move(positions=positions_box_1, box_length=box_1_length)
        for _ in range(N2):
            accepted_dmb2 = simulation.displacement_move(positions=positions_box_2, box_length=box_2_length)
        # 2) Volume exchange move
        
        if rand() < volumen_attemps_perswap*swap_attempts:
          x, cajas, npart, vmax = logic_adjustment(positions_box_1, positions_box_2, box_1_length, box_2_length)
          cajas[0][0],cajas[1][0], accepted = simulation.volume_move(x,cajas,npart, beta, dv_max)
          positions_box_1, positions_box_2, box_1_length, box_2_length = reverse_logic_adjustment(x, cajas) #Actualizamos dimensiones de cajas
        # 3) Particle transfer moves
    
        for _ in range(swap_attempts):
            pos1, pos2, L1, L2, accepted = simulation.transfer_move(pos1=positions_box_1, pos2=positions_box_2, L1= box_1_length, L2=box_2_length)

        # Record densities
        densities.append((len(pos1)/L1**3, len(pos2)/L2**3))
        v_fraction.append(((L1**3)/Vtot, (L2**3)/Vtot))
        n_fraction.append((len(pos1)/Ntot, len(pos2)/Ntot))
        if cycle % cycle_averaging == 0: #We take the mean, we dont need that much amount of data
          avg_densities.append(tuple(sum(x)/len(x) for x in zip(*densities)))
          avg_v_fraction.append(tuple(sum(x)/len(x) for x in zip(*v_fraction)))
          avg_n_fraction.append(tuple(sum(x)/len(x) for x in zip(*n_fraction)))
          #we empty the tuples for memory saving
          densities =[]
          v_fraction = []
          n_fraction=[]
        return avg_densities,avg_v_fraction, avg_n_fraction,T
    #print(densities)



# ===============================================================================
# End of script: To execute `run_gibbs_ensemble()` 
# ===============================================================================

run_gibbs_ensemble()
