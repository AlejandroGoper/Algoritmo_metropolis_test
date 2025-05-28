# Import modules
from MonteCarloSimulation import MonteCarloSimulation
from Box_module import init_positions, plot_boxes, logic_adjustment,reverse_logic_adjustment

# Import python modules
from numpy.random import rand
import math

# Simulation control parameters  
#T = 0.9
#beta = 1/T
rcut = 2.5              # Cut-off distance for tailing correction
LJ_epsilon = 1.0        # Parameter of Lennard-Jones potential
LJ_sigma = 1.0          # Parameter of Lennard-Jones potential
swap_attemps = 20
volumen_attemps = 2 
dis_attemps = 200
subcycles = swap_attemps + volumen_attemps + dis_attemps
# Number of Monte Carlo cycles
n_cycles = 4000



# =============================================================================
# Gibbs Ensemble: Monte Carlo Loop
# =============================================================================

def run_gibbs_ensemble(T):
    beta = 1/T
    # Storage variables
    densities = []
    v_fraction = []
    n_fraction = [] 
    avg_densities = []
    avg_v_fraction = []
    avg_n_fraction = []
    #Parameters need to be local
    Ntot = 256             # Total number of particles
    rho_tot = 0.3           # Total reduced density
    Vtot = Ntot / rho_tot   # Total reduced volume
    dr_max = 0.15          # Maximum step in random displacemt 
    dV_max = 0.02*Vtot
    # Box parameters (Split into two boxes)
    N1 = Ntot // 2 # Number of particles in box 1
    N2 = Ntot - N1 # Number of particles in box 2
    V1 = Vtot / 2 # Volume of box 1 
    V2 = Vtot - V1 # Volume of box 2
    acceptedV = 0
    triesV = 0
    acceptedD = 0
    triesD = 0
    acceptedS = 0
    triesS = 0
    simulation = MonteCarloSimulation(particles=Ntot, red_temperature=T,LJ_epsilon=LJ_epsilon, LJ_sigma=LJ_sigma,
                                  LJ_rcut=rcut, dr_max=dr_max, dlnV_max=dV_max)
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

    
    
    # Main cycle
    for cycle in range(n_cycles):
        for sub in range(swap_attemps + volumen_attemps + dis_attemps): 
            # 1) Displacement moves
            if (subcycles*rand() < dis_attemps):
                if (rand() > 0.5):
                    accepted_dmb1 = simulation.displacement_move(positions=positions_box_1, box_length=box_1_length)
                    if accepted_dmb1:
                        acceptedD += 1
                    triesD += 1
                else:
                    accepted_dmb2 = simulation.displacement_move(positions=positions_box_2, box_length=box_2_length)
                    if accepted_dmb2:
                        acceptedD += 1
                    triesD += 1
            # 2) Volume exchange move
        
            if subcycles*rand() < volumen_attemps:
                x, cajas, npart, vmax = logic_adjustment(positions_box_1, positions_box_2, box_1_length, box_2_length)
                cajas[0][0],cajas[1][0], accepted = simulation.volume_move(x,cajas,npart, beta, dV_max)
                if accepted:
                    acceptedV += 1
                triesV += 1
                positions_box_1, positions_box_2, box_1_length, box_2_length = reverse_logic_adjustment(x, cajas) #Actualizamos dimensiones de cajas
            # 3) Particle transfer moves
    
            if subcycles*rand() < swap_attemps:
                positions_box_1, positions_box_2, box_1_length, box_2_length , accepted = simulation.transfer_move(pos1=positions_box_1, pos2=positions_box_2, L1= box_1_length, L2=box_2_length)
                if accepted:
                    acceptedS += 1
                triesS += 1
            # Record densities
            densities.append((len(positions_box_1)/box_1_length**3, len(positions_box_2)/box_2_length**3))
            v_fraction.append(((box_1_length**3)/Vtot, (box_2_length**3)/Vtot))
            n_fraction.append((len(positions_box_1)/Ntot, len(positions_box_2)/Ntot))
        avg_densities.append(tuple(sum(x)/len(x) for x in zip(*densities)))
        avg_v_fraction.append(tuple(sum(x)/len(x) for x in zip(*v_fraction)))
        avg_n_fraction.append(tuple(sum(x)/len(x) for x in zip(*n_fraction)))
        #we empty the tuples for memory saving
        densities =[]
        v_fraction = []
        n_fraction=[]
        if (cycle % 100 == 0):
            print(f"{cycle} completed")
            if(cycle > 0):
                print(f"Se aceptó el {acceptedS*100/triesS}% de los swaps, el {acceptedV*100/triesV} de los volumen y {acceptedD*100/triesD} de los desplazamientos" )
        if (cycle % 10 == 0 & cycle != 0):
            if acceptedV/triesV > 0.55:
                dV_max *= 1.05
            elif  acceptedV/triesV < 0.45:
                dV_max *= 0.95
            acceptedV = 0
            triesV = 0
        if acceptedD/triesD > 0.55:
            dr_max *= 1.05
        elif  acceptedD/triesD < 0.45:
            dr_max *= 0.95
        acceptedD = 0
        triesD = 0
    return avg_densities,avg_v_fraction, avg_n_fraction



# ===============================================================================
# End of script: To execute `run_gibbs_ensemble()` 
# ===============================================================================

#run_gibbs_ensemble()

