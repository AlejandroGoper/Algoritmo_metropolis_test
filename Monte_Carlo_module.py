#Importamos los modulos de energia
from Energy_module import pair_energy, total_energy
from Box_module import reverse_logic_adjustment
from random import randrange, random
from numpy import round, dot, log
from numpy.random import rand
from math import exp


# =============================================================================
# Monte Carlo Move Functions
# =============================================================================

Ntot = 500              # Total number of particles
T = 1.1                 # Reduced temperature
beta = 1.0 / T           

# Move step sizes (to be tuned for ~50% acceptance)
dr_max = 0.084     # Max displacement (tuned for ~50% of acceptance)
dlnV_max = 0.01         # Max log-volume change
swap_attempts = int(0.02 * Ntot)  # ~2% of particles per cycle


def displacement_move(positions,box_length):
    """Random displacement of a single particle (NVT move)."""
    # Pick a random particle from all particles within the box. 
    
    i = randrange(len(positions)) # The index of a random particle (uniform distribution) within the array
    old_pos = positions[i].copy() # We store the position of the i-th particle as the old position
    
    # Compute local energy before move
    E_old = sum(pair_energy(dot((old_pos - positions[j] - round((old_pos-positions[j])/box_length)*box_length), 
                                   (old_pos - positions[j] - round((old_pos-positions[j])/box_length)*box_length)))
                for j in range(len(positions)) if j != i)
    
    # Propose displacement
    dr = (rand(3) * 2 - 1) * dr_max # We ensure that the displacement is uniform [-dr_max,dr_max] in all three directions.
    positions[i] = (old_pos + dr) % box_length # We ensure PBC by doing this... like PacMan 
    # Compute local energy after move
    new_pos = positions[i]
    E_new = sum(pair_energy(dot((new_pos - positions[j] - round((new_pos-positions[j])/box_length)*box_length), 
                                   (new_pos - positions[j] - round((new_pos-positions[j])/box_length)*box_length)))
                for j in range(len(positions)) if j != i)

    # =====================================================
    #  Metropolis criterion
    # =====================================================
    if random() < exp(-beta * (E_new - E_old)):
        # If an uniform random number within [0,1) is less than exp(-beta * DeltaE), the move is accepted.
        return True
    else:
        # Otherwise, we restore the old_position and return False
        positions[i] = old_pos
        return False

def volume_move(x,cajas, npart, beta, vmax):
    # Energía del estado actual
    vmin = 0.01 * vmax  # esto es para evitar que una caja colapse
    vmax_abs = 10 * vmax

    # Reversing logic
    pos1, pos2, L1, L2 = reverse_logic_adjustment(x,cajas)
    enlo1 = total_energy(pos1,L1)
    enlo2 = total_energy(pos2,L2)

    vol1 = cajas[0][0]**3
    vol2 =  cajas[1][0]**3
    vol_diff = vol2 - vol1

    # Random walk en ln(V1/V2)
    lnvn = log(vol1 / vol2) + (rand() - 0.5) * vmax/5
    v1n = vol1 * exp(lnvn) / (1 + exp(lnvn))
    v2n = vol1 + vol2 - v1n
    if v1n < vmin or v2n < vmin or v1n > vmax_abs or v2n > vmax_abs: #Si una caja colapsó rechazamos
        return cajas[0][0], cajas[1][0]
    box1n = v1n**(1/3)
    box2n = v2n**(1/3)
    # Reescalado de posiciones
    for i in range(npart):
        if x[1][i] == 0:
            factor = box1n / cajas[0][0]
        else:
            factor = box2n / cajas[1][0]
        x[0][i] *= factor

    # Energía en la nueva configuración
    pos1, pos2, L1, L2 = reverse_logic_adjustment(x,cajas)
    en1n = total_energy(pos1,L1)
    en2n = total_energy(pos2,L2)

    #Exponente del criterio de aceptacion
    arg1 = -beta * (en1n - enlo1) + (cajas[0][1] + 1) * log(v1n / vol1) / beta
    arg2 = -beta * (en2n - enlo2) + (cajas[1][1] + 1) * log(v2n / vol2) / beta

    # Regla de aceptación
    if rand() > exp(arg1 + arg2):
        # Rechazado: restaurar configuración antigua
        for i in range(npart):
            if x[1][i] == 0:
                factor = cajas[0][0] / box1n
            else:
                factor = cajas[1][0]/ box2n
            x[0][i] *= factor
        return cajas[0][0], cajas[1][0]# Retorna los valores originales
    print("Pass")
    return box1n, box2n  # Aceptado: retorna nuevos valores

def transfer_move():
    return None
