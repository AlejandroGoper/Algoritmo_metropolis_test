from numpy import round, dot
from math import pi

# Lennard-Jones parameters (reduced units) 
epsilon = 1.0
sigma = 1.0
rcut = 2.5 * sigma

# ========================================================================================
#  We work with square distances instead of the usual to improve efficiency and save time
# ========================================================================================
def pair_energy(r2):
    """Lennard-Jones pair energy for squared distance r2."""
    inv_r6 = (sigma**2 / r2)**3
    return 4 * epsilon * (inv_r6**2 - inv_r6)

# ========================================================================================
#          Compute total LJ energy with cutoff and simple tail correction.
# ========================================================================================
def total_energy(positions, box_length):
    E = 0.0
    N = len(positions)
    for i in range(N-1):
        for j in range(i+1, N):
            dr = positions[i] - positions[j]
            # Periodic boundary conditions
            dr -= round(dr / box_length) * box_length
            r2 = dot(dr, dr)
            if r2 < rcut**2:
                E += pair_energy(r2)
    # Simple mean-field tail correction
    rho = N / box_length**3
    E += (8.0/3.0) * pi * rho * N * epsilon * ( (sigma/rcut)**9 / 3 - (sigma/rcut)**3 )
    return E
