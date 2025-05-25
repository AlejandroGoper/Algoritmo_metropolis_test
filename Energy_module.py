from numpy import round, dot
from math import pi


class LennardJones:
    def __init__(self, epsilon = 1.0, sigma = 1.0, rcut = 2.5):
        # Lennard-Jones parameters (reduced units) 
        self.epsilon = epsilon
        self.sigma = sigma 
        self.rcut = rcut

    # ========================================================================================
    #  We work with square distances instead of the usual to improve efficiency and save time
    # ========================================================================================
    def pair_energy(self,r2):
        """Lennard-Jones pair energy for squared distance r2."""
        inv_r6 = (self.sigma**2 / r2)**3
        return 4 * self.epsilon * (inv_r6**2 - inv_r6)

    # ========================================================================================
    #          Compute total LJ energy with cutoff and simple tail correction.
    # ========================================================================================
    def total_energy(self, positions, box_length):
        E = 0.0
        N = len(positions)
        for i in range(N-1):
            for j in range(i+1, N):
                dr = positions[i] - positions[j]
                # Periodic boundary conditions
                dr -= round(dr / box_length) * box_length
                r2 = dot(dr, dr)
                if r2 < self.rcut**2:
                    E += self.pair_energy(r2)
        # Simple mean-field tail correction
        rho = N / box_length**3
        E += (8.0/3.0) * pi * rho * N * self.epsilon * ( (self.sigma/self.rcut)**9 / 3 - (self.sigma/self.rcut)**3 )
        return E
