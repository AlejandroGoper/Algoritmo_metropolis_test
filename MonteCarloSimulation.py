#Importamos los modulos de energia
from Energy_module import LennardJones
from Box_module import reverse_logic_adjustment
from random import randrange, random
from numpy import round, dot, log, delete, vstack
from numpy.random import rand
from math import exp


# =============================================================================
#       Monte Carlo Simulation class 
# =============================================================================


class MonteCarloSimulation(LennardJones):
    """
    A class for performing Monte Carlo moves (displacement, volume, and particle transfer)
    with parameters configurable at initialization. Inherits Lennard-Jones parameters and methods.
    """
    def __init__(self, particles=500, red_temperature=1.1,LJ_epsilon=1.0, LJ_sigma=1.0, LJ_rcut=2.5, dr_max=0.084, dlnV_max=0.01):
        self.Ntot = particles
        self.T = red_temperature    
        self.beta = 1.0 / self.T
        # Initialize Lennard-Jones base (epsilon, sigma, rcut)
        super().__init__(LJ_epsilon, LJ_sigma, LJ_rcut)
        # Move step sizes
        self.dr_max = dr_max
        self.dlnV_max = dlnV_max


    def displacement_move(self,positions,box_length):
        """Random displacement of a single particle (NVT move)."""
        # Pick a random particle from all particles within the box. 
        
        i = randrange(len(positions)) # The index of a random particle (uniform distribution) within the array
        old_pos = positions[i].copy() # We store the position of the i-th particle as the old position
        
        # Compute local energy before move
        E_old = 0.0
        for j in range(len(positions)):
            if j == i: continue
            dr = old_pos - positions[j]
            dr -= round(dr/box_length) * box_length
            r2 = dot(dr, dr)
            if r2 < self.rcut**2:
                E_old += self.pair_energy(r2)
        # Propose displacement
        dr = (rand(3) * 2 - 1) * self.dr_max # We ensure that the displacement is uniform [-dr_max,dr_max] in all three directions.
        positions[i] = (old_pos + dr) % box_length # We ensure PBC by doing this... like PacMan 
        # Compute local energy after move
        new_pos = positions[i]
        E_new = 0.0
        for j in range(len(positions)):
            if j == i: continue
            dr = new_pos - positions[j]
            dr -= round(dr/box_length) * box_length
            r2 = dot(dr, dr)
            if r2 < self.rcut**2:
                E_new += self.pair_energy(r2)
        # =====================================================
        #  Metropolis criterion
        # =====================================================
        if random() < exp(-self.beta * (E_new - E_old)):
            # If an uniform random number within [0,1) is less than exp(-beta * DeltaE), the move is accepted.
            return True
        else:
            # Otherwise, we restore the old_position and return False
            positions[i] = old_pos
            return False

    def volume_move(self,x,cajas, npart, beta, vmax):
        # Energía del estado actual
        vmin = 0.01 * vmax  # esto es para evitar que una caja colapse
        vmax_abs = vmax

        # Reversing logic
        pos1, pos2, L1, L2 = reverse_logic_adjustment(x,cajas)
        enlo1 = self.total_energy(pos1,L1)
        enlo2 = self.total_energy(pos2,L2)

        vol1 = cajas[0][0]**3
        vol2 =  cajas[1][0]**3
        vol_diff = vol2 - vol1

        # Random walk en ln(V1/V2)
        lnvn = log(vol1 / vol2) + (2*rand() - 1) * self.dlnV_max # Adjusting logic to tune the parameter for 50% Acceptance
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
        en1n = self.total_energy(pos1,L1)
        en2n = self.total_energy(pos2,L2)

        #Exponente del criterio de aceptacion
        arg1 = -self.beta * (en1n - enlo1) + (cajas[0][1] + 1) * log(v1n / vol1) / self.beta
        arg2 = -self.beta * (en2n - enlo2) + (cajas[1][1] + 1) * log(v2n / vol2) / self.beta

        # Regla de aceptación
        if rand() > exp(arg1 + arg2):
            # Rechazado: restaurar configuración antigua
            print("Rechazado")
            for i in range(npart):
                if x[1][i] == 0:
                    factor = cajas[0][0] / box1n
                else:
                    factor = cajas[1][0]/ box2n
                x[0][i] *= factor
            return cajas[0][0], cajas[1][0], False # Retorna los valores originales
        return box1n, box2n, True  # Aceptado: retorna nuevos valores

    def transfer_move(self,pos1, pos2, L1, L2):
        """
        Particle transfer move matching the paper's Eq. (7) form:
        ΔW_rev = ΔE_local + kT * [ 
            N_recv * ln((N_recv+1)/N_recv) 
        + N_don  * ln((N_don-1)/N_don)
        + ln(V_don/(N_don-1)) 
        - ln(V_recv/(N_recv+1))
        ]
        """
        # Decide donor and receiver
        if random() < 0.5:
            donor, receiver, Ld, Lr = pos2, pos1, L2, L1
            donor_is_box2 = True
        else:
            donor, receiver, Ld, Lr = pos1, pos2, L1, L2
            donor_is_box2 = False

        Nd = len(donor)
        Nr = len(receiver)
        Vd = Ld ** 3
        Vr = Lr ** 3

        # Local energy removal from donor
        idx = randrange(Nd)
        old_pos = donor[idx].copy()
        E_old = 0.0
        for j in range(Nd):
            if j == idx:
                continue
            dr = old_pos - donor[j]
            dr -= round(dr / Ld) * Ld
            r2 = dot(dr, dr)
            if r2 < self.rcut**2:
                E_old += self.pair_energy(r2)

        # Remove particle
        donor_new = delete(donor, idx, axis=0)

        # Local energy insertion into receiver
        trial = rand(3) * Lr
        E_test = 0.0
        r_overlap = 0.0001 # Condition for no overlapping
        for j in range(Nr):
            dr = trial - receiver[j]
            dr -= round(dr / Lr) * Lr
            r2 = dot(dr, dr)
            if r2 < r_overlap**2:
                print("Overlap")
                # Immediate reject: overlap too severe
                return pos1, pos2, L1, L2, False
        
            if r2 < self.rcut**2:
                E_test += self.pair_energy(r2)
        receiver_new = vstack([receiver, trial])

        # Entropic term from paper's Eq. (7)
        delta_log = (
            Nr * log((Nr + 1) / Nr)
        + Nd * log((Nd - 1) / Nd)
        + log(Vd / (Nd - 1))
        - log(Vr / (Nr + 1))
        )

        # =====================================================
        #  Metropolis criterion
        # =====================================================
        if random() < exp(-self.beta * (E_test - E_old) - delta_log):
            # Assign back to pos1/pos2 depending on donor box
            if donor_is_box2:
                return receiver_new, donor_new, Lr, Ld, True  # pos2 replaced
            else:
                return donor_new, receiver_new, Ld, Lr, True  # pos1 replaced
        else:
            # Reject: return originals
            return pos1, pos2, L1, L2, False