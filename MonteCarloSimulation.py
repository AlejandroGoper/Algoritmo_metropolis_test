#Importamos los modulos de energia
from Energy_module import LennardJones
from Box_module import reverse_logic_adjustment
from random import randrange, random
from numpy import round, dot, log, delete, vstack,exp
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

    def volume_move(self, x, cajas, npart, beta, dV_max):
        # Obtener volúmenes y número de partículas
        L1_old = cajas[0][0]
        L2_old = cajas[1][0]
        V1_old = L1_old ** 3
        V2_old = L2_old ** 3
        V_total = V1_old + V2_old
        N1 = cajas[0][1]
        N2 = cajas[1][1]
        vmin = 0.01 * V_total

        # Proponer nuevo volumen para caja 1
        deltaV = (2 * rand() - 1) * dV_max
        V1_new = V1_old + deltaV
        V2_new = V_total - V1_new

        # Verificación de volúmenes mínimos
        if V1_new < vmin or V2_new < vmin:
            return L1_old, L2_old, False

        # Reescalado de longitudes de caja
        L1_new = V1_new ** (1/3)
        L2_new = V2_new ** (1/3)

        # Reescalado de posiciones
        for i in range(npart):
            if x[1][i] == 0:
                factor = L1_new / L1_old
            else:
                factor = L2_new / L2_old
            x[0][i] *= factor

        # Calcular energías nueva y vieja
        pos1, pos2, L1_tmp, L2_tmp = reverse_logic_adjustment(x, [[L1_new, N1], [L2_new, N2]])
        U1_new = self.total_energy(pos1, L1_tmp)
        U2_new = self.total_energy(pos2, L2_tmp)

        pos1_old, pos2_old, L1_old_tmp, L2_old_tmp = reverse_logic_adjustment(x, [[L1_old, N1], [L2_old, N2]])
        U1_old = self.total_energy(pos1_old, L1_old_tmp)
        U2_old = self.total_energy(pos2_old, L2_old_tmp)

        # Regla de aceptación (Ecuación 8.3.2)
        dU = (U1_new + U2_new) - (U1_old + U2_old)
        ln_ratio = (
            N1 *log(V1_new / V1_old)
            + N2 * log(V2_new / V2_old)
            - beta * dU
        )
        acc_prob = min(1.0,exp(ln_ratio))

        if np.random.rand() < acc_prob:
            # Aceptado
            return L1_new, L2_new, True
        else:
            # Rechazado → restaurar las posiciones
            for i in range(npart):
                if x[1][i] == 0:
                    factor = L1_old / L1_new
                else:
                    factor = L2_old / L2_new
                x[0][i] *= factor
            return L1_old, L2_old, False


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
