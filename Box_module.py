from numpy import ceil, array
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# =============================================================================
#  Initialize Positions within Boxes
# =============================================================================

def init_positions(N, V):
    """Place N particles at approximate lattice positions then randomize."""
    box_length = V ** (1/3) # The length of a cube of volume V
    # Simple cubic grid
    n_side = int(ceil(N ** (1/3)))
    spacing = box_length / n_side
    positions = []
    for x in range(n_side):
        for y in range(n_side):
            for z in range(n_side):
                if len(positions) < N:
                    positions.append(array([ (x+0.5)*spacing, 
                                                (y+0.5)*spacing, 
                                                (z+0.5)*spacing ]))
    positions = array(positions)
    return positions, box_length

def plot_boxes(pos1, L1, pos2, L2):
    """
    Plot two 3D scatter subplots of particle positions in Box 1 and Box 2.
    
    Parameters:
    - pos1: (N1 x 3) array of positions in Box 1
    - L1: side length of Box 1
    - pos2: (N2 x 3) array of positions in Box 2
    - L2: side length of Box 2
    """
    fig = plt.figure(figsize=(12, 6))

    # Box 1
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.scatter(pos1[:, 0], pos1[:, 1], pos1[:, 2], marker='o')
    ax1.set_xlim(0, L1); ax1.set_ylim(0, L1); ax1.set_zlim(0, L1)
    ax1.set_title('Box 1')
    ax1.set_xlabel('x'); ax1.set_ylabel('y'); ax1.set_zlabel('z')

    # Box 2
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.scatter(pos2[:, 0], pos2[:, 1], pos2[:, 2], marker='o')
    ax2.set_xlim(0, L2); ax2.set_ylim(0, L2); ax2.set_zlim(0, L2)
    ax2.set_title('Box 2')
    ax2.set_xlabel('x'); ax2.set_ylabel('y'); ax2.set_zlabel('z')

    plt.tight_layout()
    plt.show()

import numpy as np

def logic_adjustment(pos1, pos2, L1, L2):
    """
    Convert Alex's simulation data (pos1, pos2, L1, L2) into the inputs required
    by Galo's function volume_move(x, cajas, npart, beta, vmax).

    Parameters:
    - pos1: (N1×3) array of positions in Box 1
    - pos2: (N2×3) array of positions in Box 2
    - L1, L2: side lengths of Box 1 and Box 2
    - beta: 1/(k_B T), inverse temperature

    Returns:
    - x: list [positions_array, box_index_array] where
        positions_array.shape == (N1+N2, 3)
        box_index_array.shape == (N1+N2,) entries 0 or 1
    - cajas: [[L1, N1], [L2, N2]]
    - npart: total number of particles (N1+N2)
    - vmax: total volume = V1 + V2
    """
    # Stack all positions into one (Ntot×3) array
    positions = np.vstack([pos1, pos2])
    # Create box-index array: 0 for first N1 entries, 1 for next N2 entries
    N1 = pos1.shape[0]
    N2 = pos2.shape[0]
    box_indices = np.array([0]*N1 + [1]*N2)
    # Package x as expected
    x = [positions, box_indices]
    # cajas holds [L, particle count] for each box
    cajas = [[L1, N1], [L2, N2]]
    # Total particles
    npart = N1 + N2
    # Total volume
    vmax = L1**3 + L2**3
    return x, cajas, npart, vmax

def reverse_logic_adjustment(x, cajas):
    """
    Convert the inputs used by Galo back into
    the one's written by Alex pos1, pos2, L1, L2 : expected by our total_energy function.

    Parameters:
    - x: list [positions_array, box_index_array] where
        positions_array.shape == (Ntot, 3)
        box_index_array.shape == (Ntot,) entries 0 or 1
    - cajas: [[L1, N1], [L2, N2]]

    Returns:
    - pos1: (N1×3) array of positions in Box 1
    - pos2: (N2×3) array of positions in Box 2
    - L1: side length of Box 1
    - L2: side length of Box 2
    """
    positions, box_indices = x
    L1, N1 = cajas[0]
    L2, N2 = cajas[1]

    # Extract positions belonging to each box
    pos1 = positions[box_indices == 0]
    pos2 = positions[box_indices == 1]

    return pos1, pos2, L1, L2