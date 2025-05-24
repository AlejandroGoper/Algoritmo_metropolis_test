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
