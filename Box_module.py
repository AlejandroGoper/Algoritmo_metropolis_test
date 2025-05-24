from numpy import ceil, array

# =============================================================================
#  Initialize Positions within Boxes
# =============================================================================

def init_positions(N, V):
    """Place N particles at approximate lattice positions then randomize."""
    box_length = V ** (1/3)
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
