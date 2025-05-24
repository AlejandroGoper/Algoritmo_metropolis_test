# Reduced Temperature and Density

Here we consider the expression for the reduced density in a Lennard-Jones system:

$$
\begin{aligned}
\rho^* = \frac{N \sigma^3}{V},
\end{aligned}
$$

which tells us how “crowded” our particles in the original Box are.

We use $\rho^* = 0.6$  because it sits in the liquid–vapor coexistence region for temperatures around $T^* = 1.1$. 

It’s high enough that one of our boxes will settle into a liquid‐like density and the other into a vapor‐like density once the Gibbs moves do their work.

## Initial conditions 

By setting, for example $N=500$ and splitting these particules evenly at the start, we simply give each box the same initial density, then let the volume‐ and particle‐exchange moves 
drive them toward the two coexisting densities (one high, one low).

To calculate the volume of the box, we simply use the definition of numerical density:

$$
\begin{aligned}
V = \frac{N}{\rho_{\text{total}}}.
\end{aligned}
$$
## Function displacement_move(x,cajas, npart, beta, vmax)

x is a tuple of the form [np.zeros((npart, 3)), np.zeros(npart)] the first entry is an array whith the positions of the particles and the second has a value of either 1 or 0 indicating the box. 

cajas is a tuple with the entries [np.array([caja1L, 0]), np.array([caja2L, 0])] cajas[0], cajas[1] have the lenght and number of particles of each box

npart is the number of particles, beta is 1/T where T is the reduced temperature and vmax is the total volunm

# Parameters of the Monte Carlo Moves 


## Displacement move
Every time we try a particle‐displacement move, we pick a random vector $\Delta r$ whose components lie uniformly in $\left[ -\Delta r_{\text{Max}}, \Delta r _{\text{Max}}\right]$.

If $\Delta r _{\text{Max}}$ is too large, almost all proposed moves will land in high‐energy overlap configurations and be rejected, which could affect the acceptance ratio.

In the opposite case, nearly every move is accepted, but we only explore configuration space very slowly.

We’ll typically tune $\Delta r_{\text{Max}}$ by running a short test and adjusting up or down until we hit that target (the one that gives us 50% of acceptance ratio).

## Transfer move 
```
swap_attempts = int(0.02 * Ntot)
```

This sets how many particle‐transfer moves we attempt in each Monte Carlo cycle—here, about 2% of our total particle count per cycle (that's what the article says). 

## Volume move

Pending...



