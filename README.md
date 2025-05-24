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
