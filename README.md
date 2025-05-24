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
