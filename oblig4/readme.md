# Project 4: Diffusion equation

Differences between Explicit and Implicit schemes:

- Explicit: Have the solution at $V_j$, use this to find the solution at $V_{j+1}$ as
\[ A v_{j} = v_{j+1] \]
throuh standard matrix-multiplication
- Implicit: Have the solution at $V_j$, use this to find the solution at $V_{j+1}$ as
\[ A v_{j+1} = V_{j} \]

In the example on p.313, following is done:

1. Find the solution: RHS of $Av_j = Bv_{j-1} = r$
2. Solve, using tridiag-routine $Av_j = r$ for $v_j$.

Current problem:

- *Why does the implicit scheme give 0 for the two final elements?*
- *Solve analytically for $v(x,t)$ in the differential eq., using that the boundary conditions are 0 in both ends.*
- *Remember that the solution is concave downwards below 0.*

 
