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


