# Numerical method

## Model and conventions

The supplied poster studies `u*u_x - u_yy = 0` on `[0,6] × [0,5]`. Away from `u=0`, interpreting x as an evolution variable gives `u_x = u_yy/u`; the diffusion coefficient changes sign with u. The regularized equation solved here is

```text
u*u_x - u_yy - epsilon*u_xx = f, epsilon > 0.
```

The research example has `f=0`. A nonzero source is included solely to support known-solution verification and reuse. For positive epsilon, the second-order part is elliptic. The code does not directly solve the epsilon-zero equation.

`x` and `y` are increasing uniform coordinate vectors that include their endpoints. `nx=numel(x)` and `ny=numel(y)` count **physical points, including boundaries**. Output `U` is `ny × nx`; `U(j,i)` approximates `u(x(i),y(j))`. Sparse-system numbering is `i + (j-1)*nx` (x varies fastest), implemented in MATLAB with `reshape(U.',[],1)`.

The poster writes indices from 0 through Nx/Ny, while the supplied functions take Nx/Ny as point counts. This repository uses explicit coordinate vectors to remove that ambiguity.

## Centered differences

At non-side-boundary nodes, with spacings h in x and p in y:

```text
D_x U   = (U(j,i+1) - U(j,i-1)) / (2h)
D_xx U  = (U(j,i+1) - 2U(j,i) + U(j,i-1)) / h^2
D_yy U  = (U(j+1,i) - 2U(j,i) + U(j-1,i)) / p^2
R(U)    = U .* D_x U - D_yy U - epsilon*D_xx U - f
```

This discretizes the nonconservative product `u*u_x`, as in the source project; it does not substitute a conservative numerical flux for `d(u^2/2)/dx`. No claim of shock-capturing accuracy or a proven entropy limit is made.

## Boundary conditions and ghost elimination

Left and right values are Dirichlet data of length ny. Their rows are identity equations. They take precedence at physical corners. Horizontal data are length nx and represent the positive-y derivative on **both** edges:

```text
bottom(i) = du/dy(x(i), y(1))
top(i)    = du/dy(x(i), y(end))
```

If the input is an outward-normal derivative, negate it on the bottom edge before calling the solver. The code does not apply an extra epsilon factor to these derivatives.

At the bottom, the centered Neumann relation gives `U_ghost = U(2,i) - 2p*bottom(i)`. Substitution into `-D_yy U` yields `2*(U(1,i)-U(2,i))/p^2 + 2*bottom(i)/p`. The diffusion matrix therefore has a doubled inward neighbor coefficient and the right-hand side receives **−2*bottom(i)/p**. At the top the right-hand side receives **+2*top(i)/p**, with the analogous doubled inward neighbor.

Only physical points are stored. There are no ghost-corner unknowns or placeholder corner equations. Endpoint entries of `bottom` and `top` are ignored because side Dirichlet values own those corners. Users should supply corner-compatible data for a smooth continuous solution.

The centered ghost derivative is second-order, but substituting an exact derivative into the boundary PDE stencil can produce first-order local truncation error for general smooth solutions with nonzero third y derivative. The reported smooth manufactured test displays second-order **global** solution error; that observation is not a universal accuracy proof.

## Initialization, iteration, and diagnostics

1. Assemble the sparse diffusion operator `L` and boundary-adjusted right-hand side `b`.
2. Solve `L*u0=b`, or accept a user-supplied physical-grid warm start. The original homogeneous research problem thus begins with an elliptic ramp.
3. Freeze the advecting coefficient at the current iterate: `(L + C(uk))*candidate = b`, where `C(uk)` applies `uk*D_x` at non-side nodes.
4. Backtrack from the configured relaxation toward a minimum step. Accept a reduction in the nonlinear residual, or a residual already below tolerance.
5. Require **both** relative update and absolute nonlinear residual tolerances. An initial state already satisfying the residual tolerance is accepted without a redundant solve.

The update norm is `norm(u_next-u,inf)/max(1,norm(u,inf))`. The residual is `norm((L+C(u))*u-b,inf)` and includes side-boundary equations. It is an absolute discrete-equation tolerance; it mixes PDE and Dirichlet rows with their physical units and is not a continuum error estimate. Choose a tolerance appropriate to the problem's scaling.

Each sparse solve is checked for nonfinite values and excessive scaled linear residual. Inputs are checked for finite real data, compatible dimensions, positive viscosity, uniform increasing grids, and valid iteration settings. Failure returns `info.converged=false`, a warning, and one of `maximum_iterations`, `line_search_stalled`, `nonfinite_linear_solve`, or `inaccurate_linear_solve`. Backtracking improves robustness but does not guarantee global convergence.

`info.history` records iteration, relative update, nonlinear residual, and accepted relaxation. `info.max_cell_peclet = max(abs(U))*h/(2*epsilon)` warns above 1: beyond this threshold the frozen central x stencil can lose its nonpositive off-diagonal coefficients. This is a resolution diagnostic, not a proof of convergence or accuracy. The poster's practical rule `h,p < 2*epsilon` is retained as historical context; the revised diagnostic concerns h and the actual magnitude of u.

## Example sizes and the original experiment

The poster reports `Nx=320`, `Ny=270`, decreasing epsilon from 0.1 to 0.01, at most 30 nonlinear iterations at each stage, and a relative-change threshold of `1e-3`. The precise continuation schedule and driver were not supplied, so those runs cannot be reproduced exactly from the two archived functions alone.

The new demo is intentionally modest: 121 × 31 physical points, epsilon `[0.2,0.1,0.05]`, and maximum cell Peclet 0.5 at the smallest epsilon for values bounded by 1. It uses continuation and checks convergence before plotting. To use epsilon 0.01 at this solution scale, refine x to h ≤ 0.02 (at least 301 physical x points across the length-6 domain), then assess accuracy by additional refinement. The inequality alone does not resolve every accuracy concern.
