# Source audit and validation

## Provenance

The source package consists of `Initial_guess.m`, `Nonlinear_iteration.m`, and the December 2025 LOG(M) poster. The two MATLAB files in `archive/original/` are unchanged copies. The poster and its preview are original research artifacts, not figures produced by the revised demo.

The files in `src/`, `examples/`, and `tests/`, together with the documentation, are later repository preparation. They should not be used to attribute additional work or numerical results to any particular original coauthor.

## Changes from the supplied code

| Source observation | Why it matters | Revised treatment |
|---|---|---|
| The bottom ghost equation in `Nonlinear_iteration.m` refers to the physical boundary row, although the intended centered derivative requires the next physical row. | It is not the centered stencil described in the comment. Homogeneous data and a y-independent solution can hide the error. | Eliminate the bottom ghost using the second physical row and the correct `2p` factor. |
| Four ghost-corner rows are explicitly marked placeholders and use boundary values as matrix diagonal coefficients. | A zero side value creates a zero equation coefficient; other values produce unrelated reciprocal corner values. | Remove all ghost unknowns; apply identity rows at physical side boundaries. |
| `Initial_guess.m` sets one-step differences equal to derivative data without p scaling, and uses an inward difference at the top. | Nonzero derivatives have inconsistent units and sign between edges. | Use the same positive-y derivative convention and ghost-eliminated operator for initialization and nonlinear iterations. |
| The supplied stopping test uses relative iterate change alone. | A small change does not guarantee a small nonlinear residual, especially under relaxation. | Check both update and equation residual; expose nonconvergence in the returned status. |
| No execution driver, complete viscosity sequence, or input checks accompany the two functions. | Exact historical replay and safe reuse are limited. | Add a runnable example, explicit coordinates, validation, and manufactured tests. |

These corrections address general boundary data. They do not establish that the poster's qualitative homogeneous-boundary observations are wrong. The original +1/−1 side values and nearly y-independent solution can mask several boundary defects.

## What was and was not executed

MATLAB and GNU Octave were not installed in the packaging environment. **Neither the original nor revised MATLAB solver was executed.** The MATLAB test suite is provided for the next MATLAB run. It includes a deliberate one-iteration run to verify that a stopped solve is not labeled converged.

An independent NumPy/SciPy implementation was executed locally. It assembles the stated equations and also evaluates a separate physical ghost stencil to check the residual. This validates the numerical construction and the chosen example at those parameters; it does not run, transpile, or certify the `.m` files.

To rerun the independent check in a Python environment with NumPy and SciPy:

```sh
python3 tests/reference_check.py
```

The reference script uses stricter update stopping for the demo than the MATLAB example, so its iteration counts need not match a later MATLAB run. Machine-readable output is in [reference-check-results.json](reference-check-results.json).

## Observed reference-check results

For `u=y^2+0.3y+0.2x+0.1`, `epsilon=0.4`, and `f=0.2u−2`, the 13 × 11 grid gave maximum solution error **6.34e−14**. Bottom/top positive-y derivatives are 0.3/2.3; this exercises nonzero boundary signs and spacing. The independently evaluated equation residual was **3.55e−13**.

For the smooth manufactured solution `u=0.25 sin(pi*x) cos(pi*y)+0.2x+0.1y` on the unit square:

| Points per axis | Maximum solution error | Observed order |
|---:|---:|---:|
| 9 | 3.18019e−3 | — |
| 17 | 7.90444e−4 | 2.008 |
| 33 | 1.97325e−4 | 2.002 |

The independent stencil residual was below **3.8e−12** on all three grids. A zero-data case also returned the zero solution.

For the 121 × 31 homogeneous-boundary continuation example:

| Epsilon | Accepted iterations | Independent stencil residual |
|---:|---:|---:|
| 0.2 | 25 | 8.35e−12 |
| 0.1 | 24 | 2.00e−11 |
| 0.05 | 26 | 4.02e−11 |

Vertical variation stayed below 2.2e−14, as expected for these y-independent boundary values and a symmetric solution. This supports the implementation's handling of the simple internal-layer example; it is not evidence of uniqueness, a rigorous vanishing-viscosity limit, or general two-dimensional performance.

## Next validation step

Run `tests/run_tests.m` and `examples/demo_internal_layer.m` in MATLAB, record the MATLAB version and results, and compare the returned residual and refinement orders. Further research checks should include nontrivial y-dependent boundary conditions, continuation sensitivity, stronger x refinement, and a conservative/upwind discretization comparison where appropriate.
