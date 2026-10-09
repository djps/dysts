# Rayleigh–Plesset in dysts — Phase 5: Physics caveats and validation

Oct 6, 2026 · @David Sinden

## Goal

This phase lists what limits how far the 1 MHz to 25 kHz results can be trusted, and the checks that confirm the implementation is right. Track these alongside Phase 4 rather than after it.

## Thermal behaviour (gamma)

The class uses a single polytropic exponent, `gamma`, with a default of 1.66666 (adiabatic, monatomic gas). For R0 ≈ 2 µm, the thermal Péclet number ωR0²/D_g drops from about 1 at 1 MHz to about 0.03 at 25 kHz. Over the sweep the gas therefore moves from near-adiabatic towards isothermal.

Either make `gamma` frequency-dependent (Prosperetti's effective polytropic exponent with thermal damping), or report every result at both γ = 1 and the adiabatic value.

## Cost and stiffness

The number of steps per driving period scales roughly as f0/f, because the afterbounces after each collapse have to be resolved. Expect 25 kHz runs to cost tens of times more than 1 MHz runs.

`make_trajectory` passes the analytic `_jac` to Radau (Phase 2 Step 8). At the 2 MHz default, neither numba nor the analytic Jacobian speeds things up much: 2000 steps take 7.8–8.0 s in `~/venv3.14`, against about 9 s without numba. The cost is dominated by Python overhead per right-hand-side call (`DynSys.rhs` rebuilds its parameter list every call), and by the hardcoded `rtol=1e-12, atol=1e-14` in `utils.integrate_dyn`. Those tolerances cannot be overridden: passing `rtol` through raises a duplicate-keyword error. For the sweep, either make the tolerances arguments of `integrate_dyn`, or call `solve_ivp` directly with a closure over a fixed parameter tuple and the numba `_rhs`/`_jac`. Check the effect of looser tolerances on a few Poincaré sections before relying on them. `make_trajectory_ensemble` cannot parallelise a parameter sweep, so use `multiprocessing` directly (Phase 4 Step 20).

## What the model leaves out

- **Higher-order compressibility.** The class has a first-order (R/c) correction, but it is not full Keller–Miksis or Gilmore. That matters for violent collapse at the low-frequency end (see Phase 2 Step 10).
- **Shape instabilities**, rectified diffusion and fragmentation. Small bubbles driven hard at low frequency would often not survive in reality.
- **Vapour pressure, heat and mass transfer.**

State these limits when interpreting chaotic regimes at low frequency.

## Validation

1. After the Phase 1 Step 3 fix, reproduce the linear resonance curve at small `P` and check its peak against the formula in Phase 1 Step 5.
2. Reproduce one published bifurcation diagram, using the same model variant the paper used.
3. Cross-check trajectories with a second solver. `make_trajectory(..., method="DOP853")` and `method="LSODA"` both work, because `method` is passed through to `solve_ivp`.
4. Cross-check a few Lyapunov exponents from `analysis.find_lyapunov_exponents` against an independent code such as Julia's DynamicalSystems.jl. The dysts routine uses a first-order tangent step with finite-difference Jacobians.

## Done when

- [ ] `gamma` treatment decided and applied consistently
- [ ] Linear resonance check passes
- [ ] One published bifurcation diagram reproduced
- [ ] Lyapunov exponents agree with an independent solver
- [ ] Model limits written into the results write-up
