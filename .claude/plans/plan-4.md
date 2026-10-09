# Rayleigh–Plesset in dysts — Phase 4: Frequency sweep

Oct 6, 2026 · @David Sinden

## Goal

This phase measures how the bubble's dynamics change as the driving frequency falls from 1 MHz to 25 kHz at fixed acoustic pressure.

The sweep is a separate script; it adds nothing to the database. The JSON entry holds one canonical parameter set, and the script changes attributes on one instance:

- `DynSys.rhs` reads each parameter with `getattr` on every call, so `model.omega = ...` takes effect straight away.
- `model.params` is **not** updated, because there is no custom `__setattr__`. Never read parameters back from `model.params`.
- `make_trajectory` has no `dt` or `init_cond` arguments. Set `model.dt`, `model.ic` and `model.period` instead.

## Step 16 — Frequency grid

Use 200–400 log-spaced frequencies from 1 MHz down to 25 kHz, with `omega = 2 * np.pi * f` in rad/s. Hold `P` fixed across the sweep, which means a fixed acoustic pressure.

## Step 17 — Integrate each frequency

Integrate on a uniform grid with dt = T/K, where T = 1/f and K is an integer. Every K-th sample is then a stroboscopic (Poincaré) point.

```python
import numpy as np
from dysts.flows import RayleighPlesset
from dysts.utils import integrate_dyn

model = RayleighPlesset()
K, n_trans, n_keep = 400, 200, 500          # samples/period, periods

def near_core(t, y):                        # stop if R approaches the hard core
    return y[0] - 1.05 * model.a
near_core.terminal = True

def run(f, ic):
    T = 1.0 / f
    model.omega = 2 * np.pi * f
    model.dt, model.period = T / K, T       # period is used by analysis.* helpers
    n = (n_trans + n_keep) * K
    tpts = np.arange(n) * model.dt
    traj = integrate_dyn(model, ic, tpts, dtval=model.dt,
                         method="Radau", events=near_core).T
    ok = traj.shape[0] == n                 # False if solve_ivp stopped early
    kept = traj[n_trans * K:]
    poincare = kept[::K, 0]                 # R at fixed phase
    return ok, kept, poincare, traj[-1]     # last state seeds the next run
```

- The code calls `utils.integrate_dyn` directly rather than `make_trajectory`. `make_trajectory` passes only `method` to the solver, so it cannot pass `events`. `integrate_dyn` passes any extra keyword arguments to `solve_ivp`. It returns a (d, n) array, transposed above, with no postprocessing.
- `make_trajectory(n, resample=False, postprocess=False)` gives the same grid, `arange(n) * model.dt`, if events are not needed. With a 1-D `ic` it does **not** warn when integration stops early, so check `traj.shape[0]`.
- Check for and record failed runs (`ok == False`); do not drop them silently.
- Discard about 200 or more periods as transient, then keep 500 or more.

## Step 18 — Continuation in both directions

Seed each frequency with the final state of the previous one, and run both 1 MHz → 25 kHz and 25 kHz → 1 MHz. Ω carries over unchanged, so the drive phase stays continuous when ω changes. Differences between the two passes expose coexisting attractors and hysteresis, which are common in bubble oscillators.

## Step 19 — Outputs at each frequency

| Output | How | What it shows |
| --- | --- | --- |
| Bifurcation diagram | Poincaré R/R0 values plotted against f | period-n windows, period doubling, chaos; the headline figure |
| Rmax/R0, Rmin/R0 | max and min of R over kept periods | response amplitude and compression ratio |
| Largest Lyapunov exponent | `analysis.find_lyapunov_exponents(model, ...)` after setting `model.ic`, `model.period` and `model.omega`; result in s⁻¹, so also report λ·T per drive period | chaos (positive) vs periodic |
| 0–1 test K | not in the repo; implement it (Gottwald–Melbourne) on the Poincaré samples | independent chaos check |
| Correlation dimension | `lyap.corr_dim` | attractor complexity |
| Spectrum | `utils.find_psd` on R(t), or on the radiated pressure ∝ R(RU̇ + 2U²) | subharmonics f/2, f/3, ultraharmonics, broadband content |
| Peak Mach number | max of abs(U) / `model.c` | flags model validity; unreliable above about 0.1–0.3 |

`find_lyapunov_exponents` integrates its own trajectory with `resample=True`, sized by `model.period`. It uses finite-difference Jacobians, not `_jac`, and a first-order tangent step. Convergence in `pts_per_period` has to be checked at the low-frequency end.

## Step 20 — Map the (f, P) plane

Chaos in RP lives in tongues in the frequency–pressure plane. A single fixed `P` passes in and out of chaos across a 40× frequency range, so a 1D sweep alone can mislead.

Run a coarse grid of `P` values, for example 30 kPa to 1 MPa, which includes the current default. Colour each cell by the largest Lyapunov exponent or the 0–1 test K, and mark which tongues the 1D sweep crosses.

`make_trajectory_ensemble` does not help with this. It builds fresh default instances by class name, and its `use_multiprocessing` is not implemented. Use `multiprocessing.Pool` over `P` values, or over the two continuation directions, with each worker running its own continuation.

## Done when

- [ ] Sweep script runs end to end on a coarse grid
- [ ] Forward and backward continuation both completed
- [ ] Bifurcation diagram and Lyapunov-vs-f plot produced
- [ ] Frequencies with a high Mach number or failed runs are flagged, not dropped silently
- [ ] (f, P) map produced at coarse resolution
