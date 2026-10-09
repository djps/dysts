# Rayleigh–Plesset in dysts — Phase 2: Code

Oct 6, 2026 · @David Sinden

## Goal

This phase brings the existing `RayleighPlesset` class in `dysts/flows.py` (lines 18–67) into line with the other dysts flows, and adds the plumbing the sweep needs. The conventions are: parameters passed positionally in sorted order, `@staticjit` (numba is optional; without it the code runs as plain Python), and an optional `_postprocessing` hook.

## Step 6 — Fix the RHS

Apply the Phase 1 Step 3 fix to the compressibility term. No other changes are needed to `_rhs`.

## Step 7 — Bounded phase (optional)

`Duffing` defines `_postprocessing` returning `(x, y, np.cos(z))`, so the unbounded phase is mapped to a bounded value. `RayleighPlesset` has no `_postprocessing`, so `make_trajectory` returns the raw, growing Ω. Adding one:

```python
    @staticjit
    def _postprocessing(R, U, Omega):
        return R, U, np.sin(Omega)
```

This changes the default output of `make_trajectory`, which has `postprocess=True`. It also switches `analysis.find_lyapunov_exponents` onto its coordinate-transform branch. The sweep (Phase 4) always passes `postprocess=False`, so it is not affected either way.

## Step 8 — Jacobian: rewrite, and wire it in only if wanted

The commented-out `_jac` has a different argument order from `_rhs` (`a, sigma, rho, mu, ...`). It also predates the compressibility term. It must be rewritten in sorted order before use.

No code in the repo calls `_jac`:

- `utils.integrate_dyn` calls `solve_ivp` without `jac=`, so Radau builds its own finite-difference Jacobian.
- `analysis.find_lyapunov_exponents` uses `utils.jac_fd`.

An analytic Jacobian only helps Radau if `DynSys.make_trajectory` (`dysts/base.py`) is changed to pass `jac=lambda t, y: self._jac(*y, t, *params)` to `integrate_dyn`. `integrate_dyn` already forwards extra keyword arguments to `solve_ivp`.

There is no Jacobian test in `tests/`. If `_jac` is added, also add a test comparing it with `utils.jac_fd` at a few states.

## Step 9 — Log-radius variant (optional, recommended)

Add a variant class with u = ln R as the state and `_postprocessing` returning exp(u). It needs its own JSON entry, keyed by the new class name.

As implemented in `RayleighPlessetLog`, u = ln(R / 1 m), so `_postprocessing` returns R in metres without needing R0. The solver then controls relative error in R, and R cannot become negative. It does **not** stop R from going below the hard core `a`.

## Step 10 — Full Keller–Miksis class (stretch)

The existing class already contains a first-order (R/c) correction. A full Keller–Miksis form adds (1 − U/c) and (1 + U/c) factors on the inertial and pressure terms. Below about 100 kHz, inertial collapses can reach Mach numbers where the first-order form stops being accurate. A `KellerMiksis` class, with the same SI parameters and its own JSON entry, shows how far the results depend on the model.

## Step 11 — Test reference data

`tests/test_generation.py::test_ensemble` compares 5-point Radau trajectories of every JSON system against `tests/test_data/all_trajectories.npy`. Rows follow the order of `get_attractor_list()`. Any change to the RP classes or their JSON entries changes their rows.

**Do not regenerate the whole file.** That integrates every system in the database. Recompute only the bubble rows:

```python
import numpy as np
import dysts.flows as flows
from dysts.base import get_attractor_list

path = "tests/test_data/all_trajectories.npy"
names, ref = get_attractor_list(), np.load(path, allow_pickle=True)
assert len(ref) == len(names), "a system was added or removed; insert its row with np.insert"
for name in ("RayleighPlesset", "RayleighPlessetLog"):
    sol = getattr(flows, name)().make_trajectory(5, method="Radau", resample=True)
    ref[names.index(name)] = sol[:, 0]
np.save(path, ref)
```

For day-to-day checks, run only the bubble tests: `python -m unittest tests.test_rayleigh_plesset`. `test_ensemble` (and so a bare `python -m unittest`) also integrates every system, so run it only before merging.

## Done when

- [x] RHS fix applied, and `python dysts/showRayleighPlesset.py` still runs (it needs `dvipng` for `usetex=True`; checked with usetex off)
- [x] Decision recorded on `_postprocessing`: **added**; the plotting script now passes `postprocess=False`
- [x] `_jac` rewritten (shared helpers `_rp_accel` and `_rp_accel_grad`), tested against relative-step central differences in `tests/test_rayleigh_plesset.py`, and wired into `make_trajectory` for Radau, BDF and LSODA. `Lorenz` is the only other system with a `_jac`, and it was checked too.
- [x] Log-radius variant agrees with the base class: max |ΔR| = 1.4e-7 R0 over 10 periods at 2 MHz
- [ ] Decision recorded on whether `KellerMiksis` is added
- [x] `all_trajectories.npy` regenerated once, in full (135 rows); `test_ensemble` passes. After the Phase 3 JSON changes, update only the bubble rows (Step 11).
- [ ] `python -m unittest` fully passes. The only remaining errors come from missing packages here: `pandas` (`test_datasets`) and `sdeint` (`test_trajectory_noise`).

## Notes

- Without numba, the analytic Jacobian made no measurable difference at the 2 MHz default (9.0 s vs 9.1 s for 2000 steps). Measure again at low frequency, where collapses are stiff.
