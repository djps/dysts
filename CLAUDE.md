# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A fork of `williamgilpin/dysts` (a database of ~130 chaotic ODE/DDE systems, NeurIPS 2021) that adds a `RayleighPlesset` bubble-dynamics model. The research goals (frequency-dependent bifurcation structure from 1 MHz down to ~70 kHz, then two interacting bubbles) are in `.claude/CLAUDE.md`. Phase plans are in `.claude/plans/plan-*.md`. They work with the **SI-unit** `RayleighPlesset` class in `dysts/flows.py` and its log-radius twin `RayleighPlessetLog` (state u = ln(R / 1 m)). Both use the shared `njit` helpers `_rp_accel` and `_rp_accel_grad`, so a physics change goes in one place. Each phase plan has a "Done when" checklist with the current status.

## Commands

Use the virtual environment at `~/venv3.14` (Python 3.14, numba 0.68, SciPy 1.15). Neither it nor the system `python` has pandas, sdeint, neurokit2 or nolds. The package no longer uses `pkg_resources`, and `blackmanharris` is imported from `scipy.signal.windows`, so it imports on current setuptools and SciPy.

```bash
source ~/venv3.14/bin/activate
pip install -I .                                   # install (deps: numpy, scipy, pandas; numba optional)
python -m unittest                                 # all tests; test_ensemble integrates EVERY system, so run only before merging
python -m unittest tests.test_generation.TestModels.test_trajectory   # single test
python -m unittest tests.test_rayleigh_plesset      # bubble RHS, Jacobians, log variant
python dysts/showRayleighPlesset.py                # RP phase/time plots (uses LaTeX via usetex=True)
```

`numba` is optional. Without it, `staticjit` falls back to plain Python, which is slower but behaves the same.

## Architecture

A system is defined in two places, and the two are joined by name:

1. **`dysts/flows.py`**: one `DynSys` subclass per system, with a `@staticjit` `_rhs(*state, t, *params)`. It can also have a `_jac` and a `_postprocessing` (which maps unbounded coordinates to bounded ones).
2. **`dysts/data/chaotic_attractors.json`**: metadata keyed by the **exact class name**. This holds `parameters`, `initial_conditions`, `dt`, `period`, Lyapunov spectrum and so on. Discrete maps follow the same pattern with `maps.py` and `discrete_maps.json`.

Key behaviours in `dysts/base.py` that are not obvious:

- `BaseDyn.__init__` loads the JSON entry and sets every parameter as an instance attribute. `DynSys.rhs` reads the parameters again with `getattr` on every call, in `sorted(params.keys())` order. That order is **case-sensitive ASCII**, so uppercase names come first. The `_rhs` positional arguments must match it exactly (RP: `P, Pa, R0, a, c, gamma, mu, omega, rho, sigma`).
- To change a parameter at run time, assign the attribute (`model.omega = ...`). `model.params` is **not** updated, but `rhs` uses the attribute. There is no custom `__setattr__`.
- `make_trajectory(n, method="Radau", resample=True, pts_per_period=100, ...)` has **no** `dt` or `init_cond` arguments. Set `model.dt` and `model.ic` instead.
  - With `resample=False`, the time grid is `arange(n) * model.dt`.
  - With `resample=True`, it spans `n / pts_per_period` multiples of `model.period`.
- If `model.ic` is 2-D, a trajectory is integrated for each row.
- `make_trajectory` passes `jac=` to `solve_ivp` when the system defines `_jac` and the method is Radau, BDF or LSODA (via `DynSys.jac`). Only `Lorenz` and the two RP classes define `_jac`.
- `utils.integrate_dyn` calls `solve_ivp` with hardcoded `rtol=1e-12, atol=1e-14` and uses `dt` only as `first_step`. Keep this in mind for badly scaled SI states, where R is about 1e-6.
- `make_trajectory_ensemble(n, subset=[names], **kwargs)` runs every system, or only a named subset, found through `get_attractor_list()` (the JSON keys).
- `analysis.py` holds `compute_timestep`, `find_lyapunov_exponents`, `kaplan_yorke_dimension` and related functions, which are used to fill the derived JSON fields. `lyap.py` and `equation_utils.py` contain further invariant and symbolic tools.

## Gotchas

- `tests/test_generation.py::test_ensemble` compares 5-step trajectories of **all** systems against `tests/test_data/all_trajectories.npy`. Adding, removing or reparametrising any JSON entry (including RP) changes this output, so the affected rows must be updated. Update only the bubble rows (Plan 2 Step 11); do not regenerate every system. The file was built without numba, so other systems' rows may drift past `allclose` under numba.
- `test_datasets` needs `pandas` and `test_trajectory_noise` needs `sdeint`.
- Precomputed datasets in `dysts/data/*_univariate*/*_multivariate*.json` predate the RP model and do not contain it. `load_trajectory` will fail with a KeyError for RP.
- `benchmarks/` reproduces the original paper. It has its own pinned `requirements.txt` (old darts/torch) and is not needed for the bubble work.
- `demos.ipynb` is about 16 MB. Avoid reading it whole.
