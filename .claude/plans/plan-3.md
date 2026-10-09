# Rayleigh–Plesset in dysts — Phase 3: Metadata

Oct 6, 2026 · @David Sinden

## Goal

This phase corrects the existing `RayleighPlesset` entry in `dysts/data/chaotic_attractors.json` (it is the first entry in the file), so it describes one canonical chaotic parameter set.

The JSON key must equal the class name, because `get_attractor_list()` reads the keys. The frequency sweep in Phase 4 does not add database entries. It changes `omega` on an instance at run time.

## Step 12 — Choose a published chaotic reference point

CONTRIBUTING.md asks for parameter values from published work. The current `citation` reads "diaz de la rosa or MacDonald and Gomatam (2006)", which still needs to be settled. Candidates:

- MacDonald & Gomatam (2006), *Chaotic dynamics of microbubbles in ultrasonic fields*, doi:10.1243/095440606X79596 (the source currently cited)
- Lauterborn & Parlitz, *Methods of chaos physics and their application to acoustics*, JASA 1988
- Parlitz, Englisch, Scheffczyk & Lauterborn, *Bifurcation structure of bubble oscillators*, JASA 1990

Use the chosen paper's (R0, f, P, gas, liquid) values directly in SI. Pick the paper whose model matches this class, which is RP with a hard core and a first-order c term.

## Step 13 — Fields to fix by hand

| Field | Current | Should be |
| --- | --- | --- |
| `parameters` | SI values, ω = 2π·2 MHz | values from Step 12 |
| `initial_conditions` | `[2e-6, 0, 0]` (rest state) | a point on the attractor after transients, written as [R, U, Ω] with Ω taken mod 2π |
| `nonautonomous` | false | true (as for `Duffing`) |
| `unbounded_indices` | `[]` | `[2]` (the phase) |
| `bifurcation_parameter` | null | `"omega"` |
| `description` | a note about dt | a one-line physical description |
| `citation` | two candidate sources | the single chosen reference |
| `embedding_dimension`, `delay`, `hamiltonian` | 3, false, false | unchanged |

`dt` and `period` are in seconds. `period` should equal the drive period 2π/ω (currently 5e-7 s for 2 MHz).

## Step 14 — Derived fields from repo utilities

The current derived values have to be recomputed. They predate the Phase 1 fix, and a largest exponent of 0.089 s⁻¹ is implausibly small for an SI system with MHz forcing.

| Field | Utility |
| --- | --- |
| `dt`, `period` | `analysis.compute_timestep(model)` |
| `lyapunov_spectrum_estimated`, `maximum_lyapunov_estimated` | `analysis.find_lyapunov_exponents(model, traj_length, pts_per_period=...)` (expect one exponent ≈ 0 from the phase) |
| `kaplan_yorke_dimension` | `analysis.kaplan_yorke_dimension` on that spectrum |
| `pesin_entropy` | sum of the positive exponents |
| `correlation_dimension` | `lyap.corr_dim` on a trajectory |
| `multiscale_entropy` | `analysis.mse_mv` on a trajectory (requires `neurokit2`) |

`find_lyapunov_exponents` advances the tangent space with a first-order Euler step (I + J·dt) and a finite-difference Jacobian. It needs a large `pts_per_period` to resolve collapses. Check convergence by doubling it.

## Step 15 — Test and submit

1. Install: `pip install -I .` (or `pip install -e .`).
2. Update only the bubble rows of `tests/test_data/all_trajectories.npy` (Phase 2 Step 11), then run `python -m unittest tests.test_rayleigh_plesset`. Run the full suite, which integrates every system, only before merging.
3. Optional: if contributing upstream to `williamgilpin/dysts`, post the system on that repo's issue #1, as CONTRIBUTING.md asks.

## Done when

- [ ] Reference point chosen and entered in SI
- [ ] Hand-filled fields corrected
- [ ] Derived fields recomputed after the Phase 1 fix
- [ ] `tests.test_rayleigh_plesset` passes (and the full suite before merging)
