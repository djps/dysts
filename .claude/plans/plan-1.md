# Rayleigh–Plesset in dysts — Phase 1: Formulation

Oct 6, 2026 · @David Sinden

## Goal

This phase pins down the equations that `RayleighPlesset` in `dysts/flows.py` already implements, fixes the one known error in them, and checks the default parameters. The class works in **SI units**. A single parameter, `omega` (rad/s), carries the whole 1 MHz to 25 kHz sweep.

This repo is the fork `djps/dysts`. The class was added in commit `68c811c`. Like `Duffing`, it turns periodic forcing into an autonomous system by adding a phase variable (`dotOmega = omega`).

## Step 1 — Model as implemented

The state is (R, U = Ṙ, Ω), where Ω is the driving phase. The equation is Rayleigh–Plesset with a van der Waals hard core, polytropic gas, surface tension and viscosity. It also has a first-order liquid-compressibility term of the Herring/Keller type, so it is **not** plain RP:

```latex
R\dot U + \tfrac{3}{2}U^2 = \frac{1}{\rho}\left[p_g - p_\infty(t) - \frac{2\sigma}{R} - \frac{4\rho\mu U}{R} + \frac{R}{c}\left(\dot p_g - \dot p_\infty\right)\right]

p_g = \left(P_a + \frac{2\sigma}{R_0}\right)\left(\frac{R_0^3 - a^3}{R^3 - a^3}\right)^{\gamma}, \qquad
\dot p_g = -\frac{3\gamma\, p_g R^2 U}{R^3 - a^3}, \qquad
p_\infty(t) = P_a - P\sin\Omega, \qquad \dot\Omega = \omega
```

The driving starts with rarefaction: for small Ω > 0, p∞ falls below ambient. The time derivatives of the surface-tension and viscous terms are left out of the (R/c) correction. That is a usual approximation, but it should be written up.

## Step 2 — Parameters

dysts passes parameters to `_rhs` positionally in `sorted()` order. That order is case-sensitive, so the uppercase names come first: `P, Pa, R0, a, c, gamma, mu, omega, rho, sigma`.

| Parameter | Meaning | Default in `chaotic_attractors.json` |
| --- | --- | --- |
| `P` | driving pressure amplitude (Pa) | 1.0e6 |
| `Pa` | ambient pressure (Pa) | 101325 |
| `R0` | equilibrium radius (m) | 2e-6 |
| `a` | van der Waals hard-core radius (m) | 0.1411e-6 |
| `c` | speed of sound in liquid (m/s) | 1480 |
| `gamma` | polytropic exponent | 1.66666 (monatomic, adiabatic) |
| `mu` | **kinematic** viscosity (m²/s); the code multiplies by ρ | 1.002e-6 |
| `omega` | angular driving frequency (rad/s) | 12.566e6, which is **2 MHz**, not 1 MHz |
| `rho` | liquid density (kg/m³) | 1000 |
| `sigma` | surface tension (N/m) | 0.073 |

Note that `a` is a hard-core radius. It is not a pressure ratio as in the earlier nondimensional draft.

## Step 3 — Fix the driving-derivative term

Line 42 of `dysts/flows.py` computes `(R / c) * (P * dWave + dPg)` with `dWave = -cos(Omega)`. This is wrong in two ways:

- **Units:** `P * dWave` is a pressure, but it needs to be a pressure rate. The factor `omega` is missing.
- **Sign:** the correct term is −ṗ∞ = −(d/dt)(P·Wave) = +Pω cos Ω, which equals `-P * omega * dWave`. The code has `+P * dWave`.

The corrected line is:

```python
(R / c) * (-P * omega * dWave + dPg)
```

At 2 MHz the error is a factor of about 10⁷ on that term, so all existing trajectories and derived JSON fields change once it is fixed.

## Step 4 — Scaling: stay in SI

`utils.integrate_dyn` hardcodes `rtol=1e-12` and `atol=1e-14`. With R of order 1e-6 m (and still above `a` ≈ 1.4e-7 m at collapse), atol stays about 7 orders of magnitude below R. SI units are therefore workable. In SI the frequency sweep changes only `omega`, so nondimensionalising is not needed for the sweep.

A dimensionless variant (radius scaled by R0, time by R0·√(ρ/p0)) is still an option if conditioning turns out to be a problem. If it is added, it must be a separate class with its own JSON key.

## Step 5 — Linear resonance check

Ignoring the hard core and c, the natural frequency is:

```latex
\omega_0^2 = \frac{1}{\rho R_0^2}\left[3\gamma\left(P_a + \frac{2\sigma}{R_0}\right) - \frac{2\sigma}{R_0}\right]
```

For the defaults this gives f0 ≈ 1.7 MHz (γ = 1), 2.0 MHz (γ = 1.4) and 2.25 MHz (γ = 5/3). The current 2 MHz default therefore sits near the main resonance. 1 MHz is near f0/2, and 25 kHz is about 80× below resonance.

## Done when

- [x] Step 3 sign and ω fix applied in `dysts/flows.py`
- [x] Decision recorded on the default drive frequency: **move to 1 MHz** (JSON `omega`, `period`, `dt` updated in Phase 3)
- [x] Choice of `gamma` recorded: **fixed** in the model; sweeps are run and reported at both γ = 1 and γ = 5/3 (Phase 5)
- [x] Linear resonance at small `P` checked against Step 5 (see below)

## Resonance check result

P = 1 kPa, 33 frequencies from 0.8 to 4 MHz, amplitude over the last 15 of 50 periods, with the peak refined by a parabola fit:

| γ | Step 5 formula | Simulated peak | Peak amplitude |
| --- | --- | --- | --- |
| 5/3 | 2.249 MHz | 2.239 MHz | 0.0126 R0 |
| 1 | 1.688 MHz | 1.688 MHz | 0.0197 R0 |

The small downward shift at γ = 5/3 is what damping is expected to produce.
