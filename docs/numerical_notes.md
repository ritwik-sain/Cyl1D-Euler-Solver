# Numerical Notes and Known Issues

## Numerical Thresholds

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `rhoFloor` | 10⁻²⁰ kg/cm³ | Absolute minimum density |
| `rhoVac` | m_i · n₀ · 10⁻⁹ | Vacuum threshold (matches IC density clamp); momentum zeroed below this |
| `TeFloor` | T_ieV (input) | Electron temperature floor = ambient temperature |
| `ThFloor` | T_ieV (input) | Heavy temperature floor = ambient temperature |
| `TeIzGuard` | I_E / 27 | Ionization skipped below (rate < 10⁻¹² of peak) |
| Momentum guard | \|mom\| < 10⁻²⁵ | Prevents spurious velocity from roundoff ÷ small density |
| Drag guard | \|u_P - u_N\| < 10⁻²⁰ | No relative motion |
| Thermalization guard | \|T_e - T_h\|/T_e < 10⁻⁶ | Already equilibrated |
| Conductivity guard | T < T_floor | κ negligible at floor temperature |
| CFL number | 0.4 | Stability limit for SSP-RK3 |

Temperature floors are enforced **density-locally**: ε ≥ (3/2)(ρ_local/m_i)·T_floor.
This means the energy floor scales with density — it's a temperature floor, not an energy floor.

## Advective-Pressure Flux Splitting

The standard 1D cylindrical momentum equation:

∂(ρu)/∂t + (1/r)∂(r(ρu² + p))/∂r = p/r

has a geometric source term p/r on the right. This causes well-balancing problems:
in uniform pressure, the LHS and RHS should cancel exactly, but finite differencing
introduces O(dr²) errors that create spurious velocities near the axis.

The advective-pressure split reformulates this as:

∂(ρu)/∂t + (1/r)∂(r · ρu²)/∂r + ∂p/∂r = 0

The advective flux (ρu²) gets the cylindrical divergence (1/r)∂(r·)/∂r,
while the pressure gradient is applied as a flat ∂p/∂r. No geometric source term
is needed, and uniform pressure gives exactly zero acceleration.

## Strang Splitting and the pdV Correction

The operator-split structure is: Source(dt/2) → RK3_Advection(dt) → Source(dt/2).

A known issue with Strang splitting arises when the source step modifies the
velocity field (drag) without updating the energy self-consistently. The drag
step decelerates the plasma, creating localized compression (negative ∂u/∂r).
In the subsequent advection step, this compression heats the gas through pdV work.
But the pressure feedback (which would resist the compression) doesn't act until
the next source step — a one-timestep delay.

**The pdV correction** addresses this by computing the velocity divergence change
ΔdivU = divU_post_drag - divU_pre_drag and applying the corresponding pdV work
immediately within the drag step:

ε *= (1 - (2/3) · ΔdivU · dt)

This ensures that when drag compresses the plasma, the energy response is
instantaneous rather than delayed by one timestep.

## Flux-Limited Conduction

Classical Spitzer conduction κ ∝ T^{5/2} can transport heat faster than the
free-streaming limit in regions with steep temperature gradients (e.g., shock
fronts). The harmonic-mean flux limiter:

κ_eff = κ_Sp · q_fs / (κ_Sp · |∇T| + q_fs)

smoothly transitions between:
- Classical regime (gentle gradient): κ_eff ≈ κ_Sp
- Free-streaming regime (steep gradient): κ_eff ≈ q_fs / |∇T|

The free-streaming flux is q_fs = f · n · T · v_th where f = 0.06 is standard
for laser-plasma interactions.

## Known Artifacts

### Temperature spike at the plasma front

At the leading edge of the expanding plasma, where ρ_P drops steeply to vacuum
over 2-3 cells, the HLL solver numerically diffuses ε_E slightly more than ρ_P.
Since T_e = (2/3)·ε_E·m_i/ρ_P, this creates a localized temperature spike in
cells where ρ_P has dropped but ε_E hasn't.

**Impact:** Negligible — the total thermal energy in these cells is ε_E × volume ≈ 0.
The spike does not propagate inward or affect the bulk solution.

**Mitigation:** When plotting T_e, mask values in cells where ρ_P < 10⁻⁷ · n₀ · m_i.
The `analyzeSolution` module does this automatically.

### Temperature plateau between neutral and plasma fronts

In the region between the neutral shock front (r_s) and the plasma velocity front (r_p),
drag-induced compression of the plasma velocity field creates physical heating via pdV work.
With the pdV correction active, the plateau magnitude is reduced to a physically reasonable
level consistent with the local compression ratio.

## Adaptive Timestepping

The CFL condition dt = CFL · dr / max(|u| + c_s) is recomputed every step.
For a typical plasma expansion:
- Initial: T_e ~ 10 eV → c_s ~ 40 μm/ns → dt ~ 0.001 ns
- Late time: T_e ~ 0.1 eV → c_s ~ 4 μm/ns → dt ~ 0.01 ns

This gives 3-5× fewer total steps compared to fixed timestepping with the initial dt.

## Performance Notes

For a typical run (Nr ≈ 1000, t_final = 5 ns):
- ~1500-2500 adaptive timesteps
- ~250-400 ms per step (WVM compilation on typical hardware)
- Total wall-clock: ~8-15 minutes

The dominant cost is in the source step (Steps C-E: thermalization + conduction),
which involves Nr interpreted loop iterations with inlined arithmetic.
The advection step (compiled HLL + RK3) is comparatively fast.

For parameter scans, use `ParallelTable` across independent runs for near-linear
speedup on multi-core systems.
