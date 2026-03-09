# Equations Reference

## Advection Equations (Hyperbolic System)

All equations are in 1D cylindrical coordinates (r), with velocity u in the radial direction.

### Neutral Fluid

The neutral fluid has pressure p_hN = (2/3) ε_hN [mechanical units].

**Continuity:**

∂ρ_N/∂t + (1/r) ∂(r ρ_N u_N)/∂r = 0

**Momentum** (advective-pressure split form):

∂(ρ_N u_N)/∂t + (1/r) ∂(r [ρ_N u_N² + p_hN - p_hN])/∂r + ∂p_hN/∂r = 0

The term in brackets is the advective flux (ρu² only), diverged in cylindrical form.
The pressure gradient ∂p/∂r is applied as a flat (non-cylindrical) gradient.
This avoids the geometric source term p/r that appears in the standard formulation.

**Thermal energy** (with pdV work):

∂ε_hN/∂t + (1/r) ∂(r ε_hN u_N)/∂r = -(2/3) ε_hN ∇·u_N

where ∇·u = (1/r) ∂(r u)/∂r. The energy is advected with mass-flux upwinding
(specific energy ε/ρ carried at upwind value), and the pdV term accounts for
compressive/expansive work.

### Plasma Fluid

The plasma fluid has total pressure p_tot = p_e + p_hP = (2/3)(ε_E + ε_hP) [mechanical].

**Continuity:**

∂ρ_P/∂t + (1/r) ∂(r ρ_P u_P)/∂r = 0

**Momentum:**

∂(ρ_P u_P)/∂t + (1/r) ∂(r ρ_P u_P²)/∂r + ∂p_tot/∂r = 0

**Electron thermal energy:**

∂ε_E/∂t + (1/r) ∂(r ε_E u_P)/∂r = -(2/3) ε_E ∇·u_P

**Ion thermal energy:**

∂ε_hP/∂t + (1/r) ∂(r ε_hP u_P)/∂r = -(2/3) ε_hP ∇·u_P

Each energy species does pdV work against its own partial pressure,
using the plasma velocity divergence.

---

## Source Terms

Applied via Strang splitting: S(dt/2) → Advection(dt) → S(dt/2).

### A. Collisional Ionization

**Rate coefficient** (Lotz empirical formula):

⟨σv⟩_iz = (3.015 × 10⁻⁷ / (√T_e · I_E)) · Γ(0, I_E/T_e)  [cm³/s]

where Γ(0, x) = ∫_x^∞ (e^{-t}/t) dt is the upper incomplete gamma function (= E₁(x)).

**Ionization rate:**

S_n = ⟨σv⟩_iz · n_i · n_n  [cm⁻³ ns⁻¹]

(with the 10⁻⁹ factor for ns units)

**Updates:**
- ρ_P += m_i · S_n · dt,  ρ_N -= m_i · S_n · dt  (mass transfer)
- (ρu)_P += u_N · m_i · S_n · dt  (new ions carry neutral momentum)
- ε_E -= I_E · S_n · dt  (ionization potential extracted from electrons)
- ε_hP += (3/2) S_n · T_h · dt  (new ions enter at ambient T_h)
- ε_hN -= (3/2) S_n · T_h · dt  (thermal energy leaves with neutrals)
- KE dissipation: (1/2)(u_P - u_N)² · m_i · S_n → heavy thermal (split by mass fraction)

**Guard:** Skipped when T_e < I_E/27 (rate < 10⁻¹² of peak).

### B. Ion-Neutral Drag

**Momentum transfer rate** (charge-exchange + elastic, equal masses):

v_RMS = √((u_P - u_N)² + 8 k_B (T_i + T_n) / (π m_i))  [μm/ns]

K_mt = 2.13 × 10⁻⁹ · v_RMS^{0.75}  [cm³/s]

**Implicit velocity update** (unconditionally stable):

Δu^{new} = Δu^{old} / (1 + (ν_P + ν_N) dt)

where ν_P = (1/2)(n_N) K_mt and ν_N = (1/2)(n_P) K_mt.

**KE dissipation:** KE_before - KE_after → heavy thermal energy, split by mass fraction f_P = ρ_P/(ρ_P + ρ_N).

**pdV correction:** After drag modifies velocities, the change in velocity divergence
ΔdivU = divU_post - divU_pre is computed, and the corresponding pdV work is applied:

ε *= (1 - (2/3) · ΔdivU · dt)

for each energy species (ε_E and ε_hP use plasma divU; ε_hN uses neutral divU).

### C. Electron-Ion Thermalization

**Exchange rate** (NRL formulary):

α = (3 m_e / m_i) · n_e² · ln(Λ) / (3.44 × 10⁵ · T_e^{3/2})  [ns⁻¹]

**Implicit 2×2 solve** conserving A_e · T_e + A_h · T_h:

T_e^{new} = (A_e A_h T_e + β(A_e T_e + A_h T_h)) / (A_e A_h + β(A_e + A_h))

T_h^{new} = (A_e A_h T_h + β(A_e T_e + A_h T_h)) / (A_e A_h + β(A_e + A_h))

where A_e = (3/2)n_e, A_h = (3/2)n_h, β = α · dt.

After updating T_h, the total heavy energy is split between ions and neutrals
by mass fraction: ε_hP = ε_h_tot · f_P, ε_hN = ε_h_tot · (1 - f_P).

### D. Electron Spitzer Conduction

**Spitzer conductivity:**

κ_e = 3.2 · (e/m_e) · 3.44 × 10⁵ · T_e^{5/2} / ln(Λ)  [μm² cm⁻³ ns⁻¹ eV⁻¹]

**Flux limiter** (harmonic mean):

κ_eff = κ_Sp · q_fs / (κ_Sp · |∂T/∂r| + q_fs)

where q_fs = f · n_e · T_e · v_te is the free-streaming heat flux, v_te = √(eT_e/m_e),
and f = 0.06 is the flux limiter fraction.

**Implicit tridiagonal solve** (Thomas algorithm):

(3/2)(n_e/dt) · T_e^{new} = (3/2)(n_e/dt) · T_e^{old} + (1/r) ∂/∂r(r κ_eff ∂T_e^{new}/∂r)

### E. Ion Spitzer Conduction

Same structure as electron conduction, with:

κ_ion = √(m_e/m_i) · κ_e(T_h, ln Λ)

Heat capacity uses total heavy density: C = (3/2) n_total / dt.
After solving for T_h^{new}, distribute total heavy energy by mass fraction.

---

## Coulomb Logarithm (NRL Formulary)

Three regimes:

| Condition | Formula |
|-----------|---------|
| T_i/μ ≤ T_e ≤ 10Z² | ln Λ = 23 - ln(√n_e · Z · T_e^{-3/2}) |
| T_i/μ ≤ 10Z² < T_e | ln Λ = 24 - ln(√n_e · T_e⁻¹) |
| T_e < T_i/μ | ln Λ = 30 - ln(√n_i · T_i^{-3/2} · Z²/μ) |

Floored at ln Λ = 2 for stability.
