# CylPlasmaEuler1D

**1D Cylindrical Eulerian Two-Fluid Solver for Coupled Plasma-Neutral Dynamics**

A Wolfram Mathematica solver for modeling the radial expansion of a laser-produced plasma column into a surrounding neutral gas. The code solves the full two-temperature (Te ≠ Th) fluid equations with collisional ionization, charge-exchange drag, electron-ion thermalization, and flux-limited Spitzer thermal conduction — all in 1D cylindrical geometry.

---

## Physics

Two interpenetrating fluids — a plasma (ions + electrons) and a neutral gas — are coupled through collisional processes:

- **Collisional ionization** (Lotz empirical formula): neutrals → ions + electrons
- **Ion-neutral momentum transfer** (charge-exchange + elastic drag): velocity equilibration
- **Electron-ion thermalization** (Coulomb collisions, NRL formulary): temperature equilibration
- **Thermal conduction** (Spitzer, flux-limited): separate electron and ion channels

### Assumptions
- Quasineutrality: n_e = Z · n_i
- Electrons are inertialess and co-move with ions
- Ion temperature = Neutral temperature (T_i = T_n = T_h), justified by fast ion-neutral thermalization

### Conserved Variables (7 per cell)

| Variable | Description | Units |
|----------|-------------|-------|
| ρ_N | Neutral mass density | kg/cm³ |
| (ρu)_N | Neutral momentum density | kg·μm/(cm³·ns) |
| ε_hN | Neutral thermal energy density | eV/cm³ |
| ρ_P | Ion mass density | kg/cm³ |
| (ρu)_P | Ion momentum density | kg·μm/(cm³·ns) |
| ε_E | Electron thermal energy density | eV/cm³ |
| ε_hP | Ion thermal energy density | eV/cm³ |

---

## Numerical Method

| Component | Method |
|-----------|--------|
| Spatial discretization | Finite-volume, uniform cylindrical grid |
| Riemann solver | HLL (Harten-Lax-van Leer) |
| Time integration | SSP-RK3 (3rd-order strong stability preserving) |
| Operator splitting | Strang (2nd-order): S(dt/2) → A(dt) → S(dt/2) |
| Timestepping | Adaptive CFL: dt = 0.4 · dr / max(\|u\| + c_s) |
| Conduction | Implicit tridiagonal (Thomas algorithm), flux-limited |
| Drag | Implicit analytical (unconditionally stable) |
| Thermalization | Implicit 2×2 analytical |
| Ionization | Explicit (guarded for efficiency) |

### Key Numerical Features

- **Advective-pressure flux splitting** for cylindrical momentum: eliminates the geometric source term p/r that causes well-balancing errors near the axis
- **Flux-limited Spitzer conduction**: harmonic-mean limiter caps heat flux at the free-streaming limit q_fs = f · n · T · v_th, preventing unphysical heat transport ahead of the shock
- **pdV correction in drag step**: compensates the Strang splitting error that causes spurious heating in drag-deceleration zones
- **Compiled kernels**: HLL flux loops, floor enforcement, and pdV correction compiled to WVM/C bytecode
- **Inlined physics**: all rate coefficients precomputed as constants — no external function calls during the time loop

---

## Requirements

- **Wolfram Mathematica** 12.0 or later
- Global constants must be defined before calling the solver:
  ```mathematica
  eCharge = 1.602*^-19;  (* J/eV *)
  me = 9.109*^-31;       (* electron mass, kg *)
  ```

---

## Quick Start

```mathematica
(* Load the solver and analysis modules *)
Get["path/to/intTimeSol.wl"]
Get["path/to/analyzeSolution.wl"]

(* Define global constants *)
eCharge = 1.602*^-19;
me = 9.109*^-31;
mikg = 1.673*^-27;  (* proton mass *)

(* Set up initial conditions *)
(* interpNiRInit: interpolating function for ni(r)/n0 *)
(* interpEpsEInit: interpolating function for electron energy density epsE(r) [eV/cm^3] *)

(* Run the solver *)
sol = intTimeSol[
  0.1,           (* dr [um] *)
  60,            (* rMax [um] — auto-extended if needed *)
  mikg,          (* ion mass [kg] *)
  2.4*10^18,     (* n0 [cm^-3] *)
  interpNiRInit, (* normalized ion density profile *)
  interpEpsEInit,(* electron energy density profile *)
  0.025,         (* Ti [eV] *)
  0.025,         (* Tn [eV] *)
  13.6,          (* ionization energy [eV] *)
  mikg/me,       (* ion-to-electron mass ratio *)
  1,             (* charge state Z *)
  5.0            (* simulation time [ns] *)
];

(* Analyze and visualize *)
analyzeSolution[sol, 993, 99.3, mikg, 2.4*10^18, 5.0, {1/3, 2/3, 1}]
```

---

## Solver API

### `intTimeSol`

```
intTimeSol[drBYmu, rMaxBYmu, mnkg, n0, interpNiRInit,
           interpEpsEInit, TieV, TneV, IEeV, mui, Z, tFinns,
           fFluxLim, nSnap,
           enableIoniz, enableDrag, enableTherm, enableECond, enableICond]
```

**Required arguments:**

| Argument | Type | Description |
|----------|------|-------------|
| `drBYmu` | Real | Grid spacing [μm] |
| `rMaxBYmu` | Real | Minimum domain radius [μm] |
| `mnkg` | Real | Ion/neutral mass [kg] |
| `n0` | Real | Total number density [cm⁻³] |
| `interpNiRInit` | InterpolatingFunction | Normalized ion fraction n_i(r)/n_0 |
| `interpEpsEInit` | InterpolatingFunction | Electron energy density ε_E(r) [eV/cm³] |
| `TieV` | Real | Initial ion temperature [eV] |
| `TneV` | Real | Initial neutral temperature [eV] |
| `IEeV` | Real | Ionization energy [eV] |
| `mui` | Real | Ion-to-electron mass ratio m_i/m_e |
| `Z` | Integer | Ion charge state |
| `tFinns` | Real | Simulation end time [ns] |

**Optional arguments (with defaults):**

| Argument | Default | Description |
|----------|---------|-------------|
| `fFluxLim` | 0.06 | Flux limiter fraction (0.03–0.1) |
| `nSnap` | 500 | Number of stored time snapshots |
| `enableIoniz` | True | Enable/disable ionization |
| `enableDrag` | True | Enable/disable drag |
| `enableTherm` | True | Enable/disable thermalization |
| `enableECond` | True | Enable/disable electron conduction |
| `enableICond` | True | Enable/disable ion conduction |

**Returns:** `{rules}` — a list of interpolating function replacement rules for all 7 conserved variables at each grid cell, accessible as `(variableName[cellIndex] /. rules)[time]`.

### `analyzeSolution`

```
analyzeSolution[sol, Nr, rMaxBYmu, mnkg, n0, tFinal, tPlotFracs,
                nMovieFrames, rPlotMax, nMax, TMin, TMax, exportDir]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `nMovieFrames` | 100 | Number of frames in each movie |
| `rPlotMax` | Automatic | Max radius to display [μm] |
| `nMax` | Automatic | Max density for y-axis [cm⁻³] |
| `TMin` | 0.01 | Min temperature for y-axis [eV] |
| `TMax` | Automatic | Max temperature for y-axis [eV] |
| `exportDir` | "" | Directory for MP4 export (default: NotebookDirectory) |

**Produces:**
- Printed diagnostics (mass conservation, field ranges)
- Grid of checkpoint plots (temperature, density, velocity)
- Two MP4 movies (temperature evolution, density evolution)
- Interactive `ListAnimate` widgets

---

## File Structure

```
CylPlasmaEuler1D/
├── README.md                  # This file
├── LICENSE                    # MIT License
├── .gitignore                 # Mathematica-specific ignores
├── intTimeSol.wl              # Main solver module
├── analyzeSolution.wl         # Analysis and visualization module
├── docs/
│   ├── equations.md           # Detailed equation reference
│   └── numerical_notes.md     # Numerical method details and known issues
└── examples/
    └── hydrogen_expansion.wl  # Example: H plasma column expansion
```

---

## Units

All quantities use a consistent unit system:

| Quantity | Unit |
|----------|------|
| Length | μm |
| Time | ns |
| Velocity | μm/ns (= km/s) |
| Mass density | kg/cm³ |
| Number density | cm⁻³ |
| Energy density | eV/cm³ |
| Temperature | eV |
| Pressure | eV/cm³ (multiply by e·10⁻⁶ for kg·μm²/(ns²·cm³)) |

---

## Citation

If you use this code in published work, please cite:

```
CylPlasmaEuler1D: 1D Cylindrical Two-Fluid Plasma-Neutral Solver
https://github.com/ritwik-sain/Cyl1D-Euler-Solver.git
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
