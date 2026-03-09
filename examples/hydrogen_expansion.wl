(* ============================================================================
   Example: Hydrogen Plasma Column Expansion
   
   A hot hydrogen plasma column (r ~ 5 um, Te ~ 10 eV) expands into
   a cold neutral hydrogen background (n0 = 2.4e18 cm^-3, T = 0.025 eV).
   
   This example demonstrates the basic workflow:
     1. Define global constants
     2. Set up initial condition profiles
     3. Run the solver
     4. Analyze and visualize results
   ============================================================================ *)

(* Load the solver and analysis modules — adjust paths as needed *)
Get[NotebookDirectory[] <> "intTimeSol.wl"]
Get[NotebookDirectory[] <> "analyzeSolution.wl"]

(* ================================================================ *)
(* GLOBAL CONSTANTS                                                  *)
(* ================================================================ *)
eCharge = 1.602*^-19;   (* J/eV *)
me = 9.109*^-31;        (* electron mass [kg] *)
mikg = 1.673*^-27;      (* proton mass [kg] *)

(* ================================================================ *)
(* PROBLEM PARAMETERS                                                *)
(* ================================================================ *)
n0 = 2.4*^18;           (* total number density [cm^-3] *)
TieV = 0.025;           (* ambient ion/neutral temperature [eV] *)
TneV = 0.025;           (* ambient neutral temperature [eV] *)
IEeV = 13.6;            (* hydrogen ionization energy [eV] *)
Z = 1;                  (* charge state *)
mui = mikg/me;          (* ion-to-electron mass ratio *)

(* Grid and time *)
dr = 0.1;               (* cell size [um] *)
rMax = 60.;             (* minimum domain [um] — solver auto-extends *)
tFin = 5.0;             (* simulation time [ns] *)

(* ================================================================ *)
(* INITIAL CONDITIONS                                                *)
(* Create interpolating functions for:                               *)
(*   interpNiRInit[r]: normalized ion fraction ni(r)/n0              *)
(*   interpEpsEInit[r]: electron energy density [eV/cm^3]            *)
(*                                                                   *)
(* Replace this section with your own initial conditions.            *)
(* The profiles below are a simple Gaussian example.                 *)
(* ================================================================ *)

(* Example: Gaussian plasma column *)
rColumn = 5.0;   (* column radius [um] *)
Te0 = 10.0;      (* peak electron temperature [eV] *)

(* Normalized ion fraction: Gaussian, clamped at 1e-9 *)
niProfile[r_] := Max[Exp[-r^2/(2*rColumn^2)], 1*^-9];
interpNiRInit = Interpolation[
  Table[{r, niProfile[r]}, {r, 0, 200, 0.05}],
  InterpolationOrder -> 3];

(* Electron energy density: epsE = (3/2) ni * Te *)
epsEprofile[r_] := 1.5*(n0*niProfile[r])*Te0*Exp[-r^2/(2*rColumn^2)];
interpEpsEInit = Interpolation[
  Table[{r, epsEprofile[r]}, {r, 0, 200, 0.05}],
  InterpolationOrder -> 3];

(* ================================================================ *)
(* RUN THE SOLVER                                                    *)
(* ================================================================ *)
Print["Starting solver..."];
sol = intTimeSol[dr, rMax, mikg, n0, interpNiRInit,
  interpEpsEInit, TieV, TneV, IEeV, mui, Z, tFin];

(* ================================================================ *)
(* ANALYZE AND VISUALIZE                                             *)
(* ================================================================ *)

(* Get actual grid size from solver output *)
solRules = First[sol];
Nr = Length[solRules]/7;  (* 7 rules per cell *)

analyzeSolution[sol, Nr, Nr*dr, mikg, n0, tFin,
  {1/3, 2/3, 1},    (* checkpoint times as fractions of tFin *)
  100,               (* movie frames *)
  50,                (* rPlotMax [um] *)
  Automatic,         (* nMax *)
  0.02,              (* TMin [eV] *)
  12,                (* TMax [eV] *)
  NotebookDirectory[]  (* export movies here *)
]
