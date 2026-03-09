(* ============================================================================
   intTimeSol: 1D Cylindrical Eulerian Two-Fluid Plasma-Neutral Solver
   ============================================================================
   
   PHYSICS MODEL
   =============
   Two interpenetrating fluids (plasma ions + electrons, and neutrals) in
   1D cylindrical geometry (r), coupled by collisional ionization, charge-
   exchange drag, and electron-ion thermalization.
   
   Assumptions:
     - Ion temperature = Neutral temperature (Ti = Tn = Th), justified by
       fast ion-neutral thermalization relative to the expansion timescale
     - Quasineutrality: ne = Z * ni
     - Electrons are inertialess and co-move with ions
   
   Conserved variables (7 per cell):
     rhoN  [kg/cm^3]          neutral mass density
     momN  [kg um/(cm^3 ns)]  neutral momentum density (radial)
     epsHN [eV/cm^3]          neutral thermal energy = (3/2)(rhoN/mi)Th
     rhoP  [kg/cm^3]          ion mass density
     momP  [kg um/(cm^3 ns)]  ion momentum density (radial)
     epsE  [eV/cm^3]          electron thermal energy = (3/2)(rhoP/mi)Te
     epsHP [eV/cm^3]          ion thermal energy = (3/2)(rhoP/mi)Th
   
   NUMERICAL METHOD
   ================
   Advection (hyperbolic):
     - HLL approximate Riemann solver on cell faces
     - Advective-pressure flux splitting for cylindrical momentum
       (eliminates geometric source term p/r and well-balancing errors)
     - SSP-RK3 (3rd-order strong stability preserving Runge-Kutta)
     - Mass-flux upwinding for energy advection
     - Pressure:
       * Ion momentum: pe + phP (electron + ion heavy pressure)
       * Neutral momentum: phN (neutral heavy pressure)
     - pdV work: each energy species does work against its own pressure
       using its own velocity divergence
   
   Source terms (operator-split, Strang: S(dt/2)-A(dt)-S(dt/2)):
     A. Collisional ionization (explicit, Lotz empirical formula)
     B. Ion-neutral momentum transfer drag (implicit analytical)
        + pdV correction for Strang splitting consistency
     C. Electron-ion temperature equilibration (implicit 2x2 analytical)
     D. Electron Spitzer thermal conduction (implicit tridiagonal,
        harmonic-mean flux-limited)
     E. Ion Spitzer thermal conduction (implicit tridiagonal,
        flux-limited, kappa_ion = sqrt(me/mi) * kappa_e(Th))
   
   OPTIMIZATIONS
   =============
   - Compiled (WVM/C) RHS flux computation and floor enforcement
   - All physics formulae (Spitzer kappa, Coulomb log, Lotz ionization,
     momentum transfer) inlined with precomputed constants
   - Compiled pdV correction for drag step
   - Adaptive CFL timestepping: dt = CFL * dr / max(|u|+cs)
   - Reduced floor enforcement (only after mass-transfer and at step end)
   
   CALLING CONVENTION
   ==================
   sol = intTimeSol[drBYmu, rMaxBYmu, mnkg, n0BYcm3, interpNiRInit,
                    interpEpsEInit, TieV, TneV, IEeV, mui, Z, tFinns]
   
   Required arguments:
     drBYmu         grid spacing [um]
     rMaxBYmu       minimum domain size [um] (auto-extended if needed)
     mnkg           ion/neutral mass [kg] (e.g. mH = 1.67e-27)
     n0BYcm3        total number density [cm^-3]
     interpNiRInit  interpolating function: ni(r)/n0 (normalized ion frac)
     interpEpsEInit interpolating function: epsE(r) [eV/cm^3]
     TieV           initial ion temperature [eV]
     TneV           initial neutral temperature [eV]
     IEeV           ionization energy [eV] (13.6 for hydrogen)
     mui            ion-to-electron mass ratio (mi/me)
     Z              ion charge state
     tFinns         simulation end time [ns]
   
   Optional arguments (with defaults):
     fFluxLim  (0.06)  flux limiter fraction (range 0.03-0.1)
     nSnap     (500)   number of uniformly-spaced time snapshots stored
     enableIoniz/Drag/Therm/ECond/ICond (True) source term switches
   
   Returns: {rules} where rules is a list of interpolating functions
     rhoNBYkgcm3[i][t], momN...[i][t], epsHN...[i][t],
     rhoP...[i][t], momP...[i][t], epsE...[i][t], epsHP...[i][t]
   
   UNITS
   =====
   r [um], t [ns], u [um/ns = km/s], rho [kg/cm^3],
   mom [kg um/(cm^3 ns)], energy density [eV/cm^3], T [eV]
   
   Global constants required: eCharge [J/eV = 1.602e-19], me [kg = 9.109e-31]
   ============================================================================ *)

intTimeSol[drBYmu_, rMaxBYmu_, mnkg_, n0BYcm\[Bullet]3_, interpNiRInit_,
  interpEpsEInit_, TieV_, TneV_, IEeV_, mui_, Z_, tFinns_,
  fFluxLim_:0.06, nSnap_:500,
  enableIoniz_:True, enableDrag_:True, enableTherm_:True,
  enableECond_:True, enableICond_:True] := 
Module[{Nr, rMaxEff, rCorrMaxBYmu, rArr, rFArr, eVtoMech, gamma,
  rhoFloor, rhoVac, TeFloor, ThFloor,
  rhoN, momN, epsHN, rhoP, momP, epsE, epsHP,
  rhoN1, momN1, epsHN1, rhoP1, momP1, epsE1, epsHP1,
  rhoN2, momN2, epsHN2, rhoP2, momP2, epsE2, epsHP2,
  dt, csMax, csMaxP, csMaxN, TeInit, ThInit, TeMax, TeIzGuard,
  ionCondFactor, msBYmuns, kappaConst, lotzConst, kmtVthSq, thermConst,
  TiByMui, TenZsq,
  storeIdx, nSnapshots, dtStore, tSol, rhoNSol, momNSol, epsHNSol, rhoPSol, momPSol, epsESol, epsHPSol,
  rhsFlat, rhsN, rhsP,
  compiledRHSneutral, compiledRHSplasma, compiledFloorN, compiledFloorP,
  compiledMaxCharSpeed, compiledPdVcorrection,
  applySource, enforceFloorN, enforceFloorP, inlineLnLam, compTarget},

  (* ================================================================ *)
  (* UNIT CONVERSIONS AND PHYSICAL CONSTANTS                           *)
  (* ================================================================ *)
  eVtoMech = eCharge * 10^(-6);  (* eV/cm^3 -> kg um^2/(ns^2 cm^3) *)
  gamma = 5./3;                  (* adiabatic index, monatomic ideal gas *)
  msBYmuns = 10^(-3);            (* m/s -> um/ns *)

  (* ================================================================ *)
  (* NUMERICAL THRESHOLDS                                              *)
  (* ================================================================ *)
  rhoFloor = 1*^-20;                        (* absolute min density [kg/cm^3] *)
  rhoVac = mnkg * n0BYcm\[Bullet]3 * 1*^-9; (* vacuum: matches IC density clamp *)
  TeFloor = TieV;                            (* electron T floor = ambient T *)
  ThFloor = TieV;                            (* heavy T floor = ambient T *)
  TeIzGuard = IEeV/27.;   (* skip ionization when Te < IE/27:
                              rate ~ Exp[-IE/Te] < Exp[-27] ~ 1e-12 *)
  ionCondFactor = Sqrt[me/mnkg];  (* kappa_ion/kappa_e = sqrt(me/mi) *)
  nSnapshots = Max[nSnap, 10];

  (* ================================================================ *)
  (* PRECOMPUTED CONSTANTS FOR INLINED PHYSICS                         *)
  (* ================================================================ *)
  (* Spitzer thermal conductivity [um^2 cm^-3 ns^-1 eV^-1]:
     kappa = kappaConst * Te^(5/2) / lnLambda
     From: kappa_SI = 3.2*(n*e*tau_e/me)*kB = 3.2*(e/me)*3.44e5*Te^2.5/lnLam *)
  kappaConst = N[3.2*(eCharge/me)*3.44*^5*10^3];

  (* Lotz empirical ionization rate [cm^3/s]:
     S = lotzConst * qi / (Te^0.5 * IE) * Gamma[0, IE/Te]
     where Gamma[0,x] is the upper incomplete gamma = E1(x) *)
  lotzConst = N[6.7*^-7*4.5];

  (* Momentum transfer vRMS^2 thermal contribution [um^2/ns^2 per eV]:
     For equal masses (m1=m2=m): mu=m/2, TsS=(T1+T2)/2
     vRMS^2 = (uP-uN)^2 + 8*TsS*eCharge/(Pi*mu) * 1e-6 [um^2/ns^2]
            = (uP-uN)^2 + kmtVthSq * (TieV+TneV)
     Kmt [cm^3/s] = 2.13e-9 * vRMS^0.75 *)
  kmtVthSq = N[8.*eCharge/(Pi*mnkg)*10^(-6)];

  (* Electron-ion thermalization rate constant [cm^3 ns^-1 eV^-3/2]:
     alpha = thermConst * ne^2 * lnLam / Te^1.5
     From NRL: tau_ei = 3.44e5*Te^1.5/(ne*lnLam) [s] *)
  thermConst = N[(3.*me/mnkg)*10^(-9)/3.44*^5];

  (* Coulomb logarithm regime thresholds (NRL formulary) *)
  TiByMui = TieV/mui;  (* Te threshold for regime III *)
  TenZsq = 10.*Z^2;    (* Te threshold between regimes I and II *)

  (* ================================================================ *)
  (* COMPILATION TARGET                                                *)
  (* ================================================================ *)
  compTarget = If[Length[Quiet@Check[
    Needs["CCompilerDriver`"]; CCompilers[], {}]] > 0, "C", "WVM"];

  (* ================================================================ *)
  (* INLINED COULOMB LOGARITHM (NRL Formulary)                        *)
  (* Three regimes depending on Te relative to Ti/mui and 10*Z^2.     *)
  (* Returns lnLam >= 2 always (floored for stability).               *)
  (* ================================================================ *)
  inlineLnLam[necm3_, nicm3_, TeeV_] := Module[{v},
    v = Which[
      TiByMui <= TeeV <= TenZsq,
        23. - Log[Sqrt[necm3]*Z*TeeV^(-1.5)],      (* regime I *)
      TiByMui <= TenZsq < TeeV,
        24. - Log[Sqrt[necm3]*TeeV^(-1.)],          (* regime II *)
      TeeV < TiByMui,
        30. - Log[Sqrt[nicm3]*TieV^(-1.5)*Z^2/mui], (* regime III *)
      True, 2.];
    If[v > 0., v, 2.]];

  (* ================================================================ *)
  (* DOMAIN SIZING                                                     *)
  (* Auto-extend domain if rMaxBYmu < cs_max * tFinns / 2             *)
  (* ================================================================ *)
  ThInit = TieV;
  TeMax = Max[Table[Module[{rhoPloc, epsEloc},
    rhoPloc = mnkg n0BYcm\[Bullet]3 interpNiRInit[i*drBYmu];
    epsEloc = N[interpEpsEInit[i*drBYmu]];
    If[rhoPloc < rhoVac, TeFloor, Max[(2./3)*epsEloc*mnkg/rhoPloc, TeFloor]]],
    {i, 1, Max[1, Floor[rMaxBYmu/drBYmu]]}]];
  csMaxP = Sqrt[gamma*(TeMax + ThInit)*eVtoMech/mnkg];
  csMaxN = Sqrt[gamma*ThInit*eVtoMech/mnkg];
  csMax = Max[csMaxP, csMaxN];
  rCorrMaxBYmu = csMax*tFinns/2;
  If[rMaxBYmu >= rCorrMaxBYmu,
    rMaxEff = rMaxBYmu; Nr = Round[rMaxEff/drBYmu],
    Nr = Ceiling[rCorrMaxBYmu/drBYmu]; rMaxEff = Nr*drBYmu;
    Print["*** Domain extended to ", rMaxEff, " \[Mu]m (", Nr, " cells)"]];
  rArr = Table[(i - 0.5)*drBYmu, {i, Nr}];      (* cell centers [um] *)
  rFArr = Table[1.0*i*drBYmu, {i, 0, Nr}];       (* cell faces [um], rF[0]=0 *)
  dtStore = tFinns/nSnapshots;

  Print["=== Cylindrical Plasma-Neutral Solver ==="];
  Print["Grid: ", Nr, " cells, dr = ", drBYmu, " \[Mu]m, Domain: 0 to ", rMaxEff, " \[Mu]m"];
  Print["Flux limiter f = ", fFluxLim, ", Snapshots = ", nSnapshots,
        ", Compile: ", compTarget];
  Print["Source terms: Iz=",enableIoniz," Dr=",enableDrag," Th=",enableTherm,
        " eC=",enableECond," iC=",enableICond];

  (* ================================================================ *)
  (* COMPILED FLOOR ENFORCEMENT                                        *)
  (* Density >= rhoFloor, momentum zeroed in vacuum (rho < rhoVac),    *)
  (* axis symmetry (mom[[1]]=0), energy floors density-local.          *)
  (* ================================================================ *)
  compiledFloorN = Compile[{{rhoA,_Real,1},{momA,_Real,1},{epsHA,_Real,1},
    {NrC,_Integer},{rhoFlC,_Real},{rhoVcC,_Real},{mnkgC,_Real},{ThFlC,_Real}},
    Module[{rhoO,momO,epsHO},
      rhoO=Table[Max[rhoA[[j]],rhoFlC],{j,NrC}];
      momO=Table[If[rhoA[[j]]<rhoVcC,0.,momA[[j]]],{j,NrC}]; momO[[1]]=0.;
      epsHO=Table[Max[epsHA[[j]],(3./2)*(rhoO[[j]]/mnkgC)*ThFlC],{j,NrC}];
      Join[rhoO,momO,epsHO]],CompilationTarget->compTarget,RuntimeOptions->"Speed"];
  compiledFloorP = Compile[{{rhoA,_Real,1},{momA,_Real,1},{epsEA,_Real,1},{epsHA,_Real,1},
    {NrC,_Integer},{rhoFlC,_Real},{rhoVcC,_Real},{mnkgC,_Real},{TeFlC,_Real},{ThFlC,_Real}},
    Module[{rhoO,momO,epsEO,epsHO},
      rhoO=Table[Max[rhoA[[j]],rhoFlC],{j,NrC}];
      momO=Table[If[rhoA[[j]]<rhoVcC,0.,momA[[j]]],{j,NrC}]; momO[[1]]=0.;
      epsEO=Table[Max[epsEA[[j]],(3./2)*(rhoO[[j]]/mnkgC)*TeFlC],{j,NrC}];
      epsHO=Table[Max[epsHA[[j]],(3./2)*(rhoO[[j]]/mnkgC)*ThFlC],{j,NrC}];
      Join[rhoO,momO,epsEO,epsHO]],CompilationTarget->compTarget,RuntimeOptions->"Speed"];
  enforceFloorN[rhoArr_,momArr_,epsHArr_]:=Module[{flat},
    flat=compiledFloorN[rhoArr,momArr,epsHArr,Nr,rhoFloor,rhoVac,mnkg,ThFloor];
    {flat[[1;;Nr]],flat[[Nr+1;;2Nr]],flat[[2Nr+1;;3Nr]]}];
  enforceFloorP[rhoArr_,momArr_,epsEArr_,epsHArr_]:=Module[{flat},
    flat=compiledFloorP[rhoArr,momArr,epsEArr,epsHArr,Nr,rhoFloor,rhoVac,mnkg,TeFloor,ThFloor];
    {flat[[1;;Nr]],flat[[Nr+1;;2Nr]],flat[[2Nr+1;;3Nr]],flat[[3Nr+1;;4Nr]]}];

  (* ================================================================ *)
  (* COMPILED RHS: NEUTRAL FLUID                                       *)
  (* Equations: continuity, momentum (with phN pressure via            *)
  (* advective-pressure split), energy (with pdV work).                *)
  (* Returns flat array: {drhoN, dmomN, depsHN} concatenated.          *)
  (* ================================================================ *)
  compiledRHSneutral = Compile[{{rhoA,_Real,1},{momA,_Real,1},{epsHA,_Real,1},
    {rC,_Real,1},{rFC,_Real,1},{NrC,_Integer},{drC,_Real},{gammaC,_Real},
    {eVmC,_Real},{rhoFlC,_Real},{rhoVcC,_Real}},
    Module[{u,phM,csL,flux1,flux2,fluxEH,pFA,uFA,
            rhoL,uL,pL,epsHL,csLL,rhoR,uR,pR,epsHR,csRR,
            sLL,sRR,fden,f1,f2,rhoSt,drho,dmom,depsH,divU,advF},
      u=Table[If[rhoA[[j]]<rhoVcC||Abs[momA[[j]]]<1.*^-25,0.,momA[[j]]/rhoA[[j]]],{j,NrC}];
      phM=Table[(2./3)*epsHA[[j]]*eVmC,{j,NrC}];
      csL=Table[Sqrt[Max[gammaC*phM[[j]]/Max[rhoA[[j]],rhoFlC],0.]],{j,NrC}];
      flux1=Table[0.,{NrC+1}];flux2=Table[0.,{NrC+1}];
      fluxEH=Table[0.,{NrC+1}];pFA=Table[0.,{NrC+1}];uFA=Table[0.,{NrC+1}];
      Do[
        (* Reconstruct left/right states at face j *)
        If[j==1, (* axis: ghost cell with u=0 *)
          rhoL=rhoA[[1]];uL=0.;pL=phM[[1]];epsHL=epsHA[[1]];csLL=csL[[1]];
          rhoR=rhoA[[1]];uR=u[[1]];pR=phM[[1]];epsHR=epsHA[[1]];csRR=csL[[1]],
        If[j==NrC+1, (* outer: zero-gradient *)
          rhoL=rhoA[[NrC]];uL=u[[NrC]];pL=phM[[NrC]];epsHL=epsHA[[NrC]];csLL=csL[[NrC]];
          rhoR=rhoA[[NrC]];uR=u[[NrC]];pR=phM[[NrC]];epsHR=epsHA[[NrC]];csRR=csL[[NrC]],
        (* interior *)
          rhoL=rhoA[[j-1]];uL=u[[j-1]];pL=phM[[j-1]];epsHL=epsHA[[j-1]];csLL=csL[[j-1]];
          rhoR=rhoA[[j]];uR=u[[j]];pR=phM[[j]];epsHR=epsHA[[j]];csRR=csL[[j]]]];
        (* HLL wave speeds *)
        sLL=Min[uL-csLL,uR-csRR]; sRR=Max[uL+csLL,uR+csRR];
        (* HLL flux *)
        If[sLL>=0., f1=rhoL*uL; f2=rhoL*uL*uL+pL,
        If[sRR<=0., f1=rhoR*uR; f2=rhoR*uR*uR+pR,
          fden=1./(sRR-sLL);
          f1=(sRR*rhoL*uL-sLL*rhoR*uR+sLL*sRR*(rhoR-rhoL))*fden;
          f2=(sRR*(rhoL*uL*uL+pL)-sLL*(rhoR*uR*uR+pR)+sLL*sRR*(rhoR*uR-rhoL*uL))*fden]];
        flux1[[j]]=f1; flux2[[j]]=f2;
        (* Face pressure and velocity for advective-pressure split *)
        pFA[[j]]=If[sLL>=0.,pL,If[sRR<=0.,pR,(pL+pR)/2.]];
        rhoSt=If[sLL>=0.,rhoL,If[sRR<=0.,rhoR,
          (sRR*rhoR-sLL*rhoL-(rhoR*uR-rhoL*uL))/(sRR-sLL)]];
        uFA[[j]]=If[sLL>=0.,uL,If[sRR<=0.,uR,If[rhoSt>rhoFlC,f1/rhoSt,0.]]];
        (* Energy flux: mass-flux upwinding of specific energy *)
        fluxEH[[j]]=If[f1>0.,(epsHL/Max[rhoL,rhoFlC])*f1,
                   If[f1<0.,(epsHR/Max[rhoR,rhoFlC])*f1,0.]],
      {j,1,NrC+1}];
      (* Cylindrical divergences *)
      drho=Table[-(rFC[[j+1]]*flux1[[j+1]]-rFC[[j]]*flux1[[j]])/(rC[[j]]*drC),{j,NrC}];
      advF=Table[flux2[[j]]-pFA[[j]],{j,NrC+1}];
      dmom=Table[-(rFC[[j+1]]*advF[[j+1]]-rFC[[j]]*advF[[j]])/(rC[[j]]*drC),{j,NrC}];
      dmom=Table[dmom[[j]]-(pFA[[j+1]]-pFA[[j]])/drC,{j,NrC}]; (* flat pressure gradient *)
      divU=Table[(rFC[[j+1]]*uFA[[j+1]]-rFC[[j]]*uFA[[j]])/(rC[[j]]*drC),{j,NrC}];
      depsH=Table[-(rFC[[j+1]]*fluxEH[[j+1]]-rFC[[j]]*fluxEH[[j]])/(rC[[j]]*drC),{j,NrC}];
      depsH=Table[depsH[[j]]-(2./3)*epsHA[[j]]*divU[[j]],{j,NrC}]; (* pdV work *)
      Join[drho,dmom,depsH]
    ],CompilationTarget->compTarget,RuntimeOptions->"Speed"];

  (* ================================================================ *)
  (* COMPILED RHS: PLASMA FLUID                                        *)
  (* Momentum includes pe+phP (electron + ion heavy pressure).         *)
  (* Separate pdV for epsE and epsHP using plasma velocity divergence.  *)
  (* Returns: {drhoP, dmomP, depsE, depsHP} concatenated.              *)
  (* ================================================================ *)
  compiledRHSplasma = Compile[{{rhoA,_Real,1},{momA,_Real,1},{epsEA,_Real,1},{epsHA,_Real,1},
    {rC,_Real,1},{rFC,_Real,1},{NrC,_Integer},{drC,_Real},{gammaC,_Real},
    {eVmC,_Real},{rhoFlC,_Real},{rhoVcC,_Real}},
    Module[{u,peM,phM,pTM,csL,flux1,flux2,fluxEE,fluxEH,pFA,uFA,
            rhoL,uL,pTL,epsEL,epsHL,csLL,rhoR,uR,pTR,epsER,epsHR,csRR,
            sLL,sRR,fden,f1,f2,rhoSt,drho,dmom,depsE,depsH,divU,advF},
      u=Table[If[rhoA[[j]]<rhoVcC||Abs[momA[[j]]]<1.*^-25,0.,momA[[j]]/rhoA[[j]]],{j,NrC}];
      peM=Table[(2./3)*epsEA[[j]]*eVmC,{j,NrC}];
      phM=Table[(2./3)*epsHA[[j]]*eVmC,{j,NrC}];
      pTM=Table[peM[[j]]+phM[[j]],{j,NrC}];
      csL=Table[Sqrt[Max[gammaC*pTM[[j]]/Max[rhoA[[j]],rhoFlC],0.]],{j,NrC}];
      flux1=Table[0.,{NrC+1}];flux2=Table[0.,{NrC+1}];
      fluxEE=Table[0.,{NrC+1}];fluxEH=Table[0.,{NrC+1}];
      pFA=Table[0.,{NrC+1}];uFA=Table[0.,{NrC+1}];
      Do[
        If[j==1,
          rhoL=rhoA[[1]];uL=0.;pTL=pTM[[1]];epsEL=epsEA[[1]];epsHL=epsHA[[1]];csLL=csL[[1]];
          rhoR=rhoA[[1]];uR=u[[1]];pTR=pTM[[1]];epsER=epsEA[[1]];epsHR=epsHA[[1]];csRR=csL[[1]],
        If[j==NrC+1,
          rhoL=rhoA[[NrC]];uL=u[[NrC]];pTL=pTM[[NrC]];epsEL=epsEA[[NrC]];epsHL=epsHA[[NrC]];csLL=csL[[NrC]];
          rhoR=rhoA[[NrC]];uR=u[[NrC]];pTR=pTM[[NrC]];epsER=epsEA[[NrC]];epsHR=epsHA[[NrC]];csRR=csL[[NrC]],
          rhoL=rhoA[[j-1]];uL=u[[j-1]];pTL=pTM[[j-1]];epsEL=epsEA[[j-1]];epsHL=epsHA[[j-1]];csLL=csL[[j-1]];
          rhoR=rhoA[[j]];uR=u[[j]];pTR=pTM[[j]];epsER=epsEA[[j]];epsHR=epsHA[[j]];csRR=csL[[j]]]];
        sLL=Min[uL-csLL,uR-csRR]; sRR=Max[uL+csLL,uR+csRR];
        If[sLL>=0., f1=rhoL*uL; f2=rhoL*uL*uL+pTL,
        If[sRR<=0., f1=rhoR*uR; f2=rhoR*uR*uR+pTR,
          fden=1./(sRR-sLL);
          f1=(sRR*rhoL*uL-sLL*rhoR*uR+sLL*sRR*(rhoR-rhoL))*fden;
          f2=(sRR*(rhoL*uL*uL+pTL)-sLL*(rhoR*uR*uR+pTR)+sLL*sRR*(rhoR*uR-rhoL*uL))*fden]];
        flux1[[j]]=f1; flux2[[j]]=f2;
        pFA[[j]]=If[sLL>=0.,pTL,If[sRR<=0.,pTR,(pTL+pTR)/2.]];
        rhoSt=If[sLL>=0.,rhoL,If[sRR<=0.,rhoR,
          (sRR*rhoR-sLL*rhoL-(rhoR*uR-rhoL*uL))/(sRR-sLL)]];
        uFA[[j]]=If[sLL>=0.,uL,If[sRR<=0.,uR,If[rhoSt>rhoFlC,f1/rhoSt,0.]]];
        fluxEE[[j]]=If[f1>0.,(epsEL/Max[rhoL,rhoFlC])*f1,
                   If[f1<0.,(epsER/Max[rhoR,rhoFlC])*f1,0.]];
        fluxEH[[j]]=If[f1>0.,(epsHL/Max[rhoL,rhoFlC])*f1,
                   If[f1<0.,(epsHR/Max[rhoR,rhoFlC])*f1,0.]],
      {j,1,NrC+1}];
      drho=Table[-(rFC[[j+1]]*flux1[[j+1]]-rFC[[j]]*flux1[[j]])/(rC[[j]]*drC),{j,NrC}];
      advF=Table[flux2[[j]]-pFA[[j]],{j,NrC+1}];
      dmom=Table[-(rFC[[j+1]]*advF[[j+1]]-rFC[[j]]*advF[[j]])/(rC[[j]]*drC),{j,NrC}];
      dmom=Table[dmom[[j]]-(pFA[[j+1]]-pFA[[j]])/drC,{j,NrC}];
      divU=Table[(rFC[[j+1]]*uFA[[j+1]]-rFC[[j]]*uFA[[j]])/(rC[[j]]*drC),{j,NrC}];
      depsE=Table[-(rFC[[j+1]]*fluxEE[[j+1]]-rFC[[j]]*fluxEE[[j]])/(rC[[j]]*drC),{j,NrC}];
      depsE=Table[depsE[[j]]-(2./3)*epsEA[[j]]*divU[[j]],{j,NrC}];
      depsH=Table[-(rFC[[j+1]]*fluxEH[[j+1]]-rFC[[j]]*fluxEH[[j]])/(rC[[j]]*drC),{j,NrC}];
      depsH=Table[depsH[[j]]-(2./3)*epsHA[[j]]*divU[[j]],{j,NrC}];
      Join[drho,dmom,depsE,depsH]
    ],CompilationTarget->compTarget,RuntimeOptions->"Speed"];

  (* ================================================================ *)
  (* COMPILED: MAX CHARACTERISTIC SPEED (for adaptive CFL)             *)
  (* Returns max(|u|+cs) over all cells, both species.                 *)
  (* ================================================================ *)
  compiledMaxCharSpeed = Compile[{{rhoPA,_Real,1},{momPA,_Real,1},{epsEA,_Real,1},{epsHPA,_Real,1},
    {rhoNA,_Real,1},{momNA,_Real,1},{epsHNA,_Real,1},
    {NrC,_Integer},{gammaC,_Real},{eVmC,_Real},{rhoFlC,_Real},{rhoVcC,_Real}},
    Module[{uP,uN,pTotP,phN,csP,csN,sMax},
      sMax=0.;
      Do[
        uP=If[rhoPA[[j]]<rhoVcC||Abs[momPA[[j]]]<1.*^-25,0.,momPA[[j]]/rhoPA[[j]]];
        uN=If[rhoNA[[j]]<rhoVcC||Abs[momNA[[j]]]<1.*^-25,0.,momNA[[j]]/rhoNA[[j]]];
        pTotP=(2./3)*(epsEA[[j]]+epsHPA[[j]])*eVmC;
        phN=(2./3)*epsHNA[[j]]*eVmC;
        csP=Sqrt[Max[gammaC*pTotP/Max[rhoPA[[j]],rhoFlC],0.]];
        csN=Sqrt[Max[gammaC*phN/Max[rhoNA[[j]],rhoFlC],0.]];
        sMax=Max[sMax,Abs[uP]+csP,Abs[uN]+csN],
      {j,NrC}];
      sMax
    ],CompilationTarget->compTarget,RuntimeOptions->"Speed"];

  (* ================================================================ *)
  (* COMPILED: pdV CORRECTION FOR DRAG STEP                            *)
  (* Computes the change in velocity divergence caused by drag and     *)
  (* applies the corresponding pdV work to all energy densities.       *)
  (* This corrects the Strang splitting error that otherwise causes    *)
  (* spurious heating in drag-deceleration regions.                    *)
  (* ================================================================ *)
  compiledPdVcorrection = Compile[{{uPpre,_Real,1},{uNpre,_Real,1},
    {uPpost,_Real,1},{uNpost,_Real,1},
    {epsEA,_Real,1},{epsHPA,_Real,1},{epsHNA,_Real,1},
    {rC,_Real,1},{rFC,_Real,1},{NrC,_Integer},{drC,_Real},{dtC,_Real}},
    Module[{divPpre,divPpost,divNpre,divNpost,uF,ddP,ddN,epsEO,epsHPO,epsHNO},
      divPpre=Table[0.,{NrC}];divPpost=Table[0.,{NrC}];
      divNpre=Table[0.,{NrC}];divNpost=Table[0.,{NrC}];
      Do[
        (* Plasma divU before drag *)
        uF=If[j<NrC,(uPpre[[j]]+uPpre[[j+1]])/2.,uPpre[[NrC]]];
        divPpre[[j]]=rFC[[j+1]]*uF;
        uF=If[j>1,(uPpre[[j-1]]+uPpre[[j]])/2.,0.];
        divPpre[[j]]=(divPpre[[j]]-rFC[[j]]*uF)/(rC[[j]]*drC);
        (* Plasma divU after drag *)
        uF=If[j<NrC,(uPpost[[j]]+uPpost[[j+1]])/2.,uPpost[[NrC]]];
        divPpost[[j]]=rFC[[j+1]]*uF;
        uF=If[j>1,(uPpost[[j-1]]+uPpost[[j]])/2.,0.];
        divPpost[[j]]=(divPpost[[j]]-rFC[[j]]*uF)/(rC[[j]]*drC);
        (* Neutral divU before drag *)
        uF=If[j<NrC,(uNpre[[j]]+uNpre[[j+1]])/2.,uNpre[[NrC]]];
        divNpre[[j]]=rFC[[j+1]]*uF;
        uF=If[j>1,(uNpre[[j-1]]+uNpre[[j]])/2.,0.];
        divNpre[[j]]=(divNpre[[j]]-rFC[[j]]*uF)/(rC[[j]]*drC);
        (* Neutral divU after drag *)
        uF=If[j<NrC,(uNpost[[j]]+uNpost[[j+1]])/2.,uNpost[[NrC]]];
        divNpost[[j]]=rFC[[j+1]]*uF;
        uF=If[j>1,(uNpost[[j-1]]+uNpost[[j]])/2.,0.];
        divNpost[[j]]=(divNpost[[j]]-rFC[[j]]*uF)/(rC[[j]]*drC),
      {j,NrC}];
      ddP=Table[divPpost[[j]]-divPpre[[j]],{j,NrC}];
      ddN=Table[divNpost[[j]]-divNpre[[j]],{j,NrC}];
      (* Apply pdV correction: eps *= (1 - (2/3)*DdivU*dt) *)
      epsEO=Table[epsEA[[j]]*(1.-(2./3)*ddP[[j]]*dtC),{j,NrC}];
      epsHPO=Table[epsHPA[[j]]*(1.-(2./3)*ddP[[j]]*dtC),{j,NrC}];
      epsHNO=Table[epsHNA[[j]]*(1.-(2./3)*ddN[[j]]*dtC),{j,NrC}];
      Join[epsEO,epsHPO,epsHNO]
    ],CompilationTarget->compTarget,RuntimeOptions->"Speed"];

  (* ================================================================ *)
  (* SOURCE TERMS (Strang half-step)                                   *)
  (*                                                                   *)
  (* A. Ionization: Lotz formula, explicit. Mass, momentum, energy     *)
  (*    transferred from neutral to plasma. KE dissipation -> heat.    *)
  (* B. Drag: implicit analytical. KE dissipation -> heat.             *)
  (*    Includes compiled pdV correction for splitting consistency.     *)
  (* C. Thermalization: implicit 2x2 analytical solve for Te, Th.      *)
  (*    Conserves total thermal energy (ne*Te + nh*Th = const).        *)
  (* D. Electron Spitzer conduction: implicit tridiagonal with         *)
  (*    harmonic-mean flux limiter.                                    *)
  (* E. Ion Spitzer conduction: same, with kappa_ion = sqrt(me/mi)*    *)
  (*    kappa_e(Th). Heat capacity = (3/2)*n_total (shared Th).        *)
  (*                                                                   *)
  (* Floor enforcement: after A (mass transfer safety) and after E.    *)
  (* ================================================================ *)
  applySource[rhoNIn_,momNIn_,epsHNIn_,rhoPIn_,momPIn_,epsEIn_,epsHPIn_,dtSrc_]:=
  Module[{rhoNw,momNw,epsHNw,rhoPw,momPw,epsEw,epsHPw,uN,uP,Te,Th,
          Sn,Srho,Qiz,dKEiz,
          uPpre,uNpre,Kmt,vRMS,nudP,nudN,duOld,duNew,KEbefore,KEafter,dKEdrag,
          ne,nh,nicm3,necm3,lnLam,alpha,Ae,Ah,beta,det,TeNewT,ThNewT,
          kappaFace,TeFace,ThFace,neFace,nhFace,dTdr,qSpitz,qfs,vth,
          aCoef,bCoef,cCoef,dCoef,Cplus,Cminus,C0,w,TeNewCond,ThNewCond,
          epsHtot,fP,pdvFlat,j,flat},
    rhoNw=rhoNIn;momNw=momNIn;epsHNw=epsHNIn;
    rhoPw=rhoPIn;momPw=momPIn;epsEw=epsEIn;epsHPw=epsHPIn;

    (* Precompute velocities and temperatures *)
    uP=Table[If[rhoPw[[j]]<rhoVac||Abs[momPw[[j]]]<1*^-25,0.,momPw[[j]]/rhoPw[[j]]],{j,Nr}];
    uN=Table[If[rhoNw[[j]]<rhoVac||Abs[momNw[[j]]]<1*^-25,0.,momNw[[j]]/rhoNw[[j]]],{j,Nr}];
    Te=Table[If[rhoPw[[j]]<rhoVac,TeFloor,Max[(2./3)*epsEw[[j]]*mnkg/rhoPw[[j]],TeFloor]],{j,Nr}];
    Th=Table[If[(rhoPw[[j]]+rhoNw[[j]])<rhoVac,ThFloor,
      Max[(2./3)*(epsHPw[[j]]+epsHNw[[j]])*mnkg/(rhoPw[[j]]+rhoNw[[j]]),ThFloor]],{j,Nr}];

    (* ---- A: Collisional Ionization (explicit, Lotz formula) ---- *)
    If[enableIoniz,
    Sn=Table[If[Te[[j]]<TeIzGuard,0.,
      lotzConst/(Te[[j]]^0.5*IEeV)*Quiet[N[Gamma[0,IEeV/Te[[j]]]]]*
        (rhoPw[[j]]/mnkg)*(rhoNw[[j]]/mnkg)*10^(-9)],{j,Nr}];
    Srho=mnkg*Sn; Qiz=IEeV*Sn;
    dKEiz=0.5*(uP-uN)^2*Srho/eVtoMech;
    rhoPw+=dtSrc*Srho; rhoNw-=dtSrc*Srho;
    momPw+=dtSrc*uN*Srho; momNw-=dtSrc*uN*Srho;
    epsEw-=dtSrc*Qiz;
    epsHPw+=dtSrc*1.5*Sn*Th; epsHNw-=dtSrc*1.5*Sn*Th;
    Do[fP=rhoPw[[j]]/Max[rhoPw[[j]]+rhoNw[[j]],rhoFloor];
      epsHPw[[j]]+=dtSrc*dKEiz[[j]]*fP;
      epsHNw[[j]]+=dtSrc*dKEiz[[j]]*(1.-fP),{j,Nr}];
    flat=compiledFloorP[rhoPw,momPw,epsEw,epsHPw,Nr,rhoFloor,rhoVac,mnkg,TeFloor,ThFloor];
    {rhoPw,momPw,epsEw,epsHPw}={flat[[1;;Nr]],flat[[Nr+1;;2Nr]],flat[[2Nr+1;;3Nr]],flat[[3Nr+1;;4Nr]]};
    flat=compiledFloorN[rhoNw,momNw,epsHNw,Nr,rhoFloor,rhoVac,mnkg,ThFloor];
    {rhoNw,momNw,epsHNw}={flat[[1;;Nr]],flat[[Nr+1;;2Nr]],flat[[2Nr+1;;3Nr]]}];

    (* ---- B: Ion-Neutral Drag (implicit) + pdV correction ---- *)
    If[enableDrag,
    uP=Table[If[rhoPw[[j]]<rhoVac||Abs[momPw[[j]]]<1*^-25,0.,momPw[[j]]/rhoPw[[j]]],{j,Nr}];
    uN=Table[If[rhoNw[[j]]<rhoVac||Abs[momNw[[j]]]<1*^-25,0.,momNw[[j]]/rhoNw[[j]]],{j,Nr}];
    uPpre=Table[uP[[j]],{j,Nr}]; uNpre=Table[uN[[j]],{j,Nr}];
    KEbefore=0.5*rhoPw*uP^2+0.5*rhoNw*uN^2;
    Do[If[Abs[uP[[j]]-uN[[j]]]>1*^-20,
      vRMS=Sqrt[(uP[[j]]-uN[[j]])^2+kmtVthSq*(TieV+TneV)];
      Kmt=2.13*^-9*vRMS^0.75*10^(-9);
      nudP=(0.5)*(rhoNw[[j]]/mnkg)*Kmt;nudN=(0.5)*(rhoPw[[j]]/mnkg)*Kmt;
      duOld=uP[[j]]-uN[[j]];duNew=duOld/(1.+(nudP+nudN)*dtSrc);
      uP[[j]]-=nudP*dtSrc*duNew;uN[[j]]+=nudN*dtSrc*duNew],{j,Nr}];
    momPw=rhoPw*uP;momNw=rhoNw*uN;
    KEafter=0.5*rhoPw*uP^2+0.5*rhoNw*uN^2;
    dKEdrag=(KEbefore-KEafter)/eVtoMech;
    Do[fP=rhoPw[[j]]/Max[rhoPw[[j]]+rhoNw[[j]],rhoFloor];
      epsHPw[[j]]+=dKEdrag[[j]]*fP;epsHNw[[j]]+=dKEdrag[[j]]*(1.-fP),{j,Nr}];
    (* pdV correction: compensate for Strang splitting delay *)
    pdvFlat=compiledPdVcorrection[uPpre,uNpre,uP,uN,epsEw,epsHPw,epsHNw,
      rArr,rFArr,Nr,drBYmu,dtSrc];
    epsEw=pdvFlat[[1;;Nr]];epsHPw=pdvFlat[[Nr+1;;2Nr]];epsHNw=pdvFlat[[2Nr+1;;3Nr]]];

    (* ---- C: Electron-Ion Thermalization (implicit 2x2) ---- *)
    If[enableTherm,
    Te=Table[If[rhoPw[[j]]<rhoVac,TeFloor,Max[(2./3)*epsEw[[j]]*mnkg/rhoPw[[j]],TeFloor]],{j,Nr}];
    Th=Table[If[(rhoPw[[j]]+rhoNw[[j]])<rhoVac,ThFloor,
      Max[(2./3)*(epsHPw[[j]]+epsHNw[[j]])*mnkg/(rhoPw[[j]]+rhoNw[[j]]),ThFloor]],{j,Nr}];
    Do[ne=rhoPw[[j]]/mnkg; nh=(rhoPw[[j]]+rhoNw[[j]])/mnkg; necm3=Z*ne;
      If[ne>rhoVac/mnkg && Te[[j]]>TeFloor &&
         Abs[Te[[j]]-Th[[j]]]/Max[Te[[j]],ThFloor]>1*^-6,
        lnLam=inlineLnLam[necm3,ne,Te[[j]]];
        alpha=thermConst*necm3^2*lnLam/Te[[j]]^1.5;
        Ae=1.5*ne; Ah=1.5*nh; beta=alpha*dtSrc;
        det=Ae*Ah+beta*(Ae+Ah);
        TeNewT=Max[(Ae*Ah*Te[[j]]+beta*(Ae*Te[[j]]+Ah*Th[[j]]))/det,TeFloor];
        ThNewT=Max[(Ae*Ah*Th[[j]]+beta*(Ae*Te[[j]]+Ah*Th[[j]]))/det,ThFloor];
        epsEw[[j]]=1.5*ne*TeNewT; epsHtot=1.5*nh*ThNewT;
        fP=rhoPw[[j]]/Max[rhoPw[[j]]+rhoNw[[j]],rhoFloor];
        epsHPw[[j]]=epsHtot*fP; epsHNw[[j]]=epsHtot*(1.-fP)],
    {j,Nr}]];

    (* ---- D: Electron Spitzer Conduction (flux-limited, implicit) ---- *)
    If[enableECond,
    Te=Table[If[rhoPw[[j]]<rhoVac,TeFloor,Max[(2./3)*epsEw[[j]]*mnkg/rhoPw[[j]],TeFloor]],{j,Nr}];
    aCoef=Table[0.,{Nr}];bCoef=Table[0.,{Nr}];cCoef=Table[0.,{Nr}];dCoef=Table[0.,{Nr}];
    Do[C0=1.5*(rhoPw[[j]]/mnkg)/dtSrc;
      If[j<Nr&&Te[[j]]>=TeFloor&&Te[[j+1]]>=TeFloor,
        TeFace=(Te[[j]]+Te[[j+1]])/2.;
        necm3=Z*(rhoPw[[j]]+rhoPw[[j+1]])/(2.*mnkg);
        nicm3=(rhoPw[[j]]+rhoPw[[j+1]])/(2.*mnkg);
        lnLam=inlineLnLam[necm3,nicm3,TeFace];
        kappaFace=kappaConst*TeFace^2.5/lnLam;
        neFace=necm3; dTdr=Abs[(Te[[j+1]]-Te[[j]])/drBYmu];
        vth=Sqrt[Max[eCharge*TeFace/me,0.]]*msBYmuns;
        qfs=fFluxLim*neFace*TeFace*vth; qSpitz=kappaFace*dTdr;
        kappaFace=If[qSpitz+qfs>0.,kappaFace*qfs/(qSpitz+qfs),0.];
        Cplus=rFArr[[j+1]]*kappaFace/(rArr[[j]]*drBYmu^2),Cplus=0.];
      If[j>1&&Te[[j-1]]>=TeFloor&&Te[[j]]>=TeFloor,
        TeFace=(Te[[j-1]]+Te[[j]])/2.;
        necm3=Z*(rhoPw[[j-1]]+rhoPw[[j]])/(2.*mnkg);
        nicm3=(rhoPw[[j-1]]+rhoPw[[j]])/(2.*mnkg);
        lnLam=inlineLnLam[necm3,nicm3,TeFace];
        kappaFace=kappaConst*TeFace^2.5/lnLam;
        neFace=necm3; dTdr=Abs[(Te[[j]]-Te[[j-1]])/drBYmu];
        vth=Sqrt[Max[eCharge*TeFace/me,0.]]*msBYmuns;
        qfs=fFluxLim*neFace*TeFace*vth; qSpitz=kappaFace*dTdr;
        kappaFace=If[qSpitz+qfs>0.,kappaFace*qfs/(qSpitz+qfs),0.];
        Cminus=rFArr[[j]]*kappaFace/(rArr[[j]]*drBYmu^2),Cminus=0.];
      aCoef[[j]]=-Cminus;bCoef[[j]]=C0+Cplus+Cminus;
      cCoef[[j]]=-Cplus;dCoef[[j]]=C0*Te[[j]],{j,1,Nr}];
    Do[w=aCoef[[j]]/bCoef[[j-1]];bCoef[[j]]-=w*cCoef[[j-1]];dCoef[[j]]-=w*dCoef[[j-1]],{j,2,Nr}];
    TeNewCond=Table[0.,{Nr}];TeNewCond[[Nr]]=dCoef[[Nr]]/bCoef[[Nr]];
    Do[TeNewCond[[j]]=(dCoef[[j]]-cCoef[[j]]*TeNewCond[[j+1]])/bCoef[[j]],{j,Nr-1,1,-1}];
    TeNewCond=Map[Max[#,TeFloor]&,TeNewCond];
    epsEw=1.5*(rhoPw/mnkg)*TeNewCond];

    (* ---- E: Ion Spitzer Conduction (flux-limited, implicit) ---- *)
    If[enableICond,
    Th=Table[If[(rhoPw[[j]]+rhoNw[[j]])<rhoVac,ThFloor,
      Max[(2./3)*(epsHPw[[j]]+epsHNw[[j]])*mnkg/(rhoPw[[j]]+rhoNw[[j]]),ThFloor]],{j,Nr}];
    aCoef=Table[0.,{Nr}];bCoef=Table[0.,{Nr}];cCoef=Table[0.,{Nr}];dCoef=Table[0.,{Nr}];
    Do[C0=1.5*((rhoPw[[j]]+rhoNw[[j]])/mnkg)/dtSrc;
      If[j<Nr&&Th[[j]]>=ThFloor&&Th[[j+1]]>=ThFloor,
        ThFace=(Th[[j]]+Th[[j+1]])/2.;
        necm3=Z*(rhoPw[[j]]+rhoPw[[j+1]])/(2.*mnkg);
        nicm3=(rhoPw[[j]]+rhoPw[[j+1]])/(2.*mnkg);
        lnLam=inlineLnLam[necm3,nicm3,ThFace];
        kappaFace=ionCondFactor*kappaConst*ThFace^2.5/lnLam;
        nhFace=((rhoPw[[j]]+rhoNw[[j]])+(rhoPw[[j+1]]+rhoNw[[j+1]]))/(2.*mnkg);
        dTdr=Abs[(Th[[j+1]]-Th[[j]])/drBYmu];
        vth=Sqrt[Max[eCharge*ThFace/mnkg,0.]]*msBYmuns;
        qfs=fFluxLim*nhFace*ThFace*vth; qSpitz=kappaFace*dTdr;
        kappaFace=If[qSpitz+qfs>0.,kappaFace*qfs/(qSpitz+qfs),0.];
        Cplus=rFArr[[j+1]]*kappaFace/(rArr[[j]]*drBYmu^2),Cplus=0.];
      If[j>1&&Th[[j-1]]>=ThFloor&&Th[[j]]>=ThFloor,
        ThFace=(Th[[j-1]]+Th[[j]])/2.;
        necm3=Z*(rhoPw[[j-1]]+rhoPw[[j]])/(2.*mnkg);
        nicm3=(rhoPw[[j-1]]+rhoPw[[j]])/(2.*mnkg);
        lnLam=inlineLnLam[necm3,nicm3,ThFace];
        kappaFace=ionCondFactor*kappaConst*ThFace^2.5/lnLam;
        nhFace=((rhoPw[[j-1]]+rhoNw[[j-1]])+(rhoPw[[j]]+rhoNw[[j]]))/(2.*mnkg);
        dTdr=Abs[(Th[[j]]-Th[[j-1]])/drBYmu];
        vth=Sqrt[Max[eCharge*ThFace/mnkg,0.]]*msBYmuns;
        qfs=fFluxLim*nhFace*ThFace*vth; qSpitz=kappaFace*dTdr;
        kappaFace=If[qSpitz+qfs>0.,kappaFace*qfs/(qSpitz+qfs),0.];
        Cminus=rFArr[[j]]*kappaFace/(rArr[[j]]*drBYmu^2),Cminus=0.];
      aCoef[[j]]=-Cminus;bCoef[[j]]=C0+Cplus+Cminus;
      cCoef[[j]]=-Cplus;dCoef[[j]]=C0*Th[[j]],{j,1,Nr}];
    Do[w=aCoef[[j]]/bCoef[[j-1]];bCoef[[j]]-=w*cCoef[[j-1]];dCoef[[j]]-=w*dCoef[[j-1]],{j,2,Nr}];
    ThNewCond=Table[0.,{Nr}];ThNewCond[[Nr]]=dCoef[[Nr]]/bCoef[[Nr]];
    Do[ThNewCond[[j]]=(dCoef[[j]]-cCoef[[j]]*ThNewCond[[j+1]])/bCoef[[j]],{j,Nr-1,1,-1}];
    ThNewCond=Map[Max[#,ThFloor]&,ThNewCond];
    Do[nh=(rhoPw[[j]]+rhoNw[[j]])/mnkg;epsHtot=1.5*nh*ThNewCond[[j]];
      fP=rhoPw[[j]]/Max[rhoPw[[j]]+rhoNw[[j]],rhoFloor];
      epsHPw[[j]]=epsHtot*fP;epsHNw[[j]]=epsHtot*(1.-fP),{j,Nr}]];

    (* Final floor enforcement *)
    flat=compiledFloorP[rhoPw,momPw,epsEw,epsHPw,Nr,rhoFloor,rhoVac,mnkg,TeFloor,ThFloor];
    {rhoPw,momPw,epsEw,epsHPw}={flat[[1;;Nr]],flat[[Nr+1;;2Nr]],flat[[2Nr+1;;3Nr]],flat[[3Nr+1;;4Nr]]};
    flat=compiledFloorN[rhoNw,momNw,epsHNw,Nr,rhoFloor,rhoVac,mnkg,ThFloor];
    {rhoNw,momNw,epsHNw}={flat[[1;;Nr]],flat[[Nr+1;;2Nr]],flat[[2Nr+1;;3Nr]]};
    {rhoNw,momNw,epsHNw,rhoPw,momPw,epsEw,epsHPw}];

  (* ================================================================ *)
  (* INITIAL CONDITIONS                                                *)
  (* ================================================================ *)
  rhoN = Table[mnkg n0BYcm\[Bullet]3 (1 - interpNiRInit[rArr[[i]]]), {i, Nr}];
  momN = Table[0., {i, Nr}];
  rhoP = Table[mnkg n0BYcm\[Bullet]3 interpNiRInit[rArr[[i]]], {i, Nr}];
  momP = Table[0., {i, Nr}];
  epsE = Table[N[interpEpsEInit[rArr[[i]]]], {i, Nr}];
  epsHP = Table[1.5*(rhoP[[i]]/mnkg)*ThInit, {i, Nr}];
  epsHN = Table[1.5*(rhoN[[i]]/mnkg)*ThInit, {i, Nr}];

  (* Recompute sound speed from actual grid *)
  TeInit = Table[If[rhoP[[i]]<rhoVac, TeFloor,
    Max[(2./3)*epsE[[i]]*mnkg/rhoP[[i]], TeFloor]], {i, Nr}];
  TeMax = Max[TeInit];
  csMaxP = Sqrt[gamma*(TeMax + ThInit)*eVtoMech/mnkg];
  csMaxN = Sqrt[gamma*ThInit*eVtoMech/mnkg];
  csMax = Max[csMaxP, csMaxN];

  Print["Initial Te: [", NumberForm[Min[TeInit],{3,3}], ", ",
    NumberForm[TeMax,{3,3}], "] eV, Th = ", ThInit, " eV"];
  Print["Max cs: P = ", NumberForm[csMaxP,{4,2}],
    ", N = ", NumberForm[csMaxN,{4,2}], " \[Mu]m/ns"];
  Print["Initial dt = ", NumberForm[0.4*drBYmu/csMax,{4,6}], " ns (adaptive)"];

  (* ================================================================ *)
  (* MAIN TIME LOOP: While tSim < tFinns, adaptive dt                 *)
  (*   Strang splitting: Source(dt/2) -> SSP-RK3(dt) -> Source(dt/2)  *)
  (* ================================================================ *)
  Module[{CFL = 0.4, tSim, nStep, tNextStore, tNextPrint,
          tWallStart, tWallNow, tElapsed, progressFrac = 0., sCharMax},

  (* Pre-allocate solution storage with margin *)
  tSol = Table[0., {nSnapshots + 10}];
  rhoNSol = Table[Table[0., {Nr}], {nSnapshots + 10}];
  momNSol = Table[Table[0., {Nr}], {nSnapshots + 10}];
  epsHNSol = Table[Table[0., {Nr}], {nSnapshots + 10}];
  rhoPSol = Table[Table[0., {Nr}], {nSnapshots + 10}];
  momPSol = Table[Table[0., {Nr}], {nSnapshots + 10}];
  epsESol = Table[Table[0., {Nr}], {nSnapshots + 10}];
  epsHPSol = Table[Table[0., {Nr}], {nSnapshots + 10}];

  (* Store initial state *)
  storeIdx = 1; tSol[[1]] = 0.;
  rhoNSol[[1]] = rhoN; momNSol[[1]] = momN; epsHNSol[[1]] = epsHN;
  rhoPSol[[1]] = rhoP; momPSol[[1]] = momP;
  epsESol[[1]] = epsE; epsHPSol[[1]] = epsHP;

  tSim = 0.; nStep = 0;
  tNextStore = dtStore; tNextPrint = tFinns/10;
  tWallStart = AbsoluteTime[];

  Monitor[
  While[tSim < tFinns,
    nStep++;

    (* Adaptive CFL: dt = CFL * dr / max(|u| + cs) *)
    sCharMax = compiledMaxCharSpeed[rhoP,momP,epsE,epsHP,rhoN,momN,epsHN,
                 Nr,gamma,eVtoMech,rhoFloor,rhoVac];
    sCharMax = Max[sCharMax, 1.*^-10];
    dt = CFL*drBYmu/sCharMax;
    If[tSim + dt > tFinns, dt = tFinns - tSim];

    (* Strang: first source half-step *)
    {rhoN,momN,epsHN,rhoP,momP,epsE,epsHP} =
      applySource[rhoN,momN,epsHN,rhoP,momP,epsE,epsHP,dt/2];

    (* SSP-RK3 Stage 1: U1 = U^n + dt*L(U^n) *)
    rhsFlat = compiledRHSneutral[rhoN,momN,epsHN,rArr,rFArr,Nr,drBYmu,gamma,eVtoMech,rhoFloor,rhoVac];
    rhsN = {rhsFlat[[1;;Nr]], rhsFlat[[Nr+1;;2Nr]], rhsFlat[[2Nr+1;;3Nr]]};
    rhsFlat = compiledRHSplasma[rhoP,momP,epsE,epsHP,rArr,rFArr,Nr,drBYmu,gamma,eVtoMech,rhoFloor,rhoVac];
    rhsP = {rhsFlat[[1;;Nr]], rhsFlat[[Nr+1;;2Nr]], rhsFlat[[2Nr+1;;3Nr]], rhsFlat[[3Nr+1;;4Nr]]};
    rhoN1=rhoN+dt*rhsN[[1]]; momN1=momN+dt*rhsN[[2]]; epsHN1=epsHN+dt*rhsN[[3]];
    rhoP1=rhoP+dt*rhsP[[1]]; momP1=momP+dt*rhsP[[2]]; epsE1=epsE+dt*rhsP[[3]]; epsHP1=epsHP+dt*rhsP[[4]];
    {rhoN1,momN1,epsHN1} = enforceFloorN[rhoN1,momN1,epsHN1];
    {rhoP1,momP1,epsE1,epsHP1} = enforceFloorP[rhoP1,momP1,epsE1,epsHP1];

    (* SSP-RK3 Stage 2: U2 = (3/4)U^n + (1/4)(U1 + dt*L(U1)) *)
    rhsFlat = compiledRHSneutral[rhoN1,momN1,epsHN1,rArr,rFArr,Nr,drBYmu,gamma,eVtoMech,rhoFloor,rhoVac];
    rhsN = {rhsFlat[[1;;Nr]], rhsFlat[[Nr+1;;2Nr]], rhsFlat[[2Nr+1;;3Nr]]};
    rhsFlat = compiledRHSplasma[rhoP1,momP1,epsE1,epsHP1,rArr,rFArr,Nr,drBYmu,gamma,eVtoMech,rhoFloor,rhoVac];
    rhsP = {rhsFlat[[1;;Nr]], rhsFlat[[Nr+1;;2Nr]], rhsFlat[[2Nr+1;;3Nr]], rhsFlat[[3Nr+1;;4Nr]]};
    rhoN2=(3./4)*rhoN+(1./4)*(rhoN1+dt*rhsN[[1]]);
    momN2=(3./4)*momN+(1./4)*(momN1+dt*rhsN[[2]]);
    epsHN2=(3./4)*epsHN+(1./4)*(epsHN1+dt*rhsN[[3]]);
    rhoP2=(3./4)*rhoP+(1./4)*(rhoP1+dt*rhsP[[1]]);
    momP2=(3./4)*momP+(1./4)*(momP1+dt*rhsP[[2]]);
    epsE2=(3./4)*epsE+(1./4)*(epsE1+dt*rhsP[[3]]);
    epsHP2=(3./4)*epsHP+(1./4)*(epsHP1+dt*rhsP[[4]]);
    {rhoN2,momN2,epsHN2} = enforceFloorN[rhoN2,momN2,epsHN2];
    {rhoP2,momP2,epsE2,epsHP2} = enforceFloorP[rhoP2,momP2,epsE2,epsHP2];

    (* SSP-RK3 Stage 3: U^{n+1} = (1/3)U^n + (2/3)(U2 + dt*L(U2)) *)
    rhsFlat = compiledRHSneutral[rhoN2,momN2,epsHN2,rArr,rFArr,Nr,drBYmu,gamma,eVtoMech,rhoFloor,rhoVac];
    rhsN = {rhsFlat[[1;;Nr]], rhsFlat[[Nr+1;;2Nr]], rhsFlat[[2Nr+1;;3Nr]]};
    rhsFlat = compiledRHSplasma[rhoP2,momP2,epsE2,epsHP2,rArr,rFArr,Nr,drBYmu,gamma,eVtoMech,rhoFloor,rhoVac];
    rhsP = {rhsFlat[[1;;Nr]], rhsFlat[[Nr+1;;2Nr]], rhsFlat[[2Nr+1;;3Nr]], rhsFlat[[3Nr+1;;4Nr]]};
    rhoN=(1./3)*rhoN+(2./3)*(rhoN2+dt*rhsN[[1]]);
    momN=(1./3)*momN+(2./3)*(momN2+dt*rhsN[[2]]);
    epsHN=(1./3)*epsHN+(2./3)*(epsHN2+dt*rhsN[[3]]);
    rhoP=(1./3)*rhoP+(2./3)*(rhoP2+dt*rhsP[[1]]);
    momP=(1./3)*momP+(2./3)*(momP2+dt*rhsP[[2]]);
    epsE=(1./3)*epsE+(2./3)*(epsE2+dt*rhsP[[3]]);
    epsHP=(1./3)*epsHP+(2./3)*(epsHP2+dt*rhsP[[4]]);
    {rhoN,momN,epsHN} = enforceFloorN[rhoN,momN,epsHN];
    {rhoP,momP,epsE,epsHP} = enforceFloorP[rhoP,momP,epsE,epsHP];

    (* Strang: second source half-step *)
    {rhoN,momN,epsHN,rhoP,momP,epsE,epsHP} =
      applySource[rhoN,momN,epsHN,rhoP,momP,epsE,epsHP,dt/2];

    tSim += dt;
    progressFrac = tSim/tFinns;

    (* Store at uniform time intervals *)
    If[(tSim >= tNextStore || tSim >= tFinns) && storeIdx < Length[tSol],
      storeIdx++; tSol[[storeIdx]] = tSim;
      rhoNSol[[storeIdx]] = rhoN; momNSol[[storeIdx]] = momN; epsHNSol[[storeIdx]] = epsHN;
      rhoPSol[[storeIdx]] = rhoP; momPSol[[storeIdx]] = momP;
      epsESol[[storeIdx]] = epsE; epsHPSol[[storeIdx]] = epsHP;
      tNextStore += dtStore];

    (* Print progress at 10% intervals *)
    If[tSim >= tNextPrint,
      tWallNow = AbsoluteTime[]; tElapsed = tWallNow - tWallStart;
      Print[Round[100.*progressFrac], "%  t=", NumberForm[tSim,{4,3}], " ns  step ", nStep,
        "  dt=", NumberForm[dt,{3,4}], " ns  cs=", NumberForm[sCharMax,{4,2}],
        "  Te(ax)=", NumberForm[Max[(2./3)*epsE[[1]]*mnkg/Max[rhoP[[1]],rhoFloor],TeFloor],{3,2}], " eV",
        "  Th(ax)=", NumberForm[Max[(2./3)*(epsHP[[1]]+epsHN[[1]])*mnkg/Max[rhoP[[1]]+rhoN[[1]],rhoFloor],ThFloor],{3,2}], " eV",
        "  [", NumberForm[1000.*tElapsed/nStep,{4,1}], " ms/step]"];
      tNextPrint += tFinns/10];
  ],
  (* Dynamic progress bar *)
  Module[{frac, elap, sps},
    frac = progressFrac; elap = AbsoluteTime[] - tWallStart;
    sps = If[nStep > 0, elap/nStep, 0.];
    Panel[Column[{
      Row[{ProgressIndicator[frac,{0,1},ImageSize->{400,20}], "  ",
           Style[ToString[Round[100. frac]] <> "%", Bold]}],
      Row[{"Step ", nStep, "   t=", NumberForm[tSim,{4,3}], "/",
           NumberForm[tFinns,{4,3}], " ns   dt=", NumberForm[dt,{3,4}], " ns"}],
      Row[{"Elapsed: ",
           If[elap<60, ToString[Round[elap,0.1]]<>" sec",
           If[elap<3600, ToString[Round[elap/60,0.1]]<>" min",
             ToString[Round[elap/3600,0.1]]<>" hr"]],
           "   (", NumberForm[1000. sps,{4,1}], " ms/step)"}]
    }, Spacings -> 0.3], Style["Solver Progress", Bold, 12]]]
  ]; (* end Monitor *)

  (* Ensure final state is stored *)
  If[Abs[tSol[[storeIdx]] - tFinns] > 1.*^-12 && storeIdx < Length[tSol],
    storeIdx++; tSol[[storeIdx]] = tSim;
    rhoNSol[[storeIdx]] = rhoN; momNSol[[storeIdx]] = momN; epsHNSol[[storeIdx]] = epsHN;
    rhoPSol[[storeIdx]] = rhoP; momPSol[[storeIdx]] = momP;
    epsESol[[storeIdx]] = epsE; epsHPSol[[storeIdx]] = epsHP];

  (* Timing summary *)
  tWallNow = AbsoluteTime[]; tElapsed = tWallNow - tWallStart;
  Print["\nCompleted ", nStep, " steps in ",
    If[tElapsed<60, ToString[Round[tElapsed,0.1]]<>" sec",
    If[tElapsed<3600, ToString[Round[tElapsed/60,0.1]]<>" min",
      ToString[Round[tElapsed/3600,0.01]]<>" hr"]],
    "  (", NumberForm[1000. tElapsed/nStep,{4,1}], " ms/step avg)"];
  ]; (* end Module for time loop *)

  (* ================================================================ *)
  (* TRIM AND RETURN                                                   *)
  (* ================================================================ *)
  tSol = tSol[[1;;storeIdx]];
  rhoNSol = rhoNSol[[1;;storeIdx]]; momNSol = momNSol[[1;;storeIdx]];
  epsHNSol = epsHNSol[[1;;storeIdx]];
  rhoPSol = rhoPSol[[1;;storeIdx]]; momPSol = momPSol[[1;;storeIdx]];
  epsESol = epsESol[[1;;storeIdx]]; epsHPSol = epsHPSol[[1;;storeIdx]];
  Print["Stored ", storeIdx, " snapshots"];

  (* Mass conservation check *)
  Print["\n=== Final Diagnostics ==="];
  Print["Mass N: ", ScientificForm[Total[rhoN*rArr]*drBYmu*2 Pi, 4],
    " (init: ", ScientificForm[Total[rhoNSol[[1]]*rArr]*drBYmu*2 Pi, 4], ")"];
  Print["Mass P: ", ScientificForm[Total[rhoP*rArr]*drBYmu*2 Pi, 4],
    " (init: ", ScientificForm[Total[rhoPSol[[1]]*rArr]*drBYmu*2 Pi, 4], ")"];

  (* Return interpolating functions for all 7 conserved variables *)
  {Flatten@Table[
    {rhoNBYkgcm\[Bullet]3[i] ->
       Interpolation[Transpose[{tSol, rhoNSol[[All, i]]}], InterpolationOrder -> 1],
     momNBYkgcm\[Bullet]3muns\[Bullet]1[i] ->
       Interpolation[Transpose[{tSol, momNSol[[All, i]]}], InterpolationOrder -> 1],
     epsHNBYeVcm\[Bullet]3[i] ->
       Interpolation[Transpose[{tSol, epsHNSol[[All, i]]}], InterpolationOrder -> 1],
     rhoPBYkgcm\[Bullet]3[i] ->
       Interpolation[Transpose[{tSol, rhoPSol[[All, i]]}], InterpolationOrder -> 1],
     momPBYkgcm\[Bullet]3muns\[Bullet]1[i] ->
       Interpolation[Transpose[{tSol, momPSol[[All, i]]}], InterpolationOrder -> 1],
     epsEBYeVcm\[Bullet]3[i] ->
       Interpolation[Transpose[{tSol, epsESol[[All, i]]}], InterpolationOrder -> 1],
     epsHPBYeVcm\[Bullet]3[i] ->
       Interpolation[Transpose[{tSol, epsHPSol[[All, i]]}], InterpolationOrder -> 1]},
  {i, Nr}]}
]
