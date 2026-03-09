(* ============================================================================
   analyzeSolution: Analysis and Visualization Module
   ============================================================================ *)

analyzeSolution[sol_, Nr_, rMaxBYmu_, mnkg_, n0BYcm\[Bullet]3_, tFinal_, tPlotFracs_,
  nMovieFrames_:100, rPlotMax_:Automatic, nMax_:Automatic,
  TMin_:0.01, TMax_:Automatic, exportDir_:""] :=
Module[{drBYmu, rGrid, rPlotEnd, rho0, rhoMax, TMinVal, TMaxVal, exportPath,
  solRules, rhoDispMin,
  getRhoN, getRhoP, getMomN, getMomP, getEpsE, getEpsHP, getEpsHN,
  getUN, getUP, getTe, getTh,
  tDiag, massN, massP, maxUN, maxUP,
  colTe, colTh, colRhoN, colRhoP, colUN, colUP,
  frameStyle, labelStyle, tickStyle, imgSize, plotOpts,
  lblTe, lblTh, lblRhoN, lblRhoP, lblUP, lblUN,
  lblDensAxis, lblTempAxis, lblVelAxis, lblR,
  tCheckpoints, checkGrid,
  makeTempPlot, makeDensPlot, makeVelPlot,
  tMovie, tempFrames, densFrames, tempMovie, densMovie},

  (* Setup *)
  drBYmu = rMaxBYmu/Nr;
  rGrid = Table[(i - 0.5)*drBYmu, {i, Nr}];
  rPlotEnd = If[rPlotMax === Automatic, rMaxBYmu, rPlotMax];
  rho0 = mnkg*n0BYcm\[Bullet]3;
  rhoMax = If[nMax === Automatic, 2.5*rho0, mnkg*nMax];
  TMinVal = TMin;
  TMaxVal = TMax;
  rhoDispMin = rho0*1*^-7;
  solRules = First[sol];
  exportPath = If[exportDir === "",
    Quiet@Check[NotebookDirectory[], $HomeDirectory], exportDir];
  If[!DirectoryQ[exportPath], exportPath = $HomeDirectory];

  (* Formatted labels using Subscript/Superscript *)
  lblTe = Subscript["T", "e"];
  lblTh = Subscript["T", "h"];
  lblRhoN = Row[{Subscript["\[Rho]", "N"]}];
  lblRhoP = Row[{Subscript["\[Rho]", "P"]}];
  lblUP = Subscript["u", "P"];
  lblUN = Subscript["u", "N"];
  lblR = "r [\[Mu]m]";
  lblTempAxis = "Temperature [eV]";
  lblDensAxis = Row[{"\[Rho] [kg/", Superscript["cm", "3"], "]"}];
  lblVelAxis = "Velocity [\[Mu]m/ns]";

  (* Data access *)
  getRhoN[t_] := Table[(rhoNBYkgcm\[Bullet]3[i] /. solRules)[t], {i, Nr}];
  getRhoP[t_] := Table[(rhoPBYkgcm\[Bullet]3[i] /. solRules)[t], {i, Nr}];
  getMomN[t_] := Table[(momNBYkgcm\[Bullet]3muns\[Bullet]1[i] /. solRules)[t], {i, Nr}];
  getMomP[t_] := Table[(momPBYkgcm\[Bullet]3muns\[Bullet]1[i] /. solRules)[t], {i, Nr}];
  getEpsE[t_] := Table[(epsEBYeVcm\[Bullet]3[i] /. solRules)[t], {i, Nr}];
  getEpsHP[t_] := Table[(epsHPBYeVcm\[Bullet]3[i] /. solRules)[t], {i, Nr}];
  getEpsHN[t_] := Table[(epsHNBYeVcm\[Bullet]3[i] /. solRules)[t], {i, Nr}];
  getUN[t_] := Module[{rN = getRhoN[t], mN = getMomN[t]},
    Table[If[rN[[i]] > rhoDispMin, mN[[i]]/rN[[i]], 0.], {i, Nr}]];
  getUP[t_] := Module[{rP = getRhoP[t], mP = getMomP[t]},
    Table[If[rP[[i]] > rhoDispMin, mP[[i]]/rP[[i]], 0.], {i, Nr}]];
  getTe[t_] := Module[{rP = getRhoP[t], eE = getEpsE[t]},
    Table[If[rP[[i]] < rhoDispMin, 0.025,
      Max[(2./3)*eE[[i]]*mnkg/rP[[i]], 0.025]], {i, Nr}]];
  getTh[t_] := Module[{rP = getRhoP[t], rN = getRhoN[t],
      eHP = getEpsHP[t], eHN = getEpsHN[t]},
    Table[If[(rP[[i]] + rN[[i]]) < rhoDispMin, 0.025,
      Max[(2./3)*(eHP[[i]] + eHN[[i]])*mnkg/(rP[[i]] + rN[[i]]), 0.025]], {i, Nr}]];

  (* Diagnostics *)
  Print[Style["\n\[FilledSmallSquare] SOLUTION DIAGNOSTICS", Bold, 14]];
  tDiag = Join[{0}, Table[i*tFinal/4, {i, 4}]];
  Print[Style["\nMass Conservation and Peak Velocities:", Bold]];
  Do[Module[{rN, rP, uN, uP},
    rN = getRhoN[t]; rP = getRhoP[t]; uN = getUN[t]; uP = getUP[t];
    massN = Total[rN*rGrid]*drBYmu*2*Pi;
    massP = Total[rP*rGrid]*drBYmu*2*Pi;
    maxUN = Max[Abs[uN]]; maxUP = Max[Abs[uP]];
    Print["  t = ", NumberForm[t, {3, 2}], " ns: ",
      "M(N) = ", ScientificForm[massN, 3],
      ", M(P) = ", ScientificForm[massP, 3],
      ", |uN|max = ", NumberForm[maxUN, {3, 2}],
      ", |uP|max = ", NumberForm[maxUP, {3, 2}], " \[Mu]m/ns"]],
  {t, tDiag}];
  Print[Style["\nFinal State (t = " <> ToString[tFinal] <> " ns):", Bold]];
  Module[{rN = getRhoN[tFinal], rP = getRhoP[tFinal],
      Te = getTe[tFinal], Th = getTh[tFinal]},
    Print["  rhoN: [", ScientificForm[Min[rN], 2], ", ", ScientificForm[Max[rN], 2], "]"];
    Print["  rhoP: [", ScientificForm[Min[rP], 2], ", ", ScientificForm[Max[rP], 2], "]"];
    Print["  Te:   [", NumberForm[Min[Te], {3, 3}], ", ", NumberForm[Max[Te], {3, 3}], "] eV"];
    Print["  Th:   [", NumberForm[Min[Th], {3, 3}], ", ", NumberForm[Max[Th], {3, 3}], "] eV"]];

  (* Plot styling *)
  colTe = RGBColor[0.85, 0.15, 0.15];
  colTh = RGBColor[0.15, 0.35, 0.75];
  colRhoN = RGBColor[0.20, 0.60, 0.20];
  colRhoP = RGBColor[0.60, 0.20, 0.60];
  colUN = RGBColor[0.20, 0.60, 0.20];
  colUP = RGBColor[0.60, 0.20, 0.60];
  frameStyle = Directive[Black, AbsoluteThickness[1.2]];
  labelStyle = Directive[Black, 13, FontFamily -> "Helvetica"];
  tickStyle = Directive[Black, 11];
  imgSize = 420;
  plotOpts = {Frame -> True, FrameStyle -> frameStyle,
    LabelStyle -> labelStyle, FrameTicksStyle -> tickStyle,
    ImageSize -> imgSize, ImagePadding -> {{65, 18}, {45, 30}},
    GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.85], Dashed]};

  (* Checkpoint plot generators *)
  makeTempPlot[t_] := Module[{Te, Th, tStr},
    Te = getTe[t]; Th = getTh[t];
    tStr = "t = " <> ToString[NumberForm[t, {4, 2}]] <> " ns";
    ListLogPlot[{Transpose[{rGrid, Te}], Transpose[{rGrid, Th}]},
      PlotLabel -> Style[tStr, Bold, 13],
      FrameLabel -> {lblR, lblTempAxis},
      PlotStyle -> {Directive[colTe, AbsoluteThickness[1.8]],
                    Directive[colTh, AbsoluteThickness[1.8], Dashed]},
      PlotRange -> {{0, rPlotEnd}, {TMinVal, TMaxVal}},
      PlotLegends -> Placed[LineLegend[
        {Directive[colTe, AbsoluteThickness[2]],
         Directive[colTh, AbsoluteThickness[2], Dashed]},
        {lblTe, lblTh}, LegendMarkerSize -> 20], {0.82, 0.85}],
      Joined -> True, Evaluate[plotOpts]]];

  makeDensPlot[t_] := Module[{rN, rP, tStr},
    rN = getRhoN[t]; rP = getRhoP[t];
    tStr = "t = " <> ToString[NumberForm[t, {4, 2}]] <> " ns";
    ListPlot[{Transpose[{rGrid, rN}], Transpose[{rGrid, rP}]},
      PlotLabel -> Style[tStr, Bold, 13],
      FrameLabel -> {lblR, lblDensAxis},
      PlotStyle -> {Directive[colRhoN, AbsoluteThickness[1.8]],
                    Directive[colRhoP, AbsoluteThickness[1.8], Dashed]},
      PlotRange -> {{0, rPlotEnd}, {0, rhoMax}},
      PlotLegends -> Placed[LineLegend[
        {Directive[colRhoN, AbsoluteThickness[2]],
         Directive[colRhoP, AbsoluteThickness[2], Dashed]},
        {lblRhoN, lblRhoP}, LegendMarkerSize -> 20], {0.82, 0.85}],
      Joined -> True, Evaluate[plotOpts]]];

  makeVelPlot[t_] := Module[{uN, uP, tStr},
    uN = getUN[t]; uP = getUP[t];
    tStr = "t = " <> ToString[NumberForm[t, {4, 2}]] <> " ns";
    ListPlot[{Transpose[{rGrid, uP}], Transpose[{rGrid, uN}]},
      PlotLabel -> Style[tStr, Bold, 13],
      FrameLabel -> {lblR, lblVelAxis},
      PlotStyle -> {Directive[colUP, AbsoluteThickness[1.8]],
                    Directive[colUN, AbsoluteThickness[1.8], Dashed]},
      PlotRange -> {{0, rPlotEnd}, Automatic},
      PlotLegends -> Placed[LineLegend[
        {Directive[colUP, AbsoluteThickness[2]],
         Directive[colUN, AbsoluteThickness[2], Dashed]},
        {lblUP, lblUN}, LegendMarkerSize -> 20], {0.82, 0.85}],
      Joined -> True, Evaluate[plotOpts]]];

  (* Checkpoint grid *)
  tCheckpoints = Join[{0}, tPlotFracs*tFinal];
  Print[Style["\n\[FilledSmallSquare] CHECKPOINT PLOTS", Bold, 14]];
  Print["Times: ", NumberForm[#, {4, 2}] & /@ tCheckpoints, " ns\n"];
  checkGrid = Grid[Table[{
    makeTempPlot[tCheckpoints[[i]]],
    makeDensPlot[tCheckpoints[[i]]],
    makeVelPlot[tCheckpoints[[i]]]},
    {i, Length[tCheckpoints]}],
    Spacings -> {1, 1}, Alignment -> Center];
  Print[checkGrid];

  (* Temperature movie *)
  Print[Style["\n\[FilledSmallSquare] GENERATING MOVIES", Bold, 14]];
  Print["Creating temperature movie (", nMovieFrames, " frames)..."];
  tMovie = Table[t, {t, 0, tFinal, tFinal/(nMovieFrames - 1)}];
  tempFrames = Table[Module[{Te, Th, tStr},
    Te = getTe[t]; Th = getTh[t];
    tStr = "t = " <> ToString[NumberForm[t, {4, 3}]] <> " ns";
    ListLogPlot[{Transpose[{rGrid, Te}], Transpose[{rGrid, Th}]},
      PlotLabel -> Style[tStr, Bold, 14],
      FrameLabel -> {Style[lblR, 13], Style[lblTempAxis, 13]},
      PlotStyle -> {Directive[colTe, AbsoluteThickness[2.2]],
                    Directive[colTh, AbsoluteThickness[2.2], Dashed]},
      PlotRange -> {{0, rPlotEnd}, {TMinVal, TMaxVal}},
      PlotLegends -> Placed[LineLegend[
        {Directive[colTe, AbsoluteThickness[2.5]],
         Directive[colTh, AbsoluteThickness[2.5], Dashed]},
        {Style[lblTe, 12], Style[lblTh, 12]},
        LegendMarkerSize -> 22], {0.82, 0.85}],
      Joined -> True, Frame -> True, FrameStyle -> frameStyle,
      AspectRatio -> 1, LabelStyle -> labelStyle, FrameTicksStyle -> tickStyle,
      ImageSize -> 500, ImagePadding -> {{70, 18}, {55, 40}},
      GridLines -> Automatic,
      GridLinesStyle -> Directive[GrayLevel[0.85], Dashed]]],
  {t, tMovie}];

  (* Density movie *)
  Print["Creating density movie (", nMovieFrames, " frames)..."];
  densFrames = Table[Module[{rN, rP, tStr},
    rN = getRhoN[t]; rP = getRhoP[t];
    tStr = "t = " <> ToString[NumberForm[t, {4, 3}]] <> " ns";
    ListPlot[{Transpose[{rGrid, rN}], Transpose[{rGrid, rP}]},
      PlotLabel -> Style[tStr, Bold, 14],
      FrameLabel -> {Style[lblR, 13], Style[lblDensAxis, 13]},
      PlotStyle -> {Directive[colRhoN, AbsoluteThickness[2.2]],
                    Directive[colRhoP, AbsoluteThickness[2.2], Dashed]},
      PlotRange -> {{0, rPlotEnd}, {0, rhoMax}},
      PlotLegends -> Placed[LineLegend[
        {Directive[colRhoN, AbsoluteThickness[2.5]],
         Directive[colRhoP, AbsoluteThickness[2.5], Dashed]},
        {Style[lblRhoN, 12], Style[lblRhoP, 12]},
        LegendMarkerSize -> 22], {0.82, 0.85}],
      Joined -> True, Frame -> True, FrameStyle -> frameStyle,
      AspectRatio -> 1, LabelStyle -> labelStyle, FrameTicksStyle -> tickStyle,
      ImageSize -> 500, ImagePadding -> {{70, 18}, {55, 40}},
      GridLines -> Automatic,
      GridLinesStyle -> Directive[GrayLevel[0.85], Dashed]]],
  {t, tMovie}];

  (* Export *)
  Module[{tempFile, densFile},
    tempFile = FileNameJoin[{exportPath, "temperature_evolution.mp4"}];
    densFile = FileNameJoin[{exportPath, "density_evolution.mp4"}];
    Print["Exporting to: ", exportPath];
    Print["  temperature_evolution.mp4 ..."];
    Export[tempFile, tempFrames, "FrameRate" -> 15];
    Print["  density_evolution.mp4 ..."];
    Export[densFile, densFrames, "FrameRate" -> 15];
    Print[Style["  Done.", Bold, Darker[Green]]];
    Print["  ", tempFile];
    Print["  ", densFile]];

  (* Interactive animations *)
  Print[Style["\n\[FilledSmallSquare] INTERACTIVE ANIMATIONS", Bold, 14]];
  tempMovie = ListAnimate[tempFrames, AnimationRate -> 15,
    AppearanceElements -> {"ProgressSlider", "PlayPauseButton",
      "StepLeftButton", "StepRightButton", "FasterSlowerButtons"}];
  densMovie = ListAnimate[densFrames, AnimationRate -> 15,
    AppearanceElements -> {"ProgressSlider", "PlayPauseButton",
      "StepLeftButton", "StepRightButton", "FasterSlowerButtons"}];
  Print[Panel[Column[{
    Style[Row[{"Temperature Evolution (", lblTe, " and ", lblTh, ")"}], Bold, 13],
    tempMovie}, Alignment -> Center],
    Style["Temperature", Bold, 11]]];
  Print[Panel[Column[{
    Style[Row[{"Density Evolution (", lblRhoN, " and ", lblRhoP, ")"}], Bold, 13],
    densMovie}, Alignment -> Center],
    Style["Density", Bold, 11]]];
  Print[Style["\n\[Checkmark] Analysis complete.", Bold, 14]];
]
