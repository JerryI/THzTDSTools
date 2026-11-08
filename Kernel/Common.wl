BeginPackage["JerryI`THzTDSTools`Common`"]

dropHalf;
skip10;
keepMiddle;
fourier2d;
inverseFourier2d;
resample;
filterNumerical;
iHilbert;

phaseFromKK;
phaseFromFT;
phaseFromFTKKRobust;
phaseFromFTKKGuided;

Begin["`Private`"]


dropHalf[list_List] := Drop[list, -Round[Length[list]/2]]

skip10[list_List] := Drop[list, Round[Length[list]/9]]

keepMiddle[list_List] := Drop[Drop[list, -Round[Length[list]/4]], Round[Length[list]/5]]

fourier2d[sample_NumericArray] := fourier2d[sample // Normal]
fourier2d[sample_List] := fourier2d[sample, Differences[sample[[All,1]]] // Min]
fourier2d[sample_List, step_] := With[{int = Interpolation[sample, "ExtrapolationHandler" -> {(0.) &, "WarningMessage" -> True}], resampledLength =  Length[Table[0, {i, sample[[1,1]], sample[[-1,1]], step}]]},
  Transpose[{
    Table[ i, {i, 0, (*FB[*)((1)(*,*)/(*,*)(step))(*]FB*) - (*FB[*)((1)(*,*)/(*,*)(step resampledLength))(*]FB*), (*FB[*)((1)(*,*)/(*,*)(step resampledLength))(*]FB*)}],
    Fourier[Table[int[i], {i, sample[[1,1]], sample[[-1,1]], step}]]
  }] // dropHalf
]

fourier2d[sample_List, {min_, max_, step_}] := With[{int = Interpolation[sample, "ExtrapolationHandler" -> {(0.) &, "WarningMessage" -> True}], len = Length[Table[0,  {i, min, max, step}]]},
  Transpose[{
    Table[ ( (*FB[*)((1)(*,*)/(*,*)(step len step))(*]FB*) (i - min)),  {i, min, max, step}],
    Fourier[Table[int[i], {i, min, max, step}]]
  }] // dropHalf
]


inverseFourier2d[sample_NumericArray] := inverseFourier2d[sample // Normal]
inverseFourier2d[image_List] := With[
  {
    step = Differences[image[[All, 1]]] // Abs // Min,
    min = image[[All,1]] // Min,
    max = image[[All,1]] // Max,
    int = Interpolation[image, InterpolationOrder->1]
  },
    inverseFourier2d[Table[{x,  int[x]}, {x, min, max, step}], {min, max, step}]
]

inverseFourier2d[image_List, {min_, max_, step_}] := With[
  {fourierData = image[[All, 2]]},
  {fullFourier = Join[fourierData, {0.},
    Drop[Conjugate[Reverse[fourierData]], -1]]},
  {reconstructed = InverseFourier[fullFourier]},
  Transpose@{
    Table[i, {i, 0, 1/step - 1/(step Length[reconstructed]),
      1/(step Length[reconstructed])}],
    Re[reconstructed]
  }
]

filterNumerical[sample_List] := With[{f = Select[sample, MatchQ[{__?NumberQ}]]}, If[Length[f] < 100, $Failed, f]]

resample[sample:{{{_?NumberQ,_?NumberQ}..}..}] := With[{
  min = Min[sample[[All, 1,1]]],
  max = Max[sample[[All, -1,1]]],
  step = Min[(Differences[#[[All,1]]] // Min) &/@ sample]
},
  resample[#, min, max, step] &/@ sample
];

resample[sample:{{_?NumberQ,_?NumberQ}..}] := resample[sample, sample[[1,1]], sample[[-1,1]], Differences[sample[[All,1]]] // Min]
resample[sample:{{_?NumberQ,_?NumberQ}..}, min_, max_, step_] := With[{int = Interpolation[sample, "ExtrapolationHandler" -> {(0.) &, "WarningMessage" -> False}]},
  Table[{i, int[i]}, {i, min, max, step}]
]

iHilbert[x_?VectorQ] := With[
  {fopts = FourierParameters -> {1, -1}, n = Length[x]},
  {e = Boole[EvenQ[n]]},
  Im[InverseFourier[
    Fourier[x, fopts]*
     PadRight[
      ArrayPad[ConstantArray[2, Quotient[n, 2] - e], {1, e}, 1], n], 
    fopts] ]
]



    (* Each phase function returns {frequency in THz, phase in radians} pairs. *)
phaseFromKK[spectrum_List, _, {fMin_?NumericQ, fMax_?NumericQ}] := Module[
      {freq = spectrum[[All, 1]], amplitude = Abs[spectrum[[All, 2]]],
       logAmp, phase, indices},
      (* Replace a near-zero DC coefficient by the neighboring amplitude. *)
      If[amplitude[[1]] < 10^-8 Max[amplitude], amplitude[[1]] = amplitude[[2]]];
      logAmp = Log[Map[Max[#, $MachineEpsilon] &, amplitude]];
      (* Even extension supplies negative frequencies for the Hilbert transform. *)
      phase = Take[iHilbert[Join[logAmp, Reverse[logAmp[[2 ;;]]]]], Length[freq]];
      indices = Select[Range[Length[freq]], fMin <= freq[[#]] <= fMax &];
      Transpose[{freq[[indices]], phase[[indices]]}]
    ];

phaseFromKK[spectrum_List, range:{_?NumericQ, _?NumericQ}] :=
  phaseFromKK[spectrum, Null, range];

phaseArrivalTimePs[timeTrace_List] := With[
  {peakIndex = First[Ordering[Abs[timeTrace[[All, 2]]], -1]]},
  timeTrace[[peakIndex, 1]] - timeTrace[[1, 1]]
];

restoreArrivalPhase[phasePairs_List, timeTrace_List] := Transpose[{
  phasePairs[[All, 1]],
  phasePairs[[All, 2]] + 2 Pi phasePairs[[All, 1]] phaseArrivalTimePs[timeTrace]
}];

phaseFromKK[spectrum_List, timeTrace_List,
    range:{_?NumericQ, _?NumericQ}, restoreLinearPhase_] := With[
  {result = phaseFromKK[spectrum, range]},
  If[TrueQ[restoreLinearPhase], restoreArrivalPhase[result, timeTrace], result]
];
    
phaseFromFT[spectrum_List, timeTrace_List, {fMin_?NumericQ, fMax_?NumericQ},
      restoreLinearPhase_:False] := With[
      {freq = spectrum[[All, 1]], t = spectrum[[All, 2]],
       peakTimePs = phaseArrivalTimePs[timeTrace]},
      {indices = Select[Range[Length[freq]], fMin <= freq[[#]] <= fMax &]},
      (* Remove the pulse-arrival phase before unwrapping. No fitted phase term. *)
      {wrapped = Arg[t[[indices]] Exp[-I 2 Pi freq[[indices]] peakTimePs]]},
      {unwrapped = Accumulate[
        Prepend[Mod[Differences[wrapped] + Pi, 2 Pi] - Pi, First[wrapped]]
      ]},
      {result = Transpose[{freq[[indices]], unwrapped}]},
      If[TrueQ[restoreLinearPhase], restoreArrivalPhase[result, timeTrace], result]
    ];
    
    (* Use local KK phase differences to choose each 2 Pi unwrap branch. *)
phaseFromFTKKGuided[spectrum_List, timeTrace_List,
      {fMin_?NumericQ, fMax_?NumericQ}, restoreLinearPhase_:False] := With[
      {kk = phaseFromKK[spectrum, {fMin, fMax}],
       freq = spectrum[[All, 1]], t = spectrum[[All, 2]],
       peakTimePs = phaseArrivalTimePs[timeTrace]},
      {indices = Select[Range[Length[freq]], fMin <= freq[[#]] <= fMax &]},
      {wrapped = Arg[t[[indices]] Exp[-I 2 Pi freq[[indices]] peakTimePs]]},
      {rawSteps = Differences[wrapped], kkSteps = Differences[kk[[All, 2]]]},
      {guidedSteps = rawSteps + 2 Pi Round[(kkSteps - rawSteps)/(2 Pi)]},
      {unwrapped = Accumulate[Prepend[guidedSteps, First[wrapped]]]},
      {result = Transpose[{freq[[indices]], unwrapped}]},
      If[TrueQ[restoreLinearPhase], restoreArrivalPhase[result, timeTrace], result]
    ];
    
    (* Across deep transmission nulls, compare the phase residual on both sides
       at several scales. A correction is accepted only by a two-thirds vote. *)
phaseFromFTKKRobust[spectrum_List, timeTrace_List,
      {fMin_?NumericQ, fMax_?NumericQ}, restoreLinearPhase_:False] := Module[
      {guided = phaseFromFTKKGuided[spectrum, timeTrace, {fMin, fMax}, False],
       kk = phaseFromKK[spectrum, {fMin, fMax}], freq, amplitude, indices,
       n, radius, candidates, residual, left, right, votes, shift},
      freq = guided[[All, 1]];
      n = Length[freq];
      If[n < 30,
        Return[If[TrueQ[restoreLinearPhase], restoreArrivalPhase[guided, timeTrace], guided]]
      ];
      indices = Select[Range[Length[spectrum]], fMin <= spectrum[[#, 1]] <= fMax &];
      amplitude = Abs[spectrum[[indices, 2]]];
      radius = Min[Max[6, Round[Sqrt[n]/2]], Floor[(n - 1)/4]];
      (* Relative depth, so no absolute amplitude or resonance frequency is assumed. *)
      candidates = Select[
        Range[2 radius + 1, n - 2 radius],
        amplitude[[#]] < amplitude[[# - 1]] &&
        amplitude[[#]] < amplitude[[# + 1]] &&
        amplitude[[#]] < 0.05 Median[amplitude[[# - radius ;; # + radius]]] &
      ];
      Do[
        residual = guided[[All, 2]] - kk[[All, 2]];
        votes = Table[
          left = Range[i - 2 w, i - w];
          right = Range[i + w, i + 2 w];
          Round[((Fit[Transpose[{freq[[right]], residual[[right]]}], {1, u}, u] -
                  Fit[Transpose[{freq[[left]], residual[[left]]}], {1, u}, u]) /.
                 u -> freq[[i]])/(2 Pi)],
          {w, 3, radius}
        ];
        shift = Round[Median[votes]];
        If[shift != 0 && Count[votes, shift] >= Ceiling[2 Length[votes]/3],
          guided[[i + 1 ;;, 2]] = guided[[i + 1 ;;, 2]] - 2 Pi shift
        ],
        {i, candidates}
      ];
      If[TrueQ[restoreLinearPhase], restoreArrivalPhase[guided, timeTrace], guided]
    ];
    
End[]
EndPackage[]
