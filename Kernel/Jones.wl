BeginPackage["JerryI`THzTDSTools`Jones`", {
    "JerryI`THzTDSTools`"
}]

JonesMatrix::usage = "JonesMatrix[sample, reference] constructs a frequency-dependent Jones matrix. A JonesMatrix object accepts a frequency Quantity for one matrix or a frequency-interval Quantity for its weighted average.";
OpticalAxisAnalysis::usage = "OpticalAxisAnalysis[matrix] analyzes a numeric 2x2 Jones matrix. OpticalAxisAnalysis[jones, interval] analyzes its weighted average over a frequency interval.";
PolarizationEllipse::usage = "PolarizationEllipse[matrix, scale] draws ellipses for a numeric 2x2 transmission matrix. PolarizationEllipse[jones, interval, scale] returns weighted ellipses for up to 100 frequencies.";
OpticalAxesTransmission::usage = "OpticalAxesTransmission[jones, opticalAxes] returns two {frequency, magnitude, phase} spectra using an OpticalAxisAnalysis association, two Jones vectors, or \"Circular\".";

Begin["`Private`"]

JonesMatrix::malformed = "Malformed Jones matrix: ``"
JonesMatrix::interval = "Expected a frequency interval containing sampled Jones matrices with positive total weight."
JonesMatrix[_] := (Message[JonesMatrix::malformed, "unknown pattern"]; $Failed)

validTF[trace_] := !FailureQ[TimeTraceQ[trace]]

SetAttributes[roundDeg, Listable];
roundDeg[a_] := With[{angle = Mod[Round[N[a], Pi/4], Pi]},
    Min[angle, Pi - angle]
];

MemberListQ[ar1_, ar2_] := Length[Complement[ar1, ar2]] === 0

JonesMatrix[sample:{
    ({_, _} -> _)..
}, reference: {
    ({_, _} -> _)..    
}] := Module[{s = <||>, r = <||>, temporal, data, freq, interval, weights},
    If[Length[sample] =!= Length[reference], 
        Message[JonesMatrix::malformed, "sample and reference arrays are not the same length"];
        Return[$Failed];
    ];
    If[!AllTrue[Join[sample[[All,2]], reference[[All,2]]], validTF], 
        Message[JonesMatrix::malformed, "provided spectra is not QuantityArray with {Time, 1} or {Freq, 1, 1} elements inside"];
        Return[$Failed];
    ];
    If[!AllTrue[N@Flatten[Join[sample[[All,1]], reference[[All,1]]]], NumberQ], 
        Message[JonesMatrix::malformed, "provided angles are non numerical"];
        Return[$Failed];
    ]; 
    If[Length[sample] != 4 || !MemberListQ[roundDeg/@{
        {0,0}, {90 Degree, 90 Degree}, {45 Degree, 90 Degree}, {45 Degree, 0 Degree}
    }, Map[roundDeg, sample[[All,1]]] ],
        Message[JonesMatrix::malformed, "missing angles: "<>ToString[Complement[roundDeg/@{
            {0,0}, {90 Degree, 90 Degree}, {45 Degree, 90 Degree}, {45 Degree, 0 Degree}
        }, Map[roundDeg, sample[[All,1]]] ]]];
        
        Return[$Failed];      
    ];
    If[Length[reference] != 4 || !MemberListQ[roundDeg/@{
        {0,0}, {90 Degree, 90 Degree}, {45 Degree, 90 Degree}, {45 Degree, 0 Degree}
    }, Map[roundDeg, reference[[All,1]]] ],
        Message[JonesMatrix::malformed, "missing angles: "<>ToString[Complement[roundDeg/@{
            {0,0}, {90 Degree, 90 Degree}, {45 Degree, 90 Degree}, {45 Degree, 0 Degree}
        }, Map[roundDeg, reference[[All,1]]] ]]];

        Return[$Failed];
      
    ];    

    temporal = Join[sample[[All,2]], reference[[All,2]]];
    If[AllTrue[temporal, TimeTraceQ],
        temporal = ResampleTraces[temporal];
    ];

    With[{v = roundDeg[#[[1]]]},
        s[v]  = If[TimeTraceQ[#[[2]]], TD2FD[#[[2]]], #[[2]]];
    ] &/@ Transpose[{ sample[[All,1]], Take[temporal, Length[sample]]  }];
    With[{v = roundDeg[#[[1]]]},
        r[v]  = If[TimeTraceQ[#[[2]]], TD2FD[#[[2]]], #[[2]]];
    ] &/@ Transpose[{ reference[[All,1]], Drop[temporal, Length[sample]] }];
  
    ClearAll[temporal];

    data = {
        s[{Pi/4, 0}],
        r[{Pi/4, 0}],
        r[{Pi/4, Pi/2}],
        s[{Pi/2, Pi/2}],
        r[{Pi/2, Pi/2}],
        s[{Pi/4, Pi/2}],
        s[{0, 0}],
        r[{0, 0}]
    };

    interval = GeometricMean[{
        QuantityMagnitude[FrequencyDomainConfidenceInterval[s[{0,0}]], "Terahertz"],
        QuantityMagnitude[FrequencyDomainConfidenceInterval[s[{Pi/2,Pi/2}]], "Terahertz"]
    }];

    data = Map[QuantityMagnitude[#, {"Terahertz", 1, 1}]&, data];
    
    freq = data[[1, All, 1]];
    weights = With[{max = Max[{
        data[[{1, 4, 6, 7}, All, 2]]
    }]},
        Transpose[{freq, (Table[Min[data[[{1, 4, 6, 7}, k, 2]]], {k, Length[data[[1]]]}]/max)^2}]
    ];
    
    data = Map[Map[Function[i, i[[2]] Exp[I i[[3]]]], #]&, data] // Transpose;

    ClearAll[s,r];

    With[{
        j = MapThread[{#1, jonesMatrix@@(#2)}&, {freq, data}], 
        i = interval, w = weights
    },
        ClearAll[freq, data, interval, weights];
        JonesMatrix[Interpolation[j, InterpolationOrder->0], "Terahertz", i, Interpolation[w, InterpolationOrder->0]]
    ]
]

JonesMatrix[f_InterpolatingFunction, units_String, interval_Interval, weights_InterpolatingFunction][int:Quantity[_Interval, _]] := With[{
    qi = Quiet @ Check[QuantityMagnitude[int, units], $Failed]
},
    If[!MatchQ[qi, Interval[{_?NumericQ, _?NumericQ}]],
        Message[JonesMatrix::interval]; $Failed,
        With[{
            bounds = First[qi]
        }, {
            arr = Select[extractIData[f], bounds[[1]] <= #[[1]] <= bounds[[2]] &]
        }, {
            w = weights /@ arr[[All, 1]]
        },
            If[arr === {} || Total[w] <= 0,
                Message[JonesMatrix::interval]; $Failed,
                Total[MapThread[#1 #2 &, {arr[[All, 2]], w}]]/Total[w]
            ]
        ]
    ]
]

JonesMatrix[f_InterpolatingFunction, units_String, interval_Interval, weights_InterpolatingFunction][pos_Quantity] := f[QuantityMagnitude[pos, units]]

OpticalAxisAnalysis::nojones = "Input data must be Jones matrix"
OpticalAxisAnalysis[j:{{_?NumberQ, _?NumberQ}, {_?NumberQ, _?NumberQ}}] := opticalAxisAnalysis[j]
OpticalAxisAnalysis[_] := (Message[OpticalAxisAnalysis::nojones]; $Failed);

OpticalAxisAnalysis[j: JonesMatrix[f_InterpolatingFunction, units_String, interval_Interval, weights_InterpolatingFunction]] := OpticalAxisAnalysis[j, Quantity[interval, units]]
OpticalAxisAnalysis[j: JonesMatrix[_InterpolatingFunction, _String, _Interval, _InterpolatingFunction], interval_Quantity] := With[{
    matrix = j[interval]
},
    If[matrix === $Failed, $Failed, opticalAxisAnalysis[matrix]]
]

OpticalAxesTransmission::polarizations = "Expected an OpticalAxisAnalysis association, two nonzero numeric Jones vectors, or \"Circular\"."
OpticalAxesTransmission[j_JonesMatrix, "Circular"] :=
    OpticalAxesTransmission[j, {{1, I}, {1, -I}}]

OpticalAxesTransmission[j_JonesMatrix, opticalAxes_Association] := With[{
    states = Lookup[opticalAxes, "Eigenpolarizations", {}]
},
    If[MatchQ[states, {_Association, _Association}] &&
       AllTrue[states, KeyExistsQ[#, "JonesVector"] &],
        OpticalAxesTransmission[j, states[[All, "JonesVector"]]],
        Message[OpticalAxesTransmission::polarizations]; $Failed
    ]
]

OpticalAxesTransmission[j: JonesMatrix[f_InterpolatingFunction, units_String, _Interval, _InterpolatingFunction],
                        opticalAxes:{{_?NumberQ, _?NumberQ}, {_?NumberQ, _?NumberQ}}] := With[{
    data = extractIData[f]
},
    If[AnyTrue[opticalAxes, TrueQ[Conjugate[#].# == 0] &],
        Message[OpticalAxesTransmission::polarizations]; $Failed,
        Function[vector,
            With[{
                transmission = Function[matrix,
                    Conjugate[vector].matrix.vector/(Conjugate[vector].vector)
                ] /@ data[[All, 2]]
            },
                QuantityArray[
                    Transpose[{data[[All, 1]], Abs[transmission], Arg[transmission]}],
                    {units, 1, 1}
                ]
            ]
        ] /@ opticalAxes
    ]
]

PolarizationEllipse[j: JonesMatrix[f_InterpolatingFunction, units_String, default_Interval, weights_InterpolatingFunction],
                    interval:(Automatic | _Quantity):Automatic,
                    Optional[scale_?NumericQ, 1],
                    opts:OptionsPattern[]] := With[{
    qi = QuantityMagnitude[If[interval === Automatic,
        Quantity[default, units], interval], units]
}, {
    arr = Select[extractIData[f], qi[[1,1]] <= #[[1]] <= qi[[1,2]] &]
},
    If[arr === {}, {},
        With[{sampled = If[Length[arr] > 100,
            arr[[Round[Subdivide[1, Length[arr], 99]]]], arr]},
            Function[row,
                Prepend[
                    PolarizationEllipse[row[[2]], scale, opts],
                    Opacity[Clip[weights[row[[1]]], {0, 1}]]
                ]
            ] /@ sampled
        ]
    ]
]

PolarizationEllipse[j:{{_?NumberQ, _?NumberQ}, {_?NumberQ, _?NumberQ}}, Optional[scale_?NumericQ, 1], OptionsPattern[]] := With[{
    analysis = opticalAxisAnalysis[j],
    style = OptionValue[PlotStyle]
}, {
    colors = Which[style===Automatic,
        eigenvalueColorData[analysis],

        MatchQ[{_, _}][style],
        style,

        True,
        {style, style}
    ]
},
{
    {colors[[1]], polarizationEllipse[analysis["Eigenpolarizations"][[1]], scale]},
    {colors[[2]], polarizationEllipse[analysis["Eigenpolarizations"][[2]], scale]}
}]

Options[PolarizationEllipse] = {PlotStyle->Automatic}

extractIData[ifun_] := Transpose[{Flatten[ifun["Grid"]], ifun["ValuesOnGrid"]}]


JonesMatrix /: MakeBoxes[obj: JonesMatrix[f_InterpolatingFunction, units_String, interval_Interval, weights_InterpolatingFunction], StandardForm] := Module[{
    above, below,
    report = OpticalAxisAnalysis[obj]
},
    above = { 
      {BoxForm`SummaryItem[{"Degeneracy: ", report["Degeneracy"]}]}
    };
    BoxForm`ArrangeSummaryBox[
       JonesMatrix, (* head *)
       obj,      (* interpretation *)
       Graphics[PolarizationEllipse[obj], AspectRatio->1],
       (* above and below must be in a format suitable for Grid or Column *)
       above,    (* always shown content *)
       Null (* expandable content. Currently not supported!*)
    ]
];

jonesMatrix[
  MM45_, RM45_, RP45_, M90_, R90_, MP45_, M0_, R0_
] := Module[
  {sampleTotal, referenceTotal},
  sampleTotal = {
    {M0, Sqrt[2] MP45 - M90},
    {Sqrt[2] MM45 - M0, M90}
  };
  referenceTotal = {
    {R0, Sqrt[2] RP45 - R90},
    {Sqrt[2] RM45 - R0, R90}
  };
  LinearSolve[referenceTotal, sampleTotal]
];

jonesPolarization[vector_] := Module[{
  z, s1, s2, s3, azimuth, ellipticityAngle, semiMajor, semiMinor
},
    z = N[(*FB[*)((vector)(*,*)/(*,*)((*BB[*)(Norm[vector])(*,*)(*"1:eJxTTMoPSmNmYGAo5gcSAUX5ZZkpqSn+BSWZ+XnFaYwgCRYg4ZGfk5LGBOLxAomwzNRyl9Tk/KLEkvyiYE6gSEhRYl5xQX5xajA7kKembGFkaGIN0QDSHlSakxoMNic1MQXM8MsvykWTZwUyXMtS80rAQn6lOTkAaNojwg=="*)(*]BB*)))(*]FB*)];
    s1 = (*SpB[*)Power[Abs[z[[1]]](*|*),(*|*)2](*]SpB*) - (*SpB[*)Power[Abs[z[[2]]](*|*),(*|*)2](*]SpB*);
    s2 = 2 Re[z[[1]] Conjugate[z[[2]]]];
    s3 = -2 Im[z[[1]] Conjugate[z[[2]]]];
    azimuth = Mod[((*FB[*)((1)(*,*)/(*,*)(2))(*]FB*)) ArcTan[s1, s2], \[Pi]];
    ellipticityAngle = ((*FB[*)((1)(*,*)/(*,*)(2))(*]FB*)) ArcSin[Clip[s3, {-1, 1}]];
    semiMajor = (*SqB[*)Sqrt[((*FB[*)((1)(*,*)/(*,*)(2))(*]FB*)) (1+(*SqB[*)Sqrt[(*SpB[*)Power[s1(*|*),(*|*)2](*]SpB*)+(*SpB[*)Power[s2(*|*),(*|*)2](*]SpB*)](*]SqB*))](*]SqB*);
    semiMinor = (*SqB[*)Sqrt[((*FB[*)((1)(*,*)/(*,*)(2))(*]FB*)) (1-(*SqB[*)Sqrt[(*SpB[*)Power[s1(*|*),(*|*)2](*]SpB*)+(*SpB[*)Power[s2(*|*),(*|*)2](*]SpB*)](*]SqB*))](*]SqB*);
    Association[
        "JonesVector" -> z
        ,
        "AzimuthRad" -> azimuth
        ,
        "AzimuthDeg" -> (*FB[*)((azimuth 180)(*,*)/(*,*)(\[Pi]))(*]FB*)
        ,
        "EllipticityAngleRad" -> ellipticityAngle
        ,
        "EllipticityAngleDeg" -> (*FB[*)((ellipticityAngle 180)(*,*)/(*,*)(\[Pi]))(*]FB*)
        ,
        "MinorMajorRatio" -> Abs[Tan[ellipticityAngle]]
        ,
        "SemiMajorAxis" -> semiMajor
        ,
        "SemiMinorAxis" -> semiMinor
    ]
];

opticalAxisAnalysis[jones_] := Module[{eigenvalues, eigenvectors, eigenstates, u, w, v, separation},
    {eigenvalues, eigenvectors} = Eigensystem[jones];
    eigenstates = jonesPolarization /@ eigenvectors;
    {u, w, v} = SingularValueDecomposition[jones];
    separation = Abs[eigenstates[[1]]["AzimuthDeg"] - eigenstates[[2]]["AzimuthDeg"]];
    separation = Min[separation, 180 - separation];
    Association[
        "Eigenvalues" -> eigenvalues,
        "Eigenvectors" -> eigenstates[[All, "JonesVector"]],
        "LinearizedEigenvectors" -> ({Cos[#["AzimuthRad"]], Sin[#["AzimuthRad"]]} & /@ eigenstates),
        "Eigenpolarizations" -> eigenstates
        ,
        "EigenaxisSeparationDeg" -> separation
        ,
        "SingularValues" -> Diagonal[w]
        ,
        "InputPrincipalStates" -> jonesPolarization /@ ((*BB[*)(Transpose[v])(*,*)(*"1:eJxTTMoPSmNmYGAo5gcSAUX5ZZkpqSn+BSWZ+XnFaYwgCRYg4ZGfk5LGBOLxAomwzNRyl9Tk/KLEkvyiYE6gSEhRYl5xQX5xajBISwhEKUhjUGlOajDYhNTEFFSlaIpYgQzXstS8ErCQX2lODgB64iTl"*)(*]BB*))
        ,
        "OutputPrincipalStates" -> jonesPolarization /@ ((*BB[*)(Transpose[u])(*,*)(*"1:eJxTTMoPSmNmYGAo5gcSAUX5ZZkpqSn+BSWZ+XnFaYwgCRYg4ZGfk5LGBOLxAomwzNRyl9Tk/KLEkvyiYE6gSEhRYl5xQX5xajBISwhEKUhjUGlOajDYhNTEFFSlaIpYgQzXstS8ErCQX2lODgB64iTl"*)(*]BB*))
        ,
        "RelativeAsymmetry" -> (*FB[*)(((*BB[*)(Norm[jones-(*BB[*)(Transpose[jones])(*,*)(*"1:eJxTTMoPSmNmYGAo5gcSAUX5ZZkpqSn+BSWZ+XnFaYwgCRYg4ZGfk5LGBOLxAomwzNRyl9Tk/KLEkvyiYE6gSEhRYl5xQX5xajBISwhEKUhjUGlOajDYhNTEFFSlaIpYgQzXstS8ErCQX2lODgB64iTl"*)(*]BB*)])(*,*)(*"1:eJxTTMoPSmNmYGAo5gcSAUX5ZZkpqSn+BSWZ+XnFaYwgCRYg4ZGfk5LGBOLxAomwzNRyl9Tk/KLEkvyiYE6gSEhRYl5xQX5xajA7kKembGFkaGIN0QDSHlSakxoMNic1MQXM8MsvykWTZwUyXMtS80rAQn6lOTkAaNojwg=="*)(*]BB*))(*,*)/(*,*)((*BB[*)(Norm[jones])(*,*)(*"1:eJxTTMoPSmNmYGAo5gcSAUX5ZZkpqSn+BSWZ+XnFaYwgCRYg4ZGfk5LGBOLxAomwzNRyl9Tk/KLEkvyiYE6gSEhRYl5xQX5xajA7kKembGFkaGIN0QDSHlSakxoMNic1MQXM8MsvykWTZwUyXMtS80rAQn6lOTkAaNojwg=="*)(*]BB*)))(*]FB*)
        ,
        "RelativeNonNormality" -> (*FB[*)(((*BB[*)(Norm[jones.((*BB[*)(ConjugateTranspose[jones])(*,*)(*"1:eJxTTMoPSmNmYGAo5gcSAUX5ZZkpqSn+BSWZ+XnFaYwgCRYg4ZGfk5LGBOLxAomwzNRyl9Tk/KLEkvyiYE6gSEhRYl5xQX5xajAHkKeWkpienlpkDdEB0h9UmpMaDDYoNTElWAjIcM7PyypNTyxJhWtFU80KZLiWpeaVgIX8SnNyAHmFKww="*)(*]BB*))-((*BB[*)(ConjugateTranspose[jones])(*,*)(*"1:eJxTTMoPSmNmYGAo5gcSAUX5ZZkpqSn+BSWZ+XnFaYwgCRYg4ZGfk5LGBOLxAomwzNRyl9Tk/KLEkvyiYE6gSEhRYl5xQX5xajAHkKeWkpienlpkDdEB0h9UmpMaDDYoNTElWAjIcM7PyypNTyxJhWtFU80KZLiWpeaVgIX8SnNyAHmFKww="*)(*]BB*)).jones])(*,*)(*"1:eJxTTMoPSmNmYGAo5gcSAUX5ZZkpqSn+BSWZ+XnFaYwgCRYg4ZGfk5LGBOLxAomwzNRyl9Tk/KLEkvyiYE6gSEhRYl5xQX5xajA7kKembGFkaGIN0QDSHlSakxoMNic1MQXM8MsvykWTZwUyXMtS80rAQn6lOTkAaNojwg=="*)(*]BB*))(*,*)/(*,*)((*SpB[*)Power[((*BB[*)(Norm[jones])(*,*)(*"1:eJxTTMoPSmNmYGAo5gcSAUX5ZZkpqSn+BSWZ+XnFaYwgCRYg4ZGfk5LGBOLxAomwzNRyl9Tk/KLEkvyiYE6gSEhRYl5xQX5xajA7kKembGFkaGIN0QDSHlSakxoMNic1MQXM8MsvykWTZwUyXMtS80rAQn6lOTkAaNojwg=="*)(*]BB*))(*|*),(*|*)2](*]SpB*)))(*]FB*),
        "Degeneracy" -> 1 - Clip[
          Abs[eigenvalues[[1]] - eigenvalues[[2]]]/
            (Abs[eigenvalues[[1]]] + Abs[eigenvalues[[2]]] + $MachineEpsilon),
          {0, 1}
        ]
    ]
];

polarizationEllipse[state_Association, scale_:1] := Rotate[Circle[{0, 0}, scale Re@{state["SemiMajorAxis"], state["SemiMinorAxis"]}], state["AzimuthRad"]//Re, {0, 0}]

eigenvalueColorData[axisAnalysis_Association] := Module[
  {values, separation, degeneracy, firstColor, secondColor},
  values = axisAnalysis["Eigenvalues"];
  separation = Clip[
    Abs[values[[1]] - values[[2]]]/
      (Abs[values[[1]]] + Abs[values[[2]]] + $MachineEpsilon),
    {0, 1}
  ];
  degeneracy = 1 - separation;
  firstColor = Blend[{Red, RGBColor[1,0.5,1]}, degeneracy];
  secondColor = Blend[{Blue, RGBColor[1,0.5,1]}, degeneracy];
  {firstColor, secondColor}
]


End[]
EndPackage[]
