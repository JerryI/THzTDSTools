(* Run from the project root with TestReport["Tests/Jones.wlt"]. *)

testRoot = SelectFirst[
    {Directory[], DirectoryName[Directory[]]},
    FileExistsQ[FileNameJoin[{#, "Kernel", "Jones.wl"}]] &,
    $Failed
];
If[testRoot === $Failed, Print["Run Tests/Jones.wlt from the project root or Tests directory."]; Abort[]];
Get[FileNameJoin[{testRoot, "Kernel", "Common.wl"}]];
Get[FileNameJoin[{testRoot, "Kernel", "Kernel.wl"}]];
Get[FileNameJoin[{testRoot, "Kernel", "Jones.wl"}]];

ClearAll[jonesNearQ];
jonesNearQ[a_, b_, tolerance_:10^-10] :=
    Max[Abs[Flatten[{N[a - b]}]]] < tolerance;

testJonesFrequencies = N[Range[0, 10]/10];
testJonesMatrices = Table[
    {f, {{2 + f, 0.}, {0., 1.}}}, {f, testJonesFrequencies}
];
testJonesWeights = Table[
    {f, 0.2 + 0.5 f}, {f, testJonesFrequencies}
];
testJonesWeightFunction = Interpolation[
    testJonesWeights, InterpolationOrder -> 0
];
testJones = JerryI`THzTDSTools`Jones`JonesMatrix[
    Interpolation[testJonesMatrices, InterpolationOrder -> 0],
    "Terahertz", Interval[{0.2, 0.8}], testJonesWeightFunction
];

VerificationTest[
    {
        JerryI`THzTDSTools`Jones`Private`roundDeg[
            {0, Pi/4, Pi/2, 3 Pi/4, Pi, 2 Pi, -Pi/4, -Pi/2, -Pi}],
        JerryI`THzTDSTools`Jones`Private`roundDeg[
            {{Pi, 0}, {-Pi/2, 2 Pi}}]
    },
    {{0, Pi/4, Pi/2, Pi/4, 0, 0, Pi/4, Pi/2, 0},
     {{0, 0}, {Pi/2, 0}}},
    TestID -> "Analyzer angles fold into zero through Pi over two"
]

VerificationTest[
    With[{
        analysis = JerryI`THzTDSTools`Jones`OpticalAxisAnalysis[testJones],
        narrow = JerryI`THzTDSTools`Jones`OpticalAxisAnalysis[
            testJones, Quantity[Interval[{0.2, 0.4}], "Terahertz"]],
        f = N[Range[2, 8]/10],
        fNarrow = N[Range[2, 4]/10]
    },
        AssociationQ[analysis] &&
        KeyExistsQ[analysis, "Eigenpolarizations"] &&
        jonesNearQ[Sort[analysis["Eigenvalues"]],
            {1., 2 + Total[f (0.2 + 0.5 f)]/Total[0.2 + 0.5 f]}] &&
        jonesNearQ[Sort[narrow["Eigenvalues"]],
            {1., 2 + Total[fNarrow (0.2 + 0.5 fNarrow)]/
                Total[0.2 + 0.5 fNarrow]}]
    ],
    True,
    TestID -> "Optical axis analysis uses selected weighted Jones matrices"
]

VerificationTest[
    Module[{matrixFunction, jones, spectra, axes, associated,
            complexAxes, complexData, data},
        matrixFunction = Interpolation[
            Table[{f, DiagonalMatrix[
                {2 Exp[0.3 I f], 0.7 Exp[-0.2 I f]}]},
                {f, testJonesFrequencies}],
            InterpolationOrder -> 0
        ];
        jones = JerryI`THzTDSTools`Jones`JonesMatrix[
            matrixFunction, "Terahertz", Interval[{0.2, 0.8}],
            testJonesWeightFunction
        ];
        spectra = JerryI`THzTDSTools`Jones`OpticalAxesTransmission[
            jones, {{1, 0}, {0, 1}}];
        axes = JerryI`THzTDSTools`Jones`OpticalAxisAnalysis[jones];
        associated = JerryI`THzTDSTools`Jones`OpticalAxesTransmission[
            jones, axes];
        complexAxes = JerryI`THzTDSTools`Jones`OpticalAxesTransmission[
            jones, {{1, I}, {0, 1}}];
        complexData = QuantityMagnitude[complexAxes[[1]],
            {"Terahertz", 1, 1}];
        data = QuantityMagnitude[#, {"Terahertz", 1, 1}] & /@ spectra;
        Dimensions /@ spectra === {{11, 3}, {11, 3}} &&
            jonesNearQ[data[[1, All, 2]], ConstantArray[2., 11]] &&
            jonesNearQ[data[[1, All, 3]], 0.3 testJonesFrequencies] &&
            jonesNearQ[data[[2, All, 2]], ConstantArray[0.7, 11]] &&
            jonesNearQ[data[[2, All, 3]], -0.2 testJonesFrequencies] &&
            Dimensions /@ associated === {{11, 3}, {11, 3}} &&
            jonesNearQ[
                Sort[QuantityMagnitude[#, {"Terahertz", 1, 1}][[6, 2]] & /@
                    associated], {0.7, 2.}] &&
            jonesNearQ[
                complexData[[All, 2]] Exp[I complexData[[All, 3]]],
                (2 Exp[0.3 I testJonesFrequencies] +
                    0.7 Exp[-0.2 I testJonesFrequencies])/2]
    ],
    True,
    TestID -> "Optical axes transmission returns two Euler spectra"
]

VerificationTest[
    Module[{selected, unscaled, styled, scaledRadius, baseRadius,
            denseFrequencies, dense, ellipses},
        selected = JerryI`THzTDSTools`Jones`PolarizationEllipse[
            testJones, Quantity[Interval[{0.3, 0.5}], "Terahertz"],
            2, PlotStyle -> {Green, Orange}];
        unscaled = JerryI`THzTDSTools`Jones`PolarizationEllipse[
            testJones, Quantity[Interval[{0.3, 0.5}], "Terahertz"],
            1, PlotStyle -> {Green, Orange}];
        styled = JerryI`THzTDSTools`Jones`PolarizationEllipse[
            testJones, PlotStyle -> Red];
        scaledRadius = First[Cases[selected[[1, 2]],
            Circle[_, radii_] :> radii, Infinity]];
        baseRadius = First[Cases[unscaled[[1, 2]],
            Circle[_, radii_] :> radii, Infinity]];
        denseFrequencies = N[Range[0, 1000]/1000];
        dense = JerryI`THzTDSTools`Jones`JonesMatrix[
            Interpolation[
                Table[{f, {{2., 0.}, {0., 1.}}},
                    {f, denseFrequencies}],
                InterpolationOrder -> 0],
            "Terahertz", Interval[{0., 1.}],
            Interpolation[
                Table[{f, 0.5}, {f, denseFrequencies}],
                InterpolationOrder -> 0]
        ];
        ellipses = JerryI`THzTDSTools`Jones`PolarizationEllipse[dense];
        Length[selected] == 3 &&
            jonesNearQ[selected[[1, 1, 1]], 0.35] &&
            selected[[1, 2, 1]] === Green &&
            selected[[1, 3, 1]] === Orange &&
            styled[[1, 2, 1]] === Red &&
            styled[[1, 3, 1]] === Red &&
            jonesNearQ[scaledRadius, 2 baseRadius] &&
            Length[ellipses] == 500 &&
            Length[Cases[ellipses, _Opacity, Infinity]] == Length[ellipses]
    ],
    True,
    TestID -> "Jones ellipses use weights and at most 500 frequencies"
]
