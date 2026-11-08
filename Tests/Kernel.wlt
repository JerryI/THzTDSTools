(* Run from the project root with TestReport["Tests/Kernel.wlt"]. *)

testRoot = SelectFirst[
    {Directory[], DirectoryName[Directory[]]},
    FileExistsQ[FileNameJoin[{#, "Kernel", "Kernel.wl"}]] &,
    $Failed
];
If[testRoot === $Failed, Print["Run Tests/Kernel.wlt from the project root or Tests directory."]; Abort[]];
Get[FileNameJoin[{testRoot, "Kernel", "Common.wl"}]];
Get[FileNameJoin[{testRoot, "Kernel", "Kernel.wl"}]];

ClearAll[testNearQ, testComplexSpectrum];
testNearQ[a_, b_, tolerance_:10^-10] :=
    Max[Abs[Flatten[{N[a - b]}]]] < tolerance;
testComplexSpectrum[spectrum_QuantityArray] := Module[{data},
    data = QuantityMagnitude[spectrum, {"Terahertz", 1, 1}];
    data[[All, 2]] Exp[I data[[All, 3]]]
];

testTime = QuantityArray[
    N[Table[{n/10, Sin[2 Pi 4 n/64] + 0.2 Cos[2 Pi 7 n/64]}, {n, 0, 63}]],
    {"Picoseconds", 1}
];
testSpectrum = JerryI`THzTDSTools`TD2FD[testTime];
testPolar = QuantityMagnitude[testSpectrum, {"Terahertz", 1, 1}];
testInterval = Interval[{Quantity[0.1, "Terahertz"], Quantity[2., "Terahertz"]}];

VerificationTest[
    Dimensions[testSpectrum],
    {32, 3},
    TestID -> "TD2FD returns a three-column Euler spectrum"
]

VerificationTest[
    testNearQ[
        testComplexSpectrum[testSpectrum],
        Take[Fourier[QuantityMagnitude[testTime, {"Picoseconds", 1}][[All, 2]]], 32]
    ],
    True,
    TestID -> "Magnitude and phase reconstruct the complex Fourier coefficients"
]

VerificationTest[
    testNearQ[
        QuantityMagnitude[JerryI`THzTDSTools`FD2TD[testSpectrum], {"Picoseconds", 1}][[All, 2]],
        QuantityMagnitude[testTime, {"Picoseconds", 1}][[All, 2]]
    ],
    True,
    TestID -> "TD2FD and FD2TD round trip a band-limited trace"
]

VerificationTest[
    Module[{trace, spectrum, restored},
        trace = QuantityArray[N[Table[{n/10, (-1)^n}, {n, 0, 63}]], {"Picoseconds", 1}];
        spectrum = JerryI`THzTDSTools`TD2FD[trace, "IncludeNyquist" -> True];
        restored = JerryI`THzTDSTools`FD2TD[spectrum, "IncludeNyquist" -> True];
        Dimensions[spectrum] === {33, 3} && testNearQ[
            QuantityMagnitude[restored, {"Picoseconds", 1}][[All, 2]],
            QuantityMagnitude[trace, {"Picoseconds", 1}][[All, 2]]
        ]
    ],
    True,
    TestID -> "Nyquist opt-in preserves an alternating signal"
]

VerificationTest[
    Module[{trace, restored},
        trace = QuantityArray[
            N[Table[{n/10, Sin[2 Pi 5 n/63]}, {n, 0, 62}]],
            {"Picoseconds", 1}
        ];
        restored = JerryI`THzTDSTools`FD2TD[
            JerryI`THzTDSTools`TD2FD[trace, "IncludeNyquist" -> True],
            "IncludeNyquist" -> True
        ];
        Dimensions[restored] === {64, 2} && testNearQ[
            QuantityMagnitude[restored, {"Picoseconds", 1}][[1 ;; 63, 2]],
            QuantityMagnitude[trace, {"Picoseconds", 1}][[All, 2]]
        ]
    ],
    True,
    TestID -> "Odd sample count is padded once for Nyquist opt-in"
]

testUnwrappedData = testPolar;
testUnwrappedData[[All, 3]] += 4 Pi;
testUnwrappedSpectrum = QuantityArray[testUnwrappedData, {"Terahertz", 1, 1}];

VerificationTest[
    testNearQ[
        QuantityMagnitude[testUnwrappedSpectrum, {"Terahertz", 1, 1}][[All, 3]],
        testUnwrappedData[[All, 3]]
    ] && testNearQ[
        QuantityMagnitude[JerryI`THzTDSTools`FD2TD[testUnwrappedSpectrum],
            {"Picoseconds", 1}][[All, 2]],
        QuantityMagnitude[JerryI`THzTDSTools`FD2TD[testSpectrum],
            {"Picoseconds", 1}][[All, 2]]
    ],
    True,
    TestID -> "Stored unwrapped phase survives and gives the same inverse transform"
]

VerificationTest[
    testNearQ[
        QuantityMagnitude[
            JerryI`THzTDSTools`UnwrapPhase[testUnwrappedSpectrum, testInterval, Method -> None],
            {"Terahertz", 1, 1}
        ],
        testUnwrappedData
    ],
    True,
    TestID -> "UnwrapPhase with Method None preserves the complete input spectrum"
]

testGaussian = QuantityArray[
    N[Table[{n/10, Exp[-(n - 30)^2/20]}, {n, 0, 127}]],
    {"Picoseconds", 1}
];
testGaussianSpectrum = JerryI`THzTDSTools`TD2FD[testGaussian];
testPhaseInterval = Interval[{Quantity[0.1, "Terahertz"], Quantity[1.2, "Terahertz"]}];

VerificationTest[
    Module[{before, after, inside, outside},
        before = QuantityMagnitude[testGaussianSpectrum, {"Terahertz", 1, 1}];
        after = QuantityMagnitude[
            JerryI`THzTDSTools`UnwrapPhase[testGaussian, testPhaseInterval,
                Method -> "Basic", "RemoveLinearPhase" -> True],
            {"Terahertz", 1, 1}
        ];
        inside = Select[after, 0.1 <= #[[1]] <= 1.2 &];
        outside = Select[Range[Length[after]],
            after[[#, 1]] < 0.1 || after[[#, 1]] > 1.2 &];
        Dimensions[after] === Dimensions[before] &&
            Max[Abs[inside[[All, 3]]]] < 10^-8 &&
            testNearQ[after[[All, 2]], before[[All, 2]]] &&
            testNearQ[after[[outside]], before[[outside]]]
    ],
    True,
    TestID -> "Basic phase removes pulse arrival and preserves other bins"
]

VerificationTest[
    And @@ Flatten[Table[
        Dimensions[JerryI`THzTDSTools`UnwrapPhase[input, testPhaseInterval,
            Method -> method]] === Dimensions[testGaussianSpectrum],
        {input, {testGaussian, testGaussianSpectrum}},
        {method, {None, "Basic", "KK", "KKGuided", "KKGuidedMultiscale"}}
    ]],
    True,
    TestID -> "All phase methods accept both time and Euler spectrum inputs"
]

VerificationTest[
    JerryI`THzTDSTools`FrequencyDomainConfidenceInterval[testGaussian] ===
        JerryI`THzTDSTools`FrequencyDomainConfidenceInterval[testGaussianSpectrum],
    True,
    TestID -> "Confidence interval agrees between time and frequency inputs"
]

VerificationTest[
    Module[{trimmed, expected},
        trimmed = QuantityMagnitude[
            JerryI`THzTDSTools`FrequencySpectrumTrim[testUnwrappedSpectrum,
                testInterval],
            {"Terahertz", 1, 1}
        ];
        expected = Select[testUnwrappedData, 0.1 <= #[[1]] <= 2. &];
        testNearQ[trimmed, expected]
    ],
    True,
    TestID -> "Trimming preserves magnitude and unwrapped phase"
]

VerificationTest[
    Module[{interval, trimmed},
        interval = Interval[{
            Quantity[10, 1/"Centimeters"], Quantity[20, 1/"Centimeters"]
        }];
        trimmed = QuantityMagnitude[
            JerryI`THzTDSTools`FrequencySpectrumTrim[testSpectrum, interval],
            {"Terahertz", 1, 1}
        ];
        testNearQ[trimmed[[All, 1]], {0.3125, 0.46875}]
    ],
    True,
    TestID -> "Trimming accepts inverse-centimeter bounds"
]

VerificationTest[
    {JerryI`THzTDSTools`TimeTraceQ[testTime],
     JerryI`THzTDSTools`TimeTraceQ[testSpectrum],
     JerryI`THzTDSTools`TimeTraceQ[
         QuantityArray[{{0., 1.}, {1., 2.}}, {"Terahertz", 1}]]},
    {True, False, $Failed},
    TestID -> "TimeTraceQ distinguishes time, Euler spectrum, and old format"
]

VerificationTest[
    Module[{delayed, thickness},
        delayed = QuantityArray[
            N[Table[{n/10, Boole[n == 5]}, {n, 0, 63}]],
            {"Picoseconds", 1}
        ];
        thickness = Quantity[1, "Millimeters"];
        testNearQ[
            JerryI`THzTDSTools`AveragedRefrectiveIndex[delayed, thickness],
            1.149896229
        ] && testNearQ[
            JerryI`THzTDSTools`AveragedRefrectiveIndex[
                JerryI`THzTDSTools`TD2FD[delayed], thickness],
            1.149896229
        ]
    ],
    True,
    TestID -> "Average index agrees between a delayed trace and its spectrum"
]

VerificationTest[
    Module[{reference, sample, transfer, refSpectrum, samSpectrum, expected},
        reference = QuantityArray[
            N[Table[{n/10, Boole[n == 10]}, {n, 0, 63}]],
            {"Picoseconds", 1}
        ];
        sample = QuantityArray[
            N[Table[{n/10, 0.8 Boole[n == 15]}, {n, 0, 63}]],
            {"Picoseconds", 1}
        ];
        transfer = testComplexSpectrum[JerryI`THzTDSTools`TD2FD[
            JerryI`THzTDSTools`DeconvolveTraces[sample, reference]]];
        refSpectrum = testComplexSpectrum[JerryI`THzTDSTools`TD2FD[reference]];
        samSpectrum = testComplexSpectrum[JerryI`THzTDSTools`TD2FD[sample]];
        expected = Conjugate[refSpectrum] samSpectrum/
            (Abs[refSpectrum]^2 + 10.^-5 Max[Abs[refSpectrum]^2]);
        testNearQ[Rest[transfer], Rest[expected]]
    ],
    True,
    TestID -> "Deconvolution spectrum is the regularized complex transmission"
]

VerificationTest[
    Module[{first, second, output},
        first = QuantityArray[
            N[Table[{n/10, Sin[2 Pi 4 n/64]}, {n, 0, 63}]],
            {"Picoseconds", 1}
        ];
        second = QuantityArray[
            N[Table[{n/10, Cos[2 Pi 4 n/64]}, {n, 0, 63}]],
            {"Picoseconds", 1}
        ];
        output = JerryI`THzTDSTools`LinearProjectionsToCircular[
            first, second, {0, 90 Degree}];
        Dimensions[output] === {2, 64, 2}
    ],
    True,
    TestID -> "Linear projections transform with arbitrary independent angles"
]

VerificationTest[
    Module[{c = 0.0299792458, d = 0.1, n = 1.3, k = 0.02,
            data, spectrum, thickness, estimatedN, estimatedK},
        data = N[Table[
            {j/20, (4 n/(n + 1)^2) Exp[-2 Pi (j/20) d k/c],
                2 Pi (j/20) d (n - 1)/c},
            {j, 0, 511}
        ]];
        spectrum = QuantityArray[data, {"Terahertz", 1, 1}];
        thickness = Quantity[d, "Centimeters"];
        estimatedN = QuantityMagnitude[
            JerryI`THzTDSTools`EstimateMaterialParameter[
                "n", spectrum, thickness, Method -> "Iterative",
                "Target" -> "CPU", "MovingAverage" -> False],
            {"Terahertz", 1}];
        estimatedK = QuantityMagnitude[
            JerryI`THzTDSTools`EstimateMaterialParameter[
                "\[Kappa]", spectrum, thickness, Method -> "Iterative",
                "Target" -> "CPU", "MovingAverage" -> False],
            {"Terahertz", 1}];
        Abs[estimatedN[[257, 2]] - n] < 5 10^-4 &&
            Abs[estimatedK[[257, 2]] - k] < 5 10^-4
    ],
    True,
    TestID -> "Iterative CPU solver handles spectra beyond 256 bins"
]

VerificationTest[
    Module[{c = 0.0299792458, d = 0.1, n = 1.3, k = 0.02,
            nc, fresnel, model, spectrum, thickness, corrected, timeCorrected,
            estimatedN, band, before, after},
        nc = n + I k;
        fresnel = 4 nc/(1 + nc)^2;
        model = N[Table[
            With[{f = j/20},
                With[{
                    main = fresnel Exp[I 2 Pi f d (nc - 1)/c],
                    fp = 1/(1 - ((nc - 1)/(nc + 1))^2
                        Exp[I 4 Pi f d nc/c])
                },
                    {f, Abs[main fp],
                        2 Pi f d (n - 1)/c + Arg[fresnel] + Arg[fp],
                        Abs[main]}
                ]
            ],
            {j, 0, 127}
        ]];
        spectrum = QuantityArray[model[[All, 1 ;; 3]],
            {"Terahertz", 1, 1}];
        thickness = Quantity[d, "Centimeters"];
        corrected = QuantityMagnitude[
            JerryI`THzTDSTools`DeconvolveFP[
                spectrum, thickness, "Target" -> "CPU",
                "MovingAverage" -> False],
            {"Terahertz", 1, 1}];
        timeCorrected = JerryI`THzTDSTools`DeconvolveFP[
            JerryI`THzTDSTools`FD2TD[spectrum],
            thickness, "Target" -> "CPU", "MovingAverage" -> False];
        estimatedN = QuantityMagnitude[
            JerryI`THzTDSTools`EstimateMaterialParameter[
                "n", spectrum, thickness, Method -> "Iterative",
                "SolveFP" -> True, "Target" -> "CPU",
                "MovingAverage" -> False],
            {"Terahertz", 1}];
        band = Select[Range[Length[corrected]],
            0.2 <= corrected[[#, 1]] <= 2.5 &];
        before = Max[Abs[model[[band, 2]] - model[[band, 4]]]];
        after = Max[Abs[corrected[[band, 2]] - model[[band, 4]]]];
        after < before/10 &&
            Dimensions[timeCorrected] === {256, 2} &&
            Abs[estimatedN[[21, 2]] - n] < 5 10^-4
    ],
    True,
    TestID -> "Fabry-Perot removal and iterative index recovery"
]

VerificationTest[
    Module[{c = 0.0299792458, d = 0.1, n = 1.8, k = 0.01,
            nc, fresnel, model, interval, thickness, withFP, withoutFP},
        nc = n + I k;
        fresnel = 4 nc/(1 + nc)^2;
        model = N @ Table[
            With[{f = j/100.},
                With[{main = fresnel Exp[-2 Pi f d k/c],
                      phase = 2 Pi f d (n - 1)/c + Arg[fresnel],
                      fp = 1/(1 - ((nc - 1)/(nc + 1))^2
                          Exp[I 4 Pi f d nc/c])},
                    {f, Abs[main fp], phase + Arg[fp],
                        Abs[main], phase}
                ]
            ],
            {j, 0, 511}
        ];
        interval = Interval[{
            Quantity[0.2, "Terahertz"], Quantity[4., "Terahertz"]
        }];
        thickness = Quantity[d, "Centimeters"];
        withFP = JerryI`THzTDSTools`EstimateMaterialParameter[
            "FP Reduction",
            QuantityArray[model[[All, 1 ;; 3]], {"Terahertz", 1, 1}],
            thickness, "Target" -> "CPU", "FPInterval" -> interval
        ];
        withoutFP = JerryI`THzTDSTools`EstimateMaterialParameter[
            "FP Reduction",
            QuantityArray[model[[All, {1, 4, 5}]], {"Terahertz", 1, 1}],
            thickness, "Target" -> "CPU", "FPInterval" -> interval
        ];
        NumericQ[withFP] && withFP > 2 && withoutFP === 1.
    ],
    True,
    TestID -> "FP reduction score detects a removable echo and ignores a smooth spectrum"
]

VerificationTest[
    Module[{c = 0.0299792458, d = 0.0635, index = 1.8 + 0.01 I,
            indexFunction, spectrum, data, frequencies, expected},
        indexFunction[f_?NumericQ] := index;
        spectrum = JerryI`THzTDSTools`GenerateMaterialResponse[
            indexFunction, Quantity[d, "Centimeters"],
            {Quantity[Interval[{0, 2}], "Terahertz"],
             Quantity[0.025, "Terahertz"]}];
        data = QuantityMagnitude[spectrum, {"Terahertz", 1, 1}];
        frequencies = data[[All, 1]];
        expected = (4 index/(index + 1)^2)
            Exp[I 2 Pi frequencies d (index - 1)/c];
        Dimensions[spectrum] === {81, 3} &&
            testNearQ[data[[All, 2]] Exp[I data[[All, 3]]], expected] &&
            Last[data][[3]] > 2 Pi
    ],
    True,
    TestID -> "Generated slab spectrum matches the single-pass transfer and retains delay phase"
]

VerificationTest[
    Module[{c = 0.0299792458, d = 0.08, index = 1.6 + 0.02 I,
            spectrum, data, frequencies, expected},
        spectrum = JerryI`THzTDSTools`GenerateMaterialResponse[
            Function[f, index], Quantity[d, "Centimeters"],
            {Quantity[Interval[{0, 2}], "Terahertz"],
             Quantity[0.05, "Terahertz"]}, "SolveFP" -> True];
        data = QuantityMagnitude[spectrum, {"Terahertz", 1, 1}];
        frequencies = data[[All, 1]];
        expected = (4 index/(index + 1)^2)
            Exp[I 2 Pi frequencies d (index - 1)/c]/
            (1 - ((index - 1)/(index + 1))^2
                Exp[I 4 Pi frequencies d index/c]);
        testNearQ[data[[All, 2]] Exp[I data[[All, 3]]], expected]
    ],
    True,
    TestID -> "Generated FP spectrum matches the slab multiple-reflection denominator"
]

VerificationTest[
    Module[{reference, model, output, bare, refData, modelData, outData},
        reference = QuantityArray[
            {{0., 2., 4 Pi}, {0.25, 0.4, 4 Pi + 0.3},
             {0.5, 0.6, 4 Pi + 0.8}, {0.75, 0.8, 4 Pi + 1.1},
             {1., 0.3, 4 Pi + 1.4}},
            {"Terahertz", 1, 1}];
        model = JerryI`THzTDSTools`GenerateMaterialResponse[
            Function[f, 1.5 + 0.01 I], Quantity[0.2, "Millimeters"],
            {Quantity[Interval[{0, 1}], "Terahertz"],
             Quantity[0.25, "Terahertz"]}];
        output = JerryI`THzTDSTools`GenerateMaterialResponse[
            Function[f, 1.5 + 0.01 I], Quantity[0.2, "Millimeters"],
            reference];
        bare = JerryI`THzTDSTools`GenerateMaterialResponse[
            Function[f, 1.5 + 0.01 I], Quantity[0.2, "Millimeters"],
            reference, "ConvolveWithReference" -> False];
        refData = QuantityMagnitude[reference, {"Terahertz", 1, 1}];
        modelData = QuantityMagnitude[model, {"Terahertz", 1, 1}];
        outData = QuantityMagnitude[output, {"Terahertz", 1, 1}];
        testNearQ[outData[[All, 2]], refData[[All, 2]] modelData[[All, 2]]] &&
            testNearQ[outData[[All, 3]], refData[[All, 3]] + modelData[[All, 3]]] &&
            testNearQ[QuantityMagnitude[bare, {"Terahertz", 1, 1}], modelData]
    ],
    True,
    TestID -> "Frequency reference is multiplied without losing its unwrapped phase"
]

VerificationTest[
    Module[{reference, identity, output, impulse, expectedSpectrum,
            actualSpectrum, bareSpectrum,
            thickness = Quantity[0.1, "Millimeters"]},
        reference = QuantityArray[
            N @ Table[{3 + j/10, Exp[-(j - 12)^2/8] + 0.1 (-1)^j},
                {j, 0, 63}], {"Picoseconds", 1}];
        identity = JerryI`THzTDSTools`GenerateMaterialResponse[
            Function[f, 1.], thickness, reference];
        output = JerryI`THzTDSTools`GenerateMaterialResponse[
            Function[f, 1.4 + 0.01 I], thickness, reference];
        impulse = JerryI`THzTDSTools`GenerateMaterialResponse[
            Function[f, 1.4 + 0.01 I], thickness, reference,
            "ConvolveWithReference" -> False];
        expectedSpectrum = JerryI`THzTDSTools`GenerateMaterialResponse[
            Function[f, 1.4 + 0.01 I], thickness,
            JerryI`THzTDSTools`TD2FD[reference, "IncludeNyquist" -> True]];
        bareSpectrum = JerryI`THzTDSTools`GenerateMaterialResponse[
            Function[f, 1.4 + 0.01 I], thickness,
            JerryI`THzTDSTools`TD2FD[reference, "IncludeNyquist" -> True],
            "ConvolveWithReference" -> False];
        actualSpectrum = JerryI`THzTDSTools`TD2FD[output,
            "IncludeNyquist" -> True];
        Dimensions[output] === Dimensions[reference] &&
            testNearQ[QuantityMagnitude[identity, {"Picoseconds", 1}],
                QuantityMagnitude[reference, {"Picoseconds", 1}]] &&
            QuantityMagnitude[impulse, {"Picoseconds", 1}][[1, 1]] == 0 &&
            testNearQ[testComplexSpectrum[actualSpectrum][[2 ;; -2]],
                testComplexSpectrum[expectedSpectrum][[2 ;; -2]]] &&
            testNearQ[testComplexSpectrum[
                JerryI`THzTDSTools`TD2FD[impulse,
                    "IncludeNyquist" -> True]][[2 ;; -2]],
                testComplexSpectrum[bareSpectrum][[2 ;; -2]]]
    ],
    True,
    TestID -> "Time reference convolution preserves time origin and Fourier product"
]

VerificationTest[
    Module[{c = 0.0299792458, d = 0.1, n = 1.3, k = 0.02,
            source, metadata, outputs},
        If[!TrueQ[JerryI`THzTDSTools`Private`materialGPUReadyQ[]],
            Return[True]
        ];
        source = N @ Flatten @ Table[
            With[{f = j/100.},
                {f/c, (4 n/(n + 1)^2) Exp[-2 Pi f d k/c],
                    2 Pi f d (n - 1)/c}
            ],
            {j, 0, 511}
        ];
        metadata = N[{(n - 1) d/c, d, 1., 0., 512.}];
        outputs = Table[
            JerryI`THzTDSTools`Private`materialRunGPU[
                source, metadata, 30, 1, 3],
            {3}
        ];
        SameQ @@ outputs
    ],
    True,
    TestID -> "GPU smoothing gives repeatable FP solver results"
]
