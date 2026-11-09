BeginPackage["JerryI`THzTDSTools`", {
    "JerryI`THzTDSTools`Common`"
}]

LinearProjectionsToCircular::usage =
  "LinearProjectionsToCircular[trace1, trace2, {angle1, angle2}] reconstructs the two circular-polarization time traces from linear-analyzer measurements. Equal-length lists are processed pairwise.";

ResampleTraces::usage =
  "ResampleTraces[trace] interpolates a time trace onto a uniform grid.\n ResampleTraces[{trace1, trace2, ...}] uses a common grid for all traces.";
TD2FD::usage =
  "TD2FD[trace] returns a one-sided {frequency, magnitude, phase} spectrum from a time trace. \"IncludeNyquist\" -> True includes the Nyquist bin; lists are processed elementwise.";
FD2TD::usage =
  "FD2TD[spectrum] returns a real time trace from a one-sided {frequency, magnitude, phase} spectrum. \"IncludeNyquist\" -> True expects a Nyquist bin; lists are processed elementwise.";

DeconvolveTraces::usage =
  "DeconvolveTraces[sample, reference] returns a time trace whose TD2FD spectrum is the regularized sample/reference transmission. \"Threshold\" sets the regularization strength.";

DeconvolveFP::usage =
  "DeconvolveFP[response, thickness] removes Fabry-Perot echoes using the iterative slab solver. It returns the same domain as the input: time trace or {frequency, magnitude, phase} spectrum.";

FrequencyDomainConfidenceInterval::usage =
  "FrequencyDomainConfidenceInterval[input] estimates a usable frequency Interval in terahertz from a time trace or spectrum. Lists are processed elementwise.";

FrequencySpectrumTrim::usage =
  "FrequencySpectrumTrim[input, interval] keeps only spectrum bins inside interval. Omitting interval uses FrequencyDomainConfidenceInterval; time traces are Fourier transformed first.";
BandPassFilter::usage =
  "BandPassFilter[input, center, Q] applies a Gaussian band-pass filter (Q defaults to 1). BandPassFilter[input, {{center1, Q1}, ...}] combines multiple bands. Time traces, spectra, and lists are accepted; both forms have operator forms.";
ApodizeTrace::usage =
  "ApodizeTrace[input, position, width] applies a time window centered at position. The default is HannWindow; time traces, spectra, and lists are accepted.";

TimeTraceQ::usage =
  "TimeTraceQ[input] returns True for a time trace, False for a {frequency, magnitude, phase} spectrum, and $Failed for other inputs.";

ToStandardOpticsUnits::usage =
  "ToStandardOpticsUnits[input] returns numeric data in picoseconds or terahertz from a compatible QuantityArray. Lists are processed elementwise.";

AveragedRefrectiveIndex::usage =
  "AveragedRefrectiveIndex[response, thickness] estimates the average refractive index from the response peak delay and a positive sample thickness.";
UnwrapPhase::usage =
  "UnwrapPhase[input, interval] returns a full {frequency, magnitude, phase} spectrum with phase unwrapped inside interval. Omitting interval selects it automatically; linear phase is retained by default.";

EstimateMaterialParameter::usage =
  "EstimateMaterialParameter[parameter, response, thickness] estimates an optical property or FP-reduction score. EstimateMaterialParameter[] lists the supported parameters; Method selects Basic or Iterative extraction.";

GenerateMaterialResponse::usage =
  "GenerateMaterialResponse[n, thickness, {interval, step}] generates a slab transmission spectrum. \nGenerateMaterialResponse[n, thickness, reference] returns the response in the reference domain; \"SolveFP\" and \"ConvolveWithReference\" control reflections and convolution.";

Begin["`Private`"]

LinearProjectionsToCircular::uncovUnits = "Provided QuantityArray units are not convertible to {Time, 1}.";
LinearProjectionsToCircular::angles = "Analyzer angles must be real numeric angles in radians or angular Quantities.";
LinearProjectionsToCircular::dependent = "Analyzer angles `1` and `2` describe linearly dependent projections.";
LinearProjectionsToCircular::wrongfrmt = "Expected two QuantityArray traces or equal-length lists of them, followed by a pair of analyzer angles.";

angleInRadians[angle_] := Quiet @ Check[
    N @ If[MatchQ[angle, _Quantity], QuantityMagnitude[angle, "Radians"], angle],
    $Failed
];

LinearProjectionsToCircular[first_QuantityArray, second_QuantityArray, angles:{_, _}] := Module[
    {radians, determinant, traces, sampled, x, y, circularY, times},
    radians = angleInRadians /@ angles;
    If[!VectorQ[radians, NumericQ] || !VectorQ[radians, TrueQ[Im[#] == 0] &],
        Message[LinearProjectionsToCircular::angles];
        Return[$Failed]
    ];
    determinant = Sin[radians[[2]] - radians[[1]]];
    If[Abs[determinant] < 10^-12,
        Message[LinearProjectionsToCircular::dependent, Sequence @@ angles];
        Return[$Failed]
    ];
    traces = Quiet @ Check[
        QuantityMagnitude[#, {"Picoseconds", 1}] & /@ {first, second},
        $Failed
    ];
    If[traces === $Failed,
        Message[LinearProjectionsToCircular::uncovUnits];
        Return[$Failed]
    ];
    sampled = resample[traces];
    times = sampled[[1, All, 1]];
    (* measured_i = Cos[angle_i] x + Sin[angle_i] y; invert the two projections. *)
    x = (Sin[radians[[2]]] sampled[[1, All, 2]] -
         Sin[radians[[1]]] sampled[[2, All, 2]])/determinant;
    y = (-Cos[radians[[2]]] sampled[[1, All, 2]] +
          Cos[radians[[1]]] sampled[[2, All, 2]])/determinant;
    circularY = iHilbert[y];
    (* Sqrt[2] preserves the original +45/-45-degree amplitude convention. *)
    QuantityArray[{
        Transpose[{times, Sqrt[2] (x + circularY)}],
        Transpose[{times, Sqrt[2] (x - circularY)}]
    }, {"Picoseconds", 1}]
];

LinearProjectionsToCircular[first:{__QuantityArray}, second:{__QuantityArray}, angles:{_, _}] /;
    Length[first] == Length[second] :=
    LinearProjectionsToCircular[#[[1]], #[[2]], angles] & /@ Transpose[{first, second}];

LinearProjectionsToCircular[___] := (Message[LinearProjectionsToCircular::wrongfrmt]; $Failed);

ResampleTraces::wrongfrmt = "Unknown input pattern"
ResampleTraces::uncovUnits = "Provided QuantityArray units are not convertable to {Time, 1}"
ResampleTraces[trace_] := (Message[ResampleTraces::wrongfrmt]; $Failed);

ResampleTraces[trace_QuantityArray] := With[{t = Check[QuantityMagnitude[trace, {"Picoseconds", 1}], Message[ResampleTraces::uncovUnits]; $Failed]},
    If[FailureQ[t],
        $Failed
    ,
        QuantityArray[resample[t], {"Picoseconds", 1}]
    ]
]

ResampleTraces[traces:{__QuantityArray}] := With[{ts = Check[QuantityMagnitude[#, {"Picoseconds", 1}]&/@traces, Message[ResampleTraces::uncovUnits]; $Failed]},
    If[FailureQ[ts],
        $Failed
    ,
        QuantityArray[#, {"Picoseconds", 1}] &/@ resample[ts]
    ]
]

TD2FD::uncovUnits = "Provided QuantityArray units are not convertable to {Time, 1}"
complexSpectrumToEuler[spectrum_List] :=
    ({#[[1]], Abs[#[[2]]], Arg[#[[2]]]} & /@ spectrum);

eulerSpectrumToComplex[spectrum_List] :=
    ({#[[1]], #[[2]] Exp[I #[[3]]]} & /@ spectrum);

Options[TD2FD] = {"IncludeNyquist" -> False};
TD2FD[traces:{__QuantityArray}, opts:OptionsPattern[]] := TD2FD[#, opts] & /@ traces
TD2FD[trace_QuantityArray, opts:OptionsPattern[]] := With[{t = Check[QuantityMagnitude[trace, {"Picoseconds", 1}], Message[TD2FD::uncovUnits]; $Failed]},
    If[FailureQ[t],
        $Failed
    ,
        QuantityArray[
            complexSpectrumToEuler[
                If[TrueQ[OptionValue["IncludeNyquist"]], fourier2dWithNyquist[t], fourier2d[t]]
            ],
            {"Terahertz", 1, 1}
        ]
    ]
]

FD2TD::uncovUnits = "Expected a three-column {frequency, magnitude, phase} QuantityArray with units {Terahertz, 1, 1}."
Options[FD2TD] = {"IncludeNyquist" -> False};
FD2TD[traces:{__QuantityArray}, opts:OptionsPattern[]] := FD2TD[#, opts] & /@ traces
FD2TD[trace_QuantityArray, opts:OptionsPattern[]] := If[
    TimeTraceQ[trace] =!= False,
    Message[FD2TD::uncovUnits];
    $Failed,
    With[{spectrum = eulerSpectrumToComplex[
        QuantityMagnitude[trace, {"Terahertz", 1, 1}]]},
        QuantityArray[
            If[TrueQ[OptionValue["IncludeNyquist"]],
                inverseFourier2dWithNyquist[spectrum], inverseFourier2d[spectrum]
            ],
            {"Picoseconds", 1}
        ]
    ]
];

fourier2dWithNyquist[trace_] := With[
    {sampled = resample[trace]},
    {dt = sampled[[2, 1]] - sampled[[1, 1]], values = sampled[[All, 2]]},
    {padded = If[OddQ[Length[values]], Append[values, 0.], values]},
    {n = Length[padded]},
    Transpose[{Range[0, n/2]/(n dt), Take[Fourier[padded], n/2 + 1]}]
];

inverseFourier2dWithNyquist[spectrum_] := With[
    {amplitudes = spectrum[[All, 2]],
     df = Median[Differences[spectrum[[All, 1]]]]},
    {n = 2 (Length[amplitudes] - 1),
     fullSpectrum = Join[amplitudes,
         Conjugate[Reverse[amplitudes[[2 ;; -2]]]]]},
    Transpose[{Range[0, n - 1]/(n df), Re[InverseFourier[fullSpectrum]]}]
];

DeconvolveTraces::uncovUnits = "Provided QuantityArray units are not convertible to {Time, 1}.";
DeconvolveTraces::threshold = "The Threshold option must be a positive real number.";
DeconvolveTraces::zeroRef = "The reference trace has no signal after its mean is removed.";
DeconvolveTraces[isam_QuantityArray, iref_QuantityArray, OptionsPattern[]] := Module[
    {alpha = OptionValue["Threshold"], data, sam, ref},
    If[!NumericQ[alpha] || !TrueQ[Im[alpha] == 0] || !TrueQ[alpha > 0],
        Message[DeconvolveTraces::threshold];
        Return[$Failed]
    ];
    data = Quiet @ Check[
        resample[QuantityMagnitude[#, {"Picoseconds", 1}] & /@ {isam, iref}],
        $Failed
    ];
    If[data === $Failed,
        Message[DeconvolveTraces::uncovUnits];
        Return[$Failed]
    ];
    {sam, ref} = data;
    If[Max[Abs[ref[[All, 2]] - Mean[ref[[All, 2]]]]] == 0,
        Message[DeconvolveTraces::zeroRef];
        Return[$Failed]
    ];
    QuantityArray[deconv[sam, ref, alpha], {"Picoseconds", 1}]
];

Options[DeconvolveTraces] = {"Threshold"->(10.^-5)}

deconv[sam_, ref_, alpha_:10.^-5] := With[
    {dt = Median[Differences[ref[[All, 1]]]],
     refY = ref[[All, 2]] - Mean[ref[[All, 2]]],
     samY = sam[[All, 2]] - Mean[sam[[All, 2]]]},
    (* Match TD2FD's unitary Fourier convention: TD2FD of this kernel is tHat. *)
    {nDeconv = Length[refY], refHat = Fourier[refY], samHat = Fourier[samY]},
    {power = Abs[refHat]^2},
    {tHat = Conjugate[refHat] samHat/(power + alpha Max[power])},
    Transpose[{
        dt Range[0, Ceiling[nDeconv/2] - 1 + Floor[nDeconv/2]],
        Re[InverseFourier[tHat]]
    }]
]

FrequencyDomainConfidenceInterval::uncovUnits = "Provided units don't match Time or Frequency or not QuantityArray type at all"
FrequencyDomainConfidenceInterval::short = "At least two frequency samples are required to estimate an interval.";

SetAttributes[cached, HoldRest]
SetAttributes[cachedResult, HoldRest]

cachedResult[key_, expression_] := cachedResult[key, _] = expression;

cached[keyExpr_, expr_] := With[{h = Hash[keyExpr, "SHA256"]},
    If[Length[DownValues[cachedResult]] > 256,
        (* reset cached *)
        DownValues[cachedResult] = {};
        cachedResult[key_, expression_] := cachedResult[key, _] = expression;
    ];
    cachedResult[h, expr]
];



FrequencyDomainConfidenceInterval[traces_List] := FrequencyDomainConfidenceInterval/@traces 
FrequencyDomainConfidenceInterval[trace_QuantityArray] := Module[{q},
    q = Switch[TimeTraceQ[trace],
        True, fourier2d[QuantityMagnitude[trace, {"Picoseconds", 1}]],
        False, QuantityMagnitude[trace, {"Terahertz", 1, 1}][[All, {1, 2}]],
        _, $Failed
    ];
    If[q === $Failed,
        Message[FrequencyDomainConfidenceInterval::uncovUnits];
        Return[$Failed]
    ];
    If[Length[q] < 2,
        Message[FrequencyDomainConfidenceInterval::short];
        Return[$Failed]
    ];
    frequencyDomainConfidenceInterval[q]; (* Returns wider interval for some reason. Could be a bug of contexts *)
    frequencyDomainConfidenceInterval[q]
];

FrequencyDomainConfidenceInterval[_] := (Message[FrequencyDomainConfidenceInterval::uncovUnits]; $Failed)

frequencyDomainConfidenceInterval[spectrum_] := Module[
  {d, frequencies, fallback, fit, parameters, candidates, selected, x0, \[Sigma], A},
  d = Map[Function[{pair}, {pair[[1]], Abs[pair[[2]]]^2}], spectrum];
  frequencies = d[[All, 1]];
  fallback = frequencies[[{1, Max[2, Ceiling[Length[frequencies]/2]]}]];

  fit = Quiet[Check[
    NonlinearModelFit[Drop[d, 10], A Exp[-(x0 - x)^2/\[Sigma]^2],
      {\[Sigma], A, x0}, x, ConfidenceLevel -> 0.7, Method -> "NMinimize"],
    $Failed
  ]];
  parameters = If[fit === $Failed, <||>, Quiet[Check[Association[fit["BestFitParameters"]], <||>]]];
  candidates = If[AllTrue[{x0, \[Sigma]}, KeyExistsQ[parameters, #] &],
    {Clip[parameters[x0] - Abs[parameters[\[Sigma]]], {5.0/33.356, Infinity}],
      1.2 Clip[parameters[x0] + Abs[parameters[\[Sigma]]], {30.0/33.356, Infinity}]},
    {}
  ];
  selected = If[Length[candidates] == 2 &&
      AllTrue[candidates, NumberQ[#] && TrueQ[Im[#] == 0] &] &&
      10/33.356 <= candidates[[1]] < candidates[[2]] <= 300/33.356 &&
      candidates[[2]] - candidates[[1]] >= 15/33.356,
    Select[frequencies, candidates[[1]] <= # <= candidates[[2]] &],
    {}
  ];

  Quantity[Interval[If[Length[selected] >= 2, candidates, fallback]],
    "Terahertz"]
]

FrequencySpectrumTrim::units = "Expected a two-column time trace or a three-column {frequency, magnitude, phase} QuantityArray.";
FrequencySpectrumTrim::interval = "Expected a single Interval with Quantity endpoints in terahertz or inverse centimeters.";
FrequencySpectrumTrim::empty = "The specified interval contains no frequency bins.";

frequencySpectrumTrimInput[input_QuantityArray] := Switch[TimeTraceQ[input],
    True, TD2FD[input],
    False, input,
    _, Message[FrequencySpectrumTrim::units]; $Failed
];

frequencyBoundsInTHz[interval_] := Module[{result},
    result = Quiet @ QuantityMagnitude[
        UnitConvert[interval, "Terahertz"], "Terahertz"
    ];
    If[!MatchQ[result, Interval[{_?NumericQ, _?NumericQ}]],
        result = Quiet @ QuantityMagnitude[
            UnitConvert[interval Quantity["SpeedOfLight"], "Terahertz"],
            "Terahertz"
        ]
    ];
    If[!MatchQ[result, Interval[{_?NumericQ, _?NumericQ}]],
        Return[$Failed]
    ];
    First[result]
];

frequencySpectrumTrimBounds[interval_] := With[{bounds = frequencyBoundsInTHz[interval]},
    If[bounds === $Failed, Message[FrequencySpectrumTrim::interval]];
    bounds
];

FrequencySpectrumTrim[input:{__QuantityArray}, rest___] := FrequencySpectrumTrim[#, rest] &/@ input

FrequencySpectrumTrim[input_QuantityArray] := Module[{spectrum, interval},
    spectrum = frequencySpectrumTrimInput[input];
    If[spectrum === $Failed, Return[$Failed]];
    interval = FrequencyDomainConfidenceInterval[spectrum];
    If[interval === $Failed, Return[$Failed]];
    FrequencySpectrumTrim[spectrum, interval]
];

FrequencySpectrumTrim[input_QuantityArray, interval_] := Module[
    {spectrum, bounds, rows},
    spectrum = frequencySpectrumTrimInput[input];
    If[spectrum === $Failed, Return[$Failed]];
    bounds = frequencySpectrumTrimBounds[interval];
    If[bounds === $Failed, Return[$Failed]];
    rows = Select[
        QuantityMagnitude[spectrum, {"Terahertz", 1, 1}],
        bounds[[1]] <= #[[1]] <= bounds[[2]] &
    ];
    If[rows === {},
        Message[FrequencySpectrumTrim::empty];
        Return[$Failed]
    ];
    QuantityArray[rows, {"Terahertz", 1, 1}]
];

FrequencySpectrumTrim[___] := (Message[FrequencySpectrumTrim::units]; $Failed);

BandPassFilter::input = "Expected a real two-column time trace or a three-column {frequency, magnitude, phase} QuantityArray.";
BandPassFilter::center = "The center must be a positive frequency or wavenumber Quantity.";
BandPassFilter::quality = "Q must be a positive dimensionless number or Quantity.";

bandPassPositiveRealQ[value_] :=
    NumericQ[value] && TrueQ[Im[N[value]] == 0] && TrueQ[value > 0];

bandPassCenterTHz[center_Quantity] := With[
    {frequency = Quiet @ Check[QuantityMagnitude[center, "Terahertz"], $Failed]},
    If[bandPassPositiveRealQ[frequency], frequency,
        Quiet @ Check[
            QuantityMagnitude[center Quantity["SpeedOfLight"], "Terahertz"],
            $Failed
        ]
    ]
];

bandPassQualityValue[quality_] := If[Head[quality] === Quantity,
    Quiet @ Check[QuantityMagnitude[quality, "DimensionlessUnit"], $Failed],
    quality
];

bandPassFilterSpectrum[spectrum_QuantityArray, bands_] := With[
    {data = QuantityMagnitude[spectrum, {"Terahertz", 1, 1}]},
    {frequencies = data[[All, 1]]},
    {weights = MapThread[Max,
        Exp[-4 Log[2] (#[[2]] (frequencies - #[[1]])/#[[1]])^2] & /@ bands
    ]},
    QuantityArray[
        Transpose[{frequencies, data[[All, 2]] weights, data[[All, 3]]}],
        {"Terahertz", 1, 1}
    ]
];

BandPassFilter[inputs:{__QuantityArray}, center_Quantity, quality_:1] :=
    BandPassFilter[#, center, quality] & /@ inputs;
BandPassFilter[inputs:{__QuantityArray}, bands:{{_Quantity, _}..}] :=
    BandPassFilter[#, bands] & /@ inputs;

BandPassFilter[center_Quantity, quality_:1] :=
    BandPassFilter[#, center, quality] &;
BandPassFilter[bands:{{_Quantity, _}..}] :=
    BandPassFilter[#, bands] &;

BandPassFilter[input_QuantityArray, center_Quantity, quality_:1] :=
    BandPassFilter[input, {{center, quality}}];

BandPassFilter[input_QuantityArray, bands:{{_Quantity, _}..}] := With[
    {converted = {bandPassCenterTHz[#[[1]]],
                  bandPassQualityValue[#[[2]]]} & /@ bands},
    Which[
        !AllTrue[converted[[All, 1]], bandPassPositiveRealQ],
            Message[BandPassFilter::center]; $Failed,
        !AllTrue[converted[[All, 2]], bandPassPositiveRealQ],
            Message[BandPassFilter::quality]; $Failed,
        True,
            Switch[TimeTraceQ[input],
                False, bandPassFilterSpectrum[input, converted],
                True, With[{original = QuantityMagnitude[input,
                                  {"Picoseconds", 1}]},
                    If[Length[original] < 2 ||
                       !AllTrue[original[[All, 2]], TrueQ[Im[#] == 0] &] ||
                       !AllTrue[Differences[original[[All, 1]]],
                           TrueQ[# > 0] &],
                        Message[BandPassFilter::input]; $Failed,
                        With[{filtered = bandPassFilterSpectrum[
                                  TD2FD[input, "IncludeNyquist" -> True],
                                  converted]},
                            {trace = FD2TD[filtered,
                                 "IncludeNyquist" -> True]},
                            {data = QuantityMagnitude[trace,
                                 {"Picoseconds", 1}]},
                            QuantityArray[
                                Transpose[{data[[All, 1]] + original[[1, 1]],
                                           data[[All, 2]]}],
                                {"Picoseconds", 1}
                            ]
                        ]
                    ]
                ],
                _, Message[BandPassFilter::input]; $Failed
            ]
    ]
];

BandPassFilter[___] := (Message[BandPassFilter::input]; $Failed);

ApodizeTrace::input = "Expected a two-column time trace or a three-column {frequency, magnitude, phase} QuantityArray.";
ApodizeTrace::time = "Position must be a time Quantity and width must be a positive time Quantity.";
ApodizeTrace::window = "WindowFunction must give real numeric weights on [-1/2, 1/2].";
ApodizeTrace::nyquist = "The IncludeNyquist option must be True or False.";
ApodizeTrace::spectrum = "A frequency spectrum must include DC and have strictly increasing frequency bins.";

Options[ApodizeTrace] = {
    "WindowFunction" -> HannWindow, "IncludeNyquist" -> False
};

apodizeRealQ[value_] := NumericQ[value] && TrueQ[Im[N[value]] == 0];
apodizePositiveRealQ[value_] := apodizeRealQ[value] && TrueQ[value > 0];

apodizeTimeTrace[trace_QuantityArray, positionPs_, widthPs_, window_] := Module[
    {data, scaled, weights},
    data = QuantityMagnitude[trace, {"Picoseconds", 1}];
    scaled = (data[[All, 1]] - positionPs)/widthPs;
    weights = Map[
        If[Abs[#] <= 1/2,
            Quiet @ Check[window[#], $Failed],
            0.
        ] &,
        scaled
    ];
    If[!VectorQ[weights, NumericQ] ||
       !AllTrue[weights, TrueQ[Im[N[#]] == 0] &],
        Message[ApodizeTrace::window];
        Return[$Failed]
    ];
    data[[All, 2]] = data[[All, 2]] weights;
    QuantityArray[data, {"Picoseconds", 1}]
];

ApodizeTrace[inputs:{__QuantityArray}, position_, width_,
             opts: OptionsPattern[]] :=
    ApodizeTrace[#, position, width, opts] & /@ inputs;

ApodizeTrace[input_QuantityArray, position_, width_,
             opts: OptionsPattern[]] := Module[
    {positionPs, widthPs, window = OptionValue["WindowFunction"],
     includeNyquist = OptionValue["IncludeNyquist"], domain,
     data, trace, windowed},
    positionPs = If[Head[position] === Quantity,
        Quiet @ Check[QuantityMagnitude[position, "Picoseconds"], $Failed],
        $Failed
    ];
    widthPs = If[Head[width] === Quantity,
        Quiet @ Check[QuantityMagnitude[width, "Picoseconds"], $Failed],
        $Failed
    ];
    If[!apodizeRealQ[positionPs] || !apodizePositiveRealQ[widthPs],
        Message[ApodizeTrace::time];
        Return[$Failed]
    ];
    If[!MemberQ[{True, False}, includeNyquist],
        Message[ApodizeTrace::nyquist];
        Return[$Failed]
    ];
    If[!apodizeRealQ[Quiet @ Check[window[0.], $Failed]],
        Message[ApodizeTrace::window];
        Return[$Failed]
    ];
    domain = TimeTraceQ[input];
    If[domain === $Failed,
        Message[ApodizeTrace::input];
        Return[$Failed]
    ];
    If[domain,
        Return[apodizeTimeTrace[input, positionPs, widthPs, window]]
    ];
    data = QuantityMagnitude[input, {"Terahertz", 1, 1}];
    If[Length[data] < 2 || Abs[data[[1, 1]]] > 10^-9 ||
       !AllTrue[Differences[data[[All, 1]]], TrueQ[# > 0] &],
        Message[ApodizeTrace::spectrum];
        Return[$Failed]
    ];
    trace = FD2TD[input, "IncludeNyquist" -> includeNyquist];
    If[trace === $Failed, Return[$Failed]];
    windowed = apodizeTimeTrace[trace, positionPs, widthPs, window];
    If[windowed === $Failed, Return[$Failed]];
    TD2FD[windowed, "IncludeNyquist" -> includeNyquist]
];

ApodizeTrace[___] := (Message[ApodizeTrace::input]; $Failed);

numericTraceInUnitsQ[input_, units_] := MatrixQ[
    Quiet @ Check[QuantityMagnitude[input, units], $Failed], NumericQ
];

numericEulerSpectrumQ[input_QuantityArray] := With[
    {data = Quiet @ Check[QuantityMagnitude[input, {"Terahertz", 1, 1}], $Failed]},
    MatrixQ[data, NumericQ] &&
        AllTrue[data[[All, 1]], TrueQ[Im[#] == 0] && TrueQ[# >= 0] &] &&
        AllTrue[data[[All, 2]], TrueQ[Im[#] == 0] && TrueQ[# >= 0] &] &&
        AllTrue[data[[All, 3]], TrueQ[Im[#] == 0] &]
];

TimeTraceQ[input_QuantityArray] := Which[
    MatchQ[Dimensions[input], {_?Positive, 2}] &&
        numericTraceInUnitsQ[input, {"Picoseconds", 1}], True,
    MatchQ[Dimensions[input], {_?Positive, 3}] &&
        numericEulerSpectrumQ[input], False,
    True, $Failed
];

TimeTraceQ[_] := $Failed;

ToStandardOpticsUnits[input_QuantityArray] := Switch[TimeTraceQ[input],
    True,
    QuantityMagnitude[input, {"Picoseconds", 1}],
    
    False,
    QuantityMagnitude[input, {"Terahertz", 1, 1}],
    
    _,

    With[{
        v1 = Quiet @ Check[QuantityMagnitude[input, {"Terahertz", 1}], $Failed],
        v2 = Quiet @ Check[QuantityMagnitude[input, {"Terahertz", 1/"Centimeters"}], $Failed]
    },
        If[FailureQ[v1],
            If[FailureQ[v2],
                input
            ,
                v2
            ]
        ,
            v1
        ]
    ]
]

ToStandardOpticsUnits[list_List] := ToStandardOpticsUnits /@ list
ToStandardOpticsUnits[input_] := input

UnwrapPhase::err = "Expected a numeric time trace or frequency spectrum and Method -> None, \"Basic\", \"KK\", \"KKGuided\", or \"KKGuidedMultiscale\".";
UnwrapPhase::interval = "Expected a frequency Interval with Quantity endpoints and at least two spectrum samples inside it.";
UnwrapPhase::linear = "The RemoveLinearPhase option must be True or False.";
Options[UnwrapPhase] = {Method -> "Basic", "RemoveLinearPhase" -> False};

UnwrapPhase[input_QuantityArray, opts: OptionsPattern[]] := With[
    {interval = FrequencyDomainConfidenceInterval[input]},
    If[interval === $Failed, Return[$Failed]];
    UnwrapPhase[input, interval, opts]
];

UnwrapPhase[input_QuantityArray, interval_Quantity, opts: OptionsPattern[]] := Module[
    {domain, range, data, eulerSpectrum, spectrum, indices, selected,
     timeTrace, method, removeLinearPhase, result},
    domain = TimeTraceQ[input];
    If[domain === $Failed,
        Message[UnwrapPhase::err];
        Return[$Failed]
    ];
    range = frequencyBoundsInTHz[interval];
    If[range === $Failed,
        Message[UnwrapPhase::interval];
        Return[$Failed]
    ];
    data = ToStandardOpticsUnits[input];
    eulerSpectrum = If[domain, complexSpectrumToEuler[fourier2d[data]], data];
    indices = Select[Range[Length[eulerSpectrum]],
        range[[1]] <= eulerSpectrum[[#, 1]] <= range[[2]] &];
    selected = eulerSpectrum[[indices]];
    If[Length[selected] < 2,
        Message[UnwrapPhase::interval];
        Return[$Failed]
    ];
    spectrum = eulerSpectrumToComplex[eulerSpectrum];
    method = OptionValue[Method];
    removeLinearPhase = OptionValue["RemoveLinearPhase"];
    If[!MemberQ[{True, False}, removeLinearPhase],
        Message[UnwrapPhase::linear];
        Return[$Failed]
    ];
    timeTrace = If[MemberQ[{"Basic", "KKGuided", "KKGuidedMultiscale"}, method] ||
                   (method === "KK" && !removeLinearPhase),
        If[domain, data, inverseFourier2d[spectrum]],
        Null
    ];
    result = Switch[method,
        "Basic", phaseFromFT[spectrum, timeTrace, range, !removeLinearPhase],
        "KK", If[removeLinearPhase, phaseFromKK[spectrum, range],
            phaseFromKK[spectrum, timeTrace, range, True]],
        "KKGuided", phaseFromFTKKGuided[spectrum, timeTrace, range, !removeLinearPhase],
        "KKGuidedMultiscale", phaseFromFTKKRobust[spectrum, timeTrace, range, !removeLinearPhase],
        None, selected[[All, {1, 3}]],
        _, Message[UnwrapPhase::err]; Return[$Failed]
    ];
    If[!MatchQ[result, {{_?NumericQ, _?NumericQ} ..}] ||
       Length[result] != Length[selected] ||
       Max[Abs[result[[All, 1]] - selected[[All, 1]]]] > 10^-9,
        Message[UnwrapPhase::err];
        Return[$Failed]
    ];
    eulerSpectrum[[indices, 3]] = result[[All, 2]];
    QuantityArray[eulerSpectrum, {"Terahertz", 1, 1}]
];

UnwrapPhase[___] := (Message[UnwrapPhase::err]; $Failed);

AveragedRefrectiveIndex::input = "Expected a two-column time trace or a three-column {frequency, magnitude, phase} spectrum and a positive length Quantity for thickness.";
AveragedRefrectiveIndex::zeroSignal = "The response has no peak from which to estimate a delay.";

AveragedRefrectiveIndex[input_QuantityArray, thickness_Quantity] := Module[
    {thicknessMeters, domain, trace, data, peak, lag},
    thicknessMeters = Quiet @ Check[QuantityMagnitude[thickness, "Meters"], $Failed];
    domain = TimeTraceQ[input];
    If[!NumericQ[thicknessMeters] || !TrueQ[thicknessMeters > 0] || domain === $Failed,
        Message[AveragedRefrectiveIndex::input];
        Return[$Failed]
    ];
    trace = If[domain, input, FD2TD[input]];
    If[trace === $Failed, Return[$Failed]];
    data = QuantityMagnitude[trace, {"Picoseconds", 1}];
    If[Length[data] < 2 || Max[Abs[data[[All, 2]]]] == 0,
        Message[AveragedRefrectiveIndex::zeroSignal];
        Return[$Failed]
    ];
    peak = First[Ordering[Abs[data[[All, 2]]], -1]];
    lag = data[[peak, 1]] - data[[1, 1]];
    1 + QuantityMagnitude[Quantity[lag, "Picoseconds"] Quantity["SpeedOfLight"]/thickness]
];

AveragedRefrectiveIndex[___] := (Message[AveragedRefrectiveIndex::input]; $Failed);

EstimateMaterialParameter::unkparameter = "Unknown optical parameter. Evaluate EstimateMaterialParameter[] to list the supported parameters.";
EstimateMaterialParameter::unkmethod = "Supported methods are \"Basic\" and \"Iterative\".";
EstimateMaterialParameter::iteropts = "Iterative options require positive integer NKIterations and FPIterations, Boolean SolveFP, nonnegative MovingAverage or a Boolean, and Target -> \"Automatic\", \"GPU\", or \"CPU\".";
EstimateMaterialParameter::malformed = "Expected a valid time trace or frequency spectrum and a positive thickness.";
EstimateMaterialParameter::options = "Gain must be positive and PhaseCorrection must be real.";
EstimateMaterialParameter::zero = "The transfer magnitude must be positive at every nonzero frequency used for this parameter.";
EstimateMaterialParameter::fpinterval = "FPInterval must be Automatic, All, or an Interval with endpoints in terahertz or inverse centimeters.";
EstimateMaterialParameter::fpresolution = "The selected spectrum does not resolve the predicted Fabry-Perot round-trip delay.";

EstimateMaterialParameter[] := {"\[Alpha]", "\[Kappa]", "T", "Transmission", "n", "n0", "Complex n", "FP Reduction"};

Options[EstimateMaterialParameter] = {
    Method -> "Basic", "PhaseCorrection" -> 0, "Gain" -> 1.0,
    "NKIterations" -> 30, "FPIterations" -> 8, "SolveFP" -> False,
    "MovingAverage" -> True, "Target" -> "Automatic",
    "FPInterval" -> Automatic
};

materialRealQ[value_] := NumericQ[value] && TrueQ[Im[N[value]] == 0];
materialPositiveRealQ[value_] := materialRealQ[value] && TrueQ[value > 0];

materialInputSpectrum[input_QuantityArray, unwrapPhase_] := Module[
    {domain = TimeTraceQ[input], spectrum, frequencies, data},
    spectrum = Switch[domain,
        True,
            spectrum = TD2FD[input];
            If[spectrum === $Failed, Return[$Failed]];
            If[TrueQ[unwrapPhase],
                frequencies = QuantityMagnitude[
                    spectrum, {"Terahertz", 1, 1}][[All, 1]];
                If[Length[frequencies] < 2, Return[$Failed]];
                UnwrapPhase[input, Interval[{
                    Quantity[First[frequencies], "Terahertz"],
                    Quantity[Last[frequencies], "Terahertz"]
                }], Method -> "Basic", "RemoveLinearPhase" -> False],
                spectrum
            ],
        False, input,
        _, Return[$Failed]
    ];
    If[spectrum === $Failed, Return[$Failed]];
    data = Quiet @ Check[
        QuantityMagnitude[spectrum, {"Terahertz", 1, 1}], $Failed
    ];
    If[data === $Failed || !MatrixQ[data, NumericQ] ||
       Length[data] < 2 ||
       !AllTrue[Differences[data[[All, 1]]], TrueQ[# > 0] &],
        Return[$Failed]
    ];
    data
];

materialInputData[input_QuantityArray, unwrapPhase_, positiveMagnitude_] := With[
    {data = materialInputSpectrum[input, unwrapPhase]},
    If[data === $Failed,
        Message[EstimateMaterialParameter::malformed];
        Return[$Failed]
    ];
    If[TrueQ[positiveMagnitude] &&
       AnyTrue[Select[data, #[[1]] > 0 &], #[[2]] == 0 &],
        Message[EstimateMaterialParameter::zero];
        Return[$Failed]
    ];
    data
];

materialAverageIndex[input_QuantityArray, thicknessCm_] := With[
    {n0 = AveragedRefrectiveIndex[input, Quantity[thicknessCm, "Centimeters"]]},
    If[!materialPositiveRealQ[n0],
        Message[EstimateMaterialParameter::malformed];
        Return[$Failed]
    ];
    n0
];

materialSpectrumFromData[data_List, units_, valueFunction_] := QuantityArray[
    ({#[[1]], valueFunction[#[[1]], #[[2]], #[[3]]]} & /@ data),
    units
];

materialLightSpeedCmPerPs = 0.0299792458;

estimateMaterialSpectrum["Basic", "T", input_QuantityArray, thicknessCm_,
                         gain_, phaseCorrection_] := With[
    {data = materialInputData[input, False, False]},
    If[data === $Failed, Return[$Failed]];
    materialSpectrumFromData[data, {"Terahertz", 1},
        Function[{f, amplitude, phase}, Clip[(gain amplitude)^2, {0, 1}]]]
];



estimateMaterialSpectrum["Basic", "Transmission", input_QuantityArray,
                         thicknessCm_, gain_, phaseCorrection_] :=
    estimateMaterialSpectrum["Basic", "T", input, thicknessCm, gain, phaseCorrection];

estimateMaterialSpectrum["Basic", "\[Kappa]", input_QuantityArray, thicknessCm_,
                         gain_, phaseCorrection_] := With[
    {data = materialInputData[input, False, True]},
    If[data === $Failed, Return[$Failed]];
    materialSpectrumFromData[data, {"Terahertz", 1},
        Function[{f, amplitude, phase},
            If[f == 0, 0,
                -materialLightSpeedCmPerPs Log[gain amplitude]/
                    (2 Pi f thicknessCm)]
        ]]
];

estimateMaterialSpectrum["Basic", "\[Alpha]", input_QuantityArray, thicknessCm_,
                         gain_, phaseCorrection_] := With[
    {data = materialInputData[input, False, True]},
    If[data === $Failed, Return[$Failed]];
    With[{n0 = materialAverageIndex[input, thicknessCm]},
        If[n0 === $Failed, Return[$Failed]];
        With[{fresnelPower = (4 n0/(n0 + 1)^2)^2},
            materialSpectrumFromData[data, {"Terahertz", 1/"Centimeters"},
                Function[{f, amplitude, phase},
                    If[f == 0, 0,
                        (Log[fresnelPower] - 2 Log[gain amplitude])/thicknessCm]
                ]]
        ]
    ]
];

estimateMaterialSpectrum["Basic", "n", input_QuantityArray, thicknessCm_,
                         gain_, phaseCorrection_] := With[
    {data = materialInputData[input, True, True]},
    If[data === $Failed, Return[$Failed]];
    With[{n0 = materialAverageIndex[input, thicknessCm]},
        If[n0 === $Failed, Return[$Failed]];
        materialSpectrumFromData[data, {"Terahertz", 1},
            Function[{f, amplitude, phase},
                If[f == 0, n0,
                    1 + materialLightSpeedCmPerPs (phase + phaseCorrection)/
                        (2 Pi f thicknessCm)]
            ]]
    ]
];

estimateMaterialSpectrum["Basic", "Complex n", input_QuantityArray, thicknessCm_,
                         gain_, phaseCorrection_] := With[
    {data = materialInputData[input, True, True]},
    If[data === $Failed, Return[$Failed]];
    With[{n0 = materialAverageIndex[input, thicknessCm]},
        If[n0 === $Failed, Return[$Failed]];
        materialSpectrumFromData[data, {"Terahertz", 1},
            Function[{f, amplitude, phase},
                If[f == 0, n0,
                    1 + materialLightSpeedCmPerPs
                        (phase + phaseCorrection - I Log[gain amplitude])/
                        (2 Pi f thicknessCm)]
            ]]
    ]
];

estimateMaterialSpectrum[method_, p_String, input_QuantityArray, thicknessCm_,
                         gain_, phaseCorrection_] := (
    Message[If[MemberQ[EstimateMaterialParameter[], p],
        EstimateMaterialParameter::unkmethod,
        EstimateMaterialParameter::unkparameter
    ]];
    $Failed
);

EstimateMaterialParameter[p_String][rest__] := EstimateMaterialParameter[p, rest];

EstimateMaterialParameter[p_String, input_QuantityArray, thickness_Quantity,
                          opts: OptionsPattern[]] := Module[
    {method = OptionValue[Method], gain = OptionValue["Gain"],
     phaseCorrection = OptionValue["PhaseCorrection"], thicknessCm},
    thicknessCm = Quiet @ Check[
        QuantityMagnitude[thickness, "Centimeters"], $Failed
    ];
    If[!materialPositiveRealQ[thicknessCm] || TimeTraceQ[input] === $Failed,
        Message[EstimateMaterialParameter::malformed];
        Return[$Failed]
    ];
    If[!materialPositiveRealQ[gain] || !materialRealQ[phaseCorrection],
        Message[EstimateMaterialParameter::options];
        Return[$Failed]
    ];
    If[p === "n0", Return[AveragedRefrectiveIndex[input, thickness]]];
    If[p === "FP Reduction" && method === "Basic", method = "Iterative"];
    method = materialNormalizeMethod[method,
        OptionValue["NKIterations"], OptionValue["FPIterations"],
        If[p === "FP Reduction", True, OptionValue["SolveFP"]],
        OptionValue["MovingAverage"], OptionValue["Target"]];
    If[method === $Failed, Return[$Failed]];
    If[p === "FP Reduction" &&
       MatchQ[method, {"Iterative", _, _, _, _}],
        method = Append[method, OptionValue["FPInterval"]]
    ];
    estimateMaterialSpectrum[method, p, input, thicknessCm, gain, phaseCorrection]
];

EstimateMaterialParameter[p_String, __] := (
    Message[If[MemberQ[EstimateMaterialParameter[], p],
        EstimateMaterialParameter::malformed,
        EstimateMaterialParameter::unkparameter
    ]];
    $Failed
);
materialKernelDirectory = DirectoryName[$InputFileName];

materialNormalizeMethod[method_, nk_, fp_, solveFP_, movingAverage_, target_] := Module[
    {name = method, rules = <||>, nkCount, fpCount, solve, smoothing, backend},
    If[ListQ[name] && Length[name] >= 1 &&
       MemberQ[{"Iterative", "Interative"}, First[name]],
        rules = Quiet @ Check[Association[Rest[name]], $Failed];
        If[!AssociationQ[rules],
            Message[EstimateMaterialParameter::iteropts];
            Return[$Failed]
        ];
        name = First[name]
    ];
    If[name === "Interative", name = "Iterative"];
    If[name =!= "Iterative", Return[name]];
    nkCount = Lookup[rules, "NKIterations", Lookup[rules, "Iterations", nk]];
    fpCount = Lookup[rules, "FPIterations", Lookup[rules, "FabryPerotCycles", fp]];
    solve = Lookup[rules, "SolveFP", solveFP];
    smoothing = Lookup[rules, "MovingAverage",
        Lookup[rules, "MovingAverageFilter", movingAverage]];
    backend = Lookup[rules, "Target", target];
    If[!IntegerQ[nkCount] || nkCount < 1 ||
       !IntegerQ[fpCount] || fpCount < 1 ||
       !MemberQ[{True, False}, solve] ||
       !(MemberQ[{True, False}, smoothing] ||
         (IntegerQ[smoothing] && smoothing >= 0)) ||
       !MemberQ[{"Automatic", "GPU", "CPU"}, backend],
        Message[EstimateMaterialParameter::iteropts];
        Return[$Failed]
    ];
    {"Iterative", nkCount, If[solve, fpCount, 0],
        If[smoothing === True, 1, If[smoothing === False, 0, smoothing]],
        backend}
];

materialLoadCPUBackend[library_String] := Quiet @ Check[
    With[{backend = <|
        "Load" -> LibraryFunctionLoad[library, "load", {{Real, _, "Manual"}}, Integer],
        "Set" -> LibraryFunctionLoad[library, "set",
            {Integer, {Real, _, "Constant"}}, Integer],
        "Unload" -> LibraryFunctionLoad[library, "unload", {Integer}, Integer],
        "Get" -> LibraryFunctionLoad[library, "get", {Integer}, {Real, _, Automatic}],
        "Run" -> LibraryFunctionLoad[library, "run",
            ConstantArray[Integer, 9], Integer]
    |>},
        If[AllTrue[Values[backend], Head[#] === LibraryFunction &],
            Append[backend, "Library" -> library], $Failed]
    ],
    $Failed
];

(* Exercise the 256-bin boundary and the reusable input buffer once per library. *)
materialCPUBackendQ[backend_Association] := Module[
    {src = $Failed, dest = $Failed, meta = $Failed, first, second, ok = False},
    Quiet @ Check[
        src = backend["Load"][N @ Flatten @ Table[{j, 1., 0.}, {j, 0, 257}]];
        dest = backend["Load"][ConstantArray[0., 5 258]];
        meta = backend["Load"][{0., 0.1, 1., 0., 258.}];
        If[And @@ (IntegerQ /@ {src, dest, meta}),
            ok = backend["Run"][dest, src, meta, 258, 1, 1, 0, 0, 256] === 0;
            If[ok,
                first = Normal[backend["Get"][dest]];
                ok = Length[first] === 5 258 &&
                    Abs[first[[5 256 + 3]] - 1.] < 10^-8
            ];
            If[ok,
                ok = backend["Set"][src,
                    N @ Flatten @ Table[{j, 0.8, 0.}, {j, 0, 257}]] === 0 &&
                    backend["Run"][dest, src, meta, 258, 1, 1, 0, 0, 256] === 0;
                If[ok,
                    second = Normal[backend["Get"][dest]];
                    ok = Abs[second[[5 256 + 3]] - 0.8] < 10^-8 &&
                        Abs[first[[5 256 + 3]] - 1.] < 10^-8
                ]
            ]
        ],
        ok = False
    ];
    Quiet @ Check[
        Scan[If[IntegerQ[#], backend["Unload"][#]] &,
            {src, dest, meta}],
        Null
    ];
    TrueQ[ok]
];

materialCPUBackend[] := materialCPUBackend[] = Module[
    {candidates, backend, library, sourceFile},
    sourceFile = FileNameJoin[{materialKernelDirectory, "nCPU.c"}];
    candidates = Reverse @ Sort @ FileNames[
        FileNameJoin[{DirectoryName[materialKernelDirectory],
            "LibraryResources", $SystemID <> "-v*",
            "nCPU." <> Internal`DynamicLibraryExtension[]}]
    ];
    Do[
        If[FileExistsQ[sourceFile] && FileDate[path] < FileDate[sourceFile],
            Continue[]
        ];
        With[{candidate = materialLoadCPUBackend[path]},
            If[AssociationQ[candidate] && materialCPUBackendQ[candidate],
                backend = Append[candidate, "Origin" -> "Binary"];
                Break[]
            ]
        ],
        {path, candidates}
    ];
    If[AssociationQ[backend], Return[backend]];
    library = Quiet @ Check[
        Needs["CCompilerDriver`"];
        CCompilerDriver`CreateLibrary[
            {sourceFile},
            "nCPU_" <> StringReplace[CreateUUID[], "-" -> ""],
            "TargetDirectory" -> $TemporaryDirectory
        ],
        $Failed
    ];
    If[!StringQ[library], Return[$Failed]];
    backend = materialLoadCPUBackend[library];
    If[AssociationQ[backend] && materialCPUBackendQ[backend],
        Append[backend, "Origin" -> "Built"], $Failed]
];

materialGPUReadyQ[] := materialGPUReadyQ[] = Quiet @ Check[
    Get[FileNameJoin[{materialKernelDirectory, "nGPU.wl"}]];
    TrueQ[JerryI`TDSTools`nGPU`readyQ[]] &&
        Head[JerryI`TDSTools`nGPU`run] =!= Symbol,
    False
];

(* Portable fallback when the corrected native CPU library cannot be built. *)
materialRunPure[source_List, metadata_List, nk_, smooth_, fp_] := Module[
    {data = Partition[source, 3], f, amplitude, phase, correctedAmplitude,
     correctedPhase, thickness, gain, shift, n0, coefficient, n, k, savedK,
     ratio, denominator, averagedN, averagedK, z, count, bad, clipped, pass},
    f = data[[All, 1]];
    amplitude = data[[All, 2]];
    phase = data[[All, 3]];
    correctedAmplitude = amplitude;
    correctedPhase = phase;
    {thickness, gain, shift} = metadata[[2 ;; 4]];
    n0 = 1 + 0.029979 metadata[[1]]/thickness;
    count = Length[data];
    coefficient = 0.159152/(thickness ReplacePart[f, {1 -> 1, 2 -> 1}]);
    coefficient[[1 ;; Min[2, count]]] = 0;
    n = 1 + coefficient (phase + shift);
    k = -coefficient Log[gain ReplacePart[amplitude, {1 -> 1, 2 -> 1}]];
    n[[1 ;; Min[2, count]]] = n0;
    k[[1 ;; Min[2, count]]] = 0;
    Do[
        Do[
            z = (n - I k)/(1 + n - I k)^2;
            n = 1 + coefficient (correctedPhase + shift - Arg[z]);
            k = -coefficient (Log[gain ReplacePart[correctedAmplitude,
                {1 -> 1, 2 -> 1}]] - Log[4 Abs[z]]);
            n[[1 ;; Min[2, count]]] = n0;
            k[[1 ;; Min[2, count]]] = 0,
            {nk}
        ];
        bad = Select[Range[count],
            !NumericQ[n[[#]]] || !NumericQ[k[[#]]] || n[[#]] < 1 &];
        n[[bad]] = n0;
        k[[bad]] = 0;
        If[smooth > 0 && count > smooth + 5,
            averagedN = MovingAverage[n, smooth + 1];
            averagedK = MovingAverage[k, smooth + 1];
            n[[5 ;; count - smooth - 1]] =
                averagedN[[5 ;; count - smooth - 1]];
            k[[5 ;; count - smooth - 1]] =
                averagedK[[5 ;; count - smooth - 1]]
        ];
        If[pass == 0, savedK = k];
        If[pass < fp,
            ratio = (n + I k - 1)/(n + I k + 1);
            denominator = 1 - ratio^2 Exp[I 4 Pi thickness f (n + I k)];
            correctedAmplitude = amplitude Abs[denominator];
            correctedPhase = phase + Arg[denominator];
            clipped = Select[Range[count], correctedAmplitude[[#]] > 1.2 &];
            correctedAmplitude[[clipped]] = 1.;
            correctedPhase[[clipped]] = 0.;
            correctedAmplitude[[1 ;; Min[2, count]]] =
                amplitude[[1 ;; Min[2, count]]];
            correctedPhase[[1 ;; Min[2, count]]] =
                phase[[1 ;; Min[2, count]]]
        ],
        {pass, 0, fp}
    ];
    Flatten[Transpose[{n, k, correctedAmplitude, correctedPhase, savedK}]]
];

If[!ValueQ[materialCPUWorkspace], materialCPUWorkspace = <||>];

materialReleaseCPUWorkspace[] := With[{workspace = materialCPUWorkspace},
    If[KeyExistsQ[workspace, "Size"],
        Quiet @ Check[
            workspace["Backend"]["Unload"] /@
                Lookup[workspace, {"Source", "Dest", "Meta"}],
            Null
        ]
    ];
    materialCPUWorkspace = <||>
];

materialGetCPUWorkspace[backend_, size_] := Module[
    {workspace = materialCPUWorkspace, src, dest, meta},
    If[Lookup[workspace, "Size", 0] === size &&
       Lookup[workspace, "Library", None] === backend["Library"],
        Return[workspace]
    ];
    materialReleaseCPUWorkspace[];
    {src, dest, meta} = Quiet @ Check[{
        backend["Load"][ConstantArray[0., 3 size]],
        backend["Load"][ConstantArray[0., 5 size]],
        backend["Load"][ConstantArray[0., 5]]
    }, {$Failed, $Failed, $Failed}];
    If[!And @@ (IntegerQ /@ {src, dest, meta}),
        Quiet @ Check[
            Scan[If[IntegerQ[#], backend["Unload"][#]] &,
                {src, dest, meta}], Null];
        Return[$Failed]
    ];
    materialCPUWorkspace = <|
        "Library" -> backend["Library"], "Backend" -> backend,
        "Size" -> size, "Source" -> src, "Dest" -> dest, "Meta" -> meta
    |>
];

materialRunCPU[source_List, metadata_List, nk_, smooth_, fp_] := Module[
    {backend = materialCPUBackend[], workspace, size = Length[source]/3, result},
    If[backend === $Failed,
        Return[materialRunPure[source, metadata, nk, smooth, fp]]
    ];
    workspace = materialGetCPUWorkspace[backend, size];
    If[workspace === $Failed,
        Return[materialRunPure[source, metadata, nk, smooth, fp]]
    ];
    result = Quiet @ Check[
        If[backend["Set"][workspace["Source"], source] === 0 &&
           backend["Set"][workspace["Meta"], metadata] === 0 &&
           backend["Run"][workspace["Dest"], workspace["Source"],
               workspace["Meta"], size, 1, nk, smooth, fp, 256] === 0,
            Normal[backend["Get"][workspace["Dest"]]], $Failed
        ],
        $Failed
    ];
    If[ListQ[result] && Length[result] === 5 size,
        result,
        materialReleaseCPUWorkspace[];
        materialRunPure[source, metadata, nk, smooth, fp]
    ]
];

If[!ValueQ[materialGPUWorkspace], materialGPUWorkspace = <||>];

materialReleaseGPUWorkspace[] := With[{workspace = materialGPUWorkspace},
    If[KeyExistsQ[workspace, "Size"],
        Quiet @ Check[
            JerryI`TDSTools`nGPU`unload /@
                Lookup[workspace, {"Source", "Dest"}],
            Null
        ]
    ];
    materialGPUWorkspace = <||>
];

materialGetGPUWorkspace[source_List] := Module[
    {workspace = materialGPUWorkspace, size = Length[source]/3, src, dest},
    If[Lookup[workspace, "Size", 0] =!= size,
        materialReleaseGPUWorkspace[];
        dest = Quiet @ Check[
            JerryI`TDSTools`nGPU`allocate["Float", 5 size],
            $Failed
        ];
        If[dest === $Failed, Return[$Failed]];
        src = Quiet @ Check[
            JerryI`TDSTools`nGPU`load[source, "Float"],
            $Failed
        ];
        If[src === $Failed,
            JerryI`TDSTools`nGPU`unload[dest];
            Return[$Failed]
        ];
        materialGPUWorkspace = <|
            "Size" -> size, "Source" -> src, "SourceData" -> source,
            "Dest" -> dest
        |>;
        Return[materialGPUWorkspace]
    ];
    If[workspace["SourceData"] =!= source,
        src = Quiet @ Check[
            JerryI`TDSTools`nGPU`load[source, "Float"],
            $Failed
        ];
        If[src === $Failed, Return[$Failed]];
        JerryI`TDSTools`nGPU`unload[workspace["Source"]];
        workspace["Source"] = src;
        workspace["SourceData"] = source;
        materialGPUWorkspace = workspace
    ];
    workspace
];

materialRunGPU[source_List, metadata_List, nk_, smooth_, fp_] := Module[
    {workspace, size = Length[source]/3, result},
    If[!TrueQ[materialGPUReadyQ[]], Return[$Failed]];
    workspace = materialGetGPUWorkspace[source];
    If[workspace === $Failed, Return[$Failed]];
    result = Quiet @ Check[
        JerryI`TDSTools`nGPU`run[
            workspace["Dest"], workspace["Source"],
            size, 1, nk, smooth, fp, Sequence @@ metadata[[1 ;; 4]], 256];
        Normal[JerryI`TDSTools`nGPU`get[workspace["Dest"]]],
        $Failed
    ];
    If[ListQ[result] && Length[result] === 5 size,
        result,
        materialReleaseGPUWorkspace[];
        materialGPUReadyQ[] = False;
        $Failed
    ]
];

materialNativeSolution::dc = "The iterative slab solver requires a spectrum containing the zero-frequency bin so that the pulse delay can be estimated.";
materialNativeSolution::backend = "The GPU and native CPU backends could not solve this spectrum.";

materialNativeSolution[input_QuantityArray, thicknessCm_, gain_,
                       phaseCorrection_, nk_, fp_, smooth_, target_] := Module[
    {data = materialInputData[input, True, True], n0, source, metadata, raw, rows, low},
    If[data === $Failed, Return[$Failed]];
    If[Abs[data[[1, 1]]] > 10^-9,
        Message[materialNativeSolution::dc];
        Return[$Failed]
    ];
    n0 = materialAverageIndex[input, thicknessCm];
    If[n0 === $Failed, Return[$Failed]];
    source = N @ Flatten @ Transpose[{
        data[[All, 1]]/materialLightSpeedCmPerPs,
        data[[All, 2]], data[[All, 3]]
    }];
    metadata = N[{(n0 - 1) thicknessCm/materialLightSpeedCmPerPs,
        thicknessCm, gain, phaseCorrection, Length[data]}];
    raw = If[target === "CPU", $Failed,
        materialRunGPU[source, metadata, nk, smooth, fp]];
    If[raw === $Failed, raw = materialRunCPU[source, metadata, nk, smooth, fp]];
    If[!VectorQ[raw, NumericQ] || Length[raw] != 5 Length[data],
        Message[materialNativeSolution::backend];
        Return[$Failed]
    ];
    rows = Partition[raw, 5];
    low = Min[2, Length[data]];
    rows[[1 ;; low, 1]] = n0;
    rows[[1 ;; low, 2]] = 0.;
    rows[[1 ;; low, 3]] = data[[1 ;; low, 2]];
    rows[[1 ;; low, 4]] = data[[1 ;; low, 3]];
    {data, rows}
];

materialIterativeValue["n", f_, row_] := row[[1]];
materialIterativeValue["\[Kappa]", f_, row_] := row[[2]];
materialIterativeValue["Complex n", f_, row_] := row[[1]] + I row[[2]];
materialIterativeValue["\[Alpha]", f_, row_] :=
    4 Pi f row[[2]]/materialLightSpeedCmPerPs;
materialIterativeValue["T", f_, row_, gain_] :=
    Clip[(gain row[[3]])^2, {0, 1}];

materialIterativeResult[{"Iterative", nk_, fp_, smooth_, target_}, p_,
                        input_QuantityArray, thicknessCm_, gain_,
                        phaseCorrection_] := With[
    {solution = materialNativeSolution[input, thicknessCm, gain,
        phaseCorrection, nk, fp, smooth, target]},
    If[solution === $Failed, Return[$Failed]];
    With[{data = solution[[1]], rows = solution[[2]]},
        QuantityArray[
            Table[{data[[j, 1]],
                If[p === "T",
                    materialIterativeValue["T", data[[j, 1]], rows[[j]], gain],
                    materialIterativeValue[p, data[[j, 1]], rows[[j]]]]},
                {j, Length[data]}],
            If[p === "\[Alpha]", {"Terahertz", 1/"Centimeters"},
                {"Terahertz", 1}]
        ]
    ]
];

estimateMaterialSpectrum[{"Iterative", _, 0, _, _}, "T",
                         input_QuantityArray, thicknessCm_, gain_, phaseCorrection_] :=
    estimateMaterialSpectrum["Basic", "T", input, thicknessCm, gain, phaseCorrection];

estimateMaterialSpectrum[config:{"Iterative", _, _, _, _},
                         p:("n" | "\[Kappa]" | "\[Alpha]" | "Complex n" | "T"),
                         input_QuantityArray, thicknessCm_, gain_, phaseCorrection_] :=
    materialIterativeResult[config, p, input, thicknessCm, gain, phaseCorrection];

estimateMaterialSpectrum[config:{"Iterative", _, _, _, _}, "Transmission",
                         input_QuantityArray, thicknessCm_, gain_, phaseCorrection_] :=
    estimateMaterialSpectrum[config, "T", input, thicknessCm, gain, phaseCorrection];

If[!ValueQ[materialFPIntervalCache], materialFPIntervalCache = <||>];

materialFPBounds[input_QuantityArray, Automatic] := If[
    KeyExistsQ[materialFPIntervalCache, "Input"] &&
        materialFPIntervalCache["Input"] === input,
    materialFPIntervalCache["Bounds"],
    With[{interval = FrequencyDomainConfidenceInterval[input]},
        {bounds = If[interval === $Failed, $Failed,
            frequencyBoundsInTHz[interval]]},
        materialFPIntervalCache = <|"Input" -> input, "Bounds" -> bounds|>;
        bounds
    ]
];
materialFPBounds[_, All] := All;
materialFPBounds[_, interval_] := frequencyBoundsInTHz[interval];

materialFPReduction[data_List, rows_List, thicknessCm_, bounds_] := Module[
    {indices, frequencies, raw, corrected, size, step, grid, window,
     rawTransform, correctedTransform, delayStep, lastIndex, predictedDelay,
     candidates, peak, floor, before, after},
    indices = If[bounds === All,
        Range[3, Length[data]],
        Select[Range[3, Length[data]],
            bounds[[1]] <= data[[#, 1]] <= bounds[[2]] &]
    ];
    If[Length[indices] < 16,
        Message[EstimateMaterialParameter::fpresolution];
        Return[$Failed]
    ];
    frequencies = data[[indices, 1]];
    raw = frequencies rows[[indices, 5]];
    corrected = frequencies rows[[indices, 2]];
    size = Length[frequencies];
    step = (Last[frequencies] - First[frequencies])/(size - 1);
    If[Max[Abs[Differences[frequencies] - step]] > 10^-6 step,
        grid = Subdivide[First[frequencies], Last[frequencies], size - 1];
        raw = Interpolation[Transpose[{frequencies, raw}],
            InterpolationOrder -> 1] /@ grid;
        corrected = Interpolation[Transpose[{frequencies, corrected}],
            InterpolationOrder -> 1] /@ grid
    ];
    window = Sin[Pi Range[0, size - 1]/(size - 1)]^2;
    rawTransform = Abs[Fourier[window raw]];
    correctedTransform = Abs[Fourier[window corrected]];
    delayStep = 1/(size step);
    lastIndex = Floor[size/2] + 1;
    predictedDelay = 2 thicknessCm rows[[1, 1]]/materialLightSpeedCmPerPs;
    If[predictedDelay <= 6 delayStep ||
       predictedDelay >= (lastIndex - 7) delayStep,
        Message[EstimateMaterialParameter::fpresolution];
        Return[$Failed]
    ];
    candidates = Select[Range[7, lastIndex - 5],
        Abs[(# - 1) delayStep - predictedDelay] <=
            Max[2 delayStep, 0.25 predictedDelay] &];
    If[candidates === {},
        Message[EstimateMaterialParameter::fpresolution];
        Return[$Failed]
    ];
    peak = First[MaximalBy[candidates, rawTransform[[#]] &]];
    floor = Sqrt[3] Max[
        Median[rawTransform[[peak - 5 ;; peak - 3]]],
        Median[rawTransform[[peak + 3 ;; peak + 5]]],
        Median[correctedTransform[[peak - 5 ;; peak - 3]]],
        Median[correctedTransform[[peak + 3 ;; peak + 5]]],
        $MachineEpsilon
    ];
    before = Norm[rawTransform[[peak - 1 ;; peak + 1]]];
    after = Norm[correctedTransform[[peak - 1 ;; peak + 1]]];
    If[before < 3 floor, 1., before/Max[after, floor]]
];

estimateMaterialSpectrum[{"Iterative", nk_, fp_Integer?Positive, smooth_, target_, interval_},
                         "FP Reduction", input_QuantityArray, thicknessCm_,
                         gain_, phaseCorrection_] := With[
    {bounds = materialFPBounds[input, interval]},
    If[bounds === $Failed,
        Message[EstimateMaterialParameter::fpinterval];
        Return[$Failed]
    ];
    With[{solution = materialNativeSolution[input, thicknessCm, gain,
        phaseCorrection, nk, fp, smooth, target]},
        If[solution === $Failed, Return[$Failed]];
        materialFPReduction[solution[[1]], solution[[2]],
            thicknessCm, bounds]
    ]
];

DeconvolveFP::input = "Expected a valid time trace or frequency spectrum and a positive thickness.";
DeconvolveFP::options = "Gain must be positive and PhaseCorrection must be real.";

Options[DeconvolveFP] = {
    "NKIterations" -> 30, "FPIterations" -> 8,
    "MovingAverage" -> True, "Target" -> "Automatic",
    "Gain" -> 1.0, "PhaseCorrection" -> 0
};

DeconvolveFP[input_QuantityArray, thickness_Quantity,
             opts:OptionsPattern[]] := Module[
    {domain = TimeTraceQ[input], thicknessCm, gain = OptionValue["Gain"],
     phaseCorrection = OptionValue["PhaseCorrection"], settings, solution,
     data, rows, corrected, trace, startTime},
    thicknessCm = Quiet @ Check[QuantityMagnitude[thickness, "Centimeters"], $Failed];
    If[domain === $Failed || !materialPositiveRealQ[thicknessCm],
        Message[DeconvolveFP::input];
        Return[$Failed]
    ];
    If[!materialPositiveRealQ[gain] || !materialRealQ[phaseCorrection],
        Message[DeconvolveFP::options];
        Return[$Failed]
    ];
    settings = materialNormalizeMethod["Iterative",
        OptionValue["NKIterations"], OptionValue["FPIterations"],
        True, OptionValue["MovingAverage"], OptionValue["Target"]];
    If[settings === $Failed, Return[$Failed]];
    solution = materialNativeSolution[input, thicknessCm, gain,
        phaseCorrection, Sequence @@ Rest[settings]];
    If[solution === $Failed, Return[$Failed]];
    {data, rows} = solution;
    corrected = QuantityArray[Transpose[{
        data[[All, 1]], gain rows[[All, 3]], rows[[All, 4]] + phaseCorrection
    }], {"Terahertz", 1, 1}];
    If[!domain, Return[corrected]];
    trace = FD2TD[corrected];
    If[trace === $Failed, Return[$Failed]];
    startTime = QuantityMagnitude[input, {"Picoseconds", 1}][[1, 1]];
    If[startTime == 0, Return[trace]];
    data = QuantityMagnitude[trace, {"Picoseconds", 1}];
    data[[All, 1]] = data[[All, 1]] + startTime;
    QuantityArray[data, {"Picoseconds", 1}]
];

DeconvolveFP[___] := (Message[DeconvolveFP::input]; $Failed);

GenerateMaterialResponse::input = "Expected a positive thickness, a nonnegative frequency interval and positive step, or a time trace or frequency spectrum as reference.";
GenerateMaterialResponse::index = "The refractive-index function must return numeric dimensionless values at the spectrum frequencies in terahertz.";
GenerateMaterialResponse::options = "SolveFP and ConvolveWithReference must be True or False.";

Options[GenerateMaterialResponse] = {
    "SolveFP" -> False, "ConvolveWithReference" -> True
};

materialResponseFactors[indexFunction_, thicknessCm_, frequencies_, solveFP_] := With[
    {index = N[indexFunction /@ frequencies],
     phase = 2 Pi frequencies thicknessCm/materialLightSpeedCmPerPs},
    If[!VectorQ[index, NumericQ],
        Message[GenerateMaterialResponse::index];
        $Failed,
        With[{fresnel = 4 index/(index + 1)^2,
              fp = If[solveFP,
                  1/(1 - ((index - 1)/(index + 1))^2 Exp[2 I phase index]),
                  1.]},
            {Abs[fresnel] Exp[-phase Im[index]] Abs[fp],
             phase (Re[index] - 1) + Arg[fresnel] + Arg[fp]}
        ]
    ]
];

If[!ValueQ[materialResponseReferenceCache], materialResponseReferenceCache = <||>];

materialResponseReference[input_QuantityArray] := If[
    KeyExistsQ[materialResponseReferenceCache, "Input"] &&
        materialResponseReferenceCache["Input"] === input,
    materialResponseReferenceCache,
    With[{result = Switch[Dimensions[input],
        {_Integer, 2},
            With[{data = Quiet @ Check[
                QuantityMagnitude[input, {"Picoseconds", 1}], $Failed]},
                If[!MatrixQ[data, NumericQ] || Length[data] < 2, $Failed,
                    <|"Input" -> input, "TimeDomain" -> True,
                      "Spectrum" -> fourier2dWithNyquist[data],
                      "StartTime" -> data[[1, 1]]|>
                ]
            ],
        {_Integer, 3},
            With[{data = Quiet @ Check[
                QuantityMagnitude[input, {"Terahertz", 1, 1}], $Failed]},
                If[!MatrixQ[data, NumericQ] || Length[data] < 2, $Failed,
                    <|"Input" -> input, "TimeDomain" -> False,
                      "Spectrum" -> data|>
                ]
            ],
        _, $Failed
    ]},
        If[result === $Failed, $Failed,
            materialResponseReferenceCache = result
        ]
    ]
];

materialResponseStepTHz[step_Quantity] := With[
    {direct = Quiet @ Check[QuantityMagnitude[step, "Terahertz"], $Failed]},
    If[NumericQ[direct], N[direct],
        Quiet @ Check[N @ QuantityMagnitude[
            step Quantity["SpeedOfLight"], "Terahertz"], $Failed]
    ]
];

GenerateMaterialResponse[indexFunction:(_Function | _Symbol), thickness_Quantity,
                         {interval_Quantity, step_Quantity}, opts:OptionsPattern[]] := With[
    {thicknessCm = Quiet @ Check[QuantityMagnitude[thickness, "Centimeters"], $Failed],
     bounds = N @ frequencyBoundsInTHz[interval],
     stepTHz = materialResponseStepTHz[step],
     solveFP = OptionValue["SolveFP"]},
    If[!materialPositiveRealQ[thicknessCm] ||
       !MatchQ[bounds, {_?NumericQ, _?NumericQ}] ||
       !materialPositiveRealQ[stepTHz] ||
       !TrueQ[-10^-14 <= bounds[[1]] < bounds[[2]]] ||
       !TrueQ[stepTHz <= bounds[[2]] - Max[0., bounds[[1]]]],
        Message[GenerateMaterialResponse::input];
        $Failed,
        If[!MemberQ[{True, False}, solveFP],
            Message[GenerateMaterialResponse::options];
            $Failed,
            With[{frequencies = Range[Max[0., bounds[[1]]],
                bounds[[2]], stepTHz]},
                With[{factors = materialResponseFactors[indexFunction,
                    thicknessCm, frequencies, solveFP]},
                    If[factors === $Failed, $Failed,
                        QuantityArray[Transpose[{frequencies,
                            factors[[1]], factors[[2]]}],
                            {"Terahertz", 1, 1}]
                    ]
                ]
            ]
        ]
    ]
];

GenerateMaterialResponse[indexFunction:(_Function | _Symbol), thickness_Quantity,
                         inputReference_QuantityArray, opts:OptionsPattern[]] := With[
    {thicknessCm = Quiet @ Check[QuantityMagnitude[thickness, "Centimeters"], $Failed],
     solveFP = OptionValue["SolveFP"],
     convolve = OptionValue["ConvolveWithReference"]},
    If[!materialPositiveRealQ[thicknessCm],
        Message[GenerateMaterialResponse::input];
        $Failed,
        If[!MemberQ[{True, False}, solveFP] ||
           !MemberQ[{True, False}, convolve],
            Message[GenerateMaterialResponse::options];
            $Failed,
            With[{reference = materialResponseReference[inputReference]},
                If[reference === $Failed,
                    Message[GenerateMaterialResponse::input];
                    $Failed,
                    With[{spectrum = reference["Spectrum"]},
                        With[{factors = materialResponseFactors[indexFunction,
                            thicknessCm, spectrum[[All, 1]], solveFP]},
                            If[factors === $Failed, $Failed,
                                If[reference["TimeDomain"],
                                    With[{trace = inverseFourier2dWithNyquist[
                                        Transpose[{spectrum[[All, 1]],
                                            If[convolve, spectrum[[All, 2]], 1.]
                                                factors[[1]] Exp[I factors[[2]]]}]
                                    ]},
                                        QuantityArray[Transpose[{
                                            trace[[All, 1]] + If[convolve,
                                                reference["StartTime"], 0.],
                                            trace[[All, 2]]}], {"Picoseconds", 1}]
                                    ],
                                    QuantityArray[Transpose[{spectrum[[All, 1]],
                                        If[convolve, spectrum[[All, 2]], 1.] factors[[1]],
                                        If[convolve, spectrum[[All, 3]], 0.] + factors[[2]]}],
                                        {"Terahertz", 1, 1}]
                                ]
                            ]
                        ]
                    ]
                ]
            ]
        ]
    ]
];

GenerateMaterialResponse[___] := (Message[GenerateMaterialResponse::input]; $Failed);


End[]
EndPackage[]
