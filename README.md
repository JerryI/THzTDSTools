# Time-Domain Terahertz Spectroscopy Tools 🧰

A small library for material parameter extraction from time-domain THz spectra, written in Wolfram Language, C, and OpenCL.

*Developed at University of Augsburg, Germany 🇩🇪*

🔗 Hands-on notebook with several use cases is available at [Example.html](https://jerryi.github.io/THzTDSTools/Example.html)

🔗 Jones matrix tutorial is available at [Example.html](https://jerryi.github.io/THzTDSTools/Jones.html)

*THzTDSTools* is designed mainly for processing __time-domain data in transmission optics__. It includes:

- signal (kernel) deconvolution;
- seamless conversion between the time and frequency domains;
- projections from linear to circular bases;
- extraction of material parameters (absorption, refractive index, etc.), __including multiple-beam propagation and Fresnel losses__;
- advanced 1D phase unwrapping guided by the Kramers–Kronig relation;
- optical-response modeling for thin and thick slabs;
- band-pass filtering and apodization of time traces;
- Jones matrix framework and Stokes light field analysis.

All functions are strongly typed but also support multiple input patterns: time-domain data is automatically converted to the frequency domain if needed. Iterative solvers are implemented in OpenCL, with C functions running on the CPU as a fallback and a WL implementation serving as a secondary fallback if the C version fails to load.

## GPU backends
- [x] OpenCL
- [ ] CUDA

## Supported Platforms

- WLJS Notebook (free)
- Wolfram Mathematica (paid)
- Wolfram Engine (free)
 
## Data format and units

The public processing functions use `QuantityArray` data:

| Domain | Rows | Units |
| --- | --- | --- |
| Time trace | `{time, amplitude}` | `{"Picoseconds", 1}` |
| One-sided frequency spectrum | `{frequency, magnitude, phase}` | `{"Terahertz", 1, 1}` |

Spectrum magnitudes are nonnegative. Phases are in radians and may lie outside `[-Pi, Pi]`, allowing an unwrapped phase to be stored.

Frequency-dependent material results have `{frequency, value}` rows. The absorption coefficient `\[Alpha]` is in inverse centimeters.

The package's Fourier convention gives a delayed pulse a positive linear phase. Thicknesses, time-window positions, and frequency bounds are `Quantity` values.

## Example workflow

Given measured `sample` and `reference` time traces:

```wl
Needs["JerryI`THzTDSTools`"];
response = DeconvolveTraces[sample, reference];
spectrum = UnwrapPhase[response, Method -> "KKGuidedMultiscale"];
thickness = Quantity[0.635, "Millimeters"];

index = EstimateMaterialParameter["n", spectrum, thickness,
  Method -> "Iterative", "SolveFP" -> True];
absorption = EstimateMaterialParameter["\[Alpha]", spectrum, thickness,
  Method -> "Iterative", "SolveFP" -> True];
fpFreeSpectrum = DeconvolveFP[spectrum, thickness];
```

`DeconvolveTraces` returns a time-domain transfer response. Its `TD2FD` spectrum is the regularized complex sample/reference transmission.

Keep the linear phase when unwrapping if you will invert the spectrum to time or use it to estimate the refractive index. `UnwrapPhase` retains it by default.

## Function reference

### Traces and Fourier transforms

**`LinearProjectionsToCircular[trace1, trace2, {angle1, angle2}]`**

Reconstructs two circular-polarization time traces from measurements at two linear-analyzer angles. It treats each measured projection as `Cos[angle] Ex + Sin[angle] Ey`.

The outputs are proportional to `Sqrt[2] (Ex ± Hilbert[Ey])`. The `Sqrt[2]` factor preserves the earlier ±45-degree amplitude convention.
 
Bare angles are in radians; `45 Degree` and angular `Quantity` values also work. The angles cannot differ by a multiple of 180 degrees. Equal-length lists of traces are processed pairwise.

> **Note:** This function only changes the basis of two measured linear projections. For a general polarization response — including cross-polarized components and sample/reference calibration — construct a `JonesMatrix` from the four required polarizer-analyzer configurations; see [Jones analysis](#jones-analysis).

**`ResampleTraces[trace]`**

Interpolates a trace onto a uniform time grid using its smallest step. For a list of traces, it uses the smallest step across all inputs and a common grid spanning them. Values outside an input trace are zero. Output times are in picoseconds.

**`TD2FD[trace]`**

Returns a one-sided spectrum with DC. By default, it omits the Nyquist bin. Set `"IncludeNyquist" -> True` to retain Nyquist; an odd-length resampled trace is then padded with one zero sample. Lists of traces are processed elementwise.

**`FD2TD[spectrum]`**

Reconstructs a real, even-length trace starting at time zero. It accepts unwrapped phase values and processes lists of spectra elementwise.

By default, it interpolates onto a uniform frequency grid and assumes a zero Nyquist component. Set `"IncludeNyquist" -> True` when the spectrum contains both DC and Nyquist bins.

**`DeconvolveTraces[sample, reference]`**

Resamples both time traces to a common grid, removes their means, and returns a dimensionless lag trace. The `"Threshold"` option controls regularization of the frequency-domain transfer function; its positive relative power-spectrum value defaults to `10^-5`.

FFT lags are circular. A negative-lag response appears near the end of the output rather than at a negative time. The intended response is causal, with the sample signal arriving no earlier than the reference.

### Frequency selection and filtering

**`FrequencyDomainConfidenceInterval[input]`**

Estimates a usable band from a Gaussian fit to spectral power. It accepts a time trace, a spectrum, or a list of either, and returns a frequency `Quantity` containing an `Interval` in terahertz.

If the fit is unsuitable, it falls back to a band from the first frequency to one near the midpoint.

**`FrequencySpectrumTrim[input, interval]`**

Keeps only frequency bins within `interval`. Endpoints may be in terahertz or inverse centimeters. If `interval` is omitted, the function uses `FrequencyDomainConfidenceInterval[input]`.

Time input is transformed with `TD2FD` first; the output is always a spectrum. Lists are processed elementwise. An interval containing no bins gives `$Failed`.

**`BandPassFilter[input, center, Q]`**

**`BandPassFilter[input, {{center1, Q1}, {center2, Q2}, ...}]`**

Applies one or more Gaussian frequency bands. Each center is a positive frequency or wavenumber `Quantity`; each `Q` is positive and dimensionless. The single-center form defaults to `Q = 1`.

For a band centered at `center` with quality factor `Q`, the amplitude response is:

```wl
Exp[-4 Log[2] ((f - center)/(center/Q))^2]
```

The full width at half maximum is `center/Q`. Multiple bands use the pointwise maximum of their responses, so overlapping bands never amplify the signal above the input amplitude.

Spectrum input retains its frequencies and phases. Time input returns a filtered trace with its starting time preserved; an odd sample count may gain one padded point. Lists are processed elementwise. Both forms also support operator syntax, such as `BandPassFilter[{{center1, Q1}, {center2, Q2}}][input]`.

**`ApodizeTrace[input, position, width]`**

Applies a time window centered at `position`. Both `position` and the positive `width` are time `Quantity` values. The default window is `HannWindow`.

Set `"WindowFunction" -> window` to supply a function of `(time - position)/width` on `[-1/2, 1/2]`. The output is zero outside that range.

- **Time input:** Keeps the original time grid.
- **Spectrum input:** Must contain DC. The spectrum is converted to time, windowed, and transformed back. Here `position` is relative to the reconstructed time axis starting at zero, and the output phase can be wrapped.

For a spectrum containing Nyquist, set `"IncludeNyquist" -> True` (default `False`).

Lists are processed elementwise.

### Phase unwrapping

`TD2FD` stores the principal Fourier phase in the range from $-\pi$ to $\pi$. A sampled wrapped phase does not reveal how many full $2\pi$ turns occur between frequency bins.

A nearest-neighbor unwrap can fail when the phase changes rapidly or transmission has a deep minimum. The Kramers–Kronig (KK) guided method uses the measured magnitude to help choose the phase branch.

**`UnwrapPhase[input, interval]`**

Returns a complete spectrum. It changes the phase inside the selected band while preserving the bins outside it. Time input is Fourier transformed first.

Omitting `interval` uses `FrequencyDomainConfidenceInterval[input]`. An explicit interval may have terahertz or inverse-centimeter endpoints. The available methods are:

| `Method` | Behavior |
| --- | --- |
| `"Basic"` (default) | Put each adjacent Fourier phase step in $[-\pi,\pi]$, then accumulate. |
| `"KK"` | Estimate phase from magnitude using Kramers–Kronig relations. |
| `"KKGuided"` | Use the KK estimate to choose Fourier-phase branches. |
| `"KKGuidedMultiscale"` | Also correct branch choices around deep transmission nulls. |
| `None` | Keep spectrum phase unchanged, or return wrapped Fourier phase for time input. |

With `Method -> "KKGuided"`, the library:

1. Finds the pulse-arrival time from the peak of the absolute time response and temporarily removes its linear phase, $2\pi f t_{\mathrm{peak}}$.
2. Computes a provisional KK phase from $\log |H(f)|$ using an even spectral extension and a Hilbert transform.
3. For each pair of neighboring frequencies, adds the multiple of $2\pi$ to the measured phase step that makes it closest to the corresponding KK phase step:

   $$\Delta\phi_i = \Delta\phi_i^{\mathrm{FT}} + 2\pi\,\mathrm{Round}\left(\frac{\Delta\phi_i^{\mathrm{KK}}-\Delta\phi_i^{\mathrm{FT}}}{2\pi}\right).$$

   It accumulates these steps starting from the first measured phase. Thus the KK phase guides branch selection; it does not replace the Fourier phase.
4. Restores the pulse-arrival phase unless `"RemoveLinearPhase" -> True` was requested.

`"KKGuidedMultiscale"` adds a correction near deep transmission minima. It compares the guided and KK phases on both sides of each candidate minimum at several window sizes. When those comparisons agree, it applies a $2\pi$ branch shift.

The KK estimate is approximate for finite-band data. Neither guided method guarantees a unique phase where the signal is weak or the frequency sampling is too coarse.

`"RemoveLinearPhase" -> False` is the default and retains the pulse-arrival phase term. Setting it to `True` removes that term inside the interval, so the result no longer represents the original time delay.

This option has no effect with `Method -> None`. The KK-only phase is estimated from magnitude and need not invert to the original trace.

### Material properties

**`AveragedRefrectiveIndex[response, thickness]`**

Uses the peak of the absolute time-domain response to estimate `n0 = 1 + c lag/thickness`. Spectrum input is converted with `FD2TD` first. Thickness must be a positive length `Quantity`; the estimate assumes a causal response.

**`EstimateMaterialParameter[parameter, response, thickness]`**

Accepts a deconvolved time trace or a one-sided spectrum. A spectrum supplied for index estimation should already have unwrapped phase. A time trace is unwrapped internally with the Basic method, retaining its linear phase.

`EstimateMaterialParameter[]` lists the supported parameters:

| Parameter | Result |
| --- | --- |
| `"n0"` | Peak-delay average refractive index. |
| `"n"` | Frequency-dependent refractive index. |
| `"\[Kappa]"` | Frequency-dependent extinction coefficient. |
| `"\[Alpha]"` | Absorption coefficient in inverse centimeters. |
| `"T"` or `"Transmission"` | Corrected power transmission, clipped to `[0, 1]`. |
| `"Complex n"` | `n + I \[Kappa]`. |
| `"FP Reduction"` | Dimensionless before/after FP-suppression score. |

`Method -> "Basic"` (default) estimates the material properties in one pass. Use it when internal reflections can be neglected.

Let $A(f)$ and $\phi(f)$ be the magnitude and unwrapped phase of the measured transmission, $L$ the thickness, $g$ the gain correction, and $\phi_C$ the phase correction. Here $f$ is frequency and $c$ is the speed of light. For $f>0$:

$$n(f)=1+\frac{c[\phi(f)+\phi_C]}{2\pi fL}$$

$$\kappa(f)=-\frac{c\ln[gA(f)]}{2\pi fL}$$

The absorption coefficient includes a Fresnel correction based on the peak-delay estimate $n_0$. Writing $F_0=4n_0/(n_0+1)^2$:

$$\alpha(f)=\frac{2}{L}\ln\left(\frac{F_0}{gA(f)}\right)$$

Power transmission is `T = Clip[(g A)^2, {0, 1}]`. At $f=0$, `n` and `Complex n` use $n_0$; `\[Alpha]` and `\[Kappa]` are set to zero. The defaults are `"Gain" -> 1` and `"PhaseCorrection" -> 0` radians.

`Method -> "Iterative"` runs the slab solver and needs a spectrum containing DC. The historical spelling `"Interative"` is also accepted.

Solver settings can be supplied as top-level options or in a method list such as `Method -> {"Iterative", "Iterations" -> 30}`.

| Option | Default | Meaning |
| --- | --- | --- |
| `"NKIterations"` | `30` | Index/extinction iterations. |
| `"FPIterations"` | `8` | Fabry–Perot iterations when enabled. |
| `"SolveFP"` | `False` | Include internal reflections in the iterative solve. |
| `"MovingAverage"` | `True` | Smoothing setting for the solver; also accepts a nonnegative integer. |
| `"Target"` | `"Automatic"` | Try GPU first, then CPU fallback; `"GPU"` and `"CPU"` are also accepted. |

`"FP Reduction"` always runs the iterative solver with FP removal. It compares the `f \[Kappa]` spectra before and after removal near the predicted round-trip delay.

| Score | Interpretation |
| --- | --- |
| Above 1 | FP-related structure was suppressed. |
| Below 1 | The structure became stronger. |
| Exactly 1 | No measurable reduction, or no detectable FP peak. |

This is a diagnostic, not a thickness estimate. A material feature with the same periodicity can distort the score. If the chosen frequency band cannot resolve the echo delay, the result is `$Failed`.

Set **`"FPInterval"` as a top-level option**, outside `Method`:

- `Automatic` (default) uses `FrequencyDomainConfidenceInterval`.
- `All` uses all positive-frequency bins.
- A frequency `Interval` restricts the score to a chosen measurement band.

A fine thickness scan is useful when the suppression peak is narrow:

```wl
Table[{d, EstimateMaterialParameter["FP Reduction", spectrum,
    Quantity[d, "Millimeters"], "Gain" -> 1.2,
    "FPInterval" -> Interval[{Quantity[1.5, "Terahertz"],
                              Quantity[2.5, "Terahertz"]}]]},
  {d, 0.62, 0.65, 0.0025}]
```

**`DeconvolveFP[response, thickness]`**

Removes Fabry–Perot echoes using iterative estimates of `n` and `\[Kappa]`. Thickness must be a positive length `Quantity`.

- **Time trace:** Fourier transformed and phase-unwrapped internally; returns a time trace.
- **Spectrum:** Must contain DC and is assumed to have already-unwrapped phase; returns a spectrum.

| Option | Default |
| --- | --- |
| `"NKIterations"` | `30` |
| `"FPIterations"` | `8` |
| `"MovingAverage"` | `True` |
| `"Target"` | `"Automatic"` |
| `"Gain"` | `1.0` |
| `"PhaseCorrection"` | `0` |

Gain and phase correction are included in the output. When passing that output to `EstimateMaterialParameter`, use its default gain and phase options.

**`GenerateMaterialResponse[n, thickness, {interval, step}]`**

Generates a slab-transmission spectrum as a `{frequency, magnitude, phase}` `QuantityArray`. The function `n` receives a plain numeric frequency in terahertz and returns the dimensionless complex refractive index `n + I \[Kappa]`. It may be a Wolfram `Function` or a symbol defined as a function.

Thickness must be a positive length `Quantity`.

The interval is a frequency `Quantity` containing an `Interval`, such as `Quantity[Interval[{0, 4}], "Terahertz"]`. The positive step is another frequency `Quantity`. Both accept terahertz or inverse centimeters. The grid includes the upper endpoint when it falls on a step.

The generated phase retains the slab delay and can exceed `Pi`.

`"SolveFP" -> False` (default) uses the single-pass response. Set it to `True` to include internal reflections. With complex index `m = n + I \[Kappa]`, the model is:

```wl
H0 = (4 m/(m + 1)^2) Exp[I 2 Pi f L (m - 1)/c]
FP = 1/(1 - ((m - 1)/(m + 1))^2 Exp[I 4 Pi f L m/c])
H  = H0 FP  (* when "SolveFP" -> True *)
```

Here `f` is in terahertz, `L` in centimeters, and `c` in centimeters per picosecond. The positive phase sign follows the `TD2FD` and iterative-solver convention.

**`GenerateMaterialResponse[n, thickness, reference]`**

Evaluates the same transfer function at the frequencies of `reference`. By default, `"ConvolveWithReference" -> True` multiplies the model by the reference:

- **Spectrum reference:** Returns a spectrum. The existing unwrapped phase is retained and the model phase is added.
- **Time-trace reference:** Fourier transforms the trace, multiplies by the model, and returns a trace with the original starting time.

Set `"ConvolveWithReference" -> False` to use the reference only for its sampling grid and output domain. A spectrum then returns the bare transfer spectrum; a time trace returns the impulse response starting at zero lag. This option is ignored for the `{interval, step}` form.

Time input uses the full one-sided spectrum, including Nyquist. An odd-length trace gains one padded sample. The real inverse transform discards any imaginary DC or Nyquist component of the modeled spectrum.

Convolution is circular over the finite FFT window, so allow enough time after the pulse for delayed FP echoes. Repeated calls with the same `QuantityArray` reuse the most recent reference transform.

```wl
nModel[f_] := 1.8 + 0.02 I;
thickness = Quantity[0.635, "Millimeters"];
response = GenerateMaterialResponse[nModel, thickness,
  {Quantity[Interval[{0, 4}], "Terahertz"], Quantity[0.01, "Terahertz"]},
  "SolveFP" -> True];
simulatedTrace = GenerateMaterialResponse[nModel, thickness, referenceTrace,
  "SolveFP" -> True];
```

### Jones analysis

```wl
Needs["JerryI`THzTDSTools`"];
Needs["JerryI`THzTDSTools`Jones`"];
```

**`JonesMatrix[sample, reference]`**

Constructs a calibrated, frequency-dependent Jones matrix from angle-labeled sample and reference measurements. Each argument must contain the same four polarizer-analyzer configurations:

```wl
sample = {
  {0 Degree, 0 Degree} -> sample00,
  {90 Degree, 90 Degree} -> sample9090,
  {45 Degree, 0 Degree} -> sample450,
  {45 Degree, 90 Degree} -> sample4590
};

reference = {
  {0 Degree, 0 Degree} -> reference00,
  {90 Degree, 90 Degree} -> reference9090,
  {45 Degree, 0 Degree} -> reference450,
  {45 Degree, 90 Degree} -> reference4590
};

jones = JonesMatrix[sample, reference];
```

In an angle pair `{analyzer, inputPolarizer}`, the first angle belongs to the analyzer and the second to the input polarizer. Each value may be a time trace or a one-sided frequency spectrum in the standard `QuantityArray` format. Time traces are resampled and transformed to the frequency domain automatically.

At each sampled frequency, let $M_{a,b}$ be the complex field measured at analyzer angle $a$ and input-polarizer angle $b$. The four measurements reconstruct the system transfer matrix:

$$
T=\begin{pmatrix}
M_{0,0} & \sqrt{2}M_{45,90}-M_{90,90}\\
\sqrt{2}M_{45,0}-M_{0,0} & M_{90,90}
\end{pmatrix}.
$$

The sample Jones matrix is calibrated from the left:

$$
J_{\mathrm{sample}}=T_{\mathrm{reference}}^{-1}T_{\mathrm{sample}}.
$$

The implementation evaluates this with `LinearSolve[Treference, Tsample]` instead of explicitly forming the inverse. Matrix order matters because Jones matrices generally do not commute.

Evaluate the resulting object at a frequency `Quantity` to obtain a numeric 2×2 matrix.

Evaluate it at a frequency-interval `Quantity`, such as `Quantity[Interval[{0.5, 1.5}], "Terahertz"]`, to obtain a weighted average of its sampled matrices. Before averaging, each matrix is rotated to match the common phase of the highest-weight matrix. This preserves the polarization response when propagation delay makes the phase wind across the band. The result retains the anchor matrix's common phase; use a single frequency when its absolute phase matters.

```wl
(* Inspect the calibrated matrix at one frequency. *)
matrix = jones[Quantity[1.25, "Terahertz"]];
matrix // MatrixForm

(* Average over a useful spectral band. *)
band = Quantity[Interval[{1.15, 1.5}], "Terahertz"];
averageMatrix = jones[band];
```

**`OpticalAxisAnalysis[matrix]`**

**`OpticalAxisAnalysis[jones, interval]`**

For a numeric 2×2 transmission matrix, returns an optical-axis analysis association directly. For a frequency-dependent `JonesMatrix`, first evaluates its phase-aligned, weighted average over `interval`.

Omitting `interval` uses the interval stored in `jones`.

The `"LinearizedEigenvectors"` field contains two real unit vectors along the eigenpolarization ellipses' major axes. It uses their azimuths, so these ideal linear states need not be eigenvectors of the measured matrix. The azimuth becomes unstable for nearly circular states.

For example, inspect the eigenpolarizations at one frequency and then obtain a more robust, band-averaged estimate of the linear optical axes:

```wl
analysis = OpticalAxisAnalysis[
  jones[Quantity[1.25, "Terahertz"]]
];
analysis[["Eigenpolarizations", All,
  {"JonesVector", "AzimuthDeg"}]] // TableView

bandAnalysis = OpticalAxisAnalysis[jones, band];
basis = bandAnalysis["LinearizedEigenvectors"];
```

**`OpticalAxesTransmission[jones, opticalAxes]`**

Accepts the association returned by `OpticalAxisAnalysis[jones]` or two numeric Jones vectors. `"Circular"` selects the fixed pair `{1, I}` and `{1, -I}`.

For each vector `v`, it evaluates the co-polarized complex transmission across the full frequency grid:

```wl
Conjugate[v].J[f].v/(Conjugate[v].v)
```

The result is a list of two `{frequency, magnitude, phase}` `QuantityArray` spectra. The projection formula normalizes the vectors, and `Arg` wraps the phases.

Cross-polarized conversion remains in the full Jones matrix.

```wl
axes = OpticalAxisAnalysis[jones];

eigenSpectra = OpticalAxesTransmission[jones, axes];
linearSpectra = OpticalAxesTransmission[
  jones, axes["LinearizedEigenvectors"]];
circularSpectra = OpticalAxesTransmission[jones, "Circular"];

angle = 30 Degree;
customSpectra = OpticalAxesTransmission[jones,
  {{Cos[angle], Sin[angle]}, {-Sin[angle], Cos[angle]}}];
```

The linearized axes are useful for a weakly misaligned linearly birefringent sample. Project the full Jones matrix onto them to recover the two co-polarized eigenmode spectra:

```wl
eigenSpectra = OpticalAxesTransmission[jones, basis];

ListLinePlot[
  (FrequencySpectrumTrim /@ eigenSpectra)[[All, All, {1, 2}]],
  PlotRange -> {{0.25, 1.55}, {0, 1}},
  Frame -> True,
  FrameLabel -> {"frequency (THz)", "transmission"},
  PlotLegends -> {"Polarization 1", "Polarization 2"}
]

timeTraces = FD2TD /@ eigenSpectra;
filteredTraces =
  BandPassFilter[#, Quantity[0.7, "Terahertz"], 2.5] & /@ timeTraces;

ListLinePlot[filteredTraces[[All, ;; 600]],
  PlotRange -> Full,
  Frame -> True,
  FrameLabel -> {"time (ps)", "amplitude"},
  PlotLegends -> {"Polarization 1", "Polarization 2"}
]
```

**`PolarizationEllipse[matrix, scale]`**

Draws two polarization ellipses from a numeric 2×2 transmission matrix. `scale` defaults to 1, and `PlotStyle` controls their colors.

**`PolarizationEllipse[jones, interval, scale]`**

Returns graphics primitives for the two polarization ellipses at each selected frequency. It uses the stored interval by default and samples up to 100 frequencies deterministically.

**`PolarizationEllipse[{jones1, jones2, ...}, interval, scale]`**

Returns one set of graphics primitives per Jones matrix, using the same sampled frequencies for all. By default the common band is the intersection of their stored intervals. Supply `interval` to compare a chosen band instead.

Each frequency group's opacity comes from its spectral weight. `PlotStyle` is passed to the per-frequency ellipse drawing.

```wl
Graphics[
  PolarizationEllipse[jones[Quantity[1.25, "Terahertz"]]],
  PlotRange -> 1.2,
  AspectRatio -> 1,
  Axes -> True
]

(* Overlay frequency-dependent ellipses across the stored confidence band. *)
Graphics[PolarizationEllipse[jones], PlotRange -> 1.2]

(* Compare temperatures on one frequency grid. *)
temperatures = {100., 70., 60., 50., 5.};
Graphics[#, PlotRange -> 1.2, AspectRatio -> 1] & /@
  PolarizationEllipse[d /@ temperatures]
```

### Utilities

**`AveragedRefrectiveIndex[response, thickness]`**

Estimates the mean index from the response pulse delay. For frequency spectra, it smooths the reconstructed time trace over three samples before finding the pulse, so isolated inversion spikes do not set the delay. Time-trace inputs use their original samples.

**`TimeTraceQ[input]`**

Returns `True` for a numeric, two-column time `QuantityArray`; `False` for a numeric, three-column Euler-form spectrum; and `$Failed` otherwise. Units must be compatible with the formats above.

**`ToStandardOpticsUnits[input]`**

Removes units from supported `QuantityArray` data:

| Input | Numeric rows returned |
| --- | --- |
| Time trace | `{picoseconds, amplitude}` |
| Euler-form spectrum | `{terahertz, magnitude, phase}` |
| Dimensionless material spectrum | `{terahertz, value}` |
| Absorption-coefficient spectrum | `{terahertz, inverse centimeters}` |

For example, an absorption result from `EstimateMaterialParameter` becomes numeric `{THz, cm^-1}` rows:

```wl
absorption = EstimateMaterialParameter["\[Alpha]", spectrum, thickness];
numericAbsorption = ToStandardOpticsUnits[absorption];
```

Unsupported inputs are returned unchanged. Lists are processed elementwise.

## Acknowledgments 
The work of K.V. (@JerryI) was supported by the Deutsche
Forschungsgemeinschaft (DFG, German Research Foundation)-__TRR 360-492547816__

I am deeply grateful to Dr. Joachim Deisenhofer for countless discussions, endless support, and deep insights into THz spectroscopy.

## References 
We rely on well-known methods for material parameter extraction:

- Material parameter extraction for terahertz time-domain spectroscopy using fixed-point iteration. W. Withayachumnankul, B. Ferguson, T. Rainsford, S. P. Mickan, and D. Abbott. Photonic Materials, Devices, and Applications (2005) DOI: 10.1117/12.612946.
- Phase Retrieval in Terahertz Time-Domain Measurements: A “How To” Tutorial. Peter Uhd Jepsen. Journal of Infrared, Millimeter, and Terahertz Waves (2019) DOI: 10.1007/s10762-019-00578-0.

### Published Papers

We use these tools for processing THz data in our department Experimentalphysik V. Below is a list of our publications where we used this tool to process data:

- Optical magnetoelectric effect in the polar honeycomb antiferromagnet Fe₂Mo₃O₈. K. V. Vasin, A. Strinić, F. Schilberth, S. Reschke, L. Prodan, V. Tsurkan, A. R. Nurmukhametov, M. V. Eremin, I. Kézsmárki, and J. Deisenhofer. Phys. Rev. B 110, 054401 – Published 1 August 2024.
- Magnetization reversal through an antiferromagnetic state. Somnath Ghara, Evgenii Barts, Kirill Vasin, Dmytro Kamenskyi, Lilian Prodan, Vladimir Tsurkan, István Kézsmárki, Maxim Mostovoy & Joachim Deisenhofer. Nature Communications volume 14, Article number: 5174 (2023).


## License

GNU
