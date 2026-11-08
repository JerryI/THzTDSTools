# Time-Domain Terahertz Spectroscopy Tools 🧰

A small library for material parameter extraction from time-domain THz spectra, written in Wolfram Language, C, and OpenCL.

*Developed at University of Augsburg, Germany 🇩🇪*


## Data format and units

The public processing functions use `QuantityArray` data:

| Domain | Rows | Units |
| --- | --- | --- |
| Time trace | `{time, amplitude}` | `{"Picoseconds", 1}` |
| One-sided frequency spectrum | `{frequency, magnitude, phase}` | `{"Terahertz", 1, 1}` |

Spectrum magnitudes are nonnegative and phases are numeric radians. The phase column can store values outside `[-Pi, Pi]`, so an unwrapped phase is preserved. Frequency-dependent material results have `{frequency, value}` rows; absorption coefficient `\[Alpha]` has inverse-centimeter units. Frequencies and phases must follow the package's Fourier convention, in which a delayed pulse has positive linear phase. Thickness, time-window positions, and frequency bounds are `Quantity` values.

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

`DeconvolveTraces` returns a time-domain transfer response. Its `TD2FD` spectrum is the regularized complex sample/reference transmission. Keep the linear phase when unwrapping if the spectrum will be inverted to time or used to estimate refractive index; `UnwrapPhase` does this by default.

## Function reference

### Traces and Fourier transforms

**`LinearProjectionsToCircular[trace1, trace2, {angle1, angle2}]`** reconstructs two circular-polarization time traces from measurements at two linear-analyzer angles. The measured projections are treated as `Cos[angle] Ex + Sin[angle] Ey`; the output is proportional to `Sqrt[2] (Ex ± Hilbert[Ey])`, preserving the earlier ±45-degree amplitude convention. Bare numbers are radians; `45 Degree` and angular `Quantity` values also work. The two angles cannot differ by a multiple of 180 degrees. Equal-length lists of traces are converted pairwise.

**`ResampleTraces[trace]`** interpolates onto a uniform time grid using the trace's smallest step. With a list of traces, it uses the smallest step across them and a common grid spanning all inputs. Values outside an input trace are zero. Output times are in picoseconds.

**`TD2FD[trace]`** returns a one-sided spectrum with DC and, by default, without the Nyquist bin. `"IncludeNyquist" -> True` keeps Nyquist and pads an odd-length resampled trace with one zero sample. A list of traces is processed elementwise.

**`FD2TD[spectrum]`** reconstructs a real, even-length trace starting at time zero. By default it interpolates onto a uniform frequency grid and assumes a zero Nyquist component. Use `"IncludeNyquist" -> True` when the spectrum contains both DC and Nyquist bins. Unwrapped phase values are accepted. A list of spectra is processed elementwise.

**`DeconvolveTraces[sample, reference]`** resamples both time traces to a common grid, removes their means, and returns a dimensionless lag trace. The frequency-domain transfer function is regularized by `"Threshold" -> alpha`, a positive relative power-spectrum parameter with default `10^-5`. FFT lags are circular: a negative-lag response appears near the end of the output rather than at a negative time. The intended response is causal, with the sample signal arriving no earlier than the reference.

### Frequency selection and filtering

**`FrequencyDomainConfidenceInterval[input]`** estimates a usable band from a Gaussian fit to spectral power. It accepts a time trace, spectrum, or list of either and returns an `Interval` with inverse-centimeter endpoints. If the fit is unsuitable, it falls back to the first frequency and a frequency near the midpoint.

**`FrequencySpectrumTrim[input, interval]`** keeps only frequency bins in `interval`, whose endpoints may be in terahertz or inverse centimeters. Without `interval`, it uses `FrequencyDomainConfidenceInterval[input]`. Time input is transformed with `TD2FD` first; output is always a spectrum. Lists are processed elementwise. An interval containing no bins gives `$Failed`.

**`BandPassFilter[input, center, Q]`** applies a Gaussian filter centered at a positive frequency or wavenumber `Quantity`. The positive dimensionless quality factor `Q` defaults to 1. Its amplitude response is `Exp[-4 Log[2] ((f - center)/(center/Q))^2]`, with full width at half maximum `center/Q`. Spectrum input retains frequencies and phases. Time input returns a filtered time trace with its starting time preserved; an odd sample count may gain one padded point. Lists are processed elementwise.

**`ApodizeTrace[input, position, width]`** applies `HannWindow` over a full-support width centered at `position`. Both are time `Quantity` values and width must be positive. `"WindowFunction" -> window` replaces the window with a function of `(time - position)/width` on `[-1/2, 1/2]`; output is zero outside that interval. Time input keeps its grid. Spectrum input must contain DC and is converted to time, windowed, and transformed back; `position` is then relative to the reconstructed time axis starting at zero, and the output phase can be wrapped. For a spectrum containing Nyquist, set `"IncludeNyquist" -> True` (default `False`). Lists are processed elementwise.

### Phase unwrapping

`TD2FD` stores the principal Fourier phase in the range from $-\pi$ to $\pi$. A sampled wrapped phase does not, by itself, tell us how many full $2\pi$ turns occur between frequency bins. The usual nearest-neighbor unwrap therefore becomes unreliable when phase changes rapidly or transmission has a deep minimum. This library offers a Kramers–Kronig (KK) guided method that uses the measured magnitude to help choose those turns.

**`UnwrapPhase[input, interval]`** returns a complete spectrum, changing phase inside the selected band and preserving bins outside it. Omitting `interval` uses `FrequencyDomainConfidenceInterval[input]`; an explicit interval may have terahertz or inverse-centimeter endpoints. Time input is Fourier transformed first. The methods are:

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

   $$\Delta\phi_i = \Delta\phi_i^{\mathrm{FT}} + 2\pi\,\operatorname{Round}\!\left(\frac{\Delta\phi_i^{\mathrm{KK}}-\Delta\phi_i^{\mathrm{FT}}}{2\pi}\right).$$

   It accumulates these steps starting from the first measured phase. Thus the KK phase guides branch selection; it does not replace the Fourier phase.
4. Restores the pulse-arrival phase unless `"RemoveLinearPhase" -> True` was requested.

`"KKGuidedMultiscale"` adds a correction near deep transmission minima. It compares the difference between guided and KK phases on both sides of each candidate minimum at several window sizes, then applies a $2\pi$ branch shift when the comparisons agree. The KK estimate is approximate for finite-band data, so neither method guarantees a unique phase where the signal is weak or the frequency sampling is too coarse.

`"RemoveLinearPhase" -> False` is the default and retains the pulse-arrival phase term. Setting it to `True` removes that term inside the interval, so the result no longer represents the original time delay. It has no effect with `Method -> None`. The KK-only phase is estimated from magnitude and need not invert to the original trace.

### Material properties

**`AveragedRefrectiveIndex[response, thickness]`** uses the peak of the absolute time-domain response to estimate `n0 = 1 + c lag/thickness`. Spectrum input is converted with `FD2TD` first. Thickness must be a positive length `Quantity`; the estimate assumes a causal response.

**`EstimateMaterialParameter[parameter, response, thickness]`** accepts a deconvolved time trace or a one-sided spectrum. A spectrum supplied for index estimation should already have unwrapped phase; a time trace is unwrapped internally with the Basic method and its linear phase retained. `EstimateMaterialParameter[]` lists the supported parameters:

| Parameter | Result |
| --- | --- |
| `"n0"` | Peak-delay average refractive index. |
| `"n"` | Frequency-dependent refractive index. |
| `"\[Kappa]"` | Frequency-dependent extinction coefficient. |
| `"\[Alpha]"` | Absorption coefficient in inverse centimeters. |
| `"T"` or `"Transmission"` | Corrected power transmission, clipped to `[0, 1]`. |
| `"Complex n"` | `n + I \[Kappa]`. |
| `"FP Reduction"` | Dimensionless before/after FP-suppression score. |

`Method -> "Basic"` (default) makes a single-pass estimate, suitable when internal reflections are negligible. With transfer magnitude `A`, phase `phi`, thickness `L`, gain `g`, and phase correction `phiC`, the package uses `n = 1 + c (phi + phiC)/(2 Pi f L)` and `\[Kappa] = -c Log[g A]/(2 Pi f L)`. It estimates `\[Alpha] = (Log[(4 n0/(n0 + 1)^2)^2] - 2 Log[g A])/L`, including the Fresnel transmission correction based on `n0`. `T = Clip[(g A)^2, {0, 1}]`. At zero frequency, `n` and `Complex n` use `n0`, while `\[Alpha]` and `\[Kappa]` are zero. `"Gain" -> g` defaults to 1; `"PhaseCorrection" -> phiC` defaults to 0 radians.

`Method -> "Iterative"` runs the slab solver and needs a spectrum containing DC. The historical spelling `"Interative"` is also accepted. A method list such as `Method -> {"Iterative", "Iterations" -> 30}` can hold solver settings, or use the top-level options below.

| Option | Default | Meaning |
| --- | --- | --- |
| `"NKIterations"` | `30` | Index/extinction iterations. |
| `"FPIterations"` | `8` | Fabry–Perot iterations when enabled. |
| `"SolveFP"` | `False` | Include internal reflections in the iterative solve. |
| `"MovingAverage"` | `True` | Smoothing setting for the solver; also accepts a nonnegative integer. |
| `"Target"` | `"Automatic"` | Try GPU first, then CPU fallback; `"GPU"` and `"CPU"` are also accepted. |

`"FP Reduction"` always runs the iterative solver with FP removal. It compares the pre- and post-removal `f \[Kappa]` spectra near the predicted round-trip delay: above 1 means suppression, below 1 means worsening, and exactly 1 can mean no measurable reduction **or** no detectable FP peak. This is a diagnostic, not a thickness estimate; a material feature with the same periodicity can distort it. The chosen band must resolve the echo delay or the result is `$Failed`.

Set **`"FPInterval"` as a top-level option**, outside `Method`. Its default `Automatic` uses `FrequencyDomainConfidenceInterval`, `All` uses all positive bins, and an explicit frequency `Interval` restricts the score to a reliable measurement band. A fine thickness scan is useful when the suppression peak is narrow:

```wl
Table[{d, EstimateMaterialParameter["FP Reduction", spectrum,
    Quantity[d, "Millimeters"], "Gain" -> 1.2,
    "FPInterval" -> Interval[{Quantity[1.5, "Terahertz"],
                              Quantity[2.5, "Terahertz"]}]]},
  {d, 0.62, 0.65, 0.0025}]
```

**`DeconvolveFP[response, thickness]`** removes Fabry–Perot echoes using iterative `n` and `\[Kappa]` estimates. A time trace is Fourier transformed and phase-unwrapped internally, then returned as a time trace. A spectrum must contain DC and is assumed to have already-unwrapped phase; it returns a spectrum. Thickness is a positive length `Quantity`. Options are `"NKIterations" -> 30`, `"FPIterations" -> 8`, `"MovingAverage" -> True`, `"Target" -> "Automatic"`, `"Gain" -> 1.0`, and `"PhaseCorrection" -> 0`. Gain and phase correction are included in the output, so pass it to `EstimateMaterialParameter` with default gain and phase options.

**`GenerateMaterialResponse[n, thickness, {interval, step}]`** generates the transmission through a slab in the package's `{frequency, magnitude, phase}` `QuantityArray` format. `n` is a `Function` or a symbol defined as a function; it receives each frequency as a plain number in terahertz and returns the dimensionless complex refractive index `n + I \[Kappa]`. Thickness is a positive length `Quantity`. The interval is a frequency `Quantity` containing an `Interval`, for example `Quantity[Interval[{0, 4}], "Terahertz"]`; the positive step is another frequency `Quantity`. Terahertz and inverse centimeters are accepted for both. The grid includes the upper endpoint when it lies on a step. The generated phase retains the slab delay and can exceed `Pi`.

`"SolveFP" -> False` is the default and uses the single-pass slab response `H = (4 n/(n + 1)^2) Exp[I 2 Pi f L (n - 1)/c]`, with `f` in terahertz, `L` in centimeters, and `c` in centimeters per picosecond. `"SolveFP" -> True` divides this response by `1 - ((n - 1)/(n + 1))^2 Exp[I 4 Pi f L n/c]` to include internal reflections. This sign follows the positive delay phase used by `TD2FD` and the iterative material solver.

**`GenerateMaterialResponse[n, thickness, reference]`** evaluates the same transfer function at the frequencies of `reference`. With the default `"ConvolveWithReference" -> True`, it multiplies the transfer function by the reference. A frequency spectrum returns a frequency spectrum with its existing unwrapped phase retained and the model phase added. A time trace is Fourier transformed, multiplied by the model, and returned as a time trace with the original starting time. Set `"ConvolveWithReference" -> False` to use the reference only for its sampling grid and output domain: a spectrum returns the bare transfer spectrum, while a time trace returns the impulse response starting at zero lag. The option is ignored for the `{interval, step}` form.

Time input uses the full one-sided spectrum including Nyquist; an odd-length trace gains one padded sample. Its real inverse transform discards any imaginary DC or Nyquist component of the modeled spectrum. Convolution is circular over the finite FFT window, so allow enough time after the pulse for delayed FP echoes. The most recent reference transform is reused across repeated calls with the same `QuantityArray`.

```wl
nModel[f_] := 1.8 + 0.02 I;
thickness = Quantity[0.635, "Millimeters"];
response = GenerateMaterialResponse[nModel, thickness,
  {Quantity[Interval[{0, 4}], "Terahertz"], Quantity[0.01, "Terahertz"]},
  "SolveFP" -> True];
simulatedTrace = GenerateMaterialResponse[nModel, thickness, referenceTrace,
  "SolveFP" -> True];
```

### Utilities

**`TimeTraceQ[input]`** returns `True` for a numeric two-column time `QuantityArray`, `False` for a numeric three-column Euler-form spectrum, and `$Failed` otherwise. The expected units are compatible with the formats shown above.

**`ToStandardOpticsUnits[input]`** removes units from supported `QuantityArray` data: time traces become numeric picosecond/amplitude rows, Euler-form spectra become terahertz/magnitude/phase rows, and two-column material spectra become terahertz/value rows. Unsupported inputs are returned unchanged. Lists are processed elementwise.

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

## Supported Platforms

- WLJS Notebook
- Wolfram Mathematica
- WolframScript (Wolfram Engine)

## License

GNU
