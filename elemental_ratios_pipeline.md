# Engineering Pipeline for Deriving Elemental Ratios from Spectroscopic Data: Methodology and Design Document

**Scope:** Orbital Near-Infrared (NIR) reflectance spectroscopy (0.4-3.5 $\mu$m, e.g., Moon Mineralogy Mapper M$^3$, Chandrayaan-2 IIRS, CRISM, OMEGA) with explicit treatment of **CLASS data characteristics.**

> **Terminological Clarification:** In flight terminology, CLASS = **Chandrayaan-2 Large Area Soft X-ray Spectrometer** (0.5-10 keV, X-ray fluorescence). It measures *atomic* $K_{\alpha}$ lines directly: Mg $K_{\alpha}$ 1.25 keV, Al $K_{\alpha}$ 1.49 keV, Si $K_{\alpha}$ 1.74 keV, Ca $K_{\alpha}$ 3.69 keV, Ca $K_{\beta}$ 4.01 keV, Fe $K_{\alpha}$ 6.40 keV. NIR spectrometers measure *mineralogical* crystal-field and vibrational absorptions, not elemental emission lines. This distinction is load-bearing for Section 2: **true $\text{Mg/Si}$, $\text{Al/Si}$, $\text{Ca/Si}$, $\text{Fe/Si}$ are natively CLASS-XRF products, and only mineral-chemical proxies thereof are derivable from NIR alone.** This pipeline is therefore designed as (a) NIR-only proxy pipeline, and (b) NIR+CLASS-XRF fusion pipeline for calibrated elemental ratios.

---

## 1. Raw Data Preprocessing & Calibration

Objective: Convert raw Detector Digital Numbers $DN(x,\lambda,t)$ to science-grade reflectance $R(\lambda) = I/F$ with trustworthy band center, depth, and area.

### 1.1 Instrument Calibration & Dark Current Correction

**What:**

1. Subtract dark frame $D(x,\lambda,T_{det})$ scaled to exposure time and detector temperature:

   $$DN_{1} = DN_{raw} - D(T_{det}) \cdot t_{exp} - B_{bias}$$
2. Flat-field correction: $DN_{2} = DN_{1} / F_{flat}(x,\lambda)$
3. Bad-pixel / drop-out masking + despiking (cosmic rays), linearity correction for saturation $> \sim$80% full-well.
4. Radiometric conversion: $L_{rad}(\lambda) = DN_{2} \cdot G(\lambda)$ using pre-flight + in-flight solar diffuser / lamp gain $G$.
5. Convert to $I/F$:

   $$R_{obs}(\lambda) = \frac{\pi \cdot L_{rad}(\lambda) \cdot d_{sun}^{2}}{F_{sun}(\lambda) \cdot \cos(i)}$$

   where $d_{sun}$ is heliocentric distance in AU, $F_{sun}$ is solar irradiance model (e.g., ASTM-E490, TSIS), $i$ is solar incidence angle from SPICE kernels + DEM.

**Why (physical/statistical):**

- HgCdTe / InGaAs NIR arrays have strong temperature-dependent dark current and $1/f$ bias drift. Failure to correct introduces wavelength-dependent slope that mimics continuum reddening and corrupts Band II area, biasing pyroxene chemistry.
- Flat-field non-uniformity and spectrometer `smile` cause cross-track radiometric striping. Without correction, mosaicked $\text{Mg/Si}$ proxy maps show artificial lineations.
- Linearity/saturation: saturated 1 $\mu$m shoulder flattens Band I minimum, shifting apparent band center longward $\rightarrow$ false Fe-enrichment.

**Implementation:** Temperature-interpolated dark library; per-orbit dark collects over night side; iterative sigma-clipped cosmic-ray filter in spectral-spatial domain; maintain Calibration Uncertainty Cube $\sigma_{cal}(\lambda)$.

### 1.2 Wavelength Calibration & Shift Correction

**What:**

1. Initial mapping $\lambda_0(p)$ from monochromator lab calibration (polynomial per spatial pixel).
2. In-flight refinement:
   - Atmospheric / solar Fraunhofer lines (for ground / Earth-observing): e.g., O$_2$ 0.762, 1.27 $\mu$m; solar Ca II.
   - For airless bodies: on-board laser / arc lines, or cross-correlation of observed telluric-free continuum with solar model; thermal-induced shift tracking via instrument housekeeping $T_{opt bench}$.
3. Correct smile (cross-track $\Delta\lambda$) and keystone via resampling to common grid with sinc / spline interpolation preserving integrated flux.
4. Optional sub-pixel shift: maximize cross-correlation $C(\delta) = \int S_{obs}(\lambda+\delta) \cdot S_{ref}(\lambda)d\lambda$ with $\delta$ typically $\pm$ 2-5 nm.

**Why:**

- 5 nm shift at 1 $\mu$m corresponds to $\sim$10 mol% error in olivine Fo# ($Fo = 100\cdot\text{Mg}/(\text{Mg+Fe})$) or $\sim$15 Wo error in pyroxene. Band center calibrations e.g., Adams (1974), Cloutis & Gaffey, Sanchez et al. have slopes $d\lambda/d\text{Fe} \approx$ 0.5-1.0 nm per mol% Fe.
- Thermal flexure in orbit causes seasonal 2-10 nm drift. Uncorrected, this creates spurious latitudinal compositional gradients.

**QC metric:** Report residual $\Delta\lambda_{RMS} < 0.3 \times$ spectral sampling (e.g., $<2$ nm for 10 nm sampling). Propagate into UQ.

### 1.3 Continuum Removal / Baseline Correction

**What:** Isolate absorption $R_c$-normalized features from scattering slope + thermal emission.

Pipeline stages:

1. **Thermal removal ($\lambda > 2.0$ $\mu$m for Moon/Mercury):** Model emitted radiance $B(T_{surf},\epsilon)$ using temperature from Diviner / thermal band or iterative fit:

   $$R_{true} = \frac{L_{rad} - \epsilon B(T)}{\frac{F_{sun}}{\pi d^2}\cos(i) - ...}$$

   Iterate $T$ to flatten 2.5-3.5 $\mu$m continuum or enforce Kirchhoff $\epsilon = 1-R$. Failure here mimics weak 2.8 $\mu$m OH/H$_2$O.
2. **Photometric normalization:** Lommel-Seeliger / Hapke / Akimov disk function to $i=e=g=30^{\circ}$ standard geometry.
3. **Continuum fitting per absorption complex:**
   - Convex hull (upper hull) over 0.6-2.6 $\mu$m tie points, or
   - Polynomial / log-linear continuum anchored at local maxima (e.g., 0.7, 1.5, 2.6 $\mu$m for mafics).
   - Normalized spectrum: $R_{n} = R / R_{c}$, Band Depth: $BD(\lambda_c) = 1 - R_{n}(\lambda_c)$
   - Equivalent Width: $EW = \int_{\lambda_1}^{\lambda_2} (1-R_n(\lambda))d\lambda$
   - Integrated Band Area: $IBA = \int (1-R_n)d\lambda$ for Band I, Band II.

**Why:**

- NIR continuum slope is controlled by space weathering (nanophase Fe), grain size, phase angle, not chemistry. Without removal, band depth is not comparable across scenes; $\text{Fe/Si}$ proxy would correlate with maturity (Is/FeO) rather than Fe.
- Choice of anchors is the dominant systematic for weak plagioclase 1.25 $\mu$m feature ($BD \sim$ 1-3%). Subjectivity must be tracked (see Sec. 3).

**Rule:** Never fit a single global polynomial across 0.4-3.0 $\mu$m. Fit locally per complex: 0.7-1.7 $\mu$m (Band I), 1.5-2.6 $\mu$m (Band II), 2.6-3.6 $\mu$m (hydration).

### 1.4 Atmospheric / Environmental Correction (if applicable)

- **Airless bodies (Moon, asteroids):** No telluric correction. Correct instead for: (a) exospheric scattering (negligible), (b) opposition surge and phase reddening, (c) stray light / scattered lunar dayside into nightside CLASS-XRF background.
- **Mars / Earth ground-truth:** Telluric CO$_2$ (1.43, 1.58, 2.0 $\mu$m), H$_2$O (1.4, 1.9 $\mu$m) removal via MODTRAN / ATREM ratio to bland standard (Olympus Mons / solar diffuser) observed at matched airmass:

  $$R_{surf} = R_{obs} / T_{atm}(aerosol,\tau, H_2O)$$
- **CLASS-XRF environment:** Solar state filter (only $>$B-class flare or quoted quiet-Sun flux from XSM/SXM), particle background (geotail, SEP events) subtraction, collimator response deconvolution.

**Why:** Martian 2.0 $\mu$m CO$_2$ sits directly on pyroxene Band II minimum. Residual 1% telluric creates false HCP/LCP classification $\rightarrow$ false $\text{Ca/Si}$ shift.

### 1.5 Signal-to-Noise Enhancement Preserving Line Shape

**Allowed:**

1. Savitzky-Golay (order 2-3, window 5-9 channels) — preserves 2nd moment (area) and minimum position. Verify $\Delta BD < 0.5\sigma_{noise}$.
2. Spectral binning only to Nyquist (no oversmoothing below instrument FWHM). For IIRS ($\sim$10 nm sampling, FWHM $\sim$20 nm), do not smooth $>15$ nm effective.
3. Spatial averaging with spectral homogeneity test (SAM angle $<2^{\circ}$) to boost SNR $\propto \sqrt{N}$ without mixing lithologies.
4. Wavelet (Daubechies-4) soft-thresholding for spikes; Minimum Noise Fraction (MNF) for cube-level denoising retaining first 10-20 eigenbands.
5. Mask, do not interpolate, channels with SNR $<20$ or inside deep telluric / thermal residuals.

**Forbidden:** Boxcar moving average $>3$ pixels, Gaussian smoothing broader than resolution, FFT low-pass with ringing, continuum-divided smoothing that clips minima.

**Gate:** Proceed only if after processing: SNR $>$ 50 at 1.5 $\mu$m continuum for Band I chemistry, SNR $>$ 30 at 2.0 $\mu$m for Band II; else flag as `non-quantitative, detection-only`.

> Pseudocode gate:
>
> ```text
> R, sigma = calibrate(DN, dark_lib, flat, geom, sun_model)
> lam = wavelength_calibrate(p, T_bench, smile_map)
> R_nothem, T = remove_thermal(R, lam)
> Rc = fit_continuum(R_nothem, anchors=[0.73,1.52,2.65um], method='hull')
> Rn, BD, EW = normalize(R_nothem, Rc)
> Rn_d = savitzky_golay(Rn, win=7, poly=3) if preserves_BD else Rn
> if SNR_continuum < threshold: flag QUALITY=low
> ```

---

## 2. Element and Ratio Feasibility (CLASS Data Context)

Core physics: NIR absorptions arise from **Fe$^{2+}$ crystal-field transitions, charge-transfer, and OH/H$_2$O vibrations**, not photoelectric $K$-shell ionization. Therefore NIR measures **mineral stoichiometry and Fe-site occupancy**, from which elemental ratios are *inferred* via mineral-mixing models, not directly counted.

### 2.1 Derivable (NIR Proxy) and Directly Derivable (CLASS-XRF) Ratios

| Ratio | NIR diagnostic | Physical basis | Reliability | CLASS-XRF direct lines |
|---|---|---|---|---|
| **$\text{Mg/Si}$ (mafic Mg#)** | Olivine Band I center 1.02-1.08 $\mu$m; Pyroxene Band I+II centers; BAR = Area(Band II)/Area(Band I) | Fe$^{2+}$ in M1/M2 octahedral sites: more Mg $\rightarrow$ shorter $\lambda_c$, smaller crystal-field splitting. King & Ridley: $\lambda_{ol} \approx 1.04 + 0.0012\cdot Fa$ [$\mu$m]. High-Mg orthopyroxene: BI $\sim$0.90, BII $\sim$1.85 $\mu$m; Fe-rich: BI $\sim$0.95-1.05, BII $\sim$2.0-2.3 $\mu$m. | **Semi-quantitative proxy** ($\pm$10-15 mol% Fo/En with lab cal). Requires ol/px ratio constraint. | **Direct, high reliability:** Mg $K_{\alpha}$/Si $K_{\alpha}$ flux ratio + fundamental-parameters inversion. Primary CLASS product. |
| **$\text{Fe/Si}$ / mafic abundance** | Integrated Band I+II depth/area; 1 $\mu$m $EW$; albedo + 2 $\mu$m depth | $BD \propto$ modal mafics $\times$ Fe$^{2+}$ abundance $\times$ grain size (Hapke). Low-Fe highlands: $BD_{1\mu m}<2$%; mare: 5-15%. | **Proxy for FeO wt%** ($\pm$1-2 wt% after Lucey/Gillis calibration) but degenerate with grain size, opaques, weathering. | **Direct:** Fe $K_{\alpha}$ 6.4 keV / Si. Needs strong solar flare; low quiet-Sun SNR. |
| **$\text{Ca/Si}$ (HCP vs LCP)** | Band II center + width: LCP (pigeonite/enstatite) BII $\sim$1.9-2.1 $\mu$m narrow; HCP (augite) BII $\sim$2.2-2.35 $\mu$m broad. Band I asymmetry. | Ca$^{2+}$ in M2 distorts chain silicate lattice $\rightarrow$ larger site $\rightarrow$ longer wavelength Fe$^{2+}$ transition. Adams, Cloutis calibrations: Wo% $\approx f(\lambda_{BI},\lambda_{BII})$. | **Mineralogical proxy only.** Bulk $\text{Ca/Si}$ requires plagioclase (Ca-bearing but nearly featureless) correction; NIR alone underestimates Ca in anorthosite. | **Direct:** Ca $K_{\alpha}$/Si. Good during flares. |
| **$\text{Al/Si}$ (feldspathic)** | Plagioclase 1.25-1.30 $\mu$m weak Fe$^{2+}$ band; absence of 1/2 $\mu$m bands + high albedo + Christiansen Feature (needs MIR); 2.0 $\mu$m spinel side-lobe check | Anorthite itself is NIR-transparent; traced by *lack* of mafic absorption + weak 1.3 $\mu$m shoulder. Empirical: $BD_{1\mu m} < 0.02$ + albedo $>0.25$ $\rightarrow$ FAN, $\text{Al/Si} \sim$0.6-0.9. | **Indirect / low precision.** Easily masked by $>$2% mafics or shock glass. Quantification requires fusion with XRF/GRS or MIR CF. | **Direct:** Al $K_{\alpha}$ 1.49 keV / Si. Primary CLASS product, but Al-Si line blending (ΔE $\sim$250 eV) demands high resolution ($<$200 eV FWHM) deconvolution. |
| **Hydration / OH proxies** | 2.8-3.0 $\mu$m $BD$, $EW$ | O-H stretch fundamental. Not elemental ratio but volatile indicator. | Quantitative only after rigorous thermal removal. | Not accessible to CLASS. |

**Method for NIR-to-ratio conversion (proxy chain):**

1. Extract parameters: $\lambda_{BI}$, $\lambda_{BII}$, $BD_I$, $BD_{II}$, BAR, albedo.
2. Unmix mineralogy: Linear / Hapke intimate mixture inversion against RELAB endmembers (olivine Fo$_{x}$, OPX En$_{x}$, CPX Wo$_{x}$, plagioclase An$_{x}$, glass, agglutinate) minimizing $\chi^2 = \sum (R_n-R_{model})^2/\sigma^2$.
3. Map mineral chemistry $\rightarrow$ bulk oxide via stoichiometric mass balance:

   $$\text{Mg/Si} = \frac{\sum_j f_j \cdot (\text{Mg}_j / M_j)}{\sum_j f_j \cdot (\text{Si}_j / M_j)}$$

   where $f_j$ modal fraction, $M_j$ formula mass. Report as *model-dependent* ratio.
4. If CLASS-XRF overlap exists: regress $R_{NIR-proxy}$ vs $(\text{X}/\text{Si})_{XRF}$ and apply empirical transfer $X/Si = a\cdot P_{NIR}+b$ with regional coefficients (mare vs highlands separately).

### 2.2 Non-Derivable Elements / Ratios from NIR Alone — Physical Reasoning

- **Si itself:** No NIR-active electronic or vibrational fundamental in 0.4-3.5 $\mu$m (Si-O stretch is at 8-12 $\mu$m MIR). Normalization to Si is therefore *assumed*, not measured, in NIR-only work. All NIR $\text{X/Si}$ are doubly model-dependent.
- **$\text{K/Si}$, $\text{Th/Si}$, $\text{U}$, $\text{Na}$, $\text{P}$, $\text{S}$:** Alkali / incompatible / volatile elements are hosted in feldspathic glass, KREEP mesostasis, or sulfides with no diagnostic NIR bands at orbital abundances (K$_2$O $<0.5$ wt%, Th $<5$ ppm). K-feldspar 2.4 $\mu$m overtone is far too weak and blended. Requires Gamma-Ray Spectrometer (K, Th, U) or XRF with low-energy sensitivity.
- **$\text{Ti/Si}$:** Ilmenite / Ti-rich pyroxene is spectrally bland (low albedo, blue slope, suppressed bands). Ti3+ charge transfer near 0.5 $\mu$m is non-unique vs nanophase Fe. Lucey UV/VIS ratio is a maturity-confounded proxy, not a line ratio.
- **$\text{Cr}$, $\text{Mn}$, $\text{Ni}$, REE:** Transition-metal spin-allowed bands overlap Fe$^{2+}$ and occur at $<$1000 ppm; below detection ($BD < \sigma$). REE $f$-$f$ narrow lines ($\sim$0.58, 0.80 $\mu$m Nd) require $>$100s ppm + hyperspectral SNR $>200$, not met orbitally.
- **Metallic Fe / S saturation:** Opaques (ilmenite, metal, sulfides) quench bands non-linearly (Hapke saturation). Beyond $\sim$15 wt% FeO + opaques, $BD$ saturates then *decreases* with more Fe — pipeline must flag $R_{0.75\mu m}<0.07$ as saturated, non-invertible.
- **Resolution limits:** Distinguishing Al $K_{\alpha}$ vs Si $K_{\alpha}$ in XRF, or olivine vs Fe-bearing glass in NIR (both broad 1 $\mu$m), requires $R=\lambda/\Delta\lambda >100$ and SNR thresholds rarely met in single pixels. Report upper limits, not values, when $\Delta\chi^2 <9$ between models.

**Bottom line:** Publish NIR-only $\text{Mg/Si}$, $\text{Fe/Si}$ as **mineral-chemical proxies (Mg#, FeO, BAR, Wo)** with explicit calibration; publish true elemental weight ratios **only after cross-calibration to CLASS-XRF / APXS / GRS** (Sec. 4). Do not report $\text{K/Si}$, $\text{Th/Si}$, $\text{Ti/Si}$ from NIR alone except as qualitative flags.

---

## 3. Uncertainty Quantification (UQ) & Error Propagation

### 3.1 Sources of Error (Error Budget Taxonomy)

1. **Photon / detector shot noise:** $\sigma_{shot} = \sqrt{N_e + N_{dark} + RN^2}/N_e$, propagated per channel to $\sigma_R(\lambda)$. Dominant in shadow, Band II, 3 $\mu$m.
2. **Calibration systematics:** Gain $G$ ($\sim$1-3%), solar model ($\sim$1%), photometric function ($\sim$2-5% in $R$), wavelength $\sigma_{\lambda}$ ($\sim$1-3 nm). Correlated across $\lambda$ — must use covariance, not diagonal.
3. **Thermal / photometric model residual:** $\sigma_{therm}(T\pm 5K)$ blows up $>2.3$ $\mu$m; track as inflated $\sigma_R$.
4. **Continuum placement subjectivity:** Anchor choice $\pm$20-50 nm and hull vs polynomial changes $EW$ by 5-20% for weak bands. Largest term for plagioclase $\text{Al/Si}$ proxy.
5. **Telluric / background residual:** 0.5-2% in $R_n$ inside correction bands.
6. **Smoothing bias:** Systematic shallowing $\Delta BD_{smooth}$; quantify by processing synthetic lines.
7. **Model (geological) error:** Endmember non-uniqueness, grain-size / weathering degeneracy — reported separately as *interpretation uncertainty*, not folded into measurement $\sigma$.

Maintain per-pixel vector $\sigma_R(\lambda)$ + effective correlation length $l_c$, plus $\sigma_{\lambda}$.

### 3.2 Propagation Methodology

**A. Analytical (fast, for band depths):**

For $BD = 1-R_b/R_c$:

$$\sigma_{BD}^2 = \frac{\sigma_{R_b}^2}{R_c^2} + \frac{R_b^2\sigma_{R_c}^2}{R_c^4} -2\frac{R_b}{R_c^3}\text{Cov}(R_b,R_c)$$

For $EW = \sum (1-R_i/R_{c,i})\Delta\lambda$:

$$\sigma_{EW}^2 = \mathbf{J}^T \mathbf{C}_R \mathbf{J}, \quad J_i = -\Delta\lambda/R_{c,i}$$

with full covariance $\mathbf{C}_R$ (diagonal + calibration block). For ratio $r = EW_2/EW_1$ or $BAR$:

$$\left(\frac{\sigma_r}{r}\right)^2 = \left(\frac{\sigma_1}{EW_1}\right)^2+\left(\frac{\sigma_2}{EW_2}\right)^2 -2\frac{\text{Cov}_{12}}{EW_1EW_2}$$

Band-center error from parabolic fit minimum: $\sigma_{\lambda_c} \approx \frac{\text{FWHM}}{\text{SNR}\sqrt{N_{pts}}}$ plus systematic $\sigma_{\lambda,cal}$ added in quadrature.

Transfer to chemistry via calibration slope $S = d(\text{chemistry})/d\lambda$:

$$\sigma_{Fo} = |S|\sqrt{\sigma_{\lambda_c}^2+\sigma_{cal-curve}^2}$$

**B. Monte Carlo perturbation (baseline for final products — recommended):**

```text
for k in 1..N=1000:
  R_k = R + N(0, C_R)          # correlated noise realisation
  lam_k = lam + N(0, sigma_lam)
  anchors_k = perturb(anchors, ±1 chan)  # continuum subjectivity
  Rc_k = fit_continuum(R_k, anchors_k, random choice hull/poly)
  Rn_k, params_k (BI,BII,BD,EW,BAR), chem_k = extract(R_k,Rc_k)
report median, 16/84 percentiles; full posterior for X/Si
```

- Use $N\ge 500$ (1000 for flagship maps). Vectorize; cost is trivial vs science return.
- For Hapke unmixing, MCMC (emcee) over modal fractions $f_j$ + grain size $D$ yields posterior $P(\text{Mg/Si}|R)$ directly, marginalizing degeneracies.
- Detection limit: $BD_{lim} = 3\cdot\sigma_{BD,MC}$; if $BD < BD_{lim}$, report 3$\sigma$ upper limit on $EW$ and flag ratio as limit.

Cross-check analytical vs MC on 1% random pixels; require agreement within 20% or investigate covariance misspecification.

### 3.3 Reporting Standards (Data Product Requirements)

Per-pixel / per-footprint product must contain:

- Value + $1\sigma$ standard error + 95% confidence interval (2.5/97.5 percentiles from MC), e.g., $\text{Mg/Si} = 0.42 \pm 0.05\;(1\sigma),\;95\%\,\text{CI}\,[0.32,0.52]$.
- Quality flag: `good / marginal (SNR-limited) / saturated / thermal-contaminated / detection-limit`.
- Provenance: continuum method, anchors, smoothing, calibration version, $N_{MC}$.
- Correlation note when ratio shares Si denominator: provide covariance or warn against independent-error assumption in plots.
- Global metadata: detection limits map, calibration-curve RMS (e.g., Fo $\pm$8), statement: *NIR-proxy* vs *XRF-calibrated elemental ratio*.
- Visualization: maps show value + separate relative-uncertainty ($\sigma/r$) layer; scatter plots show 68% error ellipses, not bare points.

Never report ratios to more significant figures than $\sigma$ justifies (e.g., $0.42\pm0.05$, not $0.4187$).

---

## 4. Validation, External Data, and Cross-Checks

### 4.1 Ground Truth / Correlative Flight Datasets

1. **CLASS-XRF + XSM (Chandrayaan-2):** Primary cross-validator for $\text{Mg/Si}$, $\text{Al/Si}$, $\text{Ca/Si}$, $\text{Fe/Si}$ at $\sim$12.5 km footprints. Co-register NIR proxy (100 m) averaged to CLASS footprint; regress and require $R^2>0.6$, slope $0.8-1.2$ after outlier rejection, else recalibrate. Use only flare-state CLASS (`XSM flux > threshold`, particle background quiet).
2. **APXS / ChemCam / SuperCam (landers/rovers):** Point-source absolute truth ($\pm$2-5% oxides). Use Apollo / Luna sample stations and Chang'E sites as tie points for NIR transfer function.
3. **Gamma-Ray Spectrometers (Lunar Prospector GRS, Kaguya GRS, Mars Odyssey GRS):** Independent $\text{Fe/Si}$, $\text{Th}$, $\text{K}$ at coarse ($\sim$100-300 km) scale; validates regional means, breaks NIR Fe vs maturity degeneracy.
4. **Higher-resolution orbital spectrometers:** M$^3$ vs IIRS cross-comparison; CRISM / OMEGA for Mars; Diviner CF (7-9 $\mu$m) for $\text{Al/Si}$ (plagioclase) and SiO$_2$ polymerization independent of NIR. Require band-center agreement within combined $\sigma_{\lambda}$.
5. **Topography / photometry:** LOLA / TMC DEM for incidence correction audit; Clementine UVVIS / LROC WAC for FeO-TiO$_2$ maturity cross-check.

Protocol: hold out 20% of overlap footprints as blind test; publish bias $\mu = \langle r_{NIR}-r_{XRF}\rangle$ and scatter RMS. If $|\mu| > 1\sigma_{combined}$, pipeline fails validation — do not ship.

### 4.2 Laboratory Mineral Spectra

- **RELAB (Brown), USGS splib, HOSERLab:** Endmembers for unmixing: olivine Fo$_{0-100}$, OPX En$_{30-95}$, CPX Wo$_{10-50}$, plagioclase An$_{90-98}$ + shocked variants, ilmenite, agglutinates, impact glass at multiple grain sizes (45-250 $\mu$m) and temperatures.
- Procedure: convolve lab $R=500$ spectra to flight sampling + add flight noise; verify pipeline recovers known chemistry within quoted $\sigma$. Test specifically: (a) olivine-pyroxene mixtures BAR vs ol/(ol+px), (b) plagioclase + 2% pyroxene masking test for $\text{Al/Si}$ false negatives, (c) space-weathered (laser-irradiated) series to quantify weathering bias on $\text{Fe/Si}$.
- Maintain versioned spectral library + optical constants (Hapke $n,k$) with DOI; any library change triggers full revalidation.

### 4.3 Sanity Checks & Closure Tests

1. **Synthetic injection-recovery:** Inject Gaussian / Hapke-modeled absorptions of known $\lambda_c$, $BD$, $EW$ into bland highland / noise spectra across SNR grid (10-200). Require: $|\Delta\lambda_c| < \sigma_{\lambda}$, $|\Delta BD|/BD < 10\%$ for $BD>3\sigma$, linear recovery slope $1.0\pm0.05$. Map saturation curve $BD_{obs}$ vs $BD_{true}$ to define linearity limit.
2. **Saturation / non-linearity test:** Forward-model intimate mixtures to high FeO + opaques; confirm pipeline flags quenching rather than reporting spuriously low $\text{Fe/Si}$.
3. **Continuum-sensitivity test:** Reprocess full scene with alternate continuum (hull vs 2nd-order poly, anchor $\pm1$ channel). If derived $\text{X/Si}$ shifts $>1\sigma_{MC}$, inflate reported $\sigma$ to encompass method spread.
4. **Spatial closure:** Adjacent-orbit overlap must agree within combined errors; mosaics must not show track-boundary steps in ratio maps (stripe test via cross-track median profile).
5. **Geochemical closure:** Derived oxides (from NIR+XRF fusion) must sum to $100\pm5$ wt% (allowing H, C exclusion); $\text{Mg\#}$ vs $\text{An\#}$ must fall on plausible igneous trends (e.g., FAN vs mare fields). Outliers trigger manual review.
6. **Null tests:** Process calibration-target / empty-space spectra; must return $BD$ consistent with zero. Process mature vs fresh crater rays of same unit; chemistry must be invariant while slope varies — else weathering correction has leaked into ratio.

**Deployment rule:** Pipeline version is accepted only after passing (a) lab recovery, (b) injection-recovery linearity, (c) blind XRF/Ground-truth regression, and (d) closure sum, with all metrics published in product documentation. Any change to calibration, continuum, or library increments major version and requires reprocessing + revalidation.

---

**Summary Deliverable per Observation:** L1 calibrated $I/F$ + $\sigma$ cube → L2 continuum-removed parameters ($BD$, $EW$, $\lambda_c$, BAR + MC posteriors) → L3 mineral-proxy ($\text{Mg\#}$, FeO, Wo, FAN flag) → L4 fused elemental ratios $\text{Mg/Si}$, $\text{Al/Si}$, $\text{Ca/Si}$, $\text{Fe/Si}$ with CI, detection limits, quality flags, and validation pedigree. NIR alone yields L3; L4 requires CLASS-XRF / GRS anchoring.
