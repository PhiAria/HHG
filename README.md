# HHG Argon Macroscopic Model

A modular Python model for estimating gas-phase high-harmonic generation (HHG) in argon using a cylindrically symmetric reduced macroscopic framework.

## What this repository does

This code estimates harmonic yield, per-harmonic pulse energy, and photons per pulse for odd harmonics between configurable orders (default 11 to 51). It is designed as a physically interpretable surrogate model for:

- gas-phase HHG in Ar,
- parameter scans over laser and gas conditions,
- source estimates for later beamline and monochromator modeling,
- studying qualitative effects of phase, ionization, and absorption.

The implementation lives in `HHG.py` and is structured around dataclasses plus helper functions for laser, gas, dipole, absorption, and calibration models.

## Scientific model

The code uses a reduced cylindrical model in coordinates \((r,z,t)\), assuming axial symmetry.

### 1. Driving laser

The driver is not numerically propagated through the medium. Instead, it is prescribed as a Gaussian beam in space and a Gaussian pulse in time:

\[
I(r,t;z) = I_0 \left(\frac{w_0}{w(z)}\right)^2 \exp\!\left(-\frac{2r^2}{w(z)^2}\right) \exp\!\left(-4\ln 2\, \frac{t^2}{\tau_{\mathrm{FWHM}}^2}\right)
\]

with

\[
w(z) = w_0 \sqrt{1 + (z/z_R)^2}, \qquad z_R = \pi w_0^2/\lambda_0.
\]

This captures focusing, transverse intensity variation, and pulse duration, but not self-focusing, plasma defocusing, reshaping, or depletion of the fundamental field.

### 2. Pulse energy, peak power, and peak intensity

The peak intensity is derived from average power and repetition rate rather than hardcoded:

\[
E_{\mathrm{pulse}} = \frac{P_{\mathrm{avg}}}{f_{\mathrm{rep}}}
\]

For a Gaussian temporal envelope:

\[
P_{\mathrm{peak}} = \frac{E_{\mathrm{pulse}}}{\tau_{\mathrm{FWHM}}\sqrt{\pi/(4\ln 2)}}
\]

For a Gaussian spatial beam:

\[
I_0 = \frac{2P_{\mathrm{peak}}}{\pi w_0^2}.
\]

This makes the model convenient for experimental users who know average power, repetition rate, pulse duration, and focused waist.

### 3. Gas density

The neutral argon density is estimated from the ideal gas law:

\[
N_0 = \frac{P}{k_B T}.
\]

Along the propagation axis the code supports:

- top-hat jet profile,
- Gaussian gas profile.

This gives a spatially varying source density along \(z\).

### 4. Cutoff estimate

The code uses the standard semiclassical cutoff law:

\[
E_{\mathrm{cut}} = I_p + 3.17 U_p
\]

with the ponderomotive energy approximated as

\[
U_p[\mathrm{eV}] \approx 9.337\times 10^{-14} \, I[\mathrm{W/cm^2}] \, \lambda[\mu\mathrm{m}]^2.
\]

This cutoff is used only as a guide and to apply a soft spectral rolloff in the surrogate dipole model.

### 5. Ionization model

Ionization is currently modeled with a surrogate saturation law:

\[
\eta(I) = \min\left[\left(\frac{I}{I_{\mathrm{sat}}}\right)^\gamma, \, \eta_{\max}\right].
\]

The neutral fraction entering the HHG source is:

\[
N_{\mathrm{neutral}} = N_0 (1-\eta).
\]

This captures the qualitative fact that increasing intensity eventually depletes neutrals and suppresses efficient generation, but it is not a microscopic tunnel-ionization model such as ADK or PPT.

### 6. Microscopic dipole surrogate

The code does not solve the TDSE or Lewenstein integral. Instead it uses a phenomenological dipole amplitude of the form:

\[
A_q(I,\eta) \propto I^{p_q} \exp(-\beta \eta) \times R_q
\]

where:

- \(p_q\) increases mildly with harmonic order,
- \(\exp(-\beta\eta)\) suppresses emission at high ionization,
- \(R_q\) is a soft rolloff near the semiclassical cutoff.

This is useful for scans and trend analysis, but it is not an ab initio single-atom response.

### 7. Dipole phase and trajectories

The source includes separate short- and long-trajectory contributions.

The dipole phase is modeled as:

\[
\phi_{\mathrm{dip}} = \alpha_q I.
\]

The phase coefficient \(\alpha_q\) is scaled with harmonic order and chosen separately for short and long trajectories. The full source phase also includes:

- Gouy phase,
- wavefront curvature phase,
- dipole phase.

The local source term is built schematically as:

\[
P_q \propto N_{\mathrm{neutral}} A_q \left(e^{i\phi_s} + w_l e^{i\phi_l}\right)
\]

where \(\phi_s\) and \(\phi_l\) are short- and long-trajectory phases.

This is one of the strongest features of the code, because it lets the model express physically meaningful spatial and phase-structure effects rather than using only scalar conversion efficiency.

### 8. Macroscopic buildup

For each harmonic order, the exit-plane field is assembled by coherent integration over the gas target:

\[
E_q(r,t) \propto \int P_q(r,z,t) \, e^{-\frac{z_{\mathrm{exit}}-z}{2L_{\mathrm{abs},q}}} \, dz.
\]

This means the model includes:

- coherent interference of sources along propagation,
- gas-density weighting,
- absorption of generated XUV light in the medium.

However, it does **not** solve a paraxial propagation equation for the generated harmonic field, so diffraction and reshaping of the XUV during propagation are not yet treated self-consistently.

### 9. Absorption

Absorption is included through an effective absorption length:

\[
\exp\left[-\frac{z_{\mathrm{exit}}-z}{2L_{\mathrm{abs},q}}\right].
\]

The current implementation uses a power-law scaling of absorption length with harmonic order rather than tabulated argon photoabsorption data.

### 10. Conversion to physical photons

The field integration produces a dimensionless surrogate yield:

\[
Y_q = \int \! \int |E_q(r,t)|^2 \, 2\pi r\, dr\, dt.
\]

That surrogate is mapped to physical harmonic pulse energy using a calibration constant:

\[
U_q = s\, Y_q\, E_{\mathrm{pulse}}.
\]

Then photons per pulse are estimated via

\[
N_q = \frac{U_q}{\hbar q \omega_0}.
\]

This is a practical engineering step, but it is also the main empirical normalization in the repository.

## Strengths of the repository

### Good physical insight per line of code

The code contains a lot of real HHG intuition in a compact form:

- Gaussian focusing geometry,
- gas density profile,
- cutoff scaling,
- ionization-induced suppression,
- short/long trajectory interference,
- dipole phase sensitivity,
- absorption-limited buildup.

### Clear modular structure

The code is organized into:

- dataclasses for configuration,
- reusable helper functions,
- a simulator class,
- separate summary and plotting methods.

This makes it easy to upgrade individual pieces without rewriting the whole program.

### Useful for scanning and design studies

It is well-suited for exploring:

- waist dependence,
- average power / rep-rate changes,
- gas pressure scans,
- gas-jet position shifts,
- approximate spectral content up to the 51st harmonic.

### Good basis for beamline modeling

The outputs are already close to what a monochromator model needs:

- per-harmonic photon number,
- per-harmonic pulse energy,
- selected-harmonic near-field structure.

That makes it a good starting point for a later TDCM throughput model.

## Scientific limitations

### The single-atom physics is not first-principles

The model does not compute the microscopic dipole from:

- TDSE,
- SFA/Lewenstein,
- quantitative rescattering theory.

Instead it uses a phenomenological amplitude and phase law. This is fine for trends, but not for predictive spectroscopy.

### Ionization is not ADK/PPT

The ionization model is an intensity-based saturation proxy. It does not correctly capture:

- sub-cycle ionization timing,
- wavelength dependence of tunnel rates,
- barrier suppression behavior,
- charge-state dynamics.

### No self-consistent IR propagation

The driving field is prescribed, not propagated. Therefore the model omits:

- plasma defocusing,
- Kerr self-focusing,
- neutral dispersion of the fundamental,
- pulse reshaping in the gas,
- depletion of the driver.

These effects can matter strongly near optimal phase-matching conditions.

### No explicit phase-mismatch decomposition

The code includes phase effects implicitly through Gouy, curvature, and dipole phase, but it does not explicitly solve for

\[
\Delta k = \Delta k_{\mathrm{neutral}} + \Delta k_{\mathrm{plasma}} + \Delta k_{\mathrm{geo}} + \Delta k_{\mathrm{dipole}}.
\]

As a result, it is not yet a true phase-matching model in the strict macroscopic HHG sense.

### No XUV paraxial propagation

The model integrates a source to an exit plane but does not solve propagation of the harmonic envelope in \(z\) with diffraction. Therefore:

- far-field divergence is only indirect,
- wavefront evolution is incomplete,
- collection efficiency into downstream optics is not yet reliable.

### Absolute photon numbers are calibration-dependent

The estimated photon numbers depend on the empirical `CalibrationConfig.scale` parameter. That means the code currently provides:

- good relative trends,
- plausible order-of-magnitude outputs if calibrated,
- but not fully predictive absolute flux from first principles.

## Recommended interpretation of results

Use current outputs primarily for:

- comparative studies,
- source optimization trends,
- testing sensitivity to geometry and gas parameters,
- building intuition for which harmonic orders are likely accessible.

Do **not** yet interpret the output as a high-accuracy quantitative prediction of absolute flux without calibration to experiment or a more microscopic model.

## Recommended next scientific upgrades

If the goal is realistic source estimation before modeling a TDCM monochromator, the most valuable upgrades are:

1. **Replace surrogate ionization with ADK**  
   This improves the neutral/plasma balance and makes intensity scans more trustworthy.

2. **Use tabulated argon absorption data**  
   Replace the power-law absorption length with energy-dependent Ar photoabsorption.

3. **Add explicit phase-mismatch terms**  
   Introduce neutral, plasma, geometric, and dipole contributions to \(\Delta k\).

4. **Add XUV far-field propagation**  
   Compute divergence, source size, and acceptance into downstream optics.

5. **Optionally propagate the fundamental self-consistently**  
   Needed if plasma defocusing or strong medium effects become important.

6. **Later: replace surrogate dipole with SFA/Lewenstein or lookup tables**  
   This would make spectral and phase trends more physically predictive.

## Repository structure

```text
HHG.py      Main model, helpers, configuration, simulation, and plots
README.md   Project overview and scientific notes
```

## How to run

Install dependencies:

```bash
pip install numpy matplotlib
```

Run:

```bash
python HHG.py
```

The script will:

- simulate odd harmonics from `q_min` to `q_max`,
- print a summary of laser, gas, cutoff, and photon estimates,
- display diagnostic plots.

## Key outputs

- `photons_per_pulse` for each harmonic order
- `pulse_energy_J` per harmonic
- `photon_energy_eV`
- `selected_exit_field` for a chosen harmonic
- diagnostic maps of driver intensity and effective ionization

## Intended future direction

This repository is especially well-positioned to evolve into:

1. an HHG source estimator,
2. then an XUV beam transport model,
3. then a TDCM monochromator transmission / resolution model.

That is a sensible workflow because the current code already produces the spectral and source-side quantities needed upstream of a monochromator.
