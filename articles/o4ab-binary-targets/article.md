# Reassessing binary pulsars for the O4ab narrowband search

Nathaniel Fairfield

## Abstract

The O4ab narrowband search for continuous gravitational waves omitted PSRs
J1045−4509 and J1737−0811 because of orbital uncertainty, and J1231−1411
because the ephemeris used in the search lacked required orbital parameters.
We fit the three orbits to MeerKAT pulse arrival times and test the
predictions against later observations made during O4b. With one orbital
template per pulsar, the bound on coherent power loss stays below 2% for
all three across the timing and noise models tested. The bound allows any
nonnegative weights on the O4ab time grid. At the later radio epochs, the measured
phase residuals imply less than 0.10% loss for each pulsar. All three could therefore have been included with one orbit
per target at a 3% loss tolerance. For J1400−1431, one template covers the
five published orbital errors, but the longer timing solution gives no
uncertainty for orbital drift. For J1745−0952, the historical errors also
fit within one template. The MeerKAT observations we could test end in 2021.

## 1. Which targets need more than one orbit?

The LIGO–Virgo–KAGRA (LVK) narrowband search in the first two parts of the
fourth observing run (O4ab) used a radio timing solution to correct for each
binary pulsar's motion. For J0437−4715, J1045−4509, J1737−0811 and
J1745−0952, orbital errors gave an estimated coherent observing time between
one and 90 days. LVK excluded these targets and divided the observations of
J1400−1431 into three segments, analysing the first and third. LVK omitted
J1231−1411 because the ephemeris used for that search lacked some required
orbital parameters ([LVK 2026, §4.1](#ref-lvk)). The earlier O4a targeted search had used Nançay timing
for both J1231 and J1745 and reported upper limits for both. Neither appeared
in its narrowband results ([LVK 2025, §4.2 and Tables 2–3](#ref-lvko4a)).

The selection rule asks whether the uncertainty in the Doppler shift exceeds
half a Fourier bin ([Mirasola et al. 2024](#ref-mirasola)). A small periodic phase error can cross that threshold on
a long observation while leaving most of the signal power intact. We ask how
much power a search would lose by using one orbit throughout O4ab. The count
refers only to orbits. The search still covers a band in gravitational-wave
frequency and spin-down.

[Dhurandhar & Vecchio (2001, Tables I and III)](#ref-dhurandhar) estimated
that J1045 would require 16 filters over one year at 3% mismatch. The S5
targeted search included J1045 and estimated power loss by drawing sets of
timing parameters from their joint uncertainties
([Abbott et al. 2010, §3 and Table 1](#ref-s5)). J1045's orbital period
uncertainty fell from 8.6 ms in the catalogue used for the 2001 estimate
to 7 μs in a later fit to 24 years of Parkes observations
([Reardon et al. 2021, Table A3](#ref-reardon)). We revisit the
filter count with MeerKAT timing and apply the 3% tolerance to each pulsar.

Extending the [companion study of J0437](https://github.com/thandal/gravid/tree/main/articles/j0437-one-orbital-template)
across O4ab gives a largest fractional loss of $4.5\times10^{-8}$ among
1,280 tested orbits with uniform time weights.

<div style="break-inside: avoid">

For a phase difference $\delta\phi(t_j)$ between the signal orbit and the
template, the fraction of coherent power retained at times $t_j$ is

$$
\mathcal C=\left|\frac{\sum_j w_j\exp[i\,\delta\phi(t_j)]}{\sum_j w_j}\right|^2 .
$$

</div>

A constant phase offset leaves the power unchanged. If every phase lies within
$r<\pi/2$ of a common offset, its projection along that direction is at
least $\cos r$. The power loss $1-\mathcal C$ is then at most $\sin^2 r$
for any nonnegative time weights with a positive sum. This gives a way to
bound the loss without knowing when the detectors were most sensitive.

## 2. Three orbits tested against O4b radio observations

We fit J1045, J1737 and J1231 using the 4.5-year MeerKAT Pulsar Timing
Array release ([Miles et al. 2025](#ref-miles)). The observations for these pulsars end in 2023. For each
target we fit the pulse times of arrival (TOAs) with TEMPO2, allowing an
independent dispersion correction each day ([Hobbs et al. 2006](#ref-tempo2)). The ELL1 model describes a
nearly circular orbit by its period, projected semimajor axis, epoch of
ascending node and two components of eccentricity
([Lange et al. 2001, Appendix A](#ref-lange)). We also fit rates of
change for the period, projected axis and both eccentricity components.
For J1737 we first convert the published DD orbit to ELL1, accounting for
the effect of periastron advance on orbital longitude and eccentricity.

For J1231, 15 years of Green Bank and Nançay timing already give all five
ELL1 orbital parameters and an orbital period derivative of
$(7.6\pm0.5)\times10^{-13}$. The observations end just before O4b. From the
Shapiro delay, Cromartie et al. inferred a companion mass of
$0.23^{+0.09}_{-0.06}\,M_\odot$ and an
inclination of $79.8^{+3.5}_{-4.7}$ degrees, using a prior on companion
mass from binary evolution. These are 68.3% credible intervals
([Cromartie et al. 2026, Table 1 and §4.3](#ref-cromartie)).

Our MeerKAT reference fit fixes the companion mass at $0.3\,M_\odot$ and
$\sin i=0.95$. The binary mass function then gives a pulsar mass of
$2.66\,M_\odot$. We test the effect of this choice by refitting at 24 other
combinations spanning pulsar masses of $0.8$–$3.0\,M_\odot$ and $\sin i$
from 0.2 to 0.9999. We compare all 25 fits with the same template.

### How much power can the orbital errors cost?

We use generalized least squares to retain the correlations between the
nine orbital parameters and the other fitted parameters. We compare four
descriptions of the timing noise. The first rescales the TOA errors. The
second also allows noise shared by channels in the same observation.
The other two add correlations between observations, with a power-law
spectrum or an exponential decay ([Coles et al. 2011](#ref-coles)).

We simulate 200 datasets per pulsar under its fitted power-law noise model
and refit the noise in each trial. The nominal 95% joint orbital regions
contain the input orbit in 94.5%, 90.5% and 93.0% of trials for J1045,
J1737 and J1231. We enlarge each region using a one-sided upper estimate of
the simulated 95th percentile. We vary the noise amplitudes, spectral
slope and low-frequency cutoff on a grid. We retain noise models within
13.277 of the best fit in minus twice the restricted log likelihood for
each low-frequency cutoff. Where adjacent amplitudes straddle this threshold,
we also test an amplitude just inside it. For the power-loss calculation
we raise the scale applied to the TOA errors to at least 1.1. We use the
same retained noise models for all 25 J1231 fits. The coverage of each
orbital region depends on the assumed noise.

As a check, we fit the low and high radio-frequency halves separately. The
orbital differences give nominal $p$-values of 0.68, 0.49 and 0.39 in the
same target order. Predicting each quarter of the timing span from the
other three gives the largest discrepancy for J1045 in the final quarter
(nominal orbital $p=0.058$ under the power-law noise model).

We calculate the phase at the upper edge of the search band, 1.001 times
twice the rotation frequency. At each time we take the largest excursion
allowed by the joint orbital errors, including the offset of the fit from
the reference template. We compare all fits on 8,192 equally spaced times
across O4ab, including gaps. For the fit giving the largest bound for each
pulsar, we repeat the calculation on up to 131,072 times.

<div style="break-inside: avoid">

*Table 1. Coherent power lost with one orbit per pulsar. Column 2 gives the
largest bound for the noise models tested, rounded up. It applies to any
nonnegative weights on the O4ab time grid. Column 4 uses the measured phase
residuals at the later radio epochs, with equal weight per epoch.*

| Pulsar | O4ab loss bound | O4b epochs (TOAs) | Loss at O4b epochs |
|---|---:|---:|---:|
| J1045−4509 | 0.20% | 17 (544) | 0.0069% |
| J1737−0811 | 1.0% | 16 (512) | 0.063% |
| J1231−1411 | 2.0% | 17 (544) | 0.086–0.094% |

</div>

The bounds before rounding are 0.1919%, 0.9917% and 1.9228% for the limiting
fits. Increasing the grid from 32,768 to 131,072 times changes these bounds
by less than 0.004% in relative terms. We checked the linear approximation
to the orbital phase against full TEMPO2 delays at 258 points on each
limiting error ellipsoid. The largest difference was below
$6\times10^{-6}$ rad. We also fit the second spin-frequency derivative
for all three pulsars and free parallax for J1045 and J1737. J1231 already
has a fitted parallax. With these extra parameters the largest bounds on
the 8,192-point grid are 0.208%, 0.991% and 1.864%, respectively, using
the reference mass and inclination for J1231.

### What did the later radio data show?

[Ashton et al. (2015)](#ref-ashton) used radio measurements of the Crab's
timing noise to estimate the coherent power lost in a gravitational-wave
search. Here we use the
earlier binary timing solutions to predict separate O4b observations from
the [Pulsar Portal](#ref-portal).

We convert each timing residual to gravitational-wave phase at twice the
rotation frequency, weight the channels within an epoch
by inverse variance and give each epoch equal weight. Column 4 uses the
original channel phases, including any change in dispersion since the
training observations. For J1231 the range covers all 25 fits.
For scale, introducing a 0.3-ms error in orbital delay through the projected
axis or the orbital epoch raises the sampled loss to 9–15% for J1045 and
34–37% for J1737.

To test for an error at the orbital period, we first remove dispersion at
each O4b epoch. We fit the channel residuals to a constant plus a term
proportional to the inverse square of the radio frequency. The constant
estimates the residual at infinite frequency. We propagate its uncertainty,
the uncertainty in the earlier timing fit and the noise correlations between
epochs. We then fit a quadratic trend in spin phase and test for a remaining
sine and cosine at the orbital period.

<div style="break-inside: avoid">

*Table 2. O4b residuals after removing dispersion, using each reference
orbit. The two $\chi^2$ columns use the power-law noise model, with degrees
of freedom in parentheses. The last column gives the nominal orbital-test
$p$-values across all four noise models.*

| Pulsar | Prediction $\chi^2$ (dof) | After spin trend $\chi^2$ (dof) | Orbital $p$ |
|---|---:|---:|---:|
| J1045−4509 | 49.67 (17) | 10.89 (14) | 0.684–0.811 |
| J1737−0811 | 18.24 (16) | 10.65 (13) | 0.069–0.125 |
| J1231−1411 | 79.83 (17) | 9.42 (14) | 0.952–0.988 |

</div>

The J1045 and J1231 predictions fit poorly in absolute phase. A quadratic
trend accounts for most of the discrepancy in both. J1737 leaves room for
an orbital residual of a few microseconds. These observations check the
predictions at the radio epochs. The bound over O4ab comes from the timing fits.

## 3. Two targets with incomplete timing information

**J1400−1431.** The Nançay timing solution released by LVK reports 1,372
TOAs from September 2012 through December 2024. It gives errors for five
orbital parameters and holds their secular derivatives fixed
([LVK timing release](#ref-lvkdata)). In the
first-order ELL1 delay, allowing each parameter to vary by three quoted
standard deviations bounds the power loss throughout O4ab at 0.293%, for
any time weights. Even at eight standard deviations the bound is 2.07%.
One template covers the five published errors.

An independent LOFAR dataset, released under the name J1400−1438, provides
1,468 TOAs through January 2020 ([Donner et al. 2020](#ref-donner)). Holding the Nançay orbit fixed while
refitting spin, astrometry, dispersion and station offsets describes the
LOFAR data nearly as well as freeing all nine orbital parameters
($F=0.535$ for 9 and 1,282 degrees of freedom, nominal $p=0.85$). LOFAR
alone measures the period change as $-0.024\pm1.190$ ms/yr. That uncertainty
accumulates over the years between the LOFAR observations and O4. A later
LOFAR study followed J1400 through August 2022 to measure dispersion from
the solar wind ([Susarla et al. 2024, Table 4](#ref-susarla)).

With dispersion free at each epoch, the Donner fit multiplies the quoted
TOA errors by 1.50 and gives negligible correlated noise. Twelve templates
cover the nine-parameter orbital region for this noise covariance at 3% loss
with uniform O4ab time weights. We bound the phase
variance across the linearized region and check full TEMPO2 delays for
4,151 trial orbits. The largest measured loss is 2.958%.

Three independent 49.8-MHz measurements from the Long Wavelength Array
([Kumar et al. 2025](#ref-kumar)) give O4a and O4b timing
residuals of $+298$ and $-241$ μs relative to 2019, modulo the
gravitational-wave period, assuming GPS-referenced station times. Three
epochs cannot separate orbital drift from spin, propagation and clock
errors. To count templates for the longer Nançay solution with drift free,
we would need the arrival times or the joint errors from such a fit.

**J1745−0952.** The five orbital errors published by
[Janssen et al. (2010, Table 1)](#ref-janssen) give a
three-standard-deviation bound of 2.26% under the same ELL1 calculation.
The available Portal observations span eight epochs in 2019–2021. The
folding ephemeris is labelled barycentric dynamical time (TDB). Using its
numbers as labelled gives 91% sampled loss at those epochs. Interpreting
them as barycentric coordinate time (TCB) gives 0.40%. Converting the whole
TCB model to TDB with PINT ([Luo et al. 2021](#ref-pint)) reproduces the TEMPO2 residuals to within 18 ns
after removing a constant. The historical data favour the TCB interpretation.
The [public release for LVK's O4a targeted search](https://dcc.ligo.org/G2402400/public)
contains upper limits. Its Nançay arrival times were unavailable for an O4
check of the J1745 prediction.

## 4. Conclusion

For every timing noise model tested, one orbital template per pulsar keeps
the bound on coherent power loss below 3% for J1045, J1737 and J1231 over
the full O4ab time grid.
The separate O4b radio observations support that choice for all three.
A strain search over the frequency and spin-down bands of these targets
would test for signals that the published O4ab narrowband search did not
cover.

## Data, software and acknowledgments

The [gravid repository](https://github.com/thandal/gravid) contains the
analysis code, fitted orbits and numerical results. Instructions for
reproducing Table 1 are in `articles/o4ab-binary-targets/README.md`.

We acknowledge the Pulsar Portal and its publication in preparation
([Bailes et al. 2026](#ref-bailes2026), [Portal data-use policy](https://pulsars.org.au/data-usage/)).
The October 2026 downloads contain MeerTime pulsar timing array products
for J1045, J1737 and J1745 ([Bailes et al. 2018](#ref-bailes), [2020](#ref-bailes2020),
[Miles et al. 2025](#ref-miles)) and RelBin products for J1231
([Kramer et al. 2021](#ref-kramer)). Their folding ephemerides and clock
corrections can differ from the packaged MeerKAT timing release.

<div style="break-inside: avoid">

## References

<a id="ref-s5"></a>Abbott, B. P., et al. (2010). *Searches for Gravitational Waves from Known Pulsars with Science Run 5 LIGO Data*. [ApJ 713, 671–685](https://doi.org/10.1088/0004-637X/713/1/671).

</div>

<a id="ref-ashton"></a>Ashton, G., Jones, D. I., and Prix, R. (2015). *Effect of timing noise on targeted and narrow-band coherent searches for continuous gravitational waves from pulsars*. [Physical Review D 91, 062009](https://doi.org/10.1103/PhysRevD.91.062009).

<a id="ref-bailes"></a>Bailes, M., et al. (2018). *MeerTime—the MeerKAT Key Science Program on Pulsar Timing*. [PoS MeerKAT2016, 011](https://doi.org/10.22323/1.277.0011).

<a id="ref-bailes2020"></a>Bailes, M., et al. (2020). *The MeerKAT Telescope as a Pulsar Facility: System verification and early science results from MeerTime*. [PASA 37, e028](https://doi.org/10.1017/pasa.2020.19).

<a id="ref-bailes2026"></a>Bailes, M., et al. (2026). Pulsar Portal publication, in preparation.

<a id="ref-coles"></a>Coles, W., Hobbs, G., Champion, D. J., Manchester, R. N., and Verbiest, J. P. W. (2011). *Pulsar timing analysis in the presence of correlated noise*. [MNRAS 418, 561–570](https://doi.org/10.1111/j.1365-2966.2011.19505.x).

<a id="ref-cromartie"></a>Cromartie, H. T., et al. (2026). *Shapiro Delay Measurements from Fifteen Years of PSR J1231−1411 Radio Observations*. [ApJ 999, 115](https://doi.org/10.3847/1538-4357/ae3d95).

<a id="ref-dhurandhar"></a>Dhurandhar, S. V., and Vecchio, A. (2001). *Searching for continuous gravitational wave sources in binary systems*. [Physical Review D 63, 122001](https://doi.org/10.1103/PhysRevD.63.122001).

<a id="ref-donner"></a>Donner, J. Y., Verbiest, J. P. W., Tiburzi, C., et al. (2020). *Dispersion measure variability for 36 millisecond pulsars at 150 MHz with LOFAR*. [A&A 644, A153](https://doi.org/10.1051/0004-6361/202039517). [Timing models, templates, TOAs and dispersion series](https://doi.org/10.5281/zenodo.4290012).

<a id="ref-tempo2"></a>Hobbs, G. B., Edwards, R. T., and Manchester, R. N. (2006). *TEMPO2, a new pulsar-timing package—I. An overview*. [MNRAS 369, 655–672](https://doi.org/10.1111/j.1365-2966.2006.10302.x).

<a id="ref-janssen"></a>Janssen, G. H., Stappers, B. W., Bassa, C. G., Cognard, I., Kramer, M., and Theureau, G. (2010). *Long-term timing of four millisecond pulsars*. [A&A 514, A74](https://doi.org/10.1051/0004-6361/200911728).

<a id="ref-kramer"></a>Kramer, M., et al. (2021). *The relativistic binary programme on MeerKAT: science objectives and first results*. [MNRAS 504, 2094–2114](https://doi.org/10.1093/mnras/stab375).

<a id="ref-kumar"></a>Kumar, P., Taylor, G. B., Stovall, K., Dowell, J., and White, S. M. (2025). *A Multifrequency Census of 100 Pulsars below 100 MHz with LWA: A Systematic Study of Flux Density, Spectra, Timing, Dispersion, Polarization, and Its Variation from a Decade of Observations*. [ApJ 982, 132](https://doi.org/10.3847/1538-4357/adb97a). Public folded observations: [LWA Pulsar Archive](https://lda10g.alliance.unm.edu/PulsarArchive/), retrieved 17 October 2026.

<a id="ref-lange"></a>Lange, Ch., et al. (2001). *Precision timing measurements of PSR J1012+5307*. [MNRAS 326, 274–282](https://doi.org/10.1046/j.1365-8711.2001.04606.x).

<a id="ref-lvko4a"></a>LIGO Scientific Collaboration, Virgo Collaboration and KAGRA Collaboration (2025). *Search for Continuous Gravitational Waves from Known Pulsars in the First Part of the Fourth LIGO–Virgo–KAGRA Observing Run*. [ApJ 983, 99](https://doi.org/10.3847/1538-4357/adb3a0). [O4a data products](https://dcc.ligo.org/G2402400/public).

<a id="ref-lvk"></a>LIGO Scientific Collaboration, Virgo Collaboration and KAGRA Collaboration (2026). *Narrowband Searches for Continuous Gravitational Waves from known Pulsars in the First Two Parts of the Fourth LIGO–Virgo–KAGRA Observing Run*. [ApJ 1005, 221](https://doi.org/10.3847/1538-4357/ae77fc).

<a id="ref-lvkdata"></a>LVK (2026). *O4ab narrowband data products*, [LIGO-T2600052](https://dcc.ligo.org/LIGO-T2600052/public), released analysed-target ephemerides.

<a id="ref-pint"></a>Luo, J., Ransom, S., Demorest, P., et al. (2021). *PINT: A Modern Software Package for Pulsar Timing*. [ApJ 911, 45](https://doi.org/10.3847/1538-4357/abe62f).

<a id="ref-miles"></a>Miles, M. T., Shannon, R. M., Reardon, D. J., et al. (2025). *The MeerKAT Pulsar Timing Array: the 4.5-yr data release and the noise and stochastic signals of the millisecond pulsar population*. [MNRAS 536, 1467–1488](https://doi.org/10.1093/mnras/stae2572).

<a id="ref-mirasola"></a>Mirasola, L., Leaci, P., Astone, P., et al. (2024). *New semicoherent targeted search for continuous gravitational waves from pulsars in binary systems*. [Physical Review D 110, 123043](https://doi.org/10.1103/PhysRevD.110.123043).

<a id="ref-portal"></a>The Pulsar Portal (2026). [MeerTime pulsar timing data](https://pulsars.org.au/), retrieved 16–18 October 2026. The input manifest identifies the files used.

<a id="ref-reardon"></a>Reardon, D. J., et al. (2021). *The Parkes pulsar timing array second data release: timing analysis*. [MNRAS 507, 2137–2153](https://doi.org/10.1093/mnras/stab1990).

<a id="ref-susarla"></a>Susarla, S. C., et al. (2024). *Exploring the time variability of the solar wind using LOFAR pulsar data*. [A&A 692, A18](https://doi.org/10.1051/0004-6361/202450680).
