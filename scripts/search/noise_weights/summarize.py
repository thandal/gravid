#!/usr/bin/env python3
"""Assemble the J0437 search results into one small summary.

Reads the fixed-orientation search (``scripts/j0437_upper_limit.py
--noise-weights running-median``) and the F-statistic search
(``scripts/j0437_f_statistic.py``), checks that they agree where they overlap,
and recomputes every recovery amplitude exactly from the saved sums rather
than by interpolating on an amplitude grid.
"""
import json

import numpy as np
from scipy.stats import chi2, kstest, ncx2

from _common import HERE, ROOT, cplx, statistic_coefficients

RUN = ROOT / "results/j0437_search"
# LVK (2025), ApJ 983, 99, Table 2, PSR J0437-4715, 95% upper limits on h0.
PUBLISHED = dict(bayesian=7.2e-27, f_statistic=9.0e-27, five_n_vector=5.8e-27,
                 spin_down_limit=7.95e-27)


def amplitude(threshold, a, b, c, quantile=0.95):
    """Amplitude above which `quantile` of injections stay above threshold.

    The statistic for injection j at amplitude h is a_j h^2 + b_j h + c."""
    disc = b * b - 4 * a * (c - threshold)
    h = np.where(disc > 0, (-b + np.sqrt(np.maximum(disc, 0))) / (2 * a), 0.0)
    return float(np.quantile(np.maximum(h, 0), quantile))


def expected(threshold, dof, a):
    grid = np.geomspace(2e-27, 2e-25, 600)
    eff = np.array([np.mean(ncx2.sf(threshold, dof, a * h * h)) for h in grid])
    return float(np.interp(0.95, eff, grid))


def block(dof, a, b, c, controls, interval):
    one_pct = float(chi2.isf(0.01, dof))
    return dict(
        degrees_of_freedom=dof,
        target=float(c), target_probability=float(chi2.sf(c, dof)),
        controls_at_or_above_target=int(np.sum(controls >= c)),
        control_mean=float(controls.mean()),
        control_ks_probability=float(kstest(controls, "chi2", args=(dof,)).pvalue),
        loudest_control=float(controls.max()),
        probability_of_a_louder_maximum=float(
            1 - (1 - chi2.sf(controls.max(), dof)) ** len(controls)),
        h0_95=dict(
            above_loudest_control=amplitude(controls.max(), a, b, c),
            above_loudest_control_bootstrap_90pct=interval,
            at_1pct_false_alarm=amplitude(one_pct, a, b, c),
            at_1pct_false_alarm_expected_over_noise=expected(one_pct, dof, a),
            above_observed_target=amplitude(c, a, b, c)))


fixed_run = json.loads((RUN / "result.json").read_text())
f_run = json.loads((RUN / "f_statistic.json").read_text())

coefficients = statistic_coefficients(f_run)
a4, b4, c4 = coefficients["f_statistic"]
a2, b2, c2 = coefficients["fixed_orientation"]
ys = coefficients["y_signal"]
controls4 = np.array(f_run["controls"]["two_f"])
controls2 = np.array(f_run["controls"]["fixed_orientation"])

# The two searches share the fixed-orientation statistic; they must agree.
assert np.isclose(c2, fixed_run["twoF_fixed"], rtol=1e-6)
assert np.allclose(controls2, fixed_run["background"]["central_2F"], rtol=1e-6)
assert np.allclose(ys[:, 0], cplx([row[0] for row in fixed_run["injection"]["Ys"]]), rtol=1e-5)
assert fixed_run["search_par_sha256"] == f_run["search_par_sha256"]
assert fixed_run["covariance_sha256"] == f_run["covariance_sha256"]

noise = fixed_run["noise_weights"]
captured = np.array(f_run["injection"]["captured_power_fraction"])
summary = dict(
    inputs=dict(fixed_orientation=str((RUN / "result.json").relative_to(ROOT)),
                f_statistic=str((RUN / "f_statistic.json").relative_to(ROOT))),
    noise_weights=dict(
        parameters=noise["parameters"],
        equivalent_noise_bandwidth_hz=noise["equivalent_noise_bandwidth_hz"],
        per_detector=noise["per_detector"],
        samples_total=fixed_run["n_samp"],
        samples_kept=sum(v["kept"] for v in noise["per_detector"].values())),
    fixed_orientation=block(2, a2, b2, c2, controls2,
                            f_run["injection"]["fixed_orientation"]["loudest_control"]["h0_95_interval"]),
    f_statistic=block(4, a4, b4, c4, controls4,
                      f_run["injection"]["f_statistic"]["loudest_control"]["h0_95_interval"]),
    fixed_orientation_search_as_reported=dict(
        h0_95=fixed_run["injection"]["h0_95"],
        interval=fixed_run["injection"]["sampling_uncertainty"]["h0_95_interval"],
        note="first linear crossing on that search's 16-point amplitude grid"),
    fixed_orientation_share_of_signal_power=dict(
        median=float(np.median(a2 / a4)), fifth_percentile=float(np.quantile(a2 / a4, 0.05))),
    # Injected orbits are drawn from the timing covariance; the F-statistic
    # is evaluated at the one central orbit.
    power_lost_with_one_orbit=dict(maximum=float((1 - captured).max()),
                                   median=float(np.median(1 - captured))),
    spin_down_limit=fixed_run["astro"]["h0_spindown"],
    published_lvk_o4a=PUBLISHED)
(HERE / "results/summary.json").write_text(json.dumps(summary, indent=1) + "\n")
for name in ("fixed_orientation", "f_statistic"):
    s = summary[name]
    print(f"{name}: target {s['target']:.3f} (probability {s['target_probability']:.3f}, "
          f"{s['controls_at_or_above_target']} of 40 controls at or above); controls mean "
          f"{s['control_mean']:.2f}, loudest {s['loudest_control']:.2f}")
    for key, value in s["h0_95"].items():
        print(f"   {key}: {value}")
print("power lost with one orbit:", summary["power_lost_with_one_orbit"])
