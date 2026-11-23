#!/usr/bin/env python3
"""Coherent O4a search for PSR J0437-4715 with a fixed-orientation filter.

Reads the 102 reduced (heterodyned) files of both detectors, 246 days, and
computes:
  * the matched-filter statistic 2|<m,d>|^2/<m,m>, summed coherently over
    files and detectors -- every file was heterodyned against the same base
    ephemeris, so the series share a phase reference;
  * orbital templates placed by the phase metric on the timing covariance
    (for this pulsar the metric places one);
  * the same statistic at off-source frequency offsets (controls);
  * recovery of injected signals of random orientation over an amplitude
    grid, the amplitude at which 95% exceed the loudest control, and its
    conversion to ellipticity, sensitivity depth and spin-down ratio.

Sample weights are chosen with ``--noise-weights``; see
``gravfold/noise_weights.py``.

Stages (``--stage``):
  measure    load everything, report the template count and time per template.
  onsource   alias for upperlimit.
  upperlimit target, controls and injections.

Cheap knobs for a dry run: --max-nodes, --n-offsets, --n-inj, --h0.
"""
from __future__ import annotations
import os, sys, glob, json, time, argparse, re
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np

from gravfold.upper_limit import (heterodyne_statistic_response,
                                  linear_grid_crossing,
                                  paired_grid_bootstrap_crossing)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Reduced strain: 102 files heterodyned against one base ephemeris.  BASE_PAR
# must be the converted par those files were heterodyned against, or the
# coherent sum is incoherent; it is verified byte-equal to each file's
# embedded par.
HETDIR = os.path.join(HERE, "results", "campaign_het_v2")
BASE_PAR = os.path.join(HERE, "results", "campaign_het_v2", "base_cwinpy.par")
OUTDIR = os.path.join(HERE, "results", "j0437_de440_qr")
COVARIANCE_ARTIFACT = os.path.join(
    HERE, "gravfold", "resources", "j0437_de440_binary_covariance.json")
from gravfold.orbital_metric import covariance_sha256
from gravfold.cw_phase_cache import PhaseCache
from gravfold.j0437_model import SEARCH_PAR, KEYS, binary_mean, apply_binary, load_covariance, waveform_metric

import pint.logging  # quiet PINT if pulled in transitively
try:
    pint.logging.setup(level="ERROR")
except Exception:
    pass

# J0437 physical constants for the astrophysical interpretation (Reardon et
# al. 2024, ApJL 971 L18).  The spin-down limit uses the intrinsic fdot:
# J0437's observed spin-down is dominated by the Shklovskii (kinematic) term,
# computed below from the measured proper motion and distance:
#   Pdot_kin = mu^2 * d * P / c  with  mu = 140.913 mas/yr, d = 156.98 pc
#            = 4.360e-20 s/s  (76% of Pdot_obs = 5.729e-20)
#   Pdot_int = 1.369e-20 s/s -> fdot_int = -4.13e-16 Hz/s -> h0_sd = 7.9e-27,
# close to the ~7.8e-27 quoted by the O3 searches (arXiv:2208.03932).  The
# remaining Galactic-acceleration term is ~1e-22 s/s and is neglected.
DIST_PC = 156.98            # Reardon+2024 DR3 distance [pc] (D_shk 156.98(15);
                            # the paper's NICER prior; D_parallax = 155.4(12))
_SEARCH_VALUES = {line.split()[0]: line.split()[1] for line in open(SEARCH_PAR)
                  if line.split() and not line.lstrip().startswith("#")}
F_SPIN = float(_SEARCH_VALUES["F0"])
F1_OBS = float(_SEARCH_VALUES["F1"])
PM_MAS_YR = float(np.hypot(float(_SEARCH_VALUES["PMRA"]), float(_SEARCH_VALUES["PMDEC"])))
I_ZZ = 1e38                 # canonical moment of inertia [kg m^2]
G = 6.674e-11; C = 2.998e8; PC = 3.086e16

_MU = PM_MAS_YR * 1e-3 * (np.pi / 180.0 / 3600.0) / (365.25 * 86400.0)  # rad/s
_P = 1.0 / F_SPIN
_PDOT_OBS = -F1_OBS / F_SPIN**2
_PDOT_KIN = _MU**2 * (DIST_PC * PC) * _P / C                # Shklovskii
FDOT = -(_PDOT_OBS - _PDOT_KIN) * F_SPIN**2                 # INTRINSIC [Hz/s]


from gravfold.j0437_campaign import load_chunk as _load_chunk, verify_selection

# How samples are weighted.  "file-variance" uses one plain variance per
# five-day file and drops no samples.  "running-median"
# uses gravfold.noise_weights: the transient before each stride end dropped,
# a running-median noise estimate, samples above a power cut dropped, and the
# filter's equivalent noise bandwidth, so that 2F is chi-squared with two
# degrees of freedom in Gaussian noise.
# Set once in main() before any worker is forked.
NOISE_WEIGHTS = "file-variance"


class MultiDatasetDemodulator:
    """Coherent orbit->2F over many heterodyned chunks and both detectors.

    2F = 2 |sum_i <d_i, m_i>_w|^2 / sum_i <m_i, m_i>_w, weight w_i = 1/var_i
    (per-chunk empirical variance -- studentization). The sum runs over every
    chunk of every detector: a single complex amplitude coherent across the run.
    Node templates are the SAME for every data realization, so `precompute_nodes`
    caches them once and injection recovery is then just dot products.
    """

    def __init__(self, files, base_par, studentize=True, search_par=SEARCH_PAR,
                 noise_weights=None):
        from cwinpy.parfile import PulsarParameters
        from cwinpy.signal import HeterodynedCWSimulator
        from gravfold.cwinpy_adapter import (
            _ensure_template_amplitude, _apply_orbit, theta_from_par)
        self._PP = PulsarParameters
        self._apply_orbit = apply_binary
        self._ensure_amp = _ensure_template_amplitude
        self.base_par = base_par
        self.search_par = str(search_par)
        base = _ensure_template_amplitude(PulsarParameters(base_par))
        self._base = base
        self.central_theta = binary_mean(PulsarParameters(self.search_par))
        self.search_frequency = 2*float(PulsarParameters(self.search_par)["F0"])
        # CONSOLIDATE each detector's chunks into ONE concatenated series + ONE
        # simulator. 102 barycentering simulators -> 2 collapses per-worker RAM
        # from ~5 GB to ~2 GB (the earlier --nproc 8 run OOM-thrashed on 102
        # sims/worker) and cuts model() calls 51x. Studentization is preserved
        # by carrying a PER-SAMPLE inverse-variance array: each original 5-day
        # chunk contributes 1/var over its own samples.
        by_det = {}
        for f in files:
            det, times, data, _ = _load_chunk(f)
            v = float(np.var(data)) if studentize else 1.0
            by_det.setdefault(det, []).append(
                (int(times[0]), times, data, 1.0 / max(v, 1e-300)))
        self.chunks = []            # now one entry PER DETECTOR
        self.noise_weights = noise_weights or NOISE_WEIGHTS
        self.noise_summary = {}
        for det, subs in sorted(by_det.items()):
            subs.sort(key=lambda s: s[0])                  # ascending GPS
            times = np.concatenate([s[1] for s in subs])
            data = np.concatenate([s[2] for s in subs])
            invvar = np.concatenate([np.full(s[1].size, s[3]) for s in subs])
            subs.clear()                    # release the per-file copies
            file_variance_snr = float(np.sum(invvar)**2
                                      / np.sum(invvar**2 * np.abs(data)**2))
            if self.noise_weights == "running-median":
                from gravfold.noise_weights import running_median_weights, PARAMETERS
                invvar, summary = running_median_weights(times, data, **PARAMETERS)
                summary["signal_to_noise_power_gain_over_file_variance"] = (
                    summary["signal_to_noise_power"] / file_variance_snr)
                self.noise_summary[det] = summary
            elif self.noise_weights != "file-variance":
                raise ValueError(f"unknown noise weighting: {self.noise_weights}")
            sim = HeterodynedCWSimulator(par=base, det=det, times=times)
            cache = PhaseCache(base, det, times[0], times[-1])
            self.chunks.append(dict(det=det, times=times, data=data, sim=sim, phase_cache=cache,
                                    invvar=invvar, gps=int(times[0]),
                                    n=data.size, n_sub=len(subs)))
        self.n_samp = int(sum(c["n"] for c in self.chunks))
        self._node_templates = None

    def _par_for(self, theta, freq_offset=0.0):
        """`freq_offset` is in GW Hz — the band the data, the filter passband
        and the instrumental lines all live in. The par carries SPIN
        frequencies and the model runs at freqfactor=2, so the spin shift is
        HALF the requested GW shift. Until 2026-07-14 this method added the
        offset to F0 directly, silently doubling every background offset: 4 of
        the 16 trials landed in the het filter's stopband (2F ~ 1e-4, dead
        correlators counted as noise draws), 6 aliased past the 1 Hz product's
        Nyquist back into the band, and the comb/guard exclusions were
        enforced at the wrong frequencies."""
        par = self._ensure_amp(self._apply_orbit(self._PP(self.search_par), theta))
        if freq_offset:
            F = np.array(par["F"], dtype=float)   # trap #3: F0 setter is a no-op
            F[0] += freq_offset / 2.0             # GW Hz -> spin Hz (freqfactor 2)
            par["F"] = F.tolist()
        return par

    def _model(self, chunk, par):
        carrier = 2*float(par["F0"]) - self.search_frequency
        m = chunk["phase_cache"].model(chunk["sim"], par, carrier_hz=carrier)
        if not np.any(m):
            raise RuntimeError("zero template — check par amplitude fields")
        return m

    def two_F(self, theta, freq_offset=0.0, data_override=None):
        par = self._par_for(theta, freq_offset)
        Y = 0.0 + 0.0j; NORM = 0.0
        for k, c in enumerate(self.chunks):
            m = self._model(c, par)
            d = c["data"] if data_override is None else data_override[k]
            iv = c["invvar"]                       # per-sample (studentized)
            Y += np.sum(iv * np.conj(m) * d)
            NORM += np.sum(iv * np.abs(m) ** 2)
        return 2.0 * (abs(Y) ** 2) / NORM if NORM > 0 else 0.0

    # ---- fast path: cache node templates once, reuse across data realizations
    def precompute_nodes(self, nodes):
        """Store per-node, per-chunk templates + norms so later recoveries are
        pure dot products. Memory: sum_nodes n_samp complex ~ manage node count."""
        tmpl = []
        for theta in nodes:
            par = self._par_for(theta)
            per = []
            norm = 0.0
            for c in self.chunks:
                m = self._model(c, par)
                per.append(m)
                norm += c["invvar"] * np.sum(np.abs(m) ** 2)
            tmpl.append(dict(m=per, norm=norm))
        self._node_templates = tmpl
        return self

    def two_F_nodes_cached(self, data_list):
        """2F at every cached node against the given per-chunk data list."""
        out = np.empty(len(self._node_templates))
        for j, t in enumerate(self._node_templates):
            Y = 0.0 + 0.0j
            for k, c in enumerate(self.chunks):
                Y += c["invvar"] * np.sum(np.conj(t["m"][k]) * data_list[k])
            out[j] = 2.0 * (abs(Y) ** 2) / t["norm"] if t["norm"] > 0 else 0.0
        return out


def _build_posterior_and_metric(demod, cov_scale=1.0):
    from gravfold.posterior import OrbitalPosterior
    from cwinpy.parfile import PulsarParameters
    cov, _ = load_covariance(COVARIANCE_ARTIFACT)
    post = OrbitalPosterior(demod.central_theta, cov * cov_scale, names=KEYS)
    samples = []
    for chunk in demod.chunks:
        indices = np.linspace(0, chunk["n"]-1, min(4096, chunk["n"]), dtype=int)
        samples.append((chunk["det"], chunk["times"][indices],
                        chunk["invvar"][indices] * chunk["n"] / len(indices)))
    metric = waveform_metric(demod.search_par, post.cov, samples)
    frequency = 2 * float(PulsarParameters(demod.search_par)["F0"])
    window = (min(c["times"][0] for c in demod.chunks),
              max(c["times"][-1] for c in demod.chunks))
    return post, metric, frequency, window


def _noise_record(demod):
    summary = getattr(demod, "noise_summary", None)
    mode = getattr(demod, "noise_weights", NOISE_WEIGHTS)
    record = dict(mode=mode, per_detector=summary or {})
    if mode == "running-median":
        from gravfold.noise_weights import PARAMETERS, equivalent_noise_bandwidth
        record.update(parameters=dict(PARAMETERS),
                      equivalent_noise_bandwidth_hz=equivalent_noise_bandwidth())
        for det, s in (summary or {}).items():
            print(f"[noise] {det}: kept {s['kept']:,} of {s['samples']:,} "
                  f"(stride ends {s['dropped_stride_ends']:,}, short segments "
                  f"{s['dropped_short_segments']:,}, above power cut "
                  f"{s['dropped_above_power_cut']:,}); noise power per unit "
                  f"weight {s['noise_power_per_unit_weight']:.3f}; "
                  f"signal-to-noise power x"
                  f"{s['signal_to_noise_power_gain_over_file_variance']:.1f} "
                  f"over file variance", flush=True)
    return record


def stage_measure(args, files):
    print(f"[measure] loading {len(files)} chunks ...", flush=True)
    t0 = time.time()
    demod = MultiDatasetDemodulator(files, BASE_PAR, studentize=True)
    noise_record = _noise_record(demod)
    print(f"[measure] loaded {demod.n_samp:,} samples across {len(demod.chunks)} "
          f"chunks ({time.time()-t0:.1f}s)", flush=True)
    nH = sum(c["n"] for c in demod.chunks if c["det"] == "H1")
    nL = sum(c["n"] for c in demod.chunks if c["det"] == "L1")
    print(f"[measure] H1 {nH:,} samp / L1 {nL:,} samp", flush=True)
    # one full 2F eval = the per-node cost
    tt = time.time(); F0 = demod.two_F(demod.central_theta)
    dt = time.time() - tt
    print(f"[measure] fixed-ephemeris 2F = {F0:.3f}  (per-node wall {dt:.2f}s)",
          flush=True)
    post, g, f_gw, win = _build_posterior_and_metric(demod)
    from gravfold.posterior import nodes_metric
    pts, w, info = nodes_metric(post, g, mu_max=args.mu, max_nodes=args.max_nodes)
    n = info["n_nodes"]
    print(f"[measure] window {win[0]:.0f}-{win[1]:.0f} "
          f"({(win[1]-win[0])/86400:.1f} d), f_GW={f_gw:.4f} Hz", flush=True)
    print(f"[measure] metric nodes = {n}  (mu_eff={info['effective_mu']:.3g}, "
          f"active dims={int(np.sum(np.asarray(info.get('n_per_dim',[]))>1))})",
          flush=True)
    print(f"[measure] EST on-source search = {n} nodes x {dt:.2f}s "
          f"= {n*dt/60:.1f} min", flush=True)
    print(f"[measure] EST background ({args.n_offsets} offsets) = "
          f"{args.n_offsets*n*dt/60:.1f} min", flush=True)
    print(f"[measure] EST injection (precompute {n} nodes once = {n*dt/60:.1f} "
          f"min, then {args.n_inj} inj x {len(demod.chunks)} dot ~ cheap)",
          flush=True)
    outp = args.out or os.path.join(OUTDIR, "ul_measure.json")
    payload = dict(n_chunks=len(demod.chunks), n_samp=demod.n_samp, nH=nH, nL=nL,
                   per_node_s=dt, n_nodes=int(n), mu_eff=info["effective_mu"],
                   fixed_2F=float(F0), f_gw=f_gw, window=win,
                   covariance_artifact=os.path.relpath(COVARIANCE_ARTIFACT, HERE),
                   covariance_sha256=covariance_sha256(COVARIANCE_ARTIFACT),
                   search_par_sha256=covariance_sha256(SEARCH_PAR),
                   reference_par_sha256=covariance_sha256(BASE_PAR),
                   input_identity=getattr(args, "input_identity", None),
                   noise_weights=noise_record,
                   metric=g.tolist(), nodes=np.asarray(pts).tolist())
    os.makedirs(os.path.dirname(os.path.abspath(outp)), exist_ok=True)
    tmp = outp + ".part"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, outp)
    print(f"[measure] wrote {outp}", flush=True)


# ---------------------------------------------------------------------------
# Parallel workers. Each worker builds its OWN demod after fork (sims built
# post-fork avoid LAL fork-safety issues), then answers model-heavy tasks.
# ---------------------------------------------------------------------------
_W = {}   # per-worker global state


def _winit(files, base_par, node_thetas):
    demod = MultiDatasetDemodulator(files, base_par, studentize=True)
    # full concatenated real data + per-sample inverse-variance, fixed chunk order
    d_full = np.concatenate([c["data"] for c in demod.chunks]).astype(np.complex128)
    iv_full = np.concatenate([c["invvar"] for c in demod.chunks])
    _W.update(demod=demod, d_full=d_full, iv_full=iv_full,
              slices=_chunk_slices(demod), node_thetas=np.asarray(node_thetas))


def _chunk_slices(demod):
    s, out = 0, []
    for c in demod.chunks:
        out.append((s, s + c["n"])); s += c["n"]
    return out


def _node_template_full(demod, theta, freq_offset=0.0):
    par = demod._par_for(theta, freq_offset)
    return np.concatenate([demod._model(c, par) for c in demod.chunks])


def _task_precompute(k):
    """Compute node k's full template; write to the shared memmap; return
    (NORM_k, Ydata_k)."""
    demod = _W["demod"]; iv = _W["iv_full"]; d = _W["d_full"]
    m = _node_template_full(demod, _W["node_thetas"][k]).astype(np.complex128)
    mm = np.memmap(_W["mm_path"], dtype=np.complex64, mode="r+",
                   shape=_W["mm_shape"])
    mm[k] = m.astype(np.complex64); mm.flush(); del mm
    NORM = float(np.sum(iv * np.abs(m) ** 2))
    Yd = complex(np.sum(iv * np.conj(m) * d))
    return k, NORM, Yd


def _task_background(args):
    """2F at node k, frequency offset df, against REAL data (noise-only)."""
    k, df = args
    demod = _W["demod"]; iv = _W["iv_full"]; d = _W["d_full"]
    m = _node_template_full(demod, _W["node_thetas"][k], df)
    Y, NORM, diagnostics = 0j, 0., {}
    for chunk, (lo, hi) in zip(demod.chunks, _W["slices"]):
        norm = float(np.sum(iv[lo:hi] * np.abs(m[lo:hi])**2))
        y = complex(np.sum(iv[lo:hi] * np.conj(m[lo:hi]) * d[lo:hi]))
        Y += y; NORM += norm
        diagnostics[chunk["det"]] = dict(twoF=2*abs(y)**2/norm,
                                        Y=[y.real, y.imag], norm=norm)
    return k, df, (2.0 * abs(Y)**2/NORM if NORM > 0 else 0.0), diagnostics


def _task_injection(spec):
    """Draw a random-orientation signal at a posterior-sampled orbit, unit H0;
    return the per-node cross terms Y_k^s = <m_node_k, s_unit>_w so the main
    process can scale by any h0. Reads cached node templates from the memmap."""
    j, theta_inj, cosi, psi, phi0 = spec
    demod = _W["demod"]; iv = _W["iv_full"]
    from cwinpy.parfile import PulsarParameters
    par = demod._apply_orbit(PulsarParameters(demod.search_par), theta_inj)
    par["H0"] = 1.0; par["COSIOTA"] = float(cosi)
    par["PSI"] = float(psi); par["PHI0"] = float(phi0)
    s = np.concatenate([demod._model(c, par) for c in demod.chunks])
    opt_snr2 = float(np.sum(iv * np.abs(s) ** 2))          # optimal SNR^2 at h0=1
    mm = np.memmap(_W["mm_path"], dtype=np.complex64, mode="r", shape=_W["mm_shape"])
    ivs = iv * s                                   # 427 MB, reused across nodes
    # vdot(a,b)=sum(conj(a)*b) is a fused reduction -> NO 427 MB temp per node
    # (the un-fused conj(mm[k])*ivs spiked workers to ~7 GB and OOM-thrashed).
    Ys = np.array([np.vdot(mm[k], ivs) for k in range(_W["mm_shape"][0])])
    del mm
    return j, Ys, opt_snr2


def _set_mm(mm_path, mm_shape):
    _W["mm_path"] = mm_path; _W["mm_shape"] = mm_shape


def _marginal_stats(twoF_nodes, logw):
    from scipy.special import logsumexp
    twoF_nodes = np.asarray(twoF_nodes)
    return dict(best=float(np.max(twoF_nodes)),
                logZ=float(logsumexp(0.5 * twoF_nodes + logw)))


def _astro(h0):
    d = DIST_PC * PC
    f_gw = 2.0 * F_SPIN
    eps = h0 * C**4 * d / (4 * np.pi**2 * G * I_ZZ * f_gw**2)
    h0_sd = (1.0 / d) * np.sqrt(5 * G * I_ZZ * abs(FDOT) / (2 * C**3 * F_SPIN))
    return dict(ellipticity=eps, h0_spindown=h0_sd, sd_ratio=h0 / h0_sd)


def stage_run(args, files):
    import multiprocessing as mp
    from scipy.special import logsumexp
    # 1) build once in main to place nodes + get sizes
    print(f"[run] loading {len(files)} chunks in main ...", flush=True)
    demod = MultiDatasetDemodulator(files, BASE_PAR, studentize=True)
    noise_record = _noise_record(demod)
    post, g, f_gw, win = _build_posterior_and_metric(demod)
    from gravfold.posterior import nodes_metric
    nodes, w, info = nodes_metric(post, g, mu_max=args.mu, max_nodes=args.max_nodes)
    logw = np.log(np.clip(np.asarray(w, float), 1e-300, None))
    n_nodes = len(nodes); n_samp = demod.n_samp
    print(f"[run] {n_nodes} nodes, {n_samp:,} samples, "
          f"{(win[1]-win[0])/86400:.1f} d, f_GW={f_gw:.4f}", flush=True)
    del demod  # workers rebuild their own

    outp = args.out or os.path.join(OUTDIR, "ul_result.json")
    os.makedirs(os.path.dirname(os.path.abspath(outp)), exist_ok=True)
    # Scope the large transient template cache to this result.  A shared fixed
    # path lets two smoke/science runs silently overwrite each other's memmap.
    mm_path = outp + ".node_templates.c64.mmap"

    ctx = mp.get_context("fork")
    nworkers = args.nproc if args.nproc else max(1, mp.cpu_count() - 2)
    print(f"[run] one pool, {nworkers} workers (each loads all chunks once) ...",
          flush=True)
    # THE FIXED-EPHEMERIS TEMPLATE IS NOT NODE 0. nodes_metric lays a tensor
    # grid of CELL CENTRES over the [-z, +z] box; node 0 is its most-negative
    # CORNER, and when a dimension has an even count the posterior mean is not
    # a grid point at all. Reading twoF[0] as "the fixed ephemeris" (as this
    # script did until 2026-07-13) reports a ~3-sigma corner template instead
    # of the central one — the two differ by an order of magnitude here
    # (1.25 vs 11.65). Evaluate the mean explicitly, as one extra template.
    nodes_ext = np.vstack([np.asarray(nodes, float),
                           np.asarray(post.mean, float)[None, :]])
    i_central = n_nodes                      # index of the central template
    n_ext = nodes_ext.shape[0]
    mm_shape = (n_ext, n_samp)               # the memmap holds it too
    np.memmap(mm_path, dtype=np.complex64, mode="w+", shape=mm_shape).flush()
    NORM = np.zeros(n_ext); Ydata = np.zeros(n_ext, complex)
    offs = _offsets(args.n_offsets)
    # background rows cover the BANK plus the central template: the bank rows
    # calibrate best-of-nodes; the central row gives the single-template
    # off-source distribution that the fixed-ephemeris 2F is judged against
    bg = {df: np.zeros(n_ext) for df in offs}
    bg_corrected = {}
    bg_detector = {}
    bg_response = np.array([])
    bg_stat_raw = np.array([])
    bg_central_raw = np.array([])
    p_fixed_raw = None
    rng = np.random.default_rng(20260707)
    pts_inj, _ = post.sample_mc(args.n_inj, rng)   # orbits ~ posterior
    specs = [(j, pts_inj[j], rng.uniform(-1, 1), rng.uniform(0, np.pi/2),
              rng.uniform(0, 2*np.pi)) for j in range(args.n_inj)]
    capsule = None
    if getattr(args, "input_identity", None) is not None:
        from gravfold.j0437_campaign import freeze_run
        plan = dict(input_identity=args.input_identity,
                    search_par_sha256=covariance_sha256(SEARCH_PAR),
                    reference_par_sha256=covariance_sha256(BASE_PAR),
                    covariance_sha256=covariance_sha256(COVARIANCE_ARTIFACT),
                    parameter_order=list(KEYS), mean=post.mean.tolist(),
                    nodes=np.asarray(nodes).tolist(), weights=np.asarray(w).tolist(),
                    offsets=offs, orbital_draws=pts_inj.tolist(),
                    orientation_draws=[list(spec[2:]) for spec in specs],
                    population_seed=20260707, metric=g.tolist(), window=list(win),
                    args=vars(args), phase_cache_step_seconds=300,
                    phase_inputs={os.path.basename(str(path)): dict(path=str(path), sha256=covariance_sha256(path))
                                  for path in (SEARCH_PAR, BASE_PAR, COVARIANCE_ARTIFACT)})
        capsule = freeze_run(outp, plan, COVARIANCE_ARTIFACT, HERE)
    Ys = np.zeros((args.n_inj, n_ext), complex); snr2 = np.zeros(args.n_inj)
    t0 = time.time()
    with ctx.Pool(nworkers, initializer=_winit_full,
                  initargs=(files, BASE_PAR, nodes_ext, mm_path, mm_shape)) as p:
        # 1) precompute node templates (+ the central one) -> memmap
        print(f"[run] precompute {n_nodes} node templates + the central "
              f"ephemeris -> memmap ...", flush=True)
        for k, nrm, yd in p.imap_unordered(_task_precompute, range(n_ext)):
            NORM[k] = nrm; Ydata[k] = yd
        twoF_ext = 2.0 * np.abs(Ydata) ** 2 / np.where(NORM > 0, NORM, np.inf)
        twoF_on = twoF_ext[:n_nodes]         # the BANK — marginalise over these
        twoF_central = float(twoF_ext[i_central])   # the FIXED ephemeris
        on = _marginal_stats(twoF_on, logw)
        print(f"[run] ON-SOURCE: best-node 2F={on['best']:.3f}, "
              f"logZ={on['logZ']:.3f}, fixed-ephemeris (posterior MEAN) "
              f"2F={twoF_central:.3f}  [node0 diagnostic={twoF_on[0]:.3f}]  "
              f"({time.time()-t0:.0f}s)", flush=True)

        # 2) background: off-source frequency offsets (guard band + H1 combs skipped).
        #    --n-offsets 0 skips this; it never imports a prior background.
        bg_central = np.array([])
        if offs:
            print(f"[run] background: {len(offs)} offsets x {n_ext} templates "
                  f"(bank + central) ...", flush=True)
            for k, df, val, diagnostics in p.imap_unordered(
                    _task_background,
                    [(k, df) for df in offs for k in range(n_ext)]):
                bg[df][k] = val
                bg_detector.setdefault(f"{df:+.4f}", {})[str(k)] = diagnostics
            bg_stat_raw = np.array([
                _marginal_stats(bg[df][:n_nodes], logw)["best"] for df in offs
            ])
            bg_central_raw = np.array([bg[df][i_central] for df in offs])
            # Controls farther from the heterodyne center are deterministically
            # attenuated by both the forward/backward Butterworth and the 1-s
            # block average.  Normalize them to the target-frequency response
            # before treating the controls as exchangeable background draws.
            bg_response = heterodyne_statistic_response(offs)
            bg_corrected = {
                df: bg[df] / bg_response[index]
                for index, df in enumerate(offs)
            }
            bg_stat = bg_stat_raw / bg_response
            bg_central = bg_central_raw / bg_response
            # finite-trials estimator: (k+1)/(n+1), not k/n — a zero from a
            # finite background is never evidence of p = 0
            k_exc = int(np.sum(bg_stat >= on["best"]))
            p_value = float((k_exc + 1) / (len(bg_stat) + 1))
            k_fix = int(np.sum(bg_central >= twoF_central))
            p_fixed = float((k_fix + 1) / (len(bg_central) + 1))
            k_fix_raw = int(np.sum(bg_central_raw >= twoF_central))
            p_fixed_raw = float((k_fix_raw + 1) / (len(bg_central_raw) + 1))
            print(f"[run] background best-of-nodes 2F: "
                  f"median={np.median(bg_stat):.2f}, "
                  f"p95={np.percentile(bg_stat,95):.2f}, max={bg_stat.max():.2f}; "
                  f"on-source p = ({k_exc}+1)/({len(bg_stat)}+1) = {p_value:.3f}",
                  flush=True)
            print(f"[run] central-template background 2F: "
                  f"median={np.median(bg_central):.2f}, "
                  f"max={bg_central.max():.2f}; fixed-ephemeris "
                  f"2F={twoF_central:.2f} -> p = ({k_fix}+1)/"
                  f"({len(bg_central)}+1) = {p_fixed:.3f}", flush=True)
        else:
            # No implicit reuse: historical results may have a different bank,
            # central ephemeris, or data selection, even at the same output path.
            bg_stat = np.array([]); p_value = None; p_fixed = None; offs = []
            print("[run] background SKIPPED (--n-offsets 0); p-value will be null",
                  flush=True)

        # 3) injections -> orientation-averaged h0_95%
        print(f"[run] injections: {args.n_inj} random-orientation signals ...",
              flush=True)
        for j, ys, s2 in p.imap_unordered(_task_injection, specs):
            Ys[j] = ys; snr2[j] = s2
    # Injection recovery. Per injection j, the best-of-BANK statistic against
    # (real data + h0*signal) is quadratic in h0:
    #   2F_k(h0) = 2|Ydata_k + h0 Ys_k|^2 / NORM_k,  detect = max_k 2F_k >= 2F*
    # over the n_nodes bank templates ONLY (row i_central is the fixed-ephemeris
    # diagnostic, not part of the detection statistic). Two thresholds:
    #   loudest_event      : 2F* = on-source loudest. Standard but conflated
    #                        -- eff(h0->0) ~ 0.5 because injections share the
    #                        on-source noise realization (its own loud node).
    #   loudest_background : 2F* = the loudest off-source best-of-bank value.
    #                        An EMPIRICAL threshold from n finite trials — it
    #                        was previously labelled "fixed_fap", which
    #                        overstated it: n trials bound the FAP only as
    #                        ~<1/(n+1), they do not fix it.
    # Ydata/Ys/NORM/snr2 are SAVED so any 2F* is re-derivable with no re-run.
    Yd_bank = Ydata[:n_nodes]
    NORM_bank = np.where(NORM[:n_nodes] > 0, NORM[:n_nodes], np.inf)
    thr_loud = on["best"]
    thr_bg = args.detect_threshold
    if thr_bg is None and len(bg_stat):
        thr_bg = float(bg_stat.max())
        print(f"[run] loudest-background threshold = background max = "
              f"{thr_bg:.2f} (over {len(bg_stat)} trials)", flush=True)
    h0_grid = _h0_grid(args, max(thr_loud, thr_bg or 0.0), NORM[:n_nodes], snr2)

    def detections(twostar):
        detected = np.empty((len(h0_grid), args.n_inj), dtype=bool)
        for i, h0 in enumerate(h0_grid):
            for j in range(args.n_inj):
                twoF = 2.0 * np.abs(Yd_bank + h0 * Ys[j, :n_nodes]) ** 2 / NORM_bank
                detected[i, j] = np.max(twoF) >= twostar
        return detected

    uls = {}
    detected_loud = detections(thr_loud)
    eff_loud = detected_loud.mean(axis=1)
    uls["loudest_event"] = dict(threshold=float(thr_loud),
                                h0_95=_interp_crossing(h0_grid, eff_loud, 0.95),
                                efficiency=eff_loud.tolist())
    if thr_bg:
        detected_bg = detections(thr_bg)
        eff_bg = detected_bg.mean(axis=1)
        uls["loudest_background"] = dict(
            threshold=float(thr_bg),
            h0_95=_interp_crossing(h0_grid, eff_bg, 0.95),
            efficiency=eff_bg.tolist())
    # primary UL = loudest_background when a threshold exists, else loudest_event
    primary = "loudest_background" if thr_bg else "loudest_event"
    h0_95 = uls[primary]["h0_95"]
    primary_detected = detected_bg if thr_bg else detected_loud
    bootstrap = paired_grid_bootstrap_crossing(
        h0_grid, primary_detected, target=0.95, confidence=0.90,
        n_bootstrap=20000, seed=20260825)
    astro = _astro(h0_95) if h0_95 else None
    for name, u in uls.items():
        tag = " (PRIMARY)" if name == primary else ""
        h = u["h0_95"]
        print(f"[run] UL[{name}]{tag}: 2F*={u['threshold']:.2f} -> "
              + (f"h0_95={h:.3e}" if h else "not bracketed"), flush=True)
    if astro:
        print(f"[run]   ellipticity_95 = {astro['ellipticity']:.3e}", flush=True)
        print(f"[run]   spin-down h0    = {astro['h0_spindown']:.3e}  "
              f"(sd_ratio = {astro['sd_ratio']:.3f}"
              f"{'  <-- BEATS spin-down!' if astro['sd_ratio']<1 else ''})",
              flush=True)

    out = dict(n_nodes=n_nodes, n_samp=n_samp, f_gw=f_gw, window=win,
               run_capsule=capsule,
               covariance_artifact=os.path.relpath(COVARIANCE_ARTIFACT, HERE),
               covariance_sha256=covariance_sha256(COVARIANCE_ARTIFACT),
               # self-document the run: the 2026-07-07 result was produced by
               # an invocation whose flags were never recorded — never again
               invocation=dict(argv=sys.argv[1:], args=vars(args),
                               timestamp=time.strftime("%Y-%m-%dT%H:%M:%S")),
               mu_eff=info["effective_mu"],
               search_par_sha256=covariance_sha256(SEARCH_PAR),
               reference_par_sha256=covariance_sha256(BASE_PAR),
               input_identity=getattr(args, "input_identity", None),
               noise_weights=noise_record,
               bank=dict(parameter_order=list(KEYS), mean=post.mean.tolist(),
                         nodes=np.asarray(nodes).tolist(), weights=np.asarray(w).tolist(),
                         metric=g.tolist()),
               phase_evaluation=dict(method="LAL T2 with 300-second residual-phase spline; exact frequency carrier"),
               # the fixed-ephemeris statistic, evaluated AT THE POSTERIOR
               # MEAN. `twoF_node0` is kept only to audit the pre-2026-07-13
               # artifacts, which recorded that box-corner value under the
               # `twoF_fixed` name (see ul_result.json vs ul_measure.json).
               twoF_fixed=twoF_central, twoF_node0=float(twoF_on[0]),
               twoF_fixed_p_empirical=p_fixed,
               twoF_fixed_p_empirical_raw=p_fixed_raw,
               onsource=on,
               background=dict(
                   per_detector=bg_detector,
                   convention="gw_hz",       # offsets ARE GW-frequency shifts
                   statistic_scale="target_frequency_response_corrected",
                   passband_hz=_OFF_PASSBAND_HZ, guard_hz=_OFF_GUARD_HZ,
                   comb_exclusions_hz=list(_OFF_COMBS_HZ),
                   offsets=list(offs),
                   response_model=dict(
                       kind=("forward-backward Butterworth times rectangular "
                             "block average"),
                       butterworth_order=9, butterworth_cutoff_hz=0.25,
                       averaging_seconds=1.0,
                       statistic_response_formula=(
                           "[1+(abs(df)/0.25)^18]^-2 * sinc(df*1s)^2"),
                       response_by_offset=bg_response.tolist()),
                   best_node_2F_raw=bg_stat_raw.tolist(),
                   central_2F_raw=bg_central_raw.tolist(),
                   per_node_2F_raw={f"{df:+.4f}": bg[df].tolist() for df in offs},
                   best_node_2F=bg_stat.tolist(),
                   central_2F=bg_central.tolist(),
                   per_node_2F={f"{df:+.4f}": bg_corrected[df].tolist()
                                for df in offs},
                   n_trials=len(bg_stat),
                   n_exceed_raw=(int(np.sum(bg_stat_raw >= on["best"]))
                                 if len(bg_stat_raw) else None),
                   p_value_raw=(float((np.sum(bg_stat_raw >= on["best"]) + 1)
                                      / (len(bg_stat_raw) + 1))
                                if len(bg_stat_raw) else None),
                   n_exceed=int(np.sum(bg_stat >= on["best"])) if len(bg_stat) else None,
                   p_value=p_value,          # (k+1)/(n+1) finite-trials estimator
                   p_estimator="(k+1)/(n+1)"),
               injection=dict(n_inj=args.n_inj, h0_grid=h0_grid.tolist(),
               population_seed=20260707, orbital_draws=pts_inj.tolist(),
               orientation_draws=[list(spec[2:]) for spec in specs],
               uls=uls, primary=primary, h0_95=h0_95,
               sampling_uncertainty=dict(
                   method="paired nonparametric bootstrap of injection identities",
                   estimator="first linear crossing on the stored amplitude grid",
                   h0_95_interval=bootstrap["interval"],
                   confidence=bootstrap["confidence"],
                   n_bootstrap=bootstrap["n_bootstrap"], seed=bootstrap["seed"],
                   valid_crossings=bootstrap["valid_crossings"],
                   valid_fraction=bootstrap["valid_fraction"],
                   scope=("Monte Carlo source-population sampling conditional "
                          "on the fixed on-source noise realization and observed "
                          "threshold; excludes strain calibration and repeated-"
                          "background uncertainty")),
               median_opt_snr2_at_h0_1=float(np.median(snr2)),
               # raw products: re-derive efficiency at ANY 2F* offline.
               # Rows 0..n_nodes-1 are the bank; row n_nodes is the central
               # (fixed-ephemeris) template, EXCLUDED from detection maxima.
               Ydata=[[float(z.real), float(z.imag)] for z in Ydata],
               Ys=[[[float(z.real), float(z.imag)] for z in row] for row in Ys],
               NORM=NORM.tolist(), snr2=snr2.tolist()),
               astro=astro)
    tmp = outp + ".part"
    with open(tmp, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    os.replace(tmp, outp)                      # never a half-written "result"
    print(f"[run] wrote {outp}", flush=True)
    try:
        os.remove(mm_path)
    except OSError:
        pass


def _winit_full(files, base_par, node_thetas, mm_path, mm_shape):
    _winit(files, base_par, node_thetas)
    _set_mm(mm_path, mm_shape)


# Where an off-source trial may live, in GW Hz relative to f_gw. The campaign
# products are 1 Hz complex series low-passed at a 0.25 Hz knee (9th-order
# Butterworth, applied forward+backward) and block-averaged for one second.
# The raw statistic response is the filter response times sinc(df)^2: at
# |df|=0.20 it is 0.844 (15.6% suppression), not the filter-only 0.965.
# Every completed control is divided by this known response before comparison
# with the target, while the uncorrected values remain in the result artifact.
# Thus:
#   * |df| <= 0.20 avoids the steep filter stopband and limits the size of the
#     deterministic correction (old doubled offsets read ~1e-4 there);
#   * |df| >= 0.03  clears the on-source guard band (orbital Doppler smear
#     ~ +-2(v/c)f = +-0.030 Hz for J0437);
#   * +-0.03 around each instrumental comb line (H1 heater combs at 347.175,
#     347.006, 348.114 Hz -> -0.201, -0.370, +0.738 in GW-offset space; only
#     the first is inside the passband).
_OFF_PASSBAND_HZ = 0.20
_OFF_GUARD_HZ = 0.03
_SEARCH_F0 = F_SPIN
_OFF_COMBS_HZ = tuple(f - 2*_SEARCH_F0 for f in (347.175, 347.006, 348.114))


def _offsets(n):
    """`n` symmetric off-source frequency offsets in GW Hz, all inside the
    filter passband, outside the on-source guard band, and clear of the H1
    comb lines. Deterministic: a dense +/- ladder over the allowed range,
    thinned evenly to n."""
    m = max(4 * n, 64)
    half = np.linspace(_OFF_GUARD_HZ + 0.002, _OFF_PASSBAND_HZ - 0.002, m)
    cand = np.round(np.concatenate([half, -half]), 4)
    cand = np.array([x for x in cand
                     if all(abs(x - c) > 0.03 for c in _OFF_COMBS_HZ)])
    cand = np.sort(np.round(cand, 4))
    keep = np.unique(np.round(np.linspace(0, cand.size - 1, n)).astype(int))
    return [float(x) for x in cand[keep]]


def _h0_grid(args, on_best, NORM, snr2):
    if args.h0:
        return np.array(args.h0, dtype=float)
    # centre the grid near the h0 whose median optimal SNR^2 ~ a few x on_best
    s2 = np.median(snr2) if np.median(snr2) > 0 else 1.0
    h0_star = np.sqrt(max(on_best, 4.0) / s2)
    return np.geomspace(0.15 * h0_star, 6.0 * h0_star, 16)


def _interp_crossing(x, y, target):
    return linear_grid_crossing(x, y, target)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["measure", "onsource", "upperlimit"],
                    default="measure")
    ap.add_argument("--mu", type=float, default=0.2)
    ap.add_argument("--max-nodes", type=int, default=20000)
    ap.add_argument("--n-offsets", type=int, default=16,
                    help="fresh control frequencies; 0 skips background without reuse")
    ap.add_argument("--n-inj", type=int, default=200)
    ap.add_argument("--detect-threshold", type=float, default=None,
                    help="fixed target-scale 2F* detection threshold "
                         "(e.g. the loudest off-source background). If set, this "
                         "UL is primary; else loudest-event is primary.")
    ap.add_argument("--nproc", type=int, default=1, help="worker count; default 1 to bound memory")
    ap.add_argument("--out", default=None,
                    help="result json path (default OUTDIR/ul_result.json); "
                         "point smoke runs elsewhere so they cannot clobber "
                         "the science artifact")
    ap.add_argument("--h0", type=float, nargs="*", default=None,
                    help="explicit h0 grid; default auto")
    ap.add_argument("--detectors", nargs="+", default=["H1", "L1"])
    ap.add_argument("--noise-weights",
                    choices=["file-variance", "running-median"],
                    default="file-variance",
                    help="file-variance uses one plain variance per file; "
                         "running-median drops the transient before each "
                         "stride end and samples above a power cut, and uses "
                         "a running noise estimate (gravfold/noise_weights.py)")
    ap.add_argument("--limit-chunks", type=int, default=0,
                    help="use only the first N chunks (fast dry run)")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(HETDIR, "het_*_g*.hdf5")))
    files = [f for f in files
             if re.search(r"_(%s)_2_" % "|".join(args.detectors), f)]
    if args.limit_chunks:
        # Exercise both detectors in a small smoke run rather than taking only
        # H1 from a lexicographically sorted campaign listing.
        import itertools
        groups = [[f for f in files if f"_{det}_" in os.path.basename(f)]
                  for det in args.detectors]
        interleaved = [f for row in itertools.zip_longest(*groups) for f in row if f]
        files = sorted(interleaved[:args.limit_chunks])
    if not files:
        sys.exit("no chunk files found in " + HETDIR)

    args.input_identity = verify_selection(files, BASE_PAR, smoke=bool(args.limit_chunks))
    global NOISE_WEIGHTS
    NOISE_WEIGHTS = args.noise_weights

    if args.stage == "measure":
        stage_measure(args, files)
    else:
        stage_run(args, files)


if __name__ == "__main__":
    main()
