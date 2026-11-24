"""Compare two ways to reduce six minutes of synthetic J0437 detector strain."""
import sys
from pathlib import Path

import numpy as np
from cwinpy.parfile import PulsarParameters
from cwinpy.signal import HeterodynedCWSimulator
from scipy.signal import butter, sosfiltfilt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from gravid.cw_phase_cache import PhaseCache
from gravid.j0437_model import SEARCH_PAR


def main():
    rate = 4096
    times = 1378306153.0 + np.arange(rate * 360) / rate
    seconds = times.reshape(-1, rate).mean(axis=1)
    lowpass = butter(9, 0.25, fs=rate, output="sos")
    reference = PulsarParameters(str(ROOT / "results/campaign_het_v2/base_cwinpy.par"))
    search = PulsarParameters(str(SEARCH_PAR))
    for par in (reference, search):
        par["H0"], par["COSIOTA"], par["PSI"], par["PHI0"] = 1.0, 0.3, 0.4, 0.7

    def reduce_strain(strain, carrier):
        filtered = sosfiltfilt(lowpass, strain * np.conj(carrier))
        return filtered.reshape(-1, rate).mean(axis=1)

    for detector in ("H1", "L1"):
        sim = HeterodynedCWSimulator(reference, detector, times)
        phase = sim.model(phase_only=True, freqfactor=2.0)
        carrier = np.exp(-2j * np.pi * phase)
        signal = sim.model(newpar=search, updateSSB=True, updateBSB=True)
        strain = 2 * np.real(signal * carrier)
        phase_shift = sim.model(newpar=search, updateSSB=True, updateBSB=True, phase_only=True)
        direct = reduce_strain(strain, carrier * np.exp(-2j * np.pi * phase_shift))
        cache = PhaseCache(reference, detector, seconds[0], seconds[-1])
        rephased = reduce_strain(strain, carrier) * np.conj(cache.phasor(search, seconds))

        # Exclude the edges of this short synthetic time series.
        direct, rephased = direct[90:-90], rephased[90:-90]
        relative_error = np.linalg.norm(rephased - direct) / np.linalg.norm(direct)
        print(f"{detector} relative difference {relative_error:.3g}", flush=True)
        assert relative_error < 1e-3


if __name__ == "__main__":
    main()
