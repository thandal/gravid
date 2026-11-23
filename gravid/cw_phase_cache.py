"""Interpolate the smooth heterodyne phase after removing the known carrier.

The frequency carrier is evaluated exactly at each data timestamp. Only the
slow binary/SSB residual is interpolated; validation must compare this path
with direct LAL phase evaluation over the full observing span.
"""
import numpy as np
from scipy.interpolate import CubicSpline


class PhaseCache:
    def __init__(self, reference_par, detector, start, end, step=300.0):
        from cwinpy.signal import HeterodynedCWSimulator
        self.times = np.arange(float(start)-step, float(end)+2*step, step)
        self.sim = HeterodynedCWSimulator(reference_par, detector, self.times)
        self.epoch = float(start)

    def interpolant(self, par, carrier_hz):
        phase = self.sim.model(newpar=par, updateSSB=True, updateBSB=True,
                               freqfactor=2.0, phase_only=True)
        carrier = 2*np.pi*carrier_hz*(self.times-self.epoch)
        residual = np.unwrap(np.angle(np.exp(-2j*np.pi*phase - 1j*carrier)))
        # Guard against inadequate residual sampling. A modest threshold leaves
        # ample headroom before phase unwrapping becomes ambiguous.
        if np.max(abs(np.diff(residual))) > 0.25:
            raise ValueError("phase-cache residual changes too rapidly; reduce step")
        return CubicSpline(self.times, residual, extrapolate=False)

    def phasor(self, par, times, carrier_hz=0.0):
        spline = self.interpolant(par, carrier_hz)
        phase = spline(times) + 2*np.pi*carrier_hz*(np.asarray(times)-self.epoch)
        if not np.all(np.isfinite(phase)):
            raise ValueError("phase cache does not cover requested timestamps")
        return np.exp(1j*phase)

    def model(self, direct_sim, par, carrier_hz=0.0):
        import lalpulsar
        amplitude = lalpulsar.HeterodynedPulsarGetAmplitudeModel(
            par.PulsarParameters(), 2.0, 1, 0, 0, direct_sim.gpstimes, direct_sim.resp)
        output = np.array(amplitude.data.data, dtype=complex, copy=True)
        del amplitude
        spline = self.interpolant(par, carrier_hz)
        for start in range(0, len(output), 1_000_000):
            sl = slice(start, start+1_000_000)
            # Campaign timestamps are exactly representable half seconds.
            # CWInPy stores them as long doubles; carrying that dtype into the
            # carrier exponential needlessly invokes slow complex256 trig.
            times = np.asarray(direct_sim.times[sl], dtype=np.float64)
            phase = spline(times) + 2*np.pi*carrier_hz*(times-self.epoch)
            if not np.all(np.isfinite(phase)):
                raise ValueError("phase cache does not cover requested timestamps")
            output[sl] *= np.exp(1j*phase)
        return output
