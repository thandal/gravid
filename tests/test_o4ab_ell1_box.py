from pathlib import Path

import numpy as np
import pytest


@pytest.mark.parametrize("uncertain", ["A1", "TASC"])
def test_circular_orbit_delay_bound_keeps_epoch_and_axis_units(
    uncertain, tmp_path, monkeypatch
):
    scripts = Path(__file__).resolve().parents[1] / "scripts/timing/o4ab_targets"
    monkeypatch.syspath_prepend(str(scripts))
    from historical import ell1_box

    period_days = 9.5
    axis_seconds = 8.4
    sigma = 1e-6
    values = {
        "PB": period_days,
        "A1": axis_seconds,
        "TASC": 56000,
        "EPS1": 0,
        "EPS2": 0,
    }
    model = tmp_path / "circular.par"
    model.write_text(
        "F0 300\n"
        + "".join(
            f"{name} {value} 0 {sigma if name == uncertain else 0}\n"
            for name, value in values.items()
        )
    )
    # A change in x shifts the delay amplitude directly. A change in Tasc
    # moves the orbital longitude by 2*pi*dTasc/Pb, with both times in days.
    expected = (
        3 * sigma
        if uncertain == "A1"
        else axis_seconds * 2 * np.pi * 3 * sigma / period_days
    )
    assert ell1_box(model, 3)["delay_seconds"] == pytest.approx(expected, rel=1e-12)
