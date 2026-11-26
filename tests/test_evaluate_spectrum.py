"""evaluate_spectrum.py: the user-facing folding tool reproduces the article's numbers for built-in spectra."""
import os, sys
import numpy as np
sys.path.insert(0, "analysis")
import evaluate_spectrum as ev
from spectra import SPECTRA


def _row(rows, m, t):
    return next(r for r in rows if r["material"] == m and r["areal_g_cm2"] == t)


def test_kold1972_ratios_and_limits():
    rows, lim = ev.evaluate(SPECTRA["kold1972"], ["POLYETHYLENE", "Pb"])
    assert abs(_row(rows, "POLYETHYLENE", 20)["ratio_mean_D"] - 0.580) < 0.003        # article, Table 3 (folded)
    assert lim["Al"] == "16.8" and lim["POLYETHYLENE"] == "12.2" and lim["Pb"] == "26.4"   # article, Table 4
    assert abs(_row(rows, "POLYETHYLENE", 20)["gamma"] - 3.12) < 0.02


def test_table_input_matches_function():
    F = SPECTRA["kold1972"]; E = np.geomspace(10, 10000, 150)
    Ft = ev.tabulated(E, F(E))
    a, _ = ev.evaluate(F, ["POLYETHYLENE"], [20]); b, _ = ev.evaluate(Ft, ["POLYETHYLENE"], [20])
    assert abs(_row(b, "POLYETHYLENE", 20)["BFO_D_mGy"] / _row(a, "POLYETHYLENE", 20)["BFO_D_mGy"] - 1) < 0.005


def test_band_input_matches_td_analysis():
    from gle_rule_td import band
    F1 = ev.band(1.75e8, 1.76, 5.04, 0.566)
    E = np.array([30.0, 150.0, 1000.0])
    assert np.allclose(F1(E), band(E, 1.75e8, 1.76, 5.04, 0.566))
