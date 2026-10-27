"""Checks of the response matrix (analysis/response_matrix.npz) against the stored paper results.
Only the .npz file, the PSTAR tables and the CSV summaries are needed (no raw runs, no GATE)."""
import csv
import numpy as np
import pytest
import response_matrix as rm
from spectra import SPECTRA

AREALS = [5, 10, 20, 30, 40]


@pytest.fixture(scope="module")
def npz():
    return dict(np.load(rm.OUT))


def test_matrix_contents(npz):
    assert npz["edges"].shape == (26,)
    assert npz["edges"][0] == pytest.approx(30.0) and npz["edges"][-1] == pytest.approx(5000.0)
    assert npz["z"].shape == (300,)
    for m in rm.MATS:
        for t in AREALS:
            for q in "DH":
                R = npz[f"{q}_{m}_{t}"]
                assert R.shape == (25, 300), (q, m, t)
                assert np.all(np.isfinite(R)) and np.all(R >= 0)


def test_dose_equivalent_not_below_dose(npz):
    # Q >= 1, so H >= D in every voxel that has dose (small tolerance for empty voxels)
    for m in rm.MATS:
        for t in AREALS:
            D, H = npz[f"D_{m}_{t}"], npz[f"H_{m}_{t}"]
            assert np.all(H >= D * (1 - 1e-9))


def test_fold_reference_value():
    # Koldobskiy 1989, 20 g/cm2 Al, mean dose over 0-30 g/cm2 of water (paper, Sect. 3.7)
    assert rm.fold(SPECTRA["kold1989"], "Al", 20, "D") == pytest.approx(12.58, rel=2e-3)


def _validation_rows(window):
    with open("analysis/response_matrix_validation.csv") as f:
        return [r for r in csv.DictReader(f) if r["window"] == window]


@pytest.mark.parametrize("spectrum", ["kold1956", "kold1972", "kold1989", "webber1956", "king1972"])
def test_fold_reproduces_stored_validation(spectrum):
    # Re-folding gives the stored values, and the folded mean D agrees with the direct GATE runs within 1.5 %
    rows = [r for r in _validation_rows("mean") if r["spectrum"] == spectrum]
    assert len(rows) == 2 * 7 * 3
    for r in rows:
        if r["quantity"] != "D":
            continue
        p = rm.fold(SPECTRA[spectrum], r["material"], int(r["areal"]), "D", "mean")
        assert p == pytest.approx(float(r["folded"]), rel=1e-6)
        assert abs(p / float(r["GATE"]) - 1) <= 0.015, (r["material"], r["areal"])


def test_mean_dose_equivalent_within_paper_bounds():
    # Paper: mean H within +-1.5 % for four spectra and down to -3 % for King 1972 behind 10-20 g/cm2
    for r in _validation_rows("mean"):
        if r["quantity"] != "H":
            continue
        ratio = float(r["ratio"])
        lo = 0.968 if r["spectrum"] == "king1972" else 0.985
        assert lo <= ratio <= 1.015, (r["spectrum"], r["material"], r["areal"], ratio)


def test_nasa_limit_table():
    # Table 5: minimum areal density for BFO < 250 mGy-Eq (RBE 1.5), recomputed from the matrix
    import nasa_limits_thick as nl
    with open("analysis/nasa_limits_thick.csv") as f:
        stored = {(r["spectrum"], r["material"]): r for r in csv.DictReader(f)}
    for sp, m, expected in (("kold1989", "Al", "23.8"), ("kold1989", "POLYETHYLENE", "17.6"), ("kold1956", "Al", "19.8"),
                            ("webber1956", "POLYETHYLENE", "<=5")):
        bfo = [1.5 * 10 * rm.fold(SPECTRA[sp], m, t, "D", "BFO") for t in AREALS]
        assert nl.tmin(bfo) == expected == stored[(sp, m)]["min_areal_g_cm2"]
        assert np.allclose(bfo, [float(stored[(sp, m)][f"BFO_{t}_mGyEq"]) for t in AREALS], rtol=2e-3, atol=0.06)
