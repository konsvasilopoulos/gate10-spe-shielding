"""Checks of the SPE spectra (spectra.py, analysis/gle_table.py). No GATE needed."""
import numpy as np
import pytest
import spectra
from gle_table import T2, F as gle_F

trapz = getattr(np, "trapezoid", None) or np.trapz


def test_webber1956_normalisation():
    # Phi(>E) = 1e9 exp[(239.1 - P)/100], P = rigidity in MV: Phi(>30 MeV) ~ 1e9 cm^-2 (NASA/TP-2009-215560)
    assert spectra.webber1956_integral(30.0) == pytest.approx(9.994e8, rel=1e-3)


def test_king1972_normalisation():
    # King (1974): Phi(>30 MeV) = 7.9e9 cm^-2
    assert spectra.king1972_integral(30.0) == pytest.approx(7.9e9, rel=1e-12)


def test_webber_differential_integrates_to_integral():
    E = np.geomspace(10, 2000, 20001)
    num = trapz(spectra.webber1956_differential(E), E)
    exact = spectra.webber1956_integral(10) - spectra.webber1956_integral(2000)
    assert num == pytest.approx(exact, rel=1e-4)


@pytest.mark.parametrize("name", ["webber1956", "king1972", "kold1956", "kold1972", "kold1989", "logflat"])
def test_spectrum_bins_conserve_fluence(name):
    emin, emax = (30.0, 5000.0) if name.startswith("kold") or name == "logflat" else (10.0, 2000.0)
    edges, w, tot = spectra.spectrum_bins(name, emin, emax, 100)
    F = spectra.SPECTRA[name]
    assert np.all(w > 0)
    assert tot == pytest.approx(F(emin) - F(emax), rel=1e-12)
    assert len(edges) == 101 and edges[0] == pytest.approx(emin) and edges[-1] == pytest.approx(emax)


@pytest.mark.parametrize("gle", sorted(T2))
def test_koldobskiy_band_is_continuous_and_decreasing(gle):
    # The modified Band function is continuous at the break rigidity and Phi(>E) decreases with E
    E = np.geomspace(30, 5000, 2000)
    f = gle_F(E, gle)
    assert np.all(np.isfinite(f)) and np.all(f > 0)
    assert np.all(np.diff(f) < 0)
    assert np.max(np.abs(np.diff(np.log(f)))) < 0.2      # no jump at the break


def test_spectra_module_matches_gle_table():
    # spectra.KOLDOBSKIY_GLE and gle_table.T2 must describe the same events
    E = np.array([50.0, 150.0, 500.0, 2000.0])
    for g in spectra.KOLDOBSKIY_GLE:
        assert np.allclose(spectra.koldobskiy_gle_integral(E, g), gle_F(E, g), rtol=1e-6)
