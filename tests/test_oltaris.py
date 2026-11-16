"""Checks of the GATE / HZETRN (OLTARIS) comparison against the stored tables (no GATE, no OLTARIS needed)."""
import csv
import numpy as np
import pytest
import oltaris_compare as oc


def test_oltaris_tables_are_read():
    depth, D = oc.read_table("oltaris/Al_King72/li_dose_database_mGy.dat")
    assert depth.shape == D.shape == (8,)
    assert depth[3] == pytest.approx(7.4101)                    # 20 g/cm2 Al = 7.4101 cm
    E, phi, parts = oc.read_table("oltaris/Al_Tylka89/li_flux.dat")
    assert "proton" in parts and phi.shape == (len(parts), len(E))


def test_band_function_rescaling():
    # The surface dose depends only on the input spectrum: after rescaling it must be equal for all materials
    for tag in ("GLE24", "GLE25"):
        s = [oc.oltaris(m, tag)[0][0] for m in ("Al", "PE", "Pb")]
        assert np.allclose(s, s[0], rtol=1e-4)


@pytest.mark.parametrize("spectrum", ["Webber56", "Kold72", "Tylka89"])
def test_gate_matches_hzetrn_behind_aluminium(spectrum):
    # Absorbed dose behind 20 g/cm2 Al at the BFO depth agrees within 3 %
    rows = list(csv.DictReader(open("analysis/oltaris_compare.csv")))
    r = [x for x in rows if x["spectrum"] == spectrum and x["material"] == "Al" and x["point"] == "20 + water 5 (BFO)"][0]
    assert abs(float(r["D_ratio"]) - 1) < 0.03
    F = oc.SPEC[spectrum][0] or oc.tylka_integral()
    g = oc.gate_point(F, "Al", 20, "D", (4.5, 5.5))
    assert g == pytest.approx(float(r["GATE_D_mGy"]), rel=1e-4)


def test_material_ratios_at_bfo():
    rows = {(r["spectrum"], r["point"]): r for r in csv.DictReader(open("analysis/oltaris_ratios.csv"))}
    for sp in ("Webber56", "Kold72", "Tylka89"):
        r = rows[(sp, "20 + water 5 (BFO)")]
        assert abs(float(r["PE_Al_D_GATE"]) - float(r["PE_Al_D_OLTARIS"])) < 0.02
        assert abs(float(r["Pb_Al_D_GATE"]) - float(r["Pb_Al_D_OLTARIS"])) < 0.02


@pytest.mark.parametrize("tag", ["LaRC56", "LaRC72"])
def test_larc_boundary_flux_matches_wilson_denn(tag):
    """The boundary fluxes saved by OLTARIS follow the LaRC formulas (Wilson & Denn 1976)."""
    import oltaris_larc as ol
    E, phi, parts = oc.read_table(f"oltaris/Al_{tag}/li_flux.dat")
    p = phi[parts.index("proton")]
    m = p > 1e-20                      # OLTARIS writes 1e-20 where the fluence underflows
    assert m.sum() > 90
    np.testing.assert_allclose(p[m], ol.larc_differential(tag, E[m]), rtol=1e-5)


def test_larc_integral_fluence():
    """Integral of the tabulated flux (log-log interpolation on the OLTARIS grid) against the analytic value;
    the interpolation lies slightly below a pure exponential, by 0.4 % at 100 MeV for the 1972 fit."""
    import oltaris_larc as ol
    F56, _ = ol.flux_integral("LaRC56"); F72, _ = ol.flux_integral("LaRC72")
    assert F56(100) == pytest.approx(1.5e9 * np.exp(-90 / 25) + 3e8, rel=2e-3)
    assert F72(100) == pytest.approx(6.6e8, rel=5e-3)


def test_rbe_check_hzetrn():
    """Particle-specific RBE (NCRP 132) changes the HZETRN BFO dose by less than 2 % for every case."""
    import rbe_check as rc
    for lab, tags in rc.CASES:
        for mat in ("Al", "PE", "Pb"):
            g = rc.shares(mat, tags)
            if g is None: continue
            assert sum(g.values()) == pytest.approx(1.0, abs=1e-9)
            assert abs(sum(rc.RBE[k] * g[k] for k in rc.RBE) / 1.5 - 1) < 0.02
