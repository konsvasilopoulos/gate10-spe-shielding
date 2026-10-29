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
