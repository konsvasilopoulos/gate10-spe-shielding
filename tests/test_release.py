"""Checks of the comparison of GATE releases 10.0.0 and 10.1.1 (stored direct-simulation results)."""
import csv
import release_benchmark as rb


def test_absorbed_dose_independent_of_release():
    rows = list(csv.DictReader(open("analysis/release_benchmark_hzetrn.csv")))
    old = {(r["spectrum"], r["material"], r["point"]): float(r["GATE_D_mGy"]) for r in rows if r["version"] == "v10.0.0"}
    new = {(r["spectrum"], r["material"], r["point"]): float(r["GATE_D_mGy"]) for r in rows if r["version"] == "v10.1.1"}
    for k, v in old.items():
        if k[2] in ("entrance", "0.3", "5 (BFO)"):
            assert abs(new[k] / v - 1) < 0.025, k


def test_gate_point_reads_release_data():
    d, h = rb.gate_point("v10.1.1", "kold1972", "Al", (4.5, 5.5))
    assert 100 < d < 150 and h > d                       # mGy, mSv behind 20 g/cm2 Al
