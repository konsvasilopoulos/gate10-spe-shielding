"""Threshold rule on the Tylka & Dietrich Band fits (Raukunen et al. 2018, Table 2): stored data and results."""
import csv
import numpy as np


def _rows(path):
    return list(csv.DictReader(l for l in open(path) if not l.startswith("#")))


def test_parameter_table():
    r = _rows("analysis/data/tylka_dietrich_band.csv")
    assert len(r) == 65                                        # 59 GLE components + 6 ESP components
    assert len({x["GLE"] for x in r}) == 59
    g5 = next(x for x in r if x["GLE"] == "5")
    assert (float(g5["J0"]), float(g5["g1"]), float(g5["g2"]), float(g5["R0"])) == (1.75e8, 1.76, 5.04, 0.566)


def test_band_continuity():
    import sys; sys.path.insert(0, "analysis")
    from gle_rule_td import band
    J0, g1, g2, R0 = 1.75e8, 1.76, 5.04, 0.566
    Rb = (g2 - g1) * R0; Eb = (np.sqrt((Rb * 1000) ** 2 + 938.272 ** 2) - 938.272)
    lo, hi = band(Eb * (1 - 1e-7), J0, g1, g2, R0), band(Eb * (1 + 1e-7), J0, g1, g2, R0)
    assert abs(lo / hi - 1) < 1e-5


def test_rule_on_second_reconstruction():
    r = [x for x in csv.DictReader(open("analysis/gle_rule_td.csv")) if x["material"] == "POLYETHYLENE" and x["areal"] == "20"]
    assert len(r) == 59
    g = np.array([float(x["gamma_at_Estar"]) for x in r]); f = np.array([float(x["fold_ratio"]) for x in r])
    q = f / np.array([float(x["rule_ratio"]) for x in r])
    assert np.corrcoef(g, np.log(f))[0, 1] < -0.99
    assert 0.98 <= q.min() and q.max() <= 1.05                # paper: within -1 to +4 %
