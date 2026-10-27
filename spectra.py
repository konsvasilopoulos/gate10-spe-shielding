"""
spectra.py - SPE proton spectra for the GATE source.

Webber 1956 (February 1956), as in NASA TP-2009-215560 (Heinbockel et al.):
    integral fluence  Phi(>E) = 1e9 * exp[(239.1 - P(E)) / 100]   protons / cm2
    with P(E) = sqrt(E (E + 1876))  (rigidity in MV, E kinetic energy in MeV)
    P0 = 100 MV. The normalisation gives Phi(>30 MeV) ~ 1e9 cm^-2.

King 1972 and the Koldobskiy et al. (2021) reconstructions are defined below.

Use from another script:
    from spectra import spectrum_bins
    edges, weights, fluence = spectrum_bins("webber1956", emin=10, emax=2000)

Stand-alone:  python spectra.py   -> prints basic numbers and writes spectrum_webber1956.png
"""
import numpy as np

MP2 = 1876.544  # 2 m_p c^2 (MeV)


def rigidity_MV(E):
    E = np.asarray(E, dtype=float)
    return np.sqrt(E * (E + MP2))


def webber1956_integral(E):
    """Phi(>E) in protons/cm2 for the whole event."""
    return 1e9 * np.exp((239.1 - rigidity_MV(E)) / 100.0)


def webber1956_differential(E):
    """dPhi/dE in protons/(cm2 MeV) = -dPhi(>E)/dE."""
    E = np.asarray(E, dtype=float)
    dPdE = (E + MP2 / 2) / rigidity_MV(E)
    return webber1956_integral(E) * dPdE / 100.0


def king1972_integral(E):
    """August 1972 SPE (King 1974): Phi(>E) = 7.9e9 * exp[(30 - E)/26.5] cm^-2.
    The fit was made for 30-200 MeV; here it is used for 10-1000 MeV, as is common in the literature."""
    E = np.asarray(E, dtype=float)
    return 7.9e9 * np.exp((30.0 - E) / 26.5)


# --- Koldobskiy et al. 2021 (A&A 647, A132): modified Band function per GLE ------------------------
# F(>R) = J1 (R/1GV)^-g1 exp(-R/R1)  for R < Rb,   J2 (R/1GV)^-g2 exp(-R/R2)  for R >= Rb   (cm^-2, R in GV)
# g0 = g2 - g1, R0 = R1 R2/(R2 - R1), Rb = g0 R0, J1 = J2 Rb^-g0 exp(g0).  Valid for E > 30 MeV.
# Parameters from their Table 2 (J2, g1, R1, g2, R2); R2 = None means +infinity (> 100 GV).
KOLDOBSKIY_GLE = {
    5:  (1.63e8, 1.59, 0.770, 4.84, 8.614),     # 23 Feb 1956
    24: (4.35e6, 2.02, 0.152, 11.96, None),     # 04 Aug 1972
    25: (2.60e5, 1.53, 0.039, 5.07, None),      # 07 Aug 1972
    42: (9.89e6, -1.95, 0.069, 3.39, 5.935),    # 29 Sep 1989
    43: (6.59e6, -1.21, 0.109, 4.44, 2.960),    # 19 Oct 1989
    44: (1.84e6, 0.35, 0.122, 6.07, 5.807),     # 22 Oct 1989
    45: (1.19e7, -1.77, 0.077, 2.95, 1.131),    # 24 Oct 1989
}


def koldobskiy_gle_integral(E, gle):
    J2, g1, R1, g2, R2 = KOLDOBSKIY_GLE[gle]
    R = rigidity_MV(E) / 1000.0                      # GV
    g0 = g2 - g1
    R0 = R1 if R2 is None else R1 * R2 / (R2 - R1)
    Rb = g0 * R0
    J1 = J2 * Rb ** (-g0) * np.exp(g0)
    cut2 = 0.0 if R2 is None else R / R2
    low = J1 * R ** (-g1) * np.exp(-R / R1)
    high = J2 * R ** (-g2) * np.exp(-cut2)
    return np.where(R < Rb, low, high)


def koldobskiy_sum(gles):
    return lambda E: sum(koldobskiy_gle_integral(E, g) for g in gles)


def logflat_integral(E):
    """Auxiliary spectrum for the response matrix: dN/dlnE constant, Phi(>E) = ln(5000 MeV / E).
    The fluence in each narrow window is ln(E2/E1); response = dose / window fluence."""
    return np.log(5000.0 / np.asarray(E, dtype=float))


SPECTRA = {"logflat": logflat_integral, "webber1956": webber1956_integral, "king1972": king1972_integral,
           "kold1956": koldobskiy_sum([5]),             # Feb 1956 (GLE 5)
           "kold1972": koldobskiy_sum([24, 25]),        # Aug 1972 (GLE 24 + 25)
           "kold1989": koldobskiy_sum([42, 43, 44, 45])}  # Sep-Oct 1989 (GLE 42-45)


def spectrum_bins(name="webber1956", emin=10.0, emax=2000.0, nbins=300):
    """Logarithmic bins. Weight of each bin = exact fluence inside the bin (from the integral form).
    Returns (bin edges in MeV, weights in cm^-2, total fluence emin-emax in cm^-2).
    Note: opengate expects a density per MeV; run_spectrum.py divides by the bin width."""
    integral = SPECTRA[name]
    edges = np.geomspace(emin, emax, nbins + 1)
    w = integral(edges[:-1]) - integral(edges[1:])
    return edges, w, float(w.sum())


if __name__ == "__main__":
    edges, w, tot = spectrum_bins()
    print("Webber 1956, 10-2000 MeV")
    print(f"  Total fluence:           {tot:.3e} protons/cm2")
    for e in (30, 100, 180, 300, 500, 1000):
        print(f"  Phi(>{e:>4} MeV) = {webber1956_integral(e):.3e} cm^-2   ({100*webber1956_integral(e)/tot:6.3f} % of total)")
    # check: integral of the differential = integral spectrum
    E = np.geomspace(10, 2000, 20001)
    num = np.trapezoid(webber1956_differential(E), E) if hasattr(np, "trapezoid") else np.trapz(webber1956_differential(E), E)
    print(f"  Check: int dPhi/dE = {num:.4e}, Phi(>10)-Phi(>2000) = {tot:.4e}")
    try:
        import matplotlib.pyplot as plt
        Ec = np.sqrt(edges[:-1] * edges[1:])
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].loglog(Ec, webber1956_differential(Ec)); ax[0].set_xlabel("E (MeV)"); ax[0].set_ylabel(r"d$\Phi$/dE (cm$^{-2}$ MeV$^{-1}$)")
        ax[1].loglog(Ec, webber1956_integral(Ec)); ax[1].set_xlabel("E (MeV)"); ax[1].set_ylabel(r"$\Phi$(>E) (cm$^{-2}$)")
        for a in ax: a.grid(True, which="both", alpha=0.3)
        fig.suptitle("February 1956 SPE (Webber)"); fig.tight_layout(); fig.savefig("spectrum_webber1956.png", dpi=150)
        print("  -> spectrum_webber1956.png")
    except ImportError:
        pass
