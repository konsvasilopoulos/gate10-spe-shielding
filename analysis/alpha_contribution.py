"""alpha_contribution.py - Contribution of SPE alpha particles (set 6) relative to protons, for the same spectra per nucleon.
The alpha runs use He/H = 1 at the same MeV/n, so D_total = D_p + r * D_alpha (r = actual He/H ratio).
Reports D_alpha/D_p and H_alpha/H_p (skin, BFO, mean), the increase for r = 0.03 and 0.1, and the PE/Al and Pb/Al
ratios with and without alphas. Writes analysis/alpha_contribution.csv."""
import os, csv
W = {"skin": "skin_0-0.3", "BFO": "5_gcm2", "mean": "mean_0-30"}
SP = {"kold1956": "spe_runs/materials_kold1956", "kold1972": "spe_runs/materials_kold1972", "kold1989": "spe_runs/materials_kold1989"}
MATS = ["Al", "POLYETHYLENE", "Pb"]; T = [5.0, 10.0, 20.0]; RS = [0.03, 0.10]
def load(f):
    return {(r["material"], float(r["areal_g_cm2"])): r for r in csv.DictReader(open(f))} if os.path.exists(f) else {}
rows = []
for sp, pb in SP.items():
    P = load(f"{pb}/materials_summary.csv"); A = load(f"spe_runs/alpha_{sp}/materials_summary.csv")
    if not A: print(f"[!] {sp}: alpha runs missing"); continue
    print(f"\n=== {sp}: alpha / protons at the same MeV/n (He/H = 1) ===")
    print(f"{'material':<14}{'g/cm2':>6}" + "".join(f"{'D '+w:>10}{'H '+w:>10}" for w in W) + "   increase of mean H (r=0.03 / 0.1)")
    for m in MATS:
        for t in T:
            if (m, t) not in A or (m, t) not in P: continue
            p, a = P[(m, t)], A[(m, t)]; line = f"{m:<14}{t:>6.0f}"; rec = [sp, m, t]
            for w, k in W.items():
                rd = float(a[f"D_cGy_{k}"]) / float(p[f"D_cGy_{k}"]); rh = float(a[f"H_cSv_{k}"]) / float(p[f"H_cSv_{k}"])
                line += f"{rd:>10.3f}{rh:>10.3f}"; rec += [rd, rh]
            rhm = rec[-1]; line += "   " + " / ".join(f"+{100*r*rhm:.1f} %" for r in RS)
            print(line); rows.append(rec)
    print("  Ratios to Al (mean D and H), without alphas -> with alphas (r = 0.03 / 0.1):")
    for m in ("POLYETHYLENE", "Pb"):
        for t in T:
            if not all((x, t) in A and (x, t) in P for x in (m, "Al")): continue
            out = []
            for q in ("D", "H"):
                k = f"{q}_{'cGy' if q == 'D' else 'cSv'}_{W['mean']}"
                base = float(P[(m, t)][k]) / float(P[("Al", t)][k])
                withr = [(float(P[(m, t)][k]) + r * float(A[(m, t)][k])) / (float(P[("Al", t)][k]) + r * float(A[("Al", t)][k])) for r in RS]
                out.append(f"{q}: {base:.3f} -> {withr[0]:.3f} / {withr[1]:.3f}")
            print(f"    {m:<13}{t:>3.0f} g/cm2  " + "   ".join(out))
with open("analysis/alpha_contribution.csv", "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(["spectrum", "material", "areal_g_cm2"] + [f"{q}a_over_{q}p_{x}" for x in W for q in ("D", "H")]); w.writerows(rows)
print("\nWritten: analysis/alpha_contribution.csv")
