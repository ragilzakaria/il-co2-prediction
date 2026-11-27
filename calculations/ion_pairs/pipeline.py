"""Steps between the stages of the ion-pair job (roadmap 3.5). Plain Python 3.6, no packages,
so it runs on Mahameru's default python3. Called by submit.sh from this folder:

  python3 pipeline.py prepare-a   stage A inputs: every start of starts/<system>.xyz, GFN2-xTB opt
  python3 pipeline.py select-a    per system: converged starts with unchanged bonding, sorted by
                                  energy, duplicates removed (0.1 kcal/mol and 0.15 A RMSD of the sorted
                                  interatomic distances); the lowest 3 within 3 kcal/mol go to
                                  stage B (r2SCAN-3c opt)
  python3 pipeline.py select-b    per system: lowest converged stage-B structure with unchanged
                                  bonding; stage C inputs: wB97X-D3/def2-TZVPD single points of
                                  the complex and of host and guest at the complex geometry
  python3 pipeline.py results     interaction energies -> results.tsv

"Unchanged bonding": the set of bonded atom pairs (distance below the sum of covalent radii
+ 0.45 A) is the same as in the start, so no proton moved between the ions and no ion broke up.
Systems where every structure changed bonding are listed in the select-a / select-b output.
"""
import glob
import math
import os
import re
import sys

HARTREE_KCAL = 627.509
HARTREE_KJ = 2625.4996
COVALENT = {"H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05, "Cl": 1.02}
E_WINDOW_A, MAX_B = 3.0, 3
DUP_E, DUP_RMSD = 0.1, 0.15

# Fermi smearing at 1000 K: at the default 300 K the xTB SCF of some contact pairs did not converge
# (e.g. [EMIM][DCA]); stage A only ranks structures, the energies used come from stages B and C.
XTB_OPT = "! NATIVE-GFN2-XTB Opt\n%maxcore 1000\n%scf maxiter 500 SmearTemp 1000 end\n%geom maxiter 300 end"
R2SCAN_OPT = "! r2SCAN-3c Opt MiniPrint\n%maxcore 3000\n%geom maxiter 300 end"
SP = "! wB97X-D3 def2-TZVPD def2/J RIJCOSX TightSCF\n%maxcore 2500"


def systems():
    lines = open("systems.tsv").read().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[1:] if l.strip()]


def read_frames(path):
    lines = open(path).read().splitlines()
    frames, i = [], 0
    while i < len(lines) and lines[i].strip():
        n = int(lines[i])
        frames.append([(l.split()[0], float(l.split()[1]), float(l.split()[2]), float(l.split()[3]))
                       for l in lines[i + 2:i + 2 + n]])
        i += n + 2
    return frames


def read_xyz(path):
    return read_frames(path)[0]


def write_input(path, header, title, charge, atoms):
    body = ["# " + title, header, "* xyz %d 1" % charge]
    body += ["%-2s %14.8f %14.8f %14.8f" % a for a in atoms]
    body.append("*")
    with open(path, "w") as f:
        f.write("\n".join(body) + "\n")


def dist(a, b):
    return math.sqrt((a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2)


def bonds(atoms):
    out = set()
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if dist(atoms[i], atoms[j]) < COVALENT[atoms[i][0]] + COVALENT[atoms[j][0]] + 0.45:
                out.add((i, j))
    return out


def dmat_rmsd(a, b):
    """RMSD of the sorted lists of all interatomic distances. Sorting makes it blind to swapped
    equivalent atoms (the four F of [BF4], methyl H): comparing distances atom by atom made the same
    structure look 0.6-1.0 A different, so stage B got three copies of one structure."""
    da = sorted(dist(a[i], a[j]) for i in range(len(a)) for j in range(i + 1, len(a)))
    db = sorted(dist(b[i], b[j]) for i in range(len(b)) for j in range(i + 1, len(b)))
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(da, db)) / len(da))


def finished_opt(out):
    if not os.path.exists(out):
        return None
    text = open(out, errors="ignore").read()
    if "ORCA TERMINATED NORMALLY" not in text or "THE OPTIMIZATION HAS CONVERGED" not in text:
        return None
    return float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)[-1])


def final_sp(out):
    if not os.path.exists(out):
        return None
    text = open(out, errors="ignore").read()
    if "ORCA TERMINATED NORMALLY" not in text:
        return None
    return float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)[-1])


def prepare_a():
    os.makedirs("A/inputs", exist_ok=True)
    n = 0
    for s in systems():
        for k, atoms in enumerate(read_frames("starts/%s.xyz" % s["system"])):
            write_input("A/inputs/%s__s%02d.inp" % (s["system"], k), XTB_OPT,
                        "%s start %d: GFN2-xTB optimisation" % (s["system"], k), int(s["charge"]), atoms)
            n += 1
    print("stage A: %d inputs" % n)


def select_a():
    os.makedirs("B/inputs", exist_ok=True)
    rows = ["system\tn_starts\tn_converged\tn_bonding_changed\tn_to_B\tstarts_to_B\txtb_rel_kcal"]
    for s in systems():
        starts = read_frames("starts/%s.xyz" % s["system"])
        ref = bonds(starts[0])
        cands, n_conv, n_changed = [], 0, 0
        for k in range(len(starts)):
            name = "%s__s%02d" % (s["system"], k)
            e = finished_opt("A/outputs/%s.out" % name)
            if e is None:
                continue
            n_conv += 1
            atoms = read_xyz("A/outputs/%s.xyz" % name)
            if bonds(atoms) != ref:
                n_changed += 1
                continue
            cands.append((e, k, atoms))
        cands.sort()
        chosen = []
        for e, k, atoms in cands:
            if (e - cands[0][0]) * HARTREE_KCAL > E_WINDOW_A or len(chosen) == MAX_B:
                break
            if any(abs(e - ce) * HARTREE_KCAL < DUP_E and dmat_rmsd(atoms, ca) < DUP_RMSD for ce, _, ca in chosen):
                continue
            chosen.append((e, k, atoms))
        for j, (e, k, atoms) in enumerate(chosen):
            write_input("B/inputs/%s__b%d.inp" % (s["system"], j), R2SCAN_OPT,
                        "%s: r2SCAN-3c optimisation of GFN2-xTB structure from start %d (%.2f kcal/mol above "
                        "the lowest)" % (s["system"], k, (e - cands[0][0]) * HARTREE_KCAL),
                        int(s["charge"]), atoms)
        rows.append("%s\t%d\t%d\t%d\t%d\t%s\t%s" % (
            s["system"], len(starts), n_conv, n_changed, len(chosen), ",".join(str(k) for _, k, _ in chosen),
            ",".join("%.2f" % ((e - cands[0][0]) * HARTREE_KCAL) for e, _, _ in chosen)))
        if not chosen:
            print("[NO-B] %s: %d of %d starts converged, %d changed bonding" % (
                s["system"], n_conv, len(starts), n_changed))
    open("selected_A.tsv", "w").write("\n".join(rows) + "\n")
    print("stage B: %d inputs" % len(glob.glob("B/inputs/*.inp")))


def select_b():
    os.makedirs("C/inputs", exist_ok=True)
    rows = ["system\tbest\tE_r2scan_Eh\tn_converged\tn_bonding_changed\tspread_kcal"]
    for s in systems():
        ref = bonds(read_frames("starts/%s.xyz" % s["system"])[0])
        nh, charge = int(s["n_host_atoms"]), int(s["charge"])
        host_q = 1 if s["kind"] == "pair" else -1
        guest_q = charge - host_q
        res, n_changed = [], 0
        for f in sorted(glob.glob("B/outputs/%s__b*.out" % glob.escape(s["system"]))):
            e = finished_opt(f)
            if e is None:
                continue
            atoms = read_xyz(f[:-4] + ".xyz")
            if bonds(atoms) != ref:
                n_changed += 1
                continue
            res.append((e, os.path.basename(f)[:-4], atoms))
        if not res:
            print("[NO-C] %s: no converged r2SCAN-3c structure with unchanged bonding" % s["system"])
            continue
        res.sort()
        e, best, atoms = res[0]
        sysname = s["system"]
        write_input("C/inputs/%s__complex.inp" % sysname, SP, "%s: complex, %s" % (sysname, best), charge, atoms)
        write_input("C/inputs/%s__host.inp" % sysname, SP, "%s: host at the complex geometry" % sysname,
                    host_q, atoms[:nh])
        write_input("C/inputs/%s__guest.inp" % sysname, SP, "%s: guest at the complex geometry" % sysname,
                    guest_q, atoms[nh:])
        with open("C/inputs/%s.best.xyz" % sysname, "w") as f:
            f.write("%d\n%s r2SCAN-3c %.10f\n" % (len(atoms), best, e))
            f.write("\n".join("%-2s %14.8f %14.8f %14.8f" % a for a in atoms) + "\n")
        rows.append("%s\t%s\t%.10f\t%d\t%d\t%.2f" % (sysname, best, e, len(res), n_changed,
                                                     (res[-1][0] - e) * HARTREE_KCAL))
    open("selected_B.tsv", "w").write("\n".join(rows) + "\n")
    print("stage C: %d inputs" % len(glob.glob("C/inputs/*.inp")))


def results():
    rows = ["system\tkind\thost\tguest\tE_complex_Eh\tE_host_Eh\tE_guest_Eh\tE_int_kJmol"]
    for s in systems():
        e = [final_sp("C/outputs/%s__%s.out" % (s["system"], p)) for p in ("complex", "host", "guest")]
        if None in e:
            continue
        rows.append("%s\t%s\t%s\t%s\t%.10f\t%.10f\t%.10f\t%.2f" % (
            s["system"], s["kind"], s["host"], s["guest"], e[0], e[1], e[2], (e[0] - e[1] - e[2]) * HARTREE_KJ))
    open("results.tsv", "w").write("\n".join(rows) + "\n")
    print("interaction energies: %d of %d systems" % (len(rows) - 1, len(systems())))


if __name__ == "__main__":
    {"prepare-a": prepare_a, "select-a": select_a, "select-b": select_b, "results": results}[sys.argv[1]]()
