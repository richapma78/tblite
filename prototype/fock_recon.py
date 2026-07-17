"""fock_recon.py -- the FULL AO Fock matrix, reconstructed from restart MOs + printed eps.

The H2 instrument (eht_h2.py) inverts a 2x2 -- one orbital per atom, one off-diagonal. This
generalizes: the generalized eigenproblem F C = S C eps with C^T S C = 1 inverts to

    F = S C diag(eps) C^T S

so with the restart file's MO matrix (record 2), the printed eigenvalue list, and OUR
machine-exact overlap, EVERY Fock element of ANY closed-shell molecule is measurable. That
unlocks the heteronuclear program: level-averaging rule, per-AO vs per-pair L4 scale, and
de-collinearized short-range candidates (HF vs HCl carry different CN radii).

SCOPE (measured): the oracle prints eigenvalues in a window of n_occ + 2, so the inversion is
complete only when nsao <= n_occ + 2. H2 (2/2), HF (5/5), F2 (8/8), H2O (6/6) qualify; HCl does
NOT (6 printed, nsao 10 -- chlorine's d shell starves the virtual list, and the binary exposes
no verbosity flag). The heteronuclear program therefore runs on the H/F system; HCl waits on an
eps-completion trick if it is ever needed.

Gates (declared before the first run; all four must pass on H2 + HF + F2):
  G-ortho  max|C^T S C - 1| < 1e-7, with OUR overlap in oracle AO order. This is STRONGER than
           the restart Mulliken gate: shell populations are invariant to within-shell AO
           permutations, orthonormality is not -- passing pins the p-ordering for the first time.
  G-dens   max|C diag(occ) C^T - P| < 1e-6 against the (already-gated) density record.
  G-h2     H2 at R=2.5: F_recon[0,1] = the banked 2x2 inversion value to 2e-5 (printed-eps
           precision), tying the new instrument to the old one where they overlap.
  G-sym    diatomic on z: the sigma block must not touch the pi AOs -- exactly ZERO p_x/p_y
           mixing with s rows (< 1e-6), which also self-identifies which p AO is p_z.

The MO-matrix ORIENTATION (record 2 row- vs column-major) is not assumed: both are tried and
G-ortho picks -- the choice is recorded in the result.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import constants as K  # noqa: E402
import overlap  # noqa: E402
import restart  # noqa: E402

_SYM2Z = {"h": 1, "he": 2, "c": 6, "n": 7, "o": 8, "f": 9, "si": 14, "p": 15, "s": 16, "cl": 17,
          "br": 35, "i": 53, "pd": 46}


def _pick_orient(Craw, S):
    best = None
    for tag, C in (("as-stored", Craw), ("transposed", Craw.T)):
        err = float(np.abs(C.T @ S @ C - np.eye(S.shape[0])).max())
        if best is None or err < best[1]:
            best = (tag, err, C)
    return best


def fock_ao(atoms, charge=0, uhf=0):
    """atoms in ANGSTROM (restart.converged_state convention). Closed shell: returns dict
    with F, S, C, eps (Eh), occ, ao labels, gate numbers. UKS (uhf > 0): returns the same
    plus F_a/F_b/C_a/C_b/eps_a/eps_b (per-spin full reconstructions, no symmetry assumption;
    layout decoded fifty-first push), with F = F_a for backward compatibility."""
    st = restart.converged_state(atoms, charge=charge, uhf=uhf)
    n = st["nsao"]
    eps = np.array(st["eps_ev"]) / K.EV
    occ = np.array(st["occ"])
    zs = [_SYM2Z[s.lower()] for s, _x, _y, _z in atoms]
    xyz = np.array([[x, y, z] for _s, x, y, z in atoms]) * K.BOHR
    shells, _ = overlap.build_shells(zs, xyz, charge=charge)
    S = overlap.overlap(zs, xyz, charge=charge, shells=shells, ao_order="oracle")
    assert S.shape[0] == n, f"our nsao {S.shape[0]} vs oracle {n}"
    labels = []
    for sh in shells:
        for m in range(2 * sh["l"] + 1):
            labels.append((sh["at"], zs[sh["at"]], sh["l"], m))
    if uhf == 0:
        assert len(eps) == n, f"printed eps has {len(eps)} entries, nsao {n}"
        assert len(occ) == n
        tag, ortho_err, C = _pick_orient(st["C"], S)
        dens_err = float(np.abs(C @ np.diag(occ) @ C.T - st["P"]).max())
        F = S @ C @ np.diag(eps) @ C.T @ S
        return {"F": F, "S": S, "C": C, "eps": eps, "occ": occ, "labels": labels,
                "orientation": tag, "ortho_err": ortho_err, "dens_err": dens_err,
                "state": st}
    # ---- UKS
    assert len(eps) == 2 * n, f"UKS eps has {len(eps)} entries, expected {2 * n}"
    assert "C_a" in st, "UKS restart records not found"
    out = {"S": S, "labels": labels, "state": st, "occ": occ}
    for sp, sl in (("a", slice(0, n)), ("b", slice(n, 2 * n))):
        tag, ortho_err, C = _pick_orient(st[f"C_{sp}"], S)
        e = eps[sl]
        F = S @ C @ np.diag(e) @ C.T @ S
        o = occ[sl]
        dens_err = float(np.abs(C @ np.diag(o) @ C.T - st[f"P_{sp}"]).max())
        out[f"F_{sp}"] = F
        out[f"C_{sp}"] = C
        out[f"eps_{sp}"] = e
        out[f"ortho_err_{sp}"] = ortho_err
        out[f"dens_err_{sp}"] = dens_err
    out["F"] = out["F_a"]
    out["ortho_err"] = max(out["ortho_err_a"], out["ortho_err_b"])
    out["dens_err"] = max(out["dens_err_a"], out["dens_err_b"])
    return out


def _gate():
    ok = True
    # ---- G-h2: tie to the 2x2 instrument
    R = 2.5
    r = fock_ao([("H", 0, 0, 0), ("H", 0, 0, R / K.BOHR)])
    ref_F12, ref_F11 = -0.24940460789490132, -0.2803005270651893   # banked (h2-eht-probes)
    d12, d11 = abs(r["F"][0, 1] - ref_F12), abs(r["F"][0, 0] - ref_F11)
    print(f"  H2   ortho {r['ortho_err']:.1e} ({r['orientation']})  dens {r['dens_err']:.1e}  "
          f"F12 {r['F'][0, 1]:+.7f} (banked {ref_F12:+.7f}, diff {d12:.1e})  "
          f"F11 diff {d11:.1e}", flush=True)
    ok &= r["ortho_err"] < 1e-7 and r["dens_err"] < 1e-6 and d12 < 2e-5 and d11 < 2e-5

    # ---- HF and F2: ortho + dens + sigma/pi separation
    for name, atoms in (("HF", [("H", 0, 0, 0), ("F", 0, 0, 0.92)]),
                        ("F2", [("F", 0, 0, 0), ("F", 0, 0, 1.41)])):
        r = fock_ao(atoms)
        F, labels = r["F"], r["labels"]
        hs = next(i for i, (a, _z, l, _m) in enumerate(labels) if a == 0 and l == 0)
        prows = [i for i, (a, _z, l, _m) in enumerate(labels) if a == 1 and l == 1]
        mix = sorted((abs(F[hs, i]), i) for i in prows)
        n_sigma = sum(1 for v, _i in mix if v > 1e-6)
        pz = mix[-1][1]
        print(f"  {name:4s} ortho {r['ortho_err']:.1e} ({r['orientation']})  "
              f"dens {r['dens_err']:.1e}  H(s)-p mixing: {n_sigma}/3 nonzero "
              f"(pi leak {max(mix[0][0], mix[1][0]):.1e}; p_z is AO index {pz}, "
              f"m={labels[pz][3]})  F[hs,pz] {F[hs, pz]:+.6f}", flush=True)
        ok &= r["ortho_err"] < 1e-7 and r["dens_err"] < 1e-6 and n_sigma == 1
    print("FOCK RECONSTRUCTION GATE " + ("PASSED" if ok else "FAILED"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    _gate()
