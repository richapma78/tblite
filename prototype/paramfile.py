"""paramfile.py -- parse gxtb_parameters (plain ASCII) into per-element parameters.

The g-xTB parameter file the binary reads (~/.gxtb; pinned in param/gxtb/gxtb_parameters).
Layout decoded by matching gdb-extracted values against the columns (the SI PDF has the
equations, not these numbers). Structure:

  line 0,1 : two GLOBAL rows of 10.  g1 (line 0) = [.., c0=[4], c1=[5], k1'=[6], k2'=[7], ..];
             g2 (line 1) = [.., S3=[5], S5=[6], S7=[7], omega=[8], ..]  (AES kernel + range sep).
  then a 10-line block per element Z = 1..79 (H..Au):
     offset 0 : the integer Z (header)
     offset 1 : row0, 10 values   -- element-level; [6] = Gamma (Mulliken-3rd-order hardness)
     offset 2..7 : six SHELL rows of 4 (columns = s,p,d,f):
                   shell[3] = T25  (the MFX favg base, x_l = T25[l]*ipse)
                   shell[4] = T32  (the ES Hubbard U, gp3_gam2)
     offset 8,9 : two rows of 8   -- row8[5] = the atomic core-energy increment

ipse is NOT stored here (it is derived at load; supply it via ipse=... e.g. gdb-extracted).
Everything else the MFX gamma needs (T25, c[l]=g1[4+l]) comes straight from this file, for
every element -- which is why parsing it, not extracting element-by-element, is the path to
the metals.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PARAM_FILE = os.path.join(HERE, "..", "param", "gxtb", "gxtb_parameters")
IPSE_FILE = os.path.join(HERE, "data", "gxtb_ipse.json")   # derived (not in the param file)
ELEM_LINES = 10                                 # lines per element block


def load(path=PARAM_FILE, ipse_path=IPSE_FILE):
    rows = [ln.split() for ln in open(path) if ln.strip()]
    g1 = [float(x) for x in rows[0]]
    g2 = [float(x) for x in rows[1]]
    ipse = {}
    if os.path.exists(ipse_path):               # gdb-extracted; ipse is derived at load
        ipse = {int(z): v for z, v in json.load(open(ipse_path)).items()}
    elements = {}
    i = 2
    while i + ELEM_LINES <= len(rows):
        z = int(rows[i][0])
        row0 = [float(x) for x in rows[i + 1]]
        shell = [[float(x) for x in rows[i + 2 + k]] for k in range(6)]
        row8 = [float(x) for x in rows[i + 8]]
        row9 = [float(x) for x in rows[i + 9]]
        elements[z] = {
            "ipse": ipse.get(z),                # None until extracted/derived
            "gamma": row0[6],
            "T25": shell[3],                    # [s, p, d, f]
            "T32": shell[4],
            "core": row8[5],
            "row0": row0, "shell": shell, "row8": row8, "row9": row9,
        }
        i += ELEM_LINES
    return {"g1": g1, "g2": g2, "c": [g1[4], g1[5]], "omega": g2[8], "elements": elements}


if __name__ == "__main__":
    P = load()
    print(f"parsed {len(P['elements'])} elements (Z {min(P['elements'])}..{max(P['elements'])})")
    print(f"globals: c0={P['c'][0]:.10f} c1={P['c'][1]:.10f} omega={P['omega']:.10f}")
    # gate: parsed T25 must match the gdb-extracted values
    EXTRACTED = {  # (Z): [T25_s, T25_p]  gdb-extracted at 0x3be25c0 + Z*0x20
        1: [3.6548180166], 6: [3.1541670118, 2.3521912659],
        7: [5.6891941778, 1.8783447306], 8: [8.6375010152, 1.8397908212],
        9: [5.8665628795, 2.2071162860],
    }
    worst = 0.0
    for z, t25 in EXTRACTED.items():
        for l, v in enumerate(t25):
            d = abs(P["elements"][z]["T25"][l] - v)
            worst = max(worst, d)
    print(f"T25 parsed vs gdb-extracted (H,C,N,O,F): worst |diff| = {worst:.2e}  "
          f"{'PASS' if worst < 1e-12 else 'FAIL'}")
