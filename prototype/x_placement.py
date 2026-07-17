"""x_placement.py -- pin the exchange-Fock placement at the 1e-3 level (offline).

The clean extraction (no forward fits): H2's H0 DIAGONAL is exactly -L2 (+ the measured L3
content), the mu channel is mu * (its measured response), ACP is analytic, so

    X11_req(R) = F11 - (-L2) - mu*d11_mu - L3*d11_L3 - ACP11        (1e-3 grade)
    X12_req(R) = F12 - H0fwd - mu*d12_mu - L3*d12_L3 - ACP12        (H0-forward grade)

H2's symmetry gives exact closed forms: per-spin per-atom Mulliken population m = 1/2 at
every R; P = [[1,1],[1,1]]/(1+S).

Placement candidates (v = gamma @ (m, m), the Mulliken potential):
    diag:  F11 = -w * v1                       -> w(R) = -X11_req/v1 must be R-CONSTANT
    off:   F12 = -w * S12 * (v1+v2)/2 + small  -> the remainder charted
PRE-DECLARED: the diagonal placement is PINNED iff w(R) is constant within +-2% over
R = 1..6; the off-diagonal remainder is charted, not tuned. Writes data/x-placement.json.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import params  # noqa: E402
import scf_h2  # noqa: E402

P_ = params.parse()
eH = P_["element"][1]
L2, MU, L3 = eH["shells"][0][0], eH["shells"][5][0], eH["shells"][1][0]


def main():
    ch = json.load(open(os.path.join(HERE, "data", "h2-channels.json")))
    rows = {r["R"]: r for r in
            json.load(open(os.path.join(HERE, "data", "h2-eht-probes.json")))["rows"]}
    out = []
    print("   R     X11_req    v1      w=-X11/v1  |  X12_req   -w*S*vbar   remainder")
    ws = []
    for Rs, rec in sorted(ch.items(), key=lambda t: float(t[0])):
        R = float(Rs)
        if R not in rows:
            continue
        r = rows[R]
        S, H0, A, _ES1, gam = scf_h2.build(R)
        g_off = gam[0, 1]
        v1 = 0.5 * (gam[0, 0] + g_off)
        X11 = r["F11"] - (-L2) - MU * rec["L7"]["d11"] - L3 * rec["L3"]["d11"] - A[0, 0]
        X12 = r["F12"] - H0[0, 1] - MU * rec["L7"]["d12"] - L3 * rec["L3"]["d12"] - A[0, 1]
        w = -X11 / v1
        ws.append(w)
        S12 = r["S12"]
        f12_c = -w * S12 * v1
        rem = X12 - f12_c
        out.append({"R": R, "X11_req": X11, "X12_req": X12, "v1": v1, "w": w,
                    "S12": S12, "remainder12": rem})
        print(f"  {R:4.2f}  {X11:+.5f}  {v1:+.5f}  {w:+.5f}   |  {X12:+.5f}  {f12_c:+.5f}"
              f"  {rem:+.5f}")
    wm = float(np.mean(ws))
    spread = (max(ws) - min(ws)) / abs(wm)
    print(f"\ndiag placement: w = {wm:+.5f}  spread {spread * 100:.1f}%  "
          f"({'PINNED' if spread <= 0.04 else 'NOT constant -- report'})   [1/2 = +0.50000]")
    json.dump({"points": out, "w_mean": wm, "w_spread": spread},
              open(os.path.join(HERE, "data", "x-placement.json"), "w"), indent=1)
    print("wrote data/x-placement.json")


if __name__ == "__main__":
    main()
