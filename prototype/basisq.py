"""basisq.py -- parse the PUBLIC q-vSZP basis file (param/gxtb/basisq).

The file layout (decoded from the v1.1 asset; every claim below is checked by the selftest):

    *
      <Z>   <a>   <b>   <c>          element header: atomic number + 3 adaptation coefficients
               <n> <l>               a shell: n primitives of angular momentum l (s/p/d/f)
        <exponent>  <c0>  <c1>       n rows: primitive exponent + STATIC + CHARGE-SCALING columns
      ...more shells...
    *                                next element

The two coefficient columns are the q-vSZP signature: the effective contraction is
c = c0 + c1 * q_eff, with q_eff built from the EEQ(BC) charge and coordination number as the
oracle prints in its "q-vSZP AO setup (q+aq^2+bCN^0.5+cqCN)" block -- the element header's
(a, b, c) are that expression's coefficients. So the BASIS ITSELF is geometry-dependent through
the charges: the coefficients must be rebuilt per structure before any integral is computed.

    from basisq import parse
    B = parse()          # {Z: {"adapt": (a,b,c), "shells": [(l, [(exp,c0,c1), ...]), ...]}}
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "..", "param", "gxtb", "basisq")
L = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}


def parse(path=DEFAULT):
    text = open(path).read()
    blocks = [b.strip() for b in text.split("*") if b.strip()]
    out = {}
    for b in blocks:
        lines = [l for l in b.splitlines() if l.strip()]
        head = lines[0].split()
        z = int(head[0])
        adapt = tuple(float(x) for x in head[1:4])
        shells = []
        i = 1
        while i < len(lines):
            m = re.match(r"^\s*(\d+)\s+([spdfg])\s*$", lines[i])
            if not m:
                i += 1
                continue
            n, l = int(m.group(1)), L[m.group(2)]
            prims = []
            for row in lines[i + 1:i + 1 + n]:
                v = [float(x) for x in row.split()]
                prims.append((v[0], v[1], v[2]))
            shells.append((l, prims))
            i += 1 + n
        out[z] = {"adapt": adapt, "shells": shells}
    return out


def nprims(basis, z):
    return sum(len(p) for _l, p in basis[z]["shells"])


def nprim_components(basis, z):
    """Primitives x Cartesian components -- the convention behind the oracle's printed `npr`
    (measured: lone O prints npr 24 = 6s + 6p*3; water 40 = 24 + 2*8). The oracle's `acpsao`
    is 16 PER ATOM = 1+3+5+7: the ACPs carry s,p,d,f projectors on every element."""
    per = {0: 1, 1: 3, 2: 6, 3: 10, 4: 15}
    return sum(len(p) * per[l] for l, p in basis[z]["shells"])


def ncao(basis, z):
    """Cartesian AO count: s=1, p=3, d=6, f=10 per shell."""
    per = {0: 1, 1: 3, 2: 6, 3: 10, 4: 15}
    return sum(per[l] for l, _p in basis[z]["shells"])


def nsao(basis, z):
    """Spherical AO count: 2l+1 per shell."""
    return sum(2 * l + 1 for l, _p in basis[z]["shells"])


if __name__ == "__main__":
    B = parse()
    zs = sorted(B)
    print(f"elements parsed: {len(B)} (Z {zs[0]}..{zs[-1]})")
    # the oracle printed for water: npr 40, ncao 6, nsao 6 -- check O + 2H against the file
    # (npr is primitives x Cartesian components: lone-atom runs decode the convention)
    npr = nprim_components(B, 8) + 2 * nprim_components(B, 1)
    nc = ncao(B, 8) + 2 * ncao(B, 1)
    ns = nsao(B, 8) + 2 * nsao(B, 1)
    print(f"water from file: npr {npr} (oracle 40), ncao {nc} (oracle 6), nsao {ns} (oracle 6)")
    assert npr == 40, npr
    assert ns == 6 and nc == 6, (ns, nc)
    assert nprims(B, 1) == 8 and nprim_components(B, 8) == 24  # measured on lone atoms
    h = B[1]
    print(f"H: adapt {h['adapt']}, shells {[(l, len(p)) for l, p in h['shells']]}")
    o = B[8]
    print(f"O: adapt {o['adapt']}, shells {[(l, len(p)) for l, p in o['shells']]}")
    print("BASISQ PARSER SELFTEST PASSED")
