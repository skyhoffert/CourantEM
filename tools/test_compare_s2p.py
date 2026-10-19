"""Self-test for compare_s2p.py: one network written in RI/Hz, DB/GHz and MA/MHz with
a 75-ohm reference must compare as identical; a 0.1 dB / 2 deg perturbation must
fail the openEMS rule and pass the HFSS rule. Run: python3 tools/test_compare_s2p.py"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import compare_s2p as cs  # noqa: E402

f = np.linspace(0.5e9, 40e9, 80)
bl = 2 * np.pi * f * 1.65 * 0.01 / 299792458.0
# a mismatched line section: Z = 60 ohm, 10 mm, in 50 ohm
Z, R = 60.0, 50.0
S = np.empty((len(f), 2, 2), complex)
for k, t in enumerate(bl):
    A, B, C, D = np.cos(t), 1j * Z * np.sin(t), 1j * np.sin(t) / Z, np.cos(t)
    den = A + B / R + C * R + D
    S[k] = np.array([[A + B / R - C * R - D, 2 * (A * D - B * C)], [2, -A + B / R - C * R + D]]) / den


def write(path, S, fmt, unit, R0):
    Sr = cs.renorm(S, 50.0, R0)
    sc = {"HZ": 1, "MHZ": 1e6, "GHZ": 1e9}[unit]
    with open(path, "w") as fh:
        fh.write(f"! test\n# {unit} S {fmt} R {R0}\n")
        for k in range(len(f)):
            vals = []
            for r, c in ((0, 0), (1, 0), (0, 1), (1, 1)):
                s = Sr[k, r, c]
                if fmt == "RI":
                    vals += [s.real, s.imag]
                elif fmt == "MA":
                    vals += [abs(s), np.degrees(np.angle(s))]
                else:
                    vals += [20 * np.log10(abs(s)), np.degrees(np.angle(s))]
            fh.write(" ".join(f"{v:.12e}" for v in [f[k] / sc] + vals) + "\n")


d = tempfile.mkdtemp()
pa, pb, pc, pd = (os.path.join(d, n) for n in ("a.s2p", "b.s2p", "c.s2p", "d.s2p"))
write(pa, S, "RI", "HZ", 50.0)
write(pb, S, "DB", "GHZ", 50.0)
write(pc, S, "MA", "MHZ", 75.0)
fail = 0
for p in (pb, pc):
    res, _ = cs.compare(cs.read_s2p(pa), cs.read_s2p(p), "gpu")
    worst = max(c["value"] for c in res["checks"] if c["unit"] in ("dB", "deg"))
    ok = res["pass"] and worst < 1e-6
    print(f"{'ok  ' if ok else 'FAIL'} same network, {os.path.basename(p)}: worst difference {worst:.2e}")
    fail |= not ok
P = S.copy()
P[:, 1, 0] *= 10 ** (0.1 / 20) * np.exp(1j * np.radians(2))
write(pd, P, "RI", "HZ", 50.0)
r1, _ = cs.compare(cs.read_s2p(pa), cs.read_s2p(pd), "openems")
r2, _ = cs.compare(cs.read_s2p(pa), cs.read_s2p(pd), "hfss")
ok = (not r1["pass"]) and r2["pass"]
print(f"{'ok  ' if ok else 'FAIL'} 0.1 dB / 2 deg on S21: openEMS rule fails, HFSS rule passes")
fail |= not ok
print("compare_s2p test FAILED" if fail else "compare_s2p test passed")
sys.exit(1 if fail else 0)
