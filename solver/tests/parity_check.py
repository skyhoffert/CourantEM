#!/usr/bin/env python3
"""Compare two parity outputs (native vs WebAssembly). Relative tolerance 1e-12:
the builds differ only in libm rounding of exp/sin/cos/pow."""
import sys

a = [l.split() for l in open(sys.argv[1])]
b = [l.split() for l in open(sys.argv[2])]
assert len(a) == len(b), "line count differs"
worst = 0.0
for la, lb in zip(a, b):
    assert la[0] == lb[0]
    for x, y in zip(la[1:], lb[1:]):
        x, y = float(x), float(y)
        scale = max(abs(x), abs(y), 1e-300)
        worst = max(worst, abs(x - y) / scale)
print(f"max relative difference {worst:.2e}")
sys.exit(0 if worst < 1e-12 else 1)
