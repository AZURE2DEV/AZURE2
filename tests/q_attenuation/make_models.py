#!/usr/bin/env python3
"""Models for tests/q_attenuation (finite-geometry attenuation coefficients on a particle exit).

From tests/secondary_gamma/models/grad.azr (13C(a,n)16O*(6.13), two levels, n exit d and i waves; Q = -3.914 MeV,
so g = v_cm / v_n* = 0.33-0.37 near threshold):
  plain.azr  c.m. differential test segment 2 -> 5 (isDiff 5), E 5.9-6.3 step 0.1, theta* 0-180 step 5; data line inactive
  q.azr      the same + <targetInt> Q_0..Q_39 of a uniform cone of half-angle 15 degrees (Rose: cone average of P_K)
  q1.azr     the same + 40 coefficients all 1
  q5.azr     the same + only Q_0..Q_4 of the cone (as archive fits give them)
"""
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "secondary_gamma", "models", "grad.azr")
HALF_ANGLE = 15.0


def cone_q(n, half_angle_deg):
    u, w = np.polynomial.legendre.leggauss(200)
    ca = np.cos(np.radians(half_angle_deg))
    c = 0.5 * (1 - ca) * u + 0.5 * (1 + ca)
    return [float(np.sum(w * np.polynomial.legendre.legval(c, np.eye(k + 1)[k])) / np.sum(w)) for k in range(n)]


def target_int(q):
    return ('1 "1" 10 0 0 0 0 "" 0 1 %d %s 0 "" 0 0 0.04 20 50' % (len(q), " ".join("%.15g" % v for v in q)))


def main():
    s = open(SRC).read()
    data = re.search(r"<segmentsData>\n(.*?)\n</segmentsData>", s, re.S).group(1).split("\n")
    data = ["0" + line.lstrip()[1:] for line in data]          # inactive: the test segment is what is checked
    s = re.sub(r"<segmentsData>\n.*?\n</segmentsData>", "<segmentsData>\n" + "\n".join(data) + "\n</segmentsData>", s,
               flags=re.S)
    s = s.replace("<segmentsTest>\n</segmentsTest>", "<segmentsTest>\n1  2  5  5.9  6.3  0.1  0  180  5  5\n</segmentsTest>")
    q = cone_q(40, HALF_ANGLE)
    for name, tint in (("plain", ""), ("q", target_int(q)), ("q1", target_int([1.0] * 40)), ("q5", target_int(q[:5]))):
        t = s.replace("<targetInt>\n</targetInt>", "<targetInt>\n" + (tint + "\n" if tint else "") + "</targetInt>")
        open(os.path.join(HERE, "models", name + ".azr"), "w").write(t)
    print("cone half-angle", HALF_ANGLE, "Q_0..4", ["%.6f" % v for v in q[:5]])


if __name__ == "__main__":
    main()
