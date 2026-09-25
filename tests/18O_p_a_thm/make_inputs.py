"""Regenerate the inputs of tests/18O_p_a_thm from the reproduction of
La Cognata, Spitaleri & Mukhamedzhanov, ApJ 723 (2010) 1512.

Not run by the test.  Needs the reproduction tree (THM_REPRO, default
/home/almalinux/thm_repro: model/rmat.py for the Coulomb functions and the
formal->Brune conversion, azure/build_azr.py for the level lines, digitize/ for
the digitized Fig. 4) and numpy, scipy, mpmath.

Writes, next to this file:
  18O_p_a_thm.azr            Table 3 row R (B_c = S_c(E1)) as Brune parameters,
                             two THM segments, 17 keV (c.m. sigma) folding
  data/lc723_thm_points.dat  segment 1: the digitized THM S(E) points of Fig. 4
  data/lc723_band_mid.dat    segment 2: the digitized band mid-line
  lc723_band.txt             reference curve read by check.sh
Both data files are S(E) turned into the HOES observable AZURE2 computes,
  sigma_HOES ~ (S - background) / (P_0(E) exp(2 pi eta))   (constant E/k^2 dropped),
at the proton lab energy.
"""
import os, sys
import numpy as np
REPRO = os.environ.get('THM_REPRO', '/home/almalinux/thm_repro')
sys.path.insert(0, os.path.join(REPRO, 'model'))
sys.path.insert(0, os.path.join(REPRO, 'azure'))
from rmat import CH_P, M_P, M_18O
import build_azr as B

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(REPRO, 'digitize')
up = np.loadtxt(os.path.join(D, 'lc723_band_upper.txt'))
lo = np.loadtxt(os.path.join(D, 'lc723_band_lower.txt'))
bgc = np.loadtxt(os.path.join(D, 'lc723_background.txt'))
pts = np.loadtxt(os.path.join(D, 'lc723_data.txt'))


def conv(E):   # S(E) = sigma_HOES * P_0(E) * exp(2 pi eta), up to a constant
    return np.array([CH_P.PS(e)[0] * np.exp(2 * np.pi * CH_P.keta(e)[1]) for e in E])


def lab(E):
    return E * (M_P + M_18O) / M_18O


# reference band on the grid of the reproduction (0.505-0.895 MeV, 5 keV)
Eg = np.linspace(0.505, 0.895, 79)
U = np.interp(Eg, up[:, 0], up[:, 1])
L = np.interp(Eg, lo[:, 0], lo[:, 1])
mid, half, bg, cg = 0.5 * (U + L), 0.5 * (U - L), np.polyval(bgc, Eg), conv(Eg)
with open(os.path.join(HERE, 'lc723_band.txt'), 'w') as f:
    f.write('# La Cognata et al., ApJ 723 (2010) 1512, Fig. 4: digitized red band (fits to the upper and\n'
            '# lower data limits) on a 5 keV grid; linear background S = %.6g E %+.6g MeV b digitized from\n'
            '# the dashed line; P_0(E) exp(2 pi eta) (p+18O, l=0, a=5.1 fm) turns AZURE2 sigma_HOES into\n'
            '# the S(E) shape of the figure.\n'
            '# E_cm(MeV)  S_mid(MeV b)  half_width(MeV b)  background(MeV b)  P0*exp(2pi eta)\n' % tuple(bgc))
    for r in zip(Eg, mid, half, bg, cg):
        f.write('%.4f %10.3f %9.3f %9.3f %.8e\n' % r)
with open(os.path.join(HERE, 'data', 'lc723_band_mid.dat'), 'w') as f:
    for e, m, h, b, c in zip(Eg, mid, half, bg, cg):
        f.write('%.6f 0.0 %.6e %.6e\n' % (lab(e), (m - b) / c, h / c))

Ep = pts[:, 0]
cp = conv(Ep)
err = 0.5 * (np.abs(pts[:, 2]) + np.abs(pts[:, 3]))
with open(os.path.join(HERE, 'data', 'lc723_thm_points.dat'), 'w') as f:
    for e, s, de, c in zip(Ep, pts[:, 1], err, cp):
        f.write('%.6f 0.0 %.6e %.6e\n' % (lab(e), (s - np.polyval(bgc, e)) / c, de / c))

sig_lab = 0.017 * (M_P + M_18O) / M_18O
seg = '1  1  2  %.4f  %.4f  0  180  10  1.0  1  0  0  0  0  %s  0  0'
txt = [B.CONFIG, '<levels>', B.level_lines(B.params('R', 1)), '</levels>', '<segmentsData>',
       seg % (lab(Ep[0]) - 0.001, lab(Ep[-1]) + 0.001, 'data/lc723_thm_points.dat'),
       seg % (lab(Eg[0]) - 0.001, lab(Eg[-1]) + 0.001, 'data/lc723_band_mid.dat'),
       '</segmentsData>', '<segmentsTest>', '</segmentsTest>', '<targetInt>',
       '1              "1,2"          50             1              %.6f          0              0               '
       '"" 0               0              0               0               "" 0               0 0.04 5 50' % sig_lab,
       '</targetInt>']
with open(os.path.join(HERE, '18O_p_a_thm.azr'), 'w') as f:
    f.write('\n'.join(txt) + '\n')
