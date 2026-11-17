"""Pure-python checks of the sign bookkeeping (no AZURE2 needed)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import numpy as np
from rmfit.keys import ParamKey, canonical_signs, flip, sign_distance

keys = [ParamKey("width", "1/2+#1", 1, 0, 0.5), ParamKey("width", "1/2+#1", 2, 1, 0.5),
        ParamKey("width", "1/2+#2", 1, 0, 0.5), ParamKey("width", "1/2+#2", 2, 1, 0.5),
        ParamKey("width", "3/2-#1", 1, 1, 0.5), ParamKey("width", "3/2-#1", 2, 0, 0.5),
        ParamKey("energy", "3/2-#1"), ParamKey("norm", seg=3)]
x = np.array([1., 2., 3., -4., 5., 6., 7.5, 1.0])
lv = flip(x, keys, keys[:2])                        # whole-level flip: a symmetry
pr = flip(x, keys, [keys[0], keys[2], keys[4]])     # whole-pair flip: a symmetry
assert canonical_signs(lv, keys) == canonical_signs(x, keys)
assert canonical_signs(pr, keys) == canonical_signs(x, keys)
assert sign_distance(x, lv, keys)[0] == 0 and sign_distance(x, pr, keys)[0] == 0
d, diff = sign_distance(x, flip(x, keys, [keys[3]]), keys)
assert d == 1 and diff == [keys[3]]
# a fixed pair-1 amplitude on level 1/2+#1 ties that level's gauge to pair 1:
# the free pair-1 flip is then two flips away from the original
tie = [ParamKey("width", "1/2+#1", 1, 2, 0.5)]
assert canonical_signs(pr, keys, tie) != canonical_signs(x, keys, tie)
assert sign_distance(x, pr, keys, tie)[0] == 2
# distance is symmetric and ignores zero amplitudes
z = x.copy(); z[1] = 0.0
assert sign_distance(z, flip(z, keys, [keys[1], keys[3]]), keys)[0] == 1
assert sign_distance(x, pr, keys, tie)[0] == sign_distance(pr, x, keys, tie)[0]
print("test_keys: ok")
