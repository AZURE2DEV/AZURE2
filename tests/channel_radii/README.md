# channel_radii

Per-channel radii: the `<channelRadii>` block (2026-11-10).

AZURE2 had one channel radius per particle pair (column 28 of the `<levels>` lines). The formal theory allows a
radius per channel c = (pair, l, s, J), and elastic-scattering fits in particular can want different radii for
different partial waves. The optional block, placed after `</levels>`,

```
<channelRadii>
1 1 4.0            # pair key 1, every l = 1 channel at 4.0 fm
1 2 4.5 J=5/2      # pair key 1, l = 2, J = 5/2 only
2 0 5.2 s=3/2      # pair key 2, l = 0, channel spin 3/2 only
</channelRadii>
```

gives those channels their own radius; every other channel keeps its pair's. The most specific line wins
(l+s+J > l+J > l+s > l); a channel matched by an s line and a J line but no l+s+J line is refused. Everything taken
at the channel radius follows the channel: penetrability, shift, hard-sphere phase, boundary condition, Wigner limit,
the formal <-> observed transformation, the ANC conversion and N_f integral of an external-capture final state, the
radial capture integrals and the s-wave scattering length. External capture requires the initial and final channels
of every transition to share one radius (DeBoer, 2026-11-10: with two radii it is not defined where the external
region starts); the engine refuses a model that breaks this. The intEC.dat signature includes the per-channel radii.

## Checks (`check.sh path/to/AZURE2`)

Truth without the feature: in p + 12C every J-group is a single partial wave, and a phase shift depends only on its
own channel's radius, so overriding the l = 1 radius must reproduce, wave by wave, either the unchanged model or the
model with the pair radius changed outright -- on the A- and R-matrix routes, with and without `--use-brune`.

- (a) l = 1 override: 1/2- and 3/2- equal the a = 4.0 fm model, 1/2+, 3/2+, 5/2+ the a = 5.0 fm model (identical;
  the references differ by up to 70 %);
- (b) `J=3/2`: only 3/2- moves;
- (c) every l at the pair radius == no block; every l at 4.0 fm == pair radius 4.0 fm (cross sections);
- (d) 13N with external capture: every l of pair 1 at 4.0 fm reproduces every point of the pair-radius-4.0 model
  (capture integrals, ANC conversion, target integration); l = 4 alone (no capture transition) runs and changes the
  result;
- (e) refused with an ERROR: l = 0 alone in 13N (E1 from l = 0 into the l = 1 final channel at another radius), an
  unknown pair key, a gamma pair, a non-positive radius, a malformed `s=`, a repeated line, an unterminated block,
  junk, an s line and a J line on one channel; a line that matches no channel is warned about and changes nothing;
- (f) the compound-nucleus check lists the per-channel radii.

Engines before the feature ignore the block silently: the installed dev 9b6ecb17 build fails 29 of the
checks.
