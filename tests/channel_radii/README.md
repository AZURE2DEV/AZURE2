# channel_radii

Per-channel radii: the `channel_radius` records of `<channelRadii>` (dev's `<channelRadii>` block of 2026-11-10,
written here as records; `docs/upgrade/dev_features.md` section 2).

AZURE2 had one channel radius per particle pair (column 28 of the `<levels>` lines). The formal theory allows a
radius per channel c = (pair, l, s, J), and elastic-scattering fits in particular can want different radii for
different partial waves. The optional section

```
<channelRadii>
channel_radius p12C_l1            # every l = 1 channel of p+12C at 4.0 fm
    pair    p+12C
    l       1
    radius  4.0 fm
channel_radius p12C_l2_J5_2       # l = 2, J = 5/2 only
    pair    p+12C
    l       2
    radius  4.5 fm
    J       5/2
</channelRadii>
```

gives those channels their own radius; every other channel keeps its pair's. The most specific record wins
(l+s+J, then l+J or l+s, then l); a channel named by a record with s and one with J, and by none with both, is refused.
Everything taken at the channel radius follows the channel: penetrability, shift, hard-sphere phase, boundary
condition, Wigner limit, the formal <-> observed transformation, the ANC conversion and N_f integral of an
external-capture final state, the radial capture integrals and the s-wave scattering length. External capture
requires the two channels of every transition to share one radius (DeBoer, 2026-11-10: with two radii it is not
defined where the external region starts); the set-up refuses a model that breaks this. The signature of the
external-capture integrals' file includes a channel's radius where it is not its pair's.

`AZURE2 --convert` turns dev's block (lines `<pair key> <l> <radius> [s=<s>] [J=<J>]`) into these records.

## Checks (`check.sh path/to/AZURE2`)

Truth without the feature: in p + 12C every J-group is a single partial wave, and a phase shift depends only on its
own channel's radius, so a record for l = 1 must reproduce, wave by wave, either the unchanged model or the model
with the pair radius changed outright -- on the A- and R-matrix routes, in the Brune basis and with `--no-brune`.

- (a) l = 1 record: 1/2- and 3/2- equal the a = 4.0 fm model, 1/2+, 3/2+, 5/2+ the a = 5.0 fm model (1e-9; the
  references differ);
- (b) `J 3/2`: only 3/2- moves;
- (c) every l at the pair radius == no record (identical); every l at 4.0 fm == pair radius 4.0 fm (cross sections);
- (d) 13N with external capture: every l of p+12C at 4.0 fm reproduces every point of the pair-radius-4.0 model
  (capture integrals, ANC conversion, target integration); l = 0 alone (E1 from l = 0 into the l = 1 final channel at
  another radius) is refused;
- (e) refused: a selection of channels given twice (`E_CHANNEL_RADIUS_TWICE`), a photon pair
  (`E_CHANNEL_RADIUS_PAIR_KIND`), a radius of zero, a record with s and a record with J on one channel; a record that
  names no channel is warned about and changes nothing.

Not here yet (dev has them): the session's channel form of `AZUREAPI::SetRadius`, and the per-channel radii in the
compound-nucleus check's report. THM: the vertex pieces are still taken at the pair's radius.
