# The reaction system: `System` and `JpiBlock`

An R-matrix calculation is made for a stated set of partitions, channels and levels of one or more compound nuclei.
`physics::System` holds that set as numbers, and `physics::JpiBlock` holds the levels of one total angular momentum
and parity with their reduced width amplitudes as contiguous arrays. Every later stage of `physics/` (channel
functions, level parameters, collision matrix, observables; `ARCHITECTURE.md` §3) reads them and writes nothing back.
This page states what each entry is, in which unit, and how the entries are indexed. At present the engine still
computes from `CNuc`; the System is built from a loaded `CNuc` and compared with it, number for number. Its structure
(partitions, channels, blocks without levels; `CNucSystemContent::Structure`) is what `CNuc` hands to
`level_parameters` for the boundary values and the parameter transformations (`cnuc_split.md`).

## 1. Partitions

A partition α is a pair of particles into which the compound nucleus (charge number Z_CN, mass number A_CN) can
break up, each particle in a definite state. Particle 0 is the light one, particle 1 the heavy one:

| entry | symbol | unit | meaning |
|---|---|---|---|
| `kind` | | | `Particle` (two nuclei), `Photon` (photon + final state of the compound nucleus), `Beta` (lepton + parent nucleus) |
| `nucleus` | | | index of the compound nucleus (Z_CN, A_CN) in `System::nuclei()` |
| `Z[i]` | Z_i | | charge numbers (−1 for an electron) |
| `mass_u[i]` | m_i | u | atomic masses (not mass numbers; 0 for the photon) |
| `two_j[i]` | 2 j_i | | twice the intrinsic spins |
| `parity[i]` | π_i | | intrinsic parities, ±1 |
| `g_factor[i]` | g_i | | g-factors, μ_i = g_i j_i μ_N; they enter the M1 external-capture operator |
| `excitation_MeV[i]` | e_i | MeV | excitation energies of the two particles; for a photon partition e_1 is E_x of the final state |
| `separation_MeV` | S_α | MeV | separation energy of the ground-state pair; for a β partition the decay energy Q_β |

The threshold of the partition, in excitation energy of the compound nucleus, is E_x = S_α + e_0 + e_1, and the
c.m. kinetic energy in the partition at excitation energy E_x is E_c.m. = E_x − S_α − e_0 − e_1. The reduced mass
μ_α = m_0 m_1 / (m_0 + m_1) and the Sommerfeld parameter follow from these entries; they are functions of the
System, not stored in it.

## 2. Channels

A channel c = (α, l, s) is a partition with an orbital angular momentum l and a channel spin s, the coupling
s = j_0 + j_1 of the two intrinsic spins (Lane and Thomas, Rev. Mod. Phys. 30 (1958) 257, §III). For a photon channel
l is the multipolarity L, s is the spin of the final state, and the channel is electric or magnetic according to the
parity: EL when π_λ π_f = (−1)^L, ML otherwise. A β channel is Fermi (l = 0) or Gamow–Teller (l = 1).

| entry | symbol | unit | meaning |
|---|---|---|---|
| `partition` | α | | index into `System::partitions()` |
| `kind` | | | `Particle`, `Electric`, `Magnetic`, `Fermi`, `GamowTeller` |
| `l` | l (or L) | | orbital angular momentum, or multipolarity |
| `two_s` | 2 s | | twice the channel spin |
| `radius_fm` | a_c | fm | channel radius, where the internal and external regions are matched |

A System channel is identified by (α, kind, l, s): one entry serves every J^π block that contains it, so channel
functions P_c, S_c, φ_c can later be tabulated per (α, l) and E, not per block. Today every channel of a partition has
the partition's radius; the radius is stored per channel so that an l-dependent radius later changes no layout.

## 3. J^π blocks

The level matrix is block-diagonal in total angular momentum J, parity π and compound nucleus. A block holds:

| entry | symbol | unit | meaning |
|---|---|---|---|
| `two_J`, `parity`, `nucleus` | 2J, π | | the block's quantum numbers and compound nucleus |
| `channels()[k]` | c_k | | index into `System::channels()` of the block's k-th channel, k = 0 … n_c − 1 |
| `B()[k]` | B_c | | boundary value of channel c_k in this block (dimensionless) |
| `E()[λ]` | E_λ | MeV | level energy, as excitation energy E_x of the compound nucleus, λ = 0 … n_λ − 1 |
| `gamma()[λ n_c + k]` | γ_λc | MeV^½ | reduced width amplitude of level λ in channel c_k |

The boundary value is stored with the block, not the channel, because the engine fixes it per J^π group:
B_c = S_c(E_1), the shift function at the energy of the group's first level (`CNuc::CalcBoundaryConditions`); a
block without levels has B_c = 0. γ is one row-major array: the amplitudes of one level to all channels of the block
are adjacent. A block may have no levels: then its channels scatter through their hard-sphere phases only.

What E_λ and γ_λc are depends on the parameter basis of the project, stated once by `System::basis()`:

- `Formal`: the formal R-matrix parameters with the blocks' constant B_c (Lane and Thomas);
- `Brune`: the alternative parameters of C. R. Brune, Phys. Rev. C 66, 044611 (2002), in which E_λ are the observed
  level energies (the default of AZURE2);
- `Park`: the observed-amplitude parameters of Park, Phys. Rev. C 104, 064612 (2021), on Brune's energies.

The numbers are the project's parameters after the engine's input transformation (`CNuc::TransformIn`): partial
widths Γ (eV) and ANCs (fm^−½) of the file are converted to amplitudes before the System is built. For photon and β
channels γ is in the engine's own normalization of those channels (for β channels, the amplitude B of the file). Converting between bases and to the one internal convention is the job of
`level_parameters` (`ARCHITECTURE.md` §5; `level_parameters.md`), which the System does not call: it records the basis
in its `Convention` record (`System::convention()`, `Convention.h`) and stores the parameters unconverted.

## 4. Index conventions

- All indices in `physics/` are 0-based. `CNuc` counts pairs, J-groups, levels and channels from 1.
- Spins are stored doubled, as integers (2j, 2s, 2J), so that half-integers are exact and angular-coupling code
  needs no floating-point comparison. The adapter refuses a spin that is not a multiple of ½.
- Partitions keep the order of `CNuc`'s pair list (the order of first appearance of each pair key in `<levels>`),
  blocks the order of its J-groups, block levels and block channels the order within each J-group. System channels
  are numbered in order of first appearance, block by block.
- No object holds a pointer to another: channels name their partition by index, blocks name their channels by index.
  System and JpiBlock are values; a copy is independent, and nothing can be changed after construction. The
  constructors check every index, refuse a channel whose kind does not match its partition, a block channel of
  another compound nucleus, two blocks with the same (J, π, nucleus) and a γ array of the wrong size.

## 5. Not included yet

- external-capture final states that are not R-matrix levels (`ALevel::IsInRMatrix()` false) and J^π groups that
  exist only to hold them; their ANCs and the external reduced widths;
- spin–orbit (α, l, j) channels of the hybrid-potential model, and the nuclear potential itself;
- the THM binding energy of a Trojan-horse partition (it belongs to the THM vertex, `THM.md`);
- what belongs to the fit or the data: fixed/free flags, limits, entrance pairs, decays and reaction pathways.

## 6. Checks

`physics_system_unit` (`physics/tests/system_unit_test.cpp`) checks the storage conventions and every refusal of the
constructors using `physics/` alone. `physics_standalone` configures, builds and tests `physics/` as a project of its
own, where nothing else of AZURE2 exists: an include of an engine header (`CNuc.h` was tried) fails it at
compilation. A header on the system include path compiles there (`<gsl/gsl_math.h>` did), so `repo_consistency`
refuses in `physics/` every include but `azure2/physics/` headers and the C++ standard library.
`physics_system_from_cnuc` (`tests/upgrade/system_from_cnuc_test.cpp`) loads T1, T2 and T3 in the Brune, Park and
formal bases, builds the System from `CNuc` (`src/CNucSystem.cpp`) and compares every number with the `CNuc`
accessor it came from, with exact equality; it compares the pair order, the channel count and each block's 2J, π,
channels and (Brune, Park) level energies with the `<levels>` lines of the files. A wrong accessor, an off-by-one
between the 1-based and 0-based indices, a transposed γ layout or a level attached to the wrong block fails it.
On the `feature/upgrade-1-physics-system` build (GCC 11.5) it compares 2540 numbers, all equal; each of five
deliberate errors in the adapter (one radius for all channels, e_0 and e_1 swapped, γ stored column-major, B_c of
the first channel for all, one J-group dropped) fails it with 12 to 194 mismatches.

Code: `physics/include/azure2/physics/{System,JpiBlock,Convention}.h`, `physics/src/`, adapter `include/CNucSystem.h`,
`src/CNucSystem.cpp`.
