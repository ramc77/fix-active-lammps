# fix\_active: Ten Active-Matter Models in One LAMMPS Fix

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![LAMMPS](https://img.shields.io/badge/LAMMPS-2023%2B-green.svg)](https://www.lammps.org/)

A single LAMMPS fix — `fix active` — that supplies ten distinct active-matter
self-propulsion kernels behind a uniform input-script interface. Switching
between models is a one-line change in your input deck.

## Features

- **Ten distinct kernels.** ABP, AOUP, Chiral ABP, Run-and-Tumble, Vicsek
  alignment, Active Rods, Spinners with odd viscosity, Inertial ABP,
  Motility-Induced Phase Separation, Active Nematics — each is its own
  `compute_<model>()` routine; no kernel re-uses another.
- **Single-command interface.** No multi-fix coordination, no source-tree
  modifications beyond dropping in two files.
- **MPI-correct.** Per-atom orientation, AOUP active velocity, and
  quaternion scaffold all follow the atom across MPI exchange, atom-sort
  callbacks, and ghost forward-communication. The Vicsek and nematic
  alignment kernels work under arbitrary spatial decomposition.
- **Validated.** Quantitative agreement within a few percent against the
  analytic predictions for every overdamped kernel; collective benchmarks
  reproduce the literature signatures of MIPS, the Vicsek transition, and
  active turbulence.

## Supported models

| Keyword | Model | Combine with |
|---|---|---|
| `abp`     | Active Brownian Particles                | `fix nve + fix langevin` |
| `aoup`    | Active Ornstein–Uhlenbeck Particles      | `fix nve + fix langevin` |
| `chiral`  | Chiral ABP (intrinsic angular velocity)  | `fix nve + fix langevin` |
| `rtp`     | Run-and-Tumble Particles                 | `fix nve + fix langevin` |
| `vicsek`  | Vicsek-type polar alignment              | `fix nve + fix langevin` |
| `rod`     | Active rod with anisotropic friction     | `fix nve` only (self-contained Langevin) |
| `spinner` | Chiral spinner with odd (Hall) viscosity | `fix nve + fix langevin` |
| `iabp`    | Inertial ABP                             | `fix nve` only (self-contained Langevin) |
| `mips`    | Motility-Induced Phase Separation (Cates–Tailleur quorum sensing) | `fix nve + fix langevin + pair_style` |
| `nematic` | Active nematic with apolar Q-tensor torque | `fix nve + fix langevin + pair_style` |

> **Scope.** This release applies the active force in the `xy`-plane.
> A 3D box is accepted but the active push has no `z` component (a
> warning is emitted at `init`); per-atom quaternion storage is in place
> for a future fully three-dimensional update.

## Installation

```bash
# 1. Get the source
git clone https://github.com/ramc77/fix-active-lammps.git
cd fix-active-lammps

# 2. Drop into your LAMMPS tree
cp src/fix_active.cpp src/fix_active.h $LAMMPS_ROOT/src/

# 3. Build (CMake)
cd $LAMMPS_ROOT
mkdir -p build && cd build
cmake -D BUILD_MPI=ON ../cmake
make -j8 lmp

# ... or the traditional Makefile build:
# cd $LAMMPS_ROOT/src && make mpi -j8
```

No additional LAMMPS package is required.

## Quick start

```lammps
# Minimal 2D ABP: 1000 particles, no pair interactions, isolated ensemble
units       lj
dimension   2
boundary    p p p
atom_style  atomic
atom_modify map array

region      box block 0 1000 0 1000 -0.5 0.5
create_box  1 box
create_atoms 1 random 1000 12345 box
mass        1 1.0
velocity    all set 0.0 0.0 0.0
pair_style  none

fix nve     all nve
fix lang    all langevin 0.01 0.01 1.0 11111         # provides Stokes drag + kT
fix active1 all active abp v0 1.0 Dr 0.1 gamma 1.0 seed 54321
fix enf2d   all enforce2d

compute     msd all msd com no
thermo_style custom step temp c_msd[4]
thermo      5000

timestep    0.001
run         1000000
```

Long-time MSD slope ≈ 20.04 (= 4·D\_eff with D\_eff = D\_t + v₀²/(2·D\_r)).

## Command syntax

```lammps
fix ID group-ID active <model> keyword value [keyword value ...]
```

- `ID`, `group-ID`: standard LAMMPS fields.
- `<model>`: one of the ten keywords above.
- Keyword/value pairs are order-independent. Each model consumes the
  keywords listed in [`doc/fix_active_manual.pdf`](doc/fix_active_manual.pdf);
  keywords that are not consumed by the chosen model are silently ignored.

### Selected examples

ABP — see Quick start above.

Vicsek (ordered regime):
```lammps
fix active1 all active vicsek v0 1.0 Dr 0.1 Rcut 1.0 alignment 1.0 gamma 1.0 seed 54321
```
The `alignment` parameter is in [0, 1] and interpolates between no
alignment and full Vicsek mean.

MIPS (Cates–Tailleur quorum sensing) — `Rcut` and `rho_max` are required:
```lammps
fix active1 all active mips v0 30.0 Dr 0.05 Rcut 1.5 rho_max 1.0 gamma 1.0 seed 54321
```

Spinner Hall benchmark (drag comes from `fix langevin`):
```lammps
fix nve_int  spin nve
fix drag     spin langevin 0.0 0.0 1.0 99999       # Tdamp=1 -> gamma=1
fix drive    spin addforce 1.0 0.0 0.0
fix active1  spin active spinner omega_spin 1.0 eta_odd 1.0 seed 54321
```
Steady state: ⟨v\_x⟩ = ⟨v\_y⟩ = γ·F\_ext/(γ² + (η\_odd·ω)²) = 0.5,
tan φ = η\_odd·ω/γ = 1.

For the full per-model documentation, equations of motion, and worked
examples, see [`doc/fix_active_manual.pdf`](doc/fix_active_manual.pdf).

## Repository layout

```
fix-active-lammps/
├── README.md                 # this file
├── LICENSE                   # GPL v3
├── src/
│   ├── fix_active.cpp        # LAMMPS fix implementation
│   └── fix_active.h
├── doc/
│   ├── fix_active_manual.tex # LaTeX source of the user manual
│   └── fix_active_manual.pdf # compiled manual / tutorial
└── examples/                 # ready-to-run input scripts
```

## Validation summary

All errors are measured against the closed-form analytic prediction over
the diffusive interval t ∈ [5τ\_r, 50τ\_r], well past the
ballistic-to-diffusive crossover but before finite-N noise dominates the
late-time tail of a single N = 10³ ensemble.

| Model              | Observable                | Error |
|--------------------|---------------------------|-------|
| ABP (2D)           | D\_eff                    | 1.6% |
| AOUP               | D\_eff                    | 2.8% |
| Chiral ABP         | D\_eff (Eq. for chiral)   | 2.5% |
| RTP                | D\_eff, ⟨run time⟩        | 5.1% |
| IABP               | D\_eff (overdamped limit) | 0.5% |
| Active rod         | Trajectory + persistence  | Qualitative |
| Vicsek (no align)  | D\_eff (ABP limit)        | 2.1% |
| Vicsek             | Φ vs D\_r                 | Qualitative (literature) |
| Spinner            | Hall angle tan φ          | < 0.1% |
| MIPS               | Phase separation          | Qualitative (bimodal coordination) |
| Nematic            | Active turbulence ⟨S⟩     | Qualitative (S ≪ 1) |

## Parallel computing

The per-atom orientation θ is forward-communicated to ghost atoms before
the Vicsek and nematic kernels iterate the neighbour list (declared via
`comm_forward = 1` plus overridden `pack_forward_comm` /
`unpack_forward_comm`). Without this step alignment partners across MPI
subdomain boundaries — and across periodic-image boundaries even in
serial runs — are silently dropped. Always build with MPI and run with
`mpirun -np N lmp -in ...` for production work involving the Vicsek or
nematic kernels.

The kernels are not OpenMP-threaded. For hybrid runs set
`OMP_NUM_THREADS=1` and parallelise across MPI ranks only.

## Documentation

- **[doc/fix\_active\_manual.pdf](doc/fix_active_manual.pdf)** —
  complete user manual: installation, syntax, per-model documentation,
  three worked examples, performance, troubleshooting, FAQ.

## Citation

If you use this code in published work, please cite:

```bibtex
@misc{ChandBukhari2026_fixactive,
  author       = {Ram Chand and Saqia Bukhari},
  title        = {Active Matter Model Implementation in {LAMMPS}: A Coherent {\texttt{fix active}} Framework of Ten Active Self-Propulsion Mechanisms},
  year         = {2026},
  note         = {Manuscript},
  howpublished = {\url{https://github.com/ramc77/fix-active-lammps}}
}
```

And, for the foundational LAMMPS active-matter implementation:

```bibtex
@misc{Dias2021,
  author = {C. S. Dias},
  title  = {Molecular Dynamics Simulations of Active Matter using LAMMPS},
  year   = {2021},
  note   = {\url{https://arxiv.org/abs/2102.10399}}
}
```

## Requirements

- LAMMPS (2023 or later recommended)
- C++17 compiler (g++ ≥ 9, clang++ ≥ 10, or AppleClang ≥ 13)
- MPI library (Open MPI / MPICH) for parallel runs
- Python 3.9+ with NumPy / Matplotlib for the analysis scripts

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).

## Contributing

Bug reports, pull requests, and additional active-matter kernels are
welcome. To contribute a new kernel:

1. Add a new `compute_<your_model>()` routine to `fix_active.cpp` and the
   matching enum entry / dispatch case.
2. Register any new model-specific keyword in the parser block.
3. Add a worked example under `examples/` and document the model in
   `doc/fix_active_manual.tex`.
4. Open a pull request.

## Contact

- Issues: please open a GitHub issue.
- Email: ram.chand2k11@yahoo.com

## Acknowledgments

LAMMPS development team at Sandia National Laboratories; C. S. Dias for
the foundational LAMMPS active-matter implementation; the active-matter
research community for the underlying theory.
