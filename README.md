# fix_active: Comprehensive Active Matter Simulations in LAMMPS

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![LAMMPS](https://img.shields.io/badge/LAMMPS-2023+-green.svg)](https://www.lammps.org/)

Active Matter Model Implementation in LAMMPS: A Coherent fix active Framework of Ten Active Self-Propulsion Mechanisms

## Features

- **10 Active Matter Models** in a single, unified fix command
- **No source code modification** required - just compile and use
- **Fully validated** against analytical predictions
- **Parallel-ready** - works with MPI for large-scale simulations
- **Easy model switching** via input script keywords

## Supported Models

| Model | Command | Description |
|-------|---------|-------------|
| ABP | `abp` | Active Brownian Particles |
| AOUP | `aoup` | Active Ornstein-Uhlenbeck Particles |
| Chiral | `chiral` | Circle swimmers |
| RTP | `rtp` | Run-and-Tumble (bacterial motion) |
| Vicsek | `vicsek` | Collective alignment/flocking |
| Rod | `rod` | Anisotropic active particles |
| Spinner | `spinner` | Odd viscosity / Hall effect |
| IABP | `iabp` | Inertial Active Brownian Particles |
| MIPS | `mips` | Motility-Induced Phase Separation |
| Nematic | `nematic` | Active liquid crystals |

## Quick Installation

```bash
# Clone repository
git clone https://github.com/[username]/fix_active_lammps.git
cd fix_active_lammps

# Copy source files to LAMMPS
cp src/fix_active.* /path/to/lammps/src/

# Compile LAMMPS
cd /path/to/lammps/src
make mpi
```

## Quick Start

```lammps
# Minimal ABP simulation
units lj
dimension 2
atom_style sphere
boundary p p p

region box block 0 50 0 50 -0.5 0.5
create_box 1 box
create_atoms 1 random 1000 12345 box

mass 1 1.0
set type 1 diameter 1.0

pair_style soft 1.12246
pair_coeff 1 1 100.0

# Active Brownian Particles
fix 1 all active abp v0 1.0 Dr 0.1 seed 54321
fix 2 all nve

timestep 0.001
run 100000
```

## Command Syntax

```lammps
fix ID group-ID active model [model-specific parameters]
```

### ABP Example
```lammps
fix active all active abp v0 1.0 Dr 0.1 Dt 0.01 seed 12345
```

### MIPS Example (high activity)
```lammps
fix active all active mips v0 40.0 Dr 3.0 seed 12345
```

### Vicsek Example (flocking)
```lammps
fix active all active vicsek v0 0.5 Dr 0.1 Rcut 1.5 alignment 5.0 seed 12345
```

## Repository Structure

```
fix_active_lammps/
├── README.md                 # This file
├── LICENSE                   # GPL v3 license
├── src/                      # Source code
│   ├── fix_active.cpp
│   └── fix_active.h
├── docs/                     # Documentation
│   └── USER_MANUAL.md        # Detailed user guide
├── examples/                 # Example input scripts
│   ├── abp_basic.in
│   ├── mips_phase_separation.in
│   ├── vicsek_flocking.in
│   └── ...
├── validation/               # Validation scripts and data
│   ├── abp_msd_validation.py
│   ├── mips_validation.py
│   └── ...
└── scripts/                  # Utility scripts
    └── plot_results.py
```

## Validation

All models are validated against analytical predictions:

| Model | Test | Error |
|-------|------|-------|
| ABP (2D) | Effective diffusion | < 3% |
| ABP (3D) | Effective diffusion | < 4% |
| AOUP | Effective diffusion | < 5% |
| Chiral | Reduced diffusion | < 8% |
| RTP | Run-time distribution | < 3% |
| IABP | Overdamped limit | < 2% |
| Vicsek | Order parameter | ~ 3% |
| Spinner | Hall angle | < 2% |
| MIPS | Phase separation | Qualitative |
| Nematic | Active turbulence | Qualitative |

## Documentation

- **[User Manual](docs/USER_MANUAL.md)** - Complete parameter reference and examples
- **[Examples](examples/)** - Ready-to-run input scripts
- **[Validation](validation/)** - Scripts to reproduce validation results

## Citation

If you use this code, please cite:

```bibtex
@article{fix_active2026,
  title={CActive Matter Model Implementation in LAMMPS: A Coherent \texttt{fix active}  Framework of Ten Active Self-Propulsion Mechanisms},
  author={Ram Chand},
  journal={Physical Review E},
  year={2026},
  note={In preparation}
}
```

And the foundational LAMMPS active matter work:

```bibtex
@article{dias2021molecular,
  title={Molecular Dynamics Simulations of Active Matter using LAMMPS},
  author={Dias, C. S.},
  journal={arXiv preprint arXiv:2102.10399},
  year={2021}
}
```

## Requirements

- LAMMPS (2023 or later recommended)
- C++ compiler (g++ 7.0+ or clang++ 10.0+)
- MPI library (optional, for parallel runs)
- Python 3.7+ with NumPy, Matplotlib (for analysis scripts)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-model`)
3. Commit changes (`git commit -am 'Add new model'`)
4. Push to branch (`git push origin feature/new-model`)
5. Open a Pull Request

## Contact

- **Issues**: Please open a GitHub issue for bugs or feature requests
- **Email**: ram.chand2k11@yahoo.com

## Acknowledgments

- LAMMPS development team at Sandia National Laboratories
- C. S. Dias for foundational LAMMPS active matter implementation
- Active matter research community for theoretical foundations
