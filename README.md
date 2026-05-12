# DoNOF

DoNOF is a computational chemistry code based on Natural Orbital Functional Theory (NOFT). The original code started on January 21, 2009 as PNOFID.

Current functionality includes:
- PNOF5, PNOF6, PNOF7, and GNOF
- Single-point energies, gradients, geometry optimization, numerical Hessians, transition-state optimization, and molecular dynamics
- ERPA, post-NOF corrections, and nonlinear optical properties

## Build

Requirements:
- A Fortran compiler compatible with the current `Makefile`
- `libcint`
- OpenMPI, if MPI execution is needed

Library source:
- <https://github.com/sunqm/libcint>

Build targets:

```bash
make serialg   # exe/DoNOFg.x
make ompg      # exe/DoNOFompg.x
make mpig      # exe/DoNOFmpig.x
```

## Run

Example inputs are available in `examples/`.

Available wrapper scripts:

```bash
./run_donofg filename
./run_donofompg filename
./run_donofmpig filename
./run_donofg_dyn filename
./run_donofompg_dyn filename
./run_donofmpig_dyn filename
```

Input is controlled mainly through `&INPRUN` and `&NOFINP`. Molecular-dynamics jobs also use `&INPDYN`.

## Documentation

- Online documentation: <https://donof-documentation.readthedocs.io>
- Local namelist reference: [doc/namelist.md](doc/namelist.md)
- Local PDF version: [doc/namelist.pdf](doc/namelist.pdf)

## Contact

`DoNOFsw@gmail.com`
