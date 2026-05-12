# DoNOF Input Reference

This document summarizes the current DoNOF input namelists and is intended as a practical reference for users preparing input files.

## Overview

DoNOF input is organized around three namelists:

- `&INPRUN`: run type, integral engine, ERI/RI options, electric field, and general limits
- `&NOFINP`: NOF settings, orbital optimization, perturbative options, restart behavior, and output controls
- `&INPDYN`: molecular-dynamics controls, used only when `RUNTYP='DYN'`

Most jobs use `&INPRUN`, `$DATA`, and `&NOFINP`. Dynamics jobs also use `&INPDYN`.

## Minimal Examples

### Single-point energy

```text
 &INPRUN RUNTYP='ENERGY' MULT=1 ICHARG=0 /
 $DATA
 Water
 aug-cc-pVDZ
 O 8.0  0.000000  0.000000  0.000000
 H 1.0  0.000000  0.757000  0.587000
 H 1.0  0.000000 -0.757000  0.587000
 $END
 &NOFINP IPNOF=8 /
```

### Geometry optimization

```text
 &INPRUN RUNTYP='OPTGEO' /
 $DATA
 Water optimization
 aug-cc-pVDZ
 O 8.0  0.000000  0.000000  0.000000
 H 1.0  0.000000  0.757000  0.587000
 H 1.0  0.000000 -0.757000  0.587000
 $END
 &NOFINP IPNOF=8 MOLDENGEO=1 /
```

### Transition-state optimization

```text
 &INPRUN RUNTYP='TSOPT' /
 $DATA
 Transition state guess
 aug-cc-pVDZ
 ...
 $END
 &NOFINP IPNOF=8 MOLDENGEO=1 /
```

### Molecular dynamics

```text
 &INPRUN RUNTYP='DYN' MULT=1 ICHARG=0 /
 $DATA
 Dynamics example
 cc-pVDZ
 F 9.0  -0.3000  -6.0000   0.0000
 H 1.0  -0.3707   0.0000   0.0000
 H 1.0   0.3707   0.0000   0.0000
 $END
 &NOFINP IPNOF=8 Imod=1 /
 &INPDYN dt=0.1 tmax=50.0 Vxyz=0,0.1,0,0.025,0,0,-0.025 /
```

Dynamics runs require a valid `GCF` file to initialize occupations and orbitals.

### ERPA calculation

```text
 &INPRUN RUNTYP='ENERGY' /
 $DATA
 ERPA example
 aug-cc-pVDZ
 ...
 $END
 &NOFINP IPNOF=8 ERPA=.TRUE. /
```

### Nonlinear optical properties

```text
 &INPRUN RUNTYP='ENERGY' NLOP=1 NPOINT=9 STEP=1.0d-4 /
 $DATA
 NLOP example
 aug-cc-pVDZ
 ...
 $END
 &NOFINP IPNOF=8 /
```

## `&INPRUN`

### Main job controls

- `RUNTYP='ENERGY'`
  - `ENERGY`: single-point energy
  - `GRAD`: analytic nuclear gradient
  - `OPTGEO`: geometry optimization
  - `HESS`: numerical Hessian from analytic gradients
  - `TSOPT`: first-order saddle-point search
  - `DYN`: Born-Oppenheimer on-the-fly molecular dynamics
- `MULT=1`: spin multiplicity
- `ICHARG=0`: molecular charge
- `IECP=0`
  - `0`: all-electron calculation
  - `1`: read ECP from `$ECP`
- `IEMOM=1`
  - `1`: dipole moment
  - `2`: dipole and quadrupole
  - `3`: dipole, quadrupole, and octopole

Current defaults include `ERITYP='RI'`, `USELIB=.TRUE.`, `GTYP='SPH'`, and `IPNOF=8`.

### Nonlinear optical properties

- `NLOP=0`
  - `-1`: compute `alpha`, `beta`, and `gamma`
  - `0`: disabled
  - `1`: polarizability `alpha`
  - `2`: first hyperpolarizability `beta`
  - `3`: second hyperpolarizability `gamma`
- `NPOINT=9`: number of field values used in the dyadic Romberg-Richardson procedure
- `STEP=1.0d-4`: initial electric-field step
- `ISOALPHA=0`: request isotropic and anisotropic polarizability analysis

`NLOP /= 0` is only valid with `RUNTYP='ENERGY'`.

### Coordinates and electric field

- `UNITS='ANGS'`
  - `ANGS`
  - `BOHR`
- `EVEC=0.0d0,0.0d0,0.0d0`: electric field in a.u.

### Integral engine and basis conventions

- `USELIB=.TRUE.`
  - `.TRUE.`: LIBCINT
  - `.FALSE.`: HONDO
- `GTYP='SPH'`
  - `SPH`: spherical Gaussians, valid only with `USELIB=.TRUE.`
  - `CART`: Cartesian Gaussians
- `USEHUB=.FALSE.`: Hubbard-model mode

### ERI and RI controls

- `DONTW=.TRUE.`: do not store two-electron integrals on disk
- `ERITYP='RI'`
  - `FULL`: four-center ERIs
  - `RI`: RI approximation
  - `MIX`: RI first, then full ERIs after convergence
- `CUTOFF=1.0d-9`: Schwarz screening threshold
- `RITYP='JKFIT'`
  - `JKFIT`
  - `GEN`
  - `RIFIT`
- `GEN='A2*'`: generative auxiliary basis if `RITYP='GEN'`
- `SMCD=.FALSE.`: symmetric modified Cholesky decomposition for the RI `G` matrix

### Hessian and thermochemistry

- `HSSCAL=.TRUE.`: in `OPTGEO`, compute a Hessian from analytic gradients and perform vibrational analysis
- `PROJECT=.TRUE.`: project out rotational and translational contaminants from the Hessian
- `ISIGMA=1`: rotational symmetry number

### Allocation limits

- `NATmax=100`
- `NSHELLmax=500`
- `NPRIMImax=2000`
- `NSHELLAUXmax=1000`
- `NPRIMIAUXmax=1000`

## `&NOFINP`

### Global optimization

- `MAXIT=1000`: maximum OCC-SCF iterations
- `MAXIT21=3`: initial NO-only stage when `ICOEF=21`
- `ICOEF=1`
  - `0`: optimize occupations only
  - `1`: optimize occupations and orbitals
  - `2`: optimize orbitals only
  - `21`: first run a short orbital-only stage, then switch to `ICOEF=1`
  - `3`: fragment-oriented mode
- `ISOFTMAX=1`
  - `0`: trigonometric occupation parametrization
  - `1`: softmax occupation parametrization
- `IORBOPT=2`
  - `1`: iterative diagonalization
  - `2`: ADAM
  - `3`: AdaBelief
  - `4`: YOGI
  - `5`: DEMON
  - `6`: SQP
- `IEINI=0`: stop after the initial energy evaluation
- `NO1=0`
  - `-1`: include core orbitals as fully occupied
  - `0`: all orbitals considered
  - positive integer: user-defined number of fully occupied orbitals

### Hartree-Fock guess

- `IRHF=1`
  - `0`: do not obtain HF orbitals
  - `1`: HF SCF
  - `2`: ADAM-based orbital optimization
  - `3`: iterative diagonalization
- `NCONVRHF=5`
- `MAXITRHF=100`
- `HFDAMP=.TRUE.`
- `HFEXTRAP=.TRUE.`
- `HFDIIS=.TRUE.`
- `KOOPMANS=0`

`IRHF=1` is not allowed with a nonzero electric field. In `OPTGEO`, `HESS`, and `DYN`, `IRHF` may be forced to `0`.

### Functional selection

- `IPNOF=8`
  - `5`: PNOF5
  - `6`: PNOF6
  - `7`: PNOF7
  - `8`: GNOF
- `Ista=0`
  - `0`: PNOF7
  - `1`: PNOF7s
- `Imod=0`
  - `0`: GNOF
  - `1`: GNOFm
- `HighSpin=.FALSE.`
- `NCWO=-1`
  - `1,2,3,...`
  - `-1`: automatic `NVIR/NDOC`

For `IPNOF=8`, only `Imod=0` and `Imod=1` are accepted.

### Convergence thresholds

- `NTHRESHL=4`
- `NTHRESHE=8`
- `NTHRESHEC=8`
- `NTHRESHEN=8`

These correspond to:

- `THRESHL = 10 raised to (-NTHRESHL)`
- `THRESHE = 10 raised to (-NTHRESHE)`
- `THRESHEC = 10 raised to (-NTHRESHEC)`
- `THRESHEN = 10 raised to (-NTHRESHEN)`

Notes:

- in `DYN`, if `NTHRESHE` is left at the default `8`, it is increased to `10`
- `NTHRESHEC < 8` is promoted to `8`

### Iterative-diagonalization controls

- `MAXLOOP=10`
- `SCALING=.TRUE.`
- `AUTOZEROS=.TRUE.`
- `NZEROS=0`
- `NZEROSm=5`
- `NZEROSr=2`
- `ITZITER=10`
- `DIIS=.TRUE.`
- `NTHDIIS=3`
- `NDIIS=5`
- `PERDIIS=.TRUE.`
- `DAMPING=.FALSE.`
- `EXTRAP=.FALSE.`

### ADAM-family controls

- `LR=0.01d0`
- `FACT=0.2d0`
- `BETA1=0.7d0`
- `BETA2=0.9d0`

### Perturbative and response options

- `ERPA=.FALSE.`
- `OIMP2=.FALSE.`
- `MBPT=.FALSE.`
- `NO1PT2=-1`
- `SC2MCPT=.FALSE.`
- `NEX=0`

Notes:

- `APSG` and `SC2MCPT` are only valid for `IPNOF=5`
- `OIMP2` is disabled automatically in `DYN`

### Restart and initial guesses

- `RESTART=.FALSE.`
- `INPUTGAMMA=0`
- `INPUTC=0`
- `INPUTFMIUG=0`
- `INPUTCXYZ=0`

When `RESTART=.TRUE.`, the code forces:

- `INPUTGAMMA=1`
- `INPUTC=1`
- `INPUTFMIUG=1`
- `INPUTCXYZ=1`

In `DYN`:

- the code forces `RESTART=.FALSE.`
- `INPUTGAMMA=1`
- `INPUTC=1`
- `INPUTFMIUG=1` only when `IORBOPT=1`
- `INPUTCXYZ=0`

A valid `GCF` file is still required to initialize a dynamics run.

### Output controls

- `NPRINT=0`
- `IWRITEC=0`
- `IMULPOP=0`
- `PRINTLAG=.FALSE.`
- `DIAGLAG=.FALSE.`
- `IEKT=0`
- `IAIMPAC=0`
- `IFCHK=1`
- `MOLDEN=1`
- `MOLDENGEO=0`
- `INICOND=1`
- `NOUTRDM=0`
- `NTHRESHDM=6`
- `NSQT=1`
- `NOUTCJK=0`
- `NTHRESHCJK=6`
- `NOUTTijab=0`
- `NTHRESHTijab=6`
- `APSG=.FALSE.`
- `NTHAPSG=10`

`MOLDENGEO=1` saves one Molden file per geometry visited in `OPTGEO` or `TSOPT` as `snapshot-####.mld`. It also forces `MOLDEN=1`.

In `DYN`:

- `IAIMPAC=0`
- `MOLDEN=0`
- `MOLDENGEO=0`
- `INICOND=0`
- `IFCHK=0`

### Miscellaneous controls

- `ORTHO=.TRUE.`
- `CHKORTHO=.FALSE.`
- `FROZEN=.FALSE.`
- `IFROZEN=0`
- `ICGMETHOD=1`
  - `1`: SUMSL / CG
  - `2`: NAG-based path, available only if explicitly enabled
  - `3`: LBFGS

### Excited-state ensemble controls

- `NESt=0`
- `OMEGA1=1.0d0`

## `&INPDYN`

This namelist is used only when `RUNTYP='DYN'`.

- `velflag='F'`: constant-velocity dynamics
- `dt=1.0d0`: time step in fs
- `tmax=100.0d0`: total propagation time in fs
- `Vxyz=0.0d0`: initial Cartesian velocities per atom
- `resflag='F'`: restart molecular dynamics
- `snapshot='F'`: save `snapshot-<time>.mld`
- `jumptol=1.0d-4`: threshold for treating a potential-energy jump as a real event
- `energybound='T'`: enable velocity rescaling and event guard
- `integrator='BEEVER'`
  - `BEEVER`
  - `VERLET`

### Notes

- `&INPDYN` is read only when `RUNTYP='DYN'`
- `resflag='T'` reads the last geometry and velocities from `DYNl.xyz`
- `snapshot='T'` writes Molden snapshots during the trajectory
- `energybound` accepts only `T` or `F`
- `integrator` accepts only `BEEVER` or `VERLET`
