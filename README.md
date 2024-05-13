# DoNOFsw
Donostia Natural Orbital Functional Software

DoNOF is a computational chemistry software program that stands for Donostia Natural Orbital Functional. The original code started on January 21, 2009 as PNOFID. It will run on essentially any machine with a FORTRAN 90 compiler for 64 bit processing.

DoNOF can perform computational chemistry calculations based on the Natural Orbital Functional Theory (NOFT), including PNOF5, PNOF6, PNOF7 and GNOF. Correlation corrections after PNOF calculations can be estimated by second order perturbation theories. The total spin is conserved, not just the spin projection.

The solution is established optimizing the energy functional with respect to the occupation numbers (ONs) and to the natural orbitals (NOs), separately. The constrained nonlinear programming problem for the ONs is treated under pairing restrictions as an unconstrained minimization, while the orbital optimization is carried out by a self-consistent procedure which yields the NOs automatically orthogonal. To achieve convergence, the direct inversion of the iterative subspace (DIIS) extrapolation technique is used, and a variable scale factor balances the symmetric matrix subject to the iterative diagonalizations.

You can contact us by e-mail to DoNOFsw@gmail.com

## Installation

**Requisites.** You need a FORTRAN compiler, either gfortran or ifort. Optionally, you may want to install OpenMPI for parallel execution.

Clone the code with
~~~
git clone https://github.com/DoNOF/DoNOFsw
~~~

Go inside the DoNOFsw folder (`cd DoNOFsw`) and compile with `make [option]`. For example:
~~~
make serialg # gfortran serial -> /exe/DoNOFg.x
~~~

Other options are:
~~~
make mpig    # gfortran mpi    -> /exe/DoNOFmpig.x
make mping   # gfortran mpi (for recent linux versions)
make serial  # ifort serial    -> /exe/DoNOF.x
make mpi     # ifort mpi       -> /exe/DoNOFmpi.x
~~~

The executable will be placed inside /exe

## Execution

Several input files can be found inside /examples. A basic single point calculation with the GNOF functional looks like the following:
~~~
 &INPRUN RUNTYP='ENERGY' MULT=1 ICHARG=0 ERITYP='FULL' /
 $DATA
 Water (H2O)
 cc-pVDZ
O  8.0  0.0000     0.0000    0.1173
H  1.0  0.0000     0.7572   -0.4692
H  1.0  0.0000    -0.7572   -0.4692
 $END
 &NOFINP IPNOF=8 /
~~~

If the input is placed in a file called filename.inp, it can be executed with
~~~
./run_donofg < filename     # gfortran serial
~~~
the output will be placed in filename.out.

Other options for execution are:
~~~
./run_donofmpig < filename  # gfortran mpi
./run_donof     < filename  # ifort serial
./run_donofmpi  < filename  # ifort mpi
~~~

## Capabilities

The &INPRUN and &NOFINP namelists specify the input and output, and the fundamental job options.

The functional is controlled through **IPNOF=N** in &NOFINP, with N the number of the functional. For example, INPOF=7 indicates to use PNOF7. GNOF is indicated with IPNOF=8.

Current capabilities include:
- **IRUNTYP = 1** - Single-point Energy (Default)
- **IRUNTYP = 2** - Energy + Gradients with respect to nuclear coord
- **IRUNTYP = 3** - Geometry Optimization
- **IRUNTYP = 4** - Numerical Hessian
- **IRUNTYP = 5** - Born-Oppenheimer on-the-fly molecular dynamics

Other common options include excited states calculation (**ERPA=T**) and NOF-MBPT calculations (**MBPT=T**).

A complete list of the variables can be found in the online manual through the links: https://donof-documentation.readthedocs.io/, https://donof.readthedocs.io/
