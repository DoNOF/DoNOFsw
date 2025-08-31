# ![Donostia Natural Orbital Functional Software](https://github.com/DoNOF/DoNOF-Documentation/blob/master/docs/Logo-DoNOF.jpeg)

|📫 **Contact us:**    | DoNOFsw@gmail.com                         |
|:------------------|-------------------------------------------|
|📖 **Documentation:** | https://donof-documentation.readthedocs.io|

**DoNOF** is a computational chemistry software program that stands for **Donostia Natural Orbital Functionals**. The original code started on January 21, 2009 as PNOFID. It will run on essentially any machine with a FORTRAN 90 compiler for 64 bit processing.

DoNOF can perform computational chemistry calculations based on the Natural Orbital Functional Theory (NOFT), including PNOF5, PNOF6, PNOF7 and GNOF. Correlation corrections after PNOF calculations can be estimated by second order perturbation theories. The total spin is conserved, not just the spin projection.

The solution is established optimizing the energy functional with respect to the occupation numbers (ONs) and to the natural orbitals (NOs), separately. The constrained nonlinear programming problem for the ONs is treated under pairing restrictions as an unconstrained minimization, while the orbital optimization is carried out by a self-consistent procedure which yields the NOs automatically orthogonal. To achieve convergence, the direct inversion of the iterative subspace (DIIS) extrapolation technique is used, and a variable scale factor balances the symmetric matrix subject to the iterative diagonalizations.

## Installation

**Requisites.** You need a FORTRAN compiler, either gfortran or ifort. Optionally, you may want to install OpenMPI for parallel execution.

0. Install ![libcint](https://github.com/sunqm/libcint). The following instructions are provided as example.
~~~
git clone http://github.com/sunqm/libcint.git
cd libcint
mkdir build; cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr/local/lib ..
sudo make install
~~~

1. Clone the code with
~~~
git clone https://github.com/DoNOF/DoNOFsw
~~~

2. Go inside the DoNOFsw folder (`cd DoNOFsw`) and compile with `make [option]`. For example:
~~~
make serialg # gfortran serial -> /exe/DoNOFg.x
~~~

The possible options are:
~~~
make serialg # gfortran Serial      -> /exe/DoNOFg.x
make ompg    # gfortran OpenMP      -> /exe/DoNOFompg.x
make mpig    # gfortran MPI         -> /exe/DoNOFmpig.x
make hybridg # gfortran OpenMP+MPI  -> /exe/DoNOFhybridg.x
make serial  # ifort Serial         -> /exe/DoNOF.x
make omp     # ifort OpenMP         -> /exe/DoNOFomp.x
make mpi     # ifort MPI            -> /exe/DoNOFmpi.x
make hybrid  # ifort OpenMP+MPI     -> /exe/DoNOFhybrid.x
~~~

3. The executable will be placed inside /exe

## Execution

Several input files can be found inside /examples. A basic single point calculation with the GNOF functional looks like the following:
~~~
 &INPRUN RUNTYP='ENERGY' MULT=1 ICHARG=0 USELIB=T ERITYP='FULL' /
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
./run_donofg filename     # gfortran serial
~~~
the output will be placed in filename.out.

The possible options are:
~~~
./run_donofg filename     # gfortran serial
./run_donofompg filename  # gfortran omp
./run_donofmpig filename  # gfortran mpi
./run_donof filename     # ifort serial
./run_donofomp filename  # ifort omp
./run_donofmpi filename  # ifort mpi
~~~

## Capabilities

The &INPRUN and &NOFINP namelists specify the input and output, and the fundamental job options.

The functional is controlled through **IPNOF=N** in &NOFINP, with N the number of the functional. For example:
**GNOFm**: IPNOF=8 lmod=1
**GNOF**: IPNOF=8 lmod=0
**PNOF7**: IPNOF=7
**PNOF6**: IPNOF=6
**PNOF5**: IPNOF=5

Current capabilities include:
- **RUNTYP = ENERGY** - Single-point Energy (Default)
- **RUNTYP = GRAD** - Energy + Gradients with respect to nuclear coord
- **RUNTYP = OPTGEO** - Geometry Optimization
- **RUNTYP = HESS** - Numerical Hessian
- **RUNTYP = DYN** - Born-Oppenheimer on-the-fly molecular dynamics

Other common options include excited states calculation (**ERPA=T**) and NOF-MBPT calculations (**MBPT=T**).
