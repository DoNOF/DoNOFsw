# DoNOFsw
Donostia Natural Orbital Functional Software

DoNOF is a computational chemistry software program that stands for Donostia Natural Orbital Functional. The original code started on January 21, 2009 as PNOFID. It will run on essentially any machine with a FORTRAN 90 compiler for 64 bit processing.
 
DoNOF can perform computational chemistry calculations based on the Natural Orbital Functional Theory (NOFT), including PNOF5, PNOF6 and PNOF7. Correlation corrections after PNOF calculations can be estimated by second order perturbation theories. The total spin is conserved, not just the spin projection.
 
The solution is established optimizing the energy functional with respect to the occupation numbers (ONs) and to the natural orbitals (NOs), separately. The constrained nonlinear programming problem for the ONs is treated under pairing restrictions as an unconstrained minization, while the orbital optimization is carried out by a self-consistent procedure which yields the NOs automatically orthogonal. To achieve convergence, the direct inversion of the iterative subspace (DIIS) extrapolation technique is used, and a variable scale factor balances the symmetric matrix subject to the iterative diagonalizations.
 
The &INPRUN and &NOFINP namelists specify the input and output, and the fundamental job options.

You can contact us by e-mail to DoNOFsw@gmail.com
