!======================================================================!
!                                                                      !
!                               Do N O F                               !
!                                                                      !
!       (Donostia Natural Orbital Functional Software Program)         !
!                                                                      !
!                     COPYRIGHT by Mario Piris                         !
!                                                                      !
!           Donostia International Physics Center (DIPC)               !
!            University of the Basque Country (UPV/EHU)                !
!            Basque Foundation for Science (IKERBASQUE)                !
!                                                                      !
!               GNU General Public License version 3                   !
!                                                                      !
! ==================================================================== !
!                                                                      !
!      Please inform me of any bugs, by phone at: +34 943 01 8328,     !
!             by e-mail to: mario.piris@ehu.eus, mpiris@gmail.com      !
!    or write me at: Donostia International Physics Center (DIPC),     !
!                   Manuel de Lardizabal 4, 20018 Donostia, Spain.     !
!                                                                      !
! ==================================================================== !
!                                                                      !
!                        Date: December 2020                           !
!                                                                      !
!            04/26/13 MPI module developed by Eduard Matito            !
!                                                                      !
!======================================================================*

#ifdef MPI
!
!     MPI MODULE
!
      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)
      INTEGER MY_ID, NPROCS, IERR, MASTER
      LOGICAL LMASTR
      PARAMETER (MASTER=0)
      COMMON /MPIP/ STATUS,MY_ID,NPROCS,IERR,LMASTR
#else
!
!     SERIAL MODULE
!
      LOGICAL LMASTR
      INTEGER NPROCS,MASTER,MY_ID
      COMMON /SERIAL/ MY_ID,NPROCS,LMASTR
      PARAMETER (MASTER=1)
#endif
