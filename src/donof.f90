!======================================================================!
!                                                                      !
!                               Do N O F                               !
!                                                                      !
!       (Donostia Natural Orbital Functional Software Program)         !
!                                                                      !
!                  COPYRIGHT by Mario Piris (2009)                     !
!                                                                      !
!    IPR registered under Basque Government and Spanish Ministry ECD   !
!                  Registration number 01/2020/360                     !
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
!        by e-mail to: mario.piris@ehu.eus, or write to me at:         !
!            Donostia International Physics Center (DIPC),             !
!            Manuel de Lardizabal 4, 20018 Donostia, Spain.            !
!                                                                      !
! ==================================================================== !
!                                                                      !
!                           Date: November 2024                        !
!                                                                      !
!    Program to compute the ground state properties of a molecule      !
!    in the gas phase using PNOF5 - GNOF + perturbation corrections    !
!                                                                      !
!              ( Comp. Phys. Comm. 259, 107651, 2021 )                 !
!                                                                      !
!======================================================================!

      PROGRAM DoNOF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
#include "mpip.h"                                                      
      character(8) :: date
      character(10):: time
      character(5) :: zone
      character(100):: sha
      integer,dimension(8) :: val
!
      INTEGER,ALLOCATABLE,DIMENSION(:)::IAN0,IMIN0,IMAX0,KSTART0,KATOM0
      INTEGER,ALLOCATABLE,DIMENSION(:)::IAN, IMIN, IMAX, KSTART, KATOM
      INTEGER,ALLOCATABLE,DIMENSION(:)::KTYPE0,KLOC0,INTYP0,KNG0,KMIN0
      INTEGER,ALLOCATABLE,DIMENSION(:)::KTYPE, KLOC, INTYP, KNG, KMIN
      INTEGER,ALLOCATABLE,DIMENSION(:)::KMAX0,ISH0,ITYP0
      INTEGER,ALLOCATABLE,DIMENSION(:)::KMAX, ISH, ITYP
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::Cxyz0,Cxyz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: ZAN0,ZMASS0,C10,C20
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: ZAN, ZMASS, C1, C2 
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: EX0,CS0,CP0,CD0
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: EX, CS, CP, CD
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: CF0,CG0,CH0,CI0
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: CF, CG, CH, CI
!
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: GRADS      
      DOUBLE PRECISION,DIMENSION(3) :: EVEC,DIPS
      INTEGER :: INTTYPE
      INTEGER :: CR, CM, TIMESTART, TIMEFINISH
      DOUBLE PRECISION :: RATE

      INTEGER :: SIZE_ENV, NBAS                          !LIBCINT
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: ENV0  !LIBCINT
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ATM0, BAS0   !LIBCINT
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: ENV   !LIBCINT
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ATM, BAS     !LIBCINT
!-----------------------------------------------------------------------
!     MPI initialization
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL INIMPI()
#ifdef MPI
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,MY_ID,IERR)
      IF(MY_ID > 0) GOTO 10
#endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call date_and_time(date,time,zone,val)
      write(6,1)val(5),val(6),val(2),val(3),val(1)
!-----------------------------------------------------------------------
!     Initialization for system_clock
!-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT_RATE=CR)
      CALL SYSTEM_CLOCK(COUNT_MAX=CM)
      RATE = REAL(CR)
      CALL SYSTEM_CLOCK(TIMESTART)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Header on the output file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,2)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Input namelist INPRUN variables (NINTEG=NINTMX)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!     IRUNTYP: Specifies the run calculation
!         ICH: Molecular charge  
!         MUL: Multiplicity of the electronic state
!      NINTEG: Total Number of 2e- integrals (NINTEG = NINTMX)
!      IDONTW: Do not write 2e- integrals on the disk (Unit=1)
!       IEMOM: Electrostatic moments calculation
!        EVEC: An array of the three x,y,z components of
!        IECP: Effective Core Potentials
!     IHSSCAL: Compute Hessian and vibrational analysis if IRUNTYP=3
!    IPROJECT: Project Hessian to eliminate rot/vib contaminants
!      ISIGMA: Rotational symmetric number for thermochemistry
!
!      NATmax: Maximum Number of Atoms
!   NSHELLmax: Maximum Number of Shells
!   NPRIMImax: Maximum Number of Gaussian Functions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL NAMELIST_INPRUN(IRUNTYP,ICH,MUL,NINTEG,IDONTW,IEMOM,EVEC,    &
                           ILIBCINT,IECP,IHSSCAL,IPROJECT,ISIGMA,IGTYP, &
                           NATmax,NSHELLmax,NPRIMImax,IHUBBARD)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      if(IHUBBARD==1)then         !     Hubbard Calculation
!-----------------------------------------------------------------------
       CALL HUBBARD(IRUNTYP,ICH,MUL,IDONTW)      
!-----------------------------------------------------------------------      
      else                        !     Molecular Calculation
!-----------------------------------------------------------------------
!      Initialize for the integral quadratures if ERI HONDO Calculator
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(ILIBCINT==0)CALL INIINTQUAD
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ALLOCATE(ZAN0(NATmax),Cxyz0(3,NATmax),IAN0(NATmax),IMIN0(NATmax),&
                IMAX0(NATmax),ZMASS0(NATmax),KSTART0(NSHELLmax),        &
                KATOM0(NSHELLmax),KTYPE0(NSHELLmax),KLOC0(NSHELLmax),   &
                INTYP0(NSHELLmax),KNG0(NSHELLmax),KMIN0(NSHELLmax),     &
                KMAX0(NSHELLmax),ISH0(NPRIMImax),ITYP0(NPRIMImax),      &
                C10(NPRIMImax),C20(NPRIMImax),EX0(NPRIMImax),           &
                CS0(NPRIMImax),CP0(NPRIMImax),CD0(NPRIMImax),           &
                CF0(NPRIMImax),CG0(NPRIMImax),CH0(NPRIMImax),           &
                CI0(NPRIMImax))

        SIZE_ENV = 20 + 3*NATmax + 2*NPRIMImax !LIBCINT
        ALLOCATE(ENV0(SIZE_ENV))               !LIBCINT
        ALLOCATE(ATM0(6,NATmax))               !LIBCINT
        ALLOCATE(BAS0(8,NSHELLmax))            !LIBCINT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Read in Basis Set and get initial Molecular Orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL START(IGTYP,NAT,NATmax,NBF,NBFaux,NQMT,NE,NA,NB,NSHELL,     &
                  NSHELLmax,NSHELLaux,NPRIMI,NPRIMImax,NPRIMIaux,       &
                  ZAN0,Cxyz0,IAN0,IMIN0,IMAX0,ZMASS0,KSTART0,KATOM0,    &
                  KTYPE0,KNG0,KLOC0,KMIN0,KMAX0,INTYP0,ISH0,ITYP0,      &
                  C10,C20,EX0,CS0,CP0,CD0,CF0,CG0,CH0,CI0,SIZE_ENV,ENV0,&
                  ATM0,BAS0) !LIBCINT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      NAT: Number of Atoms             
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       allocate(Cxyz(3,NAT))
       do i=1,3
        Cxyz(i,1:nat) = Cxyz0(i,1:nat)
       end do
       deallocate(Cxyz0)
       IF(NAT==1.and.IRUNTYP==3)THEN
        WRITE(6,5)
        IRUNTYP = 2
       END IF
!      
       allocate(IAN(NAT),IMIN(NAT),IMAX(NAT),ZAN(NAT),ZMASS(NAT))
       IAN(1:nat)   = IAN0(1:nat)
       IMIN(1:nat)  = IMIN0(1:nat)
       IMAX(1:nat)  = IMAX0(1:nat)
       ZAN(1:nat)   = ZAN0(1:nat)
       ZMASS(1:nat) = ZMASS0(1:nat)
       deallocate(IAN0,IMIN0,IMAX0,ZAN0,ZMASS0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      NSHELL: Total number of shells
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       allocate(KSTART(NSHELL),KATOM(NSHELL),KTYPE(NSHELL),KLOC(NSHELL))
       allocate(INTYP(NSHELL),KNG(NSHELL),KMIN(NSHELL),KMAX(NSHELL))
       KSTART(1:nshell)= KSTART0(1:nshell)
       KATOM(1:nshell) = KATOM0(1:nshell)
       KTYPE(1:nshell) = KTYPE0(1:nshell)
       KLOC(1:nshell)  = KLOC0(1:nshell) 
       INTYP(1:nshell) = INTYP0(1:nshell)
       KNG(1:nshell)   = KNG0(1:nshell)  
       KMIN(1:nshell)  = KMIN0(1:nshell) 
       KMAX(1:nshell)  = KMAX0(1:nshell) 
       deallocate(KSTART0,KATOM0,KTYPE0,KLOC0,INTYP0,KNG0,KMIN0,KMAX0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      NPRIMI: Total number of primitive exponents
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       allocate(ISH(nprimi),ITYP(nprimi),C1(nprimi),C2(nprimi))
       allocate(EX(nprimi),CS(nprimi),CP(nprimi),CD(nprimi))
       allocate(CF(nprimi),CG(nprimi),CH(nprimi),CI(nprimi))
       ISH(1:nprimi)  = ISH0(1:nprimi) 
       ITYP(1:nprimi) = ITYP0(1:nprimi) 
       C1(1:nprimi) = C10(1:nprimi)
       C2(1:nprimi) = C20(1:nprimi)
       EX(1:nprimi) = EX0(1:nprimi)
       CS(1:nprimi) = CS0(1:nprimi)
       CP(1:nprimi) = CP0(1:nprimi)
       CD(1:nprimi) = CD0(1:nprimi)
       CF(1:nprimi) = CF0(1:nprimi)
       CG(1:nprimi) = CG0(1:nprimi)
       CH(1:nprimi) = CH0(1:nprimi)
       CI(1:nprimi) = CI0(1:nprimi)      
       deallocate(ISH0,ITYP0,C10,C20,EX0,CS0,CP0,CD0,CF0,CG0,CH0,CI0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Copy LIBCINT variables
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       SIZE_ENV = 20 + 3*NAT + 2*NPRIMI + 2*NPRIMIaux !LIBCINT
       NBAS = NSHELL + NSHELLaux                      !LIBCINT
       ALLOCATE(ENV(SIZE_ENV))                        !LIBCINT
       ALLOCATE(ATM(6,NAT))                           !LIBCINT
       ALLOCATE(BAS(8,NBAS))                          !LIBCINT
       ENV(1:SIZE_ENV) = ENV0(1:SIZE_ENV)             !LIBCINT
       ATM(1:6, 1:NAT) = ATM0(1:6, 1:NAT)             !LIBCINT
       BAS(1:8, 1:NBAS) = BAS0(1:8, 1:NBAS)           !LIBCINT
       DEALLOCATE(ENV0,ATM0,BAS0)                     !LIBCINT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Header on the output file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!         NAT: Number of Atoms             
!         ICH: Charge of Molecule
!         MUL: State Multiplicity
!         NBF: Number of Basis Functions (NSQ=NBF*NBF,NBFT=NBF(NBF+1)/2)
!        NQMT: Number of linearly independent orbitals  
!          NE: Number of Electrons
!          NA: Number of Alpha electrons
!          NB: Number of Beta electrons
!        EVEC: Electric Field components
!      NSHELL: Total number of shells
!      NPRIMI: Total number of primitive exponents
!         IAN: True nuclear charge
!       IEMOM: Electrostatic moments calculation
!        IECP: Effective Core Potentials
!         ZAN: Effective Nuclear Charge
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL RUNNOFHEADER(NAT,ICH,MUL,NBF,NBFaux,NQMT,NE,NA,NB,          &
                         EVEC(1),EVEC(2),EVEC(3),NSHELL,NSHELLaux,      &
                         NPRIMI,IAN,IEMOM,IECP,IRUNTYP,Cxyz,ZAN)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      If RUNTYP = 1) ENERGY  : single-point energy
!                  2) GRADIENT: single-point energy + gradients
!                  3) OPTIMIZE: optimize the molecular geometry
!                  4) HESS    : compute numerical hessian 
!                  5) DYN     : run Born-Oppenheimer molecular dynamics
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ALLOCATE(GRADS(3*NAT))
       IF(IRUNTYP==1.or.IRUNTYP==2)THEN
        CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI, &
                      ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,   &
                      INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,   &
                      CF,CG,CH,CI,GRADS,IRUNTYP,DIPS,SIZE_ENV,ENV,ATM,  &
                      NBAS,BAS,IGTYP,1,1)
       ELSE IF(IRUNTYP==3)THEN
        CALL OPTIMIZE(NINTEG,IDONTW,NAT,ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,   &
                      KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,  &
                      ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,DIPS,SIZE_ENV, &
                      ENV,ATM,NBAS,BAS,IGTYP,GRADS,IRUNTYP,IHSSCAL,     &
                      IPROJECT,ISIGMA)
       ELSE IF(IRUNTYP==4)THEN
        CALL HESSCAL(NINTEG,IDONTW,NAT,ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,    &
                     KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,   &
                     ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,DIPS,GRADS,     &
                     IRUNTYP,IPROJECT,ISIGMA,SIZE_ENV,ENV,ATM,          &
                     NBAS,BAS,IGTYP)
       ELSE IF(IRUNTYP==5)THEN
        CALL MOLDYN (NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
                     ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,KATOM,KTYPE,   &
                     KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,  &
                     CD,CF,CG,CH,CI,GRADS,IRUNTYP,DIPS,SIZE_ENV,ENV,ATM,&
                     NBAS,BAS,IGTYP)                     
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DEALLOCATE(ZAN,ZMASS,Cxyz,GRADS,IAN,IMIN,IMAX)
       DEALLOCATE(KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH)
       DEALLOCATE(ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI)
!-----------------------------------------------------------------------
      end if
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(TIMEFINISH)
      DELTATIME = (TIMEFINISH - TIMESTART)/RATE
      WRITE(6,3)DELTATIME
!-----------------------------------------------------------------------
      call date_and_time(date,time,zone,val)
      write(6,4)val(5),val(6),val(2),val(3),val(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     MPI finalization
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,MY_ID,IERR)
      IF(MY_ID == 0) THEN
        INTTYPE = 3
        CALL MPI_BCAST(INTTYPE,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      END IF
  10  CALL MPI_FINALIZE(ierror)
#endif
!-----------------------------------------------------------------------
!     Format definitions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    1 FORMAT(/' The execution started at ',I2,'h ',I2,'min on the ',    &
                I2,'/',I2,'/',I4)
    2 FORMAT(                                                           &
      /4X,'===========================================================' &
      /4X,'!                                                         !',&                                                                     
      /4X,'!                        Do N O F                         !',&   
      /4X,'!                                                         !',&                                                                     
      /4X,'!  (Donostia Natural Orbital Functional Software Program) !',&   
      /4X,'!                                                         !',&                                                                     
      /4X,'!               COPYRIGHT by Mario Piris                  !',& 
      /4X,'!                                                         !',& 
      /4X,'!      Donostia International Physics Center (DIPC)       !',&
      /4X,'!       University of the Basque Country (UPV/EHU)        !',&
      /4X,'!       Basque Foundation for Science (IKERBASQUE)        !',&
      /4X,'!                                                         !',& 
      /4X,'!          GNU General Public License version 3           !',&                                                                     
      /4X,'!                                                         !',&                                                                     
      /4X,'!                  VERSION: November 2024                 !',&
      /4X,'!                                                         !',&
      /4X,'===========================================================')
    3 FORMAT(/,'  Elapsed real time :',F10.2,'  (Seconds)')
    4 FORMAT(/' The execution finished at ',I2,'h ',I2,'min on the ',   &
                I2,'/',I2,'/',I4)                                                                                                                       
    5 FORMAT(/1X,'Warning: For geometry optimization the number of'     &
                 'atoms must be greater,',/,10X,'than 1 '               &
                 'so RUNTYP has been set equal to GRAD, not OPTGEO')
!-----------------------------------------------------------------------
      STOP
      END
