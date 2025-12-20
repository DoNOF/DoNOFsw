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
!                           Date: October 2025                         !
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
      DOUBLE PRECISION,DIMENSION(3) :: DIPS
      INTEGER :: INTTYPE
      INTEGER :: CR, CM, TIMESTART, TIMEFINISH
      DOUBLE PRECISION :: RATE
!     LIBCINT
      INTEGER :: SIZE_ENV, NBAS
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: ENV0
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ATM0, BAS0
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: ENV
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ATM, BAS
      INTEGER :: NPRIMIecp,NSHELLecp      
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
!        NLOP: Non-linear optical property calculation (-1,0,1,2,3)
!          NP: Number of steps used in the dyadic scaling of the field
!        STEP: Initial step size for the electric field in NLOP
!    ISOALPHA: Isotropic average and the anisotropy (Raman convention)
! EFX,EFY,EFZ: The x,y,z components of the electric eield
!        IECP: Effective Core Potentials
!     IHSSCAL: Compute Hessian and vibrational analysis if IRUNTYP=3
!    IPROJECT: Project Hessian to eliminate rot/vib contaminants
!      ISIGMA: Rotational symmetric number for thermochemistry
!
!      NATmax: Maximum Number of Atoms
!   NSHELLmax: Maximum Number of Shells
!   NPRIMImax: Maximum Number of Gaussian Functions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL NAMELIST_INPRUN(IRUNTYP,ICH,MUL,NINTEG,IDONTW,IEMOM,NLOP,NP, &
                           STEP,ISOALPHA,EFX,EFY,EFZ,ILIBCINT,IECP,     &
                           IHSSCAL,IPROJECT,ISIGMA,IGTYP,NATmax,        &
                           NSHELLmax,NPRIMImax,IHUBBARD)
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
!HONDO IF(ILIBCINT==0)CALL INIINTQUAD
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
!       LIBCINT
        SIZE_ENV = 20 + 3*NATmax + 2*NPRIMImax + 1000
        ALLOCATE(ENV0(SIZE_ENV))
        ALLOCATE(ATM0(6,NATmax))
        ALLOCATE(BAS0(8,NSHELLmax + 1000))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Read in Basis Set and get initial Molecular Orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL START(IGTYP,NAT,NATmax,NBF,NBFaux,NQMT,NE,NA,NB,NSHELL,     &
                  NSHELLmax,NSHELLaux,NPRIMI,NPRIMImax,NPRIMIaux,       &
                  ZAN0,Cxyz0,IAN0,IMIN0,IMAX0,ZMASS0,KSTART0,KATOM0,    &
                  KTYPE0,KNG0,KLOC0,KMIN0,KMAX0,INTYP0,ISH0,ITYP0,      &
                  C10,C20,EX0,CS0,CP0,CD0,CF0,CG0,CH0,CI0,NPRIMIecp,    &
                  NSHELLecp,SIZE_ENV,ENV0,ATM0,BAS0)
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
       SIZE_ENV = 20 + 3*NAT + 2*NPRIMI + 2*NPRIMIaux + 2*NPRIMIecp
       NBAS = NSHELL + NSHELLaux + NSHELLecp
       ALLOCATE(ENV(SIZE_ENV))
       ALLOCATE(ATM(6,NAT))
       ALLOCATE(BAS(8,NBAS))
       ENV(1:SIZE_ENV) = ENV0(1:SIZE_ENV)
       ATM(1:6, 1:NAT) = ATM0(1:6, 1:NAT)
       BAS(1:8, 1:NBAS) = BAS0(1:8, 1:NBAS)
       DEALLOCATE(ENV0,ATM0,BAS0)
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
! EFX,EFY,EFZ: Electric Field components
!      NSHELL: Total number of shells
!      NPRIMI: Total number of primitive exponents
!         IAN: True nuclear charge
!       IEMOM: Electrostatic moments calculation
!        IECP: Effective Core Potentials
!         ZAN: Effective Nuclear Charge
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL RUNNOFHEADER(NAT,ICH,MUL,NBF,NBFaux,NQMT,NE,NA,NB,EFX,EFY,  &
                         EFZ,NSHELL,NSHELLaux,NPRIMI,IAN,IEMOM,IECP,    &
                         IRUNTYP,Cxyz,ZAN)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      If RUNTYP = 1) ENERGY  : single-point energy
!                  2) GRADIENT: single-point energy + gradients
!                  3) OPTIMIZE: optimize the molecular geometry
!                  4) HESS    : compute numerical hessian 
!                  5) DYN     : run Born-Oppenheimer molecular dynamics
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ALLOCATE(GRADS(3*NAT))
       IF(IRUNTYP==1)THEN
        if(NLOP/=0)then  ! compute non-linear optical properties
         CALL DIPNLOP(NP,STEP,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,       &
                      NSHELL,NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,      &
                      KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,    &
                      C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,      &
                      DIPS,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,NLOP,ISOALPHA)
        else
         CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,&
                       ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,  &
                       INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,  &
                       CF,CG,CH,CI,GRADS,IRUNTYP,DIPS,SIZE_ENV,ENV,ATM, &
                       NBAS,BAS,IGTYP,1,1)
        end if
       ELSE IF(IRUNTYP==2)THEN
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
       DEALLOCATE(ENV,ATM,BAS)
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
        CALL MPI_BCAST(INTTYPE,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
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
      /4X,'!                  VERSION: October 2025                  !',&
      /4X,'!                                                         !',&
      /4X,'===========================================================')
    3 FORMAT(/,'  Elapsed real time :',F15.2,'  (Seconds)')
    4 FORMAT(/' The execution finished at ',I2,'h ',I2,'min on the ',   &
                I2,'/',I2,'/',I4)                                                                                                                       
    5 FORMAT(/1X,'Warning: For geometry optimization the number of',    &
                 'atoms must be greater,',/,10X,'than 1 ',              &
                 'so RUNTYP has been set equal to GRAD, not OPTGEO')
!-----------------------------------------------------------------------
      STOP
      END

!======================================================================!
!                                                                      !
!   RUNNOFHEADER: Write NOF header on the output file                  !
!   OPENFILES: Open all general working files.                         !
!   OUTPUTBASIC: Write the basic info on the output file.              !
!   SETNO1: Determine NO1 according to true nuclear charges if NO1=-1  !
!   SETORBSPACE: Define the orbital space (NDOC,NSOC,NCWO,NVIR,...)    !
!   POINTERS: Define Pointers of the USER array for the CG subroutine  !
!                                                                      !
!======================================================================!

! RUNNOFHEADER
      SUBROUTINE RUNNOFHEADER(NATOMSn,ICHn,MULn,NBFn,NBFauxn,NQMTn,NEn, &
                              NAn,NBn,EX,EY,EZ,NSHELLn,NSHELLauxn,      &
                              NPRIMIn,IAN,IEMOMn,IECPn,IRUNTYP,Cxyz,    &
                              ZNUC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EFIELDL,RESTART,ERIACTIVATED
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPFILE_NO1PT2/NO1PT2,NEX
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      COMMON/ELPROP/IEMOM
      COMMON/ECP2/CLP(4004),ZLP(4004),NLP(4004),KFRST(1001,6),          &
      KLAST(1001,6),LMAX(1001),LPSKIP(1001),IZCORE(1001)
      COMMON/NumLinIndOrb/NQMT
      !
      INTEGER,DIMENSION(NATOMSn) :: IAN
      DOUBLE PRECISION,DIMENSION(3,NATOMSn) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NATOMSn) :: ZNUC
      CHARACTER*4,ALLOCATABLE,DIMENSION(:) :: ATMNAME
!-----------------------------------------------------------------------
!     Basic information
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NATOMS: Number of Atoms
!        ICH: Charge of Molecule
!        MUL: State Multiplicity
!        NBF: Number of Basis Functions (NSQ=NBF*NBF,NBFT=NBF(NBF+1)/2)
!       NQMT: Number of linearly independent orbitals
!         NE: Number of Electrons
!         NA >= NB
!         NA: Number of Alpha electrons
!         NB: Number of Beta electrons
!EFX,EFY,EFZ: Electric Field components
!     NSHELL: Total number of shells
!     NPRIMI: Total number of primitive exponents
!- - - - - - - - - - - - - - - - - - - - - - - - -
!     NBFaux: Number of Auxiliary Basis Functions
!  NSHELLaux: Total number of auxiliary shells
!-----------------------------------------------------------------------
      NATOMS = NATOMSn
      ICH    = ICHn
      MUL    = MULn
      NBF    = NBFn
      NQMT   = NQMTn
      NE     = NEn
      NA     = NAn
      NB     = NBn
      EFX    = EX
      EFY    = EY
      EFZ    = EZ
      NSHELL = NSHELLn
      NPRIMI = NPRIMIn
      IEMOM  = IEMOMn
      IECP   = IECPn
!- - - - - - - - - - - - - -
      NBFaux = NBFauxn
      NSHELLaux = NSHELLauxn
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NBFT: Dimension for symmetric matices
!      NSQ: Dimension for square matices
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NSQ = NBF*NBF
      NBFT = NBF*(NBF+1)/2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     State of ERIs in Nodes
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ERIACTIVATED = .FALSE.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Input namelist variables: ICOEF, MAXIT, ...
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL NAMELIST_NOFINP(IRUNTYP)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Header of NOF Calculation on the output file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL OUTPUTHEADER
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Open files (GCF,BFST,GCFe,WFN,FCHK,APSG,FRAG,CGGRAD)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL OPENFILES(IRUNTYP)
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
!     INPUTCXYZ=0: Read geometry from input file
!     INPUTCXYZ=1: Read geometry from GCF file
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
      IF(INPUTCXYZ==1) THEN
       CALL READCXYZ(ZNUC,Cxyz,NATOMS,NBF,NSQ)
       WRITE(6,2)
       ALLOCATE(ATMNAME(NATOMS))
       CALL ATOMNAMES(NATOMS,ZNUC,IZCORE,ATMNAME,Cxyz,1,1)
       DEALLOCATE(ATMNAME)
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Output Basic
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL OUTPUTBASIC
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     If NO1 = -1 calculate NO1 according to true nuclear charges (IAN)
!     NO1: Natural Orbitals with Occupation Numbers equal to one.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NO1==-1)THEN
!      If using ECPotentials (IECP/=0): NO1=0, NO1PT2=0
       IF(IECP/=0)THEN
        NO1 = 0
        NO1PT2 = 0
        WRITE(6,1)
       ELSE
        CALL SETNO1(IAN)
       ENDIF
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NO1:  Number of inactive doubly occupied orbitals (OCC=1)
!     NDOC: Number of strongly doubly occupied MOs
!     NSOC: Number of strongly singly occupied MOs
!     NDNS: Number of strongly occupied MOs (NDNS=NDOC+NSOC)
!     NCWO: Number of coupled weakly occ. MOs per strongly doubly occ.
!     NCWO*NDOC: Active orbitals in the virtual subspace
!     NO0:  Empty orbitals  (OCC=0)
!     NVIR: Number of weakly occupied MOs + empty MOs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!           NO1 | NDOC  + NSOC  |   NCWO*NDOC + NO0  = NBF
!           NO1 |      NDNS     |          NVIR      = NBF
!               | -NAC- |       |  -   NAC  - |
!                      NB      NA            NBF5
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL SETORBSPACE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Define Pointers of the USER array for the external CG subroutine
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL POINTERS
!-----------------------------------------------------------------------
    1 FORMAT(/1X,'You are using an ECP: Core Orbitals have been already'&
      ,1X,'excluded.',/1X,'NO1 and NO1PT2 = 0.')
    2 FORMAT(/1X,'Nuclear Coordinates from GCF file:'                   &
      /1X,'----------------------------------')
!-----------------------------------------------------------------------
      RETURN
      END

! OUTPUTHEADER
      SUBROUTINE OUTPUTHEADER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL DIIS,PERDIIS,DAMPING,EXTRAP,RESTART,PRINTLAG,DIAGLAG
      LOGICAL APSG,CHKORTHO,ORTHO,HighSpin,SCALING
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPSCALING/SCALING,NZEROS,NZEROSm,NZEROSr,ITZITER
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_DIIS/DIIS,PERDIIS,NDIIS,NTHDIIS,THDIIS
      COMMON/INPNOF_DAMPEXTRAP/DAMPING,EXTRAP
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_LAGRANGE/PRINTLAG,DIAGLAG
      COMMON/INPNOF_APSG/APSG,NTHAPSG,THAPSG
      COMMON/INPNOF_ORTHOGONALITY/CHKORTHO,ORTHO
      COMMON/INPNOF_RHF/CONVRHFDM,IRHF,IRHFTYP,NCONVRHF,MAXITRHF
      COMMON/INPNOF_HFID/KOOPMANS
      COMMON/INPNOF_EKT/IEKT
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/INPNOF_INICON/INICOND
      COMMON/INPNOF_ARDM/THRESHDM,NOUTRDM,NSQT,NTHRESHDM
      COMMON/INPNOF_CJK/NOUTCJK,NTHRESHCJK,THRESHCJK
      COMMON/INPNOF_Tijab/NOUTTijab,NTHRESHTijab,THRESHTijab
      COMMON/INPNOF_NTHRESH/NTHRESHL,NTHRESHE,NTHRESHEC,NTHRESHEN
      COMMON/INPNOF_COEFOPT/MAXLOOP
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPNOF_STATIC/Ista
      COMMON/INPNOF_MOD/Imod
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/INPNOF_EXSTA/OMEGA1,NESt
      INTEGER :: NESt
      DOUBLE PRECISION :: OMEGA1
!-----------------------------------------------------------------------
!     Set NZEROSm equal to NTHRESHL if NZEROSm < NTHRESHL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IORBOPT==1.and.SCALING.and.NZEROSm<NTHRESHL)THEN
       NZEROSm = NTHRESHL
       WRITE(6,'(/,6X,38A,/)')'NZEROSm has been set equal to NTHRESHL'
      ENDIF
!-----------------------------------------------------------------------
!     Write NAMELIST parameters on the Output file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,3)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write the Functional used
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(3<=IPNOF.and.IPNOF<=8)THEN
        WRITE(6,7)IPNOF
        if(IPNOF==7.and.Ista==1)WRITE(6,*)                              &
          'Static version of PNOF7 (PNOF7s):       (Ista)          1'
        if(IPNOF==8.and.Imod==1)WRITE(6,*)                              &
          'Modified version of GNOF (GNOFm):       (Imod)          1'
        if(IPNOF==8.and.Imod==2)WRITE(6,*)                              &
          'Scaled version of GNOF (GNOFs):         (Imod)          2'
      ELSE
        WRITE(6,*)'Stop Program: Select IPNOF between 3 and 8'
        CALL ABRT
      ENDIF
!     Multiplicity
      IF(MUL>1)THEN
        IF(.NOT.(IPNOF==5.or.IPNOF==7.or.IPNOF==8))THEN
        WRITE(6,'(/A44)')' Stop: IPNOF must be equal 5,7,8 for MULT>1'
        CALL ABRT
      ENDIF
      IF(HighSpin)THEN
        WRITE(6,*)                                                      &
        'High-Spin State calculation Option:     (HighSpin)      T'
      ELSE
        WRITE(6,*)                                                      &
        'Spin-Multiplet, High-Spin State Option: (HighSpin)      F'
        ENDIF
      ENDIF
!     Orbital Optimization
      IF(ICOEF==0)THEN
        WRITE(6,40)ICOEF
      ELSEIF(ICOEF==1)THEN
        WRITE(6,41)ICOEF
      ELSEIF(ICOEF==2)THEN
        WRITE(6,42)ICOEF
      ELSEIF(ICOEF==21)THEN
        WRITE(6,43)ICOEF
      ELSEIF(ICOEF==3)THEN
        WRITE(6,44)ICOEF
      ELSE
        WRITE(6,*)'Wrong ICOEF'
        STOP
      ENDIF
!
      IF(ISOFTMAX==0)THEN
        WRITE(6,50)ISOFTMAX
      ELSE IF(ISOFTMAX==1)THEN
        WRITE(6,51)ISOFTMAX
      ENDIF
!
      IF(ICOEF>0)THEN
       IF(IORBOPT==1)THEN
         WRITE(6,55)IORBOPT
       ELSE IF(2<=IORBOPT.and.IORBOPT<=5)THEN
         WRITE(6,56)IORBOPT
       ELSE IF(IORBOPT==6)THEN
         WRITE(6,57)IORBOPT
       ELSE
        WRITE(6,*)'Wrong Orb. Opt. option'
        STOP
       END IF
       IF(ICOEF==21)WRITE(6,60)MAXIT21
       WRITE(6,61)MAXIT
      END IF
!     Restart Option
      WRITE(6,*)
      IF(RESTART)THEN
        WRITE(6,*)                                                      &
        'Restart calculation Option:             (RESTART)       T'
      ELSE
        WRITE(6,*)                                                      &
        'Restart calculation Option:             (RESTART)       F'
      ENDIF
      IF(.NOT.RESTART)WRITE(6,8)INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
!
      IF(IRHF==0)THEN
       WRITE(6,81)
      ELSE IF(IRHF==1)THEN
       WRITE(6,82)
      ELSE IF(IRHF==2)THEN
       WRITE(6,83)
      ELSE IF(IRHF==3)THEN
       WRITE(6,84)KOOPMANS
      ENDIF
!
      IF(ICOEF==1.or.ICOEF==2.or.ICOEF==21)THEN
       IF(IORBOPT==1)THEN
        WRITE(6,9)NTHRESHL,NTHRESHE,NTHRESHEC,NTHRESHEN,MAXLOOP
       ELSE IF(IORBOPT==2)THEN
        WRITE(6,91)NTHRESHL,NTHRESHE,MAXLOOP
       ENDIF
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ORTHO)WRITE(6,*)                                               &
        'Orthogonalize the initial Orbitals:     (ORTHO)         T'
      IF(CHKORTHO)WRITE(6,*)                                            &
        'Check the orthonormality of MOs:        (CHKORTHO)      T'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF( IORBOPT==1 .and. (ICOEF==1.or.ICOEF==21) .and. SCALING )THEN
        WRITE(6,10)NZEROS,NZEROSm,NZEROSr,ITZITER
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF( IORBOPT==1 .and. (.not.SCALING) )THEN
        WRITE(6,*)
        WRITE(6,*)'Warning: Scaling Technique is not used !'
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF( IORBOPT==1 .and. (ICOEF==1.or.ICOEF==21) )THEN
        IF(DIIS)THEN
          WRITE(6,*)
          WRITE(6,*)                                                    &
          'DIIS technique is used in Orb. Opt.:    (DIIS)          T'
          WRITE(6,13)NTHDIIS,NDIIS
         IF(PERDIIS)THEN
          WRITE(6,*)                                                    &
          'Periodic DIIS every NDIIS:              (PERDIIS)       T'
         ELSE
          WRITE(6,*)                                                    &
          'DIIS is always applied after NDIIS      (PERDIIS)       F'
         ENDIF
        ELSE
         WRITE(6,*)                                                     &
         'DIIS Technique is not used:             (DIIS)          F'
        ENDIF
        IF(DAMPING)THEN
          WRITE(6,*)                                                    &
          'Damping of the Gen. Fock matrix:        (DAMPING)       T'
        ENDIF
        IF(EXTRAP)THEN
          WRITE(6,*)                                                    &
          'Extrapolation of the Gen. Fock matrix:  (EXTRAP)        T'
        ENDIF
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,11)NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK,MOLDEN,INICOND,IEKT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(APSG)THEN
        WRITE(6,*)                                                      &
        'Write APSG expansion coefficient file:  (APSG)          T'
        WRITE(6,14)NTHAPSG
      ELSE
        WRITE(6,*)                                                      &
        'Write APSG expansion coefficient file:  (APSG)          F'
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF( ICOEF==1 .or. ICOEF==21 )THEN
        IF(PRINTLAG)THEN
          WRITE(6,*)                                                    &
          'Output option for Lagrange Multipliers: (PRINTLAG)      T'
        ELSE
          WRITE(6,*)                                                    &
          'Output option for Lagrange Multipliers: (PRINTLAG)      F'
        ENDIF
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(DIAGLAG)THEN
        WRITE(6,*)                                                      &
        'Diagonalize matrix of Lag. Multipliers: (DIAGLAG)       T'
      ELSE
        WRITE(6,*)                                                      &
        'Diagonalize matrix of Lag. Multipliers: (DIAGLAG)       F'
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF( NOUTRDM==1 .or. NOUTRDM==2 .or. NOUTRDM==3 )THEN
        WRITE(6,12)NOUTRDM,NTHRESHDM
      ENDIF
      IF(NOUTCJK==1)WRITE(6,15)NOUTCJK,NTHRESHCJK
      IF(NOUTTijab==1)WRITE(6,16)NOUTTijab,NTHRESHTijab
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NESt>0)WRITE(6,17)NESt,OMEGA1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
!-----------------------------------------------------------------------
!     Format definitions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    3 FORMAT(/,' Input NOF Options',/,                                  &
      ' -----------------')
   40 FORMAT(                                                           &
      /1X,'Type of calculation = No Coeff. Opt.:   (ICOEF)     ',I5)
   41 FORMAT(                                                           &
      /1X,'Type of calculation = Coeff. Opt.:      (ICOEF)     ',I5)
   42 FORMAT(                                                           &
      /1X,'Type of calculation = No Occup. Opt.:   (ICOEF)     ',I5)
   43 FORMAT(                                                           &
      /1X,'Type of calculation = ICOEF2 + ICOEF1:  (ICOEF)     ',I5)
   44 FORMAT(                                                           &
      /1X,'Type of calculation = Fragment Calc.:   (ICOEF)     ',I5)
   50 FORMAT(                                                           &
       1X,'Trigonometric functions for Occ. Opt.:  (ISOFTMAX)  ',I5)
   51 FORMAT(                                                           &
       1X,'Use Softmax function for Occ. Opt.:     (ISOFTMAX)  ',I5)
   55 FORMAT(                                                           &
       1X,'Type of Orb. Opt. = Iterative Diag.:    (IORBOPT)   ',I5)
   56 FORMAT(                                                           &
       1X,'Type of Orb. Opt. = Orbital Rotations:  (IORBOPT)   ',I5)
   57 FORMAT(                                                           &
       1X,'Type of Orb. Opt. = Seq. Quad. Prog.:   (IORBOPT)   ',I5)
   60 FORMAT(                                                           &
       1X,'Maximum Number of Coef Outer Iter.:     (MAXIT21)   ',I5)
   61 FORMAT(                                                           &
       1X,'Maximum Number of Occ-Coef Outer Iter.: (MAXIT)     ',I5)
    7 FORMAT(                                                           &
      /1X,'Natural Orbital Functional Selected:    (IPNOF)     ',I5)
    8 FORMAT(                                                           &
      1X,'Restart Gamma Matrix from GCF file:     (INPUTGAMMA)',I5,    &
      /1X,'Coefficient Matrix from GCF file:       (INPUTC)    ',I5,    &
      /1X,'Diagonal Elements FMIUG from GCF file:  (INPUTFMIUG)',I5,    &
      /1X,'Cartesian Coordinates from GCF file:    (INPUTCXYZ) ',I5)
   81 FORMAT(                                                           &
      /1X,'No Hartree-Fock Calculation:            (IRHF)', 10X,'0')
   82 FORMAT(                                                           &
      /1X,'Hartree-Fock Calc. using SCF Method:    (IRHF)', 10X,'1')
   83 FORMAT(                                                           &
      /1X,'Hartree-Fock Calc. using ADAM Method:   (IRHF)', 10X,'2')
   84 FORMAT(                                                           &
      /1X,'Hartree-Fock Calc. using ID Method:     (IRHF)', 10X,'3',    &
      /1X,'Ion. Potentials by Koopmans Theorem:    (KOOPMANS)  ',I5)
    9 FORMAT(                                                           &
      /1X,'Threshold Lambda Convergence=10**(-N):  (NTHRESHL)  ',I5,    &
      /1X,'Threshold Energy Convergence=10**(-N):  (NTHRESHE)  ',I5,    &
      /1X,'Threshold Energy Convergence=10**(-N):  (NTHRESHEC) ',I5,    &
      /1X,'Threshold Energy Convergence=10**(-N):  (NTHRESHEN) ',I5,    &
      //1X,'Max. Number of Inner Coef. Iterations: (MAXLOOP)   ',I5)
   91 FORMAT(                                                           &
      /1X,'Threshold Lambda Convergence=10**(-N):  (NTHRESHL)  ',I5,    &
      /1X,'Threshold Energy Convergence=10**(-N):  (NTHRESHE)  ',I5,    &
      /1X,'Max. Number of Inner Coef. Iterations:  (MAXLOOP)   ',I5)
   10 FORMAT(                                                           &
      /1X,'Scaling Parameters:',                                        &
      /1X,'Initial Number of Zeros IN Fij:         (NZEROS)    ',I5,    &
      /1X,'Maximum Number of Zeros IN Fij:         (NZEROSm)   ',I5,    &
      /1X,'Restart Number of Zeros IN Fij:         (NZEROSr)   ',I5,    &
      /1X,'Number of Iter with constant Scaling:   (ITZITER)   ',I5)
   11 FORMAT(                                                           &
      /1X,'Output Option:                          (NPRINT)    ',I5,    &
      /1X,'Output the Coefficient Matrix:          (IWRITEC)   ',I5,    &
      /1X,'Do a Mulliken Population Analysis:      (IMULPOP)   ',I5,    &
      /1X,'Write Information into a WFN file:      (IAIMPAC)   ',I5,    &
      /1X,'Write Information into a FCHK file:     (IFCHK)     ',I5,    &
      /1X,'Write Information into a MLP file:      (MOLDEN)    ',I5,    &
      /1X,'Create Ini. Cond. (ini.xyz) file:       (INICOND)   ',I5,    &
      /1X,'Calculate IPs using Ext. Koopmans Theo: (IEKT)      ',I5)
   12 FORMAT(                                                           &
      /1X,'Print atomic RDMs to files 1DM and 2DM: (NOUTRDM)   ',I5,    &
      /1X,'Threshold DMs = 10.0**(-NTHRESHDM):     (NTHRESHDM) ',I5)
   13 FORMAT(                                                           &
      1X,'Threshold to begin DIIS = 10**(-N):      (NTHDIIS)   ',I5,    &
      /1X,'Number of considered Loops in DIIS:     (NDIIS)     ',I5)
   14 FORMAT(                                                           &
      1X,'Threshold APSG Exp. Coef. = 10**(-N):    (NTHAPSG)   ',I5)
   15 FORMAT(                                                           &
      /1X,'Print CJ12 and CK12 to file CJK:        (NOUTCJK)   ',I5,    &
      /1X,'Threshold CJKs = 10.0**(-NTHRESHCJK):   (NTHRESHCJK)',I5)
   16 FORMAT(                                                           &
      /1X,'Print OIMP2 Ampl. Tijab to file Tijab:  (NOUTTijab) ',I5,    &
      /1X,'Threshold Tijab=10.0**(-NTHRESHTijab):  (NTHRESHTijab)',I3)
   17 FORMAT(                                                           &
      /1X,'No. Excited States in w-ensemble NOFT:  (NESt)      ',I5,    &
      /1X,'Value for w1 in W = (w1,1-w1,0,...,0):  (OMEGA1)    ',F5.1)
!-----------------------------------------------------------------------
      END

! OPENFILES
      SUBROUTINE OPENFILES(IRUNTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL APSG
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_APSG/APSG,NTHAPSG,THAPSG
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/INPNOF_ARDM/THRESHDM,NOUTRDM,NSQT,NTHRESHDM
      COMMON/INPNOF_CJK/NOUTCJK,NTHRESHCJK,THRESHCJK
      COMMON/INPNOF_Tijab/NOUTTijab,NTHRESHTijab,THRESHTijab
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
!-----------------------------------------------------------------------
!     NOF Files    CONTENTS
!     ------------------------------------------------------------------
!      1  (ERI)    2e-Integrals File for 1/r12 interaction (DONTW=F)
!      2  (CGM)    Output File for CGM routines (ICGMETHOD=1,2,3)
!      3  (GCF)    GAMMA, Coefficient Matrix (C), Energies (E), FMIUG0
!      4  (BFST)   Basis Function Symbol Table
!      5           Input File
!      6           Output File
!      7  (WFN)    WFN File for AIMPAC Program
!      8  (GCFe)   GCF file corresponding to the minimum energy
!      9  (APSG)   Coefficient Matrix (C) and APSG Wavefunction of PNOF5
!     10  (FRAG)   FRAG file containing the fragment information
!     11  (CGGRAD) Output File for Optimization
!     12  (CJK)    Output File for CJ12 and CK12 (NOUTCJK)
!     13  (CND)    Output File for Non-Dynamic CK12 and MP2 amplitudes
!     14  (2DM)    Output File for atomic 2RDM (NOUTRDM)
!     15  (1DM)    Output File for atomic 1RDM (NOUTRDM)
!     16  (N2DM)   Output File for Record Number if NSQT=1
!     17  (MLD)    Output File for MOLDEN Program
!     18  (XYZ)    Output File with Geometries for MOLDEN Program
!     19  (FCHK)   Formatted checkpoint file for visualization softwares
!     20  (IRAF)   Direct File used for DIIS data in RHFCL and RHFOP
!
!     30 (DYN.xyz)  Formatted xyz trajectory file
!     31 (DYNl.xyz) Formatted xyz last trajectory file (use for restart)
!     32 (DYNen)    Formatted File with energy information
!     33 ()
!     34 (EPOTs)    Formatted File with ngcf potential energies
!
!     50 (BASIS_FILE) Opened in Subroutine ATOMS <- MOLECULE <- START
!-----------------------------------------------------------------------
!                    Open general working files
!-----------------------------------------------------------------------
      OPEN(2,FILE='CGM',STATUS='UNKNOWN',FORM='FORMATTED',              &
      ACCESS='SEQUENTIAL')
      IF(IRUNTYP<5)OPEN(3,FILE='GCF',STATUS='UNKNOWN',                  &
       FORM='FORMATTED',ACCESS='SEQUENTIAL')
      OPEN(4,FILE='BFST',STATUS='UNKNOWN',FORM='FORMATTED',             &
      ACCESS='SEQUENTIAL')
      IF(IAIMPAC==1)OPEN(7,FILE='WFN',STATUS='UNKNOWN',                 &
        FORM='FORMATTED',ACCESS='SEQUENTIAL')
      IF(IRUNTYP<5)OPEN(8,FILE='GCFe',STATUS='UNKNOWN',                 &
       FORM='FORMATTED',ACCESS='SEQUENTIAL')
      IF(APSG)OPEN(9,FILE='APSG' ,STATUS='UNKNOWN',                     &
      FORM='FORMATTED',ACCESS='SEQUENTIAL')
      IF(ICOEF==3)OPEN(10,FILE='FRAG' ,STATUS='OLD',                    &
       FORM='FORMATTED',ACCESS='SEQUENTIAL')
      IF(IRUNTYP==3.or.IRUNTYP==4)OPEN(11,FILE='CGGRAD',                &
      STATUS='UNKNOWN',FORM='FORMATTED',ACCESS='SEQUENTIAL')
      IF(NOUTCJK==1)OPEN(12,FILE='CJK',STATUS='UNKNOWN',                &
         FORM='UNFORMATTED')
      IF(NOUTTijab==1)OPEN(13,FILE='CND',STATUS='UNKNOWN',              &
           FORM='UNFORMATTED')
      IF(NOUTRDM==1.or.NOUTRDM==3)THEN
        if(NSQT==0)then
          OPEN(15,FILE='1DM',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',      &
          FORM='FORMATTED')
        else if(NSQT==1)then
          OPEN(15,FILE='1DM',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',      &
          FORM='FORMATTED')
        end if
      END IF
      IF(NOUTRDM==2.or.NOUTRDM==3)THEN
        if(NSQT==0)then
          OPEN(14,FILE='2DM',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',      &
          FORM='FORMATTED')
        else if(NSQT==1)then
          OPEN(14,FILE='2DM',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',      &
          FORM='UNFORMATTED')
          OPEN(16,FILE='N2DM',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',     &
          FORM='FORMATTED')
        end if
      END IF
      IF(MOLDEN==1)THEN
        OPEN(17,FILE='MLD',STATUS='UNKNOWN',FORM='FORMATTED',           &
        ACCESS='SEQUENTIAL')
        OPEN(18,FILE='XYZ',STATUS='UNKNOWN',FORM='FORMATTED',           &
        ACCESS='SEQUENTIAL')
      ENDIF
      IF(IFCHK==1)OPEN(19,FILE='FCHK',STATUS='UNKNOWN',                 &
       FORM='FORMATTED',ACCESS='SEQUENTIAL')
!-----------------------------------------------------------------------
      RETURN
      END

! OUTPUTBASIC
      SUBROUTINE OUTPUTBASIC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPFILE_NO1PT2/NO1PT2,NEX
!-----------------------------------------------------------------------
!     NE: Number of Electrons
!     MUL: State Multiplicity
!     NCO: Number of doubly filled molecular orbitals in HF (CLOSED)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NCO = NB
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NSOC: Number of strongly singly occupied MOs
!     NTWOPAR: 1 => Two-particle case
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NTWOPAR = 0
!     Spin-compensated
      IF( NB==(NE+MUL-1)/2 .and. NA==(NE-MUL+1)/2 )THEN
        NSOC=0
        MSpin=0
!      Two-particle case
      IF(NB==1)NTWOPAR=1
      IF(NA/=NB)THEN
        WRITE(6,1)NA,NB
        CALL ABRT
      ENDIF
!     Spin-uncompensated [ NA = NB+MUL-1 ]
      ELSE IF( NB==(NE-MUL+1)/2 .and. NA==(NE+MUL-1)/2 )THEN
        NSOC=NA-NB
        if(HighSpin)then
          MSpin=NSOC
        else
        MSpin=0
        endif
        IF(NB==0.and.NA==2)NTWOPAR=1
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NO1: Number of doubly filled NOs with occupancies equal to one
!     NO1PT2: Number of doubly filled NOs in perturbative calculations
!     Stop Program if NO1 > NB or NO1PT2 >= NB
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NO1>NB)THEN
        WRITE(6,2)
        CALL ABRT
      ENDIF
      IF(NO1PT2>=NB)THEN
        write(6,3)
        CALL ABRT
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
!-----------------------------------------------------------------------
!     Format definitions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    1 FORMAT(/'Spin compensated but NA (',I3,') not equal NB (',I3,')'  &
      /' JOB ABANDONED'/)
    2 FORMAT(/,' Error: NO1 > NB (doubly filled NOs) -> Stop Program')
    3 FORMAT(/,' Error: NO1PT2 >= NB (doubly filled NOs)->Stop Program')
!-----------------------------------------------------------------------
      END

! SETNO1
      SUBROUTINE SETNO1(IAN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      INTEGER,DIMENSION(NATOMS)::IAN
!-----------------------------------------------------------------------
!     Determine NO1 (NOs with ONs equal to 1, according to the true
!     nuclear charges if NO1=-1.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NO1= 0
      DO I=1,NATOMS
        NUCZ = IAN(I)
        IF( 1<=NUCZ.and.NUCZ<=  2)NO1i =  0          ! H-He
        IF( 3<=NUCZ.and.NUCZ<= 10)NO1i =  1          ! Li-Ne
        IF(11<=NUCZ.and.NUCZ<= 18)NO1i =  5          ! Na-Ar
        IF(19<=NUCZ.and.NUCZ<= 36)NO1i =  9          ! K-Kr
        IF(37<=NUCZ.and.NUCZ<= 49)NO1i = 18          ! Rb-In
        IF(50<=NUCZ.and.NUCZ<= 54)NO1i = 23          ! Sn-Xe
        IF(55<=NUCZ.and.NUCZ<= 71)NO1i = 27          ! Cs-Lu
        IF(72<=NUCZ.and.NUCZ<= 81)NO1i = 30          ! Hf-Tl
        IF(82<=NUCZ.and.NUCZ<= 86)NO1i = 39          ! Pb-Rn
        IF(87<=NUCZ.and.NUCZ<=109)NO1i = 43          ! Fr-Mt
        NO1 = NO1 + NO1i
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! SETORBSPACE
      SUBROUTINE SETORBSPACE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INPFILE_NO1PT2/NO1PT2,NEX
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,10)
      WRITE(6,11)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     If NBF < NE -> More orbitals have to be excluded in NO1 or NO1PT2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NBF<NE)THEN
        NDIF = NE-NBF
        IF(NO1<NDIF)THEN
          NO1 = NDIF
          WRITE(6,1)
        ENDIF
        IF(NO1PT2<NDIF)THEN
          NO1PT2 = NDIF
          WRITE(6,2)
        ENDIF
        IF(NO1==NBF)THEN
          WRITE(6,3)
          CALL ABRT
        ENDIF
        IF(NO1PT2==NBF)THEN
          WRITE(6,4)
          CALL ABRT
        ENDIF
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Frozen orbitals in perturbative calculations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NO1PT2==-1)NO1PT2 = NO1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                         Set Orbital Space
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NO1:  Number of inactive doubly occupied orbitals (OCC=1)
!     NDOC: Number of strongly doubly occupied MOs
!     NSOC: Number of strongly singly occupied MOs
!     NDNS: Number of strongly occupied MOs (NDNS=NDOC+NSOC)
!     NCWO: Number of coupled weakly occ. MOs per strongly doubly occ.
!     NCWO*NDOC: Active orbitals in the virtual subspace
!     NO0:  Empty orbitals  (OCC=0)
!     NVIR: Number of weakly occupied MOs + empty MOs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!           NO1 | NDOC  + NSOC  |   NCWO*NDOC + NO0  = NBF
!           NO1 |      NDNS     |          NVIR      = NBF
!               | -NAC- |       |  -   NAC  - |
!                      NB      NA            NBF5
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                               !  CLOSED (NB=NA=NCO,NSOC=0)
      NDOC = NB - NO1           !  NDOC = NCO - NO1, NO1 <= NCO
      NDNS = NDOC + NSOC        !  NDNS = NDOC
      NA   = NO1 + NDNS         !  NA = NB = NCO
      NVIR = NBF - NA           !  NBF - NCO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NCWO: Number of coupled weakly occ. MOs per strongly doubly occ.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NDNS/=0)THEN
        if(NDOC>0)then
          IF(NCWO/=1)THEN                      ! Extended PNOF (NCWO>1)
!- - - - - - - - - - - - - - - - - -
            if(NCWO<-1)then
              write(6,5)NCWO
              CALL ABRT
            else if(NCWO==-1)then
              NCWO = NVIR/NDOC
            else if(NCWO>NVIR/NDOC)then
              write(6,6)NCWO
              NCWO = NVIR/NDOC
            endif
!- - - - - - - - - - - - - - - - - -
          ELSE                                 ! perfect pairing (NCWO=1)
!- - - - - - - - - - - - - - - - - -
            write(6,7)
!- - - - - - - - - - - - - - - - - -
          ENDIF
        else
          NCWO = 0
        end if
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NAC: Dimension of the active natural orbital subspace
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NAC = NDOC * ( 1 + NCWO )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NBF5: Occupied Orbitals (ON /= 0), NBFT5, NSQ5
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NBF5 = NO1 + NAC + NSOC                ! NBF5 = NA  + NDOC*NCWO
      IF(NBF5>NBF)NBF5 = NBF
      NBFT5 = NBF5*(NBF5+1)/2
      NSQ5 = NBF5*NBF5
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NO0: Empty orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NO0 = NBF - NBF5
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,8)NO1,NDOC,NSOC,NCWO,NAC,NO0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
!-----------------------------------------------------------------------
!     Format definitions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    1 FORMAT(/1X,'Note NBF < NE: More orbitals have to be excluded.',   &
      1X,'NO1 has been set equal to NE-NBF')
    2 FORMAT(/1X,'Note NBF < NE: More orbitals have to be excluded.',   &
      1X,'NO1PT2 has been set equal to NE-NBF')
    3 FORMAT(/1X,'Stop: all orbitals are full occupied (NO1=NBF)')
    4 FORMAT(/1X,'Stop: all orbitals are full occupied (NO1PT2=NBF)')
    5 FORMAT(/1X,'Stop Program: Incorrect number of NCWO =',I5)
    6 FORMAT(/1X,'Your NCWO =',I5,' exceeds the maximum possible value')
    7 FORMAT(/1X,'You are doing a perfect pairing calculation: NCWO=1')
    8 FORMAT(/1X,'Inactive Doubly occupied orbitals up to NO1  =',I5,   &
      /1X,'No. considered Strongly Doubly occupied MOs  =',I5,          &
      /1X,'No. considered Strongly Singly occupied MOs  =',I5,          &
      /1X,'No. of Weakly occ. per St. Doubly occ.  MOs  =',I5,          &
      /1X,'Dimension of the Active Nat. Orb. subspace   =',I5,          &
      /1X,'Secondary Empty orbitals                     =',I5)
   10 FORMAT(/72('-'))
   11 FORMAT(/1X,'Orbital Space',/1X,'-------------')
!-----------------------------------------------------------------------
      END

! POINTERS
      SUBROUTINE POINTERS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!-----------------------------------------------------------------------
!     Define Pointers of the USER array
!-----------------------------------------------------------------------
      N1  = 1                    ! USER( N1) = RO(NBF5)
      N2  = N1  + NBF5           ! USER( N2) = CJ12(NBF5,NBF5)
      N3  = N2  + NSQ5           ! USER( N3) = CK12(NBF5,NBF5)
      N4  = N3  + NSQ5           ! USER( N4) = DR(NBF5,NBF5)
      N5  = N4  + NSQ5           ! USER( N5) = DCJ12r(NBF5,NBF5,NBF5)
      N6  = N5  + NSQ5*NBF5      ! USER( N6) = DCK12r(NBF5,NBF5,NBF5)
      N7  = N6  + NSQ5*NBF5      ! USER( N7) = QD(NBF,NBF,NBF)
      N8  = N7  + NBF*NSQ        ! USER( N8) = HCORE(NBF5)
      N9  = N8  + NBF5           ! USER( N9) = QJ(NBFT5)
      N10 = N9  + NBFT5          ! USER(N10) = QK(NBFT5)
      N11 = N10 + NBFT5          ! USER(N11) = DIPN(3)
      N12 = N11 + 3              ! USER(N12) = ADIPx(NSQ)
      N13 = N12 + NSQ            ! USER(N13) = ADIPy(NSQ)
      N14 = N13 + NSQ            ! USER(N14) = ADIPz(NSQ)
      N15 = N14 + NSQ            ! USER(N15) = DIPx(NSQ5)
      N16 = N15 + NSQ5           ! USER(N16) = DIPy(NSQ5)
      N17 = N16 + NSQ5           ! USER(N17) = DIPz(NSQ5)
      N18 = N17 + NSQ5           ! USER(N18) = QUADN(6)
      N19 = N18 + 6              ! USER(N19) = AQUADxx(NSQ)
      N20 = N19 + NSQ            ! USER(N20) = AQUADyy(NSQ)
      N21 = N20 + NSQ            ! USER(N21) = AQUADzz(NSQ)
      N22 = N21 + NSQ            ! USER(N22) = AQUADxy(NSQ)
      N23 = N22 + NSQ            ! USER(N23) = AQUADxz(NSQ)
      N24 = N23 + NSQ            ! USER(N24) = AQUADyz(NSQ)
      N25 = N24 + NSQ            ! USER(N25) = QUADxx(NSQ5)
      N26 = N25 + NSQ5           ! USER(N26) = QUADyy(NSQ5)
      N27 = N26 + NSQ5           ! USER(N27) = QUADzz(NSQ5)
      N28 = N27 + NSQ5           ! USER(N28) = QUADxy(NSQ5)
      N29 = N28 + NSQ5           ! USER(N29) = QUADxz(NSQ5)
      N30 = N29 + NSQ5           ! USER(N30) = QUADyz(NSQ5)
      N31 = N30 + NSQ5           ! USER(N31) = OCTUN(10)
      N32 = N31 + 10             ! USER(N32) = AOCTxxx(NSQ)
      N33 = N32 + NSQ            ! USER(N33) = AOCTyyy(NSQ)
      N34 = N33 + NSQ            ! USER(N34) = AOCTzzz(NSQ)
      N35 = N34 + NSQ            ! USER(N35) = AOCTxxy(NSQ)
      N36 = N35 + NSQ            ! USER(N36) = AOCTxxz(NSQ)
      N37 = N36 + NSQ            ! USER(N37) = AOCTxyy(NSQ)
      N38 = N37 + NSQ            ! USER(N38) = AOCTyyz(NSQ)
      N39 = N38 + NSQ            ! USER(N39) = AOCTxzz(NSQ)
      N40 = N39 + NSQ            ! USER(N40) = AOCTyzz(NSQ)
      N41 = N40 + NSQ            ! USER(N41) = AOCTxyz(NSQ)
      N42 = N41 + NSQ            ! USER(N42) = OCTXXX(NSQ5)
      N43 = N42 + NSQ5           ! USER(N43) = OCTYYY(NSQ5)
      N44 = N43 + NSQ5           ! USER(N44) = OCTZZZ(NSQ5)
      N45 = N44 + NSQ5           ! USER(N45) = OCTXXY(NSQ5)
      N46 = N45 + NSQ5           ! USER(N46) = OCTXXZ(NSQ5)
      N47 = N46 + NSQ5           ! USER(N47) = OCTXYY(NSQ5)
      N48 = N47 + NSQ5           ! USER(N48) = OCTYYZ(NSQ5)
      N49 = N48 + NSQ5           ! USER(N49) = OCTXZZ(NSQ5)
      N50 = N49 + NSQ5           ! USER(N50) = OCTYZZ(NSQ5)
      N51 = N50 + NSQ5           ! USER(N51) = OCTXYZ(NSQ5)
      NUSER = N51 - N1 + NSQ5
!-----------------------------------------------------------------------
      RETURN
      END

!======================================================================!
!                                                                      !
!     ENERGRAD: Calculate Energy and Gradient (IRUNTYP=1,2)            !
!                                                                      !
!======================================================================!

! ENERGRAD
      SUBROUTINE ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,    &
                          NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,   &
                          KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,   &
                          C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,     &
                          DIPS,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,NOPTCG,  &
                          IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER :: NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI
      INTEGER :: IRUNTYP,NOPTCG,IPRINTOPT
      LOGICAL RESTART,SMCD
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/USELIBCINT/ILIBCINT
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
#include "mpip.h"
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      INTEGER,DIMENSION(NAT):: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: C1,C2,EX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3*NAT) :: GRADS
      DOUBLE PRECISION,DIMENSION(3) :: DIPS
      INTEGER :: SIZE_ENV,ATM(6,NAT),NBAS,BAS(8,NBAS),IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      INTEGER(8),ALLOCATABLE,DIMENSION(:) :: IBUF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: H,S,EiHF,CHF,BUF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DIPN,QUADN,OCTUN
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DQOInt,AUX
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: AHCORE,OVERLAP
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: XINTS
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: BUFaux
!-----------------------------------------------------------------------
!     Allocate necessary arrays for 1e- and 2e- integrals + Guess
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NSQ = NBF*NBF
      NBFT = (NBF*NBF+NBF)/2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     1e Integrals: H, S ; Initial Orbitals: CHF, EiHF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(H(NBFT),S(NBFT),EiHF(NBF),CHF(NSQ))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     2e- Integrals (BUF,IBUF)
!     NINTEGtm: Maximum numbers of distinct two-electron integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IDONTW==1)THEN
        NINTEGtm = 0
        NINTEGAUXtm = 0
        if(IERITYP==1)then
          NINTEGtm = NBF*(NBF+1)*(NBF*NBF+NBF+2)/8
        else if(IERITYP==2)then
          NINTEGAUXtm = NBF*(NBF+1)/2*NBFaux
        else if(IERITYP==3)then
          NINTEGtm = NBF*(NBF+1)*(NBF*NBF+NBF+2)/8
          NINTEGAUXtm = NBF*(NBF+1)/2*NBFaux
        end if
      ELSE
        NINTEGtm = NINTEG
        NINTEGAUXtm = 0
      ENDIF
      ALLOCATE(BUF(NINTEGtm),IBUF(NINTEGtm),BUFaux(NINTEGAUXtm))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Integrals & Guess & RHF
!     Note: NINTEGt<NINTEGtm due to the CUTOFF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NSH2 = (NSHELL*NSHELL+NSHELL)/2
      ALLOCATE(XINTS(NSH2))
      CALL GuessHJKRHF(KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,EX,CS,CP,  &
                       CD,CF,CG,CH,CI,NPRIMI,Cxyz,H,S,EiHF,CHF,BUF,     &
                       IBUF,BUFaux,NSHELL,NAT,NBF,NSQ,NBFT,NINTEGtm,    &
                       NINTEGAUXtm,NINTEGt,NREC,XINTS,NSH2,IDONTW,      &
                       INPUTC,IPRINTOPT,ZAN,SIZE_ENV,ENV,ATM,NBAS,BAS,  &
                       IGTYP)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Preparing for RunNOF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IDONTW==0)REWIND(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Square Matrices AHCORE, OVERLAP
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(AHCORE(NBF,NBF),OVERLAP(NBF,NBF))
      CALL CPYTSQ(H,AHCORE,NBF)
      CALL CPYTSQ(S,OVERLAP,NBF)
      DEALLOCATE(H,S)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Multipole Moment Integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(DIPN(3),QUADN(6),OCTUN(10))
      CALL DQONuclear(DIPN,QUADN,OCTUN,Cxyz,ZAN,NAT)
!
      IF(IEMOM==1)THEN
        NVAL=3
      ELSE IF(IEMOM==2)THEN
        NVAL=3+6
      ELSE IF(IEMOM==3)THEN
        NVAL=3+6+10
      END IF
      ALLOCATE(DQOInt(NVAL*NBFT),AUX(NVAL*784))
      if(ILIBCINT==0)then
!HONDO CALL PRCALC(DQOInt,AUX,NVAL,NBFT,EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,&
!HONDO             KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,Cxyz)
      else if(ILIBCINT==1)then
       CALL PRCALCl(DQOInt,NVAL,NBFT,SIZE_ENV,ENV,NAT,ATM,NSHELL,       &
                    NBAS,BAS,IGTYP)
      end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     PNOF Calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL RunNOF(NAT,NBF,NBFT,NSHELL,NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,    &
                  KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP, &
                  C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,AHCORE,OVERLAP,CHF,EiHF,&
                  DIPN,QUADN,OCTUN,NVAL,DQOInt,NINTEG,NREC,IBUF,BUF,    &
                  BUFaux,NINTEGt,NINTEGAUXtm,IDONTW,GRADS,IRUNTYP,DIPS, &
                  XINTS,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,IPRINTOPT)
      DEALLOCATE(AHCORE,OVERLAP,CHF,EiHF,DIPN,QUADN,OCTUN,DQOInt,AUX)
      DEALLOCATE(IBUF,BUF,BUFaux,XINTS)
      NOPTCGMPI = NOPTCG
#ifdef MPI
      CALL DEACTIVATEERIs
      DO I=1,NPROCS-1
        K=0
        CALL MPI_SEND(K,1,MPI_INTEGER,I,I,MPI_COMM_WORLD,IERR)
        CALL MPI_SEND(NOPTCGMPI,1,MPI_INTEGER,I,I,MPI_COMM_WORLD,IERR)
      ENDDO
#endif
!-----------------------------------------------------------------------
      RETURN
      END

!======================================================================!
