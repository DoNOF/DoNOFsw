!======================================================================!
!                                                                      !
!               I N T E G R A L   S U B R O U T I N E S                !
!                                                                      !
!======================================================================!

! OneElecInt                                           
      SUBROUTINE OneElecInt(Cxyz,H,S,TKIN,DInteg,NBFT,IPRINTOPT,        &
                            EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,KSTART,      &
                            KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,ZAN,  &
                            NAT,IECP,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/USELIBCINT/ILIBCINT      
      LOGICAL EFIELDL
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
!      
      INTEGER :: NAT, IECP
      INTEGER :: NBFT,IPRINTOPT,NPRIMI,NSHELL
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NBFT) :: H,S,TKIN    
      DOUBLE PRECISION,DIMENSION(3*NBFT) :: DInteg
      INTEGER :: CR, CM, timestartoneE, timefinishoneE
      DOUBLE PRECISION :: RATE
!     LIBCINT
      INTEGER :: SIZE_ENV,NBAS
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NAT), BAS(8,NBAS)
      INTEGER :: IGTYP
!-----------------------------------------------------------------------
      IF(IPRINTOPT==1)WRITE(6,'(/1X,A13/,1X,13(1H-))')'1e- integrals'
!-----------------------------------------------------------------------
!     Initialization for system_clock
!-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT_RATE=CR)
      CALL SYSTEM_CLOCK(COUNT_MAX=CM)
      RATE = REAL(CR)
      CALL SYSTEM_CLOCK(timestartoneE)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
!     H, S & T integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ILIBCINT==0) THEN
!HONDO  CALL HSandT(H,S,TKIN,NBFT,EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,       &
!HONDO              KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,       &
!HONDO              ZAN,Cxyz)
      ELSE IF(ILIBCINT==1) THEN
        CALL HSandTl(H,S,TKIN,NBFT,SIZE_ENV,ENV,NAT,ATM,NSHELL,         &
                     NBAS,BAS,IGTYP)
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ECP integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF( 1<=IECP .and. IECP<=3 )THEN
        IF(ILIBCINT==0) THEN
!HONDO    CALL ECPINT(H,NBFT,EX,CS,CP,CD,CF,CG,NPRIMI,KSTART,KATOM,KNG, &
!HONDO                KLOC,KMIN,KMAX,NSHELL,Cxyz)
        ELSE IF(ILIBCINT==1) THEN
          CALL ECPINTl(H,NBFT,SIZE_ENV,ENV,NAT,ATM,NSHELL,              &
                       NBAS,BAS,IGTYP)
        END IF
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Add Electric Field Contribution to 1e Hamiltonian
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(EFIELDL)THEN
       CALL DipInt(DInteg,NBFT,EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,          &
                   KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,        &
                   Cxyz,NAT,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP)
       IF(IPRINTOPT==1)WRITE(6,'(1X,A16,3F10.5)')'Electric Field =',    &
                       EFX,EFY,EFZ
       CALL ElecFieldInt(H,DInteg,NBFT)
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL SYSTEM_CLOCK(timefinishoneE)
      DeltaToneE = (timefinishoneE - timestartoneE)/RATE
      IF(IPRINTOPT==1)                                                  &
       WRITE(6,'(1X,A22,F10.2)')'Time to do integrals =',DeltaToneE
!-----------------------------------------------------------------------
      RETURN                                                            
      END SUBROUTINE

! DipInt                                           
      SUBROUTINE DipInt(DInteg,NBFT,EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,     &
                        KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,   &
                        Cxyz,NAT,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/USELIBCINT/ILIBCINT            
      COMMON/ELPROP/IEMOM      
      COMMON/TRANSF/XP,YP,ZP
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(3*NBFT) :: DInteg
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AUX
!     LIBCINT
      INTEGER :: IGTYP
      INTEGER :: SIZE_ENV,NBAS
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NAT), BAS(8,NBAS)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Dipole integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      XP = 0.0d0                                                           
      YP = 0.0d0
      ZP = 0.0d0                                                        
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IEMSV = IEMOM                                                  
      IEMOM = 1    
      ALLOCATE(AUX(3*784))
      IF (ILIBCINT==0) THEN
!HONDO  CALL PRCALC(DInteg,AUX,3,NBFT,EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,   &
!HONDO              KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,Cxyz)
      ELSE IF(ILIBCINT==1) THEN
        CALL PRCALCl(DInteg,3,NBFT,SIZE_ENV,ENV,NAT,ATM,NSHELL,         &
                     NBAS,BAS,IGTYP)
      END IF
      IEMOM = IEMSV                                                  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DEALLOCATE(AUX)                                                                                  
      RETURN                                                            
      END                                                               

! ElecFieldInt
      SUBROUTINE ElecFieldInt(H,DInteg,NBFT)   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EFIELDL
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      DOUBLE PRECISION,DIMENSION(NBFT) :: H
      DOUBLE PRECISION,DIMENSION(3*NBFT) :: DInteg
!-----------------------------------------------------------------------
      IF(EFX/=0.0d0)CALL DAXPY(NBFT,EFX,DInteg(1       ),1,H,1)
      IF(EFY/=0.0d0)CALL DAXPY(NBFT,EFY,DInteg(1+  NBFT),1,H,1)
      IF(EFZ/=0.0d0)CALL DAXPY(NBFT,EFZ,DInteg(1+2*NBFT),1,H,1)
!----------------------------------------------------------------------- 
      RETURN                                                            
      END                                                               

! JandKl
      SUBROUTINE JandKl(BUFP2,IX2,NINTEGtm,NINTEGt,NREC,XINTS,NSH2,     &
                        IDONTW,IPRINTOPT,EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,&
                        KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,   &
                        Cxyz,NAT,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/INTFIL/NINTMX           
      COMMON/INTOPT/CUTOFF,ISCHWZ,IECP,NECP
!      
      LOGICAL SCHWRZ
      INTEGER :: NINTEGtm,NINTEGt,NREC,NSH2,IDONTW,IPRINTOPT,IGTYP
      INTEGER :: NPRIMI,NSHELL,NAT
      INTEGER(8),DIMENSION(NINTEGtm) :: IX2
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NINTEGtm) :: BUFP2
      DOUBLE PRECISION,DIMENSION(NSH2) :: XINTS
!
      INTEGER(8),ALLOCATABLE,DIMENSION(:) :: IX
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: BUFP,GHONDOl
      INTEGER :: CR, CM, timestarttwoE, timefinishtwoE
      DOUBLE PRECISION :: RATE
!     LIBCINT
      INTEGER :: SIZE_ENV, NBAS
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6, NAT), BAS(8, NBAS)
!-----------------------------------------------------------------------
!     Initialization for system_clock
!-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT_RATE=CR)
      CALL SYSTEM_CLOCK(COUNT_MAX=CM)
      RATE = REAL(CR)
      CALL SYSTEM_CLOCK(timestarttwoE)
!----------------------------------------------------------------------- 
!     Driver for 2e integrals 
!----------------------------------------------------------------------
      CALL BASCHK(LMAXIMA,KTYPE,NSHELL)
      MAXG = 4**4                                                                                  
      IF(LMAXIMA==2)MAXG =  6**4 
      IF(LMAXIMA==3)MAXG = 10**4 
      IF(LMAXIMA==4)MAXG = 15**4 
      IF(LMAXIMA==5)MAXG = 21**4 
      IF(LMAXIMA==6)MAXG = 28**4
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Debut
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL Debutl(IDONTW,IPRINTOPT,KATOM,NSHELL,Cxyz,NINTMX,NAT)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Schwarz inequality
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(GHONDOl(MAXG))
      SCHWRZ = ISCHWZ>0 
      IF(SCHWRZ)CALL ExchangeIntl(XINTS,NSH2,SIZE_ENV,ENV,NAT,ATM,      &
                                  NSHELL,NBAS,BAS,IGTYP)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     2e integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(BUFP(NINTMX),IX(NINTMX))
      NREC = 1
      CALL TwoERIl(SCHWRZ,NINTEGtm,NINTEGt,NSCHWZ,BUFP,IX,BUFP2,        &
                   IX2,XINTS,NSH2,IDONTW,IPRINTOPT,NSHELL,Cxyz,NAT,     &
                   NINTMX,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,NREC)
      DEALLOCATE(BUFP,IX,GHONDOl)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL SYSTEM_CLOCK(timefinishtwoE)
      DeltaTtwoE = (timefinishtwoE - timestarttwoE)/RATE
      IF(IPRINTOPT==1)                                                  &
       WRITE(6,'(1X,A22,F10.2)')'Time to do integrals =',DeltaTtwoE
!-----------------------------------------------------------------------
      RETURN                                                            
      END

! JandKauxl
      SUBROUTINE JandKauxl(BUFP2,NINTEGtm,IDONTW,IPRINTOPT,NBF,         &
                           EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,KSTART,KATOM, &
                           KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,Cxyz,NAT,    &
                           SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      COMMON/INTFIL/NINTMX
      INTEGER :: NINTEGtm,IDONTW,IPRINTOPT,NBF,NSHELL,NAT,IGTYP
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NINTEGtm) :: BUFP2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: Gaux
      INTEGER :: CR, CM, timestarttwoE, timefinishtwoE
      DOUBLE PRECISION :: RATE
      LOGICAL SMCD
!     LIBCINT
      INTEGER :: SIZE_ENV,NBAS
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NAT), BAS(8,NBAS)

      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
!-----------------------------------------------------------------------
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
!-----------------------------------------------------------------------
!     Initialization for system_clock
!-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT_RATE=CR)
      CALL SYSTEM_CLOCK(COUNT_MAX=CM)
      RATE = REAL(CR)
      CALL SYSTEM_CLOCK(timestarttwoE)
!-----------------------------------------------------------------------
!     Driver for 2e integrals
!-----------------------------------------------------------------------
      CALL BASCHK(LMAXIMA,KTYPE,NSHELL)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Debut
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL Debutl(IDONTW,IPRINTOPT,KATOM,NSHELL,Cxyz,NINTMX,NAT)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     2e integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL BASCHK(LAUXMAXIMA,KTYPEaux,NSHELLaux)
      MAXORI = (LMAXIMA+1)*(LMAXIMA+2)/2
      MAXORIAUX = (LAUXMAXIMA+1)*(LAUXMAXIMA+2)/2
      MAXG = MAX(MAXORI*MAXORI*MAXORIAUX,MAXORIAUX*MAXORIAUX)
!
      ALLOCATE(Gaux(MAXG))
      IF(.NOT.SMCD) THEN
       CALL AuxERIl(NINTEGtm,BUFP2,NBF,IPRINTOPT,NSHELL,NAT,            &
                    SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP)
      ELSE IF(SMCD) THEN
       CALL AuxERIModCholl(NINTEGtm,BUFP2,NBF,IPRINTOPT,NSHELL,NAT,     &
                           SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP)
      END IF
      DEALLOCATE(Gaux)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL SYSTEM_CLOCK(timefinishtwoE)
      DeltaTtwoE = (timefinishtwoE - timestarttwoE)/RATE
      IF(IPRINTOPT==1)                                                  &
       WRITE(6,'(1X,A22,F10.2)')'Time to do integrals =',DeltaTtwoE
      CALL FLUSH(6)
!-----------------------------------------------------------------------
      RETURN
      END

! Debutl
      SUBROUTINE Debutl(IDONTW,IPRINTOPT,KATOM,NSHELL,Cxyz,NINTMX,NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER :: NINTMX, NAT
      INTEGER :: IDONTW,IPRINTOPT,NSHELL,I,NBYTES,ICC
      INTEGER,DIMENSION(NSHELL) :: KATOM
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NSHELL,3) :: CO
!-----------------------------------------------------------------------
!     ERI Initializations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IPRINTOPT==1)WRITE(6,10)
      NBYTES = 16
      IF(IDONTW==1)THEN
       IF(IPRINTOPT==1)WRITE(6,20)
      ELSE
       IF(IPRINTOPT==1)WRITE(6,30)NINTMX,NBYTES
       REWIND(1)
      END IF
      DO I=1,NSHELL
       ICC = KATOM(I)
       CO(I,1)= Cxyz(1,ICC)
       CO(I,2)= Cxyz(2,ICC)
       CO(I,3)= Cxyz(3,ICC)
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Format Statements
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   10 FORMAT(/1X,'2e- integrals'/1X,13(1H-))
   20 FORMAT(' DONTW option skips storage 2e- integrals on Unit 1')
   30 FORMAT(' Storing',I8,' integrals/record on disk, using',I3,       &
             ' Bytes/integral')
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE

! OTTOINTEGR
      SUBROUTINE OTTOINTEGR(I,J,K,L,NIJ,NKL,XJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
!-----------------------------------------------------------------------
!     TRANSFORM INTO INTSCF INTEGRALS
!-----------------------------------------------------------------------
      IF(I==J)XJ=XJ*2.0
      IF(K==L)XJ=XJ*2.0
      IF(NIJ==NKL)XJ=XJ*2.0
!-----------------------------------------------------------------------
      RETURN
      END

!----------------------------------------------------------------------!
!                                                                      !
!                    M P I   S U B R O U T I N E S                     !
!                                                                      !
!      2013 Four-index transformation of the electron repulsion        !
!           integrals was parallelized by Eduard Matito                !
!                                                                      !
!          04/26/2013 module developed by Eduard Matito                !
!          06/28/2017 module modified by Ion Mitxelena                 !
!          11/15/2020 module modified by Juan Felipe Huan Lew Yee      !                    
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                      !
!   INIMPI: Initializion of the parallel or serial version of the code !
!   Abortx: Abort if OPT-MPI fatal error occurs.                       !
!   SlaveDvr: Slave driver                                             !
!   SlvInt: Calculation of the integrals in the slave nodes.           !
!   SLVHSTJ: HSTARJ subroutine adapted for Slaves                      !
!   SLVHSTK: HSTARK subroutine adapted for Slaves                      !
!   SLVFORM2JK: FORM2JK subroutine adapted for Slaves                  !
!   SLVFORMJK: FORMJK adapted for Slaves                               !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! INIMPI
      Subroutine INIMPI()
      Implicit None
#include "mpip.h"
#ifdef MPI
#ifdef _OPENMP
      !$OMP PARALLEL
      Write(6,'(A)')' This is the MPI/OpenMP version of the code'
      !$OMP END PARALLEL
#else
      Write(6,'(A)')' This is the MPI version of the code'
#endif
      CALL MPI_INIT (IERR)
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,MY_ID,IERR)
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,NPROCS,IERR)
      LMASTR=(MY_ID.EQ.MASTER)
      IF (LMASTR) Write(6,'(A,I4,A,I8)') 'Master: ',MY_ID
     
      IF (NPROCS<2) THEN
        Call Abortx ("OPT-MPI Fatal error! Use at least 2 cores!")
      ELSE IF (.NOT.LMASTR) THEN
        Call SlaveDvr()
      ENDIF

#else
#ifdef _OPENMP
      !$OMP PARALLEL
      Write(6,'(A)')' This is the OpenMP version of the code'
      !$OMP END PARALLEL
#else
      Write(6,'(A)')' This is the Serial version of the code'
#endif
      NPROCS=1
      MY_ID=0
      LMASTR=.True.
#endif
      End

! Abortx
      Subroutine Abortx(String)
      Implicit None
#include "mpip.h"
      Character*(*) String
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
#endif     
      Write(6,'(A,A)') 'Error... Aborting!',String
      Stop 'Error: Aborting!'
      End

#ifdef MPI

! SlaveDvr
      Subroutine SlaveDvr()
      Implicit None
#include "mpip.h"
      INTEGER NINTCR,INTTYPE,ID,K,IHUB
      COMMON/USEHUBBARD/IHUB
      ID=MY_ID
      K=MASTER
   10 CALL MPI_BCAST(INTTYPE,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
!HUB      CALL MPI_RECV(IHUB,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,STATUS,IERR)
      IF(INTTYPE==1) THEN
        CALL MPI_RECV(NINTCR,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,STATUS,IERR)
        Call SlvInt(NINTCR)
        GOTO 10
      ELSE IF(INTTYPE==2) THEN
        CALL MPI_RECV(NINTCR,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,STATUS,IERR)
        Call SlvIntRI(NINTCR)
        GOTO 10
      ELSE IF(INTTYPE==3) THEN
        RETURN
      END IF
      GOTO 10
      End

! SlvInt
      Subroutine SlvInt(NINTCR)
      Implicit None
#include "mpip.h"
      INTEGER NINTCR,N,K,JJ,ID,NOPT,IER,NOPTCG
      INTEGER(8) IERI
      Double Precision ERI
      ALLOCATABLE::ERI(:),IERI(:)
      ALLOCATE(IERI(NINTCR),ERI(NINTCR),STAT=IER)
      IF(IER/=0)CALL ERRORMEM(NINTCR)
      ID=MY_ID
      K=MASTER
      JJ=1
!
   10 CALL MPI_RECV(N,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,STATUS,IERR)
      IF (N<0) GOTO 20
      CALL MPI_RECV(ERI(JJ),N,MPI_REAL8,K,ID,MPI_COMM_WORLD,STATUS,IERR)
      CALL MPI_RECV(IERI(JJ),N,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,STATUS, &
                    IERR)
 
      JJ=JJ+N
      GOTO 10
!
!     Now it has all the integrals in the node, waits for wake-up signal
!     and call the pertinent subroutine (either form J or K integrals)
!
   20 CALL MPI_RECV(NOPT,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,STATUS,IERR)
      IF (NOPT==0) THEN
       CALL MPI_RECV(NOPTCG,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,         &
                     STATUS,IERR)
       IF (NOPTCG==0) GOTO 30  
       RETURN 
      ELSE IF (NOPT==1) THEN
       CALL MPI_RECV(N,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,STATUS,IERR)
       CALL SLVHSTJ(N,JJ-1,IERI,ERI)
      ELSE IF (NOPT==2) THEN
       CALL MPI_RECV(N,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,STATUS,IERR)
       CALL SLVHSTK(N,JJ-1,IERI,ERI)
      ELSE IF (NOPT==3) THEN
       CALL MPI_RECV(N,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,STATUS,IERR)
       CALL SLVFORM2JK(N,JJ-1,IERI,ERI)
      ELSE IF (NOPT==4) THEN
       CALL MPI_RECV(N,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,STATUS,IERR)
       CALL SLVFORMJK(N,JJ-1,IERI,ERI)
      ELSE IF (NOPT==5) THEN
       CALL MPI_RECV(N,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,STATUS,IERR)
       CALL SLVHSTARK3(N,JJ-1,IERI,ERI)
      ELSE IF (NOPT==6) THEN
       CALL MPI_RECV(N,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,STATUS,IERR)
       CALL SLVFFJMN1rc(N)
      ELSE
       Call Abortx ("I do not know this MPI option!")
      ENDIF    
      GOTO 20
   30 CONTINUE

      End

! SlvIntRI
      Subroutine SlvIntRI(NINTCR)
      Implicit None
#include "mpip.h"
      INTEGER NINTCR,K,JJ,ID,NOPT,IER,NOPTCG
      INTEGER NBF,NBF5,NBFaux,IAUXDIM
      Double Precision ERIaux
      ALLOCATABLE::ERIaux(:)
      ALLOCATE(ERIaux(NINTCR),STAT=IER)
      IF(IER/=0)CALL ERRORMEM(NINTCR)
      ID=MY_ID
      K=MASTER
      JJ=1
      CALL MPI_BCAST(NBF,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NBF5,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NBFaux,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      IAUXDIM = NINTCR/(NBF*(NBF+1)/2)
      CALL MPI_RECV(ERIaux(JJ),NINTCR,MPI_REAL8,K,ID,MPI_COMM_WORLD,    &
                    STATUS,IERR)
!
!     Now it has all the integrals in the node, waits for wake-up signal
!     and call the pertinent subroutine (either form J or K integrals)
!
   20 CALL MPI_RECV(NOPT,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,STATUS,IERR)
      IF (NOPT==0) THEN
       CALL MPI_RECV(NOPTCG,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,         &
                     STATUS,IERR)
       IF (NOPTCG==0) GOTO 30
       RETURN 
      ELSE IF (NOPT==1) THEN
       CALL SLVHSTJKRI(NBF,IAUXDIM,ERIaux)
      ELSE IF (NOPT==2) THEN
       CALL SLVQJKMATmRI(NBF,NBF5,IAUXDIM,ERIaux)
      ELSE IF (NOPT==6) THEN
       CALL MPI_RECV(NBF,1,MPI_INTEGER,K,ID,MPI_COMM_WORLD,STATUS,IERR)
       CALL SLVFFJMN1rc(NBF)
      ELSE
       Call Abortx ("I do not know this MPI option!")
      ENDIF
      GOTO 20
   30 DEALLOCATE(ERIaux)
      CONTINUE
!
      End
      
! SLVHSTJKRI
      SUBROUTINE SLVHSTJKRI(NBF,IAUXDIM,ERIaux)
!     HSTARJKRI adapted for Slaves 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mpip.h"
      INTEGER NBF,IAUXDIM
      DIMENSION ERIaux(NBF*(NBF+1)/2,IAUXDIM)
      DOUBLE PRECISION,DIMENSION(NBF)::C
      DOUBLE PRECISION,DIMENSION(NBF*(NBF+1)/2)::FJ,FK,FFJ,FFK
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DOUBLE PRECISION,DIMENSION(NBF)::B_IN
      DOUBLE PRECISION::B_II
!
!     Get C matrix from Master node
!
      CALL MPI_BCAST(C,NBF,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)

      FJ = 0.0d0
      FK = 0.0d0
      FFJ = 0.0d0
      FFK = 0.0d0

      !$OMP PARALLEL DO PRIVATE(M, N, MN, B_II, B_IN) REDUCTION(+:FJ, FK)
      DO K=1,IAUXDIM
        B_IN(1:NBF) = 0.0d0
        B_II = 0.0d0

        DO N=1,NBF
          DO M=1,N
            MN = M + N*(N-1)/2
            B_IN(N) = B_IN(N) + C(M)*ERIaux(MN,K)
            IF(M.NE.N) B_IN(M) = B_IN(M)+C(N)*ERIaux(MN,K)
          END DO
        END DO

        DO N=1,NBF
          B_II = B_II + C(N)*B_IN(N)
        END DO

        DO N=1,NBF
          DO M=1,N
            MN = M + N*(N-1)/2
            FJ(MN) = FJ(MN) + B_II*ERIaux(MN,K)
          END DO
        END DO

        DO N=1,NBF
          DO M=1,N
            MN = M + N*(N-1)/2
            FK(MN) = FK(MN) + B_IN(M)*B_IN(N)
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Send pieces to master
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_REDUCE(FJ,FFJ,NBF*(NBF+1)/2,MPI_REAL8,MPI_SUM,MASTER,    &
                      MPI_COMM_WORLD,IERR)                               
      CALL MPI_REDUCE(FK,FFK,NBF*(NBF+1)/2,MPI_REAL8,MPI_SUM,MASTER,    &
                      MPI_COMM_WORLD,IERR)
!-----------------------------------------------------------------------
      RETURN
      END

! SLVQJKMATmRI
      SUBROUTINE SLVQJKMATmRI(NBF,NBF5,IAUXDIM,XIJKAUX)
!     QJKMATmRI adapted for Slaves 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mpip.h"
      INTEGER NBF,NBF5,IAUXDIM
      DIMENSION XIJKAUX(NBF*(NBF+1)/2,IAUXDIM)
      DIMENSION C(NBF,NBF)
      DOUBLE PRECISION,DIMENSION(NBF5*(NBF5+1)/2)::QJ,QK,QQJ,QQK
      DIMENSION B_IN(NBF5,NBF) 
      DIMENSION B_IJ(NBF5,NBF5) 
!
!     Get C matrix from Master node
!
      CALL MPI_BCAST(C,NBF*NBF,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)


      QJ = 0.0d0
      QK = 0.0d0
      QQJ = 0.0d0
      QQK = 0.0d0

      !$OMP PARALLEL DO PRIVATE(M, N, MN, I, J, IJ, B_IN, B_IJ) REDUCTION(+:QJ,QK) 
      DO K=1,IAUXDIM

        B_IN(1:NBF5,1:NBF) = 0.0d0
        DO I=1,NBF5
          DO N=1,NBF
            DO M=1,N
              MN = M + N*(N-1)/2
              B_IN(I,N) = B_IN(I,N) + C(M,I)*XIJKAUX(MN,K)
              IF(M.NE.N) B_IN(I,M) = B_IN(I,M)+C(N,I)*XIJKAUX(MN,K)
            END DO
          END DO
        END DO

        B_IJ(1:NBF5,1:NBF5) = 0.0d0
        DO J=1,NBF5
          DO N=1,NBF
            DO I=1,J
              B_IJ(I,J) = B_IJ(I,J) + C(N,J)*B_IN(I,N)
            END DO
          END DO
        END DO

        DO J=1,NBF5
          DO I=1,J
            IJ=I+J*(J-1)/2
            QJ(IJ) = QJ(IJ) + B_IJ(I,I)*B_IJ(J,J)
          END DO
        END DO

        DO J=1,NBF5
          DO I=1,J
            IJ=I+J*(J-1)/2
            QK(IJ) = QK(IJ) + B_IJ(I,J)*B_IJ(I,J)
          END DO
        END DO

      END DO
      !$OMP END PARALLEL DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Send pieces to master
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_REDUCE(QJ,QQJ,NBF5*(NBF5+1)/2,MPI_REAL8,MPI_SUM,MASTER,  &
                      MPI_COMM_WORLD,IERR)                               
      CALL MPI_REDUCE(QK,QQK,NBF5*(NBF5+1)/2,MPI_REAL8,MPI_SUM,MASTER,  &
                      MPI_COMM_WORLD,IERR)
!-----------------------------------------------------------------------
      RETURN
      END

! SLVHSTJ
      SUBROUTINE SLVHSTJ(N,NN,IERI,ERI)
!     HSTARJ adapted for Slaves 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/USEHUBBARD/IHUB
#include "mpip.h"
      INTEGER NBF
      DIMENSION ERI(NN)
      INTEGER(8) :: IERI(NN)
      INTEGER(8) :: LABEL
      ALLOCATABLE::P(:),F(:),FF(:)
      ALLOCATE (P(N),F(N),FF(N))
!
!     Get P matrix from Master node
!
      CALL MPI_BCAST(NBF,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(P,N,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
!
      F  = 0.0d0
      FF = 0.0d0
      !$OMP PARALLEL DO PRIVATE(LABEL, I, J, K, L, XJ, NIJ, NKL) REDUCTION(+:F)
      DO M=1,NN
       LABEL = IERI(M)
       CALL LABELIJKL(LABEL,I,J,K,L)
       XJ = ERI(M)
       NIJ = I*(I-1)/2 + J
       NKL = K*(K-1)/2 + L
       IF(IHUB==0) CALL OTTOINTEGR(I,J,K,L,NIJ,NKL,XJ)
                   F(NIJ)=F(NIJ)+0.5d0*P(NKL)*XJ
       IF(NIJ/=NKL)F(NKL)=F(NKL)+0.5d0*P(NIJ)*XJ
      ENDDO
      !$OMP END PARALLEL DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Send pieces to master
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_REDUCE(F,FF,N,MPI_REAL8,MPI_SUM,MASTER,MPI_COMM_WORLD,   &
                      IERR)
!-----------------------------------------------------------------------
      DEALLOCATE(P,F,FF)
      RETURN
      END

! SLVHSTK
      SUBROUTINE SLVHSTK(N,NN,IERI,ERI)
!     HSTARK adapted for Slaves
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/USEHUBBARD/IHUB      
#include "mpip.h"
      INTEGER NBF
      INTEGER(8),DIMENSION(NN)::IERI
      INTEGER(8) :: LABEL
      DOUBLE PRECISION,DIMENSION(NN)::ERI
      ALLOCATABLE::P(:),F(:),FF(:)
      ALLOCATE (P(N),F(N),FF(N))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Get P matrix from Master node
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_BCAST(NBF,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(P,N,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)

      F = 0.0d0
      FF = 0.0d0

      !$OMP PARALLEL DO PRIVATE(LABEL, I, J, K, L, XJ, XK, NIJ, NKL, NIK, NJL, NIL, NJK) REDUCTION(+:F)
      DO M=1,NN
      !DO 1 M=1,NN
       LABEL = IERI(M)
       CALL LABELIJKL(LABEL,I,J,K,L)
       XJ = ERI(M)
       NIJ = I*(I-1)/2 + J
       NKL = K*(K-1)/2 + L

       XJ = 0.25*XJ
       IF(IHUB==0)CALL OTTOINTEGR(I,J,K,L,NIJ,NKL,XJ)

       XK = XJ
       NIK = I*(I-1)/2 + K
       NJL = MAX0(J,L)*(MAX0(J,L)-1)/2 + MIN0(J,L)
       IF(I==K.OR.J==L)XK=XK+XK
                       F(NIK)=F(NIK)+P(NJL)*XK
       IF(NIK/=NJL)    F(NJL)=F(NJL)+P(NIK)*XK
       !IF(I==J.OR.K==L)GOTO 1
       IF(I==J.OR.K==L)CYCLE
       NIL = I*(I-1)/2 + L
       NJK = MAX0(J,K)*(MAX0(J,K)-1)/2 + MIN0(J,K)
       IF(I==L.OR.J==K)XJ=XJ+XJ
                       F(NIL)=F(NIL)+P(NJK)*XJ
       IF(NIL/=NJK)    F(NJK)=F(NJK)+P(NIL)*XJ
    !1 CONTINUE
      END DO
      !$OMP END PARALLEL DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Send pieces to master
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_REDUCE(F,FF,N,MPI_REAL8,MPI_SUM,MASTER,MPI_COMM_WORLD,   &
     &                IERR)
!-----------------------------------------------------------------------
      DEALLOCATE(FF,P,F)
      RETURN
      END

! SLVFORM2JK
      SUBROUTINE SLVFORM2JK(N,NN,IERI,ERI)
!     FORM2JK adapted for Slaves
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/USEHUBBARD/IHUB      
#include "mpip.h"
      INTEGER NBF
      INTEGER(8),DIMENSION(NN)::IERI
      INTEGER(8) :: LABEL
      DOUBLE PRECISION,DIMENSION(NN)::ERI
      ALLOCATABLE::P(:),F(:),FF(:)
      ALLOCATE (P(N),F(N),FF(N))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Get P matrix from Master node
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_BCAST(NBF,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(P,N,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)

      F = 0.0d0
      FF = 0.0d0

      !$OMP PARALLEL DO PRIVATE(LABEL, XJ, XK, I, J, K, L, NIJ, NKL, NIK, NJL, NIL, NJK) REDUCTION(+:F)
      DO M=1,NN
       LABEL = IERI(M)
       CALL LABELIJKL(LABEL,I,J,K,L)
!-----------------------------------------------------------------------
!      2*J
!-----------------------------------------------------------------------
       XJ = ERI(M)
       NIJ = I*(I-1)/2 + J
       NKL = K*(K-1)/2 + L
       IF(IHUB==0) CALL OTTOINTEGR(I,J,K,L,NIJ,NKL,XJ)
                   F(NIJ)=F(NIJ)+P(NKL)*XJ
       IF(NIJ/=NKL)F(NKL)=F(NKL)+P(NIJ)*XJ
!-----------------------------------------------------------------------
!      -K
!-----------------------------------------------------------------------
       XJ = 0.25d0*XJ
       XK = XJ
       NIK = I*(I-1)/2 + K
       NJL = MAX0(J,L)*(MAX0(J,L)-1)/2 + MIN0(J,L)
       IF(I==K.OR.J==L) XK=XK+XK
                          F(NIK)=F(NIK)-P(NJL)*XK
       IF(NIK/=NJL)       F(NJL)=F(NJL)-P(NIK)*XK
       IF(I/=J.and.K/=L)THEN
        NIL = I*(I-1)/2 + L
        NJK = MAX0(J,K)*(MAX0(J,K)-1)/2 + MIN0(J,K)
        IF(I==L.OR.J==K) XJ=XJ+XJ
                           F(NIL)=F(NIL)-P(NJK)*XJ
        IF(NIL/=NJK)       F(NJK)=F(NJK)-P(NIL)*XJ
       ENDIF
      ENDDO
      !$OMP END PARALLEL DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Send pieces to master
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_REDUCE(F,FF,N,MPI_REAL8,MPI_SUM,MASTER,MPI_COMM_WORLD,   &
                      IERR)
!-----------------------------------------------------------------------
      DEALLOCATE(P,F,FF)
      RETURN
      END

! SLVFORMJK
      SUBROUTINE SLVFORMJK(N,NN,IERI,ERI)
!     FORMJK adapted for Slaves
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/USEHUBBARD/IHUB      
#include "mpip.h"
      INTEGER NBF
      INTEGER(8),DIMENSION(NN)::IERI
      INTEGER(8) :: LABEL
      DOUBLE PRECISION,DIMENSION(NN)::ERI
      ALLOCATABLE::P(:),F(:),FF(:)
      ALLOCATE (P(N),F(N),FF(N))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Get P matrix from Master node
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_BCAST(NBF,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(P,N,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)

      F = 0.0d0
      FF = 0.0d0

      DO M=1,NN
       LABEL = IERI(M)
       CALL LABELIJKL(LABEL,I,J,K,L)
!-----------------------------------------------------------------------
!      J
!-----------------------------------------------------------------------
       XJ = ERI(M)
       NIJ = I*(I-1)/2 + J
       NKL = K*(K-1)/2 + L
       IF(IHUB==0) CALL OTTOINTEGR(I,J,K,L,NIJ,NKL,XJ)
                   F(NIJ)=F(NIJ)+0.5d0*P(NKL)*XJ
       IF(NIJ/=NKL)F(NKL)=F(NKL)+0.5d0*P(NIJ)*XJ
!-----------------------------------------------------------------------
!      -K
!-----------------------------------------------------------------------
       XJ = 0.25d0*XJ
       XK = XJ
       NIK = I*(I-1)/2 + K
       NJL = MAX0(J,L)*(MAX0(J,L)-1)/2 + MIN0(J,L)
       IF(I==K.OR.J==L) XK=XK+XK
                          F(NIK)=F(NIK)-P(NJL)*XK
       IF(NIK/=NJL)       F(NJL)=F(NJL)-P(NIK)*XK
       IF(I/=J.and.K/=L)THEN
        NIL = I*(I-1)/2 + L
        NJK = MAX0(J,K)*(MAX0(J,K)-1)/2 + MIN0(J,K)
        IF(I==L.OR.J==K) XJ=XJ+XJ
                           F(NIL)=F(NIL)-P(NJK)*XJ
        IF(NIL/=NJK)       F(NJK)=F(NJK)-P(NIL)*XJ
       ENDIF
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Send pieces to master
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_REDUCE(F,FF,N,MPI_REAL8,MPI_SUM,MASTER,MPI_COMM_WORLD,   &
                      IERR)
!-----------------------------------------------------------------------
      DEALLOCATE(P,F,FF)
      RETURN
      END

! SLVHSTARK3
      SUBROUTINE SLVHSTARK3(N,NN,IERI,ERI)
!     HSTARK3 adapted for Slaves
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/USEHUBBARD/IHUB      
#include "mpip.h"
      INTEGER NBF
      INTEGER(8),DIMENSION(NN)::IERI
      INTEGER(8) :: LABEL
      DOUBLE PRECISION,DIMENSION(NN)::ERI
      ALLOCATABLE :: P1(:),F1(:),FF1(:)
      ALLOCATABLE :: P2(:),F2(:),FF2(:)
      ALLOCATABLE :: P3(:),F3(:),FF3(:)      
      ALLOCATE(P1(N),F1(N),FF1(N),P2(N),F2(N),FF2(N),P3(N),F3(N),FF3(N))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Get P matrix from Master node
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_BCAST(NBF,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(P1,N,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(P2,N,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(P3,N,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)      

      F1 = 0.0d0
      FF1 = 0.0d0
      F2 = 0.0d0
      FF2 = 0.0d0
      F3 = 0.0d0
      FF3 = 0.0d0

      !$OMP PARALLEL DO PRIVATE(LABEL, I, J, K, L, XJ, XK, NIJ, NKL, NIK, NJL, NIL, NJK)  REDUCTION(+:F1,F2,F3)
      DO M=1,NN
      !DO 1 M=1,NN
       LABEL = IERI(M)
       CALL LABELIJKL(LABEL,I,J,K,L)
       XJ = ERI(M)
       NIJ = I*(I-1)/2 + J
       NKL = K*(K-1)/2 + L

       XJ = 0.25*XJ
       IF(IHUB==0)CALL OTTOINTEGR(I,J,K,L,NIJ,NKL,XJ)

       XK = XJ
       NIK = I*(I-1)/2 + K
       NJL = MAX0(J,L)*(MAX0(J,L)-1)/2 + MIN0(J,L)
       IF(I==K.OR.J==L)XK=XK+XK
       
       F1(NIK)=F1(NIK)+P1(NJL)*XK
       F2(NIK)=F2(NIK)+P2(NJL)*XK
       F3(NIK)=F3(NIK)+P3(NJL)*XK
                       
       IF(NIK/=NJL)THEN
        F1(NJL)=F1(NJL)+P1(NIK)*XK
        F2(NJL)=F2(NJL)+P2(NIK)*XK
        F3(NJL)=F3(NJL)+P3(NIK)*XK        
       ENDIF
       
       !IF(I==J.OR.K==L)GOTO 1
       IF(I==J.OR.K==L)CYCLE
       NIL = I*(I-1)/2 + L
       NJK = MAX0(J,K)*(MAX0(J,K)-1)/2 + MIN0(J,K)
       IF(I==L.OR.J==K)XJ=XJ+XJ
       
       F1(NIL)=F1(NIL)+P1(NJK)*XJ
       F2(NIL)=F2(NIL)+P2(NJK)*XJ       
       F3(NIL)=F3(NIL)+P3(NJK)*XJ       
       
       IF(NIL/=NJK)THEN
        F1(NJK)=F1(NJK)+P1(NIL)*XJ
        F2(NJK)=F2(NJK)+P2(NIL)*XJ
        F3(NJK)=F3(NJK)+P3(NIL)*XJ        
       ENDIF
       
    !1 CONTINUE
      END DO
      !$OMP END PARALLEL DO

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Send pieces to master
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_REDUCE(F1,FF1,N,MPI_REAL8,MPI_SUM,MASTER,MPI_COMM_WORLD, &
     &                IERR)
      CALL MPI_REDUCE(F2,FF2,N,MPI_REAL8,MPI_SUM,MASTER,MPI_COMM_WORLD, &
     &                IERR)
      CALL MPI_REDUCE(F3,FF3,N,MPI_REAL8,MPI_SUM,MASTER,MPI_COMM_WORLD, &
     &                IERR)
!-----------------------------------------------------------------------
      DEALLOCATE(FF1,P1,F1,FF2,P2,F2,FF3,P3,F3)
      RETURN
      END

! SLVFFJMN1rc
      SUBROUTINE SLVFFJMN1rc(NBF)
!     FFJMN1rc adapted for Slaves 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mpip.h"
      INTEGER NBF,NBF5,NSQ
      INTEGER NA,NB,NO1,NSOC
      INTEGER LL,UL,EQPART,UNEQPART
      INTEGER i,k,ik
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::DJ,DK,F,FF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::CJ12,CK12
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::RO
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::H
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NSQ = NBF*NBF
      CALL MPI_BCAST(NBF5,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)

      ALLOCATE(DJ(NSQ,NBF5), DK(NSQ,NBF5), F(NSQ,NBF5), FF(NSQ,NBF5))
      ALLOCATE(CJ12(NBF5,NBF5), CK12(NBF5,NBF5))
      ALLOCATE(H(NBF,NBF))
      ALLOCATE(RO(NBF5))
!
      CALL MPI_BCAST(DJ,NSQ*NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DK,NSQ*NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(CJ12,NBF5*NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(CK12,NBF5*NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(RO,NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(H,NBF*NBF,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NA,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NB,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NO1,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NSOC,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
!
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,MY_ID,IERR)
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,NPROCS,IERR)

      F = 0.0D0
      FF = 0.0D0
      EQPART = INT(NBF*NBF/NPROCS)
      UNEQPART = NBF*NBF - NPROCS*EQPART
      LL = MY_ID * EQPART + 1 + MIN(MY_ID,UNEQPART)
      UL = LL + EQPART - 1
      IF(MY_ID<UNEQPART) UL = UL + 1
!-----------------------------------------------------------------------
!     Calculate Fj(m,n)
!-----------------------------------------------------------------------
      IF(NO1>1)THEN
       !$OMP PARALLEL DO PRIVATE(i,k,ik,J)
       do ik=LL,UL
        i = INT((ik-1)/nbf) + 1
        k = MOD(ik-1,nbf) + 1
        F(ik,1) = H(i,k) + SLVPRODWCWij(ik,DJ,CJ12,NBF5,NSQ)           &
                         - SLVPRODWCWij(ik,DK,CK12,NBF5,NSQ)
        DO J=NO1+1,NB
         F(ik,J) = RO(J) * ( H(i,k) + DJ(ik,J) )                       &
                 + SLVPRODWCWijq(ik,J,DJ,CJ12,NBF5,NSQ)                &
                 - SLVPRODWCWijq(ik,J,DK,CK12,NBF5,NSQ)  
        ENDDO                                                           
        if(NSOC>0)then                                                  
         DO J=NB+1,NA                                                   
          F(ik,J) = RO(J) * H(i,k)                                     &
                  + SLVPRODWCWijq(ik,J,DJ,CJ12,NBF5,NSQ)               & 
                  - SLVPRODWCWijq(ik,J,DK,CK12,NBF5,NSQ) 
         ENDDO                                                          
        end if                                                          
        DO J=NA+1,NBF5                                                  
         F(ik,J) = RO(J) * ( H(i,k) + DJ(ik,J) )                       &
                 + SLVPRODWCWijq(ik,J,DJ,CJ12,NBF5,NSQ)                & 
                 - SLVPRODWCWijq(ik,J,DK,CK12,NBF5,NSQ)
        ENDDO
       enddo
       !$OMP END PARALLEL DO
       DO J=2,NO1
        F(1:NSQ,J) = F(1:NSQ,1)
       ENDDO
      ELSE
       !$OMP PARALLEL DO PRIVATE(i,k,ik,J)
       do ik=LL,UL
        i = INT((ik-1)/NBF) + 1
        k = MOD(ik-1,NBF) + 1
        DO J=1,NB
         F(ik,J) = RO(J) * ( H(i,k) + DJ(ik,J) )                       &
                 + SLVPRODWCWijq(ik,J,DJ,CJ12,NBF5,NSQ)                &
                 - SLVPRODWCWijq(ik,J,DK,CK12,NBF5,NSQ)  
        ENDDO                                                           
        if(NSOC>0)then                                                  
         DO J=NB+1,NA                                                   
          F(ik,J) = RO(J) * H(i,k)                                     &
                  + SLVPRODWCWijq(ik,J,DJ,CJ12,NBF5,NSQ)               &
                  - SLVPRODWCWijq(ik,J,DK,CK12,NBF5,NSQ)  
         ENDDO                                                          
        end if                                                          
        DO J=NA+1,NBF5                                                  
         F(ik,J) = RO(J) * ( H(i,k) + DJ(ik,J) )                       &
                 + SLVPRODWCWijq(ik,J,DJ,CJ12,NBF5,NSQ)                &
                 - SLVPRODWCWijq(ik,J,DK,CK12,NBF5,NSQ)  
        ENDDO
       enddo
       !$OMP END PARALLEL DO
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Send pieces to master
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_REDUCE(F,FF,NSQ*NBF5,MPI_REAL8,MPI_SUM,MASTER,           &
                      MPI_COMM_WORLD,IERR)
!-----------------------------------------------------------------------
      DEALLOCATE(DJ, DK, F, FF)
      DEALLOCATE(CJ12, CK12)
      DEALLOCATE(H)
      DEALLOCATE(RO)
!-----------------------------------------------------------------------
      RETURN
      END

! SLVPRODWCWijq
      FUNCTION SLVPRODWCWijq(ij,IQ,W,CW12,NBF5,NSQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CW12
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::W
!-----------------------------------------------------------------------
!     Sum W*CW12(IQ,*) by IQP, IQ is not considered
!-----------------------------------------------------------------------
      SLVPRODWCWijq = 0.0d0
      DO IQP=1,IQ-1
       SLVPRODWCWijq = SLVPRODWCWijq + W(ij,IQP)*CW12(IQ,IQP)
      ENDDO
      DO IQP=IQ+1,NBF5
       SLVPRODWCWijq = SLVPRODWCWijq + W(ij,IQP)*CW12(IQ,IQP)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! SLVPRODWCWij
      FUNCTION SLVPRODWCWij(ij,W,CW12,NBF5,NSQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CW12
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::W
!-----------------------------------------------------------------------
!     Sum W*CW12(1,*) by IQP
!-----------------------------------------------------------------------
      SLVPRODWCWij = 0.0d0
      DO IQP=1,NBF5
       SLVPRODWCWij = SLVPRODWCWij + W(ij,IQP)*CW12(1,IQP)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

#endif

! ERRORMEM
      SUBROUTINE ERRORMEM(NMEMORY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER :: NMEMORY
      DOUBLE PRECISION :: GBNUMBER
      GBNUMBER = DFLOAT(NMEMORY*8)/(1024.0d0**3)
      WRITE(6,1)NMEMORY,GBNUMBER
    1 FORMAT(//10X,'Sorry, You need more Memory!, NMEMORY =',I20,2X,    &
                   '=',F10.2,' GB')
      STOP
      END

!----------------------------------------------------------------------!
!                                                                      !
!       2025 Use LIBCINT open source library for ERI calculation       !
!                                                                      !
!  implemented by Juan Felipe Huan Lew Yee and Jorge Martin del Campo  !
!                                                                      !
!----------------------------------------------------------------------!

! ECPINTlib
      SUBROUTINE ECPINTl(H,NBFT,SIZE_ENV,ENV,NAT,ATM,NSHELL,   &
                      NBAS,BAS,IGTYP)
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: H(NBFT)
      INTEGER, INTENT(IN) :: NBFT, SIZE_ENV, NAT, NSHELL, NBAS, IGTYP
      DOUBLE PRECISION, INTENT(IN) :: ENV(SIZE_ENV)
      INTEGER, INTENT(IN) :: ATM(6, NAT), BAS(8, NSHELL)

      INTEGER :: I, J
      INTEGER :: ISH, JSH
      INTEGER :: ERR
      INTEGER :: DI, DJ, DJEFF
      INTEGER :: LI, LJ, ID, IJ
      INTEGER(4), DIMENSION(2) :: SHLS
      INTEGER, DIMENSION(NSHELL) :: LOC, Dcgto

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: BLK

      INTEGER, EXTERNAL :: CINTcgto_cart, CINTcgto_spheric
      INTEGER, EXTERNAL :: CECPscalar_cart, CECPscalar_sph
!-----------------------------------------------------------------------
! LOC indicates the starting location of each shell in the AO basis
!-----------------------------------------------------------------------
      LOC(1) = 1
      IF(IGTYP==1) DI = CINTcgto_cart(0, BAS)
      IF(IGTYP==2) DI = CINTcgto_spheric(0, BAS)
      Dcgto(1) = DI
      DO ISH = 2,NSHELL
        LOC(ISH) = LOC(ISH-1) + DI
        IF(IGTYP==1) DI = CINTcgto_cart(ISH-1, BAS)
        IF(IGTYP==2) DI = CINTcgto_spheric(ISH-1, BAS)
        Dcgto(ISH) = DI
      END DO
!-----------------------------------------------------------------------
!                      H, S & TKIN integrals
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     I SHELL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO ISH = 1,NSHELL
        SHLS(1) = ISH - 1
        DI = Dcgto(ISH)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      J SHELL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        DO JSH = 1,ISH
          SHLS(2) = JSH - 1
          DJ = Dcgto(JSH)

          ALLOCATE(BLK(DI, DJ))
          BLK = 0.0D0

          IF(IGTYP==1) THEN
            ERR = CECPscalar_cart(BLK, SHLS, ATM, NAT, BAS, NBAS, &
                                   ENV)
          ELSE IF(IGTYP==2) THEN
            ERR = CECPscalar_sph(BLK, SHLS, ATM, NAT, BAS, NBAS, &
                                   ENV)
          END IF

          DO I = 1, DI
            LI = LOC(ISH) + I - 1
            ID = (LI * (LI - 1)) / 2
            DJEFF = MERGE(I, DJ, ISH == JSH)
            DO J = 1, DJEFF
              LJ = LOC(JSH) + J - 1
              IJ = LJ + ID
              H(IJ) =  H(IJ) + BLK(I, J)
            END DO
          END DO

          DEALLOCATE(BLK)
        END DO
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

! HSandTl
      SUBROUTINE HSandTl(H,S,TKIN,NBFT,SIZE_ENV,ENV,NAT,ATM,NSHELL,     &
                         NBAS,BAS,IGTYP)
      IMPLICIT NONE
!
      DOUBLE PRECISION, INTENT(OUT) :: H(NBFT), S(NBFT), TKIN(NBFT)
      INTEGER, INTENT(IN) :: NBFT, SIZE_ENV, NAT, NSHELL, NBAS, IGTYP
      DOUBLE PRECISION, INTENT(IN) :: ENV(SIZE_ENV)
      INTEGER, INTENT(IN) :: ATM(6, NAT), BAS(8, NSHELL)
!
      INTEGER :: I, J
      INTEGER :: ISH, JSH
      INTEGER :: ERR
      INTEGER :: DI, DJ, DJEFF
      INTEGER :: LI, LJ, ID, IJ
      INTEGER(4), DIMENSION(2) :: SHLS
      INTEGER, DIMENSION(NSHELL) :: LOC, Dcgto
!
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: SBLK, VBLK, TBLK
!
      INTEGER, EXTERNAL :: CINTcgto_spheric, CINT1e_ovlp_sph
      INTEGER, EXTERNAL :: CINT1e_kin_sph, CINT1e_nuc_sph
      INTEGER, EXTERNAL :: CINTcgto_cart, CINT1e_ovlp_cart
      INTEGER, EXTERNAL :: CINT1e_kin_cart, CINT1e_nuc_cart
!-----------------------------------------------------------------------
!     LOC indicates the starting location of each shell in the AO basis
!-----------------------------------------------------------------------
      LOC(1) = 1
      IF(IGTYP==1) DI = CINTcgto_cart(0, BAS)
      IF(IGTYP==2) DI = CINTcgto_spheric(0, BAS)
      Dcgto(1) = DI
      DO ISH = 2,NSHELL
        LOC(ISH) = LOC(ISH-1) + DI
        IF(IGTYP==1) DI = CINTcgto_cart(ISH-1, BAS)
        IF(IGTYP==2) DI = CINTcgto_spheric(ISH-1, BAS)
        Dcgto(ISH) = DI
      END DO
!-----------------------------------------------------------------------
!                      H, S & TKIN integrals
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     I SHELL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO ISH = 1,NSHELL
        SHLS(1) = ISH - 1
        DI = Dcgto(ISH)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      J SHELL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        DO JSH = 1,ISH
          SHLS(2) = JSH - 1
          DJ = Dcgto(JSH)
!
          ALLOCATE(SBLK(DI, DJ), TBLK(DI, DJ), VBLK(DI, DJ))
!
          IF(IGTYP==1) THEN
            ERR = CINT1e_ovlp_cart(SBLK,SHLS,ATM,NAT,BAS,NBAS,ENV,0_8)
            ERR = CINT1e_kin_cart(TBLK,SHLS,ATM,NAT,BAS,NBAS,ENV,0_8)
            ERR = CINT1e_nuc_cart(VBLK,SHLS,ATM,NAT,BAS,NBAS,ENV,0_8)
          ELSE IF(IGTYP==2) THEN
            ERR = CINT1e_ovlp_sph(SBLK,SHLS,ATM,NAT,BAS,NBAS,ENV,0_8)
            ERR = CINT1e_kin_sph(TBLK,SHLS,ATM,NAT,BAS,NBAS,ENV,0_8)
            ERR = CINT1e_nuc_sph(VBLK,SHLS,ATM,NAT,BAS,NBAS,ENV,0_8)
          END IF
!
          DO I = 1, DI
            LI = LOC(ISH) + I - 1
            ID = (LI * (LI - 1)) / 2
            DJEFF = MERGE(I, DJ, ISH == JSH)
            DO J = 1, DJEFF
              LJ = LOC(JSH) + J - 1
              IJ = LJ + ID
              H(IJ) =  TBLK(I, J) + VBLK(I, J)
              S(IJ) =  SBLK(I, J)
              TKIN(IJ) =  TBLK(I, J)
            END DO
          END DO
!
          DEALLOCATE(SBLK,TBLK,VBLK)
        END DO
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

! PRCALCl
      SUBROUTINE PRCALCl(XVAL,NVAL,L2,SIZE_ENV,ENV,NAT,ATM,NSHELL,      &
                         NBAS,BAS,IGTYP)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NVAL, L2, SIZE_ENV, NAT, NSHELL, NBAS
      INTEGER, INTENT(IN) :: IGTYP
      DOUBLE PRECISION, INTENT(INOUT) :: ENV(SIZE_ENV)
      INTEGER, INTENT(IN) :: ATM(6, NAT), BAS(8, NBAS)
      DOUBLE PRECISION, INTENT(OUT) :: XVAL(NVAL * L2)

      INTEGER :: IAT, ISH, JSH, I, J, K, LI, LJ, ID, IJ, NL2
      INTEGER :: DI, DJ, DJEFF, Ztot, Z, ERR
      INTEGER(4) :: SHLS(2)
      INTEGER, DIMENSION(NSHELL) :: LOC, Dcgto
      DOUBLE PRECISION :: CX, CY, CZ, RC(3)

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: WINT

      INTEGER, EXTERNAL :: CINTcgto_spheric, CINT1e_r_sph
      INTEGER, EXTERNAL :: CINTcgto_cart, CINT1e_r_cart
      COMMON/CMCoord/Xcm,Ycm,Zcm
      DOUBLE PRECISION::Xcm, Ycm, Zcm
!-----------------------------------------------------------------------
! LOC indicates the starting location of each shell in the AO basis
!-----------------------------------------------------------------------
      LOC(1) = 1
      IF(IGTYP==1) DI = CINTcgto_cart(0, BAS)
      IF(IGTYP==2) DI = CINTcgto_spheric(0, BAS)
      Dcgto(1) = DI
      DO ISH = 2,NSHELL
        LOC(ISH) = LOC(ISH-1) + DI
        IF(IGTYP==1) DI = CINTcgto_cart(ISH-1, BAS)
        IF(IGTYP==2) DI = CINTcgto_spheric(ISH-1, BAS)
        Dcgto(ISH) = DI
      END DO
!-----------------------------------------------------------------------
!                      CENTER OF MASS
!-----------------------------------------------------------------------
      !Ztot = 0
      !RC = 0.0D0
      !DO IAT=1, NAT
      !  Z = ATM(1,IAT)
      !  CX = ENV(ATM(2,IAT)+1)
      !  CY = ENV(ATM(2,IAT)+2)
      !  CZ = ENV(ATM(2,IAT)+3)
      !
      !  Ztot = Ztot + Z
      !  RC(1) = RC(1) + Z * CX
      !  RC(2) = RC(2) + Z * CY
      !  RC(3) = RC(3) + Z * CZ
      !END DO
      !RC = RC / Ztot
      ENV(2) = Xcm
      ENV(3) = Ycm
      ENV(4) = Zcm

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO ISH=1, NSHELL
        SHLS(1) = ISH - 1
        DI = Dcgto(ISH)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        DO JSH=1, ISH
          SHLS(2) = JSH - 1
          DJ = Dcgto(JSH)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!          Integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ALLOCATE(WINT(DI,DJ,3))
          IF(IGTYP==1) ERR = CINT1e_r_cart(WINT, SHLS, ATM, NAT, BAS, &
                  NBAS, ENV)
          IF(IGTYP==2) ERR = CINT1e_r_sph(WINT, SHLS, ATM, NAT, BAS,  &
                  NBAS, ENV)

          DO K=1, NVAL
            NL2 = (K-1)*L2

            DO I = 1, DI
              LI = LOC(ISH) + I - 1
              ID = (LI * (LI - 1))/2  + NL2
              DJEFF = MERGE(I, DJ, ISH == JSH)
              DO J = 1, DJEFF
                LJ = LOC(JSH) + J - 1
                IJ = LJ + ID
                XVAL(IJ) = WINT(I, J, K)
              END DO
            END DO
          END DO
          DEALLOCATE(WINT)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END DO
!-----------------------------------------------------------------------
      ENV(2:4) = 0.0D0
      RETURN
      END

! ExchangeIntl
      SUBROUTINE ExchangeIntl(XINTS,NSH2,SIZE_ENV,ENV,NAT,ATM,NSHELL,   &
                              NBAS,BAS,IGTYP)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NSH2, SIZE_ENV, NAT, NSHELL, NBAS, IGTYP
      DOUBLE PRECISION, INTENT(OUT) :: XINTS(NSH2)
      DOUBLE PRECISION, INTENT(IN) :: ENV(SIZE_ENV)
      INTEGER, INTENT(IN) :: ATM(6, NAT), BAS(8, NBAS)
      INTEGER :: ISH, JSH, ID, IJ, DI, DJ, DJEFF, I, J, ERR
      INTEGER(4) :: SHLS(4)
      DOUBLE PRECISION :: VMAX, VAL
      INTEGER, DIMENSION(NSHELL) :: LOC, Dcgto
      DOUBLE PRECISION, ALLOCATABLE :: BLK(:,:,:,:)
      INTEGER, EXTERNAL :: CINTcgto_spheric, CINT2e_sph
      INTEGER, EXTERNAL :: CINTcgto_cart, CINT2e_cart
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LOC(1) = 1
      IF(IGTYP==1) DI = CINTcgto_cart(0, BAS)
      IF(IGTYP==2) DI = CINTcgto_spheric(0, BAS)
      Dcgto(1) = DI
      DO ISH = 2,NSHELL
        LOC(ISH) = LOC(ISH-1) + DI
        IF(IGTYP==1) DI = CINTcgto_cart(ISH-1, BAS)
        IF(IGTYP==2) DI = CINTcgto_spheric(ISH-1, BAS)
        Dcgto(ISH) = DI
      END DO
!-----------------------------------------------------------------------
      DO ISH = 1,NSHELL
        ID = (ISH * (ISH - 1)) / 2
        DI = Dcgto(ISH)
        SHLS(1) = ISH - 1
        SHLS(3) = ISH - 1

        DO JSH = 1,ISH
          IJ = ID + JSH
          DJ = Dcgto(JSH)
          SHLS(2) = JSH - 1
          SHLS(4) = JSH - 1

          ALLOCATE(BLK(DI, DJ, DI, DJ))
          IF(IGTYP==1) ERR = CINT2e_cart(BLK, SHLS, ATM, NAT, BAS,    &
                  NBAS, ENV, 0_8)
          IF(IGTYP==2) ERR = CINT2e_sph(BLK, SHLS, ATM, NAT, BAS,     &
                  NBAS, ENV, 0_8)

          VMAX = 0.0D0
          DO I=1,DI
            DJEFF = MERGE(I, DJ, ISH == JSH)
            DO J=1,DJEFF
              VAL = BLK(I,J,I,J)
              IF(VAL > VMAX) VMAX = VAL
            END DO
          END DO
          XINTS(IJ) = SQRT(VMAX)
          DEALLOCATE(BLK)
        END DO
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

! TwoERIl
      SUBROUTINE TwoERIl(SCHWRZ,NINTEGtm,NINTEGt,NSCHWZ,BUFP,IX,BUFP2,  &
                         IX2,XINTS,NSH2,IDONTW,IPRINTOPT,NSHELL,Cxyz,   &
                         NAT,NINTMX,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,NREC)
      IMPLICIT NONE
      DOUBLE PRECISION :: CUTOFF
      INTEGER :: ISCHWZ,IECP,NECP
      COMMON/INTOPT/CUTOFF,ISCHWZ,IECP,NECP
!
      INTEGER, INTENT(IN) :: NINTEGtm, NSHELL, NAT, NINTMX, IGTYP
      INTEGER, INTENT(IN) :: NSH2, SIZE_ENV, NBAS
      INTEGER, INTENT(IN) :: IDONTW, IPRINTOPT
      LOGICAL, INTENT(IN) :: SCHWRZ
!
      DOUBLE PRECISION, INTENT(IN) :: Cxyz(3,NAT)
      DOUBLE PRECISION, INTENT(IN) :: XINTS(NSH2)
      DOUBLE PRECISION, INTENT(INOUT) :: BUFP(NINTMX), BUFP2(NINTEGtm)
      INTEGER(8), INTENT(INOUT) :: IX(NINTMX), IX2(NINTEGtm)
      INTEGER, INTENT(INOUT) :: NINTEGt, NSCHWZ
      DOUBLE PRECISION, INTENT(IN) :: ENV(SIZE_ENV)
      INTEGER, INTENT(IN) :: ATM(6,NAT), BAS(8,NBAS)
!
      INTEGER :: ISH, JSH, KSH, LSH, II, JJ, KK, LL, IEXCH
      INTEGER :: DI, DJ, DK, DL, IJIJ, KLKL, NREC, ICOUNT
      INTEGER :: ERR, LOC(NSHELL),Dcgto(NSHELL)
      DOUBLE PRECISION :: TEST
      LOGICAL :: SKIPA, SKIPB, SKIPC, SCHSKP
      DOUBLE PRECISION, ALLOCATABLE :: BUFPC(:)
      INTEGER(8), ALLOCATABLE :: IXC(:)
      DOUBLE PRECISION, ALLOCATABLE :: BLK(:,:,:,:)
      INTEGER(4) :: SHLS(4)
      INTEGER, EXTERNAL :: CINTcgto_spheric, CINT2e_sph
      INTEGER, EXTERNAL :: CINTcgto_cart, CINT2e_cart
      DOUBLE PRECISION,PARAMETER :: TENM12 = 1.0D-12
!-----------------------------------------------------------------------
! LOC indicates the starting location of each shell in the AO basis
!-----------------------------------------------------------------------
      LOC(1) = 1
      IF(IGTYP==1) DI = CINTcgto_cart(0, BAS)
      IF(IGTYP==2) DI = CINTcgto_spheric(0, BAS)
      Dcgto(1) = DI
      DO ISH = 2,NSHELL
        LOC(ISH) = LOC(ISH-1) + DI
        IF(IGTYP==1) DI = CINTcgto_cart(ISH-1, BAS)
        IF(IGTYP==2) DI = CINTcgto_spheric(ISH-1, BAS)
        Dcgto(ISH) = DI
      END DO
!-----------------------------------------------------------------------
!     2e- Integrals (S,P,D,F,G & L Shells)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(BUFPC(NINTMX))
      ALLOCATE(IXC(NINTMX))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NINTEGt  = 0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !$OMP PARALLEL &
      !$OMP PRIVATE(II, JJ, KK, LL, ISH, JSH, KSH, LSH, DI, DJ, DK, DL, &
      !$OMP SKIPA, SKIPB, SKIPC, IEXCH, IJIJ, KLKL, SHLS, ERR, TEST,    &
      !$OMP SCHSKP, NSCHWZ, IX, BUFP, BLK, ICOUNT)
      SCHSKP=.FALSE.
      ICOUNT = 1
      !$OMP DO SCHEDULE(DYNAMIC)
      DO II = 1, NSHELL                                ! II Shell
        DO JJ = 1, II                                  ! JJ Shell
          DO KK = 1, JJ                                ! KK Shell
            DO LL = 1, KK                              ! LL Shell

              SKIPA =  (JJ == KK)
              SKIPB = (II == KK) .or. (JJ == LL)
              SKIPC = (II == JJ) .or. (KK == LL)

              DO IEXCH = 1, 3
!- - - - - - - - - - - - - - - - - - - - - - - -
!               (II,JJ//KK,LL)
!- - - - - - - - - - - - - - - - - - - - - - - -
                IF (IEXCH == 1) THEN
                  ISH = II
                  JSH = JJ
                  KSH = KK
                  LSH = LL
!- - - - - - - - - - - - - - - - - - - - - - - -
!               (II,KK//JJ,LL)
!- - - - - - - - - - - - - - - - - - - - - - - -
                ELSE IF(IEXCH == 2) THEN
                  IF (SKIPA) CYCLE
                  ISH = II
                  JSH = KK
                  KSH = JJ
                  LSH = LL
!- - - - - - - - - - - - - - - - - - - - - - - -
!               (II,LL//JJ,KK)
!- - - - - - - - - - - - - - - - - - - - - - - -
                ELSE IF(IEXCH == 3) THEN
                  IF(SKIPB .or. SKIPC) CYCLE
                  ISH = II
                  JSH = LL
                  KSH = JJ
                  LSH = KK
                END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!               Compute 2e- Integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                IF(SCHWRZ)THEN                    ! Schwarz Inequality
                  IJIJ = (ISH*ISH-ISH)/2 + JSH
                  KLKL = (KSH*KSH-KSH)/2 + LSH
                  TEST = XINTS(IJIJ) * XINTS(KLKL)
                  SCHSKP = TEST < CUTOFF
                  IF(SCHSKP) THEN
                    NSCHWZ = NSCHWZ + 1
                    CYCLE
                  END IF
                END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!               Select integral code for ERI calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                DI = Dcgto(ISH)
                DJ = Dcgto(JSH)
                DK = Dcgto(KSH)
                DL = Dcgto(LSH)
                ALLOCATE(BLK(DI,DJ,DK,DL))

                SHLS = (/ ISH-1, JSH-1, KSH-1, LSH-1 /)
                IF(IGTYP==1) ERR = CINT2e_cart(BLK, SHLS, ATM, NAT,    &
                        BAS, NSHELL, ENV, 0_8)
                IF(IGTYP==2) ERR = CINT2e_sph(BLK, SHLS, ATM, NAT,     &
                        BAS, NSHELL, ENV, 0_8)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!               Write Label & Integral on File 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                CALL QOUTl(BUFP,IX,BUFP2,IX2,NINTEGtm,NINTMX,IDONTW,    &
                           NSHELL,BLK,DI,DJ,DK,DL,ISH,JSH,KSH,LSH,      &
                           LOC,NREC,ICOUNT,CUTOFF)
                DEALLOCATE(BLK)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              END DO
            END DO
          END DO
        END DO
      END DO
      !$OMP END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Final Record on File 1. Calculate NINTEGt.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !$OMP CRITICAL
      NINTEGt = (NREC-1)*NINTMX
      !$OMP END CRITICAL

      !$OMP BARRIER

      !$OMP CRITICAL
      CALL LASTQOUTl(BUFP,IX,BUFP2,IX2,BUFPC,IXC,NINTEGtm,NINTMX,       &
                     NINTEGt,NREC,ICOUNT,IDONTW)
      !$OMP END CRITICAL
      !$OMP END PARALLEL

      CALL FINALIZEl(BUFPC, IXC, BUFP2, IX2, NINTEGtm, NINTMX, NINTEGt, &
                     NREC, IDONTW, IPRINTOPT)
      NREC = INT(NINTEGt/NINTMX) + 1
!-----------------------------------------------------------------------
      DEALLOCATE(BUFPC,IXC)
      RETURN
      END

! QOUTl
      SUBROUTINE QOUTl(BUFP,IX,BUFP2,IX2,NINTEGtm,NINTMX,IDONTW,        &
                       NSHELL,BLK,DI,DJ,DK,DL,ISH,JSH,KSH,LSH,          &
                       LOC, NREC, ICOUNT, CUTOFF)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NINTEGtm, NINTMX, IDONTW, NSHELL
      INTEGER, INTENT(IN) :: DI, DJ, DK, DL
      INTEGER, INTENT(IN) :: ISH, JSH, KSH, LSH
      INTEGER, INTENT(IN) :: LOC(NSHELL)
      DOUBLE PRECISION, INTENT(IN) :: BLK(DI, DJ, DK, DL)
      INTEGER(8), INTENT(OUT) :: IX(NINTMX)
      DOUBLE PRECISION, INTENT(OUT) :: BUFP(NINTMX)
      INTEGER(8), INTENT(OUT) :: IX2(NINTEGtm)
      DOUBLE PRECISION, INTENT(OUT) :: BUFP2(NINTEGtm)

      INTEGER :: I, J, K, L
      INTEGER(8) :: I1, I2, I3, I4, N, LABEL
      INTEGER LOCI, LOCJ, LOCK, LOCL
      INTEGER :: IJ_COUNT, KL_COUNT, DJEFF, DKEFF
      INTEGER :: NXInteg, IJBUFi, ibuf
      DOUBLE PRECISION :: VAL
      DOUBLE PRECISION :: CUTOFF
      INTEGER :: ICOUNT, NREC
      LOGICAL :: SAME, IANDJ, KANDL
!-----------------------------------------------------------------------
!     Pack 4-indices into 1 word.
!     Write Label & integral on Unit 1 if DONTW = .False.
!-----------------------------------------------------------------------
      SAME  = ISH == KSH .and. JSH == LSH
      IANDJ = ISH == JSH
      KANDL = KSH == LSH

      LOCI = LOC(ISH)
      LOCJ = LOC(JSH)
      LOCK = LOC(KSH)
      LOCL = LOC(LSH)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IJ_COUNT = 0
      DO I = 1, DI
        DJEFF = MERGE(I, DJ, IANDJ)
        DO J = 1, DJEFF
          IJ_COUNT = IJ_COUNT + 1
          KL_COUNT = 0
          DO K = 1, DK
            DKEFF = MERGE(K, DL, KANDL)
            DO L = 1, DKEFF
              KL_COUNT = KL_COUNT + 1
              IF(SAME .AND. KL_COUNT > IJ_COUNT) CYCLE

              VAL = BLK(I, J, K, L)
              IF(ABS(VAL) < CUTOFF) CYCLE

              I1 = LOCI + I - 1
              I2 = LOCJ + J - 1
              I3 = LOCK + K - 1
              I4 = LOCL + L - 1

              IF (I1 < I2) THEN
                N = I1
                I1 = I2
                I2 = N
              END IF
              IF (I3 < I4) THEN
                N = I3
                I3 = I4
                I4 = N
              END IF
              IF(I1 < I3 .OR. (I1 == I3 .AND. I2 < I4)) THEN
                N = I1
                I1 = I3
                I3 = N
                N = I2
                I2 = I4
                I4 = N
              END IF
!
              IF(I1 == I2) VAL = VAL * 0.5D0
              IF(I3 == I4) VAL = VAL * 0.5D0
              IF(I1 == I3 .and. I2 == I4) VAL = VAL * 0.5D0
!
              LABEL = ISHFT(I1, 48) + ISHFT(I2, 32) +        &
                  ISHFT(I3, 16) + I4

              IX(ICOUNT) = LABEL
              BUFP(ICOUNT) = VAL
              ICOUNT = ICOUNT + 1

              IF(ICOUNT > NINTMX) THEN
                NXInteg = NINTMX
                !$OMP CRITICAL
                IF(IDONTW==1)THEN
                  IJBUFi = (NREC-1)*NINTMX
                  DO ibuf=1,NINTMX
                    IX2  (IJBUFi+ibuf) = IX(ibuf)
                    BUFP2(IJBUFi+ibuf) = BUFP(ibuf)
                  END DO
                ELSE
                  WRITE(1)NXInteg,IX,BUFP
                END IF
                NREC = NREC+1
                !$OMP END CRITICAL
                ICOUNT = 1
              END IF

            END DO
          END DO
        END DO
      END DO
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE

! LASTQOUTl
      SUBROUTINE LASTQOUTl(BUFP,IX,BUFP2,IX2,BUFPC,IXC,NINTEGtm,        &
                           NINTMX,NINTEGt,NREC,ICOUNT,IDONTW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER :: NINTEGtm,NINTMX,NINTEGt,IDONTW
      INTEGER :: ICOUNT,NREC
      INTEGER(8)  :: IX(NINTMX)
      DOUBLE PRECISION :: BUFP(NINTMX)
      INTEGER(8) :: IXC(NINTMX)
      DOUBLE PRECISION :: BUFPC(NINTMX)
      INTEGER(8),DIMENSION(NINTEGtm) :: IX2
      DOUBLE PRECISION,DIMENSION(NINTEGtm) :: BUFP2
      INTEGER :: COUNTER, NXInteg, IBUF, NXInteg2, IJBUFi, jbuf, IDX
!-----------------------------------------------------------------------

      COUNTER = MOD(NINTEGt,NINTMX)
      NXInteg = ICOUNT-1
      IDX = COUNTER+1

      do ibuf=1,NXInteg
        IXC(IDX) = IX(ibuf)
        BUFPC(IDX) = BUFP(ibuf)
        NINTEGt = NINTEGt + 1
        IDX = IDX + 1

        IF(IDX - 1 .EQ. NINTMX) THEN
         NXInteg2 = NINTMX
         IF(IDONTW==1)THEN
          IJBUFi = (NREC-1)*NINTMX
          do jbuf=1,NINTMX
           IX2  (IJBUFi+jbuf) = IXC(jbuf)
           BUFP2(IJBUFi+jbuf) = BUFPC(jbuf)
          end do
         ELSE
          WRITE(1)NXInteg2,IXC,BUFPC
         END IF
         NREC = NREC+1
         IDX = 1
        END IF
      end do

      RETURN

      END

! FINALIZEl
      SUBROUTINE FINALIZEl(BUFPC,IXC,BUFP2,IX2,NINTEGtm,NINTMX,NINTEGt, &
                           NREC,IDONTW,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER :: NINTEGtm,NINTMX,NINTEGt,IDONTW,IPRINTOPT,NREC
      INTEGER(8),DIMENSION(NINTMX) :: IXC
      DOUBLE PRECISION,DIMENSION(NINTMX) :: BUFPC
      INTEGER(8),DIMENSION(NINTEGtm) :: IX2
      DOUBLE PRECISION,DIMENSION(NINTEGtm) :: BUFP2
      INTEGER :: NXInteg, IJBUFi, ibuf
!-----------------------------------------------------------------------

      NXInteg = MOD(NINTEGt,NINTMX)
      IF(NXInteg>=0) THEN
       IF(IDONTW==1)THEN
        IJBUFi = (NREC-1)*NINTMX
        do ibuf=1,NXInteg
         IX2  (IJBUFi+ibuf) = IXC(ibuf)
         BUFP2(IJBUFi+ibuf) = BUFPC(ibuf)
        end do
       ELSE
        NXInteg = -NXInteg
        WRITE(1)NXInteg,IXC,BUFPC
       END IF
      ELSE
       NINTEGt = NXInteg
      ENDIF
      IF(IPRINTOPT==1)WRITE(6,10)NINTEGt,NREC
      RETURN
!-----------------------------------------------------------------------
   10 FORMAT(I20,' 2e- integrals in ',I5,' records')
!-----------------------------------------------------------------------
      END

!----------------------------------------------------------------------!
!                                                                      !
!       2020  RI Approximation for ERIs implemented by                 !
!             Juan Felipe Huan Lew Yee and Jorge Martin del Campo      !
!                                                                      !
!             ( J. Chem. Phys. 154, 064102, 2021 )                     !
!                                                                      !
!----------------------------------------------------------------------!

! AUXREADl
      SUBROUTINE AUXREADl(IGTYP,NAT,NSHELLaux,NUMaux,NATmax,ANAM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
      COMMON/EXCaux/EXaux(2000),Caux(2000)
      COMMON/BASIS_FILE/BASIS_FILE
      LOGICAL :: FILE_EXISTS
      CHARACTER(80) :: BASIS_FILE,PREFIX,AUX_FILE,AUX_FILE_1
      CHARACTER(8),DIMENSION(NATmax) :: ANAM
      DOUBLE PRECISION :: COEFICIENT
      INTEGER :: EXT_POS
      DOUBLE PRECISION,PARAMETER :: PT2953=29.53125D0
      DOUBLE PRECISION,PARAMETER :: PT1624=162.421875D0
      DOUBLE PRECISION,PARAMETER :: PT75=0.75D0
      DOUBLE PRECISION,PARAMETER :: PT187=1.875D0
      DOUBLE PRECISION,PARAMETER :: TM10=1.0D-10
      DOUBLE PRECISION,PARAMETER :: PT6562=6.5625D0
      CHARACTER(8) :: BLANK
      DATA BLANK /'        '/
      CHARACTER(8) :: CBASIS
      CHARACTER(8),DIMENSION(8) :: LABEL
      DATA LABEL/'S       ','P       ','D       ','F       ',           &
                 'G       ','H       ','I       ','L       '/
      CHARACTER(8) :: BASIS
      INTEGER :: IGTYP

      PI = 2.0d0*DASIN(1.0d0)
      PI32 = PI*SQRT(PI)
      NGAUSS=0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Open Basis Set file if exists
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CBASIS = BLANK
      EXT_POS = SCAN(TRIM(ADJUSTL(BASIS_FILE)),".", BACK= .TRUE.)  !Look for extension position
      AUX_FILE = TRIM(ADJUSTL(BASIS_FILE(1:EXT_POS-1))) !Remove extension
      IF(LEN_TRIM(AUX_FILE)>0)THEN
       CALL GETENV( 'HOME', PREFIX )
       AUX_FILE_1 = TRIM(PREFIX)//'/DoNOFsw/basis/'//                 &
                      TRIM(ADJUSTL(AUX_FILE))//"-jkfit.bas"
       INQUIRE(FILE=AUX_FILE_1,EXIST=FILE_EXISTS)
       if(FILE_EXISTS)then                                ! in DoNOFsw
        AUX_FILE = AUX_FILE_1
!      DoNOF
       else
        AUX_FILE_1= TRIM(PREFIX)//'/DoNOF/basis/'//                   &
                      TRIM(ADJUSTL(AUX_FILE))//"-jkfit.bas"
        INQUIRE(FILE=AUX_FILE_1,EXIST=FILE_EXISTS)
        if(FILE_EXISTS)then                               ! in DoNOF
         AUX_FILE = AUX_FILE_1
        else
         CALL GETENV( 'PWD', PREFIX )
         AUX_FILE_1 = TRIM(PREFIX)//"/"//                             &
                        TRIM(ADJUSTL(AUX_FILE))//"-jkfit.bas"
         INQUIRE(FILE=AUX_FILE_1,EXIST=FILE_EXISTS)
         if(FILE_EXISTS)then                              ! in pwd
          AUX_FILE = AUX_FILE_1
         else
          AUX_FILE_1 = TRIM(ADJUSTL(AUX_FILE))//"-jkfit.bas"
          INQUIRE(FILE=AUX_FILE_1,EXIST=FILE_EXISTS)
          if(FILE_EXISTS)then                             ! in given path
           AUX_FILE = AUX_FILE_1
          else                             ! auxiliar basis set file not found
           WRITE(6,*)"Basis File ",TRIM(AUX_FILE)//"-jkfit.bas"," does not exist"
           WRITE(*,*) "Plese provide the file or use auxgen auxiliary basis"
           CALL ABRT
          endif
         endif
        endif
       endif
!      Open Auxiliar Basis Set File ( Unit = 50 )
       CLOSE(50) !Close original Basis Set file
       OPEN(50,FILE=AUX_FILE,STATUS='UNKNOWN',                        &
               FORM='FORMATTED',ACCESS='SEQUENTIAL')
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Read Auxilairy Basis from the basisname-jkfit.bas
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NUMaux = 0
      NSHELLaux = 0
      DO IAT=1,NAT
       if(LEN_TRIM(AUX_FILE)>0)then
        REWIND(50)
        CALL FNDATMBASIS(ANAM(IAT),IEOF)
       endif
    2  CONTINUE
       IEOF = 0
       IERR = 0
       CALL RDCARD(50,'$DATA 6U',IEOF)
       KSIZE = -8
       CALL GSTRNG(CBASIS,KSIZE)
       READ(UNIT=CBASIS,FMT='(A8)')BASIS
       IGAUSS = IFIND('NGAUSS  ',IERR)
!!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!     Read shell information
!!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(BASIS/=BLANK)THEN
!!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        KTYP = 0
        DO I=1,7
         IF(BASIS==LABEL(I))KTYP=I
        ENDDO
        IF(KTYP==0) THEN
        WRITE(*,*) 'Stop: Illegal auxiliary basis function type ',BASIS
        CALL ABRT
        END IF
!- - - - - - -
        NSHELLaux = NSHELLaux + 1
        KSTARTaux(NSHELLaux) = NGAUSS+1
        KATOMaux(NSHELLaux) = IAT
        KTYPEaux(NSHELLaux) = KTYP
        KNGaux(NSHELLaux) = IGAUSS
        KLOCaux(NSHELLaux) = NUMaux+1
        L = KTYP-1
        NGAUSS = NGAUSS + IGAUSS
        IF(IGTYP==1) THEN
          NUMaux = NUMaux + (L + 1) * (L + 2) / 2
        ELSE IF(IGTYP==2) THEN
          NUMaux = NUMaux + 2 * L + 1
        END IF
        K1 = KSTARTaux(NSHELLaux)
        K2 = K1 + KNGaux(NSHELLaux) - 1
        ! Read exponents and coefficients of primitives
        DO K = K1,K2
         C1 = 0.0D0
         IEOF = 0
         IERR = 0
         CALL RDCARD(50,'$DATA 7U',IEOF)
         IDUM = IFIND('IDUM    ',IERR)
         IF(IERR/=0)CALL ABRT
         EXaux(K) = RFIND('ZETA    ',IERR)
         IF(IERR/=0) CALL ABRT
         C1 = RFIND('C1      ',IERR)
         IF(IERR/=0) CALL ABRT
         !CALL NORMALIZE_AUXILIAR(EXaux(K),COEFICIENT,L)
         Caux(K) = C1!*COEFICIENT
        END DO
        ! Compute contracted normalization constant
        !FACL = 0.0D0
        !DO IG = K1,K2
        ! DO JG = K1,IG
        !  EE = EXaux(IG)+EXaux(JG)
        !  FAC = EE*SQRT(EE)
        !  IF(L==0) DUM = Caux(IG)*Caux(JG)/FAC
        !  IF(L==1) DUM = 0.5D0*Caux(IG)*Caux(JG)/(EE*FAC)
        !  IF(L==2) DUM = PT75  *Caux(IG)*Caux(JG)/(EE*EE*FAC)
        !  IF(L==3) DUM = PT187 *Caux(IG)*Caux(JG)/(EE**3*FAC)
        !  IF(L==4) DUM = PT6562*Caux(IG)*Caux(JG)/(EE**4*FAC)
        !  IF(L==5) DUM = PT2953*Caux(IG)*Caux(JG)/(EE**5*FAC)
        !  IF(L==6) DUM = PT1624*Caux(IG)*Caux(JG)/(EE**6*FAC)
        !  IF(IG /= JG) THEN
        !   DUM = DUM+DUM
        !  END IF
        !  FACL = FACL+DUM
        ! END DO
        !END DO
        !IF(FACL < TM10) THEN
        ! FACL = 0.0D0
        !ELSE
        ! FACL = 1.0D0/SQRT(FACL*PI32)
        !END IF
        !DO K = K1,K2
        ! Caux(K) = Caux(K) * FACL
        !END DO
        GOTO 2
       END IF
      END DO
      CLOSE(UNIT=50)
!     Open Basis Set File ( Unit = 50 )
      OPEN(50,FILE=BASIS_FILE,STATUS='UNKNOWN',                         &
              FORM='FORMATTED',ACCESS='SEQUENTIAL')

      RETURN
      END

! NORMALIZE_AUXILIAR
      SUBROUTINE NORMALIZE_AUXILIAR(EX,COEFICIENT,L)
      IMPLICIT NONE
      INTEGER :: L
      DOUBLE PRECISION :: EX,COEFICIENT
      DOUBLE PRECISION,PARAMETER :: PT2953=29.53125D0
      DOUBLE PRECISION,PARAMETER :: PT1624=162.421875D0
      DOUBLE PRECISION,PARAMETER :: PT75=0.75D0
      DOUBLE PRECISION,PARAMETER :: PT187=1.875D0
      DOUBLE PRECISION,PARAMETER :: PT6562=6.5625D0
      DOUBLE PRECISION :: PI,PI32,EE,FACS,FACP,FACD,FACF,FACG,FACH,FACI
      PI = 2.0d0*DASIN(1.0d0)
      PI32 = PI*SQRT(PI)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Normalize Primitive Basis Functions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      EE = EX+EX
      FACS = PI32/(EE*SQRT(EE))
      FACP = 0.5D0*FACS/EE
      FACD = PT75  *FACS/(EE*EE)
      FACF = PT187 *FACS/(EE**3)
      FACG = PT6562*FACS/(EE**4)
      FACH = PT2953*FACS/(EE**5)
      FACI = PT1624*FACS/(EE**6)
      IF(L==0) COEFICIENT = 1/SQRT(FACS)
      IF(L==1) COEFICIENT = 1/SQRT(FACP)
      IF(L==2) COEFICIENT = 1/SQRT(FACD)
      IF(L==3) COEFICIENT = 1/SQRT(FACF)
      IF(L==4) COEFICIENT = 1/SQRT(FACG)
      IF(L==5) COEFICIENT = 1/SQRT(FACH)
      IF(L==6) COEFICIENT = 1/SQRT(FACI)
!-----------------------------------------------------------------------
      END

! AUXGENl
      SUBROUTINE AUXGENl(IGTYP,NAT,NPRIMI,ITYP,IMIN,IMAX,NSHELLaux,     &
                         NUMaux,EX,ZAN,Cxyz)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      INTEGER :: NPRIMI,NSHELLaux,NUMaux,IGTYP
      INTEGER,DIMENSION(NPRIMI) :: ITYP
      INTEGER,DIMENSION(NAT) :: IMIN,IMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
!
      LOGICAL SMCD
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
      COMMON/EXCaux/EXaux(2000),Caux(2000)
!
      INTEGER,DIMENSION(NAT)::LMAX
      DOUBLE PRECISION,DIMENSION(NAT)::EXMAX,EXMIN
      INTEGER::IAT,N,MINI,MAXI,lbmax,nblock,iblock,nauxis,ifa,L
!-----------------------------------------------------------------------
!     AUXGEN generates the GEN-An* auxiliary basis
!-----------------------------------------------------------------------
      do IAT=1,NAT
       MINI = IMIN(IAT)
       MAXI = IMAX(IAT)
       LMAX(IAT)  = MAXVAL(ITYP(MINI:MAXI)) - 1
       EXMAX(IAT) = MAXVAL(EX(MINI:MAXI))
       EXMIN(IAT) = MINVAL(EX(MINI:MAXI))
      end do
!
      NUMaux = 0
      NSHELLaux = 0
      do IAT=1,NAT
       lbmax = LMAX(IAT)
       zbmax = EXMAX(IAT)
       zbmin = EXMIN(IAT)
       r = 6.0 - IGEN
       nblock = 2
       if(ISTAR==1) nblock = nblock + 1
!      H and He
       if (ZAN(IAT)<=2) then
         r = r - 0.5*IGEN + 2.0
         nblock = nblock - 1
       end if
!      How many z?
       N = int(log(zbmax/zbmin)/log(r) + 0.5)
       zmin = zbmin
       zmax= zmin*r**(N-1)
       zmax= 2.0*zmax
       TMPEXP = zmax*r
       do iblock=1,nblock
        if (iblock==1) then
          nauxis = max(1,N/nblock + mod(N,nblock))
        else
          nauxis = max(1,N/nblock)
        end if
        do ifa=1,nauxis
         TMPEXP = TMPEXP/r
         if (ifa==1) TMPEXP = TMPEXP*(r+0.5*IGEN)/r
         do L=0,2*(iblock-1)
          NSHELLaux = NSHELLaux + 1
          KATOMaux(NSHELLaux) = IAT
          KTYPEaux(NSHELLaux) = L+1
          KNGaux(NSHELLaux) = 1
          KLOCaux(NSHELLaux) = NUMaux+1
          IF(IGTYP==1) THEN
            NUMaux = NUMaux + (L + 1) * (L + 2) / 2
          ELSE IF(IGTYP==2) THEN
            NUMaux = NUMaux + 2 * L + 1
          END IF
          EXaux(NSHELLaux) = TMPEXP
!         avoiding warning
          Cxyz(1:3,IAT) = Cxyz(1:3,IAT)
          !CALL NORMALIZE_AUXILIAR(TMPEXP,COEFICIENT,L)
          Caux(NSHELLaux) = 1.0D0
         end do
         if (ifa==1) TMPEXP = TMPEXP*r/(r+0.5*IGEN)
        end do
       end do
      end do
!-----------------------------------------------------------------------
      RETURN
      END

! AuxERIl
      SUBROUTINE AuxERIl(NINTEGtm, BUFP2, NBF, IPRINTOPT, NSHELL, NAT,  &
                         SIZE_ENV, ENV, ATM, NBAS, BAS, IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      INTEGER, INTENT(IN) :: NINTEGtm, NBF, IPRINTOPT, NSHELL, NAT, SIZE_ENV, NBAS
      INTEGER, INTENT(IN) :: IGTYP
      DOUBLE PRECISION, INTENT(INOUT) :: BUFP2(NINTEGtm)
      DOUBLE PRECISION, INTENT(IN) :: ENV(SIZE_ENV)
      INTEGER, INTENT(IN) :: ATM(6, NAT), BAS(8, NBAS)
!
      INTEGER :: DI,DJ,DK
      INTEGER(4) :: SHLS(4)
      INTEGER :: NBFaux, NSHELLaux
      INTEGER :: ISH, JSH, KSH
      INTEGER :: LOC(NSHELL), LOCAUX(NSHELLaux)
      INTEGER :: Dcgto(NSHELL), DAUXcgto(NSHELLaux)
      DOUBLE PRECISION, ALLOCATABLE :: GMAT(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: BLK(:,:,:)
      INTEGER :: ERR

      INTEGER, EXTERNAL :: CINTcgto_spheric, cint3c2e_sph
      INTEGER, EXTERNAL :: CINTcgto_cart, cint3c2e_cart
!-----------------------------------------------------------------------
! LOC indicates the starting location of each shell in the AO basis
!-----------------------------------------------------------------------
      LOC(1) = 1
      IF(IGTYP==1) DI = CINTcgto_cart(0, BAS)
      IF(IGTYP==2) DI = CINTcgto_spheric(0, BAS)
      Dcgto(1) = DI
      DO ISH = 2,NSHELL
        LOC(ISH) = LOC(ISH-1) + DI
        IF(IGTYP==1) DI = CINTcgto_cart(ISH-1, BAS)
        IF(IGTYP==2) DI = CINTcgto_spheric(ISH-1, BAS)
        Dcgto(ISH) = DI
      END DO

      LOCAUX(1) = 1
      IF(IGTYP==1) DI = CINTcgto_cart(NSHELL, BAS)
      IF(IGTYP==2) DI = CINTcgto_spheric(NSHELL, BAS)
      DAUXcgto(1) = DI
      DO ISH = 2,NSHELLaux
        LOCAUX(ISH) = LOCAUX(ISH-1) + DI
        IF(IGTYP==1) DI = CINTcgto_cart(NSHELL+ISH-1, BAS)
        IF(IGTYP==2) DI = CINTcgto_spheric(NSHELL+ISH-1, BAS)
        DAUXcgto(ISH) = DI
      END DO
!-----------------------------------------------------------------------
      BUFP2 = 0.0D0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Evaluate G = (P|Q) and G^{-1/2}
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(GMAT(NBFaux,NBFaux))
      GMAT = 0.0D0
      CALL METRICmatl(GMAT, SIZE_ENV, ENV, NAT, ATM, NBAS, BAS,         &
                     DAUXcgto, LOCaux, NSHELL, NSHELLaux, NBFaux, IGTYP)
      CALL PDPT_msqrtl(GMAT,NBFaux,IPRINTOPT)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     2e- Integrals (S,P,D,F,G & L Shells)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !$OMP PARALLEL PRIVATE(ISH, JSH, KSH, LQSUM, DI, DJ, DK, SHLS,    &
      !$OMP ERR, BLK, AUX, NORGH)
      !$OMP DO SCHEDULE(DYNAMIC)
      DO ISH = 1,NSHELL
        DI = Dcgto(ISH)
        SHLS(1) = ISH - 1
        DO JSH = 1,ISH
          DJ = Dcgto(JSH)
          SHLS(2) = JSH - 1
          DO KSH = 1,NSHELLaux
            DK = DAUXcgto(KSH)
            SHLS(3) = KSH - 1 + NSHELL
!- - - - - - - - - - - - - - - - - - - - - - - - - -
!           (II,JJ//KK)
!           Compute 2e- Integrals (mn|k)
!           Select integral code for ERI calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ALLOCATE(BLK(DI,DJ,DK))
            IF(IGTYP==1)ERR = cint3c2e_cart(BLK,SHLS,ATM,NAT,BAS,     &
                                            NBAS, ENV, 0_8)
            IF(IGTYP==2)ERR = cint3c2e_sph(BLK,SHLS,ATM,NAT,BAS,      &
                                           NBAS, ENV, 0_8)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!           Contract to B tensor
!- - - - - - - - - - - - - - - - - - - - - - - - - - - -
            CALL QOUT3Cl(BUFP2, GMAT, NBF, NBFaux, BLK, DI, DJ, DK,     &
                         ISH, JSH, KSH, LOC, LOCAUX, NSHELL, NSHELLaux)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - -
            DEALLOCATE(BLK)
          END DO
        END DO
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
      DEALLOCATE(GMAT)
!-----------------------------------------------------------------------
      RETURN
      END

! AuxERIModCholl
      SUBROUTINE AuxERIModCholl(NINTEGtm,BUFP2,NBF,IPRINTOPT,NSHELL,    &
                                NAT,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      INTEGER, INTENT(IN) :: NINTEGtm, NBF, IPRINTOPT, NSHELL, NAT, SIZE_ENV, NBAS
      INTEGER, INTENT(IN) :: IGTYP
      DOUBLE PRECISION, INTENT(INOUT) :: BUFP2(NBF*(NBF+1)/2,NBFaux)
      DOUBLE PRECISION, INTENT(IN) :: ENV(SIZE_ENV)
      INTEGER, INTENT(IN) :: ATM(6, NAT), BAS(8, NBAS)
!
      INTEGER :: DI,DJ,DK
      INTEGER(4) :: SHLS(4)
      INTEGER :: NBFaux, NSHELLaux
      INTEGER :: ISH, JSH, KSH
      INTEGER :: LOC(NSHELL), LOCAUX(NSHELLaux)
      INTEGER :: Dcgto(NSHELL), DAUXcgto(NSHELLaux)
      INTEGER :: ERR
      INTEGER,ALLOCATABLE :: IPIV(:)
      DOUBLE PRECISION,ALLOCATABLE :: E(:)
      DOUBLE PRECISION, ALLOCATABLE :: GMAT(:,:), BUFP(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: L(:,:), D(:,:), P(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: BLK(:,:,:)

      INTEGER, EXTERNAL :: CINTcgto_spheric, cint3c2e_sph
      INTEGER, EXTERNAL :: CINTcgto_cart, cint3c2e_cart
!-----------------------------------------------------------------------
! LOC indicates the starting location of each shell in the AO basis
!-----------------------------------------------------------------------
      LOC(1) = 1
      IF(IGTYP==1) DI = CINTcgto_cart(0, BAS)
      IF(IGTYP==2) DI = CINTcgto_spheric(0, BAS)
      Dcgto(1) = DI
      DO ISH = 2,NSHELL
        LOC(ISH) = LOC(ISH-1) + DI
        IF(IGTYP==1) DI = CINTcgto_cart(ISH-1, BAS)
        IF(IGTYP==2) DI = CINTcgto_spheric(ISH-1, BAS)
        Dcgto(ISH) = DI
      END DO

      LOCAUX(1) = 1
      IF(IGTYP==1) DI = CINTcgto_cart(NSHELL, BAS)
      IF(IGTYP==2) DI = CINTcgto_spheric(NSHELL, BAS)
      DAUXcgto(1) = DI
      DO ISH = 2,NSHELLaux
        LOCAUX(ISH) = LOCAUX(ISH-1) + DI
        IF(IGTYP==1) DI = CINTcgto_cart(NSHELL+ISH-1, BAS)
        IF(IGTYP==2) DI = CINTcgto_spheric(NSHELL+ISH-1, BAS)
        DAUXcgto(ISH) = DI
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Evaluate G = (P|Q), G = PLDL^TP^T, ModChol and get P, L, D^1/2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(BUFP(NBF*(NBF+1)/2,NBFaux))
      ALLOCATE(GMAT(NBFaux,NBFaux))
      ALLOCATE(L(NBFaux,NBFaux),D(NBFaux,NBFaux),P(NBFaux,NBFaux))
      ALLOCATE(E(NBFaux),IPIV(NBFaux))
      GMAT = 0.0D0
      CALL METRICmatl(GMAT, SIZE_ENV, ENV, NAT, ATM, NBAS, BAS,         &
                     DAUXcgto, LOCaux, NSHELL, NSHELLaux, NBFaux, IGTYP)
      CALL LDLT(GMAT,IPIV,E,NBFaux)
      CALL MODCHOL(NBFaux, GMAT, IPIV, E, 1D-10,IPRINTOPT)
      CALL GET_PLD12(GMAT,IPIV,E,NBFaux,P,L,D)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     2e- Integrals (S,P,D,F,G & L Shells)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !$OMP PARALLEL PRIVATE(ISH, JSH, KSH, LQSUM, DI, DJ, DK, SHLS,    &
      !$OMP ERR, BLK, AUX, NORGH)
      !$OMP DO SCHEDULE(DYNAMIC)
      DO ISH = 1,NSHELL
        DI = Dcgto(ISH)
        SHLS(1) = ISH - 1
        DO JSH = 1,ISH
          DJ = Dcgto(JSH)
          SHLS(2) = JSH - 1
          DO KSH = 1,NSHELLaux
            DK = DAUXcgto(KSH)
            SHLS(3) = KSH - 1 + NSHELL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!           (II,JJ//KK)
!           Compute 2e- Integrals (mn|k)
!           Select integral code for ERI calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ALLOCATE(BLK(DI,DJ,DK))
            IF(IGTYP==1) ERR = cint3c2e_cart(BLK, SHLS, ATM, NAT, BAS,&
                                             NBAS, ENV, 0_8)
            IF(IGTYP==2) ERR = cint3c2e_sph(BLK, SHLS, ATM, NAT, BAS, &
                                            NBAS, ENV, 0_8)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!           Store 3 center ERIs (mn|k)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            CALL QOUT3CModCholl(BUFP,GMAT,NBF,NBFaux,BLK,DI,DJ,DK,ISH,  &
                                JSH,KSH,LOC,LOCAUX,NSHELL,NSHELLaux)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            DEALLOCATE(BLK)
          END DO
        END DO
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
!-----------------------------------------------------------------------
!     Build b tensor, solve linear equation system LD^1/2 b = P^T(k|mn)
!-----------------------------------------------------------------------
      BUFP2 = 0.0D0
      !$OMP PARALLEL DO PRIVATE(I,J,K)
      DO I=1,NBF*(NBF+1)/2
        DO J=1,NBFaux
          DO K=1,NBFaux
            BUFP2(I,K) = BUFP2(I,K) + BUFP(I,J) * P(J,K)
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
      DEALLOCATE(BUFP)
      DO N=1,NBF
        DO M=1,N
          MN = M + N*(N-1)/2
          CALL DTRTRS('L','N','U',NBFaux,1,L,NBFaux,BUFP2(MN,1:NBFaux), &
                      NBFaux,INFO)
          CALL SOLVE_BLOCK_SYSTEM(NBFaux,GMAT,BUFP2(MN,1:NBFaux),E)
        END DO
      END DO
      DEALLOCATE(GMAT)
      DEALLOCATE(L,D,P)
      DEALLOCATE(E,IPIV)
!-----------------------------------------------------------------------
      RETURN
      END

! LDLT
      SUBROUTINE LDLT(G,ipiv,e,n)
      integer :: n
      real(8) :: G(n,n)
      integer :: ipiv(n)
      real(8) :: e(n)
      real(8) :: work(3*n-1)
      integer :: lwork
      integer :: info
!-----------------------------------------------------------------------
      lwork = 3*n-1
      ipiv(:) = 0
      call dsytrf_rk('L',n,G,n,e,ipiv,work,lwork,info)
!-----------------------------------------------------------------------
      END

! MODCHOL
      SUBROUTINE MODCHOL ( n, a, ipiv, diag, delta, iprintopt )
      implicit none
      integer    n, iprintopt
      integer ipiv(n)
      real(8) a(n,n), diag(n), delta
!
      real(8) t(2,2), eig(2), work(20), zero
      parameter (zero = 0.0d0)
      integer k, lwork, info, nchanged
!-----------------------------------------------------------------------
      lwork = 10
      k = 1
      nchanged = 0

      do while ( k.le.n )

!     Check 1x1 block
        if ( k.eq.n .or. diag(k).eq. zero ) then
          if ( a(k,k).le.delta) then
             nchanged = nchanged + 1
             a(k,k) = delta
          end if
          a(k,k) = sqrt(a(k,k))
          k = k + 1
        else

!     Check 2x2 block

          t(1,1) = a(k,k)
          t(2,1) = diag(k)
          t(1,2) = t(2,1)
          t(2,2) = a(k+1,k+1)

          call dsyev ( 'v', 'l', 2, t, 2, eig, work, lwork, info )

          if (eig(1).le.delta) then
            nchanged = nchanged + 1
            eig(1) = delta
          end if
          eig(1) = sqrt(eig(1))
          if (eig(2).le.delta) then
            nchanged = nchanged + 1
            eig(2) = delta
          end if
          eig(2) = sqrt(eig(2))
!
          a(k,k)     = t(1,1)*t(1,1)*eig(1) + t(1,2)*t(1,2)*eig(2)
          a(k+1,k+1) = t(2,1)*t(2,1)*eig(1) + t(2,2)*t(2,2)*eig(2)
          diag(k)    = t(1,1)*t(2,1)*eig(1) + t(1,2)*t(2,2)*eig(2)

!   Check off-diagonal element; if zero reset 2x2 pivot to 2(1x1) pivots

          if ( diag(k).eq.zero ) then
             ipiv(k)   = - ipiv(k)
             ipiv(k+1) = - ipiv(k+1)
          end if
!
          k = k + 2

        end if
      end do

      if(iprintopt==1.and.nchanged>0) THEN
       write(6,*)"SMCD Warning - Number of values changed from metric:",&
             nchanged
      end if

!-----------------------------------------------------------------------
      END

! GET_PLD12
      SUBROUTINE GET_PLD12(G,ipiv,e,n,P,L,D)
      integer :: n,i,j
      real(8) :: G(n,n)
      integer :: ipiv(n)
      real(8) :: e(n)
      real(8) :: L(n,n), D(n,n), vec(n)
      real(8) :: P(n,n)
!-----------------------------------------------------------------------
      P(:,:) = 0.0
      do i=1,n
        P(i,i) = 1.0
      end do
      L(:,:) = 0.0
      D(:,:) = 0.0

      do i=1,n
        do j=1,i-1
          L(i,j) = G(i,j)
        end do
        L(i,i) = 1.0
        D(i,i) = G(i,i)
        if(i.LT.n) then
          D(i+1,i) = e(i)
          D(i,i+1) = e(i)
        end if
      end do

      do i=n,1,-1
        if(abs(ipiv(i)).NE.i) then
          vec(:) = P(i,:)
          P(i,:) =P(abs(ipiv(i)),:)
          P(abs(ipiv(i)),:) = vec(:)
        end if
      end do
!-----------------------------------------------------------------------
      END

! SOLVE_BLOCK_SYSTEM
      SUBROUTINE SOLVE_BLOCK_SYSTEM( n, a, y , diag )
      implicit none
      integer n
      real(8) a(n,n), diag(n), y(n), x(n)
!
      real(8) t(2,2), zero
      parameter (zero = 0.0d0)
      integer k
!-----------------------------------------------------------------------
      k = 1

      do while ( k.le.n )

!     Check 1x1 block
        if ( k.eq.n .or. diag(k).eq. zero ) then
          x(k) = y(k)/a(k,k)
          k = k + 1
        else

!     Check 2x2 block

          t(1,1) = a(k,k)
          t(2,1) = diag(k)
          t(1,2) = t(2,1)
          t(2,2) = a(k+1,k+1)

          if(abs(t(1,1)*t(2,2)-t(1,2)*t(2,1)).le.zero)    &
                  WRITE(6,*) "Warning: no solution available"

          x(k)=(y(k)*t(2,2)-t(1,2)*y(k+1))/(t(1,1)*t(2,2)-t(1,2)*t(2,1))
        x(k+1)=(t(1,1)*y(k+1)-y(k)*t(2,1))/(t(1,1)*t(2,2)-t(1,2)*t(2,1))

          k = k + 2

        end if
      end do

      y = x
!-----------------------------------------------------------------------
      END

! METRICmatl
      SUBROUTINE METRICmatl(GMAT, SIZE_ENV, ENV, NAT, ATM, NBAS, BAS,   &
                     DAUXcgto, LOCAUX, NSHELL, NSHELLAUX, NBFaux, IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER, INTENT(IN) :: SIZE_ENV, NAT, NBAS, NSHELLAUX, IGTYP
      DOUBLE PRECISION, INTENT(IN) :: ENV(SIZE_ENV)
      INTEGER, INTENT(IN) :: ATM(6, NAT), BAS(8, NBAS)
      DOUBLE PRECISION, INTENT(OUT) :: GMAT(NBFaux,NBFaux)
      INTEGER, INTENT(IN) :: LOCAUX(NSHELLAUX)
!
      INTEGER :: DI, DK
      INTEGER(4) :: SHLS(2)
      INTEGER :: NSHELL, ISH, KSH, ERR
      INTEGER :: DAUXcgto(NSHELLAUX)

      DOUBLE PRECISION, ALLOCATABLE :: BLK(:,:)

      INTEGER(4), EXTERNAL :: CINTcgto_cart, cint2c2e_cart, cint1e_ovlp_cart
      INTEGER(4), EXTERNAL :: CINTcgto_spheric, cint2c2e_sph, cint1e_ovlp_sph
!-----------------------------------------------------------------------
!     2e- Integrals (S,P,D,F,G & L Shells)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !$OMP PARALLEL &
      !$OMP PRIVATE(II, KK, ISH, KSH, LQSUM, BLK, DI, DK, ERR, SHLS)
      !$OMP DO SCHEDULE(DYNAMIC)
      DO ISH = 1,NSHELLAUX                          ! II Shell
       SHLS(1) = ISH-1 + NSHELL
       DI = DAUXcgto(ISH)
       DO KSH = 1,ISH                               ! KK Shell
        SHLS(2) = KSH-1 + NSHELL
        DK = DAUXcgto(KSH)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Compute 2e- Integrals
!       Select integral code for ERI calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ALLOCATE(BLK(DI,DK))
        IF(IGTYP==1) ERR = cint2c2e_cart(BLK, SHLS, ATM, NAT, BAS,    &
                NBAS, ENV, 0_8)
        IF(IGTYP==2) ERR = cint2c2e_sph(BLK, SHLS, ATM, NAT, BAS,     &
                NBAS, ENV, 0_8)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Transfer integrals to G matrix
!- - - - - - - - - - - - - - - - - - - - - - -
        CALL QOUT2Cl(GMAT,NBFaux,BLK,DI,DK,LOCAUX,NSHELLAUX,ISH,KSH)
!- - - - - - - - - - - - - - - - - - - - - - -
        DEALLOCATE(BLK)
       END DO
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
!-----------------------------------------------------------------------
      RETURN
      END

! PDPT_msqrtl
      SUBROUTINE PDPT_msqrtl(GMAT,N,IPRINTOPT)
      INTEGER :: N,IPRINTOPT,I,J,NTRUNC
      DOUBLE PRECISION :: GMAT(N,N),TOL
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: VEC,AUX,AUX2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: EIG,W

      ALLOCATE(VEC(N,N),AUX(N,N),AUX2(N,N))
      ALLOCATE(EIG(N),W(N))
!
      TOL = 1.0d-12
      AUX = GMAT
      CALL DIAG(N,AUX,VEC,EIG,W)
!
      NTRUNC = 0
      AUX = 0.0d0
      DO I=1,N
       IF (EIG(I)<=TOL) THEN
        NTRUNC = NTRUNC + 1
        EIG(I) = TOL!0.0D0
        !DO J=1,N
        ! VEC(J,I) = 0.0d0
        !END DO
       !ELSE
        !AUX(I,I) = 1.0/DSQRT(EIG(I))
       END IF
       AUX(I,I) = 1.0/DSQRT(EIG(I))
      END DO

      IF(IPRINTOPT==1.AND.NTRUNC>0)                                     &
      WRITE(6,*) "RI Warning - Number of values truncated from metric:",&
                 NTRUNC

      DEALLOCATE(EIG,W)
      CALL DGEMM("N","N",N,N,N,1.0D0,VEC,N,AUX,N,0.0D0,AUX2,N)
      CALL DGEMM("N","T",N,N,N,1.0D0,AUX2,N,VEC,N,0.0D0,GMAT,N)
      DEALLOCATE(VEC,AUX,AUX2)

      END

! QOUT2Cl
      SUBROUTINE QOUT2Cl(GMAT,NBFaux,BLK,DI,DK,LOCAUX,NSHELLAUX,ISH,KSH)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: DI, DK
      INTEGER, INTENT(IN) :: NBFaux,  NSHELLAUX, ISH, KSH
      INTEGER, INTENT(IN) :: LOCAUX(NSHELLAUX)
      DOUBLE PRECISION, INTENT(INOUT) :: GMAT(NBFaux, NBFaux)
      DOUBLE PRECISION, INTENT(IN) :: BLK(DI, DK)

      INTEGER :: I, K, I1, I3, LOCI, LOCK
      LOGICAL :: SAME

!-----------------------------------------------------------------------
!     Get Auxiliary Basis Info.
!-----------------------------------------------------------------------
      SAME = (ISH == KSH)
      LOCI = LOCAUX(ISH)
      LOCK = LOCAUX(KSH)

!-----------------------------------------------------------------------
!     Store (k|l) ERIs in G matrix
!-----------------------------------------------------------------------
      DO I = 1,DI
        DO K = 1,DK
          IF(SAME.and.K>I) CYCLE
          I1 = LOCI + I - 1
          I3 = LOCK + K - 1
          GMAT(I1,I3) = BLK(I, K)
          GMAT(I3,I1) = BLK(I, K)
        END DO
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

! QOUT3Cl
      SUBROUTINE QOUT3Cl(BUFP2, GMAT, NBF, NBFAUX, BLK, DI, DJ, DK,     &
                         ISH, JSH, KSH, LOC, LOCAUX, NSHELL, NSHELLAUX)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NBF, NBFAUX
      INTEGER, INTENT(IN) :: ISH, JSH, KSH, NSHELL, NSHELLAUX
      INTEGER, INTENT(IN) :: DI, DJ, DK
      INTEGER, INTENT(IN) :: LOC(NSHELL), LOCAUX(NSHELLAUX)
      DOUBLE PRECISION, INTENT(INOUT) :: BUFP2(NBF*(NBF+1)/2, NBFAUX)
      DOUBLE PRECISION, INTENT(IN) :: GMAT(NBFAUX, NBFAUX)
      DOUBLE PRECISION, INTENT(IN) :: BLK(DI, DJ, DK)


      LOGICAL :: IANDJ
      INTEGER :: DJEFF, I
      INTEGER :: J, K, L, M, N, KK, MN, T
      INTEGER :: LOCI, LOCJ, LOCK
      DOUBLE PRECISION :: VAL
!-----------------------------------------------------------------------
!     Get Basis and Auxiliary Basis Info.
!-----------------------------------------------------------------------
      IANDJ = ISH == JSH
      LOCI = LOC(ISH)
      LOCJ = LOC(JSH)
      LOCK = LOCAUX(KSH)

!-----------------------------------------------------------------------
!     Build B tensor (Only Upper Triangular)
!      b_mn^l = sum_k (mn|k) G^{-1/2}_{kl}
!-----------------------------------------------------------------------
      DO I = 1, DI
       DJEFF = MERGE(I, DJ, IANDJ)
       DO J = 1, DJEFF
        DO K = 1, DK
          M = LOCI + I - 1
          N = LOCJ + J - 1
          KK = LOCK + K - 1
          VAL = BLK(I, J, K)

          IF (M > N) THEN
            T = M
            M = N
            N = T
          END IF
          MN = M + N*(N-1)/2

          DO L=1,NBFaux
           BUFP2(MN, L) = BUFP2(MN, L) + VAL * GMAT(KK, L)
          END DO

        END DO
       END DO
      END DO
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE

! QOUT3CModCholl
      SUBROUTINE QOUT3CModCholl(BUFP2,GMAT,NBF,NBFAUX,BLK,DI,DJ,DK,ISH, &
                                JSH,KSH,LOC,LOCAUX,NSHELL,NSHELLAUX)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NBF, NBFAUX
      INTEGER, INTENT(IN) :: ISH, JSH, KSH, NSHELL, NSHELLAUX
      INTEGER, INTENT(IN) :: DI, DJ, DK
      INTEGER, INTENT(IN) :: LOC(NSHELL), LOCAUX(NSHELLAUX)
      DOUBLE PRECISION, INTENT(INOUT) :: BUFP2(NBF*(NBF+1)/2, NBFAUX)
      DOUBLE PRECISION, INTENT(IN) :: GMAT(NBFAUX, NBFAUX)
      DOUBLE PRECISION, INTENT(IN) :: BLK(DI, DJ, DK)

      LOGICAL :: IANDJ
      INTEGER :: DJEFF, I
      INTEGER :: J, K, L, M, N, KK, MN, T
      INTEGER :: LOCI, LOCJ, LOCK
      DOUBLE PRECISION :: VAL
!-----------------------------------------------------------------------
!     Get Basis and Auxiliary Basis Info.
!-----------------------------------------------------------------------
      IANDJ = ISH == JSH
      LOCI = LOC(ISH)
      LOCJ = LOC(JSH)
      LOCK = LOCAUX(KSH)

!-----------------------------------------------------------------------
!     Build B tensor (Only Upper Triangular)
!      b_mn^l = sum_k (mn|k) G^{-1/2}_{kl}
!-----------------------------------------------------------------------
      DO I = 1, DI
       DJEFF = MERGE(I, DJ, IANDJ)
       DO J = 1, DJEFF
        DO K = 1, DK
          M = LOCI + I - 1
          N = LOCJ + J - 1
          KK = LOCK + K - 1
          VAL = BLK(I, J, K)

          IF (M > N) THEN
            T = M
            M = N
            N = T
          END IF
          MN = M + N*(N-1)/2

          BUFP2(MN,KK) = VAL

        END DO
       END DO
      END DO
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------!
