!======================================================================!
!                                                                      !
!               I N T E G R A L   S U B R O U T I N E S                !
!                                                                      !
!======================================================================!
!======================================================================!

! OneElecInt                                           
      SUBROUTINE OneElecInt(Cxyz,H,S,TKIN,DInteg,NBFT,IPRINTOPT,        &
                            EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,KSTART,      &
                            KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,ZAN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/USELIBRETA/ILIBRETA      
      COMMON/INTOPT/ISCHWZ,IECP,NECP            
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB                        
      LOGICAL EFLDL                                                  
      COMMON/EFLDC_1/EFLDL
      COMMON/EFLDC_2/EVEC(3)
!      
      INTEGER :: NBFT,IPRINTOPT,NPRIMI,NSHELL
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NBFT) :: H,S,TKIN    
      DOUBLE PRECISION,DIMENSION(3*NBFT) :: DInteg
      INTEGER :: CR, CM, timestartoneE, timefinishoneE
      DOUBLE PRECISION :: RATE
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
      if(ILIBRETA==0)then
       CALL HSandT(H,S,TKIN,NBFT,EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,        &
                   KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,        &
                   ZAN,Cxyz)
      else if(ILIBRETA==1)then
       CALL HSandTlib(H,S,TKIN,NBFT,KATOM,KLOC,KMIN,KMAX,NSHELL,        &
                      ZAN,Cxyz)      
      end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ECP integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF( 1<=IECP .and. IECP<=3 )THEN
       CALL ECPINT(H,NBFT,EX,CS,CP,CD,CF,CG,NPRIMI,KSTART,KATOM,KNG,    &
                   KLOC,KMIN,KMAX,NSHELL,Cxyz)
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Add Electric Field Contribution to 1e Hamiltonian
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(EFLDL)THEN
       CALL DipInt(DInteg,NBFT,EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,          &
                   KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,        &
                   Cxyz,NAT)
       CALL ElecFieldInt(H,DInteg,NBFT)              
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL SYSTEM_CLOCK(timefinishoneE)
      DeltaToneE = (timefinishoneE - timestartoneE)/RATE
      IF(IPRINTOPT==1)                                                  &
       WRITE(6,'(1X,A22,F10.2)')'Time to do integrals =',DeltaToneE                                                                                                    
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               

! DipInt                                           
      SUBROUTINE DipInt(DInteg,NBFT,EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,     &
                        KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,   &
                        Cxyz,NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/USELIBRETA/ILIBRETA            
      COMMON/ELPROP/IEMOM      
      COMMON/TRANSF/XP,YP,ZP
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(3*NBFT) :: DInteg
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AUX
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
      if(ILIBRETA==0)then
       CALL PRCALC(DInteg,AUX,3,NBFT,EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,    &
                   KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,Cxyz)                   
      else if(ILIBRETA==1)then
       CALL PRCALClib(DInteg,AUX,3,NBFT,KATOM,KLOC,KMIN,KMAX,NSHELL)      
      end if
      IEMOM = IEMSV                                                  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DEALLOCATE(AUX)                                                                                  
      RETURN                                                            
      END                                                               

! ElecFieldInt
      SUBROUTINE ElecFieldInt(H,DInteg,NBFT)   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      LOGICAL EFLDL                                                  
      COMMON/EFLDC_1/EFLDL
      COMMON/EFLDC_2/EVEC(3) 
      DOUBLE PRECISION,DIMENSION(NBFT) :: H
      DOUBLE PRECISION,DIMENSION(3*NBFT) :: DInteg
!-----------------------------------------------------------------------
      WRITE(6,'(1X,A16,3F10.5)')'Electric Field =',(EVEC(I),I=1,3)                     
      IF(EVEC(1)/=0.0d0)CALL DAXPY(NBFT,EVEC(1),DInteg(1       ),1,H,1)              
      IF(EVEC(2)/=0.0d0)CALL DAXPY(NBFT,EVEC(2),DInteg(1+  NBFT),1,H,1)              
      IF(EVEC(3)/=0.0d0)CALL DAXPY(NBFT,EVEC(3),DInteg(1+2*NBFT),1,H,1)              
!----------------------------------------------------------------------- 
      RETURN                                                            
      END                                                               

! JandK                                            
      SUBROUTINE JandK(BUFP2,IX2,NINTEGtm,NINTEGt,NRECO,XINTS,NSH2,     &
                       IDONTW,IPRINTOPT,EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI, &
                       KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,    &
                       Cxyz,NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/USELIBRETA/ILIBRETA      
      COMMON/INTFIL/NINTMX           
      COMMON/INTOPT/ISCHWZ,IECP,NECP            
      COMMON/RESTAR/NREC,IST,JST,KST,LST           
!      
      LOGICAL SCHWRZ
      INTEGER :: NINTEGtm,NINTEGt,NRECO,NSH2,IDONTW,IPRINTOPT
      INTEGER :: NPRIMI,NSHELL,NAT
      INTEGER,DIMENSION(NINTEGtm) :: IX2
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NINTEGtm) :: BUFP2
      DOUBLE PRECISION,DIMENSION(NSH2) :: XINTS
!
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IX
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: BUFP,GHONDOLIB
      INTEGER :: CR, CM, timestarttwoE, timefinishtwoE
      DOUBLE PRECISION :: RATE
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
      CALL Debut(IDONTW,IPRINTOPT,KATOM,NSHELL,Cxyz)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Schwarz inequality
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(GHONDOLIB(MAXG))
      SCHWRZ = ISCHWZ>0 
      IF(SCHWRZ)THEN
       if(ILIBRETA==0)then
        CALL ExchangeInt(XINTS,GHONDOLIB,NSH2,MAXG,EX,CS,CP,CD,CF,CG,   &
                         CH,CI,NPRIMI,KSTART,KATOM,KTYPE,KNG,KLOC,KMIN, &
                         KMAX,NSHELL,Cxyz,NAT)                          
       else if(ILIBRETA==1)then
        CALL ExchangeIntlib(XINTS,GHONDOLIB,NSH2,MAXG,KTYPE,KMIN,KMAX,  &
                            NSHELL)
       end if
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     2e integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(BUFP(NINTMX),IX(NINTMX))
      CALL TwoERI(SCHWRZ,NINTEGtm,NINTEGt,NSCHWZ,BUFP,IX,BUFP2,IX2,     &
                  XINTS,NSH2,GHONDOLIB,MAXG,IDONTW,IPRINTOPT,EX,CS,     &
                  CP,CD,CF,CG,CH,CI,NPRIMI,KSTART,KATOM,KTYPE,KNG,      &
                  KLOC,KMIN,KMAX,NSHELL,Cxyz,NAT)
      DEALLOCATE(BUFP,IX,GHONDOLIB)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NRECO = NREC
      CALL SYSTEM_CLOCK(timefinishtwoE)
      DeltaTtwoE = (timefinishtwoE - timestarttwoE)/RATE
      IF(IPRINTOPT==1)                                                  &
       WRITE(6,'(1X,A22,F10.2)')'Time to do integrals =',DeltaTtwoE
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               

! Debut                                            
      SUBROUTINE Debut(IDONTW,IPRINTOPT,KATOM,NSHELL,Cxyz)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/INTFIL/NINTMX           
      COMMON/RESTAR/NREC,IST,JST,KST,LST           
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
      COMMON/SHLT/SHLTOL,CUTOFF,ICOUNT
      INTEGER :: IDONTW,IPRINTOPT,NSHELL
      INTEGER,DIMENSION(NSHELL) :: KATOM
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz   
      DOUBLE PRECISION,DIMENSION(NSHELL,3) :: CO
!-----------------------------------------------------------------------
!     Initialize 2e- integral Calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                                                                              
      SHLTOL = 20*2.30258D0   
      CUTOFF = 1.0D-09
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
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
!     Normal Start
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IST = 1                                              
      JST = 1                                              
      KST = 1                                              
      LST = 1                                              
      NREC   = 1                                                        
      ICOUNT = 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                                                               
!     Format Statements
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
   10 FORMAT(/1X,'2e- integrals'/1X,13(1H-))      
   20 FORMAT(' DONTW option skips storage 2e- integrals on Unit 1')      
   30 FORMAT(' Storing',I8,' integrals/record on disk, using',I3,       &
             ' Bytes/integral')                                        
!-----------------------------------------------------------------------                                       
      RETURN                                                            
      END                                                               

! TwoERI                                            
      SUBROUTINE TwoERI(SCHWRZ,NINTEGtm,NINTEGt,NSCHWZ,BUFP,IX,BUFP2,   &
                        IX2,XINTS,NSH2,GHONDOLIB,MAXG,IDONTW,IPRINTOPT, &
                        EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,KSTART,KATOM,    &
                        KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,Cxyz,NAT)
      USE ISO_C_BINDING
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      COMMON/USELIBRETA/ILIBRETA      
      TYPE(C_PTR),DIMENSION(600)::BASLIB
      COMMON/LIBRETA/BASLIB
      LOGICAL SCHWRZ,SCHSKP,SKIPA,SKIPB,SKIPC,NPSYM                                   
      COMMON/INTFIL/NINTMX  
      COMMON/SHLEXC/NORGSH(3),NORGSP(3),IEXCH,NGTH(4)
      COMMON/RESTAR/NREC,IST,JST,KST,LST           
      COMMON/SHLNOS1/QQ4,IJKL 
!      
      INTEGER :: NINTEGtm,NINTEGt,NSCHWZ,MAXG,IDONTW,IPRINTOPT
      INTEGER :: NPRIMI,NSHELL,NAT
      INTEGER,DIMENSION(NINTEGtm) :: IX2
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DIMENSION BUFP(NINTMX),IX(NINTMX),XINTS(NSH2),GHONDOLIB(MAXG)
      DOUBLE PRECISION,DIMENSION(NINTEGtm) :: BUFP2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AUX
!-----------------------------------------------------------------------
!     2e- Integrals (S,P,D,F,G & L Shells)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CUTOFF = 1.0D-09
      ALLOCATE(AUX(49*900))                                                   
!
      IF(ILIBRETA==0)THEN
       CALL BASCHK(LMAXIMA,KTYPE,NSHELL)
       NANGM =  4                                          
       IF(LMAXIMA==2) NANGM =  6                                          
       IF(LMAXIMA==3) NANGM = 10                                          
       IF(LMAXIMA==4) NANGM = 15                                          
       IF(LMAXIMA==5) NANGM = 21                                          
       IF(LMAXIMA==6) NANGM = 28                                          
       NGTH(4) = 1                                                       
       NGTH(3) = NGTH(4) * NANGM                                         
       NGTH(2) = NGTH(3) * NANGM                                         
       NGTH(1) = NGTH(2) * NANGM                                         
       DO I=1,3                                                       
        NORGSH(I) = 0                                               
        NORGSP(I) = 0                                               
       ENDDO
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NINTEGt  = 0                                                         
      NSCHWZ= 0                                                         
      SCHSKP=.FALSE.                                                    
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO II = IST,NSHELL                            ! II Shell
       J0 = JST                                     ! JJ Shell                                            
       DO JJ = J0,II                                                 
        JST = 1                                                           
        K0 = KST                                                          
        DO KK = K0,JJ                               ! KK Shell                                                 
         KST = 1                                                           
         L0 = LST                                                          
         DO LL = L0,KK                              ! LL Shell
          LST = 1                                                           
          SKIPA =  JJ==KK                                                 
          SKIPB = (II==KK) .or. (JJ==LL)                                
          SKIPC = (II==JJ) .or. (KK==LL)                                
          NPSYM = .FALSE.                                                   
          IF(SKIPA.or.SKIPB.or.SKIPC)GO TO 1                        
          NPSYM = .TRUE.
          IH = II                                                       
          JH = JJ                                                       
          IF(JH<=IH)THEN                                              
           ID = IH                                                     
           JD = JH                                                     
          ELSE                                                           
           ID = JH                                                     
           JD = IH                                                     
          END IF                                                         
          IF(.NOT.SKIPA)                                                &
          SKIPA = (ID==II .and. JD==KK) .or. (ID==JJ .and. JD==LL)               
          IF(.NOT.SKIPB)                                                &
          SKIPB = (ID==II .and. JD==LL) .or. (ID==JJ .and. JD==KK)               
          IF(SKIPA .and. SKIPB)GO TO 2                               
          KH = KK
          IF(KH<=IH) THEN                                              
           ID = IH                                                     
           KD = KH                                                     
          ELSE                                                           
           ID = KH                                                     
           KD = IH                                                     
          END IF                                                         
          IF(.NOT.SKIPC)                                                &
          SKIPC = (ID==II .and. KD==LL) .or. (ID==JJ .and. KD==KK)               
          IF(SKIPA .and. SKIPC) GO TO 3                                
          IF(SKIPB .and. SKIPC) GO TO 4                                
          GO TO 1                                                         
    2     SKIPC = .TRUE.                                                    
          GO TO 1                                                         
    3     SKIPB = .TRUE.                                                    
          GO TO 1                                                         
    4     SKIPA = .TRUE.                                                    
!- - - - - - - - - - - - - - - - - - - - - - - -
!         (II,JJ//KK,LL)                                        
!- - - - - - - - - - - - - - - - - - - - - - - -                            
    1     CONTINUE                                                          
          IEXCH = 1                                                         
          ISH = II                                                          
          JSH = JJ                                                          
          KSH = KK                                                          
          LSH = LL                                                          
          QQ4 = 1                                                          
          IF(SKIPA .and. NPSYM) QQ4 = QQ4+1                                
          IF(SKIPB .and. NPSYM) QQ4 = QQ4+1                                
          GO TO 5                                                         
!- - - - - - - - - - - - - - - - - - - - - - - -
!         (II,KK//JJ,LL)
!- - - - - - - - - - - - - - - - - - - - - - - -                            
    6     IF (SKIPA)GO TO 7                                              
          IEXCH = 2                                                         
          ISH = II                                                          
          JSH = KK                                                          
          KSH = JJ                                                          
          LSH = LL                                                          
          QQ4 = 1                                                          
          IF (SKIPC .and. NPSYM) QQ4 = QQ4+1                               
          GO TO 5                                                         
!- - - - - - - - - - - - - - - - - - - - - - - -
!         (II,LL//JJ,KK)
!- - - - - - - - - - - - - - - - - - - - - - - -                            
    7     IF(SKIPB .or. SKIPC)GO TO 8                                   
          IEXCH = 3                                                         
          ISH = II                                                          
          JSH = LL                                                          
          KSH = JJ                                                          
          LSH = KK                                                          
          QQ4 = 1                                                          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
!         Compute 2e- Integrals                      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
    5     CONTINUE    
          IF(SCHWRZ)THEN                       ! Schwarz Inequality
           IJIJ = (ISH*ISH-ISH)/2 + JSH                                   
           KLKL = (KSH*KSH-KSH)/2 + LSH                                   
           TEST = QQ4*XINTS(IJIJ)*XINTS(KLKL)                             
           SCHSKP = TEST<CUTOFF                                        
           IF(SCHSKP)NSCHWZ = NSCHWZ + 1                                 
          END IF        
          IF(SCHSKP)GO TO 9                                              
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!         Select integral code for ERI calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          if(ILIBRETA==0)then
           LQSUM = KTYPE(ISH) + KTYPE(JSH) + KTYPE(KSH) + KTYPE(LSH) - 4     
           CALL SHELLS(1,ISH,JSH,KSH,LSH,.TRUE.,EX,CS,CP,CD,CF,CG,CH,CI,&
                       NPRIMI,KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,    &
                       NSHELL,Cxyz)                           
           CALL SHELLS(2,ISH,JSH,KSH,LSH,.TRUE.,EX,CS,CP,CD,CF,CG,CH,CI,&
                       NPRIMI,KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,    &
                       NSHELL,Cxyz)                           
           CALL IJPRIM(AUX)                                          
           NORGH = NORGSH(IEXCH)                                           
           CALL ZQOUT(GHONDOLIB,MAXG)                                
           IF(LQSUM==0) THEN                                             
            CALL S0000(GHONDOLIB(1+NORGH),AUX)                        
           ELSE                                                            
            CALL ERISPDFGHIL(GHONDOLIB(1+NORGH),AUX)                       
           END IF
          else if(ILIBRETA==1)then
!lib           CALL ERISVAL(BASLIB,ISH,JSH,KSH,LSH,GHONDOLIB)
          end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!         Write Label & Integral on File 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          if(ILIBRETA==0)then
           CALL QOUT(BUFP,IX,BUFP2,IX2,NINTEGtm,NINTMX,GHONDOLIB,IDONTW,&
                     KLOC,KMIN,KMAX,NSHELL)
          else if(ILIBRETA==1)then
           CALL QOUTlib(BUFP,IX,BUFP2,IX2,NINTEGtm,NINTMX,GHONDOLIB,    &
                        IDONTW,ISH,JSH,KSH,LSH,KLOC,KMIN,KMAX,NSHELL)
          endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    9     CONTINUE                                                          
          GO TO (6,7,8),IEXCH                                         
    8     CONTINUE                                                          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         END DO
        END DO
       END DO
      END DO
      DEALLOCATE(AUX)                                                                                                   
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Final Record on File 1. Calculate NINTEGt.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL FINAL(BUFP,IX,BUFP2,IX2,NINTEGtm,NINTMX,NINTEGt,IDONTW,      &
                 IPRINTOPT)            
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               

! FINAL                                            
      SUBROUTINE FINAL(BUFP,IX,BUFP2,IX2,NINTEGtm,NINTMX,NINTEGt,       &
                       IDONTW,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER :: NINTEGtm,NINTMX,NINTEGt,IDONTW,IPRINTOPT
      DIMENSION BUFP(NINTMX),IX(NINTMX)
      COMMON /RESTAR/ NREC,IST,JST,KST,LST          
      INTEGER,DIMENSION(NINTEGtm) :: IX2
      DOUBLE PRECISION,DIMENSION(NINTEGtm) :: BUFP2
      COMMON /SHLT  / SHLTOL,CUTOFF,ICOUNT
!-----------------------------------------------------------------------                                                                       
      IST = 1                                                           
      JST = 1                                                           
      KST = 1                                                           
      LST = 1                                                           
      NXInteg = ICOUNT-1                                                    
      IF(NXInteg>=0) THEN                                            
       IF(IDONTW==1)THEN
        IJBUFi = (NREC-1)*NINTMX
        do ibuf=1,NXInteg
         IX2  (IJBUFi+ibuf) = IX(ibuf)
         BUFP2(IJBUFi+ibuf) = BUFP(ibuf)
        end do
       ELSE
        NXInteg = -NXInteg 
        WRITE(1)NXInteg,IX,BUFP
       END IF
       NINTEGt = NINTMX*(NREC-1) + ICOUNT-1
      ELSE                                                              
       NINTEGt = NXInteg                                                     
      ENDIF                                                             
      IF(IPRINTOPT==1)WRITE(6,10)NINTEGt,NREC                                 
      RETURN                                                            
!-----------------------------------------------------------------------                                                                       
   10 FORMAT(I20,' 2e- integrals in ',I5,' records')
!-----------------------------------------------------------------------                                                                       
      END                                                               

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
!          11/15/2020 module modified by Juan Felipe Huan Lew Yee      !                                                                     !                                                                      !                                                                      !
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
      Write(6,'(A)')' This is the parallel version of the code'
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
      Write(6,'(A)')' This is the serial version of the code'
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
   10 CALL MPI_BCAST(INTTYPE,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
!HUB      CALL MPI_RECV(IHUB,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,STATUS,IERR)
      IF(INTTYPE==1) THEN
        CALL MPI_RECV(NINTCR,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,STATUS,IERR)
        Call SlvInt(NINTCR)
        GOTO 10
      ELSE IF(INTTYPE==2) THEN
        CALL MPI_RECV(NINTCR,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,STATUS,IERR)
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
      INTEGER NINTCR,N,K,JJ,IERI,ID,NOPT,IER,NOPTCG
      Double Precision ERI
      ALLOCATABLE::ERI(:),IERI(:)
      ALLOCATE(IERI(NINTCR),ERI(NINTCR),STAT=IER)
      IF(IER/=0)CALL ERRORMEM(NINTCR)
      ID=MY_ID
      K=MASTER
      JJ=1
!
   10 CALL MPI_RECV(N,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,STATUS,IERR)
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
   20 CALL MPI_RECV(NOPT,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,STATUS,IERR)
      IF (NOPT==0) THEN
       CALL MPI_RECV(NOPTCG,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,         &
                     STATUS,IERR)
       IF (NOPTCG==0) GOTO 30  
       RETURN 
      ELSE IF (NOPT==1) THEN
       CALL MPI_RECV(N,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,STATUS,IERR)
       CALL SLVHSTJ(N,JJ-1,IERI,ERI)
      ELSE IF (NOPT==2) THEN
       CALL MPI_RECV(N,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,STATUS,IERR)
       CALL SLVHSTK(N,JJ-1,IERI,ERI)
      ELSE IF (NOPT==3) THEN
       CALL MPI_RECV(N,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,STATUS,IERR)
       CALL SLVFORM2JK(N,JJ-1,IERI,ERI)
      ELSE IF (NOPT==4) THEN
       CALL MPI_RECV(N,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,STATUS,IERR)
       CALL SLVFORMJK(N,JJ-1,IERI,ERI)
      ELSE IF (NOPT==5) THEN
       CALL MPI_RECV(N,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,STATUS,IERR)
       CALL SLVHSTARK3(N,JJ-1,IERI,ERI)
      ELSE IF (NOPT==6) THEN
       CALL MPI_RECV(N,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,STATUS,IERR)
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
      CALL MPI_BCAST(NBF,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NBF5,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NBFaux,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      IAUXDIM = NINTCR/(NBF*(NBF+1)/2)
      CALL MPI_RECV(ERIaux(JJ),NINTCR,MPI_REAL8,K,ID,MPI_COMM_WORLD,    &
                    STATUS,IERR)
!
!     Now it has all the integrals in the node, waits for wake-up signal
!     and call the pertinent subroutine (either form J or K integrals)
!
   20 CALL MPI_RECV(NOPT,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,STATUS,IERR)
      IF (NOPT==0) THEN
       CALL MPI_RECV(NOPTCG,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,         &
                     STATUS,IERR)
       IF (NOPTCG==0) GOTO 30
       RETURN 
      ELSE IF (NOPT==1) THEN
       CALL SLVHSTJKRI(NBF,IAUXDIM,ERIaux)
      ELSE IF (NOPT==2) THEN
       CALL SLVQJKMATmRI(NBF,NBF5,IAUXDIM,ERIaux)
      ELSE IF (NOPT==6) THEN
       CALL MPI_RECV(NBF,1,MPI_INTEGER8,K,ID,MPI_COMM_WORLD,STATUS,IERR)
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
      DIMENSION IERI(NN),ERI(NN)
      ALLOCATABLE::P(:),F(:),FF(:)
      ALLOCATE (P(N),F(N),FF(N))
!
!     Get P matrix from Master node
!
      CALL MPI_BCAST(NBF,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
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
      INTEGER,DIMENSION(NN)::IERI
      DOUBLE PRECISION,DIMENSION(NN)::ERI
      ALLOCATABLE::P(:),F(:),FF(:)
      ALLOCATE (P(N),F(N),FF(N))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Get P matrix from Master node
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_BCAST(NBF,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
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
      INTEGER,DIMENSION(NN)::IERI
      DOUBLE PRECISION,DIMENSION(NN)::ERI
      ALLOCATABLE::P(:),F(:),FF(:)
      ALLOCATE (P(N),F(N),FF(N))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Get P matrix from Master node
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_BCAST(NBF,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
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
      INTEGER,DIMENSION(NN)::IERI
      DOUBLE PRECISION,DIMENSION(NN)::ERI
      ALLOCATABLE::P(:),F(:),FF(:)
      ALLOCATE (P(N),F(N),FF(N))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Get P matrix from Master node
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_BCAST(NBF,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
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
      INTEGER,DIMENSION(NN)::IERI
      DOUBLE PRECISION,DIMENSION(NN)::ERI
      ALLOCATABLE :: P1(:),F1(:),FF1(:)
      ALLOCATABLE :: P2(:),F2(:),FF2(:)
      ALLOCATABLE :: P3(:),F3(:),FF3(:)      
      ALLOCATE(P1(N),F1(N),FF1(N),P2(N),F2(N),FF2(N),P3(N),F3(N),FF3(N))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Get P matrix from Master node
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MPI_BCAST(NBF,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
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
      LOGICAL EFIELDL
      INTEGER NBF,NBF5,NSQ
      INTEGER NA,NB,NO1,NSOC
      INTEGER LL,UL,EQPART,UNEQPART
      INTEGER i,k,ik
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::DJ,DK,F,FF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::CJ12,CK12
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::RO
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::H,ADIPX,ADIPY,ADIPZ

      NSQ = NBF*NBF

      CALL MPI_BCAST(NBF5,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)

      ALLOCATE(DJ(NSQ,NBF5), DK(NSQ,NBF5), F(NSQ,NBF5), FF(NSQ,NBF5))
      ALLOCATE(CJ12(NBF5,NBF5), CK12(NBF5,NBF5))
      ALLOCATE(H(NBF,NBF))
      ALLOCATE(RO(NBF5))

      CALL MPI_BCAST(DJ,NSQ*NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DK,NSQ*NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(CJ12,NBF5*NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(CK12,NBF5*NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(RO,NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(H,NBF*NBF,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NA,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NB,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NO1,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NSOC,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(EFIELDL,1,MPI_LOGICAL,MASTER,MPI_COMM_WORLD,IERR)
      IF(EFIELDL) THEN
          ALLOCATE(ADIPX(NBF,NBF), ADIPY(NBF,NBF), ADIPZ(NBF,NBF))
          CALL MPI_BCAST(ADIPX,NBF*NBF,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(ADIPY,NBF*NBF,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(ADIPZ,NBF*NBF,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      END IF

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
!-----------------------------------------------------------------------
!     Including Electric Field ( Note: FEikj is constant )
!-----------------------------------------------------------------------
      IF(EFIELDL)THEN
       !$OMP PARALLEL DO PRIVATE(i,k,ik,J,FEikj)
       do ik=LL,UL
        i = INT((ik-1)/NBF) + 1
        k = MOD(ik-1,NBF) + 1
        DO J=1,NBF5
         FEikj = RO(J)*(EX*ADIPx(i,k)+EY*ADIPy(i,k)+EZ*ADIPz(i,k))
         F(ik,J) = F(ik,J) + FEikj
        ENDDO
       enddo
       !$OMP END PARALLEL DO
       DEALLOCATE(ADIPX, ADIPY, ADIPZ)
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
!                        H O N D O  Calculator                         !
!                                                                      !
!----------------------------------------------------------------------!

! HSandT                                           
      SUBROUTINE HSandT(H,S,TKIN,NBFT,EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,   &
                        KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,   &
                        ZAN,C)        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,PARAMETER :: PI212 = 1.1283791670955D0      
      LOGICAL DOBLE
      DOUBLE PRECISION,DIMENSION(NBFT) :: H,S,TKIN
      INTEGER,DIMENSION(84) :: IX,IY,IZ,JX,JY,JZ
      LOGICAL     LINEAR
      COMMON/ZMAT/LINEAR      
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB                            
      COMMON/HERMIT/H1(55),W1(55)                                        
      COMMON/ROOT/XXROOT,U(13),W(13),NROOTS    ! XX = XXROOT             
      COMMON/SHLNRM/PNRM(84)                                             
      COMMON/STV/XINTT,YINTT,ZINTT,TAA,X0X0,Y0Y0,Z0Z0,                  &       
                 XIXI,YIYI,ZIZI,XJXJ,YJYJ,ZJZJ,NINI,NJNJ 
      LOGICAL                                         IIANDJJ
      COMMON/SYMIND/II,JJ,LIT,LJT,MINI,MINJ,MAXI,MAXJ,IIANDJJ            
      INTEGER,ALLOCATABLE,DIMENSION(:)::IJX,IJY,IJZ
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
      DOUBLE PRECISION,DIMENSION(3,NAT) :: C      
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: Z,ESP1E,SBLK,TBLK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: VBLK,ZBLK,FT
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DIJ,XIN,YIN,ZIN
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: CONI,CONJ
      ALLOCATE(IJX(784),IJY(784),IJZ(784))
      ALLOCATE(Z(NBFT),ESP1E(NBFT))
      ALLOCATE(SBLK(784),TBLK(784),VBLK(784),ZBLK(784),FT(784),DIJ(784))
      ALLOCATE(XIN(343),YIN(343),ZIN(343),CONI(84),CONJ(84))
!-----------------------------------------------------------------------
      DATA IX / 1, 8, 1, 1,15, 1, 1, 8, 8, 1,                           &
               22, 1, 1,15,15, 8, 1, 8, 1, 8,                           &
               29, 1, 1,22,22, 8, 1, 8, 1,15,                           &
               15, 1,15, 8, 8,                                          &
               36, 1, 1,29,29, 8, 1, 8, 1,22,                           &
               22,15, 1,15, 1,22, 8, 8,15,15,                           &
                8,                                                      &
               43, 1, 1,36,36, 8, 1, 8, 1,29,                           &
               29,15, 1,15, 1,29, 8, 8,22,22,                           &
                1,22,22,15, 8,15, 8,15/                                  
      DATA IY / 1, 1, 8, 1, 1,15, 1, 8, 1, 8,                           &
                1,22, 1, 8, 1,15,15, 1, 8, 8,                           &
                1,29, 1, 8, 1,22,22, 1, 8,15,                           &
                1,15, 8,15, 8,                                          &
                1,36, 1, 8, 1,29,29, 1, 8,15,                           &
                1,22,22, 1,15, 8,22, 8,15, 8,                           &
               15,                                                      &
                1,43, 1, 8, 1,36,36, 1, 8,15,                           &
                1,29,29, 1,15, 8,29, 8,22, 1,                           &
               22,15, 8,22,22, 8,15,15/                                  
      DATA IZ / 1, 1, 1, 8, 1, 1,15, 1, 8, 8,                           &
                1, 1,22, 1, 8, 1, 8,15,15, 8,                           &
                1, 1,29, 1, 8, 1, 8,22,22, 1,                           &
               15,15, 8, 8,15,                                          &
                1, 1,36, 1, 8, 1, 8,29,29, 1,                           &
               15, 1,15,22,22, 8, 8,22, 8,15,                           &
               15,                                                      &
                1, 1,43, 1, 8, 1, 8,36,36, 1,                           &
               15, 1,15,29,29, 8, 8,29, 1,22,                           &
               22, 8,15, 8,15,22,22,15/                                  
      DATA JX / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,                           &
                3, 0, 0, 2, 2, 1, 0, 1, 0, 1,                           &
                4, 0, 0, 3, 3, 1, 0, 1, 0, 2,                           &
                2, 0, 2, 1, 1,                                          &
                5, 0, 0, 4, 4, 1, 0, 1, 0, 3,                           &
                3, 2, 0, 2, 0, 3, 1, 1, 2, 2,                           &
                1,                                                      &
                6, 0, 0, 5, 5, 1, 0, 1, 0, 4,                           &
                4, 2, 0, 2, 0, 4, 1, 1, 3, 3,                           &
                0, 3, 3, 2, 1, 2, 1, 2/                                  
      DATA JY / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,                           &
                0, 3, 0, 1, 0, 2, 2, 0, 1, 1,                           &
                0, 4, 0, 1, 0, 3, 3, 0, 1, 2,                           &
                0, 2, 1, 2, 1,                                          &
                0, 5, 0, 1, 0, 4, 4, 0, 1, 2,                           &
                0, 3, 3, 0, 2, 1, 3, 1, 2, 1,                           &
                2,                                                      &
                0, 6, 0, 1, 0, 5, 5, 0, 1, 2,                           &
                0, 4, 4, 0, 2, 1, 4, 1, 3, 0,                           &
                3, 2, 1, 3, 3, 1, 2, 2/                                  
      DATA JZ / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,                           &
                0, 0, 3, 0, 1, 0, 1, 2, 2, 1,                           &
                0, 0, 4, 0, 1, 0, 1, 3, 3, 0,                           &
                2, 2, 1, 1, 2,                                          &
                0, 0, 5, 0, 1, 0, 1, 4, 4, 0,                           &
                2, 0, 2, 3, 3, 1, 1, 3, 1, 2,                           &
                2,                                                      &
                0, 0, 6, 0, 1, 0, 1, 5, 5, 0,                           &
                2, 0, 2, 4, 4, 1, 1, 4, 0, 3,                           &
                3, 1, 2, 1, 2, 3, 3, 2/                                 
!-----------------------------------------------------------------------
!                      H, S & TKIN integrals
!-----------------------------------------------------------------------
      TOL = 20*2.30258D0     
      ZBLK = 0.0D+00
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     I SHELL                                               
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO III = 1,NSHELL                                           
       I = KATOM(III)                                                  
       XIXI = C(1,I)                                                    
       YIYI = C(2,I)                                                    
       ZIZI = C(3,I)                                                    
       I1 = KSTART(III)                                                
       I2 = I1+KNG(III)-1                                              
       LIT = KTYPE(III)                                                
       MINI = KMIN(III)                                                
       MAXI = KMAX(III)                                                
       LOCI = KLOC(III)-MINI                                     
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      J SHELL                                               
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO JJ = 1,III                                          
        J = KATOM(JJ)                                               
        XJXJ = C(1,J)                                                 
        YJYJ = C(2,J)                                                 
        ZJZJ = C(3,J)                                                 
        J1 = KSTART(JJ)                                             
        J2 = J1+KNG(JJ)-1                                           
        LJT = KTYPE(JJ)                                             
        MINJ = KMIN(JJ)                                             
        MAXJ = KMAX(JJ)                                             
        LOCJ = KLOC(JJ)-MINJ                                  
        NROOTS = (LIT+LJT-2)/2+1                                    
        RR = (XIXI-XJXJ)**2+(YIYI-YJYJ)**2+(ZIZI-ZJZJ)**2                       
        IIANDJJ = III == JJ                                          
        IJ = 0                                                      
        MAX = MAXJ                                                  
        DO I = MINI,MAXI                                        
         NXNX = IX(I)                                               
         NYNY = IY(I)                                               
         NZNZ = IZ(I)                                               
         IF (IIANDJJ) MAX = I                                       
         DO J = MINJ,MAX                                      
          IJ = IJ+1                                             
          IJX(IJ) = NXNX+JX(J)                                    
          IJY(IJ) = NYNY+JY(J)                                    
          IJZ(IJ) = NZNZ+JZ(J)                                    
          FT (IJ) = 2*(JX(J)+JY(J)+JZ(J)) + 3                   
         END DO
        END DO
        CALL VCLR( SBLK,1,IJ)                                       
        CALL VCLR( TBLK,1,IJ)                                       
        CALL VCLR( VBLK,1,IJ)                                       
        JGMAX = J2                                                  
        DO IG = I1,I2                                           
         AI = EX(IG)                                              
         ARRI = AI*RR                                             
         AXI = AI*XIXI                                              
         AYI = AI*YIYI                                              
         AZI = AI*ZIZI                                              
         CALL SETCONI(CONI,IG,CS,CP,CD,CF,CG,CH,CI,NPRIMI)                                    
         IF (IIANDJJ) JGMAX = IG                                    
         DO JG = J1,JGMAX                                     
          AJ = EX(JG)                                           
          AA = AI+AJ                                            
          AA1 = 1.0d0/AA                                          
          DUM = AJ*ARRI*AA1                                     
          IF(DUM<=TOL)THEN
           FAC = EXP(-DUM)                                       
           CALL SETCONI(CONJ,JG,CS,CP,CD,CF,CG,CH,CI,NPRIMI)                                 
           AX = (AXI+AJ*XJXJ)*AA1                                  
           AY = (AYI+AJ*YJYJ)*AA1                                  
           AZ = (AZI+AJ*ZJZJ)*AA1                                  
           DOBLE = IIANDJJ.and.IG/=JG                             
           MAX = MAXJ                                            
           NN = 0                                                
           DTWO = 1.0d0                                              
           IF(DOBLE)DTWO = 2.0D0                                   
           SPDIJ = CS(IG)*CP(JG)*FAC                               
           DO I = MINI,MAXI                                  
            IF (IIANDJJ) MAX = I                                 
            FACI=FAC*CONI(I)*PNRM(I)*DTWO                      
            NN1=NN+1                                           
            DO J = MINJ,MAX                                
             NN = NN+1                                       
             DIJ(NN)=FACI*CONJ(J)*PNRM(J)                    
            END DO
            IF(MINJ<=1.and.I>1.and.DOBLE)DIJ(NN1)=DIJ(NN1)*0.5D0+SPDIJ
           END DO
            TAA = SQRT(AA1)                                       
            T1 = -2.0D0*AJ*AJ*TAA                                   
            T2 = -0.5D0*TAA                                         
            X0X0 = AX                                               
            Y0Y0 = AY                                               
            Z0Z0 = AZ                                               
            IN = -7                                               
            DO I = 1,LIT                                      
             IN = IN+7                                          
             NINI = I                                             
             DO J = 1,LJT                                   
              JN = IN+J                                       
              NJNJ = J                                          
              CALL STVINT(H1,W1)                                     
              XIN(JN) = XINTT*TAA                              
              YIN(JN) = YINTT*TAA                              
              ZIN(JN) = ZINTT*TAA                              
              NJNJ = J+2                                        
              CALL STVINT(H1,W1)                                     
              XIN(JN+49) = XINTT*T1                            
              YIN(JN+49) = YINTT*T1                            
              ZIN(JN+49) = ZINTT*T1                            
              NJNJ = J-2                                        
              IF (NJNJ > 0) THEN                             
               CALL STVINT(H1,W1)                                  
              ELSE                                            
               XINTT = 0.0D+00                                  
               YINTT = 0.0D+00                                  
               ZINTT = 0.0D+00                                  
              END IF                                          
              N = (J-1)*(J-2)                                 
              DUM = N * T2                                    
              XIN(JN+98) = XINTT*DUM                           
              YIN(JN+98) = YINTT*DUM                           
              ZIN(JN+98) = ZINTT*DUM                           
              IF(LINEAR)THEN                                  
               NJNJ = J+1                                     
               CALL STVINT(H1,W1)                                  
               XIN(JN+147) = XINTT*TAA                       
               YIN(JN+147) = YINTT*TAA                       
               NJNJ = J-1                                     
               IF(NJNJ>0)THEN                          
                CALL STVINT(H1,W1)                               
               ELSE                                         
                XINTT = 0.0D+00                               
                YINTT = 0.0D+00                               
               END IF                                       
               XIN(JN+196) = XINTT*TAA*NJNJ                    
               YIN(JN+196) = YINTT*TAA*NJNJ                    
              END IF                                          
             END DO
            END DO
            DO I = 1,IJ                                       
             NXNX = IJX(I)                                        
             NYNY = IJY(I)                                        
             NZNZ = IJZ(I)                                        
             DUM   = XIN(NXNX)*YIN(NYNY)*ZIN(NZNZ)    
             DUM1X = (XIN(NXNX+49)+XIN(NXNX+98))*YIN(NYNY)*ZIN(NZNZ)    
             DUM1Y = (YIN(NYNY+49)+YIN(NYNY+98))*XIN(NXNX)*ZIN(NZNZ)    
             DUM1Z = (ZIN(NZNZ+49)+ZIN(NZNZ+98))*XIN(NXNX)*YIN(NYNY)    
             DUM1  = DUM1X + DUM1Y + DUM1Z                      
              SBLK(I) =  SBLK(I) + DIJ(I)* DUM                  
              TBLK(I) =  TBLK(I) + DIJ(I)*(DUM*AJ*FT(I)+DUM1) 
             IF(LINEAR)THEN                                     
              DUM2 = XIN(NXNX+147)*YIN(NYNY+196)                        &
                   - XIN(NXNX+196)*YIN(NYNY+147)   
              ZBLK(I) = ZBLK(I) + DIJ(I)*DUM2*ZIN(NZNZ)         
             END IF                                             
            END DO
            DUM = PI212*AA1                                    
            DO I = 1,IJ                                    
             DIJ(I) = DIJ(I)*DUM                             
            END DO
           AAX = AA*AX                                           
           AAY = AA*AY                                           
           AAZ = AA*AZ                                           
           DO IC = 1,NAT 
            ZNUC = -ZAN(IC)                                 
            CX = C(1,IC)                                    
            CY = C(2,IC)                                    
            CZ = C(3,IC)                                    
            XXROOT = AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)         
            IF (NROOTS<=3) CALL RT123                        
            IF (NROOTS==4) CALL ROOT4                        
            IF (NROOTS==5) CALL ROOT5                        
            IF (NROOTS>=6) CALL ROOT6                        
            MM = 0                                             
            DO K = 1,NROOTS                                
             UU = AA*U(K)                                    
             WW = W(K)*ZNUC                                  
             TT = 1.0d0/(AA+UU)                                
             TAA = SQRT(TT)                                  
             X0X0 = (AAX+UU*CX)*TT                             
             Y0Y0 = (AAY+UU*CY)*TT                             
             Z0Z0 = (AAZ+UU*CZ)*TT                             
             IN = -7+MM                                      
             J0 = 2                                          
             XIN(IN+8) = W1(1)                                  
             YIN(IN+8) = W1(1)
             ZIN(IN+8) = W1(1)*WW                               
             DO I = 1,LIT                                
              IN = IN+7                                    
              NINI = I                                       
              DO J = J0,LJT                            
               JN = IN+J                                 
               NJNJ = J                                    
               CALL STVINT(H1,W1)                               
               XIN(JN) = XINTT                            
               YIN(JN) = YINTT                            
               ZIN(JN) = ZINTT*WW                         
              END DO
              J0 = 1                                       
             END DO
             MM = MM+49                                      
            END DO
            DO I = 1,IJ                                    
              NXNX = IJX(I)                                     
              NYNY = IJY(I)                                     
              NZNZ = IJZ(I)                                     
              DUM = 0.0D+00                                      
              MM = 0                                          
              DO K = 1,NROOTS                             
               DUM = DUM+XIN(NXNX+MM)*YIN(NYNY+MM)*ZIN(NZNZ+MM)   
               MM = MM+49                                   
              END DO
              VBLK(I) = VBLK(I) + DUM*DIJ(I)               
             END DO
           END DO
          END IF                                             
         END DO
        END DO
        MAX = MAXJ                                                  
        NN = 0                                                      
        DO I = MINI,MAXI                                        
         LI = LOCI+I                                              
         IN = (LI*(LI-1))/2                                       
         IF (IIANDJJ) MAX = I                                       
         DO J = MINJ,MAX                                      
          LJ = LOCJ+J                                           
          JN = LJ+IN                                            
          NN = NN+1                                             
          H(JN) =  TBLK(NN) + VBLK(NN)                         
          S(JN) =  SBLK(NN)                                    
          TKIN(JN) =  TBLK(NN)                                    
          IF(LINEAR) Z(JN) = ZBLK(NN)                            
         END DO
        END DO
       END DO
      END DO
!-----------------------------------------------------------------------
      DEALLOCATE(IJX,IJY,IJZ)
      DEALLOCATE(Z,ESP1E,SBLK,TBLK,VBLK,ZBLK,FT)
      DEALLOCATE(DIJ,XIN,YIN,ZIN,CONI,CONJ)              
      RETURN                                                            
      END                                                               
      
! PRCALC                                           
      SUBROUTINE PRCALC(XVAL,WINT,NVAL,L2,EX,CS,CP,CD,CF,CG,CH,CI,      &
                 NPRIMI,KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,   &
                 Cxyz)                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      LOGICAL                                         IIANDJJ
      COMMON/SYMIND/II,JJ,LIT,LJT,MINI,MINJ,MAXI,MAXJ,IIANDJJ  
      COMMON/XYZORB/TXYZ,X00,Y00,Z00,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB                     
      LOGICAL IIIandJJJ,NORMA,DOUBLE
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DIMENSION XVAL(NVAL*L2),WINT(*)                                   
      DIMENSION DIJ(784),IJX(784),IJY(784),IJZ(784)      
      DIMENSION IX(84),IY(84),IZ(84),JX(84),JY(84),JZ(84)
      DATA IX / 1, 8, 1, 1,15, 1, 1, 8, 8, 1,                           &
               22, 1, 1,15,15, 8, 1, 8, 1, 8,                           &
               29, 1, 1,22,22, 8, 1, 8, 1,15,                           &
               15, 1,15, 8, 8,                                          &
               36, 1, 1,29,29, 8, 1, 8, 1,22,                           &
               22,15, 1,15, 1,22, 8, 8,15,15,                           &
                8,                                                      &
               43, 1, 1,36,36, 8, 1, 8, 1,29,                           &
               29,15, 1,15, 1,29, 8, 8,22,22,                           &
                1,22,22,15, 8,15, 8,15/                                  
      DATA IY / 1, 1, 8, 1, 1,15, 1, 8, 1, 8,                           &
                1,22, 1, 8, 1,15,15, 1, 8, 8,                           &
                1,29, 1, 8, 1,22,22, 1, 8,15,                           &
                1,15, 8,15, 8,                                          &
                1,36, 1, 8, 1,29,29, 1, 8,15,                           &
                1,22,22, 1,15, 8,22, 8,15, 8,                           &
               15,                                                      &
                1,43, 1, 8, 1,36,36, 1, 8,15,                           &
                1,29,29, 1,15, 8,29, 8,22, 1,                           &
               22,15, 8,22,22, 8,15,15/                                  
      DATA IZ / 1, 1, 1, 8, 1, 1,15, 1, 8, 8,                           &
                1, 1,22, 1, 8, 1, 8,15,15, 8,                           &
                1, 1,29, 1, 8, 1, 8,22,22, 1,                           &
               15,15, 8, 8,15,                                          &
                1, 1,36, 1, 8, 1, 8,29,29, 1,                           &
               15, 1,15,22,22, 8, 8,22, 8,15,                           &
               15,                                                      &
                1, 1,43, 1, 8, 1, 8,36,36, 1,                           &
               15, 1,15,29,29, 8, 8,29, 1,22,                           &
               22, 8,15, 8,15,22,22,15/                                  
      DATA JX / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,                           &
                3, 0, 0, 2, 2, 1, 0, 1, 0, 1,                           &
                4, 0, 0, 3, 3, 1, 0, 1, 0, 2,                           &
                2, 0, 2, 1, 1,                                          &
                5, 0, 0, 4, 4, 1, 0, 1, 0, 3,                           &
                3, 2, 0, 2, 0, 3, 1, 1, 2, 2,                           &
                1,                                                      &
                6, 0, 0, 5, 5, 1, 0, 1, 0, 4,                           &
                4, 2, 0, 2, 0, 4, 1, 1, 3, 3,                           &
                0, 3, 3, 2, 1, 2, 1, 2/                                  
      DATA JY / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,                           &
                0, 3, 0, 1, 0, 2, 2, 0, 1, 1,                           &
                0, 4, 0, 1, 0, 3, 3, 0, 1, 2,                           &
                0, 2, 1, 2, 1,                                          &
                0, 5, 0, 1, 0, 4, 4, 0, 1, 2,                           &
                0, 3, 3, 0, 2, 1, 3, 1, 2, 1,                           &
                2,                                                      &
                0, 6, 0, 1, 0, 5, 5, 0, 1, 2,                           &
                0, 4, 4, 0, 2, 1, 4, 1, 3, 0,                           &
                3, 2, 1, 3, 3, 1, 2, 2/                                  
      DATA JZ / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,                           &
                0, 0, 3, 0, 1, 0, 1, 2, 2, 1,                           &
                0, 0, 4, 0, 1, 0, 1, 3, 3, 0,                           &
                2, 2, 1, 1, 2,                                          &
                0, 0, 5, 0, 1, 0, 1, 4, 4, 0,                           &
                2, 0, 2, 3, 3, 1, 1, 3, 1, 2,                           &
                2,                                                      &
                0, 0, 6, 0, 1, 0, 1, 5, 5, 0,                           &
                2, 0, 2, 4, 4, 1, 1, 4, 0, 3,                           &
                3, 1, 2, 1, 2, 3, 3, 2/
!
      PARAMETER (SQRT3 = 1.73205080756887729353D+00)
      PARAMETER (SQRT5 = 2.23606797749978969641D+00)
      PARAMETER (SQRT7 = 2.64575131106459059050D+00) 
      PARAMETER (SQRT11= 3.31662479035539984911D+00)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NIJ = 784*NVAL                                                    
      TOL = 20*2.30258D0
      NORMA = .TRUE.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO III=1,NSHELL                                               
       I    = KATOM(III)                                                  
       XI   = Cxyz(1,I)                                                     
       YI   = Cxyz(2,I)                                                     
       ZI   = Cxyz(3,I)                                                     
       I1   = KSTART(III)                                                 
       I2   = I1 + KNG(III) - 1                                           
       LIT  = KTYPE(III)                                                  
       MINI = KMIN(III)                                                   
       MAXI = KMAX(III)                                                   
       LOCI = KLOC(III) - MINI                                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                      
       DO JJ=1,III                                                   
        J    = KATOM(JJ)                                                  
        XJ   = Cxyz(1,J)                                                     
        YJ   = Cxyz(2,J)                                                     
        ZJ   = Cxyz(3,J)                                                     
        J1   = KSTART(JJ)                                                 
        J2   = J1 + KNG(JJ) - 1                                           
        LJT  = KTYPE(JJ)                                                  
        MINJ = KMIN(JJ)                                                   
        MAXJ = KMAX(JJ)                                                   
        LOCJ = KLOC(JJ) - MINJ                                            
        RR = (XI-XJ)**2 + (YI-YJ)**2 + (ZI-ZJ)**2                     
        IIIandJJJ = III==JJ                                                 
        CALL VCLR(WINT,1,NIJ)                                             
        IJ = 0                                                            
        MAX = MAXJ                                                        
        DO I=MINI,MAXI                                               
         NIX = IX(I)                                                     
         NIY = IY(I)                                                     
         NIZ = IZ(I)                                                     
         IF (IIIandJJJ) MAX = I                                             
         DO J=MINJ,MAX                                             
          IJ = IJ+1                                                   
          IJX(IJ) = NIX+JX(J)                                          
          IJY(IJ) = NIY+JY(J)                                          
          IJZ(IJ) = NIZ+JZ(J)                                          
         END DO
        END DO
!
        JGMAX = J2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                        
        DO IG=I1,I2                                                   
         AI  = EX(IG)                                                      
         CSI = CS(IG)                                                      
         CPI = CP(IG)                                                      
         CDI = CD(IG)                                                      
         CFI = CF(IG)                                                      
         CGI = CG(IG)                                                      
         CHI = CH(IG)                                                      
         CII = CI(IG)                                                      
         IF(IIIandJJJ)JGMAX = IG                                             
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         DO JG=J1,JGMAX                                                
          AJ  = EX(JG)                                                      
          CSJ = CS(JG)                                                      
          CPJ = CP(JG)                                                      
          CDJ = CD(JG)                                                      
          CFJ = CF(JG)                                                      
          CGJ = CG(JG)                                                      
          CHJ = CH(JG)                                                      
          CIJ = CI(JG)                                                      
          AA  = AI + AJ                                                     
          AA1 = 1.0d0/AA                                                      
          AX  = (AI*XI + AJ*XJ)*AA1                                         
          AY  = (AI*YI + AJ*YJ)*AA1                                         
          AZ  = (AI*ZI + AJ*ZJ)*AA1                                         
          DUM = AI*AJ*RR*AA1                                                
          IF(DUM<=TOL)THEN
           FAC = EXP(-DUM)                                                   
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!          Density Factors
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           DOUBLE = IIIandJJJ.and.IG/=JG                                       
           MAX = MAXJ                                                        
           NN  = 0                                                           
           DUM1 = 0.0D+00                                                       
           DUM2 = 0.0D+00                                                       
           DO I = MINI,MAXI                                              
            IF( I==1)            DUM1 = CSI*FAC                          
            IF( I==2)            DUM1 = CPI*FAC                          
            IF( I==5)            DUM1 = CDI*FAC                          
            IF( I==8  .and.NORMA)DUM1 = DUM1*SQRT3                       
            IF( I==11)           DUM1 = CFI*FAC                          
            IF((I==14).and.NORMA)DUM1 = DUM1*SQRT5                       
            IF((I==20).and.NORMA)DUM1 = DUM1*SQRT3                       
            IF( I==21)           DUM1 = CGI*FAC                          
            IF((I==24).and.NORMA)DUM1 = DUM1*SQRT7                       
            IF((I==30).and.NORMA)DUM1 = DUM1*SQRT5/SQRT3                 
            IF((I==33).and.NORMA)DUM1 = DUM1*SQRT3                       
            IF( I==36)           DUM1 = CHI*FAC                          
            IF((I==39).and.NORMA)DUM1 = DUM1*3.0d0                       
            IF((I==45).and.NORMA)DUM1 = DUM1*SQRT7/SQRT3                 
            IF((I==51).and.NORMA)DUM1 = DUM1*SQRT3                       
            IF((I==54).and.NORMA)DUM1 = DUM1*SQRT5/SQRT3                 
            IF( I==57)           DUM1 = CII*FAC                          
            IF((I==60).and.NORMA)DUM1 = DUM1*SQRT11                      
            IF((I==66).and.NORMA)DUM1 = DUM1*3.0d0/SQRT3                 
            IF((I==72).and.NORMA)DUM1 = DUM1*SQRT3                       
            IF((I==75).and.NORMA)DUM1 = DUM1*SQRT7/(SQRT3*SQRT5)         
            IF((I==78).and.NORMA)DUM1 = DUM1*SQRT5                       
            IF((I==84).and.NORMA)DUM1 = DUM1*SQRT5/SQRT3                 
            IF(IIIandJJJ)MAX = I                                              
!---------------------------------------------------------------
            DO J = MINJ,MAX                                            
             NN = NN+1                                                   
             IF(J==1) THEN                                             
               DUM2 = DUM1*CSJ                                           
               IF(DOUBLE .and. I==1) DUM2 = DUM2 + DUM2                
               IF(DOUBLE .and. I>1) DUM2 = DUM2 + CSI*CPJ*FAC         
             ELSE IF( J==2) THEN                                       
               DUM2 = DUM1*CPJ                                           
               IF(DOUBLE) DUM2 = DUM2 + DUM2                             
             ELSE IF( J==5) THEN                                       
               DUM2 = DUM1*CDJ                                           
               IF(DOUBLE) DUM2 = DUM2 + DUM2                             
             ELSE IF((J==8).and.NORMA) THEN                             
               DUM2 = DUM2*SQRT3                                         
             ELSE IF (J==11) THEN                                      
               DUM2 = DUM1*CFJ                                           
               IF (DOUBLE) DUM2 = DUM2+DUM2                              
             ELSE IF((J==14).and.NORMA) THEN                            
               DUM2 = DUM2*SQRT5                                         
             ELSE IF((J==20).and.NORMA) THEN                            
               DUM2 = DUM2*SQRT3                                         
             ELSE IF(J==21) THEN                                       
               DUM2 = DUM1*CGJ                                           
               IF (DOUBLE) DUM2 = DUM2+DUM2                              
             ELSE IF((J==24).and.NORMA) THEN                            
               DUM2 = DUM2*SQRT7                                         
             ELSE IF((J==30).and.NORMA) THEN                            
               DUM2 = DUM2*SQRT5/SQRT3                                   
             ELSE IF((J==33).and.NORMA) THEN                            
               DUM2 = DUM2*SQRT3                                         
             ELSE IF( J==36) THEN                                      
               DUM2 = DUM1*CHJ                                           
               IF (DOUBLE) DUM2 = DUM2+DUM2                              
             ELSE IF((J==39).and.NORMA) THEN                            
               DUM2 = DUM2*3.0d0                                         
             ELSE IF((J==45).and.NORMA) THEN                            
               DUM2 = DUM2*SQRT7/SQRT3                                   
             ELSE IF((J==51).and.NORMA) THEN                            
               DUM2 = DUM2*SQRT3                                         
             ELSE IF((J==54).and.NORMA) THEN                            
               DUM2 = DUM2*SQRT5/SQRT3                                   
             ELSE IF( J==57) THEN                                      
               DUM2 = DUM1*CIJ                                           
               IF (DOUBLE) DUM2 = DUM2+DUM2                              
             ELSE IF((J==60).and.NORMA) THEN                            
               DUM2 = DUM2*SQRT11                                        
             ELSE IF((J==66).and.NORMA) THEN                            
               DUM2 = DUM2*3.0d0/SQRT3                                   
             ELSE IF((J==72).and.NORMA) THEN                            
               DUM2 = DUM2*SQRT3                                         
             ELSE IF((J==75).and.NORMA) THEN                            
               DUM2 = DUM2*SQRT7/(SQRT3*SQRT5)                           
             ELSE IF((J==78).and.NORMA) THEN                            
               DUM2 = DUM2*SQRT5                                         
             ELSE IF((J==84).and.NORMA) THEN                            
               DUM2 = DUM2*SQRT5/SQRT3                                   
             END IF                                                      
             DIJ(NN) = DUM2                                              
            END DO
!---------------------------------------------------------------
           END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
!          Integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
           CALL INTMOM(LIT,LJT,IJ,IJX,IJY,IJZ,DIJ,WINT,AA,AX,AY,AZ)        
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
          END IF
         END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        MAX = MAXJ                                                        
        DO K=1,NVAL                                                  
         NL2 = (K-1)*L2                                                 
         NN  = (K-1)*IJ                                                 
         DO I=MINI,MAXI                                            
          LI = LOCI + I                                               
          IN = LI*(LI-1)/2 + NL2                                      
          IF (IIIandJJJ) MAX = I                                          
          DO J=MINJ,MAX                                          
           LJ = LOCJ + J                                            
           JN = LJ + IN                                             
           NN = NN+1                                                
           XVAL(JN) = WINT(NN)                                      
          END DO
         END DO
        END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END DO
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               

! INIINTQUAD
      SUBROUTINE INIINTQUAD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DOUBLE PRECISION,PARAMETER :: SQRT3 = 1.73205080756887729353D0
      DOUBLE PRECISION,PARAMETER :: SQRT5 = 2.23606797749978969641D0
      DOUBLE PRECISION,PARAMETER :: SQRT7 = 2.64575131106459059050D0
      DOUBLE PRECISION,PARAMETER :: SQRT9 = 3.00000000000000000000D0
      DOUBLE PRECISION,PARAMETER :: SQRT11 = 3.31662479035539984911D0
      DOUBLE PRECISION,PARAMETER :: PI = 3.141592653589793238D0      
      COMMON/HERMIT/H1(55),W1(55)
      COMMON/RYSPAR/XASYMP(13),RTSASY(13,13),WTSASY(13,13),             &
                    NAUXS(13),MAPRYS(13),RTSAUX(55,8),WTSAUX(55,8)
      COMMON/SHLNRM/PNRM(84)     
      DOUBLE PRECISION,DIMENSION(55)::RTS,WTS,WRK
      DOUBLE PRECISION,DIMENSION(0:54)::ALPHA1,BETA1
!-----------------------------------------------------------------------
!     Set up the primitive factors for the 1e- integrals.
!     Initialize for the integral quadratures.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SQRT53 = SQRT5/SQRT3
      SQRT73 = SQRT7/SQRT3
      SQRT75 = SQRT7/SQRT5
      SQRT753 = SQRT75/SQRT3
      do i=1,84
       if(i==1.or.i==2.or.i==5.or.i==11.or.i==21.or.i==36.or.i==57)then
        PNRMi = 1.0D0
       elseif(i==8.or.i==20.or.i==33)then
        PNRMi = PNRMi*SQRT3
       elseif(i==14)then
        PNRMi = PNRMi*SQRT5
       elseif(i==24)then
        PNRMi = PNRMi*SQRT7
       elseif(i==30)then
        PNRMi = PNRMi*SQRT53
       elseif(i==39)then
        PNRMi = PNRMi*SQRT9
       elseif(i==45)then
        PNRMi = PNRMi*SQRT73
       elseif(i==51)then
        PNRMi = PNRMi*SQRT3
       elseif(i==54)then
        PNRMi = PNRMi*SQRT53
       elseif(i==60)then
        PNRMi = PNRMi*SQRT11
       elseif(i==66)then
        PNRMi = PNRMi*SQRT3
       elseif(i==72)then
        PNRMi = PNRMi*SQRT3
       elseif(i==75)then
        PNRMi = PNRMi*SQRT753
       elseif(i==78)then
        PNRMi = PNRMi*SQRT5
       elseif(i==84)then
        PNRMi = PNRMi*SQRT53
       endif
       PNRM(i) = PNRMi
      enddo
!-----------------------------------------------------------------------
!     Initialize the RYS quadrature procedure
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      QUART = 0.25D0
      EPS   = 1.0D-14
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Roots of the Hermite Polynomials of order 2n
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      XASYMP( 1)=29.0D0
      XASYMP( 2)=37.0D0
      XASYMP( 3)=43.0D0
      XASYMP( 4)=49.0D0
      XASYMP( 5)=55.0D0
      XASYMP( 6)=60.0D0
      XASYMP( 7)=65.0D0
      XASYMP( 8)=71.0D0
      XASYMP( 9)=76.0D0
      XASYMP(10)=81.0D0
      XASYMP(11)=86.0D0
      XASYMP(12)=91.0D0
      XASYMP(13)=96.0D0
!
      DO I=1,13
       N=2*I
       DO J=0,N-1
        ALPHA1(J) = 0.0D0
       ENDDO
       BETA1(0)=SQRT(PI)
       DO J=1,N-1
        BETA1(J) = J/2.0D0
       END DO
!      QL Procedure
       CALL RYSGW(N,ALPHA1,BETA1,EPS,RTS,WTS,IERR,WRK)
       DO J=1,I
        RTSASY(J,I) = RTS(I+J)*RTS(I+J)
        WTSASY(J,I) = WTS(I+J)
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Auxiliary Grids
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NAUXS( 1)=20
      NAUXS( 2)=25
      NAUXS( 3)=30
      NAUXS( 4)=30
      NAUXS( 5)=35
      NAUXS( 6)=40
      NAUXS( 7)=40
      NAUXS( 8)=40
      NAUXS( 9)=45
      NAUXS(10)=50
      NAUXS(11)=50
      NAUXS(12)=55
      NAUXS(13)=55
!
      MAPRYS( 1)=1
      MAPRYS( 2)=2
      MAPRYS( 3)=3
      MAPRYS( 4)=3
      MAPRYS( 5)=4
      MAPRYS( 6)=5
      MAPRYS( 7)=5
      MAPRYS( 8)=5
      MAPRYS( 9)=6
      MAPRYS(10)=7
      MAPRYS(11)=7
      MAPRYS(12)=8
      MAPRYS(13)=8
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Auxiliary Quadrature = Shifted Legendre
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NAUXSV=0
      IGRID=0
      DO 1 M=1,13
       NAUX=NAUXS(M)
       IF(NAUX==NAUXSV)GOTO 1
       IGRID=IGRID+1
       NAUXSV=NAUX
       DO I=0,NAUX-1
        ALPHA1(I) = 0.5D0
       ENDDO
       BETA1(0)= 1.0D0
       DO I=1,NAUX-1
        BETA1(I) = QUART/(4.0D0-(1.0D0/(I*I)))
       END DO
!      QL Procedure
       CALL RYSGW(NAUX,ALPHA1,BETA1,EPS,RTS,WTS,IERR,WRK)
       DO I=1,NAUX
        RTSAUX(I,IGRID) = RTS(I)
        WTSAUX(I,IGRID) = WTS(I)
       ENDDO
    1 CONTINUE
!-----------------------------------------------------------------------
!     Initialize Roots and Weights for Gauss-Hermite quadrature
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL SETHERMITE(H1,W1)
!-----------------------------------------------------------------------
      RETURN
      END
      
! SETHERMITE
      SUBROUTINE SETHERMITE(HTOTAL,WTOTAL)                                             
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION,DIMENSION(55),INTENT(OUT)::HTOTAL,WTOTAL
!                                                                       
      HTOTAL(1)=  0.0D0
!                                                                       
      HTOTAL(2)= -0.70710678118654752440D0                                  
      HTOTAL(3)=  0.70710678118654752440D0                                  
!                                                                       
      HTOTAL(4)= -1.22474487139158904910D0                                  
      HTOTAL(5)=  0.0D0                                                        
      HTOTAL(6)=  1.22474487139158904910D0                                  
!                                                                       
      HTOTAL(7)= -1.65068012388578455588D0                                  
      HTOTAL(8)= -0.52464762327529031788D0                                  
      HTOTAL(9)=  0.52464762327529031788D0                                  
      HTOTAL(10)=  1.65068012388578455588D0                                  
!                                                                       
      HTOTAL(11)= -2.02018287045608563293D0                                  
      HTOTAL(12)= -0.95857246461381850711D0                                  
      HTOTAL(13)=  0.0D0                                                        
      HTOTAL(14)=  0.95857246461381850711D0                                  
      HTOTAL(15)=  2.02018287045608563293D0                                  
!                                                                       
      HTOTAL(16)= -2.35060497367449222283D0                                  
      HTOTAL(17)= -1.33584907401369694971D0                                  
      HTOTAL(18)= -0.43607741192761650868D0                                  
      HTOTAL(19)=  0.43607741192761650868D0                                  
      HTOTAL(20)=  1.33584907401369694971D0                                  
      HTOTAL(21)=  2.35060497367449222283D0                                  
!                                                                       
      HTOTAL(22)= -2.65196135683523349245D0                                  
      HTOTAL(23)= -1.67355162876747144503D0                                  
      HTOTAL(24)= -0.81628788285896466304D0                                  
      HTOTAL(25)= 0.0D0                                                         
      HTOTAL(26)=  0.81628788285896466304D0                                  
      HTOTAL(27)=  1.67355162876747144503D0                                  
      HTOTAL(28)=  2.65196135683523349245D0                                  
!                                                                       
      HTOTAL(29)= -2.93063742025724401922D0                                  
      HTOTAL(30)= -1.98165675669584292585D0                                  
      HTOTAL(31)= -1.15719371244678019472D0                                  
      HTOTAL(32)= -0.38118699020732211685D0                                  
      HTOTAL(33)=  0.38118699020732211685D0                                  
      HTOTAL(34)=  1.15719371244678019472D0                                  
      HTOTAL(35)=  1.98165675669584292585D0                                  
      HTOTAL(36)=  2.93063742025724401922D0                                  
!
      HTOTAL(37)= -3.19099320178152760723D0                                  
      HTOTAL(38)= -2.26658058453184311180D0                                  
      HTOTAL(39)= -1.46855328921666793167D0                                  
      HTOTAL(40)= -0.72355101875283757332D0                                  
      HTOTAL(41)= 0.0D0                                                         
      HTOTAL(42)=  0.72355101875283757332D0                                  
      HTOTAL(43)=  1.46855328921666793167D0                                  
      HTOTAL(44)=  2.26658058453184311180D0                                  
      HTOTAL(45)=  3.19099320178152760723D0                                  
!                                                                       
      HTOTAL(46)=  -3.43615911883773760333D0                                
      HTOTAL(47)=  -2.53273167423278979641D0                                
      HTOTAL(48)=  -1.75668364929988177345D0                                
      HTOTAL(49)=  -1.03661082978951365418D0                                
      HTOTAL(50)=  -0.34290132722370460879D0                                
      HTOTAL(51)=   0.34290132722370460879D0                                
      HTOTAL(52)=   1.03661082978951365418D0                                
      HTOTAL(53)=   1.75668364929988177345D0                                
      HTOTAL(54)=   2.53273167423278979641D0                                
      HTOTAL(55)=  3.43615911883773760333D0                                
!                                                                       
      WTOTAL(1)= 1.77245385090551602730D0  ! SQRT(PI)                       
!                                                                       
      WTOTAL(2)= 8.86226925452758013649D-01                                   
      WTOTAL(3)= 8.86226925452758013649D-01                                   
!                                                                       
      WTOTAL(4)= 2.95408975150919337883D-01                                   
      WTOTAL(5)= 1.18163590060367735153D0                                   
      WTOTAL(6)= 2.95408975150919337883D-01                                   
!                                                                       
      WTOTAL(7)= 8.13128354472451771430D-02                                   
      WTOTAL(8)= 8.04914090005512836506D-01                                   
      WTOTAL(9)= 8.04914090005512836506D-01                                   
      WTOTAL(10)= 8.13128354472451771430D-02                                   
!                                                                       
      WTOTAL(11)= 1.99532420590459132077D-02                                   
      WTOTAL(12)= 3.93619323152241159828D-01                                   
      WTOTAL(13)= 9.45308720482941881226D-01                                   
      WTOTAL(14)= 3.93619323152241159828D-01                                   
      WTOTAL(15)= 1.99532420590459132077D-02                                   
!                                                                       
      WTOTAL(16)= 4.53000990550884564086D-03                                   
      WTOTAL(17)= 1.57067320322856643916D-01                                   
      WTOTAL(18)= 7.24629595224392524092D-01                                   
      WTOTAL(19)= 7.24629595224392524092D-01                                   
      WTOTAL(20)= 1.57067320322856643916D-01                                   
      WTOTAL(21)= 4.53000990550884564086D-03                                   
!                                                                       
      WTOTAL(22)= 9.71781245099519154149D-04                                   
      WTOTAL(23)= 5.45155828191270305922D-02                                   
      WTOTAL(24)= 4.25607252610127800520D-01                                   
      WTOTAL(25)= 8.10264617556807326765D-01                                   
      WTOTAL(26)= 4.25607252610127800520D-01                                   
      WTOTAL(27)= 5.45155828191270305922D-02                                   
      WTOTAL(28)= 9.71781245099519154149D-04                                   
!                                                                      
      WTOTAL(29)= 1.99604072211367619206D-04                                   
      WTOTAL(30)= 1.70779830074134754562D-02                                   
      WTOTAL(31)= 2.07802325814891879543D-01                                   
      WTOTAL(32)= 6.61147012558241291030D-01                                   
      WTOTAL(33)= 6.61147012558241291030D-01                                   
      WTOTAL(34)= 2.07802325814891879543D-01                                   
      WTOTAL(35)= 1.70779830074134754562D-02                                   
      WTOTAL(36)= 1.99604072211367619206D-04                                   
!                                                                       
      WTOTAL(37)= 3.96069772632643819046D-05                                   
      WTOTAL(38)= 4.94362427553694721722D-03                                   
      WTOTAL(39)= 8.84745273943765732880D-02                                   
      WTOTAL(40)= 4.32651559002555750200D-01                                   
      WTOTAL(41)= 7.20235215606050957124D-01                                   
      WTOTAL(42)= 4.32651559002555750200D-01                                   
      WTOTAL(43)= 8.84745273943765732880D-02                                   
      WTOTAL(44)= 4.94362427553694721722D-03                                   
      WTOTAL(45)= 3.96069772632643819046D-05                                   
!                                                                       
      WTOTAL(46)= 7.64043285523262062916D-06                                 
      WTOTAL(47)= 1.34364574678123269220D-03                                 
      WTOTAL(48)= 3.38743944554810631362D-02                                 
      WTOTAL(49)= 2.40138611082314686417D-01                                 
      WTOTAL(50)= 6.10862633735325798784D-01                                 
      WTOTAL(51)= 6.10862633735325798784D-01                                 
      WTOTAL(52)= 2.40138611082314686417D-01                                 
      WTOTAL(53)= 3.38743944554810631362D-02                                 
      WTOTAL(54)= 1.34364574678123269220D-03                                 
      WTOTAL(55)= 7.64043285523262062916D-06                                 
!-----------------------------------------------------------------------                                                                       
      RETURN                                                            
      END        

! STVINT                                           
      SUBROUTINE STVINT(HH,WW)                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION HH(28),WW(28),MINIMO(7),MAXIMO(7)  
      COMMON/STV/XINTT,YINTT,ZINTT,TAA,X0X0,Y0Y0,Z0Z0,                  &
                 XIXI,YIYI,ZIZI,XJXJ,YJYJ,ZJZJ,NINI,NJNJ      
      DATA MINIMO /1,2,4,7,11,16,22/                                       
      DATA MAXIMO /1,3,6,10,15,21,28/                                      
!-----------------------------------------------------------------------
!     Gauss-Hermite Quadrature
!-----------------------------------------------------------------------
      XINTT = 0.0D0                                                       
      YINTT = 0.0D0                                                       
      ZINTT = 0.0D0                                                       
      NPTS = (NINI+NJNJ-2)/2+1                                              
      IMIN = MINIMO(NPTS)                                                  
      IMAX = MAXIMO(NPTS)                                                  
      DO I = IMIN,IMAX                                              
       DUM = WW(I)                                                     
       PX = DUM                                                       
       PY = DUM                                                       
       PZ = DUM                                                       
       DUM = HH(I)*TAA                                                   
       IF(NINI>1) THEN                                               
        AX = DUM+X0X0-XIXI                                              
        AY = DUM+Y0Y0-YIYI                                              
        AZ = DUM+Z0Z0-ZIZI                                              
        GO TO (16,15,14,13,12,11,10),NINI                         
   10   PX = PX*AX                                                  
        PY = PY*AY                                                  
        PZ = PZ*AZ                                                  
   11   PX = PX*AX                                                  
        PY = PY*AY                                                  
        PZ = PZ*AZ                                                  
   12   PX = PX*AX                                                  
        PY = PY*AY                                                  
        PZ = PZ*AZ                                                  
   13   PX = PX*AX                                                  
        PY = PY*AY                                                  
        PZ = PZ*AZ                                                  
   14   PX = PX*AX                                                  
        PY = PY*AY                                                  
        PZ = PZ*AZ                                                  
   15   PX = PX*AX                                                  
        PY = PY*AY                                                  
        PZ = PZ*AZ                                                  
       END IF                                                         
   16  CONTINUE                                                       
       IF(NJNJ>1) THEN                                               
        BX = DUM+X0X0-XJXJ                                              
        BY = DUM+Y0Y0-YJYJ                                              
        BZ = DUM+Z0Z0-ZJZJ                                              
        GO TO (27,26,25,24,23,22,21,20),NJNJ                  
   20   PX = PX*BX                                               
        PY = PY*BY                                               
        PZ = PZ*BZ                                               
   21   PX = PX*BX                                               
        PY = PY*BY                                               
        PZ = PZ*BZ                                               
   22   PX = PX*BX                                               
        PY = PY*BY                                               
        PZ = PZ*BZ                                               
   23   PX = PX*BX                                               
        PY = PY*BY                                               
        PZ = PZ*BZ                                               
   24   PX = PX*BX                                               
        PY = PY*BY                                               
        PZ = PZ*BZ                                               
   25   PX = PX*BX                                               
        PY = PY*BY                                               
        PZ = PZ*BZ                                               
   26   PX = PX*BX                                               
        PY = PY*BY                                               
        PZ = PZ*BZ                                               
       END IF                                                         
   27  CONTINUE                                                       
       XINTT = XINTT + PX                                           
       YINTT = YINTT + PY                                           
       ZINTT = ZINTT + PZ                                           
      END DO
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               

! SETCONI
      SUBROUTINE SETCONI(CONI,IPP,CS,CP,CD,CF,CG,CH,CI,NPRIMI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION CONI(*),KARTEN(0:6)
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: CS,CP,CD,CF,CG,CH,CI
      DATA KARTEN/1,3,6,10,15,21,28/
!-----------------------------------------------------------------------
      CONI(1) = CS(IPP)
      CALL DACOPY(KARTEN(1),CP(IPP),CONI(2),1)
      CALL DACOPY(KARTEN(2),CD(IPP),CONI(5),1)
      CALL DACOPY(KARTEN(3),CF(IPP),CONI(11),1)
      CALL DACOPY(KARTEN(4),CG(IPP),CONI(21),1)
      CALL DACOPY(KARTEN(5),CH(IPP),CONI(36),1)
      CALL DACOPY(KARTEN(6),CI(IPP),CONI(57),1)
!-----------------------------------------------------------------------
      RETURN
      END
            
! DACOPY                                           
      SUBROUTINE DACOPY(N,DA,DX,INCX)                                   
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION DX(*)                                                   
!-----------------------------------------------------------------------
!     INITIALISES A VECTOR WITH A CONSTANT: DX(I) <== DA                                                
!-----------------------------------------------------------------------
      IF(N<=0)RETURN                                                  
      IF(INCX/=1)THEN
       IX = 1                                                            
       IF(INCX<0)IX=(-N+1)*INCX + 1                                 
       DO I = 1,N                                                     
        DX(IX) = DA                                                     
        IX = IX + INCX                                                  
       END DO
       RETURN                                                            
      END IF
      M = MOD(N,7)                                                      
      IF(M/=0)THEN
       DO I = 1,M                                                     
        DX(I)=DA                                                        
       END DO
       IF(N<7)RETURN
      END IF                                             
      MP1 = M + 1                                                       
      DO I = MP1,N,7                                                 
       DX(I)=DA                                                        
       DX(I + 1)=DA                                                    
       DX(I + 2)=DA                                                    
       DX(I + 3)=DA                                                    
       DX(I + 4)=DA                                                    
       DX(I + 5)=DA                                                    
       DX(I + 6)=DA                                                    
      END DO
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               

! Calling by HSandT & ERISPDFGHIL: RT123 
      SUBROUTINE RT123
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /ROOT  / X,U(13),W(13),NROOTS
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
      DATA R12,PIE4/2.75255128608411D-01, 7.85398163397448D-01/
      DATA R22,W22/ 2.72474487139158D+00, 9.17517095361369D-02/
      DATA R13/     1.90163509193487D-01/
      DATA R23,W23/ 1.78449274854325D+00, 1.77231492083829D-01/
      DATA R33,W33/ 5.52534374226326D+00, 5.11156880411248D-03/
!
      IF (X > 5.0D+00) GO TO 400
      IF (X > 1.0D+00) GO TO 280
      IF (X > 3.0D-07) GO TO 180
!     X IS APPROXIMATELY 0.0d0.         NROOTS=1,2, OR 3
      IF (NROOTS-2) 120,140,160
  120 RT1 = 0.5D+00 -X/5.0D+00
      WW1 = 1.0D+00 -X/3.0D+00
      RETURN
  140 RT1 = 1.30693606237085D-01 -2.90430236082028D-02 *X
      RT2 = 2.86930639376291D+00 -6.37623643058102D-01 *X
      WW1 = 6.52145154862545D-01 -1.22713621927067D-01 *X
      WW2 = 3.47854845137453D-01 -2.10619711404725D-01 *X
      RETURN
  160 RT1 = 6.03769246832797D-02 -9.28875764357368D-03 *X
      RT2 = 7.76823355931043D-01 -1.19511285527878D-01 *X
      RT3 = 6.66279971938567D+00 -1.02504611068957D+00 *X
      WW1 = 4.67913934572691D-01 -5.64876917232519D-02 *X
      WW2 = 3.60761573048137D-01 -1.49077186455208D-01 *X
      WW3 = 1.71324492379169D-01 -1.27768455150979D-01 *X
      RETURN
!     X = 0.0 TO 1.0                   NROOTS=1,2, OR 3
  180 IF (NROOTS == 3) GO TO 220
      F1 = ((((((((-8.36313918003957D-08*X+1.21222603512827D-06 )*X-    &
           1.15662609053481D-05 )*X+9.25197374512647D-05 )*X-           &
           6.40994113129432D-04 )*X+3.78787044215009D-03 )*X-           &
           1.85185172458485D-02 )*X+7.14285713298222D-02 )*X-           &
           1.99999999997023D-01 )*X+3.33333333333318D-01                 
      WW1 = (X+X)*F1+EXP(-X)                                             
      IF (NROOTS == 2) GO TO 200                                         
      RT1 = F1/(WW1-F1)                                                  
      RETURN                                                             
  200 RT1 = (((((((-2.35234358048491D-09*X+2.49173650389842D-08)*X-     &
           4.558315364581D-08)*X-2.447252174587D-06)*X+                 &
           4.743292959463D-05)*X-5.33184749432408D-04 )*X+              &
           4.44654947116579D-03 )*X-2.90430236084697D-02 )*X+           &
           1.30693606237085D-01                                          
      RT2 = (((((((-2.47404902329170D-08*X+2.36809910635906D-07)*X+     &
           1.835367736310D-06)*X-2.066168802076D-05)*X-                 &
           1.345693393936D-04)*X-5.88154362858038D-05 )*X+              &
           5.32735082098139D-02 )*X-6.37623643056745D-01 )*X+           &
           2.86930639376289D+00                                          
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)                    
      WW1 = WW1-WW2                                                      
      RETURN                                                             
  220 RT1 = ((((((-5.10186691538870D-10*X+2.40134415703450D-08)*X-      &
           5.01081057744427D-07 )*X+7.58291285499256D-06 )*X-           &
           9.55085533670919D-05 )*X+1.02893039315878D-03 )*X-           &
           9.28875764374337D-03 )*X+6.03769246832810D-02                 
      RT2 = ((((((-1.29646524960555D-08*X+7.74602292865683D-08)*X+      &
           1.56022811158727D-06 )*X-1.58051990661661D-05 )*X-           &
           3.30447806384059D-04 )*X+9.74266885190267D-03 )*X-           &
           1.19511285526388D-01 )*X+7.76823355931033D-01                 
      RT3 = ((((((-9.28536484109606D-09*X-3.02786290067014D-07)*X-      &
           2.50734477064200D-06 )*X-7.32728109752881D-06 )*X+           &
           2.44217481700129D-04 )*X+4.94758452357327D-02 )*X-           &
           1.02504611065774D+00 )*X+6.66279971938553D+00                 
      F2 = ((((((((-7.60911486098850D-08*X+1.09552870123182D-06 )*X-    &
           1.03463270693454D-05 )*X+8.16324851790106D-05 )*X-           &
           5.55526624875562D-04 )*X+3.20512054753924D-03 )*X-           &
           1.51515139838540D-02 )*X+5.55555554649585D-02 )*X-           &
           1.42857142854412D-01 )*X+1.99999999999986D-01                 
  240 E = EXP(-X)                                                        
      F1 = ((X+X)*F2+E)/3.0D+00                                          
      WW1 = (X+X)*F1+E                                                   
  260 T1 = RT1/(RT1+1.0D+00)                                             
      T2 = RT2/(RT2+1.0D+00)                                             
      T3 = RT3/(RT3+1.0D+00)                                             
      A2 = F2-T1*F1                                                      
      A1 = F1-T1*WW1                                                     
      WW3 = (A2-T2*A1)/((T3-T2)*(T3-T1))                                 
      WW2 = (T3*A1-A2)/((T3-T2)*(T2-T1))                                 
      WW1 = WW1-WW2-WW3                                                  
      RETURN                                                             
  280 IF (X > 3.0D+00) GO TO 340                                         
!     X = 1.0 TO 3.0                   NROOTS=1,2, OR 3                  
      Y = X-2.0D+00                                                      
      IF (NROOTS == 3) GO TO 320                                         
      F1 = ((((((((((-1.61702782425558D-10*Y+1.96215250865776D-09 )*Y-  &
           2.14234468198419D-08 )*Y+2.17216556336318D-07 )*Y-           &
           1.98850171329371D-06 )*Y+1.62429321438911D-05 )*Y-           &
           1.16740298039895D-04 )*Y+7.24888732052332D-04 )*Y-           &
           3.79490003707156D-03 )*Y+1.61723488664661D-02 )*Y-           &
           5.29428148329736D-02 )*Y+1.15702180856167D-01                 
      WW1 = (X+X)*F1+EXP(-X)                                             
      IF (NROOTS == 2) GO TO 300                                         
      RT1 = F1/(WW1-F1)                                                  
      RETURN                                                             
  300 RT1 = (((((((((-6.36859636616415D-12*Y+8.47417064776270D-11)*Y-   &
           5.152207846962D-10)*Y-3.846389873308D-10)*Y+                 &
           8.472253388380D-08)*Y-1.85306035634293D-06 )*Y+              &
           2.47191693238413D-05 )*Y-2.49018321709815D-04 )*Y+           &
           2.19173220020161D-03 )*Y-1.63329339286794D-02 )*Y+           &
           8.68085688285261D-02                                          
      RT2 = ((((((((( 1.45331350488343D-10*Y+2.07111465297976D-09)*Y-   &
           1.878920917404D-08)*Y-1.725838516261D-07)*Y+                 &
           2.247389642339D-06)*Y+9.76783813082564D-06 )*Y-              &
           1.93160765581969D-04 )*Y-1.58064140671893D-03 )*Y+           &
           4.85928174507904D-02 )*Y-4.30761584997596D-01 )*Y+           &
           1.80400974537950D+00                                          
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)                    
      WW1 = WW1-WW2                                                      
      RETURN                                                             
  320 RT1 = (((((((( 1.44687969563318D-12*Y+4.85300143926755D-12)*Y-    &
           6.55098264095516D-10 )*Y+1.56592951656828D-08 )*Y-           &
           2.60122498274734D-07 )*Y+3.86118485517386D-06 )*Y-           &
           5.13430986707889D-05 )*Y+6.03194524398109D-04 )*Y-           &
           6.11219349825090D-03 )*Y+4.52578254679079D-02                 
      RT2 = ((((((( 6.95964248788138D-10*Y-5.35281831445517D-09)*Y-     &
           6.745205954533D-08)*Y+1.502366784525D-06)*Y+                 &
           9.923326947376D-07)*Y-3.89147469249594D-04 )*Y+              &
           7.51549330892401D-03 )*Y-8.48778120363400D-02 )*Y+           &
           5.73928229597613D-01                                          
      RT3 = ((((((((-2.81496588401439D-10*Y+3.61058041895031D-09)*Y+    &
           4.53631789436255D-08 )*Y-1.40971837780847D-07 )*Y-           &
           6.05865557561067D-06 )*Y-5.15964042227127D-05 )*Y+           &
           3.34761560498171D-05 )*Y+5.04871005319119D-02 )*Y-           &
           8.24708946991557D-01 )*Y+4.81234667357205D+00                 
      F2 = ((((((((((-1.48044231072140D-10*Y+1.78157031325097D-09 )*Y-  &
           1.92514145088973D-08 )*Y+1.92804632038796D-07 )*Y-           &
           1.73806555021045D-06 )*Y+1.39195169625425D-05 )*Y-           &
           9.74574633246452D-05 )*Y+5.83701488646511D-04 )*Y-           &
           2.89955494844975D-03 )*Y+1.13847001113810D-02 )*Y-           &
           3.23446977320647D-02 )*Y+5.29428148329709D-02                 
      GO TO 240                                                          
!     X = 3.0 TO 5.0                   NROOTS =1,2, OR 3                 
  340 Y = X-4.0D+00                                                      
      IF (NROOTS == 3) GO TO 380                                         
      F1 = ((((((((((-2.62453564772299D-11*Y+3.24031041623823D-10 )*Y-  &
           3.614965656163D-09)*Y+3.760256799971D-08)*Y-                 &
           3.553558319675D-07)*Y+3.022556449731D-06)*Y-                 &
           2.290098979647D-05)*Y+1.526537461148D-04)*Y-                 &
           8.81947375894379D-04 )*Y+4.33207949514611D-03 )*Y-           &
           1.75257821619926D-02 )*Y+5.28406320615584D-02                 
      WW1 = (X+X)*F1+EXP(-X)                                             
      IF (NROOTS == 2) GO TO 360                                         
      RT1 = F1/(WW1-F1)                                                  
      RETURN                                                             
  360 RT1 = ((((((((-4.11560117487296D-12*Y+7.10910223886747D-11)*Y-    &
           1.73508862390291D-09 )*Y+5.93066856324744D-08 )*Y-           &
           9.76085576741771D-07 )*Y+1.08484384385679D-05 )*Y-           &
           1.12608004981982D-04 )*Y+1.16210907653515D-03 )*Y-           &
           9.89572595720351D-03 )*Y+6.12589701086408D-02                 
      RT2 = (((((((((-1.80555625241001D-10*Y+5.44072475994123D-10)*Y+   &
           1.603498045240D-08)*Y-1.497986283037D-07)*Y-                 &
           7.017002532106D-07)*Y+1.85882653064034D-05 )*Y-              &
           2.04685420150802D-05 )*Y-2.49327728643089D-03 )*Y+           &
           3.56550690684281D-02 )*Y-2.60417417692375D-01 )*Y+           &
           1.12155283108289D+00                                          
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)                    
      WW1 = WW1-WW2                                                      
      RETURN                                                             
  380 RT1 = ((((((( 1.44265709189601D-11*Y-4.66622033006074D-10)*Y+     &
           7.649155832025D-09)*Y-1.229940017368D-07)*Y+                 &
           2.026002142457D-06)*Y-2.87048671521677D-05 )*Y+              &
           3.70326938096287D-04 )*Y-4.21006346373634D-03 )*Y+           &
           3.50898470729044D-02                                          
      RT2 = ((((((((-2.65526039155651D-11*Y+1.97549041402552D-10)*Y+    &
           2.15971131403034D-09 )*Y-7.95045680685193D-08 )*Y+           &
           5.15021914287057D-07 )*Y+1.11788717230514D-05 )*Y-           &
           3.33739312603632D-04 )*Y+5.30601428208358D-03 )*Y-           &
           5.93483267268959D-02 )*Y+4.31180523260239D-01                 
      RT3 = ((((((((-3.92833750584041D-10*Y-4.16423229782280D-09)*Y+    &
           4.42413039572867D-08 )*Y+6.40574545989551D-07 )*Y-           &
           3.05512456576552D-06 )*Y-1.05296443527943D-04 )*Y-           &
           6.14120969315617D-04 )*Y+4.89665802767005D-02 )*Y-           &
           6.24498381002855D-01 )*Y+3.36412312243724D+00                 
      F2 = ((((((((((-2.36788772599074D-11*Y+2.89147476459092D-10 )*Y-  &
           3.18111322308846D-09 )*Y+3.25336816562485D-08 )*Y-           &
           3.00873821471489D-07 )*Y+2.48749160874431D-06 )*Y-           &
           1.81353179793672D-05 )*Y+1.14504948737066D-04 )*Y-           &
           6.10614987696677D-04 )*Y+2.64584212770942D-03 )*Y-           &
           8.66415899015349D-03 )*Y+1.75257821619922D-02                 
      GO TO 240                                                          
  400 IF (X > 15.0D+00) GO TO 560                                        
      E = EXP(-X)                                                        
      IF (X > 10.0D+00) GO TO 480                                        
!     X = 5.0 TO 10.0                  NROOTS =1,2, OR 3                 
      WW1 = (((((( 4.6897511375022D-01/X-6.9955602298985D-01)/X +       &
           5.3689283271887D-01)/X-3.2883030418398D-01)/X +              &
           2.4645596956002D-01)/X-4.9984072848436D-01)/X -              &
           3.1501078774085D-06)*E + SQRT(PIE4/X)                         
      F1 = (WW1-E)/(X+X)                                                 
      IF (NROOTS-2) 420,440,460                                          
  420 RT1 = F1/(WW1-F1)                                                  
      RETURN                                                             
  440 Y = X-7.5D+00                                                      
      RT1 = (((((((((((((-1.43632730148572D-16*Y+2.38198922570405D-16)* &
           Y+1.358319618800D-14)*Y-7.064522786879D-14)*Y-               &
           7.719300212748D-13)*Y+7.802544789997D-12)*Y+                 &
           6.628721099436D-11)*Y-1.775564159743D-09)*Y+                 &
           1.713828823990D-08)*Y-1.497500187053D-07)*Y+                 &
           2.283485114279D-06)*Y-3.76953869614706D-05 )*Y+              &
           4.74791204651451D-04 )*Y-4.60448960876139D-03 )*Y+           &
           3.72458587837249D-02                                          
      RT2 = (((((((((((( 2.48791622798900D-14*Y-1.36113510175724D-13)*Y-&
           2.224334349799D-12)*Y+4.190559455515D-11)*Y-                 &
           2.222722579924D-10)*Y-2.624183464275D-09)*Y+                 &
           6.128153450169D-08)*Y-4.383376014528D-07)*Y-                 &
           2.49952200232910D-06 )*Y+1.03236647888320D-04 )*Y-           &
           1.44614664924989D-03 )*Y+1.35094294917224D-02 )*Y-           &
           9.53478510453887D-02 )*Y+5.44765245686790D-01                 
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)                    
      WW1 = WW1-WW2                                                      
      RETURN                                                             
  460 F2 = (F1+F1+F1-E)/(X+X)                                            
      Y = X-7.5D+00                                                      
      RT1 = ((((((((((( 5.74429401360115D-16*Y+7.11884203790984D-16)*Y- &
           6.736701449826D-14)*Y-6.264613873998D-13)*Y+                 &
           1.315418927040D-11)*Y-4.23879635610964D-11 )*Y+              &
           1.39032379769474D-09 )*Y-4.65449552856856D-08 )*Y+           &
           7.34609900170759D-07 )*Y-1.08656008854077D-05 )*Y+           &
           1.77930381549953D-04 )*Y-2.39864911618015D-03 )*Y+           &
           2.39112249488821D-02                                          
      RT2 = ((((((((((( 1.13464096209120D-14*Y+6.99375313934242D-15)*Y- &
           8.595618132088D-13)*Y-5.293620408757D-12)*Y-                 &
           2.492175211635D-11)*Y+2.73681574882729D-09 )*Y-              &
           1.06656985608482D-08 )*Y-4.40252529648056D-07 )*Y+           &
           9.68100917793911D-06 )*Y-1.68211091755327D-04 )*Y+           &
           2.69443611274173D-03 )*Y-3.23845035189063D-02 )*Y+           &
           2.75969447451882D-01                                          
      RT3 = (((((((((((( 6.66339416996191D-15*Y+1.84955640200794D-13)*Y-&
           1.985141104444D-12)*Y-2.309293727603D-11)*Y+                 &
           3.917984522103D-10)*Y+1.663165279876D-09)*Y-                 &
           6.205591993923D-08)*Y+8.769581622041D-09)*Y+                 &
           8.97224398620038D-06 )*Y-3.14232666170796D-05 )*Y-           &
           1.83917335649633D-03 )*Y+3.51246831672571D-02 )*Y-           &
           3.22335051270860D-01 )*Y+1.73582831755430D+00                 
      GO TO 260                                                          
!     X = 10.0 TO 15.0                 NROOTS=1,2, OR 3                  
  480 WW1 = (((-1.8784686463512D-01/X+2.2991849164985D-01)/X -          &
           4.9893752514047D-01)/X-2.1916512131607D-05)*E + SQRT(PIE4/X)  
      F1 = (WW1-E)/(X+X)                                                 
      IF (NROOTS-2) 500,520,540                                          
  500 RT1 = F1/(WW1-F1)                                                  
      RETURN                                                             
  520 RT1 = ((((-1.01041157064226D-05*X+1.19483054115173D-03)*X -       &
           6.73760231824074D-02)*X+1.25705571069895D+00)*X + (((-       &
           8.57609422987199D+03/X+5.91005939591842D+03)/X -             &
           1.70807677109425D+03)/X+2.64536689959503D+02)/X -            &
           2.38570496490846D+01)*E + R12/(X-R12)                         
      RT2 = ((( 3.39024225137123D-04*X-9.34976436343509D-02)*X -        &
           4.22216483306320D+00)*X + (((-2.08457050986847D+03/X -       &
           1.04999071905664D+03)/X+3.39891508992661D+02)/X -            &
           1.56184800325063D+02)/X+8.00839033297501D+00)*E + R22/(X-R22) 
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)                    
      WW1 = WW1-WW2                                                      
      RETURN                                                             
  540 F2 = (F1+F1+F1-E)/(X+X)                                            
      Y = X-12.5D+00                                                     
      RT1 = ((((((((((( 4.42133001283090D-16*Y-2.77189767070441D-15)*Y- &
           4.084026087887D-14)*Y+5.379885121517D-13)*Y+                 &
           1.882093066702D-12)*Y-8.67286219861085D-11 )*Y+              &
           7.11372337079797D-10 )*Y-3.55578027040563D-09 )*Y+           &
           1.29454702851936D-07 )*Y-4.14222202791434D-06 )*Y+           &
           8.04427643593792D-05 )*Y-1.18587782909876D-03 )*Y+           &
           1.53435577063174D-02                                          
      RT2 = ((((((((((( 6.85146742119357D-15*Y-1.08257654410279D-14)*Y- &
           8.579165965128D-13)*Y+6.642452485783D-12)*Y+                 &
           4.798806828724D-11)*Y-1.13413908163831D-09 )*Y+              &
           7.08558457182751D-09 )*Y-5.59678576054633D-08 )*Y+           &
           2.51020389884249D-06 )*Y-6.63678914608681D-05 )*Y+           &
           1.11888323089714D-03 )*Y-1.45361636398178D-02 )*Y+           &
           1.65077877454402D-01                                          
      RT3 = (((((((((((( 3.20622388697743D-15*Y-2.73458804864628D-14)*Y-&
           3.157134329361D-13)*Y+8.654129268056D-12)*Y-                 &
           5.625235879301D-11)*Y-7.718080513708D-10)*Y+                 &
           2.064664199164D-08)*Y-1.567725007761D-07)*Y-                 &
           1.57938204115055D-06 )*Y+6.27436306915967D-05 )*Y-           &
           1.01308723606946D-03 )*Y+1.13901881430697D-02 )*Y-           &
           1.01449652899450D-01 )*Y+7.77203937334739D-01                 
      GO TO 260                                                          
  560 IF (X > 33.0D+00) GO TO 660                                        
!     X = 15.0 TO 33.0                 NROOTS=1,2, OR 3                  
      E = EXP(-X)                                                        
      WW1 = (( 1.9623264149430D-01/X-4.9695241464490D-01)/X -           &
           6.0156581186481D-05)*E + SQRT(PIE4/X)                         
      F1 = (WW1-E)/(X+X)                                                 
      IF (NROOTS-2) 580,600,620                                          
  580 RT1 = F1/(WW1-F1)                                                  
      RETURN                                                             
  600 RT1 = ((((-1.14906395546354D-06*X+1.76003409708332D-04)*X -       &
           1.71984023644904D-02)*X-1.37292644149838D-01)*X + (-         &
           4.75742064274859D+01/X+9.21005186542857D+00)/X -             &
           2.31080873898939D-02)*E + R12/(X-R12)                         
      RT2 = ((( 3.64921633404158D-04*X-9.71850973831558D-02)*X -        &
           4.02886174850252D+00)*X + (-1.35831002139173D+02/X -         &
           8.66891724287962D+01)/X+2.98011277766958D+00)*E + R22/(X-R22) 
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)                    
      WW1 = WW1-WW2                                                      
      RETURN                                                             
  620 F2 = (F1+F1+F1-E)/(X+X)                                            
      IF (X > 20.0D+00) GO TO 640                                        
      RT1 = ((((((-2.43270989903742D-06*X+3.57901398988359D-04)*X -     &
           2.34112415981143D-02)*X+7.81425144913975D-01)*X -            &
           1.73209218219175D+01)*X+2.43517435690398D+02)*X + (-         &
           1.97611541576986D+04/X+9.82441363463929D+03)/X -             &
           2.07970687843258D+03)*E + R13/(X-R13)                         
      RT2 = (((((-2.62627010965435D-04*X+3.49187925428138D-02)*X -      &
           3.09337618731880D+00)*X+1.07037141010778D+02)*X -            &
           2.36659637247087D+03)*X + ((-2.91669113681020D+06/X +        &
           1.41129505262758D+06)/X-2.91532335433779D+05)/X +            &
           3.35202872835409D+04)*E + R23/(X-R23)                         
      RT3 = ((((( 9.31856404738601D-05*X-2.87029400759565D-02)*X -      &
           7.83503697918455D-01)*X-1.84338896480695D+01)*X +            &
           4.04996712650414D+02)*X + (-1.89829509315154D+05/X +         &
           5.11498390849158D+04)/X-6.88145821789955D+03)*E + R33/(X-R33) 
      GO TO 260                                                          
  640 RT1 = ((((-4.97561537069643D-04*X-5.00929599665316D-02)*X +       &
           1.31099142238996D+00)*X-1.88336409225481D+01)*X -            &
           6.60344754467191D+02 /X+1.64931462413877D+02)*E + R13/(X-R13) 
      RT2 = ((((-4.48218898474906D-03*X-5.17373211334924D-01)*X +       &
           1.13691058739678D+01)*X-1.65426392885291D+02)*X -            &
           6.30909125686731D+03 /X+1.52231757709236D+03)*E + R23/(X-R23) 
      RT3 = ((((-1.38368602394293D-02*X-1.77293428863008D+00)*X +       &
           1.73639054044562D+01)*X-3.57615122086961D+02)*X -            &
           1.45734701095912D+04 /X+2.69831813951849D+03)*E + R33/(X-R33) 
      GO TO 260                                                          
!     X = 33.0 TO INFINITY             NROOTS=1,2, OR 3                  
  660 WW1 = SQRT(PIE4/X)                                                 
      IF (NROOTS-2) 680,700,720                                          
  680 RT1 = 0.5D+00/(X-0.5D+00)                                          
      RETURN                                                             
  700 IF (X > 40.0D+00) GO TO 740                                        
      E = EXP(-X)                                                        
      RT1 = (-8.78947307498880D-01*X+1.09243702330261D+01)*E + R12/(X-  &
           R12)                                                          
      RT2 = (-9.28903924275977D+00*X+8.10642367843811D+01)*E + R22/(X-  &
           R22)                                                          
      WW2 = ( 4.46857389308400D+00*X-7.79250653461045D+01)*E + W22*WW1   
      WW1 = WW1-WW2                                                      
      RETURN                                                             
  720 IF (X > 47.0D+00) GO TO 760                                        
      E = EXP(-X)                                                        
      RT1 = ((-7.39058467995275D+00*X+3.21318352526305D+02)*X -         &
           3.99433696473658D+03)*E + R13/(X-R13)                         
      RT2 = ((-7.38726243906513D+01*X+3.13569966333873D+03)*X -         &
           3.86862867311321D+04)*E + R23/(X-R23)                         
      RT3 = ((-2.63750565461336D+02*X+1.04412168692352D+04)*X -         &
           1.28094577915394D+05)*E + R33/(X-R33)                         
      WW3 = ((( 1.52258947224714D-01*X-8.30661900042651D+00)*X +        &
           1.92977367967984D+02)*X-1.67787926005344D+03)*E + W33*WW1     
      WW2 = (( 6.15072615497811D+01*X-2.91980647450269D+03)*X +         &
           3.80794303087338D+04)*E + W23*WW1
      WW1 = WW1-WW2-WW3
      RETURN
  740 RT1 = R12/(X-R12)
      RT2 = R22/(X-R22)
      WW2 = W22*WW1
      WW1 = WW1-WW2
      RETURN
  760 RT1 = R13/(X-R13)
      RT2 = R23/(X-R23)
      RT3 = R33/(X-R33)
      WW2 = W23*WW1
      WW3 = W33*WW1
      WW1 = WW1-WW2-WW3
      RETURN
      END

! Calling by HSandT & ERISPDFGHIL: ROOT4 
      SUBROUTINE ROOT4
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /ROOT  / X,U(13),W(13),NROOTS
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
      DATA R14,PIE4/1.45303521503316D-01, 7.85398163397448D-01/
      DATA R24,W24/ 1.33909728812636D+00, 2.34479815323517D-01/
      DATA R34,W34/ 3.92696350135829D+00, 1.92704402415764D-02/
      DATA R44,W44/ 8.58863568901199D+00, 2.25229076750736D-04/
!
      IF (X > 15.0D+00) GO TO 180
      IF (X > 5.0D+00) GO TO 140
      IF (X > 1.0D+00) GO TO 120
      IF (X > 3.0D-07) GO TO 100
!     X IS APPROXIMATELY 0.0d0.                   NROOTS = 4
      RT1 = 3.48198973061471D-02 -4.09645850660395D-03 *X
      RT2 = 3.81567185080042D-01 -4.48902570656719D-02 *X
      RT3 = 1.73730726945891D+00 -2.04389090547327D-01 *X
      RT4 = 1.18463056481549D+01 -1.39368301742312D+00 *X
      WW1 = 3.62683783378362D-01 -3.13844305713928D-02 *X
      WW2 = 3.13706645877886D-01 -8.98046242557724D-02 *X
      WW3 = 2.22381034453372D-01 -1.29314370958973D-01 *X
      WW4 = 1.01228536290376D-01 -8.28299075414321D-02 *X
      RETURN
!
!     X=0.0 TO 1.0                               NROOTS = 4
  100 RT1 = ((((((-1.95309614628539D-10*X+5.19765728707592D-09)*X-      &
           1.01756452250573D-07 )*X+1.72365935872131D-06 )*X-           &
           2.61203523522184D-05 )*X+3.52921308769880D-04 )*X-           &
           4.09645850658433D-03 )*X+3.48198973061469D-02                 
      RT2 = (((((-1.89554881382342D-08*X+3.07583114342365D-07)*X+       &
           1.270981734393D-06)*X-1.417298563884D-04)*X+                 &
           3.226979163176D-03)*X-4.48902570678178D-02 )*X+              &
           3.81567185080039D-01                                          
      RT3 = (((((( 1.77280535300416D-09*X+3.36524958870615D-08)*X-      &
           2.58341529013893D-07 )*X-1.13644895662320D-05 )*X-           &
           7.91549618884063D-05 )*X+1.03825827346828D-02 )*X-           &
           2.04389090525137D-01 )*X+1.73730726945889D+00                 
      RT4 = (((((-5.61188882415248D-08*X-2.49480733072460D-07)*X+       &
           3.428685057114D-06)*X+1.679007454539D-04)*X+                 &
           4.722855585715D-02)*X-1.39368301737828D+00 )*X+              &
           1.18463056481543D+01                                          
      WW1 = ((((((-1.14649303201279D-08*X+1.88015570196787D-07)*X-      &
           2.33305875372323D-06 )*X+2.68880044371597D-05 )*X-           &
           2.94268428977387D-04 )*X+3.06548909776613D-03 )*X-           &
           3.13844305680096D-02 )*X+3.62683783378335D-01                 
      WW2 = ((((((((-4.11720483772634D-09*X+6.54963481852134D-08)*X-    &
           7.20045285129626D-07 )*X+6.93779646721723D-06 )*X-           &
           6.05367572016373D-05 )*X+4.74241566251899D-04 )*X-           &
           3.26956188125316D-03 )*X+1.91883866626681D-02 )*X-           &
           8.98046242565811D-02 )*X+3.13706645877886D-01                 
      WW3 = ((((((((-3.41688436990215D-08*X+5.07238960340773D-07)*X-    &
           5.01675628408220D-06 )*X+4.20363420922845D-05 )*X-           &
           3.08040221166823D-04 )*X+1.94431864731239D-03 )*X-           &
           1.02477820460278D-02 )*X+4.28670143840073D-02 )*X-           &
           1.29314370962569D-01 )*X+2.22381034453369D-01                 
      WW4 = ((((((((( 4.99660550769508D-09*X-7.94585963310120D-08)*X+   &
           8.359072409485D-07)*X-7.422369210610D-06)*X+                 &
           5.763374308160D-05)*X-3.86645606718233D-04 )*X+              &
           2.18417516259781D-03 )*X-9.99791027771119D-03 )*X+           &
           3.48791097377370D-02 )*X-8.28299075413889D-02 )*X+           &
           1.01228536290376D-01                                          
      RETURN
!                                                                        
!     X= 1.0 TO 5.0                              NROOTS = 4              
  120 Y = X-3.0D+00                                                      
      RT1 = (((((((((-1.48570633747284D-15*Y-1.33273068108777D-13)*Y+   &
           4.068543696670D-12)*Y-9.163164161821D-11)*Y+                 &
           2.046819017845D-09)*Y-4.03076426299031D-08 )*Y+              &
           7.29407420660149D-07 )*Y-1.23118059980833D-05 )*Y+           &
           1.88796581246938D-04 )*Y-2.53262912046853D-03 )*Y+           &
           2.51198234505021D-02                                          
      RT2 = ((((((((( 1.35830583483312D-13*Y-2.29772605964836D-12)*Y-   &
           3.821500128045D-12)*Y+6.844424214735D-10)*Y-                 &
           1.048063352259D-08)*Y+1.50083186233363D-08 )*Y+              &
           3.48848942324454D-06 )*Y-1.08694174399193D-04 )*Y+           &
           2.08048885251999D-03 )*Y-2.91205805373793D-02 )*Y+           &
           2.72276489515713D-01                                          
      RT3 = ((((((((( 5.02799392850289D-13*Y+1.07461812944084D-11)*Y-   &
           1.482277886411D-10)*Y-2.153585661215D-09)*Y+                 &
           3.654087802817D-08)*Y+5.15929575830120D-07 )*Y-              &
           9.52388379435709D-06 )*Y-2.16552440036426D-04 )*Y+           &
           9.03551469568320D-03 )*Y-1.45505469175613D-01 )*Y+           &
           1.21449092319186D+00                                          
      RT4 = (((((((((-1.08510370291979D-12*Y+6.41492397277798D-11)*Y+   &
           7.542387436125D-10)*Y-2.213111836647D-09)*Y-                 &
           1.448228963549D-07)*Y-1.95670833237101D-06 )*Y-              &
           1.07481314670844D-05 )*Y+1.49335941252765D-04 )*Y+           &
           4.87791531990593D-02 )*Y-1.10559909038653D+00 )*Y+           &
           8.09502028611780D+00                                          
      WW1 = ((((((((((-4.65801912689961D-14*Y+7.58669507106800D-13)*Y-  &
           1.186387548048D-11)*Y+1.862334710665D-10)*Y-                 &
           2.799399389539D-09)*Y+4.148972684255D-08)*Y-                 &
           5.933568079600D-07)*Y+8.168349266115D-06)*Y-                 &
           1.08989176177409D-04 )*Y+1.41357961729531D-03 )*Y-           &
           1.87588361833659D-02 )*Y+2.89898651436026D-01                 
      WW2 = ((((((((((((-1.46345073267549D-14*Y+2.25644205432182D-13)*Y-&
           3.116258693847D-12)*Y+4.321908756610D-11)*Y-                 &
           5.673270062669D-10)*Y+7.006295962960D-09)*Y-                 &
           8.120186517000D-08)*Y+8.775294645770D-07)*Y-                 &
           8.77829235749024D-06 )*Y+8.04372147732379D-05 )*Y-           &
           6.64149238804153D-04 )*Y+4.81181506827225D-03 )*Y-           &
           2.88982669486183D-02 )*Y+1.56247249979288D-01                 
      WW3 = ((((((((((((( 9.06812118895365D-15*Y-1.40541322766087D-13)* &
           Y+1.919270015269D-12)*Y-2.605135739010D-11)*Y+               &
           3.299685839012D-10)*Y-3.86354139348735D-09 )*Y+              &
           4.16265847927498D-08 )*Y-4.09462835471470D-07 )*Y+           &
           3.64018881086111D-06 )*Y-2.88665153269386D-05 )*Y+           &
           2.00515819789028D-04 )*Y-1.18791896897934D-03 )*Y+           &
           5.75223633388589D-03 )*Y-2.09400418772687D-02 )*Y+           &
           4.85368861938873D-02                                          
      WW4 = ((((((((((((((-9.74835552342257D-16*Y+1.57857099317175D-14)*&
           Y-2.249993780112D-13)*Y+3.173422008953D-12)*Y-               &
           4.161159459680D-11)*Y+5.021343560166D-10)*Y-                 &
           5.545047534808D-09)*Y+5.554146993491D-08)*Y-                 &
           4.99048696190133D-07 )*Y+3.96650392371311D-06 )*Y-           &
           2.73816413291214D-05 )*Y+1.60106988333186D-04 )*Y-           &
           7.64560567879592D-04 )*Y+2.81330044426892D-03 )*Y-           &
           7.16227030134947D-03 )*Y+9.66077262223353D-03                 
      RETURN
!                                                                        
  140 IF (X > 10.0D+00) GO TO 160                                        
!     X=5.0 TO 10.0                              NROOTS = 4              
      Y = X-7.5D+00                                                      
      RT1 = ((((((((( 4.64217329776215D-15*Y-6.27892383644164D-15)*Y+   &
           3.462236347446D-13)*Y-2.927229355350D-11)*Y+                 &
           5.090355371676D-10)*Y-9.97272656345253D-09 )*Y+              &
           2.37835295639281D-07 )*Y-4.60301761310921D-06 )*Y+           &
           8.42824204233222D-05 )*Y-1.37983082233081D-03 )*Y+           &
           1.66630865869375D-02                                          
      RT2 = ((((((((( 2.93981127919047D-14*Y+8.47635639065744D-13)*Y-   &
           1.446314544774D-11)*Y-6.149155555753D-12)*Y+                 &
           8.484275604612D-10)*Y-6.10898827887652D-08 )*Y+              &
           2.39156093611106D-06 )*Y-5.35837089462592D-05 )*Y+           &
           1.00967602595557D-03 )*Y-1.57769317127372D-02 )*Y+           &
           1.74853819464285D-01                                          
      RT3 = (((((((((( 2.93523563363000D-14*Y-6.40041776667020D-14)*Y-  &
           2.695740446312D-12)*Y+1.027082960169D-10)*Y-                 &
           5.822038656780D-10)*Y-3.159991002539D-08)*Y+                 &
           4.327249251331D-07)*Y+4.856768455119D-06)*Y-                 &
           2.54617989427762D-04 )*Y+5.54843378106589D-03 )*Y-           &
           7.95013029486684D-02 )*Y+7.20206142703162D-01                 
      RT4 = (((((((((((-1.62212382394553D-14*Y+7.68943641360593D-13)*Y+ &
           5.764015756615D-12)*Y-1.380635298784D-10)*Y-                 &
           1.476849808675D-09)*Y+1.84347052385605D-08 )*Y+              &
           3.34382940759405D-07 )*Y-1.39428366421645D-06 )*Y-           &
           7.50249313713996D-05 )*Y-6.26495899187507D-04 )*Y+           &
           4.69716410901162D-02 )*Y-6.66871297428209D-01 )*Y+           &
           4.11207530217806D+00                                          
      WW1 = ((((((((((-1.65995045235997D-15*Y+6.91838935879598D-14)*Y-  &
           9.131223418888D-13)*Y+1.403341829454D-11)*Y-                 &
           3.672235069444D-10)*Y+6.366962546990D-09)*Y-                 &
           1.039220021671D-07)*Y+1.959098751715D-06)*Y-                 &
           3.33474893152939D-05 )*Y+5.72164211151013D-04 )*Y-           &
           1.05583210553392D-02 )*Y+2.26696066029591D-01                 
      WW2 = ((((((((((((-3.57248951192047D-16*Y+6.25708409149331D-15)*Y-&
           9.657033089714D-14)*Y+1.507864898748D-12)*Y-                 &
           2.332522256110D-11)*Y+3.428545616603D-10)*Y-                 &
           4.698730937661D-09)*Y+6.219977635130D-08)*Y-                 &
           7.83008889613661D-07 )*Y+9.08621687041567D-06 )*Y-           &
           9.86368311253873D-05 )*Y+9.69632496710088D-04 )*Y-           &
           8.14594214284187D-03 )*Y+8.50218447733457D-02                 
      WW3 = ((((((((((((( 1.64742458534277D-16*Y-2.68512265928410D-15)* &
           Y+3.788890667676D-14)*Y-5.508918529823D-13)*Y+               &
           7.555896810069D-12)*Y-9.69039768312637D-11 )*Y+              &
           1.16034263529672D-09 )*Y-1.28771698573873D-08 )*Y+           &
           1.31949431805798D-07 )*Y-1.23673915616005D-06 )*Y+           &
           1.04189803544936D-05 )*Y-7.79566003744742D-05 )*Y+           &
           5.03162624754434D-04 )*Y-2.55138844587555D-03 )*Y+           &
           1.13250730954014D-02                                          
      WW4 = ((((((((((((((-1.55714130075679D-17*Y+2.57193722698891D-16)*&
           Y-3.626606654097D-15)*Y+5.234734676175D-14)*Y-               &
           7.067105402134D-13)*Y+8.793512664890D-12)*Y-                 &
           1.006088923498D-10)*Y+1.050565098393D-09)*Y-                 &
           9.91517881772662D-09 )*Y+8.35835975882941D-08 )*Y-           &
           6.19785782240693D-07 )*Y+3.95841149373135D-06 )*Y-           &
           2.11366761402403D-05 )*Y+9.00474771229507D-05 )*Y-           &
           2.78777909813289D-04 )*Y+5.26543779837487D-04                 
      RETURN
!                                                                        
!     X=10.0 TO 15.0                             NROOTS = 4              
  160 Y = X-12.5D+00                                                     
      RT1 = ((((((((((( 4.94869622744119D-17*Y+8.03568805739160D-16)*Y- &
           5.599125915431D-15)*Y-1.378685560217D-13)*Y+                 &
           7.006511663249D-13)*Y+1.30391406991118D-11 )*Y+              &
           8.06987313467541D-11 )*Y-5.20644072732933D-09 )*Y+           &
           7.72794187755457D-08 )*Y-1.61512612564194D-06 )*Y+           &
           4.15083811185831D-05 )*Y-7.87855975560199D-04 )*Y+           &
           1.14189319050009D-02                                          
      RT2 = ((((((((((( 4.89224285522336D-16*Y+1.06390248099712D-14)*Y- &
           5.446260182933D-14)*Y-1.613630106295D-12)*Y+                 &
           3.910179118937D-12)*Y+1.90712434258806D-10 )*Y+              &
           8.78470199094761D-10 )*Y-5.97332993206797D-08 )*Y+           &
           9.25750831481589D-07 )*Y-2.02362185197088D-05 )*Y+           &
           4.92341968336776D-04 )*Y-8.68438439874703D-03 )*Y+           &
           1.15825965127958D-01                                          
      RT3 = (((((((((( 6.12419396208408D-14*Y+1.12328861406073D-13)*Y-  &
           9.051094103059D-12)*Y-4.781797525341D-11)*Y+                 &
           1.660828868694D-09)*Y+4.499058798868D-10)*Y-                 &
           2.519549641933D-07)*Y+4.977444040180D-06)*Y-                 &
           1.25858350034589D-04 )*Y+2.70279176970044D-03 )*Y-           &
           3.99327850801083D-02 )*Y+4.33467200855434D-01                 
      RT4 = ((((((((((( 4.63414725924048D-14*Y-4.72757262693062D-14)*Y- &
           1.001926833832D-11)*Y+6.074107718414D-11)*Y+                 &
           1.576976911942D-09)*Y-2.01186401974027D-08 )*Y-              &
           1.84530195217118D-07 )*Y+5.02333087806827D-06 )*Y+           &
           9.66961790843006D-06 )*Y-1.58522208889528D-03 )*Y+           &
           2.80539673938339D-02 )*Y-2.78953904330072D-01 )*Y+           &
           1.82835655238235D+00                                          
      WW4 = ((((((((((((( 2.90401781000996D-18*Y-4.63389683098251D-17)* &
           Y+6.274018198326D-16)*Y-8.936002188168D-15)*Y+               &
           1.194719074934D-13)*Y-1.45501321259466D-12 )*Y+              &
           1.64090830181013D-11 )*Y-1.71987745310181D-10 )*Y+           &
           1.63738403295718D-09 )*Y-1.39237504892842D-08 )*Y+           &
           1.06527318142151D-07 )*Y-7.27634957230524D-07 )*Y+           &
           4.12159381310339D-06 )*Y-1.74648169719173D-05 )*Y+           &
           8.50290130067818D-05                                          
      WW3 = ((((((((((((-4.19569145459480D-17*Y+5.94344180261644D-16)*Y-&
           1.148797566469D-14)*Y+1.881303962576D-13)*Y-                 &
           2.413554618391D-12)*Y+3.372127423047D-11)*Y-                 &
           4.933988617784D-10)*Y+6.116545396281D-09)*Y-                 &
           6.69965691739299D-08 )*Y+7.52380085447161D-07 )*Y-           &
           8.08708393262321D-06 )*Y+6.88603417296672D-05 )*Y-           &
           4.67067112993427D-04 )*Y+5.42313365864597D-03                 
      WW2 = ((((((((((-6.22272689880615D-15*Y+1.04126809657554D-13)*Y-  &
           6.842418230913D-13)*Y+1.576841731919D-11)*Y-                 &
           4.203948834175D-10)*Y+6.287255934781D-09)*Y-                 &
           8.307159819228D-08)*Y+1.356478091922D-06)*Y-                 &
           2.08065576105639D-05 )*Y+2.52396730332340D-04 )*Y-           &
           2.94484050194539D-03 )*Y+6.01396183129168D-02                 
      WW1 = (((-1.8784686463512D-01/X+2.2991849164985D-01)/X -          &
           4.9893752514047D-01)/X-2.1916512131607D-05)*EXP(-X) +        &
           SQRT(PIE4/X)-WW4-WW3-WW2                                      
      RETURN
!                                                                        
  180 WW1 = SQRT(PIE4/X)                                                 
      IF (X > 35.0D+00) GO TO 220                                        
      IF (X > 20.0D+00) GO TO 200                                        
!     X=15.0 TO 20.0                             NROOTS = 4              
      Y = X-17.5D+00                                                     
      RT1 = ((((((((((( 4.36701759531398D-17*Y-1.12860600219889D-16)*Y- &
           6.149849164164D-15)*Y+5.820231579541D-14)*Y+                 &
           4.396602872143D-13)*Y-1.24330365320172D-11 )*Y+              &
           6.71083474044549D-11 )*Y+2.43865205376067D-10 )*Y+           &
           1.67559587099969D-08 )*Y-9.32738632357572D-07 )*Y+           &
           2.39030487004977D-05 )*Y-4.68648206591515D-04 )*Y+           &
           8.34977776583956D-03                                          
      RT2 = ((((((((((( 4.98913142288158D-16*Y-2.60732537093612D-16)*Y- &
           7.775156445127D-14)*Y+5.766105220086D-13)*Y+                 &
           6.432696729600D-12)*Y-1.39571683725792D-10 )*Y+              &
           5.95451479522191D-10 )*Y+2.42471442836205D-09 )*Y+           &
           2.47485710143120D-07 )*Y-1.14710398652091D-05 )*Y+           &
           2.71252453754519D-04 )*Y-4.96812745851408D-03 )*Y+           &
           8.26020602026780D-02                                          
      RT3 = ((((((((((( 1.91498302509009D-15*Y+1.48840394311115D-14)*Y- &
           4.316925145767D-13)*Y+1.186495793471D-12)*Y+                 &
           4.615806713055D-11)*Y-5.54336148667141D-10 )*Y+              &
           3.48789978951367D-10 )*Y-2.79188977451042D-09 )*Y+           &
           2.09563208958551D-06 )*Y-6.76512715080324D-05 )*Y+           &
           1.32129867629062D-03 )*Y-2.05062147771513D-02 )*Y+           &
           2.88068671894324D-01                                          
      RT4 = (((((((((((-5.43697691672942D-15*Y-1.12483395714468D-13)*Y+ &
           2.826607936174D-12)*Y-1.266734493280D-11)*Y-                 &
           4.258722866437D-10)*Y+9.45486578503261D-09 )*Y-              &
           5.86635622821309D-08 )*Y-1.28835028104639D-06 )*Y+           &
           4.41413815691885D-05 )*Y-7.61738385590776D-04 )*Y+           &
           9.66090902985550D-03 )*Y-1.01410568057649D-01 )*Y+           &
           9.54714798156712D-01                                          
      WW4 = ((((((((((((-7.56882223582704D-19*Y+7.53541779268175D-18)*Y-&
           1.157318032236D-16)*Y+2.411195002314D-15)*Y-                 &
           3.601794386996D-14)*Y+4.082150659615D-13)*Y-                 &
           4.289542980767D-12)*Y+5.086829642731D-11)*Y-                 &
           6.35435561050807D-10 )*Y+6.82309323251123D-09 )*Y-           &
           5.63374555753167D-08 )*Y+3.57005361100431D-07 )*Y-           &
           2.40050045173721D-06 )*Y+4.94171300536397D-05                 
      WW3 = (((((((((((-5.54451040921657D-17*Y+2.68748367250999D-16)*Y+ &
           1.349020069254D-14)*Y-2.507452792892D-13)*Y+                 &
           1.944339743818D-12)*Y-1.29816917658823D-11 )*Y+              &
           3.49977768819641D-10 )*Y-8.67270669346398D-09 )*Y+           &
           1.31381116840118D-07 )*Y-1.36790720600822D-06 )*Y+           &
           1.19210697673160D-05 )*Y-1.42181943986587D-04 )*Y+           &
           4.12615396191829D-03                                          
      WW2 = (((((((((((-1.86506057729700D-16*Y+1.16661114435809D-15)*Y+ &
           2.563712856363D-14)*Y-4.498350984631D-13)*Y+                 &
           1.765194089338D-12)*Y+9.04483676345625D-12 )*Y+              &
           4.98930345609785D-10 )*Y-2.11964170928181D-08 )*Y+           &
           3.98295476005614D-07 )*Y-5.49390160829409D-06 )*Y+           &
           7.74065155353262D-05 )*Y-1.48201933009105D-03 )*Y+           &
           4.97836392625268D-02                                          
      WW1 = (( 1.9623264149430D-01/X-4.9695241464490D-01)/X -           &
           6.0156581186481D-05)*EXP(-X)+WW1-WW2-WW3-WW4                  
      RETURN
!                                                                        
!     X=20.0 TO 35.0                             NROOTS = 4              
  200 E = EXP(-X)                                                        
      RT1 = ((((((-4.45711399441838D-05*X+1.27267770241379D-03)*X -     &
           2.36954961381262D-01)*X+1.54330657903756D+01)*X -            &
           5.22799159267808D+02)*X+1.05951216669313D+04)*X + (-         &
           2.51177235556236D+06/X+8.72975373557709D+05)/X -             &
           1.29194382386499D+05)*E + R14/(X-R14)                         
      RT2 = (((((-7.85617372254488D-02*X+6.35653573484868D+00)*X -      &
           3.38296938763990D+02)*X+1.25120495802096D+04)*X -            &
           3.16847570511637D+05)*X + ((-1.02427466127427D+09/X +        &
           3.70104713293016D+08)/X-5.87119005093822D+07)/X +            &
           5.38614211391604D+06)*E + R24/(X-R24)                         
      RT3 = (((((-2.37900485051067D-01*X+1.84122184400896D+01)*X -      &
           1.00200731304146D+03)*X+3.75151841595736D+04)*X -            &
           9.50626663390130D+05)*X + ((-2.88139014651985D+09/X +        &
           1.06625915044526D+09)/X-1.72465289687396D+08)/X +            &
           1.60419390230055D+07)*E + R34/(X-R34)                         
      RT4 = ((((((-6.00691586407385D-04*X-3.64479545338439D-01)*X +     &
           1.57496131755179D+01)*X-6.54944248734901D+02)*X +            &
           1.70830039597097D+04)*X-2.90517939780207D+05)*X + (+         &
           3.49059698304732D+07/X-1.64944522586065D+07)/X +             &
           2.96817940164703D+06)*E + R44/(X-R44)                         
      IF (X <= 25.0D+00) WW4 = ((((((( 2.33766206773151D-07*X-          &
           3.81542906607063D-05)*X +3.51416601267000D-03)*X-            &
           1.66538571864728D-01)*X +4.80006136831847D+00)*X-            &
           8.73165934223603D+01)*X +9.77683627474638D+02)*X +           &
           1.66000945117640D+04/X -6.14479071209961D+03)*E + W44*WW1     
      IF (X > 25.0D+00) WW4 = (((((( 5.74245945342286D-06*X-            &
           7.58735928102351D-05)*X +2.35072857922892D-04)*X-            &
           3.78812134013125D-03)*X +3.09871652785805D-01)*X-            &
           7.11108633061306D+00)*X +5.55297573149528D+01)*E + W44*WW1    
      WW3 = (((((( 2.36392855180768D-04*X-9.16785337967013D-03)*X +     &
           4.62186525041313D-01)*X-1.96943786006540D+01)*X +            &
           4.99169195295559D+02)*X-6.21419845845090D+03)*X + ((+        &
           5.21445053212414D+07/X-1.34113464389309D+07)/X +             &
           1.13673298305631D+06)/X-2.81501182042707D+03)*E + W34*WW1     
      WW2 = (((((( 7.29841848989391D-04*X-3.53899555749875D-02)*X +     &
           2.07797425718513D+00)*X-1.00464709786287D+02)*X +            &
           3.15206108877819D+03)*X-6.27054715090012D+04)*X + (+         &
           1.54721246264919D+07/X-5.26074391316381D+06)/X +             &
           7.67135400969617D+05)*E + W24*WW1                             
      WW1 = (( 1.9623264149430D-01/X-4.9695241464490D-01)/X -           &
           6.0156581186481D-05)*E + WW1-WW2-WW3-WW4                      
      RETURN
!                                                                        
  220 IF (X > 53.0D+00) GO TO 240                                        
!     X=35.0 TO 53.0                             NROOTS = 4              
      E = EXP(-X)*(X*X)**2                                               
      RT4 = ((-2.19135070169653D-03*X-1.19108256987623D-01)*X -         &
           7.50238795695573D-01)*E + R44/(X-R44)                         
      RT3 = ((-9.65842534508637D-04*X-4.49822013469279D-02)*X +         &
           6.08784033347757D-01)*E + R34/(X-R34)                         
      RT2 = ((-3.62569791162153D-04*X-9.09231717268466D-03)*X +         &
           1.84336760556262D-01)*E + R24/(X-R24)                         
      RT1 = ((-4.07557525914600D-05*X-6.88846864931685D-04)*X +         &
           1.74725309199384D-02)*E + R14/(X-R14)                         
      WW4 = (( 5.76631982000990D-06*X-7.89187283804890D-05)*X +         &
           3.28297971853126D-04)*E + W44*WW1                             
      WW3 = (( 2.08294969857230D-04*X-3.77489954837361D-03)*X +         &
           2.09857151617436D-02)*E + W34*WW1                             
      WW2 = (( 6.16374517326469D-04*X-1.26711744680092D-02)*X +         &
          8.14504890732155D-02)*E + W24*WW1
      WW1 = WW1-WW2-WW3-WW4
      RETURN
!
!     X=47.0 TO INFINITY                         NROOTS = 4
  240 RT1 = R14/(X-R14)
      RT2 = R24/(X-R24)
      RT3 = R34/(X-R34)
      RT4 = R44/(X-R44)
      WW4 = W44*WW1
      WW3 = W34*WW1
      WW2 = W24*WW1
      WW1 = WW1-WW2-WW3-WW4
      RETURN
      END

! Calling by HSandT & ERISPDFGHIL: ROOT5 
      SUBROUTINE ROOT5
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /ROOT  / X,U(13),W(13),NROOTS
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
      DATA R15,PIE4/1.17581320211778D-01, 7.85398163397448D-01/
      DATA R25,W25/ 1.07456201243690D+00, 2.70967405960535D-01/
      DATA R35,W35/ 3.08593744371754D+00, 3.82231610015404D-02/
      DATA R45,W45/ 6.41472973366203D+00, 1.51614186862443D-03/
      DATA R55,W55/ 1.18071894899717D+01, 8.62130526143657D-06/
!
      IF (X > 15.0D+00) GO TO 180
      IF (X > 5.0D+00) GO TO 140
      IF (X > 1.0D+00) GO TO 120
      IF (X > 3.0D-07) GO TO 100
!     X IS APPROXIMATELY 0.0d0.                   NROOTS = 5
      RT1 = 2.26659266316985D-02 -2.15865967920897D-03 *X
      RT2 = 2.31271692140903D-01 -2.20258754389745D-02 *X
      RT3 = 8.57346024118836D-01 -8.16520023025515D-02 *X
      RT4 = 2.97353038120346D+00 -2.83193369647137D-01 *X
      RT5 = 1.84151859759051D+01 -1.75382723579439D+00 *X
      WW1 = 2.95524224714752D-01 -1.96867576909777D-02 *X
      WW2 = 2.69266719309995D-01 -5.61737590184721D-02 *X
      WW3 = 2.19086362515981D-01 -9.71152726793658D-02 *X
      WW4 = 1.49451349150580D-01 -1.02979262193565D-01 *X
      WW5 = 6.66713443086877D-02 -5.73782817488315D-02 *X
      RETURN
!
!     X=0.0 TO 1.0                               NROOTS = 5
  100 RT1 = ((((((-4.46679165328413D-11*X+1.21879111988031D-09)*X-      &
           2.62975022612104D-08 )*X+5.15106194905897D-07 )*X-           &
           9.27933625824749D-06 )*X+1.51794097682482D-04 )*X-           &
           2.15865967920301D-03 )*X+2.26659266316985D-02                 
      RT2 = (((((( 1.93117331714174D-10*X-4.57267589660699D-09)*X+      &
           2.48339908218932D-08 )*X+1.50716729438474D-06 )*X-           &
           6.07268757707381D-05 )*X+1.37506939145643D-03 )*X-           &
           2.20258754419939D-02 )*X+2.31271692140905D-01                 
      RT3 = ((((( 4.84989776180094D-09*X+1.31538893944284D-07)*X-       &
           2.766753852879D-06)*X-7.651163510626D-05)*X+                 &
           4.033058545972D-03)*X-8.16520022916145D-02 )*X+              &
           8.57346024118779D-01                                          
      RT4 = ((((-2.48581772214623D-07*X-4.34482635782585D-06)*X-        &
           7.46018257987630D-07 )*X+1.01210776517279D-02 )*X-           &
           2.83193369640005D-01 )*X+2.97353038120345D+00                 
      RT5 = (((((-8.92432153868554D-09*X+1.77288899268988D-08)*X+       &
           3.040754680666D-06)*X+1.058229325071D-04)*X+                 &
           4.596379534985D-02)*X-1.75382723579114D+00 )*X+              &
           1.84151859759049D+01                                          
      WW1 = ((((((-2.03822632771791D-09*X+3.89110229133810D-08)*X-      &
           5.84914787904823D-07 )*X+8.30316168666696D-06 )*X-           &
           1.13218402310546D-04 )*X+1.49128888586790D-03 )*X-           &
           1.96867576904816D-02 )*X+2.95524224714749D-01                 
      WW2 = ((((((( 8.62848118397570D-09*X-1.38975551148989D-07)*X+     &
           1.602894068228D-06)*X-1.646364300836D-05)*X+                 &
           1.538445806778D-04)*X-1.28848868034502D-03 )*X+              &
           9.38866933338584D-03 )*X-5.61737590178812D-02 )*X+           &
           2.69266719309991D-01                                          
      WW3 = ((((((((-9.41953204205665D-09*X+1.47452251067755D-07)*X-    &
           1.57456991199322D-06 )*X+1.45098401798393D-05 )*X-           &
           1.18858834181513D-04 )*X+8.53697675984210D-04 )*X-           &
           5.22877807397165D-03 )*X+2.60854524809786D-02 )*X-           &
           9.71152726809059D-02 )*X+2.19086362515979D-01                 
      WW4 = ((((((((-3.84961617022042D-08*X+5.66595396544470D-07)*X-    &
           5.52351805403748D-06 )*X+4.53160377546073D-05 )*X-           &
           3.22542784865557D-04 )*X+1.95682017370967D-03 )*X-           &
           9.77232537679229D-03 )*X+3.79455945268632D-02 )*X-           &
           1.02979262192227D-01 )*X+1.49451349150573D-01                 
      WW5 = ((((((((( 4.09594812521430D-09*X-6.47097874264417D-08)*X+   &
           6.743541482689D-07)*X-5.917993920224D-06)*X+                 &
           4.531969237381D-05)*X-2.99102856679638D-04 )*X+              &
           1.65695765202643D-03 )*X-7.40671222520653D-03 )*X+           &
           2.50889946832192D-02 )*X-5.73782817487958D-02 )*X+           &
           6.66713443086877D-02                                          
      RETURN
                                                                        
!     X=1.0 TO 5.0                               NROOTS = 5              
  120 Y = X-3.0D+00                                                      
      RT1 = ((((((((-2.58163897135138D-14*Y+8.14127461488273D-13)*Y-    &
           2.11414838976129D-11 )*Y+5.09822003260014D-10 )*Y-           &
           1.16002134438663D-08 )*Y+2.46810694414540D-07 )*Y-           &
           4.92556826124502D-06 )*Y+9.02580687971053D-05 )*Y-           &
           1.45190025120726D-03 )*Y+1.73416786387475D-02                 
      RT2 = ((((((((( 1.04525287289788D-14*Y+5.44611782010773D-14)*Y-   &
           4.831059411392D-12)*Y+1.136643908832D-10)*Y-                 &
           1.104373076913D-09)*Y-2.35346740649916D-08 )*Y+              &
           1.43772622028764D-06 )*Y-4.23405023015273D-05 )*Y+           &
           9.12034574793379D-04 )*Y-1.52479441718739D-02 )*Y+           &
           1.76055265928744D-01                                          
      RT3 = (((((((((-6.89693150857911D-14*Y+5.92064260918861D-13)*Y+   &
           1.847170956043D-11)*Y-3.390752744265D-10)*Y-                 &
           2.995532064116D-09)*Y+1.57456141058535D-07 )*Y-              &
           3.95859409711346D-07 )*Y-9.58924580919747D-05 )*Y+           &
           3.23551502557785D-03 )*Y-5.97587007636479D-02 )*Y+           &
           6.46432853383057D-01                                          
      RT4 = ((((((((-3.61293809667763D-12*Y-2.70803518291085D-11)*Y+    &
           8.83758848468769D-10 )*Y+1.59166632851267D-08 )*Y-           &
           1.32581997983422D-07 )*Y-7.60223407443995D-06 )*Y-           &
           7.41019244900952D-05 )*Y+9.81432631743423D-03 )*Y-           &
           2.23055570487771D-01 )*Y+2.21460798080643D+00                 
      RT5 = ((((((((( 7.12332088345321D-13*Y+3.16578501501894D-12)*Y-   &
           8.776668218053D-11)*Y-2.342817613343D-09)*Y-                 &
           3.496962018025D-08)*Y-3.03172870136802D-07 )*Y+              &
           1.50511293969805D-06 )*Y+1.37704919387696D-04 )*Y+           &
           4.70723869619745D-02 )*Y-1.47486623003693D+00 )*Y+           &
           1.35704792175847D+01                                          
      WW1 = ((((((((( 1.04348658616398D-13*Y-1.94147461891055D-12)*Y+   &
           3.485512360993D-11)*Y-6.277497362235D-10)*Y+                 &
           1.100758247388D-08)*Y-1.88329804969573D-07 )*Y+              &
           3.12338120839468D-06 )*Y-5.04404167403568D-05 )*Y+           &
           8.00338056610995D-04 )*Y-1.30892406559521D-02 )*Y+           &
           2.47383140241103D-01                                          
      WW2 = ((((((((((( 3.23496149760478D-14*Y-5.24314473469311D-13)*Y+ &
           7.743219385056D-12)*Y-1.146022750992D-10)*Y+                 &
           1.615238462197D-09)*Y-2.15479017572233D-08 )*Y+              &
           2.70933462557631D-07 )*Y-3.18750295288531D-06 )*Y+           &
           3.47425221210099D-05 )*Y-3.45558237388223D-04 )*Y+           &
           3.05779768191621D-03 )*Y-2.29118251223003D-02 )*Y+           &
           1.59834227924213D-01                                          
      WW3 = ((((((((((((-3.42790561802876D-14*Y+5.26475736681542D-13)*Y-&
           7.184330797139D-12)*Y+9.763932908544D-11)*Y-                 &
           1.244014559219D-09)*Y+1.472744068942D-08)*Y-                 &
           1.611749975234D-07)*Y+1.616487851917D-06)*Y-                 &
           1.46852359124154D-05 )*Y+1.18900349101069D-04 )*Y-           &
           8.37562373221756D-04 )*Y+4.93752683045845D-03 )*Y-           &
           2.25514728915673D-02 )*Y+6.95211812453929D-02                 
      WW4 = ((((((((((((( 1.04072340345039D-14*Y-1.60808044529211D-13)* &
           Y+2.183534866798D-12)*Y-2.939403008391D-11)*Y+               &
           3.679254029085D-10)*Y-4.23775673047899D-09 )*Y+              &
           4.46559231067006D-08 )*Y-4.26488836563267D-07 )*Y+           &
           3.64721335274973D-06 )*Y-2.74868382777722D-05 )*Y+           &
           1.78586118867488D-04 )*Y-9.68428981886534D-04 )*Y+           &
           4.16002324339929D-03 )*Y-1.28290192663141D-02 )*Y+           &
           2.22353727685016D-02                                          
      WW5 = ((((((((((((((-8.16770412525963D-16*Y+1.31376515047977D-14)*&
           Y-1.856950818865D-13)*Y+2.596836515749D-12)*Y-               &
           3.372639523006D-11)*Y+4.025371849467D-10)*Y-                 &
           4.389453269417D-09)*Y+4.332753856271D-08)*Y-                 &
           3.82673275931962D-07 )*Y+2.98006900751543D-06 )*Y-           &
           2.00718990300052D-05 )*Y+1.13876001386361D-04 )*Y-           &
           5.23627942443563D-04 )*Y+1.83524565118203D-03 )*Y-           &
           4.37785737450783D-03 )*Y+5.36963805223095D-03                 
      RETURN
                                                                        
  140 IF (X > 10.0D+00) GO TO 160                                        
!     X=5.0 TO 10.0                              NROOTS = 5              
      Y = X-7.5D+00                                                      
      RT1 = ((((((((-1.13825201010775D-14*Y+1.89737681670375D-13)*Y-    &
           4.81561201185876D-12 )*Y+1.56666512163407D-10 )*Y-           &
           3.73782213255083D-09 )*Y+9.15858355075147D-08 )*Y-           &
           2.13775073585629D-06 )*Y+4.56547356365536D-05 )*Y-           &
           8.68003909323740D-04 )*Y+1.22703754069176D-02                 
      RT2 = (((((((((-3.67160504428358D-15*Y+1.27876280158297D-14)*Y-   &
           1.296476623788D-12)*Y+1.477175434354D-11)*Y+                 &
           5.464102147892D-10)*Y-2.42538340602723D-08 )*Y+              &
           8.20460740637617D-07 )*Y-2.20379304598661D-05 )*Y+           &
           4.90295372978785D-04 )*Y-9.14294111576119D-03 )*Y+           &
           1.22590403403690D-01                                          
      RT3 = ((((((((( 1.39017367502123D-14*Y-6.96391385426890D-13)*Y+   &
           1.176946020731D-12)*Y+1.725627235645D-10)*Y-                 &
           3.686383856300D-09)*Y+2.87495324207095D-08 )*Y+              &
           1.71307311000282D-06 )*Y-7.94273603184629D-05 )*Y+           &
           2.00938064965897D-03 )*Y-3.63329491677178D-02 )*Y+           &
           4.34393683888443D-01                                          
      RT4 = ((((((((((-1.27815158195209D-14*Y+1.99910415869821D-14)*Y+  &
           3.753542914426D-12)*Y-2.708018219579D-11)*Y-                 &
           1.190574776587D-09)*Y+1.106696436509D-08)*Y+                 &
           3.954955671326D-07)*Y-4.398596059588D-06)*Y-                 &
           2.01087998907735D-04 )*Y+7.89092425542937D-03 )*Y-           &
           1.42056749162695D-01 )*Y+1.39964149420683D+00                 
      RT5 = ((((((((((-1.19442341030461D-13*Y-2.34074833275956D-12)*Y+  &
           6.861649627426D-12)*Y+6.082671496226D-10)*Y+                 &
           5.381160105420D-09)*Y-6.253297138700D-08)*Y-                 &
           2.135966835050D-06)*Y-2.373394341886D-05)*Y+                 &
           2.88711171412814D-06 )*Y+4.85221195290753D-02 )*Y-           &
           1.04346091985269D+00 )*Y+7.89901551676692D+00                 
      WW1 = ((((((((( 7.95526040108997D-15*Y-2.48593096128045D-13)*Y+   &
           4.761246208720D-12)*Y-9.535763686605D-11)*Y+                 &
           2.225273630974D-09)*Y-4.49796778054865D-08 )*Y+              &
           9.17812870287386D-07 )*Y-1.86764236490502D-05 )*Y+           &
           3.76807779068053D-04 )*Y-8.10456360143408D-03 )*Y+           &
           2.01097936411496D-01                                          
      WW2 = ((((((((((( 1.25678686624734D-15*Y-2.34266248891173D-14)*Y+ &
           3.973252415832D-13)*Y-6.830539401049D-12)*Y+                 &
           1.140771033372D-10)*Y-1.82546185762009D-09 )*Y+              &
           2.77209637550134D-08 )*Y-4.01726946190383D-07 )*Y+           &
           5.48227244014763D-06 )*Y-6.95676245982121D-05 )*Y+           &
           8.05193921815776D-04 )*Y-8.15528438784469D-03 )*Y+           &
           9.71769901268114D-02                                          
      WW3 = ((((((((((((-8.20929494859896D-16*Y+1.37356038393016D-14)*Y-&
           2.022863065220D-13)*Y+3.058055403795D-12)*Y-                 &
           4.387890955243D-11)*Y+5.923946274445D-10)*Y-                 &
           7.503659964159D-09)*Y+8.851599803902D-08)*Y-                 &
           9.65561998415038D-07 )*Y+9.60884622778092D-06 )*Y-           &
           8.56551787594404D-05 )*Y+6.66057194311179D-04 )*Y-           &
           4.17753183902198D-03 )*Y+2.25443826852447D-02                 
      WW4 = ((((((((((((((-1.08764612488790D-17*Y+1.85299909689937D-16)*&
           Y-2.730195628655D-15)*Y+4.127368817265D-14)*Y-               &
           5.881379088074D-13)*Y+7.805245193391D-12)*Y-                 &
           9.632707991704D-11)*Y+1.099047050624D-09)*Y-                 &
           1.15042731790748D-08 )*Y+1.09415155268932D-07 )*Y-           &
           9.33687124875935D-07 )*Y+7.02338477986218D-06 )*Y-           &
           4.53759748787756D-05 )*Y+2.41722511389146D-04 )*Y-           &
           9.75935943447037D-04 )*Y+2.57520532789644D-03                 
      WW5 = ((((((((((((((( 7.28996979748849D-19*Y-1.26518146195173D-17)&
           *Y+1.886145834486D-16)*Y-2.876728287383D-15)*Y+              &
           4.114588668138D-14)*Y-5.44436631413933D-13 )*Y+              &
           6.64976446790959D-12 )*Y-7.44560069974940D-11 )*Y+           &
           7.57553198166848D-10 )*Y-6.92956101109829D-09 )*Y+           &
           5.62222859033624D-08 )*Y-3.97500114084351D-07 )*Y+           &
           2.39039126138140D-06 )*Y-1.18023950002105D-05 )*Y+           &
           4.52254031046244D-05 )*Y-1.21113782150370D-04 )*Y+           &
           1.75013126731224D-04                                          
      RETURN
                                                                        
!     X=10.0 TO 15.0                             NROOTS = 5              
  160 Y = X-12.5D+00                                                     
      RT1 = ((((((((((-4.16387977337393D-17*Y+7.20872997373860D-16)*Y+  &
           1.395993802064D-14)*Y+3.660484641252D-14)*Y-                 &
           4.154857548139D-12)*Y+2.301379846544D-11)*Y-                 &
           1.033307012866D-09)*Y+3.997777641049D-08)*Y-                 &
           9.35118186333939D-07 )*Y+2.38589932752937D-05 )*Y-           &
           5.35185183652937D-04 )*Y+8.85218988709735D-03                 
      RT2 = ((((((((((-4.56279214732217D-16*Y+6.24941647247927D-15)*Y+  &
           1.737896339191D-13)*Y+8.964205979517D-14)*Y-                 &
           3.538906780633D-11)*Y+9.561341254948D-11)*Y-                 &
           9.772831891310D-09)*Y+4.240340194620D-07)*Y-                 &
           1.02384302866534D-05 )*Y+2.57987709704822D-04 )*Y-           &
           5.54735977651677D-03 )*Y+8.68245143991948D-02                 
      RT3 = ((((((((((-2.52879337929239D-15*Y+2.13925810087833D-14)*Y+  &
           7.884307667104D-13)*Y-9.023398159510D-13)*Y-                 &
           5.814101544957D-11)*Y-1.333480437968D-09)*Y-                 &
           2.217064940373D-08)*Y+1.643290788086D-06)*Y-                 &
           4.39602147345028D-05 )*Y+1.08648982748911D-03 )*Y-           &
           2.13014521653498D-02 )*Y+2.94150684465425D-01                 
      RT4 = ((((((((((-6.42391438038888D-15*Y+5.37848223438815D-15)*Y+  &
           8.960828117859D-13)*Y+5.214153461337D-11)*Y-                 &
           1.106601744067D-10)*Y-2.007890743962D-08)*Y+                 &
           1.543764346501D-07)*Y+4.520749076914D-06)*Y-                 &
           1.88893338587047D-04 )*Y+4.73264487389288D-03 )*Y-           &
           7.91197893350253D-02 )*Y+8.60057928514554D-01                 
      RT5 = (((((((((((-2.24366166957225D-14*Y+4.87224967526081D-14)*Y+ &
           5.587369053655D-12)*Y-3.045253104617D-12)*Y-                 &
           1.223983883080D-09)*Y-2.05603889396319D-09 )*Y+              &
           2.58604071603561D-07 )*Y+1.34240904266268D-06 )*Y-           &
           5.72877569731162D-05 )*Y-9.56275105032191D-04 )*Y+           &
           4.23367010370921D-02 )*Y-5.76800927133412D-01 )*Y+           &
           3.87328263873381D+00                                          
      WW1 = ((((((((( 8.98007931950169D-15*Y+7.25673623859497D-14)*Y+   &
           5.851494250405D-14)*Y-4.234204823846D-11)*Y+                 &
           3.911507312679D-10)*Y-9.65094802088511D-09 )*Y+              &
           3.42197444235714D-07 )*Y-7.51821178144509D-06 )*Y+           &
           1.94218051498662D-04 )*Y-5.38533819142287D-03 )*Y+           &
           1.68122596736809D-01                                          
      WW2 = ((((((((((-1.05490525395105D-15*Y+1.96855386549388D-14)*Y-  &
           5.500330153548D-13)*Y+1.003849567976D-11)*Y-                 &
           1.720997242621D-10)*Y+3.533277061402D-09)*Y-                 &
           6.389171736029D-08)*Y+1.046236652393D-06)*Y-                 &
           1.73148206795827D-05 )*Y+2.57820531617185D-04 )*Y-           &
           3.46188265338350D-03 )*Y+7.03302497508176D-02                 
      WW3 = ((((((((((( 3.60020423754545D-16*Y-6.24245825017148D-15)*Y+ &
           9.945311467434D-14)*Y-1.749051512721D-12)*Y+                 &
           2.768503957853D-11)*Y-4.08688551136506D-10 )*Y+              &
           6.04189063303610D-09 )*Y-8.23540111024147D-08 )*Y+           &
           1.01503783870262D-06 )*Y-1.20490761741576D-05 )*Y+           &
           1.26928442448148D-04 )*Y-1.05539461930597D-03 )*Y+           &
           1.15543698537013D-02                                          
      WW4 = ((((((((((((( 2.51163533058925D-18*Y-4.31723745510697D-17)* &
           Y+6.557620865832D-16)*Y-1.016528519495D-14)*Y+               &
           1.491302084832D-13)*Y-2.06638666222265D-12 )*Y+              &
           2.67958697789258D-11 )*Y-3.23322654638336D-10 )*Y+           &
           3.63722952167779D-09 )*Y-3.75484943783021D-08 )*Y+           &
           3.49164261987184D-07 )*Y-2.92658670674908D-06 )*Y+           &
           2.12937256719543D-05 )*Y-1.19434130620929D-04 )*Y+           &
           6.45524336158384D-04                                          
      WW5 = ((((((((((((((-1.29043630202811D-19*Y+2.16234952241296D-18)*&
           Y-3.107631557965D-17)*Y+4.570804313173D-16)*Y-               &
           6.301348858104D-15)*Y+8.031304476153D-14)*Y-                 &
           9.446196472547D-13)*Y+1.018245804339D-11)*Y-                 &
           9.96995451348129D-11 )*Y+8.77489010276305D-10 )*Y-           &
           6.84655877575364D-09 )*Y+4.64460857084983D-08 )*Y-           &
           2.66924538268397D-07 )*Y+1.24621276265907D-06 )*Y-           &
           4.30868944351523D-06 )*Y+9.94307982432868D-06                 
      RETURN
                                                                        
  180 IF (X > 25.0D+00) GO TO 220                                        
      IF (X > 20.0D+00) GO TO 200                                        
!     X=15.0 TO 20.0                             NROOTS = 5              
      Y = X-17.5D+00                                                     
      RT1 = (((((((((( 1.91875764545740D-16*Y+7.8357401095707D-16)*Y-   &
           3.260875931644D-14)*Y-1.186752035569D-13)*Y+                 &
           4.275180095653D-12)*Y+3.357056136731D-11)*Y-                 &
           1.123776903884D-09)*Y+1.231203269887D-08)*Y-                 &
           3.99851421361031D-07 )*Y+1.45418822817771D-05 )*Y-           &
           3.49912254976317D-04 )*Y+6.67768703938812D-03                 
      RT2 = (((((((((( 2.02778478673555D-15*Y+1.01640716785099D-14)*Y-  &
           3.385363492036D-13)*Y-1.615655871159D-12)*Y+                 &
           4.527419140333D-11)*Y+3.853670706486D-10)*Y-                 &
           1.184607130107D-08)*Y+1.347873288827D-07)*Y-                 &
           4.47788241748377D-06 )*Y+1.54942754358273D-04 )*Y-           &
           3.55524254280266D-03 )*Y+6.44912219301603D-02                 
      RT3 = (((((((((( 7.79850771456444D-15*Y+6.00464406395001D-14)*Y-  &
           1.249779730869D-12)*Y-1.020720636353D-11)*Y+                 &
           1.814709816693D-10)*Y+1.766397336977D-09)*Y-                 &
           4.603559449010D-08)*Y+5.863956443581D-07)*Y-                 &
           2.03797212506691D-05 )*Y+6.31405161185185D-04 )*Y-           &
           1.30102750145071D-02 )*Y+2.10244289044705D-01                 
      RT4 = (((((((((((-2.92397030777912D-15*Y+1.94152129078465D-14)*Y+ &
           4.859447665850D-13)*Y-3.217227223463D-12)*Y-                 &
           7.484522135512D-11)*Y+7.19101516047753D-10 )*Y+              &
           6.88409355245582D-09 )*Y-1.44374545515769D-07 )*Y+           &
           2.74941013315834D-06 )*Y-1.02790452049013D-04 )*Y+           &
           2.59924221372643D-03 )*Y-4.35712368303551D-02 )*Y+           &
           5.62170709585029D-01                                          
      RT5 = ((((((((((( 1.17976126840060D-14*Y+1.24156229350669D-13)*Y- &
           3.892741622280D-12)*Y-7.755793199043D-12)*Y+                 &
           9.492190032313D-10)*Y-4.98680128123353D-09 )*Y-              &
           1.81502268782664D-07 )*Y+2.69463269394888D-06 )*Y+           &
           2.50032154421640D-05 )*Y-1.33684303917681D-03 )*Y+           &
           2.29121951862538D-02 )*Y-2.45653725061323D-01 )*Y+           &
           1.89999883453047D+00                                          
      WW1 = (((((((((( 1.74841995087592D-15*Y-6.95671892641256D-16)*Y-  &
           3.000659497257D-13)*Y+2.021279817961D-13)*Y+                 &
           3.853596935400D-11)*Y+1.461418533652D-10)*Y-                 &
           1.014517563435D-08)*Y+1.132736008979D-07)*Y-                 &
           2.86605475073259D-06 )*Y+1.21958354908768D-04 )*Y-           &
           3.86293751153466D-03 )*Y+1.45298342081522D-01                 
      WW2 = ((((((((((-1.11199320525573D-15*Y+1.85007587796671D-15)*Y+  &
           1.220613939709D-13)*Y+1.275068098526D-12)*Y-                 &
           5.341838883262D-11)*Y+6.161037256669D-10)*Y-                 &
           1.009147879750D-08)*Y+2.907862965346D-07)*Y-                 &
           6.12300038720919D-06 )*Y+1.00104454489518D-04 )*Y-           &
           1.80677298502757D-03 )*Y+5.78009914536630D-02                 
      WW3 = ((((((((((-9.49816486853687D-16*Y+6.67922080354234D-15)*Y+  &
           2.606163540537D-15)*Y+1.983799950150D-12)*Y-                 &
           5.400548574357D-11)*Y+6.638043374114D-10)*Y-                 &
           8.799518866802D-09)*Y+1.791418482685D-07)*Y-                 &
           2.96075397351101D-06 )*Y+3.38028206156144D-05 )*Y-           &
           3.58426847857878D-04 )*Y+8.39213709428516D-03                 
      WW4 = ((((((((((( 1.33829971060180D-17*Y-3.44841877844140D-16)*Y+ &
           4.745009557656D-15)*Y-6.033814209875D-14)*Y+                 &
           1.049256040808D-12)*Y-1.70859789556117D-11 )*Y+              &
           2.15219425727959D-10 )*Y-2.52746574206884D-09 )*Y+           &
           3.27761714422960D-08 )*Y-3.90387662925193D-07 )*Y+           &
           3.46340204593870D-06 )*Y-2.43236345136782D-05 )*Y+           &
           3.54846978585226D-04                                          
      WW5 = ((((((((((((( 2.69412277020887D-20*Y-4.24837886165685D-19)* &
           Y+6.030500065438D-18)*Y-9.069722758289D-17)*Y+               &
           1.246599177672D-15)*Y-1.56872999797549D-14 )*Y+              &
           1.87305099552692D-13 )*Y-2.09498886675861D-12 )*Y+           &
           2.11630022068394D-11 )*Y-1.92566242323525D-10 )*Y+           &
           1.62012436344069D-09 )*Y-1.23621614171556D-08 )*Y+           &
           7.72165684563049D-08 )*Y-3.59858901591047D-07 )*Y+           &
           2.43682618601000D-06                                          
      RETURN
                                                                        
!     X=20.0 TO 25.0                             NROOTS = 5              
  200 Y = X-22.5D+00                                                     
      RT1 = (((((((((-1.13927848238726D-15*Y+7.39404133595713D-15)*Y+   &
           1.445982921243D-13)*Y-2.676703245252D-12)*Y+                 &
           5.823521627177D-12)*Y+2.17264723874381D-10 )*Y+              &
           3.56242145897468D-09 )*Y-3.03763737404491D-07 )*Y+           &
           9.46859114120901D-06 )*Y-2.30896753853196D-04 )*Y+           &
           5.24663913001114D-03                                          
      RT2 = (((((((((( 2.89872355524581D-16*Y-1.22296292045864D-14)*Y+  &
           6.184065097200D-14)*Y+1.649846591230D-12)*Y-                 &
           2.729713905266D-11)*Y+3.709913790650D-11)*Y+                 &
           2.216486288382D-09)*Y+4.616160236414D-08)*Y-                 &
           3.32380270861364D-06 )*Y+9.84635072633776D-05 )*Y-           &
           2.30092118015697D-03 )*Y+5.00845183695073D-02                 
      RT3 = (((((((((( 1.97068646590923D-15*Y-4.89419270626800D-14)*Y+  &
           1.136466605916D-13)*Y+7.546203883874D-12)*Y-                 &
           9.635646767455D-11)*Y-8.295965491209D-11)*Y+                 &
           7.534109114453D-09)*Y+2.699970652707D-07)*Y-                 &
           1.42982334217081D-05 )*Y+3.78290946669264D-04 )*Y-           &
           8.03133015084373D-03 )*Y+1.58689469640791D-01                 
      RT4 = (((((((((( 1.33642069941389D-14*Y-1.55850612605745D-13)*Y-  &
           7.522712577474D-13)*Y+3.209520801187D-11)*Y-                 &
           2.075594313618D-10)*Y-2.070575894402D-09)*Y+                 &
           7.323046997451D-09)*Y+1.851491550417D-06)*Y-                 &
           6.37524802411383D-05 )*Y+1.36795464918785D-03 )*Y-           &
           2.42051126993146D-02 )*Y+3.97847167557815D-01                 
      RT5 = ((((((((((-6.07053986130526D-14*Y+1.04447493138843D-12)*Y-  &
           4.286617818951D-13)*Y-2.632066100073D-10)*Y+                 &
           4.804518986559D-09)*Y-1.835675889421D-08)*Y-                 &
           1.068175391334D-06)*Y+3.292234974141D-05)*Y-                 &
           5.94805357558251D-04 )*Y+8.29382168612791D-03 )*Y-           &
           9.93122509049447D-02 )*Y+1.09857804755042D+00                 
      WW1 = (((((((((-9.10338640266542D-15*Y+1.00438927627833D-13)*Y+   &
           7.817349237071D-13)*Y-2.547619474232D-11)*Y+                 &
           1.479321506529D-10)*Y+1.52314028857627D-09 )*Y+              &
           9.20072040917242D-09 )*Y-2.19427111221848D-06 )*Y+           &
           8.65797782880311D-05 )*Y-2.82718629312875D-03 )*Y+           &
           1.28718310443295D-01                                          
      WW2 = ((((((((( 5.52380927618760D-15*Y-6.43424400204124D-14)*Y-   &
           2.358734508092D-13)*Y+8.261326648131D-12)*Y+                 &
           9.229645304956D-11)*Y-5.68108973828949D-09 )*Y+              &
           1.22477891136278D-07 )*Y-2.11919643127927D-06 )*Y+           &
           4.23605032368922D-05 )*Y-1.14423444576221D-03 )*Y+           &
           5.06607252890186D-02                                          
      WW3 = ((((((((( 3.99457454087556D-15*Y-5.11826702824182D-14)*Y-   &
           4.157593182747D-14)*Y+4.214670817758D-12)*Y+                 &
           6.705582751532D-11)*Y-3.36086411698418D-09 )*Y+              &
           6.07453633298986D-08 )*Y-7.40736211041247D-07 )*Y+           &
           8.84176371665149D-06 )*Y-1.72559275066834D-04 )*Y+           &
           7.16639814253567D-03                                          
      WW4 = (((((((((((-2.14649508112234D-18*Y-2.45525846412281D-18)*Y+ &
           6.126212599772D-16)*Y-8.526651626939D-15)*Y+                 &
           4.826636065733D-14)*Y-3.39554163649740D-13 )*Y+              &
           1.67070784862985D-11 )*Y-4.42671979311163D-10 )*Y+           &
           6.77368055908400D-09 )*Y-7.03520999708859D-08 )*Y+           &
           6.04993294708874D-07 )*Y-7.80555094280483D-06 )*Y+           &
           2.85954806605017D-04                                          
      WW5 = ((((((((((((-5.63938733073804D-21*Y+6.92182516324628D-20)*Y-&
           1.586937691507D-18)*Y+3.357639744582D-17)*Y-                 &
           4.810285046442D-16)*Y+5.386312669975D-15)*Y-                 &
           6.117895297439D-14)*Y+8.441808227634D-13)*Y-                 &
           1.18527596836592D-11 )*Y+1.36296870441445D-10 )*Y-           &
           1.17842611094141D-09 )*Y+7.80430641995926D-09 )*Y-           &
           5.97767417400540D-08 )*Y+1.65186146094969D-06                 
      RETURN
                                                                        
  220 WW1 = SQRT(PIE4/X)                                                 
      IF (X > 40.0D+00) GO TO 240                                        
!     X=25.0 TO 40.0                             NROOTS = 5              
      E = EXP(-X)                                                        
      RT1 = ((((((((-1.73363958895356D-06*X+1.19921331441483D-04)*X -   &
           1.59437614121125D-02)*X+1.13467897349442D+00)*X -            &
           4.47216460864586D+01)*X+1.06251216612604D+03)*X -            &
           1.52073917378512D+04)*X+1.20662887111273D+05)*X -            &
           4.07186366852475D+05)*E + R15/(X-R15)                         
      RT2 = ((((((((-1.60102542621710D-05*X+1.10331262112395D-03)*X -   &
           1.50043662589017D-01)*X+1.05563640866077D+01)*X -            &
           4.10468817024806D+02)*X+9.62604416506819D+03)*X -            &
           1.35888069838270D+05)*X+1.06107577038340D+06)*X -            &
           3.51190792816119D+06)*E + R25/(X-R25)                         
      RT3 = ((((((((-4.48880032128422D-05*X+2.69025112122177D-03)*X -   &
           4.01048115525954D-01)*X+2.78360021977405D+01)*X -            &
           1.04891729356965D+03)*X+2.36985942687423D+04)*X -            &
           3.19504627257548D+05)*X+2.34879693563358D+06)*X -            &
           7.16341568174085D+06)*E + R35/(X-R35)                         
      RT4 = ((((((((-6.38526371092582D-05*X-2.29263585792626D-03)*X -   &
           7.65735935499627D-02)*X+9.12692349152792D+00)*X -            &
           2.32077034386717D+02)*X+2.81839578728845D+02)*X +            &
           9.59529683876419D+04)*X-1.77638956809518D+06)*X +            &
           1.02489759645410D+07)*E + R45/(X-R45)                         
      RT5 = ((((((((-3.59049364231569D-05*X-2.25963977930044D-02)*X +   &
           1.12594870794668D+00)*X-4.56752462103909D+01)*X +            &
           1.05804526830637D+03)*X-1.16003199605875D+04)*X -            &
           4.07297627297272D+04)*X+2.22215528319857D+06)*X -            &
           1.61196455032613D+07)*E + R55/(X-R55)                         
      WW5 = (((((((((-4.61100906133970D-10*X+1.43069932644286D-07)*X -  &
           1.63960915431080D-05)*X+1.15791154612838D-03)*X -            &
           5.30573476742071D-02)*X+1.61156533367153D+00)*X -            &
           3.23248143316007D+01)*X+4.12007318109157D+02)*X -            &
           3.02260070158372D+03)*X+9.71575094154768D+03)*E + W55*WW1     
      WW4 = (((((((((-2.40799435809950D-08*X+8.12621667601546D-06)*X -  &
           9.04491430884113D-04)*X+6.37686375770059D-02)*X -            &
           2.96135703135647D+00)*X+9.15142356996330D+01)*X -            &
           1.86971865249111D+03)*X+2.42945528916947D+04)*X -            &
           1.81852473229081D+05)*X+5.96854758661427D+05)*E + W45*WW1     
      WW3 = (((((((( 1.83574464457207D-05*X-1.54837969489927D-03)*X +   &
           1.18520453711586D-01)*X-6.69649981309161D+00)*X +            &
           2.44789386487321D+02)*X-5.68832664556359D+03)*X +            &
           8.14507604229357D+04)*X-6.55181056671474D+05)*X +            &
           2.26410896607237D+06)*E + W35*WW1                             
      WW2 = (((((((( 2.77778345870650D-05*X-2.22835017655890D-03)*X +   &
           1.61077633475573D-01)*X-8.96743743396132D+00)*X +            &
           3.28062687293374D+02)*X-7.65722701219557D+03)*X +            &
           1.10255055017664D+05)*X-8.92528122219324D+05)*X +            &
           3.10638627744347D+06)*E + W25*WW1                             
      WW1 = WW1-0.01962D+00*E-WW2-WW3-WW4-WW5                            
      RETURN
                                                                        
  240 IF (X > 59.0D+00) GO TO 260                                        
!     X=40.0 TO 59.0                             NROOTS = 5              
      XXX = X**3                                                         
      E = XXX*EXP(-X)                                                    
      RT1 = (((-2.43758528330205D-02*X+2.07301567989771D+00)*X -        &
           6.45964225381113D+01)*X+7.14160088655470D+02)*E + R15/(X-R15) 
      RT2 = (((-2.28861955413636D-01*X+1.93190784733691D+01)*X -        &
           5.99774730340912D+02)*X+6.61844165304871D+03)*E + R25/(X-R25) 
      RT3 = (((-6.95053039285586D-01*X+5.76874090316016D+01)*X -        &
           1.77704143225520D+03)*X+1.95366082947811D+04)*E + R35/(X-R35) 
      RT4 = (((-1.58072809087018D+00*X+1.27050801091948D+02)*X -        &
           3.86687350914280D+03)*X+4.23024828121420D+04)*E + R45/(X-R45) 
      RT5 = (((-3.33963830405396D+00*X+2.51830424600204D+02)*X -        &
           7.57728527654961D+03)*X+8.21966816595690D+04)*E + R55/(X-R55) 
      E = XXX*E                                                          
      WW5 = (( 1.35482430510942D-08*X-3.27722199212781D-07)*X +         &
           2.41522703684296D-06)*E + W55*WW1                             
      WW4 = (( 1.23464092261605D-06*X-3.55224564275590D-05)*X +         &
           3.03274662192286D-04)*E + W45*WW1                             
      WW3 = (( 1.34547929260279D-05*X-4.19389884772726D-04)*X +         &
           3.87706687610809D-03)*E + W35*WW1                             
      WW2 = (( 2.09539509123135D-05*X-6.87646614786982D-04)*X +         &
           6.68743788585688D-03)*E + W25*WW1
      WW1 = WW1-WW2-WW3-WW4-WW5
      RETURN

!     X=59.0 TO INFINITY                         NROOTS = 5
  260 RT1 = R15/(X-R15)
      RT2 = R25/(X-R25)
      RT3 = R35/(X-R35)
      RT4 = R45/(X-R45)
      RT5 = R55/(X-R55)
      WW2 = W25*WW1
      WW3 = W35*WW1
      WW4 = W45*WW1
      WW5 = W55*WW1
      WW1 = WW1-WW2-WW3-WW4-WW5
      RETURN
      END

! Calling by HSandT & ERISPDFGHIL: ROOT6 
      SUBROUTINE ROOT6
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXAUX=55)
      DIMENSION RGRID(MXAUX),WGRID(MXAUX),P0(MXAUX),P1(MXAUX),          &
                P2(MXAUX),RTS(13),WTS(13),WRK(13),ALPHA(0:12),BETA(0:12)
      COMMON /ROOT  / XX,UF(13),WF(13),NROOTS
      COMMON /RYSPAR/ XASYMP(13),RTSASY(13,13),WTSASY(13,13),           &
                      NAUXS(13),MAPRYS(13),RTSAUX(55,8),WTSAUX(55,8)
!-----------------------------------------------------------------------
      IF(XX>=XASYMP(NROOTS)) THEN
       FACTR = 1.0d0/XX
       FACTW = SQRT(FACTR)
       DO I=1,NROOTS
         RTS(I)= FACTR * RTSASY(I,NROOTS)
         WTS(I)= FACTW * WTSASY(I,NROOTS)
       ENDDO
      ELSE
       NAUX=NAUXS(NROOTS)
       MAP=MAPRYS(NROOTS)
       DO I=1,NAUX
          T2 = RTSAUX(I,MAP)*RTSAUX(I,MAP)
          RGRID(I) = T2
          WGRID(I) = WTSAUX(I,MAP)*EXP(-XX*T2)
       ENDDO
       EPS = 1.0D-14
       CALL RYSDS(NROOTS,NAUX,RGRID,WGRID,ALPHA,BETA,IERR,P0,P1,P2)
       CALL RYSGW(NROOTS,ALPHA,BETA,EPS,RTS,WTS,IERR,WRK)
      END IF
      DO K=1,NROOTS
       DUM  = RTS(K)
       UF(K)= DUM/(1.0d0-DUM)
       WF(K)= WTS(K)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! RYSDS
      SUBROUTINE RYSDS(N,NCAP,X,W,ALPHA,BETA,IERR,P0,P1,P2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(NCAP),W(NCAP),ALPHA(N),BETA(N),P0(NCAP),              &
                P1(NCAP),P2(NCAP)

      TINY = 1.0D-40
      HUGE = 1.0D+40

      IERR=0
      IF(N<=0 .or. N>NCAP) THEN
        IERR=1
        RETURN
      END IF
      NM1=N-1

      SUM0=0.0D+00
      SUM1=0.0D+00
      DO 10 M=1,NCAP
        SUM0=SUM0+W(M)
        SUM1=SUM1+W(M)*X(M)
   10 CONTINUE
      ALPHA(1)=SUM1/SUM0
      BETA(1)=SUM0
      IF(N==1) RETURN

      DO 20 M=1,NCAP
        P1(M)=0.0D+00
        P2(M)=1.0D+00
   20 CONTINUE
      DO 40 K=1,NM1
        SUM1=0.0D+00
        SUM2=0.0D+00
        DO 30 M=1,NCAP

          IF(W(M)==0.0D+00) GOTO 30
          P0(M)=P1(M)
          P1(M)=P2(M)
          P2(M)=(X(M)-ALPHA(K))*P1(M)-BETA(K)*P0(M)

          IF(ABS(P2(M))>HUGE .or. ABS(SUM2)>HUGE) THEN
            IERR=K
            RETURN
          END IF
          T=W(M)*P2(M)*P2(M)
          SUM1=SUM1+T
          SUM2=SUM2+T*X(M)
   30   CONTINUE

        IF(ABS(SUM1)<TINY) THEN
          IERR=-K
          RETURN
        END IF
        ALPHA(K+1)=SUM2/SUM1
        BETA(K+1)=SUM1/SUM0
        SUM0=SUM1
   40 CONTINUE

      RETURN
      END
      
! RYSGW
      SUBROUTINE RYSGW(N,ALPHA,BETA,EPS,ROOTS,WEIGHT,IERR,WRK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ALPHA(N),BETA(N),ROOTS(N),WEIGHT(N),WRK(N)

      IF(N<1) THEN
        IERR=-1
        RETURN
      END IF

      IERR=0
      ROOTS(1)=ALPHA(1)
      IF(BETA(1)<0.0D+00) THEN
        IERR=-2
        RETURN
      END IF
      WEIGHT(1)=BETA(1)
      IF (N==1) RETURN

      WEIGHT(1)=1.0D+00
      WRK(N)=0.0D+00
      DO 100 K=2,N
        ROOTS(K)=ALPHA(K)
        IF(BETA(K)<0.0D+00) THEN
          IERR=-2
          RETURN
        END IF
        WRK(K-1)=SQRT(BETA(K))
        WEIGHT(K)=0.0D+00
  100 CONTINUE

      DO 240 L=1,N
        J=0

  105   DO 110 M=L,N
         IF(M==N) GO TO 120
         IF(ABS(WRK(M))<=EPS*(ABS(ROOTS(M))+ABS(ROOTS(M+1))))GOTO 120
  110   CONTINUE
  120   DP=ROOTS(L)
        IF(M==L) GO TO 240
        IF(J==30) GO TO 400
        J=J+1

        DG=(ROOTS(L+1)-DP)/(2.0D+00*WRK(L))
        DR=SQRT(DG*DG+1.0D+00)
        DG=ROOTS(M)-DP+WRK(L)/(DG+SIGN(DR,DG))
        DS=1.0D+00
        DC=1.0D+00
        DP=0.0D+00
        MML=M-L

        DO 200 II=1,MML
          I=M-II
          DF=DS*WRK(I)
          DB=DC*WRK(I)
          IF(ABS(DF)<ABS(DG)) GO TO 150
          DC=DG/DF
          DR=SQRT(DC*DC+1.0D+00)
          WRK(I+1)=DF*DR
          DS=1.0D+00/DR
          DC=DC*DS
          GO TO 160
  150     DS=DF/DG
          DR=SQRT(DS*DS+1.0D+00)
          WRK(I+1)=DG*DR
          DC=1.0D+00/DR
          DS=DS*DC
  160     DG=ROOTS(I+1)-DP
          DR=(ROOTS(I)-DG)*DS+2.0D+00*DC*DB
          DP=DS*DR
          ROOTS(I+1)=DG+DP
          DG=DC*DR-DB

          DF=WEIGHT(I+1)
          WEIGHT(I+1)=DS*WEIGHT(I)+DC*DF
          WEIGHT(I)=DC*WEIGHT(I)-DS*DF
  200   CONTINUE
        ROOTS(L)=ROOTS(L)-DP
        WRK(L)=DG
        WRK(M)=0.0D+00
        GO TO 105
  240 CONTINUE

      DO 300 II=2,N
        I=II-1
        K=I
        DP=ROOTS(I)
        DO 260 J=II,N
          IF(ROOTS(J)>=DP) GO TO 260
          K=J
          DP=ROOTS(J)
  260   CONTINUE
        IF(K==I) GO TO 300
        ROOTS(K)=ROOTS(I)
        ROOTS(I)=DP
        DP=WEIGHT(I)
        WEIGHT(I)=WEIGHT(K)
        WEIGHT(K)=DP
  300 CONTINUE
      DO 310 K=1,N
        WEIGHT(K)=BETA(1)*WEIGHT(K)*WEIGHT(K)
  310 CONTINUE
      RETURN

  400 IERR=L

      RETURN
      END

! ExchangeInt                                           
      SUBROUTINE ExchangeInt(XINTS,GHONDO,NSH2,MAXG,EX,CS,CP,CD,CF,CG,  &
                             CH,CI,NPRIMI,KSTART,KATOM,KTYPE,KNG,KLOC,  &
                             KMIN,KMAX,NSHELL,Cxyz,NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER, DIMENSION(4,3) :: IB                                                       
      COMMON/INTDEX1/IJGT(784),KLGT(784)   
      COMMON/SHLEXC/NORGSH(3),NORGSP(3),IEXCH,NGTH(4)  
      COMMON/SHLNOS1/QQ4,IJKL
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DIMENSION XINTS(NSH2),GHONDO(MAXG)
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DDIJ
      LOGICAL ISHandJSH
!-----------------------------------------------------------------------
      CALL BASCHK(LMAXIMA,KTYPE,NSHELL)
      NANGM =  4                                          
      IF(LMAXIMA==2) NANGM =  6                                          
      IF(LMAXIMA==3) NANGM = 10                                          
      IF(LMAXIMA==4) NANGM = 15                                          
      IF(LMAXIMA==5) NANGM = 21                                          
      IF(LMAXIMA==6) NANGM = 28                                          
      NGTH(4) = 1                                                       
      NGTH(3) = NGTH(4) * NANGM                                         
      NGTH(2) = NGTH(3) * NANGM                                         
      NGTH(1) = NGTH(2) * NANGM                                         
      DO I=1,3                                                       
       NORGSH(I) = 0                                               
       NORGSP(I) = 0                                               
      ENDDO                                                          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IEXCH = 1                                                         
      QQ4   = 1.0d0                                                       
      NINTEG  = 0  
!- - - - - - - - - - - - - - - - - - - - - - - -
      IJIJ = 0                                                          
      DO ISH = 1,NSHELL                                             
       DO JSH = 1,ISH                                             
        IJIJ = IJIJ+1                                               
        ALLOCATE(DDIJ(49*900))                                                      
        CALL SHELLS(1,ISH,JSH,ISH,JSH,.TRUE.,EX,CS,CP,CD,CF,CG,CH,CI,   &
                    NPRIMI,KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,&
                    Cxyz)                    
        CALL IJPRIM(DDIJ)                                           
        CALL SHELLS(2,ISH,JSH,ISH,JSH,.TRUE.,EX,CS,CP,CD,CF,CG,CH,CI,   &
                    NPRIMI,KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,&
                    Cxyz)                       
        CALL ZQOUT(GHONDO,MAXG)                                          
        IF(IJKL==1) CALL S0000(GHONDO,DDIJ)                       
        IF(IJKL>1) CALL ERISPDFGHIL(GHONDO,DDIJ)
        DEALLOCATE(DDIJ)
        VMAX = 0.0D+00                                                    
        MINI = KMIN(ISH)                                               
        MINJ = KMIN(JSH)                                               
        MAXI = KMAX(ISH)                                               
        JMAX = KMAX(JSH)                                               
        ISHandJSH=ISH==JSH                                               
        IBB = IB(1,IEXCH)                                              
        JBB = IB(2,IEXCH)                                              
        KBB = IB(3,IEXCH)                                              
        LBB = IB(4,IEXCH)                                              
        IJN = 0                                                        
        DO I=MINI,MAXI                                             
         IF(ISHandJSH) JMAX = I                                          
         DO J=MINJ,JMAX                                          
          IJN = IJN+1                                           
          NN = IJGT(IJN) + KLGT(IJN)                            
          VAL = GHONDO(NN)                                      
          IF(VAL>0.0D+00)NINTEG = NINTEG + 1                              
          IF(VAL>VMAX)VMAX=VAL                                 
         END DO
        END DO
        XINTS(IJIJ) = SQRT(VMAX)                                         
       END DO
      END DO
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               

! QOUT                                             
      SUBROUTINE QOUT(BUFP,IX,BUFP2,IX2,NINTEGtm,NINTMX,GHONDO,IDONTW,  &
                      KLOC,KMIN,KMAX,NSHELL)                            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL IANDJ,KANDL,SAME
      INTEGER,DIMENSION(NSHELL) :: KLOC,KMIN,KMAX
      DIMENSION BUFP(NINTMX),IX(NINTMX),GHONDO(*)
      COMMON /ERIOUT/ ISH,JSH,KSH,LSH,LSTRI,LSTRJ,LSTRK,LSTRL           
      COMMON /MISC  / IANDJ,KANDL,SAME                                  
      COMMON /RESTAR/ NREC,IST,JST,KST,LST                               
      INTEGER,DIMENSION(NINTEGtm) :: IX2                                 
      DOUBLE PRECISION,DIMENSION(NINTEGtm) :: BUFP2                      
      COMMON /SHLNOS/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,              &
                      MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          &
                      NIJ,IJ,KL                                    
      COMMON /SHLT  / SHLTOL,CUTOFF,ICOUNT
      SAVE IJN,KLN                                                      
!-----------------------------------------------------------------------                                                                       
!     Pack 4-indices into 1 word. 
!     Write Label & integral on Unit 1 if DONTW = .False.
!-----------------------------------------------------------------------
      SAME  = ISH == KSH .and. JSH == LSH                           
      IANDJ = ISH == JSH                                              
      KANDL = KSH == LSH                                              
      MINI = KMIN(ISH)                                                  
      MINJ = KMIN(JSH)                                                  
      MINK = KMIN(KSH)                                                  
      MINL = KMIN(LSH)                                                  
      MAXI = KMAX(ISH)                                                  
      MAXJ = KMAX(JSH)                                                  
      MAXK = KMAX(KSH)                                                  
      MAXL = KMAX(LSH)                                                  
      LOCI = KLOC(ISH)-MINI                                             
      LOCJ = KLOC(JSH)-MINJ                                             
      LOCK = KLOC(KSH)-MINK                                             
      LOCL = KLOC(LSH)-MINL          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                               
      IJN = 0                                                           
      JMAX = MAXJ                                                       
      DO I = MINI,MAXI                                              
       I_INDEX = (I-MINI)*LSTRI + 1                                   
       IF (IANDJ) JMAX = I                                            
       DO 1 J = MINJ,JMAX                                           
        IJ_INDEX = (J-MINJ)*LSTRJ + I_INDEX                         
        IJN = IJN+1                                                 
        LMAX = MAXL                                                 
        KLN = 0                                                     
        DO K =  MINK,MAXK                                       
         IJK_INDEX = (K-MINK)*LSTRK + IJ_INDEX                    
         IF (KANDL) LMAX = K                                      
         DO L = MINL,LMAX                                     
          KLN = KLN+1                                           
          IF(SAME.and.KLN>IJN)GO TO 1                 
          IJKL_INDEX = (L-MINL)*LSTRL + IJK_INDEX               
          VAL = GHONDO( IJKL_INDEX ) 
          IF(ABS(VAL)>=CUTOFF)THEN
           I1 = LOCI+I                                           
           I2 = LOCJ+J                                           
           I3 = LOCK+K                                           
           I4 = LOCL+L                                           
           IF (I1 >= I2) GO TO 100                             
           N = I1                                                
           I1 = I2                                               
           I2 = N                                                
  100      IF (I3 >= I4) GO TO 120                             
           N = I3                                                
           I3 = I4                                               
           I4 = N                                                
  120      IF (I1-I3) 140,160,180                                
  140      N = I1                                                
           I1 = I3                                               
           I3 = N                                                
           N = I2                                                
           I2 = I4                                               
           I4 = N                                                
           GO TO 180                                             
  160      IF (I2 < I4) GO TO 140                             
  180      CONTINUE                                              
!                                                                
           IF(I1 == I2) VAL = VAL*0.5D0                        
           IF(I3 == I4) VAL = VAL*0.5D0                        
           IF(I1 == I3 .and. I2 == I4) VAL = VAL*0.5D0       
!                                                                
           NPACK = ICOUNT                                        
           IPACK = I1                                            
           JPACK = I2                                            
           KPACK = I3                                            
           LPACK = I4                                            
           LABEL = ISHFT( IPACK, 48 ) + ISHFT( JPACK, 32 ) +            &
                   ISHFT( KPACK, 16 ) + LPACK                 
           IX(NPACK) = LABEL                                  
           BUFP(ICOUNT) = VAL 
           ICOUNT = ICOUNT+1                                     
           IF(ICOUNT > 0) THEN                           
            JCOUNT = ICOUNT                               
            IF(JCOUNT > NINTMX) THEN                        
             NXInteg = NINTMX
             IF(IDONTW==1)THEN
              IJBUFi = (NREC-1)*NINTMX
              do ibuf=1,NINTMX
               IX2  (IJBUFi+ibuf) = IX(ibuf)
               BUFP2(IJBUFi+ibuf) = BUFP(ibuf)
              end do
             ELSE
              WRITE(1)NXInteg,IX,BUFP
             END IF
             ICOUNT = 1                               
             NREC = NREC+1                                   
            END IF
           END IF
          END IF
         END DO
        END DO
    1  CONTINUE                                                       
      END DO
!-----------------------------------------------------------------------                                                                       
      RETURN                                                            
      END                                                               

! SHELLS                                           
      SUBROUTINE SHELLS(NELEC,ISH,JSH,KSH,LSH,FLIP,EX,CS,CP,CD,CF,CG,CH,&
                 CI,NPRIMI,KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,&
                 Cxyz)                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL       NORM
      COMMON/NORMAL/NORM
      LOGICAL     IANDJ,KANDL,SAME
      COMMON/MISC/IANDJ,KANDL,SAME
      COMMON /ERIOUT/ INU,JNU,KNU,LNU,NGTI,NGTJ,NGTK,NGTL                
      COMMON /INTDEX/ IJX(784),IJY(784),IJZ(784),IK(784),               &
                      KLX(784),KLY(784),KLZ(784)                         
      COMMON/INTDEX1/IJGT(784),KLGT(784)                                 
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      COMMON /ROOT  / XX,U(13),W(13),NROOTS                              
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NGTH(4)                  
      COMMON /SHLINF/  GA(30),CSA(30),CPA(30),CDA(30),                  &
                      CFA(30),CGA(30),CHA(30),CIA(30),                  &
                       GB(30),CSB(30),CPB(30),CDB(30),                  &
                      CFB(30),CGB(30),CHB(30),CIB(30),                  &
                       GC(30),CSC(30),CPC(30),CDC(30),                  &
                      CFC(30),CGC(30),CHC(30),CIC(30),                  &
                       GD(30),CSD(30),CPD(30),CDD(30),                  &
                      CFD(30),CGD(30),CHD(30),CID(30),                  &
                      AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,                   &
                      DX,DY,DZ,RCD,NGA,NGB,NGC,NGD                                    
      COMMON /SHLNOS/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,              &
                      MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          &
                      NIJ,IJ,KL                                          
      COMMON /SHLNOS1/QQ4,IJKL
!
      LOGICAL FLIP                                                      
      INTEGER,DIMENSION(NSHELL)::KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DIMENSION IX(84),IY(84),IZ(84),JX(84),JY(84),JZ(84),              &
                KX(84),KY(84),KZ(84),LX(84),LY(84),LZ(84)                                     
!
      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,       &
                  3,   0,   0,   2,   2,   1,   0,   1,   0,   1,       &
                  4,   0,   0,   3,   3,   1,   0,   1,   0,   2,       &
                  2,   0,   2,   1,   1,                                &
                  5,   0,   0,   4,   4,   1,   0,   1,   0,   3,       &
                  3,   2,   0,   2,   0,   3,   1,   1,   2,   2,       &
                  1,                                                    &
                  6,   0,   0,   5,   5,   1,   0,   1,   0,   4,       &
                  4,   2,   0,   2,   0,   4,   1,   1,   3,   3,       &
                  0,   3,   3,   2,   1,   2,   1,   2/                  
      DATA KX /   0,   7,   0,   0,  14,   0,   0,   7,   7,   0,       &
                 21,   0,   0,  14,  14,   7,   0,   7,   0,   7,       &
                 28,   0,   0,  21,  21,   7,   0,   7,   0,  14,       &
                 14,   0,  14,   7,   7,                                &
                 35,   0,   0,  28,  28,   7,   0,   7,   0,  21,       &
                 21,  14,   0,  14,   0,  21,   7,   7,  14,  14,       &
                  7,                                                    &
                 42,   0,   0,  35,  35,   7,   0,   7,   0,  28,       &
                 28,  14,   0,  14,   0,  28,   7,   7,  21,  21,       &
                  0,  21,  21,  14,   7,  14,   7,  14/                  
      DATA JX /   0,  49,   0,   0,  98,   0,   0,  49,  49,   0,       &
                147,   0,   0,  98,  98,  49,   0,  49,   0,  49,       &
                196,   0,   0, 147, 147,  49,   0,  49,   0,  98,       &
                 98,   0,  98,  49,  49,                                &
                245,   0,   0, 196, 196,  49,   0,  49,   0, 147,       &
                147,  98,   0,  98,   0, 147,  49,  49,  98,  98,       &
                 49,                                                    &
                294,   0,   0, 245, 245,  49,   0,  49,   0, 196,       &
                196,  98,   0,  98,   0, 196,  49,  49, 147, 147,       &
                  0, 147, 147,  98,  49,  98,  49,  98/                  
      DATA IX /   1, 344,   1,   1, 687,   1,   1, 344, 344,   1,       &
               1030,   1,   1, 687, 687, 344,   1, 344,   1, 344,       &
               1373,   1,   1,1030,1030, 344,   1, 344,   1, 687,       &
                687,   1, 687, 344, 344,                                &
               1716,   1,   1,1373,1373, 344,   1, 344,   1,1030,       &
               1030, 687,   1, 687,   1,1030, 344, 344, 687, 687,       &
                344,                                                    &
               2059,   1,   1,1716,1716, 344,   1, 344,   1,1373,       &
               1373, 687,   1, 687,   1,1373, 344, 344,1030,1030,       &
                  1,1030,1030, 687, 344, 687, 344, 687/                  
      DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,       &
                  0,   3,   0,   1,   0,   2,   2,   0,   1,   1,       &
                  0,   4,   0,   1,   0,   3,   3,   0,   1,   2,       &
                  0,   2,   1,   2,   1,                                &
                  0,   5,   0,   1,   0,   4,   4,   0,   1,   2,       &
                  0,   3,   3,   0,   2,   1,   3,   1,   2,   1,       &
                  2,                                                    &
                  0,   6,   0,   1,   0,   5,   5,   0,   1,   2,       &
                  0,   4,   4,   0,   2,   1,   4,   1,   3,   0,       &
                  3,   2,   1,   3,   3,   1,   2,   2/                  
      DATA KY /   0,   0,   7,   0,   0,  14,   0,   7,   0,   7,       &
                  0,  21,   0,   7,   0,  14,  14,   0,   7,   7,       &
                  0,  28,   0,   7,   0,  21,  21,   0,   7,  14,       &
                  0,  14,   7,  14,   7,                                &
                  0,  35,   0,   7,   0,  28,  28,   0,   7,  14,       &
                  0,  21,  21,   0,  14,   7,  21,   7,  14,   7,       &
                 14,                                                    &
                  0,  42,   0,   7,   0,  35,  35,   0,   7,  14,       &
                  0,  28,  28,   0,  14,   7,  28,   7,  21,   0,       &
                 21,  14,   7,  21,  21,   7,  14,  14/                  
      DATA JY /   0,   0,  49,   0,   0,  98,   0,  49,   0,  49,       &
                  0, 147,   0,  49,   0,  98,  98,   0,  49,  49,       &
                  0, 196,   0,  49,   0, 147, 147,   0,  49,  98,       &
                  0,  98,  49,  98,  49,                                &
                  0, 245,   0,  49,   0, 196, 196,   0,  49,  98,       &
                  0, 147, 147,   0,  98,  49, 147,  49,  98,  49,       &
                 98,                                                    &
                  0, 294,   0,  49,   0, 245, 245,   0,  49,  98,       &
                  0, 196, 196,   0,  98,  49, 196,  49, 147,   0,       &
                147,  98,  49, 147, 147,  49,  98,  98/                  
      DATA IY /   1,   1, 344,   1,   1, 687,   1, 344,   1, 344,       &
                  1,1030,   1, 344,   1, 687, 687,   1, 344, 344,       &
                  1,1373,   1, 344,   1,1030,1030,   1, 344, 687,       &
                  1, 687, 344, 687, 344,                                &
                  1,1716,   1, 344,   1,1373,1373,   1, 344, 687,       &
                  1,1030,1030,   1, 687, 344,1030, 344, 687, 344,       &
                687,                                                    &
                  1,2059,   1, 344,   1,1716,1716,   1, 344, 687,       &
                  1,1373,1373,   1, 687, 344,1373, 344,1030,   1,       &
               1030, 687, 344,1030,1030, 344, 687, 687/                  
      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,       &
                  0,   0,   3,   0,   1,   0,   1,   2,   2,   1,       &
                  0,   0,   4,   0,   1,   0,   1,   3,   3,   0,       &
                  2,   2,   1,   1,   2,                                &
                  0,   0,   5,   0,   1,   0,   1,   4,   4,   0,       &
                  2,   0,   2,   3,   3,   1,   1,   3,   1,   2,       &
                  2,                                                    &
                  0,   0,   6,   0,   1,   0,   1,   5,   5,   0,       &
                  2,   0,   2,   4,   4,   1,   1,   4,   0,   3,       &
                  3,   1,   2,   1,   2,   3,   3,   2/                  
      DATA KZ /   0,   0,   0,   7,   0,   0,  14,   0,   7,   7,       &
                  0,   0,  21,   0,   7,   0,   7,  14,  14,   7,       &
                  0,   0,  28,   0,   7,   0,   7,  21,  21,   0,       &
                 14,  14,   7,   7,  14,                                &
                  0,   0,  35,   0,   7,   0,   7,  28,  28,   0,       &
                 14,   0,  14,  21,  21,   7,   7,  21,   7,  14,       &
                 14,                                                    &
                  0,   0,  42,   0,   7,   0,   7,  35,  35,   0,       &
                 14,   0,  14,  28,  28,   7,   7,  28,   0,  21,       &
                 21,   7,  14,   7,  14,  21,  21,  14/                  
      DATA JZ /   0,   0,   0,  49,   0,   0,  98,   0,  49,  49,       &
                  0,   0, 147,   0,  49,   0,  49,  98,  98,  49,       &
                  0,   0, 196,   0,  49,   0,  49, 147, 147,   0,       &
                 98,  98,  49,  49,  98,                                &
                  0,   0, 245,   0,  49,   0,  49, 196, 196,   0,       &
                 98,   0,  98, 147, 147,  49,  49, 147,  49,  98,       &
                 98,                                                    &
                  0,   0, 294,   0,  49,   0,  49, 245, 245,   0,       &
                 98,   0,  98, 196, 196,  49,  49, 196,   0, 147,       &
                147,  49,  98,  49,  98, 147, 147,  98/                  
      DATA IZ /   1,   1,   1, 344,   1,   1, 687,   1, 344, 344,       &
                  1,   1,1030,   1, 344,   1, 344, 687, 687, 344,       &
                  1,   1,1373,   1, 344,   1, 344,1030,1030,   1,       &
                687, 687, 344, 344, 687,                                &
                  1,   1,1716,   1, 344,   1, 344,1373,1373,   1,       &
                687,   1, 687,1030,1030, 344, 344,1030, 344, 687,       &
                687,                                                    &
                  1,   1,2059,   1, 344,   1, 344,1716,1716,   1,       &
                687,   1, 687,1373,1373, 344, 344,1373,   1,1030,       &
               1030, 344, 687, 344, 687,1030,1030, 687/                 
!-----------------------------------------------------------------------                                                                       
!     PREPARE SHELL INFORMATION/FOR HONDO INTEGRATION 
      NORM = .TRUE.  
      IF(NELEC==2) GO TO 200                                          
!                                                                       
!     ----- PERMUTE ISH AND JSH SHELLS, FOR THEIR TYPE                  
!     THIS IS DONE FOR SPEED REASONS.  THE CODE GETS THE RIGHT ANSWER   
!     WITHOUT THE ANGULAR MOMENTUM FLIPPING, AND THEREFORE A CALLING    
!     ARGUMENT ALLOWS ONE DO EXACTLY THE integral BLOCK AS SPECIFIED,   
!     SHOULD THAT BE DESIRED.                                           
!                                                                       
      IANDJ = ISH == JSH                                              
      IF (KTYPE(ISH) < KTYPE(JSH)  .and.  FLIP) THEN                 
       INU = JSH                                                      
       JNU = ISH                                                      
       NGTI = NGTH(2)                                                 
       NGTJ = NGTH(1)                                                 
      ELSE                                                              
       INU = ISH                                                      
       JNU = JSH                                                      
       NGTI = NGTH(1)                                                 
       NGTJ = NGTH(2)                                                 
      END IF                                                            
!                                                                       
!     ----- ISHELL                                                      
!                                                                       
      I = KATOM(INU)                                                    
      AX = Cxyz(1,I)                                                       
      AY = Cxyz(2,I)                                                       
      AZ = Cxyz(3,I)                                                       
      I1 = KSTART(INU)                                                  
      I2 = I1+KNG(INU)-1                                                
      LIT = KTYPE(INU)                                                  
      MINI = KMIN(INU)                                                  
      MAXI = KMAX(INU)                                                  
      LOCI = KLOC(INU)-MINI                                             
      NGA = 0                                                           
      DO I = I1,I2                                                  
       NGA = NGA+1                                                    
       GA(NGA)  = EX(I)                                                
       CSA(NGA) = CS(I)                                               
       CPA(NGA) = CP(I)                                               
       CDA(NGA) = CD(I)                                               
       CFA(NGA) = CF(I)                                               
       CGA(NGA) = CG(I)                                               
       CHA(NGA) = CH(I)                                               
       CIA(NGA) = CI(I)                                               
      END DO                                                          
!                                                                       
!     ----- JSHELL                                                      
!                                                                       
      J = KATOM(JNU)                                                    
      BX = Cxyz(1,J)                                                       
      BY = Cxyz(2,J)                                                       
      BZ = Cxyz(3,J)                                                       
      J1 = KSTART(JNU)                                                  
      J2 = J1+KNG(JNU)-1                                                
      LJT = KTYPE(JNU)                                                  
      MINJ = KMIN(JNU)                                                  
      MAXJ = KMAX(JNU)                                                  
      LOCJ = KLOC(JNU)-MINJ                                             
      NGB = 0                                                           
      DO J = J1,J2                                                  
       NGB = NGB+1                                                    
       GB(NGB) = EX(J)                                                
       CSB(NGB) = CS(J)                                               
       CPB(NGB) = CP(J)                                               
       CDB(NGB) = CD(J)                                               
       CFB(NGB) = CF(J)                                               
       CGB(NGB) = CG(J)                                               
       CHB(NGB) = CH(J)                                               
       CIB(NGB) = CI(J)                                               
      END DO                                                          
      RAB = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))       
!                                                                       
!     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS                
!                                                                       
      IJ = 0                                                            
      JMAX = MAXJ                                                       
      DO I = MINI,MAXI                                              
       NX = IX(I)                                                     
       NY = IY(I)                                                     
       NZ = IZ(I)                                                     
       IF (IANDJ) JMAX = I                                            
       DO J = MINJ,JMAX                                           
        IJ = IJ+1                                                   
        IJX(IJ) = NX+JX(J)                                          
        IJY(IJ) = NY+JY(J)                                          
        IJZ(IJ) = NZ+JZ(J)                                          
        IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1                    
       END DO                                                       
      END DO                                                          
      RETURN                                                            
!     ******                                                            
!                                                                       
!        K AND L SHELL                                                  
!                                                                       
  200 CONTINUE                                                          
      KANDL = KSH == LSH                                              
      SAME = ISH == KSH .and. JSH == LSH                            
!                                                                       
!     ----- PERMUTE KSH AND LSH SHELLS, FOR THEIR TYPE                  
!                                                                       
      IF (KTYPE(KSH) < KTYPE(LSH)  .and.  FLIP) THEN                 
       KNU = LSH                                                      
       LNU = KSH                                                      
       NGTK = NGTH(4)                                                 
       NGTL = NGTH(3)                                                 
      ELSE                                                              
       KNU = KSH                                                      
       LNU = LSH                                                      
       NGTK = NGTH(3)                                                 
       NGTL = NGTH(4)                                                 
      END IF                                                            
!                                                                       
!     ----- K SHELL                                                     
!                                                                       
      K = KATOM(KNU)                                                    
      CX = Cxyz(1,K)                                                       
      CY = Cxyz(2,K)                                                       
      CZ = Cxyz(3,K)                                                       
      K1 = KSTART(KNU)                                                  
      K2 = K1+KNG(KNU)-1                                                
      LKT = KTYPE(KNU)                                                  
      MINK = KMIN(KNU)                                                  
      MAXK = KMAX(KNU)                                                  
      LOCK = KLOC(KNU)-MINK                                             
      NGC = 0                                                           
      DO K = K1,K2                                                  
       NGC = NGC+1                                                    
       GC(NGC)  = EX(K)                                                
       CSC(NGC) = CS(K)                                               
       CPC(NGC) = CP(K)                                               
       CDC(NGC) = CD(K)                                               
       CFC(NGC) = CF(K)                                               
       CGC(NGC) = CG(K)                                               
       CHC(NGC) = CH(K)                                               
       CIC(NGC) = CI(K)                                               
      END DO                                                          
!                                                                       
!     ----- LSHELL                                                      
!                                                                       
      L = KATOM(LNU)                                                    
      DX = Cxyz(1,L)                                                       
      DY = Cxyz(2,L)                                                       
      DZ = Cxyz(3,L)                                                       
      L1 = KSTART(LNU)                                                  
      L2 = L1+KNG(LNU)-1                                                
      LLT = KTYPE(LNU)                                                  
      MINL = KMIN(LNU)                                                  
      MAXL = KMAX(LNU)                                                  
      LOCL = KLOC(LNU)-MINL                                             
      NGD = 0                                                           
      DO L = L1,L2                                                  
       NGD = NGD+1                                                    
       GD(NGD) = EX(L)                                                
       CSD(NGD) = CS(L)                                               
       CPD(NGD) = CP(L)                                               
       CDD(NGD) = CD(L)                                               
       CFD(NGD) = CF(L)                                               
       CGD(NGD) = CG(L)                                               
       CHD(NGD) = CH(L)                                               
       CID(NGD) = CI(L)                                               
      END DO                                                          
      NROOTS = (LIT+LJT+LKT+LLT-2)/2                                    
      RCD = ((CX-DX)*(CX-DX) + (CY-DY)*(CY-DY) + (CZ-DZ)*(CZ-DZ))       
!                                                                       
!     ----- PREPARE INDICES FOR PAIRS OF (K,L) FUNCTIONS                
!                                                                       
      KL = 0                                                            
      LMAX = MAXL                                                       
      DO K = MINK,MAXK                                              
       NX = KX(K)                                                     
       NY = KY(K)                                                     
       NZ = KZ(K)                                                     
       IF (KANDL) LMAX = K                                            
       DO L = MINL,LMAX                                           
        KL = KL+1                                                   
        KLX(KL) = NX+LX(L)                                          
        KLY(KL) = NY+LY(L)                                          
        KLZ(KL) = NZ+LZ(L)                                          
        KLGT(KL) = NGTK*(K-MINK)+NGTL*(L-MINL)                      
       END DO                                                       
      END DO                                                          
      MAX = KL                                                          
      DO 320 I = 1,IJ                                                   
      IF (SAME) MAX = I                                                 
  320 IK(I) = MAX                                                       
      IJKL = IJ*KL                                                      
      IF (SAME) IJKL = IJ*(IJ+1)/2                                      
      RETURN                                                            
      END                                                               

! IJPRIM                                           
      SUBROUTINE IJPRIM(DDIJ)                                           
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL       NORM
      COMMON/NORMAL/NORM
      LOGICAL     IANDJ,KANDL,SAME 
      COMMON/MISC/IANDJ,KANDL,SAME
      COMMON/IJGNRL/A(900),R(900),X1(900),Y1(900),Z1(900),IJD(784)                                             
      COMMON/SHLINF/ AG(30),CSA(30),CPA(30),CDA(30),                    &
                    CFA(30),CGA(30),CHA(30),CIA(30),                    &
                     BG(30),CSB(30),CPB(30),CDB(30),                    &
                    CFB(30),CGB(30),CHB(30),CIB(30),                    &
                     CG(30),CSC(30),CPC(30),CDC(30),                    &
                    CFC(30),CGC(30),CHC(30),CIC(30),                    &
                     DG(30),CSD(30),CPD(30),CDD(30),                    &
                    CFD(30),CGD(30),CHD(30),CID(30),                    &
                    XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,                     &
                    XL,YL,ZL,RRK,NGA,NGB,NGC,NGD                                      
      COMMON/SHLNOS/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,                &
                    MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,            &
                    NIJ,IJ,KL                                          
      COMMON/SHLT/SHLTOL,CUTOFF,ICOUNT
!
      DIMENSION DDIJ(49*900)                                           
      PARAMETER (SQRT3=1.73205080756888D0,SQRT5=2.23606797749979D0,     &                       
                 SQRT7=2.64575131106459D0,SQRT11=3.3166247903553998D0)                             
!-----------------------------------------------------------------------
      MAX = MAXJ                                                        
      N = 0                                                             
      NN = 0                                                            
      NM = -2**20                                                       
      DO 180 I = MINI,MAXI                                              
         GO TO (100,100,120,120,100,120,120,100,120,120,                &
                100,120,120,100,120,120,120,120,120,100,                &
                100,120,120,100,120,120,120,120,120,100,                &
                120,120,100,120,120,                                    &
                100,120,120,100,120,120,120,120,120,100,                &
                120,120,120,120,120,100,120,120,100,120,                &
                120,                                                    &
                100,120,120,100,120,120,120,120,120,100,                &
                120,120,120,120,120,100,120,120,100,120,                &
                120,100,120,120,120,120,120,100),I                       
  100    NM = NN                                                         
  120    NN = NM                                                         
         IF (IANDJ) MAX = I                                              
         DO 170 J = MINJ,MAX                                             
            GO TO (140,140,160,160,140,160,160,140,160,160,             &
                   140,160,160,140,160,160,160,160,160,140,             &
                   140,160,160,140,160,160,160,160,160,140,             &
                   160,160,140,160,160,                                 &
                   140,160,160,140,160,160,160,160,160,140,             &
                   160,160,160,160,160,140,160,160,140,160,             &
                   160,                                                 &
                   140,160,160,140,160,160,160,160,160,140,             &
                   160,160,160,160,160,140,160,160,140,160,             &
                   160,140,160,160,160,160,160,140),J                   
  140       NN = NN+1                                                   
  160       N = N+1                                                     
            IJD(N) = NN                                                 
  170    CONTINUE                                                       
  180 CONTINUE                                                          
!                                                                       
!     ----- I PRIMITIVE                                                 
!                                                                       
      NIJ = 0                                                           
      JBMAX = NGB                                                       
      DO 540 IA = 1,NGA                                                 
         AI = AG(IA)                                                    
         ARRI = AI*RRI                                                  
         AXI = AI*XI                                                    
         AYI = AI*YI                                                    
         AZI = AI*ZI                                                    
         CSI = CSA(IA)                                                  
         CPI = CPA(IA)                                                  
         CDI = CDA(IA)                                                  
         CFI = CFA(IA)                                                  
         CGI = CGA(IA)                                                  
         CHI = CHA(IA)                                                  
         CII = CIA(IA)                                                  
!                                                                       
!        ----- J PRIMITIVE                                              
!                                                                       
         IF (IANDJ) JBMAX = IA                                          
         DO 520 JB = 1,JBMAX                                            
            AJ = BG(JB)                                                 
            AA = AI+AJ                                                  
            AAINV = 1.0d0/AA                                              
            DUM = AJ*ARRI*AAINV                                         
            IF (DUM > SHLTOL) GO TO 520                                 
            CSJ = CSB(JB)                                               
            CPJ = CPB(JB)                                               
            CDJ = CDB(JB)                                               
            CFJ = CFB(JB)                                               
            CGJ = CGB(JB)                                               
            CHJ = CHB(JB)                                               
            CIJ = CIB(JB)                                               
            NM = 49*NIJ                                                 
            NN = NM                                                     
            NIJ = NIJ+1                                                 
            R(NIJ) = DUM                                                
            A(NIJ) = AA                                                 
            X1(NIJ) = (AXI+AJ*XJ)*AAINV                                 
            Y1(NIJ) = (AYI+AJ*YJ)*AAINV                                 
            Z1(NIJ) = (AZI+AJ*ZJ)*AAINV                                 
!                                                                       
!           ----- DENSITY FACTOR                                        
!                                                                       
            DUM1 = 0.0d0                                                 
            DUM2 = 0.0d0                                                 
            DO 420 I = MINI,MAXI                                        
               GO TO (200,220,420,420,240,420,420,260,420,420,          &
                      261,420,420,262,420,420,420,420,420,263,          &
                      264,420,420,265,420,420,420,420,420,266,          &
                      420,420,267,420,420,                              &
                      268,420,420,269,420,420,420,420,420,270,          &
                      420,420,420,420,420,271,420,420,272,420,          &
                      420,                                              &
                      273,420,420,274,420,420,420,420,420,275,          &
                      420,420,420,420,420,276,420,420,277,420,          &
                      420,278,420,420,420,420,420,279),I                
  200          DUM1 = CSI*AAINV                                         
               GO TO 280                                                
  220          DUM1 = CPI*AAINV                                         
               GO TO 280                                                
  240          DUM1 = CDI*AAINV                                         
               GO TO 280                                                
  260          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  261          DUM1 = CFI*AAINV                                         
               GO TO 280                                                
  262          IF (NORM) DUM1 = DUM1*SQRT5                              
               GO TO 280                                                
  263          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  264          DUM1 = CGI*AAINV                                         
               GO TO 280                                                
  265          IF (NORM) DUM1 = DUM1*SQRT7                              
               GO TO 280                                                
  266          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
               GO TO 280                                                
  267          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  268          DUM1 = CHI*AAINV                                         
               GO TO 280                                                
  269          IF (NORM) DUM1 = DUM1*3.0d0                              
               GO TO 280                                                
  270          IF (NORM) DUM1 = DUM1*SQRT7/SQRT3                        
               GO TO 280                                                
  271          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  272          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
               GO TO 280                                                
  273          DUM1 = CII*AAINV                                         
               GO TO 280                                                
  274          IF (NORM) DUM1 = DUM1*SQRT11                             
               GO TO 280                                                
  275          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  276          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  277          IF (NORM) DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)                
               GO TO 280                                                
  278          IF (NORM) DUM1 = DUM1*SQRT5                              
               GO TO 280                                                
  279          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
!                                                                       
  280          IF (IANDJ) MAX = I                                       
               DO 400 J = MINJ,MAX                                      
                  GO TO (300,320,400,400,340,400,400,360,400,400,       &
                         361,400,400,362,400,400,400,400,400,363,       &
                         364,400,400,365,400,400,400,400,400,366,       &
                         400,400,367,400,400,                           &
                         368,400,400,369,400,400,400,400,400,370,       &
                         400,400,400,400,400,371,400,400,372,400,       &
                         400,                                           &
                         373,400,400,374,400,400,400,400,400,375,       &
                         400,400,400,400,400,376,400,400,377,400,       &
                         400,378,400,400,400,400,400,379),J             
  300             DUM2 = DUM1*CSJ                                       
                  GO TO 380                                             
  320             DUM2 = DUM1*CPJ                                       
                  GO TO 380                                             
  340             DUM2 = DUM1*CDJ                                       
                  GO TO 380                                             
  360             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  361             DUM2 = DUM1*CFJ                                       
                  GO TO 380                                             
  362             IF (NORM) DUM2 = DUM2*SQRT5                           
                  GO TO 380                                             
  363             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  364             DUM2 = DUM1*CGJ                                       
                  GO TO 380                                             
  365             IF (NORM) DUM2 = DUM2*SQRT7                           
                  GO TO 380                                             
  366             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
                  GO TO 380                                             
  367             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  368             DUM2 = DUM1*CHJ                                       
                  GO TO 380                                             
  369             IF (NORM) DUM2 = DUM2*3.0d0                           
                  GO TO 380                                             
  370             IF (NORM) DUM2 = DUM2*SQRT7/SQRT3                     
                  GO TO 380                                             
  371             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  372             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
                  GO TO 380                                             
  373             DUM2 = DUM1*CIJ                                       
                  GO TO 380                                             
  374             IF (NORM) DUM2 = DUM2*SQRT11                          
                  GO TO 380                                             
  375             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  376             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  377             IF (NORM) DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)             
                  GO TO 380                                             
  378             IF (NORM) DUM2 = DUM2*SQRT5                           
                  GO TO 380                                             
  379             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
!                                                                       
  380             NN = NN+1                                             
                  DDIJ(NN) = DUM2                                       
  400          CONTINUE                                                 
  420       CONTINUE                                                    
            IF ( .NOT. IANDJ) GO TO 520                                 
            IF (IA == JB) GO TO 520                                   
            GO TO (500,440,460,455,450,445,444),LIT                     
  440       IF (MINI == 2) GO TO 500                                  
            DDIJ(NM+2) = DDIJ(NM+2)+CSI*CPJ*AAINV                       
            GO TO 480                                                   
  444       DDIJ(NM+28) = DDIJ(NM+28)+DDIJ(NM+28)                       
            DDIJ(NM+27) = DDIJ(NM+27)+DDIJ(NM+27)                       
            DDIJ(NM+26) = DDIJ(NM+26)+DDIJ(NM+26)                       
            DDIJ(NM+25) = DDIJ(NM+25)+DDIJ(NM+25)                       
            DDIJ(NM+24) = DDIJ(NM+24)+DDIJ(NM+24)                       
            DDIJ(NM+23) = DDIJ(NM+23)+DDIJ(NM+23)                       
            DDIJ(NM+22) = DDIJ(NM+22)+DDIJ(NM+22)                       
            DDIJ(NM+21) = DDIJ(NM+21)+DDIJ(NM+21)                       
            DDIJ(NM+20) = DDIJ(NM+20)+DDIJ(NM+20)                       
            DDIJ(NM+19) = DDIJ(NM+19)+DDIJ(NM+19)                       
            DDIJ(NM+18) = DDIJ(NM+18)+DDIJ(NM+18)                       
            DDIJ(NM+17) = DDIJ(NM+17)+DDIJ(NM+17)                       
            DDIJ(NM+16) = DDIJ(NM+16)+DDIJ(NM+16)                       
  445       DDIJ(NM+15) = DDIJ(NM+15)+DDIJ(NM+15)                       
            DDIJ(NM+14) = DDIJ(NM+14)+DDIJ(NM+14)                       
            DDIJ(NM+13) = DDIJ(NM+13)+DDIJ(NM+13)                       
            DDIJ(NM+12) = DDIJ(NM+12)+DDIJ(NM+12)                       
            DDIJ(NM+11) = DDIJ(NM+11)+DDIJ(NM+11)                       
  450       DDIJ(NM+10) = DDIJ(NM+10)+DDIJ(NM+10)                       
            DDIJ(NM+9) = DDIJ(NM+9)+DDIJ(NM+9)                          
            DDIJ(NM+8) = DDIJ(NM+8)+DDIJ(NM+8)                          
            DDIJ(NM+7) = DDIJ(NM+7)+DDIJ(NM+7)                          
  455       DDIJ(NM+6) = DDIJ(NM+6)+DDIJ(NM+6)                          
            DDIJ(NM+5) = DDIJ(NM+5)+DDIJ(NM+5)                          
            DDIJ(NM+4) = DDIJ(NM+4)+DDIJ(NM+4)                          
  460       DDIJ(NM+2) = DDIJ(NM+2)+DDIJ(NM+2)                          
  480       DDIJ(NM+3) = DDIJ(NM+3)+DDIJ(NM+3)                          
  500       DDIJ(NM+1) = DDIJ(NM+1)+DDIJ(NM+1)                          
  520    CONTINUE                                                       
  540 CONTINUE
!-----------------------------------------------------------------------                                                          
      RETURN                                                            
      END                                                               

! S0000                                            
      SUBROUTINE S0000(GHONDO,DDIJ)                                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL IANDJ,KANDL,SAME
      DIMENSION GHONDO(*),DDIJ(49*900)                                 
      COMMON /IJGNRL/ A(900),R(900),X1(900),Y1(900),Z1(900),IJD(784)                                           
      COMMON /MISC  / IANDJ,KANDL,SAME                                   
      COMMON /SHLINF/  AG(30),CSA(30),CPA(30),CDA(30),                  &
                      CFA(30),CGA(30),CHA(30),CIA(30),                  &
                       BG(30),CSB(30),CPB(30),CDB(30),                  &
                      CFB(30),CGB(30),CHB(30),CIB(30),                  &
                       CG(30),CSC(30),CPC(30),CDC(30),                  &
                      CFC(30),CGC(30),CHC(30),CIC(30),                  &
                       DG(30),CSD(30),CPD(30),CDD(30),                  &
                      CFD(30),CGD(30),CHD(30),CID(30),                  &
                      XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,                   &
                      XL,YL,ZL,RRK,NGA,NGB,NGC,NGD                                    
      COMMON /SHLNOS/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,              &
                      MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          &
                      NIJ,IJ,KL                                    
      COMMON /SHLNOS1/QQ4,IJKL                         
      COMMON /SHLT  / SHLTOL,CUTOFF,ICOUNT
      PARAMETER (PI252=34.986836655250D+00, PIE4=7.85398163397448D-01)                             
!-----------------------------------------------------------------------                                                                       
!     SSSS integral for HONDO integrals          
!-----------------------------------------------------------------------                                                                       
      GGOUT = 0.0d0                                                      
      LGMAX = NGD                                                       
      DO 300 KG = 1,NGC                                                 
      BK = CG(KG)                                                       
      BRRK = BK*RRK                                                     
      BXK = BK*XK                                                       
      BYK = BK*YK                                                       
      BZK = BK*ZK                                                       
      CSK = CSC(KG)                                                     
      IF (KANDL) LGMAX = KG                                             
      DO 280 LG = 1,LGMAX                                               
      BL = DG(LG)                                                       
      BB = BK+BL                                                        
      BBINV = 1.0d0/BB                                                    
      DUM = BL*BRRK*BBINV                                               
      IF (DUM > SHLTOL) GO TO 280                                       
      BBRRK = DUM                                                       
      D2 = CSD(LG)*CSK*BBINV                                            
      IF (KANDL .and. LG /= KG) D2 = D2+D2                            
      BBX = (BXK+BL*XL)*BBINV                                           
      BBY = (BYK+BL*YL)*BBINV                                           
      BBZ = (BZK+BL*ZL)*BBINV                                           
      SUM = 0.0d0                                                        
      NN = 1                                                            
      DO 260 N = 1,NIJ                                                  
      DUM = BBRRK+R(N)                                                  
      IF (DUM > SHLTOL) GO TO 260                                       
      EXPE = EXP(-DUM)                                                  
      AA = A(N)                                                         
      AB = AA+BB                                                        
      DUM = X1(N)-BBX                                                   
      XX = DUM*DUM                                                      
      DUM = Y1(N)-BBY                                                   
      XX = DUM*DUM+XX                                                   
      DUM = Z1(N)-BBZ                                                   
      XX = DUM*DUM+XX                                                   
      X = XX*AA*BB/AB                                                   
!                                                                       
      IF (X > 5.0D+00) GO TO 160                                     
      IF (X > 1.0D+00) GO TO 120                                     
      IF (X > 3.0D-07) GO TO 100                                     
      WW1 = 1.0D+00-X/3.0D+00                                           
      GO TO 240                                                         
!                                                                       
  100 CONTINUE                                                          
      F1 = ((((((((-8.36313918003957D-08*X+1.21222603512827D-06 )*X-    &
           1.15662609053481D-05 )*X+9.25197374512647D-05 )*X-           &
           6.40994113129432D-04 )*X+3.78787044215009D-03 )*X-           &
           1.85185172458485D-02 )*X+7.14285713298222D-02 )*X-           &
           1.99999999997023D-01 )*X+3.33333333333318D-01                 
      WW1 = (X+X)*F1+EXP(-X)                                             
      GO TO 240                                                          
!                                                                        
  120 CONTINUE                                                           
      IF (X > 3.0D+00) GO TO 140                                         
      Y = X-2.0D+00                                                      
      F1 = ((((((((((-1.61702782425558D-10*Y+1.96215250865776D-09 )*Y-  &
           2.14234468198419D-08 )*Y+2.17216556336318D-07 )*Y-           &
           1.98850171329371D-06 )*Y+1.62429321438911D-05 )*Y-           &
           1.16740298039895D-04 )*Y+7.24888732052332D-04 )*Y-           &
           3.79490003707156D-03 )*Y+1.61723488664661D-02 )*Y-           &
           5.29428148329736D-02 )*Y+1.15702180856167D-01                 
      WW1 = (X+X)*F1+EXP(-X)                                             
      GO TO 240                                                          
!                                                                        
  140 CONTINUE                                                           
      Y = X-4.0D+00                                                      
      F1 = ((((((((((-2.62453564772299D-11*Y+3.24031041623823D-10 )*Y-  &
           3.614965656163D-09)*Y+3.760256799971D-08)*Y-                 &
           3.553558319675D-07)*Y+3.022556449731D-06)*Y-                 &
           2.290098979647D-05)*Y+1.526537461148D-04)*Y-                 &
           8.81947375894379D-04 )*Y+4.33207949514611D-03 )*Y-           &
           1.75257821619926D-02 )*Y+5.28406320615584D-02                 
      WW1 = (X+X)*F1+EXP(-X)                                             
      GO TO 240                                                          
!                                                                        
  160 CONTINUE                                                           
      IF (X > 15.0D+00) GO TO 200                                        
      E = EXP(-X)                                                        
      IF (X > 10.0D+00) GO TO 180                                        
      XINV = 1.0d0/X                                                      
      WW1 = (((((( 4.6897511375022D-01*XINV-6.9955602298985D-01)*XINV + &
           5.3689283271887D-01)*XINV-3.2883030418398D-01)*XINV +        &
           2.4645596956002D-01)*XINV-4.9984072848436D-01)*XINV -        &
           3.1501078774085D-06)*E + SQRT(PIE4*XINV)                      
      GO TO 240                                                          
!                                                                        
  180 CONTINUE                                                           
      XINV = 1.0d0/X                                                      
      WW1 = (((-1.8784686463512D-01*XINV+2.2991849164985D-01)*XINV      &
               -4.9893752514047D-01)*XINV-2.1916512131607D-05)*E        &
               + SQRT(PIE4*XINV)                                        
      GO TO 240                                                         
!                                                                       
  200 CONTINUE                                                          
      IF (X > 33.0D+00) GO TO 220                                    
      XINV = 1.0d0/X                                                      
      E = EXP(-X)                                                       
      WW1 = (( 1.9623264149430D-01*XINV-4.9695241464490D-01)*XINV -     &
                 6.0156581186481D-05)*E + SQRT(PIE4*XINV)                     
      GO TO 240                                                         
!                                                                       
  220 WW1 = SQRT(PIE4/X)                                                
!                                                                       
  240 CONTINUE                                                          
      SUM = SUM+DDIJ(NN)*WW1*EXPE/SQRT(AB)                
  260 NN = NN+49                                                        
      GGOUT = GGOUT+D2*SUM                                              
  280 CONTINUE                                                          
  300 CONTINUE                                                          
      GHONDO(1) = GGOUT*PI252*QQ4                                       
      RETURN                                                            
      END                                                               

! ERISPDFGHIL                                           
      SUBROUTINE ERISPDFGHIL(GHONDO,DDIJ)                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      LOGICAL       NORM
      COMMON/NORMAL/NORM
      LOGICAL     IANDJ,KANDL,SAME 
      COMMON/MISC/IANDJ,KANDL,SAME
      COMMON/DENS/DKL(784),DIJ(784)                                 
      COMMON/IJGNRL/AA(900),R(900),X1(900),Y1(900),Z1(900),IJD(784)                                             
      COMMON/ROOT/XX,U(13),W(13),NROOTS                                
      COMMON/SETINT/IN(13),KN(13),NI,NJ,NK,NL,NMAX,MMAX,                &
                    BP01,B00,B10,XCP00,XC00,YCP00,YC00,                 &
                    ZCP00,ZC00,F00,                                     &
                    DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL                        
      COMMON/SHLINF/ AG(30),CSA(30),CPA(30),CDA(30),                    &
                    CFA(30),CGA(30),CHA(30),CIA(30),                    &
                     BG(30),CSB(30),CPB(30),CDB(30),                    &
                    CFB(30),CGB(30),CHB(30),CIB(30),                    &
                     CG(30),CSC(30),CPC(30),CDC(30),                    &
                    CFC(30),CGC(30),CHC(30),CIC(30),                    &
                     DG(30),CSD(30),CPD(30),CDD(30),                    &
                    CFD(30),CGD(30),CHD(30),CID(30),                    &
                    XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,                     &
                    XL,YL,ZL,RRK,NGA,NGB,NGC,NGD                                      
      COMMON/SHLNOS/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,                &
                    MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,            &
                    NIJ,IJ,KL                                          
      COMMON/SHLNOS1/QQ4,IJKL
      COMMON/SHLT/ SHLTOL,CUTOFF,ICOUNT
!
      LOGICAL DOBLE                          
      DIMENSION IN1(13),GHONDO(*),DDIJ(*)                                                  
      PARAMETER (SQRT3=1.73205080756888D0,SQRT5=2.23606797749979D0,     &
                 SQRT7=2.64575131106459D0,PI252=34.986836655250D0,      &
                 SQRT11=3.3166247903553998D0)               
!-----------------------------------------------------------------------
      FACTOR = PI252*QQ4                                                
      NI = LIT-1                                                        
      NJ = LJT-1                                                        
      NK = LKT-1                                                        
      NL = LLT-1                                                        
      DXIJ = XI-XJ                                                      
      DYIJ = YI-YJ                                                      
      DZIJ = ZI-ZJ                                                      
      DXKL = XK-XL                                                      
      DYKL = YK-YL                                                      
      DZKL = ZK-ZL                                                      
      NMAX = NI+NJ                                                      
      MMAX = NK+NL                                                      
      MAX = NMAX+1                                                      
      DO I = 1,MAX                                                  
       N = I-1                                                        
       IF (N <= NI) IN1(I) = 343*N+1                                
       IF (N > NI) IN1(I) = 343*NI+49*(N-NI)+1                     
      END DO
      MAX = MMAX+1                                                      
      DO K = 1,MAX                                                  
       N = K-1                                                        
       IF (N <= NK) KN(K) = 7*N                                     
       IF (N > NK) KN(K) = 7*NK+N-NK                               
      END DO
      LGMAX = NGD                                                       
      DO KG = 1,NGC                        !      K Primitive                                                 
       AK = CG(KG)                                                    
       BRRK = AK*RRK                                                  
       AKXK = AK*XK                                                   
       AKYK = AK*YK                                                   
       AKZK = AK*ZK                                                   
       CSK = CSC(KG)*FACTOR                                           
       CPK = CPC(KG)*FACTOR                                           
       CDK = CDC(KG)*FACTOR                                           
       CFK = CFC(KG)*FACTOR                                           
       CGK = CGC(KG)*FACTOR                                           
       CHK = CHC(KG)*FACTOR                                           
       CIK = CIC(KG)*FACTOR                                           
       IF (KANDL) LGMAX = KG                                          
       DO LG = 1,LGMAX                     !      L Primitive
        AL = DG(LG)                                                 
        B = AK+AL                                                   
        BINV = 1.0d0/B                                                
        BBRRK = AL*BRRK*BINV                                        
        IF(BBRRK<=SHLTOL)THEN
         CSL = CSD(LG)                                               
         CPL = CPD(LG)                                               
         CDL = CDD(LG)                                               
         CFL = CFD(LG)                                               
         CGL = CGD(LG)                                               
         CHL = CHD(LG)                                               
         CIL = CID(LG)                                               
         XB = (AKXK+AL*XL)*BINV                                      
         YB = (AKYK+AL*YL)*BINV                                      
         ZB = (AKZK+AL*ZL)*BINV                                      
         BXBK = B*(XB-XK)                                            
         BYBK = B*(YB-YK)                                            
         BZBK = B*(ZB-ZK)                                            
         BXBI = B*(XB-XI)                                            
         BYBI = B*(YB-YI)                                            
         BZBI = B*(ZB-ZI)                                            
!        DENSITY FACTOR                                        
         DOBLE = KANDL.and.KG/=LG                                   
         N = 0                                                       
         MAX = MAXL                                                  
         DUM1 = 0.0d0                                                 
         DUM2 = 0.0d0                                                 
         DO K = MINK,MAXK                                        
          GO TO (140,160,220,220,180,220,220,200,220,220,               &
                 201,220,220,202,220,220,220,220,220,203,               &
                 204,220,220,205,220,220,220,220,220,206,               &
                 220,220,207,220,220,                                   &
                 208,220,220,209,220,220,220,220,220,210,               &
                 220,220,220,220,220,211,220,220,212,220,               &
                 220,                                                   &
                 213,220,220,214,220,220,220,220,220,215,               &
                 220,220,220,220,220,216,220,220,217,220,               &
                 220,218,220,220,220,220,220,219),K                      
  140     DUM1 = CSK*BINV                                          
          GO TO 220                                                
  160     DUM1 = CPK*BINV                                          
          GO TO 220                                                
  180     DUM1 = CDK*BINV                                          
          GO TO 220                                                
  200     IF (NORM) DUM1 = DUM1*SQRT3                              
          GO TO 220                                                
  201     DUM1 = CFK*BINV                                          
          GO TO 220                                                
  202     IF (NORM) DUM1 = DUM1*SQRT5                              
          GO TO 220                                                
  203     IF (NORM) DUM1 = DUM1*SQRT3                              
          GO TO 220                                                
  204     DUM1 = CGK*BINV                                          
          GO TO 220                                                
  205     IF (NORM) DUM1 = DUM1*SQRT7                              
          GO TO 220                                                
  206     IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
          GO TO 220                                                
  207     IF (NORM) DUM1 = DUM1*SQRT3                              
          GO TO 220                                                
  208     DUM1 = CHK*BINV                                          
          GO TO 220                                                
  209     IF (NORM) DUM1 = DUM1*3.0d0                              
          GO TO 220                                                
  210     IF (NORM) DUM1 = DUM1*SQRT7/SQRT3                        
          GO TO 220                                                
  211     IF (NORM) DUM1 = DUM1*SQRT3                              
          GO TO 220                                                
  212     IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
          GO TO 220                                                
  213     DUM1 = CIK*BINV                                          
          GO TO 220                                                
  214     IF (NORM) DUM1 = DUM1*SQRT11                             
          GO TO 220                                                
  215     IF (NORM) DUM1 = DUM1*SQRT3                              
          GO TO 220                                                
  216     IF (NORM) DUM1 = DUM1*SQRT3                              
          GO TO 220                                                
  217     IF (NORM) DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)                
          GO TO 220                                                
  218     IF (NORM) DUM1 = DUM1*SQRT5                              
          GO TO 220                                                
  219     IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
!                                                                  
  220     IF (KANDL) MAX = K                                       
          DO L = MINL,MAX                                      
           GO TO (240,280,340,340,300,340,340,320,340,340,              &
                  321,340,340,322,340,340,340,340,340,323,              &
                  324,340,340,325,340,340,340,340,340,326,              &
                  340,340,327,340,340,                                  &
                  328,340,340,329,340,340,340,340,340,330,              &
                  340,340,340,340,340,331,340,340,332,340,              &
                  340,                                                  &
                  333,340,340,334,340,340,340,340,340,335,              &
                  340,340,340,340,340,336,340,340,337,340,              &
                  340,338,340,340,340,340,340,339),L             
  240      DUM2 = DUM1*CSL                                       
           IF ( .NOT. DOBLE) GO TO 340                          
           IF (K > 1) GO TO 260                               
           DUM2 = DUM2+DUM2                                      
           GO TO 340                                             
  260      DUM2 = DUM2+CSK*CPL*BINV                              
           GO TO 340                                             
  280      DUM2 = DUM1*CPL                                       
           IF (DOBLE) DUM2 = DUM2+DUM2                          
           GO TO 340                                             
  300      DUM2 = DUM1*CDL                                       
           IF (DOBLE) DUM2 = DUM2+DUM2                          
           GO TO 340                                             
  320      IF (NORM) DUM2 = DUM2*SQRT3                           
           GO TO 340                                             
  321      DUM2 = DUM1*CFL                                       
           IF (DOBLE) DUM2 = DUM2+DUM2                          
           GO TO 340                                             
  322      IF (NORM) DUM2 = DUM2*SQRT5                           
           GO TO 340                                             
  323      IF (NORM) DUM2 = DUM2*SQRT3                           
           GO TO 340                                             
  324      DUM2 = DUM1*CGL                                       
           IF (DOBLE) DUM2 = DUM2+DUM2                          
           GO TO 340                                             
  325      IF (NORM) DUM2 = DUM2*SQRT7                           
           GO TO 340                                             
  326      IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
           GO TO 340                                             
  327      IF (NORM) DUM2 = DUM2*SQRT3                           
           GO TO 340                                             
  328      DUM2 = DUM1*CHL                                       
           IF (DOBLE) DUM2 = DUM2+DUM2                          
           GO TO 340                                             
  329      IF (NORM) DUM2 = DUM2*3.0d0                           
           GO TO 340                                             
  330      IF (NORM) DUM2 = DUM2*SQRT7/SQRT3                     
           GO TO 340                                             
  331      IF (NORM) DUM2 = DUM2*SQRT3                           
           GO TO 340                                             
  332      IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
           GO TO 340                                             
  333      DUM2 = DUM1*CIL                                       
           IF (DOBLE) DUM2 = DUM2+DUM2                          
           GO TO 340                                             
  334      IF (NORM) DUM2 = DUM2*SQRT11                          
           GO TO 340                                             
  335      IF (NORM) DUM2 = DUM2*SQRT3                           
           GO TO 340                                             
  336      IF (NORM) DUM2 = DUM2*SQRT3                           
           GO TO 340                                             
  337      IF (NORM) DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)             
           GO TO 340                                             
  338      IF (NORM) DUM2 = DUM2*SQRT5                           
           GO TO 340                                             
  339      IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
!                                                                
  340      N = N+1                                               
           DKL(N) = DUM2                                         
          END DO
         END DO
!
         NN = 0                                                      
         DO N = 1,NIJ                      !         I,J Primitives                                
          DUM = BBRRK+R(N)                                         
          IF(DUM<=SHLTOL)THEN
           DO I = 1,IJ                                          
            DIJ(I) = DDIJ(IJD(I)+NN)                              
           END DO
           A = AA(N)                                                
           AB = A*B                                                 
           AANDB = A+B                                              
           EXPE = EXP(-DUM)/SQRT(AANDB)                             
           RHO = AB/AANDB                                           
           XA = X1(N)                                               
           YA = Y1(N)                                               
           ZA = Z1(N)                                               
           XX = RHO*((XA-XB)*(XA-XB)+(YA-YB)*(YA-YB)+(ZA-ZB)*(ZA-ZB))             
           AXAK = A*(XA-XK)                                         
           AYAK = A*(YA-YK)                                         
           AZAK = A*(ZA-ZK)                                         
           AXAI = A*(XA-XI)                                         
           AYAI = A*(YA-YI)                                         
           AZAI = A*(ZA-ZI)                                         
           C1X = BXBK+AXAK                                          
           C2X = A*BXBK                                             
           C3X = BXBI+AXAI                                          
           C4X = B*AXAI                                             
           C1Y = BYBK+AYAK                                          
           C2Y = A*BYBK                                             
           C3Y = BYBI+AYAI                                          
           C4Y = B*AYAI                                             
           C1Z = BZBK+AZAK                                          
           C2Z = A*BZBK                                             
           C3Z = BZBI+AZAI                                          
           C4Z = B*AZAI                                             
           IF (NROOTS <= 3) CALL RT123                            
           IF (NROOTS == 4) CALL ROOT4                            
           IF (NROOTS == 5) CALL ROOT5                            
           IF (NROOTS >= 6) CALL ROOT6                            
           MM = 0                                                   
           MAX = NMAX+1                                             
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!          ERI for each root
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           DO M = 1,NROOTS                                      
            U2 = U(M)*RHO                                         
            F00 = EXPE*W(M)                                       
            DO I = 1,MAX                                      
             IN(I) = IN1(I)+MM                                  
            END DO
            DUMINV = 1.0d0/(AB+U2*AANDB)                         
            DM2INV = 0.5D0*DUMINV                               
            BP01 = (A+U2)*DM2INV                               
            B00 = U2*DM2INV                                    
            B10 = (B+U2)*DM2INV                                
            XCP00 = (U2*C1X+C2X)*DUMINV                        
            XC00 = (U2*C3X+C4X)*DUMINV                         
            YCP00 = (U2*C1Y+C2Y)*DUMINV                        
            YC00 = (U2*C3Y+C4Y)*DUMINV                         
            ZCP00 = (U2*C1Z+C2Z)*DUMINV                        
            ZC00 = (U2*C3Z+C4Z)*DUMINV                         
            CALL XYZINT                                           
            MM = MM+2401                                          
           END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!          Form (I,J//K,L) integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           CALL FormIntegrals(GHONDO)                                       
          END IF
          NN = NN + 49
         END DO
        END IF                                                  
       END DO
      END DO
!                                                                       
      RETURN                                                            
      END                                                               

! Calling by ERISPDFGHIL XYZINT                                           
      SUBROUTINE XYZINT                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL N0,N1,M0,M1,FIRST1,FIRST2,FIRST3,FIRST4                    
      COMMON /SETINT/I(13),K(13),NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX,     &
                     BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00, &
                     DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL                     
      COMMON /XYZ   / XINT(31213),YINT(31213),ZINT(31213)               
!                                                                       
      N0 = NMAX == 0                                                  
      N1 = NMAX <= 1                                                  
      M0 = MMAX == 0                                                  
      M1 = MMAX <= 1                                                  
!                                                                       
!     ----- I(0,0) -----                                                
!                                                                       
      I1 = I(1)                                                         
      XINT(I1) = 1.0d0                                                    
      YINT(I1) = 1.0d0                                                    
      ZINT(I1) = F00                                                    
      IF (N0 .and. M0) RETURN                                           
      I2 = I(2)                                                         
      K2 = K(2)                                                         
      CP10 = B00                                                        
!                                                                       
!     ----- I(1,0) -----                                                
!                                                                       
      IF (.NOT. N0) THEN                                                
        XINT(I2) = XC00                                                 
        YINT(I2) = YC00                                                 
        ZINT(I2) = ZC00*F00                                             
        IF (M0) GO TO 120                                               
      END IF                                                            
!                                                                       
!     ----- I(0,1) -----                                                
!                                                                       
      I3 = I1+K2                                                        
      XINT(I3) = XCP00                                                  
      YINT(I3) = YCP00                                                  
      ZINT(I3) = ZCP00*F00                                              
!                                                                       
!     ----- I(1,1) -----                                                
!                                                                       
      IF (.NOT. N0) THEN                                                
        I3 = I2+K2                                                      
        XINT(I3) = XCP00*XINT(I2)+CP10                                  
        YINT(I3) = YCP00*YINT(I2)+CP10                                  
        ZINT(I3) = ZCP00*ZINT(I2)+CP10*F00                              
      END IF                                                            
!                                                                       
  120 CONTINUE                                                          
      IF (.NOT. N1) THEN                                                
        C10 = 0.0d0                                                      
        I3 = I1                                                         
        I4 = I2                                                         
        DO 160 N = 2,NMAX                                               
          C10 = C10+B10                                                 
!                                                                       
!     ----- I(N,0) -----                                                
!                                                                       
          I5 = I(N+1)                                                   
          XINT(I5) = C10*XINT(I3)+XC00*XINT(I4)                         
          YINT(I5) = C10*YINT(I3)+YC00*YINT(I4)                         
          ZINT(I5) = C10*ZINT(I3)+ZC00*ZINT(I4)                         
          IF ( .NOT. M0) THEN                                           
            CP10 = CP10+B00                                             
!                                                                       
!     ----- I(N,1) -----                                                
!                                                                       
            I3 = I5+K2                                                  
            XINT(I3) = XCP00*XINT(I5)+CP10*XINT(I4)                     
            YINT(I3) = YCP00*YINT(I5)+CP10*YINT(I4)                     
            ZINT(I3) = ZCP00*ZINT(I5)+CP10*ZINT(I4)                     
          END IF                                                        
          I3 = I4                                                       
          I4 = I5                                                       
  160     CONTINUE                                                      
      END IF                                                            
      IF ( .NOT. M1) THEN                                               
        CP01 = 0.0d0                                                     
        C01 = B00                                                       
        I3 = I1                                                         
        I4 = I1+K2                                                      
        DO 220 M = 2,MMAX                                               
          CP01 = CP01+BP01                                              
!                                                                       
!     ----- I(0,M) -----                                                
!                                                                       
          I5 = I1+K(M+1)                                                
          XINT(I5) = CP01*XINT(I3)+XCP00*XINT(I4)                       
          YINT(I5) = CP01*YINT(I3)+YCP00*YINT(I4)                       
          ZINT(I5) = CP01*ZINT(I3)+ZCP00*ZINT(I4)                       
!                                                                       
!     ----- I(1,M) -----                                                
!                                                                       
          IF (.NOT. N0) THEN                                            
            C01 = C01+B00                                               
            I3 = I2+K(M+1)                                              
            XINT(I3) = XC00*XINT(I5)+C01*XINT(I4)                       
            YINT(I3) = YC00*YINT(I5)+C01*YINT(I4)                       
            ZINT(I3) = ZC00*ZINT(I5)+C01*ZINT(I4)                       
          END IF                                                        
          I3 = I4                                                       
          I4 = I5                                                       
  220   CONTINUE                                                        
      END IF                                                            
!                                                                       
!     ----- I(N,M) -----                                                
!                                                                       
      IF (.NOT. N1 .and. .NOT. M1) THEN                                 
        C01 = B00                                                       
        K3 = K2                                                         
        DO 280 M = 2,MMAX                                               
          K4 = K(M+1)                                                   
          C01 = C01+B00                                                 
          I3 = I1                                                       
          I4 = I2                                                       
          C10 = B10                                                     
          DO 260 N = 2,NMAX                                             
            I5 = I(N+1)                                                 
            XINT(I5+K4) = C10*XINT(I3+K4)+XC00*XINT(I4+K4)              &
                          +C01*XINT(I4+K3)                              
            YINT(I5+K4) = C10*YINT(I3+K4)+YC00*YINT(I4+K4)              &
                          +C01*YINT(I4+K3)                              
            ZINT(I5+K4) = C10*ZINT(I3+K4)+ZC00*ZINT(I4+K4)              &
                          +C01*ZINT(I4+K3)                              
            C10 = C10+B10                                               
            I3 = I4                                                     
            I4 = I5                                                     
  260     CONTINUE                                                      
          K3 = K4                                                       
  280   CONTINUE                                                        
      END IF                                                            
!                                                                       
!     ----- I(NI,NJ,M) -----                                            
!                                                                       
      IF (NJMAX > 0) THEN                                            
        M = 0                                                           
        I5 = I(NMAX+1)                                                  
        FIRST1 = .TRUE.                                                 
        DO 430 WHILE (FIRST1 .or. M <= MMAX)                          
          MIN = NIMAX                                                   
          KM = K(M+1)                                                   
          FIRST2 = .TRUE.                                               
          DO 360 WHILE (FIRST2 .or. MIN < NMAX)                      
            N = NMAX                                                    
            I3 = I5+KM                                                  
            FIRST3 = .TRUE.                                             
            DO 340 WHILE (FIRST3 .or. N > MIN)                       
              I4 = I(N)+KM                                              
              XINT(I3) = XINT(I3)+DXIJ*XINT(I4)                         
              YINT(I3) = YINT(I3)+DYIJ*YINT(I4)                         
              ZINT(I3) = ZINT(I3)+DZIJ*ZINT(I4)                         
              I3 = I4                                                   
              N = N-1                                                   
              FIRST3 = .FALSE.                                          
  340       END DO                                                      
            MIN = MIN+1                                                 
            FIRST2 = .FALSE.                                            
  360     END DO                                                        
          IF (NIMAX > 0) THEN                                        
            I3 = 49+KM+I1                                               
            DO 400 NJ = 1,NJMAX                                         
              I4 = I3                                                   
              DO 380 NI = 1,NIMAX                                       
                XINT(I4) = XINT(I4+294)+DXIJ*XINT(I4-49)                
                YINT(I4) = YINT(I4+294)+DYIJ*YINT(I4-49)                
                ZINT(I4) = ZINT(I4+294)+DZIJ*ZINT(I4-49)                
                I4 = I4+343                                             
  380         CONTINUE                                                  
              I3 = I3+49                                                
  400       CONTINUE                                                    
          END IF                                                        
          M = M+1                                                       
          FIRST1 = .FALSE.                                              
  430   END DO                                                          
      END IF                                                            
!                                                                       
!     ----- I(NI,NJ,NK,NL) -----                                        
!                                                                       
      IF (NLMAX > 0) THEN                                            
        I5 = K(MMAX+1)                                                  
        IA = I1                                                         
        NI = 0                                                          
        FIRST4 = .TRUE.                                                 
        DO 580 WHILE (FIRST4 .or. NI <= NIMAX)                        
          NJ = 0                                                        
          IB = IA                                                       
          FIRST1 = .TRUE.                                               
          DO 570 WHILE (FIRST1 .or. NJ <= NJMAX)                      
            MIN = NKMAX                                                 
            FIRST2 = .TRUE.                                             
            DO 530 WHILE (FIRST2 .or. MIN < MMAX)                    
              M = MMAX                                                  
              I3 = IB+I5                                                
              FIRST3 = .TRUE.                                           
              DO 520 WHILE (FIRST3 .or. M > MIN)                     
                I4 = IB+K(M)                                            
                XINT(I3) = XINT(I3)+DXKL*XINT(I4)                       
                YINT(I3) = YINT(I3)+DYKL*YINT(I4)                       
                ZINT(I3) = ZINT(I3)+DZKL*ZINT(I4)                       
                I3 = I4                                                 
                M = M-1                                                 
                FIRST3 = .FALSE.                                        
  520         END DO                                                    
              MIN = MIN+1                                               
              FIRST2 = .FALSE.                                          
  530       END DO                                                      
            IF (NKMAX > 0) THEN                                      
              I3 = IB+1                                                 
              DO 560 NL = 1,NLMAX                                       
                I4 = I3                                                 
                DO 540 NK = 1,NKMAX                                     
                  XINT(I4) = XINT(I4+6)+DXKL*XINT(I4-1)                 
                  YINT(I4) = YINT(I4+6)+DYKL*YINT(I4-1)                 
                  ZINT(I4) = ZINT(I4+6)+DZKL*ZINT(I4-1)                 
                  I4 = I4+7                                             
  540           END DO                                                  
              I3 = I3+1                                                 
  560         END DO                                                    
            END IF                                                      
            NJ = NJ+1                                                   
            IB = IB+49                                                  
            FIRST1 = .FALSE.                                            
  570     END DO                                                        
          NI = NI+1                                                     
          IA = IA+343                                                   
          FIRST4 = .FALSE.                                              
  580   END DO                                                          
      END IF                                                            
!                                                                       
      RETURN                                                            
      END                                                               

! Calling by ERISPDFGHIL FormIntegrals                                            
      SUBROUTINE FormIntegrals(GHONDO)                                          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION GHONDO(*)                                               
      COMMON /DENS  / DKL(784),DIJ(784)                                 
      COMMON /INTDEX/ IJX(784),IJY(784),IJZ(784),IK(784),               &
                      KLX(784),KLY(784),KLZ(784)              
      COMMON/INTDEX1/IJGT(784),KLGT(784)           
      COMMON /SHLNOS/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,              &
                      MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          &
                      NIJ,IJ,KL                                    
      COMMON /ROOT  / XX,U(13),W(13),NROOTS                             
      COMMON /XYZ   / XIN(31213),YIN(31213),ZIN(31213)                  
!-----------------------------------------------------------------------
      GO TO (10,20,30,40,50,60,70,80,90,100,110,120,130),NROOTS         
!     NROOTS=1                                            
   10 CONTINUE                                                          
      DO I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*(XIN(MX)*YIN(MY)*ZIN(MZ))          
      END DO
      END DO
      RETURN                                                            
!     NROOTS=2                                            
   20 CONTINUE                                                          
      DO I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                &
                  ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           &
                +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401))          
      END DO
      END DO
      RETURN                                                            
!     NROOTS=3                                            
   30 CONTINUE                                                          
      DO I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                & 
                  ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           &
                +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           &
                +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802))          
      END DO
      END DO
      RETURN                                                            
!     NROOTS=4                                            
   40 CONTINUE                                                          
      DO I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                &
                  ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           &
                +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           &
                +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           &
                +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203))          
      END DO
      END DO
      RETURN                                                            
!     NROOTS=5                                            
   50 CONTINUE                                                          
      DO I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                &
                  ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           &
                +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           &
                +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           &
                +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           &
                +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604))          
      END DO
      END DO
      RETURN                                                            
!     NROOTS=6                                            
   60 CONTINUE                                                          
      DO I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                &
                  ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           &
                +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           &
                +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           &
                +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           &
                +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           &
                +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005))          
      END DO
      END DO
      RETURN                                                            
!     NROOTS=7                                            
   70 CONTINUE                                                          
      DO I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                &
                  ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           &
                +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           &
                +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           &
                +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           &
                +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           &
                +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           &
                +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406))          
      END DO
      END DO
      RETURN                                                            
!     NROOTS=8                                            
   80 CONTINUE                                                          
      DO I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                &
                  ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           &
                +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           &
                +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           &
                +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           &
                +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           &
                +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           &
                +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           &
                +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807))          
      END DO
      END DO
      RETURN                                                            
!     NROOTS=9                                            
   90 CONTINUE                                                          
      DO I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                &
                  ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           &
                +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           &
                +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           &
                +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           &
                +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           &
                +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           &
                +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           &
                +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           &
                +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208))          
      END DO
      END DO
      RETURN                                                            
!     NROOTS=10                                           
  100 CONTINUE                                                          
      DO I = 1,IJ                                                   
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO K = 1,MAX                                                  
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                &
                  ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           &
                +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           &
                +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           &
                +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           &
                +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           &
                +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           &
                +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           &
                +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           &
                +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)           &
                +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609))          
      END DO
      END DO
      RETURN                                                            
!     NROOTS=11                                           
  110 CONTINUE                                                          
      DO I = 1,IJ                                                   
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO K = 1,MAX                                                  
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                &
                  ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           &
                +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           &
                +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           &
                +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           &
                +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           &
                +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           &
                +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           &
                +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           &
                +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)           &
                +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609)           &
                +   XIN(MX+24010)*YIN(MY+24010)*ZIN(MZ+24010))          &
                *D1*DKL(K)+GHONDO(N)                                    
      END DO
      END DO
      RETURN                                                            
!     NROOTS=12                                           
  120 CONTINUE                                                          
      DO I = 1,IJ                                                   
       D1 = DIJ(I)                                                       
       NX = IJX(I)                                                       
       NY = IJY(I)                                                       
       NZ = IJZ(I)                                                       
       N1 = IJGT(I)                                                      
       MAX = IK(I)                                                       
       DO K = 1,MAX                                                  
        MX = NX+KLX(K)                                                    
        MY = NY+KLY(K)                                                    
        MZ = NZ+KLZ(K)                                                    
        N = N1+KLGT(K)                                                    
        GHONDO(N) = GHONDO(N) + D1*DKL(K)*                              &
                    ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )         &
                  +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)         &
                  +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)         &
                  +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)         &
                  +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)         &
                  +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)         &
                  +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)         &
                  +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)         &
                  +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)         &
                  +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609)         &
                  +   XIN(MX+24010)*YIN(MY+24010)*ZIN(MZ+24010)         &
                  +   XIN(MX+26411)*YIN(MY+26411)*ZIN(MZ+26411))        &
                  *D1*DKL(K)+GHONDO(N)                                    
       END DO
      END DO
      RETURN                                                            
!     NROOTS=13                                           
  130 CONTINUE
      DO I = 1,IJ                                                   
       D1 = DIJ(I)                                                       
       NX = IJX(I)                                                       
       NY = IJY(I)                                                       
       NZ = IJZ(I)                                                       
       N1 = IJGT(I)                                                      
       MAX = IK(I)                                                       
       DO K = 1,MAX                                                  
        MX = NX+KLX(K)                                                    
        MY = NY+KLY(K)                                                    
        MZ = NZ+KLZ(K)                                                    
        N = N1+KLGT(K)                                                    
        GHONDO(N) = GHONDO(N) + D1*DKL(K)*                              &
                    ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )         &
                  +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)         &
                  +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)         &
                  +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)         &
                  +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)         &
                  +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)         &
                  +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)         &
                  +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)         &
                  +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)         &
                  +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609)         &
                  +   XIN(MX+24010)*YIN(MY+24010)*ZIN(MZ+24010)         &
                  +   XIN(MX+26411)*YIN(MY+26411)*ZIN(MZ+26411)         &
                  +   XIN(MX+28812)*YIN(MY+28812)*ZIN(MZ+28812))          
       END DO
      END DO
      RETURN
!-----------------------------------------------------------------------                                                            
      END                                                               

! ZQOUT                                            
      SUBROUTINE ZQOUT(GHONDO,MAXG)                                          
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION GHONDO(MAXG)                                               
      LOGICAL IANDJ,KANDL,SAME                                          
      COMMON /INTDEX/ IJX(784),IJY(784),IJZ(784),IK(784),               &
                      KLX(784),KLY(784),KLZ(784)              
      COMMON/INTDEX1/IJGT(784),KLGT(784)           
      COMMON /MISC  / IANDJ,KANDL,SAME                                  
      COMMON /SHLNOS/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,               &
                     MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,           &
                     NIJ,IJ,KL                                    
!-----------------------------------------------------------------------                                                                       
!     HONDO Conventional integral Output
!-----------------------------------------------------------------------  
      IJN = 0                                                           
      JMAX = MAXJ                                                       
      DO I = MINI,MAXI                                              
       IF (IANDJ) JMAX = I                                            
       DO 1 J = MINJ,JMAX                                           
          IJN = IJN+1                                                 
          N1 = IJGT(IJN)                                              
          LMAX = MAXL                                                 
          KLN = 0                                                     
        DO K =  MINK,MAXK                                       
         IF (KANDL) LMAX = K                                      
         DO L = MINL,LMAX                                     
          KLN = KLN+1                                           
          IF (SAME .and. KLN > IJN) GO TO 1                
          NN = N1+KLGT(KLN)                                     
          GHONDO(NN) = 0.0d0                                     
         END DO
        END DO
    1  CONTINUE                                                       
      END DO
!-----------------------------------------------------------------------                                                                       
      RETURN                                                            
      END                                                               

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                             ECP integrals                            ! 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

! ECPINT
      SUBROUTINE ECPINT(H0,NBFT,EX,CS,CP,CD,CF,CG,NPRIMI,KSTART,KATOM,  &
                        KNG,KLOC,KMIN,KMAX,NSHELL,Cxyz)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                        
      COMMON/ECPDIM/NCOEF1,NCOEF2,J1LEN,J2LEN,LLIM,NLIM,NTLIM,J4LEN
      COMMON/ECP1/X01,CAX,CAY,CAZ,CA,XCA,YCA,ZCA,                       &
                  X02,BAX,BAY,BAZ,BA,XBA,YBA,ZBA,                       &
                  PHASE,DAX,DAY,DAZ,DA,XDA,YDA,ZDA,XINT,KCNTR  
      !JFHLewYee: Changed NATOMS allowed dimension from 100 to 1000
      COMMON/ECP2/CLP(4004),ZLP(4004),NLP(4004),KFRST(1001,6),          &
                  KLAST(1001,6),LMAX(1001),LPSKIP(1001),IZCORE(1001)
      LOGICAL                         CANDB                            
      COMMON/ECP4/P12(3,2),R12,ACO(3),CANDB                            
      LOGICAL                                  IANDJ,NORM,NORMI,NORMJ    
      COMMON/ECPIDX/Q2,IAMIN,IAMAX,JAMIN,JAMAX,IPMIN,IPMAX,JPMIN,       &
                    JPMAX,KF1,KL1,LLMX,NPC,NPB,IANDJ,NORM,NORMI,NORMJ    
      INTEGER, DIMENSION(8192) :: IA                                     
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB      
      COMMON/MAPSHEL/MAPSHL(600,48),NT   
      
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION FPQR(25,25,25),ZLM(581)
      DIMENSION H0(NBFT),HECP(NBFT),LMF(122),LMX(582),LMY(582),LMZ(582)
      DIMENSION FP(2*11*11*11),G(28*15),COEFI(30),COEFJ(30)
      DIMENSION FQ(1),COEFQ(1)
      LOGICAL CANDA,AANDB
      DIMENSION IANG(35),MI(48)
      DATA IANG/1,3*2,6*3,10*4,15*5/
      INTEGER,ALLOCATABLE,DIMENSION(:)::JFST1,JFST2      
      INTEGER,ALLOCATABLE,DIMENSION(:,:)::LBECP1,LBECP2   
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::DCOEF1,DCOEF2,DCOEF4      
      ALLOCATE(JFST1(J1LEN+9*NCOEF1),JFST2(J2LEN+6*NCOEF2))
      ALLOCATE(LBECP1(9,NCOEF1),LBECP2(6,NCOEF2)) 
      ALLOCATE(DCOEF1(NCOEF1),DCOEF2(NCOEF2),DCOEF4(J4LEN))
!-----------------------------------------------------------------------
      FQ(1)   = 0.0d0
      COEFQ(1)= 0.0d0
      NORM = .TRUE.
      NORMI= NORM
      NORMJ= NORM
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Read ECP
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL ECCODR(DCOEF1,JFST1,LBECP1,DCOEF2,JFST2,LBECP2,              &
                  FPQR,ZLM,LMF,LMX,LMY,LMZ)
      CALL DAWT
      CALL ERRT
      CALL DAWERT
      CALL ECPINI(LMF,LMX,LMY,LMZ)
      CALL ZTAB(ZLM)
      CALL FTAB(FPQR,NLIM-1)
      CALL ECCOD3(FPQR,DCOEF4,ZLM,LMF,LMX,LMY,LMZ)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL VCLR(HECP,1,NBFT)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Triangular Index Matrix
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
      DO I = 1,8192                                                 
       IA(I) = (I*I-I)/2                                               
      ENDDO      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ISHELL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO 1 II=1,NSHELL
       DO IT=1,NT
        ID=MAPSHL(II,IT)
        IF(ID.GT.II) GO TO 1
        MI(IT)=ID
       END DO
       I1 = KSTART(II)
       I2 = I1+KNG(II)-1
       IPMIN = I1
       IPMAX = I2
       ICNTR = KATOM(II)
       IMIN = KMIN(II)
       IMAX = KMAX(II)
       LOCI = KLOC(II)-IMIN
       IAMIN = IMIN
       IAMAX = IMAX
       IIMAX = 1
       IF(IMIN.EQ.1 .AND. IMAX.EQ.4) IIMAX = 2
       DO 2 III=1,IIMAX
        IF(IIMAX.EQ.2) THEN
         IF(III.EQ.1) THEN
          IAMIN = 1
          IAMAX = 1
         ELSE
          IAMIN = 2
          IAMAX = 4
         END IF
        END IF
        NPC0= IANG(IAMAX)
        DO IG=IPMIN,IPMAX
         IF(IAMIN.LE.35) T01= CG(IG)
         IF(IAMIN.LE.20) T01= CF(IG)
         IF(IAMIN.LE.10) T01= CD(IG)
         IF(IAMIN.LE. 4) T01= CP(IG)
         IF(IAMIN.EQ. 1) T01= CS(IG)
         COEFI(IG-IPMIN+1)= T01
        END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       JSHELL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        DO 3 JJ=1,II
         Q2=1
         N2=0
         DO 5 IT=1,NT
          JD=MAPSHL(JJ,IT)
          IF(JD.GT.II) GO TO 3
          ID=MAX0(MI(IT),JD)
          JD=MIN0(MI(IT),JD)
          IF(ID.LT.II) GO TO 5
          IF(JD.LT.JJ) GO TO 5
          IF(JD.GT.JJ) GO TO 3
          N2=N2+1
    5    CONTINUE
         Q2 = NT
         Q2 = Q2/N2
         IANDJ = II.EQ.JJ
         J1 = KSTART(JJ)
         J2 = J1+KNG(JJ)-1
         JPMIN = J1
         JPMAX = J2
         JCNTR = KATOM(JJ)
         JMIN = KMIN(JJ)
         JMAX = KMAX(JJ)
         LOCJ = KLOC(JJ)-JMIN
         JAMIN = JMIN
         JAMAX = JMAX
         JJMAX = 1
         IF(JMIN.EQ.1 .AND. JMAX.EQ.4) JJMAX = 2
         DO 4 JJJ=1,JJMAX
          IF(JJMAX.EQ.2) THEN
           IF(JJJ.EQ.1) THEN
            JAMIN = 1
            JAMAX = 1
           ELSE
            IF(IANDJ .AND. IAMIN.EQ.1) GO TO 4
            JAMIN = 2
            JAMAX = 4
           END IF
          END IF
         DO JG=JPMIN,JPMAX
          IF(JAMIN.LE.35) T01= CG(JG)
          IF(JAMIN.LE.20) T01= CF(JG)
          IF(JAMIN.LE.10) T01= CD(JG)
          IF(JAMIN.LE. 4) T01= CP(JG)
          IF(JAMIN.EQ. 1) T01= CS(JG)
          COEFJ(JG-JPMIN+1)= T01
         END DO
         NPB0= IANG(JAMAX)
         IJMAX= MAX0((IAMAX-IAMIN+1)*(JAMAX-JAMIN+1),60)
         CALL VCLR(G,1,IJMAX)
         CANDB= ICNTR.EQ.JCNTR
         R12= 0.0d0
         DO M=1,3
          P12(M,1)= Cxyz(M,ICNTR)
          P12(M,2)= Cxyz(M,JCNTR)
          R12= R12+(P12(M,2)-P12(M,1))*(P12(M,2)-P12(M,1))
         END DO
         NPC= NPC0
         NPB= NPB0
         NPNP= NPC+NPB-1
         DO 6 IKCNTR=1,NAT
          KCNTR= IKCNTR
          ACO(1)= Cxyz(1,KCNTR)
          ACO(2)= Cxyz(2,KCNTR)
          ACO(3)= Cxyz(3,KCNTR)
          IF(LPSKIP(KCNTR).EQ.1)GO TO 6
          LLMX = LMAX(KCNTR)+1
          KF1 = KFRST(KCNTR,1)
          KL1 = KLAST(KCNTR,1)
         CANDA= ICNTR.EQ.KCNTR
         AANDB= KCNTR.EQ.JCNTR
         CALL ECPCBA(CANDA,AANDB,ICAB,IPOW)
         IF(ICNTR.NE.KCNTR .OR. KCNTR.NE.JCNTR)                         &
          CALL ECPPWR(IPOW,NPC0,NPB0)                                    
           IF(ICAB.EQ.1)THEN                                             
            CALL ECPAA1(NPNP,FPQR,COEFI,COEFJ,DCOEF4,G,EX,NPRIMI)                  
           ELSE IF(ICAB.EQ.2 .OR. ICAB.EQ.3) THEN                        
            CALL ECPRA2(ICAB,NPNP,FP,COEFI,COEFJ,                       &
                        DCOEF1,JFST1,LBECP1,DCOEF4,                     &
                        DCOEF2,JFST2,LBECP2,G,                          &
                        ZLM,LMF,LMX,LMY,LMZ,EX,NPRIMI)                             
           ELSE IF(ICAB.EQ.4) THEN                                       
            IC4C= 0                                                      
            CALL ECPDRA(IC4C,NPNP,FP,FQ,COEFI,COEFQ,COEFJ,DCOEF1,       &
                        JFST1,LBECP1,DCOEF2,JFST2,LBECP2,G,ZLM,LMF,     &
                        LMX,LMY,LMZ,EX,NPRIMI)
           END IF
    6     CONTINUE
          MMAX= JAMAX
          NN= 1
          DO I=IAMIN,IAMAX
           IN= IA(LOCI+I)+LOCJ
           IF(IANDJ) MMAX= I
           DO J=JAMIN,MMAX
            HECP(IN+J)= G(NN)
            NN= NN+1
           END DO
          END DO
    4    CONTINUE
    3   CONTINUE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       JSHELL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    2  CONTINUE
    1 CONTINUE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     New Core Hamiltonian
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      H0 = H0 + HECP
!-----------------------------------------------------------------------
      DEALLOCATE(LBECP1,LBECP2,JFST1,JFST2,DCOEF1,DCOEF4,DCOEF2)       
      RETURN
      END

! ECCODR
      SUBROUTINE ECCODR(DCOEF1,JFST1,LBECP1,DCOEF2,JFST2,LBECP2,        &
                        FPQR,ZLM,LMF,LMX,LMY,LMZ)                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION FPQR(25,25,25),ZLM(581),DCOEF1(*),DCOEF2(*)
      DIMENSION JFST1(*),LBECP1(9,*),JFST2(*),LBECP2(6,*),              &
                LMF(122),LMZ(581),LMX(581),LMY(581)
      COMMON /ECPDIM/ NCOEF1,NCOEF2,J1LEN,J2LEN,LLIM,NLIM,NTLIM,J4LEN
      CALL ECPINI(LMF,LMX,LMY,LMZ)
      CALL ZTAB(ZLM)
      CALL FTAB(FPQR,NLIM-1)
      CALL ECCOD1(DCOEF1,JFST1,LBECP1,FPQR,ZLM,LMF,LMX,LMY,LMZ)
      CALL ECCOD2(DCOEF2,JFST2,LBECP2,FPQR,ZLM,LMF,LMX,LMY,LMZ)
      RETURN
      END

! ECCOD1
      SUBROUTINE ECCOD1(DCOEF1,JFST1,LBECP1,FPQR,ZLM,LMF,LMX,LMY,LMZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DCOEF1(*),FPQR(25,25,25),ZLM(*)
      DIMENSION JFST1(*),LBECP1(9,*),LMF(*),LMX(*),LMY(*),LMZ(*)
      COMMON /ECPDIM/ NCOEF1,NCOEF2,J1LEN,J2LEN,LLIM,NLIM,NTLIM,J4LEN
      COMMON /GBASE / NFST(8),NX(84),NY(84),NZ(84)
      PARAMETER(SQRFPI=3.5449077018110D+00, TOL=1.0D-10)
!        -----  THREE-CENTER ONE ELECTRON INTEGRALS.            -----
      DIMENSION BINCO(28),IA(128)
      DATA BINCO/           1.0D+00,                                    &
                        1.0D+00, 1.0D+00,                               &
                     1.0D+00, 2.0D+00, 1.0D+00,                         &
                  1.0D+00, 3.0D+00, 3.0D+00, 1.0D+00,                   &
               1.0D+00, 4.0D+00, 6.0D+00, 4.0D+00, 1.0D+00,             &
            1.0D+00, 5.0D+00,10.0D+00,10.0D+00, 5.0D+00,1.0D+00,        &
         1.0D+00,6.0D+00,15.0D+00,20.0D+00,15.0D+00,6.0D+00,1.0D+00/

      J= 0
      DO 110 I=1,128
       IA(I)= J
  110 J= J+I

      JNDX = 0
      JJLST= 0
! LOOP OVER POSSIBLE I SHELLS
      DO 230 NN1=1,NLIM
         NF1 = NFST(NN1)
         NL1 = NFST(NN1+1)-1
         NN2LIM= NN1
! LOOP OVER THE I SHELL ANGUALR MOMENTUM FUNCTIONS
         DO 220 N1=NF1,NL1
            N1T = IA(N1)
            MX1 = IA(NX(N1)+1)+1
            MY1 = IA(NY(N1)+1)+1
            MZ1 = IA(NZ(N1)+1)+1
! LOOP OVER THE POSSIBLE J SHELLS (NOTE J SHELL <= I SHELL)
            DO 210 NN2=1,NN2LIM
               NF2 = NFST(NN2)
               NL2 = MIN0(NFST(NN2+1)-1,N1)
! LOOP OVER J SHELL ANGULAR MOMENTUM
               DO 200 N2=NF2,NL2
                  INDX = N2+N1T
                  MX2 = IA(NX(N2)+1)+1
                  MY2 = IA(NY(N2)+1)+1
                  MZ2 = IA(NZ(N2)+1)+1
                  LLMAX= NN1+NN2-2
! LOOP OVER LA, MU, KX...
                  DO 190 LA=0,LLMAX
                     DO 180 MU=-LA,LA
                        DO 170 KX=0,NX(N1)
                           DO 160 KY=0,NY(N1)
                              DO 150 KZ=0,NZ(N1)
                                 DO 140 KXP=0,NX(N2)
                                    DO 130 KYP=0,NY(N2)
                                       DO 120 KZP=0,NZ(N2)
                              KA = KX+KY+KZ+KXP+KYP+KZP
!
! NOTE NO NONZERO INTEGRALS FOR LA > KA
! IF KA + LA IS EVEN GO AHEAD AND CALCULATE THE ANGULAR INTEGRAL
!
                              IF(MOD(LA+KA,2).NE.1 .AND. LA.LE.KA) THEN
      DC=DCO(0,0,KX+KXP,KY+KYP,KZ+KZP,LA,MU,FPQR,ZLM,LMF,LMX,LMY,LMZ)
!
! IF THE INTEGRAL IS NON-ZERO THEN STORE IT AWAY FOR FUTURE USE
!
                                 IF(ABS(DC).LT.TOL) GO TO 120
                                 JNDX = JNDX+1
                                 LBECP1(1,JNDX)= KA
                                 LBECP1(2,JNDX)= LA
                                 LBECP1(3,JNDX)= MU
                                 LBECP1(4,JNDX)= KX
                                 LBECP1(5,JNDX)= KY
                                 LBECP1(6,JNDX)= KZ
                                 LBECP1(7,JNDX)= KXP
                                 LBECP1(8,JNDX)= KYP
                                 LBECP1(9,JNDX)= KZP
!
! THE INTEGRAL IS MULTIPLIED BY THE COMBINATION FACTORS AND BY SQRT(4PI)
! SINCE THERE IS AN EXTRA 1/SQRT(4PI) IN THE INTEGRAL
!
      DCOEF1(JNDX)= BINCO(MX1+KX) *BINCO(MY1+KY) *BINCO(MZ1+KZ) *       &
                    BINCO(MX2+KXP)*BINCO(MY2+KYP)*BINCO(MZ2+KZP)*       &
                    DC*SQRFPI
                              END IF
  120                                  CONTINUE
  130                               CONTINUE
  140                            CONTINUE
  150                         CONTINUE
  160                      CONTINUE
  170                   CONTINUE
  180                CONTINUE
  190             CONTINUE
                  JFST1(INDX)  = JJLST+1
                  JFST1(INDX+1)= JNDX+1
                  JJLST= JNDX
  200          CONTINUE
  210       CONTINUE
  220    CONTINUE
  230 CONTINUE

      IF(JNDX.GT.NCOEF1) THEN
         WRITE(6,9010) JNDX,NCOEF1
         CALL ABRT
      END IF
      NCOEF1= JNDX
      RETURN
 9010 FORMAT(/'****ECCOD1 OUT OF BOUNDS, USED ',I7,' ALLOWED ',I7)
      END

! ECCOD2
      SUBROUTINE ECCOD2(DCOEF2,JFST2,LBECP2,FPQR,ZLM,LMF,LMX,LMY,LMZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DCOEF2(*),FPQR(25,25,25),ZLM(*)
      DIMENSION JFST2(*),LBECP2(6,*),LMF(*),LMX(*),LMY(*),LMZ(*)
      COMMON /ECPDIM/ NCOEF1,NCOEF2,J1LEN,J2LEN,LLIM,NLIM,NTLIM,J4LEN
      COMMON /GBASE / NFST(8),NX(84),NY(84),NZ(84)
      PARAMETER(TOL=1.0D-10)

! BINCO -- INDEXED COMBINATIONS (B!/(A!(B-A)!)) INDEXED BY B AND A
!          AS:  (B,A): (0,0),(1,0),(1,1),(2,0),(2,1),(2,2)...
! AKA PASCAL'S TRIANGLE, ENTRIES THROUGH I FUNCTIONS (G HESSIANS)

      DIMENSION BINCO(28),IA(8)
      DATA BINCO/            1.0D+00,                                   &
                        1.0D+00, 1.0D+00,                               &
                    1.0D+00, 2.0D+00, 1.0D+00,                          &
                 1.0D+00, 3.0D+00, 3.0D+00, 1.0D+00,                    &
              1.0D+00, 4.0D+00, 6.0D+00, 4.0D+00, 1.0D+00,              &
           1.0D+00, 5.0D+00,10.0D+00,10.0D+00, 5.0D+00,1.0D+00,         &
        1.0D+00,6.0D+00,15.0D+00,20.0D+00,15.0D+00,6.0D+00,1.0D+00/

!        -----  ROUTINE FINDS FORMULA CODE FOR ONE-ELECTRON     -----
!        -----  THREE-CENTER INTEGRALS INVOLVING PROJECTION     -----
!        -----  OPERATORS.                                      -----

! THE TYPE 2 ANGULAR INTEGRAL TABLE GENERATOR

      J= 0
      DO 110 I=1,8
!        IA(I)=(I*(I-1))/2
         IA(I)= J
  110 J= J+I

      JNDX = 0
      JJLST= 0
      LLLIM= LLIM-1
      DO 200 L=0,LLLIM
         DO 190 M=L,-L,-1
            LMINDX =(L*(L+1)-M)*NTLIM
            DO 180 NN=1,NLIM
               NF = NFST(NN)
               NL = NFST(NN+1)-1
               L2MX = L+NN-1
               DO 170 N=NF,NL
                  INDX = LMINDX+N
                  MX = IA(NX(N)+1)+1
                  MY = IA(NY(N)+1)+1
                  MZ = IA(NZ(N)+1)+1
                  DO 160 LA=0,L2MX
                     DO 150 MU=-LA,LA
                        DO 140 KX=0,NX(N)
                           DO 130 KY=0,NY(N)
                              DO 120 KZ=0,NZ(N)
                                 IS = KX+KY+KZ
      DC= DCO(L,M,KX,KY,KZ,LA,MU,FPQR,ZLM,LMF,LMX,LMY,LMZ)
                                 IF(ABS(DC).LT.TOL) GO TO 120
                                 JNDX = JNDX+1
                                 LBECP2(1,JNDX)= LA
                                 LBECP2(2,JNDX)= IS
                                 LBECP2(3,JNDX)= MU
                                 LBECP2(4,JNDX)= KX
                                 LBECP2(5,JNDX)= KY
                                 LBECP2(6,JNDX)= KZ
      DCOEF2(JNDX)= BINCO(MX+KX)*BINCO(MY+KY)*BINCO(MZ+KZ)*DC
  120                         CONTINUE
  130                      CONTINUE
  140                   CONTINUE
  150                CONTINUE
  160             CONTINUE
                  JFST2(INDX)  = JJLST+1
                  JFST2(INDX+1)= JNDX+1
                  JJLST= JNDX
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE

      IF(JNDX.GT.NCOEF2) THEN
         WRITE(6,9010) JNDX,NCOEF2
         CALL ABRT
      END IF
      NCOEF2= JNDX
      RETURN
 9010 FORMAT(/'****ECCOD2 OUT OF BOUNDS, USED ',I6,' ALLOWED ',I5)
      END

! ECCOD3
      SUBROUTINE ECCOD3(FPQR,DCOEF4,ZLM,LMF,LMX,LMY,LMZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ECPDIM/NCOEF1,NCOEF2,J1LEN,J2LEN,LLIM,NLIM,NTLIM,J4LEN
      COMMON/GBASE/NFST(8),NX(84),NY(84),NZ(84)
      DOUBLE PRECISION FPQR(25,25,25),ZLM(*),DCOEF4(*)
      DIMENSION LMF(*),LMX(*),LMY(*),LMZ(*)
!-----------------------------------------------------------------------
      SQRFPI = 3.5449077018110D+00
      IIIDX = 0
      DO L=0,LLIM-1
       DO M=L,-L,-1
        LMINDX =(L*(L+1)-M)*NTLIM
        DO NN=1,NLIM
         NF = NFST(NN)
         NL = NFST(NN+1)-1
         DO N=NF,NL
          INDX = LMINDX+N
          IIIDX= MAX0(INDX,IIIDX)
          DCOEF4(INDX)= DCO(0,0,NX(N),NY(N),NZ(N),L,M,                  &
                            FPQR,ZLM,LMF,LMX,LMY,LMZ)*SQRFPI
         END DO
        END DO
       END DO
      END DO
      IF(IIIDX>J4LEN) THEN
       WRITE(6,*)'Stop: IIIDX > J4LEN'
       CALL ABRT
      END IF
!-----------------------------------------------------------------------
      RETURN
      END

! FTAB
      SUBROUTINE FTAB(FPQR,NMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER P,Q,R
      DOUBLE PRECISION FPQR(25,25,25)
!     FOR NMAX= 6, N4MAX= 4*NMAX+1=25, AND 3*N4MAX-2=73
      DIMENSION REC(73)

!        -----  ROUTINE SETS A TABLE OF F-FUNCTION VALUES.      -----
!     WHERE F MEANS THE ANGULAR INTEGRAL GIVEN BY:
!     INT[DW X^R Y^S Z^T = O FOR R, S, OR T ODD
!                        =(R-1)!!(S-1)!!(T-1)!!/(R+S+T+1)!!
!                      FOR R, S, AND T EVEN
!     NEEDS DIMENSION 4*(MAX REAL ANGULAR MOMENTUM) + 1
!     WHERE 1 IS ADDED SINCE WE REALLY NEED THE 0-TH ELEMENT
!     THUS FOR F GRADIENTS NEED G(4) SO 4*4+1=17
      PARAMETER (FPI=12.566370614359D+00)
!-----------------------------------------------------------------------
      N4MAX = 4*NMAX+1
      IF(N4MAX.GT.25) THEN
       WRITE(6,9010) NMAX
       CALL ABRT
      END IF

      DEN= 1.0d0+2.0d0
      DO P=3,3*N4MAX-2,2
       REC(P)= 1.0d0/DEN
       DEN= DEN+2.0d0
      END DO

!        -----  ZERO OUT THE TABLE.                             -----
      DO 120 P=1,N4MAX
        DO 120 Q=1,N4MAX
          DO 120 R=1,N4MAX
  120 FPQR(R,Q,P) = 0.0D+00
!        -----  RECURSIVELY GENERATE NON-ZERO ENTRIES.          -----
      FPQR(1,1,1) = FPI
      PM2=-1.0d0
      DO 150 P=1,N4MAX,2
         QM2=-1.0d0
         DO 140 Q=1,N4MAX,2
            RM2=-1.0d0
            DO 130 R=1,N4MAX,2
               IF(P.GT.1) THEN
                  FPQR(P,Q,R)= PM2*FPQR(P-2,Q,R)*REC(P-2+Q+R)
               ELSE IF(Q.GT.1) THEN
                  FPQR(P,Q,R)= QM2*FPQR(P,Q-2,R)*REC(P+Q-2+R)
               ELSE IF(R.GT.1) THEN
                  FPQR(P,Q,R)= RM2*FPQR(P,Q,R-2)*REC(P+Q+R-2)
               END IF
  130       RM2= RM2+2.0d0
  140    QM2= QM2+2.0d0
  150 PM2= PM2+2.0d0

      RETURN
 9010 FORMAT(' FTAB: NMAX TOO LARGE =',I3)
      END

! ZTAB
      SUBROUTINE ZTAB(ZLM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION ZLM(*)
!
!        -----  ROUTINE SETS UP THE REAL SPHERICAL HARMONICS    -----
!        -----  IN THE FORM OF LINEAR COMBINATIONS OF           -----
!        -----  CARTESIAN PRODUCTS.                             -----
!
      PARAMETER (THR=3.0D+00, FOR=4.0D+00,                              &
                 FIV=5.0D+00, SIX=6.0D+00, SEV=7.0D+00, EIT=8.0D+00,    &
                 F09=9.0D+00, TEN=10.0D+00, F11=11.0D+00, F13=13.0D+00, &
                 F15=15.0D+00, F17=17.0D+00, F19=19.0D+00)
      PARAMETER (ZP5=0.5D+00, FPI=12.566370614359D+00)
      DIMENSION FAP(0:17),FAT(0:17)

      FAP(0)= SQRT(1.0d0/FPI)
      FAP(1)= FAP(0)/ DSQRT(2.0d0)
      DO I=2,16,2
       FAP(I  )= FAP(I-2)*ZP5
       FAP(I+1)= FAP(I-1)*ZP5
      END DO
      DO 120 I=0,17
  120 FAT(I)= FAP(I)*THR

      SR3= SQRT(THR)
      SR5= SQRT(FIV)
      SR7= SQRT(SEV)
      S11= SQRT(F11)
      S13= SQRT(F13)
      S15= SQRT(F15)
      S17= SQRT(F17)
      S19= SQRT(F19)

      S3BYS5= SR3*SR5
      S3BYS7= SR3*SR7
      S5BYS7= SR5*SR7
      S5XS11= SR5*S11
      S7XS11= SR7*S11
      S7XS13= SR7*S13
      SR0119= SQRT(119.0D+00)
      SR0221= SQRT(221.0D+00)
      SR0247= SQRT(247.0D+00)
      SR0323= SQRT(323.0D+00)
      SR0385= S5BYS7*S11
      SR1001= S7XS11*S13
      SR3003= SR1001*SR3
      SR5005= SR1001*SR5
      SR1309= S7XS11*S17
      SR1365= S3BYS5*S7XS13
      SR2431= SR0221*S11
      S13585= SR0247*S5XS11
! L= 0, ML= 0
      ZLM(  1)= FAP(0)
! L= 1, ML=+1,0,-1
!     POWERS (X, Z, Y)
      ZLM(  2)= SR3*FAP(0)
      ZLM(  3)= ZLM(  2)
      ZLM(  4)= ZLM(  2)
! ---
! L= 2, ML=+2...-2
!     POWERS ([5]X^2+[6]Y^2,[7]XZ,[8]Z^2+[9],[10]YZ,[11]XY)
      ZLM(  5)= S15*FAP(2)
      ZLM(  6)=-ZLM(  5)
      ZLM(  7)= ZLM(  5)*2.0d0
      TEMP    = SR5*FAP(2)
      ZLM(  8)= TEMP*THR
      ZLM(  9)=-TEMP
      ZLM( 10)= ZLM(  7)
      ZLM( 11)= ZLM(  7)
! ---
! L= 3, ML=+3...-3
!     POWERS ([12]X^3+[13]XY^2,[14]X^2Z+[15]Y^2Z,[16]XZ^2+[17]X,
!             [18]Z^3+[19]Z,[20]YZ^2+[21]Y,[22]XYZ,[23]X^2Y+[24]Y^3)
      ZLM( 12)= S5BYS7*FAP(3)
      ZLM( 13)=-ZLM( 12)*THR
      ZLM( 14)= S3BYS7*SR5*FAP(2)
      ZLM( 15)=-ZLM( 14)
      TEMP    = S3BYS7*FAP(3)
      ZLM( 16)= TEMP*FIV
      ZLM( 17)=-TEMP
      TEMP    = SR7*FAP(2)
      ZLM( 18)= TEMP*FIV
      ZLM( 19)=-TEMP*THR
      ZLM( 20)= ZLM( 16)
      ZLM( 21)= ZLM( 17)
      ZLM( 22)= ZLM( 14)*2.0d0
      ZLM( 23)=-ZLM( 13)
      ZLM( 24)=-ZLM( 12)
! ---
! L= 4, ML=+4...-4
      ZLM( 25)= S5BYS7*FAT(6)
      ZLM( 26)=-ZLM( 25)*SIX
      ZLM( 27)= ZLM( 25)
      ZLM( 28)= S5BYS7*FAT(3)
      ZLM( 29)=-ZLM( 28)*THR
      TEMP    = SR5*FAT(4)
      ZLM( 30)= TEMP*SEV
      ZLM( 31)=-ZLM( 30)
      ZLM( 32)= TEMP
      ZLM( 33)=-TEMP
      TEMP    = SR5*FAT(3)
      ZLM( 34)= TEMP*SEV
      ZLM( 35)=-TEMP*THR
      TEMP    = FAT(6)
      ZLM( 36)= TEMP*FIV*SEV
      ZLM( 37)=-TEMP*FIV*SIX
      ZLM( 38)= TEMP*THR
      ZLM( 39)= ZLM( 34)
      ZLM( 40)= ZLM( 35)
      TEMP    = SR5*FAT(2)
      ZLM( 41)= TEMP*SEV
      ZLM( 42)=-TEMP
      ZLM( 43)=-ZLM( 29)
      ZLM( 44)=-ZLM( 28)
      ZLM( 45)= S5BYS7*FAT(2)
      ZLM( 46)=-ZLM( 45)
! ---
! L= 5, ML=+5...-5
      ZLM( 47)= S7XS11*FAT(7)
      ZLM( 48)=-ZLM( 47)*TEN
      ZLM( 49)= ZLM( 47)*FIV
      ZLM( 50)= SR0385*FAT(6)
      ZLM( 51)=-ZLM( 50)*SIX
      ZLM( 52)= ZLM( 50)
      TEMP    = SR0385*FAP(7)
      ZLM( 53)= TEMP*THR*THR
      ZLM( 54)=-TEMP*THR*F09
      ZLM( 55)= TEMP*THR
      ZLM( 56)=-TEMP
      TEMP    = SR0385*SR3*FAP(4)
      ZLM( 57)= TEMP*THR
      ZLM( 58)=-ZLM( 57)
      ZLM( 59)=-TEMP
      ZLM( 60)= TEMP
      TEMP    = S3BYS5*S11*FAP(6)
      ZLM( 61)= TEMP*21.0D+00
      ZLM( 62)=-TEMP*14.0D+00
      ZLM( 63)= TEMP
      TEMP    = S11*FAP(6)
      ZLM( 64)= TEMP*63.0D+00
      ZLM( 65)=-TEMP*70.0D+00
      ZLM( 66)= TEMP*15.0D+00
      ZLM( 67)= ZLM( 61)
      ZLM( 68)= ZLM( 62)
      ZLM( 69)= ZLM( 63)
      TEMP    = SR0385*SR3*FAP(2)
      ZLM( 70)= TEMP*THR
      ZLM( 71)=-TEMP
      ZLM( 72)=-ZLM( 54)
      ZLM( 73)=-ZLM( 55)
      ZLM( 74)=-ZLM( 53)
      ZLM( 75)=-ZLM( 56)
      ZLM( 76)= SR0385*FAT(2)
      ZLM( 77)=-ZLM( 76)
      ZLM( 78)= ZLM( 49)
      ZLM( 79)= ZLM( 48)
      ZLM( 80)= ZLM( 47)
! ---
! L= 6, ML=+6...-6
      ZLM( 81)= SR3003*FAP(9)
      ZLM( 82)=-ZLM( 81)*F15
      ZLM( 83)=-ZLM( 82)
      ZLM( 84)=-ZLM( 81)
      ZLM( 85)= SR1001*FAT(7)
      ZLM( 86)=-ZLM( 85)*TEN
      ZLM( 87)= ZLM( 85)*FIV
      TEMP    = S7XS13*FAT(8)
      ZLM( 88)= TEMP*F11
      ZLM( 89)=-TEMP*F11*SIX
      ZLM( 90)= ZLM( 88)
      ZLM( 91)=-TEMP
      ZLM( 92)= TEMP*SIX
      ZLM( 93)=-TEMP
      TEMP    = SR1365*FAP(7)
      ZLM( 94)= TEMP*F11
      ZLM( 95)=-TEMP*THR*F11
      ZLM( 96)= TEMP*THR*THR
      ZLM( 97)=-TEMP*THR
      TEMP    = SR1365*FAP(9)
      ZLM( 98)= TEMP*THR*F11
      ZLM( 99)=-ZLM( 98)
      ZLM(100)=-TEMP*THR*SIX
      ZLM(101)= TEMP*THR*SIX
      ZLM(102)= TEMP
      ZLM(103)=-TEMP
      TEMP    = S3BYS7*S13*FAP(6)
      ZLM(104)= TEMP*THR*F11
      ZLM(105)=-TEMP*THR*TEN
      ZLM(106)= TEMP*FIV
      TEMP    = S13*FAP(8)
      ZLM(107)= TEMP*231.0D+00
      ZLM(108)=-TEMP*315.0D+00
      ZLM(109)= TEMP*105.0D+00
      ZLM(110)=-TEMP*FIV
      ZLM(111)= ZLM(104)
      ZLM(112)= ZLM(105)
      ZLM(113)= ZLM(106)
      TEMP    = SR1365*FAP(7)
      ZLM(114)= TEMP*THR*F11
      ZLM(115)=-TEMP*THR*SIX
      ZLM(116)= TEMP
      ZLM(117)=-ZLM( 95)
      ZLM(118)=-ZLM( 94)
      ZLM(119)=-ZLM( 96)
      ZLM(120)=-ZLM( 97)
      TEMP    = S7XS13*FAT(4)
      ZLM(121)= TEMP*F11
      ZLM(122)=-ZLM(121)
      ZLM(123)=-TEMP
      ZLM(124)= TEMP
      ZLM(125)= ZLM( 87)
      ZLM(126)= ZLM( 86)
      ZLM(127)= ZLM( 85)
      TEMP    = SR3003*FAP(9)
      ZLM(128)= TEMP*SIX
      ZLM(129)=-TEMP*20.0D+00
      ZLM(130)= ZLM(128)
! ---
! L= 7, ML=-7
      ZLM(131)= S5XS11*S13*FAT(10)
      ZLM(132)= ZLM(131)*SEV*FIV
      ZLM(133)=-ZLM(131)*SEV*THR
      ZLM(134)=-ZLM(131)*SEV
! L= 7, ML=-6
      ZLM(135)= SR5005*FAT(9)
      ZLM(136)=-ZLM(135)*F15
      ZLM(137)=-ZLM(135)
      ZLM(138)=-ZLM(136)
! L= 7, ML=-5
      TEMP    = SR0385*FAT(10)
      ZLM(139)= TEMP*F13
      ZLM(140)=-ZLM(139)*TEN
      ZLM(141)= ZLM(139)*FIV
      ZLM(142)=-TEMP
      ZLM(143)=-ZLM(142)*TEN
      ZLM(144)= ZLM(142)*FIV
! L= 7, ML=-4
      TEMP    = SR0385*FAT(8)
      ZLM(145)= TEMP*F13
      ZLM(146)=-TEMP*F13*SIX
      ZLM(147)= ZLM(145)
      ZLM(148)=-TEMP*THR
      ZLM(149)= TEMP*THR*SIX
      ZLM(150)= ZLM(148)
! L= 7, ML=-3
      TEMP    = S5BYS7*FAT(10)
      ZLM(151)= TEMP*F11*F13
      ZLM(152)=-ZLM(151)*THR
      ZLM(153)=-TEMP*F11*SIX
      ZLM(154)=-ZLM(153)*THR
      ZLM(155)= TEMP*THR
      ZLM(156)=-TEMP*F09
! L= 7, ML=-2
      TEMP    = S5BYS7*FAT(9)
      ZLM(157)= TEMP*F11*F13
      ZLM(158)=-ZLM(157)
      ZLM(159)=-TEMP*F11*TEN
      ZLM(160)=-ZLM(159)
      ZLM(161)= TEMP*F15
      ZLM(162)=-ZLM(161)
! L= 7, ML=-1
      TEMP    = SR7*S15*FAP(10)
      ZLM(163)= TEMP*429.0D+00
      ZLM(164)=-TEMP*495.0D+00
      ZLM(165)= TEMP*135.0D+00
      ZLM(166)=-TEMP*FIV
! L= 7, ML= 0
      TEMP    = S15*FAP(8)
      ZLM(167)= TEMP*429.0D+00
      ZLM(168)=-TEMP*693.0D+00
      ZLM(169)= TEMP*315.0D+00
      ZLM(170)=-TEMP*35.0D+00
! L= 7, ML= 1
      ZLM(171)= ZLM(163)
      ZLM(172)= ZLM(164)
      ZLM(173)= ZLM(165)
      ZLM(174)= ZLM(166)
! L= 7, ML= 2
      ZLM(175)= ZLM(157)*2.0d0
      ZLM(176)= ZLM(159)*2.0d0
      ZLM(177)= ZLM(161)*2.0d0
! L= 7, ML= 3
      ZLM(178)=-ZLM(152)
      ZLM(179)=-ZLM(151)
      ZLM(180)=-ZLM(154)
      ZLM(181)=-ZLM(153)
      ZLM(182)=-ZLM(156)
      ZLM(183)=-ZLM(155)
! L= 7, ML= 4
      ZLM(184)= ZLM(145)*FOR
      ZLM(185)=-ZLM(184)
      ZLM(186)= ZLM(148)*FOR
      ZLM(187)=-ZLM(186)
! L= 7, ML= 5
      ZLM(188)= ZLM(141)
      ZLM(189)= ZLM(140)
      ZLM(190)= ZLM(139)
      ZLM(191)= ZLM(144)
      ZLM(192)= ZLM(143)
      ZLM(193)= ZLM(142)
! L= 7, ML= 6
      TEMP    = SR5005*FAT(7)
      ZLM(194)= TEMP*THR
      ZLM(195)=-TEMP*TEN
      ZLM(196)= ZLM(194)
! L= 7, ML= 7
      ZLM(197)=-ZLM(131)
      ZLM(198)=-ZLM(132)
      ZLM(199)=-ZLM(133)
      ZLM(200)=-ZLM(134)
! ---
! L= 8, ML=-8
      ZLM(201)= SR2431*SR5*FAT(14)
      ZLM(202)=-ZLM(201)*2.8D+01
      ZLM(203)= ZLM(201)*7.0D+01
      ZLM(204)= ZLM(202)
      ZLM(205)= ZLM(201)
! L= 8, ML=-7
      ZLM(206)= SR2431*SR5*FAT(10)
      ZLM(207)=-ZLM(206)*SEV*THR
      ZLM(208)= ZLM(206)*SEV*FIV
      ZLM(209)=-ZLM(206)*SEV
! L= 8, ML=-6
      TEMP    = SR2431*SR3*FAP(11)
      ZLM(210)= TEMP*F15
      ZLM(211)=-ZLM(210)*F15
      ZLM(212)=-ZLM(211)
      ZLM(213)=-ZLM(210)
      ZLM(214)=-TEMP
      ZLM(215)= ZLM(210)
      ZLM(216)=-ZLM(215)
      ZLM(217)= TEMP
! L= 8, ML=-5
      TEMP    = SR1001*S17*FAT(10)
      ZLM(218)= TEMP*FIV
      ZLM(219)=-ZLM(218)*TEN
      ZLM(220)= ZLM(218)*FIV
      ZLM(221)=-TEMP
      ZLM(222)= TEMP*TEN
      ZLM(223)= ZLM(221)*FIV
! L= 8, ML=-4
      TEMP    = SR1309*FAT(12)
      ZLM(224)= TEMP*F13*FIV
      ZLM(225)= ZLM(224)
      ZLM(226)=-ZLM(224)*SIX
      ZLM(227)=-TEMP*F13*2.0d0
      ZLM(228)= ZLM(227)
      ZLM(229)=-ZLM(227)*SIX
      ZLM(230)= TEMP
      ZLM(231)= TEMP
      ZLM(232)=-TEMP*SIX
! L= 8, ML=-3
      TEMP    = SR1309*S3BYS5*FAP(10)
      ZLM(233)= TEMP*F13*THR
      ZLM(234)=-ZLM(233)*THR
      ZLM(235)=-TEMP*F13*2.0d0
      ZLM(236)=-ZLM(235)*THR
      ZLM(237)= TEMP*THR
      ZLM(238)=-ZLM(237)*THR
! L= 8, ML=-2
      TEMP    = SR0119*SR5*FAT(11)
      ZLM(239)= TEMP*F11*F13
      ZLM(240)=-ZLM(239)
      ZLM(241)= ZLM(240)
      ZLM(242)= ZLM(239)
      ZLM(243)= TEMP*F11*THR
      ZLM(244)=-ZLM(243)
      ZLM(245)=-TEMP
      ZLM(246)= TEMP
! L= 8, ML=-1
      TEMP    = S17*FAT(10)
      ZLM(247)= TEMP*7.15D+02
      ZLM(248)=-TEMP*1.001D+03
      ZLM(249)= TEMP*3.85D+02
      ZLM(250)=-TEMP*3.5D+01
! L= 8, ML= 0
      TEMP    = S17*FAP(14)
      ZLM(251)= TEMP*6.435D+03
      ZLM(252)=-TEMP*1.2012D+04
      ZLM(253)= TEMP*6.930D+03
      ZLM(254)=-TEMP*1.260D+03
      ZLM(255)= TEMP*3.5D+01
! L= 8, ML= 1
      ZLM(256)= ZLM(247)
      ZLM(257)= ZLM(248)
      ZLM(258)= ZLM(249)
      ZLM(259)= ZLM(250)
! L= 8, ML= 2
      TEMP    = SR0119*SR5*FAT(9)
      ZLM(260)= TEMP*F11*F13
      ZLM(261)=-ZLM(260)
      ZLM(262)= TEMP*F11*THR
      ZLM(263)=-TEMP
! L= 8, ML= 3
      ZLM(264)=-ZLM(234)
      ZLM(265)=-ZLM(233)
      ZLM(266)=-ZLM(236)
      ZLM(267)=-ZLM(235)
      ZLM(268)=-ZLM(238)
      ZLM(269)=-ZLM(237)
! L= 8, ML= 4
      TEMP    = SR1309*FAT(8)
      ZLM(270)= TEMP*F13*FIV
      ZLM(271)=-ZLM(270)
      ZLM(272)=-TEMP*F13*2.0d0
      ZLM(273)=-ZLM(272)
      ZLM(274)= TEMP
      ZLM(275)=-TEMP
! L= 8, ML= 5
      TEMP    = SR1001*S17*FAT(10)
      ZLM(276)= TEMP*FIV*FIV
      ZLM(277)=-TEMP*FIV*TEN
      ZLM(278)= TEMP*FIV
      ZLM(279)=-TEMP*FIV
      ZLM(280)= TEMP*TEN
      ZLM(281)=-TEMP
! L= 8, ML= 6
      TEMP    = SR2431*SR3*FAP(9)
      ZLM(282)= TEMP*F15*THR
      ZLM(283)=-TEMP*F15*TEN
      ZLM(284)= ZLM(282)
      ZLM(285)=-TEMP*THR
      ZLM(286)= TEMP*TEN
      ZLM(287)= ZLM(285)
! L= 8, ML= 7
      ZLM(288)=-ZLM(209)
      ZLM(289)=-ZLM(208)
      ZLM(290)=-ZLM(207)
      ZLM(291)=-ZLM(206)
! L= 8, ML= 8
      ZLM(292)= ZLM(201)*EIT
      ZLM(293)=-ZLM(292)*SEV
      ZLM(294)=-ZLM(293)
      ZLM(295)=-ZLM(292)
! ---
! L= 9, ML=-9
      ZLM(296)= S13585*S17*FAP(15)
      ZLM(297)=-ZLM(296)*3.6D+01
      ZLM(298)= ZLM(296)*1.26D+02
      ZLM(299)=-ZLM(296)*8.4D+01
      ZLM(300)= ZLM(296)*9.0D+00
! L= 9, ML=-8
      ZLM(301)= S13585*S17*FAT(14)
      ZLM(302)=-ZLM(301)*SEV*FOR
      ZLM(303)= ZLM(301)*SEV*TEN
      ZLM(304)= ZLM(302)
      ZLM(305)= ZLM(301)
! L= 9, ML=-7
      TEMP    = S13585*FAT(15)
      ZLM(306)= TEMP*F17
      ZLM(307)=-ZLM(306)*SEV*THR
      ZLM(308)= ZLM(306)*SEV*FIV
      ZLM(309)=-ZLM(306)*SEV
      ZLM(310)=-TEMP
      ZLM(311)= TEMP*SEV*THR
      ZLM(312)=-TEMP*SEV*FIV
      ZLM(313)= TEMP*SEV
! L= 9, ML=-6
      TEMP    = SR0247*S11*S15*FAP(11)
      ZLM(314)= TEMP*F17
      ZLM(315)=-ZLM(314)*F15
      ZLM(316)=-ZLM(315)
      ZLM(317)=-ZLM(314)
      ZLM(318)=-TEMP*THR
      ZLM(319)=-ZLM(318)*F15
      ZLM(320)=-ZLM(319)
      ZLM(321)=-ZLM(318)
! L= 9, ML=-5
      TEMP    = SR0247*S11*FAT(13)
      ZLM(322)= TEMP*FIV*F17
      ZLM(323)=-ZLM(322)*TEN
      ZLM(324)= ZLM(322)*FIV
      ZLM(325)=-TEMP*FIV*SIX
      ZLM(326)=-ZLM(325)*TEN
      ZLM(327)= ZLM(325)*FIV
      ZLM(328)= TEMP
      ZLM(329)=-TEMP*TEN
      ZLM(330)= TEMP*FIV
! L= 9, ML=-4
      TEMP    = SR5005*S19*FAT(12)
      ZLM(331)= TEMP*F17
      ZLM(332)= ZLM(331)
      ZLM(333)=-ZLM(331)*SIX
      ZLM(334)=-TEMP*TEN
      ZLM(335)= ZLM(334)
      ZLM(336)=-ZLM(334)*SIX
      ZLM(337)= TEMP
      ZLM(338)= TEMP
      ZLM(339)=-TEMP*SIX
! L= 9, ML=-3
      TEMP    = S3BYS5*S7XS11*S19*FAP(13)
      ZLM(340)= TEMP*2.21D+02
      ZLM(341)=-ZLM(340)*THR
      ZLM(342)=-TEMP*F13*THR*FIV
      ZLM(343)=-ZLM(342)*THR
      ZLM(344)= TEMP*F13*THR
      ZLM(345)=-ZLM(344)*THR
      ZLM(346)=-TEMP
      ZLM(347)= TEMP*THR
! L= 9, ML=-2
      TEMP    = S5XS11*S19*FAT(11)
      ZLM(348)= TEMP*2.21D+02
      ZLM(349)=-ZLM(348)
      ZLM(350)=-TEMP*SEV*F13*THR
      ZLM(351)=-ZLM(350)
      ZLM(352)= TEMP*SEV*F13
      ZLM(353)=-ZLM(352)
      ZLM(354)=-TEMP*SEV
      ZLM(355)=-ZLM(354)
! L= 9, ML=-1
      TEMP    = SR5*S19*FAT(14)
      ZLM(356)= TEMP*2.431D+03
      ZLM(357)=-TEMP*4.004D+03
      ZLM(358)= TEMP*2.002D+03
      ZLM(359)=-TEMP*3.08D+02
      ZLM(360)= TEMP*SEV
! L= 9, ML= 0
      TEMP    = S19*FAP(14)
      ZLM(361)= TEMP*1.2155D+04
      ZLM(362)=-TEMP*2.5740D+04
      ZLM(363)= TEMP*1.8018D+04
      ZLM(364)=-TEMP*4.620D+03
      ZLM(365)= TEMP*3.15D+02
! L= 9, ML= 1
      ZLM(366)= ZLM(356)
      ZLM(367)= ZLM(357)
      ZLM(368)= ZLM(358)
      ZLM(369)= ZLM(359)
      ZLM(370)= ZLM(360)
! L= 9, ML= 2
      TEMP    = S5XS11*S19*FAT(9)
      ZLM(371)= TEMP*2.21D+02
      ZLM(372)=-TEMP*2.73D+02
      ZLM(373)= TEMP*9.1D+01
      ZLM(374)=-TEMP*SEV
! L= 9, ML= 3
      ZLM(375)=-ZLM(341)
      ZLM(376)=-ZLM(340)
      ZLM(377)=-ZLM(343)
      ZLM(378)=-ZLM(342)
      ZLM(379)=-ZLM(345)
      ZLM(380)=-ZLM(344)
      ZLM(381)=-ZLM(347)
      ZLM(382)=-ZLM(346)
! L= 9, ML= 4
      TEMP    = SR5005*S19*FAT(8)
      ZLM(383)= TEMP*F17
      ZLM(384)=-ZLM(383)
      ZLM(385)=-TEMP*TEN
      ZLM(386)=-ZLM(385)
      ZLM(387)= TEMP
      ZLM(388)=-ZLM(387)
! L= 9, ML= 5
      ZLM(389)= ZLM(324)
      ZLM(390)= ZLM(323)
      ZLM(391)= ZLM(322)
      ZLM(392)= ZLM(327)
      ZLM(393)= ZLM(326)
      ZLM(394)= ZLM(325)
      ZLM(395)= ZLM(330)
      ZLM(396)= ZLM(329)
      ZLM(397)= ZLM(328)
! L= 9, ML= 6
      TEMP    = SR0247*S3BYS5*S11*FAP(9)
      ZLM(398)= TEMP*F17*THR
      ZLM(399)=-TEMP*F17*TEN
      ZLM(400)= ZLM(398)
      ZLM(401)=-TEMP*THR*THR
      ZLM(402)= TEMP*THR*TEN
      ZLM(403)= ZLM(401)
! L= 9, ML= 7
      ZLM(404)=-ZLM(309)
      ZLM(405)=-ZLM(308)
      ZLM(406)=-ZLM(307)
      ZLM(407)=-ZLM(306)
      ZLM(408)=-ZLM(313)
      ZLM(409)=-ZLM(312)
      ZLM(410)=-ZLM(311)
      ZLM(411)=-ZLM(310)
! L= 9, ML= 8
      ZLM(412)= ZLM(301)*EIT
      ZLM(413)=-ZLM(412)*SEV
      ZLM(414)=-ZLM(413)
      ZLM(415)=-ZLM(412)
! L= 9, ML= 9
      ZLM(416)= ZLM(300)
      ZLM(417)= ZLM(299)
      ZLM(418)= ZLM(298)
      ZLM(419)= ZLM(297)
      ZLM(420)= ZLM(296)
! ---
! L=10, ML=-10
      ZLM(421)= SR0323*SR3003*FAP(17)
      ZLM(422)=-ZLM(421)*4.5D+01
      ZLM(423)= ZLM(421)*2.10D+02
      ZLM(424)=-ZLM(423)
      ZLM(425)=-ZLM(422)
      ZLM(426)=-ZLM(421)
! L=10, ML=-9
      ZLM(427)= SR0323*SR3003*SR5*FAP(15)
      ZLM(428)=-ZLM(427)*3.6D+01
      ZLM(429)= ZLM(427)*1.26D+02
      ZLM(430)=-ZLM(427)*8.4D+01
      ZLM(431)= ZLM(427)*9.0D+00
! L=10, ML=-8
      TEMP    = SR5005*S17*SR3*FAP(16)
      ZLM(432)= TEMP*F19
      ZLM(433)=-ZLM(432)*SEV*FOR
      ZLM(434)= ZLM(432)*SEV*TEN
      ZLM(435)= ZLM(433)
      ZLM(436)= ZLM(432)
      ZLM(437)=-TEMP
      ZLM(438)=-ZLM(437)*SEV*FOR
      ZLM(439)= ZLM(437)*SEV*TEN
      ZLM(440)= ZLM(438)
      ZLM(441)= ZLM(437)
! L=10, ML=-7
      TEMP    = SR5005*S17*FAT(15)
      ZLM(442)= TEMP*F19
      ZLM(443)=-ZLM(442)*SEV*THR
      ZLM(444)= ZLM(442)*SEV*FIV
      ZLM(445)=-ZLM(442)*SEV
      ZLM(446)=-TEMP*THR
      ZLM(447)=-ZLM(446)*SEV*THR
      ZLM(448)= ZLM(446)*SEV*FIV
      ZLM(449)=-ZLM(446)*SEV
! L=10, ML=-6
      TEMP    = SR5005*FAT(17)
      ZLM(450)= TEMP*3.23D+02
      ZLM(451)=-ZLM(450)*F15
      ZLM(452)=-ZLM(451)
      ZLM(453)=-ZLM(450)
      ZLM(454)=-TEMP*1.02D+02
      ZLM(455)=-ZLM(454)*F15
      ZLM(456)=-ZLM(455)
      ZLM(457)=-ZLM(454)
      ZLM(458)= TEMP*THR
      ZLM(459)=-ZLM(458)*F15
      ZLM(460)=-ZLM(459)
      ZLM(461)=-ZLM(458)
! L=10, ML=-5
      TEMP    = SR1001*FAT(13)
      ZLM(462)= TEMP*3.23D+02
      ZLM(463)=-ZLM(462)*TEN
      ZLM(464)= ZLM(462)*FIV
      ZLM(465)=-TEMP*1.70D+02
      ZLM(466)=-ZLM(465)*TEN
      ZLM(467)= ZLM(465)*FIV
      ZLM(468)= TEMP*F15
      ZLM(469)=-ZLM(468)*TEN
      ZLM(470)= ZLM(468)*FIV
! L=10, ML=-4
      TEMP    = SR5005*FAT(14)
      ZLM(471)= TEMP*3.23D+02
      ZLM(472)= ZLM(471)
      ZLM(473)=-ZLM(471)*SIX
      ZLM(474)=-TEMP*2.55D+02
      ZLM(475)= ZLM(474)
      ZLM(476)=-ZLM(474)*SIX
      ZLM(477)= TEMP*4.5D+01
      ZLM(478)= ZLM(477)
      ZLM(479)=-ZLM(477)*SIX
      ZLM(480)=-TEMP
      ZLM(481)= ZLM(480)
      ZLM(482)=-ZLM(480)*SIX
! L=10, ML=-3
      TEMP    = SR5005*FAT(13)
      ZLM(483)= TEMP*3.23D+02
      ZLM(484)=-ZLM(483)*THR
      ZLM(485)=-TEMP*3.57D+02
      ZLM(486)=-ZLM(485)*THR
      ZLM(487)= TEMP*1.05D+02
      ZLM(488)=-ZLM(487)*THR
      ZLM(489)=-TEMP*SEV
      ZLM(490)=-ZLM(489)*THR
! L=10, ML=-2
      TEMP    = SR0385*FAT(16)
      ZLM(491)= TEMP*4.199D+03
      ZLM(492)=-ZLM(491)
      ZLM(493)=-TEMP*6.188D+03
      ZLM(494)=-ZLM(493)
      ZLM(495)= TEMP*2.730D+03
      ZLM(496)=-ZLM(495)
      ZLM(497)=-TEMP*3.64D+02
      ZLM(498)=-ZLM(497)
      ZLM(499)= TEMP*SEV
      ZLM(500)=-ZLM(499)
! L=10, ML=-1
      TEMP    = SR0385*SR3*FAP(14)
      ZLM(501)= TEMP*4.199D+03
      ZLM(502)=-TEMP*7.956D+03
      ZLM(503)= TEMP*4.914D+03
      ZLM(504)=-TEMP*1.092D+03
      ZLM(505)= TEMP*6.3D+01
! L=10, ML= 0
      TEMP    = S3BYS7*FAP(16)
      ZLM(506)= TEMP*4.6189D+04
      ZLM(507)=-TEMP*1.09395D+05
      ZLM(508)= TEMP*9.0090D+04
      ZLM(509)=-TEMP*3.0030D+04
      ZLM(510)= TEMP*3.465D+03
      ZLM(511)=-TEMP*6.3D+01
! L=10, ML= 1
      ZLM(512)= ZLM(501)
      ZLM(513)= ZLM(502)
      ZLM(514)= ZLM(503)
      ZLM(515)= ZLM(504)
      ZLM(516)= ZLM(505)
! L=10, ML= 2
      ZLM(517)= ZLM(491)*2.0d0
      ZLM(518)= ZLM(493)*2.0d0
      ZLM(519)= ZLM(495)*2.0d0
      ZLM(520)= ZLM(497)*2.0d0
      ZLM(521)= ZLM(499)*2.0d0
! L=10, ML= 3
      ZLM(522)=-ZLM(484)
      ZLM(523)=-ZLM(483)
      ZLM(524)=-ZLM(486)
      ZLM(525)=-ZLM(485)
      ZLM(526)=-ZLM(488)
      ZLM(527)=-ZLM(487)
      ZLM(528)=-ZLM(490)
      ZLM(529)=-ZLM(489)
! L=10, ML= 4
      ZLM(530)= ZLM(471)*FOR
      ZLM(531)=-ZLM(530)
      ZLM(532)= ZLM(474)*FOR
      ZLM(533)=-ZLM(532)
      ZLM(534)= ZLM(477)*FOR
      ZLM(535)=-ZLM(534)
      ZLM(536)= ZLM(480)*FOR
      ZLM(537)=-ZLM(536)
! L=10, ML= 5
      ZLM(538)= ZLM(462)
      ZLM(539)= ZLM(463)
      ZLM(540)= ZLM(464)
      ZLM(541)= ZLM(465)
      ZLM(542)= ZLM(466)
      ZLM(543)= ZLM(467)
      ZLM(544)= ZLM(468)
      ZLM(545)= ZLM(469)
      ZLM(546)= ZLM(470)
! L=10, ML= 6
      ZLM(547)= ZLM(450)*SIX
      ZLM(548)=-ZLM(450)*2.0D+01
      ZLM(549)= ZLM(547)
      ZLM(550)= ZLM(454)*SIX
      ZLM(551)=-ZLM(454)*2.0D+01
      ZLM(552)= ZLM(550)
      ZLM(553)= ZLM(458)*SIX
      ZLM(554)=-ZLM(458)*2.0D+01
      ZLM(555)= ZLM(553)
! L=10, ML= 7
      ZLM(556)=-ZLM(445)
      ZLM(557)=-ZLM(444)
      ZLM(558)=-ZLM(443)
      ZLM(559)=-ZLM(442)
      ZLM(560)=-ZLM(449)
      ZLM(561)=-ZLM(448)
      ZLM(562)=-ZLM(447)
      ZLM(563)=-ZLM(446)
! L=10, ML= 8
      ZLM(564)= ZLM(432)*EIT
      ZLM(565)=-ZLM(432)*5.6D+01
      ZLM(566)=-ZLM(565)
      ZLM(567)=-ZLM(564)
      ZLM(568)= ZLM(437)*EIT
      ZLM(569)=-ZLM(437)*5.6D+01
      ZLM(570)=-ZLM(569)
      ZLM(571)=-ZLM(568)
! L=10, ML= 9
      ZLM(572)= ZLM(431)
      ZLM(573)= ZLM(430)
      ZLM(574)= ZLM(429)
      ZLM(575)= ZLM(428)
      ZLM(576)= ZLM(427)
! L=10, ML= 10
      ZLM(577)= ZLM(421)*TEN
      ZLM(578)=-ZLM(421)*1.20D+02
      ZLM(579)= ZLM(421)*2.52D+02
      ZLM(580)= ZLM(578)
      ZLM(581)= ZLM(577)
!
      RETURN
      END

! ECPINI
      SUBROUTINE ECPINI(LMF,LMX,LMY,LMZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION LMF(122),LMX(581),LMY(581),LMZ(581)
      COMMON /GBASE / NFST(8),NX(84),NY(84),NZ(84)

! LMF (FIRST) POINT TO THE FIRST PARTS OF THE
! SPHERICAL HARMONIC COEF TABLE (ZLM) AND TO THE CORRESPONDING POWERS
! OF X, Y, AND Z STORED IN LMX, LMY, LMZ (GOOD THROUGH L=10)

      DIMENSION ILMF(122),ILMX(581),ILMY(581),ILMZ(581)
      DIMENSION INFST(8),INX(84),INY(84),INZ(84)
      DATA ILMF/1, 2,3,4, 5,7,8,10,11, 12,14,16,18,20,22,23,            &
      25,28,30,34,36,39,41,43,45, 47,50,53,57,61,64,67,70,72,76,78,     &
      81,85,88,94,98,104,107,111,114,117,121,125,128,                   &
      131,135,139,145,151,157,163,167,171,175,178,184,188,194,197,      &
      201,206,210,218,224,233,239,247,251,256,260,264,270,276,282,      &
      288,292, 296,301,306,314,322,331,340,348,356,361,366,371,375,383, &
      389,398,404,412,416, 421,427,432,442,450,462,471,483,491,501,506, &
      512,517,522,530,538,547,556,564,572,577, 582/                      
      DATA ILMX/0, 1,0,0, 2,0,1,0,0,0,1, 3,1,2,0,1,1,0,0,0,0,1,2,0,     &
      4,2,0,3,1,2,0,0,2,1,1,0,0,0,0,0,1,1,2,0,3,1, 5,3,1,0,2,4,3,1,1,3, &
      2,0,2,0,1,1,1,0,0,0,0,0,0,1,1,2,2,0,0,3,1,4,2,0, 6,4,2,0,5,3,1,4, &
      2,0,4,2,0,3,1,1,3,2,0,2,0,2,0,1,1,1,0,0,0,0,0,0,0,1,1,1,2,0,2,0,3,&
      1,3,1,4,2,0,5,3,1, 7,3,5,1,6,4,0,2,5,3,1,5,3,1,4,2,0,4,2,0,3,1,3, &
      1,3,1,2,0,2,0,2,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,2,0,2,0,2,0,3,1,3,&
      1,4,2,0,4,2,0,5,3,1,0,4,2,6, 8,6,4,2,0,7,5,3,1,6,4,2,0,6,4,2,0,5, &
      3,1,5,3,1,4,0,2,4,0,2,4,0,2,3,1,3,1,3,1,2,0,2,0,2,0,2,0,1,1,1,1,0,&
      0,0,0,0,0,0,0,0,1,1,1,1,2,0,2,0,2,0,3,1,3,1,3,1,4,2,0,4,2,0,5,3,1,&
      5,3,1,6,4,2,0,7,5,3,1, 9,7,5,3,1,8,6,4,2,0,7,5,3,1,7,5,3,1,6,4,2, &
      0,6,4,2,0,5,3,1,5,3,1,5,3,1,4,0,2,4,0,2,4,0,2,3,1,3,1,3,1,3,1,2,0,&
      2,0,2,0,2,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,2,0,2,0,2,0,2,0,&
      3,1,3,1,3,1,4,2,0,4,2,0,4,2,0,5,3,1,5,3,1,6,4,2,0,6,4,2,0,7,5,3,1,&
      8,6,4,2,0, 10,8,6,4,2,0,9,7,5,3,1,8,6,4,2,0,8,6,4,2,0,7,5,3,1,7,5,&
      3,1,6,4,2,0,6,4,2,0,6,4,2,0,5,3,1,5,3,1,5,3,1,4,0,2,4,0,2,4,0,2,4,&
      0,2,3,1,3,1,3,1,3,1,2,0,2,0,2,0,2,0,2,0,1,1,1,1,1,0,0,0,0,0,0,0,  &
      0,0,0,0,1,1,1,1,1,2,0,2,0,2,0,2,0,3,1,3,1,3,1,3,1,0,2,4,0,2,4,0,2,&
      4,5,3,1,5,3,1,5,3,1,6,4,2,0,6,4,2,0,7,5,3,1,7,5,3,1,8,6,4,2,0,9,7,&
      5,3,1/                                                             
      DATA ILMY/0, 0,0,1, 0,2,0,0,0,1,1, 0,2,0,2,0,0,0,0,1,1,1,1,3,     &
      0,2,4,0,2,0,2,2,0,0,0,0,0,0,1,1,1,1,1,3,1,3, 0,2,4,4,2,0,0,2,2,0, &
      0,2,0,2,0,0,0,0,0,0,1,1,1,1,1,1,1,3,3,1,3,1,3,5, 0,2,4,6,0,2,4,0, &
      2,4,0,2,4,0,2,2,0,0,2,0,2,0,2,0,0,0,0,0,0,0,1,1,1,1,1,1,1,3,1,3,1,&
      3,1,3,1,3,5,1,3,5, 0,4,2,6,0,2,6,4,0,2,4,0,2,4,0,2,4,0,2,4,0,2,0, &
      2,0,2,0,2,0,2,0,2,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,&
      3,1,3,5,1,3,5,1,3,5,7,3,5,1, 0,2,4,6,8,0,2,4,6,0,2,4,6,0,2,4,6,0, &
      2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,0,0,0,0,&
      0,0,0,0,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,3,1,3,1,3,5,1,3,5,1,3,5,&
      1,3,5,1,3,5,7,1,3,5,7, 0,2,4,6,8,0,2,4,6,8,0,2,4,6,0,2,4,6,0,2,4, &
      6,0,2,4,6,0,2,4,0,2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,2,0,2,0,2,0,2,0,2,&
      0,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,&
      1,3,1,3,1,3,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,7,1,3,5,7,1,3,5,7,&
      1,3,5,7,9, 0,2,4,6,8,10,0,2,4,6,8,0,2,4,6,8,0,2,4,6,8,0,2,4,6,0,2,&
      4,6,0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,0,2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,&
      4,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,1,  &
      1,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,5,3,1,5,3,1,5,3,&
      1,1,3,5,1,3,5,1,3,5,1,3,5,7,1,3,5,7,1,3,5,7,1,3,5,7,1,3,5,7,9,1,3,&
      5,7,9/                                                             
      DATA ILMZ/0, 0,1,0, 0,0,1,2,0,1,0, 0,0,1,1,2,0,3,1,2,0,1,0,0,     &
      0,0,0,1,1,2,2,0,0,3,1,4,2,0,3,1,2,0,1,1,0,0, 0,0,0,1,1,1,2,2,0,0, &
      3,3,1,1,4,2,0,5,3,1,4,2,0,3,1,2,0,2,0,1,1,0,0,0, 0,0,0,0,1,1,1,2, &
      2,2,0,0,0,3,3,1,1,4,4,2,2,0,0,5,3,1,6,4,2,0,5,3,1,4,2,0,3,3,1,1,2,&
      2,0,0,1,1,1,0,0,0, 0,0,0,0,1,1,1,1,2,2,2,0,0,0,3,3,3,1,1,1,4,4,2, &
      2,0,0,5,5,3,3,1,1,6,4,2,0,7,5,3,1,6,4,2,0,5,3,1,4,4,2,2,0,0,3,3,1,&
      1,2,2,2,0,0,0,1,1,1,0,0,0,0, 0,0,0,0,0,1,1,1,1,2,2,2,2,0,0,0,0,3, &
      3,3,1,1,1,4,4,4,2,2,2,0,0,0,5,5,3,3,1,1,6,6,4,4,2,2,0,0,7,5,3,1,8,&
      6,4,2,0,7,5,3,1,6,4,2,0,5,5,3,3,1,1,4,4,2,2,0,0,3,3,3,1,1,1,2,2,2,&
      0,0,0,1,1,1,1,0,0,0,0, 0,0,0,0,0,1,1,1,1,1,2,2,2,2,0,0,0,0,3,3,3, &
      3,1,1,1,1,4,4,4,2,2,2,0,0,0,5,5,5,3,3,3,1,1,1,6,6,4,4,2,2,0,0,7,7,&
      5,5,3,3,1,1,8,6,4,2,0,9,7,5,3,1,8,6,4,2,0,7,5,3,1,6,6,4,4,2,2,0,0,&
      5,5,3,3,1,1,4,4,4,2,2,2,0,0,0,3,3,3,1,1,1,2,2,2,2,0,0,0,0,1,1,1,1,&
      0,0,0,0,0, 0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,0,0,0,0,0,3,3,3,3,1,1, &
      1,1,4,4,4,4,2,2,2,2,0,0,0,0,5,5,5,3,3,3,1,1,1,6,6,6,4,4,4,2,2,2,0,&
      0,0,7,7,5,5,3,3,1,1,8,8,6,6,4,4,2,2,0,0,9,7,5,3,1,10,8,6,4,2,0,9, &
      7,5,3,1,8,6,4,2,0,7,7,5,5,3,3,1,1,6,6,4,4,2,2,0,0,5,5,5,3,3,3,1,1,&
      1,4,4,4,2,2,2,0,0,0,3,3,3,3,1,1,1,1,2,2,2,2,0,0,0,0,1,1,1,1,1,0,0,&
      0,0,0/                                                             
      DATA INFST/1,2,5,11,21,36,57,85/                                   
      DATA INX/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,   &
      1,0,2,2,0,2,1,1, 5,0,0,4,4,1,0,1,0,3,3,2,0,2,0,3,1,1,2,2,1,       &
      6,0,0,5,5,1,0,1,0,4,4,2,0,2,0,4,1,1,3,3,0,3,3,2,1,2,1,2/,         &
           INY/0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,   &
      0,1,2,0,2,1,2,1, 0,5,0,1,0,4,4,0,1,2,0,3,3,0,2,1,3,1,2,1,2,       &
      0,6,0,1,0,5,5,0,1,2,0,4,4,0,2,1,4,1,3,0,3,2,1,3,3,1,2,2/,         &
           INZ/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1,0,0,4,0,1,0,1,   &
      3,3,0,2,2,1,1,2, 0,0,5,0,1,0,1,4,4,0,2,0,2,3,3,1,1,3,1,2,2,       &
      0,0,6,0,1,0,1,5,5,0,2,0,2,4,4,1,1,4,0,3,3,1,2,1,2,3,3,2/

      DO 110 I = 1,122
         LMF(I)= ILMF(I)
  110 CONTINUE
      DO 120 J = 1,581
         LMX(J)= ILMX(J)
         LMY(J)= ILMY(J)
         LMZ(J)= ILMZ(J)
  120 CONTINUE
      DO 130 K = 1,8
         NFST(K)= INFST(K)
  130 CONTINUE
      DO 140 L = 1,84
         NX(L) = INX(L)
         NY(L) = INY(L)
         NZ(L) = INZ(L)
  140 CONTINUE
      RETURN
      END

! DAWERT
      SUBROUTINE DAWERT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!        -----  ROUTINE ALLOCATES Parameters SPECIFYING THE     -----
!        -----  PIECEWISE CHEBYCHEV POLYNOMIAL FIT TO THE       -----
!        -----  DAWSON-ERROR FUNCTION.                          -----
!
      COMMON /DERFCM/ C(246),IFRST(40),ILAST(40),H
!
      H = 0.25D+00
! ... 0.00 < X <= 0.25  INTERVAL NO. 1, ABS.ERROR = 3.0831115438446D-
      IFRST( 1)= 1
      ILAST( 1)= 8
      C( 1)=-0.2227673357873D-15
      C( 2)=-0.3603041588494D-08
      C( 3)= 0.5641900599881D+00
      C( 4)=-0.1855716761667D-04
      C( 5)=-0.3758028089842D+00
      C( 6)=-0.2937263239086D-02
      C( 7)= 0.1586735773807D+00
      C( 8)=-0.3712362707925D-01
! ... 0.25 < X <= 0.50  INTERVAL NO. 2, ABS.ERROR = 1.6200374375330D-
      IFRST( 2)= 9
      ILAST( 2)=16
      C( 9)= 0.1654923320373D-05
      C(10)=-0.4283673776908D-04
      C(11)= 0.5646780085194D+00
      C(12)=-0.3216632920741D-02
      C(13)=-0.3626661648377D+00
      C(14)=-0.3697874666458D-01
      C(15)= 0.2103001701512D+00
      C(16)=-0.7233313577516D-01
! ... 0.50 < X <= 0.75  INTERVAL NO. 3, ABS.ERROR = 2.6527224861184D-
      IFRST( 3)=17
      ILAST( 3)=24
      C(17)=-0.4909038422764D-03
      C(18)= 0.6169678713757D-02
      C(19)= 0.5309521767077D+00
      C(20)= 0.9902446632012D-01
      C(21)=-0.5497719612405D+00
      C(22)= 0.1699289458172D+00
      C(23)= 0.8215358267937D-01
      C(24)=-0.3800879631724D-01
! ... 0.75 < X <= 1.00  INTERVAL NO. 4, ABS.ERROR = 2.1188384380366D-
      IFRST( 4)=25
      ILAST( 4)=32
      C(25)=-0.8860184581188D-02
      C(26)= 0.8182679610811D-01
      C(27)= 0.2361899484096D+00
      C(28)= 0.7408443731368D+00
      C(29)=-0.1393552470603D+01
      C(30)= 0.8398504617092D+00
      C(31)=-0.2153157176716D+00
      C(32)= 0.1898293835776D-01
! ... 1.00 < X <= 1.25  INTERVAL NO. 5, ABS.ERROR = 1.0054179711005D-
      IFRST( 5)=33
      ILAST( 5)=40
      C(33)=-0.2622864777452D-01
      C(34)= 0.2070481292278D+00
      C(35)=-0.1517729522946D+00
      C(36)= 0.1410375745664D+01
      C(37)=-0.2088589931697D+01
      C(38)= 0.1273821623264D+01
      C(39)=-0.3662054027830D+00
      C(40)= 0.4151758125850D-01
! ... 1.25 < X <= 1.50  INTERVAL NO. 6, ABS.ERROR = 1.2743583965857D-
      IFRST( 6)=41
      ILAST( 6)=48
      C(41)= 0.7652753217902D-01
      C(42)=-0.3496979906426D+00
      C(43)= 0.1142731303255D+01
      C(44)=-0.2641065409391D+00
      C(45)=-0.7870989692331D+00
      C(46)= 0.6659477061060D+00
      C(47)=-0.2082257760423D+00
      C(48)= 0.2389272621700D-01
! ... 1.50 < X <= 1.75  INTERVAL NO. 7, ABS.ERROR = 3.2651215065016D-
      IFRST( 7)=49
      ILAST( 7)=55
      C(49)= 0.4634935721533D+00
      C(50)=-0.2162775089324D+01
      C(51)= 0.4787745789535D+01
      C(52)=-0.4340108008393D+01
      C(53)= 0.1951017416635D+01
      C(54)=-0.4390388356211D+00
      C(55)= 0.3981846570969D-01
! ... 1.75 < X <= 2.00  INTERVAL NO. 8, ABS.ERROR = 8.9706020389713D-
      IFRST( 8)=56
      ILAST( 8)=63
      C(56)= 0.9892135872598D+00
      C(57)=-0.4311631180681D+01
      C(58)= 0.8556855041877D+01
      C(59)=-0.8017503990110D+01
      C(60)= 0.4106390132111D+01
      C(61)=-0.1197913275632D+01
      C(62)= 0.1884269501482D+00
      C(63)=-0.1248598098755D-01
! ... 2.00 < X <= 2.25  INTERVAL NO. 9, ABS.ERROR = 3.1388225352202D-
      IFRST( 9)=64
      ILAST( 9)=71
      C(64)= 0.4508582229747D+00
      C(65)=-0.2475021032343D+01
      C(66)= 0.5870952943693D+01
      C(67)=-0.5834804573043D+01
      C(68)= 0.3041872702481D+01
      C(69)=-0.8863349469112D+00
      C(70)= 0.1377496059452D+00
      C(71)=-0.8952617645264D-02
! ... 2.25 < X <= 2.50  INTERVAL NO.10, ABS.ERROR = 3.5207392556913D-
      IFRST(10)=72
      ILAST(10)=78
      C(72)=-0.2479150925527D+01
      C(73)= 0.6465442349154D+01
      C(74)=-0.5832144628726D+01
      C(75)= 0.2684557558935D+01
      C(76)=-0.6830155078836D+00
      C(77)= 0.9186854104822D-01
      C(78)=-0.5120121563474D-02
! ... 2.50 < X <= 2.75  INTERVAL NO.11, ABS.ERROR = 1.5842438472191D-
      IFRST(11)=79
      ILAST(11)=85
      C(79)=-0.2550246682056D+01
      C(80)= 0.6663085650406D+01
      C(81)=-0.6057008239737D+01
      C(82)= 0.2819123791215D+01
      C(83)=-0.7278168845029D+00
      C(84)= 0.9975271609922D-01
      C(85)=-0.5693960934877D-02
! ... 2.75 < X <= 3.00  INTERVAL NO.12, ABS.ERROR = 2.1390000881638D-
      IFRST(12)=86
      ILAST(12)=92
      C(86)=-0.1437322524459D+01
      C(87)= 0.4244945032011D+01
      C(88)=-0.3866819653025D+01
      C(89)= 0.1760657886042D+01
      C(90)=-0.4399493394691D+00
      C(91)= 0.5797869518089D-01
      C(92)=-0.3166941925883D-02
! ... 3.00 < X <= 3.25  INTERVAL NO.13, ABS.ERROR = 1.2176926134089D-
      IFRST(13)=93
      ILAST(13)=99
      C(93)= 0.4326982515539D-01
      C(94)= 0.1278541925538D+01
      C(95)=-0.1389457146276D+01
      C(96)= 0.6567745067205D+00
      C(97)=-0.1631568799154D+00
      C(98)= 0.2094831465123D-01
      C(99)=-0.1101922864715D-02
! ... 3.25 < X <= 3.50  INTERVAL NO.14, ABS.ERROR = 4.2543746303636D-
      IFRST(14)=100
      ILAST(14)=106
      C(100)= 0.1115286202567D+01
      C(101)=-0.7098356509943D+00
      C(102)= 0.1477452918517D+00
      C(103)= 0.2274875774917D-01
      C(104)=-0.1601072077271D-01
      C(105)= 0.2728999281923D-02
      C(106)=-0.1616695274909D-03
! ... 3.50 < X <= 3.75  INTERVAL NO.15, ABS.ERROR = 2.3193003073629D-
      IFRST(15)=107
      ILAST(15)=112
      C(107)= 0.1346851272697D+01
      C(108)=-0.1125834241251D+01
      C(109)= 0.4583884751028D+00
      C(110)=-0.1007033665315D+00
      C(111)= 0.1153266357142D-01
      C(112)=-0.5427038297057D-03
! ... 3.75 < X <= 4.00  INTERVAL NO.16, ABS.ERROR = 2.6108004647085D-
      IFRST(16)=113
      ILAST(16)=118
      C(113)= 0.1220883963156D+01
      C(114)=-0.9581085557663D+00
      C(115)= 0.3690371210149D+00
      C(116)=-0.7689811239270D-01
      C(117)= 0.8360797545174D-02
      C(118)=-0.3736140672117D-03
! ... 4.00 < X <= 4.25  INTERVAL NO.17, ABS.ERROR = 1.7570833676928D-
      IFRST(17)=119
      ILAST(17)=124
      C(119)= 0.1079863694049D+01
      C(120)=-0.7816343308043D+00
      C(121)= 0.2806819253411D+00
      C(122)=-0.5477512389825D-01
      C(123)= 0.5590565502644D-02
      C(124)=-0.2348302397877D-03
! ... 4.25 < X <= 4.50  INTERVAL NO.18, ABS.ERROR = 1.0521805648978D-
      IFRST(18)=125
      ILAST(18)=130
      C(125)= 0.9600857053756D+00
      C(126)=-0.6404875010594D+00
      C(127)= 0.2141383549717D+00
      C(128)=-0.3908625281056D-01
      C(129)= 0.3740753346938D-02
      C(130)=-0.1475725788623D-03
! ... 4.50 < X <= 4.75  INTERVAL NO.19, ABS.ERROR = 6.2718719107124D-
      IFRST(19)=131
      ILAST(19)=136
      C(131)= 0.8652746664820D+00
      C(132)=-0.5349723088368D+00
      C(133)= 0.1671594028443D+00
      C(134)=-0.2862620440719D-01
      C(135)= 0.2576075398247D-02
      C(136)=-0.9569143876433D-04
! ... 4.75 < X <= 5.00  INTERVAL NO.20, ABS.ERROR = 3.9639402871217D-
      IFRST(20)=137
      ILAST(20)=142
      C(137)= 0.7900679556534D+00
      C(138)=-0.4556957919175D+00
      C(139)= 0.1337277737706D+00
      C(140)=-0.2157594052042D-01
      C(141)= 0.1832563601783D-02
      C(142)=-0.6432286463678D-04
! ... 5.00 < X <= 5.25  INTERVAL NO.21, ABS.ERROR = 2.4784618801732D-
      IFRST(21)=143
      ILAST(21)=148
      C(143)= 0.7290302365904D+00
      C(144)=-0.3945841937680D+00
      C(145)= 0.1092502418252D+00
      C(146)=-0.1667318620030D-01
      C(147)= 0.1341496314853D-02
      C(148)=-0.4464581143111D-04
! ... 5.25 < X <= 5.50  INTERVAL NO.22, ABS.ERROR = 1.6431300764452D-
      IFRST(22)=149
      ILAST(22)=154
      C(149)= 0.6781821897825D+00
      C(150)=-0.3461064436915D+00
      C(151)= 0.9076079662805D-01
      C(152)=-0.1314681150288D-01
      C(153)= 0.1005173704471D-02
      C(154)=-0.3181374631822D-04
! ... 5.50 < X <= 5.75  INTERVAL NO.23, ABS.ERROR = 1.1222134332911D-
      IFRST(23)=155
      ILAST(23)=160
      C(155)= 0.6349563489690D+00
      C(156)=-0.3067739392218D+00
      C(157)= 0.7644326048783D-01
      C(158)=-0.1054063730153D-01
      C(159)= 0.7679506117711D-03
      C(160)=-0.2317563630640D-04
! ... 5.75 < X <= 6.00  INTERVAL NO.24, ABS.ERROR = 7.8559381222476D-
      IFRST(24)=161
      ILAST(24)=166
      C(161)= 0.5976075476333D+00
      C(162)=-0.2742701592908D+00
      C(163)= 0.6512719379650D-01
      C(164)=-0.8570613370193D-02
      C(165)= 0.5964515294181D-03
      C(166)=-0.1720313448459D-04
! ... 6.00 < X <= 6.25  INTERVAL NO.25, ABS.ERROR = 4.6252335295094D-
      IFRST(25)=167
      ILAST(25)=171
      C(167)= 0.4530821605957D+00
      C(168)=-0.1556937202255D+00
      C(169)= 0.2621052097467D-01
      C(170)=-0.2184050212463D-02
      C(171)= 0.7237572572194D-04
! ... 6.25 < X <= 6.50  INTERVAL NO.26, ABS.ERROR = 3.5389913222161D-
      IFRST(26)=172
      ILAST(26)=176
      C(172)= 0.4314398767018D+00
      C(173)=-0.1418342842259D+00
      C(174)= 0.2288200151943D-01
      C(175)=-0.1828741476402D-02
      C(176)= 0.5815166332468D-04
! ... 6.50 < X <= 6.75  INTERVAL NO.27, ABS.ERROR = 2.7426949600340D-
      IFRST(27)=177
      ILAST(27)=181
      C(177)= 0.4119346702408D+00
      C(178)=-0.1298244601295D+00
      C(179)= 0.2010878728234D-01
      C(180)=-0.1544113849832D-02
      C(181)= 0.4719618664240D-04
! ... 6.75 < X <= 7.00  INTERVAL NO.28, ABS.ERROR = 2.1502799540940D-
      IFRST(28)=182
      ILAST(28)=186
      C(182)= 0.3942461802274D+00
      C(183)=-0.1193370928285D+00
      C(184)= 0.1777693657917D-01
      C(185)=-0.1313662248549D-02
      C(186)= 0.3865502003464D-04
! ... 7.00 < X <= 7.25  INTERVAL NO.29, ABS.ERROR = 1.7034818000639D-
      IFRST(29)=187
      ILAST(29)=191
      C(187)= 0.3781177456714D+00
      C(188)=-0.1101165554215D+00
      C(189)= 0.1580007215659D-01
      C(190)=-0.1125279834923D-02
      C(191)= 0.3192276744812D-04
! ... 7.25 < X <= 7.50  INTERVAL NO.30, ABS.ERROR = 1.3624212868990D-
      IFRST(30)=192
      ILAST(30)=196
      C(192)= 0.3633408065001D+00
      C(193)=-0.1019602554247D+00
      C(194)= 0.1411174419289D-01
      C(195)=-0.9699476170226D-03
      C(196)= 0.2656330616446D-04
! ... 7.50 < X <= 7.75  INTERVAL NO.31, ABS.ERROR = 1.0990319765369D-
      IFRST(31)=197
      ILAST(31)=201
      C(197)= 0.3497439218126D+00
      C(198)=-0.9470569276499D-01
      C(199)= 0.1266017939918D-01
      C(200)=-0.8408550675085D-03
      C(201)= 0.2225784919574D-04
! ... 7.75 < X <= 8.00  INTERVAL NO.32, ABS.ERROR = 8.9381835266522D-
      IFRST(32)=202
      ILAST(32)=206
      C(202)= 0.3371845853069D+00
      C(203)=-0.8822105128463D-01
      C(204)= 0.1140456233287D-01
      C(205)=-0.7327946570967D-03
      C(206)= 0.1877023896668D-04
! ... 8.00 < X <= 8.25  INTERVAL NO.33, ABS.ERROR = 7.3210326689832D-
      IFRST(33)=207
      ILAST(33)=211
      C(207)= 0.3255431584494D+00
      C(208)=-0.8239832407718D-01
      C(209)= 0.1031237441518D-01
      C(210)=-0.6417393667562D-03
      C(211)= 0.1592339503986D-04
! ... 8.25 < X <= 8.50  INTERVAL NO.34, ABS.ERROR = 6.0371707633067D-
      IFRST(34)=212
      ILAST(34)=216
      C(212)= 0.3147180811939D+00
      C(213)=-0.7714810368023D-01
      C(214)= 0.9357439877093D-02
      C(215)=-0.5645414814808D-03
      C(216)= 0.1358301324217D-04
! ... 8.50 < X <= 8.75  INTERVAL NO.35, ABS.ERROR = 5.0113246885530D-
      IFRST(35)=217
      ILAST(35)=221
      C(217)= 0.3046231027165D+00
      C(218)=-0.7239608540250D-01
      C(219)= 0.8518560948886D-02
      C(220)=-0.4987218951555D-03
      C(221)= 0.1164632612927D-04
! ... 8.75 < X <= 9.00  INTERVAL NO.36, ABS.ERROR = 4.1842085352073D-
      IFRST(36)=222
      ILAST(36)=226
      C(222)= 0.2951839772317D+00
      C(223)=-0.6807982796391D-01
      C(224)= 0.7778392737741D-02
      C(225)=-0.4423078514719D-03
      C(226)= 0.1003385659715D-04
! ... 9.00 < X <= 9.25  INTERVAL NO.37, ABS.ERROR = 3.5125236053091D-
      IFRST(37)=227
      ILAST(37)=231
      C(227)= 0.2863364933588D+00
      C(228)=-0.6414655846973D-01
      C(229)= 0.7122648855400D-02
      C(230)=-0.3937177713169D-03
      C(231)= 0.8683627129358D-05
! ... 9.25 < X <= 9.50  INTERVAL NO.38, ABS.ERROR = 2.9645175203541D-
      IFRST(38)=232
      ILAST(38)=236
      C(232)= 0.2780246155388D+00
      C(233)=-0.6055132372497D-01
      C(234)= 0.6539470823348D-02
      C(235)=-0.3516734769846D-03
      C(236)= 0.7546893357357D-05
! ... 9.50 < X <= 9.75  INTERVAL NO.39, ABS.ERROR = 2.5142110615661D-
      IFRST(39)=237
      ILAST(39)=241
      C(237)= 0.2701996594546D+00
      C(238)=-0.5725581832360D-01
      C(239)= 0.6018987049956D-02
      C(240)=-0.3151372568482D-03
      C(241)= 0.6585092705791D-05
! ... 9.75 < X <= 10.0  INTERVAL NO.40, ABS.ERROR = 2.1422863483167D-
      IFRST(40)=242
      ILAST(40)=246
      C(242)= 0.2628187673463D+00
      C(243)=-0.5422707501601D-01
      C(244)= 0.5552907092188D-02
      C(245)=-0.2832594125266D-03
      C(246)= 0.5767453330918D-05
      RETURN
      END

! DAWT
      SUBROUTINE DAWT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!        -----  ROUTINE ALLOCATES Parameters SPECIFYING         -----
!        -----  THE PIECEWISE CHEBYSHEV POLONOMIAL FIT TO       -----
!        -----  THE DAWSON FUNCTION.                            -----
!
      COMMON /DAWFCM/ C(249),IFRST(40),ILAST(40),H
!
      H = 0.25D+00
! ... 0.00 < X <= 0.25  INTERVAL NO. 1, ABS.ERROR = 1.6571632954765D-
      IFRST( 1)= 1
      ILAST( 1)= 8
      C( 1)=-0.1587492209085D-15
      C( 2)= 0.1000000001679D+01
      C( 3)=-0.2200140634227D-06
      C( 4)=-0.6666582802511D+00
      C( 5)=-0.1406923341036D-03
      C( 6)= 0.2678610335986D+00
      C( 7)=-0.5192471402032D-02
      C( 8)=-0.6640836596489D-01
! ... 0.25 < X <= 0.50  INTERVAL NO. 2, ABS.ERROR = 3.1585400961376D-
      IFRST( 2)= 9
      ILAST( 2)=16
      C( 9)=-0.1012320887742D-04
      C(10)= 0.1000229551438D+01
      C(11)=-0.2249318275452D-02
      C(12)=-0.6542497725108D+00
      C(13)=-0.4207802992979D-01
      C(14)= 0.3555122257011D+00
      C(15)=-0.1112488455006D+00
      C(16)=-0.8490860462189D-02
! ... 0.50 < X <= 0.75  INTERVAL NO. 3, ABS.ERROR = 1.5242918038894D-
      IFRST( 3)=17
      ILAST( 3)=24
      C(17)=-0.4921356497286D-03
      C(18)= 0.1006613917835D+01
      C(19)=-0.3889997394948D-01
      C(20)=-0.5358795704274D+00
      C(21)=-0.2746372044452D+00
      C(22)= 0.6337005279825D+00
      C(23)=-0.2989429214171D+00
      C(24)= 0.4661336966923D-01
! ... 0.75 < X <= 1.00  INTERVAL NO. 4, ABS.ERROR = 9.9014130228170D-
      IFRST( 4)=25
      ILAST( 4)=32
      C(25)= 0.2895369459824D-03
      C(26)= 0.1001188430399D+01
      C(27)=-0.2450339520101D-01
      C(28)=-0.5518735978257D+00
      C(29)=-0.2745338921003D+00
      C(30)= 0.6506738753191D+00
      C(31)=-0.3141764870712D+00
      C(32)= 0.5101503644671D-01
! ... 1.00 < X <= 1.25  INTERVAL NO. 5, ABS.ERROR = 1.8118839761883D-
      IFRST( 5)=33
      ILAST( 5)=40
      C(33)= 0.4002056337725D-01
      C(34)= 0.7312836053563D+00
      C(35)= 0.7635870948288D+00
      C(36)=-0.1834201228325D+01
      C(37)= 0.9813999234300D+00
      C(38)=-0.8984573983721D-01
      C(39)=-0.7076886296272D-01
      C(40)= 0.1660415104457D-01
! ... 1.25 < X <= 1.50  INTERVAL NO. 6, ABS.ERROR = 8.9031004790741D-
      IFRST( 6)=41
      ILAST( 6)=48
      C(41)= 0.1852097215132D+00
      C(42)=-0.7878409599634D-01
      C(43)= 0.2705240106371D+01
      C(44)=-0.4425988758096D+01
      C(45)= 0.3062214704197D+01
      C(46)=-0.1094631867084D+01
      C(47)= 0.1994346954993D+00
      C(48)=-0.1461141450065D-01
! ... 1.50 < X <= 1.75  INTERVAL NO. 7, ABS.ERROR = 2.2932766796657D-
      IFRST( 7)=49
      ILAST( 7)=56
      C(49)= 0.2566795947043D+00
      C(50)=-0.4340789653964D+00
      C(51)= 0.3460647185102D+01
      C(52)=-0.5316677310217D+01
      C(53)= 0.3691359889577D+01
      C(54)=-0.1360915547303D+01
      C(55)= 0.2619748477425D+00
      C(56)=-0.2090002809252D-01
! ... 1.75 < X <= 2.00  INTERVAL NO. 8, ABS.ERROR = 5.7358562344234D-
      IFRST( 8)=57
      ILAST( 8)=64
      C(57)=-0.2953293234440D+00
      C(58)= 0.1741536362797D+01
      C(59)=-0.2179107446162D+00
      C(60)=-0.1857700863729D+01
      C(61)= 0.1737803495672D+01
      C(62)=-0.6982079451638D+00
      C(63)= 0.1369430933680D+00
      C(64)=-0.1077900614057D-01
! ... 2.00 < X <= 2.25  INTERVAL NO. 9, ABS.ERROR = 7.9634077110313D-
      IFRST( 9)=65
      ILAST( 9)=71
      C(65)=-0.1683235354314D+01
      C(66)= 0.6579307910802D+01
      C(67)=-0.7452680221211D+01
      C(68)= 0.4159630869064D+01
      C(69)=-0.1268324481886D+01
      C(70)= 0.2038507675752D+00
      C(71)=-0.1359998434782D-01
! ... 2.25 < X <= 2.50  INTERVAL NO.10, ABS.ERROR = 4.9522164147220D-
      IFRST(10)=72
      ILAST(10)=78
      C(72)=-0.1172884307439D+01
      C(73)= 0.5243170888632D+01
      C(74)=-0.5994671204809D+01
      C(75)= 0.3310821796507D+01
      C(76)=-0.9902715162510D+00
      C(77)= 0.1552556976676D+00
      C(78)=-0.1006003345052D-01
! ... 2.50 < X <= 2.75  INTERVAL NO.11, ABS.ERROR = 4.0683012514364D-
      IFRST(11)=79
      ILAST(11)=85
      C(79)= 0.2453071419396D+00
      C(80)= 0.1843675794804D+01
      C(81)=-0.2597381166786D+01
      C(82)= 0.1499064046402D+01
      C(83)=-0.4464691670097D+00
      C(84)= 0.6815297792976D-01
      C(85)=-0.4243532816569D-02
! ... 2.75 < X <= 3.00  INTERVAL NO.12, ABS.ERROR = 1.7797319173951D-
      IFRST(12)=86
      ILAST(12)=92
      C(86)= 0.1711553394682D+01
      C(87)=-0.1366851519128D+01
      C(88)= 0.3331226440003D+00
      C(89)= 0.7176989090082D-01
      C(90)=-0.5525697253082D-01
      C(91)= 0.1093749119900D-01
      C(92)=-0.7553007453680D-03
! ... 3.00 < X <= 3.25  INTERVAL NO.13, ABS.ERROR = 3.6681768733615D-
      IFRST(13)=93
      ILAST(13)=99
      C(93)= 0.2517095374277D+01
      C(94)=-0.2990903746754D+01
      C(95)= 0.1697891140418D+01
      C(96)=-0.5401196813139D+00
      C(97)= 0.9911334144514D-01
      C(98)=-0.9840574037905D-02
      C(99)= 0.4103935013215D-03
! ... 3.25 < X <= 3.50  INTERVAL NO.14, ABS.ERROR = 1.1128875598843D-
      IFRST(14)=100
      ILAST(14)=106
      C(100)= 0.2590851813134D+01
      C(101)=-0.3134054423500D+01
      C(102)= 0.1813437074513D+01
      C(103)=-0.5897752524597D+00
      C(104)= 0.1110979078173D+00
      C(105)=-0.1138103745567D-01
      C(106)= 0.4927876094977D-03
! ... 3.50 < X <= 3.75  INTERVAL NO.15, ABS.ERROR = 1.4752643551219D-
      IFRST(15)=107
      ILAST(15)=113
      C(107)= 0.2263419843285D+01
      C(108)=-0.2574296380066D+01
      C(109)= 0.1414605080120D+01
      C(110)=-0.4381747671623D+00
      C(111)= 0.7867466333846D-01
      C(112)=-0.7681612002974D-02
      C(113)= 0.3168638795614D-03
! ... 3.75 < X <= 4.00  INTERVAL NO.16, ABS.ERROR = 9.2281737806843D-
      IFRST(16)=114
      ILAST(16)=120
      C(114)= 0.1858711024930D+01
      C(115)=-0.1925962656395D+01
      C(116)= 0.9817358192129D+00
      C(117)=-0.2839951967696D+00
      C(118)= 0.4777651946157D-01
      C(119)=-0.4378302333256D-02
      C(120)= 0.1696776598692D-03
! ... 4.00 < X <= 4.25  INTERVAL NO.17, ABS.ERROR = 2.0758506025231D-
      IFRST(17)=121
      ILAST(17)=126
      C(121)= 0.1097210930469D+01
      C(122)=-0.8018697789038D+00
      C(123)= 0.2901340972873D+00
      C(124)=-0.5698500686831D-01
      C(125)= 0.5849148664856D-02
      C(126)=-0.2469443250448D-03
! ... 4.25 < X <= 4.50  INTERVAL NO.18, ABS.ERROR = 1.1124878795954D-
      IFRST(18)=127
      ILAST(18)=132
      C(127)= 0.9640388573940D+00
      C(128)=-0.6448585002328D+00
      C(129)= 0.2160731515099D+00
      C(130)=-0.3951479965599D-01
      C(131)= 0.3788248790079D-02
      C(132)=-0.1496796030551D-03
! ... 4.50 < X <= 4.75  INTERVAL NO.19, ABS.ERROR = 6.3691274476696D-
      IFRST(19)=133
      ILAST(19)=138
      C(133)= 0.8660393462895D+00
      C(134)=-0.5357754986697D+00
      C(135)= 0.1674970705868D+00
      C(136)=-0.2869722599326D-01
      C(137)= 0.2583548621624D-02
      C(138)=-0.9600615594536D-04
! ... 4.75 < X <= 5.00  INTERVAL NO.20, ABS.ERROR = 3.8777869804107D-
      IFRST(20)=139
      ILAST(20)=144
      C(139)= 0.7902176622367D+00
      C(140)=-0.4558460385033D+00
      C(141)= 0.1337881166153D+00
      C(142)=-0.2158806335901D-01
      C(143)= 0.1833781835739D-02
      C(144)=-0.6437185220420D-04
! ... 5.00 < X <= 5.25  INTERVAL NO.21, ABS.ERROR = 2.4793500585929D-
      IFRST(21)=145
      ILAST(21)=150
      C(145)= 0.7290485015677D+00
      C(146)=-0.3946016171968D+00
      C(147)= 0.1092568925649D+00
      C(148)=-0.1667445598850D-01
      C(149)= 0.1341617573053D-02
      C(150)=-0.4465044476092D-04
! ... 5.25 < X <= 5.50  INTERVAL NO.22, ABS.ERROR = 1.6431300764452D-
      IFRST(22)=151
      ILAST(22)=156
      C(151)= 0.6781842728186D+00
      C(152)=-0.3461083410440D+00
      C(153)= 0.9076148811696D-01
      C(154)=-0.1314693754430D-01
      C(155)= 0.1005185194663D-02
      C(156)=-0.3181416541338D-04
! ... 5.50 < X <= 5.75  INTERVAL NO.23, ABS.ERROR = 1.1226575225010D-
      IFRST(23)=157
      ILAST(23)=162
      C(157)= 0.6349565075405D+00
      C(158)=-0.3067740750644D+00
      C(159)= 0.7644330701296D-01
      C(160)=-0.1054064526443D-01
      C(161)= 0.7679512928007D-03
      C(162)=-0.2317565958947D-04
! ... 5.75 < X <= 6.00  INTERVAL NO.24, ABS.ERROR = 7.8470563380506D-
      IFRST(24)=163
      ILAST(24)=168
      C(163)= 0.5976070486768D+00
      C(164)=-0.2742697363056D+00
      C(165)= 0.6512705037617D-01
      C(166)=-0.8570589057672D-02
      C(167)= 0.5964494688669D-03
      C(168)=-0.1720306463540D-04
! ... 6.00 < X <= 6.25  INTERVAL NO.25, ABS.ERROR = 5.6132876125048D-
      IFRST(25)=169
      ILAST(25)=174
      C(169)= 0.5649180454987D+00
      C(170)=-0.2470090261752D+00
      C(171)= 0.5603266263223D-01
      C(172)=-0.7053466938578D-02
      C(173)= 0.4698947566794D-03
      C(174)=-0.1297991257161D-04
! ... 6.25 < X <= 6.50  INTERVAL NO.26, ABS.ERROR = 3.5389913222161D-
      IFRST(26)=175
      ILAST(26)=179
      C(175)= 0.4314398767018D+00
      C(176)=-0.1418342842259D+00
      C(177)= 0.2288200151943D-01
      C(178)=-0.1828741476402D-02
      C(179)= 0.5815166332468D-04
! ... 6.50 < X <= 6.75  INTERVAL NO.27, ABS.ERROR = 2.7426949600340D-
      IFRST(27)=180
      ILAST(27)=184
      C(180)= 0.4119346702408D+00
      C(181)=-0.1298244601295D+00
      C(182)= 0.2010878728234D-01
      C(183)=-0.1544113849832D-02
      C(184)= 0.4719618664240D-04
! ... 6.75 < X <= 7.00  INTERVAL NO.28, ABS.ERROR = 2.1502799540940D-
      IFRST(28)=185
      ILAST(28)=189
      C(185)= 0.3942461802274D+00
      C(186)=-0.1193370928285D+00
      C(187)= 0.1777693657917D-01
      C(188)=-0.1313662248549D-02
      C(189)= 0.3865502003464D-04
! ... 7.00 < X <= 7.25  INTERVAL NO.29, ABS.ERROR = 1.7034818000639D-
      IFRST(29)=190
      ILAST(29)=194
      C(190)= 0.3781177456714D+00
      C(191)=-0.1101165554215D+00
      C(192)= 0.1580007215659D-01
      C(193)=-0.1125279834923D-02
      C(194)= 0.3192276744812D-04
! ... 7.25 < X <= 7.50  INTERVAL NO.30, ABS.ERROR = 1.3624212868990D-
      IFRST(30)=195
      ILAST(30)=199
      C(195)= 0.3633408065001D+00
      C(196)=-0.1019602554247D+00
      C(197)= 0.1411174419289D-01
      C(198)=-0.9699476170226D-03
      C(199)= 0.2656330616446D-04
! ... 7.50 < X <= 7.75  INTERVAL NO.31, ABS.ERROR = 1.0990319765369D-
      IFRST(31)=200
      ILAST(31)=204
      C(200)= 0.3497439218126D+00
      C(201)=-0.9470569276499D-01
      C(202)= 0.1266017939918D-01
      C(203)=-0.8408550675085D-03
      C(204)= 0.2225784919574D-04
! ... 7.75 < X <= 8.00  INTERVAL NO.32, ABS.ERROR = 8.9381835266522D-
      IFRST(32)=205
      ILAST(32)=209
      C(205)= 0.3371845853069D+00
      C(206)=-0.8822105128463D-01
      C(207)= 0.1140456233287D-01
      C(208)=-0.7327946570967D-03
      C(209)= 0.1877023896668D-04
! ... 8.00 < X <= 8.25  INTERVAL NO.33, ABS.ERROR = 7.3210326689832D-
      IFRST(33)=210
      ILAST(33)=214
      C(210)= 0.3255431584494D+00
      C(211)=-0.8239832407718D-01
      C(212)= 0.1031237441518D-01
      C(213)=-0.6417393667562D-03
      C(214)= 0.1592339503986D-04
! ... 8.25 < X <= 8.50  INTERVAL NO.34, ABS.ERROR = 6.0371707633067D-
      IFRST(34)=215
      ILAST(34)=219
      C(215)= 0.3147180811939D+00
      C(216)=-0.7714810368023D-01
      C(217)= 0.9357439877093D-02
      C(218)=-0.5645414814808D-03
      C(219)= 0.1358301324217D-04
! ... 8.50 < X <= 8.75  INTERVAL NO.35, ABS.ERROR = 5.0113246885530D-
      IFRST(35)=220
      ILAST(35)=224
      C(220)= 0.3046231027165D+00
      C(221)=-0.7239608540250D-01
      C(222)= 0.8518560948886D-02
      C(223)=-0.4987218951555D-03
      C(224)= 0.1164632612927D-04
! ... 8.75 < X <= 9.00  INTERVAL NO.36, ABS.ERROR = 4.1842085352073D-
      IFRST(36)=225
      ILAST(36)=229
      C(225)= 0.2951839772317D+00
      C(226)=-0.6807982796391D-01
      C(227)= 0.7778392737741D-02
      C(228)=-0.4423078514719D-03
      C(229)= 0.1003385659715D-04
! ... 9.00 < X <= 9.25  INTERVAL NO.37, ABS.ERROR = 3.5125236053091D-
      IFRST(37)=230
      ILAST(37)=234
      C(230)= 0.2863364933588D+00
      C(231)=-0.6414655846973D-01
      C(232)= 0.7122648855400D-02
      C(233)=-0.3937177713169D-03
      C(234)= 0.8683627129358D-05
! ... 9.25 < X <= 9.50  INTERVAL NO.38, ABS.ERROR = 2.9645175203541D-
      IFRST(38)=235
      ILAST(38)=239
      C(235)= 0.2780246155388D+00
      C(236)=-0.6055132372497D-01
      C(237)= 0.6539470823348D-02
      C(238)=-0.3516734769846D-03
      C(239)= 0.7546893357357D-05
! ... 9.50 < X <= 9.75  INTERVAL NO.39, ABS.ERROR = 2.5142110615661D-
      IFRST(39)=240
      ILAST(39)=244
      C(240)= 0.2701996594546D+00
      C(241)=-0.5725581832360D-01
      C(242)= 0.6018987049956D-02
      C(243)=-0.3151372568482D-03
      C(244)= 0.6585092705791D-05
! ... 9.75 < X <= 10.0  INTERVAL NO.40, ABS.ERROR = 2.1422863483167D-
      IFRST(40)=245
      ILAST(40)=249
      C(245)= 0.2628187673463D+00
      C(246)=-0.5422707501601D-01
      C(247)= 0.5552907092188D-02
      C(248)=-0.2832594125266D-03
      C(249)= 0.5767453330918D-05
      RETURN
      END

! DCO
      DOUBLE PRECISION FUNCTION DCO(L,M,KX,KY,KZ,LP,MP,                 &
                                    FPQR,ZLM,LMF,LMX,LMY,LMZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION FPQR(25,25,25),ZLM(*)
      DIMENSION LMF(*),LMX(*),LMY(*),LMZ(*)
!
!        -----  ROUTINE EVALUATES ANGULAR MOMENTUM COUPLING     -----
!        -----  COEFFICIENTS.                                   -----
!
! FOR A GIVEN SET OF L AND M, THIS ROUTINE CALCULATES THE TYPE 1
! (L,M=0) OR THE TYPE 2 ANGULAR INTEGRALS. (IE THE SECOND LINE OF
! EQ 28 OR 29 IN MD'S PAPER) YOU STILL MUST DO THE SUMS.
!
      ID = L*(L+1)-M+1
      IMN= LMF(ID)
      IMX= LMF(ID+1)-1
      JD = LP*(LP+1)-MP+1
      JMN= LMF(JD)
      JMX= LMF(JD+1)-1
      SUMI= 0.0D+00
      DO 120 I=IMN,IMX
         IX = LMX(I)+KX+1
         IY = LMY(I)+KY+1
         IZ = LMZ(I)+KZ+1
         SUMJ= 0.0D+00
         DO 110 J=JMN,JMX
            JX = IX+LMX(J)
            JY = IY+LMY(J)
            JZ = IZ+LMZ(J)
  110    SUMJ= SUMJ+ZLM(J)*FPQR(JX,JY,JZ)
  120 SUMI= SUMI+ZLM(I)*SUMJ
      DCO = SUMI
      RETURN
      END

! ERRT
      SUBROUTINE ERRT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /ERRFCM/ C(142),IFRST(20),ILAST(20),H
!
!        -----  ROUTINE ALLOCATES Parameters FOR THE POLYNOMIAL -----
!        -----  FIT TO THE ERROR FUNCTION.                      -----
!
      H = 0.25D+00
! ... 0.00 < X <= 0.25  INTERVAL NO. 1, ABS.ERROR = 5.1496584774213D-
      IFRST( 1)= 1
      ILAST( 1)= 8
      C( 1)= 0.9976905682968D-15
      C( 2)= 0.1128379167616D+01
      C( 3)=-0.6825080083050D-07
      C( 4)=-0.3761237879247D+00
      C( 5)=-0.4362679036164D-04
      C( 6)= 0.1132081528735D+00
      C( 7)=-0.1608829014003D-02
      C( 8)=-0.2383745887450D-01
! ... 0.25 < X <= 0.50  INTERVAL NO. 2, ABS.ERROR = 1.0329515021112D-
      IFRST( 2)= 9
      ILAST( 2)=16
      C( 9)=-0.3265133599788D-05
      C(10)= 0.1128453123654D+01
      C(11)=-0.7236769295451D-03
      C(12)=-0.3721384745633D+00
      C(13)=-0.1348362056472D-01
      C(14)= 0.1412226404729D+00
      C(15)=-0.3539596498013D-01
      C(16)=-0.5454029355730D-02
! ... 0.50 < X <= 0.75  INTERVAL NO. 3, ABS.ERROR = 6.1888272284705D-
      IFRST( 3)=17
      ILAST( 3)=24
      C(17)=-0.1801988094647D-03
      C(18)= 0.1130782166249D+01
      C(19)=-0.1400474656066D-01
      C(20)=-0.3295493300819D+00
      C(21)=-0.9653664285517D-01
      C(22)= 0.2398163462058D+00
      C(23)=-0.1014115256923D+00
      C(24)= 0.1378314835685D-01
! ... 0.75 < X <= 1.00  INTERVAL NO. 4, ABS.ERROR = 1.6626700016786D-
      IFRST( 4)=25
      ILAST( 4)=32
      C(25)=-0.4457991000617D-03
      C(26)= 0.1133670008174D+01
      C(27)=-0.2725712863074D-01
      C(28)=-0.2961333291272D+00
      C(29)=-0.1466850258377D+00
      C(30)= 0.2847021965842D+00
      C(31)=-0.1236317583493D+00
      C(32)= 0.1848162923540D-01
! ... 1.00 < X <= 1.25  INTERVAL NO. 5, ABS.ERROR = 5.4747317790316D-
      IFRST( 5)=33
      ILAST( 5)=40
      C(33)= 0.1036667435991D-01
      C(34)= 0.1060681450430D+01
      C(35)= 0.1844323965440D+00
      C(36)=-0.6381386282176D+00
      C(37)= 0.1857801550899D+00
      C(38)= 0.9020967608584D-01
      C(39)=-0.6022465654782D-01
      C(40)= 0.9593725204468D-02
! ... 1.25 < X <= 1.50  INTERVAL NO. 6, ABS.ERROR = 2.1689317009077D-
      IFRST( 6)=41
      ILAST( 6)=47
      C(41)= 0.5096577645657D-01
      C(42)= 0.8301050352144D+00
      C(43)= 0.7461454169558D+00
      C(44)=-0.1399088183681D+01
      C(45)= 0.8049290742492D+00
      C(46)=-0.2123844940215D+00
      C(47)= 0.2202851573626D-01
! ... 1.50 < X <= 1.75  INTERVAL NO. 7, ABS.ERROR = 5.4356519285648D-
      IFRST( 7)=48
      ILAST( 7)=55
      C(48)= 0.1356983343156D+00
      C(49)= 0.4248183601853D+00
      C(50)= 0.1578617255285D+01
      C(51)=-0.2350907279473D+01
      C(52)= 0.1459137006763D+01
      C(53)=-0.4826721286933D+00
      C(54)= 0.8417719176837D-01
      C(55)=-0.6134748458862D-02
! ... 1.75 < X <= 2.00  INTERVAL NO. 8, ABS.ERROR = 1.3677947663382D-
      IFRST( 8)=56
      ILAST( 8)=63
      C(56)= 0.5650617287234D-01
      C(57)= 0.7303727164846D+00
      C(58)= 0.1073427393133D+01
      C(59)=-0.1886956498357D+01
      C(60)= 0.1203543744299D+01
      C(61)=-0.3982083398317D+00
      C(62)= 0.6867500288146D-01
      C(63)=-0.4915782383510D-02
! ... 2.00 < X <= 2.25  INTERVAL NO. 9, ABS.ERROR = 2.7604585284280D-
      IFRST( 9)=64
      ILAST( 9)=70
      C(64)=-0.7044612080369D+00
      C(65)= 0.3313423632378D+01
      C(66)=-0.2688303271381D+01
      C(67)= 0.1159938873721D+01
      C(68)=-0.2789363703341D+00
      C(69)= 0.3510668004553D-01
      C(70)=-0.1778023938338D-02
! ... 2.25 < X <= 2.50  INTERVAL NO.10, ABS.ERROR = 1.8687273950491D-
      IFRST(10)=71
      ILAST(10)=77
      C(71)=-0.8530571392369D+00
      C(72)= 0.3718753089373D+01
      C(73)=-0.3149110592528D+01
      C(74)= 0.1439416506820D+01
      C(75)=-0.3743063133443D+00
      C(76)= 0.5246810708195D-01
      C(77)=-0.3095234433810D-02
! ... 2.50 < X <= 2.75  INTERVAL NO.11, ABS.ERROR = 9.7095664841618D-
      IFRST(11)=78
      ILAST(11)=84
      C(78)=-0.6325006110506D+00
      C(79)= 0.3195352274705D+01
      C(80)=-0.2631359977752D+01
      C(81)= 0.1166146516443D+01
      C(82)=-0.2931400812425D+00
      C(83)= 0.3960487060249D-01
      C(84)=-0.2245448529720D-02
! ... 2.75 < X <= 3.00  INTERVAL NO.12, ABS.ERROR = 7.6525452641363D-
      IFRST(12)=85
      ILAST(12)=91
      C(85)=-0.1578663656307D+00
      C(86)= 0.2160403160352D+01
      C(87)=-0.1690605397678D+01
      C(88)= 0.7098560625163D+00
      C(89)=-0.1685917225356D+00
      C(90)= 0.2146466169506D-01
      C(91)=-0.1144051551819D-02
! ... 3.00 < X <= 3.25  INTERVAL NO.13, ABS.ERROR = 3.9079850466806D-
      IFRST(13)=92
      ILAST(13)=98
      C(92)= 0.3337223693570D+00
      C(93)= 0.1174814163831D+01
      C(94)=-0.8669311383990D+00
      C(95)= 0.3425820853154D+00
      C(96)=-0.7643589524863D-01
      C(97)= 0.9127119866510D-02
      C(98)=-0.4555632670720D-03
! ... 3.25 < X <= 3.50  INTERVAL NO.14, ABS.ERROR = 1.5241141682054D-
      IFRST(14)= 99
      ILAST(14)=105
      C( 99)= 0.6858505344288D+00
      C(100)= 0.5221032227938D+00
      C(101)=-0.3626430750192D+00
      C(102)= 0.1347162453894D+00
      C(103)=-0.2822370117065D-01
      C(104)= 0.3161233228942D-02
      C(105)=-0.1478642225266D-03
! ... 3.50 < X <= 3.75  INTERVAL NO.15, ABS.ERROR = 9.5496943686157D-
      IFRST(15)=106
      ILAST(15)=111
      C(106)= 0.9680982760858D+00
      C(107)= 0.4187846825891D-01
      C(108)=-0.2202918541743D-01
      C(109)= 0.5803572496370D-02
      C(110)=-0.7656564703211D-03
      C(111)= 0.4046317189932D-04
! ... 3.75 < X <= 4.00  INTERVAL NO.16, ABS.ERROR = 2.2239987629291D-
      IFRST(16)=112
      ILAST(16)=117
      C(112)= 0.9910262142894D+00
      C(113)= 0.1110389476297D-01
      C(114)=-0.5503022341054D-02
      C(115)= 0.1365301693295D-02
      C(116)=-0.1695606159046D-03
      C(117)= 0.8432380855083D-05
! ... 4.00 < X <= 4.25  INTERVAL NO.17, ABS.ERROR = 4.4764192352886D-
      IFRST(17)=118
      ILAST(17)=123
      C(118)= 0.9978511542662D+00
      C(119)= 0.2512562013763D-02
      C(120)=-0.1176282812867D-02
      C(121)= 0.2755980203801D-03
      C(122)=-0.3231375012547D-04
      C(123)= 0.1516751945019D-05
! ... 4.25 < X <= 4.50  INTERVAL NO.18, ABS.ERROR = 9.9475983006414D-
      IFRST(18)=124
      ILAST(18)=129
      C(124)= 0.9995577955752D+00
      C(125)= 0.4898449770963D-03
      C(126)=-0.2172069044718D-03
      C(127)= 0.4819046007469D-04
      C(128)=-0.5349400453269D-05
      C(129)= 0.2376735210419D-06
! ... 4.50 < X <= 4.75  INTERVAL NO.19, ABS.ERROR = 6.0396132539609D-
      IFRST(19)=130
      ILAST(19)=135
      C(130)= 0.9999195780401D+00
      C(131)= 0.8461205943888D-04
      C(132)=-0.3562696920198D-04
      C(133)= 0.7504334644182D-05
      C(134)=-0.7907161489129D-06
      C(135)= 0.3334134817123D-07
! ... 4.75 < X <= 5.00  INTERVAL NO.20, ABS.ERROR = 7.4429351570870D-
      IFRST(20)=136
      ILAST(20)=142
      C(136)= 0.8965357422439D+00
      C(137)= 0.1272754720567D+00
      C(138)=-0.6523024655166D-01
      C(139)= 0.1782843628219D-01
      C(140)=-0.2740698420287D-02
      C(141)= 0.2246825024486D-03
      C(142)=-0.7674098014832D-05
      RETURN
      END

! DAWERF
      DOUBLE PRECISION FUNCTION DAWERF(Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DERFCM/ C(246),IFRST(40),ILAST(40),H
!     EVALUATES THE HYBRID DAWSON-ERROR FUNCTION
      PARAMETER (ZP5=0.5D0, TEN=10.0D0)
      X= ABS(Y)
      IF(X>=TEN) GO TO 115

      XN = X/H
      NX = INT(XN)+1
      I1 = IFRST(NX)
      IL = ILAST(NX)
      T= C(IL)
      DO 110 K=IL-1,I1,-1
  110 T= C(K)+T*X
      DAWERF= T
      RETURN
  115 CONTINUE
      X2I = 1.0d0/(X*X)
      IL = 8
      TXT=(IL +ZP5)*X2I
      SUM= TXT
      DO 120 K=IL,1,-1
         TXT= TXT-X2I
  120 SUM= TXT*(1.0d0+SUM)
      DAWERF= X*TXT*(1.0d0+SUM)
      RETURN
      END
      
! DAWF
      DOUBLE PRECISION FUNCTION DAWF(Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DAWFCM/ C(249),IFRST(40),ILAST(40),H
!     EVALUATES THE DAWSON FUNCTION -----
      PARAMETER (ZP5=0.5D0, TEN=10.0D0)
!
      X= Y
      S= 1.0d0
      IF(X>=0.0d0) GO TO 105
      X=-Y
      S=-1.0d0
  105 IF(X>=TEN) GO TO 115
      XN = X/H
      NX = INT(XN)+1
      I1 = IFRST(NX)
      IL = ILAST(NX)
      T= C(IL)
      DO 110 K=IL-1,I1,-1
  110 T= C(K)+T*X
      DAWF  = T*S
      RETURN
!
  115 CONTINUE
      X2I = 1.0d0/(X*X)
      IL = 8
      TXT=(IL +ZP5)*X2I
      SUM= TXT
      DO 120 K=IL,1,-1
         TXT= TXT-X2I
  120 SUM= TXT*(1.0d0+SUM)
      DAWF = Y*TXT*(1.0d0+SUM)
      RETURN
      END
      
! ECPAA1
      SUBROUTINE ECPAA1(NPNP,FPQR,COEFI,COEFJ,DCOEF4,G,EX,NPRIMI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DCOEF4(*),FPQR(25,25,25)
      DIMENSION COEFI(*),COEFJ(*),G(*)
      COMMON /ECP1  / X01,CAX,CAY,CAZ,CA,XCA,YCA,ZCA,                   &
                      X02,BAX,BAY,BAZ,BA,XBA,YBA,ZBA,                   &
                      PHASE,DAX,DAY,DAZ,DA,XDA,YDA,ZDA,XINT,KCNTR        
      !JFHLewYee: Changed NATOMS allowed dimension from 100 to 1000
      COMMON /ECP2  / CLP(4004),ZLP(4004),NLP(4004),KFRST(1001,6),      &
                      KLAST(1001,6),LMAX(1001),LPSKIP(1001),IZCORE(1001)
      COMMON /ECPDIM/ NCOEF1,NCOEF2,J1LEN,J2LEN,LLIM,NLIM,NTLIM,J4LEN
      LOGICAL IANDJ,NORM,NORMI,NORMJ
      COMMON /ECPIDX/ Q2,IAMIN,IAMAX,JAMIN,JAMAX,IPMIN,IPMAX,JPMIN,     &
                      JPMAX,KF1,KL1,LLMX,NPC,NPB,IANDJ,NORM,NORMI,NORMJ                                                                        
      COMMON /GBASE / NFST(8),NX(84),NY(84),NZ(84)                       
      DIMENSION EX(NPRIMI)
      PARAMETER (SR3=1.73205080756888D+00, SR5=2.23606797749979D+00)     
      PARAMETER (SR7=2.64575131106459D+00)                               
      PARAMETER (S35=SR3*SR5, S57=SR5*SR7, S53=S57/SR3)                                                                        
      DIMENSION ECPFAC(35)                                               
      DATA ECPFAC/7*1.0d0,3*SR3,3*1.0d0,6*SR5,S35,3*1.0d0,6*SR7,        &
                  3*S53,3*S57/
!
      NBC = NPNP-1
      IF(MOD(NBC,2)==1) RETURN
      FPSUM = 0.0d0
      DO 120 IG=IPMIN,IPMAX
         X01= EX(IG)
         DO 120 JG=JPMIN,JPMAX
            X02= EX(JG)
            X12= X01+X02
            DUM= COEFI(IG-IPMIN+1)*COEFJ(JG-JPMIN+1)
            FPTEMP= 0.0d0
            DO 110 K=KF1,KL1
  110       FPTEMP= FPTEMP+FA(NBC+NLP(K),X12+ZLP(K))*CLP(K)
  120 FPSUM = FPSUM+FPTEMP*DUM
      FPSQ2 = FPSUM*Q2
      MMAX= JAMAX
      NN= 1
      DO 140 I=IAMIN,IAMAX
         DUMI= FPSQ2
         IF(NORMI) DUMI= DUMI*ECPFAC(I)
         IF(IANDJ) MMAX= I
         DO 140 J=JAMIN,MMAX
            DUMJ= DUMI
            IF(NORMJ) DUMJ= DUMJ*ECPFAC(J)
            MX= NX(I)+NX(J)+1
            MY= NY(I)+NY(J)+1
            MZ= NZ(I)+NZ(J)+1
            G(NN)= G(NN)+DUMJ*FPQR(MX,MY,MZ)
  140 NN= NN+1

      IF(LLMX<2) RETURN
      DO 300 LL=2,LLMX
         L = LL-2
         NLM1=(L*(L+1))
         KF = KFRST(KCNTR,LL)
         KL = KLAST(KCNTR,LL)
         FPSUM = 0.0d0
         DO 220 IG=IPMIN,IPMAX
            X01= EX(IG)
            DO 220 JG=JPMIN,JPMAX
               X02= EX(JG)
               X12= X01+X02
               DUM= COEFI(IG-IPMIN+1)*COEFJ(JG-JPMIN+1)
               FPTEMP= 0.0d0
               DO 210 K=KF,KL
  210          FPTEMP= FPTEMP+FA(NBC+NLP(K),X12+ZLP(K))*CLP(K)
  220    FPSUM = FPSUM+FPTEMP*DUM
         FPSQ2 = FPSUM*Q2
         MMAX= JAMAX
         NN= 1
         DO 240 I=IAMIN,IAMAX
            DUMI= FPSQ2
            IF(NORMI) DUMI= DUMI*ECPFAC(I)
            IF(IANDJ) MMAX= I
            DO 240 J=JAMIN,MMAX
               DUMJ= DUMI
               IF(NORMJ) DUMJ= DUMJ*ECPFAC(J)
               SUM= 0.0d0
               DO 230 M=-L,L
                  NLM=(NLM1-M)*NTLIM
  230          SUM= SUM+DCOEF4(NLM+I)*DCOEF4(NLM+J)
               G(NN)= G(NN)+DUMJ*SUM
  240    NN= NN+1

  300 CONTINUE
      RETURN
      END

! ECPRA2
      SUBROUTINE ECPRA2(ICAB,NPNP,FP,COEFI,COEFJ,DCOEF1,JFST1,LBECP1,   &
          DCOEF4,DCOEF2,JFST2,LBECP2,G,ZLM,LMF,LMX,LMY,LMZ,EX,NPRIMI)                             
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DOUBLE PRECISION DCOEF1(*),DCOEF2(*),DCOEF4(*),ZLM(*)
      DIMENSION FP(*),COEFI(*),COEFJ(*),JFST1(*),LBECP1(9,*),           &
                JFST2(*),LBECP2(6,*),G(*)                      
      DIMENSION LMF(*),LMX(*),LMY(*),LMZ(*)                       
      COMMON /ECP1  / X01,CAX,CAY,CAZ,CA,XCA,YCA,ZCA,                   &
                      X02,BAX,BAY,BAZ,BA,XBA,YBA,ZBA,                   &
                      PHASE,DAX,DAY,DAZ,DA,XDA,YDA,ZDA,XINT,KCNTR        
      !JFHLewYee: Changed NATOMS allowed dimension from 100 to 1000
      COMMON /ECP2  / CLP(4004),ZLP(4004),NLP(4004),KFRST(1001,6),      &
                      KLAST(1001,6),LMAX(1001),LPSKIP(1001),IZCORE(1001)                                        
      COMMON /ECP3  / ACX(7),ACY(7),ACZ(7),ABX(7),ABY(7),ABZ(7)          
      COMMON /ECPDIM/ NCOEF1,NCOEF2,J1LEN,J2LEN,LLIM,NLIM,NTLIM,J4LEN    
      LOGICAL IANDJ,NORM,NORMI,NORMJ                                     
      COMMON /ECPIDX/ Q2,IAMIN,IAMAX,JAMIN,JAMAX,IPMIN,IPMAX,JPMIN,     &
                      JPMAX,KF1,KL1,LLMX,NPC,NPB,IANDJ,NORM,NORMI,NORMJ  
      COMMON /FICMN / ALF,XI,XP0,XP1                                     
      COMMON /GBASE / NFST(8),NX(84),NY(84),NZ(84)
      
      DIMENSION EX(NPRIMI)               
      PARAMETER (SR3=1.73205080756888D+00, SR5=2.23606797749979D+00)     
      PARAMETER (SR7=2.64575131106459D+00)                               
      PARAMETER (S35=SR3*SR5, S57=SR5*SR7, S53=S57/SR3)                  
      PARAMETER (TOL=1.0D-10, FPI=12.566370614359D+00)                   
      DIMENSION ECPFAC(35),FIP(78),ZFNLM(121),CKL(11,11)                 
      DATA ECPFAC/7*1.0d0,3*SR3,3*1.0d0,6*SR5,S35,3*1.0d0,6*SR7,        &
                  3*S53,3*S57/
!     
      NPNPMX=(NPNP*(NPNP+1))/2
      CALL VCLR(FP,1,NPNPMX)

      DO 140 IG=IPMIN,IPMAX
         X01= EX(IG)
         DO 130 JG=JPMIN,JPMAX
            X02= EX(JG)
            X12= X01+X02
            DUM= COEFI(IG-IPMIN+1)*COEFJ(JG-JPMIN+1)
            IF(ICAB==2) THEN
               ALFA= X02*BA
               ALFX= ALFA*BA
            ELSE IF(ICAB==3) THEN
               ALFA= X01*CA
               ALFX= ALFA*CA
            END IF
            ALFI= 1.0d0/(ALFA+ALFA)
            XP0 = EXP(-ALFX)
            DO 120 K=KF1,KL1
               XI = 1.0d0/SQRT(X12+ZLP(K))
               ALF= ALFA*XI
               XP1= EXP(-ALFX+ALF*ALF)
               NLPK= NLP(K)
               CLPK= CLP(K)*DUM
               CALL FIECP(FIP,ALFI,NLPK,NPNP-1)
               DO 110 N=1,NPNPMX
  110          FP(N)= FP(N)+FIP(N)*CLPK
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
      CALL ZFN(ZFNLM,NPNP-1,ZLM,LMF,LMX,LMY,LMZ)
      FPIQ2 = FPI*Q2
      MMAX= JAMAX
      NN= 1
      DO 240 I=IAMIN,IAMAX
         DUMI= FPIQ2
         IF(NORMI) DUMI= DUMI*ECPFAC(I)
         IF(IANDJ) MMAX= I
         DO 240 J=JAMIN,MMAX
            DUMJ= DUMI
            IF(NORMJ) DUMJ= DUMJ*ECPFAC(J)
            IF(I>=J) THEN
               INDX= J+(I*(I-1))/2
            ELSE
               INDX= I+(J*(J-1))/2
            END IF
            JF = JFST1(INDX)
            JL = JFST1(INDX+1)-1
            DO 210 L=1,NPNP
               DO 210 K=L,NPNP
  210       CKL(K,L)= 0.0d0
            NXYZ = NX(I)+NY(I)+NZ(I)+NX(J)+NY(J)+NZ(J)
            K= I
            M= 3
            IF(ICAB==2) THEN
               K= I+J-K
               M= 3+6-M
            END IF
            IF(I>=J) M= 3+6-M
            DO 220 JJ=JF,JL
               KA = LBECP1(1,JJ)+1
               LA = LBECP1(2,JJ)+1
               MU = LBECP1(3,JJ)
               KX = LBECP1(M+1,JJ)+NX(K)
               KY = LBECP1(M+2,JJ)+NY(K)
               KZ = LBECP1(M+3,JJ)+NZ(K)
               IF((KX+KY+KZ)/=NXYZ) GO TO 220
               IX = LBECP1(4,JJ)+LBECP1(7,JJ)
               IY = LBECP1(5,JJ)+LBECP1(8,JJ)
               IZ = LBECP1(6,JJ)+LBECP1(9,JJ)
               LU = LA*(LA-1)-MU+1
               IF(ICAB==2) THEN
                  DT= ABX(KX+1-IX)*ABY(KY+1-IY)*ABZ(KZ+1-IZ)*ZFNLM(LU)
               ELSE
                  DT= ACX(KX+1-IX)*ACY(KY+1-IY)*ACZ(KZ+1-IZ)*ZFNLM(LU)
               END IF
               CKL(KA,LA)= CKL(KA,LA)+DCOEF1(JJ)*DT
  220       CONTINUE
            SUM= 0.0d0
            N= 1
            DO 230 K=1,NPNP
               DO 230 L=1,K
                  IF(ABS(CKL(K,L))>TOL) SUM= SUM+FP(N)*CKL(K,L)
  230       N= N+1
            G(NN)= G(NN)+DUMJ*SUM
  240 NN= NN+1
      IF(LLMX<2) RETURN
      DO 500 LL=2,LLMX
         CALL VCLR(FP,1,NPNPMX)
         KF = KFRST(KCNTR,LL)
         KL = KLAST(KCNTR,LL)
         DO 340 IG=IPMIN,IPMAX
            X01= EX(IG)
            DO 330 JG=JPMIN,JPMAX
               X02= EX(JG)
               X12= X01+X02
               DUM= COEFI(IG-IPMIN+1)*COEFJ(JG-JPMIN+1)
               IF(ICAB==2) THEN
                  GAMA= X02*BA
                  GAMX= GAMA*BA
               ELSE IF(ICAB==3) THEN
                  GAMA= X01*CA
                  GAMX= GAMA*CA
               END IF
               GAMI= 1.0d0/(GAMA+GAMA)
               XP0 = EXP(-GAMX)
               DO 320 K=KF,KL
                  XI = 1.0d0/SQRT(X12+ZLP(K))
                  ALF= GAMA*XI
                  XP1= EXP(-GAMX+ALF*ALF)
                  NLPK = NLP(K)
                  CLPK = CLP(K)*DUM
                  CALL FIECP(FIP,GAMI,NLPK,NPNP-1)
                  DO 310 N=1,NPNPMX
  310             FP(N)= FP(N)+FIP(N)*CLPK
  320          CONTINUE
  330       CONTINUE
  340    CONTINUE
         L2PL=(LL-2)*(LL-1)
         CALL ZFN(ZFNLM,NPNP-1,ZLM,LMF,LMX,LMY,LMZ)
         FPIQ2 = FPI*Q2
         MMAX= JAMAX
         NN= 1
         DO 440 I=IAMIN,IAMAX
            DUMI= FPIQ2
            IF(NORMI) DUMI= DUMI*ECPFAC(I)
            IF(IANDJ) MMAX= I
            DO 440 J=JAMIN,MMAX
               DUMJ= DUMI
               IF(NORMJ) DUMJ= DUMJ*ECPFAC(J)
               DO 400 L=1,NPNP
                  DO 400 K=L,NPNP
  400          CKL(K,L)= 0.0d0
               MX = NX(I)+NX(J)+1
               MY = NY(I)+NY(J)+1
               MZ = NZ(I)+NZ(J)+1
               IF(ICAB==2) THEN
                  K= I
               ELSE
                  K= J
               END IF
               KT = NX(K)+NY(K)+NZ(K)+1
               L= LL-2
               DO 420 M=L,-L,-1
                  NLM=(L2PL-M)*NTLIM
                  CT = DCOEF4(NLM+K)
                  IF(ABS(CT)<TOL) GO TO 420
                  NLM= NLM+I+J
                  JF = JFST2(NLM-K)
                  JL = JFST2(NLM-K+1)-1
                  DO 410 JJ=JF,JL
                     LA = LBECP2(1,JJ)+1
                     KA = LBECP2(2,JJ)+KT
                     MU = LBECP2(3,JJ)
                     KX = LBECP2(4,JJ)+NX(K)
                     KY = LBECP2(5,JJ)+NY(K)
                     KZ = LBECP2(6,JJ)+NZ(K)
                     LU = LA*(LA-1)-MU+1
                     IF(ICAB==2) THEN
                        DT=CT*ABX(MX-KX)*ABY(MY-KY)*ABZ(MZ-KZ)*ZFNLM(LU)
                     ELSE
                        DT=CT*ACX(MX-KX)*ACY(MY-KY)*ACZ(MZ-KZ)*ZFNLM(LU)
                     END IF
                     CKL(KA,LA)= CKL(KA,LA)+DCOEF2(JJ)*DT
  410             CONTINUE
  420          CONTINUE
               SUM= 0.0d0
               N= 1
               DO 430 K=1,NPNP
                  DO 430 L=1,K
               IF(ABS(CKL(K,L))>TOL) SUM= SUM+FP(N)*CKL(K,L)
  430          N= N+1
               G(NN)= G(NN)+DUMJ*SUM
  440    NN= NN+1
  500 CONTINUE
!
      RETURN
      END

! ECPDRA
      SUBROUTINE ECPDRA(IC4C,NPNP,FP,FQ,COEFI,COEFQ,COEFJ,DCOEF1,JFST1, &
            LBECP1,DCOEF2,JFST2,LBECP2,G,ZLM,LMF,LMX,LMY,LMZ,EX,NPRIMI)                             
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DOUBLE PRECISION DCOEF1(*),DCOEF2(*),ZLM(*)
      DIMENSION FP(*),FQ(*),COEFI(*),COEFQ(*),COEFJ(*),JFST1(*),        &
                LBECP1(9,*),JFST2(*),LBECP2(6,*),G(*)                      
      DIMENSION LMF(*),LMX(*),LMY(*),LMZ(*)                       
      COMMON /ECP1  / X01,CAX,CAY,CAZ,CA,XCA,YCA,ZCA,                   &
                      X02,BAX,BAY,BAZ,BA,XBA,YBA,ZBA,                   &
                      PHASE,DAX,DAY,DAZ,DA,XDA,YDA,ZDA,XINT,KCNTR        
      !JFHLewYee: Changed NATOMS allowed dimension from 100 to 1000
      COMMON /ECP2  / CLP(4004),ZLP(4004),NLP(4004),KFRST(1001,6),      &
                      KLAST(1001,6),LMAX(1001),LPSKIP(1001),IZCORE(1001)                                        
      COMMON /ECP3  / ACX(7),ACY(7),ACZ(7),ABX(7),ABY(7),ABZ(7)          
      LOGICAL CANDB                                                      
      COMMON /ECP4  / P12(3,2),R12,ACO(3),CANDB                          
      COMMON /ECPDIM/ NCOEF1,NCOEF2,J1LEN,J2LEN,LLIM,NLIM,NTLIM,J4LEN    
      LOGICAL IANDJ,NORM,NORMI,NORMJ                                     
      COMMON /ECPIDX/ Q2,IAMIN,IAMAX,JAMIN,JAMAX,IPMIN,IPMAX,JPMIN,     &
                      JPMAX,KF1,KL1,LLMX,NPC,NPB,IANDJ,NORM,NORMI,NORMJ  
      COMMON /ZFNCM / X,Y,Z                                              
      COMMON /FICMN / ALF,XI,XP0,XP1                                     
      COMMON /FJCMN / ALEF,BEIT,XXI,XPLS,XMNS,XP                         
      COMMON /FJNEW / XKA,XKB,GAM1,GAM2,A1,A2,C                          
      COMMON /GBASE / NFST(8),NX(84),NY(84),NZ(84)                       
      
      DIMENSION EX(NPRIMI)
      DIMENSION FIP(78),ZFNLM(121),CKLU(23,12,12)                        
      DIMENSION CKLC(11,11),CKLB(11,11),CPQ(11,11,11)                    
      DIMENSION FJPQ(11,11,11),ZFNLMC(121),ZFNLMB(121)                   
      SAVE FJPQ,ZFNLMC,ZFNLMB                                            
      PARAMETER (ABLIM=1.0D-01)                                          
      PARAMETER (SR3=1.73205080756888D+00, SR5=2.23606797749979D+00)     
      PARAMETER (SR7=2.64575131106459D+00)                               
      PARAMETER (S35=SR3*SR5, S57=SR5*SR7, S53=S57/SR3)                  
      PARAMETER (TOL=1.0D-10, FPI=12.566370614359D+00)                   
      PARAMETER (TM6=1.0D-06, ONDS4P=0.28209479177388D+00)               
      DIMENSION ECPFAC(35)                                               
      DATA ECPFAC/7*1.0d0,3*SR3,3*1.0d0,6*SR5,S35,3*1.0d0,6*SR7,        &
                  3*S53,3*S57/
!
      IF(IC4C>0) GO TO 200
      NPNPMX=(NPNP*(NPNP+1)*(2*NPNP+1))/6
      CALL VCLR(FP,1,NPNPMX)
      IF(IC4C<0)CALL VCLR(FQ,1,NPNPMX)
      DO 160 IG=IPMIN,IPMAX
         X01= EX(IG)
         DO 160 JG=JPMIN,JPMAX
            X02= EX(JG)
            X12= X01+X02
            DUM= COEFI(IG-IPMIN+1)*COEFJ(JG-JPMIN+1)
            IF(IC4C<0)DUQ= COEFQ(IG-IPMIN+1)*COEFJ(JG-JPMIN+1)
            X21= 1.0d0/X12
            Y01= X01*X21
            Y02= 1.0d0-Y01
            Y12= Y01*X02
            DAX= P12(1,1)+(P12(1,2)-P12(1,1))*Y02-ACO(1)
            DAY= P12(2,1)+(P12(2,2)-P12(2,1))*Y02-ACO(2)
            DAZ= P12(3,1)+(P12(3,2)-P12(3,1))*Y02-ACO(3)
            DA = SQRT(DAX*DAX+DAY*DAY+DAZ*DAZ)
            IF(.NOT.CANDB) THEN
               DUM= DUM*EXP(-R12*Y12)
               IF(IC4C<0)DUQ= DUQ*EXP(-R12*Y12)
            END IF
            IF(DA>=TM6) THEN
               ALFA= X12*DA
               ALFX= ALFA*DA
               ALFI= 1.0d0/(ALFA+ALFA)
               X  = DAX/DA
               Y  = DAY/DA
               Z  = DAZ/DA
               CALL ZFN(ZFNLM,NPNP-1,ZLM,LMF,LMX,LMY,LMZ)
               XP0 = EXP(-ALFX)
               DO 130 K=KF1,KL1
                  XI  = 1.0d0/SQRT(X12+ZLP(K))
                  ALF = ALFA*XI
                  XP1 = EXP(-ALFX+ALF*ALF)
                  NLPK= NLP(K)
                  CLPK= CLP(K)*DUM
                  IF(IC4C<0)CLPQ= CLP(K)*DUQ
                  CALL FIECP(FIP,ALFI,NLPK,NPNP-1)
                  N= 1
                  NN= 1
                  DO 120 KK=1,NPNP
                     DO 120 L=1,KK
                        FIPTEM= FIP(NN)*CLPK
                        IF(IC4C<0)FIPTEQ= FIP(NN)*CLPQ
                        KKLL= L*(L-1)+1+L
                        DO 110 MU=1,2*L-1
                           FP(N)= FP(N)+FIPTEM*ZFNLM(KKLL-MU)
                           IF(IC4C<0)FQ(N)= FQ(N)+FIPTEQ*ZFNLM(KKLL-MU)
  110                   N= N+1
  120             NN= NN+1
  130          CONTINUE
            ELSE
               DO 150 K=KF1,KL1
                  ZETA= X12+ZLP(K)
                  NLPK= NLP(K)
                  CLPK= CLP(K)*ONDS4P*DUM
                  IF(IC4C<0)CLPQ= CLP(K)*ONDS4P*DUQ
                  DO 140 N=1,NPNP
                     NK=(N*(N-1)*(2*N-1))/6+1
                     FATEMP= FA(N+NLPK-1,ZETA)
                     FP(NK)= FP(NK)+FATEMP*CLPK
                     IF(IC4C<0)FQ(NK)= FQ(NK)+FATEMP*CLPQ
  140             CONTINUE
  150          CONTINUE
            END IF
  160 CONTINUE
  200 CONTINUE
!
      FPIQ2 = FPI*Q2
      MMAX= JAMAX
      NN= 1
      DO 240 I=IAMIN,IAMAX
         DUMI= FPIQ2
         IF(NORMI) DUMI= DUMI*ECPFAC(I)
         IF(IANDJ) MMAX= I
         DO 240 J=JAMIN,MMAX
            DUMJ= DUMI
            IF(NORMJ) DUMJ= DUMJ*ECPFAC(J)
            IF(I>=J) THEN
               INDX= J+(I*(I-1))/2
            ELSE
               INDX= I+(J*(J-1))/2
            END IF
            JF = JFST1(INDX)
            JL = JFST1(INDX+1)-1
            NXC= NX(I)+1
            NYC= NY(I)+1
            NZC= NZ(I)+1
            NXB= NX(J)+1
            NYB= NY(J)+1
            NZB= NZ(J)+1
            DO 210 K=1,NPNP+1
               DO 210 L=1,K
                  DO 210 MU=1,2*L-1
  210       CKLU(MU,L,K)= 0.0d0

            DO 220 K=JF,JL
               KA = LBECP1(1,K)+1
               LA = LBECP1(2,K)+1
               MU = LBECP1(3,K)
               IF(I>=J) THEN
                  KX = LBECP1(4,K)
                  KY = LBECP1(5,K)
                  KZ = LBECP1(6,K)
                  KXP= LBECP1(7,K)
                  KYP= LBECP1(8,K)
                  KZP= LBECP1(9,K)
               ELSE
                  KXP= LBECP1(4,K)
                  KYP= LBECP1(5,K)
                  KZP= LBECP1(6,K)
                  KX = LBECP1(7,K)
                  KY = LBECP1(8,K)
                  KZ = LBECP1(9,K)
               END IF
      CKLU(LA+MU,LA,KA) = CKLU(LA+MU,LA,KA) + DCOEF1(K)                 &
                        * ACX(NXC-KX) *ACY(NYC-KY) *ACZ(NZC-KZ)         &
                        * ABX(NXB-KXP)*ABY(NYB-KYP)*ABZ(NZB-KZP)
  220       CONTINUE
            SUM= 0.0d0
            N= 1
            DO 230 K=1,NPNP
               DO 230 L=1,K
                  DO 230 MU=1,2*L-1
                     CKLTEM= CKLU(MU,L,K)
                     IF(IC4C<=0) THEN
                        IF(ABS(CKLTEM)>TOL) SUM= SUM+FP(N)*CKLTEM
                     ELSE
                        IF(ABS(CKLTEM)>TOL) SUM= SUM+FQ(N)*CKLTEM
                     END IF
  230       N= N+1
            G(NN)= G(NN)+DUMJ*SUM
  240 NN= NN+1
      IF(LLMX<2) RETURN
      FPISQ2= FPI*FPIQ2
      DO 500 LL=2,LLMX
         L= LL-2
         NP1  = MAX0(NPC,NPB)
         NP1PL= NP1+L
         NPCPL= NPC+L
         NPBPL= NPB+L
         NMAX =(NPNP-1)*NPCPL*(NPBPL/2+1)+6
         LTMAX= MAX0(NPNP,NP1PL)
         LEMX = MAX0(L,LTMAX/2)
         KF = KFRST(KCNTR,LL)
         KL = KLAST(KCNTR,LL)
         CALL VCLR(FP,1,NMAX)
         DO 340 IG=IPMIN,IPMAX
            X01= EX(IG)
            ALFA= X01*CA
            ALFX= ALFA*CA
            ALFI= 1.0d0/(ALFA+ALFA)
            DO 340 JG=JPMIN,JPMAX
               X02= EX(JG)
               X12= X01+X02
               BETA= X02*BA
               BETX= BETA*BA
               BETI= 1.0d0/(BETA+BETA)
               ALBE= ALFX+BETX
               XP = EXP(-ALBE)
               DUM= COEFI(IG-IPMIN+1)*COEFJ(JG-JPMIN+1)
               DO 330 K=KF,KL
                  XXI= 1.0d0/SQRT(X12+ZLP(K))
                  ALEF= ALFA*XXI
                  BEIT= BETA*XXI
                  IF(ALEF*BEIT>ABLIM) THEN
                     XKA= ALFA+ALFA
                     XKB= BETA+BETA
                     A1 = ALEF+BEIT
                     A2 = ALEF-BEIT
                     C  = X12+ZLP(K)
                     XPLS= EXP(-ALBE+(ALEF+BEIT)**2)
                     XMNS= EXP(-ALBE+(ALEF-BEIT)**2)
                     GAM1= XPLS*0.25D+00
                     GAM2= XMNS*0.25D+00
                  ELSE
                     XPLS= EXP(-ALBE+ALEF*ALEF)
                     XMNS= EXP(-ALBE+BEIT*BEIT)
                  END IF
                  NLPK = NLP(K)
                  CLPK = CLP(K)*DUM
                  CALL FJECP(FJPQ,ALFI,BETI,NLPK,NPNP,LTMAX,LEMX)
                  DO 310 IN=1,LEMX+1
  310             FP(IN)= FP(IN)+FJPQ(IN,IN,1)*CLPK
                  IF(NPNP<2) GO TO 330
                  NBEG= 1
                  N= 6+1
                  DO 320 IN=2,NPNP
                     NBEG= 1-NBEG
                     LBEG= NBEG
                     DO 320 IP=1,NPCPL
                        LBEG= 1-LBEG
                        DO 320 IQ=1+LBEG,NPBPL,2
                  FP(N)= FP(N)+FJPQ(IQ,IP,IN)*CLPK
  320             N= N+1
  330          CONTINUE
  340    CONTINUE
!
         L= LL-2
         L1MAX= MAX0(1,L+1)
         L2PL=(LL-1)*(LL-2)
         X = XCA
         Y = YCA
         Z = ZCA
         CALL ZFN(ZFNLMC,NPCPL-1,ZLM,LMF,LMX,LMY,LMZ)
         X = XBA
         Y = YBA
         Z = ZBA
         CALL ZFN(ZFNLMB,NPBPL-1,ZLM,LMF,LMX,LMY,LMZ)
!        FPISQ2= FPI*FPI*Q2
         MMAX= JAMAX
         NN= 1
         DO 440 I=IAMIN,IAMAX
            DUMI = FPISQ2
            IF(NORMI) DUMI= DUMI*ECPFAC(I)
            IF(IANDJ) MMAX= I
            DO 440 J=JAMIN,MMAX
               DUMJ= DUMI
               IF(NORMJ) DUMJ= DUMJ*ECPFAC(J)
               NXC = NX(I)+1
               NYC = NY(I)+1
               NZC = NZ(I)+1
               NXB = NX(J)+1
               NYB = NY(J)+1
               NZB = NZ(J)+1
               DO 410 M=1,2*NP1-1
                  DO 410 N=1,NPCPL
                     DO 410 K=1,NPBPL
  410          CPQ(K,N,M)= 0.0d0
               NLM=(L2PL-L)*NTLIM
               DO 420 M=L,-L,-1
                  DO 412 N=1,NP1
                     DO 412 K=1,NP1PL
                        CKLC(K,N)= 0.0d0
                        CKLB(K,N)= 0.0d0
  412             CONTINUE
                  JF = JFST2(NLM+I)
                  JL = JFST2(NLM+I+1)-1
                  DO 414 JJ=JF,JL
                     LA = LBECP2(1,JJ)+1
                     KA = LBECP2(2,JJ)+1
                     MU = LBECP2(3,JJ)
                     KX = LBECP2(4,JJ)
                     KY = LBECP2(5,JJ)
                     KZ = LBECP2(6,JJ)
                     LU = LA*(LA-1)-MU+1
                     CKLC(LA,KA)= CKLC(LA,KA)+DCOEF2(JJ)*ZFNLMC(LU)*    &
                                  ACX(NXC-KX)*ACY(NYC-KY)*ACZ(NZC-KZ)
  414             CONTINUE
                  JF = JFST2(NLM+J)
                  JL = JFST2(NLM+J+1)-1
                  DO 416 JJ=JF,JL
                     LA = LBECP2(1,JJ)+1
                     KA = LBECP2(2,JJ)+1
                     MU = LBECP2(3,JJ)
                     KX = LBECP2(4,JJ)
                     KY = LBECP2(5,JJ)
                     KZ = LBECP2(6,JJ)
                     LU = LA*(LA-1)-MU+1
                     CKLB(LA,KA)= CKLB(LA,KA)+DCOEF2(JJ)*ZFNLMB(LU)*    &
                                  ABX(NXB-KX)*ABY(NYB-KY)*ABZ(NZB-KZ)
  416             CONTINUE
                  DO 418 K=1,NP1
                     DO 418 N=1,NPCPL
                        CT= CKLC(N,K)
                        DO 418 K2=K,K+NP1-1
                           DO 418 N2=1,NPBPL
  418             CPQ(N2,N,K2)= CPQ(N2,N,K2)+CT*CKLB(N2,K2+1-K)
  420          NLM= NLM+NTLIM
               SUM= 0.0d0
               DO 422 N=1,L1MAX
  422          SUM= SUM+FP(N)*CPQ(N,N,1)
               IF(NPNP<2) GO TO 425
               JBEG= 1
               N= 6+1
               DO 424 KK=2,NPNP
                  JBEG= 1-JBEG
                  LBEG= JBEG
                  DO 424 JJ=1,NPCPL
                     LBEG= 1-LBEG
                     DO 424 M=LBEG+1,NPBPL,2
               IF(ABS(CPQ(M,JJ,KK))>TOL) SUM= SUM+FP(N)*CPQ(M,JJ,KK)
  424          N= N+1
  425          CONTINUE
               G(NN)= G(NN)+DUMJ*SUM
  440    NN= NN+1
  500 CONTINUE
!
      RETURN
      END

! ECPPWR
      SUBROUTINE ECPPWR(ISWTCH,NI,NJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /ECP1  / X01,CAX,CAY,CAZ,CA,XCA,YCA,ZCA,                   &
                      X02,BAX,BAY,BAZ,BA,XBA,YBA,ZBA,                   &
                      PHASE,DAX,DAY,DAZ,DA,XDA,YDA,ZDA,XINT,KCNTR
      COMMON /ECP3  / ACX(7),ACY(7),ACZ(7),ABX(7),ABY(7),ABZ(7)
      ACX(1)= 1.0d0
      ACY(1)= 1.0d0
      ACZ(1)= 1.0d0
      IF(ISWTCH<1) THEN
         IF(NI>1) THEN
            ACX(2)=-CAX
            ACY(2)=-CAY
            ACZ(2)=-CAZ
            IF(NI>2) THEN
               DO 110 II=3,NI
                  ACX(II)=-CAX*ACX(II-1)
                  ACY(II)=-CAY*ACY(II-1)
                  ACZ(II)=-CAZ*ACZ(II-1)
  110          CONTINUE
            END IF
         END IF
      END IF
      ABX(1)= 1.0d0
      ABY(1)= 1.0d0
      ABZ(1)= 1.0d0
      IF(ISWTCH>(-1)) THEN
         IF(NJ>1) THEN
            ABX(2)=-BAX
            ABY(2)=-BAY
            ABZ(2)=-BAZ
            IF(NJ>2) THEN
               DO 120 II=3,NJ
                  ABX(II)=-BAX*ABX(II-1)
                  ABY(II)=-BAY*ABY(II-1)
                  ABZ(II)=-BAZ*ABZ(II-1)
  120          CONTINUE
            END IF
         END IF
      END IF
      RETURN
      END

! ERRF
      DOUBLE PRECISION FUNCTION ERRF(Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /ERRFCM/ C(142),IFRST(20),ILAST(20),H
!
      X= Y
      S= 1.0d0
      IF(X>=0.0d0) GO TO 105
      X=-Y
      S=-1.0d0
  105 IF(X>=4.86D+00) GO TO 115
      XN = X/H
      NX = INT(XN)+1
      I1 = IFRST(NX)
      IL = ILAST(NX)
      T= C(IL)
      DO 110 K=IL-1,I1,-1
  110 T= C(K)+T*X
      ERRF  = T*S
      GO TO 999
  115 ERRF  = S
  999 CONTINUE
      RETURN
      END

! FA
      DOUBLE PRECISION FUNCTION FA(N,ZETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GAMMO(9),GAMME(9)
      DATA GAMMO/0.5D+00,0.5D+00,1.0D+00,3.0D+00,12.0D+00,60.0D+00,     &
                 360.0D+00,2520.0D+00,20160.0D+00/                       
      DATA GAMME/0.5D+00,0.25D+00,0.375D+00,0.9375D+00,3.28125D+00,     &
                 14.765625D+00,81.2109375D+00,527.87109375D+00,         &
                 3959.033203125D+00/
      DATA SQRPI/1.772453850905D+00/

      K= N/2
      IF(MOD(N,2)/=0) THEN
         T= GAMMO(K+1)
      ELSE
         T= GAMME(K+1)*SQRPI*SQRT(ZETA)
      END IF
      FA= T/ZETA**(K+1)
      RETURN
      END

! FIECP
      SUBROUTINE FIECP(FIP,ALFI,NLP,NPNP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FIP(*)
      COMMON /FICMN / ALF,XI,XP0,XP1
      DIMENSION FIT(19,2)
      PARAMETER (ZP5=0.5D+00)
      NEFF=NLP+NPNP
      IF(MOD(NEFF,2)/=0) THEN
         NEMX= NEFF+1
         NOMX= NEFF
      ELSE
         NEMX= NEFF
         NOMX= NEFF+1
      END IF
      IF(NEFF==0) NOMX=-1
      YI= XI*ZP5
      A1= XI*ALF
      A2= XI*YI
      IF(NEMX<0) GO TO 115
      FIT(1,1)= FSI0(0)*YI
      IF(NEMX==0) GO TO 115
      FIT(2,2)= FSI1(1)*YI*YI
      IF(NEMX==1) GO TO 115
      T1= A2
      DO 110 N=2,NEMX,2
         T2= T1-A2
         FIT(N+1,1)= A1*FIT(N  ,2)+T1*FIT(N-1,1)
         FIT(N+2,2)= A1*FIT(N+1,1)+T2*FIT(N  ,2)
  110 T1= T1+A2+A2
  115 CONTINUE
      IF(NOMX<0) GO TO 125
      FIT(1,2)= FSI1(0)*YI
      IF(NOMX==0) GO TO 125
      FIT(2,1)= FSI0(1)*YI*YI
      IF(NOMX==1) GO TO 125
      T1= A2+A2
      DO 120 N=2,NOMX,2
         T2= T1-A2-A2-A2
         FIT(N+1,2)= A1*FIT(N  ,1)+T2*FIT(N-1,2)
         FIT(N+2,1)= A1*FIT(N+1,2)+T1*FIT(N  ,1)
  120 T1= T1+A2+A2
  125 CONTINUE
      LMAX= 1
      NK=0
      NM1K=0
      DO 150 N=1,NPNP+1
!        NK=(N*(N-1))/2
         DO 130 L=1,LMAX
  130    FIP(NK+L)= FIT(N+NLP,L)
         LMAX= 2
         IF(N<=2) GO TO 145
!        NM1K=((N-1)*(N-2))/2
         T01 = ALFI+ALFI+ALFI
         DO 140 L=3,N
            FIP(NK+L)= FIP(NK-2+L)-T01*FIP(NM1K-1+L)
  140    T01 = T01 +ALFI+ALFI
  145    NM1K=NK
  150 NK=NK+N
      RETURN
      END
      
! FJFORM
      SUBROUTINE FJFORM(RHO,SGA,SGB,TAU,NMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RHO(15),SGA(15),SGB(15),TAU(15)
      COMMON /FJCMN / ALEF,BET,XI,XPLS,XMNS,XP
      COMMON /FJNEW / XKA,XKB,GAM1,GAM2,A1,A2,C
      PARAMETER (ZP5=0.5D+00)
      DIMENSION FACTI(13),B1(15),B2(15),C1(15),C2(15)
      DATA FACTI/1.0D+00,4.0D+00,1.8D+01,9.6D+01,6.0D+02,4.32D+03,      &
                 3.528D+04,3.2256D+05,3.26592D+06,3.6288D+07,           &
                 4.390848D+08,5.7480192D+09,8.09512704D+10/
      DATA SQPI/1.772453850905D+00/
!
      SQPIDC= SQPI*XI
      B1(1) = SQPIDC*GAM1
      C1(1) = SQPIDC*GAM1*ERRF(A1)
      B2(1) = SQPIDC*GAM2
      C2(1) = SQPIDC*GAM2*ERRF(A2)
      TWSQPI= SQPI+SQPI
      B1(2) = TWSQPI*GAM1*DAWERF(A1)
      C1(2) = TWSQPI*GAM1*DAWF(A1)
      B2(2) = TWSQPI*GAM2*DAWERF(A2)
      C2(2) = TWSQPI*GAM2*DAWF(A2)
      TWC= C+C
      XKAPKB= XKA+XKB
      XKAMKB= XKA-XKB
      XPP= XP
      XPM= XP
      FAP= 0.0d0
      DEN= 1.0d0
      DO 110 I=3,NMAX
         XPP= XPP*XKAPKB
         XPM= XPM*XKAMKB
!        FAP=(1+((-1)**I))*0.25D+00
!        FAM=(1-((-1)**I))*0.25D+00
         FAM= ZP5-FAP
         T01= FAP/FACTI(I-2)
         T02= FAM/FACTI(I-2)
         T03= 1.0d0/DEN
         B1(I)= XPP*T01-(TWC*B1(I-2)-XKAPKB*C1(I-1))*T03
         C1(I)= XPP*T02-(TWC*C1(I-2)-XKAPKB*B1(I-1))*T03
         B2(I)= XPM*T01-(TWC*B2(I-2)-XKAMKB*C2(I-1))*T03
         C2(I)= XPM*T02-(TWC*C2(I-2)-XKAMKB*B2(I-1))*T03
      FAP= FAM
  110 DEN= DEN+1.0d0
      DO 120 I=1,NMAX
         RHO(I)= B1(I)-B2(I)
         SGA(I)= C1(I)+C2(I)
         SGB(I)= C1(I)-C2(I)
         TAU(I)= B1(I)+B2(I)
  120 CONTINUE
      RETURN
      END
      
! FJECP
      SUBROUTINE FJECP(FJPQ,ALFI,BETI,NLP,NPNP,LMAX,LEMX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FJPQ(11,11,11)
      COMMON /FJCMN / ALEF,BET,XI,XPLS,XMNS,XP
      COMMON /FJNEW / XKA,XKB,GAM1,GAM2,A1,A2,C
      DIMENSION RHO(15),SGA(15),SGB(15),TAU(15),FJT(19,2,2)
      PARAMETER (ZP5=0.5D+00, ABLIM=1.0D-01)
      PARAMETER (ZP1=0.1D+00, FIV=5.0D+00, SIX=6.0D+00, TEN=10.0D+00)
!
      JEFF= NLP+NPNP-1
      NEMX= JEFF
      NOMX= JEFF
      IF(JEFF==0) NOMX=-1
      LOMX= LEMX
      IF(ALEF*BET>ABLIM) THEN
         XKAKB = XKA*XKB
         ONEDAB= 1.0d0/XKAKB
         ONDAB2= ONEDAB*ONEDAB
!        IF(2*LEMX+3>15) THEN
         IF(LEMX>6) THEN
            WRITE(6,*) 'LEMX TOO BIG, LEMX= ',LEMX
            CALL ABRT
         END IF
         DO 110 II=1,15
            RHO(II)= 0.0d0
            SGA(II)= 0.0d0
            SGB(II)= 0.0d0
            TAU(II)= 0.0d0
  110    CONTINUE
         CALL FJFORM(RHO,SGA,SGB,TAU,2*LEMX+3)
         FJT(1,1,1)= RHO(3)*ONEDAB
         FJT(2,1,1)= RHO(2)*ONEDAB
         FJT(1,2,1)=(SGB(3)-RHO(4)*ONEDAB*XKB)*ONEDAB
         FJT(2,2,1)=(SGB(2)-RHO(3)*ONEDAB*XKB)*ONEDAB
         FJT(1,1,2)=(SGA(3)-RHO(4)*ONEDAB*XKA)*ONEDAB
         FJT(2,1,2)=(SGA(2)-RHO(3)*ONEDAB*XKA)*ONEDAB
         FJT(1,2,2)=(RHO(5)-XKA*SGB(4)-XKB*SGA(4)+XKAKB*TAU(3))*ONDAB2
         FJT(2,2,2)=(RHO(4)-XKA*SGB(3)-XKB*SGA(3)+XKAKB*TAU(2))*ONDAB2
      ELSE
         CALL SITABL(LEMX,LOMX)
         FJT(1,1,1)= FJPS(0,0,0)
         FJT(2,1,1)= FJPS(1,0,0)
         FJT(1,2,1)= FJPS(0,1,0)
         FJT(2,2,1)= FJPS(1,1,0)
         FJT(1,1,2)= FJPS(0,0,1)
         FJT(2,1,2)= FJPS(1,0,1)
         FJT(1,2,2)= FJPS(0,1,1)
         FJT(2,2,2)= FJPS(1,1,1)
      END IF
      XA= XI*ALEF
      XB= XI*BET
      XX= XI*XI*ZP5
      IF(NEMX<0) GO TO 125
      IF(NEMX==0) GO TO 125
      IF(NEMX==1) GO TO 125
      T2= 0.0d0
      DO 120 N=2,NEMX,2
         T1= T2+XX
         T3= T2-XX-XX-XX
         FJT(N+1,1,1)= XA*FJT(N  ,2,1)+XB*FJT(N  ,1,2)+T1*FJT(N-1,1,1)
         FJT(N+1,2,2)= XA*FJT(N  ,1,2)+XB*FJT(N  ,2,1)+T3*FJT(N-1,2,2)
         FJT(N+2,2,1)= XA*FJT(N+1,1,1)+XB*FJT(N+1,2,2)+T2*FJT(N  ,2,1)
         FJT(N+2,1,2)= XA*FJT(N+1,2,2)+XB*FJT(N+1,1,1)+T2*FJT(N  ,1,2)
  120 T2= T2+XX+XX
  125 CONTINUE
      IF(NOMX<0) GO TO 135
      IF(NOMX==0) GO TO 135
      IF(NOMX==1) GO TO 135
      T2=-XX
      DO 130 N=2,NOMX,2
         T1= T2-XX
         T3= T2+XX+XX+XX
         FJT(N+1,2,1)= XA*FJT(N  ,1,1)+XB*FJT(N  ,2,2)+T2*FJT(N-1,2,1)
         FJT(N+1,1,2)= XA*FJT(N  ,2,2)+XB*FJT(N  ,1,1)+T2*FJT(N-1,1,2)
         FJT(N+2,1,1)= XA*FJT(N+1,2,1)+XB*FJT(N+1,1,2)+T3*FJT(N  ,1,1)
         FJT(N+2,2,2)= XA*FJT(N+1,1,2)+XB*FJT(N+1,2,1)+T1*FJT(N  ,2,2)
  130 T2= T2+XX+XX
  135 CONTINUE
      IF(2<=LEMX .and. LEMX<=6) THEN
         IF(ALEF*BET>ABLIM) THEN
            IF(NLP>2) THEN
               WRITE(6,*) 'NLP TOO BIG!, ',NLP,' LEMX = ',LEMX
               CALL ABRT
            END IF
         END IF
      ELSEIF(LEMX>6) THEN
         WRITE(6,*) 'ERROR LEMX = ',LEMX,' IS TOO LARGE IN FJECP'
         CALL ABRT
      END IF
      N= NLP
      DO 150 IP= 1,LEMX+1
         IF(IP<3) THEN
            FJPQ(IP,IP,1)= FJT(N+1,IP,IP)
            GO TO 150
         ELSEIF(ALEF*BET<=ABLIM) THEN
            FJPQ(IP,IP,1)= FJPS(N,IP-1,IP-1)
            GO TO 150
         END IF
         IF(IP==3) THEN
            XKA2 = XKA*XKA
            XKB2 = XKB*XKB
            XKAKB2= XKAKB*XKAKB
            DU4= ONDAB2
            QA2= 3.0D+00*XKA2
            QB2= 3.0D+00*XKB2
            T02= QA2+QB2
            T00= 9.0D+00
      DU0 = T00*RHO( 7-N)+T02*RHO( 5-N)+XKAKB2*RHO( 3-N)
      DU1 = T00*SGA( 6-N)+QA2*SGA( 4-N)
      DU2 = T00*SGB( 6-N)+QB2*SGB( 4-N)
      DU3 = T00*TAU( 5-N)
         ELSEIF(IP==4) THEN
            QA2= 1.5D+01*XKA2
            QB2= 1.5D+01*XKB2
            T02= QA2+QB2
            XRE= SIX*XKAKB2
            T00= 2.25D+02
      DU0 = T00*RHO( 9-N)+(T02*RHO( 7-N)+RHO( 5-N)*XRE)*SIX
      DU1 = T00*SGA( 8-N)+(T02+FIV* QA2)*SGA( 6-N)+SGA( 4-N)*XRE
      DU2 = T00*SGB( 8-N)+(T02+FIV* QB2)*SGB( 6-N)+SGB( 4-N)*XRE
      DU3 = T00*TAU( 7-N)+ T02*TAU( 5-N)+TAU( 3-N)*XKAKB2
         ELSEIF(IP==5) THEN
            XKA4 = XKA2*XKA2
            XKB4 = XKB2*XKB2
            QA2= 4.5D+01*XKA2
            QB2= 4.5D+01*XKB2
            T02= QA2+QB2
            T03= 1.05D+03*(XKA2+XKB2)
            QA4= 1.05D+02*XKA4
            QB4= 1.05D+02*XKB4
            T04= QA4+QB4
            XRE= TEN*XKAKB2
            T00= 1.1025D+04
      DU0 = T00*RHO(11-N)+1.05D+02*T02*RHO( 9-N)+T04*RHO( 7-N)+         &
           (2.025D+03*RHO( 7-N)+T02*RHO( 5-N)+RHO( 3-N)*XKAKB2)*XKAKB2   
      DU1 = T00*SGA(10-N)+(T03+3.675D+03*XKA2)*SGA( 8-N)+               &
            QA4*SGA( 6-N)+(4.5D+01*SGA( 6-N)+XKA2*SGA( 4-N))*XRE         
      DU2 = T00*SGB(10-N)+(T03+3.675D+03*XKB2)*SGB( 8-N)+               &
            QB4*SGB( 6-N)+(4.5D+01*SGB( 6-N)+XKB2*SGB( 4-N))*XRE
      DU3 = T00*TAU( 9-N)+T03*TAU( 7-N)+TEN*TAU( 5-N)*XRE
         ELSEIF(IP==6) THEN
            QA2= 4.2D+02*XKA2
            QB2= 4.2D+02*XKB2
            T02= QA2+QB2
            T03= 9.9225D+04*(XKA2+XKB2)
            QA4= 9.45D+02*XKA4
            QB4= 9.45D+02*XKB4
            T04= QA4+QB4
            XRE= 1.5D+01*XKAKB2
            T00= 8.93025D+05
      DU0 = T00*RHO(13-N)+4.0D+00*T03*RHO(11-N)+1.5D+01*T04*RHO( 9-N)+  &
           (1.176D+04*RHO( 9-N)+T02*RHO( 7-N)+RHO( 5-N)*XRE)*XRE         
      DU1 = T00*SGA(12-N)+(T03+7.0875D+02*QA2)*SGA(10-N)+               &
            (T04+1.4D+01*QA4)*SGA( 8-N)+                                &
           (4.41D+04*SGA( 8-N)+(T02+2.75D+00*QA2)*SGA( 6-N)+            &
            SGA( 4-N)*XRE)*XKAKB2                                        
      DU2 = T00*SGB(12-N)+(T03+7.0875D+02*QB2)*SGB(10-N)+               &
            (T04+1.4D+01*QB4)*SGB( 8-N)+                                &
            (4.41D+04*SGB( 8-N)+(T02+2.75D+00*QB2)*SGB( 6-N)+           &
            SGB( 4-N)*XRE)*XKAKB2                                        
      DU3 = T00*TAU(11-N)+T03*TAU( 9-N)+T04*TAU( 7-N)+                  &
           (1.1025D+04*TAU( 7-N)+2.5D-01*T02*TAU( 5-N)+                 &
            TAU( 3-N)*XKAKB2)*XKAKB2
         ELSEIF(IP==7) THEN
            XKA6 = XKA4*XKA2
            XKB6 = XKB4*XKB2
            QA2= 3.5D+00*XKA2
            QB2= 3.5D+00*XKB2
            T02= QA2+QB2
            T03= 6.496875D+01*T02
            QA4= 3.5D+00*XKA4
            QB4= 3.5D+00*XKB4
            T04= 3.75D-01*(QA4+QB4)
            T05= 7.7D+00*T04
            QA6= 4.8125D-02*XKA6
            QB6= 4.8125D-02*XKB6
            T06= QA6+QB6
            XRE= ZP1*XKAKB2/SIX
            T00= 5.00259375D+02
      DU0 = T00*RHO(15-N)+T03*RHO(13-N)+T05*RHO(11-N)+T06*RHO( 9-N)+    &
            (6.2015625D+03*RHO(11-N)+                                   &
            7.875D+01*T02*RHO( 9-N)+T04*RHO( 7-N)+                      &
            (7.35D+02*RHO( 7-N)+T02*RHO( 5-N)+RHO( 3-N)*XRE)*XRE)*XRE    
      DU1 = T00*SGA(14-N)+(T03-1.66753125D+02*XKB2)*SGA(12-N)+          &
            (T05-2.59875D+00*QB4)*SGA(10-N)+QA6*SGA( 8-N)+              &
           (1.65375D+04*SGA(10-N)+2.1D+02*(T02-6.25D-01*QB2)*SGA( 8-N)+ &
            QA4*SGA( 6-N)+(7.35D+02*SGA( 6-N)+QA2*SGA( 4-N))*XRE        &
           )*XRE*ZP1                                                     
      DU2 = T00*SGB(14-N)+(T03-1.66753125D+02*XKA2)*SGB(12-N)+          &
            (T05-2.59875D+00*QA4)*SGB(10-N)+QB6*SGB( 8-N)+              &
           (1.65375D+04*SGB(10-N)+2.1D+02*(T02-6.25D-01*QA2)*SGB( 8-N)+ &
            QB4*SGB( 6-N)+(7.35D+02*SGB( 6-N)+QB2*SGB( 4-N))*XRE        &
           )*XRE*ZP1                                                     
      DU3 = T00*TAU(13-N)+1.7325D+01*T02*TAU(11-N)+ZP1*T05*TAU( 9-N)+   &
           (2.1D+02*TAU( 9-N)+T02*TAU( 7-N)+3.5D+00*TAU( 5-N)*XRE       &
           )*XRE*2.1D+00
      DU4 = DU4*2.16D+05
         END IF
         DU4= DU4*ONEDAB
         FJPQ(IP,IP,1)=(DU0-DU1*XKB-DU2*XKA+DU3*XKAKB)*DU4
  150 CONTINUE
      IF(NPNP<=1) GO TO 999
      DO 220 IK=2,NPNP
         DO 160 IP=1,2
            DO 160 IQ=1,2
  160    FJPQ(IQ,IP,IK)= FJT(NLP+IK,IP,IQ)
         IF(LMAX<=2) GO TO 220
         T01= ALFI+ALFI+ALFI
         DO 170 IP=3,LMAX
            FJPQ( 1,IP,IK)= FJPQ( 1,IP-2,IK)-T01*FJPQ( 1,IP-1,IK-1)
            FJPQ( 2,IP,IK)= FJPQ( 2,IP-2,IK)-T01*FJPQ( 2,IP-1,IK-1)
  170    T01= T01 +ALFI+ALFI
         T02= BETI+BETI+BETI
         DO 180 IQ=3,LMAX
            FJPQ(IQ, 1,IK)= FJPQ(IQ-2, 1,IK)-T02*FJPQ(IQ-1, 1,IK-1)
            FJPQ(IQ, 2,IK)= FJPQ(IQ-2, 2,IK)-T02*FJPQ(IQ-1, 2,IK-1)
  180    T02= T02 +BETI+BETI

         T01= ALFI+ALFI+ALFI
         DO 210 IP=3,LMAX
            T02= BETI+BETI+BETI
            IF(IP==3) GO TO 195
            DO 190 IQ=3,IP-1
               FJPQ(IQ,IP,IK)= FJPQ(IQ,IP-2,IK)-T01*FJPQ(IQ,IP-1,IK-1)
  190       T02= T02 +BETI+BETI
  195       CONTINUE
            DO 200 IQ=IP,LMAX
               FJPQ(IQ,IP,IK)= FJPQ(IQ-2,IP,IK)-T02*FJPQ(IQ-1,IP,IK-1)
  200       T02= T02 +BETI+BETI
  210    T01= T01 +ALFI+ALFI
  220 CONTINUE

  999 CONTINUE
      RETURN
      END
      
! FJPS
      DOUBLE PRECISION FUNCTION FJPS(N,LALF,LBET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FJCMN / ALF,BET,XI,XPLS,XMNS,XP
      COMMON /FSICMN/ SI(30,7)
      PARAMETER(LIMAB=10)
      DIMENSION DFCTRL(7)
      DATA DFCTRL/1.0D+00,3.0D+00,15.0D+00,105.0D+00,945.0D+00,         &
                  1.0395D+04,1.35135D+05/
      IF((ALF-BET)>0.0D+00) THEN
         L0= LBET+1
         X = BET
      ELSE
         L0= LALF+1
         X = ALF
      END IF
      T = X**(L0-1)/DFCTRL(L0)
      X = X*X
      L1= LALF+LBET+2-L0
      L2= L0+N
      L3= L0+L0-1
      SUM = T*SI(L2,L1)
      DO 110 K=1,LIMAB
         T= T*X/(2*K*(K+K+L3))
  110 SUM = SUM+T*SI(K+K+L2,L1)
      FJPS= SUM*(0.5D+00*XI)**(N+1)
      RETURN
      END

! FSI
      DOUBLE PRECISION FUNCTION FSI(N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FICMN / ALFA,XI,XP0,XP1
      SAVE DAWFS,ERRFS
      PARAMETER (ZP5=0.5D0, OP5=1.5D0)
      PARAMETER (TP5=2.5D0, THR=3.0D0, FIV=5.0D0)
      PARAMETER (SIX=6.0D0, SP5=7.5D0, FIF=15.0D0, ALIM=0.317D0)
      DATA SQPI/1.772453850905D0/
      DATA ERRFS/0.0D0/
!     -------------
      ENTRY FSI0(N)
!     -------------
      IF(ALFA>ALIM) GO TO 010
      FSI0= FSIPS(N,0)
      GO TO 030
  010 NP1 = N+1
      IF(NP1==2) GO TO 020
      DAWFS=DAWF(ALFA)
      FSI0= SQPI*XP1*DAWFS/ALFA
      GO TO 030
  020 FSI0= ERRFS/ALFA
  030 FSI = FSI0
      RETURN
!     -------------
      ENTRY FSI1(N)
!     -------------
      IF(ALFA>ALIM) GO TO 100
      FSI1= FSIPS(N,1)
      GO TO 140
  100 NP1 = N+1
      ALF2= ALFA*ALFA
      TLF2= ALF2+ALF2
      IF(NP1-2) 110,120,130
  110 ERRFS=SQPI*XP1*ERRF(ALFA)
      FSI1=(ZP5*ERRFS-XP0*ALFA)/ALF2
      GO TO 140
  120 FSI1= SQPI*XP1*(ALFA-DAWFS)/ALF2
      GO TO 140
  130 FSI1=(2.0d0*XP0*ALFA+(TLF2-1.0d0)*ERRFS)/ALF2
  140 FSI = FSI1
      RETURN
!     -------------
      ENTRY FSI2(N)
!     -------------
      IF(ALFA>ALIM) GO TO 200
      FSI2= FSIPS(N,2)
      GO TO 250
  200 NP1 = N+1
      ALF2= ALFA*ALFA
      ALF3= ALFA*ALF2
      TLF2= ALF2+ALF2
      GO TO(210,220,230,240),NP1
  210 FSI2= ZP5*SQPI*XP1*(OP5*ALFA-(ALF2+OP5)*DAWFS)/ALF3
      GO TO 250
  220 FSI2=(THR*XP0*ALFA+(ALF2-OP5)*ERRFS)/ALF3
      GO TO 250
  230 FSI2= SQPI*XP1*(ALFA*(TLF2-THR)+THR*DAWFS)/ALF3
      GO TO 250
  240 FSI2=(2.0d0*XP0*ALFA*(TLF2-THR)+(TLF2*(TLF2-2.0)+THR)*ERRFS)/ALF3
  250 FSI = FSI2
      RETURN
!     -------------
      ENTRY FSI3(N)
!     -------------
      IF(ALFA>ALIM) GO TO 300
      FSI3= FSIPS(N,3)
      GO TO 360
  300 NP1 = N+1
      ALF2= ALFA*ALFA
      ALF3= ALFA*ALF2
      ALF4= ALFA*ALF3
      TLF2= ALF2+ALF2
      GO TO(310,320,330,340,350),NP1
  310 FSI3=(2.0d0*XP0*ALFA*(TLF2+SP5)/SIX+ZP5*(ALF2-TP5)*ERRFS)/ALF4
      GO TO 360
  320 FSI3= ZP5*SQPI*XP1*(ALFA*(TLF2-SP5)+OP5*(TLF2+FIV)*DAWFS)/ALF4
      GO TO 360
  330 FSI3=(2.0d0*XP0*ALFA*(ALF2-SP5)+(TLF2*(ALF2-THR)+SP5)*ERRFS)/ALF4
      GO TO 360
  340 FSI3= SQPI*XP1*(ALFA*(TLF2*(TLF2-FIV)+FIF)-FIF*DAWFS)/ALF4
      GO TO 360
  350 FSI3=(2.0d0*XP0*ALFA*(TLF2*(TLF2-1.0d0-THR)+FIF)+                 &
            (TLF2*(TLF2*(TLF2-THR)+SIX+THR)-FIF)*ERRFS)/ALF4
  360 FSI = FSI3
      RETURN
      END
      
! FSIPS
      DOUBLE PRECISION FUNCTION FSIPS(N,L)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FICMN / ALFA,XI,XP0,XP1
      PARAMETER (LIMA=10)
      DIMENSION FCTRL(10),DFCTRL(10)
      DATA  FCTRL/1.0d0,1.0d0,2.0D+00,6.0D+00,24.0D+00,120.0D+00,       &
                  720.0D+00,5040.0D+00,40320.0D+00,362880.0D+00/         
      DATA DFCTRL/1.0d0,1.0d0,3.0D+00,15.0D+00,105.0D+00,945.0D+00,     &
                  10395.0D+00,135135.0D+00,2027025.0D+00,34459425.0D+00/
      DATA SQPI/1.772453850905D+00/
!
      NL= N+L
      LAMBDA=NL/2
      IF(MOD(NL,2)/=0) THEN
         T  =  FCTRL(LAMBDA+1)/DFCTRL(L+2)
         FAC= 2.0d0**NL
      ELSE
         T  = DFCTRL(LAMBDA+1)/DFCTRL(L+2)
         FAC=(2.0d0**LAMBDA)*SQPI
      ENDIF
      L2= L+L+1
      X = ALFA*ALFA
      SUM= T
      DO 110 K=1,LIMA
         T =(T*X*(K+K+NL-1))/(K*(K+K+L2))
  110 SUM= SUM+T
      FSIPS = FAC*XP0*(ALFA**L)*SUM
      RETURN
      END

! RECUR
      SUBROUTINE RECUR(NMIN,NMAX,LMAX,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FSICMN/ SI(30,7)
!
      TWX= X+X
      T01=(NMIN-2.0d0)*2.0d0
      IF(LMAX<2) THEN
         DO 110 N=NMIN,NMAX,2
            T02= T01+2.0d0
            SI(N+1,1)= T02*SI(N-1,1)+TWX*SI(N  ,2)
            SI(N+2,2)= T01*SI(N  ,2)+TWX*SI(N+1,1)
  110    T01= T02+2.0d0
      ELSE IF(LMAX==3) THEN
         DO 120 N=NMIN,NMAX,2
            T02= T01+2.0d0
            SI(N+1,1)= T02*SI(N-1,1)+TWX*SI(N  ,2)
            SI(N+2,2)= T01*SI(N  ,2)+TWX*SI(N+1,1)
            SI(N+3,3)= T01*SI(N+1,3)+TWX*SI(N+2,2)
  120    T01= T02+2.0d0
      ELSE
         DO 130 N=NMIN,NMAX,2
            T02= T01+2.0d0
            SI(N+1,1)= T02*SI(N-1,1)+TWX*SI(N  ,2)
            SI(N+2,2)= T01*SI(N  ,2)+TWX*SI(N+1,1)
            SI(N+3,3)= T01*SI(N+1,3)+TWX*SI(N+2,2)
            SI(N+4,4)= T01*SI(N+2,4)+TWX*SI(N+3,3)
            SI(N+5,5)= T01*SI(N+3,5)+TWX*SI(N+4,4)
            SI(N+6,6)= T01*SI(N+4,6)+TWX*SI(N+5,5)
            SI(N+7,7)= T01*SI(N+5,7)+TWX*SI(N+6,6)
  130    T01= T02+2.0d0
      END IF
      RETURN
      END
      
! SITABL
      SUBROUTINE SITABL(LEMAX,LOMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FICMN / ALFI,XI,XP0,XP1
      COMMON /FJCMN / ALFJ,BETJ,XJ,XPLS,XMNS,XP
      COMMON /FSICMN/ SI(30,7)
!
      IF((ALFJ-BETJ)<=0.0D+00) THEN
         ALFI= BETJ
         XP1 = XMNS
      ELSE
         ALFI= ALFJ
         XP1 = XPLS
      END IF
      X   = ALFI
      XI  = XJ
      XP0 = XP
      SI(1,1)= FSI0(0)
      SI(2,2)= FSI1(1)
      IF(LEMAX<=1) GO TO 110
      SI(3,3)= FSI2(2)
      IF(LEMAX==2) GO TO 110
      SI(4,4)= FSI3(3)
      IF(LEMAX==3) GO TO 110
      SI(5,5)= FSIPS(4,4)
      IF(LEMAX==4) GO TO 110
      SI(6,6)= FSIPS(5,5)
      IF(LEMAX==5) GO TO 110
      SI(7,7)= FSIPS(6,6)
  110 CONTINUE
      CALL RECUR(2,22,LEMAX,X)
      SI(1,2)= FSI1(0)
      SI(2,1)= FSI0(1)
      SI(3,2)= FSI1(2)
      IF(LOMAX<=1) GO TO 120
      SI(4,3)= FSI2(3)
      IF(LOMAX==2) GO TO 120
      SI(5,4)= FSI3(4)
      IF(LOMAX==3) GO TO 120
      SI(6,5)= FSIPS(4,5)
      IF(LOMAX==4) GO TO 120
      SI(7,6)= FSIPS(5,6)
      IF(LOMAX==5) GO TO 120
      SI(8,7)= FSIPS(6,7)
  120 CONTINUE
      CALL RECUR(3,21,LOMAX,X)
      RETURN
      END
      
! ZFN
      SUBROUTINE ZFN(ZFNLM,LMAX,ZLM,LMF,LMX,LMY,LMZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION ZLM(*)
      DIMENSION ZFNLM(121),LMF(*),LMX(*),LMY(*),LMZ(*)
      COMMON /ZFNCM / X,Y,Z
!
      IF(LMAX<=10) THEN
         DO 130 L=0,LMAX
            ID= L*(L+1)+L+1
            DO 120 M=-L,L
               IMN= LMF(ID)
               IMX= LMF(ID+1)-1
               SUM= 0.0d0
               DO 110 I=IMN,IMX
                  DUMMY = ZLM(I)
                  IF(LMX(I)>0) DUMMY= DUMMY*(X**LMX(I))
                  IF(LMY(I)>0) DUMMY= DUMMY*(Y**LMY(I))
                  IF(LMZ(I)>0) DUMMY= DUMMY*(Z**LMZ(I))
  110          SUM= SUM+DUMMY
               ZFNLM(ID)= SUM
  120       ID= ID-1
  130    CONTINUE
      ELSE
         WRITE(6,9010) LMAX
         CALL ABRT
      END IF
      RETURN
 9010 FORMAT(' ERROR: ZFN MAX LAMDA=10, YOU REQUESTED LAMDA = ',I3)
      END

! ECPCBA
      SUBROUTINE ECPCBA(CANDA,AANDB,ICAB,IPOW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CANDA,AANDB
      COMMON /ECP1  / X01,CAX,CAY,CAZ,CA,XCA,YCA,ZCA,                   &
                      X02,BAX,BAY,BAZ,BA,XBA,YBA,ZBA,                   &
                      PHASE,DAX,DAY,DAZ,DA,XDA,YDA,ZDA,XINT,KCNTR
      LOGICAL CANDB
      COMMON /ECP4  / P12(3,2),R12,ACO(3),CANDB
      COMMON /ZFNCM / X,Y,Z

      IF(CANDA) THEN
         IF(AANDB) THEN
            ICAB = 1
         ELSE
            CAX= 0.0d0
            CAY= 0.0d0
            CAZ= 0.0d0
            CA = 0.0d0
            BAX= P12(1,2)-ACO(1)
            BAY= P12(2,2)-ACO(2)
            BAZ= P12(3,2)-ACO(3)
            BA = SQRT(BAX*BAX+BAY*BAY+BAZ*BAZ)
            X  = BAX/BA
            Y  = BAY/BA
            Z  = BAZ/BA
            ICAB= 2
            IPOW= 1
         END IF
      ELSE
         IF(AANDB) THEN
            CAX= P12(1,1)-ACO(1)
            CAY= P12(2,1)-ACO(2)
            CAZ= P12(3,1)-ACO(3)
            CA = SQRT(CAX*CAX+CAY*CAY+CAZ*CAZ)
            X  = CAX/CA
            Y  = CAY/CA
            Z  = CAZ/CA
            BAX= 0.0d0
            BAY= 0.0d0
            BAZ= 0.0d0
            BA = 0.0d0
            ICAB= 3
            IPOW=-1
         ELSE
            CAX= P12(1,1)-ACO(1)
            CAY= P12(2,1)-ACO(2)
            CAZ= P12(3,1)-ACO(3)
            CA = SQRT(CAX*CAX+CAY*CAY+CAZ*CAZ)
            XCA= CAX/CA
            YCA= CAY/CA
            ZCA= CAZ/CA
            BAX= P12(1,2)-ACO(1)
            BAY= P12(2,2)-ACO(2)
            BAZ= P12(3,2)-ACO(3)
            BA = SQRT(BAX*BAX+BAY*BAY+BAZ*BAZ)
            XBA= BAX/BA
            YBA= BAY/BA
            ZBA= BAZ/BA
            ICAB= 4
            IPOW= 0
         END IF
      END IF
      RETURN
      END

!----------------------------------------------------------------------!
!                                                                      !
!       2020 Use Libreta open source library for ERI calculation       !                                                                      !                                                                      !
!  implemented by Juan Felipe Huan Lew Yee and Jorge Martin del Campo  ! 
!                                                                      !
!----------------------------------------------------------------------!
      
! HSandTlib                                           
      SUBROUTINE HSandTlib(H,S,TKIN,NBFT,KATOM,KLOC,KMIN,KMAX,NSHELL,   &
                           ZAN,C)        
      USE ISO_C_BINDING
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION,DIMENSION(NBFT) :: H,S,TKIN
      TYPE(C_PTR),DIMENSION(600)::BASLIB
      COMMON/LIBRETA/BASLIB
      LOGICAL     LINEAR
      COMMON/ZMAT/LINEAR      
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB 
      LOGICAL IIANDJJ
      INTEGER,DIMENSION(NSHELL) :: KATOM,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
      DOUBLE PRECISION,DIMENSION(3,NAT) :: C
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: Z,SBLK,TBLK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: VBLK,ZBLK
      ALLOCATE(Z(NBFT))
      ALLOCATE(SBLK(784),TBLK(784),VBLK(784),ZBLK(784))
!-----------------------------------------------------------------------
!                      H, S & TKIN integrals
!-----------------------------------------------------------------------
      ZBLK = 0.0D+00
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     I SHELL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO III = 1,NSHELL
       I = KATOM(III)
       MINI = KMIN(III)
       MAXI = KMAX(III)
       LOCI = KLOC(III)-MINI
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      J SHELL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO JJ = 1,III
        J = KATOM(JJ)
        MINJ = KMIN(JJ)
        MAXJ = KMAX(JJ)
        LOCJ = KLOC(JJ)-MINJ
        IIANDJJ = III == JJ
!       Call to libreta
!lib        CALL sval(BASLIB,III,JJ,SBLK)
!lib        CALL tval(BASLIB,III,JJ,TBLK)
!lib        CALL vval(BASLIB,III,JJ,NAT,C(1,1:NAT),C(2,1:NAT),C(3,1:NAT),   &
!lib                  ZAN(1:NAT),VBLK)
!       avoiding warnings
        ZAN(1) = ZAN(1)
        C(1,1) = C(1,1)
!lib
!
        JMAX = MAXJ
        NN = 0
        DO I = MINI,MAXI
         LI = LOCI+I
         IN = (LI*(LI-1))/2
         IF (IIANDJJ) JMAX = I
         DO J = MINJ,JMAX
          LJ = LOCJ+J
          JN = LJ+IN
          NN = NN+1
          H(JN) =  TBLK(NN) + VBLK(NN)
          S(JN) =  SBLK(NN)
          TKIN(JN) =  TBLK(NN)
          IF(LINEAR) Z(JN) = ZBLK(NN)
         END DO
        END DO
       END DO
      END DO
!-----------------------------------------------------------------------
      DEALLOCATE(Z,SBLK,TBLK,VBLK,ZBLK)
      RETURN
      END

! PRCALClib                                           
      SUBROUTINE PRCALClib(XVAL,WINT,NVAL,L2,KATOM,KLOC,KMIN,KMAX,NSHELL)               
      USE ISO_C_BINDING              
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      TYPE(C_PTR),DIMENSION(600)::BASLIB
      COMMON/LIBRETA/BASLIB
      LOGICAL IIANDJJ
      INTEGER,DIMENSION(NSHELL) :: KATOM,KLOC,KMIN,KMAX
      DIMENSION XVAL(NVAL*L2),WINT(*)                                   
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO II=1,NSHELL
       I    = KATOM(II)
       MINI = KMIN(II)
       MAXI = KMAX(II)
       LOCI = KLOC(II) - MINI
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO JJ=1,II
        J    = KATOM(JJ)
        MINJ = KMIN(JJ)
        MAXJ = KMAX(JJ)
        LOCJ = KLOC(JJ) - MINJ
        IIandJJ = II==JJ
        IJ = 0
        JMAX = MAXJ
        DO I=MINI,MAXI
         IF (IIandJJ) JMAX = I
         DO J=MINJ,JMAX
          IJ = IJ+1
         END DO
        END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!          Integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!lib        CALL DIPVAL(BASLIB,II,JJ,WINT)
        JMAX = MAXJ
        DO K=1,NVAL
         NL2 = (K-1)*L2
         NN  = (K-1)*IJ
         DO I=MINI,MAXI
          LI = LOCI + I
          IN = LI*(LI-1)/2 + NL2
          IF (IIandJJ) JMAX = I
          DO J=MINJ,JMAX
           LJ = LOCJ + J
           JN = LJ + IN
           NN = NN+1
           XVAL(JN) = WINT(NN)
          END DO
         END DO
        END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

! ExchangeIntlib                                           
      SUBROUTINE ExchangeIntlib(XINTS,GLIBRETA,NSH2,MAXG,               &
                                KTYPE,KMIN,KMAX,NSHELL)   
      USE ISO_C_BINDING              
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER,DIMENSION(NSHELL) :: KTYPE,KMIN,KMAX
      TYPE(C_PTR),DIMENSION(600)::BASLIB
      COMMON/LIBRETA/BASLIB
      COMMON/SHLEXC/NORGSH(3),NORGSP(3),IEXCH,NGTH(4)  
      COMMON/SHLNOS1/QQ4,IJKL                                                                              
      DIMENSION XINTS(NSH2),GLIBRETA(MAXG)
      LOGICAL ISHandJSH
      INTEGER::ORII,ORIJ
!-----------------------------------------------------------------------
      CALL BASCHK(LMAXIMA,KTYPE,NSHELL)
      NANGM =  4                                          
      IF(LMAXIMA==2) NANGM =  6                                          
      IF(LMAXIMA==3) NANGM = 10                                          
      IF(LMAXIMA==4) NANGM = 15                                          
      IF(LMAXIMA==5) NANGM = 21                                          
      IF(LMAXIMA==6) NANGM = 28                                          
      NGTH(4) = 1                                                       
      NGTH(3) = NGTH(4) * NANGM                                         
      NGTH(2) = NGTH(3) * NANGM                                         
      NGTH(1) = NGTH(2) * NANGM                                         
      DO I=1,3                                                       
       NORGSH(I) = 0                                               
       NORGSP(I) = 0                                               
      ENDDO                                                          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IEXCH = 1                                                         
      QQ4   = 1.0d0                                                       
      NINTEG  = 0  
!- - - - - - - - - - - - - - - - - - - - - - - -
      IJIJ = 0
      DO ISH = 1,NSHELL
       DO JSH = 1,ISH
        IJIJ = IJIJ+1
        VMAX = 0.0D+00
        MINI = KMIN(ISH)
        MINJ = KMIN(JSH)
        MAXI = KMAX(ISH)
        JMAX = KMAX(JSH)
        ISHandJSH=ISH==JSH
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ORII = KMAX(ISH) - KMIN(ISH) + 1
        ORIJ = KMAX(JSH) - KMIN(JSH) + 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        JMAX = ORIJ
        DO I=1,ORII
         IF(ISHandJSH) JMAX = I
          DO J=1,JMAX
!lib          CALL erisval(BASLIB,ISH,JSH,ISH,JSH,GLIBRETA)
          IJKL_INDEX = (I-1)*ORIJ*ORII*ORIJ                             &
                     + (J-1)*ORII*ORIJ+(I-1)*ORIJ + (J-1)+1
          VAL = GLIBRETA(IJKL_INDEX)
          IF(VAL>0.0D+00)NINTEG = NINTEG + 1
          IF(VAL>VMAX)VMAX=VAL
         END DO
        END DO
        XINTS(IJIJ) = SQRT(VMAX)
       END DO
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

! QOUTlib
      SUBROUTINE QOUTlib(BUFP,IX,BUFP2,IX2,NINTEGtm,NINTMX,GLIBRETA,    &
                         IDONTW,ISH,JSH,KSH,LSH,KLOC,KMIN,KMAX,NSHELL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IANDJ,KANDL,SAME
      COMMON/MISC/IANDJ,KANDL,SAME
      COMMON /RESTAR/NREC,IST,JST,KST,LST
      COMMON/SHLT/SHLTOL,CUTOFF,ICOUNT
      COMMON/SHLNOS/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,MINI,MINJ,MINK, &
                    MINL,MAXI,MAXJ,MAXK,MAXL,NIJ,IJ,KL
      INTEGER :: ORII,ORIJ,ORIK,ORIL,ISH,JSH,KSH,LSH
      INTEGER,DIMENSION(NSHELL) :: KLOC,KMIN,KMAX
      INTEGER,DIMENSION(NINTEGtm) :: IX2                                 
      DIMENSION BUFP(NINTMX),IX(NINTMX),GLIBRETA(*)
      DOUBLE PRECISION,DIMENSION(NINTEGtm) :: BUFP2                      
      SAVE IJN,KLN
!-----------------------------------------------------------------------
!     Pack 4-indices into 1 word.
!     Write Label & integral on Unit 1 if DONTW = .False.
!-----------------------------------------------------------------------
      SAME  = ISH == KSH .and. JSH == LSH
      IANDJ = ISH == JSH
      KANDL = KSH == LSH
      MINI = KMIN(ISH)
      MINJ = KMIN(JSH)
      MINK = KMIN(KSH)
      MINL = KMIN(LSH)
      MAXI = KMAX(ISH)
      MAXJ = KMAX(JSH)
      MAXK = KMAX(KSH)
      MAXL = KMAX(LSH)
      LOCI = KLOC(ISH)-MINI
      LOCJ = KLOC(JSH)-MINJ
      LOCK = KLOC(KSH)-MINK
      LOCL = KLOC(LSH)-MINL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ORII = MAXI - MINI + 1
      ORIJ = MAXJ - MINJ + 1
      ORIK = MAXK - MINK + 1
      ORIL = MAXL - MINL + 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IJN = 0
      JMAX = ORIJ
      DO I = 1,ORII
       I_INDEX = (I-1)*ORIJ*ORIK*ORIL
       IF (IANDJ) JMAX = I
       DO 1 J = 1,JMAX
        IJ_INDEX = (J-1)*ORIK*ORIL + I_INDEX
        IJN = IJN+1
        LMAX = ORIL
        KLN = 0
        DO K = 1,ORIK
         IJK_INDEX = (K-1)*ORIL + IJ_INDEX
         IF (KANDL) LMAX = K
         DO L = 1,LMAX
          KLN = KLN+1
          IF(SAME.and.KLN>IJN)GO TO 1
          IJKL_INDEX = (L-1) + IJK_INDEX + 1
          VAL = GLIBRETA( IJKL_INDEX )
          IF(ABS(VAL)>=CUTOFF)THEN
           I1 = LOCI+I+MINI-1
           I2 = LOCJ+J+MINJ-1
           I3 = LOCK+K+MINK-1
           I4 = LOCL+L+MINL-1

           IF (I1 >= I2) GO TO 100
           N = I1
           I1 = I2
           I2 = N
  100      IF (I3 >= I4) GO TO 120
           N = I3
           I3 = I4
           I4 = N
  120      IF (I1-I3) 140,160,180
  140      N = I1
           I1 = I3
           I3 = N
           N = I2
           I2 = I4
           I4 = N
           GO TO 180
  160      IF (I2 < I4) GO TO 140
  180      CONTINUE
!
           IF(I1 == I2) VAL = VAL*0.5D0
           IF(I3 == I4) VAL = VAL*0.5D0
           IF(I1 == I3 .and. I2 == I4) VAL = VAL*0.5D0
!
           NPACK = ICOUNT
           IPACK = I1
           JPACK = I2
           KPACK = I3
           LPACK = I4
           LABEL = ISHFT( IPACK, 48 ) + ISHFT( JPACK, 32 ) +            &
                   ISHFT( KPACK, 16 ) + LPACK
           IX(NPACK) = LABEL
           BUFP(ICOUNT) = VAL
           ICOUNT = ICOUNT+1
           IF(ICOUNT > 0) THEN
            JCOUNT = ICOUNT
            IF(JCOUNT > NINTMX) THEN
             NXInteg = NINTMX
             IF(IDONTW==1)THEN
              IJBUFi = (NREC-1)*NINTMX
              do ibuf=1,NINTMX
               IX2  (IJBUFi+ibuf) = IX(ibuf)
               BUFP2(IJBUFi+ibuf) = BUFP(ibuf)
              end do
             ELSE
              WRITE(1)NXInteg,IX,BUFP
             END IF
             ICOUNT = 1
             NREC = NREC+1
            END IF
           END IF
          END IF
         END DO
        END DO
    1  CONTINUE
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

! UPDCOOSHELL
      SUBROUTINE UPDCOOSHELL(NSHELL,KATOM,Cxyz,NAT)     
      USE ISO_C_BINDING
      INTEGER :: NSHELL,NAT 
      INTEGER,DIMENSION(NSHELL) :: KATOM
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz 
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      LOGICAL SMCD
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      TYPE(C_PTR),DIMENSION(600)::BASLIB,AUXLIB
      COMMON/LIBRETA/BASLIB
      COMMON/LIBRETAaux/AUXLIB
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
      INTEGER :: I,IATOM
!-----------------------------------------------------------------------      
!     Update coordinates of shells
!-----------------------------------------------------------------------
      if(IERITYP==1 .or. IERITYP==3)then      ! ERITYP = FULL or MIX
       DO I=1,NSHELL
        IATOM = KATOM(I)
!       avoiding warning
        Cxyz(1:3,IATOM) = Cxyz(1:3,IATOM)
!lib        CALL updatebasis(BASLIB,I,Cxyz(1:3,IATOM))
       END DO
      end if 
      if(IERITYP==2 .or. IERITYP==3)then      ! ERITYP = RI or MIX
       DO I=1,NSHELLaux
        IATOM = KATOMaux(I)
!lib        CALL updatebasis(AUXLIB,I,Cxyz(1:3,IATOM))
       END DO      
      end if       
!-----------------------------------------------------------------------      
      RETURN
      END
            
!----------------------------------------------------------------------!
!                                                                      !
!       2020  RI Approximation for ERIs implemented by                 !
!             Juan Felipe Huan Lew Yee and Jorge Martin del Campo      !
!                                                                      !
!             ( J. Chem. Phys. 154, 064102, 2021 )                     !
!                                                                      !
!----------------------------------------------------------------------!

! AUXREAD
      SUBROUTINE AUXREAD(NAT,NSHELLaux,NUMaux,NATmax,ANAM)
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
        NUMaux = NUMaux + (L+1)*(L+2)/2
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
         CALL NORMALIZE_AUXILIAR(EXaux(K),COEFICIENT,L)
         Caux(K) = C1*COEFICIENT
        END DO
        ! Compute contracted normalization constant
        FACL = 0.0D0
        DO IG = K1,K2
         DO JG = K1,IG
          EE = EXaux(IG)+EXaux(JG)
          FAC = EE*SQRT(EE)
          IF(L==0) DUM = Caux(IG)*Caux(JG)/FAC
          IF(L==1) DUM = 0.5D0*Caux(IG)*Caux(JG)/(EE*FAC)
          IF(L==2) DUM = PT75  *Caux(IG)*Caux(JG)/(EE*EE*FAC)
          IF(L==3) DUM = PT187 *Caux(IG)*Caux(JG)/(EE**3*FAC)
          IF(L==4) DUM = PT6562*Caux(IG)*Caux(JG)/(EE**4*FAC)
          IF(L==5) DUM = PT2953*Caux(IG)*Caux(JG)/(EE**5*FAC)
          IF(L==6) DUM = PT1624*Caux(IG)*Caux(JG)/(EE**6*FAC)
          IF(IG /= JG) THEN
           DUM = DUM+DUM
          END IF
          FACL = FACL+DUM
         END DO
        END DO
        IF(FACL < TM10) THEN
         FACL = 0.0D0
        ELSE
         FACL = 1.0D0/SQRT(FACL*PI32)
        END IF
        DO K = K1,K2
         Caux(K) = Caux(K) * FACL
        END DO
        GOTO 2
       END IF
      END DO
      CLOSE(UNIT=50)
!     Open Basis Set File ( Unit = 50 )                                 
      OPEN(50,FILE=BASIS_FILE,STATUS='UNKNOWN',                        &
              FORM='FORMATTED',ACCESS='SEQUENTIAL')

      RETURN
      END

! AUXGEN
      SUBROUTINE AUXGEN(NAT,NPRIMI,ITYP,IMIN,IMAX,NSHELLaux,NUMaux,     &
                        IGEN,ISTAR,EX,ZAN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                         
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
      COMMON/EXCaux/EXaux(2000),Caux(2000)
!
      INTEGER::IGEN,N
      INTEGER,DIMENSION(NAT) :: IMIN,IMAX
      INTEGER,DIMENSION(NPRIMI) :: ITYP
      INTEGER,DIMENSION(NAT) :: LMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX
      DOUBLE PRECISION,DIMENSION(NAT) :: EXMAX,EXMIN,ZAN
      DOUBLE PRECISION :: COEFICIENT

      DO IAT=1,NAT
       MINI = IMIN(IAT)
       MAXI = IMAX(IAT)
       LMAX(IAT) = MAXVAL(ITYP(MINI:MAXI)) - 1
       EXMAX(IAT) = MAXVAL(EX(MINI:MAXI))
       EXMIN(IAT) = MINVAL(EX(MINI:MAXI))
      END DO

      NUMaux = 0
      NSHELLaux = 0
      DO IAT=1,NAT

        lbmax = LMAX(IAT)
        zbmax = EXMAX(IAT)
        zbmin = EXMIN(IAT)

        r = 6.0 - IGEN

        nblock = 2
        if(ISTAR==1) nblock = nblock + 1

!       H and He
        if (ZAN(IAT).LE.2) then
          r = r - 0.5*IGEN + 2.0
          nblock = nblock - 1
        end if

!       How many z?
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
              KSTARTaux(NSHELLaux) = NSHELLaux
              KATOMaux(NSHELLaux) = IAT
              KTYPEaux(NSHELLaux) = L+1
              KNGaux(NSHELLaux) = 1
              KLOCaux(NSHELLaux) = NUMaux+1
              NUMaux = NUMaux + (L+1)*(L+2)/2
              EXaux(NSHELLaux) = TMPEXP
              CALL NORMALIZE_AUXILIAR(TMPEXP,COEFICIENT,L)
              Caux(NSHELLaux) = COEFICIENT
            END DO
            if (ifa==1) TMPEXP = TMPEXP*r/(r+0.5*IGEN)
          END DO
       END DO
      END DO
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

! AUXGENlib
      SUBROUTINE AUXGENlib(NAT,NPRIMI,ITYP,IMIN,IMAX,NSHELLaux,         &
                           NUMaux,EX,ZAN,Cxyz)
      USE ISO_C_BINDING
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                         
!
      INTEGER :: NPRIMI,NSHELLaux,NUMaux            
      INTEGER,DIMENSION(NPRIMI) :: ITYP      
      INTEGER,DIMENSION(NAT) :: IMIN,IMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX      
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz      
!
      TYPE(C_PTR),DIMENSION(600)::AUXLIB
      COMMON/LIBRETAaux/AUXLIB
      LOGICAL SMCD
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
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
          KLOCaux(NSHELLaux) = NUMaux+1
          NUMaux = NUMaux + (L+1)*(L+2)/2
          KTYPEaux(NSHELLaux) = L+1
          KATOMaux(NSHELLaux) = IAT
!         avoiding warning
          Cxyz(1:3,IAT) = Cxyz(1:3,IAT)
!lib          CALL CREATEBASISaux(AUXLIB(NSHELLaux),TMPEXP,1.0,             &
!lib                              Cxyz(1:3,IAT),L,1)
         end do
         if (ifa==1) TMPEXP = TMPEXP*r/(r+0.5*IGEN)
        end do
       end do
      end do
!-----------------------------------------------------------------------      
      RETURN
      END

! JandKaux                
      SUBROUTINE JandKaux(BUFP2,NINTEGtm,IDONTW,IPRINTOPT,NBF,          &
                          EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,KSTART,KATOM,  &
                          KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,Cxyz,NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/USELIBRETA/ILIBRETA
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      INTEGER :: NINTEGtm,IDONTW,IPRINTOPT,NBF,NSHELL,NAT
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NINTEGtm) :: BUFP2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: Gaux
      INTEGER :: CR, CM, timestarttwoE, timefinishtwoE
      DOUBLE PRECISION :: RATE
      LOGICAL SMCD
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
      CALL Debut(IDONTW,IPRINTOPT,KATOM,NSHELL,Cxyz)
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
       if(ILIBRETA==0)then
        CALL AuxERI(NINTEGtm,BUFP2,Gaux,MAXG,NBF,IPRINTOPT,             &
                    EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,KSTART,KATOM,        &
                    KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,Cxyz,NAT)
       else if(ILIBRETA==1)then
        CALL AuxERIlib(NINTEGtm,BUFP2,Gaux,MAXG,MAXORI,NBF,             &
                       KTYPE,KLOC,NSHELL,IPRINTOPT)
       end if
      ELSE IF(SMCD) THEN
       CALL AuxERIModChol(BUFP2,Gaux,MAXG,NBF,IPRINTOPT,                &
                          EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,KSTART,KATOM,  &
                          KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,Cxyz,NAT)             
      END IF
      DEALLOCATE(Gaux)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL SYSTEM_CLOCK(timefinishtwoE)
      DeltaTtwoE = (timefinishtwoE - timestarttwoE)/RATE
      IF(IPRINTOPT==1)                                                  &
       WRITE(6,'(1X,A22,F10.2)')'Time to do integrals =',DeltaTtwoE
!-----------------------------------------------------------------------
      RETURN
      END

! AuxERI                                            
      SUBROUTINE AuxERI(NINTEGtm,BUFP2,GHONDO,MAXG,NBF,IPRINTOPT,       &
                        EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,KSTART,KATOM,    &
                        KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,Cxyz,NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
      COMMON/RESTAR/NREC,IST,JST,KST,LST
      COMMON/SHLNOS1/QQ4,IJKL
      COMMON/SHLEXC/NORGSH(3),NORGSP(3),IEXCH,NGTH(4)
!
      INTEGER :: NINTEGtm,MAXG,NBF,IPRINTOPT,NPRIMI,NSHELL,NAT
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DIMENSION GHONDO(MAXG)
      DOUBLE PRECISION,DIMENSION(NINTEGtm) :: BUFP2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::GMAT
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AUX
!-----------------------------------------------------------------------
      ALLOCATE(AUX(49*900))
      BUFP2(:) = 0.0D0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Evaluate G = (P|Q) and G^{-1/2}
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(GMAT(NBFaux,NBFaux))
      CALL METRICmat(GMAT,GHONDO,MAXG,Cxyz,NAT)
      CALL PDPT_msqrt(GMAT,NBFaux,IPRINTOPT)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     2e- Integrals (S,P,D,F,G & L Shells)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL BASCHK(LMAXIMA,KTYPE,NSHELL)
      CALL BASCHK(LMAXIMAaux,KTYPEaux,NSHELLaux)
      NANGM = 0
      IF(LMAXIMA==0) NANGM =  1
      IF(LMAXIMA==1) NANGM =  3
      IF(LMAXIMA==2) NANGM =  6
      IF(LMAXIMA==3) NANGM = 10
      IF(LMAXIMA==4) NANGM = 15
      IF(LMAXIMA==5) NANGM = 21
      IF(LMAXIMA==6) NANGM = 28
      NANGMaux = 0
      IF(LMAXIMAaux==0) NANGMaux =  1
      IF(LMAXIMAaux==1) NANGMaux =  3
      IF(LMAXIMAaux==2) NANGMaux =  6
      IF(LMAXIMAaux==3) NANGMaux = 10
      IF(LMAXIMAaux==4) NANGMaux = 15
      IF(LMAXIMAaux==5) NANGMaux = 21
      IF(LMAXIMAaux==6) NANGMaux = 28
      NGTH(4) = 1
      NGTH(3) = NGTH(4) * 1
      NGTH(2) = NGTH(3) * NANGMaux
      NGTH(1) = NGTH(2) * NANGM
!
      NORGSH = 0
      NORGSP = 0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO II = IST,NSHELL
       J0 = JST
       DO JJ = J0,II
        JST = 1
        K0 = KST
        DO KK = K0,NSHELLaux
         KST = 1
!- - - - - - - - - - - - - - - - - - - - - - - -
!        (II,JJ//KK)                                        
!- - - - - - - - - - - - - - - - - - - - - - - -                            
         IEXCH = 1
         ISH = II
         JSH = JJ
         KSH = KK
         QQ4 = 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
!         Compute 2e- Integrals (mn|k)                  
!         Select integral code for ERI calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         LQSUM = KTYPE(ISH) + KTYPE(JSH) + KTYPEaux(KSH) - 3
         CALL SHELLS3C(1,ISH,JSH,KSH,.TRUE.,EX,CS,CP,CD,CF,CG,CH,CI,    &
                       NPRIMI,KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,    &
                       NSHELL,Cxyz,NAT)
         CALL IJPRIM(AUX)
         CALL SHELLS3C(2,ISH,JSH,KSH,.TRUE.,EX,CS,CP,CD,CF,CG,CH,CI,    &
                       NPRIMI,KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,    &
                       NSHELL,Cxyz,NAT)
         NORGH = NORGSH(IEXCH)
         CALL ZQOUT(GHONDO,MAXG)
         IF(LQSUM==0) THEN
          CALL S0000(GHONDO(1+NORGH),AUX)
         ELSE
         CALL ERISPDFGHIL(GHONDO(1+NORGH),AUX)
         END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!         Contract to B tensor
!- - - - - - - - - - - - - - - - - - - - - - - - - - - -
         CALL QOUT3C(BUFP2,GHONDO,MAXG,GMAT,NBF,NBFaux,                 &
                     KLOC,KMIN,KMAX,NSHELL)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - -
        END DO
       END DO
      END DO
      DEALLOCATE(AUX)
      DEALLOCATE(GMAT)
!-----------------------------------------------------------------------
      RETURN
      END

! AuxERIlib
      SUBROUTINE AuxERIlib(NINTEGtm,BUFP2,GLIBRETA,MAXG,MAXORI,NBF,     &
                           KTYPE,KLOC,NSHELL,IPRINTOPT)
      USE ISO_C_BINDING
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER :: NINTEGtm,MAXG,MAXORI,NBF,NSHELL,IPRINTOPT
      INTEGER :: ORII,ORIJ,ORIK      
      COMMON/ORI/ORII,ORIJ,ORIK            
      COMMON/RESTAR/NREC,IST,JST,KST,LST                 
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      COMMON/SHLEXC/NORGSH(3),NORGSP(3),IEXCH,NGTH(4)      
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
!
      TYPE(C_PTR),DIMENSION(600)::BASLIB,AUXLIB
      COMMON/LIBRETA/BASLIB
      COMMON/LIBRETAaux/AUXLIB
!      
      INTEGER,DIMENSION(NSHELL) :: KTYPE,KLOC      
      DOUBLE PRECISION :: GLIBRETA(MAXG),BUFP2(NINTEGtm)
      DOUBLE PRECISION :: GMAT(NBFaux,NBFaux)
      DOUBLE PRECISION :: ERIS3C(MAXORI*MAXORI*NBFaux)
!-----------------------------------------------------------------------
!lib      CALL metric_mat(NBFaux,GMAT,NSHELLaux,AUXLIB)
      CALL PDPT_msqrt(GMAT,NBFaux,IPRINTOPT)
!
      DO II = IST,NSHELL                            ! II Shell
       J0 = JST                                     ! JJ Shell
       DO JJ = J0,II
        JST = 1
        K0 = KST
        LI = KTYPE(II)
        LJ = KTYPE(JJ)
        ORII = LI*(LI+1)/2
        ORIJ = LJ*(LJ+1)/2
        ERIS3C(1:ORII*ORIJ*NBFaux) = 0.0d0
        DO KK = K0,NSHELLaux                        ! KK Shell
         KST = 1
!- - - - - - - - - - - - - - - - - - - - - - - -
!         (II,JJ//P)
!- - - - - - - - - - - - - - - - - - - - - - - -
          ISH = II
          JSH = JJ
          KSH = KK
          QQ4 = 1
          !IF(SKIPA .and. NPSYM) QQ4 = QQ4+1
          !IF(SKIPB .and. NPSYM) QQ4 = QQ4+1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!         Compute 2e- Integrals
!         Select integral code for ERI calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          LK = KTYPEaux(KSH)
          ORIK = LK*(LK+1)/2

!lib          CALL erisval3(BASLIB,AUXLIB,ISH,JSH,KSH,GLIBRETA)

          LOCK = KLOCaux(KSH)
          N = 0
          DO I=1,ORII
           DO J=1,ORIJ
            DO L=1,NBFaux
             N = N + 1
             DO K=1,ORIK
              ERIS3C(N) = ERIS3C(N) +                                   &
              GLIBRETA((I-1)*ORIJ*ORIK+(J-1)*ORIK+K)*GMAT(LOCK+K-1,L)
             END DO
            END DO
           END DO
          END DO
!- - - - - - - - - - - - - - - - - - - - - - - -
        END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Write Label & Integral on File 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
        CALL QOUTaux(BUFP2,ERIS3C,NBF,NBFaux,ISH,JSH,KLOC,NSHELL)
       END DO
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

! AuxERIModChol
      SUBROUTINE AuxERIModChol(BUFP2,GHONDO,MAXG,NBF,IPRINTOPT,         &
                               EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,KSTART,   &
                               KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,   &
                               Cxyz,NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
      COMMON/RESTAR/NREC,IST,JST,KST,LST
      COMMON/SHLNOS1/QQ4,IJKL
      COMMON/SHLEXC/NORGSH(3),NORGSP(3),IEXCH,NGTH(4)
!
      INTEGER :: MAXG,NBF,IPRINTOPT,NPRIMI,NSHELL,NAT
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DIMENSION GHONDO(MAXG)
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::BUFP
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::GMAT,L,D,P
      DOUBLE PRECISION,DIMENSION(NBF*(NBF+1)/2,NBFaux) :: BUFP2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AUX

      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::E
      INTEGER,ALLOCATABLE,DIMENSION(:)::IPIV
      INTEGER::INFO
!-----------------------------------------------------------------------
      ALLOCATE(AUX(49*900))
!-----------------------------------------------------------------------
!     Evaluate G = (P|Q), G = PLDL^TP^T, ModChol and get P, L, D^1/2
!-----------------------------------------------------------------------
      ALLOCATE(BUFP(NBF*(NBF+1)/2,NBFaux))
      ALLOCATE(GMAT(NBFaux,NBFaux))
      ALLOCATE(L(NBFaux,NBFaux),D(NBFaux,NBFaux),P(NBFaux,NBFaux))
      ALLOCATE(E(NBFaux),IPIV(NBFaux))

      CALL METRICmat(GMAT,GHONDO,MAXG,Cxyz,NAT)
      CALL LDLT(GMAT,IPIV,E,NBFaux)
      CALL MODCHOL( NBFaux, GMAT, IPIV, E, 1D-10,IPRINTOPT)
      CALL GET_PLD12(GMAT,IPIV,E,NBFaux,P,L,D)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     2e- Integrals (S,P,D,F,G & L Shells)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL BASCHK(LMAXIMA,KTYPE,NSHELL)
      CALL BASCHK(LMAXIMAaux,KTYPEaux,NSHELLaux)
      NANGM =  0
      IF(LMAXIMA==0) NANGM =  1
      IF(LMAXIMA==1) NANGM =  3
      IF(LMAXIMA==2) NANGM =  6
      IF(LMAXIMA==3) NANGM = 10
      IF(LMAXIMA==4) NANGM = 15
      IF(LMAXIMA==5) NANGM = 21
      IF(LMAXIMA==6) NANGM = 28
      NANGMaux =  0
      IF(LMAXIMAaux==0) NANGMaux =  1
      IF(LMAXIMAaux==1) NANGMaux =  3
      IF(LMAXIMAaux==2) NANGMaux =  6
      IF(LMAXIMAaux==3) NANGMaux = 10
      IF(LMAXIMAaux==4) NANGMaux = 15
      IF(LMAXIMAaux==5) NANGMaux = 21
      IF(LMAXIMAaux==6) NANGMaux = 28
      NGTH(4) = 1
      NGTH(3) = NGTH(4) * 1
      NGTH(2) = NGTH(3) * NANGMaux
      NGTH(1) = NGTH(2) * NANGM
!
      NORGSH = 0
      NORGSP = 0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO II = IST,NSHELL
       J0 = JST
       DO JJ = J0,II
        JST = 1
        K0 = KST
        DO KK = K0,NSHELLaux
         KST = 1
!- - - - - - - - - - - - - - - - - - - - - - - -
!        (II,JJ//KK)
!- - - - - - - - - - - - - - - - - - - - - - - -
         IEXCH = 1
         ISH = II
         JSH = JJ
         KSH = KK
         QQ4 = 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!         Compute 2e- Integrals (mn|k)
!         Select integral code for ERI calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         LQSUM = KTYPE(ISH) + KTYPE(JSH) + KTYPEaux(KSH) - 3
         CALL SHELLS3C(1,ISH,JSH,KSH,.TRUE.,EX,CS,CP,CD,CF,CG,CH,CI,    &
                       NPRIMI,KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,    &
                       NSHELL,Cxyz,NAT)
         CALL IJPRIM(AUX)
         CALL SHELLS3C(2,ISH,JSH,KSH,.TRUE.,EX,CS,CP,CD,CF,CG,CH,CI,    &
                       NPRIMI,KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,    &
                       NSHELL,Cxyz,NAT)
         NORGH = NORGSH(IEXCH)
         CALL ZQOUT(GHONDO,MAXG)
         IF(LQSUM==0) THEN
          CALL S0000(GHONDO(1+NORGH),AUX)
         ELSE
         CALL ERISPDFGHIL(GHONDO(1+NORGH),AUX)
         END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!         Store 3 center ERIs (mn|k)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         CALL QOUT3CModChol(BUFP,GHONDO,MAXG,NBF,NBFaux,               &
                            KLOC,KMIN,KMAX,NSHELL)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        END DO
       END DO
      END DO
      DEALLOCATE(AUX)
!-----------------------------------------------------------------------
!     Build b tensor, solve linear equation system LD^1/2 b = P^T(k|mn)
!-----------------------------------------------------------------------
      BUFP2 = MATMUL(BUFP,P)
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

! METRICmat                                            
      SUBROUTINE METRICmat(GMAT,GHONDO,MAXG,Cxyz,NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
      COMMON/SHLEXC/NORGSH(3),NORGSP(3),IEXCH,NGTH(4)
      COMMON/RESTAR/NREC,IST,JST,KST,LST   
      COMMON/SHLNOS1/QQ4,IJKL
!
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(MAXG) :: GHONDO
      DOUBLE PRECISION,DIMENSION(NBFaux,NBFaux) :: GMAT
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AUX
!-----------------------------------------------------------------------
!     2e- Integrals (S,P,D,F,G & L Shells)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(AUX(49*900))        

      CALL BASCHK(LMAXIMAaux,KTYPEaux,NSHELLaux)
      NANGMAUX =  0                                          
      IF(LMAXIMAAUX==0) NANGMAUX =  1                                   
      IF(LMAXIMAAUX==1) NANGMAUX =  3                                   
      IF(LMAXIMAAUX==2) NANGMAUX =  6                                   
      IF(LMAXIMAAUX==3) NANGMAUX = 10                                  
      IF(LMAXIMAAUX==4) NANGMAUX = 15                                   
      IF(LMAXIMAAUX==5) NANGMAUX = 21                                   
      IF(LMAXIMAAUX==6) NANGMAUX = 28                                   
      NGTH(4) = 1                                                       
      NGTH(3) = NGTH(4) * 1                                     
      NGTH(2) = NGTH(3) * NANGMAUX                                      
      NGTH(1) = NGTH(2) * 1                                      
      DO I=1,3                                                       
       NORGSH(I) = 0                                               
       NORGSP(I) = 0                                               
      ENDDO                                                          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO II = IST,NSHELLaux                        ! II Shell
       K0 = KST                                                        
       DO KK = K0,II                               ! KK Shell     
        KST = 1                                                       
!- - - - - - - - - - - - - - - - - - - - - - - -
!       (II//KK)                                        
!- - - - - - - - - - - - - - - - - - - - - - - -                            
        IEXCH = 1                                                    
        ISH = II                                                      
        KSH = KK                                                     
        QQ4 = 1                                                       
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
!       Compute 2e- Integrals                      
!       Select integral code for ERI calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        LQSUM = KTYPEaux(ISH) + KTYPEaux(KSH) - 2    
        CALL SHELLS2C(1,ISH,KSH,Cxyz,NAT)                           
        CALL SHELLS2C(2,ISH,KSH,Cxyz,NAT)                           
        CALL IJPRIM(AUX) 
        NORGH = NORGSH(IEXCH)                                         
        CALL ZQOUT(GHONDO,MAXG)                 
        IF(LQSUM==0) THEN 
         CALL S0000(GHONDO(1+NORGH),AUX)                        
        ELSE                                                          
         CALL ERISPDFGHIL(GHONDO(1+NORGH),AUX)                       
        END IF  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Transfer integrals to G matrix
!- - - - - - - - - - - - - - - - - - - - - - -
        CALL QOUT2C(GMAT,NBFaux,GHONDO,MAXG)
!- - - - - - - - - - - - - - - - - - - - - - -
       END DO
      END DO
      DEALLOCATE(AUX) 
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               

! SHELLS2C                                           
      SUBROUTINE SHELLS2C(NELEC,ISH,KSH,Cxyz,NAT)                   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL       NORM
      COMMON/NORMAL/NORM
      LOGICAL     IANDJ,KANDL,SAME
      COMMON/MISC/IANDJ,KANDL,SAME
      COMMON/ERIOUT/INU,JNU,KNU,LNU,NGTI,NGTJ,NGTK,NGTL                
      COMMON/INTDEX/IJX(784),IJY(784),IJZ(784),IK(784),                 &
                     KLX(784),KLY(784),KLZ(784)                         
      COMMON/INTDEX1/IJGT(784),KLGT(784)
      COMMON/ROOT/XX,U(13),W(13),NROOTS
      COMMON/SHLEXC/NORGSH(3),NORGSP(3),IEXCH,NGTH(4)
      COMMON/SHLINF/ GA(30),CSA(30),CPA(30),CDA(30),                    &
                    CFA(30),CGA(30),CHA(30),CIA(30),                    &
                     GB(30),CSB(30),CPB(30),CDB(30),                    &
                    CFB(30),CGB(30),CHB(30),CIB(30),                    &
                     GC(30),CSC(30),CPC(30),CDC(30),                    &
                    CFC(30),CGC(30),CHC(30),CIC(30),                    &
                     GD(30),CSD(30),CPD(30),CDD(30),                    &
                    CFD(30),CGD(30),CHD(30),CID(30),                    &
                    AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,                     &
                    DX,DY,DZ,RCD,NGA,NGB,NGC,NGD                                      
      COMMON/SHLNOS/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,                &
                    MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,            &
                    NIJ,IJ,KL                                          
      COMMON/SHLNOS1/QQ4,IJKL                                           
!      
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
      COMMON/EXCaux/EXaux(2000),Caux(2000)
!
      DIMENSION IX(84),IY(84),IZ(84),JX(84),JY(84),JZ(84),              &
                KX(84),KY(84),KZ(84),LX(84),LY(84),LZ(84)
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      INTEGER,DIMENSION(8) :: MINF,MAXF                                  
      DATA MINF / 1, 2,  5, 11, 21, 36, 57, 1/                           
      DATA MAXF / 1, 4, 10, 20, 35, 56, 84, 4/                           
      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,       &
                  3,   0,   0,   2,   2,   1,   0,   1,   0,   1,       &
                  4,   0,   0,   3,   3,   1,   0,   1,   0,   2,       &
                  2,   0,   2,   1,   1,                                &
                  5,   0,   0,   4,   4,   1,   0,   1,   0,   3,       &
                  3,   2,   0,   2,   0,   3,   1,   1,   2,   2,       &
                  1,                                                    &
                  6,   0,   0,   5,   5,   1,   0,   1,   0,   4,       &
                  4,   2,   0,   2,   0,   4,   1,   1,   3,   3,       &
                  0,   3,   3,   2,   1,   2,   1,   2/                  
      DATA KX /   0,   7,   0,   0,  14,   0,   0,   7,   7,   0,       &
                 21,   0,   0,  14,  14,   7,   0,   7,   0,   7,       &
                 28,   0,   0,  21,  21,   7,   0,   7,   0,  14,       &
                 14,   0,  14,   7,   7,                                &
                 35,   0,   0,  28,  28,   7,   0,   7,   0,  21,       &
                 21,  14,   0,  14,   0,  21,   7,   7,  14,  14,       &
                  7,                                                    &
                 42,   0,   0,  35,  35,   7,   0,   7,   0,  28,       &
                 28,  14,   0,  14,   0,  28,   7,   7,  21,  21,       &
                  0,  21,  21,  14,   7,  14,   7,  14/                  
      DATA JX /   0,  49,   0,   0,  98,   0,   0,  49,  49,   0,       &
                147,   0,   0,  98,  98,  49,   0,  49,   0,  49,       &
                196,   0,   0, 147, 147,  49,   0,  49,   0,  98,       &
                 98,   0,  98,  49,  49,                                &
                245,   0,   0, 196, 196,  49,   0,  49,   0, 147,       &
                147,  98,   0,  98,   0, 147,  49,  49,  98,  98,       &
                 49,                                                    &
                294,   0,   0, 245, 245,  49,   0,  49,   0, 196,       &
                196,  98,   0,  98,   0, 196,  49,  49, 147, 147,       &
                  0, 147, 147,  98,  49,  98,  49,  98/                  
      DATA IX /   1, 344,   1,   1, 687,   1,   1, 344, 344,   1,       &
               1030,   1,   1, 687, 687, 344,   1, 344,   1, 344,       &
               1373,   1,   1,1030,1030, 344,   1, 344,   1, 687,       &
                687,   1, 687, 344, 344,                                &
               1716,   1,   1,1373,1373, 344,   1, 344,   1,1030,       &
               1030, 687,   1, 687,   1,1030, 344, 344, 687, 687,       &
                344,                                                    &
               2059,   1,   1,1716,1716, 344,   1, 344,   1,1373,       &
               1373, 687,   1, 687,   1,1373, 344, 344,1030,1030,       &
                  1,1030,1030, 687, 344, 687, 344, 687/                  
      DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,       &
                  0,   3,   0,   1,   0,   2,   2,   0,   1,   1,       &
                  0,   4,   0,   1,   0,   3,   3,   0,   1,   2,       &
                  0,   2,   1,   2,   1,                                &
                  0,   5,   0,   1,   0,   4,   4,   0,   1,   2,       &
                  0,   3,   3,   0,   2,   1,   3,   1,   2,   1,       &
                  2,                                                    &
                  0,   6,   0,   1,   0,   5,   5,   0,   1,   2,       &
                  0,   4,   4,   0,   2,   1,   4,   1,   3,   0,       &
                  3,   2,   1,   3,   3,   1,   2,   2/                  
      DATA KY /   0,   0,   7,   0,   0,  14,   0,   7,   0,   7,       &
                  0,  21,   0,   7,   0,  14,  14,   0,   7,   7,       &
                  0,  28,   0,   7,   0,  21,  21,   0,   7,  14,       &
                  0,  14,   7,  14,   7,                                &
                  0,  35,   0,   7,   0,  28,  28,   0,   7,  14,       &
                  0,  21,  21,   0,  14,   7,  21,   7,  14,   7,       &
                 14,                                                    &
                  0,  42,   0,   7,   0,  35,  35,   0,   7,  14,       &
                  0,  28,  28,   0,  14,   7,  28,   7,  21,   0,       &
                 21,  14,   7,  21,  21,   7,  14,  14/                  
      DATA JY /   0,   0,  49,   0,   0,  98,   0,  49,   0,  49,       &
                  0, 147,   0,  49,   0,  98,  98,   0,  49,  49,       &
                  0, 196,   0,  49,   0, 147, 147,   0,  49,  98,       &
                  0,  98,  49,  98,  49,                                &
                  0, 245,   0,  49,   0, 196, 196,   0,  49,  98,       &
                  0, 147, 147,   0,  98,  49, 147,  49,  98,  49,       &
                 98,                                                    &
                  0, 294,   0,  49,   0, 245, 245,   0,  49,  98,       &
                  0, 196, 196,   0,  98,  49, 196,  49, 147,   0,       &
                147,  98,  49, 147, 147,  49,  98,  98/                  
      DATA IY /   1,   1, 344,   1,   1, 687,   1, 344,   1, 344,       &
                  1,1030,   1, 344,   1, 687, 687,   1, 344, 344,       &
                  1,1373,   1, 344,   1,1030,1030,   1, 344, 687,       &
                  1, 687, 344, 687, 344,                                &
                  1,1716,   1, 344,   1,1373,1373,   1, 344, 687,       &
                  1,1030,1030,   1, 687, 344,1030, 344, 687, 344,       &
                687,                                                    &
                  1,2059,   1, 344,   1,1716,1716,   1, 344, 687,       &
                  1,1373,1373,   1, 687, 344,1373, 344,1030,   1,       &
               1030, 687, 344,1030,1030, 344, 687, 687/                 
      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,       &
                  0,   0,   3,   0,   1,   0,   1,   2,   2,   1,       &
                  0,   0,   4,   0,   1,   0,   1,   3,   3,   0,       &
                  2,   2,   1,   1,   2,                                &
                  0,   0,   5,   0,   1,   0,   1,   4,   4,   0,       &
                  2,   0,   2,   3,   3,   1,   1,   3,   1,   2,       &
                  2,                                                    &
                  0,   0,   6,   0,   1,   0,   1,   5,   5,   0,       &
                  2,   0,   2,   4,   4,   1,   1,   4,   0,   3,       &
                  3,   1,   2,   1,   2,   3,   3,   2/                 
      DATA KZ /   0,   0,   0,   7,   0,   0,  14,   0,   7,   7,       &
                  0,   0,  21,   0,   7,   0,   7,  14,  14,   7,       &
                  0,   0,  28,   0,   7,   0,   7,  21,  21,   0,       &
                 14,  14,   7,   7,  14,                                &
                  0,   0,  35,   0,   7,   0,   7,  28,  28,   0,       &
                 14,   0,  14,  21,  21,   7,   7,  21,   7,  14,       &
                 14,                                                    &
                  0,   0,  42,   0,   7,   0,   7,  35,  35,   0,       &
                 14,   0,  14,  28,  28,   7,   7,  28,   0,  21,       &
                 21,   7,  14,   7,  14,  21,  21,  14/                 
      DATA JZ /   0,   0,   0,  49,   0,   0,  98,   0,  49,  49,       &
                  0,   0, 147,   0,  49,   0,  49,  98,  98,  49,       &
                  0,   0, 196,   0,  49,   0,  49, 147, 147,   0,       &
                 98,  98,  49,  49,  98,                                &
                  0,   0, 245,   0,  49,   0,  49, 196, 196,   0,       &
                 98,   0,  98, 147, 147,  49,  49, 147,  49,  98,       &
                 98,                                                    &
                  0,   0, 294,   0,  49,   0,  49, 245, 245,   0,       &
                 98,   0,  98, 196, 196,  49,  49, 196,   0, 147,       &
                147,  49,  98,  49,  98, 147, 147,  98/                 
      DATA IZ /   1,   1,   1, 344,   1,   1, 687,   1, 344, 344,       &
                  1,   1,1030,   1, 344,   1, 344, 687, 687, 344,       &
                  1,   1,1373,   1, 344,   1, 344,1030,1030,   1,       &
                687, 687, 344, 344, 687,                                &
                  1,   1,1716,   1, 344,   1, 344,1373,1373,   1,       &
                687,   1, 687,1030,1030, 344, 344,1030, 344, 687,       &
                687,                                                    &
                  1,   1,2059,   1, 344,   1, 344,1716,1716,   1,       &
                687,   1, 687,1373,1373, 344, 344,1373,   1,1030,       &
               1030, 344, 687, 344, 687,1030,1030, 687/                 
!-----------------------------------------------------------------------                                                                       
!     PREPARE SHELL INFORMATION/FOR HONDO INTEGRATION        
      NORM = .FALSE.  
      
      IF(NELEC==2) GO TO 200                                          
                                                                       
      IANDJ = .FALSE.                                         
      INU = ISH                                                      
      JNU = 1                                                      
      NGTI = NGTH(1)                                                 
      NGTJ = NGTH(2)                                                 
!                                                                       
!     ----- ISHELL                                                      
!                                                                       
      I = KATOMaux(INU)                                                 
      AX = Cxyz(1,I)                                                       
      AY = Cxyz(2,I)                                                       
      AZ = Cxyz(3,I)                                          
      I1 = KSTARTAUX(INU)   !
      I2 = I1+KNGAUX(INU)-1 !
      LIT = KTYPEaux(INU)                                              
      MINI = MINF(LIT)                                                  
      MAXI = MAXF(LIT)                                                  
      LOCI = KLOCaux(INU)-MINI                                          
      NGA = 0                                                           
      DO I = I1,I2                                                  
       NGA = NGA+1    
       GA(NGA) = EXaux(I)                                             
       CSA(NGA) = Caux(I)                                            
       CPA(NGA) = Caux(I)                                           
       CDA(NGA) = Caux(I)                                           
       CFA(NGA) = Caux(I)                                            
       CGA(NGA) = Caux(I)                                           
       CHA(NGA) = Caux(I)                                           
       CIA(NGA) = Caux(I)                                           
      END DO
!                                                                       
!     ----- JSHELL (Unity Shell)                                                     
!                                                                       
      J = 1                                                   
      BX = 0                                                  
      BY = 0                                                       
      BZ = 0                                                  
      J1 = 1                                              
      J2 = 1                                             
      LJT = 1                                                  
      MINJ = 1                                               
      MAXJ = 1                                               
      LOCJ = 1                                          
      NGB = 0                                                           
      DO J = J1,J2                                                  
       NGB = NGB+1                                                    
       GB(NGB) = 0                                               
       CSB(NGB) = 1                                           
       CPB(NGB) = 1                                          
       CDB(NGB) = 1                                           
       CFB(NGB) = 1                                           
       CGB(NGB) = 1                                         
       CHB(NGB) = 1                                         
       CIB(NGB) = 1                                           
      END DO                                                          
      RAB = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))       
!                                                                       
!     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS                
!                                                                       
      IJ = 0                                                            
      JMAX = MAXJ                                                       
      DO I = MINI,MAXI                                              
       NX = IX(I)                                                     
       NY = IY(I)                                                     
       NZ = IZ(I)                                                     
       IF (IANDJ) JMAX = I                                            
       DO J = MINJ,JMAX                                           
        IJ = IJ+1                                                   
        IJX(IJ) = NX+JX(J)                                          
        IJY(IJ) = NY+JY(J)                                          
        IJZ(IJ) = NZ+JZ(J)                                          
        IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1                    
       END DO                                                       
      END DO                                                          
      RETURN                                                            
!     ******                                                            
!                                                                       
!        K AND L SHELL                                                  
!                                                                       
  200 CONTINUE                                                          
      KANDL = .FALSE.                                        
      SAME = ISH == KSH                           
!                                                                       
!     ----- PERMUTE KSH AND LSH SHELLS, FOR THEIR TYPE                  
!                                                                       
      KNU = KSH                                                      
      LNU = 1                                                      
      NGTK = NGTH(3)                                                 
      NGTL = NGTH(4)                                                 
!                                                                       
!     ----- K SHELL                                                     
!                                                                       
      K = KATOMaux(KNU)                                                 
      CX = Cxyz(1,K)                                                       
      CY = Cxyz(2,K)                                                       
      CZ = Cxyz(3,K)                                                       
      K1 = KSTARTAUX(KNU)   !
      K2 = K1+KNGAUX(KNU)-1 !
      LKT = KTYPEaux(KNU) 
      MINK = MINF(LKT)                                                  
      MAXK = MAXF(LKT)                                                  
      LOCK = KLOCaux(KNU)-MINK                                          
      NGC = 0                                                           
      DO K = K1,K2                                                  
       NGC = NGC+1                                                    
       GC(NGC) = EXaux(K)                                             
       CSC(NGC) = Caux(K)                                          
       CPC(NGC) = Caux(K)                                       
       CDC(NGC) = Caux(K)                                       
       CFC(NGC) = Caux(K)                                     
       CGC(NGC) = Caux(K)                                         
       CHC(NGC) = Caux(K)                                          
       CIC(NGC) = Caux(K)                                          
      END DO                                                          
!                                                                       
!     ----- LSHELL (Unity Shell)                                                      
!                                                                       
      L = 1                                          
      DX = 0                                                    
      DY = 0                                                  
      DZ = 0                                                      
      L1 = 1                                                 
      L2 = 1                                          
      LLT = 1                                        
      MINL = 1                                            
      MAXL = 1                                            
      LOCL = 1                                       
      NGD = 0                                                           
      DO L = L1,L2                                                  
       NGD = NGD+1                                                    
       GD(NGD) = 0                                           
       CSD(NGD) = 1                                          
       CPD(NGD) = 1                                          
       CDD(NGD) = 1                                           
       CFD(NGD) = 1                                          
       CGD(NGD) = 1                                          
       CHD(NGD) = 1                                          
       CID(NGD) = 1                                           
      END DO                                                          
      NROOTS = (LIT+LJT+LKT+LLT-2)/2                                    
      RCD = ((CX-DX)*(CX-DX) + (CY-DY)*(CY-DY) + (CZ-DZ)*(CZ-DZ))       
!                                                                       
!     ----- PREPARE INDICES FOR PAIRS OF (K,L) FUNCTIONS                
!                                                                       
      KL = 0                                                            
      LMAX = MAXL                                                       
      DO K = MINK,MAXK                                              
       NX = KX(K)                                                     
       NY = KY(K)                                                     
       NZ = KZ(K)                                                     
       IF (KANDL) LMAX = K                                            
       DO L = MINL,LMAX                                           
        KL = KL+1                                                   
        KLX(KL) = NX+LX(L)                                          
        KLY(KL) = NY+LY(L)                                          
        KLZ(KL) = NZ+LZ(L)                                          
        KLGT(KL) = NGTK*(K-MINK)+NGTL*(L-MINL)                      
       END DO                                                       
      END DO                                                      
      MAX = KL                                                          
      DO 320 I = 1,IJ                                                   
      IF (SAME) MAX = I                                                 
  320 IK(I) = MAX                                                       
      IJKL = IJ*KL        
      IF (SAME) IJKL = IJ*(IJ+1)/2 
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               

! SHELLS3C                                           
      SUBROUTINE SHELLS3C(NELEC,ISH,JSH,KSH,FLIP,EX,CS,CP,CD,CF,CG,CH,  &
                          CI,NPRIMI,KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,   &
                          KMAX,NSHELL,Cxyz,NAT)                   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL       NORM
      COMMON/NORMAL/NORM
      LOGICAL     IANDJ,KANDL,SAME
      COMMON/MISC/IANDJ,KANDL,SAME
      COMMON/ERIOUT/INU,JNU,KNU,LNU,NGTI,NGTJ,NGTK,NGTL                
      COMMON/INTDEX/IJX(784),IJY(784),IJZ(784),IK(784),                 &
                    KLX(784),KLY(784),KLZ(784)                           
      COMMON/INTDEX1/IJGT(784),KLGT(784)                                 
      COMMON/ROOT/XX,U(13),W(13),NROOTS                                 
      COMMON/SHLEXC/NORGSH(3),NORGSP(3),IEXCH,NGTH(4)                    
      COMMON/SHLINF/ GA(30),CSA(30),CPA(30),CDA(30),                    &
                    CFA(30),CGA(30),CHA(30),CIA(30),                    &
                     GB(30),CSB(30),CPB(30),CDB(30),                    &
                    CFB(30),CGB(30),CHB(30),CIB(30),                    &
                     GC(30),CSC(30),CPC(30),CDC(30),                    &
                    CFC(30),CGC(30),CHC(30),CIC(30),                    &
                     GD(30),CSD(30),CPD(30),CDD(30),                    &
                    CFD(30),CGD(30),CHD(30),CID(30),                    &
                    AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,                     &
                    DX,DY,DZ,RCD,NGA,NGB,NGC,NGD                                      
      COMMON/SHLNOS/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,                &
                    MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,            &
                    NIJ,IJ,KL                                          
      COMMON/SHLNOS1/QQ4,IJKL
!      
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
      COMMON/EXCaux/EXaux(2000),Caux(2000)
!
      LOGICAL FLIP
      INTEGER,DIMENSION(NSHELL)::KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DIMENSION IX(84),IY(84),IZ(84),JX(84),JY(84),JZ(84),              &
                KX(84),KY(84),KZ(84),LX(84),LY(84),LZ(84)
      INTEGER,DIMENSION(8) :: MINF,MAXF                                  
      DATA MINF / 1, 2,  5, 11, 21, 36, 57, 1/                           
      DATA MAXF / 1, 4, 10, 20, 35, 56, 84, 4/                           
      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,       &
                  3,   0,   0,   2,   2,   1,   0,   1,   0,   1,       &
                  4,   0,   0,   3,   3,   1,   0,   1,   0,   2,       &
                  2,   0,   2,   1,   1,                                &
                  5,   0,   0,   4,   4,   1,   0,   1,   0,   3,       &
                  3,   2,   0,   2,   0,   3,   1,   1,   2,   2,       &
                  1,                                                    &
                  6,   0,   0,   5,   5,   1,   0,   1,   0,   4,       &
                  4,   2,   0,   2,   0,   4,   1,   1,   3,   3,       &
                  0,   3,   3,   2,   1,   2,   1,   2/                  
      DATA KX /   0,   7,   0,   0,  14,   0,   0,   7,   7,   0,       &
                 21,   0,   0,  14,  14,   7,   0,   7,   0,   7,       &
                 28,   0,   0,  21,  21,   7,   0,   7,   0,  14,       &
                 14,   0,  14,   7,   7,                                &
                 35,   0,   0,  28,  28,   7,   0,   7,   0,  21,       &
                 21,  14,   0,  14,   0,  21,   7,   7,  14,  14,       &
                  7,                                                    &
                 42,   0,   0,  35,  35,   7,   0,   7,   0,  28,       &
                 28,  14,   0,  14,   0,  28,   7,   7,  21,  21,       &
                  0,  21,  21,  14,   7,  14,   7,  14/                  
      DATA JX /   0,  49,   0,   0,  98,   0,   0,  49,  49,   0,       &
                147,   0,   0,  98,  98,  49,   0,  49,   0,  49,       &
                196,   0,   0, 147, 147,  49,   0,  49,   0,  98,       &
                 98,   0,  98,  49,  49,                                &
                245,   0,   0, 196, 196,  49,   0,  49,   0, 147,       &
                147,  98,   0,  98,   0, 147,  49,  49,  98,  98,       &
                 49,                                                    &
                294,   0,   0, 245, 245,  49,   0,  49,   0, 196,       &
                196,  98,   0,  98,   0, 196,  49,  49, 147, 147,       &
                  0, 147, 147,  98,  49,  98,  49,  98/                  
      DATA IX /   1, 344,   1,   1, 687,   1,   1, 344, 344,   1,       &
               1030,   1,   1, 687, 687, 344,   1, 344,   1, 344,       &
               1373,   1,   1,1030,1030, 344,   1, 344,   1, 687,       &
                687,   1, 687, 344, 344,                                &
               1716,   1,   1,1373,1373, 344,   1, 344,   1,1030,       &
               1030, 687,   1, 687,   1,1030, 344, 344, 687, 687,       &
                344,                                                    &
               2059,   1,   1,1716,1716, 344,   1, 344,   1,1373,       &
               1373, 687,   1, 687,   1,1373, 344, 344,1030,1030,       &
                  1,1030,1030, 687, 344, 687, 344, 687/                  
      DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,       &
                  0,   3,   0,   1,   0,   2,   2,   0,   1,   1,       &
                  0,   4,   0,   1,   0,   3,   3,   0,   1,   2,       &
                  0,   2,   1,   2,   1,                                &
                  0,   5,   0,   1,   0,   4,   4,   0,   1,   2,       &
                  0,   3,   3,   0,   2,   1,   3,   1,   2,   1,       &
                  2,                                                    &
                  0,   6,   0,   1,   0,   5,   5,   0,   1,   2,       &
                  0,   4,   4,   0,   2,   1,   4,   1,   3,   0,       &
                  3,   2,   1,   3,   3,   1,   2,   2/                  
      DATA KY /   0,   0,   7,   0,   0,  14,   0,   7,   0,   7,       &
                  0,  21,   0,   7,   0,  14,  14,   0,   7,   7,       &
                  0,  28,   0,   7,   0,  21,  21,   0,   7,  14,       &
                  0,  14,   7,  14,   7,                                &
                  0,  35,   0,   7,   0,  28,  28,   0,   7,  14,       &
                  0,  21,  21,   0,  14,   7,  21,   7,  14,   7,       &
                 14,                                                    &
                  0,  42,   0,   7,   0,  35,  35,   0,   7,  14,       &
                  0,  28,  28,   0,  14,   7,  28,   7,  21,   0,       &
                 21,  14,   7,  21,  21,   7,  14,  14/                  
      DATA JY /   0,   0,  49,   0,   0,  98,   0,  49,   0,  49,       &
                  0, 147,   0,  49,   0,  98,  98,   0,  49,  49,       &
                  0, 196,   0,  49,   0, 147, 147,   0,  49,  98,       &
                  0,  98,  49,  98,  49,                                &
                  0, 245,   0,  49,   0, 196, 196,   0,  49,  98,       &
                  0, 147, 147,   0,  98,  49, 147,  49,  98,  49,       &
                 98,                                                    &
                  0, 294,   0,  49,   0, 245, 245,   0,  49,  98,       &
                  0, 196, 196,   0,  98,  49, 196,  49, 147,   0,       &
                147,  98,  49, 147, 147,  49,  98,  98/                  
      DATA IY /   1,   1, 344,   1,   1, 687,   1, 344,   1, 344,       &
                  1,1030,   1, 344,   1, 687, 687,   1, 344, 344,       &
                  1,1373,   1, 344,   1,1030,1030,   1, 344, 687,       &
                  1, 687, 344, 687, 344,                                &
                  1,1716,   1, 344,   1,1373,1373,   1, 344, 687,       &
                  1,1030,1030,   1, 687, 344,1030, 344, 687, 344,       &
                687,                                                    &
                  1,2059,   1, 344,   1,1716,1716,   1, 344, 687,       &
                  1,1373,1373,   1, 687, 344,1373, 344,1030,   1,       &
               1030, 687, 344,1030,1030, 344, 687, 687/                  
      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,       &
                  0,   0,   3,   0,   1,   0,   1,   2,   2,   1,       &
                  0,   0,   4,   0,   1,   0,   1,   3,   3,   0,       &
                  2,   2,   1,   1,   2,                                &
                  0,   0,   5,   0,   1,   0,   1,   4,   4,   0,       &
                  2,   0,   2,   3,   3,   1,   1,   3,   1,   2,       &
                  2,                                                    &
                  0,   0,   6,   0,   1,   0,   1,   5,   5,   0,       &
                  2,   0,   2,   4,   4,   1,   1,   4,   0,   3,       &
                  3,   1,   2,   1,   2,   3,   3,   2/                  
      DATA KZ /   0,   0,   0,   7,   0,   0,  14,   0,   7,   7,       &
                  0,   0,  21,   0,   7,   0,   7,  14,  14,   7,       &
                  0,   0,  28,   0,   7,   0,   7,  21,  21,   0,       &
                 14,  14,   7,   7,  14,                                &
                  0,   0,  35,   0,   7,   0,   7,  28,  28,   0,       &
                 14,   0,  14,  21,  21,   7,   7,  21,   7,  14,       &
                 14,                                                    &
                  0,   0,  42,   0,   7,   0,   7,  35,  35,   0,       &
                 14,   0,  14,  28,  28,   7,   7,  28,   0,  21,       &
                 21,   7,  14,   7,  14,  21,  21,  14/                  
      DATA JZ /   0,   0,   0,  49,   0,   0,  98,   0,  49,  49,       &
                  0,   0, 147,   0,  49,   0,  49,  98,  98,  49,       &
                  0,   0, 196,   0,  49,   0,  49, 147, 147,   0,       &
                 98,  98,  49,  49,  98,                                &
                  0,   0, 245,   0,  49,   0,  49, 196, 196,   0,       &
                 98,   0,  98, 147, 147,  49,  49, 147,  49,  98,       &
                 98,                                                    &
                  0,   0, 294,   0,  49,   0,  49, 245, 245,   0,       &
                 98,   0,  98, 196, 196,  49,  49, 196,   0, 147,       &
                147,  49,  98,  49,  98, 147, 147,  98/                  
      DATA IZ /   1,   1,   1, 344,   1,   1, 687,   1, 344, 344,       &
                  1,   1,1030,   1, 344,   1, 344, 687, 687, 344,       &
                  1,   1,1373,   1, 344,   1, 344,1030,1030,   1,       &
                687, 687, 344, 344, 687,                                &
                  1,   1,1716,   1, 344,   1, 344,1373,1373,   1,       &
                687,   1, 687,1030,1030, 344, 344,1030, 344, 687,       &
                687,                                                    &
                  1,   1,2059,   1, 344,   1, 344,1716,1716,   1,       &
                687,   1, 687,1373,1373, 344, 344,1373,   1,1030,       &
               1030, 344, 687, 344, 687,1030,1030, 687/                 
!-----------------------------------------------------------------------                                                                       
!     PREPARE SHELL INFORMATION/FOR HONDO INTEGRATION                   
      IF(NELEC==2) GO TO 200      
      NORM = .TRUE.  
!                                                                       
!     ----- PERMUTE ISH AND JSH SHELLS, FOR THEIR TYPE                  
!     THIS IS DONE FOR SPEED REASONS.  THE CODE GETS THE RIGHT ANSWER   
!     WITHOUT THE ANGULAR MOMENTUM FLIPPING, AND THEREFORE A CALLING    
!     ARGUMENT ALLOWS ONE DO EXACTLY THE integral BLOCK AS SPECIFIED,   
!     SHOULD THAT BE DESIRED.                                           
!                                                                       
      IANDJ = ISH == JSH                                              
      IF (KTYPE(ISH) < KTYPE(JSH)  .and.  FLIP) THEN                 
       INU = JSH                                                      
       JNU = ISH                                                      
       NGTI = NGTH(2)                                                 
       NGTJ = NGTH(1)                                                 
      ELSE                                                              
       INU = ISH                                                      
       JNU = JSH                                                      
       NGTI = NGTH(1)                                                 
       NGTJ = NGTH(2)                                                 
      END IF                                                            
!                                                                       
!     ----- ISHELL                                                      
!                                                                       
      I = KATOM(INU)                                                    
      AX = Cxyz(1,I)                                                       
      AY = Cxyz(2,I)                                                       
      AZ = Cxyz(3,I)                                                       
      I1 = KSTART(INU)                                                  
      I2 = I1+KNG(INU)-1                                                
      LIT = KTYPE(INU)                                                  
      MINI = KMIN(INU)                                                  
      MAXI = KMAX(INU)                                                  
      LOCI = KLOC(INU)-MINI                                             
      NGA = 0                                                           
      DO I = I1,I2                                                  
       NGA = NGA+1                                                    
       GA(NGA) = EX(I)                                                
       CSA(NGA) = CS(I)                                               
       CPA(NGA) = CP(I)                                               
       CDA(NGA) = CD(I)                                               
       CFA(NGA) = CF(I)                                               
       CGA(NGA) = CG(I)                                               
       CHA(NGA) = CH(I)                                               
       CIA(NGA) = CI(I)                                               
      END DO                                                          
!                                                                       
!     ----- JSHELL                                                      
!                                                                       
      J = KATOM(JNU)                                                    
      BX = Cxyz(1,J)                                                       
      BY = Cxyz(2,J)                                                       
      BZ = Cxyz(3,J)                                                       
      J1 = KSTART(JNU)                                                  
      J2 = J1+KNG(JNU)-1                                                
      LJT = KTYPE(JNU)                                                  
      MINJ = KMIN(JNU)                                                  
      MAXJ = KMAX(JNU)                                                  
      LOCJ = KLOC(JNU)-MINJ                                             
      NGB = 0                                                           
      DO J = J1,J2                                                  
       NGB = NGB+1                                                    
       GB(NGB) = EX(J)                                                
       CSB(NGB) = CS(J)                                               
       CPB(NGB) = CP(J)                                               
       CDB(NGB) = CD(J)                                               
       CFB(NGB) = CF(J)                                               
       CGB(NGB) = CG(J)                                               
       CHB(NGB) = CH(J)                                               
       CIB(NGB) = CI(J)                                               
      END DO                                                          
      RAB = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))       
!                                                                       
!     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS                
!                                                                       
      IJ = 0                                                            
      JMAX = MAXJ                                                       
      DO I = MINI,MAXI                                              
       NX = IX(I)                                                     
       NY = IY(I)                                                     
       NZ = IZ(I)                                                     
       IF (IANDJ) JMAX = I                                            
       DO J = MINJ,JMAX                                           
        IJ = IJ+1                                                   
        IJX(IJ) = NX+JX(J)                                          
        IJY(IJ) = NY+JY(J)                                          
        IJZ(IJ) = NZ+JZ(J)                                          
        IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1                    
       END DO                                                   
      END DO                                                         
      RETURN                                                            
!     ******                                                            
!                                                                       
!        K AND L SHELL                                                  
!                                                                       
  200 CONTINUE    
      NORM  = .FALSE.      
      KANDL = .FALSE.                                 
      SAME  = .FALSE.          

      KNU = KSH                                                      
      LNU = 1                                                      
      NGTK = NGTH(3)                                                 
      NGTL = NGTH(4)                                                 
!                                                                       
!     ----- K SHELL                                                     
!                                                                       
      K = KATOMaux(KNU)                                                 
      CX = Cxyz(1,K)                                                       
      CY = Cxyz(2,K)                                                       
      CZ = Cxyz(3,K)                                                       
      K1 = KSTARTAUX(KNU)   !
      K2 = K1+KNGAUX(KNU)-1 !
      LKT = KTYPEaux(KNU)                                               
      MINK = MINF(LKT)                                                  
      MAXK = MAXF(LKT)                                                  
      LOCK = KLOCaux(KNU)-MINK                                          
      NGC = 0                                                           
      DO K = K1,K2                                                  
       NGC = NGC+1                                                    
       GC(NGC) = EXaux(K)                                             
       CSC(NGC) = Caux(K)                                            
       CPC(NGC) = Caux(K)                                         
       CDC(NGC) = Caux(K)                                             
       CFC(NGC) = Caux(K)                                             
       CGC(NGC) = Caux(K)                                             
       CHC(NGC) = Caux(K)                                             
       CIC(NGC) = Caux(K)                                             
      END DO                                                          
!                                                                       
!     ----- LSHELL (Unity Shell)                                                     
!                                                                       
      L = 1                                                    
      DX = 0                                                       
      DY = 0                                                       
      DZ = 0                                                       
      L1 = 1                                                
      L2 = 1                                                
      LLT = 1                                                  
      MINL = 1                                                  
      MAXL = 1                                                  
      LOCL = 1                                             
      NGD = 0                                                           
      DO L = L1,L2                                                  
       NGD = NGD+1                                                    
       GD(NGD) = 0                                                
       CSD(NGD) = 1                                             
       CPD(NGD) = 1                                              
       CDD(NGD) = 1                                               
       CFD(NGD) = 1                                               
       CGD(NGD) = 1                                               
       CHD(NGD) = 1                                               
       CID(NGD) = 1                                               
      END DO                                                          
      NROOTS = (LIT+LJT+LKT+LLT-2)/2                                    
      RCD = ((CX-DX)*(CX-DX) + (CY-DY)*(CY-DY) + (CZ-DZ)*(CZ-DZ))       
!                                                                       
!     ----- PREPARE INDICES FOR PAIRS OF (K,L) FUNCTIONS                
!                                                                       
      KL = 0                                                            
      LMAX = MAXL                                                       
      DO K = MINK,MAXK                                              
       NX = KX(K)                                                     
       NY = KY(K)                                                     
       NZ = KZ(K)                                                     
       IF (KANDL) LMAX = K                                            
       DO L = MINL,LMAX                                           
        KL = KL+1                                                   
        KLX(KL) = NX+LX(L)                                          
        KLY(KL) = NY+LY(L)                                          
        KLZ(KL) = NZ+LZ(L)                                          
        KLGT(KL) = NGTK*(K-MINK)+NGTL*(L-MINL)                      
       END DO                                                       
      END DO                                                          
      MAX = KL                                                          
      DO 320 I = 1,IJ                                                   
      IF (SAME) MAX = I                                                 
  320 IK(I) = MAX                                                       
      IJKL = IJ*KL                                                      
      IF (SAME) IJKL = IJ*(IJ+1)/2

      RETURN                                                            
      END       

! PDPT_msqrt
      SUBROUTINE PDPT_msqrt(GMAT,N,IPRINTOPT)
      INTEGER :: N,IPRINTOPT,I,J,NTRUNC
      DOUBLE PRECISION :: GMAT(N,N),TOL
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: VEC,AUX,AUX2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: EIG,W

      ALLOCATE(VEC(N,N),AUX(N,N),AUX2(N,N))
      ALLOCATE(EIG(N),W(N))
!
      TOL = 1.0d-10
      AUX = GMAT
      CALL DIAG(N,AUX,VEC,EIG,W)
!
      NTRUNC = 0
      AUX = 0.0d0
      DO I=1,N
       IF (EIG(I)<=TOL) THEN
        NTRUNC = NTRUNC + 1
        EIG(I) = 0.0D0
        DO J=1,N
         VEC(J,I) = 0.0d0
        END DO
       ELSE
        AUX(I,I) = 1.0/DSQRT(EIG(I))
       END IF
      END DO

      IF(IPRINTOPT==1.AND.NTRUNC>0) THEN
        WRITE(6,*) "RI Warning - Number of values truncated from metric:",NTRUNC
      END IF

      DEALLOCATE(EIG,W)
      CALL DGEMM("N","N",N,N,N,1.0D0,VEC,N,AUX,N,0.0D0,AUX2,N)
      CALL DGEMM("N","T",N,N,N,1.0D0,AUX2,N,VEC,N,0.0D0,GMAT,N)
      DEALLOCATE(VEC,AUX,AUX2)

      END

! QOUT2C                                             
      SUBROUTINE QOUT2C(GMAT,NBFaux,GHONDO,MAXG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL IANDJ,KANDL,SAME
      COMMON/MISC  /IANDJ,KANDL,SAME                                  
      COMMON/ERIOUT/ISH,JSH,KSH,LSH,LSTRI,LSTRJ,LSTRK,LSTRL           
      COMMON/SHLNOS/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,                &
                    MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,            &
                    NIJ,IJ,KL                                          
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
!
      DOUBLE PRECISION,DIMENSION(NBFaux,NBFaux) :: GMAT
      DOUBLE PRECISION,DIMENSION(MAXG) :: GHONDO
      INTEGER,DIMENSION(8) :: MINF,MAXF
      DATA MINF / 1, 2,  5, 11, 21, 36, 57, 1/                 
      DATA MAXF / 1, 4, 10, 20, 35, 56, 84, 4/                 
!-----------------------------------------------------------------------                                                                       
!     Get Auxiliary Basis Info. 
!-----------------------------------------------------------------------
      SAME  = ISH == KSH                            
      MINI = MINF(KTYPEaux(ISH))                                        
      MINK = MINF(KTYPEaux(KSH))                                        
      MAXI = MAXF(KTYPEaux(ISH))                                        
      MAXK = MAXF(KTYPEaux(KSH))                                        
      LOCI = KLOCaux(ISH)-MINI                                         
      LOCK = KLOCaux(KSH)-MINK                                          
!-----------------------------------------------------------------------                                                                       
!     Store (k|l) ERIs in G matrix 
!-----------------------------------------------------------------------
      DO 1 I = MINI,MAXI                                              
        I_INDEX = (I-MINI)*LSTRI + 1                                   
        DO K =  MINK,MAXK                                       
          IK_INDEX = (K-MINK)*LSTRK + I_INDEX                    
          IF(SAME.and.K>I)GO TO 1                 
          VAL = GHONDO( IK_INDEX )
          I1 = LOCI+I                                           
          I3 = LOCK+K                                           
          GMAT(I1,I3) = VAL
          GMAT(I3,I1) = VAL
        END DO
    1 CONTINUE                                                       
!-----------------------------------------------------------------------                                                                       
      RETURN                                                            
      END                                                               

! QOUT3C                                     
      SUBROUTINE QOUT3C(BUFP2,GHONDO,MAXG,GMAT,NBF,NBFaux,              &
                        KLOC,KMIN,KMAX,NSHELL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL IANDJ,KANDL,SAME
      COMMON/MISC/IANDJ,KANDL,SAME                                  
      COMMON/ERIOUT/ISH,JSH,KSH,LSH,LSTRI,LSTRJ,LSTRK,LSTRL           
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
!
      INTEGER,DIMENSION(NSHELL) :: KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NBF*(NBF+1)/2,NBFaux) :: BUFP2
      DOUBLE PRECISION,DIMENSION(NBFaux,NBFaux)::GMAT
      DOUBLE PRECISION,DIMENSION(MAXG)::GHONDO
      INTEGER,DIMENSION(8) :: MINF,MAXF
      DATA MINF / 1, 2,  5, 11, 21, 36, 57, 1/                 
      DATA MAXF / 1, 4, 10, 20, 35, 56, 84, 4/      
      INTEGER :: T
!-----------------------------------------------------------------------                                                                       
!     Get Basis and Auxiliary Basis Info. 
!-----------------------------------------------------------------------
      IANDJ = ISH == JSH                                              
      MINI = KMIN(ISH)                                                  
      MINJ = KMIN(JSH)                                                  
      MINK = MINF(KTYPEaux(KSH))                                        
      MAXI = KMAX(ISH)                                                  
      MAXJ = KMAX(JSH)                                                  
      MAXK = MAXF(KTYPEaux(KSH))                                        
      LOCI = KLOC(ISH)-MINI                                             
      LOCJ = KLOC(JSH)-MINJ                                             
      LOCK = KLOCaux(KSH)-MINK                                          
!-----------------------------------------------------------------------                                                                       
!     Build B tensor (Only Upper Triangular)
!      b_mn^l = sum_k (mn|k) G^{-1/2}_{kl}
!-----------------------------------------------------------------------
      JMAX = MAXJ                                                       
      DO I = MINI,MAXI                                              
       I_INDEX = (I-MINI)*LSTRI + 1                                   
       IF (IANDJ) JMAX = I                                            
       DO J = MINJ,JMAX                                           
        IJ_INDEX = (J-MINJ)*LSTRJ + I_INDEX                         
        DO K = MINK,MAXK                                       
          IJK_INDEX = (K-MINK)*LSTRK + IJ_INDEX                    
          VAL = GHONDO( IJK_INDEX ) 
          M = LOCI+I                                           
          N = LOCJ+J                                           
          KK = LOCK+K                                           
          IF (M <= N) GO TO 100                             
          T = M                                                
          M = N                                               
          N = T                                                
  100     MN = M + N*(N-1)/2
          DO L=1,NBFaux
           BUFP2(MN,L) = BUFP2(MN,L) + VAL*GMAT(KK,L)
          END DO
        END DO
       END DO                                                       
      END DO
!-----------------------------------------------------------------------                                                                       
      RETURN 
      END                                                               

! QOUT3CModChol
      SUBROUTINE QOUT3CModChol(BUFP2,GHONDO,MAXG,NBF,NBFaux,            &
                               KLOC,KMIN,KMAX,NSHELL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IANDJ,KANDL,SAME
      COMMON/MISC/IANDJ,KANDL,SAME
      COMMON/ERIOUT/ISH,JSH,KSH,LSH,LSTRI,LSTRJ,LSTRK,LSTRL
      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
!
      INTEGER,DIMENSION(NSHELL) :: KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NBF*(NBF+1)/2,NBFaux) :: BUFP2
      DOUBLE PRECISION,DIMENSION(MAXG)::GHONDO
      INTEGER,DIMENSION(8) :: MINF,MAXF
      DATA MINF / 1, 2,  5, 11, 21, 36, 57, 1/
      DATA MAXF / 1, 4, 10, 20, 35, 56, 84, 4/
!-----------------------------------------------------------------------                                                                       
!     Get Basis and Auxiliary Basis Info. 
!-----------------------------------------------------------------------
      IANDJ = ISH == JSH
      MINI = KMIN(ISH)
      MINJ = KMIN(JSH)
      MINK = MINF(KTYPEaux(KSH))
      MAXI = KMAX(ISH)
      MAXJ = KMAX(JSH)
      MAXK = MAXF(KTYPEaux(KSH))
      LOCI = KLOC(ISH)-MINI
      LOCJ = KLOC(JSH)-MINJ
      LOCK = KLOCaux(KSH)-MINK
!-----------------------------------------------------------------------                                                                       
!     Store (mn|P) (Only Upper Triangular)
!-----------------------------------------------------------------------
      JMAX = MAXJ
      DO I = MINI,MAXI
       I_INDEX = (I-MINI)*LSTRI + 1
       IF (IANDJ) JMAX = I
       DO 1 J = MINJ,JMAX
        IJ_INDEX = (J-MINJ)*LSTRJ + I_INDEX
        DO K =  MINK,MAXK
          IJK_INDEX = (K-MINK)*LSTRK + IJ_INDEX
          VAL = GHONDO( IJK_INDEX )
          M = LOCI+I
          N = LOCJ+J
          KK = LOCK+K
          IF (M <= N) GO TO 100
          MT = M
          M = N
          N = MT
  100     MN = M + N*(N-1)/2
          BUFP2(MN,KK) = VAL
        END DO
    1  CONTINUE
      END DO
!-----------------------------------------------------------------------                                                                       
      RETURN
      END

! QOUTaux
      SUBROUTINE QOUTaux(BUFP2,GLIBRETA,NBF,NBFaux,ISH,JSH,KLOC,NSHELL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER,DIMENSION(NSHELL) :: KLOC
      DIMENSION GLIBRETA(*)
      DOUBLE PRECISION,DIMENSION(NBF*(NBF+1)/2,NBFaux) :: BUFP2
      INTEGER :: ORII,ORIJ,ORIK
      COMMON/ORI/ORII,ORIJ,ORIK
      INTEGER :: ISH,JSH
      LOGICAL IANDJ
!-----------------------------------------------------------------------
!     Pack 4-indices into 1 word.
!     Write Label & integral on Unit 1 if DONTW = .False.
!-----------------------------------------------------------------------
      IANDJ = ISH == JSH
      LOCI = KLOC(ISH)
      LOCJ = KLOC(JSH)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IJN = 0
      JMAX = ORIJ
      DO I = 1,ORII
       I_INDEX = (I-1)*ORIJ*NBFaux
       IF (IANDJ) JMAX = I
       DO J = 1,JMAX
        IJ_INDEX = (J-1)*NBFaux + I_INDEX
        IJN = IJN+1
        KLN = 0
        DO K = 1,NBFaux
         IJK_INDEX = (K-1) + IJ_INDEX + 1
         VAL = GLIBRETA( IJK_INDEX )
         I1 = LOCI+I-1
         I2 = LOCJ+J-1
         I3 = K
         IF (I1 <= I2) GO TO 100
         N = I1
         I1 = I2
         I2 = N
  100    CONTINUE
         MN = I1 + I2*(I2-1)/2
         BUFP2(MN,K) = VAL
        END DO
       END DO
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

!----------------------------------------------------------------------!
!                                                                      !
!                    Electrostatic Moment Integrals                    !     
!                                                                      !
!----------------------------------------------------------------------!

! INTMOM                                           
      SUBROUTINE INTMOM(LIT1,LJT1,IJ,IJX,IJY,IJZ,DIJ,WINT,AA,AX,AY,AZ)    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/ELPROP/IEMOM      
      COMMON/PRPINT/XINT0,XINT1,XINT2,XINT3,YINT0,YINT1,YINT2,YINT3,    &
                    ZINT0,ZINT1,ZINT2,ZINT3
      COMMON/HERMIT/H1(55),W1(55)
      COMMON/XYZORB/TXYZ,X00,Y00,Z00,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ      
      DIMENSION DIJ(*),WINT(*),IJX(*),IJY(*),IJZ(*)                     
      DIMENSION XIN(49,4),YIN(49,4),ZIN(49,4)                           
!-----------------------------------------------------------------------
!                    Electrostatic Moment Integrals 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      TXYZ  = 1.0d0/SQRT(AA)                                                 
      X00 = AX                                                           
      Y00 = AY                                                           
      Z00 = AZ                                                           
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IN = -7                                                           
      DO I=1,LIT1                                                   
       IN = IN+7                                                         
       NI = I                                                            
       DO J=1,LJT1                                                   
        JN = IN+J                                                         
        NJ = J                                                            
        CALL EMOMINT(H1,W1)                                                       
        XIN(JN,1) = XINT0*TXYZ                                               
        YIN(JN,1) = YINT0*TXYZ                                               
        ZIN(JN,1) = ZINT0*TXYZ                                               
        XIN(JN,2) = XINT1*TXYZ                                               
        YIN(JN,2) = YINT1*TXYZ                                               
        ZIN(JN,2) = ZINT1*TXYZ                                               
        IF(IEMOM>=2) THEN                                                  
         XIN(JN,3) = XINT2*TXYZ                                            
         YIN(JN,3) = YINT2*TXYZ                                            
         ZIN(JN,3) = ZINT2*TXYZ                                            
        END IF                                                            
        IF(IEMOM>=3) THEN                                                  
         XIN(JN,4) = XINT3*TXYZ                                            
         YIN(JN,4) = YINT3*TXYZ                                            
         ZIN(JN,4) = ZINT3*TXYZ                                            
        END IF                                                            
       END DO
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I=1,IJ                                                     
       NIX = IJX(I)                                                       
       NIY = IJY(I)                                                       
       NIZ = IJZ(I)                                                       
       DENS = DIJ(I)                                                      
       XIN0 = XIN(NIX,1)                                                  
       YIN0 = YIN(NIY,1)                                                  
       ZIN0 = ZIN(NIZ,1)                                                  
       XIN1 = XIN(NIX,2)                                                   
       YIN1 = YIN(NIY,2)                                                
       ZIN1 = ZIN(NIZ,2)                                                
       INDEX = I                                                       
       WINT(INDEX) = WINT(INDEX) + XIN1 * YIN0 * ZIN0 * DENS              
       INDEX = INDEX + IJ                                                
       WINT(INDEX) = WINT(INDEX) + XIN0 * YIN1 * ZIN0 * DENS             
       INDEX = INDEX + IJ                                                
       WINT(INDEX) = WINT(INDEX) + XIN0 * YIN0 * ZIN1 * DENS             
       IF(IEMOM>=2) THEN                                                  
        XIN2 = XIN(NIX,3)                                        
        YIN2 = YIN(NIY,3)                                                
        ZIN2 = ZIN(NIZ,3)                                                
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN2 * YIN0 * ZIN0 * DENS           
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN0 * YIN2 * ZIN0 * DENS           
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN0 * YIN0 * ZIN2 * DENS           
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN1 * YIN1 * ZIN0 * DENS           
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN1 * YIN0 * ZIN1 * DENS           
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN0 * YIN1 * ZIN1 * DENS           
       END IF                                                            
       IF(IEMOM>=3) THEN                                                  
        XIN3 = XIN(NIX,4)                                                
        YIN3 = YIN(NIY,4)                                                
        ZIN3 = ZIN(NIZ,4)                                                
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN3 * YIN0 * ZIN0 * DENS           
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN0 * YIN3 * ZIN0 * DENS           
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN0 * YIN0 * ZIN3 * DENS           
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN2 * YIN1 * ZIN0 * DENS           
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN2 * YIN0 * ZIN1 * DENS           
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN1 * YIN2 * ZIN0 * DENS           
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN0 * YIN2 * ZIN1 * DENS           
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN1 * YIN0 * ZIN2 * DENS           
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN0 * YIN1 * ZIN2 * DENS           
        INDEX = INDEX + IJ                                              
        WINT(INDEX) = WINT(INDEX) + XIN1 * YIN1 * ZIN1 * DENS           
       END IF                                                            
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN                                                            
      END                                                               

! EMOMINT                                           
      SUBROUTINE EMOMINT(H,WW)                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/ELPROP/IEMOM      
      COMMON/PRPINT/XINT0,XINT1,XINT2,XINT3,YINT0,YINT1,YINT2,YINT3,    &
                    ZINT0,ZINT1,ZINT2,ZINT3                              
      COMMON/XYZORB/TXYZ,X00,Y00,Z00,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
      COMMON/CMCoord/Xcm,Ycm,Zcm
      DIMENSION H(36),WW(36),MINARRAY(8),MAXARRAY(8)                               
      DATA MINARRAY /1,2,4, 7,11,16,22,29/                                   
      DATA MAXARRAY /1,3,6,10,15,21,28,36/
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
!     Gauss-Hermite Quadrature
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      XINT0 = 0.0D+00                                                      
      YINT0 = 0.0D+00                                                      
      ZINT0 = 0.0D+00                                                      
      XINT1 = 0.0D+00                                                      
      YINT1 = 0.0D+00                                                      
      ZINT1 = 0.0D+00                                                      
      XINT2 = 0.0D+00                                                      
      YINT2 = 0.0D+00                                                      
      ZINT2 = 0.0D+00                                                      
      XINT3 = 0.0D+00                                                      
      YINT3 = 0.0D+00                                                      
      ZINT3 = 0.0D+00                                                      
      NPTS = (NI + NJ + IEMOM - 2)/2 + 1                                   
      IMIN = MINARRAY(NPTS)                                                  
      IMAX = MAXARRAY(NPTS)                                                  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I = IMIN,IMAX                                              
       DUM = WW(I)                                                     
       PX = DUM                                                       
       PY = DUM                                                       
       PZ = DUM
       DUM = H(I)*TXYZ                                                   
       PTX = DUM + X00                                                 
       PTY = DUM + Y00                                                 
       PTZ = DUM + Z00                                                 
       AX = PTX - XI                                                  
       AY = PTY - YI                                                  
       AZ = PTZ - ZI                                                  
       BX = PTX - XJ                                                  
       BY = PTY - YJ                                                  
       BZ = PTZ - ZJ                                                  
       CX = PTX - Xcm                                                  
       CY = PTY - Ycm                                                  
       CZ = PTZ - Zcm                                                  
!                                                                     
       GO TO (101,102,103,104,105,106,107), NI                        
  107  PX = PX*AX                                                     
       PY = PY*AY                                                     
       PZ = PZ*AZ                                                     
  106  PX = PX*AX                                                     
       PY = PY*AY                                                     
       PZ = PZ*AZ                                                     
  105  PX = PX*AX                                                     
       PY = PY*AY                                                     
       PZ = PZ*AZ                                                     
  104  PX = PX*AX                                                     
       PY = PY*AY                                                     
       PZ = PZ*AZ                                                     
  103  PX = PX*AX                                                     
       PY = PY*AY                                                     
       PZ = PZ*AZ                                                     
  102  PX = PX*AX                                                     
       PY = PY*AY                                                     
       PZ = PZ*AZ                                                     
  101  CONTINUE                                                       
!                                                                       
       GO TO (201,202,203,204,205,206,207), NJ                        
  207  PX = PX*BX                                                     
       PY = PY*BY                                                     
       PZ = PZ*BZ                                                     
  206  PX = PX*BX                                                     
       PY = PY*BY                                                     
       PZ = PZ*BZ                                                     
  205  PX = PX*BX                                                     
       PY = PY*BY                                                     
       PZ = PZ*BZ                                                     
  204  PX = PX*BX                                                     
       PY = PY*BY                                                     
       PZ = PZ*BZ                                                     
  203  PX = PX*BX                                                     
       PY = PY*BY                                                     
       PZ = PZ*BZ                                                     
  202  PX = PX*BX                                                     
       PY = PY*BY                                                     
       PZ = PZ*BZ                                                     
  201  CONTINUE                                                       
!                                                                     
       IEMOM1 = IEMOM + 1
       GO TO (301,302,303,304),IEMOM1                                    
  304  XINT3 = XINT3 + PX*CX*CX*CX                                    
       YINT3 = YINT3 + PY*CY*CY*CY                                    
       ZINT3 = ZINT3 + PZ*CZ*CZ*CZ                                    
  303  XINT2 = XINT2 + PX*CX*CX                                       
       YINT2 = YINT2 + PY*CY*CY                                       
       ZINT2 = ZINT2 + PZ*CZ*CZ                                       
  302  XINT1 = XINT1 + PX*CX                                          
       YINT1 = YINT1 + PY*CY                                          
       ZINT1 = ZINT1 + PZ*CZ                                          
  301  XINT0 = XINT0 + PX                                             
       YINT0 = YINT0 + PY                                             
       ZINT0 = ZINT0 + PZ                                             
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN                                                            
      END                                                               
