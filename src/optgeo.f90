!======================================================================!
!                                                                      !
!  G E O M E T R Y   O P T I M I Z A T I O N   S U B R O U T I N E S   !
!                                                                      !
!======================================================================!
!                                                                      !
!   OPTLBFGS: Optimize geometries with LBFGS method (mthlib.f)         !
!   CALCOPTE: Routine to compute energy needed by OPTSUMSL             !
!   CALCOPTG: Routine to compute gradients needed by OPTSUMSL          !
!                                                                      !
!   OPTCGNAG: Optimize geometries using the CG NAG routine E04DGF      !
!   OPTSUMSL: Optimize geometries using the CG SUMSL routine           !
!   POINTERSOPT : Define Pointers of the USER array                    !
!   OPTIMIZE: Routine to carry out geometry optimization (IRUNTYP=3)   !
!   ENERGYFUN
!                                                                      !
!======================================================================!

! OPTIMIZE
      SUBROUTINE OPTIMIZE(NINTEG,IDONTW,NAT,ZAN,Cxyz,IAN,IMIN,IMAX,     &
         ZMASS,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,     &
         ISH,ITYP,C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI,DIPS,          &
         GRADS,IRUNTYP,IHSSCAL,IPROJECT,ISIGMA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL RESTART
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_CGM/ICGMETHOD
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/USELIBRETA/ILIBRETA
      COMMON/ELPROP/IEMOM
      !
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN,ZMASS
      INTEGER,DIMENSION(NAT) :: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMI)::C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT):: Cxyz
      DOUBLE PRECISION,DIMENSION(3*NAT):: GRADS
      DOUBLE PRECISION,DIMENSION(3):: DIPS
      !-----------------------------------------------------------------------
      EELEC_MIN = 1.0d20
      !     Generate an initial GCF if RESTART=F setting ICOEF=0
      IF (RESTART .eqv. .FALSE.) THEN
      ICOEFORI = ICOEF
      ICOEF = 0
      !      Update coordinates of shells if use libreta library for ERIs
      if(ILIBRETA==1)CALL UPDCOOSHELL(NSHELL,KATOM,Cxyz,NAT)
      CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
             ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,    &
             INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,   &
             CF,CG,CH,CI,GRADS,IRUNTYP,DIPS,0,0)
      RESTART = .TRUE.
      ICOEF = ICOEFORI
      ENDIF
      !     Select CG Method     
      IF(ICGMETHOD==1)THEN
      CALL OPTSUMSL(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
             ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,KATOM,KTYPE,   &
             KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP, &
             CD,CF,CG,CH,CI,IRUNTYP,GRADS,DIPS)                    
      ELSE IF(ICGMETHOD==2)THEN                                         
      CALL OPTCGNAG(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
             ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,KATOM,KTYPE,   &
             KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,    &
             CP,CD,CF,CG,CH,CI,IRUNTYP,GRADS,DIPS)                    
      ELSE IF(ICGMETHOD==3)THEN                                         
      CALL OPTLBFGS(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
             ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,    &
             INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,   &
             CF,CG,CH,CI,IRUNTYP,GRADS,DIPS,NPRINT)             
      ENDIF                                                             
      !     Compute Hessian from analytic gradients at stationary point       
      if(IHSSCAL==1)then                                                
      CALL HESSCAL(NINTEG,IDONTW,NAT,ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,     &
            KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,    &
            ITYP,C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI,DIPS,           &
            GRADS,IRUNTYP,IPROJECT,ISIGMA)
      end if
      !-----------------------------------------------------------------------
      RETURN
      END

! OPTSUMSL
      SUBROUTINE OPTSUMSL(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,    &
                          NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,   &
                          KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,     &
                          ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,IRUNTYP,   &
                          GRAD,DIPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      PARAMETER (BOHR = 0.52917724924D+00) 
      INTEGER,DIMENSION(NAT)    :: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
      INTEGER,INTENT(IN)        :: IRUNTYP
      DOUBLE PRECISION,DIMENSION(NAT)   :: ZAN,ZMASS
      DOUBLE PRECISION,DIMENSION(NPRIMI):: C1,C2,EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(3*NAT),INTENT(OUT) :: GRAD
      DOUBLE PRECISION,DIMENSION(3) :: DIPS
!      
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC_OLD,EELEC,DIF_EELEC,EELEC_MIN
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9,  &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,     &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,     &
                         NIU24,NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,NU7,     &
                         NU8,NU9,NU10,NU11,NU12,NU13,NULAST
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/GEOCONV/GEOMENERGY(200)
      COMMON/USELIBRETA/ILIBRETA 
!      
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IUSER,IV
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: USER,D,V     
      EXTERNAL CALCOPTE,CALCOPTG
!-----------------------------------------------------------------------
!     Define Pointers of the USER array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL POINTERSOPT(NAT,NSHELL,NPRIMI)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Transfer working arrays to IUSER and USER
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(IUSER(NIULAST),USER(NULAST))
!
      IUSER(NIU1)          = NINTEG
      IUSER(NIU2)          = IDONTW  
      IUSER(NIU3)          = 3
      IUSER(NIU4)          = IEMOM  
      IUSER(NIU5)          = NAT
      IUSER(NIU6)          = NBF
      IUSER(NIU7)          = NSHELL
      IUSER(NIU8)          = NPRIMI
      IUSER(NIU9 :NIU10-1) = IAN
      IUSER(NIU10:NIU11-1) = IMIN
      IUSER(NIU11:NIU12-1) = IMAX
      IUSER(NIU12:NIU13-1) = KSTART 
      IUSER(NIU13:NIU14-1) = KATOM
      IUSER(NIU14:NIU15-1) = KTYPE
      IUSER(NIU15:NIU16-1) = KLOC
      IUSER(NIU16:NIU17-1) = INTYP
      IUSER(NIU17:NIU18-1) = KNG
      IUSER(NIU18:NIU19-1) = KMIN
      IUSER(NIU19:NIU20-1) = KMAX
      IUSER(NIU20:NIU21-1) = ISH
      IUSER(NIU21:NIU22-1) = ITYP
      IUSER(NIU22)         = IRUNTYP
      IUSER(NIU23)         = 0             ! ITCG
      IUSER(NIU24)         = NBFaux      
!
      USER(NU1 :NU2-1)     = ZAN
      USER(NU2 :NU3-1)     = ZMASS 
      USER(NU3 :NU4-1)     = C1
      USER(NU4 :NU5-1)     = C2
      USER(NU5 :NU6-1)     = EX
      USER(NU6 :NU7-1)     = CS
      USER(NU7 :NU8-1)     = CP
      USER(NU8 :NU9-1)     = CD
      USER(NU9 :NU10-1)    = CF
      USER(NU10:NU11-1)    = CG
      USER(NU11:NU12-1)    = CH
      USER(NU12:NU13-1)    = CI
      USER(NU13:NULAST)    = DIPS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Initial Coordinates on File CGGRAD (Unit=11)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(11,'(/,72(1H*),/)')
      WRITE(11,*)'Initial Coordinates (Angs) for Geometry Optimization'
      WRITE(11,*)
      DO I=1,NAT
       WRITE(11,'(I5,3F15.4)')  &
             I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
      ENDDO
      WRITE(11,'(/,72(1H*),/)')      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Minimization of the total energy with respect to Cxyz
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,'(/,2X,A21,/1X,23(1H=))')'GEOMETRY OPTIMIZATION'
      WRITE(6,'(/,1X,A23,7X,A19,/)')'Call in CG Optimization',          &
                                    'Total Energy (a.u.)'
      IF(MOLDEN==1)WRITE(18,'(A18)')'[GEOMETRIES] (XYZ)'
!
      NV = 3*NAT
      LIV = 60
      LV = 71+NV*(NV+15)/2
!      
      ALLOCATE( D(NV),IV(LIV),V(LV) ) 
!      
      D(1:NV) = 1.0d0      
      IV = 0
      EMINIMA = 1.0d6
!
      CALL SUMSL(NV,D,Cxyz,CALCOPTE,CALCOPTG,IV,LIV,LV,V,IUSER,USER)
!      
      ENERGIA = EELEC + EN
      IF(EMINIMA>ENERGIA)EMINIMA = ENERGIA
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute energy at solution to print data on the output file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,'(/1X,A24,2X,F20.10,/)')'Final objective value = ',EMINIMA      
      WRITE(6,*)'New Coordinates after Conjugate Gradient Opt (Bohr)'
      WRITE(6,'(1X,51(1H-),/)')
      DO I=1,NAT
       WRITE(6,'(I5,3F15.4)')I,Cxyz(1,I),Cxyz(2,I),Cxyz(3,I)
      ENDDO
!     Update coordinates of shells if use libreta library for ERIs
      if(ILIBRETA==1)CALL UPDCOOSHELL(IUSER(NIU7),IUSER(NIU13),Cxyz,    &
                                      IUSER(NIU5))
      CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),IUSER(NIU5),    &
                    IUSER(NIU6),IUSER(NIU24),IUSER(NIU7),IUSER(NIU8),   &
                    USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),            &
                    IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),IUSER(NIU14),&
                    IUSER(NIU15),IUSER(NIU16),IUSER(NIU17),IUSER(NIU18),&
                    IUSER(NIU19),IUSER(NIU20),IUSER(NIU21),USER(NU3),   &
                    USER(NU4),USER(NU5),USER(NU6),USER(NU7),USER(NU8),  &
                    USER(NU9),USER(NU10),USER(NU11),USER(NU12),GRAD,    &                    
                    IUSER(NIU22),USER(NU13),0,1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Final Coordinates on File CGGRAD (Unit=100)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(11,'(/,72(1H*),/)')
      WRITE(11,*)'New Coordinates after Conjugate Gradient Opt (Angs)'
      WRITE(11,*)
      DO I=1,NAT
       WRITE(11,'(I5,3F15.4)')                                          &
            I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
      ENDDO
      CALL NUCDIST(NV,NAT,Cxyz)
      WRITE(11,'(/,72(1H*),/)')  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write geometry energies in [GEOCONV] section of the molden file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      IF(MOLDEN==1)WRITE(18,'(A9,/,A6)')'[GEOCONV]','energy'
      if(IUSER(NIU23)<=200)then
       DO I=1,IUSER(NIU23)
        WRITE(18,'(F20.10)')GEOMENERGY(I)
       END DO
      else
       write(6,*)'Sorry, No. of calls > 200, [GEOCONV] is not possible'
      end if
!-----------------------------------------------------------------------
      DEALLOCATE(IUSER,USER,D,IV,V)
      RETURN
      END

! POINTERSOPT
      SUBROUTINE POINTERSOPT(NAT,NSHELL,NPRIMI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9,  &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,     &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,     &
                         NIU24,NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,NU7,     &
                         NU8,NU9,NU10,NU11,NU12,NU13,NULAST
!-----------------------------------------------------------------------
!     Define Pointers of the USER array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NIU1  = 1                         ! NINTEG
      NIU2  = NIU1  + 1                 ! IDONTW 
      NIU3  = NIU2  + 1                 ! Not used
      NIU4  = NIU3  + 1                 ! IEMOM  
      NIU5  = NIU4  + 1                 ! NAT
      NIU6  = NIU5  + 1                 ! NBF
      NIU7  = NIU6  + 1                 ! NSHELL
      NIU8  = NIU7  + 1                 ! NPRIMI
      NIU9  = NIU8  + 1                 ! IAN
      NIU10 = NIU9  + NAT               ! IMIN
      NIU11 = NIU10 + NAT               ! IMAX
      NIU12 = NIU11 + NAT               ! KSTART 
      NIU13 = NIU12 + NSHELL            ! KATOM
      NIU14 = NIU13 + NSHELL            ! KTYPE
      NIU15 = NIU14 + NSHELL            ! KLOC
      NIU16 = NIU15 + NSHELL            ! INTYP
      NIU17 = NIU16 + NSHELL            ! KNG
      NIU18 = NIU17 + NSHELL            ! KMIN
      NIU19 = NIU18 + NSHELL            ! KMAX 
      NIU20 = NIU19 + NSHELL            ! ISH
      NIU21 = NIU20 + NPRIMI            ! ITYP
      NIU22 = NIU21 + NPRIMI            ! IRUNTYP
      NIU23 = NIU22 + 1                 ! ITCG
      NIU24 = NIU23 + 1                 ! NBFaux
      NIULAST = NIU24 + 1               
!
      NU1  = 1                          ! ZAN
      NU2  = NU1  + NAT                 ! ZMASS
      NU3  = NU2  + NAT                 ! C1
      NU4  = NU3  + NPRIMI              ! C2
      NU5  = NU4  + NPRIMI              ! EX   
      NU6  = NU5  + NPRIMI              ! CS
      NU7  = NU6  + NPRIMI              ! CP
      NU8  = NU7  + NPRIMI              ! CD
      NU9  = NU8  + NPRIMI              ! CF
      NU10 = NU9  + NPRIMI              ! CG
      NU11 = NU10 + NPRIMI              ! CH
      NU12 = NU11 + NPRIMI              ! CI
      NU13 = NU12 + NPRIMI              ! DIPS
      NULAST = NU13 + 3
!-----------------------------------------------------------------------
      RETURN
      END      

! CALCOPTE      
      SUBROUTINE CALCOPTE(NV,Cxyz,NF,ENERGIA,IUSER,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/USELIBRETA/ILIBRETA 
      COMMON/INPNOF_MOLDEN/MOLDEN      
      COMMON/GEOCONV/GEOMENERGY(200)
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9,  &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,     &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,     &
                         NIU24,NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,NU7,     &
                         NU8,NU9,NU10,NU11,NU12,NU13,NULAST
      !JFHLewYee: Changed NATOMS allowed dimension from 100 to 1000
      COMMON/ECP2/CLP(4004),ZLP(4004),NLP(4004),KFRST(1001,6),          &
                  KLAST(1001,6),LMAX(1001),LPSKIP(1001),IZCORE(1001)
!
      INTEGER,DIMENSION(NIULAST)         :: IUSER 
      DOUBLE PRECISION,DIMENSION(NV)     :: Cxyz
      DOUBLE PRECISION,DIMENSION(NULAST) :: USER
      PARAMETER (BOHR = 0.52917724924D+00) 
      CHARACTER*4,DIMENSION(106) :: ATMLAB
      DATA ATMLAB/'H   ','HE  ','LI  ','BE  ','B   ','C   ',            &
                  'N   ','O   ','F   ','NE  ','NA  ','MG  ',            &
                  'AL  ','SI  ','P   ','S   ','CL  ','AR  ',            &
                  'K   ','CA  ','SC  ','TI  ','V   ','CR  ',            &
                  'MN  ','FE  ','CO  ','NI  ','CU  ','ZN  ',            &
                  'GA  ','GE  ','AS  ','SE  ','BR  ','KR  ',            &
                  'RB  ','SR  ','Y   ','ZR  ','NB  ','MO  ',            &
                  'TC  ','RU  ','RH  ','PD  ','AG  ','CD  ',            &
                  'IN  ','SN  ','SB  ','TE  ','I   ','XE  ',            &
                  'CS  ','BA  ','LA  ','CE  ','PR  ','ND  ',            &
                  'PM  ','SM  ','EU  ','GD  ','TB  ','DY  ',            &
                  'HO  ','ER  ','TM  ','YB  ','LU  ','HF  ',            &
                  'TA  ','W   ','RE  ','OS  ','IR  ','PT  ',            &
                  'AU  ','HG  ','TL  ','PB  ','BI  ','PO  ',            &
                  'AT  ','RN  ','FR  ','RA  ','AC  ','TH  ',            &
                  'PA  ','U   ','NP  ','PU  ','AM  ','CM  ',            &
                  'BK  ','CF  ','ES  ','FM  ','MD  ','NO  ',            &
                  'LR  ','RF  ','X   ','BQ  '/
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: GRAD
!-----------------------------------------------------------------------
      IUSER(NIU23) = IUSER(NIU23) + 1                  ! ITCG = ITCG + 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -            
!     Update coordinates of shells if use libreta library for ERIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(ILIBRETA==1)CALL UPDCOOSHELL(IUSER(NIU7),IUSER(NIU13),Cxyz,    &
                                      IUSER(NIU5))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(GRAD(NV))
      CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),IUSER(NIU5),    &
                    IUSER(NIU6),IUSER(NIU24),IUSER(NIU7),IUSER(NIU8),   &
                    USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),            &
                    IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),IUSER(NIU14),&
                    IUSER(NIU15),IUSER(NIU16),IUSER(NIU17),IUSER(NIU18),&
                    IUSER(NIU19),IUSER(NIU20),IUSER(NIU21),USER(NU3),   &
                    USER(NU4),USER(NU5),USER(NU6),USER(NU7),USER(NU8),  &
                    USER(NU9),USER(NU10),USER(NU11),USER(NU12),GRAD,    &
                    IUSER(NIU22),USER(NU13),0,0)
      DEALLOCATE(GRAD)
!
      ENERGIA = EELEC + EN     
      WRITE(6,'(8X,I3,16X,F20.10)')NF,ENERGIA
      if(IUSER(NIU23)<=200)GEOMENERGY(nf) = ENERGIA
!
      if(MOLDEN==1)then 
       nat = IUSER(NIU5)
       write(18,'(I6)')nat
       write(18,'(I6,F20.10)')NF,ENERGIA
       do jat=1,nat
        j = (jat-1)*3
        IZNUC = INT(USER(NU1+jat-1))+IZCORE(jat)                 
        write(18,'(1X,A4,3F15.4)')ATMLAB(IZNUC),                        &
        Cxyz(1+j)*BOHR,Cxyz(2+j)*BOHR,Cxyz(3+j)*BOHR
       enddo
      endif
!-----------------------------------------------------------------------
      RETURN
      END
      
! CALCOPTG      
      SUBROUTINE CALCOPTG(NV,Cxyz,NF,GRAD,IUSER,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/USELIBRETA/ILIBRETA 

      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9,  &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,     &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,     &
                         NIU24,NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,NU7,     &
                         NU8,NU9,NU10,NU11,NU12,NU13,NULAST
      INTEGER,DIMENSION(NIULAST)         :: IUSER 
      DOUBLE PRECISION,DIMENSION(NV)     :: Cxyz
      DOUBLE PRECISION,DIMENSION(NV)     :: GRAD      
      DOUBLE PRECISION,DIMENSION(NULAST) :: USER
!-----------------------------------------------------------------------
!     Avoiding warnings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NF = NF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     Update coordinates of shells if use libreta library for ERIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(ILIBRETA==1)CALL UPDCOOSHELL(IUSER(NIU7),IUSER(NIU13),Cxyz,    &
                                      IUSER(NIU5))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),IUSER(NIU5),    &
                    IUSER(NIU6),IUSER(NIU24),IUSER(NIU7),IUSER(NIU8),   &
                    USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),            &
                    IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),IUSER(NIU14),&
                    IUSER(NIU15),IUSER(NIU16),IUSER(NIU17),IUSER(NIU18),&
                    IUSER(NIU19),IUSER(NIU20),IUSER(NIU21),USER(NU3),   &
                    USER(NU4),USER(NU5),USER(NU6),USER(NU7),USER(NU8),  &
                    USER(NU9),USER(NU10),USER(NU11),USER(NU12),GRAD,    &                    
                    IUSER(NIU22),USER(NU13),0,0)
!-----------------------------------------------------------------------
      RETURN
      END
      
! OPTCGNAG
      SUBROUTINE OPTCGNAG(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,    &
                          NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,   &
                          KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,     &
                          ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,IRUNTYP,   &
                          GRADIENT,DIPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      PARAMETER (BOHR = 0.52917724924D+00) 
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN,ZMASS
      INTEGER,DIMENSION(NAT):: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMI):: C1,C2,EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT):: Cxyz
      DOUBLE PRECISION,DIMENSION(3*NAT),INTENT(OUT):: GRADIENT
      INTEGER,INTENT(IN):: IRUNTYP
      DOUBLE PRECISION,DIMENSION(3):: DIPS
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9,  &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,     &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,     &
                         NIU24,NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,NU7,     &
                         NU8,NU9,NU10,NU11,NU12,NU13,NULAST
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/GEOCONV/GEOMENERGY(200)
      COMMON/USELIBRETA/ILIBRETA      
!                         
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IUSER,IWORK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: USER,WORK,GRADS
      EXTERNAL ENERGYFUN
!-----------------------------------------------------------------------
!     Define Pointers of the user arrays
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL POINTERSOPT(NAT,NSHELL,NPRIMI)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Transfer working arrays to IUSER and USER
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(IUSER(NIULAST),USER(NULAST))
!
      IUSER(NIU1)          = NINTEG
      IUSER(NIU2)          = IDONTW  
      IUSER(NIU3)          = 3
      IUSER(NIU4)          = IEMOM  
      IUSER(NIU5)          = NAT
      IUSER(NIU6)          = NBF
      IUSER(NIU7)          = NSHELL
      IUSER(NIU8)          = NPRIMI
      IUSER(NIU9 :NIU10-1) = IAN
      IUSER(NIU10:NIU11-1) = IMIN
      IUSER(NIU11:NIU12-1) = IMAX
      IUSER(NIU12:NIU13-1) = KSTART 
      IUSER(NIU13:NIU14-1) = KATOM
      IUSER(NIU14:NIU15-1) = KTYPE
      IUSER(NIU15:NIU16-1) = KLOC
      IUSER(NIU16:NIU17-1) = INTYP
      IUSER(NIU17:NIU18-1) = KNG
      IUSER(NIU18:NIU19-1) = KMIN
      IUSER(NIU19:NIU20-1) = KMAX
      IUSER(NIU20:NIU21-1) = ISH
      IUSER(NIU21:NIU22-1) = ITYP
      IUSER(NIU22)         = IRUNTYP
      IUSER(NIU23)         = 0             ! ITCG
      IUSER(NIU24)         = NBFaux      
!
      USER(NU1 :NU2-1)     = ZAN
      USER(NU2 :NU3-1)     = ZMASS 
      USER(NU3 :NU4-1)     = C1
      USER(NU4 :NU5-1)     = C2
      USER(NU5 :NU6-1)     = EX
      USER(NU6 :NU7-1)     = CS
      USER(NU7 :NU8-1)     = CP
      USER(NU8 :NU9-1)     = CD
      USER(NU9 :NU10-1)    = CF
      USER(NU10:NU11-1)    = CG
      USER(NU11:NU12-1)    = CH
      USER(NU12:NU13-1)    = CI
      USER(NU13:NULAST)    = DIPS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Initial Coordinates on File CGGRAD (Unit=11)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(11,'(/,72(1H*),/)')
      WRITE(11,*)'Initial Coordinates (Angs) for Geometry Optimization'
      WRITE(11,*)
      DO I=1,NAT
       WRITE(11,'(I5,3F15.4)')  &
             I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
      ENDDO
      WRITE(11,'(/,72(1H*),/)')      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Minimization of the total energy with respect to Cxyz
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,'(/,2X,A21,/1X,23(1H=))')'GEOMETRY OPTIMIZATION'
      WRITE(6,'(/,1X,A23,7X,A19,/)')'Call in CG Optimization',          &
                                    'Total Energy (a.u.)'
      IF(MOLDEN==1)WRITE(18,'(A18)')'[GEOMETRIES] (XYZ)'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Send output of E04DGF to CGM file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL X04ABF(1,2)                                           !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     Est. opt. function val.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     CALL E04DKF ('Es = ....')                                   !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Function Precision (machine precision**0.9)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Function Precision = 1.0D-6')                !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Linesearch Tolerance (0<r<1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Linesearch Tolerance = 0.1')                 !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Maximum Step Length (Default = 10**20))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Maximum Step Length  = 0.01')                !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Optimality Tolerance (deafult = relative precision**0.8)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Optimality Tolerance = 1.0D-5')              !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Verify Level (-1 = No checks, 0 = cheap test, 1 = 0 + gradients)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Verify Level = -1')                          !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Print Level (0 = No output, 1 = The final solution only)
!                 (5 = One line of summary output for each iteration)
!                 (10 = The final solution and one line for each iter.)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Print Level = 0')                            !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calling to NAG Library for using the CG method
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NV = 3*NAT
      IFAIL = 1
      EMINIMA = 1.0d6
      ALLOCATE(IWORK(NV+1),WORK(13*NV),GRADS(NV))
!      CALL E04DGF(NV,ENERGYFUN,ITER_E04DGF,ENERGIA,GRADS,Cxyz,          &  !nag
!                  IWORK,WORK,IUSER,USER,IFAIL)                             !nag
      IF(EMINIMA>ENERGIA)EMINIMA=ENERGIA
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute energy at solution to print data on the output file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(IFAIL==0)then
       WRITE(6,'(/1X,A24,2X,F20.10,/)')'Optimal solution found !',EMINIMA
      else if(IFAIL==6)then       
       WRITE(6,'(/1X,A24,2X,F20.10,/)')'Final objective value = ',EMINIMA      
      end if
      WRITE(6,*)'New Coordinates after Conjugate Gradient Opt (Bohr)'
      WRITE(6,'(1X,51(1H-),/)')
      DO I=1,NAT
       WRITE(6,'(I5,3F15.4)')I,Cxyz(1,I),Cxyz(2,I),Cxyz(3,I)
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     Update coordinates of shells if use libreta library for ERIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(ILIBRETA==1)CALL UPDCOOSHELL(IUSER(NIU7),IUSER(NIU13),Cxyz,    &
                                      IUSER(NIU5))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),IUSER(NIU5),    &
                    IUSER(NIU6),IUSER(NIU24),IUSER(NIU7),IUSER(NIU8),   &
                    USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),            &
                    IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),IUSER(NIU14),&
                    IUSER(NIU15),IUSER(NIU16),IUSER(NIU17),IUSER(NIU18),&
                    IUSER(NIU19),IUSER(NIU20),IUSER(NIU21),USER(NU3),   &
                    USER(NU4),USER(NU5),USER(NU6),USER(NU7),USER(NU8),  &
                    USER(NU9),USER(NU10),USER(NU11),USER(NU12),GRADS,   &
                    IUSER(NIU22),USER(NU13),0,1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Final Coordinates on File CGGRAD (Unit=100)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(11,'(/,72(1H*),/)')
      WRITE(11,*)'New Coordinates after Conjugate Gradient Opt (Angs)'
      WRITE(11,*)
      DO I=1,NAT
       WRITE(11,'(I5,3F15.4)')                                          &
            I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
      ENDDO
      CALL NUCDIST(NV,NAT,Cxyz)
      WRITE(11,'(/,72(1H*),/)')      
      GRADIENT = GRADS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write geometry energies in [GEOCONV] section of the molden file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      IF(MOLDEN==1)WRITE(18,'(A9,/,A6)')'[GEOCONV]','energy'
      if(IUSER(NIU23)<=200)then
       DO I=1,IUSER(NIU23)
        WRITE(18,'(F20.10)')GEOMENERGY(I)
       END DO
      else
       write(6,*)'Sorry, No. of calls > 200, [GEOCONV] is not possible'
      end if
!-----------------------------------------------------------------------
      DEALLOCATE(IUSER,USER,GRADS,IWORK,WORK)
      RETURN
      END
      
! ENERGYFUN
      SUBROUTINE ENERGYFUN(MODE,NV,Cxyz,ENERGIA,GRADS,NSTATE,IUSER,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9,  &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,     &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,     &
                         NIU24,NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,NU7,     &
                         NU8,NU9,NU10,NU11,NU12,NU13,NULAST
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC_OLD,EELEC,DIF_EELEC,EELEC_MIN
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/GEOCONV/GEOMENERGY(200)
      COMMON/USELIBRETA/ILIBRETA
      !JFHLewYee: Changed NATOMS allowed dimension from 100 to 1000
      COMMON/ECP2/CLP(4004),ZLP(4004),NLP(4004),KFRST(1001,6),          &
                  KLAST(1001,6),LMAX(1001),LPSKIP(1001),IZCORE(1001)
!
      INTEGER,DIMENSION(*) :: IUSER
      DOUBLE PRECISION,DIMENSION(*) :: USER
      DOUBLE PRECISION,DIMENSION(NV) :: Cxyz,GRADS
      PARAMETER (BOHR = 0.52917724924D+00)
      CHARACTER*4,DIMENSION(106) :: ATMLAB
      DATA ATMLAB/'H   ','HE  ','LI  ','BE  ','B   ','C   ',            &
                  'N   ','O   ','F   ','NE  ','NA  ','MG  ',            &
                  'AL  ','SI  ','P   ','S   ','CL  ','AR  ',            &
                  'K   ','CA  ','SC  ','TI  ','V   ','CR  ',            &
                  'MN  ','FE  ','CO  ','NI  ','CU  ','ZN  ',            &
                  'GA  ','GE  ','AS  ','SE  ','BR  ','KR  ',            &
                  'RB  ','SR  ','Y   ','ZR  ','NB  ','MO  ',            &
                  'TC  ','RU  ','RH  ','PD  ','AG  ','CD  ',            &
                  'IN  ','SN  ','SB  ','TE  ','I   ','XE  ',            &
                  'CS  ','BA  ','LA  ','CE  ','PR  ','ND  ',            &
                  'PM  ','SM  ','EU  ','GD  ','TB  ','DY  ',            &
                  'HO  ','ER  ','TM  ','YB  ','LU  ','HF  ',            &
                  'TA  ','W   ','RE  ','OS  ','IR  ','PT  ',            &
                  'AU  ','HG  ','TL  ','PB  ','BI  ','PO  ',            &
                  'AT  ','RN  ','FR  ','RA  ','AC  ','TH  ',            &
                  'PA  ','U   ','NP  ','PU  ','AM  ','CM  ',            &
                  'BK  ','CF  ','ES  ','FM  ','MD  ','NO  ',            &
                  'LR  ','RF  ','X   ','BQ  '/
!-----------------------------------------------------------------------
      NSTATE = 1
      MODE = 1
      IUSER(NIU23) = IUSER(NIU23) + 1                  ! ITCG = ITCG + 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     Update coordinates of shells if use libreta library for ERIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(ILIBRETA==1)CALL UPDCOOSHELL(IUSER(NIU7),IUSER(NIU13),Cxyz,    &
                                      IUSER(NIU5))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),IUSER(NIU5),    &
                    IUSER(NIU6),IUSER(NIU24),IUSER(NIU7),IUSER(NIU8),   &
                    USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),            &
                    IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),IUSER(NIU14),&
                    IUSER(NIU15),IUSER(NIU16),IUSER(NIU17),IUSER(NIU18),&
                    IUSER(NIU19),IUSER(NIU20),IUSER(NIU21),USER(NU3),   &
                    USER(NU4),USER(NU5),USER(NU6),USER(NU7),USER(NU8),  &
                    USER(NU9),USER(NU10),USER(NU11),USER(NU12),GRADS,   &
                    IUSER(NIU22),USER(NU13),0,0)
      ENERGIA = EELEC + EN     
      WRITE(6,'(8X,I3,16X,F20.10)')IUSER(NIU23),ENERGIA
      if(IUSER(NIU23)<=200)GEOMENERGY(IUSER(NIU23)) = ENERGIA
!
      if(MOLDEN==1)then 
       nat = IUSER(NIU5)
       write(18,'(I6)')nat
       write(18,'(I6,F20.10)')IUSER(NIU23),ENERGIA
       do jat=1,nat
        j = (jat-1)*3
        IZNUC = INT(USER(NU1+jat-1))+IZCORE(jat)                 
        write(18,'(1X,A4,3F15.4)')ATMLAB(IZNUC),                        &
        Cxyz(1+j)*BOHR,Cxyz(2+j)*BOHR,Cxyz(3+j)*BOHR
       enddo
      endif
!-----------------------------------------------------------------------
      RETURN
      END
      
! OPTLBFGS
      SUBROUTINE OPTLBFGS(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,    &
                          NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,   &
                          KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,      &
                          C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,IRUNTYP,        &
                          GRADIENT,DIPS,NPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC_OLD,EELEC,DIF_EELEC,EELEC_MIN
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
      COMMON/USELIBRETA/ILIBRETA                  
!      
      PARAMETER (BOHR = 0.52917724924D+00)
      INTEGER,INTENT(IN)::NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux
      INTEGER,INTENT(IN)::NSHELL,NPRIMI,IRUNTYP,NPRINT
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
      INTEGER,DIMENSION(NAT):: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMI):: C1,C2,EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT):: Cxyz
      DOUBLE PRECISION,DIMENSION(3*NAT),INTENT(OUT):: GRADIENT
      DOUBLE PRECISION,DIMENSION(3*NAT):: GRADS
      DOUBLE PRECISION,DIMENSION(3):: DIPS
      INTEGER,PARAMETER::MSAVE=7
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::W
      DOUBLE PRECISION X(3,NAT),G(3*NAT),DIAG(3*NAT)
      DOUBLE PRECISION::F,EPS,XTOL,GTOL,STPMIN,STPMAX,EELEC_MIN_LBFGS
      INTEGER::IFLAG,ICALL,N,M,MP,LP,NWORK
      INTEGER,DIMENSION(2)::IPRINT
      LOGICAL::DIAGCO
!     The driver for LBFGS must always declare LB2 as EXTERNAL
      EXTERNAL LB2
!-----------------------------------------------------------------------
      NWORK=3*NAT*(2*MSAVE +1)+2*MSAVE
      ALLOCATE(W(NWORK))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Initial Coordinates on File CGGRAD (Unit=11)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(11,'(/,72(1H*),/)')
      WRITE(11,*)'Initial Coordinates (Angs) for Geometry Optimization'
      WRITE(11,*)
      DO I=1,NAT
       WRITE(11,'(I5,3F15.4)')  &
            I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
      ENDDO
      WRITE(11,'(/,72(1H*),/)')      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Minimization of the total energy with respect to Cxyz
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,'(/,2X,A21,/1X,23(1H=))')'GEOMETRY OPTIMIZATION'
      WRITE(6,'(/,1X,A23,7X,A19,/)')'Call in CG Optimization',         &
                                    'Total Energy (a.u.)'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calling to LBFGS SUBROUTINE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     CHECK COMMON /LB3, MP SETS UNIT WHERE PRINTING OPTIMIZATION INFO, 
!     AND LP WHERE PRINTING INFO ABOUT ERRORS
      N=3*NAT ! NUMBER OF VARIABLES
      M=5     ! 0 <= M <= 7
      IPRINT(1)= 0
!     IPRINT(1)< 0 : no output is generated,
!     IPRINT(1)= 0 : output only at first and last iteration,
!     IPRINT(1)> 0 : output every IPRINT(1) iterations.      
      IPRINT(2)= 0
      IF(NPRINT==2) IPRINT(2)= 1
!     IPRINT(2)= 0 : iteration count, number of function 
!                     evaluations, function value, norm of the
!                     gradient, and steplength,
!     IPRINT(2)= 1 : same as IPRINT(2)=0, plus vector of
!                     variables and  gradient vector at the
!                     initial point,
!     IPRINT(2)= 2 : same as IPRINT(2)=1, plus vector of
!                     variables,
!     IPRINT(2)= 3 : same as IPRINT(2)=2, plus gradient vector.      
!
!     We do not wish to provide the diagonal matrices Hk0, and 
!     therefore set DIAGCO to FALSE.
      DIAGCO= .FALSE.
      EPS= 1.0D-10
      XTOL= 1.11D-16
      ICALL=0
      IFLAG=0
      X = Cxyz ! X(N) IS THE INITIAL ESTIMATE OF THE SOLUTION VECTOR
      DO
       ICALL=ICALL + 1
!      Update coordinates of shells if use libreta library for ERIs
       if(ILIBRETA==1)CALL UPDCOOSHELL(NSHELL,KATOM,Cxyz,NAT)
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
                     ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,    &
                     INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF, &
                     CG,CH,CI,GRADS,IRUNTYP,DIPS,0,0)
       WRITE(6,'(8X,I3,16X,F20.10)')ICALL,EELEC+EN
                     
!      F CONTAINS THE VALUE OF THE FUNCTION AT THE POINT X
!      G CONTAINS THE COMPONENTS OF GRADIENT AT X
       IF(ICALL==1) EELEC_MIN_LBFGS = EELEC
       IF(ICALL/=1.AND.EELEC<EELEC_MIN_LBFGS) EELEC_MIN_LBFGS = EELEC
       F = EELEC + EN
       G = GRADS
       CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
       Cxyz = X
       IF(IFLAG.LE.0) EXIT
!      We allow at most 1000 evaluations of F and G
       IF(ICALL.GT.1000) EXIT
      ENDDO
      DEALLOCATE(W)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -             
!     FINAL CALL TO COMPUTE ENERGY AND GRADS AT EQUIL. GEOMETRY
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
      IF(IFLAG==0.AND.EELEC<=EELEC_MIN_LBFGS) THEN
       WRITE(6,'(/1X,A24,/)')'Optimal solution found !'
      ELSEIF(IFLAG.LE.0.OR.EELEC>EELEC_MIN_LBFGS) THEN
       WRITE( 6,'(A38,/)')' !!! Cannot find optimal solution !!! '
       WRITE(11,'(A38,/)')' !!! Cannot find optimal solution !!! '
      ENDIF
      WRITE(6,*)'New Coordinates after the LBFGS Method (Bohr)'
      WRITE(6,'(45(1H-),/)')
      DO I=1,NAT
       WRITE(6,'(I5,3F20.10)')I,Cxyz(1,I),Cxyz(2,I),Cxyz(3,I)
      ENDDO
!     Update coordinates of shells if use libreta library for ERIs
      if(ILIBRETA==1)CALL UPDCOOSHELL(NSHELL,KATOM,Cxyz,NAT)
      CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,   &
                    ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,     &
                    INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,  &
                    CG,CH,CI,GRADS,IRUNTYP,DIPS,0,1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Final Coordinates on File CGGRAD (Unit=100)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(11,*)
      WRITE(11,*)'New Coordinates after the LBFGS Method (Angs)'
      WRITE(11,*)
      DO I=1,NAT
       WRITE(11,'(I5,3F20.10)')                                         &
            I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
      ENDDO
      CALL NUCDIST(N,NAT,Cxyz)
      GRADIENT = GRADS
!-----------------------------------------------------------------------
      RETURN
      END
