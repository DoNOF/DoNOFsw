!======================================================================!
!                                                                      !
!                F R A G M E N T   S U B R O U T I N E S               !
!                                                                      !
!                  ( Theor. Chem. Acc. 134, 151, 2015 )                !
!                                                                      !
!======================================================================!
!                                                                      !
!   FragOrbOpt: Optimize Orbitals of the Core-Fragment                 !
!   COEFW_F: COEFN = COEF*W where COEFN(N,Nf), COEF(N,N) and W(Nf,Nf)  !
!            Note. the sum index k runs within fragment kf = INDf(k)   !
!   ENERGYFrc: Calculate the electronic energy and Lagrange Multipliers!
!   FFJMNFrc: Calculate the gen. Fock matrices Fj(m,n) and Frj(m,n)    !
!   ELGF: Calculate the Lagrange Multipliers of the fragment 'f'       !              
!   PCONVE_F: Check for the symmetry of Lagrangian of the fragment 'f' !
!   FORMFMIUG_F: Form & decrease gen-Fock using a scaling factor of 'f'!
!   EELECTRr: Trace( Ct*RO*H*C + Ct*G )                                !
!   EELECTRr_EFIELD: Trace [ Ct*RO*(Ei*ADIPi)*C ] - Ei*DIPN(i)         !
!                                                                      !
!======================================================================!

! FragOrbOpt
      SUBROUTINE FragOrbOpt(ITCALL,ITLIM,AHCORE,IJKL,XIJKL,QD,COEF,RO,  &
                           CJ12,CK12,ELAG,FMIUG0,DIPN,ADIPx,ADIPy,ADIPz,&
                           ILOOP,IRUNTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL DIIS,PERDIIS,CONVGDELAG,RESTART,SCALING
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPSCALING/SCALING,NZEROS,NZEROSm,NZEROSr,ITZITER
      COMMON/CONVERGENCE/DUMEL,PCONV,CONVGDELAG
      COMMON/INPNOF_DIIS/DIIS,PERDIIS,NDIIS,NTHDIIS,THDIIS
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_THRESH/THRESHL,THRESHE,THRESHEC,THRESHEN
      COMMON/INPNOF_COEFOPT/MAXLOOP
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_FRAG/NO1f,NBFf,NBF5f
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/CONVERGESUM/SUMDIF,SUMDIF_OLD
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
!
      INTEGER,DIMENSION(NSTORE)::IJKL
      DOUBLE PRECISION,DIMENSION(NSTORE)::XIJKL
      DOUBLE PRECISION,DIMENSION(NBF)::FMIUG0
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AHCORE,COEF,ELAG
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(3)::DIPN
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
!
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::EVA,TEMP,CFM
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::FMIUG,W,COEFNEW,BFM
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:)::FK
!
      INTEGER,ALLOCATABLE,DIMENSION(:)::INO1,INDOC,INCWO,INAC,INO0,INDf
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::WFr
!-----------------------------------------------------------------------
!     Read fragment information on the FILE 10 (FRAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      REWIND(10)
      READ(10,'(I6)')NFRAG
      ALLOCATE (INO1(NFRAG),INDOC(NFRAG),INCWO(NFRAG))
      ALLOCATE (INAC(NFRAG),INO0(NFRAG))
      DO i=1,NFRAG
       READ(10,'(6I6)')INO1(i),INDOC(i),INSOC,INCWO(i),INAC(i),INO0(i)
      ENDDO 
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
!
!     NO1f:  Number of inactive doubly occupied orbitals (OCC=1)
!     NDOCf: Number of strongly occupied MOs
!     NCWOf: Number of coupled weakly occupied MOs per strongly occupied
!     NO0f: Empty orbitals  (OCC=0)
!
!     NCOf:  Number of HF occupied MOs (OCC=1)
!     NVIRf: Number of HF virtual  MOs (OCC=0)
!     NCWOf*NDOCf: Active orbitals in the virtual subspace
!     NACf: Dimension of the active natural orbital subspace
!
!           NCOf      |       NVIRf          = NBFf
!       NO1f + NDOCf  |  NCWOf*NDOCf + NO0f  = NBFf
!            |      NACf             |
!
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
      NO1f=INO1(1)
      NDOCf=INDOC(1)
      NCWOf=INCWO(1)
      NO0f=INO0(1)

      NCOf = NO1f + NDOCf
      NVIRf= NCWOf*NDOCf + NO0f
      NACf = NDOCf * ( 1 + NCWOf )
      NBFf = NO1f + NACf + NO0f

      NBF5f = NO1f + NACf
      IF(NBF5f>NBFf)NBF5f = NBFf
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Index correspondence
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (INDf(NBFf))
      DO i=1,NO1f 
       INDf(i) = i                       ! 1,NO1f
      ENDDO
      DO i=1,NACf                
       INDf(NO1f+i) = NCO-NDOCf+i        ! NCO-NDOCf+1,NCO+NCWOf*NDOCf
      ENDDO
      DO i=1,NO0f                     
       INDf(NO1f+NACf+i) = NBF-NO0f+i    ! NBF-NO0f+1,NBF
      ENDDO
!-----------------------------------------------------------------------
!                          Orbital Optimization
!-----------------------------------------------------------------------
!     First Call to FragOrbOpt
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ITCALL==1)THEN
       IF(IRUNTYP==3.or.IRUNTYP==4)NZEROS=NZEROSr      
       WRITE(6,1)
       EELEC_OLD=EELEC
       SUMDIF_OLD=0.0d0
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate Initial Electronic Energy (EELEC) and Lagrangian
!     Note: This Energy is equal to EELEC from MOLOCUPrc
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (WFr(NSQ,NBF5f))
      CALL ENERGYFrc(AHCORE,IJKL,XIJKL,QD,COEF,RO,CJ12,CK12,ELAG,       &
                     DIPN,ADIPx,ADIPy,ADIPz,INDf,WFr,0)
      CALL PCONVE_F(INDf,ELAG,DUMEL,MAXI,MAXJ,SUMDIF)
      DIF_EELEC=EELEC-EELEC_OLD
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Output of the first evaluated energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,2)ITCALL,EELEC,EELEC+EN,DIF_EELEC,DUMEL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check ITCALL and SUMDIF for changing the scaling factor
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(SCALING.and.ITCALL>2.and.ITCALL>ITLIM.and.                     &
         SUMDIF>SUMDIF_OLD)THEN
       NZEROS=NZEROS+1
       ITLIM=ITCALL+ITZITER
!      NZEROS>NZEROSm -> Restart from here with NZEROSr
       IF(NZEROS>NZEROSm)THEN
        NZEROS=NZEROSr
       ENDIF
      ENDIF
      SUMDIF_OLD=SUMDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check for the Symmetry of Lagrangian (ELAG(IQJQ)-ELAG(JQIQ))
!     for ITCALL>1 (do not stop in the first call)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      PCONV=ABS(DIF_EELEC)
      IF( ITCALL>1 .and. DUMEL<THRESHL .and. PCONV < THRESHE )THEN
       CONVGDELAG = .TRUE.
       RETURN
      ENDIF
!-----------------------------------------------------------------------
!                       START SCF-ITERATION CYCLE
!-----------------------------------------------------------------------
      ALLOCATE (FMIUG(NBFf,NBFf),W(NBFf,NBFf),EVA(NBFf),TEMP(NBFf))
      ALLOCATE (COEFNEW(NBF,NBFf),FK(MAXLOOP,NBFf,NBFf))
      ALLOCATE (BFM(MAXLOOP+1,MAXLOOP+1),CFM(MAXLOOP+1))
!
      IF(ITCALL==1.and.INPUTFMIUG==0)THEN
       MAXLP=1
      ELSE
       MAXLP=MAXLOOP
      ENDIF
!-----------------------------------------------------------------------
      ILOOP=0
      IDIIS=0
      DO LOOP=1,MAXLP
       ILOOP = ILOOP+1
       EELEC0 = EELEC
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Generalized Fock Matrix (FMIUG)
!
!      Convergent technique:
!
!      1) SCALING:   Decrease FMIUG using a scaling factor.
!                    The scaling factor varies until the number of
!                    ZEROS (.000##) is equal for all elements Fij
!      3) DIIS       Direct Inversion in the Iterative Subspace 
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(SCALING)CALL FORMFMIUG_F(INDf,FMIUG,ELAG,FMIUG0,ITCALL)
       IF(DIIS.and.DUMEL<THDIIS)THEN
        CALL FFMIUG_DIIS(NBFf,FMIUG,CFM,BFM,FK,IDIIS)
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      DIAGONALIZE SQUARE MATRIX (FMIUG) FOR REAL SYMMETRIC CASE
!      W - EIGENVECTORS, EVA - EIGENVALUES IN ALGEBRAIC DESCENDING ORDER
!      HOUSEHOLDER METHOD
!      NOTE: ONLY LOWER TRIANGLE IS USED + THIS IS DESTROYED !!!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL DIAG(NBFf,FMIUG,W,EVA,TEMP)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Move EVA -> FMIUG0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO I=1,NBFf
        IQ=INDf(i)
        FMIUG0(IQ)=EVA(I)
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      New Coefficients (COEFNEW=COEF*W), Move COEFNEW -> COEF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL COEFW_F(INDf,NBF,NBFf,COEFNEW,COEF,W)
       DO j=1,NBFf
        JQ=INDf(j)
        DO i=1,NBF
         COEF(i,JQ) = COEFNEW(i,j)
        ENDDO
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Calculate Dj(miu,niu), Jj(miu,niu), Kj(miu,niu) (j=1,nbf)
!      Keep Matrices Jj, Kj in WJj, WKj
!      Form F Matrix and keep it in WF
!      Compute G, Lagrangian Multipliers (ELAG) and one-energies (E)
!      Calculate Electronic Energy (EELEC)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL ENERGYFrc(AHCORE,IJKL,XIJKL,QD,COEF,RO,CJ12,CK12,ELAG,      &
                      DIPN,ADIPx,ADIPy,ADIPz,INDf,WFr,1)
       DELE = EELEC - EELEC0
!      Intermediate Output (Nprint=2)
       IF(NPRINT==2)WRITE(6,3)LOOP,EELEC,EELEC+EN,DELE
       CALL PCONVE_F(INDf,ELAG,DUMEL,MAXI,MAXJ,SUMDIF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Check for energy convergent solution
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF( DABS(DELE) < THRESHEC )THEN
        DEALLOCATE (COEFNEW,FMIUG,W,EVA,TEMP)
        DIF_EELEC=EELEC-EELEC_OLD
        EELEC_OLD=EELEC
        RETURN
       ENDIF
!------------------------------------------------------------------------
!                       LOOP-END OF SCF-ITERATION
!-----------------------------------------------------------------------
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Itermediate Output of the external iteration
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DIF_EELEC=EELEC-EELEC_OLD
      EELEC_OLD=EELEC
!-----------------------------------------------------------------------
!     FORMAT STATEMENTS
!-----------------------------------------------------------------------
    1 FORMAT(/2X,'RESULTS OF OCCUPATION-COEFFICIENT S.C.F. PROCEDURE'   &
            ,/1X,'===================================================', &
        //2X,'Iter',5X,'Electronic Energy',6X,'Total Energy',           &
          3X,'Energy Convergency',4X,'Max Mul-Lag Diff',/)
    2 FORMAT(I5,'.',1X,F20.10,1X,F19.10,2X,F15.10,8X,F11.6)
    3 FORMAT(2X,I3,'.',3X,F17.8,4X,F15.8,6X,F11.6)
!-----------------------------------------------------------------------
      DEALLOCATE (INO1,INDOC,INCWO,INAC,INO0,INDf,WFr)
      RETURN
      END

! COEFW_F
      SUBROUTINE COEFW_F(INDf,N,Nf,COEFN,COEF,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(Nf)::INDf
      DOUBLE PRECISION,DIMENSION(N,Nf)::COEFN
      DOUBLE PRECISION,DIMENSION(N,N)::COEF
      DOUBLE PRECISION,DIMENSION(Nf,Nf)::W
!-----------------------------------------------------------------------
!     COEFN = COEF*W
!-----------------------------------------------------------------------
      DO i=1,N
       DO j=1,Nf
        COEFN(i,j)=0.0d0
        do k=1,Nf
         kf = INDf(k)
         COEFN(i,j) = COEFN(i,j) + COEF(i,kf)*W(k,j)
        enddo
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! ENERGYFrc
      SUBROUTINE ENERGYFrc(AHCORE,IJKL,XIJKL,QD,COEF,RO,CJ12,CK12,      &
                            ELAG,DIPN,ADIPx,ADIPy,ADIPz,INDf,WFr,IE)
!-----------------------------------------------------------------------
!     Calculate the electronic energy and Lagrange Multipliers
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL EFIELDL
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INP_EFIELDL/EX,EY,EZ,EFIELDL
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_FRAG/NO1f,NBFf,NBF5f
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
!
      INTEGER,DIMENSION(NBFf)::INDf
      INTEGER,DIMENSION(NSTORE)::IJKL
      DOUBLE PRECISION,DIMENSION(NSTORE)::XIJKL
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AHCORE,COEF,ELAG
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(3)::DIPN
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5f)::WFr
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::WJj,WKj,WF,G,AUX1
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::AUX2
!-----------------------------------------------------------------------
!     Calculate Dj: QD(j,miu,niu), Jj(miu,niu), Kj(miu,niu) (j=1,NBF5) 
!-----------------------------------------------------------------------
      ALLOCATE (WJj(NSQ,NBF5),WKj(NSQ,NBF5),AUX1(NBF,NBF),AUX2(NSQ))
      DO j=1,NBF5
       CALL DENMATj(j,AUX1,COEF,NBF)
       QD(j,1:NBF,1:NBF) = AUX1(1:NBF,1:NBF)
       CALL HSTARJ(AUX2,AUX1,IJKL,XIJKL)
       WJj(1:NSQ,j) = AUX2(1:NSQ)
       CALL HSTARK(AUX2,AUX1,IJKL,XIJKL)
       WKj(1:NSQ,j) = AUX2(1:NSQ)
      ENDDO
      DEALLOCATE (AUX1,AUX2)
!-----------------------------------------------------------------------
!     Form F Matrix and keep it in WF
!     WFr: Constant part of F related to the interactions with fragments
!-----------------------------------------------------------------------
      ALLOCATE (WF(NSQ,NBF5f))
      CALL FFJMNFrc(INDf,RO,CJ12,CK12,AHCORE,WJj,WKj,WF,WFr,            &
                    ADIPx,ADIPy,ADIPz,IE)
      DEALLOCATE (WJj,WKj)
!-----------------------------------------------------------------------
!     Calculate G Matrix
!-----------------------------------------------------------------------
      ALLOCATE (G(NBF,NBF5f))
      DO jf=1,NBF5f
       IQ = INDf(jf)
       do i=1,nbf
        G(i,jf) = FC(i,IQ,WF(1,jf),COEF) + FC(i,IQ,WFr(1,jf),COEF)
       enddo
      ENDDO
      DEALLOCATE (WF)
!-----------------------------------------------------------------------
!     Lagrangian Multipliers (ELAG)
!-----------------------------------------------------------------------
      CALL ELGF(INDf,ELAG,COEF,G)
      DEALLOCATE (G)
!-----------------------------------------------------------------------
!     Calculate ELECTRONIC ENERGY
!-----------------------------------------------------------------------
!     Calculate Trace(Ct*RO*H*C+Ct*G)
      CALL EELECTRr(EELEC,AHCORE,ELAG,COEF,RO)
!     Include Nuclear Dipoles if electric field =/ 0
      IF(EFIELDL)THEN
       CALL EELECTRr_EFIELD(EELEC_EF,COEF,RO,DIPN,ADIPX,ADIPY,ADIPZ)
       EELEC=EELEC+EELEC_EF
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END

! FFJMNFrc
      SUBROUTINE FFJMNFrc(INDf,RO,CJ12,CK12,H,DJ,DK,F,Fr,               &
                          ADIPx,ADIPy,ADIPz,IE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL EFIELDL
      COMMON/INP_EFIELDL/EX,EY,EZ,EFIELDL
      COMMON/INPFILE_FRAG/NO1f,NBFf,NBF5f 
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
!
      INTEGER,DIMENSION(NBFf)::INDf
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::H,ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::DJ,DK
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5f)::F,Fr
!-----------------------------------------------------------------------
!                            Calculate Fj(m,n)
!-----------------------------------------------------------------------
      IF(NO1f>1)THEN
       do i=1,nbf
        do k=1,nbf
         ik=i+(k-1)*nbf
         F(ik,1) = H(i,k) + PRODWCWFik(ik,INDf,DJ,CJ12)                 &
                          - PRODWCWFik(ik,INDf,DK,CK12)                  
         DO jf=NO1f+1,NBF5f                                              
          j = INDf(jf)                                                   
          F(ik,jf) = RO(j) * ( H(i,k) + DJ(ik,j) )                      &
                   + PRODWCWFikq(ik,jf,INDf,DJ,CJ12)                    &
                   - PRODWCWFikq(ik,jf,INDf,DK,CK12)
         ENDDO
        enddo
       enddo
       DO jf=2,NO1f
        F(1:NSQ,jf) = F(1:NSQ,1)
       ENDDO
      ELSE
       DO jf=1,NBF5f
        j = INDf(jf)
        do i=1,nbf
         do k=1,nbf
          ik=i+(k-1)*nbf
          F(ik,jf) = RO(j) * ( H(i,k) + DJ(ik,j) )                      &
                   + PRODWCWFikq(ik,jf,INDf,DJ,CJ12)                    &
                   - PRODWCWFikq(ik,jf,INDf,DK,CK12)
         enddo
        enddo
       ENDDO
      ENDIF
!-----------------------------------------------------------------------
!                          Calculate Frj(m,n)
!-----------------------------------------------------------------------
      IF(IE==0)THEN
       IF(NO1f>1)THEN
        do i=1,nbf
         do k=1,nbf
          ik=i+(k-1)*nbf
          Fr(ik,1) = H(i,k)+PRODWCWij(ik,DJ,CJ12)-PRODWCWij(ik,DK,CK12)
          Fr(ik,1) = Fr(ik,1)-F(ik,1)
          DO jf=NO1f+1,NBF5f
           j = INDf(jf)
           Fr(ik,jf) = RO(j) * ( H(i,k) + DJ(ik,j) )                    &
                     + PRODWCWijq(ik,j,DJ,CJ12)-PRODWCWijq(ik,j,DK,CK12) 
           Fr(ik,jf) = Fr(ik,jf) - F(ik,jf)                              
          ENDDO                                                          
         enddo                                                           
        enddo                                                            
        DO jf=2,NO1f                                                     
         Fr(1:NSQ,jf) = Fr(1:NSQ,1)                                      
        ENDDO                                                            
       ELSE                                                              
        DO jf=1,NBF5f                                                    
         j = INDf(jf)                                                    
         do i=1,nbf                                                      
          do k=1,nbf                                                     
           ik=i+(k-1)*nbf                                                
           Fr(ik,jf) = RO(j) * ( H(i,k) + DJ(ik,j) )                    &
                     + PRODWCWijq(ik,j,DJ,CJ12)-PRODWCWijq(ik,j,DK,CK12)
           Fr(ik,jf) = Fr(ik,jf) - F(ik,jf)
          enddo
         enddo
        ENDDO
       ENDIF
!-----------------------------------------------------------------------
!      Including Electric Field
!      Note: FEikj is constant, so can be keeped in Fr
!-----------------------------------------------------------------------
       IF(EFIELDL)THEN
        DO jf=1,NBF5f
         j = INDf(jf)
         do i=1,nbf
         do k=1,nbf
          ik=i+(k-1)*nbf
          FEikj = ( EX*ADIPx(i,k)+EY*ADIPy(i,k)+EZ*ADIPz(i,k) ) * RO(j)
          Fr(ik,jf) = Fr(ik,jf) + FEikj
         enddo
         enddo
        ENDDO
       ENDIF
!-----------------------------------------------------------------------
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END

! ELGF
      SUBROUTINE ELGF(INDf,ELAG,C,G)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_FRAG/NO1f,NBFf,NBF5f      
      INTEGER,DIMENSION(NBFf)::INDf
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ELAG,C
      DOUBLE PRECISION,DIMENSION(NBF,NBF5f)::G
!-----------------------------------------------------------------------
!     Calculate the Lagrangian Multipliers
!-----------------------------------------------------------------------
      DO i=1,NBFf
       IQ=INDf(i)
       DO j=1,NBF5f
        JQ=INDf(j)
        ELAG(IQ,JQ) = 0.0d0
        do k=1,nbf
         ELAG(IQ,JQ) = ELAG(IQ,JQ) + C(k,IQ)*G(k,j)
        enddo
       ENDDO
       DO j=NBF5f+1,NBFf
        JQ=INDf(j)
        ELAG(IQ,JQ) = 0.0d0
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! PCONVE_F
      SUBROUTINE PCONVE_F(INDf,ELAG,DUM,MAXI,MAXJ,SUMDIF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_FRAG/NO1f,NBFf,NBF5f      
      INTEGER,DIMENSION(NBFf)::INDf
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ELAG
!-----------------------------------------------------------------------
!     MAXIMUM LAGRANGE MULTYPLIER DIFFERENCE [ELAG(IQ,JQ)-ELAG(JQ,IQ)]
!-----------------------------------------------------------------------
      DUM=0.0d0
      SUMDIF=0.0d0
      DO i=1,NBFf
       IQ=INDf(i)
       DO j=1,NBFf
        JQ=INDf(j)
        GCF=DABS(ELAG(IQ,JQ)-ELAG(JQ,IQ))
        SUMDIF=SUMDIF+GCF
        IF(GCF>DUM)THEN
         DUM=GCF
         MAXI=IQ
         MAXJ=JQ
        ENDIF
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! FORMFMIUG_F
      SUBROUTINE FORMFMIUG_F(INDf,FMIUG,ELAG,FMIUG0,ITCALL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL RESTART,SCALING
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPSCALING/SCALING,NZEROS,NZEROSm,NZEROSr,ITZITER
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPFILE_FRAG/NO1f,NBFf,NBF5f
      INTEGER,DIMENSION(NBFf)::INDf
      DOUBLE PRECISION,DIMENSION(NBF)::FMIUG0
      DOUBLE PRECISION,DIMENSION(NBFf,NBFf)::FMIUG
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ELAG
!-----------------------------------------------------------------------
!     Generalized Fock Matrix (FMIUG)
!-----------------------------------------------------------------------
      IF(itcall==1.AND.INPUTFMIUG==0)THEN          ! only for itcall==1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO i=1,NBFf
        IQ=INDf(i)
        DO j=1,i-1
         JQ=INDf(j)
         FMIUG(i,j)=(ELAG(IQ,JQ)+ELAG(JQ,IQ))/2.0        ! Nondiagonal 
         FMIUG(j,i)=FMIUG(i,j)
        ENDDO
        FMIUG(i,i)=ELAG(IQ,IQ)                            ! Diagonal
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       if(SCALING)then
        DO i=1,NBFf
         IQ=INDf(i)
         DO j=1,i-1
          JQ=INDf(j)
          FMIUG(i,j)=ELAG(IQ,JQ)-ELAG(JQ,IQ)              ! Nondiagonal 
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
!         Decrease FMIUG using a scaling factor
!         The scaling factor varies until the number of
!         ZEROS (.000##) is equal for all elements Fij
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
          CALL F01(NZEROS+9,FMIUG(i,j))
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
          FMIUG(j,i)=FMIUG(i,j)                           ! Fji=Fij
         ENDDO
         FMIUG(i,i)=FMIUG0(IQ)                            ! Diagonal
        ENDDO
       else
        DO i=1,NBFf
         IQ=INDf(i)
         DO j=1,i-1
          JQ=INDf(j)
          FMIUG(i,j)=ELAG(IQ,JQ)-ELAG(JQ,IQ)              ! Nondiagonal 
          FMIUG(j,i)=FMIUG(i,j)                           ! Fji=Fij
         ENDDO
         FMIUG(i,i)=FMIUG0(IQ)                            ! Diagonal
        ENDDO
       endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END

! EELECTRr
      SUBROUTINE EELECTRr(ENERGIA,H,ELAG,C,RO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::H,ELAG,C
!-----------------------------------------------------------------------
!     Calculate the electronic energy [ Trace( Ct*RO*H*C + Ct*G ) ]
!-----------------------------------------------------------------------
      ENERGIA = 0.0d0
      DO IQ=1,NB
       ENERGIA = ENERGIA + ELAG(IQ,IQ)
       do i=1,nbf
        ENERGIA = ENERGIA + RO(IQ)*C(i,IQ)*FC(i,IQ,H,C)
       enddo
      ENDDO
!
      if(NSOC>0)then
       if(.not.HighSpin)then
        DO IQ=NB+1,NA
         ENERGIA = ENERGIA + ELAG(IQ,IQ)
         do i=1,nbf
          ENERGIA = ENERGIA + RO(IQ)*C(i,IQ)*FC(i,IQ,H,C)
         enddo
        ENDDO
       else if(HighSpin)then
        DO IQ=NB+1,NA
         ENERGIA = ENERGIA + ELAG(IQ,IQ)
         do i=1,nbf
          ENERGIA = ENERGIA +  0.50d0*C(i,IQ)*FC(i,IQ,H,C)
         enddo
        ENDDO
       end if
      end if       
!      
      DO IQ=NA+1,NBF5
       ENERGIA = ENERGIA + ELAG(IQ,IQ)
       do i=1,nbf
        ENERGIA = ENERGIA + RO(IQ)*C(i,IQ)*FC(i,IQ,H,C)
       enddo
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! EELECTRr_EFIELD
      SUBROUTINE EELECTRr_EFIELD(EELEC_EF,C,RO,DIPN,ADIPx,ADIPy,ADIPz)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EFIELDL
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INP_EFIELDL/EX,EY,EZ,EFIELDL
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DOUBLE PRECISION,DIMENSION(3)::DIPN
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::C,ADIPx,ADIPy,ADIPz
!-----------------------------------------------------------------------
!     Calculate the electronic energy associated to the electric field
!     EELEC_EF = Trace [ Ct*RO*(Ei*ADIPi)*C ] - Ei*DIPN(i)
!-----------------------------------------------------------------------
      EELEC_EF=0.0d0
      DO IQ=1,NBF5
       do i=1,nbf
        EELEC_EF = EELEC_EF + RO(IQ)*C(i,IQ)*( EX*FC(i,IQ,ADIPx,C) +    &
                            EY*FC(i,IQ,ADIPy,C) + EZ*FC(i,IQ,ADIPz,C) )
       enddo
      ENDDO
      EELEC_EF = EELEC_EF - EX*DIPN(1) - EY*DIPN(2) - EZ*DIPN(3)
!-----------------------------------------------------------------------
      RETURN
      END
