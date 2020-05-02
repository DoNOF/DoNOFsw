!======================================================================!
!                                                                      !
!                               Do N O F                               !
!                                                                      !
!       (Donostia Natural Orbital Functional Software Program)         !
!                                                                      !
!                  COPYRIGHT by Mario Piris (2009)                     !
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
!                           Date: April 2020                           !
!                                                                      !
!    Program to compute the ground state properties of a molecule      !
!    in the gas phase using PNOF5 - PNOF7 + perturbation corrections   !
!                                                                      !
!======================================================================!
!                                                                      !
!    2013 Four-index transformation of the electron repulsion          !
!         integrals was parallelized by Eduard Matito                  !
!                                                                      !
!    2015 Electrostatic dipole, quadrupole and octupole moments        !
!         implemented by Ion Mitxelena                                 !
!                                                                      !
!    2017 NOF energy gradients implemented by Ion Mitxelena            !
!                                                                      !
!======================================================================!

!----------------------------------------------------------------------!
!                                                                      !
!   Minimization of the Electronic Energy using Conj. Gradient Method  !
!                                                                      !
!     CGOCUPSUMSLr: Prepare for calling the subroutine SUMSL           !
!     CALCOE: Compute the Occupation Energy                            !
!     CALCOG: Compute the Occupation Energy Gradients                  !
!                                                                      !
!----------------------------------------------------------------------!

! CGOCUPSUMSLr
      SUBROUTINE CGOCUPSUMSLr(NV,GAMMA,USER,ENERGY)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NV) :: GAMMA
      DOUBLE PRECISION,DIMENSION(NUSER) :: USER
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IUSER,IV      
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: D,V     
      EXTERNAL CALCOE,CALCOG
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LIV = 60
      LV = 71+NV*(NV+15)/2
      ALLOCATE( IUSER(1),IV(LIV),D(NV),V(LV) ) 
      IUSER(1) = 1
      IV = 0
      D(1:NV) = 1.0d-1       
      CALL SUMSL(NV,D,GAMMA,CALCOE,CALCOG,IV,LIV,LV,V,IUSER,USER)
      ENERGY = EELEC
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DEALLOCATE(IUSER,IV,D,V)
      RETURN
      END

! CALCOE      
      SUBROUTINE CALCOE(NV,GAMMA,NF,ENERGY,IUSER,USER)
      USE PARCOM 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(1) :: IUSER 
      DOUBLE PRECISION,DIMENSION(NV)   :: GAMMA
      DOUBLE PRECISION,DIMENSION(NUSER):: USER
!-----------------------------------------------------------------------
!     Avoiding warnings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NF = NF
      IUSER(1) = IUSER(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL OCUPACIONr(GAMMA,USER(N1),USER(N2),USER(N3),                &
                            USER(N4),USER(N5),USER(N6),NV)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      if(MSpin==0)then
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      Singlet State (S=0,Ms=0) and Multiplet States (S>0,Ms=0)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       IF(NCWO==1)THEN       ! PNOFi(1): Perfect Pairing (NCWO=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ENERGY = 0.0d0
        do in=1,NO1
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj(in,USER(N2),USER(N9))                      &
                - PRODCWQWj(in,USER(N3),USER(N10))                      
        enddo                                                           
        do i=1,NDOC                                                     
         in = NO1+i                                                     
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj(in,USER(N2),USER(N9))                      &
                - PRODCWQWj(in,USER(N3),USER(N10))                      
         in = na+ndoc-i+1                                               
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj(in,USER(N2),USER(N9))                      &
                - PRODCWQWj(in,USER(N3),USER(N10))                      
        enddo                                                           
        if(NSOC>0)then                                                  
         do i=NDOC+1,NDNS                                               
          in = NO1+i                                                    
          ENERGY = ENERGY + 2.0d0*USER(N1-1+in)*USER(N8-1+in)          &
                 + PRODCWQWj(in,USER(N2),USER(N9))                     &
                 - PRODCWQWj(in,USER(N3),USER(N10))
         enddo
        endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ELSE                  ! PNOFi(Nc): Extended PNOF (NCWO>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ENERGY = 0.0d0
        do in=1,NO1
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj(in,USER(N2),USER(N9))                      &
                - PRODCWQWj(in,USER(N3),USER(N10))                      
        enddo                                                           
        do i=1,NDOC                                                     
         in = NO1+i                                                     
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj(in,USER(N2),USER(N9))                      &
                - PRODCWQWj(in,USER(N3),USER(N10))                      
         do iw=1,ncwo-1                                                 
          in = na+ncwo*(ndoc-i)+iw                                      
          ENERGY = ENERGY + USER(N1-1+in)                              &
                 * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )      &
                 + PRODCWQWj(in,USER(N2),USER(N9))                     &
                 - PRODCWQWj(in,USER(N3),USER(N10))                     
         enddo                                                          
         in = na+ncwo*(ndoc-i)+ncwo                                     
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj(in,USER(N2),USER(N9))                      &
                - PRODCWQWj(in,USER(N3),USER(N10))                      
        enddo                                                           
        if(NSOC>0)then                                                  
         do i=NDOC+1,NDNS                                               
          in = NO1+i                                                    
          ENERGY = ENERGY + 2.0d0*USER(N1-1+in)*USER(N8-1+in)          &
                 + PRODCWQWj(in,USER(N2),USER(N9))                     &
                 - PRODCWQWj(in,USER(N3),USER(N10))
         enddo
        endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ENDIF
!- - - including Electric Field  - - - - - - - - - - - - - - - - -
       if(EFIELDL)then
        CALL DIPMOMr(USER(N11),USER(N12),USER(N13),USER(N14),          &
                     USER(N15),USER(N16),USER(N17),USER(N7),           &
                     USER(N1),DMXe,DMYe,DMZe,DMX,DMY,DMZ,DM)
        ENERGY = ENERGY - EX*DMX - EY*DMY - EZ*DMZ
       endif
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --     
      else if(MSpin>0)then
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      High-Spin Multiplet State (S>0,Ms=S)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       IF(NCWO==1)THEN       ! PNOFi(1): Perfect Pairing (NCWO=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ENERGY = 0.0d0
        do in=1,NO1
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj1(in,USER(N2),USER(N9))                     &
                - PRODCWQWj1(in,USER(N3),USER(N10))                    &
                + 2.0d0*PRODROQWj1(in,USER(N1),USER(N9))               &
                - PRODROQWj1(in,USER(N1),USER(N10))                     
        enddo                                                           
        do i=1,NDOC                                                     
         in = NO1+i                                                     
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj1(in,USER(N2),USER(N9))                     &
                - PRODCWQWj1(in,USER(N3),USER(N10))                    &
                + 2.0d0*PRODROQWj1(in,USER(N1),USER(N9))               &
                - PRODROQWj1(in,USER(N1),USER(N10))                     
         in = na+ndoc-i+1                                               
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj2(in,USER(N2),USER(N9))                     &
                - PRODCWQWj2(in,USER(N3),USER(N10))                    &
                + 2.0d0*PRODROQWj2(in,USER(N1),USER(N9))               &
                - PRODROQWj2(in,USER(N1),USER(N10))     
        enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        do i=NDOC+1,NDNS
         in = NO1+i
         ENERGY = ENERGY + USER(N1-1+in)*USER(N8-1+in)                 &
                + 0.5d0*(PRODROQWj0(in,USER(N1),USER(N9))              &
                - PRODROQWj0(in,USER(N1),USER(N10)))
        enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ELSE                  ! PNOFi(Nc): Extended PNOF (NCWO>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ENERGY = 0.0d0
        do in=1,NO1
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj1(in,USER(N2),USER(N9))                     &
                - PRODCWQWj1(in,USER(N3),USER(N10))                    &
                + 2.0d0*PRODROQWj1(in,USER(N1),USER(N9))               &
                - PRODROQWj1(in,USER(N1),USER(N10))                     
        enddo                                                           
        do i=1,NDOC                                                     
         in = NO1+i                                                     
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj1(in,USER(N2),USER(N9))                     &
                - PRODCWQWj1(in,USER(N3),USER(N10))                    &
                + 2.0d0*PRODROQWj1(in,USER(N1),USER(N9))               &
                - PRODROQWj1(in,USER(N1),USER(N10))                     
         do iw=1,ncwo-1                                                 
          in = na+ncwo*(ndoc-i)+iw                                      
          ENERGY = ENERGY + USER(N1-1+in)                              &
                 * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )      &
                 + PRODCWQWj2(in,USER(N2),USER(N9))                    &
                 - PRODCWQWj2(in,USER(N3),USER(N10))                   &
                 + 2.0d0*PRODROQWj2(in,USER(N1),USER(N9))              &
                 - PRODROQWj2(in,USER(N1),USER(N10))                    
         enddo                                                          
         in = na+ncwo*(ndoc-i)+ncwo                                     
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj2(in,USER(N2),USER(N9))                     &
                - PRODCWQWj2(in,USER(N3),USER(N10))                    &
                + 2.0d0*PRODROQWj2(in,USER(N1),USER(N9))               &
                - PRODROQWj2(in,USER(N1),USER(N10))          
        enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        do i=NDOC+1,NDNS
         in = NO1+i
         ENERGY = ENERGY + USER(N1-1+in)*USER(N8-1+in)                 &
                + 0.5d0*(PRODROQWj0(in,USER(N1),USER(N9))              &
                - PRODROQWj0(in,USER(N1),USER(N10)))
        enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ENDIF
!- - - including Electric Field  - - - - - - - - - - - - - - - - -
       if(EFIELDL)then
        CALL DIPMOMr(USER(N11),USER(N12),USER(N13),USER(N14),          &
                     USER(N15),USER(N16),USER(N17),USER(N7),           &
                     USER(N1),DMXe,DMYe,DMZe,DMX,DMY,DMZ,DM)
        ENERGY = ENERGY - EX*DMX - EY*DMY - EZ*DMZ
       endif
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      end if
!-----------------------------------------------------------------------
      EELEC = ENERGY
      RETURN
      END

! CALCOG      
      SUBROUTINE CALCOG(NV,GAMMA,NF,GRAD,IUSER,USER)
      USE PARCOM 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(1) :: IUSER 
      DOUBLE PRECISION,DIMENSION(NV) :: GAMMA,GRAD
      DOUBLE PRECISION,DIMENSION(NUSER) :: USER
!-----------------------------------------------------------------------
!     Avoiding warnings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NF = NF
      IUSER(1) = IUSER(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL OCUPACIONr(GAMMA,USER(N1),USER(N2),USER(N3),
!     &                      USER(N4),USER(N5),USER(N6),NV)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --      
      if(MSpin==0)then
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      Singlet State (S=0,Ms=0) and Multiplet States (S>0,Ms=0)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       IF(NCWO==1)THEN       ! PNOFi(1): Perfect Pairing (NCWO=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        GRAD = 0.0d0
        DO ig=1,NV
         do i=1,NDOC
          in = NO1+i
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk(nv,in,ig,USER(N5),USER(N9))   &
                             - PRODCWQWjk(nv,in,ig,USER(N6),USER(N10)) ) 
          in = na+ndoc-i+1                                               
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk(nv,in,ig,USER(N5),USER(N9))   &
                             - PRODCWQWjk(nv,in,ig,USER(N6),USER(N10)) )
         enddo
        ENDDO
!- - -  including Electric Field  - - - - - - - - - - - - - - - - -
        if(EFIELDL)then   
         DO ig=1,NV
          do i=1,NDOC
           in = NO1+i
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )                             
           in = na+ndoc-i+1                                             
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )
          enddo
         ENDDO
        endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ELSE                  ! PNOFi(Nc): Extended PNOF (NCWO>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        GRAD = 0.0d0
        DO ig=1,NV
         do i=1,NDOC
          in = NO1+i
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk(nv,in,ig,USER(N5),USER(N9))   &
                             - PRODCWQWjk(nv,in,ig,USER(N6),USER(N10)) ) 
          do iw=1,NCWO-1                                                 
           in = na+ncwo*(ndoc-i)+iw                                      
           GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)              &
                    * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )  &
                    + 2.0d0 * (PRODCWQWjk(nv,in,ig,USER(N5),USER(N9))   &
                              -PRODCWQWjk(nv,in,ig,USER(N6),USER(N10)) ) 
          enddo                                                          
          in = na+ncwo*(ndoc-i)+ncwo                                     
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk(nv,in,ig,USER(N5),USER(N9))   &
                             - PRODCWQWjk(nv,in,ig,USER(N6),USER(N10)) )
         enddo
        ENDDO
!- - -  including Electric Field  - - - - - - - - - - - - - - - -
        if(EFIELDL)then   
         DO ig=1,NV
          do i=1,NDOC
           in = NO1+i
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )                             
           do iw=1,NCWO-1                                               
            in = na+ncwo*(ndoc-i)+iw                                    
            GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)      &
                     * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)         &
                       + EZ*USER(N17-1+in) )                            
           enddo                                                        
           in = na+ncwo*(ndoc-i)+ncwo                                   
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )
          enddo
         ENDDO
        endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ENDIF
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --     
      else if(MSpin>0)then
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      High-Spin Multiplet State (S>0,Ms=S)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       IF(NCWO==1)THEN       ! PNOFi(1): Perfect Pairing (NCWO=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        GRAD = 0.0d0
        DO ig=1,NV
         do i=1,NDOC
          in = NO1+i
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk1(nv,in,ig,USER(N5),USER(N9))  &
                             - PRODCWQWjk1(nv,in,ig,USER(N6),USER(N10)))&
                   + 2.0d0 *   PRODDRQWjk1(nv,in,ig,USER(N4),USER(N9))  &
                   -           PRODDRQWjk1(nv,in,ig,USER(N4),USER(N10))  
          in = na+ndoc-i+1                                               
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk2(nv,in,ig,USER(N5),USER(N9))  &
                             - PRODCWQWjk2(nv,in,ig,USER(N6),USER(N10)))&
                   + 2.0d0 *   PRODDRQWjk2(nv,in,ig,USER(N4),USER(N9))  &
                   -           PRODDRQWjk2(nv,in,ig,USER(N4),USER(N10))
         enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ENDDO
! - - - including Electric Field - - - - - - - - - - - - - - - - - - - - 
        if(EFIELDL)then   
         DO ig=1,NV
          do i=1,NDOC
           in = NO1+i
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )                             
           in = na+ndoc-i+1                                             
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )
          enddo
         ENDDO
        endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ELSE                  ! PNOFi(Nc): Extended PNOF (NCWO>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        GRAD = 0.0d0
        DO ig=1,NV
         do i=1,NDOC
          in = NO1+i
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk1(nv,in,ig,USER(N5),USER(N9))  &
                             - PRODCWQWjk1(nv,in,ig,USER(N6),USER(N10)))&
                   + 2.0d0 *   PRODDRQWjk1(nv,in,ig,USER(N4),USER(N9))  &
                   -           PRODDRQWjk1(nv,in,ig,USER(N4),USER(N10))  
          do iw=1,NCWO-1                                                 
           in = na+ncwo*(ndoc-i)+iw                                      
           GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)              &
                    * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )  &
                    + 2.0d0 * (PRODCWQWjk2(nv,in,ig,USER(N5),USER(N9))  &
                              -PRODCWQWjk2(nv,in,ig,USER(N6),USER(N10)))&
                    + 2.0d0 *  PRODDRQWjk2(nv,in,ig,USER(N4),USER(N9))  &
                    -          PRODDRQWjk2(nv,in,ig,USER(N4),USER(N10))  
          enddo                                                          
          in = na+ncwo*(ndoc-i)+ncwo                                     
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk2(nv,in,ig,USER(N5),USER(N9))  &
                             - PRODCWQWjk2(nv,in,ig,USER(N6),USER(N10)))&
                   + 2.0d0 *   PRODDRQWjk2(nv,in,ig,USER(N4),USER(N9))  &
                   -           PRODDRQWjk2(nv,in,ig,USER(N4),USER(N10))
         enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -          
        ENDDO
! - - - including Electric Field - - - - - - - - - - - - - - - - - - - -
        if(EFIELDL)then   
         DO ig=1,NV
          do i=1,NDOC
           in = NO1+i
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )                             
           do iw=1,NCWO-1                                               
            in = na+ncwo*(ndoc-i)+iw                                    
            GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)      &
                     * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)         &
                       + EZ*USER(N17-1+in) )                            
           enddo                                                        
           in = na+ncwo*(ndoc-i)+ncwo                                   
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )
          enddo
         ENDDO
        endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ENDIF
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      end if
!-----------------------------------------------------------------------
      RETURN
      END

!----------------------------------------------------------------------!
!                                                                      !
!      SUMSL minimizes a general unconstrained objective function      !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    J E Dennis, David Gay, and R E Welsch,                            !
!    An Adaptive Nonlinear Least-squares Algorithm,                    !
!    ACM Transactions on Mathematical Software,                        !
!    Volume 7, Number 3, 1981.                                         !
!                                                                      !
!    J E Dennis, H H W Mei,                                            !
!    Two New Unconstrained Optimization Algorithms Which Use           !
!    Function and Gradient Values,                                     !
!    Journal of Optimization Theory and Applications,                  !
!    Volume 28, pages 453-482, 1979.                                   !
!                                                                      !
!    J E Dennis, Jorge More,                                           !
!    Quasi-Newton Methods, Motivation and Theory,                      !
!    SIAM Review,                                                      !
!    Volume 19, pages 46-89, 1977.                                     !
!                                                                      !
!    D Goldfarb,                                                       !
!    Factorized Variable Metric Methods for Unconstrained Optimization,!
!    Mathematics of Computation,                                       !
!    Volume 30, pages 796-811, 1976.                                   !
!                                                                      !
!----------------------------------------------------------------------!

subroutine sumsl(n,d,x,calcf,calcg,iv,liv,lv,v,uiparm,urparm)
!
! The routine uses analytic gradient and hessian approximation from
! the secant update. It interacts with subroutine sumit in an attempt
! to find an n-vector  x*  that minimizes the (unconstrained)
! objective function computed by  calcf.  (often the  x*  found is
! a local minimizer rather than a global one.)
!
! Parameters:
!
! n  (input) the number of variables on which  f  depends, i.e.,
!                  the number of components in  x.
! d  (input/output) a scale vector such that  d(i)*x(i),
!                  i = 1,2,...,n,  are all in comparable units.
!                  d can strongly affect the behavior of sumsl.
!                  finding the best choice of d is generally a trial-
!                  and-error process.  choosing d so that d(i)*x(i)
!                  has about the same value for all i often works well.
!                  the defaults provided by subroutine deflt (see iv
!                  below) require the caller to supply d.
! x........ (input/output) before (initially) calling sumsl, the call-
!                  er should set  x  to an initial guess at  x*.  when
!                  sumsl returns,  x  contains the best point so far
!                  found, i.e., the one that gives the least value so
!                  far seen for  f(x).
! calcf.... (input) a subroutine that, given x, computes f(x).  calcf
!                  must be declared external in the calling program.
!                  it is invoked by
!                    call calcf(n, x, nf, f, uiparm, urparm)
!                  when calcf is called, nf is the invocation
!                  count for calcf.  nf is included for possible use
!                  with calcg.  if x is out of bounds (e.g., if it
!                  would cause overflow in computing f(x)), then calcf
!                  should set nf to 0.  this will cause a shorter step
!                  to be attempted.  (if x is in bounds, then calcf
!                  should not change nf.)  the other parameters are as
!                  described above and below.  calcf should not change
!                  n, p, or x.
! calcg.... (input) a subroutine that, given x, computes g(x), the gra-
!                  dient of f at x.  calcg must be declared external in
!                  the calling program.  it is invoked by
!                    call calcg(n, x, nf, g, uiparm, urparm, ufaprm)
!                  when calcg is called, nf is the invocation
!                  count for calcf at the time f(x) was evaluated.  the
!                  x passed to calcg is usually the one passed to calcf
!                  on either its most recent invocation or the one
!                  prior to it.  if calcf saves intermediate results
!                  for use by calcg, then it is possible to tell from
!                  nf whether they are valid for the current x (or
!                  which copy is valid if two copies are kept).  if g
!                  cannot be computed at x, then calcg should set nf to
!                  0.  in this case, sumsl will return with iv(1) = 65.
!                  (if g can be computed at x, then calcg should not
!                  changed nf.)  the other parameters to calcg are as
!                  described above and below.  calcg should not change
!                  n or x.
! iv....... (input/output) an integer value array of length liv (see
!                  below) that helps control the sumsl algorithm and
!                  that is used to store various intermediate quanti-
!                  ties.  of particular interest are the initialization/
!                  return code iv(1) and the entries in iv that control
!                  printing and limit the number of iterations and func-
!                  tion evaluations.  see the section on iv input
!                  values below.
! liv...... (input) length of iv array.  must be at least 60.  if liv
!                  is too small, then sumsl returns with iv(1) = 15.
!                  when sumsl returns, the smallest allowed value of
!                  liv is stored in iv(lastiv) -- see the section on
!                  iv output values below.  (this is intended for use
!                  with extensions of sumsl that handle constraints.)
! lv....... (input) length of v array.  must be at least 71+n*(n+15)/2.
!                  (at least 77+n*(n+17)/2 for smsno, at least
!                  78+n*(n+12) for humsl).  if lv is too small, then
!                  sumsl returns with iv(1) = 16.  when sumsl returns,
!                  the smallest allowed value of lv is stored in
!                  iv(lastv) -- see the section on iv output values
!                  below.
! v........ (input/output) a floating-point value array of length lv
!                  (see below) that helps control the sumsl algorithm
!                  and that is used to store various intermediate
!                  quantities.  of particular interest are the entries
!                  in v that limit the length of the first step
!                  attempted (lmax0) and specify convergence tolerances
!                  (afctol, lmaxs, rfctol, sctol, xctol, xftol).
! uiparm... (input) user integer parameter array passed without change
!                  to calcf and calcg.
! urparm... (input) user floating-point parameter array passed without
!                  change to calcf and calcg.
!
!   iv input values (from subroutine deflt)
!
! iv(1)...  on input, iv(1) should have a value between 0 and 14......
!             0 and 12 mean this is a fresh start.  0 means that
!                  deflt(2, iv, liv, lv, v)
!             is to be called to provide all default values to iv and
!             v.  12 (the value that deflt assigns to iv(1)) means the
!             caller has already called deflt and has possibly changed
!             some iv and/or v entries to non-default values.
!             13 means deflt has been called and that sumsl (and
!             sumit) should only do their storage allocation.  that is,
!             they should set the output components of iv that tell
!             where various subarrays arrays of v begin, such as iv(g)
!             (and, for humsl and humit only, iv(dtol)), and return.
!             14 means that a storage has been allocated (by a call
!             with iv(1) = 13) and that the algorithm should be
!             started.  when called with iv(1) = 13, sumsl returns
!             iv(1) = 14 unless liv or lv is too small (or n is not
!             positive).  default = 12.
! iv(inith).... iv(25) tells whether the hessian approximation h should
!             be initialized.  1 (the default) means sumit should
!             initialize h to the diagonal matrix whose i-th diagonal
!             element is d(i)**2.  0 means the caller has supplied a
!             cholesky factor  l  of the initial hessian approximation
!             h = l*(l**t)  in v, starting at v(iv(lmat)) = v(iv(42))
!             (and stored compactly by rows).  note that iv(lmat) may
!             be initialized by calling sumsl with iv(1) = 13 (see
!             the iv(1) discussion above).  default = 1.
! iv(mxfcal)... iv(17) gives the maximum number of function evaluations
!             (calls on calcf) allowed.  if this number does not suf-
!             fice, then sumsl returns with iv(1) = 9.  default = 200.
! iv(mxiter)... iv(18) gives the maximum number of iterations allowed.
!             it also indirectly limits the number of gradient evalua-
!             tions (calls on calcg) to iv(mxiter) + 1.  if iv(mxiter)
!             iterations do not suffice, then sumsl returns with
!             iv(1) = 10.  default = 150.
! iv(outlev)... iv(19) controls the number and length of iteration sum-
!             mary lines printed (by itsum).  iv(outlev) = 0 means do
!             not print any summary lines.  otherwise, print a summary
!             line after each abs(iv(outlev)) iterations.  if iv(outlev)
!             is positive, then summary lines of length 78 (plus carri-
!             age control) are printed, including the following...  the
!             iteration and function evaluation counts, f = the current
!             function value, relative difference in function values
!             achieved by the latest step (i.e., reldf = (f0-v(f))/f01,
!             where f01 is the maximum of abs(v(f)) and abs(v(f0)) and
!             v(f0) is the function value from the previous itera-
!             tion), the relative function reduction predicted for the
!             step just taken (i.e., preldf = v(preduc) / f01, where
!             v(preduc) is described below), the scaled relative change
!             in x (see v(reldx) below), the step parameter for the
!             step just taken (stppar = 0 means a full newton step,
!             between 0 and 1 means a relaxed newton step, between 1
!             and 2 means a double dogleg step, greater than 2 means
!             a scaled down Cauchy step -- see subroutine dbldog), the
!             2-norm of the scale vector d times the step just taken
!             (see v(dstnrm) below), and npreldf, i.e.,
!             v(nreduc)/f01, where v(nreduc) is described below -- if
!             npreldf is positive, then it is the relative function
!             reduction predicted for a newton step (one with
!             stppar = 0).  if npreldf is negative, then it is the
!             negative of the relative function reduction predicted
!             for a step computed with step bound v(lmaxs) for use in
!             testing for singular convergence.
!                  if iv(outlev) is negative, then lines of length 50
!             are printed, including only the first 6 items listed
!             above (through reldx).
!             default = 1.
! iv(parprt)... iv(20) = 1 means print any nondefault v values on a
!             fresh start or any changed v values on a restart.
!             iv(parprt) = 0 means skip this printing.  default = 1.
! iv(prunit)... iv(21) is the output unit number on which all printing
!             is done.  iv(prunit) = 0 means suppress all printing.
!             default = standard output unit (unit 6 on most systems).
! iv(solprt)... iv(22) = 1 means print out the value of x returned (as
!             well as the gradient and the scale vector d).
!             iv(solprt) = 0 means skip this printing.  default = 1.
! iv(statpr)... iv(23) = 1 means print summary statistics upon return-
!             ing.  these consist of the function value, the scaled
!             relative change in x caused by the most recent step (see
!             v(reldx) below), the number of function and gradient
!             evaluations (calls on calcf and calcg), and the relative
!             function reductions predicted for the last step taken and
!             for a newton step (or perhaps a step bounded by v(lmaxs)
!             -- see the descriptions of preldf and npreldf under
!             iv(outlev) above).
!             iv(statpr) = 0 means skip this printing.
!             iv(statpr) = -1 means skip this printing as well as that
!             of the one-line termination reason message.  default = 1.
! iv(x0prt).... iv(24) = 1 means print the initial x and scale vector d
!             (on a fresh start only).  iv(x0prt) = 0 means skip this
!             printing.  default = 1.
!
!   (selected) iv output values
!
! iv(1)........ on output, iv(1) is a return code....
!             3 = x-convergence.  the scaled relative difference (see
!                  v(reldx)) between the current parameter vector x and
!                  a locally optimal parameter vector is very likely at
!                  most v(xctol).
!             4 = relative function convergence.  the relative differ-
!                  ence between the current function value and its lo-
!                  cally optimal value is very likely at most v(rfctol).
!             5 = both x- and relative function convergence (i.e., the
!                  conditions for iv(1) = 3 and iv(1) = 4 both hold).
!             6 = absolute function convergence.  the current function
!                  value is at most v(afctol) in absolute value.
!             7 = singular convergence.  the hessian near the current
!                  iterate appears to be singular or nearly so, and a
!                  step of length at most v(lmaxs) is unlikely to yield
!                  a relative function decrease of more than v(sctol).
!             8 = false convergence.  the iterates appear to be converg-
!                  ing to a noncritical point.  this may mean that the
!                  convergence tolerances (v(afctol), v(rfctol),
!                  v(xctol)) are too small for the accuracy to which
!                  the function and gradient are being computed, that
!                  there is an error in computing the gradient, or that
!                  the function or gradient is discontinuous near x.
!             9 = function evaluation limit reached without other con-
!                  vergence (see iv(mxfcal)).
!            10 = iteration limit reached without other convergence
!                  (see iv(mxiter)).
!            11 = STOPX returned .true. (external interrupt).  see the
!                  usage notes below.
!            14 = storage has been allocated (after a call with
!                  iv(1) = 13).
!            17 = restart attempted with n changed.
!            18 = d has a negative component and iv(dtype) <= 0.
!            19...43 = v(iv(1)) is out of range.
!            63 = f(x) cannot be computed at the initial x.
!            64 = bad parameters passed to assess (which should not
!                  occur).
!            65 = the gradient could not be computed at x (see calcg
!                  above).
!            67 = bad first parameter to deflt.
!            80 = iv(1) was out of range.
!            81 = n is not positive.
! iv(g)........ iv(28) is the starting subscript in v of the current
!             gradient vector (the one corresponding to x).
! iv(lastiv)... iv(44) is the least acceptable value of liv.  (it is
!             only set if liv is at least 44.)
! iv(lastv).... iv(45) is the least acceptable value of lv.  (it is
!             only set if liv is large enough, at least iv(lastiv).)
! iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e.,
!             function evaluations).
! iv(ngcall)... iv(30) is the number of gradient evaluations (calls on
!             calcg).
! iv(niter).... iv(31) is the number of iterations performed.
!
!   (selected) v input values (from subroutine deflt)
!
! v(bias)..... v(43) is the bias parameter used in subroutine dbldog --
!             see that subroutine for details.  default = 0.8.
! v(afctol)... v(31) is the absolute function convergence tolerance.
!             if sumsl finds a point where the function value is less
!             than v(afctol) in absolute value, and if sumsl does not
!             return with iv(1) = 3, 4, or 5, then it returns with
!             iv(1) = 6.  this test can be turned off by setting
!             v(afctol) to zero.  default = max(10**-20, machep**2),
!             where machep is the unit roundoff.
! v(dinit).... v(38), if nonnegative, is the value to which the scale
!             vector d is initialized.  default = -1.
! v(lmax0).... v(35) gives the maximum 2-norm allowed for d times the
!             very first step that sumsl attempts.  this parameter can
!             markedly affect the performance of sumsl.
! v(lmaxs).... v(36) is used in testing for singular convergence -- if
!             the function reduction predicted for a step of length
!             bounded by v(lmaxs) is at most v(sctol) * abs(f0), where
!             f0  is the function value at the start of the current
!             iteration, and if sumsl does not return with iv(1) = 3,
!             4, 5, or 6, then it returns with iv(1) = 7.  default = 1.
! v(rfctol)... v(32) is the relative function convergence tolerance.
!             if the current model predicts a maximum possible function
!             reduction (see v(nreduc)) of at most v(rfctol)*abs(f0)
!             at the start of the current iteration, where  f0  is the
!             then current function value, and if the last step attempt-
!             ed achieved no more than twice the predicted function
!             decrease, then sumsl returns with iv(1) = 4 (or 5).
!             default = max(10**-10, machep**(2/3)), where machep is
!             the unit roundoff.
! v(sctol).... v(37) is the singular convergence tolerance -- see the
!             description of v(lmaxs) above.
! v(tuner1)... v(26) helps decide when to check for false convergence.
!             this is done if the actual function decrease from the
!             current step is no more than v(tuner1) times its predict-
!             ed value.  default = 0.1.
! v(xctol).... v(33) is the x-convergence tolerance.  if a newton step
!             (see v(nreduc)) is tried that has v(reldx) <= v(xctol)
!             and if this step yields at most twice the predicted func-
!             tion decrease, then sumsl returns with iv(1) = 3 (or 5).
!             (see the description of v(reldx) below.)
!             default = machep**0.5, where machep is the unit roundoff.
! v(xftol).... v(34) is the false convergence tolerance.  if a step is
!             tried that gives no more than v(tuner1) times the predict-
!             ed function decrease and that has v(reldx) <= v(xftol),
!             and if sumsl does not return with iv(1) = 3, 4, 5, 6, or
!             7, then it returns with iv(1) = 8.  (see the description
!             of v(reldx) below.)  default = 100*machep, where
!             machep is the unit roundoff.
! v(*)........ deflt supplies to v a number of tuning constants, with
!             which it should ordinarily be unnecessary to tinker.  see
!             section 17 of version 2.2 of the nl2sol usage summary
!             (i.e., the appendix to ref. 1) for details on v(i),
!             i = decfac, incfac, phmnfc, phmxfc, rdfcmn, rdfcmx,
!             tuner2, tuner3, tuner4, tuner5.
!
!   (selected) v output values
!
! v(dgnorm)... v(1) is the 2-norm of (diag(d)**-1)*g, where g is the
!             most recently computed gradient.
! v(dstnrm)... v(2) is the 2-norm of diag(d)*step, where step is the
!             current step.
! v(f)........ v(10) is the current function value.
! v(f0)....... v(13) is the function value at the start of the current
!             iteration.
! v(nreduc)... v(6), if positive, is the maximum function reduction
!             possible according to the current model, i.e., the func-
!             tion reduction predicted for a newton step (i.e.,
!             step = -h**-1 * g,  where  g  is the current gradient and
!             h is the current hessian approximation).
!                  if v(nreduc) is negative, then it is the negative of
!             the function reduction predicted for a step computed with
!             a step bound of v(lmaxs) for use in testing for singular
!             convergence.
! v(preduc)... v(7) is the function reduction predicted (by the current
!             quadratic model) for the current step.  this (divided by
!             v(f0)) is used in testing for relative function
!             convergence.
! v(reldx).... v(17) is the scaled relative change in x caused by the
!             current step, computed as
!                  max(abs(d(i)*(x(i)-x0(i)), 1 <= i <= p) /
!                     max(d(i)*(abs(x(i))+abs(x0(i))), 1 <= i <= p),
!             where x = x0 + step.
!
!  notes
!
!   algorithm notes
!
!        this routine uses a hessian approximation computed from the
!     bfgs update (see ref 3).  only a cholesky factor of the hessian
!     approximation is stored, and this is updated using ideas from
!     ref. 4.  steps are computed by the double dogleg scheme described
!     in ref. 2.  the steps are assessed as in ref. 1.
!
!   usage notes
!
!        after a return with iv(1) <= 11, it is possible to restart,
!     i.e., to change some of the iv and v input values described above
!     and continue the algorithm from the point where it was interrupt-
!     ed.  iv(1) should not be changed, nor should any entries of iv
!     and v other than the input values (those supplied by deflt).
!        those who do not wish to write a calcg which computes the
!     gradient analytically should call smsno rather than sumsl.
!     smsno uses finite differences to compute an approximate gradient.
!        those who would prefer to provide f and g (the function and
!     gradient) by reverse communication rather than by writing subrou-
!     tines calcf and calcg may call on sumit directly.  see the com-
!     ments at the beginning of sumit.
!        those who use sumsl interactively may wish to supply their
!     own STOPX function, which should return .true. if the break key
!     has been pressed since STOPX was last invoked.  this makes it
!     possible to externally interrupt sumsl (which will return with
!     iv(1) = 11 if STOPX returns .true.).
!        storage for g is allocated at the end of v.  thus the caller
!     may make v longer than specified above and may allow calcg to use
!     elements of g beyond the first n as scratch storage.
!
!   portability notes
!
!        the sumsl distribution tape contains both single- and double-
!     precision versions of the sumsl source code, so it should be un-
!     necessary to change precisions.
!        only the function rmdcon contains machine-dependent
!     constants.  to change from one machine to another, it should
!     suffice to change the (few) relevant lines in these functions.
!        intrinsic functions are explicitly declared.  on certain com-
!     puters (e.g. univac), it may be necessary to comment out these
!     declarations.  so that this may be done automatically by a simple
!     program, such declarations are preceded by a comment having c/+
!     in columns 1-3 and blanks in columns 4-72 and are followed by
!     a comment having c/ in columns 1 and 2 and blanks in columns 3-72.
!        the sumsl source code is expressed in 1966 ansi standard
!     fortran.  it may be converted to fortran 77 by commenting out all
!     lines that fall between a line having c/6 in columns 1-3 and a
!     line having c/7 in columns 1-3 and by removing (i.e., replacing
!     by a blank) the c in column 1 of the lines that follow the c/7
!     line and precede a line having c/ in columns 1-2 and blanks in
!     columns 3-72.  these changes convert some data statements into
!     parameter statements, convert some variables from real to
!     character*4, and make the data statements that initialize these
!     variables use character strings delimited by primes instead
!     of hollerith constants.  (such variables and data statements
!     appear only in modules itsum and parck.  parameter statements
!     appear nearly everywhere.)  these changes also add save state-
!     ments for variables given machine-dependent constants by rmdcon.
!
  implicit double precision (a-h,o-z)   
  integer n, liv, lv
  integer iv(liv), uiparm(*)
  double precision d(n), x(n), v(lv), urparm(*)
!     dimension v(71 + n*(n+15)/2), uiparm(*), urparm(*)

  integer g1, iv1, nf
  double precision f
  integer nextv, nfcall, nfgcal, g, toobig, vneed

  parameter (nextv=47, nfcall=6, nfgcal=7, g=28, toobig=2, vneed=4)
  
  if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
  iv1 = iv(1)
  if (iv1 == 12 .or. iv1 == 13) iv(vneed) = iv(vneed) + n
  if (iv1 == 14) go to 10
  if (iv1 > 2 .and. iv1 < 12) go to 10
  g1 = 1
  if (iv1 == 12) iv(1) = 13
  go to 20

 10   g1 = iv(g)

 20   call sumit(d, f, v(g1), iv, liv, lv, n, v, x)
  if (iv(1) - 2) 30, 40, 50

 30   nf = iv(nfcall)
  call calcf(n, x, nf, f, uiparm, urparm)
  if (nf <= 0) iv(toobig) = 1
  go to 20

 40   call calcg(n, x, iv(nfgcal), v(g1), uiparm, urparm)
  go to 20

 50   if (iv(1) /= 14) then
        return
      end if
!
!  Storage allocation
!
  iv(g) = iv(nextv)
  iv(nextv) = iv(g) + n
  if (iv1 /= 13) go to 10

  return
end

subroutine deflt ( alg, iv, liv, lv, v )
!
!! DEFLT: supply default values to IV and V.
!
!  Discussion:
!
!   ALG = 1 means regression constants.
!   ALG = 2 means general unconstrained optimization constants.
!
  integer liv
  integer lv

  integer alg
  integer iv(liv)
  double precision v(lv)
  external vdflt
  integer miv, mv
  integer miniv(2), minv(2)
  integer algsav, covprt, covreq, dtype, hc, ierr, inith, inits
  integer ipivot, ivneed, lastiv, lastv, lmat, mxfcal, mxiter
  integer nfcov, ngcov, nvdflt, outlev, parprt, parsav, perm
  integer prunit, qrtyp, rdreq, rmat, solprt, statpr, vneed
  integer vsave, x0prt

  parameter (algsav=51, covprt=14, covreq=15, dtype=16, hc=71 )
  parameter (ierr=75, inith=25, inits=25, ipivot=76, ivneed=3 )
  parameter (lastiv=44, lastv=45, lmat=42, mxfcal=17, mxiter=18 )
  parameter (nfcov=52, ngcov=53, nvdflt=50, outlev=19, parprt=20 )
  parameter (parsav=49, perm=58, prunit=21, qrtyp=80, rdreq=57 )
  parameter (rmat=78, solprt=22, statpr=23, vneed=4, vsave=60 )
  parameter (x0prt=24)

  data miniv(1)/80/, miniv(2)/59/, minv(1)/98/, minv(2)/71/

  if ( alg < 1 .or. 2 < alg ) then
    iv(1) = 67
    return
  end if

  miv = miniv(alg)

  if ( liv < miv ) then
    iv(1) = 15
    return
  end if

  mv = minv(alg)

  if ( lv < mv ) then
    iv(1) = 16
    return
  end if

  call vdflt(alg, lv, v)
  iv(1) = 12
  iv(algsav) = alg
  iv(ivneed) = 0
  iv(lastiv) = miv
  iv(lastv) = mv
  iv(lmat) = mv + 1
  iv(mxfcal) = 200
  iv(mxiter) = 150
  iv(outlev) = 1
  iv(parprt) = 1
  iv(perm) = miv + 1
! change  
  iv(prunit) = 2
! change  
  iv(solprt) = 1
  iv(statpr) = 1
  iv(vneed) = 0
  iv(x0prt) = 1
!
!  General optimization values.
!
  if ( 2 <= alg ) then

    iv(dtype) = 0
    iv(inith) = 1
    iv(nfcov) = 0
    iv(ngcov) = 0
    iv(nvdflt) = 25
    iv(parsav) = 47
!
!  Regression values.
!
  else

    iv(covprt) = 3
    iv(covreq) = 1
    iv(dtype) = 1
    iv(hc) = 0
    iv(ierr) = 0
    iv(inits) = 0
    iv(ipivot) = 0
    iv(nvdflt) = 32
    iv(parsav) = 67
    iv(qrtyp) = 1
    iv(rdreq) = 3
    iv(rmat) = 0
    iv(vsave) = 58

  end if

  return
end

subroutine vdflt ( alg, lv, v )
!
!! VDFLT supplies default values to V.
!
!  Discussion:
!
!    alg = 1 means regression constants.
!    alg = 2 means general unconstrained optimization constants.
!
  implicit none

  integer alg, lv
  double precision v(lv)
  double precision rmdcon
  double precision machep, mepcrt, one, sqteps, three
  integer afctol, bias, cosmin, decfac, delta0, dfac, dinit, dltfdc
  integer dltfdj, dtinit, d0init, epslon, eta0, fuzz, huberc
  integer incfac, lmax0, lmaxs, phmnfc, phmxfc, rdfcmn, rdfcmx
  integer rfctol, rlimit, rsptol, sctol, sigmin, tuner1, tuner2
  integer tuner3, tuner4, tuner5, xctol, xftol

  parameter (one=1.d+0, three=3.d+0)

  parameter (afctol=31, bias=43, cosmin=47, decfac=22, delta0=44 )
  parameter ( dfac=41, dinit=38, dltfdc=42, dltfdj=43, dtinit=39 )
  parameter ( d0init=40, epslon=19, eta0=42, fuzz=45, huberc=48 )
  parameter ( incfac=23, lmax0=35, lmaxs=36, phmnfc=20, phmxfc=21 )
  parameter ( rdfcmn=24, rdfcmx=25, rfctol=32, rlimit=46, rsptol=49 )
  parameter ( sctol=37, sigmin=50, tuner1=26, tuner2=27, tuner3=28 )
  parameter ( tuner4=29, tuner5=30, xctol=33, xftol=34)

  machep = rmdcon(3)
  v(afctol) = 1.d-20

  if ( machep > 1.d-10 ) then
    v(afctol) = machep**2
  end if

  v(decfac) = 0.5d+0
  sqteps = rmdcon(4)
  v(dfac) = 0.6d+0
  v(delta0) = sqteps
  v(dtinit) = 1.d-6
  mepcrt = machep ** (one/three)
  v(d0init) = 1.d+0
  v(epslon) = 0.1d+0
  v(incfac) = 2.d+0
  v(lmax0) = 1.d+0
  v(lmaxs) = 1.d+0
  v(phmnfc) = -0.1d+0
  v(phmxfc) = 0.1d+0
  v(rdfcmn) = 0.1d+0
  v(rdfcmx) = 4.d+0
  v(rfctol) = max (1.d-10, mepcrt**2)
  v(sctol) = v(rfctol)
  v(tuner1) = 0.1d+0
  v(tuner2) = 1.d-4
  v(tuner3) = 0.75d+0
  v(tuner4) = 0.5d+0
  v(tuner5) = 0.75d+0
  v(xctol) = sqteps
  v(xftol) = 1.d+2 * machep

  if ( alg < 2 ) then
    v(cosmin) = max (1.d-6, 1.d+2 * machep)
    v(dinit) = 0.d+0
    v(dltfdc) = mepcrt
    v(dltfdj) = sqteps
    v(fuzz) = 1.5d+0
    v(huberc) = 0.7d+0
    v(rlimit) = rmdcon(5)
    v(rsptol) = 1.d-3
    v(sigmin) = 1.d-4
  else
    v(bias) = 0.8d+0
    v(dinit) = -1.0d+0
    v(eta0) = 1.0d+3 * machep
  end if

  return
end

subroutine sumit ( d, fx, g, iv, liv, lv, n, v, x)
!
!  SUMIT carries out unconstrained minimization iterations for SUMSL.
!
!  Discussion:
!
!    The routine uses double-dogleg/BFGS steps.
!
!    parameters iv, n, v, and x are the same as the corresponding
!    ones to sumsl (which see), except that v can be shorter (since
!    the part of v that sumsl uses for storing g is not needed).
!    moreover, compared with sumsl, iv(1) may have the two additional
!    output values 1 and 2, which are explained below, as is the use
!    of iv(toobig) and iv(nfgcal).  the value iv(g), which is an
!    output value from sumsl (and smsno), is not referenced by
!    sumit or the subroutines it calls.
!
!    fx and g need not have been initialized when sumit is called
!    with iv(1) = 12, 13, or 14.
!
! iv(1) = 1 means the caller should set fx to f(x), the function value
!             at x, and call sumit again, having changed none of the
!             other parameters.  an exception occurs if f(x) cannot be
!             (e.g. if overflow would occur), which may happen because
!             of an oversized step.  in this case the caller should set
!             iv(toobig) = iv(2) to 1, which will cause sumit to ig-
!             nore fx and try a smaller step.  the parameter nf that
!             sumsl passes to calcf (for possible use by calcg) is a
!             copy of iv(nfcall) = iv(6).
! iv(1) = 2 means the caller should set g to g(x), the gradient vector
!             of f at x, and call sumit again, having changed none of
!             the other parameters except possibly the scale vector d
!             when iv(dtype) = 0.  the parameter nf that sumsl passes
!             to calcg is iv(nfgcal) = iv(7).  if g(x) cannot be
!             evaluated, then the caller may set iv(nfgcal) to 0, in
!             which case sumit will return with iv(1) = 65.
!
!  Parameters:
!
! d.... scale vector.
! fx... function value.
! g.... gradient vector.
! iv... integer value array.
! liv.. length of iv (at least 60).
! lv... length of v (at least 71 + n*(n+13)/2).
! n.... number of variables (components in x and g).
! v.... floating-point value array.
! x.... vector of parameters to be optimized.
!
  integer liv
  integer lv
  integer n

  integer iv(liv)
  double precision d(n)
  double precision fx
  double precision g(n)
  double precision v(lv)
  double precision x(n)
  integer dg1, g01, i, k, l, lstgst, nwtst1, step1
  integer        temp1, w, x01, z
  double precision t
  double precision half, negone, one, onep2, zero
  logical stopx
  double precision dotprd, reldst, v2norm
  integer cnvcod, dg, dgnorm, dinit, dstnrm, dst0, f, f0, fdif
  integer gthg, gtstep, g0, incfac, inith, irc, kagqt, lmat, lmax0
  integer lmaxs, mode, model, mxfcal, mxiter, nextv, nfcall, nfgcal
  integer ngcall, niter, nreduc, nwtstp, preduc, radfac, radinc
  integer radius, rad0, reldx, restor, step, stglim, stlstg, toobig
  integer tuner4, tuner5, vneed, xirc, x0

  parameter (cnvcod=55, dg=37, g0=48, inith=25, irc=29, kagqt=33 )
  parameter ( mode=35, model=5, mxfcal=17, mxiter=18, nfcall=6 )
  parameter ( nfgcal=7, ngcall=30, niter=31, nwtstp=34, radinc=8 )
  parameter ( restor=9, step=40, stglim=11, stlstg=41, toobig=2 )
  parameter ( vneed=4, xirc=13, x0=43)

  parameter (dgnorm=1, dinit=38, dstnrm=2, dst0=3, f=10, f0=13 )
  parameter ( fdif=11, gthg=44, gtstep=4, incfac=23, lmat=42 )
  parameter ( lmax0=35, lmaxs=36, nextv=47, nreduc=6, preduc=7 )
  parameter ( radfac=16, radius=8, rad0=9, reldx=17, tuner4=29 )
  parameter ( tuner5=30)

  parameter (half=0.5d+0, negone=-1.d+0, one=1.d+0, onep2=1.2d+0, zero=0.d+0)
!
  i = iv(1)
  if (i == 1) go to 50
  if (i == 2) go to 60
!
!   check validity of iv and v input values
!
  if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
  if (iv(1) == 12 .or. iv(1) == 13) then
    iv(vneed) = iv(vneed) + n*(n+13)/2
  end if
  call parck(2, d, iv, liv, lv, n, v)
  i = iv(1) - 2
  if (i > 12) then
    return
  end if
  go to (180, 180, 180, 180, 180, 180, 120, 90, 120, 10, 10, 20), i
!
!   storage allocation
!
10    l = iv(lmat)
  iv(x0) = l + n*(n+1)/2
  iv(step) = iv(x0) + n
  iv(stlstg) = iv(step) + n
  iv(g0) = iv(stlstg) + n
  iv(nwtstp) = iv(g0) + n
  iv(dg) = iv(nwtstp) + n
  iv(nextv) = iv(dg) + n
  if (iv(1) /= 13) go to 20
     iv(1) = 14
     return
!
!   initialization
!
 20   iv(niter) = 0
  iv(nfcall) = 1
  iv(ngcall) = 1
  iv(nfgcal) = 1
  iv(mode) = -1
  iv(model) = 1
  iv(stglim) = 1
  iv(toobig) = 0
  iv(cnvcod) = 0
  iv(radinc) = 0
  v(rad0) = 0.0D+00
  if (v(dinit) >= 0.0D+00) call vscopy(n, d, v(dinit))
  if (iv(inith) /= 1) go to 40
!
!  set the initial hessian approximation to diag(d)**-2
!
     l = iv(lmat)
     call vscopy(n*(n+1)/2, v(l), zero)
     k = l - 1

     do i = 1, n
       k = k + i
       t = d(i)
       if (t <= 0.0D+00) t = one
       v(k) = t
     end do
!
!  compute initial function value
!
 40   iv(1) = 1
  return

 50   v(f) = fx
  if (iv(mode) >= 0) go to 180
  iv(1) = 2
  if (iv(toobig) == 0) then
    return
  end if
     iv(1) = 63
     go to 300
!
!   make sure gradient could be computed
!
 60   if (iv(nfgcal) /= 0) go to 70
     iv(1) = 65
     go to 300

 70   dg1 = iv(dg)
  call vvmulp(n, v(dg1), g, d, -1)
  v(dgnorm) = v2norm(n, v(dg1))

  if (iv(cnvcod) /= 0) go to 290
  if (iv(mode) == 0) go to 250
!
!   allow first step to have scaled 2-norm at most v(lmax0)
!
  v(radius) = v(lmax0)

  iv(mode) = 0
!
!  main loop
!
!   print iteration summary, check iteration limit
!
 80   call itsum(d, g, iv, liv, lv, n, v, x)
 90   k = iv(niter)
  if (k < iv(mxiter)) go to 100
     iv(1) = 10
     go to 300
!
!   update radius
!
 100  iv(niter) = k + 1
  if(k>0)v(radius) = v(radfac) * v(dstnrm)
!
!   initialize for start of next iteration
!
  g01 = iv(g0)
  x01 = iv(x0)
  v(f0) = v(f)
  iv(irc) = 4
  iv(kagqt) = -1
!
!      copy x to x0, g to g0
!
  call vcopy(n, v(x01), x)
  call vcopy(n, v(g01), g)
!
!  Check STOPX and function evaluation limit
!
 110  if ( .not. stopx ( ) ) go to 130
     iv(1) = 11
     go to 140
!
!  Come here when restarting after func. eval. limit or STOPX.
!
 120  if (v(f) >= v(f0)) go to 130
     v(radfac) = one
     k = iv(niter)
     go to 100

 130  if (iv(nfcall) < iv(mxfcal)) go to 150
     iv(1) = 9
 140     if (v(f) >= v(f0)) go to 300
!
!  in case of STOPX or function evaluation limit with
!  improved v(f), evaluate the gradient at x.
!
          iv(cnvcod) = iv(1)
          go to 240
!
!  Compute candidate step
!
 150  step1 = iv(step)
  dg1 = iv(dg)
  nwtst1 = iv(nwtstp)
  if (iv(kagqt) >= 0) go to 160
     l = iv(lmat)
     call livmul(n, v(nwtst1), v(l), g)
     v(nreduc) = half * dotprd(n, v(nwtst1), v(nwtst1))
     call litvmu(n, v(nwtst1), v(l), v(nwtst1))
     call vvmulp(n, v(step1), v(nwtst1), d, 1)
     v(dst0) = v2norm(n, v(step1))
     call vvmulp(n, v(dg1), v(dg1), d, -1)
     call ltvmul(n, v(step1), v(l), v(dg1))
     v(gthg) = v2norm(n, v(step1))
     iv(kagqt) = 0
 160  call dbdog(v(dg1), lv, n, v(nwtst1), v(step1), v)
  if (iv(irc) == 6) go to 180
!
!   check whether evaluating f(x0 + step) looks worthwhile
!
  if (v(dstnrm) <= 0.0D+00) go to 180
  if (iv(irc) /= 5) go to 170
  if (v(radfac) <= one) go to 170
  if (v(preduc) <= onep2 * v(fdif)) go to 180
!
!  Compute f(x0 + step)
!
 170  x01 = iv(x0)
  step1 = iv(step)
  call vaxpy(n, x, one, v(step1), v(x01))
  iv(nfcall) = iv(nfcall) + 1
  iv(1) = 1
  iv(toobig) = 0
  return
!
!  Assess candidate step.
!
 180  x01 = iv(x0)
  v(reldx) = reldst(n, d, x, v(x01))
  call assst(iv, liv, lv, v)
  step1 = iv(step)
  lstgst = iv(stlstg)
  if (iv(restor) == 1) call vcopy(n, x, v(x01))
  if (iv(restor) == 2) call vcopy(n, v(lstgst), v(step1))
  if (iv(restor) /= 3) go to 190
     call vcopy(n, v(step1), v(lstgst))
     call vaxpy(n, x, one, v(step1), v(x01))
     v(reldx) = reldst(n, d, x, v(x01))

 190  k = iv(irc)
  go to (200,230,230,230,200,210,220,220,220,220,220,220,280,250), k
!
!      recompute step with changed radius
!
 200     v(radius) = v(radfac) * v(dstnrm)
     go to 110
!
!   compute step of length v(lmaxs) for singular convergence test.
!
 210  v(radius) = v(lmaxs)
  go to 150
!
!   convergence or false convergence
!
 220  iv(cnvcod) = k - 4
  if (v(f) >= v(f0)) go to 290
     if (iv(xirc) == 14) go to 290
          iv(xirc) = 14
!
!  Process acceptable step.
!
 230  if (iv(irc) /= 3) go to 240
     step1 = iv(step)
     temp1 = iv(stlstg)
!
!      set  temp1 = hessian * step  for use in gradient tests
!
     l = iv(lmat)
     call ltvmul(n, v(temp1), v(l), v(step1))
     call lvmul(n, v(temp1), v(l), v(temp1))
!
!   compute gradient
!
 240  iv(ngcall) = iv(ngcall) + 1
  iv(1) = 2
  return
!
!   initializations -- g0 = g - g0, etc.
!
 250  g01 = iv(g0)
  call vaxpy(n, v(g01), negone, v(g01), g)
  step1 = iv(step)
  temp1 = iv(stlstg)
  if (iv(irc) /= 3) go to 270
!
!   set v(radfac) by gradient tests
!
!  Set  temp1 = diag(d)**-1 * (hessian*step + (g(x0)-g(x)))
!
     call vaxpy(n, v(temp1), negone, v(g01), v(temp1))
     call vvmulp(n, v(temp1), v(temp1), d, -1)
!
!  Do gradient tests
!
     if (v2norm(n, v(temp1)) <= v(dgnorm) * v(tuner4)) then
       go to 260
     end if

     if (dotprd(n, g, v(step1)) >= v(gtstep) * v(tuner5))  then
       go to 270
     end if

 260               v(radfac) = v(incfac)
!
!   update h, loop
!
 270  w = iv(nwtstp)
  z = iv(x0)
  l = iv(lmat)
  call wzbfgs(v(l), n, v(step1), v(w), v(g01), v(z))
!
!  Use the n-vectors starting at v(step1) and v(g01) for scratch.
!
  call lupdat(v(temp1), v(step1), v(l), v(g01), v(l), n, v(w), v(z))
  iv(1) = 2
  go to 80
!
!   misc. details
!
!   bad parameters to assess
!
 280  iv(1) = 64
  go to 300
!
!  Print summary of final iteration and other requested items
!
 290  iv(1) = iv(cnvcod)
  iv(cnvcod) = 0
 300  call itsum(d, g, iv, liv, lv, n, v, x)

  return
end

subroutine parck ( alg, d, iv, liv, lv, n, v )
!
!  PARCK checks parameters, prints changed values.
!
!  Discussion:
!
!    alg = 1 for regression, alg = 2 for general unconstrained opt.
!
  integer alg, liv, lv, n
  integer iv(liv)
  double precision d(n), v(lv)
  double precision rmdcon
  integer max0
  integer i, ii, iv1, j, k, l, m, miv1, miv2, ndfalt, parsv1, pu
  integer ijmp, jlim(2), miniv(2), ndflt(2)
  character*1 varnm(2), sh(2)
  character*4 cngd(3), dflt(3), vn(2,34), which(3)
  double precision big, machep, tiny, vk, vm(34), vx(34)
  integer algsav, dinit, dtype, dtype0, epslon, inits, ivneed
  integer lastiv, lastv, lmat, nextiv, nextv, nvdflt, oldn
  integer parprt, parsav, perm, prunit, vneed

  parameter (algsav=51, dinit=38, dtype=16, dtype0=54, epslon=19 )
  parameter ( inits=25, ivneed=3, lastiv=44, lastv=45, lmat=42 )
  parameter ( nextiv=46, nextv=47, nvdflt=50, oldn=38, parprt=20 )
  parameter ( parsav=49, perm=58, prunit=21, vneed=4)
  save big, machep, tiny

  data big/0.d+0/, machep/-1.d+0/, tiny/1.d+0/

     data vn(1,1),vn(2,1)/'epsl','on..'/
     data vn(1,2),vn(2,2)/'phmn','fc..'/
     data vn(1,3),vn(2,3)/'phmx','fc..'/
     data vn(1,4),vn(2,4)/'decf','ac..'/
     data vn(1,5),vn(2,5)/'incf','ac..'/
     data vn(1,6),vn(2,6)/'rdfc','mn..'/
     data vn(1,7),vn(2,7)/'rdfc','mx..'/
     data vn(1,8),vn(2,8)/'tune','r1..'/
     data vn(1,9),vn(2,9)/'tune','r2..'/
     data vn(1,10),vn(2,10)/'tune','r3..'/
     data vn(1,11),vn(2,11)/'tune','r4..'/
     data vn(1,12),vn(2,12)/'tune','r5..'/
     data vn(1,13),vn(2,13)/'afct','ol..'/
     data vn(1,14),vn(2,14)/'rfct','ol..'/
     data vn(1,15),vn(2,15)/'xcto','l...'/
     data vn(1,16),vn(2,16)/'xfto','l...'/
     data vn(1,17),vn(2,17)/'lmax','0...'/
     data vn(1,18),vn(2,18)/'lmax','s...'/
     data vn(1,19),vn(2,19)/'scto','l...'/
     data vn(1,20),vn(2,20)/'dini','t...'/
     data vn(1,21),vn(2,21)/'dtin','it..'/
     data vn(1,22),vn(2,22)/'d0in','it..'/
     data vn(1,23),vn(2,23)/'dfac','....'/
     data vn(1,24),vn(2,24)/'dltf','dc..'/
     data vn(1,25),vn(2,25)/'dltf','dj..'/
     data vn(1,26),vn(2,26)/'delt','a0..'/
     data vn(1,27),vn(2,27)/'fuzz','....'/
     data vn(1,28),vn(2,28)/'rlim','it..'/
     data vn(1,29),vn(2,29)/'cosm','in..'/
     data vn(1,30),vn(2,30)/'hube','rc..'/
     data vn(1,31),vn(2,31)/'rspt','ol..'/
     data vn(1,32),vn(2,32)/'sigm','in..'/
     data vn(1,33),vn(2,33)/'eta0','....'/
     data vn(1,34),vn(2,34)/'bias','....'/

  data vm(1)/1.0d-3/, vm(2)/-0.99d+0/, vm(3)/1.0d-3/, vm(4)/1.0d-2/
  data vm(5)/1.2d+0/, vm(6)/1.d-2/, vm(7)/1.2d+0/, vm(8)/0.d+0/
  data vm(9)/0.d+0/, vm(10)/1.d-3/, vm(11)/-1.d+0/, vm(13)/0.d+0/
  data vm(15)/0.d+0/, vm(16)/0.d+0/, vm(19)/0.d+0/, vm(20)/-10.d+0/
  data vm(21)/0.d+0/, vm(22)/0.d+0/, vm(23)/0.d+0/, vm(27)/1.01d+0/
  data vm(28)/1.d+10/, vm(30)/0.d+0/, vm(31)/0.d+0/, vm(32)/0.d+0/
  data vm(34)/0.d+0/

  data vx(1)/0.9d+0/, vx(2)/-1.d-3/, vx(3)/1.d+1/, vx(4)/0.8d+0/
  data vx(5)/1.d+2/, vx(6)/0.8d+0/, vx(7)/1.d+2/, vx(8)/0.5d+0/
  data vx(9)/0.5d+0/, vx(10)/1.d+0/, vx(11)/1.d+0/, vx(14)/0.1d+0/
  data vx(15)/1.d+0/, vx(16)/1.d+0/, vx(19)/1.d+0/, vx(23)/1.d+0/
  data vx(24)/1.d+0/, vx(25)/1.d+0/, vx(26)/1.d+0/, vx(27)/1.d+10/
  data vx(29)/1.d+0/, vx(31)/1.d+0/, vx(32)/1.d+0/, vx(33)/1.d+0/
  data vx(34)/1.d+0/

  data varnm(1)/'p'/, varnm(2)/'n'/, sh(1)/'s'/, sh(2)/'h'/
  data cngd(1),cngd(2),cngd(3)/'---c','hang','ed v'/
  data dflt(1),dflt(2),dflt(3)/'nond','efau','lt v'/
  data ijmp/33/, jlim(1)/0/, jlim(2)/24/, ndflt(1)/32/, ndflt(2)/25/
  data miniv(1)/80/, miniv(2)/59/

  pu = 0
  if (prunit <= liv) pu = iv(prunit)
  if (alg < 1 .or. alg > 2) go to 340
  if (iv(1) == 0) call deflt(alg, iv, liv, lv, v)
  iv1 = iv(1)
  if (iv1 /= 13 .and. iv1 /= 12) go to 10
  miv1 = miniv(alg)
  if (perm <= liv) miv1 = max0(miv1, iv(perm) - 1)
  if (ivneed <= liv) miv2 = miv1 + max0(iv(ivneed), 0)
  if (lastiv <= liv) iv(lastiv) = miv2
  if (liv < miv1) go to 300
  iv(ivneed) = 0
  iv(lastv) = max0(iv(vneed), 0) + iv(lmat) - 1
  iv(vneed) = 0
  if (liv < miv2) go to 300
  if (lv < iv(lastv)) go to 320
 10   if (alg == iv(algsav)) go to 30
  if (pu /= 0) write(pu,20) alg, iv(algsav)
 20 format(/39h the first parameter to deflt should be,i3,12h rather than,i3)
     iv(1) = 82
     return
 30   if (iv1 < 12 .or. iv1 > 14) go to 60
     if (n >= 1) go to 50
          iv(1) = 81
          if (pu == 0) then
            return
          end if
          write(pu,40) varnm(alg), n
 40           format(/8h /// bad,a1,2h =,i5)
          return
 50      if (iv1 /= 14) iv(nextiv) = iv(perm)
     if (iv1 /= 14) iv(nextv) = iv(lmat)
     if (iv1 == 13) then
       return
     end if
     k = iv(parsav) - epslon
     call vdflt(alg, lv-k, v(k+1))
     iv(dtype0) = 2 - alg
     iv(oldn) = n
     which(1) = dflt(1)
     which(2) = dflt(2)
     which(3) = dflt(3)
     go to 110
 60   if (n == iv(oldn)) go to 80
     iv(1) = 17
     if (pu == 0) then
       return
     end if
     write(pu,70) varnm(alg), iv(oldn), n
 70      format(/5h /// ,1a1,14h changed from ,i5,4h to ,i5)
     return

 80   if (iv1 <= 11 .and. iv1 >= 1) go to 100
     iv(1) = 80
     if (pu /= 0) write(pu,90) iv1
 90      format(/13h ///  iv(1) =,i5,28h should be between 0 and 14.)
     return

 100  which(1) = cngd(1)
  which(2) = cngd(2)
  which(3) = cngd(3)

 110  if (iv1 == 14) iv1 = 12
  if (big > tiny) go to 120
     tiny = rmdcon(1)
     machep = rmdcon(3)
     big = rmdcon(6)
     vm(12) = machep
     vx(12) = big
     vx(13) = big
     vm(14) = machep
     vm(17) = tiny
     vx(17) = big
     vm(18) = tiny
     vx(18) = big
     vx(20) = big
     vx(21) = big
     vx(22) = big
     vm(24) = machep
     vm(25) = machep
     vm(26) = machep
     vx(28) = rmdcon(5)
     vm(29) = machep
     vx(30) = big
     vm(33) = machep
 120  m = 0
  i = 1
  j = jlim(alg)
  k = epslon
  ndfalt = ndflt(alg)

  do l = 1, ndfalt
    vk = v(k)
    if (vk >= vm(i) .and. vk <= vx(i)) go to 140
      m = k
      if (pu /= 0) write(pu,130) vn(1,i), vn(2,i), k, vk,vm(i), vx(i)
 130  format(/6h ///  ,2a4,5h.. v(,i2,3h) =,d11.3,7h should, &
      11h be between,d11.3,4h and,d11.3)
 140  k = k + 1
     i = i + 1
     if (i == j) i = ijmp
  end do

  if (iv(nvdflt) == ndfalt) go to 170
     iv(1) = 51
     if (pu == 0) then
       return
     end if
     write(pu,160) iv(nvdflt), ndfalt
 160     format(/13h iv(nvdflt) =,i5,13h rather than ,i5)
     return
 170  if ((iv(dtype) > 0 .or. v(dinit) > 0.0D+00) .and. iv1 == 12) then
             go to 200
  end if

  do i = 1, n
     if (d(i) > 0.0D+00) go to 190
          m = 18
          if (pu /= 0) write(pu,180) i, d(i)
 180     format(/8h ///  d(,i3,3h) =,d11.3,19h should be positive)
 190     continue
  end do

 200  if (m == 0) go to 210
     iv(1) = m
     return

 210  if (pu == 0 .or. iv(parprt) == 0) then
        return
      end if
  if (iv1 /= 12 .or. iv(inits) == alg-1) go to 230
     m = 1
     write(pu,220) sh(alg), iv(inits)
 220 format(/22h nondefault values..../5h init,a1,14h      iv(25) =,i3)
 230  if (iv(dtype) == iv(dtype0)) go to 250
     if (m == 0) write(pu,260) which
     m = 1
     write(pu,240) iv(dtype)
 240     format(20h dtype      iv(16) =,i3)
 250  i = 1
  j = jlim(alg)
  k = epslon
  l = iv(parsav)
  ndfalt = ndflt(alg)

  do ii = 1, ndfalt
     if (v(k) == v(l)) go to 280
          if (m == 0) write(pu,260) which
 260          format(/1h ,3a4,9halues..../)
          m = 1
          write(pu,270) vn(1,i), vn(2,i), k, v(k)
 270          format(1x,2a4,5h.. v(,i2,3h) =,d15.7)
 280     k = k + 1
     l = l + 1
     i = i + 1
     if (i == j) i = ijmp
  end do

  iv(dtype0) = iv(dtype)
  parsv1 = iv(parsav)
  call vcopy(iv(nvdflt), v(parsv1), v(epslon))
  return

 300  iv(1) = 15
  if (pu == 0) then
    return
  end if
  write(pu,310) liv, miv2
 310  format(/10h /// liv =,i5,17h must be at least,i5)
  if (liv < miv1) then
    return
  end if
  if (lv < iv(lastv)) go to 320
  return

 320  iv(1) = 16
  if (pu == 0) then
    return
  end if
  write(pu,330) lv, iv(lastv)
 330  format(/9h /// lv =,i5,17h must be at least,i5)
  return

 340  iv(1) = 67
  if (pu == 0) then
    return
  end if
  write(pu,350) alg
 350  format(/10h /// alg =,i5,15h must be 1 or 2)

  return
end

subroutine itsum ( d, g, iv, liv, lv, p, v, x )

!
! ITSUM prints an iteration summary.
!
  integer liv
  integer lv
  integer p

  double precision d(p)
  double precision g(p)
  integer iv(liv)
  double precision v(lv)
  double precision x(p)
  integer alg, i, iv1, m, nf, ng, ol, pu
  character*4 model1(6), model2(6)
  double precision nreldf, oldf, preldf, reldf
  integer algsav, dstnrm, f, fdif, f0, needhd, nfcall, nfcov, ngcov
  integer ngcall, niter, nreduc, outlev, preduc, prntit, prunit
  integer reldx, solprt, statpr, stppar, sused, x0prt

  parameter (algsav=51, needhd=36, nfcall=6, nfcov=52, ngcall=30 )
  parameter ( ngcov=53, niter=31, outlev=19, prntit=39, prunit=21 )
  parameter ( solprt=22, statpr=23, sused=64, x0prt=24)
  parameter (dstnrm=2, f=10, f0=13, fdif=11, nreduc=6, preduc=7, reldx=17 )
  parameter ( stppar=5)

  data model1/'    ','    ','    ','    ','  g ','  s '/
  data model2/' g  ',' s  ','g-s ','s-g ','-s-g','-g-s'/

  pu = iv(prunit)

  if ( pu == 0 ) then
    return
  end if

  iv1 = iv(1)
  if (iv1 > 62) iv1 = iv1 - 51
  ol = iv(outlev)
  alg = iv(algsav)
  if (iv1 < 2 .or. iv1 > 15) go to 370
  if (iv1 >= 12) go to 120
  if (iv1 == 2 .and. iv(niter) == 0) go to 390
  if (ol == 0) go to 120
  if (iv1 >= 10 .and. iv(prntit) == 0) go to 120
  if (iv1 > 2) go to 10
     iv(prntit) = iv(prntit) + 1
     if (iv(prntit) < iabs(ol)) then
       return
     end if
 10   nf = iv(nfcall) - iabs(iv(nfcov))
  iv(prntit) = 0
  reldf = 0.0D+00
  preldf = 0.0D+00
  oldf = max (abs(v(f0)), abs(v(f)))
  if (oldf <= 0.0D+00) go to 20
     reldf = v(fdif) / oldf
     preldf = v(preduc) / oldf
 20   if (ol > 0) go to 60
!
!  print short summary line
!
     if (iv(needhd) == 1 .and. alg == 1) write(pu,30)
 30   format(/'   it   nf',6x,'f',7x,'reldf',3x,'preldf',3x,'reldx', &
     &  2x,'model  stppar')
     if (iv(needhd) == 1 .and. alg == 2) write(pu,40)
 40  format(/'    it   nf',7x,'f',8x,'reldf',4x,'preldf',4x,'reldx   stppar')
     iv(needhd) = 0
     if (alg == 2) go to 50
     m = iv(sused)
     write(pu,100) iv(niter), nf, v(f), reldf, preldf, v(reldx), &
       model1(m), model2(m), v(stppar)
     go to 120

 50  write(pu,110) iv(niter), nf, v(f), reldf, preldf, v(reldx), v(stppar)
     go to 120
!
!  print long summary line
!
 60   if (iv(needhd) == 1 .and. alg == 1) write(pu,70)
 70   format(/11h    it   nf,6x,1hf,7x,5hreldf,3x,6hpreldf,3x,5hreldx, &
        2x,13hmodel  stppar,2x,6hd*step,2x,7hnpreldf)
  if (iv(needhd) == 1 .and. alg == 2) write(pu,80)
 80   format(/11h    it   nf,7x,1hf,8x,5hreldf,4x,6hpreldf,4x,5hreldx, &
        3x,6hstppar,3x,6hd*step,3x,7hnpreldf)
  iv(needhd) = 0
  nreldf = 0.0D+00
  if (oldf > 0.0D+00) nreldf = v(nreduc) / oldf
  if (alg == 2) go to 90
  m = iv(sused)
  write(pu,100) iv(niter), nf, v(f), reldf, preldf, v(reldx), &
              model1(m), model2(m), v(stppar), v(dstnrm), nreldf
  go to 120

 90   write(pu,110) iv(niter), nf, v(f), reldf, preldf, &
             v(reldx), v(stppar), v(dstnrm), nreldf
 100  format(i6,i5,d10.3,2d9.2,d8.1,a3,a4,2d8.1,d9.2)
 110  format(i6,i5,d11.3,2d10.2,3d9.1,d10.2)

 120  if (iv(statpr) < 0) go to 430
  go to (999, 999, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, &
    330, 350, 520), iv1

 130  write(pu,140)
 140  format(/' x-convergence')
  go to 430

 150  write ( pu, 160 )
 160  format(/'relative function convergence')
  go to 430

 170  write(pu,180)
 180  format(/'x- and relative function convergence')
  go to 430

 190  write(pu,200)
 200  format(/'Absolute function convergence.')
  go to 430

 210  write(pu,220)
 220  format(/'Singular convergence.')
  go to 430

 230  write(pu,240)
 240  format(/'False convergence.')
  go to 430

 250  write(pu,260)
 260  format(/'Function evaluation limit.')
  go to 430

 270  write(pu,280)
 280  format(/'Iteration limit.')
  go to 430

 290  write(pu,300)
 300  format(/'STOPX')
  go to 430

 310  write(pu,320)
 320  format(/'Initial f(x) cannot be computed.')

  go to 390

 330  write(pu,340)
 340  format(/'Bad parameters to assess.')
  go to 999

 350  write(pu,360)
 360  format(/'Gradient could not be computed.')
  if (iv(niter) > 0) go to 480
  go to 390

 370  write(pu,380) iv(1)
 380  format(/'iv(1) =',i5)
  go to 999
!
!   initial call on itsum
!
 390  if (iv(x0prt) /= 0) write(pu,400) (i, x(i), d(i), i = 1, p)
 400  format(/23h     i     initial x(i),8x,4hd(i)//(1x,i5,d17.6,d14.3))
!     the following are to avoid undefined variables when the
!     function evaluation limit is 1...
!
  v(dstnrm) = 0.0D+00
  v(fdif) = 0.0D+00
  v(nreduc) = 0.0D+00
  v(preduc) = 0.0D+00
  v(reldx)  = 0.0D+00
  if (iv1 >= 12) go to 999
  iv(needhd) = 0
  iv(prntit) = 0
  if (ol == 0) go to 999
  if (ol < 0 .and. alg == 1) write(pu,30)
  if (ol < 0 .and. alg == 2) write(pu,40)
  if (ol > 0 .and. alg == 1) write(pu,70)
  if (ol > 0 .and. alg == 2) write(pu,80)
  if (alg == 1) write(pu,410) v(f)
  if (alg == 2) write(pu,420) v(f)
 410  format(/11h     0    1,d10.3)
!365  format(/11h     0    1,e11.3)
 420  format(/11h     0    1,d11.3)
  go to 999
!
!  Print various information requested on solution
!
 430  iv(needhd) = 1
  if (iv(statpr) == 0) go to 480
     oldf = max (abs(v(f0)), abs(v(f)))
     preldf = 0.0D+00
     nreldf = 0.0D+00
     if (oldf <= 0.0D+00) go to 440
          preldf = v(preduc) / oldf
          nreldf = v(nreduc) / oldf
 440     nf = iv(nfcall) - iv(nfcov)
     ng = iv(ngcall) - iv(ngcov)
     write(pu,450) v(f), v(reldx), nf, ng, preldf, nreldf
 450  format(/9h function,d17.6,8h   reldx,d17.3/12h func. evals, &
    i8,9x,11hgrad. evals,i8/7h preldf,d16.3,6x,7hnpreldf,d15.3)

     if (iv(nfcov) > 0) write(pu,460) iv(nfcov)
 460     format(/1x,i4,50h extra func. evals for covariance and diagnostics.)
     if (iv(ngcov) > 0) write(pu,470) iv(ngcov)
 470     format(1x,i4,50h extra grad. evals for covariance and diagnostics.)

 480  if (iv(solprt) == 0) go to 999
     iv(needhd) = 1
     write(pu,490)
 490  format(/22h     i      final x(i),8x,4hd(i),10x,4hg(i)/)
     do i = 1, p
          write(pu,510) i, x(i), d(i), g(i)
     end do
 510     format(1x,i5,d16.6,2d14.3)
  go to 999

 520  write(pu,530)
 530  format(/'Inconsistent dimensions.')
 999  continue

  return
end

function v2norm ( p, x )
!
!! V2NORM returns the 2-norm of the p-vector X.
!
!  Discussion:
!
!    The routine tries to avoid underflow.
!
!  Parameters:
!
  integer p

  double precision x(p)
  integer i, j
  double precision r, scale
  double precision, save :: sqteta = 0.0D+00
  double precision t, xi
  double precision rmdcon
  double precision v2norm

  v2norm = 0.0D+00

  if (p <= 0 ) then
    return
  end if

  if ( all ( x(1:p) == 0.0D+00 ) ) then
    return
  end if

  scale = 0.0D+00
  do i = 1, p
    if ( x(i) /= 0.0D+00 ) then
      scale = abs(x(i))
      exit
    end if
  end do

  if ( scale == 0.0D+00 ) then
    return
  end if

  if ( p <= i ) then
    v2norm = scale
    return
  end if

  t = 1.0D+00
  if ( sqteta == 0.0D+00 ) then
    sqteta = rmdcon(2)
  end if
!
!  sqteta is (slightly larger than) the square root of the
!  smallest positive floating point number on the machine.
!  the tests involving sqteta are done to prevent underflows.
!
  j = i + 1
  do i = j, p
    xi = abs(x(i))
    if (xi <= scale) then
      r = xi / scale
      if (r > sqteta) t = t + r*r
    else
      r = scale / xi
      if (r <= sqteta) r = 0.0D+00
      t = 1.0D+00  +  t * r*r
      scale = xi
    end if
  end do

  v2norm = scale * sqrt(t)

  return
end

function rmdcon ( k )
!
!! RMDCON returns machine dependent constants.
!
!  Discussion:
!
!    Comments below contain data statements for various machines.
!    To convert to another machine, place a c in column 1 of the
!    data statement line(s) that correspond to the current machine
!    and remove the c from column 1 of the data statement line(s)
!    that correspond to the new machine.
!
!    the constant returned depends on k...
!
!         k = 1... smallest pos. eta such that -eta exists.
!         k = 2... square root of eta.
!         k = 3... unit roundoff = smallest pos. no. machep such
!                  that 1 + machep > 1 .and. 1 - machep < 1.
!         k = 4... square root of machep.
!         k = 5... square root of big (see k = 6).
!         k = 6... largest machine no. big such that -big exists.
!
  integer k
  double precision rmdcon
  double precision big, eta, machep
  integer bigi(4), etai(4), machei(4)
  equivalence (big,bigi(1)), (eta,etai(1)), (machep,machei(1))
!
!  ibm 360, ibm 370, or xerox
!
!     data big/z7fffffffffffffff/, eta/z0010000000000000/,
!    1     machep/z3410000000000000/
!
!  data general
!
!     data big/0.7237005577d+76/, eta/0.5397605347d-78/,
!    1     machep/2.22044605d-16/
!
!  dec 11
!
!     data big/1.7d+38/, eta/2.938735878d-39/, machep/2.775557562d-17/
!
!  hp3000
!
!     data big/1.157920892d+77/, eta/8.636168556d-78/,
!    1     machep/5.551115124d-17/
!
!  honeywell
!
!     data big/1.69d+38/, eta/5.9d-39/, machep/2.1680435d-19/
!
!  dec10
!
!     data big/"377777100000000000000000/,
!    1     eta/"002400400000000000000000/,
!    2     machep/"104400000000000000000000/
!
!  burroughs
!
!     data big/o0777777777777777,o7777777777777777/,
!    1     eta/o1771000000000000,o7770000000000000/,
!    2     machep/o1451000000000000,o0000000000000000/
!
!  control data
!
!     data big/37767777777777777777b,37167777777777777777b/,
!    1     eta/00014000000000000000b,00000000000000000000b/,
!    2     machep/15614000000000000000b,15010000000000000000b/
!
!  prime
!
!     data big/1.0d+9786/, eta/1.0d-9860/, machep/1.4210855d-14/
!
!  univac
!
!     data big/8.988d+307/, eta/1.2d-308/, machep/1.734723476d-18/
!
!  vax
!
  data big/1.7d+38/, eta/2.939d-39/, machep/1.3877788d-17/
!
!  cray 1
!
!     data bigi(1)/577767777777777777777b/,
!    1     bigi(2)/000007777777777777776b/,
!    2     etai(1)/200004000000000000000b/,
!    3     etai(2)/000000000000000000000b/,
!    4     machei(1)/377224000000000000000b/,
!    5     machei(2)/000000000000000000000b/
!
!  port library -- requires more than just a data statement...
!
!     external d1mach
!     double precision d1mach, zero
!     data big/0.d+0/, eta/0.d+0/, machep/0.d+0/, zero/0.d+0/
!     if (big > 0.0D+00) go to 1
!        big = d1mach(2)
!        eta = d1mach(1)
!        machep = d1mach(4)
!1    continue
!
! end of port
!
!  body -
!
  go to (10, 20, 30, 40, 50, 60), k

 10   rmdcon = eta
  return

 20   rmdcon = sqrt(256.d+0*eta)/16.d+0
  return

 30   rmdcon = machep
  return

 40   rmdcon = sqrt(machep)
  return

 50   rmdcon = sqrt(big/256.d+0)*16.d+0
  return

 60   rmdcon = big

  return
end

subroutine vcopy ( p, y, x )
!
!! VCOPY sets y = x.
!
!  Discussion:
!
!    x and y are p-vectors
!
  implicit none

  integer p
  double precision x(p)
  double precision y(p)

  y(1:p) = x(1:p)

  return
end

subroutine vscopy ( p, y, s )
!
!! VSCOPY sets the vector Y to scalar S.
!
  implicit none
  integer p
  double precision s
  double precision y(p)

  y(1:p) = s

  return
end

subroutine vvmulp ( n, x, y, z, k )
!
!! VVMULP sets x(i) = y(i) * z(i)**k, 1 <= i <= n (for k = 1 or -1)
!
  implicit none
  integer n
  integer k
  double precision x(n)
  double precision y(n)
  double precision z(n)

  if ( k < 0 ) then
    x(1:n) = y(1:n) / z(1:n)
  else
    x(1:n) = y(1:n) * z(1:n)
  end if

  return
end

subroutine lvmul ( n, x, l, y )
!
!! LVMUL computes x = L * y.
!
!  Discussion:
!
!    L  is an  n x n  lower triangular matrix stored compactly by rows.
!    x and y may occupy the same storage.
!
  integer n

  double precision x(n), l(*), y(n)
!     dimension l(n*(n+1)/2)
  integer i, ii, ij, i0, j, np1
  double precision t

  np1 = n + 1
  i0 = n*(n+1)/2

  do ii = 1, n
    i = np1 - ii
    i0 = i0 - i
    t = 0.0D+00
    do j = 1, i
      ij = i0 + j
      t = t + l(ij)*y(j)
    end do
    x(i) = t
  end do

  return
end

subroutine livmul ( n, x, l, y )
!
!! LIVMUL solves L * x = y.
!
!  Discussion:
!
!    L is an  n x n  lower triangular
!    matrix stored compactly by rows.  x and y may occupy the same
!    storage.
!
  integer n

  double precision x(n), l(*), y(n)
  external dotprd
  double precision dotprd
  integer i, j, k
  double precision t

  do k = 1, n
    if (y(k) /= 0.0D+00 ) go to 20
    x(k) = 0.0D+00
  end do

  return

20 continue

  j = k*(k+1)/2
  x(k) = y(k) / l(j)

  if (k >= n) then
    return
  end if

  k = k + 1

  do i = k, n
     t = dotprd(i-1, l(j+1), x)
     j = j + i
     x(i) = (y(i) - t)/l(j)
  end do

  return
end

subroutine litvmu ( n, x, l, y )
!
!! LITVMU solves L' * x = y.
!
!  Discussion:
!
!    L is an  n x n  lower triangular
!    matrix stored compactly by rows.  x and y may occupy the same
!    storage.
!
  integer n

  double precision l(*)
  double precision x(n)
  double precision y(n)
  integer i, ii, ij, i0, j
  double precision xi

  x(1:n) = y(1:n)

  i0 = n*(n+1)/2

  do ii = 1, n
    i = n+1 - ii
    xi = x(i)/l(i0)
    x(i) = xi
    if ( i <= 1 ) then
      exit
    end if
    i0 = i0 - i
    if ( xi /= 0.0D+00 ) then
      do j = 1, i-1
        ij = i0 + j
        x(j) = x(j) - xi*l(ij)
      end do
    end if
  end do

  return
end

subroutine ltvmul ( n, x, l, y )
!
!! LTVMUL computes  x = (l**t)*y.
!
!  Discussion:
!
!    L is an  n x n  lower triangular matrix stored compactly by rows.
!    x and y may occupy the same storage.
!
  integer n
  double precision x(n), l(*), y(n)
!     dimension l(n*(n+1)/2)
  integer i, ij, i0, j
  double precision yi

  i0 = 0
  do i = 1, n
    yi = y(i)
    x(i) = 0.0D+00
    do j = 1, i
      ij = i0 + j
      x(j) = x(j) + yi * l(ij)
    end do
    i0 = i0 + i
  end do

  return
end

subroutine lupdat ( beta, gamma, l, lambda, lplus, n, w, z )
!
!! LUPDAT computes lplus = secant update of L.
!
!  Discussion:
!
!    this routine updates the cholesky factor  l  of a symmetric
!    positive definite matrix to which a secant update is being
!    applied -- it computes a cholesky factor  lplus  of
!    l * (i + z*w**t) * (i + w*z**t) * l**t.  it is assumed that  w
!    and  z  have been chosen so that the updated matrix is strictly
!    positive definite.
!
!    this code uses recurrence 3 of ref. 1 (with d(j) = 1 for all j)
!    to compute  lplus  of the form  l * (i + z*w**t) * q,  where  q
!    is an orthogonal matrix that makes the result lower triangular.
!    lplus may have some negative diagonal elements.
!
!  Reference:
!
!    D Goldfarb,
!    Factorized Variable Metric Methods for Unconstrained Optimization,
!    Mathematics of Computation,
!    Volume 30, pages 796-811, 1976.
!
!  Parameters:
!
!   beta = scratch vector.
!  gamma = scratch vector.
!      l (input) lower triangular matrix, stored rowwise.
! lambda = scratch vector.
!  lplus (output) lower triangular matrix, stored rowwise, which may
!             occupy the same storage as  l.
!      n (input) length of vector parameters and order of matrices.
!      w (input, destroyed on output) right singular vector of rank 1
!             correction to  l.
!      z (input, destroyed on output) left singular vector of rank 1
!             correction to  l.
!
  integer n
  double precision beta(n), gamma(n), l(*), lambda(n), lplus(*), w(n), z(n)
!     dimension l(n*(n+1)/2), lplus(n*(n+1)/2)
!
  integer i, ij, j, jj, jp1, k, nm1
  integer np1
  double precision a, b, bj, eta, gj, lj, lij, ljj, nu, s, theta, wj, zj
  double precision one

  parameter (one=1.d+0 )

  nu = one
  eta = 0.0D+00
  if (n <= 1) go to 30
  nm1 = n - 1
!
!  temporarily store s(j) = sum over k = j+1 to n of w(k)**2 in
!  lambda(j).
!
  s = 0.0D+00
  do i = 1, nm1
     j = n - i
     s = s + w(j+1)**2
     lambda(j) = s
  end do
!
!  compute lambda, gamma, and beta by goldfarb*s recurrence 3.
!
  do 20 j = 1, nm1
     wj = w(j)
     a = nu*z(j) - eta*wj
     theta = one + a*wj
     s = a*lambda(j)
     lj = sqrt(theta**2 + a*s)
     if (theta > 0.0D+00) lj = -lj
     lambda(j) = lj
     b = theta*wj + s
     gamma(j) = b * nu / lj
     beta(j) = (a - b*eta) / lj
     nu = -nu / lj
     eta = -(eta + (a**2)/(theta - lj)) / lj
 20      continue
 30   lambda(n) = one + (nu*z(n) - eta*w(n))*w(n)
!
!  update l, gradually overwriting  w  and  z  with  l*w  and  l*z.
!
  np1 = n + 1
  jj = n * (n + 1) / 2

  do k = 1, n

     j = np1 - k
     lj = lambda(j)
     ljj = l(jj)
     lplus(jj) = lj * ljj
     wj = w(j)
     w(j) = ljj * wj
     zj = z(j)
     z(j) = ljj * zj
     if (k == 1) go to 50
     bj = beta(j)
     gj = gamma(j)
     ij = jj + j
     jp1 = j + 1

     do i = jp1, n
          lij = l(ij)
          lplus(ij) = lj*lij + bj*w(i) + gj*z(i)
          w(i) = w(i) + lij*wj
          z(i) = z(i) + lij*zj
          ij = ij + i
     end do

 50      jj = jj - j

  end do

  return
end

subroutine dbdog ( dig, lv, n, nwtstp, step, v )
!
!! DBDOG: compute a double dogleg step.
!
!  Discussion:
!
!    This subroutine computes a candidate step (for use in an
!    unconstrained minimization code) by the double dogleg algorithm of
!    dennis and mei (ref. 1), which is a variation on powell*s dogleg
!    scheme (ref. 2, p. 95).
!
!    let  g  and  h  be the current gradient and hessian approxima-
!    tion respectively and let d be the current scale vector.  this
!    routine assumes dig = diag(d)**-2 * g  and  nwtstp = h**-1 * g.
!    the step computed is the same one would get by replacing g and h
!    by  diag(d)**-1 * g  and  diag(d)**-1 * h * diag(d)**-1,
!    computing step, and translating step back to the original
!    variables, i.e., premultiplying it by diag(d)**-1.
!
!  Reference:
!
!    John Dennis, Howell Mei,
!    Two New Unconstrained Optimization Algorithms Which Use
!    Function and Gradient Values,
!    Journal of Optimization Theory and Applications,
!    Volume 28, pages 453-482, 1979.
!
!    M J D Powell,
!    A Hybrid Method for Non-linear Equations,
!    in Numerical Methods for Non-linear Equations,
!    edited by Philip Rabinowitz,
!    Gordon and Breach, London, 1970.
!
!  Parameters:
!
!    dig (input) diag(d)**-2 * g -- see algorithm notes.
!      g (input) the current gradient vector.
!     lv (input) length of v.
!      n (input) number of components in  dig, g, nwtstp,  and  step.
! nwtstp (input) negative newton step -- see algorithm notes.
!   step (output) the computed step.
!      v (i/o) values array, the following components of which are
!             used here...
! v(bias)   (input) bias for relaxed newton step, which is v(bias) of
!             the way from the full newton to the fully relaxed newton
!             step.  recommended value = 0.8 .
! v(dgnorm) (input) 2-norm of diag(d)**-1 * g -- see algorithm notes.
! v(dstnrm) (output) 2-norm of diag(d) * step, which is v(radius)
!             unless v(stppar) = 0 -- see algorithm notes.
! v(dst0) (input) 2-norm of diag(d) * nwtstp -- see algorithm notes.
! v(grdfac) (output) the coefficient of  dig  in the step returned --
!             step(i) = v(grdfac)*dig(i) + v(nwtfac)*nwtstp(i).
! v(gthg)   (input) square-root of (dig**t) * (hessian) * dig -- see
!             algorithm notes.
! v(gtstep) (output) inner product between g and step.
! v(nreduc) (output) function reduction predicted for the full newton
!             step.
! v(nwtfac) (output) the coefficient of  nwtstp  in the step returned --
!             see v(grdfac) above.
! v(preduc) (output) function reduction predicted for the step returned.
! v(radius) (input) the trust region radius.  d times the step returned
!             has 2-norm v(radius) unless v(stppar) = 0.
! v(stppar) (output) code telling how step was computed... 0 means a
!             full newton step.  between 0 and 1 means v(stppar) of the
!             way from the newton to the relaxed newton step.  between
!             1 and 2 means a true double dogleg step, v(stppar) - 1 of
!             the way from the relaxed newton to the Cauchy step.
!             greater than 2 means 1 / (v(stppar) - 1) times the Cauchy
!             step.
!
  integer lv
  integer n

  double precision dig(n), nwtstp(n), step(n), v(lv)
  external dotprd, v2norm
  double precision dotprd, v2norm
  double precision cfact, cnorm, ctrnwt, ghinvg, femnsq, gnorm
  double precision nwtnrm, relax, rlambd, t, t1, t2
  double precision half, two
  integer bias, dgnorm, dstnrm, dst0, grdfac, gthg, gtstep
  integer nreduc, nwtfac, preduc, radius, stppar
  parameter (half=0.5d+0, two=2.d+0)
  parameter (bias=43, dgnorm=1, dstnrm=2, dst0=3, grdfac=45 )
  parameter ( gthg=44, gtstep=4, nreduc=6, nwtfac=46, preduc=7 )
  parameter ( radius=8, stppar=5)

  nwtnrm = v(dst0)
  rlambd = 1.0D+00
  if (nwtnrm > 0.0D+00 ) rlambd = v(radius) / nwtnrm
  gnorm = v(dgnorm)
  ghinvg = two * v(nreduc)
  v(grdfac) = 0.0D+00
  v(nwtfac) = 0.0D+00
  if (rlambd < 1.0D+00 ) go to 30
!
!  The Newton step is inside the trust region.
!
     v(stppar) = 0.0D+00
     v(dstnrm) = nwtnrm
     v(gtstep) = -ghinvg
     v(preduc) = v(nreduc)
     v(nwtfac) = -1.0D+00
     step(1:n) = -nwtstp(1:n)
     return

 30   v(dstnrm) = v(radius)
  cfact = (gnorm / v(gthg))**2
!
!  Cauchy step = -cfact * g.
!
  cnorm = gnorm * cfact
  relax = 1.0D+00 - v(bias) * ( 1.0D+00 - gnorm*cnorm/ghinvg)
  if (rlambd < relax) go to 50
!
!  Step is between relaxed Newton and full Newton steps.
!
     v(stppar) = 1.0D+00 -  (rlambd - relax) / ( 1.0D+00 - relax)
     t = -rlambd
     v(gtstep) = t * ghinvg
     v(preduc) = rlambd * ( 1.0D+00 - half*rlambd) * ghinvg
     v(nwtfac) = t
     step(1:n) = t * nwtstp(1:n)
     return

 50   if (cnorm < v(radius)) go to 70
!
!  The Cauchy step lies outside the trust region --
!  step = scaled Cauchy step.
!
     t = -v(radius) / gnorm
     v(grdfac) = t
     v(stppar) = 1.0D+00  +  cnorm / v(radius)
     v(gtstep) = -v(radius) * gnorm
  v(preduc) = v(radius)*(gnorm - half*v(radius)*(v(gthg)/gnorm)**2)
     step(1:n) = t * dig(1:n)
     return
!
!  Compute dogleg step between Cauchy and relaxed Newton
!  femur = relaxed newton step minus Cauchy step.
!
 70   ctrnwt = cfact * relax * ghinvg / gnorm
!
!  ctrnwt = inner product of Cauchy and relaxed Newton steps,
!  scaled by gnorm**-1.
!
  t1 = ctrnwt - gnorm*cfact**2
!
!  t1 = inner prod. of femur and Cauchy step, scaled by gnorm**-1.
!
  t2 = v(radius)*(v(radius)/gnorm) - gnorm*cfact**2
  t = relax * nwtnrm
  femnsq = (t/gnorm)*t - ctrnwt - t1
!
!  femnsq = square of 2-norm of femur, scaled by gnorm**-1.
!
  t = t2 / (t1 + sqrt(t1**2 + femnsq*t2))
!
!  Dogleg step  =  Cauchy step  +  t * femur.
!
  t1 = (t - 1.0D+00 ) * cfact
  v(grdfac) = t1
  t2 = -t * relax
  v(nwtfac) = t2
  v(stppar) = two - t
  v(gtstep) = t1*gnorm**2 + t2*ghinvg
  v(preduc) = -t1*gnorm * ((t2 + 1.0D+00 )*gnorm) &
                  - t2 * ( 1.0D+00 + half*t2)*ghinvg &
                   - half * (v(gthg)*t1)**2

  step(1:n) = t1 * dig(1:n) + t2 * nwtstp(1:n)

  return
end

function dotprd ( p, x, y )
!
!! DOTPRD returns the inner product of vectors X and Y.
!
  integer p

  double precision dotprd
  integer i
  double precision rmdcon
  double precision, save :: sqteta = 0.0D+00
  double precision t
  double precision x(p)
  double precision y(p)

  dotprd = 0.0D+00

  if ( sqteta == 0.0D+00 ) then
    sqteta = rmdcon(2)
  end if

  do i = 1, p

    t = max ( abs ( x(i) ), abs ( y(i) ) )
    if ( t > 1.0D+00 ) go to 10
    if (t < sqteta) go to 20
    t = (x(i)/sqteta)*y(i)
    if (abs(t) < sqteta) go to 20
 10   dotprd = dotprd + x(i)*y(i)

 20 continue

  end do

  return
end

function stopx ( )
!
!! STOPX checks to see if the BREAK key has been pressed.
!
!  Discussion:
!
!     this function may serve as the stopx (asynchronous interruption)
!     function for the nl2sol (nonlinear least-squares) package at
!     those installations which do not wish to implement a
!     dynamic stopx.
!
!     at installations where the nl2sol system is used
!     interactively, this dummy stopx should be replaced by a
!     function that returns .true. if and only if the interrupt
!     (break) key has been pressed since the last call on stopx.
!
  logical stopx

  stopx = .false.

  return
end

subroutine vaxpy ( p, w, a, x, y )
!
!! VAXPY sets w = a*x + y.
!
!  Discussion:
!
!    w, x, y = p-vectors, a = scalar
!
!  Parameters:
!
  implicit none

  integer p

  double precision a
  double precision w(p)
  double precision x(p)
  double precision y(p)

  w(1:p) = a * x(1:p) + y(1:p)

  return
end

function reldst ( p, d, x, x0 )
!
!! RELDST computes the relative difference between X and X0.
!
  integer p

  double precision reldst
  double precision d(p), x(p), x0(p)
  integer i
  double precision emax, t, xmax

  emax = 0.0D+00
  xmax = 0.0D+00

  do i = 1, p
    t = abs(d(i) * (x(i) - x0(i)))
    if (emax < t) emax = t
    t = d(i) * (abs(x(i)) + abs(x0(i)))
    if (xmax < t) xmax = t
  end do

  reldst = 0.0D+00
  if ( xmax > 0.0D+00 ) reldst = emax / xmax

  return
end

subroutine assst ( iv, liv, lv, v )
!
!! ASSST assesses a candidate step.
!
!  Discussion:
!
!    This subroutine is called by an unconstrained minimization
!    routine to assess the next candidate step.  it may recommend one
!    of several courses of action, such as accepting the step, recom-
!    puting it using the same or a new quadratic model, or halting due
!    to convergence or false convergence.  See the return code listing
!    below.
!
!  Reference:
!
!    John Dennis, David Gay, Roy Welsch,
!    An Adaptive Nonlinear Least-squares Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, 1981.
!
!    M J D Powell,
!    A Fortran Subroutine for Solving Systems of Nonlinear Algebraic Equations,
!    in Numerical Methods for Nonlinear Algebraic Equations,
!    edited by Philip Rabinowitz,
!    Gordon and Breach, London, 1970.
!
!  Parameters:
!
!  iv (i/o) integer parameter and scratch vector -- see description
!             below of iv values referenced.
!
! liv (in)  length of iv array.
!
!  lv (in)  length of v array.
!
!   v (i/o) real parameter and scratch vector -- see description
!             below of v values referenced.
!
!   iv values referenced
!
!    iv(irc) (i/o) on input for the first step tried in a new iteration,
!             iv(irc) should be set to 3 or 4 (the value to which it is
!             set when step is definitely to be accepted).  on input
!             after step has been recomputed, iv(irc) should be
!             unchanged since the previous return of assst.
!                on output, iv(irc) is a return code having one of the
!             following values...
!                  1 = switch models or try smaller step.
!                  2 = switch models or accept step.
!                  3 = accept step and determine v(radfac) by gradient
!                       tests.
!                  4 = accept step, v(radfac) has been determined.
!                  5 = recompute step (using the same model).
!                  6 = recompute step with radius = v(lmaxs) but do not
!                       evaulate the objective function.
!                  7 = x-convergence (see v(xctol)).
!                  8 = relative function convergence (see v(rfctol)).
!                  9 = both x- and relative function convergence.
!                 10 = absolute function convergence (see v(afctol)).
!                 11 = singular convergence (see v(lmaxs)).
!                 12 = false convergence (see v(xftol)).
!                 13 = iv(irc) was out of range on input.
!             return code i has precdence over i+1 for i = 9, 10, 11.
! iv(mlstgd) (i/o) saved value of iv(model).
!  iv(model) (i/o) on input, iv(model) should be an integer identifying
!             the current quadratic model of the objective function.
!             if a previous step yielded a better function reduction,
!             then iv(model) will be set to iv(mlstgd) on output.
! iv(nfcall) (in)  invocation count for the objective function.
! iv(nfgcal) (i/o) value of iv(nfcall) at step that gave the biggest
!             function reduction this iteration.  iv(nfgcal) remains
!             unchanged until a function reduction is obtained.
! iv(radinc) (i/o) the number of radius increases (or minus the number
!             of decreases) so far this iteration.
! iv(restor) (out) set to 1 if v(f) has been restored and x should be
!             restored to its initial value, to 2 if x should be saved,
!             to 3 if x should be restored from the saved value, and to
!             0 otherwise.
!  iv(stage) (i/o) count of the number of models tried so far in the
!             current iteration.
! iv(stglim) (in)  maximum number of models to consider.
! iv(switch) (out) set to 0 unless a new model is being tried and it
!             gives a smaller function value than the previous model,
!             in which case assst sets iv(switch) = 1.
! iv(toobig) (in)  is nonzero if step was too big (e.g. if it caused
!             overflow).
!   iv(xirc) (i/o) value that iv(irc) would have in the absence of
!             convergence, false convergence, and oversized steps.
!
!   v values referenced
!
! v(afctol) (in)  absolute function convergence tolerance.  if the
!             absolute value of the current function value v(f) is less
!             than v(afctol), then assst returns with iv(irc) = 10.
! v(decfac) (in)  factor by which to decrease radius when iv(toobig) is
!             nonzero.
! v(dstnrm) (in)  the 2-norm of d*step.
! v(dstsav) (i/o) value of v(dstnrm) on saved step.
!   v(dst0) (in)  the 2-norm of d times the newton step (when defined,
!             i.e., for v(nreduc) >= 0).
!      v(f) (i/o) on both input and output, v(f) is the objective func-
!             tion value at x.  if x is restored to a previous value,
!             then v(f) is restored to the corresponding value.
!   v(fdif) (out) the function reduction v(f0) - v(f) (for the output
!             value of v(f) if an earlier step gave a bigger function
!             decrease, and for the input value of v(f) otherwise).
! v(flstgd) (i/o) saved value of v(f).
!     v(f0) (in)  objective function value at start of iteration.
! v(gtslst) (i/o) value of v(gtstep) on saved step.
! v(gtstep) (in)  inner product between step and gradient.
! v(incfac) (in)  minimum factor by which to increase radius.
!  v(lmaxs) (in)  maximum reasonable step size (and initial step bound).
!             if the actual function decrease is no more than twice
!             what was predicted, if a return with iv(irc) = 7, 8, 9,
!             or 10 does not occur, if v(dstnrm) > v(lmaxs), and if
!             v(preduc) <= v(sctol) * abs(v(f0)), then assst re-
!             turns with iv(irc) = 11.  if so doing appears worthwhile,
!             then assst repeats this test with v(preduc) computed for
!             a step of length v(lmaxs) (by a return with iv(irc) = 6).
! v(nreduc) (i/o)  function reduction predicted by quadratic model for
!             newton step.  if assst is called with iv(irc) = 6, i.e.,
!             if v(preduc) has been computed with radius = v(lmaxs) for
!             use in the singular convervence test, then v(nreduc) is
!             set to -v(preduc) before the latter is restored.
! v(plstgd) (i/o) value of v(preduc) on saved step.
! v(preduc) (i/o) function reduction predicted by quadratic model for
!             current step.
! v(radfac) (out) factor to be used in determining the new radius,
!             which should be v(radfac)*dst, where  dst  is either the
!             output value of v(dstnrm) or the 2-norm of
!             diag(newd)*step  for the output value of step and the
!             updated version, newd, of the scale vector d.  for
!             iv(irc) = 3, v(radfac) = 1.0D+00 is returned.
! v(rdfcmn) (in)  minimum value for v(radfac) in terms of the input
!             value of v(dstnrm) -- suggested value = 0.1.
! v(rdfcmx) (in)  maximum value for v(radfac) -- suggested value = 4.0.
!  v(reldx) (in) scaled relative change in x caused by step, computed
!             (e.g.) by function  reldst  as
!                 max (d(i)*abs(x(i)-x0(i)), 1 <= i <= p) /
!                    max (d(i)*(abs(x(i))+abs(x0(i))), 1 <= i <= p).
! v(rfctol) (in)  relative function convergence tolerance.  if the
!             actual function reduction is at most twice what was pre-
!             dicted and  v(nreduc) <= v(rfctol)*abs(v(f0)),  then
!             assst returns with iv(irc) = 8 or 9.
! v(stppar) (in)  marquardt parameter -- 0 means full newton step.
! v(tuner1) (in)  tuning constant used to decide if the function
!             reduction was much less than expected.  suggested
!             value = 0.1.
! v(tuner2) (in)  tuning constant used to decide if the function
!             reduction was large enough to accept step.  suggested
!             value = 10**-4.
! v(tuner3) (in)  tuning constant used to decide if the radius
!             should be increased.  suggested value = 0.75.
!  v(xctol) (in)  x-convergence criterion.  if step is a newton step
!             (v(stppar) = 0) having v(reldx) <= v(xctol) and giving
!             at most twice the predicted function decrease, then
!             assst returns iv(irc) = 7 or 9.
!  v(xftol) (in)  false convergence tolerance.  if step gave no or only
!             a small function decrease and v(reldx) <= v(xftol),
!             then assst returns with iv(irc) = 12.
!
!  notes
!
!   application and usage restrictions
!
!        this routine is called as part of the nl2sol (nonlinear
!     least-squares) package.  it may be used in any unconstrained
!     minimization solver that uses dogleg, goldfeld-quandt-trotter,
!     or levenberg-marquardt steps.
!
!   algorithm notes
!
!        see (1) for further discussion of the assessing and model
!     switching strategies.  while nl2sol considers only two models,
!     assst is designed to handle any number of models.
!
!   usage notes
!
!        on the first call of an iteration, only the i/o variables
!     step, x, iv(irc), iv(model), v(f), v(dstnrm), v(gtstep), and
!     v(preduc) need have been initialized.  between calls, no i/o
!     values execpt step, x, iv(model), v(f) and the stopping toler-
!     ances should be changed.
!        after a return for convergence or false convergence, one can
!     change the stopping tolerances and call assst again, in which
!     case the stopping tests will be repeated.
!
!   history
!
!        john dennis designed much of this routine, starting with
!     ideas in (2). roy welsch suggested the model switching strategy.
!        david gay and stephen peters cast this subroutine into a more
!     portable form (winter 1977), and david gay cast it into its
!     present form (fall 1978).
!
  integer liv
  integer lv

  integer iv(liv)
  double precision v(lv)
  logical goodx
  integer i, nfc
  double precision emax, emaxs, gts, rfac1, xmax
  double precision half, one, onep2, two
  integer afctol, decfac, dstnrm, dstsav, dst0, f, fdif, flstgd, f0
  integer gtslst, gtstep, incfac, irc, lmaxs, mlstgd, model, nfcall
  integer nfgcal, nreduc, plstgd, preduc, radfac, radinc, rdfcmn
  integer rdfcmx, reldx, restor, rfctol, sctol, stage, stglim
  integer stppar, switch, toobig, tuner1, tuner2, tuner3, xctol
  integer xftol, xirc

  parameter ( half=0.5d+0, one=1.d+0, onep2=1.2d+0, two=2.d+0)
  parameter ( irc=29, mlstgd=32, model=5, nfcall=6, nfgcal=7 )
  parameter ( radinc=8, restor=9, stage=10, stglim=11, switch=12 )
  parameter ( toobig=2, xirc=13)
  parameter (afctol=31, decfac=22, dstnrm=2, dst0=3, dstsav=18 )
  parameter (f=10, fdif=11, flstgd=12, f0=13, gtslst=14, gtstep=4 )
  parameter (incfac=23, lmaxs=36, nreduc=6, plstgd=15, preduc=7 )
  parameter (radfac=16, rdfcmn=24, rdfcmx=25, reldx=17, rfctol=32 )
  parameter (sctol=37, stppar=5, tuner1=26, tuner2=27, tuner3=28 )
  parameter (xctol=33, xftol=34)

  nfc = iv(nfcall)
  iv(switch) = 0
  iv(restor) = 0
  rfac1 = one
  goodx = .true.
  i = iv(irc)

  if (i >= 1 .and. i <= 12) then
        go to (20,30,10,10,40,280,220,220,220,220,220,170), i
  end if

  iv(irc) = 13
  return
!
!  Initialize for new iteration.
!
 10   iv(stage) = 1
  iv(radinc) = 0
  v(flstgd) = v(f0)
  if (iv(toobig) == 0) go to 110
     iv(stage) = -1
     iv(xirc) = i
     go to 60
!
!  Step was recomputed with new model or smaller radius
!  first decide which
!
 20   if (iv(model) /= iv(mlstgd)) go to 30
!
!  Old model retained, smaller radius tried
!  do not consider any more new models this iteration
!
     iv(stage) = iv(stglim)
     iv(radinc) = -1
     go to 110
!
!  A new model is being tried.  decide whether to keep it.
!
 30   iv(stage) = iv(stage) + 1
!
!  Now we add the possibiltiy that step was recomputed with
!  the same model, perhaps because of an oversized step.
!
 40   if (iv(stage) > 0) go to 50
!
!  Step was recomputed because it was too big.
!
     if (iv(toobig) /= 0) go to 60
!
!  Restore iv(stage) and pick up where we left off.
!
     iv(stage) = -iv(stage)
     i = iv(xirc)
     go to (20, 30, 110, 110, 70), i

 50   if (iv(toobig) == 0) go to 70
!
!  Handle oversize step
!
  if (iv(radinc) > 0) go to 80
     iv(stage) = -iv(stage)
     iv(xirc) = iv(irc)

 60      v(radfac) = v(decfac)
     iv(radinc) = iv(radinc) - 1
     iv(irc) = 5
     iv(restor) = 1
     return

 70   if (v(f) < v(flstgd)) go to 110
!
!  The new step is a loser.  restore old model.
!
  if (iv(model) == iv(mlstgd)) go to 80
     iv(model) = iv(mlstgd)
     iv(switch) = 1
!
!  Restore step, etc. only if a previous step decreased v(f).
!
 80   if (v(flstgd) >= v(f0)) go to 110
     iv(restor) = 1
     v(f) = v(flstgd)
     v(preduc) = v(plstgd)
     v(gtstep) = v(gtslst)
     if (iv(switch) == 0) rfac1 = v(dstnrm) / v(dstsav)
     v(dstnrm) = v(dstsav)
     nfc = iv(nfgcal)
     goodx = .false.

 110  v(fdif) = v(f0) - v(f)
  if (v(fdif) > v(tuner2) * v(preduc)) go to 140
  if(iv(radinc)>0) go to 140
!
!         no (or only a trivial) function decrease
!         so try new model or smaller radius
!
     if (v(f) < v(f0)) go to 120
          iv(mlstgd) = iv(model)
          v(flstgd) = v(f)
          v(f) = v(f0)
          iv(restor) = 1
          go to 130
 120     iv(nfgcal) = nfc
 130     iv(irc) = 1
     if (iv(stage) < iv(stglim)) go to 160
          iv(irc) = 5
          iv(radinc) = iv(radinc) - 1
          go to 160
!
!  Nontrivial function decrease achieved
!
 140  iv(nfgcal) = nfc
  rfac1 = 1.0D+00
  v(dstsav) = v(dstnrm)
  if (v(fdif) > v(preduc)*v(tuner1)) go to 190
!
!  Decrease was much less than predicted -- either change models
!  or accept step with decreased radius.
!
  if (iv(stage) >= iv(stglim)) go to 150
!
!  Consider switching models
!
     iv(irc) = 2
     go to 160
!
!  Accept step with decreased radius
!
 150  iv(irc) = 4
!
!   set v(radfac) to fletcher*s decrease factor
!
 160  iv(xirc) = iv(irc)
  emax = v(gtstep) + v(fdif)
  v(radfac) = half * rfac1

  if (emax < v(gtstep)) then
    v(radfac) = rfac1 * max (v(rdfcmn),half * v(gtstep)/emax)
  end if
!
!  Do false convergence test
!
 170  if (v(reldx) <= v(xftol)) go to 180
     iv(irc) = iv(xirc)
     if (v(f) < v(f0)) go to 200
          go to 230

 180  iv(irc) = 12
  go to 240
!
!  Handle good function decrease
!
 190  if (v(fdif) < (-v(tuner3) * v(gtstep))) go to 210
!
!  Increasing radius looks worthwhile.  see if we just
!  recomputed step with a decreased radius or restored step
!  after recomputing it with a larger radius.
!
  if (iv(radinc) < 0) go to 210
  if (iv(restor) == 1) go to 210
!
!  We did not.  try a longer step unless this was a newton step.
!
     v(radfac) = v(rdfcmx)
     gts = v(gtstep)
     if (v(fdif) < (half/v(radfac) - 1.0D+00 ) * gts) then
       v(radfac) = max (v(incfac), half*gts/(gts + v(fdif)))
     end if
     iv(irc) = 4
     if (v(stppar) == 0.0D+00 ) go to 230
     if (v(dst0) >= 0.0D+00 .and. (v(dst0) < two*v(dstnrm) &
              .or. v(nreduc) < onep2*v(fdif)))  then
       go to 230
     end if
!
!  Step was not a newton step.  recompute it with a larger radius.
!
          iv(irc) = 5
          iv(radinc) = iv(radinc) + 1
!
!  Save values corresponding to good step
!
 200  v(flstgd) = v(f)
  iv(mlstgd) = iv(model)
  if (iv(restor) /= 1) iv(restor) = 2
  v(dstsav) = v(dstnrm)
  iv(nfgcal) = nfc
  v(plstgd) = v(preduc)
  v(gtslst) = v(gtstep)
  go to 230
!
!  Accept step with radius unchanged.
!
 210  v(radfac) = 1.0D+00
  iv(irc) = 3
  go to 230
!
!  Come here for a restart after convergence.
!
 220  iv(irc) = iv(xirc)
  if (v(dstsav) >= 0.0D+00 ) go to 240
     iv(irc) = 12
     go to 240
!
!  Perform convergence tests.
!
 230  iv(xirc) = iv(irc)
 240  if (iv(restor) == 1 .and. v(flstgd) < v(f0)) iv(restor) = 3
  if (abs(v(f)) < v(afctol)) iv(irc) = 10

  if (half * v(fdif) > v(preduc)) then
    return
  end if

  emax = v(rfctol) * abs(v(f0))
  emaxs = v(sctol) * abs(v(f0))
  if (v(dstnrm) > v(lmaxs) .and. v(preduc) <= emaxs) then
    iv(irc) = 11
  end if
  if (v(dst0) < 0.0D+00 ) go to 250
  i = 0

  if ((v(nreduc) > 0.0D+00 .and. v(nreduc) <= emax) .or. &
      (v(nreduc) == 0.0D+00 .and. v(preduc) == 0.0D+00 )) then
    i = 2
  end if

  if (v(stppar) == 0.0D+00 .and. v(reldx) <= v(xctol) .and. goodx) then
    i = i + 1
  end if

  if (i > 0) iv(irc) = i + 6
!
!  Consider recomputing step of length v(lmaxs) for singular
!  convergence test.
!
 250  if (iv(irc) > 5 .and. iv(irc) /= 12) then
     return
  end if

  if (v(dstnrm) > v(lmaxs)) go to 260
     if (v(preduc) >= emaxs) then
       return
     end if
          if (v(dst0) <= 0.0D+00 ) go to 270
               if (half * v(dst0) <= v(lmaxs)) then
                 return
               end if
                    go to 270
 260  if (half * v(dstnrm) <= v(lmaxs)) then
        return
      end if
  xmax = v(lmaxs) / v(dstnrm)
  if (xmax * (two - xmax) * v(preduc) >= emaxs) then
    return
  end if
 270  if (v(nreduc) < 0.0D+00 ) go to 290
!
!   recompute v(preduc) for use in singular convergence test
!
  v(gtslst) = v(gtstep)
  v(dstsav) = v(dstnrm)
  if (iv(irc) == 12) v(dstsav) = -v(dstsav)
  v(plstgd) = v(preduc)
  i = iv(restor)
  iv(restor) = 2
  if (i == 3) iv(restor) = 0
  iv(irc) = 6
  return
!
!  Perform singular convergence test with recomputed v(preduc)
!
 280  v(gtstep) = v(gtslst)
  v(dstnrm) = abs(v(dstsav))
  iv(irc) = iv(xirc)
  if (v(dstsav) <= 0.0D+00 ) iv(irc) = 12
  v(nreduc) = -v(preduc)
  v(preduc) = v(plstgd)
  iv(restor) = 3

 290  if (-v(nreduc) <= v(rfctol) * abs(v(f0))) iv(irc) = 11

  return
end

subroutine wzbfgs ( l, n, s, w, y, z )
!
!! WZBFGS compute Y and Z for LUPDAT corresponding to BFGS update.
!
!  Discussion:
!
!    When S is computed in certain ways, for example by GQTSTP or
!    DBLDOG, it is possible to save N**2/2 operations since L'*S
!    or L*L'*S is then known.
!
!    If the BFGS update to L*L' would reduce its determinant to
!    less than EPS times its old value, then this routine in effect
!    replaces Y by THETA*Y + (1-THETA)*L*L'*S, where THETA
!    (between 0 and 1) is chosen to make the reduction factor = EPS.
!
!  Parameters:
!
!    l (i/o) cholesky factor of hessian, a lower triang. matrix stored
!             compactly by rows.
!
!    n (input) order of  l  and length of  s,  w,  y,  z.
!
!    s (input) the step just taken.
!
!    w (output) right singular vector of rank 1 correction to l.
!
!    y (input) change in gradients corresponding to s.
!
!    z (output) left singular vector of rank 1 correction to l.
!
  implicit none

  integer n

  double precision dotprd
  double precision cs
  double precision cy
  double precision, parameter :: eps = 0.1D+00
  double precision epsrt
  double precision l(n*(n+1)/2)
  double precision s(n)
  double precision shs
  double precision theta
  double precision w(n)
  double precision y(n)
  double precision ys
  double precision z(n)

  call ltvmul ( n, w, l, s )
  shs = dotprd ( n, w, w )
  ys = dotprd ( n, y, s )

  if ( ys < eps * shs ) then
    theta = ( 1.0D+00 - eps ) * shs / ( shs - ys )
    epsrt = sqrt ( eps )
    cy = theta / ( shs * epsrt )
    cs = ( 1.0D+00 + ( theta - 1.0D+00 ) / epsrt ) / shs
  else
    cy = 1.0D+00 / ( sqrt ( ys ) * sqrt ( shs ) )
    cs = 1.0D+00 / shs
  end if

  call livmul ( n, z, l, y )

  z(1:n) = cy * z(1:n) - cs * w(1:n)

  return
end

!======================================================================!
!                                                                      !
!             E L E C T R O S T A T I C     M O M E N T S              !
!                                                                      !
!        Module written by Ion Mitxelena (ion.mitxelena@ehu.eus)       !
!                                                                      !
! ==================================================================== !
!                                                                      !
!           Subroutines for electrostatic moment calculations          !
!           including dipole, quadrupole and octupole moments          !
!                 ( J. Chem. Phys. 144, 204108, 2016 )                 !
!                                                                      !
!   DIPMOMr:  Electronic, nuclear and total dipole moments             !
!   QUADMOMr: Electronic, nuclear and total quadrupole moments         !
!   OCTMOMr:  Electronic, nuclear and total octupole moments           !
!                                                                      !
!   PASSDIPUSER:  Pass atomic dipole matrices to USER(13,14,15,16)     !
!   PASSUSERDIP:  Pass USER(13,14,15,16) to atomic dipole matrices     !
!   PASSQUADUSER: Pass atomic quadrupole matrices to USER(18-24)       !
!   PASSUSERQUAD: Pass USER(18-24) to atomic quadrupole matrices       !
!   PASSOCTUSER:  Pass atomic octupole matrices to USER(31-41)         !
!   PASSUSEROCT:  Pass USER(31-41) to atomic octupole matrices         !
!                                                                      !
!======================================================================!
                                                                       
!----------------------------------------------------------------------!
!          Molecular dipole, quadrupole and octupole matrices          !
!----------------------------------------------------------------------!

! DIPMOMr
      SUBROUTINE DIPMOMr(DIPN,ADIPx,ADIPy,ADIPz,DIPx,DIPy,DIPz,QD,RO, &
                         DM1e,DM2e,DM3e,DM1,DM2,DM3,DTOTAL)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::DIPN
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO,DIPx,DIPy,DIPz
      DOUBLE PRECISION,DIMENSION(NBF,NBF),INTENT(IN)::ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF),INTENT(IN)::QD
      DOUBLE PRECISION,INTENT(OUT)::DM1,DM2,DM3,DTOTAL
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::AUX
!-----------------------------------------------------------------------
!     Electron Contribution to Dipole Moment
!-----------------------------------------------------------------------
      DM1e = 0.0
      DM2e = 0.0
      DM3e = 0.0
!
      ALLOCATE (AUX(NBF,NBF))
      DO J=1,NB
       AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
       CALL TRACEm(DIPx(J),AUX,ADIPx,NBF)
       DM1e = DM1e - 2.0*RO(J)*DIPx(J)
       CALL TRACEm(DIPy(J),AUX,ADIPy,NBF)
       DM2e = DM2e - 2.0*RO(J)*DIPy(J)
       CALL TRACEm(DIPz(J),AUX,ADIPz,NBF)
       DM3e = DM3e - 2.0*RO(J)*DIPz(J)
      ENDDO
!      
      IF(NSOC>0)THEN
       if(.not.HighSpin)then
        DO J=NB+1,NA
         AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
         CALL TRACEm(DIPx(J),AUX,ADIPx,NBF)
         DM1e = DM1e - 2.0*RO(J)*DIPx(J)
         CALL TRACEm(DIPy(J),AUX,ADIPy,NBF)
         DM2e = DM2e - 2.0*RO(J)*DIPy(J)
         CALL TRACEm(DIPz(J),AUX,ADIPz,NBF)
         DM3e = DM3e - 2.0*RO(J)*DIPz(J)
        ENDDO
       else if(HighSpin)then
        DO J=NB+1,NA
         AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
         CALL TRACEm(DIPx(J),AUX,ADIPx,NBF)
         DM1e = DM1e - RO(J)*DIPx(J)
         CALL TRACEm(DIPy(J),AUX,ADIPy,NBF)
         DM2e = DM2e - RO(J)*DIPy(J)
         CALL TRACEm(DIPz(J),AUX,ADIPz,NBF)
         DM3e = DM3e - RO(J)*DIPz(J)
        ENDDO
       end if      
      ENDIF
!       
      DO J=NA+1,NBF5
       AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
       CALL TRACEm(DIPx(J),AUX,ADIPx,NBF)
       DM1e = DM1e - 2.0*RO(J)*DIPx(J)
       CALL TRACEm(DIPy(J),AUX,ADIPy,NBF)
       DM2e = DM2e - 2.0*RO(J)*DIPy(J)
       CALL TRACEm(DIPz(J),AUX,ADIPz,NBF)
       DM3e = DM3e - 2.0*RO(J)*DIPz(J)
      ENDDO
      DEALLOCATE (AUX)
!-----------------------------------------------------------------------
!     Nuclear Contribution to Dipole Moment
!-----------------------------------------------------------------------
      DM1 = DIPN(1) + DM1e
      DM2 = DIPN(2) + DM2e
      DM3 = DIPN(3) + DM3e
!-----------------------------------------------------------------------
!     Total Dipole Moment
!-----------------------------------------------------------------------
      DTOTAL = SQRT(DM1*DM1+DM2*DM2+DM3*DM3)
      RETURN
!-----------------------------------------------------------------------
      END SUBROUTINE DIPMOMr

! QUADMOMr
      SUBROUTINE QUADMOMr(QUADN,AQUADxx,AQUADyy,AQUADzz,AQUADxy, &
                          AQUADxz,AQUADyz,QUADxx,QUADyy,QUADzz,QUADxy, &
                          QUADxz,QUADyz,QD,RO,QM1e,QM2e,QM3e,QM4e,QM5e, &
                          QM6e,QM1,QM2,QM3,QM4,QM5,QM6)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION, DIMENSION(NBF,NBF,NBF), INTENT(IN) :: QD
      DOUBLE PRECISION, DIMENSION(NBF5), INTENT(IN) :: RO(NBF5)
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN) :: QUADN
      DOUBLE PRECISION, DIMENSION(NBF,NBF), INTENT(IN) :: AQUADxx,AQUADyy,AQUADzz
      DOUBLE PRECISION, DIMENSION(NBF,NBF), INTENT(IN) :: AQUADxy,AQUADxz,AQUADyz
      DOUBLE PRECISION, DIMENSION(NBF5) :: QUADxx,QUADyy,QUADzz,QUADxy,QUADxz,QUADyz
      DOUBLE PRECISION, INTENT(OUT) :: QM1,QM2,QM3,QM4,QM5,QM6
      DOUBLE PRECISION :: QM1e,QM2e,QM3e,QM4e,QM5e,QM6e
      INTEGER :: J
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::AUX
!-----------------------------------------------------------------------
!     Electron Contribution to Quadrupole Moment
!-----------------------------------------------------------------------
      QM1e=0.0d0
      QM2e=0.0d0
      QM3e=0.0d0
      QM4e=0.0d0
      QM5e=0.0d0
      QM6e=0.0d0
!
      ALLOCATE (AUX(NBF,NBF))
      DO J=1,NB
       AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
       CALL TRACEm(QUADxx(J),AUX,AQUADxx,NBF)
       QM1e = QM1e - 2.d0*RO(J)*QUADxx(J)
       CALL TRACEm(QUADyy(J),AUX,AQUADyy,NBF)
       QM2e = QM2e - 2.d0*RO(J)*QUADyy(J)
       CALL TRACEm(QUADzz(J),AUX,AQUADzz,NBF)
       QM3e = QM3e - 2.d0*RO(J)*QUADzz(J)
       CALL TRACEm(QUADxy(J),AUX,AQUADxy,NBF)
       QM4e = QM4e - 2.d0*RO(J)*QUADxy(J)
       CALL TRACEm(QUADxz(J),AUX,AQUADxz,NBF)
       QM5e = QM5e - 2.d0*RO(J)*QUADxz(J)
       CALL TRACEm(QUADyz(J),AUX,AQUADyz,NBF)
       QM6e = QM6e - 2.d0*RO(J)*QUADyz(J)
      END DO
!      
      IF(NSOC>0)THEN
       if(.not.HighSpin)then
        DO J=NB+1,NA
         AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
         CALL TRACEm(QUADxx(J),AUX,AQUADxx,NBF)
         QM1e = QM1e - 2.d0*RO(J)*QUADxx(J)
         CALL TRACEm(QUADyy(J),AUX,AQUADyy,NBF)
         QM2e = QM2e - 2.d0*RO(J)*QUADyy(J)
         CALL TRACEm(QUADzz(J),AUX,AQUADzz,NBF)
         QM3e = QM3e - 2.d0*RO(J)*QUADzz(J)
         CALL TRACEm(QUADxy(J),AUX,AQUADxy,NBF)
         QM4e = QM4e - 2.d0*RO(J)*QUADxy(J)
         CALL TRACEm(QUADxz(J),AUX,AQUADxz,NBF)
         QM5e = QM5e - 2.d0*RO(J)*QUADxz(J)
         CALL TRACEm(QUADyz(J),AUX,AQUADyz,NBF)
         QM6e = QM6e - 2.d0*RO(J)*QUADyz(J)
        END DO
       else if(HighSpin)then        
        DO J=NB+1,NA
         AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
         CALL TRACEm(QUADxx(J),AUX,AQUADxx,NBF)
         QM1e = QM1e - RO(J)*QUADxx(J)
         CALL TRACEm(QUADyy(J),AUX,AQUADyy,NBF)
         QM2e = QM2e - RO(J)*QUADyy(J)
         CALL TRACEm(QUADzz(J),AUX,AQUADzz,NBF)
         QM3e = QM3e - RO(J)*QUADzz(J)
         CALL TRACEm(QUADxy(J),AUX,AQUADxy,NBF)
         QM4e = QM4e - RO(J)*QUADxy(J)
         CALL TRACEm(QUADxz(J),AUX,AQUADxz,NBF)
         QM5e = QM5e - RO(J)*QUADxz(J)
         CALL TRACEm(QUADyz(J),AUX,AQUADyz,NBF)
         QM6e = QM6e - RO(J)*QUADyz(J)
        END DO
       end if
      END IF
!
      DO J=NA+1,NBF5
       AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
       CALL TRACEm(QUADxx(J),AUX,AQUADxx,NBF)
       QM1e = QM1e - 2.d0*RO(J)*QUADxx(J)
       CALL TRACEm(QUADyy(J),AUX,AQUADyy,NBF)
       QM2e = QM2e - 2.d0*RO(J)*QUADyy(J)
       CALL TRACEm(QUADzz(J),AUX,AQUADzz,NBF)
       QM3e = QM3e - 2.d0*RO(J)*QUADzz(J)
       CALL TRACEm(QUADxy(J),AUX,AQUADxy,NBF)
       QM4e = QM4e - 2.d0*RO(J)*QUADxy(J)
       CALL TRACEm(QUADxz(J),AUX,AQUADxz,NBF)
       QM5e = QM5e - 2.d0*RO(J)*QUADxz(J)
       CALL TRACEm(QUADyz(J),AUX,AQUADyz,NBF)
       QM6e = QM6e - 2.d0*RO(J)*QUADyz(J)
      END DO
      DEALLOCATE (AUX)
!-----------------------------------------------------------------------
!     Add Nuclear Contribution to Quadrupole Moment
!-----------------------------------------------------------------------
      QM1e = QUADN(1) + QM1e
      QM2e = QUADN(2) + QM2e
      QM3e = QUADN(3) + QM3e
      QM4e = QUADN(4) + QM4e
      QM5e = QUADN(5) + QM5e
      QM6e = QUADN(6) + QM6e
!-----------------------------------------------------------------------
!     Form Quadrupole Tensor (BUCKINHAM)
!-----------------------------------------------------------------------
      QM1 = 0.5D0*(QM1e + QM1e - QM2e - QM3e)
      QM2 = 0.5D0*(QM2e + QM2e - QM1e - QM3e)
      QM3 = 0.5D0*(QM3e + QM3e - QM1e - QM2e)
      QM4 = 1.5D0*QM4e
      QM5 = 1.5D0*QM5e
      QM6 = 1.5D0*QM6e
!-----------------------------------------------------------------------    
      RETURN
      END SUBROUTINE QUADMOMr

! OCTMOMr
      SUBROUTINE OCTMOMr(OCTUN,AOCTxxx,AOCTyyy,AOCTzzz, &
                         AOCTxxy,AOCTxxz,AOCTxyy,AOCTyyz,AOCTxzz, &
                         AOCTyzz,AOCTxyz,OCTxxx,OCTyyy,OCTzzz, &
                         OCTxxy,OCTxxz,OCTxyy,OCTyyz,OCTxzz, &
                         OCTyzz,OCTxyz,   &
                         QD,RO,OMXXXe,OMYYYe,OMZZZe,  &
                         OMXXYe,OMXXZe,OMXYYe,OMYYZe,OMXZZe, &
                         OMYZZe,OMXYZe,OXXX,OYYY,OZZZ, &
                         OXXY,OXXZ,OXYY,OYYZ,OXZZ, &
                         OYZZ,OXYZ)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION, DIMENSION(NBF,NBF,NBF), INTENT(IN) :: QD
      DOUBLE PRECISION, DIMENSION(NBF5), INTENT(IN) :: RO(NBF5)
      DOUBLE PRECISION, DIMENSION(10), INTENT(IN) :: OCTUN(10)
      DOUBLE PRECISION, DIMENSION(NBF,NBF), INTENT(IN) :: AOCTxxx,AOCTyyy,AOCTzzz
      DOUBLE PRECISION, DIMENSION(NBF,NBF), INTENT(IN) :: AOCTxyy,AOCTyyz,AOCTxzz
      DOUBLE PRECISION, DIMENSION(NBF,NBF), INTENT(IN) :: AOCTxxy,AOCTxxz
      DOUBLE PRECISION, DIMENSION(NBF,NBF), INTENT(IN) :: AOCTyzz,AOCTxyz
      DOUBLE PRECISION, DIMENSION(NBF5) :: OCTxxx,OCTyyy,OCTzzz,OCTxxy,OCTxxz
      DOUBLE PRECISION, DIMENSION(NBF5) :: OCTxyy,OCTyyz,OCTxzz,OCTyzz,OCTxyz
      DOUBLE PRECISION :: OMXXXe,OMYYYe,OMZZZe,OMXXYe,OMXXZe
      DOUBLE PRECISION :: OMXYYe,OMYYZe,OMXZZe,OMYZZe,OMXYZe
      DOUBLE PRECISION, INTENT(OUT) :: OXXX,OYYY,OZZZ,OXXY,OXXZ
      DOUBLE PRECISION, INTENT(OUT) :: OXYY,OYYZ,OXZZ,OYZZ,OXYZ
      INTEGER :: J
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::AUX
!-----------------------------------------------------------------------
!     Electron Contribution to Quadrupole Moment
!-----------------------------------------------------------------------
      OMXXXe=0.0d0
      OMYYYe=0.0d0
      OMZZZe=0.0d0
      OMXXYe=0.0d0
      OMXXZe=0.0d0
      OMXYYe=0.0d0
      OMYYZe=0.0d0
      OMXZZe=0.0d0
      OMYZZe=0.0d0
      OMXYZe=0.0d0
!
      ALLOCATE (AUX(NBF,NBF))
      DO J=1,NB
        AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
        CALL TRACEm(OCTxxx(J),AUX,AOCTxxx,NBF)
        OMXXXe = OMXXXe - 2.d0*RO(J)*OCTxxx(J)
        CALL TRACEm(OCTyyy(J),AUX,AOCTyyy,NBF)
        OMYYYe = OMYYYe - 2.d0*RO(J)*OCTyyy(J)
        CALL TRACEm(OCTzzz(J),AUX,AOCTzzz,NBF)
        OMZZZe = OMZZZe - 2.d0*RO(J)*OCTzzz(J)
        CALL TRACEm(OCTxxy(J),AUX,AOCTxxy,NBF)
        OMXXYe = OMXXYe - 2.d0*RO(J)*OCTxxy(J)
        CALL TRACEm(OCTxxz(J),AUX,AOCTxxz,NBF)
        OMXXZe = OMXXZe - 2.d0*RO(J)*OCTxxz(J)
        CALL TRACEm(OCTxyy(J),AUX,AOCTxyy,NBF)
        OMXYYe = OMXYYe - 2.d0*RO(J)*OCTxyy(J)
        CALL TRACEm(OCTyyz(J),AUX,AOCTyyz,NBF)
        OMYYZe = OMYYZe - 2.d0*RO(J)*OCTyyz(J)
        CALL TRACEm(OCTxzz(J),AUX,AOCTxzz,NBF)
        OMXZZe = OMXZZe - 2.d0*RO(J)*OCTxzz(J)
        CALL TRACEm(OCTyzz(J),AUX,AOCTyzz,NBF)
        OMYZZe = OMYZZe - 2.d0*RO(J)*OCTyzz(J)
        CALL TRACEm(OCTxyz(J),AUX,AOCTxyz,NBF)
        OMXYZe = OMXYZe - 2.d0*RO(J)*OCTxyz(J)
      END DO
!      
      IF(NSOC>0)THEN
       if(.not.HighSpin)then
        DO J=NB+1,NA
          AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
          CALL TRACEm(OCTxxx(J),AUX,AOCTxxx,NBF)
          OMXXXe = OMXXXe - 2.d0*RO(J)*OCTxxx(J)
          CALL TRACEm(OCTyyy(J),AUX,AOCTyyy,NBF)
          OMYYYe = OMYYYe - 2.d0*RO(J)*OCTyyy(J)
          CALL TRACEm(OCTzzz(J),AUX,AOCTzzz,NBF)
          OMZZZe = OMZZZe - 2.d0*RO(J)*OCTzzz(J)
          CALL TRACEm(OCTxxy(J),AUX,AOCTxxy,NBF)
          OMXXYe = OMXXYe - 2.d0*RO(J)*OCTxxy(J)
          CALL TRACEm(OCTxxz(J),AUX,AOCTxxz,NBF)
          OMXXZe = OMXXZe - 2.d0*RO(J)*OCTxxz(J)
          CALL TRACEm(OCTxyy(J),AUX,AOCTxyy,NBF)
          OMXYYe = OMXYYe - 2.d0*RO(J)*OCTxyy(J)
          CALL TRACEm(OCTyyz(J),AUX,AOCTyyz,NBF)
          OMYYZe = OMYYZe - 2.d0*RO(J)*OCTyyz(J)
          CALL TRACEm(OCTxzz(J),AUX,AOCTxzz,NBF)
          OMXZZe = OMXZZe - 2.d0*RO(J)*OCTxzz(J)
          CALL TRACEm(OCTyzz(J),AUX,AOCTyzz,NBF)
          OMYZZe = OMYZZe - 2.d0*RO(J)*OCTyzz(J)
          CALL TRACEm(OCTxyz(J),AUX,AOCTxyz,NBF)
          OMXYZe = OMXYZe - 2.d0*RO(J)*OCTxyz(J)
        END DO
       else if(HighSpin)then
        DO J=NB+1,NA
          AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
          CALL TRACEm(OCTxxx(J),AUX,AOCTxxx,NBF)
          OMXXXe = OMXXXe - RO(J)*OCTxxx(J)
          CALL TRACEm(OCTyyy(J),AUX,AOCTyyy,NBF)
          OMYYYe = OMYYYe - RO(J)*OCTyyy(J)
          CALL TRACEm(OCTzzz(J),AUX,AOCTzzz,NBF)
          OMZZZe = OMZZZe - RO(J)*OCTzzz(J)
          CALL TRACEm(OCTxxy(J),AUX,AOCTxxy,NBF)
          OMXXYe = OMXXYe - RO(J)*OCTxxy(J)
          CALL TRACEm(OCTxxz(J),AUX,AOCTxxz,NBF)
          OMXXZe = OMXXZe - RO(J)*OCTxxz(J)
          CALL TRACEm(OCTxyy(J),AUX,AOCTxyy,NBF)
          OMXYYe = OMXYYe - RO(J)*OCTxyy(J)
          CALL TRACEm(OCTyyz(J),AUX,AOCTyyz,NBF)
          OMYYZe = OMYYZe - RO(J)*OCTyyz(J)
          CALL TRACEm(OCTxzz(J),AUX,AOCTxzz,NBF)
          OMXZZe = OMXZZe - RO(J)*OCTxzz(J)
          CALL TRACEm(OCTyzz(J),AUX,AOCTyzz,NBF)
          OMYZZe = OMYZZe - RO(J)*OCTyzz(J)
          CALL TRACEm(OCTxyz(J),AUX,AOCTxyz,NBF)
          OMXYZe = OMXYZe - RO(J)*OCTxyz(J)
        END DO
       end if
      END IF
!       
      DO J=1,NBF5
        AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
        CALL TRACEm(OCTxxx(J),AUX,AOCTxxx,NBF)
        OMXXXe = OMXXXe - 2.d0*RO(J)*OCTxxx(J)
        CALL TRACEm(OCTyyy(J),AUX,AOCTyyy,NBF)
        OMYYYe = OMYYYe - 2.d0*RO(J)*OCTyyy(J)
        CALL TRACEm(OCTzzz(J),AUX,AOCTzzz,NBF)
        OMZZZe = OMZZZe - 2.d0*RO(J)*OCTzzz(J)
        CALL TRACEm(OCTxxy(J),AUX,AOCTxxy,NBF)
        OMXXYe = OMXXYe - 2.d0*RO(J)*OCTxxy(J)
        CALL TRACEm(OCTxxz(J),AUX,AOCTxxz,NBF)
        OMXXZe = OMXXZe - 2.d0*RO(J)*OCTxxz(J)
        CALL TRACEm(OCTxyy(J),AUX,AOCTxyy,NBF)
        OMXYYe = OMXYYe - 2.d0*RO(J)*OCTxyy(J)
        CALL TRACEm(OCTyyz(J),AUX,AOCTyyz,NBF)
        OMYYZe = OMYYZe - 2.d0*RO(J)*OCTyyz(J)
        CALL TRACEm(OCTxzz(J),AUX,AOCTxzz,NBF)
        OMXZZe = OMXZZe - 2.d0*RO(J)*OCTxzz(J)
        CALL TRACEm(OCTyzz(J),AUX,AOCTyzz,NBF)
        OMYZZe = OMYZZe - 2.d0*RO(J)*OCTyzz(J)
        CALL TRACEm(OCTxyz(J),AUX,AOCTxyz,NBF)
        OMXYZe = OMXYZe - 2.d0*RO(J)*OCTxyz(J)
      END DO
      DEALLOCATE (AUX)
!-----------------------------------------------------------------------
!     Nuclear Contribution to Quadrupole Moment
!-----------------------------------------------------------------------
      OMXXXe = OCTUN(1) + OMXXXe
      OMYYYe = OCTUN(2) + OMYYYe
      OMZZZe = OCTUN(3) + OMZZZe
      OMXXYe = OCTUN(4) + OMXXYe
      OMXXZe = OCTUN(5) + OMXXZe
      OMXYYe = OCTUN(6) + OMXYYe
      OMYYZe = OCTUN(7) + OMYYZe
      OMXZZe = OCTUN(8) + OMXZZe
      OMYZZe = OCTUN(9) + OMYZZe
      OMXYZe = OCTUN(10) + OMXYZe
!-----------------------------------------------------------------------
!     Form Octupole Tensor
!-----------------------------------------------------------------------
      OXXX = OMXXXe - 1.5D+00*OMXYYe - 1.5D+00*OMXZZe
      OYYY = OMYYYe - 1.5D+00*OMXXYe - 1.5D+00*OMYZZe
      OZZZ = OMZZZe - 1.5D+00*OMXXZe - 1.5D+00*OMYYZe
      OXXY = 2.0D+00*OMXXYe - 0.5D+00*OMYYYe - 0.5D+00*OMYZZe
      OXXZ = 2.0D+00*OMXXZe - 0.5D+00*OMYYZe - 0.5D+00*OMZZZe
      OXYY = 2.0D+00*OMXYYe - 0.5D+00*OMXXXe - 0.5D+00*OMXZZe
      OYYZ = 2.0D+00*OMYYZe - 0.5D+00*OMXXZe - 0.5D+00*OMZZZe
      OXZZ = 2.0D+00*OMXZZe - 0.5D+00*OMXXXe - 0.5D+00*OMXYYe
      OYZZ = 2.0D+00*OMYZZe - 0.5D+00*OMXXYe - 0.5D+00*OMYYYe
      OXYZ = 2.5D+00*OMXYZe
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE OCTMOMr

!----------------------------------------------------------------------!
!           Atomic dipole, quadrupole and octupole matrices            !
!----------------------------------------------------------------------!

! PASSDIPUSER
      SUBROUTINE PASSDIPUSER(DIPN,ADIPx,ADIPy,ADIPz,USER)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(3)::DIPN
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
!-----------------------------------------------------------------------
      CALL XtoX0(DIPN ,USER(N11),  3)
      CALL XtoX0(ADIPx,USER(N12),NSQ)
      CALL XtoX0(ADIPy,USER(N13),NSQ)
      CALL XtoX0(ADIPz,USER(N14),NSQ)
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE PASSDIPUSER

! PASSUSERDIP
      SUBROUTINE PASSUSERDIP(DIPN,ADIPx,ADIPy,ADIPz,USER)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(3)::DIPN
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
!-----------------------------------------------------------------------
      CALL XtoX0(USER(N11),DIPN ,  3)
      CALL XtoX0(USER(N12),ADIPx,NSQ)
      CALL XtoX0(USER(N13),ADIPy,NSQ)
      CALL XtoX0(USER(N14),ADIPz,NSQ)
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE PASSUSERDIP

! PASSQUADUSER
      SUBROUTINE PASSQUADUSER(QUADN,AQUADxx,AQUADyy,AQUADzz,AQUADxy, &
                              AQUADxz,AQUADyz,USER)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(6)::QUADN
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AQUADxx,AQUADyy,AQUADzz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AQUADxy,AQUADxz,AQUADyz
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
!-----------------------------------------------------------------------
      CALL XtoX0(QUADN  ,USER(N18),  6)
      CALL XtoX0(AQUADxx,USER(N19),NSQ)
      CALL XtoX0(AQUADyy,USER(N20),NSQ)
      CALL XtoX0(AQUADzz,USER(N21),NSQ)
      CALL XtoX0(AQUADxy,USER(N22),NSQ)
      CALL XtoX0(AQUADxz,USER(N23),NSQ)
      CALL XtoX0(AQUADyz,USER(N24),NSQ)
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE PASSQUADUSER

! PASSUSERQUAD
      SUBROUTINE PASSUSERQUAD(QUADN,AQUADxx,AQUADyy,AQUADzz,         &
                                    AQUADxy,AQUADxz,AQUADyz,USER)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(6)::QUADN
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AQUADxx,AQUADyy,AQUADzz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AQUADxy,AQUADxz,AQUADyz
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
!-----------------------------------------------------------------------
      CALL XtoX0(USER(N18),QUADN  ,  6)
      CALL XtoX0(USER(N19),AQUADxx,NSQ)
      CALL XtoX0(USER(N20),AQUADyy,NSQ)
      CALL XtoX0(USER(N21),AQUADzz,NSQ)
      CALL XtoX0(USER(N22),AQUADxy,NSQ)
      CALL XtoX0(USER(N23),AQUADxz,NSQ)
      CALL XtoX0(USER(N24),AQUADyz,NSQ)
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE PASSUSERQUAD

! PASSOCTUSER
      SUBROUTINE PASSOCTUSER(OCTUN,AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxy,   &
                             AOCTxxz,AOCTxyy,AOCTyyz,AOCTxzz,AOCTyzz, &
                             AOCTxyz,USER)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(10)::OCTUN
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxy,AOCTxxz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTxyy,AOCTyyz,AOCTxzz,AOCTyzz,AOCTxyz
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
!-----------------------------------------------------------------------
      CALL XtoX0(OCTUN  ,USER(N31), 10)
      CALL XtoX0(AOCTxxx,USER(N32),NSQ)
      CALL XtoX0(AOCTyyy,USER(N33),NSQ)
      CALL XtoX0(AOCTzzz,USER(N34),NSQ)
      CALL XtoX0(AOCTxxy,USER(N35),NSQ)
      CALL XtoX0(AOCTxxz,USER(N36),NSQ)
      CALL XtoX0(AOCTxyy,USER(N37),NSQ)
      CALL XtoX0(AOCTyyz,USER(N38),NSQ)
      CALL XtoX0(AOCTxzz,USER(N39),NSQ)
      CALL XtoX0(AOCTyzz,USER(N40),NSQ)
      CALL XtoX0(AOCTxyz,USER(N41),NSQ)
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE PASSOCTUSER

! PASSUSEROCT
      SUBROUTINE PASSUSEROCT(OCTUN,AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxy,   &
                             AOCTxxz,AOCTxyy,AOCTyyz,AOCTxzz,AOCTyzz, &
                             AOCTxyz,USER)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(10)::OCTUN
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxy,AOCTxxz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTxyy,AOCTyyz,AOCTxzz,AOCTyzz,AOCTxyz
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
!-----------------------------------------------------------------------
      CALL XtoX0(USER(N31),OCTUN  , 10)
      CALL XtoX0(USER(N32),AOCTxxx,NSQ)
      CALL XtoX0(USER(N33),AOCTyyy,NSQ)
      CALL XtoX0(USER(N34),AOCTzzz,NSQ)
      CALL XtoX0(USER(N35),AOCTxxy,NSQ)
      CALL XtoX0(USER(N36),AOCTxxz,NSQ)
      CALL XtoX0(USER(N37),AOCTxyy,NSQ)
      CALL XtoX0(USER(N38),AOCTyyz,NSQ)
      CALL XtoX0(USER(N39),AOCTxzz,NSQ)
      CALL XtoX0(USER(N40),AOCTyzz,NSQ)
      CALL XtoX0(USER(N41),AOCTxyz,NSQ)
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE PASSUSEROCT

!======================================================================!
!                                                                      !
!               N O F    E N E R G Y    G R A D I E N T S              !
!                                                                      !
!        Module written by Ion Mitxelena (ion.mitxelena@ehu.eus)       !
!                                                                      !
! ==================================================================== !
!                                                                      !
!           Subroutines for NOF energy gradient calculations           !
!                 ( J. Chem. Phys. 146, 014102, 2017 )                 !
!                                                                      !
!   OPTIMIZE: Routine to carry out geometry optimization               !
!   OPTCGNAG: Optimize geometries using the CG NAG routine E04DGF      !
!   ENERGYFUN: Routine to compute energy needed by OPTCGNAG            !
!   OPTLBFGS: Optimize geometries with LBFGS method (mthlib.f)         !
!   PNOFGRAD: Main driver to compute PNOF gradient at one geometry     !
!   STVDERNOF: Main driver to compute one-electron pnof gradient       !
!   JKDERNOF: Main driver to compute two-electron hf gradient          !
!                                                                      !
!   AOLAGRAN: Calculate PNOF lagrangian in AO basis                    !
!                                                                      !
!   SQUARETRIAN2: Put square matrix (HF density matrix) into diagonal  !
!              form without avoiding double counting of diagonal terms !                 
!   SQUARETRIAN3: Put square matrix (PNOF density matrix) to diagonal  !
!              form without avoiding double counting of diagonal terms !                                                    
!   TRACEs: Calculate the trace of lagrangian and one-electron         !
!           derivative integrals                                       !
!   VINTNOF: Gauss-Hermite quadrature for 1e derivative integrals      !
!   SETHERMITE: Roots and weights for Gauss-Hermite quadrature         !
!   DERINOF: Auxiliar computing 1e derivatives                         !
!   DVINTNOF: Gauss-Hermite quadrature for 1e derivative integrals     !
!   DTXYZNOF: Auxiliar computing 1e derivatives                        !
!   RT123NOF,ROOT4NOF,ROOT5NOF,ROOT6NOF: Auxiliars for roots           !
!   RYSASYNOF,RYSGWNOF,RYSDSNOF: Auxiliars for roots                   !
!                                                                      !
!   OEDHNDNOF: Set up pointers for 1e charge distribution              !
!   OEDRDNOF: Set pointes to the ij or kl charge distributions         !
!   JKDATMNOF: Select centers for derivatives                          !
!   JKDSHLNOF: Select indices for shell block                          !
!   JKDNDXNOF: Select indices for shell block                          !
!   DABCLUNOF: Obtain HF 2e density for one shell block                !
!   DABNOF: Obtain PNOF 2e density for one shell block, N**6           !
!   DABNOF2PRE: Contract CJ12 and CK12 with density matrix to later    !
!               use in DABNOF2                                         !
!   DABNOF2: Obtain PNOF 2e density for one shell block, N**5          !
!   DABNOF5: Obtain PNOF 2e density for one shell block, M*N**4, for   !
!      PNOF5 (DABNOF5PRE) and PNOF7 (DABNOF7PRE) by using separability !
!   JKDSPDNOF: Evaluate derivative integral                            !
!   JKDINVNOF: Process derivative gradient and add to total gradient   !
!   JKWRYSNOF: Compute roots and weights for quadrature                !
!   JKBCDFNOF: Compute coefficients for recursion formulae             !
!   JKGNMSNOF: Compute x,y,z integrals (2 centers)                     !
!   JKXYZSNOF: Compute x,y,z integrals (4 centers)                     !
!   JDXYZSNOF: Compute x,y,z integrals for derivatives                 !
!   DSPDFSNOF: Compute derivative integrals, in principle these        !
!          integrals are added to the gradient 'on the fly', but each  !
!          integral contribution can be stored just removing a 'return'!
!                                                                      !
!   VNNDERNOF: Nuclear repulsion contribution                          !
!   SDERNOF: Density force contribution                                !
!   HELFEYNOF: Hellmann-Feynmann contribution                          !
!   TVDERNOF: 1e contribution due to AO derivatives with respect to    !
!             nuclei motion                                            !
!                                                                      !
!======================================================================!
!-----------------------------------------------------------------------
!              DRIVERS TO CARRY OUT GEOMETRY OPTIMIZATION
!-----------------------------------------------------------------------
! OPTIMIZE
      SUBROUTINE OPTIMIZE(NINTEG,IDONTW,NAT,ZAN,Cxyz,IAN,IMIN,IMAX,    &
                 ZMASS,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,    &
                 ISH,ITYP,C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI,IRUNTYP)
      USE PARCOM            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN,ZMASS
      INTEGER,DIMENSION(NAT):: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT):: Cxyz
      DOUBLE PRECISION,DIMENSION(3*NAT):: GRADS
      DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE:: HESSIANO,DDM
      DOUBLE PRECISION,DIMENSION(3):: DIPS
!-----------------------------------------------------------------------
      EELEC_MIN = 1.0d20
!     Generate an initial GCF if RESTART=F setting ICOEF=0
      IF (RESTART .eqv. .FALSE.) THEN
       ICOEF2 = ICOEF
       ICOEF = 0
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NSHELL,NPRIMI,        &
                    ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,         &
                    KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP, &
                    CD,CF,CG,CH,CI,GRADS,IRUNTYP,DIPS,0,0)
       RESTART = .TRUE.
       ICOEF = ICOEF2
      ENDIF
      IF(ICGMETHOD==1)THEN
       CALL OPTSUMSL(NINTEG,IDONTW,IEMOM,NAT,NBF,NSHELL,NPRIMI,ZAN,    &
                     Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,KATOM,KTYPE,KLOC, &
                     INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,  &
                     CF,CG,CH,CI,IRUNTYP,GRADS,DIPS)           
      ELSE IF(ICGMETHOD==2)THEN
       CALL OPTCGNAG(NINTEG,IDONTW,IEMOM,NAT,NBF,NSHELL,NPRIMI,ZAN,    &
                     Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,KATOM,KTYPE,KLOC, &
                     INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,  &
                     CF,CG,CH,CI,IRUNTYP,GRADS,DIPS)
      ELSE IF(ICGMETHOD==3)THEN
       CALL OPTLBFGS(NINTEG,IDONTW,IEMOM,NAT,NBF,NSHELL,NPRIMI,ZAN,    &
                     Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,       &
                     INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,  &
                     CF,CG,CH,CI,IRUNTYP,GRADS,DIPS,NPRINT)
      ENDIF
!     Compute Hessian from analytic gradients at stationary point
      ALLOCATE (HESSIANO(3*NAT,3*NAT),DDM(9,NAT))
      CALL HSSNUMd(HESSIANO,3*NAT,Cxyz,GRADS,DIPS,DDM,NINTEG,IDONTW,   &
                   ZAN,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,    &
                   KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,CF,CG,    &
                   CH,CI,IRUNTYP)
!     Carry out normal mode vibrational analysis
      CALL FGMTRXd(Cxyz,HESSIANO,GRADS,3*NAT,ZAN,ZMASS,DDM)
!-----------------------------------------------------------------------
      DEALLOCATE (HESSIANO,DDM)
      RETURN
      END
! OPTSUMSL
      SUBROUTINE OPTSUMSL(NINTEG,IDONTW,IEMOM,NAT,NBF,NSHELL,NPRIMI,   &
                          ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,KATOM,   &
                          KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,  &
                          C2,EX,CS,CP,CD,CF,CG,CH,CI,IRUNTYP,GRADIENT, &
                          DIPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      PARAMETER (BOHR = 0.52917724924D+00) 
      INTEGER,DIMENSION(NAT)    :: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
      INTEGER,INTENT(IN)        :: IRUNTYP
      DOUBLE PRECISION,DIMENSION(NAT)       :: ZAN,ZMASS
      DOUBLE PRECISION,DIMENSION(NPRIMI)    :: C1,C2,EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT)     :: Cxyz
      DOUBLE PRECISION,DIMENSION(3*NAT),INTENT(OUT) :: GRADIENT
      DOUBLE PRECISION,DIMENSION(3) :: DIPS
      COMMON/EHFEN/EHF,EN,EMP2,EMP3
      COMMON/ENERGY/EELEC_OLD,EELEC,DIF_EELEC,EELEC_MIN
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9, &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,    &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,    &
                         NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,NU7,NU8,NU9,  &
                         NU10,NU11,NU12,NU13,NULAST
      COMMON/GRADIENTE/GRADS(300)               ! NATOMS=100 in donof.f   
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
      WRITE(6,'(/,1X,A23,7X,A19,/)')'Call in CG Optimization',         &
                                    'Total Energy (a.u.)'
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
      GRADIENT(1:NV) = GRADS(1:NV)
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
      CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),               &
                 IUSER(NIU5),IUSER(NIU6),IUSER(NIU7),IUSER(NIU8),      &
                 USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),              &
                 IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),               &
                 IUSER(NIU14),IUSER(NIU15),IUSER(NIU16),IUSER(NIU17),  &
                 IUSER(NIU18),IUSER(NIU19),IUSER(NIU20),IUSER(NIU21),  &
                 USER(NU3),USER(NU4),USER(NU5),USER(NU6),USER(NU7),    &
                 USER(NU8),USER(NU9),USER(NU10),USER(NU11),USER(NU12), &
                 GRADS,IUSER(NIU22),USER(NU13),0,1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Final Coordinates on File CGGRAD (Unit=100)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(11,'(/,72(1H*),/)')
      WRITE(11,*)'New Coordinates after Conjugate Gradient Opt (Angs)'
      WRITE(11,*)
      DO I=1,NAT
       WRITE(11,'(I5,3F15.4)')   &
            I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
      ENDDO
      CALL NUCDIST(NV,Cxyz)
      WRITE(11,'(/,72(1H*),/)')      
      GRADIENT(1:NV) = GRADS(1:NV)
!-----------------------------------------------------------------------
      DEALLOCATE(IUSER,USER,D,IV,V)
      RETURN
      END
! POINTERSOPT
      SUBROUTINE POINTERSOPT(NAT,NSHELL,NPRIMI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9, &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,    &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,    &
                         NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,NU7,NU8,NU9,  &
                         NU10,NU11,NU12,NU13,NULAST
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
      NIULAST = NIU23 + 1               
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
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9, &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,    &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,    &
                         NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,NU7,NU8,NU9,  &
                         NU10,NU11,NU12,NU13,NULAST
      INTEGER,DIMENSION(NIULAST) :: IUSER 
      DOUBLE PRECISION,DIMENSION(NV)         :: Cxyz
      DOUBLE PRECISION,DIMENSION(NULAST)     :: USER
      COMMON/GRADIENTE/GRADS(300)
!-----------------------------------------------------------------------
      CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),               &
                 IUSER(NIU5),IUSER(NIU6),IUSER(NIU7),IUSER(NIU8),      &
                 USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),              &
                 IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),               &
                 IUSER(NIU14),IUSER(NIU15),IUSER(NIU16),IUSER(NIU17),  &
                 IUSER(NIU18),IUSER(NIU19),IUSER(NIU20),IUSER(NIU21),  &
                 USER(NU3),USER(NU4),USER(NU5),USER(NU6),USER(NU7),    &
                 USER(NU8),USER(NU9),USER(NU10),USER(NU11),USER(NU12), &
                 GRADS,IUSER(NIU22),USER(NU13),0,0)
      ENERGIA = EELEC + EN     
      WRITE(6,'(8X,I3,16X,F20.10)')NF,ENERGIA
!-----------------------------------------------------------------------
      RETURN
      END
! CALCOPTG      
      SUBROUTINE CALCOPTG(NV,Cxyz,NF,GRAD,IUSER,USER)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9, &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,    &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,    &
                         NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,NU7,NU8,NU9,  &
                         NU10,NU11,NU12,NU13,NULAST
      INTEGER,DIMENSION(NIULAST) :: IUSER 
      DOUBLE PRECISION,DIMENSION(NV)         :: Cxyz
      DOUBLE PRECISION,DIMENSION(NV)         :: GRAD      
      DOUBLE PRECISION,DIMENSION(NULAST)     :: USER
      COMMON/GRADIENTE/GRADS(300)
!-----------------------------------------------------------------------
!     Avoiding warnings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NF = NF
      Cxyz(1) = Cxyz(1)
      IUSER(1) = IUSER(1)
      USER(1) = USER(1)
!
      GRAD(1:NV) = GRADS(1:NV)
!-----------------------------------------------------------------------
      RETURN
      END      
! OPTCGNAG
      SUBROUTINE OPTCGNAG(NINTEG,IDONTW,IEMOM,NAT,NBF,NSHELL,NPRIMI,   &
                          ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,KATOM,   &
                          KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,  &
                          C2,EX,CS,CP,CD,CF,CG,CH,CI,IRUNTYP,GRADIENT, &
                          DIPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      PARAMETER (BOHR = 0.52917724924D+00) 
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN,ZMASS
      INTEGER,DIMENSION(NAT):: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: C1,C2,EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT):: Cxyz
      DOUBLE PRECISION,DIMENSION(3*NAT),INTENT(OUT):: GRADIENT
      INTEGER,INTENT(IN):: IRUNTYP
      DOUBLE PRECISION,DIMENSION(3):: DIPS
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9, &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,    &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,    &
                         NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,NU7,NU8,NU9,  &
                         NU10,NU11,NU12,NU13,NULAST
!                         
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IUSER,IWORK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)    :: USER,WORK,GRADS
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
      WRITE(6,'(/,1X,A23,7X,A19,/)')'Call in CG Optimization',         &
                                    'Total Energy (a.u.)'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Send output of E04DGF to NAG file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!nag      CALL X04ABF(1,2)      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Maximum Step Length (Default = 10**20))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!nag      CALL E04DKF ('Maximum Step Length  = 1.0')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Optimality Tolerance (relative precision**0.8)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!nag      CALL E04DKF ('Optimality Tolerance = 1.0D-10')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Linesearch Tolerance (0<r<1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!nag      CALL E04DKF ('Linesearch Tolerance = 1.0D-3')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Function Precision (machine precision**0.9)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!nag      CALL E04DKF ('Function Precision = 1.0D-10')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Verify Level (-1 = No checks, 0 = cheap test, 1 = 0 + gradients)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!nag      CALL E04DKF ('Verify Level = -1')      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Print Level (0 = No output, 1 = The final solution only)
!                 (5 = One line of summary output for each iteration)
!                 (10 = The final solution and one line for each iter.)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!nag      CALL E04DKF ('Print Level = 5')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calling to NAG Library for using the CG method
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NV = 3*NAT
      IFAIL = 1
      EMINIMA = 1.0d6
      ALLOCATE(IWORK(NV+1),WORK(13*NV),GRADS(NV))
!nag      CALL E04DGF(NV,ENERGYFUN,ITER_E04DGF,ENERGIA,GRADS,Cxyz,         &
!nag                  IWORK,WORK,IUSER,USER,IFAIL)
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
      CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),               &
                 IUSER(NIU5),IUSER(NIU6),IUSER(NIU7),IUSER(NIU8),      &
                 USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),              &
                 IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),               &
                 IUSER(NIU14),IUSER(NIU15),IUSER(NIU16),IUSER(NIU17),  &
                 IUSER(NIU18),IUSER(NIU19),IUSER(NIU20),IUSER(NIU21),  &
                 USER(NU3),USER(NU4),USER(NU5),USER(NU6),USER(NU7),    &
                 USER(NU8),USER(NU9),USER(NU10),USER(NU11),USER(NU12), &
                 GRADS,IUSER(NIU22),USER(NU13),0,1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Final Coordinates on File CGGRAD (Unit=100)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(11,'(/,72(1H*),/)')
      WRITE(11,*)'New Coordinates after Conjugate Gradient Opt (Angs)'
      WRITE(11,*)
      DO I=1,NAT
       WRITE(11,'(I5,3F15.4)')   &
            I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
      ENDDO
      CALL NUCDIST(NV,Cxyz)
      WRITE(11,'(/,72(1H*),/)')      
      GRADIENT = GRADS
!-----------------------------------------------------------------------
      DEALLOCATE(IUSER,USER,GRADS,IWORK,WORK)
      RETURN
      END
! ENERGYFUN
      SUBROUTINE ENERGYFUN(MODE,NV,Cxyz,ENERGIA,GRADS,NSTATE,IUSER,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9, &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,    &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,    &
                         NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,NU7,NU8,NU9,  &
                         NU10,NU11,NU12,NU13,NULAST
      COMMON/EHFEN/EHF,EN,EMP2,EMP3
      COMMON/ENERGY/EELEC_OLD,EELEC,DIF_EELEC,EELEC_MIN
      INTEGER,DIMENSION(*) :: IUSER
      DOUBLE PRECISION,DIMENSION(*) :: USER
      DOUBLE PRECISION,DIMENSION(NV) :: Cxyz,GRADS
      PARAMETER (BOHR = 0.52917724924D+00)
!-----------------------------------------------------------------------
      NSTATE = 1
      MODE = 1
      IUSER(NIU23) = IUSER(NIU23) + 1                  ! ITCG = ITCG + 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),               &
                 IUSER(NIU5),IUSER(NIU6),IUSER(NIU7),IUSER(NIU8),      &
                 USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),              &
                 IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),               &
                 IUSER(NIU14),IUSER(NIU15),IUSER(NIU16),IUSER(NIU17),  &
                 IUSER(NIU18),IUSER(NIU19),IUSER(NIU20),IUSER(NIU21),  &
                 USER(NU3),USER(NU4),USER(NU5),USER(NU6),USER(NU7),    &
                 USER(NU8),USER(NU9),USER(NU10),USER(NU11),USER(NU12), &
                 GRADS,IUSER(NIU22),USER(NU13),0,0)
      ENERGIA = EELEC + EN     
      WRITE(6,'(8X,I3,16X,F20.10)')IUSER(NIU23),ENERGIA
!-----------------------------------------------------------------------
      RETURN
      END
! OPTLBFGS
      SUBROUTINE OPTLBFGS(NINTEG,IDONTW,IEMOM,NAT,NBF,NSHELL,NPRIMI,   &
                 ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,            &
                 KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,     &
                 CD,CF,CG,CH,CI,IRUNTYP,GRADIENT,DIPS,NPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (BOHR = 0.52917724924D+00)
      INTEGER,INTENT(IN)::NINTEG,IDONTW,IEMOM,NAT,NBF,NSHELL
      INTEGER,INTENT(IN)::NPRIMI,IRUNTYP,NPRINT
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
      INTEGER,DIMENSION(NAT):: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: C1,C2,EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT):: Cxyz
      DOUBLE PRECISION,DIMENSION(3*NAT),INTENT(OUT):: GRADIENT
      DOUBLE PRECISION,DIMENSION(3*NAT):: GRADS
      DOUBLE PRECISION,DIMENSION(3):: DIPS
      COMMON/EHFEN/EHF,EN,EMP2,EMP3
      COMMON/ENERGY/EELEC_OLD,EELEC,DIF_EELEC,EELEC_MIN
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
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
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NSHELL,NPRIMI,        &
                     ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,        &
                     KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP, &
                     CD,CF,CG,CH,CI,GRADS,IRUNTYP,DIPS,0,0)
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
!      We allow at most 2000 evaluations of F and G
       IF(ICALL.GT.2000) EXIT
      ENDDO
      DEALLOCATE(W)
!     FINAL CALL TO COMPUTE ENERGY AND GRADS AT EQUIL. GEOMETRY
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
      CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NSHELL,NPRIMI,         &
                   ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,          &
                   KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,   &
                   CD,CF,CG,CH,CI,GRADS,IRUNTYP,DIPS,0,1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Final Coordinates on File CGGRAD (Unit=100)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(11,*)
      WRITE(11,*)'New Coordinates after the LBFGS Method (Angs)'
      WRITE(11,*)
      DO I=1,NAT
       WRITE(11,'(I5,3F20.10)')    &
            I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
      ENDDO
      CALL NUCDIST(N,Cxyz)
      GRADIENT = GRADS
!-----------------------------------------------------------------------
      RETURN
      END
!-----------------------------------------------------------------------
!               MAIN ROUTINES TO COMPUTE PNOF GRADIENT
!-----------------------------------------------------------------------
! PNOFGRAD
      SUBROUTINE PNOFGRAD(COEF,QD,RO,ELAG,GRADS,ATMNAME,KATOM,KTYPE,   &
                          KLOC,KKMIN,KKMAX,KSTART,KNG,CX0,CY0,CZ0,ZNUC,&
                          EX1,CS,CP,CD,CF,CG,CJ12,CK12,XINTS,IPRINTOPT)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      CHARACTER*4 ATMNAME(NATOMS)
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KATOM,KTYPE,KLOC,KKMIN,KKMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KSTART,KNG
      INTEGER,INTENT(IN)::IPRINTOPT
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN)::CX0,CY0,CZ0,ZNUC
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN)::EX1,CS,CP,CD,CF,CG
      DOUBLE PRECISION,DIMENSION(NBF,NBF),INTENT(IN)::COEF
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF),INTENT(IN)::QD
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF),INTENT(IN)::ELAG
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5),INTENT(IN)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION((NSHELL*NSHELL+NSHELL)/2),INTENT(IN)::XINTS
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT):: GRADS
!
!     ONE-ELECTRON CONTRIBUTION TO THE GRADIENT
!
!      CALL STVDERNOF(COEF,QD,RO,ELAG,GRADS,ATMNAME, &
      CALL STVDERNOF(COEF,QD,RO,ELAG,GRADS,                     &
                     KATOM,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG,   &
                     CX0,CY0,CZ0,ZNUC,EX1,CS,CP,CD,CF,CG)
!
!     TWO-ELECTRON CONTRIBUTION TO THE GRADIENT
!
      IF((IPNOF==5 .or. IPNOF==7) .and. NSOC==0) THEN
       CALL JKDERNOF5(KATOM,KTYPE,KLOC,KKMIN,KKMAX,             &
                 KSTART,KNG,QD,RO,GRADS,CX0,CY0,CZ0,            &
                 EX1,CS,CP,CD,CF,CG,ATMNAME,XINTS,IPRINTOPT)
      ELSE
       CALL JKDERNOF(KATOM,KTYPE,KLOC,KKMIN,KKMAX,                 &
                     KSTART,KNG,CJ12,CK12,QD,RO,GRADS,CX0,CY0,CZ0, &
                     EX1,CS,CP,CD,CF,CG,ATMNAME,XINTS,IPRINTOPT)
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END      
! STVDERNOF
      SUBROUTINE STVDERNOF(COEF,QD,RO,ELAG,GRADS,                      &
                           KATOM,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG,    &
                           CX0,CY0,CZ0,ZNUC,EX1,CS,CP,CD,CF,CG)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KATOM,KTYPE,KLOC,KKMIN,KKMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KSTART,KNG  
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN)::CX0,CY0,CZ0,ZNUC
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN)::EX1,CS,CP,CD,CF,CG
      DOUBLE PRECISION,DIMENSION(NBF,NBF),INTENT(IN)::COEF
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF),INTENT(IN)::QD
      DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::DM2
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF),INTENT(IN)::ELAG
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT):: GRADS
      DOUBLE PRECISION,DIMENSION(NBFT)::LEPS
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00
!-----------------------------------------------------------------------
!     MAKE GRADIENT ZERO
      GRADS=ZERO
!-----------------------------------------------------------------------
!     Nuclear repulsion force
!-----------------------------------------------------------------------
      CALL VNNDERNOF(CX0,CY0,CZ0,ZNUC,GRADS)
!-----------------------------------------------------------------------
!     Density force contribution
!-----------------------------------------------------------------------      
!     GET LAGRANGIAN MATRIX (LEPS)    
      CALL AOLAGRAN(COEF,ELAG,LEPS,NBF,NBFT)
!     COMPUTE AND ADD OVERLAP DERIVATIVES
      CALL SDERNOF(KATOM,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG,  &
                CX0,CY0,CZ0,EX1,CS,CP,CD,CF,CG,LEPS,GRADS)
!-----------------------------------------------------------------------
!     One-electron hamiltonian contribution
!-----------------------------------------------------------------------
!     CONTRACT OVER OCCUPATIONS TO OBTAIN DENSITY MATRIX
      ALLOCATE(DM2(NBF,NBF))
      DM2=ZERO
      DO J=1,NBF
       DO L=1,J
        DO I=1,NBF5
          DM2(J,L)=DM2(J,L)+QD(I,J,L)*RO(I)
        ENDDO
        DM2(L,J)=DM2(J,L)
       ENDDO
      ENDDO
      DM2=DM2+DM2
!      
!     HELLMANN-FEYNMAN FORCE      
!
      CALL HELFEYNOF(KATOM,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG,  &
                  CX0,CY0,CZ0,EX1,CS,CP,CD,CF,CG,ZNUC,DM2,GRADS)
!                  
!     INTEGRAL FORCE (AO DERIVATIVE CONTRIBUTION)                  
!
      CALL TVDERNOF(KATOM,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG,  &
                  CX0,CY0,CZ0,EX1,CS,CP,CD,CF,CG,ZNUC,DM2,GRADS)
!                  
      DEALLOCATE(DM2)
!     
!     MAKE ZERO GRADIENT CORRESPONDING TO FROZEN COORDINATES
!
      IF(FROZEN) THEN
        DO I=1,200,2
         IF(IFROZEN(I).EQ.0) EXIT
         GRADS(IFROZEN(I),IFROZEN(I+1))=ZERO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!    1 FORMAT( /,5X,'-----------------------',                 &
!              /5X,' One-Electron Gradient ',                  &
!              /5X,'-----------------------',                  &
!              //2X,'Atom',10X,'EX',13X,'EY',13X,'EZ' )
!    2 FORMAT(/I4,1X,A4,3X,F10.7,5X,F10.7,5X,F10.7)
!-----------------------------------------------------------------------
      RETURN
      END   
! JKDERNOF
      SUBROUTINE JKDERNOF(KATOM,KTYPE,KLOC,KKMIN,KKMAX,                &
                          KSTART,KNG,CJ12,CK12,QD,RO,GRADS,CX0,CY0,CZ0,&
                          EX1,CS,CP,CD,CF,CG,ATMNAME,XINTS,IPRINTOPT)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
!     INPUT-OUTPUT VARIABLES OR ARGUMENTS
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KATOM,KTYPE,KLOC,KKMIN,KKMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KSTART,KNG
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5),INTENT(IN)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF),INTENT(IN)::QD
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO
      DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::P,DAAUX,DAAUX2
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN)::CX0,CY0,CZ0
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN)::EX1,CS,CP,CD,CF,CG
      DOUBLE PRECISION,DIMENSION((NSHELL*NSHELL+NSHELL)/2),INTENT(IN)::XINTS
      CHARACTER*4 ATMNAME(NATOMS)
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT)::GRADS
      INTEGER,INTENT(IN)::IPRINTOPT
!     INTERMEDIATE VARIABLES 
      INTEGER,DIMENSION(NBF)::IA
      INTEGER,DIMENSION(4,35)::IGXYZ,JGXYZ,KGXYZ,LGXYZ
      INTEGER,DIMENSION(:,:),ALLOCATABLE::NIJG,IJKLG
      INTEGER,DIMENSION(5)::LENSHL
      INTEGER::MAXFUN,NIJGDIM,MAXVEC,MINVEC,MODTYP,LVAR,LFIX
      INTEGER::MINMEM,MAXNUM,MINXYZ,MAXXYZ
      INTEGER::IDER,JDER,KDER,LDER
      INTEGER::TOTCOUNT,INVTYP,NROOTS,NIJ,NKL,KANG,I
      INTEGER::II,JJ,KK,LL,MAXLL,ISH,JSH,KSH,LSH,IIAT,JJAT,KKAT,LLAT
      INTEGER::MINJ,MAXJ,MINK,MAXK,MINL,MAXL,MINI,MAXI,NIJ0,NKL0
      INTEGER::LIT,LJT,LKT,LLTT,NUMI,NUMJ,NUMK,NUML,IIJJ,KKLL
      LOGICAL::SKIPI,SKIPJ,SKIPK,SKIPL,SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      LOGICAL::IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL,EXPNDI,EXPNDK
      DOUBLE PRECISION::DABMAX,DABCUT,CUTOFF,CUTOFF2,GMAX
      DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::DCHRG
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::DAB,GINT,FINT
      DOUBLE PRECISION,DIMENSION(84)::PNRM
      DOUBLE PRECISION::FI
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00,ONE=1.0D+00
      DOUBLE PRECISION,PARAMETER::SQRT3=1.73205080756888D+00
      DOUBLE PRECISION,PARAMETER::SQRT5=2.23606797749979D+00
      DOUBLE PRECISION,PARAMETER::SQRT7=2.64575131106459D+00
      DOUBLE PRECISION,PARAMETER::TEN=10.0D+00,TENM9=1.0D-09,TENM11=1.0D-11
      DATA LENSHL /1,4,10,20,35/
!-----------------------------------------------------------------------
!     ROUTINE BASCHK DO:
      KANG=0
      DO N=1,NSHELL
        IF(KTYPE(N).GT.KANG) KANG = KTYPE(N)
      ENDDO
      MINXYZ=(4*KANG -2 + 1)/2
      MODTYP=KANG+1
      MAXVEC= 255/3
!     ROUTINE IGETGRDVECLEN DO:
      MINVEC= MAXVEC*3+1
!
!     NOT IMPLEMENTED FOR SHELLS LARGER THAN G      
      IF(KANG.GT.5) THEN
        WRITE(6,*) "G IS THE MAXIMUM SHELL FOR GRADIENTS"
        STOP
      ENDIF      
!
!     SET NORMALIZATION CONSTANTS
!     MAXFUN=NUMBER OF FUNCTIONS WITH ANGULAR MOMENTUM EQUAL TO MAXTYP
      MAXFUN=LENSHL(KANG)
      DO I=1,MAXFUN
        PNRM(I)=ONE
      ENDDO
      SQRT53=SQRT5/SQRT3
      FI=ZERO
      DO I=1,MAXFUN
        IF (I.EQ.1.OR.I.EQ.2.OR.I.EQ.5.OR.I.EQ.11.OR.I.EQ.21) THEN
           FI = ONE
        ELSE IF (I.EQ.8.OR.I.EQ.20.OR.I.EQ.33) THEN
           FI=FI*SQRT3
        ELSE IF (I.EQ.14) THEN
           FI=FI*SQRT5
        ELSE IF (I.EQ.24) THEN
           FI=FI*SQRT7
        ELSE IF (I.EQ.30) THEN
           FI=FI*SQRT53
        END IF
        PNRM(I)=FI
      ENDDO
!
!     SET STARTING PARAMETERS
!
!     CUTOFF IS THE SCHWARZ SCREENING CUT OFF 
      CUTOFF=TENM9/TEN
      CUTOFF2=CUTOFF*0.5D+00
      ZBIG = 0.0D+00
      DO ISH=1,NSHELL
         I1=KSTART(ISH)
         I2=I1+KNG(ISH)-1
         DO IG=I1,I2
            IF(EX1(IG).GT.ZBIG) ZBIG = EX1(IG)           
         ENDDO
      ENDDO
!     DABCUT IS THE TWO PARTICLE DENSITY CUT OFF      
      DABCUT=1.0D-11
      IF(ZBIG.GT.1.0D+06) DABCUT = DABCUT/10.0D+00
      IF(ZBIG.GT.1.0D+07) DABCUT = DABCUT/10.0D+00
!     SET POINTERS FOR PARTITIONING MEMORY
      DO I=1,NBF
        IA(I) = (I*I-I)/2
      ENDDO   
!           
      NIJGDIM=0
!     CALCULATE THE LARGEST SHELL TYPE      
      DO N=1,NSHELL
!       GET NUMBER OF PRIMITIVE CHARGE DISTRIBUTIONS        
        DO NN2=1,N
          NIJGDIM=NIJGDIM+KNG(N)*KNG(NN2)
        ENDDO
      ENDDO
!      
!     MAXNUM=NUMBER OF FUNCTIONS WITH ANGULAR MOMENTUM EQUAL TO MAXTYP
      MAXNUM=((KANG)*(KANG+1))/2
!     DO AT LEAST AN L SHELL      
      IF(MAXNUM.LT.4) MAXNUM=4
      MAXNUM=(MAXNUM**4)
!      
      LVAR=     ( MODTYP**2       * MODTYP**2       )*3                 
      LVAR=LVAR+( MODTYP**2       *(MODTYP+MODTYP-1))*3                 
      LVAR=LVAR+((MODTYP+MODTYP-1)*(MODTYP+MODTYP-1))*3                 
      LVAR=LVAR+( MODTYP**2                         )*3                 
      LVAR=LVAR+((MODTYP+MODTYP-1)                  )*3                 
      LVAR=LVAR+((MODTYP+MODTYP-1)* 3               )*3                 
      LVAR=LVAR+(  3                                )*3                 
      LVAR=LVAR+(  9                                )                   
      LVAR=LVAR+(  4                                )                   
      LVAR=LVAR+(  5                                )                   
      LVAR=LVAR+( 18                                )                   
      LVAR=LVAR+(  2                                )                   
      LVAR=LVAR+(  4                                )*3                 
      LVAR=LVAR+( MODTYP**2       * MODTYP**2       )*3*14
      LFIX=(NBF*NBF + NBF)/2
      LFIX=LFIX  +(NSHELL*(NSHELL+1))/2
      LFIX=LFIX+( (NSHELL*(NSHELL+1))/2 )*2
      LFIX=LFIX+MAXNUM* 4                                               
      LFIX=LFIX+MAXNUM                                                  
      LFIX=LFIX+NIJGDIM*15
!
!     CHECK OUT NEEDED MEMORY
!     MAXXYZ=MAXIMUM NUMBER OF PRIMITIVE INTEGRALS THAT CAN BE HANDLED IN ONE VECTOR
      MAXXYZ=MAXVEC
      MINMEM=(MINXYZ*LVAR)+1+LFIX
      IF(MAXXYZ.LT.MINXYZ) THEN
        WRITE(6,3) MINMEM
        STOP
      ENDIF                                
!
      ALLOCATE(IJKLG(4,MAXNUM))
      ALLOCATE(NIJG(2,(NSHELL*(NSHELL+1))/2))
      ALLOCATE(DCHRG(15,NIJGDIM))
      ALLOCATE(DAB(MAXNUM))
      ALLOCATE(GINT(MAXNUM))
      ALLOCATE(FINT(MAXNUM))
!     SET UP 1e CHARGE DISTRIBUTION
      CALL OEDHNDNOF(DCHRG,NIJGDIM,NIJG,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG, &
                  KATOM,IA,CX0,CY0,CZ0,EX1,CS,CP,CD,CF,CG)          
      ALLOCATE(P(NBF,NBFT))
      CALL SQUARETRIAN3(QD,P,NBF,NBFT)
      P=P*0.5D+00
      ALLOCATE(DAAUX(NBF,NBFT))
      ALLOCATE(DAAUX2(NBF,NBFT))
      CALL DABNOF2PRE(CJ12,CK12,RO,P,DAAUX,DAAUX2)
!
      TOTCOUNT=1
!
!----I SHELL
!
      DO II=1,NSHELL
!      
!-----J SHELL
!
       DO JJ=1,II
        IIJJ=IA(MAX0(II,JJ))+MIN0(II,JJ)
!       SET POINTERS TO THE IJ CHARGE DISTRIBUTION        
        CALL OEDRDNOF(NIJG,NIJ,NIJ0,IIJJ)
        IF(NIJ.EQ.0) CYCLE
!        
!-----K SHELL
!
        DO KK=1,II
!        
!-----L SHELL
!
        MAXLL=KK
        IF(KK.EQ.II) MAXLL=JJ
         DO LL=1,MAXLL
!
!         IMPLEMENT INTEGRAL SCREENING HERE USING EXCHANGE INTEGRALS
!                                                                       
          IJIJ=IA(MAX0(II,JJ))+MIN0(II,JJ)                          
          KLKL=IA(MAX0(KK,LL))+MIN0(KK,LL)                          
          GMAX=(XINTS(IJIJ)*XINTS(KLKL))
!
!         COARSE SCREENING, ON JUST THE INTEGRAL VALUE         
!         ONLY WORKS IF SCHWARZ SCREENING IS ON (NATOMS>5)
          IF (GMAX.LT.CUTOFF.AND.NATOMS>5) CYCLE
!
          ISH=II
          JSH=JJ
          KSH=KK
          LSH=LL
          KKLL=IA(MAX0(KK,LL))+MIN0(KK,LL)
!         SET POINTERS TO THE KL CHARGE DISTRIBUTION                  
          CALL OEDRDNOF(NIJG,NKL,NKL0,KKLL)
          IF(NKL.EQ.0) CYCLE
!          
!         SELECT CENTERS FOR DERIVATIVES
!
          CALL JKDATMNOF(ISH,JSH,KSH,LSH,SKIPI,SKIPJ,SKIPK,SKIPL,      &
                         INVTYP,KATOM,IIAT,JJAT,KKAT,LLAT)                      
          IF(SKIPI.AND.SKIPJ.AND.SKIPK.AND.SKIPL) CYCLE
!          
!         SELECT INDICES FOR SHELL BLOCK
!
          CALL JKDSHLNOF(ISH,JSH,KSH,LSH,KTYPE,KLOC,KKMIN,KKMAX,       &
                      SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL,EXPNDI,EXPNDK   &
                      ,MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL,        &
                      IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL,              &
                      LIT,LJT,LKT,LLTT,NUMI,NUMJ,NUMK,NUML)                      
          CALL JKDNDXNOF(LIT,LJT,LKT,LLTT,NUMJ,NUMK,NUML,              &
                      SKIPI,SKIPJ,SKIPK,SKIPL,MINI,MAXI,               &
                      MINJ,MAXJ,MINK,MAXK,MINL,MAXL,NROOTS,            &
                      IGXYZ,JGXYZ,KGXYZ,LGXYZ,                         &
                      IIEQJJ,KKEQLL,IJEQKL,IJKLG,MAXNUM,               &
                      IDER,JDER,KDER,LDER)          
!                      
!         OBTAIN 2e DENSITY FOR THIS SHELL BLOCK                      
!
          CALL DABNOF2(ISH,JSH,KSH,LSH,PNRM,KKMIN,KKMAX,KLOC,          &
                      IGXYZ,JGXYZ,KGXYZ,LGXYZ,IIEQJJ,KKEQLL,IJEQKL,    &
                      IA,P,DAB,MAXNUM,DABMAX,DAAUX,DAAUX2)
!
!         FINE SCREENING, ON INTEGRAL VALUE TIMES DENSITY FACTOR
!         ONLY WORKS IF SCHWARZ SCREENING IS ON (NATOMS>5)
          IF(DABMAX*GMAX.LT.CUTOFF2.AND.NATOMS>5) CYCLE
!
!         EVALUATE DERIVATIVE INTEGRAL AND ADD TO THE GRADIENT
!
          CALL JKDSPDNOF(TOTCOUNT,GINT,FINT,DAB,DCHRG,NIJGDIM,         &
                      MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL,         &
                      IIEQJJ,KKEQLL,IJEQKL,EXPNDI,EXPNDK,              &
                      SKIPI,SKIPJ,SKIPK,SKIPL,SPI,SPJ,SPK,SPL,SPIJKL,  &
                      IJKLG,MAXNUM,INVTYP,NROOTS,MINVEC,MODTYP,MAXXYZ, &
                      NKL,NIJ,IIAT,JJAT,KKAT,LLAT,GRADS,DABMAX,DABCUT, &
                      NKL0,NIJ0,LIT,LJT,LKT,LLTT,IDER,JDER,KDER,LDER)               
!      
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      DEALLOCATE(P,IJKLG,NIJG,DCHRG,DAB,GINT,FINT,DAAUX,DAAUX2)
!     
!     MAKE ZERO GRADIENT CORRESPONDING TO FROZEN COORDINATES
!
      IF(FROZEN) THEN
        DO I=1,200,2
         IF(IFROZEN(I).EQ.0) EXIT
         GRADS(IFROZEN(I),IFROZEN(I+1))=ZERO
        ENDDO
      ENDIF      
!
!     PRINT OUT TOTAL PNOF ENERGY GRADIENT
!
      IF(IPRINTOPT==1)THEN
       WRITE(6,1)
       DO I=1,NATOMS
        WRITE(6,2)I,ATMNAME(I),GRADS(1,I),GRADS(2,I),GRADS(3,I)
       ENDDO
      ENDIF
!-----------------------------------------------------------------------
    1 FORMAT( /1X,'----------------',                                  &
              /1X,' Total Gradient ',                                  &
              /1X,'----------------',                                  &
              //9X,'Atom',7X,'Ex',10X,'Ey',10X,'Ez' )
    2 FORMAT(/1X,I4,5X,A4,F10.4,2X,F10.4,2X,F10.4)
    3 FORMAT(/,' THE MINIMUM MEMORY REQUIRED IS ',I10,' WORDS.')
!-----------------------------------------------------------------------
      RETURN
      END
! JKDERNOF5
      SUBROUTINE JKDERNOF5(KATOM,KTYPE,KLOC,KKMIN,KKMAX,               &
                           KSTART,KNG,QD,RO,GRADS,CX0,CY0,CZ0,         &
                           EX1,CS,CP,CD,CF,CG,ATMNAME,XINTS,IPRINTOPT)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
!     INPUT-OUTPUT VARIABLES OR ARGUMENTS
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KATOM,KTYPE,KLOC,KKMIN,KKMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KSTART,KNG
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF),INTENT(IN)::QD
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO
      DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::P,CJAUX,CKAUX
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN)::CX0,CY0,CZ0
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN)::EX1,CS,CP,CD,CF,CG
      DOUBLE PRECISION,DIMENSION((NSHELL*NSHELL+NSHELL)/2),INTENT(IN)::XINTS
      CHARACTER*4 ATMNAME(NATOMS)
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT)::GRADS
      INTEGER,INTENT(IN)::IPRINTOPT
!     INTERMEDIATE VARIABLES
      INTEGER,DIMENSION(NBF)::IA
      INTEGER,DIMENSION(4,35)::IGXYZ,JGXYZ,KGXYZ,LGXYZ
      INTEGER,DIMENSION(:,:),ALLOCATABLE::NIJG,IJKLG
      INTEGER,DIMENSION(5)::LENSHL
      INTEGER::MAXFUN,NIJGDIM,MAXVEC,MINVEC,MODTYP,LVAR,LFIX
      INTEGER::MINMEM,MAXNUM,MINXYZ,MAXXYZ
      INTEGER::IDER,JDER,KDER,LDER
      INTEGER::TOTCOUNT,INVTYP,NROOTS,NIJ,NKL,KANG,I
      INTEGER::II,JJ,KK,LL,MAXLL,ISH,JSH,KSH,LSH,IIAT,JJAT,KKAT,LLAT
      INTEGER::MINJ,MAXJ,MINK,MAXK,MINL,MAXL,MINI,MAXI,NIJ0,NKL0
      INTEGER::LIT,LJT,LKT,LLTT,NUMI,NUMJ,NUMK,NUML,IIJJ,KKLL
      LOGICAL::SKIPI,SKIPJ,SKIPK,SKIPL,SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      LOGICAL::IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL,EXPNDI,EXPNDK
      DOUBLE PRECISION::DABMAX,DABCUT,CUTOFF,CUTOFF2,GMAX
      DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::DCHRG
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::DAB,GINT,FINT
      DOUBLE PRECISION,DIMENSION(84)::PNRM
      DOUBLE PRECISION::FI
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00,ONE=1.0D+00
      DOUBLE PRECISION,PARAMETER::SQRT3=1.73205080756888D+00
      DOUBLE PRECISION,PARAMETER::SQRT5=2.23606797749979D+00
      DOUBLE PRECISION,PARAMETER::SQRT7=2.64575131106459D+00
      DOUBLE PRECISION,PARAMETER::TEN=10.0D+00,TENM9=1.0D-09,TENM11=1.0D-11
      DATA LENSHL /1,4,10,20,35/
!-----------------------------------------------------------------------
!     ROUTINE BASCHK DO:
      KANG=0
      DO N=1,NSHELL
        IF(KTYPE(N).GT.KANG) KANG = KTYPE(N)
      ENDDO
      MINXYZ=(4*KANG -2 + 1)/2
      MODTYP=KANG+1
      MAXVEC= 255/3
!     ROUTINE IGETGRDVECLEN DO:
      MINVEC= MAXVEC*3+1
!
!     NOT IMPLEMENTED FOR SHELLS LARGER THAN G      
      IF(KANG.GT.5) THEN
        WRITE(6,*) "G IS THE MAXIMUM SHELL FOR GRADIENTS"
        STOP
      ENDIF      
!
!     SET NORMALIZATION CONSTANTS
!     MAXFUN=NUMBER OF FUNCTIONS WITH ANGULAR MOMENTUM EQUAL TO MAXTYP
      MAXFUN=LENSHL(KANG)
      DO I=1,MAXFUN
        PNRM(I)=ONE
      ENDDO
      SQRT53=SQRT5/SQRT3
      FI=ZERO
      DO I=1,MAXFUN
        IF (I.EQ.1.OR.I.EQ.2.OR.I.EQ.5.OR.I.EQ.11.OR.I.EQ.21) THEN
           FI = ONE
        ELSE IF (I.EQ.8.OR.I.EQ.20.OR.I.EQ.33) THEN
           FI=FI*SQRT3
        ELSE IF (I.EQ.14) THEN
           FI=FI*SQRT5
        ELSE IF (I.EQ.24) THEN
           FI=FI*SQRT7
        ELSE IF (I.EQ.30) THEN
           FI=FI*SQRT53
        END IF
        PNRM(I)=FI
      ENDDO
!
!     SET STARTING PARAMETERS
!
!     CUTOFF IS THE SCHWARZ SCREENING CUT OFF 
      CUTOFF=TENM9/TEN
      CUTOFF2=CUTOFF*0.5D+00
      ZBIG = 0.0D+00
      DO ISH=1,NSHELL
         I1=KSTART(ISH)
         I2=I1+KNG(ISH)-1
         DO IG=I1,I2
            IF(EX1(IG).GT.ZBIG) ZBIG = EX1(IG)           
         ENDDO
      ENDDO
!     DABCUT IS THE TWO PARTICLE DENSITY CUT OFF      
      DABCUT=1.0D-11
      IF(ZBIG.GT.1.0D+06) DABCUT = DABCUT/10.0D+00
      IF(ZBIG.GT.1.0D+07) DABCUT = DABCUT/10.0D+00
!     SET POINTERS FOR PARTITIONING MEMORY
      DO I=1,NBF
        IA(I) = (I*I-I)/2
      ENDDO   
!           
      NIJGDIM=0
!     CALCULATE THE LARGEST SHELL TYPE      
      DO N=1,NSHELL
!       GET NUMBER OF PRIMITIVE CHARGE DISTRIBUTIONS        
        DO NN2=1,N
          NIJGDIM=NIJGDIM+KNG(N)*KNG(NN2)
        ENDDO
      ENDDO
!      
!     MAXNUM=NUMBER OF FUNCTIONS WITH ANGULAR MOMENTUM EQUAL TO MAXTYP
      MAXNUM=((KANG)*(KANG+1))/2
!     DO AT LEAST AN L SHELL      
      IF(MAXNUM.LT.4) MAXNUM=4
      MAXNUM=(MAXNUM**4)
!      
      LVAR=     ( MODTYP**2       * MODTYP**2       )*3                 
      LVAR=LVAR+( MODTYP**2       *(MODTYP+MODTYP-1))*3                 
      LVAR=LVAR+((MODTYP+MODTYP-1)*(MODTYP+MODTYP-1))*3                 
      LVAR=LVAR+( MODTYP**2                         )*3                 
      LVAR=LVAR+((MODTYP+MODTYP-1)                  )*3                 
      LVAR=LVAR+((MODTYP+MODTYP-1)* 3               )*3                 
      LVAR=LVAR+(  3                                )*3                 
      LVAR=LVAR+(  9                                )                   
      LVAR=LVAR+(  4                                )                   
      LVAR=LVAR+(  5                                )                   
      LVAR=LVAR+( 18                                )                   
      LVAR=LVAR+(  2                                )                   
      LVAR=LVAR+(  4                                )*3                 
      LVAR=LVAR+( MODTYP**2       * MODTYP**2       )*3*14
      LFIX=(NBF*NBF + NBF)/2
      LFIX=LFIX  +(NSHELL*(NSHELL+1))/2
      LFIX=LFIX+( (NSHELL*(NSHELL+1))/2 )*2
      LFIX=LFIX+MAXNUM* 4
      LFIX=LFIX+MAXNUM
      LFIX=LFIX+NIJGDIM*15
!
!     CHECK OUT NEEDED MEMORY
!     MAXXYZ=MAXIMUM NUMBER OF PRIMITIVE INTEGRALS THAT CAN BE HANDLED IN ONE VECTOR
      MAXXYZ=MAXVEC
      MINMEM=(MINXYZ*LVAR)+1+LFIX
      IF(MAXXYZ.LT.MINXYZ) THEN
        WRITE(6,3) MINMEM
        STOP
      ENDIF                                
!
      ALLOCATE(IJKLG(4,MAXNUM))
      ALLOCATE(NIJG(2,(NSHELL*(NSHELL+1))/2))
      ALLOCATE(DCHRG(15,NIJGDIM))
      ALLOCATE(DAB(MAXNUM))
      ALLOCATE(GINT(MAXNUM))
      ALLOCATE(FINT(MAXNUM))
!     SET UP 1e CHARGE DISTRIBUTION
      CALL OEDHNDNOF(DCHRG,NIJGDIM,NIJG,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG, &
                  KATOM,IA,CX0,CY0,CZ0,EX1,CS,CP,CD,CF,CG)          
      ALLOCATE(P(NBF,NBFT))  
      CALL SQUARETRIAN3(QD,P,NBF,NBFT)
      P=P*0.5D+00
      ALLOCATE(CJAUX(NBFT,NBFT))
      ALLOCATE(CKAUX(NBFT,NBFT))
      IF(IPNOF==5) CALL DABNOF5PRE(RO,P,CJAUX,CKAUX)
      IF(IPNOF==7) CALL DABNOF7PRE(RO,P,CJAUX,CKAUX)
      DEALLOCATE(P)
!
      TOTCOUNT=1
!
!----I SHELL
!
      DO II=1,NSHELL
!      
!-----J SHELL
!
       DO JJ=1,II
        IIJJ=IA(MAX0(II,JJ))+MIN0(II,JJ)
!       SET POINTERS TO THE IJ CHARGE DISTRIBUTION        
        CALL OEDRDNOF(NIJG,NIJ,NIJ0,IIJJ)
        IF(NIJ.EQ.0) CYCLE
!        
!-----K SHELL
!
        DO KK=1,II
!        
!-----L SHELL
!
        MAXLL=KK
        IF(KK.EQ.II) MAXLL=JJ
         DO LL=1,MAXLL
!
!         IMPLEMENT INTEGRAL SCREENING HERE USING EXCHANGE INTEGRALS
!                                                                       
          IJIJ=IA(MAX0(II,JJ))+MIN0(II,JJ)                          
          KLKL=IA(MAX0(KK,LL))+MIN0(KK,LL)                          
          GMAX=(XINTS(IJIJ)*XINTS(KLKL))
!
!         COARSE SCREENING, ON JUST THE INTEGRAL VALUE         
!         ONLY WORKS IF SCHWARZ SCREENING IS ON (NATOMS>5)
          IF (GMAX.LT.CUTOFF.AND.NATOMS>5) CYCLE
!
          ISH=II
          JSH=JJ
          KSH=KK
          LSH=LL
          KKLL=IA(MAX0(KK,LL))+MIN0(KK,LL)
!         SET POINTERS TO THE KL CHARGE DISTRIBUTION                  
          CALL OEDRDNOF(NIJG,NKL,NKL0,KKLL)
          IF(NKL.EQ.0) CYCLE
!          
!         SELECT CENTERS FOR DERIVATIVES
!
          CALL JKDATMNOF(ISH,JSH,KSH,LSH,SKIPI,SKIPJ,SKIPK,SKIPL,INVTYP,  &
                      KATOM,IIAT,JJAT,KKAT,LLAT)                      
          IF(SKIPI.AND.SKIPJ.AND.SKIPK.AND.SKIPL) CYCLE
!          
!         SELECT INDICES FOR SHELL BLOCK
!
          CALL JKDSHLNOF(ISH,JSH,KSH,LSH,KTYPE,KLOC,KKMIN,KKMAX,       &
                      SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL,EXPNDI,EXPNDK   &
                      ,MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL,        &
                      IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL,              &
                      LIT,LJT,LKT,LLTT,NUMI,NUMJ,NUMK,NUML)                      
          CALL JKDNDXNOF(LIT,LJT,LKT,LLTT,NUMJ,NUMK,NUML,              &
                      SKIPI,SKIPJ,SKIPK,SKIPL,MINI,MAXI,               &
                      MINJ,MAXJ,MINK,MAXK,MINL,MAXL,NROOTS,            &
                      IGXYZ,JGXYZ,KGXYZ,LGXYZ,                         &
                      IIEQJJ,KKEQLL,IJEQKL,IJKLG,MAXNUM,               &
                      IDER,JDER,KDER,LDER)          
!                      
!         OBTAIN 2e DENSITY FOR THIS SHELL BLOCK                      
!
          CALL DABNOF5(ISH,JSH,KSH,LSH,PNRM,KKMIN,KKMAX,KLOC,          &
                      IGXYZ,JGXYZ,KGXYZ,LGXYZ,IIEQJJ,KKEQLL,IJEQKL,    &
                      IA,DAB,MAXNUM,DABMAX,CJAUX,CKAUX)
!
!         FINE SCREENING, ON INTEGRAL VALUE TIMES DENSITY FACTOR
!         ONLY WORKS IF SCHWARZ SCREENING IS ON (NATOMS>5)
          IF(DABMAX*GMAX.LT.CUTOFF2.AND.NATOMS>5) CYCLE
!
!         EVALUATE DERIVATIVE INTEGRAL AND ADD TO THE GRADIENT
!
          CALL JKDSPDNOF(TOTCOUNT,GINT,FINT,DAB,DCHRG,NIJGDIM,         &
                      MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL,         &
                      IIEQJJ,KKEQLL,IJEQKL,EXPNDI,EXPNDK,              &
                      SKIPI,SKIPJ,SKIPK,SKIPL,SPI,SPJ,SPK,SPL,SPIJKL,  &
                      IJKLG,MAXNUM,INVTYP,NROOTS,MINVEC,MODTYP,MAXXYZ, &
                      NKL,NIJ,IIAT,JJAT,KKAT,LLAT,GRADS,DABMAX,DABCUT, &
                      NKL0,NIJ0,LIT,LJT,LKT,LLTT,IDER,JDER,KDER,LDER)               
!      
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      DEALLOCATE(IJKLG,NIJG,DCHRG,DAB,GINT,FINT,CJAUX,CKAUX)
!     
!     MAKE ZERO GRADIENT CORRESPONDING TO FROZEN COORDINATES
!
      IF(FROZEN) THEN
        DO I=1,200,2
         IF(IFROZEN(I).EQ.0) EXIT
         GRADS(IFROZEN(I),IFROZEN(I+1))=ZERO
        ENDDO
      ENDIF      
!
!     PRINT OUT TOTAL PNOF ENERGY GRADIENT
!
      IF(IPRINTOPT==1)THEN
       WRITE(6,1)
       DO I=1,NATOMS
        WRITE(6,2)I,ATMNAME(I),GRADS(1,I),GRADS(2,I),GRADS(3,I)
       ENDDO
      ENDIF
!-----------------------------------------------------------------------
    1 FORMAT( /1X,'----------------',                                  &
              /1X,' Total Gradient ',                                  &
              /1X,'----------------',                                  &
              //9X,'Atom',7X,'Ex',10X,'Ey',10X,'Ez' )
    2 FORMAT(/1X,I4,5X,A4,F10.4,2X,F10.4,2X,F10.4)
    3 FORMAT(/,' THE MINIMUM MEMORY REQUIRED IS ',I10,' WORDS.')
!-----------------------------------------------------------------------
      RETURN
      END      
!-----------------------------------------------------------------------
!            ROUTINES TO GET THE LAGRANGIAN MATRIX FOR PNOF
!-----------------------------------------------------------------------
! AOLAGRAN
      SUBROUTINE AOLAGRAN(C,EE,LEPS,N,NL)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN)::C,EE
      DOUBLE PRECISION,DIMENSION(NL),INTENT(OUT)::LEPS
      DOUBLE PRECISION,DIMENSION(N,N)::AUXELG,AELG
      DOUBLE PRECISION::AUX
      DOUBLE PRECISION,PARAMETER::PT5=0.5D+00,ZERO=0.0D+00
!-----------------------------------------------------------------------
!     COEF*LAMBDA*COEF
!-----------------------------------------------------------------------
      LEPS=ZERO
!     SIMETRIZATION OF LAMBDA MATRIX      
      DO I=1,NBF
       DO J=1,I
        AUXELG(I,J)=EE(I,J)
        AUXELG(J,I)=AUXELG(I,J)
       ENDDO
      ENDDO
!     SCALING AS N**3 CALCULATION, FIRST CONTRACTING OVER A SINGLE COEF
      AELG=ZERO
      DO L=1,N
       DO M=1,NBF
        DO NM=1,NBF
          AELG(M,L)=AELG(M,L)+AUXELG(M,NM)*C(L,NM)
        ENDDO
       ENDDO
      ENDDO
      KL=0
      DO K=1,N
       DO L=1,K
        KL=KL+1
        AUX=ZERO
        DO M=1,NBF
          AUX = AUX - AELG(M,L)*C(K,M)
        ENDDO
        LEPS(KL) = AUX + AUX
       ENDDO
!     DOUBLE COUNTING OF DIAGONAL IS AVOIDED IN TRACEs              
      ENDDO
      LEPS = LEPS + LEPS
!-----------------------------------------------------------------------
      RETURN
      END
!-----------------------------------------------------------------------
!                         MATHEMATICAL ROUTINES
!-----------------------------------------------------------------------
! SQUARETRIAN2
      SUBROUTINE SQUARETRIAN2(FM,F,N,NT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(N,N)::FM
      DOUBLE PRECISION,DIMENSION(NT)::F
!-----------------------------------------------------------------------
!     PUT SQUARE FM IN TRIANGULAR FORM F
!-----------------------------------------------------------------------
      NZ=0
      DO I=1,N
       DO J=1,I
        NZ=NZ+1
        F(NZ)=FM(I,J)+FM(I,J)
       ENDDO
!      DOUBLE COUNTING OF DIAGONAL TERMS IS SKIPPED IN DABCLU
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
! SQUARETRIAN3
      SUBROUTINE SQUARETRIAN3(FM,F,N,NT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(N,N,N)::FM
      DOUBLE PRECISION,DIMENSION(N,NT)::F
      NZ=0
      DO I=1,N
       DO J=1,I
        NZ=NZ+1
        F(:,NZ)=FM(:,I,J)+FM(:,J,I)
       ENDDO
!      AVOIDING DOUBLE DIAGONAL COUNTING IN DABNOF       
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END      
! TRACEs
      SUBROUTINE TRACEs(E,AA,EEPS,N,NL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN)::AA
      DOUBLE PRECISION,DIMENSION(NL),INTENT(IN)::EEPS
      DOUBLE PRECISION,INTENT(OUT)::E
      E = 0.0D+0
      KL = 0
      DO K=1,N
       DO L=1,K
        KL = KL + 1
        E = E + AA(K,L)*EEPS(KL)
       ENDDO
       E = E - AA(K,K)*EEPS(KL)*0.5d0
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
! VINTNOF
      SUBROUTINE VINTNOF(NI,NJ,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,           &
                         XINT,YINT,ZINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER,INTENT(IN)::NI,NJ
      DOUBLE PRECISION,INTENT(IN)::T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ
      DOUBLE PRECISION,INTENT(OUT)::XINT,YINT,ZINT
      INTEGER,DIMENSION(7)::IIMIN,IIMAX
      DOUBLE PRECISION,DIMENSION(55)::HTOTAL,WTOTAL
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00
!
      DATA IIMIN /1,2,4,7,11,16,22/
      DATA IIMAX /1,3,6,10,15,21,28/
!
!     ----- GAUSS-HERMITE QUADRATURE USING MINIMUM POINT FORMULA -----
!
      CALL SETHERMITE(HTOTAL,WTOTAL)
      XINT = ZERO
      YINT = ZERO
      ZINT = ZERO
      NPTS = (NI+NJ-2)/2+1
      IMIN = IIMIN(NPTS)
      IMAX = IIMAX(NPTS)
      DO I = IMIN,IMAX
         DUM = WTOTAL(I)
         PX = DUM
         PY = DUM
         PZ = DUM
         DUM = HTOTAL(I)*T
         PTX = DUM+X0
         PTY = DUM+Y0
         PTZ = DUM+Z0
         AX = PTX-XI
         AY = PTY-YI
         AZ = PTZ-ZI
         BX = PTX-XJ
         BY = PTY-YJ
         BZ = PTZ-ZJ
         IF(PX+AX .EQ. PX) AX = ZERO
         IF(PY+AY .EQ. PY) AY = ZERO
         IF(PZ+AZ .EQ. PZ) AZ = ZERO
         IF(PX+BX .EQ. PX) BX = ZERO
         IF(PY+BY .EQ. PY) BY = ZERO
         IF(PZ+BZ .EQ. PZ) BZ = ZERO
         GO TO (180,164, 162, 160,140,120,100),NI
!
  100    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  120    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  140    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  160    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  162    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  164    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  180    GO TO (320,310,300,280,260,240,220,200),NJ
!
  200    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  220    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  240    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  260    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  280    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  300    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  310    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  320    CONTINUE
         XINT = XINT+PX
         YINT = YINT+PY
         ZINT = ZINT+PZ
      ENDDO
!-----------------------------------------------------------------------  
      RETURN
      END
! DERINOF
      SUBROUTINE DERINOF(DXDI,DYDI,DZDI,X,Y,Z,LIT,LJT,AI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(5,1),INTENT(OUT)::DXDI,DYDI,DZDI
      DOUBLE PRECISION,DIMENSION(6,5),INTENT(IN)::X,Y,Z
      DOUBLE PRECISION,INTENT(IN)::AI
      INTEGER,INTENT(IN)::LIT,LJT
!
!     ----- DXDI ... -----
!
      DO J=1,LJT
        DXDI(1,J)=X(2,J)*(AI+AI)
        DYDI(1,J)=Y(2,J)*(AI+AI)
        DZDI(1,J)=Z(2,J)*(AI+AI)
      ENDDO
!
      IF(LIT.EQ.1) RETURN
!
      DO I=2,LIT
        DO J=1,LJT
          DXDI(I,J)=X(I+1,J)*(AI+AI) - X(I-1,J)*(I-1)
          DYDI(I,J)=Y(I+1,J)*(AI+AI) - Y(I-1,J)*(I-1)
          DZDI(I,J)=Z(I+1,J)*(AI+AI) - Z(I-1,J)*(I-1)
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
! DVINTNOF
      SUBROUTINE DVINTNOF(NI,NJ,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,          &
                          CX,CY,CZ,XINT,YINT,ZINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER,INTENT(IN)::NI,NJ
      DOUBLE PRECISION,INTENT(IN)::T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,CX,CY,CZ
      DOUBLE PRECISION,INTENT(OUT)::XINT,YINT,ZINT
      INTEGER,DIMENSION(7)::IIMIN,IIMAX
      DOUBLE PRECISION,DIMENSION(55)::HTOTAL,WTOTAL
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00
!
      DATA IIMIN /1,2,4,7,11,16,22/
      DATA IIMAX /1,3,6,10,15,21,28/
!
!     ----- GAUSS-HERMITE QUADRATURE USING MINIMUM POINT FORMULA -----
!
      CALL SETHERMITE(HTOTAL,WTOTAL)      
      XINT = ZERO
      YINT = ZERO
      ZINT = ZERO
      NPTS = (NI+NJ+1-2)/2+1
      IMIN = IIMIN(NPTS)
      IMAX = IIMAX(NPTS)
      DO I = IMIN,IMAX
         DUM = HTOTAL(I)*T
         PTX = DUM+X0
         PTY = DUM+Y0
         PTZ = DUM+Z0
         PX = PTX-CX
         PY = PTY-CY
         PZ = PTZ-CZ
         AX = PTX-XI
         AY = PTY-YI
         AZ = PTZ-ZI
         BX = PTX-XJ
         BY = PTY-YJ
         BZ = PTZ-ZJ
         GO TO (180,170,160,140,120,100),NI
!
  100    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  120    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  140    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  160    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  170    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  180    GO TO (320,300,280,260,240,220,200),NJ
!
  200    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  220    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  240    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  260    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  280    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  300    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  320    DUM = WTOTAL(I)
         XINT = XINT+DUM*PX
         YINT = YINT+DUM*PY
         ZINT = ZINT+DUM*PZ
      ENDDO
!-----------------------------------------------------------------------  
      RETURN
      END      
! DTXYZNOF
      SUBROUTINE DTXYZNOF(XT,YT,ZT,XS,YS,ZS,NI,NJ,AJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,INTENT(IN)::NI,NJ
      DOUBLE PRECISION,INTENT(IN)::AJ
      DOUBLE PRECISION,DIMENSION(6,5),INTENT(OUT)::XT,YT,ZT
      DOUBLE PRECISION,DIMENSION(6,7),INTENT(IN)::XS,YS,ZS
      DOUBLE PRECISION,PARAMETER::THREE=3.0D+00
!
      DO I=1,NI
        XT(I,1)=(XS(I,1  )       - XS(I,3  )*(AJ+AJ))*AJ
        YT(I,1)=(YS(I,1  )       - YS(I,3  )*(AJ+AJ))*AJ
        ZT(I,1)=(ZS(I,1  )       - ZS(I,3  )*(AJ+AJ))*AJ
      ENDDO
!
      IF(NJ.EQ.1) RETURN
!
      DO I=1,NI
        XT(I,2)=(XS(I,2  )*THREE - XS(I,4  )*(AJ+AJ))*AJ
        YT(I,2)=(YS(I,2  )*THREE - YS(I,4  )*(AJ+AJ))*AJ
        ZT(I,2)=(ZS(I,2  )*THREE - ZS(I,4  )*(AJ+AJ))*AJ
      ENDDO
!
      IF(NJ.EQ.2) RETURN
!
      DO J=3,NJ
        FACT1 = J+J-1
        FACT2 = (J-1)*(J-2)
        FACT2 = FACT2/2
        DO I=1,NI
          XT(I,J)=(XS(I,J  )*FACT1 - XS(I,J+2)*(AJ+AJ))*AJ   &
                                   - XS(I,J-2)*FACT2
          YT(I,J)=(YS(I,J  )*FACT1 - YS(I,J+2)*(AJ+AJ))*AJ   &
                                   - YS(I,J-2)*FACT2
          ZT(I,J)=(ZS(I,J  )*FACT1 - ZS(I,J+2)*(AJ+AJ))*AJ   &
                                   - ZS(I,J-2)*FACT2
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
! RT123NOF
      SUBROUTINE RT123NOF(X,NROOTS,UAUX,WAUX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION,INTENT(IN)::X
      INTEGER,INTENT(IN)::NROOTS
      DOUBLE PRECISION,DIMENSION(13),INTENT(OUT)::UAUX,WAUX
      DOUBLE PRECISION,DIMENSION(13)::U,W
!      
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
!
      DATA R12,PIE4/2.75255128608411D-01, 7.85398163397448D-01/
      DATA R22,W22/ 2.72474487139158D+00, 9.17517095361369D-02/
      DATA R13/     1.90163509193487D-01/
      DATA R23,W23/ 1.78449274854325D+00, 1.77231492083829D-01/
      DATA R33,W33/ 5.52534374226326D+00, 5.11156880411248D-03/
!-----------------------------------------------------------------------     
!      
      IF (X .GT. 5.0D+00) GO TO 400
      IF (X .GT. 1.0D+00) GO TO 280
      IF (X .GT. 3.0D-07) GO TO 180
!     X IS APPROXIMATELY ZERO.         NROOTS=1,2, OR 3
      IF (NROOTS-2) 120,140,160
  120 RT1 = 0.5D+00 -X/5.0D+00
      WW1 = 1.0D+00 -X/3.0D+00
      UAUX=U
      WAUX=W
      RETURN
  140 RT1 = 1.30693606237085D-01 -2.90430236082028D-02 *X
      RT2 = 2.86930639376291D+00 -6.37623643058102D-01 *X
      WW1 = 6.52145154862545D-01 -1.22713621927067D-01 *X
      WW2 = 3.47854845137453D-01 -2.10619711404725D-01 *X
      UAUX=U
      WAUX=W
      RETURN
  160 RT1 = 6.03769246832797D-02 -9.28875764357368D-03 *X
      RT2 = 7.76823355931043D-01 -1.19511285527878D-01 *X
      RT3 = 6.66279971938567D+00 -1.02504611068957D+00 *X
      WW1 = 4.67913934572691D-01 -5.64876917232519D-02 *X
      WW2 = 3.60761573048137D-01 -1.49077186455208D-01 *X
      WW3 = 1.71324492379169D-01 -1.27768455150979D-01 *X
      UAUX=U
      WAUX=W
      RETURN
!     X = 0.0 TO 1.0                   NROOTS=1,2, OR 3
  180 IF (NROOTS .EQ. 3) GO TO 220
      F1 = ((((((((-8.36313918003957D-08*X+1.21222603512827D-06 )*X- &
          1.15662609053481D-05 )*X+9.25197374512647D-05 )*X-        &
          6.40994113129432D-04 )*X+3.78787044215009D-03 )*X-        &
          1.85185172458485D-02 )*X+7.14285713298222D-02 )*X-        &
          1.99999999997023D-01 )*X+3.33333333333318D-01
      WW1 = (X+X)*F1+EXP(-X)
      IF (NROOTS .EQ. 2) GO TO 200
      RT1 = F1/(WW1-F1)
      UAUX=U
      WAUX=W
      RETURN
  200 RT1 = (((((((-2.35234358048491D-09*X+2.49173650389842D-08)*X-  &
          4.558315364581D-08)*X-2.447252174587D-06)*X+              &
          4.743292959463D-05)*X-5.33184749432408D-04 )*X+           &
          4.44654947116579D-03 )*X-2.90430236084697D-02 )*X+        &
          1.30693606237085D-01
      RT2 = (((((((-2.47404902329170D-08*X+2.36809910635906D-07)*X+  &
          1.835367736310D-06)*X-2.066168802076D-05)*X-              &
          1.345693393936D-04)*X-5.88154362858038D-05 )*X+           &
          5.32735082098139D-02 )*X-6.37623643056745D-01 )*X+        &
          2.86930639376289D+00
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      UAUX=U
      WAUX=W
      RETURN
  220 RT1 = ((((((-5.10186691538870D-10*X+2.40134415703450D-08)*X-   &
          5.01081057744427D-07 )*X+7.58291285499256D-06 )*X-        &
          9.55085533670919D-05 )*X+1.02893039315878D-03 )*X-        &
          9.28875764374337D-03 )*X+6.03769246832810D-02             
      RT2 = ((((((-1.29646524960555D-08*X+7.74602292865683D-08)*X+   &
          1.56022811158727D-06 )*X-1.58051990661661D-05 )*X-        &
          3.30447806384059D-04 )*X+9.74266885190267D-03 )*X-        &
          1.19511285526388D-01 )*X+7.76823355931033D-01             
      RT3 = ((((((-9.28536484109606D-09*X-3.02786290067014D-07)*X-   &
         2.50734477064200D-06 )*X-7.32728109752881D-06 )*X+        &
         2.44217481700129D-04 )*X+4.94758452357327D-02 )*X-        &
         1.02504611065774D+00 )*X+6.66279971938553D+00             
      F2 = ((((((((-7.60911486098850D-08*X+1.09552870123182D-06 )*X- &
          1.03463270693454D-05 )*X+8.16324851790106D-05 )*X-        &
          5.55526624875562D-04 )*X+3.20512054753924D-03 )*X-        &
          1.51515139838540D-02 )*X+5.55555554649585D-02 )*X-        &
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
      UAUX=U
      WAUX=W
      RETURN
  280 IF (X .GT. 3.0D+00) GO TO 340
!     X = 1.0 TO 3.0                   NROOTS=1,2, OR 3
      Y = X-2.0D+00
      IF (NROOTS .EQ. 3) GO TO 320
      F1 = ((((((((((-1.61702782425558D-10*Y+1.96215250865776D-09 )*Y- &
           2.14234468198419D-08 )*Y+2.17216556336318D-07 )*Y-         &
           1.98850171329371D-06 )*Y+1.62429321438911D-05 )*Y-         &
           1.16740298039895D-04 )*Y+7.24888732052332D-04 )*Y-         &
           3.79490003707156D-03 )*Y+1.61723488664661D-02 )*Y-         &
           5.29428148329736D-02 )*Y+1.15702180856167D-01              
      WW1 = (X+X)*F1+EXP(-X)                                           
      IF (NROOTS .EQ. 2) GO TO 300                                     
      RT1 = F1/(WW1-F1)
      UAUX=U
      WAUX=W
      RETURN
  300 RT1 = (((((((((-6.36859636616415D-12*Y+8.47417064776270D-11)*Y-  &
           5.152207846962D-10)*Y-3.846389873308D-10)*Y+               &
           8.472253388380D-08)*Y-1.85306035634293D-06 )*Y+            &
           2.47191693238413D-05 )*Y-2.49018321709815D-04 )*Y+         &
           2.19173220020161D-03 )*Y-1.63329339286794D-02 )*Y+         &
           8.68085688285261D-02                                       
      RT2 = ((((((((( 1.45331350488343D-10*Y+2.07111465297976D-09)*Y-  &
           1.878920917404D-08)*Y-1.725838516261D-07)*Y+               &
           2.247389642339D-06)*Y+9.76783813082564D-06 )*Y-            &
           1.93160765581969D-04 )*Y-1.58064140671893D-03 )*Y+         &
           4.85928174507904D-02 )*Y-4.30761584997596D-01 )*Y+         &
           1.80400974537950D+00                                       
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      UAUX=U
      WAUX=W
      RETURN
  320 RT1 = (((((((( 1.44687969563318D-12*Y+4.85300143926755D-12)*Y-   &
           6.55098264095516D-10 )*Y+1.56592951656828D-08 )*Y-          &
           2.60122498274734D-07 )*Y+3.86118485517386D-06 )*Y-          &
           5.13430986707889D-05 )*Y+6.03194524398109D-04 )*Y-          &
           6.11219349825090D-03 )*Y+4.52578254679079D-02
      RT2 = ((((((( 6.95964248788138D-10*Y-5.35281831445517D-09)*Y-    &
           6.745205954533D-08)*Y+1.502366784525D-06)*Y+                &
           9.923326947376D-07)*Y-3.89147469249594D-04 )*Y+             &
           7.51549330892401D-03 )*Y-8.48778120363400D-02 )*Y+          &
           5.73928229597613D-01
      RT3 = ((((((((-2.81496588401439D-10*Y+3.61058041895031D-09)*Y+   &
           4.53631789436255D-08 )*Y-1.40971837780847D-07 )*Y-          &
           6.05865557561067D-06 )*Y-5.15964042227127D-05 )*Y+          &
           3.34761560498171D-05 )*Y+5.04871005319119D-02 )*Y-          &
           8.24708946991557D-01 )*Y+4.81234667357205D+00
      F2 = ((((((((((-1.48044231072140D-10*Y+1.78157031325097D-09 )*Y- &
           1.92514145088973D-08 )*Y+1.92804632038796D-07 )*Y-          &
           1.73806555021045D-06 )*Y+1.39195169625425D-05 )*Y-          &
           9.74574633246452D-05 )*Y+5.83701488646511D-04 )*Y-          &
           2.89955494844975D-03 )*Y+1.13847001113810D-02 )*Y-          &
           3.23446977320647D-02 )*Y+5.29428148329709D-02               
      GO TO 240                                                        
!     X = 3.0 TO 5.0                   NROOTS =1,2, OR 3               
  340 Y = X-4.0D+00
      IF (NROOTS .EQ. 3) GO TO 380
      F1 = ((((((((((-2.62453564772299D-11*Y+3.24031041623823D-10 )*Y- &
           3.614965656163D-09)*Y+3.760256799971D-08)*Y-                &
           3.553558319675D-07)*Y+3.022556449731D-06)*Y-                &
           2.290098979647D-05)*Y+1.526537461148D-04)*Y-                &
           8.81947375894379D-04 )*Y+4.33207949514611D-03 )*Y-          &
           1.75257821619926D-02 )*Y+5.28406320615584D-02            
      WW1 = (X+X)*F1+EXP(-X)                                            
      IF (NROOTS .EQ. 2) GO TO 360                                      
      RT1 = F1/(WW1-F1)
      UAUX=U
      WAUX=W
      RETURN
  360 RT1 = ((((((((-4.11560117487296D-12*Y+7.10910223886747D-11)*Y-   &
              1.73508862390291D-09 )*Y+5.93066856324744D-08 )*Y-       &
              9.76085576741771D-07 )*Y+1.08484384385679D-05 )*Y-       &
              1.12608004981982D-04 )*Y+1.16210907653515D-03 )*Y-       &
              9.89572595720351D-03 )*Y+6.12589701086408D-02
      RT2 = (((((((((-1.80555625241001D-10*Y+5.44072475994123D-10)*Y+  &
           1.603498045240D-08)*Y-1.497986283037D-07)*Y-                &
           7.017002532106D-07)*Y+1.85882653064034D-05 )*Y-             &
           2.04685420150802D-05 )*Y-2.49327728643089D-03 )*Y+          &
           3.56550690684281D-02 )*Y-2.60417417692375D-01 )*Y+          &
           1.12155283108289D+00                                        
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)                   
      WW1 = WW1-WW2                                                     
      UAUX=U
      WAUX=W
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
  400 IF (X .GT. 15.0D+00) GO TO 560                                    
      E = EXP(-X)
      IF (X .GT. 10.0D+00) GO TO 480
!     X = 5.0 TO 10.0                  NROOTS =1,2, OR 3
      WW1 = (((((( 4.6897511375022D-01/X-6.9955602298985D-01)/X +       &
          5.3689283271887D-01)/X-3.2883030418398D-01)/X +              &
          2.4645596956002D-01)/X-4.9984072848436D-01)/X -              &
          3.1501078774085D-06)*E + SQRT(PIE4/X)                        
      F1 = (WW1-E)/(X+X)
      IF (NROOTS-2) 420,440,460
  420 RT1 = F1/(WW1-F1)
      UAUX=U
      WAUX=W
      RETURN
  440 Y = X-7.5D+00
      RT1 = (((((((((((((-1.43632730148572D-16*Y+2.38198922570405D-16)* &
             Y+1.358319618800D-14)*Y-7.064522786879D-14)*Y-            &
             7.719300212748D-13)*Y+7.802544789997D-12)*Y+              &
             6.628721099436D-11)*Y-1.775564159743D-09)*Y+              &
             1.713828823990D-08)*Y-1.497500187053D-07)*Y+              &
             2.283485114279D-06)*Y-3.76953869614706D-05 )*Y+           &
             4.74791204651451D-04 )*Y-4.60448960876139D-03 )*Y+        &
             3.72458587837249D-02                                      
      RT2 = (((((((((((( 2.48791622798900D-14*Y-1.36113510175724D-13)*Y- &
            2.224334349799D-12)*Y+4.190559455515D-11)*Y-                &
            2.222722579924D-10)*Y-2.624183464275D-09)*Y+                &
            6.128153450169D-08)*Y-4.383376014528D-07)*Y-                &
            2.49952200232910D-06 )*Y+1.03236647888320D-04 )*Y-          &
            1.44614664924989D-03 )*Y+1.35094294917224D-02 )*Y-          &
            9.53478510453887D-02 )*Y+5.44765245686790D-01               
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)                    
      WW1 = WW1-WW2
      UAUX=U
      WAUX=W
      RETURN
  460 F2 = (F1+F1+F1-E)/(X+X)
      Y = X-7.5D+00
      RT1 = ((((((((((( 5.74429401360115D-16*Y+7.11884203790984D-16)*Y-  &
           6.736701449826D-14)*Y-6.264613873998D-13)*Y+                  &
           1.315418927040D-11)*Y-4.23879635610964D-11 )*Y+               &
           1.39032379769474D-09 )*Y-4.65449552856856D-08 )*Y+            &
           7.34609900170759D-07 )*Y-1.08656008854077D-05 )*Y+            &
           1.77930381549953D-04 )*Y-2.39864911618015D-03 )*Y+            &
           2.39112249488821D-02                                          
      RT2 = ((((((((((( 1.13464096209120D-14*Y+6.99375313934242D-15)*Y-  &
          8.595618132088D-13)*Y-5.293620408757D-12)*Y-                  &
          2.492175211635D-11)*Y+2.73681574882729D-09 )*Y-               &
          1.06656985608482D-08 )*Y-4.40252529648056D-07 )*Y+            &
          9.68100917793911D-06 )*Y-1.68211091755327D-04 )*Y+            &
          2.69443611274173D-03 )*Y-3.23845035189063D-02 )*Y+            &
          2.75969447451882D-01
      RT3 = (((((((((((( 6.66339416996191D-15*Y+1.84955640200794D-13)*Y- &
           1.985141104444D-12)*Y-2.309293727603D-11)*Y+                  &
           3.917984522103D-10)*Y+1.663165279876D-09)*Y-                  &
           6.205591993923D-08)*Y+8.769581622041D-09)*Y+                  &
           8.97224398620038D-06 )*Y-3.14232666170796D-05 )*Y-            &
           1.83917335649633D-03 )*Y+3.51246831672571D-02 )*Y-            &
           3.22335051270860D-01 )*Y+1.73582831755430D+00                 
      GO TO 260                                                          
!     X = 10.0 TO 15.0                 NROOTS=1,2, OR 3
  480 WW1 = (((-1.8784686463512D-01/X+2.2991849164985D-01)/X -           &
            4.9893752514047D-01)/X-2.1916512131607D-05)*E + SQRT(PIE4/X)
      F1 = (WW1-E)/(X+X)
      IF (NROOTS-2) 500,520,540
  500 RT1 = F1/(WW1-F1)
      UAUX=U
      WAUX=W
      RETURN
  520 RT1 = ((((-1.01041157064226D-05*X+1.19483054115173D-03)*X -    &
           6.73760231824074D-02)*X+1.25705571069895D+00)*X + (((-    &
           8.57609422987199D+03/X+5.91005939591842D+03)/X -          &
           1.70807677109425D+03)/X+2.64536689959503D+02)/X -         &
           2.38570496490846D+01)*E + R12/(X-R12)
      RT2 = ((( 3.39024225137123D-04*X-9.34976436343509D-02)*X -         &
           4.22216483306320D+00)*X + (((-2.08457050986847D+03/X -        &
           1.04999071905664D+03)/X+3.39891508992661D+02)/X -             &
           1.56184800325063D+02)/X+8.00839033297501D+00)*E + R22/(X-R22) 
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      UAUX=U
      WAUX=W
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
          8.579165965128D-13)*Y+6.642452485783D-12)*Y+                &
          4.798806828724D-11)*Y-1.13413908163831D-09 )*Y+             &
          7.08558457182751D-09 )*Y-5.59678576054633D-08 )*Y+          &
          2.51020389884249D-06 )*Y-6.63678914608681D-05 )*Y+          &
          1.11888323089714D-03 )*Y-1.45361636398178D-02 )*Y+          &
          1.65077877454402D-01
      RT3 = (((((((((((( 3.20622388697743D-15*Y-2.73458804864628D-14)*Y- &
           3.157134329361D-13)*Y+8.654129268056D-12)*Y-                  &
           5.625235879301D-11)*Y-7.718080513708D-10)*Y+                  &
           2.064664199164D-08)*Y-1.567725007761D-07)*Y-                  &
           1.57938204115055D-06 )*Y+6.27436306915967D-05 )*Y-            &
           1.01308723606946D-03 )*Y+1.13901881430697D-02 )*Y-            &
           1.01449652899450D-01 )*Y+7.77203937334739D-01                 
      GO TO 260                                                          
  560 IF (X .GT. 33.0D+00) GO TO 660
!     X = 15.0 TO 33.0                 NROOTS=1,2, OR 3
      E = EXP(-X)
      WW1 = (( 1.9623264149430D-01/X-4.9695241464490D-01)/X -        &
            6.0156581186481D-05)*E + SQRT(PIE4/X)
      F1 = (WW1-E)/(X+X)
      IF (NROOTS-2) 580,600,620
  580 RT1 = F1/(WW1-F1)
      UAUX=U
      WAUX=W
      RETURN
  600 RT1 = ((((-1.14906395546354D-06*X+1.76003409708332D-04)*X -   &
           1.71984023644904D-02)*X-1.37292644149838D-01)*X + (-    &
           4.75742064274859D+01/X+9.21005186542857D+00)/X -        &
           2.31080873898939D-02)*E + R12/(X-R12)                   
      RT2 = ((( 3.64921633404158D-04*X-9.71850973831558D-02)*X -    &
           4.02886174850252D+00)*X + (-1.35831002139173D+02/X -    &
           8.66891724287962D+01)/X+2.98011277766958D+00)*E + R22/(X-R22)
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      UAUX=U
      WAUX=W
      RETURN
  620 F2 = (F1+F1+F1-E)/(X+X)
      IF (X .GT. 20.0D+00) GO TO 640
      RT1 = ((((((-2.43270989903742D-06*X+3.57901398988359D-04)*X -  &
           2.34112415981143D-02)*X+7.81425144913975D-01)*X -         &
           1.73209218219175D+01)*X+2.43517435690398D+02)*X + (-      &
           1.97611541576986D+04/X+9.82441363463929D+03)/X -          &
           2.07970687843258D+03)*E + R13/(X-R13)
      RT2 = (((((-2.62627010965435D-04*X+3.49187925428138D-02)*X -   &
          3.09337618731880D+00)*X+1.07037141010778D+02)*X -         &
          2.36659637247087D+03)*X + ((-2.91669113681020D+06/X +     &
          1.41129505262758D+06)/X-2.91532335433779D+05)/X +         &
          3.35202872835409D+04)*E + R23/(X-R23)
      RT3 = ((((( 9.31856404738601D-05*X-2.87029400759565D-02)*X -   &
          7.83503697918455D-01)*X-1.84338896480695D+01)*X +         &
          4.04996712650414D+02)*X + (-1.89829509315154D+05/X +      &
          5.11498390849158D+04)/X-6.88145821789955D+03)*E + R33/(X-R33)
      GO TO 260
  640 RT1 = ((((-4.97561537069643D-04*X-5.00929599665316D-02)*X +        &
           1.31099142238996D+00)*X-1.88336409225481D+01)*X -             &
           6.60344754467191D+02 /X+1.64931462413877D+02)*E + R13/(X-R13)
      RT2 = ((((-4.48218898474906D-03*X-5.17373211334924D-01)*X +        &
           1.13691058739678D+01)*X-1.65426392885291D+02)*X -             &
           6.30909125686731D+03 /X+1.52231757709236D+03)*E + R23/(X-R23)
      RT3 = ((((-1.38368602394293D-02*X-1.77293428863008D+00)*X +        &
         1.73639054044562D+01)*X-3.57615122086961D+02)*X -             &
         1.45734701095912D+04 /X+2.69831813951849D+03)*E + R33/(X-R33)
      GO TO 260
!     X = 33.0 TO INFINITY             NROOTS=1,2, OR 3
  660 WW1 = SQRT(PIE4/X)
      IF (NROOTS-2) 680,700,720
  680 RT1 = 0.5D+00/(X-0.5D+00)
      UAUX=U
      WAUX=W
      RETURN
  700 IF (X .GT. 40.0D+00) GO TO 740
      E = EXP(-X)
      RT1 = (-8.78947307498880D-01*X+1.09243702330261D+01)*E + R12/(X- &
           R12)                                                        
      RT2 = (-9.28903924275977D+00*X+8.10642367843811D+01)*E + R22/(X- &
           R22)                                                        
      WW2 = ( 4.46857389308400D+00*X-7.79250653461045D+01)*E + W22*WW1 
      WW1 = WW1-WW2                                                    
      UAUX=U
      WAUX=W
      RETURN
  720 IF (X .GT. 47.0D+00) GO TO 760
      E = EXP(-X)
      RT1 = ((-7.39058467995275D+00*X+3.21318352526305D+02)*X -      &
          3.99433696473658D+03)*E + R13/(X-R13)                     
      RT2 = ((-7.38726243906513D+01*X+3.13569966333873D+03)*X -     &
          3.86862867311321D+04)*E + R23/(X-R23)                    
      RT3 = ((-2.63750565461336D+02*X+1.04412168692352D+04)*X -     &
          1.28094577915394D+05)*E + R33/(X-R33)                    
      WW3 = ((( 1.52258947224714D-01*X-8.30661900042651D+00)*X +      &
          1.92977367967984D+02)*X-1.67787926005344D+03)*E + W33*WW1  
      WW2 = (( 6.15072615497811D+01*X-2.91980647450269D+03)*X +      &
           3.80794303087338D+04)*E + W23*WW1                         
      WW1 = WW1-WW2-WW3
      UAUX=U
      WAUX=W
      RETURN
  740 RT1 = R12/(X-R12)
      RT2 = R22/(X-R22)
      WW2 = W22*WW1
      WW1 = WW1-WW2
      UAUX=U
      WAUX=W
      RETURN
  760 RT1 = R13/(X-R13)
      RT2 = R23/(X-R23)
      RT3 = R33/(X-R33)
      WW2 = W23*WW1
      WW3 = W33*WW1
      WW1 = WW1-WW2-WW3
!-----------------------------------------------------------------------  
      UAUX=U
      WAUX=W
      RETURN
      END      
! ROOT4NOF
      SUBROUTINE ROOT4NOF(X,UAUX,WAUX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION,INTENT(IN)::X
      DOUBLE PRECISION,DIMENSION(13),INTENT(OUT)::UAUX,WAUX
      DOUBLE PRECISION,DIMENSION(13)::U,W
!
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
!
      DATA R14,PIE4/1.45303521503316D-01, 7.85398163397448D-01/
      DATA R24,W24/ 1.33909728812636D+00, 2.34479815323517D-01/
      DATA R34,W34/ 3.92696350135829D+00, 1.92704402415764D-02/
      DATA R44,W44/ 8.58863568901199D+00, 2.25229076750736D-04/      
!-----------------------------------------------------------------------
      IF (X .GT. 15.0D+00) GO TO 180
      IF (X .GT. 5.0D+00) GO TO 140
      IF (X .GT. 1.0D+00) GO TO 120
      IF (X .GT. 3.0D-07) GO TO 100
!     X IS APPROXIMATELY ZERO.                   NROOTS = 4
      RT1 = 3.48198973061471D-02 -4.09645850660395D-03 *X
      RT2 = 3.81567185080042D-01 -4.48902570656719D-02 *X
      RT3 = 1.73730726945891D+00 -2.04389090547327D-01 *X
      RT4 = 1.18463056481549D+01 -1.39368301742312D+00 *X
      WW1 = 3.62683783378362D-01 -3.13844305713928D-02 *X
      WW2 = 3.13706645877886D-01 -8.98046242557724D-02 *X
      WW3 = 2.22381034453372D-01 -1.29314370958973D-01 *X
      WW4 = 1.01228536290376D-01 -8.28299075414321D-02 *X
      UAUX=U
      WAUX=W
      RETURN
!
!     X=0.0 TO 1.0                               NROOTS = 4
  100 RT1 = ((((((-1.95309614628539D-10*X+5.19765728707592D-09)*X-  &
              1.01756452250573D-07 )*X+1.72365935872131D-06 )*X-   &
              2.61203523522184D-05 )*X+3.52921308769880D-04 )*X-   &
              4.09645850658433D-03 )*X+3.48198973061469D-02
      RT2 = (((((-1.89554881382342D-08*X+3.07583114342365D-07)*X+   &
             1.270981734393D-06)*X-1.417298563884D-04)*X+          &
             3.226979163176D-03)*X-4.48902570678178D-02 )*X+       &
             3.81567185080039D-01
      RT3 = (((((( 1.77280535300416D-09*X+3.36524958870615D-08)*X-  &
             2.58341529013893D-07 )*X-1.13644895662320D-05 )*X-    &
             7.91549618884063D-05 )*X+1.03825827346828D-02 )*X-    &
             2.04389090525137D-01 )*X+1.73730726945889D+00
      RT4 = (((((-5.61188882415248D-08*X-2.49480733072460D-07)*X+   &
             3.428685057114D-06)*X+1.679007454539D-04)*X+          &
             4.722855585715D-02)*X-1.39368301737828D+00 )*X+       &
             1.18463056481543D+01
      WW1 = ((((((-1.14649303201279D-08*X+1.88015570196787D-07)*X-  &
            2.33305875372323D-06 )*X+2.68880044371597D-05 )*X-     &
            2.94268428977387D-04 )*X+3.06548909776613D-03 )*X-     &
            3.13844305680096D-02 )*X+3.62683783378335D-01
      WW2 = ((((((((-4.11720483772634D-09*X+6.54963481852134D-08)*X- &
             7.20045285129626D-07 )*X+6.93779646721723D-06 )*X-     &
             6.05367572016373D-05 )*X+4.74241566251899D-04 )*X-     &
             3.26956188125316D-03 )*X+1.91883866626681D-02 )*X-     &
             8.98046242565811D-02 )*X+3.13706645877886D-01          
      WW3 = ((((((((-3.41688436990215D-08*X+5.07238960340773D-07)*X- &
             5.01675628408220D-06 )*X+4.20363420922845D-05 )*X-     &
             3.08040221166823D-04 )*X+1.94431864731239D-03 )*X-     &
             1.02477820460278D-02 )*X+4.28670143840073D-02 )*X-     &
             1.29314370962569D-01 )*X+2.22381034453369D-01
      WW4 = ((((((((( 4.99660550769508D-09*X-7.94585963310120D-08)*X+ &
             8.359072409485D-07)*X-7.422369210610D-06)*X+            &
             5.763374308160D-05)*X-3.86645606718233D-04 )*X+         &
             2.18417516259781D-03 )*X-9.99791027771119D-03 )*X+      &
             3.48791097377370D-02 )*X-8.28299075413889D-02 )*X+      &
             1.01228536290376D-01                                    
      UAUX=U
      WAUX=W
      RETURN
!
!     X= 1.0 TO 5.0                              NROOTS = 4
  120 Y = X-3.0D+00
      RT1 = (((((((((-1.48570633747284D-15*Y-1.33273068108777D-13)*Y+ &
            4.068543696670D-12)*Y-9.163164161821D-11)*Y+               &
            2.046819017845D-09)*Y-4.03076426299031D-08 )*Y+            &
            7.29407420660149D-07 )*Y-1.23118059980833D-05 )*Y+         &
            1.88796581246938D-04 )*Y-2.53262912046853D-03 )*Y+         &
            2.51198234505021D-02                                       
      RT2 = ((((((((( 1.35830583483312D-13*Y-2.29772605964836D-12)*Y- &
            3.821500128045D-12)*Y+6.844424214735D-10)*Y-               &
            1.048063352259D-08)*Y+1.50083186233363D-08 )*Y+            &
            3.48848942324454D-06 )*Y-1.08694174399193D-04 )*Y+         &
            2.08048885251999D-03 )*Y-2.91205805373793D-02 )*Y+         &
            2.72276489515713D-01                                       
      RT3 = ((((((((( 5.02799392850289D-13*Y+1.07461812944084D-11)*Y- &
            1.482277886411D-10)*Y-2.153585661215D-09)*Y+               &
            3.654087802817D-08)*Y+5.15929575830120D-07 )*Y-            &
            9.52388379435709D-06 )*Y-2.16552440036426D-04 )*Y+         &
            9.03551469568320D-03 )*Y-1.45505469175613D-01 )*Y+         &
            1.21449092319186D+00                                       
      RT4 = (((((((((-1.08510370291979D-12*Y+6.41492397277798D-11)*Y+ &
            7.542387436125D-10)*Y-2.213111836647D-09)*Y-               &
            1.448228963549D-07)*Y-1.95670833237101D-06 )*Y-            &
            1.07481314670844D-05 )*Y+1.49335941252765D-04 )*Y+         &
            4.87791531990593D-02 )*Y-1.10559909038653D+00 )*Y+         &
            8.09502028611780D+00                                       
      WW1 = ((((((((((-4.65801912689961D-14*Y+7.58669507106800D-13)*Y- &
            1.186387548048D-11)*Y+1.862334710665D-10)*Y-                &
            2.799399389539D-09)*Y+4.148972684255D-08)*Y-                &
            5.933568079600D-07)*Y+8.168349266115D-06)*Y-                &
            1.08989176177409D-04 )*Y+1.41357961729531D-03 )*Y-          &
            1.87588361833659D-02 )*Y+2.89898651436026D-01               
      WW2 = ((((((((((((-1.46345073267549D-14*Y+2.25644205432182D-13)*Y- &
            3.116258693847D-12)*Y+4.321908756610D-11)*Y-                  &
            5.673270062669D-10)*Y+7.006295962960D-09)*Y-                  &
            8.120186517000D-08)*Y+8.775294645770D-07)*Y-                  &
            8.77829235749024D-06 )*Y+8.04372147732379D-05 )*Y-            &
            6.64149238804153D-04 )*Y+4.81181506827225D-03 )*Y-            &
            2.88982669486183D-02 )*Y+1.56247249979288D-01
      WW3 = ((((((((((((( 9.06812118895365D-15*Y-1.40541322766087D-13)*  &
            Y+1.919270015269D-12)*Y-2.605135739010D-11)*Y+                &
            3.299685839012D-10)*Y-3.86354139348735D-09 )*Y+               &
            4.16265847927498D-08 )*Y-4.09462835471470D-07 )*Y+            &
            3.64018881086111D-06 )*Y-2.88665153269386D-05 )*Y+            &
            2.00515819789028D-04 )*Y-1.18791896897934D-03 )*Y+            &
            5.75223633388589D-03 )*Y-2.09400418772687D-02 )*Y+            &
            4.85368861938873D-02                                          
      WW4 = ((((((((((((((-9.74835552342257D-16*Y+1.57857099317175D-14)* &
            Y-2.249993780112D-13)*Y+3.173422008953D-12)*Y-                &
            4.161159459680D-11)*Y+5.021343560166D-10)*Y-                  &
            5.545047534808D-09)*Y+5.554146993491D-08)*Y-                  &
            4.99048696190133D-07 )*Y+3.96650392371311D-06 )*Y-            &
            2.73816413291214D-05 )*Y+1.60106988333186D-04 )*Y-            &
            7.64560567879592D-04 )*Y+2.81330044426892D-03 )*Y-            &
            7.16227030134947D-03 )*Y+9.66077262223353D-03                 
      UAUX=U
      WAUX=W
      RETURN                                                             
!
  140 IF (X .GT. 10.0D+00) GO TO 160
!     X=5.0 TO 10.0                              NROOTS = 4
      Y = X-7.5D+00
      RT1 = ((((((((( 4.64217329776215D-15*Y-6.27892383644164D-15)*Y+     &
            3.462236347446D-13)*Y-2.927229355350D-11)*Y+                 &
            5.090355371676D-10)*Y-9.97272656345253D-09 )*Y+              &
            2.37835295639281D-07 )*Y-4.60301761310921D-06 )*Y+           &
            8.42824204233222D-05 )*Y-1.37983082233081D-03 )*Y+          &
            1.66630865869375D-02                                        
      RT2 = ((((((((( 2.93981127919047D-14*Y+8.47635639065744D-13)*Y-    &
            1.446314544774D-11)*Y-6.149155555753D-12)*Y+                &
            8.484275604612D-10)*Y-6.10898827887652D-08 )*Y+             &
            2.39156093611106D-06 )*Y-5.35837089462592D-05 )*Y+          &
            1.00967602595557D-03 )*Y-1.57769317127372D-02 )*Y+          &
            1.74853819464285D-01                                        
      RT3 = (((((((((( 2.93523563363000D-14*Y-6.40041776667020D-14)*Y-   &
            2.695740446312D-12)*Y+1.027082960169D-10)*Y-                &
            5.822038656780D-10)*Y-3.159991002539D-08)*Y+                &
            4.327249251331D-07)*Y+4.856768455119D-06)*Y-                &
            2.54617989427762D-04 )*Y+5.54843378106589D-03 )*Y-          &
            7.95013029486684D-02 )*Y+7.20206142703162D-01               
      RT4 = (((((((((((-1.62212382394553D-14*Y+7.68943641360593D-13)*Y+  &
            5.764015756615D-12)*Y-1.380635298784D-10)*Y-                &
            1.476849808675D-09)*Y+1.84347052385605D-08 )*Y+             &
            3.34382940759405D-07 )*Y-1.39428366421645D-06 )*Y-          &
            7.50249313713996D-05 )*Y-6.26495899187507D-04 )*Y+          &
            4.69716410901162D-02 )*Y-6.66871297428209D-01 )*Y+          &
            4.11207530217806D+00
      WW1 = ((((((((((-1.65995045235997D-15*Y+6.91838935879598D-14)*Y-   &
            9.131223418888D-13)*Y+1.403341829454D-11)*Y-                &
            3.672235069444D-10)*Y+6.366962546990D-09)*Y-                &
            1.039220021671D-07)*Y+1.959098751715D-06)*Y-                &
            3.33474893152939D-05 )*Y+5.72164211151013D-04 )*Y-          &
            1.05583210553392D-02 )*Y+2.26696066029591D-01               
      WW2 = ((((((((((((-3.57248951192047D-16*Y+6.25708409149331D-15)*Y- &
            9.657033089714D-14)*Y+1.507864898748D-12)*Y-                &
            2.332522256110D-11)*Y+3.428545616603D-10)*Y-                &
            4.698730937661D-09)*Y+6.219977635130D-08)*Y-                &
            7.83008889613661D-07 )*Y+9.08621687041567D-06 )*Y-          &
            9.86368311253873D-05 )*Y+9.69632496710088D-04 )*Y-          &
            8.14594214284187D-03 )*Y+8.50218447733457D-02
      WW3 = ((((((((((((( 1.64742458534277D-16*Y-2.68512265928410D-15)*  &
            Y+3.788890667676D-14)*Y-5.508918529823D-13)*Y+              &
            7.555896810069D-12)*Y-9.69039768312637D-11 )*Y+             &
            1.16034263529672D-09 )*Y-1.28771698573873D-08 )*Y+          &
            1.31949431805798D-07 )*Y-1.23673915616005D-06 )*Y+          &
            1.04189803544936D-05 )*Y-7.79566003744742D-05 )*Y+          &
            5.03162624754434D-04 )*Y-2.55138844587555D-03 )*Y+          &
            1.13250730954014D-02                                        
      WW4 = ((((((((((((((-1.55714130075679D-17*Y+2.57193722698891D-16)* &
            Y-3.626606654097D-15)*Y+5.234734676175D-14)*Y-              &
            7.067105402134D-13)*Y+8.793512664890D-12)*Y-                &
            1.006088923498D-10)*Y+1.050565098393D-09)*Y-                &
            9.91517881772662D-09 )*Y+8.35835975882941D-08 )*Y-          &
            6.19785782240693D-07 )*Y+3.95841149373135D-06 )*Y-          &
            2.11366761402403D-05 )*Y+9.00474771229507D-05 )*Y-          &
            2.78777909813289D-04 )*Y+5.26543779837487D-04               
      UAUX=U
      WAUX=W
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
      WW3 = ((((((((((((-4.19569145459480D-17*Y+5.94344180261644D-16)*Y- &
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
      UAUX=U
      WAUX=W
      RETURN                                                            
!
  180 WW1 = SQRT(PIE4/X)
      IF (X .GT. 35.0D+00) GO TO 220
      IF (X .GT. 20.0D+00) GO TO 200
!     X=15.0 TO 20.0                             NROOTS = 4
      Y = X-17.5D+00
      RT1 = ((((((((((( 4.36701759531398D-17*Y-1.12860600219889D-16)*Y-   &
            6.149849164164D-15)*Y+5.820231579541D-14)*Y+                 &
            4.396602872143D-13)*Y-1.24330365320172D-11 )*Y+              &
            6.71083474044549D-11 )*Y+2.43865205376067D-10 )*Y+           &
            1.67559587099969D-08 )*Y-9.32738632357572D-07 )*Y+           &
            2.39030487004977D-05 )*Y-4.68648206591515D-04 )*Y+           &
            8.34977776583956D-03                                         
      RT2 = ((((((((((( 4.98913142288158D-16*Y-2.60732537093612D-16)*Y-   &
            7.775156445127D-14)*Y+5.766105220086D-13)*Y+                 &
            6.432696729600D-12)*Y-1.39571683725792D-10 )*Y+              &
            5.95451479522191D-10 )*Y+2.42471442836205D-09 )*Y+           &
            2.47485710143120D-07 )*Y-1.14710398652091D-05 )*Y+           &
            2.71252453754519D-04 )*Y-4.96812745851408D-03 )*Y+           &
            8.26020602026780D-02                                         
      RT3 = ((((((((((( 1.91498302509009D-15*Y+1.48840394311115D-14)*Y-   &
           4.316925145767D-13)*Y+1.186495793471D-12)*Y+                  &
           4.615806713055D-11)*Y-5.54336148667141D-10 )*Y+               &
           3.48789978951367D-10 )*Y-2.79188977451042D-09 )*Y+            &
           2.09563208958551D-06 )*Y-6.76512715080324D-05 )*Y+            &
           1.32129867629062D-03 )*Y-2.05062147771513D-02 )*Y+            &
           2.88068671894324D-01
      RT4 = (((((((((((-5.43697691672942D-15*Y-1.12483395714468D-13)*Y+   &
            2.826607936174D-12)*Y-1.266734493280D-11)*Y-                 &
            4.258722866437D-10)*Y+9.45486578503261D-09 )*Y-              &
            5.86635622821309D-08 )*Y-1.28835028104639D-06 )*Y+           &
            4.41413815691885D-05 )*Y-7.61738385590776D-04 )*Y+           &
            9.66090902985550D-03 )*Y-1.01410568057649D-01 )*Y+           &
            9.54714798156712D-01                                         
      WW4 = ((((((((((((-7.56882223582704D-19*Y+7.53541779268175D-18)*Y-  &
            1.157318032236D-16)*Y+2.411195002314D-15)*Y-                 &
            3.601794386996D-14)*Y+4.082150659615D-13)*Y-                 &
            4.289542980767D-12)*Y+5.086829642731D-11)*Y-                 &
            6.35435561050807D-10 )*Y+6.82309323251123D-09 )*Y-           &
            5.63374555753167D-08 )*Y+3.57005361100431D-07 )*Y-           &
            2.40050045173721D-06 )*Y+4.94171300536397D-05                
      WW3 = (((((((((((-5.54451040921657D-17*Y+2.68748367250999D-16)*Y+   &
            1.349020069254D-14)*Y-2.507452792892D-13)*Y+                 &
            1.944339743818D-12)*Y-1.29816917658823D-11 )*Y+              &
            3.49977768819641D-10 )*Y-8.67270669346398D-09 )*Y+           &
            1.31381116840118D-07 )*Y-1.36790720600822D-06 )*Y+           &
            1.19210697673160D-05 )*Y-1.42181943986587D-04 )*Y+           &
            4.12615396191829D-03
      WW2 = (((((((((((-1.86506057729700D-16*Y+1.16661114435809D-15)*Y+  &
            2.563712856363D-14)*Y-4.498350984631D-13)*Y+                &
            1.765194089338D-12)*Y+9.04483676345625D-12 )*Y+             &
            4.98930345609785D-10 )*Y-2.11964170928181D-08 )*Y+          &
            3.98295476005614D-07 )*Y-5.49390160829409D-06 )*Y+          &
            7.74065155353262D-05 )*Y-1.48201933009105D-03 )*Y+          &
            4.97836392625268D-02                                        
      WW1 = (( 1.9623264149430D-01/X-4.9695241464490D-01)/X -            &
            6.0156581186481D-05)*EXP(-X)+WW1-WW2-WW3-WW4
      UAUX=U
      WAUX=W
      RETURN
!
!     X=20.0 TO 35.0                             NROOTS = 4
  200 E = EXP(-X)
      RT1 = ((((((-4.45711399441838D-05*X+1.27267770241379D-03)*X -   &
            2.36954961381262D-01)*X+1.54330657903756D+01)*X -          &
            5.22799159267808D+02)*X+1.05951216669313D+04)*X + (-       &
            2.51177235556236D+06/X+8.72975373557709D+05)/X -           &
            1.29194382386499D+05)*E + R14/(X-R14)
      RT2 = (((((-7.85617372254488D-02*X+6.35653573484868D+00)*X -    &
            3.38296938763990D+02)*X+1.25120495802096D+04)*X -          &
            3.16847570511637D+05)*X + ((-1.02427466127427D+09/X +      &
            3.70104713293016D+08)/X-5.87119005093822D+07)/X +          &
            5.38614211391604D+06)*E + R24/(X-R24)
      RT3 = (((((-2.37900485051067D-01*X+1.84122184400896D+01)*X -    &
             1.00200731304146D+03)*X+3.75151841595736D+04)*X -          &
             9.50626663390130D+05)*X + ((-2.88139014651985D+09/X +      &
             1.06625915044526D+09)/X-1.72465289687396D+08)/X +          &
             1.60419390230055D+07)*E + R34/(X-R34)
      RT4 = ((((((-6.00691586407385D-04*X-3.64479545338439D-01)*X +   &
            1.57496131755179D+01)*X-6.54944248734901D+02)*X +          &
            1.70830039597097D+04)*X-2.90517939780207D+05)*X + (+       &
            3.49059698304732D+07/X-1.64944522586065D+07)/X +           &
            2.96817940164703D+06)*E + R44/(X-R44)
      IF (X .LE. 25.0D+00) WW4 = ((((((( 2.33766206773151D-07*X-      &
           3.81542906607063D-05)*X +3.51416601267000D-03)*X-          &
           1.66538571864728D-01)*X +4.80006136831847D+00)*X-          &
           8.73165934223603D+01)*X +9.77683627474638D+02)*X +         &
           1.66000945117640D+04/X -6.14479071209961D+03)*E + W44*WW1
      IF (X .GT. 25.0D+00) WW4 = (((((( 5.74245945342286D-06*X-       &
            7.58735928102351D-05)*X +2.35072857922892D-04)*X-          &
            3.78812134013125D-03)*X +3.09871652785805D-01)*X-          &
            7.11108633061306D+00)*X +5.55297573149528D+01)*E + W44*WW1 
      WW3 = (((((( 2.36392855180768D-04*X-9.16785337967013D-03)*X +   &
            4.62186525041313D-01)*X-1.96943786006540D+01)*X +          &
            4.99169195295559D+02)*X-6.21419845845090D+03)*X + ((+      &
            5.21445053212414D+07/X-1.34113464389309D+07)/X +           &
            1.13673298305631D+06)/X-2.81501182042707D+03)*E + W34*WW1
      WW2 = (((((( 7.29841848989391D-04*X-3.53899555749875D-02)*X +   &
            2.07797425718513D+00)*X-1.00464709786287D+02)*X +          &
            3.15206108877819D+03)*X-6.27054715090012D+04)*X + (+       &
            1.54721246264919D+07/X-5.26074391316381D+06)/X +           &
            7.67135400969617D+05)*E + W24*WW1
      WW1 = (( 1.9623264149430D-01/X-4.9695241464490D-01)/X -         &
            6.0156581186481D-05)*E + WW1-WW2-WW3-WW4                   
      UAUX=U
      WAUX=W
      RETURN                                                          
!                                                                     
  220 IF (X .GT. 53.0D+00) GO TO 240
!     X=35.0 TO 53.0                             NROOTS = 4
      E = EXP(-X)*(X*X)**2
      RT4 = ((-2.19135070169653D-03*X-1.19108256987623D-01)*X -  &
           7.50238795695573D-01)*E + R44/(X-R44)
      RT3 = ((-9.65842534508637D-04*X-4.49822013469279D-02)*X +  &
           6.08784033347757D-01)*E + R34/(X-R34)
      RT2 = ((-3.62569791162153D-04*X-9.09231717268466D-03)*X +   &
           1.84336760556262D-01)*E + R24/(X-R24)
      RT1 = ((-4.07557525914600D-05*X-6.88846864931685D-04)*X +    &
           1.74725309199384D-02)*E + R14/(X-R14)
      WW4 = (( 5.76631982000990D-06*X-7.89187283804890D-05)*X +     &
           3.28297971853126D-04)*E + W44*WW1
      WW3 = (( 2.08294969857230D-04*X-3.77489954837361D-03)*X +  &
           2.09857151617436D-02)*E + W34*WW1
      WW2 = (( 6.16374517326469D-04*X-1.26711744680092D-02)*X +   &
           8.14504890732155D-02)*E + W24*WW1
      WW1 = WW1-WW2-WW3-WW4
      UAUX=U
      WAUX=W
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
!-----------------------------------------------------------------------  
      UAUX=U
      WAUX=W
      RETURN
      END            
! ROOT5NOF
      SUBROUTINE ROOT5NOF(X,UAUX,WAUX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION,INTENT(IN)::X
      DOUBLE PRECISION,DIMENSION(13),INTENT(OUT)::UAUX,WAUX
      DOUBLE PRECISION,DIMENSION(13)::U,W
!
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
!
      DATA R15,PIE4/1.17581320211778D-01, 7.85398163397448D-01/
      DATA R25,W25/ 1.07456201243690D+00, 2.70967405960535D-01/
      DATA R35,W35/ 3.08593744371754D+00, 3.82231610015404D-02/
      DATA R45,W45/ 6.41472973366203D+00, 1.51614186862443D-03/
      DATA R55,W55/ 1.18071894899717D+01, 8.62130526143657D-06/
!-----------------------------------------------------------------------
      IF (X .GT. 15.0D+00) GO TO 180
      IF (X .GT. 5.0D+00) GO TO 140
      IF (X .GT. 1.0D+00) GO TO 120
      IF (X .GT. 3.0D-07) GO TO 100
!     X IS APPROXIMATELY ZERO.                   NROOTS = 5
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
      UAUX=U
      WAUX=W
      RETURN
!
!     X=0.0 TO 1.0                               NROOTS = 5
  100 RT1 = ((((((-4.46679165328413D-11*X+1.21879111988031D-09)*X-     &
           2.62975022612104D-08 )*X+5.15106194905897D-07 )*X-          &
           9.27933625824749D-06 )*X+1.51794097682482D-04 )*X-          &
           2.15865967920301D-03 )*X+2.26659266316985D-02
      RT2 = (((((( 1.93117331714174D-10*X-4.57267589660699D-09)*X+     &
           2.48339908218932D-08 )*X+1.50716729438474D-06 )*X-          &
           6.07268757707381D-05 )*X+1.37506939145643D-03 )*X-          &
           2.20258754419939D-02 )*X+2.31271692140905D-01
      RT3 = ((((( 4.84989776180094D-09*X+1.31538893944284D-07)*X-      &
           2.766753852879D-06)*X-7.651163510626D-05)*X+                &
           4.033058545972D-03)*X-8.16520022916145D-02 )*X+             &
           8.57346024118779D-01
      RT4 = ((((-2.48581772214623D-07*X-4.34482635782585D-06)*X-       &
           7.46018257987630D-07 )*X+1.01210776517279D-02 )*X-          &
           2.83193369640005D-01 )*X+2.97353038120345D+00           
      RT5 = (((((-8.92432153868554D-09*X+1.77288899268988D-08)*X+      &
           3.040754680666D-06)*X+1.058229325071D-04)*X+                &
           4.596379534985D-02)*X-1.75382723579114D+00 )*X+             &
           1.84151859759049D+01                                        
      WW1 = ((((((-2.03822632771791D-09*X+3.89110229133810D-08)*X-     &
            5.84914787904823D-07 )*X+8.30316168666696D-06 )*X-         &
            1.13218402310546D-04 )*X+1.49128888586790D-03 )*X-         &
            1.96867576904816D-02 )*X+2.95524224714749D-01
      WW2 = ((((((( 8.62848118397570D-09*X-1.38975551148989D-07)*X+    &
           1.602894068228D-06)*X-1.646364300836D-05)*X+                &
           1.538445806778D-04)*X-1.28848868034502D-03 )*X+             &
           9.38866933338584D-03 )*X-5.61737590178812D-02 )*X+          &
           2.69266719309991D-01                                    
      WW3 = ((((((((-9.41953204205665D-09*X+1.47452251067755D-07)*X-   &
           1.57456991199322D-06 )*X+1.45098401798393D-05 )*X-          &
           1.18858834181513D-04 )*X+8.53697675984210D-04 )*X-          &
           5.22877807397165D-03 )*X+2.60854524809786D-02 )*X-          &
           9.71152726809059D-02 )*X+2.19086362515979D-01
      WW4 = ((((((((-3.84961617022042D-08*X+5.66595396544470D-07)*X-   &
           5.52351805403748D-06 )*X+4.53160377546073D-05 )*X-          &
           3.22542784865557D-04 )*X+1.95682017370967D-03 )*X-          &
           9.77232537679229D-03 )*X+3.79455945268632D-02 )*X-          &
           1.02979262192227D-01 )*X+1.49451349150573D-01            
      WW5 = ((((((((( 4.09594812521430D-09*X-6.47097874264417D-08)*X+  &
           6.743541482689D-07)*X-5.917993920224D-06)*X+                &
           4.531969237381D-05)*X-2.99102856679638D-04 )*X+             &
           1.65695765202643D-03 )*X-7.40671222520653D-03 )*X+          &
           2.50889946832192D-02 )*X-5.73782817487958D-02 )*X+          &
           6.66713443086877D-02                                      
      UAUX=U
      WAUX=W
      RETURN                                                          
!
!     X=1.0 TO 5.0                               NROOTS = 5
  120 Y = X-3.0D+00
      RT1 = ((((((((-2.58163897135138D-14*Y+8.14127461488273D-13)*Y-   &
           2.11414838976129D-11 )*Y+5.09822003260014D-10 )*Y-          &
           1.16002134438663D-08 )*Y+2.46810694414540D-07 )*Y-          &
           4.92556826124502D-06 )*Y+9.02580687971053D-05 )*Y-          &
           1.45190025120726D-03 )*Y+1.73416786387475D-02               
      RT2 = ((((((((( 1.04525287289788D-14*Y+5.44611782010773D-14)*Y-  &
           4.831059411392D-12)*Y+1.136643908832D-10)*Y-                &
           1.104373076913D-09)*Y-2.35346740649916D-08 )*Y+             &
           1.43772622028764D-06 )*Y-4.23405023015273D-05 )*Y+          &
           9.12034574793379D-04 )*Y-1.52479441718739D-02 )*Y+          &
           1.76055265928744D-01                                        
      RT3 = (((((((((-6.89693150857911D-14*Y+5.92064260918861D-13)*Y+  &
           1.847170956043D-11)*Y-3.390752744265D-10)*Y-                &
           2.995532064116D-09)*Y+1.57456141058535D-07 )*Y-             &
           3.95859409711346D-07 )*Y-9.58924580919747D-05 )*Y+          &
           3.23551502557785D-03 )*Y-5.97587007636479D-02 )*Y+          &
           6.46432853383057D-01                                        
      RT4 = ((((((((-3.61293809667763D-12*Y-2.70803518291085D-11)*Y+   &
           8.83758848468769D-10 )*Y+1.59166632851267D-08 )*Y-          &
           1.32581997983422D-07 )*Y-7.60223407443995D-06 )*Y-          &
           7.41019244900952D-05 )*Y+9.81432631743423D-03 )*Y-          &
           2.23055570487771D-01 )*Y+2.21460798080643D+00               
      RT5 = ((((((((( 7.12332088345321D-13*Y+3.16578501501894D-12)*Y-  &
           8.776668218053D-11)*Y-2.342817613343D-09)*Y-                &
           3.496962018025D-08)*Y-3.03172870136802D-07 )*Y+             &
           1.50511293969805D-06 )*Y+1.37704919387696D-04 )*Y+          &
           4.70723869619745D-02 )*Y-1.47486623003693D+00 )*Y+          &
           1.35704792175847D+01                                       
      WW1 = ((((((((( 1.04348658616398D-13*Y-1.94147461891055D-12)*Y+  &
           3.485512360993D-11)*Y-6.277497362235D-10)*Y+                &
           1.100758247388D-08)*Y-1.88329804969573D-07 )*Y+             &
           3.12338120839468D-06 )*Y-5.04404167403568D-05 )*Y+          &
           8.00338056610995D-04 )*Y-1.30892406559521D-02 )*Y+          &
           2.47383140241103D-01                                        
      WW2 = ((((((((((( 3.23496149760478D-14*Y-5.24314473469311D-13)*Y+ &
           7.743219385056D-12)*Y-1.146022750992D-10)*Y+                 &
           1.615238462197D-09)*Y-2.15479017572233D-08 )*Y+              &
           2.70933462557631D-07 )*Y-3.18750295288531D-06 )*Y+           &
           3.47425221210099D-05 )*Y-3.45558237388223D-04 )*Y+           &
           3.05779768191621D-03 )*Y-2.29118251223003D-02 )*Y+           &
           1.59834227924213D-01
      WW3 = ((((((((((((-3.42790561802876D-14*Y+5.26475736681542D-13)*Y- &
           7.184330797139D-12)*Y+9.763932908544D-11)*Y-                  &
           1.244014559219D-09)*Y+1.472744068942D-08)*Y-                  &
           1.611749975234D-07)*Y+1.616487851917D-06)*Y-                  &
           1.46852359124154D-05 )*Y+1.18900349101069D-04 )*Y-            &
           8.37562373221756D-04 )*Y+4.93752683045845D-03 )*Y-            &
           2.25514728915673D-02 )*Y+6.95211812453929D-02
      WW4 = ((((((((((((( 1.04072340345039D-14*Y-1.60808044529211D-13)*  &
           Y+2.183534866798D-12)*Y-2.939403008391D-11)*Y+                &
           3.679254029085D-10)*Y-4.23775673047899D-09 )*Y+               &
           4.46559231067006D-08 )*Y-4.26488836563267D-07 )*Y+            &
           3.64721335274973D-06 )*Y-2.74868382777722D-05 )*Y+            &
           1.78586118867488D-04 )*Y-9.68428981886534D-04 )*Y+            &
           4.16002324339929D-03 )*Y-1.28290192663141D-02 )*Y+            &
           2.22353727685016D-02                                          
      WW5 = ((((((((((((((-8.16770412525963D-16*Y+1.31376515047977D-14)* &
           Y-1.856950818865D-13)*Y+2.596836515749D-12)*Y-                &
           3.372639523006D-11)*Y+4.025371849467D-10)*Y-                  &
           4.389453269417D-09)*Y+4.332753856271D-08)*Y-                  &
           3.82673275931962D-07 )*Y+2.98006900751543D-06 )*Y-            &
           2.00718990300052D-05 )*Y+1.13876001386361D-04 )*Y-            &
           5.23627942443563D-04 )*Y+1.83524565118203D-03 )*Y-            &
           4.37785737450783D-03 )*Y+5.36963805223095D-03
      UAUX=U
      WAUX=W
      RETURN
!
  140 IF (X .GT. 10.0D+00) GO TO 160
!     X=5.0 TO 10.0                              NROOTS = 5
      Y = X-7.5D+00
      RT1 = ((((((((-1.13825201010775D-14*Y+1.89737681670375D-13)*Y-   &
           4.81561201185876D-12 )*Y+1.56666512163407D-10 )*Y-          &
           3.73782213255083D-09 )*Y+9.15858355075147D-08 )*Y-          &
           2.13775073585629D-06 )*Y+4.56547356365536D-05 )*Y-          &
           8.68003909323740D-04 )*Y+1.22703754069176D-02             
      RT2 = (((((((((-3.67160504428358D-15*Y+1.27876280158297D-14)*Y-  &
           1.296476623788D-12)*Y+1.477175434354D-11)*Y+                &
           5.464102147892D-10)*Y-2.42538340602723D-08 )*Y+             &
           8.20460740637617D-07 )*Y-2.20379304598661D-05 )*Y+          &
           4.90295372978785D-04 )*Y-9.14294111576119D-03 )*Y+          &
           1.22590403403690D-01                                        
      RT3 = ((((((((( 1.39017367502123D-14*Y-6.96391385426890D-13)*Y+  &
           1.176946020731D-12)*Y+1.725627235645D-10)*Y-                &
           3.686383856300D-09)*Y+2.87495324207095D-08 )*Y+             &
           1.71307311000282D-06 )*Y-7.94273603184629D-05 )*Y+          &
           2.00938064965897D-03 )*Y-3.63329491677178D-02 )*Y+          &
           4.34393683888443D-01                                        
      RT4 = ((((((((((-1.27815158195209D-14*Y+1.99910415869821D-14)*Y+ &
           3.753542914426D-12)*Y-2.708018219579D-11)*Y-                &
           1.190574776587D-09)*Y+1.106696436509D-08)*Y+                &
           3.954955671326D-07)*Y-4.398596059588D-06)*Y-                &
           2.01087998907735D-04 )*Y+7.89092425542937D-03 )*Y-          &
           1.42056749162695D-01 )*Y+1.39964149420683D+00              
      RT5 = ((((((((((-1.19442341030461D-13*Y-2.34074833275956D-12)*Y+ &
           6.861649627426D-12)*Y+6.082671496226D-10)*Y+                &
           5.381160105420D-09)*Y-6.253297138700D-08)*Y-                &
           2.135966835050D-06)*Y-2.373394341886D-05)*Y+                &
           2.88711171412814D-06 )*Y+4.85221195290753D-02 )*Y-          &
           1.04346091985269D+00 )*Y+7.89901551676692D+00              
      WW1 = ((((((((( 7.95526040108997D-15*Y-2.48593096128045D-13)*Y+  &
           4.761246208720D-12)*Y-9.535763686605D-11)*Y+                &
           2.225273630974D-09)*Y-4.49796778054865D-08 )*Y+             &
           9.17812870287386D-07 )*Y-1.86764236490502D-05 )*Y+          &
           3.76807779068053D-04 )*Y-8.10456360143408D-03 )*Y+          &
           2.01097936411496D-01                                        
      WW2 = ((((((((((( 1.25678686624734D-15*Y-2.34266248891173D-14)*Y+ &
           3.973252415832D-13)*Y-6.830539401049D-12)*Y+                 &
           1.140771033372D-10)*Y-1.82546185762009D-09 )*Y+              &
           2.77209637550134D-08 )*Y-4.01726946190383D-07 )*Y+           &
           5.48227244014763D-06 )*Y-6.95676245982121D-05 )*Y+           &
           8.05193921815776D-04 )*Y-8.15528438784469D-03 )*Y+           &
           9.71769901268114D-02                                        
      WW3 = ((((((((((((-8.20929494859896D-16*Y+1.37356038393016D-14)*Y- &
           2.022863065220D-13)*Y+3.058055403795D-12)*Y-                  &
           4.387890955243D-11)*Y+5.923946274445D-10)*Y-                  &
           7.503659964159D-09)*Y+8.851599803902D-08)*Y-                  &
           9.65561998415038D-07 )*Y+9.60884622778092D-06 )*Y-            &
           8.56551787594404D-05 )*Y+6.66057194311179D-04 )*Y-            &
           4.17753183902198D-03 )*Y+2.25443826852447D-02
      WW4 = ((((((((((((((-1.08764612488790D-17*Y+1.85299909689937D-16)* &
           Y-2.730195628655D-15)*Y+4.127368817265D-14)*Y-                &
           5.881379088074D-13)*Y+7.805245193391D-12)*Y-                  &
           9.632707991704D-11)*Y+1.099047050624D-09)*Y-                  &
           1.15042731790748D-08 )*Y+1.09415155268932D-07 )*Y-            &
           9.33687124875935D-07 )*Y+7.02338477986218D-06 )*Y-            &
           4.53759748787756D-05 )*Y+2.41722511389146D-04 )*Y-            &
           9.75935943447037D-04 )*Y+2.57520532789644D-03                 
      WW5 = ((((((((((((((( 7.28996979748849D-19*Y-1.26518146195173D-17) &
            *Y+1.886145834486D-16)*Y-2.876728287383D-15)*Y+              &
           4.114588668138D-14)*Y-5.44436631413933D-13 )*Y+               &
           6.64976446790959D-12 )*Y-7.44560069974940D-11 )*Y+            &
           7.57553198166848D-10 )*Y-6.92956101109829D-09 )*Y+            &
           5.62222859033624D-08 )*Y-3.97500114084351D-07 )*Y+            &
           2.39039126138140D-06 )*Y-1.18023950002105D-05 )*Y+            &
           4.52254031046244D-05 )*Y-1.21113782150370D-04 )*Y+            &
           1.75013126731224D-04                                        
      UAUX=U
      WAUX=W
      RETURN                                                           
!
!     X=10.0 TO 15.0                             NROOTS = 5
  160 Y = X-12.5D+00
      RT1 = ((((((((((-4.16387977337393D-17*Y+7.20872997373860D-16)*Y+ &
           1.395993802064D-14)*Y+3.660484641252D-14)*Y-               &
           4.154857548139D-12)*Y+2.301379846544D-11)*Y-               &
           1.033307012866D-09)*Y+3.997777641049D-08)*Y-               &
           9.35118186333939D-07 )*Y+2.38589932752937D-05 )*Y-         &
           5.35185183652937D-04 )*Y+8.85218988709735D-03              
      RT2 = ((((((((((-4.56279214732217D-16*Y+6.24941647247927D-15)*Y+ &
           1.737896339191D-13)*Y+8.964205979517D-14)*Y-               &
           3.538906780633D-11)*Y+9.561341254948D-11)*Y-               &
           9.772831891310D-09)*Y+4.240340194620D-07)*Y-               &
           1.02384302866534D-05 )*Y+2.57987709704822D-04 )*Y-         &
           5.54735977651677D-03 )*Y+8.68245143991948D-02              
      RT3 = ((((((((((-2.52879337929239D-15*Y+2.13925810087833D-14)*Y+ &
           7.884307667104D-13)*Y-9.023398159510D-13)*Y-               &
           5.814101544957D-11)*Y-1.333480437968D-09)*Y-               &
           2.217064940373D-08)*Y+1.643290788086D-06)*Y-               &
           4.39602147345028D-05 )*Y+1.08648982748911D-03 )*Y-         &
           2.13014521653498D-02 )*Y+2.94150684465425D-01              
      RT4 = ((((((((((-6.42391438038888D-15*Y+5.37848223438815D-15)*Y+ &
           8.960828117859D-13)*Y+5.214153461337D-11)*Y-               &
           1.106601744067D-10)*Y-2.007890743962D-08)*Y+               &
           1.543764346501D-07)*Y+4.520749076914D-06)*Y-               &
           1.88893338587047D-04 )*Y+4.73264487389288D-03 )*Y-         &
           7.91197893350253D-02 )*Y+8.60057928514554D-01              
      RT5 = (((((((((((-2.24366166957225D-14*Y+4.87224967526081D-14)*Y+ &
           5.587369053655D-12)*Y-3.045253104617D-12)*Y-                &
           1.223983883080D-09)*Y-2.05603889396319D-09 )*Y+             &
           2.58604071603561D-07 )*Y+1.34240904266268D-06 )*Y-          &
           5.72877569731162D-05 )*Y-9.56275105032191D-04 )*Y+          &
           4.23367010370921D-02 )*Y-5.76800927133412D-01 )*Y+          &
           3.87328263873381D+00
      WW1 = ((((((((( 8.98007931950169D-15*Y+7.25673623859497D-14)*Y+  &
           5.851494250405D-14)*Y-4.234204823846D-11)*Y+               &
           3.911507312679D-10)*Y-9.65094802088511D-09 )*Y+            &
           3.42197444235714D-07 )*Y-7.51821178144509D-06 )*Y+         &
           1.94218051498662D-04 )*Y-5.38533819142287D-03 )*Y+         &
           1.68122596736809D-01                                       
      WW2 = ((((((((((-1.05490525395105D-15*Y+1.96855386549388D-14)*Y- &
           5.500330153548D-13)*Y+1.003849567976D-11)*Y-               &
           1.720997242621D-10)*Y+3.533277061402D-09)*Y-               &
           6.389171736029D-08)*Y+1.046236652393D-06)*Y-               &
           1.73148206795827D-05 )*Y+2.57820531617185D-04 )*Y-         &
           3.46188265338350D-03 )*Y+7.03302497508176D-02              
      WW3 = ((((((((((( 3.60020423754545D-16*Y-6.24245825017148D-15)*Y+ &
           9.945311467434D-14)*Y-1.749051512721D-12)*Y+                 &
           2.768503957853D-11)*Y-4.08688551136506D-10 )*Y+              &
           6.04189063303610D-09 )*Y-8.23540111024147D-08 )*Y+          &
           1.01503783870262D-06 )*Y-1.20490761741576D-05 )*Y+          &
           1.26928442448148D-04 )*Y-1.05539461930597D-03 )*Y+          &
           1.15543698537013D-02
      WW4 = ((((((((((((( 2.51163533058925D-18*Y-4.31723745510697D-17)* &
           Y+6.557620865832D-16)*Y-1.016528519495D-14)*Y+               &
           1.491302084832D-13)*Y-2.06638666222265D-12 )*Y+              &
           2.67958697789258D-11 )*Y-3.23322654638336D-10 )*Y+          &
           3.63722952167779D-09 )*Y-3.75484943783021D-08 )*Y+          &
           3.49164261987184D-07 )*Y-2.92658670674908D-06 )*Y+          &
           2.12937256719543D-05 )*Y-1.19434130620929D-04 )*Y+           &
           6.45524336158384D-04                                         
      WW5 = ((((((((((((((-1.29043630202811D-19*Y+2.16234952241296D-18)* &
           Y-3.107631557965D-17)*Y+4.570804313173D-16)*Y-               &
           6.301348858104D-15)*Y+8.031304476153D-14)*Y-                 &
           9.446196472547D-13)*Y+1.018245804339D-11)*Y-                 &
           9.96995451348129D-11 )*Y+8.77489010276305D-10 )*Y-          &
           6.84655877575364D-09 )*Y+4.64460857084983D-08 )*Y-          &
           2.66924538268397D-07 )*Y+1.24621276265907D-06 )*Y-          &
           4.30868944351523D-06 )*Y+9.94307982432868D-06
      UAUX=U
      WAUX=W
      RETURN
!
  180 IF (X .GT. 25.0D+00) GO TO 220
      IF (X .GT. 20.0D+00) GO TO 200
!     X=15.0 TO 20.0                             NROOTS = 5
      Y = X-17.5D+00
      RT1 = (((((((((( 1.91875764545740D-16*Y+7.8357401095707D-16)*Y- &
           3.260875931644D-14)*Y-1.186752035569D-13)*Y+               &
           4.275180095653D-12)*Y+3.357056136731D-11)*Y-               &
           1.123776903884D-09)*Y+1.231203269887D-08)*Y-               &
           3.99851421361031D-07 )*Y+1.45418822817771D-05 )*Y-         &
           3.49912254976317D-04 )*Y+6.67768703938812D-03              
      RT2 = (((((((((( 2.02778478673555D-15*Y+1.01640716785099D-14)*Y- &
           3.385363492036D-13)*Y-1.615655871159D-12)*Y+                &
           4.527419140333D-11)*Y+3.853670706486D-10)*Y-                &
           1.184607130107D-08)*Y+1.347873288827D-07)*Y-                &
           4.47788241748377D-06 )*Y+1.54942754358273D-04 )*Y-          &
           3.55524254280266D-03 )*Y+6.44912219301603D-02               
      RT3 = (((((((((( 7.79850771456444D-15*Y+6.00464406395001D-14)*Y- &
           1.249779730869D-12)*Y-1.020720636353D-11)*Y+               &
           1.814709816693D-10)*Y+1.766397336977D-09)*Y-               &
           4.603559449010D-08)*Y+5.863956443581D-07)*Y-               &
           2.03797212506691D-05 )*Y+6.31405161185185D-04 )*Y-         &
           1.30102750145071D-02 )*Y+2.10244289044705D-01              
      RT4 = (((((((((((-2.92397030777912D-15*Y+1.94152129078465D-14)*Y+ &
           4.859447665850D-13)*Y-3.217227223463D-12)*Y-                 &
           7.484522135512D-11)*Y+7.19101516047753D-10 )*Y+              &
           6.88409355245582D-09 )*Y-1.44374545515769D-07 )*Y+          &
           2.74941013315834D-06 )*Y-1.02790452049013D-04 )*Y+          &
           2.59924221372643D-03 )*Y-4.35712368303551D-02 )*Y+          &
           5.62170709585029D-01
      RT5 = ((((((((((( 1.17976126840060D-14*Y+1.24156229350669D-13)*Y- &
           3.892741622280D-12)*Y-7.755793199043D-12)*Y+                 &
           9.492190032313D-10)*Y-4.98680128123353D-09 )*Y-              &
           1.81502268782664D-07 )*Y+2.69463269394888D-06 )*Y+          &
           2.50032154421640D-05 )*Y-1.33684303917681D-03 )*Y+          &
           2.29121951862538D-02 )*Y-2.45653725061323D-01 )*Y+          &
           1.89999883453047D+00
      WW1 = (((((((((( 1.74841995087592D-15*Y-6.95671892641256D-16)*Y- &
          3.000659497257D-13)*Y+2.021279817961D-13)*Y+                &
          3.853596935400D-11)*Y+1.461418533652D-10)*Y-                &
          1.014517563435D-08)*Y+1.132736008979D-07)*Y-                &
          2.86605475073259D-06 )*Y+1.21958354908768D-04 )*Y-          &
          3.86293751153466D-03 )*Y+1.45298342081522D-01               
      WW2 = ((((((((((-1.11199320525573D-15*Y+1.85007587796671D-15)*Y+ &
          1.220613939709D-13)*Y+1.275068098526D-12)*Y-                &
          5.341838883262D-11)*Y+6.161037256669D-10)*Y-                &
          1.009147879750D-08)*Y+2.907862965346D-07)*Y-                &
          6.12300038720919D-06 )*Y+1.00104454489518D-04 )*Y-          &
          1.80677298502757D-03 )*Y+5.78009914536630D-02               
      WW3 = ((((((((((-9.49816486853687D-16*Y+6.67922080354234D-15)*Y+ &
          2.606163540537D-15)*Y+1.983799950150D-12)*Y-                &
          5.400548574357D-11)*Y+6.638043374114D-10)*Y-                &
          8.799518866802D-09)*Y+1.791418482685D-07)*Y-                &
          2.96075397351101D-06 )*Y+3.38028206156144D-05 )*Y-           &
          3.58426847857878D-04 )*Y+8.39213709428516D-03                
      WW4 = ((((((((((( 1.33829971060180D-17*Y-3.44841877844140D-16)*Y+ &
          4.745009557656D-15)*Y-6.033814209875D-14)*Y+                &
          1.049256040808D-12)*Y-1.70859789556117D-11 )*Y+             &
          2.15219425727959D-10 )*Y-2.52746574206884D-09 )*Y+          &
          3.27761714422960D-08 )*Y-3.90387662925193D-07 )*Y+          &
          3.46340204593870D-06 )*Y-2.43236345136782D-05 )*Y+          &
          3.54846978585226D-04                                        
      WW5 = ((((((((((((( 2.69412277020887D-20*Y-4.24837886165685D-19)* &
          Y+6.030500065438D-18)*Y-9.069722758289D-17)*Y+               &
          1.246599177672D-15)*Y-1.56872999797549D-14 )*Y+              &
          1.87305099552692D-13 )*Y-2.09498886675861D-12 )*Y+           &
          2.11630022068394D-11 )*Y-1.92566242323525D-10 )*Y+           &
          1.62012436344069D-09 )*Y-1.23621614171556D-08 )*Y+           &
          7.72165684563049D-08 )*Y-3.59858901591047D-07 )*Y+          &
          2.43682618601000D-06                                        
      UAUX=U
      WAUX=W
      RETURN                                                           
!
!     X=20.0 TO 25.0                             NROOTS = 5
  200 Y = X-22.5D+00
      RT1 = (((((((((-1.13927848238726D-15*Y+7.39404133595713D-15)*Y+  &
            1.445982921243D-13)*Y-2.676703245252D-12)*Y+                &
            5.823521627177D-12)*Y+2.17264723874381D-10 )*Y+             &
            3.56242145897468D-09 )*Y-3.03763737404491D-07 )*Y+          &
            9.46859114120901D-06 )*Y-2.30896753853196D-04 )*Y+          &
            5.24663913001114D-03                                        
      RT2 = (((((((((( 2.89872355524581D-16*Y-1.22296292045864D-14)*Y+ &
           6.184065097200D-14)*Y+1.649846591230D-12)*Y-                &
           2.729713905266D-11)*Y+3.709913790650D-11)*Y+                &
           2.216486288382D-09)*Y+4.616160236414D-08)*Y-                &
           3.32380270861364D-06 )*Y+9.84635072633776D-05 )*Y-          &
           2.30092118015697D-03 )*Y+5.00845183695073D-02               
      RT3 = (((((((((( 1.97068646590923D-15*Y-4.89419270626800D-14)*Y+ &
           1.136466605916D-13)*Y+7.546203883874D-12)*Y-                &
           9.635646767455D-11)*Y-8.295965491209D-11)*Y+                &
           7.534109114453D-09)*Y+2.699970652707D-07)*Y-                &
           1.42982334217081D-05 )*Y+3.78290946669264D-04 )*Y-          &
           8.03133015084373D-03 )*Y+1.58689469640791D-01               
      RT4 = (((((((((( 1.33642069941389D-14*Y-1.55850612605745D-13)*Y- &
           7.522712577474D-13)*Y+3.209520801187D-11)*Y-                &
           2.075594313618D-10)*Y-2.070575894402D-09)*Y+                &
           7.323046997451D-09)*Y+1.851491550417D-06)*Y-                &
           6.37524802411383D-05 )*Y+1.36795464918785D-03 )*Y-          &
           2.42051126993146D-02 )*Y+3.97847167557815D-01               
      RT5 = ((((((((((-6.07053986130526D-14*Y+1.04447493138843D-12)*Y- &
           4.286617818951D-13)*Y-2.632066100073D-10)*Y+                &
           4.804518986559D-09)*Y-1.835675889421D-08)*Y-                &
           1.068175391334D-06)*Y+3.292234974141D-05)*Y-                &
           5.94805357558251D-04 )*Y+8.29382168612791D-03 )*Y-          &
           9.93122509049447D-02 )*Y+1.09857804755042D+00               
      WW1 = (((((((((-9.10338640266542D-15*Y+1.00438927627833D-13)*Y+  &
           7.817349237071D-13)*Y-2.547619474232D-11)*Y+                &
           1.479321506529D-10)*Y+1.52314028857627D-09 )*Y+             &
           9.20072040917242D-09 )*Y-2.19427111221848D-06 )*Y+          &
           8.65797782880311D-05 )*Y-2.82718629312875D-03 )*Y+          &
           1.28718310443295D-01                                        
      WW2 = ((((((((( 5.52380927618760D-15*Y-6.43424400204124D-14)*Y-  &
           2.358734508092D-13)*Y+8.261326648131D-12)*Y+                &
           9.229645304956D-11)*Y-5.68108973828949D-09 )*Y+             &
           1.22477891136278D-07 )*Y-2.11919643127927D-06 )*Y+          &
           4.23605032368922D-05 )*Y-1.14423444576221D-03 )*Y+          &
           5.06607252890186D-02                                        
      WW3 = ((((((((( 3.99457454087556D-15*Y-5.11826702824182D-14)*Y-  &
           4.157593182747D-14)*Y+4.214670817758D-12)*Y+                &
           6.705582751532D-11)*Y-3.36086411698418D-09 )*Y+             &
           6.07453633298986D-08 )*Y-7.40736211041247D-07 )*Y+          &
           8.84176371665149D-06 )*Y-1.72559275066834D-04 )*Y+          &
           7.16639814253567D-03                                        
      WW4 = (((((((((((-2.14649508112234D-18*Y-2.45525846412281D-18)*Y+ &
           6.126212599772D-16)*Y-8.526651626939D-15)*Y+                 &
           4.826636065733D-14)*Y-3.39554163649740D-13 )*Y+              &
           1.67070784862985D-11 )*Y-4.42671979311163D-10 )*Y+           &
           6.77368055908400D-09 )*Y-7.03520999708859D-08 )*Y+           &
           6.04993294708874D-07 )*Y-7.80555094280483D-06 )*Y+           &
           2.85954806605017D-04
      WW5 = ((((((((((((-5.63938733073804D-21*Y+6.92182516324628D-20)*Y- &
           1.586937691507D-18)*Y+3.357639744582D-17)*Y-                  &
           4.810285046442D-16)*Y+5.386312669975D-15)*Y-                  &
           6.117895297439D-14)*Y+8.441808227634D-13)*Y-                 &
           1.18527596836592D-11 )*Y+1.36296870441445D-10 )*Y-           &
           1.17842611094141D-09 )*Y+7.80430641995926D-09 )*Y-           &
           5.97767417400540D-08 )*Y+1.65186146094969D-06
      UAUX=U
      WAUX=W
      RETURN
!
  220 WW1 = SQRT(PIE4/X)
      IF (X .GT. 40.0D+00) GO TO 240
!     X=25.0 TO 40.0                             NROOTS = 5
      E = EXP(-X)
      RT1 = ((((((((-1.73363958895356D-06*X+1.19921331441483D-04)*X -  &
           1.59437614121125D-02)*X+1.13467897349442D+00)*X -           &
           4.47216460864586D+01)*X+1.06251216612604D+03)*X -           &
           1.52073917378512D+04)*X+1.20662887111273D+05)*X -          &
           4.07186366852475D+05)*E + R15/(X-R15)                      
      RT2 = ((((((((-1.60102542621710D-05*X+1.10331262112395D-03)*X - &
           1.50043662589017D-01)*X+1.05563640866077D+01)*X -          &
           4.10468817024806D+02)*X+9.62604416506819D+03)*X -          &
           1.35888069838270D+05)*X+1.06107577038340D+06)*X -          &
           3.51190792816119D+06)*E + R25/(X-R25)
      RT3 = ((((((((-4.48880032128422D-05*X+2.69025112122177D-03)*X -  &
           4.01048115525954D-01)*X+2.78360021977405D+01)*X -           &
           1.04891729356965D+03)*X+2.36985942687423D+04)*X -           &
           3.19504627257548D+05)*X+2.34879693563358D+06)*X -           &
           7.16341568174085D+06)*E + R35/(X-R35)                       
      RT4 = ((((((((-6.38526371092582D-05*X-2.29263585792626D-03)*X -  &
           7.65735935499627D-02)*X+9.12692349152792D+00)*X -          &
           2.32077034386717D+02)*X+2.81839578728845D+02)*X +          &
           9.59529683876419D+04)*X-1.77638956809518D+06)*X +          &
           1.02489759645410D+07)*E + R45/(X-R45)
      RT5 = ((((((((-3.59049364231569D-05*X-2.25963977930044D-02)*X +  &
           1.12594870794668D+00)*X-4.56752462103909D+01)*X +           &
           1.05804526830637D+03)*X-1.16003199605875D+04)*X -           &
           4.07297627297272D+04)*X+2.22215528319857D+06)*X -           &
           1.61196455032613D+07)*E + R55/(X-R55)                       
      WW5 = (((((((((-4.61100906133970D-10*X+1.43069932644286D-07)*X - &
           1.63960915431080D-05)*X+1.15791154612838D-03)*X -          &
           5.30573476742071D-02)*X+1.61156533367153D+00)*X -          &
           3.23248143316007D+01)*X+4.12007318109157D+02)*X -          &
           3.02260070158372D+03)*X+9.71575094154768D+03)*E + W55*WW1
      WW4 = (((((((((-2.40799435809950D-08*X+8.12621667601546D-06)*X - &
           9.04491430884113D-04)*X+6.37686375770059D-02)*X -           &
           2.96135703135647D+00)*X+9.15142356996330D+01)*X -           &
           1.86971865249111D+03)*X+2.42945528916947D+04)*X -           &
           1.81852473229081D+05)*X+5.96854758661427D+05)*E + W45*WW1   
      WW3 = (((((((( 1.83574464457207D-05*X-1.54837969489927D-03)*X +  &
           1.18520453711586D-01)*X-6.69649981309161D+00)*X +          &
           2.44789386487321D+02)*X-5.68832664556359D+03)*X +          &
           8.14507604229357D+04)*X-6.55181056671474D+05)*X +          &
           2.26410896607237D+06)*E + W35*WW1
      WW2 = (((((((( 2.77778345870650D-05*X-2.22835017655890D-03)*X + &
           1.61077633475573D-01)*X-8.96743743396132D+00)*X +          &
           3.28062687293374D+02)*X-7.65722701219557D+03)*X +          &
           1.10255055017664D+05)*X-8.92528122219324D+05)*X +          &
           3.10638627744347D+06)*E + W25*WW1                          
      WW1 = WW1-0.01962D+00*E-WW2-WW3-WW4-WW5                         
      UAUX=U
      WAUX=W
      RETURN
!
  240 IF (X .GT. 59.0D+00) GO TO 260
!     X=40.0 TO 59.0                             NROOTS = 5
      XXX = X**3
      E = XXX*EXP(-X)
      RT1 = (((-2.43758528330205D-02*X+2.07301567989771D+00)*X -      &
           6.45964225381113D+01)*X+7.14160088655470D+02)*E + R15/(X-R15)
      RT2 = (((-2.28861955413636D-01*X+1.93190784733691D+01)*X -      &
           5.99774730340912D+02)*X+6.61844165304871D+03)*E + R25/(X-R25)
      RT3 = (((-6.95053039285586D-01*X+5.76874090316016D+01)*X -      &
           1.77704143225520D+03)*X+1.95366082947811D+04)*E + R35/(X-R35)
      RT4 = (((-1.58072809087018D+00*X+1.27050801091948D+02)*X -      &
           3.86687350914280D+03)*X+4.23024828121420D+04)*E + R45/(X-R45)
      RT5 = (((-3.33963830405396D+00*X+2.51830424600204D+02)*X -      &
           7.57728527654961D+03)*X+8.21966816595690D+04)*E + R55/(X-R55)
      E = XXX*E
      WW5 = (( 1.35482430510942D-08*X-3.27722199212781D-07)*X + &
           2.41522703684296D-06)*E + W55*WW1
      WW4 = (( 1.23464092261605D-06*X-3.55224564275590D-05)*X +  &
           3.03274662192286D-04)*E + W45*WW1
      WW3 = (( 1.34547929260279D-05*X-4.19389884772726D-04)*X +   &
           3.87706687610809D-03)*E + W35*WW1
      WW2 = (( 2.09539509123135D-05*X-6.87646614786982D-04)*X +  &
           6.68743788585688D-03)*E + W25*WW1
      WW1 = WW1-WW2-WW3-WW4-WW5
      UAUX=U
      WAUX=W
      RETURN
!
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
!-----------------------------------------------------------------------  
      UAUX=U
      WAUX=W
      RETURN
      END    
! ROOT6NOF
      SUBROUTINE ROOT6NOF(XX,NROOTS,UF,WF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      INTEGER,PARAMETER::MR=13,MXAUX=55
      DOUBLE PRECISION,PARAMETER::ONE=1.0D+00
!
      DOUBLE PRECISION,INTENT(IN)::XX
      INTEGER,INTENT(IN)::NROOTS
      DOUBLE PRECISION,DIMENSION(13),INTENT(OUT)::UF,WF
!
      DOUBLE PRECISION,DIMENSION(0:MR-1)::ALPHA,BETA
      DOUBLE PRECISION,DIMENSION(MR)::RTS,WTS,WRK
      DOUBLE PRECISION,DIMENSION(MR)::XASYMP
      DOUBLE PRECISION,DIMENSION(MR,MR)::RTSASY,WTSASY
      DOUBLE PRECISION,DIMENSION(55,8)::RTSAUX,WTSAUX
      DOUBLE PRECISION,DIMENSION(MXAUX)::RGRID,WGRID,P0,P1,P2
      INTEGER,DIMENSION(MR)::NAUXS,MAPRYS
!
!           GENERAL CASE (NROOTS=1 TO 13) OF RYS ROOTS/WEIGHTS
!
      EPS=1.0D-14
!
      CALL SETRYSNOF(XASYMP,RTSASY,WTSASY,RTSAUX,WTSAUX,NAUXS,MAPRYS)
!
      IF(XX.GE.XASYMP(NROOTS)) THEN
         CALL RYSASYNOF(NROOTS,XX,RTSASY,WTSASY,RTS,WTS)
!
      ELSE
         NAUX=NAUXS(NROOTS)
         MAP=MAPRYS(NROOTS)
         DO I=1,NAUX
            T2 = RTSAUX(I,MAP)*RTSAUX(I,MAP)
            RGRID(I) = T2
            WGRID(I) = WTSAUX(I,MAP)*EXP(-XX*T2)
         ENDDO
         CALL RYSDSNOF(NROOTS,NAUX,RGRID,WGRID,ALPHA,BETA,IERR,P0,P1,P2)
         IF (IERR.NE.0) THEN
            WRITE(6,*)"ERROR IN SUBROUTINE RYSDS, USED IN ROOT6,grad.f90"
            STOP
         ENDIF         
         CALL RYSGWNOF(NROOTS,ALPHA,BETA,EPS,RTS,WTS,IERR,WRK)
         IF (IERR.NE.0) THEN
            WRITE(6,*)"ERROR IN SUBROUTINE RYSGW, USED IN ROOT6,grad.f90"
            STOP
         ENDIF
      END IF
!
      DO K=1,NROOTS
         DUM  =RTS(K)
         UF(K)=DUM/(ONE-DUM)
         WF(K)=WTS(K)
      ENDDO
!-----------------------------------------------------------------------      
      RETURN
      END
! RYSASYNOF
      SUBROUTINE RYSASYNOF(NROOTS,XX,RTSASY,WTSASY,RTS,WTS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!      
      INTEGER,INTENT(IN)::NROOTS
      DOUBLE PRECISION,INTENT(IN)::XX
      DOUBLE PRECISION,DIMENSION(13,13),INTENT(IN)::RTSASY,WTSASY
      DOUBLE PRECISION,DIMENSION(NROOTS),INTENT(OUT)::RTS,WTS
!
      DOUBLE PRECISION,PARAMETER::ONE=1.0D+00
!
!     THE ROOTS AND WEIGHTS FOR XX = BIG ARE GIVEN BY:
!       T*T = S*S/XX    W = V/SQRT(XX)
!     WHERE S AND V ARE THE ROOTS AND WEIGHTS OF THE
!     HERMITE POLYNOMIALS, OF ORDER 2*NROOTS.
!
      FACTR = ONE/XX
      FACTW = SQRT(FACTR)
      DO I=1,NROOTS
         RTS(I)= FACTR * RTSASY(I,NROOTS)
         WTS(I)= FACTW * WTSASY(I,NROOTS)
      ENDDO
!-----------------------------------------------------------------------      
      RETURN
      END
! SETRYSNOF
      SUBROUTINE SETRYSNOF(XASYMP,RTSASY,WTSASY,RTSAUX,WTSAUX,NAUXS,MAPRYS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!      
      INTEGER,PARAMETER::MR=13, MAUX=55
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00
      DOUBLE PRECISION,DIMENSION(MAUX)::RTS,WTS,WRK
      DOUBLE PRECISION,DIMENSION(0:MAUX-1)::ALPHA,BETA
!
      DOUBLE PRECISION,DIMENSION(MR),INTENT(OUT)::XASYMP
      DOUBLE PRECISION,DIMENSION(MR,MR),INTENT(OUT)::RTSASY,WTSASY
      DOUBLE PRECISION,DIMENSION(55,8),INTENT(OUT)::RTSAUX,WTSAUX
      INTEGER,DIMENSION(MR),INTENT(OUT)::NAUXS,MAPRYS
!
!     ----- INITIALIZE THE RYS QUADRATURE PROCEDURE -----
!
      QUART = ONE/FOUR
      PI    = 3.141592653589793238D+00
      EPS   = 1.0D-14
!
!     SET UP ASYMPTOTIC ROOT COMPUTATION, BY GENERATING
!     THE ROOTS OF THE HERMITE POLYNOMIALS OF ORDER 2N.
!
      XASYMP( 1)=29.0D+00
      XASYMP( 2)=37.0D+00
      XASYMP( 3)=43.0D+00
      XASYMP( 4)=49.0D+00
      XASYMP( 5)=55.0D+00
      XASYMP( 6)=60.0D+00
      XASYMP( 7)=65.0D+00
      XASYMP( 8)=71.0D+00
      XASYMP( 9)=76.0D+00
      XASYMP(10)=81.0D+00
      XASYMP(11)=86.0D+00
      XASYMP(12)=91.0D+00
      XASYMP(13)=96.0D+00
!
      DO I=1,MR
!            NOTE THAT MAUX MUST BE AT LEAST TWICE MR, FROM NEXT LINE
         N=2*I
         DO J=0,N-1
            ALPHA(J) = ZERO
         ENDDO
         BETA(0)=SQRT(PI)
         DO J=1,N-1
            BETA(J) = J/TWO
         END DO
         CALL RYSGWNOF(N,ALPHA,BETA,EPS,RTS,WTS,IERR,WRK)
         DO J=1,I
            RTSASY(J,I) = RTS(I+J)*RTS(I+J)
            WTSASY(J,I) = WTS(I+J)
         ENDDO
      ENDDO
!
!        GENERATE AUXILIARY GRIDS, AT 8 DISTINCT POINT DENSITIES
!
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
!
!        THE AUXILIARY QUADRATURE IS TAKEN TO BE "SHIFTED LEGENDRE"
!
      NAUXSV=0
      IGRID=0
      DO M=1,MR
         NAUX=NAUXS(M)
         IF(NAUX.EQ.NAUXSV) CYCLE
         IGRID=IGRID+1
         NAUXSV=NAUX
         DO I=0,NAUX-1
            ALPHA(I) = ONE/TWO
         ENDDO
         BETA(0)= ONE
         DO I=1,NAUX-1
            BETA(I) = QUART/(FOUR-(ONE/(I*I)))
         END DO
         CALL RYSGWNOF(NAUX,ALPHA,BETA,EPS,RTS,WTS,IERR,WRK)
         DO I=1,NAUX
            RTSAUX(I,IGRID) = RTS(I)
            WTSAUX(I,IGRID) = WTS(I)
         ENDDO
      ENDDO
!-----------------------------------------------------------------------  
      RETURN
      END
! RYSGWNOF
      SUBROUTINE RYSGWNOF(N,ALPHA,BETA,EPS,ROOTS,WEIGHT,IERR,WRK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER,INTENT(IN)::N
      DOUBLE PRECISION,DIMENSION(N),INTENT(IN)::ALPHA,BETA
      DOUBLE PRECISION,INTENT(IN)::EPS
      DOUBLE PRECISION,DIMENSION(N),INTENT(OUT)::ROOTS,WEIGHT,WRK
      INTEGER,INTENT(OUT)::IERR
!
!        INPUT:  N - - THE NUMBER OF POINTS IN THE GAUSSIAN QUADRATURE
!                      FORMULA; TYPE INTEGER
!                ALPHA,BETA - - ARRAYS OF DIMENSION  N  TO BE FILLED
!                      WITH THE VALUES OF  ALPHA(K-1), BETA(K-1), K=1,2,
!                      ...,N
!                EPS - THE RELATIVE ACCURACY DESIRED IN THE NODES
!                      AND WEIGHTS
!
!        OUTPUT: ROOTS- ARRAY OF DIMENSION  N  CONTAINING THE GAUSSIAN
!                      NODES (IN INCREASING ORDER)  ROOTS(K)=X(K),
!                       K=1,2,...,N
!                WEIGHT - ARRAY OF DIMENSION  N  CONTAINING THE
!                      GAUSSIAN WEIGHTS  WEIGHT(K)=W(K), K=1,2,...,N
!                IERR- AN ERROR FLAG EQUAL TO  0  ON NORMAL RETURN,
!                      EQUAL TO  I  IF THE QR ALGORITHM DOES NOT
!                      CONVERGE WITHIN 30 ITERATIONS ON EVALUATING THE
!                      I-TH EIGENVALUE, EQUAL TO  -1  IF  N  IS NOT IN
!                      RANGE, AND EQUAL TO  -2  IF ONE OF THE BETA'S IS
!                      NEGATIVE.
!
! THE ARRAY  WRK  IS NEEDED FOR WORKING SPACE.
!
      IF(N.LT.1) THEN
        IERR=-1
        RETURN
      END IF
!
      IERR=0
      ROOTS(1)=ALPHA(1)
      IF(BETA(1).LT.0.0D+00) THEN
        IERR=-2
        RETURN
      END IF
      WEIGHT(1)=BETA(1)
      IF (N.EQ.1) RETURN
!
      WEIGHT(1)=1.0D+00
      WRK(N)=0.0D+00
      DO K=2,N
        ROOTS(K)=ALPHA(K)
        IF(BETA(K).LT.0.0D+00) THEN
          IERR=-2
          RETURN
        END IF
        WRK(K-1)=SQRT(BETA(K))
        WEIGHT(K)=0.0D+00
      ENDDO
!
      DO L=1,N
        J=0
!
! LOOK FOR A SMALL SUBDIAGONAL ELEMENT.
!
  105   DO M=L,N
          IF(M.EQ.N) EXIT
          IF(ABS(WRK(M)).LE.EPS*(ABS(ROOTS(M))+ABS(ROOTS(M+1)))) &
                   EXIT
        ENDDO
        DP=ROOTS(L)
        IF(M.EQ.L) CYCLE
        IF(J.EQ.30) GO TO 400
        J=J+1
!
! FORM SHIFT.
!
        DG=(ROOTS(L+1)-DP)/(2.0D+00*WRK(L))
        DR=SQRT(DG*DG+1.0D+00)
        DG=ROOTS(M)-DP+WRK(L)/(DG+SIGN(DR,DG))
        DS=1.0D+00
        DC=1.0D+00
        DP=0.0D+00
        MML=M-L
!
! FOR I=M-1 STEP -1 UNTIL L DO ...
!
        DO II=1,MML
          I=M-II
          DF=DS*WRK(I)
          DB=DC*WRK(I)
          IF(ABS(DF).LT.ABS(DG)) GO TO 150
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
!
! FORM FIRST COMPONENT OF VECTOR.
!
          DF=WEIGHT(I+1)
          WEIGHT(I+1)=DS*WEIGHT(I)+DC*DF
          WEIGHT(I)=DC*WEIGHT(I)-DS*DF
        ENDDO
        ROOTS(L)=ROOTS(L)-DP
        WRK(L)=DG
        WRK(M)=0.0D+00
        GO TO 105
      ENDDO
!
! ORDER EIGENVALUES AND EIGENVECTORS.
!
      DO II=2,N
        I=II-1
        K=I
        DP=ROOTS(I)
        DO J=II,N
          IF(ROOTS(J).GE.DP) CYCLE
          K=J
          DP=ROOTS(J)
        ENDDO
        IF(K.EQ.I) CYCLE
        ROOTS(K)=ROOTS(I)
        ROOTS(I)=DP
        DP=WEIGHT(I)
        WEIGHT(I)=WEIGHT(K)
        WEIGHT(K)=DP
      ENDDO
      DO K=1,N
        WEIGHT(K)=BETA(1)*WEIGHT(K)*WEIGHT(K)
      ENDDO
      RETURN
!
! SET ERROR - NO CONVERGENCE TO AN EIGENVALUE AFTER 30 ITERATIONS.
!
  400 IERR=L
!-----------------------------------------------------------------------  
      RETURN
      END
! RYSDSNOF
      SUBROUTINE RYSDSNOF(N,NCAP,X,W,ALPHA,BETA,IERR,P0,P1,P2)      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
      INTEGER,INTENT(IN)::N,NCAP
      INTEGER,INTENT(OUT)::IERR
      DOUBLE PRECISION,DIMENSION(NCAP),INTENT(IN)::X,W
      DOUBLE PRECISION,DIMENSION(N),INTENT(OUT)::ALPHA,BETA
      DOUBLE PRECISION,DIMENSION(NCAP),INTENT(OUT)::P0,P1,P2
!
!      ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE, 20, 21-62(1994)
!
! (CF. SECTION 2.1 OF W. GAUTSCHI, 'ON GENERATING ORTHOGONAL  SCI.
! POLYNOMIALS', SIAM J. STATIST. COMPUT. 3, 1982, 289-317)
!
!     (F,G)=SUM OVER K FROM 1 TO NCAP OF W(K)*F(X(K))*G(X(K)).
!
      TINYY = 1.0D-40
      HUGEE = 1.0D+40
!
      IERR=0
      IF(N.LE.0 .OR. N.GT.NCAP) THEN
        IERR=1
        RETURN
      END IF
      NM1=N-1
!
! COMPUTE THE FIRST ALPHA- AND BETA-COEFFICIENT.
!
      SUM0=0.0D+00
      SUM1=0.0D+00
      DO M=1,NCAP
        SUM0=SUM0+W(M)
        SUM1=SUM1+W(M)*X(M)
      ENDDO
      ALPHA(1)=SUM1/SUM0
      BETA(1)=SUM0
     IF(N.EQ.1) RETURN
!
! COMPUTE THE REMAINING ALPHA- AND BETA-COEFFICIENTS.
!
     DO M=1,NCAP
        P1(M)=0.0D+00
        P2(M)=1.0D+00
     ENDDO
      DO K=1,NM1
        SUM1=0.0D+00
        SUM2=0.0D+00
        DO M=1,NCAP
!
! THE FOLLOWING STATEMENT IS DESIGNED TO AVOID AN OVERFLOW CONDITION
! IN THE COMPUTATION OF  P2(M)  WHEN THE WEIGHTS  W(M)  GO TO ZERO
! FASTER (AND UNDERFLOW) THAN THE  P2(M)  GROW.
!
          IF(W(M).EQ.0.0D+00) CYCLE
          P0(M)=P1(M)
          P1(M)=P2(M)
          P2(M)=(X(M)-ALPHA(K))*P1(M)-BETA(K)*P0(M)
!
! CHECK FOR IMPENDING OVERFLOW.
!
          IF(ABS(P2(M)).GT.HUGEE .OR. ABS(SUM2).GT.HUGEE) THEN
            IERR=K
            RETURN
          END IF
          T=W(M)*P2(M)*P2(M)
          SUM1=SUM1+T
          SUM2=SUM2+T*X(M)
        ENDDO
!
! CHECK FOR IMPENDING UNDERFLOW.
!
        IF(ABS(SUM1).LT.TINYY) THEN
          IERR=-K
          RETURN
        END IF
        ALPHA(K+1)=SUM2/SUM1
        BETA(K+1)=SUM1/SUM0
        SUM0=SUM1
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
!-----------------------------------------------------------------------
!              ROUTINES RELATED TO GRADIENT COMPUTATION
!-----------------------------------------------------------------------
! OEDHNDNOF      
      SUBROUTINE OEDHNDNOF(DCHRG,NIJGDIM,NIJG,KTYPE,KLOC,KKMIN,KKMAX,  &
                        KSTART,KNG,KATOM,IA,CX0,CY0,CZ0,EX1,CS,CP,    &
                        CD,CF,CG)
      USE PARCOM 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(2,(NSHELL*(NSHELL+1))/2),INTENT(OUT)::NIJG
      DOUBLE PRECISION,DIMENSION(15,NIJGDIM),INTENT(OUT)::DCHRG
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KTYPE,KLOC,KKMIN,KKMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KSTART,KNG,KATOM
      INTEGER,DIMENSION(NBF),INTENT(IN)::IA
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN)::CX0,CY0,CZ0
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN)::EX1,CS,CP,CD,CF,CG
!
      NIJ0=0
      DO II=1,NSHELL
        DO JJ=1,II
          ISHI=II
          JSHJ=JJ
          CALL OEDSHLNOF(ISHI,JSHJ,DCHRG,NIJGDIM,NIJ,KTYPE,KLOC,KKMIN, &
                      KKMAX,KSTART,KNG,KATOM,CX0,CY0,CZ0,EX1,      &
                      CS,CP,CD,CF,CG,NIJ0)
          IIJJ=IA(MAX0(II,JJ))+MIN0(II,JJ)
          NIJG(1,IIJJ)=NIJ0
          NIJG(2,IIJJ)=NIJ
          NIJ0=NIJ0+NIJ
        ENDDO
      ENDDO
!      
      RETURN
      END
! OEDSHLNOF
      SUBROUTINE OEDSHLNOF(ISH,JSH,DCHRG,NIJGDIM,NIJ,KTYPE,KLOC,KKMIN, &
                        KKMAX,KSTART,KNG,KATOM,CX,CY,CZ,EX1,       &
                        CS,CP,CD,CF,CG,NIJ0)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
!
      DOUBLE PRECISION::RTOL,DTOL
      INTEGER,INTENT(IN)::ISH,JSH,NIJ0
      INTEGER,INTENT(OUT)::NIJ
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KTYPE,KLOC,KKMIN,KKMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KSTART,KNG,KATOM
      DOUBLE PRECISION,DIMENSION(15,NIJGDIM),INTENT(OUT)::DCHRG
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN)::CX,CY,CZ
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN)::EX1,CS,CP,CD,CF,CG
!     NOTE THAT THIS IS LIMITED TO G SHELL
      DOUBLE PRECISION,PARAMETER::ONE=1.0D+00
      DOUBLE PRECISION,DIMENSION(30)::GA,CCA,CCAS,GB,CCB,CCBS
      DOUBLE PRECISION,DIMENSION(NPRIMI,5)::CSPDFG
      LOGICAL::IIEQJJ,SPI,SPJ,EXPNDI
!      
      IIEQJJ=ISH.EQ.JSH
      RTOL=20*2.30258D+00
      DTOL=(10.0D+00)**(-20)
      CSPDFG(:,1)=CS(:)
      CSPDFG(:,2)=CP(:)
      CSPDFG(:,3)=CD(:)
      CSPDFG(:,4)=CF(:)
      CSPDFG(:,5)=CG(:)
!
!     ----- ISHELL -----
!
      I=KATOM(ISH)
      XI=CX(I)
      YI=CY(I)
      ZI=CZ(I)
      I1=KSTART(ISH)
      I2=I1+KNG(ISH)-1
      LIT=KTYPE(ISH)
      MINI=KKMIN(ISH)
      MAXI=KKMAX(ISH)
      NUMI=MAXI-MINI+1
      LOCI=KLOC(ISH)-MINI
      SPI=LIT.EQ.2.AND.MINI.EQ.1
      NGA=0
      DO I=I1,I2
        NGA=NGA+1
        GA(NGA)=EX1(I)
        CCA(NGA)=CSPDFG(I,LIT)
        IF(SPI) CCAS(NGA)=CSPDFG(I,1)/CSPDFG(I,2)
      ENDDO
!
!     ----- JSHELL -----
!
      J=KATOM(JSH)
      XJ=CX(J)
      YJ=CY(J)
      ZJ=CZ(J)
      J1=KSTART(JSH)
      J2=J1+KNG(JSH)-1
      LJT=KTYPE(JSH)
      MINJ=KKMIN(JSH)
      MAXJ=KKMAX(JSH)
      NUMJ=MAXJ-MINJ+1
      LOCJ=KLOC(JSH)-MINJ
      SPJ=LJT.EQ.2.AND.MINJ.EQ.1
      NGB=0
      DO J=J1,J2
        NGB=NGB+1
        GB(NGB)=EX1(J)
        CCB(NGB)=CSPDFG(J,LJT)
        IF(SPJ) CCBS(NGB)=CSPDFG(J,1)/CSPDFG(J,2)
      ENDDO
      RRI=((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
      EXPNDI=LIT.GE.LJT
!
!     ----- -IJ- CHARGE DISTRIBUTION -----
!
      XC=XI
      YC=YI
      ZC=ZI
      DXIJ=XI-XJ
      DYIJ=YI-YJ
      DZIJ=ZI-ZJ
      IF(.NOT.EXPNDI) THEN
       XC=XJ
       YC=YJ
       ZC=ZJ
       DXIJ=XJ-XI
       DYIJ=YJ-YI
       DZIJ=ZJ-ZI
      ENDIF
!
!     ----- - I- PRIMITIVE           -----
!
      NIJ=0
      DO IA=1,NGA
        AI=GA(IA)
        ARRI=AI*RRI
        AXI=AI*XI
        AYI=AI*YI
        AZI=AI*ZI
        CCI=CCA(IA)
!
!     ----- - J- PRIMITIVE           -----
!
        DO JB=1,NGB
          AJ=GB(JB)
          AA=AI+AJ
          AA1=ONE/AA
          DUM=AJ*ARRI*AA1
          IF(DUM.GT.RTOL) CYCLE
          DAEXPA=CCI*CCB(JB)* EXP(-DUM)*AA1
          DUM=  ABS(DAEXPA)
          IF(DUM.LE.DTOL) CYCLE
!
          NIJ=NIJ+1
          DCHRG( 1,NIJ+NIJ0)= DAEXPA
          DCHRG( 2,NIJ+NIJ0)= AA
          DCHRG( 3,NIJ+NIJ0)=(AXI+AJ*XJ)*AA1
          DCHRG( 4,NIJ+NIJ0)=(AYI+AJ*YJ)*AA1
          DCHRG( 5,NIJ+NIJ0)=(AZI+AJ*ZJ)*AA1
          DCHRG( 6,NIJ+NIJ0)= XC
          DCHRG( 7,NIJ+NIJ0)= YC
          DCHRG( 8,NIJ+NIJ0)= ZC
          DCHRG( 9,NIJ+NIJ0)= DXIJ
          DCHRG(10,NIJ+NIJ0)= DYIJ
          DCHRG(11,NIJ+NIJ0)= DZIJ
          DCHRG(12,NIJ+NIJ0)= AI+AI
          DCHRG(13,NIJ+NIJ0)= AJ+AJ
          IF(SPI) DCHRG(14,NIJ+NIJ0)=CCAS(IA)
          IF(SPJ) DCHRG(15,NIJ+NIJ0)=CCBS(JB)
        ENDDO
      ENDDO     
      RETURN
      END      
! OEDRDNOF      
      SUBROUTINE OEDRDNOF(NIJG,NIJ,NIJ0,IIJJ)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
!
      INTEGER,INTENT(IN)::IIJJ
      INTEGER,DIMENSION(2,(NSHELL*(NSHELL+1))/2),INTENT(IN)::NIJG
      INTEGER,INTENT(OUT)::NIJ0,NIJ
!
      NIJ0=NIJG(1,IIJJ)
      NIJ =NIJG(2,IIJJ)
      RETURN
      END      
! JKDATMNOF
      SUBROUTINE JKDATMNOF(II,JJ,KK,LL,SKIPI,SKIPJ,SKIPK,SKIPL,INVTYP, &
                        KATOM,IIAT,JJAT,KKAT,LLAT)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KATOM
      INTEGER,INTENT(IN)::II,JJ,KK,LL
      INTEGER,INTENT(OUT)::INVTYP,IIAT,JJAT,KKAT,LLAT
      LOGICAL,INTENT(OUT)::SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL::IANDJ,IANDK,IANDL,JANDK,JANDL,KANDL
!-----------------------------------------------------------------------      
      SKIPI=.TRUE.
      SKIPJ=.TRUE.
      SKIPK=.TRUE.
      SKIPL=.TRUE.
      IIAT=KATOM(II)
      JJAT=KATOM(JJ)
      KKAT=KATOM(KK)
      LLAT=KATOM(LL)
      IANDJ=IIAT.EQ.JJAT
      IANDK=IIAT.EQ.KKAT
      IANDL=IIAT.EQ.LLAT
      JANDK=JJAT.EQ.KKAT
      JANDL=JJAT.EQ.LLAT
      KANDL=KKAT.EQ.LLAT
      IF(.NOT.IANDJ) GO TO 500
      IF(.NOT.IANDK) GO TO 200
      IF(.NOT.IANDL) GO TO 100
!     ----- IAT = JAT = KAT = LAT ----- (IAT,IAT/IAT,IAT) -----
      INVTYP=1
      GO TO 1500
  100 CONTINUE
!     ----- IAT = JAT = KAT ; LAT ----- (IAT,IAT/IAT,LAT) -----
      SKIPL=.FALSE.
      INVTYP=2
      GO TO 1500
  200 IF(.NOT.IANDL) GO TO 300
!     ----- IAT = JAT = LAT ; KAT ----- (IAT,IAT/KAT,IAT) -----
      SKIPK=.FALSE.
      INVTYP=3
      GO TO 1500
  300 IF(.NOT.KANDL) GO TO 400
!     ----- IAT = JAT ; KAT = LAT ----- (IAT,IAT/KAT,KAT) -----
      SKIPK=.FALSE.
      SKIPL=.FALSE.
      INVTYP=4
      GO TO 1500
  400 CONTINUE
!     ----- IAT = JAT ; KAT ; LAT ----- (IAT,IAT/KAT,LAT) -----
      SKIPK=.FALSE.
      SKIPL=.FALSE.
      INVTYP=5
      GO TO 1500
  500 IF(.NOT.IANDK) GO TO 800
      IF(.NOT.IANDL) GO TO 600
!     ----- IAT = KAT = LAT ; JAT ----- (IAT,JAT/IAT,IAT) -----
      SKIPJ=.FALSE.
      INVTYP=6
      GO TO 1500
  600 IF(.NOT.JANDL) GO TO 700
!     ----- IAT = KAT ; JAT = LAT ----- (IAT,JAT/IAT,JAT) -----
      SKIPJ=.FALSE.
      SKIPL=.FALSE.
      INVTYP=7
      GO TO 1500
  700 CONTINUE
!     ----- IAT = KAT ; JAT ; LAT ----- (IAT,JAT/IAT,LAT) -----
      SKIPJ=.FALSE.
      SKIPL=.FALSE.
      INVTYP=8
      GO TO 1500
  800 IF(.NOT.IANDL) GO TO 1000
      IF(.NOT.JANDK) GO TO 900
!     ----- IAT = LAT ; JAT = KAT ----- (IAT,JAT/JAT,IAT) -----
      SKIPJ=.FALSE.
      SKIPK=.FALSE.
      INVTYP=9
      GO TO 1500
  900 CONTINUE
!     ----- IAT = LAT ; JAT , KAT ----- (IAT,JAT/KAT,IAT) -----
      SKIPJ=.FALSE.
      SKIPK=.FALSE.
      INVTYP=10
      GO TO 1500
 1000 IF(.NOT.JANDK) GO TO 1200
      IF(.NOT.JANDL) GO TO 1100
!     ----- IAT ; JAT = JAT = JAT ----- (IAT,JAT/JAT,JAT) -----
      SKIPI=.FALSE.
      INVTYP=11
      GO TO 1500
 1100 CONTINUE
!     ----- IAT ; JAT = KAT ; LAT ----- (IAT,JAT/JAT,LAT) -----
      SKIPI=.FALSE.
      SKIPL=.FALSE.
      INVTYP=12
      GO TO 1500
 1200 IF(.NOT.JANDL) GO TO 1300
!     ----- JAT = LAT ; IAT ; KAT ----- (IAT,JAT/KAT,JAT) -----
      SKIPI=.FALSE.
      SKIPK=.FALSE.
      INVTYP=13
      GO TO 1500
 1300 IF(.NOT.KANDL) GO TO 1400
!     ----- KAT = LAT ; IAT ; JAT ----- (IAT,JAT/KAT,KAT) -----
      SKIPI=.FALSE.
      SKIPJ=.FALSE.
      INVTYP=14
      GO TO 1500
 1400 CONTINUE
!     ----- IAT ; JAT ; KAT ; LAT ----- (IAT,JAT/KAT,LAT) -----
      SKIPI=.FALSE.
      SKIPJ=.FALSE.
      SKIPK=.FALSE.
      INVTYP=15
 1500 CONTINUE      
!-----------------------------------------------------------------------
      RETURN
      END
! JKDSHLNOF
      SUBROUTINE JKDSHLNOF(ISH,JSH,KSH,LSH,KTYPE,KLOC,KKMIN,KKMAX,     &
                        SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL,EXPNDI,EXPNDK &
                        ,MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL,      &
                        IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL,            &
                        LIT,LJT,LKT,LLTT,NUMI,NUMJ,NUMK,NUML)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,INTENT(IN)::ISH,JSH,KSH,LSH
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KTYPE,KLOC,KKMIN,KKMAX
      INTEGER,INTENT(OUT)::MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL
      INTEGER,INTENT(OUT)::LIT,LJT,LKT,LLTT
!     USE LLTT INSTEAD OF LLT BECAUSE LLT IS INTRINSIC FUNCTION      
      INTEGER,INTENT(OUT)::NUMI,NUMJ,NUMK,NUML
      LOGICAL,INTENT(OUT)::SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL      
      LOGICAL,INTENT(OUT)::IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      LOGICAL,INTENT(OUT)::EXPNDI,EXPNDK      
!-----------------------------------------------------------------------
      IIEQJJ=ISH.EQ.JSH
      KKEQLL=KSH.EQ.LSH
      IJEQKL=ISH.EQ.KSH.AND.JSH.EQ.LSH
      IJGTKL=MAX0(ISH,JSH).GT.MAX0(KSH,LSH)
      IJLTKL=MAX0(ISH,JSH).LT.MAX0(KSH,LSH)
!     ----- ISHELL -----
      LIT=KTYPE(ISH)
      MINI=KKMIN(ISH)
      MAXI=KKMAX(ISH)
      NUMI=MAXI-MINI+1
      LOCI=KLOC(ISH)-MINI
      SPI=LIT.EQ.2.AND.MINI.EQ.1
!     ----- JSHELL -----
      LJT=KTYPE(JSH)
      MINJ=KKMIN(JSH)
      MAXJ=KKMAX(JSH)
      NUMJ=MAXJ-MINJ+1
      LOCJ=KLOC(JSH)-MINJ
      SPJ=LJT.EQ.2.AND.MINJ.EQ.1
      SPIJ=SPI.OR.SPJ
      EXPNDI=LIT.GE.LJT
!     ----- KSHELL -----
      LKT=KTYPE(KSH)
      MINK=KKMIN(KSH)
      MAXK=KKMAX(KSH)
      NUMK=MAXK-MINK+1
      LOCK=KLOC(KSH)-MINK
      SPK=LKT.EQ.2.AND.MINK.EQ.1
!     ----- LSHELL -----
      LLTT=KTYPE(LSH)
      MINL=KKMIN(LSH)
      MAXL=KKMAX(LSH)
      NUML=MAXL-MINL+1
      LOCL=KLOC(LSH)-MINL
      SPL=LLTT.EQ.2.AND.MINL.EQ.1
      SPKL=SPK.OR.SPL
      SPIJKL=SPIJ.OR.SPKL
      EXPNDK=LKT.GE.LLTT
!-----------------------------------------------------------------------
      RETURN
      END
! JKDNDXNOF
      SUBROUTINE JKDNDXNOF(LIT,LJT,LKT,LLTT,NUMJ,NUMK,NUML,            &
                        SKIPI,SKIPJ,SKIPK,SKIPL,MINI,MAXI,       &
                        MINJ,MAXJ,MINK,MAXK,MINL,MAXL,NROOTS,    &
                        IGXYZ,JGXYZ,KGXYZ,LGXYZ,                 &
                        IIEQJJ,KKEQLL,IJEQKL,IJKLG,MAXNUM,       &
                        IDER,JDER,KDER,LDER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(4,MAXNUM)::IJKLG
      INTEGER,INTENT(IN)::LIT,LJT,LKT,LLTT,NUMJ,NUMK,NUML
      INTEGER,INTENT(IN)::MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL
      INTEGER,INTENT(OUT)::NROOTS,IDER,JDER,KDER,LDER
      INTEGER,DIMENSION(4,35),INTENT(OUT)::IGXYZ,JGXYZ,KGXYZ,LGXYZ
      LOGICAL,INTENT(IN)::SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL,INTENT(IN)::IIEQJJ,KKEQLL,IJEQKL
      INTEGER,DIMENSION(35)::IJKLX,IJKLY,IJKLZ
      INTEGER,DIMENSION(5)::IJKLN
      DATA IJKLN /   1,  4, 10, 20, 35/
      DATA IJKLX /   0,  1,  0,  0,  2,  0,  0,  1,  1,  0,      &
                     3,  0,  0,  2,  2,  1,  0,  1,  0,  1,      &
                     4,  0,  0,  3,  3,  1,  0,  1,  0,  2,      &
                     2,  0,  2,  1,  1/
      DATA IJKLY /   0,  0,  1,  0,  0,  2,  0,  1,  0,  1,      &
                     0,  3,  0,  1,  0,  2,  2,  0,  1,  1,      &
                     0,  4,  0,  1,  0,  3,  3,  0,  1,  2,      &
                     0,  2,  1,  2,  1/
      DATA IJKLZ /   0,  0,  0,  1,  0,  0,  2,  0,  1,  1,      &
                     0,  0,  3,  0,  1,  0,  1,  2,  2,  1,      &
                     0,  0,  4,  0,  1,  0,  1,  3,  3,  0,      &
                     2,  2,  1,  1,  2/
!-----------------------------------------------------------------------
      IDER=1
      JDER=1
      KDER=1
      LDER=1
      IF(SKIPI) IDER=0
      IF(SKIPJ) JDER=0
      IF(SKIPK) KDER=0
      IF(SKIPL) LDER=0
      LJTMOD=LJT + JDER
      LKTMOD=LKT + KDER
      LLTMOD=LLTT + LDER
!     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS -----
      NI=NUML*NUMK*NUMJ
      DO I=MINI,MAXI
        IGXYZ(1,I)=NI*(I-MINI)+1
      ENDDO
      LLKJT=LLTMOD*LKTMOD*LJTMOD
      DO I=1,IJKLN(LIT)
        IGXYZ(2,I)=IJKLX(I)*LLKJT+1
        IGXYZ(3,I)=IJKLY(I)*LLKJT+1
        IGXYZ(4,I)=IJKLZ(I)*LLKJT+1
      ENDDO
      NJ=NUML*NUMK
      DO J=MINJ,MAXJ
        JGXYZ(1,J)=NJ*(J-MINJ)
      ENDDO
      LLKT=LLTMOD*LKTMOD
      DO J=1,IJKLN(LJT)
        JGXYZ(2,J)=IJKLX(J)*LLKT
        JGXYZ(3,J)=IJKLY(J)*LLKT
        JGXYZ(4,J)=IJKLZ(J)*LLKT
      ENDDO
!     ----- PREPARE INDICES FOR PAIRS OF (K,L) FUNCTIONS -----
      NK=NUML
      DO K=MINK,MAXK
        KGXYZ(1,K)=NK*(K-MINK)
      ENDDO
      DO K=1,IJKLN(LKT)
        KGXYZ(2,K)=IJKLX(K)*LLTMOD
        KGXYZ(3,K)=IJKLY(K)*LLTMOD
        KGXYZ(4,K)=IJKLZ(K)*LLTMOD
      ENDDO
      NL=1
      DO L=MINL,MAXL
         LGXYZ(1,L)=NL*(L-MINL)
      ENDDO
      DO L=1,IJKLN(LLTT)
         LGXYZ(2,L)=IJKLX(L)
         LGXYZ(3,L)=IJKLY(L)
         LGXYZ(4,L)=IJKLZ(L)
      ENDDO
!     ----- PREPARE INDICES FOR (IJ/KL) -----
      IJKL=0
      DO I=MINI,MAXI
        JMAX=MAXJ
        IF(IIEQJJ) JMAX=I
        DO J=MINJ,JMAX
          KMAX=MAXK
          IF(IJEQKL) KMAX=I
          DO K=MINK,KMAX
            LMAX=MAXL
            IF(KKEQLL           ) LMAX=K
            IF(IJEQKL.AND.K.EQ.I) LMAX=J
            DO L=MINL,LMAX
              IJKL=IJKL+1
              NN=((IGXYZ(1,I)+JGXYZ(1,J))+KGXYZ(1,K))+LGXYZ(1,L)
              NX=((IGXYZ(2,I)+JGXYZ(2,J))+KGXYZ(2,K))+LGXYZ(2,L)
              NY=((IGXYZ(3,I)+JGXYZ(3,J))+KGXYZ(3,K))+LGXYZ(3,L)
              NZ=((IGXYZ(4,I)+JGXYZ(4,J))+KGXYZ(4,K))+LGXYZ(4,L)
              IJKLG(1,IJKL)=   NN
              IJKLG(2,IJKL)=3*(NX-1)+1
              IJKLG(3,IJKL)=3*(NY-1)+2
              IJKLG(4,IJKL)=3*(NZ-1)+3
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!     ----- SET NUMBER OF QUADRATURE POINTS -----
      NROOTS=(LIT+LJT+LKT+LLTT-2 + 1 )/2
!-----------------------------------------------------------------------
      RETURN
      END
! JKDSPDNOF
      SUBROUTINE JKDSPDNOF(TOTCOUNT,GINT,FINT,DAB,DCHRG,NIJGDIM,       &
                        MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL,    &
                        IIEQJJ,KKEQLL,IJEQKL,EXPNDI,EXPNDK,         &
                        SKIPI,SKIPJ,SKIPK,SKIPL,SPI,SPJ,SPK,SPL,    &
                        SPIJKL,IJKLG,MAXNUM,INVTYP,NROOTS,MINVEC,   &
                        MODTYP,MAXXYZ,NKL,NIJ,IIAT,JJAT,KKAT,LLAT,  &
                        GRADS,DABMAX,DABCUT,NKL0,NIJ0,              &
                        LIT,LJT,LKT,LLTT,IDER,JDER,KDER,LDER)          
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,INTENT(IN)::MAXNUM
      DOUBLE PRECISION,DIMENSION(MAXNUM),INTENT(IN)::DAB
      DOUBLE PRECISION,DIMENSION(15,NIJGDIM),INTENT(IN)::DCHRG
      DOUBLE PRECISION,INTENT(IN)::DABMAX,DABCUT
      LOGICAL,INTENT(IN)::IIEQJJ,KKEQLL,IJEQKL,EXPNDI,EXPNDK
      LOGICAL,INTENT(IN)::SKIPI,SKIPJ,SKIPK,SKIPL,SPI,SPJ,SPK,SPL,SPIJKL
      INTEGER,INTENT(IN)::MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL
      INTEGER,INTENT(IN)::INVTYP,NKL,NIJ,MAXXYZ,IIAT,JJAT,KKAT,LLAT
      INTEGER,INTENT(IN)::NKL0,NIJ0
      INTEGER,INTENT(IN)::LIT,LJT,LKT,LLTT,IDER,JDER,KDER,LDER
      INTEGER,INTENT(IN)::NROOTS,MINVEC
      INTEGER,INTENT(INOUT)::TOTCOUNT
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT)::GRADS
      DOUBLE PRECISION,DIMENSION(12)::FD
      DOUBLE PRECISION::XC,YC,ZC,DXIJ,DYIJ,DZIJ
      INTEGER,DIMENSION(4,MAXNUM),INTENT(IN)::IJKLG
!      
      LOGICAL::NMAXS,NMAXP,MMAXS,MMAXP
      INTEGER::NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX
!     ZERO AND FIRST DERIVATIVE INTEGRALS      
      DOUBLE PRECISION,DIMENSION(MAXNUM),INTENT(INOUT)::GINT,FINT
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::GNM,GNKL,GIJKL,FI,FJ,FK,FL,DIJ,DKL
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::AAI,AAJ,BBK,BBL,DIJSI,DKLSK,DKLSL,DIJSJ
      DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::C00,D00,F00,B00,B01,B10,RWV,ABV,CV
!
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00,PI252=34.986836655250D+00
      LOGICAL::FIRST,LAST
!-----------------------------------------------------------------------            
      ALLOCATE(GNM((MODTYP+MODTYP-1)*(MODTYP+MODTYP-1)*MAXXYZ*3))
      ALLOCATE(GNKL(MODTYP**2*(MODTYP+MODTYP-1)*MAXXYZ*3))
      ALLOCATE(GIJKL(MODTYP**2*MODTYP**2*MAXXYZ*3))
      ALLOCATE(FI(MODTYP**2*MODTYP**2*MAXXYZ*3),FJ(MODTYP**2*MODTYP**2*MAXXYZ*3))
      ALLOCATE(FK(MODTYP**2*MODTYP**2*MAXXYZ*3),FL(MODTYP**2*MODTYP**2*MAXXYZ*3))
      ALLOCATE(DIJ((MODTYP**2)*MAXXYZ*3),DKL((MODTYP+MODTYP-1)*MAXXYZ*3))
      ALLOCATE(AAI(MAXXYZ*3),AAJ(MAXXYZ*3),BBK(MAXXYZ*3),BBL(MAXXYZ*3))
      ALLOCATE(DIJSI(MAXXYZ),DKLSK(MAXXYZ),DKLSL(MAXXYZ),DIJSJ(MAXXYZ))
      ALLOCATE(C00(MAXXYZ,3),D00(MAXXYZ,3),F00(MAXXYZ,3),RWV(2,MAXXYZ),ABV(5,MAXXYZ),CV(18,MAXXYZ))
      ALLOCATE(B00((MODTYP+MODTYP-1)*MAXXYZ,3),B01((MODTYP+MODTYP-1)*MAXXYZ,3),B10((MODTYP+MODTYP-1)*MAXXYZ,3))      
!     
      NIMAX=LIT + IDER
      NJMAX=LJT + JDER
      NKMAX=LKT + KDER
      NLMAX=LLTT + LDER
      NMAX=LIT+LJT-1 + MIN0(IDER+JDER,1)
      MMAX=LKT+LLTT-1 + MIN0(KDER+LDER,1)
      NMAXS=NMAX.EQ.1
      NMAXP=NMAX.LE.2
      MMAXS=MMAX.EQ.1
      MMAXP=MMAX.LE.2
!
      DTOL=(10.0D+00)**(-20)
      DTOL=DTOL*DTOL
      Q4=PI252
      MAXG=MAXXYZ/NROOTS      
!
!     ----- PAIR OF K,L PRIMITIVES -----
!      
      FIRST=.TRUE.
      NG=0
      KLG=0
  100  KLG=KLG+1
       IF(KLG.GT.NKL) GO TO 300
       DB=DCHRG( 1,KLG+NKL0)
       BB=DCHRG( 2,KLG+NKL0)
       XB=DCHRG( 3,KLG+NKL0)
       YB=DCHRG( 4,KLG+NKL0)
       ZB=DCHRG( 5,KLG+NKL0)
       XD=DCHRG( 6,KLG+NKL0)
       YD=DCHRG( 7,KLG+NKL0)
       ZD=DCHRG( 8,KLG+NKL0)
       DXKL=DCHRG( 9,KLG+NKL0)
       DYKL=DCHRG(10,KLG+NKL0)
       DZKL=DCHRG(11,KLG+NKL0)
       Q4DB=Q4*DB
!
!     ----- PAIR OF I,J PRIMITIVES -----
!       
       IJG=0
  200  IJG=IJG+1
       IF(IJG.GT.NIJ) GO TO 100
       DA=DCHRG( 1,IJG+NIJ0)
       AA=DCHRG( 2,IJG+NIJ0)
       XA=DCHRG( 3,IJG+NIJ0)
       YA=DCHRG( 4,IJG+NIJ0)
       ZA=DCHRG( 5,IJG+NIJ0)
       Q4DBDA=Q4DB*DA
       AANDB1=(1.0D+00)/(AA+BB)        
       DUM=Q4DBDA*Q4DBDA*AANDB1
       IF(DUM.LE.DTOL) GO TO 200
       Q4DBDA=Q4DBDA*SQRT(AANDB1)
       IF(ABS(Q4DBDA*DABMAX).LT.DABCUT) GO TO 200
       RHO   =AA*BB*AANDB1
       XX=RHO*((XA-XB)**2+(YA-YB)**2+(ZA-ZB)**2)
!
       NG=NG+1
       ABV(1,NG)=AA
       ABV(2,NG)=BB
       ABV(3,NG)=RHO
       ABV(4,NG)=Q4DBDA
       ABV(5,NG)=XX
!
       XC=DCHRG( 6,IJG+NIJ0)
       YC=DCHRG( 7,IJG+NIJ0)
       ZC=DCHRG( 8,IJG+NIJ0)
       DXIJ=DCHRG( 9,IJG+NIJ0)
       DYIJ=DCHRG(10,IJG+NIJ0)
       DZIJ=DCHRG(11,IJG+NIJ0)
!
       AAI(NG)=DCHRG(12,IJG+NIJ0)
       AAJ(NG)=DCHRG(13,IJG+NIJ0)
       BBK(NG)=DCHRG(12,KLG+NKL0)
       BBL(NG)=DCHRG(13,KLG+NKL0)        
!
       IF(.NOT.MMAXS) THEN
        CV( 1,NG)=AA*(XA-XD)
        CV( 2,NG)=BB*(XB-XD)
        CV( 3,NG)=AA*(YA-YD)
        CV( 4,NG)=BB*(YB-YD)
        CV( 5,NG)=AA*(ZA-ZD)
        CV( 6,NG)=BB*(ZB-ZD)
       ENDIF
       IF(.NOT.NMAXS) THEN
        CV( 7,NG)=AA*(XA-XC)
        CV( 8,NG)=BB*(XB-XC)
        CV( 9,NG)=AA*(YA-YC)
        CV(10,NG)=BB*(YB-YC)
        CV(11,NG)=AA*(ZA-ZC)
        CV(12,NG)=BB*(ZB-ZC)
       ENDIF
       CV(13,NG)=DXIJ
       CV(14,NG)=DYIJ
       CV(15,NG)=DZIJ
       CV(16,NG)=DXKL
       CV(17,NG)=DYKL
       CV(18,NG)=DZKL
       IF(SPI) DIJSI(NG)=DCHRG(14,IJG+NIJ0)
       IF(SPJ) DIJSJ(NG)=DCHRG(15,IJG+NIJ0)
       IF(SPK) DKLSK(NG)=DCHRG(14,KLG+NKL0)
       IF(SPL) DKLSL(NG)=DCHRG(15,KLG+NKL0)        
!
       IF(NG.LT.MAXG) GO TO 200
       LAST=.FALSE.
       GO TO 310
!       
  300  CONTINUE
       LAST=.TRUE.
  310  CONTINUE
       NUMG=NG
       IF(NUMG.EQ.0) GO TO 998
!       
       IF(NROOTS.EQ.1) GO TO 480
       IF (SPI) THEN
           DO IROOT=2,NROOTS
           DO IG=1,NUMG
               DIJSI(IG+NUMG*(IROOT-1))=DIJSI(IG)
           ENDDO
           ENDDO
       ENDIF
       IF (SPJ) THEN
          DO IROOT=2,NROOTS
          DO IG=1,NUMG
             DIJSJ(IG+NUMG*(IROOT-1))=DIJSJ(IG)
          ENDDO
          ENDDO
       ENDIF
       IF (SPK) THEN
          DO IROOT=2,NROOTS
          DO IG=1,NUMG
             DKLSK(IG+NUMG*(IROOT-1))=DKLSK(IG)
          ENDDO
          ENDDO
       ENDIF
       IF (SPL) THEN
          DO IROOT=2,NROOTS
          DO IG=1,NUMG
             DKLSL(IG+NUMG*(IROOT-1))=DKLSL(IG)
          ENDDO
          ENDDO
       ENDIF      
!
  480 IF(.NOT.SKIPI) THEN
       DO IRXYZ=2,NROOTS*3
       DO IG=1,NUMG
       AAI(IG+NUMG*(IRXYZ-1))=AAI(IG)
       ENDDO
       ENDDO
      ENDIF
      IF(.NOT.SKIPJ) THEN
       DO IRXYZ=2,NROOTS*3
       DO IG=1,NUMG
       AAJ(IG+NUMG*(IRXYZ-1))=AAJ(IG)
       ENDDO
       ENDDO
      ENDIF
      IF(.NOT.SKIPK) THEN
       DO IRXYZ=2,NROOTS*3
       DO IG=1,NUMG
       BBK(IG+NUMG*(IRXYZ-1))=BBK(IG)
       ENDDO
       ENDDO
      ENDIF
      IF(.NOT.SKIPL) THEN
       DO IRXYZ=2,NROOTS*3
       DO IG=1,NUMG
       BBL(IG+NUMG*(IRXYZ-1))=BBL(IG)
       ENDDO
       ENDDO
      ENDIF
!
!     ----- COMPUTE ROOTS AND WEIGHTS FOR QUADRATURE -----
!
      CALL JKWRYSNOF(RWV,ABV,NUMG,NROOTS,NKL,NIJ)
!
!     ----- COMPUTE COEFFICIENTS FOR RECURSION FORMULAE -----
!
      CALL JKBCDFNOF(B00,B01,B10,C00,D00,F00,DIJ,DKL,  &
                  ABV,CV,RWV,NUMG,NROOTS,NKL,NIJ,   &
                  NMAXS,MMAXS)
!
!     ----- COMPUTE -X- , -Y- , -Z- INTEGRALS ( 2 CENTERS, 2-D ) -----
!
      IF(NUMG*NROOTS*3.LT.MINVEC) THEN
         CALL JKGNMSNOF(GNM,NUMG*NROOTS*3,NMAX,MMAX,   &
                     B00,B01,B10,C00,D00,F00,          &
                     NMAXS,NMAXP,MMAXS,MMAXP)
      ELSE
         WRITE(6,*) 'SCALAR MACHINES SHOULD NOT CALL -JKGNMV-'
         WRITE(6,*) 'ERROR IN SUBROUTINE JKDSPD'
         STOP
      END IF
!
!     ----- COMPUTE -X- , -Y- , -Z- INTEGRALS ( 4 CENTERS, 2-D ) -----
!
      IF(NUMG*NROOTS*3.LT.MINVEC) THEN
         CALL JKXYZSNOF(GIJKL,GIJKL,GNKL,GNKL,GNKL,GNM,      &
                    NUMG*NROOTS*3,NMAX,MMAX,NIMAX,NJMAX,     &
                    NKMAX,NLMAX,DIJ,DKL,EXPNDI,EXPNDK)
      ELSE
         WRITE(6,*) 'SCALAR MACHINES SHOULD NOT CALL -JKXYZV-'
         WRITE(6,*) 'ERROR IN SUBROUTINE JKDSPD'
         STOP
      END IF
!
!     ----- COMPUTE -X- , -Y- , -Z- INTEGRALS FOR DERIVATIVES -----
!
      IF(NUMG*NROOTS*3.LT.MINVEC) THEN
         CALL JDXYZSNOF(GIJKL,GIJKL,GIJKL,GIJKL,                       &
                     NUMG*NROOTS*3,NIMAX,NJMAX,NKMAX,NLMAX,            &
                     LIT,LJT,LKT,LLTT,AAI,AAJ,BBK,BBL,FI,FJ,FK,FL,     &
                     SKIPI,SKIPJ,SKIPK,SKIPL) 
      ELSE
         WRITE(6,*) 'SCALAR MACHINES SHOULD NOT CALL -JDXYZV-'
         WRITE(6,*) 'ERROR IN SUBROUTINE JKDSPD'
         STOP
      END IF      
!        
!       ZERO OUT FIRST TIME AROUND 
        IF(FIRST) THEN
         FD=ZERO
         FIRST=.FALSE.
        ENDIF
!
!     ----- COMPUTE DERIVATIVE INTEGRALS -----
!
      IF(NUMG*NROOTS.LT.MINVEC) THEN
         IGIJKL=MODTYP**2*MODTYP**2*MAXXYZ*3
         CALL DSPDFSNOF(TOTCOUNT,NUMG,NROOTS,IJKLG,GINT,FINT,GIJKL,   &
                     FI,FJ,FK,FL,DIJSI,DIJSJ,DKLSK,DKLSL,             &
                     DAB,SKIPI,SKIPJ,SKIPK,SKIPL,                     &
                     MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,         &
                     MAXNUM,MAXXYZ,IIEQJJ,KKEQLL,IJEQKL,              &
                     SPI,SPJ,SPK,SPL,SPIJKL,FD,IGIJKL)
      ELSE
         WRITE(6,*) 'SCALAR MACHINES SHOULD NOT CALL -DSPDFV-'
         WRITE(6,*) 'ERROR IN SUBROUTINE JKDSPD'
         STOP
      END IF
!
      IF(LAST) GO TO 998
      NG=0  
      GO TO 200
  998 CONTINUE
      IF(NUMG.EQ.0.AND.FIRST) RETURN
!      
!     ADD TO THE GRADIENT VECTOR ACCORDING TO SYMMETRY(INVTYP)
!
      CALL JKDINVNOF(INVTYP,FD,GRADS,IIAT,JJAT,KKAT,LLAT)
      DEALLOCATE(GNM,GNKL,GIJKL,FI,FJ,FK,FL,DIJ,DKL)
      DEALLOCATE(AAI,AAJ,BBK,BBL,DIJSI,DKLSK,DKLSL,DIJSJ)
      DEALLOCATE(C00,D00,F00,B00,B01,B10,RWV,ABV,CV)
!-----------------------------------------------------------------------
      RETURN
      END
! JKDINVNOF
      SUBROUTINE JKDINVNOF(INVTYP,FD,GRADS,II,JJ,KK,LL)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,INTENT(IN)::INVTYP,II,JJ,KK,LL
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT)::GRADS
      DOUBLE PRECISION,DIMENSION(3,4)::FD
!     ----- TRANSLATIONAL INVARIANCE FOR GRADIENT ELEMENTS -----
      IF (INVTYP.EQ.2) THEN
        DO IXYZ=1,3
          FD(IXYZ,1)=- FD(IXYZ,4)
        ENDDO
      ELSE IF (INVTYP.EQ.3) THEN
        DO IXYZ=1,3
          FD(IXYZ,1)=- FD(IXYZ,3)
        ENDDO
      ELSE IF (INVTYP.EQ.4.OR.INVTYP.EQ.5) THEN
        DO IXYZ=1,3
          FD(IXYZ,1)=-(FD(IXYZ,3)+FD(IXYZ,4))
        ENDDO
      ELSE IF (INVTYP.EQ.6) THEN
        DO IXYZ=1,3
          FD(IXYZ,1)=- FD(IXYZ,2)
        ENDDO
      ELSE IF (INVTYP.EQ.7) THEN
        DO IXYZ=1,3
          FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,4))
        ENDDO
      ELSE IF (INVTYP.EQ.8) THEN
        DO IXYZ=1,3
          FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,4))
        ENDDO
      ELSE IF (INVTYP.EQ.9.OR.INVTYP.EQ.10) THEN
        DO IXYZ=1,3
          FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,3))
        ENDDO
      ELSE IF (INVTYP.EQ.11) THEN
        DO IXYZ=1,3
          FD(IXYZ,2)=- FD(IXYZ,1)
        ENDDO
      ELSE IF (INVTYP.EQ.12) THEN
        DO IXYZ=1,3
          FD(IXYZ,2)=-(FD(IXYZ,1)+FD(IXYZ,4))
        ENDDO
      ELSE IF (INVTYP.EQ.13) THEN
        DO IXYZ=1,3
          FD(IXYZ,2)=-(FD(IXYZ,1)+FD(IXYZ,3))
        ENDDO
      ELSE IF (INVTYP.EQ.14) THEN
        DO IXYZ=1,3
          FD(IXYZ,3)=-(FD(IXYZ,1)+FD(IXYZ,2))
        ENDDO
      ELSE IF (INVTYP.EQ.15) THEN
        DO IXYZ=1,3
          FD(IXYZ,4)=-(FD(IXYZ,1)+FD(IXYZ,2)+FD(IXYZ,3))
        ENDDO
      ENDIF
      DO IXYZ=1,3
        GRADS(IXYZ,II)=GRADS(IXYZ,II)+FD(IXYZ,1)
        GRADS(IXYZ,JJ)=GRADS(IXYZ,JJ)+FD(IXYZ,2)
        GRADS(IXYZ,KK)=GRADS(IXYZ,KK)+FD(IXYZ,3)
        GRADS(IXYZ,LL)=GRADS(IXYZ,LL)+FD(IXYZ,4)
      ENDDO   
!-----------------------------------------------------------------------
      RETURN
      END
! DABCLUNOF
      SUBROUTINE DABCLUNOF(II,JJ,KK,LL,PNRM,                           &
                        KKMIN,KKMAX,KLOC,IGXYZ,JGXYZ,KGXYZ,LGXYZ,       &
                        IIEQJJ,KKEQLL,IJEQKL,IA,DA,DAB,MAXNUM,DABMAX)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,INTENT(IN)::II,JJ,KK,LL
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KLOC,KKMIN,KKMAX
      INTEGER,DIMENSION(4,35),INTENT(IN)::IGXYZ,JGXYZ,KGXYZ,LGXYZ
      INTEGER,DIMENSION(NBF),INTENT(IN)::IA
      LOGICAL,INTENT(IN)::IIEQJJ,KKEQLL,IJEQKL
      DOUBLE PRECISION,DIMENSION(NBFT),INTENT(IN)::DA  
      DOUBLE PRECISION,INTENT(OUT)::DABMAX
      DOUBLE PRECISION,DIMENSION(84),INTENT(IN)::PNRM
      DOUBLE PRECISION,DIMENSION(MAXNUM),INTENT(OUT)::DAB
      DOUBLE PRECISION,PARAMETER::ZER=0.0D+00,PT5=0.5D+00,F04=4.0D+00
!-----------------------------------------------------------------------
      DABMAX=ZER
!      
      MINI= KKMIN(II)
      MINJ= KKMIN(JJ)
      MINK= KKMIN(KK)
      MINL= KKMIN(LL)
      LOCI= KLOC(II)-MINI
      LOCJ= KLOC(JJ)-MINJ
      LOCK= KLOC(KK)-MINK
      LOCL= KLOC(LL)-MINL
      MAXI= KKMAX(II)
      MAXJ= KKMAX(JJ)
      MAXK= KKMAX(KK)
      MAXL= KKMAX(LL)
         DO I=MINI,MAXI
            P1I = PNRM(I)
            JMAX= MAXJ
            IF(IIEQJJ) JMAX= I
            DO J=MINJ,JMAX
               P2J = P1I*PNRM(J)
               IAJ= MAX0(LOCI+I,LOCJ+J)
               IIJ= MIN0(LOCI+I,LOCJ+J)
               KMMAX=MAXK
               IF(IJEQKL) KMMAX= I
               DO K=MINK,KMMAX
                  P3K = P2J*PNRM(K)
                  LMAX= MAXL
                  IF(KKEQLL) LMAX= K
                  IF(IJEQKL .AND. K.EQ.I) LMAX= J
                  DO L=MINL,LMAX
                     P4L= P3K*PNRM(L)
                     KAL= MAX0(LOCK+K,LOCL+L)
                     KIL= MIN0(LOCK+K,LOCL+L)
                     IN = IAJ
                     JN = IIJ
                     KN = KAL
                     LN = KIL
                     IF(IN.LT.KN .OR.(IN.EQ.KN .AND. JN.LT.LN)) THEN
                        IN = KAL
                        JN = KIL
                        KN = IAJ
                        LN = IIJ
                     ENDIF
                     IJ = IA(IN)+JN
                     IK = IA(IN)+KN
                     IL = IA(IN)+LN
                     JK = IA(MAX0(JN,KN))+MIN0(JN,KN)
                     JL = IA(JN)+LN
                     IF(JN.LT.KN) JL = IA(MAX0(JN,LN))+MIN0(JN,LN)
                     KL = IA(KN)+LN
! THE DENSITY IS BUILT ACCORDING TO THE PAPER OF DUPUIS (1978)                     
                     DF1= DA(IJ)*DA(KL)
                     DQ1= DA(IK)*DA(JL)+DA(IL)*DA(JK)
                     DF1= DF1*F04-DQ1
! AVOID DOUBLE COUNTING OF DIAGONAL TERMS NOT DONE IN SQUARETRIAN2                     
                     IF(JN.EQ.IN               ) DF1= DF1*PT5
                     IF(LN.EQ.KN               ) DF1= DF1*PT5
                     IF(KN.EQ.IN .AND. LN.EQ.JN) DF1= DF1*PT5
                     IF(DABMAX.LT. ABS(DF1)) DABMAX= ABS(DF1)
! IGXYZ AND J, K, AND L ARE SET UP IN JKDNDX
                     IJKL=IGXYZ(1,I)+JGXYZ(1,J)+KGXYZ(1,K)+LGXYZ(1,L)
                     DAB(IJKL)= DF1*P4L
                  ENDDO
               ENDDO
            ENDDO
         ENDDO   
!-----------------------------------------------------------------------
      RETURN
      END
! VNNDERNOF
      SUBROUTINE VNNDERNOF(CX0,CY0,CZ0,ZNUC,DE)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN)::CX0,CY0,CZ0,ZNUC
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT)::DE
      DOUBLE PRECISION,DIMENSION(3,NATOMS)::C
      DOUBLE PRECISION,DIMENSION(NATOMS,NATOMS)::DRG
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00,ONE=1.0D+00
!      
      C(1,:)=CX0
      C(2,:)=CY0
      C(3,:)=CZ0
!
!     ----- FORM DISTANCE MATRIX -----
!
      DRG(1,1) = ZERO
      DO K = 2,NATOMS
         DRG(K,K) = ZERO
         K1 = K-1
         DO L = 1,K1
            RKL = ZERO
            DO I = 1,3
               RKL = RKL+(C(I,K)-C(I,L))**2
            ENDDO
            DRG(K,L) = -ONE/RKL
            DRG(L,K) = SQRT(RKL)
         ENDDO
      ENDDO
!
!     ----- NUCLEAR REPULSION CONTRIBUTION TO GRADIENT -----
!
      DO KK = 1,3
         DO K = 2,NATOMS
            ZAK = ZNUC(K)
            KM1 = K-1
            DO L = 1,KM1
               ZAL = ZNUC(L)
               PKL = (C(KK,K)-C(KK,L))/DRG(L,K)
               DE(KK,K) = DE(KK,K)+PKL*DRG(K,L)*ZAK*ZAL
            ENDDO
         ENDDO
!
         NAT1 = NATOMS-1
         DO K = 1,NAT1
            ZAK = ZNUC(K)
            KP1 = K+1
            DO L = KP1,NATOMS
               ZAL = ZNUC(L)
               PKL = (C(KK,K)-C(KK,L))/DRG(K,L)
               DE(KK,K) = DE(KK,K)+PKL*DRG(L,K)*ZAK*ZAL
            ENDDO
         ENDDO
      ENDDO
!-----------------------------------------------------------------------  
      RETURN
      END
! SDERNOF
      SUBROUTINE SDERNOF(KATOM,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG,      &
                      CX0,CY0,CZ0,EX1,CS,CP,CD,CF,CG,EPS,DE)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KATOM,KTYPE,KLOC,KKMIN,KKMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KSTART,KNG
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN)::CX0,CY0,CZ0
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN)::EX1,CS,CP,CD,CF,CG
      DOUBLE PRECISION,DIMENSION(NBFT),INTENT(IN)::EPS
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT)::DE
!      
      INTEGER,DIMENSION(NBF)::IA
      INTEGER,DIMENSION(35)::IJX,IJY,IJZ
      DOUBLE PRECISION,DIMENSION(5,5)::DXS,DYS,DZS
      DOUBLE PRECISION,DIMENSION(6,5)::XS,YS,ZS
      DOUBLE PRECISION,DIMENSION(225)::DIJ
      DOUBLE PRECISION::TOL
!
      DOUBLE PRECISION,PARAMETER::SQRT3=1.73205080756888D+00
      DOUBLE PRECISION,PARAMETER::SQRT5=2.23606797749979D+00
      DOUBLE PRECISION,PARAMETER::SQRT7=2.64575131106459D+00
      DOUBLE PRECISION,PARAMETER::ONE=1.0D+00,RLN10=2.30258D+00
!
      DATA IJX / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,           &
                 4, 1, 1, 3, 3, 2, 1, 2, 1, 2,           &
                 5, 1, 1, 4, 4, 2, 1, 2, 1, 3,           &
                 3, 1, 3, 2, 2/
      DATA IJY / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,           &
                 1, 4, 1, 2, 1, 3, 3, 1, 2, 2,           &
                 1, 5, 1, 2, 1, 4, 4, 1, 2, 3,           &
                 1, 3, 2, 3, 2/
      DATA IJZ / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,           &
                 1, 1, 4, 1, 2, 1, 2, 3, 3, 2,           &
                 1, 1, 5, 1, 2, 1, 2, 4, 4, 1,           &
                 3, 3, 2, 2, 3/      
      DO I=1,NBF
        IA(I) = (I*I-I)/2
      ENDDO   
      TOL = RLN10*20
!
!     ----- I SHELL
!
      DO II = 1,NSHELL
!
      IAT = KATOM(II)
      XI = CX0(IAT)
      YI = CY0(IAT)
      ZI = CZ0(IAT)
      I1 = KSTART(II)
      I2 = I1+KNG(II)-1
      LIT = KTYPE(II)
      MINI = KKMIN(II)
      MAXI = KKMAX(II)
      LOCI = KLOC(II)-MINI
      LITDER = LIT+1     
!
!     ----- J SHELL
!
      DO JJ = 1,II
!
      JAT = KATOM(JJ)
      XJ = CX0(JAT)
      YJ = CY0(JAT)
      ZJ = CZ0(JAT)
      J1 = KSTART(JJ)
      J2 = J1+KNG(JJ)-1
      LJT = KTYPE(JJ)
      MINJ = KKMIN(JJ)
      MAXJ = KKMAX(JJ)
      LOCJ = KLOC(JJ)-MINJ
      RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      IF(II.EQ.JJ) CYCLE
!
!     ----- I PRIMITIVE
!
      DO IG = I1,I2
        AI = EX1(IG)
        ARRI = AI*RR
        AXI = AI*XI
        AYI = AI*YI
        AZI = AI*ZI
        CSI=CS(IG)
        CPI=CP(IG)
        CDI=CD(IG)
        CFI=CF(IG)
        CGI=CG(IG)
!
!     ----- J PRIMITIVE
!
        DO JG = J1,J2
          AJ = EX1(JG)
          AA = AI+AJ
          AA1 = ONE/AA
          DUM = AJ*ARRI*AA1
          IF(DUM .GT. TOL) CYCLE
          FAC = EXP(-DUM)
          CSJ = CS(JG)
          CPJ = CP(JG)
          CDJ = CD(JG)
          CFJ = CF(JG)
          CGJ = CG(JG)
          AX = (AXI+AJ*XJ)*AA1
          AY = (AYI+AJ*YJ)*AA1
          AZ = (AZI+AJ*ZJ)*AA1      
!
!     ----- DENSITY FACTOR
!
          IJ = 0
          DO I = MINI,MAXI
            IF(I.EQ.1) DUM1=CSI*FAC
            IF(I.EQ.2) DUM1=CPI*FAC
            IF(I.EQ.5) DUM1=CDI*FAC
            IF(I.EQ.8) DUM1=DUM1*SQRT3
            IF(I.EQ.11) DUM1=CFI*FAC
            IF(I.EQ.14) DUM1=DUM1*SQRT5
            IF(I.EQ.20) DUM1=DUM1*SQRT3
            IF(I.EQ.21) DUM1=CGI*FAC
            IF(I.EQ.24) DUM1=DUM1*SQRT7
            IF(I.EQ.30) DUM1=DUM1*SQRT5/SQRT3
            IF(I.EQ.33) DUM1=DUM1*SQRT3
            DO J = MINJ,MAXJ
              IF(J.EQ.1) DUM2=DUM1*CSJ
              IF(J.EQ.2) DUM2=DUM1*CPJ
              IF(J.EQ.5) DUM2=DUM1*CDJ
              IF(J.EQ.8) DUM2=DUM2*SQRT3
              IF(J.EQ.11) DUM2=DUM1*CFJ
              IF(J.EQ.14) DUM2=DUM2*SQRT5
              IF(J.EQ.20) DUM2=DUM2*SQRT3
              IF(J.EQ.21) DUM2=DUM1*CGJ
              IF(J.EQ.24) DUM2=DUM2*SQRT7
              IF(J.EQ.30) DUM2=DUM2*SQRT5/SQRT3
              IF(J.EQ.33) DUM2=DUM2*SQRT3
              IJ=IJ+1
              NN=IA(LOCI+I)+(LOCJ+J)
              DEN = EPS(NN)
              DIJ(IJ) = DUM2*DEN
            ENDDO
          ENDDO
!
!     ----- OVERLAP
!
      T = SQRT(AA1)
      X0 = AX
      Y0 = AY
      Z0 = AZ
      DO J = 1,LJT
        NJ = J
        DO I = 1,LITDER
          NI = I
          CALL VINTNOF(NI,NJ,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ &
                        ,XINT,YINT,ZINT)
          XS(I,J)=XINT*T
          YS(I,J)=YINT*T
          ZS(I,J)=ZINT*T
        ENDDO
      ENDDO
!
      CALL DERINOF(DXS,DYS,DZS,XS,YS,ZS,LIT,LJT,AI)
!
      IJ=0
      DO I=MINI,MAXI
        IX=IJX(I)
        IY=IJY(I)
        IZ=IJZ(I)
        DO J=MINJ,MAXJ
          JX=IJX(J)
          JY=IJY(J)
          JZ=IJZ(J)
          DUMX=DXS(IX,JX)* YS(IY,JY)* ZS(IZ,JZ)
          DUMY= XS(IX,JX)*DYS(IY,JY)* ZS(IZ,JZ)
          DUMZ= XS(IX,JX)* YS(IY,JY)*DZS(IZ,JZ)
          IJ=IJ+1
          DE(1,IAT)=DE(1,IAT)+(DUMX*DIJ(IJ))
          DE(2,IAT)=DE(2,IAT)+(DUMY*DIJ(IJ))
          DE(3,IAT)=DE(3,IAT)+(DUMZ*DIJ(IJ))
          DE(1,JAT)=DE(1,JAT)-(DUMX*DIJ(IJ))
          DE(2,JAT)=DE(2,JAT)-(DUMY*DIJ(IJ))
          DE(3,JAT)=DE(3,JAT)-(DUMZ*DIJ(IJ))
        ENDDO
      ENDDO
!
      ENDDO
      ENDDO
!
!     ----- END OF PRIMITIVE LOOPS -----
!
      ENDDO
      ENDDO
!
!     ----- END OF SHELL LOOPS -----
!  
!-----------------------------------------------------------------------      
      RETURN
      END
! HELFEYNOF
      SUBROUTINE HELFEYNOF(KATOM,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG,    &
                      CX0,CY0,CZ0,EX1,CS,CP,CD,CF,CG,ZAN,PM,DE)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KATOM,KTYPE,KLOC,KKMIN,KKMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KSTART,KNG
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN)::CX0,CY0,CZ0,ZAN
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN)::EX1,CS,CP,CD,CF,CG
      DOUBLE PRECISION,DIMENSION(NBF,NBF),INTENT(IN)::PM
      DOUBLE PRECISION,DIMENSION(NBFT)::P
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT)::DE
!
      LOGICAL::IANDJ,DOUBLETE
      INTEGER::NROOTS
      INTEGER,DIMENSION(NBF)::IA
      INTEGER,DIMENSION(35)::IJX,IJY,IJZ
      DOUBLE PRECISION,DIMENSION(225)::DIJ
      DOUBLE PRECISION,DIMENSION(5,5,5,2)::XIN,YIN,ZIN
      DOUBLE PRECISION,DIMENSION(13)::UROOT,WROOT
      DOUBLE PRECISION::TOL,XX,ZNUCC      
!
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00, ONE=1.0D+00
      DOUBLE PRECISION,PARAMETER::RLN10=2.30258D+00,PI212=1.1283791670955D+00
      DOUBLE PRECISION,PARAMETER::SQRT3=1.73205080756888D+00
      DOUBLE PRECISION,PARAMETER::SQRT5=2.23606797749979D+00
      DOUBLE PRECISION,PARAMETER::SQRT7=2.64575131106459D+00
!
!  THE IJX, IJY, AND IJZ ARRAYS CONTAIN THE POWERS OF THE CARTESIAN
!  GAUSSIANS PLUS 1 IN EVERY PLACE.
!
      DATA IJX/ 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,   &
                4, 1, 1, 3, 3, 2, 1, 2, 1, 2,   &
                5, 1, 1, 4, 4, 2, 1, 2, 1, 3,   &
                3, 1, 3, 2, 2/
      DATA IJY/ 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,   &
                1, 4, 1, 2, 1, 3, 3, 1, 2, 2,   &
                1, 5, 1, 2, 1, 4, 4, 1, 2, 3,   &
                1, 3, 2, 3, 2/
      DATA IJZ/ 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,   &
                1, 1, 4, 1, 2, 1, 2, 3, 3, 2,   &
                1, 1, 5, 1, 2, 1, 2, 4, 4, 1,   &
                3, 3, 2, 2, 3/      
!
!     ----- HELMANN-FEYNMAN GRADIENT TERM -----
!     INTEGRAL TYPE IS <II/H'/JJ> = <II/V'/JJ>
!                     
      DO I=1,NBF
        IA(I) = (I*I-I)/2
      ENDDO   
      TOL = RLN10*20
      CALL SQUARETRIAN2(PM,P,NBF,NBFT)
      P=P*0.5D+00
!
!     ----- I SHELL
!
      DO II = 1,NSHELL
!
      I = KATOM(II)
      XI = CX0(I)
      YI = CY0(I)
      ZI = CZ0(I)
      I1 = KSTART(II)
      I2 = I1+KNG(II)-1
      LIT = KTYPE(II)
      MINI = KKMIN(II)
      MAXI = KKMAX(II)
      LOCI = KLOC(II)-MINI
!
!     ----- J SHELL
!
      DO JJ = 1,II
!
        J = KATOM(JJ)
        XJ = CX0(J)
        YJ = CY0(J)
        ZJ = CZ0(J)
        J1 = KSTART(JJ)
        J2 = J1+KNG(JJ)-1
        LJT = KTYPE(JJ)
        MINJ = KKMIN(JJ)
        MAXJ = KKMAX(JJ)
        LOCJ = KLOC(JJ)-MINJ
        NROOTS = (LIT+LJT+1-2)/2 + 1
        RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
        IANDJ = II .EQ. JJ
!
!     ----- I PRIMITIVE
!
        DO IG = I1,I2
        AI = EX1(IG)
        ARRI = AI*RR
        AXI = AI*XI
        AYI = AI*YI
        AZI = AI*ZI
        CSI = CS(IG)
        CPI = CP(IG)
        CDI = CD(IG)
        CFI = CF(IG)
        CGI = CG(IG)
!
!     ----- J PRIMITIVE
!
        JGMAX = J2
        IF(IANDJ) JGMAX = IG
        DO JG = J1,JGMAX
          AJ = EX1(JG)
          AA = AI+AJ
          AA1 = ONE/AA
          DUM = AJ*ARRI*AA1
          IF(DUM .GT. TOL) CYCLE
          FAC = EXP(-DUM)
          CSJ = CS(JG)
          CPJ = CP(JG)
          CDJ = CD(JG)
          CFJ = CF(JG)
          CGJ = CG(JG)
          AX = (AXI+AJ*XJ)*AA1
          AY = (AYI+AJ*YJ)*AA1
          AZ = (AZI+AJ*ZJ)*AA1       
!
!     ----- DENSITY FACTOR
!
          DOUBLETE=IANDJ.AND.IG.NE.JG
          JMAX = MAXJ
          NN = 0
          DUM1 = ZERO
          DUM2 = DUM1
          DO I = MINI,MAXI
            IF(I.EQ.1) DUM1=CSI*FAC
            IF(I.EQ.2) DUM1=CPI*FAC
            IF(I.EQ.5) DUM1=CDI*FAC
            IF(I.EQ.8) DUM1=DUM1*SQRT3
            IF(I.EQ.11) DUM1=CFI*FAC
            IF(I.EQ.14) DUM1=DUM1*SQRT5
            IF(I.EQ.20) DUM1=DUM1*SQRT3
            IF(I.EQ.21) DUM1=CGI*FAC
            IF(I.EQ.24) DUM1=DUM1*SQRT7
            IF(I.EQ.30) DUM1=DUM1*SQRT5/SQRT3
            IF(I.EQ.33) DUM1=DUM1*SQRT3
!
            IF(IANDJ) JMAX = I
            DO J = MINJ,JMAX
              IF(J.EQ.1) THEN
                DUM2=DUM1*CSJ
                IF( .NOT. DOUBLETE) GO TO 350
                IF(I .GT. 1) THEN
                 DUM2 = DUM2+CSI*CPJ*FAC
                ELSE
                 DUM2 = DUM2+DUM2
                END IF
              ELSE IF(J.EQ.2) THEN
                DUM2=DUM1*CPJ
                IF(DOUBLETE) DUM2 = DUM2+DUM2
              ELSE IF(J.EQ.5) THEN
                DUM2=DUM1*CDJ
                IF(DOUBLETE) DUM2 = DUM2+DUM2
              ELSE IF(J.EQ.8) THEN
                DUM2 = DUM2*SQRT3
              ELSE IF(J.EQ.11) THEN
                DUM2=DUM1*CFJ
                IF(DOUBLETE) DUM2=DUM2+DUM2
              ELSE IF(J.EQ.14) THEN
                DUM2=DUM2*SQRT5
              ELSE IF(J.EQ.20) THEN
                DUM2=DUM2*SQRT3
              ELSE IF(J.EQ.21) THEN
                DUM2=DUM1*CGJ
                IF(DOUBLETE) DUM2=DUM2+DUM2
              ELSE IF(J.EQ.24) THEN
                DUM2=DUM2*SQRT7
              ELSE IF(J.EQ.30) THEN
                DUM2=DUM2*SQRT5/SQRT3
              ELSE IF(J.EQ.33) THEN
                DUM2=DUM2*SQRT3
              END IF
!
  350         NN = NN+1
!  
              NDUM = IA(LOCI+I)+(LOCJ+J)
              DEN = P(NDUM)
              IF(IANDJ.AND.I.EQ.J) DEN=DEN*0.5D+00
              DIJ(NN)=DUM2*DEN*PI212*AA1  
            ENDDO
          ENDDO
!
!     ..... HELLMANN-FEYNMAN TERM .....
!
          AAX = AA*AX
          AAY = AA*AY
          AAZ = AA*AZ
          DO IC = 1,NATOMS
            ZNUCC = -ZAN(IC)
            CX = CX0(IC)
            CY = CY0(IC)
            CZ = CZ0(IC)
            XX = AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
            IF(NROOTS.LE.3) CALL RT123NOF(XX,NROOTS,UROOT,WROOT)
            IF(NROOTS.EQ.4) CALL ROOT4NOF(XX,UROOT,WROOT)
            IF(NROOTS.EQ.5) CALL ROOT5NOF(XX,UROOT,WROOT)
            DO K = 1,NROOTS
              UU = AA*UROOT(K)
              WW = WROOT(K)*ZNUCC
              WW=WW*(UU+UU)
              TT = ONE/(AA+UU)
              T = SQRT(TT)
              X0 = (AAX+UU*CX)*TT
              Y0 = (AAY+UU*CY)*TT
              Z0 = (AAZ+UU*CZ)*TT
              DO J = 1,LJT
                NJ = J
                DO I = 1,LIT
                  NI = I
                  CALL VINTNOF(NI,NJ,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ &
                            ,XINT,YINT,ZINT)
                  XIN(I,J,K,1) = XINT
                  YIN(I,J,K,1) = YINT
                  ZIN(I,J,K,1) = ZINT*WW
                  CALL DVINTNOF(NI,NJ,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ &
                             ,CX,CY,CZ,XINT,YINT,ZINT)
                  XIN(I,J,K,2) = XINT
                  YIN(I,J,K,2) = YINT
                  ZIN(I,J,K,2) = ZINT*WW
                ENDDO
              ENDDO  
            ENDDO
           IJ=0
           DO I=MINI,MAXI
            IX=IJX(I)
            IY=IJY(I)
            IZ=IJZ(I)
            JMAX=MAXJ
            IF(IANDJ) JMAX=I
            DO J=MINJ,JMAX
             JX=IJX(J)
             JY=IJY(J)
             JZ=IJZ(J)
             DUMX = ZERO
             DUMY = ZERO
             DUMZ = ZERO
             DO K = 1,NROOTS
              DUMX = DUMX+XIN(IX,JX,K,2)*YIN(IY,JY,K,1)*ZIN(IZ,JZ,K,1)
              DUMY = DUMY+XIN(IX,JX,K,1)*YIN(IY,JY,K,2)*ZIN(IZ,JZ,K,1)
              DUMZ = DUMZ+XIN(IX,JX,K,1)*YIN(IY,JY,K,1)*ZIN(IZ,JZ,K,2)
             ENDDO
             IJ=IJ+1  
             DUM=DIJ(IJ)
             DE(1,IC) = DE(1,IC)+DUM*DUMX
             DE(2,IC) = DE(2,IC)+DUM*DUMY
             DE(3,IC) = DE(3,IC)+DUM*DUMZ
            ENDDO
           ENDDO
          ENDDO
      ENDDO
      ENDDO
!
!     ----- END OF *PRIMITIVE* LOOPS -----
!
      ENDDO
      ENDDO
!
!     ----- END OF *SHELL* LOOPS -----
!             
!-----------------------------------------------------------------------
      RETURN
      END      
! TVDERNOF
      SUBROUTINE TVDERNOF(KATOM,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG,     &
                      CX0,CY0,CZ0,EX1,CS,CP,CD,CF,CG,ZAN,PM,DE)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KATOM,KTYPE,KLOC,KKMIN,KKMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KSTART,KNG
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN)::CX0,CY0,CZ0,ZAN
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN)::EX1,CS,CP,CD,CF,CG
      DOUBLE PRECISION,DIMENSION(NBF,NBF),INTENT(IN)::PM
      DOUBLE PRECISION,DIMENSION(NBFT)::P
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT)::DE
!      
      INTEGER::NROOTS
      INTEGER,DIMENSION(NBF)::IA
      INTEGER,DIMENSION(35)::IJX,IJY,IJZ
      DOUBLE PRECISION,DIMENSION(5,5)::DXS,DYS,DZS,DXT,DYT,DZT
      DOUBLE PRECISION,DIMENSION(6,7)::XS,YS,ZS
      DOUBLE PRECISION,DIMENSION(6,5)::XT,YT,ZT
      DOUBLE PRECISION,DIMENSION(6,5,5)::XV,YV,ZV
      DOUBLE PRECISION,DIMENSION(5,5,5)::DXV,DYV,DZV
      DOUBLE PRECISION,DIMENSION(225)::DIJ
      DOUBLE PRECISION,DIMENSION(13)::UROOT,WROOT
      DOUBLE PRECISION::TOL
      DOUBLE PRECISION::ZNUCC
      DOUBLE PRECISION::XX
!
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00, ONE=1.0D+00,TWO=2.0D+00
      DOUBLE PRECISION,PARAMETER::RLN10=2.30258D+00,PI212=1.1283791670955D+00
      DOUBLE PRECISION,PARAMETER::SQRT3=1.73205080756888D+00
      DOUBLE PRECISION,PARAMETER::SQRT5=2.23606797749979D+00
      DOUBLE PRECISION,PARAMETER::SQRT7=2.64575131106459D+00
!
      DATA IJX / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,   &
                 4, 1, 1, 3, 3, 2, 1, 2, 1, 2,   &
                 5, 1, 1, 4, 4, 2, 1, 2, 1, 3,   &
                 3, 1, 3, 2, 2/
      DATA IJY / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,   &
                 1, 4, 1, 2, 1, 3, 3, 1, 2, 2,   &
                 1, 5, 1, 2, 1, 4, 4, 1, 2, 3,   &
                 1, 3, 2, 3, 2/
      DATA IJZ / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,   &
                 1, 1, 4, 1, 2, 1, 2, 3, 3, 2,   &
                 1, 1, 5, 1, 2, 1, 2, 4, 4, 1,   &
                 3, 3, 2, 2, 3/
!-----------------------------------------------------------------------
!
!     ----- BASIS FUNCTION DERIVATIVE CONTRIBUTIONS TO GRADIENT -----
!     INTEGRALS ARE OF TYPE <II'/H/JJ> = <II'/T+V/JJ>
!
      DO I=1,NBF
        IA(I) = (I*I-I)/2
      ENDDO
      IAZ=0
      TOL = RLN10*20
      CALL SQUARETRIAN2(PM,P,NBF,NBFT)
      P=P*0.5D+00      
!
!     ----- I SHELL
!
      DO II = 1,NSHELL
!      
        IAT = KATOM(II)
        XI = CX0(IAT)
        YI = CY0(IAT)
        ZI = CZ0(IAT)
        I1 = KSTART(II)
        I2 = I1+KNG(II)-1
        LIT = KTYPE(II)
        MINI = KKMIN(II)
        MAXI = KKMAX(II)
        LOCI = KLOC(II)-MINI
        LITDER = LIT + 1      
!
!     ----- J SHELL
!
        DO JJ = 1,NSHELL
!
          JAT = KATOM(JJ)
          XJ = CX0(JAT)
          YJ = CY0(JAT)
          ZJ = CZ0(JAT)
          J1 = KSTART(JJ)
          J2 = J1+KNG(JJ)-1
          LJT = KTYPE(JJ)
          MINJ = KKMIN(JJ)
          MAXJ = KKMAX(JJ)
          LOCJ = KLOC(JJ)-MINJ
          LJTMOD = LJT+2
          NROOTS = (LIT+LJT-1)/2 + 1
          RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2       
!
!     ----- I PRIMITIVE
!
        DO IG = I1,I2
          AI = EX1(IG)
          ARRI = AI*RR
          AXI = AI*XI
          AYI = AI*YI
          AZI = AI*ZI
          CSI=CS(IG)
          CPI=CP(IG)
          CDI=CD(IG)
          CFI=CF(IG)
          CGI=CG(IG)   
!
!     ----- J PRIMITIVE
!
          DO JG = J1,J2
            AJ = EX1(JG)
            AA = AI+AJ
            AA1 = ONE/AA
            DUM = AJ*ARRI*AA1
            IF(DUM .GT. TOL) CYCLE
            FAC = EXP(-DUM)
            CSJ = CS(JG)
            CPJ = CP(JG)
            CDJ = CD(JG)
            CFJ = CF(JG)
            CGJ = CG(JG)
            AX = (AXI+AJ*XJ)*AA1
            AY = (AYI+AJ*YJ)*AA1
            AZ = (AZI+AJ*ZJ)*AA1
!
!     ----- DENSITY FACTOR
!
            IJ = 0
            DUM1 = ZERO
            DUM2 = DUM1
            DO I=MINI,MAXI
              IF(I.EQ.1) DUM1=CSI*FAC
              IF(I.EQ.2) DUM1=CPI*FAC
              IF(I.EQ.5) DUM1=CDI*FAC
              IF(I.EQ.8) DUM1=DUM1*SQRT3
              IF(I.EQ.11) DUM1=CFI*FAC
              IF(I.EQ.14) DUM1=DUM1*SQRT5
              IF(I.EQ.20) DUM1=DUM1*SQRT3
              IF(I.EQ.21) DUM1=CGI*FAC
              IF(I.EQ.24) DUM1=DUM1*SQRT7
              IF(I.EQ.30) DUM1=DUM1*SQRT5/SQRT3
              IF(I.EQ.33) DUM1=DUM1*SQRT3
!
              DO J = MINJ,MAXJ
                IF(J.EQ.1) DUM2=DUM1*CSJ
                IF(J.EQ.2) DUM2=DUM1*CPJ
                IF(J.EQ.5) DUM2=DUM1*CDJ
                IF(J.EQ.8) DUM2=DUM2*SQRT3
                IF(J.EQ.11) DUM2=DUM1*CFJ
                IF(J.EQ.14) DUM2=DUM2*SQRT5
                IF(J.EQ.20) DUM2=DUM2*SQRT3
                IF(J.EQ.21) DUM2=DUM1*CGJ
                IF(J.EQ.24) DUM2=DUM2*SQRT7
                IF(J.EQ.30) DUM2=DUM2*SQRT5/SQRT3
                IF(J.EQ.33) DUM2=DUM2*SQRT3
!
                IJ=IJ+1            
                NN = IA(MAX0(LOCI+I,LOCJ+J))+MIN0(LOCI+I,LOCJ+J)
                DEN = P(NN)
                DIJ(IJ)=DUM2*DEN
              ENDDO
            ENDDO
!     
!     -----  KINETIC ENERGY
!
      T = SQRT(AA1)
      X0 = AX
      Y0 = AY
      Z0 = AZ
      DO J = 1,LJTMOD
        NJ =J
        DO I = 1,LITDER
          NI = I
          CALL VINTNOF(NI,NJ,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ   &
                        ,XINT,YINT,ZINT)
          XS(I,J)=XINT*T
          YS(I,J)=YINT*T
          ZS(I,J)=ZINT*T
        ENDDO
      ENDDO
      CALL DTXYZNOF(XT,YT,ZT,XS,YS,ZS,LITDER,LJT,AJ)
      CALL DERINOF(DXS,DYS,DZS,XS,YS,ZS,LIT,LJT,AI)
      CALL DERINOF(DXT,DYT,DZT,XT,YT,ZT,LIT,LJT,AI)
      IJ=0
      DO I=MINI,MAXI
        IX=IJX(I)
        IY=IJY(I)
        IZ=IJZ(I)
        DO J=MINJ,MAXJ
          JX=IJX(J)
          JY=IJY(J)
          JZ=IJZ(J)
          DUMX=DXT(IX,JX)* YS(IY,JY)* ZS(IZ,JZ)    &
              +DXS(IX,JX)* YT(IY,JY)* ZS(IZ,JZ)    &
              +DXS(IX,JX)* YS(IY,JY)* ZT(IZ,JZ)
          DUMY= XT(IX,JX)*DYS(IY,JY)* ZS(IZ,JZ)    &
              + XS(IX,JX)*DYT(IY,JY)* ZS(IZ,JZ)    &
              + XS(IX,JX)*DYS(IY,JY)* ZT(IZ,JZ)
          DUMZ= XT(IX,JX)* YS(IY,JY)*DZS(IZ,JZ)    &
              + XS(IX,JX)* YT(IY,JY)*DZS(IZ,JZ)    &
              + XS(IX,JX)* YS(IY,JY)*DZT(IZ,JZ)
          IJ=IJ+1
          DE(1,IAT)=DE(1,IAT)+ DUMX*DIJ(IJ)
          DE(2,IAT)=DE(2,IAT)+ DUMY*DIJ(IJ)
          DE(3,IAT)=DE(3,IAT)+ DUMZ*DIJ(IJ)          
        ENDDO
      ENDDO
!
!     ..... NUCLEAR ATTRACTION
!
      AAX = AA*AX
      AAY = AA*AY
      AAZ = AA*AZ
!
      DO IC = 1,NATOMS
            IF(IC.LE.NATOMS) THEN
               ZNUCC = -ZAN(IC)
               CX = CX0(IC)
               CY = CY0(IC)
               CZ = CZ0(IC)
            ENDIF
         XX = AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
         IF(NROOTS.LE.3) CALL RT123NOF(XX,NROOTS,UROOT,WROOT)
         IF(NROOTS.EQ.4) CALL ROOT4NOF(XX,UROOT,WROOT)
         IF(NROOTS.EQ.5) CALL ROOT5NOF(XX,UROOT,WROOT)
         DO K = 1,NROOTS
            UU = AA*UROOT(K)
            WW = WROOT(K)*ZNUCC
            TT = ONE/(AA+UU)
            T = SQRT(TT)
            X0 = (AAX+UU*CX)*TT
            Y0 = (AAY+UU*CY)*TT
            Z0 = (AAZ+UU*CZ)*TT
            DO J = 1,LJT
               NJ = J
               DO I = 1,LITDER
                  NI = I
                  CALL VINTNOF(NI,NJ,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ &
                        ,XINT,YINT,ZINT)
                  XV(I,J,K) = XINT
                  YV(I,J,K) = YINT
                  ZV(I,J,K) = ZINT*WW
               ENDDO
            ENDDO
            CALL DERINOF(DXV(1,1,K),DYV(1,1,K),DZV(1,1,K),    &
                      XV(1,1,K), YV(1,1,K), ZV(1,1,K),LIT,LJT,AI)
         ENDDO
         IJ=0
         DO I=MINI,MAXI
           IX=IJX(I)
           IY=IJY(I)
           IZ=IJZ(I)
           DO J=MINJ,MAXJ
             JX=IJX(J)
             JY=IJY(J)
             JZ=IJZ(J)
             DUMX=ZERO
             DUMY=ZERO
             DUMZ=ZERO
             DO K=1,NROOTS
               DUMX=DUMX+DXV(IX,JX,K)* YV(IY,JY,K)* ZV(IZ,JZ,K)
               DUMY=DUMY+ XV(IX,JX,K)*DYV(IY,JY,K)* ZV(IZ,JZ,K)
               DUMZ=DUMZ+ XV(IX,JX,K)* YV(IY,JY,K)*DZV(IZ,JZ,K)
             ENDDO
             IJ=IJ+1
             IF((IC.GT.NATOMS).AND.(IAT.EQ.IAZ)) CYCLE
             DUMINT=DIJ(IJ)*AA1*PI212
             DE(1,IAT)=DE(1,IAT)+DUMX*DUMINT
             DE(2,IAT)=DE(2,IAT)+DUMY*DUMINT
             DE(3,IAT)=DE(3,IAT)+DUMZ*DUMINT
           ENDDO
         ENDDO
      ENDDO
!
      ENDDO
      ENDDO
!
!     ----- END OF PRIMITIVE LOOPS -----
!
      ENDDO
      ENDDO
!
!     ----- END OF SHELL LOOPS -----
!
!-----------------------------------------------------------------------
      RETURN
      END
! JKWRYSNOF
      SUBROUTINE JKWRYSNOF(RWV,ABV,NUMG,NROOTS,NKL,NIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      INTEGER,INTENT(IN)::NUMG,NROOTS,NKL,NIJ
      DOUBLE PRECISION,DIMENSION(5,NKL*NIJ),INTENT(IN)::ABV
      DOUBLE PRECISION,DIMENSION(2,NUMG,NROOTS),INTENT(OUT)::RWV
      DOUBLE PRECISION,DIMENSION(13)::UROOT,WROOT
!
      DO NG=1,NUMG
         XX=ABV(5,NG)
         IF(NROOTS.LE.3) CALL RT123NOF(XX,NROOTS,UROOT,WROOT)
         IF(NROOTS.EQ.4) CALL ROOT4NOF(XX,UROOT,WROOT)
         IF(NROOTS.EQ.5) CALL ROOT5NOF(XX,UROOT,WROOT)
         IF(NROOTS.GE.6) CALL ROOT6NOF(XX,NROOTS,UROOT,WROOT)
         DO NR=1,NROOTS
           RWV(1,NG,NR)=UROOT(NR)
           RWV(2,NG,NR)=WROOT(NR)
         ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
! JKBCDFNOF
      SUBROUTINE JKBCDFNOF(B00,B01,B10,C00,D00,F00,DIJ,DKL,            &
                        ABV,CV,RWV,NUMG,NROOTS,NKL,NIJ,   &
                        NMAXS,MMAXS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!      
      LOGICAL,INTENT(IN)::NMAXS,MMAXS
      INTEGER,INTENT(IN)::NUMG,NROOTS,NKL,NIJ
      DOUBLE PRECISION,DIMENSION(5,NKL*NIJ),INTENT(IN)::ABV
      DOUBLE PRECISION,DIMENSION(18,NKL*NIJ),INTENT(IN)::CV
      DOUBLE PRECISION,DIMENSION(2,NUMG,NROOTS),INTENT(IN)::RWV
      DOUBLE PRECISION,DIMENSION(NUMG,NROOTS,3),INTENT(OUT)::B00,B01,B10,C00,D00
      DOUBLE PRECISION,DIMENSION(NUMG,NROOTS,3),INTENT(OUT)::F00,DIJ,DKL
      DOUBLE PRECISION,PARAMETER::PT5=0.5D+00, ONE=1.0D+00
!
      DO NR=1,NROOTS
        DO NG=1,NUMG
          AA =ABV(1,NG)
          BB =ABV(2,NG)
          RHO=ABV(3,NG)
          QAB=ABV(4,NG)
          UU =RHO*RWV(1,NG,NR)
          WW =    RWV(2,NG,NR)
          AAUU=AA+UU
          BBUU=BB+UU
          F00(NG,NR,1)=WW*QAB
          F00(NG,NR,2)=ONE
          F00(NG,NR,3)=ONE
          RHO0  = AA*BB/(AA+BB)
          DUM2=PT5/(AA*BB+UU*(AA+BB))
          AUDUM=AAUU*DUM2
          BUDUM=BBUU*DUM2
           UDUM=  UU*DUM2
           B00(NG,NR,1)= UDUM
           B00(NG,NR,2)= UDUM
           B00(NG,NR,3)= UDUM
           B01(NG,NR,1)=AUDUM
           B01(NG,NR,2)=AUDUM
           B01(NG,NR,3)=AUDUM
           B10(NG,NR,1)=BUDUM
           B10(NG,NR,2)=BUDUM
           B10(NG,NR,3)=BUDUM
           UDUM= UDUM+ UDUM
           IF(.NOT.MMAXS) THEN
            AUDUM=AUDUM+AUDUM
            D00(NG,NR,1)= UDUM*CV( 1,NG) + AUDUM*CV( 2,NG)
            D00(NG,NR,2)= UDUM*CV( 3,NG) + AUDUM*CV( 4,NG)
            D00(NG,NR,3)= UDUM*CV( 5,NG) + AUDUM*CV( 6,NG)
           ENDIF
           IF(.NOT.NMAXS) THEN
            BUDUM=BUDUM+BUDUM
            C00(NG,NR,1)= UDUM*CV( 8,NG) + BUDUM*CV( 7,NG)
            C00(NG,NR,2)= UDUM*CV(10,NG) + BUDUM*CV( 9,NG)
            C00(NG,NR,3)= UDUM*CV(12,NG) + BUDUM*CV(11,NG)
           ENDIF
!
        ENDDO
      ENDDO
!
      DO NR=1,NROOTS
        DO NG=1,NUMG
          DIJ(NG,NR,1)=CV(13,NG)
          DIJ(NG,NR,2)=CV(14,NG)
          DIJ(NG,NR,3)=CV(15,NG)
          DKL(NG,NR,1)=CV(16,NG)
          DKL(NG,NR,2)=CV(17,NG)
          DKL(NG,NR,3)=CV(18,NG)
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
! JKGNMSNOF
      SUBROUTINE JKGNMSNOF(GNM,NG,NMAX,MMAX,B002,B012,B102,C00,D00,    &
                           F00,NMAXS,NMAXP,MMAXS,MMAXP)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      LOGICAL,INTENT(IN)::NMAXS,NMAXP,MMAXS,MMAXP
      INTEGER,INTENT(IN)::NG,NMAX,MMAX
      DOUBLE PRECISION,DIMENSION(NG),INTENT(IN)::C00,D00,F00
      DOUBLE PRECISION,DIMENSION(NG),INTENT(INOUT)::B002,B012,B102
      DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::B00
      DOUBLE PRECISION,DIMENSION(NG,NMAX-1)::B10
      DOUBLE PRECISION,DIMENSION(NG,MMAX-1)::B01
      DOUBLE PRECISION,DIMENSION(NG,NMAX,MMAX),INTENT(OUT)::GNM
!
!     ----- G(0,0) -----
!
      DO IG=1,NG
        GNM(IG,1,1)=F00(IG)
      ENDDO
      IF(NMAXS.AND.MMAXS) RETURN
      IF(NMAXS) GO TO 30
!
!     ----- G(1,0) = C00 * G(0,0) -----
!
      DO IG=1,NG
        GNM(IG,2,1)=C00(IG)*GNM(IG,1,1)
      ENDDO
!
   30 IF(MMAXS) GO TO 60
!
!     ----- G(0,1) = D00 * G(0,0) -----
!
      DO IG=1,NG
        GNM(IG,1,2)=D00(IG)*GNM(IG,1,1)
      ENDDO
      IF(NMAXS) GO TO 60
!
!     ----- G(1,1) = B00 * G(0,0) + D00 * G(1,0) -----
!
      DO IG=1,NG
        GNM(IG,2,2)=B002(IG)*GNM(IG,1,1)+D00(IG)*GNM(IG,2,1)
      ENDDO
!
   60 MAX60=MAX0(NMAX-1,MMAX-1)
      ALLOCATE(B00(NG,MAX60))
      B00(:,1)=B002(:)
      DO M=2,MAX60
        DO IG=1,NG
          B00(IG,M)=B00(IG,M-1)+B00(IG,1)
        ENDDO
      ENDDO
!
      IF(NMAXP) GO TO 120
!
!     ----- G(N+1,0) = N * B10 * G(N-1,0) + C00 * G(N,0) -----
!
      B10(:,1)=B102(:)
      DO N=2,NMAX-1
        DO IG=1,NG
          B10(IG,N)=B10(IG,N-1)+B10(IG,1)
        ENDDO
      ENDDO
      DO N=2,NMAX-1
        DO IG=1,NG
          GNM(IG,N+1,1)=B10(IG,N-1)*GNM(IG,N-1,1)+C00(IG)*GNM(IG,N,1)
        ENDDO
      ENDDO
      IF(MMAXS) GO TO 120
!
!     ----- G(N,1) = N * B00 * G(N-1,0) + D00 * G(N,0) -----
!
      DO N=2,NMAX-1
        DO IG=1,NG
          GNM(IG,N+1,2)=B00(IG,N)*GNM(IG,N,1)+D00(IG)*GNM(IG,N+1,1)
        ENDDO
      ENDDO
!
  120 IF(MMAXP) GO TO 170
!
!     ----- G(0,M+1) = M * B01 * G(0,M-1) + D00 * G(O,M) -----
!
      B01(:,1)=B012(:)
      DO M=2,MMAX-1
        DO IG=1,NG
          B01(IG,M)=B01(IG,M-1)+B01(IG,1)
        ENDDO
      ENDDO
      DO M=2,MMAX-1
        DO IG=1,NG
          GNM(IG,1,M+1)=B01(IG,M-1)*GNM(IG,1,M-1)+D00(IG)*GNM(IG,1,M)
        ENDDO
      ENDDO
      IF(NMAXS) GO TO 170
!
!     ----- G(1,M) = M * B00 * G(0,M-1) + C00 * G(0,M) -----
!
      DO M=2,MMAX-1
        DO IG=1,NG
          GNM(IG,2,M+1)=B00(IG,M)*GNM(IG,1,M)+C00(IG)*GNM(IG,1,M+1)
        ENDDO
      ENDDO
!
  170 IF(NMAXP.OR.MMAXP) THEN
        DEALLOCATE(B00)
!        DEALLOCATE(B10)
        RETURN
      ENDIF
!
!     ----- G(N+1,M) = N * B10 * G(N-1,M  ) -----
!                    +     C00 * G(N  ,M  )
!                    + M * B00 * G(N  ,M-1)
!
      DO M=2,MMAX-1
        DO N=2,NMAX-1
          DO IG=1,NG
            GNM(IG,N+1,M+1)=B10(IG,N-1)*GNM(IG,N-1,M+1)+  &
                      C00(IG    )*GNM(IG,N  ,M+1)+        &
                      B00(IG,M  )*GNM(IG,N  ,M  )
          ENDDO
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
      DEALLOCATE(B00)
      RETURN
      END
! JKXYZSNOF
      SUBROUTINE JKXYZSNOF(GIJKL,HIJKL,GNKL,HNKL,FNKL,GNM,             &
       NG,NMAX,MMAX,NIMAX,NJMAX,NKMAX,NLMAX,DIJ,DKL,EXPNDI,EXPNDK)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL,INTENT(IN)::EXPNDI,EXPNDK
      INTEGER,INTENT(IN)::NG,NMAX,MMAX,NIMAX,NJMAX,NKMAX,NLMAX
      DOUBLE PRECISION,DIMENSION(NG*NLMAX*NKMAX,NJMAX,NIMAX),INTENT(OUT)::GIJKL
      DOUBLE PRECISION,DIMENSION(NG*NLMAX*NKMAX*NJMAX,NIMAX),INTENT(OUT)::HIJKL
      DOUBLE PRECISION,DIMENSION(NG,NLMAX,NKMAX,NMAX),INTENT(OUT)::GNKL
      DOUBLE PRECISION,DIMENSION(NG*NLMAX*NKMAX,NMAX),INTENT(IN)::HNKL
      DOUBLE PRECISION,DIMENSION(NG*NLMAX*NKMAX*NMAX),INTENT(IN)::FNKL
      DOUBLE PRECISION,DIMENSION(NG,NMAX,MMAX),INTENT(INOUT)::GNM
      DOUBLE PRECISION,DIMENSION(NG),INTENT(IN)::DIJ,DKL
!
!     ----- G(N,K,L) -----
!
      IF (.NOT.EXPNDK) THEN
!
        DO NK=1,NKMAX
          DO NL=1,NLMAX
            DO  N=1,NMAX
              DO IG=1,NG
                GNKL(IG,NL,NK,N)=GNM(IG,N,NL)
              ENDDO
            ENDDO
          ENDDO
          IF(NK.EQ.NKMAX) EXIT
          MAXXXX=MMAX-NK
          DO M=1,MAXXXX
            DO N=1,NMAX
             DO IG=1,NG
              GNM(IG,N,M)=DKL(IG)*GNM(IG,N,M)+GNM(IG,N,M+1)
             ENDDO
            ENDDO
          ENDDO
        ENDDO
!
      ELSE
!
        DO NL=1,NLMAX
          DO  NK=1,NKMAX
            DO   N=1,NMAX
              DO  IG=1,NG
                GNKL(IG,NL,NK,N)=GNM(IG,N,NK)
              ENDDO
            ENDDO
          ENDDO
          IF(NL.EQ.NLMAX) EXIT
          MAXXXX=MMAX-NL
          DO  M=1,MAXXXX
            DO  N=1,NMAX
              DO IG=1,NG
                GNM(IG,N,M)=DKL(IG)*GNM(IG,N,M)+GNM(IG,N,M+1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!
      ENDIF
!
!     ----- G(I,J,K,L) -----
!
      IF (.NOT.EXPNDI) THEN
!
        DO NI=1,NIMAX
          DO IGLKJ=1,NG*NLMAX*NKMAX*NJMAX
            HIJKL(IGLKJ,NI)=FNKL(IGLKJ)
          ENDDO
          IF(NI.EQ.NIMAX) EXIT
          MAXXXX=NMAX-NI
          DO N=1,MAXXXX
            DO NK=1,NKMAX
              DO NL=1,NLMAX
                DO IG=1,NG
                  GNKL(IG,NL,NK,N)=DIJ(IG)*GNKL(IG,NL,NK,N)+  &
                                   GNKL(IG,NL,NK,N+1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!
      ELSE
!
        DO NJ=1,NJMAX
          DO NI=1,NIMAX
            DO IGLK=1,NG*NLMAX*NKMAX
              GIJKL(IGLK,NJ,NI)=HNKL(IGLK,NI)
            ENDDO
          ENDDO
          IF(NJ.EQ.NJMAX) EXIT
          MAXXXX=NMAX-NJ
          DO N=1,MAXXXX
            DO NK=1,NKMAX
              DO NL=1,NLMAX
                DO IG=1,NG
                  GNKL(IG,NL,NK,N)=DIJ(IG)*GNKL(IG,NL,NK,N)+  &
                                   GNKL(IG,NL,NK,N+1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END
! JDXYZSNOF
      SUBROUTINE JDXYZSNOF(GI,GIJ,GIJK,GIJKL,                          &
       NG,NIMAX,NJMAX,NKMAX,NLMAX,NI,NJ,NK,NL,AAI,AAJ,AAK,AAL,         &
       FI,FJ,FK,FL,SKIPI,SKIPJ,SKIPK,SKIPL)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER,INTENT(IN)::NG,NIMAX,NJMAX,NKMAX,NLMAX,NI,NJ,NK,NL
      LOGICAL,INTENT(IN)::SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL::IS,JS,KS,LS
      DOUBLE PRECISION,DIMENSION(NG,NLMAX,NKMAX*NJMAX*NIMAX),INTENT(IN)::GIJKL
      DOUBLE PRECISION,DIMENSION(NG,NLMAX,NKMAX,NJMAX*NIMAX),INTENT(IN)::GIJK
      DOUBLE PRECISION,DIMENSION(NG,NLMAX*NKMAX,NJMAX,NIMAX),INTENT(IN)::GIJ
      DOUBLE PRECISION,DIMENSION(NG,NLMAX*NKMAX*NJMAX,NIMAX),INTENT(IN)::GI
      DOUBLE PRECISION,DIMENSION(NG),INTENT(IN)::AAI,AAJ,AAK,AAL
      DOUBLE PRECISION,DIMENSION(NG,NLMAX*NKMAX*NJMAX,NIMAX),INTENT(OUT)::FI
      DOUBLE PRECISION,DIMENSION(NG,NLMAX*NKMAX,NJMAX,NIMAX),INTENT(OUT)::FJ
      DOUBLE PRECISION,DIMENSION(NG,NLMAX,NKMAX,NJMAX*NIMAX),INTENT(OUT)::FK
      DOUBLE PRECISION,DIMENSION(NG,NLMAX,NKMAX*NJMAX*NIMAX),INTENT(OUT)::FL
!
      IS=NI.EQ.1
      JS=NJ.EQ.1
      KS=NK.EQ.1
      LS=NL.EQ.1
!
!     ----- FIRST DERIVATIVES ONLY -----
!
      IF(SKIPI) GO TO 1030
!
!     ----- -FI- ONLY -----
!
      DO LKJ=1,NLMAX*NKMAX*NJMAX
        DO IG =1,NG
           FI(IG,LKJ,1)=  GI(IG,LKJ,2)*AAI(IG)
        ENDDO
      ENDDO
      IF(IS) GO TO 1030
      DO I  =2,NI
        DO LKJ=1,NLMAX*NKMAX*NJMAX
          DO IG =1,NG
             FI(IG,LKJ,I)= GI(IG,LKJ,I+1)*AAI(IG)  &
                          -GI(IG,LKJ,I-1)*(I-1)
          ENDDO
        ENDDO
      ENDDO
!
 1030 IF(SKIPJ) GO TO 1130
!
!     ----- -FJ- ONLY -----
!
      DO I =1,NIMAX
        DO LK=1,NLMAX*NKMAX
          DO IG=1,NG
             FJ(IG,LK,1,I)=  GIJ(IG,LK,2,I)*AAJ(IG)
          ENDDO
        ENDDO
      ENDDO
      IF(JS) GO TO 1130
      DO I =1,NIMAX
        DO J =2,NJ
          DO LK=1,NLMAX*NKMAX
            DO IG=1,NG
               FJ(IG,LK,J,I)= GIJ(IG,LK,J+1,I)*AAJ(IG) &
                             -GIJ(IG,LK,J-1,I)*(J-1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
 1130 IF(SKIPK) GO TO 1230
!
!     ----- -FK- ONLY -----
!
      DO JI=1,NJMAX*NIMAX
        DO L =1,NLMAX
          DO IG=1,NG
             FK(IG,L,1,JI)=  GIJK(IG,L,2,JI)*AAK(IG)
          ENDDO
        ENDDO
      ENDDO
      IF(KS) GO TO 1230
      DO JI=1,NJMAX*NIMAX
        DO K =2,NK
          DO L =1,NLMAX
            DO IG=1,NG
               FK(IG,L,K,JI)= GIJK(IG,L,K+1,JI)*AAK(IG)  &
                             -GIJK(IG,L,K-1,JI)*(K-1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
 1230 IF(SKIPL) GO TO 1330
!
!     ----- -FL- AND -SLL- -----
!
      DO KJI=1,NKMAX*NJMAX*NIMAX
        DO IG =1,NG
           FL(IG,1,KJI)=  GIJKL(IG,2,KJI)*AAL(IG)
        ENDDO
      ENDDO
      IF(LS) GO TO 1330
      DO KJI=1,NKMAX*NJMAX*NIMAX
        DO L  =2,NL
          DO IG =1,NG
             FL(IG,L,KJI)= GIJKL(IG,L+1,KJI)*AAL(IG)   &
                          -GIJKL(IG,L-1,KJI)*(L-1)
          ENDDO
        ENDDO
      ENDDO
!
 1330 CONTINUE
!----------------------------------------------------------------------- 
      RETURN
      END
! DSPDFSNOF
      SUBROUTINE DSPDFSNOF(IIFINT,NG,NR,IJKLG,GIJKL,FIJKL,XYZ,         &
       FIXYZ,FJXYZ,FKXYZ,FLXYZ,DIJSI,DIJSJ,DKLSK,DKLSL,                &
       DAB,SKIPI,SKIPJ,SKIPK,SKIPL,MINI,MINJ,MINK,MINL,                &
       MAXI,MAXJ,MAXK,MAXL,MAXNUM,MAXXYZ,IIEQJJ,KKEQLL,                &
       IJEQKL,SPI,SPJ,SPK,SPL,SPIJKL,FD,IFI)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL::IS,JS,KS,LS
      LOGICAL::IJS,IJKS,IJKLS
      LOGICAL,INTENT(IN)::IIEQJJ,KKEQLL,IJEQKL
      LOGICAL,INTENT(IN)::SPI,SPJ,SPK,SPL,SPIJKL
      LOGICAL,INTENT(IN)::SKIPI,SKIPJ,SKIPK,SKIPL
      INTEGER,INTENT(IN)::MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      INTEGER,INTENT(IN)::NG,NR,MAXNUM,MAXXYZ,IFI
      INTEGER,INTENT(INOUT)::IIFINT
      INTEGER,DIMENSION(4,MAXNUM),INTENT(IN)::IJKLG
      DOUBLE PRECISION,DIMENSION(12),INTENT(INOUT)::FD
      DOUBLE PRECISION,DIMENSION(MAXNUM),INTENT(IN)::DAB
      DOUBLE PRECISION,DIMENSION(MAXNUM),INTENT(OUT)::GIJKL
      DOUBLE PRECISION,DIMENSION(12,MAXNUM/4),INTENT(OUT)::FIJKL
!     CAUTION, XYZ WHEN CALLING IS GIJKL, BUT HERE GIJKL STORES INTEGRALS
      DOUBLE PRECISION,DIMENSION(NG*NR,IFI-NG*NR),INTENT(IN)::XYZ
      DOUBLE PRECISION,DIMENSION(NG*NR,IFI-NG*NR),INTENT(IN)::FIXYZ,FJXYZ,FKXYZ,FLXYZ
      DOUBLE PRECISION,DIMENSION(MAXXYZ),INTENT(IN)::DIJSI,DIJSJ,DKLSK,DKLSL
      DOUBLE PRECISION,DIMENSION(NG*NR)::XY,XZ,YZ,SJ,SK,SL
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00
!
      IF(SPIJKL) GO TO 1000
!
!     ----- NO SHARED EXPONENTS ; SUM UP ( IX * IY * IZ ) -----
!
!
!     ----- GRADIENT -----
!
      IJKLN=0
      DO I=MINI,MAXI
        JMAX=MAXJ
        IF(IIEQJJ) JMAX=I
        DO J=MINJ,JMAX
          KMAX=MAXK
          IF(IJEQKL) KMAX=I
          DO K=MINK,KMAX
            LMAX=MAXL
            IF(KKEQLL           ) LMAX=K
            IF(IJEQKL.AND.K.EQ.I) LMAX=J
            DO L=MINL,LMAX
              IJKLN=IJKLN+1
              NN=IJKLG(1,IJKLN)
              NX=IJKLG(2,IJKLN)
              NY=IJKLG(3,IJKLN)
              NZ=IJKLG(4,IJKLN)
!
              DO IGR=1,NG*NR
                XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
                XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
                YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
              ENDDO
!
              IF(SKIPI) GO TO 530
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO IGR=1,NG*NR
                DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
                DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
              ENDDO
!
              IIFINT=IIFINT+3
              FD( 1)=FD( 1)+DAB(NN)*DUMFX
              FD( 2)=FD( 2)+DAB(NN)*DUMFY
              FD( 3)=FD( 3)+DAB(NN)*DUMFZ
  530         IF(SKIPJ) GO TO 550
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO IGR=1,NG*NR
                DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
                DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
              ENDDO
!
              IIFINT=IIFINT+3 
              FD( 4)=FD( 4)+DAB(NN)*DUMFX
              FD( 5)=FD( 5)+DAB(NN)*DUMFY
              FD( 6)=FD( 6)+DAB(NN)*DUMFZ
  550         IF(SKIPK) GO TO 570
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO IGR=1,NG*NR
                DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
                DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
              ENDDO
!
              IIFINT=IIFINT+3
              FD( 7)=FD( 7)+DAB(NN)*DUMFX
              FD( 8)=FD( 8)+DAB(NN)*DUMFY
              FD( 9)=FD( 9)+DAB(NN)*DUMFZ
  570         IF(SKIPL) GO TO 600
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO IGR=1,NG*NR
                DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
                DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
              ENDDO
!
              IIFINT=IIFINT+3
              FD(10)=FD(10)+DAB(NN)*DUMFX
              FD(11)=FD(11)+DAB(NN)*DUMFY
              FD(12)=FD(12)+DAB(NN)*DUMFZ
  600         CONTINUE
!
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ----- ZEROTH AND FIRST DERIVATIVE INTEGRALS -----
!
!     REMOVE NEXT RETURN FOR STORING EACH INTEGRAL CONTRIBUTION
      RETURN
!
      IJKLN=0
      DO 940 I=MINI,MAXI
        JMAX=MAXJ
        IF(IIEQJJ) JMAX=I
        DO 930 J=MINJ,JMAX
          KMAX=MAXK
          IF(IJEQKL) KMAX=I
          DO 920 K=MINK,KMAX
            LMAX=MAXL
            IF(KKEQLL           ) LMAX=K
            IF(IJEQKL.AND.K.EQ.I) LMAX=J
            DO 910 L=MINL,LMAX
            IJKLN=IJKLN+1
            NN=IJKLG(1,IJKLN)
            NX=IJKLG(2,IJKLN)
            NY=IJKLG(3,IJKLN)
            NZ=IJKLG(4,IJKLN)
!
            DO IGR=1,NG*NR
              XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
              XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
              YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
            ENDDO
            DUM=ZERO
            DO IGR=1,NG*NR
              DUM=DUM+XYZ(IGR,NX)*YZ(IGR)
            ENDDO
            GIJKL(NN)=GIJKL(NN)+DUM
!
            IF(SKIPI) GO TO 830
            DUMFX =ZERO
            DUMFY =ZERO
            DUMFZ =ZERO
            DO IGR=1,NG*NR
              DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
              DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
              DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
            ENDDO
            FIJKL( 1,NN)=FIJKL( 1,NN)+DUMFX
            FIJKL( 2,NN)=FIJKL( 2,NN)+DUMFY
            FIJKL( 3,NN)=FIJKL( 3,NN)+DUMFZ
  830       IF(SKIPJ) GO TO 850
            DUMFX =ZERO
            DUMFY =ZERO
            DUMFZ =ZERO
            DO IGR=1,NG*NR
              DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
              DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
              DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
            ENDDO
            FIJKL( 4,NN)=FIJKL( 4,NN)+DUMFX
            FIJKL( 5,NN)=FIJKL( 5,NN)+DUMFY
            FIJKL( 6,NN)=FIJKL( 6,NN)+DUMFZ
  850       IF(SKIPK) GO TO 870
            DUMFX =ZERO
            DUMFY =ZERO
            DUMFZ =ZERO
            DO IGR=1,NG*NR
              DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
              DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
              DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
            ENDDO
            FIJKL( 7,NN)=FIJKL( 7,NN)+DUMFX
            FIJKL( 8,NN)=FIJKL( 8,NN)+DUMFY
            FIJKL( 9,NN)=FIJKL( 9,NN)+DUMFZ
  870       IF(SKIPL) GO TO 900
            DUMFX =ZERO
            DUMFY =ZERO
            DUMFZ =ZERO
            DO IGR=1,NG*NR
              DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
              DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
              DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
            ENDDO
            FIJKL(10,NN)=FIJKL(10,NN)+DUMFX
            FIJKL(11,NN)=FIJKL(11,NN)+DUMFY
            FIJKL(12,NN)=FIJKL(12,NN)+DUMFZ
  900       CONTINUE
!
  910       CONTINUE
  920     CONTINUE
  930   CONTINUE
  940 CONTINUE
!
      RETURN
!
!     ----- SHARED EXPONENTS ; FORM ( IX * IY * IZ ) -----
!
 1000 CONTINUE
!
!     ----- GRADIENT -----
!
      IJKLN=0
      DO 2640 I=MINI,MAXI
        IS=SPI.AND.I.EQ.1
!
        JMAX=MAXJ
        IF(IIEQJJ) JMAX=I
        DO 2630 J=MINJ,JMAX
          JS=SPJ.AND.J.EQ.1
!
          IF(JS) THEN
             IF(IS) THEN
                DO IGR=1,NG*NR
                  SJ(IGR)=DIJSJ(IGR)*DIJSI(IGR)
                ENDDO
             ELSE
                DO IGR=1,NG*NR
                  SJ(IGR)=DIJSJ(IGR)
                ENDDO
             ENDIF
          ELSE
             IF(IS) THEN
                DO IGR=1,NG*NR
                  SJ(IGR)=DIJSI(IGR)
                ENDDO
             ENDIF
          ENDIF
          IJS=IS.OR.JS
!
          KMAX=MAXK
          IF(IJEQKL) KMAX=I
          DO 2620 K=MINK,KMAX
            KS=SPK.AND.K.EQ.1
!
            IF(KS) THEN
               IF(IJS) THEN
                  DO IGR=1,NG*NR
                    SK(IGR)=DKLSK(IGR)*SJ(IGR)
                  ENDDO
               ELSE
                  DO IGR=1,NG*NR
                    SK(IGR)=DKLSK(IGR)
                  ENDDO
               ENDIF
            ELSE
               IF(IJS) THEN
                  DO IGR=1,NG*NR
                    SK(IGR)=SJ(IGR)
                  ENDDO
               ENDIF
            ENDIF
            IJKS=IJS.OR.KS
!
            LMAX=MAXL
            IF(KKEQLL           ) LMAX=K
            IF(IJEQKL.AND.K.EQ.I) LMAX=J
            DO 2610 L=MINL,LMAX
              LS=SPL.AND.L.EQ.1
!
              IF(LS) THEN
                 IF(IJKS) THEN
                    DO IGR=1,NG*NR
                      SL(IGR)=DKLSL(IGR)*SK(IGR)
                    ENDDO
                 ELSE
                    DO IGR=1,NG*NR
                      SL(IGR)=DKLSL(IGR)
                    ENDDO
                 ENDIF
              ELSE
                 IF(IJKS) THEN
                    DO IGR=1,NG*NR
                      SL(IGR)=SK(IGR)
                    ENDDO
                 ENDIF
              ENDIF
              IJKLS=IJKS.OR.LS
!
              IJKLN=IJKLN+1
              NN=IJKLG(1,IJKLN)
              NX=IJKLG(2,IJKLN)
              NY=IJKLG(3,IJKLN)
              NZ=IJKLG(4,IJKLN)
!
              IF(IJKLS) THEN
                 DO IGR=1,NG*NR
                   XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)*SL(IGR)
                   XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)*SL(IGR)
                   YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)*SL(IGR)
                 ENDDO
              ELSE
                 DO IGR=1,NG*NR
                   XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
                   XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
                   YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
                 ENDDO
              ENDIF
!
              IF(SKIPI) GO TO 2530
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO IGR=1,NG*NR
                DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
                DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
              ENDDO
!
              IIFINT=IIFINT+3
              FD( 1)=FD( 1)+DAB(NN)*DUMFX
              FD( 2)=FD( 2)+DAB(NN)*DUMFY
              FD( 3)=FD( 3)+DAB(NN)*DUMFZ
 2530         IF(SKIPJ) GO TO 2550
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO IGR=1,NG*NR
                DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
                DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
              ENDDO
!
              IIFINT=IIFINT+3 
              FD( 4)=FD( 4)+DAB(NN)*DUMFX
              FD( 5)=FD( 5)+DAB(NN)*DUMFY
              FD( 6)=FD( 6)+DAB(NN)*DUMFZ
 2550         IF(SKIPK) GO TO 2570
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO IGR=1,NG*NR
                DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
                DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
              ENDDO
!
              IIFINT=IIFINT+3
              FD( 7)=FD( 7)+DAB(NN)*DUMFX
              FD( 8)=FD( 8)+DAB(NN)*DUMFY
              FD( 9)=FD( 9)+DAB(NN)*DUMFZ
 2570         IF(SKIPL) GO TO 2600
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO IGR=1,NG*NR
                DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
                DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
              ENDDO
!  
              IIFINT=IIFINT+3
              FD(10)=FD(10)+DAB(NN)*DUMFX
              FD(11)=FD(11)+DAB(NN)*DUMFY
              FD(12)=FD(12)+DAB(NN)*DUMFZ
 2600         CONTINUE
!
 2610       CONTINUE
 2620     CONTINUE
 2630   CONTINUE
 2640 CONTINUE
!
!     ----- ZEROTH AND FIRST DERIVATIVE INTEGRALS -----
!
!     REMOVE NEXT RETURN FOR STORING EACH INTEGRAL CONTRIBUTION
      RETURN
!
      IJKLN=0
      DO 3640 I=MINI,MAXI
        IS=SPI.AND.I.EQ.1
!
        JMAX=MAXJ
        IF(IIEQJJ) JMAX=I
        DO 3630 J=MINJ,JMAX
          JS=SPJ.AND.J.EQ.1
!
          IF(JS) THEN
             IF(IS) THEN
                DO IGR=1,NG*NR
                  SJ(IGR)=DIJSJ(IGR)*DIJSI(IGR)
                ENDDO
             ELSE
                DO IGR=1,NG*NR
                  SJ(IGR)=DIJSJ(IGR)
                ENDDO
             ENDIF
          ELSE
             IF(IS) THEN
                DO IGR=1,NG*NR
                  SJ(IGR)=DIJSI(IGR)
                ENDDO
             ENDIF
          ENDIF
          IJS=IS.OR.JS
!
          KMAX=MAXK
          IF(IJEQKL) KMAX=I
          DO 3620 K=MINK,KMAX
            KS=SPK.AND.K.EQ.1
!
            IF(KS) THEN
               IF(IJS) THEN
                  DO IGR=1,NG*NR
                    SK(IGR)=DKLSK(IGR)*SJ(IGR)
                  ENDDO
               ELSE
                  DO IGR=1,NG*NR
                    SK(IGR)=DKLSK(IGR)
                  ENDDO
               ENDIF
            ELSE
               IF(IJS) THEN
                  DO IGR=1,NG*NR
                    SK(IGR)=SJ(IGR)
                  ENDDO
               ENDIF
            ENDIF
            IJKS=IJS.OR.KS
!
            LMAX=MAXL
            IF(KKEQLL           ) LMAX=K
            IF(IJEQKL.AND.K.EQ.I) LMAX=J
            DO 3610 L=MINL,LMAX
              LS=SPL.AND.L.EQ.1
!
              IF(LS) THEN
                 IF(IJKS) THEN
                    DO IGR=1,NG*NR
                      SL(IGR)=DKLSL(IGR)*SK(IGR)
                    ENDDO
                 ELSE
                    DO IGR=1,NG*NR
                      SL(IGR)=DKLSL(IGR)
                    ENDDO
                 ENDIF
              ELSE
                 IF(IJKS) THEN
                    DO IGR=1,NG*NR
                      SL(IGR)=SK(IGR)
                    ENDDO
                 ENDIF
              ENDIF
              IJKLS=IJKS.OR.LS
!
              IJKLN=IJKLN+1
              NN=IJKLG(1,IJKLN)
              NX=IJKLG(2,IJKLN)
              NY=IJKLG(3,IJKLN)
              NZ=IJKLG(4,IJKLN)
!
              IF(IJKLS) THEN
                 DO IGR=1,NG*NR
                   XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)*SL(IGR)
                   XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)*SL(IGR)
                   YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)*SL(IGR)
                 ENDDO
              ELSE
                 DO IGR=1,NG*NR
                   XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
                   XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
                   YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
                 ENDDO
              ENDIF
!
              DUM=ZERO
              DO IGR=1,NG*NR
                DUM=DUM+XYZ(IGR,NX)*YZ(IGR)
              ENDDO
              GIJKL(NN)=GIJKL(NN)+DUM
!
              IF(SKIPI) GO TO 3530
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO IGR=1,NG*NR
                DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
                DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
              ENDDO
              FIJKL( 1,NN)=FIJKL( 1,NN)+DUMFX
              FIJKL( 2,NN)=FIJKL( 2,NN)+DUMFY
              FIJKL( 3,NN)=FIJKL( 3,NN)+DUMFZ
 3530         IF(SKIPJ) GO TO 3550
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO IGR=1,NG*NR
                DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
                DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
              ENDDO
              FIJKL( 4,NN)=FIJKL( 4,NN)+DUMFX
              FIJKL( 5,NN)=FIJKL( 5,NN)+DUMFY
              FIJKL( 6,NN)=FIJKL( 6,NN)+DUMFZ
 3550         IF(SKIPK) GO TO 3570
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO IGR=1,NG*NR
                DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
                DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
              ENDDO
              FIJKL( 7,NN)=FIJKL( 7,NN)+DUMFX
              FIJKL( 8,NN)=FIJKL( 8,NN)+DUMFY
              FIJKL( 9,NN)=FIJKL( 9,NN)+DUMFZ
 3570         IF(SKIPL) GO TO 3600
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO IGR=1,NG*NR
                DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
                DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
              ENDDO
              FIJKL(10,NN)=FIJKL(10,NN)+DUMFX
              FIJKL(11,NN)=FIJKL(11,NN)+DUMFY
              FIJKL(12,NN)=FIJKL(12,NN)+DUMFZ
 3600         CONTINUE
!
 3610       CONTINUE
 3620     CONTINUE
 3630   CONTINUE
 3640 CONTINUE
!-----------------------------------------------------------------------
      RETURN
      END
! DABNOF2
      SUBROUTINE DABNOF2(II,JJ,KK,LL,PNRM,KKMIN,KKMAX,KLOC,            &
                        IGXYZ,JGXYZ,KGXYZ,LGXYZ,IIEQJJ,KKEQLL,IJEQKL,  &
                        IA,DA,DAB,MAXNUM,DABMAX,DAAUX,DAAUX2)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,INTENT(IN)::II,JJ,KK,LL
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KLOC,KKMIN,KKMAX
      INTEGER,DIMENSION(4,35),INTENT(IN)::IGXYZ,JGXYZ,KGXYZ,LGXYZ
      INTEGER,DIMENSION(NBF),INTENT(IN)::IA
      INTEGER,INTENT(IN)::MAXNUM
      LOGICAL,INTENT(IN)::IIEQJJ,KKEQLL,IJEQKL
      DOUBLE PRECISION,DIMENSION(NBF,NBFT),INTENT(IN)::DA
      DOUBLE PRECISION,DIMENSION(NBF,NBFT),INTENT(IN)::DAAUX,DAAUX2
      DOUBLE PRECISION,INTENT(OUT)::DABMAX
      DOUBLE PRECISION,DIMENSION(84),INTENT(IN)::PNRM
      DOUBLE PRECISION,DIMENSION(MAXNUM),INTENT(OUT)::DAB
      DOUBLE PRECISION,PARAMETER::ZER=0.0D+00,PT5=0.5D+00
!-----------------------------------------------------------------------
      DABMAX=ZER         
!
!     GET 2e DENSITY FOR THIS SHELL BLOCK
!
      MINI= KKMIN(II)
      MINJ= KKMIN(JJ)
      MINK= KKMIN(KK)
      MINL= KKMIN(LL)
      LOCI= KLOC(II)-MINI
      LOCJ= KLOC(JJ)-MINJ
      LOCK= KLOC(KK)-MINK
      LOCL= KLOC(LL)-MINL
      MAXI= KKMAX(II)
      MAXJ= KKMAX(JJ)
      MAXK= KKMAX(KK)
      MAXL= KKMAX(LL)
         DO I=MINI,MAXI
            P1I = PNRM(I)
            JMAX= MAXJ
            IF(IIEQJJ) JMAX= I
            DO J=MINJ,JMAX
               P2J = P1I*PNRM(J)
               IAJ= MAX0(LOCI+I,LOCJ+J)
               IIJ= MIN0(LOCI+I,LOCJ+J)
               KMMAX=MAXK
               IF(IJEQKL) KMMAX= I
               DO K=MINK,KMMAX
                  P3K = P2J*PNRM(K)
                  LMAX= MAXL
                  IF(KKEQLL) LMAX= K
                  IF(IJEQKL .AND. K.EQ.I) LMAX= J
                  DO L=MINL,LMAX
                     P4L= P3K*PNRM(L)
                     KAL= MAX0(LOCK+K,LOCL+L)
                     KIL= MIN0(LOCK+K,LOCL+L)
                     IN = IAJ
                     JN = IIJ
                     KN = KAL
                     LN = KIL
                     IF(IN.LT.KN .OR.(IN.EQ.KN .AND. JN.LT.LN)) THEN
                        IN = KAL
                        JN = KIL
                        KN = IAJ
                        LN = IIJ
                     ENDIF
                     IJ = IA(IN)+JN
                     IK = IA(IN)+KN
                     IL = IA(IN)+LN
                     JK = IA(MAX0(JN,KN))+MIN0(JN,KN)
                     JL = IA(JN)+LN
                     IF(JN.LT.KN) JL = IA(MAX0(JN,LN))+MIN0(JN,LN)
                     KL = IA(KN)+LN
!                    CONTRACT OVER THE 2nd NO COEFFICIENT
                     DF1=ZER
                     DQ1=ZER
                     DO LQ=1,NBF5
                      DF1 = DF1 + DAAUX(LQ,KL)*DA(LQ,IJ)
                      DQ1 = DQ1 + DAAUX2(LQ,JK)*DA(LQ,IL)
                      DQ1 = DQ1 + DAAUX2(LQ,JL)*DA(LQ,IK)
                     ENDDO
!                    BUILD THE DENSITY TERM SUMMING COULOMB-LIKE
!                    AND EXCHANGE-LIKE PARTS
                     DF1 = DF1 + DF1 - DQ1
!                    AVOID DOUBLE COUNTING OF DIAGONAL TERMS                     
                     IF(JN.EQ.IN               ) DF1= DF1*PT5
                     IF(LN.EQ.KN               ) DF1= DF1*PT5
                     IF(KN.EQ.IN .AND. LN.EQ.JN) DF1= DF1*PT5
                     IF(DABMAX.LT. ABS(DF1)) DABMAX= ABS(DF1)
! IGXYZ AND J, K, AND L ARE SET UP IN JKDNDX
                     IJKL=IGXYZ(1,I)+JGXYZ(1,J)+KGXYZ(1,K)+LGXYZ(1,L)
                     DAB(IJKL)= DF1*P4L
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
!-----------------------------------------------------------------------
      RETURN
      END        
! DABNOF2PRE
      SUBROUTINE DABNOF2PRE(CJ12,CK12,RO,DA,DAAUX,DAAUX2)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5),INTENT(IN)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBFT),INTENT(IN)::DA
      DOUBLE PRECISION,DIMENSION(NBF,NBFT),INTENT(OUT)::DAAUX,DAAUX2
      DOUBLE PRECISION,PARAMETER::ZER=0.0D+00
!-----------------------------------------------------------------------
!      
!     FIRST CONTRACT OVER ONE NATURAL ORBITAL COEFFICIENT
!
      DAAUX=ZER
      DAAUX2=ZER
       DO IJ=1,NBFT
        DO LP=1,NBF5
          DAAUX(LP,IJ) = DAAUX(LP,IJ) + RO(LP)*DA(LP,IJ)
          DO LQ=1,LP-1
            DAAUX(LQ,IJ) = DAAUX(LQ,IJ) + CJ12(LP,LQ)*DA(LP,IJ)
            DAAUX2(LQ,IJ) = DAAUX2(LQ,IJ) + CK12(LP,LQ)*DA(LP,IJ)
          ENDDO
          DO LQ=LP+1,NBF5
            DAAUX(LQ,IJ) = DAAUX(LQ,IJ) + CJ12(LP,LQ)*DA(LP,IJ)
            DAAUX2(LQ,IJ) = DAAUX2(LQ,IJ) + CK12(LP,LQ)*DA(LP,IJ)
          ENDDO
        ENDDO
       ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
! DABNOF5PRE
      SUBROUTINE DABNOF5PRE(RO,DA,CJAUX,CKAUX)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Extended (Nc>1): NBF5 = NBF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO
      DOUBLE PRECISION,DIMENSION(NBF5)::BETA
      DOUBLE PRECISION,DIMENSION(NBF,NBFT),INTENT(IN)::DA
      DOUBLE PRECISION,DIMENSION(NBFT)::DENS5
      DOUBLE PRECISION,DIMENSION(NDOC,NBFT)::DAAUX2
      DOUBLE PRECISION,DIMENSION(NCO,NBFT)::DAAUX
      DOUBLE PRECISION,DIMENSION(NBFT,NBFT),INTENT(OUT)::CJAUX,CKAUX
      DOUBLE PRECISION,PARAMETER::ZER=0.0D+00
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,NBF5
        BETA(i)=DSQRT(RO(i))
      ENDDO
!     DENS5 AND DAAUX ARE RELATED WITH (NN-DELTA)
!     DAAUX2 IS RELATED TO PI MATRIX
      DENS5=ZER
      DAAUX=ZER
      DAAUX2=ZER
      CJAUX = ZER
      CKAUX = ZER
      DO IJ=1,NBFT
!     ADD HF-LIKE TERMS FROM FROZEN COEFFICIENTS      
       DO LP=1,NO1
          DENS5(IJ) = DENS5(IJ) + RO(LP)*DA(LP,IJ)
          DAAUX(LP,IJ) = DAAUX(LP,IJ) + RO(LP)*DA(LP,IJ)
       ENDDO
       DO j=1,NDOC
        jn = NO1+j
        DENS5(IJ)=DENS5(IJ)+RO(jn)*DA(jn,IJ)
        DAAUX(jn,IJ)=DAAUX(jn,IJ)+RO(jn)*DA(jn,IJ)
        DAAUX2(j,IJ)=DAAUX2(j,IJ)-BETA(jn)*DA(jn,IJ)
        DO i=NDOC+NCWO*(NDOC-j)+1,NDOC+NCWO*(NDOC-j+1)
          in = NO1+i
          DENS5(IJ)=DENS5(IJ)+RO(in)*DA(in,IJ)
          DAAUX(jn,IJ)=DAAUX(jn,IJ)+RO(in)*DA(in,IJ)
          DAAUX2(j,IJ)=DAAUX2(j,IJ)+BETA(in)*DA(in,IJ)
        ENDDO
       ENDDO
!       
       DO KL=1,IJ
        CJAUX(IJ,KL) = DENS5(IJ)*DENS5(KL)
        DO IG=1,NDOC
         IGG = IG + NO1
         CJAUX(IJ,KL) = CJAUX(IJ,KL) - DAAUX(IGG,IJ)*DAAUX(IGG,KL)
         CKAUX(IJ,KL) = CKAUX(IJ,KL) + DAAUX2(IG,IJ)*DAAUX2(IG,KL)
        ENDDO
! AACOMP        
        CKAUX(IJ,KL) = CJAUX(IJ,KL) - CKAUX(IJ,KL)
        CJAUX(IJ,KL) = CJAUX(IJ,KL) + CJAUX(IJ,KL)
! SYMMETRY IN AO INDICES         
        CKAUX(KL,IJ) = CKAUX(IJ,KL)
        CJAUX(KL,IJ) = CJAUX(IJ,KL)
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
! DABNOF7PRE
      SUBROUTINE DABNOF7PRE(RO,DA,CJAUX,CKAUX)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Extended (Nc>1): NBF5 = NBF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO
      DOUBLE PRECISION,DIMENSION(NBF5)::BETA
      DOUBLE PRECISION,DIMENSION(NBF,NBFT),INTENT(IN)::DA
      DOUBLE PRECISION,DIMENSION(NBFT)::DENS5,DENS7
      DOUBLE PRECISION,DIMENSION(NDOC,NBFT)::DAAUX2,DAAUX3
      DOUBLE PRECISION,DIMENSION(NCO,NBFT)::DAAUX
      DOUBLE PRECISION,DIMENSION(NBFT,NBFT),INTENT(OUT)::CJAUX,CKAUX
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::FIs
      DOUBLE PRECISION,PARAMETER::ZER=0.0D+00
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,NBF5
        BETA(i)=DSQRT(RO(i))
      ENDDO
      ALLOCATE(FIs(NBF5))
      FIs = 0.0d0      
      if(Ista==0)then
       DO j=NO1+1,NBF5
        FIs(j) = DSQRT( RO(j)*(1.0d0-RO(j)) )
       ENDDO
      else if(Ista==1)then
       DO j=NO1+1,NBF5
        FIs(j) = 2.0d0*RO(j)*(1.0d0-RO(j))
       ENDDO
      end if
!
      DENS5=ZER
      DENS7=ZER
      DAAUX=ZER
      DAAUX2=ZER
      DAAUX3=ZER
      CJAUX = ZER
      CKAUX = ZER
      DO IJ=1,NBFT
       DO LP=1,NO1
          DENS5(IJ) = DENS5(IJ) + RO(LP)*DA(LP,IJ)
          DAAUX(LP,IJ) = DAAUX(LP,IJ) + RO(LP)*DA(LP,IJ)
       ENDDO
       DO j=1,NDOC
        jn = NO1+j
        DENS5(IJ)=DENS5(IJ)+RO(jn)*DA(jn,IJ)
        DENS7(IJ)=DENS7(IJ)+FIs(jn)*DA(jn,IJ)
        DAAUX(jn,IJ)=DAAUX(jn,IJ)+RO(jn)*DA(jn,IJ)
        DAAUX2(j,IJ)=DAAUX2(j,IJ)-BETA(jn)*DA(jn,IJ)
        DAAUX3(j,IJ)=DAAUX3(j,IJ)+FIs(jn)*DA(jn,IJ)
        DO i=NDOC+NCWO*(NDOC-j)+1,NDOC+NCWO*(NDOC-j+1)
          in = NO1+i
          DENS5(IJ)=DENS5(IJ)+RO(in)*DA(in,IJ)
          DENS7(IJ)=DENS7(IJ)+FIs(in)*DA(in,IJ)
          DAAUX(jn,IJ)=DAAUX(jn,IJ)+RO(in)*DA(in,IJ)
          DAAUX2(j,IJ)=DAAUX2(j,IJ)+BETA(in)*DA(in,IJ)
          DAAUX3(j,IJ)=DAAUX3(j,IJ)+FIs(in)*DA(in,IJ)
        ENDDO
       ENDDO
!       
       DO KL=1,IJ
        CJAUX(IJ,KL) = DENS5(IJ)*DENS5(KL)
        CKAUX(IJ,KL) = - DENS7(IJ)*DENS7(KL)
        DO IG=1,NDOC
         IGG = IG + NO1
         CJAUX(IJ,KL) = CJAUX(IJ,KL) - DAAUX(IGG,IJ)*DAAUX(IGG,KL)
         CKAUX(IJ,KL) = CKAUX(IJ,KL) + DAAUX2(IG,IJ)*DAAUX2(IG,KL)
         CKAUX(IJ,KL) = CKAUX(IJ,KL) + DAAUX3(IG,IJ)*DAAUX3(IG,KL)
        ENDDO
! AACOMP
        CKAUX(IJ,KL) = CJAUX(IJ,KL) - CKAUX(IJ,KL)
        CJAUX(IJ,KL) = CJAUX(IJ,KL) + CJAUX(IJ,KL)
! SYMMETRY IN AO INDICES
        CKAUX(KL,IJ) = CKAUX(IJ,KL)
        CJAUX(KL,IJ) = CJAUX(IJ,KL)
       ENDDO
      ENDDO
      DEALLOCATE(FIs)
!-----------------------------------------------------------------------
      RETURN
      END
! DABNOF5
      SUBROUTINE DABNOF5(II,JJ,KK,LL,PNRM,KKMIN,KKMAX,KLOC,            &
                        IGXYZ,JGXYZ,KGXYZ,LGXYZ,IIEQJJ,KKEQLL,IJEQKL,  &
                        IA,DAB,MAXNUM,DABMAX,CJAUX,CKAUX)
      USE PARCOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,INTENT(IN)::II,JJ,KK,LL
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KLOC,KKMIN,KKMAX
      INTEGER,DIMENSION(4,35),INTENT(IN)::IGXYZ,JGXYZ,KGXYZ,LGXYZ
      INTEGER,DIMENSION(NBF),INTENT(IN)::IA
      INTEGER,INTENT(IN)::MAXNUM
      LOGICAL,INTENT(IN)::IIEQJJ,KKEQLL,IJEQKL
      DOUBLE PRECISION,DIMENSION(NBFT,NBFT),INTENT(IN)::CJAUX,CKAUX
      DOUBLE PRECISION,INTENT(OUT)::DABMAX
      DOUBLE PRECISION,DIMENSION(84),INTENT(IN)::PNRM
      DOUBLE PRECISION,DIMENSION(MAXNUM),INTENT(OUT)::DAB
      DOUBLE PRECISION,PARAMETER::ZER=0.0D+00,PT5=0.5D+00
!-----------------------------------------------------------------------
      DABMAX=ZER         
!
!     GET 2e DENSITY FOR THIS SHELL BLOCK
!
      MINI= KKMIN(II)
      MINJ= KKMIN(JJ)
      MINK= KKMIN(KK)
      MINL= KKMIN(LL)
      LOCI= KLOC(II)-MINI
      LOCJ= KLOC(JJ)-MINJ
      LOCK= KLOC(KK)-MINK
      LOCL= KLOC(LL)-MINL
      MAXI= KKMAX(II)
      MAXJ= KKMAX(JJ)
      MAXK= KKMAX(KK)
      MAXL= KKMAX(LL)
         DO I=MINI,MAXI
            P1I = PNRM(I)
            JMAX= MAXJ
            IF(IIEQJJ) JMAX= I
            DO J=MINJ,JMAX
               P2J = P1I*PNRM(J)
               IAJ= MAX0(LOCI+I,LOCJ+J)
               IIJ= MIN0(LOCI+I,LOCJ+J)
               KMMAX=MAXK
               IF(IJEQKL) KMMAX= I
               DO K=MINK,KMMAX
                  P3K = P2J*PNRM(K)
                  LMAX= MAXL
                  IF(KKEQLL) LMAX= K
                  IF(IJEQKL .AND. K.EQ.I) LMAX= J
                  DO L=MINL,LMAX
                     P4L= P3K*PNRM(L)
                     KAL= MAX0(LOCK+K,LOCL+L)
                     KIL= MIN0(LOCK+K,LOCL+L)
                     IN = IAJ
                     JN = IIJ
                     KN = KAL
                     LN = KIL
                     IF(IN.LT.KN .OR.(IN.EQ.KN .AND. JN.LT.LN)) THEN
                        IN = KAL
                        JN = KIL
                        KN = IAJ
                        LN = IIJ
                     ENDIF
                     IJ = IA(IN)+JN
                     IK = IA(IN)+KN
                     IL = IA(IN)+LN
                     JK = IA(MAX0(JN,KN))+MIN0(JN,KN)
                     JL = IA(JN)+LN
                     IF(JN.LT.KN) JL = IA(MAX0(JN,LN))+MIN0(JN,LN)
                     KL = IA(KN)+LN
                     DF1 = ZER
                     DQ1 = ZER
                     DF1 = CJAUX(IJ,KL)
                     DQ1 = CKAUX(JK,IL) + CKAUX(JL,IK)
                     DF1 = DF1 + DF1 - DQ1
!                    AVOID DOUBLE COUNTING OF DIAGONAL TERMS                     
                     IF(JN.EQ.IN               ) DF1= DF1*PT5
                     IF(LN.EQ.KN               ) DF1= DF1*PT5
                     IF(KN.EQ.IN .AND. LN.EQ.JN) DF1= DF1*PT5
                     IF(DABMAX.LT. ABS(DF1)) DABMAX= ABS(DF1)
! IGXYZ AND J, K, AND L ARE SET UP IN JKDNDX
                     IJKL=IGXYZ(1,I)+JGXYZ(1,J)+KGXYZ(1,K)+LGXYZ(1,L)
                     DAB(IJKL)= DF1*P4L
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
!-----------------------------------------------------------------------
      RETURN
      END      
!-----------------------------------------------------------------------
