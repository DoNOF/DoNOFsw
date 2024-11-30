!======================================================================!
!                                                                      !
!                N A M E L I S T   S U B R O U T I N E S               !
!                                                                      !
!======================================================================!
!                                                                      !
!   NAMELIST_INPRUN: Specifies the run calculation, multiplicity, ...  !
!   NAMELIST_NOFINP: Preset values for the NOFINP namelist variables.  !
!   NAMELIST_DYNINP: Variables for Dynamics                            ! 
!                                                                      !
!======================================================================!

! NAMELIST_INPRUN
      SUBROUTINE NAMELIST_INPRUN(ITYPRUN,ICHARG,MULT,NINTEG,            &
        IDONTW,IEMOMENTS,EVECTOR,LIBRETA,                               &
        IECPO,IHSSCAL,IPROJECT,ISIGMA,                                  &
        NATmax,NSHELLmax,NPRIMImax,IHUBBARD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
#include "mpip.h"      
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
      COMMON/INFO/IUNTRD            
      COMMON/INTFIL/NINTMX           
      LOGICAL EFLDL                                                  
      COMMON/EFLDC_1/EFLDL
      COMMON/EFLDC_2/EVEC(3) 
      COMMON/ELPROP/IEMOM      
      CHARACTER(8):: UNITS
      COMMON/CONTROL/UNITS
      COMMON/INTOPT/ISCHWZ,IECP,NECP            
      COMMON/RUNTYPE/IRUNTYP
      COMMON/USELIBRETA/ILIBRETA
      COMMON/USEHUBBARD/IHUB
      LOGICAL DONTW,USELIB,USEHUB,HSSCAL,PROJECT
      DOUBLE PRECISION,DIMENSION(3) :: EVECTOR
      CHARACTER(6) :: RUNTYP,ENERGY,GRAD,OPTGEO,HESS,DYN
      DATA ENERGY,GRAD,OPTGEO/'ENERGY','GRAD  ','OPTGEO'/
      DATA HESS,DYN/'HESS  ','DYN   '/
      CHARACTER(8) :: ANGS,BOHR
      DATA ANGS, BOHR /'ANGS    ','BOHR    '/
      CHARACTER(4) :: ERITYP,FULL,RI,MIX
      DATA FULL,RI,MIX /'FULL','RI  ','MIX '/
      CHARACTER(5) :: RITYP
      LOGICAL SMCD
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      CHARACTER(3) :: GEN
!----------------------------------------------------------------------!
!                   --- INPRUN NAMELIST VARIABLES ---                  !
!----------------------------------------------------------------------!
!
! RUNTYP           specifies the run calculation
!       = ENERGY   1)single-point energy calculation (Default)
!       = GRAD     2)energy + gradients with respect to nuclear coord
!       = OPTGEO   3)optimize the molecular geometry
!       = HESS     4)compute numerical hessian from analytic gradients
!       = DYN      5)run Born-Oppenheimer on-the-fly molecular dynamics
!
! MULT             Multiplicity of the electronic state
!       = 1        singlet (Default)
!       = 2,3,...  doublet, triplet, and so on
!
! ICHARG           Molecular charge  
!       = 0        Neutral Molecule (Default)
!
! IECP             Effective Core Potentials 
!       = 0        (Default) All electron calculation 
!       = 1        Read ECP potentials in the $ECP group
!
! IEMOM            Electrostatic moments calculation
!       = 1        calculate dipole moments (Default)
!         2        also calculate quadrupole moments
!         3        also calculate octopole moments
!
! UNITS            Distance units (any angles must be in degrees)
!       = ANGS     Angstroms (Default)
!       = BOHR     Bohr atomic units
!
! EVEC             An array of the three x,y,z components of
!                  the applied electric field, in a.u.
!                  (1 a.u. = 1 Hartree/e*bohr = 5.1422082(15)d+11 V/m)
!       = 0.0D0    (Default)
!
! USELIB           Use Libreta open source library for ERI calculation
!       = F        HONDO Calculator (Default)
!
! USEHUB           Use Hubbard Model
!       = F        (Default)
!
! DONTW            Do not write 2e- integrals on the disk (Unit=1)
!       = T        (Default)
!
! ERITYP           Typ of ERIs used in calculations
!       = FULL     4c ERIs (Default)
!       = RI       3c/2c ERIs for Resolution of the Identity (RI) App.
!       = MIX      3c/2c ERIs for Resolution of the Identity (RI) App.
!                  once converged change to 4c ERIs (FULL)
!
! RITYP            Typ of Auxiliary Basis
!       = JKFIT    Read from jkfit files (default)
!       = GEN      Use Generative Auxiliary Basis
!
! GEN              Generative Auxiliary Basis to use in RI Approx.
!                  if ERITYP = RI. Values: A2,A2*,A3,A3*,A4,A4* 
!       = A2*      (Default)
!
! SMCD             Symmetric Modified Cholesky Decomposition for the 
!                  G matrix in the RI Approximation.
!       = F        (Default)
!      
! HSSCAL           Compute Hessian from analytic gradients and carry
!                  out normal mode vibrational analysis at st. point 
!                  if RUNTYP = OPTGEO (IRUNTYP=3)
!       = T        (Default)
!
! PROJECT          Project Hessian to eliminate rot/vib contaminants
!       = T        (Default)
!
! ISIGMA           Rotational symmetric number for thermochemistry
!       = 1        There is not a center of symmetry (Default)
!       = 2        There is a center of symmetry
!                  For more info see https://cccbdb.nist.gov/thermo.asp
!
!JFHLewYee: Changed NATOMS allowed dimension from 100 to 1000
! NATmax           Maximum Number of Atoms
!       = 1001     (Default)
!
! NSHELLmax        Maximum Number of Shells
!       = 600      (Default)
!
! NPRIMImax        Maximum Number of Gaussian Functions
!       = 2000     (Default)
!
!-----------------------------------------------------------------------
      NAMELIST/INPRUN/RUNTYP,MULT,ICHARG,IECP,IEMOM,UNITS,EVEC,USELIB,  &
      USEHUB,DONTW,ERITYP,RITYP,GEN,SMCD,HSSCAL,PROJECT,&
      ISIGMA,NATmax,NSHELLmax,NPRIMImax
!-----------------------------------------------------------------------
      TI = 0.0D0                                                           
      TX = 0.0D0                                                           
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Preset values to namelist variables
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RUNTYP    = ENERGY
      MULT      = 1                                                        
      ICHARG    = 0
      IECP      = 0
      IEMOM     = 1
      UNITS     = ANGS
      EVEC      = 0.0D0     ! EVEC(1,2,3)=0
      USELIB    = .FALSE.
      USEHUB    = .FALSE.
      DONTW     = .TRUE.
      ERITYP    = FULL
      RITYP     = 'JKFIT'
      GEN       = 'A2*'
      SMCD      = .FALSE.
      HSSCAL    = .TRUE.
      PROJECT   = .TRUE.
      ISIGMA    = 1
      NATmax    = 1001     
      NSHELLmax = 600  
      NPRIMImax = 2000 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                            Read Namelist
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#ifdef MPI
      CONTINUE
#else
      REWIND(5)
#endif
      READ(5,INPRUN,ERR=1,END=1)
!-----------------------------------------------------------------------
!     Determine IRUNTYP
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                   
      IF(RUNTYP==ENERGY)THEN
        IRUNTYP = 1
      ELSE IF(RUNTYP==GRAD)THEN
        IRUNTYP = 2
      ELSE IF(RUNTYP==OPTGEO)THEN
        IRUNTYP = 3
      ELSE IF(RUNTYP==HESS)THEN
        IRUNTYP = 4
      ELSE IF(RUNTYP==DYN)THEN
        IRUNTYP = 5
      END IF
      ITYPRUN = IRUNTYP
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Use Libreta open source library for ERI calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(USELIB)THEN
        ILIBRETA = 1
      ELSE
        ILIBRETA = 0
      ENDIF
        LIBRETA = ILIBRETA
!     Stop if USELIB to use NAG library (E04DKF,E04DGF,F11JEF)
      IF(ILIBRETA==1)THEN
        WRITE(6,*)                                                       !lib
        WRITE(6,*)'Stop: To use the Libreta library you need data.cpp', &!lib
        '      and uncomment calls to the relevant routines'   !lib
        CALL ABRT                                                        !lib
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Hubbard Model
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(USEHUB)THEN
#ifdef MPI      
        WRITE(6,'(/1X,40A,/)')'For parallel calculations uncomment !HUB' !HUB
        STOP                                                             !HUB
#endif       
        IHUBBARD = 1
        if(IRUNTYP/=1)then
          IRUNTYP = 1
          WRITE(6,'(/1X,39A,/)')'Note: RUNTYP has been changed to ENERGY'
        end if
!      Forced Options
        EVEC = 0.0D0                ! No electric field
        IECP = 0                    ! Effective Core Potentials
        IEMOM = 0                   ! Electrostatic Moments
      ELSE
        IHUBBARD = 0
      END IF
      IHUB = IHUBBARD
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Typ of ERIs used in calculations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ERITYP==FULL)THEN
        IERITYP = 1
      ELSE IF(ERITYP==RI)THEN
        IERITYP = 2
      ELSE IF(ERITYP==MIX)THEN       
        IERITYP = 3
        MIXSTATE = 1                           ! 1 = RI, 2 = FULL
      ELSE
        WRITE(6,6) 
        CALL ABRT
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Typ of RIs auxiliary basis
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(RITYP=="JKFIT")THEN
        IRITYP = 1 
      ELSE IF(RITYP(1:3)=="GEN")THEN
        IRITYP = 2
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Determine Star in Auxiliary Basis if required
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(IERITYP==2 .or. IERITYP==3)then
        READ(GEN(2:2),'(I1)')IGEN
        IF(GEN(3:3)=='*') THEN
          ISTAR = 1
        ELSE IF(GEN(3:3)==' ') THEN
          ISTAR = 0
        ELSE
          ISTAR = -1
        END IF
      end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute Hessian
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(HSSCAL)THEN
        IHSSCAL = 1
      ELSE 
        IHSSCAL = 0
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Project Hessian to eliminate rot/vib contaminants
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(PROJECT)THEN
        IPROJECT = 1
      ELSE 
        IPROJECT = 0
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                   
!     Do not write 2e- integrals on the disk (Unit=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                   
      IF(DONTW)THEN
        IDONTW = 1
      ELSE 
        IDONTW = 0
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                   
!     Electric Field
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF((EVEC(1)==0.0D0).and.(EVEC(2)==0.0D0).and.(EVEC(3)==0.0D0))THEN
        EFLDL = .FALSE.
      ELSE
        EFLDL = .TRUE.
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Integral Options (TRFOPT common block)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ISCHWZ  = 0         ! Schwarz inequality off 
      NINTMX = 15000
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
!     Errors in the Input Namelist
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(UNITS==ANGS)THEN
        IUNTRD = 1                
      ELSEIF(UNITS==BOHR)THEN
        IUNTRD = -1
      ELSE
        WRITE(6,2)'UNITS ',UNITS                      
        CALL ABRT                                                      
      ENDIF                                                            
      IF(IECP<0.or.IECP>3)THEN
        WRITE(6,3)IECP                        
        CALL ABRT                                                      
      END IF
      IF(IHUB==0.and.(IEMOM<1.or.IEMOM>3))THEN   ! IEMOM = 1,2,3
        WRITE(6,4)IEMOM                       
        CALL ABRT                                                      
      END IF
      IF(IERITYP==2 .and. IDONTW==0)THEN   ! IDONTW=1 with ERITYP=RI
        WRITE(6,5)
        CALL ABRT
      END IF
      IF(IECP>0 .and. ILIBRETA==1)THEN        ! ECP only with HONDO
        WRITE(6,7)
        CALL ABRT
      END IF
      !JFHLewYee: Changed NATOMS allowed dimension from 100 to 1000
      IF(IECP>0 .and. NATmax>1001)THEN         ! due to COMMON/ECP2/
        WRITE(6,8)
        CALL ABRT
      END IF
      IF(IECP>0 .and. NSHELLmax>600)THEN      ! due to MAPSHL
        WRITE(6,9)
        CALL ABRT
      END IF
        IF(ILIBRETA==1 .and. NSHELLmax>600)THEN ! due to BASLIB
        WRITE(6,10)
        CALL ABRT
      END IF
      IF(IERITYP==2 .and. NSHELLmax>700)THEN  ! due to COMMON/NSHELaux/
        WRITE(6,11)
        CALL ABRT
        END IF
      IF(IERITYP==2 .and. NPRIMImax>2000)THEN ! due to COMMON/EXCaux/
        WRITE(6,12)
        CALL ABRT
      END IF
      !JFHLewYee: Changed NATOMS allowed dimension from 100 to 1000
      IF(IRUNTYP==5 .and. NATmax>1001)THEN ! due to COMMON/INPDYN_FLAGS/
        WRITE(6,13)
        CALL ABRT
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ICH = ICHARG                                                      
      MUL = MULT
      IEMOMENTS = IEMOM
      EVECTOR = EVEC
      NINTEG = NINTMX
      IECPO = IECP                                                   
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                        
      RETURN
!-----------------------------------------------------------------------
!     Format definitions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    1 WRITE(6,'(/2X,36A)')'Stop: Wrong INPRUN Namelist Variable'
      STOP
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    2 FORMAT(/1X,'Stop: $CONTRL KEYWORD ',A6,                           &
      ' was given an illegal value ',A8,'.'/)
    3 FORMAT(/1X,'Stop: IECP   must be between  0 and 3, not',I8,/)
    4 FORMAT(/1X,'Stop: IEMOM  must be between  1 and 3, not',I8,/) 
    5 FORMAT(/1X,'Stop: DONTW must be T with ERITYP = RI',/)
    6 FORMAT(/1X,'Stop: ERITYP must be FULL or RI',/)
    7 FORMAT(/1X,'Stop: USELIB must be F with ECP',/) 
      !JFHLewYee: Changed NATOMS allowed dimension from 100 to 1000
    8 FORMAT(/1X,'Stop: NATmax must be <= 1001 with ECP',/) 
    9 FORMAT(/1X,'Stop: NHELLmax must be <= 600 with ECP',/)
   10 FORMAT(/1X,'Stop: NHELLmax must be <= 600 with USELIB = T',/)
   11 FORMAT(/1X,'Stop: NHELLmax must be <= 700 with ERITYP = RI',/)
   12 FORMAT(/1X,'Stop: NPRIMImax must be <= 2000 with ERITYP = RI',/)
      !JFHLewYee: Changed NATOMS allowed dimension from 100 to 1000
   13 FORMAT(/1X,'Stop: NATmax must be <= 1001 with RUNTYP = DYN',/)
!-----------------------------------------------------------------------
      END

! NAMELIST_NOFINP                                                     
      SUBROUTINE NAMELIST_NOFINP(IRUNTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HFDAMP,HFEXTRAP,HFDIIS,DIIS,PERDIIS,DAMPING,EXTRAP
      LOGICAL RESTART,FROZEN,PRINTLAG,DIAGLAG,APSG,CHKORTHO,ORTHO
      LOGICAL ERPA,OIMP2,MBPT,SC2MCPT,HFID,HighSpin,SCALING,RHF
      LOGICAL AUTOLR
      DOUBLE PRECISION LR,FACT,BETA1,BETA2
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_RHF/IRHFTYP,NCONVRHF,CONVRHFDM,MAXITRHF,RHF
      COMMON/INPSCALING/SCALING,NZEROS,NZEROSm,NZEROSr,ITZITER
      COMMON/INPADAM/LR,AUTOLR,FACT,BETA1,BETA2
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_HFCONVTECH/HFDAMP,HFEXTRAP,HFDIIS
      COMMON/INPNOF_DIIS/DIIS,PERDIIS,NDIIS,NTHDIIS,THDIIS
      COMMON/INPNOF_DAMPEXTRAP/DAMPING,EXTRAP
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_FROZEN/FROZEN,IFROZEN(200)
      COMMON/INPNOF_LAGRANGE/PRINTLAG,DIAGLAG
      COMMON/INPNOF_APSG/APSG,NTHAPSG,THAPSG
      COMMON/INPNOF_ORTHOGONALITY/CHKORTHO,ORTHO
      COMMON/INPNOF_ERPA/ERPA
      COMMON/INPNOF_OIMP2/OIMP2,MBPT
      COMMON/INPNOF_SC2MCPT/SC2MCPT
      COMMON/INPNOF_HFID/HFID,NTHRESHEID,THRESHEID,MAXITID,KOOPMANS       
      COMMON/INPNOF_EINI/IEINI
      COMMON/INPNOF_EKT/IEKT
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/INPNOF_INICON/INICOND
      COMMON/INPNOF_ARDM/THRESHDM,NOUTRDM,NSQT,NTHRESHDM      
      COMMON/INPNOF_CJK/NOUTCJK,NTHRESHCJK,THRESHCJK      
      COMMON/INPNOF_Tijab/NOUTTijab,NTHRESHTijab,THRESHTijab
      COMMON/INPNOF_CGM/ICGMETHOD
      COMMON/INPNOF_NTHRESH/NTHRESHL,NTHRESHE,NTHRESHEC,NTHRESHEN
      COMMON/INPNOF_THRESH/THRESHL,THRESHE,THRESHEC,THRESHEN
      COMMON/INPNOF_COEFOPT/MAXLOOP
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPNOF_STATIC/Ista
      COMMON/INPNOF_EXSTA/NESt,OMEGA1
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/INPFILE_NO1PT2/NO1PT2,NEX
      COMMON/USEHUBBARD/IHUB
      LOGICAL SMCD
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      !
      LOGICAL AUTOZEROS
      INTEGER:: IRUNTYP
!----------------------------------------------------------------------!
!                   --- NOFINP NAMELIST VARIABLES ---                  !
!----------------------------------------------------------------------!
!
!.......... MAXIT               Maximum number of OCC-SCF iterations 
!                      = 1000   (Default)
!
!.......... MAXIT21             Maximum number of iterations for 
!                               optimizing only by NOs keeping fixed ONs 
!                               if ICOEF = 21 (first part of the Opt.)
!                      = 10     (Default)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Type of Calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.......... ICOEF               Energy Optimization with respect to NOs
!
!                      = 0      Optimize only with respect to ONs
!                      = 1      Optimize by the ONs and NOs (Default)
!                      = 2      Optimize only by NOs keeping fixed ONs
!                      = 21     ICOEF=2 with MAXIT21 and then ICOEF=1
!                      = 3      Optimize by all ONs and core-fragment 
!                               orbitals. The rest of fragment orbitals 
!                               remain frozen
!
!.......... ISOFTMAX            Use Softmax function for ON (Gamma) opt.
!                      = 1      (Default)
!
!.......... IORBOPT             Select method for NO optimization
!
!                      = 1      Iterative diagonalization (OrbOptFMIUGr)
!                      = 2      By unitary tranformations (OrbOptRot)
!                      = 3      Sequential Quadratic Program (OrbOptSQP)
!                      = 4      Adaptative Momentum (ADAM)             
!                               (Default)
!                      = 5      Gradient Descent (GD)                                  
!                      = 6      Momentum Gradient Descent (MOMENTUMGD) 
!                      = 7      Nesterov Gradient Descent (NGD)        
!                      = 8      Root Mean Square Propagation (RMSProp) 
!                      = 9      Adaptative Learning Rate (ADADELTA)    
!                      = 10     Adaptative Subgradient (ADAGRAD)       
!                      = 11     Decaying Momentum (DEMON)              
!
!.......... IEINI               Calculate only the initial energy
!                      = 0      (Default)
!
!.......... NO1                 Max. index of NOs with Occupation = 1
!                      = -1     Consider Core NOs
!                      = 0      All NOs are considered (Default)
!                      = Value  User specifies how many NOs have OCC.=1
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Hartree-Fock
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.......... RHF                 Restricted Hartree-Fock Calculation
!                      = T      (Default)
!
!.......... NCONVRHF            RHF-SCF Density Convergence Criteria
!                               CONVRHFDM=10.0**(-NCONVRHF)
!                      = 5      (Default)
!
!.......... MAXITRHF            Maximum number of RHF-SCF iterations 
!                      = 100    (Default)
!
!.......... HFDAMP              Damping of the Fock matrix
!                      = T      (Default)
!
!.......... HFEXTRAP            Extrapolation of the Fock matrix
!                      = T      (Default)
!
!.......... HFDIIS              Direct Inversion in the Iterative 
!                               Subspace in the RHF-SCF optimization
!                      = T      (Default)
!
!.......... HFID                Use the Iterative Diagonalization Method 
!                               to generate the HF Orbitals
!                      = F      (Default)
!
!.......... NTHRESHEID          Convergence of the TOTAL ENERGY
!                               THRESHEID=10.0**(-NTHRESHEID)
!                      = 5      (Default)
!
!.......... MAXITID             Maximum number of external iterations 
!                      = 100    (Default)
!
!.......... KOOPMANS            Calculate IPs using Koopmans' Theorem 
!                      = 0      (Default)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! PNOF Selection
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.......... IPNOF               Type of Natural Orbital Functional (NOF)
!                      = 5      PNOF5
!                      = 6      PNOF6
!                      = 7      PNOF7
!                      = 8      GNOF (Default)
!
!.......... Ista                Use Static version of PNOF7 
!                      = 0      PNOF7 (Default)
!                      = 1      PNOF7s
!
!.......... HighSpin            Spin-uncompensated calculation type
!                      = F      (Default) Multiplet state (Ms=0)
!                      = T      High-spin uncompensated state (Ms=S)
!
!.......... NCWO                Number of coupled weakly occupied MOs 
!                               per strongly occupied = Nc -> PNOFi(Nc)
!                      = 1      NCWO = 1
!                      = 2,3,...
!                      =-1      NCWO = NVIR/NDOC (Default)
!                               NVIR: Number of HF virtual  MOs (OCC=0)
!                               NDOC: Number of strongly occupied MOs
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Convergence Criteria in NOF calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.......... NTHRESHL            Convergence of the Lagrange Multipliers
!                               THRESHL=10.0**(-NTHRESHL)
!                      = 3      (Default)
!
!.......... NTHRESHE            Convergence of the total energy
!                               THRESHE=10.0**(-NTHRESHE)
!                      = 4      (Default)
!
!.......... NTHRESHEC           Convergence of the total energy (OrbOpt)
!                               THRESHEC=10.0**(-NTHRESHEC)
!                      = 8      (Default)
!
!.......... NTHRESHEN           Convergence of the total energy (OccOpt)
!                               THRESHEN=10.0**(-NTHRESHEN)
!                      = 10     (Default)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Options for the Orbital Optimization Program (ID Method)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.......... MAXLOOP             Maximum Iteration Number for the SCF-
!                               iteration cycle in each ITCALLs 
!                      = 30     (Default)
!
!     The straightforward iterative scheme fails to converge very 
!     often due to the values of some off-diagonal elements Fki. The 
!     latters must be suffciently small and of the same order of 
!     magnitude. A variable factor scales Fki. We establish an upper
!     bound B, in such a way that when the absolute value of the 
!     matrix element Fki is greater than B, it is scaled by a factor 
!     Cki (F'ki = Cki*Fki ), as to satisfy ABS(Fki) <= B.
!
!.......... SCALING             A variable factor scales Fki
!                      = T      (Default)
!
!.......... NZEROS              B = 10.0**(1-NZEROS). 
!                               Initial number of ZEROS in Fij. The 
!                               scaling factor varies until the number 
!                               of ZEROS (.000##) is equal for all 
!                               elements Fij.
!                      = 0      (Default)
!
!.......... NZEROSm             B = 10.0**(1-NZEROSm)
!                               Maximum number of zeros in Fij.
!                      = 5      (Default)
!
!.......... NZEROSr             B = 10.0**(1-NZEROSr)
!                               Number of zeros in Fij to restart 
!                               automatically the calculation.
!                      = 2      (Default)
!
!.......... AUTOZEROS           The code select automatically values
!                               for NZEROS,NZEROSm & NZEROSr
!                               Note: Override previously selected values
!                      = T      (Default)
!
!.......... ITZITER             Number of Iterations for constant scaling
!                      = 10     (Default)
!
!.......... DIIS                Direct Inversion in the Iterative 
!                               Subspace in the orbital optimization if 
!                               DUMEL < THDIIS every NDIIS loops
!                      = T      (Default)
!
!.......... NTHDIIS             Energy threshold to begin DIIS
!                      = 3      THDIIS = 10.0**(-NTHDIIS) (Default)
!
!.......... NDIIS               Number of considered loops to interpolate
!                               the generalized Fock matrix in the DIIS
!                      = 5      (Default)
!
!.......... PERDIIS             Periodic DIIS
!                      = T      Apply DIIS every NDIIS (Default)
!                      = F      DIIS is always applied after NDIIS
!
!.......... DAMPING             Damping of the Gen. Fock matrix (FMIUG)
!                      = F      (Default)
!
!.......... EXTRAP              Extrapolation of the Gen. Fock matrix
!                      = F      (Default)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Options for the Orbital Optimization Program (ADAM Method)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.......... LR                  Initial Learning Rate
!                      = 0.01   (Default)
!
!.......... FACT                Factor to multiply Learning Rate
!                               to reduce it after a failed ADAM     
!                      = 0.4    (Default)
!
!.......... BETA1               BETA1 Memory for first momentum
!                      = 0.7    (Default)
!
!.......... BETA2               BETA2 Memory for second momentum
!                      = 0.9    (Default)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Options for pertubative calculations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.......... ERPA                Extended Random Phase Approximation
!                     = F       (Default)
!
!.......... OIMP2               NOF - Orbital Invariant MP2
!                     = F       (Default)
!
!.......... MBPT                NOF - c - X (X=RPA, MP2, etc)
!                     = F       (Default)
!
!.......... NO1PT2              Frozen MOs in perturbative calculations
!                               Maximum index of NOs with Occupation = 1
!                      = -1     = NO1 (Default)
!                      = 0      All NOs are considered
!                      = Value  User specifies how many NOs are frozen
!
!.......... SC2MCPT             SC2-MCPT perturbation theory is used to
!                               correct the PNOF5 Energy. 
!                               2 outputs: PNOF5-SC2-MCPT and PNOF5-PT2
!                     = F       (Default)
!
!.......... NEX                 Number of excluded coupled orbitals 
!                               in the PNOF5-PT2 calculation
!                      = 0      All NOs are included (Default)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Input Options for Gamma (Occ), C and Diagonal F
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!                     --- Restart Options ---
!
!.......... RESTART             Restart from GCF file (Default=F)
!                      = F      INPUTGAMMA=0,INPUTC=0,INPUTFMIUG=0
!                      = T      INPUTGAMMA=1,INPUTC=1,INPUTFMIUG=1
!
!.......... INPUTGAMMA          Guess for GAMMA variables (ONs)
!                      = 0      Close Fermi-Dirac Distribution (Default)
!                      = 1      Input from file GCF
!
!.......... INPUTC              Guess for Coefficient matrix (NOs)
!                      = 0      Use HCORE or HF Eigenvectors (Default)
!                      = 1      Input from file GCF
!
!.......... INPUTFMIUG          Guess for Diagonal elements (FMIUG0)
!                      = 0      Use single diag. of Lagragian (Default)
!                      = 1      Input from file GCF
!
!.......... INPUTCXYZ           Nuclear Coordinates (CXYZ)
!                      = 0      Input from input file (*.inp) (Default)
!                      = 1      Input from file GCF (only if RESTART=T)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Output Options
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.......... NPRINT              Output Option (Default VALUE: 0)
!                      = 0      Short Printing (Occ,Emom,Energies)
!                      = 1      Output at initial and final iterations
!                               including MOs,Pop,APSG,Lag,IPs,DMs,CJK
!                      = 2      Output everything in each iteration
!
!.......... IWRITEC             Output Option for the Coefficient matrix
!                      = 0      No output (Default)
!                      = 1      Output the Coefficient Matrix (NOs)
!
!.......... IMULPOP             Mulliken Population Analysis
!                      = 0      Not do it (Default)
!                      = 1      Do it 
!
!.......... PRINTLAG            Output Option for Lagrange Multipliers
!                      = F      No Output (Default)
!
!.......... DIAGLAG             Diagonalize Lagrange Multipliers
!                               Print new 1e- Energies, Canonical MOs, 
!                               and new diagonal elements of the 1RDM
!                      = F      Not do it (Default)
!
!.......... IEKT                IPs by Ext. Koopmans' Theorem (EKT)
!                      = 0      Not calculate the IPs
!                      = 1      Calculate ionization potentials (IPs) 
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.......... IAIMPAC             Write information into WFN file 
!                               for the AIMPAC PROGRAM (UNIT 7)
!                      = 0      Don't write
!                      = 1      Write into WFN file (Default)
!
!.......... IFCHK               Write information into Formatted 
!                               Checkpoint (FCHK) file for visualization 
!                               software (UNIT 19)
!                      = 0      Don't write
!                      = 1      Write into FCHK file (Default)
!
!.......... MOLDEN              Write information into MLD file 
!                               for the MOLDEN PROGRAM (UNIT 17)
!                      = 0      Don't write 
!                      = 1      Write into MLD file (Default)
!
!.......... INICOND             Create initial conditions file (UNIT 33)
!                               according to normal modes with equil.
!                               geometry & ZPE velocities in Ang/fs
!                               for MOLDYN routine (IRUNTYP=3)
!                      = 0      Don't create
!                      = 1      Create ini.xyz file (Default)
!
!.......... NOUTRDM             Print OPTION for atomic RDMs 
!                      = 0      NO Output (Default)
!                      = 1      Print atomic 1-RDM in 1DM file
!                      = 2      Print atomic 2-RDM in 2DM file
!                      = 3      Print atomic RDMs in 1DM and 2DM files
!
!.......... NTHRESHDM           THRESHDM=10.0**(-NTHRESHDM)
!                      = 6      (Default)
!
!.......... NSQT                Print OPTION for 2DM file
!                      = 1      UNforMATTED (Default)
!                      = 0      forMATTED (SEE SUBROUTINE OUTPUTRDMrc)
!
!.......... NOUTCJK             Print OPTION for CJ12 and CK12
!                      = 0      NO Output (Default)
!                      = 1      Print CJ12 and CK12 in file 'CJK'
!
!.......... NTHRESHCJK          THRESHCJK=10.0**(-NTHRESHCJK)
!                      = 6      (Default)
!
!.......... NOUTTijab           Print OPTION for Tijab
!                      = 0      NO Output (Default)
!                      = 1      Print Tijab in file 'Tijab'
!
!.......... NTHRESHTijab        THRESHTijab=10.0**(-NTHRESHTijab)
!                      = 6      (Default)
!
!.......... APSG                Open an APSG file for printing the 
!                               coefficient matrix ($VEC-$END) and the
!                               expansion coefficients of the APSG
!                               generating wavefunction.
!                      = F      Output (Default)
!
!.......... NTHAPSG             Threshold for APSG expansion coefficient
!                               THAPSG = 10.0**(-NTHAPSG)
!                      = 10     (Default)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Optional Options
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.......... ORTHO               Orthogonalize the initial orbitals
!                      = F      No 
!                      = T      Yes (Default)
!
!.......... CHKORTHO            Check the Orthonormality of the MOs
!                      = F      No (Default)
!                      = T      Yes
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Options related to Frozen coordinates in gradient computation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.......... FROZEN              Is there any fixed coordinate
!                      = F      (Default)
!
!.......... IFROZEN             By pairs, what coordinate of which atom,
!                               e.g. 2,5,1,1 means "y" coordinate of
!                               atom 5 and "x" coor of atom 1 to freeze.
!                               MAXIMUM of frozen coordinates = 10
!                      = 0      (Default)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.......... ICGMETHOD           Define the Conjugate Gradient Method in
!                               OccOpt, CALTijabIsym and OPTIMIZE
!                      = 1      SUMSL: CGOCUPSUMSL,OPTSUMSL
!                               SparseSymLinearSystem_CG (Default)
!                      = 2      Use NAG subroutines:
!                               E04DGF: OPTCGNAG,CGOCUPNAG 
!                               F11JEF: SparseSymLinearSystem_NAG
!                      = 3      LBFGS: OPTLBFGS, LBFGSOCUP
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Options for Excited States
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.......... NESt                Number of Excited States considered in
!                               the ensemble NOFT W = (w1,w2,w3,...,wd)
!                      = 0      (Default)
!
!.......... OMEGA1              Value for w1 in the ensemble
!                      = 1.0    Only GS is considered (Default) 
!
!-----------------------------------------------------------------------
      NAMELIST/NOFINP/MAXIT,MAXIT21,ICOEF,ISOFTMAX,IORBOPT,IEINI,NO1,   &
                      RHF,NCONVRHF,MAXITRHF,HFDAMP,HFEXTRAP,HFDIIS,HFID,&
                      NTHRESHEID,MAXITID,KOOPMANS,IPNOF,Ista,HighSpin,  &
                      NCWO,NTHRESHL,NTHRESHE,NTHRESHEC,NTHRESHEN,       &
                      MAXLOOP,SCALING,AUTOZEROS,NZEROS,NZEROSm,NZEROSr, &
                      ITZITER,DIIS,NTHDIIS,NDIIS,PERDIIS,DAMPING,       &
                      EXTRAP,ERPA,OIMP2,NO1PT2,SC2MCPT,NEX,RESTART,     &
                      INPUTGAMMA,INPUTC,INPUTFMIUG,NPRINT,INPUTCXYZ,    &
                      IWRITEC,IMULPOP,MBPT,PRINTLAG,DIAGLAG,IEKT,       &
                      IAIMPAC,IFCHK,MOLDEN,INICOND,NOUTRDM,NTHRESHDM,   &
                      NSQT,NOUTCJK,NTHRESHCJK,NOUTTijab,NTHRESHTijab,   &
                      APSG,NTHAPSG,ORTHO,CHKORTHO,FROZEN,IFROZEN,       &
                      ICGMETHOD,NESt,OMEGA1,LR,FACT,BETA1,BETA2
!-----------------------------------------------------------------------
!     Preset values to namelist variables
!-----------------------------------------------------------------------
      MAXIT=1000
      MAXIT21=10
      
!     Type of Calculation
      ICOEF=1
      ISOFTMAX=1
      IORBOPT=4
      IEINI=0
      NO1=0
      
!     Hartree-Fock
      RHF=.TRUE.                                  ! AO Basis
      NCONVRHF=5
      MAXITRHF=100
      HFDAMP=.TRUE.
      HFEXTRAP=.TRUE.
      HFDIIS=.TRUE.
      HFID=.FALSE.                                ! MO Basis
!     For HFID=T.and.RHF=F, it may be better to take HFDAMP=HFEXTRAP=F 
      NTHRESHEID=8
      MAXITID=100
      KOOPMANS=0
      
!     PNOF Selection
      IPNOF=8
      Ista=0                                      ! PNOF7n            
      HighSpin=.FALSE.                            ! Multiplet      
      NCWO=-1
      
!     Convergence Criteria in NOF calculation
      NTHRESHL=3
      NTHRESHE=4
      NTHRESHEC=8
      NTHRESHEN=10
      
!     Options for the Orbital Optimization Program (ID Method)
      MAXLOOP=30
      SCALING=.TRUE.
      NZEROS=0
      NZEROSm=5
      NZEROSr=2      
      AUTOZEROS=.TRUE.      
      ITZITER=10
      DIIS=.TRUE.
      NTHDIIS=3
      NDIIS=5
      PERDIIS=.TRUE.
      DAMPING=.FALSE.
      EXTRAP=.FALSE.
      
!     Options for pertubative calculations
      OIMP2=.FALSE.
      MBPT=.FALSE.
      ERPA=.FALSE.
      NO1PT2=-1                  
      SC2MCPT=.FALSE.
      NEX=0      
      
!     Input Options for Gamma (Occ), C and Diagonal F
      RESTART=.FALSE.
      INPUTGAMMA=0
      INPUTC=0
      INPUTFMIUG=0
      INPUTCXYZ=0
      
!     Output Options
      NPRINT=0
      
!     for NPRINT>0
      IWRITEC=0
      IMULPOP=0
      PRINTLAG=.FALSE.
      DIAGLAG=.FALSE.  ! Only in the final Output
      IEKT=0           ! Only in the final Output
      !
      IAIMPAC=1
      IFCHK=1
      MOLDEN=1
      INICOND=1
      !
      NOUTRDM=0
      NTHRESHDM=6
      NSQT=1
!      
      NOUTCJK=0
      NTHRESHCJK=6
      NOUTTijab=0
      NTHRESHTijab=6
      !
      APSG=.FALSE.
      NTHAPSG=10
      
!     Optional Options
      ORTHO=.TRUE.
      CHKORTHO=.FALSE.
      
!     Frozen coordinates
      FROZEN=.FALSE.
      IFROZEN=0
      
!     Options for the Conjugate Gradient Method
      ICGMETHOD=1
      
!     Options for the ADAM Method
      AUTOLR=.TRUE.
      LR=0.01D0
      FACT=0.2D0
      BETA1=0.7D0
      BETA2=0.9D0

!     Excited States
      NESt=0
      OMEGA1=1.0d0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Read namelist variables
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      CONTINUE
#else
      REWIND(5)
#endif
      READ(5,NOFINP,END=1,ERR=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Note: Override AUTOZEROS to False and reduce MAXLOOP for ADAM
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IORBOPT>1) THEN
        IF(AUTOZEROS) THEN
          AUTOZEROS=.FALSE.
        END IF
      END IF
      IF(IORBOPT==4 .AND. MAXLOOP==30) THEN
        MAXLOOP=10
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Automatic selection of NZEROS                
!     Note: Override the selected values in the namelist
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(AUTOZEROS)then
        IF(RESTART)then       ! RESTART=T
          if(NTHRESHL<=3)THEN
            NZEROS  = 2
            NZEROSr = 2
            NZEROSm = 5
          else
            NZEROS  = NTHRESHL - 1
            NZEROSr = NZEROS
            NZEROSm = NTHRESHL + 2
          end if
        ELSE                  ! RESTART=F
          if(IRUNTYP/=5)then
            if(NTHRESHL>3)then
              NTHRESHL=3
              WRITE(6,'(/1X,79A/)')'Beware of high NTHRESHL with RESTART=F, NTHRESHL has been reduced to 3 !!!'
            end if
            NZEROS  = 1
            NZEROSr = 2
            NZEROSm = NTHRESHL + 2
          end if
        END IF
      end if      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NVWO=-1 if NE=2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NE==2.and.NCWO/=-1)THEN
        WRITE(6,3)
        NCWO=-1
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Force THRESHEC=1.0d-08 if larger than this value
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NTHRESHEC<8)THEN
        NTHRESHEC = 8
        WRITE(6,'(/1X,57A/)')'NTHRESHEC is very small so it has been increased to 8 !!!'
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Convergence Criteria
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CONVRHFDM   = 10.0**(-NCONVRHF) 
      THRESHEID   = 10.0**(-NTHRESHEID)
      THRESHL     = 10.0**(-NTHRESHL)
      THRESHE     = 10.0**(-NTHRESHE)
      THRESHEC    = 10.0**(-NTHRESHEC)
      THRESHEN    = 10.0**(-NTHRESHEN)
      THRESHDM    = 10.0**(-NTHRESHDM)
      THRESHCJK   = 10.0**(-NTHRESHCJK)
      THRESHTijab = 10.0**(-NTHRESHTijab)      
      THDIIS      = 10.0**(-NTHDIIS)
      THAPSG      = 10.0**(-NTHAPSG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Mandatory Options with RUNTYP=OPTGEO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IRUNTYP==3)THEN
        if(HFID.or.RHF)then
          HFID = .FALSE.
          RHF = .FALSE.        
          WRITE(6,'(/,1X,35A,/)')'!OPTGEO: HF has been set equal FALSE'
        end if
        if(OIMP2 .OR. MBPT .OR. ERPA)then      
          OIMP2 = .FALSE.
          MBPT = .FALSE.
          ERPA = .FALSE.
          NOUTTijab = 0
          WRITE(6,'(/,1X,38A,/)')'!OPTGEO: OIMP2 has been set equal FALSE'
        end if
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Mandatory Options with RUNTYP=HESS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IRUNTYP==4)THEN
        if(HFID.or.RHF)then
          HFID = .FALSE.
          RHF = .FALSE.        
          WRITE(6,'(/,1X,39A,/)')'!HESS cal.: HF has been set equal FALSE'
        end if
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Mandatory Options with RUNTYP=DYN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IRUNTYP==5)THEN
        if(HFID.or.RHF)then
          HFID = .FALSE.
          RHF = .FALSE.
          WRITE(6,'(/,1X,33A,/)')'!DYN: HF has been set equal FALSE'
        end if
        if(OIMP2)then
          OIMP2 = .FALSE.
          NOUTTijab = 0
          WRITE(6,'(/,1X,36A,/)')'!DYN: OIMP2 has been set equal FALSE'
        end if
!
        RESTART = .FALSE.
        INPUTCXYZ = 0       
        WRITE(6,'(/,1X,38A,/)')'!DYN: RESTART=FALSE, INPUTCXYZ=0'              
        INPUTGAMMA = 1
        INPUTC = 1
        INPUTFMIUG = 1
        WRITE(6,'(  1X,38A,/)')'!INPUTGAMMA=1, INPUTC=1, INPUTFMIUG=1'
        IAIMPAC=0
        MOLDEN=0
        INICOND=0
        IFCHK=0
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Restart Options
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(RESTART)THEN
        INPUTGAMMA=1
        INPUTC=1
        INPUTFMIUG=1
        INPUTCXYZ=1
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Stop if ICGMETHOD=2 to use NAG library (E04DKF,E04DGF,F11JEF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ICGMETHOD==2)THEN
        WRITE(6,*)                                                       
        WRITE(6,*)'Stop: To use the NAG library you must uncomment',&   !nag-stop 
        '      the calls to the relevant routines'                      !nag-stop
        CALL ABRT                                                       !nag-stop 
        NTHRESHEN = 16                                               
        THRESHEN = 10.0**(-NTHRESHEN)                                
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Stop if ICOEF=2 to evaluate only the initial energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF( (ICOEF==2.or.ICOEF==21) .and. IEINI==1)THEN
        WRITE(6,*)                                                       &
        'STOP: Choose ICOEF/=2,21 to evaluate the initial energy'
        CALL ABRT       
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Stop if ICOEF==3 and RI Approximation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ICOEF==3. .and. IERITYP==2)THEN
        WRITE(6,*)'Sorry: Fragment Calc. is not implemented with RI App.'
        CALL ABRT            
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     STOP if ICOEF /= 0,1,2,21,3
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(.not. (ICOEF== 0 .or. ICOEF==1 .or. ICOEF==2                   &
        .or. ICOEF==21 .or. ICOEF==3) )THEN
        WRITE(6,*)'STOP: ICOEF must be 0,1,2,21 or 3'
        CALL ABRT       
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     STOP options if IPNOF/=5 and APSG or SC2MCPT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IPNOF/=5)THEN
        IF(APSG)THEN
          WRITE(6,*)'STOP: APSG=T is only valid for PNOF5'
          CALL ABRT       
        ENDIF
        IF(SC2MCPT)THEN
          WRITE(6,*)'STOP: SC2-MCPT is only valid for PNOF5'
          CALL ABRT       
        ENDIF
      ENDIF
!
      IF(IPNOF==8)Ista=0     ! FIs = DSQRT(RO*(1.0d0-RO))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Forced options if Hubbard model is used
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IHUB==1)THEN
        NO1 = 0               ! NOs with Occupation = 1
        NO1PT2 = 0            ! Frozen MOs in perturbative calculations
        NCWO = 1              ! Number of coupled weakly occupied MOs 
        IMULPOP = 0           ! Mulliken Population Analysis
        IAIMPAC = 0           ! Write information into WFN file
        MOLDEN = 0            ! Write information into MLD file
        INICOND = 0           ! Create ini.xyz file
        IFCHK = 0             ! Write information into FCHK file
        INPUTCXYZ = 0         ! No Nuclear Coordinates (CXYZ)
!      HFID = .TRUE.         ! HF Iterative Diagonalization Method 
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
!-----------------------------------------------------------------------
!     Namelist Stop errors
!-----------------------------------------------------------------------
    1 WRITE(6,2)
      STOP
!-----------------------------------------------------------------------
!     Format definitions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    2 FORMAT(/2X,'**************************************',              &
      /2X,'*                                    *',                     &
      /2X,'*            SORRY,                  *',                     &
      /2X,'*   ERROR IN NAMELIST PARAMETERS     *',                     &
      /2X,'*                                    *',                     &
      /2X,'**************************************')               
    3 FORMAT(/1X,'!!! Warning: In the case of two electrons,'           &
      /5X,'there is only one electron pair.',                           &
      /5X,'To avoid spurious interpair contributions',                  &
      /5X,'NCWO has been set equal to -1 !!! ')
!-----------------------------------------------------------------------
      END

! NAMELIST_DYNINP                                                     
      SUBROUTINE NAMELIST_DYNINP(nat,resflag,velflag,dt,tmax,ngcf)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER :: nat,ngcf
       CHARACTER(LEN=1) :: dflag,resflag,velflag,snapshot
       DOUBLE PRECISION :: dt,tmax,Vxyz
      !JFHLewYee: Changed NATOMS allowed dimension from 100 to 1000
       COMMON/INPDYN_FLAGS/dflag(3,1001)
       COMMON/INPDYN_VELOCITY/Vxyz(3,1001)
       COMMON/INPDYN_NSHOT/snapshot
!----------------------------------------------------------------------!
!                   --- DYNINP NAMELIST VARIABLES ---                  !
!----------------------------------------------------------------------!
!                             
!...... velflag               Flag to make constant velocity dynamics
!                      = F    (Default)
!                             
!...... dt                    Time step in fs
!                      = 1    (Default)
!                             
!...... tmax                  Propagation time in fs
!                      = 100  (Default)
!                             
!...... ngcf                  number of GCF files to use in calculation
!                      = 1    (Default)
!
!...... Vxyz           = 0.0  Initial velocities per atom
!
!...... resflag               Restart MD calculation
!                      = F    (Default)
!
!...... snapshot              Save MLD file in snapshot-t.mld
!                      = F    (Default)
!                             
!-----------------------------------------------------------------------
      NAMELIST/INPDYN/velflag,dt,tmax,ngcf,Vxyz,resflag,snapshot
!-----------------------------------------------------------------------
      dflag(1:3,1:nat) = 'T'  ! Flags for selective dynamics
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                      
!     Preset values to namelist variables
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      velflag = 'F'
      dt = 1.0d0    
      tmax = 100.0d0  
      ngcf = 1 
      Vxyz(1:3,1:nat) = 0.0d0
      resflag = 'F'
      snapshot = 'F'
!-----------------------------------------------------------------------
!                            Read Namelist
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#ifdef MPI
      CONTINUE
#else
      REWIND(5)
#endif
      READ(5,INPDYN,ERR=1,END=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!               Write NAMELIST parameters on the Output file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,3)
      if(velflag=='T')then
       WRITE(6,6)
      else
       WRITE(6,7)      
      endif
      WRITE(6,8)dt 
      WRITE(6,9)tmax 
      WRITE(6,10)ngcf
      do i=1,nat
       do j=1,3
        if(Vxyz(j,i)/=0.0)WRITE(6,11)j,i,Vxyz(j,i)
       enddo
      enddo
      if(resflag=='T')then
       WRITE(6,4)
      else
       WRITE(6,5)      
      endif
      if(snapshot=='T')then
       WRITE(6,12)
      else
       WRITE(6,13)      
      endif
!      
      RETURN
!-----------------------------------------------------------------------
!     Namelist Stop errors
!-----------------------------------------------------------------------
    1 WRITE(6,2)
      STOP
!-----------------------------------------------------------------------
!     Format definitions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    2 FORMAT(/2X,'**********************************************'       &
             /2X,'*                                            *'       &
             /2X,'* SORRY, ERROR IN INPDYN NAMELIST PARAMETERS *'       &
             /2X,'*                                            *'       &
             /2X,'**********************************************')
    3 FORMAT(' Input DYN Options',/,                                    &
             ' -----------------')                                     
    4 FORMAT( ' Restart MD calculation:               (resflag)',9X,'T')    
    5 FORMAT( ' Restart MD calculation:               (resflag)',9X,'F')    
    6 FORMAT( ' Flag constant velocity dynamics:      (velflag)',9X,'T')
    7 FORMAT( ' Flag constant velocity dynamics:      (velflag)',9X,'F') 
    8 FORMAT( ' Time step in fs:                           (dt)',5X,F5.2) 
    9 FORMAT( ' Propagation time in fs:                  (tmax)',F10.2)    
   10 FORMAT( ' Number of GCF files in calculation:      (ngcf)',5X,I5)    
   11 FORMAT( ' Initial velocities per atom:',8X,'Vxyz(',I2,',',I2,')', &
                                                             4X,F6.3)
   12 FORMAT( ' Save MLD file in snapshot-t.mld:     (snapshot)',9X,'T')
   13 FORMAT( ' Save MLD file in snapshot-t.mld:     (snapshot)',9X,'F')
!-----------------------------------------------------------------------
      END

!======================================================================!
!                                                                      !
!   RUNNOFHEADER: Write NOF header on the output file                  !
!   SETORBSPACE: Define the orbital space (NDOC,NSOC,NCWO,NVIR,...)    !
!   POINTERS: Define Pointers of the USER array for the CG subroutine  !
!   OPENFILES: Open all general working files.                         ! 
!   OUTPUTBASIC: Write the basic info on the output file.              !
!   SETNO1: Determine NO1 according to true nuclear charges if NO1=-1  !
!                                                                      !
!======================================================================!

! RUNNOFHEADER
      SUBROUTINE RUNNOFHEADER(NATOMSn,ICHn,MULn,NBFn,NBFauxn,NQMTn,NEn, &
           NAn,NBn,EXn,EYn,EZn,NSHELLn,NSHELLauxn,   &
           NPRIMIn,IAN,IEMOMn,IECPn,IRUNTYP,Cxyz,ZNUC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EFIELDL,RESTART,ERIACTIVATED
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INP_EFIELDL/EX,EY,EZ,EFIELDL
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPFILE_NO1PT2/NO1PT2,NEX
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      COMMON/ELPROP/IEMOM      
      !JFHLewYee: Changed NATOMS allowed dimension from 100 to 1000
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
!       EVEC: Electric Field components (EX,EY,EZ)
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
      EX     = EXn
      EY     = EYn
      EZ     = EZn
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
      LOGICAL APSG,CHKORTHO,ORTHO,HFID,HighSpin,SCALING
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPSCALING/SCALING,NZEROS,NZEROSm,NZEROSr,ITZITER
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_DIIS/DIIS,PERDIIS,NDIIS,NTHDIIS,THDIIS
      COMMON/INPNOF_DAMPEXTRAP/DAMPING,EXTRAP      
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_LAGRANGE/PRINTLAG,DIAGLAG
      COMMON/INPNOF_APSG/APSG,NTHAPSG,THAPSG
      COMMON/INPNOF_ORTHOGONALITY/CHKORTHO,ORTHO
      COMMON/INPNOF_HFID/HFID,NTHRESHEID,THRESHEID,MAXITID,KOOPMANS       
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
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/INPNOF_EXSTA/NESt,OMEGA1      
!-----------------------------------------------------------------------
!     Set NZEROSm equal to NTHRESHL if NZEROSm < NTHRESHL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(SCALING)THEN
        IF(NZEROSm<NTHRESHL)THEN
          NZEROSm = NTHRESHL
          WRITE(6,'(/,6X,38A,/)')'NZEROSm has been set equal to NTHRESHL'
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     Write NAMELIST parameters on the Output file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,3)
      IF(ICOEF==0)THEN
        WRITE(6,4)ICOEF
      ELSEIF(ICOEF==1)THEN
        WRITE(6,50)ICOEF
      ELSEIF(ICOEF==2)THEN
        WRITE(6,51)ICOEF
      ELSEIF(ICOEF==21)THEN
        WRITE(6,52)ICOEF       
      ELSEIF(ICOEF==3)THEN
        WRITE(6,53)ICOEF
      ENDIF
      
      IF(ISOFTMAX==0)THEN
        WRITE(6,54)ISOFTMAX
      ELSE IF(ISOFTMAX==1)THEN
        WRITE(6,55)ISOFTMAX      
      ENDIF
      
      IF(ICOEF>0)THEN
        IF(IORBOPT==1)THEN
          WRITE(6,56)IORBOPT
        ELSE IF(IORBOPT==2)THEN
          WRITE(6,57)IORBOPT       
        ELSE IF(IORBOPT==3)THEN
          WRITE(6,58)IORBOPT       
        ELSE IF(IORBOPT==4)THEN
          WRITE(6,57)IORBOPT       
        ELSE IF(IORBOPT==5)THEN
          WRITE(6,57)IORBOPT       
        ELSE IF(IORBOPT==6)THEN
          WRITE(6,57)IORBOPT       
        ELSE IF(IORBOPT==7)THEN
          WRITE(6,57)IORBOPT       
        ELSE IF(IORBOPT==8)THEN
          WRITE(6,57)IORBOPT       
        ELSE IF(IORBOPT==9)THEN
          WRITE(6,57)IORBOPT       
        ELSE IF(IORBOPT==10)THEN
          WRITE(6,57)IORBOPT       
        ELSE IF(IORBOPT==11)THEN
          WRITE(6,57)IORBOPT       
        END IF
        IF(ICOEF==21)WRITE(6,59)MAXIT21
        WRITE(6,6)MAXIT       
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write the Functional used
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(3<=IPNOF.and.IPNOF<=8)THEN
        WRITE(6,7)IPNOF
        if(IPNOF==7)THEN
          if(Ista==1)WRITE(6,*)                                         &
          'Static Version of Functional PNOF7s:    (Ista)          1'
        end if
      ELSE
        WRITE(6,*)'Stop Program: Select IPNOF between 3 and 8'
        CALL ABRT       
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - -
      IF(MUL>1)THEN
        IF(.NOT.(IPNOF==5.or.IPNOF==7.or.IPNOF==8))THEN
        WRITE(6,'(/A42)')' Stop: IPNOF must be equal 5,7,8 for MULT>1'
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
!- - - - - - - - - - - - - - - - - - - - - - 
      WRITE(6,*)
      IF(RESTART)THEN
        WRITE(6,*)                                                      &
        'Restart calculation Option:             (RESTART)       T'        
      ELSE                                                               
        WRITE(6,*)                                                      &
        'Restart calculation Option:             (RESTART)       F'
      ENDIF
      IF(.NOT.RESTART)WRITE(6,8)INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(HFID)THEN
        WRITE(6,85)NTHRESHEID,MAXITID,KOOPMANS
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF( IORBOPT==1 .and. (ICOEF==1.or.ICOEF==2.or.ICOEF==21) )THEN
        WRITE(6,9)NTHRESHL,NTHRESHE,NTHRESHEC,NTHRESHEN,MAXLOOP
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
      IF(  IORBOPT==1 .and. (.not.SCALING) )THEN
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
            WRITE(6,*)                                                  &
            'Periodic DIIS every NDIIS:              (PERDIIS)       T'        
          ELSE                                                            
            WRITE(6,*)                                                  &
            'DIIS is always applied after NDIIS      (PERDIIS)       F'        
          ENDIF                                                            
        ELSE                                                              
          WRITE(6,*)                                                    &
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
    4 FORMAT(                                                           &
      /1X,'Type of calculation = No Coeff. Opt.:   (ICOEF)     ',I5)    
   50 FORMAT(                                                           &
      /1X,'Type of calculation =    Coeff. Opt.:   (ICOEF)     ',I5)    
   51 FORMAT(                                                           &
      /1X,'Type of calculation = No Occup. Opt.:   (ICOEF)     ',I5)    
   52 FORMAT(                                                           &
      /1X,'Type of calculation = ICOEF2 + ICOEF1:  (ICOEF)     ',I5) 
   53 FORMAT(                                                           &
      /1X,'Type of calculation = Fragment Calc.:   (ICOEF)     ',I5) 
   54 FORMAT(                                                           &
      1X,'Trigonometric functions for Occ. Opt.:  (ISOFTMAX)  ',I5)
   55 FORMAT(                                                           &
      1X,'Use Softmax function for Occ. Opt.:     (ISOFTMAX)  ',I5)     
   56 FORMAT(                                                           &
      1X,'Type of Orb. Opt. = Iterative Diag.:    (IORBOPT)   ',I5)    
   57 FORMAT(                                                           &
      1X,'Type of Orb. Opt. = Orbital Rotations:  (IORBOPT)   ',I5)    
   58 FORMAT(                                                           &
      1X,'Type of Orb. Opt. = Seq. Quad. Prog.:   (IORBOPT)   ',I5)      
   59 FORMAT(                                                           &
      1X,'Maximum Number of Coef Outer Iter.:     (MAXIT21)   ',I5)    
    6 FORMAT(                                                           &
      1X,'Maximum Number of Occ-Coef Outer Iter.: (MAXIT)     ',I5)    
    7 FORMAT(                                                           &
      /1X,'Natural Orbital Functional Selected:    (IPNOF)     ',I5)    
    8 FORMAT(                                                           &
      1X,'Restart Gamma Matrix from GCF file:     (INPUTGAMMA)',I5,     &
      /1X,'Coefficient Matrix from GCF file:       (INPUTC)    ',I5,    &
      /1X,'Diagonal Elements FMIUG from GCF file:  (INPUTFMIUG)',I5,    &
      /1X,'Cartesian Coordinates from GCF file:    (INPUTCXYZ) ',I5)    
   85 FORMAT(                                                           &
      /1X,'Hartree-Fock Calc. using ID Method:     (HFID)', 10X,'T',    &
      /1X,'Threshold Energy Convergence=10**(-N):  (NTHRESHEID)',I5,    &
      /1X,'Max. Number of External Iterations:     (MAXITID)   ',I5,    &
      /1X,'Ion. Potentials by Koopmans Theorem:    (KOOPMANS)  ',I5)    
    9 FORMAT(                                                           &
      /1X,'Threshold Lambda Convergence=10**(-N):  (NTHRESHL)  ',I5,    &
      /1X,'Threshold Energy Convergence=10**(-N):  (NTHRESHE)  ',I5,    &
      /1X,'Threshold Energy Convergence=10**(-N):  (NTHRESHEC) ',I5,    &
      /1X,'Threshold Energy Convergence=10**(-N):  (NTHRESHEN) ',I5,    &
      //1X,'Max. Number of Inner Coef. Iterations:  (MAXLOOP)   ',I5)    
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
      1X,'Threshold to begin DIIS = 10**(-N):     (NTHDIIS)   ',I5,     &
      /1X,'Number of considered Loops in DIIS:     (NDIIS)     ',I5)    
   14 FORMAT(                                                           &
      1X,'Threshold APSG Exp. Coef. = 10**(-N):   (NTHAPSG)   ',I5)    
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
      LOGICAL EFIELDL,HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INP_EFIELDL/EX,EY,EZ,EFIELDL
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
!     EVEC: Electric Field components
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(EX.ne.0.0d0.or.EY.ne.0.0d0.or.EZ.ne.0.0d0)THEN
        EFIELDL=.TRUE.
      ELSE
        EFIELDL=.FALSE.
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
