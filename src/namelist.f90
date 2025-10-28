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
      SUBROUTINE NAMELIST_INPRUN(ITYPRUN,ICHARG,MULT,NINTEG,IDONTW,     &
                                 IEMOMENTS,EX,EY,EZ,LIBCINT,IECPO,      &
                                 IHSSCAL,IPROJECT,ISIGMA,IGTYP,NATmax,  &
                                 NSHELLmax,NPRIMImax,IHUBBARD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
#include "mpip.h"      
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
      COMMON/INFO/IUNTRD            
      COMMON/INTFIL/NINTMX           
      LOGICAL EFIELDL
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      COMMON/ELPROP/IEMOM      
      CHARACTER(8):: UNITS
      COMMON/CONTROL/UNITS
      COMMON/INTOPT/CUTOFF,ISCHWZ,IECP,NECP
      COMMON/RUNTYPE/IRUNTYP
      COMMON/USELIBCINT/ILIBCINT
      COMMON/USEHUBBARD/IHUB
      LOGICAL DONTW,USELIB,USEHUB,HSSCAL,PROJECT
      DOUBLE PRECISION,DIMENSION(3) :: EVEC
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
      CHARACTER(4) :: GTYP, CART, SPH
      DATA CART,SPH/'CART','SPH '/
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
! USELIB           Use LIBCINT open source library for ERI calculation
!       = F        HONDO Calculator
!       = T        LIBCINT (Default)
!
! GTYP             Type of Gaussian functions
!       = CART     Cartesian (Default)
!       = SPH      Spherical (only if LIBCINT)
!
! USEHUB           Use Hubbard Model
!       = F        (Default)
!
! DONTW            Do not write 2e- integrals on the disk (Unit=1)
!       = T        (Default)
!
! ERITYP           Typ of ERIs used in calculations
!       = FULL     4c ERIs
!       = RI       3c/2c ERIs for RI Approximation (Default)
!       = MIX      3c/2c ERIs for Resolution of the Identity (RI) App.
!                  once converged change to 4c ERIs (FULL)
!
! CUTOFF           The Schwarz screening cut off for NAT>5
!       = 1.0D-9   (Default)
!
! RITYP            Typ of Auxiliary Basis
!       = JKFIT    Read from jkfit files (Default)
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
                      GTYP,USEHUB,DONTW,ERITYP,CUTOFF,RITYP,GEN,SMCD,   &
                      HSSCAL,PROJECT,ISIGMA,NATmax,NSHELLmax,NPRIMImax
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
      USELIB    = .TRUE.
      GTYP      = CART
      USEHUB    = .FALSE.
      DONTW     = .TRUE.
      ERITYP    = RI
      CUTOFF    = 1.0D-09
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
!     REWIND(5)
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
!     Use LIBCINT open source library for ERI calculation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(USELIB)THEN
! - - - - - - - - - - - - - - - -
       ILIBCINT = 1
!      Determine GTYP
       IF(GTYP==CART)THEN
        IGTYP = 1
       ELSE IF(GTYP==SPH)THEN
        IGTYP = 2
       ELSE
        WRITE(6,7)
        CALL ABRT
       END IF
! - - - - - - - - - - - - - - - -
      ELSE  ! HONDO Integrals
       ILIBCINT = 0
       IF(GTYP==CART)THEN
        IGTYP = 1
       ELSE
        WRITE(6,8)
        CALL ABRT
       END IF
      ENDIF
      LIBCINT = ILIBCINT
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
      EFX = EVEC(1)
      EFY = EVEC(2)
      EFZ = EVEC(3)
      IF(EFX.ne.0.0d0.or.EFY.ne.0.0d0.or.EFZ.ne.0.0d0)THEN
       EFIELDL=.TRUE.
      ELSE
       EFIELDL=.FALSE.
      ENDIF
!     NAMELIST_INPRUN
      EX = EVEC(1)
      EY = EVEC(2)
      EZ = EVEC(3)
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
      IF(IECP>0 .and. NATmax>1001)THEN         ! due to COMMON/ECP2/
        WRITE(6,9)
        CALL ABRT
      END IF
      IF(IECP>0 .and. NSHELLmax>600)THEN      ! due to MAPSHL
        WRITE(6,10)
        CALL ABRT
      END IF
        IF(ILIBCINT==1 .and. NSHELLmax>600)THEN ! due to BASLIB
        WRITE(6,11)
        CALL ABRT
      END IF
      IF(IERITYP==2 .and. NSHELLmax>700)THEN  ! due to COMMON/NSHELaux/
        WRITE(6,12)
        CALL ABRT
        END IF
      IF(IERITYP==2 .and. NPRIMImax>2000)THEN ! due to COMMON/EXCaux/
        WRITE(6,13)
        CALL ABRT
      END IF
      IF(IRUNTYP==5 .and. NATmax>1001)THEN ! due to COMMON/INPDYN_FLAGS/
        WRITE(6,14)
        CALL ABRT
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ICH = ICHARG                                                      
      MUL = MULT
      IEMOMENTS = IEMOM
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
    7 FORMAT(/1X,'Stop: GTYP must be CART or SPH',/) 
    8 FORMAT(/1X,'Stop: GTYP must be CART for USELIB=F',/)
    9 FORMAT(/1X,'Stop: NATmax must be <= 1001 with ECP',/) 
   10 FORMAT(/1X,'Stop: NHELLmax must be <= 600 with ECP',/)
   11 FORMAT(/1X,'Stop: NHELLmax must be <= 600 with USELIB = T',/)
   12 FORMAT(/1X,'Stop: NHELLmax must be <= 700 with ERITYP = RI',/)
   13 FORMAT(/1X,'Stop: NPRIMImax must be <= 2000 with ERITYP = RI',/)
   14 FORMAT(/1X,'Stop: NATmax must be <= 1001 with RUNTYP = DYN',/)
!-----------------------------------------------------------------------
      END

! NAMELIST_NOFINP                                                     
      SUBROUTINE NAMELIST_NOFINP(IRUNTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HFDAMP,HFEXTRAP,HFDIIS,DIIS,PERDIIS,DAMPING,EXTRAP
      LOGICAL RESTART,FROZEN,PRINTLAG,DIAGLAG,APSG,CHKORTHO,ORTHO
      LOGICAL ERPA,OIMP2,MBPT,SC2MCPT,HighSpin,SCALING
      LOGICAL AUTOLR
      DOUBLE PRECISION LR,FACT,BETA1,BETA2
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_RHF/CONVRHFDM,IRHF,IRHFTYP,NCONVRHF,MAXITRHF
      COMMON/INPSCALING/SCALING,NZEROS,NZEROSm,NZEROSr,ITZITER
      COMMON/INPADAM/LR,FACT,BETA1,BETA2,AUTOLR
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
      COMMON/INPNOF_HFID/KOOPMANS
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
      COMMON/INPNOF_MOD/Imod
      COMMON/INPNOF_EXSTA/OMEGA1,NESt 
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/INPFILE_NO1PT2/NO1PT2,NEX
      COMMON/USEHUBBARD/IHUB
      LOGICAL SMCD
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      LOGICAL EFIELDL
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      LOGICAL AUTOZEROS
      INTEGER:: IRUNTYP
      INTEGER :: NESt
      DOUBLE PRECISION :: OMEGA1
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
!                      = 3      (Default)
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
!                      = 21     ICOEF=2 & MAXIT21, then ICOEF=1
!                      = 3      Optimize by all ONs and core-fragment 
!                               orbitals. The rest of fragment orbitals 
!                               remain frozen
!
!.......... ISOFTMAX            Use Softmax function for ON (Gamma) opt.
!                      = 0      Trigonometric Parametrization for ON
!                      = 1      Softmax function (Default)
!
!.......... IORBOPT             Select method for NO optimization
!
!                      = 1      Iterative diagonalization (OrbOptFMIUGr)
!                      = 2      Adaptative Momentum (ADAM) (Default)
!                      = 3      ADABelief
!                      = 4      YOGI
!                      = 5      Decaying Momentum (DEMON)
!                      = 6      Sequential Quadratic Program (OrbOptSQP)
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
!.......... IRHF                Restricted Hartree-Fock Calculation
!                      = 0      Not obtaining HF orbitals
!                      = 1      Self Consistent Field (SCF) Method
!                               (only works with EFIELDL=.FALSE.)
!                      = 2      Orbital rotaions through ADAM (Default)
!                      = 3      Iterative Diagonalization (ID) Method
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
!..........             Calculate IPs using Koopmans' Theorem
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
!.......... Imod                Select versions of GNOFx
!                      = 0      GNOF (Default)
!                      = 1      GNOFm
!                      = 2      GNOFs
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
!                      = 4      (Default)
!
!.......... NTHRESHE            Convergence of the total energy
!                               THRESHE=10.0**(-NTHRESHE)
!                      = 8      (Default)
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
!                      = 10     (Default, recommended for ADAM method)
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
                      IRHF,NCONVRHF,MAXITRHF,HFDAMP,HFEXTRAP,HFDIIS,    &
                      KOOPMANS,IPNOF,Ista,Imod,HighSpin,NCWO,           &
                      NTHRESHL,NTHRESHE,NTHRESHEC,NTHRESHEN,            &
                      MAXLOOP,SCALING,AUTOZEROS,NZEROS,NZEROSm,NZEROSr, &
                      ITZITER,DIIS,NTHDIIS,NDIIS,PERDIIS,DAMPING,       &
                      EXTRAP,ERPA,OIMP2,NO1PT2,SC2MCPT,NEX,             &
                      RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,NPRINT,      &
                      INPUTCXYZ,IWRITEC,IMULPOP,MBPT,PRINTLAG,DIAGLAG,  &
                      IEKT,IAIMPAC,IFCHK,MOLDEN,INICOND,NOUTRDM,        &
                      NTHRESHDM,NSQT,NOUTCJK,NTHRESHCJK,NOUTTijab,      &
                      NTHRESHTijab,APSG,NTHAPSG,ORTHO,CHKORTHO,FROZEN,  &
                      IFROZEN,ICGMETHOD,NESt,OMEGA1,LR,FACT,BETA1,BETA2
!-----------------------------------------------------------------------
!     Preset values to namelist variables
!-----------------------------------------------------------------------
      MAXIT=1000
      MAXIT21=3
      
!     Type of Calculation
      ICOEF=1
      ISOFTMAX=1
      IORBOPT=2
      IEINI=0
      NO1=0
      
!     Hartree-Fock
      IRHF=2
      NCONVRHF=5
      MAXITRHF=100
      HFDAMP=.TRUE.
      HFEXTRAP=.TRUE.
      HFDIIS=.TRUE.
      KOOPMANS=0
      
!     PNOF Selection
      IPNOF=8
      Ista=0                                      
      Imod=0                      ! GNOF
      HighSpin=.FALSE.            ! Multiplet
      NCWO=-1
      
!     Convergence Criteria in NOF calculation
      NTHRESHL=4
      NTHRESHE=8
      NTHRESHEC=8
      NTHRESHEN=10
      MAXLOOP=10                  ! ADAM Method (IRHF=2)
!     MAXLOOP=30                  ! ID Method   (IRHF=3)

!     Options for the Orb. Opt. ID Method
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
!     Note: Override AUTOZEROS to False
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IORBOPT>1 .and. AUTOZEROS)AUTOZEROS=.FALSE.
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
       WRITE(6,'(/1X,57A/)')                                            &
        'NTHRESHEC is very small so it has been increased to 8 !!!'
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Convergence Criteria
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CONVRHFDM   = 10.0**(-NCONVRHF) 
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
!     IRHF=1 only works with false
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IRHF==1 .and. (EFIELDL .eqv. .TRUE.))THEN
       WRITE(6,'(/,1X,38A,/)')'IRHF=1 only works with EFIELDL=.FALSE.'
       STOP
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Mandatory Options with RUNTYP=OPTGEO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IRUNTYP==3)THEN
       if(IRHF>0)then
        IRHF = 0
        WRITE(6,'(/,1X,35A,/)')'!OPTGEO: IRHF has been set equal FALSE'
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
       if(IRHF>0)then
        IRHF = 0
        WRITE(6,'(/,1X,39A,/)')'!HESS cal.: IRHF has been set equal 0'
       end if
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Mandatory Options with RUNTYP=DYN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IRUNTYP==5)THEN
       if(IRHF>0)then
        IRHF = 0
        WRITE(6,'(/,1X,33A,/)')'!DYN: IRHF has been set equal 0'
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
       IF(IORBOPT==1)INPUTFMIUG = 1
       WRITE(6,'(  1X,38A,/)')'!INPUTGAMMA=1, INPUTC=1'
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
