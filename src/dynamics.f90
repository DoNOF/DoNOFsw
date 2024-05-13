!======================================================================!
!                                                                      !
!                    NOF-BO-MD S U B R O U T I N E S                   !
!                                                                      !
!             2021  Module implemented by Alejandro Rivero             !
!                                                                      !
!                ( J. Chem. Phys. 160, 071102, 2024 )                  !
!                                                                      !
!======================================================================!
!                                                                      !
!   MOLDYN         : Molecular Dynamics routine                        !
!   BEEVER         : Propagator routine                                !
!   PES            : Potential Energy Surface                          !
!   XYZWRITER      : Write to xyz trajectory files                     !
!   MAXNUCDIST     : Maxima Nuclear Distance                           !
!                                                                      !
!======================================================================!

! MOLDYN
SUBROUTINE MOLDYN(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,     &
                  ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,KATOM,KTYPE,KLOC, &
                  INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,CG, &
                  CH,CI,GRADS,IRUNTYP,DIPS)

 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!-----------------------------------------------------------------------
! interface  
!-----------------------------------------------------------------------
 LOGICAL EnergyFileExists
 INTEGER :: NINTEG,IDONTW,IEMOM,NAT,IRUNTYP,NBF,NBFaux,NSHELL,NPRIMI
 INTEGER,DIMENSION(NAT) :: IAN,IMIN,IMAX
 INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
 INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
 INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
 DOUBLE PRECISION,DIMENSION(NAT) :: ZAN,ZMASS 
 DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
 DOUBLE PRECISION,DIMENSION(NPRIMI):: C1,C2,EX,CS,CP,CD,CF,CG,CH,CI
 DOUBLE PRECISION,DIMENSION(3):: DIPS
 DOUBLE PRECISION,DIMENSION(3*NAT):: GRADS
!
 CHARACTER(LEN=1) :: dflag
 COMMON/INPDYN_FLAGS/dflag(3,101)
 DOUBLE PRECISION :: Vxyz
 COMMON/INPDYN_VELOCITY/Vxyz(3,101) 
 COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
 COMMON/EHFEN/EHF,EN
 COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
 COMMON/ECP2/CLP(404),ZLP(404),NLP(404),KFRST(101,6),                   &
             KLAST(101,6),LMAX(101),LPSKIP(101),IZCORE(101)
!-----------------------------------------------------------------------
! internal variables  
!-----------------------------------------------------------------------
 INTEGER :: iatom, init, ngcf
 INTEGER, parameter :: ifstopdyn = 24
 DOUBLE PRECISION,DIMENSION(NAT)  :: mass
 DOUBLE PRECISION,DIMENSION(3,NAT) :: r0,v0,a0,r1,v1,a1
 DOUBLE PRECISION :: t, tres, dt, tmax
 DOUBLE PRECISION :: Epot, Ekin, Etot
 CHARACTER(LEN=1) :: resflag,velflag
 CHARACTER(LEN=50) :: line
 CHARACTER(LEN=2) :: atom
!-----------------------------------------------------------------------
! Constants for transformation units
!-----------------------------------------------------------------------
 real(kind=8),parameter :: amtoau = 1.d0/0.52917720859d0 
 real(kind=8),parameter :: fstoau = 41.341373336561364d0
 real(kind=8),parameter :: amutoau = 1822.888485540950d0
 real(kind=8),parameter :: evtoau = 1.0d0/27.211608d0
 real(kind=8),parameter :: evau = 27.2113845  !! 1au = evau eV 
!-----------------------------------------------------------------------
! Welcome message
!-----------------------------------------------------------------------
 write(6,'(/1X,32A)')' Molecular Dynamics Calculation ' 
 write(6,'(1X,32(1H*),/)')
!-----------------------------------------------------------------------
! read dynamics setup
!-----------------------------------------------------------------------
 call NAMELIST_DYNINP(nat,resflag,velflag,dt,tmax,ngcf)
!-----------------------------------------------------------------------
! restart molecular dynamics calculation (file 31)
!-----------------------------------------------------------------------
 if (resflag == 'T') then
  rewind(31)          
  open(31,FILE='DYNl.xyz',STATUS='OLD',                                 &
          FORM='FORMATTED',ACCESS='SEQUENTIAL')
  read(31,*) line
  read(31,*) tres
  do iatom=1,nat
   read(31,*) atom,Cxyz(1,iatom), Cxyz(2,iatom), Cxyz(3,iatom),         &
                   Vxyz(1,iatom), Vxyz(2,iatom), Vxyz(3,iatom)                 
  enddo
  close(31)
! restart message
  write(6,*)'Note: t0, Cxyz & Vxyz are redefined according to DYNl.xyz'
  Cxyz = Cxyz*amtoau
 endif 
!-----------------------------------------------------------------------
 Vxyz=Vxyz*(amtoau/fstoau) ! Velocity in au/au
 dt=dt*fstoau ! dt in au
 tmax=tmax*fstoau ! tmax in au
!-----------------------------------------------------------------------
! init variables for the integrator
!-----------------------------------------------------------------------
 init=0
 distmax0 = 10.0d0*amtoau   ! Maximum internuclear distance
 mass=ZMASS*amutoau

 r0=Cxyz
 do i=1,nat
  do j=1,3
   v0(j,i)=Vxyz(j,i)  
  enddo
 enddo 
 a0=0.d0

 r1=0.d0
 v1=0.d0
 a1=0.d0 

 t=0.d0
 epot=0.d0

 if (resflag == 'T') t = tres*fstoau
!-----------------------------------------------------------------------
! init integrator message
!-----------------------------------------------------------------------
 write(6,*)
 write(6,*) '-------------------------------------------------------'
 write(6,*) 'Init integrator, time = ', t/fstoau , ' fs'
 write(6,*) '-------------------------------------------------------'
!-----------------------------------------------------------------------
! init integrator
!-----------------------------------------------------------------------
 ICOEFORI = ICOEF
 ICOEF = 21
 call beever(init,NAT,mass,dt,t,r0,v0,a0,r1,v1,a1,ngcf,dflag,velflag,   &
             NINTEG,IDONTW,IEMOM,NBF,NBFaux,NSHELL,NPRIMI,ZAN,IAN,IMIN, &
             IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP, &
             C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,DIPS)
 ICOEF = ICOEFORI             
!-----------------------------------------------------------------------
! compute potential energy (au) in t + dt
!-----------------------------------------------------------------------
 Epot = EELEC + EN
!-----------------------------------------------------------------------
! compute kinetic energy of the system in t + dt
!-----------------------------------------------------------------------
 Ekin = 0.d0

 do i=1,nat
  Ekin = Ekin + 0.5d0*mass(i)                                           &
       * ( v1(1,i)*v1(1,i)+v1(2,i)*v1(2,i)+v1(3,i)*v1(3,i) )
 enddo
!-----------------------------------------------------------------------
! compute total energy of the system
!-----------------------------------------------------------------------
 Etot = Epot + Ekin
!-----------------------------------------------------------------------
! write information of propagation in file xyz
!-----------------------------------------------------------------------
 call XYZWRITER(init,nat,zan,izcore,t,r1,v1,Epot,Ekin,Etot)
!-----------------------------------------------------------------------
! write information obout energy conservation in dynamics
!-----------------------------------------------------------------------
 inquire(FILE='endyn.dat',EXIST=EnergyFileExists)
 if(EnergyFileExists)then  
  OPEN(32,FILE='endyn.dat',STATUS='OLD',POSITION='APPEND',              &
          ACTION='WRITE',FORM='FORMATTED',ACCESS='SEQUENTIAL')
 else
  OPEN(32,FILE='endyn.dat',STATUS='NEW',ACTION='WRITE',                 &
          FORM='FORMATTED',ACCESS='SEQUENTIAL')
 end if 
 write(32,'(4(2x,f14.8))') t/fstoau,Epot,Ekin,Etot
 close(32)
!-----------------------------------------------------------------------
! begin propagation
!-----------------------------------------------------------------------            
 init = 1

 do while (abs(t) < tmax)
!-----------------------------------------------------------------------
! Propagation time message
!-----------------------------------------------------------------------
  write(6,*)
  write(6,*) '---------------------------------------------------'
  write(6,*) ' Propagation time = ', t/fstoau + dt/fstoau, ' fs'
  write(6,*) '---------------------------------------------------'
!----------------------------------------------------------------------- 
  call beever(init,nat,mass,dt,t,r0,v0,a0,r1,v1,a1,ngcf,dflag,velflag,  &
              NINTEG,IDONTW,IEMOM,NBF,NBFaux,NSHELL,NPRIMI,ZAN,IAN,IMIN,&
              IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,&
              C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,DIPS)
!-----------------------------------------------------------------------
!  compute potential energy (au)
!-----------------------------------------------------------------------
   Epot = EELEC + EN
!-----------------------------------------------------------------------
!  compute kinetic energy of the system
!-----------------------------------------------------------------------
   Ekin = 0.d0
   do i=1,nat
    Ekin = Ekin + 0.5d0*mass(i)*( v1(1,i)*v1(1,i) &
                + v1(2,i)*v1(2,i)+v1(3,i)*v1(3,i) )
   enddo
!-----------------------------------------------------------------------
!  compute total energy of the system
!-----------------------------------------------------------------------
   Etot = Epot + Ekin
!-----------------------------------------------------------------------
!  write information of propagation in file xyz
!-----------------------------------------------------------------------
   trealpart = t/fstoau - int(t/fstoau+1.0d-03)
   if ( trealpart < 0.1*dt/fstoau )  &  ! write trayectory each 1 fs
   call XYZWRITER(init,nat,zan,izcore,t,r1,v1,Epot,Ekin,Etot)      
!-----------------------------------------------------------------------
!  write information about energy conservation in dynamics
!-----------------------------------------------------------------------
   OPEN(32,FILE='endyn.dat',STATUS='OLD',POSITION='APPEND',             &
           ACTION='WRITE',FORM='FORMATTED',ACCESS='SEQUENTIAL')
   write(32,'(4(2x,f14.8))') t/fstoau,Epot,Ekin,Etot
   close(32)
!-----------------------------------------------------------------------
!  stop dynamics if Maximum internuclear distance > 10 Ang
!-----------------------------------------------------------------------
   if(t>0.1*tmax)then
    call MAXNUCDIST(nat,r1,distmax)
    if(distmax>distmax0)then
     write(6,*)'Stop: Maximum internuclear distance > 10 Ang'
     return
    end if 
   end if
!-----------------------------------------------------------------------  
 enddo
!-----------------------------------------------------------------------            
end

! MAXNUCDIST
SUBROUTINE MAXNUCDIST(nat,Cxyz,distmax)
 IMPLICIT NONE
 INTEGER,INTENT(IN) :: NAT
 DOUBLE PRECISION,INTENT(OUT) :: distmax
 DOUBLE PRECISION,DIMENSION(3,NAT),INTENT(IN) :: Cxyz
 INTEGER :: I,J
 DOUBLE PRECISION :: RR 
!-----------------------------------------------------------------------
 distmax = 0.0d0
 do I = 1,NAT
  do J = I+1,NAT
   RR = (Cxyz(1,I)-Cxyz(1,J))**2 + (Cxyz(2,I)-Cxyz(2,J))**2 +           &                       
        (Cxyz(3,I)-Cxyz(3,J))**2
   RR = DSQRT(RR)
   if(RR>distmax) distmax = RR
  end do
 end do
!-----------------------------------------------------------------------
 RETURN
END

! BEEVER
subroutine beever(init,nat,mass,dt,t,r0,v0,a0,r1,v1,a1,ngcf,dflag,      &
                  velflag,NINTEG,IDONTW,IEMOM,NBF,NBFaux,NSHELL,NPRIMI, &
                  ZAN,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,  &
                  KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,     &
                  GRADS,IRUNTYP,DIPS)  
!-----------------------------------------------------------------------
! init: flag to init propagator (in)
! nat: number of atoms in the system (in)
! mass: mass of the atoms of the system (in)
! dt : time step
! t : time of the trajectory         
! r0: position at time t-dt (in, not used), and at time t (out)
! v0: velocity at time t-dt (in, not used), and at time t (out)
! a0: acceleration at time t-dt (in), and at time t (out)
! r1: position at time t (in), and at time t+dt (out)
! v1: velocity at time t (in), and at time t+dt (out)
! a1: acceleration at time t (in), and at time t+dt (out)
! epot: potential energy at time t+dt (out)
!-----------------------------------------------------------------------
 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!-----------------------------------------------------------------------
! interface  
!-----------------------------------------------------------------------
 INTEGER, intent(in) :: init, nat, ngcf
 INTEGER :: NINTEG,IDONTW,IEMOM,IRUNTYP
 INTEGER,DIMENSION(NAT) :: IAN,IMIN,IMAX
 INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
 INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
 INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
 CHARACTER(LEN=1), dimension (3,nat), intent(in) :: dflag
 CHARACTER(LEN=1), intent(in) :: velflag
 DOUBLE PRECISION, dimension (nat), intent(in) :: mass
 DOUBLE PRECISION, dimension (3,nat), intent(inout) :: r0,v0,a0,r1,v1,a1
 DOUBLE PRECISION, intent(in)  :: dt
 DOUBLE PRECISION, intent(inout)  :: t
 DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
 DOUBLE PRECISION,DIMENSION(NPRIMI)::C1,C2,EX,CS,CP,CD,CF,CG,CH,CI
 DOUBLE PRECISION,DIMENSION(3,NAT):: GRADS
 DOUBLE PRECISION,DIMENSION(3) :: DIPS
!-----------------------------------------------------------------------
! internal variables  
!-----------------------------------------------------------------------
 INTEGER :: i
 DOUBLE PRECISION, dimension (3,NAT) :: r2,v2,a2,forces
 real(kind=8), parameter :: fstoau = 41.341373336561364d0   
!-----------------------------------------------------------------------
! init or propagation
!-----------------------------------------------------------------------
 select case (init)
!-----------------------------------------------------------------------
! init integrator
!-----------------------------------------------------------------------  
  case (0)
!-----------------------------------------------------------------------
! compute forces at initial positions r0 (t)
!-----------------------------------------------------------------------
   CALL PES(t,ngcf,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,ZAN,&
           r0,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,&
           ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,DIPS)
   forces=-GRADS

   if (dt == 0) then
    goto 10  ! special case to not allow double calculation
   endif
!-----------------------------------------------------------------------
! compute initial acceleration a0 at initial positions r0 (t)
!-----------------------------------------------------------------------
   do i=1,nat
    a0(:,i)=forces(:,i)/mass(i)
   enddo
!-----------------------------------------------------------------------
! position r1 at t + dt (Verlet algorithm)
!-----------------------------------------------------------------------
   do i=1,nat
    r1(:,i)=r0(:,i) + dt*v0(:,i) + 0.5d0*dt*dt*a0(:,i)
   enddo
!-----------------------------------------------------------------------
! compute forces at r1 (position at t + dt) 
!-----------------------------------------------------------------------
  write(6,*)
  write(6,*) '---------------------------------------------------'
  write(6,*) ' Propagation time = ', t/fstoau + dt/fstoau, ' fs'
  write(6,*) '---------------------------------------------------'

   CALL PES(t+dt,ngcf,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI, &
            ZAN,r1,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,&
        KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,DIPS)
   
   forces=-GRADS
!-----------------------------------------------------------------------
! compute acceleration a1 at t + dt
!-----------------------------------------------------------------------
   do i=1,nat
    a1(:,i)=forces(:,i)/mass(i)
   enddo
!-----------------------------------------------------------------------
! compute velocity v1 at t + dt (Verlet algorithm)
!-----------------------------------------------------------------------
   do i=1,nat
    v1(:,i)=v0(:,i) + 0.5d0*dt*(a0(:,i) + a1(:,i))
   enddo
!-----------------------------------------------------------------------
! check dynamics flags for selective dynamics
!-----------------------------------------------------------------------
   do i=1,nat
    
    if (dflag(1,i) == 'F') then 
     a0(1,i) = 0.d0
     v0(1,i) = 0.d0
     a1(1,i) = 0.d0
     v1(1,i) = 0.d0
     r1(1,i) = r0(1,i)
    endif 
    
    if (dflag(2,i) == 'F') then
     a0(2,i) = 0.d0
     v0(2,i) = 0.d0
     a1(2,i) = 0.d0
     v1(2,i) = 0.d0
     r1(2,i) = r0(2,i)
    endif 
    
    if (dflag(3,i) == 'F') then
     a0(3,i) = 0.d0
     v0(3,i) = 0.d0
     a1(3,i) = 0.d0
     v1(3,i) = 0.d0
     r1(3,i) = r0(3,i)
    endif

    if (velflag == 'T') then
     v1=v0
     a0=0.d0
     a1=0.d0
    endif 
   
   enddo
!-----------------------------------------------------------------------
! increase time in dt
!-----------------------------------------------------------------------
   t = t + dt

  case default
!-----------------------------------------------------------------------
! Equation of motion integration .....
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  new positions r2 at t + dt (Beeman algorithm)
!----------------------------------------------------------------------- 
   do i=1,nat
    r2(:,i)=r1(:,i) + dt*v1(:,i) + dt*dt*(4.0d0*a1(:,i) - a0(:,i))/6.0d0
   enddo
!-----------------------------------------------------------------------
! compute forces at r2 (position at t + dt) 
!-----------------------------------------------------------------------
   CALL PES(t+dt,ngcf,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,ZAN,&
           r2,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,&
           ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,DIPS)
   
   forces=-GRADS
!-----------------------------------------------------------------------
! compute acceleration a2 at t + dt
!-----------------------------------------------------------------------
   do i=1,nat
    a2(:,i)=forces(:,i)/mass(i)
   enddo
!-----------------------------------------------------------------------
! compute velocity v2 at t + dt (Beeman algorithm)
!-----------------------------------------------------------------------
  do i=1,nat
   v2(:,i) = v1(:,i) + dt*(2.0d0*a2(:,i)+5.0d0*a1(:,i)-a0(:,i))/6.0d0
  enddo
!-----------------------------------------------------------------------
! compute velocity v2 at t + dt (second order Adams–Moulton method)
!-----------------------------------------------------------------------
! do i=1,nat
!  v2(:,i) = v1(:,i) + dt*(5.0d0*a2(:,i)+8.0d0*a1(:,i)-a0(:,i))/12.0d0
! enddo
!-----------------------------------------------------------------------
! store values of r,v and a for the next step
!-----------------------------------------------------------------------
   r0=r1
   v0=v1
   a0=a1

   r1=r2 
   v1=v2
   a1=a2
!-----------------------------------------------------------------------
! check dynamics flags for selective dynamics
!-----------------------------------------------------------------------
   do i=1,nat
    
    if (dflag(1,i) == 'F') then 
     a0(1,i) = 0.d0
     v0(1,i) = 0.d0
     a1(1,i) = 0.d0
     v1(1,i) = 0.d0
     r1(1,i) = r0(1,i)
    endif 
    
    if (dflag(2,i) == 'F') then
     a0(2,i) = 0.d0
     v0(2,i) = 0.d0
     a1(2,i) = 0.d0
     v1(2,i) = 0.d0
     r1(2,i) = r0(2,i)
    endif 
    
    if (dflag(3,i) == 'F') then
     a0(3,i) = 0.d0
     v0(3,i) = 0.d0
     a1(3,i) = 0.d0
     v1(3,i) = 0.d0
     r1(3,i) = r0(3,i)
    endif

    if (velflag == 'T') then
     v1=v0
     a0=0.d0
     a1=0.d0
    endif 
   
   enddo  
!-----------------------------------------------------------------------
!   increase time in dt 
!-----------------------------------------------------------------------   
    t = t + dt

10 end select

end

! PES
subroutine PES(t,ngcf,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI, &
               ZAN,r0,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,  &
               KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,        &
               GRADS,IRUNTYP,DIPS)
 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!-----------------------------------------------------------------------
! interface
!-----------------------------------------------------------------------
 INTEGER ::ngcf,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,IRUNTYP
 INTEGER,DIMENSION(NAT) :: IAN,IMIN,IMAX
 INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
 INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
 INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
 DOUBLE PRECISION, intent(in)  :: t 
 DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
 DOUBLE PRECISION,DIMENSION(NPRIMI)::C1,C2,EX,CS,CP,CD,CF,CG,CH,CI
 DOUBLE PRECISION,DIMENSION(3,NAT):: r0,GRADS 
 DOUBLE PRECISION,DIMENSION(3) :: DIPS
!
 CHARACTER(LEN=1) :: snapshot
 COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
 COMMON/EHFEN/EHF,EN
 COMMON/WRTGCF/IWRTGCF
 COMMON/INPDYN_NSHOT/snapshot
 COMMON/INPNOF_MOLDEN/MOLDEN
!-----------------------------------------------------------------------
! internal variables
!-----------------------------------------------------------------------
 character (4), parameter :: seedname='GCF-'
 character (2)  :: str
 character (len=10) :: fileIn
 character (len=20) :: MLD_FILE
 character (len=8) :: str1

 INTEGER :: i,igcfmin(1)
 DOUBLE PRECISION,DIMENSION(ngcf,3,NAT):: GRADS_temp
 DOUBLE PRECISION, dimension (ngcf) :: Epot, EELEC_temp, EN_temp
 real(kind=8), parameter :: fstoau = 41.341373336561364d0  
!-----------------------------------------------------------------------
! init var
!-----------------------------------------------------------------------
 Epot(:) = 0.d0
!
 if(snapshot=='T')then
  MOLDEN = 1
  write(str1,'(F8.2)')t/fstoau
 endif 
!-----------------------------------------------------------------------
! begin loop over GCF guess
!-----------------------------------------------------------------------
 do i=1,ngcf
!-----------------------------------------------------------------------    
! Save MLD file in snapshot-t.mld
!-----------------------------------------------------------------------
  IF(snapshot=='T')THEN
   if(ngcf==1)then
    MLD_FILE = 'snapshot-'//trim(adjustl(str1))//".mld"
   else
    MLD_FILE = 'snapshot-'//trim(adjustl(str1))//'-'//char(i)//".mld"
   end if
   OPEN(17,FILE=MLD_FILE,STATUS='UNKNOWN',FORM='FORMATTED',             &
           ACCESS='SEQUENTIAL') 
  END IF
!-----------------------------------------------------------------------
! Previous GCF file is not used in BO-MD for GCF-i /= GCF-1
!-----------------------------------------------------------------------
!useGCF  if(i==1)then
!useGCF   IWRTGCF=1
!useGCF  else
!useGCF   IWRTGCF=0
!useGCF  end if
!-----------------------------------------------------------------------
! assign name to the variable fileIn
!-----------------------------------------------------------------------
  write (str,'(i2)') i
  fileIn = trim(seedname)//trim(adjustl(str))
!-----------------------------------------------------------------------
! open GCF-i in unit 3
!-----------------------------------------------------------------------
  OPEN(3,FILE=fileIn,FORM='FORMATTED',ACCESS='SEQUENTIAL')
  write(6,*)
  write(6,*)'Using GCF-',i
!-----------------------------------------------------------------------
! compute Energy and Gradients
!-----------------------------------------------------------------------
  CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,ZAN,r0,&
                IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,   &
                KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,      &
                IRUNTYP,DIPS,0,1)    ! ,0) short print on output file
!-----------------------------------------------------------------------
! store temporal variables
!-----------------------------------------------------------------------
  GRADS_temp(i,:,:) = GRADS
  EELEC_temp(i) = EELEC
  EN_temp(i) = EN
  Epot(i) = EELEC + EN
  CLOSE(3)
  IF(snapshot=='T')CLOSE(17)
 enddo
!-----------------------------------------------------------------------
! index of minimun energy solution
!-----------------------------------------------------------------------
 igcfmin = minloc(Epot)
 write(6,*)
 write(6,*)'Epot min with GCF-',igcfmin(1)
!-----------------------------------------------------------------------
! energy information multiples GCF
!-----------------------------------------------------------------------
! write(34,'(16f14.8)') t/fstoau,(Epot(j), j=1,ngcf)
!-----------------------------------------------------------------------
! passing minumun values to global variables
!-----------------------------------------------------------------------
 
 GRADS=GRADS_temp(igcfmin(1),:,:) 
 EELEC=EELEC_temp(igcfmin(1))
 EN=EN_temp(igcfmin(1))
 
end

! XYZWRITER
subroutine XYZWRITER(init,nat,zan,izcore,tau,rau,vau,Epot,Ekin,Etot)
 implicit none
!-----------------------------------------------------------------------
 logical DynFileExists,DynlFileExists
 integer,intent(in) :: init,nat
 integer,dimension (nat),intent(in) :: izcore
 real(kind=8),dimension (nat),intent(in) :: zan
 real(kind=8),intent(in) :: tau
 real(kind=8),dimension (3,nat),intent(in) :: rau,vau
 real(kind=8),intent(in) :: Epot,Ekin,Etot
!-----------------------------------------------------------------------
 character(len=4),dimension(106) :: ATMLAB
 DATA ATMLAB/'H   ','HE  ','LI  ','BE  ','B   ','C   ',                 &
             'N   ','O   ','F   ','NE  ','NA  ','MG  ',                 &
             'AL  ','SI  ','P   ','S   ','CL  ','AR  ',                 &
             'K   ','CA  ','SC  ','TI  ','V   ','CR  ',                 &
             'MN  ','FE  ','CO  ','NI  ','CU  ','ZN  ',                 &
             'GA  ','GE  ','AS  ','SE  ','BR  ','KR  ',                 &
             'RB  ','SR  ','Y   ','ZR  ','NB  ','MO  ',                 &
             'TC  ','RU  ','RH  ','PD  ','AG  ','CD  ',                 &
             'IN  ','SN  ','SB  ','TE  ','I   ','XE  ',                 &
             'CS  ','BA  ','LA  ','CE  ','PR  ','ND  ',                 &
             'PM  ','SM  ','EU  ','GD  ','TB  ','DY  ',                 &
             'HO  ','ER  ','TM  ','YB  ','LU  ','HF  ',                 &
             'TA  ','W   ','RE  ','OS  ','IR  ','PT  ',                 &
             'AU  ','HG  ','TL  ','PB  ','BI  ','PO  ',                 &
             'AT  ','RN  ','FR  ','RA  ','AC  ','TH  ',                 &
             'PA  ','U   ','NP  ','PU  ','AM  ','CM  ',                 &
             'BK  ','CF  ','ES  ','FM  ','MD  ','NO  ',                 &
             'LR  ','RF  ','X   ','BQ  '/
!-----------------------------------------------------------------------
 integer :: i,iznuc
 real(kind=8) :: tfs
 real(kind=8),dimension (3,nat) :: ram,vam
!-----------------------------------------------------------------------
! Constants for transformation units
!-----------------------------------------------------------------------
 real(kind=8),parameter :: amtoau=1.d0/0.52917720859d0 
 real(kind=8),parameter :: fstoau=41.341373336561364d0
 real(kind=8),parameter :: amutoau=1822.888485540950d0
 real(kind=8),parameter :: evtoau=1.0d0/27.211608d0 
!-----------------------------------------------------------------------
 inquire(FILE='DYN.xyz',EXIST=DynFileExists)
 if(DynFileExists)then 
  OPEN(30,FILE='DYN.xyz',STATUS='OLD',POSITION='APPEND',ACTION='WRITE', &
          FORM='FORMATTED',ACCESS='SEQUENTIAL')
 else
  OPEN(30,FILE='DYN.xyz',STATUS='NEW',ACTION='WRITE',                   &
          FORM='FORMATTED',ACCESS='SEQUENTIAL')
 end if
!
 inquire(FILE='DYNl.xyz',EXIST=DynlFileExists)
 if (DynlFileExists)OPEN(31,FILE='DYNl.xyz',STATUS='OLD',               &
                            FORM='FORMATTED',ACCESS='SEQUENTIAL')
!-----------------------------------------------------------------------
 select case(init)

  case(0)
  
!-----------------------------------------------------------------------
!  write to xyz trajectory file
!-----------------------------------------------------------------------
   tfs = tau/fstoau    
   write(30,*) nat
   write(30,'(4(a6),4(f14.8))') 't,','Epot,','Ekin,','Etot,',           &
                                tfs,Epot,Ekin,Etot
   do i=1,nat
    ram(:,i)=rau(:,i)/amtoau
    vam(:,i)=vau(:,i)*(fstoau/amtoau)
    iznuc = int(zan(i)) + izcore(i)
    write(30,'(a4, 6(1x,f14.8))')ATMLAB(iznuc),ram(:,i),vam(:,i)
   enddo

  case default

    tfs=tau/fstoau 
    write(30,*) nat
    write(30,'(4(a6),4(f14.8))') 't,','Epot,','Ekin,','Etot,',          &
                                 tfs,Epot,Ekin,Etot
    do i=1,nat
     ram(:,i)=rau(:,i)/amtoau
     vam(:,i)=vau(:,i)*(fstoau/amtoau) 
     iznuc = int(zan(i)) + izcore(i)
     write(30,'(a4, 6(1x,f14.8))')ATMLAB(iznuc),ram(:,i),vam(:,i)
    enddo
!-----------------------------------------------------------------------
!   write the last xyz backup file (DYNl.xyz)
!-----------------------------------------------------------------------
    rewind 31
    write(31,*) nat
    write(31,*) tfs
    do i=1,nat
     ram(:,i)=rau(:,i)/amtoau
     vam(:,i)=vau(:,i)*(fstoau/amtoau)
     iznuc = int(zan(i)) + izcore(i)     
     write(31,'(a4, 6(1x,f14.8))')ATMLAB(iznuc),ram(:,i),vam(:,i)
    enddo
    
  end select 
  
  close(30) 
  close(31)

end
