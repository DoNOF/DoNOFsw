SUBROUTINE MOLDYN(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,     &
                  ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,KATOM,KTYPE,KLOC, &
                  INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,CG, &
                  CH,CI,GRADS,IRUNTYP,DIPS,IZCORE,SIZE_ENV,ENV,ATM,     &
                  NBAS,BAS,IGTYP,KTYPEaux,NSHELLaux)

 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 LOGICAL EnergyFileExists
 INTEGER :: NINTEG,IDONTW,IEMOM,NAT,IRUNTYP,NBF,NBFaux,NSHELL,NPRIMI
 INTEGER,DIMENSION(NAT) :: IAN,IMIN,IMAX,IZCORE
 INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
 INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
 INTEGER,DIMENSION(NSHELLaux) :: KTYPEaux
 INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
 DOUBLE PRECISION,DIMENSION(NAT) :: ZAN,ZMASS 
 DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
 DOUBLE PRECISION,DIMENSION(NPRIMI):: C1,C2,EX,CS,CP,CD,CF,CG,CH,CI
 DOUBLE PRECISION,DIMENSION(3):: DIPS
 DOUBLE PRECISION,DIMENSION(3*NAT):: GRADS
!
 CHARACTER(LEN=1),DIMENSION(3,NAT) :: dflag
 DOUBLE PRECISION :: Vxyz
 COMMON/INPDYN_VELOCITY/Vxyz(3,100)
 CHARACTER(LEN=1) :: energybound
 COMMON/INPDYN_ENERGYBOUND/energybound
 COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
 COMMON/EHFEN/EHF,EN
 INTEGER :: SIZE_ENV,NBAS,IGTYP
 DOUBLE PRECISION :: ENV(SIZE_ENV)
 INTEGER :: ATM(6,NAT), BAS(8,NBAS)
 INTEGER :: iatom, init
 DOUBLE PRECISION,DIMENSION(NAT)  :: mass
 DOUBLE PRECISION,DIMENSION(3,NAT) :: r0,v0,a0,r1,v1,a1
 DOUBLE PRECISION :: t, tres, dt, tmax
 DOUBLE PRECISION :: Epot, Ekin, Etot
 DOUBLE PRECISION :: Epot_previous_step
 CHARACTER(LEN=1) :: resflag,velflag
 CHARACTER(LEN=6) :: integrator
 CHARACTER(LEN=50) :: line
 CHARACTER(LEN=2) :: atom
 real(kind=8),parameter :: amtoau = 1.d0/0.52917720859d0 
 real(kind=8),parameter :: fstoau = 41.341373336561364d0
 real(kind=8),parameter :: amutoau = 1822.888485540950d0
  write(6,'(/1X,32A)')' Molecular Dynamics Calculation '
  write(6,'(1X,32(1H*),/)')
  if(nat>100)then  ! due to COMMON/INPDYN_VELOCITY/Vxyz(3,100)
   write(6,*)'Stop: Molecular Dynamics is not implemented for NAT>100'
   STOP
  endif
  dflag(1:3,1:nat) = 'T'
  call NAMELIST_DYNINP(nat,resflag,velflag,dt,tmax,integrator)
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
  write(6,*)'Note: t0, Cxyz & Vxyz are redefined according to last xyz'
  Cxyz = Cxyz*amtoau
 endif 
 Vxyz=Vxyz*(amtoau/fstoau) ! Velocity in au/au
 dt=dt*fstoau ! dt in au
 tmax=tmax*fstoau ! tmax in au
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
 write(6,*)
 write(6,*) '-------------------------------------------------------'
 write(6,*) 'Init integrator, time = ', t/fstoau , ' fs'
 write(6,*) '-------------------------------------------------------'
 Epot_previous_step = 0.0d0
 call PROPAGATOR(init,NAT,mass,dt,t,r0,v0,a0,r1,v1,a1,dflag,velflag,    &
                 integrator,NINTEG,IDONTW,IEMOM,NBF,NBFaux,NSHELL,      &
                 NPRIMI,ZAN,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,      &
                 INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,     &
                 CG,CH,CI,GRADS,IRUNTYP,DIPS,IZCORE,SIZE_ENV,ENV,ATM,   &
                 NBAS,BAS,IGTYP,KTYPEaux,NSHELLaux,Epot_previous_step)
 Epot = EELEC + EN
 Ekin = 0.d0

 do i=1,nat
  Ekin = Ekin + 0.5d0*mass(i)                                           &
       * ( v1(1,i)*v1(1,i)+v1(2,i)*v1(2,i)+v1(3,i)*v1(3,i) )
 enddo
 Etot = Epot + Ekin
 Epot_previous_step = Epot
 call XYZWRITER(init,nat,zan,izcore,t,r1,v1,Epot,Ekin,Etot)
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
 init = 1

 do while (abs(t) < tmax)
  write(6,*)
  write(6,*) '---------------------------------------------------'
  write(6,*) ' Propagation time = ', t/fstoau + dt/fstoau, ' fs'
  write(6,*) '---------------------------------------------------'
  call PROPAGATOR(init,nat,mass,dt,t,r0,v0,a0,r1,v1,a1,dflag,velflag,   &
                  integrator,NINTEG,IDONTW,IEMOM,NBF,NBFaux,NSHELL,     &
                  NPRIMI,ZAN,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,     &
                  INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,    &
                  CG,CH,CI,GRADS,IRUNTYP,DIPS,IZCORE,SIZE_ENV,ENV,      &
                  ATM,NBAS,BAS,IGTYP,KTYPEaux,NSHELLaux,Epot_previous_step)
   Epot = EELEC + EN
   Ekin = 0.d0
   do i=1,nat
    Ekin = Ekin + 0.5d0*mass(i)*( v1(1,i)*v1(1,i) &
                + v1(2,i)*v1(2,i)+v1(3,i)*v1(3,i) )
   enddo
   Etot = Epot + Ekin
   Epot_previous_step = Epot
   trealpart = t/fstoau - int(t/fstoau+1.0d-03)
   if ( trealpart < 0.1*dt/fstoau )  &
   call XYZWRITER(init,nat,zan,izcore,t,r1,v1,Epot,Ekin,Etot)      
   OPEN(32,FILE='endyn.dat',STATUS='OLD',POSITION='APPEND',             &
           ACTION='WRITE',FORM='FORMATTED',ACCESS='SEQUENTIAL')
   write(32,'(4(2x,f14.8))') t/fstoau,Epot,Ekin,Etot
   close(32)
   if(t>0.1*tmax)then
    call MAXNUCDIST(nat,r1,distmax)
    if(distmax>distmax0)then
     write(6,*)'Stop: Maximum internuclear distance > 10 Ang'
     return
    end if 
   end if
 enddo
end

SUBROUTINE MAXNUCDIST(nat,Cxyz,distmax)
 IMPLICIT NONE
 INTEGER,INTENT(IN) :: NAT
 DOUBLE PRECISION,INTENT(OUT) :: distmax
 DOUBLE PRECISION,DIMENSION(3,NAT),INTENT(IN) :: Cxyz
 INTEGER :: I,J
 DOUBLE PRECISION :: RR 
 distmax = 0.0d0
 do I = 1,NAT
  do J = I+1,NAT
   RR = (Cxyz(1,I)-Cxyz(1,J))**2 + (Cxyz(2,I)-Cxyz(2,J))**2 +           &                       
        (Cxyz(3,I)-Cxyz(3,J))**2
   RR = DSQRT(RR)
  if(RR>distmax) distmax = RR
  end do
 end do
END

SUBROUTINE PROPAGATOR(init,nat,mass,dt,t,r0,v0,a0,r1,v1,a1,dflag,       &
                      velflag,integrator,NINTEG,IDONTW,IEMOM,NBF,       &
                      NBFaux,NSHELL,NPRIMI,ZAN,IAN,IMIN,IMAX,KSTART,    &
                      KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,    &
                      C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,DIPS, &
                      IZCORE,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,KTYPEaux,  &
                      NSHELLaux,EpotOldStep)
 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
 COMMON/EHFEN/EHF,EN
 INTEGER, intent(in) :: init, nat
 INTEGER :: NINTEG,IDONTW,IEMOM,IRUNTYP
 INTEGER,DIMENSION(NAT) :: IAN,IMIN,IMAX,IZCORE
 INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
 INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
 INTEGER,DIMENSION(NSHELLaux) :: KTYPEaux
 INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
 CHARACTER(LEN=1), dimension (3,nat), intent(in) :: dflag
 CHARACTER(LEN=1), intent(in) :: velflag
 CHARACTER(LEN=6), intent(in) :: integrator
 DOUBLE PRECISION, dimension (nat), intent(in) :: mass
 DOUBLE PRECISION, dimension (3,nat), intent(inout) :: r0,v0,a0,r1,v1,a1
 DOUBLE PRECISION, intent(in)  :: dt
 DOUBLE PRECISION, intent(inout)  :: t
 DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
 DOUBLE PRECISION,DIMENSION(NPRIMI)::C1,C2,EX,CS,CP,CD,CF,CG,CH,CI
 DOUBLE PRECISION,DIMENSION(3,NAT):: GRADS
 DOUBLE PRECISION,DIMENSION(3) :: DIPS
 DOUBLE PRECISION,INTENT(IN) :: EpotOldStep
 INTEGER :: SIZE_ENV,NBAS,IGTYP
 DOUBLE PRECISION :: ENV(SIZE_ENV)
 INTEGER :: ATM(6,NAT), BAS(8,NBAS)
 INTEGER :: i
 DOUBLE PRECISION, dimension (3,NAT) :: r2,v2,a2,forces
 DOUBLE PRECISION :: EpotOld, EpotNew
 real(kind=8), parameter :: fstoau = 41.341373336561364d0   
 EpotOld = 0.0d0
 EpotNew = 0.0d0
 select case (init)
  case (0)
   CALL PES(t,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,ZAN,r0,  &
            IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,  &
            ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,DIPS,  &
            IZCORE,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,KTYPEaux,NSHELLaux)
   EpotOld = EELEC + EN
   forces=-GRADS
   EpotNew = EELEC + EN

   if (dt == 0) then
    goto 10
   endif
   do i=1,nat
    a0(:,i)=forces(:,i)/mass(i)
   enddo
   do i=1,nat
    r1(:,i)=r0(:,i) + dt*v0(:,i) + 0.5d0*dt*dt*a0(:,i)
   enddo
  write(6,*)
  write(6,*) '---------------------------------------------------'
  write(6,*) ' Propagation time = ', t/fstoau + dt/fstoau, ' fs'
  write(6,*) '---------------------------------------------------'

   CALL PES(t+dt,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,ZAN,  &
            r1,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,    &
            KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,  &
            DIPS,IZCORE,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,KTYPEaux,       &
            NSHELLaux)
   
   forces=-GRADS
   EpotNew = EELEC + EN
   do i=1,nat
    a1(:,i)=forces(:,i)/mass(i)
   enddo
   do i=1,nat
    v1(:,i)=v0(:,i) + 0.5d0*dt*(a0(:,i) + a1(:,i))
   enddo
   call ENFORCE_ENERGY_BOUND(t,dt,nat,mass,r0,v0,a0,r1,v1,a1,EpotOld,   &
                             EpotNew,NINTEG,IDONTW,IEMOM,NBF,NBFaux,    &
                             NSHELL,NPRIMI,ZAN,IAN,IMIN,IMAX,KSTART,    &
                             KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,  &
                             ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,  &
                             IRUNTYP,DIPS,IZCORE,SIZE_ENV,ENV,ATM,NBAS, &
                             BAS,IGTYP,KTYPEaux,NSHELLaux)
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
   t = t + dt

  case default
   SELECT CASE(integrator)
   CASE('VERLET')
    do i=1,nat
     r2(:,i)=r1(:,i) + dt*v1(:,i) + 0.5d0*dt*dt*a1(:,i)
    enddo
   CASE('BEEVER')
    do i=1,nat
     r2(:,i)=r1(:,i) + dt*v1(:,i) + dt*dt*(4.0d0*a1(:,i)-a0(:,i))/6.0d0
    enddo
   CASE DEFAULT
    write(6,*)'Stop: Unknown DYN integrator = ',integrator
    stop
   END SELECT
   CALL PES(t+dt,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,ZAN,  &
            r2,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,    &
            KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,  &
            DIPS,IZCORE,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,KTYPEaux,       &
            NSHELLaux)
   
   forces=-GRADS
   EpotNew = EELEC + EN
   do i=1,nat
    a2(:,i)=forces(:,i)/mass(i)
   enddo
!
   SELECT CASE(integrator)
   CASE('VERLET')
    do i=1,nat
     v2(:,i) = v1(:,i) + 0.5d0*dt*(a1(:,i) + a2(:,i))
    enddo
    EpotOld = EpotOldStep
    call ENFORCE_ENERGY_BOUND(t,dt,nat,mass,r1,v1,a1,r2,v2,a2,EpotOld,  &
                              EpotNew,NINTEG,IDONTW,IEMOM,NBF,NBFaux,   &
                              NSHELL,NPRIMI,ZAN,IAN,IMIN,IMAX,KSTART,   &
                              KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH, &
                              ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS, &
                              IRUNTYP,DIPS,IZCORE,SIZE_ENV,ENV,ATM,NBAS,&
                              BAS,IGTYP,KTYPEaux,NSHELLaux)
   CASE('BEEVER')
    do i=1,nat
     v2(:,i)=v1(:,i)+dt*(2.0d0*a2(:,i)+5.0d0*a1(:,i)-a0(:,i))/6.0d0
    enddo
    EpotOld = EpotOldStep
    call ENFORCE_ENERGY_BOUND(t,dt,nat,mass,r1,v1,a1,r2,v2,a2,EpotOld,  &
                              EpotNew,NINTEG,IDONTW,IEMOM,NBF,NBFaux,   &
                              NSHELL,NPRIMI,ZAN,IAN,IMIN,IMAX,KSTART,   &
                              KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH, &
                              ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS, &
                              IRUNTYP,DIPS,IZCORE,SIZE_ENV,ENV,ATM,NBAS,&
                              BAS,IGTYP,KTYPEaux,NSHELLaux)
   CASE DEFAULT
    write(6,*)'Stop: Unknown DYN integrator = ',integrator
    stop
   END SELECT
   r0=r1
   v0=v1
   a0=a1

   r1=r2 
   v1=v2
   a1=a2
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
    t = t + dt

10 end select

end

SUBROUTINE ENFORCE_ENERGY_BOUND(t,dt,nat,mass,r_old,v_old,a_old,        &
                                r_new,v_new,a_new,E_old,E_new,          &
                                NINTEG,IDONTW,IEMOM,NBF,NBFaux,         &
                                NSHELL,NPRIMI,ZAN,IAN,IMIN,IMAX,        &
                                KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,      &
                                KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,      &
                                CD,CF,CG,CH,CI,GRADS,IRUNTYP,DIPS,      &
                                IZCORE,SIZE_ENV,ENV,ATM,NBAS,BAS,       &
                                IGTYP,KTYPEaux,NSHELLaux)
!-----------------------------------------------------------------------
! Keep each MD step within the total energy available from the previous
! state by detecting significant PES jumps and rescaling trial velocities
! stop the propagation if the new point would require negative Ekin.
!-----------------------------------------------------------------------
 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 INTEGER,INTENT(IN) :: nat
 INTEGER :: NINTEG,IDONTW,IEMOM,NBF,NBFaux
 INTEGER :: NSHELL,NPRIMI,IRUNTYP
 INTEGER,DIMENSION(NAT) :: IAN,IMIN,IMAX,IZCORE
 INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
 INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
 INTEGER,DIMENSION(NSHELLaux) :: KTYPEaux
 INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
 INTEGER :: SIZE_ENV,NBAS,IGTYP
 INTEGER :: ATM(6,NAT), BAS(8,NBAS)
 INTEGER :: i
 LOGICAL :: IsRealEvent
 LOGICAL :: CoherentResidual
 DOUBLE PRECISION,INTENT(IN) :: t, dt, E_old
 DOUBLE PRECISION,INTENT(INOUT) :: E_new
 DOUBLE PRECISION,DIMENSION(nat),INTENT(IN) :: mass, ZAN
 DOUBLE PRECISION,DIMENSION(3,nat),INTENT(IN) :: r_old,v_old
 DOUBLE PRECISION,DIMENSION(3,nat),INTENT(IN) :: a_old
 DOUBLE PRECISION,DIMENSION(3,nat),INTENT(INOUT) :: r_new,v_new
 DOUBLE PRECISION,DIMENSION(3,nat),INTENT(INOUT) :: a_new
 DOUBLE PRECISION,DIMENSION(NPRIMI) :: C1,C2,EX,CS,CP
 DOUBLE PRECISION,DIMENSION(NPRIMI) :: CD,CF,CG,CH,CI
 DOUBLE PRECISION,DIMENSION(3,NAT) :: GRADS
 DOUBLE PRECISION,DIMENSION(3) :: DIPS
 DOUBLE PRECISION :: ENV(SIZE_ENV)
 DOUBLE PRECISION,DIMENSION(3) :: grad_i
 DOUBLE PRECISION :: Kold,EtotAvail,Etrial,Ktrial,Ktarget
 DOUBLE PRECISION :: EtotTol, TinyKinetic, Scale
 DOUBLE PRECISION :: GradWork,PotentialJump,ResidualJump
 DOUBLE PRECISION :: jumptol,ResidualTol
 DOUBLE PRECISION :: SmallJumpTol,AccumResidualTol
 DOUBLE PRECISION,SAVE :: AccumResidualJump = 0.0d0
 DOUBLE PRECISION :: ResidualSign
 CHARACTER(LEN=1) :: energybound
 COMMON/INPDYN_ENERGYBOUND/energybound
 INTEGER,SAVE :: ResidualTrend = 0
 COMMON/INPDYN_EVENTTOL/jumptol
!-----------------------------------------------------------------------
 if (energybound /= 'T') return
 EtotTol = 1.0d-8
 TinyKinetic = 1.0d-20
 ResidualTol = 5.0d-4
 SmallJumpTol = 5.0d-5
 AccumResidualTol = 3.0d-4
 Kold = 0.0d0
 GradWork = 0.0d0
!
 do i=1,nat
  Kold = Kold + 0.5d0*mass(i) * ( v_old(1,i)*v_old(1,i)                 &
                                + v_old(2,i)*v_old(2,i)                 &
                                + v_old(3,i)*v_old(3,i) )
  grad_i(1) = -mass(i)*a_old(1,i)
  grad_i(2) = -mass(i)*a_old(2,i)
  grad_i(3) = -mass(i)*a_old(3,i)
  GradWork = GradWork + grad_i(1)*(r_new(1,i)-r_old(1,i))               &
                      + grad_i(2)*(r_new(2,i)-r_old(2,i))               &
                      + grad_i(3)*(r_new(3,i)-r_old(3,i))
 enddo
 EtotAvail = E_old + Kold
 PotentialJump = E_new - E_old
 ResidualJump = PotentialJump - GradWork
 CoherentResidual = .FALSE.
 if (ABS(PotentialJump) > SmallJumpTol .and.                            &
     ABS(ResidualJump) > SmallJumpTol) then
  ResidualSign = DSIGN(1.0d0,ResidualJump)
  if (PotentialJump*ResidualJump > 0.0d0) then
   CoherentResidual = .TRUE.
   if (ResidualTrend == 0) then
    ResidualTrend = NINT(ResidualSign)
    AccumResidualJump = ResidualJump
   elseif (ResidualTrend == NINT(ResidualSign)) then
    AccumResidualJump = AccumResidualJump + ResidualJump
   else
    ResidualTrend = NINT(ResidualSign)
    AccumResidualJump = ResidualJump
   endif
  else
   AccumResidualJump = 0.0d0
   ResidualTrend = 0
  endif
 else
  AccumResidualJump = 0.0d0
  ResidualTrend = 0
 endif
 Ktrial = 0.0d0
 do i=1,nat
  Ktrial = Ktrial + 0.5d0*mass(i) * ( v_new(1,i)*v_new(1,i)             &
                                    + v_new(2,i)*v_new(2,i)             &
                                    + v_new(3,i)*v_new(3,i) )
 enddo
 Etrial = E_new + Ktrial
 IsRealEvent = .FALSE.
 if (ABS(PotentialJump) > jumptol) then
  IsRealEvent = .TRUE.
 elseif (ABS(ResidualJump) > ResidualTol) then
  IsRealEvent = .TRUE.
 elseif (CoherentResidual .and.                                         &
         ABS(AccumResidualJump) > AccumResidualTol) then
  IsRealEvent = .TRUE.
 endif
 if (IsRealEvent) then
  AccumResidualJump = 0.0d0
  ResidualTrend = 0
  if (Ktrial <= TinyKinetic) then
   write(6,'(/1X,A)')                                                   &
        'Stop: event detected but trial kinetic energy is negligible'
   write(6,'(1X,A,ES20.10,A,ES20.10,A,ES20.10)')                        &
        '  E_old=',E_old,'  E_new=',E_new,'  EtotAvail=',EtotAvail
   stop
  endif
  Ktarget = EtotAvail - E_new
  if (Ktarget < -EtotTol) then
   write(6,'(/1X,A)')                                                   &
        'Stop: event step leaves no accessible kinetic energy'
   write(6,'(1X,A,ES20.10,A,ES20.10,A,ES20.10)')                        &
        '  E_old=',E_old,'  E_new=',E_new,'  EtotAvail=',EtotAvail
   write(6,'(1X,A,ES20.10)') '  Ktarget=',Ktarget
   stop
  endif
  if (Ktarget < 0.0d0) Ktarget = 0.0d0
  Scale = DSQRT(Ktarget/Ktrial)
  do i=1,nat
   v_new(:,i) = Scale*v_new(:,i)
  enddo
  return
 endif
 if (Etrial <= EtotAvail + EtotTol) return
 Ktarget = EtotAvail - E_new
 if (Ktarget < -EtotTol) then
  write(6,'(/1X,A)')                                                    &
       'Stop: step exceeds available total energy'
  write(6,'(1X,A,ES20.10,A,ES20.10,A,ES20.10)')                         &
       '  E_old=',E_old,'  E_new=',E_new,'  EtotAvail=',EtotAvail
  write(6,'(1X,A,ES20.10,A,ES20.10)')                                   &
       '  Ktrial=',Ktrial,'  EtotTrial=',Etrial
  write(6,'(1X,A,ES20.10)') '  Ktarget=',Ktarget
  stop
 endif
 if (Ktrial <= TinyKinetic) then
  write(6,'(/1X,A)')                                                    &
       'Stop: step needs velocity rescaling with negligible kinetic energy'
  write(6,'(1X,A,ES20.10,A,ES20.10,A,ES20.10)')                         &
       '  E_old=',E_old,'  E_new=',E_new,'  EtotAvail=',EtotAvail
  write(6,'(1X,A,ES20.10,A,ES20.10)')                                   &
       '  Ktrial=',Ktrial,'  EtotTrial=',Etrial
  stop
 endif
 if (Ktarget < 0.0d0) Ktarget = 0.0d0
 Scale = DSQRT(Ktarget/Ktrial)
 do i=1,nat
  v_new(:,i) = Scale*v_new(:,i)
 enddo
 return
END

SUBROUTINE PES(t,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,ZAN,  &
               r0,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,      &
               KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,        &
               GRADS,IRUNTYP,DIPS,IZCORE,SIZE_ENV,ENV,ATM,NBAS,BAS,     &
               IGTYP,KTYPEaux,NSHELLaux)
 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!-----------------------------------------------------------------------
! interface
!-----------------------------------------------------------------------
 INTEGER ::NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,IRUNTYP
 INTEGER,DIMENSION(NAT) :: IAN,IMIN,IMAX,IZCORE
 INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
 INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
 INTEGER,DIMENSION(NSHELLaux) :: KTYPEaux
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
 character (len=10) :: fileIn
 character (len=22) :: MLD_FILE
 character (len=8) :: str1

 INTEGER :: i, MOLDEN_SAVE
 real(kind=8), parameter :: fstoau = 41.341373336561364d0  
 INTEGER :: SIZE_ENV,NBAS,IGTYP
 DOUBLE PRECISION :: ENV(SIZE_ENV)
 INTEGER :: ATM(6,NAT), BAS(8,NBAS)
 do i=1,nat
  ENV(20+3*(i-1)+1) = r0(1,i)
  ENV(20+3*(i-1)+2) = r0(2,i)
  ENV(20+3*(i-1)+3) = r0(3,i) 
 enddo
!
 MOLDEN_SAVE = MOLDEN
 if(snapshot=='T')then
  MOLDEN = 1
  write(str1,'(F8.2)')t/fstoau
 endif 
 if(snapshot=='T')then
  MLD_FILE = 'snapshot-'//trim(adjustl(str1))//".mld"
  OPEN(17,FILE=MLD_FILE,STATUS='UNKNOWN',FORM='FORMATTED',             &
          ACCESS='SEQUENTIAL') 
 endif
 fileIn = 'GCF'
 OPEN(3,FILE=fileIn,FORM='FORMATTED',ACCESS='SEQUENTIAL')
 CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,ZAN,r0, &
               IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,    &
               KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,       &
               IRUNTYP,DIPS,IZCORE,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,     &
               KTYPEaux,NSHELLaux,0,1)
 CLOSE(3)
 IF(snapshot=='T')CLOSE(17)
 MOLDEN = MOLDEN_SAVE
 
end

SUBROUTINE XYZWRITER(init,nat,zan,izcore,tau,rau,vau,Epot,Ekin,Etot)
 implicit none
 logical DynFileExists,DynlFileExists
 integer,intent(in) :: init,nat
 integer,dimension (nat),intent(in) :: izcore
 real(kind=8),dimension (nat),intent(in) :: zan
 real(kind=8),intent(in) :: tau
 real(kind=8),dimension (3,nat),intent(in) :: rau,vau
 real(kind=8),intent(in) :: Epot,Ekin,Etot
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
 integer :: i,iznuc
 real(kind=8) :: tfs
 real(kind=8),dimension (3,nat) :: ram,vam
 real(kind=8),parameter :: amtoau=1.d0/0.52917720859d0 
 real(kind=8),parameter :: fstoau=41.341373336561364d0
 inquire(FILE='DYN.xyz',EXIST=DynFileExists)
 if(DynFileExists)then 
  OPEN(30,FILE='DYN.xyz',STATUS='OLD',POSITION='APPEND',ACTION='WRITE', &
          FORM='FORMATTED',ACCESS='SEQUENTIAL')
 else
  OPEN(30,FILE='DYN.xyz',STATUS='NEW',ACTION='WRITE',                   &
          FORM='FORMATTED',ACCESS='SEQUENTIAL')
 end if
!
 if (init /= 0) then
  inquire(FILE='DYNl.xyz',EXIST=DynlFileExists)
  if (DynlFileExists) then
   OPEN(31,FILE='DYNl.xyz',STATUS='OLD',POSITION='REWIND',              &
           ACTION='WRITE',FORM='FORMATTED',ACCESS='SEQUENTIAL')
  else
   OPEN(31,FILE='DYNl.xyz',STATUS='NEW',ACTION='WRITE',                 &
           FORM='FORMATTED',ACCESS='SEQUENTIAL')
  endif
 endif
 tfs=tau/fstoau 
 write(30,*) nat
 write(30,'(4(a6),4(f14.8))') 't,','Epot,','Ekin,','Etot,',            &
                              tfs,Epot,Ekin,Etot
 do i=1,nat
  ram(:,i)=rau(:,i)/amtoau
  vam(:,i)=vau(:,i)*(fstoau/amtoau)
  iznuc = int(zan(i)) + izcore(i)
  write(30,'(a4, 6(1x,f14.8))')ATMLAB(iznuc),ram(:,i),vam(:,i)
 enddo

 if (init /= 0) then
  rewind 31
  write(31,*) nat
  write(31,*) tfs
  do i=1,nat
   write(31,'(a4, 6(1x,f14.8))')ATMLAB(int(zan(i))+izcore(i)),ram(:,i),vam(:,i)
  enddo
 endif

 close(30)
 if (init /= 0) close(31)

end

