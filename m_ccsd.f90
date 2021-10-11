! Reference for Equations:
!  Stanton, et al. J. Chem. Phys 94 (6), 1991
!  Crawford, Schaefer (2000)  DOI: 10.1002/9780470125915.ch2 Eq. 134
module m_ccsd
implicit none
double precision,parameter::zero=0.0d0
double precision,parameter::half=5.0d-1
double precision,parameter::one_fourth=2.5d-1
integer,parameter::iunit=2341,iunit2=2142
logical::qnewton
integer::Nocc,Mbasis2,iter
double precision::tol8=1.0d-8
double precision,dimension(:,:),allocatable::ts,tsnew,Fae,Fmi,Fme,FockM
double precision,dimension(:,:,:,:),allocatable::td,tdnew,Wmnij,Wabef,Wmbej

private::Nocc,Mbasis2,zero,half,one_fourth,ts,tsnew,Fae,Fmi,Fme,td,tdnew,Wmnij,Wabef,Wmbej,FockM,tol8
private::ccsd_update_interm,ccsd_update_t1_t2,spin_int,taus,tau,slbasis,iter,qnewton
public::ccsd_init,ccsd_read_guess,ccsd_write_ts,ccsd_update_ts,ccsd_energy,ccsd_clean,ccsd_en_nof

contains

!!-----------------------------------------------------------
!! Public
!!-----------------------------------------------------------
subroutine ccsd_init(Mbasis,Nocc_in,qnewton_in,CCSDint,FockM_in,ERImol,ERImol2)
implicit none
! Arguments
logical,intent(in)::qnewton_in,CCSDint
integer,intent(in)::Mbasis,Nocc_in
double precision,dimension(:,:),intent(in)::FockM_in
double precision,dimension(:,:,:,:),intent(in)::ERImol
double precision,dimension(:,:,:,:),optional,intent(in)::ERImol2
! Local variables
double precision::value1
integer::i,j,a,b
! Procedures
iter=0
qnewton=qnewton_in
Mbasis2=Mbasis*2 ! Spin-less (size of the basis set) -> spin-with
Nocc=Nocc_in*2   ! Spin-less -> spin-with
write(*,*) ' ' 
write(*,'(a,i5)') ' Number of spin orbitals ',Mbasis2
write(*,'(a,i5)') ' Number of spin doubly occ. orbitals ',Nocc
write(*,*) ' ' 
allocate(ts(Nocc,Nocc+1:Mbasis2),tsnew(Nocc,Nocc+1:Mbasis2),FockM(Mbasis2,Mbasis2))
allocate(td(Nocc,Nocc,Nocc+1:Mbasis2,Nocc+1:Mbasis2),tdnew(Nocc,Nocc,Nocc+1:Mbasis2,Nocc+1:Mbasis2))
allocate(Fae(Nocc+1:Mbasis2,Nocc+1:Mbasis2),Fmi(Nocc,Nocc),Fme(Nocc,Nocc+1:Mbasis2))
allocate(Wmnij(Nocc,Nocc,Nocc,Nocc),Wabef(Nocc+1:Mbasis2,Nocc+1:Mbasis2,Nocc+1:Mbasis2,Nocc+1:Mbasis2))
allocate(Wmbej(Nocc,Nocc+1:Mbasis2,Nocc+1:Mbasis2,Nocc))
ts=zero;tsnew=zero;td=zero;tdnew=zero;Fae=zero;Fmi=zero;Fme=zero;
Wmnij=zero;Wabef=zero;Wmbej=zero;FockM=zero;
! Initialize FockM
do i=1,Mbasis2
 do j=1,Mbasis2
  if(mod(i,2)==mod(j,2)) then
   FockM(j,i)=FockM_in(slbasis(j),slbasis(i))
  endif
 enddo
enddo
! Initialize td
do i=1,Nocc
 do j=1,Nocc
  do a=Nocc+1,Mbasis2
   do b=Nocc+1,Mbasis2
    td(i,j,a,b)=spin_int(i,j,a,b,ERImol)/(FockM(i,i)+FockM(j,j)-FockM(a,a)-FockM(b,b)+tol8)
   enddo
  enddo
 enddo
enddo
! Print FockM and ERIs?
if(CCSDint) then
open(unit=iunit,form='formatted',file='fockm')
do i=1,Mbasis2
 if(abs(FockM(i,i))>tol8) then
  write(iunit,'(i5,i5,f15.10)') i,i,FockM(i,i)
 endif
enddo
close(iunit)
open(unit=iunit,form='formatted',file='terimol')
open(unit=iunit2,form='formatted',file='ferimol')
do i=1,Mbasis2
 do j=1,Mbasis2
  do a=1,Mbasis2
   do b=1,Mbasis2
    if(mod(i,2)==mod(a,2) .and. mod(j,2)==mod(b,2)) then
     if(abs(ERImol(slbasis(i),slbasis(j),slbasis(a),slbasis(b)))>tol8) then
      write(iunit,'(i5,i5,i5,i5,f15.10)') i,j,a,b,ERImol(slbasis(i),slbasis(j),slbasis(a),slbasis(b))
     endif
     if(present(ERImol2)) then
      if(abs(ERImol2(slbasis(i),slbasis(j),slbasis(a),slbasis(b)))>tol8) then
       write(iunit2,'(i5,i5,i5,i5,f15.10)') i,j,a,b,ERImol2(slbasis(i),slbasis(j),slbasis(a),slbasis(b))
      endif
     endif
    endif
   enddo
  enddo
 enddo
enddo
close(iunit)
close(iunit2)
endif
end subroutine ccsd_init

subroutine ccsd_read_guess()
implicit none
! Arguments
! Local variables
integer::i,j,a,b,istat
double precision::value1
! Procedures
td=zero
ts=zero
open(unit=iunit,form='unformatted',file='t_amps',iostat=istat,status='old')
if(istat==0) then
 do
  read(iunit,iostat=istat) i,j,a,b,value1
  if(istat==0) then 
   if(i/=0 .and. a/=0) then
    if(j/=0 .and. b/=0) then
     td(i,j,a,b)=value1
    else 
     ts(i,a)=value1
    endif
   else 
    exit 
   endif
  else
   exit
  endif
 enddo
endif
close(iunit)
end subroutine ccsd_read_guess

subroutine ccsd_write_ts(filep)
implicit none
! Arguments
character(len=50),intent(in)::filep 
! Local variables
integer::i,j,a,b
! Procedures
open(unit=iunit,form='unformatted',file=filep)
do i=1,Nocc
 do j=1,Nocc
  do a=Nocc+1,Mbasis2
   do b=Nocc+1,Mbasis2
    if(abs(td(i,j,a,b))>tol8) write(iunit) i,j,a,b,td(i,j,a,b)
   enddo
  enddo
 enddo
enddo
do i=1,Nocc
 do a=Nocc+1,Mbasis2
  if(abs(ts(i,a))>tol8) write(iunit) i,0,a,0,ts(i,a)
 enddo
enddo
write(iunit) 0,0,0,0,zero
close(iunit)
end subroutine ccsd_write_ts

subroutine ccsd_update_ts(deltaE,ERImol)
implicit none
! Arguments
double precision,intent(in)::deltaE
double precision,dimension(:,:,:,:),intent(in)::ERImol
! Local variables
character(len=50)::filep 
integer::i,j,a,b
double precision::tol4=1d-4,tol8=1d-8
! Procedures
iter=iter+1
call ccsd_update_interm(ERImol)
call ccsd_update_t1_t2(ERImol)
do i=1,Nocc
 do j=1,Nocc
  do a=Nocc+1,Mbasis2
   if(j==1) then                 
    if(abs(tsnew(i,a))<tol8) then
     tsnew(i,a)=zero             
    endif                        
   endif                         
   do b=Nocc+1,Mbasis2
    if(abs(tdnew(i,j,a,b))<tol8) then
     tdnew(i,j,a,b)=zero
    endif
   enddo
  enddo
 enddo
enddo
ts=tsnew
td=tdnew
end subroutine ccsd_update_ts

function ccsd_energy(ERImol)
implicit none
! Arguments
double precision,dimension(:,:,:,:),intent(in)::ERImol
! Local variables
integer::i,j,a,b
double precision::ccsd_energy
! Procedures
ccsd_energy=zero
do i=1,Nocc
 do a=Nocc+1,Mbasis2
  ccsd_energy=ccsd_energy+FockM(i,a)*ts(i,a)
  do j=1,Nocc
   do b=Nocc+1,Mbasis2
    ccsd_energy=ccsd_energy+one_fourth*spin_int(i,j,a,b,ERImol)*td(i,j,a,b)&
    &+half*spin_int(i,j,a,b,ERImol)*ts(i,a)*ts(j,b)
   enddo
  enddo
 enddo
enddo
end function ccsd_energy

function ccsd_en_nof(Ndocc,Nsocc,ERImol)
implicit none
! Arguments
integer,intent(in)::Ndocc,Nsocc
double precision,dimension(:,:,:,:),intent(in)::ERImol
! Local variables
integer::i,j,k,a,b
double precision::ccsd_en_nof
! Procedures
ccsd_en_nof=zero
do a=Nsocc+1,Mbasis2
 do b=Nsocc+1,Mbasis2
  do i=1,Ndocc
   if(b==Nsocc+1) then
     ccsd_en_nof=ccsd_en_nof+FockM(i,a)*ts(i,a) ! This is 0 because the FockM is diagonal.
   endif
   do j=1,Ndocc
    ccsd_en_nof=ccsd_en_nof+one_fourth*spin_int(i,j,a,b,ERImol)*td(i,j,a,b)&
    &+half*spin_int(i,j,a,b,ERImol)*ts(i,a)*ts(j,b)
   enddo
   do j=Ndocc+1,Nsocc
    ccsd_en_nof=ccsd_en_nof+half*(one_fourth*spin_int(i,j,a,b,ERImol)*td(i,j,a,b)&
    &+half*spin_int(i,j,a,b,ERImol)*ts(i,a)*ts(j,b))
   enddo
  enddo
  do i=Ndocc+1,Nsocc
   do j=1,Ndocc
    ccsd_en_nof=ccsd_en_nof+half*(one_fourth*spin_int(i,j,a,b,ERImol)*td(i,j,a,b)&
    &+half*spin_int(i,j,a,b,ERImol)*ts(i,a)*ts(j,b))
   enddo
   do j=Ndocc+1,Nsocc
    if(mod(i,2)==0) then
     k=i-1
    else
     k=i+1
    endif
    if(j/=i .and. j/=k) then
     ccsd_en_nof=ccsd_en_nof+one_fourth*(one_fourth*spin_int(i,j,a,b,ERImol)*td(i,j,a,b)&
     &+half*spin_int(i,j,a,b,ERImol)*ts(i,a)*ts(j,b))
    endif
   enddo
  enddo
 enddo
enddo
end function ccsd_en_nof

subroutine ccsd_clean()
implicit none
! Arguments
! Local variables
integer::i,iunit=1234,istat
character(len=50)::filep 
! Procedures
deallocate(ts,tsnew,td,tdnew,Fae,Fmi,Fme,FockM)
deallocate(Wmnij,Wabef,Wmbej)
end subroutine ccsd_clean

!!-----------------------------------------------------------
!! Private
!!-----------------------------------------------------------
subroutine ccsd_update_interm(ERImol)
implicit none
! Arguments
double precision,dimension(:,:,:,:),intent(in)::ERImol
! Local variables
integer::i,j,m,n,a,b,e,f
! Procedures
! Fmi and Wmnij
Fmi=zero;Wmnij=zero
do m=1,Nocc
 do i=1,Nocc
  if(i/=m) then
   Fmi(m,i)=FockM(m,i)
  endif
  do n=1,Nocc
   do j=1,Nocc
    Wmnij(m,n,i,j)=spin_int(m,n,i,j,ERImol)
    do e=Nocc+1,Mbasis2
     if(j==1) then
      Fmi(m,i)=Fmi(m,i)+ts(n,e)*spin_int(m,n,i,e,ERImol)
      if(n==1) then
       Fmi(m,i)=Fmi(m,i)+half*ts(i,e)*FockM(m,e)
      endif
     endif
     Wmnij(m,n,i,j)=Wmnij(m,n,i,j)+ts(j,e)*spin_int(m,n,i,e,ERImol)-ts(i,e)*spin_int(m,n,j,e,ERImol)
     do f=Nocc+1,Mbasis2
      Wmnij(m,n,i,j)=Wmnij(m,n,i,j)+one_fourth*tau(e,f,i,j)*spin_int(m,n,e,f,ERImol)
      if(j==1) then
       Fmi(m,i)=Fmi(m,i)+half*taus(e,f,i,n)*spin_int(m,n,e,f,ERImol)
      endif
     enddo
    enddo
   enddo
  enddo
 enddo
enddo
! Fae and Wabef
Fae=zero;Wabef=zero
do a=Nocc+1,Mbasis2
 do e=Nocc+1,Mbasis2
  if(a/=e) then
   Fae(a,e)=FockM(a,e)
  end if
  do b=Nocc+1,Mbasis2
   do f=Nocc+1,Mbasis2
    Wabef(a,b,e,f)=spin_int(a,b,e,f,ERImol)
    do m=1,Nocc
     if(b==Nocc+1) then
      if(f==Nocc+1) then
       Fae(a,e)=Fae(a,e)-half*FockM(m,e)*ts(m,a)
      endif
      Fae(a,e)=Fae(a,e)+ts(m,f)*spin_int(m,a,f,e,ERImol)
     endif
     Wabef(a,b,e,f)=Wabef(a,b,e,f)-ts(m,b)*spin_int(a,m,e,f,ERImol)+ts(m,a)*spin_int(b,m,e,f,ERImol)
     do n=1,Nocc
      Wabef(a,b,e,f)=Wabef(a,b,e,f)+one_fourth*tau(a,b,m,n)*spin_int(m,n,e,f,ERImol)
      if(b==Nocc+1) then
       Fae(a,e)=Fae(a,e)-half*taus(a,f,m,n)*spin_int(m,n,e,f,ERImol)
      endif
     enddo
    enddo
   enddo
  enddo
 enddo
enddo
! Fme and Wmbej
Fme=zero;Wmbej=zero
do m=1,Nocc
 do e=Nocc+1,Mbasis2
  Fme(m,e)=FockM(m,e)
  do b=Nocc+1,Mbasis2
   do j=1,Nocc
    Fme(m,e)=Fme(m,e)+ts(j,b)*spin_int(m,j,e,b,ERImol)
    Wmbej(m,b,e,j)=spin_int(m,b,e,j,ERImol)
    do f=Nocc+1,Mbasis2
     Wmbej(m,b,e,j)=Wmbej(m,b,e,j)+ts(j,f)*spin_int(m,b,e,f,ERImol)
    enddo
    do n=1,Nocc
     Wmbej(m,b,e,j)=Wmbej(m,b,e,j)-ts(n,b)*spin_int(m,n,e,j,ERImol)
     do f=Nocc+1,Mbasis2
      Wmbej(m,b,e,j)=Wmbej(m,b,e,j)-(half*td(j,n,f,b)+ts(j,f)*ts(n,b))*spin_int(m,n,e,f,ERImol)
     enddo
    enddo
   enddo
  enddo
 enddo
enddo
end subroutine ccsd_update_interm

subroutine ccsd_update_t1_t2(ERImol)
implicit none
! Arguments
double precision,dimension(:,:,:,:),intent(in)::ERImol
! Local variables
integer::i,j,m,n,a,b,e,f
logical::ria,rijab
! Procedures
tsnew=zero
tdnew=zero
do a=Nocc+1,Mbasis2
 do b=Nocc+1,Mbasis2
  do i=1,Nocc
   ! T1
   if(b==Nocc+1) then 
    tsnew(i,a)=FockM(i,a)
    do m=1,Nocc
     tsnew(i,a)=tsnew(i,a)-ts(m,a)*Fmi(m,i)
     do e=Nocc+1,Mbasis2
      tsnew(i,a)=tsnew(i,a)+td(i,m,a,e)*Fme(m,e)
      tsnew(i,a)=tsnew(i,a)-ts(m,e)*spin_int(m,a,i,e,ERImol)
      if(m==1) then                                         
       tsnew(i,a)=tsnew(i,a)+ts(i,e)*Fae(a,e)               
      endif                                                 
      do f=Nocc+1,Mbasis2
       tsnew(i,a)=tsnew(i,a)-half*td(i,m,e,f)*spin_int(m,a,e,f,ERImol)
      enddo
      do n=1,Nocc
       tsnew(i,a)=tsnew(i,a)-half*td(m,n,a,e)*spin_int(n,m,e,i,ERImol)
      enddo
     enddo
    enddo
    if(qnewton) then
     ria=tsnew(i,a)-ts(i,a)*(FockM(i,i)-FockM(a,a))
     if(abs(ria)>tol8) then
      tsnew(i,a)=ts(i,a)+ria/(FockM(i,i)-FockM(a,a)+tol8)
     endif
    else
     tsnew(i,a)=tsnew(i,a)/(FockM(i,i)-FockM(a,a)+tol8)
    endif
   endif
   ! T2
   do j=1,Nocc
    tdnew(i,j,a,b)=spin_int(i,j,a,b,ERImol)
    do e=Nocc+1,Mbasis2
     tdnew(i,j,a,b)=tdnew(i,j,a,b)+td(i,j,a,e)*Fae(b,e)-td(i,j,b,e)*Fae(a,e)
     do m=1,Nocc
      tdnew(i,j,a,b)=tdnew(i,j,a,b)-half*td(i,j,a,e)*ts(m,b)*Fme(m,e)+half*td(i,j,b,e)*ts(m,a)*Fme(m,e)
      if(e==Nocc+1) then                                                        
       tdnew(i,j,a,b)=tdnew(i,j,a,b)-td(i,m,a,b)*Fmi(m,j)+ td(j,m,a,b)*Fmi(m,i) 
       tdnew(i,j,a,b)=tdnew(i,j,a,b)-ts(m,a)*spin_int(m,b,i,j,ERImol)+ts(m,b)*spin_int(m,a,i,j,ERImol) 
       do n=1,Nocc                                                              
        tdnew(i,j,a,b)=tdnew(i,j,a,b)+half*tau(a,b,m,n)*Wmnij(m,n,i,j)          
       enddo                                                                    
      end if                                                                    
      tdnew(i,j,a,b)=tdnew(i,j,a,b)-half*td(i,m,a,b)*ts(j,e)*Fme(m,e)+half*td(j,m,a,b)*ts(i,e)*Fme(m,e) 
      tdnew(i,j,a,b)=tdnew(i,j,a,b)+td(i,m,a,e)*Wmbej(m,b,e,j)-ts(i,e)*ts(m,a)*spin_int(m,b,e,j,ERImol) 
      tdnew(i,j,a,b)=tdnew(i,j,a,b)-td(j,m,a,e)*Wmbej(m,b,e,i)+ts(j,e)*ts(m,a)*spin_int(m,b,e,i,ERImol) 
      tdnew(i,j,a,b)=tdnew(i,j,a,b)-td(i,m,b,e)*Wmbej(m,a,e,j)+ts(i,e)*ts(m,b)*spin_int(m,a,e,j,ERImol) 
      tdnew(i,j,a,b)=tdnew(i,j,a,b)+td(j,m,b,e)*Wmbej(m,a,e,i)-ts(j,e)*ts(m,b)*spin_int(m,a,e,i,ERIMol) 
     enddo
     tdnew(i,j,a,b)=tdnew(i,j,a,b)+ts(i,e)*spin_int(a,b,e,j,ERImol)-ts(j,e)*spin_int(a,b,e,i,ERImol)    
     do f=Nocc+1,Mbasis2                                                        
      tdnew(i,j,a,b)=tdnew(i,j,a,b)+half*tau(e,f,i,j)*Wabef(a,b,e,f)            
     enddo                                                                      
    enddo
    tdnew(i,j,a,b)=tdnew(i,j,a,b)/(FockM(i,i)+FockM(j,j)-FockM(a,a)-FockM(b,b)+tol8)
   enddo 
  enddo 
 enddo
enddo
end subroutine ccsd_update_t1_t2

function spin_int(p,q,r,s,ERImol)
implicit none
! Arguments
integer,intent(in)::p,q,r,s
double precision,dimension(:,:,:,:),intent(in)::ERImol
! Local variables
double precision::value1,value2,spin_int
! Procedures
spin_int=zero;value1=zero;value2=zero;
if(mod(p,2)==mod(r,2) .and. mod(q,2)==mod(s,2)) then
  value1=ERImol(slbasis(p),slbasis(q),slbasis(r),slbasis(s))
else
  value1=zero
endif
if(mod(p,2)==mod(s,2) .and. mod(q,2)==mod(r,2)) then
  value2=ERImol(slbasis(p),slbasis(q),slbasis(s),slbasis(r))
else
  value2=zero  
end if
spin_int=value1-value2
end function spin_int

function taus(a,b,i,j)
implicit none
! Arguments
integer,intent(in)::a,b,i,j
! Local variables
double precision::taus
! Procedures
taus=td(i,j,a,b)+half*(ts(i,a)*ts(j,b)-ts(i,b)*ts(j,a))
end function taus

function tau(a,b,i,j)
implicit none
! Arguments
integer,intent(in)::a,b,i,j
! Local variables
double precision::tau
! Procedures
tau=td(i,j,a,b)+ts(i,a)*ts(j,b)-ts(i,b)*ts(j,a)
end function tau

function slbasis(p)
implicit none
! Arguments
integer,intent(in)::p
! Local variables
integer::slbasis
! Procedures
if(mod(p,2)==0) then
  slbasis=p/2
else
  slbasis=(p+1)/2
end if
end function slbasis

end module m_ccsd
