subroutine ERPA(NOCB,NOC,NVI,NORB,CINTRA,CINTER,FOCKm,ERImol,ECd)
implicit double precision (A-H,O-Z)
common/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
common/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
integer :: NOCB,NOC,NVI,NORB
double precision :: ECd
double precision,dimension(NORB) :: CINTRA,CINTER
double precision,dimension(NORB,NORB) :: FOCKm
double precision,dimension(NOC,NBFT,NBF) :: ERImol
ECd=0.0d0

write(*,*) 'MRM, RPA IN HERE!',NOCB,NVI,NORB,NCWO
write(*,*) 'MRM, RPA IN HERE!',NOC,NBFT,NBF,NBF*NBF
do i=1,NORB
 write(*,'(*(f10.5))') FOCKm(i,:)
enddo

ECd=-0.34159d0
end subroutine ERPA
