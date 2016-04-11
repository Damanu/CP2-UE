!**********************************************************************!
!
! File: gmclj.f90
!
! Create random ("gas") initial configuration for NVT-Monte Carlo of 
! Lennard-Jonesium
!
! 13-Apr-1997 (MN)
! 17-Apr-2012
!
!**********************************************************************!

program gmclj

!**********************************************************************!

  implicit none

  integer,parameter::iout=2

  character(len=80)::fname
  integer::i,n,nran,nt,ntjob,ntprint,ntskip
  integer,dimension(:),allocatable::iran
  real::c,disp,dr,ran,rho,t
  real,dimension(:),allocatable::x,y,z

! User input

  write(unit=*,fmt="(a)",advance="no") "              n="
  read(unit=*,fmt=*) n
  write(unit=*,fmt="(a)",advance="no") "            rho="
  read(unit=*,fmt=*) rho
  write(unit=*,fmt="(a)",advance="no") "              t="
  read(unit=*,fmt=*) t
  write(unit=*,fmt="(a)",advance="no") "           disp="
  read(unit=*,fmt=*) disp
  write(unit=*,fmt="(a)",advance="no") "             dr="
  read(unit=*,fmt=*) dr
  write(unit=*,fmt="(a)",advance="no") "         ntskip="
  read(unit=*,fmt=*) ntskip
  write(unit=*,fmt="(a)",advance="no") " ntprint/ntskip="
  read(unit=*,fmt=*) ntprint
  write(unit=*,fmt="(a)",advance="no") "   ntjob/ntskip="
  read(unit=*,fmt=*) ntjob
  write(unit=*,fmt="(a)",advance="no") "          fname=[mclj_in.dat] "
  read(unit=*,fmt="(a)") fname
  if(fname=="") then
    fname="mclj_in.dat"
  end if

! Allocate arrays

  call random_seed(size=nran)

  allocate(iran(nran),x(n),y(n),z(n))

! Random positions

  c=(n/rho)**(1.0/3.0)
  do i=1,n
    call random_number(ran)
    x(i)=c*(ran-0.5)
    call random_number(ran)
    y(i)=c*(ran-0.5)
    call random_number(ran)
    z(i)=c*(ran-0.5)
  end do

! RNG seed

  call random_seed(get=iran)

! Write startup file

  nt=0

  open(unit=iout,file=fname,status="replace",action="write", &
       form="unformatted")

  write(unit=iout) n,nt,ntjob,ntprint,ntskip,iran
  write(unit=iout) disp,dr,rho,t
  write(unit=iout) x,y,z

  close(unit=iout)

! Deallocate arrays

  deallocate(iran,x,y,z)

end program gmclj

!**********************************************************************!
