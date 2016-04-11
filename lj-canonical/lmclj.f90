!**********************************************************************!
!
! File: lmclj.f90
!
! Create initial configuration (fcc lattice) for NVT-Monte Carlo of 
! Lennard-Jonesium
!
! 13-Apr-1997 (MN)
! 17-Apr-2012
!
!**********************************************************************!

program lmclj

!**********************************************************************!

  implicit none

  integer,parameter::iout=2

  character(len=80)::fname
  integer::i,ix,iy,iz,m,n,nran,nt,ntjob,ntprint,ntskip
  integer,dimension(:),allocatable::iran
  real::c,disp,dr,rho,t
  real,dimension(:),allocatable::x,y,z

! User input

  write(unit=*,fmt="(a)",advance="no") "              n="
  read(unit=*,fmt=*) n
  m=0
  do
    m=m+2
    if(m**3==2*n) then
      exit
    else if(m**3>2*n) then
      write(unit=*,fmt="(a)") " lmclj: not a magic number"
      stop
    end if
  end do

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

! Fcc lattice

  c=(n/rho)**(1.0/3.0)
  i=0
  do ix=0,m-1
    do iy=0,m-1
      do iz=0,m-1
        if(modulo(ix+iy+iz,2)==0) then
          i=i+1
          x(i)=c*((ix+0.5)/m-0.5)
          y(i)=c*((iy+0.5)/m-0.5)
          z(i)=c*((iz+0.5)/m-0.5)
        end if
      end do
    end do
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

end program lmclj

!**********************************************************************!
