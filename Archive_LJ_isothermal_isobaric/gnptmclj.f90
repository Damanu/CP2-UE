!**********************************************************************!
!
! File: gnptmclj.f90
!
! Create random ("gas") initial configuration for NpT-Monte Carlo of
! Lennard-Jonesium
!
! 29-Apr-1999 (MN)
! 17-Apr-2012
!
!**********************************************************************!

program gnptmclj

!**********************************************************************!

  implicit none

  integer,parameter::iout=2

  character(len=80)::fname
  integer::n,nt,ntjob,ntprint,ntskip
  integer::i,m,nran
  integer,dimension(:),allocatable::iran
  real::disp,dlnv,dr,p,pvm,t,v
  real::c,ran,rho
  real,dimension(:),allocatable::x,y,z

! User input

  magic: do
    write(unit=*,fmt="(a)",advance="no") "              n="
    read(unit=*,fmt=*) n
    m=0
    do
      m=m+2
      if(m**3==2*n) then     ! Check if magic number
        exit magic
      else if(m**3>2*n) then
        exit
      end if
    end do
  end do magic

  write(unit=*,fmt="(a)",advance="no") "              p="
  read(unit=*,fmt=*) p
  write(unit=*,fmt="(a)",advance="no") "              t="
  read(unit=*,fmt=*) t
  write(unit=*,fmt="(a)",advance="no") "            rho="
  read(unit=*,fmt=*) rho
  write(unit=*,fmt="(a)",advance="no") "           disp="
  read(unit=*,fmt=*) disp
  write(unit=*,fmt="(a)",advance="no") "           dv/v="
  read(unit=*,fmt=*) dlnv
  write(unit=*,fmt="(a)",advance="no") "    prob(vmove)="
  read(unit=*,fmt=*) pvm
  write(unit=*,fmt="(a)",advance="no") "             dr="
  read(unit=*,fmt=*) dr
  write(unit=*,fmt="(a)",advance="no") "         ntskip="
  read(unit=*,fmt=*) ntskip
  write(unit=*,fmt="(a)",advance="no") " ntprint/ntskip="
  read(unit=*,fmt=*) ntprint
  write(unit=*,fmt="(a)",advance="no") "   ntjob/ntskip="
  read(unit=*,fmt=*) ntjob
  write(unit=*,fmt="(a)",advance="no") &
       "          fname=[nptmclj_in.dat] "
  read(unit=*,fmt="(a)") fname
  if(fname=="") then
    fname="nptmclj_in.dat"
  end if

! RNG characteristics

  call random_seed(size=nran)

! Allocate arrays

  allocate(iran(nran),x(n),y(n),z(n))

! Random positions

  v=n/rho
  c=v**(1.0/3.0)
  do i=1,n
    call random_number(ran)
    x(i)=(ran-0.5)*c
    call random_number(ran)
    y(i)=(ran-0.5)*c
    call random_number(ran)
    z(i)=(ran-0.5)*c
  end do

! RNG seed

  call random_seed(get=iran)

! Write startup file

  nt=0

  open(unit=iout,file=fname,status="replace",action="write", &
       form="unformatted")

  write(unit=iout) n,nt,ntjob,ntprint,ntskip,iran
  write(unit=iout) disp,dlnv,dr,p,pvm,t,v
  write(unit=iout) x,y,z

  close(unit=iout)

! Deallocate arrays

  deallocate(iran,x,y,z)

end program gnptmclj

!**********************************************************************!
