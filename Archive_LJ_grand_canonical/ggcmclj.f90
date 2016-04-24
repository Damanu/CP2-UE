!**********************************************************************!
!
! File: ggcmclj.f90
!
! Create random ("gas") initial configuration for Grand Canonical (TVmu)
! Monte Carlo of Lennard-Jonesium
!
! 04-Feb-2004 (MN)
! 18-Apr-2012
!
!**********************************************************************!

program ggcmclj

!**********************************************************************!

  implicit none

  integer,parameter::iout=2

  character(len=80)::fname
  integer::i,n,nran,nt,ntjob,ntprint,ntskip
  integer,dimension(:),allocatable::iran
  real::c,disp,dr,pcre,ran,t,v,zz
  real,dimension(:),allocatable::x,y,z

! User input

  write(unit=*,fmt="(a)",advance="no") "              t="
  read(unit=*,fmt=*) t
  write(unit=*,fmt="(a)",advance="no") "              v="
  read(unit=*,fmt=*) v
  write(unit=*,fmt="(a)",advance="no") "              z="
  read(unit=*,fmt=*) zz
  write(unit=*,fmt="(a)",advance="no") "              n="
  read(unit=*,fmt=*) n

  do
    write(unit=*,fmt="(a)",advance="no") "     p(cre/del)="
    read(unit=*,fmt=*) pcre
    if(pcre<0.5) then
      exit
    end if
  end do

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

  write(unit=*,fmt="(a)",advance="no") "          fname=[gcmclj_in.dat] "
  read(unit=*,fmt="(a)") fname
  if(fname=="") then
    fname="gcmclj_in.dat"
  end if

! RNG characteristics

  call random_seed(size=nran)

! Allocate arrays

  allocate(iran(nran),x(n),y(n),z(n))

! Box length

  c=v**(1.0/3.0)
  
! Random particle positions

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
  write(unit=iout) disp,dr,pcre,t,v,zz
  write(unit=iout) x,y,z

  close(unit=iout)

! Deallocate arrays

  deallocate(iran,x,y,z)

end program ggcmclj

!**********************************************************************!
