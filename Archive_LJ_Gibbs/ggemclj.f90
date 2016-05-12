!**********************************************************************!
!
! File: ggemclj.f90
!
! Create random ("gas") initial configuration for Gibbs Ensemble (NVT)
! Monte Carlo of Lennard-Jonesium
!
! 13-May-1999 (MN)
! 19-Apr-2012
!
!**********************************************************************!

program ggemclj

!**********************************************************************!

  implicit none

  integer,parameter::iout=2

  character(len=80)::fname
  integer::i,n,n1,n2,nran,nt,ntjob,ntprint,ntskip
  integer,dimension(:),allocatable::iran
  real::c1,c2,disp1,disp2,dv,ptr,pvm,ran,t,v,v1,v2
  real,dimension(:),allocatable::x,y,z

! User input

  write(unit=*,fmt="(a)",advance="no") "             n1="
  read(unit=*,fmt=*) n1
  write(unit=*,fmt="(a)",advance="no") "             n2="
  read(unit=*,fmt=*) n2
  write(unit=*,fmt="(a)",advance="no") "             v1="
  read(unit=*,fmt=*) v1
  write(unit=*,fmt="(a)",advance="no") "             v2="
  read(unit=*,fmt=*) v2
  write(unit=*,fmt="(a)",advance="no") "              t="
  read(unit=*,fmt=*) t
  write(unit=*,fmt="(a)",advance="no") "          disp1="
  read(unit=*,fmt=*) disp1
  write(unit=*,fmt="(a)",advance="no") "          disp2="
  read(unit=*,fmt=*) disp2
  write(unit=*,fmt="(a)",advance="no") "             dv="
  read(unit=*,fmt=*) dv

  do
    write(unit=*,fmt="(a)",advance="no") "    prob(vmove)="
    read(unit=*,fmt=*) pvm
    write(unit=*,fmt="(a)",advance="no") "    prob(trans)="
    read(unit=*,fmt=*) ptr
    if(ptr+pvm<=1.0) then
      exit
    end if
  end do

  write(unit=*,fmt="(a)",advance="no") "         ntskip="
  read(unit=*,fmt=*) ntskip
  write(unit=*,fmt="(a)",advance="no") " ntprint/ntskip="
  read(unit=*,fmt=*) ntprint
  write(unit=*,fmt="(a)",advance="no") "   ntjob/ntskip="
  read(unit=*,fmt=*) ntjob
  write(unit=*,fmt="(a)",advance="no") &
       "          fname=[gemclj_in.dat] "
  read(unit=*,fmt="(a)") fname
  if(fname=="") then
    fname="gemclj_in.dat"
  end if

! Parameters

  n=n1+n2
  c1=v1**(1.0/3.0)
  c2=v2**(1.0/3.0)
  v=v1+v2

! RNG characteristics

  call random_seed(size=nran)

! Allocate arrays

  allocate(iran(nran),x(n),y(n),z(n))

! Random positions: Box 1

  do i=1,n1
    call random_number(ran)
    x(i)=(ran-0.5)*c1
    call random_number(ran)
    y(i)=(ran-0.5)*c1
    call random_number(ran)
    z(i)=(ran-0.5)*c1
  end do

! Box 2

  do i=n1+1,n
    call random_number(ran)
    x(i)=(ran-0.5)*c2
    call random_number(ran)
    y(i)=(ran-0.5)*c2
    call random_number(ran)
    z(i)=(ran-0.5)*c2
  end do

! RNG seed

  call random_seed(get=iran)

! Write startup file

  nt=0

  open(unit=iout,file=fname,status="replace",action="write", &
       form="unformatted")

  write(unit=iout) n,n1,nt,ntjob,ntprint,ntskip,iran
  write(unit=iout) disp1,disp2,dv,ptr,pvm,t,v,v1
  write(unit=iout) x,y,z

  close(unit=iout)

! Deallocate arrays

  deallocate(iran,x,y,z)

end program ggemclj

!**********************************************************************!
