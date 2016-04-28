!**********************************************************************!
!
! File: gmdlj.f90
!
! Create random ("gas") initial configuration for NVE Molecular
! Dynamics of Lennard-Jonesium
!
! 10-May-1997 (MN)
! 04-May-2012
!
!**********************************************************************!

program gmdlj

!**********************************************************************!

  implicit none

  integer,parameter::iout=2

  character(len=80)::fname
  integer::i,m,n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint
  real::c,dr,dt,ran,rho
  real,dimension(:),allocatable::vx,vy,vz,x,y,z

! User input

  magic: do
    write(unit=*,fmt="(a)",advance="no") "               n="
    read(unit=*,fmt=*) n
    m=0
    do
      m=m+2
      if(m**3==2*n) then       ! Check if magic number
        exit magic
      else if(m**3>2*n) then
        exit
      end if
    end do
  end do magic

  write(unit=*,fmt="(a)",advance="no") "             rho="
  read(unit=*,fmt=*) rho
  write(unit=*,fmt="(a)",advance="no") "              dt="
  read(unit=*,fmt=*) dt
  write(unit=*,fmt="(a)",advance="no") "              dr="
  read(unit=*,fmt=*) dr
  write(unit=*,fmt="(a)",advance="no") "         ntaskip="
  read(unit=*,fmt=*) ntaskip
  write(unit=*,fmt="(a)",advance="no") " ntprint/ntaskip="
  read(unit=*,fmt=*) ntprint
  write(unit=*,fmt="(a)",advance="no") "         ntcskip="
  read(unit=*,fmt=*) ntcskip
  write(unit=*,fmt="(a)",advance="no") "    ncor/ntcskip="
  read(unit=*,fmt=*) ncor
  write(unit=*,fmt="(a)",advance="no") "  ntorig/ntcskip="
  read(unit=*,fmt=*) ntorig
  write(unit=*,fmt="(a)",advance="no") "           ntjob="
  read(unit=*,fmt=*) ntjob

  write(unit=*,fmt="(a)",advance="no") "           fname=[mdlj_in.dat] "
  read(unit=*,fmt="(a)") fname
  if(fname=="") then
    fname="mdlj_in.dat"
  end if

! Allocate arrays

  allocate(vx(n),vy(n),vz(n),x(n),y(n),z(n))

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

! Zero velocities

  vx=0.0
  vy=0.0
  vz=0.0

! Write startup file

  nt=0

  open(unit=iout,file=fname,status="replace",action="write", &
       form="unformatted")

  write(unit=iout) n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint
  write(unit=iout) dr,dt,rho
  write(unit=iout) x,y,z
  write(unit=iout) vx,vy,vz

  close(unit=iout)

! Deallocate arrays

  deallocate(vx,vy,vz,x,y,z)

end program gmdlj

!**********************************************************************!
