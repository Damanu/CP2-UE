!**********************************************************************!
!
! File: lisomdlj.f90
!
! Create initial configuration (fcc lattice) for Gaussian isokinetic
! (NVT) Molecular Dynamics of Lennard-Jonesium
! Random velocities
!
! 03-Jun-1999 (MN)
! 05-May-2012
!
!**********************************************************************!

program lisomdlj

!**********************************************************************!

  implicit none

  integer,parameter::iout=2

  character(len=80)::fname
  integer::i,ix,iy,iz,m
  integer::n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint
  real::c,fact,ran
  real::dr,dt,rho,t
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
  write(unit=*,fmt="(a)",advance="no") "               t="
  read(unit=*,fmt=*) t
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

  write(unit=*,fmt="(a)",advance="no") &
       "           fname=[isomdlj_in.dat] "
  read(unit=*,fmt="(a)") fname
  if(fname=="") then
    fname="isomdlj_in.dat"
  end if

! Allocate arrays

  allocate(vx(n),vy(n),vz(n),x(n),y(n),z(n))

! Fcc lattice positions & random velocities

  c=(n/rho)**(1.0/3.0)

  i=0
  do ix=0,m-1
    do iy=0,m-1
      do iz=0,m-1
        if(modulo(ix+iy+iz,2)==0) then
          i=i+1
          x(i)=((ix+0.5)/m-0.5)*c
          y(i)=((iy+0.5)/m-0.5)*c
          z(i)=((iz+0.5)/m-0.5)*c
          call random_number(ran)
          vx(i)=ran-0.5
          call random_number(ran)
          vy(i)=ran-0.5
          call random_number(ran)
          vz(i)=ran-0.5
        end if
      end do
    end do
  end do

! Zero net momentum & set initial temperature

  vx=vx-sum(vx)/n
  vy=vy-sum(vy)/n
  vz=vz-sum(vz)/n

  fact=sqrt(3.0*n*t/sum(vx**2+vy**2+vz**2))
  vx=fact*vx
  vy=fact*vy
  vz=fact*vz

! Write startup file

  nt=0

  open(unit=iout,file=fname,status="replace",action="write", &
       form="unformatted")

  write(unit=iout) n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint
  write(unit=iout) dr,dt,rho,t
  write(unit=iout) x,y,z
  write(unit=iout) vx,vy,vz

  close(unit=iout)

! Deallocate arrays

  deallocate(vx,vy,vz,x,y,z)

end program lisomdlj

!**********************************************************************!
