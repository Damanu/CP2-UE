!**********************************************************************!
!
! File: zmdlj.f90
!
! NVE Molecular Dynamics of Lennard-Jonesium
! Re-initialize checkpoint file
!
! 09-May-1997 (MN)
! 04-May-2012
!
!**********************************************************************!

module getval_m

!**********************************************************************!

! Generic subroutine: read item from standard input/keep old value on
! hitting return

  implicit none

  private::getint,getreal,getstrg
  public::getval

  interface getval
    module procedure getint,getreal,getstrg
  end interface

  character(len=80),private::line

contains

!**********************************************************************!

  subroutine getint(intval)

!**********************************************************************!

    integer,intent(in out)::intval

    read(unit=*,fmt="(a)") line
    if(len_trim(adjustl(line))>0) then
      read(unit=line,fmt=*) intval
    end if

    return

  end subroutine getint

!**********************************************************************!

  subroutine getreal(realval)

!**********************************************************************!

    real,intent(in out)::realval

    read(unit=*,fmt="(a)") line
    if(len_trim(adjustl(line))>0) then
      read(unit=line,fmt=*) realval
    end if

    return

  end subroutine getreal

!**********************************************************************!

  subroutine getstrg(strgval)

!**********************************************************************!

    character(len=*),intent(in out)::strgval

    read(unit=*,fmt="(a)") line
    if(len_trim(adjustl(line))>0) then
      read(unit=line,fmt="(a)") strgval
    end if

    return

  end subroutine getstrg

end module getval_m

!**********************************************************************!

program zmdlj

!**********************************************************************!

  use getval_m

  implicit none

  integer,parameter::iin=1,iout=2

  character(len=80)::fname
  integer::n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint
  real::dr,dt,fact,rho,rhon,sk,skn
  real,dimension(:),allocatable::vx,vy,vz,x,y,z

! Read checkpoint file

  fname="mdlj_out.dat"
  write(unit=*,fmt="()")
  write(unit=*,fmt="(a)",advance="no") &
       "          infile=[mdlj_out.dat] "
  call getval(fname)

  open(unit=iin,file=fname,status="old",action="read", &
       form="unformatted",position="rewind")

! Simulation parameters

  read(unit=iin) n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint
  read(unit=iin) dr,dt,rho

! Allocate arrays

  allocate(vx(n),vy(n),vz(n),x(n),y(n),z(n))

! Positions & velocities

  read(unit=iin) x,y,z
  read(unit=iin) vx,vy,vz
  close(unit=iin)

! Save old density & kinetic energy

  rhon=rho
  sk=0.5*sum(vx**2+vy**2+vz**2)/n
  skn=sk

! User input

  write(unit=*,fmt="(a,i13)") "               n=",n
  write(unit=*,fmt="(a,f12.5,a)",advance="no") &
       "             rho=[",rho,"] "
  call getval(rhon)
  write(unit=*,fmt="(a,f12.5,a)",advance="no") &
       "             k/n=[",sk,"] "
  call getval(skn)
  write(unit=*,fmt="(a,f12.5,a)",advance="no") &
       "              dt=[",dt,"] "
  call getval(dt)
  write(unit=*,fmt="(a,f12.5,a)",advance="no") &
       "              dr=[",dr,"] "
  call getval(dr)
  write(unit=*,fmt="(a,i12,a)",advance="no") &
       "         ntaskip=[",ntaskip,"] "
  call getval(ntaskip)
  write(unit=*,fmt="(a,i12,a)",advance="no") &
       " ntprint/ntaskip=[",ntprint,"] "
  call getval(ntprint)
  write(unit=*,fmt="(a,i12,a)",advance="no") &
       "         ntcskip=[",ntcskip,"] "
  call getval(ntcskip)
  write(unit=*,fmt="(a,i12,a)",advance="no") &
       "    ncor/ntcskip=[",ncor,"] "
  call getval(ncor)
  write(unit=*,fmt="(a,i12,a)",advance="no") &
       "  ntorig/ntcskip=[",ntorig,"] "
  call getval(ntorig)
  write(unit=*,fmt="(a,i12,a)",advance="no") &
       "           ntjob=[",ntjob,"] "
  call getval(ntjob)

! Rescale positions

  if(rhon/=rho) then
    fact=(rho/rhon)**(1.0/3.0)
    rho=rhon
    x=fact*x
    y=fact*y
    z=fact*z
  end if

! Zero net momentum & rescale velocities

  vx=vx-sum(vx)/n
  vy=vy-sum(vy)/n
  vz=vz-sum(vz)/n
  if(skn/=sk) then
    fact=sqrt(skn/sk)
    vx=fact*vx
    vy=fact*vy
    vz=fact*vz
  end if

! Write new startup file

  fname="mdlj_in.dat"
  write(unit=*,fmt="(a)",advance="no") &
       "         outfile=[ mdlj_in.dat] "
  call getval(fname)

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

end program zmdlj

!**********************************************************************!
