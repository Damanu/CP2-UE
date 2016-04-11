!**********************************************************************!
!
! File: zmclj.f90
!
! NVT-Monte Carlo of Lennard-Jonesium
! Re-initialize checkpoint file
!
! 13-Apr-1997 (MN)
! 17-Apr-2012
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

program zmclj

!**********************************************************************!

  use getval_m

  implicit none

  integer,parameter::iin=1,iout=2

  character(len=80)::fname
  integer::n,nran,nt,ntjob,ntprint,ntskip
  integer,dimension(:),allocatable::iran
  real::disp,dr,fact,rho,rhon,t
  real,dimension(:),allocatable::x,y,z

! RNG characteristics

  call random_seed(size=nran)

  allocate(iran(nran))

! Read old checkpoint file

  fname="mclj_out.dat"
  write(unit=*,fmt="(a)",advance="no") "         infile=[mclj_out.dat] "
  call getval(fname)

  open(unit=iin,file=fname,status="old",action="read", &
       form="unformatted",position="rewind")

  read(unit=iin) n,nt,ntjob,ntprint,ntskip,iran
  read(unit=iin) disp,dr,rho,t

! Allocate arrays

  allocate(x(n),y(n),z(n))

! Positions

  read(unit=iin) x,y,z

  close(unit=iin)

  rhon=rho

! User input

  write(unit=*,fmt="(a,i13)") "              n=",n
  write(unit=*,fmt="(a,f12.5,a)",advance="no") "            rho=[",rhon,"] "
  call getval(rhon)
  write(unit=*,fmt="(a,f12.5,a)",advance="no") "              t=[",t,"] "
  call getval(t)
  write(unit=*,fmt="(a,f12.5,a)",advance="no") "           disp=[",disp,"] "
  call getval(disp)
  write(unit=*,fmt="(a,f12.5,a)",advance="no") "             dr=[",dr,"] "
  call getval(dr)
  write(unit=*,fmt="(a,i12,a)",advance="no") "         ntskip=[",ntskip,"] "
  call getval(ntskip)
  write(unit=*,fmt="(a,i12,a)",advance="no") " ntprint/ntskip=[",ntprint,"] "
  call getval(ntprint)
  write(unit=*,fmt="(a,i12,a)",advance="no") "   ntjob/ntskip=[",ntjob,"] "
  call getval(ntjob)

! Rescale positions

  if(rhon/=rho) then
    fact=(rho/rhon)**(1.0/3.0)
    rho=rhon
    x=fact*x
    y=fact*y
    z=fact*z
  end if

! Write new startup file

  fname="mclj_in.dat"
  write(unit=*,fmt="(a)",advance="no") "        outfile=[ mclj_in.dat] "
  call getval(fname)

  nt=0

  open(unit=iout,file=fname,status="replace",action="write", &
       form="unformatted")

  write(unit=iout) n,nt,ntjob,ntprint,ntskip,iran
  write(unit=iout) disp,dr,rho,t
  write(unit=iout) x,y,z

  close(unit=iout)

! Deallocate arrays

  deallocate(iran,x,y,z)

end program zmclj

!**********************************************************************!
