!**********************************************************************!
!
! File: zgcmclj.f90
!
! Grand Canonical (TVmu) Monte Carlo of Lennard-Jonesium
! Re-initialize startup/checkpoint file
!
! 08-Feb-2004 (MN)
! 18-Apr-2012
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

program zgcmclj

!**********************************************************************!

  use getval_m

  implicit none

  integer,parameter::iin=1,iout=2

  character(len=80)::fname
  integer::n,nt,ntjob,ntprint,ntskip
  integer::nran
  integer,dimension(:),allocatable::iran
  real::disp,dr,pcre,t,v,zz
  real::fact,v_old
  real,dimension(:),allocatable::x,y,z

! RNG characteristics

  call random_seed(size=nran)

  allocate(iran(nran))

! Read old checkpoint file

  fname="gcmclj_out.dat"
  write(unit=*,fmt="(a)",advance="no") &
       "         infile=[gcmclj_out.dat] "
  call getval(fname)

  open(unit=iin,file=fname,status="old",action="read", &
       form="unformatted",position="rewind")

  read(unit=iin) n,nt,ntjob,ntprint,ntskip,iran
  read(unit=iin) disp,dr,pcre,t,v,zz

! Allocate arrays

  allocate(x(n),y(n),z(n))

! Positions

  read(unit=iin) x,y,z

  close(unit=iin)

! User input

  write(unit=*,fmt="(a,f14.7,a)",advance="no") &
       "              t=[",t,"] "
  call getval(t)
  v_old=v
  write(unit=*,fmt="(a,f14.7,a)",advance="no") &
       "              v=[",v,"] "
  call getval(v)
  write(unit=*,fmt="(a,i15)") "              n=",n
  write(unit=*,fmt="(a,f14.7,a)",advance="no") &
       "              z=[",zz,"] "
  call getval(zz)

  do
    write(unit=*,fmt="(a,f14.7,a)",advance="no") &
         "     p(cre/del)=[",pcre,"] "
    call getval(pcre)
    if(pcre<0.5) then
      exit
    end if
  end do

  write(unit=*,fmt="(a,f14.7,a)",advance="no") &
       "           disp=[",disp,"] "
  call getval(disp)
  write(unit=*,fmt="(a,f14.7,a)",advance="no") &
       "             dr=[",dr,"] "
  call getval(dr)
  write(unit=*,fmt="(a,i14,a)",advance="no") &
       "         ntskip=[",ntskip,"] "
  call getval(ntskip)
  write(unit=*,fmt="(a,i14,a)",advance="no") &
       " ntprint/ntskip=[",ntprint,"] "
  call getval(ntprint)
  write(unit=*,fmt="(a,i14,a)",advance="no") &
       "   ntjob/ntskip=[",ntjob,"] "
  call getval(ntjob)

  fname="gcmclj_in.dat"
  write(unit=*,fmt="(a)",advance="no") &
       "        outfile=[ gcmclj_in.dat] "
  call getval(fname)

! Rescale positions

  if(v/=v_old) then
    fact=(v/v_old)**(1.0/3.0)
    x=fact*x
    y=fact*y
    z=fact*z
  end if

! Write new startup file

  nt=0

  open(unit=iout,file=fname,status="replace",action="write", &
       form="unformatted")

  write(unit=iout) n,nt,ntjob,ntprint,ntskip,iran
  write(unit=iout) disp,dr,pcre,t,v,zz
  write(unit=iout) x,y,z

  close(unit=iout)

! Deallocate arrays

  deallocate(iran,x,y,z)

end program zgcmclj

!**********************************************************************!
