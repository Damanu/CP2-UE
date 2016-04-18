!**********************************************************************!
!
! File: znptmclj.f90
!
! NpT-Monte Carlo of Lennard-Jonesium
! Re-initialize checkpoint file
!
! 01-May-1999 (MN)
! 17-Apr-2010
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

program znptmclj

!**********************************************************************!

  use getval_m

  implicit none

  integer,parameter::iin=1,iout=2

  character(len=80)::fname
  integer::n,nt,ntjob,ntprint,ntskip
  integer::nran
  integer,dimension(:),allocatable::iran
  real::disp,dlnv,dr,p,pvm,t,v
  real,dimension(:),allocatable::x,y,z

! RNG characteristics

  call random_seed(size=nran)

  allocate(iran(nran))

! Read old checkpoint file

  fname="nptmclj_out.dat"
  write(unit=*,fmt="(a)",advance="no") &
       "         infile=[nptmclj_out.dat] "
  call getval(fname)

  open(unit=iin,file=fname,status="old",action="read", &
       form="unformatted",position="rewind")

  read(unit=iin) n,nt,ntjob,ntprint,ntskip,iran
  read(unit=iin) disp,dlnv,dr,p,pvm,t,v

! Allocate arrays

  allocate(x(n),y(n),z(n))

! Positions

  read(unit=iin) x,y,z

  close(unit=iin)

! User input

  write(unit=*,fmt="(a,i16)") &
       "              n=",n
  write(unit=*,fmt="(a,f15.5,a)",advance="no") &
       "              p=[",p,"] "
  call getval(p)
  write(unit=*,fmt="(a,f15.5,a)",advance="no") &
       "              t=[",t,"] "
  call getval(t)
  write(unit=*,fmt="(a,f15.5,a)",advance="no") &
       "           disp=[",disp,"] "
  call getval(disp)
  write(unit=*,fmt="(a,f15.5,a)",advance="no") &
       "           dv/v=[",dlnv,"] "
  call getval(dlnv)
  write(unit=*,fmt="(a,f15.5,a)",advance="no") &
       "    prob(vmove)=[",pvm,"] "
  call getval(pvm)
  write(unit=*,fmt="(a,f15.5,a)",advance="no") &
       "             dr=[",dr,"] "
  call getval(dr)
  write(unit=*,fmt="(a,i15,a)",advance="no") &
       "         ntskip=[",ntskip,"] "
  call getval(ntskip)
  write(unit=*,fmt="(a,i15,a)",advance="no") &
       " ntprint/ntskip=[",ntprint,"] "
  call getval(ntprint)
  write(unit=*,fmt="(a,i15,a)",advance="no") &
       "   ntjob/ntskip=[",ntjob,"] "
  call getval(ntjob)

! Write new startup file

  fname="nptmclj_in.dat"
  write(unit=*,fmt="(a)",advance="no") &
       "        outfile=[ nptmclj_in.dat] "
  call getval(fname)

  nt=0

  open(unit=iout,file=fname,status="replace",action="write", &
       form="unformatted")

  write(unit=iout) n,nt,ntjob,ntprint,ntskip,iran
  write(unit=iout) disp,dlnv,dr,p,pvm,t,v
  write(unit=iout) x,y,z

  close(unit=iout)

! Deallocate arrays

  deallocate(iran,x,y,z)

end program znptmclj

!**********************************************************************!
