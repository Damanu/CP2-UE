!**********************************************************************!
!
! File: zgemclj.f90
!
! Gibbs Ensemble (NVT) Monte Carlo of Lennard-Jonesium
! Re-initialize startup/checkpoint file
!
! 15-May-1999 (MN)
! 19-Apr-2012
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

program zgemclj

!**********************************************************************!

  use getval_m

  implicit none

  integer,parameter::iin=1,iout=2

  character(len=80)::fname
  integer::n,n1,n2,nran,nt,ntjob,ntprint,ntskip
  integer,dimension(:),allocatable::iran
  real::disp1,disp2,dv,fact,ptr,ptrold,pvm,pvmold,t,v,v1,vold
  real,dimension(:),allocatable::x,y,z

! RNG characteristics

  call random_seed(size=nran)
  allocate(iran(nran))

! Read old checkpoint file

  fname="gemclj_out.dat"
  write(unit=*,fmt="(a)",advance="no") &
       "         infile=[gemclj_out.dat] "
  call getval(fname)

  open(unit=iin,file=fname,status="old",action="read", &
       form="unformatted",position="rewind")

  read(unit=iin) n,n1,nt,ntjob,ntprint,ntskip,iran
  read(unit=iin) disp1,disp2,dv,ptr,pvm,t,v,v1

! Allocate arrays

  n2=n-n1

  allocate(x(n),y(n),z(n))

! Positions

  read(unit=iin) x,y,z

  close(unit=iin)

! User input

  ptrold=ptr
  pvmold=pvm
  vold=v

  write(unit=*,fmt="(a,i15)") &
       "             n1=",n1
  write(unit=*,fmt="(a,i15)") &
       "             n2=",n2
  if(v<99999999.9) then
    write(unit=*,fmt="(a,f14.5,a)",advance="no") &
         "              v=[",v,"] "
  else
    write(unit=*,fmt="(a,es14.5,a)",advance="no") &
         "              v=[",v,"] "
  end if
  call getval(v)
  write(unit=*,fmt="(a,f14.5,a)",advance="no") &
       "              t=[",t,"] "
  call getval(t)
  write(unit=*,fmt="(a,f14.5,a)",advance="no") &
       "          disp1=[",disp1,"] "
  call getval(disp1)
  write(unit=*,fmt="(a,f14.5,a)",advance="no") &
       "          disp2=[",disp2,"] "
  call getval(disp2)
  write(unit=*,fmt="(a,f14.5,a)",advance="no") &
       "             dv=[",dv,"] "
  call getval(dv)

  do
    pvm=pvmold
    ptr=ptrold
    write(unit=*,fmt="(a,f14.5,a)",advance="no") &
         "    prob(vmove)=[",pvm,"] "
    call getval(pvm)
    write(unit=*,fmt="(a,f14.5,a)",advance="no") &
         "    prob(trans)=[",ptr,"] "
    call getval(ptr)
    if(ptr+pvm<=1.0) then
      exit
    end if
  end do

  write(unit=*,fmt="(a,i14,a)",advance="no") &
       "         ntskip=[",ntskip,"] "
  call getval(ntskip)
  write(unit=*,fmt="(a,i14,a)",advance="no") &
       " ntprint/ntskip=[",ntprint,"] "
  call getval(ntprint)
  write(unit=*,fmt="(a,i14,a)",advance="no") &
       "   ntjob/ntskip=[",ntjob,"] "
  call getval(ntjob)
  fname="gemclj_in.dat"
  write(unit=*,fmt="(a)",advance="no") &
       "        outfile=[ gemclj_in.dat] "
  call getval(fname)

! Rescale positions

  if(v/=vold) then
    v1=(v/vold)*v1
    fact=(v/vold)**(1.0/3.0)
    x=fact*x
    y=fact*y
    z=fact*z
  end if

! Write new startup file

  nt=0

  open(unit=iout,file=fname,status="replace",action="write", &
       form="unformatted")

  write(unit=iout) n,n1,nt,ntjob,ntprint,ntskip,iran
  write(unit=iout) disp1,disp2,dv,ptr,pvm,t,v,v1
  write(unit=iout) x,y,z

  close(unit=iout)

! Deallocate arrays

  deallocate(iran,x,y,z)

end program zgemclj

!**********************************************************************!
