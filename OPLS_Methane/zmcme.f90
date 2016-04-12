!**********************************************************************!
!
! File: zmcme.f90
!
! NVT-Monte Carlo of OPLS-AA Methane
! Re-initialize checkpoint file
!
! 23-Apr-1999 (MN)
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

program zmcme

!**********************************************************************!

  use getval_m

  implicit none

  integer,parameter::iin=1,iout=2,m=5
  real,parameter::xmc=12.0,xmh=1.0
  real,dimension(m),parameter::xmm=(/xmc,xmh,xmh,xmh,xmh/)

  character(len=80)::fname
  integer::i,n,nran,nt,ntjob,ntprint,ntskip
  integer,dimension(:),allocatable::iran
  real::da,disp,dr,fact,rho,rho_new,t,xcm,ycm,zcm
  real,dimension(:,:),allocatable::x,y,z

! RNG characteristics

  call random_seed(size=nran)

  allocate(iran(nran))

! Read old checkpoint file

  fname="mcme_out.dat"
  write(unit=*,fmt="(a)",advance="no") &
       "         infile=[mcme_out.dat] "
  call getval(fname)

  open(unit=iin,file=fname,status="old",action="read", &
       form="unformatted",position="rewind")

  read(unit=iin) n,nt,ntjob,ntprint,ntskip,iran
  read(unit=iin) da,disp,dr,rho,t

! Allocate arrays

  allocate(x(m,n),y(m,n),z(m,n))

! Positions

  read(unit=iin) x,y,z

  close(unit=iin)

  rho_new=rho

! User input

  write(unit=*,fmt="(a,i13)") "              n=",n
  write(unit=*,fmt="(a,f12.5,a)",advance="no") &
       "   rho(1/nm**3)=[",rho_new,"] "
  call getval(rho_new)
  write(unit=*,fmt="(a,f12.5,a)",advance="no") &
       "           t(K)=[",t,"] "
  call getval(t)
  write(unit=*,fmt="(a,f12.5,a)",advance="no") &
       "       disp(nm)=[",disp,"] "
  call getval(disp)
  write(unit=*,fmt="(a,f12.5,a)",advance="no") &
       "        da(deg)=[",da,"] "
  call getval(da)
  write(unit=*,fmt="(a,f12.5,a)",advance="no") &
       "         dr(nm)=[",dr,"] "
  call getval(dr)
  write(unit=*,fmt="(a,i12,a)",advance="no") &
       "         ntskip=[",ntskip,"] "
  call getval(ntskip)
  write(unit=*,fmt="(a,i12,a)",advance="no") &
       " ntprint/ntskip=[",ntprint,"] "
  call getval(ntprint)
  write(unit=*,fmt="(a,i12,a)",advance="no") &
       "   ntjob/ntskip=[",ntjob,"] "
  call getval(ntjob)

! Rescale positions

  if(rho_new/=rho) then

    fact=(rho/rho_new)**(1.0/3.0)
    rho=rho_new
    do i=1,n
      xcm=sum(xmm*x(:,i))/sum(xmm)
      ycm=sum(xmm*y(:,i))/sum(xmm)
      zcm=sum(xmm*z(:,i))/sum(xmm)
      x(:,i)=fact*xcm+(x(:,i)-xcm)
      y(:,i)=fact*ycm+(y(:,i)-ycm)
      z(:,i)=fact*zcm+(z(:,i)-zcm)
    end do

  end if

! Write new startup file

  fname="mcme_in.dat"
  write(unit=*,fmt="(a)",advance="no") &
       "        outfile=[ mcme_in.dat] "
  call getval(fname)

  nt=0

  open(unit=iout,file=fname,status="replace",action="write", &
       form="unformatted")

  write(unit=iout) n,nt,ntjob,ntprint,ntskip,iran
  write(unit=iout) da,disp,dr,rho,t
  write(unit=iout) x,y,z

  close(unit=iout)

! Deallocate arrays

  deallocate(iran,x,y,z)

end program zmcme

!**********************************************************************!
