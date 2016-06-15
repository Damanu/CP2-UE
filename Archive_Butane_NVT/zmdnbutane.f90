!**********************************************************************!
!
! File: zmdnbutane.f90
!
! NVT (Nose--Hoover) Molcular Dynamics of OPLS n-butane
! Re-initialize startup/checkpoint file
!
! 19-Apr-2009 (MN)
! 31-May-2012
!
!**********************************************************************!

module getval_m

!**********************************************************************!

! Generic subroutine: read item from standard input/keep old value on
! hitting return

  implicit none

  integer,parameter,public::double=selected_real_kind(15)

  private::getdouble,getint,getreal,getstrg
  public::getval

  interface getval
    module procedure getdouble,getint,getreal,getstrg
  end interface

  character(len=80),private::line

contains

!**********************************************************************!

  subroutine getdouble(doubleval)

!**********************************************************************!

    real(kind=double),intent(in out)::doubleval

    read(unit=*,fmt="(a)") line
    if(len_trim(adjustl(line))>0) then
      read(unit=line,fmt=*) doubleval
    end if

    return

  end subroutine getdouble

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

program zmdnbutane

!**********************************************************************!

  use getval_m

  implicit none

  integer,parameter::iin=1,iout=2
  integer,parameter::m=4,nf=7
  real(kind=double),parameter::avog=602.252_double,rg=0.0083143_double
  real(kind=double),parameter::xmch2=14.0_double,xmch3=15.0_double
  real(kind=double),dimension(m),parameter:: &
       xms=(/xmch3,xmch2,xmch2,xmch3/)

  character(len=80)::fname
  integer::i
  integer::n,nt,ntjob,ntprint,ntskip
  real(kind=double)::fact,sk,vmol,vmol_old,xc,yc,zc
  real(kind=double)::dphi,dr,dt,t,vol,tau,zeta
  real(kind=double),dimension(:,:),allocatable::x,y,z,xl,yl,zl

! Read old checkpoint file

  fname="mdnbutane_out.dat"
  write(unit=*,fmt="(a)",advance="no") &
       "         infile=[mdnbutane_out.dat] "
  call getval(fname)

  open(unit=iin,file=fname,status="old",action="read", &
       form="unformatted",position="rewind")

  read(unit=iin) n,nt,ntjob,ntprint,ntskip
  read(unit=iin) dphi,dr,dt,t,vol,tau,zeta

! Allocate arrays

  allocate(x(m,n),y(m,n),z(m,n),xl(m,n),yl(m,n),zl(m,n))

! Positions

  read(unit=iin) x,y,z
  read(unit=iin) xl,yl,zl

  close(unit=iin)

! User input

  write(unit=*,fmt="(a,i18)") &
       "              n=",n
  write(unit=*,fmt="(a,f17.7,a)",advance="no") &
       "           t(K)=[",t,"] "
  call getval(t)

  vmol=vol*avog/n
  vmol_old=vmol
  write(unit=*,fmt="(a,f17.7,a)",advance="no") &
       " vol(cm**3/mol)=[",vmol,"] "
  call getval(vmol)

  write(unit=*,fmt="(a,f17.7,a)",advance="no") &
       "        tau(ps)=[",tau,"] "
  call getval(tau)
  write(unit=*,fmt="(a,f17.7,a)",advance="no") &
       "         dt(ps)=[",dt,"] "
  call getval(dt)
  write(unit=*,fmt="(a,f17.7,a)",advance="no") &
       "             dr=[",dr,"] "
  call getval(dr)
  write(unit=*,fmt="(a,f17.7,a)",advance="no") &
       "      dphi(deg)=[",dphi,"] "
  call getval(dphi)
  write(unit=*,fmt="(a,i17,a)",advance="no") &
       "         ntskip=[",ntskip,"] "
  call getval(ntskip)
  write(unit=*,fmt="(a,i17,a)",advance="no") &
       " ntprint/ntskip=[",ntprint,"] "
  call getval(ntprint)
  write(unit=*,fmt="(a,i17,a)",advance="no") &
       "   ntjob/ntskip=[",ntjob,"] "
  call getval(ntjob)

  fname="mdnbutane_in.dat"
  write(unit=*,fmt="(a)",advance="no") &
       "        outfile=[ mdnbutane_in.dat] "
  call getval(fname)

! Rescale (center) positions

  if(vmol/=vmol_old) then

    fact=(vmol/vmol_old)**(1.0_double/3.0_double)
    vol=vmol*n/avog
    
    do i=1,n
      xc=0.5*double*(x(2,i)+x(3,i))
      yc=0.5*double*(y(2,i)+y(3,i))
      zc=0.5*double*(z(2,i)+z(3,i))
      x(:,i)=fact*xc+(x(:,i)-xc)
      y(:,i)=fact*yc+(y(:,i)-yc)
      z(:,i)=fact*zc+(z(:,i)-zc)
      xc=0.5*double*(xl(2,i)+xl(3,i))
      yc=0.5*double*(yl(2,i)+yl(3,i))
      zc=0.5*double*(zl(2,i)+zl(3,i))
      xl(:,i)=fact*xc+(xl(:,i)-xc)
      yl(:,i)=fact*yc+(yl(:,i)-yc)
      zl(:,i)=fact*zc+(zl(:,i)-zc)
    end do

  end if

! Subtract net momentum & rescale kinetic energy

  xl=xl+sum(spread(xms,2,n)*(x-xl))/(n*sum(xms))
  yl=yl+sum(spread(xms,2,n)*(y-yl))/(n*sum(xms))
  zl=zl+sum(spread(xms,2,n)*(z-zl))/(n*sum(xms))

  sk=sum(spread(xms,2,n)*((x-xl)**2+(y-yl)**2+(z-zl)**2))/dt**2

  fact=sqrt(nf*n*rg*t/sk)

  xl=x-fact*(x-xl)
  yl=y-fact*(y-yl)
  zl=z-fact*(z-zl)

! Write new checkpoint file

  nt=0
  zeta=0.0_double

  open(unit=iout,file=fname,status="replace",action="write", &
       form="unformatted")

  write(unit=iout) n,nt,ntjob,ntprint,ntskip
  write(unit=iout) dphi,dr,dt,t,vol,tau,zeta
  write(unit=iout) x,y,z
  write(unit=iout) xl,yl,zl

  close(unit=iout)

! Deallocate arrays

  deallocate(x,y,z,xl,yl,zl)

end program zmdnbutane

!**********************************************************************!
