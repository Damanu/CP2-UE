!**********************************************************************!
!
! File: gcmclj.f90
!
! Grand Canonical (TVmu) Monte Carlo of Lennard-Jonesium
!
! 04-Feb-2004 (MN)
! 19-Apr-2012
!
!**********************************************************************!

module gcmclj_glbm

!**********************************************************************!

! Parameters & global variables

  implicit none

  private

  integer,parameter,public::long=selected_int_kind(18)
  integer,parameter,public::double=selected_real_kind(15)
  integer,parameter,public::iin=1,iout=2

  integer,public::n,nt,ntjob,ntprint,ntskip
  integer,public::naccc,naccd,naccm,ndr,nmax
  integer,dimension(:),allocatable,public::iran
  real,public::disp,dr,pcre,t,v,zz
  real,public::alnmax,c,c2,c2m1,drm1,pi,r2max,rc2,su0,sw0,zv
  real,dimension(:),allocatable,public::x,y,z
  integer(kind=long),dimension(:),allocatable,public::ag
  real(kind=double),public::accrc,accrd,accrm,an,an2,au,aw

end module gcmclj_glbm

!**********************************************************************!

module gcmclj_subm

!**********************************************************************!

  use gcmclj_glbm

  implicit none

  private

  public::create,delete,getcf,means,move,putcf

contains

!**********************************************************************!

  subroutine getcf()

!**********************************************************************!

    integer::nran
    real::rc

! Maximum argument for exponential & pi

    alnmax=log(0.1*huge(1.0))
    pi=4.0*atan(1.0)

! RNG characteristics

    call random_seed(size=nran)
    allocate(iran(nran))

! Read startup/checkpoint file

    open(unit=iin,file="gcmclj_in.dat",status="old",action="read", &
         form="unformatted",position="rewind")

    read(unit=iin) n,nt,ntjob,ntprint,ntskip,iran
    read(unit=iin) disp,dr,pcre,t,v,zz

! Allocate arrays

    ndr=int(0.5*v**(1.0/3.0)/dr)

    allocate(ag(0:ndr-1),x(n),y(n),z(n))

! Positions

    read(unit=iin) x(1:n),y(1:n),z(1:n)

! Read accumulated averages/clear accumulators

    if(nt>0) then
      read(unit=iin) accrc,accrd,accrm,an,an2,au,aw
      read(unit=iin) ag(0:ndr-1)
    else
      accrc=0.0_double
      accrd=0.0_double
      accrm=0.0_double
      an=0.0_double
      an2=0.0_double
      au=0.0_double
      aw=0.0_double
      ag=0_long
    end if

    close(unit=iin)

! Parameters

    nmax=n

    c=v**(1.0/3.0)
    c2=0.5*c
    c2m1=2.0/c
    drm1=1.0/dr
    r2max=c2**2
    zv=zz*v

! Cutoff & tail corrections

    rc=0.5*c
    rc2=rc**2

    su0=2.0*pi/v*(4.0/(9.0*rc**9)-4.0/(3.0*rc**3))
    sw0=2.0*pi/v*(24.0/(3.0*rc**3)-48.0/(9.0*rc**9))

! Initialize RNG

    call random_seed(put=iran)

! Clear acceptance counters

    naccc=0
    naccd=0
    naccm=0

    return

  end subroutine getcf

!**********************************************************************!

  subroutine create()

!**********************************************************************!

    logical::accept
    integer::j
    real::dx,dy,dz,r2,rm6,ran,su,xn,yn,zn
    real,dimension(:),allocatable::xt,yt,zt

! Random position

    call random_number(ran)
    xn=(ran-0.5)*c
    call random_number(ran)
    yn=(ran-0.5)*c
    call random_number(ran)
    zn=(ran-0.5)*c

! Energy difference

    su=(2*n+1)*su0

    do j=1,n
      dx=x(j)-xn
      dx=dx-int(dx*c2m1)*c
      dy=y(j)-yn
      dy=dy-int(dy*c2m1)*c
      dz=z(j)-zn
      dz=dz-int(dz*c2m1)*c
      r2=dx**2+dy**2+dz**2
      if(r2<rc2) then
        rm6=1.0/r2**3
        su=su+(4.0*rm6-4.0)*rm6
      end if
    end do

! Acceptance test

    if(su/t>alnmax) then
      accept=.false.
    else
      if(zv*exp(-su/t)/(n+1)>=1.0) then
        accept=.true.
      else
        call random_number(ran)
        if(ran<=zv*exp(-su/t)/(n+1)) then
          accept=.true.
        else
          accept=.false.
        end if
      end if
    end if

! Accept move & append particle to arrays

    if(accept) then
      naccc=naccc+1

! Extend arrays

      if(n>=nmax) then
        allocate(xt(n),yt(n),zt(n))
        xt=x(1:n)
        yt=y(1:n)
        zt=z(1:n)
        deallocate(x,y,z)
        nmax=n+1
        allocate(x(nmax),y(nmax),z(nmax))
        x(1:n)=xt
        y(1:n)=yt
        z(1:n)=zt
        deallocate(xt,yt,zt)
      end if

      n=n+1
      x(n)=xn
      y(n)=yn
      z(n)=zn

    end if

    return

  end subroutine create

!**********************************************************************!

  subroutine move()

!**********************************************************************!

    logical::accept
    integer::i,j
    real::dx,dy,dz,r2,rm6,ran,su,xn,yn,zn

! Select random particle

    call random_number(ran)
    i=min(int(n*ran+1.0),n)

! Trial move (displacement)

    call random_number(ran)
    xn=x(i)+(ran-0.5)*disp
    xn=xn-int(xn*c2m1)*c
    call random_number(ran)
    yn=y(i)+(ran-0.5)*disp
    yn=yn-int(yn*c2m1)*c
    call random_number(ran)
    zn=z(i)+(ran-0.5)*disp
    zn=zn-int(zn*c2m1)*c

! Energy difference

    su=0.0

    do j=1,n
      if(j/=i) then

! New position

        dx=x(j)-xn
        dx=dx-int(dx*c2m1)*c
        dy=y(j)-yn
        dy=dy-int(dy*c2m1)*c
        dz=z(j)-zn
        dz=dz-int(dz*c2m1)*c
        r2=dx**2+dy**2+dz**2
        if(r2<rc2) then
          rm6=1.0/r2**3
          su=su+(4.0*rm6-4.0)*rm6
        end if

! Old position

        dx=x(j)-x(i)
        dx=dx-int(dx*c2m1)*c
        dy=y(j)-y(i)
        dy=dy-int(dy*c2m1)*c
        dz=z(j)-z(i)
        dz=dz-int(dz*c2m1)*c
        r2=dx**2+dy**2+dz**2
        if(r2<rc2) then
          rm6=1.0/r2**3
          su=su-(4.0*rm6-4.0)*rm6
        end if

      end if
    end do

! Acceptance test

    if(su<=0.0) then
      accept=.true.
    else
      if(su/t>alnmax) then
        accept=.false.
      else
        call random_number(ran)
        if(ran<=exp(-su/t)) then
          accept=.true.
        else
          accept=.false.
        end if
      end if
    end if

! Accept move

    if(accept) then
      naccm=naccm+1
      x(i)=xn
      y(i)=yn
      z(i)=zn
    end if

    return

  end subroutine move

!**********************************************************************!

  subroutine delete()

!**********************************************************************!

    logical::accept
    integer::i,j
    real::dx,dy,dz,r2,rm6,ran,su

! Empty box?

    if(n<1) then
      return
    end if

! Random particle

    call random_number(ran)
    i=min(int(n*ran+1.0),n)

! Energy difference (ok for n=1)

    su=-(2*n-1)*su0

    do j=1,n
      if(j/=i) then
        dx=x(j)-x(i)
        dx=dx-int(dx*c2m1)*c
        dy=y(j)-y(i)
        dy=dy-int(dy*c2m1)*c
        dz=z(j)-z(i)
        dz=dz-int(dz*c2m1)*c
        r2=dx**2+dy**2+dz**2
        if(r2<rc2) then
          rm6=1.0/r2**3
          su=su-(4.0*rm6-4.0)*rm6
        end if
      end if
    end do

! Acceptance test

    if(su/t>alnmax) then
      accept=.false.
    else
      if(n*exp(-su/t)/zv>=1.0) then
        accept=.true.
      else
        call random_number(ran)
        if(ran<=n*exp(-su/t)/zv) then
          accept=.true.
        else
          accept=.false.
        end if
      end if
    end if

! Accept move & eliminate particle from arrays

    if(accept) then
      naccd=naccd+1
      x(i)=x(n)
      y(i)=y(n)
      z(i)=z(n)
      n=n-1
    end if

    return

  end subroutine delete

!**********************************************************************!

  subroutine means()

!**********************************************************************!

    integer::i,ir,j
    real::dx,dy,dz,r2,rm6,su,sw

! Potential energy, virial & g(r)

    su=n**2*su0
    sw=n**2*sw0

    do i=1,n-1
      do j=i+1,n

        dx=x(j)-x(i)
        dx=dx-int(dx*c2m1)*c
        dy=y(j)-y(i)
        dy=dy-int(dy*c2m1)*c
        dz=z(j)-z(i)
        dz=dz-int(dz*c2m1)*c
        r2=dx**2+dy**2+dz**2

        if(r2<rc2) then
          rm6=1.0/r2**3
          su=su+(4.0*rm6-4.0)*rm6
          sw=sw+(24.0-48.0*rm6)*rm6
        end if

        if(r2<r2max) then
          ir=int(drm1*sqrt(r2))
          if(ir<ndr) then
            ag(ir)=ag(ir)+1_long
          end if
        end if

      end do
    end do

! Print control variables

    if(ntprint>0.and.modulo(nt,ntprint)==0) then
      write(unit=*,fmt=*) nt,naccc/max(pcre*n*ntskip,1.0), &
           naccd/max(pcre*n*ntskip,1.0), &
           naccm/max((1.0-2.0*pcre)*n*ntskip,1.0), &
           n,n/v,su/n
    end if

! Accumulate averages

    accrc=accrc+naccc/max(pcre*n*ntskip,1.0)
    accrd=accrd+naccd/max(pcre*n*ntskip,1.0)
    accrm=accrm+naccm/max((1.0-2.0*pcre)*n*ntskip,1.0)
    an=an+n
    an2=an2+n**2
    au=au+su
    aw=aw+sw

! Clear acceptance counters

    naccc=0
    naccd=0
    naccm=0

    return

  end subroutine means

!**********************************************************************!

  subroutine putcf()

!**********************************************************************!

! RNG seed

    call random_seed(get=iran)

! Write checkpoint file

    open(unit=iout,file="gcmclj_out.dat",status="replace", &
         action="write",form="unformatted")

    write(unit=iout) n,nt,ntjob,ntprint,ntskip,iran
    write(unit=iout) disp,dr,pcre,t,v,zz
    write(unit=iout) x(1:n),y(1:n),z(1:n)
    write(unit=iout) accrc,accrd,accrm,an,an2,au,aw
    write(unit=iout) ag

    close(unit=iout)

! Deallocate arrays

    deallocate(ag,iran,x,y,z)

    return

  end subroutine putcf

end module gcmclj_subm

!**********************************************************************!

program gcmclj

!**********************************************************************!

  use gcmclj_glbm
  use gcmclj_subm

  implicit none

  integer::i,j
  real::ran

! Read startup/checkpoint file & initialize variables

  call getcf()

! Do (ntskip*ntjob) passes (particle creations/deletions/moves)

  do i=1,ntjob
    nt=nt+1
    do j=1,n*ntskip
      call random_number(ran)
      if(ran<pcre) then
        call create()            ! Create particle
      else if(ran<2.0*pcre) then
        call delete()            ! Delete particle
      else
        call move()              ! Move particle
      end if
    end do
    call means()
  end do

! Write checkpoint file

  call putcf()

end program gcmclj

!**********************************************************************!
