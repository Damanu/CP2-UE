!**********************************************************************!
!
! File: nptmclj.f90
!
! NpT-Monte Carlo of Lennard-Jonesium
!
! 29-Apr-1999 (MN)
! 19-Apr-2010
!
!**********************************************************************!

module nptmclj_glbm

!**********************************************************************!

! Parameters & global variables

  implicit none

  private

  integer,parameter,public::long=selected_int_kind(18)
  integer,parameter,public::double=selected_real_kind(15)
  integer,parameter,public::iin=1,iout=2

  integer,public::n,nt,ntjob,ntprint,ntskip
  integer,public::naccp,naccv,ndr,ntryp,ntryv
  integer,dimension(:),allocatable,public::iran
  real,public::disp,dlnv,dr,p,pvm,t,v
  real,public::alnmax,c,c2,c2m1,drm1,pi,r2max,rc,rc2,su0
  real,dimension(:),allocatable,public::x,y,z
  real,dimension(:,:),allocatable,public::uatt,urep
  integer(kind=long),dimension(:),allocatable,public::ag
  real(kind=double),public::accrp,accrv,arho,au,aupv2,av,av2

end module nptmclj_glbm

!**********************************************************************!

module nptmclj_subm

  use nptmclj_glbm

  implicit none

  private

  public::getcf,means,pmove,putcf,vmove,uinit

contains

!**********************************************************************!

  subroutine getcf()

!**********************************************************************!

    integer::nran

! Maximum argument for exponential & pi

    alnmax=log(0.1*huge(1.0))
    pi=4.0*atan(1.0)

! RNG characteristics

    call random_seed(size=nran)

    allocate(iran(nran))

! Read startup/checkpoint file

    open(unit=iin,file="nptmclj_in.dat",status="old", &
         action="read",form="unformatted",position="rewind")

    read(unit=iin) n,nt,ntjob,ntprint,ntskip,iran
    read(unit=iin) disp,dlnv,dr,p,pvm,t,v

! Allocate arrays

    ndr=int(0.5*v**(1.0/3.0)/dr)

    allocate(ag(0:ndr-1),uatt(n,n),urep(n,n),x(n),y(n),z(n))

! Positions

    read(unit=iin) x,y,z

! Accumulated averages/clear accumulators

    if(nt>0) then
      read(unit=iin) accrp,accrv,arho,au,aupv2,av,av2
      read(unit=iin) ndr
      read(unit=iin) ag(0:ndr-1)
    else
      accrp=0.0_double
      accrv=0.0_double
      arho=0.0_double
      au=0.0_double
      aupv2=0.0_double
      av=0.0_double
      av2=0.0_double
      ag=0_long
    end if

    close(unit=iin)

! Box parameters

    c=v**(1.0/3.0)
    c2=0.5*c
    c2m1=1.0/c2
    r2max=c2**2
    drm1=1.0/dr

! Cutoff & tail correction

    rc=c2
    rc2=rc**2
    su0=2.0*pi*n**2/v*(4.0/(9.0*rc**9)-4.0/(3.0*rc**3))

! Initialize RNG & acceptance counters

    call random_seed(put=iran)

    ntryp=0
    ntryv=0
    naccp=0
    naccv=0

    return

  end subroutine getcf

!**********************************************************************!

  subroutine uinit()

!**********************************************************************!

    integer::i,j
    real::dx,dy,dz,r2,rm6

! Initialize pair interaction matrix
!
! urep = 4/r**12, uatt = -4/r**6

    do i=1,n-1
      urep(i,i)=0.0
      uatt(i,i)=0.0
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
          urep(i,j)=4.0*rm6**2
          uatt(i,j)=-4.0*rm6
        else
          urep(i,j)=0.0
          uatt(i,j)=0.0
        end if
        urep(j,i)=urep(i,j)
        uatt(j,i)=uatt(i,j)
      end do
    end do
    urep(n,n)=0.0
    uatt(n,n)=0.0

    return

  end subroutine uinit

!**********************************************************************!

  subroutine pmove()

!**********************************************************************!

    logical::accept
    integer::i,j
    real::dx,dy,dz,r2,ran,rm6,su,xin,yin,zin
    real,dimension(n)::uattn,urepn

    ntryp=ntryp+1

! Random particle

    call random_number(ran)
    i=min(int(ran*n+1.0),n)

! Trial move

    call random_number(ran)
    xin=x(i)+disp*(ran-0.5)
    xin=xin-int(xin*c2m1)*c
    call random_number(ran)
    yin=y(i)+disp*(ran-0.5)
    yin=yin-int(yin*c2m1)*c
    call random_number(ran)
    zin=z(i)+disp*(ran-0.5)
    zin=zin-int(zin*c2m1)*c

! Energy difference

    su=0.0
    do j=1,n
      if(j==i) then
        urepn(j)=0.0
        uattn(j)=0.0
      else
        dx=x(j)-xin
        dx=dx-int(dx*c2m1)*c
        dy=y(j)-yin
        dy=dy-int(dy*c2m1)*c
        dz=z(j)-zin
        dz=dz-int(dz*c2m1)*c
        r2=dx**2+dy**2+dz**2
        if(r2<rc2) then
          rm6=1.0/r2**3
          urepn(j)=4.0*rm6**2
          uattn(j)=-4.0*rm6
        else
          urepn(j)=0.0
          uattn(j)=0.0
        endif
        su=su+((urepn(j)+uattn(j))-(urep(i,j)+uatt(i,j)))
      end if
    end do

! Acceptance test

    if(su<=0.0) then
      accept=.true.
    else
      if(su/t<alnmax) then
        call random_number(ran)
        if(ran<=exp(-su/t)) then
          accept=.true.
        else
          accept=.false.
        end if
      else
        accept=.false.
      end if
    end if

! Accept move & update pair interaction matrix

    if(accept) then
      naccp=naccp+1
      x(i)=xin
      y(i)=yin
      z(i)=zin
      urep(i,:)=urepn
      urep(:,i)=urepn
      uatt(i,:)=uattn
      uatt(:,i)=uattn
    end if

    return

  end subroutine pmove

!**********************************************************************!

  subroutine vmove()

!**********************************************************************!

    logical::accept
    integer::i,j
    real::fact,fatt,frep,ran,rcn,su,sun0,vn

    ntryv=ntryv+1

! Volume change

    call random_number(ran)
    vn=v*exp(dlnv*(ran-0.5))

! New cutoff & tail correction

    rcn=0.5*vn**(1.0/3.0)
    sun0=2.0*pi*n**2/vn*(4.0/(9.0*rcn**9)-4.0/(3.0*rcn**3))

! "Energy" difference

    frep=(v/vn)**4-1.0
    fatt=(v/vn)**2-1.0
    su=sun0-su0
    do i=1,n-1
      do j=i+1,n
        su=su+(frep*urep(i,j)+fatt*uatt(i,j))
      end do
    end do
    su=su+p*(vn-v)-(n+1)*t*log(vn/v)

! Acceptance test

    if(su<=0.0) then
      accept=.true.
    else
      if(su/t<alnmax) then
        call random_number(ran)
        if(ran<=exp(-su/t)) then
          accept=.true.
        else
          accept=.false.
        end if
      else
        accept=.false.
      end if
    end if

! Accept

    if(accept) then
      naccv=naccv+1

! Positions

      fact=(vn/v)**(1.0/3.0)
      x=fact*x
      y=fact*y
      z=fact*z

! Box parameters

      v=vn
      c=v**(1.0/3.0)
      c2=0.5*c
      c2m1=1.0/c2
      r2max=c2**2
      ndr=min(ndr,int(c2/dr))

! Cutoff & tail correction

      rc=rcn
      rc2=rc**2
      su0=sun0

! Pair interaction matrix

      call uinit()

    end if

    return

  end subroutine vmove

!**********************************************************************!

  subroutine means()

!**********************************************************************!

    integer::i,j,k
    real::dx,dy,dz,r2,su

! Potential energy

    su=(su0+0.5*sum(urep+uatt))/n

! g(r)

    do i=1,n-1
      do j=i+1,n
        dx=x(j)-x(i)
        dx=dx-int(dx*c2m1)*c
        dy=y(j)-y(i)
        dy=dy-int(dy*c2m1)*c
        dz=z(j)-z(i)
        dz=dz-int(dz*c2m1)*c
        r2=dx**2+dy**2+dz**2
        if(r2<r2max) then
          k=int(sqrt(r2)*drm1)
          if(k<ndr) then
            ag(k)=ag(k)+1_long
          end if
        end if
      end do
    end do

! Print control variables

    if(ntprint>0.and.modulo(nt,ntprint)==0) then
      write(unit=*,fmt="(tr1,i10,5(tr1,es12.5))") nt, &
           real(naccp)/max(ntryp,1),real(naccv)/max(ntryv,1), &
           su,v,n/v
    end if

! Accumulate averages

    accrp=accrp+real(naccp)/max(ntryp,1)
    accrv=accrv+real(naccv)/max(ntryv,1)
    arho=arho+(n/v)
    au=au+su
    aupv2=aupv2+(su+p*v/n)**2
    av=av+v
    av2=av2+v**2

! Clear acceptance counters

    naccp=0
    naccv=0
    ntryp=0
    ntryv=0

    return

  end subroutine means

!**********************************************************************!

  subroutine putcf()

!**********************************************************************!

! RNG seed

    call random_seed(get=iran)

! Write checkpoint file

    open(unit=iout,file="nptmclj_out.dat",status="replace", &
         action="write",form="unformatted")

    write(unit=iout) n,nt,ntjob,ntprint,ntskip,iran
    write(unit=iout) disp,dlnv,dr,p,pvm,t,v
    write(unit=iout) x,y,z
    write(unit=iout) accrp,accrv,arho,au,aupv2,av,av2
    write(unit=iout) ndr
    write(unit=iout) ag(0:ndr-1)

    close(unit=iout)

! Deallocate arrays

    deallocate(ag,iran,uatt,urep,x,y,z)

    return

  end subroutine putcf

end module nptmclj_subm

!**********************************************************************!

program nptmclj

!**********************************************************************!

  use nptmclj_glbm
  use nptmclj_subm

  implicit none

  integer::i,j
  real::ran

! Read startup/checkpoint file & initialize pair interaction matrix

  call getcf()
  call uinit()

! Do ntskip*ntjob passes/volume changes

  do i=1,ntjob
    nt=nt+1
    do j=1,ntskip*n
      call random_number(ran)
      if(ran<pvm) then
        call vmove()     ! Volume change
      else
        call pmove()     ! Particle move
      end if
    end do
    call means()
  end do

! Write checkpoint file

  call putcf()

end program nptmclj

!**********************************************************************!
