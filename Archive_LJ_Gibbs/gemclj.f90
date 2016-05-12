!**********************************************************************!
!
! File: gemclj.f90
!
! Gibbs Ensemble (NVT) Monte Carlo of Lennard-Jonesium
!
! 15-May-1999 (MN)
! 19-Apr-2012
!
!**********************************************************************!

module gemclj_glbm

!**********************************************************************!

! Parameters & global variables

  implicit none

  private

  integer,parameter,public::long=selected_int_kind(18)
  integer,parameter,public::double=selected_real_kind(15)
  integer,parameter,public::iin=1,iout=2
  real,parameter,public::r2min=0.25**2

  integer,public::n,n1,n2,naccp1,naccp2,naccv,nacct, &
       nt,ntjob,ntprint,ntryp1,ntryp2,ntryv,ntryt1,ntryt2,ntskip
  integer,dimension(:),allocatable,public::iran
  real,public::alnmax,c1,c2,disp1,disp2,dv,pi,ptr,pvm,rc1,rc2, &
       smu1,smu2,su01,su02,t,v1,v2,v
  real,dimension(:),allocatable,public::x,y,z
  real,dimension(:),allocatable,public::uattn,urepn
  real,dimension(:,:),allocatable,public::uatt,urep
  integer(kind=long),public::accrp1,accrp2,accrt,accrv, &
       atryp1,atryp2,atryt1,atryt2,atryv,an1,an2
  real(kind=double),public::amu1,ap1,arho1,au1,av1, &
       amu2,ap2,arho2,au2,av2

end module gemclj_glbm

!**********************************************************************!

module gemclj_subm

!**********************************************************************!

  use gemclj_glbm

  implicit none

  private
  public::getcf,means,move,ptrans,putcf,uinit,vmove

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

    open(unit=iin,file="gemclj_in.dat",status="old",action="read", &
         form="unformatted",position="rewind")

    read(unit=iin) n,n1,nt,ntjob,ntprint,ntskip,iran
    read(unit=iin) disp1,disp2,dv,ptr,pvm,t,v,v1

! Allocate arrays

    n2=n-n1

    allocate(uatt(n,n),uattn(n),urep(n,n),urepn(n),x(n),y(n),z(n))

! Positions

    read(unit=iin) x,y,z

! Read accumulated averages/clear accumulators

    if(nt>0) then
      read(unit=iin) accrp1,accrp2,accrt,accrv, &
           atryp1,atryp2,atryt1,atryt2,atryv, &
           amu1,an1,ap1,arho1,au1,av1,amu2,an2,ap2,arho2,au2,av2
    else
      accrp1=0_long
      accrp2=0_long
      accrt=0_long
      accrv=0_long
      atryp1=0_long
      atryp2=0_long
      atryt1=0_long
      atryt2=0_long
      atryv=0_long
      amu1=0.0_double
      an1=0_long
      ap1=0.0_double
      arho1=0.0_double
      au1=0.0_double
      av1=0.0_double
      amu2=0.0_double
      an2=0_long
      ap2=0.0_double
      arho2=0.0_double
      au2=0.0_double
      av2=0.0_double
    end if

    close(unit=iin)

! Box parameters

    v2=v-v1
    c1=v1**(1.0/3.0)
    c2=v2**(1.0/3.0)

! Cutoff & tail corrections

    rc1=0.5*c1
    rc2=0.5*c2
    su01=2.0*pi*n1**2/v1*(4.0/(9.0*rc1**9)-4.0/(3.0*rc1**3))
    su02=2.0*pi*n2**2/v2*(4.0/(9.0*rc2**9)-4.0/(3.0*rc2**3))

! Initialize RNG & acceptance counters

    call random_seed(put=iran)

    ntryp1=0
    ntryp2=0
    ntryv=0
    ntryt1=0
    ntryt2=0
    naccp1=0
    naccp2=0
    naccv=0
    nacct=0

    smu1=0.0
    smu2=0.0

    return

  end subroutine getcf

!**********************************************************************!

  subroutine uinit()

!**********************************************************************!

    integer::i,ibox,imin,imax,j
    real::c,dx,dy,dz,f2c,r2,rcrc,rm6

! Initialize pair interaction matrix
!
! urep = 4/r**12, uatt = -4/r**6

    urep=0.0
    uatt=0.0

! Do box 1 & 2

    do ibox=1,2

      if(ibox==1) then
        imin=1
        imax=n1
        c=c1
        rcrc=rc1**2
      else
        imin=n1+1
        imax=n
        c=c2
        rcrc=rc2**2
      end if
      f2c=2.0/c

      do i=imin,imax-1
        do j=i+1,imax
          dx=x(j)-x(i)
          dx=dx-int(f2c*dx)*c
          dy=y(j)-y(i)
          dy=dy-int(f2c*dy)*c
          dz=z(j)-z(i)
          dz=dz-int(f2c*dz)*c
          r2=dx**2+dy**2+dz**2
          if(r2<rcrc) then
            rm6=1.0/r2**3
            urep(i,j)=4.0*rm6**2
            uatt(i,j)=-4.0*rm6
            urep(j,i)=urep(i,j)
            uatt(j,i)=uatt(i,j)
          end if
        end do
      end do

    end do

    return

  end subroutine uinit

!**********************************************************************!

  subroutine move()

!**********************************************************************!

    logical::accept
    integer::i,j,jmin,jmax
    real::c,disp,dx,dy,dz,f2c,r2,ran,rcrc,rm6,su,xin,yin,zin

! Random particle

    call random_number(ran)
    i=min(int(ran*n+1.0),n)

! Which box?

    if(i<=n1) then
      ntryp1=ntryp1+1     ! Box 1
      jmin=1
      jmax=n1
      c=c1
      disp=disp1
      rcrc=rc1**2
    else
      ntryp2=ntryp2+1     ! Box 2
      jmin=n1+1
      jmax=n
      c=c2
      disp=disp2
      rcrc=rc2**2
    end if
    f2c=2.0/c

! Trial move (displacement)

    call random_number(ran)
    xin=x(i)+disp*(ran-0.5)
    xin=xin-int(f2c*xin)*c
    call random_number(ran)
    yin=y(i)+disp*(ran-0.5)
    yin=yin-int(f2c*yin)*c
    call random_number(ran)
    zin=z(i)+disp*(ran-0.5)
    zin=zin-int(f2c*zin)*c

! Energy difference

    su=0.0
    do j=jmin,jmax
      if(j==i) then
        urepn(j)=0.0
        uattn(j)=0.0
      else
        dx=x(j)-xin
        dx=dx-int(f2c*dx)*c
        dy=y(j)-yin
        dy=dy-int(f2c*dy)*c
        dz=z(j)-zin
        dz=dz-int(f2c*dz)*c
        r2=dx**2+dy**2+dz**2
        if(r2<rcrc) then
          rm6=1.0/r2**3
          urepn(j)=4.0*rm6**2
          uattn(j)=-4.0*rm6
        else
          urepn(j)=0.0
          uattn(j)=0.0
        end if
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

! Accept move & update interaction matrix

    if(accept) then
      if(i<=n1) then
        naccp1=naccp1+1
      else
        naccp2=naccp2+1
      end if
      x(i)=xin
      y(i)=yin
      z(i)=zin
      urep(i,jmin:jmax)=urepn(jmin:jmax)
      urep(jmin:jmax,i)=urepn(jmin:jmax)
      uatt(i,jmin:jmax)=uattn(jmin:jmax)
      uatt(jmin:jmax,i)=uattn(jmin:jmax)
    end if

    return

  end subroutine move

!**********************************************************************!

  subroutine ptrans()

!**********************************************************************!

    logical::accept
    integer::i,iold,j,jmin,jmax
    real::c,dx,dy,dz,f2c,r2,ran,rcrc,rm6,su,xin,yin,zin

! Direction of transfer

    call random_number(ran)
    if(ran<0.5) then
      ntryt1=ntryt1+1                   ! Box 2 -> box 1 (n1 -> n1+1)
      if(n1>=n) then
        i=n+1                           ! Box 1 full
      else
        call random_number(ran)
        i=min(n1+int(ran*n2+1.0),n)     ! Delete i from box 2
        iold=n1+1                       ! Add as n1+1 to box 1
      end if
    else
      ntryt2=ntryt2+1                   ! Box 1 -> box 2 (n1 -> n1-1)
      if(n1<=0) then
        i=0                             ! Box 2 full
      else
        call random_number(ran)
        i=min(int(ran*n1+1.0),n1)       ! Delete i from box 1
        iold=n1                         ! Add as n1 to box 2
      end if
    end if

! Test particle energy

    if(i>n1) then
      jmin=1          ! Box 1
      jmax=n1
      c=c1
      rcrc=rc1**2
      su=2.0*pi*(2*n1+1)/v1*(4.0/(9.0*rc1**9)-4.0/(3.0*rc1**3))
    else
      jmin=n1+1       ! Box 2
      jmax=n
      c=c2
      rcrc=rc2**2
      su=2.0*pi*(2*n2+1)/v2*(4.0/(9.0*rc2**9)-4.0/(3.0*rc2**3))
    end if
    f2c=2.0/c

    call random_number(ran)     ! Trial position
    xin=c*(ran-0.5)
    call random_number(ran)
    yin=c*(ran-0.5)
    call random_number(ran)
    zin=c*(ran-0.5)

    do j=jmin,jmax              ! Test particle energy
      dx=x(j)-xin
      dx=dx-int(f2c*dx)*c
      dy=y(j)-yin
      dy=dy-int(f2c*dy)*c
      dz=z(j)-zin
      dz=dz-int(f2c*dz)*c
      r2=dx**2+dy**2+dz**2
      if(r2<rcrc) then
        if(r2<r2min) then     ! Overlap?
          return
        else
          rm6=1.0/r2**3
          urepn(j)=4.0*rm6**2
          uattn(j)=-4.0*rm6
          su=su+(urepn(j)+uattn(j))
        end if
      else
        urepn(j)=0.0
        uattn(j)=0.0
      end if
    end do

! Chemical potential

    if(su/t<alnmax) then
      if(i>n1) then                              ! Box 1
        if(n1<n) then
          smu1=smu1+exp(-su/t)*v1/(n1+1)
        else
          smu1=smu1+2.0*exp(-su/t)*v1/(n1+1)
        end if
      else                                       ! Box 2
        if(n2<n) then
          smu2=smu2+exp(-su/t)*v2/(n2+1)
        else
          smu2=smu2+2.0*exp(-su/t)*v2/(n2+1)
        end if
      end if
    end if

! Virtual transfer only?

    if(i<1.or.i>n) then
      return
    end if

! Complete trial "energy"

    if(i>n1) then     ! Delete i from box 2
      su=su+2.0*pi*(1-2*n2)/v2*(4.0/(9.0*rc2**9)-4.0/(3.0*rc2**3))- &
           sum(urep(i,n1+1:n)+uatt(i,n1+1:n))+ &
           t*log((n1+1)*v2/(n2*v1))
    else              ! Delete i from box 1
      su=su+2.0*pi*(1-2*n1)/v1*(4.0/(9.0*rc1**9)-4.0/(3.0*rc1**3))- &
           sum(urep(i,1:n1)+uatt(i,1:n1))+ &
           t*log((n2+1)*v1/(n1*v2))
    end if

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

! Accept transfer

    if(accept) then
      nacct=nacct+1

      x(i)=x(iold)     ! Swap positions
      y(i)=y(iold)
      z(i)=z(iold)
      x(iold)=xin
      y(iold)=yin
      z(iold)=zin

      urep(i,:)=urep(iold,:)     ! Update pair interaction matrix
      urep(:,i)=urep(i,:)
      urep(i,i)=0.0
      urepn(i)=urepn(iold)
      urep(iold,:)=urepn
      urep(:,iold)=urepn
      urep(iold,iold)=0.0
      uatt(i,:)=uatt(iold,:)
      uatt(:,i)=uatt(i,:)
      uatt(i,i)=0.0
      uattn(i)=uattn(iold)
      uatt(iold,:)=uattn
      uatt(:,iold)=uattn
      uatt(iold,iold)=0.0

      if(i>n1) then     ! Particle numbers & tail corrections
        n1=n1+1
      else
        n1=n1-1
      end if
      n2=n-n1
      su01=2.0*pi*n1**2/v1*(4.0/(9.0*rc1**9)-4.0/(3.0*rc1**3))
      su02=2.0*pi*n2**2/v2*(4.0/(9.0*rc2**9)-4.0/(3.0*rc2**3))

    end if

    return

  end subroutine ptrans

!**********************************************************************!

  subroutine vmove()

!**********************************************************************!

    logical::accept
    integer::i,j
    real::fact,fatt,frep,ran,rcn,su,su0n,v1n,v2n

    ntryv=ntryv+1

! Volume change

    call random_number(ran)
    v1n=v1+(ran-0.5)*dv
    if(v1n<=0.0.or.v1n>=v) then
      return
    end if
    v2n=v-v1n

! "Energy" difference: Box 1

    rcn=0.5*v1n**(1.0/3.0)     ! New cutoff & tail correction
    su0n=2.0*pi*n1**2/v1n*(4.0/(9.0*rcn**9)-4.0/(3.0*rcn**3))

    frep=(v1/v1n)**4-1.0
    fatt=(v1/v1n)**2-1.0
    su=su0n-su01
    do i=1,n1-1
      do j=i+1,n1
        su=su+(frep*urep(i,j)+fatt*uatt(i,j))
      end do
    end do

    su=su-t*n1*log(v1n/v1)

! Box 2

    rcn=0.5*v2n**(1.0/3.0)     ! New cutoff & tail correction
    su0n=2.0*pi*n2**2/v2n*(4.0/(9.0*rcn**9)-4.0/(3.0*rcn**3))

    frep=(v2/v2n)**4-1.0
    fatt=(v2/v2n)**2-1.0
    su=su+(su0n-su02)
    do i=n1+1,n-1
      do j=i+1,n
        su=su+(frep*urep(i,j)+fatt*uatt(i,j))
      end do
    end do

    su=su-t*n2*log(v2n/v2)

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

      fact=(v1n/v1)**(1.0/3.0)     ! Box 1
      x(1:n1)=fact*x(1:n1)
      y(1:n1)=fact*y(1:n1)
      z(1:n1)=fact*z(1:n1)
      fact=(v2n/v2)**(1.0/3.0)     ! Box 2
      x(n1+1:n)=fact*x(n1+1:n)
      y(n1+1:n)=fact*y(n1+1:n)
      z(n1+1:n)=fact*z(n1+1:n)

! Parameters

      v1=v1n               ! Box 1
      c1=v1**(1.0/3.0)
      v2=v2n               ! Box 2
      c2=v2**(1.0/3.0)

! Cutoff & tail corrections

      rc1=0.5*c1                                                 ! Box 1
      su01=2.0*pi*n1**2/v1*(4.0/(9.0*rc1**9)-4.0/(3.0*rc1**3))
      rc2=0.5*c2                                                 ! Box 2
      su02=2.0*pi*n2**2/v2*(4.0/(9.0*rc2**9)-4.0/(3.0*rc2**3))

! Pair interaction matrix

      call uinit()

    end if

    return

  end subroutine vmove

!**********************************************************************!

  subroutine means()

!**********************************************************************!

    integer::i,j
    real::c,dx,dy,dz,f2c,r2,rcrc,rho1,rho2,rm6,su1,su2,sw1,sw2

! Potential energy & virial

    su1=su01
    su2=su02
    sw1=2.0*pi*n1**2/v1*(24.0/(3.0*rc1**3)-48.0/(9.0*rc1**9))
    sw2=2.0*pi*n2**2/v2*(24.0/(3.0*rc2**3)-48.0/(9.0*rc2**9))

! Box 1

    c=c1
    f2c=2.0/c
    rcrc=rc1**2
    do i=1,n1-1
      do j=i+1,n1
        dx=x(j)-x(i)
        dx=dx-int(f2c*dx)*c
        dy=y(j)-y(i)
        dy=dy-int(f2c*dy)*c
        dz=z(j)-z(i)
        dz=dz-int(f2c*dz)*c
        r2=dx**2+dy**2+dz**2
        if(r2<rcrc) then
          su1=su1+(urep(i,j)+uatt(i,j))
          rm6=1.0/r2**3
          sw1=sw1+(24.0-48.0*rm6)*rm6
        end if
      end do
    end do
    su1=su1/max(n1,1)

! Box 2

    c=c2
    f2c=2.0/c
    rcrc=rc2**2
    do i=n1+1,n-1
      do j=i+1,n
        dx=x(j)-x(i)
        dx=dx-int(f2c*dx)*c
        dy=y(j)-y(i)
        dy=dy-int(f2c*dy)*c
        dz=z(j)-z(i)
        dz=dz-int(f2c*dz)*c
        r2=dx**2+dy**2+dz**2
        if(r2<rcrc) then
          su2=su2+(urep(i,j)+uatt(i,j))
          rm6=1.0/r2**3
          sw2=sw2+(24.0-48.0*rm6)*rm6
        end if
      end do
    end do
    su2=su2/max(n2,1)

! Print control variables

    rho1=n1/v1
    rho2=n2/v2

    if(ntprint>0.and.modulo(nt,ntprint)==0) then
      write(unit=*,fmt="(tr1,i10,8(tr1,es12.5))") nt, &
           real(naccp1)/max(ntryp1,1),real(naccp2)/max(ntryp2,1), &
           real(nacct)/max(ntryt1+ntryt2,1),real(naccv)/max(ntryv,1), &
           su1,su2,rho1,rho2
    end if

! Accumulate averages

    accrp1=accrp1+naccp1
    accrp2=accrp2+naccp2
    accrv=accrv+naccv
    accrt=accrt+nacct
    atryp1=atryp1+ntryp1
    atryp2=atryp2+ntryp2
    atryt1=atryt1+ntryt1
    atryt2=atryt2+ntryt2
    atryv=atryv+ntryv
    an1=an1+n1
    an2=an2+n2

    amu1=amu1+smu1
    ap1=ap1+(t*rho1-sw1/(3.0*v1))
    arho1=arho1+rho1
    au1=au1+su1
    av1=av1+v1

    amu2=amu2+smu2
    ap2=ap2+(t*rho2-sw2/(3.0*v2))
    arho2=arho2+rho2
    au2=au2+su2
    av2=av2+v2

! Clear acceptance counters & accumulators

    ntryp1=0
    ntryp2=0
    ntryt1=0
    ntryt2=0
    ntryv=0
    naccp1=0
    naccp2=0
    nacct=0
    naccv=0

    smu1=0.0
    smu2=0.0

    return

  end subroutine means

!**********************************************************************!

  subroutine putcf()

!**********************************************************************!

! RNG seed

    call random_seed(get=iran)

! Write checkpoint file

    open(unit=iout,file="gemclj_out.dat",status="replace", &
         action="write",form="unformatted")

    write(unit=iout) n,n1,nt,ntjob,ntprint,ntskip,iran
    write(unit=iout) disp1,disp2,dv,ptr,pvm,t,v,v1
    write(unit=iout) x,y,z
    write(unit=iout) accrp1,accrp2,accrt,accrv, &
         atryp1,atryp2,atryt1,atryt2,atryv, &
         amu1,an1,ap1,arho1,au1,av1,amu2,an2,ap2,arho2,au2,av2

    close(unit=iout)

! Deallocate arrays

    deallocate(iran,uatt,uattn,urep,urepn,x,y,z)

    return

  end subroutine putcf

end module gemclj_subm

!**********************************************************************!

program gemclj

!**********************************************************************!

  use gemclj_glbm
  use gemclj_subm

  implicit none

  integer::i,j
  real::ran

! Read startup/checkpoint file & initialize

  call getcf()
  call uinit()

! Do (ntskip*ntjob) passes (transfers/volume changes/displacements)

  do i=1,ntjob
    nt=nt+1
    do j=1,ntskip*n
      call random_number(ran)
      if(ran<ptr) then
        call ptrans()               ! Particle transfer
      else if(ran<ptr+pvm) then
        call vmove()                ! Volume change
      else
        call move()                 ! Particle displacement
      end if
    end do
    call means()
  end do

! Write checkpoint file

  call putcf()

end program gemclj

!**********************************************************************!
