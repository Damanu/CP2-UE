!**********************************************************************!
!
! File: mclj.f90
!
! NVT-Monte Carlo of Lennard-Jonesium
!
! 13-Apr-1997 (MN)
! 19-Apr-2012
!
!**********************************************************************!

module mclj_glbm

!**********************************************************************!

! Global variables

  implicit none

  private

  integer,parameter,public::long=selected_int_kind(18)
  integer,parameter,public::double=selected_real_kind(15)
  integer,parameter,public::iin=1,iout=2

  integer,public::n,nacc,ndr,nt,ntjob,ntprint,ntskip
  integer,dimension(:),allocatable,public::iran
  real,public::alnmax,c,c2,c2m1,disp,dr,drm1,pi,r2max,rc,rc2,rho, &
       su0,sw0,t
  real,dimension(:),allocatable,public::x,y,z
  real,dimension(:,:),allocatable,public::umat
  integer(kind=long),dimension(:),allocatable,public::ag
  real(kind=double),public::accr,au,au2,aw

end module mclj_glbm

!**********************************************************************!

module mclj_subm

!**********************************************************************!

  use mclj_glbm

  implicit none

  private

  public::getcf,uinit,move,means,putcf

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

    open(unit=iin,file="mclj_in.dat",status="old",action="read", &
         form="unformatted",position="rewind")

    read(unit=iin) n,nt,ntjob,ntprint,ntskip,iran
    read(unit=iin) disp,dr,rho,t

! Box parameters

    c=(n/rho)**(1.0/3.0)
    c2=0.5*c
    c2m1=1.0/c2
    r2max=c2**2
    drm1=1.0/dr
    ndr=int(c2*drm1)

! Allocate arrays

    allocate(ag(0:ndr-1),umat(n,n),x(n),y(n),z(n))

! Positions & accumulated averages

    read(unit=iin) x,y,z

    if(nt>0) then
      read(unit=iin) accr,au,au2,aw
      read(unit=iin) ag
    else
      accr=0.0_double
      au=0.0_double
      au2=0.0_double
      aw=0.0_double
      ag=0.0_long
    end if

    close(unit=iin)

! Cutoff & tail corrections

    rc=c2
    rc2=rc**2
    su0=2.0*pi*rho*n*(4.0/(9.0*rc**9)-4.0/(3.0*rc**3))
    sw0=2.0*pi*rho*n*(24.0/(3.0*rc**3)-48.0/(9.0*rc**9))

! Initialize RNG & acceptance counter

    call random_seed(put=iran)
    nacc=0

    return

  end subroutine getcf

!**********************************************************************!

  subroutine uinit()

!**********************************************************************!

    integer::i,j
    real::dx,dy,dz,r2,rm6

! Initialize pair interaction matrix

    do i=1,n-1
      umat(i,i)=0.0
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
          umat(i,j)=(4.0*rm6-4.0)*rm6
        else
          umat(i,j)=0.0
        end if
        umat(j,i)=umat(i,j)
      end do
    end do
    umat(n,n)=0.0

    return

  end subroutine uinit

!**********************************************************************!

  subroutine means()

!**********************************************************************!

    integer::i,j,k
    real::dx,dy,dz,p,r2,rm6,su,sw

! Potential energy, virial & g(r)

    su=su0
    sw=sw0
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
          su=su+umat(i,j)               ! Potential energy
          sw=sw+(24.0-48.0*rm6)*rm6     ! Virial
        end if
        if(r2<r2max) then
          k=int(sqrt(r2)*drm1)
          if(k<ndr) then
            ag(k)=ag(k)+1_long          ! g(r)
          end if
        end if
      end do
    end do
    su=su/n
    sw=sw/n

! Print control variables

    if(ntprint>0.and.modulo(nt,ntprint)==0) then
      p=rho*(t-sw/3.0)
      write(unit=*,fmt="(tr1,i10,3(tr1,es12.5))") &
           nt,real(nacc)/(n*ntskip),su,p
    end if

! Accumulate averages

    accr=accr+real(nacc)/(n*ntskip)
    au=au+su
    au2=au2+su**2
    aw=aw+sw

! Clear acceptance counter

    nacc=0

    return

  end subroutine means

!**********************************************************************!

  subroutine move()

!**********************************************************************!

    logical::accept
    integer::i,j
    real::dx,dy,dz,r2,ran,rm6,su,xin,yin,zin
    real,dimension(n)::un

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
        un(j)=0.0
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
          un(j)=(4.0*rm6-4.0)*rm6
        else
          un(j)=0.0
        end if
        su=su+(un(j)-umat(i,j))
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
      nacc=nacc+1
      x(i)=xin
      y(i)=yin
      z(i)=zin
      umat(i,:)=un
      umat(:,i)=un
    end if

    return

  end subroutine move

!**********************************************************************!

  subroutine putcf()

!**********************************************************************!

! RNG seed

    call random_seed(get=iran)

! Write checkpoint file

    open(unit=iout,file="mclj_out.dat",status="replace",action="write", &
         form="unformatted")

    write(unit=iout) n,nt,ntjob,ntprint,ntskip,iran
    write(unit=iout) disp,dr,rho,t
    write(unit=iout) x,y,z
    write(unit=iout) accr,au,au2,aw
    write(unit=iout) ag

    close(unit=iout)

! Deallocate arrays

    deallocate(ag,iran,umat,x,y,z)

    return

  end subroutine putcf

end module mclj_subm

!**********************************************************************!

program mclj

!**********************************************************************!

  use mclj_glbm
  use mclj_subm

  implicit none

  integer::i,j

! Read startup/checkpoint file & initialize

  call getcf()
  call uinit()

! Do (ntskip*ntjob) passes

  do i=1,ntjob
    nt=nt+1
    do j=1,ntskip*n
      call move()
    end do
    call means()
  end do

! Write checkpoint file

  call putcf()

end program mclj

!**********************************************************************!
