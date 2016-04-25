!**********************************************************************!
!
! File: isomdlj.f90
!
! Isokinetic (NVT) Molecular Dynamics of Lennard-Jonesium
! Leapfrog
! Gaussian thermostat
! Pair distribution function g(r)
! Velocity autocorrelation function C(t)
!
! 03-Jun-1999 (MN)
! 05-May-2012
!
!**********************************************************************!

module isomdlj_glbm

!**********************************************************************!

! Parameters & global variables

  implicit none

  private

  integer,parameter,public::long=selected_int_kind(18)
  integer,parameter,public::double=selected_real_kind(15)
  integer,parameter,public::iin=1,iout=2

  integer,public::n,ncor,ndr,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint
  integer(kind=long),dimension(:),allocatable,public::ag
  real,public::c,c2,c2m1,dr,drm1,dt,pi,r2max,rc,rc2,rho, &
       su,su0,sw,sw0,t,zeta
  real(kind=double),public::ak,au,au2,aw
  real,dimension(:),allocatable,public::fx,fy,fz,vx,vy,vz,x,y,z
  real,dimension(:,:),allocatable,public::vxt,vyt,vzt
  real(kind=double),dimension(:),allocatable,public::acf

end module isomdlj_glbm

!**********************************************************************!

module isomdlj_subm

!**********************************************************************!

  use isomdlj_glbm

  implicit none

  private

  public::force,getcf,means,move,putcf,vacf

contains

!**********************************************************************!

  subroutine getcf()

!**********************************************************************!

    integer::i

! Pi

    pi=4.0*atan(1.0)

! Read startup/checkpoint file

    open(unit=iin,file="isomdlj_in.dat",status="old",action="read", &
         form="unformatted",position="rewind")

! Simulation parameters

    read(unit=iin) n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint
    read(unit=iin) dr,dt,rho,t

! Allocate arrays

    ndr=int(0.5*(n/rho)**(1.0/3.0)/dr)

    allocate(fx(n),fy(n),fz(n),vx(n),vy(n),vz(n),x(n),y(n),z(n), &
         vxt(n,0:ncor),vyt(n,0:ncor),vzt(n,0:ncor),acf(0:ncor), &
         ag(0:ndr-1))

! Positions & velocities: r(t),v(t-dt/2)

    read(unit=iin) x,y,z
    read(unit=iin) vx,vy,vz

! Read accumulated averages/clear accumulators

    if(nt>0) then
      read(unit=iin) ak,au,au2,aw
      read(unit=iin) ag
    else
      ak=0.0_double
      au=0.0_double
      au2=0.0_double
      aw=0.0_double
      ag=0_long
    end if

! Ring buffer & velocity autocorrelation function

    if(nt/ntcskip>0) then
      do i=0,min(nt/ntcskip-1,ncor)
        read(unit=iin) vxt(:,i),vyt(:,i),vzt(:,i)
      end do
      if(nt/ntcskip>ncor) then
        read(unit=iin) acf
      end if
    else
      acf=0.0_double
    end if

    close(unit=iin)

! Box parameters

    c=(n/rho)**(1.0/3.0)
    c2=0.5*c
    c2m1=1.0/c2
    drm1=1.0/dr
    r2max=c2**2

! Cutoff & tail corrections

    rc=c2
    rc2=rc**2
    su0=2.0*pi*rho*n*(4.0/(9.0*rc**9)-4.0/(3.0*rc**3))
    sw0=2.0*pi*rho*n*(24.0/(3.0*rc**3)-48.0/(9.0*rc**9))

    return

  end subroutine getcf

!**********************************************************************!

  subroutine force()

!**********************************************************************!

    integer::i,j
    real::df,dx,dy,dz,r2,rm2,rm6

! Initialize forces

    fx=0.0
    fy=0.0
    fz=0.0

! Calculate pair forces

    do i=1,n-1
      do j=i+1,n
        dx=x(j)-x(i)
        dx=dx-int(c2m1*dx)*c   ! Periodic boundary conditions
        dy=y(j)-y(i)
        dy=dy-int(c2m1*dy)*c
        dz=z(j)-z(i)
        dz=dz-int(c2m1*dz)*c
        r2=dx**2+dy**2+dz**2
        if(r2<rc2) then         ! In range?
          rm2=1.0/r2
          rm6=rm2**3
          df=(24.0-48.0*rm6)*rm6*rm2
          fx(i)=fx(i)+df*dx
          fx(j)=fx(j)-df*dx     ! Actio = reactio
          fy(i)=fy(i)+df*dy
          fy(j)=fy(j)-df*dy
          fz(i)=fz(i)+df*dz
          fz(j)=fz(j)-df*dz
        end if
      end do
    end do

    return

  end subroutine force

!**********************************************************************!

  subroutine move()

!**********************************************************************!

    real,dimension(n)::wx,wy,wz

! Leapfrog + Gaussian thermostat

! Intermediate velocities v(t)

    wx=vx+0.5*dt*fx
    wy=vy+0.5*dt*fy
    wz=vz+0.5*dt*fz

! Friction coefficient

    zeta=sum(wx*fx+wy*fy+wz*fz)/sum(wx**2+wy**2+wz**2)

! New velocities v(t+dt/2)

    vx=vx+dt*(fx-zeta*wx)
    vy=vy+dt*(fy-zeta*wy)
    vz=vz+dt*(fz-zeta*wz)

! New positions & PBC

    x=x+dt*vx
    x=x-int(c2m1*x)*c
    y=y+dt*vy
    y=y-int(c2m1*y)*c
    z=z+dt*vz
    z=z-int(c2m1*z)*c

    return

  end subroutine move

!**********************************************************************!

  subroutine vacf()

!**********************************************************************!

    integer::k,nt1,nt2

! Velocity autocorrelation function

! Store current values in ring buffer (index nt2)

    nt2=modulo(nt/ntcskip-1,ncor+1)
    vxt(:,nt2)=vx
    vyt(:,nt2)=vy
    vzt(:,nt2)=vz

! Correlate with values at previous time steps (index nt1) 
! k = time lag

    if(nt/ntcskip>ncor.and.modulo(nt/ntcskip-1,ntorig)==0) then
      do k=0,ncor
        nt1=modulo(nt2-k,ncor+1)
        acf(k)=acf(k)+sum(vxt(:,nt1)*vxt(:,nt2)+ &
             vyt(:,nt1)*vyt(:,nt2)+vzt(:,nt1)*vzt(:,nt2))/n
      end do
    end if

    return

  end subroutine vacf

!**********************************************************************!

  subroutine means()

!**********************************************************************!

    integer::i,j,k
    real::dx,dy,dz,r2,rm6,sk

! Kinetic energy

    sk=0.5*sum(vx**2+vy**2+vz**2)

! Potential energy, virial & g(r)

    su=su0
    sw=sw0

    do i=1,n-1
      do j=i+1,n

        dx=x(j)-x(i)
        dx=dx-int(c2m1*dx)*c            ! Periodic boundary conditions
        dy=y(j)-y(i)
        dy=dy-int(c2m1*dy)*c
        dz=z(j)-z(i)
        dz=dz-int(c2m1*dz)*c
        r2=dx**2+dy**2+dz**2

        if(r2<r2max) then               ! Inscribed sphere
          k=int(sqrt(r2)*drm1)
          if(k<ndr) then
            ag(k)=ag(k)+1_long          ! g(r)
          end if
          if(r2<rc2) then               ! In range?
            rm6=1.0/r2**3
            su=su+(4.0*rm6-4.0)*rm6     ! Potential energy
            sw=sw+(24.0-48.0*rm6)*rm6   ! Virial
          end if
        end if

      end do
    end do

! Accumulate averages

    ak=ak+sk/n
    au=au+su/n
    au2=au2+(su/n)**2
    aw=aw+sw/n

! Print control variables: friction coefficient, temperature,
! potential energy & pressure

    if(ntprint>0.and.modulo(nt/ntaskip,ntprint)==0) then
      write(unit=*,fmt=*) &
           nt,zeta,sk/(1.5*n),su/n,rho*(2.0*sk-sw)/(3.0*n)
    end if

    return

  end subroutine means

!**********************************************************************!

  subroutine putcf()

!**********************************************************************!

    integer::i

! Write checkpoint file

    open(unit=iout,file="isomdlj_out.dat",status="replace", &
         action="write",form="unformatted",position="rewind")

    write(unit=iout) n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint
    write(unit=iout) dr,dt,rho,t
    write(unit=iout) x,y,z
    write(unit=iout) vx,vy,vz
    write(unit=iout) ak,au,au2,aw
    write(unit=iout) ag

    if(nt/ntcskip>0) then
      do i=0,min(nt/ntcskip-1,ncor)
        write(unit=iout) vxt(:,i),vyt(:,i),vzt(:,i)
      end do
      if(nt/ntcskip>ncor) then
        write(unit=iout) acf
      end if
    end if

    close(unit=iout)

! Deallocate arrays

    deallocate(acf,ag,fx,fy,fz,vx,vy,vz,vxt,vyt,vzt,x,y,z)

    return

  end subroutine putcf

end module isomdlj_subm

!**********************************************************************!

program isomdlj

!**********************************************************************!

  use isomdlj_glbm
  use isomdlj_subm

  implicit none

  integer::i

! Read startup/checkpoint file

  call getcf()

! Do ntjob time steps

  do i=1,ntjob
    nt=nt+1

    call force()                   ! Leapfrog
    call move()

    if(modulo(nt,ntcskip)==0) then ! Calculate VACF with sampling rate
      call vacf()                  ! 1/(ntcskip*dt), using "time origins"
    end if                         ! (ntorig*ntcskip) time steps apart

    if(modulo(nt,ntaskip)==0) then ! Calculate averages every ntaskip
      call means()                 ! time steps
    end if

  end do

! Write checkpoint file

  call putcf()

end program isomdlj

!**********************************************************************!
