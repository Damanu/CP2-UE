!**********************************************************************!
!
! File: mcme.f90
!
! NVT-Monte Carlo of OPLS-AA Methane
!
! Potential parameters and molecular geometry from:
!   G. Kaminski, E. M. Duffy, T. Matsui & W. L. Jorgensen,
!   J. Phys. Chem. 98, 13077 (1994)
!
! 23-Apr-1999 (MN)
! 19-Apr-2012
!
!**********************************************************************!
!
!        H4----------+
!       /|          /|
!      / |         / |
!     +-----------H5 |
!     |  |        |  |
!     |  |    C1  |  |
!     |  |        |  |
!     |  +--------|--H3
!     | /         | /
!     |/          |/
!     H2----------+
!
!**********************************************************************!
!
! Convention for g_AB(r):
!
!   "1"=CC   "2"=CH   "3"=HH
!
!**********************************************************************!

module mcme_glbm

!**********************************************************************!

  implicit none

  private

  integer,parameter,public::long=selected_int_kind(18)
  integer,parameter,public::double=selected_real_kind(15)
  integer,parameter,public::iin=1,iout=2,m=5,ntab=32768
  real,parameter,public::avog=602.252,atm=0.101325*(avog/1000.0), &
       qel2=0.480298**2*avog,rg=8.3143e-03
  real,parameter,public::dch=0.109,xmc=12.0,xmh=1.0
  real,parameter,public::qh=0.06,qc=-4.0*qh
  real,parameter,public::epscc=33.2*rg,epshh=15.1*rg
  real,parameter,public::sigcc=0.35,sighh=0.25
  real,parameter,public::tolbond=1.0e-05,tolbox=1.001
  real,parameter,dimension(m),public::xmm=(/xmc,xmh,xmh,xmh,xmh/)

  integer,public::n,nacc,ndr,nt,ntjob,ntprint,ntskip
  integer,dimension(:),allocatable,public::iran
  real,public::alnmax,c,c2,c22,c2m,cm,da,dcx,disp,dr,drm1,dr2m1, &
       pi,rho,rc,rc2,su0,sw0,t
  real,dimension(m),public::ufact,vfact,xm,ym,zm
  real,dimension(m,m),public::a12,c6,qq
  real,dimension(m,m,0:ntab),public::ftab,utab
  real,dimension(:),allocatable,public::xcm,ycm,zcm
  real,dimension(:,:),allocatable,public::umat,x,y,z
  integer(kind=long),dimension(:,:),allocatable,public::ag
  real(kind=double),public::accr,au,au2,aw

end module mcme_glbm

!**********************************************************************!

module mcme_subm

!**********************************************************************!

  use mcme_glbm

  implicit none

  private
  public::getcf,means,mgeom,move,pbc,putcf,tables,uinit

contains

!**********************************************************************!

  subroutine mgeom()

!**********************************************************************!

    real,dimension(m),parameter::q=(/qc,qh,qh,qh,qh/)
    real,dimension(m),parameter::eps=(/epscc,epshh,epshh,epshh,epshh/)
    real,dimension(m),parameter::sig=(/sigcc,sighh,sighh,sighh,sighh/)

    real::xcm,ycm,zcm

! Potential parameters

    qq=spread(q,2,m)*spread(q,1,m)*qel2
    a12=sqrt(spread(4.0*eps*sig**12,2,m)*spread(4.0*eps*sig**12,1,m))
    c6=sqrt(spread(4.0*eps*sig**6,2,m)*spread(4.0*eps*sig**6,1,m))

! Molecular geometry

    xm(1)=0.0                ! C1
    ym(1)=0.0
    zm(1)=0.0

    xm(2)=-dch/sqrt(3.0)     ! H2
    ym(2)=-dch/sqrt(3.0)
    zm(2)=-dch/sqrt(3.0)
    xm(3)=dch/sqrt(3.0)      ! H3
    ym(3)=dch/sqrt(3.0)
    zm(3)=-dch/sqrt(3.0)
    xm(4)=-dch/sqrt(3.0)     ! H4
    ym(4)=dch/sqrt(3.0)
    zm(4)=dch/sqrt(3.0)
    xm(5)=dch/sqrt(3.0)      ! H5
    ym(5)=-dch/sqrt(3.0)
    zm(5)=dch/sqrt(3.0)

    xcm=sum(xmm*xm)/sum(xmm) ! Subtract center of mass
    ycm=sum(xmm*ym)/sum(xmm)
    zcm=sum(xmm*zm)/sum(xmm)
    xm=xm-xcm
    ym=ym-ycm
    zm=zm-zcm

! Weights for inverse transformation

    ufact=(/0.0,-1.0,1.0,-1.0,1.0/)
    vfact=(/0.0,-1.0,1.0,1.0,-1.0/)

    return

  end subroutine mgeom

!**********************************************************************!

  subroutine getcf()

!**********************************************************************!

    integer::i,k,l,nran

! Maximum argument for exponential & pi

    alnmax=log(0.1*huge(1.0))
    pi=4.0*atan(1.0)

! RNG characteristics

    call random_seed(size=nran)

    allocate(iran(nran))

! Read startup/checkpoint file

    open(unit=iin,file="mcme_in.dat",status="old",action="read", &
         form="unformatted",position="rewind")

    read(unit=iin) n,nt,ntjob,ntprint,ntskip,iran
    read(unit=iin) da,disp,dr,rho,t

    c=(n/rho)**(1.0/3.0)
    c2=0.5*c
    drm1=1.0/dr
    ndr=int(c2*drm1)

! Allocate arrays

    allocate(ag(1:3,0:ndr-1),umat(n,n),x(m,n),y(m,n),z(m,n), &
         xcm(n),ycm(n),zcm(n))

! Positions

    read(unit=iin) x,y,z

! Averages

    if(nt>0) then
      read(unit=iin) accr,au,au2,aw
      read(unit=iin) ndr
      read(unit=iin) ag
    else
      accr=0.0_double
      au=0.0_double
      au2=0.0_double
      ag=0_long
    end if

    close(unit=iin)

! RNG seed

    call random_seed(put=iran)

! Check positions

    do i=1,n

      xcm(i)=sum(xmm*x(:,i))/sum(xmm)
      ycm(i)=sum(xmm*y(:,i))/sum(xmm)
      zcm(i)=sum(xmm*z(:,i))/sum(xmm)
      if(max(abs(xcm(i)),abs(ycm(i)),abs(zcm(i)))>tolbox*c2) then
        write(unit=*,fmt="(a,i5)") " mcme: out of box # ",i
        stop
      end if

      do k=1,m-1
        do l=k+1,m
          if(abs(sqrt((x(l,i)-x(k,i))**2+(y(l,i)-y(k,i))**2+ &
               (z(l,i)-z(k,i))**2)-sqrt((xm(l)-xm(k))**2+ &
               (ym(l)-ym(k))**2+(zm(l)-zm(k))**2))>tolbond) then
            write(unit=*,fmt="(a,3i5)") &
                 " mcme: bond length error # ",i,k,l
            stop
          end if
        end do
      end do

    end do

! Parameters

    c22=c2**2
    c2m=-c2
    cm=-c
    dcx=maxval(sqrt(xm**2+ym**2+zm**2))

! Cutoff & tail corrections

    rc=c2
    rc2=rc**2
    su0=(2.0*pi*rho*n)*(sum(a12)/(9.0*rc**9)-sum(c6)/(3.0*rc**3))
    sw0=(2.0*pi*rho*n)* &
         (6.0*sum(c6)/(3.0*rc**3)-12.0*sum(a12)/(9.0*rc**9))

! Initialize acceptance counter

    nacc=0

    return

  end subroutine getcf

!**********************************************************************!

  subroutine tables()

!**********************************************************************!

    integer::i
    real::dr2,r

! 1% safety margin

    dr2=(1.01*(rc+2.0*dcx))**2/ntab
    dr2m1=1.0/dr2

    do i=1,ntab
      r=sqrt(i*dr2)
      utab(:,:,i)=qq/r+(a12/r**6-c6)/r**6
      ftab(:,:,i)=-qq/r**3+(6.0*c6-12.0*a12/r**6)/r**8
    end do
    utab(:,:,0)=10.0*utab(:,:,1)     ! Finite values at r=0
    ftab(:,:,0)=10.0*ftab(:,:,1)

    return

  end subroutine tables

!**********************************************************************!

  subroutine pbc(dx,dy,dz,xsh,ysh,zsh)

!**********************************************************************!

    real,intent(in)::dx,dy,dz
    real,intent(out)::xsh,ysh,zsh

    if(dx>=c2) then
      xsh=c
    else if(dx<c2m) then
      xsh=cm
    else
      xsh=0.0
    end if

    if(dy>=c2) then
      ysh=c
    else if(dy<c2m) then
      ysh=cm
    else
      ysh=0.0
    end if

    if(dz>=c2) then
      zsh=c
    else if(dz<c2m) then
      zsh=cm
    else
      zsh=0.0
    end if

    return

  end subroutine pbc

!**********************************************************************!

  subroutine uinit()

!**********************************************************************!

    integer::i,ir2,j,k,l
    real::dx,dy,dz,fact,r2,xki,yki,zki,xsh,ysh,zsh

    umat=0.0

    do i=1,n-1

      do j=i+1,n

        dx=xcm(j)-xcm(i)
        dy=ycm(j)-ycm(i)
        dz=zcm(j)-zcm(i)
        call pbc(dx,dy,dz,xsh,ysh,zsh)                   ! PBC

        if((dx-xsh)**2+(dy-ysh)**2+(dz-zsh)**2<rc2) then ! COM cutoff
          do k=1,m
            xki=x(k,i)+xsh
            yki=y(k,i)+ysh
            zki=z(k,i)+zsh
            do l=1,m
              dx=x(l,j)-xki            ! Site-site interactions
              dy=y(l,j)-yki
              dz=z(l,j)-zki
              r2=dx**2+dy**2+dz**2
              fact=r2*dr2m1
              ir2=int(fact)
              fact=fact-ir2
              umat(i,j)=umat(i,j)+ &
                   ((1.0-fact)*utab(k,l,ir2)+fact*utab(k,l,ir2+1))
            end do
          end do
        end if
        umat(j,i)=umat(i,j)

      end do

    end do

    return

  end subroutine uinit

!**********************************************************************!

  subroutine move()

!**********************************************************************!

    logical::accept
    integer::i,ir2,j,k,l
    real::alpha,cosa,dx,dy,dz,fact,r2,ran,sina,su, &
         ux,uy,uz,vx,vy,vz,wx,wy,wz, &
         xcmn,ycmn,zcmn,xki,yki,zki,xsh,ysh,zsh
    real,dimension(3,3)::rmat
    real,dimension(m)::xn,yn,zn
    real,dimension(n)::ui

! Select particle

    call random_number(ran)
    i=min(int(n*ran+1.0),n)

! Translate

    call random_number(ran)
    xcmn=xcm(i)+disp*(ran-0.5)
    call random_number(ran)
    ycmn=ycm(i)+disp*(ran-0.5)
    call random_number(ran)
    zcmn=zcm(i)+disp*(ran-0.5)
    call pbc(xcmn,ycmn,zcmn,xsh,ysh,zsh)
    xcmn=xcmn-xsh
    ycmn=ycmn-ysh
    zcmn=zcmn-zsh

! Rotate

    do                                     ! Axis
      call random_number(ran)
      ux=ran-0.5
      call random_number(ran)
      uy=ran-0.5
      call random_number(ran)
      uz=ran-0.5
      if(ux**2+uy**2+uz**2<=0.25) then
        exit
      end if
    end do
    fact=1.0/sqrt(ux**2+uy**2+uz**2)
    ux=fact*ux
    uy=fact*uy
    uz=fact*uz

    call random_number(ran)           ! Angle
    alpha=da*(pi/180.0)*(ran-0.5)
    cosa=cos(alpha)
    sina=sin(alpha)

    rmat=reshape((/ &                 ! Rotation matrix
         cosa,sina*uz,-sina*uy, &
         -sina*uz,cosa,sina*ux, &
         sina*uy,-sina*ux,cosa/),(/3,3/))
    rmat=rmat+(1.0-cosa)*matmul(reshape((/ux,uy,uz/),(/3,1/)), &
         reshape((/ux,uy,uz/),(/1,3/)))

    xn=rmat(1,1)*(x(:,i)-xcm(i))+rmat(1,2)*(y(:,i)-ycm(i))+ & ! Rotate
         rmat(1,3)*(z(:,i)-zcm(i))
    yn=rmat(2,1)*(x(:,i)-xcm(i))+rmat(2,2)*(y(:,i)-ycm(i))+ &
         rmat(2,3)*(z(:,i)-zcm(i))
    zn=rmat(3,1)*(x(:,i)-xcm(i))+rmat(3,2)*(y(:,i)-ycm(i))+ &
         rmat(3,3)*(z(:,i)-zcm(i))

    ux=sum(ufact*xn)                      ! Local coordinate axes
    uy=sum(ufact*yn)
    uz=sum(ufact*zn)
    fact=1.0/sqrt(ux**2+uy**2+uz**2)
    ux=fact*ux
    uy=fact*uy
    uz=fact*uz
    vx=sum(vfact*xn)
    vy=sum(vfact*yn)
    vz=sum(vfact*zn)
    fact=ux*vx+uy*vy+uz*vz
    vx=vx-fact*ux
    vy=vy-fact*uy
    vz=vz-fact*uz
    fact=1.0/sqrt(vx**2+vy**2+vz**2)
    vx=fact*vx
    vy=fact*vy
    vz=fact*vz
    wx=uy*vz-uz*vy
    wy=uz*vx-ux*vz
    wz=ux*vy-uy*vx

! Re-assemble molecule

    xn=xcmn+xm*ux+ym*vx+zm*wx
    yn=ycmn+xm*uy+ym*vy+zm*wy
    zn=zcmn+xm*uz+ym*vz+zm*wz

! Energy difference

    su=0.0
    do j=1,n

      ui(j)=0.0

      if(j/=i) then

        dx=xcm(j)-xcmn
        dy=ycm(j)-ycmn
        dz=zcm(j)-zcmn
        call pbc(dx,dy,dz,xsh,ysh,zsh)                   ! PBC

        if((dx-xsh)**2+(dy-ysh)**2+(dz-zsh)**2<rc2) then ! COM cutoff
          do k=1,m
            xki=xn(k)+xsh
            yki=yn(k)+ysh
            zki=zn(k)+zsh
            do l=1,m
              dx=x(l,j)-xki            ! Site-site interactions
              dy=y(l,j)-yki
              dz=z(l,j)-zki
              r2=dx**2+dy**2+dz**2
              fact=r2*dr2m1
              ir2=int(fact)
              fact=fact-ir2
              ui(j)=ui(j)+ &
                   ((1.0-fact)*utab(k,l,ir2)+fact*utab(k,l,ir2+1))
            end do
          end do
        end if

        su=su+(ui(j)-umat(i,j))

      end if

    end do

! Acceptance test

    if(su<=0.0) then
      accept=.true.
    else
      if(su/(rg*t)<alnmax) then
        call random_number(ran)
        if(ran<=exp(-su/(rg*t))) then
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
      x(:,i)=xn
      y(:,i)=yn
      z(:,i)=zn
      xcm(i)=xcmn
      ycm(i)=ycmn
      zcm(i)=zcmn
      umat(:,i)=ui
      umat(i,:)=ui
    end if

    return

  end subroutine move

!**********************************************************************!

  subroutine means()

!**********************************************************************!

    integer,dimension(m,m),parameter::ig=reshape((/ &
         1,2,2,2,2, &
         2,3,3,3,3, &
         2,3,3,3,3, &
         2,3,3,3,3, &
         2,3,3,3,3/),(/m,m/))

    integer::i,ir,ir2,j,k,l
    real::dfx,dfy,dfz,dx,dy,dz,dxcm,dycm,dzcm,fact,fkl,p,r2,su,sw, &
         xki,yki,zki,xsh,ysh,zsh

! Potential energy & virial

    su=su0
    sw=sw0
    do i=1,n-1
      do j=i+1,n
        dxcm=xcm(j)-xcm(i)
        dycm=ycm(j)-ycm(i)
        dzcm=zcm(j)-zcm(i)
        call pbc(dxcm,dycm,dzcm,xsh,ysh,zsh)                    ! PBC
        if((dxcm-xsh)**2+(dycm-ysh)**2+(dzcm-zsh)**2<=rc2) then ! Cutoff
          dfx=0.0
          dfy=0.0
          dfz=0.0
          do k=1,m             ! Site-site
            xki=x(k,i)+xsh
            yki=y(k,i)+ysh
            zki=z(k,i)+zsh
            do l=1,m
              dx=x(l,j)-xki
              dy=y(l,j)-yki
              dz=z(l,j)-zki
              r2=dx**2+dy**2+dz**2
              fact=r2*dr2m1
              ir2=int(fact)
              fact=fact-ir2
              fkl=(1.0-fact)*ftab(k,l,ir2)+fact*ftab(k,l,ir2+1)
              dfx=dfx+fkl*dx
              dfy=dfy+fkl*dy
              dfz=dfz+fkl*dz
            end do
          end do
          su=su+umat(i,j)                        ! Energy
          sw=sw+(dfx*dxcm+dfy*dycm+dfz*dzcm)     ! Virial
        end if
      end do
    end do
    su=su/n
    sw=sw/n

! g_AB(r)

    do i=1,n-1
      do j=i+1,n
        do k=1,m
          do l=1,m
            dx=x(l,j)-x(k,i)
            dy=y(l,j)-y(k,i)
            dz=z(l,j)-z(k,i)
            call pbc(dx,dy,dz,xsh,ysh,zsh)
            r2=(dx-xsh)**2+(dy-ysh)**2+(dz-zsh)**2
            if(r2<c22) then
              ir=int(sqrt(r2)*drm1)
              if(ir<ndr) then
                ag(ig(k,l),ir)=ag(ig(k,l),ir)+1_long
              end if
            end if
          end do
        end do
      end do
    end do

! Print control variables

    if(ntprint>0.and.modulo(nt,ntprint)==0) then
      p=rho*(rg*t-sw/3.0)
      write(unit=*,fmt="(tr1,i10,3(tr1,es12.5))") &
           nt,real(nacc)/(n*ntskip),su,p/atm
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

  subroutine putcf()

!**********************************************************************!

! RNG seed

    call random_seed(get=iran)

! Write checkpoint file

    open(unit=iout,file="mcme_out.dat",status="replace", &
         action="write",form="unformatted")

    write(unit=iout) n,nt,ntjob,ntprint,ntskip,iran
    write(unit=iout) da,disp,dr,rho,t
    write(unit=iout) x,y,z
    write(unit=iout) accr,au,au2,aw
    write(unit=iout) ndr
    write(unit=iout) ag

    close(unit=iout)

! Deallocate arrays

    deallocate(ag,iran,umat,x,y,z,xcm,ycm,zcm)

    return

  end subroutine putcf

end module mcme_subm

!**********************************************************************!

program mcme

!**********************************************************************!

  use mcme_glbm
  use mcme_subm

  implicit none

  integer::i,j

! Read startup/checkpoint file & initialize

  call mgeom()
  call getcf()
  call tables()
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

end program mcme

!**********************************************************************!
