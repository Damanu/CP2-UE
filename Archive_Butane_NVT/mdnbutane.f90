!**********************************************************************!
!
! File: mdnbutane.f90
!
! NVT (Nose-Hoover) Molecular Dynamics of OPLS n-butane
! Simple Verlet & SHAKE algorithm
!
! Potential parameters and molecular geometry from:
!   W. L. Jorgensen, J. D. Madura, and C. J. Swenson, J. Am. Chem. Soc.
!   106, 6638 (1984)
!
! 19-Apr-2009 (MN)
! 01-Jun-2012
!
!**********************************************************************!
!
!                           CH3     Site 4
!                           /
!                          /
!                         /
!                        /
!                       /
!         CH2---------CH2           Sites 2 & 3
!         /
!        /
!       /
!      /
!     /
!   CH3                             Site 1
!
!**********************************************************************!

module mdnbutane_glbm

!**********************************************************************!

! Parameters & global variables/arrays

  implicit none

  private

  integer,parameter,public::long=selected_int_kind(18)
  integer,parameter,public::double=selected_real_kind(15)
  integer,parameter,public::iin=1,iout=2
  integer,parameter,public::itmax=200,m=4,nb=5,nf=7
  real(kind=double),parameter,public::avog=602.252_double, &
       atm=avog*0.000101325_double,cal=4.184_double,rg=0.0083143_double
  real(kind=double),parameter,public::ccc=112.0_double, &
       dcc=0.153_double,epsch2=0.118_double*cal, &
       epsch3=0.175_double*cal,sig=0.3905_double, &
       xmch2=14.0_double,xmch3=15.0_double
  real(kind=double),parameter,public::u0=2.0495_double*cal, &
       u1=-4.0495_double*cal,u2=0.315_double*cal, &
       u3=6.414_double*cal
  real(kind=double),parameter,public::tolb2=1.0e-10_double, &
       tolc=1.0e-12_double

  integer,public::ndphi,ndr
  integer,public::n,nt,ntjob,ntprint,ntskip
  integer,dimension(2,nb),public::kb
  integer(kind=long),dimension(:),allocatable,public::agcc,agss,aphi
  real(kind=double),public::c,c2,c2m1,dcx,drm1,pi,r2max,rc,rc2,rho, &
       rt,rt2,sk,su,su0,sud,sw,sw0
  real(kind=double),public::dphi,dr,dt,t,vol,tau,zeta
  real(kind=double),public::ak,ae,ae2,au,aud,aw
  real(kind=double),dimension(nb),public::d2b
  real(kind=double),dimension(m),public::xms
  real(kind=double),dimension(:),allocatable,public::xc,yc,zc
  real(kind=double),dimension(m,m),public::a12,c6
  real(kind=double),dimension(:,:),allocatable,public::fx,fy,fz, &
       x,y,z,xl,yl,zl,xn,yn,zn

end module mdnbutane_glbm

!**********************************************************************!

module mdnbutane_subm

!**********************************************************************!

  use mdnbutane_glbm

  implicit none

  private

  public::fintra,force,getcf,means,mgeom,move,putcf,shake

contains

!**********************************************************************!

  subroutine mgeom()

!**********************************************************************!

    real(kind=double),dimension(m)::xm,ym,zm
    
! Pi

    pi=4.0_double*atan(1.0_double)

! Molecular geometry

    xm(1)=(cos(ccc*pi/180.0_double)-0.5_double)*dcc
    ym(1)=-sin(ccc*pi/180.0_double)*dcc
    zm(1)=0.0_double
    xm(2)=-0.5_double*dcc
    ym(2)=0.0_double
    zm(2)=0.0_double
    xm(3)=0.5_double*dcc
    ym(3)=0.0_double
    zm(3)=0.0_double
    xm(4)=(0.5_double-cos(ccc*pi/180.0_double))*dcc
    ym(4)=sin(ccc*pi/180.0_double)*dcc
    zm(4)=0.0_double

! Site masses

    xms=(/xmch3,xmch2,xmch2,xmch3/)

! Maximum center-site distance

    dcx=maxval(sqrt(xm**2+ym**2+zm**2))

! Constraints

    kb(1,:)=(/1,1,2,2,3/)
    kb(2,:)=(/2,3,3,4,4/)

    d2b((/1,3,5/))=dcc**2
    d2b((/2,4/))=2.0_double*dcc**2*(1.0_double-cos(ccc*pi/180.0_double))

! Potential parameters

    a12=4.0_double*sqrt(spread((/epsch3,epsch2,epsch2,epsch3/),2,m)* &
         spread((/epsch3,epsch2,epsch2,epsch3/),1,m))*sig**12
    c6=4.0_double*sqrt(spread((/epsch3,epsch2,epsch2,epsch3/),2,m)* &
         spread((/epsch3,epsch2,epsch2,epsch3/),1,m))*sig**6

    return

  end subroutine mgeom

!**********************************************************************!

  subroutine getcf()

!**********************************************************************!

    integer::i,ib,k

! Read startup/checkpoint file

    open(unit=iin,file="mdnbutane_in.dat",status="old",action="read", &
         form="unformatted",position="rewind")

    read(unit=iin) n,nt,ntjob,ntprint,ntskip
    read(unit=iin) dphi,dr,dt,t,vol,tau,zeta

! Allocate arrays

    ndphi=int(180.0_double/dphi)
    ndr=int(0.5_double*vol**(1.0_double/3.0_double)/dr)

    allocate(agcc(0:ndr-1),agss(0:ndr-1),aphi(0:ndphi-1), &
         fx(m,n),fy(m,n),fz(m,n),xc(n),yc(n),zc(n), &
         x(m,n),y(m,n),z(m,n),xl(m,n),yl(m,n),zl(m,n), &
         xn(m,n),yn(m,n),zn(m,n))

! Positions

    read(unit=iin) x,y,z
    read(unit=iin) xl,yl,zl

! Averages

    if(nt>0) then
      read(unit=iin) ak,ae,ae2,au,aud,aw
      read(unit=iin) agcc
      read(unit=iin) agss
      read(unit=iin) aphi
    else
      ak=0.0_double
      ae=0.0_double
      ae2=0.0_double
      au=0.0_double
      aud=0.0_double
      aw=0.0_double
      agcc=0_long
      agss=0_long
      aphi=0_long
    end if

    close(unit=iin)

! Box parameters

    c=vol**(1.0_double/3.0_double)
    c2=0.5_double*c
    c2m1=2.0_double/c
    drm1=1.0_double/dr
    r2max=c2**2
    rho=n/vol

! Check molecule centers

    do i=1,n
      xc(i)=0.5_double*(x(2,i)+x(3,i))
      yc(i)=0.5_double*(y(2,i)+y(3,i))
      zc(i)=0.5_double*(z(2,i)+z(3,i))
      if(max(abs(xc(i)),abs(yc(i)),abs(zc(i)))-c2>tolc*c2) then
        write(unit=*,fmt=*) "mdnbutane: out of box",i,c2
        do k=1,m
          write(unit=*,fmt=*) k,x(k,i),y(k,i),z(k,i)
        end do
        stop
      end if
    end do

! Check constraints

    do i=1,n
      do ib=1,nb
        if(abs((x(kb(2,ib),i)-x(kb(1,ib),i))**2+ &
             (y(kb(2,ib),i)-y(kb(1,ib),i))**2+ &
             (z(kb(2,ib),i)-z(kb(1,ib),i))**2-d2b(ib))> &
             tolb2*d2b(ib)) then
          write(unit=*,fmt=*) "mdnbutane: bond length error", &
               i,ib,kb(1,ib),kb(2,ib)
          write(unit=*,fmt=*) sqrt((x(kb(2,ib),i)-x(kb(1,ib),i))**2+ &
               (y(kb(2,ib),i)-y(kb(1,ib),i))**2+ &
               (z(kb(2,ib),i)-z(kb(1,ib),i))**2),sqrt(d2b(ib))
          stop
        end if
      end do
    end do

! Cutoff & tail corrections

    rc=c2
    rc2=rc**2
    rt=rc-2.0_double*dcx
    rt2=rt**2

    su0=2.0_double*pi*n*(n/vol)* &
         (sum(a12)/(9.0_double*rt**9)-sum(c6)/(3.0_double*rt**3))
    sw0=2.0_double*pi*n*(n/vol)* &
         (6.0_double*sum(c6)/(3.0_double*rt**3)- &
         12.0_double*sum(a12)/(9.0_double*rt**9))

    return
    
  end subroutine getcf

!**********************************************************************!

  subroutine fintra()

!**********************************************************************!

    integer::i
    real(kind=double)::x12,y12,z12,x13,y13,z13,x23,y23,z23, &
         x24,y24,z24,x34,y34,z34
    real(kind=double)::cosphi,d123,d234,df,x123,y123,z123,x234,y234,z234
    real(kind=double)::x122,y122,z122,x132,y132,z132,x133,y133,z133, &
         x232,y232,z232,x233,y233,z233,x242,y242,z242,x243,y243,z243, &
         x342,y342,z342,x343,y343,z343

! Intramolecular (dihedral) forces & potential energy
! Called by subroutine "force"

    sud=0.0_double

    do i=1,n

! Bond vectors

      x12=x(2,i)-x(1,i)
      y12=y(2,i)-y(1,i)
      z12=z(2,i)-z(1,i)
      x13=x(3,i)-x(1,i)
      y13=y(3,i)-y(1,i)
      z13=z(3,i)-z(1,i)
      x23=x(3,i)-x(2,i)
      y23=y(3,i)-y(2,i)
      z23=z(3,i)-z(2,i)
      x24=x(4,i)-x(2,i)
      y24=y(4,i)-y(2,i)
      z24=z(4,i)-z(2,i)
      x34=x(4,i)-x(3,i)
      y34=y(4,i)-y(3,i)
      z34=z(4,i)-z(3,i)

! Normal vectors (overwritten)

      x123=y12*z23-z12*y23
      y123=z12*x23-x12*z23
      z123=x12*y23-y12*x23

      d123=sqrt(x123**2+y123**2+z123**2)

      x234=y23*z34-z23*y34
      y234=z23*x34-x23*z34
      z234=x23*y34-y23*x34

      d234=sqrt(x234**2+y234**2+z234**2)

! Dihedral angle

      cosphi=(x123*x234+y123*y234+z123*z234)/(d123*d234)

! Potential energy

      sud=sud+(u0+cosphi*(u1+cosphi*(u2+cosphi*u3)))

! Auxiliary vectors

      x122=y12*z123-z12*y123
      y122=z12*x123-x12*z123
      z122=x12*y123-y12*x123

      x132=y13*z123-z13*y123
      y132=z13*x123-x13*z123
      z132=x13*y123-y13*x123
      x133=y13*z234-z13*y234
      y133=z13*x234-x13*z234
      z133=x13*y234-y13*x234

      x232=y23*z123-z23*y123
      y232=z23*x123-x23*z123
      z232=x23*y123-y23*x123
      x233=y23*z234-z23*y234
      y233=z23*x234-x23*z234
      z233=x23*y234-y23*x234

      x242=y24*z123-z24*y123
      y242=z24*x123-x24*z123
      z242=x24*y123-y24*x123
      x243=y24*z234-z24*y234
      y243=z24*x234-x24*z234
      z243=x24*y234-y24*x234

      x342=y34*z123-z34*y123
      y342=z34*x123-x34*z123
      z342=x34*y123-y34*x123
      x343=y34*z234-z34*y234
      y343=z34*x234-x34*z234
      z343=x34*y234-y34*x234

      x123=y12*z234-z12*y234
      y123=z12*x234-x12*z234
      z123=x12*y234-y12*x234

! Distribute forces

      df=-(u1+2.0_double*u2*cosphi+3.0_double*u3*cosphi**2)

      fx(1,i)=fx(1,i)+df*(cosphi*x232/d123-x233/d234)/d123
      fy(1,i)=fy(1,i)+df*(cosphi*y232/d123-y233/d234)/d123
      fz(1,i)=fz(1,i)+df*(cosphi*z232/d123-z233/d234)/d123

      fx(2,i)=fx(2,i)+df*((x133-x342)/(d123*d234) &
           -cosphi*(x132/d123**2-x343/d234**2))
      fy(2,i)=fy(2,i)+df*((y133-y342)/(d123*d234) &
           -cosphi*(y132/d123**2-y343/d234**2))
      fz(2,i)=fz(2,i)+df*((z133-z342)/(d123*d234) &
           -cosphi*(z132/d123**2-z343/d234**2))

      fx(3,i)=fx(3,i)+df*((x242-x123)/(d123*d234) &
           -cosphi*(x243/d234**2-x122/d123**2))
      fy(3,i)=fy(3,i)+df*((y242-y123)/(d123*d234) &
           -cosphi*(y243/d234**2-y122/d123**2))
      fz(3,i)=fz(3,i)+df*((z242-z123)/(d123*d234) &
           -cosphi*(z243/d234**2-z122/d123**2))

      fx(4,i)=fx(4,i)+df*(cosphi*x233/d234-x232/d123)/d234
      fy(4,i)=fy(4,i)+df*(cosphi*y233/d234-y232/d123)/d234
      fz(4,i)=fz(4,i)+df*(cosphi*z233/d234-z232/d123)/d234

    end do

    return

  end subroutine fintra

!**********************************************************************!

  subroutine force()

!**********************************************************************!

    integer::i,j,k,l
    real(kind=double)::df,dfx,dfy,dfz,dfxc,dfyc,dfzc,dx,dy,dz, &
         dxc,dyc,dzc,r2,rm2,rm6,xki,yki,zki,xsh,ysh,zsh

! Initialize forces & potential energy/virial

    fx=0.0_double
    fy=0.0_double
    fz=0.0_double

    su=su0
    sw=sw0

! Intramolcular forces & potential energy

    call fintra()

! Intermolecular forces, potential energy & virial

    do i=1,n-1
      do j=i+1,n

! Apply periodic boundary conditions to molecule centers

        dxc=xc(j)-xc(i)
        xsh=int(dxc*c2m1)*c
        dxc=dxc-xsh
        dyc=yc(j)-yc(i)
        ysh=int(dyc*c2m1)*c
        dyc=dyc-ysh
        dzc=zc(j)-zc(i)
        zsh=int(dzc*c2m1)*c
        dzc=dzc-zsh
        if(dxc**2+dyc**2+dzc**2<rc2) then ! Center-center cutoff

          dfxc=0.0_double
          dfyc=0.0_double
          dfzc=0.0_double

! Site-site forces

          do k=1,m
            xki=x(k,i)+xsh
            yki=y(k,i)+ysh
            zki=z(k,i)+zsh
            do l=1,m
              dx=x(l,j)-xki
              dy=y(l,j)-yki
              dz=z(l,j)-zki
              r2=dx**2+dy**2+dz**2
              if(r2<rt2) then      ! Site-site cutoff
                rm2=1.0_double/r2
                rm6=rm2**3
                su=su+(a12(k,l)*rm6-c6(k,l))*rm6
                df=(6.0_double*c6(k,l)-12.0_double*a12(k,l)*rm6)*rm6*rm2
                dfx=df*dx
                dfy=df*dy
                dfz=df*dz
                fx(k,i)=fx(k,i)+dfx
                fx(l,j)=fx(l,j)-dfx
                fy(k,i)=fy(k,i)+dfy
                fy(l,j)=fy(l,j)-dfy
                fz(k,i)=fz(k,i)+dfz
                fz(l,j)=fz(l,j)-dfz
                dfxc=dfxc+dfx
                dfyc=dfyc+dfy
                dfzc=dfzc+dfz
              end if
            end do
          end do

! Virial

          sw=sw+(dxc*dfxc+dyc*dfyc+dzc*dfzc)

        end if

      end do
    end do

    return

  end subroutine force

!**********************************************************************!

  subroutine shake()

!**********************************************************************!

    logical::converged
    integer::i,ib,it
    real(kind=double),dimension(nb)::dx,dy,dz,dxn,dyn,dzn,r2,xlambda

! Iterative correction of new positions (constraint forces)
! Called by subroutine "move"

    do i=1,n

! Old bond vectors
      
      do ib=1,nb
        dx(ib)=x(kb(2,ib),i)-x(kb(1,ib),i)
        dy(ib)=y(kb(2,ib),i)-y(kb(1,ib),i)
        dz(ib)=z(kb(2,ib),i)-z(kb(1,ib),i)
        r2(ib)=dx(ib)**2+dy(ib)**2+dz(ib)**2
      end do

      xlambda=0.0_double
      converged=.false.

      do it=1,itmax

! Shake: treat constraints one by one

        do ib=1,nb

! New bond vectors
          
          dxn(ib)=xn(kb(2,ib),i)-xn(kb(1,ib),i)
          dyn(ib)=yn(kb(2,ib),i)-yn(kb(1,ib),i)
          dzn(ib)=zn(kb(2,ib),i)-zn(kb(1,ib),i)

! Lagrange multipliers

          xlambda(ib)=(dxn(ib)**2+dyn(ib)**2+dzn(ib)**2-d2b(ib)+ &
               (0.5_double*(dt**2/xms(kb(1,ib))+dt**2/xms(kb(2,ib)))* &
               xlambda(ib))**2*r2(ib))/ &
               ((dt**2/xms(kb(1,ib))+dt**2/xms(kb(2,ib)))* &
               (dxn(ib)*dx(ib)+dyn(ib)*dy(ib)+dzn(ib)*dz(ib)))

! Update new positions

          xn(kb(1,ib),i)=xn(kb(1,ib),i)+ &
               0.5_double*dt**2/xms(kb(1,ib))*xlambda(ib)*dx(ib)
          yn(kb(1,ib),i)=yn(kb(1,ib),i)+ &
               0.5_double*dt**2/xms(kb(1,ib))*xlambda(ib)*dy(ib)
          zn(kb(1,ib),i)=zn(kb(1,ib),i)+ &
               0.5_double*dt**2/xms(kb(1,ib))*xlambda(ib)*dz(ib)
          xn(kb(2,ib),i)=xn(kb(2,ib),i)- &
               0.5_double*dt**2/xms(kb(2,ib))*xlambda(ib)*dx(ib)
          yn(kb(2,ib),i)=yn(kb(2,ib),i)- &
               0.5_double*dt**2/xms(kb(2,ib))*xlambda(ib)*dy(ib)
          zn(kb(2,ib),i)=zn(kb(2,ib),i)- &
               0.5_double*dt**2/xms(kb(2,ib))*xlambda(ib)*dz(ib)

        end do

! Check for convergence

        if(all(abs(dxn**2+dyn**2+dzn**2-d2b)<tolb2*d2b)) then
          converged=.true.
          exit
        end if

      end do

      if(.not.converged) then
        write(unit=*,fmt=*) "mdnbutane: shake",nt,i
        do ib=1,nb
          write(unit=*,fmt=*) ib,kb(1,ib),kb(2,ib), &
               sqrt((x(kb(2,ib),i)-x(kb(1,ib),i))**2+ &
               (y(kb(2,ib),i)-y(kb(1,ib),i))**2+ &
               (z(kb(2,ib),i)-z(kb(1,ib),i))**2),sqrt(d2b(ib))
        end do
        stop
      end if

    end do

    return

  end subroutine shake

!**********************************************************************!

  subroutine move()

!**********************************************************************!

    integer::i
    real(kind=double)::fl,fn,xsh,ysh,zsh

! Advance positions
! Verlet & simple Nose-Hoover thermostat

    fl=1.0_double-0.5_double*zeta*dt
    fn=1.0_double/(1.0_double+0.5_double*zeta*dt)

    xn=fn*(2.0_double*x-fl*xl+spread(dt**2/xms,2,n)*fx)
    yn=fn*(2.0_double*y-fl*yl+spread(dt**2/xms,2,n)*fy)
    zn=fn*(2.0_double*z-fl*zl+spread(dt**2/xms,2,n)*fz)

! Correct positions

    call shake()

! Kinetic energy(*2)

    sk=sum(spread(xms,2,n)*((xn-x)**2+(yn-y)**2+(zn-z)**2))/dt**2

! Update friction constant

    zeta=zeta+dt/tau**2*(sk/(nf*n*rg*t)-1.0_double)

! Apply periodic boundary conditions & update positions

    do i=1,n

      xc(i)=0.5_double*(xn(2,i)+xn(3,i))
      yc(i)=0.5_double*(yn(2,i)+yn(3,i))
      zc(i)=0.5_double*(zn(2,i)+zn(3,i))
      xsh=int(xc(i)*c2m1)*c
      ysh=int(yc(i)*c2m1)*c
      zsh=int(zc(i)*c2m1)*c

      xc(i)=xc(i)-xsh
      yc(i)=yc(i)-ysh
      zc(i)=zc(i)-zsh
      xl(:,i)=x(:,i)-xsh
      yl(:,i)=y(:,i)-ysh
      zl(:,i)=z(:,i)-zsh
      x(:,i)=xn(:,i)-xsh
      y(:,i)=yn(:,i)-ysh
      z(:,i)=zn(:,i)-zsh

    end do

    return

  end subroutine move

!**********************************************************************!

  subroutine means()

!**********************************************************************!

    integer::i,iphi,ir,j,k,l
    real(kind=double)::cosphi,dx,dy,dz,r2
    real(kind=double)::x12,y12,z12,x23,y23,z23,x34,y34,z34
    real(kind=double)::d123,d234,x123,y123,z123,x234,y234,z234

! Kinetic energy

    sk=0.5_double*sk

! Dihedral angle PDF

    do i=1,n

! Bond vectors

      x12=x(2,i)-x(1,i)
      y12=y(2,i)-y(1,i)
      z12=z(2,i)-z(1,i)
      x23=x(3,i)-x(2,i)
      y23=y(3,i)-y(2,i)
      z23=z(3,i)-z(2,i)
      x34=x(4,i)-x(3,i)
      y34=y(4,i)-y(3,i)
      z34=z(4,i)-z(3,i)

! Normal vectors

      x123=y12*z23-z12*y23
      y123=z12*x23-x12*z23
      z123=x12*y23-y12*x23

      d123=sqrt(x123**2+y123**2+z123**2)

      x234=y23*z34-z23*y34
      y234=z23*x34-x23*z34
      z234=x23*y34-y23*x34

      d234=sqrt(x234**2+y234**2+z234**2)

      cosphi=(x123*x234+y123*y234+z123*z234)/(d123*d234)
      if(cosphi<-1.0_double) then
        cosphi=-1.0_double
      else if(cosphi>1.0_double) then
        cosphi=1.0_double
      end if

      iphi=min(int((180.0_double*acos(cosphi)/pi)/dphi),ndphi-1)
      aphi(iphi)=aphi(iphi)+1_long

    end do

! Pair correlation functions

    do i=1,n-1
      do j=i+1,n

! Center-center

        dx=xc(j)-xc(i)
        dx=dx-int(dx*c2m1)*c
        dy=yc(j)-yc(i)
        dy=dy-int(dy*c2m1)*c
        dz=zc(j)-zc(i)
        dz=dz-int(dz*c2m1)*c
        r2=dx**2+dy**2+dz**2
        if(r2<r2max) then
          ir=int(sqrt(r2)*drm1)
          if(ir<ndr) then
            agcc(ir)=agcc(ir)+1_long
          end if
        end if

! Site-site

        do k=1,m
          do l=1,m
            dx=x(l,j)-x(k,i)
            dx=dx-int(dx*c2m1)*c
            dy=y(l,j)-y(k,i)
            dy=dy-int(dy*c2m1)*c
            dz=z(l,j)-z(k,i)
            dz=dz-int(dz*c2m1)*c
            r2=dx**2+dy**2+dz**2
            if(r2<r2max) then
              ir=int(sqrt(r2)*drm1)
              if(ir<ndr) then
                agss(ir)=agss(ir)+1_long
              end if
            end if
          end do
        end do

      end do
    end do

! Print control variables

    if(ntprint>0.and.modulo(nt/ntskip,ntprint)==0) then
      write(unit=*,fmt=*) nt,real(zeta), &
           real(sk/(0.5_double*nf*n*rg)),real(su/n), &
           real(sud/n),real((sk+su+sud)/n), &
           real(rho/3.0_double*(sk-sw)/n/atm)
    end if

! Accumulate averages

    ak=ak+sk/n
    au=au+su/n
    aud=aud+sud/n
    ae=ae+(sk+su+sud)/n
    ae2=ae2+((sk+su+sud)/n)**2
    aw=aw+sw/n

    return

  end subroutine means

!**********************************************************************!

  subroutine putcf()

!**********************************************************************!

! Write checkpoint file

    open(unit=iout,file="mdnbutane_out.dat",status="replace", &
         action="write",form="unformatted")

    write(unit=iout) n,nt,ntjob,ntprint,ntskip
    write(unit=iout) dphi,dr,dt,t,vol,tau,zeta
    write(unit=iout) x,y,z
    write(unit=iout) xl,yl,zl
    write(unit=iout) ak,ae,ae2,au,aud,aw
    write(unit=iout) agcc
    write(unit=iout) agss
    write(unit=iout) aphi

    close(unit=iout)

! Deallocate arrays

    deallocate(agcc,agss,aphi,fx,fy,fz,xc,yc,zc,x,y,z,xl,yl,zl, &
         xn,yn,zn)

    return

  end subroutine putcf

end module mdnbutane_subm

!**********************************************************************!

program mdnbutane

!**********************************************************************!

  use mdnbutane_glbm
  use mdnbutane_subm

  implicit none

  integer::i

! Define molecular geometry & read startup/checkpoint file

  call mgeom()

  call getcf()

! Do ntjob time steps

  do i=1,ntjob
    nt=nt+1

    call force()     ! Inter- and intramolecular forces

    call move()      ! Advance positions

    if(modulo(nt,ntskip)==0) then
      call means()   ! Calculate averages every ntskip time steps
    end if

  end do

! Write checkpoint file

  call putcf()

end program mdnbutane

!**********************************************************************!
