!**********************************************************************!
!
! File: lmdnbutane.f90
!
! Create initial configuration (fcc "plastic crystal") for
! NVT (Nose-Hoover) Molecular Dynamics of OPLS n-butane
!
! Potential parameters and molecular geometry from:
!   W. L. Jorgensen, J. D. Madura, and C. J. Swenson, J. Am. Chem. Soc.
!   106, 6638 (1984)
!
! 19-Apr-2009 (MN)
! 31-May-2012
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

program lmdnbutane

!**********************************************************************!

  implicit none

  integer,parameter::double=selected_real_kind(15)
  integer,parameter::iout=2
  integer,parameter::m=4,nf=7
  real(kind=double),parameter::avog=602.252_double,rg=0.0083143_double
  real(kind=double),parameter::ccc=112.0_double,dcc=0.153_double, &
       xmch2=14.0_double,xmch3=15.0_double
  real(kind=double),dimension(m),parameter:: &
       xms=(/xmch3,xmch2,xmch2,xmch3/)

  integer::i,ix,iy,iz,nc
  integer::n,nt,ntjob,ntprint,ntskip
  real(kind=double)::c,cosphi,fact,phi,pi,ran,sinphi,sk, &
       ux,uy,uz,xc,yc,zc
  real(kind=double)::dphi,dr,dt,t,tau,vol,zeta
  real(kind=double),dimension(m)::xm,ym,zm
  real(kind=double),dimension(3,3)::rmat
  real(kind=double),dimension(:,:),allocatable::x,y,z,xl,yl,zl

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

! User input

! Fcc lattice

  write(unit=*,fmt="(a)",advance="no") "              n="
  read(unit=*,fmt=*) n
  nc=0
  do
    nc=nc+1
    if(nc**3==2*n) then
      exit
    else if(nc**3>2*n) then
      write(unit=*,fmt="(a)") " lmdnbutane: not a magic number"
      stop
    end if
  end do

! Allocate arrays

  allocate(x(m,n),y(m,n),z(m,n),xl(m,n),yl(m,n),zl(m,n))

  write(unit=*,fmt="(a)",advance="no") "           t(K)="
  read(unit=*,fmt=*) t
  write(unit=*,fmt="(a)",advance="no") "    vmol(cm**3)="
  read(unit=*,fmt=*) vol
  vol=n*vol/avog
  write(unit=*,fmt="(a)",advance="no") "        tau(ps)="
  read(unit=*,fmt=*) tau
  write(unit=*,fmt="(a)",advance="no") "         dt(ps)="
  read(unit=*,fmt=*) dt
  write(unit=*,fmt="(a)",advance="no") "         dr(nm)="
  read(unit=*,fmt=*) dr
  write(unit=*,fmt="(a)",advance="no") "      dphi(deg)="
  read(unit=*,fmt=*) dphi
  write(unit=*,fmt="(a)",advance="no") "         ntskip="
  read(unit=*,fmt=*) ntskip
  write(unit=*,fmt="(a)",advance="no") " ntprint/ntskip="
  read(unit=*,fmt=*) ntprint
  write(unit=*,fmt="(a)",advance="no") "   ntjob/ntskip="
  read(unit=*,fmt=*) ntjob

! Fcc positions, random orientations & COM velocities
! All-"trans" conformations

  c=vol**(1.0_double/3.0_double)

  i=0

  do ix=0,nc-1
    do iy=0,nc-1
      do iz=0,nc-1
        if(modulo(ix+iy+iz,2)==0) then
          i=i+1

! Random axis

          do
            call random_number(ran)
            ux=ran-0.5_double
            call random_number(ran)
            uy=ran-0.5_double
            call random_number(ran)
            uz=ran-0.5_double
            if(ux**2+uy**2+uz**2<=0.25_double) then
              exit
            end if
          end do
          fact=1.0_double/sqrt(ux**2+uy**2+uz**2)
          ux=fact*ux
          uy=fact*uy
          uz=fact*uz

! Random angle

          call random_number(ran)
          phi=2.0_double*pi*ran
          cosphi=cos(phi)
          sinphi=sin(phi)

          rmat(1,:)=(/cosphi,-sinphi*uz,sinphi*uy/)
          rmat(2,:)=(/sinphi*uz,cosphi,-sinphi*ux/)
          rmat(3,:)=(/-sinphi*uy,sinphi*ux,cosphi/)
          rmat=rmat+(1.0-cosphi)*matmul(reshape((/ux,uy,uz/),(/3,1/)), &
               reshape((/ux,uy,uz/),(/1,3/)))

! Fcc COM position
          
          xc=((ix+0.5_double)/nc-0.5_double)*c
          yc=((iy+0.5_double)/nc-0.5_double)*c
          zc=((iz+0.5_double)/nc-0.5_double)*c
          x(:,i)=xc+(rmat(1,1)*xm(:)+rmat(1,2)*ym(:)+rmat(1,3)*zm(:))
          y(:,i)=yc+(rmat(2,1)*xm(:)+rmat(2,2)*ym(:)+rmat(2,3)*zm(:))
          z(:,i)=zc+(rmat(3,1)*xm(:)+rmat(3,2)*ym(:)+rmat(3,3)*zm(:))

! Random COM velocity

          call random_number(ran)
          xl(:,i)=x(:,i)-(ran-0.5_double)
          call random_number(ran)
          yl(:,i)=y(:,i)-(ran-0.5_double)
          call random_number(ran)
          zl(:,i)=z(:,i)-(ran-0.5_double)

        end if
      end do
    end do
  end do

! Zero net momentum & set initial kinetic energy

  xl=xl+sum(spread(xms,2,n)*(x-xl))/(n*sum(xms))
  yl=yl+sum(spread(xms,2,n)*(y-yl))/(n*sum(xms))
  zl=zl+sum(spread(xms,2,n)*(z-zl))/(n*sum(xms))

  sk=sum(spread(xms,2,n)*((x-xl)**2+(y-yl)**2+(z-zl)**2))/dt**2

  fact=sqrt(nf*n*rg*t/sk)

  xl=x-fact*(x-xl)
  yl=y-fact*(y-yl)
  zl=z-fact*(z-zl)

! Write startup file

  nt=0
  zeta=0.0_double

  open(unit=iout,file="mdnbutane_in.dat",status="replace", &
       action="write",form="unformatted",position="rewind")

  write(unit=iout) n,nt,ntjob,ntprint,ntskip
  write(unit=iout) dphi,dr,dt,t,vol,tau,zeta
  write(unit=iout) x,y,z
  write(unit=iout) xl,yl,zl

  close(unit=iout)

! Deallocate arrays

  deallocate(x,y,z,xl,yl,zl)

end program lmdnbutane

!**********************************************************************!
