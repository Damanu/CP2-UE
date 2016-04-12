!**********************************************************************!
!
! File: lmcme.f90
!
! Create initial configuration (fcc "plastic crystal") for 
! NVT-Monte Carlo of OPLS-AA Methane
!
! Potential parameters and molecular geometry from:
!   G. Kaminski, E. M. Duffy, T. Matsui & W. L. Jorgensen,
!   J. Phys. Chem. 98, 13077 (1994)
!
! 23-Apr-1999 (MN)
! 18-Apr-2012
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

program lmcme

!**********************************************************************!

  implicit none

  integer,parameter::iout=2,m=5
  real,parameter::dch=0.109,xmc=12.0,xmh=1.0
  real,dimension(m),parameter::xmm=(/xmc,xmh,xmh,xmh,xmh/)

  character(len=80)::fname
  integer::i,ic,ix,iy,iz,n,nran,nt,ntjob,ntprint,ntskip
  integer,dimension(:),allocatable::iran
  real::alpha,c,cosa,da,disp,dr,fact,pi,ran,rho,sina,t, &
       ux,uy,uz,xcm,ycm,zcm
  real,dimension(3,3)::rmat
  real,dimension(m)::xm,ym,zm
  real,dimension(:,:),allocatable::x,y,z

  pi=4.0*atan(1.0)

! User input

  magic: do
    write(unit=*,fmt="(a)",advance="no") "              n="
    read(unit=*,fmt=*) n
    ic=0
    do
      ic=ic+2
      if(ic**3==2*n) then
        exit magic
      else if(ic**3>2*n) then
        exit
      end if
    end do
  end do magic

  write(unit=*,fmt="(a)",advance="no") "   rho(1/nm**3)="
  read(unit=*,fmt=*) rho
  write(unit=*,fmt="(a)",advance="no") "           t(K)="
  read(unit=*,fmt=*) t
  write(unit=*,fmt="(a)",advance="no") "       disp(nm)="
  read(unit=*,fmt=*) disp
  write(unit=*,fmt="(a)",advance="no") "        da(deg)="
  read(unit=*,fmt=*) da
  write(unit=*,fmt="(a)",advance="no") "         dr(nm)="
  read(unit=*,fmt=*) dr
  write(unit=*,fmt="(a)",advance="no") "         ntskip="
  read(unit=*,fmt=*) ntskip
  write(unit=*,fmt="(a)",advance="no") " ntprint/ntskip="
  read(unit=*,fmt=*) ntprint
  write(unit=*,fmt="(a)",advance="no") "   ntjob/ntskip="
  read(unit=*,fmt=*) ntjob
  write(unit=*,fmt="(a)",advance="no") &
       "          fname=[mcme_in.dat] "
  read(unit=*,fmt="(a)") fname
  if(fname=="") then
    fname="mcme_in.dat"
  end if

! Allocate arrays

  call random_seed(size=nran)

  allocate(iran(nran),x(m,n),y(m,n),z(m,n))

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

! Create fcc lattice with random orientations

  c=(n/rho)**(1.0/3.0)

  i=0
  do ix=0,ic-1
    do iy=0,ic-1
      do iz=0,ic-1

        if(modulo(ix+iy+iz,2)==0) then
          i=i+1

          xcm=c*((ix+0.5)/ic-0.5)     ! Center of mass
          ycm=c*((iy+0.5)/ic-0.5)
          zcm=c*((iz+0.5)/ic-0.5)

          do                            ! Rotation axis
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

          call random_number(ran)     ! Rotation angle
          alpha=2.0*pi*ran
          cosa=cos(alpha)
          sina=sin(alpha)
          
          rmat=reshape((/ &                 ! Rotation matrix
               cosa,sina*uz,-sina*uy, &
               -sina*uz,cosa,sina*ux, &
               sina*uy,-sina*ux,cosa/),(/3,3/))
          rmat=rmat+(1.0-cosa)*matmul(reshape((/ux,uy,uz/),(/3,1/)), &
               reshape((/ux,uy,uz/),(/1,3/)))
          
          x(:,i)=xcm+(rmat(1,1)*xm+rmat(1,2)*ym+rmat(1,3)*zm) ! Rotate
          y(:,i)=ycm+(rmat(2,1)*xm+rmat(2,2)*ym+rmat(2,3)*zm)
          z(:,i)=zcm+(rmat(3,1)*xm+rmat(3,2)*ym+rmat(3,3)*zm)
          
        end if

      end do
    end do
  end do

! RNG seed

  call random_seed(get=iran)

! Write startup file

  nt=0

  open(unit=iout,file=fname,status="replace",action="write", &
       form="unformatted")

  write(unit=iout) n,nt,ntjob,ntprint,ntskip,iran
  write(unit=iout) da,disp,dr,rho,t
  write(unit=iout) x,y,z

  close(unit=iout)

! Deallocate arrays

  deallocate(iran,x,y,z)

end program lmcme

!**********************************************************************!
