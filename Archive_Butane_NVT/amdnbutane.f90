!**********************************************************************!
!
! File: amdnbutane.f90
!
! NVT (Nose--Hoover) Molecular Dynamics of OPLS n-butane
! Analyze checkpoint file
!
! 03-May-2009 (MN)
! 31-May-2012
!
!**********************************************************************!

program amdnbutane

!**********************************************************************!

  implicit none

  integer,parameter::long=selected_int_kind(18)
  integer,parameter::double=selected_real_kind(15)
  integer,parameter::iin=1,iout=2
  integer,parameter::m=4,nf=7
  integer,parameter::ndphi_id=1000
  real(kind=double),parameter::avog=602.252_double, &
       atm=avog*0.000101325_double,cal=4.184_double,rg=0.0083143_double
  real(kind=double),parameter::u0=2.0495_double*cal, &
       u1=-4.0495_double*cal,u2=0.315_double*cal, &
       u3=6.414_double*cal
  real(kind=double),parameter::blue=0.0_double,green=0.0_double, &
       red=1.0_double,radius=0.3905_double

  character(len=1)::copy
  character(len=80)::fname
  integer::i,k,ndphi,ndr
  integer::n,nt,ntjob,ntprint,ntskip
  integer(kind=long),dimension(:),allocatable::agcc,agss,aphi
  real(kind=double)::aphi_id,c,c2,cosphi,cv,fact,p,phi,pi,r,r1,r2,snorm,temp
  real(kind=double)::dphi,dr,dt,t,vol,tau,zeta
  real(kind=double)::ak,ae,ae2,au,aud,aw
  real(kind=double),dimension(:,:),allocatable::x,y,z

! Pi

  pi=4.0_double*atan(1.0_double)

! Read checkpoint file

  write(unit=*,fmt="(a)",advance="no") &
       "     fname=[mdnbutane_out.dat] "
  read(unit=*,fmt="(a)") fname
  if(fname=="") then
    fname="mdnbutane_out.dat"
  end if

  open(unit=iin,file=fname,status="old",action="read", &
       form="unformatted",position="rewind")

  read(unit=iin) n,nt,ntjob,ntprint,ntskip

! Check for zero-length run

  if(nt<=0) then
    close(unit=iin)
    write(unit=*,fmt="(a)") "amdnbutane: empty file"
    stop
  end if

! Simulation parameters

  read(unit=iin) dphi,dr,dt,t,vol,tau,zeta

! Allocate arrays

  c=vol**(1.0_double/3.0_double)
  c2=0.5_double*c
  ndphi=int(180.0_double/dphi)
  ndr=int(0.5_double*c/dr)

  allocate(agcc(0:ndr-1),agss(0:ndr-1),aphi(0:ndphi-1), &
       x(m,n),y(m,n),z(m,n))

! Positions

  read(unit=iin) x,y,z
  read(unit=iin) !xl,yl,zl

! Accumulated averages, radial distribution functions &
! dihedral angle PDF

  read(unit=iin) ak,ae,ae2,au,aud,aw
  read(unit=iin) agcc
  read(unit=iin) agss
  read(unit=iin) aphi

  close(unit=iin)

! Suppress compiler warnings

  ntjob=ntjob
  ntprint=ntprint
  zeta=zeta

! Print results: simulation parameters

  write(unit=*,fmt="()")
  write(unit=*,fmt="(a,i14)") &
       "         n=",n
  write(unit=*,fmt="(a,f14.5,a)") &
       "         t=",t," (K)"
  write(unit=*,fmt="(a,f14.5,a)") &
       "       vol=",avog*vol/n," (cm**3/mol)"
  write(unit=*,fmt="(a,f14.5,a)") &
       "       tau=",tau, " (ps)"
  write(unit=*,fmt="(a,f14.5,a)") &
       "        dt=",dt, " (ps)"
  write(unit=*,fmt="(a,i14)") &
       "    ntskip=",ntskip

! Averages

  temp=ak/(0.5_double*(nt/ntskip)*nf*rg)
  p=(n/vol)/3.0_double*(ak-aw)/(nt/ntskip)/atm
  cv=n*(ae2/nt-(ae/nt)**2)/(rg*temp**2)

  write(unit=*,fmt="()")
  write(unit=*,fmt="(a,i14)") &
       "        nt=",nt
  write(unit=*,fmt="(a,es14.7,a)") &
       "       <t>=",temp," (K)"
  write(unit=*,fmt="(a,es14.7,a)") &
       " <u_inter>=",au/(nt/ntskip)," (kJ/mol)"
  write(unit=*,fmt="(a,es14.7,a)") &
       " <u_intra>=",aud/(nt/ntskip)," (kJ/mol)"
  write(unit=*,fmt="(a,es14.7,a)") &
       "       <p>=",p," (atm)"
  write(unit=*,fmt="(a,es14.7,a)") &
       "        Cv=",cv," (kJ/mol/K)"
  write(unit=*,fmt="()")

! Write g_cc(r) to file?

  write(unit=*,fmt="(a)",advance="no") &
       "   Write g_cc(r) to 'amdnbutane1.dat' [y] "
  read(unit=*,fmt="(a)") copy

  if(copy=="y".or.copy=="Y".or.copy=="") then

    open(unit=iout,file="amdnbutane1.dat",status="replace", &
         action="write",form="formatted")

    fact=3.0_double/(2.0_double*pi*n**2/vol*(nt/ntskip))
    do i=0,ndr-1
      r1=i*dr
      r2=(i+1)*dr
      r=0.5_double*(r1+r2)
      write(unit=iout,fmt="(tr1,f8.5,tr1,es12.5)") &
           r,fact*agcc(i)/(r2**3-r1**3)
    end do

    close(unit=iout)

  end if

! Write g_ss(r) to file?

  write(unit=*,fmt="(a)",advance="no") &
       "   Write g_ss(r) to 'amdnbutane2.dat' [y] "
  read(unit=*,fmt="(a)") copy

  if(copy=="y".or.copy=="Y".or.copy=="") then

    open(unit=iout,file="amdnbutane2.dat",status="replace", &
         action="write",form="formatted")

    fact=3.0_double/(2.0_double*pi*m**2*n**2/vol*(nt/ntskip))
    do i=0,ndr-1
      r1=i*dr
      r2=(i+1)*dr
      r=0.5_double*(r1+r2)
      write(unit=iout,fmt="(tr1,f8.5,tr1,es12.5)") &
           r,fact*agss(i)/(r2**3-r1**3)
    end do

    close(unit=iout)

  end if

! Write dihedral angle PDF to file?

  write(unit=*,fmt="(a)",advance="no") &
       "    Write p(phi) to 'amdnbutane3.dat' [y] "
  read(unit=*,fmt="(a)") copy

  if(copy=="y".or.copy=="Y".or.copy=="") then

    snorm=0.0_double
    do i=1,ndphi_id
      cosphi=cos((i-0.5_double)*pi/ndphi_id)
      snorm=snorm+exp(-(u0+cosphi*(u1+cosphi*(u2+cosphi*u3)))/(rg*t))
    end do
    snorm=snorm*180.0_double/ndphi_id

    open(unit=iout,file="amdnbutane3.dat",status="replace", &
         action="write",form="formatted")

    fact=1.0_double/((nt/ntskip)*n*dphi)
    do i=0,ndphi-1
      phi=(i+0.5_double)*dphi
      cosphi=cos(phi*pi/180.0_double)
      aphi_id=exp(-(u0+cosphi*(u1+cosphi*(u2+cosphi*u3)))/(rg*t))/snorm
      write(unit=iout,fmt="(tr1,f8.2,2(tr1,es12.5))") &
           phi,fact*aphi(i),aphi_id
    end do

    close(unit=iout)

  end if

! Write PDB file?

  write(unit=*,fmt="(a)",advance="no") &
       " Write PDF format to 'amdnbutane.pdb' [y] "
  read(unit=*,fmt="(a)") copy

  if(copy=="y".or.copy=="Y".or.copy=="") then

    open(unit=iout,file="amdnbutane.pdb",status="replace", &
         action="write",form="formatted")

    write(unit=iout,fmt="(a,3f9.3,3f7.2,tr2,a,tr11,a)") &
         "CRYST1",10.0_double*c,10.0_double*c,10.0_double*c, &
         90.0_double,90.0_double,90.0_double,"P1","1"

    do i=1,n
      do k=1,m
        write(unit=iout,fmt="(a,i5,tr19,3f8.3)") &
             "HETATM",m*(i-1)+k,10.0_double*(x(k,i)+c2), &
             10.0_double*(y(k,i)+c2),10.0_double*(z(k,i)+c2)
      end do
    end do

    write(unit=iout,fmt= &
         "(a,a,tr1,a,tr14,3f8.3,f6.2)") &
         "COLOR ","#####","####",red,green,blue,10.0_double*radius

    write(unit=iout,fmt="(a)") "END"

    close(unit=iout)

  end if

! Deallocate arrays

  deallocate(agcc,agss,aphi,x,y,z)

end program amdnbutane

!**********************************************************************!
