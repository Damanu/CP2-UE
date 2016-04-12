!**********************************************************************!
!
! File: amcme.f90
!
! NVT-Monte Carlo of OPLS-AA Methane
! Analyze checkpoint file
!
! 23-Apr-1999 (MN)
! 19-Apr-2012
!
!**********************************************************************!

program amcme

!**********************************************************************!

  implicit none

  integer,parameter::long=selected_int_kind(18)
  integer,parameter::double=selected_real_kind(15)
  integer,parameter::iin=1,iout=2,m=5
  real,parameter::avog=602.252,atm=0.101325*(avog/1000.0), &
       rg=8.3143e-03

  character(len=1)::copy
  character(len=80)::fname
  integer::i,k,n,ndr,nt,ntjob,ntprint,ntskip
  real::c,c2,cv,da,disp,dr,p,pi,r,r1,r2,rho,t
  real,dimension(3)::fact
  real,dimension(:,:),allocatable::x,y,z
  integer(kind=long),dimension(:,:),allocatable::ag
  real(kind=double)::accr,au,au2,aw

! Read checkpoint file

  write(unit=*,fmt="(a)",advance="no") " fname=[mcme_out.dat] "
  read(unit=*,fmt="(a)") fname
  if(fname=="") then
    fname="mcme_out.dat"
  end if

  open(unit=iin,file=fname,status="old",action="read", &
       form="unformatted",position="rewind")

  read(unit=iin) n,nt,ntjob,ntprint,ntskip!,iran

! Check for zero-length run

  if(nt<=0) then
    write(unit=*,fmt="(a)") " amcme: empty file"
    close(unit=iin)
    stop
  end if

! Simulation parameters

  read(unit=iin) da,disp,dr,rho,t

! Allocate arrays

  c=(n/rho)**(1.0/3.0)
  c2=0.5*c
  ndr=int(c2/dr)

  allocate(ag(1:3,0:ndr-1),x(m,n),y(m,n),z(m,n))

! Positions, accumulated averages & g_AB(r)

  read(unit=iin) x,y,z
  read(unit=iin) accr,au,au2,aw
  read(unit=iin) ndr
  read(unit=iin) ag

  close(unit=iin)

! Suppress compiler warnings

  ntjob=ntjob
  ntprint=ntprint

! Print results: simulation parameters

  write(unit=*,fmt="()")
  write(unit=*,fmt="(a,i12)") "     n=",n
  write(unit=*,fmt="(a,f12.5,a)") "   rho=",rho," (1/nm**3)"
  write(unit=*,fmt="(a,f12.5,a)") "     T=",t," (K)"
  write(unit=*,fmt="(a,f12.5,a)") "  disp=",disp," (nm)"
  write(unit=*,fmt="(a,f12.5,a)") "    da=",da," (deg)"
  write(unit=*,fmt="()")

! Averages

  cv=real(3.0*rg+n*(au2/nt-(au/nt)**2)/(rg*t**2))
  p=real(rho*(rg*t-aw/(3*nt)))

  write(unit=*,fmt="(a,i12,a,i5,a)") "    nt=",nt," (*",ntskip,")"
  write(unit=*,fmt="(a,es12.5)") "  accr=",accr/nt
  write(unit=*,fmt="(a,es12.5,a)") "   <U>=",au/nt," (kJ/mol)"
  write(unit=*,fmt="(a,es12.5,a)") "    Cv=",cv," (kJ/mol/K)"
  write(unit=*,fmt="(a,es12.5,a)") "     p=",p/atm," (atm)"
  write(unit=*,fmt="()")

! Write g_AB(r) to file?

  write(unit=*,fmt="(a)",advance="no") &
       " Write g_AB(r) to 'amcme.dat'? [y] "
  read(unit=*,fmt="(a)") copy

  if(copy=="y".or.copy=="Y".or.copy=="") then

    pi=4.0*atan(1.0)
    fact=(3.0/(2.0*pi*rho*n*nt))/(/1.0,8.0,16.0/)

    open(unit=iout,file="amcme.dat",status="replace", &
         action="write",form="formatted")

    do i=1,ndr-1
      r1=i*dr
      r2=(i+1)*dr
      r=0.5*(r1+r2)
      write(unit=iout,fmt="(tr1,f8.5,3(tr1,es12.5))") &
           r,fact*ag(:,i)/(r2**3-r1**3)
    end do

    close(unit=iout)

  end if

! Write PDB file?

  write(unit=*,fmt="(a)",advance="no") &
       " Write PDB format to 'amcme.pdb'? [y] "
  read(unit=*,fmt="(a)") copy

  if(copy=="y".or.copy=="Y".or.copy=="") then

    open(unit=iout,file="amcme.pdb",status="replace", &
         action="write",form="formatted")

    write(unit=iout,fmt="(a,3f9.3,3f7.2,tr2,a,tr11,a)") &
         "CRYST1",10.0*c,10.0*c,10.0*c,90.0,90.0,90.0,"P1","1"
    do i=1,n
      do k=1,m
        if(k==1) then
          write(unit=iout,fmt="(a,i5,tr2,a,tr16,3f8.3)") &
               "HETATM",m*(i-1)+k,"C",10.0*(x(k,i)+c2), &
               10.0*(y(k,i)+c2),10.0*(z(k,i)+c2)
        else
          write(unit=iout,fmt="(a,i5,tr2,a,tr16,3f8.3)") &
               "HETATM",m*(i-1)+k,"H",10.0*(x(k,i)+c2), &
               10.0*(y(k,i)+c2),10.0*(z(k,i)+c2)
        end if
      end do
    end do
    write(unit=iout,fmt="(a)") "END"

    close(unit=iout)

  end if

! Deallocate arrays

  deallocate(ag,x,y,z)

end program amcme

!**********************************************************************!
