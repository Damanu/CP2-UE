!**********************************************************************!
!
! File: amclj.f90
!
! NVT-Monte Carlo of Lennard-Jonesium
! Analyze checkpoint file
!
! 13-Apr-1997 (MN)
! 19-Apr-2012
!
!**********************************************************************!

program amclj

!**********************************************************************!

  implicit none

  character(len=1)::copy
  integer,parameter::long=selected_int_kind(18)
  integer,parameter::double=selected_real_kind(15)
  integer,parameter::iin=1,iout=2
  real,parameter::blue=0.0,green=0.0,radius=0.5,red=1.0

  character(len=80)::fname
  integer::i,n,ndr,nt,ntjob,ntprint,ntskip
  real::c,c2,cv,disp,dr,fact,p,pi,r,r1,r2,rho,t
  real,dimension(:),allocatable::x,y,z
  integer(kind=long),dimension(:),allocatable::ag
  real(kind=double)::accr,au,au2,aw

! Read checkpoint file

  write(unit=*,fmt="(a)",advance="no") " fname=[mclj_out.dat] "
  read(unit=*,fmt="(a)") fname
  if(fname=="") then
    fname="mclj_out.dat"
  end if

  open(unit=iin,file=fname,status="old",action="read", &
       form="unformatted",position="rewind")

  read(unit=iin) n,nt,ntjob,ntprint,ntskip!,iran

! Check for zero-length run

  if(nt<=0) then
    write(unit=*,fmt="(a)") " amclj: empty file"
    close(unit=iin)
    stop
  end if

! Simulation parameters

  read(unit=iin) disp,dr,rho,t

! Allocate arrays

  c=(n/rho)**(1.0/3.0)
  c2=0.5*c
  ndr=int(c2/dr)

  allocate(ag(0:ndr-1),x(n),y(n),z(n))

! Positions, accumulated averages & g(r)

  read(unit=iin) x,y,z
  read(unit=iin) accr,au,au2,aw
  read(unit=iin) ag

  close(unit=iin)

! Suppress compiler warnings

  ntjob=ntjob
  ntprint=ntprint

! Print results: simulation parameters

  write(unit=*,fmt="()")
  write(unit=*,fmt="(a,i9)") "     n=",n
  write(unit=*,fmt="(a,f9.5)") "   rho=",rho
  write(unit=*,fmt="(a,f9.5)") "     t=",t
  write(unit=*,fmt="(a,f9.5)") "  disp=",disp
  write(unit=*,fmt="()")

! Averages

  cv=real(1.5+n*(au2/nt-(au/nt)**2)/t**2)
  p=real(rho*(t-aw/(3*nt)))

  write(unit=*,fmt="(a,i12,a,i5,a)") "    nt=",nt," (*",ntskip,")"
  write(unit=*,fmt="(a,es12.5)") "  accr=",accr/nt
  write(unit=*,fmt="(a,es12.5)") " <U>/N=",au/nt
  write(unit=*,fmt="(a,es12.5)") "  Cv/N=",cv
  write(unit=*,fmt="(a,es12.5)") "     p=",p
  write(unit=*,fmt="()")

! Write g(r) to file?

  write(unit=*,fmt="(a)",advance="no") &
       " Write g(r) to 'amclj.dat'? [y] "
  read(unit=*,fmt="(a)") copy

  if(copy=="y".or.copy=="Y".or.copy=="") then

    pi=4.0*atan(1.0)
    fact=3.0/(2.0*pi*rho*n*nt)

    open(unit=iout,file="amclj.dat",status="replace",action="write", &
         form="formatted")

    do i=0,ndr-1
      r1=i*dr
      r2=(i+1)*dr
      r=0.5*(r1+r2)
      write(unit=iout,fmt="(tr1,f8.5,tr1,es12.5)") &
           r,fact*ag(i)/(r2**3-r1**3)
    end do

    close(unit=iout)

  end if

! Write PDB file?

  write(unit=*,fmt="(a)",advance="no") &
       " Write PDB format to 'amclj.pdb'? [y] "
  read(unit=*,fmt="(a)") copy

  if(copy=="y".or.copy=="Y".or.copy=="") then

    open(unit=iout,file="amclj.pdb",status="replace",action="write", &
         form="formatted")

    write(unit=iout,fmt="(a,3f9.3,3f7.2,tr2,a,tr11,a)") &
         "CRYST1",c,c,c,90.0,90.0,90.0,"P1","1"
    do i=1,n
      write(unit=iout,fmt="(a,i5,tr19,3f8.3)") &
           "HETATM",i,x(i)+c2,y(i)+c2,z(i)+c2
    end do
    write(unit=iout,fmt= &
         "(a,a,tr1,a,tr14,3f8.3,f6.2)") &
         "COLOR ","#####","####",red,green,blue,radius
    write(unit=iout,fmt="(a)") "END"

    close(unit=iout)

  end if

! Deallocate arrays

  deallocate(ag,x,y,z)

end program amclj

!**********************************************************************!
