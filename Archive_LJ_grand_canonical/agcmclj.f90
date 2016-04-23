!**********************************************************************!
!
! File: agcmclj.f90
!
! Grand Canonical (TVmu) Monte Carlo of Lennard-Jonesium
! Analyze checkpoint file
!
! 08-Feb-2004 (MN)
! 19-Apr-2012
!
!**********************************************************************!

program agcmclj

!**********************************************************************!

  implicit none

  integer,parameter::long=selected_int_kind(18)
  integer,parameter::double=selected_real_kind(15)
  integer,parameter::iin=1,iout=2
  real,parameter::blue=0.0,green=0.0,radius=0.5,red=1.0

  character(len=1)::copy
  character(len=80)::fname
  integer::n,nt,ntjob,ntprint,ntskip
  integer::i,ndr
  real::disp,dr,pcre,t,v,zz
  real::c,c2,chit,fact,p,pi,r,r1,r2,rho,xmux
  real,dimension(:),allocatable::x,y,z
  integer(kind=long),dimension(:),allocatable::ag
  real(kind=double)::accrc,accrd,accrm,an,an2,au,aw

! Read checkpoint file

  write(unit=*,fmt="(a)",advance="no") "      fname=[gcmclj_out.dat] "
  read(unit=*,fmt="(a)") fname
  if(fname=="") then
    fname="gcmclj_out.dat"
  end if

  open(unit=iin,file=fname,status="old",action="read", &
       form="unformatted",position="rewind")

  read(unit=iin) n,nt,ntjob,ntprint,ntskip!,iran

! Check for zero-length run

  if(nt<=0) then
    close(unit=iin)
    write(unit=*,fmt="(a)") " agcmclj: empty file"
    stop
  end if

! Simulation parameters

  read(unit=iin) disp,dr,pcre,t,v,zz

! Allocate arrays

  ndr=int(0.5*v**(1.0/3.0)/dr)

  allocate(ag(0:ndr-1),x(n),y(n),z(n))

! Positions & accumulated averages

  read(unit=iin) x,y,z
  read(unit=iin) accrc,accrd,accrm,an,an2,au,aw
  read(unit=iin) ag

  close(unit=iin)

! Suppress compiler warnings

  ntjob=ntjob
  ntprint=ntprint

! Print results: simulation parameters

  write(unit=*,fmt="()")
  write(unit=*,fmt="(a,f14.7)") "          t=",t
  write(unit=*,fmt="(a,f14.7)") "          v=",v
  write(unit=*,fmt="(a,f14.7)") "          z=",zz
  write(unit=*,fmt="(a,f14.7)") "       disp=",disp
  write(unit=*,fmt="(a,f14.7)") " p(cre/del)=",pcre
  write(unit=*,fmt="()")

! Averages

  rho=real(an/(v*nt))
  xmux=t*log(zz/rho)
  p=real((an*t-aw/3.0)/(v*nt))
  chit=real((an2/an-an/nt)/(rho*t))

  write(unit=*,fmt="(a,i14,a,i5,a)") &
       "         nt=",nt," (*",ntskip,")"
  write(unit=*,fmt="(a,es14.7)") "      accrc=",accrc/nt
  write(unit=*,fmt="(a,es14.7)") "      accrd=",accrd/nt
  write(unit=*,fmt="(a,es14.7)") "      accrm=",accrm/nt
  write(unit=*,fmt="(a,es14.7)") "      mu_ex=",xmux
  write(unit=*,fmt="(a,es14.7)") "        <n>=",an/nt
  write(unit=*,fmt="(a,es14.7)") "      <rho>=",rho
  write(unit=*,fmt="(a,es14.7)") "    <u>/<n>=",au/an
  write(unit=*,fmt="(a,es14.7)") "          p=",p
  write(unit=*,fmt="(a,es14.7)") "      chi_t=",chit
  write(unit=*,fmt="()")

! Write g(r) to file?

  write(unit=*,fmt="(a)",advance="no") &
       "       Write g(r) to 'agcmclj.dat'? [y] "
  read(unit=*,fmt="(a)") copy

  if(copy=="y".or.copy=="Y".or.copy=="") then

    pi=4.0*atan(1.0)
    fact=real(3.0/(2.0*pi*rho*an))

    open(unit=iout,file="agcmclj.dat",status="replace", &
         action="write",form="formatted")

    do i=0,ndr-1
      r1=i*dr
      r2=(i+1)*dr
      r=0.5*(r1+r2)
      write(unit=iout,fmt="(tr1,f8.5,tr1,es14.7)") &
           r,fact*ag(i)/(r2**3-r1**3)
    end do

    close(unit=iout)

  end if

! Write PDB file?

  write(unit=*,fmt="(a)",advance="no") &
       " Write PDB format to 'agcmclj.pdb'? [y] "
  read(unit=*,fmt="(a)") copy

  if(copy=="y".or.copy=="Y".or.copy=="") then

    c=v**(1.0/3.0)
    c2=0.5*c

    open(unit=iout,file="agcmclj.pdb",status="replace", &
         action="write",form="formatted")

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

end program agcmclj

!**********************************************************************!
