!**********************************************************************!
!
! File: amdlj.f90
!
! NVE Molecular Dynamics of Lennard-Jonesium
! Analyze checkpoint file
!
! 08-May-1997 (MN)
! 04-May-2012
!
!**********************************************************************!

program amdlj

!**********************************************************************!

  implicit none

  integer,parameter::long=selected_int_kind(18)
  integer,parameter::double=selected_real_kind(15)
  integer,parameter::iin=1,iout=2
  real,parameter::blue=0.0,green=0.0,radius=0.5,red=1.0

  character(len=1)::copy
  character(len=80)::fname
  integer::i,n,ncor,ndr,nt,ntaskip,ntcskip,ntjob,ntorig
  integer(kind=long),dimension(:),allocatable::ag
  real::c,c2,cv,dr,dt,fact,p,pi,r,r1,r2,rho,t
  real,dimension(:),allocatable::x,y,z
  real(kind=double)::ak,ak2,au,aw
  real(kind=double),dimension(:),allocatable::acf

! Read checkpoint file

  write(unit=*,fmt="(a)",advance="no") "          fname=[mdlj_out.dat] "
  read(unit=*,fmt="(a)") fname
  if(fname=="") then
    fname="mdlj_out.dat"
  end if

  open(unit=iin,file=fname,status="old",action="read", &
       form="unformatted",position="rewind")

  read(unit=iin) n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig!,ntprint

! Check for zero-length run

  if(nt<=0) then
    close(unit=iin)
    write(unit=*,fmt="(a)") " amdlj: empty file"
    stop
  end if

! Simulation parameters

  read(unit=iin) dr,dt,rho

! Allocate arrays

  ndr=int(0.5*(n/rho)**(1.0/3.0)/dr)

  allocate(acf(0:ncor),ag(0:ndr-1),x(n),y(n),z(n))

! Positions & (skip) velocities

  read(unit=iin) x,y,z
  read(unit=iin) !vx,vy,vz

! Accumulated averages

  read(unit=iin) ak,ak2,au,aw

! Pair distribution function g(r)

  read(unit=iin) ag

! Velocity autocorrelation function C(t)

  if(nt/ntcskip>ncor) then
    do i=0,ncor                       ! Skip ring buffer
      read(unit=iin) !vxt,vyt,vzt
    end do
    read(unit=iin) acf                ! C(t)
  end if

  close(unit=iin)

! Suppress compiler warnings

  ntjob=ntjob

! Print results: simulation parameters

  write(unit=*,fmt="()")
  write(unit=*,fmt="(a,i12)")   "              n=",n
  write(unit=*,fmt="(a,f12.5)") "            rho=",rho
  write(unit=*,fmt="(a,f12.5)") "             dt=",dt
  write(unit=*,fmt="(a,i12)")   "        ntaskip=",ntaskip
  write(unit=*,fmt="(a,i12)")   "        ntcskip=",ntcskip
  write(unit=*,fmt="(a,i12)")   "   ncor/ntcskip=",ncor
  write(unit=*,fmt="(a,i12)")   " ntorig/ntcskip=",ntorig

! Averages

  t=real(2.0*ak/(3.0*(nt/ntaskip)))
  cv=real(1.5/(1.0-1.5*n*(ak2/(nt/ntaskip)/(ak/(nt/ntaskip))**2-1.0)))
  p=real(rho*(t-aw/(3.0*(nt/ntaskip))))

  write(unit=*,fmt="()")
  write(unit=*,fmt="(a,i12)")    "             nt=",nt
  write(unit=*,fmt="(a,es12.5)") "            <t>=",t
  write(unit=*,fmt="(a,es12.5)") "          <k>/n=",ak/(nt/ntaskip)
  write(unit=*,fmt="(a,es12.5)") "          <u>/n=",au/(nt/ntaskip)
  write(unit=*,fmt="(a,es12.5)") "            e/n=",(ak+au)/(nt/ntaskip)
  write(unit=*,fmt="(a,es12.5)") "            <p>=",p
  write(unit=*,fmt="(a,es12.5)") "           cv/n=",cv
  write(unit=*,fmt="()")

! Write g(r) fo file?

  write(unit=*,fmt="(a)",advance="no") &
       "      Write g(r) to 'amdlj1.dat'? [y] "
  read(unit=*,fmt="(a)") copy

  if(copy=="".or.copy=="y".or.copy=="Y") then

    open(unit=iout,file="amdlj1.dat",status="replace",action="write", &
         form="formatted")

    pi=4.0*atan(1.0)
    fact=3.0/(2.0*pi*rho*n*(nt/ntaskip))
    do i=0,ndr-1
      r1=i*dr
      r2=(i+1)*dr
      r=0.5*(r1+r2)
      write(unit=iout,fmt="(tr1,f8.5,tr1,es12.5)") &
           r,fact*ag(i)/(r2**3-r1**3)
    end do

    close(unit=iout)

  end if

! Write C(t) to file?

  if(nt/ntcskip>ncor) then

    write(unit=*,fmt="(a)",advance="no") &
         "      Write C(t) to 'amdlj2.dat'? [y] "
    read(unit=*,fmt="(a)") copy

    if(copy=="".or.copy=="y".or.copy=="Y") then

      open(unit=iout,file="amdlj2.dat",status="replace", &
           action="write",form="formatted")

      do i=0,ncor
        write(unit=iout,fmt="(tr1,f8.5,2(tr1,es12.5))") i*ntcskip*dt, &
             acf(i)/((nt/ntcskip-(ncor+1))/ntorig+1), &   ! Unnormalized
             acf(i)/acf(0)                                ! Normalized
      end do

      close(unit=iout)

    end if

  end if

! Write PDB file?

  write(unit=*,fmt="(a)",advance="no") &
       " Write PDB format to 'amdlj.pdb'? [y] "
  read(unit=*,fmt="(a)") copy

  if(copy=="y".or.copy=="Y".or.copy=="") then

    open(unit=iout,file="amdlj.pdb",status="replace",action="write", &
         form="formatted")

    c=(n/rho)**(1.0/3.0)
    c2=0.5*c
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

  deallocate(acf,ag,x,y,z)

end program amdlj

!**********************************************************************!
