!**********************************************************************!
!
! File: agemclj.f90
!
! Gibbs Ensemble (NVT) Monte Carlo of Lennard-Jonesium
! Analyze checkpoint file
!
! 15-May-1999 (MN)
! 19-Apr-2012
!
!**********************************************************************!

program agemclj

!**********************************************************************!

  implicit none

  integer,parameter::long=selected_int_kind(18)
  integer,parameter::double=selected_real_kind(15)
  integer,parameter::iin=1,iout=2
  real,parameter::blue=0.0,green=0.0,radius=0.5,red=1.0

  character(len=1)::copy
  character(len=80)::fname
  integer::i,n,n1,n2,nt,ntjob,ntprint,ntskip
  real::c,disp1,disp2,dv,ptr,pvm,t,v,v1
  real,dimension(:),allocatable::x,y,z
  integer(kind=long)::accrp1,accrp2,accrt,accrv, &
       atryp1,atryp2,atryt1,atryt2,atryv,an1,an2
  real(kind=double):: amu1,ap1,arho1,au1,av1,amu2,ap2,arho2,au2,av2

! Read checkpoint file

  write(unit=*,fmt="(a)",advance="no") &
       "         fname=[gemclj_out.dat] "
  read(unit=*,fmt="(a)") fname
  if(fname=="") then
    fname="gemclj_out.dat"
  end if

  open(unit=iin,file=fname,status="old",action="read", &
       form="unformatted",position="rewind")

  read(unit=iin) n,n1,nt,ntjob,ntprint,ntskip!,iran

! Check for zero-length run

  if(nt<=0) then
    close(unit=iin)
    write(unit=*,fmt="(a)") " agemclj: empty file"
    stop
  end if

! Simulation parameters

  read(unit=iin) disp1,disp2,dv,ptr,pvm,t,v,v1

! Allocate arrays

  n2=n-n1

  allocate(x(n),y(n),z(n))

! Positions & accumulated averages

  read(unit=iin) x,y,z
  read(unit=iin) accrp1,accrp2,accrt,accrv, &
       atryp1,atryp2,atryt1,atryt2,atryv, &
       amu1,an1,ap1,arho1,au1,av1,amu2,an2,ap2,arho2,au2,av2

  close(unit=iin)

! Suppress compiler warnings

  ntjob=ntjob
  ntprint=ntprint

! Print results: simulation parameters

  write(unit=*,fmt="()")
  write(unit=*,fmt="(a,i10)") "            n1=",n1
  write(unit=*,fmt="(a,i10)") "            n2=",n2
  if(v<9999.9) then
    write(unit=*,fmt="(a,f10.5)") "             v=",v
  else
    write(unit=*,fmt="(a,es12.5)") "             v=",v
  end if
  write(unit=*,fmt="(a,f10.5)") "             t=",t
  write(unit=*,fmt="(a,f10.5)") "         disp1=",disp1
  write(unit=*,fmt="(a,f10.5)") "         disp2=",disp2
  write(unit=*,fmt="(a,f10.5)") "            dv=",dv
  write(unit=*,fmt="(a,f10.5)") "   prob(vmove)=",pvm
  write(unit=*,fmt="(a,f10.5)") "   prob(trans)=",ptr
  write(unit=*,fmt="()")

! Averages

  write(unit=*,fmt="(a,i12,a,i5,a)") &
       "            nt=",nt," (*",ntskip,")"
  write(unit=*,fmt="(a,es12.5)") &
       "        accrp1=",real(accrp1,double)/max(atryp1,1_long)
  write(unit=*,fmt="(a,es12.5)") &
       "        accrp2=",real(accrp2,double)/max(atryp2,1_long)
  write(unit=*,fmt="(a,es12.5)") &
       "         accrv=",real(accrv,double)/max(atryv,1_long)
  write(unit=*,fmt="(a,es12.5)") &
       "         accrt=",real(accrt,double)/max(atryt1+atryt2,1_long)
  write(unit=*,fmt="(a,es12.5)") &
       "          <N1>=",real(an1,double)/nt
  write(unit=*,fmt="(a,es12.5)") &
       "          <N2>=",real(an2,double)/nt
  write(unit=*,fmt="(a,es12.5)") &
       "          <V1>=",av1/nt
  write(unit=*,fmt="(a,es12.5)") &
       "          <V2>=",av2/nt
  write(unit=*,fmt="(a,es12.5)") &
       "        <rho1>=",arho1/nt
  write(unit=*,fmt="(a,es12.5)") &
       "        <rho2>=",arho2/nt
  write(unit=*,fmt="(a,es12.5)") &
       "       <U1/N1>=",au1/nt
  write(unit=*,fmt="(a,es12.5)") &
       "       <U2/N2>=",au2/nt
  write(unit=*,fmt="(a,es12.5)") &
       "          <p1>=",ap1/nt
  write(unit=*,fmt="(a,es12.5)") &
       "          <p2>=",ap2/nt
  if(amu1/=0.0_double.and.atryt1/=0_long) then
    write(unit=*,fmt="(a,es12.5)") &
         "         <mu1>=",-t*log(amu1/atryt1)
  end if
  if(amu2/=0.0_double.and.atryt2/=0_long) then
    write(unit=*,fmt="(a,es12.5)") &
         "         <mu2>=",-t*log(amu2/atryt2)
  end if
  write(unit=*,fmt="()")

! Write PDB files?

  write(unit=*,fmt="(a)",advance="no") &
       " Write PDB format to 'agemclj?.pdb'? [y] "
  read(unit=*,fmt="(a)") copy

  if(copy=="y".or.copy=="Y".or.copy=="") then

! Box 1

    c=v1**(1.0/3.0)

    open(unit=iout,file="agemclj1.pdb",status="replace",&
         action="write",form="formatted")

    write(unit=iout,fmt="(a,3f9.3,3f7.2,tr2,a,tr11,a)") &
         "CRYST1",c,c,c,90.0,90.0,90.0,"P1","1"
    do i=1,n1
      write(unit=iout,fmt="(a,i5,tr19,3f8.3)") &
           "HETATM",i,x(i)+0.5*c,y(i)+0.5*c,z(i)+0.5*c
    end do
    write(unit=iout,fmt= &
         "(a,a,tr1,a,tr14,3f8.3,f6.2)") &
         "COLOR ","#####","####",red,green,blue,radius
    write(unit=iout,fmt="(a)") "END"

    close(unit=iout)

! Box 2

    c=(v-v1)**(1.0/3.0)

    open(unit=iout,file="agemclj2.pdb",status="replace", &
         action="write",form="formatted")

    write(unit=iout,fmt="(a,3f9.3,3f7.2,tr2,a,tr11,a)") &
         "CRYST1",c,c,c,90.0,90.0,90.0,"P1","1"
    do i=n1+1,n
      write(unit=iout,fmt="(a,i5,tr19,3f8.3)") &
           "HETATM",i,x(i)+0.5*c,y(i)+0.5*c,z(i)+0.5*c
    end do
    write(unit=iout,fmt= &
         "(a,a,tr1,a,tr14,3f8.3,f6.2)") &
         "COLOR ","#####","####",red,green,blue,radius
    write(unit=iout,fmt="(a)") "END"

    close(unit=iout)

  end if

! Deallocate arrays

  deallocate(x,y,z)

end program agemclj

!**********************************************************************!
