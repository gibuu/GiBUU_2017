program testNucDLDA
  implicit none

  call testDFLDA

end program testNucDLDA


subroutine testDFLDA
  use NucDLDA
  use nucleusdefinition
  implicit none
  type(tnucleus)::nucleus
  integer::Mass,i
  character(15)::fileout


  write(*,*) 'Masse'
  read(*,*) Mass
  write(*,*) 'Dateiname'
  read(*,*) fileout
  nucleus%mass=Mass
  nucleus%charge=0

  call DFLDA(nucleus)


  open(1,file=fileout)

  do i=0,nucleus%maxIndex,1
     write(1,111) nucleus%dx*i, nucleus%densTab(i,1), nucleus%densTab(i,2)
     111  FORMAT(1X,F12.8,1x,F12.8,1x,F12.8)
  end do


  close(1)



end subroutine testDFLDA




subroutine testAapprox
  use NucDLDA
  implicit none

real, parameter::hbarc=197.326968,rho0=0.168, E0=-15.67/hbarc, M=938./hbarc, pi=3.141592654  
  real::a,b1,b2,b3,eta,ck,e0f,rhoat0,rhopat0,deltarho,Mass,ylast
  real::deltat,maxr,Massmin,ystartmin,Massin,energy,energy1
  real::e0fstart,deltae0f
  integer::steps,e0fdepth,e0fsteps,depth

  e0f=0.8
  rhoat0=0.14
  deltarho=0.0001
  rhopat0=0.
  steps=400
  deltat=0.01
  maxr=15
  depth=3
  e0fdepth=3
  e0fstart=0.65
  deltae0f=0.005
  e0fsteps=50




  b3=-1
  eta=10.8
  ck=0.3*(1/M)*(3./2.*pi**2)**(2./3.)
  a=eta/(8*M);
  call startcond(rho0,E0,ck,b3,b1,b2)
  write(*,*) 'gesuchte Masse'
  read(*,*) Massin
  
!!$  call shootinput(eta,b3,ck,a,b1,b2,E0,rhoat0,rhopat0,deltarho,steps,e0f,deltat,maxr,Massmin,ystartmin)

!!$call rhostepsearch(eta,b3,ck,a,b1,b2,E0,rhoat0,rhopat0,deltarho,steps,depth,e0f,deltat,maxr,Massmin,ystartmin)

  call Aapprox(eta,b3,ck,a,b1,b2,E0,rhoat0,rhopat0,deltarho,steps,depth,e0fdepth,deltat,maxr,Massin,Massmin,ystartmin,e0f,energy1)


!  call shoot(ystartmin,rhopat0,0.001,maxr,deltat,'diff1test01.dat',ck,a,e0*e0f,b1,b2,b3,Mass,ylast,energy)

  write(*,*) Massmin, ystartmin, e0f
  write(*,*) energy1/Massmin, betheweizs(Massmin)

!  call e0fstepsearch(eta,b3,ck,a,b1,b2,E0,rhoat0,rhopat0,deltarho,steps,depth,e0fstart,e0fsteps,deltae0f,deltat,maxr)

end subroutine testAapprox
