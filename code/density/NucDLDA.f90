!******************************************************************************
!****m* /NucDLDA
! NAME
! module NucDLDA
! PURPOSE
! This module calculates the ground state density based on a local
! density approximation without charge or symmetry energy
! AUTHOR
! Birger Steinmüller
!******************************************************************************


module NucDLDA

  implicit none
  private

  public :: Aapprox,startcond,DFLDA,getEBinding,getEParticleLaplace,DFLDAWelke, &
            startcondWelke,getEParticleLaplaceWelke


contains


!!$  !**********************************************************************
!!$  !****f* NucDLDA/betheweizs
!!$  ! NAME
!!$  ! real function betheweizs(A)
!!$  ! PURPOSE
!!$  ! This function calculates the energy per particle using only volume
!!$  ! and surface energy terms
!!$  ! INPUTS
!!$  ! * real :: A -- Mass as a real number since the routines used in this
!!$  !                module calculate it this way
!!$  ! OUTPUT
!!$  ! Energy per particle in MeV
!!$  !**********************************************************************
!!$  real function betheweizs(A)
!!$    implicit none
!!$    real,intent(in)::a
!!$    real::av,ao
!!$    av=15.84;
!!$    ao=18.33;
!!$
!!$    betheweizs=(av*A-ao*A**(2./3.))/A;
!!$
!!$  end function betheweizs

  !****************************************************************************
  !****f* NucDLDA/diffeq
  ! NAME
  ! real function diffeq(y,yp,t,ck,a,e0f,b1,b2,b3)
  ! PURPOSE
  ! This is the differential equation used for calculating the density
  ! It comes from varying the energy-density-functional
  ! INPUTS
  ! * real :: rho -- current density
  ! * real :: rhop -- first derivative of the density wrt the radius
  ! * real :: r -- current radial position
  ! * real :: ck, a, e0f, b1, b2, b3 -- constants
  ! OUTPUT
  ! Second derivative of the density wrt the radius
  !****************************************************************************

  real function diffeq(rho,rhop,r,ck,a,e0f,b1,b2,b3)

    real, intent(in)::rho,rhop,r,ck,a,e0f,b1,b2,b3

    diffeq=-2./r*rhop+5./6.*ck/a*rho**(2./3.)-e0f/(2*a)+b1/a*rho+7./6.*b2/a*rho**(4./3.)+8./6.*b3/a*rho**(5./3.)

  end function diffeq

  !****************************************************************************
  !****f* NucDLDA/energyofrho
  ! NAME
  ! real function energyofrho(rho,rhop,t,ck,a,e0f,b1,b2,b3)
  ! PURPOSE
  ! This calculates the energy density for a given density and its
  ! first derivative
  ! INPUTS
  ! * real :: rho -- current density
  ! * real :: rhop -- first derivative of the density wrt the radius
  ! * real :: ck, a, b1, b2, b3 -- constants
  ! OUTPUT
  ! Energy density in units of fm**(-1)
  !****************************************************************************

  real function energyofrho(rho,rhop,b1,b2,b3,ck,a)
    real, intent(in)::rho,rhop,b1,b2,b3,ck,a
    energyofrho=rho*(ck*rho**(2./3.)+b1*rho+b2*rho**(4./3.)+b3*rho**(5./3.))+a*rhop**2
  end function energyofrho

  !****************************************************************************
  !****f* NucDLDA/getEBinding
  ! NAME
  ! real function getEBinding(E)
  ! PURPOSE
  ! Stores the calculated binding energy of the nucleus
  ! INPUTS
  ! * real :: E -- energy
  ! OUTPUT
  ! Energy in units of MeV
  !****************************************************************************

  real function getEBinding(E)
    real, intent(in)::E
    logical, save::firstrun=.true.
    real, save::Energy=0.

    if (firstrun) then
       Energy=E
       firstrun=.false.
    end if

    getEBinding=-Energy
  end function getEBinding

!!$  !**********************************************************************
!!$  !****f* NucDLDA/getEParticle
!!$  ! NAME
!!$  ! real function getEParticle(rho,gradrho2,p2)
!!$  ! PURPOSE
!!$  ! Calculates the Energy of a Particle with
!!$  ! INPUTS
!!$  ! * real :: rho -- density
!!$  ! * real :: gradrho2 --grad(density)**2
!!$  ! * real :: p2 -- 3-momentum**2 (in GeV**2)
!!$  ! OUTPUT
!!$  ! Energy in units of MeV
!!$  !**********************************************************************
!!$
!!$  real function getEParticle(rho,gradrho2,p2GeV)
!!$    use constants, only: hbarc, mN
!!$
!!$    real , intent(in):: rho, gradrho2,p2GeV
!!$    logical, save::firstrun=.true.
!!$    real, save::b1,b2,b3,ck,a,rho1,E0,eta
!!$    real, parameter:: Mass=mN/hbarc
!!$    real::p2
!!$
!!$    if (firstrun) then
!!$       call startcond(rho1,E0,ck,b3,b1,b2,a,eta)
!!$       firstrun=.false.
!!$    end if
!!$
!!$    if (rho.eq.0.) then
!!$       getEParticle=0.
!!$       return
!!$    end if
!!$    p2=p2GeV/(hbarc**2)
!!$
!!$    getEParticle=(rho*(p2/(2*Mass)+b1*rho+b2*rho**(4./3.)+b3*rho**(5./3.))+a*gradrho2)*hbarc*1000
!!$
!!$  end function getEParticle


  !****************************************************************************
  !****s* NucDLDA/getEParticleLaplace
  ! NAME
  ! real function getEParticleLaplace(rho,gradrho2,p2)
  ! PURPOSE
  ! Calculates the Energy of a Particle with
  ! INPUTS
  ! * real :: rho -- density
  ! * real :: laplacerho --div(grad(density))
  ! * real :: p2 -- 3-momentum**2 (in GeV**2)
  ! OUTPUT
  ! Energy in units of MeV, 1st component total energy, 2nd kinetic energy,
  ! 3rd volume energy and 4th surface energy
  !****************************************************************************

  subroutine getEParticleLaplace(Etemp,rho,laplacerho,p2GeV)
    use constants, only: hbarc, mN

    real , intent(in):: rho, laplacerho,p2GeV
    real, dimension(1:4), intent(out) :: Etemp
    logical, save::firstrun=.true.
    real, save::b1,b2,b3,ck,a,rho1,E0,eta
    real, parameter::Mass=mN/hbarc
    real::p2

    if (firstrun) then
       call startcond(rho1,E0,ck,b3,b1,b2,a,eta)
       firstrun=.false.
    end if

    p2=p2GeV/(hbarc**2)

    Etemp(1)=(p2/(2*Mass)+b1*rho+b2*rho**(4./3.)+b3*rho**(5./3.)-a*laplacerho)
    Etemp(2)=(p2/(2*Mass))
    Etemp(3)=(b1*rho+b2*rho**(4./3.)+b3*rho**(5./3.))
    Etemp(4)=(-a*laplacerho)

    Etemp = Etemp*hbarc*1000

  end subroutine getEParticleLaplace


  !****************************************************************************
  !****s* NucDLDA/Aapprox
  ! NAME
  ! subroutine Aapprox(b3,ck,a,b1,b2,E0,rhostart,rhopstart,deltarho,steps,depth,e0fdepth,
  ! deltat,maxr,Massin,Massout,ystart,e0f,energy)
  ! PURPOSE
  ! This searches for the best density at the center and Lagrange multiplier
  ! e0f to get a density distribution with its Mass (Massout) close to the
  ! given Mass (Massin)
  ! INPUTS
  ! * real :: rhostart -- starting density for the search
  ! * real :: rhopstart -- first derivative of the density
  ! * real :: deltarho -- starting stepsize wrt density
  ! * integer :: depth -- searching depth in rho direction
  ! * real :: deltat -- stepsize for solving the differential equation
  ! * real :: maxr -- maximum distance from the centre
  ! * real :: Massin -- Mass of the nucleus we are looking for
  ! * real :: b3, ck, a, b1, b2, E0 -- constants
  ! * integer :: steps -- steps in density 'direction'
  ! * integer :: e0fdepth -- searching depth for e0f
  ! OUTPUT
  ! * real :: Massout -- Mass of the nucleus (calculated)
  ! * real :: ystart -- calculated starting density of this nucleus
  ! * real :: e0f -- calculated Lagrange multiplier
  ! * real :: energy -- calculated total binding energy
  !****************************************************************************

  subroutine Aapprox(b3,ck,a,b1,b2,E0,rhostart,rhopstart,deltarho,steps,depth,&
       & e0fdepth,deltat,maxr,Massin,Massout,ystart,e0f,energy)
    real, intent(in)::b3,E0,rhostart,rhopstart,deltarho,deltat,maxr,ck,a,b1,b2,Massin
    real, intent(out)::Massout,ystart,e0f,energy
    integer, intent(in)::steps,e0fdepth,depth
    real::e0fstart,deltae0f,deltatemp,maxdelta,mindelta,totalmin,maxe0f,&
         & mine0f,e0ftemp,Masstemp,ytemp,energytemp!,ylast,minA,maxA
    integer::e0fsteps,i,j

    e0fstart=0.645
    deltae0f=0.27
    e0fsteps=50

    totalmin=0
    do j=1,e0fdepth,1
       maxdelta=0
       mindelta=0
       maxe0f=0
       mine0f=0
       !maxA=0
       !minA=0
       do i=0,e0fsteps,1
          e0ftemp=e0fstart+i*deltae0f/e0fsteps
          call rhostepsearch(b3,ck,a,b1,b2,E0,rhostart,rhopstart,deltarho,&
               & steps,depth,e0ftemp,deltat,maxr,Masstemp,ytemp,energytemp)
          deltatemp=1/(Masstemp-Massin)
          if (deltatemp>maxdelta) then
             maxdelta=deltatemp
             maxe0f=e0ftemp
             !maxA=Masstemp
          end if
          if (deltatemp<mindelta) then
             mindelta=deltatemp
             mine0f=e0ftemp
             !minA=Masstemp
          end if

          if (abs(deltatemp)>totalmin) then
             totalmin=abs(deltatemp)
             e0f=e0ftemp
             Massout=Masstemp
             ystart=ytemp
             energy=energytemp
          end if
       end do

       e0fstart=min(maxe0f,mine0f)
       deltae0f=abs(maxe0f-mine0f)
    end do

  end subroutine Aapprox

  !****************************************************************************
  !****s* NucDLDA/DFLDA
  ! NAME
  ! subroutine DFLDA(nucleus)
  ! PURPOSE
  ! Initialises the density of a nucleus for a given mass
  ! INPUTS
  ! type(tnucleus)::nucleus -- Mass and Charge are important
  ! OUTPUT
  ! type(tnucleus)::nucleus -- The initialised nucleus
  !****************************************************************************

  subroutine DFLDA(nucleus)
    use nucleusDefinition

    type(tnucleus), pointer::nucleus
    real::Charge,factor,a,b1,b2,b3,ck,E0,e0f,eta,rho0,rhoat0,rhopat0,deltarho
    real::deltat,maxr,Massmin,ystartmin,Massin,energy,dummy!,e0fstart,deltae0f
    integer::steps,e0fdepth,depth!,e0fsteps

    Massin=nucleus%mass
    Charge=nucleus%charge
    factor=Charge/Massin

    call startcond(rho0,E0,ck,b3,b1,b2,a,eta)

    rhoat0=0.14
    deltarho=0.0001
    rhopat0=0.
    steps=400
    deltat=nucleus%dx
    maxr=nucleus%MaxIndex*nucleus%dx
    depth=3
    e0fdepth=3
    !e0fstart=0.65
    !deltae0f=0.005
    !e0fsteps=50

    call Aapprox(b3,ck,a,b1,b2,E0,rhoat0,rhopat0,deltarho,steps,depth, &
         & e0fdepth,deltat,maxr,Massin,Massmin,ystartmin,e0f,energy)

!    write(*,*) e0f, ystartmin

    call shootfillarray(ystartmin,rhopat0,deltat,maxr,deltat,ck,a,e0*e0f,b1,b2,b3,factor,nucleus)

    nucleus%chemPot=e0f

    dummy=getEBinding(energy)

  end subroutine DFLDA


  !****************************************************************************
  !****s* NucDLDA/rhostepsearch
  ! NAME
  ! subroutine rhostepsearch(b3,ck,a,b1,b2,E0,rhostart,rhopstart,deltarho,steps,depth,e0f,deltat,maxr,Mass,ystart,energy)
  ! PURPOSE
  ! This searches for the best density at the center for a given Lagrange
  ! mulitplier e0f to get a smooth density distribution
  ! INPUTS
  ! * real :: rhostart -- starting density for the search
  ! * real :: rhopstart -- first derivative of the density
  ! * real :: deltarho -- starting stepsize wrt density
  ! * integer :: depth -- searching depth in rho direction
  ! * real :: deltat -- stepsize for solving the differential equation
  ! * real :: maxr -- maximum distance from the centre
  ! * real :: b3, ck, a, b1, b2, E0, e0f -- constants
  ! * integer :: steps -- steps in density 'direction'
  ! OUTPUT
  ! * real :: Mass -- Mass of the nucleus (calculated)
  ! * real :: ystart -- calculated starting density of this nucleus
  ! * real :: energy -- calculated total binding energy
  !****************************************************************************

  subroutine rhostepsearch(b3,ck,a,b1,b2,E0,rhostart,rhopstart,deltarho,steps,depth,e0f,deltat,maxr,Mass,ystart,energy)
    real, intent(in)::b3,E0,rhostart,rhopstart,deltarho,e0f,deltat,maxr,ck,a,b1,b2
    real, intent(out)::Mass,ystart,energy
    integer, intent(in)::steps,depth
    integer::i
    real::Massmin,ystartmin,rhostart1,deltarho1

    rhostart1=rhostart
    deltarho1=deltarho

    do i=1,depth,1
       call shootinput(b3,ck,a,b1,b2,E0,rhostart1,rhopstart,deltarho1,steps,e0f,deltat,maxr,Massmin,ystartmin,energy)
       rhostart1=ystartmin
       deltarho1=deltarho1/steps*1.2
    end do

    Mass=Massmin
    ystart=ystartmin

  end subroutine rhostepsearch

!!$  !**********************************************************************
!!$  !****s* NucDLDA/shoot
!!$  ! NAME
!!$  ! subroutine shoot(y0,yp0,tstart,tend,delta,fileout,ck,a,e0f,b1,b2,b3,Mass,ylast,energy)
!!$  ! PURPOSE
!!$  ! A first order Euler method for solving the equation diffeq with the
!!$  ! starting parameters and prints it into a file
!!$  ! INPUTS
!!$  ! * real :: y0 -- starting density
!!$  ! * real :: yp0 -- first derivative of the density
!!$  ! * real :: tstart -- point at which to start
!!$  ! * real :: tend -- maximum distance for which the diff. eq. is solved
!!$  ! * real :: delta -- stepsize for solving the differential equation
!!$  ! * real :: maxr -- maximum distance from the centre
!!$  ! * real :: b3, ck, a, b1, b2, E0, e0f -- constants for the diff. eq.
!!$  ! * character(15) :: fileout -- the filename in which to write the output
!!$  ! OUTPUT
!!$  ! * real :: Mass -- Mass of the nucleus (calculated)
!!$  ! * real :: ylast -- density at the point after the last one which is used
!!$  ! * real :: energy -- calculated total binding energy
!!$  !**********************************************************************
!!$
!!$
!!$  subroutine shoot(y0,yp0,tstart,tend,delta,fileout,ck,a,e0f,b1,b2,b3,Mass,ylast,energy)
!!$    use constants, only: pi, hbarc
!!$
!!$    real, intent(in)::y0,yp0,tstart,tend,delta,ck,a,e0f,b1,b2,b3
!!$    real, intent(out)::Mass,ylast,energy
!!$    character(15), intent(in)::fileout
!!$    real::tnow,ynow,ypnow,tnext,ynext,ypnext
!!$
!!$    open(1,file=fileout)
!!$    tnow=tstart
!!$    ynow=y0
!!$    ypnow=yp0
!!$    Mass=0.
!!$    energy=0
!!$
!!$    do while ((tnow<(tend-delta)).and.(ypnow<=0).and.(ynow>=0))
!!$       write(1,100) tnow,ynow,ypnow
!!$       100  FORMAT(1X,F12.8,1x,F12.8,1x,F12.8)
!!$       tnext=tnow+delta
!!$       ypnext=ypnow+diffeq(ynow,ypnow,tnow,ck,a,e0f,b1,b2,b3)*delta
!!$       ynext=ynow+ypnow*delta
!!$       if (ynext>=0) then
!!$          Mass=Mass+(ynow*tnow**2+ynext*tnext**2)/2.*(tnext-tnow)
!!$          energy=energy+(energyofrho(ynow,ypnow,b1,b2,b3,ck,a)*tnow**2&
!!$               & +energyofrho(ynext,ypnext,b1,b2,b3,ck,a)*tnext**2)/2.*(tnext-tnow)
!!$       end if
!!$       ylast=ynext
!!$       tnow=tnext
!!$       ynow=ynext
!!$       ypnow=ypnext
!!$    end do
!!$    Mass=Mass*4*pi
!!$    energy=-energy*4*pi*hbarc*1000
!!$    close(1)
!!$  end subroutine shoot

  !****************************************************************************
  !****s* NucDLDA/shootfillarray
  ! NAME
  ! subroutine subroutine shootfillarray(y0,yp0,tstart,tend,delta,ck,a,e0f,b1,b2,b3,factor,nucleus)
  ! PURPOSE
  ! A first order Euler method for solving the equation diffeq with the
  ! starting parameters and fills the density array of nucleus with it
  ! INPUTS
  ! * real :: y0 -- starting density
  ! * real :: yp0 -- first derivative of the density
  ! * real :: tstart -- point at which to start
  ! * real :: tend -- maximum distance for which the diff. eq. is solved
  ! * real :: delta -- stepsize for solving the differential equation
  ! * real :: b3, ck, a, b1, b2, E0, e0f -- constants for the diff. eq.
  ! * type(tnucleus) :: nucleus -- the nucleus definition
  ! * real :: factor -- charge/Mass for rellative density of n and p
  ! OUTPUT
  ! * real :: Mass -- Mass of the nucleus (calculated)
  ! * real :: ylast -- density at the point after the last one which is used
  ! * real :: energy -- calculated total binding energy
  !****************************************************************************

  subroutine shootfillarray(y0,yp0,tstart,tend,delta,ck,a,e0f,b1,b2,b3,factor,nucleus)
    use nucleusDefinition

    real, intent(in)::y0,yp0,tstart,tend,delta,ck,a,e0f,b1,b2,b3,factor
    type(tnucleus),pointer::nucleus
    real::tnow,ynow,ypnow,tnext,ynext,ypnext
    integer::i

    tnow=tstart
    ynow=y0
    ypnow=yp0
    nucleus%densTab(0,1)=factor*ynow
    nucleus%densTab(0,2)=(1.-factor)*ynow
    i=1

    do while ((tnow<(tend-delta)).and.(ypnow<=0).and.(ynow>=0).and.(i<=nucleus%MaxIndex))
       nucleus%densTab(i,1)=factor*ynow
       nucleus%densTab(i,2)=(1.-factor)*ynow
       tnext=tnow+delta
       ypnext=ypnow+diffeq(ynow,ypnow,tnow,ck,a,e0f,b1,b2,b3)*delta
       ynext=ynow+ypnow*delta
       tnow=tnext
       ynow=ynext
       ypnow=ypnext
       i=i+1
    end do
  end subroutine shootfillarray


  !****************************************************************************
  !****s* NucDLDA/shootinput
  ! NAME
  ! subroutine shootinput(b3,ck,a,b1,b2,E0,rhoat0,rhopat0,deltarho,step,e0f,deltat,maxr,Massmin,ystartmin,energy)
  ! PURPOSE
  ! This searches a given interval of the densityfor the best density at
  ! the center for a given Lagrange mulitplier e0f to get a smooth
  ! density distribution
  ! INPUTS
  ! * real :: rhoat0 -- starting density for the search
  ! * real :: rhopat0 -- first derivative of the density
  ! * real :: deltarho -- starting stepsize wrt density
  ! * real :: deltat -- stepsize for solving the differential equation
  ! * real :: maxr -- maximum distance from the centre
  ! * real :: b3, ck, a, b1, b2, E0, e0f -- constants
  ! * integer :: step -- steps in density 'direction'
  ! OUTPUT
  ! * real :: Massmin -- Mass of the nucleus (calculated)
  ! * real :: ystartmin -- calculated starting density of this nucleus
  ! * real :: energy -- calculated total binding energy
  !****************************************************************************

  subroutine shootinput(b3,ck,a,b1,b2,E0,rhoat0,rhopat0,deltarho,step,&
       & e0f,deltat,maxr,Massmin,ystartmin,energy)
    real, intent(in)::b3,E0,rhoat0,rhopat0,deltarho,e0f,deltat,maxr,ck,a,b1,b2
    real, intent(out)::Massmin,ystartmin,energy
    integer, intent(in)::step
    real::ymin,Massout,yminout,ystart,tstart
    integer::i

    Massmin=0
    ystartmin=0
    ymin=100
    tstart=deltat

    do i=0,step,1
       ystart=rhoat0+deltarho*i
       Massout=0
       call shootnowrite(ystart,rhopat0,tstart,maxr,deltat,ck,a,e0f*e0,b1,b2,b3,Massout,yminout,energy)
       if (yminout<0) exit
       if ((yminout>=0).and.(yminout<ymin)) then
          ymin=yminout
          ystartmin=ystart
          Massmin=Massout
       end if
    end do
  end subroutine shootinput


  !****************************************************************************
  !****s* NucDLDA/shootnowrite
  ! NAME
  ! subroutine shootnowrite(y0,yp0,tstart,tend,delta,fileout,ck,a,e0f,b1,b2,b3,Mass,ylast,energy)
  ! PURPOSE
  ! Just like shoot without writing into a file
  ! INPUTS
  ! * real :: y0 -- starting density
  ! * real :: yp0 -- first derivative of the density
  ! * real :: tstart -- point at which to start
  ! * real :: tend -- maximum distance for which the diff. eq. is solved
  ! * real :: delta -- stepsize for solving the differential equation
  ! * real :: maxr -- maximum distance from the centre
  ! * real :: b3, ck, a, b1, b2, E0, e0f -- constants for the diff. eq.
  ! OUTPUT
  ! * real :: Mass -- Mass of the nucleus (calculated)
  ! * real :: ylast -- density at the point after the last one which is used
  ! * real :: energy -- calculated total binding energy
  !****************************************************************************

  subroutine shootnowrite(y0,yp0,tstart,tend,delta,ck,a,e0f,b1,b2,b3,Mass,ylast,energy)
    use constants, only: pi, hbarc

    real, intent(in)::y0,yp0,tstart,tend,delta,ck,a,e0f,b1,b2,b3
    real, intent(out)::Mass,ylast,energy
    real::tnow,ynow,ypnow,tnext,ynext,ypnext

    tnow=tstart
    ynow=y0
    ypnow=yp0
    Mass=0.
    energy=0

    do while ((tnow<(tend-delta)).and.(ypnow<=0).and.(ynow>=0))
       tnext=tnow+delta
       ypnext=ypnow+diffeq(ynow,ypnow,tnow,ck,a,e0f,b1,b2,b3)*delta
       ynext=ynow+ypnow*delta
       if (ynext>=0) then
          Mass=Mass+(ynow*tnow**2+ynext*tnext**2)/2.*(tnext-tnow)
          energy=energy+(energyofrho(ynow,ypnow,b1,b2,b3,ck,a)*tnow**2&
               & +energyofrho(ynext,ypnext,b1,b2,b3,ck,a)*tnext**2)/2.*(tnext-tnow)
       end if
       ylast=ynext
       tnow=tnext
       ynow=ynext
       ypnow=ypnext
    end do
    Mass=Mass*4*pi
    energy=-energy*4*pi*hbarc*1000
  end subroutine shootnowrite

  !****************************************************************************
  !****s* NucDLDA/startcond
  ! NAME
  ! subroutine startcond(rho0,E0,ck,b3,b1,b2,a,eta)
  ! PURPOSE
  ! calculates the constants for the potential
  ! OUTPUT
  ! * real :: rho0 -- nuclear density where the potential has a minimum
  ! * real :: E0 -- depth of the potential
  ! * real :: ck, b3 -- other constants
  ! * real :: b1, b2 -- constant for potential
  !****************************************************************************

  subroutine startcond(rho0,E0,ck,b3,b1,b2,a,eta)
    use constants, only: pi, hbarc, mN

    real, intent(out)::rho0,E0,ck,b3,b1,b2,a,eta
    real::a1,a2,a3,a4,x1,x2
    real, parameter :: M = mN/hbarc

!    nuclear matter density of 0.168
!    b3=-1
!    eta=10.8
!    rho0=0.168

    b3=2.4
    eta=12.8
    rho0=0.15

    E0=-15.84/(hbarc*1000)
    ck=0.3*(1/M)*(3./2.*pi**2)**(2./3.)
    a=eta/(8*M)
    a1=rho0
    a2=rho0**(4./3.)
    a3=3*rho0**(1./3.)
    a4=4*rho0**(2./3.)
    x1=E0-ck*rho0**(2./3.)-b3*rho0**(5./3.)
    x2=-2*ck-5*b3*rho0

    b2=(x2-a3*x1/a1)/(a4-a2*a3/a1)
    b1=(x1-a2*b2)/a1

  end subroutine startcond

  !****************************************************************************
  !****s* NucDLDA/DFLDAWelke
  ! NAME
  ! subroutine DFLDAWelke(nucleus)
  ! PURPOSE
  ! Initialises the density of a nucleus for a given mass
  ! INPUTS
  ! type(tnucleus)::nucleus -- Mass and Charge are important
  ! OUTPUT
  ! type(tnucleus)::nucleus -- The initialised nucleus
  !****************************************************************************

  subroutine DFLDAWelke(nucleus)
    use nucleusDefinition

    type(tnucleus), pointer::nucleus
    real::Charge,factor,a,b1,b2,b3,E0,e0f,rho0,rhoat0,rhopat0,deltarho,lambda,alpha,C,ck,eta
    real::deltat,maxr,Massmin,ystartmin,Massin,energy,dummy!,e0fstart,deltae0f
    integer::steps,e0fdepth,depth!,e0fsteps

    Massin=nucleus%mass
    Charge=nucleus%charge
    factor=Charge/Massin

    call startcondWelke(rho0,E0,ck,b3,b1,b2,a,eta,lambda,alpha,C)

    rhoat0=0.14
    deltarho=0.0001
    rhopat0=0.
    steps=400
    deltat=nucleus%dx
    maxr=nucleus%MaxIndex*nucleus%dx
    depth=3
    e0fdepth=3
    !e0fstart=0.65
    !deltae0f=0.005
    !e0fsteps=50

    call AapproxWelke(b3,ck,a,b1,b2,E0,lambda,alpha,rhoat0,rhopat0,deltarho,steps,depth, &
         & e0fdepth,deltat,maxr,Massin,Massmin,ystartmin,e0f,energy)

!    write(*,*) e0f, ystartmin

    call shootfillarrayWelke(ystartmin,rhopat0,deltat,maxr,deltat,ck,a,e0*e0f,b1,b2,b3,lambda,alpha,factor,nucleus)

    nucleus%chemPot=e0f

    dummy=getEBinding(energy)

  end subroutine DFLDAWelke


  real function welkeexpect(rho,lambda)
    use constants, only: pi
    real, intent(in)::rho,lambda
    real :: pf

    pf=(3./2.*pi**2*rho)**(1./3.)
    welkeexpect=3.*pf/lambda - 1./2.*lambda/pf + 1./8.*(12.*lambda/pf+&
         &lambda**3/pf**3)*log((lambda**2+4*pf**2)/lambda**2)-4*atan(2*pf/lambda)

  end function welkeexpect


  real function welkefordiffeq(rho,lambda)
    use constants, only: pi
    real, intent(in)::rho,lambda
    real :: pf,dpfdrho

    pf=(3./2.*pi**2*rho)**(1./3.)
    dpfdrho=1./3.*(3./2.)**(1./3.)/rho**(2./3.)*pi**(2./3.)

    welkefordiffeq=dpfdrho*pf/(lambda*pi**2)*(log((lambda**2+4*pf**2)/lambda**2)*2*lambda**2-&
         &8*atan(2*pf/lambda)*pf*lambda+8*pf**2)
  end function welkefordiffeq


  real function welkediff(rho,lambda)
    use constants, only: pi
    real, intent(in)::rho,lambda
    real :: pf,dpfdrho

    pf=(3./2.*pi**2*rho)**(1./3.)
    dpfdrho=1./3.*(3./2.)**(1./3.)/rho**(2./3.)*pi**(2./3.)

    welkediff=-dpfdrho*3./8.*(log((lambda**2+4*pf**2)/lambda**2)*(lambda**4+4*lambda**2*pf**2)-&
         &4*pf**2*lambda**2-8*pf**4)/(pf**4*lambda)
  end function welkediff



  !****************************************************************************
  !****s* NucDLDA/startcondWelke
  ! NAME
  ! subroutine startcondWelke(rho0,E0,ck,b3,b1,b2,a,eta)
  ! PURPOSE
  ! calculates the constants for the potential
  ! OUTPUT
  ! * real :: rho0 -- nuclear density where the potential has a minimum
  ! * real :: E0 -- depth of the potential
  ! * real :: ck, b3 -- other constants
  ! * real :: b1, b2 -- constant for potential
  !****************************************************************************
  subroutine startcondWelke(rho0,E0,ck,b3,b1,b2,a,eta,lambda,alpha,C)
    use constants, only: pi, hbarc, mN

    real, intent(out)::rho0,E0,ck,b3,b1,b2,a,eta,lambda,alpha,C
    real::a1,a2,a3,a4,x1,x2
    real, parameter :: M = mN/hbarc

    b3=2.4
    eta=12.8
    rho0=0.15

    E0=-15.84/(hbarc*1000)
    ck=0.3*(1/M)*(3./2.*pi**2)**(2./3.)
    a=eta/(8*M)
    C=-63.6
    lambda=2.13
    alpha=C/(hbarc*1000)*1./(2*pi**2)*lambda**3/rho0
    a1=rho0
    a2=rho0**(4./3.)
    a3=3*rho0**(1./3.)
    a4=4*rho0**(2./3.)
    x1=E0-ck*rho0**(2./3.)-b3*rho0**(5./3.)-alpha*welkeexpect(rho0,lambda)
    x2=-2*ck-5*b3*rho0-3*alpha*rho0**(1./3.)*welkediff(rho0,lambda)

    b2=(x2-a3*x1/a1)/(a4-a2*a3/a1)
    b1=(x1-a2*b2)/a1

  end subroutine startcondWelke


  subroutine AapproxWelke(b3,ck,a,b1,b2,E0,lambda,alpha,rhostart,rhopstart,deltarho,steps,depth,&
       & e0fdepth,deltat,maxr,Massin,Massout,ystart,e0f,energy)
    real, intent(in)::b3,E0,rhostart,rhopstart,deltarho,deltat,maxr,ck,a,b1,b2,lambda,alpha,Massin
    real, intent(out)::Massout,ystart,e0f,energy
    integer, intent(in)::steps,e0fdepth,depth
    real::e0fstart,deltae0f,deltatemp,maxdelta,mindelta,totalmin,maxe0f,&
         & mine0f,e0ftemp,Masstemp,ytemp,energytemp!,ylast,minA,maxA
    integer::e0fsteps,i,j

    e0fstart=0.645
    deltae0f=0.27
    e0fsteps=50

    totalmin=0
    do j=1,e0fdepth,1
       maxdelta=0
       mindelta=0
       maxe0f=0
       mine0f=0
       !maxA=0
       !minA=0
       do i=0,e0fsteps,1
          e0ftemp=e0fstart+i*deltae0f/e0fsteps
          call rhostepsearchWelke(b3,ck,a,b1,b2,E0,lambda,alpha,rhostart,rhopstart,deltarho,&
               & steps,depth,e0ftemp,deltat,maxr,Masstemp,ytemp,energytemp)
          deltatemp=1/(Masstemp-Massin)
          if (deltatemp>maxdelta) then
             maxdelta=deltatemp
             maxe0f=e0ftemp
             !maxA=Masstemp
          end if
          if (deltatemp<mindelta) then
             mindelta=deltatemp
             mine0f=e0ftemp
             !minA=Masstemp
          end if

          if (abs(deltatemp)>totalmin) then
             totalmin=abs(deltatemp)
             e0f=e0ftemp
             Massout=Masstemp
             ystart=ytemp
             energy=energytemp
          end if
       end do

       e0fstart=min(maxe0f,mine0f)
       deltae0f=abs(maxe0f-mine0f)
    end do

  end subroutine AapproxWelke


  subroutine rhostepsearchWelke(b3,ck,a,b1,b2,E0,lambda,alpha,rhostart,rhopstart,deltarho,steps,&
       & depth,e0f,deltat,maxr,Mass,ystart,energy)
    real, intent(in)::b3,E0,rhostart,rhopstart,deltarho,e0f,deltat,maxr,ck,a,b1,b2,lambda,alpha
    real, intent(out)::Mass,ystart,energy
    integer, intent(in)::steps,depth
    integer::i
    real::Massmin,ystartmin,rhostart1,deltarho1

    rhostart1=rhostart
    deltarho1=deltarho

    do i=1,depth,1
       call shootinputWelke(b3,ck,a,b1,b2,E0,lambda,alpha,rhostart1,rhopstart,deltarho1,steps,e0f, &
            & deltat,maxr,Massmin,ystartmin,energy)
       rhostart1=ystartmin
       deltarho1=deltarho1/steps*1.2
    end do

    Mass=Massmin
    ystart=ystartmin

  end subroutine rhostepsearchWelke


  subroutine shootinputWelke(b3,ck,a,b1,b2,E0,lambda,alpha,rhoat0,rhopat0,deltarho,step,&
       & e0f,deltat,maxr,Massmin,ystartmin,energy)
    real, intent(in)::b3,E0,rhoat0,rhopat0,deltarho,e0f,deltat,maxr,ck,a,b1,b2,lambda,alpha
    real, intent(out)::Massmin,ystartmin,energy
    integer, intent(in)::step
    real::ymin,Massout,yminout,ystart,tstart
    integer::i

    Massmin=0
    ystartmin=0
    ymin=100
    tstart=deltat

    do i=0,step,1
       ystart=rhoat0+deltarho*i
       Massout=0
       call shootnowriteWelke(ystart,rhopat0,tstart,maxr,deltat,ck,a,e0f*e0,b1,b2,b3,lambda,alpha,Massout,yminout,energy)
       if (yminout<0) exit
       if ((yminout>=0).and.(yminout<ymin)) then
          ymin=yminout
          ystartmin=ystart
          Massmin=Massout
       end if
    end do
  end subroutine shootinputWelke


  subroutine shootnowriteWelke(y0,yp0,tstart,tend,delta,ck,a,e0f,b1,b2,b3,lambda,alpha,Mass,ylast,energy)
    use constants, only: pi, hbarc

    real, intent(in)::y0,yp0,tstart,tend,delta,ck,a,e0f,b1,b2,b3,lambda,alpha
    real, intent(out)::Mass,ylast,energy
    real::tnow,ynow,ypnow,tnext,ynext,ypnext

    tnow=tstart
    ynow=y0
    ypnow=yp0
    Mass=0.
    energy=0

    do while ((tnow<(tend-delta)).and.(ypnow<=0).and.(ynow>=0))
       tnext=tnow+delta
       ypnext=ypnow+diffeqWelke(ynow,ypnow,tnow,ck,a,e0f,b1,b2,b3,lambda,alpha)*delta
       ynext=ynow+ypnow*delta
       if (ynext>=0) then
          Mass=Mass+(ynow*tnow**2+ynext*tnext**2)/2.*(tnext-tnow)
          energy=energy+(energyofrhoWelke(ynow,ypnow,b1,b2,b3,ck,a,lambda,alpha)*tnow**2&
               & +energyofrhoWelke(ynext,ypnext,b1,b2,b3,ck,a,lambda,alpha)*tnext**2)/2.*(tnext-tnow)
       end if
       ylast=ynext
       tnow=tnext
       ynow=ynext
       ypnow=ypnext
    end do
    Mass=Mass*4*pi
    energy=-energy*4*pi*hbarc*1000
  end subroutine shootnowriteWelke


  real function diffeqWelke(rho,rhop,r,ck,a,e0f,b1,b2,b3,lambda,alpha)

    real, intent(in)::rho,rhop,r,ck,a,e0f,b1,b2,b3,lambda,alpha

    diffeqWelke=-2./r*rhop+5./6.*ck/a*rho**(2./3.)-e0f/(2*a)+b1/a*rho+7./6.*b2/a*rho**(4./3.)&
         & +8./6.*b3/a*rho**(5./3.)+alpha*welkefordiffeq(rho,lambda)/(2.*a)

  end function diffeqWelke


  real function energyofrhoWelke(rho,rhop,b1,b2,b3,ck,a,lambda,alpha)
    real, intent(in)::rho,rhop,b1,b2,b3,ck,a,lambda,alpha
    energyofrhoWelke=rho*(ck*rho**(2./3.)+b1*rho+b2*rho**(4./3.)+b3*rho**(5./3.))+a*rhop**2+alpha*rho*welkeexpect(rho,lambda)
  end function energyofrhoWelke


  subroutine shootfillarrayWelke(y0,yp0,tstart,tend,delta,ck,a,e0f,b1,b2,b3,lambda,alpha,factor,nucleus)
    use nucleusDefinition

    real, intent(in)::y0,yp0,tstart,tend,delta,ck,a,e0f,b1,b2,b3,factor,lambda,alpha
    type(tnucleus),pointer::nucleus
    real::tnow,ynow,ypnow,tnext,ynext,ypnext
    integer::i

    tnow=tstart
    ynow=y0
    ypnow=yp0
    nucleus%densTab(0,1)=factor*ynow
    nucleus%densTab(0,2)=(1.-factor)*ynow
    i=1

    do while ((tnow<(tend-delta)).and.(ypnow<=0).and.(ynow>=0).and.(i<=nucleus%MaxIndex))
       nucleus%densTab(i,1)=factor*ynow
       nucleus%densTab(i,2)=(1.-factor)*ynow
       tnext=tnow+delta
       ypnext=ypnow+diffeqWelke(ynow,ypnow,tnow,ck,a,e0f,b1,b2,b3,lambda,alpha)*delta
       ynext=ynow+ypnow*delta
       tnow=tnext
       ynow=ynext
       ypnow=ypnext
       i=i+1
    end do
  end subroutine shootfillarrayWelke


  subroutine getEParticleLaplaceWelke(Etemp,rho,laplacerho,p2GeV)
    use constants, only: hbarc, mN

    real , intent(in):: rho, laplacerho,p2GeV
    real, dimension(1:5), intent(out) :: Etemp
    logical, save::firstrun=.true.
    real, save::b1,b2,b3,ck,a,rho1,E0,eta,lambda,alpha,C
    real, parameter::Mass=mN/hbarc
    real::p2

    if (firstrun) then
       call startcondWelke(rho1,E0,ck,b3,b1,b2,a,eta,lambda,alpha,C)
       firstrun=.false.
    end if

    p2=p2GeV/(hbarc**2)


!    Etemp(1)=(p2/(2*Mass)+b1*rho+b2*rho**(4./3.)+b3*rho**(5./3.)-a*laplacerho+alpha*welkeexpect(rho,lambda))*hbarc*1000

    Etemp(2)=(p2/(2*Mass))*hbarc*1000
    Etemp(3)=(b1*rho+b2*rho**(4./3.)+b3*rho**(5./3.))*hbarc*1000
    Etemp(4)=(-a*laplacerho)*hbarc*1000
    Etemp(5)=alpha*welkeexpect(rho,lambda)*hbarc*1000
    if (Etemp(5).NE.Etemp(5)) then
       Etemp(5)=0.
    end if
    if (abs(Etemp(5)).ge.100.) then
       Etemp(5)=0.
    end if

    Etemp(1)=Etemp(2)+Etemp(3)+Etemp(4)+Etemp(5)

  end subroutine getEParticleLaplaceWelke


end module NucDLDA
