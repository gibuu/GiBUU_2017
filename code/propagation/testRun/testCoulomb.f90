program testModule
  use CoulombKorrektur
  implicit none

  integer, parameter :: charge1=1
  integer, parameter :: charge2=100
  real :: mass1, mass2
  real, dimension(1:3) :: rInit,pInit,p1,r1,p2,r2
  integer :: j,i

  !Needed to test analytic Coulomb trajectories
  real, parameter :: pi=3.14159
  integer :: numSteps=100
  real, dimension(3):: r,p
  real :: masse
  real :: phi
  integer :: charge,externalCharge

  pInit=(/0.05,0.,0./)  !in GeV
  mass1=0.14
  mass2=103.

  Write(*,*) 'Initialized velocities=',pInit/sqrt(mass1**2+Sum(pInit**2))&
                                     &,pInit/sqrt(mass2**2+Sum(pInit**2))
  Do j=1,10,1
     Print *,j
     rInit=(/-2000.,float(j),0./)

     !set first particle
     r1=rInit
     p1=pInit

     !set second particle
     r2=0.
     p2=0.

     call CoulpropaTwo(r1,p1,charge1,mass1,r2,p2,charge2,mass2,2.)
     !Assume second particle to rest at origin:
     r1=rInit
     p1=pInit
     call Coulpropa(r1,p1,charge1,mass1,charge2,2.,20000.)
  End do


  ! Teste analytische Vorhersagen
  charge=charge1
  externalCharge=charge2
  masse=mass1
  Do j=1,30,1
     Do i=0,numSteps-10
        p=pInit
        r=(/-50000.,float(j),0./)

        phi=ACos(r(1)/sqrt(r(1)**2+r(2)**2))-i*(pi/2)/numSteps
        call coulAnalytic(r(1:2),p(1:2),charge,masse,externalCharge,phi)
        Write (1,*) real(r(1)),real(r(2)), real(p(1)), real(p(2))

        p=pInit
        r=(/-5000.,float(j),0./)

        phi=ACos(r(1)/sqrt(r(1)**2+r(2)**2))-i*(pi)/numSteps
        call coulAnalytic2(r(1:2),p(1:2),charge,masse,externalCharge,phi)
        Write (2,*) real(r(1)),real(r(2)), real(p(1)), real(p(2))
     End do
 end do

contains

  !***************************************************************************
  ! NAME
  ! subroutine coulAnalytic(rInit,pInit,charge,masse,externalCharge,phi)
  ! FUNCTION
  ! analytische Loesung fuer das Zentralkraftproblem in zwei dimensionen
  ! siehe Greiner, "Mechanik I", Aufgabe 34.3
  ! 1993, 6. Auflage,Verlag Harri Deutsch, Frankfurt a. M.
  !***************************************************************************
  subroutine coulAnalytic(rInit,pInit,charge,masse,externalCharge,phi)
    use constants, only: alphaQED, hbarc
    !Input
    real,dimension(2), intent(inout) :: rInit,pInit !initial position and momentum
    real, intent(in) :: masse                        !mass of particle
    integer, intent(in) :: charge                   !charge of particle
    integer, intent(in) :: ExternalCharge           !fixed charge at origin
    real,intent(inout) :: phi                       !angle
    !Variables according to Greiner
    real :: alpha
    real :: L            !Angular Momentum
    real :: r0,phi0           !radial distance and angle at t=0
    real,dimension(2) :: velocityXY !Geschwindigekeit kartesisch
    real :: velocityR,velocityPhi   !Geschwindigkeit in Zylinderkoordinaten
    real :: Energy,mass
    real :: C,A,phase,gamma
    real :: rPhi
    logical,parameter :: debug=.false.


    If (debug) then
       Write(*,*) '**In Coulanalytic'
       write(*,*) '  Momentum=',pInit
       write(*,*) '  Position=',rInit
       write(*,*) '  Mass=',masse
       write(*,*) '  Charge=',charge
       write(*,*) '  External charge=',externalCharge
    end if

    mass=masse/hbarc !in fm^-1
    pInit=pInit/hbarc !in fm^-1
    Energy=Sqrt(pInit(1)**2+pInit(2)**2+mass**2) !Energy at start of simulation in fm^-1

    !Evaluate cylindric coordinates
    r0=Sqrt(rInit(1)**2+rInit(2)**2)
    phi0=ACos(rInit(1)/r0)

    If (debug) write(*,*) 'Input: r0=',r0,'phi0=',phi0

    velocityXY(1:2)=pInit(1:2)/Energy  !no dimension
    velocityR=rInit(1)/r0*velocityXY(1)+rInit(2)/r0*velocityXY(2)
    velocityPhi=1./r0*(cos(phi0)*velocityXY(2)-sin(phi0)*velocityXY(1)) !in fm**(-1)

    gamma=1./Sqrt(1-velocityXY(1)**2-velocityXY(2)**2)
    !Angular momentum:
    L=gamma*r0**2*velocityPhi ! in fm.  It is a conserved quantity.

    If (1-(alphaQED*charge*ExternalCharge/(mass*L))**2.lt.0) then
       return
       Print *,'alpha is complex. &
            & Sqrt(',1-(alphaQED*charge*ExternalCharge/(mass*L))**2,')'
       Print *, 'charge*Externalcharge=',charge*Externalcharge
       Print *,'Angular Momentum=',L
       stop
    end if
    alpha=SQRT(1-(alphaQED*charge*ExternalCharge/(mass*L))**2)

    C=(Energy*charge*ExternalCharge*alphaQED) &
         &      /(mass**2*L**2*alpha**2)  ! units in fm^(-1)

    A=sqrt((Energy/(L*mass*alpha**2))**2-(1./(alpha*L))**2) !units fm^(-1)
    If(debug) Print *, 'A',(Energy/(L*mass*alpha**2))**2-(1./(alpha*L))**2

    phase=-ACos((1./r0+C)/A)+alpha*phi0

    rPhi=1./(-C+A*Cos(alpha*phi-phase))

    If (debug) then
       Print *, 'phase=',phase
       Print *, (1./r0+C)/A
    End if
    rInit(1)=rPhi*cos(phi)
    rInit(2)=rPhi*sin(phi)

    gamma=energy/mass-charge*externalCharge*alphaQED/mass/rPhi
    velocityR=L*A*alpha*sin(alpha*phi-phase)/gamma
    velocityPhi=L/rPhi**2/gamma

    pInit(1)=(velocityR*cos(phi)-rPhi*velocityPhi*sin(phi))*(Energy-charge*externalCharge*alphaQED/rPhi)*hbarc
    pInit(2)=(velocityR*sin(phi)+rPhi*velocityPhi*cos(phi))*(Energy-charge*externalCharge*alphaQED/rPhi)*hbarc


    If (debug) then
       Print *, 'Phi=',phi,'radius=', rPhi
       Print *, 'position in (x,y):', rInit
    end if
  end subroutine coulAnalytic



  !***************************************************************************
  ! NAME
  ! subroutine coulAnalytic2(rInit,pInit,charge,masse,externalCharge,phi)
  ! FUNCTION
  ! analytische Loesung fuer das Zentralkraftproblem in zwei dimensionen
  ! siehe Greiner, "Mechanik I", Aufgabe 34.3
  ! 1993, 6. Auflage,Verlag Harri Deutsch, Frankfurt a. M.
  ! NOTE
  ! This version can also deal with complex alpha, which occurs at
  ! small impact parameters and/or small momentum.
  !***************************************************************************
  subroutine coulAnalytic2(rInit,pInit,charge,masse,externalCharge,phi)
    use constants, only: alphaQED, hbarc
    !Input
    real,dimension(2), intent(inout) :: rInit,pInit !initial position and momentum
    real, intent(in) :: masse                        !mass of particle
    integer, intent(in) :: charge                   !charge of particle
    integer, intent(in) :: ExternalCharge           !fixed charge at origin
    real,intent(inout) :: phi                       !angle
    !Variables according to Greiner
    real :: alpha,alphaQuadrat
    real :: L            !Angular Momentum
    real :: r0,phi0           !radial distance and angle at t=0
    real,dimension(2) :: velocityXY !Geschwindigekeit kartesisch
    real :: velocityR,velocityPhi   !Geschwindigkeit in Zylinderkoordinaten
    real :: Energy,mass
    real :: B,gamma
    real :: rPhi
    real :: u0,uDot0
    logical,parameter :: debug=.false.


    If (debug) then
       Write(*,*) '**In Coulanalytic'
       write(*,*) '  Momentum=',pInit
       write(*,*) '  Position=',rInit
       write(*,*) '  Mass=',masse
       write(*,*) '  Charge=',charge
       write(*,*) '  External charge=',externalCharge
    end if

    mass=masse/hbarc !in fm^-1
    pInit=pInit/hbarc !in fm^-1
    Energy=Sqrt(pInit(1)**2+pInit(2)**2+mass**2) !Energy at start of simulation in fm^-1


    !Evaluate starting values
    r0=Sqrt(rInit(1)**2+rInit(2)**2)
    phi0=ACos(rInit(1)/r0)
    velocityXY=pInit/Energy  !no dimension
    velocityR=cos(phi0)*velocityXY(1)+sin(phi0)*velocityXY(2)
    velocityPhi=1./r0*(cos(phi0)*velocityXY(2)-sin(phi0)*velocityXY(1)) !in fm**(-1)
    gamma=1./Sqrt(1-velocityXY(1)**2-velocityXY(2)**2)


    If (debug) write(*,*) 'Input: r0=',r0,'phi0=',phi0

    !Angular momentum:
    L=gamma*r0**2*velocityPhi ! in fm.  It is a conserved quantity.

    alphaQuadrat=1.-(alphaQED*charge*ExternalCharge/(mass*L))**2
    if (alphaQuadrat.gt.0)then
       !           return
       alpha=Sqrt(alphaquadrat)
    else
       alpha=Sqrt(-alphaquadrat)  !alpha=i*alpha=beta according to my notes
       !           Print *, "alpha imaginaer.b=",rInit(2)
       !           stop
    end if

    B=(Energy*charge*ExternalCharge*alphaQED) &
         &      /(mass**2*L**2.)  ! units in fm^(-1)

    u0=1/r0+B/alphaQuadrat
    uDot0=-1./L*gamma*velocityR


    if (alphaQuadrat.gt.0)then
       rphi=1./(u0*cos(alpha*(phi-phi0))+uDot0/alpha*sin(alpha*(phi-phi0))-B/alphaQuadrat)
       gamma=energy/mass-charge*externalCharge*alphaQED/mass/rPhi
       velocityR=-L*alpha*(-u0*sin(alpha*(phi-phi0))+uDot0/alpha*cos(alpha*(phi-phi0)))/gamma
    else
       rphi=1./(u0*cosh(alpha*(phi-phi0))+uDot0/alpha*sinh(alpha*(phi-phi0))-B/alphaQuadrat)
       gamma=energy/mass-charge*externalCharge*alphaQED/mass/rPhi
       velocityR=-L*alpha*(u0*sinh(alpha*(phi-phi0))+uDot0/alpha*cosh(alpha*(phi-phi0)))/gamma
    end if

    rInit(1)=rPhi*cos(phi)
    rInit(2)=rPhi*sin(phi)

    velocityPhi=L/rPhi**2/gamma

    pInit(1)=(velocityR*cos(phi)-rPhi*velocityPhi*sin(phi))&
         &       *(Energy-charge*externalCharge*alphaQED/rPhi)*hbarc
    pInit(2)=(velocityR*sin(phi)+rPhi*velocityPhi*cos(phi))&
         &       *(Energy-charge*externalCharge*alphaQED/rPhi)*hbarc

    If (debug) then
       Print *, 'Phi=',phi,'radius=', rPhi
       Print *, 'position in (x,y):', rInit
    end if
  end subroutine coulAnalytic2

end program testModule
