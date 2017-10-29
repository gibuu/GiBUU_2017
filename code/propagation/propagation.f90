!******************************************************************************
!****m* /propagation
! NAME
! module propagation
! PURPOSE
! Module which includes the propagation of the test-particles.
!******************************************************************************
module propagation

  use CallStack, only: Traceback

  implicit none

  private

  !****************************************************************************
  !****g* propagation/Mode
  ! SOURCE
  !
  integer, save :: Mode = 2
  ! PURPOSE
  ! define the type of propagation:
  ! * 0: Cascade
  ! * 1: Euler
  ! * 2: PredictorCorrector
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/delta_E
  ! PURPOSE
  ! Delta energy in derivatives
  ! SOURCE
  !
  real, save :: delta_E=0.01
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/delta_P
  ! PURPOSE
  ! Delta Momentum in derivatives
  ! SOURCE
  !
  real, save :: delta_P=0.01
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/UseCoulombDirectly
  ! PURPOSE
  ! Whether to use coulomb force directly in propagation or not. (If switched
  ! off while coulomb is switched on in module coulomb, the effect of the
  ! coulomb potential comes in via the gradient of the potentials. With this
  ! flag you can not switch on/off coulomb, you just select, how it is
  ! treated.)
  ! SOURCE
  !
  logical,  save :: UseCoulombDirectly=.true.
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/UseHadronic
  ! PURPOSE
  ! Whether to use hadronic potentials in propagation
  ! SOURCE
  !
  logical,  save :: UseHadronic=.true.
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/dh_dp0_switch
  ! PURPOSE
  ! Switch which decides whether we use dh_dp0.
  ! SOURCE
  !
  logical,  save :: dh_dp0_switch=.true.
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/RungeKuttaOrder
  ! PURPOSE
  ! Order of Runge-Kutta in derivatives:
  ! * 1 = first order Runge-Kutta
  ! * 2 = second order Runge-Kuttay
  ! SOURCE
  !
  integer,  save :: RungeKuttaOrder=1
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/offShellInfo
  ! PURPOSE
  ! print out offShellInfo: set to .true. automatically if offShellTransport
  ! is used
  ! SOURCE
  !
  logical,save :: offShellInfo=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/offShellInfoDetail
  ! PURPOSE
  ! print out detailed offShellInfo
  ! SOURCE
  !
  logical,save :: offShellInfoDetail=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/tachyonDebug
  ! PURPOSE
  ! ...
  ! SOURCE
  !
  logical,save :: tachyonDebug=.false.
  !****************************************************************************

  logical,save :: startFlag = .true.
  logical,parameter :: warnNegativeMass2 = .false.

  !****************************************************************************
  !****s* propagation/updateVelocity
  ! NAME
  ! subroutine updateVelocity(teilchen)
  ! PURPOSE
  ! Updates velocities of a single particle or a field of particles
  ! and checks for v<c.
  ! INPUTS
  ! * type(particle),intent(inOut), dimension(:,:) :: teilchen
  ! or:
  ! * type(particle),intent(inOut), dimension(:)   :: teilchen
  ! or:
  ! * type(particle),intent(inOut) :: teilchen
  !****************************************************************************
  Interface updateVelocity
     Module Procedure updateVelocity_1, updateVelocity_field,updateVelocity_matrix
  End Interface updateVelocity


  public :: propagate
  public :: propagate_euler
  public :: propagate_cascade
  public :: updateVelocity
  public :: checkVelo
  public :: gradients

contains


  !****************************************************************************
  !****s* propagation/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads out namelist "propagation"
  ! INPUTS
  ! * (none)
  ! OUTPUT
  ! * Initializes global module variables
  !****************************************************************************
  subroutine init
    use output
    use offshellPotential

    integer :: ios

    !**************************************************************************
    !****n* propagation/Propagation
    ! NAME
    ! NAMELIST propagation
    ! PURPOSE
    ! Namelist which includes the input switches:
    ! * delta_P
    ! * delta_E
    ! * UseCoulombDirectly
    ! * UseHadronic
    ! * RungeKuttaOrder
    ! * Mode
    ! * dh_dp0_switch
    ! * offShellInfoDetail
    ! * tachyonDebug
    !**************************************************************************
    NAMELIST /propagation/ delta_P, delta_E, UseCoulombDirectly, UseHadronic, &
         RungeKuttaOrder, Mode, &
         dh_dp0_switch, offShellInfoDetail, tachyonDebug

    character(22), dimension(0:2), parameter :: cMode = (/ &
         'Cascade           ',&
         'Euler             ',&
         'PredictorCorrector' /)


    call Write_ReadingInput('propagation',0)
    rewind(5)
    read(5,nml=propagation,iostat=ios)
    if (ios>0) then
       write(*,*)
       write(*,*) 'Maybe you use an old jobcard.'
       write(*,*) 'Variabes have been renamed:'
       write(*,*) '   coulomb --> UseCoulombDirectly'
       write(*,*) '   hadronic --> Usehadronic'
       write(*,*) '   predictorCorrector --> Mode'
       write(*,*) '   DerivativeType --> RungeKuttaOrder'
    end if
    call Write_ReadingInput('propagation',0,ios)

    write(*,*) 'delta_P = ',delta_P
    write(*,*) 'delta_E = ',delta_E
    write(*,*) 'use coulomb force directly = ',UseCoulombDirectly
    write(*,*) 'hadronic potential flag    = ',UseHadronic
    write(*,*) 'RungeKuttaOrder            = ',RungeKuttaOrder
    select case (Mode)
    case (0:2)
       write(*,*) 'Propagation mode           : ',Mode,"= ",&
            cMode(Mode)
    case default
       write(*,*) 'Propagation mode           : ',Mode,"= WRONG!"
       call traceBack('wrong input value')
    end select
    write(*,*) 'dh_dp0_switch: ',dh_dp0_switch

    if (tachyonDebug) then
       write(*,*) 'Do tachyon debug!!'
    end if

    if (get_useOffShellPotentialBaryons().or.get_useOffShellPotentialMesons()) then
       offShellInfo=.true.
       if (tachyonDebug) offShellInfoDetail=.true.
       write(*,*) 'Detailed off-shell info?', offshellInfoDetail
    end if

    call Write_ReadingInput('propagation',1)

    startFlag = .false.
  end subroutine init




  !****************************************************************************
  !****s* propagation/gradients
  ! NAME
  ! subroutine gradients(part,Grad_P,Grad_R)
  ! PURPOSE
  ! determine gradients of the Hamiltonfunction
  ! INPUTS
  ! * type(particle),intent(in)  :: partIn  --
  !   particle whose gradients we want to calculate
  ! OUTPUT
  ! * real, dimension(1:3) :: Grad_P -- momentum gradient
  ! * real, dimension(0:3) :: Grad_R -- space gradients =(d/dt , d/dr(1:3) )
  !   [Note that there is no minus sign!!]
  !****************************************************************************
  subroutine gradients(partIn,Grad_P,Grad_R)

    use particleDefinition
    use coulomb, only: emfoca
    use densityModule
    use energyCalc
    use offShellPotential
    use minkowski
    use histf90
    use hist2Df90
    use idTable, only: photon
    use derivatives, only: finiteDifference

    type(particle),intent(in)  :: partIN
    real, intent(out), dimension(1:3) :: Grad_P
    real, intent(out), dimension(0:3),optional :: Grad_R

    real :: delta_R(1:3)
    real,dimension(-2:2,0:3) :: energy
    type(particle) :: part
    integer :: i,j
    real :: cpot                     ! Coulomb potential
    real, dimension(1:3) :: emForce  ! electromagnetic Force
    real :: dH_dp0           ! Derivative with respect to p_0: dH/dp_0
    type(histogram),save :: dh_dp0_hist
    type(histogram2D),save :: dh_dp0_hist2D
    integer, save :: numCalls=0
    logical, save :: DoInit=.true.
    logical :: outOfBounds!,success
    integer :: deriv_scheme

    numCalls=numCalls+1

    if (startFlag) call init

    if (get_offshell_debug()) then
       if (DoInit) then
          call createHist(dH_dp0_hist,'dH/dp_0',-2.,10.,0.1)
          call createHist2D(dH_dp0_hist2D,'Num events as function of dH/dp_0 and mass',(/0.6,-2./),(/1.3,10./),(/0.04,0.1/))
          DoInit=.false.
       end if
    end if

    ! Check Input
    if (partIn%ID <= 0) then !no particle
       write(*,*) 'ID =',partIn%ID
       call TRACEBACK('Error. Particle-ID<=0  in Gradient.')
    else if (partIn%ID.eq.photon) then
       ! Photons do not undergo any potentials
       if (present(grad_R)) grad_R=0.
       grad_P(1:3)  =  partIn%momentum(1:3)/partIn%momentum(0)
       if (Dot_Product(grad_P,grad_P).gt.1.0000000001) then
          write(*,*) 'Problem with photon in gradients'
          write(*,*) '|Velo|^2=',Dot_Product(grad_P,grad_P)
          write(*,*) 'Velo=',grad_P
          write(*,*) 'Momentum=',partIn%momentum
          return
       end if
    end if

    grad_P=0.

    if (present(grad_R)) then
       grad_R=0.
       if (UseCoulombDirectly) then ! special treatment of electromagnetic part
          cpot = emfoca(partIn%position(1:3),partIn%momentum(1:3),partIn%charge,partIn%ID,emForce)
          grad_R(1:3)=grad_R(1:3)-emForce  !Grad(Potential)=-Force
       end if
    end if

    if (UseHadronic) then !with hadronic potentials

       if (present(grad_R)) then
          ! Determine Gridsize in position space
          delta_R(1:3)=getGridSpacing()

          ! Derivative with respect to r: grad_r H
          part=partIN
          energy=0.
          do i=1,3 !loop over x,y,z
             deriv_scheme=0
             part%position=partIn%position
             do j=-RungeKuttaOrder,RungeKuttaOrder
                part%position(i)=partIN%position(i)+float(j)*delta_R(i)
                ! determine energy due to hadronic interaction
                if (treatParticleOffShell(part%Id,part%offshellparameter)) then
                   energy(j,i)=hamiltonFunc_offshell(part,outOFBounds)
                   if (outOfBounds) then
                      call OutOfBoundsMessage('r')
                      deriv_scheme=-j
                   end if
                else
                   call energyDetermination(part,&
                        ForbidCoulomb=UseCoulombDirectly,skipTestE0=.true.)
                   energy(j,i)=part%momentum(0)
                end if
             end do
             grad_R(i)=grad_R(i)+finiteDifference(energy(:,i),delta_R(i),RungeKuttaOrder,deriv_scheme)
          end do
       end if



       ! Derivative with respect to p: grad_p H
       part=partIN
       energy=0.

       do i=1,3 !loop over x,y,z
          deriv_scheme=0
          part%momentum=partIn%momentum
          do j=-RungeKuttaOrder,RungeKuttaOrder
             part%momentum(i)=partIN%momentum(i)+float(j)*delta_P
             ! determine energy due to hadronic interaction
             if (treatParticleOffShell(part%Id,part%offshellparameter)) then
                energy(j,i)=hamiltonFunc_offshell(part,outOfBounds)
                if (outOfBounds) then
                   call OutOfBoundsMessage('p')
                   deriv_scheme=-j
                end if
             else
                call energyDetermination(part,&
                     ForbidCoulomb=UseCoulombDirectly,skipTestE0=.true.)
                energy(j,i)=part%momentum(0)
             end if
          end do
          grad_P(i)=grad_P(i)+finiteDifference(energy(:,i),delta_P,RungeKuttaOrder,deriv_scheme)
       end do



       !grad_R(0)=0.

       if (treatParticleOffShell(partIn%Id,partIn%offshellparameter).and.dh_dp0_switch) then
          !********************************************************************
          ! Derivative with respect to p_0: dH/dp_0
          part=partIN
          energy=0.
          part%momentum=partIn%momentum
          deriv_scheme=0
          do j=-RungeKuttaOrder,RungeKuttaOrder
             part%momentum(0)=partIN%momentum(0)+float(j)*delta_E
             ! determine energy due to hadronic interaction
             energy(j,0)=hamiltonFunc_offshell(part,outofBounds)
             if (outOfBounds) then
                call OutOfBoundsMessage('E')
                deriv_scheme=-j
             end if
          end do
          dH_dp0=finiteDifference(energy(:,0),delta_E,RungeKuttaOrder,deriv_scheme)

          if (offShellInfoDetail) then
             write(*,*) 'dH_dp0=', dh_dp0
             if (.not.(dH_dp0>=0. .or. dH_dp0<=0.)) then
                write(*,*) "NaN in dH_dp0   ID=",partIn%Id, "     offShellParam=",partIn%offshellparameter
                write(*,*) "energy=", energy(:,0)
                stop
             end if
          end if


          ! See Oliver's notes for the case of an energy dependent Hamiltonian:
          if (present(grad_R)) grad_R(1:3)=grad_R(1:3)/(1-dH_dp0)
          grad_P=grad_P/(1-dH_dp0)


          if (get_offshell_debug()) then

             call addHist(dH_dp0_hist,dH_dp0,1.)
             call AddHist2D(dH_dp0_hist2D, (/part%mass,dH_dp0/),1.)

             if (mod(numCalls,1000).eq.0) then
                open(33,file='dH_dp0.dat')
                open(44,file='dH_dp0_2D.dat')
                call writeHist(dh_dp0_hist,33,mul=dH_dp0_hist%xbin)
                call writeHist2D_Gnuplot(dh_dp0_hist2D,44,mul=dH_dp0_hist2D%xbin(1)*dH_dp0_hist2D%xbin(2))
                close(33)
                close(44)
             end if
          end if
       end if
    else !no hadronic potentials
       !grad_R=0.
       grad_P=grad_P+partIN%momentum(1:3)/FreeEnergy(partIN)
    end if

  contains

    subroutine OutOfBoundsMessage(C)

      character*(*), intent(in) :: C
      if (offShellInfoDetail) &
           write(*,'(A,A,i4)') C,': Hamilton function off shell is out of bounds!!',j
      if (j.eq.0) then ! Out of bounds at central value
         write(*,'(A,A,I4,A,I4)') C,': Problem with outofBounds: j=',j
      else if (deriv_scheme*j.gt.0) then ! Out of bounds on both sides
         write(*,'(A,A,I4,A,I4)') C,': Problem with outofBounds: deriv_scheme='&
              & ,deriv_scheme," new deriv_scheme=",-j
      end if

    end subroutine OutOfBoundsMessage


  end subroutine gradients

  !****************************************************************************
  !****s* propagation/propagate
  ! NAME
  ! subroutine propagate(realPart, pertPart, delta_T, timeStep)
  ! PURPOSE
  ! This routine propagates the particle vectors.
  ! INPUTS
  ! * type(particle),intent(inOUT),dimension(:,:)  :: realPart
  !   -- real particle vector which should be propagated
  ! * type(particle),intent(inOUT),dimension(:,:)  :: pertPart
  !   -- perturbative particle vector which should be propagated
  ! * real, intent(in)   :: delta_T    -- time step size (fm/c)
  ! * integer,intent(in) :: timeStep   -- time step number
  ! OUTPUT
  ! * type(particle),intent(inOUT),dimension(:,:)  :: realPart
  ! * type(particle),intent(inOUT),dimension(:,:)  :: pertPart
  ! NOTES
  ! * Uses cascade, Euler or predictor-corrector time stepping.
  !****************************************************************************
  subroutine propagate(realPart, pertPart, delta_T, timeStep)
    use particleDefinition
    use inputGeneral, only: freezeRealParticles
    use densityModule, only: updateDensity
    use yukawa, only: updateYukawa
    use coulomb, only: updateCoulomb
    use offshellPotential, only: get_useOffShellPotentialMesons, get_useOffShellPotentialBaryons

    type(particle),intent(inOut),dimension(:,:) :: realPart
    type(particle),intent(inOut),dimension(:,:) :: pertPart
    real, intent(in) :: delta_T
    integer, intent(in) :: timeStep

    logical :: doReal
    ! Gradients of predictor step
    real, dimension(:,:,:), Allocatable :: gradR_real, gradP_real
    real, dimension(:,:,:), Allocatable :: gradR_pert, gradP_pert

    if (startFlag) call init

    doReal = (.not.freezeRealParticles) .or. (timeStep<=1)

    select case (Mode)
    case (0) ! === Cascade ===

       if (doReal) call propagate_cascade(realPart,delta_T)
       call propagate_cascade(pertPart,delta_T)

    case (1) ! === Euler ===

       if (doReal) call propagate_euler(realPart,delta_T)
       call propagate_euler(pertPart,delta_T)

    case (2) ! === PredictorCorrector ===

       ! Allocating and deallocating the arrays for every timestep is not
       ! optimal! This should be improved!!!

       ! Allocate fields to store the gradients
       if (doReal) then
          allocate(gradR_real(0:3, &
               lbound(realPart,dim=1):ubound(realPart,dim=1), &
               lbound(realPart,dim=2):ubound(realPart,dim=2)))
          allocate(gradP_real(1:3, &
               lbound(realPart,dim=1):ubound(realPart,dim=1), &
               lbound(realPart,dim=2):ubound(realPart,dim=2)))
       end if
       allocate(gradR_pert(0:3, &
            lbound(pertPart,dim=1):ubound(pertPart,dim=1), &
            lbound(pertPart,dim=2):ubound(pertPart,dim=2)))
       allocate(gradP_pert(1:3, &
            lbound(pertPart,dim=1):ubound(pertPart,dim=1), &
            lbound(pertPart,dim=2):ubound(pertPart,dim=2)))

       ! (1) Predictor Step: Define the predicted value of the particle at Delta_T

       if (doReal) call predictorStep(realPart,GradP_real,GradR_real,delta_T)
       call predictorStep(pertPart,GradP_pert,GradR_pert,delta_T)

       ! (2) Update potentials
       if (doReal) then
          call updateDensity(realPart)
          call updateCoulomb
          call updateYukawa(.false.)
       end if

       ! (3) Corrector step: consider also the predicted gradients
       if (doReal) then
          call correctorStep(realPart,GradP_real,GradR_real,delta_T)
          deallocate(gradR_real,gradP_real)
       end if

       call correctorStep(pertPart,GradP_pert,GradR_pert,delta_T)
       deallocate(gradR_pert,gradP_pert)

    end select

    !  ** Set masses properly:
    if (get_useOffShellPotentialBaryons() .or. get_useOffShellPotentialMesons()) then
       if (doReal) call setMass(realPart)
       call setMass(pertPart)
    end if

    if (doReal) call checkVelos(realPart)
    call checkVelos(pertPart)

  end subroutine propagate


  !****************************************************************************
  !****s* propagation/setMass
  ! NAME
  ! subroutine setMass(Parts)
  ! PURPOSE
  ! * Resets %mass according to the offshell parameter and the actual momentum and position of the particle.
  ! * if particles is NOT meant to be treated offshell, nothing is changed.
  ! INPUTS
  ! * type(particle),dimension(:,:),intent(inout) :: Parts
  !
  ! OUTPUT
  ! * type(particle),dimension(:,:),intent(inout) :: Parts
  !****************************************************************************
  subroutine setMass(Parts)
    use particleDefinition
    use offshellPotential, only:treatParticleOffShell!, getOffShellMass
    use minkowski, only: SP
    use potentialModule, only: massDetermination, scapot

    type(particle),dimension(:,:),intent(inout),target :: Parts

    integer :: iEns,iPart
    logical :: success
    type(particle), pointer :: pPart

    do iEns=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do iPart=lbound(Parts,dim=2),ubound(Parts,dim=2)

          pPart => Parts(iEns,iPart)
          if (pPart%ID.le.0) cycle

          if (.not.treatParticleOffShell(pPart%Id,pPart%offshellparameter)) cycle

          call massDetermination(pPart,success=success)

          if (.not.success) then
             write(*,*) "Can't set baremass in setMass!!!"
             write(*,*) "Particle is deleted!!!"
             call protocolFailure(pPart,'setMass1')
             pPart%ID=0
             cycle
          end if
          if (SP(pPart%momentum,pPart%momentum)<=0.01**2) then ! 10 MeV
             write(*,*)
             write(*,*) 'WARNING in setMass: p^mu p_mu too small! ',SP(pPart%momentum,pPart%momentum)
             write(*,*) 'particle is deleted!!'
             call protocolFailure(pPart,'setMass2')
             pPart%ID = 0 ! delete particle
             cycle
          end if
          if (pPart%mass<0.01) then ! 10MeV
             write(*,*)
             write(*,*) 'WARNING in setMass: baremass too small! ',pPart%mass
             write(*,*) 'particle is deleted!!'
             call protocolFailure(pPart,'setMass3')
             pPart%ID = 0 ! delete particle
             cycle
          else if (pPart%mass+scapot(pPart)>pPart%momentum(0)) then ! inv. mass > energy !
             write(*,*)
             write(*,*) 'WARNING in setMass: mass too large! ',pPart%mass
             write(*,*) 'particle is deleted!!'
             call protocolFailure(pPart,'setMass4')
             pPart%ID = 0 ! delete particle
             cycle
          end if

       end do
    end do
  end subroutine setMass


  !****************************************************************************
  !****s* propagation/predictorStep
  ! NAME
  ! subroutine predictorStep(Parts,gradP,gradR,delta_T)
  ! PURPOSE
  ! * Routine propagates a particle vector using an Euler method and returns the gradients.
  ! INPUTS
  ! * real, intent(in) :: delta_T
  !
  ! OUTPUT
  ! type(particle),intent(inOut),dimension(:,:) :: Parts -- vector at predicted position and momentum
  ! real, dimension(:,:,:),intent(out) :: gradR, gradP ! Gradients in R and P-Space
  !
  ! NOTES
  ! * Used as  predictor step in a predictor/corrector scheme
  !****************************************************************************
  subroutine predictorStep(Parts,gradP,gradR,delta_T)
    use offshellPotential, only: &
         get_useOffShellPotentialBaryons,get_useOffShellPotentialMesons
    use particleDefinition
    use minkowski, only:SP
    use output
    use IdTable, only: isBaryon

    real, dimension(1:,:,:),intent(out) :: gradP
    real, dimension(0:,:,:),intent(out) :: gradR
    type(particle),intent(inOut),dimension(:,:),target :: Parts

    real, intent(in) :: delta_T
    integer :: iEns,iPart
    type(particle), pointer :: pPart

    do iEns=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do iPart=lbound(Parts,dim=2),ubound(Parts,dim=2)

          pPart => Parts(iEns,iPart)
          if (pPart%ID < 0) exit
          if (pPart%ID == 0) cycle

          ! Evaluate dp/dt=-dH/dr and dr/dt=dH/dp
          call gradients(pPart,GradP(:,iEns,iPart),GradR(:,iEns,iPart))

          ! Set the perturbative particle vector to the predictor values:
          pPart%momentum(1:3)=pPart%momentum(1:3)-delta_T*gradR(1:3,iEns,iPart)
          pPart%momentum( 0 )=pPart%momentum( 0 )+delta_T*gradR(0,iEns,iPart)
          pPart%position(1:3)=pPart%position(1:3)+delta_T*gradP(1:3,iEns,iPart)
          pPart%velocity(1:3)=gradP(1:3,iEns,iPart)

          if (SP(pPart%momentum,pPart%momentum) < -epsilon(0.0D0)) then
             if (warnNegativeMass2) then
                write(*,*) 'problems in predictorstep, SP(p,p).lt.0:', SP(pPart%momentum,pPart%momentum)
                write(*,*) 'gradR = ', gradR(1:3,iEns,iPart)
                call writeParticle_debug(pPart)
             end if
             if (isBaryon(pPart%ID)) then
                write(*,*) 'Deleting Particle'
                call protocolFailure(pPart,'predictor')
                pPart%ID=0
             end if
          end if

       end do
    end do

    if (get_useOffShellPotentialBaryons().or.get_useOffShellPotentialMesons()) call setMass(Parts)

  end subroutine predictorStep



  !****************************************************************************
  !****s* propagation/correctorStep
  ! NAME
  ! subroutine correctorStep(Parts,gradP,gradR,delta_T)
  !
  ! PURPOSE
  ! * Routine performs a corrector step on  a particle vector previously undergoing a simple Euler method.
  !
  ! INPUTS
  ! * real, intent(in) :: delta_T
  ! * real, dimension(:,:,:),intent(out) :: gradR, gradP ! Gradients in R and P-Space of predictor step
  ! * type(particle),intent(inOut),dimension(:,:) :: Parts -- particle vector at predicted values
  !
  ! OUTPUT
  ! * type(particle),intent(inOut),dimension(:,:) :: Parts -- particle vector at t+Delta(t)
  !
  ! NOTES
  ! * Used as  corrector step in a predictor/corrector scheme
  !****************************************************************************
  subroutine correctorStep(Parts,gradP,gradR,delta_T)
    use particleDefinition
    use densityModule, only: gridsize
    use inputGeneral, only: continousBoundaries

    real, dimension(1:,:,:),intent(in) :: gradP
    real, dimension(0:,:,:),intent(in) :: gradR
    type(particle),intent(inOut),dimension(:,:),target :: Parts

    real, intent(in) :: delta_T
    real, dimension(1:3) :: Grad_P_Predictor
    real, dimension(0:3) :: Grad_R_Predictor
    integer :: i
    integer :: iEns,iPart
    logical :: checkVelo_flag
    type(particle), pointer :: pPart

    ! (3) Corrector Step: Evaluate Gradients at predicted value

    do iEns=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do iPart=lbound(Parts,dim=2),ubound(Parts,dim=2)

          pPart => Parts(iEns,iPart)
          if (pPart%ID < 0) exit
          if (pPart%ID == 0) cycle

          call gradients(pPart,Grad_P_Predictor,Grad_R_Predictor) ! Evaluate dH/dr and dH/dp

          ! (3) Do time stepping considering also the predicted gradients
          !
          !       momentum(1:3,t+delta t)=momentum(1:3,t)-delta_T*(grad_R+grad_R_Predictor)/2.
          !       momentum_predicted(1:3,t+delta t)=momentum(1:3,t)-delta_T*grad_R
          !
          ! => momentum(1:3,t+delta t)=momentum_Predicted(1:3)-delta_T*(-grad_R+grad_R_Predictor)/2.
          ! And be reminded that "Parts" is the predictor (therefore the - sign in front of gradR and gradP).
          !
          ! The same argument also holds for the position.
          !
          pPart%momentum(1:3)=pPart%momentum(1:3)-delta_T*(-gradR(1:3,iEns,iPart)+grad_R_Predictor(1:3))*0.5
          pPart%momentum( 0 )=pPart%momentum( 0 )+delta_T*(-gradR( 0, iEns,iPart)+grad_R_Predictor( 0 ))*0.5
          pPart%position(1:3)=pPart%position(1:3)+delta_T*(-gradP(1:3,iEns,iPart)+grad_P_Predictor(1:3))*0.5

          ! (4) Save velocity
          pPart%velocity(1:3)=(gradP(1:3,iEns,iPart)+grad_P_Predictor(1:3))*0.5

          if (continousBoundaries) then
             do i=1,3
                if (pPart%position(i).gt.gridsize(i)) then
                   pPart%position(i)= pPart%position(i)-2*gridsize(i)
                else if (pPart%position(i).lt.-gridsize(i)) then
                   pPart%position(i)= pPart%position(i)+2*gridsize(i)
                end if
             end do
          end if

          if (tachyonDebug) checkVelo_flag=checkVelo(pPart)
       end do
    end do

  end subroutine correctorStep



  !****************************************************************************
  !****s* propagation/propagate_euler
  ! NAME
  ! subroutine propagate_euler(Parts,delta_T)
  ! PURPOSE
  ! * Routine propagates a particle vector
  ! INPUTS
  ! * type(particle),intent(inOUT),dimension(:,:)  :: Parts  -- particle vector which should be propagated
  ! * real, intent(in) :: delta_T ! time step (fm/c)
  ! OUTPUT
  ! * type(particle),intent(inOUT),dimension(:,:)  :: Parts  -- particle vector which should be propagated
  ! NOTES
  ! * Uses simple Euler time stepping. See also files in "Documentation_Extra/propagation/".
  !****************************************************************************
  subroutine propagate_euler(Parts,delta_T)

    use particleDefinition
    use nucleusDefinition
    use densityModule, only: gridsize
    use inputGeneral, only: continousBoundaries,eventtype
    use eventtypes, only: HeavyIon, hadron
    use nucleus, only: getTarget, getProjectile

    type(particle),intent(inOut),dimension(:,:),target :: Parts
    real, intent(in) :: delta_T

    type(tNucleus), pointer :: proj, targ
    integer :: iEns,iPart
    type(particle), pointer :: pPart

    real,dimension(0:3) :: grad_R
    real,dimension(1:3) :: grad_P
    integer :: i
    logical :: checkVelo_flag

    if (startFlag) call init

    if (.not.UseHadronic .and. (eventtype==HeavyIon .or. eventtype==hadron)) then
       proj => getProjectile()
       targ => getTarget()
    end if

    do iEns=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do iPart=lbound(Parts,dim=2),ubound(Parts,dim=2)

          pPart => Parts(iEns,iPart)
          if (pPart%ID < 0) exit
          if (pPart%ID == 0) cycle

          if (UseHadronic) then

             call gradients(pPart,Grad_P,Grad_R) ! Evaluate dH/dr and dH/dp
             pPart%momentum(1:3)=pPart%momentum(1:3)-delta_T*grad_R(1:3)
             pPart%momentum( 0 )=pPart%momentum( 0 )+delta_T*grad_R( 0 )
             pPart%position(1:3)=pPart%position(1:3)+delta_T*grad_P(1:3)
             pPart%velocity(1:3)=grad_P(1:3)
             if (tachyonDebug) checkVelo_Flag=checkVelo(pPart)

          else if (pPart%perturbative .or. pPart%event(1)>=1000000 &
               .or. (eventtype/=HeavyIon .and. eventtype/=hadron) ) then ! Cascade

             call gradients(pPart,Grad_P,Grad_R) ! Evaluate dH/dr and dH/dp
             pPart%position(1:3)=pPart%position(1:3)+delta_T*grad_P(1:3)
             pPart%velocity(1:3)=grad_P(1:3)


          else  ! "Frozen" cascade (presently for eventtypes HeavyIon and hadron only !)

             if (mod(pPart%event(1),2)==1) then

                if (pPart%ID.ne.1) then
                   write(*,*) pPart%ID,pPart%event
                   call Traceback('wrong particle propagated! (1)')
                end if

                pPart%position(1:3) = pPart%position(1:3) + delta_T * targ%velocity(1:3)
                pPart%velocity(1:3) = targ%velocity(1:3)

             else if ( mod(pPart%event(1),2)==0) then

                !if(pPart%ID.ne.1) then
                !   write(*,*)pPart%ID,pPart%event
                !   call Traceback('wrong particle propagated! (2)')
                !endif

                pPart%position(1:3) = pPart%position(1:3) + delta_T * proj%velocity(1:3)
                pPart%velocity(1:3) = proj%velocity(1:3)

             else

                call gradients(pPart,Grad_P,Grad_R) ! Evaluate dH/dr and dH/dp
                pPart%position(1:3)=pPart%position(1:3)+delta_T*grad_P(1:3)
                pPart%velocity(1:3)=grad_P(1:3)

             end if

          end if

          if (continousBoundaries) then
             do i=1,3
                if (pPart%position(i).gt.gridsize(i)) then
                   pPart%position(i)= pPart%position(i)-2*gridsize(i)
                else if (pPart%position(i).lt.-gridsize(i)) then
                   pPart%position(i)= pPart%position(i)+2*gridsize(i)
                end if
             end do
          end if

       end do
    end do

  end subroutine propagate_euler

  !****************************************************************************
  !****s* propagation/propagate_cascade
  ! NAME
  ! subroutine propagate_cascade(Parts,delta_T)
  ! PURPOSE
  ! Routine propagates particles in case of cascade mode.
  ! Useful also for doing initial step in the RMF mode.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: Parts
  !   -- particles which should be propagated
  ! * real :: delta_T                            -- time step (fm/c)
  ! OUTPUT
  ! * type(particle), dimension(:,:) :: Parts
  ! NOTES
  ! Straight line trajectories; no momentum change (mean fields switched off).
  !****************************************************************************
  subroutine propagate_cascade(Parts,delta_T)
    use nucleusDefinition
    use particleDefinition
    use inputGeneral, only: eventtype
    use eventtypes, only: HeavyIon
    use inputGeneral, only: eventtype, continousBoundaries
    use nucleus, only: getTarget, getProjectile
    use densityModule, only: gridsize

    type(particle),intent(inOut),dimension(:,:),target :: Parts
    real, intent(in) :: delta_T

    type(tNucleus), pointer :: proj, targ
    integer :: iEns,iPart, i
    type(particle), pointer :: pPart

    if (startFlag) call init

    if (eventtype==HeavyIon) then
      proj => getProjectile()
      targ => getTarget()
    end if

    do iEns=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do iPart=lbound(Parts,dim=2),ubound(Parts,dim=2)

          pPart => Parts(iEns,iPart)
          if (pPart%ID < 0) exit
          if (pPart%ID == 0) cycle

          if (eventtype==HeavyIon .and. pPart%event(1)<1000000) then

             ! In case of HIC the particle is assumed to be "frozen"
             ! until it collides with other particle, i.e.
             ! it propagates with the speed of either target or projectile.

             if (pPart%ID.ne.1) then
                write(*,*) pPart%ID,pPart%event
                call Traceback('wrong particle propagated! (1)')
             end if

             if (mod(pPart%event(1),2)==0) then
                pPart%velocity(1:3) = targ%velocity(1:3)
             else
                pPart%velocity(1:3) = proj%velocity(1:3)
             end if

          else

             pPart%velocity(1:3) = pPart%momentum(1:3) / FreeEnergy( pPart )

          end if

          pPart%position(1:3) = pPart%position(1:3) + delta_T * pPart%velocity(1:3)

          if (continousBoundaries) then
             do i=1,3
                if (pPart%position(i).gt.gridsize(i)) then
                   pPart%position(i)= pPart%position(i)-2*gridsize(i)
                else if (pPart%position(i).lt.-gridsize(i)) then
                   pPart%position(i)= pPart%position(i)+2*gridsize(i)
                end if
             end do
          end if

      end do
    end do


  end subroutine propagate_cascade

  !****************************************************************************
  !****f* propagation/checkVelo
  ! NAME
  ! logical function checkVelo(part)
  ! PURPOSE
  ! Checks the velocity of a given particle
  ! INPUTS
  ! * type(particle),intent(inout) :: part
  ! NOTES
  ! * Also does some statistical stuff
  !****************************************************************************
  logical function checkVelo(part)
    use IdTable, only: rho,omegaMeson,phi,photon
    use particleDefinition
    use minkowski, only: abs4
    use offshellPotential, only: treatParticleOffShell, get_offshell_debug, get_useOffShellPotentialMesons
    use hist2Df90
    use output
    use mediumDefinition
    use mediumModule, only: mediumAt

    type(particle),intent(inout) :: part

    logical, save :: first=.true.
    integer, save :: failures=0
    integer, save :: numtry=0
    type(histogram2D),save :: hist2D_rho,hist2D_omega,hist2D_phi
    type(medium) :: med

    checkVelo=.true.
    if (part%ID.le.0 .or. part%ID==photon) return
    numtry=numtry+1

    if (startFlag) call init

    if (first .and. get_offshell_debug() .and. get_useOffShellPotentialMesons()) then
       call createHist2D(hist2D_rho,'number of tachyons as function of mass and momentum',(/0.,0./),(/2.,2./),(/0.02,0.02/))
       call createHist2D(hist2D_omega,'number of tachyons as function of mass and momentum',(/0.,0./),(/2.,2./),(/0.02,0.02/))
       call createHist2D(hist2D_phi,'number of tachyons as function of mass and momentum',(/0.,0./),(/2.,2./),(/0.02,0.02/))
       first = .false.
    end if

    if (get_offshell_debug() .and. get_useOffShellPotentialMesons()) then
       select case (part%ID)
       case (rho)
          call AddHist2D(hist2D_rho,(/abs4(part%momentum),absMom(part)/),0.,1.)
       case (omegaMeson)
          call AddHist2D(hist2D_omega,(/abs4(part%momentum),absMom(part)/),0.,1.)
       case (phi)
          call AddHist2D(hist2D_phi,(/abs4(part%momentum),absMom(part)/),0.,1.)
       end select
    end if

    if (1. - Dot_Product(part%velocity(1:3),part%velocity(1:3)) .le. 0.) then

       failures=failures+1
       checkVelo=.false.

       med = mediumAt(part%position)

       write(*,*)
       write(*,'(A,G12.5)')'Problems in CheckVelo : v = ',sqrt( Dot_Product(part%velocity(1:3),part%velocity(1:3)))
       write(*,'(A,G12.5,A,G12.5,A,G12.5)') 'm = ',abs4(part%momentum),'; pabs = ',absMom(part), &
            '; rho = ',med%densityproton+med%densityneutron
       call WriteParticle(6,0,0,part)
       write(*,*) 'Number of failures :',failures,' (',float(failures)/float(numtry)*100., '%)'
       call protocolFailure(part,'checkVelo')

       if (get_offshell_debug() .and. get_useOffShellPotentialMesons()) then
          select case (part%ID)
          case (rho)
             call AddHist2D(hist2D_rho,(/abs4(part%momentum),absMom(part)/),1.)
          case (omegaMeson)
             call AddHist2D(hist2D_omega,(/abs4(part%momentum),absMom(part)/),1.)
          case (phi)
             call AddHist2D(hist2D_phi,(/abs4(part%momentum),absMom(part)/),1.)
          end select

          if (mod(failures,100)==0) then
             call writeHist2D_Gnuplot(hist2D_rho,44,file='Tachyons2D_rho.dat')
             call writeHist2D_Gnuplot(hist2D_omega,45,file='Tachyons2D_omega.dat')
             call writeHist2D_Gnuplot(hist2D_phi,46,file='Tachyons2D_phi.dat')
          end if
       end if

       if (treatParticleOffShell(part%Id,part%offshellparameter)) then
          !ACCEPTABLE ONLY FOR OFFSHELL TRANSPORT, OTHERWISE SERIOUS ERROR!!!!!!!!!!
          write(*,*) 'this particle is now deleted!'
          part%ID=0
       else
          write(*,*) 'checkvelo: stop'
          call traceback()
          stop
       end if
    end if

  end function checkVelo


  !****************************************************************************
  !****s* propagation/checkVelos
  ! NAME
  ! subroutine checkVelos(part)
  ! PURPOSE
  ! Checks the velocities
  ! INPUTS
  ! *type(particle),dimension(:,:),intent(inout) :: part
  !****************************************************************************
  subroutine checkVelos(part)
    use particleDefinition

    type(particle),dimension(:,:),intent(inout) :: part
    integer :: iEns,iPart
    logical :: dummy

    do iEns=lbound(part,dim=1),ubound(part,dim=1)
       do iPart=lbound(part,dim=2),ubound(part,dim=2)
          if (part(iEns,iPart)%ID<=0) cycle
          dummy=checkVelo(part(iEns,iPart))
       end do
    end do

  end subroutine checkVelos


  !****************************************************************************
  ! cf. Interface updateVelocity
  !****************************************************************************
  subroutine updateVelocity_matrix(Parts)
    use particleDefinition
    use output, only: DoPR

    real,dimension(1:3)  :: grad_P
    type(particle),intent(inOut), dimension(:,:) :: Parts
    integer :: i,j
    logical :: success

    if (startFlag) call init

    if (DoPr(2)) write(*,*) 'Updating particle velocities'

    ensLoop: do i=lbound(Parts, dim=1),ubound(Parts, dim=1)
       do j=lbound(Parts, dim=2),ubound(Parts, dim=2)
          if (Parts(i,j)%ID <= 0) cycle ensLoop
          call gradients(Parts(i,j),Grad_P) ! Evaluate dH/dp
          Parts(i,j)%velocity(1:3)=grad_P
          success=checkVelo(Parts(i,j))
       end do
    end do ensLoop
  end subroutine updateVelocity_matrix



  !****************************************************************************
  ! cf. Interface updateVelocity
  !****************************************************************************
  subroutine updateVelocity_field(Parts)
    use particleDefinition

    real,dimension(1:3)  :: grad_P
    type(particle),intent(inOut), dimension(:) :: Parts
    integer :: i
    logical :: success

    if (startFlag) call init

    do i=lbound(Parts, dim=1),ubound(Parts, dim=1)
       if (Parts(i)%ID <= 0) cycle
       call gradients(Parts(i),Grad_P) ! Evaluate dH/dp
       Parts(i)%velocity(1:3)=grad_P
       success=checkVelo(Parts(i))
    end do
  end subroutine updateVelocity_field

  !****************************************************************************
  ! cf. Interface updateVelocity
  !****************************************************************************
  subroutine updateVelocity_1(part,success)
    use particleDefinition
    use output, only: writeParticle_debug

    real,dimension(1:3)  :: grad_P
    type(particle),intent(inOut) :: part
    logical,intent(out),optional :: success
    logical :: c

    if (startFlag) call init

    if (part%ID <= 0) return
    call gradients(part,Grad_P) ! Evaluate dH/dp
    part%velocity(1:3)=grad_P

    c=checkVelo(part)

    if (present(success)) then
       success = c
       if (tachyonDebug.and..not.success) then
          call writeParticle_debug(part)
          stop
       end if
    end if
  end subroutine updateVelocity_1



  subroutine protocolFailure(part,kindOfFailure)
    use particleDefinition

    character(*),   intent(in) :: kindOfFailure
    type(particle), intent(in) :: part

    logical, save :: firstTime=.true.
    character(20) :: form
    integer, save :: numFailures=0
    form='(A40,4G13.5)'

    numFailures=numFailures+1

    if (firstTime) then
       open(22,file='propa_failures.txt')
       write(22,*)
       firstTime=.false.
    else
       open(22,file='propa_failures.txt',position='append')
    end if
    write(22,'(A)')'Problems in ', kindOfFailure
    write(22,form)'|v| = '                      , sqrt(Dot_Product(part%velocity(1:3),part%velocity(1:3)))
    write(22,form)'Particle ID: '               , part%ID
    write(22,form)'Particle mass: '             , part%mass
    write(22,form)'Particle offshellparameter: ', part%offshellparameter
    write(22,form)'Particle number: '           , part%number
    write(22,form)'Particle firstevent: '       , part%firstevent
    write(22,form)'Particle velocity: '         , part%velocity
    write(22,form)'Particle momentum: '         , part%momentum
    write(22,form)'Particle perturbative?: '    , part%perturbative
    write(22,form)'Particle charge: '           , part%Charge
    write(22,form)'Particle position: '         , part%position
    write(22,*)   'Number of failures: '        , numFailures
    write(22,*)
    close(22)
  end subroutine protocolFailure


end module propagation
