!******************************************************************************
!****m* /etaDelta_resonance
! NAME
! module etaDelta_resonance
! PURPOSE
! Includes the cross sections for eta-delta scattering in the resonance regime
!
! Public routines:
! * etaDelta
!******************************************************************************
module etaDelta_resonance
  implicit none
  private

  ! Debug-flags
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! To decide wether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.

  public :: etaDelta

contains

  !****************************************************************************
  !****s* etaDelta_resonance/etaDelta
  ! NAME
  ! subroutine etaDelta (srts, teilchenIn, teilchenOut, sigmaTot, sigmaElast, plotFlag)
  !
  ! PURPOSE
  ! Evaluates eta Delta_resonance -> anything cross sections and returns also a "preevent"
  !
  ! INPUTS
  ! * real, intent(in)                          :: srts             ! sqrt(s) in the process
  ! * type(particle),dimension(1:2), intent(in) :: teilchenIn       ! colliding particles
  !
  ! Debugging:
  ! * logical, intent(in),optional              :: plotFlag         ! Switch on plotting of the  Xsections
  !
  ! RESULT
  ! * real, intent(out)                         :: sigmaTot         ! total Xsection
  ! * real, intent(out)                         :: sigmaElast       ! elastic Xsection
  !
  ! This routine does a Monte-Carlo-decision according to the partial cross sections to decide on a final state with
  ! maximal 3 final state particles. These are returned in the vector teilchenOut. The kinematics of these teilchen is
  ! only fixed in the case of a single produced resonance. Otherwise the kinematics still need to be established. The
  ! result is:
  ! * type(preEvent),dimension(1:3), intent(out)               :: teilchenOut     ! colliding particles
  !
  ! NOTES
  ! Possible final states are :
  ! * 2-particle : pi N
  !****************************************************************************
  subroutine etaDelta (srts, partIn, partOut, sigmaTot, sigmaElast, plotFlag)

    use idTable
    use particleDefinition
    use mediumDefinition
    use preEventDefinition, only: preEvent
    use twoBodyTools, only: velocity_correction, convertToAntiParticles, pcm,searchInInput
    use RMF, only: getRMF_flag
    use constants, only: mN, mPi

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Input
    real, intent(in)                              :: srts                  ! sqrt(s) in the process
    type(particle),dimension(1:2), intent(in)     :: partIn            ! colliding particles
    logical, intent(in),optional                  :: plotFlag              ! Switch on plotting of the  Xsections
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Output
    type(preEvent),dimension(1:3), intent(out) :: partOut              ! colliding particles
    real, intent(out)                          :: sigmaTot, sigmaElast     ! cross sections

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Local variables
    real,dimension(-1:1) :: piN            ! cross sections into pi N, index denotes pion charge
    real :: fluxCorrector      ! Correction of the fluxfactor due to different velocities
                               ! in the medium compared to the vacuum
    type(particle) :: eta_particle, delta_particle
    logical :: antiParticleInput, failFlag

    antiParticleINPUT=.false. ! .true. if antiparticle in the input

    ! Initialize output
    partOut(:)%ID=0                    ! ID of produced particles
    partOut(:)%charge=0                ! Charge of produced particles
    partOut(:)%antiParticle=.false.    ! Whether produced particles are particles or antiparticles
    partOut(:)%mass=0                  ! Mass of produced particles

    ! (1) Check  Input
    call searchInInput(partIn,eta,delta,eta_particle,delta_particle,failFlag)
    if (failFlag) then
       write(*,*) 'Wrong input in EtaDelta', partIn%ID
    end if

    if (eta_particle%antiParticle) then
       ! This case is not considered yet
       write(*,*) 'eta is antiparticle in "etaNuc"!!!',partIn%ID,partIn%antiparticle
       stop
    end if

    if (delta_particle%antiParticle) then
       ! Invert all particles in antiparticles
       delta_particle%Charge        =  -delta_particle%Charge
       delta_particle%antiparticle  = .false.
       eta_particle%Charge          =  -eta_particle%Charge
       antiParticleInput=.true.
    else
       antiParticleInput=.false.
    end if

    ! Correction of the fluxfactor due to different velocities in the medium compared to the vacuum
    if ( .not.getRMF_flag() ) then
      fluxCorrector=velocity_correction(partIn)
    else
      fluxCorrector=1.
    end if

    ! (2) Evaluate the cross sections
    call evaluateXsections

    ! Cutoff to kick the case out, that the cross section is zero
    if (sigmaTot.lt.1E-12) then
       sigmatot=0.
       sigmaElast=0.
       return
    end if

    ! (3) Plot them if wished
    if (present(PlotFlag).or.debugFlag) then
       if (plotFlag.or.debugFlag)  call makeOutput
    end if

    ! (4) Define final state
    call MakeDecision

    ! (5) Check Output
    if (Sum(partOut(:)%Charge).ne.delta_particle%charge+eta_particle%charge) then
       write(*,*) 'No charge conservation in pionNuc!!! Critical error' ,eta_particle%Charge, &
            & delta_particle%Charge, partOut(:)%Charge,partOut(:)%ID
       stop
    end if

    ! (6) Invert particles in antiParticles if input included antiparticles
    if (antiParticleInput) then
       if (debugFlagAnti) write(*,*) partOut
       call convertToAntiParticles(partOut)
       if (debugFlagAnti) write(*,*) partOut
    end if

  contains

    !**************************************************************************
    !****s* etaDelta/evaluateXsections
    ! NAME
    ! subroutine evaluateXsections
    !
    ! PURPOSE
    ! Evaluates eta Delta_resonance -> anything cross sections
    !
    ! NOTES
    ! There are no resonance contributions to eta Delta scattering.
    !**************************************************************************
    subroutine evaluateXsections
      use pionNucleon, only: matrixDeltaEta  ! !Matrix Element for pi N <-> eta Delta
      use clebschGordan, only: clebschSquared
      use output, only: writeparticle

      real :: pFinal, pInitial, isoz_nuk, piN_total
      integer :: pionCharge

      !************************************************************************
      ! eta Delta -> pi N
      !************************************************************************

      pFinal=pcm(srts,mPi,mN)
      pInitial=pcm(srts,eta_particle%mass,delta_particle%mass)

      if (pinitial.lt.1E-12) then
         write(*,*) 'WARNING: pInitial is zero in etaDelta', pinitial
         write(*,*) 'eta meson:'
         call writeparticle(6,0,0,eta_Particle)
         write(*,*) 'delta:'
         call writeparticle(6,0,0,delta_Particle)
         piN_total=0.
      else
         ! given by detailed balance: factor 1/2 due to (2j+1)-Terms in cross section and different
         ! spins in initial and final state
         piN_total= 0.5*matrixDeltaEta/srts**2/pinitial*pfinal
      end if

      piN=0.
      do pionCharge=-1,1
         ! Evaluate z-Component of nucleon isospin
         isoZ_nuk=delta_particle%charge-pionCharge-0.5
         if (abs(abs(isoZ_nuk)-0.5).lt.0.0001) then
            piN(pionCharge)=piN_total*clebschSquared(1.,0.5,1.5,real(pionCharge),isoZ_nuk)
         end if
      end do

      !########################################################################
      ! evaluate elastic Xsection
      !########################################################################

      sigmaElast=0.

      !########################################################################
      ! Do the flux correction for each channel
      !########################################################################

      if (fluxCorrector_flag) then
         ! We do this for each channel since they might show up seperately in the output if makeoutput is called
         sigmaElast=sigmaElast*fluxcorrector
         piN=piN*fluxcorrector
      end if

      !########################################################################
      ! Sum up everything for the total cross section
      !########################################################################

      sigmaTot=sum( piN )

    end subroutine evaluateXsections


    !333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
    !333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
    !333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333


    subroutine makeDecision
      use random, only: rn

      real :: cut!,cut2
      integer :: totalCharge
      integer :: pionCharge

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=eta_particle%Charge+delta_particle%Charge

      !########################################################################
      ! (1) Two -body final states
      !########################################################################

      ! piN production
      do pionCharge=-1,1
         if (piN(pionCharge).ge.cut) then
            partOut(1)%Id=pion
            partOut(2)%Id=nucleon

            partOut(1)%Charge=pionCharge
            partOut(2)%Charge=totalCharge-pionCharge
            return
         end if
         cut=cut-piN(pionCharge)
      end do

      ! Not event was generated:
      write(*,*) 'Error in makedecision of etaNuc', piN, cut
      stop

    end subroutine makeDecision


    !**************************************************************************
    !****s* etaDelta/makeOutput
    ! NAME
    ! subroutine makeOutput
    !
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV]
    ! .
    ! Filenames:
    ! * 'etaDelta_sigTotElast.dat'        : sigmaTot, sigmaElast
    ! * 'etaDelta_nonStrange_nuk.dat'     : non-strange meson with nucleon in final state
    !**************************************************************************
    subroutine makeOutPut()

      logical, save :: initFlag=.true.
      real :: plab
      character(len=27), parameter :: outputFile(1:2) = (/ 'etaDelta_sigTotElast.dat   ', 'etaDelta_nonStrange_nuk.dat' /)

      plab=SQRT(((srts**2-eta_particle%mass**2-delta_particle%mass**2)/2./delta_particle%mass)**2-eta_particle%mass**2)

      if (initFlag) then
         open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         write(101,*) '# srts, plab, sigmaTot, sigmaElast '
         write(102,*) '# srts, plab, piN(-1:1), etaN, pipiN '
         initFlag=.false.
      else
         open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
      end if

      write(101,'(4F9.3)') srts, plab, sigmaTot, sigmaElast
      write(102,'(5F9.3)') srts, plab, piN

      close(101)
      close(102)

    end subroutine makeOutPut

  end subroutine etaDelta


end module etaDelta_resonance
