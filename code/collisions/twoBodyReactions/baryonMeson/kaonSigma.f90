!******************************************************************************
!****m* /kaonSigma_resonance
! NAME
! module kaonSigma_resonance
!
! PURPOSE
! Includes the cross sections for kaon-sigma scattering in the resonance regime
!
! Public routines:
! * kaonSigma
! NOTES
! Resonances are included into the model of Huang et al for calculating
! the cross sections. After this, the treatment is done as for 2 -> 2 reactions.
!******************************************************************************
module kaonSigma_resonance

  implicit none
  private

  ! Debug-flags
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! To decide wether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.

  public :: kaonSigma

contains

  !****************************************************************************
  !****s* kaonSigma_resonance/kaonSigma
  ! NAME
  ! subroutine kaonSigma (srts,teilchenIN,teilchenOUT,sigmaTot,sigmaElast,plotFlag)
  !
  ! PURPOSE
  ! Evaluates kaon sigma -> anything cross sections and returns also a "preevent"
  !
  ! INPUTS
  ! * real, intent(in)                           :: srts                  ! sqrt(s) in the process
  ! * type(particle),dimension(1:2), intent(in)  :: teilchenIn            ! colliding particles
  !
  ! Debugging:
  ! * logical, intent(in),optional               :: plotFlag              ! Switch on plotting of the  Xsections
  !
  ! RESULT
  ! * real, intent(out)                          :: sigmaTot         ! total Xsection
  ! * real, intent(out)                          :: sigmaElast       ! elastic Xsection
  !
  ! This routine does a Monte-Carlo-decision according to the partial cross sections to decide on a final state with
  ! maximal 3 final state particles. These are returned in the vector teilchenOut. The kinematics of these teilchen is
  ! only fixed in the case of a single produced resonance. Otherwise the kinematics still need to be established. The
  ! result is:
  ! * type(preEvent),dimension(1:3), intent(out)               :: teilchenOut     ! colliding particles
  !
  ! NOTES
  ! Possible final states are :
  ! * 2-particle : pi N, piDelta
  !****************************************************************************
  subroutine kaonSigma (srts,partIn,partOut,sigmaTot,sigmaElast,plotFlag)

    use idTable
    use particleDefinition
    use particleProperties, only: hadron
    use mediumDefinition
    use preEventDefinition, only: preEvent
    use twoBodyTools, only: velocity_correction, convertToAntiParticles, pcm,searchInInput
    use RMF, only: getRMF_flag

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Input
    real, intent(in)                              :: srts                  ! sqrt(s) in the process
    type(particle),dimension(1:2), intent(in)     :: partIn            ! colliding particles

    logical, intent(in),optional                  :: plotFlag              ! Switch on plotting of the  Xsections

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Output
    type(preEvent),dimension(1:3), intent(out) :: partOut      ! colliding particles
    real, intent(out)                          :: sigmaTot         ! total Xsection
    real, intent(out)                          :: sigmaElast       ! elastic Xsection


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Cross sections
    real,dimension(-1:1) ::     piN                ! -> pi N, index denotes pion charge
    real,dimension(-1:1) ::     piDelta            ! -> pi Delta, index denotes pion charge



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Local variables
    real :: fluxCorrector      ! Correction of the fluxfactor due to different velocities
                               ! in the medium compared to the vacuum
    type(particle) :: kaon_particle, sigma_particle
    logical :: antiParticleInput, failFlag

    antiParticleINPUT=.false. ! .true. if antiparticle in the input

    ! Initialize output
    partOut(:)%ID=0                    ! ID of produced particles
    partOut(:)%charge=0                ! Charge of produced particles
    partOut(:)%antiParticle=.false.    ! Whether produced particles are particles or antiparticles
    partOut(:)%mass=0                  ! Mass of produced particles

    ! (1) Check  Input
    call searchInInput(partIn,kaon,sigmaResonance,kaon_particle,sigma_particle,failFlag)
    if (failFlag) then
       write(*,*) 'Wrong input in KaonNuc', partIn%ID
    end if


    if (sigma_particle%antiParticle.and.kaon_particle%antiParticle) then
       ! Both are antiparticles: s=0 scattering channel
       !
       ! Invert all particles in antiparticles
       sigma_particle%Charge        =  -sigma_particle%Charge
       sigma_particle%antiparticle  = .false.
       kaon_particle%Charge          =  -kaon_particle%Charge
       kaon_particle%antiparticle  = .false.
       antiParticleInput=.true.
    else if ((.not.(sigma_particle%antiParticle)).and.(.not.(kaon_particle%antiParticle))) then
       ! Both are no antiparticles : s=0 scattering channel
       antiParticleInput=.false.
    else
       ! Not yet implemented: S=-2,-1,1,2 scattering channels
       sigmaTot=0.
       sigmaElast=0.
       return
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
    if (Sum(partOut(:)%Charge).ne.sigma_particle%charge+kaon_particle%charge) then
       write(*,*) 'No charge conservation in pionNuc!!! Critical error' ,kaon_particle%Charge, &
            & sigma_particle%Charge, partOut(:)%Charge,partOut(:)%ID
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
    !****s* kaonSigma/evaluateXsections
    ! NAME
    ! subroutine evaluateXsections
    !
    ! PURPOSE
    ! Evaluates kaon sigma -> anything cross sections
    !
    ! NOTES
    ! There are no resonance contributions to kaon Sigma scattering. The contributions
    ! are given by Huang parametrizations of pi N(Delta) -> Kaon Sigma
    ! using detailed balance relations.
    !**************************************************************************

    subroutine evaluateXsections

      use parametrizationsBarMes, only: huang, huangd
      use clebschGordan, only: clebschSquared
      use output, only: writeparticle
      use constants, only: mN, mPi

      real,dimension(1:4) :: sigmaHuang
      real :: sigmaHuangDelta

      real :: pFinal, pInitial, detailedBalanceFactor
      integer :: pionCharge
      real ::    isoZ_delta
      !########################################################################
      ! Evaluate partial cross sections
      !########################################################################



      !************************************************************************
      ! kaon Sigma -> pi N
      !************************************************************************
      ! piN = cross section by Huang (pi N-> kaon Sigma) by detailed balance

      pFinal   = pCM(srts,mPi,mN)
      pInitial = pCM(srts,kaon_particle%mass,sigma_PARTICLE%mass)

      if (pinitial.lt.1E-12) then
         write(*,*) 'WARNING: pInitial is zero in kaonSigma', pinitial
         write(*,*) 'kaon meson:'
         call writeparticle(6,0,0,kaon_Particle)
         write(*,*) 'Sigma:'
         call writeparticle(6,0,0,sigma_Particle)
         detailedBalanceFactor= 0.
      else
         detailedBalanceFactor= (pFinal/pInitial)**2
      end if

      sigmaHuang = huang(srts) * detailedBalanceFactor

      piN=0.
      if (kaon_particle%charge.eq.0) then
         select case (sigma_particle%charge)
         case (-1)
            piN(-1)=sigmaHuang(1)
            piN(0) =0.
            piN(1) =0.
         case (0)
            piN(-1)=sigmaHuang(3)
            piN(0) =sigmaHuang(2)
            piN(1) =0.
         case (1)
            piN(-1)=0.
            piN(0) =sigmaHuang(3)
            piN(1) =sigmaHuang(4)
         end select
      else if (kaon_particle%charge.eq.1) then
         select case (sigma_particle%charge)
         case (-1)
            piN(-1)=sigmaHuang(4)
            piN(0) =sigmaHuang(3)
            piN(1) =0.
         case (0)
            piN(-1)=0.
            piN(0) =sigmaHuang(2)
            piN(1) =sigmaHuang(3)
         case (1)
            piN(-1)=0.
            piN(0) =0.
            piN(1) =sigmaHuang(1)
         end select
      else
         write(*,*) 'error in KaonSigma', kaon_particle%charge
      end if

      !************************************************************************
      ! kaon Sigma -> pi Delta
      !************************************************************************
      piDelta=0.

      ! divide out Clebsch Gordan coefficient for (pi- D++ -> S0 k) and multiply by initial clebsch Gordan
      sigmaHuangDelta = 6.*huangd(srts)*clebschSquared(0.5,1.,0.5,real(kaon_particle%charge)-0.5,real(sigma_particle%charge))
      ! Detailed balance
      pFinal   = pcm(srts,mPi,hadron(delta)%mass)
      pInitial = pcm(srts,kaon_particle%mass,sigma_particle%mass)


      if (pinitial.lt.1E-12) then
         write(*,*) 'WARNING: pInitial is zero in kaonSigma', pinitial
         write(*,*) 'kaon meson:'
         call writeparticle(6,0,0,kaon_Particle)
         write(*,*) 'Sigma:'
         call writeparticle(6,0,0,sigma_Particle)
         detailedBalanceFactor= 0.
      else
         detailedBalanceFactor= (pFinal/pInitial)**2
      end if
      ! factor of 2 due to different spin in initial and final state
      sigmaHuangDelta=sigmaHuangDelta*detailedBalanceFactor

      do pionCharge=-1,1
          isoZ_delta=kaon_particle%charge+sigma_particle%charge-pionCharge-0.5
          if (abs(isoZ_delta).lt.1.51) then
             piDelta(pionCharge)=sigmaHuangDelta*clebschSquared(1.5,1.,0.5,isoZ_delta,real(pionCharge))
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
         piDelta=piDelta*fluxcorrector
         piN=piN*fluxcorrector
      end if

      !########################################################################
      ! Sum up everything for the total cross section
      !########################################################################

      sigmaTot=sum( piN ) + sum(piDelta)

    end subroutine evaluateXsections




    !**************************************************************************
    !****s* kaonSigma/makeDecision
    ! NAME
    ! subroutine makeDecision
    !
    ! PURPOSE
    ! ...
    !
    !**************************************************************************

    subroutine makeDecision
      use random, only: rn

      real :: cut
      integer :: totalCharge, pionCharge

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=kaon_particle%Charge+sigma_particle%Charge

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

      ! pi Delta production
      do pionCharge=-1,1
         if (piDelta(pionCharge).ge.cut) then
            partOut(1)%Id=pion
            partOut(2)%Id=Delta

            partOut(1)%Charge=pionCharge
            partOut(2)%Charge=totalCharge-pionCharge
            return
         end if
         cut=cut-piDelta(pionCharge)
      end do

      ! Not event was generated:
      write(*,*) 'Error in makedecision of kaonNuc', piDelta, cut
      stop

    end subroutine makeDecision

    !**************************************************************************
    !****s* kaonSigma/makeOutput
    ! NAME
    ! subroutine makeOutput
    !
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV]
    ! .
    ! Filenames:
    ! * 'kaonSigma_sigTotElast.dat'        : sigmaTot, sigmaElast
    ! * 'kaonSigma_nonStrange.dat'         : non-strange meson in final state
    !**************************************************************************
    subroutine makeOutPut

      logical, save :: initFlag=.true.

      ! The output files
      character(30), dimension(1:6) :: outputFile
      real :: plab


      outputFile(1)='kaonSigma_sigTotElast.dat'
      outputFile(2)='kaonSigma_nonStrange.dat'

      plab=SQRT(((srts**2-kaon_particle%mass**2-sigma_particle%mass**2)/2./sigma_particle%mass)**2-kaon_particle%mass**2)


      if (initFlag) then
         open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         write(101,*) '# srts, plab, sigmaTot, sigmaElast '
         write(102,*) '# srts, plab, piN(-1:1), piDelta(-1:1) '
         initFlag=.false.
      else
         open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
      end if
      write(101,'(4F9.3)') srts, plab,sigmaTot, sigmaElast
      write(102,'(8F9.3)') srts, plab,piN, piDelta
      close(101)
      close(102)
    end subroutine makeOutPut
  end subroutine kaonSigma


end module kaonsigma_resonance
