!******************************************************************************
!****m* /pionP11_1440_resonance
! NAME
! module pionP11_1440_resonance
! PURPOSE
! Includes the cross sections for pion-P11_1440 scattering in the resonance regime.
! Implemented are the following reactions:
! * pion P11_1440 -> X
! Public routines:
! * pionNuc
!******************************************************************************
module pionP11_1440_resonance
  implicit none
  private

  ! Debug-flags
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! To decide wether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.

  public :: pionP11_1440

contains


  !****************************************************************************
  !****s* pionP11_1440_resonance/pionNuc
  ! NAME
  !  subroutine pionNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,plotFlag)
  !
  ! PURPOSE
  ! Evaluates pion P11_1440 -> anything cross sections and returns also a "preevent"
  !
  ! RESULT
  ! * real, intent(out)                                        :: sigmaTot         ! total Xsection
  ! * real, intent(out)                                        :: sigmaElast       ! elastic Xsection
  !
  ! This routine does a Monte-Carlo-decision according to the partial cross sections to decide on a final state with
  ! maximal 3 final state particles. These are returned in the vector teilchenOut. The kinematics of these teilchen is
  ! only fixed in the case of a single produced resonance. Otherwise the kinematics still need to be established. The
  ! result is:
  ! * type(preEvent),dimension(1:3), intent(out)               :: teilchenOut     ! colliding particles
  !
  ! NOTES
  ! Possible final states are :
  ! * 1-particle : baryon Resonances
  !****************************************************************************

  subroutine pionP11_1440 (srts,partIn,mediumAtColl,momLRF,partOut,sigmaTot,sigmaElast,plotFlag)

    use idTable
    use particleDefinition
    use mediumDefinition
    use preEventDefinition, only: preEvent
    use twoBodyTools, only: velocity_correction, convertToAntiParticles,searchInInput
    use RMF, only: getRMF_flag

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Input
    real, intent(in)                              :: srts                  ! sqrt(s) in the process
    type(particle),dimension(1:2), intent(in)     :: partIn            ! colliding particles
    type(medium), intent(in)                      :: mediumAtColl     ! Medium informations at the position of the collision
    logical, intent(in),optional                  :: plotFlag              ! Switch on plotting of the  Xsections
    real, intent(in) ,dimension(0:3)              :: momLRF           ! Total Momentum in LRF
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Output
    type(preEvent),dimension(1:3), intent(out) :: partOut      ! colliding particles
    real, intent(out)                          :: sigmaTot         ! total Xsection
    real, intent(out)                          :: sigmaElast       ! elastic Xsection

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Cross sections
    real, dimension(Delta:nbar) :: sigmaRes      ! partial cross sections for pi P1440 -> R
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Field to store the resonance masses
    real, dimension(Delta:nbar) :: massRes       !  Resonance masses
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Local variables
    real :: fluxCorrector      ! Correction of the fluxfactor due to different velocities
                               ! in the medium compared to the vacuum
    type(particle) :: partPion, P11_1440_particle
    logical :: antiParticleInput, failFlag


    antiParticleINPUT=.false. ! .true. if antiparticle in the input

    ! Initialize output
    partOut(:)%ID=0                    ! ID of produced particles
    partOut(:)%charge=0                ! Charge of produced particles
    partOut(:)%antiParticle=.false.    ! Whether produced particles are particles or antiparticles
    partOut(:)%mass=0                  ! Mass of produced particles

    ! (1) Check  Input
    call searchInInput(partIn,pion,P11_1440,partPion,P11_1440_particle,failFlag)
    if (failFlag) then
       write(*,*) 'Wrong input in PionNuc', partIn%ID
    end if

    if (abs(partPion%charge).gt.1) write(*,*) 'wrong pion charge in pionNuc', partPion%charge

    if (partPion%antiParticle) then
       ! This case is not considered yet
       write(*,*) 'pion is antiparticle in "pionNuc"!!!',partIn%ID,partIn%antiparticle
       stop
    end if

    if (P11_1440_particle%antiParticle) then
       ! Invert all particles in antiparticles
       P11_1440_particle%Charge        =  -P11_1440_particle%Charge
       P11_1440_particle%antiparticle  = .false.
       partPion%Charge          =  -partPion%Charge
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


    ! Cutoff to kick those case out, that the cross section is zero
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
    if (Sum(partOut(:)%Charge).ne.P11_1440_particle%charge+partPion%charge) then
       write(*,*) 'No charge conservation in pionNuc!!! Critical error' ,partPion%Charge, &
            & P11_1440_particle%Charge, partOut(:)%Charge,partOut(:)%ID
       stop
    end if

    ! (6) Invert particles in antiParticles if input included antiparticles
    if (antiParticleInput) then
       if (debugFlagAnti) write(*,*) partOut
       call convertToAntiParticles(partOut)
       if (debugFlagAnti) write(*,*) partOut
    end if

  contains

    subroutine evaluateXsections
      use resonanceCrossSections, only: barMes2resonance, barMes_R_barMes
      use idTable, only: P11_1440, pion, pion

      real, dimension(1:3) ::  position
      logical :: perturbative
      !real :: sigmaTotal_HE,sigmaElast_HE

      position=0.5*(partIn(1)%position+partIn(2)%position)
      if (partIn(1)%perturbative.or.partIn(2)%perturbative) then
         perturbative=.true.
      else
         perturbative=.false.
      end if

      !########################################################################
      ! Evaluate partial cross sections
      !########################################################################

     !*************************************************************************
      ! pion P11-1440 -> R
      !************************************************************************

      ! Full resonance contribution in the medium
      sigmaRes = barMes2resonance (pion,P11_1440,partPion%charge,P11_1440_particle%charge,.true., &
                                   mediumAtColl,momLRF,massRes, &
                                   partPion%Mass,P11_1440_particle%Mass,position,perturbative,srts)

      !########################################################################
      ! evaluate elastic Xsection
      !########################################################################

      sigmaElast=barMes_R_barMes(pion,P11_1440,pion,P11_1440,&
           & partPion%Charge,P11_1440_particle%Charge,partPion%Charge,P11_1440_particle%Charge, &
           & .false.,.false.,mediumAtColl,momLRF,&
           & partPion%Mass,P11_1440_particle%Mass,position,perturbative,srts)


      !########################################################################
      ! Do the flux correction for each channel
      !########################################################################

      if (fluxCorrector_flag) then
         ! We do this for each channel since they might show up seperately in the output if makeoutput is called
         sigmaElast=sigmaElast*fluxcorrector
         sigmaRes=sigmaRes *fluxcorrector
      end if

      !########################################################################
      ! Sum up everything for the total cross section
      !########################################################################
      ! Be careful since sigma elast is already included in the partial cross sections, therefore it is not
      ! included in the total cross section

      sigmaTot=sum (sigmaRes )

    end subroutine evaluateXsections


    subroutine makeDecision
      use random, only: rn

      real :: summe, cut, cut2
      integer :: resID, totalCharge

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=partPion%Charge+P11_1440_particle%Charge
      !########################################################################
      ! (1) Resonance production
      !########################################################################
      if (sum(sigmaRes)>=cut) then
         summe=0.
         cut2=rn()*sum(sigmaRes)
         do resId=Delta,nbar
            summe=summe+sigmaRes(resID)
            if (summe>=cut2) exit
         end do
         partOut(1)%Id=resID
         partOut(1)%Charge=totalCharge
         partOut(1)%Mass=massRes(resID)
         return
      end if
      cut=cut-sum(sigmaRes)


      ! Not event was generated:
      write(*,*) 'Error in makedecision of pionNuc', sum(sigmaRes) , cut
      stop

    end subroutine makeDecision


    !**************************************************************************
    !****s* pionP11_1440_resonance/makeOutput
    ! NAME
    ! subroutine makeOutput
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV].
    ! Filenames:
    ! * 'pionP11_1440_sigTotElast.dat'    : sigmaTot, sigmaElast
    ! * 'pionP11_1440_resProd.dat'        : Baryon resonance production
    !**************************************************************************
    subroutine makeOutPut
      logical, save :: initFlag=.true.
      real :: plab
      character(len=30), parameter :: outputFile(1:2) = (/ 'pionP11_1440_sigTotElast.dat', &
                                                           'pionP11_1440_resProd.dat    ' /)

      plab=SQRT(((srts**2-partPion%mass**2-P11_1440_particle%mass**2)/2./P11_1440_particle%mass)**2-partPion%mass**2)

      if (initFlag) then
         open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         write(101,*) '# srts, plab, sigmaTot, sigmaElast '
         write(102,*) '# srts, plab, sigmaRes(2:40)'
         initFlag=.false.
      else
         open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
      end if
      write(101,'(4F9.3)') srts, plab,sigmaTot, sigmaElast
      write(102,'(41F9.3)') srts, plab, sigmaRes(2:40)
      close(101)
      close(102)

    end subroutine makeOutPut

  end subroutine pionP11_1440


end module pionP11_1440_resonance
