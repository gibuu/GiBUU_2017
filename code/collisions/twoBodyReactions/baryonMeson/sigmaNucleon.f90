!******************************************************************************
!****m* /sigmaNucleon
! NAME
! module sigmaNucleon
! PURPOSE
! Includes the cross sections for sigma-nucleon scattering in the resonance regime.
! Public routines:
! * sigmaNuc
!******************************************************************************
module sigmaNucleon

  implicit none
  private

  ! Debug-flags
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! To decide wether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.

  public :: sigmaNuc

contains


  !****************************************************************************
  !****s* sigmaNucleon/sigmaNuc
  ! NAME
  ! subroutine sigmaNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,useHiEnergy,HiEnergySchwelle,plotFlag)
  !
  ! PURPOSE
  ! Evaluates sigma Nucleon -> anything cross sections and returns also a "preevent"
  !
  ! INPUTS
  ! * real, intent(in)                              :: srts                  ! sqrt(s) in the process
  ! * type(particle),dimension(1:2), intent(in)     :: teilchenIn            ! colliding particles
  ! * type(medium), intent(in)                      :: mediumATcollision     ! Medium informations at the position of the collision
  ! * real, intent(in) ,dimension(0:3)              :: momentumLRF           ! Total Momentum in LRF
  !
  ! High energy matching:
  ! * logical,intent(in)                            :: useHiEnergy
  ! * .true. if High-Energy cross sections are given by paramBarMesHE
  ! * real,intent(in)                               :: HiEnergySchwelle
  ! * threshold sqrt(s) for paramBarMesHE, i.e. at which energy the cross sections of paramBarMesHE are used
  !
  ! Debugging:
  ! * logical, intent(in),optional                  :: plotFlag              ! Switch on plotting of the  Xsections
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
  ! The cross sections are based upon a parametrization by Golubeva. See routine golub in
  ! parametrizationBarMes.
  ! NOTES
  ! Possible final states are :
  ! * 1-particle : baryon Resonances
  ! * 2-particle : pi N, sigma N
  !****************************************************************************
  subroutine sigmaNuc (srts,partIn,mediumAtColl,momLRF,partOut,sigmaTot,sigmaElast,&
                       useHiEnergy,HiEnergySchwelle,plotFlag)

    use idTable
    use particleDefinition
    use mediumDefinition
    use preEventDefinition, only: preEvent
    use twoBodyTools, only: velocity_correction, convertToAntiParticles, searchInInput
    use RMF, only: getRMF_flag
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Input
    real, intent(in)                              :: srts                  ! sqrt(s) in the process
    type(particle),dimension(1:2), intent(in)     :: partIn            ! colliding particles
    type(medium), intent(in)                      :: mediumAtColl     ! Medium informations at the position of the collision
    logical, intent(in),optional                  :: plotFlag              ! Switch on plotting of the  Xsections
    real, intent(in) ,dimension(0:3)              :: momLRF           ! Total Momentum in LRF
    logical,intent(in)                            :: useHiEnergy            ! .true. if High-Energy cross sections are given by paramBarMesHE
    real,intent(in)                               :: HiEnergySchwelle      ! threshold sqrt(s) for paramBarMesHE
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Output
    type(preEvent),dimension(1:3), intent(out) :: partOut      ! colliding particles
    real, intent(out)                          :: sigmaTot         ! total Xsection
    real, intent(out)                          :: sigmaElast       ! elastic Xsection

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Cross sections
    real, dimension(Delta:nbar) :: sigmaRes    ! sigma N -> R crosssection
    real,dimension(-1:1) ::     piN            ! -> pi N, index denotes pion charge
    real ::     sigmaN                         ! -> sigma N
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Field to store the resonance masses
    real, dimension(Delta:nbar) :: massRes     !  Resonance masses
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Local variables
    real :: fluxCorrector        ! Correction of the fluxfactor due to different velocities
                                 ! in the medium compared to the vacuum
    type(particle) :: sigma_particle, partNucl
    logical :: antiParticleInput, failFlag

    ! Field to store the resonance masses
          ! Resonance masses

    antiParticleINPUT=.false. ! .true. if antiparticle in the input

    ! Initialize output
    partOut(:)%ID=0                    ! ID of produced particles
    partOut(:)%charge=0                ! Charge of produced particles
    partOut(:)%antiParticle=.false.    ! Whether produced particles are particles or antiparticles
    partOut(:)%mass=0                  ! Mass of produced particles

    ! (1) Check  Input
    call searchInInput(partIn,sigmaMeson,nucleon,sigma_particle,partNucl,failFlag)
    if (failFlag) then
       write(*,*) 'Wrong input in SigmaNuc', partIn%ID
    end if

    if (sigma_particle%antiParticle) then
       ! This case is not considered yet
       write(*,*) 'sigma is antiparticle in "sigmaNuc"!!!',partIn%ID,partIn%antiparticle
       stop
    end if

    if (partNucl%antiParticle) then
       ! Invert all particles in antiparticles
       partNucl%Charge        =  -partNucl%Charge
       partNucl%antiparticle  = .false.
       sigma_particle%Charge          =  -sigma_particle%Charge
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
    if (Sum(partOut(:)%Charge).ne.partNucl%charge+sigma_particle%charge) then
       write(*,*) 'No charge conservation in pionNuc!!! Critical error' ,sigma_particle%Charge, &
            & partNucl%Charge, partOut(:)%Charge,partOut(:)%ID
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
      use resonanceCrossSections, only: barMes_R_barMes, barMes2resonance
      use mediumDefinition, only: vacuum
      !use parametrizationsBarMes, only : golub
      use idTable, only: nucleon, sigmaMeson
      use parBarMes_HighEnergy, only: paramBarMesHE
      use clebschGordan, only: clebschSquared
      use particleProperties, only: hadron
      use constants, only: mN, mPi

      real, dimension(1:3) ::  position
      logical :: perturbative
      !real, dimension(1:4) :: sigmaGolub
      !real, dimension(-1:1) :: piN_Golub
      !real :: pFinal,pInitial,detailedBalanceFactor
      integer :: pionCharge, nucCharge
      real :: elastic_Vacuum
      real, dimension(0:3) :: momentum_vacuum

      real :: sigmaTotal_HE,sigmaElast_HE ! High energy matchin

      sigmaTotal_HE = 0.
      sigmaElast_HE = 0.

      position=0.5*(partIn(1)%position+partIn(2)%position)
      if (partIn(1)%perturbative.or.partIn(2)%perturbative) then
         perturbative=.true.
      else
         perturbative=.false.
      end if

      momentum_vacuum(1:3)=partIn(1)%momentum(1:3)+partIn(2)%momentum(1:3)
      momentum_vacuum(0)=FreeEnergy(partIn(1))+FreeEnergy(partIn(2))

      !########################################################################
      ! Evaluate partial cross sections
      !########################################################################


      !call golub(srts,1,sigmaGolub) ! Results by Golubeva


      !************************************************************************
      ! sigma N -> pi N
      !************************************************************************

      ! Full resonance contribution in the vacuum
      sigmaRes = barMes2resonance (sigmaMeson,nucleon,sigma_particle%charge,partNucl%charge,.true.,vacuum, &
                                   momentum_vacuum,massRes,sigma_particle%Mass,partNucl%Mass,position,perturbative,srts)

      ! High energy matching
      piN=0.
      if (srts > mN + mPi) then
         if (useHiEnergy) then
            call paramBarMesHE(HiEnergySchwelle,sigmaMeson,nucleon,sigma_particle%charge,partNucl%charge,&
                 & mediumAtColl,sigmaTotal_HE,sigmaElast_HE)
            do pionCharge=-1,1
               nucCharge=sigma_particle%charge+partNucl%charge-pionCharge
               if ((nucCharge.eq.0).or.(nucCharge.eq.1)) then
                  piN(pionCharge)=Max(sigmaTotal_HE-sum(sigmaRes),0.)*clebschSquared(1.,0.5,0.5, &
                       & real(pionCharge),real(nucCharge)-0.5)
               end if
            end do
         end if
      end if


      !************************************************************************
      ! sigma N -> sigma N
      !************************************************************************


      elastic_vacuum=barMes_R_barMes(sigmaMeson,nucleon,sigmaMEson,nucleon,&
           & sigma_particle%Charge,partNucl%Charge,sigma_particle%Charge,partNucl%Charge,  &
           & .false.,.false.,Vacuum,momentum_vacuum,sigma_particle%Mass,partNucl%Mass, &
           & position,perturbative,srts)

      if (srts > mN + hadron(omegaMeson)%minmass) then
        sigmaN=Max(0., sigmaElast_HE - elastic_vacuum)
      else
        sigmaN=0.
      end if


      !************************************************************************
      ! sigma N -> R
      !************************************************************************

      ! Full resonance contribution in the medium
      sigmaRes = barMes2resonance(sigmaMeson,nucleon,sigma_particle%charge,partNucl%charge,.true., &
                                  mediumAtColl,momLRF,massRes, &
                                  sigma_particle%Mass,partNucl%Mass,position,perturbative,srts)

      !########################################################################
      ! evaluate elastic Xsection
      !########################################################################

      sigmaElast=barMes_R_barMes(sigmaMeson,nucleon,sigmaMeson,nucleon,&
           & sigma_particle%Charge,partNucl%Charge,sigma_particle%Charge,partNucl%Charge, &
           & .false.,.false.,mediumAtColl,momLRF,&
           & sigma_particle%Mass,partNucl%Mass,position,perturbative,srts)+sigmaN


      !########################################################################
      ! Do the flux correction for each channel
      !########################################################################

      if (fluxCorrector_flag) then
         ! We do this for each channel since they might show up seperately in the output if makeoutput is called
         sigmaElast=sigmaElast*fluxcorrector
         sigmaN=sigmaN*fluxcorrector
         piN=piN*fluxcorrector
         sigmaRes=sigmaRes *fluxcorrector
      end if

      !########################################################################
      ! Sum up everything for the total cross section
      !########################################################################
      ! Be careful since sigma elast is already included in the partial cross sections, therefore it is not
      ! included in the total cross section

      sigmaTot=sigmaN + sum( piN ) + sum (sigmaRes )

    end subroutine evaluateXsections



    subroutine makeDecision
      use random, only: rn

      real :: summe, cut, cut2
      integer :: resID, totalCharge, pionCharge

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=sigma_particle%Charge+partNucl%Charge
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

      !########################################################################
      ! (2) Two -body final states
      !########################################################################

      ! sigma N production
      if (sigmaN.ge.cut) then
         partOut(1)%Id=sigmaMeson
         partOut(2)%Id=nucleon

         partOut(1)%Charge=0
         partOut(2)%Charge=totalCharge
         return
      end if
      cut=cut-sigmaN

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
      write(*,*) 'Error in makedecision of sigmaNuc', piN, cut
      stop

    end subroutine makeDecision

    !**************************************************************************
    !****s* sigmaNuc/makeOutput
    ! NAME
    ! subroutine makeOutput
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV].
    ! Filenames:
    ! * 'sigmaN_sigTotElast.dat'        : sigmaTot, sigmaElast
    ! * 'sigmaN_resProd.dat'            : Baryon resonance production
    ! * 'sigmaN_nonStrange_nuk.dat'     : non-strange meson with nucleon in final state
    !**************************************************************************
    subroutine makeOutPut
      logical, save :: initFlag=.true.
      real :: plab
      character(len=30), parameter :: outputFile(1:3) = (/ 'sigmaN_sigTotElast.dat   ', &
                                                           'sigmaN_resProd.dat       ', &
                                                           'sigmaN_nonStrange_nuk.dat' /)

      plab=SQRT(((srts**2-sigma_particle%mass**2-partNucl%mass**2)/2./partNucl%mass)**2-sigma_particle%mass**2)

      if (initFlag) then
         open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         open(file=outputFile(3),UNIT=103,Status='Replace',Action='Write')
         write(101,*) '# srts, plab, sigmaTot, sigmaElast '
         write(102,*) '# srts, plab, sigmaRes(2:40) '
         write(103,*) '# srts, plab, piN(-1:1), sigmaN '
         initFlag=.false.
      else
         open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
         open(file=outputFile(3),UNIT=103,Status='old',Position='Append',Action='Write')
      end if
      write(101,'(4F9.3)') srts, plab,sigmaTot, sigmaElast
      write(102,'(41F9.3)') srts, plab, sigmaRes(2:40)
      write(103,'(6F9.3)') srts, plab,piN, sigmaN
      close(101)
      close(102)
      close(103)

    end subroutine makeOutPut

  end subroutine sigmaNuc


end module sigmaNucleon
