!******************************************************************************
!****m* /kaonLambda_resonance
! NAME
! module kaonLambda_resonance
! PURPOSE
! Includes the cross sections for kaon-delta scattering in the resonance regime
! Implemented are the following reactions:
!   * kaon Lambda -> X
! Public routines:
!   * kaonLambda
!******************************************************************************
module kaonLambda_resonance

  implicit none
  private

  ! Debug-flags
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! To decide wether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.

  public :: kaonLambda

contains

  !****************************************************************************
  !****m* kaonLambda_resonance/kaonLambda
  ! NAME
  ! subroutine kaonLambda(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,plotFlag)
  !
  ! PURPOSE
  ! Evaluates kaon Lambda -> anything cross sections and returns also a "preevent"
  !
  ! RESULT
  ! * real, intent(out)                                        :: sigmaTot         ! total Xsection
  ! * real, intent(out)                                        :: sigmaElast       ! elastic Xsection
  !
  ! This routine does a Monte-Carlo-decision according to the partial cross sections to decide on a final state with
  ! maximal 2 final state particles. These are returned in the vector teilchenOut. The kinematics of these teilchen is
  ! only fixed in the case of a single produced resonance. Otherwise the kinematics still need to be established. The
  ! result is:
  ! * type(preEvent),dimension(1:3), intent(out)               :: teilchenOut     ! outgoing particles
  !
  ! NOTES
  ! Possible final states are :
  ! * 1-particle : baryon Resonances
  ! * 2-particle : pi N, pi Delta
  !****************************************************************************

  subroutine kaonLambda(srts,partIn,mediumAtColl,momLRF,partOut,sigmaTot,sigmaElast,plotFlag)

    use idTable, only: nucleon, delta, pion, kaon, lambda, nbar
    use particleDefinition
    use mediumDefinition
    use preEventDefinition, only: preEvent
    use twoBodyTools, only: velocity_correction, convertToAntiParticles, pcm,searchInInput
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
    ! CrossSections named according to final states
    !real ::  lambdaKaon
    real, dimension  (-1:1)  :: piN             ! index = charge of final state pion
    real, dimension (-1:1)  :: piDelta          ! index = charge of final state kaon



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Local variables
    real :: fluxCorrector      ! Correction of the fluxfactor due to different velocities
                               ! in the medium compared to the vacuum
    type(particle) :: kaon_particle, lambda_particle
    logical :: antiParticleInput, failFlag

    ! partial cross sections for kaon Lamba -> R
    real, dimension(Delta:nbar) :: sigmaRes
    ! Field to store the resonance masses
    real , dimension(Delta:nbar) :: massRes      ! Resonance masses

    antiParticleINPUT=.false. ! .true. if antiparticle in the input

    ! Initialize output
    partOut(:)%ID=0                    ! ID of produced particles
    partOut(:)%charge=0                ! Charge of produced particles
    partOut(:)%antiParticle=.false.    ! Whether produced particles are particles or antiparticles
    partOut(:)%mass=0                  ! Mass of produced particles

    ! (1) Check  Input
    call searchInInput(partIn,kaon,lambda,kaon_particle,lambda_particle,failFlag)
    if (failFlag) then
       write(*,*) 'Wrong input in KaonNuc', partIn%ID
       stop
    end if

    if (lambda_particle%antiParticle.and.kaon_particle%antiParticle) then
       ! Both are antiparticles: s=0 scattering channel
       !
       ! Invert all particles in antiparticles
       lambda_particle%Charge        =  -lambda_particle%Charge
       lambda_particle%antiparticle  = .false.
       kaon_particle%Charge          =  -kaon_particle%Charge
       kaon_particle%antiparticle  = .false.
       antiParticleInput=.true.
    else if ((.not.(lambda_particle%antiParticle)).and.(.not.(kaon_particle%antiParticle))) then
       ! Both are no antiparticles : s=0 scattering channel
       antiParticleInput=.false.
    else
       ! Not yet implemented: S=-2,-1,1,2 scattering channels
       sigmaTot=0.
       sigmaElast=0.
       return
    end if

    if (debugFlag) then
       write(*,*) '##################################################'
       write(*,*) 'Input'
       write(*,*) 'q=',lambda_particle%charge,kaon_particle%charge
       write(*,*) 'id=', lambda_particle%id,kaon_particle%id
       write(*,*) '##################################################'
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
    if (Sum(partOut(:)%Charge).ne.kaon_particle%Charge+lambda_particle%Charge) then
       write(*,*) 'No charge conservation in kaonLambda!!! Critical error' ,kaon_particle%Charge,lambda_particle%Charge,&
            &  partOut(:)%Charge,partOut(:)%ID
       stop
    end if

    if (debugFlag) write(*,*) 'q=',partOut%charge
    if (debugFlag) write(*,*) 'id=',partOut%ID

    ! (6) Invert particles in antiParticles if input included antiparticles
    if (antiParticleInput) then
       if (debugFlagAnti) write(*,*) partOut
       call convertToAntiParticles(partOut)
       if (debugFlagAnti) write(*,*) partOut
    end if

    if (debugFlag) write(*,*) 'q=',partOut%charge
    if (debugFlag) write(*,*) 'id=',partOut%ID

  contains

    subroutine evaluateXsections
      !****s* kaonLambda/evaluateXsections
      ! NAME
      ! subroutine evaluateXsections
      ! PURPOSE
      ! Evaluates all Xsections for kaon Lambda scattering
      !************************************************************************
      use particleDefinition
      use mediumDefinition, only: vacuum
      use particleProperties, only: hadron
      use parametrizationsBarMes, only: huanglam, huanglamd
      use resonanceCrossSections, only: barMes_R_barMes, barMes2resonance
      use clebschGordan, only: clebschSquared
      use output, only: writeparticle
      use constants, only: mPi, mN

      real :: sigmaHuangLam
      real :: sigmaHuangLamd
      !logical :: background,propagated
      real,dimension(0:3)                 :: momentum_vacuum        ! Total Momentum in vacuum
      real, dimension(1:3) :: position
      logical :: perturbative

      integer:: pionCharge, deltaCharge


      real :: pInitial, pFinal, detailedBalanceFactor


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


      !************************************************************************
      ! -> pion N
      !************************************************************************

      piN=0.
      if (srts > mPi + mN) then
         pFinal   = pCM(srts,mPi,mN)
         pInitial = pCM(srts,kaon_particle%mass,lambda_particle%mass)

         if (pinitial.lt.1E-12) then
            write(*,*) 'WARNING: pInitial is zero in kaonlambda', pinitial
            write(*,*) 'kaon meson:'
            call writeparticle(6,0,0,kaon_Particle)
            write(*,*) 'Lambda:'
            call writeparticle(6,0,0,lambda_Particle)
            detailedBalanceFactor= 0.
         else
            detailedBalanceFactor= (pFinal/pInitial)**2
         end if

         ! huangLam gives : pi^{-} p -> Lambda Kaon^{0}
         sigmaHuangLam = huangLam(srts) * detailedBalanceFactor

         ! Subtract resonance contribution
         sigmaHuangLam=Max(0.,sigmaHuangLam-barMes_R_barMes(kaon,lambda,pion,nucleon, 0,0,-1,1,.false.,.true.,&
                    & Vacuum,momentum_vacuum,kaon_particle%Mass,lambda_particle%Mass, &
                    & position,perturbative,srts))

         do pionCharge=-1,1
            if ( ((kaon_particle%charge-pionCharge).eq.0).or.((kaon_particle%charge-pionCharge).eq.1)  ) then
               select case (pionCharge)
               case (-1)
                  piN(pionCharge)=sigmaHuangLam
               case (0)
                  piN(pionCharge)=0.5*sigmaHuangLam ! Clebsch-Gordan-Factor
               case (1)
                  piN(pionCharge)=sigmaHuangLam
               end select
            else
               piN(pionCharge)=0.
            end if
         end do
      end if
      if (debugFlag) write(*,*) 'piN=',piN
      !************************************************************************
      ! -> pion Delta
      !************************************************************************


      ! Simplification: neglect delta width
      piDelta=0.
      if (srts > mPi + hadron(delta)%mass) then
         pFinal   = pCM(srts,mPi,hadron(delta)%mass)
         pInitial = pCM(srts,kaon_particle%mass,lambda_particle%mass)

         if (pinitial.lt.1E-12) then
            write(*,*) 'WARNING: pInitial is zero in kaonlambda', pinitial
            write(*,*) 'kaon meson:'
            call writeparticle(6,0,0,kaon_Particle)
            write(*,*) 'Lambda:'
            call writeparticle(6,0,0,lambda_Particle)
            detailedBalanceFactor= 0.
         else
            detailedBalanceFactor= (pFinal/pInitial)**2  *2.
            ! factor 2  because of spins
         end if

         ! hunaglamd gives : delta^{++} pi- -> Lambda Kaon^{+}
         sigmaHuangLamd = huangLamd(srts) * detailedBalanceFactor

         do pionCharge=-1,1
            deltaCharge=lambda_particle%charge+kaon_particle%charge-pionCharge
            if ((deltaCharge.le.2).and.(deltaCharge.ge.-1)) then
               piDelta(pionCharge)=2.*sigmaHuangLamd*clebschSquared(1.5,1.0,0.5,real(deltaCharge)-0.5,real(pionCharge))
               ! divide out "delta^{++} pi- -> Lambda Kaon^{+}" clebsch and multiply be real clebsch
            end if
         end do
      end if

      if (debugFlag) write(*,*) 'piDelta=',piDelta
      !************************************************************************
      ! -> R
      !************************************************************************

      ! Full resonance contribution in the medium
      sigmaRes = barMes2resonance (kaon,lambda,kaon_particle%charge,lambda_particle%charge,.true.,mediumAtColl, &
                                   momLRF,massRes,kaon_particle%Mass,lambda_particle%Mass,position,perturbative,srts)

      if (debugFlag) write(*,*) 'sigmares=', sum(sigmares)
      !########################################################################
      ! evaluate elastic Xsection
      !########################################################################

      sigmaElast=barMes_R_barMes(kaon,lambda,kaon,lambda,&
           & kaon_particle%Charge,lambda_particle%Charge,kaon_particle%Charge,lambda_particle%Charge, &
           & .false.,.false.,mediumAtColl,momLRF,kaon_particle%Mass,lambda_particle%Mass, &
           & position,perturbative,srts)

      !########################################################################
      ! Do the flux correction for each channel
      !########################################################################

      ! We do this for each channel since they might show up seperately in the output if makeoutput is called
      if (fluxCorrector_flag) then
         ! Correction of the fluxfactor due to different velocities in the medium compared to the vacuum
         sigmaElast=sigmaElast*fluxcorrector
         piN=piN *fluxcorrector
         piDelta=piDelta*fluxcorrector
         sigmaRes=sigmaRes *fluxcorrector
      end if


      !########################################################################
      ! Sum up everything for the total cross section
      !########################################################################
      ! Be careful since sigma elast is already included in the partial cross sections, therefore it is not
      ! included in the total cross section

      sigmaTot=Sum(piDelta) +Sum ( piN ) + sum (sigmaRes )


    end subroutine evaluateXsections


    !**************************************************************************
    !****s* kaonLambda/makeDecision
    ! NAME
    ! subroutine MakeDecision
    ! PURPOSE
    ! Decides on the final state which is returned via teilchenOut by Monte-Carlo.
    !  * Assigns charges and ID's.
    !  * Only for resonance-production also the mass is assigned, since the mass of the resonance needed to be calculated earlier.
    ! The Monte-Carlo routine is adding up channels until the sum is exceeding x*sigma(total). x has a flat distribution in [0,1].
    ! The last added channel is then the one which is chosen for the event. After choosing a channel, the subroutine is returning to
    ! the calling routine.
    !**************************************************************************
    subroutine MakeDecision

      use random, only: rn

      real :: summe, cut, cut2
      integer :: totalCharge, resID, pionCharge

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=kaon_particle%Charge+lambda_particle%Charge
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

      if ( Sum(piN) .ge.cut) then
         partOut(1)%Id=pion
         partOut(2)%Id=nucleon
         cut2=rn()*sum(piN)
         do pioncharge=-1,1
            if (sum(piN(-1:pionCharge)).gt.cut2) then
               partOut(1)%Charge=pioncharge
               partOut(2)%Charge=totalCharge-pioncharge
               return
            end if
         end do
      end if
      cut=cut-Sum(piN)

      if ( Sum(piDelta) .ge.cut) then
         partOut(1)%Id=pion
         partOut(2)%Id=delta
         cut2=rn()*sum(piDelta)
         do pioncharge=-1,1
            if (sum(piDelta(-1:pionCharge)).gt.cut2) then
               partOut(1)%Charge=pioncharge
               partOut(2)%Charge=totalCharge-pioncharge
               return
            end if
         end do
      end if
      cut=cut-Sum(piDelta)


      write(*,*) 'Error in makeDecision : No decision made', &
           & cut,    sum(piDelta) , partOut(:)%ID,partOut(:)%Charge,sigmaTot

      Stop

    end subroutine MakeDecision


    !**************************************************************************
    !****s* kaonLambda/makeOutput
    ! NAME
    ! subroutine makeOutput
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV].
    ! Filenames:
    ! * 'kaonLambda_sigTotElast.dat'        : sigmaTot, sigmaElast
    ! * 'kaonLambda_nonstrangeProd.dat'     : non strangeness production
    ! * 'kaonLambda_resProd.dat'            : Baryon resonance production
    !**************************************************************************
    subroutine makeOutput
      logical, save :: initFlag=.true.
      real :: plab
      character(len=30), parameter :: outputFile(1:3) = (/ 'kaonLambda_sigTotElast.dat   ', &
                                                           'kaonLambda_nonstrangeProd.dat', &
                                                           'kaonLambda_resProd.dat       ' /)

      plab=SQRT(((srts**2-kaon_particle%Mass**2-lambda_particle%Mass**2)/2./lambda_particle%Mass)**2-kaon_particle%Mass**2)

      if (initFlag) then
         open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         open(file=outputFile(3),UNIT=103,Status='Replace',Action='Write')
         write(101,*) '# srts, plab, sigmaTot, sigmaElast '
         write(102,*) '# srts, plab, piN(-1:1), piDelta(-1:1)'
         write(103,*) '# srts, plab, sigmaRes(2:40)'
         initFlag=.false.
      else
         open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
         open(file=outputFile(3),UNIT=103,Status='old',Position='Append',Action='Write')
      end if
      write(101,'(4F9.3)') srts, plab,sigmaTot, sigmaElast
      write(102,'(9F9.3)') srts, plab,piN, piDelta
      write(103,'(41e12.3)') srts, plab,sigmaRes(2:40)
      close(101)
      close(102)
      close(103)
    end subroutine makeOutput

  end subroutine kaonLambda


end module kaonLambda_resonance
