!******************************************************************************
!****m* /etaNucleon
! NAME
! module etaNucleon
!
! PURPOSE
! Includes the cross sections for eta-nucleon scattering in the resonance regime
!
! Public routines:
! * etaNuc
!******************************************************************************
module etaNucleon
  implicit none
  private

  ! Debug-flags
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! To decide wether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.

  public :: etaNuc

contains

  !****************************************************************************
  !****s* etaNucleon/etaNuc
  ! NAME
  ! subroutine etaNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,useHiEnergy,HiEnergySchwelle,plotFlag)
  !
  ! PURPOSE
  ! Evaluates eta Nucleon -> anything cross sections and returns also a "preevent"
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
  ! * 2-particle : pi N, K Lambda, K Sigma
  !****************************************************************************


  subroutine etaNuc(srts,partIn,mediumAtColl,momLRF,partOut,sigmaTot,sigmaElast,&
       & useHiEnergy,HiEnergySchwelle,plotFlag)
    use idTable
    use particleDefinition
    use particleProperties, only: hadron
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
    real,dimension(-1:1) ::     piN            ! -> pi N, index denotes pion charge
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Local variables
    real :: fluxCorrector        ! Correction of the fluxfactor due to different velocities
                                 ! in the medium compared to the vacuum
    real :: s
    type(particle) :: eta_particle, partNucl
    logical :: antiParticleInput, failFlag

    ! partial cross sections for eta N -> R
    real, dimension(Delta:nbar) :: sigmaRes
    ! Field to store the resonance masses
    real, dimension(Delta:nbar) :: massRes      ! Resonance masses

    real :: lambdaKaon
    real, dimension  (0:1)  :: sigmaKaon             ! index = charge of final state kaon

    antiParticleINPUT=.false. ! .true. if antiparticle in the input

    ! Initialize output
    partOut(:)%ID=0                    ! ID of produced particles
    partOut(:)%charge=0                ! Charge of produced particles
    partOut(:)%antiParticle=.false.    ! Whether produced particles are particles or antiparticles
    partOut(:)%mass=0                  ! Mass of produced particles

    ! (1) Check  Input
    call searchInInput(partIn,eta,nucleon,eta_particle,partNucl,failFlag)
    if (failFlag) then
       write(*,*) 'Wrong input in EtaNuc', partIn%ID
    end if

    if (abs(eta_particle%charge).gt.1) write(*,*) 'wrong eta charge in etaNuc', eta_particle%charge

    if (eta_particle%antiParticle) then
       ! This case is not considered yet
       write(*,*) 'eta is antiparticle in "etaNuc"!!!',partIn%ID,partIn%antiparticle
       stop
    end if

    if (partNucl%antiParticle) then
       ! Invert all particles in antiparticles
       partNucl%Charge        =  -partNucl%Charge
       partNucl%antiparticle  = .false.
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

    s=srts**2

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
    if (Sum(partOut(:)%Charge).ne.partNucl%charge+eta_particle%charge) then
       write(*,*) 'No charge conservation in pionNuc!!! Critical error' ,eta_particle%Charge, &
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
      use parametrizationsBarMes, only: huang, huanglam
      use parBarMes_HighEnergy, only: paramBarMesHE
      use clebschGordan, only: clebschSquared
      use constants, only: mPi, mK
      use twoBodyTools, only: pCM

      real, dimension(1:3) ::  position
      logical :: perturbative
      real, dimension(0:3) :: momentum_vacuum
      real :: sigmaTotal_HE,sigmaElast_HE
      real :: p_piN, p_etaN, ratio
      real, dimension(1:4) :: sigmaHuang
      integer :: pionCharge, nucCharge

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
      ! eta N -> pi N
      !************************************************************************

      ! Full resonance contribution in the vacuum
      sigmaRes = barMes2resonance (eta,nucleon,eta_particle%charge,partNucl%charge,.true.,vacuum, &
                                   momentum_vacuum,massRes,eta_particle%Mass,partNucl%Mass,position,perturbative,srts)

      ! High energy matching
      piN=0.
      if (useHiEnergy) then
         if (srts.ge.2) then
            call paramBarMesHE(HiEnergySchwelle,eta,nucleon,eta_particle%charge,partNucl%charge,&
                 & mediumAtColl,sigmaTotal_HE,sigmaElast_HE)
            do pionCharge=-1,1
               nucCharge=eta_particle%charge+partNucl%charge-pionCharge
               if ((nucCharge.eq.0).or.(nucCharge.eq.1)) then
                  piN(pionCharge)=Max(sigmaTotal_HE-sum(sigmaRes),0.)*clebschSquared(1.,0.5,0.5,&
                       & real(pionCharge),real(nucCharge)-0.5)
               end if
            end do
         else
            piN=0.
         end if
      end if


      !************************************************************************
      ! eta N -> R
      !************************************************************************

      ! Full resonance contribution in the medium
      sigmaRes = barMes2resonance (eta,nucleon,eta_particle%charge,partNucl%charge,.true.,mediumAtColl, &
                                   momLRF,massRes,eta_particle%Mass,partNucl%Mass,position,perturbative,srts)

      !########################################################################
      ! evaluate elastic Xsection
      !########################################################################

      sigmaElast=barMes_R_barMes(eta,nucleon,eta,nucleon,&
           & eta_particle%Charge,partNucl%Charge,eta_particle%Charge,partNucl%Charge, &
           & .false.,.false.,mediumAtColl,momLRF,&
           & eta_particle%Mass,partNucl%Mass,position,perturbative,srts)

      ! Correction factor to the pion-nucleon cross sections by detailed balance
      ! (see J. Cugnon et al, PRC 40, 1822 (1989))
      ! c.m. momenta of pion-nucleon and eta-nucleon
      p_piN  = pCM(srts,mPi,partNucl%Mass)
      p_etaN = pCM(srts,eta_particle%Mass,partNucl%Mass)

      if (p_etaN.gt.1.e-06) then
         ratio=p_piN/p_etaN
      else
         ratio=1.
      end if

      !************************************************************************
      ! -> Lambda Kaon
      !************************************************************************
      lambdaKaon=0.
      if (srts > hadron(Lambda)%mass + mK) then
         ! huanglam gives : pi^{-} p -> Lambda Kaon^{0}
         lambdaKaon = 0.5 * huangLam(srts) * ratio  ! assume that sigma(eta p -> Lambda K^+) = sigma(pi^0 p -> Lambda K^+) * p_piN/p_etaN
         ! No resonance contribution
      end if

      !************************************************************************
      ! -> Sigma Kaon
      !************************************************************************
      ! sigmaKaon(0:1) : Index is charge of final state kaon
      ! sigmaHuang(1) = pi^{+}  p  ->   K^{+}  Sigma+      !
      ! sigmaHuang(2) = pi^{0}  p  ->   K^{+}  Sigma0       !
      ! sigmaHuang(3) = pi^{-}  p  ->   K^{0}  Sigma0      !
      ! sigmaHuang(4) = pi^{-}  p  ->   K^{+}   Sigma-     !
      sigmaKaon(:)=0.
      if (srts > hadron(SigmaResonance)%mass + mK) then
         sigmaHuang = huang(srts) * ratio     ! correction due to detailed balance
         if (partNucl%Charge.eq.1) then
            sigmaKaon(1)=sigmaHuang(2)     ! assume that sigma(eta p -> K^+ Sigma^0) = sigma(pi^0 p -> K^+ Sigma^0)
            sigmaKaon(0)=2.*sigmaKaon(1)  ! by isospin consideration
         else                              ! neutron Xsections by charge conjugation
            sigmaKaon(0)=sigmaHuang(2)
            sigmaKaon(1)=2.*sigmaKaon(0)
         end if
      end if

      !########################################################################
      ! Do the flux correction for each channel
      !########################################################################

      if (fluxCorrector_flag) then
         ! We do this for each channel since they might show up seperately in the output if makeoutput is called
         sigmaElast=sigmaElast*fluxcorrector
         piN=piN*fluxcorrector
         sigmaRes=sigmaRes*fluxcorrector
         lambdaKaon=lambdaKaon*fluxcorrector
         sigmaKaon=sigmaKaon*fluxcorrector
      end if

      !########################################################################
      ! Sum up everything for the total cross section
      !########################################################################
      ! Be careful since sigma elast is already included in the partial cross sections, therefore it is not
      ! included in the total cross section

      sigmaTot=sum( piN ) +sum (sigmaRes )  + lambdaKaon + sum(sigmaKaon)

    end subroutine evaluateXsections


    subroutine makeDecision
      use random, only: rn

      real :: summe, cut, cut2
      integer :: resID, totalCharge, charge, pionCharge

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=eta_particle%Charge+partNucl%Charge
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

      ! Kaon Lambda production
      if (lambdaKaon .ge.cut) then
         partOut(1)%Id=kaon
         partOut(2)%Id=Lambda
         partOut(1)%Charge=totalCharge
         partOut(2)%Charge=0
         return
      end if
      cut=cut-lambdaKaon

      ! Kaon Sigma production
      if (sum(sigmaKaon) .ge.cut) then
         partOut(1)%Id=kaon
         partOut(2)%Id=sigmaResonance
         cut2=rn()*sum(sigmaKaon)
         do charge=0,1
            if (sum(sigmaKaon(0:charge)).ge.cut2) then
               partOut(1)%Charge=charge
               partOut(2)%Charge=totalCharge-charge
               exit
            end if
         end do
         return
      end if

      ! Not event was generated:
      write(*,*) 'Error in makedecision of etaNuc', piN, cut
      stop

    end subroutine makeDecision



    !**************************************************************************
    !****s* etaNuc/makeOutput
    ! NAME
    ! subroutine makeOutput
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV].
    ! Filenames:
    ! * 'etaN_sigTotElast.dat'        : sigmaTot, sigmaElast
    ! * 'etaN_resProd.dat'            : Baryon resonance production
    ! * 'etaN_nonStrange_nuk.dat'     : non-strange meson with nucleon in final state
    ! * 'etaN_strangeProd.dat'        : Kaon and hyperon in final state
    !**************************************************************************
    subroutine makeOutPut
      logical, save :: initFlag=.true.
      real :: plab
      character(len=30), parameter :: outputFile(1:4) = (/ 'etaN_sigTotElast.dat   ', 'etaN_resProd.dat       ', &
                                                           'etaN_nonStrange_nuk.dat', 'etaN_strangeProd.dat   ' /)

      plab=SQRT(((s-eta_particle%mass**2-partNucl%mass**2)/2./partNucl%mass)**2-eta_particle%mass**2)

      if (initFlag) then
         open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         open(file=outputFile(3),UNIT=103,Status='Replace',Action='Write')
         open(file=outputFile(4),UNIT=104,Status='Replace',Action='Write')
         write(101,*) '# srts, plab, sigmaTot, sigmaElast '
         write(102,*) '# srts, plab, sigmaRes(2:40)'
         write(103,*) '# srts, plab, piN(-1:1)'
         write(104,*) '# srts, plab, lambdaKaon, sigmaKaon(0:1)'
         initFlag=.false.
      else
         open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
         open(file=outputFile(3),UNIT=103,Status='old',Position='Append',Action='Write')
         open(file=outputFile(4),UNIT=104,Status='old',Position='Append',Action='Write')
      end if
      write(101,'(4F9.3)') srts, plab,sigmaTot, sigmaElast
      write(102,'(41F9.3)') srts, plab,sigmaRes(2:40)
      write(103,'(5F9.3)') srts, plab,piN
      write(104,'(5(1x,e13.6))') srts, plab,lambdaKaon, sigmaKaon
      close(101)
      close(102)
      close(103)
      close(104)

    end subroutine makeOutPut
  end subroutine etaNuc


end module etaNucleon
