!******************************************************************************
!****m* /rhoNucleon
! NAME
! module rhoNucleon
! PURPOSE
! Includes the cross sections for rho-nucleon scattering in the resonance regime
! Public routines:
! * rhoNuc
!******************************************************************************
module rhoNucleon
  implicit none
  private

  ! Debug-flags
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! To decide wether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.

  public :: rhoNuc

contains

  !****************************************************************************
  !****s* rhoNucleon/rhoNuc
  ! NAME
  !  subroutine rhoNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,plotFlag)
  ! PURPOSE
  ! Evaluates rho Nucleon -> anything cross sections and returns a "preevent".
  ! INPUTS
  ! * real, intent(in)                              :: srts                  ! sqrt(s) in the process
  ! * type(particle),dimension(1:2), intent(in)     :: teilchenIn            ! colliding particles
  ! * type(medium), intent(in)                      :: mediumATcollision     ! Medium informations at the position of the collision
  ! * real, intent(in) ,dimension(0:3)              :: momentumLRF           ! Total Momentum in LRF
  ! High energy matching:
  ! * logical,intent(in)                            :: useHiEnergy
  ! * .true. if High-Energy cross sections are given by paramBarMesHE
  ! * real,intent(in)                               :: HiEnergySchwelle
  ! * threshold sqrt(s) for paramBarMesHE, i.e. at which energy the cross sections of paramBarMesHE are used
  ! Debugging:
  ! * logical, intent(in),optional                  :: plotFlag              ! Switch on plotting of the  Xsections
  ! RESULT
  ! * real, intent(out)                                        :: sigmaTot         ! total Xsection
  ! * real, intent(out)                                        :: sigmaElast       ! elastic Xsection
  ! * type(preEvent),dimension(1:3), intent(out)               :: teilchenOut     ! colliding particles
  ! NOTES
  ! This routine does a Monte-Carlo-decision according to the partial cross sections to decide on a final state with
  ! maximal 3 final state particles. These are returned in the vector teilchenOut. The kinematics of these teilchen is
  ! only fixed in the case of a single produced resonance. Otherwise the kinematics still need to be established.
  ! Possible final states are :
  ! * 1-particle : baryon Resonances
  ! * 2-particle : pi N, K Lambda, K Sigma
  !****************************************************************************
  subroutine rhoNuc (srts, partIn, mediumAtColl, momLRF, partOut, sigmaTot, sigmaElast, &
                     useHiEnergy, HiEnergySchwelle, plotFlag)

    use idTable
    use particleDefinition
    use mediumDefinition
    use preEventDefinition, only: preEvent
    use twoBodyTools, only: velocity_correction, convertToAntiParticles, searchInInput
    use RMF, only: getRMF_flag
    use CALLSTACK, only: TRACEBACK
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Input
    real, intent(in)                           :: srts                  ! sqrt(s) in the process
    type(particle), dimension(1:2), intent(in) :: partIn            ! colliding particles
    type(medium), intent(in)                   :: mediumAtColl     ! Medium informations at the position of the collision
    real, intent(in) ,dimension(0:3)           :: momLRF           ! Total Momentum in LRF
    logical, intent(in)                        :: useHiEnergy           ! .true. if High-Energy cross sections are given by paramBarMesHE
    real, intent(in)                           :: HiEnergySchwelle      ! threshold sqrt(s) for paramBarMesHE
    logical, intent(in), optional              :: plotFlag              ! Switch on plotting of the  Xsections
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Output
    type(preEvent),dimension(1:3), intent(out) :: partOut      ! colliding particles
    real, intent(out)                          :: sigmaTot         ! total Xsection
    real, intent(out)                          :: sigmaElast       ! elastic Xsection

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Cross sections
    real, dimension(Delta:nbar) :: sigmaRes      ! rho N -> R cross section
    real, dimension(-1:1)   :: piN           ! -> pi N, index denotes pion charge
    real                    :: lambdaKaon
    real, dimension(0:1)    :: sigmaKaon     ! index = charge of final state kaon
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Field to store the resonance masses
    real, dimension(Delta:nbar) :: massRes       !  Resonance masses
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Local variables
    real :: fluxCorrector        ! Correction of the fluxfactor due to different velocities in the medium compared to the vacuum
    real :: s
    type(particle) :: rho_particle, partNucl
    logical :: antiParticleInput, failFlag

    antiParticleINPUT=.false. ! .true. if antiparticle in the input

    ! Initialize output
    partOut(:)%ID=0                    ! ID of produced particles
    partOut(:)%charge=0                ! Charge of produced particles
    partOut(:)%antiParticle=.false.    ! Whether produced particles are particles or antiparticles
    partOut(:)%mass=0                  ! Mass of produced particles

    ! (1) Check  Input
    call searchInInput(partIn,rho,nucleon,rho_particle,partNucl,failFlag)

    if (failFlag) write(*,*) 'Wrong input in RhoNuc', partIn%ID
    if (abs(rho_particle%charge)>1) write(*,*) 'wrong rho charge in rhoNuc', rho_particle%charge

    if (rho_particle%antiParticle) then
       ! This case is not considered yet
       write(*,*) 'rho is antiparticle in "rhoNuc"!!!',partIn%ID,partIn%antiparticle
       stop
    end if

    if (partNucl%antiParticle) then
       ! Invert all particles in antiparticles
       partNucl%Charge        =  -partNucl%Charge
       partNucl%antiparticle  = .false.
       rho_particle%Charge          =  -rho_particle%Charge
       antiParticleInput=.true.
    else
       antiParticleInput=.false.
    end if

    ! Correction of the fluxfactor due to different velocities in the medium compared to the vacuum
    if (.not. getRMF_flag()) then
      fluxCorrector=velocity_correction(partIn)
    else
      fluxCorrector=1.
    end if

    s=srts**2

    ! (2) Evaluate the cross sections
    call evaluateXsections

    ! Cutoff to kick those case out, that the cross section is zero
    if (sigmaTot<1E-12) then
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
    if (Sum(partOut(:)%Charge).ne.partNucl%charge+rho_particle%charge) then
       write(*,*) 'No charge conservation in pionNuc!!! Critical error' ,rho_particle%Charge, &
                  partNucl%Charge, partOut(:)%Charge,partOut(:)%ID
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
      use mediumDefinition, only: vacuum
      use parBarMes_HighEnergy, only: paramBarMesHE
      use parametrizationsBarMes, only: huangLam, huang
      use constants, only: mPi, mK
      use particleProperties, only: hadron

      real, dimension(1:3) :: position
      real, dimension(0:3) :: momentum_vacuum
      logical :: perturbative
      real :: sigmaTotal_HE, sigmaElast_HE, sigmaHuangLam, sigma_R, sigmaHuang(1:4)
      real :: p_piN, p_rhoN, ratio

      position=0.5*(partIn(1)%position+partIn(2)%position)
      perturbative = (partIn(1)%perturbative .or. partIn(2)%perturbative)

      momentum_vacuum(1:3)=partIn(1)%momentum(1:3)+partIn(2)%momentum(1:3)
      momentum_vacuum(0)=FreeEnergy(partIn(1))+FreeEnergy(partIn(2))

      !########################################################################
      ! Evaluate partial cross sections
      !########################################################################


      !************************************************************************
      ! rho N -> pi N
      !************************************************************************

      ! Full resonance contribution in the vacuum
      sigmaRes = barMes2resonance (rho,nucleon,rho_particle%charge,partNucl%charge,.true.,vacuum, &
                                   momentum_vacuum,massRes,rho_particle%Mass,partNucl%Mass,position,perturbative,srts)

      ! High energy matching
      piN=0.
      if (useHiEnergy) then
         if (srts>=1.8) then
            call paramBarMesHE (HiEnergySchwelle,rho,nucleon,rho_particle%charge,partNucl%charge, &
                                mediumAtColl,sigmaTotal_HE,sigmaElast_HE)
            piN(rho_particle%charge)=Max(sigmaTotal_HE-sum(sigmaRes),0.)
         end if
      else
         ! Set minimal total cross section to 30 mB
         piN(rho_particle%charge)=30.-sum(sigmaRes)
      end if


      !************************************************************************
      ! rho N -> R
      !************************************************************************

      sigmaRes=0.
      ! Full resonance contribution in the medium
      sigmaRes = barMes2resonance (rho,nucleon,rho_particle%charge,partNucl%charge,.true.,mediumAtColl, &
                                   momLRF,massRes,rho_particle%Mass,partNucl%Mass,position,perturbative,srts)

      !########################################################################
      ! evaluate elastic Xsection
      !########################################################################

      sigmaElast = barMes_R_barMes (rho,nucleon,rho,nucleon, &
                                    rho_particle%Charge,partNucl%Charge,rho_particle%Charge,partNucl%Charge, &
                                    .false.,.false.,mediumAtColl,momLRF, &
                                    rho_particle%Mass,partNucl%Mass,position,perturbative,srts)

      ! Correction factor to the pion-nucleon cross sections by detailed balance
      ! (see J. Cugnon et al, PRC 40, 1822 (1989))
      ! c.m. momenta of pion-nucleon and rho-nucleon
      p_piN=(s+mPi**2-partNucl%Mass**2)**2 /(4.*s) - mPi**2
      p_piN=sqrt(max(0.,p_piN))
      p_rhoN=(s+rho_particle%Mass**2-partNucl%Mass**2)**2 /(4.*s) - rho_particle%Mass**2
      p_rhoN=sqrt(max(0.,p_rhoN))
      if (p_rhoN>1.e-06) then
         ratio=p_piN/p_rhoN
      else
         ratio=1.
      end if

      !************************************************************************
      ! -> Lambda Kaon
      !************************************************************************
      lambdaKaon=0.
      if (srts>(hadron(Lambda)%mass+mK)) then
         ! huanglam gives : pi^{-} p -> Lambda Kaon^{0}
         sigmaHuangLam = huangLam(srts) * ratio ! correction due to detailed balance
         sigma_R = barMes_R_barMes(rho,nucleon,kaon,Lambda, &
                                   -1,1,0,0,.false.,.true.,Vacuum,momentum_vacuum, &
                                   rho_particle%Mass,partNucl%Mass,position,perturbative,srts)
         select case (rho_particle%Charge)
         case (-1)
            if (partNucl%Charge==1) lambdaKaon = sigmaHuangLam - sigma_R
         case (0)
            lambdaKaon = 0.5 * (sigmaHuangLam - sigma_R)
         case (1)
            if (partNucl%Charge==0) lambdaKaon = sigmaHuangLam - sigma_R
         end select

         if (lambdaKaon<0.) then
            lambdaKaon=0.
            write(*,*) 'Problem in rhoNuc : rho N -> Lambda Kaon resonance', &
                       '  Contribution greater than elementary crossSection'
         end if
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
      if (srts>(hadron(SigmaResonance)%mass+mK)) then
         sigmaHuang = huang(srts) * ratio  ! correction due to detailed balance
         if (partNucl%Charge==1) then
            select case (rho_particle%Charge)
            case (-1)
               sigmaKaon(1)=sigmaHuang(4)
               sigmaKaon(0)=sigmaHuang(3)
            case (0)
               sigmaKaon(1)=sigmaHuang(2)
               sigmaKaon(0)=sigmaHuang(3)   ! by isospin consideration
            case (1)
               sigmaKaon(1)=sigmaHuang(1)
               sigmaKaon(0)=0.
            end select
         else ! neutron Xsections by charge conjugation
            select case (rho_particle%Charge)
            case (-1)
               sigmaKaon(1)=0.
               sigmaKaon(0)=sigmaHuang(1)
            case (0)
               sigmaKaon(1)=sigmaHuang(3)
               sigmaKaon(0)=sigmaHuang(2)
            case (1)
               sigmaKaon(1)=sigmaHuang(3)
               sigmaKaon(0)=sigmaHuang(4)
            end select
         end if
      end if

      !########################################################################
      ! Do the flux correction for each channel
      !########################################################################

      if (fluxCorrector_flag) then
         ! We do this for each channel since they might show up seperately in the output if makeoutput is called
         sigmaElast=sigmaElast*fluxcorrector
         piN=piN*fluxcorrector
         sigmaRes=sigmaRes *fluxcorrector
         lambdaKaon=lambdaKaon *fluxcorrector
         sigmaKaon=sigmaKaon*fluxcorrector
      end if

      !########################################################################
      ! Sum up everything for the total cross section
      !########################################################################
      ! Be careful since sigma elast is already included in the partial cross sections, therefore it is not
      ! included in the total cross section

      sigmaTot = sum(piN) + sum(sigmaRes) + lambdaKaon + sum(sigmaKaon)

    end subroutine evaluateXsections



    subroutine makeDecision
      use random, only: rn

      real :: summe, cut, cut2
      integer :: resID, totalCharge, pionCharge, charge

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=rho_particle%Charge+partNucl%Charge
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
         if (piN(pionCharge)>=cut) then
            partOut(1:2)%Id = (/pion,nucleon/)
            partOut(1:2)%Charge = (/pionCharge,totalCharge-pionCharge/)
            return
         end if
         cut=cut-piN(pionCharge)
      end do

      ! Lambda Kaon production
      if (lambdaKaon>=cut) then
         partOut(1:2)%Id = (/kaon,Lambda/)
         partOut(1:2)%Charge = (/totalCharge,0/)
         return
      end if
      cut=cut-lambdaKaon

      ! Sigma Kaon production
      if (sum(sigmaKaon)>=cut) then
         partOut(1:2)%Id = (/kaon,sigmaResonance/)
         cut2=rn()*sum(sigmaKaon)
         do charge=0,1
            if (sum(sigmaKaon(0:charge))>=cut2) then
               partOut(1:2)%Charge = (/charge,totalCharge-charge/)
               exit
            end if
         end do
         return
      end if

      ! No event was generated:
      write(*,*) 'Error in makedecision of rhoNuc', cut
      write(*,*) sigmaTot
      write(*,*) sigmaRes
      write(*,*)
      write(*,*) piN
      write(*,*)
      write(*,*) lambdaKaon
      write(*,*) sigmaKaon
      call TRACEBACK("")

    end subroutine makeDecision


    !**************************************************************************
    !****s* rhoNuc/makeOutput
    ! NAME
    ! subroutine makeOutput
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV].
    ! Filenames:
    ! * 'rhoN_sigTotElast.dat'        : sigmaTot, sigmaElast
    ! * 'rhoN_resProd.dat'            : Baryon resonance production
    ! * 'rhoN_nonStrange_nuk.dat'     : non-strange meson with nucleon in final state
    ! * 'rhoN_strangeProd.dat'        : Kaon and hyperon in final state
    !**************************************************************************
    subroutine makeOutPut
      logical, save :: initFlag=.true.
      real :: plab
      character(len=23), dimension(1:4), parameter :: outputFile = (/ 'rhoN_sigTotElast.dat   ', 'rhoN_resProd.dat       ', &
                                                                      'rhoN_nonStrange_nuk.dat', 'rhoN_strangeProd.dat   ' /)

      plab=SQRT(((s-rho_particle%mass**2-partNucl%mass**2)/2./partNucl%mass)**2-rho_particle%mass**2)

      if (initFlag) then
         open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         open(file=outputFile(3),UNIT=103,Status='Replace',Action='Write')
         open(file=outputFile(4),UNIT=104,Status='Replace',Action='Write')
         write(101,*) '# srts, plab, sigmaTot, sigmaElast '
         write(102,*) '# srts, plab, sum(sigmaRes), sigmaRes(2:40) '
         write(103,*) '# srts, plab, piN(-1:1)'
         write(104,*) '# srts, plab, lambdaKaon, sigmaKaon(0:1)'
         initFlag=.false.
      else
         open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
         open(file=outputFile(3),UNIT=103,Status='old',Position='Append',Action='Write')
         open(file=outputFile(4),UNIT=104,Status='old',Position='Append',Action='Write')
      end if
      write(101, '(4F9.3)') srts, plab, sigmaTot, sigmaElast
      write(102,'(42F9.3)') srts, plab, sum(sigmaRes), sigmaRes(2:40)
      write(103, '(5F9.3)') srts, plab, piN
      write(104, '(5F9.3)') srts, plab, lambdaKaon, sigmaKaon
      close(101)
      close(102)
      close(103)
      close(104)

    end subroutine makeOutPut

  end subroutine rhoNuc


end module rhoNucleon
