!******************************************************************************
!****m* /omegaNucleon
! NAME
! module omegaNucleon
! PURPOSE
! Includes the cross sections for omega-nucleon scattering in the resonance
! regime.
!******************************************************************************
module omegaNucleon
  implicit none
  private

  ! Debug-flags:
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! use the flux correction for the incoming particle velocities:
  logical, parameter :: fluxCorrector_flag=.true.

  public :: omegaNuc

contains

  !****************************************************************************
  !****s* omegaNucleon/omegaNuc
  ! NAME
  ! subroutine omegaNuc(srts,teilchenIN,mediumATcollision,momentumLRF,
  ! teilchenOUT,sigmaTot,sigmaElast,useHiEnergy,HiEnergySchwelle,plotFlag,
  ! sigmaArr)
  !
  ! PURPOSE
  ! Evaluates omega Nucleon -> anything cross sections and returns also a
  ! "preevent"
  !
  ! INPUTS
  ! * real :: srts --- sqrt(s) in the process
  ! * type(particle),dimension(1:2) :: teilchenIn --- colliding particles
  ! * type(medium) :: mediumATcollision --- Medium informations at the position
  !   of the collision
  ! * real, dimension(0:3) :: momentumLRF --- Total Momentum in LRF
  !
  ! High energy matching:
  ! * logical :: useHiEnergy ---  .true. if High-Energy cross sections are
  !   given by paramBarMesHE
  ! * real :: HiEnergySchwelle --- threshold sqrt(s) for paramBarMesHE,
  !   i.e. at which energy the cross sections of paramBarMesHE are used
  !
  ! Debugging:
  ! * logical,optional :: plotFlag --- Switch on plotting of the  Xsections
  !
  ! RESULT
  ! * real :: sigmaTot --- total Xsection
  ! * real :: sigmaElast --- elastic Xsection
  !
  ! This routine does a Monte-Carlo-decision according to the partial cross
  ! sections to decide on a final state with maximal 3 final state particles.
  ! These are returned in the vector teilchenOut. The kinematics of these
  ! teilchen is only fixed in the case of a single produced resonance.
  ! Otherwise the kinematics still need to be established. The
  ! result is:
  ! * type(preEvent),dimension(1:3) :: teilchenOut --- particles
  ! * real, dimension(6), optional :: sigmaArr -- partial cross sections
  !
  ! The cross sections are based upon a parametrization by Golubeva.
  ! See routine golub_omega in parametrizationBarMes.
  ! NOTES
  ! Possible final states are :
  ! * 1-particle : baryon resonances
  ! * 2-particle : pi N, omega N, K Lambda, K Sigma
  ! * 3-particle : pi pi N
  !****************************************************************************
  subroutine omegaNuc(srts, partIn, mediumAtColl, momLRF, partOut, sigmaTot, &
       sigmaElast, useHiEnergy, HiEnergySchwelle, plotFlag, K_Factor, sigmaArr)
    use idTable
    use particleDefinition
    use particleProperties, only: hadron
    use mediumDefinition
    use preEventDefinition, only: preEvent
    use twoBodyTools, only: velocity_correction, convertToAntiParticles, &
         pcm,searchInInput
    use RMF, only: getRMF_flag
    use constants, only: mN, mPi, mK
    use callstack, only: traceBack


    real, intent(in)                            :: srts
    type(particle), dimension(1:2), intent(in)  :: partIn
    type(medium), intent(in)                    :: mediumAtColl
    real, intent(in) ,dimension(0:3)            :: momLRF
    logical, intent(in)                         :: useHiEnergy
    real, intent(in)                            :: HiEnergySchwelle
    logical, intent(in), optional               :: plotFlag
    real, intent(in)                            :: K_Factor

    type(preEvent), dimension(1:3), intent(out) :: partOut
    real, intent(out)                           :: sigmaTot, sigmaElast
    real, dimension(6), intent(out), optional   :: sigmaArr

    ! Cross sections
    real,dimension(-1:1) :: piN, piN_R ! -> pi N, index denotes pion charge
    real :: omegaN, pipiN, lambdaKaon
    real, dimension(0:1)  :: sigmaKaon ! index = charge of final state kaon

    real :: fluxCorrector
    real :: s
    type(particle) :: partOmega, partNucl
    logical :: antiParticleInput, failFlag
    real, dimension(Delta:nbar) :: sigmaRes, massRes
    ! partial cross sections for omega N -> R, and resonance masses

    type(preEvent) :: preEv0

    antiParticleINPUT=.false. ! .true. if antiparticle in the input

    ! Initialize output
    partOut = preEv0

    ! (1) Check  Input
    call searchInInput(partIn,omegaMeson,nucleon,partOmega,partNucl,failFlag)
    if (failFlag)  call traceBack('wrong input')

    if (partOmega%antiParticle)  call traceBack('meson is anti')

    if (partNucl%antiParticle) then
       ! Invert all particles in antiparticles
       partNucl%Charge        =  -partNucl%Charge
       partNucl%antiparticle  = .false.
       partOmega%Charge          =  -partOmega%Charge
       antiParticleInput              = .true.
    else
       antiParticleInput=.false.
    end if


    ! Correction of the fluxfactor due to different velocities in the medium compared to the vacuum
    if (.not.getRMF_flag()) then
      fluxCorrector=velocity_correction(partIn)
    else
      fluxCorrector=1.
    end if

    s=srts**2

    ! (2) Evaluate the cross sections
    call evaluateXsections

    ! (2a) only fill array, no event
    if (present(sigmaArr)) then
       sigmaArr = (/ &
            omegaN, sum(piN), pipiN, &
            sum(sigmaRes), lambdaKaon, sum(sigmaKaon) /)
       return
    end if

    ! Cutoff to kick the case out, that the cross section is zero
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
    if (Sum(partOut(:)%Charge).ne.partNucl%charge+partOmega%charge) then
       write(*,*) 'No charge conservation in pionNuc!!! Critical error', partOmega%Charge, partNucl%Charge, &
                                                                         partOut(:)%Charge, partOut(:)%ID
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
      use parametrizationsBarMes, only: golub_omega, omegaN_lykasov, huang, huanglam
      use parBarMes_HighEnergy, only: paramBarMesHE
      use output, only: writeparticle

      real, dimension(1:3) ::  position
      logical :: perturbative
      real, dimension(1:2) :: sigmaGolub
      real, dimension(1:4) :: sigmaHuang
      real, dimension(-1:1) :: piN_Golub
      real :: p_piN, p_omegaN, ratio
      real :: pFinal,pInitial,detailedBalanceFactor,elast_R
      integer :: pionCharge, nucCharge
      real, dimension(0:3) :: momentum_vacuum
      real :: sigmaTotal_HE,sigmaElast_HE ! High energy matchin
      logical :: lDummy

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
      ! omega N -> pi N
      !************************************************************************
      ! piN = [cross section by Golubeva (pi^- p-> omega N) by
      ! detailed balance] - [vacuum cross section via resonances]

      if (srts > mPi + mN) then
        sigmaGolub = golub_omega (srts)            ! Results by Golubeva

        pFinal = pCM(srts, mPi, mN)
        pInitial = pCM(srts, partOmega%mass, partNucl%mass)

        if (pinitial<1E-12) then
          write(*,*) 'WARNING: pInitial is zero in omegaNuc', pinitial
          write(*,*) 'omega meson:'
          call writeparticle(6,0,0,partOmega)
          write(*,*) 'nucleon:'
          call writeparticle(6,0,0,partNucl)
          detailedBalanceFactor= 0.
        else
          detailedBalanceFactor= 1./3.*(pFinal/pInitial)**2
          ! given by detailed balance: factor 1/3 due to (2j+1)-Terms in cross
          ! section and different spins in initial and final state
        end if

        sigmaGolub(1) = 3./2. * sigmaGolub(1)
        ! 3./2. since Golub returns cross section pi^- p -> omega N,
        ! with 3./.2 we divide by the isospin-clebsch

        if (partNucl%charge==0) then
           piN_golub = (/2./3.,1./3.,0./) * sigmaGolub(1) * detailedBalanceFactor
        else if (partNucl%charge==1) then
           piN_golub = (/0.,1./3.,2./3./) * sigmaGolub(1) * detailedBalanceFactor
        else
           write(*,*) 'Error in omegaNuc', partNucl
        end if

        do pionCharge=-1,1
           nucCharge=partOmega%charge+partNucl%charge-pionCharge
           if ((nucCharge==0).or.(nucCharge==1)) then
              piN_R(pionCharge) = barMes_R_barMes(omegaMeson, nucleon, &
                   pion, nucleon, &
                   partOmega%Charge, partNucl%Charge, pionCharge, nucCharge, &
                   .false., .false., Vacuum, momentum_vacuum, partOmega%Mass, partNucl%Mass, &
                   position, perturbative, srts)
              piN(pionCharge) = max(0., piN_golub(pionCharge) - piN_R(pionCharge))
          else
              piN(pionCharge) = 0.
          end if
        end do
      else
        piN = 0.
      end if

      !************************************************************************
      ! omega N -> omega N
      !************************************************************************

      elast_R = barMes_R_barMes(omegaMeson,nucleon,omegaMeson,nucleon,&
           partOmega%Charge,partNucl%Charge,partOmega%Charge,partNucl%Charge, &
           .false.,.false.,mediumAtColl,momLRF, &
           partOmega%Mass,partNucl%Mass,position,perturbative,srts)

      ! subtract resonance contribution from (direct) omega N channel
      if (srts>mN+hadron(omegaMeson)%minmass) then
        omegaN=Max(0., omegaN_lykasov(srts,partOmega%mass,1)  - elast_R)
      else
        omegaN=0.
      end if

      !************************************************************************
      ! omega N -> R
      !************************************************************************

      ! Full resonance contribution in the medium
      sigmaRes = barMes2resonance (omegaMeson,nucleon,partOmega%charge,partNucl%charge,.true.,mediumAtColl, &
                                   momLRF,massRes,partOmega%Mass,partNucl%Mass,position,perturbative,srts)

      !************************************************************************
      ! -> pi pi N
      !************************************************************************

      ! Evaluate Omega N -> N pi pi by taking the total inelastic cross section
      ! by Lykasov (times a K-factor) and subtracting then all included
      ! inelastic cross sections.
      ! Inelastic resonance contribution = full contribution - elastic
      ! contribution

      if (srts > 2*mPi + mN) then
         pipiN=max(0.,omegaN_lykasov(srts,partOmega%mass,2)*K_Factor &
              - (Sum(sigmaRes)-elast_R)-Sum(piN))
         ! Matching to High energy region
         if (useHiEnergy) then
            call paramBarMesHE(HiEnergySchwelle,omegaMeson,nucleon,&
                 partOmega%charge,partNucl%charge, &
                 mediumAtColl,sigmaTotal_HE,sigmaElast_HE)
            sigmaTotal_HE = (sigmaTotal_HE-sigmaElast_HE)*K_factor &
                 + sigmaElast_HE
            pipiN=max(pipiN,sigmaTotal_HE-Sum(sigmaRes)-Sum(piN)-omegaN)
         end if
      else
         piPiN=0.
      end if

      !########################################################################
      ! evaluate elastic Xsection
      !########################################################################

      sigmaElast = omegaN + elast_R

      ! Correction factor to the pion-nucleon cross sections by detailed balance
      ! (see J. Cugnon et al, PRC 40, 1822 (1989))
      ! c.m. momenta of pion-nucleon and omega-nucleon
      p_piN    = pCM(s,mPi,partNucl%Mass, lDummy)
      p_omegaN = pCM(s,partOmega%Mass,partNucl%Mass, lDummy)
      if (p_omegaN.gt.1.e-06) then
         ratio=p_piN/p_omegaN
      else
         ratio=1.
      end if

      !************************************************************************
      ! -> Lambda Kaon
      !************************************************************************
      lambdaKaon=0.
      if (srts.gt.(hadron(Lambda)%mass+mK)) then
         ! huanglam gives : pi^{-} p -> Lambda Kaon^{0}
         lambdaKaon = 0.5 * huangLam(srts) * ratio
         ! assume that sigma(omega p -> Lambda K^+)
         ! = sigma(pi^0 p -> Lambda K^+) * p_piN/p_omegaN
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
      if (srts.gt.(hadron(SigmaResonance)%mass+mK)) then
         sigmaHuang = huang(srts) * ratio       ! correction due to detailed balance
         if (partNucl%Charge.eq.1) then
            ! assume that sigma(omega p -> K^+ Sigma^0) = sigma(pi^0 p -> K^+ Sigma^0)
            ! by isospin consideration
            sigmaKaon(1)=sigmaHuang(2)
            sigmaKaon(0)=2.*sigmaKaon(1)
         else
            ! neutron Xsections by charge conjugation
            sigmaKaon(0)=sigmaHuang(2)
            sigmaKaon(1)=2.*sigmaKaon(0)
         end if
      end if

      !########################################################################
      ! Do the flux correction for each channel
      !########################################################################

      if (fluxCorrector_flag) then
         ! We do this for each channel since they might show up seperately in
         ! the output if makeoutput is called
         sigmaElast=sigmaElast*fluxcorrector
         omegaN=omegaN*fluxcorrector
         piN=piN*fluxcorrector
         pipiN=pipiN*fluxcorrector
         sigmaRes=sigmaRes *fluxcorrector
         lambdaKaon=lambdaKaon*fluxcorrector
         sigmaKaon=sigmaKaon*fluxcorrector
      end if

      !########################################################################
      ! Sum up everything for the total cross section
      !########################################################################
      ! Be careful since sigma elast is already included in the partial cross
      ! sections (omegaN and sigmaRes), therefore it is not included in the
      ! total cross section.

      sigmaTot = omegaN + sum(piN) + pipiN + sum(sigmaRes) + lambdaKaon &
           + sum(sigmaKaon)

    end subroutine evaluateXsections


    subroutine makeDecision
      use random, only: rn, ranCharge

      real :: summe, cut, cut2
      integer :: resID, totalCharge, pionCharge, charge
      integer,dimension (1:3)  :: izmin,izmax,izout     ! needed for ranCharge
      logical :: ranChargeFlag

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=partOmega%Charge+partNucl%Charge
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

      ! omega N production
      if (omegaN.ge.cut) then
         partOut(1:2)%ID = (/ omegaMeson, nucleon /)
         partOut(1:2)%Charge = (/ 0, totalCharge /)
         return
      end if
      cut=cut-omegaN

      ! piN production
      do pionCharge=-1,1
         if (piN(pionCharge).ge.cut) then
            partOut(1:2)%ID = (/ pion, nucleon /)
            partOut(1:2)%Charge = (/ pionCharge, totalCharge-pionCharge /)
            return
         end if
         cut=cut-piN(pionCharge)
      end do

      ! Kaon Lambda production
      if (lambdaKaon .ge.cut) then
         partOut(1:2)%ID = (/ kaon, Lambda /)
         partOut(1:2)%Charge = (/ totalCharge, 0 /)
         return
      end if
      cut=cut-lambdaKaon

      ! Kaon Sigma production
      if (sum(sigmaKaon) .ge.cut) then
         partOut(1:2)%ID = (/ kaon, SigmaResonance /)
         cut2=rn()*sum(sigmaKaon)
         do charge=0,1
            if (sum(sigmaKaon(0:charge)).ge.cut2) then
               partOut(1:2)%Charge = (/ charge, totalCharge-charge /)
               exit
            end if
         end do
         return
      end if
      cut=cut-sum(sigmaKaon)


      !########################################################################
      ! (3) Three body final state
      !########################################################################
      ! pi pi N production

      if (pipiN.ge.cut) then
         partOut(1:3)%ID = (/pion,pion,nucleon/)
         izmin=(/-1,-1,0/)
         izmax=(/1,1,1/)
         call rancharge(izmin,izmax,totalCharge,izout,ranChargeFlag)
         if (.not.ranChargeFlag) write(*,*) 'Error in rancharge :',izmin,izmax,totalCharge,izout,ranChargeFlag
         partOut(1:3)%Charge = izout(1:3)
         return
      end if

      ! No event was generated:
      write(*,*) 'Error in makedecision of omegaNuc', &
           srts, sigmaTot, omegaN, sum(piN), pipiN, sum(sigmaRes), cut
      stop

    end subroutine makeDecision


    !**************************************************************************
    !****s* omegaNucleon/makeOutput
    ! NAME
    ! subroutine makeOutput
    !
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV]
    ! .
    ! Filenames:
    ! * 'omegaN_sigTotElast.dat'    : sigmaTot, sigmaElast
    ! * 'omegaN_resProd.dat'        : Baryon production (Resonances with ID's 1:40)
    ! * 'omegaN_nonStrange_nuk.dat' : non-strange meson with nucleon in final state
    ! * 'omegaN_strangeProd.dat'    : Kaon and hyperon in final state
    !**************************************************************************

    subroutine makeOutPut

      logical, save :: initFlag=.true.
      real :: plab
      character(30), parameter :: outputFile(1:4) = (/ &
           'omegaN_sigTotElast.dat   ', 'omegaN_resProd.dat       ', &
           'omegaN_nonStrange_nuk.dat', 'omegaN_strangeProd.dat   ' /)

      plab=SQRT(((s-partOmega%mass**2-partNucl%mass**2)/2./partNucl%mass)**2-partOmega%mass**2)

      if (initFlag) then
         open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         open(file=outputFile(3),UNIT=103,Status='Replace',Action='Write')
         open(file=outputFile(4),UNIT=104,Status='Replace',Action='Write')
         write(101,*) '# srts, plab, sigmaTot, sigmaElast '
         write(102,*) '# srts, plab, sum(sigmaRes), sigmaRes(2:40) '
         write(103,*) '# srts, plab, piN(-1:1), piN_R(-1:1), omegaN, pipiN'
         write(104,*) '# srts, plab, lambdaKaon, sigmaKaon(0:1)'
         initFlag=.false.
      else
         open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
         open(file=outputFile(3),UNIT=103,Status='old',Position='Append',Action='Write')
         open(file=outputFile(4),UNIT=104,Status='old',Position='Append',Action='Write')
      end if
      write(101, '(4F10.4)') srts, plab, sigmaTot, sigmaElast
      write(102,'(42F10.4)') srts, plab, sum(sigmaRes), sigmaRes(2:40)
      write(103,'(10F10.4)') srts, plab, piN, piN_R, omegaN, pipiN
      write(104, '(5F10.4)') srts, plab, lambdaKaon, sigmaKaon
      close(101)
      close(102)
      close(103)
      close(104)

    end subroutine makeOutPut
  end subroutine omegaNuc


end module omegaNucleon
