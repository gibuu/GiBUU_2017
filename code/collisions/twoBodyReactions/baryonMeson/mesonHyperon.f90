!******************************************************************************
!****m* /mesonHyperon
! NAME
! module mesonHyperon
! PURPOSE
! Includes the cross sections for nonstrange meson - Hyperon inelastic scattering.
! Public routines:
! * mesonY
!******************************************************************************
module mesonHyperon

  implicit none
  private

  ! Debug-flags
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! To decide wether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.

  public :: mesonY

contains

  !****************************************************************************
  !****s* mesonHyperon/mesonY
  ! NAME
  ! subroutine mesonY(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,plotFlag)
  !
  ! PURPOSE
  ! Evaluates Y pi -> Y^*, Lambda eta -> Lambda(1670),
  ! Y^* pi -> Y^*, Y pi ->  Kbar N and Y^* pi ->  Kbar N
  !                 as well as
  ! Ybar pi -> Y^*bar,  LambdaBar eta -> antiLambdaBar(1670),
  ! Y^*bar pi -> Y^*bar, Ybar pi ->  K Nbar
  ! and Y^*bar pi ->  K Nbar cross sections and returns also a "preevent"
  !
  ! INPUTS
  ! * real, intent(in)                              :: srts                  ! sqrt(s) in the process
  ! * type(particle),dimension(1:2), intent(in)     :: teilchenIn            ! colliding particles
  ! * type(medium), intent(in)                      :: mediumATcollision     ! Medium informations at the position of the collision
  ! * real, intent(in) ,dimension(0:3)              :: momentumLRF           ! Total Momentum in LRF
  !
  ! Debugging:
  ! * logical, intent(in),optional                  :: plotFlag              ! Switch on plotting of the  Xsections
  !
  ! OUTPUT
  ! * real, intent(out)                                        :: sigmaTot       ! total Xsection
  ! * real, intent(out)                                        :: sigmaElast     ! elastic Xsection
  !
  ! This routine does a Monte-Carlo-decision according to the partial cross sections to decide on a final state with
  ! maximal 2 final state particles. These are returned in the vector teilchenOut. The kinematics of these teilchen is
  ! only fixed in the case of a single produced resonance. Otherwise the kinematics still need to be established. The result is:
  ! * type(preEvent),dimension(1:3), intent(out)               :: teilchenOut    !   outgoing particles
  !
  ! The cross sections are based on parameterization by M. Effenberger.
  ! NOTES
  ! Possible final states are :
  ! * 1-particle : Y^*, Ybar^*
  ! * 2-particle : Kbar N, K Nbar
  !****************************************************************************
  subroutine mesonY(srts,partIn,mediumAtColl,momLRF,partOut,sigmaTot,sigmaElast,plotFlag)
  use idTable
  use particleDefinition
  use particleProperties, only: hadron
  use mediumDefinition
  use preEventDefinition, only: preEvent
  use twoBodyTools, only: velocity_correction,p_lab,convertToAntiParticles
  use RMF, only: getRMF_flag

  !============================================================================
  ! INPUT:
  real, intent(in)                              :: srts                  ! sqrt(s) in the process
  type(particle),dimension(1:2), intent(in)     :: partIn            ! colliding particles
  type(medium), intent(in)                      :: mediumAtColl     ! Medium informations at the position of the collision
  real, intent(in) ,dimension(0:3)              :: momLRF           ! Total Momentum in LRF
  logical, intent(in),optional                  :: plotFlag              ! Switch on plotting of the  Xsections
  !============================================================================
  ! OUTPUT:
  real, intent(out)                             :: sigmaTot, sigmaElast  ! total & elastic cross section
  type(preEvent),dimension(1:3), intent(out)    :: partOut           !   outgoing particles
  !============================================================================

  real :: sigmaRes_tot ! total cross section meson Y(Y^*) -> Resonances
  real :: sigbg        ! background cross sections pi Y(Y^*) -> Kbar N  (summed over final isospins)

  ! Local variables
  real :: fluxCorrector        ! Correction of the fluxfactor due to different velocities in the medium compared to the vacuum
  type(particle) :: meson_particle, hyperon_particle
  integer :: isrts,qmeson,qhyperon
  logical :: antiParticleInput
  integer :: totalCharge   ! total charge of colliding particles

  ! Parameters for discretisation of sqrts (same as in Effe's code):
  real, parameter :: srtmess0=1.254    !*srtmess0=m_lambda+m_pion
  real, parameter :: dmssrt=0.01       ! sqrts-bin size
  integer, parameter :: nssmmax=300    ! number of sqrts-bins

  ! Fields to store discretised background cross sections
  ! (same as in Effe's code) and respective resonance contributions
  ! to the given channels:

  !    Lambda pi^0 -> K^- p, K^0 n (sum over final isospins):
  real, dimension(0:nssmmax), save :: slaka,slaka_res

  !    Sigma pi -> Kbar N:
  !    first index = charge of incoming state = 3*q_Sigma+q_meson
  !    second index = 1 sum over final isospins
  !                 = 0 outgoing Kbar^0
  !                 = -1 outgoing K^-
  real, dimension(-4:4,-1:1,0:nssmmax), save ::  ssigka,ssigka_res

  logical, save :: flagini=.true.
  ! partial cross sections (and resonance masses) for meson Y -> R
  real, dimension(Delta:nbar) :: sigmaRes, massRes

  if (flagini) then
    ! Discretization of cross sections:
    call init
    flagini=.false.
  end if

  ! Initialize output
  partOut(:)%ID=0                    ! ID of produced particles
  partOut(:)%charge=0                ! Charge of produced particles
  partOut(:)%antiParticle=.false.    ! Whether produced particles are particles or antiparticles
  partOut(:)%mass=0.                 ! Mass of produced particles

  ! (1) Check  Input
  if ((partIn(1)%Id==pion .or. partIn(1)%Id==eta) &
      .and.partIn(2)%Id>=Lambda .and. partIn(2)%Id<=omegaResonance) then
    meson_particle   = partIn(1)
    hyperon_particle = partIn(2)
  else if ((partIn(2)%Id==pion .or. partIn(2)%Id==eta) &
           .and.partIn(1)%Id>=Lambda .and. partIn(1)%Id<=omegaResonance) then
    meson_particle   = partIn(2)
    hyperon_particle = partIn(1)
  else
    write(*,*) 'Wrong input in mesonHyperon', partIn%ID
    stop
  end if

  if (hyperon_particle%antiParticle) then
    ! Invert all particles in antiparticles:
    hyperon_particle%Charge=-hyperon_particle%Charge
    hyperon_particle%antiparticle=.false.
    meson_particle%Charge=-meson_particle%Charge
    antiParticleInput=.true.
  else
    antiParticleInput=.false.
  end if

  qmeson=meson_particle%Charge
  qhyperon=hyperon_particle%Charge

  ! Correction of the fluxfactor due to different velocities in the medium compared to the vacuum
  if (.not. getRMF_flag()) then
    fluxCorrector=velocity_correction(partIn)
  else
    fluxCorrector=1.
  end if

  ! (2) Evaluate the cross sections
  call evaluateXsections

  ! Cutoff to kick the case out, that the cross section is zero
  if (sigmaTot < 1E-12) then
    sigmatot=0.
    sigmaElast=0.
    return
  end if

  ! (3) Plot them if wished
  if (present(PlotFlag) .or. debugFlag) then
    if (plotFlag .or. debugFlag)  call makeOutput
  end if

  ! (4) Define final state
  call MakeDecision

  ! (5) Check Output
  if (Sum(partOut(:)%Charge) /= hyperon_particle%charge+meson_particle%charge) then
       write(*,*) 'No charge conservation in mesonHyperon!!! Critical error',&
                  meson_particle%Charge,hyperon_particle%Charge,&
                  partOut(:)%Charge,partOut(:)%ID
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
    !****s* mesonHyperon/evaluateXsections
    ! NAME
    ! subroutine evaluateXsections
    ! PURPOSE
    ! Evaluates Y(Y^*) pi -> Y^* and Y(Y^*) pi ->  Kbar N
    ! cross sections.
    ! NOTES
    ! There is a resonance contribution to Y(Y^*) pi scattering.
    !**************************************************************************
    subroutine evaluateXsections

    use resonanceCrossSections, only: barMes2resonance
    use constants, only: mN, mK
    use twoBodyTools, only: pCM

    real, dimension(1:3) :: position
    logical :: perturbative

    real, parameter :: m2yst_km=22., m2yst_k0=15. ! (mb GeV^2) matrix element squared for the Kbar N -> Y^* pi (see Effe's PhD, p. 239)

    real :: pinitial,pfinal,isoFactor

    position=0.5*(partIn(1)%position+partIn(2)%position)
    perturbative = (partIn(1)%perturbative .or. partIn(2)%perturbative)

    ! Resonance contribution:
    sigmaRes = barMes2resonance(meson_particle%Id,hyperon_particle%Id,qmeson,qhyperon,.true.,mediumAtColl,&
                                momLRF,massRes,meson_particle%mass,hyperon_particle%mass,position,perturbative,srts)
    sigmaRes_tot=sum(sigmaRes(:))

    sigbg=0.

    ! Non-resonant contributions:
    if (meson_particle%Id==pion .and. srts > mN + mK) then

        !  pi Y(^*) -> Kbar N:
        isrts=min(nint((srts-srtmess0)/dmssrt),nssmmax)

        select case (hyperon_particle%ID)
        case (Lambda)
          !* Lambda pi -> N Kbar:
          sigbg=slaka(isrts)
        case (SigmaResonance)
          !* Sigma pi -> N Kbar:
          sigbg=ssigka(3*qhyperon+qmeson,1,isrts)
        case (lambda_1600:sigma_1915)
          !* Y* pi -> N Kbar:
          pinitial = pCM(srts,hyperon_particle%mass,meson_particle%mass)  ! Initial c.m. momentum
          pfinal   = pCM(srts,mN,mK)                                      ! Final c.m. momentum
          if (abs(qhyperon+qmeson)==2) then
            isoFactor=0.
          else if (hadron(hyperon_particle%Id)%isoSpinTimes2==0) then
            isoFactor=1.
          else if (qhyperon+qmeson==0) then
            if (qmeson==0) then
              isoFactor=1./3.
            else
              isoFactor=5./6.
            end if
          else
            isoFactor=1./2.
          end if
          if (qhyperon+qmeson==0) then
             sigbg=isoFactor*2./(2.*hadron(hyperon_particle%Id)%spin+1.) &
                   *pfinal/pinitial/srts**2 * m2yst_km
          else
             sigbg=isoFactor*2./(2.*hadron(hyperon_particle%Id)%spin+1.) &
                   *pfinal/pinitial/srts**2 * m2yst_k0
          end if
        end select

    end if

    ! Elastic and total cross sections:
    ! NOTE: elastic Xsection is only due to the resonance contribution.
    ! Therefore, it is computed only for plotting purposes
    ! (see subroutine evaluateResContr, called from makeoutput below).
    ! Otherwise it is not needed and putted to zero to save CPU time.
    sigmaElast=0.
    sigmaTot=sigbg+sigmaRes_tot

    ! Flux correction for each channel:
    if (fluxCorrector_flag) then
      sigmaRes=sigmaRes*fluxcorrector
      sigmaRes_tot=sigmaRes_tot*fluxcorrector
      sigbg=sigbg*fluxcorrector
      sigmaElast=sigmaElast*fluxcorrector
      sigmaTot=sigmaTot*fluxcorrector
    end if

    end subroutine evaluateXsections


    !**************************************************************************
    !****s* mesonHyperon/init
    ! NAME
    ! subroutine init
    ! PURPOSE
    ! Tabulates the backgrond cross sections for  Y pi ->  Kbar N.
    !**************************************************************************
    subroutine init

    use baryonWidth, only: partialWidthBaryon, fullWidthBaryon
    use parametrizationsBarMes, only: kaonbg
    use constants, only: pi, mN, mPi, mK
    use clebschGordan, only: ClebschSquared
    use twoBodyTools, only: pCM_sqr

    real :: srts,s,momCM2_Lambda_pi,momCM2_Sigma_pi,gamma_tot
    real :: iso_res,gamma_in_Lambda_pi,gamma_in_Sigma_pi,gamma_out
    real :: SpinFactor,sigma_res_Lambda_pi,sigma_res_Sigma_pi
    real :: isoFactor
    integer :: i,k,resID

  !    Lambda pi^0 -> K^- p, K^0 n (sum over final isospins):
  !  real, dimension(0:nssmmax), save :: slaka,slaka_res

  !    Sigma pi -> Kbar N:
  !    first index = charge of incoming state = 3*q_Sigma+q_meson
  !    second index = 1 sum over final isospins
  !                 = 0 outgoing Kbar^0
  !                 = -1 outgoing K^-
  !  real, dimension(-4:4,-1:1,0:nssmmax), save ::  ssigka,ssigka_res

    slaka=0.
    slaka_res=0.
    ssigka=0.
    ssigka_res=0.

    Loop_over_srts : do i=0,nssmmax

      srts=float(i)*dmssrt+srtmess0
      if (srts <= mN + mK) cycle
      s=srts**2

      ! Lambda pi0 -> K- p, K0 n:
      slaka(i)=2.*kaonbg(srts,3,1,0)

      ! Sigma- pi0 -> K- n:
      ssigka(-3,-1,i)=kaonbg(srts,4,1,1)

      ! Sigma- pi0 -> Kbar0 N:
      ssigka(-3,0,i)=0.

      ! Sigma- pi+ -> K- p:
      ssigka(-2,-1,i)=kaonbg(srts,5,1,0)

      ! Sigma- pi+ -> Kbar0 n:
      ssigka(-2,0,i)=kaonbg(srts,4,1,0)

      ! Sigma0 pi- -> K- n:
      ssigka(-1,-1,i)=kaonbg(srts,4,1,1)

      ! Sigma0 pi- -> Kbar0 N:
      ssigka(-1,0,i)=0.

      ! Sigma0 pi0 -> K- p:
      ssigka(0,-1,i)=kaonbg(srts,6,1,0)

      ! Sigma0 pi0 -> K0 n:
      ssigka(0,0,i)=ssigka(0,-1,i)

      ! Initial c.m. momenta squared:
      momCM2_Lambda_pi = pCM_sqr(s,mPi**2,hadron(Lambda)%mass**2)
      momCM2_Sigma_pi  = pCM_sqr(s,mPi**2,hadron(SigmaResonance)%mass**2)

      ! Resonance contribution to the various channels:
      ! NOTE: if a resonance is not explicitly propagated, these contributions
      ! are added to the background fields

      ! Loop over S=-1 resonances:
      do resID=Sigma_1385,Sigma_1915

        if (.not.hadron(resId)%usedForXsections) cycle   ! Exclude resonances

        ! Total width in vacuum:
        gamma_tot=FullWidthBaryon(resID,srts)

        ! Isospin of the resonance:
        iso_res=float(hadron(resID)%isoSpinTimes2)/2.

        ! In-widths:

        if (abs(iso_res-1.).lt.0.0001) then
          gamma_in_Lambda_pi=partialwidthBaryon(resID,srts,.true.,pion,Lambda)
        else
          gamma_in_Lambda_pi=0.
        end if

        gamma_in_Sigma_pi=partialwidthBaryon(resID,srts,.true.,pion,SigmaResonance)

        ! Out-width:
        gamma_out=partialwidthBaryon(resID,srts,.false.,kaonBar,nucleon)

        SpinFactor=(2.*hadron(resId)%spin+1.)/2.

        sigma_res_Lambda_pi=SpinFactor*4.*pi/momCM2_Lambda_pi*s &
                &*gamma_in_Lambda_pi &
                &/((s-hadron(resID)%mass**2)**2+gamma_tot**2*s)*0.389

        sigma_res_Sigma_pi=SpinFactor*4.*pi/momCM2_Sigma_pi*s &
                &*gamma_in_Sigma_pi &
                &/((s-hadron(resID)%mass**2)**2+gamma_tot**2*s)*0.389


        ! Lambda pi -> Kbar N:
        if (hadron(resId)%propagated) then
          slaka_res(i)=slaka_res(i)+sigma_res_Lambda_pi*gamma_out
        else
          slaka(i)=slaka(i)+sigma_res_Lambda_pi*gamma_out
        end if

        ! Sigma- pi0 -> K- n:
        isoFactor=ClebschSquared(1.,1.,iso_res,-1.,0.)&
                &*ClebschSquared(0.5,0.5,iso_res,-0.5,-0.5)
        if (hadron(resId)%propagated) then
          ssigka_res(-3,-1,i)=ssigka_res(-3,-1,i) &
                            &+isoFactor*sigma_res_Sigma_pi*gamma_out
        else
          ssigka(-3,-1,i)=ssigka(-3,-1,i) &
                            &+isoFactor*sigma_res_Sigma_pi*gamma_out
        end if

        ! Sigma- pi0 -> Kbar0 N:
        ssigka_res(-3,0,i)=0.

        ! Sigma- pi+ -> K- p:
        isoFactor=ClebschSquared(1.,1.,iso_res,-1.,1.)&
                &*ClebschSquared(0.5,0.5,iso_res,-0.5,0.5)
        if (hadron(resId)%propagated) then
          ssigka_res(-2,-1,i)=ssigka_res(-2,-1,i) &
                            &+isoFactor*sigma_res_Sigma_pi*gamma_out
        else
          ssigka(-2,-1,i)=ssigka(-2,-1,i) &
                            &+isoFactor*sigma_res_Sigma_pi*gamma_out
        end if

        ! Sigma- pi+ -> Kbar0 n:
        ! isoFactor is the same as previous
        if (hadron(resId)%propagated) then
          ssigka_res(-2,0,i)=ssigka_res(-2,0,i) &
                            &+isoFactor*sigma_res_Sigma_pi*gamma_out
        else
          ssigka(-2,0,i)=ssigka(-2,0,i) &
                            &+isoFactor*sigma_res_Sigma_pi*gamma_out
        end if

        ! Sigma0 pi- -> K- n:
        isoFactor=ClebschSquared(1.,1.,iso_res,0.,-1.)&
                &*ClebschSquared(0.5,0.5,iso_res,-0.5,-0.5)
        if (hadron(resId)%propagated) then
          ssigka_res(-1,-1,i)=ssigka_res(-1,-1,i) &
                            &+isoFactor*sigma_res_Sigma_pi*gamma_out
        else
          ssigka(-1,-1,i)=ssigka(-1,-1,i) &
                            &+isoFactor*sigma_res_Sigma_pi*gamma_out
        end if

        ! Sigma0 pi- -> Kbar0 N:
        ssigka_res(-1,0,i)=0.

        ! Sigma0 pi0 -> K- p:
        isoFactor=ClebschSquared(1.,1.,iso_res,0.,0.)&
                &*ClebschSquared(0.5,0.5,iso_res,-0.5,0.5)
        if (hadron(resId)%propagated) then
          ssigka_res(0,-1,i)=ssigka_res(0,-1,i) &
                            &+isoFactor*sigma_res_Sigma_pi*gamma_out
        else
          ssigka(0,-1,i)=ssigka(0,-1,i) &
                            &+isoFactor*sigma_res_Sigma_pi*gamma_out
        end if

        ! Sigma0 pi0 -> K0 n:
        ! isoFactor is the same as previous
        if (hadron(resId)%propagated) then
          ssigka_res(0,0,i)=ssigka_res(0,0,i) &
                            &+isoFactor*sigma_res_Sigma_pi*gamma_out
        else
          ssigka(0,0,i)=ssigka(0,0,i) &
                            &+isoFactor*sigma_res_Sigma_pi*gamma_out
        end if

      end do
      ! End of loop over S=-1 resonances

      !*use isospin symmetrie for other channels         :
      do k=1,3
        ssigka(k,-1,i)=ssigka(-k,0,i)
        ssigka(k,0,i)=ssigka(-k,-1,i)
        ssigka_res(k,-1,i)=ssigka_res(-k,0,i)
        ssigka_res(k,0,i)=ssigka_res(-k,-1,i)
      end do

      !*isospin sums:
      do k=-3,3
        ssigka(k,1,i)=ssigka(k,-1,i)+ssigka(k,0,i)
        ssigka_res(k,1,i)=ssigka_res(k,-1,i)+ssigka_res(k,0,i)
      end do


    end do Loop_over_srts

    end subroutine init


    !**************************************************************************
    !****s* mesonHyperon/makeDecision
    ! NAME
    ! subroutine makeDecision
    ! PURPOSE
    ! Chooses randomly one of possible outgoing channels in meson-hyperon
    ! collision. Outgoing channels are: Y^* and Kbar N.
    ! Also the charges of outgoing particles are selected.
    !**************************************************************************
    subroutine makeDecision

      use random, only: rn

      real :: cut,cut2,sum,x
      integer :: resID,qkaon

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=qmeson+qhyperon
      if (abs(totalCharge)>1) then
        write(*,*) ' In mesonHyperon/makeDecision: totalCharge=',totalCharge
        stop
      end if

      if (sigmaRes_tot >= cut) then

        ! Resonance production:
        sum=0.
        do resID=Delta,nbar
          sum=sum+sigmaRes(resID)
          if (sum >= cut) exit
        end do
        partOut(1)%Id=resID
        partOut(1)%Charge=totalCharge
        partOut(1)%Mass=massRes(resID)

      else

        ! Kbar N production (background)
        partOut(1:2)%Id = (/ kaonBar, nucleon /)
        if (totalCharge==1) then
          partOut(1:2)%Charge = (/0,1/)
        else if (totalCharge==-1) then
          partOut(1:2)%Charge = (/-1,0/)
        else
          if (hyperon_particle%Id==SigmaResonance) then
            ! pi Sigma -> Kbar N:
            sum=0.
            cut2=rn()*ssigka(3*qhyperon+qmeson,1,isrts)
            do qkaon=-1,0
              sum=sum+ssigka(3*qhyperon+qmeson,qkaon,isrts)
              if (sum >= cut2) exit
            end do
            partOut(1:2)%Charge = (/ qkaon, -qkaon /)
          else
            x=rn()
            if (x<=0.5) then
              partOut(1:2)%Charge = (/-1,1/)
            else
              partOut(1:2)%Charge = 0
            end if
          end if
        end if

      end if

    end subroutine makeDecision


    !**************************************************************************
    !****s* mesonHyperon/makeOutput
    ! NAME
    ! subroutine makeOutput
    !
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV]
    ! .
    ! Filenames:
    ! * 'mesonY_sigTotElast.dat'        : sigmaTot, sigmaElast
    ! * 'mesonY_KbarN.dat'          : Kbar Nucleon outgoing channel (summed
    ! *                              over final isospins)
    ! * 'pionLambda_KbarN.dat'   : pi Lambda -> Kbar N isospin averaged
    ! * 'pionSigma_KbarN.dat'   : pi Sigma -> Kbar N isospin averaged
    !**************************************************************************
    subroutine makeOutPut

      use constants, only: mPi

      logical, save :: initFlag=.true.

      ! The output files
      character(30), dimension(1:4), parameter :: outputFile = (/ 'mesonY_sigTotElast.dat', 'mesonY_KbarN.dat      ', &
                                                                  'pionLambda_KbarN.dat  ', 'pionSigma_KbarN.dat   ' /)
      real :: plab,plab_Lambda,plab_Sigma, sigmaElast_res,sigmaKbarN_res
      integer :: isrts

      plab        = p_lab(srts,meson_particle%mass,hyperon_particle%mass)
      plab_Lambda = p_lab(srts,mPi,hadron(Lambda)%mass)
      plab_Sigma  = p_lab(srts,mPi,hadron(SigmaResonance)%mass)

      if (initFlag) then
         open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         open(file=outputFile(3),UNIT=103,Status='Replace',Action='Write')
         open(file=outputFile(4),UNIT=104,Status='Replace',Action='Write')
         write(101,*) '#   srts,    plab,   sigmaTot,   sigmaElast'
         write(102,*) '#   srts,    plab,   KbarN(bg) , KbarN(res)'
         write(103,*) '#   srts,    plab,   KbarN(bg) , KbarN(res)'
         write(104,*) '#   srts,    plab,   KbarN(bg) , KbarN(res)'
         initFlag=.false.
      else
         open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
         open(file=outputFile(3),UNIT=103,Status='old',Position='Append',Action='Write')
         open(file=outputFile(4),UNIT=104,Status='old',Position='Append',Action='Write')
      end if

      call evaluateResContr(sigmaElast_res,sigmaKbarN_res)
      write(101,'(4F9.3)') srts, plab, sigmaTot, sigmaElast_res
      write(102,'(4F9.3)') srts, plab, sigbg, sigmaKbarN_res

      isrts=min(nint((srts-srtmess0)/dmssrt),nssmmax)
      write(103,'(4F9.3)') srts, plab_Lambda, slaka(isrts), slaka_res(isrts)
      write(104,'(4F9.3)') srts, plab_Sigma, sum(ssigka(-4:4,1,isrts))/9., &
                           & sum(ssigka_res(-4:4,1,isrts))/9.
      close(101)
      close(102)
      close(103)
      close(104)
    end subroutine makeOutPut


    !**************************************************************************
    !****s* mesonHyperon/evaluateResContr
    ! NAME
    ! subroutine evaluateResContr(sigmaElast_res,sigmaKbarN_res)
    !
    ! PURPOSE
    ! Computes elastic cross section Y (Y^*) pi -> Y (Y^*) pi
    ! and Kbar N production (final isospin sum) Y (Y^*) pi -> Kbar N
    ! due to intermediate resonances
    ! OUTPUT
    ! * real, intent(out)          :: sigmaElast_res     ! resonance contribution to the elastic Xsection
    ! * real, intent(out)          :: sigmaKbarN_res  ! resonance contribution to the Kbar N production cross section summed over final isospins
    !**************************************************************************

    subroutine evaluateResContr(sigmaElast_res,sigmaKbarN_res)

    use baryonWidth, only: partialWidthBaryon, fullWidthBaryon
    use constants, only: pi
    use clebschGordan

    real, intent(out)          :: sigmaElast_res     ! resonance contribution to the elastic Xsection
    real, intent(out)          :: sigmaKbarN_res  ! resonance contribution to the Kbar N production cross section summed over final isospins
    real :: s,momCM2,gamma_tot,gamma_in,SpinFactor
    real :: iso_res,iso_meson,iso_hyperon,isoFactor
    real :: gamma_out_elastic,gamma_out_KbarN
    integer :: resID

    s=srts**2

    ! Initial c.m. momentum squared:
    momCM2=(s+meson_particle%mass**2-hyperon_particle%mass**2)**2 &
         &/(4.*s) - meson_particle%mass**2

    ! Isospins of incoming meson and hyperon:
    iso_meson=float(hadron(meson_particle%Id)%isoSpinTimes2)/2.
    iso_hyperon=float(hadron(hyperon_particle%Id)%isoSpinTimes2)/2.

    sigmaElast_res=0.
    sigmaKbarN_res=0.

    ! Loop over S=-1 resonances:
    do resID=Sigma_1385,Sigma_1915

      if (.not.hadron(resId)%usedForXsections) cycle   ! Exclude resonances

      ! Total width in vacuum:
      gamma_tot=FullWidthBaryon(resID,srts)

      ! Isospin of the resonance:
      iso_res=float(hadron(resID)%isoSpinTimes2)/2.

      ! In-width:
      gamma_in=partialwidthBaryon(resID,srts,.true.,&
               &meson_particle%Id,hyperon_particle%Id,&
               &meson_particle%mass,hyperon_particle%mass)

      ! Out-widths:
      gamma_out_elastic=partialwidthBaryon(resID,srts,.false.,&
               &meson_particle%Id,hyperon_particle%Id)
      gamma_out_KbarN=partialwidthBaryon(resID,srts,.false.,&
               &kaonBar,nucleon)

      SpinFactor=(2.*hadron(resId)%spin+1.)&
               &/(2.*hadron(hyperon_particle%Id)%spin+1.)

      ! This is valid for particles with isospin 0. and 1.:
      isoFactor=ClebschSquared(iso_hyperon,iso_meson,iso_res,&
               &float(qhyperon),float(qmeson))**2

      sigmaElast_res=sigmaElast_res+isoFactor*SpinFactor*4.*pi/momCM2*s &
                &*gamma_in*gamma_out_elastic &
                &/((s-hadron(resID)%mass**2)**2+gamma_tot**2*s)*0.389

      if (hadron(resId)%propagated) then
        ! This is valid for particles with isospin 0. and 1.:
        isoFactor=ClebschSquared(iso_hyperon,iso_meson,iso_res,&
                 &float(qhyperon),float(qmeson))
        sigmaKbarN_res=sigmaKbarN_res+isoFactor*SpinFactor*4.*pi/momCM2*s &
                  &*gamma_in*gamma_out_KbarN &
                  &/((s-hadron(resID)%mass**2)**2+gamma_tot**2*s)*0.389
      end if

    end do
    ! End of Loop over S=-1 resonances

    end subroutine evaluateResContr

  end subroutine mesonY

end module mesonHyperon
