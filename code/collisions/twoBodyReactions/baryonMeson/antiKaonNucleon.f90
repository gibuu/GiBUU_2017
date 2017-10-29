!******************************************************************************
!****m* /antiKaonNucleon
! NAME
! module antiKaonNucleon
! PURPOSE
! Includes the cross sections for Kbar(Kbar^*)-nucleon and K(K^*)-antinucleon elastic and inelastic scattering
! Public routines:
! * kaonBarNuc
!******************************************************************************
module antiKaonNucleon

  implicit none
  private

  ! To decide wether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.

  public :: kaonBarNuc

contains

  !****************************************************************************
  !****s* antiKaonNucleon/kaonBarNuc
  ! NAME
  ! subroutine kaonBarNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,plotFlag)
  !
  ! PURPOSE
  ! Evaluates Kbar(Kbar^*) N -> Y^*, Kbar N -> Kbar N, Kbar N -> Y pi and
  ! Kbar N -> Y^* pi                 as well as
  ! K(K^*) Nbar -> Ybar^*, K Nbar -> K Nbar, K Nbar -> Ybar pi,
  ! K Nbar -> Ybar^* pi cross sections
  ! and returns also a "preevent"
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
  ! only fixed in the case of a single produced resonance. Otherwise the kinematics still need to be established. The
  ! result is:
  ! * type(preEvent),dimension(1:3), intent(out)               :: teilchenOut    !   outgoing particles
  !
  ! The cross sections are based on parameterization by M. Effenberger.
  ! NOTES
  ! Possible final states are :
  ! * 1-particle : Y^*, Ybar^*
  ! * 2-particle : Kbar N, Y pi, Y^* pi,  K Nbar, Ybar pi, Ybar^* pi
  !****************************************************************************
  !****************************************************************************
  subroutine kaonBarNuc(srts,partIn,mediumAtColl,momLRF,partOut,sigmaTot,sigmaElast,plotFlag)
  use idTable
  use particleDefinition
  use particleProperties, only: hadron
  use mediumDefinition
  use preEventDefinition, only: preEvent
  use twoBodyTools, only: velocity_correction,p_lab,convertToAntiParticles
  use RMF, only: getRMF_flag

!==============================================================================
  ! INPUT:
  real, intent(in)                              :: srts                  ! sqrt(s) in the process
  type(particle),dimension(1:2), intent(in)     :: partIn            ! colliding particles
  type(medium), intent(in)                      :: mediumAtColl     ! Medium informations at the position of the collision
  real, intent(in) ,dimension(0:3)              :: momLRF           ! Total Momentum in LRF
  logical, intent(in),optional                  :: plotFlag    ! Switch on plotting of the  Xsections
!==============================================================================
  ! OUTPUT:
  real, intent(out)                             :: sigmaTot   ! total Xsection
  real, intent(out)                             :: sigmaElast ! elastic Xsection
  type(preEvent),dimension(1:3), intent(out)    :: partOut !   outgoing particles
!==============================================================================
  ! Cross sections:
    real :: sigmaRes_tot ! total cross section Kbar(Kbar^*) N -> Resonances

 !   real , dimension(:),Allocatable :: sigmaRes      ! partial cross sections for Kbar(Kbar^*) N -> R


  ! Field to store the resonance masses

!    real , dimension(:),Allocatable ::massRes      ! Resonance masses

  integer, parameter :: nchnka=4 !*number of channels for kbar n
  real, dimension(0:nchnka) :: sigbg ! background cross sections summed
                                     ! over final isospins:
                                     ! 1 --- Kbar N -> Kbar N,
                                     ! 2 --- Kbar N -> Lambda pi,
                                     ! 3 --- Kbar N -> Sigma pi,
                                     ! 4 --- Kbar N -> Y^* pi
                                     ! 0 --- sum of all above channels


  ! Local variables
  real :: fluxCorrector        ! Correction of the fluxfactor due to different velocities
                               ! in the medium compared to the vacuum
!   real :: s
  type(particle) :: kbar_particle, partNucl
  integer :: qkaon,isrts
  logical :: antiParticleInput
  integer :: totalCharge   ! total charge of colliding particles

  ! Parameters for discretisation of sqrts (same as in Effe's code):
  real, parameter :: srtmess0=1.434     !* srtmess0=m_kaon+m_nucleon    ! 1.254     !*srtmess0=m_lambda+m_pion
  real, parameter :: dmssrt=0.001       ! sqrts-bin size
  integer, parameter :: nssmmax=1500    ! number of sqrts-bins

  ! Fields to store discretised background cross sections on a proton
  ! (same as in Effe's code) and respective resonance contributions
  ! to the given channels:

  !     Kbar N -> Kbar N:
  !     first index --- charge of incoming Kbar,
  !     second index = -1 --- outgoing K^-, 0 --- outgoing Kbar^0,
  !                     1 --- sum over final isospins
  real, dimension(-1:0,-1:1,0:nssmmax), save :: skaka,skaka_res

  !     K^- p -> Lambda pi^0:
  real, dimension(0:nssmmax), save :: skapi,skapi_res

  !     Kbar N -> Sigma pi:
  !     first index --- charge of incoming Kbar,
  !     second index = -1 --- outgoing pi^-, 0 --- outgoing pi^0,
  !                    +1 --- outgoing pi^+, 2 --- sum over final isospins
  real, dimension(-1:0,-1:2,0:nssmmax), save :: skaspi2,skaspi2_res

  !     Kbar N -> Y^* pi:
  !     first index --- charge of incoming Kbar,
  !     third index = j, j=1-nsres --- partial cross section for Y^* with
  !                                    ID=1+nres+j summed over final isospins
  !                 = 0 --- sum of all partial Y^* cross sections
  real, dimension(-1:0,0:nssmmax,0:nsres), save :: skastar


  !kBar N --> Xi K:
  real, dimension(1:5) :: sigXi

  logical, save :: debug=.false.
  logical, save :: flagini=.true.

  ! partial cross sections for Kbar(Kbar^*) N -> R
  real, dimension(Delta:nbar) :: sigmaRes
  ! Field to store the resonance masses
  real, dimension(Delta:nbar) :: massRes      ! Resonance masses

  if (flagini) then
    ! Discretization of cross sections:
    call init
    flagini=.false.
  end if

  ! Initialize output
  partOut(:)%ID=0                    ! ID of produced particles
  partOut(:)%charge=0                ! Charge of produced particles
  partOut(:)%antiParticle=.false.    ! Whether produced particles are
                                         ! particles or antiparticles
  partOut(:)%mass=0.                 ! Mass of produced particles

  ! (1) Check  Input
  if (partIn(1)%Id.ge.kaon.and.partIn(1)%Id.le.kaonStarBar&
     &.and.partIn(2)%Id.eq.nucleon) then
    kbar_particle=partIn(1)  ;  partNucl=partIn(2)
  else if (partIn(2)%Id.ge.kaon.and.partIn(2)%Id.le.kaonStarBar&
     &.and.partIn(1)%Id.eq.nucleon) then
    kbar_particle=partIn(2)  ;  partNucl=partIn(1)
  else
    write(*,*) 'Wrong input in kaonBarNuc', partIn%ID
    stop
  end if

  if (kbar_particle%antiParticle) then
    ! This must not happen:
    write(*,*) 'kbar is antiparticle in "kaonBarNuc"!!!',&
               &partIn%ID,partIn%antiparticle
    stop
  end if

  if ((kbar_particle%Id==kaon.or.kbar_particle%Id==kaonStar)&
    &.and. partNucl%antiParticle) then
    ! Invert all particles in antiparticles:
    partNucl%Charge=-partNucl%Charge
    partNucl%antiparticle=.false.
    kbar_particle%Charge=-kbar_particle%Charge
    kbar_particle%Id=kbar_particle%Id+1
    antiParticleInput=.true.
  else if ((kbar_particle%Id==kaon.or.kbar_particle%Id==kaonStar)&
    &.or. partNucl%antiParticle) then
    write(*,*) 'In kaonBarNuc: K N and Kbar Nbar collisions must be treated by subroutine kaonNuc !', partIn%ID
    stop
  else
    antiParticleInput=.false.
  end if

  totalCharge=kbar_particle%Charge+partNucl%Charge

  ! Correction of the fluxfactor due to different velocities in the medium compared to the vacuum
  if ( .not.getRMF_flag() ) then
    fluxCorrector=velocity_correction(partIn)
  else
    fluxCorrector=1.
  end if

!   s=srts**2

  ! (2) Evaluate the cross sections
  call evaluateXsections

  ! Cutoff to kick the case out, that the cross section is zero
  if (sigmaTot.lt.1E-12) then
    sigmatot=0.
    sigmaElast=0.
    return
  end if

  ! (3) Plot them if wished
  if (present(PlotFlag).or.debug) then
    if (plotFlag.or.debug)  call makeOutput
  end if


  ! (4) Define final state
  call MakeDecision

  ! (5) Check Output
  if (Sum(partOut(:)%Charge).ne.partNucl%charge+kbar_particle%charge) then
    write(*,*) 'No charge conservation in kaonBarNuc!!! Critical error', &
               kbar_particle%ID, kbar_particle%Charge, partNucl%Charge, partOut(:)%Charge, partOut(:)%ID
    stop
  end if

  ! (6) Invert particles in antiParticles if input included antiparticles
  if (antiParticleInput) then
     if (debug) write(*,*) partOut
     call convertToAntiParticles(partOut)
     if (debug) write(*,*) partOut
  end if

  contains

    !**************************************************************************
    !****s* kaonBarNuc/evaluateXsections
    ! NAME
    ! subroutine evaluateXsections
    !
    ! PURPOSE
    ! Evaluates Kbar(Kbar^*) N -> Y^*, Kbar N -> Kbar N, Kbar N -> Y pi and
    ! Kbar N -> Y^* pi cross sections
    !
    ! NOTES
    ! There is a resonance contribution to Kbar(Kbar^*) N scattering.
    !**************************************************************************

    subroutine evaluateXsections
    use parametrizationsBarMes, only: sigma_KbarToXi, kaonbg
    use resonanceCrossSections, only: barMes2resonance
    use constants, only: mPi, mK
    use twoBodyTools, only: pCM

    real, dimension(1:3) :: position
    logical :: perturbative
    integer :: i
    real :: p_KN, p_KstarN, ratio

    position=0.5*(partIn(1)%position+partIn(2)%position)
    if (partIn(1)%perturbative.or.partIn(2)%perturbative) then
      perturbative=.true.
    else
      perturbative=.false.
    end if

    ! Resonance contribution:
    sigmaRes = barMes2resonance (kbar_particle%Id,partNucl%Id,kbar_particle%Charge,partNucl%Charge, &
                                 .true.,mediumAtColl,momLRF,massRes,kbar_particle%mass,partNucl%mass, &
                                 position,perturbative,srts)
    sigmaRes_tot=sum(sigmaRes(:))

    if (debug) write(*,*) ' sigmaRes_tot= ', sigmaRes_tot

    ! Non-resonant contributions:

    ratio=1.
    if (kbar_particle%Id.eq.kaonStarBar) then
       ! Correction factor to the antikaon-nucleon cross sections by detailed balance
       ! (see J. Cugnon et al, PRC 40, 1822 (1989))
       ! c.m. momenta of kBar-nucleon and kaonStarBar-nucleon
       p_KN     = pCM(srts,mK,partNucl%Mass)
       p_KstarN = pCM(srts,kbar_particle%Mass,partNucl%Mass)
       if (p_KstarN.gt.1.e-06) ratio=p_KN/p_KstarN
    end if

    isrts=max(2,min(nint((srts-srtmess0)/dmssrt),nssmmax))

    if (debug) write(*,*) ' srts, isrts: ', srts, isrts

    ! Kbar N -> X:

    qkaon=kbar_particle%Charge
    !*isospin rotation of antikaon (charge: 0->-1,-1->0):
    if (partNucl%Charge.eq.0) qkaon=-qkaon-1

    if (debug) write(*,*) ' qkaon= ', qkaon

    sigbg(:)=0.

    !*Kbar N -> Kbar N:
    if (isrts.eq.2) then
       if (qkaon.eq.-1) then
          skaka(-1,-1,isrts)=kaonbg(srts,1,1,0)
          skaka(-1,0,isrts)=kaonbg(srts,1,2,0)
          skaka(-1,1,isrts)=skaka(-1,-1,isrts)+skaka(-1,0,isrts)
       else
          skaka(0,0,isrts)=kaonbg(srts,1,1,1)
          skaka(0,-1,isrts)=0.
          skaka(0,1,isrts)=skaka(0,-1,isrts)+skaka(0,0,isrts)
       end if
    end if
    sigbg(1)=skaka(qkaon,1,isrts)

    !*Kbar N -> Lambda pi:
    if (isrts.eq.2) skapi(isrts)=kaonbg(srts,1,3,0)
    if (qkaon.eq.-1) then
       sigbg(2)=skapi(isrts)
    else
       !*factor 2 is isospin coefficient:
       sigbg(2)=2.*skapi(isrts)
    end if

    !*Kbar N -> Sigma pi:
    if (isrts.eq.2) then
       if (qkaon.eq.-1) then
              !*K- p -> Sigma+ pi-:
          skaspi2(-1,-1,isrts)=kaonbg(srts,1,4,0)
              !*K- p -> Sigma- pi+:
          skaspi2(-1,1,isrts)=kaonbg(srts,1,5,0)
              !*K- p -> Sigma0 pi0:
          skaspi2(-1,0,isrts)=kaonbg(srts,1,6,0)
              !*sum:
          skaspi2(-1,2,isrts)=skaspi2(-1,-1,isrts)+skaspi2(-1,1,isrts)+skaspi2(-1,0,isrts)
       else
              !*K0 p -> Sigma+ pi0:
          skaspi2(0,0,isrts)=kaonbg(srts,1,4,1)
              !*K0 p -> Sigma0 pi+:
          skaspi2(0,1,isrts)=skaspi2(0,0,isrts)
              !*K0 p -> ? pi-:
          skaspi2(0,-1,isrts)=0.
              !sum K0 p -> Sigma pi:
          skaspi2(0,2,isrts)=2.*skaspi2(0,1,isrts)
       end if
    end if
    sigbg(3)=skaspi2(qkaon,2,isrts)

    !*additional background for higher energies (Kbar N -> Y^* pi):
    if (srts.gt.hadron(Lambda)%mass+2.*mPi) &
       & sigbg(4)=skastar(qkaon,isrts,0)

    sigbg(0)=sum(sigbg(1:4))
    sigbg=sigbg*ratio

    if (debug) write(*,*) ' sigbg: ', sigbg(0:4)

    ! Kbar N -> Xi K
    sigXi(:) = 0.0
    if (qkaon==-1 .and. partNucl%Charge==1) then
       do i=1,5
          sigXi(i) = sigma_KbarToXi(srts,i)
       end do
    end if
    sigXi=sigXi*ratio

    if (debug) write(*,*) '  sigXi: ',  sigXi(1:5)

    ! Elastic and total cross sections:
    if (kbar_particle%Id.eq.kaonBar) then
      sigmaElast=skaka(qkaon,qkaon,isrts)+skaka_res(qkaon,qkaon,isrts)
    else
      ! NOTE: the Kbar^* N -> Kbar^* N elastic Xsection is only due to
      ! the resonance contribution.
      ! Since it is not needed we put it to zero to save CPU time.
      sigmaElast=0.
    end if
    sigmaTot=sigbg(0)+sigmaRes_tot+Sum(sigXi(:))

    if (debug) write(*,*) ' sigmaTot= ', sigmaTot
    if (debug) write(*,*) ' fluxcorrector= ', fluxcorrector

    ! Flux correction for each channel:
    if (fluxCorrector_flag) then
      sigmaRes=sigmaRes*fluxcorrector
      sigmaRes_tot=sigmaRes_tot*fluxcorrector
      sigbg=sigbg*fluxcorrector
      sigXi=sigXi*fluxcorrector
      sigmaElast=sigmaElast*fluxcorrector
      sigmaTot=sigmaTot*fluxcorrector
    end if


    end subroutine evaluateXsections



    !**************************************************************************
    !****s* kaonBarNuc/init
    ! NAME
    ! subroutine init
    !
    ! PURPOSE
    ! Discretises the backgrond cross sections for Kbar-nucleon scattering.
    !**************************************************************************

    subroutine init

    use baryonWidth, only: partialWidthBaryon, fullWidthBaryon
    use parametrizationsBarMes, only: kaonbg
    use constants, only: pi, mN, mPi, mK
    use clebschGordan, only: ClebschSquared
    use twoBodyPhaseSpace, only: Integrate_2bodyPS_resonance
    use twoBodyTools, only: pCM

    real, parameter :: m2yst_km=22., m2yst_k0=15. ! (mb GeV^2) matrix element squared for the Kbar N -> Y^* pi (see Effe's PhD, p. 239)

    real :: srts,s,momCM,gamma_tot,gamma_in,SpinFactor,sigma_res
    real :: iso_res,isoFactor,gamma_out
    real, dimension(1:5) :: ps
    integer :: i,j,resID!,k

    skaka=0.
    skaka_res=0.
    skapi=0.
    skapi_res=0.
    skaspi2=0.
    skaspi2_res=0.
    skastar=0.

    Loop_over_srts : do i=0,nssmmax

      srts=float(i)*dmssrt+srtmess0
      if (srts <= mN + mK) cycle
      s=srts**2

!*K- p -> K- p:
      skaka(-1,-1,i)=kaonbg(srts,1,1,0)
!*K- p -> K0 n:
      skaka(-1,0,i)=kaonbg(srts,1,2,0)
!*sum:
      skaka(-1,1,i)=skaka(-1,-1,i)+skaka(-1,0,i)

!*K0 p -> K0 p:
      skaka(0,0,i)=kaonbg(srts,1,1,1)
      skaka(0,-1,i)=0.
      skaka(0,1,i)=skaka(0,-1,i)+skaka(0,0,i)

!*K- p -> Lambda pi0:
      skapi(i)=kaonbg(srts,1,3,0)

!*K- p -> Sigma+ pi-:
      skaspi2(-1,-1,i)=kaonbg(srts,1,4,0)
!*K- p -> Sigma- pi+:
      skaspi2(-1,1,i)=kaonbg(srts,1,5,0)
!*K- p -> Sigma0 pi0:
      skaspi2(-1,0,i)=kaonbg(srts,1,6,0)
!*sum:
      skaspi2(-1,2,i)=skaspi2(-1,-1,i)+skaspi2(-1,1,i)+skaspi2(-1,0,i)

!*K0 p -> Sigma+ pi0:
      skaspi2(0,0,i)=kaonbg(srts,1,4,1)
!*K0 p -> Sigma0 pi+:
      skaspi2(0,1,i)=skaspi2(0,0,i)
!*K0 p -> ? pi-:
      skaspi2(0,-1,i)=0.
!sum K0 p -> Sigma pi:
      skaspi2(0,2,i)=2.*skaspi2(0,1,i)

      momCM = pCM(srts,mN,mK)  ! Initial c.m. momentum squared

      ! Resonance contribution to various channels:
      ! NOTE: if resonance is not explicitly propagated, these contributions
      ! are added to the background fields

      ! Loop over S=-1 resonances:
      do resID=Sigma_1385,Sigma_1915

        if (.not.hadron(resId)%usedForXsections) cycle   ! Exclude resonances

        ! Total width in vacuum:
        gamma_tot=FullWidthBaryon(resID,srts)

        ! In-width:
        gamma_in=partialwidthBaryon(resID,srts,.true.,kaonBar,nucleon)

        SpinFactor=(2.*hadron(resId)%spin+1.)/2.

        sigma_res=SpinFactor*4.*pi/momCM**2*s*gamma_in &
                &/((s-hadron(resID)%mass**2)**2+gamma_tot**2*s)*0.389

        ! Isospin of the resonance:
        iso_res=float(hadron(resID)%isoSpinTimes2)/2.

        !*K- p -> K- p:
        isoFactor=ClebschSquared(0.5,0.5,iso_res,-0.5,0.5)**2
        gamma_out=partialwidthBaryon(resID,srts,.false.,kaonBar,nucleon)
        if (hadron(resId)%propagated) then
          skaka_res(-1,-1,i)=skaka_res(-1,-1,i)+isoFactor*sigma_res*gamma_out
        else
          skaka(-1,-1,i)=skaka(-1,-1,i)+isoFactor*sigma_res*gamma_out
        end if

        !*K- p -> K0 n:
        if (hadron(resId)%propagated) then
          skaka_res(-1,0,i)=skaka_res(-1,0,i)+isoFactor*sigma_res*gamma_out
        else
          skaka(-1,0,i)=skaka(-1,0,i)+isoFactor*sigma_res*gamma_out
        end if

        !*K0 p -> K0 p:
        if (abs(iso_res-1.).lt.0.00001) then
          isoFactor=1.
          if (hadron(resId)%propagated) then
            skaka_res(0,0,i)=skaka_res(0,0,i)+isoFactor*sigma_res*gamma_out
          else
            skaka(0,0,i)=skaka(0,0,i)+isoFactor*sigma_res*gamma_out
          end if
        end if

        !*K- p -> Lambda pi0:
        if (abs(iso_res-1.).lt.0.00001) then
          isoFactor=ClebschSquared(0.5,0.5,iso_res,-0.5,0.5) &
                  &*ClebschSquared(0.,1.,iso_res,0.,0.)
          gamma_out=partialwidthBaryon(resID,srts,.false.,pion,Lambda)
          if (hadron(resId)%propagated) then
            skapi_res(i)=skapi_res(i)+isoFactor*sigma_res*gamma_out
          else
            skapi(i)=skapi(i)+isoFactor*sigma_res*gamma_out
          end if
        end if

        !*K- p -> Sigma+ pi-:
        isoFactor=ClebschSquared(0.5,0.5,iso_res,-0.5,0.5) &
                &*ClebschSquared(1.,1.,iso_res,1.,-1.)
        gamma_out=partialwidthBaryon(resID,srts,.false.,pion,SigmaResonance)
        if (hadron(resId)%propagated) then
          skaspi2_res(-1,-1,i)=skaspi2_res(-1,-1,i) &
                             &+isoFactor*sigma_res*gamma_out
        else
          skaspi2(-1,-1,i)=skaspi2(-1,-1,i)+isoFactor*sigma_res*gamma_out
        end if

        !*K- p -> Sigma- pi+:
        if (hadron(resId)%propagated) then
          skaspi2_res(-1,1,i)=skaspi2_res(-1,1,i)+isoFactor*sigma_res*gamma_out
        else
          skaspi2(-1,1,i)=skaspi2(-1,1,i)+isoFactor*sigma_res*gamma_out
        end if


        !*K- p -> Sigma0 pi0:
        isoFactor=ClebschSquared(0.5,0.5,iso_res,-0.5,0.5) &
                &*ClebschSquared(1.,1.,iso_res,0.,0.)
        if (hadron(resId)%propagated) then
          skaspi2_res(-1,0,i)=skaspi2_res(-1,0,i)+isoFactor*sigma_res*gamma_out
        else
          skaspi2(-1,0,i)=skaspi2(-1,0,i)+isoFactor*sigma_res*gamma_out
        end if


        if (abs(iso_res-1.).lt.0.00001) then

          !*K0 p -> Sigma+ pi0:
          isoFactor=ClebschSquared(1.,1.,iso_res,1.,0.)
          if (hadron(resId)%propagated) then
            skaspi2_res(0,0,i)=skaspi2_res(0,0,i)+isoFactor*sigma_res*gamma_out
          else
            skaspi2(0,0,i)=skaspi2(0,0,i)+isoFactor*sigma_res*gamma_out
          end if

          !*K0 p -> Sigma0 pi+:
          if (hadron(resId)%propagated) then
            skaspi2_res(0,1,i)=skaspi2_res(0,1,i)+isoFactor*sigma_res*gamma_out
          else
            skaspi2(0,1,i)=skaspi2(0,1,i)+isoFactor*sigma_res*gamma_out
          end if

        end if

      end do
      ! End of loop over S=-1 resonances

      !*isospin sums:
      skaka(-1,1,i)=skaka(-1,-1,i)+skaka(-1,0,i)
      skaka_res(-1,1,i)=skaka_res(-1,-1,i)+skaka_res(-1,0,i)
      skaka(0,1,i)=skaka(0,-1,i)+skaka(0,0,i)
      skaka_res(0,1,i)=skaka_res(0,-1,i)+skaka_res(0,0,i)
      skaspi2(-1,2,i)=skaspi2(-1,-1,i)+skaspi2(-1,0,i)+skaspi2(-1,1,i)
      skaspi2_res(-1,2,i)=skaspi2_res(-1,-1,i)+skaspi2_res(-1,0,i) &
                        &+skaspi2_res(-1,1,i)
      skaspi2(0,2,i)=skaspi2(0,-1,i)+skaspi2(0,0,i)+skaspi2(0,1,i)
      skaspi2_res(0,2,i)=skaspi2_res(0,-1,i)+skaspi2_res(0,0,i) &
                       &+skaspi2_res(0,1,i)


      !*Kbar N -> Ystar pi:
      skastar(-1,i,0)=0.
      skastar(0,i,0)=0.
      do resID=lambda_1600,sigma_1915

        ps = Integrate_2bodyPS_resonance (resID, srts, mPi, 0.)
        iso_res=float(hadron(resID)%isoSpinTimes2)/2.
        j=resID-nres-1
        if (j < 1 .or. j > nsres) then
          write(*,*) 'In antiKaonNucleon, init: j out of bounds', j, nsres
          stop
        end if

        if (abs(iso_res).lt.0.00001) then
          skastar(-1,i,j)=0.5*m2yst_km*ps(1)/momCM/s
          skastar(0,i,j)=m2yst_k0*ps(1)/momCM/s
        else
          skastar(-1,i,j)=m2yst_km*ps(1)/momCM/s
          skastar(0,i,j)=m2yst_k0*ps(1)/momCM/s
        end if
        skastar(-1,i,0)=skastar(-1,i,0)+skastar(-1,i,j)
        skastar(0,i,0)=skastar(0,i,0)+skastar(0,i,j)

      end do

    end do Loop_over_srts

    end subroutine init


    !**************************************************************************
    !****s* kaonBarNuc/makeDecision
    ! NAME
    ! subroutine makeDecision
    ! PURPOSE
    ! Chooses randomly one of possible outgoing channels in antikaon-nucleon
    ! collision. Outgoing channels are: Y^*, Kbar N, Y pi, Y^* pi.
    ! Also the charges of outgoing particles are selected.
    !**************************************************************************
    subroutine makeDecision

      use random, only: rn

      real :: cut,cut2,sum,x
      integer :: resID,qkaon2,qpion2,j

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      if (sigmaRes_tot >= cut) then

        ! Resonance production:
        sum=0.
        do resID=Delta,nbar
          sum=sum+sigmaRes(resID)
          if (sum >= cut) then
            partOut(1)%Id=resID
            partOut(1)%Charge=totalCharge
            partOut(1)%Mass=massRes(resID)
            exit
          end if
        end do

      else if (sigmaRes_tot+sigbg(1) >= cut) then

        ! Kbar N production (elastic scattering):
        sum=0.
        cut2=rn()*skaka(qkaon,1,isrts)
        do qkaon2=-1,0
          sum=sum+skaka(qkaon,qkaon2,isrts)
          if (sum >= cut2) exit
        end do
        partOut(1)%Id=kaonBar
        partOut(2)%Id=nucleon
        ! isospin rotation of antikaon:
        if (partNucl%Charge.eq.0) qkaon2=-qkaon2-1
        partOut(1)%Charge=qkaon2
        partOut(2)%Charge=totalCharge-qkaon2

      else if (sigmaRes_tot+sigbg(1)+sigbg(2) >= cut) then

        ! Lambda pi production:
        partOut(1)%Id=pion
        partOut(2)%Id=Lambda
        partOut(1)%Charge=totalCharge
        partOut(2)%Charge=0

      else if (sigmaRes_tot+sigbg(1)+sigbg(2)+sigbg(3) >= cut) then

        ! Sigma pi production:
        sum=0.
        cut2=rn()*skaspi2(qkaon,2,isrts)
        do qpion2=-1,1
          sum=sum+skaspi2(qkaon,qpion2,isrts)
          if (sum >= cut2) exit
        end do
        partOut(1)%Id=pion
        partOut(2)%Id=SigmaResonance
        ! isospin rotation of pion:
        if (partNucl%Charge.eq.0) qpion2=-qpion2
        partOut(1)%Charge=qpion2
        partOut(2)%Charge=totalCharge-qpion2

      else if (sigmaRes_tot+sigbg(0) >= cut) then

        ! Y^* pi production:

        !*select resonance:
        sum=0.
        cut2=rn()*skastar(qkaon,isrts,0)
        do j=1,nsres
          sum=sum+skastar(qkaon,isrts,j)
          if (sum >= cut2) exit
        end do
        partOut(1)%Id=pion
        partOut(2)%Id=j+nres+1

        !*determine charges:
        if (hadron(partOut(2)%Id)%isoSpinTimes2==0) then
          ! Lambda^* production:
          partOut(1)%Charge=totalCharge
          partOut(2)%Charge=0
        else
          ! Sigma^* production:
          x=rn()
          if (qkaon.eq.-1) then
            if (x.le.1./6.) then
              qpion2=0
            else if (x.le.7./12.) then
              qpion2=-1
            else
              qpion2=1
            end if
          else if (qkaon.eq.0) then
            if (x.le.0.5) then
              qpion2=1
            else
              qpion2=0
            end if
          end if
          !*isospin rotation of pion:
          if (partNucl%Charge.eq.0) qpion2=-qpion2
          partOut(1)%Charge=qpion2
          partOut(2)%Charge=totalCharge-qpion2
        end if

     else if (sigmaRes_tot+sigbg(0) + sigXi(1) >= cut) then

        partOut(1)%ID = Xi
        partOut(2)%ID = kaon
        partOut(1)%charge=-1
        partOut(2)%charge=1

     else if (sigmaRes_tot+sigbg(0) + sigXi(1) + sigXi(2) >= cut) then

        partOut(1)%ID = Xi
        partOut(2)%ID = kaon
        partOut(1)%charge=0
        partOut(2)%charge=0

     else if (sigmaRes_tot+sigbg(0) + sigXi(1) + sigXi(2) + sigXi(3) >= cut) then

        partOut(1)%ID = Xi
        partOut(2)%ID = kaon
        partOut(3)%ID = pion
        partOut(1)%charge=-1
        partOut(2)%charge=0
        partOut(3)%charge=1

     else if (sigmaRes_tot+sigbg(0) + sigXi(1) + sigXi(2) + sigXi(3) + sigXi(4) >= cut) then

        partOut(1)%ID = Xi
        partOut(2)%ID = kaon
        partOut(3)%ID = pion
        partOut(1)%charge=-1
        partOut(2)%charge=1
        partOut(3)%charge=0

     else if (sigmaRes_tot+sigbg(0) + sigXi(1) + sigXi(2) + sigXi(3) + sigXi(4) + sigXi(5) >= cut) then

        partOut(1)%ID = Xi
        partOut(2)%ID = kaon
        partOut(3)%ID = pion
        partOut(1)%charge=0
        partOut(2)%charge=1
        partOut(3)%charge=-1

      else

        ! No event was generated:
        write(*,*) 'Error in makedecision of kaonBarNuc', sigmaRes_tot,&
                   & sigbg,sigmaRes_tot+sigbg(0),sigmaTot,cut

        write(*,*) sigXi
        write(*,*) srts
        stop

      end if

    end subroutine makeDecision


    !**************************************************************************
    !****s* kaonBarNuc/makeOutput
    ! NAME
    ! subroutine makeOutput
    !
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV]
    ! .
    ! Filenames:
    ! * 'KbarN_sigTotElast.dat'  : sigmaTot, sigmaElast
    ! * 'KbarN_KmN.dat'          : K^- Nucleon outgoing channel
    ! * 'KbarN_Kbar0N.dat'       : Kbar^0 Nucleon outgoing channel
    ! * 'KbarN_LambdaPi.dat'     : Lambda pion outgoing channel
    ! * 'KbarN_SigmaPip.dat'     : Sigma pi^+ outgoing channel
    ! * 'KbarN_SigmaPi0.dat'     : Sigma pi^0 outgoing channel
    ! * 'KbarN_SigmaPim.dat'     : Sigma pi^- outgoing channel
    ! * 'KbarN_YstPi.dat'        : Y^* pion outgoing channel summed
    ! *                          : over all final hyperon resonances
    ! *                          : with M >= 1600 MeV
    ! NOTES
    ! for the Kbar^* N input channel only total resonance X section
    ! is printed out
    !**************************************************************************
    subroutine makeOutPut

      logical, save :: initFlag=.true.

      ! The output files
      character(30), dimension(1:8) :: outputFile
      real :: plab
      integer :: qkaon2,qpion2

      outputFile(1)='KbarN_sigTotElast.dat'
      outputFile(2)='KbarN_KmN.dat'
      outputFile(3)='KbarN_Kbar0N.dat'
      outputFile(4)='KbarN_LambdaPi.dat'
      outputFile(5)='KbarN_SigmaPip.dat'
      outputFile(6)='KbarN_SigmaPi0.dat'
      outputFile(7)='KbarN_SigmaPim.dat'
      outputFile(8)='KbarN_YstPi.dat'

      plab=p_lab(srts,kbar_particle%mass,partNucl%mass)

      if (initFlag) then
         open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         open(file=outputFile(3),UNIT=103,Status='Replace',Action='Write')
         open(file=outputFile(4),UNIT=104,Status='Replace',Action='Write')
         open(file=outputFile(5),UNIT=105,Status='Replace',Action='Write')
         open(file=outputFile(6),UNIT=106,Status='Replace',Action='Write')
         open(file=outputFile(7),UNIT=107,Status='Replace',Action='Write')
         open(file=outputFile(8),UNIT=108,Status='Replace',Action='Write')
         write(101,*) '#   srts,    plab,   sigmaTot,sigmaElast'
         write(102,*) '#   srts,    plab,   K^-N(bg) ,  K^-N(res)'
         write(103,*) '#   srts,    plab,   Kbar^0N(bg) , Kbar^0N(res)'
         write(104,*) '#   srts,    plab,   LambdaPi(bg) , LambdaPi(res)'
         write(105,*) '#   srts,    plab,   SigmaPip(bg) , SigmaPip(res)'
         write(106,*) '#   srts,    plab,   SigmaPi0(bg) , SigmaPi0(res)'
         write(107,*) '#   srts,    plab,   SigmaPim(bg) , SigmaPim(res)'
         write(108,*) '#   srts,    plab,   Y^*Pi'
         initFlag=.false.
      else
         open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
         open(file=outputFile(3),UNIT=103,Status='old',Position='Append',Action='Write')
         open(file=outputFile(4),UNIT=104,Status='old',Position='Append',Action='Write')
         open(file=outputFile(5),UNIT=105,Status='old',Position='Append',Action='Write')
         open(file=outputFile(6),UNIT=106,Status='old',Position='Append',Action='Write')
         open(file=outputFile(7),UNIT=107,Status='old',Position='Append',Action='Write')
         open(file=outputFile(8),UNIT=108,Status='old',Position='Append',Action='Write')
      end if

      write(101,'(4F10.4)') srts, plab, sigmaTot, sigmaElast

      if (kbar_particle%Id.eq.kaonBar) then

        if (partNucl%charge==1) then
          qkaon2=-1
          qpion2=1
        else
          qkaon2=0
          qpion2=-1
        end if
        write(102,'(4F10.4)') srts, plab, skaka(qkaon,qkaon2,isrts),&
                            & skaka_res(qkaon,qkaon2,isrts)
        write(103,'(4F10.4)') srts, plab, skaka(qkaon,-qkaon2-1,isrts),&
                            & skaka_res(qkaon,-qkaon2-1,isrts)
        if (totalCharge==0) then
          write(104,'(4F10.4)') srts, plab, skapi(isrts), skapi_res(isrts)
        else
          write(104,'(4F10.4)') srts, plab, 2.*skapi(isrts), 2.*skapi_res(isrts)
        end if
        write(105,'(4F10.4)') srts, plab, skaspi2(qkaon,qpion2,isrts), &
                            & skaspi2_res(qkaon,qpion2,isrts)
        write(106,'(4F10.4)') srts, plab, skaspi2(qkaon,0,isrts), &
                            & skaspi2_res(qkaon,0,isrts)
        write(107,'(4F10.4)') srts, plab, skaspi2(qkaon,-qpion2,isrts), &
                            & skaspi2_res(qkaon,-qpion2,isrts)
        write(108,'(4F10.4)') srts, plab, skastar(qkaon,isrts,0)

      end if

      close(101)
      close(102)
      close(103)
      close(104)
      close(105)
      close(106)
      close(107)
      close(108)
    end subroutine makeOutPut

  end subroutine kaonBarNuc

end module antiKaonNucleon
