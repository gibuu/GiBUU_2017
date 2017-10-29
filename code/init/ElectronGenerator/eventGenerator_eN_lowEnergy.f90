!******************************************************************************
!****m* /eventGenerator_eN_lowEnergy
! NAME
! module eventGenerator_eN_lowEnergy
!
! PURPOSE
! This module includes initilization routines for low energetic electron
! induced events
! (Resonance region and quasi-elastic regime, i.e. 0.7<W<2 GeV)
!
! INPUTS
! no namelist. All parameters are given as function arguments.
!******************************************************************************
module eventGenerator_eN_lowEnergy

  use CALLSTACK

  implicit none

  private

  !****************************************************************************
  !****g* eventGenerator_eN_lowEnergy/iC_XXX
  ! PURPOSE
  ! Constants defined for ease of coding: possible Channels
  !
  ! SOURCE
  !
  integer, parameter  :: iC_QE     = 1
  integer, parameter  :: iC_Res    = 2
  integer, parameter  :: iC_1Pi    = 3
  integer, parameter  :: iC_2Pi    = 4
  integer, parameter  :: iC_DIS    = 5
  integer, parameter  :: iC_VMDrho = 6
  integer, parameter  :: iC_2p2hQE = 7
  integer, parameter  :: iC_2p2hDelta = 8

  integer, parameter  :: nC        = 8
  !****************************************************************************

  logical,save :: initFlag=.true.

  public :: eventGen_eN_lowEnergy
!  PUBLIC :: iC_QE, iC_Res, iC_1Pi, iC_2Pi, iC_DIS, iC_VMDrho
  public :: init_2Pi_getBG, init_2Pi
  public :: init_DIS
  public :: checkEvent

contains

  subroutine initInput()

    initFlag=.false.

    ! --- NO CONTENT YET ---
  end subroutine initInput


  !****************************************************************************
  !****s* eventGenerator_eN_lowEnergy/eventGen_eN_lowEnergy
  ! NAME
  ! subroutine eventGen_eN_lowEnergy(eN,doC,whichRes,DoPauli,OutPart,channel,flagOK,XS,XS_Arr)
  !
  ! PURPOSE
  ! Initializes one e^- N -> X events for 0.7<W<2 GeV .
  ! The weight each event is given by (total Xsection of event=sum over all
  ! possible channels).
  ! We make a Monte Carlo decision to choose a specific channel, which is
  ! returned as "OutPart".
  !
  ! The particles are produced at the place of the nucleon target.
  !
  ! INPUTS
  ! Scattering particles:
  ! * type(electronNucleon_event)  :: eN -- The incoming electron and nucleon
  ! * logical, dimension(:)   :: doC -- switches, see below
  ! * logical :: DoPauli -- if .true., every event is checked for pauli blocking.
  !   if it is blocked, the corresponding cross section will be set to 0 and the
  !   Monte Carlo decision will neglect this event class. Thus the total cross section
  !   is reduced. Maybe all channel can be closed.
  ! * type(particle),dimension(:,:) :: realParticles -- real particle vector,
  !   only needed for pauli blocking
  !
  ! Switches (cf. iC_'module Electron_origin' and XXX):
  ! * QE     -- Switch on/off quasi-elastic events
  ! * Res    -- Switch on/off resonance production events
  ! * 1Pi    -- Switch on/off single pion production events
  ! * 2Pi    -- Switch on/off 2 pion background events
  ! * DIS    -- Switch on/off DIS events
  ! * VMDrho -- Switch on/off special treatment of gamma N -> rho0 N
  ! * 2p2hQE -- Switch on/off gamma N N --> N' N'
  ! * 2p2hDelta -- Switch on/off gamma N N --> Delta N'
  !
  ! Special Switch for resonance production:
  ! * logical, dimension(2:nres+1), intent(in)  :: whichRes
  ! With this switch special resonances can be selected, if res_flag=.true.,
  ! and whichRes is defined by the lines...
  !
  !   whichRes=.false.
  !   whichRes(Delta)=.true.
  !
  ! ... then only the Delta is included as a possible resonance channel.
  !
  ! OUTPUT
  ! * integer                      :: channel     --
  !   value according to naming scheme defined in module "Electron_origin".
  !   For all 1-Body final states the value of "channel" is given by the ID
  !   of the produced particle.
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * logical                      :: flagOK -- .true. if OutPart was created
  ! * real, OPTIONAL               :: XS -- cross section in mub
  ! * real,dimension(...),OPTIONAL :: XS_Arr -- cross sections according the
  !   internal weights (in mub, maybe also negative)
  !****************************************************************************
  subroutine eventGen_eN_lowEnergy(eN,doC,whichRes,DoPauli,realparticles,OutPart,channel,flagOK,XS,XS_Arr)
    use particleDefinition
    use eN_eventDefinition, only: electronNucleon_event
    use Electron_origin, only: origin_singlePi,origin_doublePi,origin_DIS,origin_vecmes_rho,origin_2p2hQE,origin_2p2hDelta
    use idTable, only: nres
    use monteCarlo, only: monteCarloChoose
    use pauliBlockingModule, only: checkPauli

    type(electronNucleon_event), intent(inout):: eN
    logical, dimension(nC),      intent(in)   :: doC

    logical, dimension(2:nres+1), intent(in)  :: whichRes
    logical,                      intent(in)  :: DoPauli
    type(particle), intent(in), dimension(:,:):: realParticles

    integer,        intent(out)               :: channel
    type(particle), intent(out), dimension(:) :: OutPart
    logical,        intent(out)               :: flagOK

    real, OPTIONAL, intent(out)               :: XS
    real, OPTIONAL, intent(out), dimension(nC):: XS_Arr

    integer, parameter :: nPart=5

    type(particle),dimension(1:1)      :: OutPart_Res, OutPart_QE
    type(particle),dimension(1:2)      :: OutPart_1pi, OutPart_VMDrho
    type(particle),dimension(1:2)      :: OutPart_2p2hQE, OutPart_2p2hDelta
    type(particle),dimension(1:3)      :: OutPart_2pi
    type(particle),dimension(1:nPart)  :: OutPart_DIS

    real, dimension(1:nC)          :: Weights
    real                           :: totalWeight
    logical :: flagGuess

    flagOK =.false.
    if (present(XS)) XS = 0.0
    if (present(XS_Arr)) XS_Arr = 0.0


    ! (0) Init the module

    if (initFlag) call initInput

    ! (1) Evaluate cross sections and generate for each channel a final state

    Weights = 0.0
    if (doC(iC_QE))  call init_QE(eN,OutPart_QE,weights(iC_QE))
    if (doC(iC_Res)) call init_Res(eN,whichRes,OutPart_Res,weights(iC_Res))
    if (doC(iC_1Pi)) call init_1Pi(eN,OutPart_1pi,weights(iC_1pi),-999.,-999., doC(iC_Res))
    if (doC(iC_2Pi)) call init_2Pi(eN,OutPart_2pi,weights(iC_2pi))
    if (doC(iC_DIS)) then
       call resetNumberGuess()
       call init_DIS(eN,OutPart_DIS,weights(iC_DIS))
       flagGuess = AcceptGuessedNumbers()
    end if
    if (doC(iC_VMDrho)) call init_VMDrho(eN,OutPart_VMDrho,weights(iC_VMDrho))
    if (doC(iC_2p2hQE)) call init_2p2hQE(eN,OutPart_2p2hQE,weights(iC_2p2hQE))
    if (doC(iC_2p2hDelta)) call init_2p2hDelta(eN,OutPart_2p2hDelta,weights(iC_2p2hDelta))

    ! (1a) Check Pauli:

    if (DoPauli) then
       if (abs(weights(iC_QE))>1e-5) then
          if (.not.checkPauli(OutPart_QE,realparticles)) weights(iC_QE)=0.0
       end if
       if (abs(weights(iC_Res))>1e-5) then
          if (.not.checkPauli(OutPart_Res,realparticles)) weights(iC_Res)=0.0
       end if
       if (abs(weights(iC_1Pi))>1e-5) then
          if (.not.checkPauli(OutPart_1Pi,realparticles)) weights(iC_1Pi)=0.0
       end if
       if (abs(weights(iC_2Pi))>1e-5) then
          if (.not.checkPauli(OutPart_2Pi,realparticles)) weights(iC_2Pi)=0.0
       end if
       if (abs(weights(iC_DIS))>1e-5) then
          if (.not.checkPauli(OutPart_DIS,realparticles)) weights(iC_DIS)=0.0
       end if
       if (abs(weights(iC_VMDrho))>1e-5) then
          if (.not.checkPauli(OutPart_VMDrho,realparticles)) weights(iC_VMDrho)=0.0
       end if
       if (abs(weights(iC_2p2hQE))>1e-5) then
          if (.not.checkPauli(OutPart_2p2hQE,realparticles)) weights(iC_2p2hQE)=0.0
       end if
       if (abs(weights(iC_2p2hDelta))>1e-5) then
          if (.not.checkPauli(OutPart_2p2hDelta,realparticles)) weights(iC_2p2hDelta)=0.0
       end if
    end if


    ! (2) Choose final state: Make decision for a special channel

    !write(*,'(A,1P,10e12.3)') 'weights:',weights

    OutPart%ID=0
    channel=0

    select case (monteCarloChoose(weights,totalWeight))
    case (0) ! ===== failure =====
       return ! all cross sections are zero

    case (1) ! ===== QE Event =====
       OutPart(1:1)=OutPart_QE
       channel=OutPart_QE(1)%ID

    case (2) ! ===== Resonance Event =====
       OutPart(1:1)=OutPart_Res
       channel=OutPart_Res(1)%ID

    case (3) ! ===== Pi N Event =====
       OutPart(1:2)=OutPart_1pi
       channel=origin_singlePi

    case (4) ! ===== Pi Pi N Event =====
       OutPart(1:3)=OutPart_2pi
       channel=origin_doublePi

    case (5) ! ===== DIS Event =====
       OutPart(1:nPart)=OutPart_DIS
       channel=origin_DIS

    case (6) ! ===== VMDrho Event =====
       OutPart(1:2)=OutPart_VMDrho
       channel=origin_vecmes_rho

    case (7) ! ===== 2p2h QE Event =====
       OutPart(1:2)=OutPart_2p2hQE
       channel=origin_2p2hQE

    case (8) ! ===== 2p2h Delta Event =====
       OutPart(1:2)=OutPart_2p2hDelta
       channel=origin_2p2hDelta

    case default
       call TRACEBACK('Wrong Monte-Carlo Decision')
    end select

    flagOK=.true.
    OutPart%perweight=totalWeight

    flagOK=checkEvent(eN,OutPart,channel)

    if (.not.flagOK) then
       OutPart%ID=0
       if (present(XS)) XS = 0.0
       if (present(XS_Arr)) XS_Arr = 0.0
    else
       if (present(XS)) XS = totalWeight
       if (present(XS_Arr)) XS_Arr = weights
    end if

  end subroutine eventGen_eN_lowEnergy

  !****************************************************************************
  !****f* eventGenerator_eN_lowEnergy/checkEvent
  ! NAME
  ! logical function checkEvent(eN,f,channel)
  !
  ! PURPOSE
  ! Check Charge and momentum conservation (.true. -> all is ok!)
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- incoming electron and nucleon
  ! * type(particle),dimension(:) :: f -- final state particles
  ! * integer                     :: channel -- how the event was generated,
  !   cf. iC_XXX and module Electron_origin
  ! OUTPUT
  ! * function value
  ! NOTES
  ! * if channel==2p2hQE|2p2hDelta, no checks are possible; alsways .true.
  !****************************************************************************
  logical function checkEvent(eN,f,channel)
    use particleDefinition
    use eN_eventDefinition, only: electronNucleon_event, write_electronNucleon_event
    use output
    use Electron_origin, only: origin_DIS,origin_2p2hQE,origin_2p2hDelta

    type(particle),dimension(:),intent(in) :: f
    type(electronNucleon_event) :: eN
    real, dimension(0:3)        :: mom_in,mom_out
    integer                     :: charge_in,charge_out,j
    integer, intent(in)         :: channel

    checkEvent=.true.

    ! if it was a 2p2h event, no checks are possible.
    if (channel.eq.origin_2p2hQE) return
    if (channel.eq.origin_2p2hDelta) return

    ! if it was a DIS event, we may violate energy-/momentum-conservation
    ! dramatically, since we may omit ('unknown') particles.
    if (channel.eq.origin_DIS) return

    ! Total Charge:
    ! =============
    charge_in =eN%nucleon%charge+eN%boson%charge
    charge_out=sum(f%charge, mask= f%ID.ne.0)
    if (charge_in.ne.charge_out) then
       write(*,*) 'Charge not conserved in eventGen_en_lowEnergy!!! ',channel
       write(*,*) 'IN:  ', charge_in
       write(*,*) 'OUT: ', charge_out
       call write_electronNucleon_event(eN)
       write(*,*)
       write(*,*) 'Final state particles at hadronic vertex:'
       do j=lbound(f,dim=1),ubound(f,dim=1)
          call WriteParticle_debug(f(j))
       end do
       checkEvent=.false.
       call TRACEBACK()
    end if

    ! Total Momentum:
    ! ===============

    ! do not use boson%momentum here, since this may be modified by the
    ! RemovePot routines
 !   write(*,*) 'nucleon-mom=',eN%nucleon%momentum
 !   write(*,*) 'lepton-mom_in=',eN%lepton_in%momentum
 !   write(*,*) 'lepton-mom_out=',eN%lepton_out%momentum
    mom_in =eN%nucleon%momentum+eN%lepton_in%momentum-eN%lepton_out%momentum
    do j=0,3
       mom_out(j)=sum(f%momentum(j), mask= f%ID.ne.0)
    end do
    if (maxVal(abs(mom_in-mom_out)).gt.1e-5) then
       write(*,*) 'Momentum not conserved in eventGen_en_lowEnergy!!! ',channel

       write(*,'(A,1P,5e13.4)') 'IN:  ', mom_in
       write(*,'(A,1P,5e13.4)') 'OUT: ', mom_out
       write(*,*)
       write(*,'(A,1P,5e13.4)') 'IN*: ', eN%nucleon%momentum+eN%boson%momentum
       write(*,*)

       call write_electronNucleon_event(eN)
       write(*,*)
       call WriteParticle(6, 1, f)
       write(*,*)

       write(*,*) 'Final state particles at hadronic vertex:'
       do j=lbound(f,dim=1),ubound(f,dim=1)
          if (f(j)%ID.ne.0) call WriteParticle_debug(f(j))
       end do
       checkEvent=.false.
       call TRACEBACK()
    end if

  end function checkEvent

  !****************************************************************************
  !****s* eventGenerator_eN_lowEnergy/init_2p2hQE
  ! NAME
  ! subroutine init_2p2hQE(eN,OutPart,XS)
  !
  ! PURPOSE
  ! Generate a e^- N1 N2 -> N1' N2' event.
  !
  ! NOTES
  ! just a wrapper for el_2p2h_DoQE
  !****************************************************************************
  subroutine init_2p2hQE(eN,OutPart,XS)
    use particleDefinition
    use eN_eventDefinition
    use lepton2p2h, only: lepton2p2h_DoQE

    type(electronNucleon_event)   , intent(inout) :: eN
    type(particle), dimension(:)  , intent(out)   :: OutPart
    real                          , intent(out)   :: XS

    type(electronNucleon_event) :: eN_dummy

    eN_dummy = eN
    call lepton2p2h_DoQE(eN,outPart,XS)
    call setNumber(OutPart)
    call setOutPartDefaults(OutPart, eN%nucleon)
    eN = eN_dummy
  end subroutine init_2p2hQE

  !****************************************************************************
  !****s* eventGenerator_eN_lowEnergy/init_2p2hDelta
  ! NAME
  ! subroutine init_2p2hDelta(eN,OutPart,XS)
  !
  ! PURPOSE
  ! Generate a e^- N1 N2 -> N Delta event.
  !
  ! NOTES
  ! just a wrapper for el_2p2h_DoDelta
  !****************************************************************************
  subroutine init_2p2hDelta(eN,OutPart,XS)
    use particleDefinition
    use eN_eventDefinition
    use lepton2p2h, only: lepton2p2h_DoDelta

    type(electronNucleon_event)   , intent(inout) :: eN
    type(particle), dimension(:)  , intent(out)   :: OutPart
    real                          , intent(out)   :: XS

    type(electronNucleon_event) :: eN_dummy

    eN_dummy = eN
    call lepton2p2h_DoDelta(eN,outPart,XS)
    call setNumber(OutPart)
    call setOutPartDefaults(OutPart, eN%nucleon)
    eN = eN_dummy

  end subroutine init_2p2hDelta

  !****************************************************************************
  !****s* eventGenerator_eN_lowEnergy/init_VMDrho
  ! NAME
  ! subroutine init_VMDrho(eN,OutPart,XS)
  !
  ! PURPOSE
  ! Generate a e^- N -> rho0 event.
  !
  ! NOTES
  ! This can be easily generalized to all other vector mesons.
  !****************************************************************************
  subroutine init_VMDrho(eN,OutPart,XS)
!    use offshellpotential    , only: setOffShellParameter
    use constants, only: pi, mN
    use particleDefinition
    use eN_eventDefinition
    use eN_event, only: eNeV_GetKinV!,eNeV_CheckForDIS
!    use PythiaSpecFunc       , only: Init_VM_Mass

    use photonXSections, only: calcXS_gammaN2VN
    use mediumDefinition
    use dichteDefinition
    use mediumModule, only: mediumAt,getMediumCutOff
    use RMF, only: getRMF_flag
    use IdTable, only: rho,nucleon
    use particleProperties, only: hadron
    use CollTools, only: PythiaMSTP
    use master_2Body, only: setKinematics
    use lorentzTrafo, only: lorentzCalcBeta !, lorentz
    use collisionNumbering, only: pert_numbering
    use propagation, only: updateVelocity
    use output
    use densitymodule, only: densityAt
    use ParamEP, only: CalcParamEP_R1990
    use PIL_rhoDiff, only: PIL_rhoDiffractive_PUT

    type(electronNucleon_event)   , intent(in)    :: eN
    type(particle), dimension(:)  , intent(out)   :: OutPart
    real                          , intent(out)   :: XS

    real :: nu,Q2,W,Wfree,eps,fT,epsR
    type(medium)         :: media
    real, dimension(1:4) :: sig
    real                 :: Dipole, LongEnh, mV
    real, parameter      :: PARP165 =0.5
    integer :: i
    real, dimension(1:3) :: betaToLRF
    real, dimension(1:3) :: betaToCM
    type(dichte) :: density
    type(particle), dimension (2) :: pairIN
    logical :: flagOK

    ! The following should not be a hard-wired parameter but be
    ! accessible from outside (jobcard??)
    integer, parameter :: XSscenario = 3


    flagOK = .false.
    XS = 0.0

    call eNeV_GetKinV(eN, nu,Q2,W,Wfree,eps,fT)

    ! instead of using the following parametrisation, which returns R
    ! according the TOTAL XS, we should better use something,
    ! which gives R_rho. (To Be Done!)

    call CalcParamEP_R1990(Wfree,Q2,epsR) ! returns only R
    epsR = epsR*eps

    !===== 1: Calculate cross section:

    media=mediumAt(eN%nucleon%position)

    select case (XSscenario)
    case (1) ! === PYTHIA-like ===

       call calcXS_gammaN2VN(Wfree,media,sig)

       mV = hadron(rho)%mass
       Dipole  = (mV**2/(mV**2+Q2))**2
       select case (PythiaMSTP(17))
       case (4)
          LongEnh = mV**2*Q2/((mV**2+Q2))**2
       case (5)
          LongEnh = Q2/(mV**2+Q2)
       case default
          write(*,*) 'MSTP(17)=',PythiaMSTP(17),'not implemented.'
          call TRACEBACK()
       end select
       LongEnh = 1. + eps * LongEnh * 4.*PARP165

       XS = sig(1) * Dipole * LongEnh

    case (2) ! === fit to CLAS data (total XS) ===

       ! Please note the Jacobian dnu/dW as the 'prefactor'
       XS = mN/W * (83.8*W**2-525*W+862)*(Q2+1.8)**(-3.27)

    case (3) ! fit to Pythia VMD contribution

       ! Please note the Jacobian dnu/dW as the 'prefactor'
       ! For W<2.7, there is indeed no W dependence !
       XS = mN/W * 89.8*(Q2+1.30)**(-4.33)

    end select

    fT = fT/ ( 1e3* pi/(eN%lepton_out%momentum(0)*eN%lepton_in%momentum(0)))
    XS = XS*fT ! prevent the XS from dividing by flux

    !===== 2: Generate final state:

    pairIN = (/ eN%boson, eN%nucleon /)

    call setToDefault(OutPart)

    OutPart%ID     = (/ rho, nucleon/)
    OutPart%Charge = (/ 0,   eN%nucleon%Charge /)

    OutPart%perturbative= .true.
    OutPart%perWeight   = XS ! perturbative weight

    do i=1,2
       OutPart(i)%position=eN%nucleon%position
       OutPart(i)%event=pert_numbering(eN%nucleon)
    end do

    density=densityAt(eN%nucleon%position)
    if (density%baryon(0).gt.getMediumCutOff()/100. .and. .not.getRMF_flag() ) then
      betaToLRF = lorentzCalcBeta (density%baryon, 'init_VMDrho')
    else
      betaToLRF=0.
    end if
    betaToCM = lorentzCalcBeta (pairIN(1)%momentum + pairIN(2)%momentum)

    call setKinematics(W,Wfree,betaToLRF,betaToCM,media,pairIN, OutPart,flagOK)

    if (.not.flagOK) then
       XS = 0.0
       return
    end if


    call updateVelocity(OutPart)
    call setNumber(OutPart)

!    call WriteParticle(6,1,pairIN)
!    call WriteParticle(6,2,OutPart)

    call PIL_rhoDiffractive_PUT(outPart(1)%number,&
         & .true.,epsR,outPart(2)%momentum)

!    call write_electronNucleon_event(eN)
!    call WriteParticle(6,1,OutPart)

  end subroutine init_VMDrho

  !****************************************************************************
  !****s* eventGenerator_eN_lowEnergy/init_DIS
  ! NAME
  ! subroutine init_DIS(eN,OutPart,XS)
  !
  ! PURPOSE
  ! Generate a e^- N -> DIS event.
  !
  ! NOTES
  ! This is mainly just a wrapper around DoColl_gammaN_Py.
  !
  ! The returned cross section (and also the weights of the particles)
  ! is 'flux * sigma^*',
  !****************************************************************************
  subroutine init_DIS(eN,OutPart,XS)
    use offshellpotential, only: setOffShellParameter
    use constants, only: pi
    use particleDefinition
    use eN_eventDefinition
    use eN_event, only: eNeV_GetKinV,eNeV_CheckForDIS
    use PythiaSpecFunc, only: Init_VM_Mass
    use Coll_gammaN, only: DoColl_gammaN_Py

    type(electronNucleon_event),  intent(in)  :: eN
    type(particle), dimension(:), intent(out) :: OutPart
    real,                         intent(out) :: XS

    real :: nu,Q2,W,Wfree,eps,fT
    real, dimension(0:4) :: Cross
    integer              :: EventClass
    logical              :: flagOK

    real, dimension(1:4), parameter :: rVMD = 1.0
    logical,              parameter :: DoDifr = .false.

    ! some adhoc cuts on W:
    ! do not go below 1.5, because then troubles with PYTHIA parameters
!    real, parameter :: Wcut1 = 1.60, Wcut2 = 1.65
    real, parameter :: Wcut1 = 1.50, Wcut2 = 1.50

    flagOK = .false.
    XS = 0.0

    call eNeV_GetKinV(eN, nu,Q2,W,Wfree,eps,fT)

    if (Wfree.le.Wcut1) return
    if (.not.eNeV_CheckForDIS(eN)) return ! avoid infinite loop

    call Init_VM_Mass(Wfree,eN%nucleon%position)
    call DoColl_gammaN_Py(eN,OutPart,flagOK, rVMD, DoDifr, Cross,EventClass,minW=1.0)

    fT = fT/ ( 1e3* pi/(eN%lepton_out%momentum(0)*eN%lepton_in%momentum(0)))
    if (flagOK) XS = 1000*Cross(0)
    if (Wfree.le.Wcut2) XS = XS*(Wfree-Wcut1)/(Wcut2-Wcut1) ! some adhoc W dependence

    if (flagOK) call setOffShellParameter(OutPart, flagOk)
    if (.not.flagOK) XS = 0.0

    XS = XS*fT ! prevent the XS from dividing by flux
    OutPart%perWeight=XS

    !===== set some additional fields:

    call setNumber(OutPart) ! ?????
    call setOutPartDefaults(OutPart, eN%nucleon)

  end subroutine init_DIS

  !****************************************************************************
  !****s* eventGenerator_eN_lowEnergy/init_2Pi
  ! NAME
  ! subroutine init_2Pi(eN,OutPart,XS)
  !
  ! PURPOSE
  ! Generate a e^- N -> e^- N pion pion event.
  !
  ! NOTES
  ! This is done by taking the (vacuum) 2pi-background in photoproduction
  ! and scaling it with Q2 as the total XS. (many better scalings possible!)
  !
  ! Restricted to Q^2 < 5 GeV^2 and W < 2 GeV, otherwise set to zero.
  !****************************************************************************
  subroutine init_2Pi(eN,OutPart,XS, avoidSetNumbers)
    use offShellPotential, only: setOffShellParameter
    use constants, only: pi, mN, mPi
    use particleDefinition
    use mediumDefinition
    use IdTable, only: pion,nucleon
    use eN_eventDefinition
    use eN_event, only: eNeV_GetKinV
    use random, only: rn
    use nBodyPhaseSpace, only: momenta_in_3BodyPS
    use ParamEP, only: CalcParamEP
    use energyCalc, only: energyCorrection
    use mediumModule, only: mediumAt

    type(electronNucleon_event) , intent(in)  :: eN
    type(particle), dimension(3), intent(out) :: OutPart
    real                        , intent(out) :: XS
    logical, intent(in), optional             :: avoidSetNumbers

    integer :: i
    integer :: qnuk
    real :: nu,Q2,W,Wfree,eps,fT
    real, dimension(0:3) :: ptot
    real, dimension(1:3) :: betaCMToLab = 0.0
!     type(medium)         :: mediumDUMMY    ! use as DUMMY
    type(medium)         :: mediumAtPosition

    real, dimension(0:3) :: sig2Pi !,sigRes_2pi

    real :: randomNumber
    real , dimension(1:3,1:3) :: p3       ! momenta of three particles
    real :: dum1,dum2
    logical :: flagOK, doSetNumbers

    flagOK = .false.
    XS = 0.0
    OutPart%ID=0
    OutPart%perweight=0.0

    if (eN%QSquared.ge.5.0) return

    call eNeV_GetKinV(eN, nu,Q2,W,Wfree,eps,fT)
    qnuk = eN%nucleon_free%charge

    !===== 1: Calculate cross section gamma N -> N pi pi (Q2=0):
    !===== 2: Subtract resonance cross section gamma N -> R -> N pi pi (Q2=0):
    call init_2Pi_getBG(eN%nucleon_free, Wfree, sig2pi)

    if (sig2Pi(0).lt.10e-20) return

    !===== 3: generate an event (Q2=0):

    call setToDefault(OutPart)

    OutPart(1:3)%ID=(/nucleon, pion, pion /)
    OutPart(1)%mass=mN
    OutPart(2:3)%mass=mPi

    OutPart(1:3)%antiparticle=.false.
    OutPart(1:3)%scaleCS=1.

    randomNumber=rn()*sig2Pi(0)

    if (randomNumber .lt. sig2pi(2)) then ! === one charged pion + pi^0
       OutPart(1:3)%charge = (/abs(qnuk-1), 0, qnuk-abs(qnuk-1)/)

    else if (randomNumber .lt. sig2pi(2)+sig2pi(3)) then ! === uncharged
       OutPart(1:3)%charge = (/qnuk, 0, 0/)

    else ! === double charged channel: gamma N -> N pi^+ pi^-
       OutPart(1:3)%charge = (/qnuk, 1, -1/)

    end if

    p3 = momenta_in_3BodyPS (Wfree, OutPart(1:3)%mass)

    do i=1,3
       OutPart(i)%momentum(1:3)=p3(:,i)
       OutPart(i)%momentum(0)=FreeEnergy(OutPart(i))
    end do

    !===== 4: boost the event to the gamma* N system:

    mediumAtPosition=mediumAt(eN%nucleon%position)

    ptot = eN%nucleon%momentum+eN%boson%momentum
    betaCMToLab = ptot(1:3)/ptot(0)

    call energyCorrection(W, (/0.,0.,0./), betaCMToLab, mediumAtPosition,OutPart,flagOK)

    !write(*,*) sqrts(OutPart), W
    if (.not.flagOK) then
       OutPart(1:3)%ID = 0
       return
    end if

    !===== 5: scale the cross section according Q2 dependence of sigma_tot:

    XS = sig2Pi(0)

    call CalcParamEP(Wfree,0.0,0.0, dum1)
    call CalcParamEP(Wfree,Q2, eps, dum2)
    XS = XS * dum2/dum1

    !===== 6: Take care of offshellness
    call setOffShellParameter(OutPart,flagOK)
    if (.not.flagOK) XS = 0.0

    !===== 7:  prevent the 2Pi-BG XS from dividing by flux

    XS = XS*fT/ ( 1e3* pi/(eN%lepton_out%momentum(0)*eN%lepton_in%momentum(0)))
    OutPart(1:3)%perWeight=XS

    !===== 8:  set some additional fields:

    doSetNumbers = .true.
    if (present(avoidSetNumbers)) doSetNumbers=.not.avoidSetnumbers
    if (doSetNumbers) call setNumber(OutPart) ! ?????
    call setOutPartDefaults(OutPart, eN%nucleon)

  end subroutine init_2Pi

  !****************************************************************************
  !****s* eventGenerator_eN_lowEnergy/init_2Pi_getBG
  ! NAME
  ! subroutine init_2Pi_getBG(nucleon_free, Wfree, sig2pi)
  !
  ! PURPOSE
  ! Calculate the background contribution for an e^- N -> e^- N pion pion
  ! event at Q^2=0.
  !
  ! Since te resonance cross section is only calculated up to W=2 GeV,
  ! it is also only possible to get some 2pi BG up to this value.
  !
  ! NOTES
  ! This is a helper routine for init_2Pi, but also used in the neutrino
  ! case
  !****************************************************************************
  subroutine init_2Pi_getBG(nucleon_free, Wfree, sig2pi)
    use constants, only: mN
    use particleDefinition
    use mediumDefinition
    use gamma2Pi_Xsections, only: gamma2pi
    use resProd_lepton, only: sigma_pipi_res_vac

    type(particle), intent(in)        :: nucleon_free
    real, intent(in) :: Wfree
    real, dimension(0:3), intent(out)  :: sig2Pi

    integer :: i
    real, dimension(0:3) :: sigRes_2pi
    real :: nu
    real, dimension(0:3) :: gammaMomentum, ptot
    real, dimension(1:3) :: betaCMToLab = 0.0
    type(medium)         :: mediumDUMMY    ! use as DUMMY
    real, dimension(1:3),parameter :: posDUMMY = -1000.0 ! somewhere in the no man's land, use as DUMMY

    !===== 0: check kinematics

    if (Wfree > 2.0) then
       sig2Pi = 0.0
       return
    end if

    !===== 1: Calculate cross section gamma N -> N pi pi (Q2=0):

    nu = (Wfree**2-mN**2)/(2*mN)
    gammaMomentum = (/nu,0.0,0.0,nu/)
    ptot = gammaMomentum + nucleon_free%momentum
    betaCMToLab = ptot(1:3)/ptot(0)

    call gamma2pi(nucleon_free%charge,Wfree,sig2pi,betaCMToLab,mediumDUMMY,posDUMMY)

    !===== 2: Subtract resonance cross section gamma N -> R -> N pi pi (Q2=0):

    ! since we use the free nucleon, we use the routine xxx_vac
    sigRes_2Pi=sigma_pipi_res_vac(nucleon_free,gammaMomentum)*1000.

    do i=0,3
       sig2Pi(i) = max(0.0, sig2Pi(i)-sigRes_2Pi(i))
    end do

  end subroutine init_2Pi_getBG

  !****************************************************************************
  !****s* eventGenerator_eN_lowEnergy/init_1Pi
  ! NAME
  ! subroutine init_1Pi(eN,OutPart,XS,theta_k,phi_k,modeBckGrnd)
  !
  ! PURPOSE
  ! Generate one e^- N -> e^- N pion event.
  !
  ! The weight of the each event is given by (total Xsection of event=sum
  ! over all possible pion charges).
  ! We make a Monte Carlo decision  to determine the pion charge.
  !
  ! The particles are produced at the place of the nucleon target.
  !
  ! If one or more of the input angles theta_k, phi_k are negative, then we
  ! make a Monte-Carlo decision for those angles which are negative. In this
  ! procedure, we distribute:
  ! * phi_k flat in [0,2*pi]
  ! * cos(theta_k) flat in [-1,1]
  ! Note that those angles are defined in the CM-frame of the outgoing pion and nucleon.
  !
  ! If a Monte-Carlo decision on the angle is performed, then the perweight of each event
  ! includes the following integral measure:
  ! * phi and theta decision: measure=int dphi_k dtheta_k=4*pi
  !   --> perweight=4*pi*dsigma_dOmega_pion(phi_k,theta_k)
  ! * phi integration       : measure=int dphi_k dtheta_k=2*pi
  !   --> perweight=2*pi*dsigma_dOmega_pion(phi_k,theta_k
  ! * theta integration     : measure=int dtheta_k       =2
  !   --> perweight=2   *dsigma_dOmega_pion(phi_k,theta_k
  !
  ! INPUTS
  ! * type(electronNucleon_event)  :: eN
  !   -- The underlying electron and nucleon event
  ! * real                         :: theta_k,phi_k
  !   -- Outgoing pion angles (if negative then Monte Carlo Integration).
  !   The angles are defined in the CM-Frame of the outgoing pion and nucleon
  ! * logical                      :: modeBckGrnd
  !   -- true: mode==background, false: mode=normal
  ! OUPTUT
  ! * type(particle), dimension(1:2)   :: OutPart
  ! * logical                          :: flagOK
  !****************************************************************************
  subroutine init_1Pi(eN,OutPart,XS,theta_k,phi_k,modeBckGrnd)

    use offShellPotential, only: setOffShellParameter
    use resProd_lepton, only: dSdO_fdE_fdO_k_med_res_EN
    use electronPionProd_medium_eN, only: getKinematics_eN, dSdO_fdE_fdO_k_med_eN
    use particleDefinition
    use random, only: rn, rnCos
    use constants, only: pi
    use degRad_conversion, only: degrees
    use eN_eventDefinition, only: electronNucleon_event,setVacuum
    use leptonicID, only: EM
    use monteCarlo, only: monteCarloChoose

    type(electronNucleon_event) , intent(in)  :: eN
    real                        , intent(in)  :: theta_k,phi_k
    logical                     , intent(in)  :: modeBckGrnd


    type(particle), dimension(2), intent(out) :: OutPart
    real                        , intent(out) :: XS

    ! local variables
    integer                     :: pionCharge
    real                        :: theta_k_MC, deltaMonteCarlo, phi_k_MC,totalWeight
    real, dimension(-1:1)       :: sigma, sigmaRes
    type(electronNucleon_event) :: eN_vacuum
    ! The underlying electron and nucleon event transformed to vacuum kinematics
    real, dimension(0:3,-1:1)   :: k,pf
    logical                     :: twoRoots
    integer, parameter          :: CM=2
    logical                     :: flagOK

    flagOK = .false.
    XS = 0.0
    OutPart%ID=0
    OutPart%perweight=0.0
    if (eN%QSquared.ge.5.0) return

   ! Check for Monte-Carlo-Integrations:
    deltaMonteCarlo=1.
    ! * Pion theta angle
    if (theta_k.lt.0) then
       theta_k_MC=degrees(rnCos())
       deltaMonteCarlo=deltaMonteCarlo*2.
    else
       theta_k_MC=theta_k
    end if
    ! * Pion phi angle
    if (phi_k.lt.0) then
       phi_k_MC=degrees(rn()*2*pi)
       deltaMonteCarlo=deltaMonteCarlo*2.*pi
    else
       phi_k_MC=phi_k
    end if

    ! Evaluate the cross section for all possible pion charges:
    sigma=0.

    if (modeBckGrnd) eN_vacuum=setVacuum(eN)

    pionChargeLoop: do pionCharge=eN%nucleon%charge-1,eN%nucleon%charge

       if (.not.modeBckGrnd) then

          sigma(pionCharge)= dSdO_fdE_fdO_k_med_eN(eN, pionCharge, &
               & phi_k_MC, theta_k_MC, &
               & k(:,pionCharge), pf(:,pionCharge), EM,pionNucleonSystem=CM)

       else

          sigma(pionCharge)= dSdO_fdE_fdO_k_med_eN(eN_vacuum, pionCharge, &
               & phi_k_MC, theta_k_MC, &
               & k(:,pionCharge), pf(:,pionCharge), EM,pionNucleonSystem=CM)

          if (sigma(pionCharge).lt.1E-20) then
             ! no solution to kinematical problem was found
             cycle pionChargeLoop
          end if

          ! Resonance contribution:
          sigmaRes=dSdO_fdE_fdO_k_med_res_eN(eN_vacuum,&
               & k(:,pionCharge),pf(:,pionCharge),EM,pionNucleonSystem=CM)
          ! Subtract resonance contributions:
          sigma(pionCharge)=sigma(pionCharge)-sigmares(pionCharge)

          ! Get In-Medium Kinematics for given pion angles
          ! Note: We perform this for each pion charge seperately since the
          ! pion momentum depends on the pion charge,
          ! such that the pion momentum should also depend on the pion charge
          call getKinematics_eN(eN,pionCharge,eN%nucleon%charge-pionCharge,&
               & phi_k_MC,theta_k_MC,k(:,pionCharge),pf(:,pionCharge), &
               & twoRoots,flagOK,pionNucleonSystem=CM)
          if (twoRoots) then
             call TRACEBACK('ERROR in init_1Pi. Two roots in getKinematics_eN!!!')
          end if

          if (.not.flagOK) then
             ! Kinematics could not be established in the medium
             sigma(pionCharge)=0.
          end if

       end if

    end do pionChargeLoop

    ! Check that cross section is not zero:
    if (sum(abs(sigma)).lt.10e-20) return

    sigma=sigma*deltaMonteCarlo

    ! * Choose pion charge
    pionCharge=monteCarloChoose(sigma,totalWeight)-2

    call generateEvent_1Pi(OutPart,en%nucleon, &
         & k(:,pionCharge),pf(:,pionCharge), totalWeight,pionCharge)
    flagOK =.TRUE.

    ! Take care off offshellness
    call setOffShellParameter(OutPart,flagOK)
    if (.not.flagOK) OutPart%perweight=0.0
    XS = OutPart(1)%perweight

  end subroutine init_1Pi

  !****************************************************************************
  !****s* eventGenerator_eN_lowEnergy/init_QE
  ! NAME
  ! subroutine init_QE(eN,OutPart,XS)
  !
  ! PURPOSE
  ! Initializes one e^- N -> e^- N' events. The weight of the event
  ! is given by the  Xsection for QE scattering.
  !
  ! The particles are produced at the place of the nucleon target.
  !
  ! INPUTS
  ! * type(electronNucleon_event)  :: eN -- The underlying electron and nucleon event
  ! OUPTUT
  ! * real                         :: XS - The resulting cross section (=0 for failure)
  ! * type(particle)               :: OutPart
  !****************************************************************************
  subroutine init_QE(eN,OutPart,XS)
    use offShellPotential, only: setOffShellParameter
    use quasiElastic_electron, only: dSigmadcosTheta_l_dE_l_BW_eN
    use eN_eventDefinition, only: electronNucleon_event
    use constants, only: pi
    use particleDefinition, only: particle
    use idTable, only: nucleon
    use particleProperties, only: hadron

    type(electronNucleon_event) , intent(in)  :: eN
    real                        , intent(out) :: XS
    type(particle),dimension(1:1),intent(out) :: OutPart

    real                  :: sigma, baremass_out
    real, dimension (0:3) :: pf
    logical               :: flagOK

    XS = 0.0
    flagOK = .false.
    OutPart%ID=0
    OutPart%perweight=0.0


    ! Evaluate dSigma/dcos(Theta_l)/dE_l and convert it to dSigma/dOmega_l/dE_l
    sigma=dSigmadcosTheta_l_dE_l_BW_eN(eN,pf,baremass_out)/(2*pi)

    ! Check that cross section is not zero:
    if ((sigma.lt.10E-20).or.(baremass_out.lt.min(0.5,hadron(nucleon)%minmass)))  return

    ! Generate Event
    call generateEvent_1Body(OutPart(1),eN%nucleon,pf(:),sigma,nucleon,baremass_out)
    flagOK=.true.

    ! Take care of offshellness
    call setOffShellParameter(OutPart(1),flagOK)
    if (.not.flagOK) OutPart(1)%perweight=0.0
    XS = OutPart(1)%perweight

  end subroutine init_QE


  !****************************************************************************
  !****s* eventGenerator_eN_lowEnergy/init_Res
  ! NAME
  ! subroutine init_Res(eN,OutPart,XS)
  !
  ! PURPOSE
  ! Initializes an e^- N -> Resonance event.
  ! The (per)weight of the event is given by
  ! (total Xsection of event=sum over all resonances).
  ! We make a Monte Carlo decision to determine the resonance type.
  !
  ! The particles are produced at the place of the nucleon target.
  !
  ! INPUTS
  ! * type(electronNucleon_event)  :: eN
  !   -- The underlying electron and nucleon event
  ! * logical, dimension(2:nres+1) :: useRes
  !   -- Switch on/off each resonance
  ! OUPTUT
  ! * real                         :: XS
  !   -- The resulting cross section (=0 for failure)
  ! * type(particle)               :: OutPart
  !****************************************************************************
  subroutine init_Res(eN,useRes,OutPart,XS)
    use offShellPotential, only: setOffShellParameter
    use resProd_lepton, only: dSigmadOmega_fdE_f_resProd_eN
    use eN_eventDefinition, only: electronNucleon_event
    use IdTable, only: nres, nucleon
    use particleDefinition, only: particle
    use monteCarlo, only: monteCarloChoose

    type(electronNucleon_event) , intent(in)  :: eN
    logical, dimension(2:nres+1), intent(in)  :: useRes
    real                        , intent(out) :: XS
    type(particle),dimension(1:1),intent(out) :: OutPart

    real                                      :: totalWeight
    real, dimension(nucleon+1:nucleon+nres)   :: sigma, bareMass
    real, dimension(0:3,nucleon+1:nucleon+nres)  :: pf
    integer                                   :: resID

    logical                                   :: flagOK

    flagOK = .false.
    XS = 0.0
    OutPart%ID=0
    OutPart%perweight=0.0
    if (eN%QSquared.ge.5.0) return

    sigma=0.
    ! Evaluate Xsection for each resonance
    resIDLoop: do resID=nucleon+1,nucleon+nres
       if (useRes(resID)) then
          sigma(resID)=dSigmadOmega_fdE_f_resProd_eN(eN, resID,&
               & pf(:,resID),baremass(resID))
       end if
    end do resIDLoop


    ! Check that cross section is not zero:
    if (sum(sigma).lt.10E-20)  return

    ! Generate Event
    resID=monteCarloChoose(sigma,totalWeight)+1

    call generateEvent_1Body(OutPart(1),eN%nucleon,&
         & pf(:,resID),totalWeight,resID,baremass(resID))
    flagOK=.true.

    ! Take care off offshellness
    call setOffShellParameter(OutPart(1),flagOK)
    if (.not.flagOK) OutPart(1)%perweight=0.0
    XS = OutPart(1)%perweight

  end subroutine init_Res


  !****************************************************************************
  !****s* eventGenerator_eN_lowEnergy/generateEvent_1Pi
  ! NAME
  ! subroutine generateEvent_1Pi(OutPart,initNuc,kf,pf,xSection,pionCharge)
  !
  ! PURPOSE
  ! Given the evaluated kinematics and cross sections, final state particles
  ! are initialized for pion nucleon production.
  !
  ! INPUTS
  ! * integer                       :: pionCharge -- Charge of outgoing pion
  ! * type(particle)                :: initNuc    -- initial Nucleon
  ! * real, dimension(0:3)          :: kf,pf      -- pion momentum and nucleon momentum
  ! * real                          :: xSection
  !   -- Xsection for producing this event, including e.g. d(Omega)
  ! OUPTUT
  ! * type(particle),dimension(1:2) :: OutPart -- final state particles
  !****************************************************************************
  subroutine generateEvent_1Pi(OutPart,initNuc,kf,pf,xSection,pionCharge)
    use IdTable
    use propagation, only: updateVelocity
    use particleDefinition
    use constants, only: mN, mPi

    type(particle),dimension(1:2), intent(out) :: OutPart
    integer                      , intent(in)  :: pionCharge
    type(particle)               , intent(in)  :: initNuc
    real, dimension(0:3)         , intent(in)  :: kf,pf
    real                         , intent(in)  :: xSection

    call setToDefault(OutPart)
    OutPart%ID =     (/pion,       nucleon /)
    OutPart%Charge = (/pionCharge, initNuc%charge-pionCharge/)
    OutPart%Mass =   (/mPi,        mN/)

    OutPart(1)%momentum = kf
    OutPart(2)%momentum = pf

    OutPart%perWeight    = xSection ! perturbative weight

    call updateVelocity(OutPart)
    call setNumber(OutPart) ! ???????
    call setOutPartDefaults(OutPart, initNuc)

  end subroutine generateEvent_1Pi


  !****************************************************************************
  !****s* eventGenerator_eN_lowEnergy/generateEvent_1Body
  ! NAME
  ! subroutine  generateEvent_1Body(OutPart,initNuc,pf,xSection,ID,mass)
  !
  ! PURPOSE
  ! Given the evaluated kinematics and cross sections, final state particles
  ! are initialized for 1-body final states.
  !
  ! INPUTS
  ! * integer               :: ID         -- ID of produced particle
  ! * real                  :: mass       -- mass of produced particle
  ! * type(particle)        :: initNuc    -- initial Nucleon
  ! * real, dimension(0:3)  :: pf         -- final state
  ! * real                  :: xSection
  !   -- Xsection for producing this event, including e.g. d(Omega)
  ! OUPTUT
  ! * type(particle), intent(out) :: OutPart -- finalstate particle
  !****************************************************************************
  subroutine generateEvent_1Body(OutPart,initNuc,pf,xSection,ID,mass)
    use propagation, only: updateVelocity
    use particleDefinition

    type(particle),dimension(1:1), intent(out) :: OutPart
    type(particle)               , intent(in)  :: initNuc
    real, dimension(0:3)         , intent(in)  :: pf
    real                         , intent(in)  :: xSection
    integer                      , intent(in)  :: ID
    real                         , intent(in)  :: mass

    call setToDefault(OutPart)

    OutPart(:)%ID      = ID
    OutPart(:)%charge  = initNuc%charge
    OutPart(1)%momentum= pf
    OutPart(:)%mass    = mass
    OutPart(:)%perWeight   = xSection ! perturbative weight

    call updateVelocity(OutPart)
    call setNumber(OutPart) ! ????
    call setOutPartDefaults(OutPart, initNuc)

  end subroutine generateEvent_1Body

  !****************************************************************************
  !****is* eventGenerator_eN_lowEnergy/setOutPartDefaults
  ! NAME
  ! subroutine setOutPartDefaults(OutPart, initNuc)
  ! PURPOSE
  ! set some fields of the outgoing particles to default values
  !****************************************************************************
  subroutine setOutPartDefaults(OutPart, initNuc)
    use particleDefinition
    use collisionNumbering, only: pert_numbering

    type(particle), dimension(:), intent(inout) :: OutPart
    type(particle), intent(in) :: initNuc

    integer :: i, number

    number = pert_numbering(initNuc)

    do i=1,size(OutPart)
       OutPart(i)%position=initNuc%position
       OutPart(i)%event=number   ! Important for collision term: Particles won't collide immediately with InitNuc.
       OutPart(i)%perturbative = .true.
    end do

  end subroutine setOutPartDefaults

end module eventGenerator_eN_lowEnergy
