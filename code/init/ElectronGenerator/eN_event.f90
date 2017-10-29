!******************************************************************************
!****m* /eN_event
! NAME
! module eN_event
! PURPOSE
! Handle the data type which is used to store the lepton-nucleon event.
!******************************************************************************
module eN_event

  use particleDefinition
  use eN_eventDefinition
  use constants, only: mN

  IMPLICIT NONE

  private

  public :: eNev_SetProcess ! (idProcess,idFamily)  *** 1 ***

  public :: eNev_init_sWQ   ! (sqrts,W,  Q2,flagOK) *** 2 ***
  public :: eNev_init_sxQ   ! (sqrts,xBj,Q2,flagOK) *** 2 ***
  public :: eNev_init_snQ   ! (sqrts,nu, Q2,flagOK) *** 2 ***
  public :: eNev_init_eWQ   ! (eps,  W,  Q2,flagOK) *** 2 ***
  public :: eNev_init_exQ   ! (eps,  xBj,Q2,flagOK) *** 2 ***
  public :: eNev_init_enQ   ! (eps,  nu, Q2,flagOK) *** 2 ***
  public :: eNev_init_BWQ   ! (Ebeam,W,  Q2,flagOK) *** 2 ***
  public :: eNev_init_BnQ   ! (Ebeam,nu, Q2,flagOK) *** 2 ***

  public :: eNev_init_Target! (e0,pTarget,flagOK)   *** 3 ***

  public :: eNev_init_nuStep1  ! (pTarget)
  public :: eNev_init_nuStep2  ! (Enu,IP,flagOK)
  public :: eNev_init_nuStep3a ! (Eprime,costheta)
  public :: eNev_init_nuStep3b ! (Eprime,Q2)
  public :: eNev_init_nuStep3c ! (x,y)

  public :: init_electronNucleon_event ! see below

  public :: eNev_GetKinV ! (eps,nu,Q2,W)
  public :: eNeV_Get_LightX
  public :: eNeV_Get_CostLepton

  public :: eNeV_Set_PhiLepton
  public :: eNeV_CheckForDIS

  public :: eNeV_GetLeptonCM

  public :: nuclearFluxFactor_correction

!  real,        save :: s, sqrts      ! beam-target - cm energy
!  real,        save :: W2, W         ! photon-target -cm energy
!  real,        save :: x             ! energy fraction gamma/e (<>Bjorken)
!  real,        save :: Q2            ! photon virtuality
!  real,        save :: xB            ! Bjorken-x
!  real,        save :: yB            ! Bjorken-y, light cone fraction
!  real,        save :: nu            ! Bjorken, energy of photon
!  real,        save :: eps           ! epsilon = flux_L/flux_T
!  real,        save :: fT            ! flux_T [GeV^-2]

  real,        save :: M(2)  = (/ 0.00051, mN /)    ! masses of beam and target
  real,        save :: M2(2) = (/ 0.00051, mN /)**2 ! masses of beam and target (squared)

  logical,  save :: initFlag = .true.


  !****************************************************************************
  !****g* eN_event/selectFrame
  ! SOURCE
  !
  integer,  save :: selectFrame = 2
  !
  ! PURPOSE
  ! select frame, in which the calculaton of W_free is done:
  ! * 0 = doNOT   ---  do NOT remove
  ! * 1 = CM
  ! * 2 = CALC
  ! * 3 = THRE  prescription from correct threshold behaviour, used in heavy ion collisions
  ! * 4 = NucleonRest :  boost nucleon in the rest frame, set free mN, recalculate boson momentum
  ! * 5 = THRE2 threshold with m^2: sfree=s+m^2-m*^2
  !****************************************************************************

  !****************************************************************************
  !****g* eN_event/restingNucleon
  ! SOURCE
  !
  logical, save :: restingNucleon = .true.
  !
  ! PURPOSE
  ! if this flag is .false., we use the momentum of the target nucleon in the
  ! calculation of the flux
  !****************************************************************************

contains
  !****************************************************************************
  !****is* eN_event/readInput
  ! NAME
  ! subroutine readInput
  !****************************************************************************
  subroutine readInput
    use output

    integer :: ios

    !**************************************************************************
    !****n* eN_event/eN_Event
    ! NAME
    ! NAMELIST /eN_Event/
    ! PURPOSE
    ! Namelist for stuff connected with the storage of the lepton-nucleon event:
    ! * selectFrame
    ! * restingNucleon
    !**************************************************************************
    NAMELIST /eN_event/ selectFrame,restingNucleon

    if (.not.initFlag) return

    call Write_ReadingInput('eN_event',0)
    rewind(5)
    read(5,nml=eN_event,IOSTAT=ios)
    call Write_ReadingInput("eN_event",0,ios)

    ! block non-default values:
    if (selectFrame.ne.2) then
       call notInRelease("selectFrame != 2")
    end if

    select case (selectFrame)
    case (0)
       write(*,*) 'selectFrame = ',selectFrame,'  = doNOT'
    case (1)
       write(*,*) 'selectFrame = ',selectFrame,'  = CM'
    case (2)
       write(*,*) 'selectFrame = ',selectFrame,'  = CALC'
    case (3)
       write(*,*) 'selectFrame = ',selectFrame,'  = THRE'
!    case (4)
!       write(*,*) 'selectFrame = ',selectFrame,'  = NucleonRest'
    case (5)
       write(*,*) 'selectFrame = ',selectFrame,'  = THRE2'
    case default
       write(*,*) 'selectFrame = ',selectFrame,'  = ***unknown***, stop!'
       stop
    end select
    write(*,*) 'restingNucleon = ',restingNucleon

    call Write_ReadingInput('eN_event',1)

    initFlag=.false.

  end subroutine readInput

  !****************************************************************************
  !****s* eN_event/eNev_SetProcess
  ! NAME
  ! subroutine eNev_SetProcess(e, idProcess,idFamily)
  !
  ! PURPOSE
  ! Set the process and the lepton family for the initialisation
  !
  ! INPUTS
  ! * idProcess: +-1=(anti-)EM, +-2=(anti-)CC, +-3=(anti-)NC
  ! * idFamily: 1=e-, 2=mu-, 3=tau-
  !****************************************************************************
  subroutine eNev_SetProcess(e, idProcess,idFamily)
    use IdTable
    use CallStack

    type(electronNucleon_event), intent(inout) :: e
    integer, intent(in) :: idProcess,idFamily

    type(electronNucleon_event),save :: e0 ! standard initialized eventtype

    integer, dimension(3), parameter :: arrID1 = (/electron,muon,tau/)
    integer, dimension(3), parameter :: arrID2 = (/electronNeutrino,muonNeutrino,tauNeutrino/)
    real, dimension(3), parameter :: arrMass = (/0.00051,0.10566,1.77700/)

    if (initFlag) call readInput

    if ((idFamily.lt.1).or.(idFamily.gt.3)) then
       write(*,*) 'idFamily=',idFamily
       call TRACEBACK('wrong idFamily!')
    end if

    e = e0 ! reset the whole eventtype

    select case (idProcess)
    case (-1,1) !=== EM

       call setToDefault(e%lepton_in)

       e%lepton_in%ID = arrID1(idFamily)
       e%lepton_in%mass = arrMass(idFamily)

       if (idProcess.gt.0) then
          e%lepton_in%charge = -1
       else
          e%lepton_in%charge = 1
       end if
       e%lepton_in%momentum = (/e%lepton_in%mass, 0.0, 0.0, 0.0/)

       e%lepton_out = e%lepton_in

       call setToDefault(e%boson)
       e%boson%ID = photon

    case (-2,2) !=== CC

       call setToDefault(e%lepton_in)
       e%lepton_in%ID = arrID2(idFamily)
       e%lepton_in%mass = 0.0
       e%lepton_in%charge = 0
       e%lepton_in%antiparticle = (idProcess.lt.0)
       e%lepton_in%momentum = (/0.0, 0.0, 0.0, 0.0/)

       call setToDefault(e%lepton_out)
       e%lepton_out%ID = arrID1(idFamily)
       e%lepton_out%mass = arrMass(idFamily)
       if (idProcess.gt.0) then
          e%lepton_out%charge = -1
       else
          e%lepton_out%charge = 1
       end if
       e%lepton_out%momentum = (/e%lepton_out%mass, 0.0, 0.0, 0.0/)

       call setToDefault(e%boson)
       e%boson%ID = Wboson
       e%boson%charge = -e%lepton_out%charge

    case (-3,3) !=== NC

       call setToDefault(e%lepton_in)
       e%lepton_in%ID = arrID2(idFamily)
       e%lepton_in%mass = 0.0
       e%lepton_in%charge = 0
       e%lepton_in%antiparticle = (idProcess.lt.0)
       e%lepton_in%momentum = (/0.0, 0.0, 0.0, 0.0/)

       e%lepton_out = e%lepton_in

       call setToDefault(e%boson)
       e%boson%ID = Zboson

    case default
       write(*,*) 'idProcess=',idProcess
       call TRACEBACK('wrong idProcess!')

    end select

    ! set the target: (here only for a default free nucleon)

    e%nucleon%ID = 1
    e%nucleon%mass  = mN
    e%nucleon%charge = 1
    e%nucleon%momentum  = (/e%nucleon%mass, 0.0, 0.0, 0.0/)

    e%nucleon_free = e%nucleon
    e%nucleon_free%position=999999999. ! =Vacuum


    ! set the additional variables:

    e%idProcess = idProcess
    e%idFamily  = idFamily

    ! set the masses for the internal kinematical calculations:

    M(1) = e%lepton_in%mass
    M(2) = e%nucleon%mass

    M2(1) = M(1)**2
    M2(2) = M(2)**2

  end subroutine eNev_SetProcess

  !****************************************************************************
  !****s* eN_event/eNev_init_sWQ
  ! NAME
  ! subroutine eNev_init_sWQ(e, sqrts,W,Q2, flagOK)
  !
  ! PURPOSE
  ! Initialize an instance of "type(electronNucleon_event)". Here via
  ! given sqrts,W,Q2.
  ! The target particle is a free nucleon.
  !****************************************************************************
  subroutine eNev_init_sWQ(e, sqrts,W,Q2, flagOK)
    type(electronNucleon_event), intent(inout) :: e
    real, intent(in) :: sqrts,W,Q2
    logical, intent(out) :: flagOK

    real :: nu

    nu = (W**2+Q2-M2(2))/(2*M(2))
    call eNev_init_snQ(e, sqrts,nu,Q2, flagOK)
  end subroutine eNev_init_sWQ

  !****************************************************************************
  !****s* eN_event/eNev_init_sxQ
  ! NAME
  ! suroutine eNev_init_sxQ(e, sqrts,xB,Q2, flagOK)
  !
  ! PURPOSE
  ! Initialize an instance of "type(electronNucleon_event)". Here via
  ! given sqrts,xB,Q2.
  ! The target particle is a free nucleon.
  !****************************************************************************
  subroutine eNev_init_sxQ(e, sqrts,xB,Q2, flagOK)
    type(electronNucleon_event), intent(inout) :: e
    real, intent(in) :: sqrts,xB,Q2
    logical, intent(out) :: flagOK

    real :: nu

    nu = Q2/(2*M(2)*xB)
    call eNev_init_snQ(e, sqrts,nu,Q2, flagOK)
  end subroutine eNev_init_sxQ

  !****************************************************************************
  !****s* eN_event/eNev_init_snQ
  ! NAME
  ! suroutine eNev_init_snQ(e, sqrts,nu,Q2, flagOK)
  !
  ! PURPOSE
  ! Initialize an instance of "type(electronNucleon_event)". Here via
  ! given sqrts,nu,Q2.
  ! The target particle is a free nucleon.
  ! NOTES
  ! * This is the one routine from all the same ones labeled by "..._s?Q",
  !   which does the actual work
  ! * The used kinematical relations can be derived via four-momenta squared:
  !     m^2 = (p_out)^2 = (p_in-q)^2 = m^2 - Q^2 - 2 p0 q0 + 2 pZ qZ
  !   and
  !     pX^2 = E^2 - m^2 - pZ^2
  ! * You see: this routine is only valid for EM and NC case !
  ! * Please note: the case pX^2<0 corresponds to the boundary
  !     Q^2 > 4 E (E-nu)   (for m->0)
  !   which represents the case, that the 3-momenta of lepton and photon are
  !   identical. Violation of this boundary produces unphysical kinematics.
  !
  !****************************************************************************
  subroutine eNev_init_snQ(e, sqrts,nu,Q2, flagOK)
    use minkowski, only: abs4,abs4Sq

    type(electronNucleon_event), intent(inout) :: e
    real, intent(in) :: sqrts,nu,Q2
    logical, intent(out) :: flagOK

    real :: Ebeam,pZ,pX

    flagOK = .false.

    ! in the Lab system: (Target at rest)

    Ebeam = (sqrts**2-M2(1)-M2(2))/(2*M(2))
    pZ = (2*Ebeam*nu + Q2) / (2*sqrt(nu**2+Q2))

    pX = ( Ebeam**2-M2(1)-pZ**2 )
    if (pX<0) return ! --- failure
    pX = sqrt( pX )

    ! set the incoming/outgoing lepton and the exchanged boson:

    e%lepton_in%momentum  = (/Ebeam,    pX, 0E0, pZ/)
    e%lepton_out%momentum = (/Ebeam-nu, pX, 0E0, pZ-sqrt(nu**2+Q2)/)

    e%boson%momentum = e%lepton_in%momentum-e%lepton_out%momentum


    ! calculate W etc...
    e%QSquared = -abs4Sq(e%boson%momentum)
    e%W = abs4(e%boson%momentum+e%nucleon%momentum)

    ! calculate W_free:  inv. mass for nucleon without potential

    e%nucleon_free      = e%nucleon

    e%W_free = abs4(e%boson%momentum+e%nucleon_free%momentum)

    ! calculate Wrec = inv. Mass for nucleon at rest
    e%W_rec = sqrt(mN**2+2.*mN*e%boson%momentum(0)-e%QSquared)

    flagOK = .true.
  end subroutine eNev_init_snQ

  !****************************************************************************
  !****s* eN_event/eNev_init_eWQ
  ! NAME
  ! subroutine eNev_init_eWQ(e, eps,W,Q2, flagOK)
  !
  ! PURPOSE
  ! Initialize an instance of "type(electronNucleon_event)". Here via
  ! given eps,W,Q2.
  ! The target particle is a free nucleon.
  !****************************************************************************
  subroutine eNev_init_eWQ(e, eps,W,Q2, flagOK)
    type(electronNucleon_event), intent(inout) :: e
    real, intent(in) :: eps,W,Q2
    logical, intent(out) :: flagOK

    real :: nu

    nu = (W**2+Q2-M2(2))/(2*M(2))
    call eNev_init_enQ(e, eps,nu,Q2, flagOK)
  end subroutine eNev_init_eWQ

  !****************************************************************************
  !****s* eN_event/eNev_init_exQ
  ! NAME
  ! subroutine eNev_init_exQ(e, eps,xB,Q2, flagOK)
  !
  ! PURPOSE
  ! Initialize an instance of "type(electronNucleon_event)". Here via
  ! given eps,xB,Q2.
  ! The target particle is a free nucleon.
  !****************************************************************************
  subroutine eNev_init_exQ(e, eps,xB,Q2, flagOK)
    type(electronNucleon_event), intent(inout) :: e
    real, intent(in) :: eps,xB,Q2
    logical, intent(out) :: flagOK

    real :: nu

    nu = Q2/(2*M(2)*xB)
    call eNev_init_enQ(e, eps,nu,Q2, flagOK)
  end subroutine eNev_init_exQ

  !****************************************************************************
  !****s* eN_event/eNev_init_enQ
  ! NAME
  ! subroutine eNev_init_enQ(e, eps,nu,Q2, flagOK)
  !
  ! PURPOSE
  ! Initialize an instance of "type(electronNucleon_event)". Here via
  ! given eps,nu,Q2.
  ! The target particle is a free nucleon.
  !****************************************************************************
  subroutine eNev_init_enQ(e, eps,nu,Q2, flagOK)
    type(electronNucleon_event), intent(inout) :: e
    real, intent(in) :: eps,nu,Q2
    logical, intent(out) :: flagOK

    real :: Ebeam,s

    Ebeam = (nu+sqrt((nu**2+Q2)*(1+eps)/(1-eps)))/2
    s = 2*M(2)*Ebeam +M2(1)+M2(2)
    call eNev_init_snQ(e, sqrt(s),nu,Q2, flagOK)
  end subroutine eNev_init_enQ

  !****************************************************************************
  !****s* eN_event/eNev_init_BWQ
  ! NAME
  ! subroutine eNev_init_BWQ(e, Ebeam,W,Q2, flagOK)
  !
  ! PURPOSE
  ! Initialize an instance of "type(electronNucleon_event)". Here via
  ! given Ebeam,W,Q2.
  ! The target particle is a free nucleon.
  ! NOTES
  ! This routine is used by initHiLepton.
  !****************************************************************************
  subroutine eNev_init_BWQ(e, Ebeam,W,Q2, flagOK)
    type(electronNucleon_event), intent(inout) :: e
    real, intent(in) :: Ebeam,W,Q2
    logical, intent(out) :: flagOK

    real :: s

    s = 2*M(2)*Ebeam +M2(1)+M2(2)
    call eNev_init_sWQ(e, sqrt(s),W,Q2, flagOK)
  end subroutine eNev_init_BWQ

  !****************************************************************************
  !****s* eN_event/eNev_init_BnQ
  ! NAME
  ! subroutine eNev_init_BnQ(e, Ebeam,nu,Q2, flagOK)
  !
  ! PURPOSE
  ! Initialize an instance of "type(electronNucleon_event)". Here via
  ! given Ebeam,nu,Q2.
  ! The target particle is a free nucleon.
  ! NOTES
  ! This routine is used by initHiLepton.
  !****************************************************************************
  subroutine eNev_init_BnQ(e, Ebeam,nu,Q2, flagOK)
    type(electronNucleon_event), intent(inout) :: e
    real, intent(in) :: Ebeam,nu,Q2
    logical, intent(out) :: flagOK

    real :: s

    s = 2*M(2)*Ebeam +M2(1)+M2(2)
    call eNev_init_snQ(e, sqrt(s),nu,Q2, flagOK)
  end subroutine eNev_init_BnQ

  !****************************************************************************
  !****s* eN_event/eNev_init_Target
  ! NAME
  ! subroutine eNev_init_Target(e,pTarget,flagOK)
  !
  ! PURPOSE
  ! Initialize an instance of "type(electronNucleon_event)".
  ! The purpose of this routine is to specify the target nucleon
  ! and/or to calculate W etc..
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: e -- the event to modify
  ! * type(particle), OPTIONAL :: pTarget -- target nucleon
  ! OUTPUT
  ! * type(electronNucleon_event) :: e -- the event to modify
  ! * logical :: flagOK -- .true., if everything was okay
  !
  ! NOTES
  ! You may call this routine, if all momenta are set. If yo do not give a
  ! new target nucleon, the nucleon as given in e0 is used.
  !****************************************************************************
  subroutine eNev_init_Target(e,pTarget,flagOK)
    use minkowski, only: abs4Sq
    use LorentzTrafo
    use rotation
    use CallStack
    use output

    type(electronNucleon_event), intent(inout)    :: e
    type(particle), intent(in), optional :: pTarget
    logical       , intent(out):: flagOK

    real,dimension(0:3) :: ptot,pB
    real :: theta,phi,gamma

    logical,parameter :: verbose = .false.

    if (initFlag) call readInput

    flagOK = .false.

    ! set the target, if desired:
    if (present(pTarget)) e%nucleon  = pTarget

    ! calculate W etc...
    e%QSquared = -abs4Sq(e%boson%momentum)
    ptot = e%boson%momentum+e%nucleon%momentum
    e%betacm = ptot(1:3)/ptot(0)
    e%W = abs4Sq(ptot)

    if (e%W .lt. 0.0) then
       if (verbose) then
          write(*,*)
          write(*,*) 'PROBLEM: eNev_init_Target: W^2 < 0'
          write(*,*) e%boson%momentum
          write(*,*) e%nucleon%momentum
          write(*,*) e%betacm
          write(*,*) e%W,e%QSquared
       end if
       return
    end if
    e%W = sqrt(e%W)

    gamma = 1.0 - Dot_Product(e%betacm,e%betacm)
    if (gamma .le. 0.0) then
       if (verbose) then
          write(*,*) 'PROBLEM: eNev_init_Target: gamma < 0'
          write(*,*) e%boson%momentum
          write(*,*) e%nucleon%momentum
          write(*,*) e%betacm
          write(*,*) e%W,e%QSquared
       end if
       return
    end if


    e%pcm = e%nucleon%momentum
    call lorentz(e%betacm,e%pcm)
    e%pcm=-e%pcm !for COLL_GammaN_Py

    phi = atan2(e%pcm(2),e%pcm(1))
    theta = atan2(sqrt(e%pcm(1)**2+e%pcm(2)**2),e%pcm(3))

    pB = e%lepton_in%momentum

    call lorentz(e%betacm,pB)
    pB(1:3) = rotateZY (theta, phi, pB(1:3))

    e%phiLepton = atan2(pB(2),pB(1))

!    write(*,'(A,1P,4e13.5)') 'pB (init)',pB
!    write(*,*) 'e%phiLepton = ',e%phiLepton

    ! calculate W_free:

    select case (selectFrame)
    case (doNOT)
       call RemovePot_doNOT(e)
    case (CALC)
       call RemovePot_CALC(e)
    case (CM)
       call RemovePot_CM(e)
    case (THRE)
       call RemovePot_THRE(e)
    case (NucleonRest)
       call RemovePot_NucleonRest(e)
    case (THRE2)
       call RemovePot_THRE2(e)
    case default
       write(*,*) "selectFrame=",selectFrame
       call TRACEBACK()
    end select

    flagOK = .true.
    return
  end subroutine eNev_init_Target

  !****************************************************************************
  !****s* eN_event/eNev_init_nuStep1
  ! NAME
  ! subroutine eNev_init_nuStep1(e,pTarget)
  !
  ! PURPOSE
  ! Step 1 of "neutrino init": Set the target nucleon.
  ! In Step 0 we have set the type of incoming and outgoing lepton,
  ! but not the kinematics. Therefore we can not use "init_Target" here.
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: e -- the event to modify
  ! * type(particle) :: pTarget -- target nucleon
  ! OUTPUT
  ! * type(electronNucleon_event) :: e -- the event to modify
  !
  !****************************************************************************
  subroutine eNev_init_nuStep1(e,pTarget)
    type(electronNucleon_event), intent(inout) :: e
    type(particle), intent(in) :: pTarget

    e%nucleon  = pTarget
  end subroutine eNev_init_nuStep1

  !****************************************************************************
  !****s* eN_event/eNev_init_nuStep2
  ! NAME
  ! subroutine eNev_init_nuStep2(e,Enu,IP,flagOK)
  !
  ! PURPOSE
  ! Step 2 of "neutrino init": Set the kinematic of the incoming neutrino.
  ! Additionaly, it is checked whether the energy is above threshold
  !
  ! The neutrino is aligned along the z-axis.
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: e -- the event to modify
  ! * real :: Enu -- Energy of incoming lepton
  ! * integer :: IP -- ID of process (needed for check)
  ! OUTPUT
  ! * type(electronNucleon_event) :: e -- the event to modify
  ! * logical :: flagOK -- .true., if everything was okay
  !****************************************************************************
  subroutine eNev_init_nuStep2(e,Enu,IP,flagOK)

    use idtable, only: nucleon
    use constants, only: mPi

    type(electronNucleon_event), intent(inout) :: e
    real, intent(in) :: Enu
    integer, intent(in),OPTIONAL :: IP
    logical, intent(out),OPTIONAL :: flagOK

    real :: meff, EnuThres,finalstatemass2

    e%lepton_in%momentum = (/Enu,    0.0, 0.0, Enu/)

    if (present(flagOK)) then
       flagOK = .false.

       meff=sqrtS(e%nucleon)
       if (IP.eq.nucleon) then
          finalstatemass2=(e%lepton_out%mass+meff)**2
       else
          finalstatemass2=(e%lepton_out%mass+mPi+meff)**2
       end if

       EnuThres=(finalstatemass2-meff**2)/(2*(e%nucleon%momentum(0)-e%nucleon%momentum(3)))
       if (Enu.ge.EnuThres) then
          flagOK = .TRUE.
       end if
    end if

  end subroutine eNev_init_nuStep2

  !****************************************************************************
  !****s* eN_event/eNev_init_nuStep3a
  ! NAME
  ! subroutine eNev_init_nuStep3a(e,Eprime,costheta,flagOK)
  !
  ! PURPOSE
  ! Step 3 of "neutrino init": Set the kinematic of the outgoing lepton.
  ! Now we have fixed everything (comparable to "eNev_init_Target")
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: e -- the event to modify
  ! * real :: Eprime -- Energy of scattered lepton
  ! * real :: costheta -- Angle of scattered lepton
  ! OUTPUT
  ! * type(electronNucleon_event) :: e -- the event to modify
  ! * logical :: flagOK -- .true., if everything was okay
  !
  !****************************************************************************
  subroutine eNev_init_nuStep3a(e,Eprime,costheta,flagOK)

    use random, only: rn
    use constants, only: pi
!    use minkowski, only : abs4,abs4Sq
!    use LorentzTrafo
!    use rotation
!    use output

    type(electronNucleon_event), intent(inout) :: e
    real, intent(in) :: Eprime,costheta
    logical, intent(out) :: flagOK

    real :: phi, pPrime, sintheta

    flagOK = .FALSE.

    if (Eprime.lt.e%lepton_out%mass) return ! threshold violated !

    phi=rn()*2.*pi ! random generation of azimutal angle
!!    phi=0.0 !!!! FOR TEMPORAL USE ONLY !!!!

    pPrime = sqrt(Eprime**2 - e%lepton_out%mass**2)
    sintheta = sqrt(max((1.-costheta**2),0.))

    e%lepton_out%momentum(0)=Eprime
    e%lepton_out%momentum(1)=pPrime*sintheta*cos(phi)
    e%lepton_out%momentum(2)=pPrime*sintheta*sin(phi)
    e%lepton_out%momentum(3)=pPrime*costheta

    e%boson%momentum = e%lepton_in%momentum-e%lepton_out%momentum

    call eNev_init_Target(e,flagOK=flagOK)


  end subroutine eNev_init_nuStep3a

  !****************************************************************************
  !****s* eN_event/eNev_init_nuStep3b
  ! NAME
  ! subroutine eNev_init_nuStep3b(e,Eprime,Q2,flagOK)
  !
  ! PURPOSE
  ! Step 3 of "neutrino init": Set the kinematic if the outgoing lepton.
  ! Now we have fixed everything (comparable to "eNev_init_Target")
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: e -- the event to modify
  ! * real :: Eprime -- Energy of scattered lepton
  ! * real :: Q2 -- 4-mom squared of boson
  ! OUTPUT
  ! * type(electronNucleon_event) :: e -- the event to modify
  ! * logical :: flagOK -- .true., if everything was okay
  !****************************************************************************
  subroutine eNev_init_nuStep3b(e,Eprime,Q2,flagOK)
    type(electronNucleon_event), intent(inout) :: e
    real, intent(in) :: Eprime,Q2
    logical, intent(out) :: flagOK

    real :: cost, ml_out, E_in,h1,h2,m2

    flagOK = .FALSE.

    ml_out = e%lepton_out%mass
    E_in = e%lepton_in%momentum(0)

    if (Eprime.lt.ml_out) return ! failure !!!

    m2 = ml_out**2
    h1 = 2*E_in*Eprime
    h2 = 2*E_in*sqrt(Eprime**2-m2)
    cost=-(m2+Q2-h1)/h2
    if (cost.gt. 1.0) return   ! failure !!! (Q2_min)
    if (cost.lt.-1.0) return   ! failure !!! (Q2_max)

    call eNev_init_nuStep3a(e,Eprime,cost,flagOK)

  end subroutine eNev_init_nuStep3b


  !****************************************************************************
  !****s* eN_event/eNev_init_nuStep3c
  ! NAME
  ! subroutine eNev_init_nuStep3c(e,x,y,flagOK)
  !
  ! PURPOSE
  ! Step 3 of "neutrino init": Set the kinematic if the outgoing lepton.
  ! Now we have fixed everything (comparable to "eNev_init_Target")
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: e -- the event to modify
  ! * real :: x,y -- Bjorken variables
  ! OUTPUT
  ! * type(electronNucleon_event) :: e -- the event to modify
  ! * logical :: flagOK -- .true., if everything was okay
  !
  ! NOTES
  ! we use the lorentz invariant definitions
  ! x = -(qq)/2(pq),  y = (pq)/(pk)
  ! and not the assumptions for resting target nucleon,
  ! x -> Q^2/2Mnu,  y = nu/E .
  !****************************************************************************
  subroutine eNev_init_nuStep3c(e,x,y,flagOK)
    use random, only: rn
    use minkowski, only: SP
    use constants, only: pi
    use CallStack

    type(electronNucleon_event), intent(inout) :: e
    real, intent(in) :: x,y
    logical, intent(out) :: flagOK

    real :: ml_out, E_in, phi
    real :: a,b,c,DD, qT
    real :: PK,PQ,QQ

    logical, parameter :: verbose = .false.

    flagOK = .FALSE.

    ml_out = e%lepton_out%mass
    E_in = e%lepton_in%momentum(0)

    phi=rn()*2.*pi ! random generation of azimutal angle
    !!    phi=0.0 !!!! FOR TEMPORAL USE ONLY !!!!

    PK = SP(e%lepton_in%momentum,e%nucleon%momentum)
    PQ = y * PK
    QQ = -2*x*y * PK

    if (verbose) then
       write(*,*)
       write(*,*) 'x,y     :',x,y
       write(*,*) 'PK,PQ,QQ:',PK,PQ,QQ
    end if

    a = (ml_out**2-QQ)/(2*E_in)
    b = (e%nucleon%momentum(1)*cos(phi)+e%nucleon%momentum(2)*sin(phi))/(e%nucleon%momentum(0)-e%nucleon%momentum(3))
    c = (PQ+e%nucleon%momentum(3)*a)/(e%nucleon%momentum(0)-e%nucleon%momentum(3))

    DD = (a*b)**2-a**2-2*a*c-QQ
!    if (DD.lt.0.0) write(*,*) 'DD < 0'
    if (DD.lt.0.0) return

    qT = sqrt(DD) - a*b
!    if (qT.lt.0.0) call TRACEBACK('qT < 0')
    if (qT.lt.0.0) return

!    if (-sqrt(DD) - a*b.gt.0.0) call TRACEBACK('strange: -...-... > 0 !',-1)

    e%boson%momentum = (/b*qT+c,qT*cos(phi),qT*sin(phi),b*qT+c+a/)
!    if (e%boson%momentum(0).lt.0.0) call TRACEBACK('nu < 0')
    if (e%boson%momentum(0).lt.0.0) return

    e%lepton_out%momentum = e%lepton_in%momentum-e%boson%momentum
!    if (e%lepton_out%momentum(0).lt.ml_out) call TRACEBACK('E_out < m')
    if (e%lepton_out%momentum(0).lt.ml_out) return

    call eNev_init_Target(e,flagOK=flagOK)

    if (verbose) call write_electronNucleon_event(e)


  end subroutine eNev_init_nuStep3c


  !****************************************************************************
  !****s* eN_event/init_electronNucleon_event
  ! NAME
  ! subroutine init_electronNucleon_event(e, k_in,k_out,nuc,flagOK)
  !
  ! PURPOSE
  ! Initialize an instance of "type(electronNucleon_event)". Here via
  ! incoming and outgoing lepton momenta
  ! OUTPUT
  ! * type(electronNucleon_event) :: e -- fully initialized event
  ! * logical, optional :: flagOK -- .true. if all okay
  !
  ! NOTES
  ! * This is used by Olli for the low energy lepton routines
  ! * This routine is only valid for EM and e-
  !****************************************************************************
  subroutine init_electronNucleon_event(e, k_in,k_out,nuc,flagOK)
    use particleDefinition

    type(electronNucleon_event), intent(inout) :: e
    type(particle)      , intent(in)  :: nuc
    real, dimension(0:3), intent(in)  :: k_in, k_out

    logical, optional   :: flagOK

    logical :: flag

    if (present(flagOK)) flagOK=.false.

    if (initFlag) call readInput

    ! set the incoming/outgoing lepton and the exchanged boson:

    call eNev_SetProcess(e, 1,1) ! This routine is only valid for EM and e-

    e%lepton_in%momentum  = k_in
    e%lepton_out%momentum = k_out
    e%boson%momentum      = k_in-k_out

    call eNev_init_Target(e,nuc,flag)
    if (present(flagOK)) flagOK=flag

  end subroutine init_electronNucleon_event

  !****************************************************************************
  !****is* eN_event/RemovePot_DoNOT
  ! NAME
  ! function RemovePot_DoNOT(e)
  !
  ! PURPOSE
  ! Calculate the "FREE" kinematics by assuming the same kinematics as in the
  ! InMedium case. I.e., we do not remove the potential at all.
  !****************************************************************************
  subroutine RemovePot_DoNOT(e)
    use minkowski, only: abs4
    type(electronNucleon_event) :: e
    e%nucleon_free               = e%nucleon
    e%nucleon_free%position      =999999999. ! =Vacuum if the nucleus rests at (0,0,0)
    e%W_free=abs4(e%boson%momentum+e%nucleon%momentum)
    ! calculate Wrec = inv. Mass for nucleon at rest
    e%W_rec = sqrt(mN**2+2.*mN*e%boson%momentum(0)-e%QSquared)
    e%W = abs4(e%boson%momentum+e%nucleon%momentum)
  end subroutine RemovePot_DoNOT

  !****************************************************************************
  !****is* eN_event/RemovePot_THRE
  ! NAME
  ! function RemovePot_THRE(e)
  !
  ! PURPOSE
  ! Calculate the "FREE" kinematics by setting the free W according
  !   W_free = W - m* + m ,
  ! i.e. assuring the threshold behaviour.
  ! The momenta of the 'free' target nucleon are set in the CM frame.
  !****************************************************************************
  subroutine RemovePot_THRE(e)
    use minkowski, only: abs4
    use lorentzTrafo

    type(electronNucleon_event) :: e
!!$    real :: pvec_old,  pvec, qvec, q0, costheta, Q2, a,b,c
    real :: pvec_old,  pvec

    e%nucleon_free               = e%nucleon
    e%nucleon_free%position      =999999999. ! =Vacuum if the nucleus rests at (0,0,0)
    e%W_free=abs4(e%boson%momentum+e%nucleon%momentum)-abs4(e%nucleon%momentum)+e%nucleon%mass
    ! calculate Wrec = inv. Mass for nucleon at rest
    e%W_rec = sqrt(mN**2+2.*mN*e%boson%momentum(0)-e%QSquared)
    e%W = abs4(e%boson%momentum+e%nucleon%momentum)

    ! nucleon momentum should be consistent with W_free !but this does not work in Pythia

!!$   pvec=sqrt(Dot_Product(e%nucleon%momentum(1:3),e%nucleon%momentum(1:3)))
!!$   pvec_old=pvec
!!$   qvec=sqrt(Dot_Product(e%boson%momentum(1:3),e%boson%momentum(1:3)))
!!$   costheta=Dot_Product(e%boson%momentum(1:3),e%nucleon%momentum(1:3))/pvec/qvec
!!$   q0=e%boson%momentum(0)
!!$   Q2=-abs4Sq(e%boson%momentum)
!!$   a=4*(qvec**2*costheta**2-q0**2)
!!$   c=e%W_free**2+Q2-e%nucleon%mass**2
!!$   b=2*qvec*costheta*c
!!$   c=c**2-4.*q0**2*e%nucleon%mass**2
!!$   pvec=max(  ( -b/2.+sqrt(b**2/4.-a*c) )/a, ( -b/2.-sqrt(b**2/4.-a*c) )/a )
!!$   e%nucleon_free%momentum(0) = sqrt(pvec**2+e%nucleon%mass**2)
!!$   e%nucleon_free%momentum(1:3) = e%nucleon%momentum(1:3)/pvec_old*pvec

    ! try another recipe:

    call lorentz(e%betacm,e%boson%momentum) ! boost to CM
    call lorentz(e%betacm,e%nucleon_free%momentum)
    pvec_old=sqrt(Dot_Product(e%nucleon_free%momentum(1:3),e%nucleon_free%momentum(1:3))) ! old absolute value of momentum in the CM
    e%nucleon_free%momentum(0)=max(e%W_free-e%boson%momentum(0),e%nucleon_free%mass) ! new nucleon energy in the CM
    pvec=max(0.0,sqrt(e%nucleon_free%momentum(0)**2-e%nucleon_free%mass**2)) ! new absolute value of momentum in the CM
    e%nucleon_free%momentum(1:3)=e%nucleon_free%momentum(1:3)/pvec_old*pvec ! recalculate 3-momemtum in CM
    call lorentz(-e%betacm,e%boson%momentum) ! boost back
    call lorentz(-e%betacm,e%nucleon_free%momentum)

  end subroutine RemovePot_THRE



  !****************************************************************************
  !****is* eN_event/RemovePot_THRE2
  ! NAME
  ! function RemovePot_THRE2(e)
  !
  ! PURPOSE
  ! Calculate the "FREE" kinematics by setting the free W according
  !   W_free^2 = W^2 - m*^2 + m^2 ,
  ! i.e. assuring the threshold behaviour as expected for neutrinos.
  ! The momenta of the 'free' target nucleon are set in the CM frame.
  !****************************************************************************
  subroutine RemovePot_THRE2(e)
    use minkowski, only: abs4Sq,abs4
    use lorentzTrafo

    type(electronNucleon_event) :: e
!!$    real :: pvec_old,  pvec, qvec, q0, costheta, Q2, a,b,c
    real :: pvec_old,  pvec, W2free

    e%nucleon_free               = e%nucleon
    e%nucleon_free%position      =999999999. ! =Vacuum if the nucleus rests at (0,0,0)
    W2free=abs4Sq(e%boson%momentum+e%nucleon%momentum)-abs4Sq(e%nucleon%momentum)+e%nucleon%mass**2
    if (W2free>0) then
                   e%W_free=sqrt(W2free)
                   ! calculate Wrec = inv. Mass for nucleon at rest
                   e%W_rec = sqrt(mN**2+2.*mN*e%boson%momentum(0)-e%QSquared)
                   e%W = abs4(e%boson%momentum+e%nucleon%momentum)
    else
            write(*,*) 'problems in removing potential THRE2 W2free=',W2free, '   STOP'
            stop
    end if

    ! nucleon momentum should be consistent with W_free !but this does not work in Pythia

!!$   pvec=sqrt(Dot_Product(e%nucleon%momentum(1:3),e%nucleon%momentum(1:3)))
!!$   pvec_old=pvec
!!$   qvec=sqrt(Dot_Product(e%boson%momentum(1:3),e%boson%momentum(1:3)))
!!$   costheta=Dot_Product(e%boson%momentum(1:3),e%nucleon%momentum(1:3))/pvec/qvec
!!$   q0=e%boson%momentum(0)
!!$   Q2=-abs4Sq(e%boson%momentum)
!!$   a=4*(qvec**2*costheta**2-q0**2)
!!$   c=e%W_free**2+Q2-e%nucleon%mass**2
!!$   b=2*qvec*costheta*c
!!$   c=c**2-4.*q0**2*e%nucleon%mass**2
!!$   pvec=max(  ( -b/2.+sqrt(b**2/4.-a*c) )/a, ( -b/2.-sqrt(b**2/4.-a*c) )/a )
!!$   e%nucleon_free%momentum(0) = sqrt(pvec**2+e%nucleon%mass**2)
!!$   e%nucleon_free%momentum(1:3) = e%nucleon%momentum(1:3)/pvec_old*pvec

    ! try another receipt:

    call lorentz(e%betacm,e%boson%momentum) ! boost to CM
    call lorentz(e%betacm,e%nucleon_free%momentum)
    pvec_old=sqrt(Dot_Product(e%nucleon_free%momentum(1:3),e%nucleon_free%momentum(1:3))) ! old absolute value of momentum in the CM
    e%nucleon_free%momentum(0)=max(e%W_free-e%boson%momentum(0),e%nucleon_free%mass) ! new nucleon energy in the CM
    pvec=max(0.0,sqrt(e%nucleon_free%momentum(0)**2-e%nucleon_free%mass**2)) ! new absolute value of momentum in the CM
    e%nucleon_free%momentum(1:3)=e%nucleon_free%momentum(1:3)/pvec_old*pvec ! recalculate 3-momemtum in CM
    call lorentz(-e%betacm,e%boson%momentum) ! boost back
    call lorentz(-e%betacm,e%nucleon_free%momentum)

  end subroutine RemovePot_THRE2


  !****************************************************************************
  !****is* eN_event/RemovePot_CALC
  ! NAME
  ! function RemovePot_CALC(e)
  !
  ! PURPOSE
  ! Calculate the "FREE" kinematics, while removing the potential in the CALC
  ! frame
  !****************************************************************************
  subroutine RemovePot_CALC(e)
    use minkowski, only: abs4

    type(electronNucleon_event) :: e

    e%nucleon_free               = e%nucleon
    e%nucleon_free%momentum(1:3) = e%nucleon%momentum(1:3)
    e%nucleon_free%momentum(0)   = FreeEnergy(e%nucleon)
    e%nucleon_free%position      =999999999. ! =Vacuum if the nucleus rests at (0,0,0)
    e%W_free=abs4(e%boson%momentum+e%nucleon_free%momentum)
    ! calculate Wrec = inv. Mass for nucleon at rest
    e%W_rec = sqrt(mN**2+2.*mN*e%boson%momentum(0)-e%QSquared)
    e%W = abs4(e%boson%momentum+e%nucleon%momentum)

  end subroutine RemovePot_CALC

  !****************************************************************************
  !****is* eN_event/RemovePot_CM
  ! NAME
  ! function RemovePot_CM(e)
  !
  ! PURPOSE
  ! Calculate the "FREE" kinematics, while removing the potential in the CM
  ! frame
  ! NOTES
  ! All the discussions in Juergen Lehr, PhD thesis, pages 159ff are obsolete
  !****************************************************************************
  subroutine RemovePot_CM(e)
    use minkowski, only: abs4
    use lorentzTrafo

    type(electronNucleon_event) :: e

    e%nucleon_free               = e%nucleon
    e%nucleon_free%position      =999999999

    e%nucleon_free%momentum = e%nucleon%momentum
    call lorentz(e%betacm,e%nucleon_free%momentum)
    e%nucleon_free%momentum(0)   = FreeEnergy(e%nucleon_free)
    call lorentz(-e%betacm,e%nucleon_free%momentum)

    e%W_free=abs4(e%boson%momentum+e%nucleon_free%momentum)
    ! calculate Wrec = inv. Mass for nucleon at rest
    e%W_rec = sqrt(mN**2+2.*mN*e%boson%momentum(0)-e%QSquared)
    e%W = abs4(e%boson%momentum+e%nucleon%momentum)

  end subroutine RemovePot_CM


  !****************************************************************************
  !****is* eN_event/RemovePot_NucleonRest
  ! NAME
  ! function RemovePot_NucleonRest(e)
  !
  ! PURPOSE
  ! Calculate the "FREE" kinematics, while removing the potential
  ! and adjusting the photon 4-momentum in the nucleon rest frame
  ! NOTES
  ! Recipe proposed on a group meeting May 2011
  !****************************************************************************
  subroutine RemovePot_NucleonRest(e)
    use minkowski, only: abs3, abs4, abs4Sq, SP
    use lorentzTrafo, only: lorentz

    type(electronNucleon_event) :: e

    real, dimension(1:3) :: beta ! velocity of the nucleon

    real, dimension(0:3) :: boson_momentum    ! modified momentum of the boson
    real                 :: qvectilde_abs, qvec_abs

    e%nucleon_free               = e%nucleon
    e%nucleon_free%position      =999999999

    ! velocity of nucleon in lab frame:
    beta(1:3)=e%nucleon%momentum(1:3)/e%nucleon%momentum(0)

    ! now define momentum of unbound nucleon in lab frame:
    e%nucleon_free%momentum(1:3) = 0.
    e%nucleon_free%momentum(0) = 0.938
    ! now nucleon_free is an unbound nucleon with free mass at rest

    ! boson momentum in lab frame:
    boson_momentum = e%boson%momentum
    ! now boost boson momentum to  nucleon rest frame
    call lorentz(beta,boson_momentum)

    ! calculate Wrec = inv. Mass for free nucleon at rest
    e%W_rec = sqrt(mN**2+2.*mN*e%boson%momentum(0)-e%QSquared)
    ! calculate W = inv. Mass for bound, moving nucleon
    e%W = abs4(e%boson%momentum+e%nucleon%momentum)
    ! calculate W = inv. Mass for bound, moving nucleon
    e%W_free=abs4(boson_momentum+e%nucleon_free%momentum)

  end subroutine RemovePot_NucleonRest

!  !***************************************************************************
!  !****is* eN_event/RemovePot_NucleonRest
!  ! NAME
!  ! function RemovePot_NucleonRest(e)
!  !
!  ! PURPOSE
!  ! Calculate the "FREE" kinematics, while removing the potential
!  ! and adjusting the photon 4-momentum in the nucleon rest frame
!  ! NOTES
!  ! Recipe proposed on a group meeting May 2011
!  !***************************************************************************
!  subroutine RemovePot_NucleonRest(e)
!    use minkowski, only: abs3, abs4, abs4Sq, SP
!    use lorentzTrafo, only: lorentz
!
!    type(electronNucleon_event) :: e
!
!    real, dimension(1:3) :: beta ! velocity of the nucleon
!    real, dimension(0:3) :: qtilde    ! modified momentum of the boson
!    real                 :: qvectilde_abs, qvec_abs
!
!    e%nucleon_free               = e%nucleon
!    e%nucleon_free%position      =999999999
!
!    beta(1:3)=e%nucleon%momentum(1:3)/e%nucleon%momentum(0) ! velocity of the nucleon
!
!    e%nucleon_free%momentum = e%nucleon%momentum
!    call lorentz(beta,e%nucleon_free%momentum) ! boost to the nucleon rest frame
!    e%nucleon_free%momentum(0)   = e%nucleon_free%mass ! prescription: replace bount mass with the free mass
!
!    qtilde(0)=SP(e%nucleon%momentum,e%boson%momentum)/e%nucleon_free%mass ! prescription: qtilde0 from pq=const
!
!    qvectilde_abs=sqrt(-abs4Sq(e%boson%momentum)+qtilde(0)**2) ! prescription: absolute 3-momentum from Q2=const
!
!    qvec_abs = abs3( e%boson%momentum )
!    qtilde(1:3)=e%boson%momentum(1:3)/qvec_abs*qvectilde_abs ! keep the directio of the 3-momentum, adjust the value
!
!    e%W_free=abs4(qtilde+e%nucleon_free%momentum)
!    ! calculate Wrec = inv. Mass for nucleon at rest
!    e%W_rec = sqrt(mN**2+2.*mN*e%boson%momentum(0)-e%QSquared)
!    e%W = abs4(e%boson%momentum+e%nucleon%momentum)
!
!    e%boson%momentum=qtilde  ! new boson momentum equals qtilde
!
!    ! boost nucleon and boson back to the calc frame
!    call lorentz(-beta,e%nucleon_free%momentum)
!    call lorentz(-beta,e%boson%momentum)
!
!  end subroutine RemovePot_NucleonRest
!


  !****************************************************************************
  !****s* eN_event/eNeV_GetKinV
  ! NAME
  ! subroutine eNeV_GetKinV(e, nu,Q2,W, Wfree, eps, fT)
  !
  ! PURPOSE
  ! return the kinematical variables buried in the event type.
  ! additionally calculate   :
  ! * real :: fT   -- transversal flux fT in GeV^-3
  ! * real :: eps  -- epsilon = fL/fT
  !
  ! NOTES
  ! the flux calculation was taken from PyVP:
  ! * the original formulae were implemented by Thomas Falter (PhD,
  !   eq.(3.15),(3.16)). We modified and corrected them.
  ! * we are now trying to supersede those by a "low-W/low-Q2" expression.
  ! * The equivalent photon energies are (we use the Hand convention)
  !     K = (W^2-M^2)/2M = (1-x)*nu = (1-x)*y*Ebeam  [Hand]
  !     K = sqrt[nu^2+Q^2] [Gilman]
  ! * if the global variable "restingNucleon" is true, then we assume a
  !   resting nucleon for the calculation of fT and eps. Otherwise we use
  !   the target four momentum (the returned value of nu is not influenced
  !   by this)
  ! * the returned flux fT corresponds to (c_T = 1+(1-y)^2)
  !      f_T = \frac{\alpha}{2\pi} \frac{K}{Q^2\nu^2} c_T
  !   In order to get d\sigma/dE'd\cos\theta = \Gamma \sigma^*, one has to use
  !      \Gamma = f_T \frac{E E'}{\pi}
  !****************************************************************************
  subroutine eNeV_GetKinV(e, nu,Q2,W,Wfree,eps,fT)
    use constants, only: alphaQED, twopi
    use minkowski, only: SP

    type(electronNucleon_event),intent(in) :: e
    real, intent(out) :: nu,Q2,W
    real, intent(out), optional :: Wfree, eps, fT

    real :: x,y,Ebeam,K,cL,cT, nu0

    nu  = e%boson%momentum(0)
    W   = e%W
    Q2  = e%QSquared
    if (present(Wfree)) Wfree = e%W_free
    if (present(eps).or.present(fT)) then
       if (restingNucleon) then
          nu0 = nu
          Ebeam = e%lepton_in%momentum(0)
       else
          nu0 =   SP(e%boson%momentum,e%nucleon%momentum)/M(2)
          Ebeam = SP(e%lepton_in%momentum,e%nucleon%momentum)/M(2)
       end if
       y=nu0/Ebeam
       x=Q2/(2*M(2)*nu0)
       K = (1.-x)*nu0 ! Hand convention
       cL = 2.*(1.-y-Q2/(4*Ebeam**2))/(1.+Q2/nu0**2)
       cT = cL + y**2

       if (present(fT))  fT  = alphaQED/twopi * K/(Q2*nu0**2) * cT
       if (present(eps)) eps = cL/cT

    end if

    ! this was implemented by Thomas Falter; PhD, eq.(3.15),(3.16)
    !-------------------------------------------------------------
!!$    K = (1.-x)*nu ! Hand convention
!!$    cL = 2.*(1.-y)
!!$    cT = cL + y**2
!!$    weight = alphaQED/twopi * K/(Q2*nu**2) * cT
!!$    eps = cL/cT

  end subroutine eNeV_GetKinV

  !****************************************************************************
  !****f* eN_event/eNeV_Get_LightX
  ! NAME
  ! real function eNeV_Get_LightX(e)
  !
  ! PURPOSE
  ! return the lightcone x, not Bjorken x !!
  !****************************************************************************
  real function eNeV_Get_LightX(e)
    use minkowski, only: abs4Sq

    type(electronNucleon_event),intent(in) :: e
    real :: s

    s  = abs4sq(e%lepton_in%momentum+e%nucleon_free%momentum)
    eNeV_Get_LightX = (e%W_free**2-M2(2))/(s-M2(2)+M2(1))
  end function eNeV_Get_LightX

  !****************************************************************************
  !****f* eN_event/eNeV_Get_CostLepton
  ! NAME
  ! real function eNeV_Get_CostLepton(e)
  !
  ! PURPOSE
  ! return the cosine of the angle between ncoming and outgoing lepton
  !
  ! NOTES
  ! * TODO: This is probably already defined elsewhere. Replace!
  !****************************************************************************
  real function eNeV_Get_CostLepton(e)

    type(electronNucleon_event),intent(in) :: e
    real :: a,b

    a = Dot_product(e%lepton_in%momentum(1:3), e%lepton_in%momentum(1:3))
    b = Dot_product(e%lepton_out%momentum(1:3),e%lepton_out%momentum(1:3))

    eNeV_Get_CostLepton = Dot_product(e%lepton_in%momentum(1:3), e%lepton_out%momentum(1:3))/sqrt(a*b)

  end function eNeV_Get_CostLepton

  !****************************************************************************
  !****s* eN_event/eNev_Set_PhiLepton
  ! NAME
  ! real function eNeV_Set_PhiLepton(e,phi)
  !
  ! PURPOSE
  ! rotate the lepton vectors around z-axis
  ! NOTES
  ! We assume, that before rotation, the transverse direction was the x-axis
  ! and thee y-component was zero!
  !****************************************************************************
  subroutine eNeV_Set_PhiLepton(eNev,phi)

    type(electronNucleon_event),intent(inout) :: eNev
    real, intent(in) :: phi

    real :: pT

    pT = eNev%lepton_in%momentum(1)
    eNev%lepton_in%momentum(1)=pT*cos(phi)
    eNev%lepton_in%momentum(2)=pT*sin(phi)
    pT = eNev%lepton_out%momentum(1)
    eNev%lepton_out%momentum(1)=pT*cos(phi)
    eNev%lepton_out%momentum(2)=pT*sin(phi)

  end subroutine eNeV_Set_PhiLepton

  !****************************************************************************
  !****f* eN_event/eNev_CheckForDIS
  ! NAME
  ! logical function eNeV_CheckForDIS()
  !
  ! PURPOSE
  ! Apply cuts, which are also tested in PYGAGA and may leed to infinite loops
  ! or program abortion.
  !
  ! NOTES
  ! One of the cuts in PYTHIA we commented out, so it is also not used here
  ! anymore
  !****************************************************************************
  logical function eNeV_CheckForDIS(e)
    use minkowski, only: abs4Sq
    use CollTools, only: PythiaCKIN

    type(electronNucleon_event),intent(in) :: e

    real :: x,s,PCM1,PCM3

    eNeV_CheckForDIS = .true.

    s = abs4sq(e%lepton_in%momentum+e%nucleon_free%momentum)
    x = (e%W_free**2-M2(2))/(s-M2(2)+M2(1)) ! lightcone-x

    ! This cut is not in use anymore:
!    if (e%QSquared .gt. (1-x)*(s-2*M2(2))-(2-x**2)*M2(1)/(1-x)) &
!         & eNeV_CheckForDIS = .false.

    PCM1 = s-M2(2)+M2(1)
    PCM3 = s-M2(2)-M2(1)

!    write(*,*) 'CheckForDIS:',s,x,PCM1,PCM3,(PCM1*x+e%QSquared)/PCM3


    if (((PCM1*x+e%QSquared)/PCM3).gt.PythiaCKIN(74)) &
         & eNeV_CheckForDIS = .false.

    return
  end function eNeV_CheckForDIS

  !****************************************************************************
  !****s* eN_event/eNev_GetLeptonCM
  ! NAME
  ! subroutine eNeV_GetLeptonCM(e, betacm,phi,theta,phiLepton)
  !
  ! PURPOSE
  !
  !
  ! NOTES
  ! The pcm stored in the event type is the pcm of the boson-nucleon system.
  ! Here we return the pcm of the lepton-nucleon system.
  !
  ! This is e.g. used in DoColl_nuN_Py.
  !****************************************************************************
  subroutine eNeV_GetLeptonCM(e, betacm,phi,theta,phiLepton)
    use LorentzTrafo
    use rotation
    use CallStack

    type(electronNucleon_event),intent(in) :: e

    real, dimension(3), intent(out) :: betacm
    real, intent(out) :: phi,theta,phiLepton

    real, dimension(0:3)  :: pcm           ! Lorentz-Trafo into cm-frame


    real, dimension(0:3) :: pTot,pB
    real :: gamma

    logical:: verbose = .true.

    ptot = e%lepton_in%momentum+e%nucleon%momentum
    betacm = ptot(1:3)/ptot(0)
    gamma = 1.0 - Dot_Product(betacm,betacm)
    if (gamma .le. 0.0) then
       if (verbose) then
          write(*,*) 'PROBLEM: eNev_GetLeptonPCM: gamma < 0'
          write(*,*) e%lepton_in%momentum
          write(*,*) e%boson%momentum
          write(*,*) e%nucleon%momentum
          write(*,*) betacm
          write(*,*) e%betacm
          write(*,*) e%W,e%QSquared
       end if
       call TRACEBACK('gamma < 0')
       return
    end if

    pcm = e%nucleon%momentum
    call lorentz(betacm,pcm)
    pcm=-pcm !for COLL_GammaN_Py

    phi = atan2(pcm(2),pcm(1))
    theta = atan2(sqrt(pcm(1)**2+pcm(2)**2),pcm(3))

    pB = e%lepton_out%momentum

    call lorentz(betacm,pB)
    pB(1:3) = rotateZY (theta, phi, pB(1:3))

    phiLepton = atan2(pB(2),pB(1))

  end subroutine eNeV_GetLeptonCM

  !****************************************************************************
  !****f* eN_event/nuclearFluxFactor_correction
  ! NAME
  ! real function nuclearFluxFactor_correction(p_initial,l_initial)
  ! PURPOSE
  ! Evaluates the correction factor due to considering the nuclear flux
  ! factor instead of the
  ! ones of the single nucleons. Outputs |v_e-v_n|=l^\mu p_\mu /(l_0 p_0)
  !
  ! INPUTS
  ! * real, dimension(0:3) :: p_initial, l_initial -- initial nucleon and
  !   lepton momenta
  !
  ! OUTPUT
  ! * real :: nuclearFluxFactor_correction
  !****************************************************************************
  real function nuclearFluxFactor_correction(p_initial,l_initial,withPotential)
    use Minkowski, only: SP
    real, dimension(0:3), intent(in) :: p_initial, l_initial
    logical, optional, intent(in)    :: withPotential
    logical              :: withP
    real, dimension(0:3) :: p
    !!  following the Ulrich-Kai-Alexei-Olga  discussion  on 05.04.2012
    withP=.false.
    if (present(withPotential)) withP=withPotential

    if (.not.withP) then ! redefine  p0_initial using vacuum mass
    p(0)=sqrt(  Dot_Product(p_initial(1:3),p_initial(1:3)) + mN**2  )
    p(1:3)=p_initial(1:3)
    else
    p=p_initial! otherwise use p0_initial  as normally defined including potential (== effective mass)
    end if

    nuclearFluxFactor_correction=SP(l_initial,p)/(l_initial(0)*p(0))

  end function nuclearFluxFactor_correction

end module eN_event
