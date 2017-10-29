!******************************************************************************
!****m* /lepton2p2h
! NAME
! module lepton2p2h
!
! PURPOSE
! Do all the internals for 2p2h scattering:
!
! EM:
! * ell N1 N2 --> ell' N1' N2'  == gamma* N1 N2 --> N1' N2'
! * ell N1 N2 --> ell' N Delta  == gamma* N1 N2 --> N Delta
! NC:
! * nu  N1 N2 --> nu'  N1' N2'
! * nu  N1 N2 --> nu'  N Delta
! CC:
! * nu  N1 N2 --> ell- N1' N2' (sum of hadronic charges increases by +1)
! * nu  N1 N2 --> ell- N Delta (  -- " --                              )
!
! antiEM, antiNC and antiCC are the same as EM, NC, CC.
!
! cases 1 - 3 give parametrizations for 2p2h part of structure function W1
! in terms of Q^2, no distinction for neutrinos and antineutrinos
! Cases 1 and 2 are those in: Lalakulich Gallmeister Mosel PRC86(2012)014614
! Case 3 gives a reasonable description of MiniBooNE dd neutrino data
! cases 4 - 5 give parametrization for MEC part of W1 from Christy and Bosted
! In this case also W3 is related to W1 (acc. Martini and Ericsson)
! Case 4 describes double-differential data from MiniBooNE for neutrino
! and antineutrino scattering. It also describes the dd inclusive Xsection for
! neutrinos from T2K.
!******************************************************************************
module lepton2p2h

  use particleDefinition
  use eN_eventDefinition
  use leptonicID
  use CALLSTACK, only: TRACEBACK

  implicit none

  private
  public :: lepton2p2h_DoQE
  public :: lepton2p2h_DoDelta

  !****************************************************************************
  !****g* lepton2p2h/ME_Version
  ! PURPOSE
  ! indicate the type of matrix element parametrisation
  !
  ! SOURCE
  integer, save :: ME_Version = 4
  !
  ! possible values:
  ! * 1: const ME_Norm_XX  ! const for CC  fitted to MiniBooNE is 1.8e-6
  ! * 2: constant transverse and decreasing with Enu
  ! * 3: "Dipole transverse" transverse,  fall with Q2 as 4-th power
  ! * 4: MEC from E. Christy (8/2015), with parametrization for longitudinal
  ! * 5: MEC from Bosted arXiV:1203.2262, with parametrization for longitudinal
  ! * 6: MEC additional parametrization, with parametrization for longitudinal
  ! *    not yet implemented
  !****************************************************************************


  ! The following are all tunable strength parameters for 2p2h hadronic
  ! structure functions. Default is no tuning, i.e. all parameters = 1.
  ! except for ME_Long, for which default is =0 (no longitudinal component)



  !****************************************************************************
  !****g* lepton2p2h/ME_Norm_QE
  ! PURPOSE
  ! Overall strength of 2p2h matrix element with 2N out
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_Norm_QE    = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 gives the coded strength
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/ME_Norm_Delta
  ! PURPOSE
  ! Overall strength of 2p2h matrix element with NDelta out
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_Norm_Delta = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 is a dummy value
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/ME_Mass_QE
  ! PURPOSE
  ! Cutoff-mass in some parametrizations of 2p2h matrix element for NN out
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_Mass_QE    = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 is a dummy value
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/ME_Mass_Delta
  ! PURPOSE
  ! Cutoff-mass in some parametrizations of matrix element for NDelta out
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_Mass_Delta = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 is a dummy value
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/ME_Transversity
  ! PURPOSE
  ! Parametrisation of structure functions
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_Transversity = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value = 1 turns structure function W2 on so that 2p2h is pure transverse
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/ME_LONG
  ! PURPOSE
  ! Parametrization of structure functions
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_LONG = (/0.0, 0.0, 0.0/)
  ! NOTES
  ! The value = 0 turns any longitudinal contribution to structure funct. W2 off
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/ME_W3
  ! PURPOSE
  ! Overall strength factor for structure function W3
  !
  ! only for (CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_W3 = (/0.,1.0, 1.0/)
  ! NOTES
  ! overall strength parameter for structure function W3
  !****************************************************************************

!******************************************************************************
  !****g* lepton2p2h/ME_ODW
  ! PURPOSE
  ! switch for choosing the connection between structure functions
  ! W1(electron) and W1(neutrino) and W3(neutrino)
  ! =1 for expressions from Martini et al
  ! =2 for expressions from O'Connell et al
  !
  ! only for CC,NC)
  ! SOURCE
  integer, save :: ME_ODW = 1
  ! NOTES
  ! * O'Connell et al: PR C6 (1972) 719
  ! * Martini et al: PR C80 (2009) 065501
  !****************************************************************************

!******************************************************************************
  !****g* lepton2p2h/inmedW
  ! PURPOSE
  ! Controls which inv mass W is used in parametrization of 2p2h W1
  !
  ! SOURCE
  integer, save :: inmedW = 2
  ! NOTES
  ! * inmedW = 1 : W = free, static inv. mass in 2p2h parametrization of W1
  ! * inmedW = 2 : W = inv mass for Fermi moving nucleons in potential
  ! * inmedW = 3 : W = inv mass for Fermi moving nucleons without potential
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/T
  ! PURPOSE
  ! Isospin of target nucleus, factor in Eq. (A10) in O'Connell, Donnelly,
  ! Walecka, PR C6 (1972) 719
  !
  ! SOURCE
  real, save :: T = 1
  ! NOTES
  ! This is isospin T of NN pair, for 2p2h processes
  !****************************************************************************

  logical, save:: initflag = .true.

contains

  !****************************************************************************
  !****is* lepton2p2h/readInput
  ! NAME
  ! subroutine readInput
  !****************************************************************************
  subroutine readInput

    use output

    integer :: ios, i

    !**************************************************************************
    !****n* lepton2p2h/Lepton2p2h
    ! NAME
    ! NAMELIST /Lepton2p2h/
    ! PURPOSE
    ! Includes parameters for 2p2h events:
    ! * ME_Version
    ! * ME_Norm_QE
    ! * ME_Norm_Delta
    ! * ME_Mass_QE
    ! * ME_Mass_Delta
    ! * ME_Transversity
    ! * ME_LONG
    ! * ME_W3
    ! * ME_ODW
    ! * inmedW
    ! * T
    !**************************************************************************
    NAMELIST /lepton2p2h/ ME_Version, &
                          ME_Norm_QE, ME_Norm_Delta, &
                          ME_Mass_QE, ME_Mass_Delta,ME_Transversity,ME_LONG, &
                          ME_W3,ME_ODW,inmedW,T

    if (.not.initFlag) return

    call Write_ReadingInput('lepton2p2h',0)
    rewind(5)
    read(5,nml=lepton2p2h,IOSTAT=ios)
    call Write_ReadingInput("lepton2p2h",0,ios)


    select case (ME_Version)

    case (1)
       write(*,'(A)') 'ME1  const =  4.0e-6 * ME_Norm_XX'
!  case 1 is model-I in "Lalakulich Gallmeister Mosel PRC86(2012)014614"

    case (2)
       write(*,'(A)') 'ME2  4.8e4 * 0.635 / Enu^2 * ME_Norm_XX  &
                      &  in transverse part only, decreasing with Enu'
!  case 2 is model-II from "Lalakulich Gallmeister Mosel PRC86(2012)014614"

    case (3)
       write(*,'(A)') 'ME3  dipole parameterization of the "form factor" &
                      & ME=8e4*(1+Q2/MA2)^{-4},transverse part only'
! case 3 gives a good description of MiniBooNE data with MA ~ 1.5 GeV

    case (4)
       write(*,'(A)') 'ME 4, parametrization of structure functions W1,W2,W3,&
                      & W1 from E. Christy'
       write(*,*) 'ME_NORM=', ME_NORM_QE,'ME_W3=',ME_W3,'ME_Trans=',    &
                   &   ME_Transversity, 'ME_LONG=',ME_LONG, 'inmedW=',inmedW, &
                   &   'T=',T

    case (5)
       write(*,'(A)') 'parametrization of structure functions W1,W2,W3,&
                      & W1 from Bosted'

    case (6)
       write(*,'(A)') 'parametrization of structure functions W1,W2,W3,&
                      & W1 new parametrization'

   case default
       write(*,*) 'ME_Version = ',ME_Version
       call TRACEBACK('wrong value for ME_Version')
    end select


   write(*,'(A)') 'parameters ( N N final state ):  [i=EM,CC,NC]'
   do i=1,3
      write(*,'("   A=",ES13.5,"  M=",ES13.5)') &
           & ME_Norm_QE(i),ME_Mass_QE(i),ME_Transversity(i),ME_LONG(i)
   end do
   write(*,'(A)') 'parameters ( N Delta final state ):'
   do i=1,3
      write(*,'("   A=",ES13.5,"  M=",ES13.5)') &
           & ME_Norm_Delta(i),ME_Mass_Delta(i)
   end do

    call Write_ReadingInput('lepton2p2h',1)

    initFlag=.false.

  end subroutine readInput

  !****************************************************************************
  !****s* lepton2p2h/lepton2p2h_DoQE
  ! NAME
  ! subroutine lepton2p2h_DoQE(eN,outPart,XS)
  !
  ! PURPOSE
  ! Do all the electron induced 2p2h-QE scattering gamma* N1 N2 -> N1' N2'
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- electron-Nucleon event info
  !
  ! OUTPUT
  ! * type(particle), dimension(:) :: OutPart -- the two produced nucleons
  ! * real :: XS -- the cross section
  !****************************************************************************
  subroutine lepton2p2h_DoQE(eN,outPart,XS)
!    use output, only: WriteParticle

    type(electronNucleon_event), intent(inout) :: eN
    type(particle),dimension(:), intent(inout) :: OutPart
    real, intent(out) :: XS

    logical :: flagOK

    if (initFlag) call readInput

    XS = 0.0

    call lepton2p2h_SelectN2(eN,flagOK)
    if (.not.flagOK)  return ! ==> failure

!    call write_electronNucleon_event(eN)

    call lepton2p2h_FinalState(eN,outPart,.true.,flagOK)

!!$    if (.not.flagOK) then
!!$       call write_electronNucleon_event(eN)
!!$       write(*,*) 'Failure'
!!$       stop
!!$    end if

    if (.not.flagOK) return ! ==> failure

    XS = lepton2p2h_XS(eN,outPart,.true.)
    outPart%perWeight=XS

!    if(XS < 0) write(*,*) 'XS= ',XS,'pW= ',outPart%perWeight

!    call WriteParticle(6,1,outPart)
!    stop

  end subroutine lepton2p2h_DoQE

  !****************************************************************************
  !****s* lepton2p2h/lepton2p2h_DoDelta
  ! NAME
  ! subroutine lepton2p2h_DoDelta(eN,outPart,XS)
  !
  ! PURPOSE
  ! Do all the electron induced 2p2h-QE scattering gamma* N1 N2 -> N Delta
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- electron-Nucleon event info
  !
  ! OUTPUT
  ! * type(particle), dimension(:) :: OutPart -- the two produced hadrons
  ! * real :: XS -- the cross section
  !****************************************************************************
  subroutine lepton2p2h_DoDelta(eN,outPart,XS)

    type(electronNucleon_event), intent(inout) :: eN
    type(particle),dimension(:), intent(inout) :: OutPart
    real, intent(out) :: XS

    logical :: flagOK

    if (initFlag) call readInput

    XS = 0.0

    call lepton2p2h_SelectN2(eN,flagOK)
    if (.not.flagOK) return ! ==> failure

    call lepton2p2h_FinalState(eN,outPart,.false.,flagOK)
    if (.not.flagOK) return ! ==> failure

    XS = lepton2p2h_XS(eN,outPart,.false.)
    outPart%perWeight=XS

  end subroutine lepton2p2h_DoDelta

    !**************************************************************************
  !****s* lepton2p2h/lepton2p2h_SelectN2
  ! NAME
  ! subroutine lepton2p2h_SelectN2(eN)
  !
  ! PURPOSE
  ! Finds the second nucleon for the 2p2h collision
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- electron-Nucleon event info
  !
  ! OUTPUT
  ! * type(electronNucleon_event) :: eN -- a second nucleon is added
  !
  ! NOTES
  ! * The seond particle is generated analytically, not by selecting
  !   a testparticle from the real particle vector.
  ! * This is at a very basic level. You may add more sophisticated features
  !   as eq. two-particle correlatione etc.
  ! * A threshold check  Wfree>(2*mN+1MeV) is performed
  !****************************************************************************
  subroutine lepton2p2h_SelectN2(eN,flagOK)

    use mediumDefinition
    use mediumModule, only: mediumAt
    use densitymodule, only: FermiMomAt
    use random, only: rn, rnOmega
    use constants, only: mN
    use energyCalc, only: energyDetermination
    use minkowski, only: abs4
    use lorentzTrafo, only: lorentz

    type(electronNucleon_event), intent(inout) :: eN
    logical, intent(out) :: flagOK

    type(medium)   :: media
    type(particle) :: partN2
    real :: p,pF
    type(particle), dimension(2) :: nucleon
    real, dimension(0:3) :: momentum
    integer :: i

    ! 0) Set some defaults:
    call setToDefault(partN2)
    flagOK = .false.

    partN2%ID = 1
    partN2%mass=mN

    ! 1) select charge:
    media=mediumAt(eN%nucleon%position)
    if (rn()*media%density.gt.media%densityProton) then
       partN2%charge = 0
    else
       partN2%charge = 1
    end if

    ! 2) select position:
    partN2%position = eN%nucleon%position

    ! 3) select 3-momentum:
    pF = FermiMomAt(partN2%position)
    p = pF * rn()**(1./3.)
    partN2%momentum(1:3) = p * rnOmega()
    partN2%momentum(0) = sqrt(mN**2+p**2)

    call energyDetermination(partN2)

    ! 4) change the eN information:

    eN%nucleon2 = partN2

    nucleon(1) = eN%nucleon
    nucleon(2) = eN%nucleon2
    momentum = eN%boson%momentum+nucleon(1)%momentum+nucleon(2)%momentum
    eN%betacm = momentum(1:3)/momentum(0)
    eN%W = abs4(momentum)

    ! we calculate Wfree in the CM system:
    do i=1,2
       nucleon(i)%position = 999999999
       call lorentz(eN%betacm,nucleon(i)%momentum)
       nucleon(i)%momentum(0) = FreeEnergy(nucleon(i))
       call lorentz(-eN%betacm,nucleon(i)%momentum)
    end do

    momentum = eN%boson%momentum+nucleon(1)%momentum+nucleon(2)%momentum
    eN%W_free = abs4(momentum)
! invariant mass of two-nucleon + boson system

    if (eN%W_free.le.2*mN+0.001) return ! ===> failure

    ! 5) abuse of 'offshellParameter' for storage of density,
    !    needed for the 'cross section' calculation in cases 1 - 3

    eN%nucleon2%offshellParameter = media%density
    flagOK = .true.

  end subroutine lepton2p2h_SelectN2



  !****************************************************************************
  !****s* lepton2p2h/lepton2p2h_FinalState
  ! NAME
  ! subroutine lepton2p2h_FinalState(eN,outPart,DoQE,flagOK)
  ! PURPOSE
  ! Generate the final state of the electron 2p2h event
  !****************************************************************************
  subroutine lepton2p2h_FinalState(eN,outPart,DoQE,flagOK)
    use mediumDefinition
    use mediumModule, only: mediumAt
    use collisionNumbering, only: pert_numbering
    use master_2Body, only: setKinematics
    use propagation, only: updateVelocity
    use IDtable, only: nucleon, delta
    use random, only: rn
    use particleProperties, only: hadron
    use baryonWidthMedium, only: get_MediumSwitch_coll

    type(electronNucleon_event), intent(in)    :: eN
    type(particle),dimension(:), intent(inout) :: OutPart
    logical,                     intent(in)    :: DoQE
    logical,                     intent(out)   :: flagOK

    type(particle), dimension(2) :: pairIN
    type(medium)         :: media
    real, dimension(1:3) :: betaToLRF
!     real, dimension(0:3) :: momentum
    integer :: i, ChargeIn

    flagOK = .false.

    if (size(OutPart).lt.2) call TRACEBACK('OutPart array too small.')

    pairIN = (/ eN%boson, eN%nucleon /)
    media=mediumAt(eN%nucleon%position)

    ChargeIn = eN%nucleon%Charge+eN%nucleon2%Charge

    call setToDefault(OutPart)

    if (DoQE) then !=== N N final state ===
       OutPart%ID =     (/ nucleon          , nucleon /)

       select case (eN%idProcess)
       case default ! == EM, NC, antiEM, antiNC
          OutPart%Charge = (/ eN%nucleon%Charge, eN%nucleon2%Charge /)

       case (2)      ! == CC
          select case (ChargeIn)
          case (0)
             OutPart%Charge = (/ 0, 1 /)
          case (1)
             OutPart%Charge = (/ 1, 1 /)
          case (2)
             return ! ==> failure
          case default
             call TRACEBACK('ChargeIn not allowed')
          end select

       case (-2)     ! == antiCC
          select case (ChargeIn)
          case (0)
             return ! ==> failure
          case (1)
             OutPart%Charge = (/ 0, 0 /)
          case (2)
             OutPart%Charge = (/ 0, 1 /)
          case default
             call TRACEBACK('ChargeIn not allowed')
          end select

       end select

    else           !=== N Delta final state ===
       OutPart%ID =     (/ nucleon          , delta /)

       select case (eN%idProcess)
       case default ! == EM, NC, antiEM, antiNC
          OutPart(1)%Charge = nint(rn())
          OutPart(2)%Charge = eN%nucleon%Charge+eN%nucleon2%Charge &
                            &  - OutPart(1)%Charge
       case (2)      ! == CC
          select case (ChargeIn)
          case (0)
             OutPart(1)%Charge = nint(rn())
             OutPart(2)%Charge = 1 - OutPart(1)%Charge
          case (1)
             OutPart(1)%Charge = 1
             if (rn()<1./3.) OutPart(1)%Charge = 0  ! Delta++ n : Delta+ p = 1:2
             !based on counting diagrams, but better ideas are needed
             OutPart(2)%Charge = 2 - OutPart(1)%Charge
          case (2)
             OutPart%Charge = (/ 1, 2 /)
          case default
             call TRACEBACK('ChargeIn not allowed')
          end select
!          call TRACEBACK('CC not yet implemented')
       case (-2)     ! == antiCC
          select case (ChargeIn)
          case (2)
             OutPart(1)%Charge = nint(rn())
             OutPart(2)%Charge = 1 - OutPart(1)%Charge
          case (1)
             OutPart(1)%Charge = 0
             if (rn()<3./4.) OutPart(1)%Charge = 1  ! Delta0 n : Delta- p = 4:3
             !based on counting diagrams, but better ideas are needed
             OutPart(2)%Charge =  - OutPart(1)%Charge
          case (0)
             OutPart%Charge = (/ 0, -1 /)
          case default
             call TRACEBACK('ChargeIn not allowed')
          end select
!          call TRACEBACK('antiCC not yet implemented')
       end select

       ! The following is in order to avoid problems in massass:
       if (.not.get_MediumSwitch_coll()) then
          ! minimal value: 0.938 + Delta-MinMass + epsilon
          if (eN%W_free .lt. hadron(1)%mass+hadron(2)%minmass+0.005) then
             return ! ==> failure
          end if
       end if


    end if

    OutPart%antiparticle=.false.
    OutPart%perturbative=.true.

    OutPart%perWeight=0. ! perturbative weight = XS (here only dummy)

    do i=1,2
       OutPart(i)%position=eN%nucleon%position
       OutPart(i)%event=pert_numbering(eN%nucleon)
    end do

! Now the final-state two nucleons get started with proper momentum and energy

    betaToLRF=0.
!     momentum = eN%boson%momentum+eN%nucleon%momentum+eN%nucleon2%momentum
!
!   setKinematics sets the final state for the two final nucleons in 2p2h
!

    call setKinematics (eN%W, eN%W_free, betaToLRF, eN%betacm, media, pairIn, &
                       & OutPart(1:2), flagOK, .false.)

    if (.not.flagOK) return ! ==> failure

    call updateVelocity(OutPart)

  end subroutine lepton2p2h_FinalState



  !****************************************************************************
  !****f* lepton2p2h/lepton2p2h_XS
  ! NAME
  ! real function lepton2p2h_XS(eN,outPart,DoQE)
  ! PURPOSE
  ! calculate the electron induced 2p2h-QE cross section
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- electron-Nucleon event info
  ! * type(particle),dimension(:) :: OutPart -- the outgoing particles
  ! * logical :: DoQE -- .true. for NN final state, .false. for N Delta
  ! OUTPUT
  ! * the function value
  ! NOTES
  ! * One has to give a realistic parametrization of the matrix element
  ! * If one randomly selects the position of the second particle, one
  !   has to respect this in the XS calculation (and  not
  !   multiply it with the density at the position)
  !****************************************************************************
  real function lepton2p2h_XS(eN,outPart,DoQE)

    use minkowski, only: abs4Sq
    use constants, only: pi,twopi,mN,hbarc,alphaQED,GF,coscab
    use twoBodyTools, only: pCM_sqr
    use nucleus, only: getTarget
    use nucleusdefinition

    type(tnucleus), pointer :: targetNuc
    type(electronNucleon_event), intent(in) :: eN
    type(particle),dimension(:), intent(in) :: OutPart
    logical,                     intent(in) :: DoQE

    real :: mf1_2,mf2_2  !squares of final state nucleon masses in 2p2h process
    real :: k,k1 ! absolute values of the 3-momentum of the in/outgoing leptons
    real :: sqpp,pcm2 ! sqpp = (q + p + p2)^2, pcm2 = (p_cm)^2
!    real :: d2PS      !d2PS = two-body phase space, includes factor (2pi)^4
    real :: ME        ! contraction of lepton and hadron tensor
    real :: Atarget
! NOTE: lepton tensor contains coupling constants and propagator**2
    real :: couplProp    !coupling constant times propagator^2
    real :: Q2 !,omega
    integer :: iP

    targetNuc => getTarget()
    Atarget = targetNuc%mass

    if (Atarget < 2) then
    lepton2p2h_XS = 0.
    return
    end if


    Q2=eN%QSquared                          !Q^2
!    omega=eN%boson%momentum(0)              !omega = energy transfer

    mf1_2 = abs4Sq(outPart(1)%momentum)
    mf2_2 = abs4Sq(outPart(2)%momentum)
    k1=absMom(en%lepton_out)
!    k =absMom(en%lepton_in)

!!!!!!!!!!!!!!!!
if (ME_version < 4) then
!!!!!!!!!!!!!!!!
    lepton2p2h_XS = en%lepton_out%momentum(0)/en%lepton_in%momentum(0)
    ! correct for both e and nu

    sqpp = abs4Sq(eN%boson%momentum+eN%nucleon%momentum+eN%nucleon2%momentum)
!   initial state Mandelstam s of boson + 2 nucleons
    pcm2=pCM_sqr(sqpp, mf1_2, mf2_2)
!   final state nucleon momentum squared in 2N cm system

!   dOmega = 4. * pi      ! for isotropic phase space of 2 outgoing nucleons
!   d2PS = 1.0/(4.*pi)**2  * sqrt(pcm2/sqpp) * dOmega
!   includes factor (2pi)^4, different from PDG

!    d2PS = 1.0/(4.*pi)  * sqrt(pcm2/sqpp)
!   2-body final state phase space in cm system, integrated over NN angle

    lepton2p2h_XS=lepton2p2h_XS &
         * sqrt(pcm2/(16.0*sqpp)) & ! <-- the deltas
         * 2*twopi               ! <-- the angular integration

! Now XS times probability for 2nd nucleon to be at same position
    lepton2p2h_XS=lepton2p2h_XS &
         * eN%nucleon2%offshellParameter ! <-- abuse !! = media%density

 else       ! for Bosted parametrization of MEC in W1


    iP = abs(eN%idProcess)
    select case (iP)

      case (1)
         couplProp = (4*alphaQED**2)/Q2**2 * en%lepton_out%momentum(0)**2
      case (2)
         couplProp = (GF*coscab)**2/(2*pi**2) * en%lepton_out%momentum(0)**2
      case (3)
         couplProp = 0.5*GF**2/(2*pi**2) * en%lepton_out%momentum(0)**2
         ! Coupling for NC very roughly approximated by factor 1/2
      case default
         write(*,*) 'reaction type must be EM, CC or NC'
         stop

    end select


!!!!!!!!!!!!!!!
end if
!!!!!!!!!!!!!!

! Now we have to calculate the Matrixelement:
    select case (ME_Version)

    case (1)
       ME = ME_const(eN)

    case (2)
       ME = ME_transverse(eN)*0.635/( eN%lepton_in%momentum(0) )**2
          !!!*1.17/eN%lepton_in%momentum(0)/(0.42+eN%lepton_in%momentum(0))**2

    case (3)
       ME = ME_Dipole_transverse(eN)

    case (4:)
       ME = ME_W1W2W3(eN)
!    for this case ME describes experiment (except for coupling strength)
!    and thus contains phase-space of outgoing particles

    case default
       write(*,*) 'Other cases ME_Version > 4 not yet implemented'
    end select
!!!!!!!!!!!!!!
if (ME_version < 4) then
      lepton2p2h_XS=lepton2p2h_XS* ME / (2. *eN%nucleon2%momentum(0))  &
            & * 1/(twopi**5 *8. *eN%nucleon%momentum(0))

      ! 1/GeV**2=1/1000**2/MeV**2=1/1000**2/(1/197 fm**2)=(197/1000)**2 fm**2
      !      = (197/1000)**2 * 10 mb
      ! Now the cross section is given in units of mb/GeV:

      lepton2p2h_XS = lepton2p2h_XS*hbarc**2*10.

      ! Symmetry-Factor:
      if (IsSamePart(OutPart(1),OutPart(2))) lepton2p2h_XS = lepton2p2h_XS *0.5
      if (IsSamePart(eN%nucleon,eN%nucleon2)) lepton2p2h_XS = lepton2p2h_XS *0.5

    else        ! parametrization of W1, W2, W3

     if (ME_version == 4) ME = ME * (0.145 - 0.174*Atarget**(-1./3.))/ &
     &    (0.145 - 0.174*12**(-1./3.)) * Atarget/12
!    The A-dependence here is taken from Mosel, Gallmeister, arXiv 16.06499,
!    normalized to C12, since Christy-Bosted fit was for this nucleus

     lepton2p2h_XS=couplProp * ME/Atarget
!    The scaling of ME with 1/Atarget is necessary
!    since all other cross sections are given in 1/Atarget.



    ! 1/GeV**2=1/1000**2/MeV**2=1/1000**2/(1/197 fm**2)=(197/1000)**2 fm**2
    ! = (197/1000)**2 * 10 mb
    ! Now the cross section is given in units of mb/GeV:

    lepton2p2h_XS = lepton2p2h_XS*hbarc**2*10.
!    write(*,*) 'lepton2p2h_XS = ',lepton2p2h_XS


!   No symmetry factors since the Bosted-Christy
!   X-sections are fitted to exp.

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 contains

    real function ME_const(eN)

      type(electronNucleon_event), intent(in) :: eN
      integer :: iP

      iP = abs(eN%idProcess) ! anti-... same as EM, NC and CC

      if (DoQE) then !=== N N final state ===
         ME_const=1.0e-5*ME_Norm_QE(iP)
      else           !=== N Delta final state ===
         ME_const=1.0e-5*ME_Norm_Delta(iP)
      end if

    end function ME_const


     real function ME_transverse(eN)
      use minkowski, only: metricTensor, Contract
      use leptonTensor

      type(electronNucleon_event), intent(in) :: eN
      complex, dimension(0:3,0:3) :: leptonTens, hadronTens, dummy
      integer :: mu,nu
      real :: ME
      integer :: iP

      iP = abs(eN%idProcess) ! anti-... same as EM, NC and CC

      if (DoQE) then !=== N N final state ===
         ME=4.8e4*ME_Norm_QE(iP)
      else           !=== N Delta final state ===
         ME=4.8e4*ME_Norm_Delta(iP)
      end if

      do mu=0,3
        do nu=0,3
           dummy(mu,nu)=eN%boson%momentum(mu)*eN%boson%momentum(nu)/eN%QSquared
        end do
      end do

      hadronTens = ME *( - metricTensor - dummy )
      leptonTens = leptonicTensor(eN%idProcess,eN%lepton_in%momentum,  &
                                 & eN%lepton_out%momentum)

      ME_transverse = Contract(hadronTens,leptonTens)

    end function ME_transverse



    !**************************************************************************
    !****f* lepton2p2h_XS/ME_Dipole_transverse
    ! NAME
    ! real function ME_Dipole_transverse(eN)
    ! PURPOSE
    ! calculate the 2p2h matrix element according to
    ! W_1(g_munu -q_um q_nu /Q2) * L^munu
    ! so that the contribution is only to the transverse part
    ! NOTES
    !
    ! You have full access to all incoming and outgoing particles:
    ! * eN%lepton_in  -- incoming lepton
    ! * eN%nucleon    -- incoming nucleon 1
    ! * eN%nucleon2   -- incoming nucleon 2
    !
    ! exchanged boson:
    ! * eN%boson      -- exchanged boson
    !
    ! even without considering the final state particles, you know the kind
    ! of process via 'eN%idProcess', which may take the values EM,NC,CC and
    ! also antiEM,antiNC,antiCC
    !**************************************************************************
    real function ME_Dipole_transverse(eN)
      use minkowski, only: metricTensor, Contract
      use leptonTensor

      type(electronNucleon_event), intent(in) :: eN
      complex, dimension(0:3,0:3) :: leptonTens, hadronTens, dummy
      integer :: mu,nu
      real :: ME
      integer :: iP

      iP = abs(eN%idProcess) ! anti-... same as EM, NC and CC

      if (DoQE) then !=== N N final state ===
         ME=8.0e4*ME_Norm_QE(iP)*(1. + eN%QSquared/ME_Mass_QE(iP)**2)**(-4)
      else           !=== N Delta final state ===
         ME=8.0e4*ME_Norm_Delta(iP)*(1+eN%QSquared/ME_Mass_Delta(iP)**2)**(-4)
      end if


      do mu=0,3
        do nu=0,3
           dummy(mu,nu)=eN%boson%momentum(mu)*eN%boson%momentum(nu)/eN%QSquared
        end do
      end do

      hadronTens = ME *( - metricTensor - dummy )
      leptonTens = leptonicTensor(eN%idProcess,eN%lepton_in%momentum,  &
                                 & eN%lepton_out%momentum)

      ME_Dipole_transverse = Contract(hadronTens,leptonTens)

    end function ME_Dipole_transverse



!******************************************************************************
    !****f* lepton2p2h_XS/ME_W1W2W3
    ! NAME
    ! real function ME_W1W2W3(eN)
    ! PURPOSE
    ! to calculate the 2p2h contribution to the inclusive cross sections for
    ! electrons and neutrinos, cross section depends on all 3 structure functs
    !
    ! NOTES
    !
    ! You have full access to all incoming and outgoing particles:
    ! * eN%lepton_in  -- incoming lepton
    ! * eN%nucleon    -- incoming nucleon 1
    ! * eN%nucleon2   -- incoming nucleon 2
    !
    ! exchanged boson:
    ! * eN%boson      -- exchanged boson
    !
    ! even without considering the final state particles, you know the kind
    ! of process via 'eN%idProcess', which may take the values EM,NC,CC and
    ! also antiEM,antiNC,antiCC
    !**************************************************************************
    real function ME_W1W2W3(eN)

    use particleDefinition
    use eN_eventDefinition
    use leptonicID

    type(electronNucleon_event), intent(in) :: eN

  ! real :: W1N,W2N,W3N
    real :: nuswitch= 0 ! switch for neutrino/antineutrino in structure function
  ! nuswitch = 0 for em, = +1 for neutrino, -1 for antineutrino
    real :: sinsqthetahalf,cossqthetahalf
    real :: omega, Q2   ! energy transfer, four-momentum transfer
    integer :: IP
    real :: GM0,GA0,MV,MA,GM,GA,GM2,GA2

    Q2=eN%QSquared                          !Q^2
    omega=eN%boson%momentum(0)              !omega = energy transfer

    !    vector and axial coupling constants and cutoff masses
     GM0 = 4.71
     GA0 = - 1.255
     MV =  0.84
     MA = 1.032

!   vector and axial coupling formfactors

    GM = GM0 * 1./(1 + Q2/MV**2)**2
    GA = GA0 * 1./(1 + Q2/MA**2)**2
    GM2 = GM**2
    GA2 = GA**2

    IP = abs(eN%IdProcess)
    if (IP==1) then
      nuswitch = 0
      GA2=0
    else
      nuswitch = sign(1,eN%IdProcess)
    end if

    sinsqthetahalf = 0.5*(1. - en%lepton_out%momentum(3)/k1)
    cossqthetahalf = 1 - sinsqthetahalf

    ME_W1W2W3 =  sinsqthetahalf * 2*W1(Q2,omega,GM2,GA2)                      &
       &  + cossqthetahalf * W2(Q2,omega,GM2,GA2)                             &
       &  - nuswitch*(en%lepton_in%momentum(0)+en%lepton_out%momentum(0))/mN  &
       &  * W3(Q2,omega,GM,GA) * sinsqthetahalf

    ME_W1W2W3 = ME_W1W2W3 * ME_Norm_QE(IP)

    if (IP /= 1) ME_W1W2W3 = ME_W1W2W3 * (T + 1)
  ! Factor 'T+1' here is isospin factor for neutrinos, not electrons
  ! in O'Connell, Donnelly and Walecka, PRC 6 (1972) 719, eq. A10

!!!!!!!!!!!!!!!!!!!!!!!!!
!TESTPRINTOUT
!W1N = W1(Q2,omega)
!W2N = W2(Q2,omega)
!W3N = W3(Q2,omega)
!write (*,*) 'IP =',IP,'nuswitch =',nuswitch
!write (*,*)'Enu =',en%lepton_in%momentum(0),'Q2 =', Q2,'omega =',omega
!write(*,*) 'W1=',W1N,'W2=',W2N,'W3 =',W3N,'ME =',ME_W1W2W3
!write (*,*) 'sinsqthetahalf=',sinsqthetahalf,'cossqthetahalf=',cossqthetahalf
!write(*,*) '----------------------------------------------------'
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!

! Enforce positivity constraint on structure functions by keeping matrixelement
! positive

     if (ME_W1W2W3 < 0) then
        write(*,*) 'enforce positivity constraint in 2p2h, set to 0:',ME_W1W2W3
        write(*,*) 'Q2 =',Q2,'omega =', omega
        ME_W1W2W3 = 0.
     end if

    end function ME_W1W2W3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real function W1(Q2,omega,GM2,GA2)
! Structure function W1

    real, intent(in) :: Q2,omega,GM2,GA2

      IP = abs(eN%idProcess)

 !     write (*,*) 'IP in W1 = ',IP

    select case (IP)

    case (1)                                    ! electron
      if (DoQE) then !=== N N final state ===
         W1= W1E(Q2,omega)
      else           !=== N Delta final state ===
         W1= W1E(Q2,omega)
      end if

    case (2:)                                    ! CC and NC

      if (DoQE) then !=== N N final state ===
         W1= W1NU(Q2,omega,GM2,GA2)
      else           !=== N Delta final state ===
         W1= W1NU(Q2,omega,GM2,GA2)
      end if

     case default
      write(*,*) 'ProcessID Error in 2p2h'
     end select

    end function W1

!!!!!!!!!!!!!!!!!!!!!!!!!!

    real function W1E(Q2,omega)
! Structure function for electrons, parametrization for MEC term only

        use constants, only: mN
        use nucleus, only: getTarget
        use nucleusdefinition
        use minkowski, only: abs3,abs4sq

        type(tnucleus), pointer :: targetNuc
!         type(particle), dimension(2) :: nucleon

        real ::a1,b1,c1,t1,dw2,Atarget,Wrecsq
        real, intent(in) :: Q2,omega
        real :: ENfree,pNplusQ
        real :: p0,p1,p2,p3,p4,p5,p18,p19,f

        targetNuc => getTarget()
        Atarget = targetNuc%mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Now Fermi-smearing, switch controls method
      select case (inmedW)
        case (1)
           Wrecsq = mN**2 + 2*mN*omega - Q2

        case (2)
           Wrecsq = abs4sq(eN%boson%momentum+eN%nucleon%momentum)

        case (3)
           ENfree = sqrt(abs3(eN%nucleon%momentum)**2  +mN**2)
           pNplusQ = abs3(eN%nucleon%momentum + eN%boson%momentum)
           Wrecsq = (omega + ENfree)**2 - pNplusQ**2

        case default
           write(*,*) 'wrong case for Wrecsq in 2p2h'

      end select

      if (Wrecsq <= 0.0) then
         W1E = 0.0
         return
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    select case (ME_Version)

    case (4)

!   This case returns the value for the structure function W1(MEC)
!   fitted by E.Christy to inclusive electron scattering data,
!   E. Christy, priv. comm., August 2015, fitted for C12

        a1 = 6.049*Q2**2 * exp(-1.0*Q2/1.298)/(0.314+Q2)**5.708
        b1 = 0.791 + 0.154*Q2
        c1 = 0.290
        t1 = (Wrecsq - b1)**2/c1**2/2.
        dw2 = Wrecsq + Q2 * (1. - 1./2.2) - 1.0*mN*mN
        if (dw2 < 0.0) dw2 = 0.

        if (DoQE) then !=== N N final state ===
           W1E = a1 * (exp(-1.*t1)*sqrt(dw2))/mN
        else           !=== N Delta final state ===
           W1E = a1 * (exp(-1.*t1)*sqrt(dw2))/mN
!       W1E = structure function W1 for electrons
        end if
!         write (*,*) 'W1_MEC=',W1E
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case (5)
!   This case returns the value for the structure function W1(MEC)
!   fitted by Mamyan and Bosted to inclusive electron scattering data,
!   http://arxiv.org/abs/arXiv:1203.2262, 2012


!      The following parameter-constants are from Bosted arXiV:1203.2262
          p0=0.005138
          p1=0.980710
          p2=0.046379
          p3=1.643300
          p4=6.982600
          p5=-0.226550
          p19=-0.045536
          if (4 < Atarget .and. Atarget < 21) p18 = 215
          if (20 < Atarget .and. Atarget < 51) p18 = 235
          if (50 < Atarget) p18 = 230

          f=(1. + max(0.3,Q2)/p3 )**p4 /( omega**p5          &
            &     * (1.+p18*Atarget**(1.+p19*Q2/2./mN/omega)))
!     In the Mamyan-Bosted paper eq. (10) is wrong, corrected here,
!     following Bosted code

        if (DoQE) then !=== N N final state ===
           W1E = p0*exp( -(Wrecsq-p1)**2/p2 )/f/mN
!        write (*,*) 'MEC=',MEC, '    Q2=',Q2,'   omega=',omega
      else           ! === N Delta final state ===
           W1E = p0*exp( -(Wrecsq-p1)**2/p2 )/f /mN
      end if

     case default
         write(*,*) 'ME_Version does not exist: ',ME_Version

     end select

     end function W1E

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     real function W1NU(Q2,omega,GM2,GA2)
! structure function W1 for neutrino-induced CC and NC MEC process

     use constants, only: mN

     real, intent(in) :: Q2,omega,GM2,GA2
     real :: qvec2

     qvec2 = Q2 + omega**2

  select case (ME_ODW)

   case (1)       ! case for connection between W1 and W3 from Martini et al

    W1NU= (1. + qvec2/omega**2 * GA2/GM2)* 2.0 * W1E(Q2,omega)

   case (2)        ! case for connection between W1 and W3 from O'Connell,
                  ! Donnelly, Walecka PR C6 (1972) 719
    W1NU = (1. + (2*mN)**2/qvec2 * GA2/GM2) * 2.0 * W1E(Q2,omega)

   case default
    write(*,*) 'MEW1W3 error in W1NU'
   end select

    end function W1NU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real function W2(Q2,omega,GM2,GA2)
! Structure function W2

    use constants, only: mN

    real, intent(in) :: Q2,omega,GM2,GA2
    real, parameter :: MDelta = 1.232
    real :: qvec2
    integer :: IP

    IP = abs(eN%idProcess)

    qvec2 = Q2 + omega**2
!   vector and axial coupling constants and cutoff masses

    W2 = ME_Transversity(IP) * Q2/qvec2 * W1(Q2,omega,GM2,GA2)
! W2: term necessary for purely transverse interaction, could be turned off
!     by setting ME_Transversity = 0, default = 1

    if (ME_Long(IP) > 0) &
    W2 = W2  + GA2*(MDelta - mN)**2/(2.*(Q2 + omega**2))* 1./(1 + Q2/0.3**2)**2&
         & * ME_LONG(iP) * 1.e-5
!   Structure of longitudinal term follows Martini et al (PRC 2009)
!   ME_LONG: allows to turn off longitudinal contribution to 2p2h, default = 0
!   W2: 2nd term for longitudinal response, strength function untested
    end function W2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real function W3(Q2,omega,GM,GA)
 ! Structure function W3, relevant only for neutrinos
 ! W3 is directly related to W1, either according to Martini and Ericsson,
 ! or to O'Connell, Donnelly, Walecka

    use constants, only: mN

    real, intent(in) :: Q2,omega,GM,GA
    real :: qvec2

    integer IP

    IP = abs(eN%idProcess)
    if (IP==1) then
      W3 = 0.0
      return
    end if

    qvec2 = Q2 + omega**2

 select case (ME_ODW)

   case (1)       ! case for connection between W1 and W3 from Martini et al

    W3 = 4.0 * qvec2/omega**2 * GA/GM * W1E(Q2,omega)* ME_W3(IP)
!   ME_W3(IP) = arbitrary possible strength factor for W3, default = 1

   case (2)        ! case for connection between W1 and W3 from O'Connell,
                  ! Donnelly, Walecka PR C6 (1972) 719
    W3 = 4. * (2.*mN)**2/qvec2 * GA/GM * W1E(Q2,omega) * ME_W3(IP)

   case default
    write(*,*) 'MEW1W3 error in MECNU'
   end select

    end function W3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end function lepton2p2h_XS





end module lepton2p2h
