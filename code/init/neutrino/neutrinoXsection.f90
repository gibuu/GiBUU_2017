!******************************************************************************
!****m* /neutrinoXsection
! NAME
! module neutrinoXsection
!
! PURPOSE
! This module calculates various neutrino cross sections, depending on the
! choice of the user:
! * integrated or differential or double differential cross sections
!   (according to the value of nuXsectionmode, the corresponding
!   subroutines are called in the neutrino init file)
! * EM, CC or NC (according to the value of process_ID given as input to the
!   corresponding subroutines)
! * muon, electron or tau flavor (according to the value of flavor_ID given
!   as input to the corresponding subroutines)
!******************************************************************************

module neutrinoXsection

  use particleDefinition
  use eN_eventDefinition

  implicit none
  private

  !****************************************************************************
  !****g* neutrinoXsection/debugflag
  ! SOURCE
  logical, parameter :: debugflag=.false.
  ! PURPOSE
  ! to switch on/off debug information
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/nuclear_phasespace
  ! SOURCE
  logical, save ::  nuclear_phasespace=.true.
  ! PURPOSE
  ! to change between different phasespace factors.
  ! (change only for debugging purposes):
  ! * .false. = 1/SP(k_in,p_in)  i.e. phasespace factor for each nucleon
  ! * .true.  = 1/k_in(0)*p_in(0) i.e. global phasespace factor for the nucleus
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/singlePiModel
  ! SOURCE
  integer, save :: singlePiModel=1
  ! PURPOSE
  ! to change between different models for the pion nucleon cross section:
  ! * 0 = pi N according to Nieves et al (hep-ph/0701149)
  ! * 1 = MAID-like model
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoXsection/integralPrecision
  ! SOURCE
  integer, save ::  integralPrecision=3
  ! PURPOSE
  ! precision for the Gauss integration
  ! (reduce it for nuXsectionMode.eq.0 (sigma) to e.g. 2)
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/integralPrecisionQE
  ! SOURCE
  integer, save ::  integralPrecisionQE=500
  ! PURPOSE
  ! precision for the Gauss integration over the QE peak
  ! (reduce it for nuXsectionMode.eq.0 (sigma) to e.g. 300)
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/enu
  ! SOURCE
  real, save :: enu=-10.
  ! PURPOSE
  ! neutrino energy, read in by namelist
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/delta_enu
  ! SOURCE
  real, save :: delta_enu=-10.
  ! PURPOSE
  ! value by which the neutrino energy is increased
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/Qs
  ! SOURCE
  real, save :: Qs=-10.
  ! PURPOSE
  ! momentum transfer squared
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/delta_Qs
  ! SOURCE
  real, save :: delta_Qs=-10.
  ! PURPOSE
  ! value by which Qs is increased
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/W
  ! SOURCE
  real, save :: W=-10.
  ! PURPOSE
  ! invariant mass defined as (p+q)^2
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/delta_W
  ! SOURCE
  real, save :: delta_W=-10.
  ! PURPOSE
  ! value by which W is increased
  !****************************************************************************



  !****************************************************************************
  !****g* neutrinoXsection/costheta
  ! SOURCE
  real, save :: costheta=-10.
  ! PURPOSE
  ! cosine of the angle between the neutrino (z-direction) and the
  ! outgoing lepton
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/delta_costheta
  ! SOURCE
  real, save :: delta_costheta=-10.
  ! PURPOSE
  ! value by which costheta is increased
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/elepton
  ! SOURCE
  real, save :: elepton=-10.
  ! PURPOSE
  ! energy of the outgoing lepton
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/delta_elepton
  ! SOURCE
  real, save :: delta_elepton=-10.
  ! PURPOSE
  ! value by which elepton is increased
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/invariantMassCut
  ! SOURCE
  real, save :: invariantMassCut=100.
  ! PURPOSE
  ! cut events with invariant Mass above this value (in GeV);
  ! cut pion production from Delta and DIS on Wrec = Sqrt[M^2 + 2*M*nu - Q^2]
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/invariantMassCut_BG
  ! SOURCE
  real, save :: invariantMassCut_BG=100.
  ! PURPOSE
  ! cut MAID-like background events with invariantMass_BG above this value
  ! (in GeV);
  ! cut 1pi BG on Wrec = Sqrt[M^2 + 2*M*nu - Q^2]
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/DISmassless
  ! SOURCE
  logical, save :: DISmassless = .false.
  ! PURPOSE
  ! Flag whether the PYTHIA calculation for the DIS process is done with
  ! massless or massive (Mquark~300MeV) quarks and diquarks.
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/DIScutW1
  ! SOURCE
  real, save :: DIScutW1 = 1.6
  ! PURPOSE
  ! lower W-cut for linear onset of DIS
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/DIScutW2
  ! SOURCE
  real, save :: DIScutW2 = 1.65
  ! PURPOSE
  ! upper W-cut for linear onset of DIS
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/REScutW1
  ! SOURCE
  real, save :: REScutW1 = 2.00
  ! PURPOSE
  ! lower W-cut for linear turn-off of resonances
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/REScutW2
  ! SOURCE
  real, save :: REScutW2 = 2.05
  ! PURPOSE
  ! upper W-cut for linear turn-off or resonances
  !****************************************************************************
  
  !****************************************************************************
  !****g* neutrinoXsection/mcutDIS
  ! SOURCE
  real, save :: mcutDIS = 0.71
  ! PURPOSE
  ! parameter to control Q^2 dependence of DIS
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/MC_xmax
  ! SOURCE
  real,save :: MC_xmax = 2.0
  ! PURPOSE
  !
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/DISformfakEM
  ! SOURCE
  integer,save :: DISformfakEM = 2
  ! PURPOSE
  ! Introduce an additional form factor for the DIS cross section, when
  ! processed via a photon:
  ! * 0: no form factor
  ! * 1: Q^2/(mcutDIS^2+Q^2)
  ! * 2: Q^4/(mcutDIS^2+Q^2)^2
  !
  ! In case of electron induced events, we need choose 2 in order to be
  ! compatible with Pythia's electron machinery.
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/DISformfakNCCC
  ! SOURCE
  integer,save :: DISformfakNCCC = 1
  ! PURPOSE
  ! Introduce an additional form factor for the DIS cross section, when
  ! processed via W or Z boson:
  ! * 0: no form factor
  ! * 1: Q^2/(mcutDIS^2+Q^2)
  ! * 2: Q^4/(mcutDIS^2+Q^2)^2
  !
  ! In case of electron induced events, we need choose 2 in order to be
  ! compatible with Pythia's electron machinery.
  !****************************************************************************



  real,save :: MC_x,MC_y
  real,save :: xyJacob  ! Jacobian from dx dy -> dcos(theta) dE'


  public :: MC_x,MC_y

  logical, save :: initflag=.true.

  real, parameter :: corr_factor=1.9732696817**2*10.**10 ! hbarc**2 * 10**12

  public :: Xsec_integratedSigma
  public :: Xsec_dSigmadCosThetadElepton
  public :: Xsec_dSigmadQsdElepton
  public :: Xsec_dSigmadCosTheta
  public :: Xsec_dSigmadElepton
  public :: Xsec_SigmaMC, Xsec_SigmaMC_Qs, Xsec_SigmaMC_W
  public :: SetXsecMC
  public :: XsecdCosthetadElepton  ! required in the public list for calculation of structure function F2

! for TEMPORARY use:

  public :: getEleptonLimits
  public :: SetHadronCharge
  public :: get_xsection_namelist


contains
  !****************************************************************************
  !****s* neutrinoXsection/readinput
  ! NAME
  ! subroutine readinput(nuXsectionmode)
  ! PURPOSE
  ! This subroutine reads input out of jobcard.
  ! Namelist 'nl_neutrinoxsection' is read in always.
  ! All other namelist are read depending on the value of
  ! nuXsectionmode.
  ! If a nuXsectionmode for a particular experiment is chosen, with inherent
  ! flux integration, the parameter enu must be set to some dummy, but
  ! reasonable value > 0
  !****************************************************************************
  subroutine readinput(nuXsectionmode, isExp)
    use output
    use neutrino_IDTable

    integer :: IOS
    integer, intent(in) :: nuXsectionmode
    logical, intent(in), optional :: isExp

    !**************************************************************************
    !****n* neutrinoXsection/nl_neutrinoxsection
    ! NAME
    ! NAMELIST /nl_neutrinoxsection/
    ! PURPOSE
    ! This Namelist includes:
    ! * integralPrecision
    ! * integralPrecisionQE
    ! * singlePiModel
    ! * invariantMassCut
    ! * invariantMassCut_BG
    ! * DISmassless
    ! * DIScutW1
    ! * DIScutW2
    ! * REScutW1
    ! * REScutW2
    ! * DISformfakEM
    ! * DISformfakNCCC
    ! * mcutDIS
    !**************************************************************************
    NAMELIST /nl_neutrinoxsection/ integralPrecision,&
         & integralPrecisionQE,singlePiModel,&
         & invariantMasscut,invariantMasscut_BG, &
         & DISmassless, DIScutW1, DIScutW2, REScutW1, REScutW2, &
         & DISformfakEM, DISformfakNCCC,mcutDIS

    !**************************************************************************
    !****n* neutrinoXsection/nl_integratedSigma
    ! NAME
    ! NAMELIST /nl_integratedSigma/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=integratedsigma and
    ! includes:
    ! * enu
    ! * delta_enu
    !**************************************************************************
    NAMELIST /nl_integratedSigma/ enu, delta_enu

    !**************************************************************************
    !****n* neutrinoXsection/nl_dSigmadCosThetadElepton
    ! NAME
    ! NAMELIST /nl_dSigmadCosThetadElepton/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=dSigmadCosThetadElepton
    ! and includes:
    ! * enu
    ! * costheta
    ! * elepton
    ! * delta_elepton
    !**************************************************************************
    NAMELIST /nl_dSigmadCosThetadElepton/ enu, costheta, elepton, delta_elepton

    !**************************************************************************
    !****n* neutrinoXsection/nl_dSigmadQsdElepton
    ! NAME
    ! NAMELIST /nl_dSigmadQsdElepton/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=dSigmadQsdElepton
    ! and includes:
    ! * enu
    ! * Qs
    ! * elepton
    ! * delta_elepton
    !**************************************************************************
    NAMELIST /nl_dSigmadQsdElepton/ enu, Qs, elepton, delta_elepton

    !**************************************************************************
    !****n* neutrinoXsection/nl_dSigmadQs
    ! NAME
    ! NAMELIST /nl_dSigmadQs/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=dSigmadQs and includes:
    ! * enu
    ! * Qs
    ! * delta_Qs
    !**************************************************************************
    NAMELIST /nl_dSigmadQs/ enu, Qs, delta_Qs

    !**************************************************************************
    !****n* neutrinoXsection/nl_dSigmadW
    ! NAME
    ! NAMELIST /nl_dSigmadW/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=dSigmadW and includes:
    ! * enu
    ! * W
    ! * delta_W
    !**************************************************************************
    NAMELIST /nl_dSigmadW/ enu, W, delta_W

    !**************************************************************************
    !****n* neutrinoXsection/nl_dSigmadcostheta
    ! NAME
    ! NAMELIST /nl_dSigmadcostheta/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=dSigmadcostheta and
    ! includes:
    ! * enu
    ! * costheta
    ! * delta_costheta
    !**************************************************************************
    NAMELIST /nl_dSigmadCosTheta/ enu, costheta, delta_costheta

    !**************************************************************************
    !****n* neutrinoXsection/nl_dSigmadElepton
    ! NAME
    ! NAMELIST /nl_dSigmadElepton/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=dSigmadElepton and includes:
    ! * enu
    ! * elepton
    ! * delta_elepton
    !**************************************************************************
    NAMELIST /nl_dSigmadElepton/ enu, elepton, delta_elepton

    !**************************************************************************
    !****n* neutrinoXsection/nl_SigmaMC
    ! NAME
    ! NAMELIST /nl_SigmaMC/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=dSigmaMC and
    ! includes:
    ! * enu
    ! * MC_xmax
    !**************************************************************************
    NAMELIST /nl_SigmaMC/ enu,MC_xmax

    call Write_ReadingInput('nl_neutrinoxsection',0)
    rewind(5)
    read(5,nml=nl_neutrinoxsection,IOSTAT=IOS)
    call Write_ReadingInput('nl_neutrinoxsection',0,IOS)

    write(*,*) 'integral precicions: ', integralPrecision, integralPrecisionQE

    select case (singlePiModel)
    case (0)
       write(*,*) 'pi N cross section (Delta + interfering background) according to Nieves et al'
    case (1)
       write(*,*) 'pi N background (resonances are subtracted) is taken MAID-like'
    case default
       write(*,*) 'no valid input for singlePiModel -> stop', singlePiModel
       stop
    end select

    write(*,'(a,F12.4)') ' cut         events with invariant masses above', &
         & invariantMassCut
    write(*,'(a,F12.4)') ' cut MAID BG events with invariant masses above', &
         & invariantMassCut_BG

    write(*,'(A,L2)') ' DIS with massless (di-)quarks ?',DISmassless
    write(*,'(A,2F8.3)') ' DIS W-cut = ',DIScutW1,DIScutW2
    select case (DISformfakEM)
    case (0)
       write(*,'(A)') ' DIS form factor (EM)   : 0,  = 1'
    case (1)
       write(*,'(A)') ' DIS form factor (EM)   : 1,  = Q^2/(mcutDIS^2+Q^2)'
    case (2)
       write(*,'(A)') ' DIS form factor (EM)   : 2,  = Q^4/(mcutDIS^2+Q^2)^2'
    case default
       write(*,*) 'value ',DISformfakEM,' not valid for DISformfakEM!'
       stop
    end select

    select case (DISformfakNCCC)
    case (0)
       write(*,'(A)') ' DIS form factor (CC,NC): 0,  = 1'
    case (1)
       write(*,'(A)') ' DIS form factor (CC,NC): 1,  = Q^2/(mcutDIS^2+Q^2)'
    case (2)
       write(*,'(A)') ' DIS form factor (CC,NC): 2,  = Q^4/(mcutDIS^2+Q^2)^2'
    case default
       write(*,*) 'value ',DISformfakNCCC,' not valid for DISformfakNCCC!'
       stop
    end select

    call Write_ReadingInput('nl_neutrinoxsection',1)

    select case (nuXsectionMode)

    case (integratedSigma)

       call Write_ReadingInput('nl_integratedSigma',0)
       rewind(5)
       read(5,nml=nl_integratedSigma,IOSTAT=IOS)
       call Write_ReadingInput('nl_integratedSigma',0,IOS)
       if (enu.lt.0..or.delta_enu.lt.0.) then
          write(*,*) 'input error in neutrinoXsection, enu= ', enu
          write(*,*) 'stop'
          stop
       end if
       call Write_ReadingInput('nl_integratedSigma',1)

    case (dSigmadCosThetadElepton)

       call Write_ReadingInput('nl_dSigmadCosThetadElepton',0)
       rewind(5)
       read(5,nml=nl_dSigmadCosThetadElepton,IOSTAT=IOS)
       call Write_ReadingInput('nl_dSigmadCosThetadElepton',0,IOS)
       if (enu.lt.0..or.abs(costheta).gt.1..or.elepton.lt.0..or.delta_elepton.lt.0.) then
          write(*,*) 'input error in neutrinoXsection'
          write(*,*) 'stop'
          stop
       end if
       call Write_ReadingInput('nl_dSigmadCosThetadElepton',1)

    case (dSigmadQsdElepton)

       call Write_ReadingInput('nl_dSigmadQsdElepton',0)
       rewind(5)
       read(5,nml=nl_dSigmadQsdElepton,IOSTAT=IOS)
       call Write_ReadingInput('nl_dSigmadQsdElepton',0,IOS)
       if (enu.lt.0..or.Qs.lt.0..or.elepton.lt.0..or.delta_elepton.lt.0.) then
          write(*,*) 'input error in neutrinoXsection'
          write(*,*) 'stop'
          stop
       end if
       call Write_ReadingInput('nl_dSigmadQsdElepton',1)

    case (dSigmadQs)

       call Write_ReadingInput('nl_dSigmadQs',0)
       rewind(5)
       read(5,nml=nl_dSigmadQs,IOSTAT=IOS)
       call Write_ReadingInput('nl_dSigmadQs',0,IOS)
       if (enu.lt.0. .or. Qs.lt.0. .or. delta_Qs.lt.0.) then
          write(*,*) 'input error in nl_dSigmadQs: enu or Qs or delta_Qs <0'
          write(*,*) 'stop'
          stop
       end if
       call Write_ReadingInput('nl_dSigmadQs',1)

    case (dSigmadCosTheta)

       call Write_ReadingInput('nl_dSigmadCosTheta',0)
       rewind(5)
       read(5,nml=nl_dSigmadCosTheta,IOSTAT=IOS)
       call Write_ReadingInput('nl_dSigmadCosTheta',0,IOS)
       if (enu.lt.0..or.abs(costheta).gt.1..or.delta_costheta.lt.0.) then
          write(*,*) 'input error in neutrinoXsection'
          write(*,*) 'stop'
          stop
       end if
       call Write_ReadingInput('nl_dSigmadCosTheta',1)

    case (dSigmadElepton)

       call Write_ReadingInput('nl_dSigmadElepton',0)
       rewind(5)
       read(5,nml=nl_dSigmadElepton,IOSTAT=IOS)
       call Write_ReadingInput('nl_dSigmadElepton',0,IOS)
       if (enu.lt.0..or.elepton.lt.0..or.delta_elepton.lt.0.) then
   !      write(*,*) 'enu=',enu,'elepton=',elepton,'delta_elepton=',delta_elepton
          write(*,*) 'input error in neutrinoXsection'
          write(*,*) 'stop'
          stop
       end if
       call Write_ReadingInput('nl_dSigmadElepton',1)

    case (dSigmaMC)

       if (present(isExp)) then
          if (isExp) then
             enu = 99.9 ! dummy
          end if
       end if

       call Write_ReadingInput('nl_SigmaMC',0)
       rewind(5)
       read(5,nml=nl_SigmaMC,IOSTAT=IOS)
       call Write_ReadingInput('nl_SigmaMC',0,IOS)
       if (enu.lt.0.) then
          write(*,*) 'input error in neutrinoXsection, enu= ', enu
          write(*,*) 'stop'
          stop
       end if

       write(*,'(" Enu      =",f12.3)') enu
       write(*,'(" xmax     =",f12.3)') MC_xmax

       call Write_ReadingInput('nl_SigmaMC',1)


    case (dSigmadW)

       call Write_ReadingInput('nl_dSigmadW',0)
       rewind(5)
       read(5,nml=nl_dSigmadW,IOSTAT=IOS)
       call Write_ReadingInput('nl_dSigmadW',0,IOS)
       if (enu.lt.0. .or. W.lt.0. .or. delta_W.lt.0.) then
          write(*,*) 'input error in nl_dSigmadW: enu or W or delta_W <0'
          write(*,*) 'stop'
          stop
       end if
       call Write_ReadingInput('nl_dSigmadW',1)



    case default
       write(*,*) 'error in case nuXsectionMode', nuXsectionMode
       stop

    end select

    write(*,*)
    write(*,*) 'nuXsectionMode =',sXsectionMode(nuXsectionMode)
    write(*,*)

  end subroutine readinput

  !****************************************************************************
  !****s* neutrinoXsection/Xsec_SigmaMC
  ! NAME
  ! subroutine Xsec_SigmaMC(eNev,IP,raiseFlag,raiseVal,OutPart,sig,flux_enu)
  !
  ! PURPOSE
  ! This subroutine calculates the integrated cross section depending
  ! on the input variables. It applies a MC integration technique.
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- the Lepton-Nucleon Event
  ! * integer             :: IP           -- ID of outgoing hadron/process
  ! * logical             :: raiseFlag    -- shall energy etc. be increased?
  ! * real, optional      :: flux_enu     -- if present, use this value as
  !   neutrino energy (used for experiments)
  !
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * real                         :: raiseVal -- the actual value of the
  !   running variable (i.e. elepton, enu, costheta, ...)
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !****************************************************************************
  subroutine Xsec_SigmaMC(eNev, IP, raiseFlag, raiseVal, &
       & OutPart,sig,flux_enu)

    use neutrino_IDTable

    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out

    if (initFlag) then
       call readInput(dSigmaMC, present(flux_enu))
       initFlag = .false.

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if

    end if

    raiseVal=eNev%lepton_in%momentum(0)

    ! set default output
    call setToDefault(OutPart)
    sig=0.


    if (xyJacob .eq. 0.0) return

    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    call XsecdCosthetadElepton(eNev,IP,OutPart, sig)
    ! Monte-Carlo integration over x gives * MC_xmax; for the integration over y the  the factor is 1
    ! xyJacob is Jacobian from x,y to cos(theta) dE'
    sig = sig * xyJacob * MC_xmax

  end subroutine Xsec_SigmaMC

  !****************************************************************************
  !****s* neutrinoXsection/SetXsecMC
  ! NAME
  ! subroutine SetXsecMC()
  ! PURPOSE
  ! set the values of the integration variables in the MC integration
  ! called from subroutine initNeutrino in file initNeutrino.f90
  !****************************************************************************
  subroutine SetXsecMC(eNev, flux_enu, XsectionMode)
    use random, only: rn
    use eN_event, only: eNev_init_nuStep2,eNev_init_nuStep3c
    use neutrino_IDTable
    use minkowski, only: SP, abs4Sq

    type(electronNucleon_event),  intent(inout) :: eNev
    real, intent(in) :: flux_enu
    integer, intent(in) :: XsectionMode

    real :: Ein, Eprime
    real :: PK, PP
    logical :: flagOK

    if (initFlag) then
       select case (MOD(XsectionMode,10))
       case (6)
         call readInput(dSigmaMC, (XsectionMode>10))
       case (3)
         call readInput(dSigmadQs)
       case (7)
         call readInput(dSigmadW)
       end select
       initFlag = .false.

       if (flux_enu.gt.0.0) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
    end if

    Ein = enu
    if (flux_enu.gt.0.0) Ein = flux_enu

    call eNev_init_nuStep2(eNev,Ein) ! we ignore threshold checks
    PK = SP(eNev%lepton_in%momentum,eNev%nucleon%momentum)
    PP = abs4Sq(eNev%nucleon%momentum)

    select case (MOD(XsectionMode,10))
    case (6)
      MC_x = rn()*MC_xmax
      MC_y = rn()
    case (3)
      MC_y=rn()
      MC_x=Qs/(2.*MC_y*PK)
    case (7)
      MC_y=rn()
      MC_x=1.-(W**2-PP)/(2.*MC_y*PK)
      !write(*,*) 'MC_y=',MC_y, '    MC_x=',MC_x
    end select

    call eNev_init_nuStep3c(eNev,MC_x,MC_y,flagOK)
    if (.not.flagOK) then
       xyJacob = 0.0
       !write(*,*) 'after eNev_init_nuStep3c flakOK is false. xyJacob = 0.0'
    else
       Eprime=eNev%lepton_out%momentum(0)
       ! Jacobian xyJacob from dsigma/dE1dcostheta to dsigma/dxdy:
       xyJacob = (MC_y*PK**2)/(Ein *sqrt(Eprime**2-eNev%lepton_out%mass**2)*eNev%nucleon%momentum(0))
       !write(*,*) 'xyJacob=',xyJacob
    end if

    ! MC integration is transferred to the subroutine  Xsec_SigmaMC

  end subroutine SetXsecMC





  subroutine Xsec_SigmaMC_Qs(eNev, IP, raiseFlag, raiseVal, OutPart,sig,flux_enu)

    use neutrino_IDTable

    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out


    if (initFlag) then
       call readInput(dSigmadQs)
       initFlag = .false.

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
    end if


    if (raiseFlag) then
       Qs=Qs+delta_Qs
       write(*,'(3(A,g12.5))') 'Enu=', enu, &
            & '      Qs is raised by ...', delta_Qs, '  to  Qs=', Qs
    end if
    raiseVal=Qs


    ! set default output
    call setToDefault(OutPart)
    sig=0.


    if (xyJacob .eq. 0.0) return

    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    call XsecdCosthetadElepton(eNev,IP,OutPart, sig)

    !factor /Qs *MC_x is   from dsi/x to dsi/dQ2
    sig = sig * xyJacob /Qs *MC_x

  end subroutine Xsec_SigmaMC_Qs




  subroutine Xsec_SigmaMC_W(eNev, IP, raiseFlag, raiseVal, OutPart,sig,flux_enu)

    use neutrino_IDTable
    use minkowski, only: SP


    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out
    real    :: PK

    if (initFlag) then
       call readInput(dSigmadW)
       initFlag = .false.

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
    end if


    if (raiseFlag) then
       W=W+delta_W
       write(*,'(3(A,g12.5))') 'Enu=', enu, &
            & '      W is raised by ...', delta_W, '  to  W=', W
    end if
    raiseVal=W


    ! set default output
    call setToDefault(OutPart)
    sig=0.


    if (xyJacob .eq. 0.0) return

    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    call XsecdCosthetadElepton(eNev,IP,OutPart, sig)
    !write(*,*) 'sig=',sig
    !factor W /(PK *MC_y)  is   from dsi/x to dsi/dW
    PK = SP(eNev%lepton_in%momentum,eNev%nucleon%momentum)
    !write(*,*) 'PK=',PK, '   W=',W, '   MC_y=',MC_y, '    xyJacob=',xyJacob
    sig = sig * xyJacob *W /(PK *MC_y)

  end subroutine Xsec_SigmaMC_W




















  !****************************************************************************
  !****s* neutrinoXsection/Xsec_integratedSigma
  ! NAME
  ! subroutine Xsec_integratedSigma(eNev,IP,raiseFlag,raiseVal,OutPart,sig,
  ! flux_enu)
  !
  ! PURPOSE
  ! This subroutine calculates the integrated cross section depending
  ! on the input variables:
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- the Lepton-Nucleon Event
  ! * integer             :: IP           -- ID of outgoing hadron/process
  ! * logical             :: raiseFlag    -- shall energy etc. be increased?
  ! * real, optional      :: flux_enu     -- if present, use this value as
  !   neutrino energy (used for experiments)
  !
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * real                         :: raiseVal -- the actual value of the
  !   running variable (i.e. elepton, enu, costheta, ...)
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !****************************************************************************
  subroutine Xsec_integratedSigma(eNev, IP, raiseFlag, raiseVal, &
       & OutPart,sig,flux_enu)

    use neutrino_IDTable
    use idtable, only: nucleon
    use random
    use gauss_integration
    use eN_event, only: eNev_init_nuStep2,eNev_init_nuStep3a
    !!! use lepton_kinematics_free, only: minmaxE1_costheta



    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out
    real :: costheta_max, costheta_min
    real :: elepton_max, elepton_min

    integer :: intprec,intprec1,l,j,n2,n1
    real :: eleptonint,costhetaint
    real :: sigmaximum,sigrd
    real, dimension(:),allocatable :: y,x,yy,xx
    logical :: flagOK


    if (initFlag) then
       call readInput(integratedSigma)
       initFlag = .false.

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
    end if


    ! set default output
    call setToDefault(OutPart)
    sig=0.

    if (raiseFlag.and..not.present(flux_enu)) then
       enu=enu+delta_enu
       write(*,'(a,F12.4,a,F12.4)') 'Enu is raised by ...', &
            & delta_enu, ' to ', enu
    end if

    if (present(flux_enu)) enu=flux_enu
    raiseVal=enu

    ! set neutrino 4-vector: incoming neutrino in z-direction
    ! plus threshold check for enu

    call eNev_init_nuStep2(eNev,enu,IP,flagOK)
    if (.not.flagOK) return


    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    !!! call getCosthetaLimits(costheta_min,costheta_max,enu)
    call getCosthetaLimits(costheta_min,costheta_max)
    call getEleptonLimits(IP,eNev,elepton_min,elepton_max)

    ! set integral precision
    intprec=integralPrecision
    intprec1=integralPrecision
    if (IP.eq.nucleon) intprec=integralPrecisionQE  ! due to small Breit-Wigner

    allocate (x(20*intprec1))
    allocate (y(20*intprec1))
    allocate (xx(20*intprec))
    allocate (yy(20*intprec))

    sigmaximum=0.

    call sg20r(costheta_min,costheta_max,intprec1,x,n1)
    call sg20r(elepton_min, elepton_max, intprec, xx, n2)

       do j=1,n1
       costhetaint=x(j)

       ! in this version the elepton_min depends on the costheta
       ! thus integration is performed only over kinematically allowed region
       !!! call minmaxE1_costheta(IP,enu,eNev%nucleon%mass,eNev%lepton_out%mass,costhetaint,elepton_min,elepton_max,success=flagOK)
       !!! if (.not.flagOK) cycle
       !!! if (debugflag) write(*,'(3(A,f12.4))') ' costheta=', costhetaint, '      E1min=', elepton_min, '    E1max=', elepton_max
       !!! call sg20r(elepton_min, elepton_max, intprec, xx, n2)

       yy = 0.0
       do l=1,n2
          eleptonint=xx(l)
          call eNev_init_nuStep3a(eNev,eleptonint,costhetaint,flagOK)
          if (.not.flagOK) cycle
          call XsecdCosthetadElepton(eNev,IP,OutPart,yy(l))
          if (yy(l).gt.sigmaximum) sigmaximum=yy(l)
       end do
       call rg20r(elepton_min,elepton_max,intprec,yy,y(j))
    end do
    call rg20r(costheta_min,costheta_max,intprec1,y,sig)

    ! set kinematics !!!! TODO: avoid infinite loop !!!!
    sigrd=0.
    if (sigmaximum.gt.0.) then
       do
          eleptonint=elepton_min+rn()*(elepton_max-elepton_min)
          costhetaint=costheta_min+rn()*(costheta_max-costheta_min)

          call eNev_init_nuStep3a(eNev,eleptonint,costhetaint,flagOK)
          if (.not.flagOK) cycle
          call XsecdCosthetadElepton(eNev,IP,OutPart, sigrd)
          if (sigmaximum*rn().le.sigrd) exit
       end do
    end if

    if (debugflag) write(*,*) 'sigmax', sigmaximum, 'sigrd', sigrd

    deallocate(x,y,xx,yy)

  end subroutine Xsec_integratedSigma


  !****************************************************************************
  !****s* neutrinoXsection/Xsec_dSigmadCosThetadElepton
  ! NAME
  ! subroutine Xsec_dSigmadCosThetadElepton(eNev,IP,raiseFlag,raiseVal,
  ! OutPart,sig,flux_enu)
  !
  ! PURPOSE
  ! This subroutine calculates dSigmadCosThetadElepton depending
  ! on the input variables:
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- the Lepton-Nucleon Event
  ! * integer             :: IP           -- ID of outgoing hadron/process
  ! * logical             :: raiseFlag    -- shall energy etc. be increased?
  ! * real, optional      :: flux_enu     -- if present, use this value as
  !   neutrino energy (used for experiments)
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * real                         :: raiseVal -- the actual value of the
  !   running variable (i.e. elepton, enu, costheta, ...)
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !****************************************************************************
  subroutine Xsec_dSigmadCosThetadElepton(eNev, IP, raiseFlag, raiseVal, &
       & OutPart,sig,flux_enu)

    use neutrino_IDTable
    use eN_event, only: eNev_init_nuStep2,eNev_init_nuStep3a



    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out
    real :: costheta_max, costheta_min
    real :: elepton_max, elepton_min
    logical :: flagOK

    if (initFlag) then
       call readInput(dSigmadCosThetadElepton)
       initFlag = .false.

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
       write(*,'(a,F12.4)') 'costheta= ', costheta
       write(*,'(a,F12.4)') 'elepton=  ', elepton
    end if


    !set default output
    call setToDefault(OutPart)
    sig=0.

    if (present(flux_enu)) enu=flux_enu

    if (raiseFlag) then
       elepton=elepton+delta_elepton
       write(*,*) 'Elepton is raised by ...', delta_elepton, ' to ', elepton
    end if
    raiseVal=Elepton

    ! set neutrino 4-vector: incoming neutrino in z-direction
    ! plus threshold check for enu

    call eNev_init_nuStep2(eNev,enu,IP,flagOK)
    if (.not.flagOK) return

    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    !check threshold for costheta and elepton
    call getCosthetaLimits(costheta_min,costheta_max)
    if (.not.checkLimits(costheta_min,costheta_max,costheta)) return

    call getEleptonLimits(IP,eNev,elepton_min,elepton_max)
    if (.not.checkLimits(elepton_min,elepton_max,elepton)) return

    call eNev_init_nuStep3a(eNev,elepton,costheta,flagOK)
    if (.not.flagOK) return

    call XsecdCosthetadElepton(eNev,IP,OutPart, sig)

  end subroutine Xsec_dSigmadCosThetadElepton





  !****************************************************************************
  !****s* neutrinoXsection/Xsec_dSigmadQsdElepton
  ! NAME
  ! subroutine Xsec_dSigmadQsdElepton(eNev,IP,raiseFlag,raiseVal,OutPart,
  ! sig,flux_enu)
  !
  ! PURPOSE
  ! This subroutine calculates dSigmadQsdElepton depending
  ! on the input variables:
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- the Lepton-Nucleon Event
  ! * integer             :: IP           -- ID of outgoing hadron/process
  ! * logical             :: raiseFlag    -- shall energy etc. be increased?
  ! * real, optional      :: flux_enu     -- if present, use this value as
  !   neutrino energy (used for experiments)
  !
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * real                         :: raiseVal -- the actual value of the
  !   running variable (i.e. elepton, enu, costheta, ...)
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !****************************************************************************
  subroutine Xsec_dSigmadQsdElepton(eNev, IP, raiseFlag, raiseVal, &
       & OutPart,sig,flux_enu)

    use neutrino_IDTable
    use eN_event, only: eNev_init_nuStep2,eNev_init_nuStep3b



    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out
    real :: elepton_max, elepton_min
    real :: Qs_min, Qs_max
    logical :: flagOK

    if (initFlag) then
       call readInput(dSigmadQsdElepton)
       initFlag = .false.

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
       write(*,'(a,F12.4)') 'Qs=      ', Qs
       write(*,'(a,F12.4)') 'elepton= ', elepton
    end if


    !set default output
    call setToDefault(OutPart)
    sig=0.

    if (present(flux_enu)) enu=flux_enu

    if (raiseFlag) then
       elepton=elepton+delta_elepton
       write(*,*) 'Elepton is raised by ...', delta_elepton, ' to ', elepton
    end if
    raiseVal=Elepton

    ! set neutrino 4-vector: incoming neutrino in z-direction
    ! plus threshold check for enu

    call eNev_init_nuStep2(eNev,enu,IP,flagOK)
    if (.not.flagOK) return

    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    !check threshold for Qs and elepton
    call getEleptonLimits(IP,eNev,elepton_min,elepton_max)
    if (.not.checkLimits(elepton_min,elepton_max,elepton)) return

    call getQsLimits(eNev,elepton,Qs_min,Qs_max)
    if (.not.checkLimits(Qs_min,Qs_max,Qs)) return

    call eNev_init_nuStep3b(eNev,elepton,Qs,flagOK)
    if (.not.flagOK) return

    call XsecdCosthetadElepton(eNev,IP,OutPart, sig)

    ! correct cross section for Jacobian:
    sig = sig/(2.*eNev%lepton_in%momentum(0)*sqrt(elepton**2-eNev%lepton_out%mass**2))

  end subroutine Xsec_dSigmadQsdElepton



  !****************************************************************************
  !****s* neutrinoXsection/Xsec_dSigmadcostheta
  ! NAME
  ! subroutine Xsec_dSigmadcostheta(eNev,IP,raiseFlag,raiseVal,OutPart,sig,
  ! flux_enu)
  !
  ! PURPOSE
  ! This subroutine calculates dSigmadcostheta depending
  ! on the input variables:
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- the Lepton-Nucleon Event
  ! * integer             :: IP           -- ID of outgoing hadron/process
  ! * logical             :: raiseFlag    -- shall energy etc. be increased?
  ! * real, optional      :: flux_enu     -- if present, use this value as
  !   neutrino energy (used for experiments)
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * real                         :: raiseVal -- the actual value of the
  !   running variable (i.e. elepton, enu, costheta, ...)
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !****************************************************************************
  subroutine Xsec_dSigmadcostheta(eNev, IP, raiseFlag, raiseVal, &
       & OutPart,sig,flux_enu)

    use neutrino_IDTable
    use idtable, only: nucleon
    use random
    use gauss_integration
    use eN_event, only: eNev_init_nuStep2,eNev_init_nuStep3a



    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out
    real :: costheta_max, costheta_min
    real :: elepton_max, elepton_min

    integer :: intprec,l,n2
    real :: eleptonint
    real :: sigmaximum,sigrd
    real, dimension(:),allocatable :: yy,xx
    logical :: flagOK

    if (initFlag) then
       call readInput(dSigmadcostheta)
       initFlag = .false.

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
       write(*,'(a,F12.4)') 'costheta= ', costheta
    end if


    ! set default output
    call setToDefault(OutPart)
    sig=0.

    if (present(flux_enu)) enu=flux_enu


    if (raiseFlag) then
       costheta=costheta+delta_costheta
       write(*,*) 'costheta is raised by ...', delta_costheta, ' to ', costheta
    end if
    raiseVal=costheta

    ! set neutrino 4-vector: incoming neutrino in z-direction
    ! plus threshold check for enu

    call eNev_init_nuStep2(eNev,enu,IP,flagOK)
    if (.not.flagOK) return

    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return


    call getCosthetaLimits(costheta_min,costheta_max)
    if (.not.checkLimits(costheta_min,costheta_max,costheta)) return

    call getEleptonLimits(IP,eNev,elepton_min,elepton_max)

    !set integral precision
    intprec=integralPrecision
    if (IP.eq.nucleon) intprec=integralPrecisionQE ! due to small Breit-Wigner

    allocate (xx(20*intprec))
    allocate (yy(20*intprec))

    call sg20r(elepton_min,elepton_max,intprec,xx,n2)
    sigmaximum=0.0
    yy = 0.0

    do l=1,n2
       eleptonint=xx(l)

       call eNev_init_nuStep3a(eNev,eleptonint,costheta,flagOK)
       if (.not.flagOK) cycle

       call XsecdCosthetadElepton(eNev,IP,OutPart, yy(l))
       if (yy(l).gt.sigmaximum) sigmaximum=yy(l)
    end do

    call rg20r(elepton_min,elepton_max,intprec,yy,sig)

    ! set kinematics
    if (sigmaximum.gt.0.) then
       do
          eleptonint=elepton_min+rn()*(elepton_max-elepton_min)

          call eNev_init_nuStep3a(eNev,eleptonint,costheta,flagOK)
          if (.not.flagOK) cycle

          call XsecdCosthetadElepton(eNev,IP,OutPart, sigrd)

          if (sigmaximum*rn().le.sigrd) exit
       end do
    end if

    deallocate(xx,yy)

  end subroutine Xsec_dSigmadcostheta


  !****************************************************************************
  !****s* neutrinoXsection/Xsec_dSigmadElepton
  ! NAME
  ! subroutine Xsec_dSigmadElepton(eNev,IP,raiseFlag,raiseVal,OutPart,sig,
  ! flux_enu)
  !
  ! PURPOSE
  ! This subroutine calculates dSigmadElepton depending
  ! on the input variables:
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- the Lepton-Nucleon Event
  ! * integer             :: IP           -- ID of outgoing hadron/process
  ! * logical             :: raiseFlag    -- shall energy etc. be increased?
  ! * real, optional      :: flux_enu     -- if present, use this value as
  !   neutrino energy (used for experiments)
  !
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * real                         :: raiseVal -- the actual value of the
  !   running variable (i.e. elepton, enu, costheta, ...)
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !****************************************************************************
  subroutine Xsec_dSigmadElepton(eNev, IP, raiseFlag, raiseVal, &
       & OutPart,sig,flux_enu)

    use neutrino_IDTable
    use idtable, only: nucleon
    use random
    use gauss_integration
    use eN_event, only: eNev_init_nuStep2,eNev_init_nuStep3a



    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out
    real :: costheta_max, costheta_min
    real :: elepton_max, elepton_min

    integer :: intprec,l,n2
    real :: costhetaint
    real :: sigmaximum,sigrd
    real,dimension(:),allocatable :: yy,xx
    logical :: flagOK


    if (initFlag) then
       call readInput(dSigmadElepton)
       initFlag = .false.

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
       write(*,'(a,F12.4)') 'elepton= ', elepton
    end if


    ! set default output
    call setToDefault(OutPart)
    sig=0.

    if (present(flux_enu)) enu=flux_enu

    if (raiseFlag) then
       elepton=elepton+delta_elepton
       write(*,*) 'Elepton is raised by ...', delta_elepton, ' to ', elepton
    end if
    raiseVal=Elepton

    ! set neutrino 4-vector: incoming neutrino in z-direction
    ! plus threshold check for enu

    call eNev_init_nuStep2(eNev,enu,IP,flagOK)
    if (.not.flagOK) return

    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    call getEleptonLimits(IP,eNev,elepton_min,elepton_max)
    if (.not.checkLimits(elepton_min,elepton_max,elepton)) return

    call getCosthetaLimits(costheta_min,costheta_max,enu)

    ! set integral precision
    intprec=integralPrecision
    if (IP.eq.nucleon) intprec=integralPrecisionQE  ! due to small Breit-Wigner

    allocate (xx(20*intprec))
    allocate (yy(20*intprec))

    call sg20r(costheta_min,costheta_max,intprec,xx,n2)
    sigmaximum=0.0
    yy = 0.0

    do l=1,n2
       costhetaint=xx(l)
       call eNev_init_nuStep3a(eNev,elepton,costhetaint,flagOK)
       if (.not.flagOK) cycle

       call XsecdCosthetadElepton(eNev,IP,OutPart, yy(l))
       if (yy(l).gt.sigmaximum) sigmaximum=yy(l)
    end do

    call rg20r(costheta_min,costheta_max,intprec,yy,sig)

    ! set kinematics
    if (sigmaximum.gt.0.) then
       do
          costhetaint=costheta_min+rn()*(costheta_max-costheta_min)

          call eNev_init_nuStep3a(eNev,elepton,costhetaint,flagOK)
          if (.not.flagOK) cycle

          call XsecdCosthetadElepton(eNev,IP,OutPart, sigrd)

          if (sigmaximum*rn().le.sigrd) exit
       end do
    end if

    deallocate(xx,yy)

  end subroutine Xsec_dSigmadElepton


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !****************************************************************************
  !****s* neutrinoXsection/XsecdCosthetadElepton
  ! NAME
  ! subroutine XsecdCosthetadElepton(eNev,IP,OutPart,sig)
  !
  ! PURPOSE
  ! This subroutine is the basic subroutine, which does all the job of
  ! calculating the cross section dSigma/dCost dElepton and generating
  ! a corresponding final state particle vector.
  !
  ! This routine is called by all routines called "Xsec_...".
  !
  ! INPUTS
  ! * type(electronNucleon_event)  :: eNev -- the Lepton-Nucleon Event
  ! * integer                      :: IP   -- ID of outgoing hadron/process
  !
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !
  ! NOTES
  ! For the outgoing particles OutPart, the following entries have to
  ! be set afterwards:
  ! * OutPart%firstEvent
  ! * OutPart%event(1:2)
  ! * OutPart%perturbative
  ! * OutPart%velocity
  ! * OutPart%offshellparameter
  ! * OutPart%perweight
  ! All the other values are set in this routine.
  !
  ! Returned cross section is in 10^-38 cm^2/GeV, except for EM, where the
  ! units are nb/GeV
  !****************************************************************************
  subroutine XsecdCosthetadElepton(eNev,IP,OutPart,sig)

    use constants, only: pi, mN, mPi, twopi
    use minkowski, only: SP
    use NeutrinoMatrixElement
    use spectralFunc, only:specfunc
    use ParticleProperties, only: hadron
    use leptonicID
    use idTable, only: nucleon,pion
    use neutrino_IDTable
    use singlePionProductionMAIDlike
    use singlePionProductionNHVlike, only: Nieves1piN_elepton_ct
    use Coll_nuN
    use callstack, only: traceback
    use lepton2p2h, only: lepton2p2h_DoQE,lepton2p2h_DoDelta
    use neutrino2piBack, only: DoNu2piBack
    use eN_event, only: nuclearFluxFactor_correction
    use eventGenerator_eN_lowEnergy, only: init_2Pi


    integer,             intent(in)  :: IP
    real,                intent(out) :: sig
    type(electronNucleon_event), intent(inout) :: eNev
    type(particle), dimension(:), intent(out) :: OutPart ! FinalState particles

    ! for internal purposes:

    integer :: process_ID
    real, dimension(0:3) :: k_in,k_out,p_in,p_out, pion_momentum_out
    real                 :: mass_out, ml_out
    integer              :: charge_out, pion_charge_out
    real, dimension(1:3) :: position
    integer              :: charge_in


    real :: kin_factor,invMassSQ,Wrec  !,Q2,Wrec2,W2,Wmed

    logical :: success
    real :: plep

    !set default output
    call setToDefault(OutPart)


    ! The charges are maybe already calculated earlier, but we do it here again
    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    ! set some abbreviations:
    process_ID = eNev%idProcess
    k_in  = eNev%lepton_in%momentum
    k_out = eNev%lepton_out%momentum
!    Q2 = SP(k_out - k_in,k_out - k_in)
    p_in  = eNev%nucleon%momentum
    p_out = p_in+k_in-k_out               ! momentum of outgoing particle

    position = eNev%nucleon%position
    charge_in = eNev%nucleon%charge
    ml_out = eNev%lepton_out%mass
    plep = sqrt(max((eNev%lepton_out%momentum(0)**2-ml_out**2),0.))
    ! ---------------

    sig=0.

 !  Now invariant mass cut on reconstructed Wrec = sqrt(mN**2 + 2*mN*nu - Q2)
 !   = invariant mass in incoming channel for free nucleon  at rest
 !  used for all production mechanisms

    Wrec = eNev%W_rec
    if (Wrec .gt. invariantMasscut) return


    select case (IP)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! QE and RES production
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (:31)

       ! cutting the resonances off at high W (default REScutW2=2.05 GeV)
       if (eNev%W_Free .gt. REScutW2) return


       invMassSQ=SP(p_out,p_out)
       if (invMassSQ.le.0.) return !reaction not possible

       if (-SP(p_out-p_in,p_out-p_in).lt.0.) then
          if (debugflag) write(*,*) 'neutrinoXsec, -SP(pf-pi,pf-pi).lt.0.', &
               & -SP(p_out-p_in,p_out-p_in), &
               & -SP(k_in-k_out,k_in-k_out), IP,charge_in,k_in,k_out
          return !reaction not possible
       end if

       kin_factor=plep/(pi*16.)/(SP(k_in,p_in))*specfunc(IP,charge_out,p_out,position,mass_out)

       if (nuclear_phasespace) kin_factor=kin_factor*nuclearFluxFactor_correction(p_in,k_in)

       if (mass_out.eq.0.0) return ! failure in specfunc

       !avoid the production of res. below threshold:
       if (mass_out.le.hadron(IP)%minMass) then
          if (debugflag) write(*,*) 'less than threshold -> sig=0'
          return
       end if

       sig=kin_factor*nuMaEl(process_ID,IP,charge_in,k_in,k_out,p_in,p_out,mass_out,position)
       !multiply cross section by factor to obtain results in 10^(-38) cm^2
       sig=corr_factor*sig

       ! smooth switching off the resonances
       if (eNev%W_Free .gt. REScutW1) &
            & sig = sig*(REScutW2-eNev%W_Free)/(REScutW2-REScutW1)


       ! Setting the outgoing particles:

       OutPart%position(1)=position(1)
       OutPart%position(2)=position(2)
       OutPart%position(3)=position(3)

       OutPart%formationTime = -999

       OutPart(1)%mass = mass_out
       OutPart(1)%ID   = IP
       OutPart(1)%Charge = Charge_Out
       OutPart(1)%antiparticle=.false.
       OutPart(1)%Momentum = p_out



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! single pi N
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (onePionCH_n,onePionCH_p)

       ! cutting the background off at high W (default REScutW2=2.05 GeV)
       if (eNev%W_Free .gt. REScutW2) return

       ! invariant mass cut for pion BG contribution
       if (Wrec .gt. invariantMasscut_BG) return

       ! remarks: the routines to be called here should return:
       ! - dsigma/(dElepton dcostheta_lepton) in units of GeV**4 for CC,EM,NC
       ! - p_out, mass_out
       ! - pion_momentum_out

       select case (singlePiModel)

       case (0) !HNV model

          !if (debugflag) write(*,*) 'singlePiModel nieves started'
          call Nieves1piN_elepton_ct( process_ID,k_in,k_out,.true., &
               & p_in,position,charge_in, &
               & p_out,charge_out,pion_momentum_out,pion_charge_out,sig)

          mass_out=mN

       case (1) !MAID like

          !check on outgoing particles
          invMassSQ=SP(p_out,p_out)
          if (invMassSQ.le.0.) return !reaction not possible

          sig=MAIDlike_singlePi(eNev,charge_out,pion_charge_out, &
               & p_out, pion_momentum_out,nuclear_phasespace)

          mass_out=mN

       case default

          write(*,*) 'wrong choice for singlePiModel -> stop',singlePiModel
          call traceback()

       end select

       !multiply cross section by factor to obtain results in 10^(-38) cm^2
       sig=corr_factor*sig

       ! smooth switching off the  background
       if (eNev%W_Free .gt. REScutW1) &
            & sig = sig*(REScutW2-eNev%W_Free)/(REScutW2-REScutW1)

       ! Setting the outgoing particles:

       OutPart%position(1)=position(1)
       OutPart%position(2)=position(2)
       OutPart%position(3)=position(3)

       OutPart%formationTime = -999

       OutPart(1:2)%ID = (/nucleon,pion/)
       OutPart(1:2)%Mass=(/mass_out,mPi/)
       OutPart(1:2)%Charge=(/charge_out,pion_charge_out/)
       OutPart(1:2)%antiparticle=.false.

       OutPart(1)%Momentum = p_out
       OutPart(2)%Momentum = pion_momentum_out

    case default
       write(*,*) 'strange IP, STOP:',IP
       call traceback()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (DIS_CH)

       if (eNev%W_Free .lt. DIScutW1) return

       !check on outgoing particles
       invMassSQ=SP(p_out,p_out)
       if (invMassSQ.le.0.) return !reaction not possible

       OutPart%position(1)=position(1)
       OutPart%position(2)=position(2)
       OutPart%position(3)=position(3)

       call DoColl_nuN_Py(eNev,OutPart,success, DISmassless, sig)

       ! correction for the nuclear phase space
       if (nuclear_phasespace) sig=sig*nuclearFluxFactor_correction(p_in,k_in)

       ! aditional form factor:
       select case (eNev%idProcess)
       case (-1,1) !=== EM
          select case (DISformfakEM)
          case (1)
             sig = sig * eNev%Qsquared/(eNev%Qsquared+mcutDIS**2)
          case (2)
             sig = sig * (eNev%Qsquared/(eNev%Qsquared+mcutDIS**2))**2
          end select
       case default !=== CC, NC
          select case (DISformfakNCCC)
          case (1)
             sig = sig * eNev%Qsquared/(eNev%Qsquared+mcutDIS**2)
          case (2)
             sig = sig * (eNev%Qsquared/(eNev%Qsquared+mcutDIS**2))**2
          end select
       end select

       ! returned XS is dsigma/dE'dcost in mb/GeV
       if (.not.success) sig=0.0
       sig = sig * 1e11 ! cross section in 10^-38 cm^2 (from mbar=10^{-27} to 10^{-38})

       if (eNev%W_Free .lt. DIScutW2) &
            & sig = sig*(eNev%W_Free-DIScutW1)/(DIScutW2-DIScutW1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 2p2h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (QE2p2h)


       call lepton2p2h_DoQE(eNev,outPart,sig)
       ! returned XS is dsigma/dE'dOmega in mb/GeV/A
       sig = sig * twopi * 1.e11 ! cross section dsigma/dE'dcost in 10^-38 cm^2

    case (Delta2p2h)
       call lepton2p2h_DoDelta(eNev,outPart,sig)
       ! returned XS is dsigma/dE'dOmega in mb/GeV/A
       sig = sig * twopi * 1.e11 ! cross section dsigma/dE'dcost in 10^-38 cm^2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! two pion background
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (twoPion)

       ! mass cut:
       if (eNev%W_Free .lt. DIScutW1) return

       invMassSQ=SP(p_out,p_out)
       if (invMassSQ.le.0.) return ! reaction not possible

       select case (abs(process_ID))
       case (1)
          call init_2Pi(eNev,OutPart,sig, .true.)
       case (2)
          call DoNu2piBack(eNev,outPart,sig)
       end select

       ! returned XS is dsigma/dE'dOmega in mb/GeV
       sig = sig * twopi * 1e11 ! cross section dsigma/dE'dcost in 10^-38 cm^2

       if (process_ID.eq.antiCC) sig=sig/2.    !decrease xsec for antineutrinos

    end select

    if (process_ID.eq.EM) sig=1e-5*sig    !cross section in nanobarn for el-m reactions

!    write(*,*) '>> XsecdCosthetadElepton',IP,sig

  end subroutine XsecdCosthetadElepton



  !****************************************************************************
  !****f* neutrinoXsection/SetHadronCharge
  ! NAME
  ! logical function SetHadronCharge(eNev,IP,Q_R,Q_pi)
  !
  ! PURPOSE
  ! This function sets the charge of the outgoing hadrons depending on the
  ! reaction process. If successful, SetHadronCharge=.true., if not (when
  ! the reaction is not possible) SetHadronCharge=.false.
  !
  ! INPUTS
  ! * type(electronNucleon_event)  :: eNev -- the Lepton-Nucleon Event
  ! * integer                      :: IP   -- ID of reaction
  !
  ! OUTPUT
  ! * integer                      :: Q_R  -- charge of outgoing baryon
  ! * integer                      :: Q_pi -- charge of outgoing pion
  ! * function value -- indicates possible failures
  !****************************************************************************
  logical function SetHadronCharge(eNev,IP,Q_R,Q_pi)
    use neutrino_IDTable
    use leptonicID
    use ParticleProperties, only: hadron

    type(electronNucleon_event), intent(in) :: eNev
    integer, intent(in) :: IP
    integer, intent(out) :: Q_R, Q_pi

    !set default output
    Q_R=0
    Q_pi=0

    SetHadronCharge=.true.

    select case (IP)

    case (1:31)
       !=======================================================================
       !===== QE and RES production
       !=======================================================================
       Q_pi = 0

       select case (eNev%idProcess)
       case default
          ! == EM, NC, antiEM, antiNC
          Q_R = eNev%nucleon%charge
       case (2)
          ! == CC
          Q_R = eNev%nucleon%charge + 1
          if ((Q_R.eq.2).and.(hadron(IP)%isoSpinTimes2.ne.3)) &
               & SetHadronCharge=.false. ! reaction not possible
       case (-2)
          ! == antiCC
          Q_R = eNev%nucleon%charge - 1
          if ((Q_R.eq.-1).and.(hadron(IP)%isoSpinTimes2.ne.3)) &
               & SetHadronCharge=.false. ! reaction not possible

       end select

    case (onePionCH_n)
       !=======================================================================
       !===== pi + n production
       !=======================================================================

       !CC:    nu + n -> l- + pi+ + n  ;  nu~ + n -> l+  + pi- + n
       !CC:    nu + p ->  ---          ;  nu~ + p -> l+  + pi0 + n
       !NC/EM: nu + n -> nu + pi0 + n  ;  nu~ + n -> nu~ + pi0 + n
       !NC/EM: nu + p -> nu + pi+ + n  ;  nu~ + p -> nu~ + pi+ + n

       Q_R = 0
       select case (eNev%idProcess)
       case default ! == EM, NC, antiEM, antiNC
          Q_pi = eNev%nucleon%charge
       case (2)    ! == CC
          Q_pi = eNev%nucleon%charge+1
          if (Q_pi.eq.2) SetHadronCharge=.false. ! reaction not possible
       case (-2)   ! == antiCC
          Q_pi = eNev%nucleon%charge-1
       end select

    case (onePionCH_p)
       !=======================================================================
       !===== pi + p production
       !=======================================================================

       !CC:    nu + n -> l- + pi0 + p  ;  nu~ + n ->  ---
       !CC:    nu + p -> l- + pi+ + p  ;  nu~ + p -> l+  + pi- + p
       !NC/EM: nu + n -> nu + pi- + p  ;  nu~ + n -> nu~ + pi- + p
       !NC/EM: nu + p -> nu + pi0 + p  ;  nu~ + p -> nu~ + pi0 + p

       Q_R = 1
       select case (eNev%idProcess)
       case default ! == EM, NC, antiEM, antiNC
          Q_pi = eNev%nucleon%charge-1
       case (2)    ! == CC
          Q_pi = eNev%nucleon%charge
       case (-2)   ! == antiCC
          Q_pi = eNev%nucleon%charge-2
          if (Q_pi.eq.-2) SetHadronCharge=.false. ! reaction not possible
       end select

    case (DIS_CH,QE2p2h,Delta2p2h)
       !=======================================================================
       !===== DIS, 2p2h
       !=======================================================================
       !... everything as the defaults

    case (twoPion)
       !=======================================================================
       !===== 2 pion backround
       !=======================================================================
       !... everything as the defaults

    case default
       !=======================================================================
       !===== DEFAULT
       !=======================================================================
       write(*,*) 'wrong IP in SetHadronCharge:',IP,' STOP!'
       stop

    end select

  end function SetHadronCharge


  !****************************************************************************
  !****s* neutrinoXsection/getCosthetaLimits
  ! NAME
  ! subroutine getCosthetaLimits(costheta_min,costheta_max)
  !
  ! PURPOSE
  ! This subroutine returns the limits in costheta.
  !
  ! OUTPUT
  ! * real             :: costheta_min
  ! * real             :: costheta_max
  !****************************************************************************
  subroutine getCosthetaLimits(costheta_min,costheta_max,Enu)

    real, intent(out) :: costheta_min,costheta_max
    real, intent(in), optional :: Enu
    costheta_min=-1.
    costheta_max=1.
    if (present(Enu)) then
        ! practical cut for high neutrino energies - the formula is pure "educated guess"
        ! the reason for cut is physical --- above Q^2>4 everything dies,
        ! so for high neutrino energies only forward scattering is possible
        if (Enu.ge.2) costheta_min=1.-4./Enu/Enu
    end if

  end subroutine getCosthetaLimits


  !****************************************************************************
  !****s* neutrinoXsection/getEleptonLimits
  ! NAME
  ! subroutine getEleptonLimits(IP,eNev,elepton_min,elepton_max)
  !
  ! PURPOSE
  ! This subroutine returns the limits in elepton depending on the neutrino
  ! energy.
  !
  ! INPUTS
  ! * type(electronNucleon_event)  :: eNev -- the Lepton-Nucleon Event
  ! * integer                      :: IP   -- ID of outgoing hadron/process
  !
  ! OUTPUT
  ! * real             :: elepton_min
  ! * real             :: elepton_max
  !****************************************************************************
  subroutine getEleptonLimits(IP,eNev,elepton_min,elepton_max)
    use idtable, only: nucleon,delta
    use constants, only: mPi

    integer, intent(in) :: IP
    real, intent(out) :: elepton_min,elepton_max
    type(electronNucleon_event), intent(in) :: eNev

    real :: enu, ml_out

    enu = eNev%lepton_in%momentum(0)
    ml_out = eNev%lepton_out%mass

    elepton_min=ml_out
    if (IP.eq.nucleon) elepton_max=enu
    if (IP.ge.delta) elepton_max=enu-mPi
! IP=35 = 2p2h process
    if (IP.eq.35) elepton_max=enu
  end subroutine getEleptonLimits


  !****************************************************************************
  !****s* neutrinoXsection/getQsLimits
  ! NAME
  ! subroutine getQsLimits(eNev,elepton,Qs_min,Qs_max)
  !
  ! PURPOSE
  ! This subroutine returns the limits in Qs depending on the neutrino energy
  ! and the energy of the outgoing lepton.
  !
  ! INPUTS
  ! * type(electronNucleon_event)  :: eNev -- the Lepton-Nucleon Event
  ! * real                         :: elepton -- energy of outgoing lepton
  !
  ! OUTPUT
  ! * real             :: Qs_min
  ! * real             :: Qs_max
  !****************************************************************************
  subroutine getQsLimits(eNev,elepton,Qs_min,Qs_max)

    type(electronNucleon_event), intent(in) :: eNev
    real, intent(in) :: elepton
    real, intent(out) :: Qs_min,Qs_max

    real :: enu,ml_out

    enu = eNev%lepton_in%momentum(0)
    ml_out = eNev%lepton_out%mass

    !! Attention! you should have checked before, that elepton>ml_out
    !! is guaranteed. Therefore we could skip the "sqrt(max(" stuff !!

    Qs_max=-ml_out**2+2.*enu*(elepton+sqrt(max((elepton**2-ml_out**2),0.)))
    Qs_min=-ml_out**2+2.*enu*(elepton-sqrt(max((elepton**2-ml_out**2),0.)))
  end subroutine getQsLimits

  !****************************************************************************
  !****f* neutrinoXsection/checkLimits
  ! NAME
  ! logical function checkLimits(X_min,X_max,X)
  !
  ! PURPOSE
  ! This function checks whether X is out of bounds.
  ! If not, checkLimits=.true., if yes checkQsLimits=.false.
  !
  ! INPUTS
  ! * real  :: X_min,X_max,X
  !****************************************************************************
  logical function checkLimits(X_min,X_max,X)

    real, intent(in) :: X_min,X_max,X

    checkLimits = .false.
    if (X.lt.X_min) return
    if (X.gt.X_max) return
    checkLimits = .true.

  end function checkLimits



  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************




  !****************************************************************************
  !****s* neutrinoXsection/get_xsection_namelist
  ! NAME
  ! subroutine get_xsection_namelist
  !
  ! PURPOSE
  ! This subroutine returns variables that are set in the according xsection namlist.
  !
  ! OUTPUT
  ! * logical,optional :: Gdebugflag,Gnuclear_phasespace
  ! * integer,optional :: GsinglePiModel,GintegralPrecision,GintegralPrecisionQE
  ! * real,optional :: Genu, Gdelta_enu, GQs,Gdelta_Qs, Gcostheta, Gdelta_costheta, Gdelta_elepton, GinvariantMasscut
  !
  !****************************************************************************

  subroutine get_xsection_namelist(XsectionMode,Gdebugflag,Gnuclear_phasespace,GsinglePiModel,GintegralPrecision, &
                            & GintegralPrecisionQE,Genu,Gdelta_enu,GQs,Gdelta_Qs,Gcostheta,Gdelta_costheta, &
                            & Gelepton,Gdelta_elepton,GinvariantMasscut)

    integer, intent(in), optional ::  XsectionMode
    logical,optional,intent(out) :: Gdebugflag,Gnuclear_phasespace
    integer,optional,intent(out) :: GsinglePiModel,GintegralPrecision,GintegralPrecisionQE
    real,optional,intent(out) :: Genu,Gdelta_enu,GQs,Gdelta_Qs,Gcostheta,Gdelta_costheta,Gelepton,Gdelta_elepton,GinvariantMasscut

    if (present (XsectionMode)) then
        if (initflag) then
           call readInput(MOD(XsectionMode,10),(XsectionMode>10))
           initflag = .false.
        end if
    end if

    if (present(Gdebugflag)) Gdebugflag=debugflag
    if (present(Gnuclear_phasespace)) Gnuclear_phasespace=nuclear_phasespace
    if (present(GsinglePiModel)) GsinglePiModel=singlePiModel
    if (present(GintegralPrecision))GintegralPrecision =integralPrecision
    if (present(GintegralPrecisionQE))GintegralPrecisionQE =integralPrecisionQE
    if (present(Genu)) Genu=enu
    if (present(Gdelta_enu))Gdelta_enu= delta_enu
    if (present(GQs)) GQs =Qs
    if (present(Gdelta_Qs))Gdelta_Qs =delta_Qs
    if (present(Gcostheta))Gcostheta =costheta
    if (present(Gdelta_costheta))Gdelta_costheta =delta_costheta
    if (present(Gelepton))Gelepton =elepton
    if (present(Gdelta_elepton))Gdelta_elepton =delta_elepton
    if (present(GinvariantMasscut))GinvariantMasscut =invariantMasscut
  end subroutine get_xsection_namelist

end module neutrinoXsection
