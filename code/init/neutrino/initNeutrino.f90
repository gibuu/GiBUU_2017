!******************************************************************************
!****m* /initNeutrino
! NAME
! module initNeutrino
!
! PURPOSE
! This module is the main module for neutrino-induced reactions. It handles
! the calculation of various neutrino cross sections, depending on the
! choice of the user for each entry in the particle vector and finally, it
! sets the produced output in the pertubativeVector.
!******************************************************************************
module initNeutrino

  use neutrino_IDTable
  use Electron_origin, only: origin_singlePi, origin_doublePi, &
       origin_DIS, origin_2p2hQE, origin_2p2hDelta

  implicit none

  private

  public :: init_neutrino
  public :: getFirstEventRange
  public :: getNeutrinoInfo
  public :: nuXsectionMode, process_ID, flavor_ID
  public :: nuExp
  public :: includeQE, includeDELTA, includeRES, include1pi, includeDIS, &
       & include2p2hQE, include2p2hDelta, include2pi
  public :: neutrinoInit_getRealRun
  public :: cleanup

  public :: get_init_namelist
  public :: get_runtime_vars
  public :: max_finalstate_ID, max_Hist,includeHist,K2Hist, &
       numberOfExperiments,&
       OscLength,Osc



  !****************************************************************************
  !****g* initNeutrino/process_ID
  ! SOURCE
  integer, save :: process_ID = 2
  ! PURPOSE
  ! Determine the process (cf. module leptonicID):
  ! * 1 = EM
  ! * 2 = CC
  ! * 3 = NC
  ! * -1 = antiEM
  ! * -2 = antiCC
  ! * -3 = antiNC
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/flavor_ID
  ! SOURCE
  integer, save :: flavor_ID = 2
  ! PURPOSE
  ! Determine the lepton flavor:
  ! * 1 = electron
  ! * 2 = muon
  ! * 3 = tau
  !****************************************************************************

  character*(*), dimension(3), parameter :: sProcess = (/"EM","CC","NC"/)
  character*(*), dimension(3), parameter :: sFamily  = (/"e  ","mu ","tau"/)

  !****************************************************************************
  !****g* initNeutrino/nuXsectionMode
  ! SOURCE
  integer, save :: nuXsectionMode = 0
  ! PURPOSE
  ! To choose which kind of Xsection is calculated. All values set in
  ! module neutrino_IDTable.f90
  !
  ! possible values:
  ! * 0 = integratedSigma: required input: enu
  ! * 1 = dSigmadCosThetadElepton: required input: enu, costheta, elepton
  ! * 2 = dSigmadQsdElepton: required input: enu, Qs, elepton
  ! * 3 = dSigmadQs: required input: enu, Qs
  ! * 4 = dSigmadCosTheta: required input: enu, costheta
  ! * 5 = dSigmadElepton: required input: enu, elepton
  ! * 6 = dSigmaMC: required input: enu
  ! * 7 = dSigmadW: required input: enu, W
  !
  ! calculation for specific experiments taking into account the flux
  ! (choose your favorite experiment with flag nuExp):
  ! * 10 = EXP_dSigmadEnu
  ! * 11 = EXP_dSigmadCosThetadElepton
  ! * 12 = EXP_dSigmadQsdElepton
  ! * 13 = EXP_dSigmadQs
  ! * 14 = EXP_dSigmadCosTheta
  ! * 15 = EXP_dSigmadElepton
  ! * 16 = EXP_dSigmaMC
  ! * 17 = EXP_dSigmadW
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/nuExp
  ! SOURCE
  integer, save :: nuExp = 0
  ! PURPOSE
  ! * 0 = no specific experiment
  ! * 1 = MiniBooNE neutrino flux (in neutrino mode = positive polarity)
  ! * 2 = ANL
  ! * 3 = K2K
  ! * 4 = BNL
  ! * 5 = MiniBooNE antienutrino flux (in antineutrino mode = negative polarity)
  ! * 6 = MINOS muon-neutrino  in neutrino mode
  ! * 7 = MINOS muon-antineutrino  in neutrino mode
  ! * 8 = NOVA neutrino (medium energy NuMI, 14 mrad off-axis), FD
  ! * 9 = T2K neutrino off-axix 2.5 degrees ( at ND280 detector )
  ! * 10 = uniform distribution from Eflux_min to Eflux_max
  !       (see namelist nl_neutrino_energyFlux in the module expNeutrinoFluxes)
  ! * 11 = MINOS muon-neutrino  in antineutrino mode
  ! * 12 = MINOS muon-antineutrino  in antineutrino mode
  ! * 13 = MINERvA muon neutrino, old flux
  ! * 14 = MINERvA muon antineutrino, old flux
  ! * 15 = LBNF/DUNE neutrino in neutrino mode
  ! * 16 = LBNF/DUNE antineutrino in antineutrino mode
  ! * 17 = LBNO neutrino in neutrino mode
  ! * 18 = NOMAD
  ! * 19 = BNB nue          BNB= Booster Neutrino Beam
  ! * 20 = BNB nuebar
  ! * 21 = BNB numu
  ! * 22 = BNB numubar
  ! * 23 = NOvA ND
  ! * 24 = T2K on axis
  ! * 25 = MINERvA, 2016 flux
  integer, parameter :: numberOfExperiments=25
  !****************************************************************************


  character*(*), dimension(0:numberOfExperiments), parameter ::  sExp  = (/ &
        "no specific experiment ", &
        "MiniBooNE nu           ", "ANL                    ", &
        "K2K                    ", &
        "BNL                    ", "MiniBooNE barnu        ", &
        "MINOS nu numode        ", "MINOS barnu numode     ", &
        "NOvA FD                ", &
        "T2K OffAxis 2.5deg     ", "uniform distribution   ", &
        "MINOS nu barnumode     ", "MINOS barnu barnumode  ", &
        "MINERvA nu numode      ", "MINERvA barnu barnumode", &
        "LBNF-DUNE nu           ", "LBNF-DUNE barnu        ", &
        "LBNO nu numode         ", "NOMAD                  ", &
        "BNB nue                ", "BNB nuebar             ", &
        "BNB numu               ", "BNB numubar            ", &
        "NOvA ND                ", "T2K on axis            ", &
        "MINERvA, 2016 flux     " /)


  real, dimension(0:numberOfExperiments), parameter :: OscLength = &
  (/ 0., 0.541, 0., 250., 0., 0.541, 735., 735., 810., 295., 0., 735., 735., &
     0.5, 0.5, 1300., 1300., 2300., 0.6262,0.,0.,0.,0.,0.,0.,0. /)
  ! oscillation length for various experiments in kilometers

  logical, dimension(0:numberOfExperiments), parameter:: Osc = &
  (/ .FALSE.,.FALSE.,.FALSE.,.TRUE.,.FALSE.,.FALSE.,.TRUE.,.TRUE.,.TRUE.,&
     .TRUE.,.FALSE.,.TRUE.,.TRUE.,.FALSE.,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,&
     .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE. /)
  ! OSC is true for oscillation experiments, false otherwise
  !


  !****************************************************************************
  !****g* initNeutrino/debugFlag
  ! SOURCE
  logical, parameter :: debugFlag = .false.
  ! PURPOSE
  ! To switch on debugging information
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/includeQE
  ! SOURCE
  logical, save :: includeQE = .true.
  ! PURPOSE
  ! include QE scattering
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/includeDELTA
  ! SOURCE
  logical, save :: includeDELTA = .true.
  ! PURPOSE
  ! include Delta excitation
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/includeRES
  ! SOURCE
  logical, save :: includeRES = .true.
  ! PURPOSE
  ! include excitation of higher resonances
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/include1pi
  ! SOURCE
  logical, save :: include1pi = .false.
  ! PURPOSE
  ! include one-pion cross section
  ! see neutrinoXsection.f90 for details: there
  ! one might choose between different models and
  ! also whether it is taken as background or as
  ! total cross section
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/include2pi
  ! SOURCE
  logical, save :: include2pi = .false.
  ! PURPOSE
  ! include 2 pion background channel
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/includeDIS
  ! SOURCE
  logical, save :: includeDIS = .false.
  ! PURPOSE
  ! include DIS contribution
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/include2p2hQE
  ! SOURCE
  logical, save :: include2p2hQE = .false.
  ! PURPOSE
  ! include 2p2h QE contribution
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/include2p2hDelta
  ! SOURCE
  logical, save :: include2p2hDelta = .false.
  ! PURPOSE
  ! include 2p2h Delta contribution
  !****************************************************************************


  !****************************************************************************
  !****g* initNeutrino/realRun
  ! SOURCE
  logical, save :: realRun = .false.
  ! PURPOSE
  ! Do not initialize the final state particles as perturbative particles but
  ! as real ones.
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/printAbsorptionXS
  ! SOURCE
  logical, save :: printAbsorptionXS = .false.
  ! PURPOSE
  ! flag to produce output about inclusive (absorption) cross sections
  !****************************************************************************

  real, save :: raiseVal


  !****************************************************************************
  !****g* initNeutrino/max_finalstate_ID
  ! SOURCE
  integer,parameter :: max_finalstate_ID=37
  ! * Parameter determines the reaction mechanism and kind of final states
  ! * Final states are numbered (often IP) by
  ! * 1: nucleon (QE)
  ! * 2-31: non-strange baryon resonance (as in IdTable)
  ! * 32: pi neutron-background  (e.g. nu + n -> mu + pi+ + n)
  ! * 33: pi proton-background   (e.g. nu + n -> mu + pi0 + p)
  ! * 34: DIS
  ! * 35: 2p2h QE
  ! * 36: 2p2h Delta
  ! * 37: two pion background
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/sigmacut
  ! SOURCE
  real, save :: sigmacut=10e-4
  ! PURPOSE
  ! events with a cross section smaller than this value are skipped.
  !****************************************************************************


  logical, save :: initFlag = .true.
  logical, save :: readinputflag = .true.
  integer, save :: first=0

  !resulting cross sections that can be accessed by analysis routines
  real,dimension(0:max_finalstate_ID)      :: sigabsArr
  real,dimension(0:max_finalstate_ID),save :: sigabsArrFinal

  logical, dimension(1:max_finalstate_ID) :: includeK
  integer, dimension(1:max_finalstate_ID),parameter :: K2Hist = (/&
  ! 1=QE, 2=Delta, 3=highRES, 4=1piBG, 5=DIS, 6=2p2hQE, 7=2p2hDelta, 8=2pi
       & 1,2,3,3,3,3,3,3,3,3,&
       & 3,3,3,3,3,3,3,3,3,3,&
       & 3,3,3,3,3,3,3,3,3,3,&
       & 3,4,4,5,6,7,8/)
  ! K2Hist controls output of various reaction components.
  ! The first 31 elements stand for the number of nucleon resonances taken
  ! into account. The many '3' have the effect that all higher-lying
  ! resonances beyond the Delta are for the analysis being summed into
  ! one effective higher-lying resonance.

  integer, parameter :: max_Hist = 8 ! number of different histograms

  logical, dimension(0:max_Hist), save :: includeHist

  ! needed for checkEvent:
  integer,dimension(1:max_finalstate_ID),parameter :: eOrigin = (/&
       1,2,3,4,5,6,7,8,9,10, &
       11,12,13,14,15,16,17,18,19,20, &
       21,22,23,24,25,26,27,28,29,30, &
       31,origin_singlePi,origin_singlePi,origin_DIS, &
       origin_2p2hQE, origin_2p2hDelta, origin_doublePi /)

contains

  !****************************************************************************
  !****s* initNeutrino/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! This subroutine reads input out of jobcard
  ! from namelist 'neutrino_induced'.
  !****************************************************************************
  subroutine readInput
    use output
    use esample

    integer :: ios

    !**************************************************************************
    !****n* initNeutrino/neutrino_induced
    ! NAME
    ! NAMELIST /neutrino_induced/
    ! PURPOSE
    ! This Namelist includes:
    ! * process_ID
    ! * flavor_ID
    ! * nuXsectionMode
    ! * nuExp
    ! * includeQE
    ! * includeDELTA
    ! * includeRES
    ! * include1pi
    ! * include2pi
    ! * includeDIS
    ! * include2p2hQE
    ! * include2p2hDelta
    ! * sigmacut
    ! * realRun
    ! * printAbsorptionXS
    !**************************************************************************
    NAMELIST /neutrino_induced/  process_ID,flavor_ID,nuXsectionMode,nuExp, &
         includeQE,includeDELTA,includeRES,include1pi,includeDIS,&
         include2p2hQE, include2p2hDelta, include2pi, &
         sigmacut, realRun, printAbsorptionXS

    if (.not.readinputflag) return

    call Write_ReadingInput('neutrino_induced',0)
    rewind(5)
    read(5,nml=neutrino_induced,IOSTAT=ios)
    call Write_ReadingInput("neutrino_induced",0,ios)

    ! block not-implemented values:
    if (include2pi) then
       write(*,*) "2pi bg for neutrinos not implemented, turned off"
       include2pi = .FALSE.
    end if
    write(*,*) 'include2pi =', include2pi

    select case (flavor_ID)
    case (electron,muon,taulepton)
       write(*,'(A,A,i5)') ' ...flavor = ',sFamily(flavor_ID)
    case default
       write(*,'(A,A,i5)') ' ...flavor = ','***unknown*** ', flavor_ID
       stop
    end select

    select case (process_ID)
    case (1:3)
       write(*,'(A,A,i5)') ' ...Process: ',sProcess(process_ID)
    case (-3:-1)
       write(*,'(A,A,i5)') ' ...Process: anti-',sProcess(-process_ID)
    case default
       write(*,'(A,A,i5)') ' ...Process: ','***unknown*** ',process_ID
       stop
    end select

    write(*,'(a,2I2)') ' nuXsectionMode,nuExp:',nuXsectionMode,nuExp



    write(*,*)
    write(*,'(a,L2)')' ...include QE        : ',includeQE
    write(*,'(a,L2)')' ...include Delta     : ',includeDELTA
    write(*,'(a,L2)')' ...include higher RES: ',includeRES
    write(*,'(a,L2)')' ...include 1pi       : ',include1pi
    write(*,'(a,L2)')' ...include 2pi       : ',include2pi
    write(*,'(a,L2)')' ...include DIS       : ',includeDIS
    write(*,'(a,L2)')' ...include 2p2h QE   : ',include2p2hQE
    write(*,'(a,L2)')' ...include 2p2h Delta: ',include2p2hDelta


    if (nuXsectionMode.ge.10) then
       select case (nuExp)
       case (1:numberOfExperiments)
          write(*,*) '##### calculation is done for the ',trim(sExp(nuExp)),&
               & ' experiment #####'
       case default
          write(*,*) 'combination nuXsectionMode.ge.10.and.nuExp makes no sense -> STOP', &
             & nuexp,nuXsectionmode
          stop
       end select
    end if

    if (nuExp.gt.0.and.nuXsectionMode.lt.10) then
       write(*,*) 'combination nuExp.gt.0.and.nuXsectionMode.lt.10 makes no sense -> STOP', &
          & nuexp,nuXsectionmode
       stop
    end if


    if (nuExp.eq.1   .and. process_ID.le.0) then
       write(*,*) 'combination of MiniBooNE neutrino flux and antineutrino run makes no sense &
          &  -> STOP', nuexp, process_ID
       stop
    end if

    if (nuExp.eq.5   .and. process_ID.ge.0) then
       write(*,*) 'combination of MiniBooNE antineutrino flux and neutrino run makes no sense &
          & -> STOP', nuexp, process_ID
       stop
    end if


    if ((nuExp.eq.6 .or. nuExp.eq.11)  .and. process_ID.le.0) then
       write(*,*) 'combination of MINOS neutrino flux and antineutrino run makes no sense &
          & -> STOP', nuexp, process_ID
       stop
    end if

    if ((nuExp.eq.7  .or. nuExp.eq.12)  .and. process_ID.ge.0) then
       write(*,*) 'combination of MINOS antineutrino flux and neutrino run makes no sense &
          & -> STOP', nuexp, process_ID
       stop
    end if

    if (realRun) write(*,*) '#### REAL RUN ####'

    call Write_ReadingInput('neutrino_induced',1)

    readinputflag = .false.


  end subroutine readInput

  !****************************************************************************
  !****************************************************************************
  subroutine cleanUp
    use formfactors_A_main, only: cleanupMAID => cleanup

    call cleanupMAID
  end subroutine cleanUp


  !****************************************************************************
  !****s* initNeutrino/init_neutrino
  ! NAME
  ! subroutine init_neutrino(realParticles,pertParticles,raiseFlagIn,
  ! num_runs_sameEnergy,targetNuc)
  !
  ! PURPOSE
  ! This subroutine initializes a neutrino event on each nucleon in the
  ! realparticles vector (given to the routine). The resulting particles
  ! are set into the pertParticles vector. The reaction process is
  ! determined by the values read in by 'readInput':
  ! * process_ID
  ! * flavor_ID
  ! * nuXsectionMode
  !
  ! The user might also choose which contributions should be included
  ! (QE, DELTA and/or higher resonances) and whether the calculation should
  ! be done for a specific experiment.
  !
  ! INPUTS
  ! * type(particle), dimension(:,:) :: realParticles
  ! * integer                        :: num_runs_sameEnergy
  ! * logical                        :: raiseFlagIn -- if .true. then the
  !   energy etc is raised by raiseValue
  ! * type(tnucleus), pointer        :: targetNuc
  !
  ! OUTPUT
  ! * type(particle), dimension(:,:) :: pertParticles
  !
  !****************************************************************************
  subroutine init_neutrino(realParticles,pertParticles,raiseFlagIn,&
       & num_runs_sameEnergy,targetNuc)
    use particleDefinition
    use random, only: rn
    use idtable
    use pauliBlockingModule, only: checkPauli
    use neutrinoXsection
    use propagation, only: gradients,updateVelocity
    use collisionNumbering, only: pert_numbering,real_numbering
    use insertion, only: setIntoVector
    use inputGeneral, only: fullEnsemble
    use output, only: WriteParticleVector,WriteParticle_debug,IntToChar
    use histf90
    use hist2Df90
    use neutrinoInfoStorage
    use neutrinoProdInfo, only: neutrinoProdInfo_Init,neutrinoProdInfo_Store
    use offShellPotential
    use ExpNeutrinofluxes
    use esample
    use eN_eventDefinition
    use eN_event
    use eventGenerator_eN_lowEnergy, only: checkEvent
    use CallStack
    use Coll_nuN, only: CalcXY
    use nucleusDefinition
    use NuclearPDF, only: SetNuclearPDFA
    use monteCarlo, only: MonteCarloChoose
    use constants, only: mN


    type(particle), dimension(:,:),intent(inOut) :: realParticles
    type(particle), dimension(:,:),intent(inOut) :: pertParticles
    integer, intent(in) :: num_runs_sameEnergy
    logical, intent(in) :: raiseFlagIn
    type(tnucleus), pointer :: targetNuc

    logical :: raiseFlag

    integer :: numtry=0
    integer,dimension(1:max_finalstate_ID)       :: numberofsuccess=0
    integer,dimension(1:max_finalstate_ID), save :: numberofsuccessfinal=0
    real :: sigtot!,sig
    real, dimension(1:max_finalstate_ID) ::  sigma=0.

    type(particle),dimension(1:max_finalstate_ID,1:2),target :: OutPart
    type(particle),dimension(1:20),target :: OutPartDIS
    type(particle),dimension(1:3), target :: OutPart2pi
    type(particle),dimension(:), pointer :: pOutPart
    type(particle),dimension(:), allocatable :: finalstate
    integer :: i, j, k

    real :: totalWeight
    logical :: setflag
    integer, save :: numberofcalls=0
    integer :: countfluxcutoff
    integer :: countsigmacutoff
    type(histogram), save :: energyInit

    real :: flux_enu
    integer :: firstEvent
    integer :: whichReal_nucleon,numNucleons
    integer :: failuresV=0,failuresO=0
    logical :: success
    real,dimension(1:3)  :: grad_P
    real :: fak1,fak2


    type(electronNucleon_event) :: eNeV0,eNev1
    type(electronNucleon_event),dimension(1:max_finalstate_ID) :: eNev
    integer :: number
    logical :: NumbersAlreadySet
    logical,save :: MCmode

    type(histogram), save :: hSigmaMC_qz(0:max_Hist)
    type(histogram), save :: hSigmaMC_nu(0:max_Hist)
    type(histogram), save :: hSigmaMC_Q2(0:max_Hist)
    type(histogram), save :: hSigmaMC_X(0:max_Hist)
    type(histogram), save :: hSigmaMC_Xrec(0:max_Hist)
    type(histogram2D), save :: hSigmaMC_nuQ2(0:max_Hist)
    type(histogram2D), save :: hSigmaMC_EprimeCost(0:max_Hist)
    type(histogram2D), save :: hSigmaMC_XY(0:max_Hist)
!
!  There are 3 different definitions for the invariant mass W:
!  1. a boson with four-momentum q hits Fermi-moving nucleon in potential well: W
!  2. a boson with four-momentum q hits Fermi-moving nucleon without potential: W_free
!  3. a boson with four-momentum q hits free nucleon: Wrec
!
    type(histogram), save :: hSigmaMC_W(0:max_Hist)
    type(histogram), save :: hSigmaMC_Wrec(0:max_Hist)
    type(histogram), save :: hSigmaMC_Wfree(0:max_Hist)

    real :: Q2max,numax,W2max,cost
    integer :: iHist

    logical :: flagDUMMY

    write(*,*) '################### NEUTRINO INIT STARTS #######################'

    if (initFlag) then
       call readInput
       call DoInit
       initFlag=.false.
    end if

    if (fullEnsemble.and.realRun) then
       write(*,*) 'FullEnsemble+Real particles in final state not yet implemented!!! STOP'
       stop
    end if


    numberofcalls=numberofcalls+1
    if (numberofcalls.gt.num_runs_sameEnergy) then
       numberofcalls=1
       sigabsArrFinal = 0.
       numberofsuccessfinal=0
    end if

    sigabsArr = 0.
    numberofsuccess=0
    countsigmacutoff=0
    countfluxcutoff=0
    failuresV=0
    failuresO=0
    first=0

    !loop to determine numtry (number of testteilchen)
    numtry=0
    do i = lbound(realParticles,dim=1),ubound(realParticles,dim=1)
       do j = lbound(realParticles,dim=2),ubound(realParticles,dim=2)
          if (realParticles(i,j)%ID.ne.nucleon) cycle
          numtry=numtry+1
       end do
    end do

    raiseflag = raiseflagin
    call neutrinoProdInfo_Init(numtry)

    call SetNuclearPDFA(targetNuc%mass)

    select case (nuExp)
    case (1)
       if (nuXsectionMode.eq.EXP_dSigmadQs&
            & .or.nuXsectionMode.eq.EXP_dSigmadEnu&
            & .or.nuXsectionMode.eq.EXP_dSigmaMC) then
          call neutrinoInfoStorage_Init(numtry)
       end if
    case (3)
       if (nuXsectionMode.eq.EXP_dSigmadEnu &
            & .or. nuXsectionMode.eq.EXP_dSigmaMC) then
          call neutrinoInfoStorage_Init(numtry)
       end if
    end select

    ! set the overall kinematics (most of it as dummy):
    call eNev_SetProcess(eNev0, process_ID,flavor_ID)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!  ENSEMBLE LOOP  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Now a loop over all ensembles starts
! realParticles: in runs with frozen nuclear configuration these are the
! target nucleons, 1st index gives ensemble, 2nd index gives particle
! identity, p or n
!

    loopOverEnsemble: do i = lbound(realParticles,dim=1),ubound(realParticles,dim=1)
       if (realRun) then
          numNucleons=0
          do j = lbound(realParticles,dim=2),ubound(realParticles,dim=2)
             if (realParticles(i,j)%ID.ne.nucleon) cycle
             numNucleons=numNucleons+1
          end do
          if (numNucleons.gt.0) then
             ! Choose randomly one nucleon in the ensemble to make
             ! the collsion with:
             whichReal_nucleon=1+int(rn()*numNucleons)
          else
             cycle loopOverEnsemble
          end if
          numNucleons=0
       end if


       !print out status info
       if (mod(i,50)==0) write(*,*) 'now starting ensemble ',i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!   BOUND NUCLEON LOOP   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Now a loop over all bound nucleons in the target nucleus starts
!
       loopVector: do j = lbound(realParticles,dim=2),ubound(realParticles,dim=2)

          if (realParticles(i,j)%ID.ne.nucleon) cycle
          if (realRun) then
             ! Only allow for a collision with the chosen nucleon.
             numNucleons=numNucleons+1
             if (.not.(numNucleons.eq.whichReal_nucleon)) cycle
          end if


          flux_enu = -99.9
          if (nuExp.gt.0) then
             select case (nuExp)
             case (1)
                flux_enu=MiniBooNEenergy()
             case (2)
                flux_enu=ANLenergy()
             case (3)
                flux_enu=K2Kenergy()
             case (4)
                flux_enu=BNLenergy()
             case (5)
                flux_enu=MiniBooNEenergyBARNU()
             case (6)
                flux_enu=MINOSenergyNU_fluxNU()
             case (7)
                flux_enu=MINOSenergyBARNU_fluxNU()
             case (8)
                flux_enu=NOVAenergy_FD(Flavor_ID,Process_ID)
             case (9)
                flux_enu=T2Kenergy_ND(Flavor_ID,Process_ID)
             case (10)
                flux_enu=uniformFlux()
             case (11)
                flux_enu=MINOSenergyNU_fluxBARNU()
             case (12)
                flux_enu=MINOSenergyBARNU_fluxBARNU()
             case (13)
                flux_enu=MINERVAenergyNU()
             case (14)
                flux_enu=MINERVAenergyBARNU()
             case (15)
                flux_enu=DUNEenergyNU()
             case (16)
                flux_enu=DUNEenergyBARNU()
             case (17)
                flux_enu=LBNOenergyNU()
             case (18)
                flux_enu=NOMADenergyNU()
             case (19)
                flux_enu=BNBenergyNUe()
             case (20)
                flux_enu=BNBenergyNUebar()
             case (21)
                flux_enu=BNBenergyNUmu()
             case (22)
                flux_enu=BNBenergyNUmubar()
             case (23)
                flux_enu=NOVAenergy_ND(Flavor_ID,Process_ID)
             case (24)
                flux_enu=T2Konaxisenergy(Flavor_ID,Process_ID)
             case (25)
                flux_enu=MINERvAenergy(Flavor_ID,Process_ID)
             case default
                write(*,*) 'Experiment does not exist'
                stop
             end select

             if (flux_enu.lt.Enu_lower_cut .or. flux_enu.gt.Enu_upper_cut) then
                countfluxcutoff=countfluxcutoff+1
                cycle
             end if

             call AddHist(energyInit,flux_enu,1.)
          end if

!
! now determine target nucleon on which primary interaction takes place
! Step 1 of "neutrino init": Set the target nucleon.
!
          eNev1 = eNev0
          call eNev_init_nuStep1(eNev1,realParticles(i,j))

          sigma=0.
!
! now compute cross section in MC mode
! SetXsecMC in module neutrinoXsection
!
          if (MCmode) call SetXsecMC(eNev1,flux_enu,nuXsectionMode)

          !          call write_electronNucleon_event(eNev1)

          call resetNumberGuess()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!   LOOP OVER REACTION MECHANISMS and FINAL STATES   !!!!!!!!!!!!!!

!calculate cross sections

!max_finalstate_ID set earlier in this file
!
!first all cross sections are calculated, then a MC decision is made
!which process takes place and is further propagated


          particle_ID_loop: do k=1, max_finalstate_ID

             if (.not.includeK(k)) cycle

             eNev(k) = eNev1

             select case (k)
             case (DIS_CH)
                pOutPart => OutPartDIS(:)
             case (twoPion)
                pOutPart => OutPart2pi(:)
             case default
                pOutPart => OutPart(k,:)
             end select


             select case (nuXsectionMode)

             case (integratedSigma)
                call Xsec_integratedSigma(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

             case (dSigmadCosThetadElepton)
                call Xsec_dSigmadCosThetadElepton(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))
! if(k == 35) write(*,*) 'sigma(35) in InitNeutrino 1 =', sigma(k)

             case (dSigmadQsdElepton)
                call Xsec_dSigmadQsdElepton(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

             case (dSigmadQs)
                call Xsec_SigmaMC_Qs(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

             case (dSigmadCosTheta)
                call Xsec_dSigmadCosTheta(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

             case (dSigmadElepton)
                call Xsec_dSigmadElepton(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

             case (dSigmaMC)
                call Xsec_SigmaMC(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

             case (dSigmadW)
                call Xsec_SigmaMC_W(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

                !experimental cross sections

             case (EXP_dSigmadEnu)
                call Xsec_integratedSigma(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)

             case (EXP_dSigmadCosThetadElepton)
                call Xsec_dSigmadCosThetadElepton(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)

             case (EXP_dSigmadQsdElepton)
                call Xsec_dSigmadQsdElepton(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)

             case (EXP_dSigmadQs)
                call Xsec_SigmaMC_Qs(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)

             case (EXP_dSigmadCosTheta)
                call Xsec_dSigmadCosTheta(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)

             case (EXP_dSigmadElepton)
                call Xsec_dSigmadElepton(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)

             case (EXP_dSigmaMC)
                call Xsec_SigmaMC(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)

             case (EXP_dSigmadW)
                call Xsec_SigmaMC_W(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)


             case default
                write(*,*) 'nuXsectionMode=',nuXsectionMode
                call TRACEBACK('error in case nuXsectionMode')

             end select

             raiseflag=.false.

             ! Pauli blocking
             if (sigma(k).ne.0.) then
                if (.not.checkPauli(pOutPart,realParticles)) sigma(k)=0.
             end if

          end do particle_ID_loop

          if (realRun) sigma=sigma*real(size(realParticles,2))/11.
! number 11 comes from main.f90 where lengthReal is
! set to lengthReal=targetNuc%mass+targetNuc%mass*10

          sigtot=sum(sigma)
          if (debugflag) write(*,'(A,g13.5)') ' sigtot= ', sigtot

!!$   write(*,'(A,g13.5)') ' sigtot= ', sigtot
!!$   write(*,'(2g13.5)') sigma(1:2)

!     check if total cross section smaller than fixed cutoff,
!     if so, then throw cross section away

          if (sigtot.lt.sigmacut) then
             countsigmacutoff=countsigmacutoff+1
             if (debugflag) write(*,'(A,i7.3)') ' countsigmacutoff=',&
           & countsigmacutoff
             if (debugflag) write(*,'(2(A,g12.5))') &
           & 'In initNeutrino.f90: sigtot=',sigtot,' is less than sigmacut=',&
           & sigmacut
             cycle loopVector
          end if

          firstEvent=getFirstEvent()

! event selection:
          k = MonteCarloChoose(sigma,totalWeight)
! totalWeight=Perweight which should be assigned to the chosen channel

          select case (k)
          case (0)
             write(*,*) 'Problem initNeutrino: no event generated:', sigma
             write(*,*) 'Stop'
             stop
          case (DIS_CH)
             pOutPart => OutPartDIS(:)
          case (twoPion)
             pOutPart => OutPart2pi(:)
          case default
             pOutPart => OutPart(k,:)
          end select

          allocate(finalstate(size(pOutPart)))
          finalstate = pOutPart

          finalState%perturbative=(.not.realRun)
          finalState%firstEvent=FirstEvent           ! Number of first event,
                                                     ! stays constant during run.
          finalState%perweight=totalWeight/float(numtry)
          finalState%history=0

          numberofsuccess(k)=numberofsuccess(k)+1

          ! all what's to do with xsections:

          fak1 = sigtot/float(numtry)
          fak2 = 1./float(num_runs_sameEnergy)

          call neutrinoProdInfo_Store(firstEvent, k, fak1,&
               eNev(k)%lepton_in%momentum,&
               eNev(k)%lepton_out%momentum,&
               eNev(k)%boson%momentum,&
               eNev(k)%nucleon_free%momentum,&
               eNev(k)%nucleon_free%charge)

          iHist = K2Hist(k)

          call AddHist2D(hSigmaMC_nuQ2(0),hSigmaMC_nuQ2(iHist),&
               & (/eNev(k)%boson%momentum(0),eNev(k)%Qsquared/),&
               & fak1)
          call AddHist(hSigmaMC_nu(0),hSigmaMC_nu(iHist),&
               & eNev(k)%boson%momentum(0),&
               & fak1)
          call AddHist(hSigmaMC_qz(0),hSigmaMC_qz(iHist),&
               & sqrt(eNev(k)%boson%momentum(0)**2+eNev(k)%Qsquared),&
               & fak1)
          call AddHist(hSigmaMC_Q2(0),hSigmaMC_Q2(iHist),&
               & eNev(k)%Qsquared,&
               & fak1)
          call AddHist(hSigmaMC_Xrec(0),hSigmaMC_Xrec(iHist),&
               & eNev(k)%Qsquared/(2*mN*eNev(k)%boson%momentum(0)),&
               & fak1)

          if (eNev(k)%W > 0) then
           call AddHist(hSigmaMC_W(0),hSigmaMC_W(iHist),&
                  & eNev(k)%W,&
                  & fak1)
          end if

          if (eNev(k)%W_free > 0) then
  !        Write(*,*) 'Wfree', eNev(k)%W_free
             call AddHist(hSigmaMC_Wfree(0),hSigmaMC_Wfree(iHist),&
                  & eNev(k)%W_free,&
                  & fak1)
          end if

          if (eNev(k)%W_rec > 0) then
             call AddHist(hSigmaMC_Wrec(0),hSigmaMC_Wrec(iHist),&
                  & eNev(k)%W_rec,&
                  & fak1)
          end if


          cost = eNeV_Get_CostLepton(eNev(k))

          call AddHist2D(hSigmaMC_EprimeCost(0),hSigmaMC_EprimeCost(iHist),&
               & (/eNev(k)%lepton_out%momentum(0),cost/),&
               & fak1)

          if (.not.MCmode) then
             call CalcXY(eNev(k)%lepton_in%momentum, &
                  & eNev(k)%lepton_out%momentum,&
                  & eNev(k)%nucleon%momentum, MC_x,MC_y, flagDUMMY)
          end if

          call AddHist2D(hSigmaMC_XY(0),hSigmaMC_XY(iHist),&
               & (/MC_x,MC_y/),&
               & fak1)
          call AddHist(hSigmaMC_X(0),hSigmaMC_X(iHist),&
               & MC_x,&
               & fak1)


          sigabsArr(0)  = sigabsArr(0)  + fak1
          sigabsArr(1:) = sigabsArr(1:) + sigma(1:)/float(numtry)



          ! calculate the finalstate offshellness

          select case (k)
          case default ! ==== single-pi, DIS, 2p2h, two-pi-backg

             finalstate%offshellparameter=0.
             call updateVelocity(finalstate)

          case (:31)    ! ==== RES or QE

             finalstate%offshellparameter=0.
             if (get_useOffShellPotentialBaryons()) then
                finalstate%offshellparameter=getOffShellParameter(k,   &
                     & finalstate(1)%mass,finalstate(1)%momentum, &
                     & RealParticles(i,j)%position,success)
                if (.not.success) then
                   if (debugflag) write(*,*) 'offshell parameter > max_offshellparameter &
                                            & in initNeutrino => sig=0'
                   failuresO=failuresO+1
                   deallocate(finalstate)
                   cycle loopVector
                end if
             end if

             call gradients(finalstate(1),grad_P) ! Evaluate dH/dp
             finalstate(1)%velocity=grad_P
             if (1. - Dot_Product(grad_P(1:3),grad_P(1:3)) .le. 0.) then
                write(*,'(A,5G13.5)')'problems in initNeutrino: &
                     & velocity**2 greater or equal 1.', &
                     & k, Dot_Product(grad_P(1:3) &
                     & ,grad_P(1:3)), finalstate(1)%mass,&
                     & finalstate(1)%offshellparameter,sigma(k)/float(numtry)
                if (debugflag) then
                   call WriteParticle_debug(finalstate(1))
                   write(*,*)
                   write(*,*) '...this funny particle is now deleted'
                end if
                failuresV=failuresV+1
                deallocate(finalstate)
                cycle loopVector
             end if

          end select

          ! now we can give them the first-event-number:

          if (.not.checkEvent(eNev(k),finalstate,eOrigin(k))) then
             call TRACEBACK('conservations violated.')
          end if


          if (realRun) then
             number = real_numbering()
          else
             number = pert_numbering(realParticles(i,j))
          end if
          finalState%event(1) = number
          finalState%event(2) = number

          ! set the particles in the particle vector:

          if (k.ne.DIS_CH) call resetNumberGuess()
          NumbersAlreadySet = AcceptGuessedNumbers()

          if (fullEnsemble) then
             if (realRun) stop 'RealRun+fullEnsemble not yet implemented: initNeutrino'
             call setIntoVector(finalState,pertParticles,setFlag,NumbersAlreadySet)
          else
             if (realRun) then

                stop 'RealRun seems to be corrupted: initNeutrino'

!!! ATTENTION !!! The following seems to be untested !!!!!!!!

                ! Delete the scattering partner
                realParticles(i,j)%ID=0
                ! Set new particles into real particle vector
                realParticles(i:i,:)%perweight=finalState(1)%perweight      ! ????????
                realParticles(i:i,:)%firstEvent=finalState(1)%firstEvent    ! ????????
                call setIntoVector(finalState,realParticles(i:i,:),setFlag,NumbersAlreadySet)
                ! ????????
             else
                call setIntoVector(finalState,pertParticles(i:i,:),setFlag,NumbersAlreadySet)
             end if
          end if

          deallocate(finalstate)

          if (.not.setFlag) then
             write(*,*) 'error setIntoVector in initNeutrino'
             if (.not.realRun) then
                call WriteParticleVector('pert',pertParticles)
             else
                call WriteParticleVector('real',realParticles)
             end if
             stop
          end if

       end do loopVector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!   END OF BOUND NUCLEON LOOP   !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end do loopOverEnsemble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!   END OF LOOPS OVER ENSEMBLE   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






    call SetNuclearPDFA(1)

    write(*,'(40("-"))')
    write(*,'(a,I10)') ' failures v>c            :',failuresV
    write(*,'(a,I10)') ' failures due to offshell:',failuresO
    write(*,'(a,I10)') ' numtry                  =',numtry
    write(*,'(a,I10)') ' numberofsuccess         =',sum(numberofsuccess)
    write(*,'(a,I10)') ' countsigmacutoff        =', countsigmacutoff
    write(*,'(a,I10)') ' countfluxcutoff         =', countfluxcutoff
    write(*,'(a,F12.4)')         ' sigabs       =', sigabsArr(0)
    write(*,'(a,F12.4,a,I10,a)') ' sigabsQE     =', sigabsArr(1), '(#',numberofsuccess(1),')'
    write(*,'(a,F12.4,a,I10,a)') ' sigabsDelta  =', sigabsArr(2), '(#',numberofsuccess(2),')'
    write(*,'(a,F12.4,a,I10,a)') ' sigabshighRes=', sum(sigabsArr(3:31)), &
                                                     & '(#',sum(numberofsuccess(3:31)), ')'
    write(*,'(a,F12.4,a,I10,a)') ' sigabs1pi (n)=', sigabsArr(32), '(#',numberofsuccess(32), ')'
    write(*,'(a,F12.4,a,I10,a)') ' sigabs1pi (p)=', sigabsArr(33), '(#',numberofsuccess(33), ')'
    write(*,'(a,F12.4,a,I10,a)') ' sigabs2pi    =', sigabsArr(37), '(#',numberofsuccess(37), ')'
    write(*,'(a,F12.4,a,I10,a)') ' sigabsDIS    =', sigabsArr(34), '(#',numberofsuccess(34),')'
    write(*,'(a,F12.4,a,I10,a)') ' sigabs2p2h   =', sum(sigabsArr(35:36)), &
                                                     & '(#',sum(numberofsuccess(35:36)),')'
    write(*,'(40("-"))')

    sigabsArrFinal = sigabsArrFinal+sigabsArr
    numberofsuccessfinal=numberofsuccess+numberofsuccessfinal

    call DoWrite

    !**************************************************************************
    !****o* initNeutrino/neutrino_initialized_energyFlux.dat
    ! NAME
    ! file neutrino_initialized_energyFlux.dat
    ! PURPOSE
    ! This file provides the incoming neutrino flux
    ! Only energies at which actual sampling has taken place are in this file
    !
    ! Columns:
    ! * #1: Enu (in GeV)
    ! * #2: flux (in 1/GeV), for the calculations of X-sections relevant is only
    !   the normalized flux
    ! * #3: number of events
    ! * #4: should be ignored
    !**************************************************************************
    if (nuExp.gt.0) then
       call writeHist(energyInit,10, &
            & file='neutrino_initialized_energyFlux.dat')
    end if

    if (printAbsorptionXS) then

    !**************************************************************************
    !****o* initNeutrino/neutrino.NuQ2planeXS.ZZZ.dat
    ! NAME
    ! file neutrino.NuQ2planeXS.ZZZ.dat,   ZZZ=000 - 008
    ! PURPOSE
    ! This files provides absorption cross section in the
    ! nu(transferred energy)-Q^2 plane.
    ! Format is suitable as input to gnuplots "splot" command
    !
    ! The units are 10^{-38} cm^2/GeV^3 for process_ID=CC, NC and
    ! 10^{-33} cm^2/GeV^3  for EM
    !
    ! ZZZ = ...:
    ! * 000: sum over all channels
    ! * 001: QE cross section
    ! * 002: Delta
    ! * 003: highRES
    ! * 004: 1pi
    ! * 005: DIS
    ! * 006: 2p2h QE
    ! * 007: 2p2h Delta
    ! * 008: 2pi background
    !
    ! Columns:
    ! * #1: nu (transfered energy),  GeV
    ! * #2: Q^2, GeV^2
    ! * #3: cross section
    ! * #4: number of entries that lead to cross section
    ! * #5: should be ignored
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.Q2NuplaneXS.ZZZ.dat
    ! NAME
    ! file neutrino.Q2NuplaneXS.ZZZ.dat
    ! PURPOSE
    ! Nearly the same as  neutrino.NuQ2planeXS.ZZZ.dat, but:
    ! * #1: Q^2, GeV^2
    ! * #2: nu,  GeV
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.EprimeCostplaneXS.ZZZ.dat
    ! NAME
    ! file neutrino.EprimeCostplaneXS.ZZZ.dat,   ZZZ=000 - 008
    ! PURPOSE
    ! This files provides absorption cross section in the
    ! Eprime(energy of the outgoing lepton)-costheta(cos of the angle between
    ! the incoming and outgoing leptons) plane.
    ! Format is suitable as input to gnuplots "splot" command
    !
    ! The units are 10^{-38} cm^2/GeV for process_ID=CC, NC and
    ! 10^{-33} cm^2/GeV  for EM
    !
    ! ZZZ is the origin (the first interaction vertex) of the event:
    ! * 000: sum over all origins
    ! * 001: QE cross section
    ! * 002: Delta
    ! * 003: highRES
    ! * 004: 1pi
    ! * 005: DIS
    ! * 006: 2p2h QE
    ! * 007: 2p2h Delta
    ! * 008: 2pi background
    !
    ! Columns:
    ! * #1: Eprime (energy of the outgoing lepton),  GeV
    ! * #2: costheta(cos of the angle between the incoming and outgoing leptons)
    ! * #3: cross section
    ! * #4: number of entries that lead to cross section
    ! * #5: should be ignored
    !**************************************************************************


    !**************************************************************************
    !****o* initNeutrino/neutrino.XYplaneXS.ZZZ.dat
    ! NAME
    ! file neutrino.XYplaneXS.ZZZ.dat
    ! ZZZ=000 - 008  is the origin (the first interaction vertex) of the event:
    !                (see description in  neutrino.EprimeCostplaneXS.ZZZ.dat)
    ! PURPOSE
    ! This files provides absorption cross section in the
    ! x_{Bjorken}-y(=nu/Enu) plane.
    ! Format is suitable as input to gnuplots "splot" command
    !
    ! The units are 10^{-38} cm^2 for process_ID=CC, NC and
    ! 10^{-33} cm^2  for EM
    !
    ! Columns:
    ! * #1: x_Bjorken
    ! * #2: y (=nu/Enu)
    ! * #3: cross section
    ! * #4: number of entries that lead to cross section
    ! * #5: should be ignored
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.NuXS.ZZZ.dat
    ! NAME
    ! file neutrino.NuXS.ZZZ.dat
    ! ZZZ=000 - 008  is the origin (the first interaction vertex) of the event:
    !                (see description in  neutrino.EprimeCostplaneXS.ZZZ.dat)
    ! PURPOSE
    ! This files provides absorption cross section versus nu(=transfered energy)
    ! Format is suitable as input to gnuplots "splot" command
    !
    ! The units are 10^{-38} cm^2/GeV for process_ID=CC, NC and
    ! 10^{-33} cm^2/GeV  for EM
    !
    ! Columns:
    ! * #1: nu (transfered energy), GeV
    ! * #2: cross section
    ! * #3: number of entries that lead to cross section
    ! * #4: should be ignored
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.Q2XS.ZZZ.dat
    ! NAME
    ! file neutrino.NuQ2.ZZZ.dat
    ! PURPOSE
    ! Similar to neutrino.NuXS.ZZZ.dat, but versus Q2
    !
    ! Coulmns:
    ! * #1: Q2 , GeV^2
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.XXS.ZZZ.dat
    ! NAME
    ! file neutrino.XXS.ZZZ.dat
    ! PURPOSE
    ! Similar to neutrino.NuXS.ZZZ.dat, but versus x_Bjorken
    !
    ! Columns:
    ! * #1: x_Bjorken
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.XrecXS.ZZZ.dat
    ! NAME
    ! file neutrino.XrecXS.ZZZ.dat
    ! PURPOSE
    ! Similar to neutrino.NuXS.ZZZ.dat, but versus reconstructed x_Bjorken,
    ! defined for nucleon at rest: Xrec = Q^2/(2*mN*nu)
    !
    ! Columns:
    ! * #1: free x_Bjorken
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.WrecXS.ZZZ.dat
    ! NAME
    ! file neutrino.WrecXS.ZZZ.dat
    ! PURPOSE
    ! Similar to neutrino.NuXS.ZZZ.dat, but versus Wrec=qsrt(mN^2+2*mN*nu-Q2) ,
    ! where mN=0.938 GeV = nucleon mass (see constants.f90), i.e. W
    ! reconstructed for nucleon at rest
    !
    ! Columns:
    ! * #1: Wrec , GeV
    !**************************************************************************


    !**************************************************************************
    !****o* initNeutrino/neutrino.WfreeXS.ZZZ.dat
    ! NAME
    ! file neutrino.WfreeXS.ZZZ.dat
    ! PURPOSE
    ! Similar to neutrino.NuXS.ZZZ.dat, but versus Wfree, i.e. W value for
    ! boson with four-momentum q and Fermi-moving nucleon, without potential
    !
    ! Columns:
    ! * #1: Wfree , GeV
    !**************************************************************************


    !**************************************************************************
    !****o* initNeutrino/neutrino.WXS.ZZZ.dat
    ! NAME
    ! file neutrino.WXS.ZZZ.dat
    ! PURPOSE
    ! Similar to neutrino.NuXS.ZZZ.dat, but versus W, i.e. W value for
    ! boson with four-momentum q and Fermi-moving nucleon in the potential
    !
    ! Columns:
    ! * #1: Wfree , GeV
    !**************************************************************************

    do iHist=0,max_Hist
       if (.not.includeHist(iHist)) cycle
       call WriteHist2D_Gnuplot(hSigmaMC_nuQ2(iHist),10, &
            & mul=1.0/numberofcalls, add=1e-20,&
            & file='neutrino.NuQ2planeXS.'//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)
       call WriteHist2D_Gnuplot(hSigmaMC_nuQ2(iHist),10, &
            & mul=1.0/numberofcalls, add=1e-20,&
            & file='neutrino.Q2NuplaneXS.'//trim(intToChar(iHist))//'.dat',&
            & dump=.false., SwapXY=.true.)
       call WriteHist2D_Gnuplot(hSigmaMC_EprimeCost(iHist),10, &
            & mul=1.0/numberofcalls, add=1e-20,&
            & file='neutrino.EprimeCostplaneXS.'//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)
       call WriteHist2D_Gnuplot(hSigmaMC_XY(iHist),10, &
            & mul=1.0/numberofcalls, add=1e-20,&
            & file='neutrino.XYplaneXS.'//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)
       call WriteHist(hSigmaMC_nu(iHist),10, &
            & mul=1.0/numberofcalls, add=1e-20,&
            & file='neutrino.NuXS.'//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)
       call WriteHist(hSigmaMC_Q2(iHist),10, &
            & mul=1.0/numberofcalls, add=1e-20,&
            & file='neutrino.Q2XS.'//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)
       call WriteHist(hSigmaMC_qz(iHist),10, &
            & mul=1.0/numberofcalls, add=1e-20,&
            & file='neutrino.qzXS.'//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)
       call WriteHist(hSigmaMC_X(iHist),10, &
            & mul=1.0/numberofcalls, add=1e-20,&
            & file='neutrino.XXS.'//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)
       call WriteHist(hSigmaMC_Xrec(iHist),10, &
            & mul=1.0/numberofcalls, add=1e-20,&
            & file='neutrino.XrecXS.'//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)
       call WriteHist(hSigmaMC_Wrec(iHist),10, &
            & mul=1.0/numberofcalls, add=1e-20,&
            & file='neutrino.WrecXS.'//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)
       call WriteHist(hSigmaMC_Wfree(iHist),10, &
            & mul=1.0/numberofcalls, add=1e-20,&
            & file='neutrino.WfreeXS.'//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)
       call WriteHist(hSigmaMC_W(iHist),10, &
            & mul=1.0/numberofcalls, add=1e-20,&
            & file='neutrino.WXS.'//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)
    end do

    end if

    write(*,*) '################### NEUTRINO INIT FINISHED #######################'

  contains
    !**************************************************************************
    !****s* init_neutrino/DoInit
    ! NAME
    ! subroutine DoInit
    ! PURPOSE
    ! Doing the actual initialzation steps
    !**************************************************************************
    subroutine DoInit()

      ! numbers given in array nuMaxArr are upper boundaries for energy-transfer
      real, dimension(1:numberOfExperiments), parameter :: nuMaxArr = (/&
           & 2.5, 3.0, 4.0, 5.0, 2.5, & ! MiniBooNE-nu,ANL,K2K,BNL,MiniBooNE-barnu
           & 30.0, 30.0, 15.0, 3.5, &   ! MINOS-nu-numode, MINOS-barnu-numode, NOvA, T2K OA2.5,
           & 3.0, 30.0,  &              ! uniform distr, MINOS-nu-barnumode,
           & 30.0,      &               ! MINOS-barnu-barnumode
           & 30.0, 30.0, 30.0, 30.0, &  ! MINERvA-numu, MINERvA-antinumu,DUNE-nu,DUNE-barnu
           & 30.0, 300.0,      &        ! LBNO-nu, NOMAD
           & 7.5, 7.5, 7.5, 7.5, &      ! BNB-nue,BNB-nuebar,BNBnumu,BNBnumubar
           & 15., 20., 20. /)           ! NOvA, T2K, MINERvA 2016


      !---------------------------------------------------------------------
      ! setting up some arrays for switching on/off the channels
      !---------------------------------------------------------------------
      includeHist = (/.true.,includeQE, includeDELTA, includeRES, &
           include1pi, includeDIS, include2p2hQE, include2p2hDelta, &
           include2pi /)

      includeK = .false.
      do k=1,max_finalstate_ID

         select case (k) ! === check for inclusion of process
         case (nucleon)
            if (.not.includeQE) cycle
         case (delta)
            if (.not.includeDELTA) cycle
         case (P11_1440:F37_1950)
            if (.not.includeRES) cycle
         case (onePionCH_n,onePionCH_p)
            if (.not.include1pi) cycle
         case (DIS_CH)
            if (.not.includeDIS) cycle
         case (QE2p2h)
            if (.not.include2p2hQE) cycle
         case (delta2p2h)
            if (.not.include2p2hDelta) cycle
         case (twoPion)
            if (.not.include2pi) cycle
         end select

         select case (k) ! === check for inclusion of resonance
         case (S11_2090,D13_2080,G17_2190,P11_2100,P13_1900,F15_2000,S31_1900,D33_1940, &
               & D35_1930,D35_2350,P31_1750,F35_1750)
            cycle
         end select

         includeK(k) = .true.
      end do

      !---------------------------------------------------------------------
      ! initialize some output files
      !---------------------------------------------------------------------
      open(10,File='neutrino_info.dat')
      if (process_ID.lt.0) then
         write(10,*)'# process_ID : anti-',sProcess(-process_ID)
      else
         write(10,*)'# process_ID : ',sProcess(process_ID)
      end if
      write(10,*)'# flavor_ID  : ',sFamily(flavor_ID)
      write(10,*)'# '
      write(10,'(a,L2)')' #...include QE        : ',includeQE
      write(10,'(a,L2)')' #...include Delta     : ',includeDELTA
      write(10,'(a,L2)')' #...include higher RES: ',includeRES
      write(10,'(a,L2)')' #...include 1pi       : ',include1pi
      write(10,'(a,L2)')' #...include 2pi       : ',include2pi
      write(10,'(a,L2)')' #...include DIS       : ',includeDIS
      write(10,'(a,L2)')' #...include 2p2h QE   : ',include2p2hQE
      write(10,'(a,L2)')' #...include 2p2h Delta: ',include2p2hDelta
      write(10,*)'# '
      write(10,*)'# nuXsectionMode: ',sXsectionMode(MOD(nuXsectionMode,10))
      write(10,*)'# '
      if (nuExp.gt.0) write(10,*)'##### calculation is done for the ',&
           & trim(sExp(nuExp)),' experiment #####'
      close(10)

      if (printAbsorptionXS) then

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section.dat
      ! NAME
      ! file neutrino_absorption_cross_section.dat
      ! PURPOSE
      ! The file is produced in the runs with eventtype=5=neutrino .
      !
      ! The file shows the absorption cross section for lepton
      ! ( neutrino or charged lepton) scattering for the sum of all channels
      ! which were set to TRUE in the namelist "neutrino_induced"
      ! ( QE+Delta+highRES+1pi+DIS if includeQE, includeDELTA, includeRES,
      ! include1pi, includeDIS were TRUE)
      !
      ! For process_ID=CC and NC  the units 10^{-38} cm^2 for integrated xsec
      ! (10^{-38)cm^2/GeV  for  dsigma/dElepton,  10^{-38)cm^2/GeV^2  for
      ! dsigma/dQ^2, and so on)                                               =the
      ! For process_ID=EM the units are nanobars=10^{-33}cm^2
      !
      ! Columns:
      ! * #1: variable which was raised
      !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode, Elepton for
      !   nuXsectionMode=2=dSigmadQsdElepton  and so on)
      ! * #2: cross section
      !************************************************************************
      open(10,File='neutrino_absorption_cross_section.dat')
      write(10,*) '#  raiseVal, xsection'
      close(10)

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_ALL.dat
      ! NAME
      ! file neutrino_absorption_cross_section_ALL.dat
      ! PURPOSE
      ! More detailed information than in
      ! neutrino_absorption_cross_section.dat:
      !
      ! Columns:
      ! * #1: variable which was raised (the same as
      !   neutrino_absorption_cross_section.dat)
      !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode, Elepton for
      !   nuXsectionMode=2=dSigmadQsdElepton )
      ! * #2: cross section, sum over all channels  (the same as
      !   neutrino_absorption_cross_section.dat)
      ! * #3: cross section for QE events (the same as column 2 in
      !   neutrino_absorption_cross_section_QE.dat)
      ! * #4: cross section for DELTA events (the same as column 2 in
      !   neutrino_absorption_cross_section_Delta.dat)
      ! * #5: cross section for hihgRES events (the same as column 2 in
      !   neutrino_absorption_cross_section_highRES.dat)
      ! * #6: cross section for 1pi events (the same as column2+column3 in
      !   neutrino_absorption_cross_section_1pi.dat)
      ! * #7: cross section for DIS events (the same as column 2 in
      !   neutrino_absorption_cross_section_DIS.dat)
      ! * #8: cross section for 2p2h QE events
      ! * #9: cross section for 2p2h Delta events
      ! * #10: cross section for 2 pion background events
      !
      !************************************************************************
      open(10,File='neutrino_absorption_cross_section_ALL.dat')
      write(10,*) '# 1:var 2:sum 3:QE 4:Delta 5:highRES 6:1pi 7:DIS 8:2p2h-QE &
                  & 9:2p2h-Delta 10:2pi'
      close(10)

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_QE.dat
      ! NAME
      ! file neutrino_absorption_cross_section_QE.dat
      ! PURPOSE
      ! The same structure as neutrino_absorption_cross_section.dat,
      ! but only for QE events (=the first interaction act was quasielastic
      ! or elastis scattering)
      !************************************************************************
      if (includeQE) then
         open(10,File='neutrino_absorption_cross_section_QE.dat')
         write(10,*) '#  raiseVal, xsection'
         close(10)
      end if

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_Delta.dat
      ! NAME
      ! file neutrino_absorption_cross_section_Delta.dat
      ! PURPOSE
      ! The Delta production events (=the first interaction was
      ! production of the Delta resonance)
      !************************************************************************
      if (includeDELTA) then
         open(10,File='neutrino_absorption_cross_section_Delta.dat')
         write(10,*) '#  raiseVal, xsection'
         close(10)
      end if

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_highRES.dat
      ! NAME
      ! file neutrino_absorption_cross_section_highRES.dat
      ! PURPOSE
      ! For highRES production events (=the first interaction was
      ! production of any resonance beyond Delta)
      !
      ! Columns:
      ! * #1: variable which was raised (the same as
      !   neutrino_absorption_cross_section.dat)
      !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode, Elepton for
      !   nuXsectionMode=2=dSigmadQsdElepton )
      ! * #2: cross section, sum over all higher resonances beyond the Delta
      ! * #3 - 31: contribution of individual nucleon resonances beyond the Delta
      !   Individual resonance numbers in Module IdTable
      !   3: P11(1440), 4: S11(1535), ..., 7: D13(1520), ... etc
      !
      !************************************************************************
      if (includeRES) then
         open(10,File='neutrino_absorption_cross_section_highRES.dat')
         write(10,*) '#  raiseVal, sum, single contribution xsection'
         close(10)
      end if

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_1pi.dat
      ! NAME
      ! file neutrino_absorption_cross_section_1pi.dat
      ! PURPOSE
      ! Nearly the same structure as neutrino_absorption_cross_section.dat,
      ! but only for nonresonant 1-pion production events (=the first
      ! interaction was production of 1-pion final state)
      !
      ! Columns:
      ! * #2: cross section  for the channel with neutron in the final state
      !   ("final" here is "after the first interaction act", that is before
      !   final state interactions,
      !   e.g. nu n \to mu- n pi^0)
      ! * #3: cross section  for the channel with proton in the final state,
      !   e.g. nu n \to mu- p pi^-)
      !************************************************************************
      if (include1pi) then
         open(10,File='neutrino_absorption_cross_section_1pi.dat')
         write(10,*) '#  raiseVal, pi neutron-backgr, pi proton-backgr'
         close(10)
      end if

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_DIS.dat
      ! NAME
      ! file neutrino_absorption_cross_section_DIS.dat
      ! PURPOSE
      ! The same structure as neutrino_absorption_cross_section.dat,
      ! but only for DIS production events (=the first interaction was DIS)
      !************************************************************************
      if (includeDIS) then
         open(10,File='neutrino_absorption_cross_section_DIS.dat')
         write(10,*) '#  raiseVal, xsection'
         close(10)
      end if

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_numbers.dat
      ! NAME
      ! file neutrino_absorption_cross_section_numbers.dat
      ! PURPOSE
      ! Shows the number of test particles which interacted (here we mean
      ! the first interaction) via various channels
      !
      ! Columns:
      ! * #1: variable which was raised (the same as
      !   neutrino_absorption_cross_section.dat)
      !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode, Elepton for
      !   nuXsectionMode=2=dSigmadQsdElepton )
      ! * #2: QE channel
      ! * #3: Delta production
      ! * #4: P_{11}(1440) production
      ! * #5: S_{11}(1535) production
      ! * #6: D_{13}(1520) production
      ! * #7: neutron and pi+  in the final state
      ! * #8: proton and pi0 in the final state
      ! * #9: number of test particles that underwent any kind of interaction
      ! * #10: number of test particles (= number of tries)
      !************************************************************************
      open(10,File='neutrino_absorption_cross_section_numbers.dat')
      write(10,*) '#  raiseVal, numberofQE, numberofDELTA, numberP11,', &
           &' numberS11, numberD13,', &
           & '# number_npi+, number_ppi0, numberTOT, numberTestParticles'
      close(10)

      end if

      !---------------------------------------------------------------------
      ! Allocate some histograms
      !---------------------------------------------------------------------

      MCmode = ((nuXsectionMode.eq.dSigmaMC)&
           & .or.(nuXsectionMode.eq.EXP_dSigmaMC) &
           & .or.(nuXsectionMode.eq.dSigmadQs) &
           & .or.(nuXsectionMode.eq.EXP_dSigmadQs) &
           & .or.(nuXsectionMode.eq.dSigmadW) &
           & .or.(nuXsectionMode.eq.EXP_dSigmadW) )



      if (nuExp.gt.0) then
         call CreateHist(energyInit, 'initialized energy',0.,150.,0.02)
      end if

      ! this is to estimate the "reasonable" numax in the neutrino
      ! experiments in order to use it for initiallisation of histograms
      ! for differential absorption cross sections
      if (nuExp.gt.0) then
         numax = nuMaxArr(nuExp)
      else
         !  the changes are made because numax=-10. was always returned here.
         !  This  was because the namelist nl_neutrinoxsection is not read yet
         !  and the default value is used.
         call get_xsection_namelist(XsectionMode=nuXsectionMode,Genu=numax)
      end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! In these 2D histograms the last 3 2d arrays determine the binning.
! For example, in the EprimeCost histogram the input
! (/0.,-1./),(/numax,1.0/),(/0.01,0.1/) stands for 0 < Eprime < numax with a
! binwidth of 0.01 and -1 cost < +1 with a bindwidth of 0.1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) 'In initialization of the histograms  for absorption xsec:'
      write(*,*) '    numax=', numax
      numax = max(numax,1e-10)
      Q2max = 2*mN*numax+mN**2       ! rough estimate
      W2max = Q2max                  ! rough estimate
      Q2max = min(10.0,Q2max)
      do iHist=0,8
         if (.not.includeHist(iHist)) cycle
         call CreateHist2D(hSigmaMC_nuQ2(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs nu and Q2 ', &
              & (/0.,0./),(/numax,Q2max/),(/0.1,0.1/))
         call CreateHist2D(hSigmaMC_EprimeCost(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs Eprime and cost ', &
              & (/0.,-1./),(/numax,1.0/),(/0.01,0.1/))
         call CreateHist2D(hSigmaMC_XY(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs x and y ', &
              & (/0.,0./),(/2.1,2.1/),(/1.0/100,1.0/100/))
         call CreateHist(hSigmaMC_nu(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs nu ', &
              & 0.,numax,0.01)
         call CreateHist(hSigmaMC_qz(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs qz=sqrt(nu^2+Q2) ', &
              & 0.,numax,0.01)
         call CreateHist(hSigmaMC_Q2(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs Q2 ', &
              & 0.,Q2max,0.01)
         call CreateHist(hSigmaMC_X(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs X_Bjorken ', &
              & 0.,2.0,0.01)
         call CreateHist(hSigmaMC_Xrec(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs X_Bjorken reconstructed = &
              & Q2/2./mN/nu', 0.,2.0,0.01)
         call CreateHist(hSigmaMC_Wrec(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs Wrec = &
              & sqrt(mN^2 +2*mN*nu - Q2) ', &
              & 0.,sqrt(W2max),0.02)
         call CreateHist(hSigmaMC_Wfree(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs Wfree', &
              & 0.,sqrt(W2max),0.02)
         call CreateHist(hSigmaMC_W(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs W', &
              & 0.,sqrt(W2max),0.02)
      end do

    end subroutine DoInit


    !**************************************************************************
    !****s* init_neutrino/DoWrite
    ! NAME
    ! subroutine DoWrite
    ! PURPOSE
    ! Doing the write out of most data
    !**************************************************************************
    subroutine DoWrite()

      if (.not.printAbsorptionXS) return

      open(10,File='neutrino_absorption_cross_section.dat',&
           & position='append')
      if (numberofcalls.ne.1) backspace(10)
      write(10,'(10g13.5)') raiseVal, sigabsArrFinal(0)/real(numberofcalls)
      close(10)

      open(10,File='neutrino_absorption_cross_section_ALL.dat',&
           & position='append')
      if (numberofcalls.ne.1) backspace(10)
      write(10,'(30g13.5)') raiseVal, sigabsArrFinal(0)/real(numberofcalls),&
           & sigabsArrFinal(1)/real(numberofcalls), &
           & sigabsArrFinal(2)/real(numberofcalls), &
           & sum(sigabsArrFinal(3:31))/real(numberofcalls), &
           & sum(sigabsArrFinal(32:33))/real(numberofcalls), &
           & sigabsArrFinal(34)/real(numberofcalls), &
           & sigabsArrFinal(35)/real(numberofcalls), &
           & sigabsArrFinal(36)/real(numberofcalls), &
           & sigabsArrFinal(37)/real(numberofcalls)
      close(10)

!!$      write(113,'(10g13.5)') raiseVal, sigabsArrFinal(0)/real(numberofcalls),&
!!$           & sigabsArrFinal(1)/real(numberofcalls), &
!!$           & sigabsArrFinal(2)/real(numberofcalls), &
!!$           & sum(sigabsArrFinal(3:31))/real(numberofcalls), &
!!$           & sum(sigabsArrFinal(32:33))/real(numberofcalls), &
!!$           & sigabsArrFinal(34)/real(numberofcalls)

      if (includeQE) then
         open(10,File='neutrino_absorption_cross_section_QE.dat',&
              & position='append')
         if (numberofcalls.ne.1) backspace(10)
         write(10,'(10g13.5)') raiseVal, sigabsArrFinal(1)/real(numberofcalls)
         close(10)
      end if

      if (includeDELTA) then
         open(10,File='neutrino_absorption_cross_section_Delta.dat',&
              & position='append')
         if (numberofcalls.ne.1) backspace(10)
         write(10,'(10g13.5)') raiseVal, sigabsArrFinal(2)/real(numberofcalls)
         close(10)
      end if

      if (includeRES) then
         open(10,File='neutrino_absorption_cross_section_highRES.dat',&
              & position='append')
         if (numberofcalls.ne.1) backspace(10)
         write(10,'(31g13.5)') raiseVal, &
              & sum(sigabsArrFinal(3:31))/real(numberofcalls),&
              & sigabsArrFinal(3:31)/real(numberofcalls)
         close(10)
      end if

      if (include1pi) then
         open(10,File='neutrino_absorption_cross_section_1pi.dat',&
              & position='append')
         if (numberofcalls.ne.1) backspace(10)
         write(10,'(10g13.5)') raiseVal, &
              & sigabsArrFinal(32:33)/real(numberofcalls)
         close(10)
      end if

      if (includeDIS) then
         open(10,File='neutrino_absorption_cross_section_DIS.dat',&
              & position='append')
         if (numberofcalls.ne.1) backspace(10)
         write(10,'(10g13.5)') raiseVal, sigabsArrFinal(34)/real(numberofcalls)
         close(10)
      end if

      open(10,File='neutrino_absorption_cross_section_numbers.dat',&
           & position='append')
      if (numberofcalls.ne.1) backspace(10)
      write(10,'(10g13.5)') raiseVal, numberofsuccessfinal(1)/numberofcalls,&
           & numberofsuccessfinal(2)/numberofcalls, &
           & numberofsuccessfinal(3)/numberofcalls,&
           & numberofsuccessfinal(4)/numberofcalls, &
           & numberofsuccessfinal(7)/numberofcalls,&
           & numberofsuccessfinal(32)/numberofcalls, &
           & numberofsuccessfinal(33)/numberofcalls,&
           & sum(numberofsuccessfinal)/numberofcalls,numtry
      close(10)

    end subroutine DoWrite


   end subroutine init_neutrino

  !****************************************************************************
  !****************************************************************************
  integer function getFirstEvent()

    first=first+1
    getFirstEvent=first
  end function getFirstEvent

  !****************************************************************************
  !****************************************************************************
  function getFirstEventRange()

    integer, dimension (1:2) :: getFirstEventRange
    getFirstEventRange=(/1,first/)
  end function getFirstEventRange

  !****************************************************************************
  !****************************************************************************
  subroutine getNeutrinoInfo(raiseVal_)

    real, intent(out) :: raiseVal_
    raiseVal_=raiseVal
  end subroutine getNeutrinoInfo

  !****************************************************************************
  !****************************************************************************
  logical function neutrinoInit_getRealRun()

    if (readinputflag) then
       call readInput
    end if
    neutrinoInit_getRealRun=realRun
  end function neutrinoInit_getRealRun

  !****************************************************************************
  !****s* initNeutrino/get_init_namelist
  ! NAME
  ! subroutine get_init_namelist(process_ID, flavor_ID,
  ! nuXsectionMode, nuExp, debugflag, includeQE, includeDELTA, includeRES,
  ! include1pi, realRun)
  !
  ! PURPOSE
  ! This subroutine returns any entry of the neutrino init namelist.
  !
  ! OUTPUT
  ! * logical, optional :: debugflag,includeQE,includeDELTA,
  !   includeRES,include1pi,realRun
  ! * integer, optional :: process_ID,flavor_ID,nuXsectionMode,nuExp
  !
  !****************************************************************************
  subroutine get_init_namelist(Gprocess_ID,Gflavor_ID, &
       & GnuXsectionMode, &
       & GnuExp,Gdebugflag,GincludeQE,GincludeDELTA,GincludeRES,Ginclude1pi,&
       & GrealRun, outLepton_ID, outLepton_charge)

    logical, optional, intent(out) :: Gdebugflag, &
         & GincludeQE,GincludeDELTA,GincludeRES,Ginclude1pi,GrealRun
    integer, optional, intent(out) :: Gprocess_ID,Gflavor_ID, &
         & GnuXsectionMode, GnuExp, outLepton_ID, outLepton_charge

    integer::     leptonout_ID, leptonout_charge

    if (present(Gflavor_ID)) Gflavor_ID=flavor_ID
    if (present(Gprocess_ID)) Gprocess_ID=process_ID
    if (present(GnuXsectionMode)) GnuXsectionMode=nuXsectionMode
    if (present(GnuExp)) GnuExp=nuExp
    if (present(Gdebugflag)) Gdebugflag=debugflag
    if (present(GincludeQE)) GincludeQE=includeQE
    if (present(GincludeDELTA)) GincludeDELTA=includeDELTA
    if (present(GincludeRES)) GincludeRES=includeRES
    if (present(Ginclude1pi)) Ginclude1pi=include1pi
    if (present(GrealRun)) GrealRun=realRun


     if (abs(process_ID).eq.2 .or. abs(process_ID).eq.1) then
          leptonout_ID= 900+flavor_ID
          leptonout_charge=-sign(1, process_ID)
     else if (abs(process_ID).eq.3) then
          leptonout_ID= sign(910+flavor_ID, process_ID)
          leptonout_charge=0
     end if

    if (present(outLepton_ID)) outLepton_ID=leptonout_ID
    if (present(outLepton_charge)) outLepton_charge=leptonout_charge


  end subroutine get_init_namelist

  !****************************************************************************
  !****s* initNeutrino/get_runtime_vars
  ! NAME
  ! subroutine get_runtime_vars(sigabsArrFinal,sigabsArr)
  !
  ! PURPOSE
  ! This subroutine returns variables that are changed with every init.
  !
  ! OUTPUT
  ! * real,dimension(0:max_finalstate_ID),optional :: sigabsArrFinal, sigabsArr
  !
  !****************************************************************************
  subroutine get_runtime_vars(GsigabsArrFinal,GsigabsArr)

    real,dimension(0:max_finalstate_ID),optional,intent(out) :: GsigabsArrFinal,GsigabsArr
    if (present(GsigabsArrFinal)) GsigabsArrFinal=sigabsArrFinal
    if (present(GsigabsArr)) GsigabsArr=sigabsArr
  end subroutine get_runtime_vars

end module initNeutrino
