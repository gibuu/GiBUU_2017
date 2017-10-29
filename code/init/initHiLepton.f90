!******************************************************************************
!****m* /initHiLepton
! NAME
! module initHiLepton
!
! PURPOSE
! This module is for high-energetic lepton-nucleus collisions.
! The kinematical setup is such that the exchanged virtual photon runs in
! z direction, while the leptons (incoming/outgoing) move in the x-z plane.
!
! INPUTS
! The namelists "HiLeptonNucleus" and "HiPhotonKinematics".
!******************************************************************************
module initHiLepton

  use eN_eventDefinition, only: electronNucleon_event

  implicit none
  private

  public :: InitHiLeptonInduced
  public :: AccWeight
  public :: GetiExperiment
  public :: GetPhotonKin
  public :: GetEnergies
  public :: HiLepton_getRealRun

  !****************************************************************************
  !****g* initHiLepton/iExperiment
  ! SOURCE
  !
  integer,save :: iExperiment=0
  !
  ! PURPOSE
  ! choice of experiment, detector and energy
  !
  ! possible values are:
  ! *  0: no experiment/fixed kinematics
  ! *  1: Hermes, 27GeV, D,N,Kr
  ! *  2: Hermes, 27GeV, Ne
  ! *  3: Hermes, 27GeV, H
  ! *  4: JLAB, 12GeV
  ! *  5: JLAB,  5GeV
  ! *  6: EMC, 100GeV
  ! *  7: EMC, 120GeV
  ! *  8: EMC, 200GeV
  ! *  9: EMC, 280GeV
  ! * 10: Hermes, 12GeV
  ! * 11: Hermes, 27GeV, arXiv:0704.3270
  ! * 12: Mainz, Yoon: Ebeam=1.5GeV
  ! * 13: Hermes, 27GeV, arXiv:0704.3712 (pT-broadening)
  ! * 14: JLAB,  5GeV, rho0 experiment
  ! * 15: JLAB,  4GeV, rho0 experiment
  ! * 16: EIC, E_e and E_A given explicit (3+30,11+30,4+100)
  ! * 17: no detector, total cross section, Ebeam
  ! * 18: E665, 470GeV
  ! * 19: CLAS/JLAB, 12GeV RunGroupA optimized 10.6 GeV
  ! * 20: CLAS/JLAB, 12GeV RunGroupA theoterical
  !
  ! please note:
  ! The entry "iExperiment == 0" replaces the old HiPhoton event type.
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/iDetector
  ! SOURCE
  !
  integer,save :: iDetector = -1
  !
  ! PURPOSE
  ! This sets the treatment of the detector:
  ! * -1 : not valid/not initialized/use default
  ! *  0 : no detector, as AccFlag=.false.
  ! *  1 : HERMES, full efficiency
  ! *  2 : EMC, full efficiency
  ! *  3 : CLAS, only cuts (th_e=12°..50°, th_hadron=6°..143°)
  ! *  4 : CLAS, full efficiency + cuts as for 5GeV
  ! *  5 : CLAS, electron: cuts (th_e=12°..50°),
  !              hadrons: efficiency+cuts as for 5GeV
  ! * 90 : full acceptance
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/shadow
  ! SOURCE
  !
  logical,save :: shadow=.true.
  !
  ! PURPOSE
  ! flag: Consider shadowing or not
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/minimumMomentum
  ! SOURCE
  !
  real,save    :: minimumMomentum=0.1
  !
  ! PURPOSE
  ! minimal momentum considered. (in GeV)
  !****************************************************************************


  !****************************************************************************
  !****g* initHiLepton/ModusCalcFluxNorm
  ! SOURCE
  !
  logical,save :: ModusCalcFluxNorm=.false.
  !
  ! PURPOSE
  ! if this flag is true, than we do not really generate events. We only
  ! select nu and Q2 according an equal distribution and plot the flux (and
  ! the flux multiplied with AccWeight).
  ! Normally we choose nu and Q2 according flux*Accweight via von-Neumann-
  ! rejection method (where we loose access to the absolute normalisation).
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/realRun
  ! SOURCE
  !
  logical,save :: realRun=.false.
  !
  ! PURPOSE
  ! Flag to indicate, whether we produce real or perturbative particles.
  ! NOTES
  ! run with real particles untested !!!
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/DoStatistics
  ! SOURCE
  !
  logical,save :: DoStatistics=.false.
  !
  ! PURPOSE
  ! switch on/off statistical output of init routines
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/Ebeam
  ! SOURCE
  !
  real, save :: Ebeam
  !
  ! PURPOSE
  ! electron beam energy [GeV]
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/user_ymax
  ! SOURCE
  real :: user_ymax = -99.9
  ! PURPOSE
  ! user given value for ymax, overrides default value if reasonable
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/user_numin
  ! SOURCE
  real :: user_numin = -99.9 ! GeV
    ! PURPOSE
  ! user given value for numin, overrides default value if reasonable
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/user_numax
  ! SOURCE
  real :: user_numax = -99.9 ! GeV
  ! PURPOSE
  ! user given value for numax, overrides default value if reasonable
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/user_costmin
  ! SOURCE
  real :: user_costmin = -99.9
  ! PURPOSE
  ! user given value for costmin, overrides default value if reasonable
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/user_costmax
  ! SOURCE
  real :: user_costmax =  99.9
  ! PURPOSE
  ! user given value for costmax, overrides default value if reasonable
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/user_smin
  ! SOURCE
  real :: user_smin = -99.9
  ! PURPOSE
  ! user given value for smin, overrides default value if reasonable
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/user_xBmin
  ! SOURCE
  real :: user_xBmin = -99.9
  ! PURPOSE
  ! user given value for xBmin, overrides default value if reasonable
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/user_qsqmin
  ! SOURCE
  real :: user_qsqmin = -99.9
  ! PURPOSE
  ! user given value for qsqmin, overrides default value if reasonable
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/user_maxw
  ! SOURCE
  real :: user_maxw = -99.9
  ! PURPOSE
  ! user given value for maxw, overrides default value if reasonable
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/earlyPauli
  ! SOURCE
  !
  logical,save :: earlyPauli=.false.
  !
  ! PURPOSE
  ! Flag to indicate, whether we should check Pauli blocking already during
  ! generation or only at the end.
  !
  ! if .false. (default), events will be generated in a first stage without
  ! Pauli blocking. This is then tested afterwards. If the generated event
  ! is blocked, it will be redone! Thus Pauli blocking does *not* change the
  ! total cross section, only the relative strength will be reshuffled.
  !
  ! if .true., then blocked events will be excluded from the Monte Carlo
  ! decision and the total cross section will be reduced.
  !
  ! NOTES
  ! The behaviour, if no event at all is possible, is at the moment a little
  ! bit unpredictable ;)
  !****************************************************************************

  logical,save :: initFlag=.true.


  !****************************************************************************
  ! some global (module) variables:

  real, save :: numin,numax
  real, save :: costmin,costmax
  real, save :: ymax
  real, save :: Q2min,Q2max
  real, save :: smin
  real, save :: xBmin
  real, save :: maxw
  logical,parameter :: AccFlag=.true.
  real, save :: weight1,weight2

  Type(electronNucleon_event), save :: eNev_InitData

  real, parameter, dimension(0:5,0:20) :: maxwArr = RESHAPE ( (/ &
       & 0. , 0.     , 0.     , 0.     , 0.     , 0.     ,& !  0: 0
       & 0. , 4.57e-1, 0.     , 0.     , 0.     , 0.     ,& !  1: 1
       & 0. , 0.     , 0.     , 0.     , 0.     , 0.     ,& !  2: 1
       & 0. , 4.57e-1, 0.     , 0.     , 0.     , 0.     ,& !  3: 1
       & 0. , 0.     , 0.     , 6.20e-3, 0.     , 0.     ,& !  4: 3
       & 0. , 0.     , 0.     , 0.     , 6.60e-4, 0.     ,& !  5: 4
       & 0. , 0.     , 1.15e-2, 0.     , 0.     , 0.     ,& !  6: 2
       & 0. , 0.     , 1.30e-2, 0.     , 0.     , 0.     ,& !  7: 2
       & 0. , 0.     , 1.10e-2, 0.     , 0.     , 0.     ,& !  8: 2
       & 0. , 0.     , 8.50e-3, 0.     , 0.     , 0.     ,& !  9: 2
       & 0. , 4.10e-3, 0.     , 0.     , 0.     , 0.     ,& ! 10: 1
       & 0. , 4.57e-1, 0.     , 0.     , 0.     , 0.     ,& ! 11: 1
       & 0. , 0.     , 0.     , 0.     , 0.     , 0.     ,& ! 12: 0
       & 0. , 4.57e-1, 0.     , 0.     , 0.     , 0.     ,& ! 13: 1
       & 0. , 0.     , 0.     , 1.93e-3, 1.19e-3, 1.80e-3,& ! 14: 4
       & 0. , 0.     , 0.     , 0.     , 8.15e-4, 0.     ,& ! 15: 4
       & 0. , 0.     , 0.     , 0.     , 0.     , 0.     ,& ! 16: 0
       & 1.e-4 , 0.  , 0.     , 0.     , 0.     , 0.     ,& ! 17: 0
       & 5.5e-2, 0.  , 0.     , 0.     , 0.     , 0.     ,& ! 18: 0
       & 0. , 0.     , 0.     , 6.20e-3, 0.     , 0.     ,& ! 19: 3
       & 0. , 0.     , 0.     , 6.20e-3, 0.     , 0.     &  ! 20: 3
       &/) , (/6,21/) )

  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/EIC_Ee
  ! SOURCE
  !
  real, save :: EIC_Ee = -99.9
  !
  ! PURPOSE
  ! the electron beam energy, if iExperiment=EIC
  !****************************************************************************

  !****************************************************************************
  !****g* initHiLepton/EIC_EA
  ! SOURCE
  !
  real, save :: EIC_EA = -99.9
  !
  ! PURPOSE
  ! the hadron beam energy, if iExperiment=EIC
  !****************************************************************************

contains

  !****************************************************************************
  !****f* initHiLepton/HiLepton_getRealRun
  ! NAME
  ! logical function HiLepton_getRealRun()
  ! PURPOSE
  ! return the value of realRun
  !****************************************************************************
  logical function HiLepton_getRealRun()
    HiLepton_getRealRun = realRun
  end function HiLepton_getRealRun

  !****************************************************************************
  !****s* initHiLepton/InitHiLeptonInduced
  ! NAME
  ! subroutine InitHiLeptonInduced(rParts,pParts,targetNuc)
  !
  ! PURPOSE
  ! This routine initializes the total given perturbative particle vector
  ! "pParts" by calling "genHiPhotonEvent" for every nucleon given in "rParts"
  ! and inserting the output of every successful event into "pParts".
  !
  ! For every incoming "photon" a new set of the kinematic variables
  ! according the constraints of the given experiment/detector is choosen.
  !
  ! Particle weights are set by the total XS of the generating event.
  !
  ! Outgoing particles are in a system, where the (intermediate) PHOTON
  ! defines the z-axis, not the (incoming) LEPTON.
  ! The scattering happens in a plane (except fermi motion), where the
  ! second momentum component (i.e. y-direction) vanishes.
  !
  ! INPUTS
  ! * type(particle),dimension(:,:) :: rParts
  ! * type(tNucleus),pointer        :: targetNuc
  !
  ! OUTPUT
  ! * type(particle),dimension(:,:) :: pParts
  !
  !****************************************************************************
  subroutine InitHiLeptonInduced(rParts,pParts,targetNuc)

    use particleDefinition
    use nucleusDefinition
    use baryonPotentialModule, only: getNoPertPot_baryon
    use checks, only: ChecksCallEnergy
    use CollHistory, only: CollHist_ClearArray,CollHist_UpdateHist
    use constants, only: mN
    use eN_event, only: eNeV_GetKinV,eNev_init_Target
    use energyCalc, only: updateEnergies
    use eventGenerator_eN_HiEnergy, only: eventGen_eN_HiEnergy
    use EventInfo_HiLep, only: EventInfo_HiLep_Init,EventInfo_HiLep_Store, &
         EventInfo_HiLep_Dump
    use hadronFormation, only: forceInitFormation
    use histf90
    use hist2Df90
    use insertion, only: FindLastUsed, setIntoVector
    use ParamEP, only: CalcParamEP_R1990
    use pauliBlockingModule, only: checkPauli
    use PythiaSpecFunc, only: Init_VM_Mass
    use random, only: rn
    use shadowing, only: AeffCalc
    use sourceAnalysis, only: getSMM_Flag
    use output, only: DoPr, WriteParticle, intTochar

    type(particle), dimension(:,:),intent(inout),TARGET :: rParts
    type(particle), dimension(:,:),intent(inout),TARGET :: pParts
    type(tNucleus), pointer :: targetNuc

    ! Input variables according jobcard:

    integer :: i
    integer :: iEns,nEns, iPart,nPart

    real, dimension(4) :: scaleVMD
    real :: Q2, nu, eps, W, Wfree
    real :: phiLep, R

    integer, parameter :: nOutPart = 100
    type(particle),dimension(nOutPart) :: outPart
    type(particle)                     :: target

    real    :: XS_tot
    real, dimension(0:4) :: XS_Arr
    real, dimension(8)   :: XS_Arr_low
    logical :: flagOK
    logical :: NumbersAlreadySet
    integer :: nTries,nTries2
    integer, save :: nCalls = 0

    real,save                :: Sum_XS
    real,dimension(0:4),save :: Sum_XS_Arr
    real,dimension(8),save   :: Sum_XS_Arr_low

    type(histogram2D), save :: h2D_shad, h2D_Pos
    type(histogram2D), save :: h2D_nuQ2_Eps,h2D_nuQ2_EpsR
    type(histogram2D), save :: h2D_nuQ2_Flux,h2D_nuQ2_FluxW
    type(histogram2D), dimension(0:5), save :: h2D_nuQ2
!!$    type(histogram), save :: h_RhoTP, h_RhoTF ! FOR TEST PURPOSES

    type(histogram2D), save :: h2D_WvsWfree,h2D_WvsW0
    type(histogram), save :: h_pFermi

    real, dimension(2) :: P1,P2,dP

    integer, dimension(2,nOutPart) :: posOut
    type(particle),save :: partPhoton
    type(particle), POINTER :: pPart
    integer :: channel

    type(electronNucleon_event), save :: eNev_RunData
    real :: ptot(0:3)
    real,parameter :: add0=1e-20
    real :: mul0,mul1
    integer :: iiPart=0, iNr

    type(particle),dimension(:,:),pointer :: pArr
    real :: bT,bZ

    if (DoPr(2)) write(*,*) '**Initializing high-energy lepton induced reactions'

    if (initFlag) then
       call initInput ! Read input
       call forceInitFormation

       call CreateHist2D(h2D_shad,"ShadFak",(/-10.,0./),(/10.,10./),(/.1,.1/),.true.)
       call CreateHist2D(h2D_Pos,"Interaction pos",(/-10.,0./),(/10.,10./),(/.1,.1/),.true.)

       P1 = (/numin,Q2min/)
       P2 = (/numax,Q2max/)
       dP = (P2-P1)/100
       if (iExperiment.eq.0) dP = (P2-P1)/1

       call CreateHist2D(h2D_nuQ2_Flux, "nuQ2: Flux",P1,P2,dP,.true.)
       call CreateHist2D(h2D_nuQ2_FluxW,"nuQ2: Flux (weights)",P1,P2,dP,.true.)

       call CreateHist2D(h2D_nuQ2(0),"nuQ2: considered",P1,P2,dP,.true.)
!!$       call CreateHist2D(h2D_nuQ2(1),"nuQ2: 103 [phi;el]",P1,P2,dP,.true.)
!!$       call CreateHist2D(h2D_nuQ2(2),"nuQ2: 300 [Lambda K]",P1,P2,dP,.true.)
!!$       call CreateHist2D(h2D_nuQ2(3),"nuQ2: 400 [Sigma  K]",P1,P2,dP,.true.)
!!$       call CreateHist2D(h2D_nuQ2(4),"nuQ2: 999 [misc]",P1,P2,dP,.true.)

!!$       call CreateHist(h_RhoTP, "t_P (Rho0)", 0.0,20.0,0.1) ! FOR TEST PURPOSES
!!$       call CreateHist(h_RhoTF, "t_F (Rho0)", 0.0,20.0,0.1)

       call CreateHist2D(h2D_WvsWfree, "W vs. Wfree", (/0.,0./),(/4.,4./),(/0.02,0.02/) )
       call CreateHist2D(h2D_WvsW0,    "W vs. W0", (/0.,0./),(/4.,4./),(/0.02,0.02/) )
       call CreateHist(h_pFermi, "p_Fermi", 0., 1.0, 0.01)

       call CreateHist2D(h2D_nuQ2_Eps, "nuQ2: eps",P1,P2,dP,.true.)
       call CreateHist2D(h2D_nuQ2_EpsR, "nuQ2: eps*R",P1,P2,dP,.true.)

       Sum_XS = 0.0
       Sum_XS_Arr = 0.0
       Sum_XS_Arr_low = 0.0


       call setToDefault(partPhoton)
       partPhoton%ID = 999

       initFlag=.false.
    end if

    if (DoPr(2)) write(*,*) '**'
    nCalls = nCalls+1

    nEns=size(rParts,dim=1)

    nPart = FindLastUsed(rParts(1,:))

    mul0 = 1./nCalls
    mul1 = 1./(nCalls*nEns*nPart)

    call EventInfo_HiLep_Init(nEns,nPart)

    call CollHist_ClearArray()

    if ((realRun).or.(.not.getNoPertPot_baryon())) then
       call updateEnergies(rParts)
    end if

    call ChecksCallEnergy(-99.9,rParts)

    weight1 = 1.0 ! Only for h2D_nuQ2_FluxW !!!
    weight2 = 1.0 ! Only for h2D_nuQ2_FluxW !!!

    if (realRun) then
       pArr => rParts
    else
       pArr => pParts
    end if

    ! setting last collision time of all target nucleons to the past:

    rParts%lastCollisionTime = -9.9 ! DUMMY value

    do iEns=1,nEns
       if (realRun) iiPart = int(rn()*nPart)+1 ! int just cuts the digits!

!       write(*,*) 'iiPart:',iEns,iiPart

       ParticleLoop : do iPart=1,nPart
          if (realRun) then
             if (iPart.ne.iiPart) cycle
          end if

          target = rParts(iEns,iPart)
          nTries2 = 0

101       continue

          !... determine photon kinematics
          nTries = 0
          do
             nTries = nTries+1

             if (iExperiment.eq.0) then
                eNev_RunData = eNev_InitData
             else
                call photonFluxMC(eNev_RunData,target)
             end if

!             call write_electronNucleon_event(eNev_RunData)

             call eNev_init_Target(eNev_RunData,target,flagOK)
             if (flagOK) exit

             if (nTries.gt.1000) then
                write(*,*) 'FluxLoop: more than 1000 tries. Skip target!'
                cycle ParticleLoop
             end if
          end do

!!$          call WriteParticle(6,iEns,iPart,target)
!!$          call write_electronNucleon_event(eNev_RunData, .false.)

          call eNeV_GetKinV(eNev_RunData, nu,Q2,W,Wfree,eps)
          phiLep = atan2(eNev_RunData%lepton_in%momentum(2),eNev_RunData%lepton_in%momentum(1))

          call CalcParamEP_R1990(W,Q2,R)
          call AddHist2D(h2D_nuQ2_Eps, (/nu,Q2/), 1.0,eps)
          call AddHist2D(h2D_nuQ2_EpsR, (/nu,Q2/), 1.0,eps*R)

          call AddHist2D(h2D_nuQ2_Flux, (/nu,Q2/), 1.0)
          call AddHist2D(h2D_nuQ2_FluxW, (/nu,Q2/), weight2, weight1)

          ptot = eNev_RunData%boson%momentum +(/mN,0.,0.,0./)

          call AddHist2D(h2D_WvsWfree, (/W,Wfree/), 1.0)
          call AddHist2D(h2D_WvsW0,    (/W,sqrt(ptot(0)**2-sum(ptot(1:3)**2)) /), 1.0)
          call AddHist(h_pFermi, absMom(target), 1.0)

          if (ModusCalcFluxNorm) cycle ParticleLoop ! do not really generate events

          ! save W & position for calls of VM_Mass (Otherwise VM_Mass doesn't know about srts!!)
          call Init_VM_Mass(Wfree,target%position)

!          write(*,'(i3,A,1P,10g12.5)') ipart,'  nu,Q2,Ebeam,eps,W = ',nu,Q2,Ebeam,eps,W

          !...shadowing
          if (shadow) then
             call AeffCalc(targetNuc,target%position,nu,Q2,scaleVMD)
          else
             scaleVMD=1.
          end if


          !( save some of the parameters already here...)
          iNr = iPart
          if (realRun) iNr = -1
          call EventInfo_HiLep_Store(iEns,iNr,0.0,nu,Q2,eps,0,Ebeam)

          nTries = 0
          RetryLoop: do
             nTries = nTries+1

             if (nTries.gt.1000) then
                write(*,*) 'RetryLoop: more than 1000 tries. Skip target!'
                cycle ParticleLoop
             end if

             call ResetNumberGuess()

             call eventGen_eN_HiEnergy(eNev_RunData,iEns*1000+iPart,scaleVMD, &
                  earlyPauli,rParts, &
                  outPart,channel,flagOK,XS_tot,XS_Arr,XS_Arr_low )

             if (abs(XS_tot).lt.1e-10) then
!                write(*,*) 'RetryLoop: XS_tot=0.'

                nTries2 = nTries2+1
                if (nTries2.lt.1000) goto 101
                write(*,*) 'RetryLoop: 1000 times XS_tot=0.'
                cycle ParticleLoop
!!$                write(*,*) 'RetryLoop: XS_tot=0. Skip target!'
!!$                cycle ParticleLoop
             end if

             if (.not.flagOK) cycle RetryLoop

             outPart%perWeight = XS_tot

             !...test charge, baryon, energy and momentum conservation
             flagOK = checkConservation(eNev_RunData,outPart,channel)
             if (.not.flagOK) then
                if (DoPr(1)) write(*,*) 'Redo event: conservation! ',iEns,iPart
                cycle RetryLoop
             end if

             !...check Pauli Blocking
             flagOK = checkPauli(outPart, rParts)
             if (.not.flagOK) then
                if (DoPr(1)) write(*,*) 'Redo event: pauli blocked! ',iEns,iPart
                cycle RetryLoop
             end if

             exit RetryLoop
          end do RetryLoop

          !...Give the particles their (final) number

          NumbersAlreadySet = AcceptGuessedNumbers()

          posOut = 0
          if (realRun) then
             outPart%perturbative = .false.
             call setIntoVector(outPart(:),pArr(iEns:iEns,:), &
                  & flagOK,NumbersAlreadySet,positions=posOut(:,:))
             rParts(iEns,iPart)%ID = 0
             rParts(iEns:iEns,:)%perweight = XS_tot
             rParts(iEns:iEns,:)%firstevent = iEns*1000+iPart
          else
             call setIntoVector(outPart(:),pArr(iEns:iEns,:), &
                  & flagOK,NumbersAlreadySet,positions=posOut(:,:))
          end if
          posOut(1,:) = iEns

          call CollHist_UpdateHist((/partPhoton,target/), outPart, &
                 & (/0,0/), posOut, XS_tot)

          iNr = iPart
          if (realRun) iNr = -1
          call EventInfo_HiLep_Store(iEns,iNr,XS_tot,nu,Q2,eps,channel,Ebeam,phiLep,W,target%position)

          !...Delete particles with small momenta !!!

          do i=1,nOutPart
             if (posOut(2,i).lt.1) cycle
             pPart => pArr(posOut(1,i),posOut(2,i))
             if (pPart%ID <= 0) cycle
             if (DOT_PRODUCT(pPart%momentum(1:3),pPart%momentum(1:3)) .lt. minimumMomentum**2) then
                if (DoPr(1)) write(*,*) '~~~ Deleting particle (momentum too small).'
                if (DoPr(1)) call WriteParticle(6,posOut(1,i),posOut(2,i),pPart)
                pPart%ID = 0
             end if
          end do

!!$          !... FOR TEST PURPOSES:
!!$
!!$          do i=1,nOutPart
!!$             if (posOut(2,i).lt.1) cycle
!!$             pPart => pArr(posOut(1,i),posOut(2,i))
!!$
!!$             if (pPart%ID <= 0) cycle
!!$
!!$             if (pPart%ID .ne. 103) cycle ! only rho
!!$             if (pPart%charge .ne. 0) cycle ! only rho0
!!$             if (pPart%momentum(0) .lt. 0.9*nu) cycle ! only rho0 with z>0.9
!!$
!!$             call AddHist(h_RhoTP, pPart%productionTime,XS_tot)
!!$             call AddHist(h_RhoTF, pPart%formationTime,XS_tot)
!!$
!!$          end do

          !... Some statistics per target nucleon:

          bZ = target%position(3)
          bT = sqrt(target%position(1)**2+target%position(2)**2)
          call AddHist2D(h2D_shad,(/bZ,bT/),XS_Arr(1),scaleVMD(1)*XS_Arr(1))
          call AddHist2D(h2D_Pos, (/bZ,bT/),XS_tot)

          call AddHist2D(h2D_nuQ2(0),(/nu,Q2/), XS_tot)
!!$          select case(HiPhotonEventType)
!!$          case (103)
!!$             call AddHist2D(h2D_nuQ2(1),(/nu,Q2/), XS_tot)
!!$          case (300)
!!$             call AddHist2D(h2D_nuQ2(2),(/nu,Q2/), XS_tot)
!!$          case (400)
!!$             call AddHist2D(h2D_nuQ2(3),(/nu,Q2/), XS_tot)
!!$          case DEFAULT
!!$             call AddHist2D(h2D_nuQ2(4),(/nu,Q2/), XS_tot)
!!$          end select


          Sum_XS = Sum_XS + XS_tot
          Sum_XS_Arr = Sum_XS_Arr + XS_Arr
          Sum_XS_Arr_low = Sum_XS_Arr_low + XS_Arr_low

!!$          write(*,*) 'channel = ',channel
!!$          write(*,*) XS_Arr_low




       end do ParticleLoop
    end do

    if (DoPr(2)) write(*,*) '^^'

    W = sqrt(max(0., mN**2+2*mN*nu-Q2))
!    write(165,'(20ES14.5)') nu,Q2,W,eps, Sum_XS*mul1, Sum_XS_Arr_low*mul1, Sum_XS_Arr*mul1

    if (DoStatistics) then

       call WriteHist2D_Gnuplot(h2D_shad,DoAve=.true.,MaxVal=-1.0,file='initHiLep.Shadowing.dat',dump=.true.)

       call WriteHist2D_Gnuplot(h2D_Pos,add=add0,mul=mul0,file='initHiLep.Pos.dat',dump=.true.)

       call WriteHist2D_Gnuplot(h2D_nuQ2_Eps,DoAve=.true.,MaxVal=-1.0,file='initHiLep.nuQ2.Eps.dat',dump=.true.)
       call WriteHist2D_Gnuplot(h2D_nuQ2_EpsR,DoAve=.true.,MaxVal=-1.0,file='initHiLep.nuQ2.EpsR.dat',dump=.true.)

       call WriteHist2D_Gnuplot(h2D_nuQ2_flux,add=add0,mul=mul0,file='initHiLep.nuQ2.Flux.dat',dump=.true.)
       call WriteHist2D_Gnuplot(h2D_nuQ2_fluxW,add=add0,mul=mul0,file='initHiLep.nuQ2.FluxW.dat',dump=.true.)
       do i=0,4
          call WriteHist2D_Gnuplot(h2D_nuQ2(i),add=add0,mul=mul0,file='initHiLep.nuQ2.'//intTochar(i)//'.dat',dump=.true.)
       end do

!!$       call WriteHist(h_RhoTP,add=add0,mul=mul0,file='initHiLep.RhoTP.dat',dump=.true.) ! FOR TEST PURPOSES
!!$       call WriteHist(h_RhoTF,add=add0,mul=mul0,file='initHiLep.RhoTF.dat',dump=.true.)

       call WriteHist2D_Gnuplot(h2D_WvsWfree,add=add0,mul=mul1,file='initHiLep.WvsWfree.dat',dump=.true.)
       call WriteHist2D_Gnuplot(h2D_WvsW0,add=add0,mul=mul1,file='initHiLep.WvsW0.dat',dump=.true.)
       call WriteHist(h_pFermi,add=add0,mul=mul1,file='initHiLep.pFermi.dat',dump=.true.)


       if (getSMM_Flag()) then
          call EventInfo_HiLep_Dump()
       end if

    end if ! DoStatistics

  end subroutine InitHiLeptonInduced

  !****************************************************************************
  !****s* initHiLepton/initInput
  ! NAME
  ! subroutine initInput
  !
  ! PURPOSE
  ! read initialisation of "InitHiLeptonInduced" from input file,
  ! namelist "HiLeptonNucleus"
  !
  !****************************************************************************
  subroutine initInput

    use checks, only: ChecksSwitchRealRun
    use eventGenerator_eN_HiEnergy, only: ReadHiGammaNucleus
    use Dilepton_Analysis, only: Dilep_Init
    use output, only: writeFileDocu, Write_ReadingInput
    use callstack, only: traceback

    !**************************************************************************
    !****n* initHiLepton/HiLeptonNucleus
    ! NAME
    ! NAMELIST /HiLeptonNucleus/
    ! PURPOSE
    ! Namelist for initHiLepton includes:
    ! * iExperiment
    ! * shadow
    ! * minimumMomentum
    ! * ModusCalcFluxNorm
    ! * iDetector
    ! * EIC_Ee
    ! * EIC_EA
    ! * realRun
    ! * DoStatistics
    ! * user_numin
    ! * user_numax,
    ! * user_costmin
    ! * user_costmax
    ! * user_ymax
    ! * user_smin
    ! * user_xBmin
    ! * user_qsqmin
    ! * user_maxw
    ! * earlyPauli
    !**************************************************************************
    NAMELIST /HiLeptonNucleus/ iExperiment,shadow, &
         minimumMomentum, ModusCalcFluxNorm, iDetector, &
         EIC_Ee, EIC_EA, Ebeam, realRun, &
         DoStatistics, &
         user_numin, user_numax, &
         user_costmin, user_costmax, &
         user_ymax, &
         user_smin, &
         user_xBmin, &
         user_qsqmin, &
         user_maxw, &
         earlyPauli

    character(40), dimension(0:20) :: NN
    character(40), dimension(0: 5) :: ND
    integer :: ios

    NN( 0)= 'no experiment/fixed kinematics'
    NN( 1)= 'Hermes:27, D,N,Kr'
    NN( 2)= 'Hermes:27, Ne'
    NN( 3)= 'Hermes:27, H'
    NN( 4)= 'JLAB:12'
    NN( 5)= 'JLAB:5'
    NN( 6)= 'EMC:100'
    NN( 7)= 'EMC:120'
    NN( 8)= 'EMC:200'
    NN( 9)= 'EMC:280'
    NN(10)= 'Hermes:12'
    NN(11)= 'Hermes:27, 0704.3270'
    NN(12)= 'Mainz, Yoon: Ee=1.5GeV'
    NN(13)= 'Hermes:27, 0704.3712'
    NN(14)= 'JLAB:5, rho0'
    NN(15)= 'JLAB:4, rho0'
    NN(16)= 'EIC'
    NN(17)= 'total cross section'
    NN(18)= 'E665:470'
    NN(19)= 'CLAS/JLAB12 RunA'
    NN(20)= 'CLAS/JLAB12 RunA'
    ND( 0)= 'no detector'
    ND( 1)= 'Hermes, full efficiency'
    ND( 2)= 'EMC, full efficiency'
    ND( 3)= 'CLAS, only cuts'
    ND( 4)= 'CLAS, full efficiency'
    ND( 5)= 'CLAS, electron: cuts; hadrons: eff.'

    call ReadHiGammaNucleus

    iExperiment = -1
    iDetector   = -1
    Ebeam = -99.9

    call Write_ReadingInput('HiLeptonNucleus',0)
    rewind(5)
    read(5,nml=HiLeptonNucleus,iostat=ios)

    if (ios>0) then
       write(*,*)
       write(*,*) '>>>>> parameters moved to "HiGammaNucleus" !?!'
    end if
    call Write_ReadingInput('HiLeptonNucleus',0,ios)

    if (iExperiment.eq.-1) &
         & call TRACEBACK('ERROR: you must provide iExperiment!!!')

    if (iDetector.eq. -1) then ! value not set; set to default
       select case (iExperiment)
       case (0)
          iDetector = 0 ! full acceptance
       case (1:3,10,11,13)
          iDetector = 1 ! HERMES
       case (6:9)
          iDetector = 2 ! EMC
       case (4,19,20)
          iDetector = 3 ! CLAS
       case (5,14,15)
          iDetector = 4 ! CLAS
       case (12)
          iDetector = 0 ! Mainz
       case (16)
          iDetector = 0 ! EIC
       case (17)
          iDetector = 0 ! total cross section
       case (18)
          iDetector = 0 ! E665
       end select
    end if

    if (.not.AccFlag) iDetector = 0

    write(*,*) '  iExperiment =',iExperiment, ' >',trim(NN(iExperiment)),'<'
    write(*,*) '  iDetector   =',iDetector,   ' >',trim(ND(iDetector)),'<'
    write(*,*) '  useShadowing       =',shadow
    write(*,*) '  minimumMomentum    =',minimumMomentum
    write(*,*) '  realRun            =',realRun
    write(*,*) '  DoStatistics       =',DoStatistics
    write(*,*) '  EarlyPauli         =',EarlyPauli
    write(*,*)

    if (realRun) then
       write(*,*)
       write(*,*) '  ATTENTION: real run !!!'
       write(*,*)
       call ChecksSwitchRealRun(.true.)
    end if

    select case (iExperiment)
    case (16)
       if (EIC_Ee.le.0.0) &
            & call TRACEBACK('ERROR: you must provide EIC_Ee!!!')
       if (EIC_EA.le.0.0) &
            & call TRACEBACK('ERROR: you must provide EIC_EA!!!')
       write(*,*) '  EIC:',EIC_Ee,' + ',EIC_EA
    case (17)
       if (Ebeam .lt. 0.0) &
            & call TRACEBACK('ERROR: you must provide Ebeam!!!')
       write(*,*) '  Ebeam:',Ebeam
    end select

    if (ModusCalcFluxNorm) then
       write(*,*)
       write(*,*) '  ATTENTION: We run the code in the CalcFluxNorm mode!'
       write(*,*)
    end if

    if (iExperiment.eq.0) then
       call initFixedKin
    else
       call initExperiment
    end if

    call Dilep_Init (Ebeam)  ! initialize the dilepton-analysis module

    if (DoStatistics) then
       call writeFileDocu('initHiLep.Shadowing.dat','shadowing faktor for (bZ,bT)')
       call writeFileDocu('initHiLep.Pos.dat','interaction position: XS(bZ,bT)')
       call writeFileDocu('initHiLep.nuQ2.Eps[R].dat','eps(nu,Q2) [weighted with R=sigma_L/sigma_T]')
       call writeFileDocu('initHiLep.nuQ2.Flux[W].dat','')
       call writeFileDocu('initHiLep.nuQ2.n.dat, n=0..4','')
!!$       call writeFileDocu('initHiLep.RhoTP.dat','production time')
!!$       call writeFileDocu('initHiLep.RhoTF.dat','formation time')
       call writeFileDocu('initHiLep.WvsWfree.dat','W vs. Wfree')
       call writeFileDocu('initHiLep.WvsW0.dat','W vs. W0')
       call writeFileDocu('initHiLep.pFermi.dat','fermi momentum')
    end if


    call Write_ReadingInput('HiLeptonNucleus',1)

  end subroutine initInput

  !****************************************************************************
  !****s* initHiLepton/initExperiment
  ! NAME
  ! subroutine initExperiment
  !
  ! PURPOSE
  ! set kinematical constraints according the given experiment
  !****************************************************************************
  subroutine initExperiment

    use constants, only: melec,mN
    use lorentzTrafo, only: lorentzCalcBeta, lorentz
    use callstack, only: traceback

    real :: Q2max1,Q2max2,qsqmin
    real, dimension(0:3) :: pE,pA
    real, dimension(1:3) :: beta

    smin=4.
    ymax=0.85

    select case (iExperiment)

    case (1:3,11,13) !=== HERMES@27GeV

       Ebeam=27.6
       xBmin=0.02
       costmin=cos(0.2202)
       costmax=cos(0.04)
       qsqmin=1.

       select case (iExperiment)
       case (1) ! D,N,Kr
          numin=7.
       case (2) ! Ne
          numin=2.
       case (3) ! H
          numin=5.4
          xBmin=0.
       case (11) ! final paper
          numin=6.
       case (13) ! pT-broadening
          smin = 10.
          numin=5.4
       end select

    case (10) !=== HERMES@12GeV

       Ebeam=12.0
       xBmin=0.
       costmin=cos(0.2202)
       costmax=cos(0.04)
       qsqmin=0.5
       numin=2.

    case (4:5,14:15,19:20) !=== JLab

       costmax=cos(0.1047) !   6°
       costmin=cos(2.4958) ! 142°
       xBmin=0.
       numin=2.
       qsqmin=1.
       select case (iExperiment)
       case (4)!JLab after 12 GeV Upgrade
          Ebeam=12.
       case (5)!JLab currently 5 GeV
          Ebeam=5.
       case (14)!JLab currently 5 GeV, rho0 experiment
          costmax=cos(0.1745) ! 10°
          costmin=cos(1.0472) ! 60°
          Ebeam=5.
          numin=2.0
          qsqmin=0.8
          ymax=0.97
       case (15)!JLab currently 4 GeV, rho0 experiment
          Ebeam=4.
          numin=2.05
          qsqmin=0.6
          ymax=0.97
       case (19)!JLab after 12 GeV Upgrade RunGroupA optimized
          Ebeam=10.6
       case (20)!JLab after 12 GeV Upgrade RunGroupA theoretical
          Ebeam=11.
       end select

    case (6:9) !=== EMC

       costmin=-1.
       xBmin=0.
       select case (iExperiment)
       case (6)
          Ebeam=100.
          costmax=cos(0.016)
          numin=10.
          qsqmin=2.
       case (7)
          Ebeam=120.
          costmax=cos(0.016)
          numin=10.
          qsqmin=2.
       case (8)
          Ebeam=200.
          costmax=cos(0.014)
          numin=30.
          qsqmin=2.
       case (9)
          Ebeam=280.
          costmax=cos(0.014)
          numin=50.
          qsqmin=5.
       end select

    case (12) !=== Mainz, Yoon, Ebeam=1.5GeV

       Ebeam = 1.5
!       Ebeam = 4.5

       smin = 1.4 ! GeV^2

       costmin=-1.0
       costmax= 1.0

       ymax = 1.0-1e-3
       xBmin = 0.0

       numin = 0.0
!       numax=Ebeam*ymax

!       qsqmin = 0.0
       qsqmin = 0.01

    case (16) !=== EIC

       pE=(/EIC_Ee,0.0,0.0,sqrt(EIC_Ee**2-melec**2)/)
       pA=(/EIC_EA,0.0,0.0,-sqrt(EIC_EA**2-mN**2)/)

       write(*,*)
       write(*,*) ' Lab system:'
       write(*,'("    ",A,1P,4e13.5)') 'pE = ',pE
       write(*,'("    ",A,1P,4e13.5)') 'pA = ',pA

       beta = lorentzCalcBeta(pA)
       call lorentz(beta,pE)
       call lorentz(beta,pA)

       write(*,*) ' fixed target system:'
       write(*,'("    ",A,1P,4e13.5)') 'pE = ',pE
       write(*,'("    ",A,1P,4e13.5)') 'pA = ',pA
       write(*,*)

       Ebeam = pE(0)

       write(*,*) '  ==> shifting to fixed target with Ebeam=',Ebeam
       write(*,*)

       costmin=-1.0
       costmax= 1.0

       xBmin = 0.0

       numin = 5.0
       ymax = 1.0-1e-3
!       numax=Ebeam*ymax

       numin = 0.1*Ebeam ! == ymin = 0.1
       ymax = 0.8

!       qsqmin = 0.0
       qsqmin = 0.1

    case (17) !=== total cross section

       ! Ebeam is already given via the jobcard
       smin = mN**2 ! GeV^2

       costmin=-1.0
       costmax= 1.0

       ymax = 1.0-1e-3
       xBmin = 0.0

       numin = 1e-3
!       numax=Ebeam*ymax

!       qsqmin = 0.0
       qsqmin = 0.01

    case (18) !=== E665

       Ebeam = 470.

       smin = 1.4 ! GeV^2

       costmin=-1.0
       costmax= 1.0

       ymax = 0.7
       xBmin = 0.0

       numin = 20.0
!       numax=Ebeam*ymax

       qsqmin = 0.1

    case default

       write(*,*) 'Case not implemented:',AccFlag,iExperiment
       call TRACEBACK()

    end select

    numax=Ebeam*ymax

    !***** checking for user-override: *****

    if (user_ymax>0) then
       write(*,'("OVERRIDE: ",A10," = ",0P,f13.5," --> ",f13.5)') "ymax",ymax,user_ymax
       ymax = user_ymax
       write(*,'("        : ",A10," = ",1P,g13.5," --> ",g13.5)') "numax",numax,Ebeam*ymax
       numax=Ebeam*ymax
    end if

    if (user_numin>0) then
       write(*,'("OVERRIDE: ",A10," = ",1P,e13.5," --> ",e13.5)') "numin",numin,user_numin
       numin = user_numin
    end if

    if (user_numax>0) then
       if (user_ymax>0) then
          write(*,*) "Error: Avoid setting 'user_numax' and 'user_ymax' simultanously!"
          call TRACEBACK()
       end if

       write(*,'("OVERRIDE: ",A10," = ",1P,e13.5," --> ",e13.5)') "numax",numax,user_numax
       numax = user_numax
       write(*,'("        : ",A10," = ",0P,f13.5," --> ",f13.5)') "ymax",ymax,numax/Ebeam
       ymax = numax/Ebeam
    end if

    if (user_costmin>-1.0) then
       write(*,'("OVERRIDE: ",A10," = ",0P,f13.5," --> ",f13.5)') "costmin",costmin,user_costmin
       costmin = user_costmin
    end if

    if (user_costmax<1.0) then
       write(*,'("OVERRIDE: ",A10," = ",1P,e13.5," --> ",e13.5)') "costmax",costmax,user_costmax
       costmax = user_costmax
    end if

    if (user_smin>0) then
       write(*,'("OVERRIDE: ",A10," = ",1P,e13.5," --> ",e13.5)') "smin",smin,user_smin
       smin = user_smin
    end if

    if (user_xBmin>0) then
       write(*,'("OVERRIDE: ",A10," = ",1P,e13.5," --> ",e13.5)') "xBmin",xBmin,user_xBmin
       xBmin = user_xBmin
    end if

    if (user_qsqmin>0) then
       write(*,'("OVERRIDE: ",A10," = ",1P,e13.5," --> ",e13.5)') "qsqmin",qsqmin,user_qsqmin
       qsqmin = user_qsqmin
    end if

    !***** write out some sample values, as given here:  *****
    write(*,*) "!   user_ymax    = ",ymax
    write(*,*) "!   user_numin   = ",numin
    write(*,*) "!   user_numax   = ",numax
    write(*,*) "!   user_costmin = ",costmin
    write(*,*) "!   user_costmax = ",costmax
    write(*,*) "!   user_smin    = ",smin
    write(*,*) "!   user_xBmin   = ",xBmin
    write(*,*) "!   user_qsqmin  = ",qsqmin
    write(*,*) "!   user_maxw    = ",maxw


    !***** Calculate the boundaries: *****

    Q2min=max(qsqmin,2.*Ebeam*(Ebeam-numax)*(1.-costmax))
    Q2max1=2.*ebeam*(Ebeam-numin)*(1.-costmin) !detector geometry
    Q2max2=2.*mN*numax+mN**2-smin !W cut
    Q2max=min(Q2max1,Q2max2)

    write(*,*)
    write(*,*) '   Ebeam          = ',Ebeam
    write(*,*) '   nu   (min,max) = ',numin,numax
    write(*,*) '   cost (min,max) = ',costmin,costmax
    write(*,*) '   Q2   (min,max) = ',Q2min, Q2max
!!$    write(*,*) '   Q2min          = ',Q2min
!!$    write(*,*) '   Q2min          = ',qsqmin,2.*Ebeam*(Ebeam-numax)*(1.-costmax)
!!$    write(*,*) '   Q2max          = ',Q2max
!!$    write(*,*) '   Q2max          = ',Q2max1,Q2max2
    write(*,*)

    if (costmin>costmax)  call TRACEBACK('costmin > costmax')
    if (Q2min.eq.0.0) call TRACEBACK('Q2min = 0')
    if (Q2min>Q2max)  call TRACEBACK('Q2min > Q2max')
    if (numin>numax)  call TRACEBACK('numin > numax')
    if (numin.eq.0.0)  call TRACEBACK('numin = 0')

    maxw = maxwArr(iDetector,iExperiment)
    if (user_maxw>0) then
       write(*,'("OVERRIDE: ",A10," = ",1P,e13.5," --> ",e13.5)') "maxw",maxw,user_maxw
       maxw = user_maxw
    end if

    write(*,*) '   maxw = ', maxw
    write(*,*)

  end subroutine initExperiment

  !****************************************************************************
  !****s* initHiLepton/initFixedKin
  ! NAME
  ! subroutine initFixedKin
  !
  ! PURPOSE
  ! set kinematical constraints according to no experiment but to given
  ! fixed kinematics. The namelist 'HiPhotonKinematics' is read.
  !
  ! NOTES
  ! you also have to set e.g. numin,numax etc. in order to ensure correct
  ! treatment of histograms
  !****************************************************************************
  subroutine initFixedKin

    use eN_eventDefinition, only: write_electronNucleon_event
    use eN_event, only: eNeV_GetKinV,eNev_SetProcess, &
         eNev_init_eWQ,eNev_init_exQ,eNev_init_enQ, &
         eNev_init_sWQ,eNev_init_sxQ,eNev_init_snQ, &
         eNev_init_BWQ,eNev_init_BnQ
    use ParamEP, only: CalcParamEP_R1990
    use output, only: Write_ReadingInput
    use callstack, only: traceback

    !**************************************************************************
    !****g* initFixedKin/nu
    ! SOURCE
    real :: nu = -99.9
    ! PURPOSE
    ! Photon energy [GeV]
    !**************************************************************************

    !**************************************************************************
    !****g* initFixedKin/Q2
    ! SOURCE
    real :: Q2 = -99.9
    ! PURPOSE
    ! transfer four momentum squared [GeV^2]
    !**************************************************************************

    !**************************************************************************
    !****g* initFixedKin/eps
    ! SOURCE
    real :: eps =-99.9
    ! PURPOSE
    ! Photon polarisation [1]
    !**************************************************************************

    !**************************************************************************
    !****g* initFixedKin/srts
    ! SOURCE
    real :: srts = -99.9
    ! PURPOSE
    ! sqrt(s) of electron nucleon system [GeV]
    !**************************************************************************

    !**************************************************************************
    !****g* initFixedKin/W
    ! SOURCE
    real :: W = -99.9
    ! PURPOSE
    ! sqrt(s) of photon nucleon system [GeV]
    !**************************************************************************

    !**************************************************************************
    !****g* initFixedKin/xBj
    ! SOURCE
    real :: xBj = -99.9
    ! PURPOSE
    ! Bjorken x [1]
    !**************************************************************************

    integer :: ios
    logical :: flagOK
    real :: R

    !**************************************************************************
    !****n* initHiLepton/HiPhotonKinematics
    ! NAME
    ! NAMELIST /HiPhotonKinematics/
    ! PURPOSE
    ! Namelist for initHiLepton in the case of iExperiment=0 includes:
    ! * nu
    ! * Q2
    ! * eps
    ! * srts
    ! * W
    ! * xBj
    ! * Ebeam
    ! NOTES
    ! you have to give a valid combination of three of them.
    !**************************************************************************
    NAMELIST /HiPhotonKinematics/ nu,Q2,eps,srts,W,xBj,Ebeam

    iDetector = 90
    Ebeam = -99.9

    call Write_ReadingInput('HiPhotonKinematics',0)
    rewind(5)
    read(5,nml=HiPhotonKinematics,iostat=ios)
    call Write_ReadingInput('HiPhotonKinematics',0,ios)

    call eNev_SetProcess(eNev_InitData, 1,1)  ! set to EM and electron

    flagOK = .false.
    if (eps .gt. 0) then
       if (W .gt. 0) then
          write(*,*) 'eNev_init_eWQ(eps,W,Q2)'
          write(*,*) '       ',eps,W,Q2
          call eNev_init_eWQ(eNev_InitData, eps,W,Q2, flagOK)

       else if (xBj .gt. 0) then
          write(*,*) 'eNev_init_exQ(eps,xBj,Q2)'
          write(*,*) '       ',eps,xBj,Q2
          call eNev_init_exQ(eNev_InitData, eps,xBj,Q2, flagOK)

       else if (nu .gt. 0) then
          write(*,*) 'eNev_init_enQ(eps,nu,Q2)'
          write(*,*) '       ',eps,nu,Q2
          call eNev_init_enQ(eNev_InitData, eps,nu,Q2, flagOK)

       else
          call TRACEBACK('you must provide W or xBj or nu!')
       end if
    else if (srts .gt. 0) then
       if (W .gt. 0) then
          write(*,*) 'eNev_init_sWQ(srts,W,Q2)'
          write(*,*) '       ',srts,W,Q2
          call eNev_init_sWQ(eNev_InitData, srts,W,Q2, flagOK)

       else if (xBj .gt. 0) then
          write(*,*) 'eNev_init_sxQ(srts,xBj,Q2)'
          write(*,*) '       ',srts,xBj,Q2
          call eNev_init_sxQ(eNev_InitData, srts,xBj,Q2, flagOK)

       else if (nu .gt. 0) then
          write(*,*) 'eNev_init_snQ(srts,nu,Q2)'
          write(*,*) '       ',srts,nu,Q2
          call eNev_init_snQ(eNev_InitData, srts,nu,Q2, flagOK)

       else
           call TRACEBACK('you must provide W or xBj or nu!')
       end if
    else if (Ebeam .gt. 0) then
       if (W .gt. 0) then
          write(*,*) 'eNev_init_BWQ(Ebeam,W,Q2)'
          write(*,*) '       ',Ebeam,W,Q2
          call eNev_init_BWQ(eNev_InitData, Ebeam,W,Q2, flagOK)

       else if (xBj .gt. 0) then
          write(*,*) 'eNev_init_BxQ(Ebeam,xBj,Q2)'
          write(*,*) '       ',Ebeam,xBj,Q2
!          call eNev_init_BxQ(eNev_InitData, Ebeam,xBj,Q2, flagOK)
          call TRACEBACK(' --- not implemented! STOP')

       else if (nu .gt. 0) then
          write(*,*) 'eNev_init_BnQ(Ebeam,nu,Q2)'
          write(*,*) '       ',Ebeam,nu,Q2
          call eNev_init_BnQ(eNev_InitData, Ebeam,nu,Q2, flagOK)

       else
          call TRACEBACK('you must provide W or xBj or nu!')
       end if
    else
       call TRACEBACK('you must provide eps or srts or Ebeam!')
    end if

    if (.not.flagOK) call TRACEBACK('kinematics not allowed!')

    call write_electronNucleon_event(eNev_InitData, .false., .true.)

    call eNeV_GetKinV(eNev_InitData, nu,Q2,W,eps=eps)
    Ebeam = eNev_InitData%lepton_in%momentum(0)
    write(*,*) '  nu   =',nu
    write(*,*) '  Q2   =',Q2
    write(*,*) '  W    =',W
    write(*,*) '  eps  =',eps
    write(*,*) '  Ebeam=',Ebeam
    call CalcParamEP_R1990(W,Q2,R)
    write(*,*) '  R1990=',R


    call Write_ReadingInput('HiPhotonKinematics',1)

    numin = nu - 0.1
    numax = nu + 0.1
    Q2min = Q2 - 0.1
    Q2max = Q2 + 0.1
    ymax = 1.0
    costmin = -1.0
    costmax =  1.0
    smin = 0.0


  end subroutine initFixedKin

  !****************************************************************************
  !****f* initHiLepton/GetiExperiment
  ! NAME
  ! integer function GetiExperiment()
  !
  ! PURPOSE
  ! return value of variable "iExperiment"
  !
  ! needed e.g. in analysis routines
  !****************************************************************************
  integer function GetiExperiment()
    GetiExperiment = iExperiment
  end function GetiExperiment

  !****************************************************************************
  !****s* initHiLepton/GetEnergies
  ! NAME
  ! real GetEnergies(Ebeam_,EIC_Ee_,EIC_EA_)
  !
  ! PURPOSE
  ! return value of variables Ebeam and (if given) EIC_Ee, EIC_EA
  !
  ! needed e.g. in analysis routines
  !****************************************************************************
  subroutine GetEnergies(Ebeam_,EIC_Ee_,EIC_EA_)
    real, intent(out) :: Ebeam_
    real, intent(out),OPTIONAL :: EIC_Ee_, EIC_EA_
    Ebeam_ = Ebeam
    if (present(EIC_Ee_)) EIC_Ee_=EIC_Ee
    if (present(EIC_EA_)) EIC_EA_=EIC_EA
  end subroutine GetEnergies


  !****************************************************************************
  !****s* initHiLepton/GetPhotonKin
  ! NAME
  ! subroutine GetPhotonKin(nu,Q2,W)
  !
  ! PURPOSE
  ! in the case of Experiment==0, this returns the given kinematics
  !
  ! needed e.g. in analysis routines
  !****************************************************************************
  subroutine GetPhotonKin(nu,Q2,W)
    use eN_event, only: eNeV_GetKinV
    real, intent(out):: nu,Q2,W
    call eNeV_GetKinV(eNev_InitData, nu,Q2,W)
  end subroutine GetPhotonKin


!!$  !*************************************************************************
!!$  !****s* initHiLepton/determineHiPhoton
!!$  ! NAME
!!$  ! subroutine determineHiPhoton(p,nu,Q2,eps,thetaLep,phiLep,W,Wfree,pcm,betacm)
!!$  !
!!$  ! PURPOSE
!!$  ! randomly generate all photon kinematics by Monte Carlo according flux
!!$  ! and detector constraints
!!$  !
!!$  ! INPUTS
!!$  ! * real,dimension(0:3) :: p -- momentum of target nucleon
!!$  !
!!$  ! OUTPUT
!!$  ! * real                :: nu,Q2,eps -- photon variables
!!$  ! * real                :: thetaLep,phiLep -- angles scattered lepton
!!$  ! * real                :: W,Wfree -- W of process
!!$  ! * real,dimension(0:3) :: pcm -- boost vector
!!$  ! * real,dimension(3)   :: betacm -- boost vector
!!$  !*************************************************************************
!!$  subroutine determineHiPhoton(p,nu,Q2,eps,thetaLep,phiLep,W,Wfree,pcm,betacm)
!!$
!!$!    use InitHiPhoton, only : determineSrtsfree
!!$    use minkowski, only : abs4
!!$    use CallStack
!!$    implicit none
!!$
!!$    real,dimension(0:3),intent(in)  :: p
!!$    real,               intent(out) :: nu,Q2,eps,thetaLep,phiLep,W, Wfree
!!$    real,dimension(0:3),intent(out) :: pcm
!!$    real,dimension(3),  intent(out) :: betacm
!!$
!!$    logical :: flagOK
!!$
!!$    flagOK=.false.
!!$    do while(.not.flagOK)
!!$       call photonFluxMC(p,nu,Q2,eps,thetaLep,phiLep)
!!$       call TRACEBACK("determineSrtsfree is replaced")
!!$!!!!       call determineSrtsfree(p,nu,Q2,Wfree,pcm,betacm,flagOK)
!!$!       if (Wfree.lt.HighEnergyThreshold) flagOK=.false.
!!$    end do
!!$
!!$    W = abs4( p + (/nu, 0.,0., sqrt(nu**2+Q2)/) )
!!$
!!$  end subroutine determineHiPhoton

  !****************************************************************************
  !****s* initHiLepton/photonFluxMC
  ! NAME
  ! subroutine photonFluxMC(eNev,pTarget)
  !
  ! PURPOSE
  ! Choose randomly some nu, Q2 and eps values corresponding to "iExperiment".
  ! Then calculate the flux for this given kinematics. Take the flux and
  ! some weight (called AccWeight) representing (averaged) probabilities
  ! of detecting the scattered lepton in the detector in order to accept
  ! or reject the choice accoring Monte Carlo techniques.
  !
  ! INPUTS
  ! * type(particle)  :: pTarget -- target nucleon to be considered
  !
  ! OUTPUT
  ! * type(electronNucleon_event) :: eNev -- The incoming electron
  !
  ! NOTES
  ! * since we changed the definition according the MC selection, the
  !   values given for "maxW" are not valid anymore.
  !   We have done some preliminary runs for some configurations to
  !   readjust he values. All non-checked values are set to zero.
  !   You may get a lot of "adjust maxW"-messages in
  !   first runs. See the file "initHiLep.AdjustMaxW.txt" for the largest
  !   value used so far, readjust the values in the code, recompile, restart.
  ! * depending on the flag "eN_event/restingNucleon" the flux is calculated
  !   as average over the fermi momentum or using the momentum of the
  !   target nucleon for fT and epsilon.
  !****************************************************************************
  subroutine photonFluxMC(eNev,pTarget)

    use particleDefinition
    use constants, only: pi, mN
    use eN_event, only: eNeV_GetKinV,eNev_SetProcess,eNev_init_BnQ, &
         eNeV_Set_PhiLepton
    use vector, only: theta_in
    use random, only: rn
    use output, only: DoPr

    type(particle), intent(in)  :: pTarget
    type(electronNucleon_event),intent(out) :: eNev

    real :: nu,Q2,W,W2, xB, thetaLepScatt,pLepScatt
    real :: wtest
    logical :: flagOK

    logical, save :: lWriteAdjustMaxW = .true.

    if (lWriteAdjustMaxW) then ! reset the value in the file to zero
       call WriteAdjustMaxW(-99.9)
       lWriteAdjustMaxW=.false.
    end if


    fluxLoop : do
       !in the following the lepton mass is neglected

       ! we select nu according f(nu)=c/nu

       nu = numin*(numax/numin)**rn()

       ! we select Q2 according f(Q2)=c/Q2
       ! Therefore we have to multiply wtest by Q2 (see below)
       ! Unfortunately we forgot to multiply also by log(Q2max/Q2min)=1/c
       ! in order to ensure the normalisation. If we repair this,
       ! we have to readjust maxw for all experiments.

       Q2 = Q2min*(Q2max/Q2min)**rn()

       W2 = 2.*mN*nu+mN**2-Q2

       if (W2.lt.smin) then !W-cut (includes xB<1 for smin>0.938**2)
!          write(*,*) 'W-cut: ', nu,Q2,sqrt(2.*mN*nu+mN**2-Q2)
          cycle fluxLoop
       end if

       if (iExperiment.eq.5) then
          if (W2.gt.10.24) cycle fluxloop ! experimental cut for JLab @ 5GeV
       end if

       xB=Q2/(2.*mN*nu)
       if (xB.lt.xBmin) then  !xB-cut
          cycle fluxLoop
       end if

       call eNev_SetProcess(eNev, 1,1)  ! set to EM and electron
       call eNev_init_BnQ(eNev, Ebeam,nu,Q2, flagOK)
       if (.not.flagOK) cycle fluxLoop
       call eNeV_Set_PhiLepton(eNev,rn()*2.*pi)

       eNev%nucleon%momentum = pTarget%momentum

       call eNeV_GetKinV(eNev, nu,Q2,W,fT=wtest)

!       call write_electronNucleon_event(eNev, .false., .true.)

       wtest = wtest * nu * log(numax/numin) ! because f(nu)=c/nu
       wtest = wtest * Q2 * log(Q2max/Q2min) ! because f(Q2)=1/Q2

       weight1 = wtest

       if (AccFlag) then
          pLepScatt=sqrt(Dot_product(eNev%lepton_out%momentum(1:3),eNev%lepton_out%momentum(1:3)))
          thetaLepScatt=theta_in(eNev%lepton_in%momentum(1:3),eNev%lepton_out%momentum(1:3),2)
          wtest=wtest*AccWeight(iDetector,0,pLepScatt,thetaLepScatt)
       end if

       weight2 = wtest

       if (wtest.gt.maxw) then
          if (DoPr(4)) write(*,*) 'in photonFluxMC: adjust maxw',wtest
          maxw = wtest
          call WriteAdjustMaxW(wtest)
       end if

       if (.not.ModusCalcFluxNorm) then
          if (wtest.lt.rn()*maxw) cycle fluxLoop
       end if


       exit fluxLoop
    end do fluxLoop



  contains

    subroutine WriteAdjustMaxW(r)
      use output, only: writeFileDocu

      real, intent(in) :: r

      real, save :: rmax = -99.9

      if (r > rmax) then
         open(121, file="initHiLep.AdjustMaxW.txt", status='unknown')
         rewind(121)
         write(121,*) r
         write(121,*)
         close(121)

         if (rmax < 0.0) then ! it's the first time
            call writeFileDocu('initHiLep.AdjustMaxW.txt',&
                 &'maximal value, please adjust input value')
         end if
         rmax = r + 1e-10

      end if

    end subroutine WriteAdjustMaxW

  end subroutine photonFluxMC




  !****************************************************************************
  !****f* initHiLepton/AccWeight
  ! NAME
  ! real function AccWeight(iDetector,part_type,momentum,theta)
  !
  ! PURPOSE
  ! Return the probability to detect the scattered lepton or produced
  ! hadrons with a given detector symmetry.
  !
  ! The detctor is rotated around the lepton beam axis and only the angle
  ! between particle vector and beam axis is taken into account.
  ! Therefore not only boolean values (1: detected, 0:not detected) are
  ! possible but all values inbetween.
  !
  ! INPUTS
  ! * integer :: iDetector -- choice of detector (see below)
  ! * integer :: part_type -- electron, pos./neg.hadron
  !     (only used for CLAS detector)
  ! * real    :: momentum -- particle momentum
  ! * real    :: theta -- angle of particle vector to beam axis
  !
  ! OUTPUT
  ! probability to detect the particle
  !
  ! NOTES
  ! The variable "iDetector" can take the values:
  ! * 1 = HERMES
  ! * 2 = EMC
  ! * 3 = CLAS
  ! * 4 = CLAS@5GeV
  ! * 5 = CLAS (mixture): electron: as 3, hadrons: as 4
  !
  ! The variable "part_type" has the values:
  ! *  0 = electron
  ! *  1 = hadron, positive charged
  ! * -1 = hadron, negative charged
  !
  !****************************************************************************
  real function AccWeight(iDetector,part_type,momentum,theta)

    use constants, only: pi
    use callstack, only: traceback

    integer,intent(in) :: iDetector, part_type
    real,intent(in):: momentum, theta

    real,parameter :: th1H=0.04, th2H=0.14, th3H=0.1746, th4H=0.2202
    real,parameter :: th1E=0.087, th2E=0.140, th3E=0.165
    real,parameter :: th1C =0.1047, th2C =2.4958 ! = 6°,143°
!    real,parameter :: th1Ce=0.2443, th2Ce=0.8727 ! = 14°,50°
    real,parameter :: th1Ce=0.2094, th2Ce=0.8727 ! = 12°,50°
    real,parameter :: pihalf=pi/2. ! PI/2
    real :: th

    th=abs(theta)

    select case (iDetector)

    case (1) ! ===== HERMES

       if (th.lt.th1H) then
          AccWeight = 0.
       else if (th.lt.th2H) then
          AccWeight = pihalf-asin(th1H/th)
       else if (th.lt.th3H) then
          AccWeight = asin(th2H/th)-asin(th1H/th)
       else if (th.lt.th4H) then
          AccWeight = asin(th2H/th)-acos(0.17/th)
       else
          AccWeight = 0.
       end if
       AccWeight = AccWeight/pihalf

    case (2) ! ===== EMC

       if (th.lt.th1E) then
          AccWeight = 1.
       else if (th.lt.th2E) then
          AccWeight = asin(th1E/th)
       else if (th.lt.th3E) then
          AccWeight = asin(th1E/th)-acos(th2E/th)
       else
          AccWeight = 0.
       end if
       AccWeight = AccWeight/pihalf

    case (3) ! ===== CLAS

       select case (part_type)
       case (0) ! --- electron
          if (th.lt.th1Ce) then
             AccWeight = 0.
          else if (th.lt.th2Ce) then
             AccWeight = 1.
          else
             AccWeight = 0.
          end if
       case (-1,1) ! --- hadron
          if (th.lt.th1C) then
             AccWeight = 0.
          else if (th.lt.th2C) then
             AccWeight = 1.
          else
             AccWeight = 0.
          end if
       end select

    case (4) ! ===== CLAS @ 5GeV

       AccWeight=CLAS(part_type,momentum,theta)

    case (5) ! ===== CLAS (mixture)

       select case (part_type)
       case (0) ! --- electron
          if (th.lt.th1Ce) then
             AccWeight = 0.
          else if (th.lt.th2Ce) then
             AccWeight = 1.
          else
             AccWeight = 0.
          end if
       case default
          AccWeight=CLAS(part_type,momentum,theta)
       end select

    case (0,90) ! ===== full acceptance

       AccWeight=1.0

    case default


       write(*,*) 'AccWeight: iDetector=',iDetector
       call TRACEBACK()

    end select

  contains

    !**************************************************************************
    !****f* AccWeight/CLAS
    ! NAME
    ! real FUNCTION CLAS(part_type,p,theta_rad)
    ! PURPOSE
    ! detector efficiency of CLAS @ 5 GeV
    ! NOTES
    ! * W.Brooks, private communication
    !**************************************************************************
    real FUNCTION CLAS(part_type,p,theta_rad)

      use constants, only: pi

      integer,intent(in) :: part_type
      real,intent(in) :: p,theta_rad
      REAL ::t_current,acc!,phi
      REAL :: t_max
      REAL :: phi0_el, phi0_nh, phi0_ph
      REAL :: theta0_el, theta0_nh, theta0_ph
      REAL :: thetas_el, thetas_nh, thetas_ph
      REAL :: p_shift, cel_ex, pel_ex
      REAL :: ch_ex,theta_cut
      REAL :: theta_min, theta_max,delta_phi, exp
      INTEGER :: electron,pos_hadron, neg_hadron
      REAL :: d2r,theta

      data t_max/3375./
      data phi0_el/30./
      data phi0_nh/25./
      data phi0_ph/25./
      data theta0_el/12.5/
      data theta0_nh/10./
      data theta0_ph/5./
      data thetas_el/15./
      data thetas_nh/15./
      data thetas_ph/25./
      data theta_max/50./
      data p_shift/0.15/
      data pel_ex/0.333/
      data cel_ex/0.35/
      data theta_cut/75./
      data ch_ex/0.3/
      data electron/0/
      data pos_hadron/1/
      data neg_hadron/-1/
      data d2r/0.01754533/

      data t_current/2250./

      theta=theta_rad*180./pi
      Acc = 0.0

      if (part_type.EQ.electron) then
         theta_min = theta0_el+thetas_el/(p*t_max/t_current+p_shift)
         if (theta.gt.theta_min.and.theta.lt.theta_max) then
            exp = cel_ex*(p*t_max/t_current)**pel_ex
            delta_phi = phi0_el*sin((theta-theta_min)*d2r)**exp
            Acc=delta_phi/30.     ! fraction of phi solid angle
         end if
      else if (part_type.EQ.pos_hadron) then
         theta_min=8.
         if (theta.gt.theta_min) then
            exp=(p*t_max/t_current/5.)**(1./8.)
            delta_phi = phi0_ph*sin((theta-theta_min)*d2r)**exp
            Acc=delta_phi/30.     ! fraction of phi solid angle
            if (p*(t_max/t_current).ge.1.and.theta.gt.139.) then
               Acc=0.
            else if (p*(t_max/t_current).lt.1.and.theta.gt.135.) then
               Acc=0.
            else if (p.lt.0.1) then ! this is for pions
               Acc=0.
               !c          else if (p.lt.0.3) then ! this is for protons
               !c             Acc=0.
            end if
         end if
      else if (part_type.EQ.neg_hadron) then
         theta_min = theta0_nh+thetas_nh/(p*t_max/t_current+p_shift)
         if (theta.gt.theta_min) then
            exp = ch_ex*(p*t_max/t_current)**pel_ex
            delta_phi = phi0_nh*sin((theta-theta_min)*d2r)**exp
            Acc=delta_phi/30.     ! fraction of phi solid angle
            if (theta.gt.139.) then
               Acc=0.
            end if
         end if
      end if

!c      vb1 = acc*sin(theta/57.3)/2. ! fraction of theta solid angle

!c     I comment out the above. I define the acceptance as the probability
!c     of an event in a particular bin getting accepted. It must range from
!c     0 to 1. It does not have to equal 4pi when integrated over.

      CLAS = acc

    end FUNCTION CLAS

  end function AccWeight

  !****************************************************************************
  !****f* initHiLepton/checkConservation
  ! NAME
  ! logical function checkConservation(eNev,outPart,channel)
  !
  ! PURPOSE
  ! Perform checks about conservation of:
  ! * Charge
  ! * Baryon number
  ! * Energy
  ! on the given event.
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev
  ! * type(particle),dimension(:) :: outPart
  ! * integer                     :: channel
  !
  ! OUTPUT
  ! * .true. if all checks okay
  !
  ! NOTES
  ! * if it was a 2p2h event, no checks are done. This could be improved,
  !   since we know all incoming particles stored in eNev. To be done!
  ! * To be done: Check what happens for DIS events !!!
  !****************************************************************************
  logical function checkConservation(eNev,outPart,channel)

    use particleDefinition
    use baryonPotentialModule, only: getNoPertPot_baryon
    use Electron_origin, only: origin_2p2hQE,origin_2p2hDelta
    use output, only: DoPr

    type(electronNucleon_event),intent(in)   :: eNev
    type(particle),dimension(:),intent(inout):: outPart
    integer,                    intent(in)   :: channel

    integer :: i
    type(particle)        :: outSUM
    integer               :: outSUM_nBaryon
    real, dimension (0:3) :: inMomentum

    checkConservation = .TRUE.

!    if (DoToyModel) then
!       return
!    endif

    ! if it was a 2p2h event, no checks are possible.
    if (channel.eq.origin_2p2hQE) return
    if (channel.eq.origin_2p2hDelta) return

    checkConservation = .FALSE.

    call setToDefault(outSUM)
    outSUM_nBaryon = 0

    inMomentum = eNev%boson%momentum + eNev%nucleon_free%momentum


    do i=1,size(outPart)
       if (outPart(i)%ID <= 0) cycle
       outSUM%charge  = outSUM%charge  + outPart(i)%charge
       outSUM%momentum= outSUM%momentum+ outPart(i)%momentum
       if (outPart(i)%ID .lt. 100) then
          if (outPart(i)%antiparticle) then
             outSUM_nBaryon = outSUM_nBaryon - 1
          else
             outSUM_nBaryon = outSUM_nBaryon + 1
          end if
       end if
    end do

    !...Charge:

    if (outSUM%charge .ne. eNev%nucleon_free%charge) then
       if (DoPr(2)) write(*,'(A,2i5)') 'Problem in initHiPhoton: charge not conserved',&
            outSUM%charge,eNev%nucleon%charge
       return
    end if

    !...Baryon number:

    if (outSUM_nBaryon .ne. 1) then
       if (DoPr(2)) write(*,'(A,2i5)') 'Problem in initHiPhoton: baryon number not conserved',&
            outSUM_nBaryon, eNev%nucleon%antiparticle
       return
    end if

    if ((.not.realRun).and.(getNoPertPot_baryon())) then

       !...energy:

       if ( abs(outSUM%momentum(0)-inMomentum(0)) .gt. 1e-2*inMomentum(0)) then
          if (DoPr(2)) write(*,'(A,1P,2e14.5)') 'Problem in initHiPhoton: energy not conserved',&
               outSUM%momentum(0),inMomentum(0)
          return
       end if

!!$       !...momentum:
!!$       do i=1,3
!!$          if ( abs(outSUM%momentum(i)-inMomentum(i)) .gt. 1e-2*inMomentum(i)) then
!!$             if (DoPr(2)) write(*,'(A,1P,2e14.5)') 'Problem in initHiPhoton: mom not conserved',&
!!$                  outSUM%momentum(i),inMomentum(i)
!!$             return
!!$          end if
!!$       end do

    end if

    checkConservation = .TRUE.

  end function checkConservation



end module initHiLepton
