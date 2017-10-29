!******************************************************************************
!****m* /HiLeptonAnalysis
! NAME
! module HiLeptonAnalysis
!
! PURPOSE
! Do the Analysis stuff for HiLepton
! INPUTS
! (no input)
! OUTPUT
! cf. "writeBinning"
! NOTES
! the old routines were replaced by "HistMP" stuff etc.
!******************************************************************************
module HiLeptonAnalysis

  use histf90
  use hist2Df90
  use histMPf90
  use histMC
  use AnaEventDefinition
  use AnaEvent
  use preEventDefinition
  use PreEvListDefinition
  use PreEvList
  use CollHistory
  use ValueRangeDefinition
  use InvMassesAnalysis

  implicit none
  private

  !****************************************************************************
  !****g* HiLeptonAnalysis/DoTimes
  ! SOURCE
  !
  logical, save :: DoTimes       = .false.
  !
  ! PURPOSE
  ! switch on/off: reporting of times
  !****************************************************************************

  !****************************************************************************
  !****g* HiLeptonAnalysis/DoOutChannels
  ! SOURCE
  !
  logical, save :: DoOutChannels = .false.
  !
  ! PURPOSE
  ! switch on/off: reporting of all final state channels
  !****************************************************************************

  !****************************************************************************
  !****g* HiLeptonAnalysis/DoInvMasses
  ! SOURCE
  !
  logical, save :: DoInvMasses   = .false.
  !
  ! PURPOSE
  ! switch on/off: reporting of pairwise-invariant-masses
  !****************************************************************************

  !****************************************************************************
  !****g* HiLeptonAnalysis/DoFindRho0
  ! SOURCE
  !
  logical, save :: DoFindRho0    = .false.
  !
  ! PURPOSE
  ! switch on/off: reconstructing rho0 from final pions
  !****************************************************************************

  !****************************************************************************
  !****g* HiLeptonAnalysis/DoClasie
  ! SOURCE
  !
  logical, save :: DoClasie      = .false.
  !
  ! PURPOSE
  ! switch on/off: Do pion analysis as Clasie et al., arXiv:0701.1481
  !****************************************************************************

  !****************************************************************************
  !****g* HiLeptonAnalysis/DoBrooks
  ! SOURCE
  !
  logical, save :: DoBrooks      = .false.
  !
  ! PURPOSE
  ! switch on/off: Do pi+ pT2 spectra for Brooks delta pT2
  !****************************************************************************

  !****************************************************************************
  !****g* HiLeptonAnalysis/DoMorrow
  ! SOURCE
  !
  logical, save :: DoMorrow      = .false.
  !
  ! PURPOSE
  ! switch on/off: Do pion analysis as Morrow et al., Morrow:2008ek
  !****************************************************************************

  !****************************************************************************
  !****g* HiLeptonAnalysis/DoMandelT
  ! SOURCE
  !
  logical, save :: DoMandelT      = .false.
  !
  ! PURPOSE
  ! switch on/off: Do pion analysis with Mandelstam t.
  !****************************************************************************

  !****************************************************************************
  !****g* HiLeptonAnalysis/DoCentralN
  ! SOURCE
  !
  logical, save :: DoCentralN     = .false.
  !
  ! PURPOSE
  ! switch on/off: Do centrality analysis with slow nucleons
  !****************************************************************************

  !****************************************************************************
  !****ig* HiLeptonAnalysis/DoEventAdd
  ! SOURCE
  !
  logical, save :: DoEventAdd  = .true.
  logical, save :: DoEventAdd2 = .true.
  !
  ! PURPOSE
  ! Decide, whether we have to do an EventArr analysis or not. These flags
  ! are not directly eccessible, but computed from other values:
  ! * DoEventAdd  = DoOutChannels.or.DoInvMasses.or.DoFindRho0
  ! * DoEventAdd2 = DoOutChannels.or.DoMorrow
  !****************************************************************************

  !****************************************************************************
  !****g* HiLeptonAnalysis/DoFSIsqrts
  ! SOURCE
  !
  logical, save :: DoFSIsqrts = .false.
  ! PURPOSE
  ! switch on/off: Estimate potential/future final state interactions
  !
  ! Plot the sqrt(s) distribution of potential final state interactions of
  ! perturbative particles with the nucleus (real) particles).
  ! (The interactions do not happen, this is calculated before every
  ! propagation.)
  ! In order to select the particle class for which one wants to report the
  ! FSI, please change directly the code.
  !****************************************************************************

  !****************************************************************************
  !****g* HiLeptonAnalysis/DoClassifyFirst
  ! SOURCE
  !
  logical, save :: DoClassifyFirst = .false.
  ! PURPOSE
  ! Classifying 'FirstEvent' into some classes
  !
  ! Needs DoEventAdd.
  !****************************************************************************

  !****************************************************************************
  !****g* HiLeptonAnalysis/DoLeptonKinematics
  ! SOURCE
  !
  logical, save :: DoLeptonKinematics = .false.
  ! PURPOSE
  ! switch on/off lepton kinematics output
  !****************************************************************************

  !****************************************************************************
  !****g* HiLeptonAnalysis/DoHadronKinematics
  ! SOURCE
  !
  logical, save :: DoHadronKinematics = .false.
  ! PURPOSE
  ! switch on/off hadron kinematics output
  !****************************************************************************

  !****************************************************************************
  !****it* HiLeptonAnalysis/tAverage
  ! NAME
  ! type tAverage
  ! PURPOSE
  ! Store general averages, as eq. <Q^2>, <W> etc.
  ! SOURCE
  !
  type tAverage
     real :: sum0,sum,W,Q2,nu
  end type tAverage
  !****************************************************************************

  logical, save :: FlagReadInput = .true.


  integer,parameter :: IDmin=-1, IDmax=21
  integer,parameter :: IZmin=-1, IZmax=1
  integer,parameter :: nAccWeight=2

  real, parameter :: add0 = 1e-20 ! dummy add for histograms

  real,save :: NLeptons

  type(histogramMP), save :: hMP_zH,hMP_nu,hMP_Q2,hMP_pT2 ! identified hadrons
  type(histogramMP), save :: hMP_pT2zH,hMP_pT2nu,hMP_pT2Q2,hMP_pT2pT2 ! <pT2> of identified hadrons

  type(histogram),   save :: hCH_zH,hCH_nu,hCH_Q2,hCH_pT2 ! charged hadrons

  type(histogram),   save :: hLep_Q2,hLep_nu
  type(histogram),dimension(2),save :: hKin_zH_Q2,hKin_nu_Q2 ! <Q2>(zH), <Q2>(nu)
  type(histogram),dimension(2),save :: hKin_zH_nu,hKin_nu_zH ! <nu>(zH), <zH>(nu)


  type(ValueRange) :: EhR,nuR,zHR,pTR,Q2R

  type(histogram2D), dimension(-1:1) :: H2DpTAvePion
  type(histogram2D), dimension(-1:1) :: H2DpTPionZH,H2DpTPionNU,H2DpTPionQ2

  type(histogram2D),dimension(0:9) :: H2DleptonXS

  type(histogramMP),dimension(3), save :: hMP_nu_zh,hMP_Q2_zh,hMP_pT2_zh
  type(histogramMP),dimension(3), save :: hMP_zH_nu,hMP_Q2_nu,hMP_pT2_nu
  type(histogramMP),dimension(2), save :: hMP_nu_pT,hMP_zH_pT,hMP_Q2_pT

  integer, save :: iExperiment ! store the value of the JobCard 'HiLeptonNucleus'
  integer, save :: iDetector

  type(tAverage), save :: GlobalAverage

  type(tAnaEvent), dimension(:,:), allocatable, save :: EventArr
  type(tAnaEvent), dimension(:,:), allocatable, save :: EventArr0

  type(tPreEvListEntry), dimension(:,:), allocatable, save, TARGET :: PreEvArr,PreEvArr0


  integer,save :: nRuns=0

  type(histogram), dimension(:,:), allocatable, save :: hJLAB5_R_nu,hJLAB5_R_Q2,hJLAB5_R_zH
  real, dimension(:,:,:,:), allocatable, save :: hJLAB5_pT2_ARR

  type(histogram2D), dimension(2), save :: H2D_CollHistPT

  type(histogramMP), save :: hMP_EIC_zH, hMP_EIC_nu, hMP_EIC_Q2
  type(histogramMP), dimension(1:5), save :: hMP_EIC_zH_Q2, hMP_EIC_nu_Q2

  type(histogram), save   :: HHH ! DUMMY for use by others

  !************ Histograms for 'DoTimes':

  type(histogram),  save :: HH1(1:3)
  type(histogram),  save :: HH2(0:7)
  type(histogramMP),save :: HHMP(1:2),HHMPlead(1:2)
  type(histogramMP),save :: HHMPnu(1:2),HHMPnulead(1:2)

  !************ Histograms for 'DoOutChannels':

  ! ---none---

  !************ Histograms for 'DoInvMasses':

  ! ---stored in the module---

  !************ Histograms for 'DoFindRho0':

  ! proc:
  ! * 2001 : VMD, 2002 : direct, 2003 : anomalous, 2004 : DIS, 2000 : misc

  type(histogram),  save :: hRho0MV, hRho0MX, hRho0DE, hRho0DecTime
  type(histogram2D), save :: H2D_rho0NuQ2,H2D_rho0NuQ2_proc(0:4)
  type(histogram2D), save :: H2D_rho0WQ2, H2D_rho0WQ2_proc(0:4)
  type(histogram), save :: hRho0zH, hRho0zH_proc(0:4)
  type(histogram), save :: hRho0t, hRho0t_proc(0:4)
  type(histogram), save :: hRho0sigmaW, hRho0sigmaW_proc(0:4)
  type(histogram), save :: hRho0Mom,hRho0W,hRho0Q2,hRho0Theta

  !************ Histograms for 'DoClasie':

  type(histogram),dimension(2),save :: histMX


  !************ Histograms for 'DoMorrow':

  type(histogramMC),dimension(3), save :: hMorrow
  type(histogramMC),dimension(2), save :: hMorrowT

  !************ Histograms for 'DoFSIsqrts':

  type(histogram),save :: hFSIsqrts

  !************ Histograms for 'DoBrooks':

  type(histogramMC), save :: hBrooks
  real, dimension(:,:,:,:), allocatable :: ArrBrooks

  !************ Histograms for 'DoClassifyFirst':

  type(histogramMP), dimension(0:3),  save :: hMP_nleadPT
  type(histogramMP), dimension(0:8),  save :: hMP_zHClass, hMP_pT2Class
  type(histogramMP), dimension(-1:10),  save :: hMP_zH_Generation

  !************ Histograms for 'DoMandelT':

  type(histogramMP), save :: hMP_MandelstamT

  !************ Histograms for 'DoCentralN':

  type(histogram2D), save :: H2D_CN_b(0:1), H2D_CN_bT(0:1), H2D_CN_bZ(0:1)
  type(histogram2D), save :: H2D_CN_bZE(0:1)
  type(histogram), save   :: H_CN_b, H_CN_bT, H_CN_bZ
  type(histogram2D), save :: H2D_CN_pTz(0:1)

  !****************************************************************************

!  integer, parameter :: iPartSet = 1 ! Set of Particles in HistMP
  integer, parameter :: iPartSet = 2 ! Set of Particles in HistMP

  character(3), dimension(-1:1), parameter :: piName  = (/'pi-','pi0','pi+'/)
  character(2), dimension(0:1),  parameter :: nucName = (/'N0','N+'/)

  public :: DoHiLeptonAnalysis
  public :: HiLeptonAnalysisPerTime

contains



  !****************************************************************************
  !****s* HiLeptonAnalysis/DoHiLeptonAnalysis
  ! NAME
  ! subroutine DoHiLeptonAnalysis(realPart,pertPart,MassNum,finalFlag,beforeRUN)
  !
  ! PURPOSE
  ! This is the routine, which does the actual analysis for HiLepton
  ! calculations
  !
  ! INPUTS
  ! * type(particle),dimension(:,:) :: realPart -- the real particle vector
  ! * type(particle),dimension(:,:) :: pertPart -- the perturbative particle vector
  ! * integer                       :: MassNum  -- the size of the nucleus
  ! * logical                       :: finalFlag -- flag, whether it is the final call
  ! * logical,             optional :: beforeRUN -- flag, whether this routine is
  !   called before all timesteps (i.e. directly after init) or at the end of a run
  !
  ! OUTPUT
  ! (none)
  !
  ! NOTES
  ! If this routine is called with beforeRUN=.true., then production and
  ! formation times of the particles produced in the first event are
  ! reported. Otherwise, after all timesteps, particle spectra etc are
  ! stored.
  !****************************************************************************
  subroutine DoHiLeptonAnalysis(realPart,pertPart,MassNum,finalFlag,beforeRUN)
    use particleDefinition
    use EventInfo_HiLep, only: EventInfo_HiLep_Get
    use initHiLepton, only: GetiExperiment,HiLepton_getRealRun
    use histf90
    use hist2Df90
    use collisionReporter

    use PIL_FormInfo
    use PIL_rho0Dec

    use Rho0Finder
    use particlePointerListDefinition
    use particlePointerList, only: ParticleList_CLEAR,ParticleList_INIT!,ParticleList_Print

    use output
    use CallStack
    use minkowski, only: abs4,abs4Sq

    use PauliBlockingModule, only: WriteBlockMom
    use idTable, only: pion

    type(particle),dimension(:,:),intent(in),   target :: realPart
    type(particle),dimension(:,:),intent(inout),target :: pertPart
    integer,intent(in) :: MassNum
    logical,intent(in) :: finalFlag
    logical,intent(in),optional :: beforeRUN

    integer :: i,j,iEns,nEns,iNuc,ii,nPart
    integer :: EventType
    real :: nu,Q2,epsilon,Weight
    real :: phi_Lepton,Ebeam


    logical,save :: startFlag=.true.

    integer,Allocatable,save :: nPertPart(:)
    integer, save       :: nPertPartMax = 0

    type(particle) :: Part
    type(particle), POINTER :: pPart

    real :: w1, w0
    integer :: iInt, iIntErr, iProc

    logical :: DoBeforeRUN

    type(tPreEvList), save :: PreEvList0,PreEvList
    type(tPreEvList), save :: PreEvListRhoOrig
!    type(tPreEvListEntry), save :: PreEvListEntry

    type(tParticleListNode), POINTER :: pNode


    real :: AccWeight, hWeight
    real :: zH
    real :: pT2   ! pT^2 according photon momentum axis
    integer :: nLead,h
    real :: W,Wcalc
    real :: dum
    real :: mul0 ! parameters for histograms
    integer :: iNr
    integer :: iQ2, ixB, izH
    real :: r

    type(particle),dimension(:,:),pointer :: pArr

    mul0 = 1./(max(nLeptons,1.0))

    DoBeforeRUN = .FALSE.
    if (present(beforeRUN)) DoBeforeRUN = beforeRUN

    write(*,*) '------------------- Do High Energy Lepton Analysis ---------------------'

    if (HiLepton_getRealRun()) then
       pArr => realPart
    else
       pArr => pertPart
    end if
    nEns = size(pertPart,dim=1)
    nPart = size(pArr,dim=2)

    ! **********initialization******************************
    if (startFlag) then
       write(*,*) '***HiLeptonAnalysis: Initializing...'
       startFlag=.false.

       if (FlagReadInput) call readInput()

       call WriteFileNames()

       nRuns=0
       iExperiment = getiExperiment()
       call resetBinning()

       allocate(nPertPart(nEns))

       if (DoTIMES) call DoTimes_INIT()

       call GlobalAverage_ZERO()

       if (DoEventAdd) then
          allocate(EventArr(nEns,MassNum))
          allocate(EventArr0(nEns,MassNum))
          allocate(PreEvArr(nEns,MassNum))
          allocate(PreEvArr0(nEns,MassNum))
          do i=1,nEns
             do j=1,MassNum
                call ParticleList_INIT(EventArr(i,j)%particleList)
                call event_CLEAR(EventArr(i,j))
                call ParticleList_INIT(EventArr0(i,j)%particleList)
                call event_CLEAR(EventArr0(i,j))
             end do
          end do
       end if


       if (DoInvMasses) call InvMasses_INIT

       if (DoClasie) then
!          call CreateHist(histMX(1),   'MissMass noFSI',      0.9*MassNum,1.0*MassNum,0.0005*MassNum)
!          call CreateHist(histMX(2),   'MissMass FSI',        0.9*MassNum,1.0*MassNum,0.0005*MassNum)
          dum = 0.92*MassNum
          call CreateHist(histMX(1),   'MissMass noFSI',      dum,dum+3.0,0.001)
          call CreateHist(histMX(2),   'MissMass FSI',        dum,dum+3.0,0.001)
       end if

       if (DoClassifyFirst) then
          do i=0,3
             call CreateHistMP(hMP_nleadPT(i), "d<pT2>(chgd.pi)/dzH "//trim(IntToChar(i)), 0.0, 1.2, 0.01, 2)
          end do
       end if

       if (DoFSIsqrts) call CreateHist(hFSIsqrts, 'sqrt(s) of possible FSI', 0.0,5.0,0.02)

       if (DoMorrow) then
          call CreateHistMC(hMorrow(1), 'invMass (pi+ pi-)', 0.0, 2.0, 0.01, 5)
          call CreateHistMC(hMorrow(2), 'invMass (p pi+)', 1.0, 3.0, 0.01, 5)
          call CreateHistMC(hMorrow(3), 'invMass (p pi-)', 1.0, 3.0, 0.01, 5)

          call CreateHistMC(hMorrowT(1), '-t (pi+ pi-)', 0.0, 3.5, 0.01, 5)
          call CreateHistMC(hMorrowT(2), '-t (rho0)', 0.0, 3.5, 0.01, 5)

          hMorrow%xDesc='Mass [GeV]'
          do i=1,3
             hMorrow(i)%yDesc=(/ "p rho0     ","Delta++ pi-",&
                  &              "Delta0 pi+ ","PS         ","Misc       " /)
          end do
          hMorrowT%xDesc='-t [GeV]'
          hMorrowT(1)%yDesc=(/ "p rho0     ","Delta++ pi-",&
                &              "Delta0 pi+ ","PS         ","Misc       " /)
          hMorrowT(2)%yDesc=(/ "VMD      ", "direct   ", &
               &               "anomalous", "DIS      ", "Misc     " /)

       end if

       if (DoBrooks) then
          call CreateHistMC(hBrooks, 'dN(pi+)/dpT2', 0.0,2.0,0.02, 8)
          hBrooks%xDesc='pT2 [GeV2]'
          hBrooks%yDesc=(/ "Q2=1-2, nu=2-3, z=0.5-0.6", &
               &           "Q2=1-2, nu=2-3, z=0.6-0.7", &
               &           "Q2=1-2, nu=3-4, z=0.5-0.6", &
               &           "Q2=1-2, nu=3-4, z=0.6-0.7", &
               &           "Q2=2-3, nu=2-3, z=0.5-0.6", &
               &           "Q2=2-3, nu=2-3, z=0.6-0.7", &
               &           "Q2=2-3, nu=3-4, z=0.5-0.6", &
               &           "Q2=2-3, nu=3-4, z=0.6-0.7" /)

          allocate( ArrBrooks(6,5,10,0:1) ) ! Q, xB, zH
          ArrBrooks = 0.0
       end if

       if (DoMandelT) then
          call CreateHistMP(hMP_MandelstamT,"dN_id/(N_e d|t|)", 0.0, 10.0, 0.02, iPartSet)
       end if

       if (DoCentralN) then
          do i=0,1
             call CreateHist2D(H2D_CN_b(i),"p vs b, "//nucName(i), &
                  & (/0.0,0.0/), (/5.0,2.0/), (/0.1,0.02/) )
             call CreateHist2D(H2D_CN_bT(i),"p vs bT, "//nucName(i), &
                  & (/0.0,0.0/), (/5.0,2.0/), (/0.1,0.02/) )
             call CreateHist2D(H2D_CN_bZ(i),"p vs bZ, "//nucName(i), &
                  & (/-5.0,0.0/), (/5.0,2.0/), (/0.1,0.02/) )
             call CreateHist2D(H2D_CN_bZE(i),"Ekin vs bZ, "//nucName(i), &
                  & (/-5.0,-0.1/), (/5.0,1.0/), (/0.1,0.005/) )

             call CreateHist2D(H2D_CN_pTz(i),"pT vs z, "//nucName(i), &
                  & (/0.0,0.0/), (/2.0,1.1/), (/0.02,0.02/) )
          end do
          call CreateHist(H_CN_b,"b",0.,5.,0.1)
          call CreateHist(H_CN_bT,"bT",0.,5.,0.1)
          call CreateHist(H_CN_bZ,"bZ",-5.,5.,0.1)
       end if

       write(*,*) '***HiLeptonAnalysis: Initializing... [END]'
    end if

    ! **********Do the actual Run **************************

    if (DoBeforeRUN) then
       if (DoTIMES) call DoTimes_CALC()

       do i=1,nEns
          do j=1,MassNum
             iNr = j
             if (HiLepton_getRealRun()) iNr = -1
             if (EventInfo_HiLep_Get(i,iNr,Weight,nu,Q2,epsilon,EventType)) then
                call DoBinningLepton(nu,Q2,Weight)

                call AddHist2D(H2DleptonXS(0), (/Q2,nu/), Weight) ! = total cross section
                if (EventType.lt.3000) then
                   call AddHist2D(H2DleptonXS(1+EventType/1000), (/Q2,nu/), Weight)
                   if (EventType.gt.2000.and.EventType.lt.3000) &
                        & call AddHist2D(H2DleptonXS(3+mod(EventType,1000)), (/Q2,nu/), Weight)
                end if
                call GlobalAverage_ADD(Weight, Q2, nu)

             end if
             if (HiLepton_getRealRun()) exit
          end do
       end do

       if (DoEventAdd) then

          ! repeating here 3 times the double loop structure
          ! is necessary: (iEns,iNuc) .ne. (i,j) !!!

          do i=1,nEns
             do j=1,MassNum
                call event_CLEAR(EventArr0(i,j))
             end do
          end do

          do i=1,nEns
             do j=1,size(pArr,dim=2)
                if (pArr(i,j)%Id <  0) exit
                if (pArr(i,j)%Id <= 0) cycle

                if (pArr(i,j)%lastCollisionTime <  0) cycle

                iEns = pArr(i,j)%firstEvent / 1000
                iNuc = mod(pArr(i,j)%firstEvent,1000)
                pPart => pArr(i,j)
                call event_ADD(EventArr0(iEns,iNuc), pPart)
             end do
          end do

          do i=1,nEns
             do j=1,MassNum
                if (CreateSortedPreEvent(EventArr0(i,j),PreEvArr0(i,j)%preE)) then
                   PreEvArr0(i,j)%weight = EventArr0(i,j)%particleList%first%V%perweight
                   if (DoOutChannels) call PreEvList_INSERT(PreEvList0,PreEvArr0(i,j))
                end if
             end do
          end do

          if (DoOutChannels) then
             open(141,file='OutChannels.INIT.dat', status='unknown')
             rewind(141)
             call PreEvList_Print(141,PreEvList0,mul0,withLN=.true.)
             close(141)
          end if
       end if

       do i=1,nEns
          do j=1,nPart
             pPart => pArr(i,j)
             if (pPart%Id <  0) exit
             if (pPart%Id <= 0) cycle

             iNr = pPart%firstEvent
             if (HiLepton_getRealRun()) iNr = -1

             if (.not.(EventInfo_HiLep_Get(i,iNr,Weight,nu,Q2,epsilon,EventType))) &
                  & exit
             if (.not.(EventInfo_HiLep_Get(i,iNr,Weight,nu,Q2,epsilon,EventType))) &
                  & call TRACEBACK()

             zH = pPart%momentum(0)/nu
             pT2 = pPart%momentum(1)**2 + pPart%momentum(2)**2

             if (DoClassifyFirst) then
                nLead = 0
                if (PIL_FormInfo_GET( pPart%number,h)) nLead=mod(h/100,10)
                call AddHistMP(hMP_nleadPT(nLead),pPart,&
                     & zH,pPart%perweight,pPart%perweight*pT2)
             end if

             if (DoClasie) then
                if ((pPart%Id==pion).and.(pPart%charge==1)) then
                   nRuns=nRuns+1 ! temporary
                   call CalculateExclPiP(1)
                   nRuns=nRuns-1 ! temporary, see above
                end if
             end if

             if (DoFSIsqrts) call CalcFSIsqrts()

          end do
       end do

       return
    end if



    nRuns=nRuns+1 ! number of runs with same energy

    if (DoEventAdd) then
       do i=1,nEns
          do j=1,MassNum
             call event_CLEAR(EventArr(i,j))
          end do
       end do
    end if


    nPertPart = 0

    ensembleLoop : do i=1,nEns
       indexLoop : do j=1,nPart
          pPart => pArr(i,j)
          if (pPart%Id <  0) exit
          if (pPart%Id <= 0) cycle

          iNr = pPart%firstEvent
          if (HiLepton_getRealRun()) iNr = -1

          nPertPart(i) =  nPertPart(i)+1

          if (.not.(EventInfo_HiLep_Get(i,iNr,Weight,nu,Q2,epsilon,EventType,&
               & Ebeam,phi_Lepton))) then

             write(*,*) 'not a photon induced particle'
             call WriteParticle(6)
             call WriteParticle(6,i,j,pPart)
             if (HiLepton_getRealRun()) then
                cycle
             else
                call TRACEBACK()
             end if
          end if
          if (abs(pPart%perWeight-Weight).gt.1e-8) then
             write(*,*) 'Error: perWeight,Weight=',pPart%perWeight,Weight
             write(*,*)
             write(*,*) nu,Q2,epsilon,EventType,Ebeam,phi_Lepton
             call WriteParticle(6)
             call WriteParticle(6,i,j,pPart)
             call TRACEBACK()
          end if

          call DoBinningHadron(pPart,nu,Q2,Ebeam,phi_Lepton,i,j)

          if (pPart%lastCollisionTime <  0) cycle

          if (DoEventAdd) then
             iEns = pPart%firstEvent / 1000
             iNuc = mod(pPart%firstEvent,1000)
             call event_ADD(EventArr(iEns,iNuc), pPart)
          end if

          if (DoClasie) then
             if ((pPart%Id==pion).and.(pPart%charge==1)) then
                call CalculateExclPiP(2)
             end if
          end if

       end do indexLoop

       if (nPertPart(i).gt.nPertPartMax) then
          nPertPartMax = nPertPart(i)
       end if

    end do ensembleLoop

    if (DoFindRho0) then
       call FindRho0Analysis()
    end if

    if (DoMandelT) then
       call MandelTAnalysis()
    end if

    if (DoCentralN) then
       call CentralNAnalysis()
    end if


    call writeBinning()

    if (DoTIMES) call DoTimes_Write()

    write(*,*) 'Number of Particles per Ensemble: [initialized:',&
         & size(pArr,dim=2),']'
    write(*,*) '                                      Maximum :',&
         & nPertPartMax


    if (DoInvMasses) then
       do i=1,nEns
          do j=1,MassNum
             if (EventInfo_HiLep_Get(i,j,Weight,nu,Q2,epsilon,EventType)) then
!               write(*,'("#### EVENT ",i5,i3," ##### ",i7)') i,j,EventType
!               call ParticleList_Print(EventArr(i,j))
                call InvMasses_FillEvent(EventArr(i,j))
             end if
          end do
       end do
       call InvMasses_Write (mul0)
    end if

    if (DoEventAdd2) then
       do i=1,nEns
          do j=1,MassNum
             if (CreateSortedPreEvent(EventArr(i,j),PreEvArr(i,j)%preE)) then
                PreEvArr(i,j)%weight = EventArr(i,j)%particleList%first%V%perweight
                if (DoOutChannels) call PreEvList_INSERT(PreEvList,PreEvArr(i,j))
             end if
          end do
       end do
    end if

    if (DoOutChannels) then
       open(141,file='OutChannels.FINAL.dat', status='unknown')
       rewind(141)
       call PreEvList_Print(141,PreEvList,mul0)
       close(141)
    end if

    if (DoMorrow) then
       call MorrowAnalysis()
       do i=1,3
          call WriteHistMC(hMorrow(i), add=add0,mul=mul0,&
               & file='HiLep.Morrow.'//trim(intToChar(i))//'.dat')
       end do
       do i=1,2
          call WriteHistMC(hMorrowT(i), add=add0,mul=mul0,&
               & file='HiLep.MorrowT.'//trim(intToChar(i))//'.dat')
       end do
    end if

    if (DoBrooks) call BrooksAnalysis()


    if (DoFindRho0) then
       call WriteHist(hRho0MV,add=add0,mul=mul0,file='HiLep.Rho0_MV.dat',dump=.true.)
       call WriteHist(hRho0MX,add=add0,mul=mul0,file='HiLep.Rho0_MX.dat',dump=.true.)
       call WriteHist(hRho0DE,add=add0,mul=mul0,file='HiLep.Rho0_DE.dat',dump=.true.)
       call WriteHist(hRho0DecTime,add=add0,mul=mul0,file='HiLep.Rho0_DecTime.dat',dump=.true.)
       call WriteHist(hRho0Mom,add=add0,mul=mul0,file='HiLep.Rho0_Mom.dat',dump=.true.)
       call WriteHist(hRho0Theta,add=add0,mul=mul0,file='HiLep.Rho0_Theta.dat',dump=.true.)

       open(141,file='HiLep.OutChannels.RhoOrig.dat', status='unknown')
       rewind(141)
       call PreEvList_Print(141,PreEvListRhoOrig,mul0)
       close(141)

       ! Due to historical reasons, the following file names had a "JLAB" in it.
       ! Some analysis tools may be confused by this.
       ! You may replace "Rho0_" instead of "JLABrho" and vice versa.


       ! --- nu Q2 ---

       call WriteHist2D_Gnuplot(h2D_rho0NuQ2, mul=mul0, add=add0,&
            & file='HiLep.Rho0_NuQ2.dat',dump=.true.)
       call WriteHist2D_Gnuplot(h2D_rho0NuQ2, DoAve=.true.,&
            & file='HiLep.Rho0_NuQ2.AveMis.dat',dump=.false.)
       call IntegrateHist2D(h2D_rho0NuQ2,HHH,2)
       call WriteHist(HHH,add=add0,mul=mul0,&
            & file="HiLep.Rho0_NuQ2.Int2.dat",dump=.true.)
       call AverageHist2D(h2D_rho0NuQ2,HHH,2)
       call WriteHist(HHH,file="HiLep.Rho0_NuQ2.Ave2.dat",dump=.true.)

       do i=0,4
          call WriteHist2D_Gnuplot(h2D_rho0NuQ2_proc(i), mul=mul0, add=add0,&
            & file='HiLep.Rho0_NuQ2_proc'//trim(intToChar(i))//'.dat',dump=.true.)
          call WriteHist2D_Gnuplot(h2D_rho0NuQ2_proc(i), DoAve=.true.,&
            & file='HiLep.Rho0_NuQ2_proc'//trim(intToChar(i))//'.AveMis.dat',dump=.false.)
          call IntegrateHist2D(h2D_rho0NuQ2_proc(i),HHH,2)
          call WriteHist(HHH,add=add0,mul=mul0,&
               & file='HiLep.Rho0_NuQ2_proc'//trim(intToChar(i))//'.Int2.dat',dump=.true.)
          call AverageHist2D(h2D_rho0NuQ2_proc(i),HHH,2)
          call WriteHist(HHH,&
               & file='HiLep.Rho0_NuQ2_proc'//trim(intToChar(i))//'.Ave2.dat',dump=.true.)
       end do

       ! --- W Q2 ---

       call WriteHist(hRho0W,add=add0,mul=mul0,&
            & file="HiLep.Rho0_W.dat",dump=.true.)
       call WriteHist(hRho0Q2,add=add0,mul=mul0,&
            & file="HiLep.Rho0_Q2.dat",dump=.true.)

       call WriteHist2D_Gnuplot(h2D_rho0WQ2, mul=mul0, add=add0,&
            & file='HiLep.Rho0_WQ2.dat',dump=.true.)
       call WriteHist2D_Gnuplot(h2D_rho0WQ2, DoAve=.true.,&
            & file='HiLep.Rho0_WQ2.AveMis.dat',dump=.false.)
       call IntegrateHist2D(h2D_rho0WQ2,HHH,1)
       call WriteHist(HHH,add=add0,mul=mul0,&
            & file="HiLep.Rho0_WQ2.Int1.dat",dump=.true.)
       call AverageHist2D(h2D_rho0WQ2,HHH,1)
       call WriteHist(HHH,file="HiLep.Rho0_WQ2.Ave1.dat",dump=.true.)

       call IntegrateHist2D(h2D_rho0WQ2,HHH,2)
       call WriteHist(HHH,add=add0,mul=mul0,&
            & file="HiLep.Rho0_WQ2.Int2.dat",dump=.true.)
       call AverageHist2D(h2D_rho0WQ2,HHH,2)
       call WriteHist(HHH,file="HiLep.Rho0_WQ2.Ave2.dat",dump=.true.)

       do i=0,4
          call WriteHist2D_Gnuplot(h2D_rho0WQ2_proc(i), mul=mul0, add=add0,&
            & file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.dat',dump=.true.)
          call WriteHist2D_Gnuplot(h2D_rho0WQ2_proc(i), DoAve=.true.,&
            & file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.AveMis.dat',dump=.false.)

          call IntegrateHist2D(h2D_rho0WQ2_proc(i),HHH,1)
          call WriteHist(HHH,add=add0,mul=mul0,&
               & file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.Int1.dat',dump=.true.)
          call AverageHist2D(h2D_rho0WQ2_proc(i),HHH,1)
          call WriteHist(HHH,&
               & file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.Ave1.dat',dump=.true.)

          call IntegrateHist2D(h2D_rho0WQ2_proc(i),HHH,2)
          call WriteHist(HHH,add=add0,mul=mul0,&
               & file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.Int2.dat',dump=.true.)
          call AverageHist2D(h2D_rho0WQ2_proc(i),HHH,2)
          call WriteHist(HHH,&
               & file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.Ave2.dat',dump=.true.)
       end do

       ! --- zH ---
       ! --- t  ---

       call WriteHist(hRho0zH,mul=mul0, add=add0,&
            & file='HiLep.Rho0_zH.dat',dump=.true.)
       call WriteHist(hRho0t,mul=mul0, add=add0,&
            &file='HiLep.Rho0_t.dat',dump=.true.)
       do i=0,4
          call WriteHist(hRho0zH_proc(i),mul=mul0, add=add0,&
               &file='HiLep.Rho0_zH_proc'//trim(intToChar(i))//'.dat',dump=.true.)
          call WriteHist(hRho0t_proc(i),mul=mul0, add=add0,&
               &file='HiLep.Rho0_t_proc'//trim(intToChar(i))//'.dat',dump=.true.)
       end do

       ! --- sigmaW ---

       call WriteHist(hRho0sigmaW,mul=mul0, add=add0,&
            &file='HiLep.Rho0_sigmaW.dat',dump=.true.)
       call WriteHist(hRho0sigmaW,mul=mul0, add=add0, DoAve=.true.,&
            &file='HiLep.Rho0_sigmaW.Ave.dat')
       do i=0,4
          call WriteHist(hRho0sigmaW_proc(i),mul=mul0, add=add0,&
               &file='HiLep.Rho0_sigmaW_proc'//trim(intToChar(i))//'.dat',dump=.true.)
          call WriteHist(hRho0sigmaW_proc(i),mul=mul0, add=add0, DoAve=.true.,&
               &file='HiLep.Rho0_sigmaW_proc'//trim(intToChar(i))//'.Ave.dat')
       end do

    end if

!    call DoingCollHistAnalysis

    if (DoClassifyFirst) then
       do i=0,3
          call WriteHistMP(hMP_nleadPT(i), add=add0,mul=mul0,iColumn=1,&
               & file='HiLep.nleadPT.N.'//trim(intToChar(i))//'.dat',dump=.true.) ! noAcc
          call WriteHistMP(hMP_nleadPT(i), DoAve=.true.,&
               & file='HiLep.nleadPT.PT2.'//trim(intToChar(i))//'.dat') ! noAcc
       end do
    end if

    if (DoClasie) then
       call WriteHist(histMX(1),mul=mul0,add=add0,file='HiLep.Clasie.1.dat',dump=.true.)
       call WriteHist(histMX(2),mul=mul0,add=add0,file='HiLep.Clasie.2.dat',dump=.true.)
    end if

    if (DoBrooks) then
       call WriteHistMC(hBrooks,'HiLep.Brooks.dat',add=add0,mul=mul0,dump=.true.)
       open(141,file='HiLep.Brooks.detailed.dat', status='unknown')
       rewind(141)
       write(141,'(A)') "# <pT2> for Q2,xB,zH bins"
       write(141,'(A)') "#    Q2 = 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4"
       write(141,'(A)') "#    xB = 0.1-0.19, 0.19-0.28, 0.28-0.37, 0.37-0.46, 0.46-0.55"
       write(141,'(A)') "#    zH = 0-0.1, ... 0.9-1"
       write(141,'(A)') "#"
       write(141,'(A)') "# iQ2,ixB,izH, <pT2>, sum(w), sum(pT2*w)"

       do iQ2=1,6
          do ixB=1,3
             do izH=1,10
                if (ArrBrooks(iQ2,ixB,izH,0)>0) then
                   r = ArrBrooks(iQ2,ixB,izH,1)/ArrBrooks(iQ2,ixB,izH,0)
                else
                   r = 99.9
                end if
                write(141,'(3i3.0,1P,3e13.4)') iQ2,ixB,izH, r, ArrBrooks(iQ2,ixB,izH,:)
             end do
          end do
       end do
       close(141)
    end if

    if (DoFSIsqrts) then
       call WriteHist(hFSIsqrts,mul=mul0,add=add0,file='HiLep.FSIsqrts.dat',dump=.true.)
    end if

    if (DoMandelT) then
       ! ===== Mandelstam-t-spectra =====

       call WriteHistMP(hMP_MandelstamT,add=add0,mul=mul0,iColumn=3,&
            &file='HiLep.MandelstamT.idH.excl.dat')
       call WriteHistMP(hMP_MandelstamT,add=add0,mul=mul0,iColumn=1,&
            &file='HiLep.MandelstamT.idH.incl.dat',dump=.true.)

    end if

    if (DoCentralN) then
       call WriteHist(H_CN_b, mul=mul0, add=add0,&
            & file='HiLep.CentralN_b.C_CN_b.dat',dump=.true.)
       call WriteHist(H_CN_bT, mul=mul0, add=add0,&
            & file='HiLep.CentralN_b.C_CN_bT.dat',dump=.true.)
       call WriteHist(H_CN_bZ, mul=mul0, add=add0,&
            & file='HiLep.CentralN_b.C_CN_bZ.dat',dump=.true.)

       do i=0,1
          call WriteHist2D_Gnuplot(H2D_CN_b(i), mul=mul0, add=add0,&
               & file='HiLep.CentralN_b.'//nucName(i)//'.dat',dump=.true.)
          call WriteHist2D_Gnuplot(H2D_CN_b(i), MaxVal=0.0,H2=H_CN_b,&
               & file='HiLep.CentralN_b.norm.'//nucName(i)//'.dat')

          call WriteHist2D_Gnuplot(H2D_CN_bT(i), mul=mul0, add=add0,&
               & file='HiLep.CentralN_bT.'//nucName(i)//'.dat',dump=.true.)
          call WriteHist2D_Gnuplot(H2D_CN_bT(i), MaxVal=0.0,H2=H_CN_bT,&
               & file='HiLep.CentralN_bT.norm.'//nucName(i)//'.dat')

          call WriteHist2D_Gnuplot(H2D_CN_bZ(i), mul=mul0, add=add0,&
               & file='HiLep.CentralN_bZ.'//nucName(i)//'.dat',dump=.true.)
          call WriteHist2D_Gnuplot(H2D_CN_bZ(i), MaxVal=0.0,H2=H_CN_bZ,&
               & file='HiLep.CentralN_bZ.norm.'//nucName(i)//'.dat')
          call IntegrateHist2D(H2D_CN_bZ(i),HHH,2)
          call WriteHist(HHH, mul=mul0, add=add0,&
               & file='HiLep.CentralN_bZ.'//nucName(i)//'.y.dat',dump=.true.)

          call WriteHist2D_Gnuplot(H2D_CN_bZE(i), mul=mul0, add=add0,&
               & file='HiLep.CentralN_bZ_Ekin.'//nucName(i)//'.dat',dump=.true.)
          call WriteHist2D_Gnuplot(H2D_CN_bZE(i), MaxVal=0.0,H2=H_CN_bZ,&
               & file='HiLep.CentralN_bZ_Ekin.norm.'//nucName(i)//'.dat')
          call IntegrateHist2D(H2D_CN_bZE(i),HHH,2)
          call WriteHist(HHH, mul=mul0, add=add0,&
               & file='HiLep.CentralN_bZ_Ekin.'//nucName(i)//'.y.dat',dump=.true.)

          call WriteHist2D_Gnuplot(H2D_CN_pTz(i), mul=mul0, add=add0,&
               & file='HiLep.CentralN.pTz.'//nucName(i)//'.dat',dump=.true.)
       end do
    end if

    call WriteBlockMom(mul0)

    call cR_Write(nRuns)

    call GlobalAverage_Write(6, nRuns,nEns,massNum)
    open(141,file='HiLep.GlobalAverages.dat', status='unknown')
    rewind(141)
    call GlobalAverage_Write(141, nRuns,nEns,massNum)
    close(141)

    write(*,*) 'Histograms written.'

    if (finalFlag) then
       startFlag=.true.
       if (DoTIMES) call DoTimes_Finish()
    end if


  contains

    !**************************************************************************

    subroutine DoTimes_INIT()

      use initHiLepton, only: GetEnergies

      real :: Ebeam

      call CreateHist(HH1(1), 'time 1'   ,0.0,1.1,0.02)
      call CreateHist(HH1(2), 'time 2'   ,0.0,1.1,0.02)
      call CreateHist(HH1(3), 'timediff' ,0.0,1.1,0.02)

      call CreateHist(HH2(0), 'ALL'   ,0.0,1.1,0.02)
      call CreateHist(HH2(1), 'QQ string '   ,0.0,1.1,0.02)
      call CreateHist(HH2(2), 'QgQ string'   ,0.0,1.1,0.02)
      call CreateHist(HH2(3), 'gg string '   ,0.0,1.1,0.02)
      call CreateHist(HH2(4), 'Cluster1  '   ,0.0,1.1,0.02)
      call CreateHist(HH2(5), 'Cluster2  '   ,0.0,1.1,0.02)
      call CreateHist(HH2(6), 'DokuLine  '   ,0.0,1.1,0.02)
      call CreateHist(HH2(7), 'MISC      '   ,0.0,1.1,0.02)

      call CreateHistMP(HHMP(1), 'time 1 vs z'   ,0.0,1.1,0.02, 2)
      call CreateHistMP(HHMP(2), 'time 2 vs z'   ,0.0,1.1,0.02, 2)
      call CreateHistMP(HHMPlead(1), 'time 1 vs z (only lead)' ,0.0,1.1,0.02, 2)
      call CreateHistMP(HHMPlead(2), 'time 2 vs z (only lead)' ,0.0,1.1,0.02, 2)

      call GetEnergies(Ebeam)

      call CreateHistMP(HHMPnu(1), 'time 1 vs nu'   ,0.0,Ebeam,Ebeam/200, 2)
      call CreateHistMP(HHMPnu(2), 'time 2 vs nu'   ,0.0,Ebeam,Ebeam/200, 2)

      call CreateHistMP(HHMPnulead(1), 'time 1 vs nu (only lead)', &
           &0.0,Ebeam,Ebeam/200, 2)
      call CreateHistMP(HHMPnulead(2), 'time 2 vs nu (only lead)', &
           &0.0,Ebeam,Ebeam/200, 2)

    end subroutine DoTimes_INIT

    !**************************************************************************

    subroutine DoTimes_CALC()

      integer :: internalN
!      real :: gamma
      real :: rDum1, rDum2

      real, parameter :: zH1=0.6,zH2=0.8

      do i=1,nEns
         do j=1,size(pertPart,dim=2)
            if (pertPart(i,j)%ID <  0) exit
            if (pertPart(i,j)%ID <= 0) cycle
            Part = pertPart(i,j)

            rDum1 = Part%productionTime
            rDum2 = Part%formationTime

            if (EventInfo_HiLep_Get(i,Part%firstEvent,Weight,nu,Q2,epsilon,EventType,Ebeam)) then
            else
               write(*,*) 'Ooops,L111'
            end if

            zH = Part%momentum(0)/nu
            !velo = sqrt(Part%velocity(1)**2+Part%velocity(2)**2+Part%velocity(3)**2)

!!$            gamma = Part%momentum(0)/Part%mass

!!$            rDum1=rDum1/gamma
!!$            rDum2=rDum2/gamma

!!$            rDum1=rDum1/Part%momentum(0)
!!$            rDum2=rDum2/Part%momentum(0)

            w0 = Part%perweight

            if (rDum2 .ge. 0.0) then
               call AddHist(hh1(1),zH,w0,w0*rDum1)
               call AddHist(hh1(2),zH,w0,w0*rDum2)
               call AddHist(hh1(3),zH,w0,w0*(rDum2-rDum1))

               call AddHistMP(hhMP(1),Part,zH,w0,w0*rDum1)
               call AddHistMP(hhMP(2),Part,zH,w0,w0*rDum2)

               if ((zH.ge.zH1).and.(zH.le.zH2)) then
                  call AddHistMP(hhMPnu(1),Part,nu,w0,w0*rDum1)
                  call AddHistMP(hhMPnu(2),Part,nu,w0,w0*rDum2)
               end if

               if (rDum1.le.0.0) then
                  call AddHistMP(hhMPlead(2),Part,zH,w0,w0*rDum2)
                  if ((zH.ge.zH1).and.(zH.le.zH2)) then
                     call AddHistMP(hhMPnulead(2),Part,nu,w0,w0*rDum2)
                  end if
               end if
            end if


            if (PIL_FormInfo_GET(Part%number,internalN)) then
               ! write(*,*) 'ok'
            else
               ! write(*,*) 'fail'
            end if

            iInt = mod(internalN/10,10)
            if (iInt==0) iInt = 7
            iIntErr = mod(internalN,10)

            w1 = 0.0
            if (iIntErr.ne.0) w1=w0

            call AddHist(hh2(0),hh2(iInt),zH,w0,w1)

         end do
      end do

    end subroutine DoTimes_CALC

    !**************************************************************************

    subroutine DoTimes_Write()
      use output

      do i=1,3
         call WriteHist(hh1(i),add=add0,mul=1./nRuns,DoAve=.true.,&
              &file='DoTIMES.'//trim(intToChar(140+i))//'.dat')
      end do

      do i=0,7
         call WriteHist(hh2(i),add=add0,mul=1./nRuns,&
              &file='DoTIMES.'//trim(intToChar(150+i))//'.dat')
      end do

      do i=1,2
         ! write averages:
         call WriteHistMP(hhMP(i),DoAve=.true.,&
              &file='DoTIMES.MP'//trim(intToChar(140+i))//'.dat')
         call WriteHistMP(hhMPnu(i),DoAve=.true.,&
              &file='DoTIMES.MPnu'//trim(intToChar(140+i))//'.dat')

         call WriteHistMP(hhMPlead(i),DoAve=.true.,&
              &file='DoTIMES.MP'//trim(intToChar(140+i))//'.lead.dat')
         call WriteHistMP(hhMPnulead(i),DoAve=.true.,&
              &file='DoTIMES.MPnu'//trim(intToChar(140+i))//'.lead.dat')

         ! write corresponding particle spectrum:
         call WriteHistMP(hhMP(i),add=add0,mul=mul0,&
              &file='DoTIMES.MP'//trim(intToChar(140+i))//'.N.dat')
         call WriteHistMP(hhMPnu(i),add=add0,mul=mul0,&
              &file='DoTIMES.MPnu'//trim(intToChar(140+i))//'.N.dat')

         call WriteHistMP(hhMPlead(i),add=add0,mul=mul0,&
              &file='DoTIMES.MP'//trim(intToChar(140+i))//'.lead.N.dat')
         call WriteHistMP(hhMPnulead(i),add=add0,mul=mul0,&
              &file='DoTIMES.MPnu'//trim(intToChar(140+i))//'.lead.N.dat')
      end do

    end subroutine DoTimes_Write

    !**************************************************************************

    subroutine DoTimes_Finish()
      do i=1,3
         call RemoveHist(hh1(i))
      end do

      do i=0,7
         call RemoveHist(hh2(i))
      end do
    end subroutine DoTimes_Finish


    !**************************************************************************
    !****s* DoHiLeptonAnalysis/CalculateExclPiP
    ! NAME
    ! subroutine CalculateExclPiP
    ! PURPOSE
    ! Calculate the Cuts according Clasie et al.,arXiv:0701.1481
    ! for exclusive pi plus production (transparency).
    !
    ! INPUTS
    ! * integer :: iHist -- The histogram histMX(iHist) is used
    ! OUTPUT
    ! * The current event is added to histogram histMX(iHist).
    !
    ! NOTES
    ! we use the global variables:
    ! * nu,Q2,pPart,nRuns
    ! * histMX
    ! * PreEvArr
    ! The files fort.1173..1178 are written.
    !**************************************************************************
    subroutine CalculateExclPiP(iHist)
      use ParticleProperties, only: PartName
      use history, only: history_getgeneration
      use constants, only: mN

      integer, intent(in) :: iHist

      real, dimension (0:3) :: vecQ, vecMiss
      real :: MassCut, MassMiss, kTg, ph
      real, parameter :: dph = 0.1

      integer :: iEns, iNuc
      real :: weight

      type(tPreEvListEntry), POINTER :: pV
      integer :: j
      character*(15), dimension(8) :: AA


      if (pPart%ID /= pion) return
      if (pPart%charge.ne.1) return

      iEns = pPart%firstEvent / 1000
      iNuc = mod(pPart%firstEvent,1000)
      weight = pPart%perweight

      vecQ = (/nu, 0., 0., sqrt(nu**2+Q2) /)
      vecMiss = vecQ - pPart%momentum

      select case (MassNum)
      case (1)
         vecMiss(0) = vecMiss(0) + mN
         MassCut=99.9

      case (2)
         vecMiss(0) = vecMiss(0) + 2.01355*0.931494
         if (Q2.eq.1.1) then
            MassCut=2.00
         else if (Q2.eq.2.15) then
            MassCut=2.00
         else if (Q2.eq.3.00) then
            MassCut=2.025
         else if (Q2.eq.3.91) then
            MassCut=2.04
         else if (Q2.eq.4.69) then
            MassCut=2.08
         else
            write(*,*) 'Clasie: wrong Q2:',Q2
            return
         end if

      case (12)
         vecMiss(0) = vecMiss(0) + 12.0107*0.931494
         if (Q2.eq.1.1) then
            MassCut=11.35
         else if (Q2.eq.2.15) then
            MassCut=11.375
         else if (Q2.eq.3.00) then
            MassCut=11.400
         else if (Q2.eq.3.91) then
            MassCut=11.400
         else if (Q2.eq.4.69) then
            MassCut=11.425
         else
            write(*,*) 'Clasie: wrong Q2:',Q2
            return
         end if

      case (27)
         vecMiss(0) = vecMiss(0) + 26.9800*0.931494
         if (Q2.eq.1.1) then
            MassCut=25.275
         else if (Q2.eq.2.15) then
            MassCut=25.325
         else if (Q2.eq.3.00) then
            MassCut=25.35
         else if (Q2.eq.3.91) then
            MassCut=25.38
         else if (Q2.eq.4.69) then
            MassCut=25.40
         else
            write(*,*) 'Clasie: wrong Q2:',Q2
            return
         end if

      case (63)
         vecMiss(0) = vecMiss(0) + 63.5460*0.931494
         if (Q2.eq.1.1) then
            MassCut=59.35
         else if (Q2.eq.2.15) then
            MassCut=59.40
         else if (Q2.eq.3.00) then
            MassCut=59.40
         else if (Q2.eq.3.91) then
            MassCut=59.45
         else if (Q2.eq.4.69) then
            MassCut=59.50
         else
            write(*,*) 'Clasie: wrong Q2:',Q2
            return
         end if

      case (197)
         vecMiss(0) = vecMiss(0) + 196.9237*0.931494
         if (Q2.eq.1.1) then
            MassCut=183.57
         else if (Q2.eq.2.15) then
            MassCut=183.63
         else if (Q2.eq.3.00) then
            MassCut=183.63
         else if (Q2.eq.3.91) then
            MassCut=183.67
         else if (Q2.eq.4.69) then
            MassCut=183.74
         else
            write(*,*) 'Clasie: wrong Q2:',Q2
            return
         end if

      case default
         write(*,*) 'Clasie: wrong M:',MassNum
         return
      end select

      MassMiss = sqrt(vecMiss(0)**2-Dot_Product(vecMiss(1:3),vecMiss(1:3)))

      ! Select the forward pion production events:

      kTg = sqrt(Dot_Product(vecMiss(1:3),vecMiss(1:3)))
      ph = absmom(pPart)

      if (vecQ(3)-kTg.le.ph-dpH) return
      if (vecQ(3)-kTg.ge.ph+dpH) return


      if (MassMiss .ge. MassCut) then
         call AddHist(histMX(iHist),MassMiss,0.,weight)
         return
      end if
      call AddHist(histMX(iHist),MassMiss,weight,weight)

      pV => PreEvArr(iEns,iNuc)

      do j=1,8
         AA(j) = PartName(pV%preE(j)%ID,pV%preE(j)%charge,pV%preE(j)%antiparticle)
      end do
!      write(6,'(A,8A16)') 'Clasie..Primary Event:',AA
      write(1172+iHist,'(A,8A16)') 'Clasie..Primary Event:',AA

!      call WriteParticle(6,iEns,iNuc,pPart)
      call WriteParticle(1174+iHist,iEns,iNuc,pPart)

      write(1176+iHist,*) nRuns, pPart%firstEvent, &
           &history_getGeneration(pPart%history),pPart%history

    end subroutine CalculateExclPiP


    !**************************************************************************

!     subroutine DoingCollHistAnalysis
!       use CollHistory
!       use output
!
!       integer cc
!
!       call CollHist_WriteList(1./(1000*1000.*nLeptons)) ! 1000 as a dummy
!
!       do i=1,nEns
!          do j=1,size(pertPart,dim=2)
!             if(pertPart(i,j)%Id <  0) exit
!             if(pertPart(i,j)%Id <= 0) cycle
!
!             if(pertPart(i,j)%Id /= pion) cycle
!             if(pertPart(i,j)%charge .ne. 1) cycle
!
!             cc = CollHist_ClassifyHist(i,j)
! !            write(*,*) 'CCC:',cc
!             if (cc.eq.0) then
!                call WriteParticle(6,i,j,pertPart(i,j))
!                call CollHist_WriteHistParticle(6,i,j)
!             endif
!
!          end do
!       end do
!     end subroutine DoingCollHistAnalysis

    !**************************************************************************

    integer function SelectChannel(evtype)
      integer, intent(IN) :: evtype

      integer, parameter, dimension(17) :: iCH = (/&
           &  101, 102, 103, 104,&
           &  201, 202, 203, 204,&
           &  300, 400, 500,1000,&
           & 2001,2002,2003,2004,&
           & 5000 &
           &/)

      integer :: i

      if (evtype.gt.5000) then
         SelectChannel = 17
         return
      end if
      do i=1,17
         if (iCH(i).eq.evtype) then
            SelectChannel = i
            return
         end if
      end do


      SelectChannel = 0
      write(*,*) 'SelectChannel: strange one:',evtype
      call TRACEBACK()

    end function SelectChannel

    !**************************************************************************

    subroutine CalcFSIsqrts
      use twoBodyTools, only: sqrtS_Free
      use idTable, only: rho

      integer :: iEnsR, iPartR
      real :: srts

      if (pPart%Id /= rho) return
      if (pPart%charge.eq.0) return
      if (zH < 0.9) return

      do iEnsR=1,size(RealPart,dim=1)
         do iPartR=1,size(RealPart,dim=2)
            srts=sqrtS_Free( (/pPart,RealPart(iEnsR,iPartR)/) )
            call AddHist(hFSIsqrts, srts,pPart%perWeight)
         end do
      end do

    end subroutine CalcFSIsqrts


    !**************************************************************************

    subroutine MorrowAnalysis
      use idTable, only: rho

      logical,save :: DoMorrowInit = .true.

      type(particle), dimension(1,1:3)  :: finalState
      type(tPreEvListEntry), dimension(0:4),save :: PreCompArr
      integer :: ii, iiEntry, iEntry, iH1,iH2,iH, iProc
      logical :: flagOK
      real :: mass, MandelstamT
      type(tParticleListNode), POINTER :: pNode1,pNode2

      real :: nu,Q2,epsilon,Weight
      integer :: EventType
      real, dimension(0:3) :: pPhoton

      if (DoMorrowInit) then

         finalstate%ID = 0
         finalstate%antiparticle = .false.

         !=== p rho0 : entry 1
         finalstate(1,1:2)%ID     =(/1,rho/)
         finalstate(1,1:2)%charge =(/1,  0/)
         if (.not.CreateSortedPreEvent(finalState(1,:),PreCompArr(1)%preE)) call TRACEBACK()

         !=== Delta++ pi- : entry 2
         finalstate(1,1:2)%ID     =(/2,pion/)
         finalstate(1,1:2)%charge =(/2, -1/)
         if (.not.CreateSortedPreEvent(finalState(1,:),PreCompArr(2)%preE)) call TRACEBACK()

         !=== Delta0 pi+ : entry 3
         finalstate(1,1:2)%ID     =(/2,pion/)
         finalstate(1,1:2)%charge =(/0,  1/)
         if (.not.CreateSortedPreEvent(finalState(1,:),PreCompArr(3)%preE)) call TRACEBACK()

         !=== p pi+ pi- : entry 0 and 4
         finalstate(1,1:3)%ID     =(/1,pion,pion/)
         finalstate(1,1:3)%charge =(/1,  1, -1/)
         if (.not.CreateSortedPreEvent(finalState(1,:),PreCompArr(0)%preE)) call TRACEBACK()
         if (.not.CreateSortedPreEvent(finalState(1,:),PreCompArr(4)%preE)) call TRACEBACK()

         DoMorrowInit = .false.
      end if


      do i=1,nEns
         do j=1,MassNum

            !=== 1: Decide whether it is a two pion event:

            flagOK=.true.

            do ii=1,4
               if (PreEvArr(i,j)%preE(ii)%mass.ne.PreCompArr(0)%preE(ii)%mass) &
                    & flagOK=.false. ! attention: abuse of mass
            end do
            if (.not.flagOK) cycle ! This is not a p pi+ pi- event

            if (.not.EventInfo_HiLep_Get(i,j,Weight,nu,Q2,epsilon,EventType)) call TRACEBACK()
            pPhoton = (/nu,0.0,0.0,sqrt(nu**2+Q2)/)

            !=== 2a: Select origin process:
            ! iProc=1..5 = VMD, direct, anomalous, DIS, Misc

            select case (EventType)
            case (2001:2004)
               iProc = EventType-2000
            case default
               iProc = 5 ! this is the dummy value
            end select

            !=== 2b: Select origin channel:
            ! iEntry=1..5 = p rho0, D++ pi-, D0 pi+, PS, Misc

            iEntry = 5 ! this is the dummy value
            do iiEntry=1,4
               flagOK=.true.
               do ii=1,4
                  if (PreEvArr0(i,j)%preE(ii)%mass.ne.PreCompArr(iiEntry)%preE(ii)%mass) &
                       & flagOK=.false. ! attention: abuse of mass
               end do
               if (flagOK) then
                  iEntry = iiEntry
                  exit
               end if
            end do

            !=== 3: Calculate the 3 invariant masses and t:

            MandelstamT = 0.0

            pNode1 => EventArr(i,j)%particleList%first
            do
              if (.not. associated(pNode1)) exit
              iH1 = MorrowPartClass(pNode1%V%ID,pNode1%V%charge)

              pNode2=>pNode1%next
              do
                 if (.not. associated(pNode2)) exit
                 iH2 = MorrowPartClass(pNode2%V%ID,pNode2%V%charge)

                 mass = abs4(pNode1%V%momentum+pNode2%V%momentum)

                 select case (iH1+iH2)
                 case (6)
                    iH = 1 ! pi+ pi-
                    MandelstamT = abs4Sq(pPhoton-(pNode1%V%momentum+pNode2%V%momentum))
                 case (3)
                    iH = 2 ! p pi+
                 case (5)
                    iH = 3 ! p pi-
                 case default
                    write(*,*) 'wrong iH:',iH1,iH2
                    write(*,*) EventArr(i,j)%particleList%nEntries
                    call TRACEBACK()
                    ! This error means: it is not a p pi pi event,
                    ! but plus some additional particle
                 end select

                 call AddHistMC(hMorrow(iH),mass,iEntry,Weight)

                 pNode2 => pNode2%next
              end do ! pNode2

              pNode1 => pNode1%next
           end do ! pNode1

           call AddHistMC(hMorrowT(1),-MandelstamT,iEntry,Weight)
           if (iEntry==1) call AddHistMC(hMorrowT(2),-MandelstamT,iProc,Weight)

         end do ! j
      end do ! i

    end subroutine MorrowAnalysis

    integer function MorrowPartClass(iD,iQ)
      use idTable, only: nucleon, pion
      integer, intent(in) :: iD,iQ

      MorrowPartClass = 0
      select case (iD)
      case (nucleon)
         select case (iQ)
         case (1)
            MorrowPartClass = 1 ! == proton
         end select
      case (pion)
         select case (iQ)
         case (1)
            MorrowPartClass = 2 ! == pi+
         case (-1)
            MorrowPartClass = 4 ! == pi-
         end select
      end select

    end function MorrowPartClass

    !**************************************************************************

    subroutine BrooksAnalysis

      use constants, only: mN

      real :: zH, pT2, xB
      integer :: i,j,iCh
      integer :: iQ2, ixB, izH

      do i=1,nEns
         do j=1,size(pertPart,dim=2)
            if (pertPart(i,j)%Id <  0) exit
            if (pertPart(i,j)%Id /= pion) cycle   ! only pion
            if (pertPart(i,j)%charge .ne. 1) cycle ! only positive

            if (.not.(EventInfo_HiLep_Get(i,pertPart(i,j)%firstEvent,Weight,nu,Q2,epsilon,EventType))) then
               call TRACEBACK()
            end if

            zH = pertPart(i,j)%momentum(0)/nu
            xB = Q2/(2*mN*nu)
            pT2 = pertPart(i,j)%momentum(1)**2 + pertPart(i,j)%momentum(2)**2


            iQ2 = 0
            if (Q2 < 1.5) then
               iQ2 = 1
            else if (Q2 < 2.0) then
               iQ2 = 2
            else if (Q2 < 2.5) then
               iQ2 = 3
            else if (Q2 < 3.0) then
               iQ2 = 4
            else if (Q2 < 3.5) then
               iQ2 = 5
            else if (Q2 < 4.0) then
               iQ2 = 6
            end if

            izH = int(zH/0.1)+1

            ixB = 0
            if (xB < 0.19) then
               ixB = 1
            else if (xB < 0.28) then
               ixB = 2
            else if (xB < 0.37) then
               ixB = 3
            else if (xB < 0.46) then
               ixB = 4
            else if (xB < 0.55) then
               ixB = 5
            end if

            if (iQ2 > 0 .and. ixB > 0) then
               ArrBrooks(iQ2,ixB,izH,0) = ArrBrooks(iQ2,ixB,izH,0)+Weight
               ArrBrooks(iQ2,ixB,izH,1) = ArrBrooks(iQ2,ixB,izH,1)+(pT2*Weight)
            end if


            if (zH<0.5) cycle
            if (zH>0.7) cycle

            if (Q2<1.0) cycle
            if (Q2>3.0) cycle

            if (nu<2.0) cycle
            if (nu>4.0) cycle

            iCh = 0
            if (Q2>2) iCh = iCh + 4
            if (nu>3) iCh = iCh + 2
            if (zH>0.6) iCh = iCh + 1

            call AddHistMC(hBrooks,pT2,iCh+1,Weight)


         end do
      end do
    end subroutine BrooksAnalysis

    !**************************************************************************

    subroutine FindRho0Analysis
      use idTable, only: rho
      use ParticleProperties, only: hadron
      use constants, only: mN

      type(tParticleList) :: PartsOut
      real, dimension(:), allocatable :: Probs
      logical :: AccFlag,CutFlag
      real :: MandelstamT

      do i=1,nEns
         do j=1,MassNum
            if (.not.EventInfo_HiLep_Get(i,j,Weight,nu,Q2,epsilon,EventType,W=Wcalc)) cycle

!            write(*,'("#### EVENT ",i5,i3," ##### ",i7)') i,j,EventType
!            call ParticleList_Print(EventArr(i,j)%particleList)
!!            stop

            ! In the case of a decaying rho0, we first have to recombine them from the
            ! pions of the decay. If the rho0 was selected not to decay, we directly use
            ! the rho0 from the event for the analysis.

            if (hadron(rho)%stability .ne. 0) then

               allocate(Probs(EventArr(i,j)%particleList%nEntries))

               pNode => EventArr(i,j)%particleList%first
               Probs = 0
               ii = 1
               do while (associated(pNode))
                  call checkCuts(pNode%V,nu,Q2,Ebeam,phi_Lepton,AccFlag,AccWeight)

                  if (iExperiment.eq.18) then
                     if (pNode%V%momentum(0).lt.10.0) AccWeight=0.0
                  end if

                  Probs(ii) = AccWeight
                  pNode=>pNode%next
                  ii=ii+1
               end do

               call FindRho0(EventArr(i,j)%particleList,PartsOut,&
                    & (/mN+nu,0.0,0.0,sqrt(nu**2+Q2)/),Probs)

!               call ParticleList_Print(PartsOut)

               pNode => PartsOut%first

            else
               pNode => EventArr(i,j)%particleList%first
            end if

            do while (associated(pNode))

               ! first check, whether it is really a rho0:
               AccFlag = .true.
               if (pNode%V%charge .ne. 0) AccFlag = .false.
               if (pNode%V%ID /= rho) AccFlag = .false.
               if (.not.AccFlag) then
                  pNode=>pNode%next
                  cycle
               end if

               CutFlag = .true.

               select case (iExperiment)
               case (14,15) !  JLAB,  4/5GeV (CLAS, rho0)
!                  if (abs(pNode%V%position(2)).gt.0.1) AccFlag = .false. ! <-- Delta E; not used acc. Kawtar

                  MandelstamT = abs4Sq((/nu,0.0,0.0,sqrt(nu**2+Q2)/)-pNode%V%momentum)
                  if (MandelstamT.lt.-0.4) CutFlag = .false.
                  if (MandelstamT.gt.-0.1) CutFlag = .false.

                  if (pNode%V%momentum(0).lt.0.9*nu) CutFlag = .false.

               case default

                  ! for TEMORARY USE !!!!!!!!!!!!!!
                  MandelstamT = abs4Sq((/nu,0.0,0.0,sqrt(nu**2+Q2)/)-pNode%V%momentum)
!                     if (MandelstamT.lt.-0.5) CutFlag = .false.
!                     if (MandelstamT.lt.-0.4) CutFlag = .false.
!                     if (MandelstamT.gt.-0.1) CutFlag = .false.

                  if (pNode%V%momentum(0).lt.0.9*nu) CutFlag = .false.
!                     if (pNode%V%momentum(0).lt.0.8*nu) CutFlag = .false.

               end select

               ! Default behaviour:
!               AccFlag = CutFlag

               if (AccFlag) then

!!$                      write(*,*) '%%%%%%%%%%%%%%%%%%%% found rho0:'
!!$                      call WriteParticle(6,1,1,pNode%V)
!!$                      write(*,*)
!!$                      write(*,*) Probs(1:ii-1)
!!$                      call ParticleList_Print(EventArr(i,j)%particleList)

                  w0 = pNode%V%perWeight
                  w1 = w0
                  if (pNode%V%scaleCS.lt.1.0) w1=0 ! do not count 'half' rhos
                                                   ! abuse: scaleCS==prob.decay
                  if (AccFlag.neqv.CutFlag) w1=0


                  call AddHist(hRho0MV,pNode%V%mass,w0,w1)
                  call AddHist(hRho0MX,pNode%V%position(1),w0,w1)      ! abuse: pos(1)==MX !!!
                  call AddHist(hRho0DE,pNode%V%position(2),w0,w1)      ! abuse: pos(2)==DE !!!
                  call AddHist(hRho0DecTime,pNode%V%position(3),w0,w1) ! abuse: pos(3)==DecayTime !!!
                  call AddHist(hRho0Theta,pNode%V%offShellParameter,w0,w1) ! abuse: ...==thetaDecay !!!

                  call AddHist(hRho0Mom,pNode%V%momentum(0),w0,w1)


                  hWeight = PreEvArr0(i,j)%weight
                  PreEvArr0(i,j)%weight = w0 * 1000
                  call PreEvList_INSERT(PreEvListRhoOrig,PreEvArr0(i,j))
                  PreEvArr0(i,j)%weight = hWeight

                  W = sqrt(mN**2-Q2+2*mN*nu) ! != Wcalc!
                  zH = pNode%V%momentum(0)/nu

                  select case (EventType)
                  case (2001:2004)
                     iProc = EventType-2000
                  case default
                     iProc = 0
                  end select

                  call AddHist2D(h2D_rho0NuQ2,h2D_rho0NuQ2_proc(iProc), (/nu,Q2/), w0,w1)
                  call AddHist2D(h2D_rho0WQ2,h2D_rho0WQ2_proc(iProc),  (/W,Q2/),  w0,w1)
                  call AddHist(hRho0zH,hRho0zH_proc(iProc),     zH,          w0,w1)
                  call AddHist(hRho0t,hRho0t_proc(iProc),     -MandelstamT, w0,w1)
                  call AddHist(hRho0sigmaW,hRho0sigmaW_proc(iProc), Wcalc, w1,w1*w1)

                  call AddHist(hRho0W,W,w0,w1)
                  call AddHist(hRho0Q2,Q2,w0,w1)

               end if

               pNode=>pNode%next
            end do


            if (hadron(rho)%stability /= 0) then
               call ParticleList_CLEAR(PartsOut,.true.)
               deallocate(Probs)
            end if

         end do
      end do

    end subroutine FindRho0Analysis

    !**************************************************************************

    subroutine MandelTAnalysis()

      real, dimension(0:3) :: pPhoton
      real :: MandelstamT
      integer :: i,j
      real :: w2

      do i=1,nEns
         do j=1,MassNum
            if (.not.EventInfo_HiLep_Get(i,j,Weight,nu,Q2,epsilon,EventType,W=Wcalc)) cycle

            pPhoton = (/nu,0.0,0.0,sqrt(nu**2+Q2)/)
            pNode => EventArr(i,j)%particleList%first

            w2 = 0.0
            if (EventArr(i,j)%particleList%nEntries .le. 2) w2 = Weight ! 'exclusive' event

            do while (associated(pNode))
               pPart => pNode%V

               MandelstamT = abs4Sq(pPhoton-pPart%momentum)
               call AddHistMP(HMP_MandelstamT, pPart, -MandelstamT,  Weight,w2)

               pNode=>pNode%next
            end do
         end do
      end do

    end subroutine MandelTAnalysis

    !**************************************************************************
    subroutine CentralNAnalysis()
      use vector, only: absVec
      !use dichteDefinition
      !use densitymodule
      use constants, only: mN

      integer :: i,j,iNr
      real, dimension(3) :: pos
      real :: b,bT,bZ
      real :: p,pT, Ekin
      real :: z
!      real,parameter :: densMax = 0.16
!       type(dichte) :: dens
      logical :: isUnbound

      do i=1,nEns
         do j=1,MassNum
            iNr = j
            if (HiLepton_getRealRun()) iNr = -1
            if (.not.EventInfo_HiLep_Get(i,iNr,Weight,nu,Q2,epsilon,EventType,W=Wcalc,pos=pos)) cycle

            b = absVec(pos(1:3))
            bT = absVec(pos(1:2))
            bZ = pos(3)

            call AddHist(H_CN_b, b, Weight)
            call AddHist(H_CN_bT, bT, Weight)
            call AddHist(H_CN_bZ, bZ, Weight)

            pNode => EventArr(i,j)%particleList%first
            do while (associated(pNode))
               pPart => pNode%V

               if ((pPart%ID .eq. 1).and.(.not.pPart%antiparticle)) then

!                  dens = densityAt(pPart%position)
!                  isUnbound = (dens%baryon(0) .lt. densMax/100)

                  isUnbound = .true. ! accept all

                  if (isUnbound) then

                     p = absVec(pPart%momentum(1:3))
                     pT = absVec(pPart%momentum(1:2))
                     z = pPart%momentum(0)/nu
                     Ekin = pPart%momentum(0)-mN

                     call AddHist2D(H2D_CN_b(pPart%charge), (/b,p/), Weight)
                     call AddHist2D(H2D_CN_bT(pPart%charge), (/bT,p/), Weight)
                     call AddHist2D(H2D_CN_bZ(pPart%charge), (/bZ,p/), Weight)
                     call AddHist2D(H2D_CN_bZE(pPart%charge), (/bZ,Ekin/), Weight)

                     call AddHist2D(H2D_CN_pTz(pPart%charge), (/pT,z/), Weight)

                  end if
               end if

               pNode=>pNode%next
            end do
            if (HiLepton_getRealRun()) exit
         end do
      end do

    end subroutine CentralNAnalysis
  end subroutine DoHiLeptonAnalysis


  !****************************************************************************
  !****s* HiLeptonAnalysis/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "HiLepton_Analysis"
  !****************************************************************************
  subroutine readInput

    use output

    integer :: ios

    !**************************************************************************
    !****n* HiLeptonAnalysis/HiLepton_Analysis
    ! NAME
    ! NAMELIST HiLepton_Analysis
    ! PURPOSE
    ! Includes the switches:
    ! * DoTimes
    ! * DoOutChannels
    ! * DoInvMasses
    ! * DoFindRho0
    ! * DoClasie
    ! * DoMorrow
    ! * DoBrooks
    ! * DoMandelT
    ! * DoClassifyFirst
    ! * DoFSIsqrts
    ! * DoCentralN
    ! * DoLeptonKinematics
    ! * DoHadronKinematics
    !**************************************************************************
    NAMELIST /HiLepton_Analysis/ DoTimes, DoOutChannels, DoInvMasses, &
         & DoFindRho0, DoClasie, DoMorrow,DoBrooks,DoClassifyFirst, &
         & DoFSIsqrts, DoMandelT, DoCentralN, &
         & DoLeptonKinematics, DoHadronKinematics

    call Write_ReadingInput('HiLepton_Analysis',0)
    rewind(5)
    read(5,nml=HiLepton_Analysis,IOSTAT=ios)
    call Write_ReadingInput('HiLepton_Analysis',0,ios)

    write(*,*) 'DoLeptonKinematics :',DoLeptonKinematics
    write(*,*) 'DoHadronKinematics :',DoHadronKinematics
    write(*,*)
    write(*,*) 'DoTimes         :',DoTimes
    write(*,*) 'DoOutChannels   :',DoOutChannels
    write(*,*) 'DoInvMasses     :',DoInvMasses
    write(*,*) 'DoFindRho0      :',DoFindRho0
    write(*,*) 'DoClasie        :',DoClasie
    write(*,*) 'DoMorrow        :',DoMorrow
    write(*,*) 'DoBrooks        :',DoBrooks
    write(*,*) 'DoFSIsqrts      :',DoFSIsqrts
    write(*,*) 'DoMandelT       :',DoMandelT
    write(*,*) 'DoCentralN      :',DoCentralN
    write(*,*)
    write(*,*) 'DoClassifyFirst :',DoClassifyFirst
    write(*,*)

    DoEventAdd  = DoOutChannels.or.DoInvMasses.or.DoFindRho0.or.DoClassifyFirst.or.DoMorrow.or.DoClasie.or.DoMandelT.or.DoCentralN
    DoEventAdd2 = DoOutChannels.or.DoMorrow

    write(*,*) 'DoEventAdd    =',DoEventAdd
    write(*,*) 'DoEventAdd2   =',DoEventAdd2

    call Write_ReadingInput('HiLepton_Analysis',1)
    FlagReadInput = .false.

  end subroutine readInput




  !****************************************************************************
  !****s* HiLeptonAnalysis/resetBinning
  ! NAME
  ! subroutine resetBinning
  ! PURPOSE
  ! Reset all Histogramms used.
  !
  ! Allocates all arrays and sets necessary parameters.
  ! INPUTS
  ! ---
  ! OUTPUT
  ! ---
  !****************************************************************************
  subroutine resetBinning
    use output, only: IntToChar
    use initHiLepton, only: GetPhotonKin,GetEnergies
    use constants, only: mN

    integer :: i
    integer :: iNu,iQ2,iZH, nNu,nQ2,nZH
    real :: nu,Q2,W
    real :: Ebeam_,EIC_Ee_,EIC_EA_

    type(ValueRange) :: hh_Q2R, hh_nuR, hh_WR, hh_zHR

    character(5), dimension(2),    parameter :: AccName = (/'noAcc','inAcc'/)

    character(13), dimension(3), parameter :: BinNameZ = (/'[zh=0.2..0.4]','[zh=0.4..0.7]','[zh=0.7..1.2]'/)
    character(11), dimension(3), parameter :: BinNameN = (/'[nu= 6..12]',  '[nu=12..17]',  '[nu=17..24]'/)
    character( 9), dimension(2), parameter :: BinNameP = (/'[pT2<0.7]',    '[pT2>0.7]'/)

    pTR = ValueRange(0.0,5.0, 0.05)
    zHR = ValueRange(0.0,1.0, 0.02)
    Q2R = ValueRange(0.0,30.0, 1.0)

    select case (iExperiment)

    case (0) ! no experiment/fixed kinematics
       iDetector = 90

       call GetPhotonKin(nu,Q2,W)
       nuR = ValueRange( nu-0.1, nu+0.1, 0.1)
       Q2R = ValueRange( Q2-0.1, Q2+0.1, 0.1)


    case (1) ! HERMES: D,N.Kr
       !-------------------------------------------------------
       iDetector = 1
       nuR = ValueRange(7.0, 25.0, 1.0)
       Q2R = ValueRange(0.0, 18.0, 0.5)
       zHR%l=0.2     !minimum zh of hadron for binning (except. zh-binning)
       EHR%l=1.4

    case (2) ! HERMES: Ne
       !-------------------------------------------------------
       iDetector = 1
       nuR = ValueRange(2.0, 24.0, 1.0)
       Q2R = ValueRange(0.0, 18.0, 0.5)
       zHR%l=0.2     !minimum zh of hadron for binning (except. zh-binning)
       EHR%l=1.4

    case (11) ! HERMES: final paper
       !-------------------------------------------------------
       iDetector = 1
       nuR = ValueRange(6.0, 24.0, 1.0)
       Q2R = ValueRange(0.0, 18.0, 0.5)
       zHR%l=0.2     !minimum zh of hadron for binning (except. zh-binning)
       zHR%u=1.2
       EHR%l=1.4

    case (13) ! HERMES: pT-broadening
       !-------------------------------------------------------
       iDetector = 1
       nuR = ValueRange(5.0, 24.0, 1.0)
       Q2R = ValueRange(0.0, 18.0, 0.5)
       zHR%l=0.2     !minimum zh of hadron for binning (except. zh-binning)
       zHR%u=1.2
       EHR%l=1.4

    case (10) ! HERMES@12GeV
       !-------------------------------------------------------
       iDetector = 1
       nuR = ValueRange(2.0, 11.0, 0.2)
       Q2R = ValueRange(0.0,  6.0, 0.2)
       zHR%l=0.2     !minimum zh of hadron for binning (except. zh-binning)
       EHR%l=1.4

    case (4,19,20) ! JLAB@12GeV (CLAS)
       !-------------------------------------------------------
       iDetector = 3
       nuR = ValueRange(2.0, 11.0, 0.2)
       Q2R = ValueRange(0.0, 16.0, 0.2)
       zHR%l=0.2     !minimum zh of hadron for binning (except. zh-binning)
       EHR%l=1.4

    case (5) ! JLAB@5GeV (CLAS)
       !-------------------------------------------------------
       iDetector = 4
       nuR = ValueRange(2.0, 4.6, 0.1)
       Q2R = ValueRange(0.0, 5.0, 0.1)
       zHR%l=0.2     !minimum zh of hadron for binning (except. zh-binning)
!       EHR%l=1.4
       EHR%l=0.0

    case (14) ! JLAB@5GeV (CLAS, rho0)
       !-------------------------------------------------------
       iDetector = 4
       nuR = ValueRange(2.0, 4.6, 0.1)
       Q2R = ValueRange(0.0, 5.0, 0.1)
       zHR%l=0.2     !minimum zh of hadron for binning (except. zh-binning)
!       EHR%l=1.4
       EHR%l=0.0

    case (15) ! JLAB@4GeV (CLAS, rho0)
       !-------------------------------------------------------
       iDetector = 4
       nuR = ValueRange(2.0, 4.6, 0.1)
       Q2R = ValueRange(0.0, 5.0, 0.1)
       zHR%l=0.2     !minimum zh of hadron for binning (except. zh-binning)
!       EHR%l=1.4
       EHR%l=0.0

    case (6) ! EMC: 100
       !-------------------------------------------------------
       iDetector = 2
       nuR = ValueRange(10.0,  85.0, 2.0)
       Q2R = ValueRange( 0.0, 160.0, 2.0)
       zHR%l=0.0
       EHR%l=3.0

    case (8) ! EMC: 200
       !-------------------------------------------------------
       iDetector = 2
       nuR = ValueRange(30.0, 170.0,  5.0)
       Q2R = ValueRange( 0.0, 320.0,  2.0)
       zHR%l=0.0
       EHR%l=3.0

    case (9) ! EMC: 280
       !-------------------------------------------------------
       iDetector = 2
       nuR = ValueRange(50.0, 240.0,  5.0)
       Q2R = ValueRange( 0.0, 450.0,  2.0)
       zHR%l=0.0
       EHR%l=5.0

    case (12) ! Mainz, Yoon: 1.5
       !-------------------------------------------------------
       iDetector = 90
       nuR = ValueRange( 0.0, 2.0,  0.02)
       Q2R = ValueRange( 0.0, 20.0, 0.02)

    case (16) ! EIC
       !-------------------------------------------------------
       iDetector = 90

       call GetEnergies(Ebeam_,EIC_Ee_,EIC_EA_)

       nuR = ValueRange( 0.0, Ebeam_,  Ebeam_/100)
       zHR%l=0.4
       Q2R = ValueRange(0.0,30.0, 0.2)

    case (17) ! total cross section
       !-------------------------------------------------------
       iDetector = 90

       call GetEnergies(Ebeam_)

       nuR = ValueRange( 0.0, Ebeam_,  Ebeam_/100)
       Q2=min(10.0,2*mN*Ebeam_+mN**2)
       Q2R = ValueRange(0.0,Q2, Q2/100)

    case (18) ! E665
       !-------------------------------------------------------
       iDetector = 90

       call GetEnergies(Ebeam_)

       nuR = ValueRange( 0.0, Ebeam_,  Ebeam_/100)
       Q2=min(10.0,2*mN*Ebeam_+mN**2)
       Q2R = ValueRange(0.0,Q2, Q2/100)

    case default
       !-------------------------------------------------------
       write(*,*) 'Binning not prepared for iExp=',iExperiment
       stop
       !-------------------------------------------------------
    end select

    hh_Q2R = Q2R
    hh_nuR = nuR
    hh_WR  = ValueRange( 1.5, 3.2,  0.02)
    hh_zHR = ValueRange( 0.0, 1.2,  0.005)

    nLeptons=0.

    call CreateHist(hLep_Q2, "N_e(Q2), <nu>(Q2)", Q2R%l,Q2R%u,Q2R%d)
    call CreateHist(hLep_nu, "N_e(nu), <Q2>(nu)", 0.0  ,nuR%u,nuR%d)

    do i=1,2
       call CreateHist(hKin_zH_Q2(i), "<Q2>(zH) ("//trim(AccName(i))//")", zHR%l,zHR%u,zHR%d)
       call CreateHist(hKin_nu_Q2(i), "<Q2>(nu) ("//trim(AccName(i))//")", nuR%l,nuR%u,nuR%d)
       call CreateHist(hKin_zH_nu(i), "<nu>(zH) ("//trim(AccName(i))//")", zHR%l,zHR%u,zHR%d)
       call CreateHist(hKin_nu_zH(i), "<zH>(nu) ("//trim(AccName(i))//")", nuR%l,nuR%u,nuR%d)
    end do


    call CreateHist(hCH_zH, "dN_ch/(N_e dzH)",  0.0,  1.2,  zhR%d)
    call CreateHist(hCH_nu, "dN_ch/(N_e dnu)",  0.0,  nuR%u,nuR%d)
    call CreateHist(hCH_Q2, "dN_ch/(N_e dQ2)",  Q2R%l,Q2R%u,Q2R%d)
    call CreateHist(hCH_pT2,"dN_ch/(N_e dpT2)", pTR%l,pTR%u,pTR%d)

    call CreateHistMP(hMP_zH, "dN_id/(N_e dzH)",  0.0,  1.2,  zhR%d, iPartSet)
    call CreateHistMP(hMP_nu, "dN_id/(N_e dnu)",  0.0,  nuR%u,nuR%d, iPartSet)
    call CreateHistMP(hMP_Q2, "dN_id/(N_e dQ2)",  Q2R%l,Q2R%u,Q2R%d, iPartSet)
    call CreateHistMP(hMP_pT2,"dN_id/(N_e dpT2)", pTR%l,pTR%u,pTR%d, iPartSet)

    call CreateHistMP(hMP_pT2zH, "<pT2>_id(zH)",  0.0,  1.2,  zhR%d, iPartSet)
    call CreateHistMP(hMP_pT2nu, "<pT2>_id(nu)",  0.0,  nuR%u,nuR%d, iPartSet)
    call CreateHistMP(hMP_pT2Q2, "<pT2>_id(Q2)",  Q2R%l,Q2R%u,Q2R%d, iPartSet)
    call CreateHistMP(hMP_pT2pT2,"<pT2>_id(pT2)", pTR%l,pTR%u,pTR%d, iPartSet)

    do i=-1,1
       call CreateHist2D(H2DpTAvePion(i),"<pT2>(nu,zH) "//piName(i), &
            &(/nuR%l,zHR%l/), (/nuR%u,zhR%u/), (/nuR%d,zhR%d/) , .true.)

       call CreateHist2D(H2DpTPionZH(i),"pT2 vs zH "//piName(i), &
            &(/0.0,pTR%l/), (/1.1,pTR%u/), (/0.02,pTR%d/) , .true.)
       call CreateHist2D(H2DpTPionNU(i),"pT2 vs nu "//piName(i), &
            &(/0.0,pTR%l/), (/nuR%u,pTR%u/), (/nuR%d,pTR%d/) , .true.)
       call CreateHist2D(H2DpTPionQ2(i),"pT2 vs Q2 "//piName(i), &
            &(/Q2R%l,pTR%l/), (/Q2R%u,pTR%u/), (/Q2R%d,pTR%d/) , .true.)

    end do

    call CreateHist2D(H2DleptonXS(0), "XS (total)", &
         &(/Q2R%l,nuR%l/), (/Q2R%u,nuR%u/), (/Q2R%d,nuR%d/), .true.)
    call CreateHist2D(H2DleptonXS(1), "XS (low energy)", &
         &(/Q2R%l,nuR%l/), (/Q2R%u,nuR%u/), (/Q2R%d,nuR%d/), .true.)
    call CreateHist2D(H2DleptonXS(2), "XS (FRITIOF)", &
         &(/Q2R%l,nuR%l/), (/Q2R%u,nuR%u/), (/Q2R%d,nuR%d/), .true.)
    call CreateHist2D(H2DleptonXS(3), "XS (PYTHIA)", &
         &(/Q2R%l,nuR%l/), (/Q2R%u,nuR%u/), (/Q2R%d,nuR%d/), .true.)

    call CreateHist2D(H2DleptonXS(4), "XS (PYTHIA VMD)", &
         &(/Q2R%l,nuR%l/), (/Q2R%u,nuR%u/), (/Q2R%d,nuR%d/), .true.)
    call CreateHist2D(H2DleptonXS(5), "XS (PYTHIA direct)", &
         &(/Q2R%l,nuR%l/), (/Q2R%u,nuR%u/), (/Q2R%d,nuR%d/), .true.)
    call CreateHist2D(H2DleptonXS(6), "XS (PYTHIA GVMD)", &
         &(/Q2R%l,nuR%l/), (/Q2R%u,nuR%u/), (/Q2R%d,nuR%d/), .true.)
    call CreateHist2D(H2DleptonXS(7), "XS (PYTHIA DIS)", &
         &(/Q2R%l,nuR%l/), (/Q2R%u,nuR%u/), (/Q2R%d,nuR%d/), .true.)

    if (DoClassifyFirst) then
       do i=0,8
          call CreateHistMP(hMP_zHClass(i), "dN_id/(N_e dzH) "//trim(IntToChar(i)),  0.0,  1.2,  0.01, iPartSet)
          call CreateHistMP(hMP_pT2Class(i), "d<pT2>/dzH "//trim(IntToChar(i)),  0.0,  1.2,  0.01, iPartSet)
       end do

       call CreateHistMP(hMP_zH_Generation(-1), "dN_id/(N_e dzH) Generation ALL",  0.0,  1.2,  0.01, iPartSet)
       do i=0,10
          call CreateHistMP(hMP_zH_Generation(i), "dN_id/(N_e dzH) Generation "//intToChar(i),  0.0,  1.2,  0.01, iPartSet)
       end do
    end if

    if (CollHist_GetDoCollHistory()) then
       call CreateHist2D(H2D_CollHistPT(1), "d<pT2>(chgd.pi)/dzH (classified)", (/-16.5,0.0/), (/16.5,1.2/), (/1.0,0.01/))
       call CreateHist2D(H2D_CollHistPT(2), "d<pT2>(chgd.pi)/dzH (classified,100)", (/-16.5,0.0/), (/16.5,1.2/), (/1.0,0.01/))
    end if

    select case (iExperiment)
    case (11)  ! HERMES (final paper)
       do i=1,3
          call CreateHistMP(hMP_nu_zh(i), "dN_id/(N_e dnu)"//BinNameZ(i),  0.0,  nuR%u,nuR%d, 1)
          call CreateHistMP(hMP_Q2_zh(i), "dN_id/(N_e dQ2)"//BinNameZ(i),  Q2R%l,Q2R%u,Q2R%d, 1)
          call CreateHistMP(hMP_pT2_zh(i),"dN_id/(N_e dpT2)"//BinNameZ(i), pTR%l,pTR%u,pTR%d, 1)

          call CreateHistMP(hMP_zH_nu(i), "dN_id/(N_e dzH)"//BinNameN(i),  0.0,  1.2,  zhR%d, 1)
          call CreateHistMP(hMP_Q2_nu(i), "dN_id/(N_e dQ2)"//BinNameN(i),  Q2R%l,Q2R%u,Q2R%d, 1)
          call CreateHistMP(hMP_pT2_nu(i),"dN_id/(N_e dpT2)"//BinNameN(i), pTR%l,pTR%u,pTR%d, 1)
       end do
       do i=1,2
          call CreateHistMP(hMP_nu_pT(i), "dN_id/(N_e dnu)"//BinNameP(i), 0.0,  nuR%u,nuR%d, 1)
          call CreateHistMP(hMP_zH_pT(i), "dN_id/(N_e dzH)"//BinNameP(i), 0.0,  1.2,  zhR%d, 1)
          call CreateHistMP(hMP_Q2_pT(i), "dN_id/(N_e dQ2)"//BinNameP(i), Q2R%l,Q2R%u,Q2R%d, 1)
       end do

    case (5) ! JLAB,  5GeV (CLAS)

       if (allocated(hJLAB5_R_nu)) then
          write(*,*) 'hJLAB5_R_nu allready allocated. stop!'
          stop
       end if

       call BinCutsJLAB5_R(-99.0,99.0,99.0, nNu,nQ2,nZH)

       if (max(nNu,nQ2,nZH).gt.9) then
          write(*,*) 'in resetBinning: BinCutsJLAB5_R > 9. stop!'
          stop
       end if

       allocate(hJLAB5_R_nu(nQ2,nZH))
       do iQ2=1,nQ2
          do iZH=1,nZH
             call CreateHist(hJLAB5_R_nu(iQ2,iZH), &
                  & "dN_pi/(N_e dnu)[bin(Q2)="//Achar(iQ2+48)//",bin(zH)="//Achar(iZH+48)//"]", &
                  & 0.0,  nuR%u,nuR%d)
          end do
       end do

       allocate(hJLAB5_R_Q2(nNu,nZH))
       do iNu=1,nNu
          do iZH=1,nZH
             call CreateHist(hJLAB5_R_Q2(iNu,iZH), &
                  & "dN_pi/(N_e dQ2)[bin(nu)="//Achar(iNu+48)//",bin(zH)="//Achar(iZH+48)//"]", &
                  & Q2R%l,Q2R%u,Q2R%d)
          end do
       end do

       allocate(hJLAB5_R_ZH(nQ2,nNu))
       do iQ2=1,nQ2
          do iNu=1,nNu
             call CreateHist(hJLAB5_R_zH(iQ2,iNu), &
                  & "dN_pi/(N_e dzH)[bin(Q2)="//Achar(iQ2+48)//",bin(nu)="//Achar(iNu+48)//"]", &
                  &  0.0,  1.2,  zhR%d)
          end do
       end do

       call BinCutsJLAB5_R(-99.0,99.0,99.0, nNu,nQ2,nZH)

       allocate(hJLAB5_pT2_ARR(nNu,nQ2,nZH,0:2))
       hJLAB5_pT2_ARR = 0.0

    case (14,15) !  JLAB,  4/5GeV (CLAS, rho0)

       hh_Q2R = ValueRange( 0.6,  4.0,  0.05) ! we use here a different binning
       hh_nuR = ValueRange( 2.05, 4.85, 0.05)
       hh_WR  = ValueRange( 1.8,  3.2,  0.02)

    case (16) ! EIC

       call CreateHistMP(hMP_EIC_zH, "dN/dzH", 0.0,  1.2,  zhR%d, 3)
       call CreateHistMP(hMP_EIC_nu, "dN/dnu", nuR%l,nuR%u,nuR%d, 3)
       call CreateHistMP(hMP_EIC_Q2, "dN/dnu", Q2R%l,Q2R%u,Q2R%d, 3)
       do iQ2=1,5
          call CreateHistMP(hMP_EIC_nu_Q2(iQ2), "dN/dnu [bin(Q2)="//Achar(iQ2+48)//"]", nuR%l,nuR%u,nuR%d, 3)
          call CreateHistMP(hMP_EIC_zH_Q2(iQ2), "dN/dzH [bin(Q2)="//Achar(iQ2+48)//"]", 0.0,  1.2,  zhR%d, 3)
       end do

    case (18) ! E665
       hh_zHR = ValueRange(0.95,1.05, 0.001)
       hh_WR  = ValueRange( 1.8,20.0,   0.1)

    end select

    if (DoFindRho0) then
       call CreateHist(hRho0MV, "rho0, M_V", 0.0,2.0,0.01)
       call CreateHist(hRho0MX, "rho0, M_X", 0.0,5.0,0.05)
       call CreateHist(hRho0DE, "rho0, Delta E", -5.0,20.0,0.1)
       call CreateHist(hRho0DecTime, "rho0, Decay Time", 0.0,20.0,0.1)
       call CreateHist(hRho0Mom, "rho0, P", 0.0,10.0,0.02)
       call CreateHist(hRho0Theta, "rho0, theta", 0.0,180.0,2.0)

       call CreateHist(hRho0W, "rho0, W",  hh_WR%l,hh_WR%u,hh_WR%d)
       call CreateHist(hRho0Q2, "rho0, Q2",  hh_Q2R%l,hh_Q2R%u,hh_Q2R%d)

       call CreateHist2D(h2D_rho0NuQ2, "dN_rho0/N_e dnu dQ2 (with cuts)", &
            & (/hh_nuR%l,hh_Q2R%l/), (/hh_nuR%u,hh_Q2R%u/),(/hh_nuR%d,hh_Q2R%d/),.true.)
       call CreateHist2D(h2D_rho0WQ2, "dN_rho0/N_e dW dQ2 (with cuts)", &
            & (/hh_WR%l,hh_Q2R%l/), (/hh_WR%u,hh_Q2R%u/),(/hh_WR%d,hh_Q2R%d/),.true.)
       call CreateHist(hRho0zH, "dN_rho0/N_e dzH (with cuts)", hh_zHR%l,hh_zHR%u,hh_zHR%d)

       call CreateHist(hRho0t, "dN_rho0/N_e d|t| (with cuts)", 0.0,  1.0,  0.005)
       call CreateHist(hRho0sigmaW, "dsigma_rho0/N_e dW (with cuts)", hh_WR%l, hh_WR%u,  hh_WR%d)

       do i=0,4
          call CreateHist2D(h2D_rho0NuQ2_proc(i), "dN_rho0/N_e dnu dQ2 proc="//Achar(i+48)//" (with cuts)", &
             & (/hh_nuR%l,hh_Q2R%l/), (/hh_nuR%u,hh_Q2R%u/),(/hh_nuR%d,hh_Q2R%d/),.true.)
          call CreateHist2D(h2D_rho0WQ2_proc(i), "dN_rho0/N_e dW dQ2 proc="//Achar(i+48)//" (with cuts)", &
               & (/hh_WR%l,hh_Q2R%l/), (/hh_WR%u,hh_Q2R%u/),(/hh_WR%d,hh_Q2R%d/),.true.)
          call CreateHist(hRho0zH_proc(i), "dN_rho0/N_e dzH proc="//Achar(i+48)//" (with cuts)", &
               & hh_zHR%l,hh_zHR%u,hh_zHR%d)
          call CreateHist(hRho0t_proc(i), "dN_rho0/N_e d|t| proc="//Achar(i+48)//" (with cuts)", &
               & 0.0,  1.0,  0.005)
          call CreateHist(hRho0sigmaW_proc(i), "dsigma_rho0/N_e dW proc="//Achar(i+48)//" (with cuts)",&
               & hh_WR%l, hh_WR%u,  hh_WR%d)
       end do
    end if

  end subroutine resetBinning

  !****************************************************************************
  !****s* HiLeptonAnalysis/DoBinningLepton
  ! NAME
  ! subroutine DoBinningLepton(nu,Q2,Weight)
  ! PURPOSE
  ! By calling this routine, the histogramms according the leptons are filled.
  ! INPUTS
  ! * real :: nu -- nu value of Event
  ! * real :: Q2 -- Q**2 value of Event
  ! * real :: Weight -- perturbative weight of Event
  ! OUTPUT
  ! ---
  !****************************************************************************
  subroutine DoBinningLepton(nu,Q2,Weight)

    real,intent(in) :: nu,Q2,Weight

!    write(*,*) 'DoBinningLepton: ',nu,Q2,Weight

    if (nu.ge.nuR%l.and.nu.le.nuR%u) then
       nLeptons=nLeptons+Weight
    end if

    call AddHist(hLep_Q2, Q2, Weight,Weight*nu)
    call AddHist(hLep_nu, nu, Weight,Weight*Q2)

  end subroutine DoBinningLepton


  !****************************************************************************
  !****s* HiLeptonAnalysis/DoBinningHadron
  ! NAME
  ! subroutine DoBinningHadron(Part,nu,Q2,Ebeam,phi_Lepton,iPart)
  ! PURPOSE
  ! By calling this routine, the histogramms according the hadron are filled.
  ! INPUTS
  ! * type(particle) :: Part -- the particle
  ! * real :: nu -- nu value of Event
  ! * real :: Q2 -- Q**2 value of Event
  ! * real :: Ebeam -- Energy of the Lepton Beam
  ! * real :: phi_Lepton -- Angle of scattered Lepton
  ! * integer :: iPart -- number of particle in pertPart
  ! OUTPUT
  ! ---
  !****************************************************************************
  subroutine DoBinningHadron(Part,nu,Q2,Ebeam,phi_Lepton,iEns,iPart)
    use particleDefinition
    use IdTable
    use history
    use PIL_rho0Dec
    use output

    type(particle),intent(in) :: Part
    real,intent(in) :: nu,Q2,Ebeam,phi_Lepton
    integer, intent(in) :: iEns,iPart

    integer :: ID, IZ
    logical :: acceptFlag
    real :: Accweight,zh, W0,W1,W2, pp
    real :: pT2   ! pT^2 according photon momentum axis
    real :: pT2L  ! pT^2 according lepton momentum axis
    integer :: i, iClass, iEns1, iNuc, iGen
!    integer, dimension(1:3) :: parents

    call checkCuts(Part,nu,Q2,Ebeam,phi_Lepton,acceptFlag,AccWeight, pT2L)
    if (.not.acceptFlag) return

    ID=Part%ID
    IZ=Part%charge
    zh=Part%momentum(0)/nu
    pp = absMom(Part)

    pT2 = Part%momentum(1)**2 + Part%momentum(2)**2 ! pT^2 according photon axis
!      write(*,*) 'pT2:',pT2,pT2L
!      pT2 = pT2L                                      ! pT^2 according lepton axis


    W0 = Part%perWeight
    W1 = W0 * AccWeight ! weight for charged particles
    W2 = W1             ! weight for identified spectra

    select case (iExperiment)

    case (1) ! Hermes (H,N,Kr)

       select case (ID)
       case (nucleon)
          if (pp <  4.0) W2 = 0.0
          if (pp > 15.0) W2 = 0.0
       case default
          if (pp <  2.5) W2 = 0.0
          if (pp > 15.0) W2 = 0.0
       end select

    case (2) ! Hermes (Ne)

       select case (ID)
       case (nucleon)
          if (pp <  4.0) W2 = 0.0
          if (pp > 15.0) W2 = 0.0
       case (pion)
          if (pp <  0.6) W2 = 0.0
          if (pp > 15.0) W2 = 0.0
       case (kaon,kaonBar)
          if (pp <  2.0) W2 = 0.0
          if (pp > 15.0) W2 = 0.0
       case default
          if (pp <  2.5) W2 = 0.0
          if (pp > 15.0) W2 = 0.0
       end select

    case (11) ! Hermes (final paper)

       if (pp > 15.0) W2 = 0.0
       if (pp <  2.0) W2 = 0.0
       if (id.eq.pion .and. iz.eq.0 .and. pp <  2.5) W2 = 0.0

    case (13) ! Hermes (pT-broadening)

       select case (ID)
       case (pion,kaon,kaonBar)
          if (pp <  2.0) W2 = 0.0
          if (pp > 15.0) W2 = 0.0
       case default
          if (pp <  1.0) W2 = 0.0
       end select

    case (10) ! Hermes@12GeV

       if (pp <  1.0) W2 = 0.0

    end select

    !-----------------------------------------------------------------
    ! determine average kinematics (zH)
    !-----------------------------------------------------------------

    call AddHist(hKin_zH_Q2(1),zH, W0,Q2*W0) ! without acceptance
    call AddHist(hKin_zH_Q2(2),zH, W1,Q2*W1) ! with acceptance

    call AddHist(hKin_zH_nu(1),zH, W0,nu*W0) ! without acceptance
    call AddHist(hKin_zH_nu(2),zH, W1,nu*W1) ! with acceptance

    !-----------------------------------------------------------------
    ! determine hadron sepectra (zH)
    !-----------------------------------------------------------------

    if (iZ.ne.0) then
       call AddHist(hCH_zH, zH, W1,W0)
    end if

    call AddHistMP(HMP_zH,   Part, zH, W2,W0)
    call AddHistMP(HMP_pT2zH,Part, zH, W2,W2*pT2)

    if (ID==pion) then
       call AddHist2D(H2DpTPionZH(IZ), (/zH,pT2/),W2)
    end if

    iEns1 = Part%firstEvent / 1000
    iNuc = mod(Part%firstEvent,1000)

!    iClass = 0
!    if (Part%history.eq.0) then
!       iClass = 1 ! direct production
!    else if (history_1Body(Part%history)) then
!       iClass = 2 ! Decay
!       call history_getParents(Part%history,parents)
!       if (parents(1).eq.103) then
!          iClass = 3 ! Decay of rho
!          if (PIL_rho0Dec_Get(Part%number,i)) then
!             iClass = 4 ! Decay of diffr. rho0
!          endif
!       endif
!    else
!       iClass = 5 ! something else
!    endif

    if (DoClassifyFirst) then
       call AddHistMP(HMP_zHClass(0),     Part, zH,W2)
       call AddHistMP(HMP_pT2Class(0),    Part, zH,W2,W2*pT2)
       iClass = ClassifyFirstEvent(iEns1,iNuc)
       if (iClass.gt.0) then
          call AddHistMP(HMP_zHClass(iClass), Part, zH,W2)
          call AddHistMP(HMP_pT2Class(iClass),Part, zH,W2,W2*pT2)
       end if

       iGen = history_getGeneration(Part%history)
       call AddHistMP(hMP_zH_Generation(-1),Part, zH, W2)
       if (iGen.lt.11) then
          call AddHistMP(hMP_zH_Generation(iGen),Part, zH, W2)
       end if
    end if

    select case (iExperiment)
    case (5)
       call DoBinning_5
    case (16)
       call DoBinning_16
    end select

    if (CollHist_GetDoCollHistory()) then
       if ((ID==pion).and.(IZ/=0)) then

          iClass = CollHist_ClassifyHist(iEns,iPart)

!!$          if (iClass.eq.0) then
!!$             write(*,*) '?????????',iEns,iPart
!!$             call WriteParticle(6,iEns,iPart,Part)
!!$             call CollHist_WriteHistParticle(6,iEns,iPart)
!!$          else
!!$             write(*,*) '!!!!!!!!!',iEns,iPart
!!$          end if

          iGen = 1
          if (abs(iClass).ge.100) then
             iGen = 2
             iClass = iClass/100
          end if

          call AddHist2D(H2D_CollHistPT(iGen), (/iClass*1.0, zH/), W2,W2*pT2)

       end if
    end if

    if (zh.lt.zHR%l) return ! <<<<<<<<<<<<<<<<<<<<<<<< minimal zH !

    !-----------------------------------------------------------------
    ! determine average kinematics (nu)
    !-----------------------------------------------------------------

    call AddHist(hKin_nu_Q2(1),nu, W0,Q2*W0) ! without acceptance
    call AddHist(hKin_nu_Q2(2),nu, W1,Q2*W1) ! with acceptance

    call AddHist(hKin_nu_zH(1),nu, W0,zH*W0) ! without acceptance
    call AddHist(hKin_nu_zH(2),nu, W1,zH*W1) ! with acceptance

    !-----------------------------------------------------------------
    ! determine hadron sepectra (nu,Q2,pT2)
    !-----------------------------------------------------------------

    if (iZ.ne.0) then
       call AddHist(hCH_nu,  nu,  W1,W0)
       call AddHist(hCH_Q2,  Q2,  W1,W0)
       call AddHist(hCH_pT2, pT2, W1,W0)
    end if

    call AddHistMP(HMP_nu, Part, nu,  W2,W0)
    call AddHistMP(HMP_Q2, Part, Q2,  W2,W0)
    call AddHistMP(HMP_pT2,Part, pT2, W2,W0)

    call AddHistMP(HMP_pT2nu, Part, nu,  W2,W2*pT2)
    call AddHistMP(HMP_pT2Q2, Part, Q2,  W2,W2*pT2)
    call AddHistMP(HMP_pT2pT2,Part, pT2, W2,W2*pT2)


    if (ID==pion) then
       call AddHist2D(H2DpTAvePion(IZ), (/nu,zH/),W2,W2*pT2 )

       call AddHist2D(H2DpTPionNU(IZ), (/nu,pT2/),W2)
       call AddHist2D(H2DpTPionQ2(IZ), (/Q2,pT2/),W2)

    end if

    if (iExperiment.eq.11) then ! HERMES (final paper)
       ! please note: we have already zH > 0.2. Is this correct???

       call DoBinning_11

    end if

  contains

    !-----------------------------------------------------------------
    subroutine DoBinning_5
      integer :: iNu,iQ2,iZH

      if (ID/=pion) return ! only Pions
      if (IZ.ne. 1)  return ! only positive charged

      call BinCutsJLAB5_R(nu,Q2,zH, iNu,iQ2,iZH)

      if (min(iQ2,iZH)>0) call AddHist(hJLAB5_R_nu(iQ2,iZH), nu,W2,W0)
      if (min(iNu,iZH)>0) call AddHist(hJLAB5_R_Q2(iNu,iZH), Q2,W2,W0)
      if (min(iQ2,iNu)>0) call AddHist(hJLAB5_R_zH(iQ2,iNu), zH,W2,W0)

      call BinCutsJLAB5_pT2(nu,Q2,zH, iNu,iQ2,iZH)

      if (min(iNu,iQ2,iZH)>0) hJLAB5_pT2_ARR(iNu,iQ2,iZH,:)= hJLAB5_pT2_ARR(iNu,iQ2,iZH,:)+(/1.0,W2,W2*pT2/)

    end subroutine DoBinning_5

    !-----------------------------------------------------------------
    subroutine DoBinning_11

      if (zH.lt.0.4) then
         i=1
      else if (zH.lt.0.7) then
         i=2
      else
         i=3
      end if
      call AddHistMP(HMP_nu_zH(i), Part, nu,  W2,W0)
      call AddHistMP(HMP_Q2_zH(i), Part, Q2,  W2,W0)
      call AddHistMP(HMP_pT2_zH(i),Part, pT2, W2,W0)

      if (nu.lt.12.0) then
         i=1
      else if (nu.lt.17.0) then
         i=2
      else
         i=3
      end if
      call AddHistMP(HMP_zH_nu(i), Part, zH,  W2,W0)
      call AddHistMP(HMP_Q2_nu(i), Part, Q2,  W2,W0)
      call AddHistMP(HMP_pT2_nu(i),Part, pT2, W2,W0)

      if (pT2 < 0.7) then
         i=1
      else
         i=2
      end if
      call AddHistMP(HMP_nu_pT(i), Part, nu,  W2,W0)
      call AddHistMP(HMP_zH_pT(i), Part, zH,  W2,W0)
      call AddHistMP(HMP_Q2_pT(i), Part, Q2,  W2,W0)


    end subroutine DoBinning_11

    !-----------------------------------------------------------------
    subroutine DoBinning_16
      integer :: iQ2

      iQ2=0
      if (Q2.lt.2.0) then
         iQ2 = 1
      else if (Q2.lt.4.0) then
         iQ2 = 2
      else if (Q2.lt.6.0) then
         iQ2 = 3
      else if (Q2.lt.10.0) then
         iQ2 = 4
      else
         iQ2 = 5
      end if

      call AddHistMP(HMP_EIC_zH, Part, zH,  W2,W0)
      call AddHistMP(HMP_EIC_zH_Q2(iQ2), Part, zH, W2,W0)

      if (zh.lt.zHR%l) return

      call AddHistMP(HMP_EIC_nu, Part, nu,  W2,W0)
      call AddHistMP(HMP_EIC_nu_Q2(iQ2), Part, nu,  W2,W0)

      call AddHistMP(HMP_EIC_Q2, Part, Q2,  W2,W0)

    end subroutine DoBinning_16
    !-----------------------------------------------------------------

  end subroutine DoBinningHadron

  !****************************************************************************
  !****is* HiLeptonAnalysis/writeBinning
  ! NAME
  ! subroutine writeBinning()
  ! PURPOSE
  ! write out all histogramms
  !
  ! INPUTS
  ! (none)
  !
  ! OUTPUT
  ! several files are (re-)written.
  !
  ! NOTES
  !
  !****************************************************************************

  subroutine writeBinning()
    use output

    integer :: i
    real :: mul0

    !type(histogram)   :: HHH

    mul0 = 1./nLeptons

    if (DoLeptonKinematics) then

       ! ===== Lepton Kinematics =====

       call WriteHist(hLep_Q2,add=add0,&
            &file='HiLep.lep.Q2.kinematics.dat',dump=.true.)
       call WriteHist(hLep_nu,add=add0,&
            &file='HiLep.lep.nu.kinematics.dat',dump=.true.)

       call WriteHist(hLep_Q2,add=add0,DoAve=.true.,&
            &file='HiLep.lep.nuAve_Q2.kinematics.dat')
       call WriteHist(hLep_nu,add=add0,DoAve=.true.,&
            &file='HiLep.lep.Q2Ave_nu.kinematics.dat')


       ! -----

       open(140,file='HiLep.NuQ2planeXS.SYS'//'.dat', status='unknown')
       rewind(140)
       do i=0,9
          call WriteHist2D_SYSTEM(H2DleptonXS(i),140)

          call WriteHist2D_Gnuplot(H2DleptonXS(i),mul=mul0, add=add0,&
               &file='HiLep.NuQ2planeXS.'//trim(intToChar(i))//'.dat',dump=.true.)
       end do
       close(140)

    end if

    if (DoHadronKinematics) then

       ! ===== nu-spectra =====

       call WriteHistMP(hMP_nu,add=add0,H2=hLep_nu,iColumn=3,&
            &file='HiLep.nu.idH.noAcc.dat') ! noAcc
       call WriteHistMP(hMP_nu,add=add0,H2=hLep_nu,iColumn=1,&
            &file='HiLep.nu.idH.Acc.dat',dump=.true.) ! Acc

       call WriteHist(hCH_nu,add=add0,H2=hLep_nu,&
            &file='HiLep.nu.chH.dat',dump=.true.) ! Acc, NoAcc

       call WriteHist(hKin_nu_Q2(1),add=add0,DoAve=.true.,&
            &file='HiLep.nu.Q2.kinematics.noAcc.dat',dump=.true.) ! <Q2>(nu), noAcc
       call WriteHist(hKin_nu_Q2(2),add=add0,DoAve=.true.,&
            &file='HiLep.nu.Q2.kinematics.Acc.dat',dump=.true.) ! <Q2>(nu), Acc

       call WriteHist(hKin_nu_zH(1),add=add0,DoAve=.true.,&
            &file='HiLep.nu.zH.kinematics.noAcc.dat',dump=.true.) ! <zH>(nu), noAcc
       call WriteHist(hKin_nu_zH(2),add=add0,DoAve=.true.,&
            &file='HiLep.nu.zH.kinematics.Acc.dat',dump=.true.) ! <zH>(nu), Acc

       call WriteHistMP(hMP_pT2nu,DoAve=.true.,&
         &file='HiLep.AvePT2.nu.dat',dump=.true.)

       ! ===== Q2-spectra =====

       call WriteHistMP(hMP_Q2,add=add0,H2=hLep_Q2,iColumn=3,&
            &file='HiLep.Q2.idH.noAcc.dat') ! noAcc
       call WriteHistMP(hMP_Q2,add=add0,H2=hLep_Q2,iColumn=1,&
            &file='HiLep.Q2.idH.Acc.dat',dump=.true.) ! Acc

       call WriteHist(hCH_Q2,add=add0,H2=hLep_Q2,&
            &file='HiLep.Q2.chH.dat',dump=.true.) ! Acc, NoAcc

       call WriteHistMP(hMP_pT2Q2,DoAve=.true.,&
            &file='HiLep.AvePT2.Q2.dat',dump=.true.)

       ! ===== zH-spectra =====

       call WriteHistMP(hMP_zH,add=add0,mul=mul0,iColumn=3,&
            &file='HiLep.zH.idH.noAcc.dat') ! noAcc
       call WriteHistMP(hMP_zH,add=add0,mul=mul0,iColumn=1,&
            &file='HiLep.zH.idH.Acc.dat',dump=.true.) ! Acc

       call WriteHist(hCH_zH,add=add0,mul=mul0,&
            &file='HiLep.zH.chH.dat',dump=.true.) ! Acc, NoAcc

       call WriteHist(hKin_zH_Q2(1),add=add0,DoAve=.true.,&
            &file='HiLep.zH.Q2.kinematics.noAcc.dat',dump=.true.) ! <Q2>(zH), noAcc
       call WriteHist(hKin_zH_Q2(2),add=add0,DoAve=.true.,&
            &file='HiLep.zH.Q2.kinematics.Acc.dat',dump=.true.) ! <Q2>(zH), Acc

       call WriteHist(hKin_zH_nu(1),add=add0,DoAve=.true.,&
            &file='HiLep.zH.nu.kinematics.noAcc.dat',dump=.true.) ! <nu>(zH), noAcc
       call WriteHist(hKin_zH_nu(2),add=add0,DoAve=.true.,&
            &file='HiLep.zH.nu.kinematics.Acc.dat',dump=.true.) ! <nu>(zH), Acc

       call WriteHistMP(hMP_pT2zH,DoAve=.true.,&
            &file='HiLep.AvePT2.zH.dat',dump=.true.)

    end if

!    rewind(56)
!    call WriteHistMP(hMP_zHClass(0),56,add=add0,mul=mul0,iColumn=1) ! Acc
!    do i=1,5
!       rewind(560+i)
!       call WriteHistMP(hMP_zHClass(i),560+i,add=add0,mul=mul0,iColumn=1) ! Acc
!    end do

    if (DoClassifyFirst) then
       do i=0,8
          call WriteHistMP(hMP_zHClass(i),add=add0,mul=mul0,iColumn=1,&
               & file='HiLep.zHclass.'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc
          call WriteHistMP(hMP_pT2Class(i),DoAve=.true.,&
               & file='HiLep.pT2class.'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc
       end do

       call WriteHistMP(hMP_zH_Generation(-1),add=add0,mul=mul0,iColumn=1,&
            & file='HiLep.zH_Generation.'//intToChar(999)//'.dat',dump=.true.)
       do i=0,10
          call WriteHistMP(hMP_zH_Generation(i),add=add0,mul=mul0,iColumn=1,&
               & file='HiLep.zH_Generation.'//intToChar(i)//'.dat',dump=.true.)
       end do
    end if


    if (DoHadronKinematics) then
       ! ===== pT2-spectra =====

       call WriteHistMP(hMP_pT2,add=add0,mul=mul0,iColumn=3,&
            &file='HiLep.pT2.idH.noAcc.dat') ! noAcc
       call WriteHistMP(hMP_pT2,add=add0,mul=mul0,iColumn=1,&
            &file='HiLep.pT2.idH.Acc.dat',dump=.true.) ! Acc

       call WriteHist(hCH_pT2,add=add0,mul=mul0,&
            &file='HiLep.pT2.chH.dat',dump=.true.) ! Acc, NoAcc

       call WriteHistMP(hMP_pT2pT2,DoAve=.true.,&
            &file='HiLep.AvePT2.pT2.dat',dump=.true.)

       ! =====================

       do i=-1,1
          call WriteHist2D_Gnuplot(H2DpTAvePion(i),DoAve=.true.,maxval=0.0,&
               &file='HiLep.pT2Ave_Plane.'//piName(i)//'.dat',dump=.true.)

          call WriteHist2D_Gnuplot(H2DpTPionZH(i), mul=mul0, add=add0,&
               &file='HiLep.pT2_zH.'//piName(i)//'.dat',dump=.true.)

          call WriteHist2D_Gnuplot(H2DpTPionNU(i), H2=hLep_nu, add=add0,&
               &file='HiLep.pT2_nu.'//piName(i)//'.dat',dump=.true.)

          call WriteHist2D_Gnuplot(H2DpTPionQ2(i), H2=hLep_Q2, add=add0,&
               &file='HiLep.pT2_Q2.'//piName(i)//'.dat',dump=.true.)
       end do


       call WriteHist2D_Gnuplot(H2D_CollHistPT(1), DoAve=.true., mul=mul0, add=add0,&
            & file='HiLep.CollHistPT2.1.dat',dump=.true.) ! Acc
       call WriteHist2D_Gnuplot(H2D_CollHistPT(2), DoAve=.true., mul=mul0, add=add0,&
            & file='HiLep.CollHistPT2.2.dat',dump=.true.) ! Acc



       select case (iExperiment)
       case (11) ! HERMES (final paper)

          call WriteBinning_11

       case (5) ! JLAB@5GeV (CLAS)

          call WriteBinning_5

       case (16) ! EIC

          !       call WriteHistMP(hMP_EIC_zH,add=add0,mul=mul0,iColumn=3,&
          !            &file='HiLep.EIC_zH.idH.noAcc.dat') ! noAcc
          call WriteHistMP(hMP_EIC_zH,add=add0,mul=mul0,iColumn=1,&
               &file='HiLep.EIC_zH.idH.Acc.dat',dump=.true.) ! Acc

          !       call WriteHistMP(hMP_EIC_nu,add=add0,H2=hLep_nu,iColumn=3,&
          !            &file='HiLep.EIC_nu.idH.noAcc.dat') ! noAcc
          call WriteHistMP(hMP_EIC_nu,add=add0,H2=hLep_nu,iColumn=1,&
               &file='HiLep.EIC_nu.idH.Acc.dat',dump=.true.) ! Acc

          !       call WriteHistMP(hMP_EIC_Q2,add=add0,H2=hLep_Q2,iColumn=3,&
          !            &file='HiLep.EIC_Q2.idH.noAcc.dat') ! noAcc
          call WriteHistMP(hMP_EIC_Q2,add=add0,H2=hLep_Q2,iColumn=1,&
               &file='HiLep.EIC_Q2.idH.Acc.dat',dump=.true.) ! Acc

          do i=1,5
             call WriteHistMP(hMP_EIC_zH_Q2(i),add=add0,mul=mul0,iColumn=1,&
                  &file='HiLep.EIC_zH.idH.Q2_'//Achar(i+48)//'.Acc.dat',dump=.true.) ! Acc
             call WriteHistMP(hMP_EIC_nu_Q2(i),add=add0,H2=hLep_nu,iColumn=1,&
                  &file='HiLep.EIC_nu.idH.Q2_'//Achar(i+48)//'.Acc.dat',dump=.true.) ! Acc
          end do


       end select

    end if


    contains
      !-----------------------------------------------------------------
      subroutine WriteBinning_5
        integer :: iNu,iQ2,iZH, nNu,nQ2,nZH

        call BinCutsJLAB5_R(-99.0,99.0,99.0, nNu,nQ2,nZH)

        do iQ2=1,nQ2
           do iZH=1,nZH
              call WriteHist(hJLAB5_R_nu(iQ2,iZH), add=add0,mul=mul0, &
                   & file='HiLep.JLAB5_R.nu_'//Achar(iQ2+48)//'_'//Achar(iZH+48)//'.dat',dump=.true.)
           end do
        end do

        do iNu=1,nNu
           do iZH=1,nZH
              call WriteHist(hJLAB5_R_Q2(iNu,iZH), add=add0,mul=mul0, &
                   & file='HiLep.JLAB5_R.Q2_'//Achar(iNu+48)//'_'//Achar(iZH+48)//'.dat',dump=.true.)
           end do
        end do

        do iQ2=1,nQ2
           do iNu=1,nNu
              call WriteHist(hJLAB5_R_zH(iQ2,iNu), add=add0,mul=mul0, &
                   & file='HiLep.JLAB5_R.zH_'//Achar(iQ2+48)//'_'//Achar(iNu+48)//'.dat',dump=.true.)
           end do
        end do

        call BinCutsJLAB5_pT2(-99.0,99.0,99.0, nNu,nQ2,nZH)

        open(141,file='HiLep.JLAB5_pT2.ARR_Q2_nu_zH.dat')
        rewind(141)
        do iQ2=1,nQ2
           do iNu=1,nNu
              do iZH=1,nZH
                 write(141,'(3i2,0P,f11.5,1P,3e12.4,0P)') iQ2,iNu,iZH,&
                      & hJLAB5_pT2_ARR(iNu,iQ2,iZH,2)/(hJLAB5_pT2_ARR(iNu,iQ2,iZH,1)+add0),&
                      & hJLAB5_pT2_ARR(iNu,iQ2,iZH,:)
              end do
           end do
        end do
        close(141)


      end subroutine WriteBinning_5

      !-----------------------------------------------------------------
      subroutine WriteBinning_11

        do i=1,3
           call WriteHistMP(hMP_nu_zh(i),add=add0,mul=mul0,iColumn=1,&
                &file='HiLep.dN_id.nu.zH.'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc
           call WriteHistMP(hMP_Q2_zh(i),add=add0,mul=mul0,iColumn=1,&
                &file='HiLep.dN_id.Q2.zH.'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc
           call WriteHistMP(hMP_pT2_zh(i),add=add0,mul=mul0,iColumn=1,&
                &file='HiLep.dN_id.pT2.zH.'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc

           call WriteHistMP(hMP_zH_nu(i),add=add0,mul=mul0,iColumn=1,&
                &file='HiLep.dN_id.zH.nu.'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc
           call WriteHistMP(hMP_Q2_nu(i),add=add0,mul=mul0,iColumn=1,&
                &file='HiLep.dN_id.Q2.nu.'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc
           call WriteHistMP(hMP_pT2_nu(i),add=add0,mul=mul0,iColumn=1,&
                &file='HiLep.dN_id.pT2.nu.'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc
        end do

        do i=1,2
           call WriteHistMP(hMP_nu_pT(i),add=add0,mul=mul0,iColumn=1,&
                &file='HiLep.dN_id.nu.pT2.'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc
           call WriteHistMP(hMP_zH_pT(i),add=add0,mul=mul0,iColumn=1,&
                &file='HiLep.dN_id.zH.pT2.'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc
           call WriteHistMP(hMP_Q2_pT(i),add=add0,mul=mul0,iColumn=1,&
                &file='HiLep.dN_id.Q2.pT2.'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc
        end do


      end subroutine WriteBinning_11
      !-----------------------------------------------------------------

  end subroutine writeBinning


  !****************************************************************************
  !****is* HiLeptonAnalysis/BinCutsJLAB5_R
  ! NAME
  ! subroutine BinCutsJLAB5_R(nu,Q2,zH, iNu,iQ2,iZH)
  ! PURPOSE
  ! select the histogramm to bin the values of R
  !****************************************************************************

  subroutine BinCutsJLAB5_R(nu,Q2,zH, iNu,iQ2,iZH)

    real, intent(IN)     :: nu,Q2,zH
    integer, intent(OUT) :: iNu,iQ2,iZH

    integer, parameter :: nNu=6, nQ2=6, nZH=6
    real, dimension(nNu), parameter :: cutsNu = (/2.2, 3.0, 3.3, 3.5, 4.0, 4.3/)
    real, dimension(nQ2), parameter :: cutsQ2 = (/1.0, 1.25,1.5,1.85, 2.4, 4.2/)
    real, dimension(nZH), parameter :: cutsZH = (/0.0, 0.2, 0.4, 0.6, 0.8, 1.0/)

    integer :: i

    if (nu<0.0) then ! return max histogram number
       iNu = nNu-1
       iQ2 = nQ2-1
       iZH = nZH-1
       return
    end if

    iNu = 0
    do i=1,nNu
       if (nu.lt.cutsNu(i)) exit
       iNu = iNu+1
    end do
    if (iNu.eq.nNu) iNu=0 ! failure

    iQ2 = 0
    do i=1,nQ2
       if (Q2.lt.cutsQ2(i)) exit
       iQ2 = iQ2+1
    end do
    if (iQ2.eq.nQ2) iQ2=0 ! failure

    iZH = 0
    do i=1,nZH
       if (ZH.lt.cutsZH(i)) exit
       iZH = iZH+1
    end do
    if (iZH.eq.nZH) iZH=0 ! failure

  end subroutine BinCutsJLAB5_R

  !****************************************************************************
  !****is* HiLeptonAnalysis/BinCutsJLAB5_pT2
  ! NAME
  ! subroutine BinCutsJLAB5_pT2(nu,Q2,zH, iNu,iQ2,iZH)
  ! PURPOSE
  ! select the histogramm to bin the values of R
  !****************************************************************************

  subroutine BinCutsJLAB5_pT2(nu,Q2,zH, iNu,iQ2,iZH)

    real, intent(IN)     :: nu,Q2,zH
    integer, intent(OUT) :: iNu,iQ2,iZH

    integer, parameter :: nNu=4, nQ2=4, nZH=5
    real, dimension(nNu), parameter :: cutsNu = (/2.0, 3.0, 4.0, 5.0/)
    real, dimension(nQ2), parameter :: cutsQ2 = (/1.0, 2.0, 3.0, 4.0/)
    real, dimension(nZH), parameter :: cutsZH = (/0.4, 0.5, 0.6, 0.7, 0.8/)

    integer :: i

    if (nu<0.0) then ! return max histogram number
       iNu = nNu-1
       iQ2 = nQ2-1
       iZH = nZH-1
       return
    end if

    iNu = 0
    do i=1,nNu
       if (nu.lt.cutsNu(i)) exit
       iNu = iNu+1
    end do
    if (iNu.eq.nNu) iNu=0 ! failure

    iQ2 = 0
    do i=1,nQ2
       if (Q2.lt.cutsQ2(i)) exit
       iQ2 = iQ2+1
    end do
    if (iQ2.eq.nQ2) iQ2=0 ! failure

    iZH = 0
    do i=1,nZH
       if (ZH.lt.cutsZH(i)) exit
       iZH = iZH+1
    end do
    if (iZH.eq.nZH) iZH=0 ! failure

  end subroutine BinCutsJLAB5_pT2

  !****************************************************************************
  !****is* HiLeptonAnalysis/checkCuts
  ! NAME
  ! subroutine checkCuts(Part,nu,Q2,Ebeam,phi_Lepton,acceptFlag,weight)
  ! PURPOSE
  ! This routine applies the experimental cuts and returns some weight,
  ! which indicate the probability, that the parton is detected.
  !
  ! As a side effect, it calculates the transverse momentum of the
  ! particle in respect to the photon direction.
  ! INPUTS
  ! * type(particle) :: Part -- particle to consider
  ! * real :: nu -- nu value of Event
  ! * real :: Q2 -- Q**2 value of Event
  ! * real :: Ebeam -- Energy of the Lepton Beam
  ! * real :: phi_Lepton -- Angle of scattered Lepton
  ! OUTPUT
  ! * logical :: acceptFlag -- true, if particle is detected at all
  ! * real :: weight -- weight of acceptance
  ! * real, OPTIONAL :: pT2 -- transverse momentum according lepton beam
  ! NOTES
  ! phi_Lepton indicates the azimuthal angle of the scattered lepton
  ! around the axis, which is determined by the momentum of the
  ! incoming lepton.
  ! This is normally not the z-axis. (see below)
  !
  ! A particle can be detected, but, e.g. due to rotational systematics
  ! around the lepton beam axis, only in 50% of all cases. If the details
  ! of the detector are not important, you can then give the particle
  ! the acceptance weight 0.5.
  !
  ! Be aware: The events are normally in a system, where the PHOTON
  ! determines the z-axis.
  !****************************************************************************
  subroutine checkCuts(Part,nu,Q2,Ebeam,phi_Lepton,acceptFlag,weight, pT2)
    use particleDefinition
    use initHiLepton, only: AccWeight

    type(particle),intent(in) :: Part
    logical,intent(out) :: acceptFlag
    real,intent(in) :: nu,Q2,Ebeam,phi_Lepton
    real,intent(out) :: weight
    real,intent(out),optional :: pT2

    real :: pBeam,pPhoton,alpha,theta, h
!    real :: Eprime,pPrime
!    real,parameter :: mLepton=0. ! neglect lepton mass just as in initLepton


    weight=0.
    acceptFlag=.false.
    if (Part%momentum(0).lt.EHR%l) return

    acceptFlag=.true.

    pPhoton=sqrt(nu**2+Q2)
    !Eprime=Ebeam-nu
    !pPrime=sqrt(Eprime**2-mLepton**2)
    !pBeam=sqrt(Ebeam**2-mLepton**2)
    pBeam=Ebeam
    !alpha=acos((pPhoton**2+pBeam**2-pPrime**2)/(2.*pBeam*pPhoton)) ! angle between photon momentum and beam axis
    ! The following is the same:
    alpha=acos((2*Ebeam*nu+Q2)/(2*pBeam*pPhoton)) ! angle between photon momentum and beam axis


    call labrotate(alpha,phi_Lepton,Part%momentum,theta,pT2=h)
    if (present(pT2))  pT2 =  h

    weight=AccWeight(iDetector,sign(1,Part%Charge),absMom(Part),theta) !theta is angle between hadron momentum and beam axis



  end subroutine checkCuts


  !****************************************************************************
  !****is* HiLeptonAnalysis/labrotate
  ! NAME
  ! subroutine labrotate(alpha,phi_Lepton,p,theta,phi,pT2)
  ! PURPOSE
  ! Because AccWeight is an rotation-averaged value, we have
  ! now to choose the "phi" angle.
  ! INPUTS
  ! * alpha      : angle to rotate around y-axis
  ! * phi_Lepton : angle to rotate arund z-axis
  ! * p          : vector to rotate
  ! OUTPUT
  ! * theta  : angle
  ! * phi    : angle
  ! * pT2    : transverse momentum (squared) acording lepton axis
  !
  ! NOTES
  ! First the vector is rotated around the y-axis by the angle
  ! alpha (i.e. the angle between incoming lepton and photon).
  ! Secondly, this vector is rotated around the z-axis by the angle
  ! phi_Lepton.
  !****************************************************************************
  subroutine labrotate(alpha,phi_Lepton,p,theta,phi,pT2)

    real,intent(in) :: alpha,phi_Lepton
    real,dimension(0:3),intent(in) :: p
    real,intent(out) :: theta
    real,intent(out),optional :: phi
    real,intent(out),optional :: pT2

    real :: cp,sp,ca,sa,plx,ply,plz,ptot

    ptot=sqrt(p(1)**2+p(2)**2+p(3)**2)

    cp=cos(phi_Lepton)
    sp=sin(phi_Lepton)

    ca=cos(alpha)
    sa=sin(alpha)

    plx= p(1)*ca*cp +p(2)*sp + p(3)*sa*cp
    ply=-p(1)*ca*sp +p(2)*cp - p(3)*sa*sp
    plz=-p(1)*sa             + p(3)*ca

    theta=acos(plz/ptot)

    if (present(phi)) phi=atan2(ply,plx)
    if (present(pT2)) pT2=plx**2 + ply**2

  end subroutine labrotate


  !****************************************************************************
  !****s* HiLeptonAnalysis/HiLeptonAnalysisPerTime
  ! NAME
  ! subroutine HiLeptonAnalysisPerTime(Time, pParts)
  ! PURPOSE
  ! produce statistical output after every time step
  !
  ! INPUTS
  ! * real           :: Time        -- time of time step (in fm)
  ! * type(particle) :: pParts(:,:) -- particle vector of perturbative parts
  !
  ! OUTPUT
  ! internal variables are changed and some output to predefined files.
  !
  !****************************************************************************
  subroutine HiLeptonAnalysisPerTime(Time, pParts)
    use particleDefinition
    use hist2Df90
    use EventInfo_HiLep

    real,    intent(in)       :: Time
    type(particle),intent(in) :: pParts(:,:)

    logical, save :: DoInit = .TRUE.

    type(histogram2D), save :: hist_Time, hist_Pos

    integer :: i,j
    real :: r, Eh, zH, w, ws
    real :: weight,nu,Q2,epsilon
    integer :: EventType

    if (DoInit) then

       call CreateHist2D(hist_Time, "XS vs. zH and TIME",&
            & (/0.0,0.0/), (/1.1,30.0/), (/0.02,0.5/)  , .true.)
       call CreateHist2D(hist_Pos , "XS vs. zH and POS",&
            & (/0.0,0.0/), (/1.1,30.0/), (/0.02,0.5/)  , .true.)

       DoInit = .FALSE.
    end if


    do j=1,Size(pParts,dim=2)
       do i=1,Size(pParts,dim=1)
          if (pParts(i,j)%ID > 0) then

             Eh = pParts(i,j)%momentum(0)

             if (EventInfo_HiLep_Get(i,pParts(i,j)%firstEvent,Weight,nu,Q2,epsilon,EventType)) then
             else
                write(*,*) 'Ooops, L1087'
                stop
             end if

             zH = Eh /nu

             w  = pParts(i,j)%perWeight
             ws = w * pParts(i,j)%ScaleCS

             call AddHist2D(hist_Time, (/zh,time/), ws, w)

             r = sqrt(DOT_PRODUCT(pParts(i,j)%position,pParts(i,j)%position))
             call AddHist2D(hist_Pos, (/zh,r/), ws, w)
          end if
       end do
    end do

    call WriteHist2D_Gnuplot(hist_Time, file='HiLeptonAnalysisPerTime.611.dat')
    call WriteHist2D_Gnuplot(hist_Pos,  file='HiLeptonAnalysisPerTime.612.dat')

  end subroutine HiLeptonAnalysisPerTime

  !****************************************************************************

  subroutine GlobalAverage_ZERO()

    GlobalAverage%sum0 = 0
    GlobalAverage%sum  = 0
    GlobalAverage%W    = 0
    GlobalAverage%Q2   = 0
    GlobalAverage%nu   = 0

  end subroutine GlobalAverage_ZERO

  !****************************************************************************

  subroutine GlobalAverage_ADD(weight, Q2, nu)
    use constants, only: mN

    real, intent(in) :: weight, Q2, nu

    real :: W

    W = sqrt(max(mN**2+2*mN*nu-Q2, 0.))

    GlobalAverage%sum0 = GlobalAverage%sum0 + 1.0
    GlobalAverage%sum  = GlobalAverage%sum  + weight
    GlobalAverage%W    = GlobalAverage%W    + weight*W
    GlobalAverage%Q2   = GlobalAverage%Q2   + weight*Q2
    GlobalAverage%nu   = GlobalAverage%nu   + weight*nu

  end subroutine GlobalAverage_ADD

  !****************************************************************************

  subroutine GlobalAverage_Write(iFile,nRuns,nEns,massNum)
    integer, intent(in) :: iFile
    integer, intent(in) :: nRuns,nEns,massNum

    if (GlobalAverage%sum.le.0.) return ! no entries

    write(iFile,'(A)') '======== GlobalAverages ========== '
    write(iFile,'(A,g14.5)') ' W  [GeV]    :      ',GlobalAverage%W/GlobalAverage%sum
    write(iFile,'(A,g14.5)') ' Q2 [GeV^2]  :      ',GlobalAverage%Q2/GlobalAverage%sum
    write(iFile,'(A,g14.5)') ' nu [GeV]    :      ',GlobalAverage%nu/GlobalAverage%sum
    write(iFile,'(A)') '================================== '
    write(iFile,'(A,g14.5)') ' sigma_tot   [mub]: ',GlobalAverage%sum/(1.0*nRuns*nEns)
    write(iFile,'(A,g14.5)') ' sigma_tot/A [mub]: ',GlobalAverage%sum/(1.0*nRuns*nEns*massNum)
    write(iFile,'(A)') '================================== '
    write(iFile,'(A,g14.5)') ' N_lepton   : ',GlobalAverage%sum0/(1.0*nRuns*nEns)
    write(iFile,'(A,g14.5)') ' N_lepton/A : ',GlobalAverage%sum0/(1.0*nRuns*nEns*massNum)
    write(iFile,'(A)') '================================== '
  end subroutine GlobalAverage_Write

  !****************************************************************************

  subroutine WriteFileNames
    use output

    logical, save :: first = .true.

    if (.not.first) return
    first = .false.

    call writeFileDocu('HiLep.GlobalAverages.dat',&
         &'some global averages and cross sections')

    if (DoOutChannels) then
       call writeFileDocu('HiLep.OutChannels.INIT.dat',&
            &'List of event channels (before interactions and particle decays)')
       call writeFileDocu('HiLep.OutChannels.FINAL.dat',&
            &'List of event channels (after interactions and particle decays)')
    end if

    if (DoTimes) then
       call writeFileDocu('DoTIMES.14n.dat (n=1..3)',&
            &'averaged formation times vs zH (Prod, Form, Form-Prod)')
       call writeFileDocu('DoTIMES.15n.dat (n=0..7)',&
            &'fragmentation type vs zH (cf. eArr)')
       call writeFileDocu('DoTIMES.MP.14n[.lead][.N].dat (n=1..2)',&
            &'')
       call writeFileDocu('DoTIMES.MP.14n[.lead][.N].dat (n=1..2)',&
            &'')

    end if

    if (DoLeptonKinematics) then
       call writeFileDocu('HiLep.lep.XX.kinematics.dat; XX=nu,Q2',&
            &'lepton kinematics: N_e(XX)')
       call writeFileDocu('HiLep.lep.nuAve_Q2.kinematics.dat',&
            &'kinematics (via leptons): <nu>(Q2)')
       call writeFileDocu('HiLep.lep.Q2Ave_nu.kinematics.dat',&
            &'kinematics (via leptons): <Q2>(nu)')


       call writeFileDocu('HiLep.NuQ2planeXS.00n.dat,  n=0..9 (cf.HiLep.NuQ2planeXS.SYS.dat)',&
            &'cross section as function of (nu,Q2)-plane; total and different types')
    end if

    if (DoHadronKinematics) then
       call writeFileDocu('HiLep.XX.YY.kinematics.[no]Acc.dat; XX,YY=nu,Q2,zH',&
            &'kinematics (via hadrons): <YY>(XX) , with/without acceptance')
!!$       call writeFileDocu('HiLep.XX.kinematics.dat; XX=nu,Q2,zH,pT2',&
!!$            &'kinematics (via hadrons): --- not used ---')

       call writeFileDocu('HiLep.XX.idH.[no]Acc.dat; XX=nu,Q2,zH,pT2',&
            &'identified hadrons, dN_id/(N_e dXX), with/without acceptance')

       call writeFileDocu('HiLep.XX.chH.dat; XX=nu,Q2,zH,pT2',&
            &'charged hadrons, dN_ch/(N_e dXX), with & without acceptance')

       call writeFileDocu('HiLep.AvePT2.XX.dat; XX=nu,Q2,zH,pT2',&
            &'identified hadrons, d<pT2>/(N_e dXX), with & without acceptance')


!    call writeFileDocu('fort.56',&
!         &'identified hadrons, dN_id/(N_e dzH), all particles')
!    call writeFileDocu('fort.561',&
!         &'identified hadrons, dN_id/(N_e dzH), all particles, direct production')
!    call writeFileDocu('fort.562',&
!         &'identified hadrons, dN_id/(N_e dzH), all particles, via decay')
!    call writeFileDocu('fort.563',&
!         &'identified hadrons, dN_id/(N_e dzH), all particles, via decay of rho')
!    call writeFileDocu('fort.564',&
!         &'identified hadrons, dN_id/(N_e dzH), all particles, via decay of diffractive rho0')


       call writeFileDocu('HiLep.pT2Ave_Plane.piX.dat; X=+,0,-',&
            &'<pT2>(nu,zH)')
       call writeFileDocu('HiLep.pT2_YY.piX.dat; YY=nu,zH; X=+,0,-',&
            &'d2N_pi/(N_e dYY dpT2)')

       select case (iExperiment)
       case (16) ! EIC
          call writeFileDocu('HiLep.EIC_XX.idH.Acc.dat, XX=nu,Q2,zH',&
               &'identified hadrons, dN/dXX')
          call writeFileDocu('HiLep.EIC_XX.idH.Q2_n.Acc.dat, XX=nu,zH, n=1..5',&
               &'identified hadrons, dN/dXX [bin(Q2)=...]')
       end select

    end if

    ! DoInvMasses = okay. Done in the module

    if (DoFindRho0) then
       call writeFileDocu('HiLep.Rho0_MV.dat','...')
       call writeFileDocu('HiLep.Rho0_MX.dat','...')
       call writeFileDocu('HiLep.Rho0_DE.dat','...')
       call writeFileDocu('HiLep.Rho0_DecTime.dat','...')
       call writeFileDocu('HiLep.Rho0_Mom.dat','...')
       call writeFileDocu('HiLep.Rho0_Theta.dat','...')
       call writeFileDocu('HiLep.OutChannels.RhoOrig.dat','...')
       call writeFileDocu('HiLep.Rho0_*.dat','...')
       ! to be continued ...
    end if

    if (DoClasie) then
       call writeFileDocu('HiLep.Clasie.n.dat; n=1,2','...')
    end if

    if (DoBrooks) then
       call writeFileDocu('HiLep.Brooks.dat','pi+ dN/dpT2, <pT2>; see header in file')
       call writeFileDocu('HiLep.Brooks.detailed.dat','<pT2>; see header in file')
    end if

    if (DoMorrow) then
       call writeFileDocu('HiLep.Morrow.n.dat; n=1..3','...')
       call writeFileDocu('HiLep.MorrowT.n.dat; n=1..2','...')
    end if

    if (DoMandelT) then
       call writeFileDocu('HiLep.MandelstamT.idH.excl.dat',&
            'dN_id/(N_e d|t|), exclusive')
       call writeFileDocu('HiLep.MandelstamT.idH.incl.dat',&
            'dN_id/(N_e d|t|), inclusive')
    end if

    if (DoCentralN) then
       ! to be done
    end if

    if (DoFSIsqrts) then
       call writeFileDocu('HiLep.FSIsqrts.dat',&
            'sqrt(s) of possible FSI')
    end if


  end subroutine WriteFileNames


  !****************************************************************************

  integer function ClassifyFirstEvent(iEns,iNuc)
    use idTable, only: nucleon, pion

    integer, intent(in) :: iEns,iNuc

    type(tPreEvListEntry), POINTER :: pV
    !integer :: nPart, j
    !character*(15), dimension(8) :: AA

    ClassifyFirstEvent = 0

    if (.not.DoClassifyFirst) return

    pV => PreEvArr0(iEns,iNuc)

!!$    do j=1,8
!!$       AA(j) = PartName(pV%preE(j)%ID,pV%preE(j)%charge,pV%preE(j)%antiparticle)
!!$    enddo
!!$    write(6,'(f12.5,8A16)') pV%weight * 1.0,AA


    if (PV%preE(1)%ID == 0) then
    else if (PV%preE(2)%ID == 0) then
    else if (PV%preE(3)%ID == 0) then
       if (PV%preE(1)%ID==nucleon .and. PV%preE(2)%ID==pion) then
          ClassifyFirstEvent=1
       else if (PV%preE(1)%ID==nucleon .and. PV%preE(2)%ID>pion) then
          ClassifyFirstEvent=2
       else if (PV%preE(1)%ID>nucleon .and. PV%preE(2)%ID==pion) then
          ClassifyFirstEvent=3
       else if (PV%preE(1)%ID>nucleon .and. PV%preE(2)%ID>pion) then
          ClassifyFirstEvent=4
       end if
    else if (PV%preE(4)%ID == 0) then
       if (PV%preE(1)%ID==nucleon .and. PV%preE(2)%ID==pion .and. PV%preE(3)%ID==pion) then
          ClassifyFirstEvent=5
       else if (PV%preE(1)%ID==nucleon .and. PV%preE(2)%ID==pion .and. PV%preE(3)%ID>pion) then
          ClassifyFirstEvent=6
       else if (PV%preE(1)%ID>nucleon .and. PV%preE(2)%ID==pion .and. PV%preE(3)%ID==pion) then
          ClassifyFirstEvent=7
       else if (PV%preE(1)%ID>nucleon .and. PV%preE(2)%ID==pion .and. PV%preE(3)%ID>pion) then
          ClassifyFirstEvent=8
       end if
    end if



    return
  end function ClassifyFirstEvent



end module HiLeptonAnalysis
