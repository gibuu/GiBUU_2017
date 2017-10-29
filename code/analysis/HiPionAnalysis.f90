!******************************************************************************
!****m* /HiPionAnalysis
! NAME
! module HiPionAnalysis
!
! PURPOSE
! This module is for the analysis of high energetic pion-nucleus collisions.
!
! INPUTS
! (none)
!
! NOTES
! This stuff is close to module HiLeptonAnalysis.
!
!******************************************************************************
module HiPionAnalysis

  use histf90
  use hist2Df90
  use histMPf90
  use hist2DMPf90
  use AnaEventDefinition
  use AnaEvent

  implicit none
  private

  !****************************************************************************
  !****g* HiPionAnalysis/Enable
  ! SOURCE
  !
  logical, save :: Enable = .true.
  ! PURPOSE
  ! If .true. the HiPion analysis will be performed, otherwise not.
  !****************************************************************************

  !****************************************************************************
  !****g* HiPionAnalysis/EnablePerTime
  ! SOURCE
  !
  logical, save :: EnablePerTime = .false.
  ! PURPOSE
  ! If .true. the HiPion analysis per timestep will be performed, otherwise
  ! not.
  !****************************************************************************

  !****************************************************************************
  !****g* HiPionAnalysis/DoHarp
  ! SOURCE
  !
  logical, save :: DoHarp = .false.
  ! PURPOSE
  ! switch on/off: Analysis for the HARP experiment
  !****************************************************************************

  !****************************************************************************
  !****g* HiPionAnalysis/DoBlobel
  ! SOURCE
  !
  logical, save :: DoBlobel = .false.
  ! PURPOSE
  ! switch on/off: Analysis according Blobel et al.
  !****************************************************************************

  !****************************************************************************
  !****g* HiPionAnalysis/DoInvMasses
  ! SOURCE
  !
  logical, save :: DoInvMasses   = .false.
  !
  ! PURPOSE
  ! switch on/off: reporting of pairwise-invariant-masses
  !****************************************************************************

  !****************************************************************************
  !****g* HiPionAnalysis/DoOutChannels
  ! SOURCE
  !
  logical, save :: DoOutChannels = .false.
  !
  ! PURPOSE
  ! switch on/off: reporting of all final state channels
  !****************************************************************************

  !****************************************************************************
  !****g* HiPionAnalysis/DoSimpleKin
  ! SOURCE
  !
  logical, save :: DoSimpleKin = .false.
  ! PURPOSE
  ! switch on/off: Analysis for simple kinematics: pZ-, pT-spectra etc.
  !****************************************************************************

  !****************************************************************************
  !****ig* HiPionAnalysis/DoEventAdd
  ! SOURCE
  !
  logical, save :: DoEventAdd   = .true.
  logical, save :: DoEventAdd2  = .true.
  !
  ! PURPOSE
  ! Decide, whether we have to do an EventArr analysis or not. These flags
  ! are not directly accessible, but computed from other values:
  ! * DoEventAdd  = DoOutChannels.or.DoInvMasses
  ! * DoEventAdd2 = DoOutChannels.or.DoInvMasses
  !****************************************************************************

  !************ Histograms for 'DoHarp':

  type(histogramMP),   save :: HMP_Harp
  type(histogram2DMP), save :: H2DMP_Harp
  type(histogram2DMP), save :: H2DMP_HarpHi,H2DMP_HarpHi2,H2DMP_HarpHi3
  type(histogram2DMP), save :: H2DMP_HarpHi2a
  type(histogram2DMP), save :: H2DMP_Harp_CDP

  !************ Histograms for 'DoBlobel':

  type(histogram2DMP), save :: H2DMP_Blobel

  !************ Histograms for 'DoInvMasses':

  ! ---stored in the module---

  !************ Histograms for 'DoSimpleKin':

  type(histogram),save :: hist1
  type(histogram),save :: histPT_pi0,histPT_pip,histPT_pim
  type(histogram2D),save :: H2D_AllPion,H2D_AllOther
  type(histogramMP),save :: HMP_pT
  type(histogramMP),save :: HMP_y2, HMP_y3, HMP_y4, HMP_y5

  !****************************************************************************


  integer, save :: nEventArr = 10000
  type(tAnaEvent), dimension(:), allocatable, save :: EventArr, EventArr0

  logical, save :: FlagReadInput = .true.

  public :: DoHiPionAnalysis
  public :: HiPionAnalysisPerTime

contains


  !****************************************************************************
  !****s* HiPionAnalysis/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "HiPion_Analysis".
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n* HiPionAnalysis/HiPion_Analysis
    ! NAME
    ! NAMELIST /HiPion_Analysis/
    ! PURPOSE
    ! Includes the switches:
    ! * Enable
    ! * EnablePerTime
    ! * DoSimpleKin
    ! * DoHarp
    ! * DoBlobel
    ! * DoInvMasses
    ! * DoOutChannels
    !**************************************************************************

    integer :: ios
    NAMELIST /HiPion_Analysis/ Enable,EnablePerTime, &
         & DoSimpleKin,DoHarp,DoBlobel,DoInvMasses,DoOutChannels

    call Write_ReadingInput('HiPion_Analysis',0)
    rewind(5)
    read(5,nml=HiPion_Analysis,IOSTAT=IOS)
    call Write_ReadingInput('HiPion_Analysis',0,IOS)
    write(*,*) 'Enable HiPion analysis :', Enable
    write(*,*) 'Enable perTime         :', EnablePerTime
    write(*,*)
    write(*,*) 'DoSimpleKin   :',DoSimpleKin
    write(*,*) 'DoHarp        :',DoHarp
    write(*,*) 'DoBlobel      :',DoBlobel
    write(*,*) 'DoOutChannels :',DoOutChannels
    write(*,*) 'DoInvMasses   :',DoInvMasses
    write(*,*)

    DoEventAdd  = DoOutChannels.or.DoInvMasses
    DoEventAdd2 = DoOutChannels.or.DoInvMasses

    write(*,*) 'DoEventAdd    =',DoEventAdd
    write(*,*) 'DoEventAdd2   =',DoEventAdd2

    call Write_ReadingInput('HiPion_Analysis',1)
    FlagReadInput = .false.


  end subroutine readInput



  !****************************************************************************
  !****s* HiPionAnalysis/DoHiPionAnalysis
  ! NAME
  ! subroutine DoHiPionAnalysis(pertPart,finalFlag, beforeRUN)
  !
  ! PURPOSE
  ! This is the routine, which does the actual analysis for HiLepton
  ! calculations
  !
  ! INPUTS
  ! * type(particle),dimension(:,:) :: pertPart -- the perturbative particle vector
  ! * logical                       :: finalFlag -- flag, whether it is the final call
  ! * logical,             optional :: beforeRUN -- flag, whether this routine is
  !   called before all timesteps (i.e. directly after init) or at the end of a run
  !
  ! OUTPUT
  ! (none)
  !****************************************************************************
  subroutine DoHiPionAnalysis(pertPart,finalFlag, beforeRUN)

    use particleDefinition
    use particlePointerListDefinition
    use InitHiPion, only: ConfigParticles, GetBetaFak, getTotalPerweight
    use collisionNumbering, only: pert_firstnumbering12
    use collisionReporter, only: CR_write
    use CollHistory, only: CollHist_ClassifyHist
    use constants, only: pi
    use preEventDefinition
    use PreEvListDefinition
    use preEvList, only: CreateSortedPreEvent, PreEvList_INSERT, PreEvList_Print
    use InvMassesAnalysis, only: InvMasses_INIT, InvMasses_FillEvent, InvMasses_Write
    use output, only: DoPR

    type(particle), intent(in),dimension(:,:), target  :: pertPart
    logical, intent(in)          :: finalFlag
    logical, intent(in),optional :: beforeRUN

    integer, save :: numberRuns=0
    logical, save :: startFlag=.true.
    integer :: iEnsemble, iPart, i,ii,j

    type(particle)          :: Part
    type(particle), pointer :: pPart

    real :: sigmaTot, sigmaAbs
    real, save:: SumSigmaTot, SumSigmaAbs

    integer,Allocatable,save :: nPertPart(:)
    integer,            save :: nPertPartMax = 0

    real :: ptot,pT, yAct, yCM, theta, w,w1, mul0, totalPerWeight

!    real, parameter :: betafak = -7.001922 ! pi-+p, 515 GeV
!      beta = -(P(1,3)+P(2,3))/(P(1,4)+P(2,4)),  fak = log((1+beta)/(1-beta))
!      vgl Pythia, DoPiN.F

    real, save :: betafak
    real, parameter :: yBin1 = -0.75, yBin2= 0.75

    logical :: DoBeforeRUN

    type(tPreEvList), save :: PEList, PEList0
    type(tPreEvListEntry), save :: PreEvListEntry

    type(tParticleListNode),Pointer  :: pNode

    real, dimension(-1:1,-21:21) :: PionOrigin


    if (FlagReadInput) call readInput
    if (.not. enable) return

    DoBeforeRUN = .FALSE.
    if (present(beforeRUN)) DoBeforeRUN = beforeRUN

    write(*,*) '------------------- Perform High Energy Pion Analysis ---------------------'


    ! **********initialization******************************
    if (startflag) then

       write(*,*) '***DoHiPionAnalysis: Initializing...'

       startflag = .FALSE.

       numberRuns = 0
       SumSigmaTot = 0
       SumSigmaAbs = 0

       allocate(nPertPart(size(pertPart,dim=1)))

       if (DoSimpleKin) then

          call CreateHist(hist1,'p_Z, all',-10.,510.,5.)
          call CreateHist(histPT_pi0,'pT, pi0',0.,10.,0.2)
          call CreateHist(histPT_pip,'pT, pip',0.,10.,0.2)
          call CreateHist(histPT_pim,'pT, pim',0.,10.,0.2)


!          call CreateHist2D(H2D_AllPion, 'H2D_AllPion', (/-20.,0./), (/520.,10./), (/20.,0.2/), .true.) ! p_Z, p_T
          call CreateHist2D(H2D_AllPion,  'H2D_AllPion',  (/-10.,0./), (/10.,10./), (/0.2,0.2/), .true.) ! y_L, p_T
          call CreateHist2D(H2D_AllOther, 'H2D_AllOther', (/-10.,0./), (/10.,10./), (/0.2,0.2/), .true.) ! y_L, p_T


          call CreateHistMP(HMP_pT, 'pT', 0.,10.,0.2, 2)

          call CreateHistMP(HMP_y2, 'y (pT=2)', -5.,5.,0.2, 2)
          call CreateHistMP(HMP_y3, 'y (pT=3)', -5.,5.,0.2, 2)
          call CreateHistMP(HMP_y4, 'y (pT=4)', -5.,5.,0.2, 2)
          call CreateHistMP(HMP_y5, 'y (pT=5)', -5.,5.,0.2, 2)
       end if

       if (DoEventAdd) call AllocateEventArr()
       if (DoInvMasses) call InvMasses_INIT

       betafak = GetBetaFak()
       write(*,*) 'betafak=',betafak

       if (DoHarp) then
          call CreateHistMP(HMP_Harp,"p, forward",0.1,1.0,0.05, 1)

          call CreateHist2DMP(H2DMP_Harp,"p vs theta",&
               &(/0.350,0.1/), (/2.150,1.0/), (/0.200,0.05/) , 1, .true.)
          call CreateHist2DMP(H2DMP_HarpHi,"p vs theta, Hi",&
               &(/0.030,0.0/), (/0.300,20.0/), (/0.030,0.500/) , 1, .true.)
          call CreateHist2DMP(H2DMP_HarpHi2,"p vs theta, Hi2",&
               &(/0.000,0.0/), (/0.300,20.0/), (/0.050,0.500/) , 1, .true.)
          call CreateHist2DMP(H2DMP_HarpHi2a,"p vs theta, Hi2",&
               &(/0.000,0.0/), (/0.300,20.0/), (/0.050,0.100/) , 1, .true.)
          call CreateHist2DMP(H2DMP_HarpHi3,"p vs theta, Hi3",&
               &(/0.000,0.0/), (/0.150,20.0/), (/0.025,0.250/) , 1, .true.)

          call CreateHist2DMP(H2DMP_Harp_CDP,"pT vs theta, CDP",&
               &(/0.,0.0/), (/125.,1.5/), (/5.,0.02/) , 1, .true.)
       end if

       if (DoBlobel) then
          call CreateHist2DMP(H2DMP_Blobel, "pT vs y*",&
               &(/-0.05,-4.05/), (/1.4,4.05/), (/0.1,0.1/) , 1, .true.)
       end if

       PionOrigin = 0.0

       write(*,*) '***DoHiPionAnalysis: Initializing... [END]'
    end if

    ! **********Do the actual Run **************************

!    call QYGIVE(' PARP( 91)=')
!    call QYGIVE(' PARJ( 21)=')
!    stop

    if (DoEventAdd) then
       ii = pert_firstnumbering12(DoNotInc=.true.)
       if (ii.gt.nEventArr) then
          do i=1,nEventArr
             call event_CLEAR(EventArr0(i))
             call event_CLEAR(EventArr(i))
          end do
          deallocate(EventArr0)
          deallocate(EventArr)

          ! set new size:
          nEventArr = nint(ii*1.1) ! 10% security buffer

          write(*,*) '...Reallocating EventArr: --> ',nEventArr

          call AllocateEventArr()
       end if
    end if


    if (DoBeforeRUN) then
       if (DoEventAdd) then

          do i=1,nEventArr
             call event_CLEAR(EventArr0(i))
          end do

          do i=1,size(pertPart,dim=1)
             do j=1,size(pertPart,dim=2)
                if (pertPart(i,j)%Id <  0) exit
                if (pertPart(i,j)%Id <= 0) cycle
                if (pertPart(i,j)%firstEvent .eq.0) cycle ! particle did not interact !

                pPart => pertPart(i,j)
                call event_ADD(EventArr0(pPart%firstEvent), pPart)
             end do
          end do

          do i=1,nEventArr
             if (CreateSortedPreEvent(EventArr0(i),PreEvListEntry%preE)) then
                PreEvListEntry%weight = EventArr0(i)%particleList%first%V%perweight
                if (DoOutChannels) call PreEvList_INSERT(PEList0,PreEvListEntry)
             end if
          end do

          if (DoOutChannels) then
             !*****************************************************************
             !****o* HiPionAnalysis/OutChannels.INIT.dat
             ! NAME
             ! file OutChannels.INIT.dat
             ! PURPOSE
             ! This file reports the produced events at initialization time.
             ! It can be enabled by the switch 'DoOutChannels'.
             ! Entries:
             ! * Column #1: Cross section of channel in mb.
             ! * Column #2-X: Final state particles.
             ! NOTES
             ! See also OutChannels.FINAL.dat.
             !*****************************************************************
             open(141,file='OutChannels.INIT.dat', status='unknown')
             rewind(141)
             call PreEvList_Print(141,PEList0,1./(numberRuns+1))
             close(141)
          end if
       end if

       return
    end if


    numberRuns=numberRuns+1
    totalPerWeight = getTotalPerweight()

    SigmaTot = totalPerWeight

    nPertPart = 0

    if (DoEventAdd2) then
       do i=1,nEventArr
          call event_CLEAR(EventArr(i))
       end do
    end if

    do iEnsemble = 1,size(pertPart,dim=1)
       PartLoop:do iPart = 1,size(pertPart,dim=2)
          Part = pertPart(iEnsemble,iPart)
          if (Part%ID <= 0) cycle PartLoop

          w = Part%perweight

          ii = Part%firstEvent
          if (ii.eq.0) then ! particle did not interact !
             SigmaTot = SigmaTot - w
             cycle PartLoop
          end if

          pPart => pertPart(iEnsemble,iPart)
          if (DoEventAdd2) call event_ADD(EventArr(ii), pPart)

          nPertPart(iEnsemble) =  nPertPart(iEnsemble)+1


!          write(*,*) Part%momentum
          ! hadron: pT, y|_Lab
          pT = sqrt(Part%momentum(1)**2+Part%momentum(2)**2)
!          yAct = 0.5 * log((Part%momentum(0)+Part%momentum(3))/(Part%momentum(0)-Part%momentum(3)))
          yAct = rapidity(Part)

!          write(*,*) 'yAct :', yAct,Part%momentum(3),Part%momentum(3)/Part%momentum(0)

!          stop

          theta = 0.0
          if (DoHarp) then

             theta = atan2(pT,Part%momentum(3))
             ptot = absMom(Part)
             if ((theta.gt.0.350).and.(theta.lt.1.55)) &
                  & call AddHistMP(HMP_Harp, Part, ptot,w)

             call AddHist2DMP(H2DMP_Harp, Part,(/theta,ptot/), w)
             call AddHist2DMP(H2DMP_HarpHi, Part,(/theta,ptot/), w)
             call AddHist2DMP(H2DMP_HarpHi2, Part,(/theta,ptot/), w)
             call AddHist2DMP(H2DMP_HarpHi2a, Part,(/theta,ptot/), w)
             call AddHist2DMP(H2DMP_HarpHi3, Part,(/theta,ptot/), w)

             call AddHist2DMP(H2DMP_Harp_CDP, Part,(/theta*180./pi,pT/), w)

!                if (absMom(Part)<0.150) then
             if ((absMom(Part)>5.0).and.((theta.gt.0.030).and.(theta.lt.0.060))) then
!                   call CollHist_WriteHistParticle(6,iEnsemble,iPart)

                ii = max(-21,min(21,CollHist_ClassifyHist(iEnsemble,iPart)))
                PionOrigin(Part%charge,ii) = PionOrigin(Part%charge,ii)+Part%perweight
             end if

          end if


          yCM = yAct + 0.5*betaFak

          w1=0
          if (pT>0) w1=w/pT ! avoid division by zero

          if (DoBlobel) then
             call AddHist2DMP(H2DMP_Blobel, Part, (/pT,yCM/), w1)
          end if

          if (DoSimpleKin) then
             call AddHist(hist1,Part%momentum(3),w)

             if (Part%ID == 101) then
                call AddHist2D(H2D_AllPion, (/yAct, pT/), w)
             else
                call AddHist2D(H2D_AllOther, (/yAct, pT/), w)
             end if

             if (pT<2.0-0.2) then
             else if (pT<2.0+0.2) then
                call AddHistMP(HMP_y2, Part, yCM, w/(2*0.2))
             else if (pT<3.0-0.2) then
             else if (pT<3.0+0.2) then
                call AddHistMP(HMP_y3, Part, yCM, w/(2*0.2))
             else if (pT<4.0-0.2) then
             else if (pT<4.0+0.2) then
                call AddHistMP(HMP_y4, Part, yCM, w/(2*0.2))
             else if (pT<5.0-0.2) then
             else if (pT<5.0+0.2) then
                call AddHistMP(HMP_y5, Part, yCM, w/(2*0.2))
             end if

             if ((yCM.ge.yBin1).and.(yCM.le.yBin2)) then
                w1 = w1/((yBin2-yBin1)*6.28)

                call AddHistMP(HMP_pT, Part, pT, w1)  ! E ds/d^3p !!!

                select case (Part%ID)
                case (101)
                   select case (Part%charge)
                   case (1)
                      call AddHist(histPT_pip, pT, w1)  ! E ds/d^3p !!!
                   case (0)
                      call AddHist(histPT_pi0, pT, w1)  ! E ds/d^3p !!!
                   case (-1)
                      call AddHist(histPT_pim, pT, w1)  ! E ds/d^3p !!!
                   end select
                end select
             end if

          end if



       end do PartLoop
    end do

    SumSigmaTot = SumSigmaTot + SigmaTot

    SigmaAbs = SigmaTot

    if (DoEventAdd2) then
       do i=1,nEventArr
          select case (EventArr(i)%particleList%nEntries)
          case (1)
             if (EventArr(i)%numberParticles(1,ConfigParticles(1)%Charge)==1) then
                pNode => EventArr(i)%particleList%first
                if (associated(pNode)) then
                   SigmaAbs = SigmaAbs - pNode%V%perweight
                else
                   write(*,*) 'OOps, HiPionAnalysis, 339. STOP'
                   stop
                end if
             end if
          case (2)
             if ((EventArr(i)%numberParticles(1,ConfigParticles(1)%Charge)==1).and. &
                  & (SUM(EventArr(i)%numberParticles(7,-2:2))==1)) then
                pNode => EventArr(i)%particleList%first
                if (associated(pNode)) then
                   SigmaAbs = SigmaAbs - pNode%V%perweight
                else
                   write(*,*) 'OOps, HiPionAnalysis, 350. STOP'
                   stop
                end if
             end if

          end select

       end do

    end if

    SumSigmaAbs = SumSigmaAbs + SigmaAbs

    do iEnsemble=1,size(pertPart,dim=1)
      if (nPertPart(iEnsemble)>nPertPartMax) nPertPartMax = nPertPart(iEnsemble)
    end do
    write(*,'(A,i5,A,i5)') 'Number of pert Particles per Ensemble: ', nPertPartMax,' / ',size(pertPart,dim=2)

    mul0 = 1./numberRuns

    write(*,*) 'sigma_tot = ',SumSigmaTot*mul0
    if (DoEventAdd2) write(*,*) 'sigma_abs = ',SumSigmaAbs*mul0

    if (DoSimpleKin) then
       call WriteHist(hist1,101,mul=mul0,add=1e-20,file="KIN.100")

!       call WriteHist_Spline(hist1,101,mul=mul0,add=1e-20,file="KIN.1100")
!       call WriteHist_BSpline(hist1,101,mul=mul0,add=1e-20,file="KIN.1101")

       call WriteHist(histPT_pip,101,mul=mul0,add=1e-20,file="KIN.101")
       call WriteHist(histPT_pi0,101,mul=mul0,add=1e-20,file="KIN.102")
       call WriteHist(histPT_pim,101,mul=mul0,add=1e-20,file="KIN.103")

       call WriteHist2D_Gnuplot(H2D_AllPion, 101, mul=mul0,add=1e-20,file="KIN.310")
       call WriteHist2D_Gnuplot(H2D_AllOther, 101, mul=mul0,add=1e-20,file="KIN.311")

       call WriteHistMP(HMP_pT,101,mul=mul0,add=1e-20,file="KIN.201")
       call WriteHistMP(HMP_y2,101,mul=mul0,add=1e-20,file="KIN.202")
       call WriteHistMP(HMP_y3,101,mul=mul0,add=1e-20,file="KIN.203")
       call WriteHistMP(HMP_y4,101,mul=mul0,add=1e-20,file="KIN.204")
       call WriteHistMP(HMP_y5,101,mul=mul0,add=1e-20,file="KIN.205")
    end if


    if (DoHarp) then
       call WriteHist2DMP_Gnuplot(H2DMP_Harp,  mul=mul0,add=1e-20,&
               & file='HiPion.Harp.xxx.dat',dump=.true.)
       call WriteHist2DMP_Gnuplot(H2DMP_HarpHi,  mul=mul0,add=1e-20,&
               & file='HiPion.Harp.Hi.xxx.dat',dump=.true.)
       call WriteHist2DMP_Gnuplot(H2DMP_HarpHi2,  mul=mul0,add=1e-20,&
               & file='HiPion.Harp.Hi2.xxx.dat',dump=.true.)
       call WriteHist2DMP_Gnuplot(H2DMP_HarpHi2a,  mul=mul0,add=1e-20,&
               & file='HiPion.Harp.Hi2a.xxx.dat',dump=.true.)
       call WriteHist2DMP_Gnuplot(H2DMP_HarpHi3,  mul=mul0,add=1e-20,&
               & file='HiPion.Harp.Hi3.xxx.dat',dump=.true.)

       call WriteHist2DMP_Gnuplot(H2DMP_Harp_CDP,  mul=mul0,add=1e-20,&
               & file='HiPion.Harp.CDP.xxx.dat',dump=.true.)


       call WriteHistMP(HMP_Harp, 101, mul=mul0,add=1e-20,&
               & file='HiPion.Harp.AllForward.dat',dump=.true.)

       open(101,file='HiPion.Harp.PionOrigin.dat',status='unknown')
       rewind(101)
       write(101,'(i3,1P,3e13.5,"  ",3e13.5)') 0,PionOrigin(-1:1,0)*mul0, 0.0,0.0,0.0
       do i=1,21
          write(101,'(i3,1P,3e13.5,"  ",3e13.5)') i,&
               & PionOrigin(-1:1,i)*mul0, PionOrigin(-1:1,-i)*mul0
       end do
       write(101,*)
       write(101,*)
       write(101,'("   ",1P,3e13.5,"  ",3e13.5)') sum(PionOrigin,dim=2)*mul0
       write(101,*)
       write(101,*)
       write(101,'(i3,2P,3f13.5,"  ",3f13.5)') 0,PionOrigin(-1:1,0)/sum(PionOrigin,dim=2), 0.0,0.0,0.0
       do i=1,21
          write(101,'(i3,2P,3f13.5,"  ",3f13.5)') i,&
               & PionOrigin(-1:1,i)/sum(PionOrigin,dim=2), PionOrigin(-1:1,-i)/sum(PionOrigin,dim=2)
       end do
       close(101)
    end if

    if (DoBlobel) then
       call WriteHist2DMP_Gnuplot(H2DMP_Blobel, 101,  mul=mul0,add=1e-20,&
            & file='HiPion.Blobel.xxx.dat',dump=.true.)
    end if

    if (DoInvMasses) then
       do i=1,nEventArr
          call InvMasses_FillEvent(EventArr(i))
       end do

       call InvMasses_Write(mul=mul0)
    end if

    if (DoEventAdd) then
       do i=1,nEventArr
          if (CreateSortedPreEvent(EventArr(i),PreEvListEntry%preE)) then
             PreEvListEntry%weight = EventArr(i)%particleList%first%V%perweight
             if (DoOutChannels) call PreEvList_INSERT(PEList,PreEvListEntry)
          end if
       end do

       if (DoOutChannels) then
          !********************************************************************
          !****o* HiPionAnalysis/OutChannels.FINAL.dat
          ! NAME
          ! file OutChannels.FINAL.dat
          ! PURPOSE
          ! This file reports the produced events at the end of the simulation.
          ! It can be enabled by the switch 'DoOutChannels'.
          ! Entries:
          ! * Column #1: Cross section of channel in mb.
          ! * Column #2-X: Final state particles.
          ! NOTES
          ! See also OutChannels.INIT.dat.
          !********************************************************************
          open(141,file='OutChannels.FINAL.dat', status='unknown')
          rewind(141)
          call PreEvList_Print(141,PEList,mul0)
          close(141)
       end if
    end if


    call cR_Write(numberRuns)


    if (DoPR(0)) write(*,*) 'Histograms written.'


    if (finalFlag) then
!!$          do i=1,Size(pertPart,dim=2)
!!$             if (pertPart(1,i)%ID.gt.0) then
!!$                if  (pertPart(1,i)%FirstEvent.eq.0) then
!!$                   write(*,*) 'xxx: ',pertPart(1,i)%number
!!$                endif
!!$             endif
!!$          enddo
    end if


!   contains
    !**************************************************************************

!     subroutine calcMomenta
!       use constants, only: mN, mPi
!
!       double precision srtS, theta, phi, beta(3)
!
!       write(*,*) '~~~~~'
!       write(*,*) '~~~~~ [] calcMomenta; start...'
!
!
!       call MP_WriteVersion(6)
!       call MP_Set3(1, mPi, 0d0,0d0,515d0) ! pion
!       call MP_Set3(2,  mN, 0d0,0d0,  0d0) ! Nucleon
!       call MP_Copy(1,3)
!       call MP_Copy(2,4)
!
!       call MP_CalcROBO(1,2, srts, theta, phi, beta)
!       !         write(*,*) '~~~~~ [] calcMomenta; ...'
!       !         write(*,*) '~~~~~                 sqrts, theta, phi, beta=',sqrts, theta, phi, beta
!
!       call MP_ROBO_inv(3,4, theta, phi, beta(1),beta(2),beta(3))
!
!       call MP_Write(6,1,4)
!
!       write(*,*) '~~~~~ [] calcMomenta; finished!'
!       write(*,*) '~~~~~'
!
!     end subroutine calcMomenta

    !**************************************************************************


  end subroutine DoHiPionAnalysis


  !****************************************************************************
  !****s* HiPionAnalysis/HiPionAnalysisPerTime
  ! NAME
  ! subroutine HiPionAnalysisPerTime(iTime, Time, pParts)
  ! PURPOSE
  ! produce statistical output after every time step
  !
  ! INPUTS
  ! * integer        :: iTime       -- number of timestep
  ! * real           :: Time        -- time of time step (in fm)
  ! * type(particle) :: pParts(:,:) -- particle vector of perturbative parts
  !
  ! OUTPUT
  ! internal variables are changed and some output to predefined files.
  !
  !****************************************************************************
  subroutine HiPionAnalysisPerTime(iTime, Time, pParts)
    use particleDefinition
    use hist2Df90
    use InitHiPion, only: GetBetaFak
    use Dilepton_Analysis, only: Dilep_Decays
    use inputGeneral, only: numTimeSteps

    integer, intent(in)       :: iTime
    real,    intent(in)       :: Time
    type(particle),intent(in) :: pParts(:,:)

    logical, save :: DoInit = .TRUE.

    type(histogram2D), save :: hist_pZ_Time, hist_pT_Time

    integer :: i,j
    real :: pZ,pT, w, ws, yAct,yCM
    !real :: weight,nu,Q2,epsilon
    !integer :: EventType

    real, save :: betafak
    real, parameter :: yBin1 = -0.5, yBin2= 0.5

    if (DoInit) then

       betafak = GetBetaFak()
       write(*,*) 'betafak=',betafak

       call CreateHist2D(hist_pZ_Time, "XS vs. pZ and TIME",&
            & (/-50.0,0.0/), (/550.0,30.0/), (/10.0,0.5/), .true.)
       call CreateHist2D(hist_pT_Time , "XS vs. pT and TIME",&
            & (/0.0,0.0/), (/15.0,30.0/), (/0.2,0.5/), .true. )



       DoInit = .FALSE.
    end if

    if (EnablePerTime) then

       do j=1,Size(pParts,dim=2)
          do i=1,Size(pParts,dim=1)
             if (pParts(i,j)%ID > 0) then

                w  = pParts(i,j)%perWeight
                ws = w * pParts(i,j)%ScaleCS

                pZ = pParts(i,j)%momentum(3)

                call AddHist2D(hist_pZ_Time, (/pZ,time/), ws, w)

                pT = sqrt(pParts(i,j)%momentum(1)**2+pParts(i,j)%momentum(2)**2)

                yAct = 0.5 * log((pParts(i,j)%momentum(0)+pParts(i,j)%momentum(3)) &
                     & / (pParts(i,j)%momentum(0)-pParts(i,j)%momentum(3)))

                yCM = yAct + 0.5*betaFak

                if (yCM.lt.yBin1) cycle
                if (yCM.gt.yBin2) cycle


                call AddHist2D(hist_pT_Time, (/pT,time/), ws, w)
             end if
          end do
       end do

       call WriteHist2D_Gnuplot(hist_pZ_Time, file="HiPionAnalysisPerTime_pZ.dat")
       call WriteHist2D_Gnuplot(hist_pT_Time, file="HiPionAnalysisPerTime_pT.dat")
    end if

    call Dilep_Decays(Time,pParts,iTime/numTimeSteps)   ! do Dilepton analysis

  end subroutine HiPionAnalysisPerTime

  !****************************************************************************


!   subroutine WriteFileNames
!     use output, only: writeFileDocu
!
!     logical, save :: first = .true.
!
!     if (.not.first) return
!     first = .false.
!
!     if (DoOutChannels) then
!        call writeFileDocu('HiPion.OutChannels.INIT.dat',&
!             &'List of event channels (before interactions and particle decays)')
!        call writeFileDocu('HiPion.OutChannels.FINAL.dat',&
!             &'List of event channels (after interactions and particle decays)')
!     end if
!
!     if (DoHarp) then
!        call writeFileDocu('HiPion.Harp.xxx.dat',&
!             &'---to be explained---')
!        call writeFileDocu('HiPion.Harp.Hi.xxx.dat',&
!             &'---to be explained---')
!        call writeFileDocu('HiPion.Harp.Hi2.xxx.dat',&
!             &'---to be explained---')
!        call writeFileDocu('HiPion.Harp.Hi3.xxx.dat',&
!             &'---to be explained---')
!
!        call writeFileDocu('HiPion.Harp.AllForward.dat',&
!             &'---to be explained---')
!
!        call writeFileDocu('HiPion.Harp.PionOrigin.dat',&
!             &'this is not a histogram ---to be explained---')
!     end if
!
!     if (DoBlobel) then
!        call writeFileDocu('HiPion.Blobel.xxx.dat',&
!             &'---to be explained---')
!     end if
!
!     if (DoSimpleKin) then
!        call writeFileDocu('KIN.100','p_Z for all particles')
!
!        call writeFileDocu('KIN.10[123]','p_T for pi+,pi0,pim')
!        call writeFileDocu('KIN.310','(yL,pT) for all Pions')
!        call writeFileDocu('KIN.311','(yL,pT) for all other but Pions')
!
!        call writeFileDocu('KIN.201','p_T for all particles')
!        call writeFileDocu('KIN.20[2-5]','y with pT=2,3,4,5 for all particles')
!
!     end if
!
!
!
!   end subroutine WriteFileNames

  !****************************************************************************

  subroutine AllocateEventArr()

    use particlePointerList, only: ParticleList_INIT

    integer :: i

    allocate(EventArr(nEventArr))
    allocate(EventArr0(nEventArr))
    do i=1,nEventArr
       call ParticleList_INIT(EventArr0(i)%particleList)
       call event_CLEAR(EventArr0(i))
       call ParticleList_INIT(EventArr(i)%particleList)
       call event_CLEAR(EventArr(i))
    end do
  end subroutine AllocateEventArr

end module HiPionAnalysis
