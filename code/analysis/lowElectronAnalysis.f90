!******************************************************************************
!****m* /lowElectronAnalysis
! NAME
! lowElectronAnalysis
! PURPOSE
! Module with all routines necessary to make output for low energy electron
! induced processes.
!******************************************************************************
module lowElectronAnalysis

  implicit none
  private

  Public:: lowElectron_Analyze


  !****************************************************************************
  !****g* lowElectronAnalysis/dE_switch
  ! SOURCE
  !
  logical, save :: dE_switch=.false.
  ! PURPOSE
  ! If .true. then also dSigma/dE is produced, if false not..
  !****************************************************************************


  !****************************************************************************
  !****g* lowElectronAnalysis/dOmega_switch
  ! SOURCE
  !
  logical, save :: dOmega_switch=.false.
  ! PURPOSE
  ! If .true. then also dSigma/dOmega is produced, if false not..
  !****************************************************************************


  logical, save :: initflag=.true.

contains

  !****************************************************************************
  !****s* lowElectronAnalysis/readinput
  ! NAME
  ! subroutine readinput
  ! INPUTS
  ! NONE
  ! OUTPUT
  ! NONE
  ! PURPOSE
  ! This subroutine reads out the jobcard "electronAnalysis".
  ! Only called once to initialize the module.
  !****************************************************************************
  subroutine readinput
    use output, only: Write_ReadingInput

    integer :: IOS

    !**************************************************************************
    !****n* lowElectronAnalysis/lowElePhoto_Analysis
    ! NAME
    ! NAMELIST /lowElePhoto_Analysis/
    ! PURPOSE
    ! This namelist for module "lowElectronAnalysis" includes:
    ! * dOmega_switch
    ! * dE_switch
    !**************************************************************************
    NAMELIST /lowElePhoto_Analysis/ dOmega_switch,dE_switch

    call Write_ReadingInput('lowElePhoto_Analysis',0)
    rewind(5)
    read(5,nml=lowElePhoto_Analysis,IOSTAT=IOS)
    write(*,*) 'dSigma/dOmega output?', dOmega_switch
    write(*,*) 'dSigma/dE     output?', dE_switch
    call Write_ReadingInput('lowElePhoto_Analysis',0,IOS)

  end subroutine readinput


  !****************************************************************************
  !****s* lowElectronAnalysis/lowElectron_Analyze
  ! NAME
  ! subroutine lowElectronAnalysisAnalysis(particles,finalFlag)
  ! INPUTS
  ! * type(particle), intent(in),dimension(:,:)  :: particles        ! Particles which shall be analyzed
  ! * logical, intent(in) :: finalFlag                               ! if .true. than the final output
  !   for a series of calls will be done
  ! USES
  ! * AnaEvent
  ! NOTES
  ! * This subroutine produces output for electron-nucleus scattering
  !****************************************************************************
  subroutine lowElectron_Analyze  (particles,finalFlag)
    use particleDefinition
    use lowElectron, only: le_get_FirstEventRange,le_get_Energy_li,le_get_Energy_lf,writeOriginXS
    use Electron_origin, only: lE_isResEvent
    use eventTypes, only: RealPhoton, LoLepton
    use histf90
    use hist2Df90
    use degRad_conversion, only: radian
    use output, only: intToChar
    use initLowPhoton, only: energy_gamma
    use AnaEventDefinition
    use AnaEvent
    use ZeroPionAnalysis
    use inputGeneral, only: input_eventType => eventtype

    type(particle), intent(in),dimension(:,:) ,target :: particles
    logical, intent(in) :: finalFlag

    ! Local variables:
    integer, dimension (1:2) :: firstEvents
    type(tAnaEvent), Allocatable, dimension(:) :: events ! A list of all events
    type(tAnaEvent), Allocatable, dimension(:) :: events_res ! A list of resonance induced events
    type(tAnaEvent), Allocatable, dimension(:) :: events_bg ! A list of background induced events
    type(particle)         , POINTER :: particlePointer
    integer :: i,j,first

    real, parameter :: ekinMax=2. ! Maximal kinetic energy for dsigma/dEkin
    real, parameter :: dEkin=0.01 ! Delta(eKin) for dsigma/dEKin
    real  :: dPhi   ! Delta(phi) for dsigma/dOmega
    real  :: dTheta ! Delta(theta) for dsigma/dOmega

    integer,save :: numberOfCalls=0
    integer,save :: numberOfFinals=0
    integer,save :: numberOfCallsFinals=0

    logical, dimension(1:numStableParts) :: printflags

    ! In these hists we save the information which the "AnaEvent" subroutines are returning us
    real , dimension(1:dimSigma,1:2),save :: sigma,sigma_res,sigma_bg,check
    real , dimension(1:dimSigma,1:2),save :: sigma_0pions,sigma_0pions_res,sigma_0pions_bg


    ! used for histograms for each energy of the initial lepton
    type(histogram),save,  &
         & dimension(1:numStableParts,-2:2) :: dE_hists, dE_hists_1X, dE_hists_2X, dE_hists_MULTI
!         & dimension(lBound(particleIDs,dim=1):uBound(particleIDs,dim=1),-2:2) :: dE_hists

    ! used for histograms sumed over energies of the initial lepton
    type(histogram),save,  &
         & dimension(1:numStableParts,-2:2) :: dE_hists_ALL, dE_hists_1X_ALL, dE_hists_2X_ALL, dE_hists_MULTI_ALL


    type(histogram),save,  &
         & dimension(1:numStableParts,-2:2) :: dTheta_hists, dPhi_hists
!         & dimension(lBound(particleIDs,dim=1):uBound(particleIDs,dim=1),-2:2) :: dTheta_hists,dPhi_hists

    type(histogram2D),save, &
         & dimension(1:numStableParts,-2:2) :: dOmega_hists
!         & dimension(lBound(particleIDs,dim=1):uBound(particleIDs,dim=1),-2:2) :: dOmega_hists

    real :: energyLabel
    integer :: check_index

    dPhi=radian(10.)   ! Delta(phi) for dsigma/dOmega
    dTheta=radian(10.) ! Delta(theta) for dsigma/dOmega


    if (initflag) then
       call readinput
       printFlags=.false.
       printFlags(1)=.true. ! nucleon
       printFlags(2)=.true. ! eta
       printFlags(7)=.true. ! pion
       call set_particleIDs_flag(printFlags)
       initflag=.false.
       numberOfCalls=0
       numberOfFinals=0
       numberOfCallsFinals=0

       !***********************************************************************
       !****o* lowElectronAnalysis/lowPhotoEle_sigma.dat
       ! NAME
       ! file lowPhotoEle_sigma.dat
       !
       ! PURPOSE
       ! The file is produced in the runs with eventtype=3=RealPhoton  and  eventtype=4=LoLepton.
       !
       ! Cross sections (in microbarns) for electron or photon induced events for  preselected  final states
       !
       ! Columns:
       ! * #1: variable which was raised
       !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode, Elepton for nuXsectionMode=2=dSigmadQsdElepton  and so on)
       ! * #2-#99:  xsec for preselected states
       !   see description in AnaEvent.f90, subroutine event_sigma  OR in the output file sigma.dat
       !   In each channel the outgoing lepton is presupposed (unless explicitely stated otherwise)
       !
       !   The description of some columns is also given in documentation to neutrino_total_Xsection_multiplicities.dat
      !************************************************************************
       open(111,file='lowPhotoEle_sigma.dat')
       write(111,*) '#'
       write(111,*) '# Columns:'
       write(111,*) '# photon energy [GeV], Xsections (same order as in sigma.dat) [microbarn/GeV/A]'
       write(111,'(A)',advance='no') '#     001'
       do i=2,dimSigma
          write(111,'(A15)',advance='no') intToChar(i)
       end do
       write(111,*)
       close(111)
       open(111,file='lowPhotoEle_sigma_0pions.dat')
       write(111,*) '#'
       write(111,*) '# Columns:'
       write(111,*) '# photon energy [GeV], Xsections (same order as in sigma_0pions.dat) [microbarn/GeV/A]'
       write(111,'(A)',advance='no') '#     001'
       do i=2,dimSigma
          write(111,'(A15)',advance='no') intToChar(i)
       end do
       write(111,*)
       close(111)

       !***********************************************************************
       !****o* lowElectronAnalysis/lowPhotoEle_sigma.res.dat
       ! NAME
       ! file lowPhotoEle_sigma.res.dat
       ! PURPOSE
       ! Cross sections for electron or photon induced events.
       ! The same as lowPhotoEle_sigma.dat, but
       ! Only such Events are included, where the initial scattering event was a resonance excitation.
       !***********************************************************************
       open(111,file='lowPhotoEle_sigma.res.dat')
       write(111,*) '#'
       write(111,*) '# Columns:'
       write(111,*) '# photon energy [GeV], Xsections (same order as in sigma.dat)  [microbarn/GeV/A]'
       close(111)

       !***********************************************************************
       !****o* lowElectronAnalysis/lowPhotoEle_sigma.bg.dat
       ! NAME
       ! file lowPhotoEle_sigma.bg.dat
       ! PURPOSE
       ! Cross sections for electron or photon induced events.
       ! The same as lowPhotoEle_sigma.dat, but
       ! Only such Events  are included, where the initial scattering event was a background event.
       !***********************************************************************
       open(111,file='lowPhotoEle_sigma.bg.dat')
       write(111,*) '#'
       write(111,*) '# Columns:'
       write(111,*) '# photon energy [GeV], Xsections (same order as in sigma.dat)  [microbarn/GeV/A]'
       close(111)

    end if

    numberOfCalls=numberOfCalls+1
    numberOfCallsFinals=numberOfCallsFinals+1


    write(*,'(A,2I9)') '** In lowElectron_Analyze    numberofCalls   numberofFinals', numberofCalls,numberofFinals

    ! (1) Setting up the particles into the events
    ! This is done with the help of %firstIndex:
    ! Particles stemming from the same event get in the init the same %firstIndex entry.
    ! During the run %firstinit stays constant and is inherited during collisions.

    firstEvents=le_get_FirstEventRange()
    write(*,*) firstevents
    allocate(Events(firstEvents(1):firstEvents(2)))
    allocate(Events_res(firstEvents(1):firstEvents(2)))
    allocate(Events_bg(firstEvents(1):firstEvents(2)))
    do i=firstEvents(1),firstEvents(2)
       call event_init(events(i))
       call event_init(events_res(i))
       call event_init(events_bg(i))
    end do


    do i=lbound(particles,dim=1),ubound(particles,dim=1)
       do j=lbound(particles,dim=2),ubound(particles,dim=2)
          if (particles(i,j)%ID.le.0) cycle
          first=particles(i,j)%firstEvent
          if (first.le.0) cycle ! If this analysis routine is used with realParticles as input (e.g. called by LowPhotonAnalysis), then there are spectator nucleons with first=0, which should not be added to the event.
!          write(*,*) first
          particlePointer=>particles(i,j)
          ! Add particle to the event with its firstEvent index.
          if (lE_isResEvent(first)) then
             call event_add(events_res(first),particlePointer)
          else
             call event_add(events_bg(first),particlePointer)
          end if
          call event_add(events(first),particlePointer)
       end do
    end do

    ! (2) Use the list "events" to evaluate cross sections

    select case (input_eventtype)
    case (LoLepton)
       ! electron induced events
       energyLabel= le_get_energy_li()-le_get_energy_lf()
    case (RealPhoton)
       ! photon induced events
       energyLabel= energy_gamma
    end select


    ! if(numberOfCalls.eq.1) then we initialize the histograms
    if (dE_switch) call event_dSigma_dE(events,0.,ekinMax,dEkin,'diff_'//trim(intToChar(numberofFinals)), &
         & numberOfCalls, &
         & dE_hists,(numberOfCalls.eq.1),sameFileNameIn=.true., &
         & histsMulti=dE_hists_Multi,hists1X=dE_hists_1X,hists2X=dE_hists_2X)

    if (dOmega_switch) call event_dSigma_dOmega(events,dTheta,dPhi,'diff_'//trim(intToChar(numberofFinals)), &
         & numberOfCalls, &
         & dTheta_hists, dPhi_hists, dOmega_hists,(numberOfCalls.eq.1),sameFileNameIn=.true.)


    ! sumed over all energies of the incoming lepton
    if (dE_switch) call event_dSigma_dE(events,0.,ekinMax,dEkin,'diff_ALL', &
         & numberOfCallsFinals, &
         & dE_hists_ALL,(numberOfCalls.eq.1  .and. numberofFinals.eq.0),sameFileNameIn=.true., &
         & histsMulti=dE_hists_Multi_ALL, hists1X=dE_hists_1X_ALL, hists2X=dE_hists_2X_ALL)



    call event_sigma(events,     sigma,     (numberOfCalls.eq.1), numberOfCalls, identifier=energyLabel)
    call event_sigma(events_res, sigma_res, (numberOfCalls.eq.1), numberOfCalls, identifier=energyLabel)
    call event_sigma(events_bg,  sigma_bg,  (numberOfCalls.eq.1), numberOfCalls, identifier=energyLabel)

    ! final states with 0 pions and 2 pions  (title 0pions for historical reasons)
    call event_sigma_0pions(events,     sigma_0pions,    (numberOfCalls.eq.1), numberOfCalls, identifier=energyLabel)
    call event_sigma_0pions(events_res, sigma_0pions_res,(numberOfCalls.eq.1), numberOfCalls, identifier=energyLabel)
    call event_sigma_0pions(events_bg,  sigma_0pions_bg, (numberOfCalls.eq.1), numberOfCalls, identifier=energyLabel)




    check=sigma_bg+sigma_res-sigma
    do check_index=lbound(check,dim=1)+1,ubound(check,dim=1)
       if (abs(check(check_index,1)).gt.0.0001) then
          write(*,*) 'error in lowElectronAnalysis'
          write(*,*) check_index, check(check_index,1),sigma_bg(check_index,1),sigma_res(check_index,1),-sigma(check_index,1)
          write(*,*) numberofCalls
          stop
       end if
    end do


    if (finalFlag) then
       open(111,file='lowPhotoEle_sigma.dat',position='append')
       if (numberOfCalls.gt.1) then
          write(111,'(1P,300E15.5)') energyLabel,sigma(2:,1)/float(numberOfCalls), &
               & sqrt(max(0.,(sigma(:,2)-sigma(:,1)**2/float(numberOfCalls))/float(numberOfCalls-1)/float(numberOfCalls)))
       else
          write(111,'(1P,300E15.5)') energyLabel,sigma(2:,1)
       end if
       close(111)
       ! Resonance contribution
       open(111,file='lowPhotoEle_sigma.res.dat',position='append')
       if (numberOfCalls.gt.1) then
          write(111,'(1P,300E15.5)') energyLabel,sigma_res(2:,1)/float(numberOfCalls) &
               & ,   sqrt(max(0.,(sigma_res(:,2)-sigma_res(:,1)**2/float(numberOfCalls)) &
               & /float(numberOfCalls-1)/float(numberOfCalls)))
       else
          write(111,'(1P,300E15.5)') energyLabel,sigma_res(2:,1)
       end if
       close(111)
       ! Background contribution
       open(111,file='lowPhotoEle_sigma.bg.dat',position='append')
       if (numberOfCalls.gt.1) then
          write(111,'(1P,300E15.5)') energyLabel,sigma_bg(2:,1)/float(numberOfCalls) &
               & ,   sqrt(max(0.,(sigma_bg(:,2)-sigma_bg(:,1)**2/float(numberOfCalls)) &
               & /float(numberOfCalls-1)/float(numberOfCalls)))
       else
          write(111,'(1P,300E15.5)') energyLabel,sigma_bg(2:,1)
       end if
       close(111)
       !!
       ! The same for final states with 0 and 2 pions
       open(111,file='lowPhotoEle_sigma_0pions.dat',position='append')
       if (numberOfCalls.gt.1) then
          write(111,'(1P,300E15.5)') energyLabel,sigma_0pions(2:,1)/float(numberOfCalls) &
               & ,   sqrt(max(0.,(sigma_0pions(:,2)-sigma_0pions(:,1)**2/float(numberOfCalls))&
               & /float(numberOfCalls-1)/float(numberOfCalls)))
       else
          write(111,'(1P,300E15.5)') energyLabel,sigma_0pions(2:,1)
       end if
       close(111)
       ! Resonance contribution
       open(111,file='lowPhotoEle_sigma_0pions.res.dat',position='append')
       if (numberOfCalls.gt.1) then
          write(111,'(1P,300E15.5)') energyLabel,sigma_0pions_res(2:,1)/float(numberOfCalls) &
               & ,   sqrt(max(0.,(sigma_0pions_res(:,2)-sigma_0pions_res(:,1)**2/float(numberOfCalls)) &
               & /float(numberOfCalls-1)/float(numberOfCalls)))
       else
          write(111,'(1P,300E15.5)') energyLabel,sigma_0pions_res(2:,1)
       end if
       close(111)
       ! Background contribution
       open(111,file='lowPhotoEle_sigma_0pions.bg.dat',position='append')
       if (numberOfCalls.gt.1) then
          write(111,'(1P,300E15.5)') energyLabel,sigma_0pions_bg(2:,1)/float(numberOfCalls) &
               & ,   sqrt(max(0.,(sigma_0pions_bg(:,2)-sigma_0pions_bg(:,1)**2/float(numberOfCalls)) &
               & /float(numberOfCalls-1)/float(numberOfCalls)))
       else
          write(111,'(1P,300E15.5)') energyLabel,sigma_0pions_bg(2:,1)
       end if
       close(111)

       call writeOriginXS( 1.0/float(numberOfCalls) )

    end if

    ! (3) Clear the list "events" to clear the memory
    do i=firstEvents(1),firstEvents(2)
       call event_clear(events(i))
       call event_clear(events_res(i))
       call event_clear(events_bg(i))
    end do
    deallocate(events)
    deallocate(events_res)
    deallocate(events_bg)

    if (finalFlag) then
       numberOfCalls=0
       numberOfFinals= numberOfFinals+1
    end if

  end subroutine lowElectron_Analyze
end module lowElectronAnalysis
