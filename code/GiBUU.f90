!******************************************************************************
!****p* /GiBUU
! NAME
! program GiBUU
!
! PURPOSE
! This is the main file of the Giessen Boltzmann-Uehling-Uhlenbeck (GiBUU) code.
! It steers the whole simulation, including initialization, time evolution
! (propagation & collisions) and analysis.
!
! COPYRIGHT
! (C) 2005-2014 The GiBUU Team (see full list of authors below)
!
! The GiBUU code is licensed under the GNU General Public License (GPL v2).
! See accompanying LICENSE file or http://www.gnu.org/licenses/gpl-2.0.html.
!
! AUTHOR
! Academic Supervisor:
! * Ulrich Mosel
!
! Alphabetical list of authors:
! * Oliver Buss
! * Thomas Falter
! * Theodoros Gaitanos
! * Kai Gallmeister
! * David Kalok
! * Murat Kaskulov
! * Alexei Larionov
! * Olga Lalakulich
! * Ivan Lappo-Danilewski
! * Tina Leitner
! * Ulrich Mosel
! * Birger Steinmueller
! * Janus Weil
!
! Email:
! * gibuu@projects.hepforge.org
!
! Postal address:
! * Institut fuer Theoretische Physik I
! * Justus-Liebig-Universitaet Giessen
! * Heinrich-Buff-Ring 16
! * D-35392 Giessen
!
! Telephone:
! * +49 (0)641 99 33344   (U. Mosel)
!
! SEE ALSO
! The full physics included in the code is described in this review paper:
! * O. Buss et al., Phys. Rept. 512 (2012) 1,
!   http://inspirehep.net/record/912923
! Additional information and documentation can be found on the website:
! * https://gibuu.hepforge.org
!
! BUGS
! Please report bugs and suggestions for improvements to
! gibuu@projects.hepforge.org.
!******************************************************************************
program GiBUU

  use checks, only: ChecksSetDefaulSwitches
  use collisionNumbering, only: writeCountedEvents
  use inputGeneral, only: eventType, delta_T, &
       num_Energies, num_runs_SameEnergy, current_run_number, readinputGeneral
  use nucleusDefinition, only: tNucleus
  use output, only: header, chapter, subchapter, doPR, &
       timeMeasurement, printTime, writeParticleVector
  use particleDefinition, only: particle, setToDefault, setNumbersToDefault
  use particleProperties, only: initParticleProperties
  use random, only: setRandom
  use statistics, only: splashInfo
  use version, only: printVersion
  use CallStack, only: traceback

  implicit none


  !****************************************************************************
  !****ig* GiBUU/realParticles
  ! SOURCE
  !
  type(particle),Allocatable :: realParticles(:,:)
  !
  ! PURPOSE
  ! The particle vector: real particles
  !
  ! NOTES
  ! * First index : ensemble
  ! * Second index : position within ensemble
  !****************************************************************************

  !****************************************************************************
  !****ig* GiBUU/pertParticles
  ! SOURCE
  !
  type(particle),Allocatable :: pertParticles(:,:)
  !
  ! PURPOSE
  ! The particle vector: perturbative particles
  !
  ! NOTES
  ! * First index : ensemble
  ! * Second index : position within ensemble
  !****************************************************************************

  !****************************************************************************
  !****ig* GiBUU/targetNuc
  ! SOURCE
  !
  type(tNucleus),pointer :: targetNuc => NULL()
  !
  ! PURPOSE
  ! Target nucleus
  !****************************************************************************

  !****************************************************************************
  !****ig* GiBUU/projectileNuc
  ! SOURCE
  !
  type(tNucleus),pointer :: projectileNuc => NULL()   ! Projectile nucleus
  !
  ! PURPOSE
  ! Projectile nucleus
  !****************************************************************************


  integer :: j,k
  logical :: raiseEnergy
  real, save :: delta_T_max   !field to store initial delta_T from inputGeneral


  ! Formats for output:
  character(200), parameter :: format1 = '(79("="),/,5("=")," ",A," ",i5," / ",i5,/,79("="),/)'
  character(200), parameter :: format2 = '(79("#"),/,5("#")," ",A," ",i8," / ",i8,/,79("#"),/)'
  character(200), parameter :: format3 = '(5("=")," ",A," ",i5,"/",i5,"  ",5("#")," ",A," ",i9,"/",i9)'

  call printTime("START")
  write(*,header) 'BUU simulation: start'

  !============= Output the version of the code
  call printVersion

  !============= Initialize the database entries
  write(*,chapter) 'Init database:...'
  call readInputGeneral
  call initParticleProperties
  write(*,chapter) 'Init database: finished'
  call ChecksSetDefaulSwitches(EventType)

  raiseEnergy=.false.
  delta_T_max=delta_T

  energyLoop: do j=1,num_Energies         ! loop over different energies
     write(*,format1) 'Energy loop :',j,num_Energies

     subsequentRunsLoop: do k=1,num_runs_SameEnergy ! loop over subsequent runs
        write(*,format2) 'Run loop :',k,num_runs_SameEnergy

        !============= Output the present status:
        open(123,file='main.run',status='unknown')
        write(123,format3) 'Energy loop :',j,num_Energies, &
             'Run loop :',k,num_runs_SameEnergy
        close(123)

        !============= Take care of random numbers:
        call setRandom()
        current_run_number = current_run_number + 1

        !============= Initialize configuration
        write(*,subchapter) 'Init starting configuration'
        call initConfig(raiseEnergy)

        !============= Do some analysis after the init and before the run.
        call analysis(.false.,.true.)

        !============= Transport the hadronic matter
        write(*,subchapter) 'Do RUN'
        call run

        !============= Final analysis
        write(*,subchapter) 'Analysing output of the run'
        call analysis(k.eq.num_runs_sameEnergy,.false.)
        raiseEnergy=.false.

        !============= Print out statistical informations
        write(*,*) '(FINAL) number of counted Events:'
        call writeCountedEvents(0)

        !============= Sets info module back to normal
        call splashInfo

     end do subsequentRunsLoop
     raiseEnergy=.true.
  end do energyLoop

  call finalCleanup()

  write(*,header) 'BUU simulation: finished'
  call printTime("STOP")


contains


  !****************************************************************************
  !****s* GiBUU/initConfig
  ! NAME
  ! subroutine initConfig(energyRaiseFlag)
  !
  ! PURPOSE
  ! * Initializes starting configuration for particles in each run.
  ! * Here we have some "select case(eventType)":
  !   for each eventType different subroutines take care of the proper
  !   initialization.
  ! * The vectors realParticles and pertParticles are allocated.
  !
  ! INPUTS
  ! * logical :: energyRaiseFlag -- if (true) then we raise the bombarding
  !   energy (used to perform an energy scan)
  !
  !****************************************************************************
  subroutine initConfig(energyRaiseFlag)

    use inputGeneral, only: length_perturbative, length_real, numEnsembles, &
         printParticleVectors

    use initInABox, only: initializeInABox, BoostToEps
    use initInABoxDelta, only: InitInABoxDelta_init
    use initBox, only: initializeBox
    use initPion, only: InitPionInduced
    use initHiPion, only: InitHiPionInduced
    use initLowPhoton, only: initialize_lowPhoton, lowPhotonInit_getRealRun
    use lowElectron, only: init_lowElectron
    use initHiLepton, only: InitHiLeptonInduced
    use initHeavyIon, only: initHeavyIonCollision
    use initElementary, only: initElementaryCollision
    use initHadron, only: initHadronInduced
    use initNeutrino, only: init_Neutrino, neutrinoInit_getRealRun
    use initTransportGivenParticle, only: init_transportGivenParticle
    use initExternal, only: initializeExternal, ExternalIsPerturbative

    use nucleus, only: getTarget, getProjectile
    use densityModule, only: updateDensity, updateRMF, storeFields
    use RMF, only: getRMF_flag
    use yukawa, only: updateYukawa
    use coulomb, only: updateCoulomb
    use energyCalc, only: updateEnergies
    use eventtypes
    use checks, only: ChecksCallEnergy
    use propagation, only: propagate_cascade, updateVelocity
    use insertion, only: GarbageCollection
    use PILCollected, only: PILCollected_ZERO
    use initNucleus_in_PS, only: initNucPhaseSpace, shiftBack
    use groundStateAnalysis, only: countMass
    use baryonPotentialModule, only: HandPotentialToDensityStatic
    use deuterium_PL, only: deuteriumPL_assign


    logical, intent(in) :: energyRaiseFlag
    integer :: lengthReal = 0  ! max number of real particles per ensemble
    integer :: lengthPert = 0  ! max number of pert particles per ensemble
    logical, save :: first = .true.

    if (first) then

      select case (eventType)
      case (elementary, Box)
         if (.not. associated(targetNuc)) then
            allocate(targetNuc) ! dummy target, only %velocity is used
         end if

      case default
         targetNuc => getTarget()    ! set up target resting at 0
         lengthReal = targetNuc%mass ! Real particles fixed by nucleons in target

         if (targetNuc%ReAdjustForConstBinding) then
            write(*,*) 'we now have to readjust the density'
            call handPotentialToDensityStatic(targetNuc)
         end if

      end select

      if (eventType==HeavyIon) projectileNuc => getProjectile()    ! set up projectile resting at 0

      !...Defining maximal lenghts of particle vectors: (Default)

      select case (eventType)
      case (elementary)
         lengthPert =   1
         lengthReal = 100
      case (HeavyIon,ExternalSource)
         lengthPert =   10
         lengthReal = 3000
      case (LoPion)
         lengthPert =  300
      case (RealPhoton)
         if (lowPhotonInit_getRealRun()) then
            lengthReal = targetNuc%mass*11
            lengthPert = 1
         else
            lengthPert = max(100,15*targetNuc%mass)
         end if
      case (LoLepton)
         lengthPert = 10000
      case (Neutrino)
         if (neutrinoInit_getRealRun()) then
            lengthReal = targetNuc%mass*11
            lengthPert = 1
         else
            lengthPert =  max(100,10*targetNuc%mass)
         end if
      case (HiPion)
         lengthPert = 10000
      case (HiLepton)
         lengthPert = 5000
      case (InABox)
         lengthPert =  1
         lengthReal = 5000
      case (InABox_pion,inABox_delta)
         lengthPert =  150
         lengthReal = 5000
      case (Box)
         lengthPert =  1
         lengthReal = 5000
      case (groundState)
         lengthPert =   1
         lengthReal = 500
      case (transportGivenParticle)
         lengthPert = 100
      case (hadron)
         lengthPert =  1
         lengthReal = targetNuc%mass + 100
      end select

      !...Defining maximal lenghts of particle vectors: (From Input)

      if (length_perturbative >= 0) lengthPert = length_perturbative
      if (length_real >= 0)         lengthReal = length_real

      !...Allocate the vectors

      allocate(realparticles(1:numEnsembles,1:lengthReal))
      allocate(pertparticles(1:numEnsembles,1:lengthPert))

      if (doPr(1)) then
         write(*,'(79("#"))')
         write(*,'(79("#"))')
         write(*,'("##              ","#Ensemble","*(","len(Pert)","+","len(Real)",")",'// &
              &  ' "   * ",i3," Bytes                ","##")') sizeof(realParticles(1,1))
         write(*,'("## MEMORY USAGE=",i9,"*(",i9,"+",i9,") = ",i12," <-> ",f6.1," MB ##")') &
              numEnsembles,lengthPert,lengthReal, numEnsembles*(lengthReal+lengthPert), &
              real(numEnsembles*(lengthReal+lengthPert)*sizeof(realParticles(1,1)))/(1024**2)
         write(*,'(79("#"))')
         write(*,'(79("#"))')
      end if

      first=.false.

    else

      call setToDefault(realParticles)
      call setToDefault(pertParticles)

    end if

    call setNumbersToDefault

    call GarbageCollection(pertParticles,.true.)
    call GarbageCollection(realParticles)

    call PILCollected_ZERO ! reset all particle info lists


    !...Do the Init

    select case (eventType)

    case (elementary)

       call initElementaryCollision(realparticles,energyRaiseFlag)

    case (HeavyIon)

       call initHeavyIonCollision(targetNuc,projectileNuc)

       ! Set up test-particles which represent the nuclei:
       call initNucPhaseSpace(realparticles,projectileNuc)
       call initNucPhaseSpace(realparticles,targetNuc)
       ! shift back coordinates along z-axis (see notes in initNucPhaseSpace):
       call shiftBack(realparticles)

       if (.not.getRMF_flag()) then

          call updateDensity(realParticles)
          call updateCoulomb
          call updateYukawa(.true.)

       else

          call updateRMF(realParticles)
          if (doPr(2)) write(*,*) 'Do delta_T propag to get baryon 4-current'
          call propagate_cascade(realParticles,delta_T)
          call storeFields
          call updateRMF(realParticles)
          call updateCoulomb

       end if

    case (LoPion)

       call setUpTarget(.false.)
       call initPionInduced(pertParticles,energyRaiseFlag,targetNuc)

    case (RealPhoton) ! 3= photon nucleus collision

       call setUpTarget()
       call deuteriumPL_assign(realParticles)
       call initialize_lowPhoton(realParticles, pertParticles, &
            energyRaiseFlag, targetNuc)

    case (LoLepton)

       call setUpTarget()
       call init_lowElectron(realParticles, pertParticles, &
            energyRaiseFlag, targetNuc)

    case (Neutrino)

       call setUpTarget()
       call init_Neutrino(realParticles, pertParticles, &
            energyRaiseFlag, num_runs_sameEnergy, targetNuc)

    case (HiPion)

       call setUpTarget(.false.)
       call initHiPionInduced(pertParticles,realParticles,targetNuc)

    case (HiLepton)

       call setUpTarget(.false.)
       call ChecksCallEnergy(-99.9,realParticles)
       call initHiLeptonInduced(realParticles,pertParticles,targetNuc)

    case (InABox,InABox_pion,InABox_delta)

       ! set up the box with nucleons:
       call initializeInABox(realParticles)

       call updateDensity(realParticles)
       call updateCoulomb
       call updateYukawa(.true.)
       call updateVelocity(realParticles)

       ! Set up test-particles which shall be propagated in the box:

       select case (eventType)
       case (inABox)
          call boostToEps(realParticles)
       case (inABox_pion)
          call initPionInduced(pertParticles,energyRaiseFlag)
       case (inABox_delta)
          call initInABoxDelta_init(pertParticles,.true.)
       case default
          call traceback('error inAbox')
       end select

    case (Box)

       call initializeBox(realParticles)

    case (ExternalSource)

       call initializeExternal(realParticles,pertParticles)

       if (ExternalIsPerturbative()) then
          call setUpTarget()
       else

          call updateDensity(realParticles)
          call updateCoulomb
          if (.not.getRMF_flag()) then
             call updateYukawa(.true.)
          else
             call updateRMF(realParticles)
             call storeFields
          end if
       end if

    case (groundState)

       ! Set up test-particles which represent the nucleus:
       call initNucPhaseSpace(realparticles,targetNuc)

       call countMass(realparticles)

       if (.not.getRMF_flag()) then

          call updateDensity(realParticles)
          call updateCoulomb
          call updateYukawa(.true.)

       else

          call updateRMF(realParticles)
          call storeFields
          call updateCoulomb

       end if

    case (hadron)

       call initNucPhaseSpace(realparticles,targetNuc)

       if (getRMF_flag()) then
          call updateRMF(realParticles)
          call updateCoulomb
       end if

       ! Initialise a hadron:
       call initHadronInduced(realParticles,pertParticles)

       if (.not.getRMF_flag()) then

          call updateDensity(realParticles)
          call updateCoulomb
          call updateYukawa(.true.)

       else

          call updateRMF(realParticles)
          call storeFields
          call updateCoulomb

       end if

    case (transportGivenParticle)
       ! Set up test-particles which represent the nucleus:
       call setUpTarget(.false.)

       call init_transportGivenParticle(realParticles,pertParticles)

    case default
       write(*,*) 'No valid eventtype:',eventtype,'STOP'
       call traceback()
    end select

    ! clean up after init routines
    call GarbageCollection(pertParticles,.true.)
    call GarbageCollection(realParticles)

    call deuteriumPL_assign(realParticles)

    if (.not. getRMF_flag()) then
       !... update momentum(0) after updating the mean fields
       call updateEnergies(realParticles)
       call updateEnergies(pertParticles)
       ! Update velocities
       call updateVelocity(pertParticles)
       call updateVelocity(realParticles)
    end if

    if (printParticleVectors) then
       write(*,'(A)') "The initialized particle vectors can be found in the files *_Init_*.dat"
       call writeParticleVector('RealParticles_Init',realParticles)
       call writeParticleVector('PertParticles_Init',pertParticles)
    end if

    ! For safety a final garbage collection:
    call GarbageCollection(pertParticles,.true.)
    call GarbageCollection(realParticles)

    !=== Elementary checks ===

    call ChecksCallEnergy(0.0,realParticles)

  end subroutine initConfig

  !****************************************************************************
  !****s* GiBUU/setUpTarget
  ! NAME
  ! subroutine setUpTarget(doUpdateEnergies)
  !
  ! PURPOSE
  ! Initialize the test particles representing the target nucleus
  !
  ! INPUTS
  ! * logical, OPTIONAL :: doUpdateEnergies -- Flag to indicate, whether
  !   'updateEnergies' should be called
  ! OUTPUT
  ! The particle vector realparticles is modified, also global arrays
  ! connected with this.
  !****************************************************************************
  subroutine setUpTarget(doUpdateEnergies)

    use initNucleus_in_PS, only: initNucPhaseSpace
    use densityModule, only: updateDensity
    use coulomb, only: updateCoulomb
    use yukawa, only: updateYukawa
    use propagation, only: updateVelocity
    use energyCalc, only: updateEnergies

    logical, intent(in), optional :: doUpdateEnergies

    logical :: doUE

    doUE = .true.
    if (present(doUpdateEnergies)) doUE = doUpdateEnergies

    call initNucPhaseSpace(realparticles,targetNuc)
    call updateDensity(realParticles)
    call updateCoulomb
    call updateYukawa(.true.)
    call updateVelocity(realParticles)
    if (doUE) call updateEnergies(realParticles)

  end subroutine setUpTarget


  !****************************************************************************
  !****s* GiBUU/run
  ! NAME
  ! subroutine run
  !
  ! PURPOSE
  ! Propagate a given realParticle and pertParticle vector
  ! according to the BUU equations.
  !****************************************************************************
  subroutine run

    use povray, only: povray_output
    use RMF, only: getRMF_flag
    use yukawa, only: updateYukawa
    use densitymodule, only: updateDensity, gridSpacing, gridSize
    use coulomb, only: updateCoulomb
    use propagation, only: propagate
    use propagation_RMF, only: propagate_RMF
    use checks, only: ChecksCallAll, ChecksCallEnergy, evaluateTimeStep, &
         CheckGridSize !,ChecksCallOccupied
    use collisionTerm, only: collideMain, forceDecays
    use insertion, only: GarbageCollection
    use thermoDynamics, only: updateTemperature
    use energyCalc, only: updateEnergies
    use hadronFormation, only: formation
    use collisionNumbering, only: writeCountedEvents, &
         nullCountedEvents, nulln_participants

    use eventtypes

    use HeavyIonAnalysis, only: doHeavyIonAnalysisTime, HeavyIon_evol
    use sourceAnalysis, only: doSourceAnalysis, getSMM_Flag, resetSMMvalues, &
         stopGiBUU
    use HiLeptonAnalysis, only: HiLeptonAnalysisPerTime
    use HiPionAnalysis, only: HiPionAnalysisPerTime
    use hadronAnalysis, only: doHadronAnalysisTime
    use BoxAnalysis, only: doBoxAnalysisTime
    use transportGivenParticleAnalysis, only: transportGivenParticle_analyze
    use InABoxAnalysis, only: doInABoxAnalysisTime
    use InABoxAnalysisPion, only: InABoxAnalysisPion_count
    use InABoxAnalysisDelta, only: InABoxAnalysisDelta_count
    use Dilepton_Analysis, only: Dilep_Decays
    use FreezeoutAnalysis, only: doFreezeoutAnalysisPerTime
    use EventOutputAnalysis, only: doEventOutput

    use radiativeDeltaDecay, only: doRadiativeDeltaDecay
    use inputGeneral, only: numTimeSteps, printParticleVectorTime, &
         variableTimeStep, time_max, povray_switch, freezeRealParticles, &
         checkGridSize_Flag, doFragmentNucleons
    use deuterium_PL, only: deuteriumPL_assign
    use FragmentNucleons, only: doAddFragmentNucleons
    use analyzeSpectra, only: doAnalyzeSpectra

    integer :: timeStep
    real :: time,delta_T_new

    time=0.

    call nullCountedEvents(0)
    call nulln_participants

    !loop over time steps
    if (numTimeSteps==0) then
       select case (eventtype)
       case (RealPhoton,HiPion,HiLepton)
          call Dilep_Decays(0.,pertParticles,1)
       case (HeavyIon,hadron)
          call Dilep_Decays(0.,realParticles,1)
       end select

       if (eventtype==neutrino.or.eventtype==transportGivenParticle) &
            call doRadiativeDeltaDecay(0,pertParticles)

       call deuteriumPL_assign(realParticles)
    end if

    if (eventtype==Hadron .or. eventType==groundState .or. eventType==ExternalSource) then
       if ( getSMM_Flag() ) then !determine properties of fragmenting source(s)
          call doSourceAnalysis(realParticles,time,delta_T,.false.,targetNuc)
       end if
    end if

    if (printParticleVectorTime) then
       select case (eventType)
       case (HeavyIon,Hadron,groundState)
          call doHeavyIonAnalysisTime(realParticles, time)
       end select
    end if

    ! event output at initialization (before time evolution)
    call doEventOutput(realparticles, pertParticles, 0)

    select case (eventType)
    case (inABox)
      call doInABoxAnalysisTime(realParticles, 0)
    case (inABox_pion)
      call InABoxAnalysisPion_count(pertParticles,time)
    case (inABox_delta)
      call InABoxAnalysisDelta_count(pertParticles,time)
    end select

    TimeStep=0

    PhaseSpaceEvolution : do while(time<time_max-delta_T/2.)

       if (doPR(1)) call timeMeasurement(.true.) ! Reset stopWatch

       if (variableTimeStep) then
          select case (eventtype)
          case (elementary,Heavyion,hadron)
             ! real-real collisions
             call evaluateTimeStep(1,1.2,delta_T_max,time,delta_T_new)
          case default
             ! pert-real collisions
             call evaluateTimeStep(2,1.,delta_T_max,time,delta_T_new)
          end select
          delta_T=delta_T_new
       end if

       TimeStep=TimeStep+1

       !***********************************************************************
       ! Reduce time step in some time-window:
       !if (5. <= time .and. time <= 15.) then
       !  delta_T=delta_T_max/10.
       !else
       !  delta_T=delta_T_max
       !end if
       !***********************************************************************

       time=time+delta_T

       write(*,'(79("*"),/,5("*")," TimeStep ",i4,": ",f7.3," fm"/,79("*"))') timeStep,time

       !=== Do some analysis/statistics on "per timestep" basis:

       if (povray_switch) call povray_output(timestep, realParticles, pertParticles)

       !=== Some "Per-Time-Step - Analysis"

       select case (eventtype)

       case (RealPhoton)
          call Dilep_Decays(time, pertParticles, timeStep/numTimeSteps)

       case (neutrino)
          call doRadiativeDeltaDecay(TimeStep,pertParticles)

       case (HiLepton)
          call HiLeptonAnalysisPerTime(time,pertParticles)
          if ( getSMM_Flag() ) then
             call doSourceAnalysis(realParticles,time,delta_T,.false., &
                  targetNuc)
             ! stop BUU-run after on-set of equilibration:
             if (stopGiBUU) exit PhaseSpaceEvolution
          end if
          call Dilep_Decays(time, pertParticles, timeStep/numTimeSteps)

       case (HiPion)
          call HiPionAnalysisPerTime(timestep,time,pertParticles)

       case (HeavyIon)
          if ( getSMM_Flag() ) then
             call doSourceAnalysis(realParticles,time,delta_T,.false., &
                  targetNuc,projectileNuc)
             ! stop BUU-run after on-set of equilibration:
             if (stopGiBUU) exit PhaseSpaceEvolution
          end if
          call Dilep_Decays(time, realParticles, timeStep/numTimeSteps)

       case (Hadron,groundState,ExternalSource)
          if ( getSMM_Flag() ) then
             call doSourceAnalysis(realParticles,time,delta_T,.false., &
                  targetNuc)
             ! stop BUU-run after on-set of equilibration:
             if (stopGiBUU) exit PhaseSpaceEvolution
          end if
          call Dilep_Decays(time, realParticles, timeStep/numTimeSteps)

       case (transportGivenParticle)
          call transportGivenParticle_analyze(pertParticles,timestep)
          call doRadiativeDeltaDecay(TimeStep,pertParticles)

       case (InABox)
          call doInABoxAnalysisTime(realParticles, timestep)

       case (Box)
          call doBoxAnalysisTime(realParticles, timestep)

       end select

       call doFreezeoutAnalysisPerTime(timestep,time, &
            pertParticles,realParticles)

       ! event output at particular timestep:
       call doEventOutput(realparticles, pertParticles, timestep)

       if (doPR(1)) call timeMeasurement()

       !=== Elementary checks ===

       call ChecksCallAll(timestep,time,realParticles,pertParticles)

       if (doPR(1)) call timeMeasurement()

       !=== Garbage Collection ===
       !         Write(*,*) '**Do Garbage Collection'
       call GarbageCollection(pertParticles, .true.)
       call GarbageCollection(realParticles)

       !=== propagation ===
       write(*,*) '**Do Propagation'
       if (.not. getRMF_flag()) then
         call propagate(realParticles, pertParticles, delta_T, TimeStep)
       else
         call propagate_RMF(realParticles, pertParticles, delta_T, TimeStep)
      end if

       !=== updating mean fields ===
       if (.not.freezeRealParticles) then
           if (.not.getRMF_flag()) then
              call updateDensity(realParticles)
              call updateCoulomb
              call updateYukawa(.false.)
              !=== update momentum(0) after updating the mean fields ===
              call updateEnergies(realParticles)
           end if
           call updateTemperature(realParticles)
       end if

       if (.not.getRMF_flag()) call updateEnergies(pertParticles,.true.)

       if (doPR(1)) call timeMeasurement()

       !=== collisions ===

       write(*,*) '**Do Collisions'
       call formation(pertParticles,realParticles,time,.false.)
       call nullCountedEvents(1)
       call nulln_participants
       call collideMain(pertParticles,realParticles,time)

       call deuteriumPL_assign(realParticles)

       if (doPR(1)) call timeMeasurement()

       !=== checks ===

       call ChecksCallEnergy(time,realParticles)

       ! checks only if particles start to escape out of grid:
       if (checkGridSize_Flag) then
          call checkGridSize(realParticles,time,time_max,eventType, &
               gridSpacing,gridSize)
       end if

       if (doPr(2)) call writeCountedEvents(1,time)
       if (doPr(1)) call writeCountedEvents(2)

       select case (eventtype)
       case (inABox_pion)
          call InABoxAnalysisPion_count(pertParticles,time)
       case (inABox_delta)
          call InABoxAnalysisDelta_count(pertParticles,time)
       case (HeavyIon,groundState,ExternalSource)
          call HeavyIon_evol(realparticles, time, timestep)
       case (hadron)
          call HeavyIon_evol(realparticles, time, timestep)
          call doHadronAnalysisTime(realparticles, time, .false.)
       end select

       call doAnalyzeSpectra(realparticles, pertParticles, timestep,&
            time, .false.)

       !=== Print the particle vector frequently ===

       if (printParticleVectorTime) then

          select case (eventType)
          case (HeavyIon,Hadron,groundState,ExternalSource)
             call doHeavyIonAnalysisTime(realParticles, time)
          end select

       end if

       if (doPR(1)) call timeMeasurement() ! Print stopWatch

    end do PhaseSpaceEvolution

    ! === Forcing decays of all particles at the end ===
    call formation(pertParticles,realParticles,time,.true.)
    call forceDecays(pertParticles,realParticles,time)

    if (.not.getRMF_flag()) then
       ! === update momentum(0) after the forced decays ===
       if (.not.freezeRealParticles) call updateEnergies(realParticles)
       call updateEnergies(pertParticles,.true.)
    end if

    !=== Elementary checks ===

    call ChecksCallEnergy(time,realParticles)
    call ChecksCallAll(timestep,time,realParticles,pertParticles)

    ! === Add nucleons from Fragmentation

    if (doFragmentNucleons) then
       call doAddFragmentNucleons(realParticles,targetNuc)
    end if

    ! === Again "Per-Time-Step"-Analysis after the forced decays

    select case (eventtype)
    case (HeavyIon)
       if (getSMM_Flag()) then
          ! Reset variables needed for sourceAnalysis for the next run
          call doSourceAnalysis(realParticles,time,delta_T,.true., &
               targetNuc,projectileNuc)
          call resetSMMvalues
       end if

    case (Hadron)
       if (getSMM_Flag()) then
          ! Reset variables needed for sourceAnalysis for the next run
          call doSourceAnalysis(realParticles,time,delta_T,.true., &
               targetNuc)
          call resetSMMvalues
       end if
       call doHadronAnalysisTime(realparticles, time, .true.)

    case (ExternalSource)
       if (getSMM_Flag()) then
          ! Reset variables needed for sourceAnalysis for the next run
          call doSourceAnalysis(realParticles,time,delta_T,.true., &
               targetNuc)
          call resetSMMvalues
       end if

    case (neutrino,transportGivenParticle)
       call doRadiativeDeltaDecay(numTimeSteps+1,pertParticles)

    case (HiLepton)
       if (getSMM_Flag()) then
          ! Reset variables needed for sourceAnalysis for the next run
          call doSourceAnalysis(realParticles,time,delta_T,.true., &
               targetNuc)
          call resetSMMvalues
       end if

    end select

    call doFreezeoutAnalysisPerTime(timestep,time,pertParticles,realParticles)

!     call ChecksCallOccupied(realParticles,pertParticles,'At the end of the run: ')
    if (doPR(1)) call timeMeasurement() ! Print stopWatch


  end subroutine run



  !****************************************************************************
  !****s* GiBUU/analysis
  ! NAME
  ! subroutine analysis(finalizeFlag,beforeRUN)
  ! PURPOSE
  ! Analyze output of run.
  ! INPUTS
  ! * logical :: finalizeFlag -- This flag signalises that the last run of
  !   one energy was taking place.
  ! * logical :: beforeRUN    -- Flag to indicate, whether this routine
  !   is called before or after "run" is called. Makes it possible to produce
  !   some analysis-output of the particle vector directly after its init.
  !****************************************************************************
  subroutine analysis(finalizeFlag,beforeRUN)
    use eventtypes
    use inputGeneral, only: FinalCoulombCorrection, printParticleVectors, &
         numTimeSteps
    use CoulombKorrektur, only: CoulombPropagation
    use pionXsection, only: pionXsectionAnalysis
    use HiLeptonAnalysis, only: doHiLeptonAnalysis
    use HiPionAnalysis, only: doHiPionAnalysis
    use InABoxAnalysisPion, only: InABoxAnalysisPion_eval
    use lowphotonanalysis, only: analyze_Photon
    use HeavyIonAnalysis, only: doHeavyIonAnalysis
    use ElementaryAnalysis, only: doElementaryAnalysis
    use lowElectronAnalysis, only: lowElectron_Analyze
    use yScalingAnalysis, only: yScaling_Analyze
    use neutrinoAnalysis, only: neutrino_Analyze
    use initLowPhoton, only: lowPhotonInit_getRealRun
    use initNeutrino, only: neutrinoInit_getRealRun
    use history, only: history_print
    use Dilepton_Analysis, only: Dilep_write_CS
    use radiativeDeltaDecay, only: radiativeDeltaDecay_write_CS
    use EventOutputAnalysis, only: doEventOutput
    use analyzeSpectra, only: doAnalyzeSpectra

    logical, intent(in) :: finalizeFlag, beforeRUN
    integer :: i,j
    logical :: first_historyPrint

    if (beforeRUN) then
       select case (eventType)
       case (HiPion)
          call doHiPionAnalysis(pertParticles,finalizeFlag,beforeRUN)
       case (HiLepton)
          call doHiLeptonAnalysis(realparticles,pertParticles, &
               targetNuc%mass,finalizeFlag,beforeRUN)
       end select

       return ! leave this routine !!!
    end if

    !********** Final Coulomb propagation

    select case (eventType)
    case (LoPion,RealPhoton,neutrino,HiPion,HiLepton)
       if (FinalCoulombCorrection) call CoulombPropagation(pertParticles)
    end select

    !********** Analysis

    select case (eventType)
    case (elementary)
       write(*,*) '   Main : Calling elementary collision analysis routine'
       call doElementaryAnalysis(realparticles, finalizeFlag)
       write(*,*) '   Main : Finished elementary collision analysis routine'
       call doHeavyIonAnalysis(realparticles,pertParticles,finalizeFlag)

    case (HeavyIon,hadron,ExternalSource)
       write(*,*) '   Main : Calling heavy ion induced analysis routine'
       call doHeavyIonAnalysis(realparticles,pertParticles,finalizeFlag)
       call Dilep_write_CS()

    case (LoPion)
       write(*,*) '   Main : Calling pion induced analysis routine'
       call  pionXsectionAnalysis(pertParticles,finalizeFlag)

    case (RealPhoton)
       if (printParticleVectors) then
          write(*,'(A)') "The particle vectors before analysis can be found in the files *_preAnalysis_*.dat"
          call writeParticleVector('RealParticles_preAnalysis',realParticles)
          call writeParticleVector('PertParticles_preAnalysis',pertParticles)
       end if
       call Dilep_write_CS()
       write(*,*) '   Main : Calling photon induced analysis routine'
       if (lowPhotonInit_getRealRun()) then
          call analyze_Photon(realParticles,finalizeFlag)
       else
          call analyze_Photon(pertParticles,finalizeFlag)
       end if

    case (LoLepton)
       call lowElectron_Analyze(pertParticles,finalizeFlag)
       call yScaling_Analyze(pertParticles,finalizeFlag,eventType)

    case (neutrino)
       write(*,*) '   Main : Calling neutrino induced analysis routine'
       if (neutrinoInit_getRealRun()) then
          call neutrino_Analyze(realParticles,finalizeFlag,num_runs_sameEnergy)
       else
          call yScaling_Analyze(pertParticles,finalizeFlag,eventType)
          call neutrino_Analyze(pertParticles,finalizeFlag,num_runs_sameEnergy)
          call radiativeDeltaDecay_write_CS(pertParticles)
       end if

    case (HiPion)
       call Dilep_write_CS()
       call doHiPionAnalysis(pertParticles,finalizeFlag)

    case (HiLepton)
       call doHiLeptonAnalysis(realparticles,pertParticles, &
            targetNuc%mass,finalizeFlag)
       call Dilep_write_CS()

    case (inABox_pion)
       if (finalizeFlag) call InABoxAnalysisPion_eval()

    case (groundState)
       call doHeavyIonAnalysis(realparticles,pertParticles,finalizeFlag)

    case (transportGivenParticle)
       call radiativeDeltaDecay_write_CS(pertParticles)

    end select


    if (PrintParticleVectors) then
       call writeParticleVector('RealParticles_Final',realParticles)
       call writeParticleVector('PertParticles_Final',pertParticles)

       first_historyPrint=.true.
       open(100,file="PertParticles_Final_collHistory.dat")
       do i=1,size(PertParticles,dim=1)
          do j=1,size(PertParticles,dim=2)
             if (Pertparticles(i,j)%ID > 0) then
                call history_print(i,PertParticles(i,j),100,first_historyPrint)
                first_historyPrint=.false.
             end if
          end do
       end do
       close(100)
       first_historyPrint=.true.
       open(100,file="RealParticles_Final_collHistory.dat")
       do i=1,size(realParticles,dim=1)
          do j=1,size(realParticles,dim=2)
             if (realparticles(i,j)%ID > 0) then
                call history_print(i,realParticles(i,j),100)
                first_historyPrint=.false.
             end if
          end do
       end do
       close(100)

    end if

    ! final event output after time evolution
    call doEventOutput(realparticles, pertParticles, numTimeSteps+1)

    call doAnalyzeSpectra(realparticles, pertParticles, numTimeSteps+1,&
         0.0, finalizeFlag)

  end subroutine analysis


  !****************************************************************************
  !****s* GiBUU/finalCleanup
  ! NAME
  ! subroutine finalCleanup
  ! PURPOSE
  ! Clean up all "global" memory. Deallocates particle vectors and
  ! projectile/target nuclei and calls cleanup routines in several modules.
  !****************************************************************************
  subroutine finalCleanup
    use eventtypes
    use baryonWidth, only: cleanupBaryon => cleanUp
    use mesonWidth, only: cleanupMeson => cleanUp
    use densityModule, only: cleanupDensity => cleanup
    use yukawa, only: cleanupYukawa => cleanup
    use thermoDynamics, only: cleanupThermo => cleanup
    use parametrizationsBarMes, only: cleanupBarMes => cleanup
    use initLowPhoton, only: cleanupLowPhoton => cleanup
    use PILCollected, only: PILCollected_DeAllocate
    use selfenergy_baryons, only: cleanupRealParts => cleanup
    use neutrinoAnalysis, only: cleanupNeutrinoAna => cleanup
    use initNeutrino, only: cleanupNeutrinoInit => cleanup
    use baryonWidthMedium_tables, only: cleanupBaryonMedium => cleanup
    use mesonWidthMedium_tables, only: cleanupMesonMedium => cleanup
    use coulomb, only: cleanupCoulomb => cleanup
    use volumeElements, only: cleanupVE => cleanup
    use hadronAnalysis, only: cleanupHadronAna => cleanup

    write(*,*)
    write(*,*) "Cleaning Up All Allocated Memory!"
    write(*,*)
    ! Particle Vectors (real & perturbative):
    if (allocated(realParticles)) deallocate(realParticles, pertParticles)
    ! Target & Projectile Nuclei:
    if (associated(projectileNuc)) deallocate(projectileNuc)
    if (associated(targetNuc)) deallocate(targetNuc)
    ! General stuff:
    ! (Particle Properties, Widths, Cross Sections, Potentials, Densities, ...
    call cleanupBaryon
    call cleanupMeson
    call cleanupDensity
    call cleanupYukawa
    call cleanupThermo
    call cleanupBarMes
    call PILCollected_DeAllocate
    call cleanupRealParts
    call cleanupBaryonMedium
    call cleanupMesonMedium
    call cleanupCoulomb
    call cleanupVE
    ! Inits, Analysis, etc. (specific to different eventTypes):
    select case (eventType)
    case (RealPhoton)
      call cleanupLowPhoton
    case (neutrino)
      call cleanupNeutrinoInit
      call cleanupNeutrinoAna
    case (hadron)
      call cleanupHadronAna
    end select
  end subroutine finalCleanup


end program GiBUU
