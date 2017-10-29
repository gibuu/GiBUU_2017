!******************************************************************************
!****m* /inputGeneral
! NAME
! module inputGeneral
! PURPOSE
! Collects some general parameters into the namelist 'input'
! and provides a routine for reading the values from a JobCard.
!******************************************************************************
module inputGeneral

  implicit none


  !****************************************************************************
  !****g* inputGeneral/povray_switch
  ! PURPOSE
  ! Switch for generating Povray-Output
  ! SOURCE
  !
  logical        ,save :: povray_switch=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/freezeRealParticles
  ! PURPOSE
  ! Switch for not propagating real particles
  ! SOURCE
  !
  logical        ,save :: freezeRealParticles=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/FinalCoulombCorrection
  ! PURPOSE
  ! Switch for Coulomb correction at the end of each run of the
  ! outgoing particles
  ! SOURCE
  !
  logical        ,save :: FinalCoulombCorrection=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/length_perturbative
  ! PURPOSE
  ! Length of perturbative particle vector (per ensemble). If negative,
  ! it will be determined automatically by event type.
  ! SOURCE
  !
  integer, save :: length_perturbative = -1
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/length_real
  ! PURPOSE
  ! Length of real particle vector (per ensemble). If negative,
  ! it will be determined automatically by event type.
  ! SOURCE
  !
  integer, save :: length_real = -1
  !****************************************************************************


  !****************************************************************************
  !****g* inputGeneral/eventtype
  ! PURPOSE
  ! Switch for the type of event
  !
  ! possible values: see module eventtypes
  !
  ! SOURCE
  !
  integer        ,save :: eventtype=3
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/fullEnsemble
  ! PURPOSE
  ! Switch for the type of simulation:
  ! * .false.=parallel ensembles
  ! * .true.=full ensemble
  ! See also "localEnsemble".
  ! SOURCE
  !
  logical, save :: fullEnsemble = .false.
  !****************************************************************************


  !****************************************************************************
  !****g* inputGeneral/localEnsemble
  ! PURPOSE
  ! Switch for the type of simulation:
  ! * .false. = parallel or full ensembles (depending on the value of the fullEnsemble switch).
  ! * .true. = fullEnsemble with "local collisionCriteria", see Lang/Babovsky et al., J. Comput. Phys. 106 (1993) 391-396.
  ! Setting localEnsemble = .true. will implicitly set fullEnsemble = .true. (disregarding its value in the jobcard).
  ! SOURCE
  !
  logical,save :: localEnsemble = .false.
  !****************************************************************************


  !****************************************************************************
  !****g* inputGeneral/numEnsembles
  ! PURPOSE
  ! Number of parallel ensembles
  ! SOURCE
  !
  integer        ,save :: numEnsembles=300
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/numTimeSteps
  ! PURPOSE
  ! Number of time steps
  ! SOURCE
  !
  integer        ,save :: numTimeSteps=100
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/delta_T
  ! PURPOSE
  ! time difference for time stepping
  ! SOURCE
  !
  real        ,save :: delta_T=0.2
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/variableTimeStep
  ! PURPOSE
  ! Switch for using of variable time step:
  ! * .false.= use constant time step delta_T (see above).
  ! * .true.= use time step computed from the frequency of collisions.
  !   In this case the input delta_T is used as the maximum
  !   allowed time step.
  ! SOURCE
  !
  logical        ,save :: variableTimeStep=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/checkGridSize_Flag
  ! PURPOSE
  ! Switch for checking if particles escape out of grid.
  ! possible values:
  ! * .false.= no check.
  ! * .true. = check is performed, and a warning flag is printed out,
  !   in case that particles are outside of the grid.
  ! * check valid only for real particles.
  ! SOURCE
  !
  logical        ,save :: checkGridSize_Flag=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/time_max
  ! PURPOSE
  ! Maximum time until which the time evolution will be computed
  ! in the case of variableTimeStep = .true.
  ! SOURCE
  !
  real        ,save :: time_max=30.
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/num_energies
  ! PURPOSE
  ! Number of different energies for energy scans
  ! SOURCE
  !
  integer, save :: num_energies=1
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/num_runs_sameEnergy
  ! PURPOSE
  ! Number of runs with the same energy in the initialization.
  ! SOURCE
  !
  integer, save :: num_runs_sameEnergy=1
  !****************************************************************************


  !****************************************************************************
  !****g* inputGeneral/current_run_number
  ! PURPOSE
  ! Counter to keep track of the number of the current run.
  ! This is incremented for each new run, no matter if the
  ! energy was increased or not.
  ! SOURCE
  !
  integer, save :: current_run_number = 0
  !****************************************************************************


  !****************************************************************************
  !****g* inputGeneral/printParticleVectors
  ! PURPOSE
  ! Switch to turn on the printing of the particle vector at
  ! the start and end of a run.
  ! SOURCE
  !
  logical, save :: printParticleVectors=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/printParticleVectorTime
  ! PURPOSE
  ! * Switch to turn on the printing of the particle vector
  !   as function of time.
  ! * Useful for event classes using real particles (HeavyIon,Hadron).
  ! * See also 'timeForOutput' and 'timeSequence'.
  ! SOURCE
  !
  logical, save :: printParticleVectorTime = .false.
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/printParticleVectorsFormat
  ! PURPOSE
  ! Select the format for printing the particle vectors.
  ! Possible values are:
  ! * 1: ASCII
  ! * 2: binary
  ! SOURCE
  !
  integer, save :: printParticleVectorsFormat = 1
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/timeForOutput
  ! PURPOSE
  ! * Time (fm/c) after which the particle vector
  !   is printed during run (see also variable "timeSequence").
  ! * valid only if printParticleVectorTime = .true.
  ! SOURCE
  !
  Real, save :: timeForOutput=50.
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/timeSequence
  ! PURPOSE
  ! * Time sequence (fm/c) of time dependent printing of the particle vector
  ! * valid only if printParticleVectorTime = .true.
  ! SOURCE
  !
  Real, save :: timeSequence=10.
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/path_To_Input
  ! PURPOSE
  ! Path to input files. This switch needs to be set to the local path of the
  ! 'buuinput' directory, which contains various input files for GiBUU.
  ! SOURCE
  !
  character(300), save :: path_To_Input = '../../buuinput'
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/continousBoundaries
  ! SOURCE
  !
  logical, save :: continousBoundaries=.false.
  !
  ! PURPOSE
  ! * Switch to turn on continous boundary conditions.
  ! * Implications for density and propagation.
  ! * This means that particles are propagated according to continous
  !   boundaries. A particle leaving the grid will move back in from the
  !   opposite side. The densities are carefully constructed such that places
  !   at the opposite side contribute to places on the near side.
  ! * What is still missing is the full implementation in collision criteria,
  !   this is not done yet for the two body collisions!
  !   Be careful therefore with the 2-Body-collisions at the edges.
  !   A particle at one edge does not see its scattering partner at the
  !   opposite edge.
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/LRF_equals_CALC_frame
  ! PURPOSE
  ! * Switch to turn on the assumption that calculation frame and LRF frame coincide
  ! * Only useful for reactions close to ground state !!!
  ! SOURCE
  !
  logical, save :: LRF_equals_CALC_frame=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* inputGeneral/DoFragmentNucleons
  ! PURPOSE
  ! * Switch to turn on/off adding of nucleons stemming from fragmentation
  !   of bound clusters.
  ! SOURCE
  !
  logical, save :: DoFragmentNucleons=.false.
  !****************************************************************************


contains


  !****************************************************************************
  !****s* inputGeneral/readInputGeneral
  ! NAME
  ! subroutine readInputGeneral
  ! PURPOSE
  ! Reads input in jobcard out of namelist "input".
  !****************************************************************************
  subroutine readInputGeneral
    use output, only: Write_ReadingInput,doPrLevel,DoPrLevelDefault, &
         setPrintParticleVectorsFormat
    use eventTypes, only: HeavyIon, HeavyIon
    use CallStack, only: Traceback

    use nucleus, only: getTarget
    use nucleusdefinition
    type(tnucleus), pointer :: targetNuc
    real ::Atarget

    !**************************************************************************
    !****n* inputGeneral/input
    ! NAME
    ! NAMELIST /input/
    ! PURPOSE
    ! Includes the input switches:
    ! * path_To_Input
    ! * numEnsembles
    ! * eventtype
    ! * fullEnsemble
    ! * localEnsemble
    ! * delta_T
    ! * numTimeSteps
    ! * variableTimeStep
    ! * time_max
    ! * num_energies
    ! * num_runs_sameEnergy
    ! * checkGridSize_Flag
    ! * continousBoundaries
    ! * FinalCoulombCorrection
    ! * length_perturbative
    ! * length_real
    ! * freezeRealParticles
    ! * printParticleVectors
    ! * printParticleVectorTime
    ! * printParticleVectorsFormat
    ! * timeForOutput
    ! * timeSequence
    ! * DoPrLevel
    ! * povray_switch
    ! * LRF_equals_CALC_frame
    ! * DoFragmentNucleons
    !**************************************************************************
    NAMELIST /input/ path_To_Input, numEnsembles, eventtype, fullEnsemble, localEnsemble, &
                     delta_T, numTimeSteps, variableTimeStep, time_max, &
                     num_energies, num_runs_sameEnergy, &
                     checkGridSize_Flag, continousBoundaries, FinalCoulombCorrection, &
                     length_perturbative, length_real, freezeRealParticles, &
                     printParticleVectors, printParticleVectorTime, timeForOutput, timeSequence, &
                     printParticleVectorsFormat, &
                     DoPrLevel, povray_switch, LRF_equals_CALC_frame, DoFragmentNucleons

    integer :: i,ios
    character(*), parameter :: format1 = '(" Simulation type is ",A)'
    logical :: ex

    call Write_ReadingInput("input",0)
    rewind(5,IOSTAT=ios)
    if (ios /= 0) then
      write(*,*)
      write(*,*) "Rewind error: Maybe you forgot to provide a jobcard?"
      write(*,*) "The input jobcard file needs to be piped into the GiBUU executable like this:"
      write(*,*) "./GiBUU.x < jobcard"
      write(*,*) "Furthermore, you can redirect the output of GiBUU into a log file if desired:"
      write(*,*) "./GiBUU.x < jobcard > log.txt"
      write(*,*)
      stop
    end if
    read(5,nml=input,IOSTAT=ios)
    call Write_ReadingInput("input",0,ios)
    if (ios /= 0) then
       call Traceback('Error in  namelist "input" : This namelist is crucial. STOP!')
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set numTimeSteps = 0 for target mass A = 1  to avoid crash for numTimeSteps>0
    TargetNuc =>  getTarget()
    Atarget = targetNuc%mass
    if (Atarget == 1) numTimeSteps = 0
    write(*,*) 'number of time steps set to 0 since target = H'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (localEnsemble) then
       write(*,format1) '"LOCAL ENSEMBLES" (="FULL ENSEMBLE" + local coll.crit.)'
       write(*,*)
       if (.not.fullensemble) then
          write(*,'(A)') ' WARNING!!! setting fullEnsemble=.true. (localEnsemble=.true. implies this)'
          write(*,*)
          fullEnsemble=.true.
       end if

    else
       if (fullEnsemble) then
          write(*,format1) '"FULL ENSEMBLE" (Kodama coll.crit.)'
       else
          write(*,format1) '"PARALLEL ENSEMBLES"'
       end if
   end if

    write(*,*) 'Eventtype    =',eventtype
    write(*,*) 'Number of ensembles   =',NumEnsembles
    write(*,*) 'Use variable time step:', variableTimeStep
    if (variableTimeStep) then
      write(*,*) '   Maximum time:', time_max
    else
      time_max = float(numTimeSteps)*delta_T
    end if
    write(*,*) 'Time step size=',delta_T
    write(*,*) 'Check grid dimensions :',checkGridSize_Flag
    write(*,*) 'Number time steps/run =',numTimeSteps
    write(*,*) 'Number run/energy     =',num_runs_sameEnergy
    write(*,*) 'Number of energies    =',num_energies
    write(*,*)
    write(*,*) 'Print particle vector at start/end of a run: ',printParticleVectors
    write(*,*) 'Print particle vector at different times   : ',printParticleVectorTime
    if (printParticleVectorTime) then
       if (TimeForOutput < 0. .or. TimeForOutput > time_max) then
          write(*,'(A,F7.2,A)') '   Not an appropriate starting time for output!',TimeForOutput,' [fm/c]'
          write(*,*) '  Set starting time to its default value (50 fm/c)'
          TimeForOutput = 50.
       end if
       if (TimeSequence < delta_T .or. TimeSequence > time_max) then
          write(*,'(A,F7.2,A)') '   Not appropriate value for timeSequence!',TimeSequence,' [fm/c]'
          write(*,*) '   Set timeSequence to its default value (10 fm/c)'
          TimeSequence = 10.
       end if
       write(*,'(A,F7.2,A)') '   time for output: ',TimeForOutput,' [fm/c]'
       write(*,'(A,F7.2,A)') '   time sequence  : ',timeSequence,' [fm/c]'
    end if
    if (printParticleVectors.or.printParticleVectorTime) then
       write(*,*) 'Print particle vector format: ',printParticleVectorsFormat
    end if
    call setPrintParticleVectorsFormat(printParticleVectorsFormat)

    write(*,*) 'Use continous boundary conditions: ',continousBoundaries
    write(*,*) 'Final Coulomb correction at end of run=',FinalCoulombCorrection

    call ExpandPath(path_to_input)

    write(*,'(2A)') ' Path to input files= ',trim(path_to_Input)

    ! Check whether the directory exists, note: inquire does not work on directories!
    Inquire(file=trim(path_to_Input)//'/baryonWidthVacuum.dat.bz2',exist=ex)
    if (.not. ex) then
       write(*,*) 'Directory does not exist: ', path_to_input
       call Traceback('Stop!')
    end if

    ! Try to open dummy file to check for write access
    open(333,file=trim(path_to_Input)//'/directory.test',iostat=ios)
    if (ios /= 0) write(*,*) 'WARNING: No rights to write to', path_to_input
    close(333)

    if (freezeRealParticles) write(*,*) 'WARNING: REAL PARTICLES ARE FROZEN!!!!!'

    if (length_real >= 0) then
       write(*,*) 'Length of real         vector=', length_real
    else
       write(*,*) 'Length of real         vector= (hard wired in the code)'
    end if
    if (length_perturbative >= 0) then
       write(*,*) 'Length of perturbative vector=', length_perturbative
    else
       write(*,*) 'Length of perturbative vector= (hard wired in the code)'
    end if

    write(*,*)
    do i=-10,5
       if (DoPrLevel(i).neqv.DoPrLevelDefault(i)) then
          write(*,*) 'DoPr(',i,') changed: ',DoPrLevelDefault(i), &
          '-->',DoPrLevel(i)
       end if
    end do
    write(*,*)
    write(*,*) 'Generate povray output: ',povray_switch
    write(*,*)

    if (LRF_equals_CALC_frame) then
       write(*,'(A)') ' WARNING: In the following we will assume that calculation frame = LRF frame!!!'
       if (eventtype==HeavyIon) then
          write(*,'(A)') ' STOP: Serious error. It does not make sense to use this assumption for heavy ion runs!!'
          call Traceback('Use LRF_equals_CALC_frame=.false. !!!')
       end if
    end if

    write(*,*)
    if (DoFragmentNucleons) then
       write(*,'(A)') ' Nucleons of fragmenting sources: added at the end'
    else
       write(*,'(A)') ' Nucleons of fragmenting sources: none'
    end if
    write(*,*)

    call Write_ReadingInput("input",1)

  end subroutine readInputGeneral



  subroutine ExpandPath(p)
    character(len=*),intent(inout) :: p
    character(len=100) :: home
    integer :: stat

    if (p(1:1)/='~') return

    call get_environment_variable(name="HOME",value=home,status=stat)
    if (stat/=0) then
      write(*,*) "Environment variable $HOME not set!"
      return
    end if

    p = trim(home) // p(2:)

  end subroutine


end module inputGeneral
