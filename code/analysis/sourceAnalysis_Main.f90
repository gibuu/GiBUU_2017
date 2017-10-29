!******************************************************************************
!****m* /sourceAnalysis
! NAME
! module sourceAnalysis
!
! PURPOSE
! The main module which determines parameters for (afterwards) statistical
! fragmentation and stopping GiBUU-run.
! NOTES
! * The major files which are neccessary for the afterwards statistical
!   fragmentation have the structure "Source<realToChar(time)>fmc.dat", and
!   they store the major information event-by-event.
! * The file "SourceEvol.dat" provides with additional information on time
!   evolution of thermodynamical properties at the center of the sources,
!   mass, charge numbers and total energy.
!******************************************************************************
module sourceAnalysis

  private

  !****************************************************************************
  !****n* sourceAnalysis/SMM_input
  ! NAME
  ! NAMELIST /SMM_input/
  ! PURPOSE
  ! Includes input switches:
  ! * SMM_Flag
  ! * rho_cutoff
  ! * spectator_cutoff
  ! * A_cutoff
  ! * SelectionMethod
  ! * betaChoice
  ! * MaxTimePrinting
  ! * DetailedHyperonOutput
  ! * hyperSource
  !****************************************************************************

  !****************************************************************************
  !****g* sourceAnalysis/SMM_Flag
  ! SOURCE
  !
  logical, SAVE :: SMM_Flag         = .false.
  ! PURPOSE
  ! if .true. then source analysis is switched on
  !****************************************************************************

  !****************************************************************************
  !****g* sourceAnalysis/rho_cutoff
  ! SOURCE
  !
  real,    SAVE :: rho_cutoff       = 100.
  ! PURPOSE
  ! density cutoff (in units of the saturation density "rhoNull")
  ! which defines "emitting" particles
  !****************************************************************************

  !****************************************************************************
  !****g* sourceAnalysis/spectator_cutoff
  ! SOURCE
  !
  real,    SAVE :: spectator_cutoff = 1.
  ! PURPOSE
  ! min. value of number of collisions which defines
  ! "spectator"-matter
  !****************************************************************************


  !****************************************************************************
  !****g* sourceAnalysis/A_cutoff
  ! SOURCE
  !
  integer,    SAVE :: A_cutoff = 2
  ! PURPOSE
  ! min. value of the source mass number
  !****************************************************************************


  !****************************************************************************
  !****g* sourceAnalysis/SelectionMethod
  ! SOURCE
  !
  integer, SAVE :: SelectionMethod  = 0
  ! PURPOSE
  ! defines the selection method of spectators and fireball.
  ! Can be used in high energy Hadron-Nucleus events.
  !****************************************************************************

  !****************************************************************************
  !****g* sourceAnalysis/betaChoice
  ! SOURCE
  !
  integer, SAVE ::  betaChoice = 0
  ! PURPOSE
  ! Defines the way to calculate the source velocity in RMF mode.
  ! Has no influence in calculations with Skyrme potential.
  !****************************************************************************

  !****************************************************************************
  !****g* sourceAnalysis/DetailedHyperonOutput
  ! SOURCE
  !
  logical, SAVE :: DetailedHyperonOutput = .true.
  ! PURPOSE
  ! print more informations for Hyperons and pions.
  !****************************************************************************


  !****************************************************************************
  !****g* sourceAnalysis/hyperSource
  ! SOURCE
  !
  logical, SAVE :: hyperSource = .false.
  ! PURPOSE
  ! If true, the Lambda and Sigma0 hyperons will be included
  ! into source
  !****************************************************************************


  !****************************************************************************
  !****g* sourceAnalysis/stopGiBUU
  ! SOURCE
  !
  logical, SAVE :: stopGiBUU = .false.
  ! PURPOSE
  ! Flag which indicates when to stop the GiBUU run.
  !****************************************************************************

  !****************************************************************************
  !****g* sourceAnalysis/iprint
  ! SOURCE
  !
  Integer, SAVE :: iprint = 0
  ! PURPOSE
  ! Indicates how many times one saves the results. When the results
  ! have been printed MaxTimePrinting-times, the variable "StopGiBUU" is set
  ! to true and the actual GiBUU run stops.
  !****************************************************************************

  !****************************************************************************
  !****g* sourceAnalysis/MaxTimePrinting
  ! SOURCE
  !
  Integer, SAVE :: MaxTimePrinting = 10
  ! PURPOSE
  ! Indicates how many times the results are printed into files.
  ! NOTES
  ! Set MaxTimePrinting to a very big value, i.e. 1000, if you wish that
  ! the BUU-run developes until time=time_max.
  !****************************************************************************

  !****************************************************************************
  !****g* sourceAnalysis/isut
  ! SOURCE
  !
  Integer, SAVE :: isut = 0
  ! PURPOSE
  ! Actual number of the GiBUU run.
  !****************************************************************************

  !****************************************************************************
  !****g* sourceAnalysis/time_out
  ! SOURCE
  !
  real,    save :: time_out=0.0
  ! PURPOSE
  ! Printing control.
  !****************************************************************************

  !****************************************************************************
  !****g* sourceAnalysis/itSteps
  ! SOURCE
  !
  integer, save :: itSteps=0
  ! PURPOSE
  ! Printing control.
  !****************************************************************************


  !****************************************************************************
  !****g* sourceAnalysis/stossParameter
  ! SOURCE
  !
  Real, SAVE :: stossParameter=1000.
  ! PURPOSE
  ! * Impact parameter [fm] of the reaction (HeavyIon or Hadron eventTypes).
  !   The impact parameter is taken from the initialization modules
  !   "initHeavyIon" or "initHadron" and copied to the local variable
  !   "stossParameter".
  ! * The (unrealistic) default value is neccessary for control that the
  !   impact parameter has been extracted correctly.
  !****************************************************************************

  public :: DoSourceAnalysis,getSMM_Flag,resetSMMvalues
  public :: stopGiBUU

contains

  !****************************************************************************
  !****s* sourceAnalysis/DoSourceAnalysis
  ! NAME
  ! subroutine DoSourceAnalysis
  !
  ! PURPOSE
  ! The main routine which decides when to STOP BUU and switch to
  ! Statistical Multifragmentation.
  ! INPUTS
  ! * type(particle),dimension(:,:) :: realPV -- real particle vector
  ! * real                          :: time   -- actual time of simulation
  ! * type(tNucleus)                :: targetNuc -- target properties
  ! * type(tNucleus),optional       :: projectileNuc -- projectile properties
  ! * logical                       :: FinalFlag -- print source(s) info
  !   after forced decays (important for high energy runs)
  !
  ! NOTES
  ! * This analysis is valid only for real particles.
  ! * Determination of fragmenting source(s) as function of time
  !   (controlled by the variable "timeSequence").
  ! * Determination of physical properties at the center of the source(s),
  !   e.g. pressure components, density, degree of equilibration.
  ! * GiBUU-run does not immediatly stop after onset of equilibration.
  !   First source(s) info is printed out at several times, and then
  !   GiBUU-run is terminated.
  ! * For high energy runs analysis routine is called again after
  !   the forced decays.
  ! * Statistical fragmentation is performed afterwards using an extra
  !   program (see workingCode/testRun/auswerteTools/clusters/smm/smm_Main.f90).
  !****************************************************************************
  subroutine DoSourceAnalysis(realPV,time,delta_T,FinalFlag,targetNuc,projectileNuc)

    use particleDefinition
    use nucleusDefinition
    use initHeavyIon, only: b_HI => b
    use initHadron, only: b_had => b
    use RMF, only: getRMF_flag
    use inputGeneral, only: timeForOutput,timeSequence
    use determineSource, only: Get_FragmentingSource,deallocate_source,&
         & Get_InitialPosX
    use sourceProperties

    implicit none

    !-----------------------------------------------------------------------
    ! Input variables
    !-----------------------------------------------------------------------
    type(particle),dimension(:,:),intent(in) :: realPV ! real particleVector
    real,                         intent(in) :: time,delta_T
    type(tNucleus),               pointer    :: targetNuc
    type(tNucleus),optional,      pointer    :: projectileNuc
    logical,                      intent(in) :: FinalFlag
    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    integer, dimension(Size(realPV,dim=1),Size(realPV,dim=2)) :: sourceType
    integer, dimension(1:Size(realPV,dim=1))                  :: NumSources
    integer :: numEnsemples,numParticles

    logical, save :: impact_Flag=.true.
    logical, save :: FirstTimeCall=.true.
    !-----------------------------------------------------------------------
    ! Set the impact parameter of the actual run. Unfortunatly, the
    ! impact parameter between heavyIon and Hadron eventTypes is defined
    ! differently. For this reason we have two different variables from
    ! the corresponding init-modules (see the use declarations) for the
    ! same physical quantity.
    !-----------------------------------------------------------------------
    if (impact_Flag) then
       if (present(projectileNuc)) then
          stossParameter = b_HI
       else
          stossParameter = b_had
       end if
       impact_Flag = .false.
    end if
    if (stossParameter==1000.) then
       write(*,*) 'Module sourceAnalysis_Main, routine DoSourceAnalysis:'
       write(*,*) 'Wrong input for the impact parameter: ',stossParameter
       write(*,*) '!!! Termination of the program !!!'
       STOP
    end if
    !-----------------------------------------------------------------------
    numEnsemples = Size(realPV,dim=1)
    numParticles = Size(realPV,dim=2)


    if (present(projectileNuc) .and. SelectionMethod==2) then
       if (FirstTimeCall) then
          call Get_InitialPosX(NumEnsemples,NumParticles,realPV)
          FirstTimeCall = .false.
       end if
    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (FinalFlag) then

       ! determine again source(s) after forced decays:
       if ( present(projectileNuc) ) then

          call Get_FragmentingSource(numEnsemples,numParticles,&
               &  SelectionMethod,betaChoice,hyperSource, &
               &  realPV,Spectator_cutoff,rho_cutoff,A_cutoff,stossParameter, &
               &  Numsources,sourceType, &
               &  targetNuc,projectileNuc)
       else

          call Get_FragmentingSource(numEnsemples,numParticles,&
               &  SelectionMethod,betaChoice,hyperSource, &
               &  realPV,Spectator_cutoff,rho_cutoff,A_cutoff,stossParameter, &
               &  Numsources,sourceType,targetNuc)
       end if

       if ( getRMF_Flag() ) then
          call sourceProperties_Main(time,NumEnsemples,NumSources,realPV)
       end if

       ! overwrite output files after forced decays:
       call WriteSourceInfo(time,stossParameter,isut,NumSources,hyperSource,FinalFlag)
       call printParticleVector(isut,numEnsemples,numParticles,&
            &                   time,sourceType,realPV)

       ! deallocate fields:
       call deallocate_source
       if ( getRMF_Flag() ) call deallocate_sourceFields

       RETURN

    end if
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (time.ge.timeForOutput) then

       Time_out = Time_out + delta_T
       itSteps  = itSteps  + 1

       if (itSteps==1 .or. abs(Time_out-timeSequence)<1.e-6) then

          Time_out = 0.0

          !-----------------------------------------------------------------
          !(1) determine source(s)
          !-----------------------------------------------------------------
          if ( present(projectileNuc) ) then

             call Get_FragmentingSource(numEnsemples,numParticles,&
                  &  SelectionMethod,betaChoice,hyperSource, &
                  &  realPV,Spectator_cutoff,rho_cutoff,A_cutoff,stossParameter, &
                  &  Numsources,sourceType, &
                  &  targetNuc,projectileNuc)
          else

             call Get_FragmentingSource(numEnsemples,numParticles,&
                  &  SelectionMethod,betaChoice,hyperSource, &
                  &  realPV,Spectator_cutoff,rho_cutoff,A_cutoff,stossParameter, &
                  &  Numsources,sourceType,targetNuc)
          end if

          !-----------------------------------------------------------------
          ! (2) determine local quantities of the source(s):
          ! * calculation of pressure and density fields on the grid.
          ! * determination of the pressure & density fields
          !   at the center of the source(s).
          ! * Determination of the average equilibration.
          ! * Determination of internal energy of the source(s).
          ! * Printing of source info as function of time.
          !-----------------------------------------------------------------
          if ( getRMF_Flag() ) then
             call sourceProperties_Main(time,NumEnsemples,NumSources,realPV)
          end if
          call WriteSourceInfo(time,stossParameter,isut,NumSources,hyperSource)
          call printParticleVector(isut,numEnsemples,numParticles,&
               &                   time,sourceType,realPV)

          !-----------------------------------------------------------------
          !(4) Deallocate "type(quelle) TheSource" and its properties
          !    before the next call at time= t + timeSequence
          !-----------------------------------------------------------------
          call deallocate_source
          if ( getRMF_Flag() ) call deallocate_sourceFields

       end if

       !--------------------------------------------------------------------
       !(4) Stop GiBUU-run after particlevector printed out 10 times after
       !    onset of equilibration.
       !--------------------------------------------------------------------
       if (iprint==MaxTimePrinting) then
          iprint = iprint + 1
          stopGiBUU = .true.
       end if

    end if

  end subroutine DoSourceAnalysis


  !****************************************************************************
  !****s* sourceAnalysis/resetSMMvalues
  ! NAME
  ! subroutine resetSMMvalues
  !
  ! PURPOSE
  ! * Reset global variables to default values for the next subsequent run.
  ! * called in main.f90
  !****************************************************************************
  subroutine resetSMMvalues

    isut         = isut + 1 !next subsequent run

    !reset default values for the next subsequent run:
    Time_out    = 0.0      ! reset to its default value
    itSteps     = 0
    iprint      = 0        !initialize again number of printings
    stopGiBUU   = .false.  !initialize again switch to STOP BUU-run

  end subroutine resetSMMvalues


  !****************************************************************************
  !****f* sourceAnalysis/getSMM_Flag
  ! NAME
  ! logical function getSMM_flag()
  ! PURPOSE
  ! Return the value of SMM_flag. Reads NAMELIST 'SMM_input' from jobCard if
  ! necessary
  !****************************************************************************
  function getSMM_flag() Result (flag)
    use output
    use constants, only: singlePrecision
    implicit none
    integer :: ios
    logical :: flag
    logical, save :: init_getSMM_flag=.true.

    NAMELIST /SMM_input/ SMM_flag, rho_cutoff, spectator_cutoff, A_cutoff, &
         & SelectionMethod,betaChoice,MaxTimePrinting,DetailedHyperonOutput, &
         & hyperSource

    if (init_getSMM_flag) then

       call Write_ReadingInput('SMM_input',0)
       rewind(5)
       read(5,nml=SMM_input,iostat=ios)
       call Write_ReadingInput('SMM_input',0,ios)
       if (SMM_Flag) then
          write(*,*) ' Set SMM_flag to', SMM_flag,'.'
          write(*,*) ' Set rho_cutoff to', rho_cutoff,'.'
          write(*,*) ' Set spectator_cutoff to', spectator_cutoff,'.'
          write(*,*) ' Set A_cutoff to', A_cutoff,'.'
          write(*,*) ' Set SelectionMethod to', SelectionMethod,'.'
          write(*,*) ' Set betaChoice to', betaChoice,'.'
          write(*,*) ' Max. Number of printing', MaxTimePrinting,'.'
          write(*,*) ' DetailedHyperonOutput ',DetailedHyperonOutput,'.'
          write(*,*) ' hyperSource ', hyperSource
          write(*,*) ' Single/double Precision ?', singlePrecision,'.'
       end if
       call Write_ReadingInput('SMM_input',1)

       init_getSMM_flag = .false.

    end if

    flag=SMM_flag

  end function getSMM_flag


  !****************************************************************************
  !****s* sourceAnalysis/printParticleVector
  ! NAME
  ! subroutine printParticleVector
  ! PURPOSE
  ! Prints out the particle vector.
  !****************************************************************************
  subroutine printParticleVector(isut,NumEns,NumPart,time,&
       &                         sourceType,realPV,FinalFlag)
    use particleDefinition
    use output, only: realTochar
    implicit none
    integer,           intent(in) :: isut,NumEns,NumPart
    real,              intent(in) :: time
    logical, optional, intent(in) :: FinalFlag
    type(particle), dimension(:,:),                intent(in) :: realPV
    integer,        dimension(1:NumEns,1:NumPart), intent(in) :: sourceType

    integer :: i,j, ID
    logical :: flag

    if (.not.present(FinalFlag)) then

       open(103,file='auauSMM'//realTochar(time)//'fmc.dat',position='Append')
       if (isut==0) then
          if (DetailedHyperonOutput) then
             write(103,*) NumEns, 2
          else
             write(103,*) NumEns, 1
          end if
          write(103,*) '# time = ',time,' fm/c'
          write(103,*)
       end if

       do i = 1,size(realPV,dim=1) ! ensembles
          do j = 1,size(realPV,dim=2) ! particles
             ID = realPV(i,j)%ID
             if (ID == 0) cycle
             if (ID < 0) exit
             if (sourceType(i,j).ne.999) cycle !print out only emitted particles (mesons included)

             flag = DetailedHyperonOutput .and. (ID==32 .or. ID==33 .or. ID==101)

             if (realPV(i,j)%antiparticle) ID = -ID

             if ( flag ) then !detailed output for hyperons & pions
                write(103,49) &
                     & realPV(i,j)%number, &
                     & ID, &
                     & realPV(i,j)%productionTime, &
                     & realPV(i,j)%lastCollisionTime, &
                     & realPV(i,j)%history, &
                     & realPV(i,j)%charge, &
                     & realPV(i,j)%mass, &
                     & realPV(i,j)%position(1:3),&
                     & realPV(i,j)%momentum(1:3),i,isut+1

             else       ! standard output for all other emitted particles
                write(103,50) ID, &
                     & realPV(i,j)%charge, &
                     & realPV(i,j)%mass, &
                     & realPV(i,j)%position(1:3), &
                     & realPV(i,j)%momentum(1:3),i,isut+1
             end if

          end do
       end do

       close(103)

    else

       open(105,file='auauFile_noresonances.dat',position='Append')

       do i = 1,size(realPV,dim=1) ! ensembles
          do j = 1,size(realPV,dim=2) ! particles
             ID = realPV(i,j)%ID
             if (ID == 0) cycle
             if (ID < 0) exit
             if (realPV(i,j)%antiparticle) ID = -ID

             write(105,51) ID, &
                  & sourceType(i,j), &
                  & realPV(i,j)%charge, &
                  & realPV(i,j)%mass, &
                  & realPV(i,j)%position(1:3), &
                  & realPV(i,j)%momentum(1:3),i,isut+1, time
          end do
       end do

       close(105)

    end if


49  format(i8,1x,i3,1x,2(1x,f7.2),1x,i10,1x,i3,1x,f7.4,6(1x,f9.3),1x,i5,1x,i3)
50  format(i4,1x,i3,1x,f7.4,6(1x,f9.3),1x,i5,1x,i3)
51  format(2(1x,i4),1x,i3,1x,f7.4,6(1x,f9.3),1x,i5,1x,i3,1x,f7.3)

  end subroutine printParticleVector

end module sourceAnalysis
