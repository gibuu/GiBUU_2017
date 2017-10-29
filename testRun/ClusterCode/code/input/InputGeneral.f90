!***************************************************************************
!****m* /inputGeneral
! NAME
! module inputGeneral
! PURPOSE
! This is the main module for reading some general parameters concerning 
! the whole cluster code
!***************************************************************************
Module InputGeneral

  !*************************************************************************
  !****g* InputGeneral/TheModel
  ! SOURCE
  !
  Integer,save :: TheModel=0
  !
  ! PURPOSE
  ! Defines the model applied for nuclear fragmentation
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputGeneral/SubEvents
  ! SOURCE
  !
  Integer,save :: SubEvents = 0
  !
  ! PURPOSE
  ! Number of sub-sequent runs at same beam energy
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputGeneral/NumEnsemples
  ! SOURCE
  !
  Integer,save :: NumEnsemples = 0
  !
  ! PURPOSE
  ! Number of ensemples (=test-particles per nucleon)
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputGeneral/Realparticles
  ! SOURCE
  !
  Integer,save :: RealParticles = 0
  !
  ! PURPOSE
  ! Suggested maximal number of (real) particles per event.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputGeneral/BUU_DataFile
  ! SOURCE
  !
  Character(60),save :: BUU_DataFile='BUU_File'
  !
  ! PURPOSE
  ! GiBUU Data File to be analyzed
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputGeneral/PathToBUUInput
  ! SOURCE
  !
  character(100), save :: pathToBUUInput='~/'
  !
  ! PURPOSE
  ! Link to the directory of GiBUU files to be analyzed
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/Get_Hyp
  ! SOURCE
  !
  Logical,save :: Get_Hyp = .false.
  !
  ! PURPOSE
  ! Switch for Coalescence between hyperons and clusters. 
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/Get_GiBUUvec
  ! SOURCE
  !
  Logical,save :: Get_GiBUUvec = .false.
  !
  ! PURPOSE
  ! Switch to read the GiBUU particle vector. 
  ! NOTES
  ! Needed when: 
  ! (1) emitted particles are included
  ! (2) production of hypernuclei is included
  ! in the analysis (Get_Hyp=.true.).
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/ALADIN_Flag
  ! SOURCE
  !
  Logical,save :: ALADIN_Flag = .false.
  !
  ! PURPOSE
  ! if true, then ALADIN/INDRA-data for Spectator-Fragmentation only are 
  ! analyzed. Fireball matter is not included in the analysis. 
  ! if false, then FOPI-data for fireballs could be analyzed. 
  ! This part is still under construction. 
  ! NOTES
  ! This flag is needed only for Heavy-Ion runs. If false, then fireball 
  ! matter is analyzed, where mostly FOPI-date exist. Note, that this option 
  ! works only for "Get_GiBUUvec=.true.", because the particle vector is 
  ! needed to calculate radial flow profiles. 
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/CHARMS_Flag
  ! SOURCE
  !
  Logical,save :: CHARMS_Flag = .false.
  !
  ! PURPOSE
  ! if true, then CHARMS-data for Spectator-Fragmentation are 
  ! analyzed, otherwise FOPI-data.
  !
  !*************************************************************************



!==========================================================================*
CONTAINS
!==========================================================================*

  !*************************************************************************
  !****s* inputGeneral/Get_InputGeneral
  ! NAME
  ! subroutine Get_InputGeneral
  ! PURPOSE
  ! Reads input in jobcard out of namelist "inputGeneral"
  !***********************************************************************

  subroutine Get_InputGeneral
    implicit none
    integer :: ios
    !-----------------------------------------------------------------------
    NAMELIST /inputGeneral/TheModel,SubEvents,NumEnsemples,RealParticles, & 
         & Get_Hyp,ALADIN_Flag,CHARMS_Flag, Get_GiBUUvec, & 
         & pathToBUUInput, BUU_DataFile

    write(*,*) '----- Reading Namelist "inputGeneral": Start -----'
    rewind(5)
    read(5,nml=inputGeneral,IOSTAT=ios)
    !-----------------------------------------------------------------------
    if(ios.ne.0) then
       write(*,*) 'Error in  namelist "inputGeneral" : This namelist is crucial. STOP!'
       stop
    end if
    !-----------------------------------------------------------------------

    if ( (Get_Hyp) .and. (.not.Get_GiBUUvec) ) then
       write(*,*) 
       write(*,*) 'Get_InputGeneral: you need the GiBUU particle vector '
       write(*,*) 'for production of hypernuclei!!!'
       write(*,*) 'Get_GiBUUvec = ',Get_GiBUUvec
       write(*,*) 'Change Get_GiBUUvec to true...'
       Get_GiBUUvec = .true.
       write(*,*) 'actual value: Get_GiBUUvec = ',Get_GiBUUvec
    endif

    if (SubEvents .eq. 0) then
       write(*,*) 'SubEvents .eq. 0 ! STOP'
       stop
    end if

    if ( (.not. ALADIN_Flag) .and. (.not.Get_GiBUUvec) ) then
       write(*,*) 'Get_InputGeneral: you need the GiBUU particle vector '
       write(*,*) 'for the analysis of fireball matter in HIC-runs!!!'
       write(*,*) 'STOP'
       STOP
    end if

    write(*,*) '----- Reading Namelist "inputGeneral": END   -----'
    write(*,*)

  !*************************************************************************
  end subroutine Get_InputGeneral !*****************************************
  !*************************************************************************


end module InputGeneral
