!******************************************************************************
!****m* /history
! NAME
! module history
! PURPOSE
! Here the %history field of type(particle) is dealt with. This field is used
! to store a history of the collisions of a particle.
!
! In particular, we save the particle IDs of those particles which produced
! the particle ('parents'). E.g. if a pion was produced in AB -> X pion,
! then we store the IDs of the two incoming particles "A" and "B".
!
! Furthermore, we define a so-called "generation". Initially all particles
! have generation=0. Particles which are created by particles of generation x
! and y get a value for the generation of generation=max( x, y )+1.
! Therefore this "generation" is a measure for the number of interactions
! which were necessary to produce a particle. For perturbative particles in
! ground state calculations,"generation" gives directly the number of
! interactions which were necessary to produce a particle.
!******************************************************************************
module history

  implicit none
  private

  integer, parameter :: million=1000000
  integer, parameter :: nnn=1          ! Identifier for a (N N N) interaction
  integer, parameter :: nnDelta=2      ! Identifier for a (N N Delta) interaction
  integer, parameter :: nnPion=3       ! Identifier for a (N N Pion) interaction
  integer, parameter :: some3Body=4    ! Identifier for any other 3Body interaction being not NNN, NN Delta or NN Pion.


  !****************************************************************************
  !****g* history/IncGeneration_Decay
  ! SOURCE
  !
  logical, save :: IncGeneration_Decay = .true.
  ! PURPOSE
  ! This flag determines whether we will increase the stored
  ! 'generation' of the daughter particles in a resonance decay.
  !****************************************************************************

  !****************************************************************************
  !****g* history/IncGeneration_Elastic
  ! SOURCE
  !
  logical, save :: IncGeneration_Elastic = .true.
  ! PURPOSE
  ! This flag determines whether we will increase the stored 'generation'
  ! of particles in an elastic collision. Setting it to .false. will also
  ! prevent elastic collisions from showing up as parents in the history.
  !****************************************************************************


  !****************************************************************************
  !****s* history/setHistory
  ! NAME
  ! interface setHistory
  !
  ! PURPOSE
  ! Evaluates the value of %history for a particle which was produced in
  ! a interaction, which was either a resonance decay, 2 body interaction or 3
  ! body interaction.
  !
  ! NOTES
  ! We use the following rule to determine the value of %history:
  ! * For 1-Body processes : history=1.000.000*generation+id1
  ! * For 2-Body processes : history=1.000.000*generation+1.000*id2+id1
  ! * For 3-Body processes : history=-(1.000.000*generation+x) where x is
  !   defined by the identifiers given as parameters in the module header.
  !
  ! This scheme is employed to store this information in a compact way into one integer,
  ! and to save memory by doing this.
  !****************************************************************************
  interface setHistory
    module procedure setHistory1, setHistory2, setHistory3
  end interface

  public :: setHistory, history_getParents, history_getGeneration, history_print

  logical, save :: initFlag = .true.

contains

  !****************************************************************************
  !****s* history/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "History".
  !****************************************************************************
  subroutine readInput
    !use output           ! using output not possible (circular)
    integer :: ios

    !**************************************************************************
    !****n* history/History
    ! NAME
    ! Namelist /History/
    ! PURPOSE
    ! Namelist for history includes:
    ! * IncGeneration_Decay
    ! * IncGeneration_Elastic
    !**************************************************************************
    NAMELIST /history/ IncGeneration_Decay, IncGeneration_Elastic

1000 FORMAT (/,'------ Init "',A,'": reading...')
1001 FORMAT ('------ Init "',A,'": reading finished.',/)
1002 FORMAT (/,'--- !!!!! ERROR while reading namelist "',A,'" !!!!! STOPPING !!',/)
1003 FORMAT ('--- *** namelist "',A,'" missing in jobcard. Values-->Default!')

    !call Write_ReadingInput('history',0)
    write(*,1000) 'history'

    rewind(5)
    read(5,nml=history,iostat=ios)

    !call Write_ReadingInput('history',0,ios)
    if (ios > 0) then
      write(*,1002) 'history'
      stop
    else if (ios < 0) then
      write(*,1003) 'history'
    end if

    write(*,*) 'Increase Generation in Decay:        ', IncGeneration_Decay
    write(*,*) 'Increase Generation in Elastic Coll: ', IncGeneration_Elastic

    !call Write_ReadingInput('history',1)
    write(*,1001) 'history'

    initFlag = .false.

  end subroutine readInput



  subroutine setHistory1 (p1, final)
    use particleDefinition
    type(particle), intent(in) :: p1
    type(particle), dimension(:) :: final
    integer :: id1, generation

    if (initFlag) call readInput

    id1 = p1%ID
    if (p1%antiparticle) id1 = id1 + 200

    generation = abs(p1%history)/million
    if (IncGeneration_Decay) generation = generation + 1

    final%history = id1 + generation*million

  end subroutine setHistory1



  subroutine setHistory2 (p, final)
    use particleDefinition
    use twoBodyTools, only: IsElastic
    type(particle), dimension(2), intent(in) :: p
    type(particle), dimension(:) :: final
    integer :: id1, id2, generation

    if (initFlag) call readInput

    if (final(3)%ID<=0 .and. IsElastic(p,final(1:2)) .and. .not.IncGeneration_Elastic) then
      if (isSamePart(p(1),final(1))) then
        final(1)%history = p(1)%history
        final(2)%history = p(2)%history
      else
        final(1)%history = p(2)%history
        final(2)%history = p(1)%history
      end if
      return
    end if

    generation = max (abs(p(1)%history), abs(p(2)%history)) / million + 1

    id1 = p(1)%ID
    if (p(1)%antiparticle) id1 = id1 + 200

    id2 = p(2)%ID
    if (p(2)%antiparticle) id2 = id2 + 200

    final%history = id1 + 1000*id2 + generation*million

  end subroutine setHistory2



  subroutine setHistory3 (p1, p2, p3, final)
    use idTable, only: nucleon, delta, pion
    use particleDefinition
    type(particle), intent(in) :: p1,p2,p3
    type(particle), dimension(:) :: final
    integer :: generation, id

    generation = max(abs(p1%history),abs(p2%history),abs(p3%history))/million + 1

    if ((p1%ID==nucleon).and.(p2%ID==nucleon).and.(p3%ID==nucleon)) then
      id = -nnn      ! ---- NNN channel
    else if ( ((p1%ID==nucleon).and.(p2%ID==nucleon).and.(p3%ID==delta)) .or. &
              ((p1%ID==nucleon).and.(p3%ID==nucleon).and.(p2%ID==delta)) .or. &
              ((p2%ID==nucleon).and.(p3%ID==nucleon).and.(p1%ID==delta)) ) then
      id = -nnDelta  ! ---- NNDelta channel
    else if ( ((p1%ID==nucleon).and.(p2%ID==nucleon).and.(p3%ID==pion)) .or. &
              ((p1%ID==nucleon).and.(p3%ID==nucleon).and.(p2%ID==pion)) .or. &
              ((p2%ID==nucleon).and.(p3%ID==nucleon).and.(p1%ID==pion)) ) then
      id = -nnPion   ! ---- NNpi channel
    else
      id = -some3Body
    end if

    final%history = id - generation*million

  end subroutine setHistory3



  !****************************************************************************
  !****s* history/history_getParents
  ! NAME
  ! function history_getParents (history) result (parents)
  !
  ! PURPOSE
  ! Analyzes a given value of %history of a particle and returns an array of integers which include
  ! the IDs of the particles which produced the latter particle.
  !
  ! INPUTS
  ! * integer, intent(in) :: history -- value of %history of the regarded particle
  !
  ! OUTPUT
  ! * integer, dimension(1:3) :: parents -- parents of the regarded particle
  !****************************************************************************
  function history_getParents (history) result (parents)
    use idTable, only: nucleon, Delta, pion
    integer, intent(in) :: history
    integer, dimension(1:3) :: parents

    integer :: generation, dummy

    parents = 0
    generation = history_getGeneration(history)

    if (history > 0) then
       ! decays and two-body collisions
       parents(2)=(history-million*generation)/1000
       parents(1)=history-million*generation-1000*parents(2)
       if (parents(2)>0 .and. parents(1)>parents(2)) then
          ! sort: smaller ID first
          dummy=parents(2)
          parents(2)=parents(1)
          parents(1)=dummy
       end if
    else
       ! three-body collisions
       select case (-history-million*generation)
       case (nnn)
          parents(1:3) = nucleon
       case (nnDelta)
          parents(1:3) = (/nucleon,nucleon,Delta/)
       case (nnPion)
          parents(1:3) = (/nucleon,nucleon,pion/)
       end select
    end if
  end function history_getParents


  !****************************************************************************
  !****f* history/history_getGeneration
  ! NAME
  ! integer function history_getGeneration(history)
  !
  ! PURPOSE
  ! Analyzes a given value of %history of a particle and returns the value
  ! of generation of the latter particle.
  !
  ! INPUTS
  ! * integer, intent(in) :: history -- value of %history of the regarded particle
  !
  !****************************************************************************
  integer function history_getGeneration(history)
    integer, intent(in) :: history

    history_getGeneration = abs(history)/million

  end function history_getGeneration


!!$  !*************************************************************************
!!$  !****f* history/history_1Body
!!$  ! NAME
!!$  ! logical function history_1Body(history)
!!$  !
!!$  ! PURPOSE
!!$  ! Analyzes a given value of %history of a particle and returns .true. if the
!!$  ! particle was produced in a resonance decay.
!!$  !
!!$  ! INPUTS
!!$  ! * integer, intent(in) :: history -- value of %history of the regarded particle
!!$  !
!!$  !*************************************************************************
!!$  logical function history_1Body(history)
!!$    integer, intent(in) :: history
!!$    integer :: generation
!!$    generation=history_getGeneration(history)
!!$    if((history-million*generation.lt.1000).and.(history-million*generation.gt.0)) then
!!$       history_1Body=.true.
!!$    else
!!$       history_1Body=.false.
!!$    end if
!!$  end function history_1Body
!!$
!!$
!!$  !*************************************************************************
!!$  !****f* history/history_2Body
!!$  ! NAME
!!$  ! logical function history_2Body(history)
!!$  !
!!$  ! PURPOSE
!!$  ! Analyzes a given value of %history of a particle and returns .true. if the
!!$  ! particle was produced in 2 Body interaction.
!!$  !
!!$  ! INPUTS
!!$  ! * integer, intent(in) :: history -- value of %history of the regarded particle
!!$  !
!!$  !*************************************************************************
!!$  logical function history_2Body(history)
!!$    integer, intent(in) :: history
!!$    integer :: generation
!!$    generation=history_getGeneration(history)
!!$    if((history-million*generation.lt.million).and.(history-million*generation.gt.1000)) then
!!$       history_2Body=.true.
!!$    else
!!$       history_2Body=.false.
!!$    end if
!!$  end function history_2Body
!!$
!!$
!!$  !*************************************************************************
!!$  !****f* history/history_3Body
!!$  ! NAME
!!$  ! logical function history_3Body(history)
!!$  !
!!$  ! PURPOSE
!!$  ! Analyzes a given value of %history of a particle and returns .true. if the
!!$  ! particle was produced in 3 Body interaction.
!!$  !
!!$  ! INPUTS
!!$  ! * integer, intent(in) :: history -- value of %history of the regarded particle
!!$  !
!!$  !*************************************************************************
!!$  logical function history_3Body(history)
!!$    integer, intent(in) :: history
!!$    !integer :: generation
!!$    !generation=history_getGeneration(history)
!!$    if(history.lt.0) then
!!$       history_3Body=.true.
!!$    else
!!$       history_3Body=.false.
!!$    end if
!!$  end function history_3Body
!!$

  !****************************************************************************
  !****f* history/history_print
  ! NAME
  ! subroutine history_print(ensemble,p, iFile)
  !
  ! PURPOSE
  ! Print the history of a given particle to channel "ifile"
  !
  ! INPUTS
  ! * integer, intent(in) :: iFile, ensemble
  ! * type(particle),intent(in)  :: p
  !
  !****************************************************************************
  subroutine history_print(ensemble,p, iFile,initFlag)
    use particleDefinition
    integer, intent(in) :: iFile, ensemble
    type(particle),intent(in)  :: p
    integer :: generation, parents(3)
    logical, optional :: initFlag

    generation = history_getGeneration (p%history)
    parents = history_getParents (p%history)
    if (present(initFlag)) then
       if (initFlag)   write(iFile,'(A)') 'ensemble,ID, charge, history, generation, parents(1:3), firstEvent'
    end if

    write(iFile,'(12I10)') ensemble,p%ID, p%charge, p%history, generation, parents,p%firstEvent

  end subroutine history_print


end module history
