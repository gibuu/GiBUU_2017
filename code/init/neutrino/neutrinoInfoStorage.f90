!******************************************************************************
!****m* /neutrinoInfoStorage
! NAME
! module neutrinoInfoStorage
!
! PURPOSE
! This module stores information about the initial neutrino event.
! (structure taken from EventInfo_HiLep)
!
!******************************************************************************
module neutrinoInfoStorage

  implicit none
  private

  !****************************************************************************
  !****t* neutrinoInfoStorage/tneutrinoInfoStorage
  ! SOURCE
  !
  type tneutrinoInfoStorage
     logical :: flagOK=.False.
     real    :: value=0.
     real    :: value_rec=0.
     real    :: value_rec_2=0.
  end type tneutrinoInfoStorage
  !
  ! PURPOSE
  ! This holds all information we want to store connected to
  ! neutrino induced reactions.
  !****************************************************************************

  type(tneutrinoInfoStorage),save,dimension(:),allocatable :: EXP_nuInfo

  public :: neutrinoInfoStorage_Init,neutrinoInfoStorage_Clear !,neutrinoInfoStorage_Store,neutrinoInfoStorage_Get

contains

  !****************************************************************************
  !****s* neutrinoInfoStorage/neutrinoInfoStorage_Init
  ! NAME
  ! subroutine neutrinoInfoStorage_Init(NumInitialEvents)
  !
  ! PURPOSE
  ! allocate memory and reset the corresponding arrays.
  !
  ! INPUTS
  ! * integer :: NumInitialEvents -- number of possible initial events
  !
  ! OUTPUT
  ! ---
  !****************************************************************************
  subroutine neutrinoInfoStorage_Init(NumInitialEvents)
    integer,intent(in) :: NumInitialEvents


    if (allocated(EXP_nuInfo)) then
       deallocate(EXP_nuInfo)
    end if
    allocate(EXP_nuInfo(1:NumInitialEvents))
    EXP_nuInfo%flagOK=.false.
    EXP_nuInfo%value=0.
    EXP_nuInfo%value_rec=0.
    EXP_nuInfo%value_rec_2=0.

  end subroutine NeutrinoInfoStorage_Init

  !****************************************************************************
  !****s* neutrinoInfoStorage/neutrinoInfoStorage_clear
  ! NAME
  ! subroutine neutrinoInfoStorage_clear
  !
  ! PURPOSE
  ! if necessary clear allocated memory.
  !
  ! OUTPUT
  ! ---
  !****************************************************************************
  subroutine neutrinoInfoStorage_clear
    !logical,save :: initFlag=.true.

    if (allocated(EXP_nuInfo)) then
       deallocate(EXP_nuInfo)
    end if

  end subroutine NeutrinoInfoStorage_clear


  !****************************************************************************
  !****s* neutrinoInfoStorage/neutrinoInfoStorage_Store
  ! NAME
  ! subroutine neutrinoInfoStorage_Store(i,value,value_rec)
  !
  ! PURPOSE
  ! Store the event info connected with number "i":
  !
  !
  ! INPUTS
  ! * integer :: i -- actual number of event
  ! * real    :: value -- neutrino kinematic variables
  ! * real    :: value_rec, value_rec_2  -- reconstructed variables
  !
  ! OUTPUT
  ! ---
  !****************************************************************************
!   subroutine neutrinoInfoStorage_Store(i,value,value_rec,value_rec_2)
!     integer,intent(in)          :: i
!     real,   intent(in)          :: value,value_rec
!     real, intent(in), optional :: value_rec_2
!
!     EXP_nuInfo(i)%flagOK=.true.
!     EXP_nuInfo(i)%value=value
!     EXP_nuInfo(i)%value_rec=value_rec
!     if(present(value_rec_2)) EXP_nuInfo(i)%value_rec_2=value_rec_2
!
!   end subroutine NeutrinoInfoStorage_Store


  !****************************************************************************
  !****f* neutrinoInfoStorage/neutrinoInfoStorage_Get
  ! NAME
  ! logical function neutrinoInfoStorage_Get(i,value,value_rec)
  !
  ! PURPOSE
  ! Get the event info stored connected with number "i".
  !
  ! INPUTS
  ! * integer :: i -- actual number of event
  !
  ! OUTPUT
  ! * real    :: value -- neutrino kinematic variables
  ! * real    :: value_rec -- reconstructed variables
  !
  !****************************************************************************
!   logical function neutrinoInfoStorage_Get(i,value,value_rec,value_rec_2)
!     integer,intent(in)           :: i
!     real,   intent(out)          :: value,value_rec
!     real, intent(out), optional  :: value_rec_2
!
!     neutrinoInfoStorage_Get = .FALSE.
!
!     if (.not.ALLOCATED(EXP_nuInfo)) return
!
!     value=EXP_nuInfo(i)%value
!     value_rec=EXP_nuInfo(i)%value_rec
!
!     if(present(value_rec_2)) value_rec_2=EXP_nuInfo(i)%value_rec_2
!
!     neutrinoInfoStorage_Get = EXP_nuInfo(i)%flagOK
!
!
!   end function NeutrinoInfoStorage_Get

end module NeutrinoInfoStorage
