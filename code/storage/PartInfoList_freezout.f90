!******************************************************************************
!****m* /PIL_freezeout
! NAME
! module PIL_freezeout
! PURPOSE
! Provide some storage method for the freeze out position of some particles
!
! This is closely connected to the module PILIndex
!
! INPUTS
! (none)
!
! NOTES
! * "PIL" stands for "PartInfoList"
!******************************************************************************
module PIL_freezeout

  use PILIndex

  implicit none
  private

  !****************************************************************************
  !****t* PIL_freezeout/tValueEntry
  ! NAME
  ! type tValueEntry
  ! PURPOSE
  ! a container of the information to be stored
  ! SOURCE
  !
  type tValueEntry
     real, dimension(0:3) :: val = 0
     integer :: history = 0
     logical :: escaped = .false.
  end type
  !****************************************************************************


  !****************************************************************************
  !****g* PIL_freezeout/IndexList
  ! PURPOSE
  ! The list, were the particle numbers are connected with the (physical)
  ! storage position
  !
  ! SOURCE
  !
  type(tIndexList), save :: IndexList
  !****************************************************************************


  !****************************************************************************
  !****g* PIL_freezeout/ValueList
  ! PURPOSE
  ! The list, were the information is stored
  !
  ! SOURCE
  !
  type(tValueEntry), allocatable, save :: ValueList(:)
  !****************************************************************************


  public :: PIL_freezeout_PUT, PIL_freezeout_GET
  public :: PIL_freezeout_DeAllocate
  public :: PIL_freezeout_ZERO
!   PUBLIC :: PIL_freezeout_Print

contains

  !****************************************************************************
  !****s* PIL_freezeout/PIL_freezeout_DeAlloc
  ! NAME
  ! subroutine PIL_freezeout_DeAlloc
  ! PURPOSE
  ! Deallocate the memory for this list and the corresponding index list.
  !****************************************************************************
  subroutine PIL_freezeout_DeAllocate()
    call PILIndex_DeAllocate(IndexList)
    if (allocated(ValueList)) deallocate(ValueList)
  end subroutine PIL_freezeout_DeAllocate


  !****************************************************************************
  !****s* PIL_freezeout/PIL_freezeout_ZERO
  ! NAME
  ! subroutine PIL_freezeout_ZERO()
  ! PURPOSE
  ! Reset the list by setting the counter of stored information to 0.
  ! No allocation or deallocation of memory happens.
  !****************************************************************************
  subroutine PIL_freezeout_ZERO()
    IndexList%nEntry = 0
  end subroutine PIL_freezeout_ZERO


  !****************************************************************************
  !****s* PIL_freezeout/PIL_freezeout_PUT
  ! NAME
  ! subroutine PIL_freezeout_PUT(number,r)
  ! PURPOSE
  ! Store the information "r" connected with particle "number" in the list.
  ! INPUTS
  ! * integer :: number -- the (unique) particle number
  ! * real    :: r      -- the 4-position to store
  ! * integer :: hist   -- the history of the particle
  ! * logical :: escaped-- flag to indicate, whether particle has 'escaped'
  ! OUTPUT
  ! none
  !****************************************************************************
  subroutine PIL_freezeout_PUT(number,r,hist,escaped)
    integer, intent(IN) :: number
    real, dimension(0:3), intent(IN) :: r
    integer, intent(IN) :: hist
    logical, intent(IN) :: escaped

    integer :: iEntry

    iEntry = PILIndex_PUT(IndexList, number,"freezeout")

!    write(*,*) '###PUT: ',number,r,iEntry

    if (iEntry>0) then      ! everything is ok
       ValueList(iEntry)%val = r
       ValueList(iEntry)%history = hist
       ValueList(iEntry)%escaped = escaped
    else
       call PIL_freezeout_Allocate() ! do (re)allocate
       ValueList(-iEntry)%val = r
       ValueList(-iEntry)%history = hist
       ValueList(-iEntry)%escaped = escaped
    end if
  end subroutine PIL_freezeout_PUT


  !****************************************************************************
  !****f* PIL_freezeout/PIL_freezeout_GET
  ! NAME
  ! logical function PIL_freezeout_GET(number,r,hist,escaped)
  ! PURPOSE
  ! Get the stored information of particle "number"
  ! INPUTS
  ! * integer :: number -- the (unique) particle number
  ! OUTPUT
  ! * real    :: r      -- the 4-position to store
  ! * integer :: hist   -- the history of the particle
  ! * logical :: escaped-- flag to indicate, whether particle has 'escaped'
  ! * the (logical) return value signals, whether inforation about
  !   this particle was found in the list or not.
  !****************************************************************************
  logical function PIL_freezeout_GET(number,r,hist,escaped)
    integer, intent(IN)  :: number
    real, dimension(0:3), intent(OUT) :: r
    integer, intent(OUT) :: hist
    logical, intent(OUT) :: escaped

    integer :: iEntry

    iEntry = PILIndex_FIND(IndexList,number)

    ! ATTENTION: iEntry is the line of information in the IndexList.
    ! The information connected to the particle is stored in the line
    ! IndexList%Entry(iEntry) in the array ValueList !!!!!

    if (iEntry > 0) then
       r = ValueList(IndexList%Entry(iEntry))%val
       hist = ValueList(IndexList%Entry(iEntry))%history
       escaped = ValueList(IndexList%Entry(iEntry))%escaped
       PIL_freezeout_GET = .TRUE.
    else
       r = 0.0
       hist = 0
       escaped = .false.
       PIL_freezeout_GET = .FALSE.
    end if

  end function PIL_freezeout_GET


  !****************************************************************************
  !****is* PIL_freezeout/PIL_freezeout_Allocate
  ! NAME
  ! subroutine PIL_freezeout_Allocate
  ! PURPOSE
  ! Do the allocation and reallocation of the value vector.
  ! The new size is taken from the size of the IndexList vectors.
  ! NOTES
  ! For security one should insert here checks, whether the memory allocations
  ! failed and stop execution in these cases.
  !****************************************************************************
  subroutine PIL_freezeout_Allocate
    integer :: n0, n1,i
    type(tValueEntry), allocatable :: L0(:)

    n1 = size(IndexList%PartNumber) ! new size

    if (.not.allocated(ValueList)) then
       allocate(ValueList(n1))
       return
    end if

    n0 = size(ValueList)            ! old size

    allocate(L0(n0))
    do i=1,n0
       L0(i)%val = ValueList(i)%val
       L0(i)%history = ValueList(i)%history
       L0(i)%escaped = ValueList(i)%escaped
    end do
    deallocate(ValueList)
    allocate(ValueList(n1))
    do i=1,n0
       ValueList(i)%val = L0(i)%val
       ValueList(i)%history = L0(i)%history
       ValueList(i)%escaped = L0(i)%escaped
    end do
    deallocate(L0)

  end subroutine PIL_freezeout_Allocate


  !****************************************************************************
  !****s* PIL_freezeout/PIL_freezeout_Print
  ! NAME
  ! subroutine PIL_freezeout_Print(file)
  ! PURPOSE
  ! Print the list to file
  !****************************************************************************
!   subroutine PIL_freezeout_Print(file)
!     integer, intent(IN) :: file
!
!     integer :: i
!
!     write(file,*) '****** PIL_freezeout:',IndexList%nEntry
!     do i=1,IndexList%nEntry
!        write(file,'(i8.0,i8.0,3g12.5)') i,IndexList%PartNumber(i),ValueList(IndexList%Entry(i))%val
!     end do
!
!   end subroutine PIL_freezeout_Print



end module PIL_freezeout
