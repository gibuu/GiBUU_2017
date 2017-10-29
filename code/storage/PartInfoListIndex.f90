!******************************************************************************
!****m* /PILIndex
! NAME
! module PILIndex
! PURPOSE
! Provide a multi usable list connecting the (unique) particle number with
! the index number of some additional information in some other array.
!
! While the list in the 'other array' is filled entry by entry, the list
! provided here is sorted by the index number in increasing order.
! This is done in order to speed up finding information about a particle
! from the 'other array' and whether information about this specific
! particle is actually stored in the array.
!
! This implementation speeds up getting information at the expense of
! inserting new information.
!
! INPUTS
! (none)
! NOTES
! * ...maybe one could change the access to PIL_nEntry0
! * "PIL" stands for "PartInfoList"
!******************************************************************************
module PILIndex
  use CallStack, only: Traceback

  implicit none
  private

  !****************************************************************************
  !****t* PILIndex/tIndexList
  ! NAME
  ! type tIndexList
  ! PURPOSE
  ! store all information of a IndexList
  ! SOURCE
  !
  type, public :: tIndexList
     integer :: nEntry = -999              ! number of used entries
     integer :: nHole  = 0                 ! number of deleted entries
     integer, allocatable :: PartNumber(:) ! array of particle numbers
     integer, allocatable :: Entry(:)      ! array of indices
  end type tIndexList
  !****************************************************************************

  !****************************************************************************
  !****g* PILIndex/PIL_nEntry0
  ! PURPOSE
  ! default size of arrays (at the beginning)
  !
  ! SOURCE
  !
  integer,parameter :: PIL_nEntry0 = 1000
  !****************************************************************************

  public :: PILIndex_FIND, PILIndex_PUT, PILIndex_DeAllocate

contains

  !****************************************************************************
  !****f* PILIndex/PILIndex_PUT
  ! NAME
  ! integer function PILIndex_PUT(List, number)
  ! PURPOSE
  ! find the place, where some information concerning a particle with "number"
  ! should be put into some info vector.
  !
  ! INPUTS
  ! * type(tIndexList) :: List   -- The List to put the entry in
  ! * integer          :: number -- The particle number to insert
  ! * character(*),intent(in),optional :: CallName -- Name of the List
  !
  ! OUTPUT
  ! * List is modified
  ! * return value gives the position where the info to be stored.
  !   If this value is negative, some allocation/reallocation of the
  !   infovector has to be done, then the info can be stored at the position
  !   given by the absolute value of the return value.
  !****************************************************************************
  integer function PILIndex_PUT(List, number, CallName)

    type(tIndexList), intent(INOUT) :: List
    integer,          intent(IN)    :: number
    character(*),intent(in),optional :: CallName
    integer :: i1
    logical :: flagAlloc
    integer :: newEntry, newPos

    flagAlloc = .false.

    !=== first call: ===
!    if (List%nEntry == -999) then
    if (.not.allocated(List%PartNumber)) then
       call PILIndex_Allocate(List)
       List%nEntry = 1
       newEntry = 1
       newPos = 1
       List%PartNumber(newEntry) = number
       List%Entry(newEntry)      = newPos
       PILIndex_PUT              = -newPos ! sign because of allocate
       return
    end if

    !=== empty list: ===
    if (List%nEntry == 0) then
       List%nEntry = 1
       newEntry = 1
       newPos = 1
       List%PartNumber(newEntry) = number
       List%Entry(newEntry)      = newPos
       PILIndex_PUT              = newPos
       return
    end if

    !=== number larger than largest stored number --> add at the end: ===
    ! (This is a duplicate of code from below; just as a shortcut
    ! Maybe we should delete it for easier maintanance !!!)
    !
    if (number > List%PartNumber(List%nEntry)) then
       if (List%nEntry >= size(List%PartNumber)) then
          if (present(CallName)) then
             write(*,*) 'PILIndex: reallocate ',size(List%PartNumber),'->... [',CallName,']'
          else
             write(*,*) 'PILIndex: reallocate ',size(List%PartNumber),'->...'
          end if
          call PILIndex_Allocate(List)
          flagAlloc = .true.
       end if
       List%nEntry = List%nEntry + 1
       newEntry = List%nEntry
       if (List%nHole==0) then
          newPos = List%nEntry
       else
          newPos = List%Entry(size(List%Entry)-List%nHole)
          List%nHole = List%nHole -1
       end if
       List%PartNumber(newEntry) = number
       List%Entry(newEntry)      = newPos
       if (flagAlloc) then
          PILIndex_PUT  = -newPos
       else
          PILIndex_PUT  = newPos
       end if
       return
    end if

!    call PILIndex_Print(List,6)

    !=== search for entry position: ===
    i1 = PILIndex_FIND(List, number)
    if (i1 > 0) then !=== entry found !!!
       PILIndex_PUT = List%Entry(i1)
       return
    end if

    !=== now we have to insert the info at position (-i1+1) !!!

    !=== allocate new space, if needed: ===
    if (List%nEntry >= size(List%PartNumber)) then
       if (present(CallName)) then
          write(*,*) 'PILIndex: reallocate ',size(List%PartNumber),'->... [',CallName,']'
       else
          write(*,*) 'PILIndex: reallocate ',size(List%PartNumber),'->...'
       end if
       call PILIndex_Allocate(List)
       flagAlloc = .true.
    end if

    !=== shift entry information by copying it: ===
    List%PartNumber(-i1+2:List%nEntry+1) = List%PartNumber(-i1+1:List%nEntry)
    List%Entry(-i1+2:List%nEntry+1) = List%Entry(-i1+1:List%nEntry)

    !=== insert the info at position (-i1+1) ===
    List%nEntry = List%nEntry + 1
    newEntry = -i1+1
    if (List%nHole==0) then
       newPos = List%nEntry
    else
       newPos = List%Entry(size(List%Entry)-List%nHole)
       List%nHole = List%nHole -1
    end if
    List%PartNumber(newEntry) = number
    List%Entry(newEntry)      = newPos
    if (flagAlloc) then
       PILIndex_PUT  = -newPos
    else
       PILIndex_PUT  = newPos
    end if

!   call PILIndex_Print(List,6)

  end function PILIndex_PUT


  !****************************************************************************
  !****f* PILIndex/PILIndex_FIND
  ! NAME
  ! integer function PILIndex_FIND(List, number)
  ! PURPOSE
  ! find the place, where the "number" is in the list
  ! or where it should be stored
  !
  ! INPUTS
  ! * type(tIndexList) :: List   -- The List
  ! * integer          :: number -- The particle number to find
  !
  ! OUTPUT
  ! if the return value R is > 0:
  ! here the entry "number" was found in the list
  !
  ! if the return value R <= 0:
  ! This means, that the entry at position (-R) is smaller than number,
  ! while the entry at position (-R)+1 is already larger.
  !
  ! NOTES
  ! the return value R==0 is possible. This means of course, that already
  ! the first entry is larger than number.
  ! But you should never (!) try to access some entries at position R==0,
  ! because this undershoots the lower bound of the arrays.
  !****************************************************************************
  integer function PILIndex_FIND (List, number)

    type(tIndexList), intent(IN) :: List
    integer,          intent(IN) :: number

    integer :: jL, jU, jM

    if (List%nEntry<=0) then
       ! No entry available
       PILIndex_FIND=0
       return
    end if

    jL = 0
    jU = List%nEntry + 1

    if (number > List%PartNumber(jU-1)) then
       PILIndex_FIND = -jU
       return
    else if (number < List%PartNumber(jL+1)) then
       PILIndex_FIND = 0
       return
    end if

    do
       if (jU-jL <= 1) exit
       jM = (jU+jL)/2
       if (number .ge. List%PartNumber(jM)) then
          jL = jM
       else
          jU = jM
       end if
    end do

    if (number == List%PartNumber(1) ) then
       PILIndex_FIND = 1
    else if (number == List%PartNumber(List%nEntry) ) then
       PILIndex_FIND = List%nEntry
    else
       PILIndex_FIND = -jL
       if (jL<1) return
       if (number==List%PartNumber(jL)) PILIndex_FIND = jL
    end if

  end function PILIndex_FIND


  !****************************************************************************
  !****s* PILIndex/PILIndex_Print
  ! NAME
  ! subroutine PILIndex_Print(List, file)
  ! PURPOSE
  ! just print the list to a file stream
  !
  ! INPUTS
  ! * type(tIndexList) :: List -- The List
  ! * integer          :: file -- The file number
  !****************************************************************************
!   subroutine PILIndex_Print(List, file)
!
!     type(tIndexList), intent(IN) :: List
!     integer,          intent(IN) :: file
!
!     integer :: i
!
!     do i=1,List%nEntry
!        write(file,'(i5.0,2i7.0)') i,List%PartNumber(i),List%Entry(i)
!     end do
!   end subroutine PILIndex_Print


  !****************************************************************************
  !****s* PILIndex/PILIndex_Allocate
  ! NAME
  ! subroutine PILIndex_Allocate(List)
  ! PURPOSE
  ! Do all allocation and reallocatiomn stuff connected with "List"
  !
  ! INPUTS
  ! * type(tIndexList) :: List -- The List
  ! OUTPUT
  ! * if the arrays of "list" were not allocated before, they are allocated with
  !   some default size
  ! * if the arrays were allocated, they are reallocated with a size 1.3 times larger
  !   than the original size
  ! NOTES
  ! I think to remeber having read that the factor 1.3 is a good compromise
  ! between the effort of copying around all data and how often this happens.
  ! Of course this implies that one starts with some reasonable value.
  ! Maybe one could improve. (Kai)
  !
  ! For security one should insert here checks, whether the memory allocations
  ! failed and stop execution in these cases.
  !
  ! It would be a good idea to check whether the compiler knows the intrinsic
  ! MOVE_ALLOC subroutine and use this instead of copying the contents twice.
  !****************************************************************************
  subroutine PILIndex_Allocate(List)

    type(tIndexList), intent(INOUT) :: List

    integer, allocatable :: L0(:)
    integer :: n0, n1

    !=== first call: ===
    if (.not.allocated(List%PartNumber)) then
       allocate(List%PartNumber(PIL_nEntry0))
       allocate(List%Entry(PIL_nEntry0))
       List%nEntry = 0
       List%nHole  = 0
       return
    end if

    if (List%nHole>0) call Traceback("Allocate not necessary with nHole>0")

    n0 = size(List%PartNumber)
    n1 = int ( n0 * 1.3 )

    allocate(L0(n0))

    L0 = List%PartNumber

    deallocate(List%PartNumber)
    allocate(List%PartNumber(n1))

    List%PartNumber(1:n0) = L0(:)

    L0 = List%Entry

    deallocate(List%Entry)
    allocate(List%Entry(n1))

    List%Entry(1:n0) = L0

    deallocate(L0)
  end subroutine PILIndex_Allocate

  !****************************************************************************
  !****s* PILIndex/PILIndex_DeAllocate
  ! NAME
  ! subroutine PILIndex_DeAllocate(List)
  ! PURPOSE
  ! DeAllocate all memory connected wit this index list.
  !****************************************************************************
  subroutine PILIndex_DeAllocate(List)

    type(tIndexList), intent(INOUT) :: List

    if (allocated(List%PartNumber)) deallocate(List%PartNumber)
    if (allocated(List%Entry)     ) deallocate(List%Entry)
    List%nEntry = 0
    List%nHole  = 0
  end subroutine PILIndex_DeAllocate

  !****************************************************************************
  !****f* PILIndex/PILIndex_DELETE
  ! NAME
  ! integer function PILIndex_DELETE(List, number)
  ! PURPOSE
  ! find the place, where the "number" is in the list,
  ! delete it and return the index for others list to follow.
  !
  ! INPUTS
  ! * type(tIndexList) :: List   -- The List
  ! * integer          :: number -- The particle number to find
  !
  ! OUTPUT
  ! if the return value R is > 0:
  ! here the entry "number" was found in the list
  !
  ! if the return value R <= 0:
  ! the entry was not found and not deleted
  !
  ! NOTES
  ! Attention, this routine returns iH=List%Entry(i1), i.e. directy the
  ! index for ValueList(...). (Contrary to the _FIND routine, where
  ! you have to calculate the index by looking it up in %Entry(...). But here
  ! this entry is already deleted!)
  !****************************************************************************
  integer function PILIndex_DELETE(List, number)

    type(tIndexList), intent(INOUT) :: List
    integer,          intent(IN) :: number

    integer :: i1, iH

    i1 = PILIndex_FIND(List, number)
    PILIndex_DELETE = i1
    if (i1 <= 0) return

    iH = List%Entry(i1)
    PILIndex_DELETE = iH

    !=== delete entry information by shifting the stuff above: ===
    List%PartNumber(i1:List%nEntry-1) = List%PartNumber(i1+1:List%nEntry)
    List%Entry(i1:List%nEntry-1) = List%Entry(i1+1:List%nEntry)
    List%nEntry = List%nEntry - 1
    List%nHole  = List%nHole  + 1

    ! a safe place to store the free slot:
    List%Entry(size(List%Entry)-List%nHole) = iH

    ! List totally cleaned up?
    if (List%nEntry==0) then
       List%nHole = 0 ! we may forget all the 'holes'
       write(*,*) 'list totally cleaned.'
    end if


  end function PILIndex_DELETE

end module PILIndex
