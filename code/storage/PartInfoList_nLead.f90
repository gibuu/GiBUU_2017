!******************************************************************************
!****m* /PIL_nLead
! NAME
! module PIL_nLead
! PURPOSE
! Provide some storage method for a real valued information connected
! to a particle characterized by its unique particleNumber.
!
! This is closely connected to the module PILIndex
!
! The actual stored value is the scaled cross section divided by the hard
! scale Q2.
!
! INPUTS
! (none)
!
! NOTES
! * "PIL" stands for "PartInfoList"
! * unfortunately this module does NOT store the number of leading quarks,
!   as the name could suggest; for this purpose see PIL_FormInfo.
!   The neame of this module should better be "..._Pedest..."
!******************************************************************************
module PIL_nLead

  use PILIndex

  implicit none
  private

  !****************************************************************************
  !****t* PIL_nLead/tValueEntry
  ! NAME
  ! type tValueEntry
  ! PURPOSE
  ! a container of the information to be stored
  ! SOURCE
  !
  type tValueEntry
     real :: val
  end type
  !****************************************************************************


  !****************************************************************************
  !****g* PIL_nLead/IndexList
  ! PURPOSE
  ! The list, were the particle numbers are connected with the (physical)
  ! storage position
  !
  ! SOURCE
  !
  type(tIndexList), save :: IndexList
  !****************************************************************************


  !****************************************************************************
  !****g* PIL_nLead/ValueList
  ! PURPOSE
  ! The list, were the information is stored
  !
  ! SOURCE
  !
  type(tValueEntry), allocatable, save :: ValueList(:)
  !****************************************************************************


  public :: PIL_nLead_PUT, PIL_nLead_GET
  public :: PIL_nLead_DeAllocate
  public :: PIL_nLead_ZERO
!   PUBLIC :: PIL_nLead_Print

contains

  !****************************************************************************
  !****s* PIL_nLead/PIL_nLead_DeAlloc
  ! NAME
  ! subroutine PIL_nLead_DeAlloc
  ! PURPOSE
  ! Deallocate the memory for this list and the corresponding index list.
  !****************************************************************************
  subroutine PIL_nLead_DeAllocate()
    call PILIndex_DeAllocate(IndexList)
    if (allocated(ValueList)) deallocate(ValueList)
  end subroutine PIL_nLead_DeAllocate


  !****************************************************************************
  !****s* PIL_nLead/PIL_nLead_ZERO
  ! NAME
  ! subroutine PIL_nLead_ZERO()
  ! PURPOSE
  ! Reset the list by setting the counter of stored information to 0.
  ! No allocation or deallocation of memory happens.
  !****************************************************************************
  subroutine PIL_nLead_ZERO()
    IndexList%nEntry = 0
  end subroutine PIL_nLead_ZERO


  !****************************************************************************
  !****s* PIL_nLead/PIL_nLead_PUT
  ! NAME
  ! subroutine PIL_nLead_PUT(number,r)
  ! PURPOSE
  ! Store the information "r" connected with particle "number" in the list.
  ! INPUTS
  ! * integer :: number -- the (unique) particle number
  ! * real    :: r      -- the information to store
  ! OUTPUT
  ! none
  !****************************************************************************
  subroutine PIL_nLead_PUT(number,r)
    integer, intent(IN) :: number
    real,    intent(IN) :: r

    integer :: iEntry

    iEntry = PILIndex_PUT(IndexList, number,"nLead")

!    write(*,*) '###PUT: ',number,r,iEntry

    if (iEntry>0) then      ! everything is ok
       ValueList(iEntry)%val = r
    else
       call PIL_nLead_Allocate() ! do (re)allocate
       ValueList(-iEntry)%val = r
    end if
  end subroutine PIL_nLead_PUT


  !****************************************************************************
  !****f* PIL_nLead/PIL_nLead_GET
  ! NAME
  ! logical function PIL_nLead_GET(number,r)
  ! PURPOSE
  ! Get the stored information of particle "number"
  ! INPUTS
  ! * integer, intent(IN)  :: number -- the (unique) particle number
  ! OUTPUT
  ! * real,    intent(OUT) :: r -- the stored information (or 0.)
  ! * the (logical) return value signals, whether inforation about
  !   this particle was found in the list or not.
  !****************************************************************************
  logical function PIL_nLead_GET(number,r)
    integer, intent(IN)  :: number
    real,    intent(OUT) :: r

    integer :: iEntry

    iEntry = PILIndex_FIND(IndexList,number)

    ! ATTENTION: iEntry is the line of information in the IndexList.
    ! The information connected to the particle is stored in the line
    ! IndexList%Entry(iEntry) in the array ValueList !!!!!

    if (iEntry > 0) then
       r = ValueList(IndexList%Entry(iEntry))%val
       PIL_nLead_GET = .TRUE.
    else
       r = 0.0
       PIL_nLead_GET = .FALSE.
    end if

  end function PIL_nLead_GET


  !****************************************************************************
  !****is* PIL_nLead/PIL_nLead_Allocate
  ! NAME
  ! subroutine PIL_nLead_Allocate
  ! PURPOSE
  ! Do the allocation and reallocation of the value vector.
  ! The new size is taken from the size of the IndexList vectors.
  ! NOTES
  ! For security one should insert here checks, whether the memory allocations
  ! failed and stop execution in these cases.
  !****************************************************************************
  subroutine PIL_nLead_Allocate
    integer :: n0, n1
    type(tValueEntry), allocatable :: L0(:)

    n1 = size(IndexList%PartNumber) ! new size

    if (.not.allocated(ValueList)) then
       allocate(ValueList(n1))
       return
    end if

    n0 = size(ValueList)            ! old size

    allocate(L0(n0))
    L0(:)%val = ValueList(:)%val
    deallocate(ValueList)
    allocate(ValueList(n1))
    ValueList(1:n0)%val = L0(1:n0)%val
    deallocate(L0)

  end subroutine PIL_nLead_Allocate


  !****************************************************************************
  !****is* PIL_nLead/PIL_nLead_Print
  ! NAME
  ! subroutine PIL_nLead_Print(file)
  ! PURPOSE
  ! Print the list to file
  !****************************************************************************
!   subroutine PIL_nLead_Print(file)
!     integer, intent(IN) :: file
!
!     integer :: i
!
!     write(file,*) '****** PIL_nLead:',IndexList%nEntry
!     do i=1,IndexList%nEntry
!        write(file,'(i8.0,i8.0,g12.5)') i,IndexList%PartNumber(i),ValueList(IndexList%Entry(i))%val
!     end do
!
!   end subroutine PIL_nLead_Print


end module PIL_nLead
