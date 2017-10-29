!******************************************************************************
!****m* /PIL_mesonMom
! NAME
! module PIL_mesonMom
! PURPOSE
! Provide some storage method for a real valued information connected
! to a particle characterized by its unique particleNumber.
!
! This is closely connected to the module PILIndex
!
! This module has a doubled purpose:
! * It should do a practical work:
!   It stores information connected to the cross section evolution
!   in the hadronization process.
! * By copy'n'paste you can transfer this module to other information
!   storage tasks. The necessary changes are on a query'n'replace
!   level.
!
! The routines ...PUT and ...GET could be overwritten by routines using
! the explicit type definition of tValueEntry, which would make sense
! for more complicated informations.
!
! INPUTS
! (none)
!
! NOTES
! * "PIL" stands for "PartInfoList"
!******************************************************************************
module PIL_mesonMom

  use PILIndex

  implicit none
  private

  !****************************************************************************
  !****t* PIL_mesonMom/tValueEntry
  ! NAME
  ! type tValueEntry
  ! PURPOSE
  ! a container of the information to be stored
  ! SOURCE
  !
  type tValueEntry
     real, dimension(1:3) :: val
  end type
  !****************************************************************************


  !****************************************************************************
  !****g* PIL_mesonMom/IndexList
  ! PURPOSE
  ! The list, were the particle numbers are connected with the (physical)
  ! storage position
  !
  ! SOURCE
  !
  type(tIndexList), save :: IndexList
  !****************************************************************************


  !****************************************************************************
  !****g* PIL_mesonMom/ValueList
  ! PURPOSE
  ! The list, were the information is stored
  !
  ! SOURCE
  !
  type(tValueEntry), allocatable, save :: ValueList(:)
  !****************************************************************************


  public :: PIL_mesonMom_PUT, PIL_mesonMom_GET
  public :: PIL_mesonMom_DeAllocate
  public :: PIL_mesonMom_ZERO
!   PUBLIC :: PIL_mesonMom_Print

contains

  !****************************************************************************
  !****s* PIL_mesonMom/PIL_mesonMom_DeAlloc
  ! NAME
  ! subroutine PIL_mesonMom_DeAlloc
  ! PURPOSE
  ! Deallocate the memory for this list and the corresponding index list.
  !****************************************************************************
  subroutine PIL_mesonMom_DeAllocate()
    call PILIndex_DeAllocate(IndexList)
    if (allocated(ValueList)) deallocate(ValueList)
  end subroutine PIL_mesonMom_DeAllocate


  !****************************************************************************
  !****s* PIL_mesonMom/PIL_mesonMom_ZERO
  ! NAME
  ! subroutine PIL_mesonMom_ZERO()
  ! PURPOSE
  ! Reset the list by setting the counter of stored information to 0.
  ! No allocation or deallocation of memory happens.
  !****************************************************************************
  subroutine PIL_mesonMom_ZERO()
    IndexList%nEntry = 0
  end subroutine PIL_mesonMom_ZERO


  !****************************************************************************
  !****s* PIL_mesonMom/PIL_mesonMom_PUT
  ! NAME
  ! subroutine PIL_mesonMom_PUT(number,r)
  ! PURPOSE
  ! Store the information "r" connected with particle "number" in the list.
  ! INPUTS
  ! * integer :: number -- the (unique) particle number
  ! * real    :: r      -- the information to store
  ! OUTPUT
  ! none
  !****************************************************************************
  subroutine PIL_mesonMom_PUT(number,r)
    integer, intent(IN) :: number
    real, dimension(1:3),   intent(IN) :: r

    integer :: iEntry

    iEntry = PILIndex_PUT(IndexList, number,"mesonMom")

!    write(*,*) '###PUT: ',number,r,iEntry

    if (iEntry>0) then      ! everything is ok
       ValueList(iEntry)%val = r
    else
       call PIL_mesonMom_Allocate() ! do (re)allocate
       ValueList(-iEntry)%val = r
    end if
  end subroutine PIL_mesonMom_PUT


  !****************************************************************************
  !****f* PIL_mesonMom/PIL_mesonMom_GET
  ! NAME
  ! logical function PIL_mesonMom_GET(number,r)
  ! PURPOSE
  ! Get the stored information of particle "number"
  ! INPUTS
  ! * integer, intent(IN)  :: number -- the (unique) particle number
  ! OUTPUT
  ! * real,    intent(OUT) :: r -- the stored information (or 0.)
  ! * the (logical) return value signals, whether inforation about
  !   this particle was found in the list or not.
  !****************************************************************************
  logical function PIL_mesonMom_GET(number,r)
    integer, intent(IN)  :: number
    real, dimension(1:3),   intent(OUT) :: r

    integer :: iEntry

    iEntry = PILIndex_FIND(IndexList,number)

    ! ATTENTION: iEntry is the line of information in the IndexList.
    ! The information connected to the particle is stored in the line
    ! IndexList%Entry(iEntry) in the array ValueList !!!!!

    if (iEntry > 0) then
       r = ValueList(IndexList%Entry(iEntry))%val
       PIL_mesonMom_GET = .TRUE.
    else
       r = 0.0
       PIL_mesonMom_GET = .FALSE.
    end if

  end function PIL_mesonMom_GET


  !****************************************************************************
  !****is* PIL_mesonMom/PIL_mesonMom_Allocate
  ! NAME
  ! subroutine PIL_mesonMom_Allocate
  ! PURPOSE
  ! Do the allocation and reallocation of the value vector.
  ! The new size is taken from the size of the IndexList vectors.
  ! NOTES
  ! For security one should insert here checks, whether the memory allocations
  ! failed and stop execution in these cases.
  !****************************************************************************
  subroutine PIL_mesonMom_Allocate
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
    end do
    deallocate(ValueList)
    allocate(ValueList(n1))
    do i=1,n0
       ValueList(i)%val = L0(i)%val
    end do
    deallocate(L0)

  end subroutine PIL_mesonMom_Allocate


  !****************************************************************************
  !****s* PIL_mesonMom/PIL_mesonMom_Print
  ! NAME
  ! subroutine PIL_mesonMom_Print(file)
  ! PURPOSE
  ! Print the list to file
  !****************************************************************************
!   subroutine PIL_mesonMom_Print(file)
!     integer, intent(IN) :: file
!
!     integer :: i
!
!     write(file,*) '****** PIL_mesonMom:',IndexList%nEntry
!     do i=1,IndexList%nEntry
!        write(file,'(i8.0,i8.0,3g12.5)') i,IndexList%PartNumber(i),ValueList(IndexList%Entry(i))%val
!     end do
!
!   end subroutine PIL_mesonMom_Print



end module PIL_mesonMom
