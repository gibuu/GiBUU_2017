!******************************************************************************
!****m* /PIL_omegaDec
! NAME
! module PIL_omegaDec
! PURPOSE
!
! Provide some storage method for additional information connected
! to a particle, characterized by its unique particleNumber.
!
! This is used in order to store the density at decay point for each
! omega decay (into pi0 gamma).
!
! NOTES
! * "PIL" stands for "PartInfoList"
!******************************************************************************
module PIL_omegaDec

  use PILIndex
  implicit none
  private

  !****************************************************************************
  !****t* PIL_omegaDec/decayInfo
  ! NAME
  ! type decayInfo
  ! PURPOSE
  ! a container of the information to be stored
  ! SOURCE
  type decayInfo
     real :: dens
  end type
  !****************************************************************************


  !****************************************************************************
  !****ig* PIL_omegaDec/IndexList
  ! SOURCE
  type(tIndexList), save :: IndexList
  ! PURPOSE
  ! The list, were the particle numbers are connected with the (physical)
  ! storage position
  !****************************************************************************


  !****************************************************************************
  !****ig* PIL_omegaDec/ValueList
  ! SOURCE
  type(decayInfo), save, allocatable :: ValueList(:)
  ! PURPOSE
  ! The list, were the information is stored
  !****************************************************************************


  public :: PIL_omegaDec_Put, PIL_omegaDec_Get
  public :: PIL_omegaDec_Deallocate
  public :: PIL_omegaDec_Zero

contains

  !****************************************************************************
  !****s* PIL_omegaDec/PIL_omegaDec_Deallocate
  ! NAME
  ! subroutine PIL_omegaDec_Deallocate
  ! PURPOSE
  ! Deallocate the memory for this list and the corresponding index list.
  !****************************************************************************
  subroutine PIL_omegaDec_Deallocate()
    call PILIndex_DeAllocate(IndexList)
    if (allocated(ValueList)) deallocate(ValueList)
  end subroutine


  !****************************************************************************
  !****s* PIL_omegaDec/PIL_omegaDec_Zero
  ! NAME
  ! subroutine PIL_omegaDec_Zero()
  ! PURPOSE
  ! Reset the list by setting the counter of stored information to 0.
  ! No allocation or deallocation of memory happens.
  !****************************************************************************
  subroutine PIL_omegaDec_Zero()
    IndexList%nEntry = 0
  end subroutine


  !****************************************************************************
  !****s* PIL_omegaDec/PIL_omegaDec_Put
  ! NAME
  ! subroutine PIL_omegaDec_Put (number, d)
  ! PURPOSE
  ! Store the information connected with particle "number" in the list.
  ! INPUTS
  ! * integer :: number -- the (unique) particle number
  ! * real    :: d      -- density at decay point
  ! OUTPUT
  ! none
  !****************************************************************************
  subroutine PIL_omegaDec_Put (number, d)
    integer, intent(in) :: number
    real, intent(in)    :: d

    integer :: iEntry

    iEntry = PILIndex_PUT(IndexList, number,"rho0Dec")

!    write(*,*) '###PUT: ',number,r,iEntry

    if (iEntry>0) then      ! everything is ok
       ValueList(iEntry)%dens = d
    else
       call PIL_omegaDec_Allocate() ! do (re)allocate
       ValueList(-iEntry)%dens = d
    end if
  end subroutine PIL_omegaDec_Put


  !****************************************************************************
  !****f* PIL_omegaDec/PIL_omegaDec_Get
  ! NAME
  ! logical function PIL_omegaDec_Get (number, d)
  ! PURPOSE
  ! Get the stored information of particle "number".
  ! INPUTS
  ! * integer :: number -- the (unique) particle number
  ! OUTPUT
  ! * real :: d -- the stored density
  ! * The (logical) return value signals, whether inforation about
  !   this particle was found in the list or not.
  !****************************************************************************
  logical function PIL_omegaDec_Get (number, d)
    integer, intent(in) :: number
    real, intent(out)   :: d

    integer :: iEntry

    iEntry = PILIndex_Find(IndexList,number)

    ! ATTENTION: iEntry is the line of information in the IndexList.
    ! The information connected to the particle is stored in the line
    ! IndexList%Entry(iEntry) in the array ValueList !!!!!

    if (iEntry > 0) then
       d = ValueList(IndexList%Entry(iEntry))%dens
       PIL_omegaDec_Get = .TRUE.
    else
       d = -1.
       PIL_omegaDec_Get = .FALSE.
    end if

  end function PIL_omegaDec_Get


  !****************************************************************************
  !****is* PIL_omegaDec/PIL_omegaDec_Allocate
  ! NAME
  ! subroutine PIL_omegaDec_Allocate
  ! PURPOSE
  ! Do the allocation and reallocation of the value vector.
  ! The new size is taken from the size of the IndexList vectors.
  ! NOTES
  ! For security one should insert here checks, whether the memory allocations
  ! failed and stop execution in these cases.
  !****************************************************************************
  subroutine PIL_omegaDec_Allocate
    integer :: n0, n1
    type(decayInfo), allocatable :: L0(:)

    n1 = size(IndexList%PartNumber) ! new size

    if (.not.allocated(ValueList)) then
       allocate(ValueList(n1))
       return
    end if

    n0 = size(ValueList)            ! old size

    allocate(L0(n0))
    L0(:) = ValueList(:)
    deallocate(ValueList)
    allocate(ValueList(n1))
    ValueList(1:n0) =  L0(1:n0)
    deallocate(L0)

  end subroutine PIL_omegaDec_Allocate


end module PIL_omegaDec
