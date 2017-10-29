!******************************************************************************
!****m* /PIL_FormInfo
! NAME
! module PIL_FormInfo
! PURPOSE
! cf. module PIL_nLead
!
! Provide some storage method for a integer valued information connected
! to a particle characterized by its unique particleNumber.
!
! This is closely connected to the module PILIndex
!
! The stored information is a flag connected with string fragmentation:
!   EArr(2,...) + 10 * EArr(3,...) + 100 * nLead
! (cf. GetJetSetVec.f)
!
! EArr(1,...):
! * 0: hadron has not been processed
! * 1: everything should be ok (in principle)
! * 3: next to last hadron in fragmentation
! * 4: last hadron in fragmentation
!
! EArr(2,...):
! * 1: problems with ProdTime 1 ; additive
! * 2: problems with ProdTime 2 ; additive
! * 4: problems with FormTime   ; additive
! Values can be tested by the intrinsic f77-function IAnd(...,1|2|4).
!
! EArr(3,...):
! * 1: hadron from QQ-String
! * 2: hadron from QQ-String with internal Gluons
! * 3: hadron from pure gluonic GG-string
! * 4: hadron from Cluster-Decay -> 1 hadron
! * 5: hadron from Cluster-Decay -> 2 hadrons
! * 6: hadron from Doku-Line
!
! EArr(4,...): number of the string
!
! EArr(5,...): Rank of the particle, like it was calculated
!
! EArr(6,...): rank as the minmal number of
! particles to the left or right (plus 1). (reverse ordering is considered)
!
! NOTES
! * "PIL" stands for "PartInfoList"
!******************************************************************************
module PIL_FormInfo

  use PILIndex
  implicit none
  private

  !****************************************************************************
  !****t* PIL_FormInfo/FormationVals
  ! NAME
  ! type FormationVals
  ! PURPOSE
  ! a container of the information to be stored
  ! SOURCE
  !
  type FormationVals
     integer :: val
  end type FormationVals
  !****************************************************************************


  !****************************************************************************
  !****ig* PIL_FormInfo/IndexList
  ! PURPOSE
  ! The list, were the particle numbers are connected with the (physical)
  ! storage position
  !
  ! SOURCE
  !
  type(tIndexList),save :: IndexList
  !****************************************************************************


  !****************************************************************************
  !****ig* PIL_FormInfo/ValueList
  ! PURPOSE
  ! The list, were the information is stored
  !
  ! SOURCE
  !
  type(FormationVals),save,  allocatable :: ValueList(:)
  !****************************************************************************


  public :: PIL_FormInfo_PUT, PIL_FormInfo_GET
  public :: PIL_FormInfo_DeAllocate
  public :: PIL_FormInfo_ZERO
!   PUBLIC :: PIL_FormInfo_Print

contains

  !****************************************************************************
  !****s* PIL_FormInfo/PIL_FormInfo_DeAlloc
  ! NAME
  ! subroutine PIL_FormInfo_DeAlloc
  ! PURPOSE
  ! Deallocate the memory for this list and the corresponding index list.
  !****************************************************************************
  subroutine PIL_FormInfo_DeAllocate()
    call PILIndex_DeAllocate(IndexList)
    if (allocated(ValueList)) deallocate(ValueList)
  end subroutine PIL_FormInfo_DeAllocate

  !****************************************************************************
  !****s* PIL_FormInfo/PIL_FormInfo_ZERO
  ! NAME
  ! subroutine PIL_FormInfo_ZERO()
  ! PURPOSE
  ! Reset the list by setting the counter of stored information to 0.
  ! No allocation or deallocation of memory happens.
  !****************************************************************************
  subroutine PIL_FormInfo_ZERO()
    IndexList%nEntry = 0
  end subroutine PIL_FormInfo_ZERO


  !****************************************************************************
  !****s* PIL_FormInfo/PIL_FormInfo_PUT
  ! NAME
  ! subroutine PIL_FormInfo_PUT(number,r)
  ! PURPOSE
  ! Store the information "r" connected with particle "number" in the list.
  ! INPUTS
  ! * integer :: number -- the (unique) particle number
  ! * integer    :: r      -- the information to store
  ! OUTPUT
  ! none
  !****************************************************************************
  subroutine PIL_FormInfo_PUT(number,r)
    integer, intent(IN) :: number
    integer, intent(IN) :: r

    integer :: iEntry

    iEntry = PILIndex_PUT(IndexList, number,"FormInfo")

!    write(*,*) '###PUT: ',number,r,iEntry

    if (iEntry>0) then      ! everything is ok
       ValueList(iEntry)%val = r
    else
       call PIL_FormInfo_Allocate() ! do (re)allocate
       ValueList(-iEntry)%val = r
    end if
  end subroutine PIL_FormInfo_PUT


  !****************************************************************************
  !****f* PIL_FormInfo/PIL_FormInfo_GET
  ! NAME
  ! logical function PIL_FormInfo_GET(number,r)
  ! PURPOSE
  ! Get the stored information of particle "number"
  ! INPUTS
  ! * integer, intent(IN)  :: number -- the (unique) particle number
  ! OUTPUT
  ! * integer,    intent(OUT) :: r -- the stored information (or .false.)
  ! * the (logical) return value signals, whether inforation about
  !   this particle was found in the list or not.
  !****************************************************************************
  logical function PIL_FormInfo_GET(number,r)
    implicit none
    integer, intent(IN)  :: number
    integer, intent(OUT) :: r

    integer :: iEntry

    iEntry = PILIndex_FIND(IndexList,number)

    ! ATTENTION: iEntry is the line of information in the IndexList.
    ! The information connected to the particle is stored in the line
    ! IndexList%Entry(iEntry) in the array ValueList !!!!!

    if (iEntry > 0) then
       r = ValueList(IndexList%Entry(iEntry))%val
       PIL_FormInfo_GET = .TRUE.
    else
       r = 0
       PIL_FormInfo_GET = .FALSE.
    end if

  end function PIL_FormInfo_GET


  !****************************************************************************
  !****is* PIL_FormInfo/PIL_FormInfo_Allocate
  ! NAME
  ! subroutine PIL_FormInfo_Allocate
  ! PURPOSE
  ! Do the allocation and reallocation of the value vector.
  ! The new size is taken from the size of the IndexList vectors.
  ! NOTES
  ! For security one should insert here checks, whether the memory allocations
  ! failed and stop execution in these cases.
  !****************************************************************************
  subroutine PIL_FormInfo_Allocate
    implicit none

    integer :: n0, n1
    type(FormationVals), allocatable :: L0(:)

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

  end subroutine PIL_FormInfo_Allocate


  !****************************************************************************
  !****is* PIL_FormInfo/PIL_FormInfo_Print
  ! NAME
  ! subroutine PIL_FormInfo_Print(file)
  ! PURPOSE
  ! Print the list to file
  !****************************************************************************
!   subroutine PIL_FormInfo_Print(file)
!     implicit none
!     integer, intent(IN) :: file
!
!     integer :: i
!
!     write(file,*) '****** PIL_FormInfo:',IndexList%nEntry
!     do i=1,IndexList%nEntry
!        write(file,'(i8.0,i8.0,g12.5)') i,IndexList%PartNumber(i),ValueList(IndexList%Entry(i))%val
!     end do
!
!   end subroutine PIL_FormInfo_Print

end module PIL_FormInfo
