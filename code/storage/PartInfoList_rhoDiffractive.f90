!******************************************************************************
!****m* /PIL_rhoDiffractive
! NAME
! module PIL_rhoDiffractive
! PURPOSE
! cf. module PIL_nLead
!
! Provide some storage method for a real valued information connected
! to a particle characterized by its unique particleNumber.
!
! This is closely connected to the module PILIndex
!
! NOTES
! * "PIL" stands for "PartInfoList"
!******************************************************************************
module PIL_rhoDiff

  use PILIndex

  !****************************************************************************
  !****t* PIL_rhoDiffractive/diffractiveInfo
  ! NAME
  ! type diffractiveInfo
  ! PURPOSE
  ! a container of the information to be stored
  ! SOURCE
  !
  type diffractiveInfo
     logical :: val
     real    :: valEpsR
     real,dimension(0:3) :: valRecoil
  end type diffractiveInfo
  !****************************************************************************

  private

  !****************************************************************************
  !****ig* PIL_rhoDiffractive/IndexList
  ! PURPOSE
  ! The list, were the particle numbers are connected with the (physical)
  ! storage position
  !
  ! SOURCE
  !
  type(tIndexList),save :: IndexList
  !****************************************************************************



  !****************************************************************************
  !****ig* PIL_rhoDiffractive/ValueList
  ! PURPOSE
  ! The list, were the information is stored
  !
  ! SOURCE
  !
  type(diffractiveInfo),save,  allocatable :: ValueList(:)
  !****************************************************************************


  public :: PIL_rhoDiffractive_PUT, PIL_rhoDiffractive_GET
  public :: PIL_rhoDiffractive_DeAllocate
  public :: PIL_rhoDiffractive_ZERO
!   PUBLIC :: PIL_rhoDiffractive_Print

contains
  !****************************************************************************
  !****s* PIL_rhoDiffractive/PIL_rhoDiffractive_DeAlloc
  ! NAME
  ! subroutine PIL_rhoDiffractive_DeAlloc
  ! PURPOSE
  ! Deallocate the memory for this list and the corresponding index list.
  !****************************************************************************
  subroutine PIL_rhoDiffractive_DeAllocate()
    call PILIndex_DeAllocate(IndexList)
    if (allocated(ValueList)) deallocate(ValueList)
  end subroutine PIL_rhoDiffractive_DeAllocate

  !****************************************************************************
  !****s* PIL_rhoDiffractive/PIL_rhoDiffractive_ZERO
  ! NAME
  ! subroutine PIL_rhoDiffractive_ZERO()
  ! PURPOSE
  ! Reset the list by setting the counter of stored information to 0.
  ! No allocation or deallocation of memory happens.
  !****************************************************************************
  subroutine PIL_rhoDiffractive_ZERO()
    IndexList%nEntry = 0
  end subroutine PIL_rhoDiffractive_ZERO


  !****************************************************************************
  !****s* PIL_rhoDiffractive/PIL_rhoDiffractive_PUT
  ! NAME
  ! subroutine PIL_rhoDiffractive_PUT(number,r,epsR,recoil)
  ! PURPOSE
  ! Store the information "r" connected with particle "number" in the list.
  ! INPUTS
  ! * integer :: number -- the (unique) particle number
  ! * logical    :: r      -- the information to store
  ! * real :: epsR -- epsilon * R
  ! * real, dimension(0:3) :: recoil -- momentum vector of target recoil
  ! OUTPUT
  ! none
  !****************************************************************************
  subroutine PIL_rhoDiffractive_PUT(number,r,epsR,recoil)
    implicit none
    integer, intent(IN) :: number
    logical, intent(IN) :: r
    real, intent(IN) :: epsR
    real, dimension(0:3),intent(IN) :: recoil

    integer :: iEntry

    iEntry = PILIndex_PUT(IndexList, number,"rhoDiffractive")

!    write(*,*) '###PUT: ',number,r,iEntry

    if (iEntry>0) then      ! everything is ok
       ValueList(iEntry)%val = r
       ValueList(iEntry)%valEpsR = epsR
       ValueList(iEntry)%valRecoil = recoil
    else
       call PIL_rhoDiffractive_Allocate() ! do (re)allocate
       ValueList(-iEntry)%val = r
       ValueList(-iEntry)%valEpsR = epsR
       ValueList(-iEntry)%valRecoil = recoil
    end if
  end subroutine PIL_rhoDiffractive_PUT


  !****************************************************************************
  !****f* PIL_rhoDiffractive/PIL_rhoDiffractive_GET
  ! NAME
  ! logical function PIL_rhoDiffractive_GET(number,r,epsR,recoil)
  ! PURPOSE
  ! Get the stored information of particle "number"
  ! INPUTS
  ! * integer  :: number -- the (unique) particle number
  ! OUTPUT
  ! * logical :: r -- the stored information (or .false.)
  ! * real, OPTIONAL :: epsR -- epsilon * R
  ! * real, dimension(0:3), OPTIONAL :: recoil -- momentum vector of target recoil
  ! * the (logical) return value signals, whether inforation about
  !   this particle was found in the list or not.
  !****************************************************************************
  logical function PIL_rhoDiffractive_GET(number,r,epsR,recoil)
    implicit none
    integer, intent(IN)  :: number
    logical, intent(OUT) :: r
    real, intent(OUT), optional :: epsR
    real, dimension(0:3),intent(OUT),OPTIONAL :: recoil

    integer :: iEntry

    iEntry = PILIndex_FIND(IndexList,number)

    ! ATTENTION: iEntry is the line of information in the IndexList.
    ! The information connected to the particle is stored in the line
    ! IndexList%Entry(iEntry) in the array ValueList !!!!!

    if (iEntry > 0) then
       r = ValueList(IndexList%Entry(iEntry))%val
       if (present(recoil)) recoil = ValueList(IndexList%Entry(iEntry))%valRecoil
       if (present(epsR)) epsR = ValueList(IndexList%Entry(iEntry))%valEpsR
       PIL_rhoDiffractive_GET = .TRUE.
    else
       r = .false.
       if (present(recoil)) recoil = 0.0
       if (present(epsR)) epsR = 0.0
       PIL_rhoDiffractive_GET = .FALSE.
    end if

  end function PIL_rhoDiffractive_GET


  !****************************************************************************
  !****is* PIL_rhoDiffractive/PIL_rhoDiffractive_Allocate
  ! NAME
  ! subroutine PIL_rhoDiffractive_Allocate
  ! PURPOSE
  ! Do the allocation and reallocation of the value vector.
  ! The new size is taken from the size of the IndexList vectors.
  ! NOTES
  ! For security one should insert here checks, whether the memory allocations
  ! failed and stop execution in these cases.
  !****************************************************************************
  subroutine PIL_rhoDiffractive_Allocate
    implicit none

    integer :: n0, n1
    type(diffractiveInfo), allocatable :: L0(:)

    n1 = size(IndexList%PartNumber) ! new size

    if (.not.allocated(ValueList)) then
       allocate(ValueList(n1))
       return
    end if

    n0 = size(ValueList)            ! old size

    allocate(L0(n0))
    L0(:)%val = ValueList(:)%val
    L0(:)%valEpsR = ValueList(:)%valEpsR
    L0(:)%valRecoil(0) = ValueList(:)%valRecoil(0)
    L0(:)%valRecoil(1) = ValueList(:)%valRecoil(1)
    L0(:)%valRecoil(2) = ValueList(:)%valRecoil(2)
    L0(:)%valRecoil(3) = ValueList(:)%valRecoil(3)
    deallocate(ValueList)
    allocate(ValueList(n1))
    ValueList(1:n0)%val = L0(1:n0)%val
    ValueList(1:n0)%valEpsR = L0(1:n0)%valEpsR
    ValueList(1:n0)%valRecoil(0) = L0(1:n0)%valRecoil(0)
    ValueList(1:n0)%valRecoil(1) = L0(1:n0)%valRecoil(1)
    ValueList(1:n0)%valRecoil(2) = L0(1:n0)%valRecoil(2)
    ValueList(1:n0)%valRecoil(3) = L0(1:n0)%valRecoil(3)

    deallocate(L0)

  end subroutine PIL_rhoDiffractive_Allocate


  !****************************************************************************
  !****is* PIL_rhoDiffractive/PIL_rhoDiffractive_Print
  ! NAME
  ! subroutine PIL_rhoDiffractive_Print(file)
  ! PURPOSE
  ! Print the list to file
  !****************************************************************************
!   subroutine PIL_rhoDiffractive_Print(file)
!     implicit none
!     integer, intent(IN) :: file
!
!     integer :: i
!
!     write(file,*) '****** PIL_rhoDiffractive:',IndexList%nEntry
!     do i=1,IndexList%nEntry
!        write(file,'(i8.0,i8.0,g12.5)') i,IndexList%PartNumber(i),ValueList(IndexList%Entry(i))%val
!     end do
!
!   end subroutine PIL_rhoDiffractive_Print

end module PIL_rhoDiff
