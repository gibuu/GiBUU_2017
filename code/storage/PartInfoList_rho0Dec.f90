!******************************************************************************
!****m* /PIL_rho0Dec
! NAME
! module PIL_rho0Dec
! PURPOSE
! cf. module PIL_nLead
!
! Provide some storage method for a integer valued information connected
! to a particle characterized by its unique particleNumber.
!
! This is used in order to store for every charged pion the unique
! number of the opposite charged pion in a rho0 decay.
!
! This is closely connected to the module PILIndex
!
! NOTES
! * "PIL" stands for "PartInfoList"
!******************************************************************************
module PIL_rho0Dec

  use PILIndex
  implicit none
  private

  !****************************************************************************
  !****t* PIL_rho0Dec/decayInfo
  ! NAME
  ! type decayInfo
  ! PURPOSE
  ! a container of the information to be stored
  ! SOURCE
  !
  type decayInfo
     integer :: val
     real    :: tC, tP, tF
  end type decayInfo
  !****************************************************************************


  !****************************************************************************
  !****ig* PIL_rho0Dec/IndexList
  ! SOURCE
  !
  type(tIndexList),save :: IndexList
  ! PURPOSE
  ! The list, were the particle numbers are connected with the (physical)
  ! storage position
  !****************************************************************************


  !****************************************************************************
  !****ig* PIL_rho0Dec/ValueList
  ! SOURCE
  !
  type(decayInfo),save,  allocatable :: ValueList(:)
  ! PURPOSE
  ! The list, were the information is stored
  !****************************************************************************


  public :: PIL_rho0Dec_PUT, PIL_rho0Dec_GET
  public :: PIL_rho0Dec_DeAllocate
  public :: PIL_rho0Dec_ZERO
!   PUBLIC :: PIL_rho0Dec_Print

contains

  !****************************************************************************
  !****s* PIL_rho0Dec/PIL_rho0Dec_DeAlloc
  ! NAME
  ! subroutine PIL_rho0Dec_DeAlloc
  ! PURPOSE
  ! Deallocate the memory for this list and the corresponding index list.
  !****************************************************************************
  subroutine PIL_rho0Dec_DeAllocate()
    call PILIndex_DeAllocate(IndexList)
    if (allocated(ValueList)) deallocate(ValueList)
  end subroutine PIL_rho0Dec_DeAllocate

  !****************************************************************************
  !****s* PIL_rho0Dec/PIL_rho0Dec_ZERO
  ! NAME
  ! subroutine PIL_rho0Dec_ZERO()
  ! PURPOSE
  ! Reset the list by setting the counter of stored information to 0.
  ! No allocation or deallocation of memory happens.
  !****************************************************************************
  subroutine PIL_rho0Dec_ZERO()
    IndexList%nEntry = 0
  end subroutine PIL_rho0Dec_ZERO


  !****************************************************************************
  !****s* PIL_rho0Dec/PIL_rho0Dec_PUT
  ! NAME
  ! subroutine PIL_rho0Dec_PUT(number,r,tC,tP,tF)
  ! PURPOSE
  ! Store the information  connected with particle "number" in the list.
  ! INPUTS
  ! * integer :: number -- the (unique) particle number
  ! * integer :: r      -- the information to store
  ! * real    :: tC,tP,tF  -- collision, production and formation time
  ! OUTPUT
  ! none
  !****************************************************************************
  subroutine PIL_rho0Dec_PUT(number,r,tC,tP,tF)

    integer, intent(IN) :: number
    integer, intent(IN) :: r
    real, intent(IN)    :: tC,tP,tF

    integer :: iEntry

    iEntry = PILIndex_PUT(IndexList, number,"rho0Dec")

!    write(*,*) '###PUT: ',number,r,iEntry

    if (iEntry>0) then      ! everything is ok
       ValueList(iEntry)%val = r
       ValueList(iEntry)%tC  = tC
       ValueList(iEntry)%tP  = tP
       ValueList(iEntry)%tF  = tF
    else
       call PIL_rho0Dec_Allocate() ! do (re)allocate
       ValueList(-iEntry)%val = r
       ValueList(-iEntry)%tC  = tC
       ValueList(-iEntry)%tP  = tP
       ValueList(-iEntry)%tF  = tF
    end if
  end subroutine PIL_rho0Dec_PUT


  !****************************************************************************
  !****f* PIL_rho0Dec/PIL_rho0Dec_GET
  ! NAME
  ! logical function PIL_rho0Dec_GET(number,r,tC,tP,tF)
  ! PURPOSE
  ! Get the stored information of particle "number"
  ! INPUTS
  ! * integer  :: number -- the (unique) particle number
  ! OUTPUT
  ! * integer         :: r -- the stored information (or 0)
  ! * real, OPTIONAL  :: tC,tP,tF -- the stored collsion,production and
  !   formation times
  ! * the (logical) return value signals, whether inforation about
  !   this particle was found in the list or not.
  !****************************************************************************
  logical function PIL_rho0Dec_GET(number,r,tC,tP,tF)
    integer, intent(IN)  :: number
    integer, intent(OUT) :: r
    real, intent(OUT),OPTIONAL    :: tC,tP,tF

    integer :: iEntry

    iEntry = PILIndex_FIND(IndexList,number)

    ! ATTENTION: iEntry is the line of information in the IndexList.
    ! The information connected to the particle is stored in the line
    ! IndexList%Entry(iEntry) in the array ValueList !!!!!

    if (iEntry > 0) then
       r = ValueList(IndexList%Entry(iEntry))%val
       if (present(tC)) tC = ValueList(IndexList%Entry(iEntry))%tC
       if (present(tP)) tP = ValueList(IndexList%Entry(iEntry))%tP
       if (present(tF)) tF = ValueList(IndexList%Entry(iEntry))%tF
       PIL_rho0Dec_GET = .TRUE.
    else
       r = 0
       if (present(tC)) tC = -99.9
       if (present(tP)) tP = -99.9
       if (present(tF)) tF = -99.9
       PIL_rho0Dec_GET = .FALSE.
    end if

  end function PIL_rho0Dec_GET


  !****************************************************************************
  !****is* PIL_rho0Dec/PIL_rho0Dec_Allocate
  ! NAME
  ! subroutine PIL_rho0Dec_Allocate
  ! PURPOSE
  ! Do the allocation and reallocation of the value vector.
  ! The new size is taken from the size of the IndexList vectors.
  ! NOTES
  ! For security one should insert here checks, whether the memory allocations
  ! failed and stop execution in these cases.
  !****************************************************************************
  subroutine PIL_rho0Dec_Allocate
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

  end subroutine PIL_rho0Dec_Allocate


  !****************************************************************************
  !****is* PIL_rho0Dec/PIL_rho0Dec_Print
  ! NAME
  ! subroutine PIL_rho0Dec_Print(file)
  ! PURPOSE
  ! Print the list to file
  !****************************************************************************
!   subroutine PIL_rho0Dec_Print(file)
!     integer, intent(IN) :: file
!
!     integer :: i
!
!     write(file,*) '****** PIL_rho0Dec:',IndexList%nEntry
!     do i=1,IndexList%nEntry
!        write(file,'(i8.0,i8.0,i8, 3f11.5)') i,IndexList%PartNumber(i),&
!             & ValueList(IndexList%Entry(i))%val, &
!             & ValueList(IndexList%Entry(i))%tC, &
!             & ValueList(IndexList%Entry(i))%tP, &
!             & ValueList(IndexList%Entry(i))%tF
!     end do
!
!   end subroutine PIL_rho0Dec_Print

end module PIL_rho0Dec
