!******************************************************************************
!****m* /MassAssInfoDefinition
! NAME
! module MassAssInfoDefinition
!
! PURPOSE
! Provide information necessary for an improved mass selection according
! a relativistic Breit Wigner with mass dependent width. Also for in-medium
! widths
!******************************************************************************
module MassAssInfoDefinition

  implicit none
  private

  !****************************************************************************
  !****m* MassAssInfoDefinition/tMassAssInfo
  ! SOURCE
  !
  type, public :: tMassAssInfo
     logical :: IsStable
     logical :: IsMomDep
     real :: Mass0
     real :: Gamma0
     real,dimension(:), allocatable :: W ! = weight
     real,dimension(:), allocatable :: Q ! = Qmax
     real,dimension(:), allocatable :: M ! = M_min, M_max
     real,dimension(:), allocatable :: Y ! = Y_min, Y_max
  end type tMassAssInfo
  !
  ! PURPOSE
  ! store the info concerning the mass slices
  !
  ! The stored information is:
  ! * IsStable: flag if particle is stable
  ! * IsMomDep: flag if width is momentum dependent
  ! * Mass0,Gamma0: pole mass and width at the pole mass
  ! * Q: The maximum of the rejection weight for each bin
  !   (This is the ratio of the relativistic Breit-Wigner over the
  !   non-relativistic one)
  ! * W: The weight the selected bin has to be chosen
  ! * M: The mass boundaries of every bin
  ! * Y: The value 2*atan(2*(M-M0)/G) at the boundaries of every bin
  !****************************************************************************

  !****************************************************************************
  !****t* MassAssInfoDefinition/tMassAssInfoArray
  ! SOURCE
  !
  type, public :: tMassAssInfoArray
     real,    dimension(:,:), allocatable :: W ! = weight
     real,    dimension(:,:), allocatable :: Q ! = Qmax
     real,    dimension(:,:), allocatable :: M ! = M_min, M_max
     real,    dimension(:,:), allocatable :: Y ! = Y_min, Y_max
     integer, dimension(:),   allocatable :: n1, n2 ! = used bins
  end type tMassAssInfoArray
  !
  ! PURPOSE
  ! The storage of density dependent information in the meson and baryon
  ! modules.
  !
  ! see tMassAssInfo for further information
  !****************************************************************************

  !****************************************************************************
  !****g* MassAssInfoDefinition/UseMassAssInfo
  ! SOURCE
  !
  logical, save :: UseMassAssInfo  = .true.
  ! PURPOSE
  ! This switch indicates, whether we want to use the whole MassAssInfo
  ! machinery or stick to the old prescription of mass assignment.
  !
  ! You may set this switch via the jobcard. Anyhow, if your selection
  ! of switches for baryon and medium switches leads to cases which are
  ! not yet implemented, this flag is set to false automatically.
  !****************************************************************************


  logical, save :: initFlag=.true.

  public :: ResetMassAssInfo
  public :: SetMassAssInfo_Stable
  public :: SetMassAssInfo_FromArray
  public :: SetUpperBin
  public :: MassAssInfoQ
  public :: Get_UseMassAssInfo
  public :: Set_UseMassAssInfo

contains

  !****************************************************************************
  !****s* MassAssInfoDefinition/ResetMassAssInfo
  ! NAME
  ! subroutine ResetMassAssInfo(MAI)
  ! PURPOSE
  ! Reset the fields of the type
  ! INPUTS
  ! * type(tMassAssInfo) :: MAI -- the Info to reset
  ! OUTPUT
  ! * type(tMassAssInfo) :: MAI -- the Info resetted
  !****************************************************************************
  subroutine ResetMassAssInfo(MAI)
    type(tMassAssInfo), intent(inout) :: MAI

    MAI%IsStable = .false.
    MAI%IsMomDep = .false.
    MAI%Mass0    = 0.0
    MAI%Gamma0   = 0.0
    if (allocated(MAI%W)) deallocate(MAI%W)
    if (allocated(MAI%Q)) deallocate(MAI%Q)
    if (allocated(MAI%M)) deallocate(MAI%M)
    if (allocated(MAI%Y)) deallocate(MAI%Y)

  end subroutine ResetMassAssInfo

  !****************************************************************************
  !****s* MassAssInfoDefinition/SetMassAssInfo_Stable
  ! NAME
  ! subroutine SetMassAssInfo_Stable(MAI)
  ! PURPOSE
  ! Ensure correct init, if the particle is considered stable
  ! INPUTS
  ! * MAI%Mass0 and MAI%Gamma0 have to be set
  ! OUTPUT
  ! * MAI
  !****************************************************************************
  subroutine SetMassAssInfo_Stable(MAI)

    type(tMassAssInfo), intent(inout) :: MAI
    integer, parameter :: nB = 1

    if (allocated(MAI%W)) deallocate(MAI%W)
    if (allocated(MAI%Q)) deallocate(MAI%Q)
    if (allocated(MAI%M)) deallocate(MAI%M)
    if (allocated(MAI%Y)) deallocate(MAI%Y)
    allocate(MAI%W(nB))
    allocate(MAI%Q(nB))
    allocate(MAI%M(nB+1))
    allocate(MAI%Y(nB+1))

    MAI%IsStable = .true.
    MAI%IsMomDep = .false.
    MAI%W = 1.0
    MAI%Q = 1.0
    MAI%M = MAI%Mass0
    MAI%Y = 0.0

  end subroutine SetMassAssInfo_Stable

  !****************************************************************************
  !****s* MassAssInfoDefinition/SetMassAssInfo_FromArray
  ! NAME
  ! subroutine SetMassAssInfo_FromArray(MAI, MAIarr, iG, mix)
  ! PURPOSE
  ! Set the information 'MAI' by taken the values out of a stored array.
  ! INPUTS
  ! * type(tMassAssInfo) :: MAI -- the Info to be set
  ! * type(tMassAssInfoArray) :: MAIarr -- the array where the info has to
  !   be taken from
  ! * integer :: iG -- which bin for the width to use
  ! * real :: mixG -- weight for interpolating with next width bin
  ! OUTPUT
  ! * type(tMassAssInfo) :: MAI -- the Info
  ! NOTES
  ! at the moment, the interpolation is implemented to be linear in all
  ! components
  !****************************************************************************
  subroutine SetMassAssInfo_FromArray(MAI, MAIarr, iG, mixG)
    use CALLSTACK, only: TRACEBACK

    type(tMassAssInfo), intent(inout) :: MAI
    type(tMassAssInfoArray), intent(in) :: MAIarr
    integer, intent(in) :: iG
    real, intent(in) :: mixG

    integer :: iB1,iB2,nB

    if (.not.Get_UseMassAssInfo()) then
       call TRACEBACK('arrays can not be initialized')
    end if

    if (allocated(MAI%W)) deallocate(MAI%W)
    if (allocated(MAI%Q)) deallocate(MAI%Q)
    if (allocated(MAI%M)) deallocate(MAI%M)
    if (allocated(MAI%Y)) deallocate(MAI%Y)

    if (mixG>0.0) then ! we mix this bin and the next one
       iB1 = min(MAIarr%n1(iG),MAIarr%n1(iG+1))
       iB2 = max(MAIarr%n2(iG),MAIarr%n2(iG+1))
    else
       iB1 = MAIarr%n1(iG)
       iB2 = MAIarr%n2(iG)
    end if

    nB = iB2-iB1 + 1
    allocate(MAI%W(nB))
    allocate(MAI%Q(nB))
    allocate(MAI%M(nB+1))
    allocate(MAI%Y(nB+1))

    if (mixG>0.0) then
       ! linear interpolation of all values:

       MAI%W = (1-mixG)*MAIarr%W(iG,iB1:iB2) + mixG*MAIarr%W(iG+1,iB1:iB2)
       MAI%Q = (1-mixG)*MAIarr%Q(iG,iB1:iB2) + mixG*MAIarr%Q(iG+1,iB1:iB2)
       MAI%M = (1-mixG)*MAIarr%M(iG,iB1:iB2+1) + mixG*MAIarr%M(iG+1,iB1:iB2+1)
       MAI%Y = (1-mixG)*MAIarr%Y(iG,iB1:iB2+1) + mixG*MAIarr%Y(iG+1,iB1:iB2+1)
    else
       MAI%W = MAIarr%W(iG,iB1:iB2)
       MAI%Q = MAIarr%Q(iG,iB1:iB2)
       MAI%M = MAIarr%M(iG,iB1:iB2+1)
       MAI%Y = MAIarr%Y(iG,iB1:iB2+1)
    end if


  end subroutine SetMassAssInfo_FromArray

  !****************************************************************************
  !****s* MassAssInfoDefinition/SetUpperBin
  ! NAME
  ! subroutine SetUpperBin(MAI,maxmass,gamma,nBin)
  ! PURPOSE
  ! Change the information stored in MAI such, that the maximal mass value
  ! is maxmass. The mass bins above the mass cut are deleted. For the bin
  ! containing the cut the upper boundary is set to the max value. Also
  ! The weight W (and Y) is adjusted. If the Q value stored for the bin
  ! was negative, also this value is adjusted (negative Q indicates, that
  ! the maximum was at the upper boundary).
  ! INPUTS
  ! * type(tMassAssInfo) :: MAI -- the Info to be set
  ! * real :: maxmass -- the upper limit of the mass range
  ! * real :: gamma -- the width at maxmass
  ! OUTPUT
  ! * integer :: nBin -- number of bins remaining
  ! * MAI is changed
  !****************************************************************************

  subroutine SetUpperBin(MAI,maxmass,gamma,nBin)
!    use CALLSTACK, only : TRACEBACK

    type(tMassAssInfo), intent(inout) :: MAI
    real, intent(in) :: maxmass
    real, intent(in) :: gamma
    integer, intent(out) :: nBin

    nBin = size(MAI%W)

    if (MAI%IsStable) return

    if (maxmass .lt. MAI%M(1)) then
       ! do something
    end if
    if (maxmass .gt. MAI%M(nBin+1)) then
       ! do something
    end if

!!$    write(*,*) '-----------------------'
!!$    write(*,*) '-----------------------'
!!$    write(*,'(A10,20f13.4)') 'BinM',MAI%M
!!$    write(*,'(A10,7(" "),20f13.4)') 'BinW',MAI%W(1:nBin)/sum(MAI%W(1:nBin))
!!$    write(*,'(A10,7(" "),20f13.4)') 'BinQ',MAI%Q(1:nBin)
!!$    write(*,*) '-----------------------'
    do
       if (nBin.eq.1) exit
       if (MAI%M(nBin).lt.maxmass) exit

       MAI%W(nBin) = 0.0
       nBin=nBin-1

    end do

!    if ((MAI%M(nBin).gt.maxmass).or.(MAI%M(nBin+1).lt.maxmass)) then
!       call TRACEBACK('something strange')
!    end if

    ! adjust W

    MAI%M(nBin+1) = maxmass
    MAI%Y(nBin+1) = 2*atan2(2*(maxmass-MAI%mass0),MAI%gamma0)
    MAI%W(nBin) = MAI%Y(nBin+1)-MAI%Y(nBin)

    ! adjust Q (if possible)

    if (MAI%Q(nBin).lt.0.0) then
       MAI%Q(nBin) = MassAssInfoQ(MAI%mass0,MAI%gamma0,maxmass,gamma)
    end if

    MAI%Q = abs(MAI%Q)

!!$    write(*,'(A10,20f13.4)') 'BinM',MAI%M
!!$    write(*,'(A10,7(" "),20f13.4)') 'BinW',MAI%W(1:nBin)/sum(MAI%W(1:nBin))
!!$    write(*,'(A10,7(" "),20f13.4)') 'BinQ',MAI%Q(1:nBin)
!!$    write(*,*) '-----------------------'
!!$    write(*,*)

  end subroutine SetUpperBin


  !****************************************************************************
  !****f* MassAssInfoDefinition/MassAssInfoQ
  ! NAME
  ! real function MassAssInfoQ(Mass0,Gamma0,Mass,Gamma)
  !
  ! PURPOSE
  ! calculate the weight Q=A/B, while A is the relativistic and B is the
  ! nonrelativistic Breit-Wigner.
  ! INPUTS
  ! * real :: Mass0 -- the pole mass
  ! * real :: Gamma0 -- the width at the pole mass
  ! * real :: Mass -- the actual mass
  ! * real :: Gamma -- the actual width
  !****************************************************************************
  real function MassAssInfoQ(Mass0,Gamma0,Mass,Gamma)
    real, intent(in) :: Mass0,Gamma0,Mass,Gamma

    real :: A,B
    real :: mass2,mass02,gamma024

    mass2  = mass**2
    mass02 = mass0**2
    gamma024 = gamma0**2/4

    A = mass2*gamma0*gamma / ((mass2-mass02)**2+mass2*gamma**2)
    B = gamma024 / ((mass-mass0)**2+gamma024)
    MassAssInfoQ = A/B

  end function MassAssInfoQ

  !****************************************************************************
  !****is* MassAssInfoDefinition/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "MassAssInfo".
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    integer :: ios

    !**************************************************************************
    !****n* MassAssInfoDefinition/MassAssInfo
    ! NAME
    ! NAMELIST MassAssInfo
    ! includes the parameters:
    ! * UseMassAssInfo
    !**************************************************************************
    NAMELIST /MassAssInfo/ UseMassAssInfo

    call Write_ReadingInput('MassAssInfo',0)
    rewind(5)
    read(5,nml=MassAssInfo,IOSTAT=IOS)
    call Write_ReadingInput('MassAssInfo',0,IOS)

    write(*,*) 'UseMassAssInfo = ',UseMassAssInfo

    call Write_ReadingInput('MassAssInfo',1)

    initFlag=.false.
  end subroutine readInput

  !****************************************************************************
  !****f* MassAssInfoDefinition/Get_UseMassAssInfo
  ! NAME
  ! logical function Get_UseMassAssInfo()
  ! PURPOSE
  ! return the value of UseMassAssInfo
  !****************************************************************************
  logical function Get_UseMassAssInfo()
    if (initFlag) call readInput
    Get_UseMassAssInfo = UseMassAssInfo
  end function Get_UseMassAssInfo

  !****************************************************************************
  !****f* MassAssInfoDefinition/Set_UseMassAssInfo
  ! NAME
  ! subroutine Set_UseMassAssInfo(value)
  ! PURPOSE
  ! Override the value of UseMassAssInfo
  !****************************************************************************
  subroutine Set_UseMassAssInfo(value)
    logical, intent(in) :: value
    if (initFlag) call readInput
    if (UseMassAssInfo .neqv. value) then
       write(*,*) 'WARNING: Overriding value read form jobcard. New value:'
       write(*,*) '  UseMassAssInfo = ',value
    end if
    UseMassAssInfo = value
  end subroutine Set_UseMassAssInfo

end module MassAssInfoDefinition
