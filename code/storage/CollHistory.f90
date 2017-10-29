!******************************************************************************
!****m* /CollHistory
! NAME
! module CollHistory
! PURPOSE
! Here we can store for debugging purposes the full collision history of
! every (perturbative) particle from the particle vector
!
! Every collision type is stored in a list of type tPreEvList
! (In fact 1-Body, 2-Body and 3-Body collisions are stored in three
! independent lists).
! We identify the collision type by its position in the list. Additionally
! we encode the "list" by modifying the iEntry value:
! * 1-Body: iEntry = iEntry     (range: 1 ... 10000)
! * 2-Body: iEntry = -iEntry
! * 3-Body: iEntry = iEntry + 10000
! Therefore for every collision type the (modified) iEntry represents an
! unique code number.
!
! For each particle in the perturbative particle vector we store now in an
! additional array the code numbers of the corresponding collition types.
! This array has the same dimensions as the perturbative particle vector +
! one additional dimension -- the subsequent collision type codes.
!
! In a collision we have therefore to delete the array entry of the incoming
! perturbative particle and to create for every final state particle an
! entry in the array which has one more entry in the third dimension than
! the array of the incoming particle. (We use the information
! from insertion/setIntoVector, where the particles were inserted.)
!
! Doing a Garbage Collection on the perturbative particle vector also
! reshuffles the entries. This has to be respected.
!
! On could think of enhancing the actual stored information by maybe also
! storing...:
! * the HiEnergy type of the collision (this would be an additional
!   field in the PreEvListEntry)
! * the time of the collision (this would be an additional field in
!   CollHistEntry)
!
!
!******************************************************************************

module CollHistory
  use preEventDefinition
  use PreEvListDefinition

  implicit none

  private

  !****************************************************************************
  !****it* CollHistory/CollHistEntry
  ! SOURCE
  !
  type CollHistEntry
     integer, allocatable :: Entries(:)
  end type CollHistEntry
  !
  ! PURPOSE
  ! The type to store the collision history. The entries are the codes for
  ! the collsion types. The size of the array corresponds to the generation.
  !****************************************************************************

  !****************************************************************************
  !****ig* CollHistory/CollHistArray
  ! SOURCE
  !
  type(CollHistEntry), allocatable :: CollHistArray(:,:)
  !
  ! PURPOSE
  ! The global array of the size of the perturbative particle vector, which
  ! holds all Collision History information
  !****************************************************************************

  !****************************************************************************
  !****ig* CollHistory/CollList
  ! SOURCE
  !
  type(tPreEvList), save :: CollList1, CollList2, CollList3
  !
  ! PURPOSE
  ! The list to store the 1Body, 2Body and 3Body collision types.
  !****************************************************************************

  logical, save :: initFlag = .true.

  logical, parameter :: verb = .false. ! flag for debugging verbosity


  !****************************************************************************
  !****ig* CollHistory/DoCollHistory
  ! SOURCE
  !
  logical, save :: DoCollHistory = .false.
  !
  ! PURPOSE
  ! Flag to switch on/off the whole Collision History machinery.
  !
  ! You may set this variable via your jobcard, namelist "collHistory"
  !****************************************************************************

  integer, save :: nEns=0, nPart=0

  !****************************************************************************
  !****n* CollHistory/collHistory
  ! NAME
  ! NAMELIST collHistory
  ! PURPOSE
  ! Includes the switches:
  ! * DoCollHistory
  !****************************************************************************

  public :: CollHist_UpdateHist
  public :: CollHist_SetSize
  public :: CollHist_DoGBC
  public :: CollHist_WriteList
  public :: CollHist_WriteHistParticle
  public :: CollHist_ClearArray
  public :: CollHist_ClassifyHist
  public :: CollHist_GetDoCollHistory

contains
  !****************************************************************************
  !****f* CollHistory/CollHist_GetDoCollHistory
  ! NAME
  ! logical function CollHist_GetDoCollHistory
  ! PURPOSE
  ! return value of DoCollHistory
  !****************************************************************************
  logical function CollHist_GetDoCollHistory()
    CollHist_GetDoCollHistory = DoCollHistory
  end function CollHist_GetDoCollHistory


  !****************************************************************************
  !****is* CollHistory/DoInit
  ! NAME
  ! subroutine DoInit
  ! PURPOSE
  ! Do the init of this module.
  ! INPUTS
  ! * via namelist "collHistory"
  ! OUTPUT
  ! * memory is allocated, no explicit output
  !****************************************************************************
  subroutine DoInit
    use output, only: Write_ReadingInput

    integer :: ios

    NAMELIST /collHistory/ DoCollHistory

    ! Read namelist
    rewind(5)
    call Write_ReadingInput('collHistory',0)
    read(5,nml=collHistory,iostat=ios)
    call Write_ReadingInput('collHistory',0,ios)

    write(*,*) 'Doing Collision History:',DoCollHistory

    call Write_ReadingInput('collHistory',1)

    ! Allocate memory, if desired

    if (DoCollHistory) then
       if (nEns.eq.0 .or. nPart.eq.0) then
          write(*,*) 'Call CollHist_SetSize before. stop !'
          stop
       end if
       write(*,*) 'Allocating CollHistArray(',nEns,',',nPart,')'
       allocate(CollHistArray(nEns,nPart))
    end if

    initFlag = .false.

  end subroutine DoInit

  !****************************************************************************
  !****s* CollHistory/CollHist_SetSize
  ! NAME
  ! subroutine CollHist_SetSize(sizes,localflag)
  ! PURPOSE
  ! Prepare the init of this module by storing the size of the perturbtive
  ! particle vector.
  !
  ! INPUTS
  ! * integer, dimension(2) :: sizes -- corresponds to ubound(PertParticles)
  ! OUTPUT
  ! * (none)
  !
  ! NOTES
  ! This routine has to been called before DoInit in order to set the size
  ! of the array.
  !****************************************************************************
  subroutine CollHist_SetSize(sizes)
    use CallStack, only: traceback

    integer, dimension(2), intent(in) :: sizes

    if (initFlag) then
       nEns = sizes(1)
       nPart = sizes(2)
    else
       if (.not.DoCollHistory) return
       if (nEns .ne. sizes(1) .or. nPart .ne. sizes(2)) then
          write(*,*) 'CollHist_SetSize,ooops:',nEns,sizes(1),nPart,sizes(2)
          call traceback()
          stop
       end if
    end if

  end subroutine CollHist_SetSize

  !****************************************************************************
  !****s* CollHistory/CollHist_ClearArray
  ! NAME
  ! subroutine CollHist_ClearArray()
  ! PURPOSE
  ! Clears the used CollHistArray. Necessary after every run/in every Init
  ! call
  !
  ! INPUTS
  ! * (none)
  ! OUTPUT
  ! * (none)
  !
  ! NOTES
  ! The lists of collisions types are unchanged.
  !****************************************************************************
  subroutine CollHist_ClearArray()
    integer i,j

    if (.not.allocated(CollHistArray)) return
    do i=1,nEns
       do j=1,nPart
          if (allocated(CollHistArray(i,j)%Entries)) &
               & deallocate(CollHistArray(i,j)%Entries)
       end do
    end do
  end subroutine CollHist_ClearArray

  !****************************************************************************
  !****s* CollHistory/CollHist_UpdateHist
  ! NAME
  ! subroutine CollHist_UpdateHist(partIn, partOut, posIn, posOut, weight)
  ! PURPOSE
  ! This is the major routine doing all the work for keeping the history
  ! of all particles up to date:
  ! * store the collision type in the corresponding list (if necessary)
  ! * remove history for incoming particle
  ! * enhance history for all (stored) outgoing particles
  !
  ! INPUTS
  ! * type(particle), dimension(:) :: partIn -- incoming particles
  ! * type(particle), dimension(:) :: partOut-- outgoing particles
  ! * integer, dimension(2)   :: posIn  -- (iEns,iPart) of incoming particles
  ! * integer, dimension(:,:) :: posOut -- (iEns,iPart) of outgoing particles
  !   (the array has to be as (2,:))
  ! * real :: weight -- perweight of the incoming perturbative particle
  ! OUTPUT
  ! * internal arrays and lists modified.
  !****************************************************************************
  subroutine CollHist_UpdateHist(partIn, partOut, posIn, posOut, weight)
    use particleDefinition
    use PreEvList, only: CreateSortedPreEvent, PreEvList_INSERT, PreEvList_PrintEntry
    !use CallStack

    type(particle), dimension(:), intent(in) :: partIn
    type(particle), dimension(:), intent(in) :: partOut
    integer, dimension(2),   intent(in) :: posIn
    integer, dimension(:,:), intent(in) :: posOut ! = (2,:)
    real, intent(in) :: weight

    integer :: nPartIn, nPartOut, iEntry, EntrySize, i
    type(preEvent), allocatable, dimension(:) :: preEvIn, preEvOut
    type(tPreEvListEntry) :: Entry
    type(CollHistEntry) :: DummyEntry

    if (initFlag) call DoInit
    if (.not.DoCollHistory) return


    ! Create PreEvent and find it in the list:

    if (.not.CreateSortedPreEvent(partIn,preEvIn)) then
       write(*,*) 'Problems in CollHist_UpdateHist: preEvIn. STOP!'
       stop
    end if
    if (.not.CreateSortedPreEvent(partOut,preEvOut)) then
       write(*,*) 'Problems in CollHist_UpdateHist: preEvOut. STOP!'
       stop
    end if

    nPartIn = size(preEvIn)
    nPartOut = size(preEvOut)

!    write(*,*) 'nPartIn,nPartOut:',nPartIn,nPartOut
!    call TRACEBACK(STRING="just write stack:",USER_EXIT_CODE=-1)

    allocate(Entry%preE(nPartIn+nPartOut))
    Entry%preE(1:nPartIn) = preEvIn(1:nPartIn)
    Entry%preE(nPartIn+1:nPartIn+nPartOut) = preEvOut(1:nPartOut)
    Entry%weight = weight

    select case (nPartIn)
    case (1)
       call PreEvList_INSERT(CollList1,Entry,iEntry)

       if (iEntry.eq.10000) then
          write(*,*) 'Problems in CollHist_UpdateHist: iEntry=10000'
          stop
       end if

    case (2)
       call PreEvList_INSERT(CollList2,Entry,iEntry)
       iEntry = -iEntry

    case (3)
       call PreEvList_INSERT(CollList3,Entry,iEntry)
       iEntry = iEntry+10000

    case default
       write(*,*) 'something fishy'
       stop
    end select

    if (verb) then
       write(*,*) '...iEntry = ',iEntry
       call PreEvList_PrintEntry(6,Entry,1.0,20,iBreak=nPartIn,sBreak="===>")
    end if


    ! Update the histories of the final particles:

    EntrySize = 0
    if ((posIn(1).gt.0).and.(posIn(2).gt.0)) then
       if (allocated(CollHistArray(posIn(1),posIn(2))%Entries)) &
          & EntrySize = size(CollHistArray(posIn(1),posIn(2))%Entries)
    end if
    allocate(DummyEntry%Entries(EntrySize+1))

    if (EntrySize.gt.0) then
       DummyEntry%Entries(1:EntrySize) = CollHistArray(posIn(1),posIn(2))%Entries
       deallocate(CollHistArray(posIn(1),posIn(2))%Entries)

       if (verb) write(*,*) 'CollHist: delete ',posIn(1),posIn(2)
    else
       if (verb) write(*,*) 'CollHist: use    ',posIn(1),posIn(2)
    end if
    EntrySize = EntrySize+1
    DummyEntry%Entries(EntrySize) = iEntry

    do i=1,nPartOut
       if (posOut(1,i).gt.0 .and.posOut(2,i).gt.0) then
          if (verb) write(*,*) 'CollHist: SET    ',posOut(1,i),posOut(2,i),i

          if (allocated(CollHistArray(posOut(1,i),posOut(2,i))%Entries))&
               & deallocate(CollHistArray(posOut(1,i),posOut(2,i))%Entries)
          allocate(CollHistArray(posOut(1,i),posOut(2,i))%Entries(EntrySize))
          CollHistArray(posOut(1,i),posOut(2,i))%Entries = DummyEntry%Entries
       end if

    end do

    deallocate(DummyEntry%Entries)


  end subroutine CollHist_UpdateHist


  !****************************************************************************
  !****s* CollHistory/CollHist_DoGBC
  ! NAME
  ! subroutine CollHist_DoGBC(iEns,iPart1,iPart2)
  ! PURPOSE
  ! This routine replays modifications done by insertion/GarbageCollection:
  ! The information of (iEns,iPart2) has to be moved to (iEns,iPart1)
  ! INPUTS
  ! * integer :: iEns,iPart1,iPart2 -- ensemble and particle position of
  !   enties to move
  ! OUTPUT
  ! * internal arrays modified.
  !****************************************************************************
  subroutine CollHist_DoGBC(iEns,iPart1,iPart2)
    integer, intent(in) :: iEns,iPart1,iPart2
    integer :: nsize=0

    if (.not.DoCollHistory) return
    if (verb) write(*,*) 'CollHist: DoGBC ',iEns,iPart1,'<-',iEns,iPart2

    if (allocated(CollHistArray(iEns,iPart1)%Entries)) &
         & deallocate(CollHistArray(iEns,iPart1)%Entries)

    if (allocated(CollHistArray(iEns,iPart2)%Entries)) then
       nsize = size(CollHistArray(iEns,iPart2)%Entries)
       allocate(CollHistArray(iEns,iPart1)%Entries(nsize))
       CollHistArray(iEns,iPart1)%Entries = CollHistArray(iEns,iPart2)%Entries
       deallocate(CollHistArray(iEns,iPart2)%Entries)
    end if
  end subroutine CollHist_DoGBC


  !****************************************************************************
  !****s* CollHistory/CollHist_WriteList
  ! NAME
  ! subroutine CollHist_WriteList(fak)
  ! PURPOSE
  ! Write all the collision type lists to files.
  ! INPUTS
  ! * real :: fak -- faktor to multiply the weights with
  ! OUTPUT
  ! * files written
  !****************************************************************************
  subroutine CollHist_WriteList(fak)
    use PreEvList

    real, intent(in) :: fak

    if (.not.DoCollHistory) return

    rewind(1001)
    rewind(1002)
    rewind(1003)

    call PreEvList_Print(1001,CollList1,fak,20,iBreak=1,sBreak="===>",withLN=.true.)
    call PreEvList_Print(1002,CollList2,fak,20,iBreak=2,sBreak="===>",withLN=.true.)
    call PreEvList_Print(1003,CollList3,fak,20,iBreak=3,sBreak="===>",withLN=.true.)

  end subroutine CollHist_WriteList

  !****************************************************************************
  !****s* CollHistory/CollHist_WriteHistParticle
  ! NAME
  ! subroutine CollHist_WriteHistParticle(iFile,iEns,iPart)
  ! PURPOSE
  ! Write the history of a specific particle to a file
  ! INPUTS
  ! * integer :: ifile -- file to be used
  ! * integer :: iEns, iPart -- coordinates of the particle in the particle
  !   vector
  ! OUTPUT
  ! * information written to file
  !****************************************************************************
  subroutine CollHist_WriteHistParticle(iFile,iEns,iPart)
    use PreEvList, only: PreEvList_PrintEntry

    integer, intent(in) :: iFile,iEns,iPart

    integer :: i, EntrySize, iEntry, nIn
    logical :: flag
    type(tPreEvListEntry) :: V

    if (.not.DoCollHistory) return

    if (allocated(CollHistArray(iEns,iPart)%Entries)) then
       EntrySize = size(CollHistArray(iEns,iPart)%Entries)
    else
       EntrySize = 0
    end if

!    if (EntrySize .gt. 0) then
!       write(*,*) '>>',EntrySize,CollHistArray(iEns,iPart)%Entries
!    else
!       write(*,*) '>>',EntrySize
!    endif

    do i=1,EntrySize
       iEntry = CollHistArray(iEns,iPart)%Entries(i)

       flag = CollHist_GetV(iEntry, V, nIn)

       if (.not.flag) then
          write(*,*) 'CollHist_WriteHistParticle: Entry not found in list. stop'
          write(*,*) iEntry,nIn
          stop
       end if

       call PreEvList_PrintEntry(iFile,V,1.0,20,iBreak=nIn,sBreak="===>",iLN=i)

    end do

    write(*,*) 'Classify: ',CollHist_ClassifyHist(iEns,iPart)

  end subroutine CollHist_WriteHistParticle

  !****************************************************************************
  !****f* CollHistory/CollHist_GetV
  ! NAME
  ! logical function CollHist_GetV(iEntry, V, nIn)
  ! PURPOSE
  ! with Entry code as input return the corresponding History entry
  ! INPUTS
  ! * integer :: iEntry -- history code
  ! OUTPUT
  ! * type(tPreEvListEntry) :: V -- The PreEvList
  ! * integer :: nIn -- number of incoming particles (1-Body,2-Body,3-Body)
  ! * function value: .true. on success
  !****************************************************************************
  logical function CollHist_GetV(iEntry, V, nIn)
    use PreEvList, only: PreEvList_GET

    integer, intent(IN) :: iEntry
    integer, intent(OUT):: nIn
    type(tPreEvListEntry),intent(inOUT) :: V

    CollHist_GetV = .false.
    nIn = 0

    if (iEntry .lt. 0) then
       nIn = 2
       CollHist_GetV = PreEvList_GET(CollList2, V, -iEntry)
    else if (iEntry .gt. 10000) then
       nIn = 3
       CollHist_GetV = PreEvList_GET(CollList3, V, iEntry-10000)
    else
       nIn = 1
       CollHist_GetV = PreEvList_GET(CollList1, V, iEntry)
    end if

  end function CollHist_GetV

  !****************************************************************************
  !****f* CollHistory/CollHist_ClassifyHist
  ! NAME
  ! integer function CollHist_ClassifyHist(iEns,iPart)
  ! PURPOSE
  ! ...
  ! INPUTS
  ! * integer :: iEns, iPart -- coordinates of the particle in the particle
  !   vector
  ! OUTPUT
  ! * function value
  ! NOTES
  ! This is here a very special version for N+gamma -> pi+- + X
  !****************************************************************************
  integer function CollHist_ClassifyHist(iEns,iPart)

    integer, intent(in) :: iEns,iPart

    integer :: EntrySize, iEntry, nIn, i
    logical :: flag
    type(tPreEvListEntry), dimension(2) :: V
    integer, dimension(2) :: iClass

    CollHist_ClassifyHist = 0

    if (.not.DoCollHistory) return

    if (allocated(CollHistArray(iEns,iPart)%Entries)) then
       EntrySize = size(CollHistArray(iEns,iPart)%Entries)
    else
       EntrySize = 0
    end if

    select case (EntrySize)
    case (0)
       return

    case (1)
       iEntry = CollHistArray(iEns,iPart)%Entries(1)
       flag = CollHist_GetV(iEntry, V(1), nIn)
       iClass(1) = 0
       select case (nIn)
       case (1)
          iClass(1) = CollHist_Classify1Body(V(1))
       case (2)
          if (V(1)%preE(2)%ID.eq.999) iClass(1) = 1
       end select

       CollHist_ClassifyHist = iClass(1)

    case default
       iEntry = CollHistArray(iEns,iPart)%Entries(1)
       flag = CollHist_GetV(iEntry, V(1), nIn)
       if (nIn.ne.2) return
!       if (V(1)%preE(2)%ID.ne.999) return

       iClass(1) = 0
       do i=EntrySize,2,-1

          iEntry = CollHistArray(iEns,iPart)%Entries(i)
          flag = CollHist_GetV(iEntry, V(2), nIn)
          select case (nIn)
          case (1)
             iClass(1) = CollHist_Classify1Body(V(2))
          case (2)
             iClass(1) = CollHist_Classify2Body(V(2))
          end select
          CollHist_ClassifyHist = iClass(1)
          if (i.lt.EntrySize) CollHist_ClassifyHist = CollHist_ClassifyHist*100
          if (i.gt.2) CollHist_ClassifyHist = -CollHist_ClassifyHist

          if (CollHist_ClassifyHist.ne.0) return
       end do


    end select


  end function CollHist_ClassifyHist

  !****************************************************************************

  integer function CollHist_Classify1Body(V)

    type(tPreEvListEntry),intent(in) :: V

    integer :: i
    logical :: isPi

    CollHist_Classify1Body = 0 ! not charcterized

    isPi = .false.
    do i=2,size(V%preE)
       if (V%preE(i)%ID .eq. 101) isPi = .true.
    end do
    if (.not.isPi) return

    select case (V%preE(1)%ID)
    case (103)
       CollHist_Classify1Body = 2 ! rho -> pi + X
    case (105,107)
       CollHist_Classify1Body = 3 ! V(w/o rho) -> pi + X
    case (101,102,104,106,108:121)
       CollHist_Classify1Body = 4 ! m -> pi + X
    case (2)
       CollHist_Classify1Body = 5 ! Delta -> pi + X
    case (1,3:61)
       CollHist_Classify1Body = 6 ! B -> pi + X
    end select

  end function CollHist_Classify1Body

  !****************************************************************************

  integer function CollHist_Classify2Body(V)

    type(tPreEvListEntry),intent(in) :: V

    integer :: i
    logical :: isPi

    CollHist_Classify2Body = 0 ! not charcterized

    isPi = .false.
    do i=3,size(V%preE)
       if (V%preE(i)%ID .eq. 101) isPi = .true.
    end do
    if (.not.isPi) return

    isPi = .false.
    do i=1,2
       if (V%preE(i)%ID .eq. 101) isPi = .true.
    end do
    if (isPi) return

    i=2
    if (V%preE(2)%ID.eq.1) i=1
    select case (V%preE(i)%ID)
    case (103)
       CollHist_Classify2Body = 12 ! N + rho -> pi + X
    case (105,107)
       CollHist_Classify2Body = 13 ! N + V(w/o rho) -> pi + X
    case (102,104,106,108:121)
       CollHist_Classify2Body = 14 ! N + m -> pi + X
    case (2)
       CollHist_Classify2Body = 15 ! N + Delta -> pi + X
    case (1,3:61)
       CollHist_Classify2Body = 16 ! N + B -> pi + X
    case default
       CollHist_Classify2Body = 11 ! A+B -> pi + X, A,B != pi
    end select

  end function CollHist_Classify2Body


end module CollHistory
