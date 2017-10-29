!******************************************************************************
!****m* /PreEvList
! NAME
! module PreEvList
! PURPOSE
! Here the handling of "type tPreEvList" are given.
!******************************************************************************
module PreEvList
  use PreEvListDefinition
  use AnaEventDefinition

  implicit none
  private

  public :: CreateSortedPreEvent
  public :: ComparePreEvent
  public :: PreEvList_INIT, PreEvList_CLEAR, PreEvList_INSERT
  public :: PreEvList_Print, PreEvList_PrintEntry
  public :: PreEvList_GET


  !****************************************************************************
  !****f* PreEvList/CreateSortedPreEvent
  ! NAME
  ! logical function CreateSortedPreEvent(E, PreE, chargeZero)
  ! logical function CreateSortedPreEvent(Parts, PreE, chargeZero)
  !
  ! PURPOSE
  ! This routine creates a sorted list of particles on the basis of
  ! type "preEvent" out of the given tEvent.
  !
  ! INPUTS
  ! * logical, OPTIONAL :: chargeZero -- whether to set the charge to 0
  ! * type(tAnaEvent) ::              E     -- The event
  ! or
  ! * type(particle), dimension(:) :: Parts -- The Particles
  !
  ! OUTPUT
  ! * type(preEvent), dimension(:) :: PreE -- the list of particles
  ! * return value flags success or not
  !
  ! NOTES
  ! Only a maximum number of particles is considered, at the moment: 40.
  !
  !****************************************************************************
  Interface CreateSortedPreEvent
      Module Procedure CreateSortedPreEvent_E,CreateSortedPreEvent_P
   end Interface

contains
  !****************************************************************************
  ! cf. interface "CreateSortedPreEvent" :
  !****************************************************************************

  !-------------------------------------------------------------------------
  logical function CreateSortedPreEvent_E(E, preE, chargeZero)
    use sorting, only: indexx
    use particleDefinition
    use particlePointerListDefinition

    type(tAnaEvent), intent(in) :: E
    type(preEvent), allocatable, dimension(:), intent(out) :: preE
    logical, intent(in), optional :: chargeZero

    integer, parameter :: nPartMax0=40
    real, dimension(1:nPartMax0) :: idCodes
    integer, dimension(1:nPartMax0) :: ind
    type(tParticleListNode),dimension(1:nPartMax0) :: pPartx
    type(tParticleListNode), POINTER  :: pNode
    integer :: i,nPart
    logical :: doZero


    CreateSortedPreEvent_E = .false.

    nPart = E%particleList%nEntries
    if (nPart<1 .or. nPart>nPartMax0) return

    doZero = .false.
    if (present(chargeZero)) then
       doZero = chargeZero
    end if

    i = 1
    pNode => E%particleList%first
    do
       if (.not. associated(pNode)) exit
       pPartx(i)%V => pNode%V

       if (doZero) then
          idCodes(i) = pNode%V%ID*100 + (5)*10 + 2
       else
          idCodes(i) = pNode%V%ID*100 + (pNode%V%charge+5)*10 + 2
       end if
       if (pNode%V%antiparticle) idCodes(i) = idCodes(i) - 1

       i = i+1
       pNode => pNode%next
    end do

    call indexx(idCodes(1:nPart),ind(1:nPart))

    allocate(preE(1:nPart))
    do i=1,nPart
       preE(i)%ID = pPartx(ind(i))%V%ID
       if (doZero) then
          preE(i)%charge = 0
       else
          preE(i)%charge = pPartx(ind(i))%V%charge
       end if
       preE(i)%antiparticle = pPartx(ind(i))%V%antiparticle
       preE(i)%mass = idCodes(ind(i)) ! abuse of mass for storage of code
    end do

    CreateSortedPreEvent_E = .true.

  end function CreateSortedPreEvent_E

  !-------------------------------------------------------------------------

  logical function CreateSortedPreEvent_P(Parts, preE, chargeZero)
    use sorting, only: indexx
    use particleDefinition

    type(particle), dimension(:), intent(in) :: Parts
    type(preEvent), allocatable, dimension(:), intent(out) :: preE
    logical, intent(in), optional :: chargeZero

    integer :: i, nPart
    integer, parameter :: nPartMax0=40
    real, dimension(1:nPartMax0) :: idCodes
    integer, dimension(1:nPartMax0) :: ind
    logical :: doZero

    CreateSortedPreEvent_P = .false.

    doZero = .false.
    if (present(chargeZero)) then
       doZero = chargeZero
    end if


    nPart = 0
    do i=1,size(Parts,dim=1)
       if (Parts(i)%ID <= 0) exit
       nPart = i
       if (i>nPartMax0) exit

       if (doZero) then
          idCodes(i) = Parts(i)%ID*100 + (5)*10 + 2
       else
          idCodes(i) = Parts(i)%ID*100 + (Parts(i)%charge+5)*10 + 2
       end if
       if (Parts(i)%antiparticle) idCodes(i) = idCodes(i) - 1
    end do

    if (nPart==0 .or. nPart>nPartMax0) return

    call indexx(idCodes(1:nPart),ind(1:nPart))

    allocate(preE(1:nPart))
    do i=1,nPart
       preE(i)%ID     = Parts(ind(i))%ID
       if (doZero) then
          preE(i)%charge = 0
       else
          preE(i)%charge = Parts(ind(i))%charge
       end if
       preE(i)%antiparticle = Parts(ind(i))%antiparticle
       preE(i)%mass   = idCodes(ind(i)) ! abuse of mass for storage of code
    end do

    CreateSortedPreEvent_P = .true.

  end function CreateSortedPreEvent_P
  !-------------------------------------------------------------------------

  !****************************************************************************
  !****f* PreEvList/ComparePreEvent
  ! NAME
  ! integer function ComparePreEvent
  ! PURPOSE
  ! return -1,0,1 as ordering parameter of preEvents
  !****************************************************************************
  integer function ComparePreEvent(preE1, preE2)

    type(preEvent), dimension(:), intent(in) :: preE1, preE2

    integer :: i

    if (size(preE1) < size(preE2)) then
       ComparePreEvent = -1
       return
    end if
    if (size(preE1) > size(preE2)) then
       ComparePreEvent = 1
       return
    end if

    do i=1,size(preE1)
       if (preE1(i)%mass < preE2(i)%mass) then
          ComparePreEvent = -1
          return
       end if
       if (preE1(i)%mass > preE2(i)%mass) then
          ComparePreEvent = 1
          return
       end if
    end do

    ComparePreEvent = 0

  end function ComparePreEvent


  !****************************************************************************
  !****s* PreEvList/PreEvList_INIT
  ! NAME
  ! subroutine PreEvList_INIT(L)
  ! PURPOSE
  ! Initialize the List
  ! (call only at start; to reset the list please call PreEvList_CLEAR)
  !
  ! INPUTS
  ! * type(tPreEvList) :: L -- The List
  !
  ! OUTPUT
  ! (none)
  !****************************************************************************
  subroutine PreEvList_INIT(L)
    type(tPreEvList) :: L

    NULLIFY(L%first,L%last)
    L%nEntries=0
  end subroutine PreEvList_INIT

  !****************************************************************************
  !****s* PreEvList/PreEvList_CLEAR
  ! NAME
  ! subroutine PreEvList_CLEAR(L)
  ! PURPOSE
  ! Reset the List: Delete all Nodes and Entries and re-init the pointers
  !
  ! INPUTS
  ! * type(tPreEvList) :: L -- The List
  !
  ! OUTPUT
  ! (none)
  !****************************************************************************
  subroutine PreEvList_CLEAR(L)
    type(tPreEvList) :: L

    type(tPreEvListNode), POINTER :: pNode,pNodeP

    if (L%nEntries>0) then
       pNodeP => L%first
       do
          if (.not. associated(pNodeP)) exit
          pNode => pNodeP%next
          deallocate(pNodeP%V%preE)
          deallocate(pNodeP%V)
          deallocate(pNodeP)
          pNodeP=>pNode
       end do
    end if
    call PreEvList_INIT(L)

  end subroutine PreEvList_CLEAR

  !****************************************************************************
  !****s* PreEvList/PreEvList_INSERT
  ! NAME
  ! subroutine PreEvList_INSERT(L, V, ipos)
  ! PURPOSE
  ! Insert the entry "V" into the list "L", if it is not already in.
  ! In this case, just sum up the value of "weight" of the corresponding
  ! node
  ! INPUTS
  ! * type(tPreEvList)      :: L -- The List
  ! * type(tPreEvListEntry) :: V -- The Entry
  ! OUTPUT
  ! * integer, OPTIONAL :: ipos -- the position in the list
  ! * The list is modified
  !
  ! NOTES
  ! one could think about building up a sorted list instead of just
  ! appending new entries at the end.
  !****************************************************************************
  subroutine PreEvList_INSERT(L, V, ipos)
    type(tPreEvList) :: L
    type(tPreEvListEntry) :: V
    integer, OPTIONAL, intent(OUT) :: ipos

    type(tPreEvListNode), POINTER :: pNode
    integer :: i,j,n
    type(tPreEvListEntry), POINTER :: Vnew

    pNode => L%first
    j=1
    n = size(V%preE)
    nodeLoop: do
       if (.not. associated(pNode)) exit

       if (size(pNode%V%preE) .ne. n) then
          pNode=>pNode%next
          j=j+1
          cycle nodeLoop
       end if
       do i=1,n
          if (pNode%V%preE(i)%mass .ne. V%preE(i)%mass) then ! attention: abuse of mass
             pNode=>pNode%next
             j=j+1
             cycle nodeLoop
          end if
       end do

       ! we have found it:

       pNode%V%weight = pNode%V%weight + V%weight
       if (present(ipos)) ipos=j

       return
    end do nodeLoop

    ! now we have to append it:

    allocate(Vnew)
    allocate(Vnew%preE(n))
    Vnew%weight = V%weight
    Vnew%PreE   = V%PreE
    allocate(pNode)
    NULLIFY(pNode%next)
    pNode%V => Vnew

    if (.not. associated(L%first)) then
       L%first => pNode
    else
       L%last%next => pNode
    end if
    L%last => pNode

    L%nEntries =  L%nEntries+1
    if (present(ipos)) ipos=L%nEntries

  end subroutine PreEvList_INSERT

  !****************************************************************************
  !****s* PreEvList/PreEvList_GET
  ! NAME
  ! logical function PreEvList_GET(L, V, ipos)
  ! PURPOSE
  ! Insert the entry "V" into the list "L", if it is not already in.
  ! In this case, just sum up the value of "weight" of the corresponding
  ! node
  ! INPUTS
  ! * type(tPreEvList)      :: L -- The List
  ! * type(tPreEvListEntry) :: V -- The Entry
  ! * integer :: ipos -- the position in the list
  ! OUTPUT
  ! * type(tPreEvListEntry) :: V -- The Entry
  ! * function value : .true. at success
  !
  !****************************************************************************
  logical function PreEvList_GET(L, V, ipos)
    type(tPreEvList), intent(IN) :: L
    type(tPreEvListEntry), intent(INOUT) :: V
    integer, intent(IN) :: ipos

    type(tPreEvListNode), POINTER :: pNode
    integer :: j,n

    PreEvList_GET = .false.

    pNode => L%first
    j=1
    nodeLoop: do
       if (.not. associated(pNode)) exit

       if (j.eq.ipos) exit

       pNode=>pNode%next
       j=j+1
       cycle nodeLoop

    end do nodeLoop

    if (j.eq.ipos) then
       n = size(pNode%V%preE)

       if (allocated(V%preE)) then
          if (size(V%preE).ne.n) then
             deallocate(V%preE)
             allocate(V%preE(n))
          end if
       else
          allocate(V%preE(n))
       end if
       V%preE = pNode%V%preE
       V%weight = pNode%V%weight
       PreEvList_GET = .true.
    end if


  end function PreEvList_GET


  !****************************************************************************
  !****s* PreEvList/PreEvList_Print
  ! NAME
  ! subroutine PreEvList_Print (iFile, L, fak, n, iBreak, sBreak, withLN, doSort)
  ! PURPOSE
  ! Write the list "L" to file "iFile". Multiply the written weights
  ! by the factor "fak". As a side effect, the list is being sorted (before printing).
  ! INPUTS
  ! * integer          :: iFile -- The file number
  ! * type(tPreEvList) :: L -- The List
  ! * real             :: fak -- the factor to multiply the weights with
  ! * integer,OPTIONAL :: n -- number of columns to print
  ! * integer,OPTIONAL :: iBreak -- some string is inserted after entry nnn
  ! * character*(*),OPTIONAL :: sBreak -- string to insert
  ! * logical, OPTIONAL :: withLN -- print with line numbers
  ! * logical, OPTIONAL :: doSort -- switch sorting on/off
  ! OUTPUT
  ! witten to file "iFile"
  ! NOTES
  ! The format is not very clean.
  !****************************************************************************
  subroutine PreEvList_Print (iFile, L, fak, n, iBreak, sBreak, withLN, doSort)
    use ParticleProperties, only: isStrange

    integer,          intent(IN) :: iFile
    type(tPreEvList), intent(inout) :: L
    real,             intent(IN) :: fak
    integer,OPTIONAL, intent(IN) :: n, iBreak
    character*(*),OPTIONAL, intent(IN) :: sBreak
    logical,OPTIONAL, intent(IN) :: withLN
    logical,OPTIONAL, intent(IN) :: doSort

    type(tPreEvListNode), POINTER :: pNode
    integer :: j,ii
    real,dimension(0:1,0:3) :: Sum
    logical :: flagS,doLN,doS
    real :: W

    doLN = .false.
    if (present(withLN)) doLN = withLN

    doS = .true.
    if (present(doSort)) doS = doSort

    write(iFile,*) "PreEvList_Print: nEntries     = ", L%nEntries

    if (doS) call PreEvList_Sort(L)  ! sort before printing

    Sum = 0
    pNode => L%first
    ii = 1
    do
       if (.not. associated(pNode)) exit

       if (doLN) then
          call PreEvList_PrintEntry(iFile,pNode%V,fak,n,iBreak,sBreak,ii)
       else
          call PreEvList_PrintEntry(iFile,pNode%V,fak,n,iBreak,sBreak)
       end if

       flagS = .false.
       do j=1,size(pNode%V%preE)
          flagS = flagS .or. isStrange(pNode%V%preE(j)%ID)
       end do
       W = pNode%V%weight * fak

       Sum(0:1,0) = Sum(0:1,0) + (/1.0,W/)
       if (flagS) Sum(0:1,1) = Sum(0:1,1) + (/1.0,W/)

       pNode=>pNode%next
       ii = ii+1
    end do

    write(iFile,'(A,f6.0,f9.3)') "total  : ", Sum(0:1,0)
    write(iFile,'(A,f6.0,f9.3)') "strange: ", Sum(0:1,1)

  end subroutine PreEvList_Print


  !****************************************************************************
  !****s* PreEvList/PreEvList_PrintEntry
  ! NAME
  ! subroutine PreEvList_PrintEntry(iFile,E,fak,n,iBreak,sBreak,iLN,doWeight,noAdvance)
  ! PURPOSE
  ! Print a tPreEvListEntry to output.
  ! INPUTS
  ! * integer          :: iFile -- The file number
  ! * type(tPreEvListEntry) :: E -- The List-Entry
  ! * real             :: fak -- the factor to multiply the weights with
  ! * integer,OPTIONAL :: n -- number of columns to print
  ! * integer,OPTIONAL :: iBreak -- some string is inserted after entry nnn
  ! * character*(*),OPTIONAL :: sBreak -- string to insert
  ! * integer, OPTIONAL :: iLN -- line number to print
  ! * logical, OPTIONAL :: doWeight -- print weight or not
  ! * logical, OPTIONAL :: noAdvance -- print line feed or not
  ! OUTPUT
  ! witten to file "iFile"
  !
  ! NOTES
  ! The format is not very clean.
  !****************************************************************************
  subroutine PreEvList_PrintEntry(iFile,E,fak,n,iBreak,sBreak,iLN,doWeight,noAdvance)
    use ParticleProperties, only: PartName

    integer,          intent(IN) :: iFile
    type(tPreEvListEntry),  intent(IN)  :: E
    real,             intent(IN) :: fak
    integer,OPTIONAL, intent(IN) :: n
    integer,OPTIONAL, intent(IN) :: iBreak
    character*(*),OPTIONAL, intent(IN) :: sBreak
    integer, OPTIONAL, intent(IN) :: iLN
    logical, OPTIONAL, intent(IN) :: doWeight
    logical, OPTIONAL, intent(IN) :: noAdvance

    character*(15), dimension(20) :: AA
    integer :: nAA
    integer :: j, iiB
    character*(100) :: BUF,ssB
    character*(1000):: BUF2
    real :: W
    logical :: doLN, doW, noAdv

    nAA = 8 ! default value
    if (present(n)) nAA = min(n,20)

    doLN = .false.
    if (present(iLN)) doLN = .true.

    doW = .true.
    if (present(doWeight)) doW = doWeight

    noAdv = .false.
    if (present(noAdvance)) noAdv = noAdvance

    iiB = 0 ! default value
    if (present(iBreak)) iiB = iBreak

    if (iiB.gt.0) then
       ssB = ' : '
       if (present(sBreak)) ssB = trim(sBreak)
       write(BUF,1001) iiB,trim(ssB)
    else
       write(BUF,1000)
    end if


1000 format ("(20A16)")
1001 format ("(",i2,"A16,'",A,"',20A16)")

!    write(*,*) 'PreEvList_PrintEntry: BUF=#',trim(BUF),'#'

    nAA = min(nAA,size(E%preE))
    AA = ""
    do j=1,nAA
       AA(j) = PartName(E%preE(j)%ID,E%preE(j)%charge,E%preE(j)%antiparticle)
    end do
    W = E%weight * fak

    if (doLN) write(iFile, '(i9)', advance='no') iLN

    if (doW) write(iFile, '(f12.5)', advance='no') W

    write(BUF2,BUF) AA

    if (noAdv) then
       write(iFile,'(A)',advance='no') trim(BUF2)
    else
       write(iFile,'(A)') trim(BUF2)
    end if

  end subroutine PreEvList_PrintEntry


  !****************************************************************************
  !****s* PreEvList/PreEvList_Sort
  ! NAME
  ! subroutine PreEvList_Sort (L)
  ! PURPOSE
  ! Sort the list 'L' by weights (using the 'mergesort' algorithm).
  ! INPUTS
  ! * type(tPreEvList) :: L -- the unsorted list
  ! OUTPUT
  ! * type(tPreEvList) :: L -- the sorted list
  ! NOTE
  ! This routine represents a rather direct translation of the C code from
  ! http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html
  !****************************************************************************
  subroutine PreEvList_Sort (L)
    type(tPreEvList), intent(inout) :: L

    type(tPreEvListNode), pointer :: head, tail, p, q, e
    integer :: insize, nmerges, psize, qsize, i

    head => L%first
    if (.not. associated(head)) return

    insize = 1

    do while (.true.)
      p => head
      head => NULL()
      tail => NULL()

      nmerges = 0  ! count number of merges we do in this pass

      do while (associated(p))
         nmerges=nmerges+1  ! there is a merge to be done
         ! step `insize' places along from p
         q => p
         psize = 0
         do i = 1, insize
            psize=psize+1
            q => q%next
            if (.not. associated(q)) exit
         end do

         ! if q hasn't fallen off end, we have two lists to merge
         qsize = insize

         ! now we have two lists: merge them
         do while (psize > 0 .or. (qsize > 0 .and. associated(q)))

            ! decide whether next element of merge comes from p or q
            if (psize == 0) then
               ! p is empty, e must come from q
               e => q; q => q%next; qsize=qsize-1
            else if (qsize == 0 .or. .not. associated(q)) then
               ! q is empty, e must come from p
               e => p; p => p%next; psize=psize-1
            else if (p%V%weight>q%V%weight) then
               ! first element of p is larger, e must come from p
               e => p; p => p%next; psize=psize-1
            else
               ! first element of q is larger, e must come from q
               e => q; q => q%next; qsize=qsize-1
            end if

            ! add the next element to the merged list
            if (associated(tail)) then
               tail%next => e
            else
               head => e
            end if
            tail => e
         end do

         ! now p has stepped `insize' places along, and q has too
         p => q
      end do

      tail%next => NULL()

      ! if we have done only one merge, we're finished
      if (nmerges <= 1) exit

      ! otherwise repeat, merging lists twice the size
      insize = insize * 2
    end do

    L%first => head
    L%last  => tail

  end subroutine PreEvList_Sort


end module PreEvList
