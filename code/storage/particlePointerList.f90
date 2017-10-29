!******************************************************************************
!****m* /particlePointerList
! NAME
! module particlePointerList
!
! PURPOSE
! This module defines all necesary types for storing pointers to particles.
! This includes lists of particles and lists to lists and ... ;)
!
! INPUTS
! (none)
!******************************************************************************
module particlePointerList

  use particleDefinition
  use particlePointerListDefinition

  implicit none

  private

  public :: ParticleList_Print
  public :: ParticleList_getParticle
  public :: ParticleList_INIT, ParticleList_CLEAR
  public :: ParticleList_APPEND
  public :: ParticleList_PREPEND
!   public :: ParticleList_REMOVE

contains

  !****************************************************************************
  !****s* particlePointerList/ParticleList_Print
  ! NAME
  ! subroutine ParticleList_Print(L)
  !
  ! PURPOSE
  ! Loop over the List and print every particle
  !
  ! INPUTS
  ! * type(tParticleList) :: L -- The List
  !
  ! OUTPUT
  ! written on stdout
  !****************************************************************************
  subroutine ParticleList_Print(L)
    use output, only: WriteParticle

    type(tParticleList) :: L
    type(tParticleListNode), POINTER :: pNode

    write(*,*) "ParticleList_Print: ", L%nEntries

    pNode => L%first
    do
       if (.not. associated(pNode)) exit

       call WriteParticle(6,0,0,pNode%V)
!       write(*,*) pNode%V%A, pNode%V%B

       pNode=>pNode%next
    end do

  end subroutine ParticleList_Print


  !****************************************************************************
  !****s* particlePointerList/ParticleList_getParticle
  ! NAME
  ! logical function ParticleList_getParticle(L, ID, n, P, antiparticle)
  !
  ! PURPOSE
  ! Loop over the List "L" and find the "n"-th particle in the list with
  ! %ID="ID", "%antiparticle=antiparticle" and
  ! %charge="charge".  This particle "P" is then returned.
  !
  ! INPUTS
  ! * type(tParticleList) :: L -- The List of particles
  ! * integer :: ID     -- ID of particle which shall be returned
  ! * integer :: charge -- charge of particle which shall be returned
  ! * integer :: n      -- We return the n-th particle in the list with %ID=ID
  ! * logical :: antiparticle -- .false. if we search for a particle, .true. for an antiparticle
  ! * logical, OPTIONAL :: weightNonZero -- if .true. only count particles with perweight > 0
  !
  ! OUTPUT
  ! * type(particle) :: P -- n-th particle with wished ID
  ! * logical        :: success -- True if it was possible to find n-Particles
  !   with the wished ID, False otherwise
  !****************************************************************************
  function ParticleList_getParticle(L, ID, charge, n, P, antiparticle, weightNonZero) result (success)

    type(tParticleList),intent(in) :: L
    integer,intent(in) :: ID, charge, n
    type(particle), intent(out) ::  P
    logical, intent(in) :: antiparticle
    logical, intent(in), optional :: weightNonZero
    logical :: success

    type(tParticleListNode), POINTER :: pNode
    integer :: foundIDs
    logical :: checkWeight

    checkWeight =.false.
    if (present(weightNonZero)) checkWeight = weightNonZero

    ! Default return values
    success = .false.
    call SetToDefault(p)

    ! Search particle:
    pNode => L%first
    foundIDs = 0
    do
       if (.not. associated(pNode)) exit

       if (pNode%V%ID == ID) then
          if (pNode%V%charge == charge) then
             if (pNode%V%antiparticle.eqv.antiparticle) then
                if (checkWeight) then
                   if (pNode%V%perWeight > 1e-20) then
                      foundIDs = foundIDs + 1
                   end if
                else
                   foundIDs = foundIDs + 1
                end if
                if (foundIDs == n) then
                   ! n particles with the wanted ID are found, and the n-th particle is returned.
                   p = pNode%V
                   success = .true.
                   return
                end if
             end if
          end if
       end if
       pNode => pNode%next
    end do

  end function ParticleList_getParticle




  !****************************************************************************
  !****s* particlePointerList/ParticleList_INIT
  ! NAME
  ! subroutine ParticleList_INIT(L)
  !
  ! PURPOSE
  ! Initialize the List
  ! (call only at start; to reset the list please call ParticleList_CLEAR)
  !
  ! INPUTS
  ! * type(tParticleList) :: L -- The List
  !
  ! OUTPUT
  ! (none)
  !****************************************************************************
  subroutine ParticleList_INIT (L)

    type(tParticleList) :: L

    NULLIFY(L%first,L%last)
    L%nEntries=0

  end subroutine ParticleList_INIT



  !****************************************************************************
  !****s* particlePointerList/ParticleList_CLEAR
  ! NAME
  ! subroutine ParticleList_CLEAR(L,all)
  !
  ! PURPOSE
  ! Reset the List: Delete all Nodes and re-init the pointers
  !
  !
  ! INPUTS
  ! * type(tParticleList) :: L -- The List
  ! * logical, OPTIONAL   :: all -- if present and true, also the particle is deallocated
  !
  ! OUTPUT
  ! (none)
  !****************************************************************************
  subroutine ParticleList_CLEAR(L,all)

    type(tParticleList) :: L
    logical,optional :: all

    type(tParticleListNode), POINTER :: pNode,pNodeP

    if (L%nEntries>0) then
       pNodeP => L%first
       do
          if (.not. associated(pNodeP)) exit
          pNode => pNodeP%next
          if (present(all)) then
             if (all) deallocate(pNodeP%V)
          end if
          deallocate(pNodeP)
          pNodeP=>pNode
       end do
    end if
    call ParticleList_INIT(L)

  end subroutine ParticleList_CLEAR



  !****************************************************************************
  !****s* particlePointerList/ParticleList_APPEND
  ! NAME
  ! subroutine ParticleList_APPEND(L,V)
  !
  ! PURPOSE
  ! Append the particle (which V points at) at the end of the list.
  !
  !
  ! INPUTS
  ! * type(tParticleList)     :: L -- The List
  ! * type(particle), POINTER :: V -- The particle to add
  !
  ! OUTPUT
  ! (none)
  !****************************************************************************
  subroutine ParticleList_APPEND(L, V)

    type(tParticleList)              :: L
    type(particle)         , POINTER :: V

    type(tParticleListNode), POINTER :: pNode

    allocate(pNode)
    NULLIFY(pNode%next)
    pNode%V => V

    if (.not. associated(L%first)) then
       L%first => pNode
    else
       L%last%next => pNode
    end if
    L%last => pNode

    L%nEntries = L%nEntries+1

  end subroutine ParticleList_APPEND


  !****************************************************************************
  !****s* particlePointerList/ParticleList_PREPEND
  ! NAME
  ! subroutine ParticleList_PREPEND(L,V)
  !
  ! PURPOSE
  ! Prepend the particle (which V points at) at the beginning of the list.
  !
  !
  ! INPUTS
  ! * type(tParticleList)     :: L -- The List
  ! * type(particle), POINTER :: V -- The particle to add
  !
  ! OUTPUT
  ! (none)
  !****************************************************************************
  subroutine ParticleList_PREPEND(L, V)

    type(tParticleList)              :: L
    type(particle)         , POINTER :: V

    type(tParticleListNode), POINTER :: pNode

    allocate(pNode)
    NULLIFY(pNode%next)
    pNode%V => V

    if (.not. associated(L%last)) then
       L%last => pNode
    else
       pNode%next => L%first
    end if
    L%first => pNode

    L%nEntries = L%nEntries+1

  end subroutine ParticleList_PREPEND


  !****************************************************************************
  !****f* particlePointerList/ParticleList_REMOVE
  ! NAME
  ! logical function ParticleList_REMOVE(L,iEntry,V)
  !
  ! PURPOSE
  ! Remove the iEntry-th particle from the List.
  ! The pointer to this particle is returned, the memory of the node is
  ! freed.
  !
  ! NOTES
  ! This routine is fastest if one removes the first particle. For all
  ! higher indices, the routine has to loop over the entries of the list
  ! until it has found the right number.
  !
  ! INPUTS
  ! * type(tParticleList)     :: L -- The List
  ! * integer                 :: iEntry -- The number of the particle to
  !   remove
  ! * type(particle), POINTER :: V -- The particle to add
  !
  ! OUTPUT
  ! * type(particle), POINTER :: V -- The particle removed
  ! * return value -- TRUE at success (index was in allowed range)
  !****************************************************************************
!   logical function ParticleList_REMOVE(L,iEntry, V)
!
!     type(tParticleList)              :: L
!     type(particle)         , POINTER :: V
!     integer                          :: iEntry
!
!     type(tParticleListNode), POINTER :: pNode,pNodeP
!     integer :: i
!
!     NULLIFY(V)
!     ParticleList_REMOVE = .FALSE.
!
!     if (iEntry < 1) return
!     if (iEntry > L%nEntries) return
!
!     pNode => L%first
!     i = 2
!
!     if (L%nEntries == 1) then
!        V => pNode%V
!        DEALLOCATE(pNode)
!        NULLIFY(L%first,L%last)
!        L%nEntries = 0
!        ParticleList_REMOVE = .TRUE.
!        return
!     end if
!
!     if (iEntry == 1) then
!        V => pNode%V
!        L%first => pNode%next
!        DEALLOCATE(pNode)
!        L%nEntries = L%nEntries-1
!        ParticleList_REMOVE = .TRUE.
!        return
!     end if
!
!     do
!        if (i == iEntry) exit ! now pNode point to the precessor
!        i = i+1
!        pNode => pNode%next
!     end do
!
!     pNodeP => pNode
!     pNode => pNode%next
!
!     V => pNode%V
!     pNodeP%next => pNode%next
!
!     DEALLOCATE(pNode)
!     if (iEntry == L%nEntries) L%last => pNodeP
!     L%nEntries = L%nEntries-1
!     ParticleList_REMOVE = .TRUE.
!     return
!
!   end function ParticleList_REMOVE

end module particlePointerList
