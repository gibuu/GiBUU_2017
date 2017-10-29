!******************************************************************************
!****m* /PreEvListDefinition
! NAME
! module PreEvListDefinition
! PURPOSE
! Here the type definitions
! * type tPreEvListEntry
! * type tPreEvListNode
! * type tPreEvList
! are given.
!
! Routines to work with these type are defined elsewhere, as e.g. in
! module PreEvList
!******************************************************************************
module PreEvListDefinition

  use preEventDefinition

  implicit none

  !****************************************************************************
  !****t* PreEvListDefinition/tPreEvListEntry
  ! SOURCE
  !
  type tPreEvListEntry
     type(preEvent), allocatable :: preE(:)
     real :: weight
  end type tPreEvListEntry
  ! PURPOSE
  ! store information about a kind of event
  !****************************************************************************

  !****************************************************************************
  !****t* PreEvListDefinition/tPreEvListNode
  ! SOURCE
  !
  type tPreEvListNode
     type(tPreEvListEntry), POINTER :: V      => null()
     type(tPreEvListNode), POINTER :: next    => null()
  end type tPreEvListNode
  ! PURPOSE
  ! node in order to create a single connected list of "tPreEvListEntry"
  ! elements.
  !****************************************************************************

  !****************************************************************************
  !****t* PreEvListDefinition/tPreEvList
  ! SOURCE
  !
  type tPreEvList
     type(tPreEvListNode), POINTER :: first  => null()
     type(tPreEvListNode), POINTER :: last   => null()
     integer :: nEntries
  end type tPreEvList
  ! PURPOSE
  ! store a single connected list of "tPreEvListEntry" elements.
  !****************************************************************************

end module PreEvListDefinition
