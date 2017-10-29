!******************************************************************************
!****m* /particlePointerListDefinition
! NAME
! module particlePointerListDefinition
!
! PURPOSE
! This module defines all necesary types for storing pointers to particles.
! This includes lists of particles and lists to lists and ... ;)
!
! Routines to handle these type are defined elsewhere.
!
! INPUTS
! (none)
!******************************************************************************
module particlePointerListDefinition

  use particleDefinition, only: particle

  !****************************************************************************
  !****t* particlePointerListDefinition/tParticleListNode
  ! PURPOSE
  ! Define the node of a (single linked) particle list:
  ! * pointer to type(particle)
  ! * pointer to (next) tParticleListNode
  !
  ! The data is not stored directly but connected via a "pointer" to some
  ! entry elsewhere, typically in some "particle vector".
  !
  ! SOURCE
  !
  type tParticleListNode
     sequence
     type(particle)         , POINTER :: V    => null()
     type(tParticleListNode), POINTER :: next => null()
  end type tParticleListNode
  !****************************************************************************

  !****************************************************************************
  !****t* particlePointerListDefinition/tParticleList
  ! PURPOSE
  ! Store some elements of a (single linked) list of particles.
  !
  ! The list is build up via a chain of nodes according "tParticleListNode"
  !
  ! NOTES
  ! Keeping the number of entries in a seperate value is not really
  ! necessary for managing the list, but is included in aspect of the
  ! usage of this type: e.g. one maybe wants to remove the fifth entry out
  ! of nine, while the value "5" was choosen by a random generator.
  !
  ! SOURCE
  !
  type tParticleList
     sequence
     type(tParticleListNode), POINTER :: first => null()
     type(tParticleListNode), POINTER :: last  => null()
     integer :: nEntries = 0
  end type tParticleList
  !****************************************************************************


end module particlePointerListDefinition
