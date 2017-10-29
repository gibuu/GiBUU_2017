!******************************************************************************
!****m* /preEventDefinition
! NAME
! Module  preEventDefinition
! PURPOSE
! Here the type(preEvent) is defined.
! Routines handling this type are defined elsewhere.
!******************************************************************************
module preEventDefinition

  implicit none

  !****************************************************************************
  !****t*  preEventDefinition/preEvent
  ! SOURCE
  !
  Type preEvent
     real    :: mass = 0.
     Integer :: ID = 0
     Integer :: charge = 0
     logical :: antiparticle = .false.
  End Type preEvent
  !
  ! PURPOSE
  ! This is the minimal information about a particle.
  ! An array of this type is used by XSectionMaster etc. in order to store
  ! the outgoing particles.
  !
  ! NOTES
  ! We also use this type for listings of all possible out--channels.
  ! It would be worthwile to rethink the naming.
  !****************************************************************************

end module preEventDefinition
