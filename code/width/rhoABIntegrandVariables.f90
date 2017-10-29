!******************************************************************************
!****m* /rhoABIntegrandVariables
! NAME
! module rhoABIntegrandVariables
! NOTES
! This module declares some global variables for the integration routines
! ('rho_AB_Integrand') of baryon and meson decays into unstable particles.
! It is necessary to NOT declare 'rho_AB_Integrand' as  an internal function,
! since we want to have it as an argument in procedure call.
! Therefore we have to introduce this module as an interface of global
! variables between "semiStableFinalState" and "rho_ab_integrand".
!******************************************************************************
module rhoABIntegrandVariables

  implicit none

  Integer  :: idUnstable_copy  ! ID of unstable decay product
  Integer  :: L_copy           ! angular momentum of final state

  real :: massStable, massUnstable, gammaUnstable ! properties of the decay products
  real :: srts

end module rhoABIntegrandVariables
