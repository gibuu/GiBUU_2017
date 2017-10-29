!******************************************************************************
!****m* /Coll_gammaN_exclPi
! NAME
! module Coll_gammaN_exclPi
!
! PURPOSE
! This module contains all routines necessary to produce /virtual) photon
! induced exclusive pion production event.
!******************************************************************************
module Coll_gammaN_exclPi

  IMPLICIT NONE

  private

  public :: DoColl_gammaN_exclPi

contains

  !****************************************************************************
  !****s* Coll_gammaN_exclPi/DoColl_gammaN_exclPi
  ! NAME
  ! subroutine DoColl_gammaN_exclPi(eNev,ExclPiCharge,flagOK,outPart, XS_tot)
  !
  ! PURPOSE
  ! generate a exclusive pion production event
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev     -- electron nucleon kinematics
  ! * integer :: ExclPiCharge -- charge of pion to be produced
  !
  ! OUTPUT
  ! * logical   :: FlagOK  -- .true. if event was successfull
  ! * real      :: XS_Tot  -- total cross section
  ! * type(particle), dimension(:)  :: OutPart -- Final State particle vector
  ! NOTES
  ! * The returned Cross Section 'XS_Tot' is
  !     dsigma/dE'dOmega (in mub/MeV)
  !   in the target rest frame
  !****************************************************************************
  subroutine DoColl_gammaN_exclPi(eNev,ExclPiCharge,flagOK,outPart, XS_tot)
    use particleDefinition
    use eN_eventDefinition, only: electronNucleon_event

    type(electronNucleon_event),  intent(in)  :: eNev
    integer,                      intent(in)  :: ExclPiCharge
    type(particle), dimension(:), intent(out) :: OutPart
    real,                         intent(out) :: XS_Tot
    logical,                      intent(out) :: FlagOK


    integer :: NucCharge

    ! reset all output:

    OutPart%ID = 0
    XS_tot = 0
    flagOK = .false.

    ! check charge:
    NucCharge = eNev%nucleon%charge - ExclPiCharge
    if ((NucCharge.lt.0).or.(NucCharge.gt.1)) return


    ! now do some stuff...


  end subroutine DoColl_gammaN_exclPi




end module Coll_gammaN_exclPi
