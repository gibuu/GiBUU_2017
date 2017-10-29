!******************************************************************************
!****m* /PILCollected
! NAME
! module PILCollected
! PURPOSE
! Provides some routines cleaning up or resetting PIL modules. We call
! this module "Collected", because it just provides abbreviational calls
! for all other modules: Call it here, and it will be called elsewhere!
!
! The lists 'collected' here are:
! * PIL_FormInfo
! * PIL_mesonMom
! * PIL_nLead
! * PIL_rho0Dec
! * PIL_rhoDiff
!
! Even if you do not use one of the modlues listed here, you may call the
! restoring routines provided by this module.
!******************************************************************************
module PILCollected

  use PIL_FormInfo
  use PIL_mesonMom
  use PIL_nLead
  use PIL_rho0Dec
  use PIL_rhoDiff
  use PIL_omegaDec
  use PIL_freezeout
  implicit none

  private

  public :: PILCollected_DeAllocate
  public :: PILCollected_ZERO

contains

  !****************************************************************************
  !****s* PILCollected/PILCollected_DeAlloc
  ! NAME
  ! subroutine PILCollected_DeAlloc
  ! PURPOSE
  ! Deallocate the memory of all lists
  !****************************************************************************
  subroutine PILCollected_DeAllocate()

    call PIL_FormInfo_DeAllocate
    call PIL_mesonMom_DeAllocate
    call PIL_nLead_DeAllocate
    call PIL_rho0Dec_DeAllocate
    call PIL_rhoDiffractive_DeAllocate
    call PIL_omegaDec_Deallocate
    call PIL_freezeout_Deallocate

  end subroutine PILCollected_DeAllocate


  !****************************************************************************
  !****s* PILCollected/PILCollected_ZERO
  ! NAME
  ! subroutine PILCollected_ZERO()
  ! PURPOSE
  ! Reset all lists by setting the counter of stored information to 0.
  ! No allocation or deallocation of memory happens.
  !****************************************************************************
  subroutine PILCollected_ZERO()

    call PIL_FormInfo_ZERO
    call PIL_mesonMom_ZERO
    call PIL_nLead_ZERO
    call PIL_rho0Dec_ZERO
    call PIL_rhoDiffractive_ZERO
    call PIL_omegaDec_Zero
    call PIL_freezeout_Zero

  end subroutine PILCollected_ZERO


end module PILCollected
