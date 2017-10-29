!******************************************************************************
!****m* /degRad_conversion
! NAME
! module degRad_conversion
!
! PURPOSE
! This module defines functions to convert degrees to radians and vice versa
!
!******************************************************************************


module degRad_conversion
Private

Public :: radian, degrees

contains
  !****************************************************************************
  !****f* degRad_conversion/radian
  ! NAME
  ! real function radian(x)
  ! PURPOSE
  ! Converts input angle from degrees to radian
  ! INPUTS
  ! * real  :: x
  ! OUTPUT
  ! * real  :: radian
  !****************************************************************************
  real function radian(x)
    ! Converts input angle from degree to radians
    use constants, only: pi
    implicit none
    real, intent(in) :: x ! input in units of degree
    radian=x/180.*pi
  end function radian

  !****************************************************************************
  !****f* degRad_conversion/degrees
  ! NAME
  ! real function degrees(x)
  ! PURPOSE
  ! Converts input angle from radians to degrees
  ! INPUTS
  ! * real  :: x
  ! OUTPUT
  ! * real  :: degrees
  !****************************************************************************
  real function degrees(x)
    ! Converts input angle from radians to degrees
    use constants, only: pi
    implicit none
    real, intent(in) :: x ! input in units of sr
    degrees=x*180./pi
  end function degrees


end module degRad_conversion
