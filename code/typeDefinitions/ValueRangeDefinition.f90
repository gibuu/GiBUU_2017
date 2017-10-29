!******************************************************************************
!****m* /ValueRangeDefinition
! NAME
! module ValueRangeDefinition
! PURPOSE
! Defines a type which holds 3 values: min, max and delta.
!
! NOTES
! * used in HiLeptonAnalysis
! * may be used in all histogram definitions etc.
!******************************************************************************
module ValueRangeDefinition

  !****************************************************************************
  !****t* ValueRangeDefinition/ValueRange
  ! NAME
  ! type ValueRange
  ! PURPOSE
  ! Store a minimal value, a maximal value and an increment in one structure.
  ! SOURCE
  !
  type ValueRange
     real :: l ! = lower bound
     real :: u ! = upper bound
     real :: d ! = increment, delta
  end type ValueRange
  !
  ! NOTES
  ! you can set values by e.g.:
  ! * type(ValueRange) :: xRange
  ! * type(ValueRange) :: yRange = ValueRange(5.0,7.0,0.2)
  ! * type(ValueRange) :: zRange
  ! * xRange%l=3
  ! * zRange = ValueRange(3.0,7.0,0.2)
  !****************************************************************************

end module ValueRangeDefinition
