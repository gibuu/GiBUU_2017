!******************************************************************************
!****m* /leptonicID
! NAME
! module leptonicID
!
! PURPOSE
! define some constants for characterizing lepton induced reactions
!******************************************************************************

module leptonicID
  !****************************************************************************
  !****t* leptonicID/Constants
  ! PURPOSE
  ! values for process_ID; reactions induced by leptons are >0, while
  ! reactions induced by anti-leptons are <0
  ! SOURCE
  !
  integer, parameter :: EM = 1
  integer, parameter :: CC = 2
  integer, parameter :: NC = 3

  integer, parameter :: antiEM = -1
  integer, parameter :: antiCC = -2
  integer, parameter :: antiNC = -3


  !values for charge
  integer, parameter :: proton = 1
  integer, parameter :: neutron = 0
  !****************************************************************************
end module leptonicID
