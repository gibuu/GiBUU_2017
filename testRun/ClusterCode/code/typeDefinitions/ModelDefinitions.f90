!***************************************************************************
!****m* /ModelDefinitions
! NAME
! module ModelDefinitions
! PURPOSE
! Declare constants which define the different models of fragmentation:
! * 1 -- Phase space coalescence model
! * 2 -- statistical multifragmentation model
! NOTES
! * Model 1 contains also optionally an evaporation procedure of excited 
!   fragments, however, Fermi break-up is not taken into account.
! * Model 2 contains everything! Fermi break-up and sequential evaporation. 
!   This Statistical Multifragmentation Model (SMM) has been provided by 
!   Alexandre Botvina (GSI).
!***************************************************************************
module ModelDefinitions

  integer, parameter :: Coalescence = 1
  integer, parameter :: SMM = 2

end module ModelDefinitions
