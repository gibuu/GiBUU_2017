
!***************************************************************************
!****m* /photon_flux
! NAME
! module photon_flux
!
! PURPOSE
! This module defines functions for evaluating the photon_flux
!***************************************************************************

module photon_flux

  implicit none

contains 

  !*************************************************************************
  !****f* photon_flux/gamma
  ! NAME
  ! real function gamma(Ebeam,q,theta)
  ! PURPOSE
  ! * Calculates the photon flux of a given electron beam
  ! * Returns Gamma of Hand convention [in rest frame of nucleon!]
  ! INPUTS
  ! * real, dimension(0:3),intent(in) :: q -- photon 4-momemtum
  ! * real, intent(in)  :: theta           -- electron scattering angle
  ! * real, intent(in)  :: Ebeam           -- electron beam energy
  !*************************************************************************
  real function gamma(Ebeam,q,theta)
    use minkowski, only : SP
    use degRad_conversion, only : radian
    use constants, only : pi, alphaQED, mN

    real, dimension(0:3),intent(in) :: q
    real, intent(in)  :: theta ! electron scattering angle
    real, intent(in)  :: Ebeam ! beam energy
    
    real :: K, epsilon,QS,s

    QS=-SP(q,q)
    
    s=(mN+q(0))**2-dot_product(q(1:3),q(1:3))
    K=(s-mN**2)/2./mN

    epsilon=1./(1.+2.*dot_product(q(1:3),q(1:3))/QS*(tan(radian(theta)/2.))**2)

    Gamma=alphaQED/2./pi**2*(ebeam-q(0))/ebeam * K/QS /(1-epsilon)
    
  end function gamma


  
  


end module photon_flux
