!******************************************************************************
!****m* /propagators
! NAME
! module propagators
! PURPOSE
! * Provides propagators, i.e. two-point functions
!******************************************************************************
module propagators

  implicit none
  private

  public :: propagator_3_2_vac

contains

  !****************************************************************************
  !****f* propagators/propagator_3_2_vac
  ! NAME
  ! function propagator_3_2_vac(resID,p_res) result(propagator)
  !
  ! PURPOSE
  ! * Provides a spin=3/2 propagator G^{mu,nu} of a resonance in vacuum. Note that G^{mu nu} is a 4x4 matrix.
  !
  ! INPUTS
  ! * real, dimension(0:3),intent(in) :: p_res  ! resonance 4-momentum
  ! * integer :: resID                          ! ID of the resonance
  !
  ! OUTPUT
  ! complex,dimension(0:3,0:3,0:3,0:3), intent(out) ::propagator
  ! First two indices =mu,nu,  i.e.   (G^{mu,nu} )_ij=propagator(mu,nu,i,j)
  !****************************************************************************
  function propagator_3_2_vac(resID,p_res) result(propagator)
    use baryonWidth, only:FullWidthBaryon
    use particleProperties, only: hadron
    use minkowski, only: SP
    use spinProjector, only: spin32proj_tensor
    use constants, only: ii

    real, dimension(0:3), intent(in)    :: p_res
    integer             , intent(in)    :: resID

    complex,dimension(0:3,0:3,0:3,0:3)  ::propagator

    real :: Gamma,s, sqrtS

    s=SP(p_res,p_res)
    if (s.lt.0) then
       write(*,'(/,A,/,A,G15.5)') 'Error in propagator_3_2_vac', 'Mass^2 of resonance less than 0! mass^2=', s
       STOP 'stop'
    end if

    sqrtS=sqrt(s)
    ! Retrieving vacuum width:
    Gamma=FullWidthBaryon(resID,sqrtS)

    ! The propagator
    propagator= spin32proj_tensor(resID,p_res)/(S-cmplx(hadron(resID)%mass**2) + ii*sqrtS*Gamma)

  end function propagator_3_2_vac


end module propagators
