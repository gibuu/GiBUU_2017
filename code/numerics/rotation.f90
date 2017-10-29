!******************************************************************************
!****m* /rotation
! NAME
! module rotation
!
! PURPOSE
! This module defines rotation matrices and evaluates the spherical
!coordinates of a vector.
!******************************************************************************
module rotation
  implicit none
  private

  public :: rotateYZ, get_phi_theta, rotateZY, rotateTo, rotateFrom

contains


  !****************************************************************************
  !****f* rotation/rotateTo
  ! NAME
  ! function rotateTo(targ, v_in)
  ! PURPOSE
  !
  ! Rotate some vector with the same matrix as one would use to rotate the
  ! z-axis onto a vector 'targ'.
  !
  ! INPUTS
  ! * real, dimension(1:3) :: targ -- Target vector
  ! * real, dimension(1:3) :: v_in -- Vector before rotation
  ! OUTPUT
  ! * real, dimension(1:3) :: v_out -- Vector after rotation
  !
  ! USAGE
  ! * calling
  !     b = rotateTo( a, (/0.,0.,1./) )
  !   b is the unit-vector in direction of a
  !****************************************************************************
  function rotateTo(targ, v_in) result (v_out)
    real, intent(in),  dimension(1:3) :: targ, v_in
    real, dimension(1:3) :: v_out

    real :: theta, phi

    call get_phi_theta(targ(1:3), theta, phi)
    v_out(1:3) = rotateYZ(theta, phi, v_in(1:3))

  end function rotateTo

  !****************************************************************************
  !****f* rotation/rotateFrom
  ! NAME
  ! function rotateFrom(targ, v_in)
  ! PURPOSE
  !
  ! Rotates some vector with the same matrix as one would use to rotate the
  ! vector 'targ' onto the z-axis
  !
  ! INPUTS
  ! * real, dimension(1:3) :: targ -- Target vector
  ! * real, dimension(1:3) :: v_in -- Vector before rotation
  ! OUTPUT
  ! * real, dimension(1:3) :: v_out -- Vector after rotation
  !
  ! USAGE
  ! * calling
  !     b = rotateFrom( a, a )
  !   b is (/ 0., 0., |a| /)
  ! * calling
  !     b = rotateTo( a, (/ 0., 0., 1. /) )
  !     c = rotateFrom( a, b )
  !   c is (/ 0., 0., 1. /) again.
  !****************************************************************************
  function rotateFrom(targ, v_in) result (v_out)
    real, intent(in),  dimension(1:3) :: targ, v_in
    real, dimension(1:3) :: v_out

    real :: theta, phi

    call get_phi_theta(targ(1:3), theta, phi)
    v_out(1:3) = rotateYZback(theta, phi, v_in(1:3))

  end function rotateFrom


  !****************************************************************************
  !****f* rotation/rotateYZ
  ! NAME
  ! function rotateYZ(theta, phi, v_in, cosTheta)
  ! PURPOSE
  ! Does first a rotation around the y-axis with angle theta, then around
  ! z with angle phi.
  !
  ! Rotates the z-axis on a vector which has the spherical coordinates theta
  ! and phi.
  !
  ! INPUTS
  ! * real                 :: phi   -- Angle around z-axis
  ! * real                 :: theta -- Angle around y-Axis
  ! * real, dimension(1:3) :: v_in  -- Vector before rotation
  ! * real, OPTIONAL       :: cosTheta -- cos(theta)
  ! OUTPUT
  ! * real, dimension(1:3) :: v_out -- Vector after rotation
  ! NOTES
  ! if optional parameter cosTheta is given, cosT and sinT are calculated
  ! using this parameter instead of theta
  !****************************************************************************
  function rotateYZ(theta, phi, v_in, cosTheta) result (v_out)

    real, intent(in) :: theta, phi
    real, intent(in), dimension(1:3) :: v_in
    real, intent(in), optional :: cosTheta
    real, dimension(1:3) :: v_out

    real :: sinT,cosT,cosP,sinP

    if (present(cosTheta)) then
       cosT = cosTheta
       sinT = sqrt(max(1.-cosT**2, 0.))
    else
       cosT=cos(theta)
       sinT=sin(theta)
    end if

    cosP=cos(phi)
    sinP=sin(phi)

    v_out(1)=  cosT*cosP*v_in(1) -sinP*v_in(2) +sinT*cosP*v_in(3)
    v_out(2)=  cosT*sinP*v_in(1) +cosP*v_in(2) +sinT*sinP*v_in(3)
    v_out(3)= -sinT     *v_in(1)               +cosT     *v_in(3)

  end function rotateYZ

  !****************************************************************************
  !****f* rotation/rotateYZback
  ! NAME
  ! function rotateYZback(theta, phi, v_in, cosTheta)
  ! PURPOSE
  ! Does first a rotation around the z-axis with angle -theta, then around
  ! y with angle -phi.
  !
  ! Is the back rotation performed by rotateYZ
  !
  ! INPUTS
  ! * real                 :: phi   -- Angle around z-axis
  ! * real                 :: theta -- Angle around y-Axis
  ! * real, dimension(1:3) :: v_in  -- Vector before rotation
  ! * real, OPTIONAL       :: cosTheta -- cos(theta)
  ! OUTPUT
  ! * real, dimension(1:3) :: v_out -- Vector after rotation
  ! NOTES
  ! if optional parameter cosTheta is given, cosT and sinT are calculated
  ! using this parameter instead of theta
  !****************************************************************************
  function rotateYZback(theta, phi, v_in, cosTheta) result (v_out)

    real, intent(in) :: theta, phi
    real, intent(in), dimension(1:3) :: v_in
    real, intent(in), optional :: cosTheta
    real, dimension(1:3) :: v_out

    real :: sinT,cosT,cosP,sinP

    if (present(cosTheta)) then
       cosT = cosTheta
       sinT = sqrt(max(1.-cosT**2, 0.))
    else
       cosT=cos(theta)
       sinT=sin(theta)
    end if

    cosP=cos(phi)
    sinP=sin(phi)

    v_out(1)=  cosT*cosP*v_in(1) +cosT*sinP*v_in(2) -sinT*v_in(3)
    v_out(2)= -sinP     *v_in(1) +cosP     *v_in(2)
    v_out(3)=  sinT*cosP*v_in(1) +sinT*sinP*v_in(2) +cosT*v_in(3)


  end function rotateYZback

  !****************************************************************************
  !****f* rotation/rotateZY
  ! NAME
  ! function rotateZY(theta, phi, v_in, cosTheta)
  ! PURPOSE
  ! Does first rotation around z and thereafter around y.
  !
  ! As a result, e.g. the vector which had originally the spherical coordinates
  ! theta and phi is now rotated on the z-Axis. It's basically the inverse of
  ! rotateYZ.
  !
  ! INPUTS
  ! * real                 :: phi   -- Angle around z-axis [radian]
  ! * real                 :: theta -- Angle around y-Axis [radian]
  ! * real, dimension(1:3) :: v_in  -- Vector before rotation
  ! * real, OPTIONAL       :: cosTheta -- cos(theta)
  ! OUTPUT
  ! * real, dimension(1:3) :: v_out -- Vector after rotation
  ! NOTES
  ! if optional parameter cosTheta is given, cosT and sinT are calculated
  ! using this parameter instead of theta
  !****************************************************************************
  function rotateZY(theta, phi, v_in, cosTheta) result (v_out)

    real, intent(in) :: theta, phi
    real, intent(in), dimension(1:3) :: v_in
    real, intent(in), optional :: cosTheta
    real, dimension(1:3) :: v_out

    real :: sinT,cosT,cosP,sinP
    real, dimension(1:3,1:3) :: RotMatrix

    if (present(cosTheta)) then
       cosT = cosTheta
       sinT = sqrt(max(1.-cosT**2, 0.))
    else
       cosT=cos(theta)
       sinT=sin(theta)
    end if

    cosP=cos(phi)
    sinP=sin(phi)

    RotMatrix = reshape((/ cosT*cosP , cosT*sinP , -sinT , &
                           -sinP     , cosP      , 0.    , &
                           sinT*cosP , sinT*sinP , cosT  /) , &
                           (/3,3/) ,order=(/2,1/))

    v_out = Reshape(MatMul(RotMatrix,reshape(v_in,(/3,1/))),(/3/))

  end function rotateZY



  !****************************************************************************
  !****s* rotation/get_phi_theta
  ! NAME
  ! subroutine get_phi_theta(vector, theta, phi)
  ! PURPOSE
  ! Translate the input "vector" into spherical coordinates.
  ! Returns phi and theta of vector:
  ! * x-axis defines phi=0
  ! * z-axis defines theta=0
  ! INPUTS
  ! * real, dimension(1:3) :: vector --
  ! OUTPUT
  ! * real :: phi    -- phi   [radian], range: 0 <= phi   <= 2*pi
  ! * real :: theta  -- theta [radian], range: 0 <= theta <= pi
  !****************************************************************************
  subroutine get_phi_theta(vector, theta, phi)
    use constants, only: pi

    real, dimension(1:3), intent(in) :: vector
    real, intent(out) :: phi, theta
    real :: v2,cosT

    v2 = dot_product(vector,vector)
    if (v2>0) then
       cosT=vector(3)/sqrt(v2)
       theta=acos(cosT)
       phi=atan2(vector(2),vector(1))
       if (phi<0) phi=phi+2.*pi
    else
       phi=0.
       theta=0.
    end if

  end subroutine get_phi_Theta


end module rotation
