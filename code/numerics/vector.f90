!******************************************************************************
!****m* /vector
! NAME
! module vector
!
! PURPOSE
! This module defines functions for vector handling
!
!******************************************************************************



module vector

  public :: theta_in, absVec, crossProduct, sphericalVector

contains


  !****************************************************************************
  !****f* vector/sphericalVector
  ! NAME
  ! function sphericalVector(theta,phi,rabs) result(r)
  ! PURPOSE
  ! Returns cartesian vector with spherical coordinates (theta,phi,rAbs)
  ! INPUTS
  ! * real,intent (in) :: theta [degrees], phi [degrees]
  ! * real,intent (in) :: rabs
  ! OUTPUT
  ! * real, dimension(1:3) :: r
  !****************************************************************************
  function sphericalVector(theta,phi,rabs) result(r)
    use degRad_conversion
    implicit none
    real,intent (in) :: theta, phi,rabs
    real, dimension(1:3) :: r
    r= rabs * (/sin(Radian(theta))*cos(radian(phi)),sin(Radian(theta))*sin(radian(phi)),cos(Radian(theta)) /)
  end function  sphericalVector

  !****************************************************************************
  !****f* vector/sphericalVector_radian
  ! NAME
  ! function sphericalVector_radian(theta,phi,rabs) result(r)
  ! PURPOSE
  ! Returns cartesian vector with spherical coordinates (theta,phi,rAbs)
  ! INPUTS
  ! * real,intent (in) :: theta [radian], phi [radian]
  ! * real,intent (in) :: rabs
  ! OUTPUT
  ! * real, dimension(1:3) :: r
  !****************************************************************************
  function sphericalVector_radian(theta,phi,rabs) result(r)
    implicit none
    real,intent (in) :: theta, phi,rabs
    real, dimension(1:3) :: r
    r= rabs * (/sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta) /)
  end function  sphericalVector_radian




  !****************************************************************************
  !****f* vector/theta_in
  ! NAME
  ! real function theta_in(a,b,switch)
  ! PURPOSE
  ! Calculates angle in between vector a and b.
  ! INPUTS
  ! * real ,dimension (1:3),intent(in) :: a,b
  ! * integer,intent(in), optional :: switch
  ! OUTPUT
  ! * real  :: theta_in
  ! NOTES
  ! * If no switch is given or switch=1, then units : [degrees]
  ! * If switch=2, then units : [radian]
  !****************************************************************************
  real function theta_in(a,b,switch)
    use degRad_conversion, only: degrees
    implicit none
    real ,dimension (1:3),intent(in) :: a,b
    integer,intent(in), optional :: switch
    real :: abs_ab


    abs_ab=sqrt(Dot_product(a,a)*Dot_product(b,b))

    if (present(switch)) then
       select case (switch)
       case (1)
          theta_in=degrees(acos(Dot_Product(a,b)/abs_ab))
       case (2)
          theta_in=acos(Dot_Product(a,b)/abs_ab)
       case default
          write(*,*) 'Invalid switch in theta_in, module vector:',switch
          write(*,*) 'STOP'
          stop
       end select
    else
       theta_in=degrees(acos(Dot_Product(a,b)/abs_ab))
    end if
  end function theta_in

  !****************************************************************************
  !****f* vector/absVec
  ! NAME
  ! real function absVec(a)
  ! PURPOSE
  ! Calculates the absolute value of a vector a.
  ! INPUTS
  ! * real ,dimension (:),intent(in) :: a
  !****************************************************************************
  real function absVec(a)
    implicit none
    real ,dimension (:),intent(in) :: a

    absVec=sqrt(Dot_product(a,a))

  end function absVec



  !****************************************************************************
  !****f* vector/crossProduct
  ! NAME
  ! function crossProduct(a,b)
  ! PURPOSE
  ! Calculates cross product of two 3-vectors: c= a x b = sum_{j,k} \epsilon_{ijk} a_j b_k
  ! INPUTS
  ! real ,dimension (1:3),intent(in) :: a,b
  ! OUTPUT
  ! real ,dimension (1:3)            :: crossProduct
  !****************************************************************************
  function crossProduct(a,b)
    implicit none
    real ,dimension (1:3),intent(in) :: a,b
    real ,dimension (1:3)            :: crossProduct

    crossProduct(1)=a(2)*b(3)-a(3)*b(2)
    crossProduct(2)=a(3)*b(1)-a(1)*b(3)
    crossProduct(3)=a(1)*b(2)-a(2)*b(1)


  end function crossProduct



end module vector
