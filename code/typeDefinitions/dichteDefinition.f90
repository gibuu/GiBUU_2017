!******************************************************************************
!****m* /dichteDefinition
! PURPOSE
! Defines the type "dichte".
!******************************************************************************
module dichteDefinition
  implicit none
  private

  !****************************************************************************
  !****t* dichteDefinition/dichte
  ! SOURCE
  !
  type, public :: dichte
     real,dimension(0:3) :: baryon = 0.   ! full baryon flux  in units 1/fm**3
     real,dimension(0:3) :: proton = 0.   ! full proton flux  in units 1/fm**3
     real,dimension(0:3) :: neutron= 0.   ! full neutron flux in units 1/fm**3
  end type dichte
  !****************************************************************************
  !****************************************************************************

  !Interfaces:
  Interface operator(+)! "+" for type(density)
     module procedure plus
  end Interface

  Interface operator(*)! "*" for type(density) and a real
     module procedure times1
!      module procedure times2
  end Interface

  public :: operator(+), operator(*)

contains

  function plus(x,y)     !plus for densities
    type(dichte),intent (in) :: x,y
    type(dichte)  :: plus
    plus%baryon=x%baryon+y%baryon
    plus%neutron=x%neutron+y%neutron
    plus%proton=x%proton+y%proton
  end function plus

  function times1(x,r)    !scalar multiplication for densities
    type(dichte),intent (in) :: x
    real,intent(in) :: r
    type(dichte)  :: times1
    times1%baryon=x%baryon*r
    times1%neutron=x%neutron*r
    times1%proton=x%proton*r
  end function times1

!   function times2(r,x)    !scalar multiplication for densities
!     type(dichte),intent (in) :: x
!     real,intent(in) :: r
!     type(dichte)  :: times2
!     times2%baryon=x%baryon*r
!     times2%neutron=x%neutron*r
!     times2%proton=x%proton*r
!   end function times2

end module dichteDefinition
