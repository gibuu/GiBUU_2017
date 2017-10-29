!******************************************************************************
!****m* /derivatives
! NAME
! module derivatives
! PURPOSE
! Module which includes derivative and finite difference routines.
!******************************************************************************
module derivatives

  implicit none
  private

  public :: finiteDifference, derivative

contains

  !****************************************************************************
  !****f* derivatives/finiteDifference
  ! NAME
  ! function finiteDifference(y,dx,order,scheme) result(derivative)
  ! PURPOSE
  ! Evaluates the derivatives using finite differences.
  ! INPUTS
  ! * real   ,intent(in), dimension(-2:2) :: y ! function value at y(x-2*dx),...,y(x+2*dx)
  ! * real   ,intent(in) :: dx                 ! delta x
  ! * integer,intent(in) :: order              ! 1=first order, 2=second order
  ! * integer,intent(in) :: scheme             ! negative=backward difference,0=central difference,positive=forward difference
  ! OUTPUT
  ! * real :: derivative
  !****************************************************************************
  function finiteDifference(y,dx,order,scheme) result(derivative)
    real :: derivative
    real   ,intent(in), dimension(-2:2) :: y
    real   ,intent(in) :: dx
    integer,intent(in) :: order
    integer,intent(in) :: scheme

    select case (order)
    case (1)
       ! First order derivative
       if (scheme.lt.0) then
          ! Backward difference
          derivative=(y(0)-y(-1))/dx
       else if (scheme.gt.0) then
          ! Forward difference
          derivative=(y(1)-y(0))/dx
       else
          ! Central difference
          derivative=(y(1)-y(-1))/dx/2.
       end if
    case (2)
       ! Second order derivative
       if (scheme.lt.0) then
          ! Backward difference
          derivative= (3.*y(0)-4.*y(-1)+y(-2))/(2.*dx)
       else if (scheme.gt.0) then
          ! Forward difference
          derivative= (-3.*y(0)+4.*y(1)-y(2))/(2.*dx)
       else
          ! Central difference
          derivative=(y(-2)-8.*y(-1)+8.*y(1)-y(2))/(12.*dx)
       end if
    case default
       write(*,*) 'Order not yet implemented in derivative', order
       stop
    end select

  end function finiteDifference


  !****************************************************************************
  !****f* derivatives/derivative
  ! NAME
  ! real function derivative(f,x,dx,order,scheme)
  ! PURPOSE
  ! Evaluates the derivative of "f" at position "x".
  ! INPUTS
  ! * real, external     :: f
  ! * real,intent(in)    :: x
  ! * real,intent(in)    :: dx  -- step size
  ! * integer,intent(in) :: order        --  1=first order, 2=second order
  ! * integer,intent(in) :: scheme       --  negative=backward difference scheme,0=central difference,positive=forward difference
  !****************************************************************************
  real function derivative(f,x,dx,order,scheme)
    real, external     :: f
    real,intent(in)    :: x,dx
    integer,intent(in) :: order, scheme
    integer :: i
    real, dimension(-2:2) :: y
    y=0.
    if (scheme.lt.0) then
       ! backward scheme
       do i=-order,0
          y(i)=f(x+float(i)*dx)
       end do
    else if (scheme.gt.0) then
       ! forward scheme
       do i=0,order
          y(i)=f(x+float(i)*dx)
       end do
    else
       ! central scheme
       do i=-order,order
          if (i.eq.0) cycle
          y(i)=f(x+float(i)*dx)
       end do
    end if
    derivative=finiteDifference(y,dx,order,scheme)
  end  function derivative


end module derivatives
