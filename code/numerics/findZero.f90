!******************************************************************************
!****m* /findZero
! NAME
! module findZero
!
! PURPOSE
! * This module includes routines to find roots of functions, i.e. the x for which f(x)=0.
!
! METHODS
! * newton_findZero
!******************************************************************************
module findZero


contains



  !****************************************************************************
  !****f* findZero/newton_findZero
  ! NAME
  ! function newton_findZero(f, start,accuracy,maxSteps,success) result(zero)
  !
  ! PURPOSE
  ! * This function evaluates the root of the function "f", i.e. the x for which f(x)=0.
  ! * Uses a simple Newton algorithm. Not suited for several roots!!
  !
  ! INPUTS
  ! * real, external :: f                      -- external function
  ! * real, intent(in) :: start, accuracy      -- starting value and accuracy for the zero
  ! * integer, intent(in) :: maxSteps          -- maximal number of iterations
  !
  ! Output
  ! * The return value of the function is the root itself
  ! * logical, intent(out) :: success         -- "true" if the search was succesful, else "false"
  !
  !****************************************************************************
  function newton_findZero(f, start,accuracy,maxSteps,success) result(zero)
    implicit none
    real, external :: f
    real, intent(in) :: start, accuracy
    integer, intent(in) :: maxSteps
    logical, intent(out) :: success
    real :: zero
    real, dimension(1:2) :: x
    real :: y
    character(100) :: name='module findZero: function newton_findZero'
    real :: dx,d
    integer :: counter


    zero=0.
    success=.false.


    ! Set starting values
    x(1)=start

    do counter=1,maxSteps
       y=f(x(1))
       ! Adapt sensible values for dx via the value of x(1)
       dx=abs(x(1)*10E-5)
       ! Perform Newton Algorithm
       d= deriv(x(1))
       x(2)=x(1)-y/d
       y=f(x(2))
       if (abs(y).lt.accuracy) then
          success=.true.
          zero=x(2)
          return
       end if
       if (abs(d).lt.epsilon(d)) then
          ! Derivative is zero
          write(*,*) name
          write(*,*) 'Error: derivative=0! No SUCCESS!!! '
          write(*,*) 'derivative=',d
          write(*,*) 'counter=',counter
          return
       end if
       x(1)=x(2)
    end do

    write(*,*) name
    write(*,*) 'Error: counter > maxSteps ! No SUCCESS!!! '


  contains
    real function deriv(x)
      implicit none
      real :: x
      deriv= (f(x+dx)-f(x-dx))/(2.*dx)
    end function deriv
  end function newton_findZero





  !****************************************************************************
  !****f* findZero/bisection_findZero
  ! NAME
  ! function bisection_findZero(f, a,b,accuracy,maxSteps,success) result(zero)
  !
  ! PURPOSE
  ! * This function evaluates the root of the function "f" on the interval (a,b) , i.e. the x for which f(x)=0.
  ! * Uses a simple bisection algorithm. Not suited for several roots!!
  ! * Works best if f(a)*f(b)<0. If not, then the computation time is increasing.
  !
  ! INPUTS
  ! * real, external :: f                      -- external function
  ! * real, intent(in) :: a,b                  -- defines the search range a<x<b
  ! * real, intent(in) :: accuracy             -- accuracy for the zero
  ! * integer, intent(in) :: maxSteps          -- maximal number of iterations
  !
  ! Output
  ! * The return value of the function is the root itself
  ! * logical, intent(out) :: success         -- "true" if the search was succesful, else "false"
  !
  !****************************************************************************
  function bisection_findZero(f, a,b,accuracy,maxSteps,success) result(zero)
    implicit none
    real, external :: f
    real, intent(in) :: a,b, accuracy
    integer, intent(in) :: maxSteps
    logical, intent(out) :: success
    real :: zero
    real, dimension(1:3) :: x,y
    character(60) :: name='module findZero: function bisection_findZero'
    !integer :: counter
    integer :: n,i
    logical :: successNegative
    logical :: debug=.false.
    zero=0.
    success=.false.

    if (debug) write(*,*) 'x,y,a,b=',x,y,a,b


    ! Set starting values
    x(1)=a
    x(2)=b
    y(1)=f(x(1))
    y(2)=f(x(2))

    if (debug) then
       write(*,*) 'x1=',x
       write(*,*) 'y1=',y
    end if


    ! First establish that y(1)*y(2)<0
    successNegative=.false.
    if (y(1)*y(2).gt.0) then
       ! Slowly increase the grid spacings to find a sign change in f(x):
       findNegative_loop: do n=2, maxSteps,3 ! n= Number of grid points
          do i=1,n-1  ! loop over all gridpoints
             x(3)=x(1)+(x(2)-x(1))/float(n)*float(i)
             y(3)=f(x(3))
             if (y(1)*y(3).lt.0) then
                ! The function changed sign in between x=x(1)+(x(2)-x(1))/float(n)*float(i)
                ! and x(1)+(x(2)-x(1))/float(n)*float(i-1)
                ! => Root must be in the middle
                successNegative=.true.
                x(1)=x(3)-(x(2)-x(1))/float(n)
                x(2)=x(3)
                y(1)=f(x(1))
                y(2)=f(x(2))
             end if
             if ( successNegative ) exit findNegative_loop
          end do
       end do findNegative_loop
    else
       successNegative=.true.
    end if


    if (.not. successNegative ) then
       write(*,*) name
       write(*,*) 'Could not find y(1)*y(2)<0!!'
       return
    end if

    ! Do bisection method
    findRoot_loop: do  i=1, maxSteps
       x(3)=(x(2)+x(1))/2.
       y(3)=f(x(3))
       if (debug) then
          write(*,'(A,3E15.4)') 'x im loop (1)=',x
          write(*,'(A,3E15.4)') 'y im loop (1)=',y
       end if
       if (y(3)*y(1).lt.0) then
          ! root in between x(1) and x(3)
          x(2)=x(3)
          y(2)=f(x(2))
       else if (y(3)*y(2).lt.0) then
          ! root in between x(3) and x(2)
          x(1)=x(3)
          y(1)=f(x(1))
       else
          write(*,*) name
          write(*,*) 'Error in bisection algorithm'
          stop
       end if
       if (debug) then
          write(*,'(A,3E15.4)') 'x im loop (2)=',x
          write(*,'(A,3E15.4)') 'y im loop (2)=',y
       end if
       if (abs(y(3)).lt.accuracy) then
          success=.true.
          zero=x(3)
          return
       end if
    end do findRoot_loop
    write(*,*) name
    write(*,*) 'steps > maxsteps. Could not reach accuracy!!'
  end function bisection_findZero


end module findZero
