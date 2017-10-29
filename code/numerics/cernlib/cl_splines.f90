!******************************************************************************
!****m* /cl_splines
! NAME
! module cl_splines
!
! PURPOSE
! This module implements cubic splines according to cernlib (E211). See also
! Documentation for the routine "E211" in cernlib documentation.
!******************************************************************************
module cl_splines
  implicit none
  private
  character(30) :: modulename='cl_splines.f90'


  !****************************************************************************
  !****t* cl_splines/tspline
  ! PURPOSE
  ! Type which contains all information of the splines.
  ! SOURCE
  !
  type tspline
     real, dimension(:)  , Allocatable :: X
     real, dimension(:,:), Allocatable :: Y
     real, dimension(:,:), Allocatable :: A,B,C,D
     integer :: splineMode=1
     logical :: isInitialized=.false.
  end type tspline
  !
  !****************************************************************************

  public :: tspline, cl_initSpline, cl_spline, cl_error, cl_cleanupSpline

contains


  !****************************************************************************
  !****f* cl_splines/cl_initSpline
  ! NAME
  ! function cl_initSpline(x,y,splineMode) result(s)
  !
  ! PURPOSE
  ! Initializes a spline.
  !
  ! INPUTS
  ! * real, dimension(:),intent(in) :: x,y  -- x and y values of the data points. X and y must have same size!
  ! * integer,intent(in), optional :: splineMode --  1=second derivative of spline vanishes at end-points,
  !   2=second derivative is constant in between first (last) and second (last but one) point
  !
  ! OUTPUT
  ! * type(tspline) :: s -- an initialized spline
  !
  !****************************************************************************
  function cl_initSpline(x,y,splineMode) result(s)
    implicit none
    type(tspline) :: s
    real, dimension(:),intent(in) :: x,y
    integer,intent(in), optional :: splineMode

    integer :: n,i

    ! Interface to original cernlib subroutine:
    Interface
       SUBROUTINE RCSPLN(N,X,M,Y,NDIM,MODE,A,B,C,D)
         DIMENSION X(0:*),Y(0:NDIM,*)
         DIMENSION A(NDIM,*),B(NDIM,*),C(NDIM,*),D(NDIM,*)
       end SUBROUTINE RCSPLN
    end Interface

    if (size(x,dim=1).ne.size(y,dim=1)) then
       write(*,*) 'Error in module ',modulename
       write(*,*) 'size(x) not equal size(y). STOP!'
       stop
    end if

    if (s%isInitialized) then
       write(*,*) 'Error in module ',modulename
       write(*,*) 'Spline is already initialized. STOP!'
       stop
    end if


    do i=lbound(x,dim=1),ubound(x,dim=1)-1
       if (x(i).ge.x(i+1)) then
          write(*,*) 'Error in module ',modulename
          write(*,*) 'Vector x in non-monotonic', x
          write(*,*) 'STOP!'
          stop
       end if
    end do

    n=size(x,dim=1)

    allocate(s%x(0:n-1))
    allocate(s%y(0:n-1,1:1))

    allocate(s%A(1:n-1,1:1))
    allocate(s%B(1:n-1,1:1))
    allocate(s%C(1:n-1,1:1))
    allocate(s%D(1:n-1,1:1))

    s%x(:)=x
    s%y(:,1)=y

    if (present(splineMode)) s%splineMode = splineMode

    call RCSPLN(n-1,s%x,1,s%y,n-1,s%splineMode,s%a,s%b,s%c,s%d)

    s%isInitialized=.true.

  end function cl_initSpline


  !****************************************************************************
  !****f* cl_splines/cl_cleanupSpline
  ! NAME
  ! subroutine cl_cleanupSpline(s)
  !
  ! PURPOSE
  ! Deallocates all the memory that the spline occupied.
  !
  ! INPUTS
  ! * type(tspline) :: s -- an initialized spline
  !
  !****************************************************************************
  subroutine cl_cleanupSpline(s)
    implicit none
    type(tspline) :: s
    if (.not. s%isInitialized) return
    deallocate(s%x, s%y, s%A, s%B, s%C, s%D)
    s%isInitialized = .false.
  end subroutine


  !****************************************************************************
  !****f* cl_splines/cl_spline
  ! NAME
  ! real function cl_spline(s,x,successFlag,errorCode)
  !
  ! PURPOSE
  ! Evaluates the function value of a spline "s" at position "x" .
  !
  ! INPUTS
  ! * type(tspline) :: s             -- initialized spline
  ! * real   , intent(in)  :: x      -- position at which the spline shall be evaluated
  !
  ! OUTPUT
  ! * real :: value of the spline at x
  !
  ! * logical, intent(out) :: successFlag !.true. if spline was evaluated within bounds
  !
  ! * integer, intent(out) :: errorCode
  !
  ! NOTES
  ! * errorCode= 0 -> no error
  ! * errorCode=-1 -> x is lower   than any data point used for the spline
  ! * errorCode= 1 -> x is greater than any data point used for the spline
  !
  !****************************************************************************
  real function cl_spline(s,x,successFlag,errorCode)
    implicit none
    type(tspline) :: s
    real   , intent(in)  :: x
    logical, intent(out) :: successFlag
    integer, intent(out) :: errorCode

    integer :: i

    successFlag=.false.


    if (.not.s%isInitialized) then
       write(*,*) 'Error in module ',modulename
       write(*,*) 'Spline is not yet initialized. STOP!'
       stop
    end if

    indexLoop: do i=lbound(s%x,dim=1)+1,ubound(s%x,dim=1)
       if (x.ge.s%x(i-1).and. x.le.s%x(i)) then
          successFlag=.true.
          errorCode=0
          cl_spline=eval(i)
          exit indexLoop
       end if
    end do indexLoop


    if (.not.successFlag) then
       if (x.ge.s%x(ubound(s%x,dim=1))) then
          errorCode=1
          cl_spline=eval(ubound(s%a,dim=1))
       else if (x.le.s%x(lbound(s%x,dim=1))) then
          errorCode=-1
          cl_spline=eval(lbound(s%a,dim=1))
       else
          write(*,*) 'Unexpected error in module ',x,s%x
          write(*,*) 'STOP!'
          stop
       end if
    end if
  contains
    real function eval(k)
      implicit none
      integer, intent(in) :: k
      eval=s%a(k,1)+s%b(k,1)*(x-s%x(k-1))+s%c(k,1)*(x-s%x(k-1))**2+s%d(k,1)*(x-s%x(k-1))**3
    end function eval
  end function cl_spline

  !****************************************************************************
  !****s* cl_splines/cl_error
  ! NAME
  !    subroutine cl_error(error,name,x)
  !
  ! PURPOSE
  ! Implements output of error messages
  !
  ! INPUTS
  ! * integer     , intent(in) :: error -- type of error
  ! * character(*), intent(in) :: name  -- identifier (Where did the error happen)
  ! * real        , intent(in) :: x     -- value where spline should be evaluated
  !
  ! NOTES
  ! * errorCode= 0 -> no error
  ! * errorCode=-1 -> x is lower   than any data point used for the spline
  ! * errorCode= 1 -> x is greater than any data point used for the spline
  !
  !****************************************************************************
  subroutine cl_error(error,name,x)
    implicit none
    integer     , intent(in) :: error
    character(*), intent(in) :: name
    real        , intent(in) :: x
    select case (error)
    case (0)
       write(*,*) 'No error'
    case (-1)
       write(*,*) 'Underflow of spline in', name
       write(*,*) 'x=', x
    case (1)
       write(*,*) 'Overflow of spline in', name
       write(*,*) 'x=', x
    end select
  end subroutine cl_error

end module cl_splines
