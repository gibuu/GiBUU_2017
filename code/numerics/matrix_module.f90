!******************************************************************************
!****m* /matrix
! NAME
! module matrix
!
! PURPOSE
! This module defines functions which are connected to matrices: Matrix multiplication, ....
!******************************************************************************
module matrix_module

  implicit none
  private

  ! two-dimensional unit matrix
  real, dimension(1:2,1:2), parameter :: unit2 = reshape((/ 1.,0., &
                                                            0.,1. /),(/2,2/))

  ! four-dimensional unit matrix
  real, dimension(0:3,0:3), parameter :: unit4 = reshape((/ 1.,0.,0.,0., &
                                                            0.,1.,0.,0., &
                                                            0.,0.,1.,0., &
                                                            0.,0.,0.,1. /),(/4,4/))

  public:: matrixMult,unit2,unit4,dagger,trace,printMatrix


  !****************************************************************************
  !****s* matrix/matrixMult
  ! NAME
  ! function matrixMult(a,b,c,...) result(matrix)
  ! PURPOSE
  ! Evaluates a x b x c x ...   for matrices
  ! INPUTS
  ! * a,b,c : matrices a,b,c -- this interface excepts a number of matrices between 2 and 9
  ! * The matrices should be complex
  ! OUTPUT
  ! * complex, dimension(0:3,0:3) :: matrix
  !****************************************************************************
  Interface matrixMult
     Module Procedure matrixMult_C2,matrixMult_C3,matrixMult_C4,matrixMult_C5,matrixMult_C6,matrixMult_C7, &
                      matrixMult_C8,matrixMult_C9,matrixMult_C10 !,matrixMult_C2_vec
  End Interface

contains


!   function matrixMult_C2_vec(a,b) result(matrix)
!     complex, dimension(0:3) :: matrix
!     complex, dimension(0:3,0:3),intent(in) :: a
!     complex, dimension(0:3),intent(in) :: b
!     integer :: i,k!,j
!     do i=0,3
!        matrix(i)=0.
!        do k=0,3
!           matrix(i)=matrix(i)+a(i,k)*b(k)
!        end do
!     end do
!   end function matrixMult_C2_vec


  function matrixMult_C2(a,b) result(matrix)
    complex, dimension(0:3,0:3) :: matrix
    complex, dimension(0:3,0:3),intent(in) :: a,b
    matrix=MatMul(a,b)
  end function matrixMult_C2

  function matrixMult_C3(a,b,c) result(matrix)
    complex, dimension(0:3,0:3) :: matrix
    complex, dimension(0:3,0:3),intent(in) :: a,b,c
    matrix=MatMul(a,MatMul(b,c))
  end function matrixMult_C3

  function matrixMult_C4(a,b,c,d) result(matrix)
    complex, dimension(0:3,0:3) :: matrix
    complex, dimension(0:3,0:3),intent(in) :: a,b,c,d
    matrix=matMul(a,MatrixMult_C3(b,c,d))
  end function matrixMult_C4

  function matrixMult_C5(a,b,c,d,e) result(matrix)
    complex, dimension(0:3,0:3) :: matrix
    complex, dimension(0:3,0:3),intent(in) :: a,b,c,d,e
    matrix=MatMul(a,MatrixMult_C4(b,c,d,e))
  end function matrixMult_C5

  function matrixMult_C6(a,b,c,d,e,f) result(matrix)
    complex, dimension(0:3,0:3) :: matrix
    complex, dimension(0:3,0:3),intent(in) :: a,b,c,d,e,f
    matrix=MatMul(a,MatrixMult_C5(b,c,d,e,f))
  end function matrixMult_C6

  function matrixMult_C7(a,b,c,d,e,f,g) result(matrix)
    complex, dimension(0:3,0:3) :: matrix
    complex, dimension(0:3,0:3),intent(in) :: a,b,c,d,e,f,g
    matrix=MatMul(a,MatrixMult_C6(b,c,d,e,f,g))
  end function matrixMult_C7

  function matrixMult_C8(a,b,c,d,e,f,g,h) result(matrix)
    complex, dimension(0:3,0:3) :: matrix
    complex, dimension(0:3,0:3),intent(in) :: a,b,c,d,e,f,g,h
    matrix=MatMul(a,MatrixMult_C7(b,c,d,e,f,g,h))
  end function matrixMult_C8

  function matrixMult_C9(a,b,c,d,e,f,g,h,i) result(matrix)
    complex, dimension(0:3,0:3) :: matrix
    complex, dimension(0:3,0:3),intent(in) :: a,b,c,d,e,f,g,h,i
    matrix=MatMul(a,MatrixMult_C8(b,c,d,e,f,g,h,i))
  end function matrixMult_C9

  function matrixMult_C10(a,b,c,d,e,f,g,h,i,j) result(matrix)
    complex, dimension(0:3,0:3) :: matrix
    complex, dimension(0:3,0:3),intent(in) :: a,b,c,d,e,f,g,h,i,j
    matrix=MatMul(a,MatrixMult_C9(b,c,d,e,f,g,h,i,j))
  end function matrixMult_C10



  !****************************************************************************
  !****f* matrix/dagger
  ! NAME
  ! function dagger(a) result(a_dagger)
  ! PURPOSE
  ! Evaluates a^Dagger=Transpose(Conj(a))
  ! INPUTS
  ! * complex, dimension(0:3,0:3) :: a
  ! OUTPUT
  ! * complex, dimension(0:3,0:3) :: a_dagger
  !****************************************************************************
  function dagger(a) result(a_dagger)
    complex, dimension(0:3,0:3),intent(in) :: a
    complex, dimension(0:3,0:3) :: a_dagger
    integer :: mu,nu
    a_dagger=transpose(a)
    do mu=0,3
       do nu=0,3
          a_dagger(mu,nu) =conjg(a_dagger(mu,nu))
       end do
    end do
  end function dagger


  !****************************************************************************
  !****f* matrix/trace
  ! NAME
  ! complex function trace(a)
  ! PURPOSE
  ! Evaluates Trace (a_{mu nu})
  ! INPUTS
  ! * real, dimension(0:3,0:3) :: a
  ! OUTPUT
  ! * complex
  !****************************************************************************
  function trace(a) result(c)
    complex :: c
    complex, intent(in),  dimension(0:3,0:3) :: a
    integer :: mu
    c=0.
    do mu=0,3
       c=c+a(mu,mu)
    end do
  end function trace


  !****************************************************************************
  !****s* matrix/printMatrix
  ! NAME
  ! subroutine printMatrix(A)
  ! PURPOSE
  ! * Prints the input matrix
  ! INPUTS
  ! * real, dimension(0:3,0:3) :: A
  ! OUTPUT
  ! * NONE
  !****************************************************************************
  subroutine printMatrix(A)
    integer :: j
    complex, intent(in),  dimension(0:3,0:3) :: a
    do j=0,3
       write(*,'(4("(",2G18.9,")"))') A(j,:)
    end do
  end subroutine printMatrix

end module matrix_module
