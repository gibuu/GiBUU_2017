!******************************************************************************
!****m* /sorting
! NAME
! module sorting
!
! PURPOSE
! This module defines some routines for sorting purposes
!******************************************************************************
module sorting

  implicit none
  private

  public :: indexx

contains

  !****************************************************************************
  !****s* sorting/indexx
  ! NAME
  ! SUBROUTINE indexx(arr,index)
  !
  ! PURPOSE
  ! Indexes an array arr, i.e., outputs the array index of length N such that
  ! arr(index(j)) is in ascending order for j = 1, 2, . . . ,N.
  ! The input quantity arr is not changed.
  !
  ! (Uses SORTTF from CernLib; cf. SORTZV)
  !
  ! INPUTS
  ! * real, dimension(:) :: arr
  !
  ! OUTPUT
  ! * integer, dimension(:) :: index
  !
  ! NOTES
  ! the integer array "index" has to be provided as input. The size may not
  ! be smaller than the size of "arr"
  !
  ! USAGE
  ! Assume you have to arrays A and B and you want to access the information
  ! first in an unsorted way and secondly sorted according the values of
  ! array A.
  !
  ! Example:
  !
  !    do i=1,n
  !      write(*,*) A(i),B(i)         ! unsorted
  !    enddo
  !    call indexx(A,ii)
  !    do i=1,n
  !      write(*,*) A(ii(i)),B(ii(i)) ! sorted according A
  !    enddo
  !
  !****************************************************************************
  SUBROUTINE indexx(arr,index)

    REAL, DIMENSION(:), INTENT(IN) :: arr
    INTEGER, DIMENSION(:), INTENT(OUT) :: index

    INTEGER :: n,i

    n = size(index)
    if (n .ne. size(arr)) then
       write(*,*) 'indexx: sizes must be equal:',n,size(arr)
       stop
    end if

    do i=1,n
       index(i) = i
    end do

    call SORTTF(arr,index,n)

  contains

    !==========================================================================
    SUBROUTINE SORTTF (A,INDEX,N1)

!      DIMENSION A(N1),INDEX(N1)
      real A(*)
      integer index(*),N1

      real AI
      integer N, I1,I2,I22,I222,I3,I33


      N = N1
      do 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      if (I2) 3,3,2
    2 I22 = INDEX(I2)
      if (AI.LE.A (I22)) GO TO 3
      INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      if (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      if (I2.LE.N) I22= INDEX(I2)
      if (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      if (A(I22)-A(I222)) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 if (AI-A(I22)) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
   END SUBROUTINE SORTTF
   !===========================================================================


  END SUBROUTINE indexx

!******************************************************************************

end module sorting
