program test

implicit none

! call splineTest
call clebschTest

contains


  subroutine splineTest
    use spline
    real, dimension(1:5) :: y, y2
    real, dimension(1:5) :: x
    real :: value, res
    integer :: i
    x=(/ 1.,2.,12.,30.,50. /)
    y=(/ 1.,3.,126.,200.,-50.1 /)
    call spline2(x,y,y2)
    Do i=lBound(y,dim=1),uBound(y,dim=1)
       write(201,*) x(i), y(i)
    End do
    Do i=-100,100
       value=0.6*i
       call splint2(x,y,y2,value,res)
       write(200,*) value, res
    End do
  end subroutine splineTest


  subroutine clebschTest
    use clebschGordan, only: clebschSquared
    use random, only: rn

    real :: i1,i2,i3,iz1,iz2,iz1_a,iz2_a
    integer :: i

    do i=1,100000
    i1=NINT(rn()*5)/2.
    i2=NINT(rn()*5.)/2.
    i3=NINT(rn()*5.)/2.
    do iz1=-Abs(2*i1),Abs(2*i1),2  !loop over all possible z-Values of the isospin
      do iz2=-Abs(2*i2),Abs(2*i2),2
         iz1_a=iz1/2.
         iz2_a=iz2/2.
         Write(*,*) clebschSquared(i1,i2,i3,iz1_a,iz2_a)
      end do
    End do
    End do
  End subroutine clebschTest

end program test
