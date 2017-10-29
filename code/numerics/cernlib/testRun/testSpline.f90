program test

use cl_splines, only: tspline, cl_initSpline, cl_spline
! use spline

implicit none
integer, parameter :: N=6
integer, parameter :: mode=0
real, dimension(0:N) :: X,Y  !,Y2

type(tspline) :: s,t
real :: xs,ys,yt, ys_old
integer :: error
logical :: flag
integer :: i

x=(/-4.,-1.,2.,3.,7.,8.,10./)

y=(/1.,2.,-20.,6.,4.,31.,1./)

do i=0,n
   write(30,*) x(i), y(i)
end do

! call spline2(x,y,y2)

s=cl_initSpline(x,y)
t=cl_initSpline(x,y,0)


do i=0,1000
   xs=-30+float(i)*0.1
   ys=cl_spline(s,xs,flag,error)
   yt=cl_spline(t,xs,flag,error)
!    call splint2(x,y,y2,xs,ys_old)
   write(*,'(4E20.4,L,I8)') xs,ys,yt,ys_old,flag,error
   write(20,'(4E20.4,L,I8)') xs,ys,yt,ys_old,flag,error
end do


end program test
