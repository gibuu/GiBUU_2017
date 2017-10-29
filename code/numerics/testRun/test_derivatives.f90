

!!$ Testing program for the derivatives module
program test
use derivatives
implicit none
real :: dx=0.1
real :: x
real, dimension(-2:2) :: y
integer :: scheme, order

do scheme=-1,1
   do order=1,2
      x=0.
      do 
         x=x+dx
         y=(/sin(x-2.*dx),sin(x-dx),sin(x),sin(x+dx),sin(x+2.*dx)/)
         write(order*100+scheme+1,'(4G20.7)') x, finiteDifference(y,dx,order,scheme), cos(x),derivative(sin,x,dx,order,scheme)
         if (x.gt.10.) exit
      end do
   end do
end do
end program test
