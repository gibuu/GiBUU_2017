program test
call tester
end program test


subroutine tester()
use interpolation
implicit none

real, dimension(-1:0,1:2) :: table
integer :: j,i
real, dimension (1:2) :: X, A1,A2
table(-1,1)=1
table(0,2)=2
table(-1,2)=3
table(0,1)=4

A1=(/-1,0/)
A2=(/1,2/)


do i=-20,20
   x(1)=float(i)*0.2 
   do j=0,10
  x(2)=float(j)*0.2 
  
write(12,'(5F9.3)')  x(1), x(2),inter2(X,A1,A2,Table), table(min(0,max(-1,nint(x(1)))),min(2,max(1,nint(x(2)))))
end do
write(12,*) 
end do
end subroutine tester
