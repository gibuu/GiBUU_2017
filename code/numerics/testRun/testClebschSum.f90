program TestClebschSum

use ClebschGordan, only: CG

implicit none

real :: j1,j2,j,m

real, parameter :: jmax = 50.

j1=0.
do while (j1<=jmax)
  j2=0.
  do while (j2<=jmax)
     j=abs(j1-j2)
     do while (j<=j1+j2)
       m=-j
       do while (m<=j)
         call checkSum(j1,j2,j,m)
         m=m+1
       end do
       j=j+1.
     end do
    j2=j2+0.5
  end do
  j1=j1+0.5
end do

contains

  subroutine checkSum(j1,j2,j,m)
    implicit none
    real,intent(in) :: j1,j2,j,m
    real :: m1,s,c
    s=0.
    m1 = max(-j1,m-j2)
    do while (m1<=min(j1,m+j2))
      c = CG(nint(2*j1),nint(2*j2),nint(2*j),nint(2*m1),nint(2*(m-m1)))**2
      s = s + c
      m1 = m1 + 1.
    end do
    if (abs(s-1.)>1E-6) then
      print '(A,4G12.5)',"sum rule violated:",j1,j2,j,m
      print *,"sum = ",s
      stop
    else
      print '(A,4G12.5)',"sum rule fulfilled:",j1,j2,j,m
    end if
  end subroutine

end program 
