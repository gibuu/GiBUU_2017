program testGamma

use minkowski, only: gamma,gamma_old

complex, dimension(0:3,0:3) :: g

integer :: i,j,k

do k=0,11
  if (k==4) cycle
  g = gamma_old(k)
  do i=0,3
  do j=0,3
    write(*,"(3I3,6F12.6)") i,j,k,gamma(i,j,k),g(i,j),gamma(i,j,k)-g(i,j)
    if (gamma(i,j,k)-g(i,j)/=0.) write(*,*) "ERROR!"
  end do
  end do
  write(*,*)
end do


end program 
