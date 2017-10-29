program test
  
  use inputGeneral
  implicit none


  call readinputGeneral
  call init_database


call sub



end program


subroutine sub
use degRad_conversion
use photon_flux
implicit none

real :: ebeam=0.855
real :: efinal=0.2
real :: theta=30 
real, dimension(0:3) :: q,pi,pf

pi=(/ebeam,0.,0.,ebeam/)
pf=(/efinal,0.,efinal*sin(radian(theta)),efinal*cos(radian(theta))/)

q=pi-pf

write(*,*) gamma(Ebeam,q,theta)

end subroutine sub
