program testFF
use particleProperties, only: initParticleProperties
use FF_QE_nucleonScattering
implicit none
real, parameter :: dq=0.005
integer :: i
real :: F1,F2, Qsquared,dummy,GE,GM
call initParticleProperties

open(10,file='protonFF.dat')
open(11,file='neutronFF.dat')
do i=0,400
   QSquared=dq*i
   call formfactors_QE(QSquared,2,1,F1,F2,dummy,dummy,GE,GM)
   write(10,'(5F12.5)') QSquared, F1, F2,GE,GM
   call formfactors_QE(QSquared,2,0,F1,F2,dummy,dummy,GE,GM)
   write(11,'(5F12.5)') QSquared, F1, F2,GE,GM
end do 


end program testFF
