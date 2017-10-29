! Test program to compare the different form factor paramerizations.
!
! (to compile this, the "parametrization"-flag in FF_QE_nucleonScattering must temporarily be declared "public")

program testFF_para
use FF_QE_nucleonScattering, only: formfactors_QE, parametrization
implicit none
real, parameter :: dq=0.005
integer :: i,j
real :: F1,F2, Qsquared,dummy,GE(0:3),GM(0:3)
call init_Database

open(10,file='proton_GE.dat')
open(11,file='proton_GM.dat')
open(12,file='neutron_GE.dat')
open(13,file='neutron_GM.dat')

do i=0,400
   QSquared=dq*i

   do j=0,3
     parametrization=j
     call formfactors_QE(QSquared,2,1,F1,F2,dummy,dummy,GE(j),GM(j))
   end do

   write(10,'(5F12.5)') QSquared, GE
   write(11,'(5F12.5)') QSquared, GM

   do j=0,3
     parametrization=j
     call formfactors_QE(QSquared,2,0,F1,F2,dummy,dummy,GE(j),GM(j))
   end do

   write(12,'(5F12.5)') QSquared, GE
   write(13,'(5F12.5)') QSquared, GM
end do 


end program testFF_para
