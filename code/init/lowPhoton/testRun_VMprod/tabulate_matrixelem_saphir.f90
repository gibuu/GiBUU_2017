program tab_me_saphir

use IdTable, only: nucleon,omegaMeson
use particleProperties, only: initParticleProperties
use twoBodyPhaseSpace, only: Integrate_2bodyPS_resonance
use spline, only: bsplint2
use inputGeneral, only: readInputGeneral,path_To_Input
use constants, only: mN

implicit none

real, dimension(1:23) :: x,y
real, dimension(0:100) :: M,s
real, parameter :: srts_min=1.72,srts_max=2.45
real :: srts,k,E_g,sig,ps(1:5)
integer :: i

call readInputGeneral
call initParticleProperties

! read cross section data from file
open(202,file=trim(path_To_Input)//"/gammaN_omegaN_saphir.dat")
Do i=lBound(y,dim=1),uBound(y,dim=1)
  read(202,*) x(i),y(i)
end do
close(202)

! calculate matrix element
do i=lBound(M,dim=1),uBound(M,dim=1)
  srts = srts_min+i*(srts_max-srts_min)/(uBound(M,dim=1)-lBound(M,dim=1))
  s(i)=srts
  k = (srts**2-mN**2)/(2*srts)
  ! get exp. cross section
  E_g = (srts**2-mN**2)/(2*mN)
  sig = bsplint2(x,y,E_g)
  ! phase space
  ps = Integrate_2bodyPS_resonance(omegaMeson,srts,mN,0.)
  ! matrix element
  M(i)=k*srts**2*sig/ps(1)
end do

! write matrix element table to file
open(203,file=trim(path_To_Input)//"/gammaN_omegaN_ME_saphir.dat")
do i=lBound(M,dim=1),uBound(M,dim=1)
  write(203,*) s(i),M(i)
end do
close(203)


end program 
