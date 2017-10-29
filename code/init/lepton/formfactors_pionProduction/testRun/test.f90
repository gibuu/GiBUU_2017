program test
use formfactors_A_main
use particleProperties, only: initParticleProperties
use inputGeneral, only: readinputGeneral
implicit none
real :: QS

call readinputGeneral
call initParticleProperties
QS=0.
do 
QS=Qs+0.01
write(10,'(13E15.5)') QS,getA(-1,1,90.,1.8,QS)
end do
end program
