program testSpecFuncMes

use particleDefinition
use inputGeneral, only: readInputGeneral
use particleProperties, only: initParticleProperties
use spectralFuncMesons
use selfEnergy_Mesons, only: getDispersion

implicit none

integer,parameter :: N = 3000
real, parameter :: dm = 0.001

integer :: i, cien
real :: sf
type(particle) :: p

call readInputGeneral
call initParticleProperties

if (getDispersion()) then
  cien = 100
else
  cien = 0
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

p%position = (/ 999., 999., 999. /) ! vacuum

p%ID = 103
do i=0,N
  p%momentum(0) = i*dm
  sf = specFuncMes (p)
  write(p%ID+cien,*) p%momentum(0),sf
end do

p%ID = 105
do i=0,N
  p%momentum(0) = i*dm
  sf = specFuncMes (p)
  write(p%ID+cien,*) p%momentum(0),sf
end do

p%ID = 107
do i=0,N
  p%momentum(0) = i*dm
  sf = specFuncMes (p)
  write(p%ID+cien,*) p%momentum(0),sf
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

p%position = (/ 0., 0., 0. /) ! medium

p%ID = 103
do i=0,N
  p%momentum(0) = i*dm
  sf = specFuncMes (p)
  write(1000+p%ID+cien,*) p%momentum(0),sf
end do

p%ID = 105
do i=0,N
  p%momentum(0) = i*dm
  sf = specFuncMes (p)
  write(1000+p%ID+cien,*) p%momentum(0),sf
end do

p%ID = 107
do i=0,N
  p%momentum(0) = i*dm
  sf = specFuncMes (p)
  write(1000+p%ID+cien,*) p%momentum(0),sf
end do

end program
