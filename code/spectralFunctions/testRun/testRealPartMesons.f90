program testRealPartMesons

use constants, only: rhoNull
use inputGeneral
use particleProperties, only: initParticleProperties, hadron
use selfEnergy_mesons
use derivatives, only: finiteDifference
use mediumDefinition

implicit none

integer,parameter :: ID = 105
real, parameter :: dm = 0.001
integer, parameter :: N = 3000

real :: m,im(0:N),re(0:N),m0,y,c
integer :: i

type(medium) :: med

call initParticleProperties
call readInputGeneral

!call tabulate_realPart()

m0 = hadron(ID)%mass

! (1) vacuum
med%density = 0.
med%densityProton = 0.
med%densityNeutron = 0.
med%useMedium = .false.

do i=0,N
  m = i*dm
  im(i) = get_imagPart (ID, m, med)
  re(i) = get_realPart (ID, m, med)
end do

do i=2,N-2
  m = i*dm
  y = ( m**2 - m0**2 ) / (im(i))
  c = ( y*finiteDifference(im(i-2:i+2),dm,2,0) + finiteDifference(re(i-2:i+2),dm,2,0) ) / (2.*m)
  write (ID,'(5G13.5)') m,im(i),re(i),y,c
end do

! (2) medium
med%density = rhoNull
med%densityProton = rhoNull/2.
med%densityNeutron = rhoNull/2.
med%useMedium = .true.

do i=0,N
  m = i*dm
  im(i) = get_imagPart (ID, m, med)
  re(i) = get_realPart (ID, m, med)
end do

do i=2,N-2
  m = i*dm
  y = ( m**2 - m0**2 ) / (im(i))
  c = ( y*finiteDifference(im(i-2:i+2),dm,2,0) + finiteDifference(re(i-2:i+2),dm,2,0) ) / (2.*m)
  write (ID+100,'(5G13.5)') m,im(i),re(i),y,c
end do

end program 
