program test_nBodyPhaseSpace

use IdTable
use particleProperties
use nBodyPhaseSpace

implicit none

integer :: i
real :: srts
real, dimension(1:6) :: mass
real, dimension(1:5) :: phaseSpace

  call initParticleProperties

      open(1,file='test_nBodyPhaseSpace.dat',status='unknown')
      write(1,*)'# Phase space volumes for different outgoing channels'
      write(1,*)'# srts:   NN:    NNpi:   NN2pi:  NN3pi:  NN4pi:'

      do i = 1,300

        srts = 2.*baryon(nucleon)%mass + 0.01*i

        mass(1) = baryon(nucleon)%mass
        mass(2) = baryon(nucleon)%mass
        mass(3) = meson(pion)%mass
        mass(4) = meson(pion)%mass
        mass(5) = meson(pion)%mass
        mass(6) = meson(pion)%mass

        call integrate_nBodyPS(srts,mass(1:2),phaseSpace(1))
        call integrate_nBodyPS(srts,mass(1:3),phaseSpace(2))
        call integrate_nBodyPS(srts,mass(1:4),phaseSpace(3))
        call integrate_nBodyPS(srts,mass(1:5),phaseSpace(4))
        call integrate_nBodyPS(srts,mass(1:6),phaseSpace(5))

        write(1,'(1x,6(e13.6,1x))') srts, phaseSpace

      enddo

end program test_nBodyPhaseSpace


       
        
