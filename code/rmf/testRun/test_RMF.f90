program test_RMF

use inputGeneral, only: readInputGeneral
use particleProperties, only: initParticleProperties
use RMF, only: getRMF_parSet, walecka, m_nucleon
use constants, only : pi, hbarc

implicit none

real :: mass_seed, rhobar, shift, rhos, energyDensity, press, U, p_fermi
real :: ScalarPotential, VectorPotential
integer :: i
character(len=50) :: f

call readInputGeneral
call initParticleProperties

write(f,'(A,i1,A)') 'RMF_set',getRMF_parSet(),'.dat'

open(1,file=trim(f),status='unknown')
write(1,*)'# rhobar, fm^-3:    p_fermi, GeV:   rhos, fm^-3:   m^*, GeV:'
write(1,*)'# (cont.)  E/A-m_N, GeV:  p, GeV^-3:  S, GeV:  V, GeV:   U, GeV:'

mass_seed = m_nucleon

do i=0,400

   rhobar = i*0.001

   p_fermi = hbarc*(1.5*pi**2*rhobar)**0.333333

   call walecka (rhobar, shift, em0=mass_seed, rhoscalar=rhos,     &
                 endens=energyDensity, pressure=press,             &
                 S=ScalarPotential, V=VectorPotential, potential=U)

   write(1,'(9(e13.6,1x))') rhobar, p_fermi, rhos, m_nucleon-shift,&
                            energyDensity/rhobar-m_nucleon, press,&
                            ScalarPotential, VectorPotential, U

   mass_seed = m_nucleon - shift

end do

close(1)

end program test_RMF
