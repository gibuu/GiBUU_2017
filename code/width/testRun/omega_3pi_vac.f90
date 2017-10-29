! This test case calculates the cross section for e+e- -> omega -> 3pi using the width parametrization from GiBUU,
! in order to compare to the data shown in Klingl NPA650 (1999).

program test_omega_3pi
  use inputGeneral
  use particleProperties
  use mesonWidth
  use IdTable, only: pion, omegaMeson, photon
  use constants, only: melec
  use twoBodyTools, only: pCM
  implicit none

  integer :: i, ID
  real, parameter :: dm = 0.001
  real :: mass, M0, width_tot, width_3pi, width_ee, width_p0g, p_ab, sig, BR

  call InitParticleProperties
  call readinputGeneral

  ID = omegaMeson
  M0 = meson(ID)%mass
  Write(*,*)      "Testing Vacuum Decay Width of Meson #",ID
  do i=0,5000
     mass=i*dm
     width_tot = fullWidthMeson(ID,mass)                                ! omega -> X
     width_3pi = partialwidthMeson(ID,mass,pion,pion,pion,+1,-1,0)      ! omega -> 3pi
     width_ee  = 0.767e-06 * M0**4 / mass**3                            ! omega -> e+e-
     p_ab = pCM(mass,melec,melec)
     ! cross section for e+e- -> omega -> 3pi, cf. Effenberger (2.52), leaving out constant factors
     sig = mass**2 * width_ee * width_3pi / p_ab**2 / ((mass**2-M0**2)**2+mass**2*width_tot**2)
     write(ID+100,'(6F15.9)') mass, width_tot, width_3pi, width_ee, p_ab, sig
  end do 

  Write(*,*)      "Testing Branching Ratio omega -> pi0 gamma"
  do i=0,5000
     mass=i*dm
     width_tot = fullWidthMeson(ID,mass)                                ! omega -> X
     width_p0g = partialwidthMeson(ID,mass,pion,photon)      ! omega -> pi0 gamma
     if (width_tot>0.) then
       BR = width_p0g/width_tot
     else
       BR = 0.
     end if
     write(ID+200,'(4F15.9)') mass, width_tot, width_p0g, BR
  end do 

end
