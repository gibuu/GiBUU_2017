program test_rho_dil
  use inputGeneral, only: readinputGeneral
  use particleProperties, only: hadron, InitParticleProperties
  use mesonWidth, only: fullWidthMeson
  use mesonWidthVacuum, only: dileptonWidth
  use IdTable, only: rho
  use constants, only: pi

  implicit none

  integer, parameter :: ID = rho
  real, parameter :: dm = 0.001
  
  integer :: i
  real :: mass, M0, width_had, width_dil, width_tot, sf1, sf2, dil1, dil2

  call readinputGeneral
  call InitParticleProperties

  M0 = hadron(ID)%mass
  Write(*,*)      "Testing Vacuum Decay Width of Meson #",ID
  do i=0,2000
     mass=i*dm
     width_had = fullWidthMeson(ID,mass)
     width_dil = dileptonWidth(ID,mass)
     width_tot = width_had+width_dil
     sf1 = width_had/pi*mass/((mass**2-M0**2)**2 + mass**2*width_had**2)
     sf2 = width_tot/pi*mass/((mass**2-M0**2)**2 + mass**2*width_tot**2)
     dil1 = sf1 * width_dil/width_had
     dil2 = sf2 * width_dil/width_tot
     write(33,'(8F15.9)') mass,width_had,width_dil,width_tot,sf1,sf2,dil1,dil2
  end do

end
