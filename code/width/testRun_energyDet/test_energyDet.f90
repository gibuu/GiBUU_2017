program test
  use idTable
  use inputGeneral, only: readinputGeneral
  use particleProperties, only: initParticleProperties
  implicit none

  write(*,*) "**********************************" 
  write(*,*) "**********************************" 
  write(*,*)
  Write(*,*)      "Testing energy determination"
  write(*,*)
  write(*,*) "**********************************" 
  write(*,*) "**********************************" 
  write(*,*)
  Write(*,*)      "Calling'init_Database'"
  call initParticleProperties
  call readinputGeneral

  call tester


end program test

subroutine tester
  use particleDefinition
  use energycalc, only : energyDetermination
  use IdTable, only: rho
  implicit none
  type(particle) :: p
  real, dimension(1:3) :: beta
  real :: error
  p%ID=rho
  p%charge=0
  p%mass= 0.6375E+00
  p%momentum=(/ 0.8660E+00 , 0.7168E-01,  0.9806E-01,  0.5772E+00 /)
  p%perturbative=.true.
  p%offshellParameter= -0.4922E+00
  beta=(/  0.1453E+00, -0.2453E-01, -0.4565E+00 /)
  p%position=(/-0.227E+01 ,    -0.187E+01  ,    0.181E+01/)
  p%velocity=(/ -0.162E+00  ,   -0.563E-01  ,    0.887E+00 /)
  call energyDetermination(p,beta,100,error)
  write(*,*) 'energy=',p%momentum(0)
  write(*,*) 'error= ',error
end subroutine tester
