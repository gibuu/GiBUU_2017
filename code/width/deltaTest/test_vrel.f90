program test
use inputGeneral
use particleProperties, only: initParticleProperties

real, dimension(1:3) :: p
real :: mass=1.2
real :: rho=0.17

call initParticleProperties
call readinputGeneral



p=(/0.,0.,0./)
call tester

p=(/1.,0.,0./)
call tester

p=(/1.,2.,0./)
call tester

p=(/1.,2.,4./)
call tester

p=(/111.,2.,4./)
call tester

p=(/1111.,2.,121212124./)
call tester


contains 
  subroutine tester
    use relativeVelocity
    implicit none
    real :: quad, riemann
    quad=vrel_times_rho(rho,mass,p,.true.)
    riemann=vrel_times_rho(rho,mass,p,.false.)
    write(*,*) quad/(rho*0.197**3), riemann/(rho*0.197**3)
  end subroutine tester
  
end program test
