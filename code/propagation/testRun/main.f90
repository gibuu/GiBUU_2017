program main

  use particleDefinition
  use idTable, only: nucleon
  use particleProperties, only: initParticleProperties
  use propagation, only: propagate_euler
  use inputGeneral, only: readinputgeneral, delta_T
  implicit none

  type(particle),dimension(1:1,1:1) :: teilchenPert

  type(particle) :: teilchen

  call initParticleProperties
  call readinputgeneral

 teilchen%Mass=   1.00696913621373  
 teilchen%Momentum=   (/1.01105783827159,   8.864022276520209E-002 , -0.308415047504890   ,    1.386015716455711E-002/)
 teilchen%OFfshellparameter=   6.84038642412723
 teilchen%Position= (/-0.666866762400002   ,    0.619665061175609 ,  0.165761574601560 /)
 teilchen%velocity= (/-3.744196450168186E-003 ,-0.187235768511816, -8.805597163916382E-003 /)
 teilchen%Id=nucleon
 
 teilchenPert=teilchen

 write(*,*) teilchenPert(1,1)%momentum
 
 call Propagate_euler(teilchenPert,delta_T)
 
 If(Dot_Product(teilchen%velocity,teilchen%velocity).ge.1) then
    write(*,*) ' Error : Particles faster than 1:',Dot_Product(teilchen%velocity,teilchen%velocity)
 end if

end program main


