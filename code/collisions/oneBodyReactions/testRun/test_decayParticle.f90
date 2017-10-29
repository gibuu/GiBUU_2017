program test_decayParticle
use master_1body
use IdTable
use particleProperties, only: initParticleProperties
use particleDefinition
use minkowski
implicit none

type(particle) :: resonance
type(particle),dimension(1:20) :: finalState
logical :: pauliFlag, collisionFlag
real :: time
integer :: i,j
real, dimension(0:3) :: momTot

call initParticleProperties


write(*,*) '**************************************************************************'

Write(*,*) 'Testing delta Bar'

time=0.
resonance%ID=delta
resonance%antiParticle=.true.
resonance%Charge=-2

resonance%Mass=   1.19628101249960    
resonance%Momentum=  (/ 0.125499491E+01   ,      -0.413693032E+00, -0.150027894E+00,-0.516731763E-01/)
resonance%Position=  (/ 0.189334479E+01   ,      -0.345616075E+01, -0.147158862E+01 /)
resonance%Velocity=  (/ -0.392547837E+00  ,       -0.142500801E+00 ,        -0.490951234E-01 /)
resonance%Offshellparameter=  0.000000000000000E+000 





Do j=1,100
   finalState(1:4)%ID=0
   write(*,*)
   write(*,*) '***Run:' , j
   call decayParticle(resonance,finalState,collisionFlag,pauliFlag,.false.,time)
   !Write(*,*) 'Flags=' , CollisionFlag, pauliFlag

   If(collisionFlag) then
      Write(*,*) finalState(1:4)%ID
      Write(*,*) finalState(1:4)%antiParticle
      Write(*,*) finalState(1:4)%Charge
      momtot=0.
      Do i=1,4
         if(finalState(i)%ID.gt.0)  momtot=momtot+finalState(i)%Momentum
      end do
      !  exit
      write(*,'(A4,5F15.5)') 'In' ,resonance%momentum, abs4(resonance%momentum)
      write(*,'(A4,5F15.5)') 'Out',momtot, abs4(momTot)
   end if
End do


end program test_decayParticle
