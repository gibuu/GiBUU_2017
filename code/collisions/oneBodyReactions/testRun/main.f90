program main
use master_1body
use IdTable
use particleProperties, only: hadron, initParticleProperties
use particleDefinition

implicit none

type(particle) :: resonance
type(particle),dimension(1:20) :: finalState
logical :: pauliFlag, collisionFlag
real :: time, mom
integer :: i,j

call initParticleProperties

Write(*,*) 'Testing dStar'

time=0.
resonance%ID=dStar
resonance%antiParticle=.true.
resonance%Charge=-1
resonance%mass=hadron(dSStar_plus)%mass+0.1
mom=0.5

resonance%momentum(0)=Sqrt(resonance%mass**2+mom**2)
resonance%momentum(1)=mom


Do j=1,100
write(*,*)
write(*,*) '***Run:' , j
call decayParticle(resonance,finalState,collisionFlag,pauliFlag,.false.,time)
!Write(*,*) 'Flags=' , CollisionFlag, pauliFlag

If(collisionFlag) then
   Write(*,*) finalState(1:4)%ID
   Write(*,*) finalState(1:4)%antiParticle
   Write(*,*) finalState(1:4)%Charge
   Do i=1,4
      Write(*,*) finalState(i)%Momentum
   end do
!   exit
end if
End do



write(*,*) '**************************************************************************'

Write(*,*) 'Testing delta Bar'

time=0.
resonance%ID=delta
resonance%antiParticle=.true.
resonance%Charge=-2
resonance%mass=hadron(delta)%mass+0.1
mom=0.5

resonance%momentum(0)=Sqrt(resonance%mass**2+mom**2)
resonance%momentum(1)=mom


Do j=1,100
write(*,*)
write(*,*) '***Run:' , j
call decayParticle(resonance,finalState,collisionFlag,pauliFlag,.false.,time)
!Write(*,*) 'Flags=' , CollisionFlag, pauliFlag

If(collisionFlag) then
   Write(*,*) finalState(1:4)%ID
   Write(*,*) finalState(1:4)%antiParticle
   Write(*,*) finalState(1:4)%Charge
   Do i=1,4
      Write(*,*) finalState(i)%Momentum
   end do
 !  exit
end if
End do


end program main
