program main
use master_1body
use IdTable
use particleProperties
use particleDefinition

implicit none

type(particle),dimension(1:1) :: resonance
type(particle),dimension(1:20) :: finalState
logical :: pauliFlag, collisionFlag
real :: time, mom
integer :: i,j
integer :: counter(1:3)
call init_dataBase

Write(*,*) 'etaPrime'

time=0.
resonance%ID=etaPrime
resonance%antiParticle=.false.
resonance%Charge=0
resonance%mass=meson(etaPrime)%mass
mom=0.5

resonance%momentum(0)=Sqrt(resonance%mass**2+mom**2)
resonance%momentum(1)=mom
counter=0.

Do j=1,10000
   write(*,*)
   write(*,*) '***Run:' , j
   call decayParticle(resonance,finalState,collisionFlag,time,pauliFlag,.true.)
   Write(*,*) 'Flags=' , CollisionFlag, pauliFlag
   If(collisionFlag) then
      Write(*,*) finalState(1:4)%ID
      Write(*,*) finalState(1:4)%antiParticle
      Write(*,*) finalState(1:4)%Charge
      Do i=1,4
         Write(*,*) finalState(i)%Momentum
      end do
      if((finalState(1)%ID.eq.pion).and.(finalState(2)%ID.eq.pion).and.(finalState(3)%ID.eq.eta)) then
         if((finalState(1)%charge.eq.0).and.(finalState(2)%charge.eq.0).and.(finalState(3)%charge.eq.0)) then
            counter(1)=counter(1)+1
         else
            counter(2)=counter(2)+1
         end if
         !   exit
      else if(finalState(3)%ID.eq.0) then
         counter(3)=counter(3)+1
      end if
   end if
End do

write(*,*) counter/10000.





end program main
