program main
use IDtable
use particleDefinition
use particleProperties
use energycalc
use master_1body
implicit none
type(particle),dimension(1:1) :: teilchen
type(particle),dimension(1:20) :: finalstate
real, dimension(1:3) :: betaToLRF
logical :: collisionflag,pauliflag

real :: time

call init_Database

betaToLRF=0.
teilchen(1)%ID=Delta
teilchen(1)%momentum(1:3)=(/12.,3.,4./)
teilchen(1)%mass=1.212
teilchen(1)%momentum(0)=freeEnergy(teilchen(1))
teilchen(1)%perturbative=.true.
betaToLrf(1)=0.8
!Print * , teilchen%momentum
!call energyDetermination(teilchen, betaToLRF)
!Print * , teilchen%momentum

time=20.
call decayParticle(teilchen,finalState,collisionFlag,time,pauliFlag,.true.)
Write(*,*) SQRTS(finalState(1),finalState(2))

end program main


