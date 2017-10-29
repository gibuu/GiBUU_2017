program main

  use inputGeneral
  implicit none 

  call init_database
  call readInputGeneral

  call testMomConserv

end program main



subroutine testMomConserv

  use collisionTerm
  use particleDefinition
  use particleProperties
  use IdTable
  use inputGeneral
  use timing
  implicit none

  type(particle), dimension (:,:),Allocatable :: teilchenPert
  type(particle), dimension (:,:),Allocatable :: teilchenReal
  type(particle), dimension (1:2) :: pair
  integer :: i,j,chargeVorher,chargeNachher
  integer, parameter :: m=10000 ! number ensembles
  real,dimension(0:3) :: momTot_vorher,momTot_nachher

  write(*,*)
  write( *,*) '******************************************************'
  write( *,*) 'Testing collisionTerm'
  write( *,*) '******************************************************'


  momTot_vorher=0.
  momTot_nachher=0.
  chargeNachher=0
  chargeVorher=0


  Allocate(teilchenPert(1:m,1:1))
  Allocate(teilchenReal(1:m,1:100))

  pair%Id=(/nucleon,nucleon/)
  pair%charge=(/1,1/)

  pair(1)%event=(/1,1/)
  pair(2)%event=(/2,2/)

  pair%mass=(/baryon(nucleon)%mass,baryon(nucleon)%mass/)

  pair(1)%position=(/0.1,0.,0./)
  pair(2)%position=(/0.,0.,0./)

  pair(1)%momentum(1:3)=(/0.,-0.3,0./)
  pair(2)%momentum(1:3)=(/0.,27.3,0./)

  pair(1)%momentum(0)=sqrt(pair(1)%mass**2+Dot_Product(pair(1)%momentum(1:3),pair(1)%momentum(1:3)))
  pair(2)%momentum(0)=sqrt(pair(2)%mass**2+Dot_Product(pair(2)%momentum(1:3),pair(2)%momentum(1:3)))

  pair(1)%velocity=pair(1)%momentum(1:3)/pair(1)%momentum(0)
  pair(2)%velocity=pair(2)%momentum(1:3)/pair(2)%momentum(0)

  Print*, 'Velocities=', pair(1)%velocity,pair(2)%velocity

  Do i=1,m
     teilchenPert(i,:)%ID=0
     teilchenReal(i,1)=pair(1)
     teilchenReal(i,2)=pair(2)
     teilchenReal(i,1)%perturbative=.false.
  End do

  call setNumber(teilchenReal)
  call setNumber(teilchenPert)

  Print *, 'positions'
  Print *, pair(1)%position
  Print *, pair(2)%position
  Print *, 'momenta'
  Print *, pair(1)%momentum
  Print *, pair(2)%momentum
  Print *, 'Wurzel(s)=',sqrts(pair(1),pair(2))
  Print *, 'Total momentum=',pair(1)%momentum+pair(2)%momentum
  print * ,'***********************'

  Print *, 'Testing the collisions'
  !Do i =1,10
  call timeMeasurement(.true.)

  Print *, 'Before the collisions, Printing real vector'
  Print *, '*************'
  Print *, 'Id,Charge,Perturbative, antiparticle, perweight,scaleCS, collTime, prodTime, inFormation, number,event'
  do i=lBound(teilchenReal,dim=1),uBound(teilchenReal,dim=1)
     do j=lBound(teilchenReal,dim=2),uBound(teilchenReal,dim=2)
        call writeParticleOut(teilchenReal(i,j),1)
        if (teilchenReal(i,j)%ID > 0) then
           momTot_vorher=momTot_vorher+teilchenReal(i,j)%momentum
           chargeVorher=chargeVorher+teilchenReal(i,j)%charge
        end if
     end do
  end do
  Print *, ' '
  Print *, ' '
  Print *, ' '
  Print *, ' '
  call collideMain(teilchenPert,teilchenReal,0.3)
  call timeMeasurement(.false.)


  Print *, 'Finished the collisions, Printing real vector'
  Print *, '*************'
  Print *, 'Id,Charge,Perturbative, antiparticle, perweight,scaleCS, collTime, prodTime, inFormation, number,event'
  do i=lBound(teilchenReal,dim=1),uBound(teilchenReal,dim=1)
     do j=lBound(teilchenReal,dim=2),uBound(teilchenReal,dim=2)
        If(teilchenReal(i,j)%ID > 0) then
           call writeParticleOut(teilchenReal(i,j),2)
           momTot_nachher=momTot_nachher+teilchenReal(i,j)%momentum
           chargeNachher=chargeNachher+teilchenReal(i,j)%charge
        end if
     end do
  end do


  Print *, 'vorher', momTot_vorher,chargeVorher
  Print *, 'nachher', momTot_nachher,chargeNachher



end subroutine testMomConserv


subroutine writeParticleOut(teilchen,j)
  use particleDefinition
  type(particle) :: teilchen
  integer :: j
  write(10+j,'(2I3,2L,4F6.3,L,3I8)')&
       & teilchen%ID,teilchen%Charge,teilchen%perturbative,teilchen%antiparticle, teilchen%perweight,teilchen%scaleCS, &
       & teilchen%lastCollisionTime, teilchen%productionTime, &
       & teilchen%in_Formation,teilchen%number,teilchen%event 
end subroutine writeParticleOut

