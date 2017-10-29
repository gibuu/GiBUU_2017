program main

  use inputGeneral
  implicit none 

  call init_database
  call readInput

  !  call testCollisionTerm
  !  call testXsectionMaster
  call test3Body

end program main


subroutine testXsectionMaster
  use collisionTerm
  use particleDefinition
  use particleProperties
  use IdTable
  use inputGeneral
  use timing
  use mediumDefinition
  use master_2Body
  use preEventDefinition

  implicit none
  logical :: HiEnergyFlag, plotFlag
  real :: srts, sigmaTot, sigmaElast,p
  real, dimension (0:3) :: momentumLRF
  integer :: i
  type(medium) :: mediumATcollision
  type(particle),dimension(1:2) :: teilchenIN
  type(preEvent),dimension(1:20) :: teilchenOUT

  write(*,*)
  write( *,*) '******************************************************'
  write( *,*) 'Testing Xsections'
  write( *,*) '******************************************************'

  Open(100,File='XSectionOutPut.dat')
  write(100,*) '# Collider:',teilchenIN(1)%ID,teilchenIN(1)%charge,'Target:',teilchenIN(2)%ID,teilchenIN(2)%charge
  write(100,*) '# p of Collider(Target at rest), sqrt(s),sigmaTot,sigmaElast'
  Do i=1,1000
     p=0.05*i

     teilchenIN%ID=(/pion,nucleon/)
     teilchenIN%mass=(/meson(teilchenIN(1)%ID)%mass,baryon(nucleon)%mass/)
     teilchenIN%Charge=(/-1,1/)
     mediumATCollision%UseMedium=.false.

     teilchenIN(1)%momentum(1:3)=(/p,0.,0./)
     teilchenIN(2)%momentum(1:3)=(/0.,0.,0./)
     teilchenIN(1)%momentum(0)=SQRT(teilchenIN(1)%mass**2+Dot_Product(teilchenIN(1)%momentum(1:3),teilchenIN(1)%momentum(1:3)))
     teilchenIN(2)%momentum(0)=SQRT(baryon(nucleon)%mass**2+Dot_Product(teilchenIN(1)%momentum(1:3),teilchenIN(1)%momentum(1:3)))

     momentumLRF=teilchenIn(1)%momentum+teilchenIN(2)%momentum
     srts=sqrtS(teilchenIN(1),teilchenIN(2))
     write(*,*) 'srts=', srts

     call XsectionMaster(srts,teilchenIn,mediumATcollision,momentumLRF,teilchenOut,sigmaTot,sigmaElast,HiEnergyFlag,plotFlag)
     Write(*,*) p, sigmaTot, sigmaElast, HiEnergyFlag
     Write(100,'(4F8.4,L)') p, srts, sigmaTot, sigmaElast, HiEnergyFlag
  end do

end subroutine testXsectionMaster







!******************************************



subroutine testCollisionTerm

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
  integer :: i,j
  integer, parameter :: m=3 ! number ensembles


  write(*,*)
  write( *,*) '******************************************************'
  write( *,*) 'Testing collisionTerm'
  write( *,*) '******************************************************'


  Allocate(teilchenPert(1:m,1:35))
  Allocate(teilchenReal(1:m,1:1))

  pair%Id=(/nucleon,nucleon/)
  pair%charge=(/1,1/)
  pair(1)%event=(/-999,-999/)
  pair(2)%event=(/-999,-999/)

  pair%mass=(/baryon(nucleon)%mass,baryon(nucleon)%mass/)

  pair(1)%position=(/0.1,0.,0./)
  pair(2)%position=(/0.,0.,0./)

  pair(1)%momentum(1:3)=(/0.,0.,0./)
  pair(2)%momentum(1:3)=(/0.,1.3,0./)

  pair(1)%momentum(0)=sqrt(pair(1)%mass**2+Dot_Product(pair(1)%momentum(1:3),pair(1)%momentum(1:3)))
  pair(2)%momentum(0)=sqrt(pair(2)%mass**2+Dot_Product(pair(2)%momentum(1:3),pair(2)%momentum(1:3)))

  pair(1)%velocity=pair(1)%momentum(1:3)/pair(1)%momentum(0)
  pair(2)%velocity=pair(2)%momentum(1:3)/pair(2)%momentum(0)
  Print*, 'Velocities=', pair(1)%velocity,pair(2)%velocity

  Do i=1,m
     teilchenReal(i,1)=pair(1)
     teilchenReal(i,1)%perturbative=.false.

     teilchenPert(i,(/1,5,10,15,20,25,30,35/))=pair(2)
     teilchenPert(i,:)%perturbative=.true.  
     teilchenPert(i,:)%perweight=float(i)
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

  Print *, 'Before the collisions, Printing perturbative vector'
  Print *, '*************'
  Print *, 'Id,Charge,Perturbative, antiparticle, perweight,scaleCS, collTime, prodTime, inFormation, number,event'
  do i=lBound(teilchenPert,dim=1),uBound(teilchenPert,dim=1)
     do j=lBound(teilchenPert,dim=2),uBound(teilchenPert,dim=2)
        call writeParticleOut(teilchenPert(i,j))
     end do
  end do

  Print *, 'Before the collisions, Printing real vector'
  Print *, '*************'
  Print *, 'Id,Charge,Perturbative, antiparticle, perweight,scaleCS, collTime, prodTime, inFormation, number,event'
  do i=lBound(teilchenReal,dim=1),uBound(teilchenReal,dim=1)
     do j=lBound(teilchenReal,dim=2),uBound(teilchenReal,dim=2)
        call writeParticleOut(teilchenReal(i,j))
     end do
  end do
  Print *, ' '
  Print *, ' '
  Print *, ' '
  Print *, ' '
  call collideMain(teilchenPert,teilchenReal,0.3)
  call timeMeasurement(.false.)
  !end do
  Print *, 'Finished the collisions, Printing perturbative vector'
  Print *, '*************'
  Print *, 'Id,Charge,Perturbative, antiparticle, perweight,scaleCS, collTime, prodTime, inFormation, number,event'
  do i=lBound(teilchenPert,dim=2),uBound(teilchenPert,dim=1)
     do j=lBound(teilchenPert,dim=2),uBound(teilchenPert,dim=2)
        call writeParticleOut(teilchenPert(i,j))
     end do
  end do

  Print *, 'Finished the collisions, Printing real vector'
  Print *, '*************'
  Print *, 'Id,Charge,Perturbative, antiparticle, perweight,scaleCS, collTime, prodTime, inFormation, number,event'
  do i=lBound(teilchenReal,dim=1),uBound(teilchenReal,dim=1)
     do j=lBound(teilchenReal,dim=2),uBound(teilchenReal,dim=2)
        call writeParticleOut(teilchenReal(i,j))
     end do
  end do


end subroutine testCollisionTerm
subroutine writeParticleOut(teilchen)
  use particleDefinition
  type(particle) :: teilchen
  write(*,'(2I3,2L,4F6.3,L,3I8)')&
       & teilchen%ID,teilchen%Charge,teilchen%perturbative,teilchen%antiparticle, teilchen%perweight,teilchen%scaleCS, &
       & teilchen%lastCollisionTime, teilchen%productionTime, &
       & teilchen%in_Formation,teilchen%number,teilchen%event 
end subroutine writeParticleOut

!*************************************************************************************************


subroutine test3Body

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
  integer :: i,j,steps
  integer, parameter :: m=1 ! number ensembles
  integer :: numberPions_before,numberPions_after

  real, save :: elab=0.005
  integer, save :: charge=1
  

  write(*,*)
  write( *,*) '******************************************************'
  write( *,*) 'Testing collisionTerm'
  write( *,*) '******************************************************'
  
  call readinputtest

  
  do steps=1,1
     elab=elab+0.005

  Allocate(teilchenPert(1:m,1:20))
  Allocate(teilchenReal(1:m,1:4))

  ! Set up real particles
  do i=1,m
     Do j=1,4
        teilchenReal(i,j)%Id=nucleon
        teilchenReal(i,j)%mass=baryon(nucleon)%mass
        teilchenReal(i,j)%momentum=(/baryon(nucleon)%mass,0.,0.,0./)
        teilchenReal(i,j)%velocity=teilchenReal(i,j)%momentum(1:3)/teilchenReal(i,j)%momentum(0)
        teilchenReal(i,j)%event=0
        teilchenReal(i,j)%perturbative=.false.
        teilchenReal(i,j)%antiparticle=.false.
     end do
     teilchenReal(i,1)%charge=1
     teilchenReal(i,1)%position=(/0.,0.,1./)
     teilchenReal(i,2)%charge=1
     teilchenReal(i,2)%position=(/0.,1.,0./)
     teilchenReal(i,3)%charge=0
     teilchenReal(i,3)%position=(/1.,0.,0./)
     teilchenReal(i,4)%charge=0
     teilchenReal(i,4)%position=(/-1.,0.,0./)
  end do

  ! Set up perturbative particles
  numberPions_before=0
  do i=1,m
     Do j=lBound(teilchenPert,dim=2),NINT(float(uBound(teilchenPert,dim=2))/3.)
        teilchenPert(i,j)%Id=pion
        teilchenPert(i,j)%mass=meson(pion)%mass
        teilchenPert(i,j)%momentum(0)=baryon(nucleon)%mass+elab
        teilchenPert(i,j)%momentum(1:3)=(/0.,0.,SQRT((meson(pion)%mass+elab)**2-meson(pion)%mass**2) /)
        teilchenPert(i,j)%velocity=teilchenPert(i,j)%momentum(1:3)/teilchenPert(i,j)%momentum(0)
        teilchenPert(i,j)%event=-999
        teilchenPert(i,j)%perturbative=.true.
        teilchenPert(i,j)%antiparticle=.false.
        teilchenPert(i,j)%charge=charge
        teilchenPert(i,j)%position=(/0.,0.,0./)
        teilchenPert(i,j)%perweight=float(i)
        numberPions_before=numberPions_before+1
     end do
  end do

  call setNumber(teilchenReal)
  call setNumber(teilchenPert)

  If(m.lt.10) then
     Print *, 'Wurzel(s)=',sqrts(teilchenPert(1,1),teilchenReal(1,1),teilchenReal(1,2))
     Print *, 'Total momentum=',teilchenPert(1,1)%momentum+teilchenReal(1,1)%momentum+teilchenReal(1,2)%momentum
     print * ,'***********************'

     Print *, 'Testing the collisions'
     call timeMeasurement(.true.)
     Print *, 'Before the collisions, Printing perturbative vector'
     Print *, '*************'
     Print *, 'Id,Perturbative, antiparticle, perweight,scaleCS, collTime, prodTime, inFormation, number,event'
     do i=lBound(teilchenPert,dim=1),uBound(teilchenPert,dim=1)
        do j=lBound(teilchenPert,dim=2),uBound(teilchenPert,dim=2)
           call writeParticleOut(teilchenPert(i,j))
        end do
     end do

     Print *, 'Before the collisions, Printing real vector'
     Print *, '*************'
     Print *, 'Id,Perturbative, antiparticle, perweight,scaleCS, collTime, prodTime, inFormation, number,event'
     do i=lBound(teilchenReal,dim=1),uBound(teilchenReal,dim=1)
        do j=lBound(teilchenReal,dim=2),uBound(teilchenReal,dim=2)
           call writeParticleOut(teilchenReal(i,j))
        end do
     end do
     Print *, ' '
     Print *, ' '
     Print *, ' '
     Print *, ' '
  end if
  call collideMain(teilchenPert,teilchenReal,0.5)
  call timeMeasurement(.false.)

  If(m.lt.10) then
     Print *, 'Finished the collisions, Printing perturbative vector'
     Print *, '*************'
     Print *, 'Id,Perturbative, antiparticle, perweight,scaleCS, collTime, prodTime, inFormation, number,event'
     do i=lBound(teilchenPert,dim=2),uBound(teilchenPert,dim=1)
        do j=lBound(teilchenPert,dim=2),uBound(teilchenPert,dim=2)
           call writeParticleOut(teilchenPert(i,j))
        end do
     end do

     Print *, 'Finished the collisions, Printing real vector'
     Print *, '*************'
     Print *, 'Id,Perturbative, antiparticle, perweight,scaleCS, collTime, prodTime, inFormation, number,event'
     do i=lBound(teilchenReal,dim=1),uBound(teilchenReal,dim=1)
        do j=lBound(teilchenReal,dim=2),uBound(teilchenReal,dim=2)
           call writeParticleOut(teilchenReal(i,j))
        end do
     end do
  end if

  numberPions_after=0
  ! Count pions
  do i=lBound(teilchenPert,dim=2),uBound(teilchenPert,dim=1)
     do j=lBound(teilchenPert,dim=2),uBound(teilchenPert,dim=2)
        if (teilchenPert(i,j)%ID.eq.pion) numberPions_after=numberPions_after+1
     end do
  end do

  write(*,*) 
  write(*,*) '# Kinetic energy of pion:', elab
  write(*,*) '# Pions before:', numberPions_before
  write(*,*) '# Pions after:', numberPions_after
  write(*,*) '=> Decay rate:', (log(float(numberPions_before))-log(float(numberPions_after)))/delta_T*197.
  write(111,*) elab,(log(float(numberPions_before))-log(float(numberPions_after)))/delta_T*197.

  DeAllocate(teilchenPert)
  DeAllocate(teilchenReal)
  end do

contains

  subroutine readInputtest
    use output
    implicit none

    NAMELIST /test/ elab,charge

    call Write_ReadingInput('test',0)
    rewind(5)
    read(5,nml=test)
    write(*,*) '  Kinetic energy of pion',elab
    write(*,*) '  Charge of pion', charge

    call Write_ReadingInput('test',1)


  end subroutine readInputtest


end subroutine test3Body





