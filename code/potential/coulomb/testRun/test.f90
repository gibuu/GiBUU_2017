program testCoulomb
!program to test the Coulomb routines
use particleDefinition
use idTable, only: nucleon
use coulomb, only: updateCoulomb, emfoca

implicit none
type(particle),dimension(1:1,1:3) :: teilchenVektor
real :: x,y,z,cpot
real,dimension(1:3) ::emForce,r,p
integer :: i !,j
real, parameter :: elmcon=0.0014398  !=e^2*0.2 Gev fm=1/137*0.2 GeV fm

teilchenVektor(1,1)%Charge=1
teilchenVektor(1,1)%position=(/0.,0.4,0./)
teilchenVektor(1,1)%momentum=(/0.14,0.,0.,0./)
teilchenVektor(1,1)%ID=nucleon

!!$!teilchenVektor(2,1)%Charge=0
!!$teilchenVektor(2,1)%position=(/100.,0.,0./)
!!$teilchenVektor(2,1)%momentum=(/0.14,0.,0.,0./)
!!$teilchenVektor(2,1)%ID=1


write(*,*) 'update and initializing Coulomb'
call updateCoulomb  !(teilchenVektor,.true.)
write(*,*) 'Testing Coulomb'
Do i=-25,25
   !Do j=-25,25
      x=i*2.
      y=0!j*2.
      z=0.
      r=(/x,y,z/)
      p=(/0.,0.,0./)
      cpot = emfoca(r,p,1,nucleon,emForce)
      Write(11,'(4F15.7)') x,y,cpot,1./sqrt(x**2+y**2)*elmcon
   !End do
      !write(11,*) 
End do


end program testCoulomb
