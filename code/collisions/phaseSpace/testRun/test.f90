program test
use particleProperties, only: baryon, meson, initParticleProperties
implicit none
!!$Print *,"Test 3 body phase space"
!!$call test3Body
!!$
!!$Print *
!!$Print *
!!$Print *
!!$Print *
!!$
!!$Print *,"Test 4 body phase space"
!!$call test4Body
!!$
!!$
!!$Print *
!!$Print *
!!$Print *
!!$Print *
!!$
Print *, 'Initializing databench for particles'
call initParticleProperties
!!$
!!$
Print *,"Test 3 body phase space integration"
call test3Body_Integration
!!$
!!$Print *,"Test 2 body phase space integration"
!!$call test2Body_Integration

!Print *,"Test final State: massAss"
!call testMassAss
!write(*,*) '***********************************************'
!write(*,*) '***********************************************'
!write(*,*) '***********************************************'
!write(*,*) '***********************************************'
!write(*,*) '***********************************************'
!write(*,*) '***********************************************'
!call testAssMass



Print *
Print *
Print *
Print *



contains

subroutine testMassAss
  use finalStateModule
  use MediumDefinition
  use IDTable
  use particleDefinition
  use particleProperties
    implicit none
    real                                     :: srts                             ! sqrt(s)
    type(medium)                           :: medium_AtCollision   ! medium information : density, temperature,...
    integer, dimension(1:2)             :: idOUT                         ! ids of produced particles
    integer, dimension(1:2)             :: chargeOUT                ! charges of produced particles
    real, dimension(1:2)                  :: spotOUT                    ! scalar potential of produced particles
    type(particle), dimension (1:2)  :: pairIn                         ! incoming particles a and b
    real, dimension(1:3)                 :: betaToLRF                 ! beta for boost to LRF
    real, dimension(1:3)                 :: betaToCM                  ! beta for boost to CM-Frame

    type(particle), dimension (1:2)  :: pairOut ! outgoing particles c and d
    logical                                        :: flag       ! .true. if it was possible to determine a final state
    real :: plab
    integer :: i
    real, parameter :: deltaP=0.05
   medium_AtCollision%useMedium=.false.

  idOut=(/nucleon, phi/)
  chargeOut=(/0,0/)
  spotOut=0.

  pairIN(1)%ID=nucleon
  pairIN(2)%ID=pion
  pairIN(1)%mass=baryon(nucleon)%mass
  pairIN(2)%mass=meson(pion)%mass
  pairIN(1)%charge=0
  pairIN(2)%Charge=0
  pairIN(1)%momentum(0)=baryon(nucleon)%mass
  pairIN(1)%momentum(1:3)=0.


  Do i=1,100

     write(*,*)
     write(*,*)
     write(*,*) '***********************************************************************'
     plab=0.1+i*deltaP

     pairIN(2)%momentum(0)=sqrt(pairIN(2)%mass**2+pLab**2)
     pairIN(2)%momentum(1:3)=(/0.,0.,pLab/)
     srts=SQRT((Sum(pairIN%momentum(0)))**2-Dot_Product(pairIN(2)%momentum(1:3),pairIN(2)%momentum(1:3)))
     write(*,*) 'srts vorher=', srts

     If(srts.gt.3.) exit

     betaToLRF=0.
     betaToCM=(/0.,0.,pLab/)/(pairIN(1)%momentum(0)+pairIN(2)%momentum(0))

     Write(*,*) pairIN(:)%ID
     Write(*,'(3F8.3)') pairIN(1)%Momentum(1:3)
     Write(*,'(3F8.3)') pairIN(2)%Momentum(1:3)
     Write(*,*) 'Masses IN' , pairIN(1)%Mass , pairIN(2)%Mass
     Write(*,*) 'Charges IN' ,pairIN(1)%Charge, pairIN(2)%Charge

     write(*,*) 'P vorher=', pLab

     call massass(srts,medium_AtCollision,pairIn,pairout,betaToLRF,betaToCM,flag)

     If (flag) then
        Write(*,*) pairOut(:)%ID
        Write(*,'(3F8.3)') pairOut(1)%Momentum(1:3)
        Write(*,'(3F8.3)') pairOut(2)%Momentum(1:3)
        Write(*,*) 'Masses OUT' , pairOut(1)%Mass ,  pairOut(2)%Mass
        Write(*,*) 'Charges OUT', pairOut(1)%Charge , pairOut(2)%Charge

        pairOut(1)%Momentum(0)=sqrt( pairOut(1)%mass**2+Dot_Product(pairOut(1)%Momentum,pairOut(1)%Momentum))
        pairOut(2)%Momentum(0)=sqrt( pairOut(2)%mass**2+Dot_Product(pairOut(2)%Momentum,pairOut(2)%Momentum))

        betaToCM=-betaToCM
        call lorentz(betaTOCM,pairOUT(1)%momentum)
        call lorentz(betaTOCM,pairOUT(2)%momentum)

        srts=sqrt(Sum(pairOUT%momentum(0))**2-Dot_Product(pairOUT(2)%momentum(1:3)+pairOUT(1)%momentum(1:3), &
                                                          pairOUT(2)%momentum(1:3)+pairOUT(1)%momentum(1:3)))
        write(*,*) 'P nachher=', pairOUT(2)%momentum(1:3)+pairOUT(1)%momentum(1:3)
        write(*,*) 'srts nachher=' , srts
     else
        Write(*,*) 'Failed to find solution!!!'
     end if
  End do

  ! ********************************************************************************************************************************

end subroutine testMassAss

subroutine testAssMass
  use finalStateModule
  use MediumDefinition
  use IDTable
  use particleDefinition
  use particleProperties
  implicit none
  real                                     :: srts                             ! sqrt(s)
  type(medium)                           :: medium_AtCollision   ! medium information : density, temperature,...
  integer, dimension(1:3)             :: idOUT                         ! ids of produced particles
  integer, dimension(1:3)             :: chargeOUT                ! charges of produced particles
  real, dimension(1:3)                  :: spotOUT                    ! scalar potential of produced particles
  type(particle), dimension (1:2)  :: pairIn                         ! incoming particles a and b
  real, dimension(1:3)                 :: betaToLRF                 ! beta for boost to LRF
  real, dimension(1:3)                 :: betaToCM                  ! beta for boost to CM-Frame

  type(particle), dimension (1:3)  :: pairOut ! outgoing particles c and d
  logical                                        :: flag       ! .true. if it was possible to determine a final state
  real :: plab
  integer :: i
  real, parameter :: deltaP=0.05
  medium_AtCollision%useMedium=.false.

  idOut=(/nucleon, phi, pion /)
  chargeOut=(/0,0,0/)
  spotOut=0.

  pairIN(1)%ID=nucleon
  pairIN(2)%ID=pion
  pairIN(1)%mass=baryon(nucleon)%mass
  pairIN(2)%mass=meson(pion)%mass
  pairIN(1)%charge=0
  pairIN(2)%Charge=0
  pairIN(1)%momentum(0)=baryon(nucleon)%mass
  pairIN(1)%momentum(1:3)=0.


  Do i=1,100

     write(*,*)
     write(*,*)
     write(*,*) '***********************************************************************'
     plab=0.1+i*deltaP

     pairIN(2)%momentum(0)=sqrt(pairIN(2)%mass**2+pLab**2)
     pairIN(2)%momentum(1:3)=(/0.,0.,pLab/)
     srts=SQRT((Sum(pairIN%momentum(0)))**2-Dot_Product(pairIN(2)%momentum(1:3),pairIN(2)%momentum(1:3)))
     write(*,*) 'srts vorher=', srts

     If(srts.gt.3.) exit

     betaToLRF=0.
     betaToCM=(/0.,0.,pLab/)/(pairIN(1)%momentum(0)+pairIN(2)%momentum(0))

     Write(*,*) pairIN(:)%ID
     Write(*,'(3F8.3)') pairIN(1)%Momentum(1:3)
     Write(*,'(3F8.3)') pairIN(2)%Momentum(1:3)
     Write(*,*) 'Masses IN' , pairIN(1)%Mass , pairIN(2)%Mass
     Write(*,*) 'Charges IN' ,pairIN(1)%Charge, pairIN(2)%Charge
     write(*,*) 'P vorher=', pLab


     call assMass(srts,medium_AtCollision,pairIn,pairOut,spotOut,betaToLRF,betaToCM,flag)

     If (flag) then
        Write(*,*) pairOut(:)%ID
        Write(*,'(3F8.3)') pairOut(1)%Momentum(1:3)
        Write(*,'(3F8.3)') pairOut(2)%Momentum(1:3)
        Write(*,'(3F8.3)') pairOut(3)%Momentum(1:3)
        Write(*,*) 'Masses OUT' , pairOut(1)%Mass ,  pairOut(2)%Mass,  pairOut(3)%Mass
        Write(*,*) 'Charges OUT', pairOut(1)%Charge , pairOut(2)%Charge,   pairOut(3)%Mass

        pairOut(1)%Momentum(0)=sqrt( pairOut(1)%mass**2+Dot_Product(pairOut(1)%Momentum(1:3),pairOut(1)%Momentum(1:3)))
        pairOut(2)%Momentum(0)=sqrt( pairOut(2)%mass**2+Dot_Product(pairOut(2)%Momentum(1:3),pairOut(2)%Momentum(1:3)))
        pairOut(3)%Momentum(0)=sqrt( pairOut(3)%mass**2+Dot_Product(pairOut(3)%Momentum(1:3),pairOut(3)%Momentum(1:3)))

        betaToCM=-betaToCM
        call lorentz(betaTOCM,pairOUT(1)%momentum)
        call lorentz(betaTOCM,pairOUT(2)%momentum)
        call lorentz(betaTOCM,pairOUT(3)%momentum)

        srts = sqrt(Sum(pairOUT%momentum(0))**2 - &
                    Dot_Product(pairOUT(3)%momentum(1:3)+pairOUT(2)%momentum(1:3)+pairOUT(1)%momentum(1:3), &
                                pairOUT(3)%momentum(1:3)+pairOUT(2)%momentum(1:3)+pairOUT(1)%momentum(1:3)))
        write(*,*) 'P nachher=', pairOUT(2)%momentum(1:3)+pairOUT(1)%momentum(1:3)+pairOUT(3)%momentum(1:3)
        write(*,*) 'srts nachher=' , srts
     else
        Write(*,*) 'Failed to find solution!!!'
     end if
  End do


end subroutine testAssMass

!**********************************



subroutine test3body
use nBodyPhaseSpace, only : momenta_in_3BodyPS
implicit none

real :: m1,m2,m3,srts
real, dimension(3,3) :: p3
real,parameter :: deltaS=0.1
integer :: i
m1=0.12
m2=0.9
m3=1.2
Print *,"Masses", m1, m2, m3
Do i=0,100
   srts=i*deltaS
   call momenta_in_3BodyPS(srts,(/m1,m2,m3/),p3)
End do


m1=0.12
m2=0.000000000000001
m3=1.2
Print *,"Masses", m1, m2, m3
Do i=0,100
   srts=i*deltaS
   call momenta_in_3BodyPS(srts,(/m1,m2,m3/),p3)
End do


m1=0.12
m2=0.9
m3=300.
Print *,"Masses", m1, m2, m3
Do i=0,100
   srts=i*deltaS
   call momenta_in_3BodyPS(srts,(/m1,m2,m3/),p3)
End do
end subroutine test3body




subroutine test4body
use nBodyPhaseSpace, only : momenta_in_4BodyPS
implicit none

real :: m1,m2,m3,m4,srts
real, dimension(3,4) :: p4
real,parameter :: deltaS=0.1
integer :: i

m4=0.12

m1=0.12
m2=0.9
m3=1.2
Print *,"Masses", m1, m2, m3
Do i=0,100
   srts=i*deltaS
   call momenta_in_4BodyPS(srts,(/m1,m2,m3,m4/),p4)

End do


m1=0.12
m2=0.000000000000001
m3=1.2
Print *,"Masses", m1, m2, m3
Do i=0,100
   srts=i*deltaS
   call momenta_in_4BodyPS(srts,(/m1,m2,m3,m4/),p4)
End do


m1=0.12
m2=0.9
m3=300.
Print *,"Masses", m1, m2, m3
Do i=0,100
   srts=i*deltaS
   call momenta_in_4BodyPS(srts,(/m1,m2,m3,m4/),p4)

End do
end subroutine test4body


!**************************************************


subroutine test3body_Integration
use threeBodyPhaseSpace
use idTable
implicit none

real :: m1,m2,m3,srts
real,parameter :: deltaS=0.025
integer :: i
real :: ps,scalarPotential
real , dimension(1:2) :: psR
integer :: resonanceID

m1=0.12
m2=0.9
!m3=baryon(delta)%mass
m3=meson(sigmaMeson)%mass

!resonanceId=delta
resonanceId=sigmaMeson

scalarPotential=0.0 
Print *,"Masses", m1, m2, m3
Do i=0,1000
   srts=i*deltaS
   call Integrate_3bodyPS(ps,srts,m1,m2,m3)
   call Integrate_3bodyPS_Resonance(psR,srts,m1,m2,resonanceID,scalarPotential)
   Write(701,'(4F12.4)') srts,ps,psR
End do


m1=0.0000000000001
m2=0.000000000000001

Print *,"Masses", m1, m2, m3
Do i=0,1000
   srts=i*deltaS
   call Integrate_3bodyPS(ps,srts,m1,m2,m3)
   call Integrate_3bodyPS_Resonance(psR,srts,m1,m2,resonanceID,scalarPotential)
   Write(702,'(4F12.4)') srts,ps,psR
End do


m1=0.12
m2=8.9

Print *,"Masses", m1, m2, m3
Do i=0,1000
   srts=i*deltaS
   call Integrate_3bodyPS(ps,srts,m1,m2,m3)
   call Integrate_3bodyPS_Resonance(psR,srts,m1,m2,resonanceID,scalarPotential)
   Write(703,'(4F12.4)') srts,ps,psR
End do
end subroutine test3body_Integration



!**************************************************


subroutine test2body_Integration
use twoBodyPhaseSpace, only: Integrate_2bodyPS_resonance, nnRR
use idTable
implicit none

real :: m1,srts
real,parameter :: deltaS=0.025
integer :: i
real ::scalarPotential
real , dimension(1:5) :: ps
integer :: ID(1:2)
real :: integral


Id = (/delta,delta/)

m1=baryon(delta)%mass

scalarPotential=0.0 
Print *,"Masses", m1
Do i=0,1000
   srts=i*deltaS
   call  Integrate_2bodyPS_resonance(ps,ID(1),srts,m1,scalarPotential)
   integral = nnRR(srts,ID)
   Write(901,'(9F10.4)') srts,ps,integral
End do

Id(1)=SigmaStar
m1=0.1
scalarPotential=0.0 
Print *,"Masses", m1
Do i=0,1000
   srts=i*deltaS
   call  Integrate_2bodyPS_resonance(ps,ID(1),srts,m1,scalarPotential)
   Write(902,'(7F12.4)') srts,ps
End do

Id(1)=p11_1710
m1=0.1
scalarPotential=0.0 
Print *,"Masses", m1
Do i=0,1000
   srts=i*deltaS
   call  Integrate_2bodyPS_resonance(ps,ID(1),srts,m1,scalarPotential)
   Write(903,'(7F12.4)') srts,ps
End do

end subroutine test2body_Integration

end program test
