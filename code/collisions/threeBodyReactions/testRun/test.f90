program test

  use particleProperties, only: initParticleProperties
  implicit none

  call initParticleProperties

  call testMatrixElements
  !call tester
  ! call testGamma

contains

  subroutine testMatrixElements
    use NNPION_NN, only : matrix_NN_NNPion

    integer :: i 
    real :: srts
    real, dimension (1:4) :: matrixElements

    do i=1,2000
       srts=1.8+i*0.001
       matrixElements = matrix_NN_NNPion(srts, .true.)
    end do

  end subroutine testMatrixElements


  !************************************************************************************************


!   subroutine testgamma
!     use NNPION_NN, only : gamma_NNPion_NN,gamma_NNPion_NN_output
! 
!     integer :: i ,charge,charge1,charge2
!     real :: srts
!     real, dimension (1:4) :: gamma
!     real :: mPion, Enucleon(1:2), ePion
!     real :: gamma1
!     integer, dimension(1:2) :: chargeNucleon
! 
! 
!     ENUCLEON=0.938
!     mPion=0.138
!     do charge=-1,1
!        do i=1,200
!           ePion=mPion+i*0.002
!           srts=sqrt((enucleon(1)+enucleon(2)+ePion)**2-(ePion**2-mPion**2))
!           call gamma_NNPion_NN_output(srts,ePION,0.085,0.085,charge, gamma,.true.)
!        end do
!     end do
! 
!     
!     ENUCLEON=0.938
!     do charge=-1,1  ! Loop over pion charge
!        do i=1,200
!           ePion=mPion+i*0.002
!           srts=sqrt((enucleon(1)+enucleon(2)+ePion)**2-(ePion**2-mPion**2))
!           gamma1=0
!           write(*,*) 'i=',i
!           do charge1=0,1 ! Loop over first nucleon charge
!              do charge2=charge1,1 ! Loop over second nucleon charge
!                 chargeNucleon=(/charge1,charge2/)
!                 gamma1=gamma1+gamma_NNPion_NN(srts,ePION,enucleon,0.085,0.085,charge,chargeNucleon)
!              end do
!           end do
!           write(charge+10,*) epion-mPion, gamma1*1000
!        end do
!     end do
! 
!   end subroutine testgamma


  !************************************************************************************************


  subroutine tester
    use particleDefinition
    use constants, only: mN, mPi

    type(particle) :: p1,p2,p3
    real :: srts, pcm
    integer :: i

    p1%mass=mN
    p2%mass=mN
    p3%mass=mPi
    p2%momentum(1:3)=(/0.1,0.,0./)
    p1%momentum(1:3)=(/0.,0.,0./)
    p3%momentum(1:3)=(/0.,0.,0./)

    p1%momentum(0)=p1%mass
    p2%momentum(0)=sqrt(p2%mass**2+Dot_Product(    p2%momentum(1:3),    p2%momentum(1:3)))

    do i=1,100
       p3%momentum(0)=p3%mass+0.005*i
       p3%momentum(1)=sqrt((p3%momentum(0))**2-p3%mass**2)
       srts=sqrtS(p1,p2,p3)
       pcm=SQRT(Max(0.,srts**2/4.-mN**2))
       write(*,'(3(A,F8.3))') 'SQRT(s)=', srts,' PCM=',pcm,'elab=',0.005*i
       write(10,*) srts,pcm,0.025*i
    end do


  end subroutine tester


end program
