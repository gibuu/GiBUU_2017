program test
  use inputGeneral, only: readInputGeneral
  use particleProperties, only: initParticleProperties
  implicit none

  call readInputGeneral
  call initParticleProperties

  call test_vac

!   call test_oset

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! test the Delta vacuum width
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine test_vac
    use constants, only: pi
    use particleProperties, only: hadron
    use IDTable, only: Delta
    use baryonWidth, only: FullWidthBaryon
    integer :: i
    real :: m,w,sf

    open(99,file="DeltaVacWidth.dat")
    do i=0,150
      m = hadron(Delta)%minmass + i*0.01
      w = FullWidthBaryon (Delta, m)
      sf = w/pi*m / ( (m**2-hadron(Delta)%mass**2)**2 + (m*w)**2)
      write (99,*) m,w,sf
    end do
    close(99)
  end subroutine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! test the Oset in-medium width
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine test_oset
    use deltaWidth, only: deloset
    real :: mass
    real :: imsig2,imsig3,imsigq
    real :: gcoll1,gcoll2,gcoll1a
    integer :: i

    do i=0,100
      mass = 1. + i*0.01
      call deloset(mass,0.17,imsig2,imsig3,imsigq)
      gcoll1=2.*(imsig2+imsig3)
      gcoll1a=2.*imsig2
      gcoll2=2.*imsigq
      write(100,'(7ES15.5)') mass,imsig2,imsig3,imsigq,gcoll1,gcoll1a,gcoll2
    end do

  end subroutine


end program test
