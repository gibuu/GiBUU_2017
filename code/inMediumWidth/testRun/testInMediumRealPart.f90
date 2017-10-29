program test
use particleProperties, only: initParticleProperties
use inputGeneral, only: readInputGeneral
implicit none

call initParticleProperties
call readInputGeneral

call testerImagPart
call testerRealPart

contains

subroutine testerRealPart
  ! Test real part of self energy
  use output, only: intToChar
  use selfenergy_baryons, only : get_RealPart
  use particleProperties, only: hadron
  use mediumDefinition

  integer :: i,j,k,particleID
  real,parameter :: rho=0.15
  real :: pabs,mass,E
  type(medium) :: med

  med%useMedium=.true.

  do particleID=1,10
!  do particleID=2,2
     densLoop: do i=4,4
        open(23,file='RealPart_'//IntToChar(particleID)//'_dens'//IntToChar(i)//'.dat')
        write(*,*) "i",i
        med%densityNeutron=rho/4.*float(i)/2.
        med%densityProton=rho/4.*float(i)/2.
        med%density=med%densityNeutron+med%densityProton
        write((i+1)*100+particleId,*) '#rho=',med%density
        write((i+1)*100+particleId,*) '#ID=',particleID
        mass=0.940
        massenLoop: do j=0,50
           write(*,*) "j",j
           mass=hadron(particleID)%minmass-0.2+float(j)*0.02
           impulsLoop: do k=0,20 
              pabs=float(k)*0.2
              E=sqrt(pabs**2+mass**2)
              write(23,'(5E18.6)') mass, pabs, get_RealPart(particleID,pAbs,mass,med)
           end do impulsLoop
           write(23,*)
        end do massenLoop
        close(23)
     end do densLoop
  end do

end subroutine testerRealPart



subroutine testerImagPart
  ! Test imaginary part of self energy
  use selfenergy_baryons, only : selfenergy_imag
  use output, only: intToChar
  use mediumDefinition

  integer :: i,j,k,particleID
  real,parameter :: rho=0.17
  real :: pabs,mass,E
  type(medium) :: med

  med%useMedium=.true.

  do particleID=1,20
     open(23,file='ImagPart_'//IntToChar(particleID)//'.dat')
     densLoop: do i=4,4
        write(*,*) "i",i
        med%densityNeutron=rho/4.*float(i)/2.
        med%densityProton=rho/4.*float(i)/2.
        med%density=med%densityNeutron+med%densityProton
        write((i+1)*100+particleId,*) '#rho=',med%density
        write((i+1)*100+particleId,*) '#ID=',particleID
        mass=0.940
        massenLoop: do j=0,200
           write(*,*) "j",j
           mass=0.+float(j)*0.04
           impulsLoop: do k=0,20 
              pabs=float(k)*0.2
              E=sqrt(pabs**2+mass**2)
              write(23,'(5E18.6)') mass, pabs, selfenergy_imag(particleID,pAbs,E,med)
           end do impulsLoop
           write(23,*)
        end do massenLoop
     end do densLoop
     close(23)
  end do

end subroutine testerImagPart


end program test
