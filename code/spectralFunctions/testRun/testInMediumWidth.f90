program test
use inputGeneral
use particleProperties, only: initParticleProperties

call initParticleProperties
call readInputGeneral
call tester

end program test


subroutine tester
  use output
  use mediumDefinition
  use baryonWidthMedium,only:WidthBaryonMedium
  use IdTable, only : nres
  implicit none
  integer :: i,j,k,particleID
  real,parameter :: rho=0.16
!  real,parameter :: rho=0.001
  real :: rhoN,rhoP,pabs,mass
  type(medium) :: med
  real, dimension(0:3) :: momentumLRF
  real :: dummy
  particleID=1
  do i=1,4
     write(*,*) "i",i
     Open(57,File='Width_'//intToChar(particleID)//'_rho_'//intToChar(i)//'.dat')
     rhoN=rho/4.*float(i)/2.
     rhoP=rho/4.*float(i)/2.
     call setMedium()
     write(57,*) '#rho=',rhoN+rhoP
     write(57,*) '#ID=',particleID
     mass=0.940
     do j=0,70
        write(*,*) "j",j
        mass=0.8+float(j)*0.005
        do k=0,34 
           pabs=float(k)*0.03
           momentumLRF=(/sqrt(mass**2+pabs**2),pabs,0.,0./)
           write(57,'(3E18.6)') pabs,mass, WidthBaryonMedium(particleID,mass,momentumLRF,med)
        end do
        write(57,*)
     end do
     close(57)
  end do

!  stop

  do particleID=2,31
     do i=0,4
        Open(57,File='Width_'//intToChar(particleID)//'_rho_'//intToChar(i)//'.dat')
        rhoN=rho/4.*float(i)/2.
        rhoP=rho/4.*float(i)/2.
        call setMedium()
        write(57,*) '#rho=',rhoN+rhoP
        write(57,*) '#ID=',particleID
        do j=0,100
        write(*,*) "j",j
           mass=0.8+float(j)*0.02
           do k=0,34 
              pabs=float(k)*0.03
              momentumLRF=(/sqrt(mass**2+pabs**2),pabs,0.,0./)
              write(57,'(3E18.6)') pabs,mass, WidthBaryonMedium(particleID,mass,momentumLRF,med)
           end do
           write(57,*)
        end do
        close(57)
     end do
  end do


  ! TIMINGS
  write(*,*) 'Start timing'
  call timeMeasurement()
  do j=0,700
     !write(*,*) "j",j
     mass=0.8+float(j)*0.0005
     do k=0,30000
        pabs=float(k)*0.00003
        momentumLRF=(/sqrt(mass**2+pabs**2),pabs,0.,0./)
        dummy= WidthBaryonMedium(23,mass,momentumLRF,med)
     end do

  end do
  call timeMeasurement()



  contains 
    subroutine setMedium()
      implicit none
      med%useMedium=.true.
      med%densityProton=rhoP
      med%densityNeutron=rhoP
      med%density=rhoP+rhoN
      med%temperature=0.
    end subroutine setMedium

end subroutine tester
