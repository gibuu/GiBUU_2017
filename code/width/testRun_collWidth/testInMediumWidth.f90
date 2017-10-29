program test
use inputGeneral
use particleProperties, only: initParticleProperties

call initParticleProperties
call readInputGeneral
call tester

end program test


subroutine tester
  use output
  use baryonWidthMedium_tables, only : get_inMediumWidth
  use IdTable, only : nres
  implicit none
  integer :: i,j,k,particleID
  real,parameter :: rho=0.16
!  real,parameter :: rho=0.001
  real :: rhoN,rhoP,pabs,mass
  real, dimension(0:3) :: momentum

  particleID=1
  do i=1,4
     Open(57,File='Width_'//inttochar(particleID)//'_rho_'//intToChar(i)//'.dat')
     write(*,*) "i",i
     rhoN=rho/4.*float(i)/2.
     rhoP=rho/4.*float(i)/2.
     write(57,*) '#rhoN', rhoN
     write(57,*) '#rhoP', rhoP

     write(57,*) '#rho=',rhoN+rhoP
     write(57,*) '#ID=',particleID
     mass=0.938
     do k=0,100 
        pabs=float(k)*0.01
        momentum=(/sqrt(mass**2+pabs**2),pabs,0.,0./)
        write(57,'(5E18.6)') pabs,mass, get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,1)&
             , get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,2), get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,3)
     end do
     close(57)
  end do
 

  do particleID=2,2     
     do i=0,4
        Open(57,File='Width_'//IntToChar(particleID)//'_rho_'//intToChar(i)//'.dat')
        rhoN=rho/4.*float(i)/2.
        rhoP=rho/4.*float(i)/2.
        write(57,*) '#rho=',rhoN+rhoP
        write(57,*) '#ID=',particleID
        do j=0,100
        write(*,*) "j",j
           mass=0.8+float(j)*0.02
           do k=0,34 
              pabs=float(k)*0.03
              momentum=(/sqrt(mass**2+pabs**2),pabs,0.,0./)
              write(57,'(5G18.6)') pabs,mass, get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,1)&
                   & , get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,2),&
                   & get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,3)

           end do
           write(57,*)
        end do
     end do
     close(57)
  end do


end subroutine tester
