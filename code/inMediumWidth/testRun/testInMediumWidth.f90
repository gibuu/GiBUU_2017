program test
use inputGeneral, only: readInputGeneral
use particleProperties, only: initParticleProperties

call initParticleProperties
call readInputGeneral
call tester

end program test


subroutine tester
  use output
  use deltaWidth
  use baryonWidthMedium_tables, only : get_inMediumWidth
  use IdTable, only : nres,delta
  use potentialModule, only : massDetermination
  use mediumDefinition
  implicit none
  integer :: i,j,k,particleID
  real,parameter :: rho=0.17
!  real,parameter :: rho=0.001
  real :: pabs,mass,baremass
  real, dimension(0:3) :: momentum
  type(medium) :: med

  med%useMedium=.true.

  particleID=1
  Open(57,File='Width_001.dat')
  do i=4,4
     write(*,*) "i",i
     med%densityNeutron=rho/4.*float(i)/2.
     med%densityProton=rho/4.*float(i)/2.
     med%density=med%densityNeutron+med%densityProton
     write((i+1)*100+particleId,*) '#rho=',med%density
     write((i+1)*100+particleId,*) '#ID=',particleID
     mass=0.940
     do j=0,70
        write(*,*) "j",j
        mass=0.8+float(j)*0.005
        do k=0,30 
           pabs=float(k)*0.03
           momentum(0)=sqrt(mass**2+pabs**2)
           momentum(1)=pabs
           momentum(2:3)=0.
           call massDetermination(particleID,momentum,med,baremass)
           write(57,'(5E18.6)') pabs,mass, &
                                get_inMediumWidth(particleID,momentum,baremass,med%densityNeutron,med%densityProton,1), &
                                get_inMediumWidth(particleID,momentum,baremass,med%densityNeutron,med%densityProton,2), &
                                get_inMediumWidth(particleID,momentum,baremass,med%densityNeutron,med%densityProton,3)
        end do
        write(57,*)
     end do
  end do
  close(57)
!  stop

  do particleID=2,31
     Open(57,File='Width_'//IntToChar(particleID)//'.dat')
     do i=4,4
        med%densityNeutron=rho/4.*float(i)/2.
        med%densityProton=rho/4.*float(i)/2.
        med%density=med%densityNeutron+med%densityProton
        write((i+1)*100+particleId,*) '#rho=',med%density
        write((i+1)*100+particleId,*) '#ID=',particleID
        write((i+1)*100+particleId,*) '# pabs, mass, coll width, free width, full width'
        do j=0,100
        write(*,*) "j",j
           mass=0.8+float(j)*0.02
           do k=0,30 
              pabs=float(k)*0.03
              momentum(0)=sqrt(mass**2+pabs**2)
              momentum(1)=pabs
              momentum(2:3)=0.
              call massDetermination(particleID,momentum,med,baremass)
              write(57,'(5E18.6)') pabs,mass, &
                                   get_inMediumWidth(particleID,momentum,baremass,med%densityNeutron,med%densityProton,1), &
                                   get_inMediumWidth(particleID,momentum,baremass,med%densityNeutron,med%densityProton,2), &
                                   get_inMediumWidth(particleID,momentum,baremass,med%densityNeutron,med%densityProton,3)

           end do
           write(57,*)
        end do
     end do
     close(57)
  end do




  Open(57,File='Width_DeltaOset.dat')
  do i=4,4
     med%densityNeutron=rho/4.*float(i)/2.
     med%densityProton=rho/4.*float(i)/2.
     med%density=med%densityNeutron+med%densityProton
     write((i+1)*100+particleId,*) '#rho=',med%density
     write((i+1)*100+particleId,*) '#ID=',particleID
     write((i+1)*100+particleId,*) '# pabs, mass, coll width, free width, full width'
     do j=0,100
        write(*,*) "j",j
        mass=0.8+float(j)*0.02
        do k=0,30 
           pabs=float(k)*0.03
           momentum(0)=sqrt(mass**2+pabs**2)
           momentum(1)=pabs
           momentum(2:3)=0.
           call massDetermination(particleID,momentum,med,baremass)
           write(57,'(6E18.6)') pabs,mass, get_inMediumWidth(delta,momentum,baremass,med%densityNeutron,med%densityProton,1), &
                                           get_inMediumWidth(delta,momentum,baremass,med%densityNeutron,med%densityProton,2), &
                                           get_inMediumWidth(delta,momentum,baremass,med%densityNeutron,med%densityProton,3), &
                                           delta_fullWidth(baremass,(/pabs,0.,0./),med%density)
        end do
        write(57,*)
     end do
  end do
  close(57)










end subroutine tester
