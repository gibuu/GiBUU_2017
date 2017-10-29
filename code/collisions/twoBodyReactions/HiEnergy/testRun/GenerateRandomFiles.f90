program test
  use output
  use particleDefinition
  use particleProperties
  use Coll_Fritiof
  use Coll_Pythia
  !use rhoMassParameter, only : srtFreeRhoMass
  use hadronFormation, only : forceInitFormation
  use random

  implicit none


  integer :: i,N
  real :: r

  real qyr

  call forceInitFormation
  call InitParticleProperties

  N = rn(0)*100000*qyr(0)

  write(*,*) N
  write(*,*) qyr(0)

  do i=1,N
     r = qyr(0)
  enddo

  call DumpRANDOM



end program test


