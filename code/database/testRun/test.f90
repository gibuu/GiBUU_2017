program test
  use particleProperties
  use IdTable
  use particleDefinition

  implicit none
  integer :: i
  type(particle) :: teilchen
  write(*,*) "**********************************"
  write(*,*) "**********************************"
  write(*,*)
  Write(*,*)      "Testing the database "
  write(*,*)
  write(*,*) "**********************************"
  write(*,*) "**********************************"
  write(*,*)
  Write(*,*)      "Calling'init_Database'"
  call init_Database
  write(*,*) "Testing minimal masses. Writing results to fort.100"
  Do i=lbound(meson,dim=1),ubound(meson,dim=1)
     write(100,*) i, minimalMass(i)
  End do
  Do i=lbound(baryon,dim=1),ubound(baryon,dim=1)
     write(100,*) i, minimalMass(i)
  end do

  write(*,*) "**********************************"

  write(*,*) "**Testing IsMeson"


  write(*,*) 'ISMESON(100)', isMeson(100)
  write(*,*) 'ISMESON(120)', isMeson(120)
  write(*,*) 'ISMESON(130)', isMeson(130)
  write(*,*) 'ISMESON(200)', isMeson(140)
  write(*,*) 'ISBARYON(1)',  isBaryon(1)
  write(*,*) 'ISBARYON(3)',  isBaryon(3)
  write(*,*) 'ISBARYON(100)',isBaryon(100)

  write(*,*) "**********************************"

  write(*,*) "**Testing ValidCharge"


  teilchen%id=nucleon
  teilchen%charge=1
  teilchen%antiParticle=.true.
  write(*,*) 'ValidCharge(',teilchen%ID,teilchen%charge,teilchen%antiparticle,')',ValidCharge(teilchen)
  teilchen%antiParticle=.false.
  write(*,*) 'ValidCharge(',teilchen%ID,teilchen%charge,teilchen%antiparticle,')',ValidCharge(teilchen)

  teilchen%id=kaonBar
  teilchen%charge=1
  teilchen%antiParticle=.true.
  write(*,*) 'ValidCharge(',teilchen%ID,teilchen%charge,teilchen%antiparticle,')',ValidCharge(teilchen)
  teilchen%antiParticle=.false.
  write(*,*) 'ValidCharge(',teilchen%ID,teilchen%charge,teilchen%antiparticle,')',ValidCharge(teilchen)




  write(*,*) "**********************************"
  Write(*,*)      "Finished Testing"
  write(*,*) "**********************************"



end program test
