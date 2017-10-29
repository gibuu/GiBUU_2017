! Testing the module parametrizationBarMes

program test_param
  use inputGeneral, only: readInputGeneral
  use version, only: PrintVersion
  use particleProperties, only: initParticleProperties
  use parametrizationsBarMes, only :  piN_to_strangeBaryon_kaon_pion

  implicit none

  real :: sqrts
  real, dimension(1:7) :: sigma_piMinus_p
  real, dimension(1:4) :: sigma_piPlus_p

  call PrintVersion

  call readInputGeneral           ! in module inputGeneral
  call initParticleProperties


  open(10,file="piN_to_strangeBaryon_kaon_pion.dat")
  sqrts=0.
  write(10,*) '# sqrt(s), Cross sections pi^-p (1:7),Cross sections pi^+p (1:4)'
  do
     call  piN_to_strangeBaryon_kaon_pion(sqrtS,sigma_piMinus_p,sigma_piPlus_p)
     write(10,'(12F17.9)') sqrtS,sigma_piMinus_p,sigma_piPlus_p
     sqrts=sqrts+0.01
     if(sqrts.gt.10.) exit
  end do
  close(10)
  write(*,*) 'FINISHED TESTING piN_to_strangeBaryon_kaon_pion!!'


end program test_param
