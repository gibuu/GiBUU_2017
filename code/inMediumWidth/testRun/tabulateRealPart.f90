program tabulateRealPart
  use inputGeneral, only: readInputGeneral
  use version, only: PrintVersion
  use particleProperties, only: initParticleProperties
  use selfenergy_baryons, only: get_RealPart
  use mediumDefinition, only: vacuum
  use IdTable, only: nucleon
  implicit none

  real :: pabs,mass,E,dummy

  call initParticleProperties
  call readInputGeneral

  call PrintVersion

  mass=1.
  pabs=.5
  E=sqrt(pabs**2+mass**2)
  !the values of the above variables are totally irrelevant
  !important is only, that get_realPart is called at all
  !once it has started, it will tabulate (cf. jobcard_realPartTabulate for the correct switches) 

  dummy=get_RealPart(nucleon,pAbs,mass,vacuum)

end program tabulateRealPart
