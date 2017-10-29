program test
  use inputGeneral, only: readinputGeneral
  use particleProperties, only: initParticleProperties
  use version, only: PrintVersion
  use inMediumWidth

  call PrintVersion

  call readinputGeneral
  call initParticleProperties

  call tabulate_inMediumWidth_mesons
  !call tabulate_inMediumWidth_baryons

end program test
