program printWidth
  use particleProperties
  use idTable
  use inputGeneral

  implicit none

  call readinputGeneral
  call initParticleProperties

  call print

end program printWidth


subroutine print
  use mediumDefinition
  use mesonWidthMedium,  only : GetMassAssInfo_Meson, WidthMesonMedium

   use fgsl, only: kn => fgsl_sf_bessel_kcn

  real,dimension(0:3)  :: momLRF


  integer :: ID, iM
  real :: mass

  ID = 103
  do iM=1,500
     mass = iM*0.01

     write(115,*) mass,WidthMesonMedium(ID,mass,momLRF,vacuum),kn(2,mass/0.15),kn(2,mass/0.21)


  end do

end subroutine print
