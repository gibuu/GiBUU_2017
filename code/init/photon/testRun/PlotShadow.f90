program PlotShadow
  use nucleusDefinition
  use nucleus
  use shadowing

  implicit none

  type(tNucleus),pointer :: Nuc
  real,dimension(1:3) :: position
  real :: nu,Q2
  real,dimension(4) :: shadfac
  real :: z, dz, r

  integer:: i,nZ


  Nuc => getTarget()

  nZ = 50
  nu = 10.0
  Q2 = 1.!1.5

  write(*,*) Nuc%MaxDist, Nuc%Radius

  call WriteNucleusStaticDens("DensTab.dat", Nuc)


  z = -Nuc%MaxDist
  dz = -z/nZ
  do i=-nZ,nZ
     position = (/0.0, 0.0, z/)
     r = sqrt(position(1)**2+position(2)**2+position(3)**2)
     
     call AeffCalc(Nuc,position,nu,Q2,shadfac)

     write(91,*) z, shadfac(1), NucleusStaticDens(Nuc,r,0)

     z = z+dz

  enddo
  


end program PlotShadow
