program Compton

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!! Computes the cross section for Real Compton Scattering (RCS) via resonances: !!!!!
  !!!!! gamma N -> R -> gamma N                                                      !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use IDTable, only: nucleon,nres
  use particleDefinition
  use particleProperties, only: initParticleProperties
  use constants, only: mN
  use resProd_lepton, only: sigma_resProd
  use baryonWidth, only: FullWidthBaryon, baryonWidth_gammaN

  implicit none

  real, parameter :: delta_E = 0.01
  integer, parameter :: res_max = 10

  integer :: i,res_ID
  real :: E_gamma,p_gamma(0:3),mass_res,gamma_photon
  real,dimension(2:nres+1) :: sigma1,sigma2,BR
  type(particle) :: nuc

  call initParticleProperties

  call setToDefault (nuc)
  nuc%ID = nucleon
  nuc%charge = 1
  nuc%mass = mN
  nuc%momentum = (/mN,0.,0.,0./)

  open (80,file="gammaN_R.dat")
  open (81,file="RCS.dat")
  open (82,file="BR.dat")

  do i=10,200
    E_gamma = i*delta_E ! photon energy
    p_gamma = (/E_gamma,0.,0.,E_gamma/)  ! real photon along z-axis
    do res_ID=2,nres+1  ! loop over resonances
      sigma1 (res_ID) = sigma_resProd(nuc,res_ID,p_gamma,mass_res)*1E3 ! cross section in microbarn
      if (sigma1 (res_ID) > 0.) then
        gamma_photon = baryonWidth_gammaN(res_ID,mass_res,0.,nuc%charge)
        BR (res_ID) = gamma_photon/(FullWidthBaryon(res_ID,mass_res)+gamma_photon) ! branching ratio into N+gamma
        sigma2 (res_ID) = sigma1 (res_ID) * BR (res_ID)
      end if
    end do
    write (80,'(33G12.5)') E_gamma,mass_res,sum(sigma1),sigma1 ! output to file
    write (81,'(33G12.5)') E_gamma,mass_res,sum(sigma2),sigma2
    write (82,'(33G12.5)') E_gamma,mass_res,sum(BR),BR
  end do

  close (80)
  close (81)
  close (82)

end program
