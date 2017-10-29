program main
  ! Test the 5-fold differential cross section by integrating over dOmega(pion). 
  ! The resulting dSigma/dOmega(electron)/dE(electron) can be compared to the 
  ! experimental results of O'Connell et al. (PRL 53(17))
  ! Writes result to file 'dsigma_dOmegadE.dat'.
  call init_database

  call test_piN_new
contains


  subroutine test_piN_new
  ! Test the 5-fold differential cross section by integrating over dOmega(pion). 
  ! The resulting dSigma/dOmega(electron)/dE(electron) can be compared to the 
  ! experimental results of O'Connell et al. (PRL 53(17))
  ! Writes result to file 'dsigma_dOmegadE.dat'.

  use constants, only :pi
  use electronNuc_to_Delta, only : elepi

  implicit none

  ! The O'Connell experiment: PRL 53(17)
  real,parameter :: initial_energy=0.73
  real,parameter :: thetaElec_rad=37.1/180.*pi

  integer, parameter :: maxSteps=100               ! Number of photon energies
  integer, parameter :: maxSteps_phi=50          ! Number of integration steps in phi(pion)
  integer, parameter :: maxSteps_dCosTheta=50   ! Number of integration steps in cosTheta(pion)

  real     :: final_energy
  real     :: phiPion_rad,thetaPion_rad
  real     :: cs,total
  integer  :: i,j,k
  real     :: dcosTheta, dPhi

  ! Step sizes for dOmega(pion) integration:
  dcosTheta=2./float(maxsteps_dCosTheta)
  dPhi=2.*pi/float(maxsteps_phi)

  write(*,*) 'dPhi=',dphi
  write(*,*) 'dCos(Theta)=',dCosTheta
  write(*,*)
  open(11,file='dsigmadOmegadE.dat')

  thetaPion_rad=10./180.*pi
  phiPion_rad=60./180.*pi
  final_energy=0.3
  
  ! Sum both proton contributions:
!  call elepi(1,1,initial_energy,final_energy,thetaElec_rad,thetaPion_rad,phiPion_rad,cs)
!  stop

  
  energyLoop: do i=1,maxsteps-1
     final_Energy=initial_energy/float(maxsteps)*float(i)
     write(*,*) i,'/',maxsteps-1,'q_0=',-final_Energy+initial_energy

     ! Integrate over phi and cos(theta) of the outgoing pion
     total=0.
     phiLoop: do j=1,maxsteps_phi
        phiPion_rad=(float(j)-0.5)*dphi
        cosThetaLoop: do k=1,maxsteps_dCostheta
           thetaPion_rad=acos(-1.+(float(k)-0.5)*dcosTheta)
           ! Sum both proton contributions:
           call elepi(1,initial_energy,final_energy,thetaElec_rad,thetaPion_rad,phiPion_rad,cs)
           total=total+cs
           call elepi(2,initial_energy,final_energy,thetaElec_rad,thetaPion_rad,phiPion_rad,cs)
           total=total+cs
        end do cosThetaLoop
     end do phiLoop
     total=total*dphi*dcosTheta
     write(11,'(6E15.3)')initial_energy-final_energy,  total
  end do energyLoop
  close(11)

end subroutine test_piN_new










end program main
