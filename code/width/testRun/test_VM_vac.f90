module deriv

  implicit none

  integer :: ID

contains

  real function gamma_(m)
    use mesonWidth
    real, intent(in) :: m
    gamma_ = m * fullWidthMeson(ID,m)
  end function

end module


program test_VM_vac
  use inputGeneral, only: readinputGeneral
  use particleProperties, only: InitParticleProperties, hadron
  use mesonWidth, only: partialwidthMeson, fullWidthMeson
  use mesonWidthVacuum, only: dileptonWidth
  use IdTable
  use derivatives, only: derivative
  use deriv, only: ID, gamma_
  implicit none

  real, save :: Gamma_coll_rho   = 0.150 ! in GeV
  real, save :: Gamma_coll_omega = 0.150
  real, save :: Gamma_coll_phi   = 0.030

  integer :: i!, ID
  real, parameter :: dm = 0.001
  real :: mass, width(0:4), M0, Y, C, C_med

  write(*,*) "**********************************"
  write(*,*) "**********************************"
  write(*,*)
  Write(*,*)      "Testing the routines which generate the vacuum width of the vector mesons!"
  write(*,*)
  write(*,*) "**********************************"
  write(*,*) "**********************************"
  write(*,*)

  call readinputGeneral
  call InitParticleProperties

  ID = rho
  M0 = hadron(ID)%mass
  Write(*,*)      "Testing Vacuum Decay Width of Meson #",ID
  do i=0,5000
     mass=i*dm
     width(0)=fullWidthMeson(ID,mass)
     width(1)=partialwidthMeson(ID,mass,pion,pion) ! rho -> 2pi
     width(2)=dileptonWidth(ID,mass)
     Y = ( mass**2 - M0**2 ) / (mass*width(0)) ! off-shell parameter
     C = Y/(2*mass) * derivative(gamma_,mass,dm,1,0)
     C_med = (mass**2-M0**2)/(2*mass**2) * Gamma_coll_rho/(width(0)+Gamma_coll_rho)
     write(ID,'(6F15.9)') mass,width(0:2),C,C_med
  end do

  ID = omegaMeson
  M0 = hadron(ID)%mass
  Write(*,*)      "Testing Vacuum Decay Width of Meson #",ID
  do i=0,5000
     mass=i*dm
     width(0)=fullWidthMeson(ID,mass)
     width(1)=partialwidthMeson(ID,mass,pion,pion,pion,+1,-1,0) ! omega -> 3pi
     width(2)=partialwidthMeson(ID,mass,pion,photon) ! omega -> pi0 gamma
     width(3)=partialwidthMeson(ID,mass,pion,pion)   ! omega -> 2pi
     width(4)=dileptonWidth(ID,mass)
     Y = ( mass**2 - M0**2 ) / (mass*width(0)) ! off-shell parameter
     C = Y/(2*mass) * derivative(gamma_,mass,dm,1,0)
     C_med = (mass**2-M0**2)/(2*mass**2) * Gamma_coll_omega/(width(0)+Gamma_coll_omega)
     write(ID,'(8F15.9)') mass,width(0:4),C, C_med
  end do

  ID = phi
  M0 = hadron(ID)%mass
  Write(*,*)      "Testing Vacuum Decay Width of Meson #",ID
  do i=0,5000
     mass=i*dm
     width(0)=fullWidthMeson(ID,mass)
     width(1)=partialwidthMeson(ID,mass,kaon,kaonBar) ! phi -> 2K
     width(2)=partialwidthMeson(ID,mass,rho,pion)     ! phi -> rho pi0
     width(3)=partialwidthMeson(ID,mass,pion,pion,pion,+1,-1,0) ! phi -> 3pi
     width(4)=dileptonWidth(ID,mass)
     Y = ( mass**2 - M0**2 ) / (mass*width(0)) ! off-shell parameter
     C = Y/(2*mass) * derivative(gamma_,mass,dm,1,0)
     C_med = (mass**2-M0**2)/(2*mass**2) * Gamma_coll_phi/(width(0)+Gamma_coll_phi)
     write(ID,'(8F15.9)') mass,width(0:4),C,C_med
  end do

end program test_VM_vac
