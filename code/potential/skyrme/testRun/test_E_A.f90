
program main
  ! Testrun for skyrme potential
  ! * can be used to extract new parameter sets!!
  use skyrme, only: evaluate_skyrmeParameters, E_div_A
  use particleProperties, only: initParticleProperties
  implicit none

  real, save :: U_0  =-0.13
  real, save :: p_0  = 0.8             !
  real, parameter :: E_bin=-0.016           ! GeV
  real, parameter :: rho_0=0.16*0.197**3   ! GeV^3
  real, parameter :: compressibility= 0.220 ! GeV
  integer :: i,j
  real :: lambda, rho,A,B,C,tau,eA
  logical :: success
  !call evaluateParameter(U_0,p_0,E_bind,K,rho_0,A,B,tau,C,lambda)

  call initParticleProperties
  
  
  do j=1,3
     Select case(j)
     case(1)
        U_0=-0.09
        p_0=0.8
     case(2)
        U_0=-0.1
        p_0=0.6
     case(3)
        U_0=-0.11
        p_0=0.45
     end select
     call evaluate_skyrmeParameters(rho_0,p_0, u_0,E_bin, compressibility, A,B,C,tau,lambda,success)

     do i=0,500
        rho=0.001*0.197**3+float(i)*0.001*0.197**3
        eA= E_div_A(rho, rho_0, A,B,C, tau, lambda)
        write(10+j,*) rho,rho/0.197**3,eA
     end do
  end do


end program main
