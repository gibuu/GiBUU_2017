program main
  ! Testrun for skyrme potential
  ! * can be used to extract new parameter sets!!
  use particleProperties, only: initParticleProperties
  use skyrme, only: evaluate_skyrmeParameters, lambdaRoot, check
  use output, only: line
  implicit none

  real, save :: U_0  =-0.13
  real, save :: p_0  = 0.8             !
  real, parameter :: E_bin=-0.016           ! GeV
  real, parameter :: rho_0=0.16*0.197**3   ! GeV^3
  real, parameter :: compressibility= 0.220 ! GeV
  integer :: i,j
  real :: lambda,toSolve,A,B,C,tau
  logical :: success
  !call evaluateParameter(U_0,p_0,E_bind,K,rho_0,A,B,tau,C,lambda)

  call initParticleProperties


  open(100,file="lambda.dat")
  write(100,*) "#lambda[GeV], lambda [1/fm], f(lambda)"
  do i=1,1000
     lambda=float(i)*0.001
     tosolve=lambdaRoot(lambda,rho_0,p_0, u_0,E_bin, compressibility)
     write(*,*) lambda, toSolve
     write(100,*) lambda, lambda/0.197,toSolve
  end do
  close(100)

  write(*,*)
  write(*,*)
  write(*,*) line
  write(*,*) line
  write(*,*) ' Getting new parameters:'
  write(*,*) line
  write(*,*) line


  open(33,file='parameterSpace.dat')
  write(33,'(7A15,A18)') 'p_0',' U_0',' A',' B','C','tau', 'lambda', 'success (I) und (L)'
  do i=0,500
     p_0=0.3+float(i)*0.002
     do j=0,75
        U_0=-0.15+float(j)*0.002
        call evaluate_skyrmeParameters(rho_0,p_0, u_0,E_bin, compressibility, A,B,C,tau,lambda,success)
        write(*,'(7E15.4,I3)') p_0, U_0, A, B,C,tau, lambda, success
        write(33,'(7E15.4,I10,L8)') p_0, U_0, A, B,C,tau, lambda, success, success
     end do
     write(33,*)
  end do
  close(33)


contains

  subroutine checkOld()

    real, parameter :: U_0  =-0.075
    real, parameter :: p_0  = 0.8             !
    real, parameter :: E_bin=-0.016           ! GeV
    real, parameter :: rho_0=0.168*0.197**3   ! GeV^3
    real, parameter :: compressibility= 0.290 ! GeV

    write(*,*)
    write(*,*)
    write(*,*) line
    write(*,*) line
    write(*,*) ' Check old parameters (EQS=5)'
    write(*,*) line
    write(*,*) line
    call check(rho_0, p_0, E_bin, compressibility,-0.029253,0.057248,-0.063516, 1.76, 2.13*0.197)

    write(*,*)
    write(*,*)
    write(*,*) line
    write(*,*) line
    write(*,*) ' Check old parameters (EQS=3)'
    write(*,*) line
    write(*,*) line
    call check(rho_0, p_0, E_bin, 0.215,-0.28699993, 0.2336517, 0.,    1.225707, 2.13*0.197)

  end subroutine checkOld


end program main
