!******************************************************************************
!****m* /relativeVelocity
! NAME
! module relativeVelocity
! PURPOSE
! Evaluates relative velocity of particle and a fermi sea of nucleons
!******************************************************************************
module relativeVelocity
  implicit none
  private
  real, save :: p_delta,e_delta !Globals to be used to pass additional input to vrel_integrand

  public ::  vrel_times_rho

contains

  !****************************************************************************
  !****f* relativeVelocity/vrel_times_rho
  ! NAME
  ! real function vrel_times_rho(rho,mass,p,fastIN)
  ! PURPOSE
  ! * Evaluates the average relative velocity of a particle with a fermi sea of protons or neutrons.
  ! * Returns  (average velocity x proton/neutron density)
  ! * See appendix A.4 of oliver's thesis for details
  ! INPUTS
  ! real, intent(in) :: rho              ! Density in fm^-3
  ! real, intent(in) :: mass             ! Mass of particle
  ! real, intent(in),dimension(1:3) :: p ! Momentum of particle
  ! OUTPUT
  ! average velocity x density [GeV^3]
  !****************************************************************************
  real function vrel_times_rho(rho,mass,p,fastIN)
    use vector, only: absVec
    use tabulation
    implicit none
    real, intent(in) :: rho              ! Density in fm^-3
    real, intent(in) :: mass             ! Mass of particle
    real, intent(in),dimension(1:3) :: p ! Momentum of particle
    logical, optional :: fastIN
    logical, parameter :: tabulate=.false. ! the tabulation is slower than the real calculation :(
    real(4), dimension(1:3) :: x,min, max, delta
    type(table_3dim),save :: table
    integer :: outOfBounds
    logical, save :: firstTime=.true.
    if (tabulate) then
       if (firstTime) then
          min=(/0.,1.,0./)
          max=(/0.1,2.,2./)
          delta=(/0.01,0.02,0.02/)
          write(*,*) ' Initializing vrel_times_rho...'
          table=init_table_3dim(min,max,delta,func)
          write(*,*) ' ... done.'
          firstTime=.false.
       end if
       x=(/rho,mass,absVec(p)/)
       vrel_times_rho= getValue(x,table,outOfBounds)
       if (outOfBounds.ne.0)  vrel_times_rho=calculate_vrel_times_rho(rho,mass,absVec(p),fastIN)
    else
       vrel_times_rho=calculate_vrel_times_rho(rho,mass,absVec(p),fastIN)
    end if
  end function vrel_times_rho


  real function func(a,b,c)
    implicit none
    real,intent(in) :: a,b,c
    func = calculate_vrel_times_rho (a,b,c)
  end function func


  real function calculate_vrel_times_rho(rho,mass,p,fastIN)
    use constants, only: pi
    use densityModule, only: fermiMomentum_noIsospin
    use quadpack, only: qag
    implicit none

    real, intent(in) :: rho              ! Density in fm^-3
    real, intent(in) :: mass             ! Mass of particle
    real, intent(in) :: p                ! Absolute Momentum of particle
    logical, optional :: fastIN

    logical,parameter :: fastDefault=.true. ! Switch for integration method
    logical :: fast
    real :: pf
    real :: p_nuc, delta_p_nuc
    integer,parameter :: numSteps=100000
    integer :: i
    real :: integral
    real :: abserr,epsAbs,epsRel
    integer :: neval, ier

    if (present(fastIN)) then
       fast=fastIN
    else
       fast=fastDefault
    end if

    p_delta=p

    if (p_delta.lt.1E-20) then
       ! If the particle rests, then the average relative velocity is 0.
       calculate_vrel_times_rho=0.
       return
    end if

    E_delta=sqrt(mass**2+p_delta**2)

    ! Fermi momentum:
    pf=fermiMomentum_noIsospin(rho)
    if (fast) then
       ! Use Quadpack
       epsAbs=0.05*p_delta/e_delta !=5%*velocity of the particle
       epsRel=0.05 ! 5% relative error
       call qag (vrel_integrand, 0., pf,epsAbs, epsRel, 1, integral, abserr, neval, ier )
       if (ier.ne.0) then
          write(*,*) 'Problem with qag in vrel:', abserr, neval, ier
       end if
    else
       ! Simple Riemann integral (to check the other method!)
       delta_p_nuc=pf/float(numSteps)
       integral=0.
       do i=1,numSteps
          p_nuc=float(i)*delta_p_nuc
          integral=integral+vrel_integrand(p_nuc)
       end do
       integral=integral*delta_p_nuc
    end if
    calculate_vrel_times_rho=(2./((2.*pi)**3)*2./3.*pi*E_delta/p_delta)*integral
  end function calculate_vrel_times_rho

  !****************************************************************************
  !****f* relativeVelocity/vrel_integrand
  ! NAME
  ! real function vrel_integrand(p_nuc)
  ! PURPOSE
  ! * Integrand for vrel_times_rho
  !****************************************************************************
  real function vrel_integrand(p_nuc)
    use constants, only: mN
    implicit none
    real, intent(in) :: p_nuc
    real :: e_nuc
    real :: a,b,c
    e_nuc=sqrt(mN**2+p_nuc**2)

    a=p_delta**2/e_delta**2
    b=p_nuc**2/e_nuc**2
    c=2.*p_delta*p_nuc/e_delta/e_nuc

    vrel_integrand=p_nuc*e_nuc*((max(a+b+c,1.E-12))**(3./2.)-(max(a+b-c,1.E-12))**(3./2.))

  end function vrel_integrand

end module relativeVelocity
