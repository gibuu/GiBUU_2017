!******************************************************************************
!****m* /skyrme
! NAME
! module skyrme
!
! PURPOSE
! * This module includes various routines for the Skyrme mean field potential
!
! NOTES
! * For documentation see appendix A.4 of Oliver's thesis
!
!******************************************************************************
module skyrme
  implicit none

  private

  public :: evaluate_skyrmeParameters, U

  public :: E_div_A, lambdaRoot, check

contains

  !****************************************************************************
  !****s* skyrme/evaluate_skyrmeParameters
  ! NAME
  ! subroutine evaluate_skyrmeParameters (rhoNull,pNull, uNull,bindingEnergy, compressibility,A,B,C,tau,lambda,success)
  ! PURPOSE
  ! Evaluates the Skyrme parameters A,B,C,tau,lambda for given input
  ! * rhoNull -- saturation density in GeV^3
  ! * pNull   -- U(rho=rhoNull,p=p_0)=0  -> p_0 is root of potential in GeV
  ! * uNull   -- U(rho=rhoNull,p=0  )=uNull  -> depth of potential at p=0 in GeV
  ! * bindingEnergy   -- Nuclear matter binding energy in GeV
  ! * compressibility -- Nuclear matter compressibility in GeV
  ! OUTPUT
  ! *  real, intent(out)     :: A,B,C,lambda ! Potential parameters, units: GeV
  ! *  real, intent(out)     :: tau          ! Potential parameter , no units
  ! *  logical, intent(out)  :: success      ! true if parameters could be determined
  ! NOTES
  ! See appendix A.4 of Oliver's thesis for details
  !****************************************************************************
  subroutine evaluate_skyrmeParameters (rhoNull,pNull, uNull,bindingEnergy, compressibility,A,B,C,tau,lambda,success)
    use findZero, only: bisection_findZero
    use constants, only: hbarc

    real, intent(in)  :: rhoNull,pNull, uNull,bindingEnergy, compressibility
    real, intent(out) :: A,B,C,tau,lambda
    logical,intent(out) :: success

    !integer :: i
    real :: dummy,lambdaPole,lambda1,lambda2
    logical :: success1,success2


    ! Initializing the variables with some dummy value for lambda
    dummy=lambdaRoot (3.,rhoNull,pNull, uNull,bindingEnergy, compressibility)
    dummy=lambdaRootPole (3.,rhoNull,pNull, uNull,bindingEnergy)
    ! Searching the root of the function:
    ! Assume lambda > pole of lambdaRoot(x)
    ! 1) Search pole
    lambdaPole=bisection_findZero(lambdaRootPole_dummy,0.0001,5.,0.0001,100,success)
    if (success) then
       ! 2) Search lambda before pole and after pole
       ! Before:
       write(*,*) 'Searching root before the pole. lambda_pole=', lambdaPole
       lambda1=bisection_findZero(lambdaRoot_dummy,0.001,lambdaPole-0.001,0.0001,100,success1)
       ! After:
       write(*,*) 'Searching root after the pole. lambda_pole=', lambdaPole
       lambda2=bisection_findZero(lambdaRoot_dummy,lambdaPole+0.001,30.*hbarc,0.0001,100,success2)
       if (success1.and.success2) then
          write(*,*) 'TWO solutions!'
          write(*,*) 'Lambda',lambda1,lambda2
          success=.false.
          stop
       else if (success1) then
          lambda=lambda1
          success=.true.
       else if (success2) then
          lambda=lambda2
          success=.true.
       else
          success=.false.
       end if
    else
       !write(*,*) "Can't find pole"
       lambda=bisection_findZero(lambdaRoot_dummy,0.001,30.*hbarc,0.0001,1000,success)
    end if


    if (success) then
       call getParameters(rhoNull,pNull, uNull, compressibility,A,B,C,tau,lambda)
       call check(rhoNull, pNull, bindingEnergy, compressibility,A,B,C, tau, lambda)
    else
       A=0.
       B=0.
       C=0.
       tau=0.
       lambda=0.
    end if
  end subroutine evaluate_skyrmeParameters




  !****************************************************************************
  !****f* skyrme/U
  ! NAME
  ! real function U(rho, p, rhoNull, A,B,C, tau, lambda)
  ! PURPOSE
  ! Returns the single-particle potential in units of GeV.
  ! INPUTS
  ! * rhoNull -- saturation density in GeV^3
  ! * rho     -- baryon density in GeV^3
  ! * p       -- momentum in GeV
  ! * A,B,C,lambda -- potential parameters in GeV
  ! * tau -- potential parameters (no units)
  ! NOTES
  ! See appendix A.4 of Oliver's thesis for details
  !****************************************************************************
  real function U(rho, p, rhoNull, A,B,C, tau, lambda)
    real, intent(in) :: rho, p, rhoNull, A,B,C, tau, lambda
    U=A*rho/rhoNull+B*(rho/rhoNull)**tau+2.*C/rhoNull*f2(lambda,rho,p)
  end function U



  !****************************************************************************
  !****f* skyrme/E_div_A
  ! NAME
  ! real function E_div_A(rho, rhoNull, A,B,C, tau, lambda)
  ! PURPOSE
  ! Returns the nuclear binding energy in GeV
  ! INPUTS
  ! * rhoNull -- saturation density in GeV^3
  ! * rho     -- baryon density in GeV^3
  ! * A,B,C,lambda -- potential parameters in GeV
  ! * tau -- potential parameters (no units)
  ! NOTES
  ! See appendix A.4 of Oliver's thesis for details
  !****************************************************************************
  real function E_div_A(rho, rhoNull, A,B,C, tau, lambda)
    use constants, only: mN
    real, intent(in) :: rho, rhoNull, A,B,C, tau, lambda
    E_div_A = 1./rho*(ekin(rho)+A/2.*rho**2/rhoNull+B/(tau+1.)*rho**(tau+1.)/rhoNull**tau+C/rhoNull*f1(lambda,rho)) &
              - mN
  end function E_DIV_A


  !****************************************************************************
  !****f* skyrme/dRho_E_div_A
  ! NAME
  ! real function dRho_E_div_A(rho, rhoNull, A,B,C, tau, lambda)
  ! PURPOSE
  ! Returns the rho-derivative of the nuclear binding energy
  ! INPUTS
  ! * rhoNull -- saturation density in GeV^3
  ! * rho     -- baryon density in GeV^3
  ! * A,B,C,lambda -- potential parameters in GeV
  ! * tau -- potential parameters (no units)
  ! NOTES
  ! See appendix A.4 of Oliver's thesis for details
  !****************************************************************************
  real function dRho_E_div_A(rho, rhoNull, A,B,C, tau, lambda)
    use derivatives
    use constants, only: hbarc

    real, intent(in) :: rho, rhoNull, A,B,C, tau, lambda

    real :: delta_rho=0.001*hbarc**3
    integer :: i, scheme
    real, dimension (-2:2)  ::  deriv_field

    scheme=0
    do i=-2,2
       if (rho+float(i)*delta_rho.gt.0) then
          deriv_field(i)=E_div_A(rho+float(i)*delta_rho, rhoNull, A,B,C, tau, lambda)
       else
          deriv_field(i)=0.
          scheme=999 ! -> Must use forward differencing
       end if
    end do

    dRho_E_div_A=finiteDifference(deriv_field,delta_rho,1,scheme)

  end function DRHO_E_DIV_A


  !****************************************************************************
  !****f* skyrme/dRhoSquared_E_div_A
  ! NAME
  ! real function dRhoSquared_E_div_A(rho, rhoNull, A,B,C, tau, lambda)
  ! PURPOSE
  ! Returns the second rho-derivative of the nuclear binding energy
  ! INPUTS
  ! * rhoNull -- saturation density in GeV^3
  ! * rho     -- baryon density in GeV^3
  ! * A,B,C,lambda -- potential parameters in GeV
  ! * tau -- potential parameters (no units)
  ! NOTES
  ! See appendix A.4 of Oliver's thesis for details
  !****************************************************************************
  real function dRhoSquared_E_div_A(rho, rhoNull, A,B,C, tau, lambda)
    use derivatives
    use constants, only: hbarc

    real, intent(in) :: rho, rhoNull, A,B,C, tau, lambda

    real :: delta_rho=0.001*hbarc**3
    integer :: i, scheme
    real, dimension (-2:2)  ::  deriv_field

    scheme=0
    do i=-2,2
       if (rho+float(i)*delta_rho.gt.0) then
          deriv_field(i)=dRho_E_div_A(rho+float(i)*delta_rho, rhoNull, A,B,C, tau, lambda)
       else
          deriv_field(i)=0.
          scheme=999 ! -> Must use forward differencing
       end if
    end do

    dRhoSquared_E_div_A=finiteDifference(deriv_field,delta_rho,1,scheme)


  end function DRHOSQUARED_E_DIV_A




  !****************************************************************************
  !****f* skyrme/pf
  ! NAME
  ! real function pf(rho)
  ! PURPOSE
  ! Returns the Fermi momentum in GeV
  ! INPUTS
  ! * rho     -- baryon density in GeV^3
  !****************************************************************************
  real function pf(rho)
    use constants, only: pi
    real, intent (in) :: rho ! Baryon density in units GeV^3
    pf=(3./2.*pi**2*rho)**(1./3.)
  end function pf



  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  ! f1
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

  ! f1= Integral d3p' d3p f(r,p) f(r,p') /(1+(p-p')**2/Lambda**2)
  real function f1(lambda,rho)
    use constants, only: pi

    real, intent(in) :: lambda
    real, intent(in) :: rho       ! Units GeV^3

    real :: pfermi
    pfermi=pf(rho)

    ! See Welke PRC 38(5), 2101, Formula 5.6
    f1=32.*(pi**2)/3.*pfermi**4*lambda**2* (  &
         &   3./8.-lambda/(2.*pfermi)*atan(2.*pfermi/lambda)-lambda**2/(16.*pfermi**2) &
         & + (3./16.*lambda**2/pfermi**2+1./64.*lambda**4/pfermi**4)*log(1.+4.*pfermi**2/lambda**2) &
         &     ) !/(2.*pi)**6
    f1=f1*16./(2.*pi)**6    !f1=g^2 * int d^3p/(2pi)^3 ...
  end function f1

  ! f1= Integral d3p' d3p f(r,p) f(r,p') /(1+(p-p')**2/Lambda**2)
  ! Integration by Riemann sums
  ! * Used to check f1 from above.
!   real function integral_f1(lambda,rho)
!     use constants, only : pi
!
!     real, intent(in) :: lambda
!     real, intent(in) :: rho       ! Units GeV^3
!
!     real,parameter :: dp         =0.001
!     real,parameter :: dp_prime   =0.001
!     real,parameter :: dcos_theta =0.001
!
!     real :: pfermi, p, p_prime, cos_theta
!     integer :: i,j,k
!     integer:: N_steps_p_prime, N_steps_p, N_steps_cos_theta
!
!     pfermi=pf(rho)
!
!
!     N_steps_p        =Nint( pfermi/dp       )
!     N_steps_p_prime  =Nint( pfermi/dp_prime )
!     N_steps_cos_theta=Nint( 2./dcos_theta   )
!
!
!     integral_f1=0.
!
!     do i=0,N_steps_p-1
!        p=(float(i)+0.5)*dp
!        do j=0,N_steps_p_prime-1
!           p_prime=(float(j)+0.5)*dp_prime
!           do k=0,N_steps_cos_theta-1
!              cos_theta=-1.+(float(k)+0.5)*dcos_theta
!              integral_f1=integral_f1+p**2*p_prime**2 /(1.+(p**2+p_prime**2-2.*p*p_prime*cos_theta)/Lambda**2)
!           end do
!        end do
!     end do
!     integral_f1=   integral_f1   * dp * dp_prime * dcos_theta
!     integral_f1=   integral_f1   * 4.*pi*2.*pi
!   end function integral_f1



  ! df_1/drho
  real function df1_drho(lambda,rho)
    use derivatives
    use constants, only: hbarc

    real, intent(in) :: lambda
    real, intent(in) :: rho       ! Units GeV^3
    real :: delta_rho=0.001*hbarc**3
    integer :: i, scheme
    real, dimension (-2:2)  ::  f1_field

    scheme=0
    do i=-2,2
       if (rho+float(i)*delta_rho.gt.0) then
          f1_field(i)=f1(lambda,rho+float(i)*delta_rho)
       else
          f1_field(i)=0.
          scheme=999 ! -> Must use forward differencing
       end if
    end do

    df1_drho=finiteDifference(f1_field,delta_rho,2,scheme)

  end function df1_drho




  ! d^2f_1/drho^2
  real function df1_drhoSquared(lambda,rho)
    use derivatives
    use constants, only: hbarc

    real, intent(in) :: lambda

    real, intent(in) :: rho       ! Units GeV^3
    real :: delta_rho=0.001*hbarc**3
    integer :: i, scheme
    real, dimension (-2:2)  ::  f1_field

    scheme=0
    do i=-2,2
       if (rho+float(i)*delta_rho.gt.0) then
          f1_field(i)=df1_drho(lambda,rho+float(i)*delta_rho)
       else
          f1_field(i)=0.
          scheme=999 ! -> Must use forward differencing
       end if
    end do

    df1_drhoSquared=finiteDifference(f1_field,delta_rho,2,scheme)


  end function df1_drhoSquared



  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  ! f2
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************


  ! f2= Integral d3p' f(r,p') /(1+(p-p')**2/Lambda**2)
  real function f2(lambda,rho,p)
    use constants, only: pi

    real, intent(in) :: lambda
    real, intent(in) :: p
    real, intent(in) :: rho       ! Units GeV^3

    real :: pfermi, logApprox,y

    pfermi=pf(rho)

    if (abs(2.*pfermi*p/(pfermi**2+lambda**2)).lt.1E-6) then
       ! Make Taylor expansion for the logarithm to find the x-> 0 limit of the first term:
       ! We define
       y=2.*pfermi/(pfermi**2+p**2+lambda**2)
       ! and realize that this term becomes small:
       ! y'=p*y .
       ! So we can make a Taylor expansion in y':
       ! log((1+y')/(1-y'))=2y'+2y'^3/3+2y'^5/5+...
       ! Thus the p-> 0 limit of the first term is given
       ! 1/p log( ((p+pfermi)**2+lambda**2)/((p-pfermi)**2+lambda**2)    )  for small y =
       ! = 2/p (y'+y'^3/3+y'^5/5+..)
       logApprox=2.*(y+y**3*p**2+y**5*p**4)

       f2=pi*lambda**3*(&
            & (pfermi**2+lambda**2-p**2)/(2.*lambda)*logApprox &
            & +2.*pfermi/lambda-2.*(atan((p+pfermi)/lambda)-atan((p-pfermi)/lambda)) &
            & )
    else
       ! See Welke PRC 38(5), 2101, Formula 5.6
       f2=pi*lambda**3*(&
            & (pfermi**2+lambda**2-p**2)/(2.*p*lambda)*log( ((p+pfermi)**2+lambda**2)/((p-pfermi)**2+lambda**2)    ) &
            & +2.*pfermi/lambda-2.*(atan((p+pfermi)/lambda)-atan((p-pfermi)/lambda)) &
            & )
    end if
    f2=f2*4./(2.*pi)**3

  end function f2



  ! f2= Integral d3p' f(r,p') /(1+(p-p')**2/Lambda**2)
  ! Integration by Riemann sums
  ! * Used to check f2 from above.
!   real function integral_f2(lambda,rho,p)
!     use constants, only : pi
!
!     real, intent(in) :: lambda
!     real, intent(in) :: p
!     real, intent(in) :: rho       ! Units GeV^3
!
!
!     real,parameter :: dp_prime   =0.001
!     real,parameter :: dcos_theta =0.01
!
!     real :: pfermi, p_prime, cos_theta
!     integer :: j,k
!     integer:: N_steps_p_prime, N_steps_cos_theta
!
!     pfermi=pf(rho)
!     N_steps_p_prime  =Nint( pfermi/dp_prime )
!     N_steps_cos_theta=Nint( 2./dcos_theta   )
!
!     integral_f2=0.
!
!     do j=0,N_steps_p_prime-1
!        p_prime=(float(j)+0.5)*dp_prime
!        do k=0,N_steps_cos_theta-1
!           cos_theta=-1.+(float(k)+0.5)*dcos_theta
!           integral_f2=integral_f2+p_prime**2 /(1.+(p**2+p_prime**2-2.*p*p_prime*cos_theta)/Lambda**2)
!        end do
!     end do
!
!     integral_f2=   integral_f2   * dp_prime * dcos_theta
!     integral_f2=   integral_f2   * 2.*pi
!
!     integral_f2=integral_f2*4./(2.*pi)**3
!
!   end function integral_f2


  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  ! f3
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************


  ! f3= 2*f1/rhoNull**2 -2/rhoNull*d/drho (f1) +d^2/drho^2 (f1)
  real function f3(lambda,rho,rhoNull)

    real, intent(in) :: lambda
    real, intent(in) :: rho       ! Units GeV^3
    real, intent(in) :: rhoNull   ! Units GeV^3

    ! See Welke PRC 38(5), 2101, Formula 5.6
    f3= 2./rhoNull**2 * f1(lambda,rho)-2./rhoNull *df1_drho(lambda,rho) +df1_drhoSquared(lambda,rho)

  end function f3


  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  ! f4
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************



  ! f4= d/drho (f1) - f1/rhoNull
  real function f4(lambda,rho,rhoNull)

    real, intent(in) :: lambda
    real, intent(in) :: rho       ! Units GeV^3
    real, intent(in) :: rhoNull   ! Units GeV^3

    ! See Welke PRC 38(5), 2101, Formula 5.6
    f4= df1_drho(lambda,rho) -f1(lambda,rho)/rhoNull

  end function f4



  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  ! f5
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************



  ! f5= f1/rhoNull**2-2*f2(lambda,p=p0)/rhoNull
  real function f5(lambda,rho,rhoNull,pNull)

    real, intent(in) :: lambda
    real, intent(in) :: rho       ! Units GeV^3
    real, intent(in) :: rhoNull   ! Units GeV^3
    real, intent(in) :: pNull   ! Units GeV^3

    f5=f1(lambda,rho)/rhoNull**2-f2(lambda,rho,pNull)/rhoNull

  end function f5

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  ! f6
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************


  ! f6= f4/rhoNull**2-2*f2(lambda,p=p0)/rhoNull
  real function f6(lambda,rho,rhoNull,pNull)

    real, intent(in) :: lambda
    real, intent(in) :: rho       ! Units GeV^3
    real, intent(in) :: rhoNull   ! Units GeV^3
    real, intent(in) :: pNull   ! Units GeV^3

    f6=f4(lambda,rho,rhoNull)/rhoNull**2-f2(lambda,rho,pNull)/rhoNull**2

  end function f6


  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  ! f7
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

  ! f7=2./rhoNull*( f2(p=0.) -  f2(p=pNull))
  real function f7(lambda,rho,rhoNull,pNull)

    real, intent(in) :: lambda
    real, intent(in) :: rho       ! Units GeV^3
    real, intent(in) :: rhoNull   ! Units GeV^3
    real, intent(in) :: pNull   ! Units GeV^3

    !write(11,*) rhoNull
    !write(11,*) rhoNull,f2(lambda,rho,0.)
    !write(11,*) rhoNull,f2(lambda,rho,0.),f2(lambda,rho,pNull)
    f7=2./rhoNull*( f2(lambda,rho,0.)-  f2(lambda,rho,pNull))

  end function f7


  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  ! The constants
  ! See A.4 of Oliver's thesis for details
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

  real function c1(rhoNull,bindingEnergy)
    use constants, only: mN
    real, intent(in) :: rhoNull , bindingEnergy
    c1=bindingEnergy-ekin(rhoNull)/rhoNull+mN
  end function c1

  real function c2(compressibility,rhoNull)
    real, intent(in) :: compressibility, rhoNull
    c2=compressibility/9./rhoNull**2-d_drhoSquared_ekin_div_rho(rhoNull)
  end function c2

  real function c3(rhoNull)
    real, intent(in) :: rhoNull
    c3=ekin(rhoNull)/rhoNull**2-1./rhoNull*d_drho_ekin(rhoNull)
  end function c3



  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  ! The kinetic energy contribution
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************


  !****************************************************************************
  !****s* skyrme/ekin
  ! NAME
  ! real function ekin(rho)
  ! PURPOSE
  ! Evaluates the kinetic energy density
  ! INPUTS
  ! real, intent(in) :: rho ! density in GeV^3
  ! NOTES
  ! * See appendix A.4 of Oliver's thesis for details
  !****************************************************************************
  real function ekin(rho)
    use constants, only: pi, mN

    real, intent(in) :: rho
    real :: p_f

    p_f=pf(rho)

    ekin=4.* 4.*pi*( &
         &  p_f/4.*(p_f**2+mN**2)**(3./2.)  &
         & - mN**2/8.*(p_f*sqrt(p_f**2+mN**2)+mN**2*log( (p_f+sqrt(p_f**2+mN**2))/mN ) )  &
         & ) /(2.*pi)**3

  end function ekin


  ! Only for checking ekin(rho)
!   real function ekin_integral(rho)
!     use constants, only: pi, mN
!
!     real, intent(in) :: rho
!
!     real,parameter :: dp         =0.001
!     real    :: pfermi, p
!     integer :: i
!     integer :: N_steps_p
!
!
!     pfermi=pf(rho)
!     N_steps_p        =Nint( pfermi/dp)
!
!     ekin_integral=0.
!
!     do i=0,N_steps_p-1
!        p=(float(i)+0.5)*dp
!        ekin_integral= ekin_integral+p**2*sqrt(p**2+mN**2)
!     end do
!
!     ekin_integral=ekin_integral*2.* 4.*pi*dp /(2.*pi)**3
!
!   end function ekin_integral


  !****************************************************************************
  !****s* skyrme/d_drho_ekin
  ! NAME
  ! real function d_drho_ekin(rho)
  ! PURPOSE
  ! Evaluates the rho-derivative of the kinetic energy density
  ! INPUTS
  ! real, intent(in) :: rho ! density in GeV^3
  ! NOTES
  ! * See appendix A.4 of Oliver's thesis for details
  !****************************************************************************
  real function d_drho_ekin(rho)
    use derivatives
    use constants, only: hbarc

    real, intent(in) :: rho
    real :: delta_rho=0.001*hbarc**3
    integer :: i, scheme
    real, dimension (-2:2)  ::  deriv_field

    scheme=0
    do i=-2,2
       if (rho+float(i)*delta_rho.gt.0) then
          deriv_field(i)=ekin(rho+float(i)*delta_rho)
       else
          deriv_field(i)=0.
          scheme=999 ! -> Must use forward differencing
       end if
    end do

    d_drho_ekin=finiteDifference(deriv_field,delta_rho,2,scheme)

  end function d_drho_ekin

  !****************************************************************************
  !****s* skyrme/d_drhoSquared_ekin
  ! NAME
  ! real function d_drhoSquared_ekin(rho)
  ! PURPOSE
  ! Evaluates the 2nd rho-derivative of the kinetic energy density
  ! INPUTS
  ! real, intent(in) :: rho ! density in GeV^3
  ! NOTES
  ! * See appendix A.4 of Oliver's thesis for details
  !****************************************************************************
  real function d_drhoSquared_ekin(rho)
    use derivatives
    use constants, only: hbarc

    real, intent(in) :: rho
    real :: delta_rho=0.001*hbarc**3
    integer :: i, scheme
    real, dimension (-2:2)  ::  deriv_field

    scheme=0
    do i=-2,2
       if (rho+float(i)*delta_rho.gt.0) then
          deriv_field(i)=d_drho_ekin(rho+float(i)*delta_rho)
       else
          deriv_field(i)=0.
          scheme=999 ! -> Must use forward differencing
       end if
    end do

    d_drhoSquared_ekin=finiteDifference(deriv_field,delta_rho,2,scheme)

  end function d_drhoSquared_ekin


  ! d/drho(Ekin/rho)
!   real function d_drho_ekin_div_rho(rho)
!     real, intent(in) :: rho
!     d_drho_ekin_div_rho=d_drho_ekin(rho)/rho-ekin(rho)/rho**2
!   end function d_drho_ekin_div_rho


  ! d^2/drho^2(Ekin/rho)
  real function d_drhoSquared_ekin_div_rho(rho)
    real, intent(in) :: rho
    d_drhoSquared_ekin_div_rho=1./rho*d_drhoSquared_ekin(rho)-2./rho**2*d_drho_ekin(rho)+2*ekin(rho)/rho**3
  end function d_drhoSquared_ekin_div_rho



  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  ! Solving for LAMBDA
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

  real function lambdaRoot(lambda,rhoNull_in,pNull_in, uNull_in,bindingEnergy_in, compressibility_in)

    real, intent(in) :: lambda
    real, intent(in),optional :: rhoNull_in,pNull_in, uNull_in,bindingEnergy_in, compressibility_in
    real, save  :: rhoNull_save,pNull_save, uNull_save,bindingEnergy_save, compressibility_save
    real, save  :: rhoNull,pNull, uNull,bindingEnergy, compressibility
    logical , save :: initialized=.false.

    if (present(rhoNull_in).and.present(pNull_in).and.present(uNull_in)&
         & .and.present(bindingEnergy_in).and.present(compressibility_in)) then
       initialized=.true.
       rhoNull_save=rhoNull_in
       pNull_save  =pNull_in
       uNull_save  =uNull_in
       bindingEnergy_save=bindingEnergy_in
       compressibility_save=compressibility_in
    else if (.not.initialized) then
       write(*,*) 'Variables not initialized in eq! STOP!'
       stop
    end if

    rhoNull=rhoNull_save
    pNull  =pNull_save
    uNull  =uNull_save
    bindingEnergy=bindingEnergy_save
    compressibility=compressibility_save

    lambdaRoot=1.+1/rhoNull/(c3(rhoNull)-uNull *f6(lambda,rhoNull,rhoNull,pNull)/f7(lambda,rhoNull,rhoNull,pNull)) &
         &  *(c1(rhoNull,bindingEnergy)-uNull*f5(lambda,rhoNull,rhoNull,pNull)/f7(lambda,rhoNull,rhoNull,pNull))


    !write(99,'(6E20.5)') lambda, lambdaRoot , c3(rhoNull)-uNull *f6(lambda,rhoNull,rhoNull,pNull)/f7(lambda,rhoNull,rhoNull,pNull), &
    !     &  (c1(rhoNull,bindingEnergy)-uNull*f5(lambda,rhoNull,rhoNull,pNull)/f7(lambda,rhoNull,rhoNull,pNull)) ,lambdaRootPole(lambda,rhoNull,pNull, uNull,bindingEnergy),lambdaRootPole(lambda)


  end function lambdaRoot


  real function lambdaRootPole(lambda,rhoNull_in,pNull_in, uNull_in,bindingEnergy_in)

    real, intent(in) :: lambda
    real, intent(in),optional :: rhoNull_in,pNull_in,uNull_in,bindingEnergy_in
    real, save  :: rhoNull_save,pNull_save, uNull_save,bindingEnergy_save
    real, save  :: rhoNull,pNull, bindingEnergy, uNull
    logical , save :: initialized=.false.

    if (present(rhoNull_in).and.present(pNull_in).and.present(uNull_in).and.present(bindingEnergy_in)) then
       initialized=.true.
       rhoNull_save=rhoNull_in
       pNull_save  =pNull_in
       uNull_save  =uNull_in
       bindingEnergy_save=bindingEnergy_in
    else if (.not.initialized) then
       write(*,*) 'Variables not initialized in eq! STOP!'
       stop
    end if

    rhoNull=rhoNull_save
    pNull  =pNull_save
    uNull  =uNull_save
    bindingEnergy=bindingEnergy_save
    lambdaRootPole=c3(rhoNull)-uNull *f6(lambda,rhoNull,rhoNull,pNull)/f7(lambda,rhoNull,rhoNull,pNull)
  end function lambdaRootPole


  ! This dummy function is needed for root searching in the lambda variable
  real function lambdaRoot_dummy(lambda)
    real, intent(in) :: lambda
    lambdaRoot_dummy=lambdaRoot(lambda)
  end function lambdaRoot_dummy


  ! This dummy function is needed for root searching in the lambda variable
  real function lambdaRootPole_dummy(lambda)
    real, intent(in) :: lambda
    lambdaRootPole_dummy=lambdaRootPole(lambda)
  end function lambdaRootPole_dummy


  !****************************************************************************
  !****s* skyrme/getParameters
  ! NAME
  ! subroutine getParameters(rhoNull,pNull, uNull, compressibility,A,B,C,tau,lambda)
  ! PURPOSE
  ! Evaluates A,B,C and tau if lambda and rhoNull, pNull, uNull, compressibility are given
  ! INPUTS
  ! real, intent(in) :: lambda
  ! real, intent(in)  :: rhoNull,pNull, uNull, compressibility
  ! OUTPUT
  ! real, intent(out) :: A,B,C,tau
  ! NOTES
  ! * All inputs and outputs in multiples of GeV
  ! * See appendix A.4 of Oliver's thesis for details
  !****************************************************************************
  subroutine getParameters(rhoNull,pNull, uNull, compressibility,A,B,C,tau,lambda)

    real, intent(in) :: lambda
    real, intent(in)  :: rhoNull,pNull, uNull, compressibility
    real, intent(out) :: A,B,C,tau

    tau=1./2./(c3(rhoNull)-uNull*f6(lambda,rhoNull,rhoNull,pNull)/f7(lambda,rhoNull,rhoNull,pNull)) &
         & *(c2(compressibility,rhoNull)-uNull/rhoNull**2*f3(lambda,rhoNull,rhoNull)/f7(lambda,rhoNull,rhoNull,pNull))*rhoNull

    B=rhoNull/(tau/(tau+1.)-0.5)*(c3(rhoNull)-uNull*f6(lambda,rhoNull,rhoNull,pNull)/f7(lambda,rhoNull,rhoNull,pNull))

    C=uNull/f7(lambda,rhoNull,rhoNull,pNull)

    A=-B-2.*C/rhoNull*f2(lambda,rhoNull,pNull)

  end subroutine getParameters



  !****************************************************************************
  !****s* skyrme/check
  ! NAME
  ! subroutine check(rhoNull, pNull, bindingEnergy, compressibility, A,B,C, tau, lambda)
  ! PURPOSE
  ! Checks whether  A,B,C, tau, lambda are valid parameters for given input
  ! rhoNull, pNull, bindingEnergy, compressibility
  ! INPUTS
  ! * real, intent(in) :: rhoNull, pNull, bindingEnergy, compressibility, A,B,C, tau, lambda
  ! * all units in multiples of GeV!
  ! NOTES
  ! See appendix A.4 of Oliver's thesis for details
  !****************************************************************************
  subroutine check(rhoNull, pNull,bindingEnergy, compressibility, A,B,C, tau, lambda)
    use constants, only: hbarc

    real, intent(in) :: rhoNull, pNull,bindingEnergy, compressibility, A,B,C, tau, lambda

    write(*,'(A20,F15.8)') 'A [GeV] =', A
    write(*,'(A20,F15.8)') 'B [GeV] =', B
    write(*,'(A20,F15.8)') 'C [GeV] =', C
    write(*,'(A20,F15.8)') 'Tau =', tau
    write(*,'(A20,F15.8)') 'Lambda [GeV]=', Lambda
    write(*,'(A20,F15.8)') 'Lambda [1/fm]=', Lambda/hbarc

    write(*,'(A20,E15.4,"=",E15.4,"??")') 'E/A',  E_div_A(rhoNull, rhoNull, A,B,C, tau, lambda), bindingEnergy
    write(*,'(A20,E15.4,"=",E15.4,"??")') 'd(E/A)/dRho'    ,  dRho_E_div_A(rhoNull, rhoNull, A,B,C, tau, lambda), 0.
    write(*,'(A20,E15.4,"=",E15.4,"??")') 'd^2(E/A)/dRho^2',  dRhoSquared_E_div_A(rhoNull, rhoNull, A,B,C, tau, lambda)&
         & , compressibility/9./rhoNull**2.
    write(*,'(A20,E15.4,"=",E15.4,"??")') 'U(p=0)'      ,  U(rhoNull, 0.   , rhoNull, A,B,C, tau, lambda)
    write(*,'(A20,E15.4,"=",E15.4,"??")') 'U(p=p_0)'    ,  U(rhoNull, pNull, rhoNull, A,B,C, tau, lambda)
  end subroutine check

end module skyrme
