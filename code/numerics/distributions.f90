!******************************************************************************
!****m* /distributions
! NAME
! module distributions
! PURPOSE
! This module contains various distribution functions:
! * Gaussian
! * Breit-Wigner, relativistic and non-relativistic
! * Novosibirsk
! * Woods-Saxon
! * Fermi
! * Blatt-Weisskopf
!******************************************************************************
module distributions
  implicit none
  private

  public :: gauss, novo
  public :: bW, RelBW   ! bW_norm
  public :: woods_saxon, Fermi
  public :: BlattWeisskopf
  public :: markusPostFormFactor

contains

  !****************************************************************************
  !****f* distributions/gauss
  ! NAME
  ! real function gauss (x, x0, w)
  ! PURPOSE
  ! Returns the value of the Gaussian distribution at x.
  ! INPUTS
  ! * real :: x    --- the argument (where the function is evaluated)
  ! * real :: x0   --- pole parameter
  ! * real :: w    --- width parameter
  !****************************************************************************
  real function gauss (x, x0, w)
    use constants, only: pi
    real,intent(in) :: x,x0,w
    gauss = exp(-0.5*((x-x0)/w)**2) / (w*sqrt(2*pi))
  end function gauss


  !****************************************************************************
  !****f* distributions/novo
  ! NAME
  ! real function novo (x, x0, s, t)
  ! PURPOSE
  ! Returns the value of the Novosibirsk distribution at x.
  ! INPUTS
  ! * real :: x    --- the argument (where the function is evaluated)
  ! * real :: x0   --- pole parameter
  ! * real :: s    --- width parameter sigma
  ! * real :: t    --- skewness parameter tau
  !****************************************************************************
  real function novo (x, x0, s, t)
    use constants, only: pi
    real,intent(in) :: x,x0,s,t
    real :: qx
    qx = 1. + (x-x0)/s * sinh(t*sqrt(log(4.)))/sqrt(log(4.))
    if (qx>0.) then
      novo = sinh(t*sqrt(log(4.))) * exp(-1.5*t**2) / (s*t*sqrt(2.*pi*log(4.))) * exp( -0.5*(log(qx)/t)**2 + t**2)
    else
      novo = 0.
    end if
  end function


  !****************************************************************************
  !****f* distributions/woods_saxon
  ! NAME
  ! real function woods_saxon (x, f0, x0, dx)
  ! PURPOSE
  ! Returns the value of the Woods-Saxon distribution at x.
  ! INPUTS
  ! * real :: x     --- the argument (where the function is evaluated)
  ! * real :: f0    --- normalization constant, approximately equals the value at x=0
  ! * real :: x0    --- 'radius' parameter
  ! * real :: dx    --- 'diffuseness' parameter
  !****************************************************************************
  real function woods_saxon (x, f0, x0, dx)
    real, intent(in) :: x, f0, x0, dx
    woods_saxon = f0 / (1. + exp((x-x0)/dx))
  end function woods_saxon


  !****************************************************************************
  !****f* distributions/Fermi
  ! NAME
  ! real function Fermi (E, mu, T)
  ! PURPOSE
  ! Returns the value of the Fermi distribution for the energy E,
  ! chemical potential mu and temperature T.
  ! INPUTS
  ! * real :: E     --- energy in GeV
  ! * real :: mu    --- chemical potential in GeV
  ! * real :: T     --- temperature in GeV
  !****************************************************************************
  real function Fermi (E, mu, T)
    real, intent(in) :: E, mu, T
    Fermi = 1. / (1. + exp((E-mu)/T))
  end function Fermi


  !****************************************************************************
  !****f* distributions/bW
  ! NAME
  ! real function bW(x,pole,width)
  ! PURPOSE
  ! Returns the value of the non-relativistic Breit-Wigner distribution at x
  ! with a given pole and width.
  ! INPUTS
  ! * real :: x, pole, width
  !****************************************************************************
  real function bW(x,pole,width)
    use constants, only: pi
    real, intent(in) :: x, pole, width
    bW=width/( (2*pi)*(  (x-pole)**2 + width**2/4. ))
  end function bW

!!$  !***********************************************************************
!!$  !****f* distributions/bW_norm
!!$  ! NAME
!!$  ! real function bW_norm(pole,width,a,b)
!!$  ! PURPOSE
!!$  ! Returns the Integral of the non-relativistic Breit-Wigner distribution
!!$  ! with a given pole and width in the interval (a,b).
!!$  ! INPUTS
!!$  ! * real :: pole, width, a,b
!!$  !***********************************************************************
!!$  real function bW_norm(pole,width,a,b)
!!$    use constants, only : pi
!!$    real, intent(in) :: a,b
!!$    real, intent(in) ::  pole, width
!!$
!!$    bW_norm=1/pi * ( atan(2*(b-pole)/width) - atan(2*(a-pole)/width)  )
!!$
!!$!    call qromb(f,a,b,bW_norm)
!!$!  contains
!!$!    real function f(x)
!!$!      real,intent(in) :: x
!!$!      f=bW(x,pole,width)
!!$!    end function f
!!$  end function bW_norm


  !****************************************************************************
  !****f* distributions/RelBW
  ! NAME
  ! real function RelBW(x,pole,width)
  ! NOTES
  ! Returns the value of the relativistic Breit-Wigner distribution at x
  ! with a given pole and width.
  ! INPUTS
  ! * real :: x       --- invariant mass sqrt(p_nu*p^nu) [GeV]
  ! * real :: pole    --- pole mass [GeV]
  ! * real :: width   --- width [GeV]
  ! NOTES
  ! Evaluates the spectral function of a particle with given mass and gamma
  ! at p^nu p_nu= x**2.
  ! Normalized to Integral(RelBW) from p**2=(0,infinity)GeV**2 = 1
  !****************************************************************************
  real function RelBW (x, pole, width)
    use constants, only: pi
    real, intent(in) :: x, pole, width
    real, parameter:: epsilon=0.0000000001

    if (x>epsilon) then
       RelBW=2./pi*(x**2*width)/((x**2-pole**2)**2+x**2*width**2)
    else if (x<-epsilon) then
       write(*,*) 'x in RelBW less than zero: ',x
       write(*,*) width, pole
       stop
    else
       RelBW=0.
    end if

  end function RelBW


  !****************************************************************************
  !****f* distributions/BlattWeisskopf
  ! NAME
  ! real function BlattWeisskopf(x,l)
  ! NOTES
  ! Returns the function value of the Blatt-Weisskopf-Functions,
  ! which govern the momentum dependence of the width of a resonance
  ! decaying into AB.
  ! See e.g. Effenberger Dr. thesis, page 28
  ! INPUTS
  ! * real :: x -- = p_ab *R  with
  !   p_ab = Relative Momentum of outgoing particles AB and
  !   R    = Interaction-Radius
  ! * integer :: l -- angular Momentum of outgoing particles AB
  !****************************************************************************
  real function BlattWeisskopf(x,l)
    real, intent(in)    :: x
    integer, intent(in) :: l

    select case (l)
    case (0)
       BlattWeisskopf = 1
    case (1)
       BlattWeisskopf = x/sqrt(1+x**2)
    case (2)
       BlattWeisskopf = x**2/sqrt(9+3*x**2+x**4)
    case (3)
       BlattWeisskopf = x**3/sqrt(225+45*x**2+6*x**4+x**6)
    case (4)
       BlattWeisskopf = x**4/sqrt(11025+1575*x**2+135*x**4+10*x**6+x**8)
    case default
       write(*,*) 'Wrong Angular Momentum in BlattWeisskopf:',l
       write(*,*) 'Critical Error! Stop'
       BlattWeisskopf = 0
       stop
    end select

  end function BlattWeisskopf


  !****************************************************************************
  !****f* distributions/markusPostFormFactor
  ! NAME
  ! real function markusPostFormFactor (mass, pole, srts0, lambda)
  ! PURPOSE
  ! See Post's Phd thesis, page 35, eq. (3.22).
  !****************************************************************************
  real function markusPostFormFactor (mass, pole, srts0, lambda)

    real, intent(in) :: mass, pole, lambda
    real, intent(in) :: srts0                  ! Decay threshold

    real :: L4, p2, s0
    L4 = lambda**4
    p2 = pole**2
    s0 = srts0**2

    markusPostFormFactor = (L4+0.25*(s0-p2)**2) / (L4+(mass**2-0.5*(s0+p2))**2)

  end function markusPostFormFactor


end module distributions
