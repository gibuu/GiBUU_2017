!******************************************************************************
!****m* /winkel_tools
! PURPOSE
! Contains the slope b-parameters and other routines used for sampling
! the angular distributiuons in two-body scattering.
!******************************************************************************
module winkel_tools

  implicit none
  private

  public :: Cugnon_bpp, Cugnon_bnp
  public :: tSlope_Effenb, tSlope_PYTHIA
  public :: tSlope_EL_pbarp, tSlope_CEX_pbarp
  public :: dsigdt_Regge, dsigdt_LE

contains


  !****************************************************************************
  !****f* winkel_tools/Cugnon_bpp
  ! NAME
  ! real function Cugnon_bpp(plab)
  ! NOTES
  ! b-coefficients from new Cugnon parameterization
  ! INPUTS
  ! * plab -- lab. momentum (GeV/c)
  !****************************************************************************
  real function Cugnon_bpp(plab)
    real :: plab
    if ( plab.lt.2. ) then
       Cugnon_bpp = 5.5*plab**8 / (7.7+plab**8)
    else
       Cugnon_bpp = 5.334 + 0.67*(plab-2.)
    end if
  end function Cugnon_bpp

  !****************************************************************************
  !****f* winkel_tools/Cugnon_bnp
  ! NAME
  ! real function Cugnon_bnp(plab)
  ! NOTES
  ! b-coefficients from new Cugnon parameterization
  ! INPUTS
  ! * plab --lab. momentum (GeV/c)
  !****************************************************************************
  real function Cugnon_bnp(plab)
    real :: plab
    if ( plab.lt.0.225 ) then
       Cugnon_bnp = 0.
    else if ( plab.lt.0.6 ) then
       Cugnon_bnp = 16.53*(plab-0.225)
    else if ( plab.lt.1.6 ) then
       Cugnon_bnp = -1.63*plab + 7.16
    else
       Cugnon_bnp = Cugnon_bpp(plab)
    end if
  end function Cugnon_bnp

  !****************************************************************************
  !****f* winkel_tools/tSlope_Effenb
  ! NAME
  ! real function tSlope_Effenb(sqrts)
  ! NOTES
  ! b-coefficient for pp elastic scattering parameterization
  ! from J. Cugnon et al., NPA 352, 505 (1981)
  ! INPUTS
  ! * sqrts -- inv. energy (GeV)
  !****************************************************************************
  real function tSlope_Effenb(sqrts)
    use constants, only: mN
    real, intent(in) :: sqrts
    real :: h
    h = (3.65*(sqrts-2*mN))**6
    tSlope_Effenb = 6*h/(1.0+h)
  end function tSlope_Effenb

  !****************************************************************************
  !****f* winkel_tools/tSlope_PYTHIA
  ! NAME
  ! real function tSlope_PYTHIA(sqrts, ID1, ID2)
  ! NOTES
  ! b-coefficient for pp elastic scattering,
  ! taken from PYTHIA model, see T. Falter et al., PRC 70, 054609 (2004)
  ! INPUTS
  ! * sqrts -- inv. energy (GeV)
  ! * ID1, ID2 --- ids of colliding particles
  !****************************************************************************
  real function tSlope_PYTHIA(sqrts, ID1, ID2)
    real, intent(in):: sqrts
    integer, intent(in):: ID1, ID2

    real :: b1, b2
    real, parameter :: eps=0.0808

    b1 = tSlope_PYTHIA_b(ID1)
    b2 = tSlope_PYTHIA_b(ID2)

    tSlope_PYTHIA = 2*(b1+b2)+4*sqrts**(2*eps)-4.2
    return

  contains

    real function tSlope_PYTHIA_b(ID)
      integer ID

      if (ID.le.100) then        ! Proton
         tSlope_PYTHIA_b = 2.3
      else if (ID.eq.109) then   ! J/psi
         tSlope_PYTHIA_b = 0.23
      else
         tSlope_PYTHIA_b = 1.4
      end if
      return
    end function tSlope_PYTHIA_b

  end function tSlope_PYTHIA

  !****************************************************************************
  !****f* winkel_tools/tSlope_EL_pbarp
  ! NAME
  ! real function tSlope_EL_pbarp(sqrts)
  ! NOTES
  ! b-coefficient for pbar p elastic scattering,
  ! taken from L.A. Kondratyuk and M.G. Sapozhnikov, Sov. J. Nucl. Phys. 46, 56 (1987).
  ! INPUTS
  ! * sqrts -- inv. energy (GeV)
  !****************************************************************************
  real function tSlope_EL_pbarp(sqrts)
    use constants, only: hbarc, mN
    real, parameter :: A=0.67, B=0.35
    real :: sqrts, qcm
    qcm=sqrt(sqrts**2/4.-mN**2)
    tSlope_EL_pbarp=(A+B*hbarc/qcm)**2/hbarc**2
  end function tSlope_EL_pbarp

  !****************************************************************************
  !****f* winkel_tools/tSlope_CEX_pbarp
  ! NAME
  ! real function tSlope_CEX_pbarp(sqrts)
  ! NOTES
  ! b-coefficient for pbar p --> nbar n charge exchange,
  !   parameterization of available HEP data at plab=0.7-35 GeV/c by A.L.
  ! INPUTS
  ! * sqrts -- inv. energy (GeV)
  !****************************************************************************
  real function tSlope_CEX_pbarp(sqrts)
    use constants, only: mN
    real :: sqrts, plab
    plab=sqrt( (sqrts**2/(2.*mN))**2-sqrts**2 )
    tSlope_CEX_pbarp= 11.*exp(-0.23*plab)+8.*plab**2.2/(254.+plab**2.2)
  end function tSlope_CEX_Pbarp


  !****************************************************************************
  !****f* winkel_tools/dsigdt_LE
  ! NAME
  ! function dsigdt_LE (sqrts, prcm) result(costheta)
  ! PURPOSE
  ! Choose LambdaBar polar angle for Nbar N -> LambdaBar Lambda.
  ! INPUTS
  ! * sqrts -- inv. energy (GeV)
  ! * prcm  -- c.m. momentum of colliding particles (GeV/c)
  ! OUTPUT
  ! * costheta -- cosine of the polar angle
  ! NOTES
  ! The fit by Theodoros Gaitanos is used: dSigDt = a0 + a1*exp(a2*(tmin-t)).
  ! The angle is chosen with respect to the Nbar c.m. momentum.
  ! Valid at p_lab < 1.8 GeV/c.
  !****************************************************************************
  function dsigdt_LE (sqrts, prcm) result(costheta)

    use random, only: rn
    use IdTable, only: Lambda
    use particleProperties, only: hadron

    real, intent(in) :: sqrts, prcm
    real :: costheta

    real :: exc,a0,a1,a2,prcm_out,t_min,t_max,dt,sig1,sig2,x,t
    real, parameter :: aa = -4.66944, bb = 1.49962

    exc = sqrts-2.*hadron(Lambda)%mass
    a0 = (19.1946+132202*exc**2.5)/(1.+42592*exc**3.5)
    a1 = (277.044*exc**0.25-3145.39*exc+4.99314e+06*exc**2.8)/(1.+19075.7*exc**3)
    if (exc<=0.0033) then
       a2 = -(-67357.8*exc+257.692)
    else
       !--> fit ohne die letzten 2 Punkte:
       a2 = -aa*alog(bb*exc)
    end if
    prcm_out = sqrts**2/4. - hadron(Lambda)%mass**2
    prcm_out = sqrt( max(0.,prcm_out) )
    t_min=-(prcm-prcm_out)**2
    t_max=-(prcm+prcm_out)**2
    dt=t_min-t_max
    sig1=a0*dt
    sig2=a1/a2*(exp(a2*dt)-1.)
    x = rn()*(sig1+sig2)
    if (x<=sig1) then
       x= rn()
       t=t_max+x*dt
    else
       x= rn()
       t=-alog(x*exp(-a2*t_min)+exp(-a2*t_max)*(1.-x))/a2
    end if
    costheta= (2.*t-t_max-t_min)/(t_min-t_max)

  end function dsigdt_LE


  !****************************************************************************
  !****f* winkel_tools/dsigdt_Regge
  ! NAME
  ! function dsigdt_Regge (sqrts, prcm, n) result (costheta)
  ! PURPOSE
  ! Choose the antihyperon polar angle for Nbar N -> Ybar Y.
  ! INPUTS
  ! * sqrts -- inv. energy (GeV)
  ! * prcm  -- c.m. momentum of colliding particles (GeV/c)
  ! * n = 1 --- choose LambdaBar polar angle for Nbar N -> LambdaBar Lambda
  !     = 2 --- choose LambdaBar polar angle for Nbar N -> LambdaBar Sigma
  !             or Lambda polar angle for        Nbar N -> SigmaBar Lambda
  ! OUTPUT
  ! * costheta -- cosine of the polar angle
  ! NOTES
  ! The angle is chosen with respect to the Nbar c.m. momentum for outgoing LambdaBar
  ! or with respect to the N c.m. momentum for outgoing Lambda.
  ! Used phenomenological Regge-like fit from B. Sadoulet, CERN-HERA 69-2 :
  ! dSigDt \propto sigma*a*exp(a*t)*c**(2.*alpha) + tau*b**2/4.*sqrt(-t)*exp(-b*sqrt(-t))*c**(2.*beta)
  ! Valid at p_lab > 1.8 GeV/c
  !****************************************************************************
  function dsigdt_Regge (sqrts, prcm, n) result (costheta)
    use random, only: rn
    use IdTable, only: Lambda,SigmaResonance
    use particleProperties, only: hadron
    use constants, only: mN

    real, intent(in) :: sqrts, prcm   ! [GeV]
    integer, intent(in) :: n
    real :: costheta

    !**** Parameters for the pbar p -> LambdaBar Lambda  (1)
    !****           and      pbar p -> SigmaBar^0 Lambda (2) cross sections:
    real, dimension(1:2), parameter :: sigma=(/808.,295./), a=(/8.5,8.5/), alpha=(/0.28,0.41/), &
                                       tau=(/1280.,580./), b=(/3.5,3.5/), beta=(/0.35,0.73/)
    real :: s,prcm_out,ta,tb,t_min,t_max,t,sum_m2,c,exp_min,exp_max,f_min,f_max,f0,sig1,sig2
    real :: x,srtt,fun,deriv
    integer :: i

    s=sqrts**2
    select case (n)
    case (1)
       prcm_out=s/4.-hadron(Lambda)%mass**2
       prcm_out = sqrt( max(0.,prcm_out) )
       ta=-prcm**2-prcm_out**2
       sum_m2=2.*(mN**2+hadron(Lambda)%mass**2)
    case (2)
       prcm_out=(s+hadron(SigmaResonance)%mass**2-hadron(Lambda)%mass**2)**2/(4.*s) &
                -hadron(SigmaResonance)%mass**2
       prcm_out= sqrt( max(0.,prcm_out) )
       ta= mN**2 + hadron(Lambda)%mass**2 &
            - 2.*sqrt(prcm**2+mN**2)*sqrt(prcm_out**2+hadron(Lambda)%mass**2)
       sum_m2=2.*mN**2+hadron(Lambda)%mass**2+hadron(SigmaResonance)%mass**2
    case default
       write(*,*) ' Wrong mode in dsigdt_Regge : ', n
       stop
    end select
    tb=2.*prcm*prcm_out

    t_min=ta+tb
    t_max=ta-tb

    !  write(*,*)' t_min, t_max : ', t_min, t_max

    c=(2.*s-sum_m2)/(17.48-sum_m2)
    exp_min=exp(a(n)*t_min)
    exp_max=exp(a(n)*t_max)
    f_min=f(sqrt(max(0.,-t_min)))
    f_max=f(sqrt(max(0.,-t_max)))

    !  write(*,*)' f_min, f_max :', f_min, f_max
    !  stop

    sig1= sigma(n) * c**(2.*alpha(n)) * ( exp_min - exp_max )
    sig2=tau(n) * b(n)**2/4. * c**(2.*beta(n)) * ( f_min - f_max )
    x = rn()*(sig1+sig2)
    if (x<=sig1) then
       x= rn()
       t=alog( x*exp_min + (1.-x)*exp_max ) / a(n)
    else
       x= rn()
       t=0.5*(t_min+t_max)
       f0=f_max-x*(f_max-f_min)
       do i=1,10
          srtt=sqrt(max(0.,-t))
          fun=f(srtt)-f0
          deriv=dfdt(srtt)
          if (deriv/=0.) t=t-fun/deriv
          if (t>=0.) t=t_min
          if (abs(fun)<=1.e-06) exit
       end do
       if (abs(fun)>1.e-06) then
          write(*,*) ' Problem in dsigdt_Regge, fun: ', fun
          write(*,*) ' sqrts, prcm, n: ', sqrts, prcm, n
          write(*,*) ' x, t, t_min,t_max: ', x, t, t_min,t_max
       end if
    end if

    costheta=max(-1.,min(1.,(2.*t-t_max-t_min)/(t_min-t_max)))

  contains

    real function f(x)
      real, intent(in) :: x
      f=2./b(n)**3*((b(n)*x+1.)**2+1.)*exp(-b(n)*x)
    end function f

    real function dfdt(x)
      real, intent(in) :: x
      dfdt=x*exp(-b(n)*x)
    end function dfdt

  end function dsigdt_Regge


end module winkel_tools
