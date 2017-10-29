module DeltaDalitz_integration

  implicit none
  private
  public :: Gamma_DeltaDalitz_integrated

contains

  ! this function gives the total width of Delta->e+e-N (integrated over dilepton mass)
  function Gamma_DeltaDalitz_integrated (m_D) result (gamma)
    use quadpack, only: qags
    use constants, only: mN, melec
    use Dilepton_Analysis, only: dGamma_dM_DeltaDalitz

    real, intent(in) :: m_D
    real :: gamma

    real,   parameter :: qag_absError=0.      ! absolute error=0  -> relative error gives the accuracy
    real,   parameter :: qag_relError=0.001   ! relative error
    real              :: qag_error_estimation
    integer           :: qag_neval,qag_error_code
    real :: a,b

    a = 2.*melec  ! lower integration boundary
    b = m_D - mN  ! upper integration boundary
    call qags (integrand,a,b,qag_absError,qag_relError,gamma,qag_error_estimation,qag_neval,qag_error_code)
    if (qag_error_code/=0 .and. qag_error_estimation/gamma>qag_relError*3.) then
      write(*,*) 'Error with QAGS in calcXS_omega_saphir:',qag_error_code
      write(*,*) 'result: ',gamma,' +- ',qag_error_estimation,' (',qag_error_estimation/gamma*100.,'%)'
      write(*,*) 'neval=',qag_neval
      write(*,*) 'm_D=',m_D
    end if
    
  contains
  
    real function integrand (m_ee)
      real, intent(in) :: m_ee
      integrand = dGamma_dM_DeltaDalitz (m_D, m_ee, 1)
    end function

  end function

end module

!*******************************************************************************


program testDileptonDelta

  ! This is a testcase for the routine 'Gamma0_DeltaDalitz'.
  ! It checks the real-photon limit and reproduces fig. 6 from
  ! Wolf et al, Nucl. Phys. A517 (1990), 615-638.

  use constants, only: pi,alphaQED,mN,mPi,melec
  use IdTable, only: delta
  use inputGeneral, only: readInputGeneral
  use ParticleProperties, only: hadron, initParticleProperties
  use Dilepton_Analysis, only: Gamma0_DeltaDalitz_Wolf, Gamma0_DeltaDalitz_Krivo, Gamma0_DeltaDalitz_Ernst, &
                               DeltaWidth_gammaN, Dilep_Init
  use baryonWidth, only: baryonWidth_gammaN
  use DeltaDalitz_integration, only: Gamma_DeltaDalitz_integrated

  implicit none

  integer, parameter :: N = 5
  real, parameter :: dm = 0.001

  real :: m_ee,m_delta(N),dGdM_W(N),dGdM_K(N),dGdM_E(N),dGdM_HT0(N),dGdM_HT1(N),m,G0,fac
  integer :: i,j

  call readInputGeneral
  call initParticleProperties
  
  call Dilep_Init (1.)

  ! check Gamma(0)
  G0 = Gamma0_DeltaDalitz_Wolf(hadron(delta)%mass,0.)
  print *, "Gamma(0)=", G0, " (Wolf)"
  G0 = Gamma0_DeltaDalitz_Krivo(hadron(delta)%mass,0.)
  print *, "Gamma(0)=", G0, " (Krivo.)"
  G0 = Gamma0_DeltaDalitz_Ernst(hadron(delta)%mass,0.)
  print *, "Gamma(0)=", G0, " (Ernst)"
  G0 = baryonWidth_gammaN(Delta,hadron(delta)%mass,0.,0)
  print *, "Gamma(0)=", G0, " (matrix,charge 0)"
  G0 = baryonWidth_gammaN(Delta,hadron(delta)%mass,0.,1)
  print *, "Gamma(0)=", G0, " (matrix,charge 1)"

  ! write dGamma/dM to file
  m_ee = dm
  m_delta = (/1.23,1.43,1.63,1.83,2.03/)
  open (911,file="dGamma_dM_DeltaDalitz.dat")
  do i=1,1100
    fac = 2.*alphaQED/(3.*pi*m_ee) * (1.+2.*melec**2/m_ee**2) * sqrt(1.-4.*melec**2/m_ee**2)
    do j=1,N
      if (m_ee+mN>m_delta(j)) then
        dGdM_W(j) = 0.
        dGdM_K(j) = 0.
        dGdM_E(j) = 0.
        dGdM_HT0(j) = 0.
        dGdM_HT1(j) = 0.
      else
        dGdM_W(j) = fac * Gamma0_DeltaDalitz_Wolf (m_delta(j), m_ee)
        dGdM_K(j) = fac * Gamma0_DeltaDalitz_Krivo (m_delta(j), m_ee)
        dGdM_E(j) = fac * Gamma0_DeltaDalitz_Ernst (m_delta(j), m_ee)
        dGdM_HT0(j) = fac * baryonWidth_gammaN (Delta,m_delta(j), m_ee, 0)
        dGdM_HT1(j) = fac * baryonWidth_gammaN (Delta,m_delta(j), m_ee, 1)
      end if
    end do
    write (911,'(26G12.5)') m_ee,dGdM_W,dGdM_K,dGdM_E,dGdM_HT0,dGdM_HT1
    m_ee = m_ee + dm
  end do
  close (911)

  ! write Gamma to file
  open (911,file="Gamma_DeltaDalitz.dat")
  m = mN + 0.01
  do i=1,120
    write (911,'(4G12.5)') m, Gamma_DeltaDalitz(m), Gamma_DeltaDalitz_integrated (m), DeltaWidth_gammaN (m, 0., 1)
    m = m + 0.01
  end do
  close (911)

  call do_dGammadM_rhoN

contains

  ! produce dGamma/dM for Delta -> rho N -> e+e-N
  subroutine do_dGammadM_rhoN
    use particleDefinition
    use histMC
    use mediumDefinition, only: vacuum
    use minkowski, only: abs4
    use master_1Body, only: decayParticle_rhoN
    use mesonWidthVacuum, only: dileptonWidth
    use mesonWidthMedium, only: WidthMesonMedium
    use output, only: timeMeasurement

    type(histogramMC) :: hist
    type(particle) D, fs(1:2)
    real :: gam, rmass, rmom(0:3)
    character(len=10) :: str
    integer :: Nevt = 20000000
    integer :: NN

    D%ID = Delta
    call CreateHistMC (hist, "dGamma/dM", 0., 1.2, 0.01, 5)

    do j=1,N
      write(str,'(f4.2,a)') m_delta(j)," GeV"
      hist%yDesc(j) = str
    end do

    print *, "do_dGammadM_rhoN: number of events = ", Nevt
    
    NN = Nevt

    do j = N,1,-1   ! Delta mass loop
      print *, "*** m = ",m_delta(j),NN," ***"
      D%mass = m_delta(j)
      D%momentum = (/ D%mass, 0., 0., 0. /)
      call TimeMeasurement (.true.)
      do i=1,NN  ! event loop
        gam = decayParticle_rhoN (D, fs)
        rmom = fs(1)%momentum
        rmass = abs4(rmom)
        gam = gam * dileptonWidth (fs(1)%ID,rmass) / WidthMesonMedium(fs(1)%ID,rmass,rmom,vacuum)
        call AddHistMC (hist, rmass, j, gam/float(NN))
      end do
      call TimeMeasurement()

      call WriteHistMC (hist, "dGamma_dM_rhoN.dat")
      NN = NN / 2.5
    end do

  end subroutine


  ! Gamma(Delta->e+e-N) - analytically integrated over dilepton mass (but wrong?!?)
  real function Gamma_DeltaDalitz(m_D)
    use constants, only: pi,alphaQED,melec,mN

    real, intent(in) :: m_D
    real :: b,c,x_m,x0,N,k
    real,parameter  :: G2 = 3.029**2  ! form factor squared at M = 0

    x0 = 4. * melec**2
    x_m = (m_D-mN)**2
    b = 2. * (m_D**2 + mN**2)
    c = (m_D**2 - mN**2)
    k = - 1. / (4.*m_D**2*mN**2)
    N = alphaQED**2 / (48.*pi) * G2 * (m_D+mN)**2 / (m_D**3*mN**2)

    Gamma_DeltaDalitz = N * ( - sqrt(c)*(x_m+b/4.) - (b*x_m+1./k)/2.*log(m_D/mN) + x_m*sqrt(c)*log(c/(x0*m_D*mN)) ) 

  end function

end program
