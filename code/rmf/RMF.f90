!******************************************************************************
!****m* /RMF
! NAME
! module RMF
! PURPOSE
! Includes all information about relativistic mean-field potential for baryons and mesons.
! NOTES
! * When hyperon coupling is scaled by the well known factor of 2/3, the kaons are scaled
!   by the factor 1/3 in order to compensate the missing self energy between incoming and
!   outgoing channel. This is because the threshold condition, e.g. \pi N->YK
!   sqrt(s*)>m*_y+m*_k, in the medium assumes no changes in the self energy between
!   initial and final states (see Alexei's paper on three-body collisions).
!   This prescription should better not be used at energies
!   near the kaon-production threshold, since in this method the kaon potential
!   is not consistent within Chiral Perturbation Theory or One-Boson-Exchange models.
! * The same non-trivial feature appears if the baryon self energies depend on isospin.
!   Presently no isospin-dependent part is included in the baryon fields.
! * Going beyond this simple approximation means to explicitly include different
!   threshold conditions for all channels considered in the collision term.
!******************************************************************************
module RMF

  implicit none
  private

  !****************************************************************************
  !****g* RMF/RMF_flag
  ! SOURCE
  !
  logical, save :: RMF_flag = .false.
  ! PURPOSE
  ! If .true. then use relativistic mean fields.
  !****************************************************************************


  !****************************************************************************
  !****g* RMF/N_set
  ! SOURCE
  !
  integer, save :: N_set = 1
  ! PURPOSE
  ! Select which parameter set to use:
  ! * 1 --- NL1 from G.A. Lalazissis et al., PRC 55, 540 (1997),  (K=211.29 MeV, m*/m=0.57)
  ! * 2 --- NL3 from G.A. Lalazissis et al., PRC 55, 540 (1997),  (K=271.76 MeV, m*/m=0.60)
  ! * 3 --- NL2 set from A. Lang et al., NPA 541, 507 (1992),     (K=210 MeV,    m*/m=0.83)
  ! * 4 --- NLZ2 set from M. Bender et al., PRC 60, 34304 (1999), (K=172 MeV,    m*/m=0.583)
  ! * 5 --- NL3* set from G.A. Lalazissis, private communication, (K=258.28 MeV, m*/m=0.594)
  ! * 6 --- Same as N_set=3, but including the rho meson.
  ! * 7 --- NL1 set from S.J. Lee et al., PRL 57, 2916 (1986),    (K=212 MeV,    m*/m=0.57)
  ! * 8 --- NL2 set from S.J. Lee et al., PRL 57, 2916 (1986),    (K=399 MeV,    m*/m=0.67)
  ! * 9 --- Set I from B. Liu et al., PRC 65, 045201 (2002),      (K=240 MeV,    m*/m=0.75)
  !****************************************************************************


  !****************************************************************************
  !****g* RMF/grad_flag
  ! SOURCE
  !
  logical, save, public :: grad_flag = .false.
  ! PURPOSE
  ! If .true. then include space derivatives of the fields.
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/lorentz_flag
  ! SOURCE
  !
  logical, save, public :: lorentz_flag = .true.
  ! PURPOSE
  ! If .false. then the space components of the omega field are put to zero.
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/fourMomDen_flag
  ! SOURCE
  !
  logical, save, public :: fourMomDen_flag = .false.
  ! PURPOSE
  ! If .true. then compute the four-momentum density field
  ! (not used in propagation).
  !****************************************************************************


  !******* Modification factors for the coupling constants with mean meson fields:

  !****************************************************************************
  !****g* RMF/fact_pbar
  ! SOURCE
  real, save :: fact_pbar    = 1.
  ! PURPOSE
  ! Modification factor for the antiproton coupling constants.
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/fact_hyp
  ! SOURCE
  real, save :: fact_hyp     = 1.
  ! PURPOSE
  ! Modification factor for the hyperon coupling constants.
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/fact_antihyp
  ! SOURCE
  real, save :: fact_antihyp = 1.
  ! PURPOSE
  ! Modification factor for the antihyperon coupling constants.
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/fact_Xi
  ! SOURCE
  real, save :: fact_Xi      = 1.
  ! PURPOSE
  ! Modification factor for the Xi and XiStar coupling constants.
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/fact_antiXi
  ! SOURCE
  real, save :: fact_antiXi  = 1.
  ! PURPOSE
  ! Modification factor for the antiXi and antiXiStar coupling constants.
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/fact_kaon
  ! SOURCE
  real, save :: fact_kaon    = 0.
  ! PURPOSE
  ! Modification factor for the Kaon and antikaon coupling constants.
  !****************************************************************************



  !****************************************************************************
  !****g* RMF/kaonpot_flag
  ! SOURCE
  !
  logical, save, public :: kaonpot_flag = .false.
  ! PURPOSE
  ! This switch turns on the Kaon potential in RMF mode.
  !****************************************************************************


  !**** variables related to the kaon-nucleon potential in RMF mode:
  real, parameter, public :: gs_kaon = 52.0291   ! Sigma_KN/f_pi^2 in GeV^-1 (with Sigma_KN=0.450 GeV and f_pi=0.093 GeV)
  real, parameter, public :: gv_kaon = 72.2627   ! 3/(8f_pi*^2) in GeV^-2 (with f_pi=sqrt(0.6)*f_pi)

  !**** Mean field parameters in the notations of Ref.
  !**** G.A. Lalazissis et al., PRC 55, 540 (1997):

  real, public, save :: m_nucleon       ! -- nucleon mass (GeV),
  real, public, save :: m_sigma ! -- sigma-meson mass (GeV),
  real, public, save :: m_omega ! -- omega-meson mass (GeV),
  real, public, save :: m_rho           ! -- rho-meson mass (GeV),
  real, public, save :: g_sigma ! -- sigma-nucleon coupling constant,
  real, public, save :: g_omega ! -- omega-nucleon coupling constant,
  real, public, save :: g_rho=0. ! -- rho-nucleon coupling constant,
  real, public, save :: g_2     ! -- coefficient at sigma^3 in the Lagrangian (GeV),
  real, public, save :: g_3     ! -- coefficient at sigma^4 in the Lagrangian

  real, public, save :: a_1, a_2, a_3, a_4, a_5, a_6, a_7    ! auxiliary parameters

  logical, save :: initFlag=.true.


  public :: walecka, f, fprime, fshift
  public :: getRMF_flag, getRMF_parSet, ModificationFactor


contains


  !****************************************************************************
  !****f* RMF/ModificationFactor
  ! NAME
  ! real function ModificationFactor(Id,antiFlag)
  ! PURPOSE
  ! Returns the modification factor of the RMF coupling constants for a given particle.
  ! INPUTS
  ! * integer, intent(in) :: Id         ! Id of particle
  ! * logical, intent(in) :: antiFlag   ! if .true. the particle is an antiparticle
  !****************************************************************************
  pure real function ModificationFactor(Id,antiFlag)
    use IdTable

    integer, intent(in) :: Id         ! Id of particle
    logical, intent(in) :: antiFlag   ! if .true. the particle is an antiparticle

    if ( nucleon <= Id .and. Id <= F37_1950 .and. antiFlag ) then ! Nonstrange antibaryon
      ModificationFactor=fact_pbar
    else if ( Lambda  <= Id .and. Id <= sigma_1915 ) then  ! Baryon with s=-1
      if (.not.antiFlag) then
         ModificationFactor=fact_hyp
      else
         ModificationFactor=fact_antihyp
      end if
    else if ( Xi  <= Id .and. Id <= XiStar ) then  ! Baryon with s=-2
      if (.not.antiFlag) then
         ModificationFactor=fact_Xi
      else
         ModificationFactor=fact_antiXi
      end if
    else if (id==Kaon .or. id==kaonBar) then
       ModificationFactor=fact_kaon
    else
      ModificationFactor=1.
    end if
  end function ModificationFactor


  logical function getRMF_flag()
    if (initFlag) call init
    getRMF_flag = RMF_flag
  end function getRMF_flag


  integer function getRMF_parSet()
    if (initFlag) call init
    getRMF_parSet = N_set
  end function getRMF_parSet


  !****************************************************************************
  !****s* RMF/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads input switches. Initializes the mean field parameters.
  !****************************************************************************
  subroutine init

      use constants, only: hbarc
      use output, only: Write_ReadingInput

      integer :: ios

      !************************************************************************
      !****n* RMF/RMF_input
      ! NAME
      ! NAMELIST /RMF_input/
      ! PURPOSE
      ! Includes the following input switches:
      ! * RMF_flag
      ! * N_set
      ! * grad_flag
      ! * lorentz_flag
      ! * fourMomDen_flag
      ! * kaonpot_flag
      ! * fact_pbar
      ! * fact_hyp
      ! * fact_antihyp
      ! * fact_Xi
      ! * fact_antiXi
      ! * fact_kaon
      !************************************************************************
      NAMELIST /RMF_input/ RMF_flag, N_set, grad_flag, lorentz_flag, fourMomDen_flag, kaonpot_flag, &
                           fact_pbar, fact_hyp, fact_antihyp, fact_Xi, fact_antiXi, fact_Kaon

      call Write_ReadingInput('RMF_input',0)
      rewind(5)
      read(5,nml=RMF_input,iostat=ios)
      call Write_ReadingInput('RMF_input',0,ios)
      write(*,*) 'Set RMF_flag to', RMF_flag,'.'
      if (RMF_flag) then
         write(*,*) 'Set N_set to', N_set,'.'
         write(*,*) 'Set grad_flag to', grad_flag,'.'
         write(*,*) 'Set lorentz_flag to', lorentz_flag,'.'
         write(*,*) 'Set kaonpot_flag to', kaonpot_flag,'.'
         write(*,*) 'Set fourMomDen_flag to', fourMomDen_flag,'.'
         write(*,*) 'Set fact_pbar to', fact_pbar,'.'
         write(*,*) 'Set fact_hyp to', fact_hyp,'.'
         write(*,*) 'Set fact_antihyp to', fact_antihyp,'.'
         write(*,*) 'Set fact_Xi to', fact_Xi,'.'
         write(*,*) 'Set fact_antiXi to', fact_antiXi,'.'
         write(*,*) 'Set fact_kaon to', fact_kaon,'.'
      end if
      call Write_ReadingInput('RMF_input',1)

      initFlag = .false.

      if (RMF_flag) then

         select case (N_set)

         case (1)
            ! NL1 set from G.A. Lalazissis et al., PRC 55, 540 (1997)
            ! (K=211.29 MeV, m*/m=0.57):
            m_nucleon=0.938   ! -- nucleon mass (GeV),
            m_sigma=0.492250  ! -- sigma-meson mass (GeV),
            m_omega=0.795359  ! -- omega-meson mass (GeV),
            m_rho=0.763000    ! -- rho-meson mass (GeV),
            g_sigma=10.138    ! -- sigma-nucleon coupling constant,
            g_omega=13.285    ! -- omega-nucleon coupling constant,
            g_rho=4.976       ! -- rho-nucleon coupling constant,
            g_2=-12.172*hbarc ! -- coefficient at sigma^3 in the Lagrangian (GeV),
            g_3=-36.265 ! -- coefficient at sigma^4 in the Lagrangian

         case (2)
            ! NL3 set from G.A. Lalazissis et al., PRC 55, 540 (1997)
            ! (K=271.76 MeV, m*/m=0.60) :
            m_nucleon=0.939   ! -- nucleon mass (GeV),
            m_sigma=0.508194  ! -- sigma-meson mass (GeV),
            m_omega=0.782501  ! -- omega-meson mass (GeV),
            m_rho=0.763000    ! -- rho-meson mass (GeV),
            g_sigma=10.217    ! -- sigma-nucleon coupling constant,
            g_omega=12.868    ! -- omega-nucleon coupling constant,
            g_rho=4.474       ! -- rho-nucleon coupling constant,
            g_2=-10.431*hbarc ! -- coefficient at sigma^3 in the Lagrangian (GeV),
            g_3=-28.885 ! -- coefficient at sigma^4 in the Lagrangian

         case (3)
            ! NL2 set from A. Lang et al., NPA 541, 507 (1992)
            ! (K=210 MeV, m*/m=0.83) :
            m_nucleon=0.938   ! -- nucleon mass (GeV),
            m_sigma=0.5505 ! -- sigma-meson mass (GeV),
            m_omega=0.7833  ! -- omega-meson mass (GeV),
            m_rho=0.763000    ! -- rho-meson mass (GeV),
            g_sigma=8.50    ! -- sigma-nucleon coupling constant,
            g_omega=7.54    ! -- omega-nucleon coupling constant,
            g_rho=0.      ! -- rho-nucleon coupling constant,
            g_2=-50.37*hbarc ! -- coefficient at sigma^3 in the Lagrangian (GeV),
            g_3=-6.26 ! -- coefficient at sigma^4 in the Lagrangian

         case (4)
            ! NLZ2 set from M. Bender et al., PRC 60, 34304 (1999)
            ! (K=172 MeV, m*/m=0.583) :
            m_nucleon=0.9389   ! -- nucleon mass (GeV),
            m_sigma=0.493150  ! -- sigma-meson mass (GeV),
            m_omega=0.7800  ! -- omega-meson mass (GeV),
            m_rho=0.763000    ! -- rho-meson mass (GeV),
            g_sigma=10.1369    ! -- sigma-nucleon coupling constant,
            g_omega=12.9084   ! -- omega-nucleon coupling constant,
            g_rho=4.55627       ! -- rho-nucleon coupling constant,
            g_2=-13.7561*hbarc ! -- coefficient at sigma^3 in the Lagrangian (GeV),
            g_3=-41.4013 ! -- coefficient at sigma^4 in the Lagrangian

         case (5)
            ! NL3* set from G.A. Lalazissis, private communication.
            ! (K=258.28 MeV, m*/m=0.594) :
            m_nucleon=0.939   ! -- nucleon mass (GeV),
            m_sigma=0.5026  ! -- sigma-meson mass (GeV),
            m_omega=0.7826  ! -- omega-meson mass (GeV),
            m_rho=0.763000    ! -- rho-meson mass (GeV),
            g_sigma=10.0944    ! -- sigma-nucleon coupling constant,
            g_omega=12.8065   ! -- omega-nucleon coupling constant,
            g_rho=4.5748       ! -- rho-nucleon coupling constant,
            g_2=-10.8093*hbarc ! -- coefficient at sigma^3 in the Lagrangian (GeV),
            g_3=-30.1486 ! -- coefficient at sigma^4 in the Lagrangian

         case (6)
            ! Same as N_set=3, but including the rho meson.
            m_nucleon=0.938   ! -- nucleon mass (GeV),
            m_sigma=0.5505 ! -- sigma-meson mass (GeV),
            m_omega=0.7833  ! -- omega-meson mass (GeV),
            m_rho=0.763000    ! -- rho-meson mass (GeV),
            g_sigma=8.50    ! -- sigma-nucleon coupling constant,
            g_omega=7.54    ! -- omega-nucleon coupling constant,
            g_rho=4.271      ! -- rho-nucleon coupling constant,
            g_2=-50.37*hbarc ! -- coefficient at sigma^3 in the Lagrangian (GeV),
            g_3=-6.26 ! -- coefficient at sigma^4 in the Lagrangian

         case (7)
            ! NL1 set from S.J. Lee et al., PRL 57, 2916 (1986)
            ! (K=212 MeV, m*/m=0.57) :
            m_nucleon=0.938   ! -- nucleon mass (GeV),
            m_sigma=0.49225   ! -- sigma-meson mass (GeV),
            m_omega=0.795359  ! -- omega-meson mass (GeV),
            m_rho=0.763000    ! -- rho-meson mass (GeV),
            g_sigma=10.138    ! -- sigma-nucleon coupling constant,
            g_omega=13.285    ! -- omega-nucleon coupling constant,
            g_rho=4.976       ! -- rho-nucleon coupling constant,
            g_2=-12.172*hbarc ! -- coefficient at sigma^3 in the Lagrangian (GeV),
            g_3=-36.265       ! -- coefficient at sigma^4 in the Lagrangian

         case (8)
            ! NL2 set from S.J. Lee et al., PRL 57, 2916 (1986)
            ! (K=399 MeV, m*/m=0.67) :
            m_nucleon=0.938   ! -- nucleon mass (GeV),
            m_sigma=0.50489   ! -- sigma-meson mass (GeV),
            m_omega=0.78000   ! -- omega-meson mass (GeV),
            m_rho=0.763000    ! -- rho-meson mass (GeV),
            g_sigma=9.111     ! -- sigma-nucleon coupling constant,
            g_omega=11.493    ! -- omega-nucleon coupling constant,
            g_rho=5.507       ! -- rho-nucleon coupling constant,
            g_2=-2.304*hbarc  ! -- coefficient at sigma^3 in the Lagrangian (GeV),
            g_3=13.783        ! -- coefficient at sigma^4 in the Lagrangian

         case (9)
            ! Set I from B. Liu et al., PRC 65, 045201 (2002)
            ! (K=240 MeV, m*/m=0.75)
            m_nucleon = 0.939                         ! -- nucleon mass (GeV)
            m_sigma   = 0.550                         ! -- sigma-meson mass (GeV)
            m_omega   = 0.783                         ! -- omega-meson mass (GeV)
            m_rho     = 0.763                         ! -- rho-meson mass (GeV)
            g_sigma = m_sigma * sqrt(10.33) / hbarc   ! -- sigma-nucleon coupling constant
            g_omega = m_omega * sqrt(5.42)  / hbarc   ! -- omega-nucleon coupling constant
            g_rho   = m_rho   * sqrt(0.95)  / hbarc   ! -- rho-nucleon coupling constant
            g_2 = -0.033  * g_sigma**3 * hbarc        ! -- coefficient at sigma^3 in the Lagrangian (GeV)
            g_3 = -0.0048 * g_sigma**4                ! -- coefficient at sigma^4 in the Lagrangian

         case default
            write(*,*) ' In RMF: invalid value for N_set !!!', N_set
            stop

         end select

         a_1=(g_sigma/m_sigma)**2*hbarc**3
         a_2=g_2*g_sigma/m_sigma**2
         a_3=g_3*g_sigma/m_sigma**2
         a_4=2.*g_2/m_sigma**2
         a_5=3.*g_3/m_sigma**2
         a_6=(g_omega/m_omega)**2*hbarc**3
         a_7=(g_rho/m_rho)**2*hbarc**3

      end if

  end subroutine init


  !****************************************************************************
  !****s* RMF/walecka
  ! NAME
  ! subroutine walecka(rhobar,shift,em0,rhoscalar,endens,S,V,potential)
  ! PURPOSE
  ! Determine the mass shift of the nucleon in equilibrated
  ! isospin symmetric nuclear matter at zero temperature
  ! within Walecka model with nonlinear sigma-coupling.
  ! INPUTS
  ! * real, intent(in) :: rhobar ! -- baryon density (fm^-3),
  ! * real, optional, intent(in) :: em0 ! -- starting value of mass for iterations (GeV),
  ! OUTPUT
  ! * real, intent(out) :: shift ! = m - m^* -- mass shift (GeV),
  ! * real, optional, intent(out) :: rhoscalar ! -- scalar density (fm^-3),
  ! * real, optional, intent(out) :: endens    ! -- energy density (GeV/fm^3),
  ! * real, optional, intent(out) :: pressure  ! -- pressure (GeV/fm^3),
  ! * real, optional, intent(out) :: S         ! -- scalar potential (GeV),
  ! * real, optional, intent(out) :: V         ! -- vector potential (GeV),
  ! * real, optional, intent(out) :: potential ! -- Schroedinger equivalent potential (GeV)
  !****************************************************************************
  subroutine walecka(rhobar,shift,em0,rhoscalar,endens,pressure,S,V,potential)

    use constants, only: hbarc, pi, mN

    real, intent(in) :: rhobar ! -- baryon density (fm^-3),
    real, intent(out) :: shift ! = m - m^* -- mass shift (GeV),
    real, optional, intent(in) :: em0 ! -- starting value of mass for iterations (GeV),
    real, optional, intent(out) :: rhoscalar ! -- scalar density (fm^-3)
    real, optional, intent(out) :: endens ! -- energy density (GeV/fm^3)
    real, optional, intent(out) :: pressure  ! -- pressure (GeV/fm^3),
    real, optional, intent(out) :: S         ! -- scalar potential (GeV),
    real, optional, intent(out) :: V         ! -- vector potential (GeV),
    real, optional, intent(out) :: potential ! -- potential at zero momentum (GeV)

    real :: pf, dmstm, sigma, a, rhos, drhos, fun, derfun, U_s, U_v
    integer :: niter
    logical :: flagit
    logical, save :: debug=.false.

    if (initFlag) call init

    if (rhobar < 0.001) then
      shift = 0.
      if (present(rhoscalar)) rhoscalar = 0.
      if (present(endens)) endens = 0.
      return
    end if

    pf = (1.5*pi**2*rhobar)**0.333333*hbarc

    ! Here we reset m_nucleon using the default value for the
    ! nucleon mass:
    m_nucleon = mN

    if (present(em0)) then
      dmstm = em0 - m_nucleon
    else
      dmstm = -fshift(rhobar)
    end if

    ! Test: **********************************
    !shift = fshift(rhobar)
    !return
    !**************************************************************************

    flagit = .true.
    niter = 0

    Loop_over_iterations : do while(flagit)

       niter = niter + 1

       sigma = dmstm/g_sigma    ! mean value of the sigma-meson field

       a = (dmstm + m_nucleon)/pf

       rhos = rhobar*f(a)        ! scalar density

       drhos = rhobar/pf*fprime(a) ! d rhos / d mst

       fun = dmstm + a_1*rhos + a_2*sigma**2 + a_3*sigma**3   ! we want to have fun = 0.

       derfun = 1. + a_1*drhos + a_4*sigma + a_5*sigma**2  ! d fun / d mst

       if (derfun.ne.0.) dmstm = dmstm - fun/derfun

       if ( abs(fun) <= 1.e-04 ) then
         flagit = .false.
       else if ( niter == 10 ) then
         write(*,*) ' In walecka: bad convergence after 100 iterations:',&
                   & rhobar,dmstm,abs(fun)
         stop
       end if

       if (debug) then
         write(*,'(a26,1x,f4.2,1x,i2,1x,e13.6,1x,e13.6)') &
              'rhobar, niter, dmstm, fun:', rhobar, niter, dmstm, fun
       end if

    end do Loop_over_iterations

    shift = -dmstm
    a = (dmstm + m_nucleon)/pf

    if (present(rhoscalar)) rhoscalar=rhobar*f(a)

    if (present(endens)) then
      ! Compute also the energy density and pressure:
      sigma = dmstm/g_sigma
      endens = ( 2.*pf**4/pi**2*g(a) &
            &+   0.5*(m_sigma*sigma)**2 + g_2*sigma**3/3. &
            &+   g_3*sigma**4/4. )/hbarc**3 &
            &+ 0.5*a_6*rhobar**2
      pressure = rhobar*sqrt(pf**2+(dmstm + m_nucleon)**2) &
              &-( 2.*pf**4/pi**2*g(a) &
              &+  0.5*(m_sigma*sigma)**2 + g_2*sigma**3/3. &
              &+  g_3*sigma**4/4. )/hbarc**3 &
              &+ 0.5*a_6*rhobar**2
    end if

    U_s=g_sigma*sigma
    U_v=a_6*rhobar

    if (present(S)) S=U_s
    if (present(V)) V=U_v
    if (present(potential)) potential = U_s + U_v + (U_s**2-U_v**2)/(2.*m_nucleon)

  end subroutine walecka


  !****************************************************************************
  !****f* RMF/fshift
  ! NAME
  ! real function fshift(rho)
  ! PURPOSE
  ! Fit of the nucleon mass shift m - m* for the various RMF parameter sets.
  ! INPUTS
  ! * real, intent(in) :: rho  -- baryon density (fm**-3)
  ! OUTPUT
  ! * real :: fshift  -- m - m* (GeV)
  ! NOTES
  ! This is a very rough fit which is only good to provide
  ! the starting value for iterations in walecka.
  ! The density rho must be in the interval from 0 up to 12*rhoNull.
  !****************************************************************************
  real function fshift(rho)

    use constants, only: rhoNull

    real, intent(in) :: rho

    ! constants for each parameter set (index = N_set)
    real, parameter :: A(1:9) = (/ 1.03099, 0.783659, 0.201221, 1.01769, 0.75293, 0.201221, 1.00994, 0.599025, 0.359 /)
    real, parameter :: B(1:9) = (/ 1.52439, 1.41926,  0.888537, 1.57845, 1.63993, 0.888537, 1.67593, 1.16527,  1.105 /)

    if (initFlag) call init

    fshift = m_nucleon * ( 1.-1./(1.+A(N_set)*(rho/rhoNull)**B(N_set)) )

  end function fshift


  !****************************************************************************
  !****f* RMF/f
  ! NAME
  ! real function f
  ! PURPOSE
  ! Computes analytically the expression
  !  3*a*\int_0^1 dx x^2/\sqrt(x^2+a^2)
  ! INPUTS
  ! * real, intent(in) :: a  -- dimensionless parameter equal to m^*/p_F
  !****************************************************************************
  real function f(a)

    real, intent(in) :: a
    real :: tmp

    tmp = sqrt(1.+a**2)
    f = 1.5*a * ( tmp - 0.5*a**2*log((tmp + 1.)/(tmp - 1.)) )

  end function f


  !****************************************************************************
  !****f* RMF/fprime
  ! NAME
  ! real function fprime(a)
  ! PURPOSE
  ! Computes analytically the derivative
  ! of function f(a) with respect to a.
  ! INPUTS
  ! * real, intent(in) :: a  -- dimensionless parameter equal to m^*/p_F
  !****************************************************************************
  real function fprime(a)

    real, intent(in) :: a
    real :: tmp

    tmp = sqrt(1.+a**2)
    fprime = f(a)/a + 3.*a**2*( 1./tmp - log((tmp+1.)/a) )

  end function fprime


  !****************************************************************************
  !****f* RMF/g
  ! NAME
  ! real function g(a)
  ! PURPOSE
  ! Computes analytically the expression
  !  \int_0^1 dx x^2*\sqrt(x^2+a^2)
  ! INPUTS
  ! * real, intent(in) :: a  -- dimensionless parameter equal to m^*/p_F
  !****************************************************************************
  real function g(a)

    real, intent(in) :: a
    real :: tmp

    tmp = sqrt(1.+a**2)
    g = ( tmp**3 + tmp - 0.5*a**4*log((tmp + 1.)/(tmp - 1.)) )/8.

  end function g


end module RMF
