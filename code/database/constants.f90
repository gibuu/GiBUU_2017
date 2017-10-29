!******************************************************************************
!****m* /constants
! NAME
! module constants
! PURPOSE
! Includes constants which are useful all over the code.
! SOURCE
!
module constants
  implicit none
  public

  ! math
  real, parameter    :: pi    = 3.14159265358979323846264338327950288419716939937510
  real, parameter    :: twopi = 2*pi
  complex, parameter :: ii    = (0.,1.)                       ! the imaginary unit
  complex, parameter :: zero  = (0.,0.)                       ! complex zero

  ! unit conversion
  real, parameter :: hbarc = 0.197326968                      ! in GeV fm
  real, parameter :: GeVSquared_times_mb = 0.1*hbarc**(-2)    ! = 1/(0.389 GeV^2 mb)

  ! couplings and mixing angles
  real, parameter :: alphaQED=0.007297                        ! fine-structure constant alpha
  real, parameter :: electronChargeSQ=alphaQED*4.*pi          ! e**2 = 4*pi*alpha
  real, parameter :: GF=1.16637e-5                            ! Fermi's constant in GeV^-2
  real, parameter :: coscab=0.9745                            ! cosine of Cabbibo mixing angle
  real, parameter :: sinsthweinbg=0.2228                      ! sin_theta_Weinberg squared!!!
  real, parameter :: f_pi=0.093                               ! pion weak decay constant in GeV
  real, parameter :: g_A=1.26                                 ! axial nucleon coupling used in nonlinear sigma model

  ! particle masses (in GeV):
  ! leptons
  real, parameter :: melec = 0.00051099892                    ! electron mass
  real, parameter :: mmuon = 0.105658369                      ! muon mass
  real, parameter :: mtau  = 1.77699                          ! tau mass
  ! hadrons
  real, parameter :: mPi   = 0.138                            ! pion mass
  real, parameter :: mN    = 0.938                            ! nucleon mass
  real, parameter :: mK    = 0.496                            ! Kaon mass

  ! nuclear physics
  real, parameter :: rhoNull = 0.168                          ! normal nuclear density in fm**(-3)
  real, parameter :: massu   = 0.931494028                    ! mass of atomic unit in GeV (pdg.lbl.gov 2009)

  ! numerics
  integer, parameter :: singlePrecision = 4                   ! number of bytes for single precision variables

!******************************************************************************


end module constants
