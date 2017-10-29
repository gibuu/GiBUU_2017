!******************************************************************************
!****m* /barBar_barBarMes
! NAME
! module barBar_barBarMes
! PURPOSE
! This module implements the calculation of baryon baryon -> baryon baryon meson
! cross sections.
!******************************************************************************
module barBar_barBarMes
  implicit none
  private

  !****************************************************************************
  !****g* barBar_barBarMes/NNpi_BG
  ! SOURCE
  !
  integer, save :: NNpi_BG = 2
  ! PURPOSE
  ! Switch for the N N -> N N pi background:
  ! * 0 = no BG
  ! * 1 = BG according to Teis
  ! * 2 = BG according to Buss (improves threshold behavior, default)
  ! * 3 = BG according to Weil
  !****************************************************************************

  !****************************************************************************
  !****g* barBar_barBarMes/NNV_BG
  ! SOURCE
  !
  logical, save :: NNV_BG = .true.
  ! PURPOSE
  ! Incude a N N -> N N V background term, where V=omega,phi (in addition to
  ! possible resonance contributions).
  !****************************************************************************

  !****************************************************************************
  !****g* barBar_barBarMes/isofac_omega
  ! SOURCE
  !
  real, save :: isofac_omega = 1.
  ! PURPOSE
  ! Isospin enhancement factor for p n -> p n omega, relative to
  ! p p -> p p omega. Data indicate that this is around 2, while theory predicts
  ! even larger values (up to 5).
  ! Reference: Barsov et al., EPJ A21 (2004) 521-527.
  !****************************************************************************

  !****************************************************************************
  !****g* barBar_barBarMes/isofac_phi
  ! SOURCE
  !
  real, save :: isofac_phi = 1.
  ! PURPOSE
  ! Isospin enhancement factor for p n -> p n phi, relative to p p -> p p phi.
  ! Theory predicts values of 3-4, cf.:
  ! Kaptari, Kaempfer, Eur.Phys.J. A23 (2005) 291-304.
  !****************************************************************************


  public :: NN_NNpi_direct, sig_NNV, chooseCharge_NNpiomega
  public :: Np_NRplus_NNPion, sstatepi

  logical, save :: first = .true.

contains


  subroutine readInput

    use output, only: Write_ReadingInput

    integer :: ios

    !**************************************************************************
    !****n* barBar_barBarMes/barBar_barBarMes
    ! NAME
    ! NAMELIST barBar_barBarMes
    ! PURPOSE
    ! Includes the switches:
    ! * NNpi_BG
    ! * NNV_BG
    ! * isofac_omega
    ! * isofac_phi
    !**************************************************************************
    NAMELIST /barBar_barBarMes/ NNpi_BG, NNV_BG, isofac_omega, isofac_phi

    call Write_ReadingInput('barBar_barBarMes',0)
    rewind(5)
    read(5,nml=barBar_barBarMes,iostat=ios)
    call Write_ReadingInput('barBar_barBarMes',0,ios)
    write(*,*) 'NNpi_BG      = ', NNpi_BG
    write(*,*) 'NNV_BG       = ', NNV_BG
    write(*,*) 'isofac_omega = ', isofac_omega
    write(*,*) 'isofac_phi   = ', isofac_phi
    call Write_ReadingInput('barBar_barBarMes',1)

    first = .false.

  end subroutine readInput


  !****************************************************************************
  !****f* barBar_barBarMes/NN_NNpi_direct
  ! NAME
  ! function NN_NNpi_direct (srtfree, bg) result(Sigma_NNPion)
  ! PURPOSE
  ! Returns the background contributions to NN-> NN pion (via non-propagated resonances, plus non-resonant part).
  ! The pion cross sections contain the following channels:
  ! * 1: p p -> p p pi^0 ; which is equivalent to n n -> n n pi^0
  ! * 2: p n -> n n pi^+ ; which is equivalent to  p n -> p p pi^-
  ! * 3: p n -> p n pi^0 ;
  ! * 4: p p -> p n pi^+ ; which is equivalent to n n -> p n pi^-
  ! INPUTS
  ! * real, intent(in) :: srtfree             --- sqrt(s) in GeV
  ! OUTPUT
  ! * real, dimension (1:4) :: Sigma_NNPION   --- Cross sections in mB
  !****************************************************************************
  function NN_NNpi_direct (srtfree) result(Sigma_NNPion)
    use constants, only: mN, mPi

    real, intent(in) :: srtfree
    real, dimension (1:4) :: Sigma_NNPion  ! 1: p p -> p p pi^0 ; which is equivalent to n n -> n n pi^0
                                           ! 2: p n -> n n pi^+ ; which is equivalent to  p n -> p p pi^-
                                           ! 3: p n -> p n pi^0 ;
                                           ! 4: p p -> p n pi^+ ; which is equivalent to n n -> p n pi^-

    real :: sigmaDelta, sigma1_2, sigma3_2, sigmaAdd(4)
    real, parameter :: isofac12(1:4) = (/ 1./3., 1./3., 1./3., 2./3. /)
    real, parameter :: isofac32(1:4) = (/ 2./3., 1./3., 4./3., 10./3. /)

    if (first) call readInput

    sigma_NNPion = 0.

    ! Evaluate background cross sections for NN -> NN Pion
    if (srtfree > 2*mN + mPi) then
      sigmaAdd = sstatepi (srtFree)                ! non-resonant background

      ! pp: contribution of non-progagated resonances
      call Np_NRplus_NNPion (srtFree, .true., 1, sigmaDelta, sigma1_2, sigma3_2)
      Sigma_NNPion(1) = isofac12(1) * sigma1_2 + isofac32(1) * (sigmaDelta+sigma3_2) + sigmaAdd(1)
      Sigma_NNPion(4) = isofac12(4) * sigma1_2 + isofac32(4) * (sigmaDelta+sigma3_2) + sigmaAdd(4)

      ! pn: contribution of non-progagated resonances
      call Np_NRplus_NNPion (srtFree, .true., 0, sigmaDelta, sigma1_2, sigma3_2)
      Sigma_NNPion(2) = isofac12(2) * sigma1_2 + isofac32(2) * (sigmaDelta+sigma3_2) + sigmaAdd(2)
      Sigma_NNPion(3) = isofac12(3) * sigma1_2 + isofac32(3) * (sigmaDelta+sigma3_2) + sigmaAdd(3)
    end if

  end function NN_NNpi_direct


  !****************************************************************************
  !****s* barBar_barBarMes/Np_NRplus_NNPion
  ! NAME
  ! subroutine Np_NRplus_NNPion (srtS, background, ch, sigmaDelta, sigma1_2, sigma3_2)
  ! PURPOSE
  ! Evaluates the cross sections for N p -> N Resonance(+) -> N N pion with vacuum assumptions. It returns
  ! the cross sections sigma1_2, which is summed over all I=1/2 resonances in the intermediate state, and
  ! sigma3_2, which is summed over all I=3/2 resonances in the intermediate state.
  ! INPUTS
  ! * real, intent(in) :: srts    ! sqrt(s)
  ! * logical, intent(in) :: background
  !   .true. : do not consider any propagated resonances, therefore the background contribution is returned
  !   .false. : Consider all resonances, Therefore the full resonance contribution is returned
  ! * integer, intent(in) :: ch   ! charge of incoming nucleon
  ! OUTPUT
  ! * real, intent(out) :: sigmaDelta ! the cross section contribution of the Delta resonance
  ! * real, intent(out) :: sigma1_2   ! the cross section summed over all I=1/2 resonances
  ! * real, intent(out) :: sigma3_2   ! the cross section summed over all I=3/2 resonances (except the Delta)
  !****************************************************************************
  subroutine Np_NRplus_NNPion (srtS, background, ch, sigmaDelta, sigma1_2, sigma3_2)

    use IdTable, only: nucleon, Delta, F37_1950
    use barBar_BarBar, only: NN_NRes
    use particleProperties, only: hadron
    use dimi, only: dimiIntegrated
    use constants, only: mN
    use twoBodyTools, only: pCM

    real, intent(in) :: srts
    logical, intent(in) :: background  ! .true.  : do not consider any propagated resonances, therefore the background contribution is returned
                                       ! .false. : consider all resonances, therefore the full resonance contribution is returned
    integer ,intent(in) :: ch          ! charge of incoming nucleon
    real, intent(out) :: sigmaDelta, sigma1_2, sigma3_2

    real :: p_cm
    integer :: resID

    ! Initialize output
    sigmaDelta = 0.
    sigma1_2   = 0.
    sigma3_2   = 0.

    if (srts <= 2.*mN) return

    p_cm = pCM (srts, mN, mN)

    ! The Delta contribution
    if (.not.(background.and.hadron(Delta)%propagated) .and. (p_cm > 1E-9)) then       ! If Delta is not switched off
      ! factor 1/3 because dimixsec contains the xsection for p p->n D++. And Clebsch(N p->N D+)/(p p->n D++) = (1/4) / (3/4) = 1/3
      sigmaDelta = dimiIntegrated(srts) / (3.*p_cm)
    end if

    ! All other resonances
    do resID = Delta+1, F37_1950
       if (.not. hadron(resID)%usedForXsections) cycle           ! Exclude resonances from the Xsections
       if (background .and. hadron(resID)%propagated) cycle      ! Exclude propagated resonances for background studies
       if (hadron(resID)%IsoSpinTimes2 == 1) then ! Isospin 1/2 resonances
          sigma1_2 = sigma1_2 + NN_NRes(srts,p_cm,(/nucleon,resID/),(/1,1/),(/1,ch/),2)
       else if (hadron(resID)%IsoSpinTimes2 == 3) then ! Isospin 3/2 resonances
          ! factor 1/4 because the NN_NRES gives the result summed over all final states. And Clebsch(N p->N R+)=1/4 for resonance R with I=3/2.
          sigma3_2 = sigma3_2 + NN_NRes(srts,p_cm,(/nucleon,resID/),(/1,1/),(/1,ch/),2) / 4.
       end if
    end do

  end subroutine Np_NRplus_NNPion


  !****************************************************************************
  !****f* barBar_barBarMes/sstatepi
  ! NAME
  ! function sstatepi (srts, choice) result(sigma)
  ! PURPOSE
  ! Function for calculation of N N -> N N pi background (non-resonant), for the following channels:
  ! * 1: p p -> p p pi0
  ! * 2: p n -> n n pi+   =   p n -> p p pi-
  ! * 3: p n -> p n pi0
  ! * 4: p p -> p n pi+
  ! INPUTS
  ! * real, intent(in) :: srts                --- sqrt(s)
  ! * integer, intent(in), optional :: choice --- choice of parametrization: 1=Teis,2=Buss,3=Weil;
  !                                               if not given, the value is taken from the namelist switch 'NNpi_BG'
  ! OUTPUT
  ! * real :: sigma(4)
  !****************************************************************************
  function sstatepi (srts, choice) result(sigma)

    real, intent(in) :: srts
    integer, intent(in), optional :: choice
    real :: sigma(4)

    real :: x(5,4) , xfit
    real, parameter :: s0 = 2.015
    real :: par1(11),par2(7),par3(7),par4(10)
    integer :: param

    if (present(choice)) then
      param = choice
    else
      param = NNpi_BG
    end if
    select case (param)

    case (1) !!! Teis parametrization

       x(:,1) = (/ 61.32865, 6.1777919, 1.520738, 3.481812,      2.503206      /)
       x(:,2) = (/ 24.94817, 1.928099,  3.300971, 2.4472713e-03, 0.8496530     /)
       x(:,3) = (/ 7.25069,  2.306925,  0.883237, 3.641924,      6.6390312e-05 /)

       xfit = (srts-s0)*5
       if (xfit > 0.) then
          sigma(1) = x(1,1) * xfit**x(2,1) * exp(-(x(3,1)*xfit**x(4,1) +x(5,1)*xfit))
          sigma(2) = x(1,2) * xfit**x(2,2) * exp(-(x(3,2)*xfit**x(4,2) +x(5,2)*xfit))
          sigma(3) = x(1,3) * xfit**x(2,3) * exp(-(x(3,3)*xfit**x(4,3) +x(5,3)*xfit))
          sigma(4) = 2*sigma(1)
       else
          sigma(1:4)=0.
       end if

    case (2) !!! Buss parametrization

       par1(:) = (/ 0.00421147, 0.987631, -71.7949, 2073.54,  -24143.3, 130013., &
                                -319134.,   -16.1583, 306499.,  0.154771, 0.000823241 /)
       par2(:) = (/ 0.214637, -0.264214, 114.609, -835.136, 1615.66,  18.485,  -54.2815 /)
       par3(:) = (/ 10.6202,  -178.313,  958.58,  -40.171,  76.1355,  336.747, -425.802 /)
       par4(:) = (/ 2.50696,  -38.0018,  186.721, 14.5268,  -17.8739, 245.919, -319.183, 0.111811, 7.78858, -0.752894 /)

       xfit = srts-s0
       if (xfit > 0.) then
          sigma(1) = exp(par1(8)*xfit) * (xfit*par1(1)+xfit**2*par1(2)+xfit**3*par1(3)+xfit**4*par1(4)+xfit**5*par1(5) &
                     +xfit**6*par1(6)+xfit**7*par1(7)+xfit**8*par1(9)) / (par1(11)+(xfit-par1(10))**2)
          sigma(2) = exp(par2(6)*xfit+par2(7)*xfit**2) * (xfit*par2(1)+xfit**2*par2(2)+xfit**3*par2(3) &
                     + xfit**4*par2(4)+xfit**5*par2(5))
          sigma(3) = exp(par3(6)*xfit+par3(7)*xfit**1.2) * (par3(1)+xfit*par3(2)+xfit**2*par3(3)) / (par3(5)+(xfit-par3(4))**2)
          sigma(4) = exp(par4(6)*xfit+par4(7)*xfit**1.2) * (par4(1)+xfit*par4(2)+xfit**2*par4(3) &
                     +xfit**3*par4(4)+xfit**4*par4(5)) / (par4(9)*xfit**par4(10)+(xfit-par4(8))**2)
       else
          sigma(1:4) = 0.
       end if

    case (3) !!! Weil parametrization

       x(:,1) = (/ 14.4301, 6.32665,   19.8245,  1.28221,  -17.3336      /)
       x(:,2) = (/ 36.8751,  4.08966,  16.8907,  0.997063, -12.6746      /)
       x(:,3) = (/ 7.25069,  2.306925, 0.883237, 3.641924, 6.6390312e-05 /)
       x(:,4) = (/ 6.17963,  2.22284,  4.96358,  1.66218,  -3.65066      /)

       xfit = (srts-s0)*5
       if (xfit > 0.) then
          sigma(:) = x(1,:) * xfit**x(2,:) * exp(-(x(3,:)*xfit**x(4,:) + x(5,:)*xfit))
       else
          sigma(:) = 0.
       end if

    case default   !!! no background

      sigma = 0.

    end select

  end function sstatepi


  !****************************************************************************
  !****f* barBar_barBarMes/sig_NNV
  ! NAME
  ! real function sig_NNV (srts, charge_in)
  ! PURPOSE
  ! Calculates cross sections for
  !   1) N N -> N N omega
  !   2) N N -> N N phi
  !   3) N N -> N N pi omega
  ! in mb as function of sqrt(s), according to the parametrization (eq. 15) from
  ! Sibirtsev, Nucl. Phys. A 604 (1996) 455.
  ! INPUTS
  ! * real, intent(in) :: srts           --- sqrt(s) in GeV
  ! * integer, intent(in) :: charge_in   --- total incoming charge (2=pp, 1=pn, 0=nn)
  ! OUTPUT
  ! * real :: sig(1:3)                   --- cross section in mb for the three channels
  ! NOTE
  ! The parameters for channel (1) have been taken from Sibirtsev.
  ! The first parameter 'a' has been updated slightly by M. Abdel-Bary et al., Phys. Lett. B 647 (2007) 351.
  ! The parameters for channel (2) are from E.Ya. Paryev, J.Phys.G G36 (2009) 015103.
  ! The parameters for the third channel have been fitted to results of Pythia 6.4.25 (tuned).
  !****************************************************************************
  function sig_NNV (srts, charge_in) result (sig)

    use particleProperties, only: hadron
    use idTable, only: omegaMeson, phi
    use constants, only: mN, mPi

    real, intent(in) :: srts
    integer, intent(in) :: charge_in
    real :: sig(1:3)

    real :: srts0(1:3),x
    integer :: i
    ! Parameters (cf. Sibirtsev, table 1)
    real, parameter :: a(1:3) = (/5.3, 0.01, 1.20 /)
    real, parameter :: b(1:3) = (/2.3, 1.26, 1.55 /)
    real, parameter :: c(1:3) = (/2.4, 1.66, 1.47 /)

    sig = 0.

    if (first) call readInput
    if (.not. NNV_BG) return

    srts0(1) = 2.*mN+hadron(omegaMeson)%mass
    srts0(2) = 2.*mN+hadron(phi)%mass
    srts0(3) = 2.*mN+mPi+hadron(omegaMeson)%mass

    do i=1,3
      x = (srts0(i)/srts)**2
      if (x<1) sig(i) = a(i) * (1.-x)**b(i) * x**c(i)
    end do

    if (charge_in == 1) then
      sig(1) = sig(1) * isofac_omega
      sig(2) = sig(2) * isofac_phi
    end if

  end function sig_NNV


  !****************************************************************************
  !****f* barBar_barBarMes/chooseCharge_NNpiomega
  ! NAME
  ! subroutine chooseCharge_NNpiomega (chIn, chOut)
  ! PURPOSE
  ! Choose final state charges for N N -> N N pi omega,
  ! using Clebsch-Gordan factors which assume the pion is produced via a Delta.
  ! INPUTS
  ! * integer, dimension(1:2), intent(in) :: chIn   --- initial state charges (N,N)
  ! OUTPUT
  ! * integer, dimension(1:3) :: chOut              --- final state charges (N,N,pi)
  !****************************************************************************
  function chooseCharge_NNpiomega (chIn) result (chOut)
    use random, only: rn

    integer, dimension(1:2), intent(in) :: chIn    ! initial state charges (N,N)
    integer, dimension(1:3) :: chOut  ! final state charges (N,N,pi)

    real :: x

    x = 6*rn()

    if (chIn(1) == chIn(2)) then
      ! pp / nn
      chOut(1) = chIn(1)
      if (x<1.) then
        ! 1/6: N N pi0
        chOut(2) = chIn(2)
        chOut(3) = 0
      else
        ! 5/6: p n pi+/-
        chOut(2) = 1-chIn(2)
        chOut(3) = sum(chIn)-sum(chOut(1:2))
      end if
    else
      ! pn
      if (x<1.) then
        ! 1/6: n n pi+
        chOut = (/0,0,1/)
      else if (x<2.) then
        ! 1/6: p p pi-
        chOut = (/1,1,-1/)
      else
        ! 4/6: p n pi0
        chOut = (/1,0,0/)
      end if
    end if

  end function chooseCharge_NNpiomega



end module barBar_barBarMes
