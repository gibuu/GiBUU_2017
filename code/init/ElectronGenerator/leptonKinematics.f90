!******************************************************************************
!****m* /leptonKinematics
! NAME
! module leptonKinematics
!
! PURPOSE
! This module includes routines, which deal with kinematics of
! lepton-nucleon scattering.
!
!******************************************************************************

module leptonKinematics
  use particleDefinition
  use eN_eventDefinition
  implicit none

  private

  public :: &
       & evaluate_QSquared, &
       & evaluate_theta_lf,&
       & evaluate_epsilon,&
       & evaluate_W, &
       & evaluate_energy_lf,&
       & buildElectrons

! evaluate_energy_lf(Qsquared,theta_lf,energy_li,W)
! evaluate_theta_lf(QSquared,energy_li,energy_lf)

! evaluate_epsilon(theta_lf,energy_li,energy_lf)
! evaluate_W(theta_lf,energy_li,energy_lf)
! evaluate_QSquared(theta_lf,energy_li,energy_lf)


contains

  !****************************************************************************
  !****s* leptonKinematics/evaluate_energy_lf
  ! function evaluate_energy_lf(Qsquared,theta_lf,energy_li,W) result(energy_lf)
  !
  ! PURPOSE
  ! * Evaluates final energy of a lepton
  ! * Either (QSquared, energy_li, W) must be given or (theta_lf, energy_li, W)
  ! * Assumes leptons to be massless
  !
  ! INPUTS
  ! * real, intent (in),optional :: theta_lf          -- lepton scattering angle in degrees
  ! * real, intent (in),optional :: energy_li         -- initial lepton energy in GeV
  ! * real, intent (in),optional :: W                 -- center of mass energy at the hadronic vertex assuming resting nucleon
  ! * real, intent (in),optional :: Qsquared          -- photon momentum transfer Q^2
  !
  ! NOTES
  ! * ATTENTION: All inputs are optional, so you MUST use identifiers in the call-statement to specify your input!
  !   E.g.:    x=evaluate_energy_lf  (Qsquared=a,energy_li=b,W=c)
  !
  !
  ! OUTPUT
  ! *  real :: energy_lf  -- final lepton energy
  !****************************************************************************
  function evaluate_energy_lf(Qsquared,theta_lf,energy_li,W) result(energy_lf)
    use degRad_conversion, only: radian
    use constants, only: mN
    use callStack

    real, intent (in),optional :: Qsquared,theta_lf,energy_li,W
    real :: energy_lf

    if (present(Qsquared).and.present(energy_li).and.present(W)) then
       ! Evaluate energy_lf via Q^2, energy_li and W
       energy_lf=1./(2.*mn)*(-W**2-Qsquared+mn**2+2.*mn*energy_li)

    else if (present(theta_lf).and.present(energy_li).and.present(W)) then
       ! Evaluate energy_lf via theta_lf, energy_li and W
       energy_lf=(mN**2+2.*mN*energy_li-W**2)/(2.*mN+2.*energy_li*(1-cos(radian(theta_lf))))
    else
       write(*,*) 'Not enough input to evaluate_energy_lf!'
       write(*,*) present(Qsquared), present(theta_lf), present(energy_li), present(W)
       write(*,*) 'STOP!'
       call traceback()
       stop
    end if
  end function evaluate_energy_lf



  !****************************************************************************
  !****s* leptonKinematics/evaluate_epsilon
  ! real function  evaluate_epsilon(theta_lf,energy_li,energy_lf)
  !
  ! PURPOSE
  ! * Evaluates epsilon for given final and initial lepton energy and given lepton scattering angle.
  ! * Assumes leptons to be massless
  !
  ! INPUTS
  ! * real, intent (in) :: theta_lf          -- lepton scattering angle in degrees
  ! * real, intent (in) :: energy_li         -- initial lepton energy in GeV
  ! * real, intent (in) :: energy_lf         -- final lepton energy in GeV
  !
  ! OUTPUT
  ! *  real function evaluate_epsilon  -- epsilon [dimensionless]
  !****************************************************************************
  real function evaluate_epsilon(theta_lf,energy_li,energy_lf)
    use degRad_conversion, only: radian

    real, intent (in) :: theta_lf,energy_li,energy_lf
    real :: QS,vec_q_squared

    QS=evaluate_QSquared(theta_lf,energy_li,energy_lf)

    ! We use vec(q)=(energy_lf * sin(theta), 0., energy_li-energy_lf*cos(theta) )
    ! =>
    vec_q_squared=energy_lf**2+energy_li**2-2.*energy_li*energy_lf*cos(radian(theta_lf))

    ! \epsilon=\left(1+2\frac{\vec{q}\;^2}{Q^2} \tan ^2\left(\frac{\theta_f}{2} \right)  \right)^{-1}
    ! =>
    evaluate_epsilon=    1./(1.+2. * vec_q_squared/QS * (tan(radian(theta_lf/2.)))**2)

  end function evaluate_epsilon



  !****************************************************************************
  !****s* leptonKinematics/evaluate_W
  ! real function evaluate_W(theta_lf,energy_li,energy_lf)
  !
  ! PURPOSE
  ! * Evaluates W for given final and initial lepton energy and given lepton scattering angle.
  ! * Assumes resting nucleon and leptons to be massless
  !
  ! INPUTS
  ! * real, intent (in) :: theta_lf          -- lepton scattering angle in degrees
  ! * real, intent (in) :: energy_li         -- initial lepton energy in GeV
  ! * real, intent (in) :: energy_lf         -- final lepton energy in GeV
  !
  ! OUTPUT
  ! *  real function evaluate_W  -- W in GeV
  !****************************************************************************
  real function evaluate_W(theta_lf,energy_li,energy_lf)
    use degRad_conversion, only: radian
    use constants, only: mN
    real, intent (in) :: theta_lf,energy_li,energy_lf
    evaluate_W=sqrt(mn**2+2.*mn*(energy_li-energy_lf)-2.*energy_li*energy_lf*(1.-cos(radian(theta_lf))))
  end function evaluate_W


  !****************************************************************************
  !****s* leptonKinematics/evaluate_QSquared
  ! real function evaluate_QSquared(theta_lf,energy_li,energy_lf)
  !
  ! PURPOSE
  ! * Evaluates QSquared for given final and initial lepton energy and given lepton scattering angle.
  ! * Assumes leptons to be massless
  !
  ! INPUTS
  ! * real, intent (in) :: theta_lf          -- lepton scattering angle in degrees
  ! * real, intent (in) :: energy_li         -- initial lepton energy in GeV
  ! * real, intent (in) :: energy_lf         -- final lepton energy in GeV
  !
  ! OUTPUT
  ! *  real function evaluate_QSquared  -- Q^2 in GeV^2
  !****************************************************************************
  real function evaluate_QSquared(theta_lf,energy_li,energy_lf)
    use degRad_conversion, only: radian
    real, intent (in) :: theta_lf,energy_li,energy_lf
    evaluate_QSquared=2.*energy_li*energy_lf*(1.-cos(radian(theta_lf)))
  end function evaluate_QSquared


  !****************************************************************************
  !****s* leptonKinematics/evaluate_theta_lf
  ! real function evaluate_theta_lf(QSquared,energy_li,energy_lf)
  !
  ! PURPOSE
  ! * Evaluates lepton scattering angle for given QSquared and  given final and initial lepton energy.
  ! * Assumes leptons to be massless
  !
  ! INPUTS
  ! * real, intent (in) :: QSquared          -- virtual photon Q^2 in GeV^2
  ! * real, intent (in) :: energy_li         -- initial lepton energy in GeV
  ! * real, intent (in) :: energy_lf         -- final lepton energy in GeV
  !
  ! OUTPUT
  ! *  real function evaluate_theta_lf     -- lepton scattering angle in degrees
  !****************************************************************************
  real function evaluate_theta_lf(QSquared,energy_li,energy_lf)
    use degRad_conversion, only: degrees
    real, parameter :: eps=0.000001 ! =1 eV
    real, intent (in) :: Qsquared,energy_li,energy_lf
    if (energy_li.lt.eps.or.energy_lf.lt.eps) then
       write(*,*) 'Error in leptonKinematics/evaluate_theta_lf'
       write(*,*) 'Energies too small!!!', energy_li, energy_lf
       write(*,*) 'STOP!'
       stop
    end if
    evaluate_theta_lf=degrees( acos(1.-QSquared/(2.*energy_li*energy_lf)) )
  end function evaluate_theta_lf



  !****************************************************************************
  !****s* leptonKinematics/buildElectrons
  ! subroutine buildElectrons(p_in,p_out,theta,phi,energy_in,energy_out,direction_in)
  !
  ! PURPOSE
  ! * Calculate 4-momenta of an electron pair for given scattering angles theta,phi and given energies
  !   energy_in,energy_out. Direction_in denotes a space-vector pointing in the direction of the
  !   incoming electron.
  !
  ! INPUTS
  ! * real, intent(in)                  :: theta_lf [radian],phi_MC [radian],energy_li [GeV],energy_lf [GeV]
  ! * real, dimension(1:3), intent(in)  :: direction_in
  !
  ! OUTPUT
  ! * real, dimension(0:4), intent(out) :: p_in, p_out  [all in GeV]
  !****************************************************************************
  subroutine buildElectrons(p_in,p_out,theta,phi,energy_in,energy_out,direction_in)
    use rotation, only: rotateTo
    use vector, only: sphericalVector_radian,absVec

    real, dimension(0:3), intent(out) :: p_in, p_out
    real, intent(in)                  :: theta,phi,energy_in,energy_out
    real, dimension(1:3), intent(in)  :: direction_in
    real, dimension(1:3)              :: dummy,unitVector,unitVector_rotated

    p_in(0)=energy_in
    p_out(0)=energy_out

    ! Unit vector in direction of incoming electron
    unitVector=direction_in/absVec(direction_in)

    ! Create Unit vector in direction of outgoing electron in a system where incoming electron lies along z-axis
    dummy=sphericalVector_radian(theta,phi,1.)
    ! Rotate this unitVector such that incoming electron lies along unitVector
    unitVector_rotated = rotateTo (unitVector, dummy)

    p_in(1:3) =unitVector         * p_in(0)
    p_out(1:3)=unitVector_rotated * p_out(0)

    ! Check theta:
    !    if(abs(acos(Dot_Product(p_in(1:3),p_out(1:3))&
    !         & /AbsVec(p_out(1:3))/AbsVec(p_in(1:3)))-theta).gt.1E-2) then
    !       write(*,*) 'Theta',acos(Dot_Product(p_in(1:3),p_out(1:3))/AbsVec(p_out(1:3))/AbsVec(p_in(1:3)))
    !       write(*,*) 'Real theta',theta
    !       stop 'Error in buildElectrons'
    !    end if
  end subroutine buildElectrons




end module leptonKinematics
