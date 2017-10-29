!******************************************************************************
!****m* /PhotonPionProduction_medium
! NAME
! module PhotonPionProduction_medium
! PURPOSE
! * Evaluates the photon nucleon ->  pion nucleon cross section in the medium.
! * Fermi motion of the incoming nucleon, as well as the potentials of the in and outgoing nucleons
!   are taken into account.
! * Always assumes that incoming photon travels in z-direction
!******************************************************************************
module PhotonPionProduction_medium

  use minkowski, only: pair => SP
  use electronPionProduction_kine, only: get_k_abs, pionPot

  implicit none

  private

  public :: dSigmadOmega_k_med,getKinematics_photon
  logical, save :: debug=.false.

contains


  !****************************************************************************
  !****f* PhotonPionProduction_medium/dSigmadOmega_k_med
  ! NAME
  ! function dSigmaddOmega_k_med(init_Nuc,pionCharge,energy_li,phi_k,theta_k,q,k,pf) Result(dSigma)
  !
  ! PURPOSE
  ! * Evaluates the gamma nucleon -> pion nucleon cross section as
  !   2-fold differential dSigma/dOmega(pion).
  ! * This is done using the formulas which include explicitly H_mu nu g^mu nu.
  ! * The angles of the outgoing pion are measured relative to the photon momentum "q"
  ! * Result in units of mb/sr .
  !
  ! INPUTS
  ! * type(particle),intent(in)   :: init_Nuc       ! incoming nucleon in the frame where electron runs in +z direction
  ! * integer, intent(in)         :: pionCharge     ! charge of outgoing pion
  ! * real, intent(in)            :: phi_k, theta_k ! pion scattering angles in units of DEGREE with respect to momentum transfer q
  ! * real, intent(in), dimension(0:3) :: q         ! photon 4-momentum
  !
  ! RESULT
  ! * real                              :: dSigma   ! dsigma/dOmega(final pion)
  ! * real, intent(out), dimension(0:3) :: pf       ! outgoing nucleon 4-momentum
  ! * real, intent(out), dimension(0:3) :: k        ! outgoing pion  4-momentum
  !****************************************************************************
  function dSigmadOmega_k_med(init_Nuc,pionCharge,phi_k,theta_k,q,k,pf,success) Result(dSigma)
    use degRad_conversion, only: degrees
    use constants, only: pi, mN, mPi, hbarc
    use electronPionProduction_kine, only: getV_out, get_dV_Pi_dk
    use particleDefinition
    use minkowski, only: SP

    real ::dSigma ! dsigma/dOmega(final electron)/dE(final Electron)/dOmega(final pion)

    type(particle),intent(in)   :: init_Nuc       ! incoming nucleon
    integer, intent(in)         :: pionCharge     ! charge of outgoing pion
    real, intent(in)            :: phi_k, theta_k ! pion scattering angles in units of degree
    real, intent(in), dimension(0:3) :: q         ! virtual photon 4-momentum


    real, intent(out), dimension(0:3) :: pf       ! outgoing nucleon 4-momentum
    real, intent(out), dimension(0:3) :: k        ! outgoing pion  4-momentum

    real,dimension(0:3) :: pin          ! incoming electron 4-momentum

    logical :: success

    real :: mf_n ! mass of final nucleon   =sqrt(pf^mu pf_mu)
    real :: mi_n ! mass of initial nucleon =sqrt(pin^mu pin_mu)

    integer :: charge_nucOut

    real :: kvec_abs !,qvec_abs,
    real :: pfvec_abs,Pf_freeEnergy, V_out,dV_out, dV_pi_dk
    !real :: thetaCM,

    real :: matrixElementSquared

    success=.true.

    if (debug) then
       write(*,*)
       write(*,'(A)')'+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(*,'(A)') ' Photon induced pion production'
       write(*,'(A)')'+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(*,*)
       write(*,'(A)')'###############################################'
       write(*,'(A)') 'Input:'
       write(*,*)
       write(*,'(A)')       'Initial nucleon'
       write(*,'(A,4F9.5)') ' * Momentum :',init_nuc%momentum
       write(*,'(A,3F9.5)') ' * Position :',init_nuc%position
       write(*,'(A,I9)')    ' * Charge   :',init_nuc%charge
       write(*,'(A,F9.5)')  ' * Mass     :',init_nuc%mass
       write(*,'(A,I8)')    ' * ID       :',init_nuc%ID
       write(*,'(A,L8)')    ' * perturbative      :',init_nuc%perturbative
       write(*,*)
       write(*,'(A)')       'Initial Photon'
       write(*,'(A,4F9.5)') ' * momentum :',q
       write(*,*)
       write(*,'(A)')       'Final pion'
       write(*,'(A,F9.5)')  ' * Phi of pion   :',phi_k
       write(*,'(A,F9.5)')  ' * Theta of pion :',theta_k
       write(*,*)
       write(*,'(A)')'###############################################'
    end if

    charge_nucOut=init_Nuc%charge-pionCharge

    call getKinematics_photon(pionCharge,init_Nuc,phi_k,theta_k,q,k,pf,success)


    if (.not.success) then
       ! No solution to energy/momentum conservation
       dSigma=0.
!       write(*,'(A,4E14.5)') 'no success', q
       return
    else
       pin=init_Nuc%momentum
       ! If W too small, then return:
       ! The value 1.096 is dictated by the MAID input grid which has the lowest point at W=1.1 and dW=0.01
       if (sqrt(pair(pin+q,pin+q)).lt.1.096) then
          dSigma=0.
          return
       end if

       ! Evaluate full final state kinematics and cross section
       if (debug) then
          write(*,*)
          write(*,'(A)')'###############################################'
          write(*,'(A)')'###############################################'
          write(*,'(A)')'********** Incoming :'
          write(*,'(A,4F9.5)') 'q=',q
          write(*,'(A,4F9.5)') 'pi=',pin
          write(*,'(A,4F9.5)') 'Total=' , pi+q
          write(*,*)
          write(*,'(A)')'***********Outgoing :'
          write(*,'(A,4F9.5)') 'pf=',pf
          write(*,'(A,4F9.5)') 'k=', k
          write(*,'(A,4F9.5)') 'Total=' , k+pf
          write(*,*)
          write(*,'(A,4F9.5)') 's        =' , pair(pin+q,pin+q)
          write(*,'(A,4F9.5)') 't        =' , pair(k-q,k-q)
          write(*,*)
          write(*,'(A,4F9.5)') 'W=sqrt(s)=' , sqrt(pair(pin+q,pin+q))
          write(*,'(A,4F9.5)') 'Q^2      =' , -pair(q,q)
          write(*,'(A)')'###############################################'
          write(*,'(A)')'###############################################'
          write(*,*)
          write(*,'(A,F9.5)') 'theta_pion=',degrees(acos(Dot_product(k(1:3),q(1:3)) &
               & /sqrt(Dot_product(q(1:3),q(1:3)))/sqrt(Dot_product(k(1:3),k(1:3)))))
          write(*,*)
          write(*,*) 'Checks (inserting the masses to calculate the energy):'
          write(*,'(A,4F9.5)') 'k=',sqrt(0.140**2+Dot_Product(k(1:3),k(1:3))), k(1:3)
       end if

       Kvec_ABS=sqrt(Dot_Product(k(1:3),k(1:3)))
       !Qvec_ABS=sqrt(Dot_Product(q(1:3),q(1:3)))
       Pfvec_ABS=sqrt(Dot_Product(pf(1:3),pf(1:3)))
       Pf_freeEnergy=sqrt(Dot_Product(pf(1:3),pf(1:3))+mN**2)

       mf_N=sqrt(pair(pf,pf))
       mi_N=sqrt(pair(pin,pin))

       ! Matrix element
       matrixElementSquared=matrixElement_photon(pin,pf,k,q,pionCharge,charge_nucOut)

       if (debug) write(*,*) 'mi,mf',mf_n,mi_n


       !       call getV_out(V_out,dV_out,pf,init_Nuc,pionCharge)
       !       if(debug) write(*,*) 'After getV_out',V_out,dV_out

       ! dsigma in units of 1/GEV**3


       call getV_out(V_out,dV_out,pf,init_Nuc,pionCharge)
       if (pionPot) call get_dV_Pi_dk(dV_Pi_dk,kvec_abs,init_Nuc,pionCharge)
       if (debug) write(*,*) 'After getV_out',V_out,dV_out

       ! dsigma in units of 1/GEV**3
       if (pionPot) then
          dSigma=mf_N*mi_N/(4.*pf(0)*k(0)*sqrt((SP(q,pin))**2))/((2.*pi)**2)  * kVec_abs**2 *matrixElementSquared  &
               & /abs( &
               & kvec_abs/sqrt(kvec_abs**2+mPi**2)+dV_Pi_dk       &
               & +(kvec_abs-Dot_product(pin(1:3)+q(1:3),k(1:3))/kvec_abs)/pf(0) &
               &  * (1.+1./pfVec_ABS*(2.*pfVec_ABS/Pf_freeEnergy*V_out+2.*Pf_freeEnergy*dV_out+2.*V_out*dV_out)) &
               & )
       else
          ! all potential terms neglected, only for debugging
          dSigma=mf_N*mi_N/(4.*pf(0)*k(0)*sqrt((SP(q,pin))**2))/((2.*pi)**2)  * kVec_abs**2 *matrixElementSquared   &
               & /abs( &
               & kvec_abs/k(0)+(kvec_abs-Dot_product(q(1:3),k(1:3))/kvec_abs)/pf(0) &
               & )
       end if



       ! 1/GeV**2=1/1000**2/MeV**2=1/1000**2/(1/197 fm**2)=(197/1000)**2 fm**2= (197/1000)**2 * 10 mb
       dSigma=dsigma*hbarc**2.*10
    end if


  end function dSigmadOmega_k_med




  !****************************************************************************
  !****s* PhotonPionProduction_medium/getKinematics_photon
  ! NAME
  ! subroutine getKinematics_photon(pionCharge,initNuc,phi_k,theta_k,q,k,pf,success)
  ! PURPOSE
  ! Evaluates the full kinematics for a gamma nucleon -> pion nucleon reaction.
  !
  ! The angles of the outgoing pion are measured relative to the momentum transfer "q". All angles in degree.
  ! INPUTS
  ! * real,intent(in) :: phi_k, theta_k             ! pion scattering angles in units of degree
  ! * type(particle), intent(in) :: initNuc         ! initial nucleon
  ! RESULT
  ! * real, dimension(0:3), intent(out) :: k,q,pf
  ! * 4-momenta of pion, photon and final nucleon in [GeV]
  ! * logical, intent(out) :: success ! Flag shows whether the kinematics could be established.
  !   .true.=success, .false.=no success
  !****************************************************************************
  subroutine getKinematics_photon(pionCharge,initNuc,phi_k,theta_k,q,k,pf,success)
    use particleDefinition
    use degRad_conversion, only: degrees, radian
    use rotation, only: rotateTo

    type(particle)       , intent(in)  :: initNuc
    integer              , intent(in)  :: pionCharge
    real                 , intent(in)  :: phi_k,theta_k
    real, dimension(0:3) , intent(in)  :: q
    real, dimension(0:3) , intent(out) :: k,pf
    logical, intent(out) :: success ! whether the kinematics could be established

    real :: k_abs
    !integer :: m
    !real :: mpion
    logical,parameter :: debug_kin=.false.
    real, dimension(1:3) :: k_unit

    success=.false.

    ! Define unit vector in direction of outgoing pion:
    ! Here z-axis is chosen in the direction of q(1:3)
    k_unit(1:3)=(/sin(radian(theta_k))*cos(radian(phi_k)),sin(radian(theta_k))*sin(radian(phi_k)),cos(radian(theta_k))/)

    ! Rotate this pion vector to a system where q is defining the z-axis
    k_unit = rotateTo (q(1:3), k_unit)

    if (debug) write(*,'(A,F9.5)') 'theta_pion=',degrees(acos(Dot_product(k_unit,q(1:3)) &
         & /sqrt(Dot_product(q(1:3),q(1:3)))))

    ! Solve energy momentum conservation to get absolute value of k(1:3):
    k_abs=get_K_abs(k_unit(1),k_unit(2),k_unit(3),q,initNuc,pionCharge,initNuc%charge-pionCharge,success,k)

    if (success) then
       ! Check solutions
       if (debug_kin) write(*,*) "k_abs=", k_abs
       if (k_abs.lt.0.0000001) then
          write(*,*) 'Warning: Outgoing pion momentum 0!'
          write(*,*) "k_abs=", k_abs
          !          write(*,*) 'stop'
          !          stop
       end if
       pf=q+initNuc%momentum-k
    else
       k=0
       pf=0
    end if
  end subroutine getKinematics_photon

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

  !****************************************************************************
  ! Evaluates the Matrixelement in the CM-Frame
  real function  matrixElement_photon(pi_lab,pf_lab,k_lab,q_lab,charge_pionOut,charge_nucOut)
    use degRad_conversion, only: degrees
    use lorentzTrafo

    integer :: charge_pionOut,charge_nucOut
    real, dimension(0:3),intent(in) :: pi_lab, pf_lab, k_lab,q_lab ! in- and outgoing four-momenta in lab frame
    real, dimension(0:3) :: pi, pf, k,q ! in- and outgoing four-momenta in cm frame
    !real  :: Mi, Mf ! masses of nucleons
    real,dimension(1:3) :: betaToCM


    if (debug) then
       write(*,*)
       write(*,'(A)')'+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(*,'(A)') ' Photon induced pion production: The Matrix Element'
       write(*,'(A)')'+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(*,*)
       write(*,'(A)')'###############################################'
       write(*,'(A)') 'Input (lab frame):'
       write(*,*)
       write(*,'(A,4F9.5)') 'nucleon   incoming :',pi_lab
       write(*,'(A,4F9.5)') 'nucleon   outgoing :',pf_lab
       write(*,'(A,4F9.5)') 'pion      outgoing :',k_lab
       write(*,*)
       write(*,'(A,2F9.5)') 'incoming nucleon mass:',sqrt(pair(pi_lab,pi_lab))
       write(*,'(A,2F9.5)') 'ougoing  nucleon mass:',sqrt(pair(pf_lab,pf_lab))
    end if

    !mi=sqrt(pair(pi_lab,pi_lab))
    !mf=sqrt(pair(pf_lab,pf_lab))

    ! Transforming everything to CM-Frame of the hadronic vertex:

    betatoCM = lorentzCalcBeta (q_lab+pi_lab, 'matrixElement')

    q =q_lab
    pi=pi_lab
    pf=pf_lab
    k =k_lab

    call lorentz(betaToCM, pi, 'matrixElement')
    call lorentz(betaToCM, pf, 'matrixElement')
    call lorentz(betaToCM, q , 'matrixElement')
    call lorentz(betaToCM, k , 'matrixElement')

    if (debug) then
       write(*,'(A)')'###############################################'
       write(*,'(A)') 'Converted to CM frame of hadronic vertex:'
       write(*,*)
       write(*,'(A,4F9.5)') 'nucleon   incoming :',pi
       write(*,'(A,4F9.5)') 'nucleon   outgoing :',pf
       write(*,'(A,4F9.5)') 'pion      outgoing :',k
       write(*,'(A,4F9.5)') 'virtual   photon   :',q
       write(*,*)
       write(*,'(A,2F9.5)') 'incoming nucleon mass:',sqrt(pair(pi,pi))
       write(*,'(A,2F9.5)') 'ougoing  nucleon mass:',sqrt(pair(pf,pf))
    end if

    matrixElement_photon=contraction()

  contains

    real function contraction()
      use hadronTensor_npi, only: h_munu
      use formfactors_A_main, only: getA
      use constants, only: mN

      complex,dimension(1:6) :: A
      !complex :: hadron
      real :: theta,s_Vacuum
      real, dimension(0:3) :: pi_Vacuum

      ! This is meant to solve the problem that acos(1) can sometimes give NAN since 1 is 1.000000000000001 or so.
      theta=degrees(Acos( &
           &    Min( &
           &    Max(&
           &    Dot_product(q(1:3),k(1:3))&
           &      /sqrt(Dot_product(q(1:3),q(1:3)))/sqrt(Dot_product(k(1:3),k(1:3))) &
           &    ,-0.999999999999999999) &
           &    , 0.999999999999999999)&
           ))

      pi_Vacuum(1:3)=pi(1:3)
      pi_Vacuum(0)=sqrt(mN**2+Dot_product(pi_Vacuum(1:3),pi_Vacuum(1:3)))
      s_Vacuum=pair(q+pi_Vacuum,q+pi_Vacuum)

      A=getA(charge_pionOut,charge_nucOut,theta,s_Vacuum,-pair(q,q))
!      A=getA(0,1,theta,s_Vacuum,-pair(q,q))

      contraction=-1./2.*(real(h_munu(0,0,pi,pf,k,q,A))-real(h_munu(1,1,pi,pf,k,q,A))  &
           & -real(h_munu(2,2,pi,pf,k,q,A))-real(h_munu(3,3,pi,pf,k,q,A) ))

    end function contraction

  end function matrixElement_photon



end module PhotonPionProduction_medium
