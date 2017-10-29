!******************************************************************************
!****m* /electronPionProduction_medium_eN
! NAME
! module  electronPionProduction_medium_eN
! PURPOSE
! * Evaluates the electron nucleon -> electron pion nucleon cross section as
!   5-fold differential dSigma/dOmega(final e^-) /dE(final e^-) /dOmega(pion)
!   in the medium.
! * Fermi motion of the incoming nucleon, as well as the potentials of the in
!   and outgoing nucleons are taken into account.
! * While the module "electronPionProduction_medium" always assumes that the
!   incoming electron travels in z-direction, this module deals with
!   arbitrary momenta
!******************************************************************************
module electronPionProd_medium_eN
  use electronPionProduction_kine, only: pionPot, getKinematics_eN
  use minkowski, only: pair => SP

  implicit none

  private

  public :: dSdO_fdE_fdO_k_med_eN, getKinematics_eN

  logical, save :: debug=.false.


contains

  !****************************************************************************
  !****f* electronPionProduction_medium_eN/dSdO_fdE_fdO_k_med_eN
  ! NAME
  ! function dSdO_fdE_fdO_k_med_eN(eN, pionCharge, phi_k, theta_k, k, pf, processID,pionNucleonSystem)
  !
  ! PURPOSE
  ! * Evaluates the electron nucleon -> electron pion nucleon cross section as
  !   5-fold differential
  !   dSigma/dOmega(final e^-) /dE(final e^-) /dOmega(pion).
  ! * This is done using the formulas which include explicitly H_mu nu L^mu nu.
  ! * The angles of the outgoing pion are measured relative to the momentum
  !   transfer "q"
  ! * Result in units of mcb/MeV=10^3 mcb/GeV=mb/GeV
  !
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN       --
  !   underlying electron-nucleon event
  ! * integer                     :: pionCharge     -- charge of outgoing pion
  ! * real                        :: phi_k, theta_k --
  !   pion scattering angles in units of degree
  ! * integer           , optional:: processID      --
  !   EM=Electromagnetic (see module leptonicID for definition)
  ! * integer           , optional:: pionNucleonSystem     --
  !   1 = pion angles phi_k and theta_k are defined in lab,
  !   2 = they are defined in the CM frame of  the incoming nucleon and photon
  !
  !
  ! RESULT
  ! * real ::dSigma --
  !   dsigma/dOmega(final electron)/dE(final Electron)/dOmega(final pion)
  ! * real, dimension(0:3) :: pf       -- outgoing nucleon 4-momentum
  ! * real, dimension(0:3) :: k        -- outgoing pion  4-momentum
  !****************************************************************************
  function dSdO_fdE_fdO_k_med_eN(eN, pionCharge, phi_k, theta_k, k, pf, processID,pionNucleonSystem) Result(dSigma)
    use degRad_conversion, only: degrees, radian
    use constants, only: pi, mN, mPi, hbarc
    use particleDefinition
    use leptonicID
    use vector, only: absVec, theta_IN
    use eN_eventDefinition, only:electronNucleon_event,write_electronNucleon_event
    use electronPionProduction_kine, only: getV_out, get_dV_Pi_dk
    use lorentzTrafo, only: lorentz


    implicit none

    real ::dSigma ! dsigma/dOmega(final electron)/dE(final Electron)/dOmega(final pion)

    type(electronNucleon_event) , intent(in)  :: eN
    integer                     , intent(in)  :: pionCharge
    real                        , intent(in)  :: phi_k, theta_k
    integer           , optional, intent(in)  :: processID
    integer           , optional, intent(in)  :: pionNucleonSystem

    real        , intent(out), dimension(0:3) :: pf
    real        , intent(out), dimension(0:3) :: k

    real, dimension(0:3) :: lf       ! outgoing electron 4-momentum
    real, dimension(0:3) :: q        ! virtual photon 4-momentum
    real, dimension(0:3) :: li       ! incoming electron 4-momentum
    real, dimension(0:3) :: pin      ! incoming electron 4-momentum


    logical :: success,twoRoots

    real :: mf_n ! mass of final nucleon   =sqrt(pf^mu pf_mu)
    real :: mi_n ! mass of initial nucleon =sqrt(pin^mu pin_mu)

    integer :: charge_nucOut
    integer :: process_ID
    real, dimension(1:3) :: electron_velocity
    real :: kvec_abs, pfvec_abs, lfvec_ABS, Pf_freeEnergy
    real :: dV_out, V_out,relativeVelocity
    real :: dV_Pi_dk, jacobian, kcm_abs, theta_k_lab
    real, dimension (0:3) :: kcm
    real, dimension (1:3) :: betaToCM

    if (debug) then
       call write_electronNucleon_event(eN,DoShort_=.true.)

       write(*,'(A)')'###############################################'
       write(*,*)
       write(*,'(A)')      'Final pion   :'
       write(*,'(A,F9.5)') 'Phi of pion       :',phi_k
       write(*,'(A,F9.5)') 'Theta of pion     :',theta_k
       write(*,*)
       write(*,'(A)')'###############################################'
    end if

    process_ID=EM
    if (present(processID))process_ID=processID

    charge_nucOut=eN%nucleon%charge-pionCharge
    if (process_ID.eq.CC) charge_nucOut=charge_nucOut+1  !for negatively charged outgoing lepton


    li=eN%lepton_in%momentum
    lf=eN%lepton_out%momentum
    q =eN%boson%momentum
    if (present(pionNucleonSystem)) then
       call getKinematics_eN(eN,pionCharge,charge_nucOut,phi_k,theta_k,&
            & k,pf,twoRoots,success,pionNucleonSystem)
       if (success) then
          theta_k_lab=theta_In( k(1:3),en%boson%momentum(1:3))
          if (pionNucleonSystem.ne.2) then
             ! Check
             if (abs((theta_k_lab-theta_k)/theta_k).gt.1E-5) then
                write(*,*) 'Error in getKinematics_EN'
                write(*,*) theta_k_lab , theta_k, pionNucleonSystem
                STOP 'in module electronPionProduction_medium_eN.f90'
             end if
          end if
       end if
    else
       call getKinematics_eN(eN,pionCharge,charge_nucOut,phi_k,theta_k,&
            & k,pf,twoRoots,success)
    end if
    if (debug)  write(*,*) 'After getKinematics:', success


    if (.not.success) then
       ! No solution to energy/momentum conservation
       dSigma=0.
       return
    else
       pin=eN%nucleon%momentum
       ! If W too small, then return:
       if (eN%W_free.lt.1.076) then
          if (debug) write(*,*) 'W_free.lt.1.076',eN%W_free
          dSigma=0.
          return
       end if
       ! Evaluate full final state kinematics and cross section
       if (debug) then
          write(*,*)
          write(*,'(A)')'###############################################'
          write(*,'(A)')'###############################################'
          write(*,'(A)')'********** Incoming :'
          write(*,'(A,4F9.5)') 'li=',li
          write(*,'(A,4F9.5)') 'pi=',pin
          write(*,'(A,4F9.5)') 'Total=' , li+pin
          write(*,*)
          write(*,'(A)')'********** Virtual Photon :'
          write(*,'(A,4F9.5)') 'q=',q
          write(*,*)
          write(*,'(A)')'***********Outgoing :'
          write(*,'(A,4F9.5)') 'lf=',lf
          write(*,'(A,4F9.5)') 'pf=',pf
          write(*,'(A,4F9.5)') 'k=',k
          write(*,'(A,4F9.5)') 'Total=' , k+lf+pf
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

       lfvec_ABS=sqrt(Dot_Product(lf(1:3),lf(1:3)))
       Kvec_ABS=sqrt(Dot_Product(k(1:3),k(1:3)))
       Pfvec_ABS=sqrt(Dot_Product(pf(1:3),pf(1:3)))
       Pf_freeEnergy=sqrt(Dot_Product(pf(1:3),pf(1:3))+mN**2)

       mf_N=sqrt(pair(pf,pf))
       mi_N=sqrt(pair(pin,pin))

       electron_velocity=eN%lepton_in%momentum(1:3)/eN%lepton_in%momentum(0)
       relativeVelocity=absVec(electron_velocity-eN%nucleon%velocity)

       if (debug) write(*,*) 'mi,mf,relativeVelocity',mf_n,mi_n,relativevelocity

       call getV_out(V_out,dV_out,pf,eN%nucleon,pionCharge)
       if (pionPot) call get_dV_Pi_dk(dV_Pi_dk,kvec_abs,eN%nucleon,pionCharge)
       if (debug) write(*,*) 'After getV_out',V_out,dV_out


       ! dsigma in units of 1/GEV**3
       if (pionPot) then
          dSigma=mf_N*mi_N/relativeVelocity/((2.*pi)**5)/2. &
               & *kVec_abs**2*lfvec_ABS/k(0)/li(0)/pf(0)/pin(0)*&
               & matrixElement_eN(pin,pf,li,lf,k,q,pionCharge,charge_nucOut,process_ID,eN%W_free) &
               & /abs( &
               & kvec_abs/sqrt(kvec_abs**2+mPi**2)+dV_Pi_dk &
               & +(kvec_abs-Dot_product(pin(1:3)+q(1:3),k(1:3))/kvec_abs)/pf(0) &
               &  * (1.+1./pfVec_ABS*(2.*pfVec_ABS/Pf_freeEnergy*V_out+2.*Pf_freeEnergy*dV_out+2.*V_out*dV_out)) &
               & )
       else
          dSigma=mf_N*mi_N/relativeVelocity/((2.*pi)**5)/2. &
               & *kVec_abs**2*lfvec_ABS/k(0)/li(0)/pf(0)/pin(0)*&
               & matrixElement_eN(pin,pf,li,lf,k,q,pionCharge,charge_nucOut,process_ID,eN%W_free) &
               & /abs( &
               & kvec_abs/k(0)+(kvec_abs-Dot_product(pin(1:3)+q(1:3),k(1:3))/kvec_abs)/pf(0) &
               &  * (1.+1./pfVec_ABS*(2.*pfVec_ABS/Pf_freeEnergy*V_out+2.*Pf_freeEnergy*dV_out+2.*V_out*dV_out)) &
               & )
       end if


       if (debug) then
          write(*,*) ' kvec_abs/k(0)',  kvec_abs/k(0)
          write(*,*) '+ ... =',(kvec_abs-Dot_product(pin(1:3)+q(1:3),k(1:3))/kvec_abs)/pf(0) &
               &  * (1.+1./pfVec_ABS*(2.*pfVec_ABS/Pf_freeEnergy*V_out+2*Pf_freeEnergy*dV_out+V_out*dV_out))
          write(*,*) '1>> ... ?', 1./pfVec_ABS*(2.*pfVec_ABS/Pf_freeEnergy*V_out+2*Pf_freeEnergy*dV_out+2.*V_out*dV_out)
       end if

       ! 1/GeV**2=1/1000**2/MeV**2=1/1000**2/(1/197 fm**2)=(197/1000)**2 fm**2= (197/1000)**2 * 10 mb
       dSigma=dsigma*hbarc**2.*10

       if (present(pionNucleonSystem)) then
          if (pionNucleonSystem.eq.2) then
             ! Transform dOmega_pion(Calculation Frame) to dOmega_pion(Center of Mass Frame of photon nucleon system)
             ! boost pion to CM system
             kcm=k
             betaToCM=(pin(1:3)+q(1:3))/(pin(0)+q(0))
             call lorentz(betaToCM, kcm)
             kcm_abs=absVec(kcm(1:3))
             jacobian=Dot_product(k(1:3),k(1:3))*sqrt(pair(pin+q,pin+q))/&
                  & kCM_abs/abs(   (Sqrt(pair(pin,pin))+q(0))*sqrt(Dot_Product(k(1:3),k(1:3)))&
                  & -k(0)*sqrt(Dot_Product(q(1:3),q(1:3)))* cos(radian(theta_k_lab)))
             dsigma=dsigma/jacobian
          end if
       end if


       if (twoRoots) then
          ! There were two solutions for the kinematics from which we picked randomly one. So we
          ! now need to multiply the cross section by 2, such that we get on average sigma(kinematics1)+sigma(kinematics2)
          dsigma=dsigma*2.
       end if
       if (debug) write( *,*) 'dsigma=',dsigma

    end if
  end function dSdO_fdE_fdO_k_med_eN


  !****************************************************************************
  ! Evaluates the Matrixelement in the CM-Frame
  real function  matrixElement_EN(pi_lab,pf_lab,li_lab,lf_lab,k_lab,q_lab,charge_pionOut,charge_nucOut,processID,W_free)
    use degRad_conversion, only: degrees
    use lorentzTrafo
    use leptonicID
    implicit none

    real   , intent(in) :: W_free
    integer, intent(in) :: charge_pionOut,charge_nucOut
    real, dimension(0:3),intent(in) :: pi_lab, pf_lab, li_lab,lf_lab,k_lab,q_lab ! in- and outgoing four-momenta in lab frame
    integer, intent(in),optional :: processID
    real, dimension(0:3) :: pi, pf, li,lf,k,q ! in- and outgoing four-momenta in cm frame
    !real  :: Mi, Mf ! masses of nucleons
    !real  :: mie, mfe ! masses of electrons
    real,dimension(1:3) :: betaToCM
    !real :: u_ms, s_ms, t_ms  ! Mandelstam variables
    !real :: theta
    !real :: mpion
    !real :: kcm
    !integer :: process_ID

    !process_ID=EM
    !if(present(processID)) process_ID=processID

    !mpion=mPi

    if (debug) then
       write(*,*)
       write(*,'(A)')'+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(*,'(A)') ' Electron induced pion production: The Matrix Element'
       write(*,'(A)')'+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(*,*)
       write(*,'(A)')'###############################################'
       write(*,'(A)') 'Input (lab frame):'
       write(*,*)
       write(*,'(A,4F9.5)') 'electron  incoming :',li_lab
       write(*,'(A,4F9.5)') 'electron  outgoing :',lf_lab
       write(*,'(A,4F9.5)') 'nucleon   incoming :',pi_lab
       write(*,'(A,4F9.5)') 'nucleon   outgoing :',pf_lab
       write(*,'(A,4F9.5)') 'pion      outgoing :',k_lab
       write(*,*)
       write(*,'(A,2F9.5)') 'incoming nucleon mass:',sqrt(pair(pi_lab,pi_lab))
       write(*,'(A,2F9.5)') 'ougoing  nucleon mass:',sqrt(pair(pf_lab,pf_lab))
    end if

    ! Transforming everything to CM-Frame of the hadronic vertex:

    betatoCM = lorentzCalcBeta (q_lab+pi_lab, 'matrixElement')

    li=li_lab
    lf=lf_lab
    q =q_lab
    pi=pi_lab
    pf=pf_lab
    k =k_lab
    call lorentz(betaToCM, li, 'matrixElement')
    call lorentz(betaToCM, lf, 'matrixElement')
    call lorentz(betaToCM, pi, 'matrixElement')
    call lorentz(betaToCM, pf, 'matrixElement')
    call lorentz(betaToCM, q , 'matrixElement')
    call lorentz(betaToCM, k , 'matrixElement')

    if (debug) then
       write(*,'(A)')'###############################################'
       write(*,'(A)') 'Converted to CM frame of hadronic vertex:'
       write(*,*)
       write(*,'(A,4F9.5)') 'electron  incoming :',li
       write(*,'(A,4F9.5)') 'electron  outgoing :',lf
       write(*,'(A,4F9.5)') 'nucleon   incoming :',pi
       write(*,'(A,4F9.5)') 'nucleon   outgoing :',pf
       write(*,'(A,4F9.5)') 'pion      outgoing :',k
       write(*,'(A,4F9.5)') 'virtual   photon   :',q
       write(*,*)
       write(*,'(A,2F9.5)') 'incoming nucleon mass:',sqrt(pair(pi,pi))
       write(*,'(A,2F9.5)') 'ougoing  nucleon mass:',sqrt(pair(pf,pf))
    end if

   matrixElement_eN=contraction(W_free)
  contains

    real function contraction(W_free)
      use minkowski, only: metricTensor
      use leptonTensor, only: l_munu
      use hadronTensor_npi, only: h_munu
      use formfactors_A_main, only: getA

      implicit none
      real, intent(in)       :: W_free
      complex,dimension(1:6) :: A
      integer :: mu,nu
      real :: contraction_c
      real :: theta,s_Vacuum
      real :: s_Vacuum_belowGrid
      logical :: belowGrid


      ! This is meant to solve the problem that acos(1)
      ! can sometimes give NAN since 1 is 1.000000000000001 or so.
      theta=degrees(Acos( &
           &    Min( &
           &    Max(&
           &    Dot_product(q(1:3),k(1:3))&
           &      /sqrt(Dot_product(q(1:3),q(1:3)))/sqrt(Dot_product(k(1:3),k(1:3))) &
           &    ,-0.999999999999999999) &
           &    , 0.999999999999999999)&
           ))


      s_Vacuum=W_free**2

      ! The value 1.096 is dictated by the MAID input grid which has the
      ! lowest point at W=1.1 and dW=0.01
      ! if W below this value, we do a linear extrapolation to the
      ! pion-nucleon-threshold at 1.076 GeV
      belowGrid=.false.
      if (sqrt(s_Vacuum).lt.1.096) then
         belowGrid=.true.
         s_Vacuum_belowGrid=s_Vacuum
         s_Vacuum=1.096**2
      end if


      select case (processID)
      case (EM)

         A=getA(charge_pionOut,charge_nucOut,theta,s_Vacuum,-pair(q,q))

      case (CC)

         select case (charge_pionout)
         case (1)
            A=sqrt(2.)*getA(0,1,theta,s_Vacuum,-pair(q,q))&
                 & -getA(-1,1,theta,s_Vacuum,-pair(q,q))
         case (0)
            A=-getA(0,0,theta,s_Vacuum,-pair(q,q))-  &
              & sqrt(2.)*getA(-1,1,theta,s_Vacuum,-pair(q,q))&
              & +getA(0,1,theta,s_Vacuum,-pair(q,q))
         end select

      case (NC)
         write(*,*) 'MAidlike BG not yet implemented -> stop'
         stop
      end select

      ! extrapolate to lower W if necessary
      if (belowGrid) then
         A=A/0.04344*(s_Vacuum_belowGrid-1.076**2)
         !0.04344=1.096**2-1.076**2, i.e. last Maid gridpoint
         ! minus pion-nucleon threshold
      end if

      contraction_c=0.
      do mu=0,3
         do nu=0,3
            contraction_c=contraction_c+ &
                 & metricTensor(mu,mu)*metricTensor(nu,nu)*l_munu(mu,nu,li,lf)&
                 & *real(h_munu(mu,nu,pi,pf,k,q,A) )
         end do
      end do
      contraction=Contraction_c

    end function contraction

  end function matrixElement_EN



end module electronPionProd_medium_eN
