!******************************************************************************
!****m* /singlePionProductionNHVlike
! NAME
! module singlePionProductionNHVlike
! PURPOSE
! ...
!******************************************************************************
module singlePionProductionNHVlike


  private

  !****************************************************************************
  !****g* singlePionProductionNHVlike/debug_HNV
  ! SOURCE
  integer, parameter ::  debug_HNV=0
  ! PURPOSE
  ! level of input debug information for HNV 1-pion xsec
  ! 0 - no output, 1 - basic messages, 2 - more details
  !****************************************************************************

  !****************************************************************************
  !****g* singlePionProductionNHVlike/integrate_over
  ! SOURCE
  integer, save ::  integrate_over=2      !1=costhetaPi , 2=Epi  3=over cosThetaPi_star_qz in CM frame
  ! PURPOSE
  ! possible values:
  ! * 1 = costhetaPi
  ! * 2 = Epi
  ! * 3 = over cosThetaPi_star_qz in CM frame
  !
  ! which 3-pl differential cross sectio to use for integration:
  ! * 1= dsigma/dcostheta/dElepton/dcosThetaPion was originally used and
  !   works for nuclei.
  !   disadvantage: for some cosThetaPion there are two solutions for Epi,
  !   this  leads to fluctuations on the cross section
  ! * 2= dsigma/dcostheta/dElepton/dEPion  has an advantage, that for a given pion energy there is
  !   only one solution for  the angle between the resonance and pion momenta.
  !   so the integration is simpler and results should be smoother
  !
  ! NOTES
  ! * for 1 : the only option checked for nucleus
  ! * for 2 : code works better and faster,
  !   gives significantly smoother results below Delta peak.
  !   disadvantage: now for the free nucleon only, TO DO : nuclei
  !****************************************************************************


  logical ,save :: pionPot=.true.

  logical, save :: initflag=.true.

  public :: Nieves1piN_elepton_ct_Epi
  public :: Nieves1piN_elepton_ct_ctPi_PhiPi
  public :: Nieves1piN_elepton_ct_ctPi
  public :: Nieves1piN_elepton_ct
  public :: Nieves_MaEl_1pi

contains


  subroutine readinput_HNV()
    use output
    implicit none
    integer :: IOS

    !**************************************************************************
    !****n* singlePionProductionNHVlike/nl_singlePionProductionNHVlike
    ! NAME
    ! NAMELIST /nl_singlePionProductionNHVlike/
    ! PURPOSE
    ! This Namelist includes:
    ! * integrate_over  !1=costhetaPi(must be for nuclei), 2=Epi, 3=over cosThetaPi_star_qz in CM frame
    !**************************************************************************
    NAMELIST /nl_singlePionProductionNHVlike/ integrate_over

    call Write_ReadingInput('nl_singlePionProductionNHVlike',0)
    rewind(5)
    read(5,nml=nl_singlePionProductionNHVlike,IOSTAT=IOS)
    call Write_ReadingInput('nl_singlePionProductionNHVlike',0,IOS)
    write(*,*) 'integrate_over=: ', integrate_over
    call Write_ReadingInput('nl_singlePionProductionNHVlike',1)


  end subroutine readinput_HNV




  !****************************************************************************
  !****s* singlePionProductionNHVlike/Nieves1piN_elepton_ct_Epi
  ! NAME
  ! subroutine Nieves1piN_elepton_ct_Epi(process_ID,k_in,k_out,p_in,position, charge_in,  &
  !     & p_out,charge_out,p_pion_out,pion_charge_out,Epi,sig)
  !
  ! PURPOSE
  ! This subroutine is called by Nieves1piN_elepton_ct,
  ! and does the actual calculation of the 4-fold differential x-sec dSigma/dcostheta dElepton.dcosThetaPion.dPhiPion
  ! It is only used in the case:  include1pi=.true. and singlePiModel=0   (Nieves model for 1-pion x-sec)
  !
  ! INPUTS
  !  integer             :: process_ID   -- CC, NC or EM
  !  real, dimension(0:3):: k_in         -- momentum of incoming lepton
  !  real, dimension(0:3):: k_out        -- momentum of outgoing lepton
  !  real, dimension(0:3):: p_in         -- momentum of incoming hadron
  !  real, dimension(1:3):: position     -- position of incoming hadron
  !  integer             :: charge_in    -- charge of incoming hadron
  !  integer             :: charge_out   -- charge of outgoing nucleon
  !  integer             :: pion_charge_out  -- charge of outgoing pion
  !  real                :: Epi      -- pion energy in lab frame
  !
  !
  ! OUTPUT
  ! real, dimension(0:3) :: p_out        -- momentum of outgoing nucleon (directions)
  ! real, dimension(0:3) :: p_pion_out -- momentum of outgoing pion (directions generated randomly)
  ! real                 :: sig          -- calculated cross section in 1/GeV^3
  !
  ! NOTES
  ! This cross section has an advantage, because for a given pion energy there is only one solution for
  ! the angle between the resonance and pion momenta.
  !
  !****************************************************************************
  subroutine  Nieves1piN_elepton_ct_Epi( process_ID,k_in,k_out,firstTime,p_in,position, charge_in, &
       & p_out,charge_out,p_pion_out,pion_charge_out,Epi,sig)

    use random, only: rn
    use constants, only: pi
    use minkowski, only: SP, abs4
    use particleDefinition
    use degRad_conversion
    use rotation
    use ParticleProperties, only: hadron
    use lepton_kinematics_free, only: CosThetaPiW_W_Epi
    use constants, only: mN, mPi

    implicit none

    integer,     intent(in)          :: process_ID
    real, dimension(0:3), intent(in) :: k_in, k_out
    logical,            intent(in)   :: firstTime   ! guarantees that lepton variables are calculated only once
    real, dimension(0:3), intent(in) :: p_in
    real, dimension(1:3), intent(in) :: position
    integer,              intent(in) :: charge_in
    integer,              intent(in) :: charge_out
    integer,              intent(in) :: pion_charge_out
    real,                 intent(in) :: Epi         ! pion energy in lab frame

    real, dimension(0:3),intent(out) :: p_out, p_pion_out
    real,                intent(out) :: sig

    real    ::  k1lepton2
    real, dimension(0:3) :: q, W
    real    :: thetaW, phiW
    real    :: costhetaPiW, sinthetaPiW, phiPiW, p_pion_abs

    real, dimension(1:3) ::    pion_direction,  pion_direction_relative_to_W

    logical     :: success

    real    :: kinemat_coeff

    if (initFlag) then
       call readinput_HNV()
       initFlag=.false.
    end if

    ! values for the incoming and outgoing leptons
    k1lepton2=Dot_Product(k_out(1:3),k_out(1:3))  ! absolute value of the 3-momentum of outgoing lepton squared
    q=k_in-k_out

    ! direction of the "p_in+q" 3-vector
    W=p_in+q
    ! if below thereshold , reaction not possible, cross section = 0
    if (SP(W,W).lt.(hadron(2)%minmass**2)) then
       sig=0.
       if (debug_HNV.ge.1) write(*,'(2(A,g12.5))') 'invarian mass =', (sqrt(SP(p_in+q,p_in+q))) ,'   is below threshold, sig=0'
       return
    end if
    call get_phi_Theta(W(1:3),thetaW,phiW)

    if (debug_HNV.ge.3) write(*,'((A,4g12.5))') 'In Nieves1piN_elepton_ct_Epi:   W=', W
    if (debug_HNV.ge.3) write(*,'(3(A,g12.5))') 'In Nieves1piN_elepton_ct_Epi:   thetaW=', thetaW,  '    phiW=',phiW


    ! ******************************
    ! solve delta-function for  kinematics
    ! ******************************
    call CosThetaPiW_W_Epi(W,Epi,mN,mpi,costhetaPiW,success)
    if (debug_HNV.ge.2) write(*,'(4(A,g12.5))') 'In Nieves1piN_elepton_ct_Epi:      Epi=',Epi, '     costhetaPiW=', costhetaPiW

    if (.not.success) then
       ! No solution to energy/momentum conservation
       sig=0.
       if (debug_HNV.ge.1) write(*,'(A)') 'No solution found for momentum conservation, sig=0'
       return
    else


       ! ******************************
       ! outgoing kinematics and cross section
       ! ******************************

       !Epi=sqrt(p_pion_abs**2+mpi**2)
       !p_pion_out=(/ Epi , p_pion_abs*sinThetaPi*cos(PhiPi), p_pion_abs*sinThetaPi*sin(PhiPi), p_pion_abs*cosThetaPi /)


       ! define p_pion_out
       sinthetaPiW=sqrt((1.-costhetaPiW)*(1.+costhetaPiW))
       phiPiW=rn()*2.*pi                                       ! this is effectively 1-point MonteCarlo integration

       pion_direction_relative_to_W=(/sinthetaPiW*cos(phiPiW),sinthetaPiW*sin(phiPiW),costhetaPiW/)
       pion_direction = rotateYZ (thetaW, phiW, pion_direction_relative_to_W)
       if (debug_HNV.ge.3) write(*,'(A,3g12.4)') 'pion_direction_W=', pion_direction_relative_to_W
       if (debug_HNV.ge.3) write(*,'(A,g12.4)') '|pion_direction_W|=', &
            & dot_product(pion_direction_relative_to_W,pion_direction_relative_to_W)
       if (debug_HNV.ge.3) write(*,'(A,3g12.4)') 'pion_direction=', pion_direction
       if (debug_HNV.ge.3) write(*,'(A,g12.4)') '|pion_direction|=', dot_product(pion_direction,pion_direction)
       if (debug_HNV.ge.3) write(*,'(A,g12.4)') '|costhetaPiW|=', &
            &  dot_product(W(1:3),pion_direction)/sqrt(dot_product(W(1:3),W(1:3)))


       p_pion_abs=sqrt((Epi-mpi)*(Epi+mpi))
       p_pion_out=(/Epi,p_pion_abs*pion_direction(1:3)/)
       if (debug_HNV.ge.3) write(*,'(A,4g12.4)') 'In Nieves1piN_elepton_ct_Epi:  p_pion_out=', p_pion_out
       if (debug_HNV.ge.3) write(*,'(A,g12.4)') 'In Nieves1piN_elepton_ct_Epi:  sqrt(p_pion_out^2)=', abs4(p_pion_out)

       p_out=W-p_pion_out
       if (debug_HNV.ge.3) write(*,'(A,4g12.4)') 'In Nieves1piN_elepton_ct_Epi:  p_out=', p_out
       if (debug_HNV.ge.3) write(*,'(A,g12.4)') 'In Nieves1piN_elepton_ct_Epi:  sqrt(p_out^2)=', abs4(p_out)

       kinemat_coeff=1./((2.*pi)**3)/32./SP(p_in,k_in)*sqrt(k1lepton2)/sqrt(Dot_Product(W(1:3),W(1:3)))

       if (debug_HNV.ge.2) write(*,'(A,g12.5)') 'In Nieves1piN_elepton_ct_Epi:  kinemat_coeff=', kinemat_coeff

       sig=kinemat_coeff*REAL( Nieves_MaEl_1pi(process_ID, k_in,k_out, p_in, position, &
            & charge_in, p_out,charge_out, p_pion_out,pion_charge_out) )

       if (debug_HNV.ge.2) write(*,'(1(A,g12.5))') 'In Nieves1piN_elepton_ct_ctPi      sig=', sig
    end if

    if (debug_HNV.ge.1)  write(*,*) ' '

  end subroutine  Nieves1piN_elepton_ct_Epi

















  !****************************************************************************
  !****s* singlePionProductionNHVlike/Nieves1piN_elepton_ct_ctPi_PhiPi
  ! NAME
  ! subroutine Nieves1piN_elepton_ct_ctPi_PhiPi(process_ID,k_in,k_out,p_in,position, charge_in,  &
  !     & p_out,charge_out,p_pion_out,pion_charge_out,cosThetaPi,PhiPi,sig)
  !
  ! PURPOSE
  ! This subroutine is called by Nieves1piN_elepton_ct,
  ! and does the actual calculation of the 4-fold differential x-sec dSigma/dcostheta dElepton.dcosThetaPion.dPhiPion
  ! It is only used in the case:  include1pi=.true. and singlePiModel=0   (Nieves model for 1-pion x-sec)
  !
  ! INPUTS
  !  integer             :: process_ID   -- CC, NC or EM
  !  real, dimension(0:3):: k_in         -- momentum of incoming lepton
  !  real, dimension(0:3):: k_out        -- momentum of outgoing lepton
  !  real, dimension(0:3):: p_in         -- momentum of incoming hadron
  !  real, dimension(1:3):: position     -- position of incoming hadron
  !  integer             :: charge_in    -- charge of incoming hadron
  !  integer             :: charge_out   -- charge of outgoing nucleon
  !  integer             :: pion_charge_out  -- charge of outgoing pion
  !  real                :: cosThetaPi   -- cos of the polar angle of the outgoing pion
  !  real                :: PhiPi        -- azimuthal angle of the outgoing pion
  !
  !
  ! OUTPUT
  ! real, dimension(0:3) :: p_out        -- momentum of outgoing nucleon
  ! real, dimension(0:3) :: p_pion_out -- momentum of outgoing pion
  ! real                 :: sig          -- calculated cross section in 1/GeV^3
  !
  !****************************************************************************
  subroutine  Nieves1piN_elepton_ct_ctPi_PhiPi( process_ID,k_in,k_out,firstTime,p_in,position, charge_in, &
       & p_out,charge_out,p_pion_out,pion_charge_out,cosThetaPi,PhiPi,sig)

    use constants, only: pi, mN, mPi
    use minkowski, only: SP, abs4, abs4Sq
    use IdTable, only: nucleon
    use particleDefinition
    use ParticleProperties, only: hadron
    use electronPionProduction_kine, only: getV_out, get_dV_Pi_dk, get_k_abs_improved
    use lepton_kinematics_free, only: lowBoundary_for2costhetaPiSolutions
    use degRad_conversion


    implicit none

    integer, intent(in) :: process_ID
    real, dimension(0:3), intent(in) :: k_in, k_out
    logical, intent(in) :: firstTime ! guarantees that lepton variables are calculated only once
    real, dimension(0:3), intent(in) :: p_in
    real, dimension(1:3), intent(in) :: position
    integer,              intent(in) :: charge_in
    integer,              intent(in) :: charge_out
    integer,              intent(in) :: pion_charge_out
    real,                 intent(in) :: cosThetaPi, PhiPi ! cos of the polar angle of the outgoing pion, azimuthal pion angle

    real, dimension(0:3),intent(out) :: p_out, p_pion_out
    real,                intent(out) :: sig


    real, save  :: k1lepton2
    integer     :: numRoots, r
    real, dimension(0:3) :: q
    real, dimension(1:3) :: pion_direction, betaTOCM
    real, dimension(1:2,0:3) :: p_pion_out_root
    real :: p_pion_abs, p_out_abs, p_out_freeEnergy, p_pion_freeEnergy
    real :: sinThetaPi

    type(particle) :: nucleon_in
    logical :: success


    real :: V_out, dV_out, dV_Pi_dk, kinemat_coeff
    real :: EWmin, numin, Q2min

    if (initFlag) then
       call readinput_HNV()
       initFlag=.false.
    end if


    if (firstTime.eqv..false.) then ! guarantees that lepton masses are calculated only once
    else

       ! values for the incoming and outgoing leptons
       k1lepton2=Dot_Product(k_out(1:3),k_out(1:3))  ! absolute value of the 3-momentum of outgoing lepton squared

       q=k_in-k_out

       ! define unit vector in pion direction
       sinThetaPi=sqrt((1.-cosThetaPi)*(1.+cosThetaPi))
       pion_direction=(/ sinThetaPi*cos(PhiPi), sinThetaPi*sin(PhiPi), cosThetaPi /)

    end if


    ! ******************************
    ! set up incoming nucleon
    ! ******************************
    call setToDefault(nucleon_in)
    nucleon_in%position=position
    nucleon_in%momentum=p_in
    nucleon_in%mass=mN                  ! baremass mass
    !nucleon_in%momentum(0)=freeEnergy(nucleon_in)! E=sqrt(p(1:3)^2+m_0^2)
    nucleon_in%ID=nucleon
    nucleon_in%charge =charge_in
    nucleon_in%antiparticle=.false.
    nucleon_in%perturbative=.true.

    ! ******************************
    ! solve delta-function for  kinematics
    ! ******************************


    betaTOCM(1:3)=(q(1:3)+nucleon_in%momentum(1:3))/(q(0)+nucleon_in%momentum(0))

    if (sqrt(Dot_Product(betaTOCM,betaTOCM))>1.) then
       sig=0
       return
    end if

    call get_k_abs_improved(NumRoots,pion_direction(1),pion_direction(2),pion_direction(3),q,&
         & nucleon_in,pion_charge_out,charge_out,success,p_pion_out_root,betaTOCM)


    if (debug_HNV.ge.2) then
       call lowBoundary_for2costhetaPiSolutions(abs4Sq(p_in+q),p_in(0),mpi,EWmin,numin,Q2min)
       write(*,'(3(A,g12.5))') 'W=', abs4(p_in+q), '    EW=', (p_in(0)+q(0)), '   EWmin=',EWmin
       write(*,'(A,3g12.5,A,g12.5,A,I5)') 'betaTOCM=', betaTOCM, &
            & '    abs(betaTOCM)=', (sqrt(Dot_Product(betaTOCM,betaTOCM)))    ,'   NumRoots=', NumRoots
    end if

    !if (debug_HNV.ge.1) write(*,'(//5(A,g11.5))') ' After getKinematics: Enu=',enu, '    E1=',e1lepton, &
    ! & '    cosThetaPi=',cosThetaPi, '      PhiPi=', PhiPi,   '  p_pion_abs=',p_pion_abs

    !if (debug_HNV.ge.3) write(*,'(//5(A,4g13.5))') ' Check getKinematics:  p_pion_d+p_out_d=',(p_pion_d+p_out_d),  &
    ! & '     sould be equal to q_d+p_in=',(q_d+p_in)



    if (.not.success) then
       ! No solution to energy/momentum conservation
       sig=0.
       if (debug_HNV.ge.1) write(*,'(2(A,g12.5))') 'No solution found for momentum conservation, sig=0'
       return
    else
       ! If W too small, then return:
       if (SP(p_in+q,p_in+q).lt.(hadron(2)%minmass**2)) then
          sig=0.
          if (debug_HNV.ge.1) write(*,'(2(A,g12.5))') 'invarian mass =', (sqrt(SP(p_in+q,p_in+q))) ,'   is below threshold, sig=0'
          return
       end if


       ! ******************************
       ! outgoing kinematics and cross section
       ! ******************************

       !Epi=sqrt(p_pion_abs**2+mpi**2)
       !p_pion_out=(/ Epi , p_pion_abs*sinThetaPi*cos(PhiPi), p_pion_abs*sinThetaPi*sin(PhiPi), p_pion_abs*cosThetaPi /)

       ! summing up over root for absolute pion momentum for a given cosThetaPi
       sig=0
       do r=1,numRoots

          p_pion_out=p_pion_out_root(r,:)
          p_pion_abs=sqrt(Dot_Product(p_pion_out(1:3),p_pion_out(1:3)))
          p_pion_freeEnergy=sqrt(p_pion_abs**2+mpi**2)

          p_out=p_in+q-p_pion_out
          p_out_abs=sqrt(Dot_Product(p_out(1:3),p_out(1:3)))
          p_out_freeEnergy=sqrt(p_out_abs**2+mN**2)

          if (debug_HNV.ge.2) write(*,'(A,g12.5,A,I5,2(A,g12.5))') 'cosThetaPi=', cosThetaPi, &
               & '  Solutions: r=',r, '   p_pion_abs=',p_pion_abs, '   p_out_abs=',p_out_abs


          call getV_out(V_out,dV_out,p_out,nucleon_in,pion_charge_out)
          if (debug_HNV.ge.3) write(*,'(A,2(A,g12.5))') 'In Nieves1piN_elepton_ct_ctPi:  getV_out  has been called, ', &
               & '   V_out=',V_out, '    dV_pit=',dV_out

          if (pionPot) then
             call get_dV_Pi_dk(dV_Pi_dk,p_pion_abs,nucleon_in,pion_charge_out)
             if (debug_HNV.ge.3) write(*,'(A,A,g12.5)') 'In Nieves1piN_elepton_ct_ctPi:  Pion potential has been called', &
                  & '  dV_Pi_dk=', dV_Pi_dk
          else
             dV_Pi_dk=0
          end if

          ! dsigma in units of 1/GEV**3
          !!       if(pionPot) then
          !!         sig=0
          !!       else
          !          Eq. (5.13) in Oliver's thesis
          ! 2008-10-13
          ! extra 1./16. and NO ml^2*mN*mN1  because of my normalization of the fermion spinors
          ! 1/(mN*Enu*relativ_velocity) is replace by INVARIANT 1./SP(p_in,k_in)
          ! 1./(2.*pi)**4   instead of 1./(2.*pi)**5 because  integrations over phi_l is done
          ! e1lepton is replaced with sqrt(k1lepton2), which is correct for outgoing leptons with nonzero masses

          kinemat_coeff=1./((2.*pi)**4)/2./16. &
               & *p_pion_abs**2*sqrt(k1lepton2)/p_pion_out(0)/SP(p_in,k_in)/p_out(0)&
               & /abs( &
               & p_pion_abs/p_pion_freeEnergy + dV_Pi_dk &
               & + (p_pion_abs-Dot_product(p_in(1:3)+q(1:3),p_pion_out(1:3))/p_pion_abs)/p_out(0) &
               &  * (1.+2./p_out_abs*(p_out_abs/p_out_freeEnergy*V_out+p_out_freeEnergy*dV_out+V_out*dV_out)) &
               & )


          if (debug_HNV.ge.3) then
             write(*,'(A,4g12.5,A,g12.5)') 'In Nieves1piN_elepton_ct_ctPi_PhiPi:  k_in=', k_in, '   abs4(k_in)=', abs4(k_in)
             write(*,'(A,4g12.5,A,g12.5)') '                                      p_in=', p_in, '   abs4(p_in)=', abs4(p_in)
             write(*,'(A,4g12.5,A,g12.5)') '                                      k_out=', k_out, '   abs4(k_out)=', abs4(k_out)
             write(*,'(A,4g12.5,A,g12.5)') '                                      p_out=', p_out, '   abs4(p_out)=', abs4(p_out)
             write(*,'(A,4g12.5,A,g12.5)') '                                      p_pion_out=', p_pion_out, &
                  & '   abs4(p_pion_out)=', abs4(p_pion_out)
             write(*,'(5(A,g12.5))') 'p_pion_abs=', p_pion_abs, '   p_out_abs=',p_out_abs, &
                  & '    p_out_freeEnergy=', p_out_freeEnergy, '    p_pion_freeEnergy=', p_pion_freeEnergy,&
                  & '    V_out=',V_out, '   dV_out=',dV_out
          end if

          if (debug_HNV.ge.2) write(*,'(A,g12.5)') 'In Nieves1piN_elepton_ct_ctPi:  kinemat_coeff=', kinemat_coeff

          sig=sig+kinemat_coeff*REAL( Nieves_MaEl_1pi(process_ID, k_in,k_out, p_in, position, &
               & charge_in, p_out,charge_out, p_pion_out,pion_charge_out) )
          !!       end if
       end do !summing up over root for absolute pion momentum for a given cosThetaPi



       if (debug_HNV.ge.1) write(*,'(A,g12.5)') 'In Nieves1piN_elepton_ct_ctPi   sig=', sig
    end if
    if (debug_HNV.ge.1)  write(*,*) ' '

  end subroutine  Nieves1piN_elepton_ct_ctPi_PhiPi



  !****************************************************************************
  !****s* singlePionProductionNHVlike/Nieves1piN_elepton_ct_ctPi
  ! NAME
  ! subroutine Nieves1piN_elepton_ct_ctPi(process_ID,k_in,k_out,p_in,position, charge_in,  &
  !     & p_out,charge_out,p_pion_out,pion_charge_out,cosThetaPi,sig)
  !
  ! PURPOSE
  ! This subroutine  is called by Nieves1piN_elepton_ct
  ! does the calculation of the 3-pl differential x-sec dSigma/dcostheta dElepton.dcosThetaPion
  ! It is only used in the case:  include1pi=.true. and singlePiModel=0   (Nieves model for 1-pion x-sec)
  !
  ! INPUTS
  !  integer             :: process_ID   -- CC, NC or EM
  !  real, dimension(0:3):: k_in         -- momentum of incoming lepton
  !  real, dimension(0:3):: k_out        -- momentum of outgoing lepton
  !  real, dimension(0:3):: p_in         -- momentum of incoming hadron
  !  real, dimension(1:3):: position     -- position of incoming hadron
  !  integer             :: charge_in    -- charge of incoming hadron
  !  integer             :: charge_out   -- charge of outgoing nucleon
  !  integer             :: pion_charge_out  -- charge of outgoing pion
  !  real        :: cosThetaPi   -- cos of the polar angle of the outgoing pion
  !
  !
  ! OUTPUT
  ! real, dimension(0:3) :: p_out        -- momentum of outgoing nucleon
  ! real, dimension(0:3) :: p_pion_out -- momentum of outgoing pioninclude_W_dist=.false.
  ! real                 :: sig          -- calculated cross section in 1/GeV^3
  !
  !****************************************************************************
  subroutine  Nieves1piN_elepton_ct_ctPi( process_ID,k_in,k_out,firstTime,p_in,position, charge_in, &
       & p_out,charge_out,p_pion_out,pion_charge_out,cosThetaPi,sig)
    use constants, only:pi
    use random, only: rn

    implicit none
    integer,              intent(in) :: process_ID
    real, dimension(0:3), intent(in) :: k_in, k_out
    logical,              intent(in) :: firstTime ! guarantees that lepton variables are calculated only once
    real, dimension(0:3), intent(in) :: p_in
    real, dimension(1:3), intent(in) :: position
    integer,              intent(in) :: charge_in
    integer,              intent(in) :: charge_out
    integer,              intent(in) :: pion_charge_out
    real,                 intent(in) :: cosThetaPi ! cos of the polar angle of the outgoing pion

    real, dimension(0:3),intent(out) :: p_out, p_pion_out
    real,                intent(out) :: sig


    real :: PhiPi_int

    if (initFlag) then
       call readinput_HNV()
       initFlag=.false.
    end if

    PhiPi_int= -rn()*2.*pi

    call  Nieves1piN_elepton_ct_ctPi_PhiPi( process_ID,k_in,k_out,firstTime,p_in,position, charge_in, &
         &p_out, charge_out, p_pion_out, pion_charge_out,cosThetaPi,PhiPi_int,sig)
    sig=sig*2*pi

  end subroutine  Nieves1piN_elepton_ct_ctPi









  !****************************************************************************
  !****s* singlePionProductionNHVlike/Nieves1piN_elepton_ct_ctPi_star
  ! NAME
  ! subroutine Nieves1piN_elepton_ct_ctPi_star(process_ID,k_in,k_out,p_in,position, charge_in,  &
  !     & p_out,charge_out_star,p_pion_out_star,pion_charge_out,cosThetaPi_star_qz,sig)
  !
  ! PURPOSE
  ! calculates 3-fold differential cross section dsigma/d E' dcostheta_lepton d costheta_pion_star
  ! for pion angle in the CM frame with respect to the z-axis defined by the photon 3-momentum
  !
  ! Advantage: nothing should depend on the azimuthal angle of pion, so integration is pure 2*pi
  !                 (without assumption of 1-point Monte-Carlo integration)
  !
  ! It is only used in the case:  include1pi=.true. and singlePiModel=0   (Nieves model for 1-pion x-sec)
  !
  ! Useful for CM pion-angle distribution and comparison with electroproduction data like Egiyan PRC 73
  !
  ! INPUTS
  !  integer             :: process_ID   -- CC, NC or EM
  !  real, dimension(0:3):: k_in         -- momentum of incoming lepton
  !  real, dimension(0:3):: k_out        -- momentum of outgoing lepton
  !  real, dimension(0:3):: p_in         -- momentum of incoming hadron
  !  real, dimension(1:3):: position     -- position of incoming hadron
  !  integer             :: charge_in    -- charge of incoming hadron
  !  integer             :: charge_out   -- charge of outgoing nucleon
  !  integer             :: pion_charge_out  -- charge of outgoing pion
  !  real        :: cosThetaPi_star         -- angle of the outgoing pion with respect to 3-momentum of q in the CM frame
  !
  !
  ! OUTPUT
  ! real, dimension(0:3) :: p_out_star_qz        -- momentum of outgoing nucleon in the CM frame
  ! real, dimension(0:3) :: p_pion_star_qz       -- momentum of outgoing pion in the CM frame
  ! real                 :: sig                  -- calculated cross section in 1/GeV^3
  !
  !****************************************************************************
  subroutine  Nieves1piN_elepton_ct_ctPi_star( process_ID, k_in,k_out, firstTime, p_in,position, charge_in, &
       & p_out_star_qz,charge_out,p_pion_star_qz,pion_charge_out,cosThetaPi_star_qz,sig)


    use constants, only: pi, mN, mPi
    use minkowski, only: SP, abs4
    use IdTable, only: nucleon
    use particleDefinition
    use random, only: rn
    use lorentzTrafo, only: lorentz
    use rotation, only: rotateTo

    implicit none

    integer,     intent(in)          :: process_ID
    real, dimension(0:3), intent(in) :: k_in, k_out
    logical,              intent(in) :: firstTime   ! guarantees that lepton variables are calculated only once
    real, dimension(0:3), intent(in) :: p_in
    real, dimension(1:3), intent(in) :: position
    integer,              intent(in) :: charge_in
    integer,              intent(in) :: charge_out
    integer,              intent(in) :: pion_charge_out
    real,                 intent(in) :: cosThetaPi_star_qz
    real, dimension(0:3),intent(out) :: p_out_star_qz, p_pion_star_qz
    real,                intent(out) :: sig


    real, dimension(0:3) :: q, W, qCM,    p_pion, p_out,  p_pion_star, p_out_star
    real, dimension(1:3) :: betaTOCM, pion_unit, pion_unit_rot
    real    :: k1lepton2, W2,  p_pion_abs_star, EN1_star, Epi_star
    real    :: sinThetaPi_star_qz, phi_pion_star

    type(particle)  :: nucleon_in

    real    :: kinemat_coeff

    if (initFlag) then
       call readinput_HNV()
       initFlag=.false.
    end if


    ! ******************************
    ! set up lepton kinematics
    ! ******************************

    k1lepton2=Dot_Product(k_out(1:3),k_out(1:3))  ! absolute value of the 3-momentum of outgoing lepton squared
    q=k_in-k_out
    W=p_in+q
    W2=SP(W,W)

    ! ******************************
    ! set up incoming nucleon
    ! ******************************
    call setToDefault(nucleon_in)
    nucleon_in%position=position
    nucleon_in%momentum=p_in
    nucleon_in%mass=mN
    !nucleon_in%momentum(0)=freeEnergy(nucleon_in)! E=sqrt(p(1:3)^2+m_0^2)
    nucleon_in%ID=nucleon
    nucleon_in%charge =charge_in
    nucleon_in%antiparticle=.false.
    nucleon_in%perturbative=.true.

    ! ******************************
    ! set kinematics in the CM frame
    ! ******************************

    betaTOCM(1:3)=W(1:3)/W(0)

    ! 4-vector of the vector boson in the CM frame
    qCM=q
    call lorentz(betaTOCM, qCM, ' getKinematics_eN')

    ! polar angle of pion
    sinThetaPi_star_qz=sqrt((1.-cosThetaPi_star_qz)*(1.+cosThetaPi_star_qz))
    ! azimuthal angle of pion is choosen randomly
    phi_pion_star=2.*pi*rn()
    ! unit pion vector along z-star-axis
    pion_unit=(/sinThetaPi_star_qz*cos(phi_pion_star), sinThetaPi_star_qz*cos(phi_pion_star), cosThetaPi_star_qz/)
    ! Rotate this pion vector to a system where the CM q is defining the z-axis
    pion_unit_rot = rotateTo (qCM(1:3), pion_unit)


    ! pion energy in the CM frame
    Epi_star=(W2-mN**2+mpi**2)/2./sqrt(W2)
    ! absolute value of pion 3-momentum in the CM frame
    p_pion_abs_star=sqrt((Epi_star-mpi)*(Epi_star+mpi))
    ! outgoing nucleon energy in the CM frame
    EN1_star=sqrt(p_pion_abs_star**2+mN**2)

    ! these two are needed for matrix element calculation
    ! 4-vector of the outgoin pion in the CM frame
    p_pion_star=(/Epi_star,p_pion_abs_star*pion_unit_rot(1:3)/)
    ! 4-vector of the outgoin nucleon in the CM frame
    p_out_star=(/EN1_star,-p_pion_abs_star*pion_unit_rot(1:3)/)


    ! these two are needed for the output
    ! 4-vector of the outgoin pion in the CM frame where q is in z-direction
    p_pion_star_qz=(/Epi_star,p_pion_abs_star*pion_unit(1:3)/)
    ! 4-vector of the outgoin nucleon in the CM frame where q is in z-direction
    p_out_star_qz=(/EN1_star,-p_pion_abs_star*pion_unit(1:3)/)


    if (debug_HNV.ge.2) then
       write(*,'(A,4g12.5,A,g12.5)') 'W=', W, '     abs4(W)=', abs4(W)
       write(*,'(A,3g12.5,A,g12.5)') 'betaTOCM=', betaTOCM, '    abs(betaTOCM)=', (sqrt(Dot_Product(betaTOCM,betaTOCM)))
       write(*,'(A,4g12.5)') 'qCM=', qCM
       write(*,'(A,4g12.5)') 'p_pion_star=', p_pion_star
       write(*,'(A,4g12.5)') 'cos of the angle between qCM and p_pion_star=', &
            & ( Dot_Product(qCM(1:3),p_pion_star(1:3))/p_pion_abs_star/sqrt(Dot_Product(qCM(1:3),qCM(1:3))) )
       write(*,'(A,4g12.5)') '       to be compared with ', cosThetaPi_star_qz
    end if



    ! ******************************
    ! set kinematics in the lab frame   - needed to evaluate the matrix element
    ! ******************************

    ! 4-vector of the outgoin pion in the lab frame
    p_pion=p_pion_star
    call lorentz(-betaTOCM, p_pion, ' Nieves1piN_elepton_ct_star')
    ! 4-vector of the outgoin nucleon in the lab frame
    p_out=p_out_star
    call lorentz(-betaTOCM, p_out,  ' Nieves1piN_elepton_ct_star')



    ! ******************************
    ! cross section
    ! ******************************

    kinemat_coeff=1./((2.*pi)**3)/32./mN/k_in(0)*sqrt(k1lepton2)*p_pion_abs_star/sqrt(W2)

    sig=kinemat_coeff*REAL( Nieves_MaEl_1pi(process_ID, k_in,k_out, p_in, position, &
         & charge_in, p_out,charge_out, p_pion,pion_charge_out) )

    if (debug_HNV.ge.3) write(*,'(A,g12.5)') 'In Nieves1piN_elepton_ct_ctPi   sig=', sig

  end subroutine  Nieves1piN_elepton_ct_ctPi_star





  !****************************************************************************
  !****s* singlePionProductionNHVlike/Nieves1piN_elepton_ct
  ! NAME
  ! subroutine Nieves1piN_elepton_ct(process_ID,k_in,k_out,firstTime,p_in,position, charge_in,  &
  !     & p_out,charge_out,p_pion_out,pion_charge_out,sig)
  !
  ! PURPOSE
  ! This subroutine is called by XsecdCosthetadElepton,
  ! and does the calculation of the 2-bl differential x-sec dSigma/dcostheta dElepton.
  ! It is only used in the case:  include1pi=.true. and singlePiModel=0   (Nieves model for 1-pion x-sec)
  !
  ! INPUTS
  !  integer             :: process_ID   -- CC, NC or EM
  !  real, dimension(0:3):: k_in         -- momentum of incoming lepton
  !  real, dimension(0:3):: k_out        -- momentum of outgoing lepton
  !  logical,        :: firstTime    -- if true than lepton variables are calculated only once
  !  real, dimension(0:3):: p_in         -- momentum of incoming hadron
  !  real, dimension(1:3):: position     -- position of incoming hadron
  !  integer             :: charge_in    -- charge of incoming hadron
  !  integer             :: charge_out   -- charge of outgoing nucleon
  !  integer             :: pion_charge_out  -- charge of outgoing pion
  !
  ! OUTPUT
  ! real, dimension(0:3) :: p_out        -- momentum of outgoing nucleon, angles being generated randomly
  ! real, dimension(0:3) :: p_pion_out -- momentum of outgoing pion, angles being generated randomly
  ! real                 :: sig          -- calculated cross section in 1/GeV^3
  !
  !****************************************************************************
  subroutine  Nieves1piN_elepton_ct( process_ID, k_in,k_out, firstTime, p_in,position, charge_in, &
       & p_out,charge_out,p_pion_out,pion_charge_out,sig)

    use random, only: rn
    use lepton_kinematics_free, only: minmaxEpi_W
    use constants, only: mN, mPi

    implicit none

    integer,              intent(in) :: process_ID
    real, dimension(0:3), intent(in) :: k_in, k_out
    logical,              intent(in)   :: firstTime   ! guarantees that lepton variables are calculated only once
    real, dimension(0:3), intent(in) :: p_in
    real, dimension(1:3), intent(in) :: position
    integer,              intent(in) :: charge_in
    integer,              intent(in) :: charge_out
    integer,              intent(in) :: pion_charge_out

    real, dimension(0:3),intent(out) :: p_out, p_pion_out
    real,                intent(out) :: sig

    real :: cosThetaPi_int
    real    :: Epi_int, Epi_min, Epi_max
    real, dimension(0:3)    :: W
    logical :: success_kinemLimits

    if (initFlag) then
       call readinput_HNV()
       initFlag=.false.
    end if


    if (integrate_over.eq.1) then

       cosThetaPi_int= -1. + rn()*2.
       call  Nieves1piN_elepton_ct_ctPi( process_ID, k_in,k_out,firstTime,p_in,position, charge_in, &
            &p_out, charge_out, p_pion_out, pion_charge_out,cosThetaPi_int,sig)
       sig=sig*2.


       ! checked for free nucleon only, should be modified for the bound nucleon
    else if (integrate_over.eq.3) then
       cosThetaPi_int= -1. + rn()*2.
       ! the pion angle in this case is in CM frame where q is along z-axis
       ! the outgoing 4-momenta p_out and p_pion_out are in the same frame
       call Nieves1piN_elepton_ct_ctPi_star( process_ID, k_in,k_out, firstTime, p_in,position, charge_in, &
            & p_out,charge_out,p_pion_out,pion_charge_out,cosThetaPi_int,sig)
       sig=sig*2.

       ! checked for free nucleon only, should be modified for the bound nucleon
    else if (integrate_over.eq.2) then

       W=p_in+k_in-k_out
       call minmaxEpi_W(W,mN,mpi,Epi_min,Epi_max,success_kinemLimits)

       if (.not.success_kinemLimits) then
          sig=0
          p_out=(/mN,0.,0.,0./)
          p_pion_out=(/mpi,0.,0.,0./)
          return
       end if

       Epi_int=Epi_min+rn()*(Epi_max-Epi_min)
       call  Nieves1piN_elepton_ct_Epi( process_ID, k_in,k_out,firstTime,p_in,position, charge_in, &
            & p_out, charge_out, p_pion_out, pion_charge_out,Epi_int,sig)

       sig=sig*(Epi_max-Epi_min)

    else
       stop 'STOP. Integrate_over should be either 1 or 2 or 3'
    end if


  end subroutine  Nieves1piN_elepton_ct




  !****************************************************************************
  !****f* singlePionProductionNHVlike/Nieves_MaEl_1pi
  ! NAME
  ! complex function  Nieves_MaEl_1pi(process_ID,  k_in, k_out,  p_in, position, charge_in,
  ! p_out, charge_out,  ppi_out, pion_charge_out)
  !
  ! PURPOSE
  ! This function is called by Nieves1piN_elepton_ct_ctPi,
  ! and does the actual calculation of the matrix element of 1-pion production according to the
  ! Hernandez-Nieves-Valverde model.
  ! It is only used in the case:  include1pi=.true. and singlePiModel=0   (Nieves model for 1-pion x-sec)
  !
  ! INPUTS
  !  integer             :: process_ID   -- CC, NC or EM
  !  real, dimension(0:3):: k_in         -- momentum of incoming lepton
  !  real, dimension(0:3):: k_out        -- momentum of outgoing lepton
  !  real, dimension(0:3):: p_in         -- momentum of incoming hadron
  !  real, dimension(1:3):: position     -- position of incoming hadron
  !  integer             :: charge_in    -- charge of incoming hadron
  !  real, dimension(0:3):: p_out        -- momentum of outgoing nucleon
  !  integer             :: charge_out   -- charge of outgoing nucleon
  !  real, dimension(0:3):: ppi_out        -- momentum of outgoing pion
  !  integer             :: pion_charge_out  -- charge of outgoing pion
  !
  ! OUTPUT
  ! Matrix element as a complex variable. Matrix element is real by definition, so at the end of the
  ! function the spetial warning is given if the imaginary part exceeds 10^{-6}
  !
  !****************************************************************************
  complex function  Nieves_MaEl_1pi(process_ID,  k_in, k_out,  p_in, position, charge_in, &
       & p_out, charge_out,  ppi_out, pion_charge_out)

    use leptonTensor
    use NievesHadronTensor
    use minkowski, only: Contract

    implicit none

    integer, intent(in) :: process_ID !CC, NC, em
    real, dimension(0:3), intent(in) :: k_in, k_out, p_in, p_out, ppi_out
    integer, intent(in) :: charge_in, charge_out, pion_charge_out
    real, dimension(1:3), intent(in) :: position



    complex, dimension(0:3,0:3) :: leptonTens, hadronTens
    logical :: success_hadronTens
    real, dimension(0:3) :: q

    if (initFlag) then
       call readinput_HNV()
       initFlag=.false.
    end if

    ! leptonic tensor
    ! factors with GF or alpha_QED sit here
    leptonTens=leptonicTensor(process_ID,k_in,k_out)
    if (debug_HNV.ge.3) write(*,'(A,g12.4,g12.4)') 'In Nieves_MaEl_1pi:   leptonTens is calculated '

    ! hadronic tensor and contraction
    q=k_in-k_out
    success_hadronTens = NievesHadronTensor_1pi(process_ID,q,p_in, position, charge_in, &
         & p_out,charge_out,ppi_out,pion_charge_out,hadronTens)
    if (success_hadronTens) then
       if (debug_HNV.ge.3) write(*,'(A,g12.4,g12.4)') 'In Nieves_MaEl_1pi:   hadronTens  is calculated'

       Nieves_MaEl_1pi=Contract(hadronTens,leptonTens)
       if (debug_HNV.ge.1) write(*,'(A,g12.4,g12.4)') 'In Nieves_MaEl_1pi:    (Matrix element)^2  is :', Nieves_MaEl_1pi

    else
       Nieves_MaEl_1pi=0
       if (debug_HNV.ge.2) write(*,'(A)') 'In Nieves_MaEl_1pi:  matrix element is put to 0'
    end if

    ! check that matrix element has no imaginary part
    if (abs(AIMAG(Nieves_MaEl_1pi)).gt.0.0000001) then
       write(*,'(A,2g12.5)') 'ERROR: In Nieves_MaEl_1pi:    (Matrix element)^2  is imaginary: ', Nieves_MaEl_1pi
       stop
    end if

  end function Nieves_MaEl_1pi






end module singlePionProductionNHVlike
