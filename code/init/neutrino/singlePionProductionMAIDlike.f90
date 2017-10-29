!******************************************************************************
!****m* /singlePionProductionMAIDlike
! NAME
! module singlePionProductionMAIDlike
! PURPOSE
! ...
!******************************************************************************
module singlePionProductionMAIDlike

  implicit none
  private

  !****************************************************************************
  !****g* singlePionProductionMAIDlike/elpionData
  ! SOURCE
  !character(len=3) :: elpionData='BNL'
  ! PURPOSE
  ! Switch to use fit to ANL or BNL elementary pion production data
  ! Use: elpionData = ANL or elpionData = BNL
  ! switch used to set parameters for 1 pion BG and for Delta resonance:
  !real,save :: b_proton_pinull, b_neutron_piplus, MADelta
  !
  !****************************************************************************

  !****************************************************************************
  !****g* singlePionProductionMAIDlike/b_proton_pinull
  ! SOURCE
  real, save :: b_proton_pinull=3.
  ! PURPOSE
  ! 3. is  tuned to ANL, 6. is tuned to BNL
  !****************************************************************************

  !****************************************************************************
  !****g* singlePionProductionMAIDlike/b_neutron_piplus
  ! SOURCE
  real, save :: b_neutron_piplus=1.5
  ! PURPOSE
  ! 1.5 is tuned to ANL, 3. is tuned to BNL
  !****************************************************************************

  logical, save :: initFlag=.true.

  public :: MAIDlike_singlePi

contains

 subroutine readInput
    use output

    integer :: ios
    !**************************************************************************
    !****n* singlePionProductionMAIDlike/neutrino_MAIDlikeBG
    ! NAME
    ! NAMELIST neutrino_MAIDlikeBG
    ! PURPOSE
    ! Includes parameters:
    ! * b_proton_pinull
    ! * b_neutron_piplus
    !**************************************************************************
    NAMELIST /neutrino_MAIDlikeBG/  b_proton_pinull, b_neutron_piplus
    call Write_ReadingInput('neutrino_MAIDlikeBG',0)
    rewind(5)
    read(5,nml=neutrino_MAIDlikeBG,IOSTAT=ios)
    call Write_ReadingInput("neutrino_MAIDlikeBG",0,ios)
    write(*,'(a,2F12.4)') ' b_proton_pinull, b_neutron_piplus:', b_proton_pinull, b_neutron_piplus
    call Write_ReadingInput('neutrino_MAIDlikeBG',1)
 end subroutine readInput


  !****************************************************************************
  !****f* singlePionProductionMAIDlike/MAIDlike_singlePi
  ! NAME
  ! real function MAIDlike_singlePi(eNev, charge_out,pion_charge_out,
  ! p_out,pion_momentum_out,nuclear_phasespace)
  ! PURPOSE
  ! ...
  ! INPUTS
  ! ...
  ! OUTPUT
  ! ...
  !****************************************************************************
  real function MAIDlike_singlePi(eNev, charge_out,pion_charge_out,  &
       & p_out,pion_momentum_out,nuclear_phasespace)

    use electronPionProd_medium_EN
    use particleDefinition
    use IdTable, only: nucleon
    use random, only: rn, rnCos
    use constants, only: alphaQED, coscab, GF, mN, mPi, pi, hbarc
    use degRad_conversion, only: degrees
    use resProd_lepton, only: dSdO_fdE_fdO_k_med_res_EN
    use minkowski, only: SP
    use leptonicID
    use eN_eventDefinition, only: electronNucleon_event
    use eN_event, only: eNev_init_Target, nuclearFluxFactor_correction


    type(electronNucleon_event),intent(in) :: eNev
    integer,              intent(in) :: charge_out
    integer,              intent(in) :: pion_charge_out

    logical,              intent(in) :: nuclear_phasespace

    real, dimension(0:3), intent(out) :: p_out
    real, dimension(0:3), intent(out) :: pion_momentum_out

    type(electronNucleon_event) :: eN
    real ,dimension(0:3) :: pf,k  !pf --- outgoing nucleon, k -- pion momentum
    real :: sigmaPion
    real,dimension(-1:1)  :: sigmaRes
    real :: phi_k_MC, theta_k_MC!, theta_lf_MC
    real :: deltaMonteCarlo
    type(particle) :: targetNuc
    logical:: success, twoRoots
    real :: sigmaVEC, sigmaAXint
    real :: coupling, Qs
    integer :: process_ID
    real, dimension(0:3) :: k_in, p_in, k_out


    if (initFlag) then
       call readInput
       initFlag=.false.
    end if



    ! set some abbreviations:
    process_ID = eNev%idProcess
    k_in  = eNev%lepton_in%momentum
    k_out = eNev%lepton_out%momentum
    p_in  = eNev%nucleon%momentum
    ! ---------------

    MAIDlike_singlePi=0.
    p_out=0.
    pion_momentum_out=0.
    sigmaAXInt=0.
    coupling=1.

    if (process_ID.eq.CC.and.((pion_charge_out+charge_out).eq.2)) return  !no background for isopin3/2 channel in CC
    if (process_ID.le.0) return !no background for anti-leptons
    if (process_ID.eq.NC) return !no background for NC

    ! Check for Monte-Carlo-Integrations:
    deltaMonteCarlo=4.*pi*2.*pi   !4pi from pion angles, 2pi from lepton_phi

    theta_k_MC=degrees(rnCos())
    phi_k_MC=degrees(rn()*2*pi)


    ! Construct particle which is free and has the 3-momentum of the considered nucleon
    call setToDefault(targetNuc)
    targetNuc%momentum=p_in
    targetNuc%charge=eNev%nucleon%charge
    targetNuc%id=nucleon
    targetNuc%mass=mN        ! m=m_0
    targetNuc%momentum(0)=freeEnergy(targetNuc)! E=sqrt(p(1:3)^2+m_0^2)
    targetNuc%position(:)=1000.                ! Far away the nucleus
    targetNuc%perturbative=.true.

    eN = eNev ! make a copy with a free nucleon
    call eNev_init_Target(eN,targetNuc,success)


    !threshold cut
    if (eN%W_free.le.(mN+mPi)) return !less than threshold - reaction not possible

    !cuts required by MAID due to limited grid size
    if ((eN%W_free.ge.2.).or.(-SP(k_in-k_out,k_in-k_out).ge.4.9)) return


    ! Full pion production cross section:
    sigmaPion=dSdO_fdE_fdO_k_med_eN(eN, pion_Charge_out, phi_k_MC, theta_k_MC,&
                                   & k, pf, process_ID,pionNucleonSystem=2)
    ! the last parameter pionNucleonSystem=2 is necessary to correctly evaluate
    ! the nucleon structure functions in the cm frame, pointed out by Leitner,
    ! not in the lab frame, as used by Lalakulich and Paschos

    if (nuclear_phasespace) sigmaPion=sigmaPion*nuclearFluxFactor_correction(targetNuc%momentum &
         & ,k_out+k+pf-targetNuc%momentum)

    if (sigmaPion.lt.1E-20) return ! no solution to kinematical problem was found


    ! Resonance contribution:
    sigmaRes=dSdO_fdE_fdO_k_med_res_EN(eN,k,pf,process_ID,pionNucleonSystem=2)


    ! Subtract resonance contributions:
    sigmaVEC=sigmaPion-sigmares(pion_charge_out)
    sigmaVEC=sigmaVEC*0.1/hbarc**2


    if (process_ID.eq.CC.and.charge_out.eq.neutron) then
       sigmaAXInt=b_neutron_piplus*sigmaVEC
    else if (process_ID.eq.CC.and.charge_out.eq.proton) then
       sigmaAXInt=b_proton_pinull*sigmaVEC
    end if


    !rescale coupling for CC
    if (process_ID.eq.CC) then
       Qs=-SP(k_in-k_out,k_in-k_out)
       coupling=GF**2/2.*coscab**2/((4.*pi*alphaQED)**2/Qs**2/2.)
    end if

    MAIDlike_singlePi=(sigmaVEC+sigmaAXInt)*deltaMonteCarlo*coupling



    ! Evaluate kinematics in the medium:
    call setToDefault(targetNuc)
    targetNuc%momentum=p_in
    targetNuc%charge=eNev%nucleon%charge
    targetNuc%id=nucleon
    targetNuc%mass=mN
    targetNuc%position=eNev%nucleon%position
    targetNuc%perturbative=.true.

    call eNev_init_Target(eN,targetNuc,success)

    ! Evaluate kinematics in the medium:
    call getKinematics_eN(eN,pion_Charge_out,charge_out,phi_k_MC,theta_k_MC,k, &
                         & pf,twoRoots,success,pionNucleonSystem=2)
    ! the last parameter pionNucleonSystem=2 is necessary to correctly evaluate
    ! the nucleon structure functions in the cm frame, pointed out by Leitner,
    ! not in the lab frame, as used by Lalakulich and Paschos

    if (.not.success) then
       ! Kinematics could not be established in the medium
       MAIDlike_singlePi=0.0
       return
    else if (tworoots) then
       write(*,*) 'Error in getKinematics_en! Two solutions!!'
       stop 'singlePionProductionMAIDlike.f90'
    end if

    p_out=pf
    pion_momentum_out=k


  end function MAIDlike_singlePi


end module singlePionProductionMAIDlike
