!******************************************************************************
!****m* /quasiElastic_electron
! NAME
! module quasiElastic_electron
! PURPOSE
! * Evaluates cross sections for gamma N -> N'
! * For details see the notes about this in the work of Oliver Buss
!******************************************************************************
module quasiElastic_electron
  implicit none

  private

  public :: dSigmadcosTheta_l_dE_l_BW_eN
  public :: matrixElement

  logical, save :: initFlag=.true.

  !****************************************************************************
  !****g* quasiElastic_electron/vac_kinematics
  ! SOURCE
  !
  logical,save :: vac_kinematics=.false.
  !
  ! PURPOSE
  ! Evaluates the cross section with matrix element in medium, but kinematical
  ! stuff not.
  !****************************************************************************

  !****************************************************************************
  !****g* quasiElastic_electron/benharMethod
  ! SOURCE
  !
  logical,save :: benharMethod=.false.
  !
  ! PURPOSE
  ! Evaluates the in-medium cross sections as Benhar does.
  !****************************************************************************

  !****************************************************************************
  !****g* quasiElastic_electron/benharMethod_sim
  ! SOURCE
  !
  logical,save :: benharMethod_sim=.false.
  !
  ! PURPOSE
  ! Evaluates the in-medium cross sections similar as Benhar does.
  !****************************************************************************


  !****************************************************************************
  !****g* quasiElastic_electron/final_free
  ! SOURCE
  !
  logical,save :: final_free=.false.
  !
  ! PURPOSE
  ! Evaluates the cross section assuming the final particle is free.
  !
  ! not used up to now
  !****************************************************************************

  !****************************************************************************
  !****g* quasiElastic_electron/simple_BW
  ! SOURCE
  !
  logical,save :: simple_BW=.false.
  !
  ! PURPOSE
  ! * If .true. then we use for the nucleon a simple Breit-Wigner with the
  ! width defined in "width".
  ! * If .false. then we use the module "spectralFunc".
  !****************************************************************************

  !****************************************************************************
  !****g* quasiElastic_electron/width
  ! SOURCE
  !
  real,save :: width=0.001
  !
  ! PURPOSE
  ! Width of the BW used for the QE cross section.
  !****************************************************************************





contains

  !****************************************************************************
  !****s* quasiElastic_electron/initInput
  ! NAME
  ! subroutine initInput
  !
  ! PURPOSE
  ! Reads in job card.
  !****************************************************************************
  subroutine initInput
    use output

    integer :: ios

    NAMELIST /quasiElastic_elec/ benharMethod,benharMethod_sim,vac_kinematics,final_free,width,simple_BW

    call Write_ReadingInput('quasiElastic_elec',0)

    rewind(5)
    read(5,nml=quasiElastic_elec,IOSTAT=ios)
    call Write_ReadingInput("quasiElastic_elec",0,ios)

    write(*,*) 'Use Benhar method for in-medium QE-scattering?        ', benharMethod
    write(*,*) 'Use method similar Benhar for in-medium QE-scattering?', benharMethod_sim
    write(*,*) 'Use method where final nucleon is assumed free?       ', final_free
    if (benharMethod.and.benharMethod_sim.and.final_free) then
       write(*,*) " Use either Benhard's method or a similar one, but not two at the same time! STOP"
       stop
    end if
    write(*,*) 'Vaccuum kinematics in the cross section pre factor?   ', vac_kinematics

    write(*,*)
    write(*,*) 'Use the real spectral function=',.not.simple_BW
    if (simple_BW) write(*,*) '  => Width of QE-peak smearing=',width

    call Write_ReadingInput('quasiElastic_elec',1)

  end subroutine initInput


  !****************************************************************************
  !****f* quasiElastic_electron/dSigmadcosTheta_l_dE_l_BW_eN
  ! NAME
  ! function dSigmadcosTheta_l_dE_l_BW_eN(eN,pf,nuc_bareMass) result(sigma)
  !
  ! PURPOSE
  ! * Evaluates cross sections for e N -> e' N'
  ! * For details see the notes about this in the work of Oliver Buss
  ! * Evaluates dSigma/dcos(theta_lepton)/dE_lepton by replacing the delta
  !   function delta(p'^2-m'^2) by a Breit-Wigner function or a real spectral
  !   function
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN  -- underlying electron nucleon event
  !
  ! OUTPUT
  ! * real, dimension(0:3) :: pf   -- outgoing nucleon 4-momenta
  ! * real                 :: nuc_bareMass        --  bare mass of nucleon
  ! * The cross section "dSigmadcosTheta_l_dE_l_BW_eN" is given in units
  !   of mb/GeV=mcb/MeV
  !
  !****************************************************************************
  function dSigmadcosTheta_l_dE_l_BW_eN(eN,pf,nuc_bareMass,sigmaPS) result(sigma)
    use particleDefinition
    use Idtable, only: nucleon
    use minkowski, only: SP, abs3,abs4
    use constants, only: twopi,mN,hbarc
    use degRad_conversion
    use distributions, only: bw
    use spectralFunc, only: specFunc
    use eN_eventDefinition, only: electronNucleon_event
    use potentialModule, only: potential_LRF, scapot

    type(electronNucleon_event) , intent(in) :: eN
    real, dimension(0:3),intent(out)         :: pf
    real, intent(out)                        :: nuc_bareMass
    real, intent(out), optional              :: sigmaPS

    ! local variables
    real    :: sigma
    type(particle) :: targetNuc
    real, dimension(0:3) :: lin,pin,lf
    real                 :: pi_vec_abs  !, pf_vec_abs
    real :: dummy,mi,mf
    real :: kinematics
    integer, save :: zerocount
    logical :: flagOk

    if (initFlag) then
       call initInput()
       zerocount=0
       initFlag=.false.
    end if

    ! Initial lepton: Assume lepton in z-direction
    lin=eN%lepton_in%momentum

    ! Initial nucleon
    pin=eN%nucleon%momentum
    targetNuc=en%nucleon

    if (benharMethod.or.benharMethod_sim) then
       pin(0)=freeEnergy(targetNuc)
    end if
    mi        = abs4(pin)
    pi_vec_abs= abs3(pin)

!    write(*,*) 'mi = ',mi,pi_vec_abs
!    write(*,*) '.. = ',mi-scalarPotential_nucleon(pi_vec_abs,targetNuc%charge,targetNuc%position),targetNuc%charge

    ! Final lepton
    lf=en%lepton_out%momentum

    ! Final Nucleon: Momentum conservation
    pf=pin+lin-lf
!     pf_vec_abs=abs3(pf)

    mf=abs4(pf,flagOk)
    if (.not.flagOk) then
       zerocount=zerocount +1
       write(*,*) ' Negative baryon pf^mu pf_mu in dSigmadcosTheta_l_dE_l_BW '
       write(*,*) '   counter = ', zerocount
       sigma=0.
       return
    end if

    ! Evaluate the bare mass of the nucleon:
    dummy=scapot(1,targetNuc%charge,pf,targetNuc%position,nuc_bareMass)

    if (benharMethod.or.benharMethod_sim) then
       nuc_bareMass=mf
       pf(0)=pin(0)+lin(0)-lf(0)+potential_LRF(1,targetNuc%charge,pin,targetNuc%position)
       if (benharMethod_sim) &
            & pf(0)=pf(0)-potential_LRF(1,targetNuc%charge,pf,targetNuc%position)
    end if


!    write(15,*)nuc_bareMass,bw(nuc_bareMass,mN,width)

    if (vac_kinematics) then
       pi_vec_abs=sqrt(Dot_product(pin(1:3),pin(1:3)))
       kinematics = mN*lf(0)/(twopi*SP(lin,(/sqrt(mN**2+pi_vec_abs**2),pin(1:3)/)))
       write(99,*) kinematics, mi*lf(0)/(twopi*SP(lin,pin))
    else
       kinematics = mi*lf(0)/(twopi*SP(lin,pin))
    end if


    if (simple_bw) then
       sigma=kinematics &
            & * matrixElement(pin,pf,lin,lf,targetNuc%charge)* bw(nuc_bareMass,mN,width)
       if (present(sigmaPS)) sigmaPS = kinematics*bw(nuc_bareMass,mN,width)*hbarc**2*10.
    else
       sigma=kinematics &
            & * matrixElement(pin,pf,lin,lf,targetNuc%charge) &
            & * 2*mf*specFunc(nucleon,targetNuc%charge,pf,targetNuc%position,nuc_bareMass)
       if (present(sigmaPS)) &
            sigmaPS = kinematics&
            & *2*mf*specFunc(nucleon,targetNuc%charge,pf,targetNuc%position,nuc_bareMass) &
            & *hbarc**2*10.
    end if

    ! Converting to units
    ! 1/GeV**2=1/1000**2/MeV**2=1/1000**2/(1/197 fm**2)=(197/1000)**2 fm**2= (197/1000)**2 * 10 mb

    sigma=sigma*hbarc**2*10. ! Now the cross section is given in units of mb/GeV

  end function dSigmadcosTheta_l_dE_l_BW_eN



  !****************************************************************************
  !****f* quasiElastic_electron/matrixElement
  ! NAME
  ! real function matrixElement(pin,pf,lin,lf,nuc_charge)
  !
  ! PURPOSE
  ! * Evaluates matrix element for e N -> e' N'
  ! * For details see the notes about this in the work of Oliver Buss
  !
  ! INPUTS
  ! * integer              :: nucCharge   -- nucleon charge
  ! * real, dimension(0:3) :: pin,pf      -- incoming and outgoing nucleon 4-momentum
  ! * real, dimension(0:3) :: lin,lf      -- incoming and outgoing lepton 4-momentum
  !
  ! NOTES
  ! Note that the lepton tensor does not include the factor of m_e^2,
  ! therefore it is not included in this matrixelement.
  !****************************************************************************
  real function matrixElement(pin,pf,lin,lf,nuc_charge)
    use minkowski
    use hadronTensor_QE, only: H_munu_QE
    use leptonTensor, only: l_munu

    real, dimension(0:3), intent(in) :: pin,pf,lin,lf
    integer             , intent(in) :: nuc_Charge
    integer :: mu,nu

    matrixElement=0
    do mu=0,3
       do nu=0,3
          matrixElement=matrixElement &
               & + metricTensor(mu,mu)*metricTensor(nu,nu)&
               &   * l_munu(mu,nu,lin,lf)*H_munu_QE(mu,nu,pin,pf,nuc_charge)
       end do
    end do
  end function matrixElement

end module quasiElastic_electron
