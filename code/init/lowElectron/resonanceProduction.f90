!******************************************************************************
!****m* /resProd_lepton
! NAME
! module resProd_lepton
! PURPOSE
! * Evaluates cross sections for gamma N -> R
! * For details see the notes about this in the work of Oliver Buss
!******************************************************************************
module resProd_lepton
  implicit none
  private

  public :: dSigmadOmega_fdE_f_resProd_eN
  public :: dSdO_fdE_fdO_k_med_res_eN
  public :: sigma_resProd
  public :: dSdOmega_k_med_res
!  public :: sigma_pipi_res
  public :: sigma_pipi_res_vac, sigma_barMes_res_vac

  !****************************************************************************
  !****g* resProd_lepton/debug
  ! SOURCE
  logical, parameter :: debug=.false.
  ! PURPOSE
  !
  !****************************************************************************

  logical, save :: initFlag=.true.

contains

  !****************************************************************************
  !****f* resProd_lepton/dSigmadOmega_fdE_f_resProd_eN
  ! NAME
  ! function  dSigmadOmega_fdE_f_resProd_eN(eN,resID,pout, bareMass,processID) result(xSection)
  !
  ! PURPOSE
  ! * Evaluates cross sections for e N -> R
  ! * For details see the notes about this in the work of Oliver Buss
  !
  ! INPUTS
  ! * type(electronNucleon_event)  :: eN          -- underlying electron nucleon event
  ! * integer                      :: resID       -- resonance ID
  ! * integer, optional            :: processID
  !
  ! OUTPUT
  ! * real                         :: bareMass    -- bare resonance mass
  ! * real, dimension(0:3)         :: pout        -- resonance momentum
  ! * real :: xSection     -- dsigma/dOmega/dE in units of mb/(GeV sr)
  !
  ! NOTES
  ! * Enhances dSigmadOmega_fdE_f_resProd by allowing arbitrary electron
  !   momentum directions.
  !****************************************************************************
  function dSigmadOmega_fdE_f_resProd_eN(eN,resID,pout, bareMass,processID) result(xSection)
    use leptonicID, only: EM, CC
    use spectralFunc, only: specFunc
    use constants, only: GeVSquared_times_mb, pi
    use minkowski, only: contract,SP
    use hadronTensor_ResProd
    use leptonTensor
    use degRad_conversion
    use particleDefinition
    use eN_eventDefinition, only: electronNucleon_event

    type(electronNucleon_event) , intent(in)  :: eN
    integer                     , intent(in)  :: resID
    integer, optional           , intent(in)  :: processID
    real                        , intent(out) :: bareMass
    real, dimension(0:3)        , intent(out) :: pout

    ! local variables
    integer                    :: process_ID
    real, dimension(0:3)       :: pin, lin,q,lout
    complex,dimension(0:3,0:3) :: hadronTensor,leptonTens
    complex                    :: matrixElement_Squared
    real                       :: kinematics, Xsection
    real                       :: spec

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    ! ********** Kinematics ******************

    ! Incoming nucleon
    pin=eN%nucleon%momentum

    ! Initial lepton: Assume lepton in z-direction
    lin=eN%lepton_in%momentum

    !Outgoing lepton
    lout=eN%lepton_out%momentum

    q=lin-lout
    pout=pin+q

    ! ********** Cross section ******************

    process_ID=EM

    if (present(processID)) then
       process_ID=processID
       if (processID.eq.CC) process_ID=99   !needed in formfactors
    end if

    leptonTens=l_munu_matrix(lin,lout)
    spec=specFunc(resID,en%nucleon%Charge,pout,en%nucleon%position,baremass)
    if (hadronTensor_R(pin,pout,resID,en%nucleon%Charge,process_ID,hadronTensor,baremass) ) then
       matrixElement_Squared=Contract( leptonTens, hadronTensor)
    else
       matrixElement_Squared=0.
       xsection=0.
       return
    end if

    if (abs(AIMAG(matrixElement_Squared)).gt.0.0000001) then
       write(*,*) 'ERROR: (Matrix element)^2  is imaginary:', matrixElement_Squared
    end if


    kinematics=1./SP(lin,pin)*sqrt(Dot_Product(lout(1:3),lout(1:3)))/(32.*pi**2)

    Xsection=kinematics  *REAL(matrixElement_Squared)*spec


    ! ********** UNIT CONVERSION ******************

    ! Convert from 1/(sr GeV^3) to mb/(GeV sr)
    ! 1/(GeV^2)=1/(1000 MeV)^2=1/(1000/(197 fm)^2) =0.197**2 *10 mb
    Xsection=Xsection/GeVSquared_times_mb

    ! Convert from mb/(GeV sr) to mb/(MeV sr)
    !  Xsection=Xsection/1000.

  end function dSigmadOmega_fdE_f_resProd_eN

  !****************************************************************************
  !****f* resProd_lepton/sigma_pipi_res_vac
  ! NAME
  ! function  sigma_pipi_res_vac(targetNuc,q) result(sigma)
  !
  ! PURPOSE
  !
  ! Evaluates the resonance contribution to double-pion production in
  ! gamma N -> R -> N 2Pi scattering.
  ! The return value is sigma in [mb]. Converts target nucleon first to
  ! vacuum nucleon!!!
  ! dsigma/dOmega_electron/dE_electron/dOmega_pion
  !
  ! Assumptions:
  ! * No interferences among resonances.
  !
  ! INPUTS
  ! * type(particle)         :: targetNuc -- Target nucleon
  ! * real, dimension (0:3)  :: q         -- Virtual photon 4-momentum
  ! OUTPUT
  ! * sigma(0:3) Cross sections for gamma+nucleon->nucleon+2Pi production
  ! * sigma(1) -> nucleon piMinus piPlus
  ! * sigma(2) -> nucleon piPlus piNull or nucleon piMinus piNull
  ! * sigma(3) -> nucleon piNull piNull
  ! * sigma(0) Total Xsection into nucleon+2 Pions
  !****************************************************************************
  function sigma_pipi_res_vac(targetNuc,q) result(sigma_pipi)
    use particleDefinition
    use constants, only: mN

    type(particle), intent(in)        :: targetNuc
    real, dimension (0:3), intent(in) :: q
    real,dimension(0:3)               :: sigma_pipi

    type(particle)                    :: targetNuc_VAC     ! Target nucleon, converted to vacuum

    targetNuc_vac=targetNuc

    targetNuc_vac%mass=mN
    targetNuc_vac%position=1000.
    targetNuc_vac%momentum(0)=freeEnergy(targetNuc)

    sigma_pipi=sigma_pipi_res(targetNuc_vac,q)

  end function Sigma_pipi_res_vac

  !****************************************************************************
  !****f* resProd_lepton/sigma_pipi_res
  ! NAME
  ! function  sigma_pipi_res(targetNuc,q) result(sigma)
  !
  ! PURPOSE
  !
  ! Evaluates the resonance contribution to double-pion production in
  ! gamma N -> R -> N 2Pi scattering.
  ! The return value is sigma in [mb].
  !
  ! Assumptions:
  ! * No interferences among resonances.
  !
  ! INPUTS
  ! * type(particle)         :: targetNuc -- Target nucleon
  ! * real, dimension (0:3)  :: q         -- Virtual photon 4-momentum
  ! OUTPUT
  ! * sigma(0:3) Cross sections for gamma+nucleon->nucleon+2Pi production
  ! * sigma(1) -> nucleon piMinus piPlus
  ! * sigma(2) -> nucleon piPlus piNull or nucleon piMinus piNull
  ! * sigma(3) -> nucleon piNull piNull
  ! * sigma(0) Total Xsection into nucleon+2 Pions
  !****************************************************************************
  function sigma_pipi_res(targetNuc,q) result(sigma_pipi)
    use particleDefinition
    use idTable, only: nres, nucleon,pion,P11_1440,rho,sigmaMeson,Delta
    use baryonWidth, only: partialWidthBaryon, FullWidthBaryon
    use clebschGordan, only:clebschSquared
    use particleProperties, only: hadron
    use CALLSTACK, only: TRACEBACK

    type(particle), intent(in)         :: targetNuc
    real, dimension (0:3) , intent(in) :: q
    real,dimension(0:3)                :: sigma_pipi

    real, dimension (nucleon+1:nucleon+nres,0:3) :: sigma ! dsigma/dOmega_electron/dE_electron


    real, dimension (1:3,1:4) :: clebsch
    integer             :: pionCharge_1, pionCharge_2    ! Charges of outgoing pions
    integer :: resID
    !real :: branchingRatio
    integer :: i
    real :: sigma_tot
    real :: res_iz, res_I, nuc_iz
    real :: mass, baremass, fullWidth
    real :: c1,c2,c3,c4
    real, dimension(1:4) :: gamma_Out

    if (sqrt(Dot_Product(targetNuc%position,targetNuc%position)).lt.50) then
       write(*,*) 'In medium not yet implemented in Sigma_pipi_res. STOP!'
       write(*,*) 'position=',targetNuc%position
       call TRACEBACK()
       stop
    end if

    do i=0,3
       sigma(:,i)=0.
    end do


    resIDLoop: do resID=nucleon+1,nucleon+nres
!       pf_res=targetNuc%momentum+q

       sigma_tot=sigma_resProd(targetNuc,resId,q,baremass)
       mass=bareMass

       fullWidth=FullWidthBaryon(resID,mass)
       if (fullWidth.lt.1E-10) then
          sigma(resID,:)=0.
          cycle resIdLoop
       end if

       ! Get gamma N -> pi pi N by summing just over the ...
       ! gamma N -> pion Delta, gamma N -> N rho, gamma N -> N sigma, gamma N -> pion P11_1440
       ! ... channels. All resonances but the P_11(1440) decay fully into N pi.
       ! We correct for the P_11(1440) later.
       gamma_Out(1)=partialwidthBaryon(resID,mass,.false.,pion,Delta)        /fullWidth
       gamma_Out(2)=partialwidthBaryon(resID,mass,.false.,rho,nucleon)       /fullWidth
       gamma_Out(3)=partialwidthBaryon(resID,mass,.false.,sigmaMeson,nucleon)/fullWidth
       ! Here we make the simplifying assumption that the decay ratio of P11_1440 to NPi is constant in mass:
       gamma_Out(4)=partialwidthBaryon(resID,mass,.false.,pion,P11_1440)     /fullWidth*hadron(P11_1440)%decays(1)

       res_I=float(hadron(resId)%isoSpinTimes2)/2.
       res_Iz=targetNuc%charge-1./2.

       clebsch(1,:)=0.
       clebsch(2,:)=0.
       clebsch(3,:)=0.
       loop1: do pionCharge_1=-1,1
          loop2: do pionCharge_2=-1,1
             if (abs(pionCharge_1+pionCharge_2).gt.1) cycle loop2

             ! R ->  pion Delta -> pion pion N
             nuc_iz=res_iZ-real(pionCharge_2+pionCharge_1)
             if (abs(nuc_iz).lt.0.6) then
                c1=       clebschSquared(1.,1.5,res_I,real(pionCharge_1),res_iZ-real(pionCharge_1)) &
                  &    *  clebschSquared(1.,0.5,1.5  ,real(pionCharge_2),nuc_iz)
             else
                c1=0.
             end if

             ! R ->  rho N -> pion pion N
             nuc_iz=res_iZ-real(pionCharge_2+pionCharge_1)
             if (abs(pionCharge_1+pionCharge_2).le.1.and.abs(nuc_iz).lt.0.6) then
                c2=       clebschSquared(1.,0.5,res_I,real(pionCharge_1+pionCharge_2),nuc_iz) &
                     & *  clebschSquared(1.,1.,1.,real(pionCharge_2),real(pionCharge_1))
             else
                c2=0.
             end if

             ! R ->  sigma N -> pion pion N
             if (abs(pionCharge_1+pionCharge_2).le.0) then
                c3= clebschSquared(1.,1.,0.,real(pionCharge_2),real(pionCharge_1))
             else
                c3=0.
             end if

             nuc_iz=res_iZ-real(pionCharge_2+pionCharge_1)
             ! R ->  pion P_11(1440) -> pion pion N
             if (abs(nuc_iZ).lt.0.6.and.abs(res_iZ-real(pionCharge_1)).lt.0.6) then
                c4=       clebschSquared(1.,0.5,res_I,real(pionCharge_1),res_iZ-real(pionCharge_1)) &
                     & *  clebschSquared(1.,0.5,0.5  ,real(pionCharge_2),nuc_iz)
             else
                c4=0.
             end if

             if (((pionCharge_1.eq.1).and.(pionCharge_2.eq.-1)).or.((pionCharge_1.eq.-1).and.(pionCharge_2.eq.1))) then
               clebsch(1,1)=clebsch(1,1)+c1
               clebsch(1,2)=clebsch(1,2)+c2
               clebsch(1,3)=clebsch(1,3)+c3
               clebsch(1,4)=clebsch(1,4)+c4
             else if (((pionCharge_1.eq.1).and.(pionCharge_2.eq.0)).or.((pionCharge_1.eq.0).and.(pionCharge_2.eq.1)) &
                  & .or.((pionCharge_1.eq.-1).and.(pionCharge_2.eq.0)).or.((pionCharge_1.eq.0).and.(pionCharge_2.eq.-1))) then
               clebsch(2,1)=clebsch(2,1)+c1
               clebsch(2,2)=clebsch(2,2)+c2
               clebsch(2,3)=clebsch(2,3)+c3
               clebsch(2,4)=clebsch(2,4)+c4
             else if ((pionCharge_1.eq.0).and.(pionCharge_2.eq.0)) then
               clebsch(3,1)=clebsch(3,1)+c1
               clebsch(3,2)=clebsch(3,2)+c2
               clebsch(3,3)=clebsch(3,3)+c3
               clebsch(3,4)=clebsch(3,4)+c4
              end if
           end do loop2
        end do loop1

        !write(*,'(5E15.4)') sigma_tot, gamma_Out
       do i=1,4
          sigma(resID,1:3)=sigma(resID,1:3)+sigma_tot*gamma_out(i)*clebsch(1:3,i)
       end do
       sigma(resID,0)=sum(sigma(resID,1:3))

    end do resIDLoop

    do i=0,3
       sigma_pipi(i)=sum(sigma(:,i))
    end do
  end function Sigma_pipi_res


  !****************************************************************************
  !****f* resProd_lepton/sigma_barMes_res_vac
  ! NAME
  ! function  sigma_barMes_res_vac(targetNuc,q,IDbar,IDmes) result (sigma)
  !
  ! PURPOSE
  !
  ! Evaluates the resonance contribution of gamma N -> R -> B m^0 scattering
  ! (where X may be a nucleon or Delta, while m^0 is a neutral meson).
  ! The return value is sigma in [mb]. Converts target nucleon first to
  ! vacuum nucleon!!!
  !
  ! Assumptions:
  ! * No interferences among resonances.
  !
  ! INPUTS
  ! * type(particle)         :: targetNuc -- Target nucleon
  ! * real, dimension (0:3)  :: q         -- Virtual photon 4-momentum
  ! * integer                :: IDbar     -- ID of produced baryon N (nucleon or Delta)
  ! * integer                :: IDmes     -- array containing the IDs of produced mesons m
  ! OUTPUT
  ! * sigma -- Cross sections for gamma N -> R -> B m^0 production, for all mesons you asked for
  !****************************************************************************
  function  sigma_barMes_res_vac(targetNuc,q,IDbar,IDmes) result (sigma)
    use particleDefinition
    use constants, only: mN

    type(particle), intent(in)         :: targetNuc     ! Target nucleon
    real, dimension (0:3) , intent(in) :: q             ! Virtual photon 4-momentum
    integer, intent(in)                :: IDbar, IDmes(:)
    real                               :: sigma(lbound(IDMes,dim=1):ubound(IDMes,dim=1))

    type(particle)                     :: targetNuc_VAC     ! Target nucleon, converted to vacuum

    targetNuc_vac=targetNuc

    targetNuc_vac%mass=mN
    targetNuc_vac%position=1000.
    targetNuc_vac%momentum(0)=freeEnergy(targetNuc)

    sigma = sigma_barMes_res(targetNuc_vac,q,IDbar,IDmes)

  end function Sigma_barMes_res_vac



!******************************************************************************
  !****f* resProd_lepton/sigma_barMes_res
  ! NAME
  ! function  sigma_barMes_res(targetNuc,q,IDbar,IDmes) result(sigma_VM)
  !
  ! PURPOSE
  !
  ! Evaluates the resonance contribution of gamma N -> R -> B m^0 scattering
  ! (where B may be a nucleon or Delta, while m^0 is a neutral meson).
  ! The return value is sigma in [mb]. The cross section is calculated separately
  ! for all mesons which are passed in IDmes.
  !
  ! Assumptions:
  ! * No interferences among resonances.
  !
  ! INPUTS
  ! * type(particle)         :: targetNuc -- Target nucleon
  ! * real, dimension (0:3)  :: q         -- Virtual photon 4-momentum
  ! * integer                :: IDbar     -- ID of produced baryon B (nucleon or Delta)
  ! * integer                :: IDmes     -- array containing the IDs of produced mesons m
  ! OUTPUT
  ! * sigma -- Cross sections for gamma N -> R -> B m^0 production, for all mesons you asked for
  !****************************************************************************
  function  sigma_barMes_res(targetNuc,q,IDbar,IDmes) result(sigma_VM)
    use particleDefinition
    use idTable, only: nres,nucleon
    use baryonWidth, only: partialWidthBaryon, FullWidthBaryon
    use particleProperties, only: hadron
    use clebschGordan, only: clebschSquared
    use CALLSTACK, only: TRACEBACK

    type(particle), intent(in)         :: targetNuc     ! Target nucleon
    real, dimension (0:3) , intent(in) :: q             ! Virtual photon 4-momentum
    integer, intent(in)                :: IDbar, IDmes(:)
    real                               :: sigma_VM(lbound(IDmes,dim=1):ubound(IDmes,dim=1))

    real, dimension (nucleon+1:nucleon+nres,lbound(IDmes,dim=1):ubound(IDmes,dim=1))   :: sigma
    integer :: resID,i
    real :: sigma_tot, mass, baremass, fullWidth,res_I,res_Iz,gamma_Out

    if (sqrt(Dot_Product(targetNuc%position,targetNuc%position)).lt.50) then
       write(*,*) 'In medium not yet implemented in Sigma_VM_res. STOP!'
       write(*,*) 'position=',targetNuc%position
       call TRACEBACK()
       stop
    end if

    sigma = 0.

    ! loop over intermediate resonances
    do resID=nucleon+1,nucleon+nres

      sigma_tot=sigma_resProd(targetNuc,resId,q,baremass)
      mass=bareMass

      fullWidth=FullWidthBaryon(resID,mass)
      if (fullWidth.lt.1E-10) cycle

      res_I=float(hadron(resId)%isoSpinTimes2)/2.
      res_Iz=targetNuc%charge-1./2.

      ! loop over meson final states
      do i = lbound(IDmes,dim=1),ubound(IDmes,dim=1)
        ! Branching ratios into V N
        gamma_Out = partialwidthBaryon(resID,mass,.false.,IDmes(i),IDbar)/fullWidth &
                    * clebschSquared(hadron(IDmes(i))%isospinTimes2/2.,hadron(IDbar)%isospinTimes2/2.,res_I,0.,res_Iz)

        sigma(resID,i) = sigma_tot * gamma_out
      end do

    end do

    ! sum over all resonances
    sigma_VM(:) = sum(sigma,dim=1)

  end function Sigma_barMes_res

  !****************************************************************************
  !****f* resProd_lepton/dSdOmega_k_med_res
  ! NAME
  ! function  dSdOmega_k_med_res(targetNuc,q,k,pf) result(sigma_dOmega)
  !
  ! PURPOSE
  !
  ! Evaluates the resonance contribution to pion production in
  ! gamma R->eNPi scattering. The return value
  ! is dsigma/dOmega(pion) in [mub/Sr].
  !
  ! Assumptions:
  ! * No interferences among resonances.
  ! * Isotropic decay of the resonance in its rest-frame.
  ! * In the vacuum.
  !
  ! INPUTS
  ! * type(particle)         :: targetNuc -- Target nucleon
  ! * real, dimension (0:3)  :: q         -- Virtual photon 4-momentum
  ! * real, dimension (0:3)  :: k         -- pion 4-momentum
  ! * real, dimension (0:3)  :: pf        -- Outgoing nucleon 4-momentum
  ! OUTPUT
  ! * logical                :: success -- flag
  ! * real,dimension(-1:1)   :: sigma_dOmega  -- dsigma/dOmega_pion; Index: PionCharge
  !****************************************************************************
  function  dSdOmega_k_med_res(targetNuc,q,k,pf,success) result(sigma_dOmega)
    use constants, only: pi
    use particleDefinition
    use mediumDefinition, only: vacuum
    use idTable, only: nres, nucleon,pion
    use baryonWidthMedium, only: partialWidthBaryonMedium, WidthBaryonMedium
    use clebschGordan, only:clebschSquared
    use particleProperties, only: hadron
    use CALLSTACK, only: TRACEBACK

    type(particle), intent(in)         :: targetNuc
    real, dimension (0:3) , intent(in) :: q
    real, dimension (-1:1,0:3) , intent(in) :: pf
    real, dimension (-1:1,0:3) , intent(in) :: k
    real,dimension(-1:1)  :: sigma_dOmega  ! dsigma/dOmega_electron/dE_electron/dOmega_pion
    real, dimension (-1:1,nucleon+1:nucleon+nres) :: sigma         ! dsigma/dOmega_electron/dE_electron

    real, dimension (0:3) :: pf_res!,q_res
    integer :: pionCharge    ! Charge of outgoing pion
    integer :: resID
    real :: branchingRatio,mass,fullWidth
    integer :: i
    real :: sigma_tot, baremass
    logical :: success

    if (sqrt(Dot_Product(targetNuc%position,targetNuc%position))<=50.) then
       write(*,*) 'In medium not yet implemented in dSdOmega_k_med_res. STOP!'
       write(*,*) 'position=',targetNuc%position
       call TRACEBACK()
       stop
    end if

    do i=-1,1
       sigma(i,:)=0.
    end do

    success=.true.

    resIDLoop: do resID=nucleon+1,nucleon+nres
       pf_res=targetNuc%momentum+q

       sigma_tot=sigma_resProd(targetNuc,resId,q,baremass)
       mass=bareMass

       fullWidth=WidthBaryonMedium(resID,mass,pf_res,vacuum)
       if (fullWidth.lt.1E-10) then
          sigma(:,resID)=0.
          cycle resIdLoop
       end if

       branchingRatio=partialWidthBaryonMedium(resID,mass,.false.,pion,nucleon,pf_res,vacuum)&
            &        /fullWidth

       sigma_tot=sigma_tot*branchingRatio/4./pi
       do pionCharge=targetNuc%charge-1,targetNuc%charge
          sigma(pionCharge,resID)=sigma_tot*dOmegaCM_dOmega(pf_res,k(pionCharge,:),pf(pionCharge,:),success) &
               & *clebschSquared(1.,0.5,float(hadron(resId)%isoSpinTimes2)/2.,&
               & real(pionCharge),real(targetNuc%charge)-1./2.-real(pionCharge))
          if (.not.success) then
             write(*,*)
             write(*,*) "k=",k(pionCharge,:)
             write(*,*) "pf=", pf(pionCharge,:)
             write(*,*) "pf of res=", pf_res(:)
             write(*,*) "resID, pion charge ", resID,pionCharge
             sigma_dOmega=0.
             return
          end if
       end do
    end do resIDLoop

    do i=-1,1
       sigma_dOmega(i)=sum(sigma(i,:))
    end do
  end function dSdOmega_k_med_res




  !****************************************************************************
  !****f* resProd_lepton/dSdO_fdE_fdO_k_med_res_EN
  ! NAME
  ! function dSdO_fdE_fdO_k_med_res_EN(eN,k,pf,processID) result(sigma_dOmega)
  !
  ! PURPOSE
  !
  ! Evaluates the resonance contribution to pion production in
  ! eN->eR->eNPi scattering. The return value
  ! is dsigma/dOmega(electron)/dE(electron)/dOmega(pion) in [mb/GeV/Sr**2].
  !
  ! Assumptions:
  ! * No interferences among resonances.
  ! * Isotropic decay of the resonance in its rest-frame.
  ! * In the vacuum.
  !
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN  -- electron-nucleon scattering event
  ! * real, dimension (0:3)       :: pf  -- Outgoing nucleon 4-momentum
  ! * real, dimension (0:3)       :: k   -- Outgoing pion 4-momentum
  ! * integer, optional           :: processID -- See module leptonicID for
  !   usage
  ! * integer, optional           :: pionNucleonSystem --
  !   If this parameter is set to 1, then we evaluate dOmega_pion in the
  !   calculation frame. If it's 2 then it is evaluated
  !   in the cm frame of the outgoing pion  and nucleon.
  ! OUTPUT
  ! * real,dimension(-1:1):: sigma_dOmega --
  !   dsigma/dOmega_electron/dE_electron/dOmega_pion; Index: PionCharge
  !
  ! NOTES
  ! * Enhances dSdO_fdE_fdO_k_med_res by allowing arbitrary electron
  !   momentum directions
  !****************************************************************************
  function  dSdO_fdE_fdO_k_med_res_EN(eN,k,pf,processID_IN,pionNucleonSystem) result(sigma_dOmega)
    use particleDefinition
    use particleProperties, only: hadron
    use constants, only: pi
    use mediumDefinition, only: vacuum
    use idTable, only: nres, nucleon,pion
    use baryonWidthMedium, only: partialWidthBaryonMedium, WidthBaryonMedium
    use clebschGordan, only: clebschSquared
    use leptonicID, only: EM
    use eN_eventDefinition, only: electronNucleon_event,write_electronNucleon_event
    use CALLSTACK, only: TRACEBACK

    ! Input
    type(electronNucleon_event), intent(in)       :: eN
    real, dimension (0:3)      , intent(in)       :: pf
    real, dimension (0:3)      , intent(in)       :: k
    integer, optional          , intent(in)       :: processID_IN
    integer, optional          , intent(in)       :: pionNucleonSystem
    ! Result:
    real,dimension(-1:1)                          :: sigma_dOmega  ! dsigma/dOmega_electron/dE_electron/dOmega_pion

    real, dimension (-1:1,nucleon+1:nucleon+nres)   :: sigma         ! dsigma/dOmega_electron/dE_electron

    real, dimension (0:3) :: pf_res
    integer               :: pionCharge    ! Charge of outgoing pion
    integer               :: resID
    real                  :: branchingRatio,mass,fullWidth
    integer               :: i, processID
    real                  :: sigma_tot
    real                  :: baremass
    logical               :: piN_inCM_frame


    if (present(pionNucleonSystem)) then
       piN_inCM_frame=(pionNucleonSystem.eq.2)
    end if

    if (sqrt(Dot_Product(eN%nucleon%position,eN%nucleon%position))<=50.) then
       write(*,*) 'In medium not yet implemented in dSdO_fdE_fdO_k_med_res_EN. STOP!'
       write(*,*) 'position=',eN%nucleon%position
       call write_electronNucleon_event(eN,.FALSE.,.FALSE.)
       call TRACEBACK()
       stop
    end if

    sigma(:,:)=0.

    processID=EM
    if (present(processID_IN)) processID=processID_IN

    resIDLoop: do resID=nucleon+1,nucleon+nres
       sigma_tot=dSigmadOmega_fdE_f_resProd_eN(eN,resID,pf_res,baremass,processID)

       mass=bareMass

       fullWidth=WidthBaryonMedium(resID,mass,pf_res,vacuum)
       if (fullWidth.lt.1E-10) then
          sigma(:,resID)=0.
          cycle resIdLoop
       end if

       branchingRatio=partialWidthBaryonMedium(resID,mass,.false.,pion,nucleon,pf_res,vacuum)&
            &        /fullWidth

       if (piN_inCM_frame) then
          sigma_tot=sigma_tot*branchingRatio/4./pi
       else
          sigma_tot=sigma_tot*branchingRatio*dOmegaCM_dOmega(pf_res,k,pf)/4./pi
       end if
       do pionCharge=eN%nucleon%charge-1,eN%nucleon%charge
          sigma(pionCharge,resID)=sigma_tot *   &
               & clebschSquared(1.,0.5,hadron(resId)%isoSpinTimes2/2.,               &
               & real(pionCharge),real(eN%nucleon%charge)-1/2.-real(pionCharge))
       end do
    end do resIDLoop

    do i=-1,1
       sigma_dOmega(i)=sum(sigma(i,:))
    end do
  end function dSdO_fdE_fdO_k_med_res_eN


  !****************************************************************************
  !****f* resProd_lepton/dOmegaCM_dOmega
  ! NAME
  ! real function dOmegaCM_dOmega()
  ! PURPOSE
  ! Evaluates the Jacobian for dOmega_CM(pion)/dOmega_lab(pion)
  !****************************************************************************
  real function dOmegaCM_dOmega(pf_res,k,pf,success)
    use minkowski, only: abs4
    use twoBodyTools, only: pcm
    use constants, only: mN, mPi

    real :: abs_k_vec,abs_pf_vec!,abs_q_vec
    real :: cos_theta_k
    real :: kcm
    real, dimension (0:3),intent(in) :: k,pf,pf_res
    logical, optional :: success

    if (present(success)) success=.true.
    kcm= pCM(abs4(pf_res), mN, mPi)
    if (kcm.lt.1E-10) then
       write(*,*) 'WARNING: Trying to produce resting pion in CM=', k
       write(*,*) 'kcm=', kcm
       dOmegaCM_dOmega=0.
       if (present(success)) success=.false.
       return
    end if

    abs_k_vec=sqrt(Dot_Product(k(1:3),k(1:3)))
    !abs_q_vec=sqrt(Dot_Product(q(1:3),q(1:3)))
    abs_pf_vec=sqrt(Dot_Product(pf(1:3),pf(1:3)))

    ! Make case study for the evaluation of cos(theta)
    if (abs_k_vec.lt.1E-10) then
       write(*,*) 'WARNING: Trying to produce resting pion=', k
       write(*,'(A,4E15.3)') 'k=', k
       !write(*,'(A,4E15.3)') 'energy_li, energy_lf=', energy_li, energy_lf
       dOmegaCM_dOmega=0.
       if (present(success)) success=.false.
       return
    else if (abs_pf_vec.lt.1E-10) then
       dOmegaCM_dOmega=abs4(pf_res)*abs_K_vec**2/kCM/(abs_K_vec*pf(0))

    else
       cos_theta_k=Dot_Product(k(1:3),pf(1:3))/abs_k_vec/abs_pf_vec
       dOmegaCM_dOmega=abs4(pf_res)*abs_K_vec**2/kCM/(abs_K_vec*pf(0)-abs_pf_vec*k(0)*cos_theta_k)
    end if

  end function dOmegaCM_dOmega




  !****************************************************************************
  !****f* resProd_lepton/sigma_resProd
  ! NAME
  ! function sigma_resProd(targetNucleon,resID,q,baremass) result(xSection)
  !
  ! PURPOSE
  !
  ! Evaluates the cross section for gamma N -> R scattering. The return value
  ! is sigma in [mb].
  !
  ! Assumptions:
  ! * No interferences among resonances.
  !
  ! INPUTS
  ! * type(particle)         :: targetNuc -- Target nucleon
  ! * real, dimension (0:3)  :: q         -- Virtual photon 4-momentum
  ! * integer                :: resId     -- ID of resonance
  ! OUTPUT
  ! *  real                  :: baremass  --
  !    bare mass of resonance (mass without scalar potential)
  ! *  real                  :: xSection  -- sigma
  !****************************************************************************
  function sigma_resProd(targetNucleon,resID,q,baremass_res) result(xSection)
    use leptonicID, only: EM
    use spectralFunc, only: specFunc
    use constants, only: GeVSquared_times_mb, pi
    use minkowski, only: SP
    use hadronTensor_ResProd
    use degRad_conversion
    use particleDefinition
    use constants, only: electronChargeSQ

    type(particle)      , intent(in)     :: targetNucleon
    integer             , intent(in)     :: resID             ! resonance ID
    real, dimension(0:3), intent(in)     :: q                 ! Photon momentum
    real                , intent(out)    :: baremass_res          ! Bare mass of resonance

    real                                 :: xSection
    real, dimension(0:3)                 :: pin, pout
    complex,dimension(0:3,0:3)           :: hadronTensor
    complex                              :: matrixElement_Squared
    real                                 :: kinematics
    logical,parameter                    :: debug_this=.false.
    real                                 :: spec

    if (debug_this) then
       write(*,'(A,4E15.4)')      'q=',q
       write(*,'(A,I8)')          'resID=',resID
       write(*,'(A,4E15.4)')      'p=',targetNucleon%momentum
    end if

    baremass_res=0.
    xSection=0.

    pin =targetNucleon%momentum
    pout=pin+q

    ! ********** Cross section ******************

    SPEC=specFunc(resID,targetNucleon%Charge,pout,targetNucleon%position,baremass_res)

    if (hadronTensor_R(pin,pout,resID,targetNucleon%Charge,EM,hadronTensor,baremass_res) ) then
       matrixElement_Squared=1./2.*electronChargeSQ*(-hadronTensor(0,0)+hadronTensor(1,1)+hadronTensor(2,2)+hadronTensor(3,3))
       if (debug_this) then
          write(*,'(A,4E15.3)') 'M=',Real(hadronTensor(0,0)),Real(hadronTensor(1,1)),Real(hadronTensor(2,2)),Real(hadronTensor(3,3))
          write(*,*) real(matrixElement_Squared), &
                     1./2.*(-real(hadronTensor(0,0))+real(hadronTensor(1,1))+real(hadronTensor(2,2))+real(hadronTensor(3,3)))
       end if
    else
       baremass_res=0.
       xsection=0.
       if (debug_this) then
          write(*,'(A,E15.4)') 'Sigma=',xsection
          write(*,*)
       end if
       return
    end if

    if (abs(AIMAG(matrixElement_Squared))>1E-7) write(*,*) 'ERROR: (Matrix element)^2  is imaginary:', matrixElement_Squared

    kinematics=1./(4.*abs(SP(q,pin)))

    Xsection=kinematics *REAL(matrixElement_Squared)*2.*pi * SPEC

    ! ********** UNIT CONVERSION ******************
    ! Convert from 1/GeV**2 to mb
    Xsection=Xsection/GeVSquared_times_mb

    if (Xsection<-1E-9) then
       write(*,*) 'Xsection less than zero'
       write(*,*) targetNucleon%momentum
       write(*,*) targetNucleon%charge
       write(*,*) resID
       write(*,*) q
       stop
    end if

    if (debug_this) then
       write(*,'(A,E15.4)') 'Sigma=',xsection
       write(77,'(3E15.4)') q(0),xsection,baremass_res
       write(*,*)
    end if
  end function sigma_resProd



  subroutine readInput

!!$    use output
!!$
!!$    integer :: ios
!!$    !*************************************************************************
!!$    !****n* resProd_lepton/resonanceProd_electron
!!$    ! NAME
!!$    ! NAMELIST /resonanceProd_electron/
!!$    ! PURPOSE
!!$    ! Namelist for moduleresProd_lepton includes:
!!$    ! * debug
!!$    !*************************************************************************
!!$    NAMELIST /resonanceProd_electron/ debug
!!$
!!$    call Write_ReadingInput('resonanceProd_electron',0)
!!$    rewind(5)
!!$    read(5,nml=resonanceProd_electron,IOSTAT=ios)
!!$    call Write_ReadingInput("resonanceProd_electron",0,ios)
!!$
!!$    write(*,*) 'debug?', debug
!!$
!!$    call Write_ReadingInput('resonanceProd_electron',1)

  end subroutine readInput


end module resProd_lepton
