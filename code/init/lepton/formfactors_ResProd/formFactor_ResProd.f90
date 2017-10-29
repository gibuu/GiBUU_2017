!******************************************************************************
!****m* /formFactor_ResProd
! NAME
! module formFactor_ResProd
!
! PURPOSE
! Provides the form factors for electro-weak resonance excitation.
!
! INPUTS
! Via FF_ResProd (namelist "input_FF_ResProd" in the Jobcard) one might choose
! whether the form factors are calculated from MAID's helicity amplitudes directly
! or whether the fit of Lalakulich (PRD 74, 014009 (2006)) is used.
!
!******************************************************************************
module formFactor_ResProd
  use leptonicID

  implicit none

  private

  public:: getFormfactor_Res, ff_resProd_setCutoff
!   public :: width_gammaNR

  logical, save :: initFlag=.true.

  !****************************************************************************
  !****g* formFactor_ResProd/debugFlag
  ! SOURCE
  logical, parameter :: debugFlag=.false.
  ! PURPOSE
  ! Switch for debug information
  !****************************************************************************


  !****************************************************************************
  !****g* formFactor_ResProd/FF_ResProd
  ! SOURCE
  !
  integer, save :: FF_ResProd=0
  ! PURPOSE
  ! With FF_ResProd (namelist "input_FF_ResProd" in the Jobcard) one can choose
  ! how the form factors are calculated:
  ! * 0: MAID's helicity amplitudes (Luis' helicity expressions - CM frame)
  ! * 1: fit of Lalakulich (PRD 74, 014009 (2006))
  ! * 2: MAID's helicity amplitudes (Lalakulich's helicity expressions - LAB frame)
  !****************************************************************************



  !****************************************************************************
  !****g* formFactor_ResProd/DeltaAxFF
  ! SOURCE
  !
  integer, save :: DeltaAxFF=1
  ! PURPOSE
  ! choose between different axial form factors for the Delta:
  ! * 1: Adler
  ! * 2: Paschos
  ! * 3: dipol
  !****************************************************************************

  !****************************************************************************
  !****g* formFactor_ResProd/MA
  ! SOURCE
  !
  real, save :: MA=0.95   ! 0.95 tuned to ANL  ! 1.3 tuned to BNL
  ! PURPOSE
  ! delta resonance axial mass parameter.
  !****************************************************************************

  !****************************************************************************
  !****g* formFactor_ResProd/aDelta
  ! SOURCE
  !
  real, save :: aDelta=-0.25
  ! PURPOSE
  ! fit parameter for C_5^A (Adler)
  !****************************************************************************

  !****************************************************************************
  !****g* formFactor_ResProd/bDelta
  ! SOURCE
  !
  real, save :: bDelta=0.04
  ! PURPOSE
  ! fit parameter for C_5^A (Adler)
  !****************************************************************************

  !****************************************************************************
  !****g* formFactor_ResProd/cDelta
  ! SOURCE
  !
  real, save :: cDelta=3.
  ! PURPOSE
  ! fit parameter for C_5^A (Paschos)
  !****************************************************************************

  !****************************************************************************
  !****g* formFactor_ResProd/DeltaCouplrelErr
  ! SOURCE
  !
  real, save :: DeltaCouplrelErr=0.
  ! PURPOSE
  ! error in percent for C_5^A(0) for the Delta
  !****************************************************************************

  !****************************************************************************
  !****g* formFactor_ResProd/HNV_axialFF
  ! SOURCE
  !
  logical, save :: HNV_axialFF=.false.
  !
  ! PURPOSE
  ! With .true. or .false. HNV_axialFF (namelist "input_FF_ResProd" in the
  ! Jobcard)
  ! one can choose which axial form factors to use for Delta-resonance:
  ! * .true. is Hernandez-Nieves-Valverde fit with C5A=0.867,MA=0.985 (PRD 76)
  ! * .false. is as it was used by Lalakulich et al in PRD 74
  !****************************************************************************


  !****************************************************************************
  !****g* formFactor_ResProd/nenner_C5A_Lalakulich
  ! SOURCE
  !
  real, save :: nenner_C5A_Lalakulich=3.0
  !
  ! PURPOSE
  ! Factor wich appear in the Lalakulich parameterization of the axial C_5^A form factor
  ! 3.0 was fitted to BNL and  used in Lalakulich PRD71 and PRD 74
  ! fit of ANL gave 0.5
  !****************************************************************************


  !****************************************************************************
  !****g* formFactor_ResProd/refit_barnu_axialFF
  ! SOURCE
  !
  logical, save :: refit_barnu_axialFF=.false.
  !
  ! PURPOSE
  !
  ! With .true. refit_barnu_axialFF (namelist "input_FF_ResProd" in the Jobcard)
  ! means that the axial form factors are refitted to explain the low value
  ! of antineutrino cross section ( exper data  Bolognese PLB81,393 (1979) )
  !****************************************************************************



  !****************************************************************************
  !****g* formFactor_ResProd/W_cutOff_switch
  ! SOURCE
  !
  logical, save :: W_cutOff_switch=.false.
  ! PURPOSE
  ! Switch to include a W-dependent cut-off function for the vector form
  ! factor of the Delta:
  ! * false = excluded
  ! * true  = included
  !****************************************************************************

  !****************************************************************************
  !****g* formFactor_ResProd/W_cutOff_switchAll
  ! SOURCE
  !
  logical, save :: W_cutOff_switchAll=.false.
  ! PURPOSE
  ! Switch to include a W-dependent cut-off function for the vector and the
  ! axial form factor of all resonances:
  ! * false = excluded
  ! * true  = included
  ! NOTES
  ! we assume the same dependence as for the Delta vector form factor
  !****************************************************************************

  !****************************************************************************
  !****g* formFactor_ResProd/W_cutOff_lambda
  ! SOURCE
  !
  real, save :: W_cutOff_lambda=1.071
  ! PURPOSE
  ! Value for lambda in the W-dependent cut-off function.
  !****************************************************************************

  !****************************************************************************
  !****g* formFactor_ResProd/vector_FF_switch
  ! SOURCE
  !
  logical, save :: vector_FF_switch=.true.
  ! PURPOSE
  ! Switch to turn off the vector form factors:
  ! * false = off
  ! * true  = on
  !****************************************************************************

  !****************************************************************************
  !****g* formFactor_ResProd/axial_FF_switch
  ! SOURCE
  !
  logical, save :: axial_FF_switch=.true.
  ! PURPOSE
  ! Switch to turn off the axial form factors:
  ! * false = off
  ! * true  = on
  !****************************************************************************


contains
  subroutine readInputResProdFormFactors
    use output
    use callstack, only: traceback
 !   use singlePionProductionMAIDlike
 !   use readresparam, only: MADelta

    integer :: ios
    !**************************************************************************
    !****n* formFactor_ResProd/input_FF_ResProd
    ! NAME
    ! NAMELIST /input_FF_ResProd/
    ! PURPOSE
    ! This Namelist for module "formFactor_ResProd" includes:
    ! * FF_ResProd
    ! * aDelta
    ! * bDelta
    ! * cDelta
    ! * DeltaAxFF
    ! * HNV_axialFF
    ! * nenner_C5A_Lalakulich
    ! * refit_barnu_axialFF
    ! * W_cutOff_lambda
    ! * W_cutOff_switch
    ! * vector_FF_switch
    ! * axial_FF_switch
    ! * W_cutOff_switchAll
    ! * DeltaCouplrelErr
    ! * MA
    !**************************************************************************
    NAMELIST /input_FF_ResProd/ FF_ResProd,aDelta,bDelta,cDelta,&
         & DeltaAxFF,HNV_axialFF, &
         & nenner_C5A_Lalakulich, refit_barnu_axialFF, &
         & W_cutOff_lambda,W_cutOff_switch,vector_FF_switch,axial_FF_switch,&
         & W_cutOff_switchAll,DeltaCouplrelErr,MA


    call Write_ReadingInput('input_FF_ResProd',0)
    rewind(5)
    read(5,nml=input_FF_ResProd,IOSTAT=ios)
    call Write_ReadingInput("input_FF_ResProd",0,ios)

    select case (FF_ResProd)
    case (0)
       write(*,*) 'form factors for resonance production: MAID with Luis helicity amplitudes'
    case (1)
       write(*,*) 'form factors for resonance production: Lalakulich Fit'
    case (2)
       write(*,*) 'form factors for resonance production: MAID with Lalakulich helicity amplitudes'
    case default
       write(*,*) 'form factors for resonance production: something strange', FF_ResProd
       call traceback()
    end select

    write(*,*) 'delta axial mass: ', MA
    write(*,*) 'assumed error in delta C5A(0): ', DeltaCouplrelErr
    select case (deltaAxFF)
    case (1)
       write(*,'(a,2F12.4)') ' using Adlers form for delta C_5^A with fit parameters a, b: ', aDelta, bDelta
    case (2)
       write(*,'(a,F12.4)') ' using Paschos form for delta C_5^A with fit parameters c: ', cDelta
    case (3)
       write(*,'(a,F12.4)') ' using dipole form for delta C_5^A'
    case default
       write(*,*) 'wrong input for deltaAxFF:',deltaAxFF
       call traceback()
    end select

    if (W_cutOff_switch) then
       write(*,*) 'W-dependence in vector form factor of Delta is INCLUDED!'
       write(*,'(a,F12.4)') ' -> Value for the cut-off factor lambda:', W_cutOff_lambda
    end if
    if (W_cutOff_switchAll) then
       write(*,*) 'W-dependence in all form factors of all resonances is INCLUDED!'
       write(*,'(a,F12.4)') ' -> Value for the cut-off factor lambda:', W_cutOff_lambda
    end if
    if (.not.W_cutOff_switch.and..not.W_cutOff_switchAll) then
       write(*,*) 'W-dependence in form factors is EXCLUDED!'
    end if

    if (.not.vector_FF_switch) write(*,*) 'WARNING: vector FF are switched off!!!'
    if (.not.axial_FF_switch)  write(*,*) 'WARNING: axial FF are switched off!!!'

    if (HNV_axialFF) then
       if (FF_ResProd.eq.1) then
          write(*,*) 'use HNV axial form factor'
       else
          write(*,*) 'NOTE: HNV_axialFF is set to true but has no effect for FF_ResProd.ne.1',FF_ResProd
       end if
    end if


    if (refit_barnu_axialFF .and. HNV_axialFF) then
       call traceback('Only one of the parameters HNV_axialFF  _OR_  refit_barnu_axialFF can be true, not both. STOP')
    end if

    call Write_ReadingInput('input_FF_ResProd',1)
  end subroutine readInputResProdFormFactors


  !****************************************************************************
  !****f* formFactor_ResProd/getFormfactor_Res
  ! NAME
  ! function getFormfactor_Res (Qs, bare_mass, resID, targetCharge, process, FF_set)
  !
  ! PURPOSE
  ! This function serves as an interface to the subroutines which calculate
  ! the form factors.
  !
  ! INPUTS
  ! * real      :: Qs            -- momentum transfer
  ! * real      :: bare_mass     -- bare mass of the resonance (no potentials included!)
  ! * integer   :: resID         -- ID of resonance (see IDtable)
  ! * integer   :: targetCharge  -- charge of target (nucleon)
  ! * integer   :: process       -- EM, CC or NC (or anti if negative)
  ! * logical   :: FF_set        -- controll flag whether form factors
  !                                 have been set successfully
  !
  ! OUTPUT
  ! real, dimension(1:8)  :: getFormfactor_Res --
  ! first 4 entries contain the vector form factors,
  ! last 4 entries contain the axial form factors
  !****************************************************************************
  function getFormfactor_Res (Qs, bare_mass, resID, targetCharge, process, FF_set) result (formfactorfield)
    use IDTable, only: Delta
    use particleProperties, only: hadron
    use distributions, only: markusPostFormfactor

    real, intent(in) :: Qs, bare_mass
    integer, intent(in) :: resID, process, targetCharge
    logical, intent(out) :: FF_set             !true if form factors set, false if form factors could not be set
    real, dimension(1:8) :: formfactorfield    !convention: 4 x FF_V, 4 x FF_A

    !*** Read input:
    if (initFlag) then
       call readInputResProdFormFactors
       initFlag=.false.
    end if

    formfactorfield=0.

    FF_set=.true.

    if (FF_ResProd.eq.0) then
       call vectorformfactor_MAID_CM(formfactorfield(1:4),Qs,resID,targetCharge,process,FF_set)
       if (process.ne.EM) call axialformfactor(formfactorfield(5:8),Qs,resID,targetCharge,process,FF_set)

    else if (FF_ResProd.eq.1) then
       call vectorformfactor_Lalakulich(formfactorfield(1:4),Qs,resID,targetCharge,process,FF_set)
       if (process.ne.EM) call axialformfactor_Lalakulich(formfactorfield(5:8),Qs,resID,targetCharge,process,FF_set)

    else if (FF_ResProd.eq.2) then
       call vectorformfactor_MAIDLala(formfactorfield(1:4),Qs,resID,targetCharge,process,FF_set)
       if (process.ne.EM) call axialformfactor(formfactorfield(5:8),Qs,resID,targetCharge,process,FF_set)

    end if

    ! W-dependence for the vector form factors of the Delta
    if (W_cutoff_switch .and. (resID==Delta) .and. bare_mass>0) then
      formfactorfield(1:4) = formfactorfield(1:4) &
                           * markusPostFormfactor(bare_mass,hadron(resID)%mass,hadron(resID)%minmass,W_cutoff_lambda)
    end if
    ! W-dependence for all form factors of all resonances
    if (W_cutoff_switchAll .and. bare_mass>0) then
      formfactorfield(1:8) = formfactorfield(1:8) &
                           * markusPostFormfactor(bare_mass,hadron(resID)%mass,hadron(resID)%minmass,W_cutoff_lambda)
    end if

    if (.not.vector_FF_switch) formfactorfield(1:4)=0.
    if (.not.axial_FF_switch)  formfactorfield(5:8)=0.


  end function getFormfactor_Res





  !****************************************************************************
  !****s* formFactor_ResProd/vectorformfactor_MAID_CM
  ! NAME
  ! subroutine vectorformfactor_MAID_CM(vecformfactor,Qs,resID,targetCharge,process,FF_set)
  !
  ! PURPOSE
  ! returns the vector form factors calculated using the MAID helicity
  ! amplitudes, helicity amplitudes from Luis are used (CM frame!!)
  !
  ! INPUTS
  ! * real      :: Qs            -- momentum transfer
  ! * integer   :: resID         -- ID of resonance (see IDtable)
  ! * integer   :: targetCharge  -- charge of target (nucleon)
  ! * integer   :: process       -- EM, CC or NC (or anti if negative)
  ! * logical   :: FF_set        -- controll flag whether form factors
  !                                 have been set successfully
  !
  ! OUTPUT
  ! * real, dimension(1:4):: vecformfactor  --   contains the vector FF
  !****************************************************************************
  subroutine vectorformfactor_MAID_CM(vecformfactor,Qs,resID,targetCharge,process,FF_set)
    use IDtable
    use particleProperties, only: hadron
    use constants, only: sinsthweinbg
    use output, only: write_initstatus

    real,  intent(in) :: Qs
    integer, intent(in) :: resID, process, targetCharge
    real, dimension(1:4), intent(out) :: vecformfactor
    logical, intent(inout) :: FF_set             !true if form factors set, false if form factors could not be set

    real,dimension(0:1) :: C3,C4,C5
    integer, parameter :: QS_max_index=400 ! Number of grid points in Q^2 for tabulation
    real, parameter    :: QS_max=2.          ! Maximal Q^2 for tabulating the form factors
    real, save         :: QS_delta
    real,dimension(0:1,nucleon+1:nucleon+nres,0:QS_max_index),save :: C3_table,C4_table,C5_table

    logical, dimension(nucleon+1:nucleon+nres),save :: FF_SET_table
    integer :: QSquared_index
    logical ,save :: firstTime=.true.
    logical, parameter :: useTable=.true.

    vecformfactor=0.

    if (resid.eq.nucleon) then
       ff_set=.false.
       return
    end if

    if (firstTime.and.useTable) then
       firstTime=.false.
       QS_delta=QS_max/float(QS_max_index)
       call write_initstatus('Tabulating form factors for resonance production using CM helicity amplitudes',0)
       call tabulate(FF_SET_table,C3_table,C4_table,C5_table)
       call write_initstatus('Tabulating form factors for resonance production using CM helicity amplitudes',1)
    end if

    if (QS>0. .and. QS<QS_max .and. useTable) then
       ! Use look-up table
       QSquared_index=Nint(QS/QS_delta)
       C3=C3_table(0:1,resID,QSquared_index)
       C4=C4_table(0:1,resID,QSquared_index)
       C5=C5_table(0:1,resID,QSquared_index)
       FF_SET=FF_SET_table(resID)
    else
       ! Calculate explicitly
       call MAID_CM(C3,C4,C5,Qs,resID,FF_set)
    end if

    if (.not.FF_SET) then
       return
    end if


    select case (process)

    case (EM,antiEM)
       vecformfactor(1)=C3(targetCharge)
       vecformfactor(2)=C4(targetCharge)
       vecformfactor(3)=C5(targetCharge)
       FF_set=.true.

    case (CC,antiCC)
       select case (hadron(resID)%isoSpinTimes2)
       case (3)
          !see appendix of Luis notes for the minus
          vecformfactor(1)=-C3(targetCharge)
          vecformfactor(2)=-C4(targetCharge)
          vecformfactor(3)=-C5(targetCharge)
       case (1)
          vecformfactor(1)=C3(proton)-C3(neutron)
          vecformfactor(2)=C4(proton)-C4(neutron)
          vecformfactor(3)=C5(proton)-C5(neutron)
       end select
       FF_set=.true.

    case (NC,antiNC)
       select case (hadron(resID)%isoSpinTimes2)
       case (3)
          !see appendix of Luis notes for the minus
          vecformfactor(1)=-(1.-2.*sinsthweinbg)*C3(targetCharge)
          vecformfactor(2)=-(1.-2.*sinsthweinbg)*C4(targetCharge)
          vecformfactor(3)=-(1.-2.*sinsthweinbg)*C5(targetCharge)
       case (1)
          if (targetcharge.eq.proton) then
             vecformfactor(1)=(0.5-2.*sinsthweinbg)*C3(proton)-0.5*C3(neutron)
             vecformfactor(2)=(0.5-2.*sinsthweinbg)*C4(proton)-0.5*C4(neutron)
             vecformfactor(3)=(0.5-2.*sinsthweinbg)*C5(proton)-0.5*C5(neutron)
          else
             vecformfactor(1)=(0.5-2.*sinsthweinbg)*C3(neutron)-0.5*C3(proton)
             vecformfactor(2)=(0.5-2.*sinsthweinbg)*C4(neutron)-0.5*C4(proton)
             vecformfactor(3)=(0.5-2.*sinsthweinbg)*C5(neutron)-0.5*C5(proton)
          end if
       end select
       FF_set=.true.

    case default
       if (debugflag) write(*,*) 'strange processID -> STOP ',process
       stop

    end select


  contains


    subroutine tabulate(FF_SET_table,C3_table,C4_table,C5_table)
      real,dimension(0:1,nucleon+1:nucleon+nres,0:QS_max_index),intent(out) :: C3_table,C4_table,C5_table
      logical, dimension(nucleon+1:nucleon+nres),intent(out) :: FF_SET_table
      real, dimension(0:1) :: C3,c4,c5
      real :: QSquared
      integer :: ID
      logical :: flag
      integer :: i
      resonanceLoop: do Id=nucleon+1,nucleon+nres
         QSquaredLoop: do i=0,QS_max_index
            QSquared=float(i)*QS_Delta
            call MAID_CM(C3,C4,C5,QSquared,ID,flag)
            C3_table(0:1,ID,i)=C3
            C4_table(0:1,ID,i)=C4
            C5_table(0:1,ID,i)=C5
            FF_SET_table(ID)=flag
         end do QSquaredLoop
      end do resonanceLoop

    end subroutine tabulate

    subroutine MAID_CM(C3,C4,C5,Qs,resID,FF_set)
      use constants, only: pi,alphaQED,mN
      use particleProperties, only: hadron
      use helicityAmplitudes

      real               , intent(in)    :: Qs
      integer            , intent(in)    :: resID
      logical            , intent(inout) :: FF_set             !true if form factors set, false if form factors could not be set
      real,dimension(0:1), intent(out)   :: C3,C4,C5
      real,dimension(0:1) :: A12,A32,S12

      real,dimension(0:1) :: F1,F2
      real :: ms, q2!,mu

      FF_set=.true.

      q2=-Qs

      ms=hadron(resID)%mass
      !      mu=ms+mn

      F1=0.
      F2=0.
      C3=0.
      C4=0.
      C5=0.

      if (((mn + ms)**2 - q2).lt.0..or.((mn - ms)**2 - q2).lt.0.) then
         vecFormfactor=0.
         FF_set=.false.
         return
      end if

      call get_helicityAmplitudes(proton,resID,Qs,A12(proton),A32(proton),S12(proton))
      call get_helicityAmplitudes(neutron,resID,Qs,A12(neutron),A32(neutron),S12(neutron))

      select case (resID)


      case (Delta,P13_1720,F15_1680,F35_1905,F37_1950)

         C3=(((3*A12 + Sqrt(3.)*A32)*ms)/(Sqrt(3.)*mn*Sqrt(pi)*((mn + ms)**2 - q2)*    &
              &  Sqrt((alphaQED*(-(mn - ms)**2 + q2))/(mn**5*(mn - ms)*(mn + ms)))))

         C4=-(Sqrt(0.6666666666666666)*(Sqrt(6.)*A32*((mn - ms)**2 - q2)*(mn**2 - mn*ms +       &
              & ms**2 - q2)*Sqrt((alphaQED*(-(mn + ms)**2 + q2))/(mn*(mn - ms)*(mn + ms))) -   &
              & 3*mn*ms*(Sqrt(2.)*A12*((mn - ms)**2 - q2)*Sqrt((alphaQED*(-(mn + ms)**2 + q2))/ &
              & (mn*(mn - ms)*(mn + ms))) + 2*mn*ms*Sqrt((alphaQED*(-(mn - ms)**2 + q2))/      &
              & (mn**5*(mn - ms)*(mn + ms)))*(mn**2 - ms**2 + q2)*S12)))/ (Sqrt(pi)*           &
              & ((mn - ms)**2 - q2)**2*((mn + ms)**2 - q2)*Sqrt((alphaQED*(-(mn - ms)**2 + q2))&
              & /(mn**5*(mn - ms)*(mn + ms)))*Sqrt((alphaQED*(-(mn + ms)**2 + q2))/(mn*        &
              & (mn - ms)*(mn + ms))))

         C5=((Sqrt(0.6666666666666666)*ms**2*(-3*Sqrt(2.)*A12*((mn - ms)**2 - q2)*            &
              & Sqrt((alphaQED*(-(mn + ms)**2 + q2))/(mn*(mn - ms)*(mn + ms))) + Sqrt(6.)*A32* &
              & ((mn - ms)**2 - q2)*Sqrt((alphaQED*(-(mn + ms)**2 + q2))/(mn*(mn - ms)*(mn +   &
              & ms))) + 6*mn**2*Sqrt((alphaQED*(-(mn - ms)**2 + q2))/(mn**5*(mn - ms)*(mn + ms &
              & )))*(-mn**2 + ms**2 + q2)*S12))/(Sqrt(pi)*((mn - ms)**2 - q2)**2*((mn + ms)**2 &
              & - q2)*Sqrt((alphaQED*(-(mn - ms)**2 + q2))/(mn**5*(mn - ms)*(mn + ms)))*       &
              & Sqrt((alphaQED*(-(mn + ms)**2 + q2))/(mn*(mn - ms)*(mn + ms)))))


      case (P11_1440,P31_1910)

         F1=(2*(A12*Sqrt((alphaQED*ms)/(-mn**2 + ms**2))*((mn - ms)**2 - q2)*Sqrt((mn + ms)**2 &
              & - q2) + 2*Sqrt(2.)*mn**2.5*ms**1.5*(mn + ms)*Sqrt((alphaQED*(-(mn - ms)**2 +   &
              & q2))/(mn**5*(mn - ms)*(mn + ms)))*S12))/ (Sqrt((alphaQED*ms)/(-2*mn**2 + 2*    &
              & ms**2))*Sqrt(pi)*((mn - ms)**2 - q2)*((mn + ms)**2 - q2)**1.5*Sqrt((alphaQED*  &
              & (-(mn - ms)**2 + q2))/(mn**5*(mn - ms)*(mn + ms))))

         F2=(A12*(mn + ms)*Sqrt((alphaQED*ms)/(-mn**2 + ms**2))*((mn - ms)**2 - q2)*Sqrt((mn   &
              & + ms)**2 - q2) + 2*Sqrt(2.)*mn**2.5*ms**1.5*q2*Sqrt((alphaQED*(-(mn - ms)**2 +  &
              & q2))/(mn**5*(mn - ms)*(mn + ms)))*S12)/ (mn*Sqrt((alphaQED*ms)/(-2*mn**2 + 2*  &
              & ms**2))*Sqrt(pi)*Sqrt((mn + ms)**2 - q2)* Sqrt((alphaQED*(-(mn - ms)**2 + q2)) &
              & /(mn**5*(mn - ms)*(mn + ms)))*(mn**4 + (ms**2 - q2)**2 - 2*mn**2*(ms**2 + q2)))

         !         C3=F1*mu**2/(2*mn)**2
         !         C4=F2*mu/(2*mn)
         C3=F1
         C4=F2


      case (D13_1520,D33_1700,D15_1675)

         C3=((3*A12 - Sqrt(3.)*A32)*mn*ms)/(Sqrt(3.)*Sqrt(pi)*((mn - ms)**2 - q2)*Sqrt((       &
              & alphaQED*(-(mn + ms)**2 + q2))/(mn*(mn - ms)*(mn + ms))))

         C4=(Sqrt(0.6666666666666666)*(Sqrt(6.)*A32*mn**2*(mn**2 + mn*ms + ms**2 - q2)*((mn +   &
              & ms)**2 - q2)*Sqrt((alphaQED*(-(mn - ms)**2 + q2))/(mn**5*(mn - ms)*(mn + ms))) &
              & + 3*ms*(-(Sqrt(2.)*A12*mn**3*((mn + ms)**2 - q2)*Sqrt((alphaQED*(-(mn - ms)**2 &
              & + q2))/(mn**5*(mn - ms)*(mn + ms)))) + 2*ms*(-mn**2 + ms**2 - q2)*Sqrt((       &
              & alphaQED*(-(mn + ms)**2 + q2))/(mn*(mn - ms)*(mn + ms)))*S12)))/(Sqrt(pi)*((mn &
              & - ms)**2 - q2)*((mn + ms)**2 - q2)**2*Sqrt((alphaQED*(-(mn - ms)**2 + q2))/(   &
              & mn**5*(mn - ms)*(mn + ms)))* Sqrt((alphaQED*(-(mn + ms)**2 + q2))/(mn*(mn -    &
              & ms)*(mn + ms))))

         C5=-((Sqrt(0.6666666666666666)*ms**2*(3*Sqrt(2.)*A12*mn**2*((mn + ms)**2 - q2)*       &
              & Sqrt((alphaQED*(-(mn - ms)**2 + q2))/(mn**5*(mn - ms)*(mn + ms))) + Sqrt(6.)*  &
              & A32*mn**2*((mn + ms)**2 - q2)*Sqrt((alphaQED*(-(mn - ms)**2 + q2))/(mn**5*(mn  &
              & - ms)*(mn + ms))) + 6*(-mn**2 + ms**2 + q2)*Sqrt((alphaQED*(-(mn + ms)**2 +    &
              & q2))/(mn*(mn - ms)*(mn + ms)))*S12))/(Sqrt(pi)*((mn - ms)**2 - q2)*((mn + ms   &
              & )**2 - q2)**2* Sqrt((alphaQED*(-(mn - ms)**2 + q2))/(mn**5*(mn - ms)*(mn +     &
              & ms)))* Sqrt((alphaQED*(-(mn + ms)**2 + q2))/(mn*(mn - ms)*(mn + ms)))))


      case (S11_1535,S31_1620,S11_1650)

         F1=(2*mn**2*(Sqrt(2.)*A12*Sqrt(-((alphaQED*((mn - ms)**2 - q2)*((mn + ms)**2 - q2)**2)/ &
              & (mn*(mn - ms)*ms**2*(mn + ms)))) + 4*(mn - ms)*Sqrt((alphaQED*(-(mn + ms)**2 +   &
              & q2))/(mn*(mn - ms)*(mn + ms)))*S12))/(Sqrt(pi)*((mn - ms)**2 - q2)*Sqrt(-((      &
              & alphaQED*((mn - ms)**2 - q2)*((mn + ms)**2 - q2)**2)/(mn*(mn - ms)*ms**2*(mn +   &
              & ms))))*Sqrt((alphaQED*(-(mn + ms)**2 + q2))/(mn*(mn - ms)*(mn + ms))))

         F2=-((Sqrt(2.)*mn*(A12*(mn - ms)*Sqrt(-((alphaQED*((mn - ms)**2 - q2)*((mn + ms)**2 -   &
              & q2)**2)/ (mn*(mn - ms)*ms**2*(mn + ms)))) +  2*Sqrt(2.)*q2*Sqrt((alphaQED*(-(mn  &
              & + ms)**2 + q2))/(mn*(mn - ms)*(mn + ms)))*S12))/ (Sqrt(pi)*((mn - ms)**2 - q2)*  &
              & Sqrt(-((alphaQED*((mn - ms)**2 - q2)*((mn + ms)**2 - q2)**2)/(mn*(mn - ms)*ms**2 &
              & *(mn + ms))))*Sqrt((alphaQED*(-(mn + ms)**2 + q2))/(mn*(mn - ms)*(mn + ms)))))

         ! C3=F1*mu**2/(2*mn)**2
         ! C4=F2*mu/(2*mn)
         C3=F1
         C4=F2


      case default

         if (debugflag) write(*,*) 'formfactors of resonance ', resID, ' not yet included'
         FF_set=.false.

      end select


    end subroutine MAID_CM


  end subroutine vectorformfactor_MAID_CM




  !****************************************************************************
  !****s* formFactor_ResProd/vectorformfactor_MAIDLala
  ! NAME
  ! subroutine vectorformfactor_MAIDLala(vecformfactor,Qs,resID,targetCharge,process,FF_set)
  !
  ! PURPOSE
  ! returns the vector form factors calculated using the MAID helicity
  ! amplitudes, helicity amplitudes from Lalakulich are used
  !
  ! INPUTS
  ! * real      :: Qs            -- momentum transfer
  ! * integer   :: resID         -- ID of resonance (see IDtable)
  ! * integer   :: targetCharge  -- charge of target (nucleon)
  ! * integer   :: process       -- EM, CC or NC (or anti if negative)
  ! * logical   :: FF_set        -- controll flag whether form factors
  !                                 have been set successfully
  !
  ! OUTPUT
  ! * real, dimension(1:4):: vecformfactor --    contains the vector FF
  !****************************************************************************
  subroutine vectorformfactor_MAIDLala(vecformfactor,Qs,resID,targetCharge,process,FF_set)
    use minkowski, only: SP
    use IDtable
    use particleproperties, only: hadron
    use constants, only: pi, alphaQED
    use helicityAmplitudes, only: get_helicityAmplitudes

    real,  intent(in) :: Qs
    integer, intent(in) :: resID, process, targetCharge
    real, dimension(1:4), intent(out) :: vecformfactor
    logical, intent(inout) :: FF_set             !true if form factors set, false if form factors could not be set

    real,dimension(0:1) :: C3,C4,C5
    integer, parameter :: QS_max_index=400 ! Number of grid points in Q^2 for tabulation
    real, parameter    :: QS_max=2.          ! Maximal Q^2 for tabulating the form factors
    real, save         :: QS_delta
    real,dimension(0:1,nucleon+1:nucleon+nres,0:QS_max_index),save :: C3_table,C4_table,C5_table

    logical, dimension(nucleon+1:nucleon+nres),save :: FF_SET_table
    integer :: QSquared_index
    logical ,save :: firstTime=.true.
    logical, parameter :: useTable=.true.

    vecformfactor=0.

    if (resid.eq.nucleon) then
       ff_set=.false.
       return
    end if

    if (firstTime.and.useTable) then
       firstTime=.false.
       QS_delta=QS_max/float(QS_max_index)
       write(*,*)
       write(*,*) 'Tabulating Form Factors for resonance production using Lalakulich helicity amplitudes ...'
       call tabulate(FF_SET_table,C3_table,C4_table,C5_table)
       write(*,*) '... done'
       write(*,*)
    end if

    if (QS.lt.QS_max.and.useTable) then
       ! Use look-up table
       QSquared_index=Nint(QS/QS_delta)
       C3=C3_table(0:1,resID,QSquared_index)
       C4=C4_table(0:1,resID,QSquared_index)
       C5=C5_table(0:1,resID,QSquared_index)
       FF_SET=FF_SET_table(resID)
    else
       ! Calculate explicitly
       call MAIDLala(C3,C4,C5,Qs,resID,FF_set)
    end if

    if (.not.FF_SET) then
       return
    end if


    select case (process)

    case (EM,antiEM)
       vecformfactor(1)=C3(targetCharge)
       vecformfactor(2)=C4(targetCharge)
       vecformfactor(3)=C5(targetCharge)
       FF_set=.true.

    case (CC,antiCC)
       select case (hadron(resID)%isoSpinTimes2)
       case (3)
          vecformfactor(1)=-C3(targetCharge)
          vecformfactor(2)=-C4(targetCharge)
          vecformfactor(3)=-C5(targetCharge)
       case (1)
          vecformfactor(1)=C3(proton)-C3(neutron)
          vecformfactor(2)=C4(proton)-C4(neutron)
          vecformfactor(3)=C5(proton)-C5(neutron)
       end select
       FF_set=.true.

    case (NC,antiNC)
       write(*,*) 'NC not yet implemented'
       FF_set=.false.

    case default
       if (debugflag) write(*,*) 'strange processID -> STOP ',process
       stop

    end select


  contains


    subroutine tabulate(FF_SET_table,C3_table,C4_table,C5_table)
      use IdTable, only: nucleon, nres
      real,dimension(0:1,nucleon+1:nucleon+nres,0:QS_max_index),intent(out) :: C3_table,C4_table,C5_table
      logical, dimension(nucleon+1:nucleon+nres),intent(out) :: FF_SET_table
      real, dimension(0:1) :: C3,c4,c5
      real :: QSquared
      integer :: ID
      logical :: flag
      integer :: i
      resonanceLoop: do Id=nucleon+1,nucleon+nres
         QSquaredLoop: do i=0,QS_max_index
            QSquared=float(i)*QS_Delta
            call MAIDLala(C3,C4,C5,QSquared,ID,flag)
            C3_table(0:1,ID,i)=C3
            C4_table(0:1,ID,i)=C4
            C5_table(0:1,ID,i)=C5
            FF_SET_table(ID)=flag
         end do QSquaredLoop
      end do resonanceLoop

    end subroutine tabulate

    subroutine MAIDLala(C3,C4,C5,Qs,resID,FF_set)
      use particleProperties, only: hadron
      use constants, only: mN

      real               , intent(in)    :: Qs
      integer            , intent(in)    :: resID
      logical            , intent(inout) :: FF_set             !true if form factors set, false if form factors could not be set
      real,dimension(0:1), intent(out)   :: C3,C4,C5
      real,dimension(0:1) :: A12,A32,S12

      real :: MR, mu,q0,qz,ppr0,qtimesp,qtimesppr
      real :: sqrtN, sqrt2


      real  :: W           ! invariant mass of resonance

      real, dimension(0:3) :: q,pin,pf

      FF_set=.true.

      W=hadron(resID)%mass
      call getKine_lab(W,QS,q,pin,pf)

      MR=hadron(resID)%mass
      mu=MR+mN

      q0=q(0)
      qz=q(3)
      ppr0=pf(0)

      qtimesp=SP(q,pin)
      qtimesppr=SP(q,pf)

      if ((pi * alphaQED * 2. * mN * (ppr0+MR))/(mN*(W**2-mN**2)).lt.0) then
         vecFormfactor=0.
         FF_set=.false.
         return
      end if
      sqrtN=sqrt((pi * alphaQED * 2. * mN * (ppr0+MR))/(mN*(W**2-mN**2)))

      sqrt2=sqrt(2.)

      C3=0.
      C4=0.
      C5=0.

      call get_helicityAmplitudes(proton,resID,Qs,A12(proton),A32(proton),S12(proton))
      call get_helicityAmplitudes(neutron,resID,Qs,A12(neutron),A32(neutron),S12(neutron))

      select case (resID)

      case (Delta,P13_1720,F15_1680,F35_1905,F37_1950)

         !compared to Lalakulich, there is an overall minus sign less here
         !her eq. 4.4 holds only for isospin 3/2 resonances, we include this minus sign
         !explicitely later in select case (processID)

         C3=((Sqrt(3.)*A12 + A32)*MR)/(2.*qz*SqrtN)

         C4=-(mN*(3*Sqrt2*A12*MR*(mN*(mN + MR)*(mN + q0) - MR*qtimesp)*qz - Sqrt(6.)*A32*    &
              & (mN*((mN - MR)*MR + 2*mN*ppr0)*(mN + q0) + MR**2*qtimesp)*qz -              &
              & 6*mN*MR*(MR + ppr0)*qtimesp*S12))/(2.*Sqrt(6.)*qz**2*SqrtN*(mN*(mN + q0)*   &
              & qtimesppr - qtimesp*W**2))

         C5=-(mN*(3*Sqrt2*mN*MR*(MR + ppr0)*qtimesppr*S12 + 3*A12*MR*qz*(MR*qtimesppr -      &
              & (mN + MR)*W**2) + Sqrt(3.)*A32*qz*(MR**2*qtimesppr + ((mN - MR)*MR +        &
              & 2*mN*ppr0)*W**2)))/(2.*Sqrt(3.)*qz**2*SqrtN*(mN*(mN + q0)*qtimesppr - qtimesp*W**2))



      case (P11_1440,P31_1910)

         C3=(mu**2*(MR + ppr0)*(A12*qz + sqrt2*(mN + MR)*S12))/(sqrt2*sqrtN*((mN + MR)**2 + Qs)*qz**2)

         C4=(mu*(MR + ppr0)*(A12*(mN + MR)*qz - sqrt2*Qs*S12))/(sqrt2*sqrtN*((mN + MR)**2 + Qs)*qz**2)

         ! Different definition of electromagnetic currents
         C3=C3/mu**2*(2*mn)**2
         C4=C4/mu*(2*mn)



      case (D13_1520,D33_1700,D15_1675)

         !compared to Lalakulich, there is an overall minus sign less here
         !her eq. 4.4 holds only for isospin 3/2 resonances, we include this minus sign
         !explicitely later in select case (processID)

         C3=-((-3*A12 + Sqrt(3.)*A32)*MR*(MR + ppr0))/(2.*Sqrt(3.)*sqrtN*qz**2)

         C4=(mN*(3*sqrt2*A12*MR*(MR + ppr0)*(mN*(mN - MR)*(mN + q0) + MR*qtimesp) -             &
              &  Sqrt(6.)*A32*(MR*(MR + ppr0)*(mN*(mN - MR)*(mN + q0) + MR*qtimesp) +               &
              &  2*mN**2*(mN + q0)*qz**2) - 6*mN*MR*qtimesp*qz*S12))/                              &
              &  (2.*Sqrt(6.)*sqrtN*((-(mN*(mN + 2*q0)) + Qs)*qtimesp + mN*(mN + q0)*qtimesppr)*qz**2)

         C5=(mN*(-3*sqrt2*A12*MR*(MR + ppr0)*((mN - MR)*(mN**2 + 2*mN*q0 - Qs) + MR*qtimesppr) + &
              &  Sqrt(6.)*A32*(MR*(MR + ppr0)*((mN - MR)*(mN**2 + 2*mN*q0 - Qs) + MR*qtimesppr) +    &
              &  2*mN*(mN**2 + 2*mN*q0 - Qs)*qz**2) + 6*mN*MR*qtimesppr*qz*S12))/                   &
              &  (2.*Sqrt(6.)*sqrtN*((-(mN*(mN + 2*q0)) + Qs)*qtimesp + mN*(mN + q0)*qtimesppr)*qz**2)



      case (S11_1535,S31_1620,S11_1650)

         C3=(mu**2*(A12*qz + sqrt2*(mN - MR)*S12))/(sqrt2*sqrtN*((mN - MR)**2 + Qs)*qz)

         C4=(mu*(A12*(-mN + MR)*qz + sqrt2*Qs*S12))/(sqrt2*sqrtN*((mN - MR)**2 + Qs)*qz)

         ! Different definition of electromagnetic currents
         C3=C3/mu**2*(2*mn)**2
         C4=C4/mu*(2*mn)

      case default

         if (debugflag) write(*,*) 'formfactors of resonance ', resID, ' not yet included'
         FF_set=.false.

      end select


    end subroutine MAIDLala


  end subroutine vectorformfactor_MAIDLala




  !****************************************************************************
  !****s* formFactor_ResProd/vectorformfactor_Lalakulich
  ! NAME
  ! subroutine vectorformfactor_Lalakulich(vecformfactor,Qs,resID,targetCharge,process,FF_set)
  !
  ! PURPOSE
  ! returns the vector form factors according to a fit by Lalakulich et al.
  ! PRD 74, 014009 (2006).
  !
  ! INPUTS
  ! * real      :: Qs            -- momentum transfer
  ! * integer   :: resID         -- ID of resonance (see IDtable)
  ! * integer   :: targetCharge  -- charge of target (nucleon)
  ! * integer   :: process       -- EM, CC or NC (or anti if negative)
  ! * logical   :: FF_set        -- controll flag whether form factors
  !                                 have been set successfully
  !
  ! OUTPUT
  ! * real, dimension(1:4):: vecformfactor --    contains the vector FF
  !****************************************************************************
  subroutine vectorformfactor_Lalakulich(vecformfactor,Qs,resID,targetCharge,process,FF_set)
    use IDtable
    use particleproperties, only: hadron
    use constants, only: sinsthweinbg, mN

    real, intent(in) :: Qs
    integer, intent(in) :: resID, process, targetCharge
    real, dimension(1:4), intent(out) :: vecformfactor
    logical, intent(inout) :: FF_set             !true if form factors set, false if form factors could not be set

    real,dimension(0:1) :: C3,C4,C5

    real, parameter :: MV=0.84
    real :: MV2, DV
    real :: mu,mr

    MR=hadron(resID)%mass
    mu=MR+mN

    vecformfactor=0.
    C3=0.
    C4=0.
    C5=0.

    MV2=MV**2

    DV=(1.+Qs/MV2)**2

    select case (resID)

    case (Delta) ! cf. eq. (4.20) of Lalakulich

       C3(proton)=2.13/DV/(1.+Qs/4./MV2)
       C3(neutron)=C3(proton)

       C4(proton)=-1.51/DV/(1.+Qs/4./MV2)
       C4(neutron)=C4(proton)

       C5(proton)=0.48/DV/(1.+Qs/0.776/MV2)
       C5(neutron)=C5(proton)


    case (P11_1440) ! cf. eq. (4.24) of Lalakulich

       C3(proton)=2.3/DV/(1.+Qs/4.3/MV2)
       C3(neutron)=-C3(proton)

       C4(proton)=-0.76/DV*(1.-2.8*log(1.+Qs))
       C4(neutron)=-C4(proton)

       ! Different definition of electromagnetic currents
       C3=C3/mu**2*(2*mn)**2
       C4=C4/mu*(2*mn)


    case (D13_1520) ! cf. eq. (4.8) of Lalakulich

       C3(proton)=2.95/DV/(1.+Qs/8.9/MV2)
       C3(neutron)=-1.13/DV/(1.+Qs/8.9/MV2)

       C4(proton)=-1.05/DV/(1.+Qs/8.9/MV2)
       C4(neutron)=0.46/DV/(1.+Qs/8.9/MV2)

       C5(proton)=-0.48/DV
       C5(neutron)=-0.17/DV

    case (S11_1535) ! cf. eq. (4.33) of Lalakulich

       C3(proton)=2.0/DV/(1.+Qs/1.2/MV2)*(1.+7.2*log(1.+Qs))
       C3(neutron)=-C3(proton)

       C4(proton)=0.84/DV*(1.+0.11*log(1.+Qs))
       C4(neutron)=-C4(proton)

       ! Different definition of electromagnetic currents
       C3=C3/mu**2*(2*mn)**2
       C4=C4/mu*(2*mn)

    case default

       if (debugflag) write(*,*) 'formfactors of resonance ', resID, ' not yet included'
       FF_set=.false.


    end select

    if (.not.FF_SET) then
       return
    end if


    select case (process)

    case (EM,antiEM)
       vecformfactor(1)=C3(targetCharge)
       vecformfactor(2)=C4(targetCharge)
       vecformfactor(3)=C5(targetCharge)
       FF_set=.true.

    case (CC,antiCC)
       select case (hadron(resID)%isoSpinTimes2)
       case (3)
          vecformfactor(1)=C3(targetCharge)
          vecformfactor(2)=C4(targetCharge)
          vecformfactor(3)=C5(targetCharge)
       case (1)
          vecformfactor(1)=C3(neutron)-C3(proton)  !Lalakulich uses a different definition here than the rest of the world
          vecformfactor(2)=C4(neutron)-C4(proton)
          vecformfactor(3)=C5(neutron)-C5(proton)
          !vecformfactor(1)=C3(proton)-C3(neutron)   ! conventional here and then positive axial FF
          !vecformfactor(2)=C4(proton)-C4(neutron)
          !vecformfactor(3)=C5(proton)-C5(neutron)
       end select
       FF_set=.true.

    case (NC,antiNC)
       select case (hadron(resID)%isoSpinTimes2)
       case (3)
          vecformfactor(1)=(1.-2.*sinsthweinbg)*C3(targetCharge)
          vecformfactor(2)=(1.-2.*sinsthweinbg)*C4(targetCharge)
          vecformfactor(3)=(1.-2.*sinsthweinbg)*C5(targetCharge)
       case (1)
          if (targetcharge.eq.proton) then
             vecformfactor(1)=(0.5-2.*sinsthweinbg)*C3(proton)-0.5*C3(neutron)
             vecformfactor(2)=(0.5-2.*sinsthweinbg)*C4(proton)-0.5*C4(neutron)
             vecformfactor(3)=(0.5-2.*sinsthweinbg)*C5(proton)-0.5*C5(neutron)
          else
             vecformfactor(1)=(0.5-2.*sinsthweinbg)*C3(neutron)-0.5*C3(proton)
             vecformfactor(2)=(0.5-2.*sinsthweinbg)*C4(neutron)-0.5*C4(proton)
             vecformfactor(3)=(0.5-2.*sinsthweinbg)*C5(neutron)-0.5*C5(proton)
          end if
       end select
       FF_set=.true.

    case default
       if (debugflag) write(*,*) 'strange processID -> STOP ',process
       stop

    end select


  end subroutine vectorformfactor_Lalakulich





  !****************************************************************************
  !****s* formFactor_ResProd/axialformfactor
  ! NAME
  ! subroutine axialformfactor(axialformfactor,Qs,resID,targetCharge,process,FF_set)
  !
  ! PURPOSE
  ! returns the axial form factors, see Tina's notes and Mathematica
  ! notebook 'axialff_coupling.nb' for the calculation of the axial couplings
  !
  ! INPUTS
  ! * real      :: Qs            -- momentum transfer
  ! * integer   :: resID         -- ID of resonance (see IDtable)
  ! * integer   :: targetCharge  -- charge of target (nucleon)
  ! * integer   :: process       -- EM, CC or NC (or anti if negative)
  ! * logical   :: FF_set        -- controll flag whether form factors
  !                                 have been set successfully
  !
  ! OUTPUT
  ! * real, dimension(1:4):: axialff --    contains the axial FF
  !****************************************************************************
  subroutine axialformfactor(axialff,Qs,resID,targetCharge,process,FF_set)
    use idtable
    use particleproperties, only: hadron
    use constants, only: mN, mPi

    real, intent(in) :: Qs
    integer, intent(in) :: resID, process, targetCharge
    real, dimension(1:4), intent(out) :: axialff
    logical, intent(inout) :: FF_set             !true if form factors set, false if form factors could not be set
    real :: C3A,C4A,C5A,C6A
    real :: C5A0,C3A0

    real :: DA

    real, parameter :: MAresSq=1.


    axialff=0.
    FF_set=.true.

    C3A=0.
    C4A=0.
    C5A=0.
    C6A=0.


    DA=1./(1.+Qs/MAresSq)**2


    if (resid.eq.nucleon) then
       ff_set=.false.
       return
    end if


    select case (resID)

    case (P11_1440,P31_1910)

       if (resID.eq.P11_1440) C3A0=-0.520198
       if (resID.eq.P31_1910) C3A0=0.0758086

       C3A=C3A0*DA

       C4A=C3A*(hadron(resID)%mass+mN)*mN/(Qs+mPi**2)


    case (S11_1535,S31_1620,S11_1650)

       if (resID.eq.S11_1535) C3A0=-0.226645
       if (resID.eq.S31_1620) C3A0=0.0519172
       if (resID.eq.S11_1650) C3A0=-0.250395

       C3A=C3A0*DA

       C4A=C3A*(hadron(resID)%mass-mN)*mN/(Qs+mPi**2)



    case (Delta,P13_1720,F15_1680,F35_1905,F37_1950)

       if (resID.eq.Delta)    C5A0=1.17466
       if (resID.eq.P13_1720) C5A0=-0.293934
       if (resID.eq.F15_1680) C5A0=-0.431952
       if (resID.eq.F35_1905) C5A0=0.148606
       if (resID.eq.F37_1950) C5A0=0.235433

       C3A=0.

       C4A=0.

       !modified dipol form
       if (resID.eq.Delta) then
          if (deltaAxFF.eq.1) DA=1./(1.+Qs/MA**2)**2*(1.+aDelta*Qs/(bDelta+Qs))
          if (deltaAxFF.eq.2) DA=1./(1.+Qs/MA**2)**2/(1.+Qs/cDelta/MA**2)
          if (deltaAxFF.eq.3) DA=1./(1.+Qs/MA**2)**2
       end if

       if (resID.eq.Delta.and.abs(DeltaCouplrelErr).gt.0) C5A0=C5A0+DeltaCouplrelErr/100.*C5A0


       C5A=C5A0*DA

       C6A=mN**2*C5A/(mPi**2+Qs)


       if (resID.eq.Delta) C4A=-C5A/4.   !Adler model



    case (D13_1520,D33_1700,D15_1675)

       if (resID.eq.D13_1520) C5A0=-2.15421
       if (resID.eq.D33_1700) C5A0=0.839608
       if (resID.eq.D15_1675) C5A0=-1.38016

       C3A=0.

       C4A=0.

       C5A=C5A0*DA

       C6A=mN**2*C5A/(mPi**2+Qs)


    case default

       if (debugflag) write(*,*) 'formfactors of resonance ', resID, ' not yet included'
       FF_set=.false.

    end select


    if (.not.FF_SET) then
       return
    end if


    select case (process)

    case (EM,antiEM)
       write(*,*) 'for EM process you should not be in the axialformfactor routine -> STOP'
       stop

    case (CC,antiCC)
       axialff(1)=C3A
       axialff(2)=C4A
       axialff(3)=C5A
       axialff(4)=C6A
       FF_set=.true.

    case (NC,antiNC)
       select case (hadron(resID)%isoSpinTimes2)
       case (3)
          axialff(1)=C3A
          axialff(2)=C4A
          axialff(3)=C5A
       case (1)
          if (targetcharge.eq.proton) then
             axialff(1)=0.5*C3A
             axialff(2)=0.5*C4A
             axialff(3)=0.5*C5A
          else
             axialff(1)=-0.5*C3A
             axialff(2)=-0.5*C4A
             axialff(3)=-0.5*C5A
          end if
       end select
       axialff(4)=0.
       FF_set=.true.

    case default
       if (debugflag) write(*,*) 'strange processID -> STOP ',process
       stop

    end select


  end subroutine axialformfactor







  !****************************************************************************
  !****s* formFactor_ResProd/axialformfactor_Lalakulich
  ! NAME
  ! subroutine axialformfactor_Lalakulich(axialformfactor,Qs,resID,targetCharge,process,FF_set)
  !
  ! PURPOSE
  ! returns the axial form factors according to a fit by Lalakulich et al.
  ! PRD 74, 014009 (2006).
  !
  ! INPUTS
  ! * real      :: Qs            -- momentum transfer
  ! * integer   :: resID         -- ID of resonance (see IDtable)
  ! * integer   :: targetCharge  -- charge of target (nucleon)
  ! * integer   :: process       -- EM, CC or NC (or anti if negative)
  ! * logical   :: FF_set        -- controll flag whether form factors
  !                                 have been set successfully
  !
  ! OUTPUT
  ! * real, dimension(1:4):: axialformfactor --    contains the axial FF
  !****************************************************************************
  subroutine axialformfactor_Lalakulich(axialff,Qs,resID,targetCharge,process,FF_set)
    use IDtable
    use particleproperties, only: hadron
    use constants, only: mN, mPi

    real, intent(in) :: Qs
    integer, intent(in) :: resID, process, targetCharge
    real, dimension(1:4), intent(out) :: axialff
    logical, intent(inout) :: FF_set             !true if form factors set, false if form factors could not be set

    real :: C3A,C4A,C5A,C6A

    real,parameter :: MAp=1.05

    real :: MA2, DA

    axialff=0.
    C3A=0.
    C4A=0.
    C5A=0.
    C6A=0.

    MA2=MAp**2

    DA=(1.+Qs/MA2)**2


    select case (resID)

    case (Delta)

       if (HNV_axialFF) then
          C5A=0.867/(1.+Qs/0.985**2)**2/(1.+Qs/3./0.985**2)
       else
          C5A=1.2/DA/(1.+Qs/nenner_C5A_Lalakulich/MA2)
       end if

       C3A=0.

       if (refit_barnu_axialFF) then
          C4A=0
       else
          C4A=-C5A/4.
       end if


       C6A=mN**2*C5A/(mPi**2+Qs)


    case (P11_1440) ! cf. eq. (4.30) of Lalakulich

       C3A=-0.51/DA/(1.+Qs/3./MA2)

       C4A=C3A*(hadron(P11_1440)%mass+mN)*mN/(Qs+mPi**2)


    case (D13_1520) ! cf. eq. (4.9), (4.10) of Lalakulich

       C3A=0.

       C4A=0.

       C5A=-2.1/DA/(1.+Qs/3./MA2)

       C6A=mN**2*C5A/(mPi**2+Qs)


    case (S11_1535) ! cf. eq. (4.34) of Lalakulich

       C3A=-0.21/DA/(1.+Qs/3./MA2)

       C4A=C3A*(hadron(S11_1535)%mass-mN)*mN/(Qs+mPi**2)


    case default

       if (debugflag) write(*,*) 'formfactors of resonance ', resID, ' not yet included'
       FF_set=.false.


    end select

    if (.not.FF_SET) then
       return
    end if


    select case (process)

    case (EM,antiEM)
       write(*,*) 'for EM process you should not be in the axialformfactor routine -> STOP'
       stop

    case (CC,antiCC)
       axialff(1)=C3A
       axialff(2)=C4A
       axialff(3)=C5A
       axialff(4)=C6A

       if (hadron(resID)%Spin.lt.1) then
          axialff(1)=-C3A
          axialff(2)=-C4A    !we use a different sign in the axial spin 1/2 current than Lalakulich
       end if

       FF_set=.true.

    case (NC,antiNC)
       select case (hadron(resID)%isoSpinTimes2)
       case (3)
          axialff(1)=C3A
          axialff(2)=C4A
          axialff(3)=C5A
       case (1)
          if (targetcharge.eq.proton) then
             axialff(1)=-0.5*C3A  ! "-" sign because in PRD74 the vector FF was defines as neutron-proton
             axialff(2)=-0.5*C4A  ! correspondingly the axial form factors  were set negative
             axialff(3)=-0.5*C5A
          else
             axialff(1)=0.5*C3A   ! "-" sign because in PRD74 the vector FF was defines as neutron-proton
             axialff(2)=0.5*C4A   ! correspondingly the axial form factors  were set negative
             axialff(3)=0.5*C5A
          end if
       end select
       axialff(4)=0.
       FF_set=.true.


    case default
       if (debugflag) write(*,*) 'strange processID -> STOP ',process
       stop

    end select

  end subroutine axialformfactor_Lalakulich








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  subroutine getKine_lab(W,QS,q,pi,pf)
    ! Establishes final and initial nucleon momentum and the momentum exchange q
    ! Assumptions: Initial nucleon at rest, vec(q) in direction of z-axis
    ! Input : W=sqrt(s)
    ! Input : QS=-SP(q,q)
    !
    use constants, only: mN
    real, dimension(0:3),intent(out) :: q,pi,pf
    real :: W, QS
    pi(1:3)=0.
    pi(0)=mN
    q(0)=(W**2-mN**2+QS)/(2.*mN)
    q(1:2)=0.
    q(3)=sqrt(QS+q(0)**2)
    pf=pi+q
  end subroutine getKine_lab


  ! Drechsel et al., nucl-th

!   real function width_gammaNR(id,W,k_w_out,k_r_out)
!     use particleProperties, only: hadron
!     use constants, only: mN
!
!     integer, intent(in) :: ID
!     real   , intent(in) :: W
!     real :: k_W, k_R,f_gammaN
!     real,optional :: k_w_out,k_r_out
!     real, parameter :: X_R=0.500
!     integer :: n
!
!     k_W=(              W**2-mN**2)/(2*W)
!     k_R=(hadron(id)%mass**2-mN**2)/(2*hadron(id)%mass)
!
!     n=2
!
!     f_gammaN=( (k_W/k_R)**n)*(X_R**2+k_R**2)/(X_R**2+k_W**2)
!     width_gammaNR=f_gammaN**2
!
!     if(present(k_w_out)) k_w_out=k_w
!     if(present(k_r_out)) k_r_out=k_r
!
!   end function width_gammaNR


  !****************************************************************************
  !****s* formFactor_ResProd/ff_resProd_setCutoff
  ! NAME
  ! subroutine ff_resProd_setCutoff(x)
  ! INPUTS
  ! real, intent(in) :: x  -- new value for lambda in the W-dependence of the form factor
  ! PURPOSE
  ! Sets a new cut-off for the W-dependent form factor, in particular it sets W_cutoff_lambda=x.
  !****************************************************************************
  subroutine ff_resProd_setCutoff(x)
    real, intent(in) :: x
    ! Read input
    if (initFlag) then
       call readInputResProdFormFactors
       initFlag=.false.
    end if
    ! Reset lambda
    W_cutoff_lambda=x
    write(*,'(A,/,A,G13.5)') 'WARNING: New value for lambda in module formFactor_ResProd!' &
         & ,' -> W_cutoff_lambda=', W_cutoff_lambda
  end subroutine ff_resProd_setCutoff


end module formFactor_ResProd
