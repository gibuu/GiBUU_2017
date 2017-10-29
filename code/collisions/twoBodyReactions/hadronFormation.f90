!******************************************************************************
!****m* /hadronFormation
! NAME
! module hadronFormation
! PURPOSE
! This module contains routines handling the formation process of particles
! INPUTS
! namelist "hadronformation"
!
! NOTES
! 1) Old Concept:
!
! Between Creation time and tauForma*gamma, the XS is scaled by
! #leading quarks/#quarks
! (i.e. 0/2,1/2,2/2 for mesons and 0/3,1/3,2/3,3/3 for baryons).
! After tauForma*gamma, the original XS is restored.
! tauForma is in the order of 0.8fm/c.
!
! 2) New Concept: (GetJetSetVec)
!
! We have 4 times:
! - (hard) interaction (INT),
! - production of first quark (P1),
! - production of second quark (P2),
! - formation of hadron (F). (i.e. first meeting point of YoYo-Lines)
! According to input, 2 times are selected. Before time 1, the XS is zero.
! Between the two times, the cross section evolves according to a given power
! of the time variable. If zero is given, the XS is scaled by a given factor.
! After time 2 the XS is set to 1
!
! Unfortunately, GetJetSetVec works not for all particles but produces
! some "funny" numbers for some particles (e.g. tau_F < tau_P1, tau_P2)
! [approx 10% of particles]
! In this case, Concept 1) is used for these particles (!!! ATTENTION !!!)
!
! 3) Enhanced New Concept
!
! as concept 2), but the cross section do not start at 0 for leading
! particles, but at a value ~ 1/Q2, where Q2 has to be defined
! process dependent.
!******************************************************************************

!******************************************************************************
!****n* hadronFormation/hadronformation
! NAME
! NAMELIST hadronformation
! PURPOSE
! Namelist for module "hadronFormation" includes:
! * tauProda
! * tauForma
! * tauFormaFak
! * useJetSetVec
! * powerCS
! * useTimeFrom
! * useTimeTo
! * GuessDiffrTimes
! * useJetSetVec_Q
! * useJetSetVec_R
! * pedestalCS
!******************************************************************************
module hadronFormation

  implicit none
  private

  !****************************************************************************
  !****g* hadronFormation/useJetSetVec
  ! SOURCE
  !
  logical,save, public :: useJetSetVec = .TRUE.
  ! PURPOSE
  ! Flag to select fragmentation time estimates:
  ! * false -> old concept 1)
  ! * true -> new concepts 2) and 3)
  ! NOTES
  ! select false in case of calculations on a nucleon (speed!).
  !****************************************************************************

  !****************************************************************************
  !****g* hadronFormation/useJetSetVec_Q
  ! SOURCE
  !
  logical,save, public :: useJetSetVec_Q = .TRUE.
  ! PURPOSE
  ! if useJetSetVec, then also use Q2 as measure for XS-pedestal, i.e.
  ! select concept 3) instead of concept 2)
  !****************************************************************************

  !****************************************************************************
  !****g* hadronFormation/useJetSetVec_R
  ! SOURCE
  !
  logical,save, public :: useJetSetVec_R = .TRUE.
  ! PURPOSE
  ! if not useJetSetVec_Q, then use rLead as measure for XS-pedestal
  !****************************************************************************

  !****************************************************************************
  !****g* hadronFormation/tauProda
  ! SOURCE
  !
  real,save :: tauProda    = 0.5
  ! PURPOSE
  ! in formation time concept 2) and 3) for "error particles":
  ! production time of non-leading in rest frame of hadron (in fm)
  !****************************************************************************

  !****************************************************************************
  !****g* hadronFormation/tauForma
  ! SOURCE
  !
  real,save :: tauForma    = 0.8
  ! PURPOSE
  ! in formation time concept 1) and in concept 2),3) for "error particles":
  ! formation time in rest frame of hadron (in fm)
  !****************************************************************************

  !****************************************************************************
  !****g* hadronFormation/tauFormaFak
  ! SOURCE
  !
  real,save :: tauFormaFak = 1.0
  ! PURPOSE
  ! in formation time concept 1):
  ! scale factor for constituent quark model,
  ! rescales #(lead quarks)/#quarks
  !****************************************************************************

  !****************************************************************************
  !****g* hadronFormation/powerCS
  ! SOURCE
  !
  real, save    :: powerCS=1.0
  ! PURPOSE
  ! in formation time concept 2):
  ! power of 't' (constant, linear, quadratic)
  !****************************************************************************

  !****************************************************************************
  !****g* hadronFormation/useTimeFrom
  ! SOURCE
  !
  integer, save :: useTimeFrom = 1
  ! PURPOSE
  ! in formation time concept 2):
  ! encode time XS starts to evolve: 1: tP_min, 2: tP_max, 3: tF
  !****************************************************************************

  !****************************************************************************
  !****g* hadronFormation/useTimeTo
  ! SOURCE
  !
  integer, save :: useTimeTo   = 3
  ! PURPOSE
  ! in formation time concept 2):
  ! encode time XS stops to evolve: 1: tP_min, 2: tP_max, 3: tF
  !****************************************************************************

  !****************************************************************************
  !****g* hadronFormation/pedestalCS
  ! SOURCE
  !
  real, save    :: pedestalCS=0.0
  ! PURPOSE
  ! in formation time concept 2):
  ! encode time XS stops to evolve: 1: tP_min, 2: tP_max, 3: tF
  !****************************************************************************

  !****************************************************************************
  !****g* hadronFormation/GuessDiffrTimes
  ! SOURCE
  !
  logical, save :: GuessDiffrTimes = .TRUE.
  ! PURPOSE
  ! if true, then the times for diffractive particles are treated like them
  ! of all other particles, otherwise particles from "diffractive" events
  ! hadronize immediately.
  !****************************************************************************

  logical, save :: initFlag=.true.

  public :: forceInitFormation, formation, SetJSVFormation

contains

  !****************************************************************************
  !****s* hadronFormation/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "hadronFormation"
  !****************************************************************************
  !****************************************************************************
  subroutine readInput
    use output

    NAMELIST /hadronFormation/ tauProda,tauForma,tauFormaFak, &
         useJetSetVec, powerCS, useTimeFrom, useTimeTo, &
         GuessDiffrTimes, useJetSetVec_Q, useJetSetVec_R, pedestalCS

    common /DataGJV/ Arr(3,4,200),EArr(6,200),verb,AtOrigin
    ! Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
    ! EArr(6,nArrMax),     ! errFlag, rank
    ! verb,                ! verbosity
    ! AtOrigin             ! treatment of outmost prod points
    double precision Arr
    integer EArr
    integer verb
    logical AtOrigin
    save /DataGJV/

    integer useTimeCode

    character*(*), dimension (0:3), parameter :: tN = (/'t_int ','tP_min','tP_max','tF    '/)

    integer :: ios
    real :: h_powerCS

    verb = 0
    AtOrigin = .TRUE.

    ! save some defaults:
    h_powerCS = powerCS
    powerCS = -99.9

    call Write_ReadingInput('hadronFormation',0)
    rewind(5)
    read(5,nml=hadronFormation,IOSTAT=ios)
    call Write_ReadingInput('hadronFormation',0,ios)

    write(*,*) 'Use JetSetVec = ',useJetSetVec
    write(*,*) 'Use Q2        = ',useJetSetVec_Q
    write(*,*)

    if (powerCS.eq.-99.9) then
       powerCS = h_powerCS
       if (.not.useJetSetVec) powerCS = 0.0
    end if

    if (useJetSetVec_Q) then
       if (pedestalCS < 0.0001) then
          pedestalCS = 1.0
          write(*,*) '            --- Value of pedestalCS changed !!!'
       end if
       write(*,1001) ' ... power of time : ',powerCS,'   pedestal= r_Lead/Q^2 * ',pedestalCS
    else if (useJetSetVec_R) then
       if (pedestalCS < 0.0001) then
          pedestalCS = 1.0
          write(*,*) '            --- Value of pedestalCS changed !!!'
       end if
       write(*,1001) ' ... power of time : ',powerCS,'   pedestal= r_Lead * ',pedestalCS
    else
       write(*,1001) ' ... power of time : ',powerCS,'   pedestal=',pedestalCS
    end if

    if (useJetSetVec) then
       useTimeCode= 10*useTimeFrom + useTimeTo
       select case (useTimeCode)
       case (0,11,22,33)
          write(*,*) '... evolution     : ',tN(useTimeFrom)
       case (1,2,3,12,13,23)
          write(*,*) '... evolution     : ',tN(useTimeFrom),' ... ',tN(useTimeTo)
       case default
          write(*,*) 'STOP: useTimeFrom and useTimeTo not valid!'
          stop
       end select
       write(*,*) '... verb,AtOrigin = ',verb,AtOrigin
       write(*,*) '... Guess Times for Diffractive particles: ',GuessDiffrTimes
    else
       write(*,*) '  ... evolution     : ',tN(0),' ... ','gamma*tauForma'
    end if



    write(*,*)
    write(*,1000) ' Production Propertime tauProda= ',tauProda,' fm'
    write(*,1000) ' Formation  Propertime tauForma= ',tauForma,' fm'
    write(*,1000) ' ...ReScale Pre-tauForma-XS by:  ',tauFormaFak
    if (tauforma.lt.1e-3) then
       write(*,*) 'tauforma too small!!!!'
       stop
    end if


1000 FORMAT(A,f7.3,A)
1001 FORMAT(A,f7.3,A,f7.3)


    call SetSwitchPYSTRF(useJetSetVec)
!    call SetSwitchLUSTRF(.FALSE.)
    call SetSwitchLUSTRF(useJetSetVec)

    call Write_ReadingInput('hadronFormation',1)

  end subroutine readInput

  !****************************************************************************
  !****s* hadronFormation/forceInitFormation
  ! NAME
  ! subroutine forceInitFormation
  ! PURPOSE
  ! Forces to read the JobCard (if not done already) to initialize this module.
  !****************************************************************************
  subroutine forceInitFormation
    if (initFlag) then
       call readInput
       initFlag=.false.
    end if
  end subroutine forceInitFormation

  !****************************************************************************
  !****s* hadronFormation/formation
  ! NAME
  ! subroutine formation(PertPart,RealPart,time,finalFlag)
  ! PURPOSE
  ! Loop over all particles in PertPart and RealPart and call
  ! for them the subroutine "performFormation"
  !****************************************************************************
  subroutine formation(PertPart,RealPart,time,finalFlag)
    use particleDefinition

    type(particle),intent(inOUT),dimension(:,:) :: RealPart ! real particles
    type(particle),intent(inOUT),dimension(:,:) :: PertPart ! perturbative particles
    real, intent(in) :: time
    logical, intent(in) :: finalFlag

    integer :: numberEnsembles,size_real_dim2,size_pert_dim2,ensemble,index

    ! Initialize at first call
    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    numberEnsembles=size(RealPart,dim=1)

!!$    If (size(PertPart,dim=1).ne.size(RealPart,dim=1)) then
!!$       write(*,*) 'Number of ensembles in real and perturbative particle vectors do not fit'
!!$       Write(*,*) 'Real :' , size(RealPart,dim=1)
!!$       Write(*,*) 'Perturbative :' , size(PertPart,dim=1)
!!$       Write(*,*) 'Critical Error in FormationMain! Stop!'
!!$       Stop
!!$    end if

    size_real_dim2=size(RealPart,dim=2)
    size_pert_dim2=size(PertPart,dim=2)

    ensemble_loop : do ensemble=1,numberEnsembles             !loops over ensemble
       !--- real particle formation ---
       index_loop_real : do index=1,size_real_dim2
          if (RealPart(ensemble,index)%ID <  0) exit index_loop_real
          if (RealPart(ensemble,index)%ID <= 0) cycle index_loop_real
          call performFormation(RealPart(ensemble,index),time,finalFlag)
       end do index_loop_real
       !--- pert. particle formation ---
       index_loop_pert : do index=1,size_pert_dim2
          if (PertPart(ensemble,index)%ID <  0) exit index_loop_pert
          if (PertPart(ensemble,index)%ID <= 0) cycle index_loop_pert
          call performFormation(PertPart(ensemble,index),time,finalFlag)
       end do index_loop_pert
    end do ensemble_loop
  end subroutine formation

  !****************************************************************************
  !****s* hadronFormation/performFormation
  ! NAME
  ! subroutine performFormation(Part,time,finalFlag)
  !
  ! PURPOSE
  ! Set the Flag "in_Formation" and the scaling of the CrossSection
  ! "ScaleCS" according to the "formation time concept": cf. top of the module
  !
  ! NOTES
  ! scaleCS should be set to something like
  !   (#leading quarks/#quarks) * tauformaFak
  !****************************************************************************
  subroutine performFormation(Part,time,finalFlag)
    use particleDefinition

    type(particle),intent(inOUT) :: Part
    real,intent(in) :: time
    logical,intent(in) :: finalFlag
    real :: tForma, timefak, CS0

    if (finalflag) Part%in_Formation=.false. ! Force to hadronize NOW

    if (.not.Part%in_Formation) then
       Part%scaleCS=1.
       return
    end if

    if (time.lt.Part%productionTime) then
       Part%scaleCS = 0
       return
    end if

    ! If there was a problem in GetJetSetVec, we have set the
    ! FormationTime of the Particle to something <0.
    ! Therefore, also the reported production times are not in use.

    if (useJetSetVec.and.Part%formationTime.ge.0.0) then
       tForma = Part%formationTime
    else
       ! This is called for particles with:
       ! useJetSetVec = .FALSE. or formationTime < 0,
       ! but only if in_Formation=.TRUE. !!!
       tForma=Part%productionTime+freeEnergy(Part)/Part%mass*tauForma
    end if

    if (time.lt.tForma) then

       CS0 = pedestalCS
       if (useJetSetVec_Q.or.useJetSetVec_R) call CalcSC0_UseQ()

       if (powerCS==0.0) then
          Part%scaleCS = CS0
       else
          timefak = (time-Part%productionTime)/(tForma-Part%productionTime)
          Part%scaleCS = CS0+(1.0-CS0)*timefak**powerCS
       end if
    else
       Part%scaleCS=1.
       Part%in_Formation=.false.
       if (Part%formationTime.lt.0.0)  Part%formationTime = time
    end if

  contains
    subroutine CalcSC0_UseQ()
      use PIL_nLead

      real    :: r

      if (PIL_nLead_GET(Part%number, r)) then
         CS0 = r * pedestalCS
!         write(*,*) '...GET:',Part%number, r
      else
         CS0 = 0.0
      end if
    end subroutine CalcSC0_UseQ

  end subroutine performFormation

  !****************************************************************************
  !****s* hadronFormation/SetJSVFormation
  ! NAME
  ! subroutine SetJSVFormation(Part,iPart, HardScaleQ2)
  !
  ! PURPOSE
  ! Calculates and sets the "Production/Formation Time" to be stored
  ! with the particle (in case of GetJetSetVec).
  !
  ! INPUTS
  ! * type(particle) :: Part
  ! * integer        :: iPart
  ! * real           :: HardScaleQ2 -- the Q2 value of the production
  !   process connected to this particle
  ! * integer        :: nLead -- number of leading quarks (0...3)
  !
  ! OUTPUT
  ! * type(particle) :: Part
  !****************************************************************************
  subroutine SetJSVFormation(Part,iPart, HardScaleQ2, nLead)

    use particleDefinition
    use IdTable, only: isMeson, isPhoton
    use PIL_FormInfo, only: PIL_FormInfo_PUT

    type(particle),intent(inOUT) :: Part
    integer, intent(in)          :: iPart
    real, intent(in)             :: HardScaleQ2
    integer, intent(in)          :: nLead

    common /DataGJV/ Arr(3,4,200),EArr(6,200),verb,AtOrigin
    ! Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
    ! EArr(6,nArrMax),     ! errFlag, rank
    ! verb,                ! verbosity
    ! AtOrigin             ! treatment of outmost prod points
    double precision Arr
    integer EArr
    integer verb
    logical AtOrigin
    save /DataGJV/

    integer :: iSmaller
    real :: time, gamma

    integer :: h, nL, nL0
    real :: rL

    if (isPhoton(Part%ID)) then
       Part%in_Formation=.false.
       Part%scaleCS=1.
       return
    end if

    nL0 = 3
    if (isMeson(Part%ID)) nL0 = 2
    nL = max(min(nLead,nL0),0)
    rL = real(nL)/nL0

    Part%in_Formation=.true.
    Part%scaleCS=rL*tauFormaFak

    if (.not.useJetSetVec) then

       h = 100*nL
       call PIL_FormInfo_PUT(Part%number, h)

    else

       if (EArr(1,iPart).eq.0) return ! hadron has not been processed

       h = EArr(2,iPart) + 10*EArr(3,iPart) + 100*nL
       call PIL_FormInfo_PUT(Part%number, h)

       ! this is just a bad way to get the "time" when the event happened:
       time = Part%productionTime

       if ((EArr(3,iPart)==6).and. .not.GuessDiffrTimes) then
          Part%in_Formation=.false.
          Part%scaleCS=1.
          return
       end if

       select case (EArr(2,iPart))
       case (0) ! ----- no problems with times:

          call CalcTimeOrder

          select case (useTimeFrom)
          case (0)
             Part%productionTime = time                       ! interaction time
          case (1)
             Part%productionTime = time + Arr(iSmaller,4,iPart)   ! smaller time
          case (2)
             Part%productionTime = time + Arr(3-iSmaller,4,iPart) ! larger time
          case (3)
             Part%productionTime = time + Arr(3,4,iPart)          ! formation time
          end select

          select case (useTimeTo)
          case (0)
             Part%formationTime = time                       ! interaction time
          case (1)
             Part%formationTime = time + Arr(iSmaller,4,iPart)   ! smaller time
          case (2)
             Part%formationTime = time + Arr(3-iSmaller,4,iPart) ! larger time
          case (3)
             Part%formationTime = time + Arr(3,4,iPart)          ! formation time
          end select

       case default ! ----- some severe problems with times:

          gamma = freeEnergy(Part)/Part%mass

          select case (useTimeFrom)
          case (0)
             Part%productionTime = time                  ! interaction time
          case (1)
             if (Part%scaleCS.gt.0.0) then
                Part%productionTime = time
             else
                Part%productionTime = time + gamma*tauProda
             end if
          case (2)
             Part%productionTime = time + gamma*tauProda
          case (3)
             Part%productionTime = time + gamma*tauForma ! formation time
          end select

          select case (useTimeTo)
          case (0)
             Part%formationTime = time                  ! interaction time
          case (1)
             if (Part%scaleCS.gt.0.0) then
                Part%formationTime = time
             else
                Part%formationTime = time + gamma*tauProda
             end if
          case (2)
             Part%formationTime = time + gamma*tauProda
          case (3)
             Part%formationTime = time + gamma*tauForma ! formation time
          end select

       end select
    end if

    ! both times set correctly, set also XS-scaling:

    !    write(*,*) 'in SetJSVFormation:', Part%scaleCS

    if (useJetSetVec_Q) then
       call SetQ_UseQ(HardScaleQ2)
    else if (useJetSetVec_R) then
       call SetQ_UseQ(1.0)
    end if

    call performFormation(Part,time,.FALSE.)

  contains

    !**************************************************************************
    !****s* SetJSVFormation/CalcTimeOrder
    ! NAME
    ! subroutine CalcTimeOrder
    ! PURPOSE
    ! set "iSmaller" to 1 or 2 according 4D points Arr(i,1:4,iPart)
    !
    ! Time ordering is defined to be according: tau
    !
    !**************************************************************************
    subroutine CalcTimeOrder
      if (Arr(1,4,iPart)**2-Arr(1,1,iPart)**2-Arr(1,2,iPart)**2-Arr(1,3,iPart)**2 .lt. &
           & Arr(2,4,iPart)**2-Arr(2,1,iPart)**2-Arr(2,2,iPart)**2-Arr(2,3,iPart)**2) then
         iSmaller = 1
      else
         iSmaller = 2
      end if
    end subroutine CalcTimeOrder

    !**************************************************************************
    !****s* SetJSVFormation/SetQ_UseQ
    ! NAME
    ! subroutine SetQ_UseQ(ScaleQ2)
    ! PURPOSE
    ! Get the Q2 information stored on a event basis and connect it with
    ! the leading particle.
    !
    ! NOTES
    ! * If this routine is called for a nonleading particle : return
    ! * If this routine is called with ScaleQ2==0       : return
    !**************************************************************************
    subroutine SetQ_UseQ(ScaleQ2)
      use PIL_nLead
      use inputGeneral
      use eventtypes

      real, intent(in) :: ScaleQ2

      !logical :: scatterFlag
      !integer :: EType
      !real    :: nu,Q2,eps,Weight

      COMMON/PYINT1/MINT(400),VINT(400)
      integer MINT
      double precision VINT
      SAVE /PYINT1/


      if (Part%scaleCS == 0) return ! not a leading particle
      if (ScaleQ2 == 0) return  ! Q2 value not valid

      call PIL_nLead_PUT(Part%number, Part%scaleCS/ScaleQ2)

!      write(*,*) 'SetQ_UseQ(): Q2 = ',ScaleQ2

    end subroutine SetQ_UseQ


  end subroutine SetJSVFormation

  !****************************************************************************



end module hadronFormation
