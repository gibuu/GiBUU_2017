!******************************************************************************
!****m* /CollTools
! PURPOSE
! This module contains routines needed by the modules DoColl_...
!******************************************************************************
module CollTools

  IMPLICIT NONE
  private

  public :: SetSomeDefaults_PY, SetSomeDefaults_PY_Hermes, &
       SetSomeDefaults_FR, &
       SetVectorFromPYJETS, &
       ConvertInPartPythia, ConvertInPartFritiof, &
       ReduceCharge, CorrectChargeVector, &
       ConvertToAnti, &
       PythiaCKIN, PythiaMSTP, &
       CheckUndecayedString


  !****************************************************************************
  !****s* CollTools/ConvertToAnti
  ! NAME
  ! subroutine ConvertToAnti(part)
  !
  ! subroutine ConvertToAnti(part, partNEW)
  !
  ! subroutine ConvertToAnti(parts)
  ! PURPOSE
  ! Convert particle(s) to anti-particle(s).
  ! INPUTS
  ! * type(particle) :: part
  ! * type(particle), dimension(:) :: parts
  ! OUTPUT
  ! This depends on input/output combinations:
  !
  ! Case 1:
  ! * type(particle) :: part
  !   (The input particle is replaced by the new values)
  !
  ! Case 2:
  ! * type(particle) :: partNEW
  !   (The input values are not modified; all information is transfered to
  !   a new copy and this copy is modified)
  !
  ! Case 3:
  ! * type(particle), dimension(:) :: parts
  !   (For every particle in the vector, Case 1) is called)
  !
  ! NOTES
  ! This function is overloaded.
  !****************************************************************************
  Interface ConvertToAnti
     Module Procedure ConvertToAnti_1, ConvertToAnti_11, ConvertToAnti_2
  end Interface


contains


  !****************************************************************************
  !****s* CollTools/SetSomeDefaults_PY
  ! NAME
  ! subroutine SetSomeDefaults_PY
  ! PURPOSE
  ! Set some default values for (new) PYTHIA/JETSET.
  !****************************************************************************
  subroutine SetSomeDefaults_PY
    use ID_translation, only: SetBruteForceMasses_PY

    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    integer MSTU,MSTJ
    double precision PARU,PARJ
    SAVE /PYDAT1/

    COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
    integer KCHG
    double precision PMAS,PARF,VCKM
    SAVE /PYDAT2/

    COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
    integer MDCY,MDME,KFDP
    double precision BRAT
    SAVE /PYDAT3/

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
    integer MSEL,MSELPD,MSUB,KFIN
    double precision CKIN
    SAVE /PYSUBS/

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

    logical GetSwitchPythiaHermes

    logical, parameter :: lBruteForceMasses = .TRUE.

    MSEL = 1  ! 1: default, 2: + elastic & diffractive

    MSTJ(21) = 0              ! particle decay on/off

    if (MSTP(182).lt.400) then
       MSTU(12) = 0              ! writing of Pythia title page, v6.2
    else
       MSTU(12) = 12345          ! writing of Pythia title page, v6.4
    end if
    MSTU(16) = 2              ! needed by GetLeading
    MSTU(21) = 1              ! don't stop after errors
    MSTU(22) = 0              ! #errors to be written
    MSTU(25) = 0              ! don't write warnings
    MSTU(26) = 0              ! #warnings to be written

    PARP(2) = 0.1d0           ! lowest c.m. energy
    CKIN(3) = 0d0             ! min pThat (allow low pT events)
    PARP(111) = 0d0           ! otherwise problems at low sqrts
    MSTP(111) = 1             ! master switch fragmentation/decay
    MSTP(122) = 0             ! switch init and max XS print-out

    PARP(104)=0.000d0         ! Min.energy for XS definition

                              ! forbid the decay of:
    MDCY(102,1)=0             ! pi0
!    MDCY(113,1)=0             ! K0
!    MDCY(112,1)=0             ! K_S0

    MDME(556:560,1) = 0       ! switch off rho_0-decay-channels
                              ! (works only with v6.206)
                              ! necessary, but should not be

!    MSTP(142) = 0             ! use reweighted events

!    PARP( 91)=0.44            ! width intrinsic kT

    CKIN( 77)=1.50            ! Wmin (see also CKIN(3),CKIN(5))
!    CKIN(  5)=0.75            ! min pThat (singular processes)

! COMMENTS ON CKIN(5):
! the minimal W in gamma*N is calculated as max(CKIN(77),2*CKIN(3),2*CKIN(5)).
! unfortunately changing CKIN(5) may change the overall yields.
! We just want to allow DIS events for W<2GeV, but W>CKIN(77)
! without changing the B+B and B+M behaviour.


    if (GetSwitchPythiaHermes()) call SetSomeDefaults_PY_Hermes

    if (lBruteForceMasses) call SetBruteForceMasses_PY

    call SetSomeDefaults_PY_Job


  end subroutine SetSomeDefaults_PY


  !****************************************************************************
  !****s* CollTools/SetSomeDefaults_PY_Hermes
  ! NAME
  ! subroutine SetSomeDefaults_PY_Hermes
  ! PURPOSE
  ! Set some default values for (new) PYTHIA/JETSET as tuned by the HERMES
  ! collab. Meant to be called after a call to SetSomeDefaults_PY.
  ! NOTES
  ! To use this, this routine has to be called immediately before
  ! a call to InitPythia (when used for "Virtual Photons").
  !
  ! Parameters-setting for own calculation are not set here, but in
  ! "InitPythia". Some "useless" paramtere changings are excluded or
  ! marked with "(???)". No support for MSTP(199)=1 ("use of radgen")!!!
  !****************************************************************************
  subroutine SetSomeDefaults_PY_Hermes

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    integer MSTU,MSTJ
    double precision PARU,PARJ
    SAVE /PYDAT1/

    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
    integer MSEL,MSELPD,MSUB,KFIN
    double precision CKIN
    SAVE /PYSUBS/

!... PYTHIA-PARAMETER:

    MSTP(  1)=2               ! #generations  (???)
    MSTP( 17)=6               ! r for longitudinal photons
                              ! =6 not provided by standard PYTHIA.
                              ! make sure to link the right routines
    MSTP( 20)=0               ! suppression of resolved XS
    MSTP( 41)=1               ! handling Resonance Decays (???)

    MSTP( 61)=0               ! master: ISR (QCD/QED)
    MSTP( 71)=0               ! master: FSR (QCD/QED)
    MSTP( 81)=0               ! master: multiple interactions

    MSTP( 92)=4               ! energy splitting in remnant
    MSTP(101)=1               ! structure of diffr. system

!    MSTP(199)=1               ! use RADGEN

    PARP( 18)=0.17            ! scale GVMD <-> VMD
    PARP( 62)=0.50            ! cut-off (Q evolution) space-like showers
    PARP( 65)=0.50            ! cut-off (min. energy) space-like showers

    PARP( 91)=0.44            ! width intrinsic kT
    PARP( 93)=2.00            ! cut-off intrinsic kT

    PARP( 99)=0.44            ! width intrinsic kT (photon)
    PARP(100)=2.00            ! cut-off intrinsic kT (photon)

    PARP(102)=0.50            ! start of diffr.states mass spectrum
    PARP(103)=0.50            ! (cf. diffr.states mass spectrum)

    PARP(104)=0.30            ! Min.energy for XS definition

    PARP(161)=2.69            ! VMD-Coupling: rho
    PARP(162)=24.6            ! VMD-Coupling: omega
    PARP(163)=18.8            ! VMD-Coupling: phi
!    PARP(164)=                ! VMD-Coupling: J/Psi

    PARP(165)=0.33            ! factor for transversal XS

!... JETSET-PARAMETER:

    PARJ(  1)=0.03            ! Prob(qq)/Prob(q)
    PARJ(  2)=0.12            ! Prob(s)/Prob(u)
    PARJ(  3)=0.25            ! (P(us)/P(ud))/(P(s)/P(d))

    PARJ( 11)=0.25            ! Prob(Spin=1) of light meson
    PARJ( 12)=0.30            ! Prob(Spin=1) of strange meson

    PARJ( 21)=0.38            ! pT of Hadron in Fragm.: width
    PARJ( 23)=0.03            !    -"-                : (1)
    PARJ( 24)=2.50            !    -"-                : (2)

    PARJ( 33)=0.20            ! cut-off: String-Fragm <-> Cluster-Decay

    PARJ( 41)=1.13            ! LUND-FF: Parameter a
    PARJ( 42)=0.37            !   -"-  : Parameter b
    PARJ( 45)=0.80            !   -"-  : change Par. a for qq

    MSTJ( 12)=1               ! Switch baryon prod model

    MSTJ( 45)=4               ! #flavour in shower
    MSTJ(112)=4               ! #flavour in alphaS (nominal)
    MSTJ(113)=4               !     -"-            (min)
    MSTJ(114)=4               !     -"-            (max)

    CKIN(  1)=1.00            ! sqrt(s^hat): min

  end subroutine SetSomeDefaults_PY_Hermes


  !****************************************************************************
  !****s* CollTools/SetSomeDefaults_FR
  ! NAME
  ! subroutine SetSomeDefaults_FR
  ! PURPOSE
  ! Set some default values for FRITIOF and (old) PYTHIA/JETSET.
  !****************************************************************************
  subroutine SetSomeDefaults_FR
    use ID_translation, only: SetBruteForceMasses_FR

    COMMON/FRPARA1/KFR(20),VFR(20)
    integer KFR
    real VFR
    SAVE /FRPARA1/

    COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    integer MSTU,MSTJ
    real PARU,PARJ
    SAVE /LUDAT1/

    COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
    integer KCHG
    real PMAS,PARF,VCKM
    SAVE /LUDAT2/

    COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
    integer MDCY,MDME,KFDP
    real BRAT
    SAVE /LUDAT3/

    COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    real PARP,PARI
    SAVE /RYPARS/

    integer LUCOMP ! prototype

    logical, parameter :: lBruteForceMasses = .TRUE.


    !...PYTHIA/JETSET:

    MSTJ(21) = 0              ! particle decay on/off
    MSTU(12) = 0              ! writing of Pythia title page
                              ! (also affects FRITIOF!)
    MSTU(16) = 2              ! needed for GetLeading
    MSTU(21) = 1              ! don't stop after errors
    MSTU(22) = 0              ! #errors to be written
    MSTU(25) = 0              ! don't write warnings
    MSTU(26) = 0              ! #warnings to be written

    PARP(2) = 10.d0           ! lowest c.m. energy

                              ! forbid the decays of:
    MDCY(LUCOMP(111),1) = 0   ! pi0
!    MDCY(LUCOMP(311),1) = 0   ! K0

    !...FRITIOF:

    KFR( 2) = 1               ! multiple gluon emission
    KFR(11) = 0               ! write out init-message
    KFR( 7) = 1               ! PYTHIA usage

    if (lBruteForceMasses) call SetBruteForceMasses_FR

  end subroutine SetSomeDefaults_FR


  !****************************************************************************
  !****s* CollTools/SetVectorFromPYJETS
  ! NAME
  ! subroutine SetVectorFromPYJETS(Part, HardScaleQ2)
  ! PURPOSE
  ! Take the information from the PYJETS common block and set the output
  ! particle array Part accordingly.
  !
  ! INPUTS
  ! * type(particle),dimension(:) :: Part -- particle vector
  ! * real,                       :: HardScaleQ2 --
  !   Q2 of the production process
  !
  ! OUTPUT
  ! * integer, OPTIONAL :: iDiffrRho -- if given and a diffractive rho
  !   was generated, this returns the index
  ! * "Part" chaged
  !****************************************************************************
  subroutine SetVectorFromPYJETS (Part, HardScaleQ2, iDiffrRho)

    use hadronformation, only: SetJSVFormation
    use particledefinition
!    use propagation, only: checkVelo
    use CALLSTACK, only: TRACEBACK
    use ID_translation, only: KFtoBUU

    type(particle), dimension(:), intent(inout) :: Part
    real, intent(in)                            :: HardScaleQ2
    integer, intent(OUT), optional              :: iDiffrRho

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

    integer :: i1,i2
    integer :: ID0, IZ0
    !logical :: c

    if (N>size(Part)) then
       write(*,*) 'Error in SetVectorFromPYJETS: N>size', N, size(Part)
       write(*,*) 'You have to increase the size of the particle vector'
       write(*,*) 'in the calling routines!'
       call PYLIST(2)
       call TRACEBACK()
    end if

    call resetNumberGuess()

    if (present(iDiffrRho)) iDiffrRho = 0

    i2 = 0
    do i1=1,N
       call KFtoBUU (K(i1,2), ID0, IZ0)
       if (ID0.eq.0) cycle ! "unknown particles wont be propagated"
       i2 = i2+1

       Part(i2)%ID            = abs(ID0)
       Part(i2)%antiparticle  = (ID0.lt.0)
       Part(i2)%charge        = IZ0
       Part(i2)%mass          = P(i1,5)
       Part(i2)%momentum(1:3) = P(i1,1:3)
       Part(i2)%momentum(0)   = P(i1,4)

       Part(i2)%velocity = Part(i2)%momentum(1:3)/Part(i2)%momentum(0)
       !c=CheckVelo(Part(i2))

       call setNumberGuess(Part(i2))

!       Part(i2)%position = ...pos(coll)... ! already set
!       Part(i2)%productionTime = time      ! already set
!       Part(i2)%in_Formation = .TRUE.      ! already set

       if (present(iDiffrRho)) then
          if ((Part(i2)%ID.eq.103).and.(Part(i2)%charge.eq.0)&
               &.and.(K(i1,4).eq.2)) iDiffrRho = i2
       end if

       Part(i2)%scaleCS=1.
       Part(i2)%formationTime = -999

       call SetJSVFormation(Part(i2),i1,HardScaleQ2,K(i1,4)) ! if JetSetVec !!!

!!$ not set here:
!!$      real    :: lastCollisionTime=0.
!!$      real    :: offshellParameter=0.
!!$      real, dimension(1:2) :: formationTime=0.
!!$      real    :: perWeight=0.
!!$      real    :: coulombPotential=0.
!!$      integer   ::  number
!!$      integer,dimension(1:2) :: event=0
!!$      integer :: firstEvent=0
!!$      logical :: perturbative=.false.

    end do

  end subroutine SetVectorFromPYJETS



  !****************************************************************************
  !****f* CollTools/ConvertInPartPythia
  ! NAME
  ! function ConvertInPartPythia (ID) result (IDc)
  ! PURPOSE
  ! Reduce variety of incoming hadrons to some simpler manyfold
  ! (special for Pythia).
  !****************************************************************************
  function ConvertInPartPythia (ID) result (IDc)

    integer, intent(in) :: ID
    integer :: IDc

    integer :: aID

    aID = abs(ID)
    IDc = aID

    if (aID < 100) then ! BARYON

       select case (aID)
       case (2:31)
          IDc = 1   ! Nucleon
       case (34, 43:45, 50:52)
          IDc = 33  ! Sigma
       case (35:42, 46:49)
          IDc = 32  ! Lambda
       case (54)
          IDc = 53  ! Xi
       case (56:)
          write(*,*) 'ConvertInPartPythia: no charm! ',ID
          IDc = 0
       end select

       if (ID<0) IDc = -IDc

    else               ! MESON

       select case (ID)
       case (101:109,122)
          IDc = 101  ! pion
       case (112)
          IDc = 110  ! Kaon
       case (113)
          IDc = 111  ! Kbar
       case (114:121)
          write(*,*) 'ConvertInPartPythia: no charm! ',ID
          IDc = 0
       end select
    end if

  end function ConvertInPartPythia


  !****************************************************************************
  !****f* CollTools/ConvertInPartFritiof
  ! NAME
  ! function ConvertInPartFritiof (ID) result (IDc)
  ! PURPOSE
  ! Reduce variety of incoming hadrons to some simpler manyfold
  ! (special for Fritiof).
  !****************************************************************************
  function ConvertInPartFritiof (ID) result (IDc)

    integer, intent(in) :: ID
    integer :: IDc

    IDc = ID

    if (ID < 100) then ! BARYON

       select case (ID)
       case (1:31)
          IDc = 1
       case (34, 43:45, 50:52)
          IDc = 33
       case (35:42, 46:49)
          IDc = 32
       case (54)
          IDc = 53
       case (56:)
!          write(*,*) 'ConvertInPartFritiof: no charm! ',ID
       end select

    else               ! MESON

       select case (ID)
       case (101:109)
          IDc = 101
       case (112)
          IDc = 110
       case (113)
          IDc = 111
       case (114:)
!          write(*,*) 'ConvertInPartFritiof: no charm! ',ID
       end select
    end if

  end function ConvertInPartFritiof


  !****************************************************************************
  !****f* CollTools/ReduceCharge
  ! NAME
  ! function ReduceCharge(ID,IZ) result (IZc)
  ! PURPOSE
  ! Reduce charge of incoming Delta particles
  ! (e.g. Delta++ is mapped onto a Delta+).
  !****************************************************************************
  function ReduceCharge(ID,IZ) result (IZc)

    integer, intent(in)  :: ID, IZ
    integer :: IZc

    IZc = IZ
    if (abs(ID)<=31) then ! N, Delta, ...
      if (ID>0) then
        if (IZ==-1) then
          IZc = 0
        else if (IZ==2) then
          IZc = 1
        end if
      else
        if (IZ==1) then
          IZc = 0
        else if (IZ==-2) then
          IZc = -1
        end if
      end if
    end if

  end function ReduceCharge


  !****************************************************************************
  !****s* CollTools/CorrectChargeVector
  ! NAME
  ! subroutine CorrectChargeVector(Part,DeltaQ)
  ! PURPOSE
  ! Correct the Charge of the outgoing particles of the whole event.
  ! If deltaQ>0 one has to add positive charges.
  ! This is necessary, if ReduceCharge changed the charge of an
  ! incoming particle.
  !****************************************************************************
  subroutine CorrectChargeVector(Part,DeltaQ)
    use ID_translation, only: BUUKFDeltaQ
    use particleDefinition

    type(particle),dimension(:),intent(inout)   :: Part   ! particles
    integer, intent(inout)                      :: DeltaQ

    integer :: i,Q,dQMax, ID
    integer, ALLOCATABLE :: dQArray(:)

    allocate(dQArray(1:size(Part)))

    dQMax = 0
    do i=1,size(Part)
       if (Part(i)%ID <= 0) cycle

       ID = Part(i)%ID
       if (Part(i)%antiparticle) ID=-ID

       dQArray(i) = BUUKFDeltaQ(deltaQ, ID, Part(i)%charge)
       dQMax = dQMax + dQArray(i)
    end do

    if (abs(dQMax).lt.abs(deltaQ)) then
!       write (*,*) 'ERROR: CorrectChargeEvent: dQmax<deltaQ!'
       return
    end if

    Q = sign(1,deltaQ)

    do
       do i=1,size(Part)
          if (dQArray(i).ne.0) then
            Part(i)%charge = Part(i)%charge + Q
            dQArray(i)     = dQArray(i)     - Q
            deltaQ         = deltaQ         - Q
            if (deltaQ.eq.0) return
          end if
       end do
    end do

    deallocate(dQArray)

  end subroutine CorrectChargeVector


  !****************************************************************************
  !****s* CollTools/SetSomeDefaults_PY_Job
  ! NAME
  ! subroutine SetSomeDefaults_PY_Job
  !
  ! PURPOSE
  ! Read in the namelist "pythia" from the jobcard.
  !
  ! NOTES
  ! * only at the first call, the namelist is actually read in, in further
  !   calls the stored modifications are "replayed"
  ! * in order to keep the names in the namelist as usual, the arrays used
  !   in the common blocks are renamed.
  ! * not ALL Pythia parameters can be set, please see NAMELIST /pythia/
  !****************************************************************************
  subroutine SetSomeDefaults_PY_Job

    COMMON/PYDAT1/MSTU_D(200),PARU_D(200),MSTJ_D(200),PARJ_D(200)
    integer MSTU_D,MSTJ_D
    double precision PARU_D,PARJ_D
    SAVE /PYDAT1/

    COMMON/PYDAT2/KCHG_D(500,4),PMAS_D(500,4),PARF_D(2000),VCKM_D(4,4)
    integer KCHG_D
    double precision PMAS_D,PARF_D,VCKM_D
    SAVE /PYDAT2/

    COMMON/PYDAT3/MDCY_D(500,3),MDME_D(8000,2),BRAT_D(8000),KFDP_D(8000,5)
    integer MDCY_D,MDME_D,KFDP_D
    double precision BRAT_D
    SAVE /PYDAT3/

    COMMON/PYPARS/MSTP_D(200),PARP_D(200),MSTI_D(200),PARI_D(200)
    integer MSTP_D,MSTI_D
    double precision PARP_D,PARI_D
    SAVE /PYPARS/

    COMMON/PYSUBS/MSEL_D,MSELPD,MSUB_D(500),KFIN_D(2,-40:40),CKIN_D(200)
    integer MSEL_D,MSELPD,MSUB_D,KFIN_D
    double precision CKIN_D
    SAVE /PYSUBS/

    !**************************************************************************
    !****g* CollTools/MSTU
    ! SOURCE
    integer,dimension(200),  save :: MSTU
    ! PURPOSE
    ! Pythia array
    !**************************************************************************

    !**************************************************************************
    !****g* CollTools/MSTJ
    ! SOURCE
    integer,dimension(200),  save :: MSTJ
    ! PURPOSE
    ! Pythia array
    !**************************************************************************

    !**************************************************************************
    !****g* CollTools/MSTP
    ! SOURCE
    integer,dimension(200),  save :: MSTP
    ! PURPOSE
    ! Pythia array
    !**************************************************************************

    !**************************************************************************
    !****g* CollTools/MSTI
    ! SOURCE
    integer,dimension(200),  save :: MSTI
    ! PURPOSE
    ! Pythia array
    !**************************************************************************

    !**************************************************************************
    !****g* CollTools/PARU
    ! SOURCE
    real,   dimension(200),  save :: PARU
    ! PURPOSE
    ! Pythia array
    !**************************************************************************

    !**************************************************************************
    !****g* CollTools/PARJ
    ! SOURCE
    real,   dimension(200),  save :: PARJ
    ! PURPOSE
    ! Pythia array
    !**************************************************************************

    !**************************************************************************
    !****g* CollTools/PARP
    ! SOURCE
    real,   dimension(200),  save :: PARP
    ! PURPOSE
    ! Pythia array
    !**************************************************************************

    !**************************************************************************
    !****g* CollTools/PARI
    ! SOURCE
    real,   dimension(200),  save :: PARI
    ! PURPOSE
    ! Pythia array
    !**************************************************************************

    !**************************************************************************
    !****g* CollTools/MSEL
    ! SOURCE
    integer,                 save :: MSEL
    ! PURPOSE
    ! Pythia variable
    !**************************************************************************

    !**************************************************************************
    !****g* CollTools/CKIN
    ! SOURCE
    real,   dimension(200),  save :: CKIN
    ! PURPOSE
    ! Pythia array
    !**************************************************************************

    !**************************************************************************
    !****g* CollTools/PMAS
    ! SOURCE
    real,   dimension(500,4),save :: PMAS
    ! PURPOSE
    ! Pythia array
    !**************************************************************************

    !**************************************************************************
    !****g* CollTools/MDCY
    ! SOURCE
    integer,dimension(500,3),save :: MDCY
    ! PURPOSE
    ! Pythia array
    !**************************************************************************

    integer, save :: &
         & nMSTU,nMSTJ,nMSTP,nMSTI,nPARU,nPARJ,nPARP,nPARI,nCKIN,nMSEL,nPMAS,nMDCY
    integer,dimension(200), save :: &
         & iMSTU,iMSTJ,iMSTP,iMSTI,iPARU,iPARJ,iPARP,iPARI,iCKIN
    integer,dimension(500), save :: iPMAS, iMDCY

    logical, save :: NoChanges = .true.
    logical, save :: DoInit = .true.


    if (DoInit) call Init
    if (NoChanges) return

    call PlayChanges


  contains

    !**************************************************************************
    !****s* SetSomeDefaults_PY_Job/Init
    ! NAME
    ! subroutine Init
    !
    ! PURPOSE
    ! Read in the namelist "pythia" from the jobcard
    ! and find the differences to the original values.
    !**************************************************************************
    subroutine Init
      use output
      use hadronFormation, only: forceInitFormation

      integer :: i,ios

      !************************************************************************
      !****n* CollTools/pythia
      ! NAME
      ! NAMELIST /pythia/
      ! PURPOSE
      ! In addition to the default changes done by GiBUU, you can set additional
      ! parameters from the arrays:
      ! * MSEL
      ! * MSTU
      ! * MSTJ
      ! * MSTP
      ! * MSTI
      ! * PARU
      ! * PARJ
      ! * PARP
      ! * PARI
      ! * CKIN
      ! * PMAS
      ! * MDCY
      ! cf. PYTHIA documentation for details.
      !
      ! NOTES
      ! Only changes to the default values (i.e. values at the stage, when the
      ! namelist is read in) will be stored and replayed every time when
      ! necessary.
      !************************************************************************
      NAMELIST /pythia/ MSTU,MSTJ,MSTP,MSTI,PARU,PARJ,PARP,PARI,CKIN,MSEL,PMAS,MDCY

      ! (1) set to default:

      MSTU = MSTU_D
      MSTJ = MSTJ_D
      MSTP = MSTP_D
      MSTI = MSTI_D
      PARU = PARU_D
      PARJ = PARJ_D
      PARP = PARP_D
      PARI = PARI_D
      CKIN = CKIN_D
      MSEL = MSEL_D
      PMAS = PMAS_D
      MDCY = MDCY_D

      ! (2) read in namelist:

      call Write_ReadingInput('pythia',0)
      rewind(5)
      read(5,nml=pythia,iostat=ios)
      call Write_ReadingInput('pythia',0,ios)

      ! (3) find changes:

      nMSTU = 0
      nMSTJ = 0
      nMSTP = 0
      nMSTI = 0
      nPARU = 0
      nPARJ = 0
      nPARP = 0
      nPARI = 0
      nCKIN = 0
      nMSEL = 0
      nPMAS = 0
      nMDCY = 0

      do i=1,200
         if (MSTU(i).ne.MSTU_D(i)) then
            if (i.eq.12) then
               write(*,*) 'Changes of MSTU(12) are ignored'
            else
               nMSTU = nMSTU+1
               iMSTU(nMSTU) = i
            end if
         end if

         if (MSTJ(i).ne.MSTJ_D(i)) then
            nMSTJ = nMSTJ+1
            iMSTJ(nMSTJ) = i
         end if

         if (MSTP(i).ne.MSTP_D(i)) then
            nMSTP = nMSTP+1
            iMSTP(nMSTP) = i
         end if

         if (MSTI(i).ne.MSTI_D(i)) then
            nMSTI = nMSTI+1
            iMSTI(nMSTI) = i
         end if

         if (PARU(i).ne.PARU_D(i)) then
            nPARU = nPARU+1
            iPARU(nPARU) = i
         end if

         if (PARJ(i).ne.PARJ_D(i)) then
            nPARJ = nPARJ+1
            iPARJ(nPARJ) = i
         end if

         if (PARP(i).ne.PARP_D(i)) then
            nPARP = nPARP+1
            iPARP(nPARP) = i
         end if

         if (PARI(i).ne.PARI_D(i)) then
            nPARI = nPARI+1
            iPARI(nPARI) = i
         end if

         if (CKIN(i).ne.CKIN_D(i)) then
            nCKIN = nCKIN+1
            iCKIN(nCKIN) = i
         end if
      end do

      if (MSEL.ne.MSEL_D) then
         nMSEL = 1
      end if

      do i=1,500
         if (SUM(PMAS(i,:)-PMAS_D(i,:))/=0) then
            nPMAS = nPMAS+1
            iPMAS(nPMAS) = i
         end if
      end do

      do i=1,500
         if (MDCY(i,1).ne.MDCY_D(i,1)) then
            nMDCY = nMDCY+1
            iMDCY(nMDCY) = i
         end if
         if ((MDCY(i,2).ne.MDCY_D(i,2)).or.((MDCY(i,3).ne.MDCY_D(i,3)))) then
            write(*,*) 'Entries MDCY(:,2) or MDCY(:,3) will be ignored!'
         end if

      end do

      noChanges= MAX(nMSTU,nMSTJ,nMSTP,nMSTI,nPARU,nPARJ,nPARP,nPARI,nCKIN,nMSEL,nPMAS,nMDCY)==0

      if (NoChanges) then
         write(*,*) 'No Changes indicated by namelist'
      else
         call WriteChanges
      end if


      call Write_ReadingInput('pythia',1)
      DoInit = .false.

      call PYLOGO

      call forceInitFormation

    end subroutine Init


    !**************************************************************************
    !****s* SetSomeDefaults_PY_Job/WriteChanges
    ! NAME
    ! subroutine WriteChanges
    !
    ! PURPOSE
    ! Write the changes to stdout.
    !**************************************************************************
    subroutine WriteChanges

      integer :: i

!      call WriteValues()

      if (nMSEL>0) then
         write(*,1001) 'MSEL',MSEL_D,MSEL
      end if

      do i=1,nMSTU
         write(*,1002) 'MSTU',iMSTU(i),MSTU_D(iMSTU(i)),MSTU(iMSTU(i))
      end do
      do i=1,nMSTJ
         write(*,1002) 'MSTJ',iMSTJ(i),MSTJ_D(iMSTJ(i)),MSTJ(iMSTJ(i))
      end do
      do i=1,nMSTP
         write(*,1002) 'MSTP',iMSTP(i),MSTP_D(iMSTP(i)),MSTP(iMSTP(i))
      end do
      do i=1,nMSTI
         write(*,1002) 'MSTI',iMSTI(i),MSTI_D(iMSTI(i)),MSTI(iMSTI(i))
      end do

      do i=1,nPARU
         write(*,1003) 'PARU',iPARU(i),PARU_D(iPARU(i)),PARU(iPARU(i))
      end do
      do i=1,nPARJ
         write(*,1003) 'PARJ',iPARJ(i),PARJ_D(iPARJ(i)),PARJ(iPARJ(i))
      end do
      do i=1,nPARP
         write(*,1003) 'PARP',iPARP(i),PARP_D(iPARP(i)),PARP(iPARP(i))
      end do
      do i=1,nPARI
         write(*,1003) 'PARI',iPARI(i),PARI_D(iPARI(i)),PARI(iPARI(i))
      end do

      do i=1,nCKIN
         write(*,1003) 'CKIN',iCKIN(i),CKIN_D(iCKIN(i)),CKIN(iCKIN(i))
      end do

      do i=1,nPMAS
         write(*,1004) 'PMAS',iPMAS(i),PMAS_D(iPMAS(i),:),PMAS(iPMAS(i),:)
      end do

      do i=1,nMDCY
         write(*,1005) 'MDCY',iMDCY(i),MDCY_D(iMDCY(i),1),MDCY(iMDCY(i),1)
      end do

1001  FORMAT ("  ",A,"       "," : ",i14," -> ",i14)
1002  FORMAT ("  ",A,"(",i3,")  "," : ",i14," -> ",i14)
1003  FORMAT ("  ",A,"(",i3,")  "," : ",f14.5," -> ",f14.5)
1004  FORMAT ("  ",A,"(",i3,")  "," : ",4f14.5," -> ",4f14.5)
1005  FORMAT ("  ",A,"(",i3,",1)"," : ",i14," -> ",i14)

    end subroutine WriteChanges


    !**************************************************************************
    !****s* SetSomeDefaults_PY_Job/PlayChanges
    ! NAME
    ! subroutine WriteChanges
    !
    ! PURPOSE
    ! Play the changes.
    !**************************************************************************
    subroutine PlayChanges

      integer :: i

!!$      write(*,*) 'Play Changes:'
!!$      call WriteChanges

      if (nMSEL>0)  MSEL_D = MSEL


      do i=1,nMSTU
         MSTU_D(iMSTU(i))=MSTU(iMSTU(i))
      end do
      do i=1,nMSTJ
         MSTJ_D(iMSTJ(i))=MSTJ(iMSTJ(i))
      end do
      do i=1,nMSTP
         MSTP_D(iMSTP(i))=MSTP(iMSTP(i))
      end do
      do i=1,nMSTI
         MSTI_D(iMSTI(i))=MSTI(iMSTI(i))
      end do

      do i=1,nPARU
         PARU_D(iPARU(i))=PARU(iPARU(i))
      end do
      do i=1,nPARJ
         PARJ_D(iPARJ(i))=PARJ(iPARJ(i))
      end do
      do i=1,nPARP
         PARP_D(iPARP(i))=PARP(iPARP(i))
      end do
      do i=1,nPARI
         PARI_D(iPARI(i))=PARI(iPARI(i))
      end do

      do i=1,nCKIN
         CKIN_D(iCKIN(i))=CKIN(iCKIN(i))
      end do

      do i=1,nPMAS
         PMAS_D(iPMAS(i),:)=PMAS(iPMAS(i),:)
      end do

      do i=1,nMDCY
         MDCY_D(iMDCY(i),1)=MDCY(iMDCY(i),1)
      end do

    end subroutine PlayChanges


    !**************************************************************************
    !****s* SetSomeDefaults_PY_Job/WriteValues
    ! NAME
    ! subroutine WriteValues
    !
    ! PURPOSE
    ! Write all the values to stdout.
    !**************************************************************************
!     subroutine WriteValues
!
!       integer :: i
!
!       write(*,1001) 'MSEL',MSEL
!
!       do i=1,200
!          write(*,1002) 'MSTU',i,MSTU(i)
!       end do
!       do i=1,200
!          write(*,1002) 'MSTJ',i,MSTJ(i)
!       end do
!       do i=1,200
!          write(*,1002) 'MSTP',i,MSTP(i)
!       end do
!       do i=1,200
!          write(*,1002) 'MSTI',i,MSTI(i)
!       end do
!
!       do i=1,200
!          write(*,1003) 'PARU',i,PARU(i)
!       end do
!       do i=1,200
!          write(*,1003) 'PARJ',i,PARJ(i)
!       end do
!       do i=1,200
!          write(*,1003) 'PARP',i,PARP(i)
!       end do
!       do i=1,200
!          write(*,1003) 'PARI',i,PARI(i)
!       end do
!
!       do i=1,200
!          write(*,1003) 'CKIN',i,CKIN(i)
!       end do
!
! !!$      do i=1,500
! !!$         write(*,1004) 'PMAS',i,PMAS(i),PMAS(iPMAS(i),:)
! !!$      end do
! !!$
! !!$      do i=1,500
! !!$         write(*,1005) 'MDCY',i,MDCY(i),MDCY(iMDCY(i),1)
! !!$      end do
!
!
!       stop
!
! 1001  FORMAT ("  ",A,"       "   ," = ",i14)
! 1002  FORMAT ("  ",A,"(",i3,")  "," = ",i14)
! 1003  FORMAT ("  ",A,"(",i3,")  "," = ",f14.5)
! !1004  FORMAT ("  ",A,"(",i3,")  "," = ",4f14.5)
! !1005  FORMAT ("  ",A,"(",i3,",1)"," = ",i14)
!
!     end subroutine WriteValues

  end subroutine SetSomeDefaults_PY_Job


  !****************************************************************************
  ! cf. interface ConvertToAnti
  !****************************************************************************
  subroutine ConvertToAnti_1(part)
    use particleDefinition
    use IdTable, only: isMeson, isBaryon, getAntiMeson, photon

    type(particle),intent(inout)   :: part

    if (part%ID==0) return
    if (part%ID==photon) return

    if (isMeson(part%ID)) then
       call getAntiMeson(part%ID,part%Charge, part%ID,part%Charge)
    else if (isBaryon(part%ID)) then
       part%Charge       = -part%Charge
       part%antiparticle = .not.part%antiparticle
    else
       write(*,*) 'OOps',part%ID,part%Charge,part%antiparticle
       stop
    end if
  end subroutine ConvertToAnti_1

  !-------------------------------------------------------------------------

  subroutine ConvertToAnti_11(part,partNEW)
    use particleDefinition
    type(particle),intent(in)   :: part
    type(particle),intent(out)  :: partNEW

    partNEW = part
    call ConvertToAnti_1(partNEW)
  end subroutine ConvertToAnti_11

  !-------------------------------------------------------------------------

  subroutine ConvertToAnti_2(parts)
    use particleDefinition
    type(particle),dimension(:),intent(inout)   :: parts
    integer :: i

    do i=lbound(parts,dim=1),ubound(parts,dim=1)
       if (parts(i)%ID > 0) call ConvertToAnti_1(parts(i))
    end do
  end subroutine ConvertToAnti_2

  !****************************************************************************
  !****************************************************************************


  !****************************************************************************
  !****f* CollTools/PythiaCKIN
  ! NAME
  ! real function PythiaCKIN(i)
  ! PURPOSE
  ! Return the value from the common block.
  !****************************************************************************
  real function PythiaCKIN(i)
    integer :: i

    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
    integer MSEL,MSELPD,MSUB,KFIN
    double precision CKIN
    SAVE /PYSUBS/

    PythiaCKIN = CKIN(i)

  end function PythiaCKIN


  !****************************************************************************
  !****f* CollTools/PythiaMSTP
  ! NAME
  ! integer function PythiaMSTP(i)
  ! PURPOSE
  ! Return the value from the common block.
  !****************************************************************************
  integer function PythiaMSTP(i)
    integer :: i

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    PythiaMSTP = MSTP(i)

  end function PythiaMSTP


  !****************************************************************************
  !****s* CollTools/CheckUndecayedString
  ! NAME
  ! subroutine CheckUndecayedString()
  ! PURPOSE
  ! Check whether some undecayed quark etc survived in the JETSET output.
  ! Should not happen, but occurs sometimes.
  ! OUTPUT
  ! * MINT(51) = 2, if problem occured
  !****************************************************************************
  subroutine CheckUndecayedString()

    !...common blocks:

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

    COMMON/PYINT1/MINT(400),VINT(400)
    integer MINT
    double precision VINT
    SAVE /PYINT1/

    integer :: i
    logical :: err

    if (MINT(51).eq.2) return

    err = .false.
    do i=1,N
       if (K(i,1).lt.10) then ! stable particle
          if (abs(K(i,2)).lt.10) err = .true. ! quark
          if (K(i,2).eq.21) err = .true. ! gluon
       end if
    end do

    if (err) then
       MINT(51) = 2
       write(*,*) 'MESSAGE: Jetset problem, Event rejected!'
       call PYLIST(2)
    end if

  end subroutine CheckUndecayedString


end module CollTools
