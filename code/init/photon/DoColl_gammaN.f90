!******************************************************************************
!****m* /Coll_gammaN
! NAME
! module Coll_gammaN
!
! PURPOSE
! ...
!******************************************************************************
module Coll_gammaN

  use particleDefinition
  use CollTools
  use hadronFormation
  use PyVP
  use eN_eventDefinition

  implicit none
  private

  logical, save :: DoColl_gammaN_verbose = .false.

  public :: DoColl_gammaN_Py, DoColl_gammaN_Fr
  public :: DoColl_gammaN_Toy
  public :: DoColl_gammaN_verbose

contains

  !****************************************************************************
  !****s* Coll_gammaN/DoColl_gammaN_Py
  ! NAME
  ! subroutine DoColl_gammaN_Py(eNev,outPart,flagOK, rVMD, DoDifr, Cross,EventClass,MinW)
  !
  ! PURPOSE
  ! generate a high energy Photon event with PYTHIA
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev     -- electron nucleon kinematics
  ! * real, dimension(1:4)        :: rVMD     -- scaling of VMD XS (rho, omega, phi, J/psi)
  ! * logical                     :: DoDifr   -- Do diffractive events or not
  ! * real, OPTIONAL              :: minW     -- if given, return immediately if W<minW
  !
  ! OUTPUT
  ! * type(particle),dimension(:) :: outPart  -- outgoing particles
  ! * real, dimension(0:4)        :: Cross    -- cross sections for the different
  !   event classes VMD, direct, anomalous, DIS and their sum
  ! * integer                     :: EventClass -- reported class of this event
  !   (cf. PYTHIA encoding of events)
  ! * logical                     :: flagOK   -- .TRUE., if everything okay
  !
  ! NOTES
  ! The returned cross section is sigma_T(y,Q^2) + eps * sigma_L(y,Q^2) in mb.
  !
  ! This follows the definitions of Christy&Bosted, so one has to
  ! multiply this returned cross section with
  !   \alpha/(2\pi) K/(Q^2 \nu^2) (E E')/\pi cT
  ! with cT ~ 1+(1-y)^2 in order to get the value of
  !   d\sigma/(dE'd\cos\theta) in mb/GeV
  !
  ! PYTHIA fails primarily below W=2GeV. Therefore it does not make sense to
  ! reduce the hard wired cut in this routine. You may increase it with
  ! the parameter MinW; if you *REALLY* know what you are doing, you may
  ! also decrease it this way, as e.g. for DIS events below W=2GeV.
  !
  !****************************************************************************
  subroutine DoColl_gammaN_Py(eNev,outPart,flagOK, rVMD, DoDifr, Cross,EventClass,MinW)

    use CollGetLeading
    use PIL_rhoDiff
    use output
    use eN_eventDefinition
    use eN_event

    type(electronNucleon_event),intent(in)   :: eNev
    type(particle),dimension(:),intent(inout):: outPart  ! outgoing particles
    logical,                    intent(out)  :: flagOK
    real, dimension(1:4),       intent(in)   :: rVMD
    logical,                    intent(in)   :: DoDifr
    real, dimension(0:4),       intent(out)  :: Cross
    integer,                    intent(out)  :: EventClass
    real,                       intent(IN), optional :: minW


    ! common blocks

    COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
    integer KCHG
    double precision PMAS,PARF,VCKM
    SAVE /PYDAT2/

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    COMMON/PYINT1/MINT(400),VINT(400)
    integer MINT
    double precision VINT
    SAVE /PYINT1/

    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
    integer MSEL,MSELPD,MSUB,KFIN
    double precision CKIN
    SAVE /PYSUBS/

    ! for temprorary use:

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

    ! internal variables:

    real :: theta, phi
    integer :: iEv,nEv
    integer :: i, iDiffrRho!,nPion
    integer :: MSTP03,MSTP81,MSTP82
    real :: CKIN05

    real, dimension(0:3) :: recoil

    real  :: nu,Q2,W,Wfree


    ! set default output

    outPart%ID = 0
    Cross = 0.
    flagOK = .FALSE.
    EventClass = 0

    ! set up PYTHIA



    call eNeV_GetKinV(eNev, nu,Q2,W,Wfree)


    ! only high energy photons:
    ! (we allow for minW < W < 2GeV, but not for W < 2GeV alone !!!)

    if (.not.present(minW)) then
       if (Wfree.lt.2.0) then
          write(*,*) 'Error in DoColl_gammaN (W < 2 GeV) : W=',Wfree
!          return
          stop
       end if
    else
       if (minW.lt.2.0) then
          if (Wfree.lt.minW) then
             write(*,*) 'Error in DoColl_gammaN (W < minW <2 GeV) : W=',Wfree,minW
             stop
          end if
       else
          if (Wfree.lt.minW) then
             write(*,*) 'Error in DoColl_gammaN (W < minW >2 GeV) : W=',Wfree,minW
             return
          end if
       end if
    end if


    call SetSomeDefaults_PY
    call ScaleVMD(rVMD)


    ! adjust the pThat cut off according W:

    CKIN05 = CKIN(5)
    if (2*CKIN(5).ge.Wfree) CKIN(5) = Wfree/2-1e-3

    ! set some values, Pythia 6.4 can not live without !!!!

    MSTP03 = MSTP(3)
    MSTP81 = MSTP(81)
    MSTP82 = MSTP(82)

    if (MSTP(182).ge.400) then ! Pythia v6.4
       MSTP(81) = 1      ! selection of main routine
       MSTP(82) = 1      ! structure of multiple interactions
       MSTP(3)  = 1      ! selection of Lambda in alpha_S
    end if

    ! calculate cross sections:

    if (DoDifr) then
       nEv = 100              ! KOG: check number!
    else
       nEv = 100              ! KOG: check number!
    end if


    call InitPythia(eNev,.FALSE.,0,DoDifr)

    do iEv=1,nEv
       if (useJetSetVec) call GetJetsetVecINIT
       call PYEVNT

!!$       if (MINT(51).eq.2) then
!!$          write(*,*) 'N=',N
!!$          call PYGIVE('VINT(1)=')
!!$          call PYGIVE('VINT(3)=')
!!$          call PYGIVE('VINT(4)=')
!!$          call PYGIVE('VINT(67)=')
!!$          call PYGIVE('VINT(68)=')
!!$          call PYGIVE('VINT(69)=')
!!$          call PYGIVE('VINT(70)=')
!!$          call PYGIVE('MINT(17)=')
!!$          call PYGIVE('MINT(18)=')
!!$          call PYGIVE('MINT(101)=')
!!$          call PYGIVE('MINT(102)=')
!!$          call PYGIVE('MINT(103)=')
!!$          call PYGIVE('MINT(104)=')
!!$          call PYGIVE('MINT(107)=')
!!$          call PYGIVE('MINT(108)=')
!!$          call PYGIVE('MINT(122)=')
!!$
!!$       endif




    end do

    call CollectXS_Class(Cross, 0)

    ! generate THE EVENT:

    call InitPythia(eNev,.TRUE.,0,DoDifr)

    if (useJetSetVec) call GetJetSetVecINIT

!!$    select case (inPart%number)
!!$    case (100010)
!!$       verb = 1
!!$    case default
!!$       verb = 0
!!$    end select

    call PYEVNT

    CKIN(5) = CKIN05 ! reset PYTHIA defaults
    MSTP(3) = MSTP03
    MSTP(81)= MSTP81
    MSTP(82)= MSTP82

    if (MINT(51).eq.2) then
!!$       write(*,*) 'N=',N
!!$       call PYGIVE('VINT(1)=')
!!$       call PYGIVE('VINT(3)=')
!!$       call PYGIVE('VINT(4)=')
!!$       call PYGIVE('VINT(67)=')
!!$       call PYGIVE('VINT(68)=')
!!$       call PYGIVE('VINT(69)=')
!!$       call PYGIVE('VINT(70)=')
!!$       call PYGIVE('MINT(17)=')
!!$       call PYGIVE('MINT(18)=')
!!$       call PYGIVE('MINT(101)=')
!!$       call PYGIVE('MINT(102)=')
!!$       call PYGIVE('MINT(103)=')
!!$       call PYGIVE('MINT(104)=')
!!$       call PYGIVE('MINT(107)=')
!!$       call PYGIVE('MINT(108)=')
!!$       call PYGIVE('MINT(122)=')

       N = 0
       return ! -> FAILURE
    end if

    call MarkLepton

!    call PYLIST(2)

    if (useJetSetVec) then

!!$       if (verb>0) then
!!$          call SFREPS_Write(1)
!!$          call SFREPS_Write(2)
!!$          call SFREPS_Write(3)
!!$       end if
!!$       write(*,*) '##### Part=',inPart%number
!!$       call PYGIVE('MSTP(14)=')
!!$       call PYGIVE('MSTI(1)=')
!!$       call PYGIVE('MINT(1)=')
!!$       call PYGIVE('MSTI(9)=')
!!$       call PYGIVE('MINT(103)=')
!!$       call PYGIVE('MINT(107)=')
!!$       call PYGIVE('MINT(122)=')
!!$       call PYGIVE('MINT(123)=')
!!$       call PYGIVE('MINT(17)=')
!!$       call PYGIVE('MINT(18)=')
!!$
!!$
!!$       call PYLIST(2)

       call GetJetsetVec(.TRUE.)
       call GetJetsetVecCheckT(-1d-5)
       call GetJetsetVecPYEDIT

    end if

    call GetLeading_PY         ! find leading particles

!    call PYLIST(2)
!    stop

    if (DoColl_gammaN_verbose) call PYLIST(2)

    !...Rotate and boost the whole event to final system:

    phi = atan2(eNev%pcm(2),eNeV%pcm(1))
    theta = atan2(sqrt(eNev%pcm(1)**2+eNev%pcm(2)**2),eNev%pcm(3))
    call PYROBO(1,N, 0.0,eNev%phiLepTON, 0.0,0.0,0.0)
    call PYROBO(1,N, theta,phi, eNev%betacm(1),eNev%betacm(2),eNev%betacm(3))

    if (DoColl_gammaN_verbose) call PYLIST(2)

!    call PYLIST(2)

    call PYEDIT(1)            ! clean up event list

    if (DoColl_gammaN_verbose) call PYLIST(2)
    if (useJetSetVec) then

       call GetJetsetVecPYROBO(0.0,  eNev%phiLepTON, 0.0,0.0,0.0)
       call GetJetsetVecPYROBO(theta,phi, eNev%betacm(1),eNev%betacm(2),eNev%betacm(3))

       if (DoColl_gammaN_verbose) then
          write(*,*) '##### Part=',eNev%nucleon_free%number
          call GetJetSetVec_List(6,1,N)
       end if
!       stop
    end if

    !...Copy Particles to ouput-vector

    call SetVectorFromPYJETS(outPart, Q2, iDiffrRho)

    if (DoColl_gammaN_verbose) call WriteParticle(6,1,outPart)
    if (DoColl_gammaN_verbose) call WriteParticle(6,1,0,eNev%nucleon_free)

    if (iDiffrRho/=0) then

!!$!       call WriteParticle(6,1,0,inPart)

!!$       write(*,*) '###############'
!!$       write(*,*) 'N=',N,MSTI(9)
!!$       call WriteParticle(6,0,iDiffrRho,outPart(iDiffrRho))
!!$       if (N.eq.2) call WriteParticle(6,0,3-iDiffrRho,outPart(3-iDiffrRho))
!!$       write(*,*) '###############'

       !...Adjust gamma*N -> rho0N VMD-XS by a 'fit by eye' at low W:

       if (N.eq.2) then
          if ((Wfree.ge.2.0).and.(Wfree.lt.2.3)) then
             Cross(0) = Cross(0)*(1.0-0.5*(2.3-Wfree))
             Cross(1) = Cross(1)*(1.0-0.5*(2.3-Wfree))
          end if
       end if

       recoil = 0
       do i=1,N
          recoil(1:3) = recoil(1:3) + P(i,1:3)
          recoil(0)   = recoil(0)   + P(i,4)
       end do
       recoil = recoil  - outPart(iDiffrRho)%momentum(0:3)
!!$       write(*,*) recoil
!!$       write(*,*) PARI(113)

       call PIL_rhoDiffractive_PUT(outPart(iDiffrRho)%number,&
            & .true.,real(PARI(113)-1.0),recoil)

!             stop
    end if

    !...remember event class of this event:


    select case (MSTP(14))
    case (30)
       EventClass = MSTI(9)
    case (26)
       EventClass = 4
    case default
       write(*,*) 'wrong MSTP(14). STOP!', MSTP(14)
       stop
    end select

    flagOK = .TRUE.

!!$    if (N.eq.3) then
!!$       nPion = 0
!!$       do i=1,3
!!$          if (K(i,2).eq.111) nPion=nPion+1
!!$          if (abs(K(i,2)).eq.211) nPion=nPion+1
!!$       enddo
!!$
!!$       if (nPion.eq.2) then
!!$          call PYLIST(2)
!!$          call GetJetSetVec_List(6,1,N)
!!$          write(*,*)
!!$       endif
!!$    end if

  end subroutine DoColl_gammaN_Py


  !****************************************************************************
  !****s* Coll_gammaN/DoColl_gammaN_Fr
  ! NAME
  ! subroutine DoColl_gammaN_Fr(inPart,outPart,flagOK, W,Q2,eps, pcm,beta, iTyp)
  !
  ! PURPOSE
  ! generate a high energy Photon event with FRITIOF (=VMD!)
  !
  ! With FRITIOF we can only do VMD events. "iTyp" selects the type of
  ! vector meson.
  !
  ! INPUTS
  ! * type(particle)              :: inPart   -- incoming nucleon
  ! * real                        :: W        -- incoming photon (W)
  ! * real                        :: Q2       -- incoming photon (Q^2)
  ! * real                        :: eps      -- incoming photon (epsilon)
  ! * real, dimension(0:3)        :: pcm      -- boost vector
  ! * real, dimension(1:3)        :: beta     -- boost vector
  ! * integer                     :: iTyp     -- 1:4 = rho, omega, phi, J/Psi
  !
  ! OUTPUT
  ! * type(particle),dimension(:) :: outPart  -- outgoing particles
  ! * logical                     :: flagOK   -- .TRUE., if everything okay
  !
  ! NOTES
  ! With Q2=Q^2 > 0 we have a projectile with imaginary mass.
  ! The VMD prescription holds for Q^2=0 with m_V=0.
  !
  ! Here, with Q2>0,
  ! the mass of the projectile (mass1) can be choosen to be =0 or
  ! to be massArr(iTyp) = /0.76850, 0.78194, 1.01940, 3.09688/.
  !
  ! Events of the kind V N -> pi0 pi0 R where both pi0 come from a cluster->2
  ! decay are excluded (by returning "flagOK=.false.")
  !****************************************************************************
  subroutine DoColl_gammaN_Fr(inPart,outPart,flagOK, W,Q2,eps, pcm,beta, iTyp)

    use CollGetLeading

    type(particle),             intent(in)   :: inPart   ! incoming nucleon
    type(particle),dimension(:),intent(inout):: outPart  ! outgoing particles
    logical,                    intent(out)  :: flagOK
    real,                       intent(in)   :: W
    real,                       intent(in)   :: Q2
    real,                       intent(in)   :: eps
    real, dimension(0:3),       intent(in)   :: pcm
    real, dimension(1:3),       intent(in)   :: beta
    integer,                    intent(in)   :: iTyp

!...common blocks:

    COMMON/FRCODES/IPT(2),PACD(27),NNUC(27),NPROT(27),KCD(27),RO1(27,2),EXMA(9,2)
    integer IPT, NNUC, NPROT, KCD
    real RO1, EXMA
    character*4 PACD
    SAVE /FRCODES/

    COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    integer MSTU,MSTJ
    real PARU,PARJ
    SAVE /LUDAT1/

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

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

!...prototypes:

    real ULMASS

!...internal variables:

    ! rho, omega, phi, J/psi:
    integer, dimension(4), parameter :: iKF1Arr = (/     113,     223,     333,     443 /)
    real,    dimension(4), parameter :: mass1Arr= (/ 0.76850, 0.78194, 1.01940, 3.09688 /)

    integer :: iKF1, iKF2
    real    :: mass1, mass2, plab
    integer :: iTry
    real    :: theta, phi

    integer :: nPion,i
    integer, dimension(3) :: iPion

!...select incoming particles:

!    N = 0


    flagOK = .FALSE.
    outPart%ID = 0 ! reset outgoing particles

    if ((iTyp.lt.1) .or. (iTyp.gt.4)) then
       write(*,*) 'DoColl_GammaN_Fr: iTyp=',iTyp,' <1 or >4. STOP.'
       stop
    end if

    iKF1 = iKF1Arr(iTyp)    ! Beam is vector meson

    select case (InPart%charge)
    case (1)
       iKF2 = 2212          ! Target is proton
    case (0)
       iKF2 = 2112          ! Target is neutron
    case default
       write(*,*) 'wrong charge. STOP!'
       stop
    end select


     mass1 = mass1Arr(iTyp)   ! Meson Mass = Peak Mass !!!
!     mass1 = 0.0               ! Meson Mass = 0 !!!
    mass2 = ULMASS(iKF2)

    if (W.lt. mass1+mass2) then
       write(*,*) 'DoColl_GammaN_Fr: W=',W,' < m_1+m2 =',mass1,'+',mass2,'=',mass1+mass2,'. STOP.'
       stop
    end if

!...calculate Plab for fixed target event:

!    plab = ((W**2-mass1Arr(iTyp)**2-mass2**2)/(2.0*mass2))**2 - mass1Arr(iTyp)**2
!   if (plab.lt.0.0) then
!       write(*,*) 'DoColl_GammaN_Fr: plab^2=',plab,' < 0 . STOP. [try]'
!       return
!       !stop
!    endif


    plab = ((W**2-mass1**2-mass2**2)/(2.0*mass2))**2 - mass1**2
    if (plab.lt.0.0) then
       write(*,*) 'DoColl_GammaN_Fr: plab^2=',plab,' < 0 . STOP.'
       stop
    end if
    plab = sqrt(plab)

!...set up PYTHIA/JETSET/FRITIOF::

    call SetSomeDefaults_FR

! -- Projectile: Vector meson --

    KCD(1) = iKF1
    NPROT(1) = 0
    NNUC(1) = 1
    EXMA(1,1) = mass1             ! Fritiof sets mass to default !!!
    EXMA(1,2) = mass1

! -- Target: proton/neutron --

    KCD(2) = iKF2
    NPROT(2) = InPart%charge  ! proton/neutron
    NNUC(2) = 1
    EXMA(2,1) = mass2         ! Fritiof sets mass to default !!!
    EXMA(2,2) = mass2


!    write(*,*) '===',mass1,mass2


!...Start the Event Loop

    iTry = 0
    do
       outPart%ID = 0 ! reset outgoing particles

       iTry=iTry+1
       if (iTry.ge.100) then
          write(*,*) 'DoColl_GammaN_Fr: itry=',iTry
          return
       end if

       if (useJetSetVec) call GetJetsetVecINIT

!...Generate THE EVENT:

       call FRHADHAD('CMS','NEW1','NEW2',W, mass1,mass2)

       if (MSTU(24).ne.0) then
          if (iTRY.eq.99) write(*,*) 'DoColl_GammaN_Fr: MSTU(24).ne.0)',MSTU(24)
          cycle
       end if

       if (useJetSetVec) then
          call GetJetsetVec(.FALSE.)
!          call PYLIST(2)
!          call GetJetSetVec_List(6,1,N)

          call GetJetsetVecCheckT(-1d-5)

          call GetJetsetVecPYEDIT
       end if

       call GetLeading_FR ! find leading particles
       call LUEDIT(1)     ! clean up event list

!...Rotate and boost the whole event to final system:

       phi = atan2(pcm(2),pcm(1))
       theta = atan2(sqrt(pcm(1)**2+pcm(2)**2),pcm(3))
       call LUROBO(theta, phi, beta(1),beta(2),beta(3))

       if (useJetSetVec) then
          call GetJetsetVecPYROBO(theta,phi, beta(1),beta(2),beta(3))

!         call PYLIST(2)
!         call GetJetSetVec_List(6,1,N)
!         stop
       end if

!...Copy Particles to ouput-vector

       call SetVectorFromPYJETS(outPart, Q2)

!...exit the loop
       exit

    end do


    flagOk = .TRUE.

    if (N.eq.3) then
       nPion = 0
       do i=1,3
          if (K(i,2).eq.111) then ! its pi0
             nPion=nPion+1
             iPion(nPion) = i
          end if
       end do

       if (nPion.eq.2) then
          if (EArr(3,iPion(1)).eq.5 .and. EArr(3,iPion(2)).eq.5) then
             flagOk = .FALSE.

!!$             call PYLIST(2)
!!$             call GetJetSetVec_List(6,1,N)
!!$             write(*,*)
          end if
       end if
    end if


  end subroutine DoColl_gammaN_Fr

  !****************************************************************************
  !****s* Coll_gammaN/DoColl_gammaN_Toy
  ! NAME
  ! subroutine DoColl_gammaN_Toy(eNev,flagOK,outPart)
  !
  ! PURPOSE
  ! generate a high energy Photon event according to a toy model
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev     -- electron nucleon kinematics
  !
  ! OUTPUT
  ! * type(particle),dimension(:) :: outPart  -- outgoing particles
  ! * logical                     :: flagOK   -- .TRUE., if everything okay
  !
  ! NOTES
  ! * This generates a single pion according to a very simple toy model.
  !   There are no conservations: baryon number, energy, momentum etc is
  !   violated
  !
  ! * The full glory of "type(electronNucleon_event)" is not used up to now
  !****************************************************************************
  subroutine DoColl_gammaN_Toy(eNev,flagOK,outPart)

    use constants, only: twoPi, mPi, mN
    use random, only: rn, rnExp
    use eN_event

    type(electronNucleon_event),intent(in)   :: eNev
    type(particle),dimension(:),intent(inout):: outPart
    logical,                    intent(out)  :: flagOK

    real                          :: W,Wfree
    real                          :: Q2
    real                          :: eps

    real :: pT,pT2,phi,z,pz,nu,EE,pp
    integer :: i!, iQ
    real :: phiB, thetaB ! boost angles

    double precision MP_P ! prototype

    flagOK = .FALSE.

    call eNeV_GetKinV(eNev, eps,nu,Q2,W,Wfree)

    outPart%ID = 0 ! reset outgoing particles
    outPart%number = 0
    outPart%antiparticle=.false.
    outPart%perWeight = 1.0

    outPart(1)%ID = 101
    outPart(1)%mass = mPi

!    outPart(2)%ID = 1
!    outPart(2)%mass = mN

    ! select charges:
!!$
!!$    if (rn(0).gt.0.5) then
!!$       iQ = 1
!!$    else
!!$       iQ = 0
!!$    endif
!!$    outPart(1)%charge = inPart%charge - iQ
!!$!    outPart(2)%charge = iQ

    outPart(1)%charge = 1

    ! select momenta according MC distributions:

    do
       pT2 = rnExp(-3.0)
       phi = twopi*rn()
       z = rn()

       pz = (z*nu)**2-outPart(1)%mass**2-pT2
       if (pz.ge.0.0) exit
    end do
    pz = sqrt(pz)
    pT = sqrt(pT2)

    ! set momenta in nucleon rest frame:

    call MP_Set4(11, 0.0, 0.0, 0.0, sqrt(nu**2+Q2), nu )
    call MP_Set3(12,  mN, 0.0, 0.0, 0.0)
    call MP_Set3(13, outPart(1)%mass,  pT*cos(phi),  pT*sin(phi), pZ)
!    call MP_Set3(14, outPart(2)%mass, -pT*cos(phi), -pT*sin(phi), sqrt(nu**2+Q2)-pZ)
!    call MP_write(6,11,14)

    ! boost to cm frame:

    EE = (Wfree**2+Q2+mN**2)/(2*Wfree)
    pp = sqrt(EE**2-mN**2)

    call MP_ROBO(11,14, 0.0,0.0, 0.0,0.0,-pp/EE)

!    call MP_write(6,11,14)

    ! boost to final system:

    phiB   = atan2(eNev%pcm(2),eNev%pcm(1))
    thetaB = atan2(sqrt(eNev%pcm(1)**2+eNev%pcm(2)**2),eNev%pcm(3))

    call MP_ROBO(11,14, thetaB, phiB, eNev%betacm(1),eNev%betacm(2),eNev%betacm(3) )

!    call MP_write(6,11,14)
!    stop

    ! set momenta into particle vector

    do i=1,3
       outPart(1)%momentum(i) = real(MP_P(13,i))
!       outPart(2)%momentum(i) = real(MP_P(14,i))
    end do

    outPart(1)%momentum(0) = outPart(1)%mass
    outPart(1)%momentum(0) = sqrt(DOT_PRODUCT(outPart(1)%momentum,outPart(1)%momentum))

!    outPart(2)%momentum(0) = outPart(2)%mass
!    outPart(2)%momentum(0) = sqrt(DOT_PRODUCT(outPart(2)%momentum,outPart(2)%momentum))


    flagOK = .TRUE.

  end subroutine DoColl_gammaN_Toy

end module Coll_gammaN
