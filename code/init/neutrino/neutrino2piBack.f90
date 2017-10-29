!******************************************************************************
!****m* /neutrino2piBack
! NAME
! module neutrino2piBack
!
! PURPOSE
! provide some routines for generating a 2 pion background for neutrino
! induced events
!******************************************************************************
module neutrino2piBack
  use particleDefinition
  use eN_eventDefinition

  implicit none



  public :: DoNu2piBack

  real, save:: Norm2pi = 1.0

  real, dimension(0:1, -3:3, 3), save:: RelFak2pi ! (qnuk,eN%idProcess,iProc)

  private

  logical, save:: initflag = .true.

  integer, dimension(0:1, -3:3, 3, 3), save:: Proc2piCharge

contains

  !****************************************************************************
  !****is* neutrino2piBack/readInput
  ! NAME
  ! subroutine readInput
  !****************************************************************************
  subroutine readInput
    use output

    integer :: ios

    NAMELIST /neutrino2piBack/ &
         & Norm2pi

    if (.not.initFlag) return

    ! CC:

    RelFak2pi( 0, 2, 1:3 ) = (/ 1./3., 1./3., 1./3. /)
    Proc2piCharge(0, 2, 1, 1:3) = (/0, 1, 0/)
    Proc2piCharge(0, 2, 2, 1:3) = (/1, 0, 0/)
    Proc2piCharge(0, 2, 3, 1:3) = (/1, 1,-1/)

    RelFak2pi( 1, 2, 1:3 ) = (/ 1./3., 1./3., 0. /)
    Proc2piCharge(1, 2, 1, 1:3) = (/0, 1, 1/)
    Proc2piCharge(1, 2, 2, 1:3) = (/1, 1, 0/)
    Proc2piCharge(1, 2, 3, 1:3) = (/9, 9, 9/) ! does not exist


    ! antiCC:

    RelFak2pi( 0,-2, 1:3 ) = (/ 1./3., 1./3., 0. /)
    Proc2piCharge(0,-2, 1, 1:3) = (/0,-1, 0/)
    Proc2piCharge(0,-2, 2, 1:3) = (/1,-1,-1/)
    Proc2piCharge(0,-2, 3, 1:3) = (/9, 9, 9/) ! does not exist

    RelFak2pi( 1,-2, 1:3 ) = (/ 1./3., 1./3., 1./3. /)
    Proc2piCharge(1,-2, 1, 1:3) = (/0, 0, 0/)
    Proc2piCharge(1,-2, 2, 1:3) = (/0, 1,-1/)
    Proc2piCharge(1,-2, 3, 1:3) = (/1,-1, 0/)



    call Write_ReadingInput('neutrino2piBack',0)
    rewind(5)
    read(5,nml=neutrino2piBack,IOSTAT=ios)
    call Write_ReadingInput("neutrino2piBack",0,ios)

    call Write_ReadingInput('neutrino2piBack',1)

    initFlag=.false.
  end subroutine readInput


  subroutine DoNu2piBack(eN,outPart,XS)

    use constants, only: pi, mN, mPi, alphaQED, GF
    use IdTable, only: pion,nucleon
    use eventGenerator_eN_lowEnergy, only: init_2Pi_getBG
    use eN_event, only: eNeV_GetKinV
    use ParamEP, only: CalcParamEP
    use nBodyPhaseSpace, only: momenta_in_3BodyPS
    use energyCalc, only: energyCorrection
    use mediumDefinition
    use mediumModule, only: mediumAt
    use monteCarlo, only: monteCarloChoose
    use offShellPotential, only: setOffShellParameter

    type(electronNucleon_event), intent(inout) :: eN
    type(particle),dimension(:), intent(inout) :: OutPart
    real, intent(out) :: XS

    real :: nu,Q2,W,Wfree,eps,fT
    integer :: i,qnuk, iProc
    real, dimension(0:3) :: sig2Pi
    real :: dum1,dum2
    real, dimension(1:3,1:3) :: p3       ! momenta of three particles
    logical :: flagOK
    real, dimension(0:3) :: ptot
    real, dimension(1:3) :: betaCMToLab = 0.0
    type(medium)         :: mediumAtPosition
    real :: sumWeight
    real :: Q2fak

    if (initFlag) call readInput

    XS = 0.0
    if (eN%QSquared.ge.5.0) return

    call eNeV_GetKinV(eN, nu,Q2,W,Wfree,eps,fT)

!    write(*,*) 'Wfree, Q2: ', Wfree, Q2


    if (Wfree.gt.1.7) return

    call init_2Pi_getBG(eN%nucleon_free, Wfree, sig2pi)

    if (sig2Pi(0).lt.10e-20) return

    !===== Scale the event acording Q2:

    XS = sig2pi(0)

    call CalcParamEP(Wfree,0.0,0.0, dum1)
    call CalcParamEP(Wfree,Q2, eps, dum2)
!    call CalcParamEP(Wfree,Q2, 0.0, dum2)
    XS = XS * dum2/dum1

    !==========================================================================
    XS = XS*fT/ ( 1e3* pi/(eN%lepton_out%momentum(0)*eN%lepton_in%momentum(0)))

    XS = XS * (GF/(2*pi*alphaQED))**2 ! scaling according DIS cross sections
    XS = XS*Norm2pi

    Q2fak = Q2/(1+Q2/(80.399**2)) ! mW=80.399 GeV

!    XS = XS*Q2 ! DUMMY !!!
    XS = XS*(pi**2/2)*Q2fak**2  ! DUMMY !!!

    !===== 2: Select Process
    qnuk = eN%nucleon_free%charge
    iProc = monteCarloChoose(RelFak2pi(qnuk,eN%idProcess,1:3), sumWeight)
    if (iProc.eq.0) then
       write(*,*) 'neutrino 2piBack.f90:  Failure!'
       stop
    end if
    OutPart(1:3)%charge = Proc2piCharge(qnuk,eN%idProcess,iProc,1:3)
    XS = XS*sumWeight

    !===== 3: generate an event (Q2=0):

    OutPart(1:3)%ID=(/nucleon, pion, pion /)
    OutPart(1)%mass=mN
    OutPart(2:3)%mass=mPi

    OutPart(1:3)%antiparticle=.false.
    OutPart(1:3)%scaleCS=1.

    p3 = momenta_in_3BodyPS (Wfree, OutPart(1:3)%mass)

    do i=1,3
       OutPart(i)%momentum(1:3)=p3(i,:)
       OutPart(i)%momentum(0)=FreeEnergy(OutPart(i))
    end do

    !===== 4: boost the event to the gamma* N system:

    mediumAtPosition=mediumAt(eN%nucleon%position)

    ptot = eN%nucleon%momentum+eN%boson%momentum
    betaCMToLab = ptot(1:3)/ptot(0)

    call energyCorrection(W, (/0.,0.,0./), betaCMToLab, mediumAtPosition,OutPart,flagOK)

    !write(*,*) sqrts(OutPart), W
    if (.not.flagOK) then
       OutPart(1:3)%ID = 0
       return
    end if

    !===== 6: Take care of offshellness
    call setOffShellParameter(OutPart,flagOK)
    if (.not.flagOK) XS = 0.0

    OutPart(1:3)%perWeight=XS


!    write(*,*) ' ----> XS = ',XS


  end subroutine DoNu2piBack




end module neutrino2piBack
