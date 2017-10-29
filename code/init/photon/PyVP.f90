!******************************************************************************
!****m* /PyVP
! NAME
! module PyVP
!
! PURPOSE
! [Py]thia and [V]irtual [P]hotons:
! some routines for gamma* p with PYTHIA 6.x
!
! NOTES
! * KG, 12.12.07:
!   The lines concerning fT and eps marked with "!######" should be
!   checked and be consistent with the definitions in initHiLepton.f90
! * This was originally PyVP.F, 18.10.2001...25.4.2005, and is now
!   rewritten in F90
! * during the rewrite in favour of Type(electronNucleon_event) many of
!   its functionality has been moved somewhere or was obsolete
!******************************************************************************
module PyVP

  implicit none
  private

! --- INITIALISATION:---

  public :: InitPythia   ! (DoFrag,Print1,DoDifr)

! --- GET CROSS SECTION:---

  public :: CollectXS_process     ! (N,iXS,XS, file)
  public :: CollectXS_class       ! (XS, file)
!!$  PUBLIC :: CollectNG_class     ! (NG, file)

! --- SCALE VMD COUPLINGS:---

  public :: ScaleVMD            ! (xVMD)

! --- SOME ROUTINES:---

  public :: MarkLepton

  public :: SetPYTHIAthresh, GetPYTHIAthresh

!!$  PUBLIC :: CalcFlux


! --- Former Common Block variables: ---

  real, parameter :: Par1 = 0.00116171 ! alpha/2 pi
  real, parameter :: Par2 = 8.91227    ! 1/(4*pi^2*alpha*(hc)^2) [GeV^-2 mb^-1]

  real,        save :: PYTHIAthresh=2.0  ! if W<this, Pythia is only doing DIS


contains


  !****************************************************************************
  !****s* PyVP/SetPYTHIAthresh
  ! NAME
  ! subroutine SetPYTHIAthresh(thresh)
  ! PURPOSE
  ! Set the value of module variable PYTHIAthresh
  ! NOTES
  ! cf. also eventGenerator_eN_HiEnergy/PYTHIAthresh
  !****************************************************************************
  subroutine SetPYTHIAthresh(thresh)
    real,intent(in) :: thresh
    PYTHIAthresh = thresh
  end subroutine SetPYTHIAthresh


  !****************************************************************************
  !****f* PyVP/getPYTHIAthresh
  ! NAME
  ! real function GetPYTHIAthresh()
  ! PURPOSE
  ! return the value of module variable PYTHIAthresh
  ! NOTES
  ! cf. also eventGenerator_eN_HiEnergy/PYTHIAthresh
  !****************************************************************************
  real function GetPYTHIAthresh()
    GetPYTHIAthresh = PYTHIAthresh
  end function GetPYTHIAthresh


  !****************************************************************************
  !****s* PyVP/InitPythia
  ! NAME
  ! subroutine InitPythia(eNev,DoFrag,Print1,DoDifr)
  !
  ! PURPOSE
  ! STEP (4) of the initialization: Initialize Pythia
  !
  ! Some default parameters and CKIN-parameters are set
  ! and PYINIT will be called.
  !
  ! After calling this function, you are able to call PYEVNT without any
  ! additional initialisation (everything else should have been done already
  ! by calling InitPythiaCalc etc.)
  !
  ! INPUTS
  ! * type(electronNucleon_event):: eNev -- the electron nucleon event
  ! * logical :: DoFrag -- ???
  ! * integer :: Print1 -- ???
  ! * logical :: DoDifr -- ???
  !
  ! NOTES
  ! * BEWARE: the cuts on the PYTHIA-values x and Q2 cannot be exact.
  !****************************************************************************
  subroutine InitPythia(eNev,DoFrag,Print1,DoDifr)
    use eN_eventDefinition
    use eN_event

    use lorentzTrafo
    use rotation

    type(electronNucleon_event),intent(in)   :: eNev
    logical, intent(in) :: DoFrag,DoDifr
    integer, intent(in) :: Print1


    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    integer MSTU,MSTJ
    double precision PARU,PARJ
    SAVE /PYDAT1/

    COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
    integer MDCY,MDME,KFDP
    double precision BRAT
    SAVE /PYDAT3/

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
    integer MSEL,MSELPD,MSUB,KFIN
    double precision CKIN
    SAVE /PYSUBS/

    integer :: KFi
    character*20 Buf
    character*20,save :: cBeam,cTarget ! Pythia init strings
    real, dimension(0:3) :: pB,pT
    real :: phi,theta, x,Q2
    real, dimension(1:3) :: beta

    KFi = 11               ! Beam   is e-
    call PYNAME(KFi,Buf)
    write(cBeam,101) Buf(1:10)
101 FORMAT ('gamma/',A)


    if (eNev%nucleon_free%charge.gt.0) then
       KFi = 2212          ! Target is proton
    else
       KFi = 2112          ! Target is neutron
    end if
    call PYNAME(KFi,Buf)
    write(cTarget,102) Buf(1:10)
102 FORMAT (A)

    ! set some default

    if (eNev%W_free.lt.PYTHIAthresh) then
       MSTP(14) = 26          ! structure of incoming photon: DIS
    else
       MSTP(14) = 30          ! structure of incoming photon: All
    end if

    if (DoDifr) then
       MSEL = 2               ! 1: default, 2: + elastic & diffractive
    else
       MSEL = 1
    end if

    if (DoFrag) then
       MSTP(111) = 1          ! master switch fragmentation/decay
       MSTP(122) = 0          ! switch init and max XS print-out
    else
       MSTP(111) = 0          ! master switch fragmentation/decay
       MSTP(122) = Print1     ! switch init and max XS print-out
    end if

    ! set kinematical cuts

    x = eNeV_Get_LightX(eNev)
    Q2 = eNev%QSquared

    CKIN(61) = x
    CKIN(62) = x * 1.001d0

    CKIN(65) = Q2
    CKIN(66) = Q2 * 1.001d0

    ! set incoming vectors
!!$    call write_electronNucleon_event(eNev, .false.)

    ! please note: we need the electron vectors in the x-z plane,
    ! therefore we have to append a additional rotation with phi
    ! around the z-axis in order to set the y-component to 0 !!!


    pB = eNev%lepton_in%momentum
    pT = eNev%nucleon_free%momentum

    phi = atan2(eNev%pcm(2),eNeV%pcm(1))
    theta = atan2(sqrt(eNev%pcm(1)**2+eNev%pcm(2)**2),eNev%pcm(3))
    beta = eNev%betacm

    call lorentz(beta,pB)
    call lorentz(beta,pT)

    pB(1:3) = rotateZY (theta, phi, pB(1:3))
    pT(1:3) = rotateZY (theta, phi, pT(1:3))

    pB(1:3) = rotateZY (0.0, eNev%phiLepton, pB(1:3))
    pT(1:3) = rotateZY (0.0, eNev%phiLepton, pT(1:3))

    P(1,1:3) = pB(1:3)
    P(2,1:3) = pT(1:3)


!    write(*,'(A,1P,4e13.5)') 'pB (for Pythia)',pB

    ! Initialize Pythia

    call PYINIT('3MOM', cBeam, cTarget, eNev%lepton_in%momentum(0))

  end subroutine InitPythia



  !****************************************************************************
  !****s* PyVP/CollectXS_process
  ! NAME
  ! subroutine CollectXS_process(N,iXS,XS, file)
  !
  ! PURPOSE
  ! Collect the cross sections for single processes
  !
  ! The array XS will be filled with the cross sections of
  ! the processes listed in iXS. N is the maximal number of listable
  ! processes (as input), but also the number of processes as output.
  ! If file!=0 also output is written to file.
  !
  ! (adapted from PYSTAT)
  !
  ! INPUTS
  ! * integer :: N -- ???
  ! * integer :: iXS(*) -- ???
  ! * integer :: file -- if >0, output is written to specific file
  !
  ! OUTPUT
  ! * real :: XS(*) -- ???
  !
  !****************************************************************************
  subroutine CollectXS_process(N,iXS,XS, file)
    integer file, N, iXS(*)
    double precision XS(*)


    COMMON/PYINT1/MINT(400),VINT(400)
    integer MINT
    double precision VINT
    SAVE /PYINT1/

    COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
    integer NGENPD,NGEN
    double precision XSEC
    SAVE /PYINT5/

    COMMON/PYINT6/PROC(0:500)
    CHARACTER PROC*28
    SAVE /PYINT6/

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
    integer MSEL,MSELPD,MSUB,KFIN
    double precision CKIN
    SAVE /PYSUBS/

    integer i,i1

    if (MINT(121).GT.1) CALL PYSAVE(5,0) ! ?????

    do i=0,500
       if (i.eq.0 .or. MSUB(i).eq.1) then
          i1=i1+1
          if (i1.gt.N) stop
          iXS(i1) = i
          XS(i1)  = XSEC(i,3)

          if (file.gt.0) write(file,5200) iXS(i1),PROC(iXS(i1)),XS(i1)

       end if
    end do
    N = i1

5200 FORMAT (I3,1X,A28,': ',1P,D10.3,' mb')

  end subroutine CollectXS_process


  !****************************************************************************
  !****s* PyVP/CollectXS_class
  ! NAME
  ! subroutine CollectXS_class(XS, file)
  !
  ! PURPOSE
  ! Collect the cross sections for classes
  !
  ! The array XS will be filled with the cross sections of
  ! VMD, direct, anomalous, DIS.
  ! If file!=0 also output is written to file.
  !
  ! This routine returns
  !   \sigma^*(y,Q^2) = \sigma_T(y,Q^2) + \epsilon * \sigma_L(y,Q^2)
  ! given by
  !   d\sigma/dy dQ^2 = flux^T \sigma^*(y,Q^2)
  ! with
  !   flux^T = \alpha/(2 \pi) ( (1+(1-y)^2)/y 1/Q^2 - 2 m_e^ y/Q^4 )
  ! as given in [Friberg:2000ra,eq.(2.1)ff]
  !
  ! return values are given im mb.
  !
  ! (adapted from PYSTAT)
  !
  ! INPUTS
  ! * integer :: file -- if >0, output is written to specific file
  !
  ! OUTPUT
  ! * see above
  !
  ! NOTES
  ! The original XS reported by PYTHIA is dsigma/dx dQ2 (x not Bjorken-x!),
  ! since the MC variables are x and Q2. We transform this to dsigma/dy dQ2
  ! by multiplying with the Jacobian. We divide by the flux as used by
  ! PYTHIA. This results in some sigma^(gamma*)(y,Q2), which differs from
  ! the usual definition of sigma_T + eps * sigma_L=dsigma/dE'dOmega / flux
  ! by the Jacobian (dy dQ2)<->(dE'dOmega).
  !
  ! The cross section is calculated with Hand convention.
  !****************************************************************************
  subroutine CollectXS_class(XS, file)
    integer, intent(in) :: file
    real ::  XS(0:4)

    COMMON/PYINT1/MINT(400),VINT(400)
    integer MINT
    double precision VINT
    SAVE /PYINT1/

    COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
    integer NGENPD,NGEN
    double precision XSEC
    SAVE /PYINT5/

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
    integer MSEL,MSELPD,MSUB,KFIN
    double precision CKIN
    SAVE /PYSUBS/

    integer i
    double precision fac, fT_h

    CHARACTER PROGP4(0:4)*28
    DATA PROGP4/&
         &'SUM                         ',&
         &'VMD * hadron                ',&
         &'direct * hadron             ',&
         &'anomalous * hadron          ',&
         &'DIS * hadron                '/

    if (file.gt.0) write(file,*) '=== xSect === (corrected for cuts)'

    fT_h = 2*Par1*VINT(319)/(VINT(320)*VINT(307))
    if (MSTP(16).EQ.0) then
       fT_h = fT_h/VINT(309)
    else
       fT_h = fT_h/VINT(305)
    end if

    !    now fT_h is exactly eq.(2.2) in Friberg:2000ra

    if (MINT(121).eq.4.and.MSTP(14).eq.30) then
       XS(0) = 0d0
       fac = 1d0/( (CKIN(62)-CKIN(61))*(CKIN(66)-CKIN(65)) )
       do i=1,4
          call PYSAVE(3,i)
          XS(i) = XSEC(0,3)*fac/fT_h
          XS(0) = XS(0)+XS(i)
          if (file.gt.0) write(file,5200) PROGP4(i),XS(i),ft_h
       end do
       if (file.gt.0) write(file,5200) PROGP4(0),XS(0)
    else if (MINT(121).eq.1.and.MSTP(14).eq.26) then
       fac = 1d0/( (CKIN(62)-CKIN(61))*(CKIN(66)-CKIN(65)) )
       do i=1,3
          XS(i) = 0d0
       end do
       XS(4) = XSEC(0,3)*fac/fT_h
       XS(0) = XS(4)

       if (file.gt.0) then
          do i=1,4
             write(file,5200) PROGP4(i),XS(i),ft_h
          end do
          write(file,5200) PROGP4(0),XS(0),ft_h
       end if

    else
       write(*,*) 'Ooops! Stop.'
       stop

    end if

5200 FORMAT (A28,': ',1P,D10.3,' mb  ',0P,e13.4)

  end subroutine CollectXS_class

!!$  !***************************************************************************
!!$  !****s* PyVP/CollectNG_class
!!$  ! NAME
!!$  ! subroutine CollectNG_class(NG, file)
!!$  !
!!$  ! PURPOSE
!!$  ! Collect the cross sections for classes
!!$  !
!!$  ! The array NG will be filled with the number of gen. events of
!!$  ! VMD, direct, anomalous, DIS.
!!$  ! If file!=0 also output is written to file.
!!$  !
!!$  ! cf. "CollectXS_class"
!!$  !
!!$  ! (adapted form PYSTAT)
!!$  !
!!$  ! INPUTS
!!$  ! * integer :: file -- if >0, output is written to specific file
!!$  !
!!$  ! OUTPUT
!!$  ! * see above
!!$  !***************************************************************************
!!$  subroutine CollectNG_class(NG, file)
!!$    integer file
!!$    double precision NG(0:4,0:2)
!!$
!!$    COMMON/PYINT1/MINT(400),VINT(400)
!!$    integer MINT
!!$    double precision VINT
!!$    SAVE /PYINT1/
!!$
!!$    COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
!!$    integer NGENPD,NGEN
!!$    double precision XSEC
!!$    SAVE /PYINT5/
!!$
!!$    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
!!$    integer MSTP,MSTI
!!$    double precision PARP,PARI
!!$    SAVE /PYPARS/
!!$
!!$    integer i,j
!!$
!!$    CHARACTER PROGP4(0:4)*28
!!$    DATA PROGP4/&
!!$         &'SUM                         ',&
!!$         &'VMD * hadron                ',&
!!$         &'direct * hadron             ',&
!!$         &'anomalous * hadron          ',&
!!$         &'DIS * hadron                '/
!!$
!!$    if (MINT(121).eq.4.and.MSTP(14).eq.30) then
!!$       NG(0,0) = 0d0
!!$       NG(0,1) = 0d0
!!$       NG(0,2) = 0d0
!!$
!!$       do i=1,4
!!$          call PYSAVE(3,i)
!!$          NG(i,0) = NGEN(0,3) + 1d-20
!!$          NG(0,0) = NG(0,0)+NG(i,0)
!!$       enddo
!!$
!!$       do i=1,4
!!$          call CalcMCError(NG(0,0),NG(i,0),NG(i,1),NG(i,2))
!!$       enddo
!!$
!!$
!!$       if (file.gt.0) then
!!$          write(file,*) '=== NGEN === (with stat. errors)'
!!$          do i=1,4
!!$             write(file,5200) PROGP4(i),(NG(i,j),j=0,2)
!!$          enddo
!!$          write(file,5200) PROGP4(0),(NG(0,j),j=0,2)
!!$       endif
!!$    else
!!$       write (*,*) 'CollectNG_class:Ooops! Stop.'
!!$       stop
!!$    endif
!!$
!!$5200 FORMAT (A28,': ',1P,D10.3,' [+ ',D10.3,',- ',D10.3,']')
!!$
!!$  end subroutine CollectNG_class


  !****************************************************************************
  !****s* PyVP/ScaleVMD
  ! NAME
  ! subroutine ScaleVMD(xVMD)
  !
  ! PURPOSE
  ! Scale the VMD coupling constants:
  ! f_V^2/4pi -> f_V^2/4pi / x
  !
  ! INPUTS
  ! * real :: xVMD(4) -- ???
  ! OUTPUT
  ! * ...
  !****************************************************************************
  subroutine ScaleVMD(xVMD)
    real, intent(in) :: xVMD(4)


    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    double precision PARP0(4)
    save PARP0
    logical Init
    data Init /.TRUE./
    save Init


    integer i

    !      data PARP0 /2.20d0, 23.6d0, 18.4d0, 11.5d0/

    if (Init) then
       do i=1,4
          PARP0(i) = PARP(160+i)
       end do
       !         write(*,*) '...ScaleVMD stored:',PARP0
       Init = .FALSE.
    end if


    do i=1,4
       PARP(160+i) = PARP0(i) / xVMD(i)
    end do
  end subroutine ScaleVMD

  !****************************************************************************
  !****s* PyVP/MarkLepton
  ! NAME
  ! subroutine MarkLepton
  ! PURPOSE
  ! Mark lepton in PYTHIA output as "DOCUMENTATION LINE"
  ! This ensures, that the initial electron will be removed by PYEDIT
  !
  !****************************************************************************
  subroutine MarkLepton

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    integer KF0,i

    KF0 = K(1,2)               ! the incoming lepton

    do i=MSTI(4)+1,N          ! skip documentation
       if (K(i,3).eq.3) then
          if (K(i,2).ne.KF0 .or. K(K(i,3),2).ne.KF0) then
             write(*,*) 'Ooops, what an event. Stop!'
             call PYLIST(1)
             stop
          end if
          K(i,1) = 21
          goto 100
       end if
    end do
100 continue
  end subroutine MarkLepton


!!$c=================================================================
!!$c Save Photon
!!$c saves the four momentum of the Photon in the MomPart-Array Line 4
!!$
!!$      subroutine SavePhoton
!!$
!!$      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
!!$      integer N,NPAD,K
!!$      double precision P,V
!!$      SAVE /PYJETS/
!!$
!!$      if (K(4,2).ne.22) then
!!$         write(*,*) 'Oops, where is the photon?'
!!$         call PYLIST(2)
!!$         stop
!!$      endif
!!$
!!$      call MP_Set4(4, P(4,5), P(4,1), P(4,2), P(4,3), P(4,4))
!!$
!!$      end

!!$c=================================================================
!!$c Calculate the t = (p_i - p_gamma)^2 value
!!$c
!!$c The photon 4-momentum has to be saved in MomPart Array Line 4
!!$
!!$      double precision function CalcT(iLine)
!!$      integer iLine
!!$
!!$      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
!!$      integer N,NPAD,K
!!$      double precision P,V
!!$      SAVE /PYJETS/
!!$
!!$      double precision MP_P
!!$
!!$      CalcT = (P(iLine,4)-MP_P(4,4))**2
!!$     $     - (P(iLine,1)-MP_P(4,1))**2
!!$     $     - (P(iLine,2)-MP_P(4,2))**2
!!$     $     - (P(iLine,3)-MP_P(4,3))**2
!!$
!!$      return
!!$      end

!!$c=================================================================
!!$c Is it an elastic event?
!!$c
!!$c returns true if it is an elastic event (means: gamma N -> V N)
!!$c sets iMeson to 1 (rho), 2 (omega), 3 (phi), 4 (J/Psi)
!!$c sets lMeson to the line of the meson
!!$
!!$      logical function IsItElastic(iMeson,lMeson)

!!$      integer iMeson, lMeson
!!$
!!$      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
!!$      integer N,NPAD,K
!!$      double precision P,V
!!$      SAVE /PYJETS/
!!$
!!$      double precision Par1,Par2
!!$      parameter (
!!$     $     Par1 = 0.00116171d0, ! alpha/2 pi
!!$     $     Par2 = 8.91227d0)    ! 1/(4*pi^2*alpha*(hc)^2) [GeV^-2 mb^-1]
!!$
!!$      COMMON /InitPyt/
!!$     $     s, sqrts,    ! beam-target - cm energy
!!$     $     W2, W,       ! photon-target -cm energy
!!$     $     x,               ! energy fraction gamma/e (<>Bjorken)
!!$     $     Q2,              ! photon virtuality
!!$     $     xB,              ! Bjorken-x
!!$     $     yB,              ! Bjorken-y, light cone fraction
!!$     $     nu,              ! Bjorken, energy of photon
!!$     $     eps,             ! epsilon = flux_L/flux_T
!!$     $     fT,              ! flux_T [GeV^-2]
!!$     $     IP_KF(2),            ! The KF codes
!!$     $     M(2), M2(2), ! masses of beam and target (and squared)
!!$     $     cBeam,cTarget ! Pythia init strings
!!$      integer
!!$     $     IP_KF
!!$      double precision
!!$     $     s,sqrts, W2,W, x, Q2,
!!$     $     xB, yB, nu, eps, fT,
!!$     $     M, M2
!!$      character*20 cBeam,cTarget
!!$      SAVE /InitPyt/
!!$
!!$      iMeson = 0
!!$      lMeson = 0
!!$      IsItElastic = .false.
!!$      if (N.gt.2) then
!!$         return
!!$      endif
!!$      if (K(1,2).eq.IP_KF(2)) then
!!$         lMeson = 2
!!$      elseif (K(2,2).eq.IP_KF(2)) then
!!$         lMeson = 1
!!$      else
!!$         return
!!$      endif
!!$
!!$      if (K(lMeson,2).eq.113) then
!!$         iMeson = 1
!!$      elseif (K(lMeson,2).eq.223) then
!!$         iMeson = 2
!!$      elseif (K(lMeson,2).eq.333) then
!!$         iMeson = 3
!!$      elseif (K(lMeson,2).eq.443) then
!!$         iMeson = 4
!!$      endif
!!$
!!$      IsItElastic = (iMeson.ne.0)
!!$
!!$      return
!!$      end

!!$c=================================================================
!!$c Check the Event, if a particle (forbiddenly) decayed
!!$c
!!$c necessary, because PYTHIA sometimes dacays rho0 (and others?), while
!!$c all particle decay has been explicitly forbidden
!!$
!!$      subroutine CheckEventDecay
!!$
!!$      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
!!$      integer N,NPAD,K
!!$      double precision P,V
!!$      SAVE /PYJETS/
!!$
!!$      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
!!$      integer MSTP,MSTI
!!$      double precision PARP,PARI
!!$      SAVE /PYPARS/
!!$
!!$      integer i
!!$
!!$      do i=MSTI(4)+1,N          ! skip documentation
!!$         if (K(i,1).eq.11) then
!!$            if (K(i,2).lt.10) goto 100 ! quark
!!$            if (K(i,2).eq.21) goto 100 ! gluon
!!$            if (K(i,2).eq.91) goto 100 ! cluster
!!$            if (K(i,2).eq.92) goto 100 ! string
!!$            if (K(i,2).eq.1103) goto 100 ! dd_1
!!$            if (K(i,2).eq.2101) goto 100 ! ud_0
!!$            if (K(i,2).eq.2103) goto 100 ! ud_1
!!$            if (K(i,2).eq.2203) goto 100 ! uu_1
!!$
!!$            write(*,*) 'Ooops, a decayed particle in line ',i,'!!!'
!!$            call PYLIST(2)
!!$            stop
!!$ 100        continue
!!$         endif
!!$      enddo
!!$      end

!!$c=================================================================
!!$c The Conversion-Factor F_2 = F2Fac * sigma [FriSjo00]
!!$c
!!$
!!$      double precision function F2Fac()
!!$
!!$      double precision Par1,Par2
!!$      parameter (
!!$     $     Par1 = 0.00116171d0, ! alpha/2 pi
!!$     $     Par2 = 8.91227d0)    ! 1/(4*pi^2*alpha*(hc)^2) [GeV^-2 mb^-1]
!!$
!!$      COMMON /InitPyt/
!!$     $     s, sqrts,    ! beam-target - cm energy
!!$     $     W2, W,       ! photon-target -cm energy
!!$     $     x,               ! energy fraction gamma/e (<>Bjorken)
!!$     $     Q2,              ! photon virtuality
!!$     $     xB,              ! Bjorken-x
!!$     $     yB,              ! Bjorken-y, light cone fraction
!!$     $     nu,              ! Bjorken, energy of photon
!!$     $     eps,             ! epsilon = flux_L/flux_T
!!$     $     fT,              ! flux_T [GeV^-2]
!!$     $     IP_KF(2),            ! The KF codes
!!$     $     M(2), M2(2), ! masses of beam and target (and squared)
!!$     $     cBeam,cTarget ! Pythia init strings
!!$      integer
!!$     $     IP_KF
!!$      double precision
!!$     $     s,sqrts, W2,W, x, Q2,
!!$     $     xB, yB, nu, eps, fT,
!!$     $     M, M2
!!$      character*20 cBeam,cTarget
!!$      SAVE /InitPyt/
!!$
!!$      F2Fac = Par2 * Q2**2 * (1d0-xB)/(Q2+4d0*M2(2)*xB**2)
!!$      return
!!$      end

!!$c=================================================================
!!$c Calculate the Monte Carlo Error
!!$c
!!$c Given the total Number of events Ntot and the number of interesting
!!$c events N, this sets the absolute (negative and positive) errors
!!$c dN_p and dN_m.
!!$c
!!$c cf. Year 01, eq.{47.4)ff
!!$c
!!$      subroutine CalcMCError(Ntot, N, dN_p, dN_m)
!!$      double precision Ntot, N, dN_p, dN_m
!!$      double precision h0,h1,h2
!!$
!!$      h0 = Ntot/(Ntot+1)
!!$      h1 = N+0.5d0
!!$      h2 = sqrt(0.25d0+N-N**2/Ntot)
!!$
!!$      dN_p = h0*(h1+h2) - N
!!$      dN_m = N - min(N, h0*(h1-h2))
!!$
!!$      end
!!$

!!$  !*************************************************************************
!!$  !****s* PyVP/CalcFlux
!!$  ! NAME
!!$  ! subroutine CalcFlux(Ebeam,Q2,y,x,weight,eps)
!!$  !
!!$  ! PURPOSE
!!$  ! calculate for given kinematics the flux; used by the MC routine which
!!$  ! selects the photon parameters.
!!$  !
!!$  ! INPUTS
!!$  ! * real :: Ebeam -- Energy of lepton beam
!!$  ! * real :: Q2    -- virtuality
!!$  ! * real :: y     -- energy fraction
!!$  ! * real :: x     -- Bjorken-x
!!$  ! OUTPUT
!!$  ! * real :: weight -- transversal flux fT in GeV^-3
!!$  ! * real :: eps    -- epsilon = fL/fT
!!$  !
!!$  ! NOTES
!!$  ! * the original formulae were implemented by Thomas Falter (PhD,
!!$  !   eq.(3.15),(3.16)). We modified and corrected them.
!!$  ! * we are now trying to supersede those by a "low-W/low-Q2" expression.
!!$  ! * The equivalent photon energies are (we use the Hand convention)
!!$  !     K = (W^2-M^2)/2M = (1-x)*nu = (1-x)*y*Ebeam  [Hand]
!!$  !     K = sqrt[nu^2+Q^2] [Gilman]
!!$  ! * The MC random variables are nu and Q2. Unfortunatly, in older code versions,
!!$  !   the returned flux corresponds to the definition of d\sigma/dy dQ^2 (y=nu/E)
!!$  !   and is therefore only correct up to a factor E^2. ?????????????
!!$  !
!!$  !*************************************************************************
!!$  subroutine CalcFlux(Ebeam,Q2,y,x,weight,eps)
!!$
!!$    use constants
!!$
!!$    real,intent(in)  :: Ebeam,Q2,y,x
!!$    real,intent(out) :: weight,eps
!!$
!!$    real :: K,nu,cL,cT
!!$
!!$    nu = Ebeam*y
!!$
!!$    ! this was implemented by Thomas Falter; PhD, eq.(3.15),(3.16)
!!$    !-------------------------------------------------------------
!!$! !$    K = (1.-x)*nu ! Hand convention
!!$! !$    cL = 2.*(1.-y)
!!$! !$    cT = cL + y**2
!!$! !$    weight = alphaQED/twopi * K/(Q2*nu**2) * cT
!!$! !$    eps = cL/cT
!!$
!!$    ! this is the low-W expansion of the above ones: (m_e=0)
!!$    !-----------------------------------------------
!!$    ! (cf. Falter PhD [appendix D], Tytgat PhD)
!!$
!!$    K = (1.-x)*nu ! Hand convention
!!$    cL = 2.*(1.-y-Q2/(4*Ebeam**2))/(1.+Q2/nu**2)
!!$    cT = cL + y**2
!!$    weight = alphaQED/twopi * K/(Q2*nu**2) * cT
!!$    eps = cL/cT
!!$
!!$
!!$  end subroutine CalcFlux


end module PyVP
