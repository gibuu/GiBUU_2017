!******************************************************************************
!****m* /ParamEP
! NAME
! module ParamEP
!
! PURPOSE
! This module defines routines which return parametrizations of the
! electron proton cross section.
!
! INPUTS
! The Namelist "paramEP" in the Jobcard.
!
! NOTES
! At the moment we include the parametrizations:
! * F. W. Brasse et al.,
!   ``Parametrization Of The Q**2 Dependence Of Virtual Gamma P Total
!   Cross-Sections In The Resonance Region,''
!   Nucl. Phys. B {\bf 110}, 413 (1976).
! * M. E. Christy and P. E. Bosted,
!   ``Empirical Fit to Precision Inclusive Electron-Proton Cross Sections
!   in the Resonance Region,''
!   Phys.Rev. C81 (2010) 055213
!
! In addition, we also include here the the ALLM parametrization for high W
! values (W>1.75 GeV):
! * H. Abramowicz, E. M. Levin, A. Levy and U. Maor,
!   ``A Parametrization of sigma-T (gamma* p) above the resonance region
!   Q**2 >= 0,'' Phys. Lett. B {\bf 269} (1991) 465.
! Here the authors claim, that this provides a smooth continuation of the
! parametrization by Brasse et al.. The newest version is also implemented:
! * H. Abramowicz and A. Levy,
!   ``The ALLM parameterization of sigma(tot)(gamma* p): An update,''
!   arXiv:hep-ph/9712415.
!
! We also provide the parametrizations for R=sigma_L/sigma_T by:
! * L.W.Whitlow et al.,
!   ``A Precise extraction of R = sigma-L / sigma-T from a global analysis
!   of the SLAC deep inelastic e p and e d scattering cross-sections,''
!   Phys.Lett.B250:193-198,1990.
! * V.Tvaskis et al.,
!   ``Longitudinal-transverse separations of structure functions at low Q**2
!   for hydrogen and deuterium,''
!   Phys.Rev.Lett.98:142301,2007.
!   PhD thesis, http://www1.jlab.org/Ul/Publications/documents/thesis_V_Tvaskis.pdf
!******************************************************************************
module ParamEP

  implicit none
  private

  public :: CalcParamEP
  public :: CalcParamEP_ALLM
  public :: CalcParamEP_ALLM97
  public :: CalcParamEP_R1990
  public :: Flux_Bosted

  logical,save :: initFlag=.true.

  !****************************************************************************
  !****g* ParamEP/useParam
  ! SOURCE
  integer, save :: useParam = 2
  ! PURPOSE
  ! select, which parametrization to use:
  ! * 1: Brasse
  ! * 2: Bosted
  !****************************************************************************

  real, save, dimension(3,56,4) :: Bra1
  real, save, dimension(3,56,6) :: Bra2


contains

  !****************************************************************************
  !****s* ParamEP/initInput
  ! NAME
  ! subroutine initInput
  !
  ! PURPOSE
  ! Reads in job card, checks the settings of the input parameters and also
  ! reads the data arrays
  !****************************************************************************
  subroutine initInput
    use output
    use inputGeneral, only: path_to_input

    integer :: ios

    character(100) :: filename
    integer :: i
    !**************************************************************************
    !****n* ParamEP/paramEP
    ! NAME
    ! NAMELIST /paramEP/
    ! PURPOSE
    ! Namelist for module ParamEP includes:
    ! * useParam
    !**************************************************************************
    NAMELIST /paramEP/ useParam

    call Write_ReadingInput('paramEP',0)

    rewind(5)
    read(5,nml=paramEP,IOSTAT=ios)
    call Write_ReadingInput("paramEP",0,ios)

    ! Reading Brasse:

    filename=trim(path_to_input)//'/electronNucleon/Brasse/brasse.dat'
    open(77,file=filename,iostat=ios,status='old')
    call iosCheck(ios,fileName)
    filename=trim(path_to_input)//'/electronNucleon/Brasse/brasse_0906.dat'
    open(78,file=filename,iostat=ios,status='old')
    call iosCheck(ios,fileName)
    filename=trim(path_to_input)//'/electronNucleon/Brasse/brasse_lt06.dat'
    open(79,file=filename,iostat=ios,status='old')
    call iosCheck(ios,fileName)
    filename=trim(path_to_input)//'/electronNucleon/Brasse/brerror.dat'
    open(87,file=filename,iostat=ios,status='old')
    call iosCheck(ios,fileName)
    filename=trim(path_to_input)//'/electronNucleon/Brasse/brerror_0906.dat'
    open(88,file=filename,iostat=ios,status='old')
    call iosCheck(ios,fileName)
    filename=trim(path_to_input)//'/electronNucleon/Brasse/brerror_lt06.dat'
    open(89,file=filename,iostat=ios,status='old')
    call iosCheck(ios,fileName)

    do i=1,56
       read(77,*) Bra1(1,i,:)
       read(78,*) Bra1(2,i,:)
       read(79,*) Bra1(3,i,:)
       read(87,*) Bra2(1,i,:)
       read(88,*) Bra2(2,i,:)
       read(89,*) Bra2(3,i,:)
    end do

    Bra2 = Bra2 * 1e-4

    close(77)
    close(78)
    close(79)
    close(87)
    close(88)
    close(89)

    initFlag = .false.

    call Write_ReadingInput('paramEP',1)

  contains

    subroutine iosCheck(ios,name)
      integer :: ios
      character(*) :: name

      if (ios.ne.0) then
         write(*,'(2A)') 'Error in opening input file: ',name
         write(*,'(A,I5)') 'I/O status=', ios
         write(*,'(A)') 'STOP!'
         stop
      end if

    end subroutine iosCheck

  end subroutine initInput

  !****************************************************************************
  !****s* ParamEP/CalcParamEP
  ! NAME
  ! subroutine CalcParamEP(W,Q2,eps, XS, XSerr)
  !
  ! PURPOSE
  ! Calculate the XS (and its error) according the selected Parametrization
  !
  ! INPUTS
  ! * real                        :: W        -- incoming photon (W)
  ! * real                        :: Q2       -- incoming photon (Q^2)
  ! * real                        :: eps      -- incoming photon (epsilon)
  !
  ! OUTPUT
  ! * real                        :: XS       -- cross section
  ! * real,OPTIONAL               :: XSerr    -- error
  !****************************************************************************
  subroutine CalcParamEP(W,Q2,eps, XS, XSerr)
    use CallStack

    real, intent(in)            :: W, Q2, eps
    real, intent(out)           :: XS
    real, intent(out), optional :: XSerr

    if (initFlag) call initInput

    select case (useParam)
    case (1)
       if (present(XSerr)) then
          call CalcParamEP_Brasse(W,Q2,eps, XS, XSerr)
       else
          call CalcParamEP_Brasse(W,Q2,eps, XS)
       end if
    case (2)
       call CalcParamEP_Bosted(W,Q2,eps, XS)
       if (present(XSerr)) XSerr = 0.0 ! error not possible here
    case default
       call TRACEBACK("wrong useParam")
    end select

  end subroutine CalcParamEP

  !****************************************************************************
  !****s* ParamEP/CalcParamEP_Brasse
  ! NAME
  ! subroutine CalcParamEP_Brasse(W,Q2,eps, XS, XSerr)
  !
  ! PURPOSE
  ! Calculate the XS (and its error) according the Brasse Parametrization:
  ! * F.~W.~Brasse et al.,
  !   ``Parametrization Of The Q**2 Dependence Of Virtual Gamma P Total
  !   Cross-Sections In The Resonance Region,''
  !   Nucl. Phys. B {\bf 110}, 413 (1976)
  !
  ! INPUTS
  ! * real                        :: W        -- incoming photon (W)
  ! * real                        :: Q2       -- incoming photon (Q^2)
  ! * real                        :: eps      -- incoming photon (epsilon)
  !
  ! OUTPUT
  ! * real                        :: XS       -- cross section
  ! * real,OPTIONAL               :: XSerr    -- error
  !
  ! NOTES
  ! The range of validity is:
  ! * eps = 0...0.6...0.9...1.0
  ! * Q2 = 0 ... ??? GeV^2
  ! * W = 1.1 ... 2.0 GeV
  !
  ! The returned cross section is
  !    \sigma^* = \sigma_T+\epsilon\sigma_L
  !             = \frac{1}{\Gamma} \frac{d\sigma}{dE' d\Omega}
  ! Unfortunately, the authors do not giv a definition of Gamma.
  !
  ! An interpolation between different W bins would smoothen the results.
  !****************************************************************************
  subroutine CalcParamEP_Brasse(W,Q2,eps, XS, XSerr)
    use CallStack, only: TRACEBACK
    use constants, only: mN


    real, intent(in)            :: W, Q2, eps
    real, intent(out)           :: XS
    real, intent(out), optional :: XSerr

    integer :: iSet,iW
    real :: dipol,vecq,vecq0,logsig, dsda,dsdb,dsdc
    real, parameter :: d = 3.0

    ! reset output:

    XS = 0.0
    if (present(XSerr)) XSerr = 0.0

    ! select epsilon value:

    iSet = 0
    if (eps.ge.0.9) then
       iSet = 1
    else if (eps.gt.0.6) then
       iSet = 2
    else if (eps.ge.0.0) then
       iSet = 3
    end if

    if (iSet.eq.0) then
       write(*,*) 'eps=',eps
       call TRACEBACK("wrong eps")
    end if

    ! search entry:

    if (W.lt.1.775) then
       iW = nint( (W-1.110)/0.015 )+1
    else
       iW = nint( (W-1.790)/0.020 )+46
    end if

    if (iW.lt.1) return
    if (iW.gt.56) return

    ! calculate XS:

    dipol=1./(1.+Q2/0.71)**2

    vecq = sqrt(((Bra1(iSet,iW,1)**2-mN**2+Q2)/(2*mN))**2+Q2)
    vecq0= (Bra1(iSet,iW,1)**2-mN**2)/(2*mN)
    logsig = Bra1(iSet,iW,2)+Bra1(iSet,iW,3)*log(vecq/vecq0) &
         & + Bra1(iSet,iW,4)*abs(log(vecq/vecq0))**d
    XS = dipol**2*exp(logsig)

    ! calculate XSerr:

    if (present(XSerr)) then
       dsda=XS
       dsdb=XS*log(vecq/vecq0)
       dsdc=XS*abs(log(vecq/vecq0))**d

       XSerr = Bra2(iSet,iW,1)*dsda**2 &
            & + Bra2(iSet,iW,2)*dsdb**2 &
            & + Bra2(iSet,iW,3)*dsdc**2 &
            & - 2*Bra2(iSet,iW,4)*dsda*dsdb &
            & + 2*Bra2(iSet,iW,5)*dsda*dsdc &
            & - 2*Bra2(iSet,iW,6)*dsdb*dsdc
       XSerr = sqrt(abs(XSerr))
    end if

!    write(*,'(1P,4e14.4)') W,Bra1(iSet,iW,1),XS

  end subroutine CalcParamEP_Brasse

  !****************************************************************************
  !****f* ParamEP/Flux_Bosted
  ! NAME
  ! real function Flux_Bosted(W,Q2,eps)
  !
  ! PURPOSE
  ! return the value of equation (4) in Phys.Rev. C81 (2010) 055213
  !
  ! Multiplying this return value with the value given by
  ! CalcParamEP_Bosted() yields
  !     \frac{d\sigma}{dE' d\Omega}
  !
  ! INPUTS
  ! * real                        :: W        -- incoming photon (W)
  ! * real                        :: Q2       -- incoming photon (Q^2)
  ! * real                        :: eps      -- incoming photon (epsilon)
  !
  ! OUTPUT
  ! * function value -- flux in 1/GeV
  !****************************************************************************
  real function Flux_Bosted(W,Q2,eps)
     use constants, only: twopi,alphaQED,mN

     real, intent(in)            :: W, Q2, eps
     real :: Ebeam,nu

     nu = (W**2+Q2-mN**2)/(2*mN)
     Ebeam = (nu+sqrt((nu**2+Q2)*(1+eps)/(1-eps)))/2
     Flux_Bosted = alphaQED*(Ebeam-nu)*(W**2-mN**2)/(twopi**2*Q2*mN*Ebeam*(1-eps))
  end function Flux_Bosted

  !****************************************************************************
  !****s* ParamEP/CalcParamEP_Bosted
  ! NAME
  ! subroutine CalcParamEP_Bosted(W,Q2,eps, XS)
  !
  ! PURPOSE
  ! Calculate the XS according the Bosted Parametrization:
  ! * M. E. Christy and P. E. Bosted,
  !   ``Empirical Fit to Precision Inclusive Electron-Proton Cross Sections
  !   in the Resonance Region,'' Phys.Rev. C81 (2010) 055213
  !
  ! INPUTS
  ! * real                        :: W        -- incoming photon (W)
  ! * real                        :: Q2       -- incoming photon (Q^2)
  ! * real                        :: eps      -- incoming photon (epsilon)
  !
  ! OUTPUT
  ! * real                        :: XS       -- cross section
  !
  ! NOTES
  ! This is a wrapper function around code provided by P.Bosted.
  !
  ! The returned cross section is
  !    \sigma^* = \sigma_T+\epsilon\sigma_L
  !             = \frac{1}{\Gamma} \frac{d\sigma}{dE' d\Omega}
  ! with
  !    \Gamma = \frac{\alpha E' (W^2-M^2)}{(2\pi)^2 Q^2 M E (1-\epsilon)}
  !
  ! The returned values are given in mub.
  !
  ! The range of validity is:
  ! * eps = 0...1.0
  ! * Q2 = 0 ... 10 GeV^2
  ! * W = 1.1 ... 3.0 GeV
  !****************************************************************************
  subroutine CalcParamEP_Bosted(W,Q2,eps, XS)

    real, intent(in)            :: W, Q2, eps
    real, intent(out)           :: XS

    real(8) F1,R,sigt,sigl

    call christy507(dble(W**2),dble(Q2), F1,R,sigt,sigl)
    XS = sigt+eps*sigl

  contains

    SUBROUTINE christy507(W2,Q2,F1,R,sigt,sigl)
      !c   M.E. Christy and P.E. Bosted, ``Empirical Fit to Precision
      !c    Inclusive Electron-Proton Cross Sections in the Resonance Region'',
      !c    (arXiv:0712.3731). To be submitted to Phys. Rev. C.

      real(8) w2,q2,xval1(50),xvall(50),xval(100)
      real(8) mp,mp2,pi,alpha,xb,sigT,sigL,F1,R  !,F2,FL
      integer i !,npts,sf

      mp = .9382727
      mp2 = mp*mp
      pi = 3.141593
      alpha = 1./137.036


      data xval / &
           & 0.12298D+01,0.15304D+01,0.15057D+01,0.16980D+01,0.16650D+01,&
           & 0.14333D+01,0.13573D+00,0.22000D+00,0.82956D-01,0.95782D-01,&
           & 0.10936D+00,0.37944D+00,0.77805D+01,0.42291D+01,0.12598D+01,&
           & 0.21242D+01,0.63351D+01,0.68232D+04,0.33521D+05,0.25686D+01,&
           & 0.60347D+00,0.21240D+02,0.55746D-01,0.24886D+01,0.23305D+01,&
           & -.28789D+00,0.18607D+00,0.63534D-01,0.19790D+01,-.56175D+00,&
           & 0.38964D+00,0.54883D+00,0.22506D-01,0.46213D+03,0.19221D+00,&
           & 0.19141D+01,0.24606D+03,0.67469D-01,0.13501D+01,0.12054D+00,&
           & -.89360D+02,0.20977D+00,0.15715D+01,0.90736D-01,-.38495D-02,&
           & 0.10362D-01,0.19341D+01,0.38000D+00,0.34187D+01,0.14462D+00,&
           & 0.00000D+00,0.00000D+00,0.00000D+00,0.00000D+00,0.00000D+00,&
           & 0.00000D+00,0.00000D+00,0.00000D+00,0.00000D+00,0.00000D+00,&
           & 0.00000D+00,0.00000D+00,0.29414D+02,0.19910D+02,0.22587D+00,&
           & 0.00000D+00,0.00000D+00,0.38565D+04,0.65717D+00,0.00000D+00,&
           & 0.15792D+03,0.97046D+02,0.31042D+00,0.00000D+00,0.42160D+01,&
           & 0.38200D-01,0.12182D+01,0.00000D+00,0.13764D+02,0.31393D+00,&
           & 0.29997D+01,0.00000D+00,0.55124D+01,0.53743D-01,0.13091D+01,&
           & 0.00000D+00,0.86746D+02,0.40864D-04,0.40294D+01,0.31285D+01,&
           & 0.33403D+00,0.49623D+01,0.00000D+00,0.00000D+00,0.11000D+02,&
           & 0.18951D+01,0.51376D+00,0.00000D+00,0.42802D+01,0.00000D+00 /


      do i=1,50
         xval1(i) = xval(i)
         xvalL(i) = xval(50+i)
         if (i.LE.12) xvalL(i) = xval1(i)
      end do
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)


      xb = q2/(w2+q2-mp2)

      call resmod507(1,w2,q2,xval1,sigT)
      call resmod507(2,w2,q2,xvalL,sigL)

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
!       FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = sigL/sigT

    end SUBROUTINE christy507

    SUBROUTINE RESMOD507(sf,w2,q2,xval,sig)

      real(8) W,w2,q2,mp,mp2,xb,sig,xval(50),mass(7),width(7)!,xth(4),mpi2
      real(8) height(7),rescoef(6,4)!,sig_del,sig_21,sig_22,sig_31,sig_32
      real(8) nr_coef(3,4),sigr(7),wdif(2),sig_nr!,sig_4
      real(8) mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      real(8) petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      real(8) eetacmr(7),epicm,epi2cm,eetacm,br(7,3),ang(7)!,spin(7)
      real(8) pgam(7),pwid(7,3),x0(7),q20,h_nr(3)!,dip,mon
      real(8) sig_res,t,xpr(2),m0!,sig_4L,sigtemp,slope
      real(8) sigrsv(7),sig_nrsv
      INTEGER i,j,num,sf!,l
      real(8) sig_mec
      !logical first/.true./
      common/tst2/sigrsv,sig_nrsv,sig_mec


      mp = 0.9382727
      mpi = 0.135
!       mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + 2.*mpi)

      m0 = 0.125
      if (sf.EQ.2) m0 = xval(49)

      if (sf.EQ.1) then
         q20 = 0.05
      else
         q20 = 0.125
      end if


      !CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.0       !!!  P33(1232)
      br(2,1) = 0.45      !!!  S11(1535)
      br(3,1) = 0.65      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.4       !!!  S11(1650)
      br(6,1) = 0.65      !!!  P11(1440) roper
      br(7,1) = 0.50      !!!  F37(1950)

      !CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.45      !!!  S11(1535)
      br(3,3) = 0.0       !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.1       !!!  S11(1650)
      br(6,3) = 0.0       !!!  P11(1440) roper
      br(7,3) = 0.0       !!!  F37(1950)

      !CCCC  2-pion branching ratios  CCCC

      do i=1,7
         br(i,2) = 1.-br(i,1)-br(i,3)
      end do


      !CCCC   Meson angular momentum   CCCC



      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper
      ang(7) = 3.       !!!  F37(1950)

      do i=1,7     !!!  resonance damping parameter  !!!
         x0(i) = 0.215
         !c        x0(i) = xval(50)
      end do

      x0(1) = 0.15
      x0(1) = xval(50)

      do i=1,7
         br(i,2) = 1.-br(i,1)-br(i,3)
      end do


!       dip = 1./(1.+q2/0.71)**2.             !!!  Dipole parameterization  !!!
!       mon = 1./(1.+q2/0.71)**1.

      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+mpi+mpi)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)


      t = log(log((q2+m0)/0.330**2)/log(m0/0.330**2))

      !CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0d0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0d0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0d0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6              !!!  Read in resonance masses     !!!
         num = num + 1
         mass(i) = xval(i)
      end do
      do i=1,6              !!!  Read in resonance widths     !!!
         num = num + 1
         intwidth(i) = xval(num)
         width(i) = intwidth(i)
      end do

      if (sf.EQ.2) then      !!!  Put in 4th resonance region  !!!
         mass(7) = xval(43)
         intwidth(7) = xval(44)
         width(7) = intwidth(7)
      else
         mass(7) = xval(47)
         intwidth(7) = xval(48)
         width(7) = intwidth(7)
      end if

      do i=1,7
         kr(i) = (mass(i)**2-mp2)/2./mp
         kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
         epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
         ppicmr(i) = SQRT(MAX(0d0,(epicmr(i)**2 - mpi**2)))
         epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
         ppi2cmr(i) = SQRT(MAX(0d0,(epi2cmr(i)**2 - (2.*mpi)**2)))
         eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
         petacmr(i) =  SQRT(MAX(0d0,(eetacmr(i)**2 - meta**2)))

         !CCC   Calculate partial widths   CCC

         pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)&
              &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
         !c         !!!  1-pion decay mode


         pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)&
              &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))&
              &         **(ang(i)+2)   !!!  2-pion decay mode

         pwid(i,2) = W/mass(i)*pwid(i,2)


         pwid(i,3) = 0.          !!!  eta decay mode


         if (i.EQ.2.OR.i.EQ.5) then
            pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)&
                 & *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
            !c         !!!  eta decay only for S11's
         end if



         pgam(i) = (kcm/kcmr(i))**2*&
              &  (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

         pgam(i) = intwidth(i)*pgam(i)

         width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      end do

      !CCC    End resonance kinematics and Widths calculations   CCC


      !CCC    Begin resonance Q^2 dependence calculations   CCC


      do i=1,6
         do j=1,4
            num = num + 1
            rescoef(i,j)=xval(num)
         end do

         if (sf.EQ.1) then

            height(i) = rescoef(i,1)*&
                 &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/&
                 &          (1.+q2/0.91)**rescoef(i,4)

         else

            height(i) = rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)&
                 &                             *exp(-1.*rescoef(i,3)*q2)

         end if


         height(i) = height(i)*height(i)

      end do


      if (sf.EQ.2) then      !!!  4th resonance region  !!!

         height(7) = xval(45)*q2/(1.+xval(46)*q2)*exp(-1.*xval(47)*q2)

      else
         height(7) = xval(49)/(1.+q2/0.91)**1.

      end if
      height(7) = height(7)*height(7)



      !CCC    End resonance Q^2 dependence calculations   CCC


      do i=1,3               !!!  Non-Res coefficients  !!!
         do j=1,4
            num = num + 1
            nr_coef(i,j)=xval(num)
         end do
      end do


      !CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7
         sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. &
              &              + (mass(i)*width(i))**2.)
         sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
         if (sf.eq.1) sigrsv(i) = sigr(i)
         sig_res = sig_res + sigr(i)
      end do

      sig_res = sig_res*w


      !CCC    Finish resonances / start non-res background calculation   CCC


      sig_nr = 0.

      if (sf.EQ.1) then

         do i=1,2

            h_nr(i) = nr_coef(i,1)/     &
                 &       (q2+nr_coef(i,2))** &
                 &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
            sig_nr = sig_nr +h_nr(i)*(wdif(1))**(float(2*i+1)/2.)
         end do

         sig_nr = sig_nr*xpr(1)
         sig_nrsv = sig_nr

      else if (sf.EQ.2) then

         do i=1,1
            sig_nr = sig_nr + nr_coef(i,1)* &
                 &      (1.-xpr(i))**(nr_coef(i,3)+nr_coef(i,2)*t) &
                 &               /(1.-xb)/(q2+q20) &
                 & *(q2/(q2+q20))**nr_coef(i,4)*xpr(i)**(xval(41)+xval(42)*t)
         end do

      end if


      sig = sig_res + sig_nr



!1000  format(8f12.5)

      RETURN
    END SUBROUTINE RESMOD507


  end subroutine CalcParamEP_Bosted

  !****************************************************************************
  !****s* ParamEP/CalcParamEP_ALLM
  ! NAME
  ! subroutine CalcParamEP_ALLM(W,Q2, XS)
  !
  ! PURPOSE
  ! Calculate the XS according the ALLM Parametrization:
  ! * H. Abramowicz, E. M. Levin, A. Levy and U. Maor,
  !   ``A Parametrization of sigma-T (gamma* p) above the resonance region
  !   Q**2 >= 0,'' Phys. Lett. B {\bf 269} (1991) 465.
  !
  ! INPUTS
  ! * real                        :: W        -- incoming photon (W)
  ! * real                        :: Q2       -- incoming photon (Q^2)
  !
  ! OUTPUT
  ! * real                        :: XS       -- cross section (in mub)
  !
  ! NOTES
  ! The range of validity is:
  ! * W > 1.75 GeV
  ! * Q2 = 0...2000GeV^2
  !
  ! This routine returns sigma_tot == sigma_L+sigma_T
  !****************************************************************************
  subroutine CalcParamEP_ALLM(W,Q2, XS)

    real, intent(in)            :: W, Q2
    real, intent(out)           :: XS

    real :: x, W2, t, xPom,xReg, CPom,CReg, aPom,aReg, bPom,bReg, F2Pom,F2Reg

    real, parameter :: C = 0.11220d0      ! 4*Pi^2*alpha*(hbar c)^2
    real, parameter :: M2 = 0.879844d0    ! M_proton^2


    real :: M02,MPom2,MReg2,L2,Q02, &
         &     P_aReg(3),P_bReg(3),P_CReg(3), &
         &     P_aPom(3),P_bPom(3),P_CPom(3)

    real :: f1,f2, a1,a2,a3

    f1(t,a1,a2,a3) = a1+a2*t**a3
    f2(t,a1,a2,a3) = a1+(a1-a2)*(1E0/(1E0+t**a3)-1E0)

    data M02,MPom2,MReg2,L2,Q02 / &
!         &     0.30508d0,10.67564d0,0.20623d0,0.06527d0,0.27799d0/
         &     0.30508d0,10.67564d0,0.20623d0,0.06527d0,0.27001d0/

    data P_aReg,P_bReg,P_CReg / &
!         &     0.60408d0, 0.17353d0, 1.61812d0, &
!         &     1.26066d0, 1.83624d0, 0.81141d0, &
!         &     0.67639d0, 0.49027d0, 2.66275d0/
         &     0.60408d0, 0.08144d0, 2.18355d0, &
         &     1.21849d0, 1.82378d0, 0.71033d0, &
         &     0.67639d0, 0.52072d0, 1.90782d0/

    data P_aPom,P_bPom,P_CPom / &
!         &     -0.04503d0, -0.36407d0, 8.17091d0, &
!         &      0.49222d0,  0.52116d0, 3.55115d0, &
!         &      0.26550d0,  0.04856d0, 1.04682d0/
         &     -0.04503d0, -0.41254d0, 7.54998d0, &
         &      0.34131d0,  0.58400d0, 3.25216d0, &
         &      0.26550d0,  0.10525d0, 5.55143d0/

!      W2 = M2+Q2*(1d0/x-1d0)    ! for x as input
    W2 = W**2                 ! for W as input
    x = Q2/(Q2+W2-M2)         ! for W as input
    xPom = 1d0/(1d0+(W2-M2)/(Q2+MPom2))
    xReg = 1d0/(1d0+(W2-M2)/(Q2+MReg2))
    t = log(log((Q2+Q02)/L2)/log(Q02/L2))

    CReg = f1(t,P_CReg(1),P_CReg(2),P_CReg(3))
    aReg = f1(t,P_aReg(1),P_aReg(2),P_aReg(3))
    bReg = f1(t,P_bReg(1),P_bReg(2),P_bReg(3))

    CPom = f2(t,P_CPom(1),P_CPom(2),P_CPom(3))
    aPom = f2(t,P_aPom(1),P_aPom(2),P_aPom(3))
    bPom = f1(t,P_bPom(1),P_bPom(2),P_bPom(3))

    F2Pom = CPom * xPom**aPom * (1-x)**bPom
    F2Reg = CReg * xReg**aReg * (1-x)**bReg

    XS = 1000* C/((Q2+M02)*(1-x)) * &
         &     (1d0+(4d0*M2*Q2)/(Q2+W2-M2)**2)*(F2Pom+F2Reg)

  end subroutine CalcParamEP_ALLM


  !****************************************************************************
  !****s* ParamEP/CalcParamEP_ALLM97
  ! NAME
  ! subroutine CalcParamEP_ALLM97(W,Q2, XS)
  !
  ! PURPOSE
  ! Calculate the XS according the ALLM97 Parametrization:
  ! * H. Abramowicz and A. Levy,
  !   ``The ALLM parameterization of sigma(tot)(gamma* p): An update,''
  !   arXiv:hep-ph/9712415.
  !
  ! INPUTS
  ! * real                        :: W        -- incoming photon (W)
  ! * real                        :: Q2       -- incoming photon (Q^2)
  !
  ! OUTPUT
  ! * real                        :: XS       -- cross section (in mub)
  ! NOTES
  ! The range of validity is:
  ! * W > sqrt(3) GeV
  ! * Q2 = 0...2000GeV^2
  !
  ! This routine returns sigma_tot == sigma_L+sigma_T
  !****************************************************************************
  subroutine CalcParamEP_ALLM97(W,Q2, XS)

    real, intent(in)            :: W, Q2
    real, intent(out)           :: XS

    real :: x, W2, t, xPom,xReg, CPom,CReg, aPom,aReg, bPom,bReg, F2Pom,F2Reg

    real, parameter :: C = 0.11220d0      ! 4*Pi^2*alpha*(hbar c)^2
    real, parameter :: M2 = 0.879844d0    ! M_proton^2


    real :: M02,MPom2,MReg2,L2,Q02, &
         &     P_aReg(3),P_bReg(3),P_CReg(3), &
         &     P_aPom(3),P_bPom(3),P_CPom(3)

    real :: f1,f2, a1,a2,a3

    f1(t,a1,a2,a3) = a1+a2*t**a3
    f2(t,a1,a2,a3) = a1+(a1-a2)*(1E0/(1E0+t**a3)-1E0)

    data M02,MPom2,MReg2,L2,Q02 / &
         &     0.31985d0, 49.457d0, 0.15052d0, 0.06527d0, 0.52544d0/

    data P_aReg,P_bReg,P_CReg / &
         &     0.58400d0, 0.37888d0, 2.60630d0, &
         &     0.01147d0, 3.75820d0, 0.49338d0, &
         &     0.80107d0, 0.97307d0, 0.31985d0/

    data P_aPom,P_bPom,P_CPom / &
         &     -0.08080d0, -0.44812d0, 1.1709d0, &
         &      0.36292d0,  1.89170d0, 1.8439d0, &
         &      0.28067d0,  0.22291d0, 2.1979d0/

!      W2 = M2+Q2*(1d0/x-1d0)    ! for x as input
    W2 = W**2                 ! for W as input
!    x = 1d0/(1d0+(W2-M2)/Q2)  ! for W as input
    x = Q2/(Q2+W2-M2)         ! for W as input
    xPom = 1d0/(1d0+(W2-M2)/(Q2+MPom2))
    xReg = 1d0/(1d0+(W2-M2)/(Q2+MReg2))
    t = log(log((Q2+Q02)/L2)/log(Q02/L2))

    CReg = f1(t,P_CReg(1),P_CReg(2),P_CReg(3))
    aReg = f1(t,P_aReg(1),P_aReg(2),P_aReg(3))
    bReg = f1(t,P_bReg(1),P_bReg(2),P_bReg(3))

    CPom = f2(t,P_CPom(1),P_CPom(2),P_CPom(3))
    aPom = f2(t,P_aPom(1),P_aPom(2),P_aPom(3))
    bPom = f1(t,P_bPom(1),P_bPom(2),P_bPom(3))

    F2Pom = CPom * xPom**aPom * (1-x)**bPom
    F2Reg = CReg * xReg**aReg * (1-x)**bReg

    XS = 1000* C/((Q2+M02)*(1-x)) * &
         &     (1d0+(4d0*M2*Q2)/(Q2+W2-M2)**2)*(F2Pom+F2Reg)

  end subroutine CalcParamEP_ALLM97

  !****************************************************************************
  !****s* ParamEP/CalcParamEP_R1990
  ! NAME
  ! subroutine CalcParamEP_R1990(W,Q2, R)
  !
  ! PURPOSE
  ! Calculate R=sigma_L/sigma_T according:
  ! * L.W.Whitlow et al.,
  !   ``A Precise extraction of R = sigma-L / sigma-T from a global analysis
  !   of the SLAC deep inelastic e p and e d scattering cross-sections,''
  !   Phys.Lett.B250:193-198,1990.
  !
  ! NOTES
  ! We corrected a typo in the constant b1 from the paper.
  !****************************************************************************
  subroutine CalcParamEP_R1990(W,Q2, R)

    real, intent(in)            :: W, Q2
    real, intent(out)           :: R

    real :: x,QQ
    real, parameter :: M2 = 0.879844d0    ! M_proton^2
    real, parameter :: Lam=0.2, b1=0.0635,b2=0.5747,b3=-0.3534

    QQ=max(0.35,Q2)

    x = 1d0/(1d0+(W**2-M2)/QQ)  ! for W as input

    R = b1*(1+12*QQ*0.125**2/((QQ+1)*(0.125**2+x**2)))/log(QQ/Lam**2) &
         & + b2/QQ &
         & + b3/(QQ**2+0.3**2)

  end subroutine CalcParamEP_R1990

end module ParamEP
