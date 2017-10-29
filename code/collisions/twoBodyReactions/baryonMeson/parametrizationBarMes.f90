!******************************************************************************
!****m* /parametrizationsBarMes
! NAME
! module parametrizationsBarMes
!
! PURPOSE
! Includes routines which are parametrizations of "baryon meson -> X" data.
!******************************************************************************
module parametrizationsBarMes

  use cl_splines
  implicit none

  private

  type(tspline),save :: s1,s2,s3,s4

  public :: golub_omega, golub_phi
  public :: omegaN_lykasov
  public :: piN_elastic, piN_chargeExchange
  public :: sigelkp, sigCEXkp, sigCEXk0
  public :: kppbg, kpnbg, kaonbg
  public :: huang, huangd, huanglamd, huanglam
  public :: sibirtpi
  public :: JPsiN
  public :: pin_to_strangebaryon_kaon_pion
  public :: pin_to_strangebaryon_kaon_pion_matrix, sigma_KbarToXi
  public :: cleanup


contains


  subroutine cleanup
    call cl_cleanupSpline(s1)
    call cl_cleanupSpline(s2)
    call cl_cleanupSpline(s3)
    call cl_cleanupSpline(s4)
  end subroutine


  !****************************************************************************
  !****f* parametrizationsBarMes/piN_elastic
  ! NAME
  ! real function piN_elastic(charge,srts,success)
  !
  ! PURPOSE
  ! Returns the cross section for elastic pion proton scattering.
  !
  ! These are own brute force fits to data for plab=0..20 GeV, i.e. for sqrt(s)<6.2 GeV !!!
  !
  ! INPUTS
  ! * integer ::  charge   -- pion charge
  ! * integer ::  srts     -- SQRT(s) in GeV
  !
  ! OUTPUT
  ! * logical :: success -- .true. if result is useful
  ! * funcion value :: sigma(pi N -> pi N) in mb
  !
  ! NOTES
  ! * Fitted to data in the range of plab=0..20 GeV (only pi- and pi+ !!).
  ! * Brute force fits to data
  !****************************************************************************
  real function piN_elastic(charge,srts,success)
    use twoBodyTools, only: p_lab
    use constants, only: mN,mPi

    integer,intent(in) :: charge
    real  :: plab
    real,intent(in) :: srts
    logical, intent(out) :: success

    success=.true.

    if (srts.le.mPi+mN) then
       piN_elastic=0.
       return
    end if
    plab=p_lab(srts,mPi,mN)

    select case (charge)
    case (-1)
       if (plab.lt.0.53) then
          piN_elastic=piMinus_1(plab)
       else if (plab.lt.2.0) then
          piN_elastic=piMinus_2(plab)
       else if (plab.lt.20.) then
          piN_elastic=piMinus_3(plab)
       else
          piN_elastic=piMinus_3(20.)
          success=.false.
       end if
    case (1)
       if (plab.lt.0.35) then
          piN_elastic=piPlus_1(plab)
       else if (plab.lt.0.95) then
          piN_elastic=piPlus_2(plab)
       else if (plab.lt.2.0) then
          piN_elastic=piPlus_2a(plab)
       else if (plab.lt.20.) then
          piN_elastic=piPlus_3(plab)
       else
          piN_elastic=piPlus_3(20.)
          success=.false.
       end if
    case default
       piN_elastic=0.
       success=.false.
    end select
    piN_elastic=max(0.,piN_elastic)

  contains
    real function piPlus_1(x)
      real,intent(in) :: x
      real,parameter ::w     = 0.0377622 !     +/- 0.3358       (889.3%)
      real,parameter ::h     = 0.209973  !     +/- 2.861        (1362%)
      real,parameter ::a     = -20.2246  !     +/- 15.58        (77.05%)
      real,parameter ::b     = 216.724   !     +/- 811.5        (374.5%)
      piPlus_1=h/((x-0.301)**2+w**2)+a+b*x
    end function piPlus_1

    real function piPlus_2(x)
      real,intent(in) :: x
      real,parameter :: a0     = 1387.3    !       +/- 159.1        (11.47%)
      real,parameter :: a1     = -8045.61  !       +/- 1048         (13.02%)
      real,parameter :: a2     = 18300.7   !       +/- 2576         (14.08%)
      real,parameter :: a3     = -19012.1  !       +/- 2810         (14.78%)
      real,parameter :: a4     = 7553.98   !       +/- 1151         (15.23%)
      real,parameter :: a5     = -5460.35  !       +/- 889.1        (16.28%)
      !      piPlus_2=a0+(a1-0.5)*x+(a2-0.5)*x**2+(a3-0.5)*x**3+(a4-0.5)*x**4 +(a5-0.5)*(x-0.5)**5

      ! since the fit function is a little bit strange (c.f. (x-0.5)**5)
      ! the coefficients in the Horner scheme are a little bit compilcated
      real,parameter :: b0 = 0.015625 + a0 - 0.03125 * a5
      real,parameter :: b1 = -0.65625 + a1 + 0.3125 * a5
      real,parameter :: b2 = 0.125 + a2 - 1.25 * a5
      real,parameter :: b3 = -1.75 + a3 + 2.5 * a5
      real,parameter :: b4 = 0.75 + a4 - 2.5 * a5
      real,parameter :: b5 = -0.5 + a5
      piPlus_2=b0+x*(b1+x*(b2+x*(b3+x*(b4+x*(b5)))))

    end function piPlus_2

    real function piPlus_2a(x)
      real,intent(in) :: x
      real,parameter :: s0     = 141.954  !        +/- 151          (106.4%)
      real,parameter :: s1     = -492.093 !        +/- 418.4        (85.02%)
      real,parameter :: s2     = 629.658  !        +/- 426.9        (67.8%)
      real,parameter :: s3     = -326.288 !        +/- 190.5        (58.39%)
      real,parameter :: s4     = 59.9584  !        +/- 31.43        (52.42%)
      piPlus_2a=s0+(s1-0.5)*x+(s2-0.5)*x**2+(s3-0.5)*x**3+(s4-0.5)*x**4
    end function piPlus_2a


    real function piPlus_3(x)
      real, intent(in) :: x
      real,parameter :: c0      = 2.98991   !       +/- 0.0513       (1.716%)
      real,parameter :: c1      = 14.565
      piPlus_3=c0+c1/x
    end function piPlus_3

    real function piMinus_1(x)
      real,intent(in) :: x
      real,parameter :: w       = 0.0656649
      real,parameter :: h       = 0.0912381
      real,parameter :: a       = 47.8669
      real,parameter :: alpha   = 2.4718
      piMinus_1=h/((x-0.301)**2+w**2)+a*x**alpha
    end function piMinus_1

    real function piMinus_2(x)
      real,intent(in) :: x
      real,parameter ::w    = 0.0614504 !       +/- 0.01331      (21.66%)
      real,parameter ::h    = 0.0317369 !       +/- 0.01659      (52.29%)
      real,parameter ::w2   = 0.0802949 !       +/- 0.01178      (14.68%)
      real,parameter ::h2   = 0.0901384 !       +/- 0.03118      (34.59%)
      real,parameter ::w3   = 4.09971   !       +/- 370.4        (9034%)
      real,parameter ::h3   = 9.37678   !       +/- 2237         (2.386e+04%)
      real,parameter ::w4   = 1.24181   !       +/- 3.483        (280.5%)
      real,parameter ::h4   = 18.2444   !       +/- 177.8        (974.4%)
      piMinus_2=h/((x-0.72)**2+w**2)+h2/((x-1.026)**2+w2**2)+h3/((x-0.869)**2+w3**2)+h4/((x-1.2)**2+w4**2)
    end function piMinus_2

    real function piMinus_3(x)
      real, intent(in) :: x
      real,parameter ::  c0              = 3.2356
      real,parameter ::  c1              = 11.7764
      piMinus_3=c0+c1/x
    end function piMinus_3
  end function piN_elastic


  !****************************************************************************
  !****f* parametrizationsBarMes/piN_chargeExchange
  ! NAME
  ! real function piN_chargeExchange(srts)
  !
  ! PURPOSE
  ! Returns the cross section for pion nucleon charge exchange, i.e.:
  ! * p pi- <-> n pi0
  ! * n pi+ <-> p pi0
  !
  ! Fitted to data in the range of plab=0..1.4 GeV (i.e. sqrt(s)<1.9 GeV)
  !
  ! INPUTS
  ! * integer ::  srts     -- sqrt(s) in GeV
  !
  ! OUTPUT
  ! * function value       -- charge exchange Xsection in mb
  !
  ! NOTES
  ! See buuinput/pionNucleon/chargeExchange/fit.nb
  !****************************************************************************
  real function piN_chargeExchange(srts)
    use twoBodyTools, only: p_lab
    use constants, only: mN, mPi

    real, intent(in) :: srts
    real :: x ! =momentum of pion in restframe of nucleon
    logical,save :: initTable=.true.
    real, parameter :: binSize=0.003
    real, parameter :: min=0.8
    real, parameter :: maximum=3
    integer, parameter :: maxIndex = NINT((maximum-min)/binSize)
    real, dimension(0:maxIndex),save :: table
    integer :: i
    logical, parameter :: useTable=.false.
    real :: srts_Dummy

    if (useTable) then
       if (initTable) then
          do i=0,maxIndex
             srts_dummy=mPi+mN+float(i)*binsize
             x=p_lab(srts_dummy,mPi,mN)
             table(i)= max(0.,0.17774869654285305 /(0.03299564550037992   + (-0.5206120292039138 + x)**2) &
                            + 0.07425990408062914 /(0.005773792650686387  + (-0.3391240794763685 + x)**2) &
                            +  0.11108734870135036/(0.0029855031054111065 + (-0.2712812392094471 + x)**2) &
                            + 8.48700351456715*x - 4.084626239576959*x**2)
          end do
          initTable=.false.
       end if

       i=NINT((srts-mPi-mN)/binsize)
       if (i<0) then
!          write(*,*) 'Warning: srts too low in piN_chargeExchange'
          piN_chargeExchange=0.
       else if (i<=maxIndex) then
             piN_chargeExchange=table(i)
       else
          x=p_lab(srts,mPi,mN)
          piN_chargeExchange= 0.17774869654285305 /(0.03299564550037992   + (-0.5206120292039138 + x)**2) &
                            + 0.07425990408062914 /(0.005773792650686387  + (-0.3391240794763685 + x)**2) &
                            +  0.11108734870135036/(0.0029855031054111065 + (-0.2712812392094471 + x)**2) &
                            + 8.48700351456715*x - 4.084626239576959*x**2
          piN_chargeExchange=max(piN_chargeExchange,0.)
       end if
    else
       x=p_lab(srts,mPi,mN)
       piN_chargeExchange= 0.17774869654285305 /(0.03299564550037992   + (-0.5206120292039138 + x)**2) &
                         + 0.07425990408062914 /(0.005773792650686387  + (-0.3391240794763685 + x)**2) &
                         +  0.11108734870135036/(0.0029855031054111065 + (-0.2712812392094471 + x)**2) &
                         + 8.48700351456715*x - 4.084626239576959*x**2
       piN_chargeExchange=max(piN_chargeExchange,0.)
    end if
  end function piN_chargeExchange


  !****************************************************************************
  !****s* parametrizationsBarMes/pionquasibg
  ! NAME
  ! subroutine pionquasibg(nukcharge,pioncharge,sqrts,elastic,chargeEx)
  !
  ! PURPOSE
  ! Subroutine evaluates quasielastic background for low energy pi N->pi N
  !
  ! INPUTS
  ! * integer,intent(in) ::  nukcharge,pioncharge --- charges
  ! * real,intent(in) ::sqrts                --- SQRT(s)
  ! OUTPUT
  ! * real,intent(out) ::  elastic, chargeEx --background crosssection in mB,
  !   elastic and charge exchange
  !****************************************************************************
!   subroutine pionquasibg(nukcharge,pioncharge,sqrtsOriginal,elastic,chargeEx)
!
!     integer,intent(in) :: nukcharge,pioncharge
!     real,intent(in) :: sqrtsOriginal
!     real,intent(out) ::  elastic, chargeEx
!
!     !local
!     real, parameter :: limit=1.15712
!     real, parameter :: nullpunkt=1.18
!     real  ::    dummy
!     real ::     sqrts
!
!     elastic=0.
!     chargeEx=0.
!
!     sqrts=sqrtsOriginal
!
!     If (sqrts.le.nullpunkt) then
!        ! Evaluate Elastic = Kanal C
!        If (pioncharge.eq.0) then     !Kanal D
!           elastic=0.
!        else if (pioncharge+nukcharge.eq.0 .or.pioncharge+nukcharge.eq.1) then !Kanal C
!           If (sqrts.le.limit) then
!              elastic=kanalc(sqrts)
!           else
!              dummy=sqrts-limit
!              sqrts=limit
!              elastic=kanalc(sqrts)*(1-dummy/(nullpunkt-limit))
!           end if
!        else if (pioncharge+nukcharge.eq.2.or.pioncharge+nukcharge.eq.-1) then !Kanal A=Isospin=3/2 channel
!           elastic=0.
!        end if
!
!        sqrts=sqrtsOriginal
!
!
!        ! Evaluate Charge Exchange = Kanal B
!        If((pioncharge+nukCharge.ne.2).and.(pioncharge+nukCharge.ne.-1)) then
!           !          If (sqrts.le.limit) then
!           chargeEx=kanalb(sqrts)
!           !          else
!           !             dummy=sqrts-limit
!           !             sqrts=limit
!           !             chargeEx=kanalb(sqrts)*(1-dummy/(nullpunkt-limit))
!           !          end if
!           if(sqrts.gt.1.4) chargeEx=0.
!        end if
!     end if
!
!     If (elastic.le.0) then !No negative crosssections!!!
!        elastic=0.
!     end if
!     If (chargeEx.le.0) then !No negative crosssections!!!
!        chargeEx=0.
!     end if
!
!
!     return
!
!   contains
!
!     function kanala(sqrts) RESULT (Ergebnis)
!       ! Elastic background for N pi -> N pi
!       ! Corresponds to C_4 at page 53 of Oliver's diploma thesis
!       ! Set to zero when I realized, that there are data for the Isospin=3/2 channel and that we don't need any background there.
!       ! This fit was based on the results of kanalc and kanalb which are exact. So I estimated the C_4  Xsection by a non-coherence assumption
!       ! of I=1/2 and I=3/2. This failed!!!
!       real sqrts
!       real, parameter :: a=1.0292E6
!       real, parameter :: b=-1.07355
!       real, parameter :: c=-0.966852
!       real, parameter :: d=1.39615
!       real, parameter :: e=-2.36214
!       real  Ergebnis
!       !Ergebnis=a*(b+sqrts)*(c+sqrts)*(d+e*sqrts+sqrts**2)
!       Ergebnis=0.
!     end function kanala
!
!     function kanalb(sqrts)  RESULT (Ergebnis)
!       ! Charge exchange cross section background for N pi -> N pi
!       ! Corresponds to C_2 at page 53 of Oliver's diploma thesis
!       real sqrts
!       real, parameter :: a=234585.0
!       real, parameter :: b=-1.06751
!       real, parameter :: c=-1.0073
!       real, parameter :: d=1.39469
!       real, parameter :: e=-2.36109
!       real  Ergebnis
!       Ergebnis=a*(b+sqrts)*(c+sqrts)*(d+e*sqrts+sqrts**2)
!     end function kanalb
!
!     function kanalc(sqrts)  RESULT (Ergebnis)
!       ! Elastic background for N pi -> N pi
!       ! Corresponds to C_3 at page 53 of Oliver's diploma thesis
!       real sqrts
!       real, parameter :: a=126088.0
!       real, parameter :: b=-1.19405
!       real, parameter :: c=-1.15482
!       real, parameter :: d=1.18733
!       real, parameter :: e=-2.17777
!       real  Ergebnis
!       Ergebnis=a*(b+sqrts)*(c+sqrts)*(d+e*sqrts+sqrts**2)
!     end function kanalc
!
!
!     function kanald(sqrts)  RESULT (Ergebnis)
!       ! Elastic background for N pi -> N pi
!       ! Corresponds to C_1 at page 53 of Oliver's diploma thesis
!       ! set to zero when I realized, that there are data for the Isospin=3/2 channel and that we don't need any background there
!       ! This fit was based on the results of kanalc and kanalb which are exact. So I estimated the C_1
!       ! Xsection by a non-coherence assumption of I=1/2 and I=3/2. This failed!!!
!       real sqrts
!       real, parameter :: a=271316.0
!       real, parameter :: b=-1.07289
!       real, parameter :: c=-0.902965
!       real, parameter :: d=1.40657
!       real, parameter :: e=-2.37095
!       real  Ergebnis
!       !Ergebnis=a*(b+sqrts)*(c+sqrts)*(d+e*sqrts+sqrts**2)
!       Ergebnis=0.
!     end function kanald
!
!   end subroutine pionquasibg

  !****************************************************************************
  !****f* parametrizationsBarMes/sigelkp
  ! NAME
  ! real function sigelkp(p,n)
  ! PURPOSE
  ! Calculate:
  ! * K+ p -> K+ p for n=1
  ! * K+ n -> K+ n for n=2
  ! INPUTS
  ! * real :: p -- kaon momentum in the rest frame of the proton
  ! * integer :: n -- tells which cross section is calculated (see above)
  !****************************************************************************
  real function sigelkp(p,n)
    integer :: i, n
    real :: a(0:4),p

    data (a(i),i=0,4) /10.508,-3.716,1.845,-0.764,0.508/

    select case (n)
    case (1)
       ! Direct fit to the data on K+ p -> K+ p:
       sigelkp=(a(0)+a(1)*p+a(2)*p**2)/(1.+a(3)*p+a(4)*p**2)
    case (2)
       ! Obtained by subtracting charge-exchange and inelastic contributions
       ! from the total K+ n cross section. The data on sig_tot and sig_inel
       ! are from B.R. Martin, NPB 94, 413 (1975)
       ! (see also the review of C.B. Dover and G.E. Walker, Phys. Rep. 89, 1 (1982)).
       ! sig_cex is parameterized in function sigCEXkp of this module.
       if (p.ge.0. .and. p.le.0.8) then
          sigelkp=5.+3.75*p
       else if (p.le.1.13) then
          sigelkp=6.4+2.*p
       else
          sigelkp=261.2/(exp((p-0.217)/0.220)+1.)+4.53
       end if
     case default
       write(*,*) ' in sigelkp: wrong n= ', n
       stop
     end select

  end function sigelkp

  !****************************************************************************
  !****f* parametrizationsBarMes/sigCEXkp
  ! NAME
  ! real function sigCEXkp(srts)
  ! PURPOSE
  ! Calculate
  ! * K+ n -> K0 p
  ! INPUTS
  ! * real :: srts -- invariant energy (GeV)
  !****************************************************************************
  real function sigCEXkp(srts)

    real, intent(in) :: srts
    real :: s

    s=srts**2

    ! Parameterization of the HEP data:
    if (s > 2.062097) then
       sigCEXkp=331.399*(s/2.062096-1.)**1.77226*(2.062096/s)**6.68116
    else
       sigCEXkp=0.
    end if

  end function sigCEXkp

  !****************************************************************************
  !****f* parametrizationsBarMes/sigCEXk0
  ! NAME
  ! real function sigCEXk0(srts)
  ! PURPOSE
  ! Calculate
  ! * K0 p -> K+ n
  ! INPUTS
  ! * real :: srts -- invariant energy (GeV)
  !****************************************************************************
  real function sigCEXk0(srts)

    real, intent(in) :: srts
    real, parameter :: mKp=0.4937, mK0=0.4977, mp=0.9383, mn=0.9396
    real, parameter :: mKp2=mKp**2, mK02=mK0**2, mp2=mp**2, mn2=mn**2
    real :: s, qf2, qi2

    s=srts**2

    qi2=(s+mK02-mp2)**2/(4.*s)-mK02
    qf2=max(0.,(s+mKp2-mn2)**2/(4.*s)-mKp2)

    if (qi2 > 1.e-06) then
       sigCEXk0=sigCEXkp(srts)*qf2/qi2
    else
       sigCEXk0=sigCEXkp(srts)
    end if

  end function sigCEXk0



  !****************************************************************************
  !****f* parametrizationsBarMes/kppbg
  ! NAME
  ! real function kppbg(plab2)
  ! PURPOSE
  ! Calculate
  ! * K+ p -> p K pion
  ! INPUTS
  ! * real :: plab2
  !****************************************************************************
  real function kppbg(plab2)
    use output

    logical, save :: initFlag=.true.
    integer n
    parameter (n=17)
    real,save ::  plab(n),sig(n)!,sig2(n)
    real plab2,sigma
    logical :: successFlag
    integer :: errorType


    if (initFlag) then
       call Write_InitStatus("kppBg",0)
       call kppbgini
       initFlag=.false.
       call Write_InitStatus("kppBg",1)
    end if

    if (plab2.lt.plab(1)) then
       sigma=0.
    else if (plab2.gt.plab(n)) then
       sigma=sig(n)
    else
       sigma=cl_spline(s1,plab2,successFlag,errorType)
       if (.not.successFlag) call cl_error(errorType,'kppbg',plab2)
    end if
    kppbg=sigma
    return

  contains

    subroutine kppbgini
      use inputGeneral, only: path_to_input
      integer i
      open(13,file=trim(path_to_input)//'/kpp_tot_bg_spl.dat',status='unknown')
      do i=1,n
         read(13,*)plab(i),sig(i)
      end do
      close(13)
      s1=cl_initSpline(plab,sig)
    end subroutine kppbgini

  end function kppbg


  !****************************************************************************
  !****f* parametrizationsBarMes/kpnbg
  ! NAME
  ! real function kpnbg(plab2)
  ! PURPOSE
  ! Calculate:
  ! * K0 p -> n K pion
  ! INPUTS
  ! * real :: plab2
  !****************************************************************************
  real function kpnbg(plab2)
    use output

    integer n
    parameter (n=24)
    real,save ::  plab(n),sig(n)!,sig2(n)
    real plab2,sigma
    logical :: initFlag=.true.
    logical :: successFlag
    integer :: errorType

    if (initFlag) then
       call Write_InitStatus("kpnBg",0)
       call kpnbgini
       initFlag=.false.
       call Write_InitStatus("kpnBg",1)
    end if

    if (plab2.lt.plab(1)) then
       sigma=0.
    else if (plab2.gt.plab(n)) then
       sigma=sig(n)
    else
       sigma=cl_spline(s2,plab2,successFlag,errorType)
       if (.not.successFlag) call cl_error(errorType,'kpnbg',plab2)
    end if
    kpnbg=sigma
    return

  contains

    subroutine kpnbgini
      use inputGeneral, only: path_to_input
      integer i
      open(13,file=trim(path_to_input)//'/kpn_tot_bg_spl.dat',status='unknown')
      do i=1,n
         read(13,*)plab(i),sig(i)
      end do
      close(13)
      s2=cl_initSpline(plab,sig)
    end subroutine kpnbgini

  end function kpnbg


  !****************************************************************************
  !****f* parametrizationsBarMes/kaonbg
  ! NAME
  ! real function kaonbg(srts,ichin,ichout,isospinc)
  ! PURPOSE
  ! Calculate K- p and K- n cross section im mb
  ! INPUTS
  ! * real, intent(in) :: srts        -- invariant energy (GeV)
  ! * integer, intent(in) :: ichin    -- initial channel
  ! * integer, intent(in) :: ichout   -- final channel
  ! * integer, intent(in) :: isospinc -- switch for isospin of incoming nucleon:
  !   0 -> K- p, 1 -> K- n
  ! NOTES
  ! Enumeration of incoming and outgoing channels is as follows:
  ! * 1 --- K- p     or  K- n
  ! * 2 --- Kbar0 n
  ! * 3 --- Lambda pi0
  ! * 4 --- Sigma+ pi-  or Sigma- pi0 (Sigma0 pi-)
  ! * 5 --- Sigma- pi+
  ! * 6 --- Sigma0 pi0
  !****************************************************************************
  real function kaonbg(srts,ichin,ichout,isospinc)

    use ParticleProperties, only: hadron
    use IDTable, only: Lambda, sigmaResonance
    use output
    use constants, only: mN, mPi, mK
    use twoBodyTools, only: pcm

    real, intent(in) :: srts        ! invariant energy (GeV)
    integer, intent(in) :: ichin    ! initial channel
    integer, intent(in) :: ichout   ! final channel
    integer, intent(in) :: isospinc ! switch for isospin of incoming nucleon
                                    ! 0 -> K- p
                                    ! 1 -> K- n

    real :: mass1,mass2,mass3,mass4,pinitial,pfinal,pfinalo,plab,s
    integer :: ichannel

    if (isospinc.ne.0 .and. isospinc.ne.1) then
       write(*,*) ' wrong input in kaonbg, isospinc = ', isospinc
       stop
    end if

    if (min(ichin,ichout).ne.1) then
       kaonbg=0
       return
    end if

    if (ichin.le.2) then
       mass1=mN
       mass2=mK
    else if (ichin.eq.3) then
       mass1=hadron(lambda)%mass
       mass2=mPi
    else
       mass1=hadron(sigmaResonance)%mass
       mass2=mPi
    end if

    if (mass1+mass2.ge.srts) then
       kaonbg=0
       return
    end if

    if (ichout.le.2) then
       mass3=mN
       mass4=mK
    else if (ichout.eq.3) then
       mass3=hadron(lambda)%mass
       mass4=mPi
    else
       mass3=hadron(sigmaResonance)%mass
       mass4=mPi
    end if

    if (mass3+mass4.ge.srts) then
       kaonbg=0
       return
    end if

    pinitial=pcm(srts,mass1,mass2)
    pfinal=pcm(srts,mass3,mass4)
    s=srts**2
    plab=sqrt((s-mK**2-mN**2)**2/(4.*mN**2)-mK**2)
    !*pfinalo appears in the matrix element and was the final momentum
    !*when Kbar N was initial, therefore:
    if (ichin.le.2) then
       pfinalo=pfinal
    else if (ichout.le.2) then
       pfinalo=pinitial
    else
       write(*,*) 'problems in kaonbg ichin,ichout',ichin,ichout
       stop
    end if

    ichannel=max(ichin,ichout)

    if (isospinc.eq.0) then

       !*K- p in initial or final state

       select case (ichannel)

       case (1)  ! K- p <-> K- p

          if (srts.lt.1.669) then
             kaonbg=235.4/s*(0.0144/(0.0144+pfinalo**2))**0.78
          else if (srts.lt.2.038) then
             kaonbg=44.9*(s/2.7556-1.)**0.23*(2.7556/s)**3.13-3.78
          else
             kaonbg=382.6*(s/3.24-1.)**2.25*(3.24/s)**7.93+3.76
          end if

       case (2)  ! K- p <-> Kbar0 n

          if (srts.lt.1.56) then
             kaonbg=62./s*(0.0324/(0.0324+pfinalo**2))**1.8
          else if (srts.lt.1.66) then
             kaonbg=76.5-410.3*plab+756.5*plab**2-464.5*plab**3
          else if (srts.lt.1.67) then
             kaonbg=0.
          else if (srts.lt.1.75) then
             kaonbg=0.81*exp(-(srts-1.71)**2/0.00036)
          else if (srts.lt.1.95) then
             kaonbg=0.86*exp(-(srts-1.85)**2/0.0027)
          else if (srts.lt.2.04) then
             kaonbg=0.
          else
             kaonbg=43.1*(s/4.-1.)**2.3*(4./s)**7.8
          end if

       case (3)  ! K- p <-> pi0 Lambda

         if (srts.lt.1.521) then
            kaonbg=617.7*pfinal/pinitial/s*(5.625e-03/(5.625e-03+pfinalo**2))**1.49
         else
            kaonbg=617.5*pfinal/pinitial/s*(4.e-04/(4.e-04+pfinalo**2))**0.8
         end if

       case (4)  ! K- p <-> Sigma+ Pi-

          if (srts.lt.1.51) then
             kaonbg=600.*pfinal/pinitial/s*(3.6e-05/(3.6e-05+pfinalo**2))**0.48
          else if (srts.lt.1.67) then
             kaonbg=(15.8-21.1*plab) * (pfinal/pfinalo)**2
          else if (srts.lt.2.15) then
             kaonbg=(0.57*exp(-(srts-1.71)**2/0.0019)+0.72*exp(-(srts-1.98)**2/0.082)) &
                   & * (pfinal/pfinalo)**2
          else
             kaonbg=1.33/plab**1.66 * (pfinal/pfinalo)**2
          end if

       case (5)  ! K- p <-> Sigma- Pi+

          if (srts.lt.1.62) then
             kaonbg=4999.*pfinal/pinitial/s*(8.41e-04/(8.41e-04+pfinalo**2))**1.43
          else if (srts.lt.1.94) then
             kaonbg=(1.85*exp(-(srts-1.65)**2/0.0038)+0.34*exp(-(srts-1.84)**2/0.0028)) &
                   & * (pfinal/pfinalo)**2
          else
             kaonbg=0.
          end if

       case (6)  ! K- p <-> Sigma0 Pi0

          if (srts.lt.3.5) then
             kaonbg=2640.*pfinal/pinitial/s*(1.936e-03/(1.936e-03+pfinalo**2))**1.66
          else
             kaonbg=0.
          end if

       case default

          write(*,*) ' in antiKaonNucleon: wrong channel'
          write(*,*) ' ichin,ichout,isospinc :', ichin,ichout,isospinc
          stop

       end select

    else

       !*K- n in initial or final state

       select case (ichannel)

       case (1)  ! K- n <-> K- n

          if (srts.lt.1.61) then
             kaonbg=4.381
          else
             kaonbg=325.*(s/2.56-1.)**1.28*(2.56/s)**8.84+3.31
          end if

       case (4)  ! K- n <-> Sigma- Pi0 (Sigma0 Pi-)

          if (srts.lt.1.55) then
             ! from K- p -> Sigma- Pi+ :
             kaonbg=4999.*pfinal/pinitial/s*(8.41e-04/(8.41e-04+pfinalo**2))**1.43
          else if (srts.lt.1.638) then
             kaonbg=(45.197-27.331*srts) * (pfinal/pfinalo)**2
          else if (srts.lt.2.006) then
             kaonbg=0.72*exp(-(srts-1.74)**2/0.019)
          else
             kaonbg=8.*(s/4.-1.)**1.2*(4./s)**7.2 * (pfinal/pfinalo)**2
          end if

       case default

          write(*,*) ' in antiKaonNucleon: wrong channel'
          write(*,*) ' ichin,ichout,isospinc :', ichin,ichout,isospinc
          stop

       end select

    end if

  end function kaonbg


  !****************************************************************************
  !****f* parametrizationsBarMes/kaonbgt
  ! NAME
  ! real function kaonbgt(plab2)
  ! PURPOSE
  ! Calculate K- p  total cross section im mb
  !****************************************************************************
!   real function kaonbgt(plab2)
!     use output
!
!     integer n
!     !      real b2,error1,error2
!     parameter(n=44)
!
!     real, save ::  plab(n),sig(n)!,sig2(n)
!
!     real plab2,sigma
!
!     logical,save :: initFlag=.true. ! checks wether this routine is initialized
!
!     logical :: successFlag
!     integer :: errorType
!
!     If(initFlag) then
!        call Write_InitStatus("kaonBgT",0)
!        call kaonbgtini
!        initFlag=.false.
!        call Write_InitStatus("kaonBgT",1)
!     end if
!
!     if(plab2.lt.plab(1)) then
!        sigma=0.
!     else if(plab2.gt.plab(n)) then
!        sigma=sig(n)
!     else
!        sigma=cl_spline(s4,plab2,successFlag,errorType)
!        if(.not.successFlag) call cl_error(errorType,'kppbg',plab2)
!     end if
!     kaonbgt=sigma
!     return
!
!   contains
!
!     subroutine kaonbgtini
!       use inputGeneral, only : path_to_input
!       integer i
!       real error1,error2,b2
!
!       open (13,file=trim(path_to_input)//'/kmp_tot_bg.dat',  status='unknown')
!       do i=1,n
!          read(13,*)plab(i),b2,sig(i),error1,error2
!       end do
!       close(13)
!       s4=cl_initSpline(plab,sig)
!     end subroutine kaonbgtini
!   end function  kaonbgt


  !****************************************************************************
  !****f* parametrizationsBarMes/huang
  ! NAME
  ! function huang(srts) result(sig)
  ! PURPOSE
  ! cross sections for pi p -> K Sigma in mb:
  ! * sig(1) = pi^{+}  p  ->   K^{+} Sigma+
  ! * sig(2) = pi^{0}  p  ->   K^{+} Sigma0 or pi^{-}  p  ->   K^{0} Sigma0
  ! * sig(3) = pi^{-}  p  ->   K^{0} Sigma0
  ! * sig(4) = pi^{-}  p  ->   K^{+} Sigma-
  ! NOTES
  ! See Effenberger PHD, page 236
  !****************************************************************************
  function huang(srts) result(sig)

    real, intent(in)  :: srts
    real :: sig(4)

    if (srts <= 1.688) then
       sig = 0.
       return
    end if

    ! pi+ p -> Sigma+ K+
    SIG(1)=0.03591*(SRTS-1.688)**0.9541/((SRTS-1.89)**2+0.01548) +0.1594*(SRTS-1.688)**0.01056/((SRTS-3.0)**2+0.9412)
    ! pi0 p -> Sigma0 K+
    SIG(2)=0.003978*(SRTS-1.688)**0.5848/((SRTS-1.74)**2+0.00667)   +0.04709*(SRTS-1.688)**2.165/((SRTS-1.905)**2+0.006358)
    ! pi- p -> Sigma0 K0
    SIG(3)=0.05014*(SRTS-1.688)**1.2878/((SRTS-1.73)**2+0.006455)
    ! pi- p -> Sigma- K+
    SIG(4)=0.009803*(SRTS-1.688)**0.6021/((SRTS-1.742)**2+0.006583)+0.006521*(SRTS-1.688)**1.4728/((SRTS-1.94)**2+0.006248)

    !*      do k=1,4
    !*         sig(k)=sig(k)*2.
    !*      end do

  end function huang


  !****************************************************************************
  !****f* parametrizationsBarMes/huanglam
  ! NAME
  ! real function huanglam(srts)
  ! PURPOSE
  ! cross section for pi- p -> Lambda K0 in mb
  !****************************************************************************
  real function huanglam(srts)
    real, intent(in) :: srts
    real :: s0
    if (srts < 1.612) then
       huanglam = 0.
    else
      S0=1.612
      huanglam = 0.007665 * (srts-S0)**0.1341 / ((srts-1.72)**2+0.007826)
      !*      huanglam=huanglam*2.
    end if
  end function huanglam


  !****************************************************************************
  !****f* parametrizationsBarMes/huanglamd
  ! NAME
  ! real function huanglamd(srts)
  ! PURPOSE
  ! cross section for pi- Delta++ -> Lambda K+ in mb
  !
  ! Clebsch for (pi^- D++ -> I=1/2 I_z=1/2)   is 1/2
  !
  ! Therefore : Xsection( pi- D++ -> Lambda k+)=1/2 Xsection( pion Delta -> Lambda kaon)
  !
  ! => Xsection( pion Delta -> R(I=1/2)-> Sigma kaon) = 2 * Xsection( pi- D++ -> S0 k+)
  !****************************************************************************
  real function huanglamd (srts)
    real, intent(in) :: srts
    if (srts <= 1.612) then
       huanglamd = 0.
    else
      huanglamd = 0.006545 * (srts-1.612)**0.7866 / ((srts-1.72)**2+0.004852)
      !*      huanglamd=huanglamd*2.
    end if
  end function huanglamd


  !****************************************************************************
  !****f* parametrizationsBarMes/huangd
  ! NAME
  ! real function huangd(srts)
  ! PURPOSE
  ! cross section for pi- Delta++ -> Sigma0 K+ in mb
  !
  ! Assumption: all processes happen via I=1/2 resonance
  !
  ! Clebsch for (pi^- D++ -> I=1/2 I_z=1/2)   is 1/2
  !
  ! Clebsch for (k^+ S0 -> I=1/2 I_z=1/2)     is 1/3
  !
  ! Therefore : Xsection( pi- D++ -> S0 k+)=1/2*1/3* Xsection( pion Delta -> R(I=1/2)-> Sigma kaon)
  !
  ! => Xsection( pion Delta -> R(I=1/2)-> Sigma kaon) = 6 * Xsection( pi- D++ -> S0 k+)
  !****************************************************************************
  real function huangd (srts)
    real, intent(in) :: srts
    !*pi- D++ -> S0 k+
    if (srts <= 1.688) then
       huangd = 0.
    else
      huangd = 0.004959 * (srts-1.688)**0.7785 / ((srts-1.725)**2+0.008147)
      !*      huangd=huangd*2.
    end if
  end function huangd


  !****************************************************************************
  !****f* parametrizationsBarMes/sibirtpi
  ! NAME
  ! real function sibirtpi(srts)
  ! INPUTS
  ! real :: srts    --- sqrt(s)
  ! PURPOSE
  ! cross section for  pi N -> N K Kbar
  !****************************************************************************
  real function sibirtpi(srts)
    use constants, only: mN, mK
    real,intent(in) :: srts
    real :: s0,ss,ss0

    S0 = mN + 2*mK
    if (srts <= S0) then
       sibirtpi = 0.
    else
      SS=srts**2
      SS0=S0**2
      !*cross section for pi- p -> n K0 k0bar,therefore factor 2.
      sibirtpi = 2. * 1.121 * (1.-SS0/SS)**1.86 * (SS0/SS)**2
    end if
  END function sibirtpi



  !****************************************************************************
  !****f* parametrizationsBarMes/golub_omega
  ! NAME
  ! function golub_omega (srts) result(sigma)
  ! INPUTS
  ! * real :: srts = sqrt(s) in GeV
  ! * real :: m_V = mass of vector meson V in GeV, only used for sigma(3) and sigma(4)   (unused!)
  ! NOTES
  ! Calculates the cross sections (in mb) for interactions of omegas with nucleons.
  ! * sigma(1): pi- p -> V n
  ! * sigma(2): pi- p -> V pi N
  ! * sigma(3): V N -> V N         (unused!)
  ! * sigma(4): V N -> inelastic   (unused!)
  !
  ! Original parameterizations are taken from:
  ! * J. Cugnon et al., PRC 41, 1701 (1990) for sigma(pi^- p -> omega n);
  ! * A. Sibirtsev et al., Z. Phys. A 358, 357 (1997) for sigma(pi^- p -> phi n), sigma(pi^+ p -> omega X)
  !   and sigma(pi^+ p -> phi X).
  ! * Some of these parameterizations are also used in Golubeva NPA625 (1997) 832.
  ! * See also Effenberger Diss. (A.20)...(A.25)
  !
  ! This routine needs a careful check yet (see subtractions between various cross sections below)!!!!
  !****************************************************************************
  function golub_omega (srts) result (sigma)
    use particleProperties, only: hadron
    use idTable, only: omegaMeson
    use constants, only: mN, mPi
    use twoBodyTools, only: p_lab

    real, intent(in) :: srts !, m_V
    real :: sigma(2)

    real :: ppi0,ppi,x,mass1 ! ,plab

       mass1=hadron(omegaMeson)%mass

       if (srts > mass1 + mN) then
          ! threshold value for pion momentum in lab frame
          ppi0 = p_lab(mass1+mN,mPi,mN)
          ppi  = p_lab(srts,mPi,mN)
          ! Parameterization of sigma(pi^- p -> omega n) from J. Cugnon et al., PRC 41, 1701 (1990):
          sigma(1)=13.76*(ppi-ppi0)/(ppi**3.33-1.07)
       else
          sigma(1)=0.
       end if

       if (srts > mass1 + mPi + mN) then
          x=srts**2/(mass1+mN)**2
          ! sigma(pi^+ p -> omega X) - sigma(pi^- p -> omega n), where
          ! sigma(pi^+ p -> omega X) is taken from A. Sibirtsev et al., Z. Phys. A 358, 357 (1997):
          sigma(2)=max(4.8*(x-1)**1.47*x**(-1.26)-sigma(1),0.)
       else
          sigma(2)=0
       end if

       !if(srts.gt.mass1+mN) then
!           plab = p_lab(srts,m_V,mN)
!           sigma(3) = 20./(1.+plab)
!           sigma(4) = 11. + 9./plab-sigma(3)
       !else
       !   sigma(3)=20.
       !   sigma(4)=100.
       !end if

  end function golub_omega


  !****************************************************************************
  !****f* parametrizationsBarMes/golub_phi
  ! NAME
  ! function golub_phi (srts, m_V) result(sigma)
  ! INPUTS
  ! * real :: srts = sqrt(s) in GeV
  ! * real :: m_V = mass of vector meson V in GeV, only used for sigma(3) and sigma(4)
  ! NOTES
  ! Calculates the cross sections (in mb) for interactions of phis with nucleons.
  ! * sigma(1): pi- p -> V n
  ! * sigma(2): pi- p -> V pi N
  ! * sigma(3): V N -> V N
  ! * sigma(4): V N -> inelastic
  !
  ! Original parameterizations are taken from:
  ! * A. Sibirtsev et al., Z. Phys. A 358, 357 (1997) for sigma(pi^- p -> phi n), sigma(pi^+ p -> omega X)
  !   and sigma(pi^+ p -> phi X).
  ! * Some of these parameterizations are also used in Golubeva NPA625 (1997) 832.
  ! * See also Effenberger Diss. (A.26)...(A.31)
  !
  ! This routine needs a careful check yet (see subtractions between various cross sections below)!!!!
  !****************************************************************************
  function golub_phi (srts, m_V) result(sigma)
    use particleProperties, only: hadron
    use idTable, only: phi
    use constants, only: pi, mN, mPi
    use twoBodyTools, only: pCM, p_lab

    real, intent(in) :: srts, m_V
    real :: sigma(4)

    real :: ppi,pphi,x,plab,mass1

       mass1=hadron(phi)%mass

       if (srts > mass1 + mN) then
          ppi  = pCM(srts,mPi,mN)
          pphi = pCM(srts,mass1,mN)
          ! Parameterization of sigma(pi^- p -> phi n) from A. Sibirtsev et al., Z. Phys. A 358, 357 (1997):
          sigma(1)=0.00588/2.*pi**2/srts*0.99**2/((srts-1.8)**2+0.99**2/4.)*pphi/ppi**2
       else
          sigma(1)=0.
       end if

       sigma(2)=0
       if (srts > mass1 + mPi + mN) then
          x=srts**2/(mass1+mN)**2
          if (x > 2.3) then
             ! sigma(pi^+ p -> phi X) - sigma(pi^- p -> phi n), where
             ! sigma(pi^+ p -> phi X)  is taken from A. Sibirtsev et al., Z. Phys. A 358, 357 (1997)
             sigma(2)=0.09*(x-1.3-1.)**2.54*(x-1.3)**(-2.1)-sigma(1)
          end if
       end if
       !sigma(2)=0.

       !if(srts.gt.mass1+mN) then
          plab = p_lab(srts,m_V,mN)
          sigma(3) = 10./(1.+plab)
          sigma(4) = 5. + 4.5/plab-sigma(3)
       !else
       !   sigma(3)=20.
       !   sigma(4)=100.
       !end if

  end function golub_phi


  !****************************************************************************
  !****f* parametrizationsBarMes/omegaN_lykasov
  ! NAME
  ! real function omegaN_lykasov(srts,m,ichannel)
  ! INPUTS
  ! * real,    intent(in) :: srts --- sqrt(s) in GeV
  ! * real,    intent(in) :: m --- mass of the omega meson in GeV
  ! * integer, intent(in) :: ichannel --- 1=elastic, 2=inelastic
  ! NOTES
  ! calculates omega-nucleon cross section in mb:
  ! * ichannel=1 : elastic
  ! * ichannel=2 : inelastic
  !
  ! References:
  ! * Lykasov,Cassing,Sibirtsev,Rzyanin, Eur. Phys. J. A6 (1999) 71, nucl-th/9811019.
  ! * see also Muehlich Diss. (9.1),(9.2)
  !****************************************************************************
  real function omegaN_lykasov(srts,m,ichannel)
    use constants, only: mN

    real,    intent(in) :: srts
    real,    intent(in) :: m
    integer, intent(in) :: ichannel
    real :: pabs

    pabs = sqrt(max(((srts**2-mN**2-m**2)/2./mN)**2-m**2,1e-06))

    if (ichannel == 1) then
      ! omega N (elastic)
      omegaN_lykasov = 5.4 + 10.*exp(-0.6*pabs)
    else if (ichannel == 2) then
      ! omega N (inelastic)
      omegaN_lykasov = 20. + 4./pabs
    end if

  end function omegaN_lykasov


  !****************************************************************************
  !****s* parametrizationsBarMes/JPsiN
  ! NAME
  ! subroutine JPsiN(srts,sigma)
  ! INPUTS
  ! * real :: srts = sqrt(s) in GeV
  ! NOTES
  ! subroutine for calculation of cross sections for interactions of
  ! a J/Psi with a nucleon
  ! OUTPUT
  ! cross sections in mb:
  ! * sigma(1)=J N -> J N
  ! * sigma(2)=J N -> Lambda_c Dbar
  ! * sigma(3)=J N -> Lambda_c D*bar
  ! * sigma(4)=J N -> N D Dbar
  !
  ! Parameterizations of the boson-exchange model calculations
  ! from A. Sibirtsev et al., PRC 63, 044906 (2001)
  !****************************************************************************
  subroutine JPsiN(srts,sigma)
    use particleProperties, only: hadron
    use idTable, only: JPsi,dBar,dStarBar,Lambda_cPlus
    use constants, only: mN
    use twoBodyTools, only: p_lab

    real,   intent(in) :: srts
    real, dimension(1:4), intent(out):: sigma

    real :: mJ,mDbar,mDsbar,mLc,plab
    real :: s,s0,srts0,sig_Dexch,sig_Dsexch

    mJ=hadron(JPsi)%mass
    mDbar=hadron(dBar)%mass
    mDsbar=hadron(dStarBar)%mass
    mLc=hadron(Lambda_cPlus)%mass

    plab=p_lab(srts,mJ,mN)
    s=srts**2

    sigma=0.

    !******* J N -> J N (elastic) ********************************
    !******* fit of the optical model result from A. Sibirtsev
    if (plab.lt.68.7) then
       sigma(1)=1.22*exp(-plab**2/19.46)+0.49*exp(-(plab-68.7)**2/9086.)
    else
       sigma(1)=0.49
    end if

    !******* J N -> Lambda_c Dbar: ********************************
    srts0=mLc+mDbar
    if (srts.gt.srts0) then
       s0=srts0**2
       sig_Dexch=5.477*(s/s0-1.)**0.869*(s0/s)**6.79
       if (srts.lt.4.595) then
          sig_Dsexch=82.83*(s/s0-1.)**1.97*(s0/s)**8.92
       else
          sig_Dsexch=303.4*(srts-srts0)**1.45*exp(-1.34*srts)+0.52
       end if
       sigma(2)=sig_Dexch+sig_Dsexch
    end if

    !******* J N -> Lambda_c D*bar: *******************************
    srts0=mLc+mDsbar
    if (srts.gt.srts0) then
       s0=srts0**2
       sigma(3)=5.2*(s/s0-1.)**0.65*(s0/s)**4.7
    end if

    !******* J N -> N D Dbar: *************************************
    srts0=mN+2.*mDbar
    if (srts.gt.srts0) then
       s0=srts0**2
       sigma(4)=5.85*(s/s0-1.)**1.67*(s0/s)**1.74
    end if

  end subroutine JPsiN


  !****************************************************************************
  !****s* parametrizationsBarMes/weideta
  ! NAME
  ! subroutine weideta(srts,sig)
  ! INPUTS
  ! * real :: srts =sqrt(s)
  ! OUTPUT
  ! * real sig :  parameterization of pi- p -> eta n from Weidmann diploma thesis
  !****************************************************************************
!   subroutine weideta(srts,sig)
!     use particleProperties, only: hadron
!     use idTable, only: eta
!     use constants, only: mN
!
!     real srt0,srts,sig
!
!     srt0=mN+hadron(eta)%mass
!     if(srts.le.srt0) then
!        sig=0.
!     else if(srts.lt.1.5906) then
!        sig=13.07*(srts-srt0)**0.25288
!     else
!        sig=0.14486*(srts-srt0)**(-1.452426)
!     end if
!   end subroutine weideta



!   subroutine twopiback(srts,sigplus,sigminus,signull,sigminull,  signuplus)
!     use constants, only: mN, mPi
!
!     real,intent(in)::srts
!     real,intent(out)::sigplus,sigminus,signull,sigminull,signuplus
!     real::plab,a,b1,b2,b3,b4,plabs
!
!     if(srts>(mN+2.*mPi)) then
!
!        plabs=sqrt(max(((srts**2-mN**2-mPi**2)/(2.*mN))**2- mPi**2,0.))*1000. ! MeV
!
!        if(plabs<=1000.) then  !bis 1 GeV
!           plab=min(plabs,500.)
!
!           a=578.71924
!           b1=-5.4295
!           b2=0.01603
!           b3=-1.44403E-05
!           sigplus=max(a+b1*plab+b2*plab**2+b3*plab**3,0.)/1000. !mb
!
!           a=-2161.32444
!           b1=24.39403
!           b2=-0.09118
!           b3=1.13403E-04
!           sigminus=max(a+b1*plab+b2*plab**2+b3*plab**3,0.)/1000. !mb
!
!           a=-9463.22072
!           b1=98.92348
!           b2=-0.35753
!           b3=4.84901E-04
!           b4=-1.37707E-07
!           signull=max(a+b1*plab+b2*plab**2+b3*plab**3+b4*plab**4,0.)/1000. !mb
!
!           a=-1029.9821
!           b1=10.19464
!           b2=-0.03398
!           b3=3.82246E-05
!           sigminull=max(a+b1*plab+b2*plab**2+b3*plab**3,0.)/1000. !mb
!
!           a=-2226.21419
!           b1=26.85298
!           b2=-0.1201
!           b3=2.34545E-04
!           b4=-1.66873E-07
!           signuplus=max(a+b1*plab+b2*plab**2+b3*plab**3+   b4*plab**4,0.)/1000. !mb
!
!           if(plabs>=500.) then
!              sigplus=max(sigplus*((1000.-plabs)/500.),0.)
!              sigminus=max(sigminus*((1000.-plabs)/500.),0.)
!              signull=max(signull*((1000.-plabs)/500.),0.)
!              sigminull=max(sigminull*((1000.-plabs)/500.),0.)
!              signuplus=max(signuplus*((1000.-plabs)/500.),0.)
!           end if
!
!        else
!           sigplus=0.
!           sigminus=0.
!           signull=0.
!           sigminull=0.
!           signuplus=0.
!        end if
!
!     else
!        sigplus=0.
!        sigminus=0.
!        signull=0.
!        sigminull=0.
!        signuplus=0.
!     end if
!
!     return
!   end subroutine twopiback


  !****************************************************************************
  !****s* parametrizationsBarMes/piN_to_strangeBaryon_kaon_pion
  ! NAME
  ! subroutine piN_to_strangeBaryon_kaon_pion(sqrtS,crossSections_piMinus_p,crossSections_piPlus_p)
  ! INPUTS
  ! * real :: srts =sqrt(s)
  ! RESULT
  ! * real, dimension(1:4) :: crossSections_piPlus_p
  ! * real, dimension(1:7) :: crossSections_piMinus_p
  !
  !
  ! crossSections_piPlus_p gives the cross section for pi^+ p -> ...:
  ! *  1 : Lambda K^+ pi^+
  ! *  2 : Sigma^0 K^+ pi^+
  ! *  3 : Sigma^+ K^+ pi^0
  ! *  4 : Sigma^+ K^0 pi^+
  !
  !
  ! crossSections_piMinus_p gives the cross section for pi^- p -> ...:
  ! *  1 : Lambda K^0 pi^0
  ! *  2 : Lambda K^+ pi^-
  ! *  3 : Sigma^0 K^0 pi^0
  ! *  4 : Sigma^0 K^+ pi^-
  ! *  5 : Sigma^- K^+ pi^0
  ! *  6 : Sigma^- K^0 pi^+
  ! *  7 : Sigma^+ K^0 pi^-
  !
  ! All cross sections in units of mb!
  ! NOTES
  ! This parametrization is based on a fit by K. Gallmeister (12/2008)
  ! * Region in which the fit was performed: p_lab<3.5 GeV.
  ! * Data error bars have not been considered in the fit.
  !****************************************************************************
  subroutine piN_to_strangeBaryon_kaon_pion(sqrtS,crossSections_piMinus_p,crossSections_piPlus_p)

    real, intent(in) :: sqrts
    real, dimension(1:4), intent(out) :: crossSections_piPlus_p
    real, dimension(1:7), intent(out) :: crossSections_piMinus_p
    ! Fit parameters by K. Gallmeister
    real, parameter  :: s0_lambdaKPi =  1.850
    real, parameter  :: s0_SigmaKPi  =  1.923
    real, parameter, dimension(1:7)  :: ai_piMinus_p=(/0.169 ,0.140 , 0.100, 0.0724, 0.0520, 0.117, 0.0514 /)
    real, parameter, dimension(1:4)  :: ai_piPlus_p =(/0.217 ,0.0426, 0.126, 0.0887 /)
    real ::  f_lambdaKPi, f_sigmaKPi

    ! Check thresholds
    if (sqrts.gt.s0_lambdaKPi) then
       f_lambdaKPi  = f((sqrtS/ s0_lambdaKPi)**2)
    else
       f_lambdaKPi  = 0.
    end if

    if (sqrts.gt.s0_SigmaKPi) then
       f_sigmaKPi   = f((sqrtS/ s0_sigmaKPi)**2)
    else
       f_sigmaKPi   = 0.
    end if

    ! Construct cross sections for pi^+ p -> ...
    crossSections_piPlus_p(1)     = f_lambdaKPi
    crossSections_piPlus_p(2:4)   = f_sigmaKPi
    !  * Normalization:
    crossSections_piPlus_p        = crossSections_piPlus_p    * ai_piPlus_p

    ! Construct cross sections for pi^- p -> ...
    crossSections_piMinus_p(1:2)  = f_lambdaKPi
    crossSections_piMinus_p(3:7)  = f_sigmaKPi
    !  * Normalization:
    crossSections_piMinus_p       = crossSections_piMinus_p   * ai_piMinus_p

    contains
      real function f(x)
        real, intent(in) :: x
        ! Fit parameters by K. Gallmeister
        real, parameter  :: A            = 86.027
        real, parameter  :: B            =  2.197
        real, parameter  :: C            =  7.363
        real, parameter  :: eps=1E-5
        if (x.gt.1.+eps) then
           f=A*(x -1.)**B * x**(-C)
        else
           f=0.
        end if
      end function f
  end subroutine piN_to_strangeBaryon_kaon_pion

  !****************************************************************************
  !****s* parametrizationsBarMes/piN_to_strangeBaryon_kaon_pion_matrix
  ! NAME
  ! subroutine piN_to_strangeBaryon_kaon_pion_matrix(sqrtS,pionCharge,nukCharge,matrix_lambda,matrix_sigma)
  !
  ! INPUTS
  ! * real :: srts =sqrt(s)
  ! * integer, intent(in) :: pionCharge     -- charge of incoming pion
  ! * integer, intent(in) ::  nucCharge     -- charge of incoming nucleon
  !
  ! RESULT
  ! Cross section for kaon pion Lambda production :
  ! * real, dimension(0:1,-1:1)  , intent(out) :: matrix_lambda ! 1st Index: kaon charge, 2nd: pion charge
  !
  ! Cross section for kaon pion Sigma production (Charge of Sigma is fixed by total charge):
  ! * real, dimension(0:1,-1:1)  , intent(out) :: matrix_sigma  ! 1st Index: kaon charge, 2nd: pion charge
  !
  ! All cross sections in units of mb!
  !
  ! NOTES
  ! * This parametrization is based on the subroutine piN_to_strangeBaryon_kaon_pion
  ! * No Xsections for pi^0 induced reactions!!
  !****************************************************************************
  subroutine piN_to_strangeBaryon_kaon_pion_matrix(sqrtS,pionCharge,nucCharge,matrix_lambda,matrix_sigma)

    real   , intent(in) :: sqrts
    integer, intent(in) :: pionCharge
    integer, intent(in) ::  nucCharge
    real, dimension(0:1,-1:1)  , intent(out) :: matrix_lambda ! 1st Index: kaon charge, 2nd: pion charge
    real, dimension(0:1,-1:1)  , intent(out) :: matrix_sigma  ! 1st Index: kaon charge, 2nd: pion charge
    integer, parameter :: proton  = 1
    integer, parameter :: neutron = 0

    real, dimension(1:4) :: sigma_piPlus_p
    real, dimension(1:7) :: sigma_piMinus_p

    call piN_to_strangeBaryon_kaon_pion(sqrtS,sigma_piMinus_p,sigma_piPlus_p)


    ! pion nukleon -> kaon pion lambda
    matrix_lambda=0.

    select case (nucCharge)
    case (proton)
       select case (pionCharge)
       case (-1)
          matrix_lambda(0,0) = sigma_piMinus_p(1)
          matrix_lambda(1,-1)= sigma_piMinus_p(2)
       case (1)
          matrix_lambda(1,1) = sigma_piPlus_p(1)
       end select
    case (neutron)
       select case (pionCharge)
       case (-1)
          !pi^- n -> lambda K^0 pi^-      <=> pi^+ p -> lambda^0 k^+ pi^+
          matrix_lambda(0,-1) = sigma_piPlus_p(1)
       case (1)
          !pi^+ n -> lambda K^+ pi^0      <=> pi^- p -> lambda k^0 pi^0
          !pi^+ n -> lambda K^0 pi^+      <=> pi^- p -> lambda k^+ pi^-
          matrix_lambda(1,0) = sigma_piMinus_p(1)
          matrix_lambda(0,1) = sigma_piMinus_p(2)
       end select
    end select


    ! pion nukleon -> kaon pion sigma
    matrix_sigma =0.

    select case (nucCharge)
    case (proton)
       select case (pionCharge)
       case (-1)
          matrix_sigma(0, 0)  = sigma_piMinus_p(3)
          matrix_sigma(1,-1) = sigma_piMinus_p(4)
          matrix_sigma(1, 0) = sigma_piMinus_p(5)
          matrix_sigma(0, 1) = sigma_piMinus_p(6)
          matrix_sigma(0,-1) = sigma_piMinus_p(7)
       case (1)
          matrix_sigma(1,1) = sigma_piPlus_p(2)
          matrix_sigma(1,0) = sigma_piPlus_p(3)
          matrix_sigma(0,1) = sigma_piPlus_p(4)
       end select
    case (neutron)
       select case (pionCharge)
       case (-1)
          !pi^- n -> K^0 pi^- Sigma^0 <=>  pi^+ p -> K^+ pi^+ Sigma^0 -> Channel 2
          !pi^- n -> K^0 pi^0 Sigma^- <=>  pi^+ p -> K^+ pi^0 Sigma^+ -> Channel 3
          !pi^- n -> K^+ pi^- Sigma^- <=>  pi^+ p -> K^0 pi^+ Sigma^+ -> Channel 4
          matrix_sigma(0,-1) = sigma_piPlus_p(2)
          matrix_sigma(0, 0) = sigma_piPlus_p(3)
          matrix_sigma(1,-1) = sigma_piPlus_p(4)
       case (1)
          !pi^+ n -> K^+ pi^0 Sigma^0 <=>  pi^- p -> Sigma^0 K^0 pi^0 -> Channel 3
          !pi^+ n -> K^0 pi^+ Sigma^0 <=>  pi^- p -> Sigma^0 K^+ pi^- -> Channel 4
          !pi^+ n -> K^0 pi^0 Sigma^+ <=>  pi^- p -> Sigma^- K^+ pi^0 -> Channel 5
          !pi^+ n -> K^+ pi^- Sigma^+ <=>  pi^- p -> Sigma^- K^0 pi^+ -> Channel 6
          !pi^+ n -> K^+ pi^+ Sigma^- <=>  pi^- p -> Sigma^+ K^0 pi^- -> Channel 7
          matrix_sigma(1, 0)  = sigma_piMinus_p(3)
          matrix_sigma(0, 1)  = sigma_piMinus_p(4)
          matrix_sigma(0, 0)  = sigma_piMinus_p(5)
          matrix_sigma(1,-1)  = sigma_piMinus_p(6)
          matrix_sigma(1, 1)  = sigma_piMinus_p(7)
       end select
    end select




  end subroutine piN_to_strangeBaryon_kaon_pion_matrix


  !****************************************************************************
  !****s* parametrizationsBarMes/sigma_KbarToXi
  ! NAME
  ! real function sigma_KbarToXi(srts,ichannel)
  ! INPUTS
  ! * real :: srts =sqrt(s)
  ! * integer :: ichannel = outgoing channel
  ! NOTES
  ! Considered channels:
  ! * (ichannel = 1) K^- p --> Xi^- K^+
  ! * (ichannel = 2) K^- p --> Xi^0 K^0
  ! * (ichannel = 3) K^- p --> Xi^- K^0 \pi^+
  ! * (ichannel = 4) K^- p --> Xi^- K^+ \pi^0
  ! * (ichannel = 5) K^- p --> Xi^0 K^+ \pi^-
  !****************************************************************************
  real function sigma_KbarToXi(srts,ichannel)
    use particleProperties, only: hadron
    use IdTable, only: Xi
    use constants, only: mN, mPi, mK

    real,    intent(in) :: srts
    integer, intent(in) :: ichannel

    real :: m_hyp, m_ka, m_nuc, m_pi, plab, s0

    !Fit Function: f(x)=a0*((x/x0)**2-1)**a1*((x0/x)**2)**a2 [mb]

    m_hyp = hadron(Xi)%mass
    m_ka  = mK
    m_nuc = mN
    m_pi  = mPi

    plab = sqrt( (srts**2-(m_ka-m_nuc)**2)*(srts**2-(m_ka+m_nuc)**2) )/(2.*m_nuc)

    select case (ichannel)

    case (1) !K^{-} Proton --> Xi^{-} K^{+}
       s0 = m_hyp+m_ka
       if (srts < (s0+0.0001) .or. (plab < 1.042 ) ) then
          sigma_KbarToXi = 0.0
       else
          sigma_KbarToXi = 3.88563*((plab/1.042)**2-1.)**2.30195*((1.042/plab)**2)**4.48525
       end if

    case (2) !K^{-} Proton --> Xi^{0} K^{0}
       s0 = m_hyp+m_ka
       if (srts < (s0+0.0001) .or. (plab < 1.042 ) ) then
          sigma_KbarToXi = 0.0
       else
          sigma_KbarToXi = 3.9*((plab/1.042)**2-1.)**3*((1.042/plab)**2)**5.5
       end if

    case (3) !K^{-} Proton --> Xi^{-} K^{0} \pi^+
       s0 = m_hyp+m_ka+m_pi
       if (srts < (s0+0.0001) .or. (plab < 1.03535) ) then
          sigma_KbarToXi = 0.0
       else
          sigma_KbarToXi = 50.5*((plab/1.03535)**2-1.)**9.9*((1.03535/plab)**2)**12.5
       end if

    case (4) !K^{-} Proton --> Xi^{-} K^{+} \pi^0
       s0 = m_hyp+m_ka+m_pi
       if (srts < (s0+0.0001) .or. (plab < 1.035) ) then
          sigma_KbarToXi = 0.0
       else
          sigma_KbarToXi = 17.6354*((plab/1.035)**2-1.)**8.9688*((1.035/plab)**2)**11.4948
       end if

    case (5) !K^{-} Proton --> Xi^{0} K^{+} \pi^-
       s0 = m_hyp+m_ka+m_pi
       if (srts < (s0+0.0001) .or. (plab < 1.035) ) then
          sigma_KbarToXi = 0.0
       else
          sigma_KbarToXi = 64664.5*((plab/1.035)**2-1.)**18.4351*((1.035/plab)**2)**24.5073
       end if

    case default

       write(*,*) 'module parametrizationsBarMes/function sigma_KbarToXi: '
       write(*,*) 'wrong channel! ichannel = ',ichannel
       write(*,*) 'STOP'
       STOP

    end select


  end function sigma_KbarToXi


end module parametrizationsBarMes
