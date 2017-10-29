! execute this program by calling:
! ./TryNeutrino.x < /dev/null
! or
! ./TryNeutrino.x < jobTryNu


program TryNeutrino
  use output
  use particleDefinition
  use particleProperties
  use Coll_Pythia
  use CollTools
  !use rhoMassParameter, only : srtFreeRhoMass
  use hadronFormation, only : forceInitFormation

  implicit none

  real, save :: Eel=-1.0   ! sqrt(s) = 1.66

!  real, save :: Eel=-5.0
!  real, save :: Eel=-1.66  ! sqrt(s) = 2
!  real, save :: Eel=-1.80  ! sqrt(s) = 2.063
!  real, save :: Eel=-1.90  ! sqrt(s) = 2.108
!  real, save :: Eel=-2.00  ! sqrt(s) = 2.152
!  real, save :: Eel=-4.33  ! sqrt(s) = 3
!  real, save :: Eel=-8.07  ! sqrt(s) = 4
!  real, save :: Eel=-18.7  ! sqrt(s) = 6
!  real, save :: Eel=-33.64 ! sqrt(s) = 8

  real,save :: srts2 ! = 0.938**2+2*0.938*Eel

  integer, save :: nEv = -100000
!  integer, save :: nEv = -1000000
!  integer, save :: nEv = -10000000


  call forceInitFormation
  call InitParticleProperties

  call ReadInit

  call SetSomeDefaults_PY

  call DoAnaEstimate
  call DoAnaEstimate2
!  call DoPythia

  write(*,*) 'ok'


contains
  subroutine DoPythia

    use histf90
    use hist2Df90
    use CollTools
    IMPLICIT NONE

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

    COMMON/PYINT1/MINT(400),VINT(400)
    integer MINT
    double precision VINT
    SAVE /PYINT1/


    type(histogram) :: hQ2, hW, hEl

    type(histogram2D) :: h2d, h2dXY,h2dnuQ2,h2dTauZ,h2dWth,h2dXY_py,h2X,h2Y

    integer :: iEv,nEv100,nEv10, i

    real :: fak, QQ, costheta, x, y
    real :: W, sqrts, theta, phi, beta(3)

    real MP_P,ml_out, MP_ScalProd4,MP_CosAngle3


    nEv100 = nEv / 100
    nEv10 = nEv / 10


    write(*,*) 'START OF DOPYTHIA'

    call CreateHist(hQ2, "Q2", 0.0, 5.0, 0.1)
    call CreateHist(hW,  "W",  0.0, 5.0, 0.1)
    call CreateHist(hEl,  "Eprime",  0.0, 5.0, 0.1)

!    call CreateHist2D(h2d, "Eprime vs cos(theta)", (/0.0,0.0/),(/1.0,25.0/),(/0.01,0.1/))
    call CreateHist2D(h2d, "Eprime vs cos(theta)", (/0.0,0.0/),(/1.0,5.0/),(/0.01,0.01/))
    call CreateHist2D(h2dXY, "x vs y", (/0.0,0.0/),(/1.0,1.0/),(/0.01,0.01/))
    call CreateHist2D(h2dnuQ2, "nu vs Q2", (/0.0,0.0/),(/20.0,35.0/),(/0.2,0.2/))
    call CreateHist2D(h2dTauZ, "tau vs z", (/0.0,-1.0/),(/1.0,1.0/),(/0.01,0.01/))
    call CreateHist2D(h2dWth, "W vs cos(theta)", (/0.0,0.0/),(/1.0,5.0/),(/0.01,0.01/))

    call CreateHist2D(h2dXY_py, "x vs y (pythia vals)", (/0.0,0.0/),(/1.0,1.0/),(/0.01,0.01/))
    call CreateHist2D(h2X, "x in x-y plane", (/0.0,0.0/),(/1.0,1.0/),(/0.01,0.01/))
    call CreateHist2D(h2Y, "y in x-y plane", (/0.0,0.0/),(/1.0,1.0/),(/0.01,0.01/))

    call PYGIVE('MSEL=0') ! selction of subprocesses
    call PYGIVE('MSUB(10)=1') !

    call PYGIVE('MSTP(111) = 0') ! master switch fragmentation/decay

!    call PYGIVE('MSTP(21) = 3') ! neutral current, Z0 only
    call PYGIVE('MSTP(21) = 5') ! charged current only


    call PYGIVE('CKIN( 1) = 0.1') !                                     ***2***
!    call PYGIVE('CKIN(31) = 0.1') ! srts_min in 3,4body final state    ***2***
    call PYGIVE('CKIN(39) = 0.1') ! W2_min in DIS (loose cut)          ***2***
!    call PYGIVE('CKIN(77) = 0.1') ! W_min for photon-hadron            ***2***

    call PYGIVE('MSTP( 61)=0')               ! master: ISR (QCD/QED)
    call PYGIVE('MSTP( 71)=0')               ! master: FSR (QCD/QED)

!    call PYGIVE('CKIN(3) = 0.1d0')    ! min pThat (allow low pT events)
!    call PYGIVE('CKIN(3) = 0.010d0')  ! min pThat (allow low pT events)
!    call PYGIVE('CKIN(3) = 0.1d0')    ! min pThat (allow low pT events)
!    call PYGIVE('CKIN(5) = 0.75d0')   ! min pThat (singular processes)
!    call PYGIVE('CKIN(5) = 0.050d0')  ! min pThat (singular processes)

!    call PYGIVE('CKIN(5) = 0.25d0')   ! min pThat (singular processes) ***1***


!    call PYGIVE('MSTJ(17)=3') ! number of attemts                        ***3***
!    call PYGIVE('PARJ(32)=0.1') ! min.energy of string(+quark masses)    ***3***
!    call PYGIVE('PARJ(33)=0.5') ! border string fragm / 2 hadron cluster ***3***
!    call PYGIVE('PARJ(34)=1.0') ! -"-                                    ***3***
!    call PYGIVE('PARJ(36)=0.3') ! -"-                                    ***3***

!=======================

!!$    call PYGIVE('PMAS(1,1)=0.010') ! Mass of d quark
!!$    call PYGIVE('PMAS(2,1)=0.010') ! Mass of u quark
!!$
!!$    call PYGIVE('PMAS(150,1)=0.020') ! Mass of diquark
!!$    call PYGIVE('PMAS(152,1)=0.020') ! Mass of diquark
!!$    call PYGIVE('PMAS(153,1)=0.020') ! Mass of diquark
!!$    call PYGIVE('PMAS(156,1)=0.020') ! Mass of diquark

!=======================


!    call PYGIVE('PMAS(154,1)=1.200') ! Mass of n0
!    call PYGIVE('PMAS(157,1)=1.200') ! Mass of p+

    call PYGIVE('MSTP(32)= 5') ! Q2 definition: -t

!!    call PYGIVE('CKIN(6)= 0.0005') ! singular if masses below this cut
!    call PYGIVE('CKIN(6)= 0.1') ! singular if masses below this cut
    call PYGIVE('CKIN(6)= 0.000') ! singular if masses below this cut   ***4***

!    call PYGIVE('CKIN(27)= -0.95') ! cut on cos(theta_hat)
!    call PYGIVE('CKIN(28)=  0.95') ! cut on cos(theta_hat)



!    call PYGIVE('PARP(91)=0.00') ! width intrinsic kT      ***5***
    call PYGIVE('PARP(91)=0.44') ! width intrinsic kT


!    call PYINIT('CMS','nu_e','p', sqrt(srts2))
    call PYINIT('CMS','nu_e','n', sqrt(srts2))

    do iEv=1,nEv

       CALL WriteStatusMC(iEv,nEv,nEv100,nEv10)

       call PYEVNT
       if (MINT(51).eq.2) cycle
       if (iEv.lt.11) call PYLIST(2)


       if (iEv.lt.11) then
          call PYGIVE('VINT( 1)=') ! sqrt(s)
          write(*,*) '== tau:'
          call PYGIVE('VINT(11)=') ! taumin
          call PYGIVE('VINT(21)=') ! tau
          call PYGIVE('VINT(31)=') ! taumax
          write(*,*) '== y:'
          call PYGIVE('VINT(12)=') ! ystmin
          call PYGIVE('VINT(22)=') ! yst
          call PYGIVE('VINT(32)=') ! ystmax
          write(*,*) '== cos(theta):'
          call PYGIVE('VINT(13)=') ! theta: ctnmin
          call PYGIVE('VINT(14)=') ! theta: ctpmin
          call PYGIVE('VINT(23)=') ! theta: cth
          call PYGIVE('VINT(33)=') ! theta: ctnman
          call PYGIVE('VINT(34)=') ! theta: ctpmax
          write(*,*) '== x_1,2, Q2:'
          call PYGIVE('VINT(41)=') ! x_1
          call PYGIVE('VINT(42)=') ! x_2
          call PYGIVE('VINT(52)=') ! inner Q2
          call PYGIVE('VINT(54)=') ! outer Q2
          write(*,*) '== shat,that,uhat:'
          call PYGIVE('VINT(44)=')
          call PYGIVE('VINT(45)=')
          call PYGIVE('VINT(46)=')

          call PYGIVE('MINT(43)=')
       endif

       if (VINT(14).ne.0.0) then
          call PYGIVE('VINT(14)=') ! taumax
       endif

       call MP_Set4(1, P(1,5),P(1,1),P(1,2),P(1,3),P(1,4))
       call MP_Set4(2, P(2,5),P(2,1),P(2,2),P(2,3),P(2,4))

       call MP_CalcROBO(2,2, sqrts, theta, phi, beta)

       call MP_Set4(3, P(9,5),P(9,1),P(9,2),P(9,3),P(9,4))

       call MP_ROBO_INV(1,3,theta, phi, beta(1),beta(2),beta(3))

       call MP_Copy(1,4)
       call MP_AddMom(3,4, -1.0)

       if (iEv.lt.11) call MP_write(6,1,4)


!       call AddHist(hQ2, PARI(22), 1.0)
!       call AddHist(hW, PARI(13), 1.0)


!       if (K(9,2).ne.K(1,2)) then
!          write(*,*) ' ooops'
!          call PYLIST(2)
!       endif


       QQ = -MP_ScalProd4(4,4)

!       write(*,*) QQ, PARI(22) ! <-- not equal !!!

       W = (P(1,4)+P(2,4)-P(9,4))**2
       do i=1,3
          W = W - (P(1,i)+P(2,i)-P(9,i))**2
       enddo
       W = sqrt(W)

!       write(*,*) W
!       call AddHist(hW, PARI(13), 1.0)
       call AddHist(hW, W, 1.0)

       call AddHist(hQ2, QQ, 1.0)
       call AddHist(hEl, MP_P(3,4), 1.0)

       ml_out=0.

       costheta = (ml_out**2+QQ-2.*MP_P(1,4)*MP_P(3,4))/(-2.*MP_P(1,4)*sqrt(MP_P(3,4)**2-ml_out**2))

!       write(*,*) MP_CosAngle3(1,3),costheta ! <-- is equal !!

       call AddHist2d(h2d, (/costheta,MP_P(3,4)/), 1.0)
       call AddHist2d(h2dWth, (/costheta,W/), 1.0)

       X = QQ/(2*0.938*MP_P(4,4))
       Y = MP_P(4,4)/MP_P(1,4)

       call AddHist2d(h2dXY, (/x,y/), 1.0)
       call AddHist2d(h2dXY_py, (/VINT(21),VINT(22)/), 1.0)

       call AddHist2d(h2X, (/x,y/), 1.0,VINT(21))
       call AddHist2d(h2Y, (/x,y/), 1.0,VINT(22))


       call AddHist2d(h2dnuQ2, (/MP_P(4,4),QQ/), 1.0)

!       write(*,*) PARI(33),PARI(34)


       if (iEv.lt.11) then
          write(*,*) 'x,y,costheta:',x,y,costheta,QQ
       end if

       call AddHist2d(h2dTauZ, (/VINT(21),VINT(23)/), 1.0)

    end do

    call PYSTAT(1)

    fak = PARI(2)*1e11
!    fak = PARI(2)/((CKIN(36)-CKIN(35)))
!    fak = PARI(2)/((CKIN(36)-CKIN(35))*(CKIN(40)-CKIN(39)))

    call WriteHist(hQ2,101,mul=fak,file="Neutrino.h.Q2.dat")
    call WriteHist(hW, 102,mul=fak,file="Neutrino.h.W.dat")
    call WriteHist(hEl, 103,mul=fak,file="Neutrino.h.El.dat")

    call WriteHist2D_Gnuplot(h2d,201,mul=fak,add=1e-20,file="Neutrino.h2D.theta_Eprime.dat")
    call WriteHist2D_Gnuplot(h2dXY,202,mul=fak,add=1e-20,file="Neutrino.h2D.X_Y.dat")
    call WriteHist2D_Gnuplot(h2dnuQ2,203,mul=fak,add=1e-20,file="Neutrino.h2D.nu_Q2.dat")
    call WriteHist2D_Gnuplot(h2dTauZ,204,mul=fak,add=1e-20,file="Neutrino.h2D.tau_z.dat")
    call WriteHist2D_Gnuplot(h2dWth,205,mul=fak,add=1e-20,file="Neutrino.h2D.theta_W.dat")

    call WriteHist2D_Gnuplot(h2dXY_py,206,mul=fak,add=1e-20,file="Neutrino.h2D.X_Y_Py.dat")
    call WriteHist2D_Gnuplot(h2X,206,DoAve=.true.,file="Neutrino.h2D.X_Y_xAve.dat")
    call WriteHist2D_Gnuplot(h2Y,206,DoAve=.true.,file="Neutrino.h2D.X_Y_yAve.dat")

  end subroutine DoPythia

  !====================================================

  real function AnaXS_def(iC,Ein,x,y)
    implicit none
    real, intent(in) :: Ein,x,y
    integer, intent(iN) :: iC

    real :: nu,Q2,srts2,sigma0
    real, dimension(-25:25) :: xq

    real, parameter :: GF=1.166E-5 ! in GeV-2
    real, parameter :: hc2 = 0.389 ! in GeV2 mb


    srts2 = 0.938**2+2*0.938*Ein
    sigma0 = GF**2*srts2*hc2/(2*3.14)  ! in mb !!!
    nu = y*Ein
    Q2 = 2*0.938*nu * x
    select case(iC)
    case (0)
       call pypdfl(2112, x,Q2,xq)
    case (1)
       call pypdfl(2212, x,Q2,xq)
    case default
       write(*,*) 'wrong charge of target. stop!'
       stop
    end select

    AnaXS_def = 2*sigma0 * (xq(1)+ (1-y)**2*xq(-2)) ! x(d + (1-y)**2ubar)
    return
  end function AnaXS_def

  !====================================================

  real function AnaXS_Aivazis(iC,Ein,x,y)
    implicit none
    real, intent(in) :: Ein,x,y
    integer, intent(iN) :: iC

    real :: nu,Q2,srts2,sigma0,chi,eta,coshpsi,dp,dm
    real :: h1,h2,h3, c1,c2,c3
    real, dimension(-25:25) :: xq

    real, parameter :: GF=1.166E-5 ! in GeV-2
    real, parameter :: hc2 = 0.389 ! in GeV2 mb

    real, parameter :: mq(6) = (/0.03,0.03,99.,99.,99.,99./) ! quark mass in GeV

    integer :: iq1,iq2 ! 1,2,... = d,u,...
    integer :: il      ! the incoming lepton
    integer :: nl      ! no of spin states of incoming lepton
    real :: LmnWmn

    AnaXS_Aivazis = 0.0

    il = 1
    nl = 1

    srts2 = 0.938**2+2*0.938*Ein
    sigma0 = GF**2*srts2*hc2/(2*3.14)  ! in mb !!!
    nu = y*Ein
    Q2 = 2*0.938*nu * x


    iq1 = 1
    iq2 = 2

    !... some kinematical variables

    eta = 2*x/(1d0+sqrt(1d0+(2*0.938*x)**2/Q2))
    chi = eta*((Q2-mq(iq1)**2+mq(iq2)**2)+Delta(-Q2,mq(iq1)**2,mq(iq2)**2))/(2*Q2)
    coshpsi=((eta*0.938)**2-Q2+2*eta*(srts2-0.938**2))/((eta*0.938)**2+Q2)
    dp = (coshpsi+1)/2d0
    dm = (coshpsi-1)/2d0

!    write(*,*) x,y
!    write(*,*) eta,chi

    select case(iC)
    case (0)
       call pypdfl(2112, chi,Q2,xq)
    case (1)
       call pypdfl(2212, chi,Q2,xq)
    case default
       write(*,*) 'wrong charge of target. stop!'
       stop
    end select

    !... factors used in LmunuWmunu:

    h1 = (Q2*chi**2*dm+(mq(iq1)*eta)**2*dp+eta*chi*Q2)*(Q2*chi**2*dp+(mq(iq1)*eta)**2*dm)/(2*eta*chi)**2
    h2 = (Q2*chi**2*dp+(mq(iq1)*eta)**2*dm-eta*chi*Q2)*(Q2*chi**2*dm+(mq(iq1)*eta)**2*dp)/(2*eta*chi)**2
    h3 = mq(iq1)*mq(iq2)*Q2/2

    !... now the couplings:

    c1 = (gRa(iq1)*gRl(il))**2+(gLa(iq1)*gLl(il))**2
    c2 = (gRa(iq1)*gLl(il))**2+(gLa(iq1)*gRl(il))**2
    c3 = gRa(iq1)*gLa(iq1)*(gRl(il)**2+gRl(il)**2)

    !LmnWmn = 2**6/(nl*Q2*Delta(-Q2,mq(iq1)**2,mq(iq2)**2)) * (c1*h1+c2*h2-c3*h3)

    !... approx: (summed over d and ubar)
    LmnWmn = 2**6/(nl*Q2*Delta(-Q2,mq(iq1)**2,mq(iq2)**2)*chi) * c1 * &
         & (h1*xq(1) + h2*xq(-2))

    !please note: pypdfl returns x*f(x,Q2)

!    write(*,*) c1,h1
!    write(*,*) c2,h2
!    write(*,*) c3,h3
!    stop

    AnaXS_Aivazis = (4*3.14*Q2)/(2*(srts2-0.938**2))*(0.938*Ein*y)/(8*3.14**2) * LmnWmn ! now: dsigma/dxdy

!    AnaXS_def = 2*sigma0 * (xq(1)+ (1-y)**2*xq(-2)) ! x(d + (1-y)**2ubar)
    return

  end function AnaXS_Aivazis

  real function Delta(a,b,c)
    real, intent(in) :: a,b,c
    Delta=(a**2+b**2+c**2-2*(a*b+a*c+b*c))
    if (Delta.lt.0.0) then
       write(*,*) 'Problem in delta. stop.'
       stop
    endif
    Delta=sqrt(Delta)
    return
  end function Delta

  real function gLa(iq)
    integer, intent(in) :: iq
    gLa = 0.0
    if (iq.gt.0) then
       gLa = 1.0 ! here we have to add V_{ab} !!!
    else if (iq.lt.0) then
       gLa = 0.0
    else
       write(*,*) 'iq=0 in gLa. stop!'
       stop
    endif
    return
  end function gLa

  real function gRa(iq)
    integer, intent(in) :: iq
    gRa = gLa(-iq)
    return
  end function gRa

  real function gLl(il)
    integer, intent(in) :: il
    gLl = 0.0
    if (il.gt.0) then
       gLl = 1.0
    else if (il.lt.0) then
       gLl = 0.0
    else
       write(*,*) 'il=0 in gLl. stop!'
       stop
    endif
    return
  end function gLl

  real function gRl(il)
    integer, intent(in) :: il
    gRl = gLl(-il)
    return
  end function gRl


  !====================================================


  subroutine DoAnaEstimate

    use CollTools
    IMPLICIT NONE

    integer :: iX, iY
    real :: x,y, Q2, nu, Eprime,cost, sigma0,sigma1,sigma2


    call SetSomeDefaults_PY

    do iX=2,99,2
       x = iX*0.01
       do iY=2,99,2
          y = iY*0.01

          nu = y*Eel
          Q2 = 2*0.938*nu * x
          Eprime = Eel - nu
!          cost = 1.0 -2*y ! Halzen Martin, wrong !!!!!
          cost = 1- Q2/(2*Eel*Eprime)

          if (cost.lt.-1) then
             write(122,'(2F13.5,4E13.5,5e13.5)') x,y,Q2,nu,Eprime,cost,&
               & 1e-20,1e-20,1e-20,1e-20,1e-20
             write(222,'(2F13.5,4E13.5,5e13.5)') x,y,Q2,nu,Eprime,cost,&
               & 1e-20,1e-20,1e-20,1e-20,1e-20
             cycle
          end if

          sigma1 = AnaXS_def(1,Eel,x,y) ! nu p
          sigma2 = AnaXS_def(0,Eel,x,y) ! nu n

          write(122,'(2F13.5,4E13.5,5e13.5)') x,y,Q2,nu,Eprime,cost,&
               & 0d0, sigma1,sigma2, 2*0.938*nu*Eel, 0.938*nu/(Eel-nu)

          sigma1 = AnaXS_Aivazis(1,Eel,x,y) ! nu p
          sigma2 = AnaXS_Aivazis(0,Eel,x,y) ! nu n

          write(222,'(2F13.5,4E13.5,5e13.5)') x,y,Q2,nu,Eprime,cost,&
               & 0d0, sigma1,sigma2, 2*0.938*nu*Eel, 0.938*nu/(Eel-nu)
       end do
       write(122,*)
       write(222,*)

    end do

  end subroutine DoAnaEstimate

  !====================================================

  subroutine DoAnaEstimate2

    use CollTools
    IMPLICIT NONE

    integer :: iE, ict
    real :: x,y, Q2, nu, Eprime,cost, sigma1,sigma2, fak

    logical:: skip =.false.


    do ict=0,99
       cost = ict*0.01

!!$    do ict=5,9
!!$       cost = ict*0.1

!!$       do iE=1,99
!!$          Eprime = iE*0.01*Eel
       do iE=0,200
          Eprime = iE*0.01
          nu = Eel-Eprime
          y = nu/Eel

          Q2 = 2*Eel*Eprime*(1-cost)
          x = Q2/(2*0.938*nu)

          skip = .false.
          if (x.gt.1.0) skip = .true.
          if (x.le.0.0) skip = .true.

          if (skip) then
             sigma1=0
             sigma2=0
          else
             sigma1 = AnaXS_def(1,Eel,x,y) ! nu p
             sigma2 = AnaXS_def(0,Eel,x,y) ! nu n
          endif

          fak = (0.938*Eel*y)/(Eprime)

          write(124,'(2F13.5,1P,4E13.5,5e13.5)') x,y,Q2,nu,Eprime,cost,&
               & 0d0, sigma1,sigma2
          write(125,'(2F13.5,1P,4E13.5,5e13.5)') x,y,Q2,nu,Eprime,cost,&
               & 0d0, sigma1/fak,sigma2/fak

          if (skip) then
             sigma1=0
             sigma2=0
          else
             sigma1 = AnaXS_Aivazis(1,Eel,x,y) ! nu p
             sigma2 = AnaXS_Aivazis(0,Eel,x,y) ! nu n
          endif

          write(224,'(2F13.5,1P,4E13.5,5e13.5)') x,y,Q2,nu,Eprime,cost,&
               & 0d0, sigma1,sigma2
          write(225,'(2F13.5,1P,4E13.5,5e13.5)') x,y,Q2,nu,Eprime,cost,&
               & 0d0, sigma1/fak,sigma2/fak

       enddo
       write(124,*)
       write(125,*)
       write(224,*)
       write(225,*)
    end do


  end subroutine DoAnaEstimate2

  !====================================================

  subroutine WriteStatusMC(iEV, nEV, nEV_100, nEv_10)
    IMPLICIT NONE
    integer iEV, nEV, nEV_100, nEv_10


    integer n_Perc

    if (nEV.lt.100) return  ! dont write status
    if (mod(iEV,nEV_100) .eq. 0) then
       n_Perc = iEV/nEV_100

         write(*,*) 'TryNeutrino : ',n_Perc,'%'

       open(12,file='TryNeutrino.run',status='unknown')
       write(12,*) 'TryNeutrino.run : ',n_Perc,'%'
       close(12)

    endif

  end subroutine WriteStatusMC

  subroutine ReadInit
    use output
    IMPLICIT NONE


    integer :: ios
    NAMELIST /TryNeutrino/ Eel,nEv

    call Write_ReadingInput('TryNeutrino',0)
    rewind(5)
    read(5,nml=TryNeutrino,iostat=ios)
    call Write_ReadingInput('TryNeutrino',0,ios)

    if (Eel<0.0) then
       Eel = -Eel
       write(*,*) ' Eel = ',Eel,' (default)'
    else
       write(*,*) ' Eel = ',Eel
    endif
    if (nEv<0) then
       nEv = -nEv
       write(*,*) ' nEv = ',nEv,'   = 10^',log10(1d0*nEv),' (default)'
    else
       write(*,*) ' nEv = ',nEv,'   = 10^',log10(1d0*nEv)
    endif

    call Write_ReadingInput('TryNeutrino',1)

    srts2 = 0.938**2+2*0.938*Eel

  end subroutine ReadInit


end program TryNeutrino
