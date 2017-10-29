!******************************************************************************
!****m* /pionPot
! NAME
! module pionPot
! PURPOSE
! Includes several functions to determine the pion potential.
!******************************************************************************
module pionPot

  implicit none
  private

  public :: pionPot_Main

contains


  !****************************************************************************
  !****f* pionPot/pionPot_Main
  ! NAME
  ! real function pionPot_Main (masse, momentum, rhoProton, rhoNeutron, charge, potentialSwitch)
  ! PURPOSE
  ! Determine the scalar potential for pions (used in mesonpotential).
  ! Potential is defined as 0th component of a vector potential in the LRF.
  ! The potential should not be used above a momentum of 250 MeV.
  !****************************************************************************
  real function pionPot_Main (masse, momentum, rhoProton, rhoNeutron, charge, potentialSwitch)

    ! input
    real,    intent(in)                 :: masse
    real,    intent(in), dimension(1:3) :: momentum
    real,    intent(in)                 :: rhoProton, rhoNeutron
    integer, intent(in)                 :: charge, potentialSwitch

    ! local
    real :: rho, absMomentum
    real, parameter :: maxMom = 0.25   ! maximum momentum in GeV

    ! initialize
    rho=rhoProton+rhoNeutron
    absMomentum=SQRT(momentum(1)**2+momentum(2)**2+momentum(3)**2)


    select case (potentialSwitch)
    case (1)
       if (absMomentum>maxMom) then
          !        write(*,*) 'Problem in file pionPot_Main.f90 subroutine pionPot_Main'
          !          write(*,*) 'Problem with pionPotential. Momentum too high:' ,absMomentum
          !          write(*,*) 'If you want to propagate pions with such high momentum, then you need to switch off their potential!'
          !          write(*,*) 'Critical error! STOP'
          pionPot_Main=0.
          !
          !          stop
       else
         ! potential in Hamiltonian H=SQRT(p**2+m**2)+V(Nuk)+V(Coulomb)
         pionPot_Main = pionPot_Oset (momentum, rhoProton, rhoNeutron, charge, .true.)
       end if
    case (2)
       ! potential in Hamiltonian H=SQRT(p**2+m**2)+V(Nuk)+V(Coulomb)
       pionPot_Main = pionPot_Kapusta (absMomentum, rho)
    case (3)
       if (absMomentum>maxMom) then
          !          write(*,*) 'Problem in file pionPot_Main.f90 subroutine pionPot_Main'
          !            write(*,*) 'Problem with pionPotential. Momentum too high:' ,absMomentum
          !            write(*,*) 'If you want to propagate pions with such high momentum, then you need to switch off their potential!'
          write(*,*) 'Critical error! STOP'
          stop
       else
         pionPot_Main = pionPot_DeltaHole (momentum, rhoProton, rhoNeutron, masse)
       end if
    case (4)
       if (absMomentum>maxMom) then
          !            write(*,*) 'Problem in file pionPot_Main.f90 subroutine pionPot_Main'
          !              write(*,*) 'Problem with pionPotential. Momentum too high:' ,absMomentum
          !              write(*,*) 'If you want to propagate pions with such high momentum, then you need to switch off their potential!'
          !          pionPot_Main=0.
          pionPot_Main = pionPot_DeltaHole ((/maxMom,0.,0./), rhoProton, rhoNeutron, masse)
          !
          !              write(*,*) 'Critical error! STOP'
          !             stop
       else
          pionPot_Main = pionPot_Smooth (momentum, rhoProton, rhoNeutron, masse, charge)
       end if
    case default
       pionPot_Main=0.
    end select

  end function pionpot_Main


  !****************************************************************************
  !****f* pionPot/pionPot_DeltaHole
  ! NAME
  ! function pionPot_DeltaHole (p, rhop, rhon, masse) result (Vopt)
  ! PURPOSE
  ! Evaluates potential for pion according to simple Delta-Hole-Model
  ! See Diploma-Thesis O.Buss, pages 18-22.
  ! INPUTS
  ! * real   ::  rhop,rhon -- proton and neutron density
  ! * real,dimension(1:3), intent(in)    ::  p-- !momentum of pion in GEV
  ! * real, intent(in)    ::  masse -- mass of pion in GEV
  ! RESULT
  ! * real, intent(out)   :: Vopt -- in GEV
  !****************************************************************************
  function pionPot_DeltaHole (p, rhop, rhon, masse) result (Vopt)
    use constants, only: mN, hbarc

    ! input
    real, dimension(1:3), intent(in) :: p
    real, intent(in) :: rhop, rhon, masse
    ! output
    real :: Vopt

    ! local
    double precision :: qtot, rho, maPion
    double precision, parameter :: mDelta=1.232
    double precision, parameter :: fDelta=2.0
    double precision, parameter :: gprime=0.62

    qtot=DBLE(sqrt(p(1)**2+p(2)**2+p(3)**2))
    rho=DBLE(rhop+rhon)
    maPion=Dble(masse)

    Vopt=Real(Zpion(qtot,rho)*pion(qtot,rho)+Zloch(qtot,rho)*loch(qtot,rho))-SQRT(masse**2+qtot**2)

  contains

    Function pion(q,rho) Result (Ergebnis)   !Ergebnis in GEV
      double precision q,rho
      double precision Epion,EnukDel,Edelta,RhoGEV
      double precision Ergebnis,C
      rhoGEV=rho*hbarc**3
      Epion=SQRT(mapion**2+q**2)
      Edelta=SQRT(mDelta**2+q**2)-mN
      C=8./9.*rhoGEV*(fDelta/mapion)**2 !units: GEV
      ENukDel=SQRT(Edelta*(Edelta+gprime*C))
      Ergebnis=SQRT(0.5*(ENukDel**2+Epion**2-SQRT((ENukDel**2-Epion**2)**2+4*q**2*C*Edelta)))
    end function pion

    Function loch(q,rho) Result (Ergebnis) !Ergebnis in GEV
      double precision q,rho
      double precision Epion,EnukDel,Edelta,RhoGEV
      double precision Ergebnis,C
      rhoGEV=rho*hbarc**3
      Epion=SQRT(mapion**2+q**2)
      Edelta=SQRT(mDelta**2+q**2)-mN
      C=8./9.*rhoGEV*(fDelta/mapion)**2 !units: GEV
      ENukDel=SQRT(Edelta*(Edelta+gprime*C))
      Ergebnis=SQRT(0.5*(ENukDel**2+Epion**2+SQRT((ENukDel**2-Epion**2)**2+4*q**2*C*Edelta)))
    end function loch

    Function Zpion(q,rho) Result (Ergebnis)   !Ergebnis einheitenlos
      double precision q,rho
      double precision Epion,EnukDel,Edelta,RhoGEV
      double precision Ergebnis,C
      rhoGEV=rho*hbarc**3
      Epion=SQRT(mapion**2+q**2)
      Edelta=SQRT(mDelta**2+q**2)-mN
      C=8./9.*rhoGEV*(fDelta/mapion)**2 !units: GEV
      ENukDel=SQRT(Edelta*(Edelta+gprime*C))
      Ergebnis=1./2.*(1+(ENukDel**2-Epion**2)/SQRT((ENukDel**2-Epion**2)**2+4*q**2*C*Edelta))
    end function Zpion

    Function Zloch(q,rho) Result (Ergebnis) !Ergebnis einheitenlos
      double precision q,rho
      double precision Epion,EnukDel,Edelta,RhoGEV
      double precision Ergebnis,C
      rhoGEV=rho*hbarc**3
      Epion=SQRT(mapion**2+q**2)
      Edelta=SQRT(mDelta**2+q**2)-mN
      C=8./9.*rhoGEV*(fDelta/mapion)**2 !units: GEV
      ENukDel=SQRT(Edelta*(Edelta+gprime*C))
      Ergebnis=1./2.*(1-(ENukDel**2-Epion**2)/SQRT((ENukDel**2-Epion**2)**2+4*q**2*C*Edelta))
    end function Zloch

  end function pionPot_DeltaHole


  !****************************************************************************
  !****f* pionPot/pionPot_Kapusta
  ! NAME
  ! function pionPot_Kapusta (k, rho) result (VOpt)
  ! PURPOSE
  ! Evaluates potential for pion according to Kapusta.
  ! INPUTS
  ! * real   ::  rho -- density in fm^-3
  ! * real, intent(in) :: k -- absolute momentum in GeV
  ! RESULT
  ! * real, intent(out):: VOpt -- optical potential as oth component of vector potential in GeV
  !****************************************************************************
  function pionPot_Kapusta (k, rho) result (VOpt)
    use constants, only: mPi, rhoNull

    real, intent(in) :: k, rho
    real :: VOpt

    real, parameter :: alpha=0.154
    real :: omega,omegaNeu,mNull,kNull, x

    omega=SQRT(k**2+mPi**2)
    x=exp(-alpha*rho/rhonull)
    mNull=(1+6.5*(1.-x**10))*mPi
    kNull=SQRT((1.-x)**2+2.*mNull/mPi*(1.-x))*mPi

    omegaNeu=SQRT((abs(k)-kNull)**2+mNull**2)-SQRT(kNull**2+mNull**2)+mPi

    Vopt=omegaNeu-omega
  end function pionPot_Kapusta


  !****************************************************************************
  !****f* pionPot/pionPot_Smooth
  ! NAME
  ! function pionPot_Smooth (momentum, rhop, rhon, masse, charge) result (VOpt)
  ! PURPOSE
  ! Evaluates the pion potential by linear interpolation between Oset potential
  ! (around Ekin=0) and a simple Delta-hole potential. Valid up to Ekin=130MeV.
  ! INPUTS
  ! * real   ::  rhop,rhon !proton and neutron density
  ! * real, dimension(1:3),intent(in) :: momentum  ! in GeV
  ! * real, intent(in)    ::  masse                ! mass of pion in GEV
  ! RESULT
  ! * real, intent(out)   :: Vopt                  ! potential in GEV
  !****************************************************************************
  function pionPot_Smooth (momentum, rhop, rhon, masse, charge) result (VOpt)
    real, dimension(1:3), intent(in) :: momentum
    real, intent(in) :: rhop, rhon, masse
    integer, intent(in) :: charge

    real :: Vopt

    real, parameter :: oben=0.14
    real, parameter :: unten=0.08
    real :: qtot

    logical:: DeltaFlag
    real, dimension(1:3) :: mom

    real :: u,o,uPrime,oPrime!,uDo,oDo ! Variables for spline

    real, parameter :: deltaP=0.005 ! Fuer Ableitung
    integer :: k
    real, dimension(0:3) :: wert

    !**************************************************************************

    DeltaFlag=.true.
    qtot=sqrt(momentum(1)**2+momentum(2)**2+momentum(3)**2)

    if (rhop+rhon<0.001) then
       Vopt=0
       return
    end if

    if (qtot<unten) then
       Vopt = pionPot_Oset (momentum, rhop, rhon, charge, Deltaflag)
    else if (qtot>oben) then
       Vopt = pionPot_DeltaHole (momentum, rhop, rhon, masse)
    else
       !Werte an unterer Schranke
       do k=0,3
          mom(1)=unten-float(k)*DeltaP
          mom(2:3)=0.
          wert(k) = pionPot_Oset (mom, rhop, rhon, charge, Deltaflag)
       end do
       u=wert(0)
       uPrime=(3*wert(0)-4*wert(1)+wert(2))/2./DeltaP
       !uDo=(2*wert(0)-5*wert(1)+4*wert(2)-wert(3))/(DeltaP**2)

       !Werte an oberer Schranke
       do k=0,3
          mom(1)=oben+float(k)*DeltaP
          mom(2:3)=0
          wert(k) = pionPot_DeltaHole(mom,rhop,rhon,masse)
       end do
       o=wert(0)
       oPrime=(-3*wert(0)+4*wert(1)-wert(2))/2./DeltaP
       !oDo=(2*wert(0)-5*wert(1)+4*wert(2)-wert(3))/(DeltaP**2)
       !Vopt=splineFunc2(u,o,uPrime,oPrime,uDo,oDo,qtot)
       Vopt=splineFunc(u,o,uPrime,oPrime,qtot)
       !Vopt=u+(o-u)*(qtot-unten)/(oben-unten)
    end if

  Contains

    Function splineFunc(u,o,uPrime,oPrime,x) Result (Ergebnis)
      Real :: u,o,uPrime,oPrime,x,Ergebnis
      Real :: a,b,c,d
      !Ableitung stetig
      a=(-2*o+2*u+(OPrime+uPrime)*(oben-unten))
      a=a/((oben-unten)**3)
      b=(3*o-3*u)*(oben+unten)-(oben-unten)*((oPrime+2*uPrime)*oben+(2*oPrime+uPrime)*unten)
      b=b/((oben-unten)**3)
      c=uPrime*oben*(oben**2+oben*unten-2*unten**2)-unten*(6*o*oben  &
           &         -6*u*oben+oPrime*(unten**2+oben*unten-2*oben**2))
      c=c/((oben-unten)**3)
      d=u*oben**2*(oben-3*unten)+unten*(uPrime*oben**2*(unten-oben)  &
           &         +unten*(o*(3*oben-unten)+oPrime*oben*(unten-oben)))
      d=d/((oben-unten)**3)
      Ergebnis=a*(x**3)+b*(x**2)+c*x+d
    end function splineFunc

!     Function splineFunc2(u,o,uPrime,oPrime,uDo,oDo,x) Result(Ergebnis)
!       Real :: u,o,uPrime,oPrime,uDo,oDo,x,Ergebnis
!       Real :: a,b,c,d,e,f
!       Real :: obenS  !oben reSkaliert
!
!       obenS=oben-unten  !Reskalieren auf unten=0
!       x=x-unten
!
!       !Ableitung differenzierbar
!       a=12.*o-12.*u+obenS*((oDo-uDo)*obenS-6.*(oPrime+uPrime))
!       a=a/2./(obenS**5)
!       b=30.*obenS*(u-o)+obenS*(14*oPrime*obenS+16*uPrime*obenS      &
!            &                           -(obenS**2)*(2*oDo-3*uDo))
!       b=b/2./(obenS**5)
!       c=20*(obenS**2)*(u-o)              &
!            &      +(obenS**3)*(8*oPrime+12*uPrime-oDO*obenS+3*uDo*obenS)
!       c=-c/2./(obenS**5)
!       d=uDo/2.
!       e=uPrime
!       f=u
!       Ergebnis=a*(x**5)+b*(x**4)+c*(x**3)+d*(x**2)+e*x+f
!     end function splineFunc2

  end function pionPot_Smooth


  !****************************************************************************
  !****f* pionPot/pionPot_Oset
  ! NAME
  ! function pionPot_Oset (momentum, rhop, rhon, charge, Deltaflag, Decay) result (VOpt)
  ! PURPOSE
  ! Evaluates real part of the pion selfenergy and pion decay width
  ! according to Oset, Salcedo et al, NPA 554 (pages 554 and following).
  ! All calculations are done in MeV and then converted to GeV in the end!
  ! NOTES
  ! This routine calls two subroutines, which evaluate s-Wave and P-Wave of the
  ! selfenergy. These Subroutines evaluate the full selfenergy.
  !****************************************************************************
  function pionPot_Oset (momentum, rhop, rhon, charge, Deltaflag, Decay) result (VOpt)
    use constants, only: pi

    !input
    real, dimension(1:3),intent(in) :: momentum  !Dreier-Impuls
    real, intent(in)                :: rhop,rhon !Proton-,Neutrondichte
    integer, intent(in)             :: charge    !Ladung
    logical, intent(in)             :: DeltaFlag !.true.=decay-constant with direct delta terms
                                                 !.false.=decay-constant without direct delta terms
    ! output
    real, intent(out), optional :: Decay    ! Decay width
    real                        :: VOpt     ! Potential

    !local
    real, parameter :: maPion=138.
    real, parameter :: mNukleon=938.5
    real, parameter :: mDelta=1232.
    real, parameter :: mRes=mDelta-mNukleon
    real, parameter :: fmMEv=197.327053
    real, parameter :: rhoNull=0.17*fmMEv**3
    real            :: qcm,epsilon,E,V!,DeltaE         !kinematics
    real            :: kF,s                            !kinematics
    real            :: ImSwave, ImPwave                !absorbtive Part of Potential
    complex         :: PwavePart,SwavePart             !P and S wave of potential
    real            :: pD,nD                           !Proton and neutron Density
    real            :: dummy,Edummy
    real            :: T                               !kinematics
    logical         :: bigFlag,toobigflag
    real,parameter  :: oben=10000.0   !hier ist Re[potential] null gesetzt("Kuenstlich"!!)
    real,parameter  :: unten=10000.0  !hier ist obere Gueltigkeitsgrenze des potentials

    E=SQRT((momentum(1)**2+momentum(2)**2+momentum(3)**2)*1000.**2+maPion**2)
    !Energy of pion in MEV, 1000**2 is MEV<->GEV conversion

    bigflag=.false.
    toobigflag=.false.

    if (E>unten) then
       Edummy=E
       E=unten
       bigFlag=.true.
       if (Edummy>oben) then
          toobigflag=.true.
       end if
    end if
    pD=rhop*fmMev**3         !Proton Density in MEV**3
    nD=rhon*fmMev**3         !Neutron Density in MEV**3

    V=0.                               !Don't consider Coulomb Potential to evaluate strong potential
    epsilon=E/mNukleon
    kF=(3./2.*pi**2.*(pd+nd))**(1./3.)                              !Fermi momentum
    s=mNukleon**2+maPion**2+2*E*(mNukleon+3./5.*kF**2/2./mNukleon)  ! meanvalue of s

    ! Meanvalue of Center of Mass Momentum:
    qcm=((((s-mNukleon**2-maPion**2)/2.)**2-mNukleon**2*maPion**2)/s)**(1./2.)

    T=(E-maPion-V)/maPion

    !Call routines, which evaluate s- and p-wave

    select case (charge)
    case (-1)  ! pi-
       call swave(ImSwave,SwavePart)     !Evaluates S-wave
       call pwave(ImPwave,PwavePart)     !Evaluates P-Wave
       VOpt=Real(PwavePart+SwavePart)/2./E/1000.
       if (present(Decay)) Decay=(ImPwave+ImSwave)/E/(-1000.)
       !VOpt=Real(PwavePart)/2./E/1000.
       !Decay=(ImPwave)/E/(-1000.)

    case (1)  ! pi+
       !Rotating Isospin, because Oset's paper is made for PiMinus
       dummy=pD
       pD=nD
       nD=dummy
       call swave(ImSwave,SwavePart)     !Evaluates S-wave
       call pwave(ImPwave,PwavePart)     !Evaluates P-Wave
       VOpt=Real(PwavePart+SwavePart)/2./E/1000.
       if (present(Decay)) Decay=(ImPwave+ImSwave)/E/(-1000.)

    case (0)  ! pi0
       call swave(ImSwave,SwavePart)     !Evaluates S-wave
       call pwave(ImPwave,PwavePart)     !Evaluates P-Wave
       VOpt=Real(PwavePart+SwavePart)
       if (present(Decay)) Decay=(ImPwave+ImSwave)/E

       !Rotating Isospin, because Oset's paper is made for PiMinus
       dummy=pD
       pD=nD
       nD=dummy
       call swave(ImSwave,SwavePart)     !Evaluates S-wave
       call pwave(ImPwave,PwavePart)     !Evaluates P-Wave
       VOpt=(VOpt+Real(PwavePart+SwavePart))/2./E/1000./2.
       if (present(Decay)) Decay=((Decay+(ImPwave+ImSwave)/E)/2.)/(-1000.)

    case default ! ERROR ERROR
       write(*,*) "Problems with charges in oset.f"
       stop
    end select

    ! RE[Pot] stetig gegen Null gehen lassen an Grenze des Gueltigkeitsbereichs
    if (bigFlag) then
       Vopt=Vopt*(1.-(Edummy-unten)/(oben-unten))
    end if
    if (toobigFlag) then
       Vopt=0.
       if (present(Decay)) Decay=0.
    end if


  contains


    !  ** interne Prozeduren **********************************************

    subroutine swave(ImSwave,SwavePart)

      !************************************************************************
      !Evaluates Swave part of selfenergy of pion
      !The names of the variabels are according to NPA 554, pages 554ff
      !************************************************************************

      !output
      real, intent(out)   :: ImSwave    !Only absorptive part of S-wave
      complex,intent(out) :: SwavePart  !Full S-wave
      !local
      real, parameter       :: ImBO=0.041*maPion**(-4)
      real, parameter       :: bNull=-0.013/maPion
      real, parameter       :: bOne=-0.092/maPion
      real, parameter       :: DeltabNull=-0.0053/maPion
      real, parameter       :: DeltabOne=-0.013/maPion
      real, parameter       :: DeltaImBO=0.0064*maPion**(-4)
      real, dimension(1:14) :: bOneQ,bNullQ
      complex               :: V1,V3
      Real                  :: V2,V1abs, V1Quasi, V1Real
      integer               :: k
      Real                  :: bOneQdummy,bNullQDummy

      V1Abs=-4.*pi*ImBO*(1+epsilon/2.)*2.*(pD**2.+nD*pD) !Absorptive Part of Swave

      V1Quasi=-4.*pi*(1.+epsilon/2.)*(pD+nD)**2*IMBNullQ(T) !Quasielastic part of SWave
      V1Real=-4.*pi*((1.+epsilon)*(bNull-3./2./pi*kF*(1.+epsilon)*  &
           &          (bNull**2+2.*bOne**2))*(pD+nD)*(1.+T*maPion/100.)  &
           &           +(1.+epsilon)*bOne*(nD-pD))
      V1=CMPLX(V1Real,V1Abs+V1Quasi)
      ! Evaluate V2

      bNullQ(1)= 6.0
      bNullQ(2)= 13.0
      bNullQ(3)= 21.0
      bNullQ(4)= 30.0
      bNullQ(5)= 39.0
      bNullQ(6)= 50.0
      bNullQ(7)= 62.0
      bNullQ(8)= 75.0
      bNullQ(9)= 88.0
      bNullQ(10)= 103.0
      bNullQ(11)= 119.0
      bNullQ(12)= 136.0
      bNullQ(13)= 155.0
      bNullQ(14)= 174.0

      bOneQ(1)= -2.0
      bOneQ(2)= -4.0
      bOneQ(3)= -5.0
      bOneQ(4)= -5.0
      bOneQ(5)= -5.0
      bOneQ(6)= -5.0
      bOneQ(7)= -4.0
      bOneQ(8)= -2.0
      bOneQ(9)= 1.0
      bOneQ(10)= 4.0
      bOneQ(11)= 8.0
      bOneQ(12)= 13.0
      bOneQ(13)= 18.0
      bOneQ(14)= 25.0
      do k=1,14
         bNullQ(k)=bNullQ(k)/maPion/10000.
         bOneQ(k)=bOneQ(k)/maPion/10000.
      end do
!!$    *        k=ANINT(T*maPion/5.)
!!$    *        if(k.le.0) then
!!$    *          k=1
!!$    *        end if
!!$    *        if(k.ge.14) then
!!$    *          k=14
!!$    *        end if
!!$    *        V2=-4.*pi*(bNullQ(k)*(pD+nD)+bOneQ(k)
!!$    *     &     *(nD-pD))

      k=INT(T*maPion/5.)
      if (k.le.0) then
         k=1
      end if
      if (k.gt.13) then
         k=14
         write(*,*) "Energy too big in Oset",k,T
         stop
      end if
      bNullQDummy=bNullQ(k)+(T*maPion/5.-k)*(bNullQ(k+1)-bNullQ(k))
      bOneQDummy=bOneQ(k)+(T*maPion/5.-k)*(bOneQ(k+1)-bOneQ(k))
      V2=-4.*pi*(bNullQdummy*(pD+nD)+bOneQdummy*(nD-pD))



      ! Evaluate V3
      V3=-4.*pi*CMPLX((1.+epsilon)*(DeltaBNull*(pD+nD)                 &
           &       +DeltaBOne*(nD-pD)),DeltaIMBO*(1+epsilon/2.)*2.          &
           &       *(pD**2.+pD*nD))

      ImSwave=V1abs
      SwavePart=V1+CMPLX(0.,V2)+V3

    end subroutine swave
    !**************************************************************************
    !**************************************************************************


    subroutine pwave(ImPwave,PwavePart)

      !************************************************************************
      !Evaluates Pwave part of selfenergy of pion
      !The names of the variabels are according to NPA 554, pages 554ff
      !************************************************************************

      !output
      real, intent(out)   ::ImPwave
      complex, intent(out) ::PwavePart
      !local
      real, parameter :: fStarSquare=0.36*4.*pi
      real, parameter :: fSquare=0.08*4.*pi
      real, parameter :: ReCoNR=0.498*maPion**(-6)
      real, parameter :: gPrime=0.63
      real            :: x
      real            :: alphaN,alphaC,renormp
      complex         :: SigmaDelP,SigmaDelN!,SigmaDelRho
      complex         :: SigmaQA3Rho,p1,p2,p3,SigmaA3
      complex         :: alphaD, alphaT,alphaDAbs
      complex         :: renormpFull

      ! Nonresonant terms:
      alphaN=1./4./pi*2.*fSquare/(maPion**2.)*((1.+epsilon/2.)**2.)/E*(nD-pD)

      ! Direct Delta/Hole
      x=2.*pD/rhoNull
      SigmaDelP=CMPLX(-53.*x,ImSigmaDelta(x,T))
      x=2.*nD/rhoNull
      SigmaDelN=CMPLX(-53.*x,ImSigmaDelta(x,T))
      x=(pD+nD)/rhoNull
      !SigmaDelRho=CMPLX(-53.*x,ImSigmaDelta(x,T))
      SigmaQA3Rho=CMPLX(0,ImQA3Delta(x,T))
      SigmaA3=CMPLX(0,-ca3(T)*x**alphaa3(T))

      alphaD=CMPLX(nD,0.)*(CMPLX(SQRT(s)                                         &
           &          -mDelta,0.)-SigmaDelP+CMPLX(0.,GammaDelta(qcm,kF,s)/2.)         &
           &          -SigmaQA3Rho)**(-1.)                                            &
           &          +CMPLX(pD/3.,0.)                                                &
           &          *(CMPLX(SQRT(s)-mDelta,0.)-2./3.*SigmaDelP                      &
           &            -1./3.*SigmaDelN+CMPLX(0.,GammaDelta(qcm,kF,s)                &
           &            /2.)-SigmaQA3Rho)**(-1.)

      ! Crossed Delta/Hole
      alphaC=pD*(-SQRT(s)+mNukleon-mRes+53.*2.*nD/rhoNull)**(-1.)                &
           &     +nD/3.*(-SQRT(s)+mNukleon-mRes+2./3.*53.*2.*nD/rhoNull+1./3.         &
           &     *53.*2.*pD/rhoNull)**(-1.)

      ! Lot of second order stuff
      alphaT=CMPLX(ReCoNR*(pD+nD)**2.,ImCONR(T)*rhoNull                          &
           &          /ImSigmaDeltaHat(T)/anor(T)                                     &
           &          *(app(T)*pD*AIMAG(SigmaDelP)+anp(T)*nD*AIMAG(SigmaDelP)         &
           &          +apn(T)*pD*AIMAG(SigmaDelN)))

      P1=(alphaD+CMPLX(alphaC,0))*                                               &
           &         (-2./3./4./pi*fStarSquare/maPion**2)                             &
           &         +CMPLX(alphaN,0)+alphaT

      alphaDAbs=CMPLX(nD,0.)*(abs(CMPLX(SQRT(s)                                &
           &          -mDelta,0.)-SigmaDelP+CMPLX(0.,GammaDelta(qcm,kF,s)/2.)       &
           &          -SigmaQA3Rho)**2)**(-1.)*                                     &
           &          (-CMPLX(0.,GammaDelta(qcm,kF,s)                               &
           &            /2.)+SigmaA3+SigmaDelP)                                     &
           &          +CMPLX(pD/3.,0.)                                              &
           &          *(abs(CMPLX(SQRT(s)-mDelta,0.)-2./3.*SigmaDelP                &
           &            -1./3.*SigmaDelN+CMPLX(0.,GammaDelta(qcm,kF,s)              &
           &            /2.)-SigmaQA3Rho)**2)**(-1.)*                               &
           &          (-CMPLX(0.,GammaDelta(qcm,kF,s)                               &
           &            /2.)+SigmaA3+2./3.*SigmaDelP                                &
           &            +1./3.*SigmaDelN)


      Call PwaveRest(p2,p3)

      if (DeltaFlag) then
         renormP=AIMAG(alphaT-2./3./4./pi*fStarSquare/maPion**2                &
              &             *alphaDAbs)/(1+16*pi**2*gPrime**2*(abs(p1+p2+p3))**2    &
              &           +4.*pi*gPrime*2.*Real(p1+p2+p3)) !Absorptive part
      else
         renormP=AIMAG(alphaT)/(1+16*pi**2*gPrime**2                               &
              &           *(abs(p1+p2+p3))**2                                           &
              &           +4.*pi*gPrime*2.*Real(p1+p2+p3)) !nonresonant absorbtive part
      end if

      renormpFull=(P1+P2+P3)/(CMPLX(1.,0)+4.*pi*gPrime*(p1+p2+p3)) !full expression

      ImPwave=-4.*pi*mNukleon**2./s*(1.-epsilon/2.)*renormP*(E**2-mapion**2) !Absorptive part of Pwave
      PWavePart=-4.*pi*mNukleon**2./s*(1.-epsilon/2.)*renormpFull*(E**2-mapion**2.) !Full P-Wave

    end subroutine pwave


    !**************************************************************************

    subroutine PwaveRest(p2,p3)

      !************************************************************************
      !Evaluates quasielastic part of the pwave and real part of pwave
      !The names of the variabels are according to NPA 554, pages 554ff
      !************************************************************************

      !output
      Complex, intent(out)   :: p2,p3

      !local
      Real,parameter :: DeltacNull=0.03*maPion**(-3)
      Real,parameter :: DeltacOne=0.1*maPion**(-3)
      Real,parameter :: DeltaImCNull=0.107*maPion**(-6)
      Real           :: ImSigmaDelP, ImSigmaDelN,x
      Real,parameter :: deltaT=5.
      Real           :: cNullQ(14),cOneQ(14)
      Integer        :: i,k
      Real           :: cOneQdummy,cNullQDummy

      cNullQ(1)= 1.0
      cNullQ(2)= 3.0
      cNullQ(3)= 8.0
      cNullQ(4)= 14.0
      cNullQ(5)= 23.0
      cNullQ(6)= 34.0
      cNullQ(7)= 46.0
      cNullQ(8)= 61.0
      cNullQ(9)= 78.0
      cNullQ(10)= 97.0
      cNullQ(11)= 118.0
      cNullQ(12)= 141.0
      cNullQ(13)= 166.0
      cNullQ(14)= 192.0

      cOneQ(1)= 0.0
      cOneQ(2)= 1.0
      cOneQ(3)= 3.0
      cOneQ(4)= 6.0
      cOneQ(5)= 9.0
      cOneQ(6)= 14.0
      cOneQ(7)= 20.0
      cOneQ(8)= 27.0
      cOneQ(9)= 36.0
      cOneQ(10)= 46.0
      cOneQ(11)= 57.0
      cOneQ(12)= 69.0
      cOneQ(13)= 81.0
      cOneQ(14)= 94.0

      do k=1,14
         cNullQ(k)=CNullQ(k)/(maPion**3.)/10000.
         cOneQ(k)=cOneQ(k)/(maPion**3.)/10000.
      end do


!!$    *      i=ANINT(T*maPion/deltaT)
!!$
!!$    *      if(i.le.0) then
!!$    *         i=1
!!$    *      end if
!!$    *      if(i.ge.14) then
!!$    *         i=14
!!$    *      end if

      i=INT(T*maPion/deltaT)

      if (i.le.0) then
         i=1
      end if
      if (i.gt.13) then
         write(*,*) "Energy too big in oset",i,T
         stop
      end if

      x=2.*pD/rhoNull
      ImSigmaDelP=ImSigmaDelta(x,T)
      x=2.*nD/rhoNull
      ImSigmaDelN=ImSigmaDelta(x,T)

      cNullQDummy=cNullQ(i)+(T*maPion/deltaT-i)*(cNullQ(i+1)-cNullQ(i))
      cOneQDummy=cOneQ(i)+(T*maPion/deltaT-i)*(cOneQ(i+1)-cOneQ(i))

      p2=CMPLX(0.,cNullQDummy*(pD+nD)+cOneQDummy*(nD-pD))

      !    *      p2=CMPLX(0.,cNullQ(i)*(pD+nD)+cOneQ(i)*(nD-pD))
      p3=CMPLX((1.+epsilon)*(DeltaCNull*(pD+nD)+DeltaCOne*(nD-pD))     &
           &   ,DeltaImCNull*rhoNull/IMSigmaDeltaHat(T)*(pD*ImSigmaDelN     &
           &    +nD*ImSigmaDelP))


    end subroutine PWaveRest


    ! *******************intern function declarations****************************************************
    !  * Functions are named according to NPA 554, pages 554ff****************************************

    Function ImBNullQ(T) Result(Ergebnis)
      Real T, Ergebnis
      Ergebnis=T*(-14.+190.8*T-182.3*T**2.)*maPion**(-4.)/1000.
    end function ImBNullQ

    Function ImCONR(T) Result(Ergebnis)
      Real T, Ergebnis
      Ergebnis=(0.36-0.497*T+0.52*T**2.)*maPion**(-6.)
    end function ImCONR

    Function anor(T) RESULT (Ergebnis)
      Real T,Ergebnis
      Ergebnis=1.3-3.22*T+2.59*T**2.
    end function anor

    Function anp(T) RESULT (Ergebnis)
      Real T,Ergebnis
      Ergebnis=1.735-4.07*T+3.11*T**2.
    end function anp

    Function apn(T) RESULT (Ergebnis)
      Real T,Ergebnis
      Ergebnis=0.295-0.75*T+0.613*T**2.
    end function apn

    Function app(T) RESULT (Ergebnis)
      Real T,Ergebnis
      Ergebnis=0.57-1.617*T+1.436*T**2.
    end function app

    Function ImSigmaDeltaHat(T) RESULT (Ergebnis)
      Real T,Ergebnis
      Ergebnis=-38.3*(1.-0.85*T+0.54*T**2.)
    end function ImSigmaDeltaHat

    Function a(T) RESULT (Ergebnis)
      Real T,Ergebnis
      Ergebnis=2.72-4.07*T+3.07*T**2
    end function a
    ! *******************************************
    Function cq(T) RESULT (Ergebnis)
      Real T,Ergebnis
      Ergebnis=T*(20.2-8.58*T+0.702*T**2.)
    end function cq

    Function alphaq (T) RESULT (Ergebnis)
      Real T,Ergebnis
      Ergebnis=1.+T*(-0.309-0.315*T+0.151*T**2.)
    end function alphaq
    !**************************************************************************
    Function ca3(T) RESULT (Ergebnis)
      Real T,Ergebnis
      if (T.ge.85/maPion) then
         Ergebnis=T*(-7.08+27.4*T-9.49*T**2)
      else
         Ergebnis=T*maPion*3.7/85.
      end if
    end function ca3

    Function alphaa3 (T) RESULT (Ergebnis)
      Real T,Ergebnis
      Ergebnis=1.+T*(0.984-0.512*T+0.1*T**2)
    end function alphaa3
    !**************************************************************************
    Function Ione(qtilte) Result (Ergebnis)
      Real qtilte, Ergebnis
      if (qtilte.ge.1) then
         Ergebnis=1.-1./5./(qTilte**2)
      else
         Ergebnis=1.+qtilte-1.-(qtilte**3.)/5.
      end if
    end function Ione

    Function Itwo(qtilte) Result (Ergebnis)
      Real qtilte, Ergebnis
      if (qtilte.ge.1) then
         Ergebnis=1-3./5./(qtilte**2)-4./21./(qtilte**6)+18./35./(qtilte**4)
      else
         Ergebnis=33./35.*qtilte-23./105.*qtilte**3
      end if
    end function Itwo

    !**************************************************************************

    Function ImSigmaDelta(x,T) Result (Ergebnis)
      Real x,T,Ergebnis
      Ergebnis=ImSigmaDeltaHat(T)/a(T)*ATAN(x*a(T))
    end function ImSigmaDelta

    Function ImQA3Delta(x,T) Result (Ergebnis)
      Real x,T,Ergebnis
      Ergebnis=-cq(T)*x**alphaq(T)-ca3(T)*x**alphaa3(T)
    end function ImQA3Delta

    Function GammaDelta(qcm,kF,s) Result (Ergebnis)
      Real qcm,s,Ergebnis,qtilte,kF,GamFree
      real, parameter :: fStarSquare=0.36*4.*pi

      qtilte=qcm/kF
      GamFree=2./3.*fStarSquare/(maPion**2.)/4./pi*qcm**3*mNukleon/SQRT(s)
      Ergebnis=(Ione(qtilte)+Itwo(qtilte))/2.*GamFree
    end function GammaDelta

    !**************************************************************************

  end function pionPot_Oset


end module pionPot
