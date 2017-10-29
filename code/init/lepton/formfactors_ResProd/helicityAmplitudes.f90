!******************************************************************************
!****m* /helicityAmplitudes
! NAME
! module helicityAmplitudes
! PURPOSE
! * Provides the helicity amplitudes of the MAID analysis
!   -> http://www.kph.uni-mainz.de/MAID/maid2003/maid2003.html
! * Modified version of a source code provided by L. Tiator
! * Citation: D. Drechsel,  S.S. Kamalov, L. Tiator; Nucl. Phys. A645 (1999) 145-174
!******************************************************************************
module helicityAmplitudes
  private

  public :: get_helicityAmplitudes,writeOutFORMFACTOR
  public :: HP33_MAID05  ! Delta(1232) form factor

  real, parameter :: X3P33=1
  real, parameter :: X1P33=1
  real, parameter :: XSP33=1
  real, parameter :: X1S31=1
  real, parameter :: XSS31=1
  real, parameter :: X1D33=1
  real, parameter :: X3D33=1
  real, parameter :: XSD33=1
  real, parameter :: X1P11p=1
  real, parameter :: XSP11p=1
  real, parameter :: X1S11p=1
  real, parameter :: XSS11p=1
  real, parameter :: X1P11n=1
  real, parameter :: XSP11n=1
  real, parameter :: X1S11n=1
  real, parameter :: XSS11n=1
  real, parameter :: X1D13p=1
  real, parameter :: X3D13p=1
  real, parameter :: XSD13p=1
  real, parameter :: X1F15p=1
  real, parameter :: X3F15p=1
  real, parameter :: XSF15p=1
  real, parameter :: X1D13n=1
  real, parameter :: X3D13n=1
  real, parameter :: XSD13n=1
  real, parameter :: X1F15n=1
  real, parameter :: X3F15n=1
  real, parameter :: XSF15n=1
  real, parameter :: X1S2p=1
  real, parameter :: XSS2p=1
  real, parameter :: X1S2n=1
  real, parameter :: XSS2n=1
  real, parameter :: X1P13p=1
  real, parameter :: X3P13p=1
  real, parameter :: XSP13p=1
  real, parameter :: X1P31=1
  real, parameter :: XSP31=1
  real, parameter :: X1P13n=1
  real, parameter :: X3P13n=1
  real, parameter :: XSP13n=1
  real, parameter :: X1D15p=1
  real, parameter :: X3D15p=1
  real, parameter :: XSD15p=1
  real, parameter :: X1D15n=1
  real, parameter :: X3D15n=1
  real, parameter :: XSD15n=1
  real, parameter :: X1F35=1
  real, parameter :: X3F35=1
  real, parameter :: XSF35=1
  real, parameter :: X1F37=1
  real, parameter :: X3F37=1
  real, parameter :: XSF37=1


contains

  subroutine writeOutFORMFACTOR()
    implicit none
    integer :: ig,iso
    real(8) :: Q2G
    do  IG=1,3
       Q2G=(IG-1)*0.5
       do  ISO=1,2
          CALL HEL_OUT(ISO,Q2G)
       end do
       write(6,*)
    end do
    STOP
  END subroutine writeOutFORMFACTOR



  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  !
  ! MAID 2003
  !
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************


  subroutine HP33(Q2G,qcm,qcm0,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    Fq=exp(-0.21*Q2G)/(1+Q2G/0.71)**2*(qcm/qcm0)
    !c *****************************************************************
    Phi_R=0.
    A10=-140.250*X1P33
    A30=-265.437*X3P33
    S10=  27.577*XSP33

    A1= A10*(1.+ 0.0214486*Q2G)*Fq
    A3= A30*(1.- 0.0065431*Q2G)*Fq
    S1= S10*(1.+ 0.0166667*Q2G**3)*Fq

    AE= (1./2.)*(A3/sqrt(3.) -A1)
    AM=-(1./2.)*(Sqrt(3.0)*A3+A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R,  X1P33, A1,X3P33, A3, XSP33, S1
123 format( 5x,'P33(1232):',1x,F6.2,3(1x,2(F8.3,1x)))
    return
  end subroutine HP33

  subroutine HP11(ISO,PhiR,Q2G,AM,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)

    PI=3.1415926536D0
    Phi_R = -15.48
    PhiR=Phi_R*Pi/180.
    COSR=COS(PhiR)
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    !c ***************************************************************
    A10 = -74.585/COSR*X1P11p
    S10 = 32.995/COSR*XSP11p
    A1= A10*(1.- 0.95*Q2G-0.95*Q2G**4)*exp(-1.55*Q2G)
    S1= S10*(1.+ 7.001*Q2G)*exp(-2.1*Q2G)
    GO TO 20
    !c ***************************************
10  A10 = 49.926/COSR*X1P11n
    S10 = -40./COSR*XSP11n
    A1= A10*(1.+ 0.7*Q2G)*exp(-1.77*Q2G)
    S1=S10*(1.+ 2.97843*Q2G)*exp(-1.55*Q2G)
    !c **************************************************
20  AM=A1
    AS=-sqrt(2.)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3)  write(6,123)  Phi_R, X1P11p, A1, XSP11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4)  write(6,123)  Phi_R, X1P11n, A1, XSP11n, S1
123 format( 5x,'P11(1440):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HP11

  subroutine HD13(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)

    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 32.
    PhiR=Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-25.7796/COSR*X1D13p
    A30=140.439/COSR*X3D13p
    S10=-27.14/COSR*XSD13p

    A1= A10*(1.+6.2193*Q2G)*exp(-1.08576*Q2G)
    A3= A30*(1.+2.2357*Q2G)*exp(-2.54217*Q2G)
    S1= S10*(1.+ 8.429*Q2G)*exp(-3.4534*Q2G)
    GO TO 20
    !c ***************************************************
10  CONTINUE
    Phi_R = 19.
    PhiR=Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-81.02/COSR*X1D13n
    A30=-140.6/COSR*X3D13n
    S10= 12.85/COSR*XSD13n

    A1= A10*(1.-2.31468*Q2G)*exp(-1.55*Q2G)
    A3= A30*(1.+1.0903*Q2G)*exp(-1.75*Q2G)

    S1= S10*(1.+7.13437*Q2G)*exp(-1.57380*Q2G)
    !c ***************************************************
20  AE=-(1./2.)*(sqrt(3.)*A3+A1)
    AM=-(1./2.)*(A3/Sqrt(3.0)-A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3)  write(6,123)  Phi_R, X1D13p, A1, X3D13p, A3, XSD13p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4)  write(6,123)  Phi_R, X1D13n, A1, X3D13n, A3, XSD13n, S1
123 format( 5x,'D13(1520):',1x,F6.2,3(1x,2(F8.3,1x)))
    return
  end subroutine HD13

  subroutine HD33(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)

    PI=3.1415926536D0
    Phi_R = 61.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=65.5267/COSR*X1D33
    A30=103.0903/COSR*X3D33
    S10=1./COSR*XSD33
    !c *************************************************

    A1= A10*(1.+3.83921*Q2G)*exp(-1.77207*Q2G)
    A3= A30*(1.+1.9722*Q2G)*exp(-2.2*Q2G)
    S1= S10/(1.+Q2G/0.71)**2

    AE=-(1./2.)*(A3*sqrt(3.)+A1)
    AM=-(1./2.)*(A3/sqrt(3.)-A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1D33, A1, X3D33, A3, XSD33, S1
123 format( 5x,'D33(1700):',1x,F6.2,3(1x,2(F8.3,1x)))

    Return
  End subroutine HD33

  subroutine HS11f(ISO,PhiR,Q2G,AE,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !    COMMON/par11p/ X1P11p,XSP11p,X1S11p,XSS11p
    !    COMMON/par11n/ X1P11n,XSP11n,X1S11n,XSS11n
    !c *****************************************************************
    PI=3.1415926536D0
    Phi_R=8.193
    PhiR=Phi_R*PI/180.
    COSR=COS(PhiR)

    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10

    A10=72.482/COSR*X1S11p
    S10=-15.70/COSR*XSS11p

    A1=A10*(1.+1.34470*Q2G)*exp(-0.75879*Q2G)
    S1=S10*(1.+2.8261*Q2G)*exp(-0.73735*Q2G)

    GO TO 20
    !c **************************************************

10  A10=-41.94/COSR*X1S11n
    S10=28.18/COSR*XSS11n

    A1=A10*(1.+4.33783*Q2G)*exp(-1.68723*Q2G)
    S1=S10*(1.+0.35874*Q2G)*exp(-1.55*Q2G)

    !c ****************************************************
20  AE=-A1
    AS=-sqrt(2.)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3)  write(6,123)  Phi_R, X1S11p, A1, XSS11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4)  write(6,123) Phi_R,  X1S11n, A1, XSS11n, S1
123 format( 5x,'S11(1535):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HS11f

  subroutine HS11s(ISO,PhiR,Q2G,AE,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !   COMMON/secS11/ X1S11p,XSS11p, X1S11n,XSS11n
    !c *****************************************************************
    PI=3.1415926536D0
    Phi_R=6.961
    PhiR=Phi_R*PI/180.
    COSR=COS(PhiR)

    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10

    A10=32.0913/COSR*X1S11p
    S10=-3.489/COSR*XSS11p

    A1=A10*(1.+1.45359*Q2G)*exp(-0.6167*Q2G)
    S1=S10*(1.+2.878*Q2G)*exp(-0.75879*Q2G)

    GO TO 20
    !c **************************************************

10  A10=26.32/COSR*X1S11n
    S10=10./COSR*XSS11n

    A1=A10*(1.+0.13305*Q2G)*exp(-1.55*Q2G)
    S1=S10*(1.-0.5*Q2G)*exp(-1.55*Q2G)

    !c ********************************************
20  AE=-A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123) Phi_R, X1S11p, A1, XSS11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123) Phi_R,  X1S11n, A1, XSS11n, S1
123 format( 5x,'S11(1650):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HS11s


  subroutine HS31(Q2G,PhiR,AE,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !    COMMON/param3/ X3P33,X1P33,XSP33,X1S31,XSS31,X1D33,X3D33,XSD33
    !c *****************************************************************
    PI=3.1415926536D0
    Phi_R=22.54
    PhiR=Phi_R*PI/180.
    COSR=COS(PhiR)

    A10 = 65.2101/COSR*X1S31
    S10 =1./COSR*XSS31

    A1=A10*(1.+1.5*Q2G)*exp(-3.*Q2G)
    S1=S10/(1.+Q2G/0.71)**2
    AE=-A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R,  X1S31, A1, XSS31, S1
123 format( 5x,'S31(1620):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))
    return
  end subroutine HS31


  subroutine HF15(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !   COMMON/parDFp/ X1D13p,X3D13p,XSD13p,X1F15p,X3F15p,XSF15p
    !  COMMON/parDFn/ X1D13n,X3D13n,XSD13n,X1F15n,X3F15n,XSF15n
    !c ***************************************************************
    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 10.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-24.2131/COSR*X1F15p
    A30=132.2624/COSR*X3F15p
    S10=-25.65/COSR*XSF15p

    A1= A10*(1.+3.72353*Q2G)*exp(-1.2*Q2G)
    A3= A30*(1.+1.53671*Q2G)*exp(-2.22357*Q2G)
    S1= S10*(1.+5.5987*Q2G)*exp(-1.55*Q2G)
    GO TO 20
    !c ***************************************************************
10  CONTINUE
    Phi_R = 15.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=24.49/COSR*X1F15n
    A30=-33.70/COSR*X3F15n
    S10=1./COSR*XSF15n

    A1= A10*(1.+0.*Q2G)*exp(-1.2*Q2G)
    A3= A30*(1.+4.*Q2G)*exp(-1.75*Q2G)
    S1= S10/(1.+Q2G/0.71)**2
    !c *******************************************************
20  AE=-(1./3.)*( A3*sqrt(2.)+A1)
    AM=-(1./3.)*( A3/sqrt(2.)-A1)
    AS=-(1./3.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3)  write(6,123) Phi_R,  X1F15p, A1, X3F15p, A3, XSF15p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4)  write(6,123) Phi_R,  X1F15n, A1, X3F15n, A3, XSF15n, S1
123 format( 5x,'F15(1680):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF15

  subroutine HP31(Q2G,PhiR,AM,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    ! COMMON/parPPp/ X1P13p,X3P13p,XSP13p,X1P31,XSP31

    PI=3.1415926536D0
    Phi_R=35.
    PhiR = Phi_R*PI/180.
    COSR=COS(PhiR)
    A10=34.03/COSR*X1P31
    S10=1./COSR*XSP31

    A1= A10/(1.+Q2G/0.71)**2
    S1= S10/(1.+Q2G/0.71)**2

    AM= A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1P31, A1, XSP31, S1
123 format( 5x,'P31(1910):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    Return
  End subroutine HP31

  subroutine HD15(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !    COMMON/parD15/ X1D15p,X3D15p,XSD15p,X1D15n,X3D15n,XSD15n

    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 20.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)
    A10=21.397/COSR*X1D15p
    A30=22.600/COSR*X3D15p
    S10=1./COSR*XSD15p

    A1=A10/(1.+Q2G/0.71)**2
    A3= A30/(1.+Q2G/0.71)**2
    S1=S10/(1.+Q2G/0.71)**2
    GO TO 20
    !c ********************************************
10  CONTINUE
    Phi_R = 0.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)
    A10=-61.45/COSR*X1D15n
    A30=-74.31/COSR*X3D15n
    S10=-1./COSR*XSD15n

    A1=A10/(1.+Q2G/0.71)**2
    A3= A30/(1.+Q2G/0.71)**2
    S1=S10/(1.+Q2G/0.71)**2


    !c *************************************************
20  AE= (1./3.)*(A3/sqrt(2.)-A1)
    AM=-(1./3.)*(A3*sqrt(2.)+A1)
    AS=-(1./3.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3)  write(6,123) Phi_R,  X1D15p, A1, X3D15p, A3, XSD15p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4)  write(6,123) Phi_R,  X1D15n, A1, X3D15n, A3, XSD15n, S1
123 format( 5x,'D15(1675):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HD15


  subroutine HF35(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !   COMMON/parF3/ X1F35,X3F35,XSF35,X1F37,X3F37,XSF37

    PI=3.1415926536D0
    Phi_R =40.
    PhiR =Phi_R*Pi/180.
    COSR=COS(PhiR)
    A10=11.092/COSR*X1F35
    A30=-16.63/COSR*X3F35
    S10=1./COSR*XSF35

    A1= A10/(1.+Q2G/0.71)**2
    A3= A30/(1.+Q2G/0.71)**2
    S1= S10/(1.+Q2G/0.71)**2

    AE=-(1./3.)*(A3*sqrt(2.)+A1)
    AM=-(1./3.)*(A3/sqrt(2.)-A1)
    AS=-(1./3.)*S1*sqrt(2.)

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1F35, A1, X3F35, A3, XSF35, S1
123 format( 5x,'F35(1905):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF35

  subroutine HF37(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !    COMMON/parF3/ X1F35,X3F35,XSF35,X1F37,X3F37,XSF37

    PI=3.1415926536D0
    Phi_R = 30.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-67.55/COSR*X1F37
    A30=-87.21/COSR*X3F37
    S10=1./COSR*XSF37

    A1=A10/(1.+Q2G/0.71)**2
    A3=A30/(1.+Q2G/0.71)**2
    S1=S10/(1.+Q2G/0.71)**2

    AE= (1./4.)*(A3*3./sqrt(15.)-A1)
    AM=-(1./4.)*(A3*5./sqrt(15.)+A1)
    AS=-(1./4.)*Sqrt(2.)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1F37, A1, X3F37, A3, XSF37, S1
123 format( 5x,'F37(1950):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF37

  subroutine HP13(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !   COMMON/parPPp/ X1P13p,X3P13p,XSP13p,X1P31,XSP31
    !  COMMON/parPPn/ X1P13n,X3P13n,XSP13n

    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 0.
    PhiR = Phi_R*PI/180.
    COSR = COS(PhiR)

    A10=54.89/COSR*X1P13p
    A30=-31.69/COSR*X3P13p
    S10=-53.03/COSR*XSP13p

    A1= A10*(1.+4.24137*Q2G)*EXP(-1.55*Q2G)
    A3= A30*(1.+3.9974*Q2G)*EXP(-1.55*Q2G)
    S1= S10*(1.+2.45819*Q2G)*EXP(-1.55*Q2G)
    GO TO 20
    !c ***************************************************
10  CONTINUE
    Phi_R = 0.
    PhiR = Phi_R*PI/180.
    COSR = COS(PhiR)

    A10=16.68/COSR*X1P13n
    A30=-74.93/COSR*X3P13n
    S10=-1./COSR*XSP13n

    A1= A10*(1.+4.24137*Q2G)*EXP(-1.55*Q2G)
    A3= A30*(1.+0.9974*Q2G)*EXP(-1.55*Q2G)
    S1= S10/(1.+Q2G/0.71)**2

    !c ****************************************************

20  AE= (1./2.)*(A3/sqrt(3.)-A1)
    AM=-(1./2.)*(A3*sqrt(3.)+A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123) Phi_R, X1P13p, A1, X3P13p, A3, XSP13p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123) Phi_R, X1P13n, A1, X3P13n, A3, XSP13n, S1
123 format( 5x,'P13(1720):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HP13

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************




  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  !
  ! MAID 2005
  !
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

  subroutine HP33_MAID05(Q2G,qcm,qcm0,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    ! ***************************************************************
    Fq=exp(-0.21*Q2G)/(1+Q2G/0.71)**2*(qcm/qcm0)
    ! *****************************************************************
    Phi_R=0.
    !C ******  fit parameters changed to XE and XM, 19 April 2005
    !C ***     changes onlyn inside of HP33
    !C ***     outside XMP33 is still calles X3P33 and XEP33 -> X1P33
    XMP33=X3P33
    XEP33=X1P33

    A10=-137.445 *X1P33
    A30=-260.128 *X3P33
    S10=  27.577 *XSP33

    A1= A10*(1.+ 0.0214486*Q2G)*Fq
    A3= A30*(1.- 0.0065431*Q2G)*Fq
    S1= S10*(1.+ 0.0166667*Q2G**3)*Fq

    AE= (1./2.)*(A3/sqrt(3.) -A1)*XEP33
    AM=-(1./2.)*(Sqrt(3.0)*A3+A1)*XMP33
    AS=-(1./2.)*Sqrt(2.0)*S1

    A1=-(1./2.)*(3*AE+AM)
    A3= (1./2.)*Sqrt(3.)*(AE-AM)

    if (IPRN.EQ.0) return
    write(6,123) Phi_R,  X1P33, A1,X3P33, A3, XSP33, S1
123 format( 5x,'P33(1232):',1x,F6.2,3(1x,2(F8.3,1x)))
    return
  end subroutine HP33_MAID05

  subroutine HP11_MAID05(ISO,PhiR,Q2G,AM,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    PI=3.1415926536D0
    Phi_R = -15.48
    PhiR=Phi_R*Pi/180.
    COSR=COS(PhiR)
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    !c ***************************************************************
    A10 = -59.134/COSR  *X1P11p
    S10 = 32.995/COSR   *XSP11p
    A1= A10*(1.- 0.95*Q2G -1.007*Q2G**4)*exp(-1.51*Q2G)
    S1= S10*(1.+ 7.001*Q2G)*exp(-2.1*Q2G)
    GO TO 20
    !c ***************************************
10  A10 = 52.137/COSR *X1P11n
    S10 = -40./COSR   *XSP11n
    A1= A10*(1.+0.9450*Q2G)*exp(-1.77*Q2G)
    S1=S10*(1.+ 2.97843*Q2G)*exp(-1.55*Q2G)
    !c **************************************************
20  AM=A1
    AS=-sqrt(2.)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123)  Phi_R, X1P11p, A1, XSP11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123)  Phi_R, X1P11n, A1, XSP11n, S1
123 format( 5x,'P11(1440):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HP11_MAID05

  subroutine HD13_MAID05(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !       COMMON/parDFp/ X1D13p,X3D13p,XSD13p,X1F15p,X3F15p,XSF15p
    !       COMMON/parDFn/ X1D13n,X3D13n,XSD13n,X1F15n,X3F15n,XSF15n
    !c ***************************************************************
    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    XED13p=1
    XMD13p=1

    Phi_R = 32.
    PhiR=Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-23.2016/COSR *X1D13p
    A30=136.2258/COSR *X3D13p
    S10=-27.14/COSR   *XSD13p

    A1= A10*(1.+6.84123*Q2G)*exp(-1.08576*Q2G)
    A3= A30*(1.+2.2357*Q2G)*exp(-2.54217*Q2G)
    S1= S10*(1.+6.06888*Q2G)*exp(-3.4534*Q2G)
    AE=-(1./2.)*(sqrt(3.)*A3+A1)   * XED13p
    AM=-(1./2.)*(A3/Sqrt(3.0)-A1)  * XMD13p
    AS=-(1./2.)*Sqrt(2.0)*S1
    A1=(3*AM-AE)/2.0
    A3=-Sqrt(3.0)*(AE+AM)/2.0
    GO TO 20
    !c ***************************************************
10  CONTINUE
    Phi_R = 19.
    PhiR=Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-72.362/COSR   *X1D13n
    A30=-145.620/COSR  *X3D13n
    S10= 12.85/COSR    *XSD13n

    A1= A10*(1.-1.20549*Q2G)*exp(-1.55*Q2G)
    A3= A30*(1.+0.24466*Q2G)*exp(-1.75*Q2G)
    S1= S10*(1.+7.13437*Q2G)*exp(-1.57380*Q2G)
    AE=-(1./2.)*(sqrt(3.)*A3+A1)
    AM=-(1./2.)*(A3/Sqrt(3.0)-A1)
    AS=-(1./2.)*Sqrt(2.0)*S1
    !c ***************************************************
20  CONTINUE

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123)  Phi_R, X1D13p, A1, X3D13p, A3, XSD13p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123)  Phi_R, X1D13n, A1, X3D13n, A3, XSD13n, S1
123 format( 5x,'D13(1520):',1x,F6.2,3(1x,2(F8.3,1x)))
    return
  end subroutine HD13_MAID05

  subroutine HD33_MAID05(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !       COMMON/param3/ X3P33,X1P33,XSP33,X1S31,XSS31,X1D33,X3D33,XSD33
    !c*****************************************************************
    PI=3.1415926536D0
    Phi_R = 61.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=109.631/COSR *X1D33
    A30=101.914/COSR *X3D33
    S10=1./COSR *XSD33
    !c *************************************************

    A1= A10*(1.+1.906628*Q2G)*exp(-1.77207*Q2G)
    A3= A30*(1.+1.9722*Q2G)*exp(-2.2*Q2G)
    S1= S10*exp(-2.0*Q2G)

    AE=-(1./2.)*(A3*sqrt(3.)+A1)
    AM=-(1./2.)*(A3/sqrt(3.)-A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1D33, A1, X3D33, A3, XSD33, S1
123 format( 5x,'D33(1700):',1x,F6.2,3(1x,2(F8.3,1x)))

    Return
  End subroutine HD33_MAID05

  subroutine HS11f_MAID05(ISO,PhiR,Q2G,AE,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !COMMON/par11p/ X1P11p,XSP11p,X1S11p,XSS11p
    !COMMON/par11n/ X1P11n,XSP11n,X1S11n,XSS11n
    !c *****************************************************************
    PI=3.1415926536D0
    Phi_R=8.193
    PhiR=Phi_R*PI/180.
    COSR=COS(PhiR)

    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10

    A10=65.751/COSR   *X1S11p
    S10=-15.70/COSR   *XSS11p

    A1=A10*(1.+1.61364*Q2G)*exp(-0.75879*Q2G)
    S1=S10*(1.+2.8261*Q2G)*exp(-0.73735*Q2G)

    GO TO 20
    !c **************************************************

10  A10=-50.148/COSR *X1S11n
    S10=28.18/COSR   *XSS11n

    A1=A10*(1.+2.86297*Q2G)*exp(-1.68723*Q2G)
    S1=S10*(1.+0.35874*Q2G)*exp(-1.55*Q2G)

    !c ****************************************************
20  AE=-A1
    AS=-sqrt(2.)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123)  Phi_R, X1S11p, A1, XSS11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123) Phi_R,  X1S11n, A1, XSS11n, S1
123 format( 5x,'S11(1535):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HS11f_MAID05

  subroutine HS11s_MAID05(ISO,PhiR,Q2G,AE,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !       COMMON/secS11/ X1S11p,XSS11p, X1S11n,XSS11n
    !c *****************************************************************
    PI=3.1415926536D0
    Phi_R=6.961
    PhiR=Phi_R*PI/180.
    COSR=COS(PhiR)

    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10

    A10=33.0210/COSR *X1S11p
    S10=-3.489/COSR  *XSS11p

    A1=A10*(1.+1.45359*Q2G)*exp(-0.6167*Q2G)
    S1=S10*(1.+2.878*Q2G)*exp(-0.75879*Q2G)

    GO TO 20
    !c **************************************************

10  A10=9.186/COSR *X1S11n
    S10=10./COSR   *XSS11n

    A1=A10*(1.+0.13305*Q2G)*exp(-1.55*Q2G)
    S1=S10*(1.-0.5*Q2G)*exp(-1.55*Q2G)

    !c ********************************************
20  AE=-A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123) Phi_R, X1S11p, A1, XSS11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123) Phi_R,  X1S11n, A1, XSS11n, S1
123 format( 5x,'S11(1650):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HS11s_MAID05

  subroutine HS31_MAID05(Q2G,PhiR,AE,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !       COMMON/param3/ X3P33,X1P33,XSP33,X1S31,XSS31,X1D33,X3D33,XSD33
    !c *****************************************************************
    PI=3.1415926536D0
    Phi_R=22.54
    PhiR=Phi_R*PI/180.
    COSR=COS(PhiR)

    A10 = 60.6258/COSR*X1S31
    S10 =1./COSR*XSS31

    A1=A10*(1.+1.5*Q2G)*exp(-3.*Q2G)
    S1=S10*exp(-2.*Q2G)

    AE=-A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R,  X1S31, A1, XSS31, S1
123 format( 5x,'S31(1620):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))
    return
  end subroutine HS31_MAID05

  subroutine HF15_MAID05(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !       COMMON/parDFp/ X1D13p,X3D13p,XSD13p,X1F15p,X3F15p,XSF15p
    !       COMMON/parDFn/ X1D13n,X3D13n,XSD13n,X1F15n,X3F15n,XSF15n
    !c ***************************************************************
    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 10.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-24.7443/COSR   *X1F15p
    A30=132.2624/COSR   *X3F15p
    S10=-25.65/COSR     *XSF15p

    A1= A10*(1.+3.72353*Q2G)*exp(-1.2*Q2G)
    A3= A30*(1.+1.352305*Q2G)*exp(-2.22357*Q2G)
    S1= S10*(1.+4.47896*Q2G)*exp(-1.55*Q2G)
    GO TO 20
    !c ***************************************************************
10  CONTINUE
    Phi_R = 15.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=26.94/COSR   *X1F15n
    A30=-37.07/COSR  *X3F15n
    S10=1./COSR*XSF15n

    A1= A10*(1.+0.001*Q2G)*exp(-1.2*Q2G)
    A3= A30*(1.+3.*Q2G)*exp(-1.75*Q2G)
    S1= S10 *exp(-1.55*Q2G)
    !c *******************************************************
20  AE=-(1./3.)*( A3*sqrt(2.)+A1)
    AM=-(1./3.)*( A3/sqrt(2.)-A1)
    AS=-(1./3.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123) Phi_R,  X1F15p, A1, X3F15p, A3, XSF15p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123) Phi_R,  X1F15n, A1, X3F15n, A3, XSF15n, S1
123 format( 5x,'F15(1680):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF15_MAID05

  subroutine HP31_MAID05(Q2G,PhiR,AM,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !COMMON/parPPp/ X1P13p,X3P13p,XSP13p,X1P31,XSP31

    PI=3.1415926536D0
    Phi_R=35.
    PhiR = Phi_R*PI/180.
    COSR=COS(PhiR)
    A10=14.786/COSR*X1P31
    S10=1./COSR*XSP31

    A1= A10*(1.+0.*Q2G)*exp(-2.*Q2G)
    S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)

    AM= A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1P31, A1, XSP31, S1
123 format( 5x,'P31(1910):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    Return
  End subroutine HP31_MAID05

  subroutine HD15_MAID05(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !COMMON/parD15/ X1D15p,X3D15p,XSD15p,X1D15n,X3D15n,XSD15n

    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 20.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)
    A10=14.356/COSR  *X1D15p
    A30=20.322/COSR  *X3D15p
    S10=1./COSR*XSD15p


    A1= A10*(1.+0.1*Q2G)*exp(-2.*Q2G)
    A3= A30*(1.+0.1*Q2G)*exp(-2.*Q2G)
    S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)
    GO TO 20
    !c ********************************************
10  CONTINUE
    Phi_R = 0.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)
    A10=-61.738/COSR *X1D15n
    A30=-83.868/COSR *X3D15n
    S10=-1./COSR*XSD15n

    A1= A10*(1.+0.01*Q2G)*exp(-2.*Q2G)
    A3= A30*(1.+0.01*Q2G)*exp(-2.*Q2G)
    S1= S10*(1.+0.01*Q2G)*exp(-2.*Q2G)


    !c *************************************************
20  AE= (1./3.)*(A3/sqrt(2.)-A1)
    AM=-(1./3.)*(A3*sqrt(2.)+A1)
    AS=-(1./3.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123) Phi_R,  X1D15p, A1, X3D15p, A3, XSD15p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123) Phi_R,  X1D15n, A1, X3D15n, A3, XSD15n, S1
123 format( 5x,'D15(1675):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HD15_MAID05

  subroutine HF35_MAID05(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !COMMON/parF3/ X1F35,X3F35,XSF35,X1F37,X3F37,XSF37

    PI=3.1415926536D0
    Phi_R =40.
    PhiR =Phi_R*Pi/180.
    COSR=COS(PhiR)
    A10=13.934/COSR*X1F35
    A30=-21.427/COSR*X3F35
    S10=1./COSR*XSF35

    A1= A10*(1.+0.*Q2G)*exp(-2.*Q2G)
    A3= A30*(1.+0.*Q2G)*exp(-2.*Q2G)
    S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)

    AE=-(1./3.)*(A3*sqrt(2.)+A1)
    AM=-(1./3.)*(A3/sqrt(2.)-A1)
    AS=-(1./3.)*S1*sqrt(2.)

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1F35, A1, X3F35, A3, XSF35, S1
123 format( 5x,'F35(1905):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF35_MAID05

  subroutine HF37_MAID05(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !COMMON/parF3/ X1F35,X3F35,XSF35,X1F37,X3F37,XSF37

    PI=3.1415926536D0
    Phi_R = 30.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-81.06/COSR*X1F37
    A30=-104.65/COSR*X3F37
    S10=1./COSR*XSF37

    A1= A10*(1.+0.*Q2G)*exp(-2.*Q2G)
    A3= A30*(1.+0.*Q2G)*exp(-2.*Q2G)
    S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)

    AE= (1./4.)*(A3*3./sqrt(15.)-A1)
    AM=-(1./4.)*(A3*5./sqrt(15.)+A1)
    AS=-(1./4.)*Sqrt(2.)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1F37, A1, X3F37, A3, XSF37, S1
123 format( 5x,'F37(1950):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF37_MAID05

  subroutine HP13_MAID05(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !COMMON/parPPp/ X1P13p,X3P13p,XSP13p,X1P31,XSP31
    !COMMON/parPPn/ X1P13n,X3P13n,XSP13n

    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 0.
    PhiR = Phi_R*PI/180.
    COSR = COS(PhiR)

    A10=73.002/COSR   *X1P13p
    A30=-11.465/COSR  *X3P13p
    S10=-53.03/COSR   *XSP13p

    A1= A10*(1.+1.891651*Q2G)*EXP(-1.55*Q2G)
    A3= A30*(1.+15.9896*Q2G)*EXP(-1.55*Q2G)
    S1= S10*(1.+2.45819*Q2G)*EXP(-1.55*Q2G)
    GO TO 20
    !c ***************************************************
10  CONTINUE
    Phi_R = 0.
    PhiR = Phi_R*PI/180.
    COSR = COS(PhiR)

    A10=-2.904/COSR  *X1P13n
    A30=-30.972/COSR *X3P13n
    S10=-1./COSR   *XSP13n

    A1= A10*(1.+12.72411*Q2G)*EXP(-1.55*Q2G)
    A3= A30*(1.+4.987*Q2G)*EXP(-1.55*Q2G)
    S1= S10*EXP(-1.55*Q2G)
    !c ****************************************************

20  AE= (1./2.)*(A3/sqrt(3.)-A1)
    AM=-(1./2.)*(A3*sqrt(3.)+A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123) Phi_R, X1P13p, A1, X3P13p, A3, XSP13p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123) Phi_R, X1P13n, A1, X3P13n, A3, XSP13n, S1
123 format( 5x,'P13(1720):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HP13_MAID05



  !****************************************************************************
  !****************************************************************************
  !****************************************************************************


  SUBROUTINE HEL_OUT(ISO,Q2G)

    ! Q2G -- Q^2 in GeV
    ! ISO -- 1=proton
    !        2=neutron

    IMPLICIT real(8)(A-H,O-Z)
    DOUBLE PRECISION mi, kgcm0
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,2) Q2G
2   FORMAT(5x,'proton e.m. helicity amplitudes at Q^2=',F7.3,2x,'in units 10^-3/sqrt(GeV)')
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,3) Q2G
3   FORMAT(5x,'neutron e.m. helicity amplitudes at Q^2=',F7.3,2x,'in units 10^-3/sqrt(GeV)')
    write(6,4)
4   format(17X,'Phi_R',6x,'X1',5x,'A1/2',8x,'X3',5x,'A3/2',8x,'XS', 5x,'S1/2')

    W0=1.232
    mi=0.9382723
    kgcm0 = (W0*W0-mi*mi)/2./W0
    egcm = (W0*W0-Q2G-mi*mi)/2./W0
    qcm=sqrt(egcm**2+Q2G)

    CALL HP33 (Q2G,qcm,real(kgcm0,8),DE,DM,DS,A1,A3,S1,1)
    CALL HP11 (ISO,PhiR,Q2G,DM,DS,A1,S1,1)
    CALL HD13 (ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
    CALL HS11f(ISO,PhiR,Q2G,DM,DS,A1,S1,1)
    CALL HS31 (Q2G,PhiR,DM,DS,A1,S1,1)
    CALL HS11s(ISO,PhiR,Q2G,DM,DS,A1,S1,1)
    CALL HD15 (ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
    CALL HF15 (ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
    CALL HD33 (Q2G,PhiR,DE,DM,DS,A1,A3,S1,1)
    CALL HP13 (ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
    CALL HF35 (Q2G,PhiR,DE,DM,DS,A1,A3,S1,1)
    CALL HP31 (Q2G,PhiR,DM,DS,A1,S1,1)
    CALL HF37 (Q2G,PhiR,DE,DM,DS,A1,A3,S1,1)

    RETURN
  END SUBROUTINE HEL_OUT


  !****************************************************************************
  !****s* helicityAmplitudes/get_helicityAmplitudes
  ! NAME
  ! subroutine get_helicityAmplitudes(targetCharge,resonanceID,QSquared,A1,A3,S1,MAID_version_in)
  ! PURPOSE
  ! * Interface for the MAID helicity amplitudes. Returns the helicity amplitudes for a given resonance.
  !
  ! INPUTS
  ! * real, intent(in) :: QSquared        ! QSquared in GeV**2
  ! * integer, intent(in) :: targetCharge ! Charge of target nucleon
  ! * integer, intent(in) :: resonanceID  ! ID of the resonance according to GiBUU
  ! * integer,optional, intent(in) :: MAID_version_in
  !   MAID_version = 1 : use MAID 2003
  !   MAID_version = 2 : use MAID 2005
  ! OUTPUT
  ! * real, intent(out):: A1,A3,S1 ! In units of GeV^(-1/2)
  !
  ! NOTES
  ! * In this routine I am converting from real to real(8) since the MAID routines are made for
  !   real(8) input and output.
  !****************************************************************************

  subroutine get_helicityAmplitudes(targetCharge,resonanceID,QSquared_IN,A1_OUT,A3_OUT,S1_OUT,MAID_version_in)
    use IDTABLE
    implicit none
    real, intent(in) :: QSquared_IN        ! QSquared in GeV**2
    integer, intent(in) :: targetCharge ! Charge of target nucleon
    integer, intent(in) :: resonanceID  ! ID of the resonance
    real, intent(out):: A1_OUT,A3_OUT,S1_OUT
    real(8) :: A1,A3,S1
    real(8) :: QSquared        ! QSquared in GeV**2

    integer :: iso, printInfos
    real(8) :: egcm,qcm,kgcm0,DE,DM,DS,phiR
    real(8), parameter :: W0=1.232
    real(8), parameter :: mi=0.9382723

    integer,optional, intent(in) :: MAID_version_in

    integer :: MAID_version   !MAID_version = 1 : use MAID 2003
                              !MAID_version = 2 : use MAID 2005

    if (present(MAID_version_in)) then
       if (MAID_version_in.eq.1.or.MAID_version_in.eq.2) then
          MAID_version=MAID_version_in
       else
          MAID_version = 2
       end if
    else
       MAID_version = 2
    end if


    !Checks
    if (.not.(targetCharge.eq.0.or.targetCharge.eq.1)) then
       write(*,*) 'Wrong nucleon charge in get_helicityAmplitudes', targetCharge
    end if

    QSquared=QSquared_In

    A1=0.
    A3=0.
    S1=0.

    iso=2-targetCharge

    printInfos=0
    ! 0= no printout
    ! 1= printout

    if (MAID_version.eq.1) then

       select case (resonanceID)
       case (Delta)
          kgcm0 = (W0*W0-mi*mi)/2./W0
          egcm = (W0*W0-QSquared-mi*mi)/2./W0
          qcm=sqrt(egcm**2+QSquared)
          CALL HP33 (QSquared,qcm,kgcm0,DE,DM,DS,A1,A3,S1,printInfos)
       case (P11_1440)
          CALL HP11 (ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case (D13_1520)
          CALL HD13 (ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (S11_1535)
          CALL HS11f(ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case ( S31_1620)
          CALL HS31 (QSquared,PhiR,DM,DS,A1,S1,printInfos)
       case (S11_1650)
          CALL HS11s(ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case (D15_1675)
          CALL HD15 (ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (F15_1680)
          CALL HF15 (ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (D33_1700)
          CALL HD33 (QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       case (P13_1720)
          CALL HP13 (ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (F35_1905)
          CALL HF35 (QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       case (P31_1910)
          CALL HP31 (QSquared,PhiR,DM,DS,A1,S1,printInfos)
       case (F37_1950)
          CALL HF37 (QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       end select

    else if (MAID_version.eq.2) then

       select case (resonanceID)
       case (Delta)
          kgcm0 = (W0*W0-mi*mi)/2./W0
          egcm = (W0*W0-QSquared-mi*mi)/2./W0
          qcm=sqrt(egcm**2+QSquared)
          CALL HP33_MAID05 (QSquared,qcm,kgcm0,DE,DM,DS,A1,A3,S1,printInfos)
       case (P11_1440)
          CALL HP11_MAID05 (ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case (D13_1520)
          CALL HD13_MAID05 (ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (S11_1535)
          CALL HS11f_MAID05 (ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case ( S31_1620)
          CALL HS31_MAID05 (QSquared,PhiR,DM,DS,A1,S1,printInfos)
       case (S11_1650)
          CALL HS11s_MAID05(ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case (D15_1675)
          CALL HD15_MAID05 (ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (F15_1680)
          CALL HF15_MAID05 (ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (D33_1700)
          CALL HD33_MAID05 (QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       case (P13_1720)
          CALL HP13_MAID05 (ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (F35_1905)
          CALL HF35_MAID05 (QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       case (P31_1910)
          CALL HP31_MAID05 (QSquared,PhiR,DM,DS,A1,S1,printInfos)
       case (F37_1950)
          CALL HF37_MAID05 (QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       end select

    else
       write(*,*) 'MAID version ', MAID_version, ' makes no sense -> STOP'
       stop
    end if


    ! Converting to GeV^(-1/2)
    A1_OUT=A1/1000.
    A3_OUT=A3/1000.
    S1_OUT=S1/1000.

    !we use different sign in S_1/2 amplitude, thus:
    S1_OUT=-S1_OUT


  end subroutine get_helicityAmplitudes
end module helicityAmplitudes
