!******************************************************************************
!****m* /NDeltaFF_Ramalho
! NAME
! module NDeltaFF_Ramalho
! FUNCTION
! Calculate the Delta -> N gamma* transition form factor.
! "Model 4" - new pion cloud parametrization (private communication).
! Author: Gilberto Ramalho.
!******************************************************************************
module NDeltaFF_Ramalho

  implicit none
  private

  public :: NDeltaSL  ! NDelta FF in spacelike region
  public :: NDeltaTL  ! NDelta FF in timelike region

  real(8), parameter :: PI = 3.14159265358979323846264338327950288419716939937510d0

  ! mass parameters
  real(8), parameter :: mpi  = 0.138d0
  real(8), parameter :: mN   = 0.939d0
  real(8), parameter :: mD0  = 1.232d0
  real(8), parameter :: mrho = 0.775d0

  ! decay widths
  real(8), parameter :: drho = 0.1491d0   ! rho decay width
  real(8), parameter :: dmx  = 0.5964d0   ! effective decay width (4*drho)

  ! parameters of the quark current
  real(8), parameter :: cp = 4.1600d0
  real(8), parameter :: cm = 1.1150d0
  real(8), parameter :: dp = -0.6860d0
  real(8), parameter :: dm = -0.6860d0

  real(8), parameter :: kp = 1.6410d0
  real(8), parameter :: km = 1.8240d0
  real(8), parameter :: lamb = 1.2100d0

  ! parameters of the (radial) wave function
  real(8), parameter :: N0 = 3.35696291311   ! normalization constant
  real(8), parameter :: beta1 = 0.0490d0
  real(8), parameter :: beta2 = 0.7170d0

  real(8), parameter :: alf1 = 0.33660d0
  real(8), parameter :: alf2 = 0.33660d0      ! model from PRD 80, 013008 (2009)
  real(8) :: NS

  ! pion cloud parameters: model from PRD 78, 114017 (2008)
  real(8), parameter :: lpiD = 0.4410d0
  real(8), parameter :: LLpiD = 1.530D0/mN**2  ! in natural units

  logical, save :: first = .true.  ! indicates whether initialization is required

contains

  !****************************************************************************
  ! Subroutine INIT
  ! initialize Delta normalization constant
  !****************************************************************************
  SUBROUTINE Init
    NS = Normaliza()
    first = .false.
  END SUBROUTINE


  !****************************************************************************
  ! NDeltaTL
  ! calculate NDelta FF in timelike region q2t > 0
  !
  ! input:
  ! * q2t = q^2
  ! * W   = Delta mass
  ! output:
  ! * GM2  = abs(G_M*)**2
  !****************************************************************************
  real(8) FUNCTION NDeltaTL(q2t,W) result(GM2)
    real(8), intent(in)  :: q2t,W

    real(8) :: Q2,GMa,GMb,GMpiA,GMpiB,ax,bx,GMre,GMim

    if (first) call Init

    Q2=-q2t

    ax = 1d0  ! bare contribution       (1=use it, 0=don't)
    bx = 1d0  ! pion cloud contribution (1=use it, 0=don't)

    if ( q2t > (W-mN)**2.0 ) then

      GMpiA=0d0
      GMpiB=0d0

    else

      CALL DeltaTL3(Q2,NS,W,GMa,GMb)

      if (abs(Q2) > 1d-6) then
        CALL PionCloudR4(Q2,GMpiA,GMpiB)
      else
        CALL PionCloudR4S(GMpiA,GMpiB)
      end if

    end if

    GMre=ax*GMa + bx*GMpiA
    GMim=ax*GMb + bx*GMpiB
    GM2=GMre**2.0 + GMim**2.0

  END FUNCTION


  !****************************************************************************
  ! NDeltaSL
  ! calculate NDelta FF in spacelike region q2t < 0
  !
  ! input:
  ! * q2t = q^2
  ! * W   = Delta mass
  ! output:
  ! * GM2 = abs(G_M*)**2
  !****************************************************************************
  real(8) FUNCTION NDeltaSL(q2t,W) result(GM2)
    real(8), intent(in)  :: q2t,W

    real(8) :: Q2,GMb,GMpi,GM,ax,bx

    if (first) call Init

    ax = 1d0  ! bare contribution       (1=use it, 0=don't)
    bx = 1d0  ! pion cloud contribution (1=use it, 0=don't)

    if (q2t > 0d0) then
      write(*,*)
      write(*,*) 'wrong sign for q2t'
      stop
    end if

    if (abs(q2t) < 1d-4) then
      Q2=1d-4
    else
      Q2=-q2t
    end if

    GMb  = DeltaFF0(Q2,NS,W)
    GMpi = PionCloudR(Q2)

    GM=ax*GMb + bx*GMpi
    GM2=GM*GM

  END FUNCTION


  !****************************************************************************
  ! Time-like region
  ! gR=gR(Q2)
  !****************************************************************************
  SUBROUTINE DeltaTL3(Q2,NS,W,GMa,GMb)
    real(8), intent(in)  :: Q2, NS, W
    real(8), intent(out) :: GMa,GMb

    integer, parameter :: nmax=400
    real(8) :: mD
    real(8) :: q2vec,qvec,EN,omega
    integer :: Nk,Nz,ik,iz
    real(8) :: GMS
    real(8) :: mdel,mdel2,dm2,sm2,sm,dmp
    real(8) :: intfac,fac,GD3,F1m,F2m,jm
    real(8), dimension(nmax) :: Xk,WX,U,WU
    real(8) :: cut,Q2m
    real(8) :: k,dk,Es,z,dz,k2,z2,intS
    real(8) :: chi1,psi1,chi2,psi2
    real(8) :: mrho2,Rfactor,gR
    real(8) :: Rfactor2A,Rfactor2B,mh2,gmh
    real(8)::  Rfactor3A,Rfactor3B,Rfactor2
    real(8) :: RfactorA,RfactorB,F1mA,F1mB,F2mA,F2mB
    real(8) :: jmA,jmB,Q2R

    gmh=4d0*(-Q2/(4d0*mN**2.0-Q2))**2.0*dmx

    ! Q2R threshold for gR(Q2)
    Q2R = 4*mpi**2

    gR=0d0

    if ( Q2 < - Q2R) then
        gR=((-Q2R-Q2)/(-Q2R+mrho**2.0))**1.5*mrho**2.0/(-Q2)*drho
    end if

    mD=W

    q2vec=(W**2.0+ mN**2.0+ Q2)**2.0/(4d0*W**2.0)-mN**2.0
    qvec=sqrt(q2vec)

    omega=(W**2.0-mN**2.0-Q2)/(2d0*W)

    En=sqrt(mN**2.0+q2vec)

    ! Constants  ----------------------------------------------------
    mdel=mD/mN        ! Delta mass (nucleon unities)
    mdel2=mdel*mdel

    dm2=mdel*mdel-1d0   ! md**2-mn**2
    sm2=mdel*mdel+1d0   ! md**2+mn**2

    dmp=mdel-1d0
    sm=mdel+1d0

    ! Integral factor
    intfac=.5d0/(2d0*PI)**2.0
    GD3=3d0/(1d0+Q2/0.71d0)**2.0

    ! Integral coeficient
    fac=1d0/sqrt(3d0)

    ! Integration parameters ----------------------------------------

    Nk=200     ! k grid
    Nz=300     ! z grid

    CALL GAUSS( 0.,1.,Xk,WX,Nk)
    CALL GAUSS(-1.,1.,U,WU,Nz)

    cut=0.5d0

    Q2m=Q2/mN**2     ! Normalized momenta

    ! current functions

    mrho2=(mrho/mN)**2.0

    Rfactor=mrho**2/(mrho**2+Q2)

    RfactorA = mrho**2.0*(mrho**2+Q2) / ( (mrho**2+Q2)**2.0 + mrho**2*gR**2)
    RfactorB = mrho**2.0*mrho*gR / ( (mrho**2+Q2)**2.0 + mrho**2*gR**2)

    mh2=4d0*mN**2.0

    Rfactor2=mh2/((mh2+Q2)**2.0+mh2*gmh**2.0)

    Rfactor3A=(mh2+Q2)*Rfactor2
    Rfactor3B=sqrt(mh2)*gmh*Rfactor2

    Rfactor2A = Rfactor2**2.0 * ( (mh2+Q2)**2.0 - Mh2*gmh**2.0)
    RFactor2B = Rfactor2**2.0 * ( 2d0*(mh2+Q2)*sqrt(mh2)*gmh )

    F1m = lamb + (1d0-lamb)*Rfactor + cm*.25d0*Q2m/(1d0+.25d0*Q2m)**2.0
    F2m = km * ( dm*Rfactor + (1d0-dm)/(1d0+.25d0*Q2m))

    ! Real part

    F1mA = lamb + (1d0-lamb)*RfactorA + cm*Rfactor2A*Q2/mh2
    F2mA = km * ( dm*RfactorA + (1d0-dm)*Rfactor3A )

    ! Imaginary part

    F1mB = lamb*0 + (1d0-lamb)*RfactorB + cm*Rfactor2B*Q2/mh2
    F2mB = km * ( dm*RfactorB + (1d0-dm)*Rfactor3B )

    ! Initiate integration

    intS=0d0

    do ik=1,Nk
        k=cut*Xk(ik)/(1d0-Xk(ik))
        k2=k*k
        dk=cut/(1d0-Xk(ik))**2.0*WX(ik)
        Es=sqrt(1d0+k2)
        do iz=1,Nz
          z=U(iz)
          z2=z*z
          dz=WU(iz)

          ! Wave functions (Delta rest frame)
          chi1=2d0*Es
          psi1=1d0/((alf1-2d0+chi1)*(alf2-2d0+chi1)**2.0)

          chi2=2d0*(En*Es-qvec*k*z)/mN
          psi2=1d0/((beta1-2d0+chi2)*(beta2-2d0+chi2))

          intS=intS+psi1*psi2*k2*dk*dz/Es

        end do   ! End of iz cicle
    end do      ! End of ik cicle

    intS=  N0*NS*intfac*intS

    jm=F1m+sm*F2m/2d0     ! F1m (normalized to 1)

    jmA=F1mA+sm*F2mA/2d0
    jmB=F1mB+sm*F2mB/2d0

    ! S-state Form Factors
    GMS=8d0*jm/(3d0*sm)*intS*fac

    GMa=8d0*jmA/(3d0*sm)*intS*fac
    GMb=8d0*jmB/(3d0*sm)*intS*fac

  end subroutine


  !****************************************************************************
  ! Routine DeltaFF0
  ! Input: Q2,NS,mD
  !
  ! Output: GM
  !
  ! Uses Delta rest frame coordenates
  ! SS transition (frame invariant): intS
  !****************************************************************************
  real(8) FUNCTION DeltaFF0(Q2,NS,mD) result(GM)
    real(8), intent(in) :: Q2, NS, mD

    integer, parameter :: nmax=400
    integer :: Nk,Nz,ik,iz
    real(8) :: GMS
    real(8) :: mdel,mdel2,dm2,sm2,sm,dmp
    real(8) :: intfac,fac,GD3,F1m,F2m,jm
    real(8), dimension(nmax) :: Xk,WX,U,WU
    real(8) :: cut,Q2m,qD,qD2,En
    real(8) :: k,dk,Es,z,dz,k2,z2,intS
    real(8) :: chiN,chiD,psiN,psiD
    real(8) :: mrho2

    ! Constants  ----------------------------------------------------
    mdel=mD/mN        ! Delta mass (nucleon unities)
    mdel2=mdel*mdel

    dm2=mdel*mdel-1d0   ! md**2-mn**2
    sm2=mdel*mdel+1d0   ! md**2+mn**2

    dmp=mdel-1d0
    sm=mdel+1d0

    ! Integral factor
    intfac=.5d0/(2d0*PI)**2.0
    GD3=3d0/(1d0+Q2/0.71d0)**2.0


    ! Integral coeficient
    fac=1d0/sqrt(3d0)

    ! Integration parameters
    Nk=200     ! k grid
    Nz=300     ! z grid

    CALL GAUSS( 0.,1.,Xk,WX,Nk)
    CALL GAUSS(-1.,1.,U,WU,Nz)

    cut=0.5d0

    Q2m=Q2/mN**2     ! Normalized momenta

    ! Delta rest frame variables
    qD2=.25d0*(sm*sm+Q2m)*(dmp*dmp+Q2m)/mdel2
    qD=sqrt(qD2)

    En=(mdel2+1d0+Q2m)/(2d0*mdel)

    ! current functions

    mrho2=(mrho/mN)**2.0

    ! F1m     (No pion cloud)
    F1m = lamb + (1d0-lamb)/(1d0+Q2m/mrho2) + cm*.25d0*Q2m/(1d0+.25d0*Q2m)**2.0
    ! F2m
    F2m = km * ( dm/(1d0+Q2m/mrho2) + (1d0-dm)/(1d0+.25d0*Q2m) )

    ! Initiate integration

    intS=0d0

    do ik=1,Nk
        k=cut*Xk(ik)/(1d0-Xk(ik))
        k2=k*k
        dk=cut/(1d0-Xk(ik))**2.0*WX(ik)
        Es=sqrt(1d0+k2)
        do iz=1,Nz
          z=U(iz)
          z2=z*z
          dz=WU(iz)

          ! Wave function parameters
          ! Delta Rest Frame
          chiN=2d0*En*Es+2d0*k*qD*z
          chiD=2d0*Es

          ! Nucleon wave function
          psiN=1d0/((beta1-2d0+chiN)*(beta2-2d0+chiN))

          ! Delta wave functions
          ! S-state

          psiD=1d0/((alf1-2d0+chiD)*(alf2-2d0+chiD)**2.0)

          intS=intS+psiN*psiD*k2*dk*dz/Es

        end do   ! End of iz cicle
    end do      ! End of ik cicle

    intS=  N0*NS*intfac*intS

    jm=F1m+sm*F2m/2d0     ! F1m (normalized to 1)

    ! S-state Form Factors; no pion cloud
    GMS=8d0*jm/(3d0*sm)*intS*fac

    GM=GMS

  END FUNCTION DeltaFF0


  !****************************************************************************
  ! Pion cloud in TL
  !****************************************************************************
  SUBROUTINE PionCloudR4(Q2,GMpiA,GMpiB)
    real(8), intent(in) :: Q2
    real(8), intent(out) :: GMpiA,GMpiB

    real(8) :: Q2m
    real(8) :: Rfactor,cut2,RfactorA,RfactorB,fpi
    real(8) :: t1,t2,mpi2,mrho2
    real(8) :: a,b,gm,RfactorC
    real(8) :: at,bt
    real(8) :: gmD,q2R,GD,GD2t

    at=0.83334330616420560d0
    bt=0.17809674271214604d0

    Q2m=Q2/mN**2.0

    mpi2=mpi**2.0
    mrho2=mrho**2.0

    q2R=4d0*mpi2

    ! cut_pi width defined in units mN
    gm=4d0*(-Q2m/(LLpiD-Q2m))**2.0*dmx/mN

    ! term 1 -- diagram (a)

    t1=Q2*log(-Q2/mpi2)*(bt)/pi
    t2=Q2*(bt)

    Rfactor=at**2.0/( (at**2.0 + Q2 + t1)**2.0 + t2**2.0)

    RfactorA=Rfactor*(at**2.0 + Q2 + t1)
    RfactorB=Rfactor*t2

    RfactorC=( LLpiD/((LLpiD+Q2m)**2+ LLpiD*gm**2) )**2

    fpi=lpiD*RfactorC

    a= (LLpiD+Q2m)**2.0 - LLpiD*gm**2.0
    b= 2d0*(LLpiD+Q2m)*sqrt(LLpiD)*gm

    GMpiA=fpi*(a*RfactorA - b*RfactorB)
    GMpiB=fpi*(b*RfactorA + a*RfactorB)

    ! term 2 -- diagram (a)

    cut2=.9d0

    ! effective width
    gmD=4d0*(-Q2/(cut2**2.0-Q2))**2.0*dmx

    GD=cut2**2.0/( (cut2 + Q2)**2.0 + cut2*gmD**2.0)

    GD2t=GD*GD

    ! Final result

    GMpiA= 3d0*(GMpiA + lpiD*GD2t)/2d0
    GMpiB= 3d0*GMpiB/2d0     ! no imaginary contribution from GD2t

  END SUBROUTINE


  !****************************************************************************
  ! Pion cloud in TL/SL  (limit Q2=0d0)
  !****************************************************************************
  SUBROUTINE PionCloudR4S(GMpiA,GMpiB)
    real(8), intent(out) :: GMpiA,GMpiB
    GMpiA=3d0*lpiD
    GMpiB=0d0
  END SUBROUTINE


  !****************************************************************************
  ! Pion cloud Space-like
  !****************************************************************************
  real(8) FUNCTION PionCloudR(Q2) result(GMpi)
    real(8), intent(in) :: Q2

    real(8) :: GD,Q2m,GDp
    real(8) :: at,bt,cut2

    at= 0.83439499288356689d0
    bt= 0.17821529978939199d0

    cut2=0.9d0

    GD=1d0/(1d0+Q2/cut2)**2.0

    GDp = at**2.0 / (at**2.0 + Q2 + bt*Q2/pi*log(Q2/mpi**2.0))

    Q2m=Q2/mN**2.0

    GMpi=lpiD*(LLpiD/(LLpiD+Q2m))**2*GDp

    GMpi= 3d0*(GMpi + lpiD*GD**2.0)/2d0
  END FUNCTION


  !****************************************************************************
  ! Routine Normalize
  ! Evaluate normalization constants NS
  ! INPUT: alf1,alf2
  !****************************************************************************
  real(8) FUNCTION Normaliza()
    integer, parameter :: nmax=400
    integer :: Nk,Nz,ik,iz
    real(8) :: mdel,mdel2,dm2,sm2,sm,dm
    real(8) :: intfac
    real(8), dimension(nmax) :: Xk,WX,U,WU
    real(8) :: cut,qD,qD2,En
    real(8) :: k,dk,Es,z,dz,k2,z2
    real(8) :: intS,chiD,psiD

    ! Constants  ----------------------------------------------------
    mdel=mD0/mN        ! Delta mass (nucleon unities)
    mdel2=mdel*mdel

    dm2=mdel2-1d0   ! md**2-mn**2
    sm2=mdel2+1d0   ! md**2+mn**2

    dm=mdel-1d0
    sm=mdel+1d0

    ! Integral factor
    intfac=.5d0/(2d0*PI)**2.0
    !GD3=3d0/(1d0+Q2/0.71d0)**2.0

    ! Integration parameters ----------------------------------------

    Nk=200     ! k grid
    Nz=300     ! z grid

    CALL GAUSS( 0.,1.,Xk,WX,Nk)
    CALL GAUSS(-1.,1.,U,WU,Nz)

    cut=0.5d0

    ! Delta rest frame variables Q2=0
    qD2=dm2/(2d0*mdel)
    qD=sqrt(qD2)
    En=sm2/(2d0*mdel)

    ! Initiate integration

    intS=0d0

    do ik=1,Nk
        k=cut*Xk(ik)/(1d0-Xk(ik))
        k2=k*k
        dk=cut/(1d0-Xk(ik))**2.0*WX(ik)
        Es=sqrt(1d0+k2)
        do iz=1,Nz
          z=U(iz)
          z2=z*z
          dz=WU(iz)

          ! Wave function parameters
          ! Delta Rest Frame
          chiD=2d0*Es

          ! Delta wave function initial state
          psiD=1d0/((alf1-2d0+chiD)*(alf2-2d0+chiD)**2.0)

          intS=intS+(psiD*psiD)*k2*dk*dz/Es

        end do   ! End of iz cicle
    end do      ! End of ik cicle

    intS= intfac*intS

    Normaliza = 1d0/sqrt(intS)

  END FUNCTION



  !****************************************************************************
  !     Gaussian method: calculate coordinates and weights
  !****************************************************************************
  SUBROUTINE GAUSS(X1,X2,X,W,N)
    real(8), INTENT(IN) :: X1, X2
    INTEGER, INTENT(IN) :: N
    real(8), DIMENSION(N), INTENT(OUT) :: X, W

    INTEGER :: I, J, M
    real(8) :: P1, P2, P3, PP, XL, XM, Z, Z1
    real(8), PARAMETER :: EPS = 3.D-12

    M=(N+1)/2
    XM=0.5D0*(X2+X1)
    XL=0.5D0*(X2-X1)
    do 12 I=1,M
      Z=COS(3.14159265358979D0*(I-.25D0)/(N+.5D0))
1     CONTINUE
        P1=1.D0
        P2=0.D0
        do 11 J=1,N
          P3=P2
          P2=P1
          P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11      CONTINUE
        PP=N*(Z*P1-P2)/(Z*Z-1.D0)
        Z1=Z
        Z=Z1-P1/PP
      if (abs(Z-Z1).GT.EPS) GOTO 1
      X(I)=XM-XL*Z
      X(N+1-I)=XM+XL*Z
      W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
      W(N+1-I)=W(I)
12  CONTINUE

  END SUBROUTINE

end module NDeltaFF_Ramalho
