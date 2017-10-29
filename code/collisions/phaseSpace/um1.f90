!******************************************************************************
!****m* /um1
! NAME
! module um1
!
! PURPOSE
! Analytical Phase Volume Calculation
! cf.: G.I. Kopyolv, Nucl.Phys. 36(1962)425
!******************************************************************************
module um1

  implicit none
  private

  public :: SN

  !          --- CONSTANTS FOR THE PHASE VOLUME EVALUATION ---
  !          --- G.I. KOPYLOV'70 , P.477

  real, dimension(21), parameter :: DZET = (/&
     &           0.00000, 0.34200, 0.51913, 0.61706,&
     &           0.68052, 0.72541, 0.75899, 0.78513, 0.80609,&
     &           0.82328, 0.83764, 0.84983, 0.86031, 0.86941,&
     &           0.87739, 0.88445, 0.89074, 0.89637, 0.90146,&
     &           0.90606, 0.91025/)
  real, dimension(21), parameter :: CK = (/&
     &            0.0000,  &
     &            3.1416,  16.3400,  43.4590,  78.2010, 106.6800,&
     &          117.3900, 108.3400,  86.1630,  60.2200,  37.5480,&
     &           21.1360,  10.8440,   5.1127,   2.2294,   0.90436,&
     &            0.34295,  0.12210,  0.040966, 0.012996, 0.0039103/)


contains
  !****************************************************************************
  !****f* um1/sn
  ! NAME
  ! real function SN(AMN,NP,AMI)
  ! INPUTS
  ! * real, dimension(21) :: AMI -- Masses of the particles in the initial order
  ! * integer             :: NP -- number of particles (2<=NP<=20)
  ! * real                :: AMN -- Total Energy
  !****************************************************************************
  real function SN(AMN,NP,AMI)

    real, dimension(21) :: AMI
    integer             :: NP
    real                :: AMN

    real :: DMU(21),DTAU(21),DMN(21),DMI(21),DPN(21),DEN(21),DLEV,DPRO
    real :: AMW(21),AMF(21)
    integer :: i,INDC,IP,IWAY,j,k,l,n,n1,NM

    real :: AMX



    N=NP
    DMN(N)=DBLE(AMN)
    DMU(N)=0.D+00
    do I=1,N
       AMW(I)=AMI(I)
       DMU(N)=DMU(N)+DBLE(AMW(I))
    end do
    if (DMU(N).GE.DMN(N)) GOTO 100
    !
    !           --- THE BEST ARRANGEMENT IS SYMMETRICAL ONE ---
    !
    do I=1,N
       CALL PMAXN(AMW,N,AMX,NM)
       AMW(NM)=0.
       INDC=2*MOD(I,2)-1
       IP=INT(N/2.+1.)+INT(I/2.)*INDC
       AMF(IP)=AMX
!            PRINT 20,IP,AMF(IP)
! 20         FORMAT(8X,' IP=',I4,6X,'AMF=',F12.5)
    end do
    !
    !              --- MAIN CALCULATIONS ---
    !
    N1=N-1
    DPRO=1.D+00
    do J=1,N
       DMI(J)=DBLE(AMF(J))
    end do
    do J=1,N1
       DMU(J)=0.D+00
       do L=1,J
          DMU(J)=DMU(J)+DMI(L)
       end do
    end do
    IWAY=-1
    DLEV=1.0D-04*DMU(N)
    DTAU(N)=AMN-DMU(N)
    if (DTAU(N).LE.DLEV) IWAY=1
    do J=1,N1
       K=N-J+1
       DTAU(K-1)=DTAU(K)*DZET(K-1)
       DMN(K-1)=DTAU(K-1)+DMU(K-1)
       if (IWAY.GT.0) GOTO 7
       DEN(K)=(DMN(K)**2-DMN(K-1)**2+DMI(K)**2)/2.D+00/DMN(K)
       DPN(K)=SQRT(DEN(K)**2-DMI(K)**2)
       GOTO 8
7      DPN(K)=SQRT((2.D+00*DMI(K)*DMU(K-1))*(DTAU(K)-DTAU(K-1)) /DMU(K))
! 8         PRINT 800,K,DPN(K),DTAU(K)
! 800      FORMAT(3X,'K=',I3,4X,'DPN(K)=',E10.3,4X,'DTAU(K)=',E10.3)
8      DPRO=DPRO*DPN(K)
    end do
!         PRINT 700,DPRO,CK(N),DTAU(N)
! 700     FORMAT(3X,'DPRO=',E10.3,' CK(N)=',E10.3,' DTAU(N)=',E10.3)
    SN=CK(N)*(DTAU(N)**(N-2))*DPRO/DMN(N)
    RETURN

100 SN=0.
    RETURN
  end function SN


  SUBROUTINE PMAXN(A,N,AMX,NM)
    ! --  TO FIND THE GREATEST VALUE 'AMX' AND ITS NUMBER 'NM'           --
    ! --  IN THE ARRAY 'A(N)'.  IN THE CASE OF EQUAL ELEMENTS TO ASSIGN: --
    ! --  AMX=A(1), NM=1

    real,dimension(*), intent(in) :: A
    integer, intent(in) :: N
    real, intent(out) :: AMX
    integer, intent(out) :: NM

    integer :: i

    NM=1
    AMX=A(1)
    do I=1,N
       if (A(I).LE.AMX) cycle
       AMX=A(I)
       NM=I
    end do
    RETURN
  END SUBROUTINE PMAXN


end module um1
