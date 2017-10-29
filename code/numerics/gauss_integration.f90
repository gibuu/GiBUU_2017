!******************************************************************************
!****m* /gauss_integration
! NAME
! module gauss_integration
! PURPOSE
! This module includes the routines needed for Gauss integration.
!******************************************************************************
module gauss_integration

  implicit none
  private

  public :: sg20r, rg20r, sg64r, rg64r

contains

  !****************************************************************************
  ! SUBROUTINAS DE INTEGRACION NUMERICA POR GAUSS(SINGLE PRECISION)
  !****************************************************************************
  SUBROUTINE SG20R(A,B,N,X,NP)

    real, intent(in) :: A,B
    integer, intent(in) :: N
    integer, intent(out) :: NP
    real, intent(out),dimension(:) :: X

    real,dimension(10),parameter :: Y=  &
         & (/.9931285991,.9639719272,.9122344282,.8391169718, &
         & .7463319064,.6360536807,.5108670019,.3737060887, &
         & .2277858511,.0765265211/)

    integer :: I1,I2,I,J,J1,J2
    real :: DINT,DELT,ORIG,DORIG

    NP=20*N
    DINT=(B-A)/float(N)
    DELT=DINT*0.5
    ORIG=A-DELT
    I1=-20
    do I=1,N
       ORIG=ORIG+DINT
       DORIG=ORIG+ORIG
       I1=I1+20
       I2=I1+21
       do J=1,10
          J1=I1+J
          J2=I2-J
          X(J1)=ORIG-DELT*Y(J)
          X(J2)=DORIG-X(J1)
       end do
    end do
  END SUBROUTINE SG20R


  SUBROUTINE RG20R(A,B,N,CF,CRES)
    real, intent(in) :: A,B
    integer, intent(in) :: N
    real, intent(in),dimension(:) :: CF
    real, intent(out) :: CRES
    integer :: I1,I,I2,J,J1,J2
    real :: CR

    real,dimension(10),parameter :: W= &
         & (/.0176140071,.0406014298,.0626720483,.0832767415, &
         & .1019301198,.1181945319,.1316886384,.1420961093,.1491729864,  &
         & .1527533871/)

    CR=0.
    I1=-20
    do I=1,N
       I1=I1+20
       I2=I1+21
       do J=1,10
          J1=I1+J
          J2=I2-J
          CR=CR+W(J)*(CF(J1)+CF(J2))
       end do
    end do
    CRES=CR*0.5*(B-A)/float(N)
  END SUBROUTINE RG20R





  !****************************************************************************
  ! the same as before but with 64 gauss points instead of 20
  !****************************************************************************


  SUBROUTINE SG64R(A,B,N,X,NP)

    real, intent(in) :: A,B
    integer, intent(in) :: N
    integer, intent(out) :: NP
    real, intent(out),dimension(:) :: X

    real,dimension(32),parameter :: Y=  &
         & (/ 0.999305041735772, 0.996340116771955, 0.991013371476744, 0.983336253884626, &
              & 0.973326827789911, 0.961008799652054, 0.946411374858403, 0.929569172131940, &
              & 0.910522137078503, 0.889315445995114, 0.865999398154093, 0.840629296252580, &
              & 0.813265315122798, 0.783972358943341, 0.752819907260532, 0.719881850171611, &
              & 0.685236313054233, 0.648965471254657, 0.611155355172393, 0.571895646202634, &
              & 0.531279464019895, 0.48940314570705,  0.446366017253464, 0.402270157963992, &
              & 0.357220158337668, 0.311322871990211, 0.264687162208767, 0.217423643740007, &
              & 0.169644420423993, 0.121462819296121, 0.072993121787799, 0.0243502926634244   /)

    integer :: I1,I2,I,J,J1,J2
    real :: DINT,DELT,ORIG,DORIG

    NP=64*N
    DINT=(B-A)/float(N)
    DELT=DINT*0.5
    ORIG=A-DELT
    I1=-64
    do I=1,N
       ORIG=ORIG+DINT
       DORIG=ORIG+ORIG
       I1=I1+64
       I2=I1+65
       do J=1,32
          J1=I1+J
          J2=I2-J
          X(J1)=ORIG-DELT*Y(J)
          X(J2)=DORIG-X(J1)
       end do
    end do
  END SUBROUTINE SG64R



  SUBROUTINE RG64R(A,B,N,CF,CRES)
    real, intent(in) :: A,B
    integer, intent(in) :: N
    real, intent(in),dimension(:) :: CF
    real, intent(out) :: CRES
    integer :: I1,I,I2,J,J1,J2
    real :: CR

    real,dimension(32),parameter :: W= &
         & (/0.00178328072169643,0.00414703326056247,0.00650445796897836,0.00884675982636395, &
         &   0.0111681394601311, 0.0134630478967186, 0.0157260304760247, 0.0179517157756973,  &
         &   0.0201348231535302, 0.0222701738083833, 0.0243527025687109, 0.0263774697150547,  &
         &   0.0283396726142595, 0.0302346570724025, 0.0320579283548516, 0.0338051618371416,  &
         &   0.0354722132568824, 0.0370551285402400, 0.0385501531786156, 0.0399537411327203,  &
         &   0.0412625632426235, 0.0424735151236536, 0.0435837245293235, 0.0445905581637566,  &
         &   0.0454916279274181, 0.0462847965813144, 0.0469681828162100, 0.0475401657148303,  &
         &   0.0479993885964583, 0.0483447622348030, 0.0485754674415034, 0.0486909570091397     /)

    CR=0.
    I1=-64
    do I=1,N
       I1=I1+64
       I2=I1+65
       do J=1,32
          J1=I1+J
          J2=I2-J
          CR=CR+W(J)*(CF(J1)+CF(J2))
       end do
    end do
    CRES=CR*0.5*(B-A)/float(N)
  END SUBROUTINE RG64R


end module gauss_integration
