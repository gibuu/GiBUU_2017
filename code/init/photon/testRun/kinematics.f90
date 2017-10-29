program kinematics

  use constants
  use hist2df90

  implicit none

  type(histogram2D) :: H2D1,H2D2

  real :: eps,W,Q2,nu

  double precision sqrtsLab, thetaLab, phiLab, betaLab(3)
  double precision MP_P ! prototype
  double precision MP_SqrtS ! prototype

  eps = 0.99
  W = 2.50
  Q2 = 1.00
  nu = 4.5

  call CreateHist2D(H2D1,"pT2 vz zH (1)", &
       &(/0.0,0.0/), (/1.0,1.8/), (/0.01,0.01/) , .true.)
  call CreateHist2D(H2D2,"pT2 vz zH (2)", &
       &(/0.0,0.0/), (/1.0,1.8/), (/0.01,0.01/) , .true.)

  call InitPythiaSetMasses(11,2212)   ! Target is proton

!  call InitPythiaSetKinV3(eps,W,Q2)
  call InitPythiaSetKinV5(eps,nu,Q2)
  call InitPythiaCalc

  call purify()

  call InitPythiaDump(6)

  call MP_Write(6, 1,4)

  ! 1: electron
  ! 2: nucleon
  ! 3: photon
  ! 4: scattered electron

  call MP_CalcROBO(2,2, sqrtsLab, thetaLab, phiLab, betaLab)
  call MP_ROBO_Inv(1,4, thetaLab,phiLab,betaLab(1),betaLab(2),betaLab(3))
  call MP_Write(6, 1,4)

  nu = MP_P(3,4)
  W  = MP_SqrtS(2,3)

  write(*,*) ">>>>>>>> nu, W : ",nu, W


!  call scenario1(101,0.938) ! Nucleon
!  call scenario1(102,1.232) ! Delta
  call scenario1(103,1.076) ! pi+N

  call WriteHist2D_Gnuplot(H2D1,141,mul=1., add=1e-3)

  call scenario2(111,0.770) ! rho

  call WriteHist2D_Gnuplot(H2D2,142,mul=1., add=1e-20)

contains
  subroutine purify()

  double precision MP_P,MP_M ! prototype

  call MP_Set3(2, MP_M(2), -1e-20, 0.0, MP_P(2,3))
  call MP_Set4(3, MP_M(3),  1e-20, 0.0, MP_P(3,3),MP_P(3,4))

  end subroutine purify


  subroutine scenario1(iFile,mA)

    implicit none

    integer :: iFile

    double precision MP_SqrtS ! prototype
    double precision MP_P ! prototype

    real :: srts,mA,mB, pA, theta
    integer :: i1

    integer, parameter :: n1 = 1000

    real :: zH, pT2,weight
    real :: theta1

    double precision sqrtsRho, thetaRho, phiRho, betaRho(3)

!    mA = 0.938 ! Nucleon
!    mA = 1.232 ! Delta
    mB = 0.138


    srts = MP_SqrtS(2,3)

    pA = sqrt((srts**2-mA**2-mB**2)**2-4*mA**2*mB**2)/(2*srts)


    weight = 1.0/(pi*n1)

    do i1=0,n1
       theta1 = (pi/n1)*i1


       call MP_Set3(5,mA, 0.0, 0.0, pA)
       call MP_Set3(6,mB, 0.0, 0.0,-pA)

       call MP_ROBO(5,6, theta1,0.0, 0.0,0.0,0.0)


!       call MP_Write(6, 5,6)

       call MP_ROBO_Inv(5,6, thetaLab,phiLab,betaLab(1),betaLab(2),betaLab(3))

       call MP_Write(6, 5,6)

       zH = MP_P(6,4)/nu
       pT2= MP_P(6,1)**2+MP_P(6,2)**2

       write(iFile,*) zH,sqrt(pT2)

       call AddHist2D(H2D1, (/zH,pT2/),weight*sin(theta1))

    enddo

  end subroutine scenario1

  subroutine scenario2(iFile,mA)

    implicit none

    integer :: iFile

    double precision MP_SqrtS ! prototype
    double precision MP_P ! prototype

    real :: srts,mA,mB,mC, pA, theta, pC
    integer :: i1,i2,i3
    integer,parameter :: n1=200,n2=200,n3=400
    real :: theta1,theta2,phi2


    real :: zH, pT2,weight

    double precision sqrtsRho, thetaRho, phiRho, betaRho(3)


    mB = 0.938
    mC = 0.138

    srts = MP_SqrtS(2,3)

    pA = sqrt((srts**2-mA**2-mB**2)**2-4*mA**2*mB**2)/(2*srts)
    pC = sqrt((mA**2  -mC**2-mC**2)**2-4*mC**2*mC**2)/(2*mA)

    weight = 1.0/(pi**3 * n1*n2*2*n3)

    do i1=0,n1

       call MP_Set3(5,mA, 0.0,0.0, pA)
       call MP_Set3(6,mB, 0.0,0.0,-pA)

       theta1 = (pi/n1)*i1

       call MP_ROBO(5,6, theta1,0.0, 0.0,0.0,0.0)

       call MP_Write(6, 5,6)

       call MP_CalcROBO(5,5, sqrtsRho, thetaRho, phiRho, betaRho)

!       call MP_ROBO_Inv(5,5, thetaRho,phiRho,betaRho(1),betaRho(2),betaRho(3))

       do i2= 0,n2
          theta2 = (pi/n2)*i2
          do i3 = 0,n3
             phi2 = (2*pi/n3)*i3

             call MP_Set3(7,mC, 0.0, 0.0, pC)
             call MP_Set3(8,mC, 0.0, 0.0,-pC)
             call MP_ROBO(7,8, theta2,phi2, 0.0,0.0,0.0)

             call MP_ROBO(7,8, thetaRho,phiRho,betaRho(1),betaRho(2),betaRho(3))
             
!             call MP_Write(6, 1,8)
             
             call MP_ROBO_Inv(5,8, thetaLab,phiLab,betaLab(1),betaLab(2),betaLab(3))
             
!             call MP_Write(6, 1,8)


             zH = MP_P(7,4)/nu
             pT2= MP_P(7,1)**2+MP_P(7,2)**2
!             write(iFile,*) zH,sqrt(pT2)
!             write(iFile,*) zH,pT2

             call AddHist2D(H2D2, (/zH,pT2/),weight*sin(theta1)*sin(theta2))


          enddo


       enddo
 


    end do

  end subroutine scenario2

end program kinematics
