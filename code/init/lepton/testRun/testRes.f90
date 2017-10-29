program test
  use inputGeneral, only: readInputGeneral
  use particleProperties, only: initParticleProperties
!   use minkowski
!   use IDTABLE
  implicit none
!   INTEGER :: ID

  call readInputGeneral
  call initParticleProperties
!   call printMatrices()

  call testingAll()

!  DO ID=NUCLEON,NUCLEON+nres
 ! call testing(ID)
  !  call testing(S11_1535)
  ! call testing(P11_1440)
  !call testing(D13_1520)
  !END DO
end program test


subroutine testing(resID)
  use leptonicID
  use spectralFunc, only : specFunc
  use constants, only : GeVSquared_times_mb, pi, mN
  use minkowski, only : contract,SP
  use hadronTensor_ResProd, only: hadronTensor_R
  use particleProperties, only: hadron
  use leptonTensor
  use degRad_conversion
  implicit none
  real, dimension(0:3) :: lin,lout
  real, dimension(0:3) :: q
  real, dimension(0:3) :: pin,pout
  complex,dimension(0:3,0:3) :: hadronTensor,leptonTens
  complex :: matrixElement_Squared
  real :: kinematics, Xsection
  integer :: i
  integer, parameter :: maxSteps=80

  real                , save :: theta_lepton_out=37.1 ! in degrees
  real                , save      :: e_in=0.73
  real                            :: e_out
  integer             , save      :: nukCharge=1
  integer             , intent(in):: resID
  real, dimension(1:3), parameter :: position=(/100.,100.,100./)

  real :: e,de,integral

  call readInput

  pin(0)=mN
  pin(1:3)=0.


  If(.false.) then
     integral=0
     dE=3./float(maxSteps)
     do i=10,maxSteps
        E=dE*float(i)
        pOut=(/E,0.,0.,0./)
        integral=integral+ 2.*E*specFunc(resID,nukCharge,pout,position)
     end do
     integral=integral*dE
     write(111,*) resId, integral
  end If

  e_out_loop: do i=1,maxSteps
     e_out=float(i)*e_in/float(maxSteps)

     ! ********** Kinematics ******************

     ! Initial lepton: Assume lepton in z-direction 
     lin (0)=e_in 
     lin(3)=e_in
     lin(1:2)=0.

     !lepton: Assume phi=0. 
     lout(0)=E_out
     lout(1)=E_out*sin(radian(theta_lepton_out))
     lout(2)=0.
     lout(3)=E_out*cos(radian(theta_lepton_out))

     q=lin-lout
     pout=pin+q

     ! ********** Cross section ******************
     leptonTens=l_munu_matrix(lin,lout)
     if(hadronTensor_R(pin,pout,resID,nukCharge,EM,hadronTensor,hadron(resID)%mass) ) then
        matrixElement_Squared=Contract( leptonTens, hadronTensor)
     else
        matrixElement_Squared=0. 
     end if

     if(abs(AIMAG(matrixElement_Squared)).gt.0.0000001) then
        write(*,*) 'ERROR: (Matrix element)^2  is imaginary:', matrixElement_Squared
     end if


     kinematics=1./SP(lin,pin)*lout(0)/(32.*pi**2)

     Xsection=kinematics  *REAL(matrixElement_Squared)*  specFunc(resID,nukCharge,pout,position)

     ! ********** UNIT CONVERSION ******************

     ! Convert from 1/(sr GeV^3) to mb/(GeV sr)
     ! 1/(GeV^2)=1/(1000 MeV)^2=1/(1000/(197 fm)^2) =0.197**2 *10 mb
     Xsection=Xsection/GeVSquared_times_mb
     
     ! Convert from mb/(GeV sr) to mb/(MeV sr)
     Xsection=Xsection/1000.

!     write(*,*) E_out, Xsection, -SP(q,q)
     write(10+resID,'(5E15.4)') E_out, Xsection, -SP(q,q),specFunc(resID,nukCharge,pout,position), REAL(matrixElement_Squared)
  end do e_out_loop

  contains



    subroutine readInput

      use output

      implicit none
      integer :: ios

      NAMELIST /test/ e_in,nukCharge,theta_lepton_out

      rewind(5)
      read(5,nml=test,IOSTAT=ios)
      call Write_ReadingInput("test",0,ios)

      write(*,*) 'Incoming energy=', e_in
      write(*,*) 'NukCharge=', nukCharge
      write(*,*) 'theta_lepton_out=', theta_lepton_out


      call Write_ReadingInput('test',1)

    end subroutine readInput


    real function gamma_flux(lin,lout)
      ! not yet finished
      use constants, only : alphaQED, pi, mN
      use minkowski
      implicit none
      real, dimension(0:3) :: lin,lout,q,p
      real :: K, epsilon,W,theta

      p=0
      p(0)=mN
      
      q=lin-lout
      
      W=abs4(p+q)

      K=(W**2-mN**2)/mN
      theta=acos(Dot_product(lin(1:3),lout(1:3))/sqrt(Dot_product(lin(1:3),lin(1:3)))/sqrt(Dot_product(lout(1:3),lout(1:3))))

      epsilon=1./(1.+2.*Dot_product(q(1:3),q(1:3))/(-SP(q,q))*(tan(theta/2.))**2)

      Gamma_flux=alphaQED/2/pi**2*lin(0)/lout(0)*K/(-SP(q,q))/(1.-epsilon)
      


    end function gamma_flux

  
end subroutine testing


subroutine testingAll()
  use leptonicID
  use spectralFunc, only : specFunc
  use constants, only : GeVSquared_times_mb, pi, mN
  use minkowski, only : contract,SP
  use IdTable, only : nucleon,nres
  use hadronTensor_ResProd, only: hadronTensor_R
  use particleProperties, only: hadron
  use leptonTensor
  use degRad_conversion
  implicit none
  real, dimension(0:3) :: lin,lout
  real, dimension(0:3) :: q
  real, dimension(0:3) :: pin,pout
  complex,dimension(0:3,0:3) :: hadronTensor,leptonTens
  complex :: matrixElement_Squared
  real :: kinematics, Xsection
  integer :: i
  integer, parameter :: maxSteps=80

  real                , save :: theta_lepton_out=37.1 ! in degrees
  real                , save      :: e_in=0.73
  real                            :: e_out
  integer             , save      :: nukCharge=1
  integer             :: resID
  real, dimension(1:3), parameter :: position=(/100.,100.,100./)

  real :: sigTot

  call readInput

  pin(0)=mN
  pin(1:3)=0.


  e_out_loop: do i=1,maxSteps
     sigTot=0.
     e_out=float(i)*e_in/float(maxSteps)
     res_loop: do resId=nucleon,nucleon+nres


        ! ********** Kinematics ******************

        ! Initial lepton: Assume lepton in z-direction 
        lin (0)=e_in 
        lin(3)=e_in
        lin(1:2)=0.

        !lepton: Assume phi=0. 
        lout(0)=E_out
        lout(1)=E_out*sin(radian(theta_lepton_out))
        lout(2)=0.
        lout(3)=E_out*cos(radian(theta_lepton_out))

        q=lin-lout
        pout=pin+q

        ! ********** Cross section ******************
        leptonTens=l_munu_matrix(lin,lout)
        if(hadronTensor_R(pin,pout,resID,nukCharge,EM,hadronTensor,hadron(resID)%mass) ) then
           matrixElement_Squared=Contract( leptonTens, hadronTensor)
        else
           matrixElement_Squared=0. 
        end if

        if(abs(AIMAG(matrixElement_Squared)).gt.0.0000001) then
           write(*,*) 'ERROR: (Matrix element)^2  is imaginary:', matrixElement_Squared
        end if


        kinematics=1./SP(lin,pin)*lout(0)/(32.*pi**2)

        Xsection=kinematics  *REAL(matrixElement_Squared)*  specFunc(resID,nukCharge,pout,position)

        ! ********** UNIT CONVERSION ******************

        ! Convert from 1/(sr GeV^3) to mb/(GeV sr)
        ! 1/(GeV^2)=1/(1000 MeV)^2=1/(1000/(197 fm)^2) =0.197**2 *10 mb
        Xsection=Xsection/GeVSquared_times_mb

        ! Convert from mb/(GeV sr) to mb/(MeV sr)
        Xsection=Xsection/1000.

        !     write(*,*) E_out, Xsection, -SP(q,q)
        write(10+resID,'(5E15.4)') E_out, Xsection, -SP(q,q),specFunc(resID,nukCharge,pout,position), REAL(matrixElement_Squared)
        sigTot=sigTot+Xsection
     end do res_loop
     write(100,*) E_out,sigTot
  end do e_out_loop

contains



  subroutine readInput

    use output

    implicit none
    integer :: ios

    NAMELIST /test/ e_in,nukCharge,theta_lepton_out

    rewind(5)
    read(5,nml=test,IOSTAT=ios)
    call Write_ReadingInput("test",0,ios)

    write(*,*) 'Incoming energy=', e_in
    write(*,*) 'NukCharge=', nukCharge
    write(*,*) 'theta_lepton_out=', theta_lepton_out


    call Write_ReadingInput('test',1)

  end subroutine readInput

end subroutine testingAll





