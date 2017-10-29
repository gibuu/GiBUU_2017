!*******************************************************************************************************
!****p* test_med/
! NAME
! program test_med
! PURPOSE
! * Basic test of the in-medium X section for e^- n -> e^- n pi
! * Results have to be compared to MAID results for virtual photon cross sections in the vacuum
! * See http://www.kph.uni-mainz.de/MAID/maid2003
!************************************************************************************


program test_med
  use inputGeneral, only: readinputGeneral
  implicit none
  call init_database
  call readinputGeneral

  call test_CM_MonteCARLO

end program test_med


!*******************************************************************************************************
!****s* test_CM_MonteCarlo/
!
! PURPOSE
! Calculates dSigma/dOmega_f/dE_f/dOmega_k in the nucleon rest-frame. 
!
! INPUTS
! * Namelist Xsection : energy_li,energy_lf,theta_lf,theta_k,charge_pionOut,charge_nucOut
!
! OUTPUT
! * File fort.111
!************************************************************************************
subroutine test_CM_MonteCArlo
  use electronPionProd_medium_eN
  use vector, only : absVec, theta_in
  use eN_eventDefinition , only : electronNucleon_event, CM, CALC,write_electronNucleon_event
  use eN_event           , only : init_electronNucleon_event
  use rotation           , only : get_phi_Theta
  use leptonKinematics   , only : buildElectrons
  use particleDefinition
  use particleProperties, only : baryon
  use degRad_conversion  , only : degrees, radian
  use IdTable,only : nucleon,pion  
  use leptonicID, only :EM
  use random
  use minkowski, only : abs4
  use constants, only : pi
  implicit none
  integer :: pionCharge=-1
  real    :: phi=-10.
  real    :: theta=16.
  real    :: energy_li=3.595
  real    :: energy_lf=2.795


  !real    :: theta=20.
  !  real    :: energy_li=4.4
  !  real    :: energy_lf=2.9

  type(electronNucleon_event) :: eN
  type(particle)              :: nuc
  real, dimension(0:3)        :: electron_in,electron_out,k,pf
  real :: phi_k, theta_k, xsection,theta_qk,dOmega,integral_cm,integral_lab,xsection2, integral_gen, integral_nonEN
  integer :: i,j,precision,noSuccess_cm,noSuccess_lab, noSuccess_gen, noSuccess_nonEN
  integer :: num_loops_start=3
  integer :: num_loops,total_Loops
  character(100) :: f

  real, external :: crossSec_nonEN


  f='(A30,E12.5,A13,I4,A15,F9.4,A15, I7)'

  call buildElectrons(electron_in,electron_out,radian(theta),radian(phi),energy_li,energy_lf,(/0.,0.,1./))

  ! Set which creates trouble (taken from Kai's testRun)
 electron_In =(/          2.341807425E+01 ,  1.407124727E+01 ,  0.000000000E+00 ,  1.871913998E+01 /)
 electron_Out=(/          2.098728747E+01 ,  1.407124727E+01 ,  0.000000000E+00 ,  1.557132735E+01/)


  ! Set up resting nucleon
  nuc%ID=nucleon
  nuc%momentum(0)=baryon(nucleon)%mass
  nuc%mass=baryon(nucleon)%mass
  nuc%momentum(1:3)=0.
  nuc%charge=0
  ! (2) Set up electron-nucleon event
  eN=init_electronNucleon_event(electron_in,electron_out,nuc,CALC)
  call write_electronNucleon_event(eN)


  ! (3) Get Pion Production Xsection in LAB
  total_Loops=0

  integral_cm=0.
  integral_lab=0.
  integral_gen=0.
  integral_nonEN=0.
  noSuccess_lab=0
  noSuccess_cm=0
  noSuccess_gen=0
  noSuccess_nonEN=0

  precision_loop: do precision=1,20

     num_loops=num_loops_start*2**precision
     total_loops=total_loops+num_loops
     dOmega=360.*2./float(total_loops)
     dOmega=4*pi/float(total_loops)


     if (.true.) then
        lab_Loop : do j=1,num_loops
           ! Monte-Carlo
           phi_k=rn()*360.
           theta_k=degrees(rnCos())

           Xsection=dSdO_fdE_fdO_k_med_eN(eN, pionCharge, phi_k, theta_k, k, pf, EM,pionNucleonSystem=1)

           ! Check
           !Xsection2=dSdO_fdE_fdO_k_med_eN(eN, pionCharge, phi_k, theta_k, k, pf, EM)
           !if(abs(xsection-xsection2).gt.1E-8) then
           !   write(*,*) xsection-xsection2
           !   stop 'test_cm_monteCarlo'
           !end if
           !write(*,*) 'k=',k
           if(abs(xsection).gt.1E-20) then
              call get_phi_Theta(k,phi, theta)
              theta_qk=acos(Dot_product(en%photon(1:3),k(1:3))/absVec(k(1:3))/absVec(en%photon(1:3)))
              !write(100,*) theta_qk,  Xsection
              integral_lab=integral_lab+Xsection
           else
              noSuccess_lab=noSuccess_lab+1
           end if
        end do lab_Loop

        !write(*,'(A,E15.5,A)') 'No SUCCESS =', noSuccess /(float(num_loops))*100.,'%'
        !write(*,*) 'Integral in LAB=', integral,' Precision=',precision
        !write(*,*) 'Num points=',num_loops
        write(*,f) 'Integral  LAB=', integral_lab*dOmega,' Precision=',precision, 'No SUCCESS =', &
             & noSuccess_lab /(float(total_loops))*100.," %, # points=",total_loops
     end if
     if(.true.) then
        CM_Loop : do j=1,num_loops
           phi_k=rn()*360.
           theta_k=degrees(rnCos())

           Xsection=dSdO_fdE_fdO_k_med_eN(eN, pionCharge, phi_k, theta_k, k, pf, EM,pionNucleonSystem=2)
           if(abs(xsection).gt.1E-20) then
              call get_phi_Theta(k,phi, theta)
              theta_qk=acos(Dot_product(en%photon(1:3),k(1:3))/absVec(k(1:3))/absVec(en%photon(1:3)))
              !write(100,*) theta_qk,  Xsection
              integral_cm=integral_cm+Xsection
           else
              noSuccess_cm=noSuccess_cm+1
           end if
        end do CM_Loop
        !write(*,'(A,E15.5,A)') 'No SUCCESS =', noSuccess /(float(num_loops))," %"
        write(*,f) 'Integral  CM=', integral_cm*dOmega,' Precision=',precision, 'No SUCCESS =', &
             & noSuccess_cm /(float(total_loops))*100.," %, # points=",total_loops
     end if
     if (.false.) then
        generator_Loop : do j=1,num_loops
           ! Monte-Carlo
           !phi_k=rn()*360.
           !theta_k=degrees(rnCos())

           Xsection=crossSec_generator(eN,.true.,pionCharge)
           ! Check
           !Xsection2=dSdO_fdE_fdO_k_med_eN(eN, pionCharge, phi_k, theta_k, k, pf, EM)
           !if(abs(xsection-xsection2).gt.1E-8) then
           !   write(*,*) xsection-xsection2
           !   stop 'test_cm_monteCarlo'
           !end if
           !write(*,*) 'k=',k
           if(abs(xsection).gt.1E-20) then
              !call get_phi_Theta(k,phi, theta)
              !theta_qk=acos(Dot_product(en%photon(1:3),k(1:3))/absVec(k(1:3))/absVec(en%photon(1:3)))
              !write(100,*) theta_qk,  Xsection
              integral_gen=integral_gen+Xsection
           else
              noSuccess_gen=noSuccess_gen+1
           end if
        end do generator_Loop

        !write(*,'(A,E15.5,A)') 'No SUCCESS =', noSuccess /(float(num_loops))*100.,'%'
        !write(*,*) 'Integral in LAB=', integral,' Precision=',precision
        !write(*,*) 'Num points=',num_loops
        write(*,f) 'Integral for Generator=', integral_gen/(float(total_loops)),' Precision=',precision, 'No SUCCESS =', &
             & noSuccess_gen /(float(total_loops))*100.," %, # points=",total_loops
     end if
     if (.true.) then
        nonEN_Loop : do j=1,num_loops
           ! Monte-Carlo
           phi_k=rn()*360.
           theta_k=degrees(rnCos())

           Xsection=crossSec_nonEN(eN,pionCharge,phi_k,theta_k)
           if(abs(xsection).gt.1E-20) then
              !call get_phi_Theta(k,phi, theta)
              !theta_qk=acos(Dot_product(en%photon(1:3),k(1:3))/absVec(k(1:3))/absVec(en%photon(1:3)))
              !write(100,*) theta_qk,  Xsection
              integral_nonEN=integral_nonEN+Xsection
           else
              noSuccess_nonEN=noSuccess_nonEN+1
           end if
        end do nonEN_Loop

        write(*,f) 'Integral for old non-EN=', integral_nonEN*dOmega,&
             & ' Precision=',precision, 'No SUCCESS =', &
             & noSuccess_nonEN /float(total_loops)*100.," %, # points=",total_loops
     end if



     write(*,*)
  end do precision_loop


contains 


  function crossSec_generator(eN,doRes,pionCharge) result(XS)
    use master_1body, only : decayParticle
    use idTable, only : nres,pion,nucleon
    use eventGenerator_eN_lowEnergy, only : eventGen_eN_lowEnergy
    use particleProperties, only : isBaryon
    implicit none
    type(electronNucleon_event), intent(in) :: eN
    logical                    , intent(in) :: doRes
    real                                    :: XS
    integer                    , intent(in) :: pionCharge
    integer                :: channel
    type(particle), dimension(1:10) :: finalState
    real, dimension(1:5) :: XS_Arr
    logical :: success,collisionFlag,pauliFlag
    integer :: i
    logical, dimension(2:nres+1),parameter  :: whichRes=.true.
    type(particle) , dimension(1:1):: resonance
    real :: time


    call eventGen_eN_lowEnergy(eN,.false.,doRes,whichRes,.true.,.false.,.false.,finalState,channel,success,XS,XS_Arr)      
    if(doRes.and.isBaryon(finalState(1)%ID).and.finalState(1)%ID.ne.nucleon) then
       resonance=finalState(1)
       time=0.
       call decayParticle(resonance,finalState,collisionFlag,time,pauliFlag,finalFlag=.true.)
    end if

    do i=lbound(finalState,dim=1), ubound(finalState,dim=1)
       if(finalState(i)%ID.eq.pion.and. finalState(i)%charge.ne.pionCharge) then
          XS=0.
       end if
    end do
  end function crossSec_generator

end subroutine test_CM_MonteCArlo



function crossSec_nonEN(eN,pionCharge,phi_k,theta_k) result(XS)
  use vector, only : theta_in
  use eN_eventDefinition
  !use electronPionProduction_medium, only : dSigmadOmega_fdE_fdOmega_k_med
  use degRad_conversion, only : degrees
  implicit none
  type(electronNucleon_event), intent(in) :: eN
  real                                    :: XS
  integer                    , intent(in) :: pionCharge
  real, intent(in)                        :: phi_k,theta_k

  real :: theta_lf
  real, dimension(0:3) :: q,lf,k,pf

  theta_lf=theta_in(en%electron_in(1:3),en%electron_out(1:3))
  XS=0.! dSigmadOmega_fdE_fdOmega_k_med(en%nucleon,pionCharge,en%electron_in(0),en%electron_out(0), &
    !   theta_lf,phi_k, theta_k,q,lf,k,pf)
  !write(*,*) theta_lf, XS
end function crossSec_nonEN

