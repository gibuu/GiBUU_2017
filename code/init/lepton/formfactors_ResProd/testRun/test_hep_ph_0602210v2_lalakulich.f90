!***************************************************************************************************
!****p* /test_hep_ph_0602210v2_lalakulich
! NAME
! program test_hep_ph_0602210v2_lalakulich
! PURPOSE
! Checks mainly plots of Lalakulich et al. in hep-ph/0602210v2
! and the impact of inserting our form factors in their formulas.
!
! See helcityAmplitudes_test.f90 for a consistency test of our helicity amplitudes.
!***************************************************************************************************

program test_hep_ph_0602210v2_lalakulich

  use particleProperties, only: initParticleProperties
  use leptonicID
  implicit none

  call initParticleProperties
  call testingP11_1440
  call testingD13_1520
  call testingS11_1535
  call testingdelta


contains


subroutine getKine_lab(W,QS,q,pi,pf)
  ! Establishes final and initial nucleon momentum and the momentum exchange q
  ! Assumptions: Initial nucleon at rest, vec(q) in direction of z-axis
  ! Input : W=sqrt(s)
  ! Input : QS=-SP(q,q)
  !
  use constants, only: mN

  real, dimension(0:3) :: q,pi,pf
  real :: W, QS
  pi(1:3)=0.
  pi(0)=mN
  q(0)=(W**2-mN**2+QS)/(2.*mN)
  q(1:2)=0.
  q(3)=sqrt(QS+q(0)**2)
  pf=pi+q
end subroutine getKine_lab


subroutine testingP11_1440
  ! Checks plot 5 of Lalakulich hep-ph/0602210v2 
  use formFactor_ResProd
  use IdTable, only : P11_1440
  use leptonicID,only : EM
  use constants, only : pi, alphaQED, mN
  use particleProperties, only: hadron
  use helicityAmplitudes

  real, dimension(1:8) :: f
  logical :: ff_set
  real :: QS,W,A12,A32,S12,mu,N,mR
  integer :: i,targetCharge
  real, dimension(0:3) :: q,pin,pf
  real :: A32maid1,A12maid1,S12maid1,A32maid2,A12maid2,S12maid2

  mR=hadron(P11_1440)%mass
  mu=mN+mR

  do targetCharge=0,1
     do i=0,100
        W=mR
        QS=float(i)*6./100.

        call getKine_lab(W,QS,q,pin,pf)

        f=getFormfactor_Res(Qs,mR,P11_1440,targetCharge,EM,FF_set)

        N=pi*alphaQED/(mN*(W**2-mN**2))*2*mN*(pf(0)+W)

        A32=0.
        A12=sqrt(N)*sqrt(2.)*q(3)/(pf(0)+W)*(F(1)/mu**2*QS+F(2)/mu*(W+mN))
        S12=sqrt(N)*q(3)**2/(pf(0)+W)*(F(1)/mu**2*(W+mN)-F(2)/mu)


        call get_helicityAmplitudes(targetcharge,P11_1440,Qs,A12maid1,A32maid1,S12maid1,1)  !MAID 2003
        call get_helicityAmplitudes(targetcharge,P11_1440,Qs,A12maid2,A32maid2,S12maid2,2)  !MAID 2005

        write(10+targetCharge,'(10E15.4)') QS,A32,A12,S12,A32maid1,A12maid1,S12maid1,A32maid2,A12maid2,S12maid2
     end do
  end do


end subroutine testingP11_1440


subroutine testingS11_1535
  ! Checks plot 5 of Lalakulich hep-ph/0602210v2 
  use formFactor_ResProd
  use IdTable, only : S11_1535
  use leptonicID,only : EM
  use constants, only : pi, alphaQED, mN
  use particleProperties, only: hadron
  use helicityAmplitudes

  real, dimension(1:8) :: f
  logical :: ff_set
  real :: QS,W,A12,A32,S12,mu,N,mr
  integer :: i,targetCharge
  real, dimension(0:3) :: q,pin,pf
  real :: A32maid1,A12maid1,S12maid1,A32maid2,A12maid2,S12maid2

  mR=hadron(S11_1535)%mass
  mu=mN+mR

  do targetCharge=0,1
     do i=0,100
        W=mR
        QS=float(i)*6./100.

        call getKine_lab(W,QS,q,pin,pf)

        f=getFormfactor_Res(Qs,mR,S11_1535,targetCharge,EM,FF_set)

        N=pi*alphaQED/(mN*(W**2-mN**2))*2*mN*(pf(0)+W)

        A32=0.
        A12=sqrt(2.*N)*(F(1)/mu**2*QS+F(2)/mu*(mR-mN))
        S12=sqrt(N)*q(3)*(-F(1)/mu**2*(mR-mN)+F(2)/mu)


        call get_helicityAmplitudes(targetcharge,S11_1535,Qs,A12maid1,A32maid1,S12maid1,1) !MAID 2003
        call get_helicityAmplitudes(targetcharge,S11_1535,Qs,A12maid2,A32maid2,S12maid2,2)  !MAID 2005

        write(40+targetCharge,'(10E15.4)') QS,A32,A12,S12,A32maid1,A12maid1,S12maid1,A32maid2,A12maid2,S12maid2
     end do
  end do

end subroutine testingS11_1535



subroutine testingD13_1520
  ! Checks plot 2 of Lalakulich hep-ph/0602210v2 
  use formFactor_ResProd
  use IdTable, only : D13_1520
  use leptonicID,only : EM
  use constants, only : pi, alphaQED, mN
  use particleProperties, only: hadron
  use minkowski, only : SP
  use helicityAmplitudes

  real, dimension(1:8) :: f
  logical :: ff_set
  real :: QS,W,A12,A32,S12,mR,N
  real :: A32maid1,A12maid1,S12maid1,A32maid2,A12maid2,S12maid2
  integer :: i,targetCharge
  real, dimension(0:3) :: q,pin,pf

  mR=hadron(D13_1520)%mass

  do targetCharge=0,1
     do i=0,100
        W=mR
        QS=float(i)*6./100.

        call getKine_lab(W,QS,q,pin,pf)

        f=getFormfactor_Res(Qs,mR,D13_1520,targetCharge,EM,FF_set)

        N=pi*alphaQED/(mN*(W**2-mN**2))*2*mN*(pf(0)+W)

        A32=-sqrt(N)*(F(1)/mN*(mR-mN) + F(2)/mN**2*SP(q,pf) +F(3)/mN**2*SP(q,pin))
        A12=-sqrt(N/3.)*(F(1)/mN*(mR-mN-2.*mN/mR*q(3)**2/(pf(0)+mR)) + F(2)/mN**2*SP(q,pf) + F(3)/mN**2 *SP(q,pin))
        S12=-sqrt(2.*N/3.)*q(3)/mR*(-F(1)/mN*mR + F(2)/mN**2*(Qs-2.*mN*q(0)-mN**2) - F(3)/mN*(q(0)+mN))

        call get_helicityAmplitudes(targetcharge,D13_1520,Qs,A12maid1,A32maid1,S12maid1,1) !MAID 2003
        call get_helicityAmplitudes(targetcharge,D13_1520,Qs,A12maid2,A32maid2,S12maid2,2)  !MAID 2005


        write(20+targetCharge,'(10E15.4)') QS,A32,A12,S12,A32maid1,A12maid1,S12maid1,A32maid2,A12maid2,S12maid2

     end do
  end do

end subroutine testingD13_1520


subroutine testingdelta
  ! Checks plot 2 of Lalakulich hep-ph/0602210v2 
  use formFactor_ResProd
  use IdTable, only : delta
  use leptonicID,only : EM
  use constants, only : pi, alphaQED, mN
  use particleProperties, only: hadron
  use minkowski, only : SP
  use helicityAmplitudes

  real, dimension(1:8) :: f
  logical :: ff_set
  real :: QS,W,A12,A32,S12,mR,N
  real :: A32maid1,A12maid1,S12maid1,A32maid2,A12maid2,S12maid2
  integer :: i,targetCharge
  real, dimension(0:3) :: q,pin,pf

  mR=hadron(delta)%mass

  do targetCharge=0,1
     do i=0,100
        W=mR
        QS=float(i)*6./100.

        call getKine_lab(W,QS,q,pin,pf)

        f=getFormfactor_Res(Qs,mR,delta,targetCharge,EM,FF_set)

        N=pi*alphaQED/(mN*(W**2-mN**2))*2*mN*(pf(0)+W)

        A32=(-1.)*(-sqrt(N))*q(3)/(pf(0)+mR)* (F(1)/mN*(mR+mN) + F(2)/mN**2*SP(q,pf) +F(3)/mN**2*SP(q,pin))
        A12=(-1.)*sqrt(N/3.)*(F(1)/mN*(mR+mN-2.*mN/mR*(pf(0)+mR)) + F(2)/mN**2*SP(q,pf) + F(3)/mN**2 *SP(q,pin))*q(3)/(pf(0)+mR)
        S12=(-1.)*sqrt(2.*N/3.)*q(3)**2/(pf(0)+mR)/mR*(F(1)/mN*mR + F(2)/mN**2*W**2 + F(3)/mN*(q(0)+mN))

        call get_helicityAmplitudes(targetcharge,delta,Qs,A12maid1,A32maid1,S12maid1,1) !MAID 2003
        call get_helicityAmplitudes(targetcharge,delta,Qs,A12maid2,A32maid2,S12maid2,2)  !MAID 2005


        write(30+targetCharge,'(10E15.4)') QS,A32,A12,S12,A32maid1,A12maid1,S12maid1,A32maid2,A12maid2,S12maid2

     end do
  end do

end subroutine testingdelta


end program test_hep_ph_0602210v2_lalakulich