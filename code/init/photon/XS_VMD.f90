!******************************************************************************
!****m* /XS_VMD
! NAME
! module XS_VMD
!
! PURPOSE
! This module contains the routines for cross section calculations
! of VMD and GVMD processes.
!
! INPUTS
! (none -- see NOTES)
!
! NOTES
! * The routine vmd uses the shadowing scaling by using PARP(161)..PARP(164)
! * The routines vmd and gvmd return 0 for a cross section, if srts<0.938+m_V.
!   This could (and should) be improved.
!******************************************************************************
module XS_VMD

  implicit none
  private
  public :: vmd, gvmd !, eleformf

  real,dimension(4),parameter :: mV  = (/ 0.76850, 0.78194, 1.01940, 3.09688 /) ! Masses of Mesons

!...internal parameters:

  real, parameter :: alpha = 7.29927e-3
  real, parameter :: epsilon = 0.0808
  real, parameter :: eta = -0.4525

  real,parameter :: a=0.5, N=3.0

contains

  !****************************************************************************
  !****s* XS_VMD/vmd
  ! NAME
  ! subroutine vmd(srts,Q2,eps,sigma,useVM)
  ! PURPOSE
  ! calculate the VMD cross sections for (rho,omega,phi,J/psi)
  ! (using the PYTHIA parameters for scaleVMD)
  !
  ! INPUTS
  ! * real :: srts -- the W value
  ! * real :: Q2   -- the Q2 value
  ! * real :: eps  -- the epsilon value of the photon
  ! * logical, dimension(4), OPTIONAL :: useVM -- true, if the corresponding
  !   meson should be considered; i.e. if false, the XS of the VM is set to 0
  !
  ! OUTPUT
  ! * real, dimension(0:4) :: sigma -- the XS
  !
  ! NOTES
  ! * Attention, we respect the vacuum thresholds, i.e. W>0.938+m_V.
  !   This should be improved!
  ! * Ref: C.Friberg, T.Sjöstrand, JHEP 09 (2000) 010
  !        sigma^{Vp}       : eq. (2.15)
  !        sigma^{\gamma^*p}: eq. (2.24) (for eps=0)
  !                         : eq. (2.28) (for eps=1)
  ! * uses the interpolation (1+eps*r_i) between eq.(2.24) and (2.28)
  !   with r_1: eq.(2.26) [=default] or r_2: eq.(2.27) as set in
  !   MSTP(17)=4,5 in /PYPARS/.
  ! * CfV(4) are the coupling constants [f_V^2/4pi] encoded in PARP(160+i).
  ! * PYTHIA uses mV == mrho. Therefore we have the second column in the
  !   output array XS. If one wants to compare the output of this routine
  !   with the Monte Carlo output, then one realizes:
  !      XS(1) from CollectXS_class() <=> XS(0,2)*(W2/(W2+Q2))**3 from here
  !****************************************************************************
  subroutine vmd(srts,Q2,eps,sigma,useVM)
    use constants, only: mN

    real :: srts,Q2,eps
    real,intent(out):: sigma(0:4)
    logical, dimension(4), intent(IN),OPTIONAL :: useVM

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    real,dimension(4),parameter :: CXSa = (/ 13.63, 13.63, 10.01, 0.970 /)
    real,dimension(4),parameter :: CXSb = (/ 31.79, 31.79, -1.52,-0.146 /)
!     real,dimension(4),parameter :: CmV  = (/ 0.76850, 0.78194, 1.01940, 3.09688 /) ! Masses of Mesons

    real XS_Vp(4),test,XS(0:4,2)
    real W2, Dipole(4), LongEnh(4), PARP165
    integer i
    real CfV(4)
    integer MSTP17

    data PARP165 /0.5/

!...get coupling constants and MSTI17:

    do i=1,4
       CfV(i) = PARP(160+i)   ! f_V^2/4pi
    end do
    MSTP17 = MSTP(17)

!...calculate sigma^{Vp}(s) (and slope parameter)

    W2 = srts**2
    do i=1,4
       XS_Vp(i) = CXSa(i)*W2**epsilon+CXSb(i)*W2**eta
       if (i>1 .and. srts<mN+mV(i)) XS_Vp(i) = 0.0 ! vacuum thresholds
       if (present(useVM).and.(.not.(useVM(i))))  XS_Vp(i) = 0.0
    end do

!...calculate the Dipole factor, the longitudal photon enh. factor

    do i=1,4
       Dipole(i) = (mV(i)**2/(mV(i)**2+Q2))**2
       if (MSTP17.eq.4) then
          LongEnh(i) = mV(i)**2*Q2/((mV(i)**2+Q2))**2
       else if (MSTP17.eq.5) then
          LongEnh(i) = Q2/(mV(i)**2+Q2)
       else if (MSTP17.eq.6) then ! only possible with HERMES parameters
          LongEnh(i) = 0.25* (Q2/mV(i)**2)**0.61
       else
          write(*,*) 'VMDParam: MSTP(17) <> (4,5,6). stop'
          stop
       end if
       LongEnh(i) = 1. + eps * LongEnh(i) * 4.*PARP165
    end do

!...calculate the sum (with Q2 dependence, different mV--handling):
    XS(0,1) = 0.
    XS(0,2) = 0.
    do i=1,4
       XS(i,1) = alpha*XS_VP(i)/CfV(i) * Dipole(i)*LongEnh(i) ! m_V
       XS(i,2) = alpha*XS_VP(i)/CfV(i) * Dipole(1)*LongEnh(1) ! m_rho

       XS(0,1) = XS(0,1) + XS(i,1)
       XS(0,2) = XS(0,2) + XS(i,2)
    end do

!...calculate VMD cross sections with PYTHIA formfactor
    sigma(0)=XS(0,2)*(W2/(W2+Q2))**3*1000. !mub
    test=0.
    do i=1,4
       sigma(i)=XS(i,2)*(W2/(W2+Q2))**3*1000. !mub
       test=test+sigma(i)
    end do

    if (abs(test-sigma(0)).gt.1.e-5) then
       write(*,*) 'problems in vmd',sigma(0),test
       stop
    end if

    return
  end subroutine vmd


  !****************************************************************************
  !****s* XS_VMD/gvmd
  ! NAME
  ! subroutine gvmd(srts,Q2,eps,siggvmd)
  ! PURPOSE
  ! calculate the GVMD cross sections
  !****************************************************************************
  subroutine gvmd(srts,Q2,eps,siggvmd,useVM)
    use parBarMes_HighEnergy
    use constants, only: alphaQED,pi, mN

    real,               intent(in)  :: srts,Q2,eps
    real,dimension(0:4),intent(out) :: siggvmd
    logical, dimension(4), intent(IN),OPTIONAL :: useVM

    real,dimension(4)               :: sigvec
    integer :: i
    real :: ea,formfac,kv,sum
    real,dimension(0:1) :: sum1,k
    real,dimension(4) :: eq2

    eq2(1)=0.5*(5./9.) !(u+d)/2
    eq2(2)=0.5*(5./9.) !(u-d)/2
    eq2(3)=1./9. !s
    eq2(4)=4./9. !c

    k(0)=0.5
    k(1)=1.9*(srts**2/1e06)**epsilon
    kv=0.4

    ea=eps*a
    kv=kv**2

    formfac=(srts**2/(Q2+srts**2))**3

    do i=0,1
       k(i)=k(i)**2
       sum1(i)=((3+2*ea)*Q2**2+24*(1.+ea)*Q2*k(i)+48*k(i)**2)/(Q2+4*k(i))**3
    end do

    sum=alphaQED/pi*4./3.*kv*(sum1(0)-sum1(1))*formfac
    call paramBarMesHE_v(srts,sigvec)

    do i=1,4 !...respect vacuum thresholds (attention: vacuum !!!!)
       if (i>1 .and. srts<mN+mV(i)) sigvec(i) = 0.0 ! vacuum thresholds
       if (present(useVM).and.(.not.(useVM(i))))  sigvec(i) = 0.0
    end do

    siggvmd(0)=0.
    do i=1,4
       siggvmd(i)=eq2(i)*sum*sigvec(i)*1000. !mub
       siggvmd(0)=siggvmd(0)+siggvmd(i)
    end do

  end subroutine gvmd


  !****************************************************************************
  !****s* XS_VMD/Eleformf
  ! NAME
  ! subroutine eleformf(srts,Q2,eps,formfac)
  ! PURPOSE
  ! calculate a Formfactor
  !****************************************************************************
!   subroutine eleformf(srts,Q2,eps,formfac)
!     real,             intent(in)  :: srts,Q2,eps
!     real,dimension(4),intent(out) :: formfac
!
!     real :: s
!
!     s=srts**2
!     formfac = (1.0 + eps * a*4*mv**2*Q2/(mv**2+Q2)**2) * (mv**2/(Q2+mv**2))**2 * (s/(Q2+s))**n
!
!   end subroutine eleformf



end module XS_VMD
