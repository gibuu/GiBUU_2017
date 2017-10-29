!******************************************************************************
!****m* /hadronTensor_Npi
! NAME
! module hadronTensor_Npi
! PURPOSE
! * Evaluates the hadron tensor for gamma N -> N pi
! * For details see the notes about this in the work of Oliver Buss
!******************************************************************************

module hadronTensor_NPi


contains

  complex function cMult(a,b)
    ! Evaluates  a times b*
    implicit none
    complex :: a,b
    real :: r,i
    r=Real(a)*Real(b)+AImag(a)*AImag(b)
    i=Real(b)*AImag(a)-Real(a)*AImag(b)
    cMult=Cmplx(r,i)
  end function cMult

  complex function H_munu(mu,nu,pin,pout,k,q,A_maid)
    use minkowski, only: metricTensor,SP
    implicit none
    integer,intent(in) :: nu, mu                ! Lorentz indices
    real, intent(in),dimension(0:3) :: pin,pout ! incoming and outgoing nucleon momentum
    real, intent(in),dimension(0:3) :: q        ! incoming virtual photon momentum
    real, intent(in),dimension(0:3) :: k        ! outgoing pion momentum
    complex, intent(in),dimension (1:6) :: A_maid    ! Form factors

    real, dimension(0:3) :: pave
    real :: mf,mi
    !complex ::H_munu_complex
    complex     :: a
    integer :: alpha, beta


    if ((max(mu,nu).gt.3).or.(min(mu,nu).lt.0)) then
       write(*,*) 'Error in hadronic tensor: indices are not well defined'
       write(*,*) mu,nu
    end if

    if (Dot_Product(pin+q-pout-k,pin+q-pout-k).gt.0.0001) then
       write(*,*) 'Error in hadronic tensor: momentum is not conserved'
       write(*,*) pin
       write(*,*) '+', q,  '=>'
       write(*,*) pout
       write(*,*) '+', k
    end if

    mi=sqrt(SP(pin,pin))
    mf=sqrt(SP(pout,pout))

    PAve=(pout+pin)/2.

    a=-A_Maid(1)+2*mi*A_maid(4)

    H_munu=cMult(a,a)*s1(mu,nu)+cmult(a,c(nu))*s3(mu)+cmult(c(mu),a)*s3(nu)+cmult(c(nu),c(mu))*s9()

    ! Single beta contractions
    do beta=0,3
       H_munu=H_munu+cMult(a,b(nu,beta))*s2(mu,beta) &
            & +cMult(c(mu),b(nu,beta))*s6(beta) &
            & +cMult(b(mu,beta),a)*s2(nu,beta)+cMult(b(mu,beta),c(nu))*s6(beta)
    end do

    ! beta and alpha contractions
    do beta=0,3
       do alpha=0,3
          H_munu=H_munu+cMult(b(mu,alpha),b(nu,beta))*s5(alpha,beta)
       end do
    end do
    H_munu=-H_munu/(8.*mi*mf)

  contains
    real function s1(mu,nu)
      implicit none
      integer, intent(in) :: mu,nu
      s1=4.*SP(q,q)*(  pout(mu)*pin(nu)+pout(nu)*pin(mu)-metricTensor(mu,nu)*SP(pout,pin)  ) &
           & -8.*SP(pin,q)*(  pout(mu)*q(nu)+pout(nu)*q(mu)-metricTensor(mu,nu)*SP(pout,q)    ) &
           & +4.* mi*mf*SP(q,q)*MetricTensor(mu,nu)
    end function s1

    !**************************************************************************

    real function s3(mu)
      implicit none
      integer, intent(in) :: mu
      s3=-4.*(pout(mu)*SP(pin,q)+q(mu)*SP(pin,pout)-pin(mu)*SP(q,pout)) + 4*mf*mi*q(mu)
    end function s3

    !**************************************************************************

    real function s2(mu,beta)
      implicit none
      integer, intent(in) :: mu,beta
      s2=-4.*mi*(pout(mu)*q(beta)+q(mu)*pout(beta)-metricTensor(mu,beta)*SP(pout,q)) &
           & +4.*mf*(-pin(mu)*q(beta)+q(mu)*pin(beta)+metricTensor(mu,beta)*SP(pin,q))
    end function s2

    !**************************************************************************

    real function s9()
      implicit none
      s9=4.*mf*mi-4.*SP(pout,pin)
    end function s9

    !**************************************************************************

    real function s6(mu)
      implicit none
      integer, intent(in) :: mu
      s6=-4.*mi*pout(mu)+4.*mf*pin(mu)
    end function s6

    !**************************************************************************

    real function s5(alpha,beta)
      implicit none
      integer, intent(in) :: alpha,beta
      s5=4.*mi*mf*metricTensor(alpha,beta)-4.*(pout(alpha)*pin(beta)+pout(beta)*pin(alpha)&
           &   -SP(pout,pin)*metricTensor(alpha,beta))
    end function s5

    !**************************************************************************


    complex function  c(mu)
      implicit none
      integer, intent(in) :: mu
      c=(A_maid(1)-2*mi*A_maid(4))*q(mu) &
           & +2*A_maid(2)*(pave(mu)*SP(q,k-q/2.)-SP(pave,q)*(k(mu)-q(mu)/2.)) &
           & +A_maid(5)*(q(mu)*SP(k,q)-SP(q,q)*k(mu))
    end function c

    !**************************************************************************

    complex function  b(mu,alpha)
      implicit none
      integer, intent(in) :: mu,alpha
      complex :: term
      b=kronecker_delta(mu,alpha)*(-A_maid(3)*SP(k,q)-2*A_maid(4)*SP(q,pave)+A_maid(6)*SP(q,q))
      term=A_maid(3)*k(mu)+2*A_maid(4)*pave(mu)-A_maid(6)*q(mu)
      if (alpha.eq.0) then
         b=b+q(alpha)*term
      else
         b=b-q(alpha)*term
      end if
    end function b

  end function H_munu

  !****************************************************************************
  !****************************************************************************

  real function kronecker_delta(i,j)
    integer, intent(in) :: i,j
    if (i.eq.j) then
       kronecker_delta=1.
    else
       kronecker_delta=0.
    end if
  end function kronecker_delta

end module hadronTensor_NPi
