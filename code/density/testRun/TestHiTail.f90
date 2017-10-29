program TestHiTail
  use random
  use histf90


  implicit none

  integer :: i
  integer, parameter :: n = 1000000
  type(histogram) :: hhh
  real :: x
  real, dimension(1:3) :: p


  call CreateHist(hhh,"HiTail distribution", 0.0,5.0,0.02)

  do i=1,n
     x = chooseAbsMomentum()

!!$     do   ! Monte Carlo distribution of momentum in sphere of radius pFermi 
!!$        p(1)=1.-2.*rn(0) 
!!$        p(2)=1.-2.*rn(0) 
!!$        p(3)=1.-2.*rn(0) 
!!$        x = sqrt(dot_product(p,p))
!!$        if (x.le.1.) exit 
!!$     end do


     call AddHist(hhh,x,1.0)
  enddo

  call WriteHist(hhh,101,mul = 1.0/n, add=1e-20)



contains
  real function chooseAbsMomentum()
    implicit none
    
!    logical, parameter :: HiTail = .false.
    logical, parameter :: HiTail = .true.

    real, parameter :: v1 = 0.85  ! n(p) for p=0..1
    real, parameter :: v2 = 0.15  ! n(p) for p==1.000001
!    real, parameter :: A = 3.1784 ! slope for exp(-A p)
    real, parameter :: A = 2.3    ! slope for exp(-A p)/p^2

    real, parameter :: h1 = v1/3, h2 = h1 + v2/A
    real :: h,u

    if (.not.HiTail) then
       chooseAbsMomentum = rn(0)**(1./3.)
    else
       h = v2*exp(A) ! scaling of hiTail
       u = rn(0)

       chooseAbsMomentum = (h2/h1*u)**(1./3.)

       if (chooseAbsMomentum.gt.1.0) then
          chooseAbsMomentum = -log(h2*A/h*(1-u))/A
       endif
    endif
  end function chooseAbsMomentum
  



end program TestHiTail
