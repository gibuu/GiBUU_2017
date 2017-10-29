program testRandom
  implicit none

!   call test_rnCos
!   call test_rn_trueFalse

  call test_rnPower(0.)

contains

  subroutine test_rn_trueFalse
    use random, only: rn_trueFalse
    integer,dimension(1:2) :: hits
    integer :: i

    do i=1,10000000
      if(rn_trueFalse()) then
        hits(1)=hits(1)+1
      else
        hits(2)=hits(2)+1
      end if
    end do
    write(*,*) hits
  end subroutine test_rn_trueFalse


  subroutine test_rnCos
    use histf90
    use random, only: rnCos
    type(histogram) :: h,cosh
    integer :: i
    real :: x,y
    integer :: n=1000000

    call createHist(H, "x",-4.,4.,0.01)
    call createHist(cosH, "cos(x)",-4.,4.,0.01)
    do i=1,n
      x=rnCos()
      y=cos(x)
      call AddHist(h, x,1.)
      call AddHist(cosh, y,1.)
    end do
    call writeHist(h,11,0.,1./float(n))
    call writeHist(cosh,12,0.,1./float(n))
  end subroutine test_rnCos


  subroutine test_rnPower (p)
    use random, only: rnPower
    use histf90
    
    real, intent(in) :: p
    
    type(histogram) :: H
    real :: x
    integer :: i
    integer, parameter :: n = 1000000

    call createHist (H, "x", 1., 2., 0.01)
    do i=1,n
      x = rnPower (p, 1., 2.)
      call AddHist (H, x, 1./float(n))
    end do
    call writeHist (H, file="rnPower.dat")
  end subroutine test_rnPower

end program testRandom
