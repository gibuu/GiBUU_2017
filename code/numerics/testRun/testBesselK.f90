program testBesselK
  use besselK

  implicit none

!  call test1
!  call test2

  call testAll

contains
  subroutine test1

    use fgsl, only: kn => fgsl_sf_bessel_kcn

    integer :: ix
    real :: x

    do ix=1,500
       x = ix * 0.01
       write(111,*) x, besselK1(x), kn(1,x)
    end do

  end subroutine test1

  subroutine test2

    use fgsl, only: kn => fgsl_sf_bessel_kcn

    integer :: ix
    real :: x

    do ix=1,500
       x = ix * 0.01
       write(112,*) x, besselK2(x), kn(2,x)
    end do

  end subroutine test2

  subroutine testAll

    use fgsl, only: kn => fgsl_sf_bessel_kcn, in => fgsl_sf_bessel_icn

    integer :: ix
    real :: x

    do ix=1,5000
       x = ix * 0.01
       write(110,*) x, &
            besselI0(x), in(0,x), &
            besselI1(x), in(1,x), &
            besselK0(x), kn(0,x), &
            besselK1(x), kn(1,x), &
            besselK2(x), kn(2,x)
    end do

  end subroutine testAll

end program testBesselK
