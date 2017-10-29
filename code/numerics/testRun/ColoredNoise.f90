program ColoredNoise

  ! This program generates a list of random numbers as a function of time,
  ! where the correation function is given by an exponential.
  !
  ! Taken from
  ! J.Schmidt, A.Meistrenko, H.van Hees, Z.Xu, C.Greiner
  ! Phys Rev E 91 (2015) 032125
  !
  ! The parameter M could be reduced, since the region deltaT*[-M...M] only
  ! has to cover the range, where G(t) does not vanish
  !
  ! please note: the fourier transform of a real, even function is
  ! also real and even, thus we do not need to use complex numbers

  implicit none

  real, parameter :: correlD = 1.2
  real, parameter :: correlTau = 7.5
  integer, parameter :: nFile = 100 ! number of files to generate

  call generate

contains

  real function correl(t)
    real, intent(in) :: t
    correl = correlD/(2*correlTau)*exp(-abs(t)/correlTau)
  end function correl

  subroutine generate

    use constants, only: pi
    use random, only: rnGauss
    use output, only: intToChar4
!    use fgsl, only: kn => fgsl_sf_bessel_kcn

    real, parameter :: deltaT = 0.1
    integer, parameter :: M = 10000
    integer, parameter :: N = 50000


    real, dimension(-m:m) :: S, G
    real, dimension(-m:m+n) :: xiw
!    real, dimension(0:n) :: xi


    integer :: i,j, iFile
    real :: omegan,tm,rh
!!$    complex :: ch

    write(*,*) '=== eq. (A2) ==='

    do i=-M,M
       omegan = 2*pi*i/((2*M+1)*deltaT)

!!$       ch = (0.,0.)
!!$       do j=-M,M
!!$          tm = j*deltaT
!!$          ch = ch + CMPLX(cos(omegan*tm),sin(omegan*tm))*correl(tm)
!!$       end do
!!$       ch = ch*deltaT
!!$
!!$!       rh = correlD/(1+(correlTau*omegan)**2)
!!$!       write(*,*) i,ch,rh,REAL(ch)/rh
!!$
!!$       S(i) = REAL(ch)

       rh = 0.
       do j=-M,M
          tm = j*deltaT
          rh = rh + cos(omegan*tm)*correl(tm)
       end do
       rh = rh*deltaT

!       rh = correlD/(1+(correlTau*omegan)**2)
!       write(*,*) i,ch,rh,REAL(ch)/rh

       S(i) = rh
    end do


    write(*,*) '=== eq. (A3) ==='

    do j=-M,M
       tm = j*deltaT

!!$       ch = (0.,0.)
!!$       do i=-M,M
!!$          omegan = 2*pi*i/((2*M+1)*deltaT)
!!$
!!$          ch = ch + CMPLX(cos(omegan*tm),-sin(omegan*tm))*sqrt(S(i))
!!$       end do
!!$       ch = ch/((2*M+1)*deltaT)
!!$
!!$       !       write(*,*) j, ch
!!$       G(j) = REAL(ch)

       rh = 0.
       do i=-M,M
          omegan = 2*pi*i/((2*M+1)*deltaT)

          rh = rh + cos(omegan*tm)*sqrt(S(i))
       end do
       rh = rh/((2*M+1)*deltaT)

       !       write(*,*) j, ch
       G(j) = rh

!       if (j.ne.0) then
!          rh = sqrt(correlD)/(pi*correlTau)*kn(0,abs(tm)/correlTau)
!       else
!          rh = 0.0
!       end if
!       write(112,*) j,tm,G(j),rh
    end do

    do iFile=1,nFile
       write(*,*) 'iFile = ',iFile
       write(*,*) '=== eq. (A4) ==='

       do j=-M,M+N
          xiw(j) = rngauss(1.0,0.0)/sqrt(deltaT)
          !       write(*,*) j,xiw(j)
       end do

       write(*,*) '=== eq. (A5) ==='

       open(113,file="noise."//intToChar4(iFile-1)//".dat", status="unknown")

       do j=0,N
          tm = j*deltaT
          rh = 0
          do i=-M,M
             rh = rh + G(i)*xiw(i+j)
          end do
          rh = rh*deltaT

          !       xi(j) = rh

          write(113,*) j, rh
       end do

       close(113)
    end do

  end subroutine generate

end program ColoredNoise
