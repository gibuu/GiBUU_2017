program testMaxwell
  use randomMaxwell

  implicit none

  call check1
  call check2


contains

  subroutine check1
    use histf90

    integer :: ir, nr=1000000
    real :: r
    type(histogram) :: HH

    call CreateHist(HH, 'vals', 0.0, 10.0, 0.1)

    call initMaxwell(0.135, 0.500)

    do ir=1,nr

       r = rnMaxwell()
       call AddHist(HH, r, 1.0)
    end do

    call WriteHist(HH, 114, mul=1.0/nr)

  end subroutine check1

  ! plots the sqrt(s) distribution of two particles
  subroutine check2

    use histf90
    use random
    use fgsl, only: kn => fgsl_sf_bessel_kcn

    real, parameter :: mass = 0.001
!    real, parameter :: mass = 0.135
    real, parameter :: T = 0.150

    integer :: ir, nr=1000000
    real :: p1,p2,srts
    real,dimension(0:3) :: mom1,mom2,momTot
    type(histogram) :: HH

    integer :: im
    real :: m,x


    call CreateHist(HH, 'dN/dsrts', 0.0, 4.0, 0.02)

    call initMaxwell(mass, T)

    do ir=1,nr

       p1 = rnMaxwell()
       p2 = rnMaxwell()

       mom1(1:3) = p1*rnOmega()
       mom2(1:3) = p2*rnOmega()

       mom1(0) = sqrt(p1**2 + mass**2)
       mom2(0) = sqrt(p2**2 + mass**2)

       momTot = mom1+mom2

       srts = sqrt( momTot(0)**2-sum(momTot(1:3)**2) )

       call AddHist(HH, srts, 1.0)
    end do

    call WriteHist(HH, 115, mul=1.0/nr)


    do im=1,500
       m = im * 0.01
       x = m/T
       write(116,*) m, x**3/(32.0*T)*( 2*kn(2,x) + x*kn(1,x) ), &
            x**3/(32.0*T)*( 2*kn(2,x) ), x**3/(32.0*T)*( x*kn(1,x) )
    end do


  end subroutine check2


end program testMaxwell
