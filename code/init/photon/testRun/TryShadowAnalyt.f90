program TryShadowAnalyt

  real, parameter :: mRho = 0.770, pi=3.14
  real :: nu, Q2
  real :: kV, qV, sigma, alpha
  real :: rho0, A, R, z, zI


  complex :: css, S, hc
  integer :: i1,i2,N

  N = 100

  nu = 10.0
  Q2 = 0.0!1.5

  kV = sqrt(nu**2-mRho**2) ! [GeV]
  qV = sqrt(Q2+nu**2)-kV   ! [GeV]

  sigma = 20.8*(1+0.766/sqrt(kV)) ! [mb]
  alpha = -0.7667/(sqrt(kV)+0.766)! [1]

  A = 4.0
  R = 2.12 ! [fm]

!  A = 208. !4.0
!  R = 7.915! 2.12 ! [fm]


  rho0 = 3*A/(4*pi*R**3) ! [fm^-3]

  css = sigma*cmplx(1,-alpha)/(2*10.) ! [fm^2] 


  write(*,*) '-----'
  write(*,*) kV
  write(*,*) qV
  write(*,*) sigma
  write(*,*) alpha
  write(*,*) css
  write(*,*) rho0
  write(*,*) '-----'
  

  S = cmplx(0.,0.)
  do i1=-N,N
     z = i1*R/N
     S = cmplx(0.,0.)
     do i2=-N,i1
        zI=i2*R/N
        hc = rho0*css*exp(cmplx(0,qV*(zI-z)/0.197))*exp(-css*II(zI,z))
        S = S + hc * R/N
     end do
!     write(*,*) z,S,1.-S
     write(*,*) z,1.-S,(abs(1.-S))**2
     write(192,*) z,(abs(1.-S))**2
  enddo

  

contains
  real function II(zI,z)
    implicit none
    real :: zI,z

    II = 0.
    if (zI.gt.R) return
    if (z.lt.-R) return
    II = rho0*(min(z,R)-max(zI,-R))
    return

  end function II



end program TryShadowAnalyt
