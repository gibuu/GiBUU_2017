! Evaluates average Fermi Momentum for an Oxygen nucleus
program test
  call fold_centers_to_matter()
  !call average_p_f()
end program test

subroutine average_p_f
  use densityModule, only : fermiMomentum_sym
  use densityStatic, only : densityLuis
  use constants, only : pi
  implicit none
  real, parameter :: dr=0.01
  integer i
  real :: integral, norm,r
  real :: weight,rhop_mat,rhon_mat,rho,rhop_cent,rhon_cent

  integral=0.
  norm=0.
  do i=0,10000
     r=float(i)*dr
     call densityLuis(r,6,12,rhop_mat,rhon_mat)
     rho=rhop_mat+rhon_mat
     weight=rho * r**2*dr
     integral=integral+ fermiMomentum_sym(rho) *weight
     norm=norm+weight
  end do
  write(*,'(5F16.5)') norm, integral, norm*4.*pi, integral*4.*pi, integral/norm


end subroutine average_p_f


subroutine fold_centers_to_matter
  use densityModule, only : fermiMomentum_sym
  use densityStatic, only : densityLuis
  use constants, only : pi
  implicit none
  real, parameter :: dr=0.01
  integer, parameter :: numPoints_phi=100
  integer, parameter :: numPoints_theta=100
  real :: dPhi, dTheta
  integer :: i,j,k
  real :: integral, norm,rPrime,r
  real :: weight,rhop_mat,rhon_mat,rho,rhop_cent,rhon_cent

  ! Plot Luis density profile for both density of matter and centers
  open(101,file='DensLuis.dat')
  write(101,'(5A12)') '#r','rhop_mat','rhon_mat','rhop_cent','rhon_cent'
  do i=0,100
     r=float(i)*0.1
     call densityLuis(r,6,12,rhop_mat,rhon_mat)
     call densityLuis(r,6,12,rhop_cent,rhon_cent,center_in=.true.)
     write(101,'(5E12.4)') r,rhop_mat,rhon_mat,rhop_cent,rhon_cent
  end do
  close(101)


  open(101,file='DensLuis_folding.dat')
  do j=1,100
     r=float(j)*0.1
     integral=0.
     norm=0.
     do i=0,1000
        rPrime=float(i)*dr
        call densityLuis(rPrime,6,12,rhop_cent,rhon_cent,center_in=.true.)
        weight=rPrime*dr*(gauss(r+rPrime)-gauss(r-rPrime))
        integral=integral+ rhop_cent * weight
        norm    =norm    + weight
     end do
     write(101,'(4E12.4)') r,integral,norm,integral/norm
     write(*,'(4E12.4)') r,integral,norm,integral/norm
  end do
  close(101)
  ! Folding
  contains
    real function gauss(x)
      implicit none
      real, intent(in) :: x
      real, parameter :: width=0.5
      gauss=exp(-x**2/2./width**2)
    end function gauss
end subroutine fold_centers_to_matter
