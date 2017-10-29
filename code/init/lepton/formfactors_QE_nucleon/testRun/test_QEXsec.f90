program test
  use constants, only : pi
  use particleProperties, only: initParticleProperties
  implicit none
  real :: dsigma_dOmega,e_in,theta,QSquared
  integer :: i

  call initParticleProperties

  ! Evaluates the cross section for quasielastic scattering off a free nucleon
  e_in=3.114
  open(10,file='QE_dOmega_xsection.dat')
  do i=0,180
     theta=pi/180.*float(i)
     call xsec(e_in,theta,dsigma_dOmega, QSquared)
     write(10,'(5E14.4)') e_in, theta,float(i), QSquared, dsigma_dOmega
  end do

contains

  subroutine xsec(e_in,theta,dsigma_dOmega,QSquared)
    use FF_QE_nucleonScattering
    use constants, only: alphaQED, mN

    real, intent(in) :: e_in, theta
    real,intent(out) :: dsigma_dOmega, Qsquared

    real :: F1,F2,dummy,e_out,dsigma_dOmega_Mott

    e_out=    e_in/(1.+2*e_in/mN*sin(theta/2.)**2)
    QSquared= 4.*e_in*e_out*sin(theta/2.)**2

    dsigma_dOmega_Mott=alphaQED**2 /(4.*e_in**2*sin(theta/2.)**4)    *e_out/e_in  *cos(theta/2.)**2

    call formfactors_QE(QSquared,2,1,F1,F2,dummy,dummy) 

    dsigma_dOmega=dsigma_dOmega_Mott*(F1**2+QSquared/4./mN**2*(F2**2+2.*(F1+F2)**2*tan(theta/2.)**2))

    ! Converting to units 
    ! 1/GeV**2=1/1000**2/MeV**2=1/1000**2/(1/197 fm**2)=(197/1000)**2 fm**2= (197/1000)**2 * 10 mb

    dsigma_dOmega=dsigma_dOmega*0.197**2*10.

  end subroutine xsec


end program test
