program test

use quasiElastic_electron, only : dsigmadt
implicit none


!real, parameter     :: E_in=0.73
real, parameter     :: E_in=1.148
!real, parameter     :: E_in=3.114
integer, parameter  :: nucCharge=1
integer :: theta_degree
real :: Qsquared,dsigmadOmega,dsdt

call init_database()

open(200,File="QE_xsection.dat.1148")
write(200,*) '# Incoming electron energy=',e_in
write(200,*) '# Nucleon charge=',1
write(200,*) '# Qsquared, theta(degree), dsigma/dt, dsigma/dOmega (mb)'
do theta_degree=5,170,1
   dsdt=dSigmadt(E_in,float(theta_degree),nucCharge,QSquared,dsigmadOmega)
   write(200,'(4E12.3)') Qsquared,float(theta_degree), dsdt  ,dsigmadOmega
end do

close(2)
!call hadTensor_qe_test()

contains



subroutine hadTensor_qe_test()
    use minkowski
    use  hadronTensor_QE, only : H_munu_QE
    use leptonTensor    , only : l_munu

    implicit none
    real, dimension(0:3) :: pin,pf,q
    integer              :: nuc_Charge
    integer :: mu,nu,alpha
    real :: matrixElement




    pin(1:3)=(/1.,3.,0.2/)
    pin(0)=sqrt(0.938**2+Dot_product(pin(1:3),pin(1:3)))
    pf(1:3)=(/0.2,1.,0.4/)
    pf(0)=sqrt(0.938**2+Dot_product(pf(1:3),pf(1:3)))
    q=pf-pin
    do nuc_charge=0,1
       matrixElement=0
       do mu=0,3
          do nu=0,3
             do alpha=0,3
                if(nu.ne.alpha) cycle
                matrixElement=matrixElement &
                     & + q(alpha)*metricTensor(alpha,nu)*H_munu_QE(mu,nu,pin,pf,nuc_charge)
             end do
          end do
!          write(99,*) mu, nuc_charge,matrixElement
       end do
    end do



  end subroutine hadTensor_qe_test
end program test
