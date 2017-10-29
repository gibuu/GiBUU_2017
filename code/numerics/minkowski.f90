!******************************************************************************
!****m* /minkowski
! NAME
! module minkowski
!
! PURPOSE
! This module defines functions which are connected to Relativity:
! Metric Tensor, Scalar Product, Gamma matrices, ...
!
! NOTES
! * Uses the "mostly -" metric.
!******************************************************************************
module minkowski
  use matrix_module, only: unit4
  use constants, only: ii

  implicit none
  private

  ! ATTENTION:
  ! The reshape command generates 'transposed' arrays: e.g:
  ! A(2,2) = reshape((/1,2, &
  !                    3,4/),(/2,2/)
  ! leads to
  ! A(1,:) --> 1, 3
  ! A(2,:) --> 2, 4


  !****************************************************************************
  ! the Pauli matrices:

  complex, dimension(1:2,1:2), public, parameter :: sigma1 =    reshape((/0., 1., 1., 0./),(/2,2/))
  complex, dimension(1:2,1:2), public, parameter :: sigma2 = ii*reshape((/0., 1.,-1., 0./),(/2,2/))
  complex, dimension(1:2,1:2), public, parameter :: sigma3 =    reshape((/1., 0., 0.,-1./),(/2,2/))

  !****************************************************************************

  complex, dimension(0:3,0:3), parameter :: zero4 = 0.

  complex, dimension(0:3,0:3), public, parameter :: gamma0 = reshape((/ 1., 0., 0., 0., &
                                                                        0., 1., 0., 0., &
                                                                        0., 0.,-1., 0., &
                                                                        0., 0., 0.,-1. /),(/4,4/))

  complex, dimension(0:3,0:3), public, parameter :: gamma1 = reshape((/ 0., 0., 0.,-1., &
                                                                        0., 0.,-1., 0., &
                                                                        0., 1., 0., 0., &
                                                                        1., 0., 0., 0. /),(/4,4/))

  complex, dimension(0:3,0:3), public, parameter :: gamma2 = ii*reshape((/ 0., 0., 0.,-1., &
                                                                           0., 0., 1., 0., &
                                                                           0., 1., 0., 0., &
                                                                          -1., 0., 0., 0. /),(/4,4/))

  complex, dimension(0:3,0:3), public, parameter :: gamma3 = reshape((/ 0., 0.,-1., 0., &
                                                                        0., 0., 0., 1., &
                                                                        1., 0., 0., 0., &
                                                                        0.,-1., 0., 0. /),(/4,4/))

  complex, dimension(0:3,0:3), public, parameter :: gamma5 = reshape((/ 0., 0., 1., 0., &
                                                                        0., 0., 0., 1., &
                                                                        1., 0., 0., 0., &
                                                                        0., 1., 0., 0. /),(/4,4/))

  complex, dimension(0:3,0:3), parameter :: gamma6 = (unit4-gamma5)/2.
  complex, dimension(0:3,0:3), parameter :: gamma7 = (unit4+gamma5)/2.

  complex, dimension(0:3,0:3), parameter :: gamma8 = reshape((/ 0., 0.,-1., 0., &
                                                                0., 0., 0.,-1., &
                                                                1., 0., 0., 0., &
                                                                0., 1., 0., 0. /),(/4,4/))

  complex, dimension(0:3,0:3), parameter :: gamma9 = reshape((/ 0., 1., 0., 0., &
                                                                1., 0., 0., 0., &
                                                                0., 0., 0.,-1., &
                                                                0., 0.,-1., 0. /),(/4,4/))

  complex, dimension(0:3,0:3), parameter :: gamma10 = ii*reshape((/ 0., 1., 0., 0., &
                                                                   -1., 0., 0., 0., &
                                                                    0., 0., 0.,-1., &
                                                                    0., 0., 1., 0. /),(/4,4/))

  complex, dimension(0:3,0:3), parameter :: gamma11 = reshape((/ 1., 0., 0., 0., &
                                                                 0.,-1., 0., 0., &
                                                                 0., 0.,-1., 0., &
                                                                 0., 0., 0., 1. /),(/4,4/))

  !****************************************************************************
  !****g* minkowski/gamma
  ! NAME
  ! complex, dimension(0:3,0:3,0:11), parameter :: gamma
  ! PURPOSE
  ! Represents the gamma matrices gamma0-gamma3, gamma5, gamma6-gamma11.
  !****************************************************************************
  complex, dimension(0:3,0:3,0:11), public, parameter :: gamma = reshape( (/gamma0,gamma1,gamma2,gamma3,zero4,gamma5, &
                                                                            gamma6,gamma7,gamma8,gamma9,gamma10,gamma11/), &
                                                                          (/4,4,12/))


  !****************************************************************************
  !****************************************************************************

  complex, dimension(0:3,0:3), parameter :: sigma4_01 = ii* reshape((/ 0., 0., 0., 1., &
                                                                       0., 0., 1., 0., &
                                                                       0., 1., 0., 0., &
                                                                       1., 0., 0., 0. /),(/4,4/))
  complex, dimension(0:3,0:3), parameter :: sigma4_02 =     reshape((/ 0., 0., 0.,-1., &
                                                                       0., 0., 1., 0., &
                                                                       0.,-1., 0., 0., &
                                                                       1., 0., 0., 0. /),(/4,4/))
  complex, dimension(0:3,0:3), parameter :: sigma4_03 = ii* reshape((/ 0., 0., 1., 0., &
                                                                       0., 0., 0.,-1., &
                                                                       1., 0., 0., 0., &
                                                                       0.,-1., 0., 0. /),(/4,4/))
  complex, dimension(0:3,0:3), parameter :: sigma4_12 =     reshape((/ 1., 0., 0., 0., &
                                                                       0.,-1., 0., 0., &
                                                                       0., 0., 1., 0., &
                                                                       0., 0., 0.,-1. /),(/4,4/))
  complex, dimension(0:3,0:3), parameter :: sigma4_13 = ii* reshape((/ 0.,-1., 0., 0., &
                                                                       1., 0., 0., 0., &
                                                                       0., 0., 0.,-1., &
                                                                       0., 0., 1., 0. /),(/4,4/))
  complex, dimension(0:3,0:3), parameter :: sigma4_23 =     reshape((/ 0., 1., 0., 0., &
                                                                       1., 0., 0., 0., &
                                                                       0., 0., 0., 1., &
                                                                       0., 0., 1., 0. /),(/4,4/))

  complex, dimension(0:3,0:3,0:3,0:3), public, parameter :: sigma4par =     &
       reshape((/     zero4,-sigma4_01,-sigma4_02,-sigma4_03, &
                  sigma4_01,     zero4,-sigma4_12,-sigma4_13, &
                  sigma4_02, sigma4_12,     zero4,-sigma4_23, &
                  sigma4_03, sigma4_13, sigma4_23,     zero4 /), (/4,4,4,4/))


  !****************************************************************************
  !****g* minkowski/metricTensor
  ! PURPOSE
  ! The metric tensor
  ! SOURCE
  !
  real, dimension(0:3,0:3), public, parameter :: metricTensor = reshape((/ 1., 0., 0., 0., &
                                                                           0.,-1., 0., 0., &
                                                                           0., 0.,-1., 0., &
                                                                           0., 0., 0.,-1. /),(/4,4/))
  !****************************************************************************

  public :: SP,abs4,abs4Sq,abs3,op_ang
  public :: slashed,slashed5,contract,sigma4,tilde,levi_civita



  !****************************************************************************
  !****s* minkowski/contract
  ! NAME
  ! real function Contract(a,b)
  ! PURPOSE
  ! Evaluates a^(mu nu) b_(mu nu)
  ! INPUTS
  ! * a,b : matrices a^(mu nu) b^(mu nu)
  ! * The matrices can be real or complex (therefore we use an interface)
  ! OUTPUT
  ! * real
  !****************************************************************************
  Interface Contract
     Module Procedure ContractCC  !, ContractRR, ContractCR, ContractRC
  End Interface


  !****************************************************************************
  !****f* minkowski/abs4
  ! NAME
  ! real function abs4(a)
  ! real function abs4(a, flagOk)
  !
  ! PURPOSE
  ! Absolute value of a 4-vector.
  ! INPUTS
  ! * real,dimension(0:3)  :: a
  ! OUTPUT
  ! * real  ::    abs4=sqrt(a(0)*a(0)-a(1)*a(1)-a(2)*a(2)-a(3)*a(3))
  !****************************************************************************
  Interface abs4
     Module Procedure abs4a, abs4b
  end Interface abs4

contains


  !****************************************************************************
  !****f* minkowski/op_ang
  ! NAME
  ! real function op_ang(p1,p2)
  ! PURPOSE
  ! Computes the opening angle between the spatial components of two 4-vectors.
  ! INPUTS
  ! * real,dimension(0:3)  :: p1,p2
  ! OUTPUT
  ! * opening angle in degrees [0...180]
  !****************************************************************************
  real function op_ang(p1,p2)
    use constants, only: pi
    real, dimension(0:3), intent(in) :: p1,p2
    op_ang = acos( dot_product(p1(1:3),p2(1:3)) / (abs3(p1)*abs3(p2)) ) * 180./pi
  end function op_ang


  !****************************************************************************
  !****f* minkowski/abs3
  ! NAME
  ! function abs3(a)
  ! PURPOSE
  ! Absolute value of the spatial components of a 4-vector.
  ! INPUTS
  ! * real,dimension(0:3)  :: a
  ! OUTPUT
  ! * real  ::    abs3=sqrt(a(1)**2+a(2)**2+a(3)**2)
  !****************************************************************************
  real function abs3(a)
    real, dimension(0:3), intent(in) :: a
    abs3 = sqrt(Dot_Product(a(1:3),a(1:3)))
  end function abs3


  !****************************************************************************
  ! cf. interface "abs4" :
  !****************************************************************************
  real function abs4a(a)
    real, dimension(0:3), intent(in) :: a
    abs4a = sqrt(SP(a,a))
  end function abs4a
  !-------------------------------------------------------------------------
  real function abs4b(a,flagOk)
    real, dimension(0:3), intent(in) :: a
    logical, intent(out) :: flagOk
    abs4b = SP(a,a)
    if (abs4b.ge.0) then
       flagOk = .true.
       abs4b=sqrt(abs4b)
    else
       flagOk = .false.
       abs4b = 0
    end if
  end function abs4b

  !****************************************************************************
  !****f* minkowski/abs4Sq
  ! NAME
  ! real function abs4Sq(a)
  ! PURPOSE
  ! Absolute value squared of a 4-Vector.
  ! INPUTS
  ! * real,dimension(0:3)  :: a ! four vector
  ! OUTPUT
  ! * real  ::    abs4Sq=a(0)*a(0)-a(1)*a(1)-a(2)*a(2)-a(3)*a(3)
  !****************************************************************************
  real function abs4Sq(a)
    real, dimension(0:3), intent(in) :: a
    abs4Sq=SP(a,a)
  end function abs4Sq


  !****************************************************************************
  !****f* minkowski/SP
  ! NAME
  ! function SP(a,b)
  ! PURPOSE
  ! Scalar Product for 4-Vectors, "mostly -" metric
  ! INPUTS
  ! * real,dimension(0:3)  :: a,b ! four vectors
  ! OUTPUT
  ! * real  ::    SP=a(0)*b(0)-a(1)*b(1)-a(2)*b(2)-a(3)*b(3)
  !****************************************************************************
  function SP(a,b)
    real :: SP
    real, dimension(0:3), intent(in) :: a,b
    SP=a(0)*b(0)-Dot_Product(a(1:3),b(1:3))
  end function SP


  !****************************************************************************
  !****f* minkowski/ContractCC
  ! NAME
  ! real function ContractCC(a,b)
  ! PURPOSE
  ! Evaluates a^(mu nu) b_(mu nu)
  ! INPUTS
  ! * complex,dimension(0:3,0:3)  :: a,b ! matrices a^(mu nu) b^(mu nu)
  ! OUTPUT
  ! * real
  !****************************************************************************
  real function ContractCC(a,b)
    complex, dimension(0:3,0:3), intent(in) :: a,b
    integer :: mu,nu
    contractCC=0.
    do mu=0,3
       do nu=0,3
          contractCC=contractCC+real(a(mu,nu)*b(mu,nu)*metricTensor(mu,mu)*metricTensor(nu,nu))
       end do
    end do
  end function CONTRACTCC


  !****************************************************************************
  !****f* minkowski/ContractRR
  ! NAME
  ! real function ContractRR(a,b)
  ! PURPOSE
  ! Evaluates a^(mu nu) b_(mu nu)
  ! INPUTS
  ! * real,dimension(0:3,0:3)  :: a,b ! matrices a^(mu nu) b^(mu nu)
  ! OUTPUT
  ! * real
  !****************************************************************************
!   real function ContractRR(a,b)
!     real, dimension(0:3,0:3), intent(in) :: a,b
!     integer :: mu,nu
!     contractRR=0.
!     do mu=0,3
!        do nu=0,3
!           contractRR=contractRR+a(mu,nu)*b(mu,nu)*metricTensor(mu,mu)*metricTensor(nu,nu)
!        end do
!     end do
!   end function CONTRACTRR

!!$  !*************************************************************************
!!$  !****f* minkowski/ContractCR
!!$  ! NAME
!!$  ! real function ContractCR(a,b)
!!$  ! PURPOSE
!!$  ! Evaluates a^(mu nu) b_(mu nu)
!!$  ! INPUTS
!!$  ! complex, dimension(0:3,0:3), intent(in) :: a ! matrix a(mu nu)
!!$  ! real, dimension(0:3,0:3), intent(in)    :: b ! matrix b^(mu nu)
!!$  ! OUTPUT
!!$  ! * real
!!$  !*************************************************************************
!   real function ContractCR(a,b)
!     implicit none
!     complex, dimension(0:3,0:3), intent(in) :: a
!     real, dimension(0:3,0:3), intent(in) :: b
!     integer :: mu,nu,alpha,beta
!     contractCR=0.
!     do mu=0,3
!        do nu=0,3
!           do alpha=0,3
!              if(alpha.ne.mu) cycle
!              do beta=0,3
!                 if(beta.ne.nu) cycle
!                 contractCR=contractCR+real(a(mu,nu)*b(alpha,beta)*metricTensor(mu,alpha)*metricTensor(nu,beta))
!              end do
!           end do
!        end do
!     end do
!   end function CONTRACTCR


!!$  !*************************************************************************
!!$  !****f* minkowski/ContractRC
!!$  ! NAME
!!$  ! real function ContractRC(a,b)
!!$  ! PURPOSE
!!$  ! Evaluates a^(mu nu) b_(mu nu)
!!$  ! INPUTS
!!$  ! real, dimension(0:3,0:3), intent(in) :: a ! matrix a(mu nu)
!!$  ! complex, dimension(0:3,0:3), intent(in)    :: b ! matrix b^(mu nu)
!!$  ! OUTPUT
!!$  ! * real
!!$  !*************************************************************************
!   real function ContractRC(a,b)
!     implicit none
!     real, dimension(0:3,0:3), intent(in) :: a
!     complex, dimension(0:3,0:3), intent(in) :: b
!     !integer :: mu,nu,alpha,beta
!     contractRC=ContractCR(b,a)
!   end function CONTRACTRC


  !****************************************************************************
  !****f* minkowski/sigma4
  ! NAME
  ! complex function sigma4(a)
  ! PURPOSE
  ! Returns sigma^(mu nu)=i/2 [gamma^mu, gamma^nu]
  ! INPUTS
  ! OUTPUT
  ! * complex, dimension(0:3,0:3) :: matrix
  !****************************************************************************
  function sigma4(mu,nu) result(matrix)
    complex, dimension(0:3,0:3) :: matrix
    integer,intent(in) :: mu, nu

!    matrix=ii/2.*(MatMul(gamma(:,:,mu),gamma(:,:,nu))-MatMul(gamma(:,:,nu),gamma(:,:,mu)))
    matrix = sigma4par(:,:,mu,nu)

  end function sigma4


  !****************************************************************************
  !****f* minkowski/slashed
  ! NAME
  ! function slashed(p) result(matrix)
  ! PURPOSE
  ! Evaluates gamma^mu*p_mu
  ! INPUTS
  ! * real, dimension(0:3) :: p
  ! OUTPUT
  ! * complex, dimension(0:3,0:3) :: matrix
  !****************************************************************************
  function slashed(p) result(matrix)
    real, intent(in),  dimension(0:3) :: p
    complex, dimension(0:3,0:3) :: matrix
    integer :: mu
    matrix=0.
    do mu=0,3
       matrix=matrix+p(mu)*metricTensor(mu,mu)*gamma(:,:,mu)
    end do
  end function slashed


  !****************************************************************************
  !****f* minkowski/slashed5
  ! NAME
  ! function slashed5(p) result(matrix)
  ! PURPOSE
  ! Evaluates gamma^mu*p_mu*gamma_5
  ! INPUTS
  ! * real, dimension(0:3) :: p
  ! OUTPUT
  ! * complex, dimension(0:3,0:3) :: matrix
  !****************************************************************************
  function slashed5(p) result(matrix)
    real, intent(in),  dimension(0:3) :: p
    complex, dimension(0:3,0:3) :: matrix
    integer :: mu
    matrix=0.
    do mu=0,3
       matrix=matrix+p(mu)*metricTensor(mu,mu)*gamma(:,:,mu+8)
    end do
  end function slashed5



  !****************************************************************************
  !****f* minkowski/tilde
  ! NAME
  ! function tilde(a) result(a_tilde)
  ! PURPOSE
  ! Evaluates a^tilde=gamma_0 a^dagger gamma_0
  ! INPUTS
  ! * complex, dimension(0:3,0:3) :: a
  ! OUTPUT
  ! * complex, dimension(0:3,0:3) :: a_tilde
  !****************************************************************************
  function tilde(a) result(a_tilde)
    use matrix_module, only: dagger, MatrixMult
    complex, dimension(0:3,0:3),intent(in) :: a
    complex, dimension(0:3,0:3) :: a_tilde
    a_tilde = MatrixMult(gamma0,dagger(a),gamma0)
  end function tilde


  !****************************************************************************
  !****f* minkowski/levi_civita
  ! NAME
  ! integer function levi_civita (i, j, k, l)
  ! PURPOSE
  ! Calculates the fully antisymmetric Levi-Civita tensor \epsilon_{ijkl}
  ! in four dimensions.
  ! INPUTS
  ! * integer :: i, j, k, l   --- tensor indices
  ! OUTPUT
  ! * tensor value for the given indices
  !****************************************************************************
  integer function levi_civita (i, j, k, l)
    integer, intent(in) :: i, j, k, l
    levi_civita = (i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l) / 12
  end function


end module minkowski
