!******************************************************************************
!****m* /leptonTensor
! NAME
! module leptonTensor
! PURPOSE
! * Provides the lepton-tensor for the lepton vertex in electro-weak processes
! * NOTE: lepton tensor contains coupling constants and propagator**2
!******************************************************************************

module leptonTensor
  private

  public :: l_munu,l_munu_Matrix, electronMass,leptonicTensor

  !Electron mass in GeV
  real, parameter :: electronMass=0.00051099892

contains
  !****************************************************************************
  !****f* leptonTensor/l_munu
  ! NAME
  ! real function l_munu(mu,nu,lin,lout)
  !
  ! PURPOSE
  ! * Provides the electromagnetic lepton tensor
  ! * Evaluates  L_mu,nu * m_e^2 =1/2 e^2/q^4 1/4 Tr[ {GS(lin)+me} gamma_mu {GS(l_f) + me}  gamma_nu ] * m_e^2
  !   according to Bjorken-Drell Notation
  ! * NOTE: lepton tensor contains coupling constants and propagator**2
  !
  ! INPUTS
  ! * real, dimension(0:3),intent(in) :: lin,lout  ! incoming and outgoing lepton 4-momentum
  ! * integer ::  mu , nu ! Indices of the tensor
  !
  ! OUTPUT
  ! *   real function l_munu(mu,nu,lin,lout)
  !****************************************************************************
  real function l_munu(mu,nu,lin,lout)

    use minkowski, only: metricTensor,SP
    use constants, only: pi
    implicit none
    real, dimension(0:3),intent(in) :: lin,lout
    real, parameter :: electronChargeSQ = 1./137.04*4.*pi ! =e^2
    integer :: nu, mu
    real, dimension(0:3) :: q

    q=lin-lout ! momentum of virtual photon

    l_munu=electronChargeSQ/(SP(q,q))**2* &
         & (electronMass**2*metricTensor(mu,nu)+lin(mu)*lout(nu)+lin(nu)*lout(mu)-SP(lin,lout)*metricTensor(mu,nu))/2.
  end function l_munu


  !****************************************************************************
  !****f* leptonTensor/l_munu_matrix
  ! NAME
  ! function l_munu_matrix(lin,lout) result(L)
  !
  ! PURPOSE
  ! * Provides the electromagnetic lepton tensor
  ! * Evaluates  L_mu,nu * (2m_e)^2 e^2 =1/2 e^4/q^4 Tr[ {GS(lin)+me} gamma_mu {GS(l_f) + me}  gamma_nu ] * m_e^2
  !   according to Bjorken-Drell Notation
  !
  ! INPUTS
  ! * real, dimension(0:3),intent(in) :: lin,lout  ! incoming and outgoing lepton 4-momentum
  !
  ! OUTPUT
  ! * real, dimension(0:3,0:3) :: L
  !****************************************************************************
  function l_munu_matrix(lin,lout) result(L)

    use constants, only: pi
    implicit none
    real, dimension(0:3,0:3)        :: L
    real, dimension(0:3),intent(in) :: lin,lout
    integer ::  mu , nu ! Indices of the tensor
    real, parameter :: electronChargeSQ = 1./137.04*4.*pi

    do mu=0,3
       do nu=0,3
          L(mu,nu)=l_munu(mu,nu,lin,lout)*4.*electronChargeSQ
       end do
    end do
  end function l_munu_matrix





  !****************************************************************************
  !****f* leptonTensor/leptonicTensor
  ! NAME
  ! function leptonicTensor(process_ID,flavor_ID,k_in,k_out) result(matrix)
  !
  ! PURPOSE
  ! * Provides the full electro-weak lepton tensor for lepton spinor normalization   baru.u = 2m
  ! * NOTE: lepton tensor contains coupling constants and propagator**2
  !
  ! INPUTS
  ! * integer, intent(in)             :: process_ID  -- EM, CC or NC (or anti if negative)
  ! * real, dimension(0:3),intent(in) :: k_in,k_out  -- incoming and outgoing lepton 4-momentum
  !
  ! OUTPUT
  ! * real, dimension(0:3,0:3) :: matrix
  !****************************************************************************
  function leptonicTensor(process_ID,k_in,k_out) result(matrix)
    use minkowski, only: gamma,gamma5,SP,slashed
    use matrix_module, only: unit4,trace,matrixmult
    use leptonicID
    use constants, only: pi,alphaQED,GF,coscab
    implicit none

    integer, intent(in)             :: process_ID
    real, dimension(0:3),intent(in) :: k_in,k_out

    complex, dimension(0:3,0:3)    :: matrix

    complex, dimension(0:3,0:3) :: projector_in,projector_out,j_mu,j_nu
    real, dimension(0:3)        :: q
    real                        :: mi,mf
    real                        :: Qs
    integer                     :: mu,nu
    real                        :: coupling
    real                        :: a


    if (process_ID.gt.0) then
       a=1.
    else !anti-particle
       a=-1.
    end if
    if (process_ID.eq.EM) a=0.

    q=k_in-k_out
    Qs=-SP(q,q)
    mi=sqrt(max(SP(k_in,k_in),0.))
    mf=sqrt(max(SP(k_out,k_out),0.))

    projector_out=slashed(k_out)+mf*unit4
    projector_in=slashed(k_in)+mi*unit4

    matrix=0.

    do mu=0,3
       do nu=0,3

          select case (process_ID)

          case (NC,antiNC)
             coupling=GF**2/2.

          case (CC,antiCC)
             coupling=GF**2/2.*coscab**2

          case (EM,antiEM)
             coupling=(4.*pi*alphaQED)**2/Qs**2*1./2.    !1/2 due to lepton spin averaging

            case default
             write(*,*) 'problem with processID in leptonTensor -> STOP'
             stop
          end select

          j_mu=MatMul(gamma(:,:,mu),unit4-a*gamma5)
          j_nu=MatMul(gamma(:,:,nu),unit4-a*gamma5)

          matrix(mu,nu)=coupling*Trace(MatrixMult(projector_out,j_mu,projector_in,j_nu))

       end do
    end do

  end function leptonicTensor

end module leptonTensor
