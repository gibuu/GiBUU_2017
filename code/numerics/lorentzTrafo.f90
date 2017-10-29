!******************************************************************************
!****m* /lorentzTrafo
! NAME
! module lorentzTrafo
! PURPOSE
! Implements lorentz transformation.
!******************************************************************************
module lorentzTrafo
  implicit none
  private

  public :: lorentz
  public :: lorentzCalcBeta
  public :: BoostTensor
  public :: eval_sigmaBoost
!!$  Public :: BoostToCM


  !****************************************************************************
  !****f* lorentzTrafo/lorentzCalcBeta
  ! NAME
  ! function lorentzCalcBeta (Mom, CallName) result (beta)
  ! function lorentzCalcBeta (mom3, mass, CallName) result (beta)
  ! PURPOSE
  ! calculate the beta-Vector for a Lorentz-Boost.
  ! checks whether a valid vector results
  ! INPUTS
  ! * real,dimension(0:3)    :: Mom
  ! * character(40),optional :: CallName
  ! or:
  ! * real, dimension(1:3)   :: mom3
  ! * real                   :: mass
  ! * character(40),optional :: CallName
  ! RESULT
  ! * real, dimension(1:3) :: beta
  !****************************************************************************
  interface lorentzCalcBeta
     module procedure lorentzCalcBeta4,lorentzCalcBeta3
  end interface

contains


!!$  !************************************************************************
!!$  !****s* lorentzTrafo/BoostToCM
!!$  ! NAME
!!$  ! subroutine BoostToCM(a,b)
!!$  ! PURPOSE
!!$  ! Boosts vectors a and b to their common center of momentum frame.
!!$  ! INPUTS
!!$  ! * real, dimension(0:3), intent(inout) :: a,b  -- vectors a,b in some frame
!!$  ! RESULT
!!$  ! * real, dimension(0:3), intent(inout) :: a,b  -- vectors a,b in CM frame
!!$  !************************************************************************
!!$  subroutine BoostToCM(a,b)
!!$    real, dimension(0:3), intent(inout) :: a,b
!!$    real, dimension(1:3) :: beta
!!$    real, dimension(0:3) :: total
!!$    total=a+b
!!$    call lorentzCalcBeta(beta,total, 'BoostToCM')
!!$    CALL lorentz(beta,a,'BoostToCM')
!!$    CALL lorentz(beta,b,'BoostToCM')
!!$  end subroutine BoostToCM



  !****************************************************************************
  !****s* lorentzTrafo/lorentz
  ! NAME
  ! subroutine lorentz(beta,fourVector,CallName)
  ! PURPOSE
  ! performs Lorentz transformation of fourVector into a
  ! system which is traveling with the velocity beta(1:3)
  ! INPUTS
  ! * real,dimension(0:3),intent(inout) :: fourVector
  ! * real,dimension(1:3),intent(in)    :: beta
  ! * character(*),intent(in),optional  :: CallName
  ! RESULT
  ! * real,dimension(0:3),intent(inout) :: fourVector
  !****************************************************************************
  subroutine lorentz (beta, fourVector, CallName)
    use callstack, only: traceback

    real,dimension(0:3),intent(inout) :: fourVector
    real,dimension(1:3),intent(in)    :: beta
    character(*),optional,intent(in)  :: CallName

    real :: gamma, betaFour

    !Evaluate gamma

    gamma = Dot_Product(beta,beta)
    if (gamma<1e-20) return ! do nothing
    if (gamma>=1.0) then
       write(*,*)
       if (present(CallName)) then
          write(*,*) 'lorentz called by ',CallName
       end if
       write(*,*) '(1-beta**2) in lorentz less or equal zero: ', 1.0-gamma
       write(*,*) 'beta=',  beta
       write(*,*) 'Stop program'
       call TRACEBACK()
    end if
    gamma = 1.0/sqrt(1.0-gamma)

    !beta*fourVector(1:3)
    betaFour = Dot_Product(beta, fourVector(1:3))
    !do transformation
    fourVector(1:3)= fourVector(1:3) + gamma*beta(1:3)*(gamma/(gamma+1.)*betaFour-fourvector(0))
    fourVector(0)= gamma*(fourVector(0)-betaFour)
    return
  end subroutine lorentz


  !****************************************************************************
  ! cf. interface lorentzCalcBeta
  !****************************************************************************
  function lorentzCalcBeta4 (Mom, CallName) result (beta)
    use minkowski, only: abs4Sq
    use callstack, only: traceback

    real, dimension(0:3), intent(in)   :: Mom
    character(*), intent(in), optional :: CallName
    real, dimension(1:3) :: beta

    character(80) :: Name
    real, parameter :: eps = 1e-10 ! accuracy

    if (present(CallName)) then
       Name = 'lorentzCalcBeta  (called by '//CallName//')'
    else
       Name = 'lorentzCalcBeta'
    end if

    if (abs4Sq(Mom)<eps*Mom(0)**2) then
       write(*,*) Name,': not a valid boost vector!'
       write(*,*) Mom,Mom(0)**2-Dot_Product(Mom(1:3),Mom(1:3))
       call traceBack()
       stop
    end if

    beta(1:3) = Mom(1:3)/Mom(0)

  end function lorentzCalcBeta4
  !---------------------------------------------------------------------------
  function lorentzCalcBeta3 (mom3, mass, CallName) result (beta)
    real, dimension(1:3), intent(in)   :: mom3
    real, intent(in)                   :: mass
    character(*), intent(in), optional :: CallName
    real, dimension(1:3) :: beta

    real,dimension(0:3)  :: Mom

    Mom(0) = sqrt(mass**2+sum(mom3**2))
    Mom(1:3) = mom3(1:3)
    beta = lorentzCalcBeta4 (Mom, CallName)
  end function lorentzCalcBeta3



  !****************************************************************************
  !****s* lorentzTrafo/BoostTensor
  ! NAME
  ! subroutine BoostTensor(u,e,e0)
  ! PURPOSE
  ! Lorentz-Transformation of a tensor e(0:3,0:3) into a system moving
  ! with 4-velocity u(0,3).
  ! INPUTS
  ! * real,dimension(0:3),intent(in)     :: u -- boost 4-velocity
  ! * real,dimension(0:3,0:3),intent(in) :: e -- tensor in CF
  ! RESULT
  ! * real,dimension(0:3,0:3),intent(in) :: e0 -- tensor in LRF
  ! NOTES
  ! Lorentztransformation des Energie-Impuls-Tensor in das Ruhesystem
  !****************************************************************************
  subroutine BoostTensor(u,e,e0)

    real, dimension(0:3),     intent(in) :: u
    real, dimension(0:3,0:3), intent(in) :: e
    real, dimension(0:3,0:3), intent(out) :: e0

    real, dimension(0:3,0:3) :: d,lambda
    real, dimension(1:3)     :: beta

    integer :: i,j,k,l
    real :: qbeta,gamma

    do i=1,3
       do j=1,3
          d(i,j) = 0.0
       end do
    end do

    do i=1,3
       d(i,i) = 1.0
    end do

    gamma = u(0)

    qbeta = 0.0
    do i=1,3
       beta(i) = -u(i)/u(0)
       qbeta = qbeta + beta(i)**2
    end do

    lambda(0,0) = gamma
    do i=1,3
       lambda(0,i) = gamma*beta(i)
       lambda(i,0) = lambda(0,i)
    end do

    do i=1,3
       do j=1,3
          lambda(i,j) = d(i,j)+(gamma-1.)*beta(i)*beta(j)/(qbeta+0.000001)
       end do
    end do

    do k=0,3
       do l=0,3

          e0(k,l) = 0.0

          do i=0,3
             do j=0,3
                e0(k,l) = e0(k,l) + lambda(i,k)*lambda(j,l)*e(i,j)
             end do
          end do

       end do
    end do

  end subroutine BoostTensor



  !****************************************************************************
  !****f* lorentzTrafo/eval_sigmaBoost
  ! NAME
  ! real function eval_sigmaBoost(mom1,mom2)
  !
  ! PURPOSE
  ! This function calculates the boost factor for a cross section. Given a
  ! cross section "sigma" which is defined the rest frame of particles A or B,
  ! the cross section "sigmaPrime" in a frame where both a moving
  ! is given by:
  ! * sigmaPrime=sigma* eval_sigmaBoost
  !
  ! INPUTS
  ! * real , dimension(0:3),intent(in) :: mom1,mom2 -- 4-momenta of the colliding pair
  !
  ! OUTPUT
  ! * Boost factor
  !
  ! NOTES
  ! * For derivation confer Oliver's Phd thesis (Appendix G)
  ! * Note that the formula given in Effenberger's diploma thesis is only
  !   very approximate!!
  !****************************************************************************
  real function eval_sigmaBoost (mom1, mom2)
    use particleDefinition
    use minkowski, only: abs3, abs4
    use callStack, only: TRACEBACK

    real, dimension(0:3), intent(in) :: mom1, mom2

    real, dimension(0:3,1:2) :: mom
    real, dimension(1:2) :: value
    real, dimension(1:3) :: BETA
    integer :: i,j,numTries
    real :: factor
    logical,dimension(1:2) :: isMassLess
    real, parameter :: epsilon=1.E-8

    mom(:,1)=mom1
    mom(:,2)=mom2

    ! Check whether one of the particles is massless
    isMassless(1)= (abs4(mom1)<epsilon)
    isMassless(2)= (abs4(mom2)<epsilon)

    if (isMassless(1).and.isMassless(2)) then
       write(*,*) 'error in eval_sigmaBoost! Both momenta are massless!!!'
       write(*,*) mom1
       write(*,*) mom2
       write(*,*) 'STOP !!'
       call TRACEBACK()
       stop
    end if

    !**************************************************************************
    ! Boost to restframes of particles (first to the one of particle 1, then
    ! to the one of  2) and evaluate the boost factor in both frames.
    !
    ! In principle this procedure is an overkill: we could also simply choose
    ! one of the two rest frames since the result should not dependend on this
    ! choice! However, we do it right now as a numerical check!
    !**************************************************************************
    factor=1./abs3(mom(0,1)*mom(:,2)-mom(0,2)*mom(:,1))
    numTries=0
    value=0.
    do i=1,2
       if (isMassless(i)) cycle ! Cycle since we can't boost to a rest frame of a massless particle
       j=3-i
       mom(:,1)=mom1
       mom(:,2)=mom2
       beta=mom(1:3,i)/mom(0,i)
       call lorentz(beta,mom(:,j))
       value(i)=factor*abs4(mom(:,i))*abs3(mom(:,j))
       numTries=numTries+1
    end do

    eval_sigmaBoost=sum(value)/float(numTries)

    ! Numtries is greater 1 if both particles have mass. Then we make a
    ! consistency check:
    if (numTries>1) then
       ! check on numerics:
       if (abs(eval_sigmaBoost-value(1))>0.001) then
          ! Both entries in value should be identical!!!
          write(*,*) 'Critical Error in eval_sigmaBoost', value
          STOP
       end if
    end if
  end function eval_sigmaBoost


end module lorentzTrafo
