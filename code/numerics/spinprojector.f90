!******************************************************************************
!****m* /spinProjector
! NAME
! module spinProjector
! PURPOSE
! * Provides the spin projector for spin 3/2 resonances
!******************************************************************************

module spinProjector

  implicit none
  private

  public :: spin32proj, spin32proj_tensor

contains

  !****************************************************************************
  !****f* spinProjector/spin32proj
  ! NAME
  ! function spin32proj(resID,mu,nu,p_res)
  !
  ! PURPOSE
  ! * provides the spin projector for spin 3/2 resonances
  !   as defined e.g. in PRC 73, 065502 (2006), eq. (21)
  ! * return value is contravariant: spin32proj^(mu nu)
  !
  ! INPUTS
  ! * real, dimension(0:3),intent(in) :: p_res  ! resonance 4-momentum
  ! * integer ::  mu , nu                       ! Lorentz indices
  ! * integer :: resID                          ! ID of the resonance
  ! * logical, optional :: useInvMass           ! if true, use invariant mass in projector, if false, use pole mass (needed for u-channel)
  !
  ! OUTPUT
  ! *  spin32proj is complex with dimension(0:3,0:3)
  !****************************************************************************
  function spin32proj(resID,mu,nu,p_res,useInvMass) result(getspin32proj)
    use minkowski, only: abs4, gamma, slashed, metricTensor
    use matrix_module, only: unit4
    use particleproperties, only: hadron

    complex, dimension(0:3,0:3) :: getspin32proj
    real, dimension(0:3),intent(in) :: p_res
    integer, intent(in) :: nu, mu, resID
    logical, intent(in), optional :: useInvMass

    real :: W
    complex, dimension(0:3,0:3) :: spinproj1, spinproj2

    logical :: useInvariantMass    !if true, use invariant mass in projector, if false, use pole mass

    if ( present(useInvMass) ) then
       useInvariantMass=useInvMass
    else
       useInvariantMass=.true.
    end if

    if (mu.lt.0.or.mu.gt.3.or.nu.lt.0.or.nu.gt.3) then
       write(*,*) 'something strange with mu, nu in spin32proj -> STOP', mu, nu
       stop
    end if

    if (useInvariantMass) then
       W=abs4(p_res)
    else
       W=hadron(resID)%mass
    end if

    spinproj1= -(slashed(p_res) + W*unit4)

    spinproj2= metricTensor(mu,nu)*unit4 - 2./3.*p_res(mu)*p_res(nu)/W**2*unit4   &
         &     + 1./3.*(p_res(mu)*gamma(:,:,nu) - p_res(nu)*gamma(:,:,mu))/W - 1./3.* matmul(gamma(:,:,mu),gamma(:,:,nu))


    getspin32proj=matmul(spinproj1,spinproj2)

  end function spin32proj


  function spin32proj_tensor(resID,p_res) result(getspin32proj)
    use minkowski, only: abs4, gamma, slashed, metricTensor
    use matrix_module, only: unit4
    use particleproperties, only: hadron

    complex, dimension(0:3,0:3,0:3,0:3) :: getspin32proj
    real, dimension(0:3),intent(in) :: p_res
    integer, intent(in) :: resID
    integer  :: nu, mu

    real :: W
    complex, dimension(0:3,0:3) :: spinproj1, spinproj2

    logical :: useInvariantMass=.true.    !if true, use invariant mass in projector, if false, use pole mass

    if (useInvariantMass) then
       W=abs4(p_res)
    else
       W=hadron(resID)%mass
    end if
    spinproj1= -(slashed(p_res) + W*unit4)
    do mu=0,3
       do nu=0,3

          spinproj2= metricTensor(mu,nu)*unit4 - 2./3.*p_res(mu)*p_res(nu)/W**2*unit4   &
               &     + 1./3.*(p_res(mu)*gamma(:,:,nu) - p_res(nu)*gamma(:,:,mu))/W - 1./3.* matmul(gamma(:,:,mu),gamma(:,:,nu))

          getspin32proj(mu,nu,:,:)=matmul(spinproj1,spinproj2)

       end do
    end do
  end function spin32proj_tensor


end module spinProjector
