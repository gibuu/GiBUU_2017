module initInABoxDelta

  implicit none
  private

  real, parameter :: mom=0.1
  real, parameter ::  mass=1.2
  integer, parameter ::  charge=-1

  public :: InitInABoxDelta_init

contains

  subroutine InitInABoxDelta_init(part,flag)
    ! Initialisiert Deltas, die in der x-y Ebene mit z=0 starten und dann mit gegebenem p_z in die z-Richtung loslaufen.
    use particleDefinition
    use random, only: rn
    use idTable, only: delta
!     use deltaWidth, only: deloset

    integer :: i,j
    type(particle), dimension(:,:) :: part
    logical, intent(in) :: flag
!     real :: imsig2,imsig3,imsigq,rho

!!$    open(100,file="width_oset.dat")
!!$    do i=0,100
!!$       rho=0.002*float(i)
!!$       call  deloset(mass,rho,imsig2,imsig3,imsigq)
!!$       write(100,'(5G18.4)') rho,imsig2,imsig3, imsigq,2*(imsig2+imsig3+imsigq)
!!$    end do
!!$    close(100)

    do i=lbound(part,dim=1),ubound(part,dim=1)
       call setToDefault(part(i,:))
       do j=1,50
          part(i,j)%ID=delta
          part(i,j)%charge=charge
          part(i,j)%mass=mass
          part(i,j)%momentum(1:3) = (/ 0., 0., mom /)
          part(i,j)%momentum(0)   = sqrt(mass**2 + mom**2)
          part(i,j)%position(1:3) = (/ rn()*2., rn()*2., 0. /)
          part(i,j)%perturbative = flag
       end do
    end do
  end subroutine InitInABoxDelta_init


end module initInABoxDelta
