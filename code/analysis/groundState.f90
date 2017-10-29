!******************************************************************************
!****m* /groundStateAnalysis
! NAME
! module groundStateAnalysis
! PURPOSE
! ...
!******************************************************************************
module groundStateAnalysis


contains


  subroutine countMass(teilchen)
    use particleDefinition
    use histf90
    use idTable, only: nucleon
    implicit none
    type(particle), dimension(:,:) :: teilchen
    integer :: i,j
    real :: pAbs
    type(histogram) :: histMass, histMomentum, histMomentum_p, histMomentum_n

    call CreateHist(histMass, 'Massen',0.,1.,0.001)
    call CreateHist(histMomentum, 'Momentum',0.,1.,0.001)
    call CreateHist(histMomentum_p, 'Momentum_p',0.,1.,0.001)
    call CreateHist(histMomentum_n, 'Momentum_n',0.,1.,0.001)

    do i=lbound(teilchen,dim=1),ubound(teilchen,dim=1)
       do j=lbound(teilchen,dim=2),ubound(teilchen,dim=2)
          if (teilchen(i,j)%ID.ne.nucleon) cycle
          !          write(*,*)teilchen(i,j)%mass
          call AddHist(histMass,teilchen(i,j)%mass,1.)

          pAbs=sqrt(Dot_Product(teilchen(i,j)%momentum(1:3),teilchen(i,j)%momentum(1:3) ))
          if (teilchen(i,j)%charge.eq.0) then
             call AddHist(histMomentum_n,pAbs,1.)

          else
             call AddHist(histMomentum_p,pAbs,1.)
          end if
          call AddHist(histMomentum,pAbs,1.)
       end do
    end do
    open(11,file='masses.dat')
    open(12,file='momentum.dat')
    open(13,file='momentum_n.dat')
    open(14,file='momentum_p.dat')
    call WriteHist(histMass,11)
    call WriteHist(histMomentum,12)
    call WriteHist(histMomentum_n,13)
    call WriteHist(histMomentum_p,14)
    close(11)
    close(12)
    close(13)
    close(14)


  end subroutine countMass


end module groundStateAnalysis
