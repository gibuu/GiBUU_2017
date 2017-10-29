program main
  use GridOrdering

  implicit none

  integer :: i


  call GridOrdering_Init



  do i=lbound(nDistance,dim=1),ubound(nDistance,dim=1)
     if (nDistance(i,2)<0) cycle
     call GridOrdering_RandomizeRadius(i)

  end do



end program main
