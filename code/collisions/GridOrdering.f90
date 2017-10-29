!******************************************************************************
!****m* /GridOrdering
! NAME
! module Gridordering
!
! PURPOSE
! This module defines some routines for accesing in a 3D grid the
! closest neighbours etc.
!
! INPUTS
! ---
!******************************************************************************
module GridOrdering

  implicit none

  private

  !****************************************************************************
  !****g* GridOrdering/nX
  ! SOURCE
  !
  integer, parameter :: nX = 5
  ! PURPOSE
  ! Grid elements [-nX:nX, -nX:nX, -nX:nX] are considered
  !****************************************************************************

  !****************************************************************************
  !****g* GridOrdering/nR
  ! SOURCE
  !
  integer, parameter :: nR = (2*nX+1)**3
  ! PURPOSE
  ! The number of elements: (-nX:nX)**3
  !****************************************************************************

  !****************************************************************************
  !****g* GridOrdering/DeltaV
  ! SOURCE
  !
  integer, dimension(nR,3), public :: DeltaV
  ! PURPOSE
  ! Store the 3D Delta of the indizes
  !****************************************************************************

  !****************************************************************************
  !****g* GridOrdering/iDeltaV
  ! SOURCE
  !
  integer, dimension(nR), public   :: iDeltaV
  ! PURPOSE
  ! The index to access the information stored in "DeltaV":
  ! access via DeltaV(iDeltaV(i)) returns a sorted list
  !
  !****************************************************************************

  !****************************************************************************
  !****g* GridOrdering/nDistance
  ! SOURCE
  !
  integer, dimension(0:nX**2,2), public :: nDistance
  ! PURPOSE
  ! ...
  !****************************************************************************



  public :: GridOrdering_Init
  public :: GridOrdering_RandomizeRadius


contains
  !****************************************************************************
  !****s* GridOrdering/GridOrdering_Init
  ! NAME
  ! subroutine GridOrdering_Init
  ! PURPOSE
  ! ...
  !****************************************************************************
  subroutine GridOrdering_Init
    use sorting

    implicit none

    real, dimension(nR) :: r
    integer :: ix,iy,iz, ii, ii0,ii1

    ii = 0
    do ix=-nX,nX
       do iy=-nX,nX
          do iz=-nX,nX
             ii=ii+1
             ! Distance^2: via Distance of cell centers:
!             r(ii) = ix**2+iy**2+iz**2

             ! Distance^2: via minimal Distance of cell corners:

             r(ii) = max(abs(ix)-1,0)**2 &
                  & + max(abs(iy)-1,0)**2 &
                  & + max(abs(iz)-1,0)**2
             if (ix.ne.0.or.iy.ne.0.or.iz.ne.0) r(ii) = r(ii)+1


             DeltaV(ii,1:3) = (/ix,iy,iz/)
          end do
       end do
    end do

    call indexx(r, iDeltaV)

    nDistance(:,1) = 0
    nDistance(:,2) =-1
    ii0 = 0
    ii1 = 1
    do ii=1,nR
       if (r(iDeltaV(ii)) > ii0) then
          nDistance(ii0,1)=ii1
          nDistance(ii0,2)=ii-1
          ii0 = int(r(iDeltaV(ii)))
          if (ii0 > nX**2) exit
          ii1=ii
       end if
    end do


!!$    write(*,*) '===='
!!$    do ii=1,nR
!!$       if (r(iDeltaV(ii))<=nX**2) then
!!$          write(*,'(i5,f12.0,3i5)') ii,r(iDeltaV(ii)),&
!!$               & DeltaV(iDeltaV(ii),1),&
!!$               & DeltaV(iDeltaV(ii),2),&
!!$               & DeltaV(iDeltaV(ii),3)
!!$       end if
!!$    end do
!!$    write(*,*) '===='
!!$    do ii=0,nX**2
!!$       write(*,'(i4,2i6,i12)') ii,nDistance(ii,1),nDistance(ii,2),&
!!$            &nDistance(ii,2)-nDistance(ii,1)+1
!!$    enddo
!!$    write(*,*) '===='

  end subroutine GridOrdering_Init

  !****************************************************************************
  !****s* GridOrdering/GridOrdering_RandomizeRadius
  ! NAME
  ! subroutine GridOrdering_RandomizeRadius(iRadius)
  ! PURPOSE
  ! ...
  !****************************************************************************
  subroutine GridOrdering_RandomizeRadius(iRadius)
    use random

    implicit none

    integer, intent(IN) :: iRadius

    integer :: i,i1,i2,n, j, h

    i1 = nDistance(iRadius,1)
    i2 = nDistance(iRadius,2)

    n = i2-i1+1
    if (n<2) return

    do i=i1,i2
       j = min(i2, i1 + int(n * rn()))

!       write(*,*) 'RAN:',i1,i2,i,j


       h = iDeltaV(i)
       iDeltaV(i) = iDeltaV(j)
       iDeltaV(j) = h


    end do



!!$    write(*,*) '===='
!!$    do i=i1,i2
!!$       write(*,'(i5,f12.0,3i5)') i,float(iRadius),&
!!$            & DeltaV(iDeltaV(i),1),&
!!$            & DeltaV(iDeltaV(i),2),&
!!$            & DeltaV(iDeltaV(i),3)
!!$    end do
!!$    write(*,*) '===='



  end subroutine GridOrdering_RandomizeRadius

end module GridOrdering
