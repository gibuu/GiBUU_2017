!***************************************************************************
!****m* /Output
! NAME
! module Output
!
! PURPOSE
! This module contains a subroutine which prints out the fragment vector.
!***************************************************************************
module Output

  implicit none
  PRIVATE

  PUBLIC :: WriteGlobalVector

  contains

    subroutine WriteGlobalVector(iev,FragmentVector)
      use InputGeneral,     only : SubEvents
      use InputCoalescence, only : Out4
      use typeDefinitions,  only : cluster

      type(cluster), dimension(:), intent(in) :: FragmentVector
      integer, intent(in)                     :: iev

      integer              :: i,idp,Z,A,Y
      real                 :: Mass
      real, dimension(1:3) :: xx,pp
      character(2)         :: LS

      if (iev==1) open(Unit=4,file=Out4)

      do i=1,size(FragmentVector,dim=1)
         if (FragmentVector(i)%ID == 0) cycle
         idp    = FragmentVector(i)%ID
         Z      = FragmentVector(i)%ChargeNumber
         A      = FragmentVector(i)%MassNumber
         Y      = FragmentVector(i)%HypNumber
         LS     = FragmentVector(i)%HypType
         Mass   = FragmentVector(i)%Mass*0.19733
         xx(:)  = FragmentVector(i)%position(:)
         pp(1:3)= FragmentVector(i)%momentum(1:3)
         write(4,10) idp,Z,A,Y,LS,Mass,xx,pp,iev
      end do

      if (iev==SubEvents) close(Unit=4)

10    format(i4,3(1x,i3),1x,a2,1x,f5.3,3(1x,f8.3),3(1x,f8.3),1x,i4)

    end subroutine WriteGlobalVector

end module Output
