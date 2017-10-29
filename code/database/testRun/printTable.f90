program printTable
  use inputGeneral, only: readInputGeneral
  use particleProperties, only: initParticleProperties, hadron, nDecays

  implicit none

  integer :: ID, j

  call readInputGeneral
  call initParticleProperties

  do ID=1,121

     do j= 1, nDecays
       write(*,'(I3,6I4,6F12.5)') ID,hadron(ID)%decaysID(j),hadron(ID)%decays(j)
       write(120,'(I3,6I4,6F12.5)') ID,hadron(ID)%decaysID(j),hadron(ID)%decays(j)
     end do

  end do


end program printTable
