! this program simulates pi+N -> pi+pi+N and pi+N -> pi+pi+pi+N
!
! start it with by typing: time ./testnpions.x < /dev/null

program testnpions
  use particleProperties
  use nBodyPhaseSpace, only : momenta_in_3BodyPS, momenta_in_4BodyPS
  use lorentzTrafo
  use histf90
  implicit none

  real :: pLab, srts 
  real, parameter :: mp=0.138, mN=0.938
  integer :: iMC,i
  integer, parameter :: nMC=1000000
  real, dimension(3,3) :: p3
  real, dimension(4,3) :: p4
  real, dimension(1:3) :: beta
  real, dimension(0:3) :: x
  real, dimension(1:4) :: mm
  type(histogram),save :: hist2,hist3

  pLab = 2.0
  mm=(/mN,mp,mp,mp/)

  call CreateHist(hist2,'dN/dE, 2 pion',0.0,5.,0.01)
  call CreateHist(hist3,'dN/dE, 3 pion',0.0,5.,0.01)



  beta = (/0.0,0.0,-pLab/sqrt(mp**2+pLab**2)/)

  srts = sqrt(mN**2+mP**2+2*mN*sqrt(mp**2+pLab**2))

  ! check boost philosophy:

  x(1:2) = 0
  x(3) = -sqrt( (srts**2-mN**2-mp**2)**2-(2*mN*mp)**2 )/(2*srts)
  x(0) = sqrt(mN**2+x(3)**2)

  beta = x(1:3)/x(0)
  
  write(*,'(1P,4e14.5)') x
  call lorentz(beta,x)
  write(*,'(1P,4e14.5)') x ! this should now be the Nucleon at rest
  write(*,*) '========'



  do iMC=1,nMC
     call momenta_in_3BodyPS(srts,(/mN,mp,mp/),p3)
     
     do i=2,3
        x(1:3) = p3(:,i)
        x(0)   = sqrt(mm(i)**2+sum(x(1:3)**2))

!        write(*,'(1P,4e14.5)') x
        call lorentz(beta,x)
!        write(*,'(1P,4e14.5)') x
!        write(*,*)

        call AddHist(hist2,x(0),1.0)

     enddo

     call momenta_in_4BodyPS(srts,(/mN,mp,mp,mp/),p4)
     
     do i=2,4
        x(1:3) = p4(:,i)
        x(0)   = sqrt(mm(i)**2+sum(x(1:3)**2))

!        write(*,'(1P,4e14.5)') x
        call lorentz(beta,x)
!        write(*,'(1P,4e14.5)') x
!        write(*,*)

        call AddHist(hist3,x(0),1.0)

     enddo
  end do

  call WriteHist(hist2,141,mul=1./(nMC), add=1e-20)
  call WriteHist(hist3,142,mul=1./(nMC), add=1e-20)


end program testnpions
