program printMesMes

  use constants, only: mK
  use IdTable
  use inputGeneral, only: readInputGeneral, path_to_input
  use mesonMeson_Xsections
  use mesonWidth, only: fullWidthMeson
  use particleProperties, only: initParticleProperties, hadron
  use twoBodyTools, only: pCM
  use version, only: printVersion

  implicit none

  integer :: iQ,iSrts
  real :: srts,pinitial2,const,dummy
  real, dimension(21) :: msigbg
  logical :: useWidth
  real :: mKK 

  call PrintVersion

  call readInputGeneral
  call initParticleProperties

  dummy = fullWidthMeson(103, 0.9) ! Dummy call for enforce init

  call mesMes_Tabulate

  iQ = 0
!  usewidth = .false.
  usewidth = .true.
  const = 2.0
  mKK = hadron(kaonStar)%mass


  do iSrts=1,200 ! 80 ! 40
     srts = iSrts*0.02 ! 0.05 ! 0.1
     if (srts < 2*mK) cycle
     pinitial2 = pCM(srts,mK,mK)**2

!     call kkbar_cross(iQ,srts,pinitial2,const,msigbg)     
     call kkbar_cross(iQ,srts,pinitial2,const,msigbg,useWidth)

     write(21,'(99g12.5)') srts,msigbg,sum(msigbg)

     if (srts < mK+mKK) cycle
     pinitial2 = pCM(srts,mK,mKK)**2

     msigbg = 0.0

!     call kstarkbar_cross(iQ,srts,pinitial2,const,msigbg)
     call kstarkbar_cross(iQ,srts,pinitial2,const,msigbg,useWidth)

     write(22,'(99g12.5)') srts,msigbg,sum(msigbg)

  end do


end program printMesMes
