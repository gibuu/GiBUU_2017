program Check_gammaN_rhoN

  use inputGeneral
  use output
  use version
  use particleDefinition
  use particleProperties
  use hadronFormation
  use CollTools
  use VMMassPythia
  use mediumDefinition

  use photonXSections
  use XS_VMD

  IMPLICIT NONE

  type(medium) :: med
  real, dimension(1:4) :: sig
  real, dimension(0:4) :: sigma

  real, dimension(1:4) :: XS

  integer :: iSrts,nSrts
  real :: Srts,Srts1,Srts2,dSrts
  integer :: i


  call PrintVersion
  call readInputGeneral

  call forceInitFormation
  call InitParticleProperties

  Srts1 = 1.6
  Srts2 = 4.2
  dSrts = 0.02
  nSrts = (Srts2-Srts1)/dSrts

  do iSrts=0,nSrts
     Srts = Srts1+iSrts*dSrts

     do i=1,3
       call setIParam(i) 
       call calcXS_gammaN2VN(srts,med,sig)
       XS(i) = sig(1)
    end do

    call vmd(srts,0.0,0.0,sigma)
    XS(4) = sigma(1)

    write(121,'(f7.4,1P,20e13.5)') srts,XS

  end do





end program Check_gammaN_rhoN
