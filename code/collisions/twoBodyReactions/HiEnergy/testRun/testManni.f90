program test
  use output
  use particleDefinition
  use particleProperties
  use Coll_Manni
  use Coll_Pythia
  use VMMassParameter, only : srtFreeVMMass
  use hadronFormation, only : forceInitFormation

  use hist2Df90

  implicit none


  type(particle),dimension(2):: inPart   ! incoming particles
  type(particle),dimension(30):: outPart  ! outgoing particles
  real                       :: srtS
  real, dimension(0:3)       :: pcm
  real, dimension(1:3)       :: beta
  logical                    :: flagOK
  
  real                       :: weight, sigma, pL,pT
  character(3), dimension(-1:1), parameter :: piName  = (/'pi-','pi0','pi+'/)
  integer :: i,iEv,nEv

  type(histogram2D),dimension(-1:1) :: h2D

  call forceInitFormation
  call InitParticleProperties

!  call DoPythia

  
  nEv = 10000000

  inPart%ID = (/103, 1 /)
  inPart%charge = (/0, 0 /)
  inPart%antiparticle = (/.False., .False. /)
  inPart%mass = (/0.770, 0.938/)

  beta = 0
  pcm = (/1, 0, 0, 0/)
  srtS = 2.500


  srtFreeVMMass=srts


  call WriteParticle(6,99,1,inPart(1))
  call WriteParticle(6,99,2,inPart(2))

  do i=-1,1
     call CreateHist2D(h2D(i),"pT vs pL "//piName(i), &
            &(/-2.0,0.0/), (/2.0,2.0/), (/0.02,0.02/) , .true.)
  enddo

  OutPart%productionTime = 0.0
  OutPart%in_Formation = .TRUE.

  sigma = DoColl_ManniSigma(inPart,srtS)

  do iEv=1,nEv

     call DoColl_Manni(inPart,outPart,flagOK, srtS,pcm,beta)

!     call QYLIST(2)
  
     do i=1,30
        if (OutPart(i)%ID .eq. 1) then
           pL = OutPart(i)%momentum(3)
           pT = sqrt(OutPart(i)%momentum(1)**2+OutPart(i)%momentum(2)**2)
           call AddHist2D(h2D(OutPart(i)%charge),(/pL, pT/),1.0/pT,1.0)

        endif
     enddo

  end do

  do i=-1,1
     call WriteHist2D_Gnuplot(H2D(i),102+i,add=1e-20,mul=1./nEv)
  end do

  write(*,*) 'ok'


  
  
contains 
  subroutine DoPythia
    use VMMassParameter, only : srtFreeVMMass
    IMPLICIT NONE

    integer :: iEv

    write(*,*) 'START OF DOPYTHIA'

    call QYGIVE('MSTP(111) = 0') ! master switch fragmentation/decay
!    call QYGIVE('MSTP(125) = 2')


    call QYINIT('CMS','p','n', 20d0)
    srtFreeVMMass = 20d0

    do iEv=1,10
       call QYEVNT
       call QYLIST(2)
    end do
    
    stop

  end subroutine DoPythia
  

end program test


