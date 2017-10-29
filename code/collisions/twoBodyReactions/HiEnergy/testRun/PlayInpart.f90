program test
  use output
  use particleDefinition
  use particleProperties
  use Coll_Fritiof
  use Coll_Pythia
  use VMMassParameter, only : srtFreeVMMass
  use hadronFormation, only : forceInitFormation
  use ID_translation, only: KFfromBUU

  implicit none


  type(particle),dimension(2):: inPart   ! incoming particles
  type(particle),dimension(30):: outPart  ! outgoing particles
  real                       :: srtS
  real, dimension(0:3)       :: pcm
  real, dimension(1:3)       :: beta
  logical                    :: flagOK
  
  real                       :: weight

  integer :: i

  integer :: id1,id2,iz1,iz2,kf

  call forceInitFormation
  call InitParticleProperties

  beta = 0
  pcm = (/1, 0, 0, 0/)
  srtS = 10.0

!  do id1=1,61
!     do iz1=-1,2
  do id1=101,113
     do iz1=-1,1
        call KFfromBUU(KF, ID1,IZ1)
        if (KF==0) cycle

        do id2=1,61
           do iz2=-1,2

              call KFfromBUU(KF, ID2,IZ2)
              if (KF==0) cycle

              write(*,*) id1,iz1,id2,iz2
!==========
              inPart%ID     = (/ID1,ID2/)

              inPart%charge = (/IZ1,IZ2/)
              inPart%antiparticle = (/.False., .False. /)

!              inPart%charge = (/IZ1,-IZ2/)
!              inPart%antiparticle = (/.False., .TRUE. /)

              call DoColl_Pythia(inPart,outPart,flagOK, srtS,pcm,beta,.FALSE.,weight)

!==========

           end do
        end do
     end do
  end do


!  call DoPythia

  
!  inPart%ID = (/34, 2 /)
!  inPart%charge = (/1, -2 /)
!  inPart%antiparticle = (/.False., .TRUE. /)
!  inPart%charge = (/-1, 2 /)
!  inPart%antiparticle = (/.TRUE.,.FALSE. /)
 
  
!  inPart%ID = (/1, 1 /)
!  inPart%charge = (/1, 0 /)
!  inPart%antiparticle = (/.False., .False.. /)

!  inPart%ID = (/34, 34 /)
!  inPart%charge = (/2, -2 /)
!  inPart%antiparticle = (/.False., .TRUE. /)


  inPart%ID = (/110, 34 /)
  inPart%charge = (/1, -1 /)
  inPart%antiparticle = (/.False., .False. /)

  inPart%mass = (/0.496, 1.371/)

  beta = 0
  pcm = (/1, 0, 0, 0/)
!  srtS = 3.732
!  srtS = 20.000
  srtS = 2.3466

  srtFreeVMMass=srts


  call WriteParticle(6,99,1,inPart(1))
  call WriteParticle(6,99,2,inPart(2))
  

  OutPart%productionTime = 0.0
  OutPart%in_Formation = .TRUE.

  call DoColl_Pythia(inPart,outPart,flagOK, srtS,pcm,beta,.FALSE.,weight)

  call QYLIST(2)
  
  do i=1,30
     if (OutPart(i)%ID > 0) call WriteParticle(6,99,i,outPart(i))
  enddo


!  call DoColl_Fritiof(inPart,outPart,flagOK, srtS,pcm,beta)
!  call DoColl_Fritiof(inPart,outPart,flagOK, srtS,pcm,beta)

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


