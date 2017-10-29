program test
  use output
  use particleDefinition
  use particleProperties
  use Coll_BaB

  implicit none


  type(particle),dimension(2):: inPart   ! incoming particles
  type(particle),dimension(30):: outPart  ! outgoing particles
  real                       :: srtS
  real, dimension(0:3)       :: pcm
  real, dimension(1:3)       :: beta
  logical                    :: flagOK
  
  call InitParticleProperties
  
  inPart%ID = (/34, 2 /)
  inPart%charge = (/1, -2 /)
  inPart%antiparticle = (/.False., .TRUE. /)
!  inPart%charge = (/-1, 2 /)
!  inPart%antiparticle = (/.TRUE.,.FALSE. /)
 
  
!  inPart%ID = (/1, 2 /)
!  inPart%charge = (/1, -2 /)
!  inPart%antiparticle = (/.False., .TRUE. /)

!  inPart%ID = (/34, 34 /)
!  inPart%charge = (/2, -2 /)
!  inPart%antiparticle = (/.False., .TRUE. /)


  beta = 0
  pcm = (/1, 0, 0, 0/)
  srtS = 2.732

  call WriteParticle(6,99,1,inPart(1))
  call WriteParticle(6,99,2,inPart(2))
  


  call DoColl_BaB(inPart,outPart,flagOK, srtS,pcm,beta)
  call DoColl_BaB(inPart,outPart,flagOK, srtS,pcm,beta)

  write(*,*) 'ok'
  
  



end program test
