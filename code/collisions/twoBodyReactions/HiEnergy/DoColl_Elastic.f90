!******************************************************************************
!****m* /Coll_Elastic
! PURPOSE
! Implement a + b -> X elastic processes/reactions (for the HiEnergy part)
!******************************************************************************
module Coll_Elastic

  private

  public :: DoColl_Elast

contains

  !****************************************************************************
  !****s* Coll_Elastic/DoColl_Elast
  ! NAME
  ! subroutine DoColl_Elast(inPart,outPart,flagOK,sqrtS,pcm,beta,AngDist)
  !
  ! PURPOSE
  ! Perform a collision of particles given in "inPart" with energy "sqrtS" and
  ! return outgoing particles in "outPart".
  !
  ! "pcm" and "beta" are vectors used for Boost and Rotation of the event.
  !
  ! if "flagOK" is false, no event happened, the output in "outPart" should
  ! be neglected!
  !
  ! Different angular (tSlope) distributions are possible:
  ! * isotropic
  ! * tSlope_Effenb(sqrts)     --- pp elastic scattering parameterization from
  !   J. Cugnon et al., NPA 352, 505 (1981)
  ! * tSlope_PYTHIA(sqrts, ID1, ID2) --- taken from PYTHIA model, see
  !   T. Falter et al., PRC 70, 054609 (2004)
  !
  ! INPUTS
  ! * type(particle),dimension(:) :: inPart   -- incoming particles
  ! * real                        :: sqrtS    -- energy of ollision
  ! * real, dimension(0:3)        :: pcm      -- boost-vector
  ! * real, dimension(1:3)        :: beta     -- boost vector
  ! * integer,                    :: AngDist  -- angular distribution
  !
  ! OUTPUT
  ! * type(particle),dimension(:) :: outPart  ! outgoing particles
  ! * logical                     :: flagOK   ! event okay ?
  !
  ! NOTES
  ! cf. DoColl_Pythia
  !
  ! in order to understand the meaning of "pcm" and "beta":
  ! The (Pythia-)event is done in the restframe of the two given particles.
  ! Then a call to PYROBO according
  !       phi = atan2(pcm(2),pcm(1))
  !       theta = atan2(sqrt(pcm(1)**2+pcm(2)**2),pcm(3))
  !       call PYROBO(1,N, theta,phi, beta(1),beta(2),beta(3))
  ! is performed in order to transform the system into the desired
  ! (Lab-) system.
  !****************************************************************************
  subroutine DoColl_Elast(inPart,outPart,flagOK,sqrtS,pcm,beta,AngDist)

    use constants
    use particleDefinition, only: particle
    use random
    use RMF, only: getRMF_flag
    use IdTable, only: nucleon
    use winkel_tools

    IMPLICIT NONE

    type(particle),dimension(:),intent(in)   :: inPart   ! incoming particles
    type(particle),dimension(:),intent(inout):: outPart  ! outgoing particles
    real,                       intent(in)   :: sqrtS
    real, dimension(0:3),       intent(in)   :: pcm
    real, dimension(1:3),       intent(in)   :: beta
    integer,                    intent(in)   :: AngDist
    logical,                    intent(out)  :: flagOK

    real :: pA, pcmL
    real :: ctheta, stheta, phi, cphi, sphi ! random variables
    real :: h1,h2,h3,h4
    real :: phiB, thetaB ! boost angles
    integer :: i,j
    logical :: flag

    double precision MP_P ! prototype


    flagOK = .TRUE.

    ! set momentum: pA and sqrt(pcm(1)**2+pcm(2)**2+pcm(3)**2)
    ! shoud be equal in vacuum,
    ! but in RMF mode pA and sqrt(pcm(1)**2+pcm(2)**2+pcm(3)**2)
    ! are not equal, since pcm(0:3) is the kinetic four-momentum.
    ! Therefore, in RMF mode we have to use pA as the analog
    ! of the vacuum c.m. momentum.

    pA = 0.5/sqrtS * sqrt( (sqrtS**2-(inPart(1)%mass+inPart(2)%mass)**2) &
                           & * (sqrtS**2-(inPart(1)%mass-inPart(2)%mass)**2) )

    if ( .not.getRMF_flag() ) then
       pcmL = sqrt(pcm(1)**2+pcm(2)**2+pcm(3)**2)
    else
       pcmL = pA
    end if

!    if (abs(pA-pcmL)/pA.gt.1e-3) then
!       write(*,*) '(DoColl_Elast) Warning: abs(pA-|pcm|)/pA>1e-3 !',pA,pcmL
!    endif


    ! set random variables: ctheta, phi

    if (AngDist.eq.1) then
       ctheta = 2*rn()-1
    else

       if ( inPart(1)%Id.eq.nucleon .and. inPart(2)%Id.eq.nucleon &
            & .and. (inPart(1)%antiParticle.neqv.inPart(2)%antiParticle) ) then
          h1 = tSlope_EL_pbarp(sqrts)
       else
          if (AngDist.eq.2) then
             h1 = tSlope_Effenb(sqrts)
          else if (AngDist.eq.3) then
             h1 = tSlope_Pythia(sqrts, inPart(1)%ID,inPart(2)%ID)
          else
             write(*,*) '(COLL_El): wrong AngDist! stop!'
             stop
          end if
       end if

       h1 = -2*pcmL**2 * h1   ! -dt/d(ctheta)
       h3 = rn()             ! RAN

       if (h1 < -100) then
          h4 = log(0.+h3)/h1
       else
          h4 = log((1.0-h3)*exp(2*h1)+h3)/h1
                              ! = RAN(0...4*pA^2) accord. exp(-A*x)
       end if
       ctheta = 1.0-h4
       if (abs(ctheta).gt.1.0) then
          ctheta = 2*h3-1
       end if
    end if

    phi = twopi*rn()

! set other variables:

    stheta = sqrt(1.0-ctheta**2)

    cphi = cos(phi)
    sphi = sin(phi)

    h1 = pcmL * stheta*cphi
    h2 = pcmL * stheta*sphi
    h3 = pcmL * ctheta

! boost to final system:

    call MP_Set3(1, inPart(1)%mass,  h1,  h2,  h3 )
    call MP_Set3(2, inPart(2)%mass, -h1, -h2, -h3 )

    phiB   = atan2(pcm(2),pcm(1))
    thetaB = atan2(sqrt(pcm(1)**2+pcm(2)**2),pcm(3))

    call MP_ROBO(1,2, thetaB, phiB, beta(1),beta(2),beta(3) )

! set outgoing particles:

    flag = outPart(1)%perturbative ! remember this flag

    outPart%ID = 0 ! reset all particles

    outPart(1:2) = inPart(1:2) ! set all variables to the ingoing values
    outPart(1:2)%number = 0
    outPart(1:2)%perturbative = flag

    do i=1,2                   ! set momenta according scattering
       do j=1,3
          outPart(i)%momentum(j) = real(MP_P(i,j))
       end do

       outPart(i)%momentum(0) = outPart(i)%mass
       outPart(i)%momentum(0) = sqrt(DOT_PRODUCT(outPart(i)%momentum,outPart(i)%momentum))
    end do

  end subroutine DoColl_Elast

end module Coll_Elastic
