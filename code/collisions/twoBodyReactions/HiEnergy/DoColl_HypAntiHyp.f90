!******************************************************************************
!****m* /Coll_HypAntiHyp
! PURPOSE
! Implement Baryon + AntiBaryon -> Hyperon + AntiHyperon  processes (for the HiEnergy part)
!******************************************************************************
module Coll_HypAntiHyp

  IMPLICIT NONE
  private

  public :: DoColl_YYbar

contains


  !****************************************************************************
  !****s* Coll_HypAntiHyp/DoColl_YYbar
  ! NAME
  ! subroutine DoColl_YYbar(inPart,outPart,flagOK,sqrtS,pcm,beta,outChannel)
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
  ! INPUTS
  ! * type(particle),dimension(:) :: inPart   -- incoming particles
  ! * real                        :: sqrtS    -- energy of collision
  ! * real, dimension(0:3)        :: pcm      -- c.m. momentum of 1-st particle
  ! * real, dimension(1:3)        :: beta     -- boost vector
  ! * integer,                    :: outChannel  -- outgoing channel = 1 - LambdaBar+Lambda
  ! *                                                                = 2 - Lambda+SigmaBar or LambdaBar+Sigma
  ! *                                                                = 3 - Xi+XiBar
  !
  ! OUTPUT
  ! * type(particle),dimension(:) :: outPart  ! outgoing particles
  ! * logical                     :: flagOK   ! event okay ?
  !
  ! NOTES
  ! This routine is applicable for baryon-antibaryon collisions only.
  !****************************************************************************

  subroutine DoColl_YYbar(inPart,outPart,flagOK,sqrtS,pcm,beta,outChannel)

    use constants, only: twoPi
    use particleDefinition, only: particle
    use particleProperties, only: hadron
    use random, only: rn
    use IdTable, only: Lambda, SigmaResonance, Xi
    use RMF, only: getRMF_flag
    use winkel_tools, only: dsigdt_Regge

    type(particle),dimension(:),intent(in)     :: inPart   ! incoming particles
    type(particle),dimension(:),intent(inout)  :: outPart  ! outgoing particles
    real,                       intent(in)   :: sqrtS
    real, dimension(0:3),       intent(in)   :: pcm
    real, dimension(1:3),       intent(in)   :: beta
    integer,                    intent(in)   :: outChannel
    logical,                    intent(out)  :: flagOK

    real :: x, s, prcm, pcmOut
    real :: ctheta, stheta, phi, cphi, sphi ! random variables
    real :: h1,h2,h3
    real :: phiB, thetaB ! boost angles
    integer :: nAnti, totCharge, i, j

    double precision MP_P ! prototype


    if (inPart(1)%antiparticle.eqv.inPart(2)%antiparticle) then
       write(*,*) 'DoColl_YYbar can not be used for Bar+Bar or AntiBar+Antibar collisions !'
       stop
    end if

    flagOK = .TRUE.

    outPart(1:2)%antiparticle=inPart(1:2)%antiparticle
    if ( inPart(1)%antiparticle ) then
       nAnti=1
    else
       nAnti=2
    end if

    totCharge= sum(inPart(1:2)%charge)

    s=sqrtS**2

    if ( .not.getRMF_flag() ) then
       prcm= sqrt(dot_product(pcm(1:3),pcm(1:3)))
    else
       prcm = ( s + inPart(1)%mass**2 - inPart(2)%mass**2 )**2/(4.*s) &
               & - inPart(1)%mass**2
       prcm = sqrt( max(0.,prcm) )
    end if


    select case (outChannel)

    case (1)   ! LambdaBar+Lambda

       if (totCharge /= 0) then
          write(*,*) 'In DoColl_YYbar: wrong charge of initial state'
          write(*,*) 'for Lambda+LambdaBar production : ', totCharge
          stop
       end if

       outPart(1:2)%Id=Lambda
       outPart(1:2)%charge=0
       outPart(1:2)%mass=hadron(Lambda)%mass

       ctheta = dsigdt_Regge (sqrtS, prcm, 1)

    case (2)   ! Lambda+SigmaBar or LambdaBar+Sigma

       if (abs(totCharge) > 1) then
          write(*,*) 'In DoColl_YYbar: wrong charge of initial state'
          write(*,*) 'for  Lambda SigmaBar / LambdaBar Sigma production : ', totCharge
          stop
       end if

       x = rn()
       if (x < 0.5) then
          outPart(nAnti)%ID=Lambda  !LambdaBar Sigma
          outPart(nAnti)%charge=0
          outPart(nAnti)%mass=hadron(Lambda)%mass
          outPart(3-nAnti)%ID=SigmaResonance
          outPart(3-nAnti)%charge=totCharge
          outPart(3-nAnti)%mass=hadron(SigmaResonance)%mass
       else
          outPart(nAnti)%ID=SigmaResonance   !Lambda SigmaBar
          outPart(nAnti)%charge=totCharge
          outPart(nAnti)%mass=hadron(SigmaResonance)%mass
          outPart(3-nAnti)%ID=Lambda
          outPart(3-nAnti)%charge=0
          outPart(3-nAnti)%mass=hadron(Lambda)%mass
       end if

       ctheta = dsigdt_Regge (sqrtS, prcm, 2)

    case (3)   ! Xi+XiBar

       outPart(1:2)%Id=Xi
       if ( totCharge == 0 ) then
          x = rn()
          if (x < 0.5) then                ! Xi^- + XiBar^+
             outPart(nAnti)%charge=1
          else                             ! Xi^0 + XiBar^0
             outPart(nAnti)%charge=0
          end if
       else if ( totCharge == -1) then     ! Xi^- + XiBar^0
          outPart(nAnti)%charge=0
       else if ( totCharge == 1) then      ! Xi^0 + XiBar^+
          outPart(nAnti)%charge=1
       else
          write(*,*) 'In DoColl_YYbar: wrong charge of initial state'
          write(*,*) 'for Xi+XiBar production : ', totCharge
          stop
       end if

       outPart(3-nAnti)%charge= totCharge - outPart(nAnti)%charge

       outPart(1:2)%mass=hadron(Xi)%mass

       ctheta = 2*rn()-1   ! isotropic angular distribution

    end select

    if (sqrtS < sum(outPart(1:2)%mass)) then
       flagOK=.false.
       outPart%Id=0
       return
    end if

! c.m. momentum of outgoing particles:
    pcmOut=sqrt((s-outPart(1)%mass**2+outPart(2)%mass**2)**2/(4.*s)-outPart(2)%mass**2)

! set random variables: phi

    phi = twopi*rn()

! set other variables:

    stheta = sqrt(1.0-ctheta**2)

    cphi = cos(phi)
    sphi = sin(phi)

    h1 = pcmOut * stheta*cphi
    h2 = pcmOut * stheta*sphi
    h3 = pcmOut * ctheta

! boost to final system:

    call MP_Set3(1, outPart(1)%mass,  h1,  h2,  h3 )
    call MP_Set3(2, outPart(2)%mass, -h1, -h2, -h3 )

    phiB   = atan2(pcm(2),pcm(1))
    thetaB = atan2(sqrt(pcm(1)**2+pcm(2)**2),pcm(3))

    call MP_ROBO(1,2, thetaB, phiB, beta(1),beta(2),beta(3) )

    do i=1,2                   ! set momenta according scattering
       do j=1,3
          outPart(i)%momentum(j) = real(MP_P(i,j))
       end do

       outPart(i)%momentum(0) = outPart(i)%mass
       outPart(i)%momentum(0) = sqrt(DOT_PRODUCT(outPart(i)%momentum,outPart(i)%momentum))
    end do

    outPart(1:2)%number = 0

  end subroutine DoColl_YYbar

end module Coll_HypAntiHyp
