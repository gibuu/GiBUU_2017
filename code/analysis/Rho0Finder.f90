!******************************************************************************
!****m* /Rho0Finder
! NAME
! module Rho0Finder
!
! PURPOSE
! This module includes routines, which try to identify rho0 from
! charged particles in a particle vetor.
!
! INPUTS
! (none)
!******************************************************************************
module Rho0Finder
  use particlePointerListDefinition
  use particlePointerList

  implicit none

  private

  public :: FindRho0

  logical,  save :: initFlag = .true.

  !--- steering parameters:
  logical, save :: TwoTracks = .false.
  logical, save :: OnlyPion  = .true.
  logical, save :: TotChargeZero = .true.
  logical, save :: UseProbs  = .false.
  logical, save :: UseProbsBool  = .false.
  !---

contains

  !****************************************************************************
  !****s* Rho0Finder/readInput
  ! NAME
  ! subroutine readInput
  !
  ! PURPOSE
  ! Read in namelist
  !****************************************************************************
  subroutine readInput
    use output

    implicit none
    integer :: ios

    NAMELIST /Rho0Finder/ TwoTracks,OnlyPion,TotChargeZero,UseProbs,useProbsBool

    if (.not.initFlag) return

    call Write_ReadingInput('Rho0Finder',0)
    rewind(5)
    read(5,nml=Rho0Finder,IOSTAT=ios)
    call Write_ReadingInput('Rho0Finder',0,ios)

    write(*,*) 'TwoTracks     = ',TwoTracks
    write(*,*) 'OnlyPion      = ',OnlyPion
    write(*,*) 'TotChargeZero = ',TotChargeZero
    write(*,*) 'UseProbs      = ',UseProbs
    write(*,*) 'UseProbsBool  = ',UseProbsBool

    call Write_ReadingInput('Rho0Finder',1)
    initFlag=.false.

  end subroutine readInput



  !****************************************************************************
  !****s* Rho0Finder/FindRho0
  ! NAME
  ! subroutine FindRho0(PartsIn,PartsOut,momIn,ProbIn)
  !
  ! PURPOSE
  ! This routine implements a rho0 finding process as in experiments.
  ! Given the list of detected particles in a specific event, it returns
  ! a particle list of (possible) reconstructed rho0.
  !
  ! If the particles have a detection probability less than 1, the
  ! probability of the rho0 is the product of the two decay candidates.
  ! One has the option to multiply this probability
  ! with 1 minus the detection probability of every other
  ! particle.
  !
  ! INPUTS
  ! * type(tParticleList) :: PartsIn -- particles produced in the event
  ! * real, dimension(0:3) :: momIn -- momentum used for calculation of
  !   missing mass, e.g. p+q, i.e. photon + proton momenta
  ! * real, dimension(:), OPTIONAL :: ProbIn -- Array of detection
  !   probabilities (AccWeight)
  !
  ! OUTPUT
  ! * type(tParticleList) :: PartsOut -- list of possible rho0
  !
  ! NOTES
  ! * At the moment, this routine just implements the method used by
  !   the Hermes experiment, cf. M.Tytgat, PhD thesis (Diffractive
  !   production of rho0 and omega mesons at Hermes), p.116f.
  !   A more general definition with some steering parameters (maybe
  !   via jobCard) is anticipated.
  ! * This routine allocates memory for the rho0. Therefore it is necessary
  !   that you clean up the particle list in the calling routine after use
  !   via
  !        call ParticleList_CLEAR(PartsOut,.true.)
  ! * In order to complete the Tytgat analysis, an additional cut on -t is
  !   necessary
  ! * In order to get comparable results as in experiments, Kaons should
  !   be instable in calculations and decay into pions.
  ! * We abuse many of the variables of the type to store information
  !   different to their intention. We list them in the following:
  ! * Abuse: x-position --> missing mass M_X
  ! * Abuse: y-position --> Delta E
  ! * Abuse: z-position --> DecayTime
  ! * Abuse: scaleCS --> probability, that the pions really come out
  !   of a decay (1.0, if both pions are stored in
  !   PIL_rho0Dec and 0.5, if only one pion is stored and 0, if we
  !   constructed a rho0, which does not stem from a decay [or we forgot
  !   to switch on storeRho0Info in collisionterm]).
  ! * Abuse: offShellParameter --> Decay angle theta
  !   (cf. CalculateDecayAngle)
  !****************************************************************************

  subroutine FindRho0(PartsIn,PartsOut,momIn,ProbIn)

    use PIL_rho0Dec
    use CallStack
    use output
    use constants, only: mN

    implicit none
    type(tParticleList), intent(IN) :: PartsIn
    real, dimension(0:3), intent(IN) :: momIn
    type(tParticleList), intent(OUT) :: PartsOut
    real, dimension(:), intent(IN), optional :: ProbIn

    real, dimension(:), allocatable :: Prob
    type(particle), POINTER :: Part

    integer :: nParts, i1,i2,i3
    real :: probRho, mRho2, mX2, DeltaE
    real, dimension(0:3) :: momRho,momX
    type(tParticleListNode), POINTER :: pNode1,pNode2,pNode3
    integer :: nr2

    real :: tC,tP,tF
    real :: thetaDecay

    if (initFlag) call readInput

    call ParticleList_INIT(PartsOut)

    nParts=PartsIn%nEntries
    if (nParts.eq.0) return

    allocate(Prob(nParts))
    Prob = 1.0

    if ((useProbs).and.(present(ProbIn))) then
       if (size(ProbIN).lt.nParts) then
          write(*,*) 'ERROR:size(ProbIn) too small!',size(ProbIN),nParts
          stop
       end if
       Prob = ProbIn(1:nParts)
       if (useProbsBool) then
          do i1=1,nParts
             if (Prob(i1).gt.0) Prob(i1)=1.0
          end do
       end if
    end if

    pNode1 => PartsIn%first
    i1 = 1
    do while (associated(pNode1))

       if (Prob(i1).eq.0.0) goto 101
       if (pNode1%V%charge .eq.0) goto 101
       if ((OnlyPion).and.(pNode1%V%ID.ne.101)) goto 101


       pNode2 => pNode1%next
       i2 = i1+1
       do while (associated(pNode2))

          if (Prob(i2).eq.0.0) goto 102
          if (pNode2%V%charge .eq.0) goto 102
          if ((OnlyPion).and.(pNode2%V%ID.ne.101)) goto 102

          if ((TotChargeZero).and.(pNode1%V%charge+pNode2%V%charge .ne.0)) &
               & goto 102

          probRho=Prob(i1)*Prob(i2)

          if (TwoTracks) then
             pNode3 => PartsIn%first
             i3 = 1
             do while (associated(pNode3))
                if (i3.eq.i1) goto 103
                if (i3.eq.i2) goto 103

                probRho=probRho*(1.0-Prob(i3))

103             pNode3 => pNode3%next
                i3 = i3+1
             end do ! Loop3
          end if

          if (ProbRho.eq.0.0) goto 102

          momRho = pNode1%V%momentum+pNode2%V%momentum
          mRho2 = momRho(0)**2-Sum(momRho(1:3)**2)

!          if (mRho2 .lt. 0.6**2) goto 102
!          if (mRho2 .gt. 1.0**2) goto 102

          momX = momIn-momRho
          mX2 = momX(0)**2-Sum(momX(1:3)**2)

          DeltaE=(mX2-mN**2)/(2*mN)

!          if (DeltaE.gt.0.6) goto 102

          call CalculateDecayAngle(thetaDecay)


          allocate(Part)
          Part%momentum = momRho
          Part%mass = sqrt(mRho2)
          Part%perWeight = pNode1%V%perWeight * probRho
          Part%ID = 103
          Part%charge = 0
          Part%event = pNode1%V%event
          Part%firstEvent = pNode1%V%firstEvent
          Part%position(1) = sqrt(max(0.0,mX2)) ! abuse !!!
          Part%position(2) = DeltaE    ! abuse !!!

          Part%offShellParameter=thetaDecay ! abuse !!!

          Part%scaleCS = 0.0 ! abuse !!!
          if (PIL_rho0Dec_GET(pNode1%V%number,nr2,tC,tP,tF)) then
             if (pNode2%V%number .eq. nr2) then
                Part%scaleCS = 1.0
             else
                Part%scaleCS = 0.5
             end if
          else if (PIL_rho0Dec_GET(pNode2%V%number,nr2,tC,tP,tF)) then
             if (pNode1%V%number .eq. nr2) then
                write(*,*) 'oops, why not earlier?'
                call TRACEBACK()
                Part%scaleCS = 1.0
             else
                Part%scaleCS = 0.5
             end if
          end if

          Part%lastCollisionTime=-99.9
          Part%productionTime   =-99.9
          Part%formationTime    =-99.9
          Part%position(3) = 0    ! abuse !!!
          if (Part%scaleCS.eq.1.0) then
             Part%lastCollisionTime= tC
             Part%productionTime   = tP
             Part%formationTime    = tF
             Part%position(3) = pNode2%V%productionTime    ! abuse !!!
          end if

          call ParticleList_APPEND(PartsOut,Part)

102       pNode2=>pNode2%next
          i2 = i2+1
       end do ! Loop2

101    pNode1=>pNode1%next
       i1 = i1+1
    end do ! Loop1

    deallocate(Prob)

  contains

    !**************************************************************************
    !****s* FindRho0/CalculateDecayAngle
    ! NAME
    ! subroutine CalculateDecayAngle(theta)
    ! PURPOSE
    ! Calculate the decay angle theta in the "s-channel helicity system", i.e.
    ! the angle between the positive z-axis and the pi+ momentum
    ! in the rest frame of the rho (-z axis defined by recoil particle)
    !
    ! The angle phi (i.e. the azimuthal angle of the pi+ in respect of the
    ! photon-recoil plane) is not calculated
    ! INPUTS
    ! * momIn
    ! * pNode1%V%momentum,pNode2%V%momentum
    ! * momRho
    ! OUTPUT
    ! * real :: theta -- the angle (in degrees)
    !**************************************************************************
    subroutine CalculateDecayAngle(theta)
!       use rotation, only : get_phi_Theta, rotateZY
      use lorentzTrafo, only: lorentz, lorentzCalcBeta
      use vector, only: theta_in

      real, intent(out) :: theta

      real, dimension(0:3) :: momRecoil,momPi
      real, dimension(3) :: beta
!      real :: theta_rot,phi_rot

      momRecoil = momIn - momRho

      if (pNode1%V%charge.gt.0) then
         momPi = pNode1%V%momentum
      else
         momPi = pNode2%V%momentum
      end if

      beta = lorentzCalcBeta (momRho)
      call lorentz(beta,momRecoil)
      call lorentz(beta,momPi)

      ! if we only calculate theta, we do not have to rotate
      ! the recoil momentum to be the -z didrection.
      ! we directly calculate theta as the angle between
      ! -recoil and pi+

!      call get_phi_Theta(-momRecoil(1:3),theta_Rot,phi_rot)
!      call rotateZY(theta_Rot,phi_Rot,momRecoil(1:3),momRecoil(1:3))
!      call rotateZY(theta_Rot,phi_Rot,momPi(1:3),momPi(1:3))

      theta = theta_in(-momRecoil(1:3),momPi(1:3))

!      write (*,*) 'theta:',theta
!      write (*,*)

    end subroutine CalculateDecayAngle

  end subroutine FindRho0


end module Rho0Finder
