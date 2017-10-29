!******************************************************************************
!****m* /initElementary
! NAME
! module initElementary
! PURPOSE
! Implement the init routines for collisions of "elementary"
! (non-perturbative) particles.
!******************************************************************************
module initElementary

  use particleDefinition, only: particle

  IMPLICIT NONE

  private

  !****************************************************************************
  !****g* elementary/impactParameter
  ! SOURCE
  !
  real, save :: impactParameter = -1. ! [fm]
  ! PURPOSE
  ! * >=0: this is the actual impact parameter
  ! * <0 : Impact parameter integration up to an automatically determined b_max.
  !   The actual impact parameter is randomly sampled in the interval
  !   [0.,b_max] with a proper geometrical weight.
  !****************************************************************************

  !****************************************************************************
  !****g* elementary/particleId
  ! SOURCE
  !
  integer,dimension(2), save :: particleId=(/1,1/)
  ! PURPOSE
  ! Id of particles
  !****************************************************************************

  !****************************************************************************
  !****g* elementary/particleAnti
  ! SOURCE
  !
  logical,dimension(2), save :: particleAnti=(/.false.,.false./)
  ! PURPOSE
  ! if .true. then particles are antiparticles
  !****************************************************************************

  !****************************************************************************
  !****g* elementary/particleCharge
  ! SOURCE
  !
  integer,dimension(2), save :: particleCharge=(/0,0/)
  ! PURPOSE
  ! Charge of the particles
  !****************************************************************************

  !****************************************************************************
  !****g* elementary/particleMass
  ! SOURCE
  !
  real,   dimension(2), save :: particleMass=(/-1.0,-1.0/)
  ! PURPOSE
  ! Mass of particles (if <0, then mass is set to default values)
  !****************************************************************************

  !****************************************************************************
  !****g* elementary/srtsRaiseFlag
  ! SOURCE
  !
  logical, save :: srtsRaiseFlag=.false.
  ! PURPOSE
  ! * if .true. then the srts stepping is done
  ! * if .false. then the ekin_lab stepping is done
  !****************************************************************************

  !****************************************************************************
  !****g* elementary/ekin_lab
  ! SOURCE
  !
  real, save :: ekin_lab=1. ! [GeV]
  ! PURPOSE
  ! kin. energy of first particle in the rest frame of second particle
  ! (starting value for the energy scan: the number of different energies
  ! is set by parameter num_Energies in the namelist "input")
  !****************************************************************************

  !****************************************************************************
  !****g* elementary/delta_ekin_lab
  ! SOURCE
  !
  real, save :: delta_ekin_lab=0.03 ! [GeV]
  ! PURPOSE
  ! kin. energy step for energy scan
  !****************************************************************************

  !****************************************************************************
  !****g* elementary/srts
  ! SOURCE
  !
  real, save :: srts=3. ! [GeV]
  ! PURPOSE
  ! invariant energy
  ! (starting value for the energy scan)
  !****************************************************************************

  !****************************************************************************
  !****g* elementary/delta_srts
  ! SOURCE
  !
  real, save :: delta_srts=1. ! [GeV]
  ! PURPOSE
  ! srts step for srts scan
  !****************************************************************************

  !****************************************************************************
  !****g* elementary/debug
  ! SOURCE
  !
  logical, parameter :: debug=.false.
  ! PURPOSE
  ! if .true., additional printouts are done for debugging
  !****************************************************************************

  !Working variables:
  real, save :: p_lab                      ! Laboratory momentum (of 1-st particle in the r.f. of 2-nd)
  real, save :: siggeo                     ! Geometrical cross section (mbarn)
  type(particle), save :: projectile, target ! 1-st and 2-nd particles, respectively

  public :: initElementaryCollision
  public :: impactParameter, srts, p_lab, particleID, target, projectile, siggeo, &
            particleCharge, particleAnti, ekin_lab, particlemass

contains

  !****************************************************************************
  !****s* initElementary/initElementaryCollision
  ! NAME
  ! subroutine initElementaryCollision(teilchen,energyRaiseFlag)
  ! PURPOSE
  ! Provides initial conditions for the collision of two elementary
  ! particles.
  ! INPUTS
  ! * type(particle),dimension(:,:),intent(inout) :: teilchen
  !   -- array of particles to store the two colliding elementary particles.
  ! * logical, intent(in) :: energyRaiseFlag
  !   -- if .true., then the lab. energy of the first particle is increased
  !   by delta_ekin_lab.
  ! NOTES
  ! The user has to provide the namelist 'elementary', which contains
  ! the impact parameter and lab. energy of the projectile and id's of
  ! colliding particles. If the input impact parameter is negative, then
  ! the actual impact parameter is chosen by Monte-Carlo between
  ! 0 and abs(input impact parameter).
  !****************************************************************************
  subroutine initElementaryCollision(teilchen,energyRaiseFlag)

    use insertion, only: GarbageCollection

    type(particle),dimension(:,:),intent(inout) :: teilchen
    logical, intent(in) :: energyRaiseFlag

    logical, save :: initFlag=.true.
    integer :: j
    write(*,*)
    write(*,*) '**Initializing  elementary events'

    if (initFlag) call initInput

    if (energyRaiseFlag) then

      if (srtsRaiseFlag) then
        srts = srts + delta_srts
      else
        ekin_lab=ekin_lab+delta_ekin_lab
      end if

    end if

    do j=1,size(teilchen,dim=1) ! Loop over all Ensembles

       call setKinematics
       call setPosition

    end do

    projectile=Teilchen(1,1)
    target=Teilchen(1,2)

    call GarbageCollection(teilchen)

  contains

    !**************************************************************************
    !****s* initElementaryCollision/initInput
    ! NAME
    ! subroutine initInput
    ! PURPOSE
    ! Reads input out of jobcard. Namelist 'elementary'.
    !**************************************************************************
    subroutine initInput

      use output, only: Write_ReadingInput
      use IdTable, only: isHadron
      use ParticleProperties, only: hadron, PartName

      integer :: ios,i

      character*15, dimension(2) :: name

      !************************************************************************
      !****n* initElementary/elementary
      ! NAME
      ! NAMELIST /elementary/
      ! PURPOSE
      ! Includes parameters for initialization of a collision between two
      ! elementary particles:
      ! * impactParameter
      ! * particleId
      ! * particleAnti
      ! * particleCharge
      ! * srtsRaiseFlag
      ! * ekin_lab
      ! * delta_ekin_lab
      ! * srts
      ! * delta_srts
      !************************************************************************
      NAMELIST /elementary/  impactParameter, &
                             particleId, particleAnti, particleCharge, &
                             srtsRaiseFlag, ekin_lab, delta_ekin_lab, &
                             srts, delta_srts

      call Write_ReadingInput('Elementary',0)
      rewind(5)
      read(5,nml=elementary,IOSTAT=ios)
      call Write_ReadingInput('Elementary',0,ios)

      do i=1,2
         name(i) = PartName(particleID(i),particleCharge(i),particleAnti(i))

         if (isHadron(particleId(i))) then
            if (ParticleMass(i)<0.) ParticleMass(i)=hadron(particleId(i))%mass
         else
            write(*,*) ' In initElementary: particle ',i,' is wrong'
            stop
         end if

      end do

      write(*,*) '  Impact Parameter = ',impactParameter
      write(*,*) '  Particle: ID     = ',particleId
      write(*,*) '            Anti   = ',particleAnti
      write(*,*) '            charge = ',particleCharge
      write(*,*) '            mass   = ',particleMass
      write(*,*)
      write(*,*) '  >> ',trim(name(1)),' <-> ',trim(name(2)),' <<'
      write(*,*)
      write(*,*) '  srtsRaiseFlag=', srtsRaiseFlag
      if (srtsRaiseFlag) then
         write(*,*) '  srts      =', srts
         write(*,*) '  delta_srts=', delta_srts
      else
         write(*,*) '  Kinetic Energy of 1-st particle in lab frame=',ekin_lab
         write(*,*) '  Delta(Energy) for energy scans              =',delta_ekin_lab
      end if

      do i=1,2
         if (iand(hadron(particleId(i))%stability,1).ne.0) then
            write(*,'(A,i2,A)') 'Particle ',i,' is not stable during time steps. stop !!!!'
            write(*,'(A,i3,A,i2,A,i2,A)') 'you may consider setting "stabilityflag(',&
                 & particleId(i),')" [=',hadron(particleId(i))%stability,&
                 &'] to ',hadron(particleId(i))%stability-1,&
                 &' in the jobcard (namelist "ModifyParticles")'
            stop
         end if
      end do

      call Write_ReadingInput('Elementary',1)

      initFlag=.false.

    end subroutine initInput

    !**************************************************************************
    !****s* initElementaryCollision/setKinematics
    ! NAME
    ! subroutine setKinematics
    ! PURPOSE
    ! Sets basic kinematics of the elmentary particle collision.
    !**************************************************************************
    subroutine setKinematics
      use particleDefinition, only: setToDefault, setNumber

      !integer :: i
      !character*15, dimension(2) :: name

      call setToDefault(teilchen(j,1:2)) !set teilchen to its default values

      Teilchen(j,1:2)%ID=particleId(1:2)
      Teilchen(j,1:2)%antiparticle=particleAnti(1:2)
      Teilchen(j,1:2)%charge=particleCharge(1:2)
      Teilchen(j,1:2)%mass=particleMass(1:2)

      if (srtsRaiseFlag) then ! use srts for initialisation

        ekin_lab = ( srts**2 - (Teilchen(j,1)%mass+Teilchen(j,2)%mass)**2 ) &
                & / (2.*Teilchen(j,2)%mass)
      end if

      Teilchen(j,1)%momentum(0)=ekin_lab+Teilchen(j,1)%mass
      !1-st particle is initialized moving in positive z-direction:
      p_lab=Sqrt(Teilchen(j,1)%momentum(0)**2-Teilchen(j,1)%mass**2)
      Teilchen(j,1)%momentum(1:3)=(/0.,0.,p_lab/)
      Teilchen(j,2)%momentum(0)=Teilchen(j,2)%mass
      !2-nd particle is initialized at rest:
      Teilchen(j,2)%momentum(1:3)=0.

      if (.not.srtsRaiseFlag) then

          srts = Sqrt((Teilchen(j,1)%momentum(0)+Teilchen(j,2)%mass)**2 &
             & - p_lab**2)

      end if

      !assume vacuum dispersion relation:
      Teilchen(j,1)%velocity(1:3)=&
           & teilchen(j,1)%momentum(1:3)/teilchen(j,1)%momentum(0)
      Teilchen(j,2)%velocity(1:3)=0.

      teilchen(j,1)%event=1
      teilchen(j,2)%event=2

      ! Give the particles their unique numbers
      call setNumber(teilchen(j,1))
      call setNumber(teilchen(j,2))

      if (debug) then
         write(*,*) 'Ids:',teilchen(j,1:2)%Id
         write(*,*) 'Charges:',teilchen(j,1:2)%charge
         write(*,*) 'Antiparticles ?:',teilchen(j,1:2)%antiparticle
         write(*,*) 'Masses:',teilchen(j,1:2)%mass
         write(*,*) '4-momentum of 1-st particle',teilchen(j,1)%momentum
         write(*,*) 'Velocity of 1-st particle:',teilchen(j,1)%velocity
         write(*,*) '4-momentum of 2-nd particle',teilchen(j,2)%momentum
         write(*,*) 'Velocity of 2-nd particle:',teilchen(j,2)%velocity
       end if

    end subroutine setKinematics

    !**************************************************************************
    !****s* initElementaryCollision/setPosition
    ! NAME
    ! subroutine setPosition
    ! PURPOSE
    ! Sets positions of the elementary particles.
    ! NOTES
    ! * If impactParameter is choosen to be less than zero, then the
    !   impact parameter is choosen by a Monte-Carlo decision.
    ! * The second particle is always initialized at zero (already done
    !   by subroutine setToDefault called in setKinematics).
    !**************************************************************************
    subroutine setPosition

      use random, only: rn
      use inputGeneral, only: delta_T
      use XsectionRatios, only: getRatio
      use constants, only: pi

      real :: b, z, ratio
      logical, parameter :: AutoAdjustB = .true.

      z = -teilchen(j,1)%velocity(3)*delta_T

      b = abs(impactParameter)
      if (impactParameter < 0. .and. AutoAdjustB) b =  AdjustImpact(teilchen(j,1:2))

      ratio = getRatio((/teilchen(j,1),teilchen(j,2)/)) ! Medium modification factor
      siggeo = 10.*pi * b**2  * ratio

      if (impactParameter < 0.) b = b * sqrt(ratio) * sqrt(rn())

      teilchen(j,1)%position=(/b,0.,z/)

      if (debug) then
         write(*,*) 'Position of 1-st particle:',teilchen(j,1)%position
         write(*,*) 'Position of 2-nd particle:',teilchen(j,2)%position
      end if

    end subroutine setPosition

    !**************************************************************************
    real function AdjustImpact(pair)
      use particleDefinition, only: sqrts
      use mediumDefinition, only: vacuum
      use preEventDefinition
      use master_2Body, only: XsectionMaster
      use constants, only: pi

      type(particle),dimension(1:2)    :: pair
      real, dimension(0:7) :: sigs
      real :: srts
      real, dimension(0:3) :: momLRF
      type(preEvent), dimension(1:4)         :: finalState
      logical :: HiFlag

      srts=sqrts(pair(1),pair(2))
      momLRF=pair(2)%momentum+pair(1)%momentum

      call XsectionMaster(srts,pair,vacuum,momLRF,finalState,sigs,HiFlag)

      AdjustImpact = sqrt(sigs(0)/(10.*pi))

    end function AdjustImpact

  end subroutine initElementaryCollision



end module initElementary
