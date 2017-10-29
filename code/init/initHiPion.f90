!******************************************************************************
!****m* /initHiPion
! NAME
! module initHiPion
!
! PURPOSE
! This module is for (high-energy) pion-nucleus and proton-nucleus collisions.
! In fact also any other hadron can be used as projectile.
! The eventtype 'HiPion' uses perturbative projectile particles.
!
! INPUTS
! The namelist "HiPionNucleus".
!
! NOTES
! following input parameters are changed/renamed/deleted/new:
! * numberPions -> nTestparticles
! * DoProton is removed
! * Charge -> ProjectileCharge
! * ProjectileID is introduced
! * ProjectileAnti is introduced
!******************************************************************************
module initHiPion

  use particleDefinition
  use idTable, only: pion
  implicit none

  private

  public :: initHiPionInduced
  public :: getTotalPerweight, getBetaFak
  public :: ConfigParticles

  !****************************************************************************
  !****g* initHiPion/ekin_lab
  ! SOURCE
  !
  real,save      :: ekin_lab = -99.9
  !
  ! PURPOSE
  ! Kinetic energy of projectile particles in lab frame [GeV].
  !****************************************************************************

  !****************************************************************************
  !****g* initHiPion/p_lab
  ! SOURCE
  !
  real,save      :: p_lab = -99.9
  !
  ! PURPOSE
  ! Momentum of projectile particles in lab frame [GeV/c].
  !****************************************************************************

  !****************************************************************************
  !****g* initHiPion/ProjectileCharge
  ! SOURCE
  !
  integer,save   :: ProjectileCharge = 0
  !
  ! PURPOSE
  ! Charge of projectile particles.
  !****************************************************************************

  !****************************************************************************
  !****g* initHiPion/ProjectileID
  ! SOURCE
  !
  integer,save   :: ProjectileID = pion
  !
  ! PURPOSE
  ! ID of projectile particles.
  !****************************************************************************

  !****************************************************************************
  !****g* initHiPion/ProjectileAnti
  ! SOURCE
  !
  logical,save   :: ProjectileAnti = .false.
  !
  ! PURPOSE
  ! Antiparticle flag of projectile particles.
  !****************************************************************************


  !****************************************************************************
  !****g* initHiPion/nTestparticles
  ! SOURCE
  !
  integer,save   :: nTestparticles = 200
  !
  ! PURPOSE
  ! Number of projectile testparticles per ensemble.
  !****************************************************************************

  !****************************************************************************
  !****g* initHiPion/minimumMomentum
  ! SOURCE
  !
  real,   save   :: minimumMomentum = 1.0
  !
  ! PURPOSE
  ! Minimal momentum of particles (in GeV) produced in the init routines.
  ! Only particles with an absolute momentum larger than this will be
  ! further propagated.
  !****************************************************************************

  !****************************************************************************
  !****g* initHiPion/useHermesPythiaPars
  ! SOURCE
  !
  logical,save :: useHermesPythiaPars = .false.
  !
  ! PURPOSE
  ! flag: Use "PYTHIA tuning done by HERMES collab"
  !****************************************************************************

  !****************************************************************************
  !****g* initHiPion/distance
  ! SOURCE
  !
  real,save      :: distance = 15.
  !
  ! PURPOSE
  ! Distance in z-direction from the nucleus center in fm,
  ! where the projectiles are initialzed.
  ! If negative, the distance will be chosen automatically.
  !****************************************************************************

  !****************************************************************************
  !****g* initHiPion/impact_parameter
  ! SOURCE
  !
  real,save      :: impact_parameter = 0.
  !
  ! PURPOSE
  ! Impact parameter of the projectiles in fm.
  ! If positive (or zero), this fixed value is used for all projectiles.
  ! If negative, the impact parameter is chosen by Monte Carlo,
  ! so that the projectiles are distributed over a certain disk.
  ! Cf. 'setPosition'.
  !****************************************************************************

  !****************************************************************************
  !****g* initHiPion/DoPerturbativeInit
  ! SOURCE
  !
  logical, save :: DoPerturbativeInit = .false.
  !
  ! PURPOSE
  ! If this flag is set to .true., the first collision of the projectile
  ! particles with a nucleon in the target nucleus will be done in this
  ! init routine (at timestep 0). This enables you to treat the first
  ! (hard) collision different from those in the FSI.
  !
  ! If this flag is set to .false., the projectile particles have to be
  ! propagated onto the nucleus as in the default trasport treatment.
  !
  ! See documentation of 'initHiPionInducedCollide' and
  ! 'initHiPionInducedCollideFull' for further information.
  !****************************************************************************

  !****************************************************************************
  !****g* initHiPion/DoOnlyOne
  ! SOURCE
  !
  logical, save :: DoOnlyOne = .true.
  !
  ! PURPOSE
  ! If the first interaction of beam and target particles is treated already
  ! here in the init (cf. DoPerturbativeInit), you may select whether a beam
  ! particle may interact only once (flag set to .true.) or with all other
  ! target nucleons (flag set to .false.).
  !
  ! See documentation of 'initHiPionInducedCollide' and
  ! 'initHiPionInducedCollideFull' for further information.
  !****************************************************************************

  real, save :: betaCM2Lab(1:3)

  !****************************************************************************
  !****g* initHiPion/betafak
  ! SOURCE
  !
  real,   save   :: betafak = 0.
  !
  ! PURPOSE
  ! This variable holds 2 times the value of the rapidity of the center of
  ! mass of the colliding pion-nucleon system (modulo sign).
  !   beta    = -(P(1,3)+P(2,3))/(P(1,4)+P(2,4))
  !   betafak = log((1+beta)/(1-beta))
  ! Here
  !   E_pion  = M_pion + Ekin_lab
  !   beta    = - sqrt(E_pion**2 - m_Pion**2) / (E_pion+M_nucleon)
  ! The rapidity of a particle is given by
  !   y_lab = 0.5 * log( (E+p_z)/(E-p_z) )   [lab frame]
  !   y_cm  = y_lab + 0.5 * betafak          [cm frame]
  !
  ! This is used in "HiPionAnalysis"
  !****************************************************************************

  !****************************************************************************
  !****g* initHiPion/totalPerweight
  ! SOURCE
  !
  real,save      :: totalPerweight=0.
  !
  ! PURPOSE
  ! This is the sum over the perweights of all initialized beam particles.
  ! Is needed in order to calculate total cross sections,
  !****************************************************************************

  !****************************************************************************
  !****g* initHiPion/ConfigParticles
  ! SOURCE
  !
  type(particle),dimension(2), save :: ConfigParticles
  !
  ! PURPOSE
  ! Store beam and target particles as dummy
  !****************************************************************************


  !****************************************************************************
  !****g* initHiPion/NucCharge
  ! SOURCE
  !
  integer, save :: NucCharge = -1
  !
  ! PURPOSE
  ! Select charge state of nucleons to scatter on. If this value is >=0,
  ! then we only scatter on nucleons with the respective charge,
  ! i.e. only on neutrons if NucCharge==0 and only on protons if NucCharge==1.
  ! Useful e.g. for selecting only pn events in a pd collision.
  !****************************************************************************


contains


  !****************************************************************************
  !****f* initHiPion/getBetaFak
  ! NAME
  ! real function getBetaFak()
  ! PURPOSE
  ! return the value of the variable "betafak"
  !****************************************************************************
  real function getBetaFak()
    getBetaFak = betafak
  end function getBetaFak

  !****************************************************************************
  !****f* initHiPion/getTotalPerweight
  ! NAME
  ! real function getTotalPerweight()
  ! PURPOSE
  ! This function returns the total perweight of the pions as they were
  ! initialized at the last call of initHiPionInduced.
  !****************************************************************************
  real function getTotalPerweight()
    getTotalPerweight = totalPerweight
  end function getTotalPerweight

!   !*************************************************************************
!   !****f* initHiPion/getNumberPions
!   ! NAME
!   ! integer function getNumberPions()
!   ! PURPOSE
!   ! This function returns the number of pions.
!   !*************************************************************************
!   integer function getNumberPions()
!     getNumberPions = numberPions
!   end function getNumberPions

  !****************************************************************************
  !****s* initHiPion/initHiPionInduced
  ! NAME
  ! subroutine initHiPionInduced(pertPart,realPart,targetNuc)
  ! PURPOSE
  ! This routine initializes high energy pions in pion-nucleus scattering.
  ! INPUTS
  ! * type(particle),dimension(:,:) :: pertPart,realPart -- particle vectors
  ! * type(tNucleus),pointer        :: targetNuc -- the target nucleus
  ! OUTPUT
  ! pertPart is modified
  !****************************************************************************
  subroutine initHiPionInduced(pertPart,realPart,targetNuc)

    use idTable, only: nucleon
    use nucleusDefinition
    use collisionNumbering, only: pert_Numbering,pert_firstNumbering12
    use hadronFormation, only: forceInitFormation
    use inputGeneral, only: fullEnsemble, localEnsemble, numEnsembles
    use energyCalc, only: energyDetermination
    use output, only: DoPr, timeMeasurement
    use CollHistory, only: CollHist_ClearArray
    use baryonPotentialModule, only: getNoPertPot_baryon
    use collisionTerm, only: readInputCollTerm => readInput


    type(particle),dimension(:,:),intent(inout) :: pertPart,realPart
    type(tNucleus),pointer                      :: targetNuc

    ! internal variables:
    integer  :: index            ! index of pion in vector pertPart
    integer  :: i,j,iEns,iPart   ! loop indizes
    logical, save :: initFlag=.true.
    logical :: lDummy

   if (DoPR(1)) then
     write(*,*)
     write(*,*) '**Initializing pions for HiPion induced events...'
     call timeMeasurement(.true.) ! Reset stopWatch
   end if

    if (initFlag) then
      call initInput
      call forceInitFormation
      call readInputCollTerm ! initialize enrgyCheck
      initFlag=.false.
    end if

    i = pert_firstnumbering12(.true.)

    call CollHist_ClearArray()

    totalPerweight=0.

    do j=1,numEnsembles ! Loop over all Ensembles
       index=1
       do i=1,nTestparticles
          do while(pertPart(j,index)%Id > 0) ! Find free place in particle vector
             index=index+1
             if (index.gt.size(pertPart,dim=2)) then
                write(*,*) 'Particle vector too small in initHiPion'
                write(*,*) 'Size=',size(pertPart,dim=2)
                write(*,*) 'Ensemble: ',j
                write(*,*) 'Number testparticles per ensemble=',nTestparticles
                stop
             end if
          end do

          call setKinematics(pertPart(j,index))
          call setPosition
          call setNumber(pertPart(j,index)) ! Give the particle its unique number

       end do
    end do

    if (.not.getNoPertPot_baryon()) then
       if (DoPr(2)) write(*,*) 'Updating energies of Particles'
       do iPart=1,size(realPart,dim=2)
          do iEns=1,numEnsembles
             if (realPart(iEns,iPart)%ID > 0) then
                call energyDetermination(realPart(iEns,iPart))
             end if
          end do
       end do
    end if

    if (DoPR(1)) then
      write(*,*)
      write(*,*) '**Initializing pions for HiPion induced events. finished!'
      write(*,*)
      call timeMeasurement()
    end if

    if (DoPerturbativeInit) then
       if (FullEnsemble) then
          lDummy = localEnsemble
          localEnsemble = .false.
          call initHiPionInducedCollideFull(pertPart,realPart)
          localEnsemble = lDummy
       else
          call initHiPionInducedCollide(pertPart,realPart)
       end if
    end if


!    stop

    !**************************************************************************

  contains

    !**************************************************************************
    !****s* initHiPionInduced/initInput
    ! NAME
    ! subroutine initInput
    ! PURPOSE
    ! Reads input out of jobcard. Namelist 'HiPionNucleus'.
    !**************************************************************************
    subroutine initInput
      use output, only: Write_ReadingInput, writeParticle
      use Dilepton_Analysis, only: Dilep_Init, Dilep_UpdateProjectile
      use constants, only: mN
      use callStack, only: traceback

      !************************************************************************
      !****n* initHiPion/HiPionNucleus
      ! PURPOSE
      ! Namelist for initHiPion includes:
      ! * distance
      ! * impact_parameter
      ! * ProjectileCharge
      ! * ProjectileID
      ! * ProjectileAnti
      ! * nTestparticles
      ! * ekin_lab
      ! * p_lab
      ! * DoPerturbativeInit
      ! * DoOnlyOne
      ! * minimumMomentum
      ! * useHermesPythiaPars
      ! * NucCharge
      !************************************************************************
      NAMELIST /HiPionNucleus/ distance, impact_parameter, nTestparticles, &
           ProjectileCharge, ProjectileID, ProjectileAnti, &
           ekin_lab, p_lab, &
           DoPerturbativeInit, DoOnlyOne, &
           minimumMomentum, useHermesPythiaPars, NucCharge

      real :: rdum

      call Write_ReadingInput('HiPionNucleus',0)
      rewind(5)
      read(5,nml=HiPionNucleus)

      if ((ekin_lab < 0.0).and.(p_lab < 0.0)) then
         call traceback("You either have to give ekin_lab or p_lab")
      else if ((ekin_lab > 0.0).and.(p_lab > 0.0)) then
         call traceback("You either have to give ekin_lab or p_lab, not both")
      end if

      write(*,*) '  Impact Parameter =',impact_parameter
      write(*,*) '  Distance         =',distance
      write(*,*) '  Kinetic Energy of projectile in lab frame =',ekin_lab
      write(*,*) '  Momentum of projectile in lab frame       =',p_lab
      write(*,*) '  Number of testparticles per ensemble: ',nTestparticles
      write(*,*) '  Projectile:    ID      = ',ProjectileID
      write(*,*) '                 charge  = ',ProjectileCharge
      write(*,*) '                 anti    = ',ProjectileAnti
      write(*,*) '  Do perturbative init            :',DoPerturbativeInit
      write(*,*) '  Do perturbative init - only One :',DoOnlyOne
      write(*,*) '  Minimum Momentum of init parts  :',minimumMomentum
      write(*,*) '  useHermesPythiaPars             :',useHermesPythiaPars
      write(*,*) '  NucCharge                       :',NucCharge

      call SetSwitchPythiaHermes(useHermesPythiaPars)

      !*** reset Distance:

      if (Distance.lt.0) then
         write(*,*) '... Distance<0 --> Distance=targetNuc%MaxDist-Distance'
         Distance=targetNuc%MaxDist-Distance
         write(*,*) '  Distance        =',distance
      end if


      !***Check wether pion is initialized inside the nucleus

      if (.not.doPerturbativeInit) then
         rdum = sqrt(distance**2+(max(0.0,impact_parameter))**2)

         if (rdum.lt.max(targetNuc%radius+targetNuc%surface,targetNuc%MaxDist)) then
            write(*,*)
            write(*,*) 'WARNING: Pions are initialized inside the nucleus! (1)'
            distance = max(targetNuc%radius+targetNuc%surface,targetNuc%MaxDist)
            write(*,*) 'distance reset to: ',distance
         end if
      end if

      !***Write Kinematics

      write(*,*)
      write(*,*) '--- Beam, Target: ---'

      call setKinematics(ConfigParticles(1))

      ConfigParticles(2)%ID = nucleon
      ConfigParticles(2)%charge = 1
      ConfigParticles(2)%mass = mN
      ConfigParticles(2)%momentum(0) = ConfigParticles(2)%mass

      call WriteParticle(6)
      call WriteParticle(6,0,1,ConfigParticles(1))
      call WriteParticle(6,0,2,ConfigParticles(2))

      write(*,*)
      write(*,*) '  Sqrt(S) = ',sqrtS(ConfigParticles)

!... generate a table ...
!      ekin_lab = 0.0
!      do
!         ekin_lab = ekin_lab + 0.1
!         call setKinematics
!         write(*,'(3f12.3)') ekin_lab,pertPart(1,1)%momentum(3),sqrtS(pertPart(1,1),X)
!         if (ekin_lab.gt.3.0) stop
!      end do


      call Write_ReadingInput('HiPionNucleus',1)

      call SetBetaFak

      call Dilep_Init(ekin_lab,nTestparticles,beta=betaCM2Lab)
      call Dilep_UpdateProjectile(ekin_lab,ConfigParticles(1)%momentum)

    end subroutine initInput


    !**************************************************************************
    !****s* initHiPionInduced/setKinematics
    ! NAME
    ! subroutine setKinematics
    ! PURPOSE
    ! Sets basic kinematics of the pions.
    !**************************************************************************
    subroutine setKinematics(part)

      use ParticleProperties, only: hadron, validCharge
      use callStack, only: traceback

      type(particle) :: part

      part%ID = ProjectileID
      part%charge = ProjectileCharge
      part%antiparticle = ProjectileAnti
      part%mass = hadron(ProjectileID)%mass

      if (.not.validCharge(part)) then
         call traceback("charge of projectile not valid!")
      end if

      part%perturbative=.true.
      part%productionTime=0.

      if (ekin_lab > 0.0) then
         part%momentum(0)=ekin_lab+part%mass
         part%momentum(1:3)=(/0.,0.,sqrt(part%momentum(0)**2-part%mass**2)/)
         p_lab = part%momentum(3)
      else
         part%momentum(0)=sqrt(p_lab**2+part%mass**2)
         part%momentum(1:3)=(/0.,0.,p_lab/)
         ekin_lab=part%momentum(0)-part%mass
      end if
      part%velocity(1:3)=part%momentum(1:3)/part%momentum(0)
      part%event(1:2)=pert_Numbering()

    end subroutine setKinematics

    !**************************************************************************
    !****s* initHiPionInduced/setPosition
    ! NAME
    ! subroutine setPosition
    !
    ! PURPOSE
    ! Sets positions of the pions.
    !
    ! If Impact_Parameter is choosen to be less than zero, than the impact
    ! parameter is choosen by a Monte-Carlo-decision. This is made such that
    ! the pion is initialized on a disk of  radius "bmax_Innerdisk" or on a
    ! ring which surrounds the inner disk and has an outer radius of
    ! "bmaxOuterRing".
    ! The probability to be on the inner disk is given by "pInnerDisk".
    ! The inner disk and the outer ring are separetely populated by a constant
    ! number density of pions.
    ! One distinguishes between inner disk and outer ring to have the possibility
    ! to have different  population densities. Assumed one would only have one
    ! disk, then most of the particles would be initialized with high
    ! impact-parameter where only few reactions take place.
    !**************************************************************************
    subroutine setPosition

      use random, only: rn
      use inputGeneral, only: fullEnsemble, numEnsembles
      use constants, only: pi

      real :: bmax_OuterRing              ! maximal Radius of outer ring
      real :: bmax_InnerDisk              ! Radius of inner ring
      real, parameter :: pInnerDisk=0.7   ! probability for initialization on inner ring
!      real :: minimalDistance=2.52        ! sqrt(maximal crossection of pion and nucleus/pi)
      real :: minimalDistance=1.25        ! sqrt(maximal crossection of pion and nucleus/pi)
      real :: phi
      real, parameter :: ratioRadius=1.8  ! bmax_Outerring=ratioRadius*nuclearRadius+...
      integer :: nTestparticlesTot
      real :: randomNumber, impact
      logical, save :: flag = .true.
      real :: radius, PerWeight1, PerWeight2

      nTestparticlesTot=nTestparticles*numEnsembles

      if (impact_parameter.ge.0.) then
         pertPart(j,index)%position=(/impact_Parameter,0.,-distance/)
         pertPart(j,index)%perweight=1./float(nTestparticlesTot)
      else     ! Monte Carlo decision to have impact parameter integration in the end
               ! maximum impact parameter of outer ring:
!         if(fullEnsemble) then
         if (fullEnsemble.and. (.not.localEnsemble)) then
            minimalDistance=minimalDistance/sqrt(float(numEnsembles))
         end if

         radius = max(targetNuc%radius, 1.0)

         bmax_OuterRing=ratioRadius*radius + minimaldistance
         bmax_InnerDisk=            radius + minimaldistance

         PerWeight1 = pi*bmax_InnerDisk**2/pInnerDisk &
                 &/float(nTestparticlesTot)*10  ! in mB (factor 10 due to fm**2 to mb conversion)
         PerWeight2 = pi*(bmax_OuterRing**2-bmax_InnerDisk**2)/(1.-pInnerDisk) &
                 &/float(nTestparticlesTot)*10  ! in mB (factor 10 due to fm**2 to mb conversion)

         if (flag) then
            write(*,*) '  Radius of outer ring  :',bmax_OuterRing
!            write(*,*) '            -----> sigma_max=',(bmax_OuterRing**2*31.4),'mb'
            write(*,*) '  Radius of inner circle:',bmax_InnerDisk
            write(*,*) '  perweight for pion in inner circle :', PerWeight1
            write(*,*) '  perweight for pion in outer ring   :', PerWeight2
            flag=.false.
         end if

         randomNumber=rn()
         phi=rn()*2*pi
         impact=rn()
         if (randomNumber.le.pInnerDisk) then ! impact parameter within nuclear radius

            radius = sqrt(impact)*bmax_InnerDisk
            pertPart(j,index)%position(1)=radius*cos(phi)
            pertPart(j,index)%position(2)=radius*sin(phi)
            pertPart(j,index)%position(3)=-distance
            pertPart(j,index)%perweight=PerWeight1

         else                                ! impact parameter not within nuclear radius

            radius = sqrt(impact*(bmax_OuterRing**2-bmax_InnerDisk**2)+bmax_InnerDisk**2)
            pertPart(j,index)%position(1)=radius*cos(phi)
            pertPart(j,index)%position(2)=radius*sin(phi)
            pertPart(j,index)%position(3)=-distance
            pertPart(j,index)%perweight=PerWeight2
         end if
      end if

      totalPerweight=totalPerweight+pertPart(j,index)%perweight

    end subroutine setPosition


    !**************************************************************************
    !****s* initHiPionInduced/setBetaFak
    ! NAME
    ! subroutine setBetaFak
    !
    ! PURPOSE
    ! set the variable "BetaFak"
    !**************************************************************************
    subroutine setBetaFak

      use constants, only: mN
      use ParticleProperties, only: hadron

      real :: E1,E2,m1,beta

      m1 = hadron(ProjectileID)%mass

      E1 = ekin_lab + m1
      E2 = mN

      beta = - sqrt(E1**2-m1**2)/(E1+E2)
      betaCM2Lab = beta * (/0.,0.,1./)
      betafak = log((1+beta)/(1-beta))

    end subroutine setBetaFak


  end subroutine initHiPionInduced


  !****************************************************************************
  !****s* initHiPion/initHiPionInducedCollide
  ! NAME
  ! subroutine initHiPionInducedCollide(pertPart,realPart)
  !
  ! PURPOSE
  ! if all (perturbative) pions and (real) nucleons are set,
  ! this routine performs for every pion its first collision.
  ! All collisions are done at time=zero.
  !
  ! The nucleons are sorted by its z-position. Then for a given pion
  ! it is checked whether a collision happens with the first nucleon, with
  ! the secons and so on.
  ! If the flag "DoOnlyOne" is true, only the first possible collision is
  ! performed [sigma ~ A^(2/3)], otherwise all possible collisions are
  ! considered [sigma ~ A].
  !
  ! Pions without any possible collisions are set somewhere in space with a
  ! large z-coordinate, otherwise the pion is deleted from the particle vector.
  !
  ! INPUTS
  ! * type(particle),dimension(:,:) :: pertPart,realPart -- particle vectors
  !
  ! OUTPUT
  ! pertPart is modified
  !****************************************************************************
  subroutine initHiPionInducedCollide(pertPart,realPart)

    use inputGeneral, only: numEnsembles
    use master_2Body, only: collide_2body
    use collisionTools, only: finalCheck
    use collisionNumbering, only: pert_numbering, pert_firstnumbering, reportEventNumber
    use output, only: DoPr, timeMeasurement, writeParticle
    use sorting, only: indexx
    use VolumeElements, only: VolumeElements_CLEAR, VolumeElements_SETUP_Real, VolumeElements_SETUP_Pert, VolumeElements_Statistics
    use pauliBlockingModule, only: checkPauli
    use CollHistory, only: CollHist_UpdateHist
    use insertion, only: setIntoVector, GarbageCollection

    type(particle),dimension(:,:),intent(inout),TARGET :: pertPart
    type(particle),dimension(:,:),intent(inout)        :: realPart

    integer :: iEnsemble,i,j,ii,HiEnergyType,number,nColl,nNucleon
    integer, parameter :: maxout = 100
    type(particle),dimension(1:2)      :: pair       ! incoming particles
    type(particle),dimension(1:maxout) :: finalState ! final state of particles
    integer, dimension(2,maxout) :: posOut
    logical :: flagOK, HiEnergyFlag, setFlag,NumbersAlreadySet
    real    :: time
    type(particle), POINTER :: pPart
    real,    dimension(size(realPart,dim=2)) :: zPos
    integer, dimension(size(realPart,dim=2)) :: zPosIndex
    real, dimension(0:20,2) :: SumEW
    real, dimension(2)      :: sSumEW
    real                    :: ssSumEW

    write(*,*)
    write(*,*) '**HiPionInduced: Doing Collisions, replacing pert particles!!!'
    write(*,*)

    call GarbageCollection(pertPart,.true.)
    call GarbageCollection(realPart)

    call VolumeElements_CLEAR
    call VolumeElements_SETUP_Pert(pertPart)
    call VolumeElements_SETUP_Real(realPart)
    call VolumeElements_Statistics

    SumEW = 0.
    time = 0.
    nNucleon = Size(realPart,dim=2)

    do iEnsemble=1,numEnsembles

       if (DoOnlyOne) then
          zPos(:) = realPart(iEnsemble,:)%position(3)
          call indexx(zPos,zPosIndex)
       end if


       pPartLoop: do i=1,nTestparticles
          pair(2) = pertPart(iEnsemble,i)

          nColl = 0

!!$          write(99,*) 'pertPart:'
!!$          call WriteParticle(99,iEnsemble,i,pair(2))


          rPartLoop: do j=1,nNucleon
             if (DoOnlyOne) then
                pair(1) = realPart(iEnsemble,zPosIndex(j))
             else
                pair(1) = realPart(iEnsemble,j)
             end if
             if (pair(1)%ID <= 0) cycle rPartLoop

             if (NucCharge >= 0) then
               if (pair(1)%charge /= NucCharge) cycle rPartLoop
             end if

!!$             write(99,*) 'realPart:'
!!$             call WriteParticle(99,iEnsemble,j,pair(1))

             call binSrts(sqrtS(pair),pair(2)%perweight)

             pair(2)%position(3) = pair(1)%position(3) ! force a collision

             call collide_2body (pair, finalState, time, flagOK, HiEnergyFlag, HiEnergyType)
             if (.not.flagOK) cycle rPartLoop

             flagOk = finalCheck (pair, finalState, HiEnergyFlag, 'initHiPion')
             if (.not. flagOk) cycle rPartLoop

             !...check Pauli Blocking
             flagOK = checkPauli(finalState, realPart)
             if ( DoPR(1) .and. (.not.flagOK) ) write(*,*) iEnsemble, 'Redo event: pauli blocked'
             if ( .not.flagOK ) cycle rPartLoop

             !...This was a successful event !

             NumbersAlreadySet = AcceptGuessedNumbers()

             nColl = nColl + 1

!             write(*,*) 'HiEnergyType=',HiEnergyType
!             write(99,*) 'HiEnergyType=',HiEnergyType

             ! Insert finalState into pertPart-vector (starting at nPionXX)

             number=pert_numbering(pair(1))
             finalState%event(1)=number
             finalState%event(2)=number

             finalState%lastCollisionTime=time ! = 0 !!! Attention!!!
             finalState%perweight = pair(2)%perweight
             finalState%firstEvent=pert_firstnumbering(pair(1),pair(2))
             finalState%perturbative=.true.

             posOut = 0
             call setIntoVector(finalState,pertPart(iEnsemble:iEnsemble,:),setFlag,NumbersAlreadySet,positions=posOut(:,:))
             posOut(1,:) = iEnsemble

             call CollHist_UpdateHist(pair, finalState, (/0,0/), posOut, finalState(1)%perweight)

             call ReportEventNumber (pair, finalState, finalState(1)%event, time, 212, HiEnergyType)

             !...Delete particles with small momenta
             do ii=1,maxout
                if (posOut(2,ii).lt.1) cycle
                pPart => pertPart(posOut(1,ii),posOut(2,ii))
                if (pPart%ID <= 0) cycle
                if (absMom(pPart) < minimumMomentum) then
                   if (DoPr(1)) write(*,*) '~~~ Deleting particle (momentum too small).'
                   if (DoPr(1)) call WriteParticle(6,posOut(1,ii),posOut(2,ii),pPart)
                   pPart%ID = 0
                end if

!                if (rapidity(pPart) < 5.5) then
!                   if (DoPR(1)) write(*,*) '~~~ Deleting particle (rapidity < 5.5)'
!                   pPart%ID = 0
!                end if

             end do


!199          continue

             pair%firstEvent = 0 ! they can rescatter as often as desired !!!

             if (DoOnlyOne) exit rPartLoop

          end do rPartLoop


          if (nColl.gt.0) then
             ! Erase Pion from pertPart-Vector:
             pertPart(iEnsemble,i)%ID = 0
          else
             ! Move Pion far out of the interaction zone:
             pertPart(iEnsemble,i)%position(3) = 200.0

             ! Erase Pion from pertPart-Vector:
!             pertPart(iEnsemble,i)%ID = 0
!             totalPerweight=totalPerweight-pertPart(iEnsemble,i)%perweight
!             write(*,*) 'Erasing non interacted pion!'
          end if


          SumEW(nColl,1) = SumEW(nColl,1) + 1
          SumEW(nColl,2) = SumEW(nColl,2) + pair(2)%perWeight


       end do pPartLoop

!       write(*,*) 'perturbative HiPion: iEnsemble = ',iEnsemble,' finished.'
    end do

    call binSrts()

    SumEW = SumEW / numEnsembles ! normalize to "per Ensemble"
    sSumEW = SUM(SumEW,dim=1)

    write(*,*)
    if (DoPR(1)) then
      write(*,*) 'sSumEW                   =',sSumEW
      write(*,*) 'sSumEW(1:2)-SumEW(0,1:2) =',sSumEW(1:2)-SumEW(0,1:2)
    end if
    write(*,*) '=====> sigmaTot   = ',(sSumEW(2)-SumEW(0,2))*numEnsembles

    if (nNucleon>1) then
      write(*,*) '=====> sigmaTot/A = ',(sSumEW(2)-SumEW(0,2))*numEnsembles/real(nNucleon)

      ssSumEW = 0.
      do i=1,20
        ssSumEW = ssSumEW + i * SumEW(i,2)
        if (sumEW(i,1)>0.) write(*,'(i3,3f12.5)') i,SumEW(i,1),SumEW(i,2),SumEW(i,2)/sSumEW(2)
      end do
      ssSumEW = ssSumEW / (sSumEW(2)-SumEW(0,2))
      write(*,*) '-->', ssSumEW
    end if

    if (DoPR(1)) then
      write(*,*)
      write(*,*) '**HiPionInduced: Doing Collisions, replacing pert particles!!! Finished.'
      write(*,*)
      call timeMeasurement()
    end if

!!$    do i=1,Size(pertPart,dim=2)
!!$       if (pertPart(1,i)%ID.gt.0) then
!!$          if  (pertPart(1,i)%FirstEvent.eq.0) then
!!$             write(*,*) 'xxx: ',pertPart(1,i)%number
!!$          endif
!!$       endif
!!$    enddo


  end subroutine initHiPionInducedCollide


  !****************************************************************************
  !****s* initHiPion/initHiPionInducedCollideFull
  ! NAME
  ! subroutine initHiPionInducedCollideFull(pertPart,realPart)
  !
  ! PURPOSE
  !
  ! same as initHiPionInducedCollide, but now for full (and) local
  ! ensemble calculations
  !
  ! if all (perturbative) pions and (real) nucleons are set,
  ! this routine performs for every pion its first collision.
  ! All collisions are done at time=zero.
  !
  ! The nucleons are sorted by its z-position. Then for a given pion
  ! it is checked whether a collision happens with the first nucleon, with
  ! the second and so on.
  ! If the flag "DoOnlyOne" is true, only the first possible collision is
  ! performed [sigma ~ A^(2/3)], otherwise all possible collisions are
  ! considered [sigma ~ A].
  !
  ! Pions without any possible collisions are set somewhere in space with a
  ! large z-coordinate, otherwise the pion is deleted from the particle vector.
  !
  ! INPUTS
  ! * type(particle),dimension(:,:) :: pertPart,realPart -- particle vectors
  !
  ! OUTPUT
  ! pertPart is modified
  !****************************************************************************
  subroutine initHiPionInducedCollideFull(pertPart,realPart)

    use inputGeneral, only: numEnsembles
    use master_2Body, only: collide_2body
    use collisionTools, only: finalCheck
    use collisionNumbering, only: pert_numbering, pert_firstnumbering, ReportEventNumber
    use output, only: DoPr, timeMeasurement, writeParticle
    use sorting, only: indexx
    use VolumeElements, only: VolumeElements_CLEAR, VolumeElements_SETUP_Real, VolumeElements_SETUP_Pert, VolumeElements_Statistics
    use pauliBlockingModule, only: checkPauli
    use CollHistory, only: CollHist_UpdateHist
    use insertion, only: setIntoVector, GarbageCollection

    type(particle),dimension(:,:),intent(inout),TARGET :: pertPart
    type(particle),dimension(:,:),intent(inout)        :: realPart

    integer :: iEnsemble,i, nNucleon,ii,HiEnergyType,number,nColl,iEnsR, jR
    integer :: j=0
    integer, parameter :: maxout = 100
    type(particle),dimension(1:2)      :: pair       ! incoming particles
    type(particle),dimension(1:maxout) :: finalState ! final state of particles
    integer, dimension(2,maxout) :: posOut
    logical :: flagOK, HiEnergyFlag, setFlag, NumbersAlreadySet
    real    :: time
    type(particle), POINTER :: pPart
    ! here we have size = numEnsembles*nParticles!!!
    real,    dimension(size(realPart)) :: zPos
    integer, dimension(size(realPart)) :: zPosIndex
    real, dimension(0:20,2) :: SumEW
    real, dimension(2)      :: sSumEW
    real                    :: ssSumEW

    write(*,*)
    write(*,*) '**HiPionInduced: Doing FULL Collisions, replacing pert particles!!!'
    write(*,*)

    call GarbageCollection(pertPart,.true.)
    call GarbageCollection(realPart)

    call VolumeElements_CLEAR
    call VolumeElements_SETUP_Pert(pertPart)
    call VolumeElements_SETUP_Real(realPart)
    call VolumeElements_Statistics

    SumEW = 0.
    nNucleon  = Size(realPart,dim=2)
    time = 0.

    if (DoOnlyOne) then
       zPos(:) = RESHAPE(realPart(:,:)%position(3), (/size(realPart)/) )
       call indexx(zPos,zPosIndex)
    else
       zPosIndex = (/ (j, j=1,size(zPosIndex)) /)
    end if

    do iEnsemble=1,numEnsembles

       pPartLoop: do i=1,nTestparticles
          pair(2) = pertPart(iEnsemble,i)

          nColl = 0

!!$          write(99,*) 'pertPart:'
!!$          call WriteParticle(99,iEnsemble,i,pair(2))

          rPartLoop: do j=1,Size(realPart)

             iEnsR = int((zPosIndex(j)-1)/nNucleon)+1
             jR    = mod((zPosIndex(j)-1),nNucleon)+1

             pair(1)  = realPart(iEnsR,jR)
             if (pair(1)%ID <= 0) cycle rPartLoop

             if (NucCharge >= 0) then
               if (pair(1)%charge /= NucCharge) cycle rPartLoop
             end if

!!$             write(99,*) 'realPart:'
!!$             call WriteParticle(99,iEnsemble,j,pair(1))

             call binSrts(sqrtS(pair),pair(2)%perweight)

             pair(2)%position(3) = pair(1)%position(3) ! force a collision

             call collide_2body (pair, finalState, time, flagOK, HiEnergyFlag, HiEnergyType)
             if (.not.flagOK) cycle rPartLoop

             flagOk = finalCheck (pair, finalState, HiEnergyFlag, 'initHiPion')
             if (.not. flagOk) cycle rPartLoop

             !...check Pauli Blocking
             flagOK = checkPauli(finalState, realPart)
             if (DoPR(1).and.(.not.flagOK)) write(*,*) iEnsemble, 'Redo event: pauli blocked'
             if (.not.flagOK) cycle rPartLoop

             !...This was a successful event !

             NumbersAlreadySet = AcceptGuessedNumbers()

             nColl = nColl + 1

!             write(*,*) 'HiEnergyType=',HiEnergyType
!             write(99,*) 'HiEnergyType=',HiEnergyType

             ! Insert finalState into pertPart-vector

             number=pert_numbering(pair(1))
             finalState%event(1)=number
             finalState%event(2)=number

             finalState%lastCollisionTime=time ! = 0 !!! Attention!!!
             finalState%perweight = pair(2)%perweight
             finalState%firstEvent=pert_firstnumbering(pair(1),pair(2))
             finalState%perturbative=.true.

             posOut = 0
             call setIntoVector(finalState,pertPart,setFlag,NumbersAlreadySet,positions=posOut(:,:))

             call CollHist_UpdateHist(pair, finalState, (/0,0/), posOut, finalState(1)%perweight)

             call ReportEventNumber (pair, finalState, finalState(1)%event, time, 212, HiEnergyType)


             !...Delete particles with small momenta
             do ii=1,maxout
                if (posOut(2,ii).lt.1) cycle
                pPart => pertPart(posOut(1,ii),posOut(2,ii))
                if (pPart%ID <= 0) cycle
                if (absMom(pPart) < minimumMomentum) then
                   if (DoPr(1)) write(*,*) '~~~ Deleting particle (momentum too small).'
                   if (DoPr(1)) call WriteParticle(6,posOut(1,ii),posOut(2,ii),pPart)
                   pPart%ID = 0
                end if

!                if (rapidity(pPart) < 5.5) then
!                   if (DoPR(1)) write(*,*) '~~~ Deleting particle (rapidity < 5.5)'
!                   pPart%ID = 0
!                end if

             end do

             pair%firstEvent = 0 ! they can rescatter as often as desired !!!

             if (DoOnlyOne) exit rPartLoop

          end do rPartLoop


          if (nColl.gt.0) then
             ! Erase Pion from pertPart-Vector:
             pertPart(iEnsemble,i)%ID = 0
          else
             ! Move Pion far out of the interaction zone:
             pertPart(iEnsemble,i)%position(3) = 200.0

             ! Erase Pion from pertPart-Vector:
!             pertPart(iEnsemble,i)%ID = 0
!             totalPerweight=totalPerweight-pertPart(iEnsemble,i)%perweight
!             write(*,*) 'Erasing non interacted pion!'
          end if


          SumEW(nColl,1) = SumEW(nColl,1) + 1
          SumEW(nColl,2) = SumEW(nColl,2) + pair(2)%perWeight


       end do pPartLoop

!       write(*,*) 'perturbative HiPion: iEnsemble = ',iEnsemble,' finished.'
    end do

    call binSrts()

    SumEW = SumEW / numEnsembles ! normalize to "per Ensemble"
    sSumEW = SUM(SumEW,dim=1)

    write(*,*)
    if (DoPR(1)) then
      write(*,*) 'sSumEW                   =',sSumEW
      write(*,*) 'sSumEW(1:2)-SumEW(0,1:2) =',sSumEW(1:2)-SumEW(0,1:2)
    end if
    write(*,*) '=====> sigmaTot   = ',(sSumEW(2)-SumEW(0,2))*numEnsembles

    if (nNucleon>1) then
      write(*,*) '=====> sigmaTot/A = ',(sSumEW(2)-SumEW(0,2))*numEnsembles/real(nNucleon)

      ssSumEW = 0.
      do i=1,20
        ssSumEW = ssSumEW + i * SumEW(i,2)
        if (sumEW(i,1)>0.) write(*,'(i3,3f12.5)') i,SumEW(i,1),SumEW(i,2),SumEW(i,2)/sSumEW(2)
      end do
      ssSumEW = ssSumEW / (sSumEW(2)-SumEW(0,2))
      write(*,*) '-->', ssSumEW
    end if

    if (DoPR(1)) then
      write(*,*)
      write(*,*) '**HiPionInduced: Doing FULL Collisions, replacing pert particles!!! Finished.'
      write(*,*)
      call timeMeasurement()
    end if

!!$    do i=1,Size(pertPart,dim=2)
!!$       if (pertPart(1,i)%ID.gt.0) then
!!$          if  (pertPart(1,i)%FirstEvent.eq.0) then
!!$             write(*,*) 'xxx: ',pertPart(1,i)%number
!!$          endif
!!$       endif
!!$    enddo


  end subroutine initHiPionInducedCollideFull


  !****************************************************************************
  !****s* initHiPion/binSrts
  ! NAME
  ! subroutine binSrts (srts, weight)
  !
  ! PURPOSE
  ! Create a histogram of the sqrt(s) of the first collision.
  ! If called without argument, the histogram is writen to
  ! "initHiPion_binSrts.dat"
  !
  ! INPUTS
  ! * real, optional: srts   -- the value to add
  ! * real, optional: weight -- weight of the value
  !****************************************************************************
  subroutine binSrts (srts, weight)
    use histf90

    real, intent(in), optional :: srts, weight

    logical, save :: binsrts_initflag=.true.
    type(histogram), save :: hSRTS

    ! parameters for the binning:
    ! (please adjust corresponding your needs)
    real, parameter :: xmin = 2.0, xmax = 5.5, xbin = 0.005

    if (binsrts_initflag) then
       call CreateHist(hSRTS,"srts of first coll", xmin, xmax, xbin)
       binsrts_initflag = .false.
    end if

    if (present(srts) .and. present(weight)) then
       call AddHist(hSRTS, srts, weight)
    else
       call WriteHist(hSRTS, file="initHiPion_binSrts.dat")
    end if

  end subroutine binSrts


end module initHiPion
