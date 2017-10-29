!******************************************************************************
!****m* /masterNBody
! NAME
! module masterNBody
! PURPOSE
! Implements N-body collisions.
!******************************************************************************
module masterNBody

  implicit none
  private

  public :: check_for_Nbody
  public :: generate_3body_collision
  public :: Nbody_analysis

  logical, parameter :: debug=.false.

contains

  !****************************************************************************
  !****s* masterNBody/check_for_Nbody
  ! NAME
  ! subroutine check_for_Nbody(index,time,coll_time,parts,n_found,ind_found,ind_min,sigma,gamma)
  ! PURPOSE
  ! searches for a third particle (or more particles) in proximity of the colliding pair
  ! INPUTS
  ! * integer, dimension(1:2), intent(in) :: index -- indices of colliding pair
  ! * real, intent(in) :: time -- current time in calc. frame
  ! * real, intent(in) :: coll_time -- collision time instant in calc. frame
  ! * type(particle), dimension(:) :: parts -- all particles from the ensemble of the colliding pair
  !
  ! OUTPUT
  ! * integer, intent(out) :: n_found -- number of particles found in vicinity of the colliding pair
  ! * integer, dimension(:) :: ind_found -- array of indices of the found particles
  ! * integer, intent(out) :: ind_min -- index of the closest particle to the collision point
  ! * real, intent(out) :: sigma, gamma -- total cross section (mbarn) and gamma-factor in CM frame of colliding particles
  !
  ! NOTES
  ! Only baryons are looked for in an interaction volume around colliding
  ! pair of baryons 1 and 2. The interaction volume is an ellipsoid centered in
  ! the center-of-mass of 1 and 2. The ellipsoid is symmetric with respect
  ! to the collision axis and compressed by gamma-factor along this axis.
  ! The transverse half-axis of the ellipsoid is sqrt(sigma/pi/10.),
  ! where sigma is the total cross section of 1 and 2 (mbarn).
  !****************************************************************************
  subroutine check_for_Nbody(index, time, coll_time, parts, n_found, ind_found, ind_min, sigma, gamma)

  use particleDefinition
  use constants, only: pi
!  use IdTable, only: isMeson, isBaryon
  use lorentzTrafo
  use XsectionRatios, only: getSigmaTotal, getSigmaScreened !,getRatio

  integer, dimension(1:2), intent(in) :: index
  real, intent(in) :: time
  real, intent(in) :: coll_time
  type(particle), dimension(:) :: parts
  integer, intent(out) :: n_found
  integer, dimension(:) :: ind_found
  real, intent(out)  :: sigma, gamma

  real :: rdist2,r02,dx2,dx2_min
  real, dimension(1:3) :: dx,dx_scaled,beta,beta_use
  integer :: nmax,i,ind_min
  type(particle) :: part1,part2
  real, dimension(0:3) :: dummy,coll_4point,particle_4point,particle_momentum
  real :: beta_use2,gamma_use,dx_times_beta,mstar1,mstar2,mstar_i

  integer, save :: imode=1
  ! 0 --- interaction volume is divided by the average gamma-factor of
  !       the colliding particles 1 and 2 in their c.m. frame
  ! 1 --- interaction volume is divided by the gamma-factor of
  !       the 3-d particle in the cm frame of 1 and 2

  logical, save :: debug=.false.

  n_found=0
  nmax=size(ind_found)

  !Copy colliding particles to working variables part1 and part2:
  part1=parts(index(1))
  part2=parts(index(2))

  if (debug) then
    write(*,*) 'In check_for_Nbody:', part1%Id,part2%Id
    write(*,*) 'time, coll_time:',time,coll_time
  end if

  !Cross section:
  !sigma = xsection(part1,part2) * getRatio((/part1,part2/))
  sigma = getSigmaTotal((/part1,part2/))

  ! In-medium screening:
  sigma = getSigmaScreened((/part1,part2/),sigma)

  if (sigma==0.) return

  r02=sigma/pi/10.

  !Evaluate positions of colliding particles at coll_time in the computatinal
  !frame:
  part1%position(1:3)=part1%position(1:3)&
                    &+(coll_time-time)*part1%velocity(1:3)
  part2%position(1:3)=part2%position(1:3)&
                    &+(coll_time-time)*part2%velocity(1:3)

  !Evaluate CM position at coll_time in the computational frame:
  coll_4point(1:3) = ( part1%momentum(0)*part1%position(1:3) &
                     &+part2%momentum(0)*part2%position(1:3) ) &
                  &/( part1%momentum(0) + part2%momentum(0) )

  !Evaluate the velocity of the CM-frame of colliding particles:
  beta(1:3) = ( part1%momentum(1:3)+part2%momentum(1:3) ) &
          & / ( part1%momentum(0) + part2%momentum(0) )

  !Lorentz transformations to the CM frame:

  ! 1) collision 4-point:
  coll_4point(0)=coll_time
  call lorentz(beta,coll_4point,'check_for_Nbody 1')

  ! 2) 4-momentum of the 1-st colliding particle:
  dummy=part1%momentum
  call lorentz(beta,dummy,'check_for_Nbody 2')
  part1%momentum=dummy

  ! 3) 4-momentum of the 2-nd colliding particle:
  dummy=part2%momentum
  call lorentz(beta,dummy,'check_for_Nbody 3')
  part2%momentum=dummy

  ! Dirac masses:
  mstar1 = part1%momentum(0)**2 - dot_product( part1%momentum(1:3),&
                                                 & part1%momentum(1:3) )
  mstar1 = sqrt( max(0.,mstar1) )
  mstar2 = part2%momentum(0)**2 - dot_product( part2%momentum(1:3),&
                                                 & part2%momentum(1:3) )
  mstar2 = sqrt( max(0.,mstar2) )

  if (imode.eq.0) then
    ! Use average gamma-factor in CM frame:
    gamma_use=0.5*(part1%momentum(0)/mstar1 &
           &     + part2%momentum(0)/mstar2)
    !*** Test:
    !   gamma_use=1.

    ! Use velocity of the 1-st colliding particle in CM frame:
    beta_use(1:3)=part1%momentum(1:3)/part1%momentum(0)
    beta_use2=dot_product(beta_use(1:3),beta_use(1:3))
    if (debug) write(*,*) ' beta_use2:', beta_use2
  end if

  gamma=0.
  dx2_min=r02
  ind_min=0

  particles_loop : do i=1,size(parts)

    if (parts(i)%ID < 0) exit particles_loop
    if (parts(i)%ID <= 0) cycle particles_loop

!   do not count mesons:
!    if (isMeson(parts(i)%ID)) cycle particles_loop

!   avoid autocorrelations:
    if (i==index(1) .or. i==index(2)) cycle particles_loop

!   Lorentz transformation of the i-th particle 4-position to the
!   CM frame of colliding particles:
    particle_4point=(/time,parts(i)%position(1:3)/)
    call lorentz(beta,particle_4point,'check_for_Nbody 4')

!   Lorentz transformation of the i-th particle 4-momentum to the
!   CM frame of colliding particles:
    particle_momentum(0:3)=parts(i)%momentum(0:3)
    call lorentz(beta,particle_momentum,'check_for_Nbody 5')

!   Compute position of the i-th particle at the collision time
!   in the CM frame of colliding particles:
    particle_4point(1:3)=particle_4point(1:3)&
                        &+(coll_4point(0)-particle_4point(0))&
                        &*particle_momentum(1:3)/particle_momentum(0)


    dx(1:3)=particle_4point(1:3)- coll_4point(1:3)

!   In the case of imode.eq.0 we determine whether the i-th particle is inside
!   an interaction volume, which is an ellipsoid centered in the collision point in
!   CM frame compressed along the collision axis by gamma-factor. The transverse
!   half-axis of the ellipsoid is r0=sqrt(sigma/pi/10.)

    if (imode.eq.1) then

!     Determine whether the i-th particle is inside an interaction volume,
!     which is an ellipsoid centered in the collision point in CM frame
!     compressed along the velocity of the i-th particle by gamma-factor of the i-th particle.
!     The transverse half-axis of the ellipsoid is r0=sqrt(sigma/pi/10.)

      mstar_i = particle_momentum(0)**2 - dot_product( particle_momentum(1:3),&
                                                     & particle_momentum(1:3) )
      mstar_i = sqrt( max(0.,mstar_i) )
      if (mstar_i.eq.0.) then
        write(*,*) 'In check_for_Nbody, mstar_i=0'
        stop
      end if
      gamma_use= particle_momentum(0)/mstar_i
      beta_use(1:3)=particle_momentum(1:3)/particle_momentum(0)
      beta_use2=dot_product(beta_use(1:3),beta_use(1:3))

    end if

    dx_times_beta=dot_product(dx(1:3),beta_use(1:3))
    dx_scaled=dx+(gamma_use-1.)*dx_times_beta/beta_use2*beta_use

    rdist2=dot_product(dx_scaled,dx_scaled)

    if (rdist2 < r02) then
      n_found=n_found+1
      if (n_found > nmax) then
        write(*,*) ' in check_for_Nbody: n_found > nmax', n_found,nmax
        write(*,*) ' Increase nmax !!!!'
        stop
      end if
      ind_found(n_found)=i
      dx2=dot_product(dx,dx)
      if (dx2.lt.dx2_min) then
        dx2_min=dx2
        ind_min=i
        gamma=gamma_use
      end if
    end if

  end do particles_loop

  end subroutine check_for_Nbody


!   real function xsection(part1,part2)
! !****s* masterNBody/xsection
! ! PURPOSE
! ! computes the total cross section (mb) for the scattering of part1 on part2
!
!   use particleDefinition
!   use particleProperties, only: isMeson, isBaryon, hadron
!   use twoBodyTools, only: sqrtS_free
!   use constants, only: mN
!
!   type(particle), intent(in) :: part1,part2 ! colliding particles
!   real :: stringfactor,srts,plab
!
! ! Currently only baryon-baryon collisions are included:
!   if (isMeson(part1%ID) .or. isMeson(part2%ID)) then
!     xsection= 0.
!     return
!   end if
!
! ! Simple recipie: *********
!   xsection= 40.
!   return
! !**************************
!
! ! Don't consider antibaryons:
!   if (part1%antiparticle .or. part2%antiparticle) then
!     xsection= 0.
!     return
!   end if
!
!   stringFactor=1.
!   if (part1%in_Formation) stringFactor=stringFactor*part1%scaleCS
!   if (part2%in_Formation) stringFactor=stringFactor*part2%scaleCS
!   if (stringFactor.lt.1E-8) then
!     xsection= 0.
!     return
!   end if
!
! ! Use proton-proton total cross section (mb):
!   srts=sqrtS_free((/part1,part2/))
!   if (srts.lt.2.6) then
!     if (hadron(part1%ID)%strangeness.ne.0 &
!        .or. hadron(part2%ID)%strangeness.ne.0 &
!        .or. hadron(part1%ID)%charm.ne.0 &
!        .or. hadron(part2%ID)%charm.ne.0 &
!        .or. part1%antiparticle .or. part2%antiparticle) then
!        xsection= 0.
!        return
!     end if
!   end if
!
!   plab= sqrt(srts**4/(4*mN**2)-srts**2)
!
!   if (plab < 0.4) then
!     xsection= 34.*(plab/0.4)**(-2.104)
!   else if (plab < 0.8) then
!     xsection= 23.5 + 1000.*(plab-0.7)**4
!   else if (plab < 1.5) then
!     xsection= 23.5 + 24.6/(1.+exp(-(plab-1.2)/0.10))
!   else if (plab < 3.) then
!     xsection= 41. + 60.*(plab-0.9)*exp(-1.2*plab)
!   else
!     xsection= 48.0+0.522*log(plab)**2-4.51*log(plab)
!   end if
!
! ! Cutoffs:
!   if (max(part1%ID,part2%ID) > 1) then
!     xsection= min(xsection,80.4) ! Resonance scattering
!   else
!     xsection= min(xsection,55.)  ! Nucleon-nucleon scattering
!   end if
!
!   xsection= xsection*stringfactor
!
! !  write(*,*) 'id1,id2,srts,plab,sigma:',part1%ID,part2%ID,srts,plab,xsection
!
!   end function xsection


  !****************************************************************************
  !****s* masterNBody/generate_3body_collision
  ! NAME
  ! subroutine generate_3body_collision(parts,finalState,collFlag,time,HiEnergyFlag,HiEnergyType)
  ! PURPOSE
  ! simulates a 3-body collision of baryons
  ! INPUTS
  ! * type(particle), dimension(1:3) ,intent(in) :: parts --- Incoming particles
  ! * real, intent(in) :: time ! current time
  ! OUTPUT
  ! * type(particle), dimension(:) :: finalState --- outgoing particles
  ! * logical, intent(out) :: collFlag --- .true. if collision happens
  ! * logical, intent(out) :: HiEnergyFlag  --- .true. if fritiof was used
  ! * integer, intent(out) :: HiEnergyType  --- 0:LowEnergy, 1:Fritiof, 2:Pythia
  ! NOTES
  ! Incoming particles 1 and 2 are actually colliding.
  ! The incoming particle 3 is the closest one to the c.m. of 1 and 2,
  ! which is determined before by the subroutine 'check_for_Nbody'.
  ! The energy and momentum of the particles 1,2 and 3 in their c.m.
  ! system are redistributed in such a way that the particle 3 stops in that
  ! system and gives its energy to the relative motion of 1 and 2.
  ! Then a 2-body collision of 1 and 2 is simulated in a usual way.
  !****************************************************************************
  subroutine generate_3body_collision(parts, finalState, collFlag, time, HiEnergyFlag, HiEnergyType)

  use particleDefinition, only: particle, setToDefault
  use lorentzTrafo, only: lorentzCalcBeta, lorentz
  use random, only: rnOmega
  use densitymodule, only: densityAt
  use mediumDefinition
  use dichtedefinition
  use preEventDefinition
  use master_2Body, only: XsectionMaster, XsectionCutOff, setKinematics, setKinematicsHiEnergy
  use twoBodyStatistics, only: sqrts_distribution
  use IdTable, only: isMeson, isBaryon, pion
  use propagation, only: updateVelocity, checkVelo
  use nBodyPhaseSpace, only: momenta_in_3BodyPS
  use RMF, only: getRMF_flag
  use energyCalc, only: energyCorrection
  use XsectionRatios, only: accept_event
  use mediumModule, only: mediumAt,getMediumCutOff
  use Annihilation, only: annihilate

  type(particle), dimension(1:3) ,intent(in) :: parts
  type(particle), dimension(:) :: finalState
  real, intent(in) :: time

  logical, intent(out) :: collFlag
  logical, intent(out) :: HiEnergyFlag
  integer, intent(out) :: HiEnergyType

  real, dimension(0:3) :: momentum_calc, momentum_calc_12, momentum_sum, momentum_final
  real, dimension(0:3) :: momLRF_12, impuls
  real, dimension(1:3) :: betacm,betacm_12, betaToLRF, position, betaDummy
  real :: sqrtsStar, mstar1, mstar2, mstar3, sqrts_12_ini, sqrts_12, sqrtS_corr, sqrtS_final
  real, dimension(0:7) :: sigs
  real :: tmp

  type(dichte) :: density
  type(medium) :: mediumAtColl
  type(particle) :: part1,part2,particle3
  type(preEvent),dimension(1:4) :: chosenEvent

  integer :: maxID, k, nAttempts, i

  logical :: successFlag, potFailure

  ! Copy incoming particles to working fields:
  part1=parts(1)
  part2=parts(2)
  particle3=parts(3)

!  if (debug) write(*,*) 'In generate_3body_collision, Ids, charges, masses:',&
!            &part1%Id,part2%Id,particle3%Id,&
!            &part1%charge,part2%charge,particle3%charge,&
!            &part1%mass,part2%mass,particle3%mass


  ! Total 4-momentum in the calculational frame:
  momentum_calc(0:3)=part1%momentum(0:3)+part2%momentum(0:3)+particle3%momentum(0:3)

!  if (debug) write(*,*) 'Total 4-momentum before collision:', momentum_calc

  ! Define the center-of-mass velocity:
  betacm = lorentzCalcBeta (momentum_calc, 'generate_3body_collision_1')

  ! Transform particle 4-momenta to their c.m. frame:
  impuls=part1%momentum
  call lorentz(betacm,impuls,'generate_3body_collision_1')
  part1%momentum=impuls
  impuls=part2%momentum
  call lorentz(betacm,impuls,'generate_3body_collision_2')
  part2%momentum=impuls
  impuls=particle3%momentum
  call lorentz(betacm,impuls,'generate_3body_collision_3')
  particle3%momentum=impuls

  ! Define total CM energy:
  sqrtsStar = part1%momentum(0) + part2%momentum(0) + particle3%momentum(0)

  ! Determine Dirac masses:
  mstar1 = part1%momentum(0)**2 - dot_product( part1%momentum(1:3),&
                                                 & part1%momentum(1:3) )
  mstar1 = sqrt( max(0.,mstar1) )
  mstar2 = part2%momentum(0)**2 - dot_product( part2%momentum(1:3),&
                                                 & part2%momentum(1:3) )
  mstar2 = sqrt( max(0.,mstar2) )
  mstar3 = particle3%momentum(0)**2 - dot_product( particle3%momentum(1:3),&
                                                 & particle3%momentum(1:3) )
  mstar3 = sqrt( max(0.,mstar3) )



  ! Reshuffle the 3-momenta of incoming particles:
  call threeBodyMomenta(1)

  ! Boost all particles back to the computational frame (needed for setKinematics/setKinematicsHiEnergy):
  betaDummy=-betacm
  impuls=part1%momentum
  call lorentz(betaDummy,impuls,'generate_3body_collision_4')
  part1%momentum=impuls
  impuls=part2%momentum
  call lorentz(betaDummy,impuls,'generate_3body_collision_5')
  part2%momentum=impuls
  impuls=particle3%momentum
  call lorentz(betaDummy,impuls,'generate_3body_collision_6')
  particle3%momentum=impuls


  ! Check energy-momentum conservation:
  momentum_sum(0:3)=part1%momentum(0:3)+part2%momentum(0:3)+particle3%momentum(0:3)
  if (abs(momentum_calc(0)-momentum_sum(0)).gt.1.e-06) then
    write(*,*) ' Energy is not conserved in generate_3body_collision:',momentum_calc(0),momentum_sum(0)
    stop
  else
    tmp=dot_product(momentum_calc(1:3)-momentum_sum(1:3),momentum_calc(1:3)-momentum_sum(1:3))
    if (tmp.gt.1.e-06) then
      write(*,*) ' Momentum is not conserved in generate_3body_collision:', tmp
      stop
    end if
  end if

  ! Total momentum of 1 and 2 in calc. frame after modifications:
  momentum_calc_12=part1%momentum+part2%momentum
  betacm_12 = lorentzCalcBeta (momentum_calc_12, 'generate_3body_collision_2')


  ! Check:
  !  sqrts_12=sqrts_12_ini
  !  momentum_calc_12(0:3)=parts(1)%momentum(0:3)+parts(2)%momentum(0:3)
  !  betacm_12 = lorentzCalcBeta (momentum_calc_12, 'generate_3body_collision_3')

  ! Read out medium at collision point in LRF:
  position=(part1%position+part2%position)/2.
  density=densityAt(position)

  mediumAtColl=mediumAt(density,position)

  ! Set total momentum of 1 and 2 in LRF. If density not negligible then boost is needed.
  momLRF_12=momentum_calc_12
  if ( density%baryon(0).gt.getMediumCutOff()/100. .and. .not.getRMF_flag() ) then
     betaToLRF = lorentzCalcBeta (density%baryon, 'generate_3body_collision_4')
     call lorentz(betaToLRF,momLRF_12,'generate_3body_collision_7')
  else
     betaToLRF=0.
  end if

  if (debug) write(*,*) 'sqrts_12:', sqrts_12

  ! Determine the total and elastic cross sections for the collision of 1 and 2
  ! and also the "preevent" (in case of low energy collision):

  if (.not.getRMF_flag()) then
      ! Here we assume that calculation is done without mean field at all
      sqrtS_corr = sqrts_12
  else
      sqrtS_corr = sqrts_12 - mstar1 - mstar2 + part1%mass + part2%mass
  end if


  call XsectionMaster(sqrtS_corr,(/part1,part2/),mediumAtColl,momLRF_12,&
       chosenEvent, sigs, HiEnergyFlag)

  if (debug) write(*,*) 'Xsections:', sigs(0:1), HiEnergyFlag

  ! We assume that the collision will always happen provided that the total
  ! cross section is not zero. The geometrical collision criterion is not
  ! checked here, since it has been already checked before in collide_2body.

  if (sigs(0).lt.XsectionCutOff) then
    collFlag=.false.
    return
  end if

  ! Initialize output:
  call setToDefault(finalState)

  finalState(1:4)%ID=chosenEvent(1:4)%ID
  finalState(1:4)%charge=chosenEvent(1:4)%charge
  finalState(1:4)%antiParticle=chosenEvent(1:4)%antiParticle
  finalState(1:4)%mass=chosenEvent(1:4)%mass

  finalState%productionTime = time
  finalState%formationTime  = time

  ! Set positions and kinematics of outgoing particles in collision of 1 and 2:

  if (.not.HiEnergyFlag) then

     if ( isBaryon(part1%Id) .and. isBaryon(part2%Id) .and. &
            & (part1%antiParticle.neqv.part2%antiParticle) .and. &
            & finalState(1)%Id.eq.pion .and. finalState(2)%Id.eq.pion .and.&
            & finalState(3)%Id.eq.pion ) then   ! Simulate annihilation

          if (part1%antiParticle) then
             call annihilate(part1,part2,time,finalState,collFlag,HiEnergyType)
          else
             call annihilate(part2,part1,time,finalState,collFlag,HiEnergyType)
          end if

          if (.not.collFlag) return

          call ResetPosition ! This also sets maxId

          HiEnergyFlag=.true.   ! Bad trick to avoid charge check after annihilation
          ! (unfortunately charge is not always conserved by annihilate)

     else

          call ResetPosition ! This also sets maxId

          HiEnergyType=0

          call setKinematics(sqrts_12,sqrtS_corr,betaToLRF,betacm_12,mediumAtColl,&
                             (/part1,part2/),finalState(1:maxID),collFlag)

          if (.not.collFlag) return

     end if



  else

     if ( .not.getRMF_flag() ) then

         call setToDefault(finalState)

         call setKinematicsHiEnergy(sqrts_12,sqrts_12, sigs, &
              betacm_12,(/part1,part2/),time,finalState,collFlag,&
              HiEnergyType)

         if (.not.collFlag) return

         call ResetPosition

     else

         nAttempts = 0

         do k = 1,10 ! Loop over attempts to simulate a Fritiof event with energy conservation

              nAttempts = nAttempts + 1

              call setToDefault(finalState)

              call setKinematicsHiEnergy(sqrts_12,sqrtS_corr, sigs, &
                   betacm_12,(/part1,part2/),time,finalState,collFlag,&
                   HiEnergyType)

              if (.not.collFlag) return

              call ResetPosition

              do i=1,maxId

                ! This is needed because energyCorrection takes finalState in the CM frame
                ! and returns it in the calculational frame:
                impuls = finalState(i)%momentum
                call lorentz(betacm_12,impuls, 'generate_3body_collision_8')
                finalState(i)%momentum = impuls

              end do

              call energyCorrection(sqrts_12, betaToLRF, betacm_12, mediumAtColl, &
                                   &finalState(1:maxId), successFlag, potFailure)

              if (potFailure) then
                write(*,*) ' In generate_3body_collision: potentialFailure happen'
                collFlag = .false.
                return
              end if

              if (successFlag) exit

         end do

         if ( .not.successFlag ) then

                write(*,*) ' In generate_3body_collision: energy correction FAILED:'
                write(*,*) ' Colliding particles:', part1%Id, part2%Id, &
                                                  & part1%antiparticle, part2%antiparticle
                write(*,*) ' Corrected sqrtS:', sqrtS_corr
                write(*,*) ' wished sqrtS:', sqrts_12
                write(*,*) ' Bare masses:', part1%mass, part2%mass
                write(*,*) ' Effective masses:', mstar1, mstar2
                write(*,*) ' Final particles:', finalState(1:maxId)%Id, finalState(1:maxId)%antiparticle
                write(*,*) ' Final bare masses:',  finalState(1:maxId)%mass
                write(*,*) ' Sum of final bare masses:', sum(finalState(1:maxId)%mass)
                do k = 0,3
                  momentum_final(k) = Sum(finalState(1:maxId)%momentum(k))
                end do
                sqrtS_final = sqrt(  momentum_final(0)**2 &
                               & - dot_product(momentum_final(1:3),momentum_final(1:3)) )
                write(*,*) ' final sqrtS:', sqrtS_final
                collFlag = .false.
                return

         end if

     end if

  end if

  ! Accept or reject the event randomly, using the in-medium cross section ratios:

  if ( .not.accept_event((/part1,part2/),finalState,sqrts_12_ini) ) then
     collFlag = .false.
     return
  end if


  ! Fill histograms for stat. analysis:
  call sqrts_distribution((/part1,part2/),3)

  ! Add 3-d colliding particle to the list of outgoing particles:
  ! (It is assumed that the particle 3 interacted elastically with
  !  1 and 2. Thus its productionTime and formationTime are not changed.)
  finalState(maxID+1)=particle3
  finalState(maxID+1)%number=0

  finalState(1:maxID+1)%lastCollisionTime=time

!  if (debug) write(*,*) 'Total 4-momentum of final particles:',&
!            &sum(finalState(:)%momentum(0)),sum(finalState(:)%momentum(1)),&
!            &sum(finalState(:)%momentum(2)),sum(finalState(:)%momentum(3))


  if ( .not.getRMF_flag() ) then
     call updateVelocity(finalState)
  else
     do i = 1,maxId+1
       finalstate(i)%velocity = finalState(i)%momentum(1:3)/finalState(i)%momentum(0)
       if (.not. checkVelo(finalState(i))) then
          write(*,*) ' In generate_3body_collision: checkVelo failed!'
          collFlag = .false.
          return
       end if
     end do
  end if


  if (debug) then
    do k=1,maxId+1
      if (finalState(k)%id <= 0) exit
      write(*,*) 'Id, antiparticle?, charge:',finalState(k)%Id,finalState(k)%antiparticle,&
                               &finalState(k)%charge
    end do
  end if

  contains

    !**************************************************************************
    !****s* generate_3body_collision/threeBodyMomenta
    ! NAME
    ! subroutine threeBodyMomenta(imode)
    ! PURPOSE
    ! Chooses randomly momenta of incoming three particles in their common CM
    ! frame.
    ! Also free energies of particles and the new sqrts_12 of colliding
    ! particles 1 and 2 are computed.
    !
    ! INPUTS
    ! * integer, intent(in) :: imode --- see below
    !
    ! possible values of imode:
    ! * 1: 3-body phase space sampling,
    ! * 2: 3-d particle is stopped, momenta of 1-st and 2-nd particles
    !   are randomly rotated,
    ! * 3: 3-d particle is stopped,  momenta of 1-st and 2-nd particles
    !   do not change their direction.
    !**************************************************************************
    subroutine threeBodyMomenta(imode)

    integer, intent(in) :: imode
    real, dimension(3,3) :: p3      ! first index : momentum, second: particle
    real, dimension(0:3) :: momentum_12,momentum_1
    real, dimension(1:3) :: beta_12
    real :: s_12,qcm_12
    real :: tmp!,phi,costh,sinth

    ! Total 4-momentum of 1 and 2 in the cm frame of 1,2 and 3 before modifications:
    momentum_12(0:3)=part1%momentum(0:3)+part2%momentum(0:3)

    ! Sqrt(s) of colliding pair before modifications:
    sqrts_12_ini=momentum_12(0)**2-dot_product(momentum_12(1:3),momentum_12(1:3))
    sqrts_12_ini=sqrt(max(0.,sqrts_12_ini))

    select case (imode)
    case (1)

      p3 = momenta_in_3BodyPS (sqrtsStar, (/mstar1,mstar2,mstar3/))

      part1%momentum(1:3)=p3(1:3,1)
      part2%momentum(1:3)=p3(1:3,2)
      particle3%momentum(1:3)=p3(1:3,3)

      part1%momentum(0) = sqrt( mstar1**2 + dot_product( part1%momentum(1:3), &
                                                           & part1%momentum(1:3) ) )
      part2%momentum(0) = sqrt( mstar2**2 + dot_product( part2%momentum(1:3), &
                                                           & part2%momentum(1:3) ) )
      particle3%momentum(0) = sqrt( mstar3**2 + dot_product( particle3%momentum(1:3), &
                                                           & particle3%momentum(1:3) ) )

      momentum_12(0:3)=part1%momentum(0:3)+part2%momentum(0:3)
      sqrts_12=momentum_12(0)**2-dot_product(momentum_12(1:3),momentum_12(1:3))
      sqrts_12=sqrt(max(0.,sqrts_12))

      return

    case (2,3)

      ! Maximum possible invariant energy of the 1-st and 2-nd particle
      ! assuming that the 3-d particle is stopped in the common c.m. frame:
      sqrts_12=sqrtsStar-mstar3

      if (debug) write(*,*) 'sqrts_12_ini,sqrts_12:', sqrts_12_ini,sqrts_12

      s_12=sqrts_12**2

      ! Compute maximum possible c.m. momentum of the 1-st and 2-nd particle:
      qcm_12=(s_12+mstar1**2-mstar2**2)**2/4./s_12-mstar1**2
      qcm_12=sqrt(max(0.,qcm_12))

    end select

    select case (imode)
    case (2)

      ! Random choice of the direction:
      momentum_1(1:3)=rnOmega()

    case (3)

      ! Cm velocity of 1 and 2 in the cm frame of 1,2 and 3:
      beta_12 = lorentzCalcBeta (momentum_12, 'generate_3body_collision_5')

      ! Lorentz trafo of 1 to the cm frame of 1 and 2:
      momentum_1(0:3)=part1%momentum(0:3)
      call lorentz(beta_12,momentum_1,'generate_3body_collision_8')

      ! Normalize momentum_1 to unity:
      tmp=Dot_Product(momentum_1(1:3),momentum_1(1:3))
      if (tmp.gt.0.) then
        momentum_1(1:3)=momentum_1(1:3)/sqrt(tmp)
      else
        ! Random choice of the direction:
        momentum_1(1:3)=rnOmega()
      end if

    end select

    ! Determine the new momenta of all three particles in their common CM frame,
    ! 3-d particle is assumed to be at rest (imode=2 or 3)
    part1%momentum(1:3)=qcm_12*momentum_1(1:3)
    part2%momentum(1:3)=-part1%momentum(1:3)
    particle3%momentum(1:3)=0.

    part1%momentum(0) = sqrt( mstar1**2 + dot_product( part1%momentum(1:3), &
                                                         & part1%momentum(1:3) ) )
    part2%momentum(0) = sqrt( mstar2**2 + dot_product( part2%momentum(1:3), &
                                                         & part2%momentum(1:3) ) )
    particle3%momentum(0) = sqrt( mstar3**2 + dot_product( particle3%momentum(1:3), &
                                                         & particle3%momentum(1:3) ) )

    end subroutine threeBodyMomenta


    subroutine ResetPosition
      ! Doku?????
      logical :: flag1,flag2
      maxId= ubound(finalState,dim=1)
      do k=lbound(finalState,dim=1),ubound(finalState,dim=1)
         if (finalState(k)%Id.eq.0) then
            maxID=k-1
            exit
         end if
      end do
      flag1= .false.
      flag2= .false.
      do k=lbound(finalState,dim=1),maxID
         if ( .not.flag1 .and. &
            ( (isBaryon(part1%ID) .and. isBaryon(finalState(k)%ID) .and. &
               (part1%antiparticle.eqv.finalState(k)%antiparticle))&
            .or. &
            (isMeson(part1%ID) .and. isMeson(finalState(k)%ID)) ) ) then
           finalState(k)%position = part1%position
           flag1= .true.
         else if ( .not.flag2 .and. &
            ( (isBaryon(part2%ID) .and. isBaryon(finalState(k)%ID) .and. &
               (part2%antiparticle.eqv.finalState(k)%antiparticle))&
            .or. &
            (isMeson(part2%ID) .and. isMeson(finalState(k)%ID)) ) ) then
           finalState(k)%position = part2%position
           flag2 = .true.
         else
           finalState(k)%position = position
         end if
      end do
    end subroutine ResetPosition



  end subroutine generate_3body_collision


  !****************************************************************************
  !****s* masterNBody/Nbody_analysis
  ! NAME
  ! subroutine Nbody_analysis(index,time,parts,n_found,ind_found,sigma,gamma,flag)
  ! PURPOSE
  ! Performs statistical analysis of many-body collisions using the output
  ! information from subroutine check_for_Nbody.
  ! INPUTS
  ! * integer, dimension(1:2) :: index --- indices of colliding pair
  ! * real :: time --- current time in calc. frame
  ! * type(particle), dimension(:) :: parts --- all particles from the ensemble of the colliding pair
  ! * integer :: n_found --- number of particles found in vicinity of the colliding pair
  ! * integer, dimension(:) :: ind_found --- array of indices of the found particles
  ! * real :: sigma --- total cross section (mbarn)
  ! * real :: gamma --- gamma-factor in CM frame of colliding particles
  ! * logical, optional :: flag --- if .true. --- do output
  !****************************************************************************
  subroutine Nbody_analysis(index, time, parts, n_found, ind_found, sigma, gamma, flag)

  use particleDefinition
  use inputGeneral, only: time_max,delta_T

  integer, dimension(1:2), intent(in) :: index ! indices of colliding pair
  real, intent(in) :: time  ! current time in calc. frame
  type(particle), dimension(:) :: parts ! all particles from the ensemble of the colliding pair
  integer, intent(in) :: n_found ! number of particles found in vicinity of the colliding pair
  integer, dimension(:) :: ind_found  ! array of indices of the found particles
  real, intent(in) :: sigma, gamma ! total cross section (mbarn) and gamma-factor in CM frame of colliding particles
  logical, optional, intent(in) :: flag ! if .true. --- do output

! Statistical variables (do not influence dynamics):
  integer, parameter :: Nmax=400 ! max possible value of n_found
  integer, save, dimension(0:1,2:Nmax+2) :: num_Nbody_collisions ! number of N-body collisions (0 -- at a given time step, 1 -- time integrated)
  integer, parameter :: numQ=100 ! number of Q-bins
  real, parameter :: dQ=0.1 ! Q-bin
  integer, save, dimension(numQ,2:Nmax+2) :: num_Nbody_collisions_vs_Q_pair, num_Nbody_collisions_vs_Q_cluster
  real, save, dimension(0:1,2:Nmax+2) :: Q_pair_vs_N, Q_cluster_vs_N

! Working variables:
  real, save :: Q_pair_tot, Q_cluster_tot
  real, dimension(0:3) :: p_pair,p_cluster
  real :: sqrts_pair,sqrts_cluster,Q_pair,Q_cluster,tmp
  integer :: i,ibin
  type(particle) :: part1,part2
  real, save :: sigma_average,gamma_average,srts_average
  real, save :: time_prev=0.
  logical, save :: switch=.true.

  if (switch) then
    open(36,file='check_for_Nbody.dat',status='unknown')
    write(36,*)'# Number of N-body collisions at a given time step vs time'
    write(36,*)'# time:  delta_T:  2-body:  3-body:  .... '
    open(37,file='Q_pair_vs_time.dat',status='unknown')
    write(37,*)'# Kin energy per baryon for colliding pair vs time'
    write(37,*)'# time:  2-body:  3-body:  .... '
    open(38,file='Q_cluster_vs_time.dat',status='unknown')
    write(38,*)'# Kin energy per baryon for all cluster vs time'
    write(38,*)'# time:  2-body:  3-body:  .... '
    open(39,file='sigma_gamma_srts_vs_time.dat',status='unknown')
    write(39,*)'# Cross section, gamma-factor and srts vs time'
    write(39,*)'# time:  sigma:  gamma:  srts:'
    num_Nbody_collisions(:,:)=0
    num_Nbody_collisions_vs_Q_pair(:,:)=0
    num_Nbody_collisions_vs_Q_cluster(:,:)=0
    Q_pair_vs_N(:,:)=0.
    Q_cluster_vs_N(:,:)=0.
    Q_pair_tot=0.
    Q_cluster_tot=0.
    switch=.false.
  end if

  if (time.ne.time_prev) then
    num_Nbody_collisions(0,:)=0
    Q_pair_vs_N(0,:)=0.
    Q_cluster_vs_N(0,:)=0.
    sigma_average=0.
    gamma_average=0.
    srts_average=0.
    time_prev=time
  end if

  if (present(flag)) then

    if (flag) then

      write(36,'(f7.3,1x,f7.3,100(1x,i5))') time,delta_T,&
                                         &num_Nbody_collisions(0,2:20)
      write(37,'(f7.3,100(1x,f6.2))')  time,&
               & Q_pair_vs_N(0,:)/float(max(1,num_Nbody_collisions(0,:)))
      write(38,'(f7.3,100(1x,f6.2))') time,&
               & Q_cluster_vs_N(0,:)/float(max(1,num_Nbody_collisions(0,:)))
      tmp=float(sum(num_Nbody_collisions(0,3:Nmax+2)))
      write(39,'(4(f7.3,1x))') time,&
               &sigma_average/max(1.,tmp),gamma_average/max(1.,tmp),&
               &srts_average/max(1.,tmp)
      num_Nbody_collisions(1,:)= num_Nbody_collisions(1,:) &
                            &  + num_Nbody_collisions(0,:)
      Q_pair_vs_N(1,:)= Q_pair_vs_N(1,:) + Q_pair_vs_N(0,:)
      Q_cluster_vs_N(1,:)= Q_cluster_vs_N(1,:) + Q_cluster_vs_N(0,:)

      if (time.ge.time_max) then

         open(36,file='check_for_Nbody_sum.dat',status='unknown')
         write(36,*)'# Integrated number of Nbody collisions:'
         write(36,*)'# N:    Number of collisions:'
         do i=2,ubound(num_Nbody_collisions,dim=2)
           write(36,*)i,num_Nbody_collisions(1,i)
         end do

         open(36,file='num_Nbody_collisions_vs_Q_pair.dat',status='unknown')
         write(36,*)'# Integrated number of Nbody collisions vs Q_pair:'
         write(36,*)'# Q_pair:    2-body:  3-body:  .... '
         do i=1,size(num_Nbody_collisions_vs_Q_pair,dim=1)
           Q_pair= (float(i)-0.5)*dQ
           write(36,'(f7.3,401i7)') Q_pair, num_Nbody_collisions_vs_Q_pair(i,:)
         end do

         open(36,file='Q_pair_vs_N.dat',status='unknown')
         write(36,*)'# Kin energy per particle vs cluster size'
         write(36,*)'# accouniting the colliding pair only'
         write(36,*)'# N:  Q_pair:'
         do i=2,ubound(Q_pair_vs_N,dim=2)
           write(36,*) i,Q_pair_vs_N(1,i)&
                      & /float(max(1,num_Nbody_collisions(1,i)))
         end do
         write(36,*)'# Average kin energy:', &
                  Q_pair_tot/float(max(1,sum(num_Nbody_collisions(1,:))))

         open(36,file='num_Nbody_collisions_vs_Q_cluster.dat',status='unknown')
         write(36,*)'# Integrated number of Nbody collisions vs Q_cluster:'
         write(36,*)'# Q_cluster:    2-body:  3-body:  .... '
         do i=1,size(num_Nbody_collisions_vs_Q_cluster,dim=1)
           Q_cluster= (float(i)-0.5)*dQ
           write(36,'(f7.3,401i7)') Q_cluster, num_Nbody_collisions_vs_Q_cluster(i,:)
         end do

         open(36,file='Q_cluster_vs_N.dat',status='unknown')
         write(36,*)'# Kin energy per particle vs cluster size'
         write(36,*)'# N:  Q_cluster:'
         do i=2,ubound(Q_cluster_vs_N,dim=2)
           write(36,*) i,&
           Q_cluster_vs_N(1,i)/float(max(1,num_Nbody_collisions(1,i)))
         end do
         write(36,*)'# Average kin energy:', &
                    Q_cluster_tot/float(max(1,sum(num_Nbody_collisions(1,:))))

      end if

      return

    end if

  end if

  if (n_found > Nmax) then
    write(*,*) 'In Nbody_analysis: n_found > Nmax', n_found,Nmax
    write(*,*) 'Increase Nmax !!!!'
    stop
  end if

  num_Nbody_collisions(0,2+n_found)=num_Nbody_collisions(0,2+n_found) + 1
  if (n_found.gt.0) then
    sigma_average=sigma_average+sigma
    gamma_average=gamma_average+gamma
  end if

! Copy colliding particles to working variables part1 and part2:
  part1=parts(index(1))
  part2=parts(index(2))

  p_pair= part1%momentum + part2%momentum
  p_cluster= p_pair
  if (n_found > 0) then
    do i=1,n_found
      p_cluster=p_cluster+parts(ind_found(i))%momentum
    end do
  end if
  sqrts_pair= sqrt(p_pair(0)**2 - dot_product(p_pair(1:3),p_pair(1:3)))
  if (n_found.gt.0)  srts_average=srts_average+sqrts_pair
  sqrts_cluster= sqrt(p_cluster(0)**2 &
               &    - dot_product(p_cluster(1:3),p_cluster(1:3)))
  Q_pair= sqrts_pair - part1%mass - part2%mass
  Q_cluster= sqrts_cluster - part1%mass - part2%mass
  if (n_found > 0) then
    do i=1,n_found
      Q_cluster=  Q_cluster - parts(ind_found(i))%mass
    end do
  end if
  Q_pair= Q_pair/2.
  Q_cluster= Q_cluster/float(2+n_found)
! Cluster distribution over kin energy per particle
! accounting for the colliding pair only:
  ibin= int(Q_pair/dQ) + 1
  if (ibin >= 1 .and. ibin <= numQ) then
    num_Nbody_collisions_vs_Q_pair(ibin,2+n_found)= &
                    & num_Nbody_collisions_vs_Q_pair(ibin,2+n_found) + 1
  end if
! Cluster distribution over kin energy per particle:
  ibin= int(Q_cluster/dQ) + 1
  if (ibin >= 1 .and. ibin <= numQ) then
    num_Nbody_collisions_vs_Q_cluster(ibin,2+n_found)= &
                    & num_Nbody_collisions_vs_Q_cluster(ibin,2+n_found) + 1
  end if
! Average kin energy per particle:
! accounting for the colliding pair only:
  Q_pair_vs_N(0,2+n_found)= Q_pair_vs_N(0,2+n_found) + Q_pair
  Q_pair_tot= Q_pair_tot + Q_pair
! Average kin energy per particle:
  Q_cluster_vs_N(0,2+n_found)= Q_cluster_vs_N(0,2+n_found) + Q_cluster
  Q_cluster_tot= Q_cluster_tot + Q_cluster

  end subroutine Nbody_analysis

end module masterNBody
