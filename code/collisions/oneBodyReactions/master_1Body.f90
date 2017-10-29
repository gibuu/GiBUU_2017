!******************************************************************************
!****m* /master_1Body
! NAME
! module master_1Body
! PURPOSE
! Implements all decays (a -> X).
!******************************************************************************
module master_1Body
  implicit none
  private


  !****************************************************************************
  !****g* master_1Body/correctEnergy
  ! SOURCE
  !
  logical,save :: correctEnergy=.true.
  ! PURPOSE
  ! Scale final state momenta to fulfill energy and momentum conservation.
  ! If .false. energy conservation is violated
  !****************************************************************************


  !****************************************************************************
  !****g* master_1Body/debug
  ! SOURCE
  !
  logical, parameter :: debug = .false.
  ! PURPOSE
  ! If .true., additional debugging information will be printed out.
  !****************************************************************************


  !****************************************************************************
  !****g* master_1Body/gammaCutOff
  ! SOURCE
  !
  real, parameter :: gammaCutOff = 1E-4
  ! PURPOSE
  ! If the decay width is lower than this value, than we treat this particle
  ! to be stable during propagation;
  ! particle is forced to decay in last time step, if gammaDecay > 0.
  !****************************************************************************


  !****************************************************************************
  !****g* master_1Body/StableInFormation
  ! SOURCE
  !
  logical, save :: StableInFormation = .true.
  ! PURPOSE
  ! Particles during its formation time are considered to be stable or not.
  !****************************************************************************


  !****************************************************************************
  !****g* master_1Body/omegaDecayMediumInfo
  ! SOURCE
  !
  logical, save :: omegaDecayMediumInfo = .false.
  ! PURPOSE
  ! Write out information about all decaying omega mesons to a file
  ! called "omegaMediumInfo.dat" (decay point, momentum, density, etc).
  !****************************************************************************


  !****************************************************************************
  !****g* master_1Body/omegaDecay_restriction
  ! SOURCE
  !
  integer, save :: omegaDecay_restriction = 0
  ! PURPOSE
  ! This switch, like omegaDecayMediumInfo, helps to analyze omega -> pi0 gamma decays.
  ! It will only have an effect for omegaDecayMediumInfo = .true.
  ! Possible values:
  ! * 0 = none (default)
  ! * 1 = vacuum ( rho < 0.1 rho0)
  ! * 2 = medium ( rho > 0.1 rho0)
  ! With the default value (0), all omega decays are carried out as usual.
  ! For the value 1, the decay products are only kept, if the decay happens in the vacuum (i.e. at rho < 0.1 * rho0).
  ! For the value 2, the decay products are only kept, if the decay happens in the medium (i.e. at rho > 0.1 * rho0).
  ! If the density does not meet these conditions, the decay products are simply removed and will not be put in the particle vector
  ! (and thus they will not appear in the analysis).
  !****************************************************************************


  logical, save :: initFlag=.true.


  public :: decayParticle, decayParticle_rhoN, assignCharge


contains


  !****************************************************************************
  !****s* master_1Body/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "master_1body".
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n* master_1Body/master_1body
    ! NAME
    ! NAMELIST /master_1body/
    ! PURPOSE
    ! Includes the switches:
    ! * correctEnergy
    ! * StableInFormation
    ! * omegaDecayMediumInfo
    ! * omegaDecay_restriction
    !**************************************************************************
    NAMELIST /master_1Body/ correctEnergy, StableInFormation, omegaDecayMediumInfo, omegaDecay_restriction

    integer :: ios

    call Write_ReadingInput('master_1Body',0)
    rewind(5)
    read(5,nml=master_1Body,IOSTAT=ios)
    call Write_ReadingInput('master_1Body',0,ios)

    write(*,*) 'correctEnergy          =', correctEnergy
    write(*,*) 'StableInFormation      =', StableInFormation
    write(*,*) 'omegaDecayMediumInfo   =', omegaDecayMediumInfo
    write(*,*) 'omegaDecay_restriction =', omegaDecay_restriction

    call Write_ReadingInput('master_1Body',1)
  end subroutine readInput


  !****************************************************************************
  !****s* master_1Body/decayParticle
  ! NAME
  ! subroutine decayParticle(resonanceIN, finalState, collisionFlag, pauliFlag, finalFlag, time, gammaOut)
  ! PURPOSE
  ! Evaluates for a given particle the final states of a decay process.
  ! First it checks wether the decay takes place in this time step,
  ! then it evaluates the final state if the decay criteria is fulfilled.
  !
  ! Treatment of the antiParticles:
  ! If the incoming resonance is an anti-particle, then we charge conjugate the
  ! incoming resonance and promote it to a particle. In the end we will
  ! do the charge conjugation again and promote all particle to antiparticles.
  ! If the particle is still in its formation period (%in_formation=.true.)
  ! then we do not allow it to decay!
  !
  ! INPUTS
  ! * type(particle),dimension(1:1) :: resonanceIN    -- incoming particle
  ! * logical                       :: finalFlag    -- if .true. then decay is forced to happen
  ! * real                          :: time         -- current time in fermi [does not influence any MC decisions; only used for analysis]
  !
  ! OUTPUT
  ! * type(particle),dimension(:)   :: finalState    -- produced final state
  ! * logical                       :: collisionFlag -- true if decay took place
  ! * logical                       :: pauliFlag     -- Set to .true. if Pauli blocking
  !   is already considered in decay decision
  ! * real, optional                :: gammaOut      -- full width [GeV]
  !
  !****************************************************************************
  subroutine decayParticle(resonanceIN, finalState, collisionFlag, pauliFlag, finalFlag, time, gammaOut)

    use IDTable, only: photon, rho, omegaMeson, isMeson, isBaryon, getAntiMeson
    use particleDefinition
    use particleProperties, only: hadron, nDecays, get_rho_dilep
    use MediumDefinition
    use decayChannels, only: Decay2BodyBaryon, Decay2BodyMeson, Decay3BodyMeson
    use inputGeneral, only: delta_T
    use baryonWidthMedium, only: decayWidthBaryonMedium
    use mesonWidthMedium, only: decayWidthMesonMedium
    use propagation, only: updateVelocity, checkVelo
    use constants, only: hbarc, rhoNull, mPi
    use random, only: rn
    use RMF, only: getRMF_flag
    use offShellPotential, only: getOffShellParameter
    use energyCalc, only: energyDetermination
    use history, only: setHistory
    use output, only: DoPR

    type(particle), intent(in)   :: resonanceIN
    type(particle), dimension(:) :: finalState
    logical, intent(out)         :: collisionFlag, pauliFlag
    logical, intent(in)          :: finalFlag
    real, intent(in)             :: time
    real, optional,intent(out)   :: gammaOut

    logical :: success, AntiParticleFlag
    type(medium) :: MediumAtPosition
    real :: gammaDecay,gamma,wahrscheinlichkeit, zufall
    real, dimension(1:3) :: betaToLRF
    real, dimension(0:3) :: momLRF
    integer :: i, anzahl, dId, antiCharge, antiID, stability, L
    type(particle)  :: resonance         !  used to copy the incoming resonance
    real, dimension(1:nDecays) :: width  ! field to hold the widths for the decay channels
    integer, save :: countspacelike=0

    collisionFlag=.false.
    pauliFlag=.false.

    resonance=resonanceIN ! Copy the input to temporary variable

    if (StableInFormation .and. resonance%in_Formation .and. .not.finalFlag) then
       ! Resonance is still in its formation period, therefore it can not decay
       return
    end if

    stability = hadron(resonance%ID)%stability
    if (finalFlag) then
       if (iand(stability,2) == 0) return
    else
       if (iand(stability,1) == 0) return
    end if

    ! (1) Initialize switches at first call
    if (initFlag) then
       call ReadInput
       initFlag=.false.
    end if

    ! (2) Evaluate medium informations
    if (.not. getRMF_flag()) then
       call getMomentum_and_Medium(resonance,momLRF, betaToLRF, mediumAtPosition)
    else
       momLRF = resonance%momentum
       betaToLRF = 0.
       mediumAtPosition%useMedium = .false.
    end if

    if (resonance%antiParticle) then
       antiParticleFlag=.true.
       ! Convert antiParticle to particle:
       resonance%antiparticle=.false.
       resonance%charge=-resonance%charge
    else
       antiParticleFlag=.false.
    end if

    ! (3) Evaluate the decay width of the particle
    if (isMeson(resonance%ID)) then
      width = decayWidthMesonMedium (resonance%ID, resonance%mass, resonance%charge, pauliFlag)
    else if (isBaryon(resonance%ID)) then
      call decayWidthBaryonMedium (resonance%ID,resonance%mass,momLRF,mediumATposition,width,pauliFlag)
    else
       write(*,*) 'Error in master_1Body/decayParticle'
       write(*,*) 'Resonance is no meson and no baryon:', resonance%ID
       write(*,*) 'Nonsense!!! Stop!!!'
       stop
    end if

    gammaDecay = Sum(width)
    if (present(gammaOut)) gammaOut=gammaDecay

    if (.not.finalFlag .and. gammaDecay<gammaCutOff) then
       collisionFlag=.false.
       return
    end if

    ! (4) Apply time criteria to decide on the decay of the particle
    if (dot_product(resonance%velocity(1:3),resonance%velocity(1:3)) > 1.) then
       countspacelike=countspacelike +1
       write(*,*) 'velocity greater than c in decayParticle! happened', countspacelike, ' times!'
       collisionFlag=.false.
       return
    end if
    !gamma=resonance%momentum(0)/resonance%mass   ! <-- not good with mean field
    gamma = 1./sqrt( 1. - dot_product(resonance%velocity(1:3),resonance%velocity(1:3)) )

    wahrscheinlichkeit=delta_T/gamma * gammaDecay/hbarc ! decay probability
    if (wahrscheinlichkeit > 1.0) then
       write(*,*) 'wahrsch > 1.0: ',wahrscheinlichkeit
    end if

    ! normally, for long timesteps, the decay prob is given by an exponential,
    ! since multiple events could happen. But in this code, we assume e.g. for
    ! 2body collisions etc, that the time step size is small, such that only
    ! one event is possible per timestep.
    ! You should *not* use the following expression:
!    wahrscheinlichkeit = 1.-exp(-wahrscheinlichkeit)


    zufall=rn()

    if (debug) then
       write(*,*) 'ID=', resonance%ID
       write(*,*) 'Gamma=', gamma
       write(*,*) 'DecayWidth=', gammaDecay
       write(*,*) 'Decay probability=', wahrscheinlichkeit
       write(*,*) 'Zufallszahl=', zufall
    end if

    if (finalFlag) then
       ! If "Finalflag" is set then we force the decay to happen,
       ! if there is some decay width available.
       if (gammaDecay > 0.) then
         zufall = 0.
       else if (resonance%ID==rho .and. resonance%charge==0 .and. get_rho_dilep() .and. resonance%mass<=2*mPi) then
         ! force rho0 to decay into dileptons (which we do not put in the particle vector)
         finalState(:)%ID = 0
         finalState(:)%charge = 0
         finalState(1)%momentum = resonance%momentum  ! trick to satisfy mom. cons. in finalCheck
         finalState(2)%momentum = 0.
         finalState(3)%momentum = 0.
         collisionFlag = .true.
         return
       else
         if (DoPR(3)) write(*,'(A,2I4,2ES13.5)') 'warning in DecayParticle: forced decay not possible', &
                                                 resonance%ID, resonance%charge, resonance%mass, gammaDecay   !, width
         !stop
       end if
    end if

    if (zufall>=wahrscheinlichkeit) then     !   monte-carlo decision
       ! Decay does not take place
       collisionFlag=.false.
       return
    end if

    ! Initialize output
    collisionFlag=.true.
    call setToDefault(finalstate)

    ! (5) Decay takes place: Determine the decay channel
    do
       zufall=rn()
       if (zufall>0.) exit
    end do
    ! Particles inherit perturbative switch: Important for the kinematics:
    finalState%perturbative=resonance%perturbative

    width=width/gammaDecay
    wahrscheinlichkeit=0.
    dId = 0
    do i=1,nDecays
       wahrscheinlichkeit=wahrscheinlichkeit+width(i)
       if (wahrscheinlichkeit>=zufall) then
          dId = hadron(resonance%ID)%decaysID(i)
          exit
       end if
    end do

    L = 0

    select case (dId)
    case (0)    ! no decay
       write(*,*) "error in decayParticle: no decay", resonance%ID, resonance%mass, gammaDecay, width
       stop
    case (1:)   ! 2Body Decay
       ! Set the position
       do i=1,2
          finalState(i)%position=resonance%position
       end do

       finalState(3:)%ID=0
       if (isBaryon(resonance%ID)) then
         finalState(1:2)%ID = Decay2BodyBaryon(dId)%ID(1:2)
         L = Decay2BodyBaryon(dId)%angularMomentum
       else
         finalState(1:2)%ID = Decay2BodyMeson(dId)%ID(1:2)
         L = Decay2BodyMeson(dId)%angularMomentum
       end if
       call assignCharge(finalState(1:2),resonance%ID,resonance%charge)
       anzahl = 2
    case (:-1)  ! 3Body Decay
       ! Set the position
       do i=1,3
          finalState(i)%position=resonance%position
       end do

       finalState(4:)%ID=0
       if (isBaryon(resonance%ID)) then
          write(*,*) 'Error in master_1Body : No three body decays for baryons implemented'
          write(*,*) dId, zufall, wahrscheinlichkeit
          write(*,*) 'critical error. Stop'
          stop
       else
          finalState(1:3)%ID=Decay3BodyMeson(-dId)%ID(1:3)
          finalState(1:3)%Charge=Decay3BodyMeson(-dId)%Charge(1:3)
       end if
       anzahl = 3
    end select

    ! Convert particles in final state to AntiParticles ,
    ! if there was an antiparticle in the incoming channel:
    if (AntiParticleFlag) then
       do i=1,anzahl
          if (isBaryon(finalState(i)%ID)) then
             finalState(i)%charge=-finalState(i)%charge
             finalState(i)%antiParticle=.true.
          else if (isMeson(finalState(i)%ID)) then
             call getAntiMeson(finalState(i)%ID,finalState(i)%charge,antiID,antiCharge)
             finalState(i)%ID=antiID
             finalState(i)%charge=antiCharge
          else
             write(*,*) 'Error in master_1Body. Final StateParticle is no Meson and no Baryon:', finalState%ID
             write(*,*) resonance%ID,resonanceIN%ID
             write(*,*) 'Stop'
             stop
          end if
       end do
    end if

    ! Set the kinematics of the final state
    ! * Must be done after AntiParticle conversion, since particles and antiparticles have different potentials!
    collisionFlag = setKinematics (resonance, finalState(1:anzahl), mediumAtPosition, betaToLRF, L)
    if (.not.collisionFlag) return

    ! Update velocities
    if (.not. getRMF_flag()) then
       do i = 1,anzahl
          if (finalstate(i)%id == 0) cycle
          call energyDetermination(finalstate(i),check=.true.)
          finalstate(i)%offshellparameter=getOffShellParameter(finalstate(i)%ID,  &
               finalstate(i)%Mass,finalstate(i)%momentum,finalstate(i)%position,success)
          if (.not. success) then
             collisionFlag = .false.
             return
          end if
          if (finalState(i)%ID /= photon) then
             call updateVelocity(finalState(i),success)
             if (.not. success) then
                write(*,*) 'Master_1Body(1): velocity ge 1: collisionFlag -> FALSE'
                collisionFlag = .false.
                return
             end if
          end if
       end do
    else
       do i = 1,anzahl
         finalstate(i)%velocity = finalState(i)%momentum(1:3)/finalState(i)%momentum(0)
         if (.not. checkVelo(finalState(i))) then
            write(*,*) 'Master_1Body(2): velocity ge 1: collisionFlag -> FALSE'
            collisionFlag = .false.
            return
         end if
       end do
    end if

    ! Label event by eventNumber
    finalState%event(1) = resonance%event(1)
    finalState%event(2) = resonance%event(2)

    finalState%perturbative = resonance%perturbative
    finalState%perweight    = resonance%perweight
    finalState%firstEvent   = resonance%firstEvent
    call setHistory (resonance, finalState)

    finalState%lastCollisionTime = time
    finalState%productionTime    = time
    finalState%formationTime     = time

    ! omega -> pi0 gamma
    if (omegaDecayMediumInfo .and. (resonance%ID == omegaMeson) .and. (dId == 6)) then
      call omegaMediumInfo (resonance, mediumAtPosition%density, time)
      if (omegaDecay_restriction==1 .and. (mediumATposition%density/rhoNull>0.1)) finalState%ID = 0
      if (omegaDecay_restriction==2 .and. (mediumATposition%density/rhoNull<=0.1)) finalState%ID = 0
    end if

!!$    write(*,*) '==='
!!$    call WriteParticle(6,0,1,resonanceIN(1))
!!$    if (abs(rapidity(resonanceIN(1))-3.5).lt.0.75) &
!!$               & write(*,*) 'MidRap:'
!!$    do i=1,size(finalState)
!!$       if (finalstate(i)%ID <= 0) cycle
!!$       call WriteParticle(6,1,i,finalstate(i))
!!$
!!$       if (finalstate(i)%ID==101 .and. finalstate(i)%charge==0) then
!!$          if (abs(rapidity(finalstate(i))-3.5).lt.0.75) &
!!$               & write(*,*) 'MidRap:',rapidity(finalstate(i)),&
!!$               & sqrt(finalstate(i)%momentum(1)**2+finalstate(i)%momentum(2)**2)
!!$       end if
!!$    end do

  end subroutine decayParticle


  !****************************************************************************
  !****f* master_1Body/decayParticle_rhoN
  ! NAME
  ! function decayParticle_rhoN (resonanceIN, finalState) result (gammaDecay)
  ! PURPOSE
  ! Performs the decay of a baryon resonance into a rho0-N final state and
  ! returns the decay width associated with that process.
  !
  ! INPUTS
  ! * type(particle) :: resonanceIN       -- incoming resonance
  !
  ! OUTPUT
  ! * type(particle) :: finalState(1:2)   -- produced final state
  ! * real           :: gammaDecay        -- decay width into rho0-N [GeV]
  !****************************************************************************
  function decayParticle_rhoN (resonanceIN, finalState) result (gammaDecay)

    use IDTable, only: isMeson, isBaryon, getAntiMeson, rho, nucleon
    use particleDefinition
    use particleProperties, only: hadron, nDecays
    use MediumDefinition
    use decayChannels, only: Decay2BodyBaryon
    use baryonWidthMedium, only: decayWidthBaryonMedium
    use random, only: rn
    use RMF, only: getRMF_flag
    use clebschGordan, only: CG
    use callstack, only: traceback
    use minkowski, only: abs4

    type(particle), intent(in) :: resonanceIN
    type(particle)             :: finalState(1:2)
    real                       :: gammaDecay

    logical :: success, antiParticleFlag
    type(medium) :: MediumAtPosition
    real :: wahrscheinlichkeit, zufall, betaToLRF(1:3), momLRF(0:3)
    integer :: i, dId, antiCharge, antiID
    type(particle) :: resonance          !  used to copy the incoming resonance
    real, dimension(1:nDecays) :: width  ! field to hold the partial widths for the different decay channels

    resonance = resonanceIN   ! Copy the input to temporary variable

    ! (1) Initialize switches at first call
    if (initFlag) then
      call ReadInput
      initFlag=.false.
    end if

    ! (2) Evaluate medium informations
    if (.not. getRMF_flag()) then
      call getMomentum_and_Medium (resonance, momLRF, betaToLRF, mediumAtPosition)
    else
      momLRF = resonance%momentum
      betaToLRF = 0.
      mediumAtPosition%useMedium = .false.
    end if

    if (resonance%antiParticle) then
      antiParticleFlag=.true.
      ! Convert antiParticle to particle:
      resonance%antiparticle=.false.
      resonance%charge=-resonance%charge
    else
      antiParticleFlag=.false.
    end if

    ! (3) Evaluate the decay width of the particle
    call decayWidthBaryonMedium (resonance%ID, resonance%mass, momLRF, mediumAtPosition, width, success)

    ! only keep rho0-N width
    do i=1,nDecays
      dId = hadron(resonance%ID)%decaysID(i)
      if (dId<9 .or. dId>12) then
        width(i) = 0.
      else
        ! apply isospin factor for rho0 decay
        width(i) = width(i) * CG(hadron(rho)%isospinTimes2,hadron(nucleon)%isospinTimes2, &
                                 hadron(resonance%Id)%isospinTimes2,0,resonance%charge*2-1)**2
      end if
    end do

    gammaDecay = Sum(width)

    ! Initialize output
    call setToDefault(finalstate)

    ! (5) Decay takes place: Determine the decay channel
    do
      zufall = rn()
      if (zufall>0.) exit
    end do

    ! Particles inherit perturbative switch (important for the kinematics)
    finalState%perturbative = resonance%perturbative

    width = width / gammaDecay
    wahrscheinlichkeit = 0.
    dId = 0
    do i=1,nDecays
      wahrscheinlichkeit = wahrscheinlichkeit + width(i)
      if (wahrscheinlichkeit >= zufall) then
        dId = hadron(resonance%ID)%decaysID(i)
        exit
      end if
    end do

    if (dID==0) then
      write(*,*) "problems in decayParticle_rhoN!"
      write(*,*) width(:)
      write(*,*) resonance%ID, resonance%charge, antiParticleFlag
      write(*,*) resonance%mass, abs4(resonance%momentum)
      call traceback()
    end if

    do i=1,2
      finalState(i)%position = resonance%position
    end do
    finalState(1:2)%ID     = Decay2BodyBaryon(dId)%ID(1:2)
    finalState(1:2)%charge = (/ 0, resonance%charge /)

    ! Convert particles in final state to AntiParticles,
    ! if there was an antiparticle in the incoming channel:
    if (antiParticleFlag) then
      do i=1,2
        if (isBaryon(finalState(i)%ID)) then
          finalState(i)%charge=-finalState(i)%charge
          finalState(i)%antiParticle=.true.
        else if (isMeson(finalState(i)%ID)) then
          call getAntiMeson(finalState(i)%ID,finalState(i)%charge,antiID,antiCharge)
          finalState(i)%ID=antiID
          finalState(i)%charge=antiCharge
        else
          write(*,*) 'Error in master_1Body. Final StateParticle is no Meson and no Baryon:', finalState%ID
          write(*,*) resonance%ID,resonanceIN%ID
          write(*,*) 'Stop'
          stop
        end if
      end do
    end if

    ! Set the kinematics of the final state
    ! * Must be done after AntiParticle conversion, since particles and antiparticles have different potentials!
    success = setKinematics (resonance, finalState(1:2), mediumAtPosition, betaToLRF, Decay2BodyBaryon(dId)%angularMomentum)

!     if (.not. success) then
!       print *, "error in decayParticle_rhoN: setKinematics failed!"
!       print *, resonance%ID, resonance%charge, resonance%mass
!       print *, finalState(1:2)%ID, finalState(1:2)%charge
!       stop
!     end if

  end function decayParticle_rhoN



  subroutine omegaMediumInfo (part, dens, time)
    use particleDefinition
    use constants, only: rhoNull
    use inputGeneral, only: current_run_number, num_Runs_sameEnergy, num_Energies
    use minkowski, only: abs4
    use histf90
    use PIL_omegaDec, only: PIL_omegaDec_Put

    type(particle), intent(in) :: part
    real, intent(in) :: dens
    real, intent(in) :: time
    logical, save :: init = .true.
    type(histogram), save :: dsigma_dm, decDensity
    real, parameter :: massres_sigma = 0.016
    !real, parameter :: massres_tau   = 0.023

    !**************************************************************************
    !****o* master_1Body/omegaMediumInfo.dat
    ! NAME
    ! file omegaMediumInfo.dat
    ! PURPOSE
    ! This file contains informations about omega mesons at decay time
    ! (event number, perweight, 4-momentum, position, bare mass, density at decy point, time, etc).
    !**************************************************************************
    open(66, file = "omegaMediumInfo.dat", position = 'append')

    if (init) then
      call CreateHist (dsigma_dm, 'dsigma/dm (pi0 gamma) in microbarn/GeV/A', 0., 1., 0.001)
      call CreateHist (decDensity,'density at decay point [rho/rho0]',0.,1.,0.01)

      rewind(66)
      write(66,'(A)')   "### This file contains the medium information for all omega -> pi0 gamma decays:"
      write(66,'(A,A)') "### run, event, perweight, 4-momentum [GeV], coordinates (x,y,z) [fm], ", &
                        "vacuum mass [GeV], density at decay vertex (rho/rho_0), time [fm]"
      write(66,*)
      init = .false.
    end if

    write(66,'(2I7,11ES15.7)') current_run_number, part%firstEvent, part%perweight, part%momentum, part%position, part%mass, &
                               dens/rhoNull, time
    close(66)

    call PIL_omegaDec_Put (part%firstEvent, dens)

    if (dens/rhoNull<0.1) then
      ! "vacuum" contribution
      call addHist (dsigma_dm, abs4(part%momentum), y = part%perweight/float(num_Runs_sameEnergy*num_Energies), y2 = 0.)
    else
      ! "in-medium" contribution
      call addHist (dsigma_dm, abs4(part%momentum), y = 0., y2 = part%perweight/float(num_Runs_sameEnergy*num_Energies))
    end if
    call addHist (decDensity, dens/rhoNull, part%perweight/float(num_Runs_sameEnergy*num_Energies))

    call WriteHist       (dsigma_dm, file="omegaMediumInfo_dsigma_dm.dat")
    call WriteHist_Gauss (dsigma_dm, "omegaMediumInfo_dsigma_dm_gauss.dat", massres_sigma)
    !call WriteHist_Novo  (dsigma_dm, "omegaMediumInfo_dsigma_dm_novo.dat",  massres_sigma, massres_tau)
    call WriteHist (decDensity, file='omegaDecayDensity.dat')

  end subroutine



  !****************************************************************************
  !****s* master_1Body/assignCharge
  ! NAME
  ! subroutine assignCharge(outPart,inID,inCharge)
  ! PURPOSE
  ! Random Choice of charges for a 2-body final-state in a given decay of a
  ! resonance "inID" of charge "inCharge".  The charges
  ! are choosen for hadronic decays according to the Clebsch-Gordan
  ! coefficients.
  ! For weak decays the charges are distributed according to their weak decay
  ! channels.
  ! INPUTS
  ! * type(particle), dimension(1:2) :: outPart  -- pair resulting of a decay
  ! * integer                        :: inID     -- Id of decaying particle
  ! * integer                        :: inCharge -- Charge of decaying particle
  ! OUTPUT
  ! * type(particle), dimension(1:2) :: outPart  -- pair resulting of a decay
  ! NOTES
  ! Be careful: variables named isospin can also contain isospin*2 to convert to
  ! integer values.
  !****************************************************************************
  subroutine assignCharge(outPart,inID,inCharge)
    use particleProperties, only: hadron
    use random, only: rn
    use particleDefinition
    use clebschGordan, only: clebschSquared
    use IdTable, only: photon, dsStar_plus, dsStar_minus, DMeson, dBar, dStar, dStarBar, isMeson, isBaryon
    use callstack, only: traceBack

    type(particle), dimension(1:2),intent(inOUT) :: outPart
    integer, intent(in) :: inID
    integer, intent(in) :: inCharge

    integer :: inIsospin_times2,inStrange, inCharm
    integer :: inIsospin_z_times2,outIsospin_zMax_times2, outIsospin_zMin_times2
    integer,dimension(1:2) :: outIsospin_times2
    integer,dimension(1:2) :: outStrange,outCharm, iz
    real :: xrn, prob
    real :: iso1,iso2,iso3,iso_z1,iso_z2  ! real Isospin (not multiplied by 2!!)
    integer :: j

    !**************************************************************************
    ! Take first care of weak decays, which violate isospin.
    ! Possible decay channels which are of weak nature:
    ! ds*^+ => X
    ! ds*^- => X
    ! d*^0 => X
    ! Y => photon+X (Includes ds*->photon ds)
    !**************************************************************************

    if (outPart(1)%Id.eq.photon) then
       outPart(1)%charge=0
       outPart(2)%charge=inCharge
       return
    end if
    if (outPart(2)%Id.eq.photon) then
       outPart(2)%charge=0
       outPart(1)%charge=inCharge
       return
    end if

    ! Test whether Ds*->Ds pion is taking place:

    if (inID==dSStar_minus .or. inID==dSStar_plus) then
       if (outPart(1)%ID == inID-2) then
          outPart(2)%charge=0
          outPart(1)%charge=inCharge
       else if (outPart(2)%ID == inID-2) then
          outPart(1)%charge=0
          outPart(2)%charge=inCharge
       else
          write(*,*) 'assignCharge:',inID,outPart(1)%ID,outPart(2)%ID
          call traceBack('Strange outgoing particle in dsStar decay')
       end if
       return
    end if

    ! Test whether D*(0)  is decaying (decays violating isospin rules):
    ! a) D*(0)-> D(0) Pi(0)
    ! b) D*(0)-> D(0) gamma (is already treated above with the photon case)
    ! NOTE :  The decays of the charged D* mesons [ D*(+-)-> D Pi] conserve Isospin    !

    if ((inID==dStar .or. inID==dStarBar).and.(Incharge==0)) then
       if (outPart(1)%ID == dMeson.or.outPart(1)%ID == dbar) then
          outPart(1:2)%charge = 0
       else if (outPart(2)%ID == dMeson.or.outPart(2)%ID == dbar) then
          outPart(1:2)%charge = 0
       else
          write(*,*) 'assignCharge:',inID,outPart(1)%ID,outPart(2)%ID
          call traceBack('Strange outgoing particle in dStar decay')
       end if
       return
    end if

    !**************************************************************************
    ! No weak decay, therefore use isospin arguments to distribute charges:
    !**************************************************************************


    ! Define isoSpin, strangeness and charm of the resonance
    inIsospin_times2=hadron(inID)%isospinTimes2
    inStrange=hadron(inID)%strangeness
    inCharm=hadron(inID)%charm

    ! Define isoSpin, strangeness and charm of the final state particles
    do j=1,2
       outIsospin_times2(j)=hadron(outPart(j)%ID)%isospinTimes2
       outStrange(j)=hadron(outPart(j)%ID)%strangeness
       outCharm(j)=hadron(outPart(j)%ID)%charm
    end do

    ! Checks
    if ((Sum(outStrange).ne.inStrange).or. (Sum(outCharm).ne.inCharm)) then
       write(*,*) 'assignCharge:',inID,outPart(1)%ID,outPart(2)%ID
       write(*,*) 'problems with strangeness or charm in master_1Body'
       write(*,*) inIsospin_times2,outIsospin_times2
       write(*,*) inStrange,outStrange
       write(*,*) inCharm,outCharm
       call traceBack()
    end if

    if (IsBaryon(inID)) then
       inIsospin_z_times2=2*inCharge-inStrange-inCharm-1
    else
       inIsospin_z_times2=2*inCharge-inStrange-inCharm
    end if

    ! Define bounds for z-component of isospin for the first decay product :
    outIsospin_zMax_times2=min(inIsospin_z_times2+outIsospin_times2(2),outIsospin_times2(1))
    outIsospin_zMin_times2=max(inIsospin_z_times2-outIsospin_times2(2),-outIsospin_times2(1))

    if ( outIsospin_zMax_times2.lt.outIsospin_zMin_times2) then
       write(*,*) 'assignCharge:',inID,outPart(1)%ID,outPart(2)%ID
       write(*,*) 'problems on master_1Body,   isospin _Out _zMax.lt.isospin _Out _zMin',  &
            & outIsospin_zMax_times2, outIsospin_zMin_times2
       write(*,*) 'outPart%ID: ',outPart%ID
       write(*,*) 'ID,IZ: ',inID, inCharge

       call traceBack()
    end if

    xrn=rn()
    iz(1)=outIsospin_zMin_times2
    prob=0.
    ! Convert to real isospin : divide by 2
    iso1=float(outIsospin_times2(1))/2.
    iso2=float(outIsospin_times2(2))/2.
    iso3=float(inIsospin_times2)/2.
    isospinLoop : do
       if (iz(1).gt.outIsospin_zmax_times2) then
          write(*,*) 'assignCharge:',inID,outPart(1)%ID,outPart(2)%ID
          write(*,*) 'master_1Body: problems in charge assignment'
          write(*,*) outPart%ID
          write(*,*) inID,inCharge,inStrange,inCharm
          write(*,*) iz,outIsospin_zmax_times2
          call traceBack()
       end if
       iz(2)=inIsospin_z_times2-iz(1)

       !Evaluate Clebsch-Gordon's :
       ! Convert to real isospin : divide by 2
       iso_z1=float(iz(1))/2.
       iso_z2=float(iz(2))/2.
       prob=prob+clebschSquared(iso1,iso2,iso3,iso_z1,iso_z2)

       if (prob.ge.xrn) then
          exit isoSpinLoop
       else
          iz(1)=iz(1)+2
       end if
    end do isospinLoop

    do j=1,2
       if (isMeson(outPart(j)%ID)) then
          outPart(j)%charge=int(float(iz(j)+outStrange(j)+outCharm(j))/2.)
       else
          outPart(j)%charge=int(float(iz(j)+outStrange(j)+outCharm(j)+1)/2.)
       end if
    end do

  end subroutine assignCharge

  !****************************************************************************
  !****s* master_1Body/getMomentum_and_Medium
  ! NAME
  ! subroutine getMomentum_and_Medium(resonance,momentumLRF, betatToLRF, mediumAtDecay)
  ! PURPOSE
  ! Evaluates for a given particle the medium at its position and transforms its
  ! momentum into the  LRF.
  ! INPUTS
  ! * type(particle)        :: resonance      -- considerer resonance particle
  ! OUTPUT
  ! * type(medium)          :: mediumAtDecay
  ! * real, dimension(0:3)  :: momentumLRF    -- momentum in the LRF
  ! * real, dimension(1:3)  :: betaToLRF      -- beta for boost to LRF
  !****************************************************************************
  subroutine getMomentum_and_Medium(resonance,momLRF, betaToLRF, mediumAtDecay)

    use densitymodule, only: densityAt
    use particleDefinition
    use mediumDefinition
    use dichtedefinition
    use lorentzTrafo, only: lorentz
    use mediumModule, only: mediumAt,getMediumCutOff

    type(particle), intent(in) :: resonance
    type(medium),intent(out)            :: mediumAtDecay
    real,dimension(0:3) ,intent(out)  :: momLRF
    real, dimension(1:3),intent(out)   :: betaToLRF

    type(dichte) :: density
    real,dimension(1:3)  :: position

    ! (1) Read out medium at collision point in LRF
    position=resonance%position
    density=densityAt(position)
    mediumAtDecay=mediumAt(density,position)

    ! (2) Define total momentum in LRF. If density not negligible then boost is needed.
    momLRF= resonance%momentum(0:3)
    if (density%baryon(0).gt.getMediumCutOff()/100.) then
       betaToLRF=density%baryon(1:3)/density%baryon(0)
       call lorentz(betaToLRF, momLRF, 'master_1Body(1)')
    else
       betaToLRF=0.
    end if

  end subroutine getMomentum_and_Medium


  !****************************************************************************
  !****f* master_1Body/setKinematics
  ! NAME
  ! function setKinematics (resonance, finalState, mediumAtCollision, betaToLRF, L) result(collisionFlag)
  ! PURPOSE
  ! Evaluates the kinematics for the "finalState" particles.
  ! INPUTS
  ! * type(particle)              :: resonance         -- resonance which decays into "FinalState"
  ! * type(medium)                :: mediumAtCollision -- medium information
  ! * real,dimension(1:3)         :: betaToLRF         -- beta for boost to Local Rest Frame
  ! * integer                     :: L                 -- angular momentum in final state
  ! * type(particle),dimension(:) :: finalState        -- vector of final-state particles
  ! OUTPUT
  ! * type(particle),dimension(:) :: finalState        -- vector of final-state particles
  ! * logical                     :: collisionFlag     -- "true" if kinematics could be set, "false" if not
  ! NOTES
  ! It's important that the IDs and charges of the "finalState"
  ! particles are already set when calling this subroutine.
  !
  ! Only kinematics including masses of this finalState will be set.
  !****************************************************************************
  function setKinematics (resonance, finalState, mediumAtColl, betaToLRF, L) result(collisionFlag)

    use mediumDefinition
    use particleDefinition
    use finalStateModule, only: assMass, massAss
    use energyCalc, only: energyCorrection
    use lorentzTrafo, only: lorentz, lorentzCalcBeta
    use IdTable, only: isMeson
    use offShellPotential, only: treatParticleOffShell
    use minkowski, only: abs4
    use callstack, only: traceback

    type(particle), intent(in)       :: resonance
    type(medium), intent(in)         :: mediumAtColl
    real, dimension(1:3), intent(in) :: betaToLRF
    integer, intent(in)              :: L
    type(particle), dimension(:)     :: finalState
    logical                          :: collisionFlag

    integer :: i,j
    real, parameter, dimension (1:3) :: spotOut3 = 0.
    logical :: flag, successFlag
    integer, parameter :: maxCorrectLoop=10     ! maximal number of iterations for energy correction
    real :: srts, srts_vacuum, betaToCM(1:3)
    type(particle), dimension(1:2) :: pair ! Incoming particles

    ! Usually massass and assmass expect a pair of initial state particles, here there is a single initial state particle. Therefore:
    pair(1)=resonance
    pair(2)%ID=0

    collisionFlag=.true.

    ! CHECK INPUT
    do i=lbound(finalState,dim=1),ubound(finalState,dim=1)
      if (finalState(i)%Id <= 0) then
        write(*,*) 'SetKinematics[1]: Particle ID == 0!!!'
        write(*,*) resonance%ID, finalState(:)%ID
        call traceback()
      end if
    end do

    ! Set boost velocities:

    flag=.true.

    ! Set Energy and momenta
    srts=abs4(resonance%momentum)
    srts_Vacuum=resonance%mass

    betaToCM = lorentzCalcBeta (resonance%momentum, 'master_1Body/setKinematics')

    select case (size(finalState,dim=1))
    case (2)
       energyCorrectLoop : do i=1, maxCorrectLoop
          ! two body final state.
          if (debug) write(*,*) 'Vor massAss'
          call massass(srts_vacuum,mediumAtColl,pair,finalState,betaToLRF,betaToCM,L,flag)
          if (.not. flag) then
             write(*,*) 'SetKinematics[1]: Impossible to find final state [2]'
             write(*,*) resonance%Id, finalState%ID, srts
             !stop
             successFlag=.false.
             cycle
          end if
          if (correctEnergy) then
             if (debug) write(*,*) '****Correcting Energy*********************'
             if (debug) write(*,*) 'wished srts=', srts
             if (debug) write(*,*) 'initial srts=', sqrts(Pair(1),Pair(2))
             call energyCorrection(srts,betaToLRF,betaToCM, mediumAtColl, finalState, successFlag)
             if (debug) write(*,*) 'final srts=', sqrts(finalState(1),finalState(2))
!            If (debug) write (*,*) 'initial srts=', sqrts(Pair(1),Pair(2))
             if (successFlag) exit energyCorrectLoop
          else
             ! just boost to Calculation frame
             successFlag=.true.
             do j=1,2
                call lorentz(-betaToCM,finalState(j)%momentum(0:3), 'master_1Body(2)')
             end do
             exit energyCorrectLoop
          end if
       end do energyCorrectLoop

       if (.not.successFlag) then
          write(*,'(A,I4,A,2I4,A,2ES12.4)') 'SetKinematics[1]: Energy correction failed. [2]', &
                                            resonance%ID,' ->',finalState(1:2)%Id,' @',resonance%mass,srts
          ! Kill Event
          collisionFlag=.false.
          finalState(1:2)%ID=0
          return
       end if

    case (3)
       do i=1, maxCorrectLoop
          ! three body final state. Neglecting possible scalar potentials, since spotOut=0
          call assMass(srts_vacuum,mediumAtColl,pair,finalState,spotOut3,betaToLRF,betaToCM,flag)
          if (.not. flag) then
             write(*,*) 'SetKinematics[1]: Impossible to find final state [3]'
             write(*,*) resonance%Id, finalState%ID, srts
             if (isMeson(resonance%ID) .and. treatParticleOffShell(resonance%ID,resonance%offShellParameter)) then
               ! Kill Event
               collisionFlag=.false.
               finalState(1:3)%ID = 0
               return
             else
               stop
             end if
          end if
          if (correctEnergy) then
             call energyCorrection(srts,betaToLRF,betaToCM, mediumAtColl, finalState, successFlag)
             if (successFlag) exit
          else
             ! just boost to Calculation frame
             successFlag=.true.
             do j=1,3
                call lorentz(-betaToCM,finalState(j)%momentum(0:3), 'master_1Body(3)')
             end do
          end if
       end do

       if (.not.successFlag) then
          write(*,'(A,I4,A,3I4,A,2ES12.4)') 'SetKinematics[1]: Energy correction failed. [3]', &
                                            resonance%ID,' ->',finalState(1:3)%Id,' @',resonance%mass,srts
          collisionFlag=.false.
          ! Kill Event
          finalState(1:3)%ID = 0
          return
       end if

    case default
       write(*,*) 'Error in setKinematics: No treatment of more than a three-particle-final state implemented yet!.' &
         & ,'Also single-particle Final-state not implemented'
       write(*,*) size(finalState), finalState(:)%ID
    end select

    ! set velocities in vacuum approximation

!    Do i=lBound(finalState,dim=1), uBound(finalState, dim=1)
!       finalstate(i)%velocity=finalState(i)%momentum(1:3)/finalState(i)%momentum(0)
!    End do
    finalState%scaleCS=1.
    finalState%in_formation=.false.

  end function setKinematics


end module master_1Body
