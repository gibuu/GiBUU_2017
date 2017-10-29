!******************************************************************************
!****m* /potentialModule
! NAME
! module potentialModule
! PURPOSE
! Does administration for the hadronic potentials in the code.
!******************************************************************************
module potentialModule

  implicit none
  private

  !****************************************************************************
  !****s* potentialModule/massDetermination
  ! NAME
  ! subroutine massDetermination(particleID,momentumLRF,med,baremass,verbose,success,pos)
  !
  ! subroutine massDetermination(teilchen,betaToCF,verbose,success)
  !
  ! PURPOSE
  ! This subroutine determines the bare mass of a particle with a given
  ! 4-momentum considering a momentum-dependent scalar potential that is
  ! defined in the local rest frame.
  !
  ! NOTES
  ! We keep p(0:3) constant!
  ! To solve   p_0=sqrt(p**2+m**2)+V_LRF_0 (p(1:3)) in the LRF frame => solution for m!
  !
  ! Attention: Does not respect Coulomb potential !!!
  !
  ! INPUTS
  ! * integer              :: ID             -- ID of particle to consider
  ! * real, dimension(0:3) :: momentumLRF    -- momentum of particle in LRF
  ! * type(medium)         :: med            -- medium information
  ! * logical, OPTIONAL    :: verbose        -- flag to print warnings (Default: .true.)
  !
  ! or:
  ! * type(particle)       :: teilchen -- Particle whose energy should be calculated.
  ! * real,dimension(3), OPTIONAL :: betaToCF -- Velocity of "Calc Frame" in the
  !   frame where the particle's momentum is defined
  ! * logical, OPTIONAL    :: verbose        -- flag to print warnings (Default: .true.)
  !
  ! OUTPUT
  ! * real                 :: baremass
  ! * logical, OPTIONAL    :: success
  !
  ! or:
  ! * type(particle)       :: teilchen -- Particle whose energy should be calculated.
  ! * logical, OPTIONAL    :: success
  !
  !****************************************************************************
  Interface massDetermination
     Module Procedure massDetermination_part,massDetermination_noPart
  End Interface

  !****************************************************************************
  !****f* potentialModule/potential_LRF
  ! NAME
  ! real function potential_LRF(teilchen, addCoulomb)
  !
  ! real function potential_LRF(ID,IQ,mom,pos, addCoulomb)
  !
  ! real function potential_LRF(teilchen,density, addCoulomb)
  !
  ! real function potential_LRF(teilchen,rhoB, addCoulomb)
  !
  ! real function potential_LRF(teilchen,med, addCoulomb)
  !
  ! PURPOSE
  ! Evaluates the hadronic potential of particles in the local rest frame (LRF).
  ! It's considered to be the 0th component of a vector potential.
  ! If 'addCoulomb' is given and true, it also adds the Coulomb potential.
  !
  ! INPUTS
  ! * type(particle)    :: teilchen  -- particle (already boosted to LRF)
  ! * logical, OPTIONAL :: addCoulomb -- flag, whether to add Coulomb
  !   (default: cf. addCoulombDefault)
  !
  ! or:
  !
  ! * integer                 :: ID  -- Id of baryon
  ! * integer                 :: IQ  -- Charge of baryon
  ! * real, dimension(0:3)    :: mom -- momentum of baryon in LRF
  ! * real, dimension(1:3)    :: pos -- position of baryon
  !
  ! If additional arguments are given, they reduce to define the medium:
  ! * med%density,med%densityProton,med%densityNeutron = rhoB,rhoB/2,rhoB/2
  ! * med%density,med%densityProton,med%densityNeutron = med
  ! In these cases, the position of the particle is not set properly.
  !
  ! Calling the routine with an optional 'density' argument is as in the first case, but
  ! giving 'density' also as an optional argument to the internal routine 'mediumAt'.
  !
  ! OUTPUT
  ! function value
  !****************************************************************************
  Interface potential_LRF
     Module Procedure &
          & potential_LRF_1,potential_LRF_2,potential_LRF_3,&
          & potential_LRF_4,potential_LRF_5, &
          & potential_LRF_1C,potential_LRF_2C,potential_LRF_3C,&
          & potential_LRF_4C,potential_LRF_5C
  End Interface

  !****************************************************************************
  !****f* potentialModule/scaPot
  ! NAME
  ! real function scapot(part,baremass_out,success, addCoulomb)
  !
  ! real function scapot(ID,IQ,mom,pos,baremass_out,success, addCoulomb)
  !
  ! PURPOSE
  ! Returns the scalar potential of a particle.
  ! If 'addCoulomb' is given and true, it includes the Coulomb potential in the
  ! calculation.
  !
  ! INPUTS
  ! * type(particle)          :: part
  ! * logical, OPTIONAL :: addCoulomb -- flag, whether to add Coulomb
  !   (default: cf. addCoulombDefault)
  !
  ! or:
  !
  ! * integer                 :: ID  -- Id of baryon
  ! * integer                 :: IQ  -- Charge of baryon
  ! * real, dimension(0:3)    :: mom -- momentum of baryon in LRF
  ! * real, dimension(1:3)    :: pos -- position of baryon
  !
  ! OUTPUT
  ! * real, optional          :: bareMass_out
  ! * logical, optional       :: success
  !
  ! NOTES
  ! Does not respect coulomb potential !
  !****************************************************************************
  Interface scapot
     Module Procedure scapot1, scapot2
  end Interface scapot


  logical, parameter :: addCoulombDefault = .true.

  logical :: debugFlag=.false.

  public :: potential_LRF
  public :: trueEnergy
  public :: massDetermination
  public :: scapot

contains

  !****************************************************************************
  ! cf. Interface 'potential_LRF'
  !****************************************************************************
  real function potential_LRF_1(teilchen)
    use particleDefinition

    type(particle),intent(in) :: teilchen

    potential_LRF_1=potential_LRF_1C(teilchen,addCoulombDefault)

  end function potential_LRF_1
  !-------------------------------------------------------------------------
  real function potential_LRF_1C(teilchen, addCoulomb)
    use particleDefinition
    use mesonPotentialModule, only: MesonPotential
    use baryonPotentialModule, only: BaryonPotential
    use mediumDefinition
    use mediumModule, only: mediumAt
    use callStack
    use IDTable, only: photon, isMeson, isBaryon
    use coulomb, only: emfoca

    type(particle),intent(in) :: teilchen
    logical,intent(in) :: addCoulomb

    type(medium) :: med

    potential_LRF_1C=0.

    ! The position of the particle is well defined, therefore we use the position to deduce the density
    med = mediumAt(teilchen%position)

    if (isBaryon(teilchen%ID)) then !All baryons
       potential_LRF_1C=BaryonPotential(teilchen,med,.false.)
    else if (isMeson(teilchen%ID)) then !All Mesons
       potential_LRF_1C=MesonPotential(teilchen,med)
    else if (teilchen%ID/=photon) then
       write(*,*) 'Funny particle in potential_LRF_1. ID=',teilchen%ID
       call traceback()
    end if

    if (addCoulomb) potential_LRF_1C = potential_LRF_1C &
         + emfoca(teilchen%position,teilchen%momentum(1:3),teilchen%charge,teilchen%ID)

  end function potential_LRF_1C
  !-------------------------------------------------------------------------
  real function potential_LRF_2(ID,IQ,mom,pos)
    integer, intent(in) :: ID
    integer, intent(in) :: IQ
    real, dimension(0:3), intent(in) :: mom
    real, dimension(1:3), intent(in) :: pos

    potential_LRF_2 = potential_LRF_2C(ID,IQ,mom,pos,addCoulombDefault)

  end function potential_LRF_2
  !-------------------------------------------------------------------------
  real function potential_LRF_2C(ID,IQ,mom,pos, addCoulomb)
    use particleDefinition
    use constants, only: mN
    use callStack, only: TRACEBACK

    integer, intent(in) :: ID
    integer, intent(in) :: IQ
    real, dimension(0:3), intent(in) :: mom
    real, dimension(1:3), intent(in) :: pos
    logical,intent(in) :: addCoulomb

    type(particle) :: part

    ! The notes marked with '<===' are included while converting the
    ! functionality of "potential_nucleon" into this routine.
    ! Is the mass and/or p(0) important here at all ????

    if (id.ne.1) call TRACEBACK("Use only with ID = 1") ! <=== necessary ???

    call SetToDefault(part)
    part%ID=ID
    part%charge=IQ
    part%antiparticle=.false.
    part%mass=mN ! <=== only with ID=1; necessary ???
    part%momentum(1:3)=mom(1:3) ! <=== ???
    part%position=pos
    part%antiparticle=.false.
    part%perturbative=.true.

    potential_LRF_2C = potential_LRF_1C(part, addCoulomb)

  end function potential_LRF_2C
  !-------------------------------------------------------------------------
  real function potential_LRF_3(teilchen,rhoB)
    use particleDefinition

    type(particle),intent(in) :: teilchen
    real, intent(in) :: rhoB

    potential_LRF_3 = potential_LRF_3C(teilchen,rhoB, addCoulombDefault)

  End function potential_LRF_3
  !-------------------------------------------------------------------------
  real function potential_LRF_3C(teilchen,rhoB, addCoulomb)
    use particleDefinition
    use mediumDefinition

    type(particle),intent(in) :: teilchen
    real, intent(in) :: rhoB
    logical,intent(in) :: addCoulomb

    type(medium) :: med

    ! The position of the particle is not well defined,
    ! only rhoB is given as input
    med%densityProton=rhoB/2.
    med%densityNeutron=rhoB/2.
    med%density=rhoB
    med%useMedium=.true.

    potential_LRF_3C = potential_LRF_4C(teilchen,med,addCoulomb)

  end function potential_LRF_3C
  !-------------------------------------------------------------------------
  real function potential_LRF_4(teilchen,med)
    use particleDefinition
    use mediumDefinition

    type(particle),intent(in) :: teilchen
    type(medium),intent(in) :: med

    potential_LRF_4=potential_LRF_4C(teilchen,med, addCoulombDefault)

  end function potential_LRF_4
  !-------------------------------------------------------------------------
  real function potential_LRF_4C(teilchen,med, addCoulomb)
    use particleDefinition
    use mesonPotentialModule, only: MesonPotential
    use baryonPotentialModule, only: BaryonPotential
    use mediumDefinition
    use callStack
    use IDTable, only: photon, isMeson, isBaryon

    type(particle),intent(in) :: teilchen
    type(medium),intent(in) :: med
    logical,intent(in) :: addCoulomb

    potential_LRF_4C=0.

    if (addCoulomb) then
       call traceback("Can not calculate Coulomb, since position not given.")
    end if

    if (isBaryon(Teilchen%ID)) then
       potential_LRF_4C=BaryonPotential(teilchen,med,.true.)
    else if (isMeson(Teilchen%ID)) then
       potential_LRF_4C=MesonPotential(teilchen,med)
    else if (Teilchen%ID/=photon) then
       write(*,*) 'Funny particle in potential_LRF_4. ID=',Teilchen%ID
       call traceback()
    end if

  end function potential_LRF_4C
  !-------------------------------------------------------------------------
  real function potential_LRF_5(teilchen,density)
    use particleDefinition
    use dichteDefinition

    type(particle),intent(in) :: teilchen
    type(dichte),intent(in) :: density

    potential_LRF_5=potential_LRF_5C(teilchen,density, addCoulombDefault)

  end function potential_LRF_5
  !-------------------------------------------------------------------------
  real function potential_LRF_5C(teilchen,density, addCoulomb)
    use particleDefinition
    use mesonPotentialModule, only: MesonPotential
    use baryonPotentialModule, only: BaryonPotential
    use mediumDefinition
    use mediumModule, only: mediumAt
    use callStack
    use IDTable, only: photon, isMeson, isBaryon
    use dichteDefinition
    use coulomb, only: emfoca

    type(particle),intent(in) :: teilchen
    type(dichte),intent(in) :: density
    logical,intent(in) :: addCoulomb

    type(medium) :: med

    potential_LRF_5C=0.

    ! The position of the particle is well defined, therefore we use the
    ! position to deduce the density
    med = mediumAt(density,teilchen%position)

    if (isBaryon(Teilchen%ID)) then !All baryons
       potential_LRF_5C=BaryonPotential(teilchen,med,.false.)
    else if (isMeson(Teilchen%ID)) then !All Mesons
       potential_LRF_5C=MesonPotential(teilchen,med)
    else if (Teilchen%ID/=photon) then
       write(*,*) 'Funny particle in potential_LRF_5. ID=',Teilchen%ID
       call traceback()
    end if

    if (addCoulomb) potential_LRF_5C = potential_LRF_5C &
         + emfoca(teilchen%position,teilchen%momentum(1:3),teilchen%charge,teilchen%ID)

  end function potential_LRF_5C
  !****************************************************************************

  !****************************************************************************
  ! cf. interface scapot
  !****************************************************************************
  real function scapot1(part,baremass_out,success, addCoulomb)
    use particleDefinition
    use minkowski, only: abs4,abs4Sq
    use densitymodule, only: densityAt
    use dichteDefinition

    type(particle),intent(in) :: part
    logical, intent(out), optional :: success
    real,    intent(out), optional :: bareMass_out !bare mass of the particle
    logical, intent(in), optional :: addCoulomb

    real :: vecpot, invMass, baremass
    type(dichte) :: dens

    integer, save :: zerocount=0
    logical :: flagOK, addC
    real, dimension(0:3) :: p

    ! Default return values:
    scapot1=0.
    if (present(success)) success=.false.
    if (present(baremass_out)) baremass_out=0.

    addC = addCoulombDefault
    if (present(addCoulomb)) addC = addCoulomb

    vecpot=potential_LRF(part, addC)

    p = part%momentum
    invMass=abs4(p,flagOK)
    if ((.not.flagOK).or.(invMass.le.1E-3)) return ! ==> failure

    p(0) = p(0)-vecpot
    baremass = abs4(p,flagOK)
    if (.not.flagOK) then
       zerocount=zerocount+1
       if (debugflag) then
          write(*,*) 'final state mass less than zero!',abs4Sq(p),&
               ' count = ',zerocount
          write(*,*) 'pot=',vecpot
          write(*,*) 'mom=',part%momentum(1:3)
          write(*,*) 'pos=',part%position(1:3)
          dens=densityAt(part%position(1:3))
          write(*,*) 'dens=',dens%baryon,dens%proton,dens%neutron
       end if
       return ! ==> failure
    end if

    scapot1=invMass-baremass

    if (present(success)) success=.true.
    if (present(baremass_out)) baremass_out=baremass

  end function scapot1
  !-------------------------------------------------------------------------
  real function scapot2(ID,IQ,mom,pos, baremass_out,success, addCoulomb)
    use particleDefinition
    use baryonPotentialmodule, only: getPotentialEQSType
    use minkowski, only: abs4, abs4Sq
    use densitymodule, only: densityAt
    use dichteDefinition

    integer, intent(in) :: ID
    integer, intent(in) :: IQ
    real, dimension(0:3), intent(in) :: mom
    real, dimension(1:3), intent(in) :: pos
    logical, intent(out), optional :: success
    real,    intent(out), optional :: bareMass_out !bare mass of the particle
    logical, intent(in), optional :: addCoulomb

    real :: vecpot, invMass, baremass
    type(dichte) :: dens
    type(particle) :: part

    integer, save :: zerocount=0
    logical :: flagOK, addC
    real, dimension(0:3) :: p

    ! Default return values:
    scapot2=0.
    if (present(success)) success=.false.
    if (present(baremass_out)) baremass_out=0.

    !define outgoing particle
    call setToDefault(part)
    part%ID=ID
    part%charge=IQ
    part%momentum=mom
    part%position=pos
    part%antiparticle=.false.
    part%perturbative=.true.


    ! determine bare mass of the particle,
    ! assume calculation frame = rest frame of nucleus = local rest frame

    if (getPotentialEQSType().eq.0) then
       scapot2=0.
       if (present(success)) success=.true.
       if (present(baremass_out)) baremass_out=abs4(mom)
       return
    end if

    addC = addCoulombDefault
    if (present(addCoulomb)) addC = addCoulomb

    vecpot=potential_LRF(part, addC)

    p = mom
    invMass=abs4(p,flagOK)
    if ((.not.flagOK).or.(invMass.le.1E-3)) return ! ==> failure

    p(0) = p(0)-vecpot
    baremass = abs4(p,flagOK)
    if (.not.flagOK) then
       zerocount=zerocount+1
       if (debugflag) then
          write(*,*) 'final state mass less than zero!',abs4Sq(p),&
               ' count = ',zerocount
          write(*,*) 'pot=',vecpot
          write(*,*) 'mom=',part%momentum(1:3)
          write(*,*) 'pos=',part%position(1:3)
          dens=densityAt(part%position(1:3))
          write(*,*) 'dens=',dens%baryon
          write(*,*) 'dens=',dens%proton
          write(*,*) 'dens=',dens%neutron
       end if
       return ! ==> failure
    end if

    scapot2=invMass-baremass

    if (present(success)) success=.true.
    if (present(baremass_out)) baremass_out=baremass

  end function scapot2

  !****************************************************************************
  !****s* potentialModule/trueEnergy
  ! NAME
  ! function trueEnergy(part)
  ! PURPOSE
  ! * (1) Boosts particle to "Local Rest Frame" (LRF).
  ! * (2) Evaluates energy-rearrangement terms.
  ! * (3) Boosts this back to calculation frame.
  ! INPUTS
  ! * type(particle) :: part -- Particle whose energy is to be calculated
  ! NOTES
  ! The LRF is the frame in which the baryon current vanishes.
  ! If the density is very small, then no boost takes place, and the free
  ! energy is assumed.
  ! For all particles besides the nucleon, the free 1-particle energy is
  ! equal to the true energy. But for the nucleon we need to substract the
  ! rearrangement terms due to the potential.
  !****************************************************************************
  real function trueEnergy(partIn, addCoulomb)
    use particleDefinition
    use mediumDefinition
    use densitymodule, only: boostToLRF
    use mediumModule, only: mediumAt
    use baryonPotentialModule, only: rearrangementPotential
    use IdTable, only: isMeson, isBaryon
    use coulomb, only: emfoca

    type(particle), intent(in) :: partIn
    logical,intent(in),OPTIONAL :: addCoulomb

    logical :: doC
    type(medium) :: med
    type(particle) :: part
    !real,parameter :: densityCutOff=1E-08

    part=partIn
    doC = .false.
    if (present(addCoulomb)) doC = addCoulomb

    if (isBaryon(Part%Id)) then
       med = mediumAt(part%position(1:3))  ! evaluate density
       call boostToLRF(part,1)  ! boost from calculation frame to LRF
       part%momentum(0) = part%momentum(0) &
            + rearrangementPotential(part, med)
       if (doC) part%momentum(0) = part%momentum(0) &
            - emfoca(part%position,part%momentum(1:3),part%charge,part%ID)/2
       call boostToLRF(part,2)  ! boost from LRF to calculation frame

    else if (isMeson(Part%Id)) then
       if (doC .and. part%ID.ne.0) then
          call boostToLRF(part,1)  ! boost from calculation frame to LRF
          part%momentum(0) = part%momentum(0) &
               - emfoca(part%position,part%momentum(1:3),part%charge,part%ID)/2
          call boostToLRF(part,2)  ! boost from LRF to calculation frame
       end if
    else
       trueEnergy=0.
       write(*,*) 'Funny particle in trueEnergy. ID=',Part%ID
       stop
    end if

    trueEnergy=part%momentum(0)

  end function trueEnergy


  !****************************************************************************
  ! cf. Interface 'massDetermination'
  !****************************************************************************
  subroutine massDetermination_noPart(ID,momLRF,med,baremass,verbose,success,pos)
    use particleDefinition
    use mediumDefinition
    use IdTable, only: isMeson
    use minkowski, only: abs4
    use selfEnergy_mesons, only: get_realPart

    integer, intent(in ) :: ID
    real, dimension(0:3) :: momLRF
    type(medium),intent(in) :: med
    real, intent(out)    :: baremass
    logical, optional    :: verbose,success
    real, dimension(1:3),intent(in),OPTIONAL :: pos

    real :: bareMassSquared
    type(particle) :: teilchen
    logical :: verbose_
    real :: rp, pot

    verbose_=.true.
    if (present(verbose)) verbose_=verbose

    if (present(success)) success=.false.

    ! Evaluate scalar potential in LRF:
    call setToDefault(teilchen)
    teilchen%ID      =  ID
    teilchen%momentum=  momLRF
    teilchen%perturbative = .false.

    if (present(pos)) teilchen%position = pos

    if (isMeson(teilchen%ID)) then
      rp = get_realPart(ID, abs4(momLRF), med)
    else
      rp = 0.
    end if

    ! p_0=sqrt(p**2+m**2)+V_LRF_0 (p(1:3)) =>

    if (present(pos)) then
       pot = potential_LRF(teilchen)
    else
       pot = potential_LRF(teilchen,med,.false.)
    end if
    bareMassSquared=(teilchen%momentum(0)-pot)**2-AbsMom(teilchen)**2-rp

    if (bareMassSquared.lt.0) then
       if (verbose_) then
          write(*,*) 'Problem: mass**2 less than zero in massDetermination_noPart'
          write(*,'(A,I3,A,4E15.4,A,2F9.3,A,F9.4)') 'ID:', ID, &
               & 'p=',momLRF,'Dens=',med%DensityProton,med%DensityNeutron,&
               'pot=',pot
       end if
       baremass=0.
       if (present(success)) success=.false.
    else
       bareMass=sqrt(max(0.,bareMassSquared))
       if (present(success)) success=.true.
    end if

  end subroutine massDetermination_noPart

  !-------------------------------------------------------------------------
  subroutine massDetermination_part(partIn, betaToCF, verbose, success)
    use IdTable, only: isHadron, isMeson
    use particleProperties, only: hadron
    use particleDefinition
    use densityModule, only: boostToLRF
    use lorentzTrafo, only: lorentz
    use minkowski, only: abs4
    use mediumDefinition
    use mediumModule, only: mediumAt
    use selfEnergy_mesons, only: get_realPart
    use callstack, only: traceBack

    type(particle), intent(inOut) :: partIn
    real,    optional, intent(in), dimension(1:3) :: betaToCF
    logical, optional, intent(in)  :: verbose
    logical, optional, intent(out) :: success

    logical            :: verbose_
    type(particle)     :: teilchen
    logical, parameter :: debugFlag=.false.
    type(medium) :: med
    real :: rp, h

    if (DebugFlag) write(*,*) '**In MassDetermination'

    verbose_=.true.
    if (present(verbose)) verbose_=verbose

    if (present(success)) success=.false.

    ! Check input
    if (AbsMom(partIn)>1000.) then
       call errorMessage
       if (isHadron(teilchen%ID)) then
          partIn%mass=hadron(teilchen%ID)%mass
       else
          partIn%mass=0.
       end if
       return
    end if

    ! Evaluate scalar potential in LRF:
    ! Boost a copy of the particle to LRF.
    teilchen=partIn

    ! (1.1) Boost to calculation frame
    if (present(betaToCF)) then  !boost particle to calculation frame
       call lorentz(betaToCF, teilchen%momentum, 'massDetermination')
    end if

    ! (1.2) Boost from calculation frame to LRF
    call boostToLRF(teilchen,1)

    ! p_0=sqrt(p**2+m**2)+V_LRF_0 (p(1:3)) =>

    h = (teilchen%momentum(0)-potential_LRF(teilchen))**2-AbsMom(teilchen)**2

    if (h<0.) then
       if (verbose_) write(*,*) 'WARNING in massDetermination: negative mass',teilchen%ID,h
       if (.not.present(success)) then
          write(*,*) 'Deleting particle ! teilchenIn%ID=0 !'
          partIn%ID=0 ! delete particle
       end if
    else
       if (isMeson(teilchen%ID)) then
          med = mediumAt(teilchen%position)
          rp = get_realPart(teilchen%ID, abs4(teilchen%momentum), med)
       else
          rp = 0.
       end if
       partIn%mass = sqrt(max(0.,h-rp))
       if (present(success)) success=.true.
    end if

  contains

    subroutine errorMessage
      use densitymodule, only: densityAt
      use dichteDefinition

      type(dichte) :: dens

      write(*,*)
      write(*,*) 'WARNING: Particle with huge momentum in MassDetermination'
      write(*,*) 'ID:           ', partIn%Id
      write(*,*) 'Charge:       ', partIn%charge
      write(*,*) 'Bare mass:    ', partIn%mass
      write(*,*) 'Eff. mass:    ', abs4(partIn%momentum)
      write(*,*) 'Momentum:     ', partIn%momentum
      write(*,*) 'Position:     ', partIn%position
      write(*,*) 'Perturbative: ', partIn%perturbative
      write(*,*) 'Offshell par: ', partIn%offshellParameter
      dens = densityAt(partIn%position)
      write(*,*) 'Density:      ', dens%baryon, dens%proton, dens%neutron
      if (present(betaToCF)) write(*,*) 'beta to calculation frame: ', betaToCF
      write(*,*)

      call traceBack('problem in massDetermination: particle with huge momentum!')

    end subroutine errorMessage

  end subroutine massDetermination_part



end module potentialModule
