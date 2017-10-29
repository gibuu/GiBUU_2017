!******************************************************************************
!****m* /initTransportGivenParticle
! NAME
! module initTransportGivenParticle
!
! PURPOSE
! Initialisation of the particle vector for eventtype "TransportGivenParticle".
!******************************************************************************
module initTransportGivenParticle
  implicit none
  private

  !****************************************************************************
  !****g* initTransportGivenParticle/particle_ID
  ! SOURCE
  integer, save :: particle_ID=1
  ! PURPOSE
  ! Determines what kind of particle is initialized (see idTable)
  !****************************************************************************


  !****************************************************************************
  !****g* initTransportGivenParticle/charge
  ! SOURCE
  integer, save :: charge=1
  ! PURPOSE
  ! Determines what charge
  !****************************************************************************


  !****************************************************************************
  !****g* initTransportGivenParticle/position
  ! SOURCE
  real,dimension(1:3), save :: position=(/0.,0.,0./)
  ! PURPOSE
  ! Determines the position.
  !****************************************************************************


  !****************************************************************************
  !****g* initTransportGivenParticle/threemomentum
  ! SOURCE
  real,dimension(1:3), save :: threemomentum=(/0.,0.,1./)
  ! PURPOSE
  ! Determines the three-momentum.
  !****************************************************************************


  !****************************************************************************
  !****g* initTransportGivenParticle/mass
  ! SOURCE
  real, save :: mass=-1.
  ! PURPOSE
  ! Determines the mass (if negative, choose mass according to spectral function).
  !****************************************************************************


  !****************************************************************************
  !****g* initTransportGivenParticle/maxmass
  ! SOURCE
  real, save :: maxmass = 1.5
  ! PURPOSE
  ! Determines the maximum mass (if mass is chosen according to spectral function).
  !****************************************************************************


  !****************************************************************************
  !****g* initTransportGivenParticle/perweight
  ! SOURCE
  real, save :: perweight=1.
  ! PURPOSE
  ! Determines the weight.
  !****************************************************************************


  !****************************************************************************
  !****g* initTransportGivenParticle/frequency
  ! SOURCE
  integer, save :: frequency = 10
  ! PURPOSE
  ! after this amount of time steps a new output file is generated
  !****************************************************************************


  !****************************************************************************
  !****g* initTransportGivenParticle/initRandomRadiativeDelta
  ! SOURCE
  logical, save :: initRandomRadiativeDelta = .false.
  ! PURPOSE
  ! intented use: radiativeDelta decay.
  ! chooses position,threemomentum,mass of Delta randomly; charge is choosen either 0 or 1
  !****************************************************************************


  integer, save :: first

  logical, save :: initFlag=.true.

  public :: init_transportGivenParticle,particle_ID,getFrequency


contains

  !****************************************************************************
  !****s* initTransportGivenParticle/readinput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! This subroutine reads input out of jobcard
  ! from namelist 'TransportGivenParticle'.
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    integer :: ios

    !**************************************************************************
    !****n* initTransportGivenParticle/TransportGivenParticle
    ! NAME
    ! NAMELIST /TransportGivenParticle/
    ! PURPOSE
    ! This Namelist includes:
    ! * particle_ID
    ! * charge
    ! * position
    ! * threemomentum
    ! * mass
    ! * maxmass
    ! * perweight
    ! * frequency
    ! * initRandomRadiativeDelta
    !**************************************************************************
    NAMELIST /TransportGivenParticle/ particle_ID, charge, position, threemomentum, mass, maxmass, perweight, &
                                      frequency, initRandomRadiativeDelta

    call Write_ReadingInput('TransportGivenParticle',0)
    rewind(5)
    read(5,nml=TransportGivenParticle,IOSTAT=ios)
    call Write_ReadingInput("TransportGivenParticle",0,ios)
    call Write_ReadingInput('TransportGivenParticle',1)

    write(*,*) 'particle to be transported:'
    write(*,*) 'id: ', particle_ID

    if (initRandomRadiativeDelta .and. particle_ID/=2) then
       write(*,*) 'initRandomRadiativeDelta is true, but particle.ne.Delta -> STOP'
       stop
    end if

    if (initRandomRadiativeDelta) then
       write(*,*) 'choose position,threemomentum,mass of Delta randomly; charge is choosen either 0 or 1'
    else
       write(*,*) 'charge: ',charge
       write(*,*) 'position and three-momentum: ', position, threemomentum
       write(*,*) 'mass (if negative, use SF): ',mass
    end if

    write(*,*) 'max. mass: ', maxmass

    write(*,*) 'weight: ', perweight
    write(*,*) 'frequency: ', frequency

  end subroutine readInput


  integer function getFrequency()
    if (initFlag) then
       call readInput
       initFlag=.false.
    end if
    getFrequency = frequency
  end function getFrequency


  !****************************************************************************
  !****s* initTransportGivenParticle/init_transportGivenParticle
  ! NAME
  ! subroutine init_transportGivenParticle(realParticles,pertParticles)
  !
  ! PURPOSE
  ! This subroutine initializes a particle.
  !
  ! INPUTS
  ! * type(particle), dimension(:,:) :: realParticles
  !
  ! OUTPUT
  ! * type(particle), dimension(:,:) :: pertParticles
  !****************************************************************************
  subroutine init_transportGivenParticle(realParticles,pertParticles)
    use particleDefinition
    use random, only: rn,rnOmega
    use idtable, only: nucleon, isMeson, isBaryon
    use constants, only: pi
    use ParticleProperties, only: hadron
    use energyCalc, only: energyDetermination
    use propagation, only: updateVelocity
    use collisionNumbering, only: pert_numbering
    use insertion, only: setIntoVector
    use inputGeneral, only: fullEnsemble
    use output, only: WriteParticleVector
    use offShellPotential, only: getOffShellParameter
    use spectralFunc, only: specFunc
    use spectralFuncMesons, only: specFuncMes
    use hist2Df90
    use mediumDefinition
    use mediumModule, only: mediumAt
    use nucleusDefinition
    use nucleus, only: getTarget
    use RadiativeDeltaStorage, only: radiativeDeltaStorage_clear, radiativeDeltaStorage_Init, radiativeDeltaStorage_Store
    use transportGivenParticleAnalysis, only: transportGivenParticle_analyze

    type(particle), dimension(:,:),intent(inOut) :: realParticles
    type(particle), dimension(:,:),intent(inOut) :: pertParticles

    type(particle),dimension(1:1) :: finalstate
    integer :: i,j,n

    real, save :: minmass

    real :: baremass_out,polewidth
    real :: specFunc_RD,specFunc_max
    !real, dimension(0:3) :: momentum,momentum_dummy
    !    real, dimension(1:3) :: position_dummy

    real :: invmass, m !,width
    integer :: numEvents,numevents_failure

    logical :: setflag

    type(histogram2D) :: velos

    logical :: success
    type(medium) :: med

    type(tnucleus),pointer :: TargetNuc
    real :: radius,mom
    real, save :: extNucRadius
    type(particle) :: part
    real :: SF

    if (initFlag) then
       call readInput
       initFlag=.false.

       med = mediumAt(position)

       if (isMeson(particle_ID) .and. med%useMedium) then
          minmass = 0.
       else
          minmass = hadron(particle_ID)%minmass
       end if

       targetNuc => getTarget()
       extNucRadius=targetNuc%radius*1.2  !give 20% extra

    end if


    first=0

    call CreateHist2D(velos, 'velocities',(/0.,0.6/),(/10.,1.3/),(/0.05,0.01/))


    if (initRandomRadiativeDelta) then
       numEvents=0
       loopOverEnsemble0: do i = lbound(realParticles,dim=1),ubound(realParticles,dim=1)
          loopVector0: do j = lbound(realParticles,dim=2),ubound(realParticles,dim=2)
             if (realParticles(i,j)%ID.ne.nucleon) cycle
             numEvents=numEvents+1
          end do loopVector0
       end do loopOverEnsemble0
       call radiativeDeltaStorage_clear
       call radiativeDeltaStorage_Init(numEvents)
    end if

    numEvents=0
    numEvents_failure=0

    loopOverEnsemble: do i = lbound(realParticles,dim=1),ubound(realParticles,dim=1)
       loopVector: do j = lbound(realParticles,dim=2),ubound(realParticles,dim=2)
          if (realParticles(i,j)%ID.ne.nucleon) cycle

          call setToDefault(finalState)
          finalState%antiparticle=.false.
          finalState%perturbative=.true.
          finalState%productionTime=0.
          finalState%lastCollisionTime=0.
          finalState%formationTime=0.
          finalState%scaleCS=1.
          finalState%in_Formation=.false.
          finalState%event(1)=pert_numbering(realParticles(i,j))
          finalState%firstEvent=getFirstEvent()
          finalstate%history=0

          finalstate%id=particle_id

          numEvents=numEvents+1

          if (initRandomRadiativeDelta) then

             !choose charge between 1 and 0
             finalstate%charge=nint(rn())

             !choose radius
             radius=extNucRadius*rn()
             finalState(1)%position=radius*rnOmega()

             !choose threemomentum between 0.1 and 2 GeV
             mom= 0.1+1.9*rn()
             finalstate(1)%momentum(1:3)=mom*rnOmega()

             ! choose mass
             m=minmass+rn()*(maxmass-minmass)
             finalstate%mass=m

             call radiativeDeltaStorage_Store(finalState(1)%firstEvent,radius,mom,m,finalstate(1)%charge)

          else

             finalstate(1)%momentum(1:3)=threemomentum
             finalState(1)%position=position
             finalstate%charge=charge

             if (mass<0.) then
                !choose mass according to spectral function
                if (isBaryon(particle_ID)) then
                   polewidth=hadron(particle_ID)%width
                   if (particle_ID.eq.nucleon) polewidth=0.001
                   specfunc_max=1.5/(pi*hadron(particle_ID)%mass*polewidth)   !very very rough estimation -> CHECK
                   n=0
                   do
                      n=n+1
                      ! choose m**2 from a flat distribution (instead of m),
                      ! since specFunc is normalized in m**2
                      finalstate%mass=sqrt(minmass**2+rn()*(maxmass**2-minmass**2))
                      finalstate%offshellparameter=0.
                      call energyDetermination(finalstate(1))
                      specFunc_RD=specFunc(finalstate(1)%ID,finalstate(1)%charge,finalstate(1)%momentum, &
                                  finalstate(1)%position,bareMass_out)
                      if (specfunc_max<=specFunc_RD) then
                         write(*,*) 'problem in initTransportgivenParticle: specfunc_max.le.specFunc_RD ',specfunc_max,specFunc_RD
                         stop
                      end if
                      if (specfunc_max*rn()<=specFunc_RD) exit
                      !if(n.eq.1000) write(*,*) 'finding mass problem',i,j
                   end do
                   !write(*,*) 'successfull after ',n,' tries'
                   if (abs(finalstate(1)%mass-baremass_out)>=1.E-3) then
                      write(*,*) 'masses should be equal, stop'
                      stop
                   end if
                else if (isMeson(particle_ID)) then
                   polewidth=hadron(particle_ID)%width
                   specfunc_max=1.5/(pi*hadron(particle_ID)%mass*polewidth)   !very very rough estimation -> CHECK
                   n=0
                   do
                      n=n+1
                      ! choose m**2 from a flat distribution (instead of m),
                      ! since specFuncMes is normalized in m**2
                      finalstate%mass=sqrt(minmass**2+rn()*(maxmass**2-minmass**2))
                      finalstate%offshellparameter=0.
                      call energyDetermination(finalstate(1))
                      specFunc_RD=specFuncMes(finalstate(1),bareMass_out)
                      if (specfunc_max<=specFunc_RD) then
                         write(*,*) 'problem in initTransportgivenParticle: specfunc_max.le.specFunc_RD ',specfunc_max,specFunc_RD
                         stop
                      end if
                      if (specfunc_max*rn()<=specFunc_RD) exit
                      !if(n.eq.1000) write(*,*) 'finding mass problem',i,j
                   end do
                   !write(*,*) 'successfull after ',n,' tries'
                   if (abs(finalstate(1)%mass-baremass_out)>=1.E-3) then
                      write(*,*) 'masses should be equal, stop'
                      stop
                   end if
                else
                   write(*,*) 'particle is neither baryon nor meson, stop'
                   stop
                end if
             else
                finalstate%mass=mass
             end if
          end if  !initRandomRadiativeDelta

          finalstate%offshellparameter=0.
          call energyDetermination(finalstate(1))

          finalstate%offshellparameter=getOffShellParameter(finalstate(1)%ID,  &
               finalstate(1)%Mass,finalstate(1)%momentum,finalstate(1)%position,success)

          if (.not.success) then
             ! REJECT event
             finalstate%id=0
             numEvents_failure=numEvents_failure+1
             cycle loopVector
          end if
          !if flag useOffShellPotential is .false. then this returns 0.

          finalState%perweight=perweight

          call updateVelocity(finalState(1),success)
          if (.not.success) then
             write(*,*) 'problems in inittransportgivenparticle: velocity**2 greater or equal 1.'
             call AddHist2D(velos, (/finalstate(1)%offshellparameter,finalstate(1)%mass/),1.)
             finalstate%id=0
             numEvents_failure=numEvents_failure+1
             cycle loopVector
          end if

          call energyDetermination(finalstate(1))

          if (fullEnsemble) then
             call setIntoVector(finalState,pertParticles,setFlag)
          else
             call setIntoVector(finalState,pertParticles(i:i,:),setFlag)
          end if

          if (.not.setFlag) then
             write(*,*) 'error setIntoVector in initTransportGivenParticles'
             call WriteParticleVector('pert',pertParticles)
             stop
          end if
          cycle loopVector
          write(*,*) 'problem in init transport given particle, no event generated, Stop'
          stop
       end do loopVector
    end do loopOverEnsemble

    call transportGivenParticle_analyze (pertParticles, 0, particle_ID, frequency)

    part%ID = particle_ID
    part%charge = charge

    ! calculate vacuum SF
    part%position = 999. ! vacuum
    open(10,file='TGP_SF_vac.dat')
    write(10,*) '# vac SF'
    invMass=hadron(particle_ID)%minmass
    do
      invMass=invMass+0.001
      if (invMass>maxmass) exit
      part%momentum(1:3) = threemomentum(1:3)
      part%momentum(0)=sqrt(invMass**2+Dot_Product(threemomentum,threemomentum))
      if (isBaryon(particle_ID)) then
        SF = specFunc (particle_ID,charge,part%momentum,part%position)
      else if (isMeson(particle_ID)) then
        SF = specFuncMes (part)
      end if
      write(10,*) invMass,2.*invMass*SF
    end do
    close(10)

    ! calculate medium SF
    part%position = 0. ! medium
    open(10,file='TGP_SF_med.dat')
    write(10,*) '# med SF'
    invMass = 0.
    do
      invMass = invMass+0.001
      if (invMass>maxmass) exit
      part%momentum(1:3) = threemomentum(1:3)
      part%momentum(0)=sqrt(invMass**2+Dot_Product(threemomentum,threemomentum))
      if (isBaryon(particle_ID)) then
        SF = specFunc (particle_ID,charge,part%momentum,part%position)
      else if (isMeson(particle_ID)) then
        SF = specFuncMes (part)
      end if
      write(10,*) invMass,2*invMass*SF
    end do
    close(10)

    write(*,*) 'end of init transport given particle'
    open(199,file='velocities_greater_1.dat')
    call WriteHist2D_Gnuplot(velos,199,mul=velos%xbin(1)*velos%xbin(2))
    close(199)
    call RemoveHist2D(velos)

    write(*,*) 'Events, failed events:',numEvents, numEvents_failure
    open(199,file='velocities_greater_1.txt')
    write(199,'(A,I12,A,I12,A,G10.3,A)') 'Events, failed events [total, %]:',numEvents, '[',numEvents_failure,','&
         &  ,float(numEvents_failure)/float(numEvents)*100., '%]'
    close(199)
  end subroutine init_transportGivenParticle


  integer function getFirstEvent()
    first=first+1
    getFirstEvent=first
  end function getFirstEvent


end module initTransportGivenParticle
