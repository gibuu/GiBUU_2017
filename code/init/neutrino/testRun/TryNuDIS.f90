program TryNuDIS

  use constants
  use inputGeneral
  use particleDefinition
  use Coll_gammaN
  use VMMassPythia
  use hadronFormation, only : forceInitFormation
  use Coll_Pythia
  use CollTools

  use eN_eventDefinition
  use eN_event
  use eventGenerator_eN_lowEnergy
  use eventGenerator_eN_HiEnergy
  use idTable,only :nres

  use ParamEP
  use output
  use minkowski, only: abs4,abs4Sq

  use collisionTerm, only: collideMain
  use preEventDefinition
  use PreEvList, only: CreateSortedPreEvent
  use lowElectron_origin, only : origin_singlePi,origin_doublePi,origin_DIS
  use insertion, only: GarbageCollection

  use Coll_nuN

  implicit none

  Type(electronNucleon_event), save :: eNev0
  Type(electronNucleon_event), save :: eNev
  type(particle) :: TargetNuc
  type(particle), dimension(1,1) :: realPart
  type(particle),dimension(1,1:25)  :: finalState

  real :: cross
  logical :: flagOK
  integer:: iEV, nEV=100

  ! Parameters to play with:

  real :: Ebeam = 1.0
  logical :: doCC = .true.
!  logical :: doCC = .false.
  logical :: doMassless=.true.
!  logical :: doMassless=.false.
  integer :: iGeneration = 2




  call readinputGeneral
  call init_database
  call forceInitFormation
!  call InitParticleProperties

  call SetSomeDefaults_PY

  !...Set up the target nucleon:

  call setToDefault(TargetNuc)
  TargetNuc%ID = 1
  TargetNuc%charge = 1
  TargetNuc%mass = 0.938
  TargetNuc%momentum = (/0.938, 0.0, 0.0, 0.0 /)
  TargetNuc%Position = 9999999d0


  !...Set up the initial kinematics:

  eNev0%electron_in = (/Ebeam,0d0, 0d0,Ebeam/)
  eNev0%electron_out= eNev0%electron_in ! !!!! DUMMY !!!!
  eNev0%photon      = eNev0%electron_in-eNev0%electron_out

  eNev0%nucleon%ID = 1
  eNev0%nucleon%mass  = 0.938
  eNev0%nucleon%charge = 1
  eNev0%nucleon%momentum  = (/0.938,    0d0, 0d0, 0d0/)
  eNev0%nucleon_free%position=999999999.

  eNev0%QSquared = -abs4Sq(eNev0%photon)
  eNev0%W = abs4(eNev0%photon+eNev0%nucleon%momentum)
    
  eNev0%nucleon_free      = eNev0%nucleon

  eNev0%W_free = abs4(eNev0%photon+eNev0%nucleon_free%momentum)

  call write_electronNucleon_event(eNev0,.FALSE.,.TRUE.)


  !...Do some Events:

  do iEV=1,nEV

     eNev = eNev0
     call eNev_init_Target(eNev,TargetNuc,flagOK)
     if (.not.flagOK) cycle

     call DoColl_nuN_Py(eNev,finalState(1,:),flagOK, doCC, iGeneration, doMassless, cross)

     if (flagOK) then
        call write_electronNucleon_event(eNev,.FALSE.,.FALSE.)
        call WriteParticle(6,1,finalState(1,:))
        write(*,*) 'cross section: (in mb)',cross
     else
        write(*,*) '---Event failed.---'
     endif

  end do

  write(*,*) 'done.'

end program TryNuDIS
