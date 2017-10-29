!******************************************************************************
!****m* /lowElectron
! NAME
! module lowElectron
!
! PURPOSE
! This module includes a initilization routine for low energetic electron induced events
! (Resonance region)
!
! INPUTS
! The Namelist "LowElectron" in the Jobcard.
!******************************************************************************
module lowElectron
  use idTable, only: nres

  implicit none

  private

  !****************************************************************************
  !****g* lowElectron/runType
  ! SOURCE
  integer, save :: runType=1
  !
  ! PURPOSE
  ! * If runType=1, then we make runs at some fixed angle defined by lowElectron/theta_lf.
  ! * If runType=2, then we make runs at some fixed QSquared defined by lowElectron/QSquared
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/inputType
  ! SOURCE
  integer, save :: inputType=1
  !
  ! PURPOSE
  ! Decides which set of variables is used to determine the final electron
  ! energy energy_lf and the step size delta_energy_lf:
  ! * If inputType=1, then we use directly energy_lf and delta_energy_lf as input
  ! * If inputType=2, then we use W_min and W_max as input. For this we assume the
  !   nucleon to be at rest to calculate energy_lf out of W.
  ! * If inputType=3, then we use energy_lf_min and energy_lf_max as input.
  !****************************************************************************


  integer, parameter :: constantTheta=1
  integer, parameter :: constantQSquared=2

  !****************************************************************************
  !****g* lowElectron/theta_lf
  ! SOURCE
  real,save :: theta_lf=10.
  !
  ! PURPOSE
  ! Theta scattering angle of outgoing electron with respect to the incoming one.
  ! Only relevant of runType=1.
  !****************************************************************************


  !****************************************************************************
  !****g* lowElectron/phi_lf
  ! SOURCE
  real,save :: phi_lf=-10.
  !
  ! PURPOSE
  ! Phi scattering angle of outgoing electron with respect to the incoming one.
  ! If less than 0, then we do a Monte-Carlo-Integration over phi!
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/QSquared
  ! SOURCE
  real,save :: QSquared=0.5
  !
  ! PURPOSE
  ! QSquared of virtual photon.
  ! Only relevant of runType=2.
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/energy_li
  ! SOURCE
  real,save :: energy_li=1.2
  !
  ! PURPOSE
  ! Energy of incoming electron in GeV.
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/energy_lf
  ! SOURCE
  real,save :: energy_lf=0.8
  !
  ! PURPOSE
  ! Energy of final state electron in GeV.
  ! * Only used if inputType=1
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/delta_energy_lf
  ! SOURCE
  real,save :: delta_energy_lf=0.8
  !
  ! PURPOSE
  ! delta(Energy) of final state electron in GeV for energy scans.
  ! * Only used if inputType=1
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/W_min
  ! SOURCE
  real,save :: W_min=0.9
  !
  ! PURPOSE
  ! Minimal W at the hadronic vertex assuming a resting nucleon
  ! * Only used if inputType=2
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/W_max
  ! SOURCE
  real,save :: W_max=1.9
  !
  ! PURPOSE
  ! Maximal W at the hadronic vertex assuming a resting nucleon
  ! * Only used if inputType=2
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/energy_lf_min
  ! SOURCE
  real,save :: energy_lf_min=0.1
  !
  ! PURPOSE
  ! Minimal energy_lf
  ! * Only used if inputType=3
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/energy_lf_max
  ! SOURCE
  real,save :: energy_lf_max=0.1
  !
  ! PURPOSE
  ! Maximal energy_lf
  ! * Only used if inputType=3
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/Do_QE
  ! SOURCE
  logical,save :: Do_QE=.true.
  !
  ! PURPOSE
  ! Switch for including or excluding Quasi-Elastic (QE) processes
  !****************************************************************************



  !****************************************************************************
  !****g* lowElectron/minMass_QE
  ! SOURCE
  real,save :: minMass_QE=0.3
  !
  ! PURPOSE
  ! Minimal mass of a nucleon in QE event. Prevents super-luminous nucleons
  ! when embedded in a Skyrme potential.
  !****************************************************************************


  !****************************************************************************
  !****g* lowElectron/Do_1Pi
  ! SOURCE
  logical,save :: Do_1Pi=.true.
  !
  ! PURPOSE
  ! Switch for including or excluding direct Single pion production processes. If the
  ! resonances are included (Do_Res=.true.) then only the background part is included.
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/Do_2Pi
  ! SOURCE
  logical,save :: Do_2Pi=.true.
  !
  ! PURPOSE
  ! Switch for including or excluding direct Double pion production processes. If the
  ! resonances are included (Do_Res=.true.) then only the background part is included.
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/Do_DIS
  ! SOURCE
  logical,save :: Do_DIS=.true.
  !
  ! PURPOSE
  ! Switch for including or excluding deeply inelastic scattering (DIS) events. Only
  ! relevant for W>1.4-1.5 GeV.
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/Do_2p2hQE
  ! SOURCE
  logical,save :: Do_2p2hQE=.false.
  !
  ! PURPOSE
  ! Switch for including or excluding event according
  ! gamma* N1 N2 -> N1' N2'
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/Do_2p2hDelta
  ! SOURCE
  logical,save :: Do_2p2hDelta=.false.
  !
  ! PURPOSE
  ! Switch for including or excluding event according
  ! gamma* N1 N2 -> N' Delta
  !****************************************************************************


  !****************************************************************************
  !****g* lowElectron/Do_Res
  ! SOURCE
  logical,save :: Do_Res=.true.
  !
  ! PURPOSE
  ! Switch for including or excluding resonance production processes
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/useRes
  ! SOURCE
  logical,save,dimension(2:nres+1) :: useRes=.true.
  !
  ! PURPOSE
  ! Switch for including/excluding specific resonances
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/onlyDelta
  ! SOURCE
  logical,save  :: onlyDelta=.false.
  !
  ! PURPOSE
  ! Switch for including only delta resonance
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/nuclearTarget_corr
  ! SOURCE
  logical , save :: nuclearTarget_corr=.true.
  !
  ! PURPOSE
  ! * If the input is a nuclear targer, then the target nucleus is at rest
  !   and we calculate the cross section for
  !   nuclear target: use flux with respect to the nucleus.
  ! * Use .false. only for debugging.
  !****************************************************************************


  !****************************************************************************
  !****g* lowElectron/minEnergy_1pi
  ! SOURCE
  real, save :: minEnergy_1pi=0.16
  !
  ! PURPOSE
  ! Minimal q_0 such that pion production processes are considered.
  !****************************************************************************

  !****************************************************************************
  !****g* lowElectron/low_Ele_incXsection
  ! SOURCE
  !
  real, save :: low_Ele_incXsection = 0
  ! PURPOSE
  ! Collects at event generagtion the inclusive cross section for each run,
  ! which can then be easily processed by analysis routines
  !****************************************************************************

  logical,save :: initFlag=.true.

  real, save :: originXS(-2:9) ! collect cross section according origin_...


  public :: init_lowElectron
  public :: le_get_FirstEventRange
  public :: getfirstEvent, resetFirstEvent
  public :: le_get_energy_li, le_get_energy_lf
  public :: le_get_theta_lf
  public :: le_get_QSquared
  public :: low_Ele_incXsection
  public :: writeOriginXS


  ! Integer to take care of the firstEvent variable
  integer, save :: first=0

contains

  !****************************************************************************
  !****s* lowElectron/init_lowElectron
  ! NAME
  ! subroutine init_lowElectron(realParticles,pertParticles,raiseEnergy,targetNuc)
  !
  ! PURPOSE
  ! Initializes e^- A events.
  !
  ! INPUTS
  ! * logical :: raiseEnergy -- if true than we increase energy_lf by delta_energy_lf
  ! * type(particle), dimension(:,:) :: realParticles  -- real particle vector
  ! * type(particle), dimension(:,:) :: pertParticles  -- perturbative particle vector
  ! * type(tNucleus), pointer        :: targetNuc      -- target Nucleus
  ! OUPTUT
  ! * type(particle), dimension(:,:) :: pertParticles  -- perturbative particle vector
  !****************************************************************************
  subroutine init_lowElectron(realParticles,pertParticles,raiseEnergy,targetNuc)
    use nucleusDefinition
    use PauliBlockingModule, only: checkPauli
    use particleDefinition
    use IdTable, only: nucleon, nres
    use random, only: rn
    use insertion, only: setIntoVector
    use output, only: WriteParticleVector, subchapter!, writeparticle
    use degRad_conversion, only: radian
    use vector, only: absVec
    use inputGeneral, only: fullEnsemble
    use eN_eventDefinition, only: electronNucleon_event
    use eN_event, only: init_electronNucleon_event, nuclearFluxFactor_correction
    use leptonKinematics, only: evaluate_QSquared,evaluate_theta_lf,evaluate_epsilon,evaluate_W, buildElectrons
    use Electron_origin
    use eventGenerator_eN_lowEnergy, only: eventGen_eN_lowEnergy

    type(particle), intent(in)   , dimension(:,:) :: realParticles    ! real particles
    type(particle), intent(inout), dimension(:,:) :: pertParticles    ! perturbative particles

    type(tNucleus), pointer                       :: targetNuc        ! Target Nucleus
    logical, intent(in)                           :: raiseEnergy      ! If true than we increase energy_lf by delta_energy_lf


    logical,parameter :: debug_local=.false. ! Switch on local debugging

    type(particle), dimension(1:10) :: finalState
    type(electronNucleon_event) :: eN
    logical :: nuclearTarget, successFlag, Do_1Pi_switch
    integer :: numtries,numEvents, i,j, channel
    real :: phi, xsection
    real, dimension(0:3) :: electron_in,electron_out, electron_out_total
    logical :: flagOK
    logical, dimension(8) :: doC

    write(*,*)
    write(*,subchapter) 'Initializing electron induced events'


    call resetFirstEvent
    call le_whichOrigin(0,size(realParticles,dim=1)*size(realParticles,dim=2))

    if (initFlag) then
       call initInput

       originXS = 0.0

       !***********************************************************************
       !****o* lowElectron/lowElectron_protocol.log
       ! NAME
       ! file lowElectron_protocol.log
       ! PURPOSE
       ! Logs the kinematics being used in the electron nucleus initialization
       !***********************************************************************
       open(40,file='lowElectron_protocol.log')
       write(40,'(6A20)') 'energy_li [GeV]','energy_lf [GeV]',&
            & 'QSquared [GeV^2]','theta_lf [degrees]', 'W [GeV]', 'epsilon [dimensionless]'
       close(40)

       !***********************************************************************
       !****o* lowElectron/originXS.dat
       ! NAME
       ! file originXS.dat
       ! PURPOSE
       ! cross section seperated according origin, see Electron_origin.f90
       !
       ! given vales are 1/A d\sigma/dE´d\Omega given in mb/sr GeV
       !***********************************************************************
       open(40,file='originXS.dat')
       write(40,'(A)') '# XS is 1/A d\sigma/dE´d\Omega given in mb/sr GeV'
       write(40,'(12A20)') '# energy_lf [GeV]',&
            & 'events', &
            & 'total', 'QE', 'Res', '1pi', '2pi', 'DIS', &
            & 'VMDrho', '2p2h QE', '2p2h Delta', &
            & 'only Delta'
       close(40)


    end if

    ! Check that target nucleus is at rest
    if (targetNuc%mass.gt.1) then
       if (absVec(targetNuc%velocity).lt.1E-4) then
          nuclearTarget=nuclearTarget_corr
       else
          write(*,*) 'Error: Target nucleus is in motion in init_lowElectron. velocity=',targetNuc%velocity
          write(*,*) 'CRITICAL ERROR -> STOP'
          stop
       end if
    else
       nuclearTarget=.false.
    end if

    if (raiseEnergy) then
       energy_lf=energy_lf+delta_energy_lf
       ! CHECK that further energy increase is sensible:
       if (energy_lf.gt.energy_li) then
          write(*,*) 'Incoming electron energy smaller than outgoing'
          write(*,*) 'Energy incoming=',energy_li
          write(*,*) 'Energy outgoing=',energy_lf
          write(*,*) 'STOP'
          stop
       end if
       originXS = 0.0
    end if

    write(*,*) 'Electron energy incoming=',energy_li
    write(*,*) 'Electron energy outgoing=',energy_lf

    originXS(-2) = energy_lf

    select case (runType)
    case (constantQSquared)
       write(*,*) 'The Q^2 of the photon is chosen constant!'
       write(*,*) '* Q^2 (GeV^2) :', QSquared
       theta_lf=evaluate_theta_lf(QSquared,energy_li,energy_lf)
    case (constantTheta)
       write(*,*) 'The electron scattering angle is chosen constant!'
       write(*,*) '* Electron scattering angle    (degree) :', theta_lf
       if (theta_lf.lt.0)     write(*,*) '   -> Random theta (degree) (MC Integration)!!!'
       QSquared=evaluate_QSquared(theta_lf,energy_li,energy_lf)
    case default
       write(*,*) 'Wrong runType in lowElectron.f90', runType
       write(*,*) 'STOP!'
       stop
    end select

    ! Write kinematics to file
    open(40,file='lowElectron_protocol.log',position='Append')
    write(40,'(1P,6E20.5)') energy_li,energy_lf,QSquared,theta_lf,&
         & evaluate_W(theta_lf,energy_li,energy_lf),evaluate_epsilon(theta_lf,energy_li,energy_lf)
    close(40)

    numEvents=countParticles(realParticles,nucleon) !numEvents = NumNucleons*NumEnsembles

    if (debug_local) then
       open(100,file="ElectronMomenta.dat")
       write(100,'(2A8,2A12,2A48)') "i","j","theta","phi","p_in","p_out"
       open(101,file="W_free.dat")
       write(101,'(2A12)') "W", "W_free"
       electron_out_total=0.
       numTries=0
    end if

    low_Ele_incXsection = 0

    ! Loop over all target particles

    ensembleLoop: do j=lbound(realParticles,dim=2),ubound(realParticles,dim=2)
       indexLoop: do i=lbound(realParticles,dim=1),ubound(realParticles,dim=1)

          if (realParticles(i,j)%ID.ne.nucleon) cycle ! only nucleon targets allowed

          ! (1) Set up electron
          if (phi_lf.lt.0) then
             phi=rn()*360.
          else
             phi=phi_lf
          end if
          call buildElectrons(electron_in,electron_out,radian(theta_lf),radian(phi),energy_li,energy_lf,(/0.,0.,1./))

          ! (2) Set up electron-nucleon event

          call init_electronNucleon_event(eN, electron_in,electron_out,realParticles(i,j),flagOK)
          if (.not.flagOK) cycle indexLoop


          !call write_electronNucleon_event(eN)
          !stop

          if (debug_local) then
             ! Protocol electrons
             write(100,'(2I8,10E12.4)') i,j, theta_lf,phi,electron_in,electron_out
             electron_out_total=electron_out_total+electron_out
             numTries=numTries+1
             ! Protocol W
             write(101,'(2E12.5)') eN%W, eN%W_free
          end if

          ! (3) Generate Final State
          if (eN%boson%momentum(0).lt.minEnergy_1pi) then
             Do_1Pi_switch=.false.
          else
             Do_1Pi_switch=Do_1Pi
          end if

          doC = (/Do_QE,Do_Res,Do_1Pi_switch,Do_2Pi,Do_DIS, .false.,Do_2p2hQE, Do_2p2hDelta/)

          call eventGen_eN_lowEnergy(eN,doC,useRes,.false.,realParticles, &
               finalState,channel,successFlag,xsection)
          if (.not.successFlag) cycle indexLoop

          if (channel.eq.nucleon) then
             ! Cut on too low mass nucleons
             if (finalState(1)%mass.lt.minMass_QE) then
                finalState%ID=0
                successFlag=.false.
                cycle indexLoop
             end if
          end if

          successFlag = checkPauli(finalState,realParticles)
          if (.not.successFlag) cycle indexLoop

          if (nuclearTarget) then
             ! Nuclear flux correction
             finalState%perweight=finalState%perweight*nuclearFluxFactor_correction(realParticles(i,j)%momentum,electron_in)
             xsection= xsection*nuclearFluxFactor_correction(realParticles(i,j)%momentum,electron_in)
          end if

          finalState%perweight=finalState%perweight/float(numEvents) !resulting sigma.dat will have "cross section/nucleon"
          xsection= xsection/float(numEvents) !numEvents = NumNucleons*NumEnsembles
          low_Ele_incXsection = low_Ele_incXsection + xsection

          originXS(-1) = originXS(-1) + 1.0 ! number of events
          originXS(0) = originXS(0) + xsection ! total cross section

          select case (channel)
          case (nucleon) ! ===== QE Event =====
             originXS(1) = originXS(1) + xsection

          case (2) ! ===== Delta Event =====
             originXS(2) = originXS(2) + xsection
             originXS(9) = originXS(9) + xsection

          case (3:nres+1) ! ===== Resonance Event =====
             originXS(2) = originXS(2) + xsection

          case (origin_singlePi) ! ===== Pi N Event =====
             originXS(3) = originXS(3) + xsection

          case (origin_doublePi) ! ===== Pi Pi N Event =====
             originXS(4) = originXS(4) + xsection

          case (origin_DIS) ! ===== DIS Event =====
             originXS(5) = originXS(5) + xsection

          case (origin_vecmes_rho) ! ===== VMDrho Event =====
             originXS(6) = originXS(6) + xsection

          case (origin_2p2hQE) ! ===== 2p2h QE Event =====
             originXS(7) = originXS(7) + xsection

          case (origin_2p2hDelta) ! ===== 2p2h Delta Event =====
             originXS(8) = originXS(8) + xsection

          case default
             write(*,*) 'strange channel: ',channel

          end select


          ! Store the production channel via firstEvent
          finalState%firstEvent=getFirstEvent()

!!$          if (finalState(1)%firstEvent .eq.81) then
!!$             write(*,*) 'channel =',channel
!!$             call WriteParticle(6,99,finalState)
!!$             stop
!!$          end if

          call le_whichOrigin(1,finalState(1)%firstEvent,channel)
          ! (4) Set into vector

          if (fullEnsemble) then
             call setIntoVector(finalState,pertParticles,successFlag,.true.)
          else
             call setIntoVector(finalState,pertParticles(i:i,:),successFlag,.true.)
          end if

          if (.not.successFlag) then
             write(*,*) 'setIntoVector failed in init_singlePi. See pert* files! STOP! '
             call WriteParticleVector('pert',pertParticles)
             stop
          end if
       end do indexLoop
     end do ensembleLoop

     if (debug_local.and.numTries.gt.0) then
        ! Protocol averaged electron momenta
        write(100,*) '# Average electron momentum', electron_out_total/float(numTries)
     end if
     close(100); close(101)

    write(*,*)
    write(*,subchapter) '... finished initializing electron induced events'
    !stop
  end subroutine init_lowElectron


  !****************************************************************************
  !****f* lowElectron/le_get_energy_lf
  ! NAME
  ! real function le_get_energy_lf()
  !
  ! PURPOSE
  ! Returns energy_lf
  !****************************************************************************
  real function le_get_energy_lf()
    if (initFlag) call initInput
    le_get_energy_lf=energy_lf
  end function le_get_energy_lf

  !****************************************************************************
  !****f* lowElectron/le_get_energy_li
  ! NAME
  ! real function le_get_energy_li()
  !
  ! PURPOSE
  ! Returns energy_li
  !****************************************************************************
  real function le_get_energy_li()
    if (initFlag) call initInput
    le_get_energy_li=energy_li
  end function le_get_energy_li

  !****************************************************************************
  !****f* lowElectron/le_get_theta_lf
  ! NAME
  ! real function le_get_theta_lf()
  !
  ! PURPOSE
  ! Returns theta_lf
  !****************************************************************************
  real function le_get_theta_lf()
    if (initFlag) call initInput
    le_get_theta_lf=theta_lf
  end function le_get_theta_lf


  !****************************************************************************
  !****s* lowElectron/initInput
  ! NAME
  ! subroutine initInput
  !
  ! PURPOSE
  ! Reads in job card and checks the settings of the input parameters.
  !****************************************************************************
  subroutine initInput
    use IdTable, only: delta
    use output
    use inputGeneral, only: num_Energies
    use leptonKinematics, only: evaluate_theta_lf,evaluate_energy_lf

    integer :: ios

    !**************************************************************************
    !****n* lowElectron/LowElectron
    ! NAME
    ! NAMELIST /LowElectron/
    ! PURPOSE
    ! This Namelist includes:
    ! * runType
    ! * inputType
    ! * theta_lf
    ! * phi_lf
    ! * energy_li
    ! * energy_lf
    ! * energy_lf_min
    ! * energy_lf_max
    ! * delta_energy_lf
    ! * W_min
    ! * W_max
    ! * QSquared
    ! * Do_QE
    ! * Do_1Pi
    ! * Do_Res
    ! * Do_2Pi
    ! * Do_DIS
    ! * Do_2p2hQE
    ! * Do_2p2hDelta
    ! * minMass_QE
    ! * minEnergy_1pi
    ! * onlyDelta
    ! * nuclearTarget_corr
    !
    ! NOTES
    ! For runType=1 (fixed angle):
    ! * real,save :: theta_lf
    ! For runType=2 (fixed Qsquared):
    ! * real,save :: QSquared
    !
    ! For inputType=1 (energy_lf as input, fixed step size):
    ! * real,save :: energy_li
    ! * real,save :: energy_lf
    ! * real,save :: delta_energy_lf
    !
    ! For inputType=2 (W as input):
    ! * real,save :: W_min
    ! * real,save :: W_max
    ! * real,save :: delta_energy_lf
    !
    ! For inputType=3 (Energy_lf as input, step size according to number of runs):
    ! * real,save :: energy_lf_min
    ! * real,save :: energy_lf_max
    !
    !
    ! For Quasi-elastic scattering:
    ! * Do_QE
    ! * minMass_QE
    !
    ! For Pion production:
    ! * Do_1Pi
    ! * minEnergy_1pi
    !
    ! For Resonance excitations:
    ! * Do_Res
    ! * onlyDelta
    !
    ! For 2Pi background:
    ! * Do_2Pi
    !
    ! For DIS events:
    ! * Do_DIS
    !
    ! For 2p2h excitations:
    ! * Do_2p2hQE
    ! * Do_2p2hDelta
    !
    !
    ! Technical issues:
    ! *  nuclearTarget_corr
    !
    !**************************************************************************

    NAMELIST /lowElectron/ energy_li,energy_lf,theta_lf,delta_energy_lf,&
         & Do_QE, minMass_QE, &
         & Do_1Pi,minEnergy_1pi,&
         & Do_2Pi,Do_DIS, &
         & Do_2p2hQE, Do_2p2hDelta, &
         & Do_Res,onlyDelta,nuclearTarget_corr, &
         & runType, inputType, &
         & QSquared, energy_lf_min,energy_lf_max,w_min,w_max,phi_lf

    initFlag=.false.

    call Write_ReadingInput('lowElectron',0)

    rewind(5)
    read(5,nml=lowElectron,IOSTAT=ios)

    if (ios>0) then
       write(*,*) "Please note: "
       write(*,*) "we have renamed some switches:"
       write(*,*) "  QE_scattering   --> Do_QE"
       write(*,*) "  singlePi        --> Do_1Pi"
       write(*,*) "  doublePi_switch --> Do_2Pi"
       write(*,*) "  resonances      --> Do_Res"
       write(*,*) "  dis_switch      --> Do_DIS"
       write(*,*) "  minimalMass_QE  --> minMass_QE"
       write(*,*) "  singlePi_minEnergy --> minEnergy_1pi"
    end if

    call Write_ReadingInput("lowElectron",0,ios)

    write(*,*) '() Incoming Electron:'
    write(*,*) '   * Energy incoming:',energy_li

    ! Evaluate energy_lf and step size delta_energy_lf
    select case (inputType)
       case (1)
          write(*,*) '     We use energy_lf and delta_energy_lf as input!'
       case (2)
          write(*,*) '     We use W_min and W_max as input!'
          write(*,*) '   * W_min :', W_min
          write(*,*) '   * W_max :', W_max
          select case (runType)
          case (constantQSquared)
             energy_lf=evaluate_energy_lf(Qsquared=Qsquared,energy_li=energy_li,W=W_max)
             delta_energy_lf=(evaluate_energy_lf(Qsquared=QSquared,energy_li=energy_li,W=W_min)-energy_lf)/float(num_Energies-1)
          case (constantTheta)
             energy_lf=evaluate_energy_lf(theta_lf=theta_lf,energy_li=energy_li,W=W_max)
             delta_energy_lf=(evaluate_energy_lf(theta_lf=theta_lf,energy_li=energy_li,W=W_min)-energy_lf)/float(num_Energies-1)
          case default
             write(*,*) 'Wrong runType in lowElectron.f90', runType
             write(*,*) 'STOP!'
             stop
          end select
       case (3)
          write(*,*) '     We use energy_lf_min and energy_lf_max as input!'
          write(*,*) '   * energy_lf_min:' , energy_lf_min
          write(*,*) '   * energy_lf_max:' , energy_lf_max
          energy_lf=energy_lf_min
          if (num_Energies.eq.1) write(*,*) 'ATTENTION: num_Energies=1 !!!!!'
          delta_energy_lf=(energy_lf_max-energy_lf_min)/float(num_Energies-1)

       case default
          write(*,*) 'Wrong inputType in lowElectron.f90', inputType
          write(*,*) 'STOP!'
          stop
    end select

    write(*,*) '   * Energy outgoing:',energy_lf
    write(*,*) '   * Energy step size:',delta_energy_lf

    select case (runType)
    case (constantQSquared)
       write(*,*) '   The Q^2 of the photon is chosen constant!'
       write(*,*) '   * Q^2 (GeV^2) :', QSquared
       theta_lf=evaluate_theta_lf(QSquared,energy_li,energy_lf)
    case (constantTheta)
       write(*,*) '   The electron scattering angle  is chosen constant!'
       write(*,*) '   * Electron scattering angle    (degree) :', theta_lf
       if (theta_lf.lt.0)     write(*,*) '   -> Random theta (degree) (MC Integration)!!!'
    case default
       write(*,*) 'Wrong runType in lowElectron.f90', runType
       write(*,*) 'STOP!'
       stop
    end select
    write(*,*) '   Phi scattering angle (degree):' , phi_lf
    if (phi_lf.lt.0)     write(*,*) '   -> Random phi (degree) (MC Integration)!!!'

    if (.not.nuclearTarget_corr) write(*,*) '() WARNING: The nuclear target correction is switched off!!!!!!'
    write(*,*)
    write(*,*) '() Switches:'
    write(*,*) '      Do_QE        :',Do_QE
    write(*,*) '      Do_Res       :',Do_Res
    write(*,*) '      Do_1Pi       :',Do_1Pi
    write(*,*) '      Do_2Pi       :',Do_2Pi
    write(*,*) '      Do_DIS       :',Do_DIS
    write(*,*) '      Do_2p2hQE    :',Do_2p2hQE
    write(*,*) '      Do_2p2hDelta :',Do_2p2hDelta
    write(*,*)

    if (Do_1Pi.and.(.not.Do_Res)) then
       ! Only direct pion production
       write(*,*) '() Direct single pion production included.    Do_1Pi=',Do_1Pi
       write(*,*) '   Outgoing pion:'
       write(*,*) '    -> Random theta (degree) (MC Integration)!!!'
       write(*,*) '    -> Random phi   (degree) (MC Integration)!!!'
       write(*,*) '() Single pion production excluded if q_0<',  minEnergy_1pi

    else if (Do_Res) then
       write(*,*) '() Resonance production included.             Do_Res=',Do_Res
       if (onlyDelta) then
          write(*,*) '  -> No resonance besides the Delta is considered!!!!!!!!'
          useRes=.false.
          useRes(delta)=.true.
       end if
       if (Do_1Pi) then
          write(*,*) '() Direct single pion production included. Do_1Pi=',Do_1Pi
       else
          write(*,*) '() Direct single pion production EXCLUDED! Do_1Pi=',Do_1Pi
      end if
    else
       ! No resonances and no direct pion production
       write(*,*) '() Resonance production EXCLUDED!             Do_Res=',Do_Res
       write(*,*) '() Direct single pion production EXCLUDED!    Do_1Pi=',Do_1Pi
    end if


    if (Do_QE) then
       write(*,*) '() QE-scattering  included.                   Do_QE =', Do_QE
       write(*,*) '   Minimal Mass of a nucleon in QE-event =', minMass_QE
    else
       write(*,*) '() QE-scattering  EXCLUDED!                   Do_QE =',Do_QE
    end if

    if (Do_2Pi) then
       write(*,*) '() 2Pi bg-production  included.               Do_2pi=',Do_2Pi
    else
       write(*,*) '() 2Pi bg-production  EXCLUDED!               Do_2pi=',Do_2Pi
    end if
    if (Do_DIS) then
       write(*,*) '() DIS included.                              Do_DIS=',Do_DIS
    else
       write(*,*) '() DIS EXCLUDED!                              Do_DIS=',Do_DIS
    end if
    if (Do_2p2hQE) then
       write(*,*) '() 2p2h QE included.                       Do_2p2hQE=',Do_2p2hQE
    else
       write(*,*) '() 2p2h QE EXCLUDED!                       Do_2p2hQE=',Do_2p2hQE
    end if

    call Write_ReadingInput('lowElectron',1)

  end subroutine initInput

  !****************************************************************************
  !****s* lowElectron/resetFirstEvent
  ! NAME
  ! subroutine resetFirstEvent
  ! PURPOSE
  ! Resets the firstEvent counter to zero when starting a new run.
  !****************************************************************************
  subroutine resetFirstEvent
    first = 0
  end subroutine


  !****************************************************************************
  !****if* lowElectron/getFirstEvent
  ! NAME
  ! integer function getFirstEvent()
  ! PURPOSE
  ! Each event gets another entry in %firstEvent. This function is used to set firstEvent.
  !****************************************************************************
  integer function getFirstEvent()
    first=first+1
    getFirstEvent=first
  end function getFirstEvent


  !****************************************************************************
  !****f* lowElectron/le_get_FirstEventRange
  ! NAME
  ! function le_get_FirstEventRange()
  ! PURPOSE
  ! Each event gets another entry in %firstEvent. This entry is an integer value between 1 and number of events.
  ! This function returns the upper and lower limits of the range of %firstEvent.
  ! OUTPUT
  ! integer, dimension (1:2) :: getFirstEventRange
  !****************************************************************************
  function le_get_FirstEventRange()
    integer, dimension (1:2) :: le_get_FirstEventRange
    le_get_FirstEventRange=(/1,first/)
  end function le_get_FirstEventRange

  !****************************************************************************
  !****f* lowElectron/le_get_QSquared
  ! NAME
  ! real function le_get_Qsquared()
  ! PURPOSE
  ! Return the Q^2 value of the electron
  ! OUTPUT
  ! real :: le_get_QSquared()
  !****************************************************************************
  real function le_get_QSquared()
    le_get_QSquared= QSquared
  end function le_get_QSquared

  !****************************************************************************
  !****s* lowElectron/writeOriginXS
  ! NAME
  ! subroutine writeOriginXS
  ! PURPOSE
  ! write line to file 'originXS.dat'
  !****************************************************************************
  subroutine writeOriginXS(mulfak)

    real, intent(IN) :: mulfak

    open(40,file='originXS.dat',position='Append')
    write(40,'(1P,12E20.5)') originXS(-2:-1),originXS(0:9)*mulfak
    close(40)

  end subroutine writeOriginXS

end module lowElectron
