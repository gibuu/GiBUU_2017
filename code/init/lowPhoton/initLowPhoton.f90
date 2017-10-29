!******************************************************************************
!****m* /initLowPhoton
! NAME
! module initLowPhoton
! PURPOSE
! Includes initialization of photon induced events at low energies.
!******************************************************************************
module initLowPhoton
  use particleDefinition
  use histMC
  implicit none
  private

  !****************************************************************************
  !****g* initLowPhoton/pascalTwoPi
  ! SOURCE
  !
  logical, save :: pascalTwoPi=.true.
  ! PURPOSE
  ! To switch on Pascal Muehlich's event generator for the gamma N -> N pi pi
  ! reaction:
  ! * .true.  : Pascal's prescription
  ! * .false. : phase space
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/equalDistribution_twoPi
  ! SOURCE
  !
  logical, save :: equalDistribution_twoPi=.false.
  ! PURPOSE
  ! All two pion production channels are populated with the same probability.
  ! This enhances the production of pi^0 pi0 pairs.
  ! Here, one has been careful to set perweight in a careful manner.
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/twoPi
  ! SOURCE
  !
  logical, save :: twoPi = .true.
  ! PURPOSE
  ! Switch for the direct 2pi production (gamma N -> N pi pi).
  ! If resonances=.true., then it is only a background, else the full cross
  ! section.
  ! Only works up to sqrt(s) = 2.1 GeV and is set to zero above.
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/vecMes
  ! SOURCE
  !
  logical, save :: vecMes(1:3) = .false.
  ! PURPOSE
  ! Switch for the production of neutral vector mesons: gamma N -> V N
  ! (1=rho^0,2=omega,3=phi)
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/vecMes_Delta
  ! SOURCE
  !
  logical, save :: vecMes_Delta(1:3) = .false.
  ! PURPOSE
  ! Switch for the production of neutral vector mesons: gamma N -> V Delta
  ! (1=rho^0,2=omega,3=phi)
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/resonances
  ! SOURCE
  logical, save :: resonances = .false.
  !
  ! PURPOSE
  ! Switch for including or excluding resonance production processes.
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/singlePi
  ! SOURCE
  logical, save :: singlePi = .false.
  !
  ! PURPOSE
  ! Switch for including direct single pion production. If resonances=.true.
  ! then it's only a background, else the full cross section.
  ! Only works up to sqrt(s) = 2.0 GeV and is set to zero above.
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/pi0Eta
  ! SOURCE
  logical, save :: pi0Eta = .false.
  !
  ! PURPOSE
  ! Switch for including direct pi^0 eta production.
  ! Only works up to sqrt(s) = 2.5 GeV and is set to zero above.
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/pi0Eta=_factor_neutron
  ! SOURCE
  real, save :: pi0eta_factor_neutron=1.
  !
  ! PURPOSE
  ! Scaling the gamma n->n pi^0 eta cross section.
  ! We assume
  !    sigma(gamma n->n pi^0 eta)
  !    = pi0eta_factor_neutron * sigma(gamma p->p pi^0 eta)
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/fritiof
  ! SOURCE
  logical, save :: fritiof = .false.
  !
  ! PURPOSE
  ! Switch for including FRITIOF events, cf. 'transitionEvent'.
  ! Will not give any contributions below sqrt(s) = 1.7 GeV.
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/includeResonance
  ! SOURCE
  logical,save,dimension(1:100) :: includeResonance=.true.
  !
  ! PURPOSE
  ! Switch for including/excluding specific resonances,
  ! e.g. includeResonance(21)=.false. excludes the D33(1700).
  !****************************************************************************


  logical,save,allocatable,dimension(:) :: useRes
  real,   save,allocatable,dimension(:) :: save_sigma_res


  !****************************************************************************
  !****g* initLowPhoton/onlyDelta
  ! SOURCE
  logical,save  :: onlyDelta=.false.
  !
  ! PURPOSE
  ! Switch for including only the Delta resonance.
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/onlyDelta_not
  ! SOURCE
  logical,save  :: onlyDelta_not=.false.
  !
  ! PURPOSE
  ! Switch for excluding only the Delta resonance.
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/debugFlag
  ! SOURCE
  !
  logical, save :: debugFlag=.false.
  ! PURPOSE
  ! To switch on debugging information.
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/energy_gamma
  ! SOURCE
  !
  real, save :: energy_gamma = 0.
  ! PURPOSE
  ! Energy of incoming photon in nucleus rest frame (in GeV).
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/delta_energy
  ! SOURCE
  !
  real, save :: delta_energy=0.
  ! PURPOSE
  ! Increase of energy for energy scans.
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/energy_weighting
  ! SOURCE
  !
  integer, save :: energy_weighting = 0
  ! PURPOSE
  ! Determines the relative weight of different photon energies in energy scans
  ! Possible values:
  ! * 0 = flat distribution (all energies are weighted equal)
  ! * 1 = exponential distr. (energies are weighted ~ 1/E)
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/noNucs_twoPi
  ! SOURCE
  !
  logical, save :: noNucs_twoPi=.false.
  ! PURPOSE
  ! Do not propagate the nucleons which are produced in a 2Pi event.
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/realRun
  ! SOURCE
  !
  logical, save :: realRun=.false.
  ! PURPOSE
  ! Do not initialize the final state particles as perturbative particles but
  ! as real ones.
  !****************************************************************************


  !****************************************************************************
  !****g* initLowPhoton/nuclearTarget_corr
  ! SOURCE
  logical , save :: nuclearTarget_corr=.true.
  !
  ! PURPOSE
  ! If the input is a nuclear target, then the target nucleus is
  ! assumed to be at rest and we calculate the cross section for a nuclear
  ! target.
  ! Therefore, we use the flux factor with respect to the nucleus.
  ! NOTES
  ! * Use ".false." only for debugging.
  !****************************************************************************

  !****************************************************************************
  !****g* initLowPhoton/makeHist_mass_energy
  ! SOURCE
  logical, save :: makeHist_mass_energy=.false.
  !
  ! PURPOSE
  ! Protocol the bare mass and the absolute momentum of excited resonances
  !
  ! Output is written to files named
  !   "Mass_and_energy_resID_'resID'_'energy of photon'_MeV.dat"
  !****************************************************************************


  !****************************************************************************
  !****g* initLowPhoton/useXsectionBoost
  ! SOURCE
  logical, save ::  useXsectionBoost=.true.
  !
  ! PURPOSE
  ! boost Xsections from nucleon rest frame to lab frame (where necessary)
  ! NOTES
  ! * Up to now only implemented for 2pi,
  !   not necessary for single pi and resonances
  ! * Check vecMes production !!!!
  !****************************************************************************


  logical, save :: gammaN_VN_formation = .true.


  !****************************************************************************
  !****g* initLowPhoton/selectFrame
  ! SOURCE
  integer, save :: selectFrame = 1
  !
  ! PURPOSE
  ! Select the frame in which 'sqrt(s)_free' is calulated:
  ! 1 = CM frame (default)
  ! 2 = LAB frame
  !
  ! In the chosen frame, one removes the 0th component of the potential from
  ! the nucleon momentum.
  ! Then this modified momentum is used to calculate sqrt(s)_free.
  !****************************************************************************


  !****************************************************************************
  !****g* initLowPhoton/printXS
  ! SOURCE
  logical, save :: printXS = .false.
  !
  ! PURPOSE
  ! Write out the photoproduction cross sections to files
  ! "lowPhoton_XS.dat" and "lowPhoton_XS_res.dat".
  !****************************************************************************

  logical, save :: initFlag = .true.

  ! Used to store coordinates of production points
  real, dimension(:,:),allocatable :: coordinate

  type(particle),save :: gamma ! the photon

  type(histogramMC),save:: vm_mass, vm_energy, vm_mom, vm_ekin, omega_excit_func

  ! for energy scans: save minimum and maximum energies, since 'energy_gamma'
  ! will be overwritten:
  real, save :: Emin, Emax

  public :: initialize_lowPhoton
  public :: energy_gamma
  public :: getCoordinate
  public :: gammaN_VN_formation
  public :: energyWeight
  public :: lowPhotonInit_getRealRun

  public :: getTwoPi
  public :: getVecMes

  public :: getSinglePi
  public :: getResonances
  public :: getEnergyGamma
  public :: getDeltaEnergy
  public :: getpi0eta

  public :: cleanup

contains

  !****************************************************************************
  !****s* initLowPhoton/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input out of jobcard. Namelist 'low_photo_induced'.
  !****************************************************************************
  subroutine readInput
    use particleProperties, only: hadron
    use IDTable, only: delta,lambda,photon,nbar
    use output, only: Write_ReadingInput, line
    use inputGeneral, only: num_Energies
    use Dilepton_Analysis, only: Dilep_Init, Dilep_UpdateProjectile
    use hadronFormation, only: forceInitFormation
    use collisionTerm, only: readInputCollTerm => readInput

    !**************************************************************************
    !****n* initLowPhoton/low_photo_induced
    ! NAME
    ! NAMELIST low_photo_induced
    ! PURPOSE
    ! Includes input parameters
    !
    ! Kinematical parameters:
    ! * energy_gamma
    ! * delta_energy
    ! * energy_weighting
    !
    ! Single Pion production:
    ! * singlePi
    !
    ! Two pion production:
    ! * twoPi
    ! * pi0eta
    ! * pi0eta_factor_neutron
    ! * pascalTwoPi
    ! * noNucs_twoPi
    ! * equalDistribution_twoPi
    !
    ! Resonance production:
    ! * resonances
    ! * onlyDelta
    ! * onlyDelta_not
    ! * includeResonance
    !
    ! Vector meson production:
    ! * vecMes
    ! * vecMes_Delta
    ! * gammaN_VN_formation
    ! * fritiof
    !
    ! Technical:
    ! * realRun
    ! * debugFlag
    ! * nuclearTarget_corr
    ! * useXsectionBoost
    ! * selectFrame
    ! * printXS
    !**************************************************************************
    NAMELIST /low_photo_induced/ &
         energy_gamma, delta_energy, energy_weighting, realRun, debugFlag, &
         singlePi, twoPi, pascalTwoPi, equalDistribution_twoPi, noNucs_twoPi, &
         resonances, onlyDelta, onlyDelta_not, includeResonance, &
         vecMes, vecMes_Delta, gammaN_VN_formation, &
         pi0eta, pi0eta_factor_neutron, fritiof, &
         nuclearTarget_corr, useXsectionBoost, selectFrame, printXS

    integer :: ios, i

    allocate(useRes(1:nbar))
    includeResonance=.true.
    useRes=.true.

    call Write_ReadingInput('low_photo_induced',0)

    rewind(5)
    read(5,nml=low_photo_induced,IOSTAT=ios)

    call Write_ReadingInput('low_photo_induced',0,ios)

    ! This construction is necessary since one can not use allocatable arrays
    ! in namelists:
    useRes=includeResonance(lbound(useRes,dim=1):ubound(useRes,dim=1))

    ! compute minimum and maximum energy
    Emin = energy_gamma
    Emax = Emin + delta_energy*(num_energies-1)

    if (realRun) then
       write(*,*) 'Final state particles are REAL particles!'
    else
       write(*,*) 'Final state particles are PERTURBATIVE particles!'
    end if
    write(*,*)
    if (.not.nuclearTarget_corr) then
       write(*,*) '() WARNING: The nuclear target correction is switched off!!!!!!'
       write(*,*) line
    end if
    write(*,*) 'Photon energy        =',energy_gamma
    write(*,*) 'Delta(Energy) per run=',delta_energy
    write(*,*) 'Energy weighting     =',energy_weighting
    write(*,*) 'Debugging :',debugFlag
    write(*,*) '* Use 1 pion production                  :',singlePi
    if (.not.singlePi) write(*,*) '   -> No direct 1Pi production!'
    write(*,*) '* Use 2 pion production                  :',twoPi
    if (twoPi) then
       if (resonances .and. pascalTwoPi) then
          write(*,*) '  NOTE: Pascal Muehlichs 2Pi not possible with resonances!'
          write(*,*) '        set pascalTwoPi=.false.!!!'
          pascalTwoPi=.false.
       end if

       write(*,*) '  Use Pascal Muehlichs 2Pi-Prescription  :',pascalTwoPi
       write(*,*) '  Use equal probability for 2pi channels :',equalDistribution_twoPi
       write(*,*) '  Do not propagate final state nucleons  :',noNucs_twoPi
    else
       write(*,*) '   -> No direct 2Pi production!'
    end if

    write(*,*) '* Use vector meson production (V N)      :',vecMes
    write(*,*) '* Use vector meson production (V Delta)  :',vecMes_Delta
    if (vecMes(1).or.vecMes(2).or.vecMes(3).or.vecMes_Delta(1).or.vecMes_Delta(2).or.vecMes_Delta(3)) &
         write(*,*) '  gammaN_VN_formation: ',gammaN_VN_formation
    write(*,*) '* Use pi^0 eta production                :',pi0eta
    if (.not.pi0eta)     write(*,*) '   -> No direct pi^0 eta production!'
    write(*,*) '* Use FRITIOF events                     :',fritiof

    write(*,*) '* Use Resonance model                    :',resonances
    if (resonances.and.onlyDelta) then
       write(*,*) '  -> No resonance besides the Delta is considered!!!!!!!!'
       useRes=.false.
       useRes(delta)=.true.
    else if (resonances.and.onlyDelta_not) then
       write(*,*) '  -> All resonance besides the Delta is considered!!!!!!!!'
       useRes(delta)=.false.
       ! Switch off exotic resonance:
       useRes(Lambda:)=.false.
    else if (resonances) then
       do i=1,nbar
          if (.not.useRes(i)) then
             write(*,*) " Warning: ", hadron(i)%name, " is not included!"
          end if
       end do
       ! Switch off exotic resonance:
       useRes(Lambda:)=.false.
    else
       ! Switch off exotic resonance:
       useRes(Lambda:)=.false.
    end if

    write(*,*)
    write(*,*) 'use Xsection Boost :', useXsectionBoost
    write(*,*) 'selectFrame = ', selectFrame
    write(*,*)

    !**************************************************************************
    !****o* initLowPhoton/init_LowPhoton.log
    ! NAME
    ! file init_LowPhoton.log
    ! PURPOSE
    ! Logfile for the lowPhoton init
    !**************************************************************************
    open(40, file="init_LowPhoton.log")
    write(40,*) '#Single Pion production:'
    write(40,*) 'singlePi',singlePi
    write(40,*)
    write(40,*) '#Two Pion production:'
    write(40,*) 'pascalTwoPi', pascalTwoPi
    write(40,*) 'equalDistribution_twoPi',equalDistribution_twoPi
    write(40,*) 'twoPi', twoPi
    write(40,*) 'noNucs_twoPi',noNucs_twoPi
    write(40,*)
    write(40,*) '#Vector meson production'
    write(40,*) 'Vecmes',vecMes
    write(40,*) 'Vecmes_Delta',vecMes_Delta
    write(40,*) 'gammaN_VN_formation',gammaN_VN_formation
    write(40,*)
    write(40,*) '#Resonance production'
    write(40,*) 'resonances',resonances
    write(40,*) 'onlyDelta',onlyDelta
    write(40,*) 'onlyDelta_not',onlyDelta_not
    do i=lbound(useRes,dim=1),ubound(useRes,dim=1)
       write(40,*) hadron(i)%name, " included? ",useRes(i)
    end do
    write(40,*)
    write(40,*) '# Other switches'
    write(40,*) 'nuclearTarget_corr',nuclearTarget_corr
    write(40,*) 'realRun', realRun
    write(40,*) 'useXsectionBoost',useXsectionBoost
    close(40)

    call Write_ReadingInput('low_photo_induced',1)

    ! Definition of the photon:
    gamma%ID=photon
    gamma%Charge=0
    gamma%momentum=(/energy_gamma,0.,0.,energy_Gamma/)

    call Dilep_Init(Emax)
    call Dilep_UpdateProjectile(energy_gamma,gamma%momentum)


    if (fritiof) call forceInitFormation

    call readInputCollTerm ! initialize enrgyCheck

    initFlag = .false.

  end subroutine readInput

  !****************************************************************************
  !****s* initLowPhoton/cleanUp
  ! NAME
  ! subroutine cleanUp
  ! PURPOSE
  ! cleanup allocated array and histograms
  !****************************************************************************
  subroutine cleanUp
    use Electron_origin, only: cleanupLEOrigin => cleanup
    use formfactors_A_main, only: cleanupMAID => cleanup
    use gamma2Pi_Xsections, only: cleanup2pi => cleanup
    if (allocated(useRes)) deallocate(useRes)
    if (allocated(save_sigma_res)) deallocate(save_sigma_res)
    if (allocated(coordinate)) deallocate(coordinate)
    call RemoveHistMC(vm_mass)
    call RemoveHistMC(vm_energy)
    call RemoveHistMC(vm_mom)
    call RemoveHistMC(vm_ekin)
    call RemoveHistMC(omega_excit_func)
    call cleanupLEOrigin
    if (singlePi) call cleanupMAID
    if (twoPi) call cleanup2pi
  end subroutine

  !****************************************************************************
  !****f* initLowPhoton/getResonances
  ! NAME
  ! logical function getResonances()
  ! PURPOSE
  ! return variable 'resonances'
  !****************************************************************************
  logical function getResonances()
    if (initFlag) call readInput
    getResonances=resonances
  end function getResonances

  !****************************************************************************
  !****f* initLowPhoton/getSinglePi
  ! NAME
  ! logical function getSinglePi()
  ! PURPOSE
  ! return variable 'singlePi'
  !****************************************************************************
  logical function getSinglePi()
    if (initFlag) call readInput
    getsinglePi=singlePi
  end function getSinglePi

  !****************************************************************************
  !****f* initLowPhoton/getPi0eta
  ! NAME
  ! logical function getPi0eta()
  ! PURPOSE
  ! return variable 'pi0eta'
  !****************************************************************************
  logical function getpi0eta()
    if (initFlag) call readInput
    getpi0eta=pi0eta
  end function getpi0eta

  !****************************************************************************
  !****f* initLowPhoton/getTwoPi
  ! NAME
  ! logical function getTwoPi()
  ! PURPOSE
  ! return variable 'twopi'
  !****************************************************************************
  logical function getTwoPi()
    if (initFlag) call readInput
    getTwoPi=twoPi
  end function getTwoPi

  !****************************************************************************
  !****f* initLowPhoton/getVecMes
  ! NAME
  ! logical function getVecMes()
  ! PURPOSE
  ! return logical 'or' of vecMes() and vecMes_Delta()
  !****************************************************************************
  logical function getVecMes()
    if (initFlag) call readInput
    getVecMes = vecMes(1) .or. vecMes(2) .or. vecMes(3) .or. vecMes_Delta(1) .or. vecMes_Delta(2) .or. vecMes_Delta(3)
  end function getVecMes

  !****************************************************************************
  !****f* initLowPhoton/lowPhotonInit_getRealRun
  ! NAME
  ! logical function lowPhotonInit_getRealRun()
  ! PURPOSE
  ! return variable 'realRun'
  !****************************************************************************
  logical function lowPhotonInit_getRealRun()
    if (initFlag) call readInput
    lowPhotonInit_getRealRun=realRun
  end function lowPhotonInit_getRealRun

  !****************************************************************************
  !****f* initLowPhoton/getEnergyGamma
  ! NAME
  ! real function getEnergyGamma()
  ! PURPOSE
  ! return variable 'energy_gamma'
  !****************************************************************************
  real function getEnergyGamma()
    if (initFlag) call readInput
    getEnergyGamma=energy_gamma
  end function getEnergyGamma

  !****************************************************************************
  !****f* initLowPhoton/getDeltaEnergy
  ! NAME
  ! real function getDeltaEnergy()
  ! PURPOSE
  ! return variable 'delta_energy'
  !****************************************************************************
  real function getDeltaEnergy()
    if (initFlag) call readInput
    getDeltaEnergy=delta_Energy
  end function getDeltaEnergy


  !****************************************************************************
  !****s* initLowPhoton/energyWeight
  ! NAME
  ! real function energyWeight(E)
  ! PURPOSE
  ! Computes a weighting factor which depends on the photon energy.
  ! Currently this weighting can be flat or exponential,
  ! cf. 'energy_weighting'.
  ! INPUTS
  ! * real, intent(in) :: E  --  the current photon energy
  ! NOTES
  ! All weighting functions implemented here should be normalized to the same
  ! value, i.e. doing a summation over all energies should yield a value of
  ! num_Energies, like in the flat case.
  !****************************************************************************
  real function energyWeight(E)
    use inputGeneral, only: num_Energies
    real, intent(in) :: E
    real, save :: norm = -1.
    integer :: i

    if (norm < 0. .and. energy_weighting == 1) then
      ! compute normalization factor
      norm = 0.
      do i = 0,(num_Energies-1)
        norm = norm + 1/(Emin + i*delta_energy)
      end do
      norm = num_Energies / norm
    end if

    select case (energy_weighting)
    case (0) ! flat
      energyWeight = 1.
    case (1) ! exponential
      energyWeight = norm / E
    case default
      write(*,*) "energy_weighting has an illegal value: ",energy_weighting
      stop
    end select
  end function


  !****************************************************************************
  !****s* initLowPhoton/initialize_lowPhoton
  ! NAME
  ! subroutine initialize_lowPhoton(realParticles, pertParticles, raiseFlag,
  ! targetNuc)
  !
  ! PURPOSE
  ! * Initializes a photon event on each nucleon in the realparticles vector.
  ! * The resulting particles are set into the pertParticles vector.
  ! * The assigned perweight is given for each single event by
  !   sigma_Total/(number of ensembles)/(number of nucleons per ensemble).
  ! * Pauli blocking is accounted for.
  !
  ! INPUTS
  ! * type(particle), dimension(:,:) :: realParticles
  ! * type(particle), dimension(:,:) :: pertParticles
  ! * logical :: raiseFlag  -- If .true. then the photon energy is raised by
  !   delta_energy
  ! * type(tNucleus), pointer :: targetNuc  -- target Nucleus
  !
  ! OUTPUT
  ! * type(particle), dimension(:,:) :: realParticles
  ! * type(particle), dimension(:,:) :: pertParticles

  ! NOTES
  ! * For the 1Pi Background, we allow for negative perweights. Henceforth,
  !   the Monte-Carlo sampling is done according to
  !      sum ( abs (cross sections) ).
  !   Then the perweight for a special channel "i" is given by
  !      sig(i)/abs(sig(i)) * sum ( abs (cross sections) ).
  ! * cross sections in [\mu b].
  !****************************************************************************
  subroutine initialize_lowPhoton(realParticles, pertParticles, raiseFlag, targetNuc)

    use nucleusDefinition
    use output, only:realToChar4
    use particleProperties, only: hadron
    use collisionTools, only: finalCheck
    use pauliBlockingModule, only: checkPauli
    use insertion, only: setIntoVector
    use IdTable, only: nbar, nucleon, Delta, pion, rho, omegaMeson, phi, photon
    use inputGeneral, only: fullensemble, numEnsembles
    use collisionNumbering, only: pert_numbering,real_numbering
    use gamma2Pi_Xsections, only: gamma2pi
    use energyCalc, only: energyDetermination, energyCorrection
    use mediumDefinition
    use mediumModule, only: mediumAt
    use constants, only: pi
    use random, only: rn, rn_openInterval,rnOmega_angles
    use resProd_lepton, only: sigma_resProd,sigma_pipi_res_vac,sigma_barmes_res_vac
    use deuterium_PL, only: deuterium_pertOrigin
    use photonPionProduction_medium, only:  dSigmadOmega_k_med
    use lowElectron, only: getFirstEvent, resetFirstEvent
    use Electron_origin
    use minkowski, only: SP,abs4,abs3
    use vector, only: absVec
    use offShellPotential, only: getOffShellParameter
    use propagation, only: updateVelocity
    use hist2Df90
    use photonXSections, only: calcXS_gammaN2VN, calcXS_gammaN2VDelta
    use lorentzTrafo, only: eval_sigmaBoost, lorentz
    use inputGeneral, only: numEnsembles,num_Runs_SameEnergy,current_run_number
    use history, only: setHistory
    use ClebschGordan
    use Dilepton_Analysis, only: Dilep_BH, Dilep_UpdateProjectile
    use pi0eta_photoproduction
    use XS_VMD, only: vmd,gvmd
    use constants, only: mN, mPi

    type(particle), dimension(:,:),intent(inOut),target :: realParticles
    type(particle), dimension(:,:),intent(inOut) :: pertParticles
    logical, intent(in) :: raiseFlag
    type(tNucleus),pointer :: targetNuc

    type(histogram2D), save, dimension(1:31) :: mass_energy_hist
    integer :: i, j, k, resID, iii, ID, ch, index
    real :: meff, pCM_Squared, srtFree, srts, phi_pi, theta_pi
    real, dimension(1:3) :: betaCMToLab,betaToLRF
    logical :: flag, setFlag, success, nuclearTarget
    real, dimension(0:3) :: sig2Pi,sigRes_2pi
    real, dimension(0:8) :: sigVM,sigRes_VM
    real, dimension(-1:1) :: sig1Pi
    real, dimension(0:4) :: sigVMD,XSvmd,XSgvmd
    real :: sigtot,sig,sigFr,sig_pi0eta
    type(particle), dimension(1:20) :: finalstate
    type(medium) :: mediumAtPosition
    real, dimension(1:nbar) :: mass_res , sigma_res
    integer :: numNucleons, whichReal_Nucleon, pionCharge, whichOrigin, nout
    real, dimension(-1:1,0:3) :: k_pi,pf_pi
    ! For immediate analysis of the cross sections and debugging:
    real, dimension(-1:1),save :: save_sig1Pi=0.
    real, save :: save_sigFr=0.
    integer, save :: numEvents_res=0, numEvents_1pi=0, numEvents_Fr=0
    integer, save :: failuresV=0, failuresO=0, failuresM=0
    real    :: fluxCorrection, BR_2pi, IS, integral
    !logical, save ::    once_Hist=.true.
    real, dimension(1:2) :: x1,x2,bin
    real, dimension(0:3) :: p_free ! Nucleon momentum without potential
    real, dimension(0:3) :: p_cm
    type(particle), pointer :: pPart

    if (.not.allocated(save_sigma_res)) then
       allocate(save_sigma_res(1:nbar))
       save_sigma_res=0.
    end if

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

    call resetFirstEvent
    call le_whichOrigin(0,size(realParticles,dim=1)*size(realParticles,dim=2))

    if (makeHist_mass_energy) then
       x1=(/0.,0.8/)
       x2=(/1.5,2.5/)
       bin=(/0.02,0.03/)
       do iii=lbound(mass_energy_hist,dim=1),ubound(mass_energy_hist,dim=1)
          call CreateHist2D(mass_energy_hist(iii), 'Momentum and mass of resonances',x1,x2,Bin,.true.)
       end do
    end if
    if (initFlag) call readInput
    if (raiseFlag) then
       energy_gamma=energy_gamma+delta_energy
       numEvents_1pi=0
       numEvents_res=0
       numEvents_Fr=0
       save_sig1Pi=0.
       save_sigma_res=0.
       save_sigFr=0.
    end if

    if (getVecMes() .or. fritiof) then
      if (.not. vm_mass%initialized) then
        ! initialize histogram for VM mass distribution
        call CreateHistMC(vm_mass,'VM Mass Spectra',0.,2.,0.001,3)
        vm_mass%xDesc='Mass [GeV]'
        vm_mass%yDesc(1:3) = (/ "rho  ", "omega", "phi  " /)
        call CreateHistMC(vm_energy,'VM Energy Spectra',0.,2.,0.01,3)
        vm_energy%xDesc='Energy [GeV]'
        vm_energy%yDesc(1:3) = vm_mass%yDesc(1:3)
        call CreateHistMC(vm_mom,'VM Momentum Spectra',0.,2.,0.01,3)
        vm_mom%xDesc='Momentum [GeV]'
        vm_mom%yDesc(1:3) = vm_mass%yDesc(1:3)
        call CreateHistMC(vm_ekin,'VM Kinetic Energy Spectra',0.,2.,0.01,3)
        vm_ekin%xDesc='Kinetic Energy [GeV]'
        vm_ekin%yDesc(1:3) = vm_mass%yDesc(1:3)
      end if
    end if

    if ((fritiof .or. vecMes(2) .or. vecMes_Delta(2)) .and. .not. omega_excit_func%initialized .and. delta_energy>0.) then
      call CreateHistMC (omega_excit_func, "omega excitation function (i.e. omega production cross section in microbarn/A)", &
                         Emin-delta_Energy/2., Emax+delta_Energy/2, delta_energy, 3)
      omega_excit_func%xDesc='E_gamma [GeV]'
      omega_excit_func%yDesc(1:3) = (/ "VN     ","VDelta ", "Fritiof"/)
    end if


    if (fullEnsemble.and.realRun) then
       write(*,*) 'FullEnsemble+Real particles in final state not yet implemented!!! STOP'
       stop
    end if

    write(*,'(A,I6,A,G15.7,A,G15.7)') &
         ' ** Low photon init: Run #', current_run_number, &
         ', photon energy = ',energy_gamma, &
         ', energy weight = ',energyWeight(energy_gamma)

    ! Generating the events
    integral=0.


    ! Definition of the photon:
    gamma%ID=photon
    gamma%Charge=0
    gamma%momentum=(/energy_gamma,0.,0.,energy_Gamma/)

    call Dilep_UpdateProjectile(energy_gamma,gamma%momentum)

    ensLoop: do i = lbound(realParticles,dim=1),ubound(realParticles,dim=1)
       if (realRun) then
          numNucleons=0
          do j = lbound(realParticles,dim=2),ubound(realParticles,dim=2)
             if (realParticles(i,j)%ID.ne.nucleon) cycle
             numNucleons=numNucleons+1
          end do
          if (numNucleons.gt.0) then
             ! Choose randomly one nucleon in  the ensemble to make the collsion with:
             whichReal_nucleon=1+int(rn()*numNucleons)
          else
             cycle ensLoop
          end if
          numNucleons=0
       end if

       partLoop: do j = lbound(realParticles,dim=2),ubound(realParticles,dim=2)
          if (realParticles(i,j)%ID.ne.nucleon) cycle

          if (realRun) then
             ! Only allow for a collision with the chosen nucleon.
             numNucleons=numNucleons+1
             if (.not.(numNucleons.eq.whichReal_nucleon)) cycle
          end if

          Nullify(deuterium_pertOrigin)
          deuterium_pertOrigin=>realParticles(i,j)

          pPart => realParticles(i,j)

          !********************************************************************
          ! Set up kinematics
          !********************************************************************
          meff=abs4(pPart%momentum)
          srts=sqrts(gamma,pPart)
          if (debugFlag) write(*,*) 'SQRT(s)=',srts
          pCM_Squared=(srts**2-meff**2)**2/(4.*srts**2)
          betaCMtoLab=-(pPart%momentum(1:3)+gamma%momentum(1:3))/(pPart%momentum(0)+energy_gamma)

          ! Nucleon momentum without (0th component of) potential:
          p_free(1:3)=pPart%momentum(1:3)
          p_free(0)  =sqrt(mN**2+Dot_product(p_free(1:3),p_free(1:3)))

          ! calculate srtfree (frame-dependent)
          select case (selectFrame)
          case (1)  ! calculate srtfree in CM frame
            srtfree=sqrt(mN**2+pCM_Squared)+sqrt(pCM_Squared)
          case (2)  ! calculate srtfree in LAB frame
            srtfree = abs4(p_free+gamma%momentum)
          case default
            write(*,*) "selectFrame = ",selectFrame
            stop
          end select

          if (debugFlag.or.twoPi.or.getVecMes()) then
             if ((i.eq.ubound(realParticles,dim=1)).and.(j.eq.uBound(realParticles,dim=2))) then
                call binSrts(srts,srtFree,.true.) ! Call and print results
             else if (realRun.and.(i.eq.ubound(realParticles,dim=1))) then
                call binSrts(srts,srtFree,.true.) ! Call and print results
             else
                call binSrts(srts,srtFree,.false.) ! Call w/o print results
             end if
          end if

          mediumAtPosition=mediumAt(pPart%position)

          !********************************************************************
          ! Add up cross sections (in microbarn)
          !********************************************************************

          sigTot=0.
          sig2pi=0.
          sigVM=0.
          sigRes_VM=0.
          sigma_res=0.
          sig1Pi=0.
          sig_pi0eta=0.
          sigVMD=0.
          sigFr=0.

          !********************************************************************
          !*** gamma N -> V N (direct photoproduction of vector mesons. V=rho^0,omega,phi)
          !********************************************************************
          if (vecMes(1) .or. vecMes(2) .or. vecMes(3)) then
            call calcXS_gammaN2VN(srtfree, mediumAtPosition, sigVM(1:4))
            if (resonances) then
              ! Subtract resonance contribution from V X channels (only important for rho)
              sigRes_VM(1:3) = sigma_barmes_res_vac(pPart,gamma%momentum,nucleon,(/rho,omegaMeson,phi/))*1000.
              do index=1,3
                sigVM(index)=max(0.,sigVM(index)-sigRes_VM(index))
              end do
            end if
            if (.not. vecMes(1)) sigVM(1) = 0.
            if (.not. vecMes(2)) sigVM(2) = 0.
            if (.not. vecMes(3)) sigVM(3) = 0.
          end if

          !********************************************************************
          !*** gamma N -> V Delta (direct photoproduction of vector mesons. V=rho^0,omega,phi)
          !********************************************************************
          if (vecMes_Delta(1) .or. vecMes_Delta(2) .or. vecMes_Delta(3)) then
            call calcXS_gammaN2VDelta(srtfree, sigVM(5:8), mediumAtPosition)
            if (resonances) then
              ! Subtract resonance contribution from V X channels (only important for rho)
              sigRes_VM(5:7) = sigma_barmes_res_vac(pPart,gamma%momentum,Delta,(/rho,omegaMeson,phi/))*1000.
              do index=5,7
                sigVM(index)=max(0.,sigVM(index)-sigRes_VM(index))
              end do
            end if
            if (.not. vecMes_Delta(1)) sigVM(5) = 0.
            if (.not. vecMes_Delta(2)) sigVM(6) = 0.
            if (.not. vecMes_Delta(3)) sigVM(7) = 0.
          end if

          sigVM(0) = SUM(sigVM(1:8))
          sigTot=sigTot+sigVM(0)

          !********************************************************************
          !*** gamma N -> pi pi + N
          !********************************************************************
          if (twoPi .and. srtfree<2.1) then
             !if(debugFlag) write(*,*)'vor 2pi-Querschnitt'
             call gamma2pi(pPart%charge,srtfree,sig2pi,betaCMToLab,mediumAtPosition,pPart%position)
             if (useXsectionBoost) sig2pi=sig2pi*eval_sigmaBoost(gamma%momentum,p_free)
             !if(debugflag) write(*,*)'nach 2pi-Querschnitt'
             if (resonances) then
                ! Subtract resonance contribution
                sigRes_2Pi=sigma_pipi_res_vac(realParticles(i,j),gamma%momentum)*1000.
                do index=lbound(sig2pi,dim=1),ubound(sig2pi,dim=1)
                   sig2Pi(index)=max(0.,sig2pi(index)-sigRes_2pi(index))
                end do
             end if
             if (getVecMes()) then
               ! Subtract vector meson contribution
               ! Problem: BR(V->2pi) is mass dependent, but the mass is not fixed.
               ! => use pole mass as an approximation.
               do index=1,3                                ! loop over vector mesons
                  ID = 101+2*index                         ! ID = 103,105,107
                  BR_2pi = hadron(ID)%Decays(1)        ! branching ratio (V->2pi) at pole mass
                  IS = float(hadron(ID)%isoSpinTimes2)/2.   ! isospin of vector meson
                  ! total :
                  sig2Pi(0)=max(0.,sig2Pi(0)-sigVM(index)*BR_2pi)
                  ! pi+ pi- :
                  sig2Pi(1)=max(0.,sig2Pi(1)-sigVM(index)*BR_2pi*(ClebschSquared(1.,1.,IS,1.,-1.)+ClebschSquared(1.,1.,IS,-1.,1.)))
                  ! pi+ pi0 and pi- pi0 have no contributions from neutral vector mesons
                  ! pi0 pi0 :
                  sig2Pi(3)=max(0.,sig2Pi(3)-sigVM(index)*BR_2pi*ClebschSquared(1.,1.,IS,0.,0.))
               end do
             end if
             sigTot=sigTot+sig2Pi(0)
          end if

          !********************************************************************
          !*** gamma N -> pi0 eta + N
          !********************************************************************
          if (pi0eta .and. srtfree<2.5) then
             if (pPart%charge==1) then
                ! Proton Xsection
                sig_pi0eta=sigma_gamma_p_to_p_pi0_eta(srtfree)
             else
                ! Neutron Xsection
                ! We just scale the proton Xsection since we have no model for this Xsection.
                sig_pi0eta=pi0eta_factor_neutron*sigma_gamma_p_to_p_pi0_eta(srtfree)
             end if
             sigTot=sigTot+sig_pi0eta
          end if

          !********************************************************************
          !*** Resonance production
          !********************************************************************
          mass_res =0.
          if (resonances) then
             numEvents_res=numEvents_res+1
             do resID=1,nbar
                if (useRes(resID)) &
                     sigma_res(resId)=sigma_resProd(pPart,resID,gamma%momentum,mass_res(resID))*1000.
             end do
             sigTot=sigTot+sum(sigma_res)
             save_sigma_res=save_sigma_res+sigma_res
          end if

          !********************************************************************
          !*** Direct 1-Pion production
          !********************************************************************
          if (singlePi) then
            if (.not.resonances) then
              numEvents_1pi=numEvents_1pi+1
              ! NOT BACKGROUND
              ! Choose by random theta and phi
              call rnOmega_angles(theta_pi,phi_pi)
              do pionCharge=pPart%charge-1,pPart%charge
                 sig1pi(pionCharge)=dSigmadOmega_k_med(pPart,pionCharge,phi_pi,theta_pi,&
                      gamma%momentum,k_pi(pionCharge,:),pf_pi(pionCharge,:), &
                      success)*4.*pi*1000.
              end do
              save_sig1pi=save_sig1pi+sig1pi
            else if (srtfree<2.0) then
              ! MAID input is only valid up to 2 GeV
              numEvents_1pi=numEvents_1pi+1
              ! Single Pi only as BACKGROUND
              call rnOmega_angles(theta_pi,phi_pi)
              do pionCharge=pPart%charge-1,pPart%charge
                sig1pi=sigma_1Pi_BG(pPart,phi_pi,theta_pi,gamma%momentum,k_pi,pf_pi)*4*pi*1000.
              end do
              save_sig1pi=save_sig1pi+sig1pi
            end if

            sigTot=sigTot+SUM(abs(sig1pi(:)))
          end if

          !********************************************************************
          !*** FRITIOF
          !********************************************************************
          if (fritiof .and. srtfree>mN+hadron(rho)%mass+0.01) then
            call vmd(srtfree,0.,0.,XSvmd)
            call gvmd(srtfree,0.,0.,XSgvmd)
            sigVMD = XSvmd + XSgvmd
            ! When FRITIOF is enabled, the total cross section is given by VMD.
            ! After subtracting the cross sections for the exclusive processes, the rest is left for FRITIOF.
            sigFr = max(sigVMD(0)-sigTot,0.)
            sigTot = sigTot + sigFr
            numEvents_Fr = numEvents_Fr + 1
            save_sigFr = save_sigFr + sigFr
          end if

          if (nuclearTarget) then

             ! nuclear flux correction:
             ! Multiply by relative velocity of photon and nucleon and divide
             ! by relative velocity of target and photon (=1)

             fluxCorrection=SP(gamma%momentum,pPart%momentum)/(pPart%momentum(0)*gamma%momentum(0))
             sigTot      =sigTot       *fluxCorrection
             sig2pi      =sig2pi       *fluxCorrection
             sigVM       =sigVM        *fluxCorrection
             sigma_res   =sigma_res    *fluxCorrection
             sig1Pi      =sig1Pi       *fluxCorrection
             sig_pi0eta  =sig_pi0eta   *fluxCorrection
          end if

          if (sigTot<1e-08) cycle
          integral=integral+sigTot

          !********************************************************************
          ! calculate Bethe-Heitler cross section (for Dilepton Analysis)
          !********************************************************************
          call Dilep_BH(gamma%momentum,pPart)


          !********************************************************************
          ! Choose Channel
          !********************************************************************

          ! Find random number which is not 1 or 0:
          sig=rn_openInterval()*sigTot

          !if(debugFlag) write(*,*) 'Vor choose channel', sig

          call setToDefault(finalState)
          if (twoPi .and. sig<sig2Pi(0)) then
             !---
             !--- 2pi production
             !---
             call twoPiProduction(pPart%charge,pPart%momentum,pPart%position, &
                  & energy_Gamma, sig2Pi, sigtot, srtFree, betaCMToLab, &
                  & mediumAtPosition, finalState(1:3), flag)
             whichOrigin=origin_doublePi

          else if (getVecMes() .and. sig<sig2Pi(0)+sigVM(0)) then
             !---
             !--- vector meson production
             !---
             call vecmesProduction(pPart, gamma, sigVM, sigTot, srts, srtFree,&
                  & betaCMToLab, mediumAtPosition, finalState(1:2), flag)
             if (finalState(1)%ID == Delta) then
                whichOrigin = origin_vecmes_Delta
             else
                select case (finalState(2)%ID)
                case (rho)
                   whichOrigin=origin_vecmes_rho
                case (omegaMeson)
                   whichOrigin=origin_vecmes_omega
                case (phi)
                   whichOrigin=origin_vecmes_phi
                case default
                   whichOrigin=origin_vecmes
                end select
             end if

          else if (resonances .and. sig<sig2Pi(0)+sigVM(0)+sum(sigma_res)) then
             !---
             !--- resonances
             !---
             sig=sig-sig2Pi(0)-sigVM(0)
             resLoop: do resID=1,nbar
                if (sig.le.sigma_res(resID)) then
                   ! Initialize finalstate with ID=resID and momentum p+q
                   finalState(1)%ID=resID
                   finalState(1)%momentum=realParticles(i,j)%momentum+gamma%momentum
                   finalState(1)%mass=mass_res(resID)
                   finalState(1)%charge=realParticles(i,j)%charge
                   finalState(:)%perweight=sigTot
                   whichOrigin=resID
                   if (resID.lt.ubound(mass_energy_hist,dim=1).and.makeHist_mass_energy) &
                        & call AddHist2D(mass_energy_hist(resID), &
                        & (/ sqrt(Dot_Product(finalState(1)%momentum(1:3),finalState(1)%momentum(1:3))), &
                        & mass_res(resID)/),1.,sigma_res(resID))
                   exit resLoop
                end if
                sig=sig-sigma_res(resID)
                if (resID.eq.nbar) then
                   write(*,*) 'Critical error in resonance production!',sigma_res,sig
                   write(*,*) 'STOP in initLowPhoton'
                   stop
                end if
             end do resLoop
             flag=.true.

          else if (singlePi .and. sig<sig2Pi(0)+sigVM(0)+sum(sigma_res)+sum(abs(sig1Pi))) then
             !---
             !--- single pion production
             !---
             sig=sig-sig2Pi(0)-sigVM(0)-sum(sigma_res)
             pionLoop: do pionCharge=-1,1
                if (sig.le.abs(sig1Pi(pionCharge))) then
                   ! Initialize pi N finalstates:
                   ! (1) PION:
                   finalState(1)%ID=pion
                   finalState(1)%momentum=k_pi(pionCharge,:)
                   finalState(1)%mass=mPi
                   finalState(1)%charge=pionCharge
                   ! (2) NUCLEON:
                   finalState(2)%ID=nucleon
                   finalState(2)%momentum=pf_pi(pionCharge,:)
                   finalState(2)%mass=mN
                   finalState(2)%charge=realParticles(i,j)%charge-pionCharge
                   finalState(:)%perweight=sigTot*sig1Pi(pionCharge)/abs(sig1Pi(pionCharge))
                   exit pionLoop
                end if
                sig=sig-abs(sig1Pi(pionCharge))
                if (pionCharge.eq.1) then
                   write(*,*) 'Critical error in single pion production!',sig1Pi,sig
                   write(*,*) 'STOP in initLowPhoton'
                   stop
                end if
             end do pionLoop
             flag=.true.
             !if(debugFlag)  write(*,*) 'pion charge=', finalState(1)%charge
             whichOrigin=origin_singlePi

          else if (pi0eta .and. sig<sig2Pi(0)+sigVM(0)+sum(sigma_res)+sum(abs(sig1Pi)+sig_pi0eta)) then
             !---
             !--- pi^0 eta production
             !---
             call event_gamma_p_to_p_pi0_eta(pPart,gamma,sigTot,srts,srtFree, &
                  & betaCMToLab,mediumAtPosition,finalState(1:3),flag,.not.realRun)
             whichOrigin=origin_pi0eta
          else if (fritiof) then
             !---
             !--- inclusive events with FRITIOF
             !---
             p_cm = pPart%momentum
             call lorentz(-betaCMToLab,p_cm) ! p_cm now is the nucleon momentum in the center-of-mass frame
             p_cm(1:3)=-p_cm(1:3)

             call transitionEvent(pPart, flag, finalState, srtFree, sigVMD, p_cm, -betaCMToLab, nout)

             do k=1,nout
                ! This is needed because energyCorrection takes finalState in the CM frame
                ! and returns it in the calculational frame:
                call lorentz(-betaCMToLab,finalState(k)%momentum)
             end do

             betaToLRF = 0.
             call energyCorrection(srts,betaToLRF,-betaCMToLab, mediumAtPosition, finalState(1:nout), flag)

             finalState(:)%perweight=sigTot
             whichOrigin = origin_fritiof

          else
             !---
             !--- ERROR
             !---
             write(*,*) 'Error in initLowPhoton.f90: Failed to determine production channel!', energy_gamma, srts
             write(*,*) singlePi, twoPi, vecMes, pi0eta, resonances
             write(*,*) sig, sigTot, sig1Pi(-1:1), sig2Pi(0), sigVM(0), sig_pi0eta
             write(*,*) sigma_res
             stop
          end if

          if (.not.flag) cycle

          !********************************************************************
          ! Check Pauli Blocking etc
          !********************************************************************

          ! set position
          do k=lbound(finalState,dim=1), ubound(finalState,dim=1)
             finalState(k)%position=pPart%position
          end do

          finalState(:)%perturbative=(.not.realRun)

          ! adjust perweights
          ! for real runs, we hit only one nucleon per ensemble,
          ! for perturbative runs, we hit each nucleon
          finalState(:)%perweight = finalState(:)%perweight / float(numEnsembles)*energyWeight(energy_gamma)
          if (.not.realRun) &
               finalState(:)%perweight = finalState(:)%perweight / float(targetNuc%mass)

          flag = checkPauli(finalState,realParticles)
          if (.not.flag) cycle

          !********************************************************************
          ! Determine omega excitation function
          !********************************************************************

          if ((vecmes(2) .or. vecmes_Delta(2) .or. fritiof) .and. delta_energy>0.) then
             do k=lbound(finalState,dim=1), ubound(finalState,dim=1)
                if (finalState(k)%ID == omegaMeson) then
                   if (whichorigin == origin_fritiof) then
                      ch = 3
                   else if (whichorigin == origin_vecMes_Delta) then
                      ch = 2
                   else
                      ch = 1
                   end if
                   ! add perweight to histogram (without weighting by energy)
                   call AddHistMC(omega_excit_func, energy_gamma, ch, &
                        finalState(k)%perweight / (float(num_Runs_SameEnergy)*energyWeight(energy_gamma)) )
                end if
             end do
          end if

          !********************************************************************
          ! Set into vector
          !********************************************************************

          do k=lbound(finalState,dim=1), ubound(finalState,dim=1)
             if (noNucs_TwoPi) then
                if (finalstate(k)%ID.eq.nucleon) finalState(k)%ID=0
             end if
             if (finalState(k)%ID<=0) cycle
             call energyDetermination(finalState(k))
             finalstate(k)%offshellparameter = &
                  getOffShellParameter(finalstate(k)%ID,  &
                  finalstate(k)%Mass,finalstate(k)%momentum, &
                  finalstate(k)%position,success)

             ! check offshell parameter
             if (.not.success) then
                failuresO=failuresO+1
                write(*,*) 'PROBLEM with offshell parameter in initLowPhoton, reaction not possible'
                write(*,*) 'Number of failures: ',failuresO
                cycle partLoop
             end if
             ! check p^mu p_mu
             if (abs4(finalState(k)%momentum)<0.) then
                failuresM=failuresM+1
                write(*,*) 'PROBLEM in initLowPhoton: m<0'
                write(*,*) 'Number of failures: ',failuresM
                cycle partLoop
             end if
             ! check velocity
             call updateVelocity(finalState(k),success)
             if (.not.success) then
                failuresV=failuresV+1
                write(*,'(A,4G12.5)')'PROBLEM in initLowPhoton: v>1'
                write(*,'(A,I5)')     'Number of failures: ',failuresV
                cycle partLoop
             end if
          end do

          if (.not.noNucs_twoPi) then
             flag = finalCheck((/pPart,gamma/), finalState, whichOrigin==origin_fritiof, 'initLowPhoton')
             if (.not. flag) cycle
          end if

          ! create unique particle number
          do k=lbound(finalState,dim=1), ubound(finalState,dim=1)
             if (finalState(k)%ID>0) call setNumber(finalState(k))
          end do

          !  Label event by eventNumber:
          if (.not.realRun) then
             finalState%event(1)=pert_numbering(realParticles(i,j))
             finalState%event(2)=pert_numbering(realParticles(i,j))
          else
             finalState%event(1)=real_numbering()
             finalState%event(2)=real_numbering()
          end if
          finalState%lastCollisionTime=0.
          ! This is important for the later analysis: f
          ! irstEvent must be unique for every event:
          finalState%firstEvent=getFirstEvent()
          ! Store the info about the origin of the event
          call le_whichOrigin(1,finalState(1)%firstEvent,whichorigin)

          call setHistory((/realParticles(i,j),gamma/), finalState)

          if (getvecmes() .or. fritiof) then
            ! put vector meson masses into histogram
            do k=lbound(finalState,dim=1), ubound(finalState,dim=1)
              select case (finalState(k)%ID)
                case (rho)
                call AddHistMC(vm_mass,abs4(finalState(k)%momentum),1, &
                               finalState(k)%perweight/float(num_Runs_SameEnergy*numEnsembles))
                call AddHistMC(vm_energy,finalState(k)%momentum(0),1, &
                               finalState(k)%perweight/float(num_Runs_SameEnergy*numEnsembles))
                call AddHistMC(vm_mom,abs3(finalState(k)%momentum),1, &
                               finalState(k)%perweight/float(num_Runs_SameEnergy*numEnsembles))
                call AddHistMC(vm_ekin,finalState(k)%momentum(0)-abs4(finalState(k)%momentum),1, &
                               finalState(k)%perweight/float(num_Runs_SameEnergy*numEnsembles))
              case (omegaMeson)
                call AddHistMC(vm_mass,abs4(finalState(k)%momentum),2, &
                               finalState(k)%perweight/float(num_Runs_SameEnergy*numEnsembles))
                call AddHistMC(vm_energy,finalState(k)%momentum(0),2, &
                               finalState(k)%perweight/float(num_Runs_SameEnergy*numEnsembles))
                call AddHistMC(vm_mom,abs3(finalState(k)%momentum),2, &
                               finalState(k)%perweight/float(num_Runs_SameEnergy*numEnsembles))
                call AddHistMC(vm_ekin,finalState(k)%momentum(0)-abs4(finalState(k)%momentum),2, &
                               finalState(k)%perweight/float(num_Runs_SameEnergy*numEnsembles))
                !!! write out mass distribution vs. sqrt(s) to reproduce Muehlich fig. 9.14.
!                 write(666,'(i12,5G13.6)') finalState(k)%number, srts, srtFree, abs4(finalState(k)%momentum), &
!                                           abs3(finalState(k)%momentum), finalState(k)%momentum(0)
                call omegaProdInfo (finalState(k), mediumAtPosition)
              case (phi)
                call AddHistMC(vm_mass,abs4(finalState(k)%momentum),3, &
                               finalState(k)%perweight/float(num_Runs_SameEnergy*numEnsembles))
                call AddHistMC(vm_energy,finalState(k)%momentum(0),3, &
                               finalState(k)%perweight/float(num_Runs_SameEnergy*numEnsembles))
                call AddHistMC(vm_mom,abs3(finalState(k)%momentum),3, &
                               finalState(k)%perweight/float(num_Runs_SameEnergy*numEnsembles))
                call AddHistMC(vm_ekin,finalState(k)%momentum(0)-abs4(finalState(k)%momentum),3, &
                               finalState(k)%perweight/float(num_Runs_SameEnergy*numEnsembles))
              end select
            end do
          end if

          if (twoPi) call saveCoordinate(finalState(lbound(finalState,dim=1))%firstEvent,realParticles(i,j)%position)

          !If(DebugFLAG) write(*,*) FinalState%ID
          !If(DebugFLAG) write(*,*) FinalState%momentum(0)-mN

          ! Create final state particles
          if (fullensemble) then
            if (realRun) stop 'RealRun+fullEnsemble not yet implemented: initLowPhoton'
            call setIntoVector(finalState, pertParticles, setFlag, .true.)
            ! (4) Check that setting into perturbative particle vector worked out
            if (.not.setFlag) then
                write(*,*) 'Perturbative particle vector too small!'
                write(*,*) size(finalState),  lbound(pertParticles), ubound(pertParticles)
                write(*,*) finalState%ID
                stop ' initLowPhoton'
            end if
          else
            if (realRun) then
                ! Delete the scattering partner
                realParticles(i,j)%ID=0
                ! Set new particles into real particle vector
                realParticles(i:i,:)%perweight=finalState(1)%perweight
                realParticles(i:i,:)%firstEvent=finalState(1)%firstEvent
                call setIntoVector(finalState,realParticles(i:i,:), setFlag, .true.)
                ! (4) Check that setting into real particle vector worked out
                if (.not.setFlag) then
                  write(*,*) 'Real particle vector too small!'
                  write(*,*) size(finalState),  lbound(realParticles), ubound(realParticles)
                  write(*,*) finalState%ID
                  stop ' initLowPhoton'
                end if
            else
                call setIntoVector(finalState,pertParticles(i:i,:), setFlag, .true.)
                ! (4) Check that setting into perturbative particle vector worked out
                if (.not.setFlag) then
                  write(*,*) 'Perturbative particle vector too small!'
                  write(*,*) size(finalState),  lbound(pertParticles), ubound(pertParticles)
                  write(*,*) finalState%ID
                  stop ' initLowPhoton'
                end if
            end if
          end if

       end do partLoop
    end do ensLoop

    Nullify(deuterium_pertOrigin)

    !**************************************************************************
    ! analysis
    !**************************************************************************

    ! write out histogram: momentum and mass of resonances
    if (makeHist_mass_energy) then
       do resID=lbound(mass_energy_hist,dim=1),ubound(mass_energy_hist,dim=1)
          if (useRes(resID)) then
             open(87,file='Mass_and_energy_resID_'//trim(hadron(resID)%name) //'__E_gamma_'//&
                  & realTochar4(1000.*energy_gamma)//'MeV.dat')
             call WriteHist2D_Gnuplot(mass_energy_hist(resID),87)
             call RemoveHist2D(mass_energy_hist(resID))
             close(87)
          end if
       end do
    end if

    ! print total cross section
    if (realRun) then
      integral = integral / float(numEnsembles)
    else
      integral = integral / float(numEnsembles*targetNuc%mass)
    end if
    write(*,'(A,G12.4,A)') 'After lowPhoton init... (sigma_tot=', integral,')'

    ! write out histogram: mass distribution of vector mesons
    if (getvecMes() .or. fritiof) then
      call WriteHistMC(vm_mass,  'VMmass_init.dat',   .false.)
      call WriteHistMC(vm_energy,'VMenergy_init.dat', .false.)
      call WriteHistMC(vm_mom,   'VMmom_init.dat',    .false.)
      call WriteHistMC(vm_ekin,  'VMekin_init.dat',   .false.)
    end if


    !**************************************************************************
    !****o* initLowPhoton/OmegaExcitFunc.dat
    ! NAME
    ! file OmegaExcitFunc.dat
    ! PURPOSE
    ! Contains the omega excitation function, i.e. the energy-dependent
    ! inclusive omega photo-production cross section. In contrast to
    ! OmegaExcitFunc_pi0gamma.dat, this does not include FSI effects.
    !**************************************************************************

    ! write out histogram: omega excitation function
    if (vecMes(2) .or. vecMes_Delta(2) .or. fritiof) &
         call WriteHistMC(omega_excit_func, 'OmegaExcitFunc.dat',.false.)

    if (printXS) call printOutXS

  contains


    !**************************************************************************
    !****is* initLowPhoton/printOutXS
    ! NAME
    ! subroutine printOutXS
    !
    ! PURPOSE
    ! ...
    !**************************************************************************
    subroutine printOutXS

      logical,save :: firstTime=.true.
      real :: onepi(-1:1),res(1:nbar),fr

      if (firstTime) then
        open(300,file='lowPhoton_XS.dat')
        open(301,file='lowPhoton_XS_res.dat')
        firsttime=.false.
      else
        open(300,file='lowPhoton_XS.dat',position='Append')
        open(301,file='lowPhoton_XS_res.dat',position='Append')
      end if

      if (numEvents_1pi>0) then
        onepi = save_sig1pi/float(numEvents_1pi)
      else
        onepi = 0.
      end if

      if (numEvents_res>0) then
        res = save_sigma_res(1:nbar)/float(numEvents_res)
      else
        res = 0.
      end if

      if (numEvents_Fr>0) then
        fr = save_sigFr/float(numEvents_Fr)
      else
        fr = 0.
      end if

      !************************************************************************
      !****o* initLowPhoton/lowPhoton_XS.dat
      ! NAME
      ! file lowPhoton_XS.dat
      ! PURPOSE
      ! Contains various photoproduction cross sections per nucleon
      ! (without FSI).
      ! PURPOSE
      ! Columns:
      !  * 1    : E_gamma [Gev]
      !  * 2    : sqrt(s) [GeV]
      !  * 3    : total cross section [microbarn/A]
      !  * 4    : sum of all resonances
      !  * 5-7  : exclusive pi-, pi0, pi+
      !  * 8-11 : 2pi
      !  * 12-20: vector mesons
      !  * 21   : pi0 eta
      !  * 22   : Fritiof
      !************************************************************************
      write(300,'(26G13.6)') energy_gamma, srts, integral, sum(res), &
           onepi(-1:1), sig2pi(0:3), sigVM(0:8), sig_pi0eta, fr

      !************************************************************************
      !****o* initLowPhoton/lowPhoton_XS_res.dat
      ! NAME
      ! file lowPhoton_XS_res.dat
      ! PURPOSE
      ! Contains the resonance photoproduction cross sections per nucleon
      ! (without FSI).
      ! PURPOSE
      ! Columns:
      !  * 1 : E_gamma [Gev]
      !  * 2 : sqrt(s) [GeV]
      !  * 3-63: resonance cross sections [microbarn/A]
      !************************************************************************
      write(301,'(63G13.6)') energy_gamma, srts, res(1:nbar)

      close(300)
      close(301)

    end subroutine

  end subroutine initialize_lowPhoton

  !****************************************************************************
  !****s* initLowPhoton/saveCoordinate
  ! NAME
  ! subroutine saveCoordinate(firstEventNumber,x)
  ! PURPOSE
  ! Used to save production coordinates of particles
  ! INPUTS
  ! * integer :: firstEventNumber -- firstEventnumber associated to position x
  ! * real, dimension(1:3) ::  x -- position
  !****************************************************************************
  subroutine saveCoordinate(firstEventNumber,x)

    real, dimension(1:3), intent(in) ::  x
    integer, intent(in) :: firstEventNumber
    real, dimension(:,:),allocatable :: save_coordinate
    logical, save :: initFlag_coord=.true.
    integer :: i

    if (initFlag_coord) then
       allocate(coordinate(0:100000,1:3)) ! 2.29 MB
       initFlag_coord=.false.
    end if

    do
       if ((firstEventNumber.le.ubound(coordinate,dim=1)).and.(firstEventNumber.ge.lbound(coordinate,dim=1))) then
          coordinate(firstEventNumber,:)=x
          exit
       else
          if (size(coordinate,dim=1).gt.10000000) then
             write(*,*) 'Vector coordinate gets too big in saveCoordinate',firstEventNumber,lbound(coordinate,dim=1),&
                  &uBound(coordinate,dim=1)
             stop
          end if
          ! Save field, enlarge it and then write saved field back to first entries:
          allocate(save_coordinate(lbound(coordinate,dim=1):uBound(coordinate,dim=1),1:3))
          do i=lbound(coordinate,dim=1),ubound(coordinate,dim=1)
             save_coordinate(i,:)=coordinate(i,:)
          end do
          deallocate(coordinate)
          allocate (coordinate(lbound(save_coordinate,dim=1):uBound(save_coordinate,dim=1)*3,1:3))
          do i=lbound(save_coordinate,dim=1),ubound(save_coordinate,dim=1)
             coordinate(i,:)=save_coordinate(i,:)
          end do
          deallocate(save_coordinate)
       end if
    end do
  end subroutine saveCoordinate

  !****************************************************************************
  !****s* initLowPhoton/getCoordinate
  ! NAME
  ! subroutine getCoordinate(firstEventNumber,x)
  ! PURPOSE
  ! Used to return production coordinates of particles
  ! INPUTS
  ! * integer :: firstEventNumber -- firstEventnumber associated to position x
  ! OUTPUT
  ! * real, dimension(1:3) ::  x -- position
  !****************************************************************************
  subroutine getCoordinate(firstEventNumber,x)

    real, dimension(1:3), intent(out) ::  x
    integer, intent(in) :: firstEventNumber

    if (.not.allocated(coordinate)) then
       write(*,*) 'ERROR in getcoordinate: Vector coordinate is not allocated'
       x=0.
       return
    end if

    if ((firstEventNumber.le.ubound(coordinate,dim=1)).and.(firstEventNumber.ge.lbound(coordinate,dim=1))) then
       x=coordinate(firstEventNumber,:)
    else
       x=0.
       write(*,*) 'WARNING : getCoordinate is out of bounds',firstEventNumber,lbound(coordinate,dim=1),ubound(coordinate,dim=1)
    end if

  end subroutine getCoordinate


  !****************************************************************************
  !****s* initLowPhoton/omegaProdInfo
  ! NAME
  ! subroutine omegaProdInfo(part, med)
  ! PURPOSE
  ! ...
  !****************************************************************************
  subroutine omegaProdInfo(part, med)
    use particleDefinition
    use mediumDefinition
    use constants, only: rhoNull
    use inputGeneral, only: current_run_number, num_Runs_sameEnergy, &
         num_Energies
    use histf90

    type(particle), intent(in) :: part
    type(medium), intent(in) :: med

    logical, save :: init = .true.
    type(histogram), save :: prodDensity

    !**************************************************************************
    !****o* initLowPhoton/omegaProdInfo.dat
    ! NAME
    ! file omegaProdInfo.dat
    ! PURPOSE
    ! This file contains informations about omega mesons at production time
    ! (event number, perweight, 4-momentum, position, density at production
    ! point, etc).
    !**************************************************************************

    open(66, file = "omegaProdInfo.dat", position = 'append')

    if (init) then
      call CreateHist(prodDensity,'density at production point [rho/rho0]',0.,1.,0.01)

      rewind(66)
      write(66,'(A)')   "### This file contains information for all omega photoproduction events:"
      write(66,'(A,A)') "### run, event, perweight, 4-momentum [GeV], coordinates (x,y,z) [fm], density at prod. point [rho/rho0]"
      write(66,*)
      init = .false.
    end if

    write(66,'(2I7,9ES15.7)') current_run_number,part%firstEvent,part%perweight,part%momentum,part%position,med%density/rhoNull
    close(66)

    !**************************************************************************
    !****o* initLowPhoton/omegaProdDensity.dat
    ! NAME
    ! file omegaProdDensity.dat
    ! PURPOSE
    ! This file contains a histogram of the density at the production point
    ! of omega mesons.
    !**************************************************************************

    call addHist (prodDensity, med%density/rhoNull, part%perweight/float(num_Runs_sameEnergy*num_Energies))
    call WriteHist (prodDensity, file='omegaProdDensity.dat')

  end subroutine omegaProdInfo



  !****************************************************************************
  !****s* initLowPhoton/twoPiProduction
  ! NAME
  ! subroutine twoPiProduction(qnuk, momNuk, posNuk, egamma, sig2Pi, sigtot,
  ! srtFree, betaToCalcFrame, mediumAtPosition, finalState, successFlag)
  ! PURPOSE
  ! Initialize a photon nucleon -> nucleon pion pion event.
  !
  ! The resulting particles are returned via the field "finalState".
  !
  ! The assigned perweight is given for each single event by
  !   sigma_Total/(number of ensembles)/(number of nucleons per ensemble).
  !
  ! INPUTS
  ! * integer :: qnuk -- charge of nucleon
  ! * real :: egamma -- energy of gamma
  ! * real, dimension(0:3) :: sig2pi -- Xsections for two-pi-production,
  !   according to gamma2pi
  ! * real, dimension(0:3) :: momNuk -- momentum of nucleon in incoming channel
  ! * real, dimension(1:3) :: posNuk -- position of nucleon in incoming channel
  ! * real :: srtFree -- sqrt(s) in vacuum prescription
  ! * real, dimension(1:3) :: betaToCalcFrame -- velocity of calc frame in cm
  !   frame, necessary to boost from CM frame to calc-frame
  ! * real :: sigtot -- total cross section, used to set perweight
  ! OUTPUT
  ! * logical :: successFlag -- whether event could be generated
  ! * type(particle), dimension(1:3) :: finalState -- final state particles
  !****************************************************************************
  subroutine twoPiProduction(qnuk, momNuk, posNuk, egamma, sig2Pi, sigtot, &
       srtFree, betaToCalcFrame, mediumAtPosition, finalState, successFlag)

    use eventGenerator_twoPi
    use nBodyPhaseSpace, only: momenta_in_3BodyPS
    use IdTable, only: nucleon, pion
    use random, only: rn
    use lorentzTrafo, only: lorentz
    use mediumDefinition
    use minkowski, only: abs4
    use energyCalc, only: energyCorrection
    use constants, only: mN, mPi

    type(particle), dimension(1:3), intent(out) :: finalState
    logical, intent(out) :: successFlag
    integer, intent(in) :: qnuk
    real, intent(in) :: egamma
    real, dimension(0:3),intent(in) :: sig2pi
    real, dimension(0:3),intent(in) :: momNuk
    real, dimension(1:3),intent(in) :: posNuk
    real, intent(in) :: srtFree
    real, dimension(1:3), intent(in) :: betaToCalcFrame
    real, intent(in) :: sigtot

    real :: x

    real :: probability_uncharged,probability_singleCharge
    !integer ::nukcharge, pionCharge_1, pionCharge_2
    integer :: rch
    real , dimension(0:3) :: pi1, pi2 ! incoming momenta
    real , dimension(0:3) :: po1, po2,po3 ! outgoing momenta
    real , dimension(1:3,1:3) :: p3       ! momenta of three particles
    real , dimension(1:3) :: betaToLRF

    type(medium),intent(in) :: mediumAtPosition

    logical, parameter :: energyCorrection_flag=.true. ! Correct the final energies for potential
    logical :: potFailure
    real :: srtS_inMed
    integer :: numTries, index
    integer, parameter :: max_numTries = 50 ! Maximal number to find an event for the energy correction procedure (i.e. including potentials)

    !**************************************************************************
    ! Choose channel
    !**************************************************************************
    integer, save :: num_failed=0, num_total=0

    successFlag=.false.

    x=rn()

    probability_uncharged=   sig2pi(3)/sig2pi(0)  ! Uncharged channel
    probability_singleCharge=sig2pi(2)/sig2pi(0)  ! Single Charged Channel

    if (equalDistribution_twoPi) then
       probability_uncharged=   1./3. ! Uncharged channel
       probability_singleCharge=1./3. ! Single Charged Channel
    end if

    finalState(1:3)%ID=(/nucleon, pion, pion /)
    finalState(1)%mass=mN
    finalState(2:3)%mass=mPi
    finalState%perturbative = .not. realRun

    if (x.lt.probability_singleCharge) then ! one charged pion + pi^0
       finalState(1)%charge=abs(qnuk-1)
       finalState(2)%charge=0
       finalState(3)%charge=qnuk-finalState(1)%charge-finalState(2)%charge
       if (equalDistribution_twoPi) then
          ! Divide out artificial probability factor (1/3) and multiply by real probability for this channel which is given by sig2pi(2)/sig2pi(0)
          finalState%perweight=sig2pi(2)/sig2pi(0) / (1./3.) * sigtot
       else
          finalState%perweight=sigTot
       end if
       if (qnuk.eq.1) then
          rCh=2
       else
          rCh=5
       end if

    else if (x.lt.(probability_singleCharge+probability_UnCharged)) then ! uncharged channel
       finalState(1)%charge=qnuk
       finalState(2)%charge=0
       finalState(3)%charge=0
       if (equalDistribution_twoPi) then
          ! Divide out artificial probability factor (1/3) and multiply by real probability for this channel which is given by sig2pi(3)/sig2pi(0)
          finalState%perweight=sig2pi(3)/sig2pi(0) / (1./3.) * sigtot
       else
          finalState%perweight=sigtot
       end if
       if (qnuk.eq.1) then
          rCh=3
       else
          rCh=6
       end if

    else ! Double charge channel gamma Nucleon -> nucleon pi^+ pi^-
       finalState(1)%charge=qnuk
       finalState(2)%charge=1
       finalState(3)%charge=-1
       if (equalDistribution_twoPi) then
          ! Divide out artificial probability factor (1/3) and multiply by real probability for this channel which is given by sig2pi(1)/sig2pi(0)
          finalState%perweight=sig2pi(1)/sig2pi(0) / (1./3.) * sigtot
       else
          finalState%perweight=sigtot
       end if
       if (qnuk.eq.1) then
          rCh=1
       else
          rCh=4
       end if
    end if

    !**************************************************************************
    !Auswuerfeln der Impulse
    !**************************************************************************
    successFlag=.false.
    do index=1,3
       finalState(index)%position=posNuk
    end do
    correctLoop:  do numTries=0,max_numTries
       if (PascalTwoPi.and.(srtfree.lt.1.55)) then
          !  rch=1: gamma p -> pi+ pi- p
          !  rch=2: gamma p -> pi+ pi0 n
          !  rch=3: gamma p -> pi0 pi0 p
          !  rch=4: gamma n -> pi+ pi- n
          !  rch=5: gamma n -> pi- pi0 p
          !  rch=6: gamma n -> pi0 pi0 n

          ! Momentum of photon
          pi1(0)=egamma
          pi1(1)=0.
          pi1(2)=0.
          pi1(3)=egamma

          ! momentum of nucleon
          pi2(1:3)=momNuk(1:3)
          pi2(0)=sqrt(mN**2+momNuk(1)**2+momNuk(2)**2+momNuk(3)**2)

          if (debugFlag) then
             write(*,*) 'vor eventgenerator'
             write(*,*) '.srtfree=',srtfree
          end if
          ! Check threshold:
          if (srtfree.le.mN+2*mPi) return
          if (abs4(pi1+pi2).le.mN+2*mPi) return

          call eventGenerator(rch,pi1,pi2,po1,po2,po3,mediumAtPosition,posNuk,energyCorrection_Flag)

          finalState(1)%momentum(1:3)=po1(1:3)
          finalState(2)%momentum(1:3)=po2(1:3)
          finalState(3)%momentum(1:3)=po3(1:3)
          finalState(1)%momentum(0)=FreeEnergy(finalState(1))
          finalState(2)%momentum(0)=FreeEnergy(finalState(2))
          finalState(3)%momentum(0)=FreeEnergy(finalState(3))



          if (debugFlag) write(*,*) 'nach eventgenerator'
          if (energyCorrection_flag) then
             srts_inMed=abs4(momnuk+(/egamma,0.,0.,egamma/))
             betaToLRF=0.
             if (debugFlag)write(*,*) 'sum momenta', finalstate(1)%momentum+finalstate(2)%momentum+finalstate(3)%momentum
             if (debugFlag)write(*,*) 'Do energyCorrection (1)!'
             call energyCorrection (srtS_inMed, betaToLRF, -betaToCalcFrame, mediumAtPosition, &
                                    finalState(1:3), successFlag, potFailure)
             if (debugFlag)write(*,*) 'sum momenta', finalstate(2)%momentum+finalstate(1)%momentum+finalstate(3)%momentum
             if (debugFlag)write(*,*) successFlag
             if (successFlag) exit correctLoop
          else
             successFlag=.true.
             exit correctLoop
          end if

       else


          p3 = momenta_in_3BodyPS (srtfree, (/mN,mPi,mPi/))

          finalState(1)%momentum(1:3)=p3(:,1)
          finalState(2)%momentum(1:3)=p3(:,2)
          finalState(3)%momentum(1:3)=p3(:,3)

          finalState(1)%momentum(0)=FreeEnergy(finalState(1))
          finalState(2)%momentum(0)=FreeEnergy(finalState(2))
          finalState(3)%momentum(0)=FreeEnergy(finalState(3))

          if (energyCorrection_flag) then
             ! Correct the energy for the potentials
             if (debugFlag)write(*,*) 'Do energyCorrection (2)!'

             srts_inMed=abs4(momnuk+(/egamma,0.,0.,egamma/))
             betaToLRF=0.
             call energyCorrection(srtS_inMed, betaToLRF, -betaToCalcFrame, mediumAtPosition, &
                                   finalState(1:3), successFlag, potFailure)
             if (successFlag) exit correctLoop
          else
             ! Just Boost to Calculation frame
             finalState(1)%momentum(0)=sqrt(mN**2+Dot_Product(p3(:,1),p3(:,1)))
             call lorentz(betatoCalcFrame  ,finalState(1)%momentum )

             finalState(2)%momentum(0)=sqrt(mPi**2+Dot_Product(p3(:,2),p3(:,2)))
             call lorentz(betatoCalcFrame  ,finalState(2)%momentum )

             finalState(3)%momentum(0)=sqrt(mPi**2+Dot_Product(p3(:,3),p3(:,3)))
             call lorentz(betatoCalcFrame ,finalState(3)%momentum )
             successFlag=.true.
             exit correctLoop
          end if
       end if
    end do correctLoop
    num_Total=num_Total+1
    if (.not. successFlag) then
       num_failed=num_failed+1
       if (Mod(num_failed,100).eq.0) then
          write(*,*) 'Warning: No success in creating final state for 2Pi events', &
               & energyCorrection_flag
          write(*,'(A,F9.3,A,2I10)') 'failures=',float(num_failed)/float(num_Total)*100., '%, ', num_failed, num_total
       end if
    end if

  end subroutine twoPiProduction

  !****************************************************************************
  !****is* initLowPhoton/binSrts
  ! NAME
  ! subroutine binSrts(srts,srtsFree,PrintFlag)
  !
  ! PURPOSE
  ! ...
  !****************************************************************************
  subroutine binSrts(srts,srtsFree,PrintFlag)
    real, intent(in) :: srtsFree, srts
    logical, intent(in) :: PrintFlag
    real, parameter :: minSrts = 1.0
    real, parameter :: maxSrts = 3.0
    integer, parameter :: binsize = 400
    integer,dimension(0:binSize,1:2),save :: srtsArray = 0
    real :: deltaS
    integer :: index
    integer :: i
    integer,save :: number =0

    number=number+1

    deltaS=(maxSrts-minSrts)/float(binSize)

    index=Min(binSize,Max(NINT((srtsFree-minSrts)/deltaS),0))
    srtsArray(index,1)=srtsArray(index,1)+1
    index=Min(binSize,Max(NINT((srts-minSrts)/deltaS),0))
    srtsArray(index,2)=srtsArray(index,2)+1


    if (PrintFlag) then
       open(11,File='Srts_photoInit.dat')
       write(11,*) '# srts respectively srtsfree [GeV], frequency srtsFree, frequency srts'
       do i=0,binSize
          write(11,'(3F10.4)') minSrts+float(i)*deltaS,float(srtsArray(i,1:2))/float(number)/deltaS
       end do
       close(11)
    end if
  end subroutine binSrts


  !****************************************************************************
  !****s* initLowPhoton/sigma_1pi_bg
  ! NAME
  ! subroutine sigma_1pi_bg(targetNuc,phi_pi,theta_pi,q,k_pi,pf_pi)
  !
  ! PURPOSE
  ! * Evaluated bg cross section for gamma N -> N pion events.
  ! * The cross section is the difference of the total pion production cross
  !   section in the vacuum and the resonance contributions in the vacuum.
  ! * The final state momenta are evaluated in the medium
  !****************************************************************************
  function sigma_1pi_bg(targetNuc,phi,theta,q,k,pf) result(sigma)

    use photonPionProduction_medium
    use IdTable, only: nucleon
    use constants, only: mN
    use output, only: line
    use histf90
    use resProd_lepton, only: dSdOmega_k_med_res

    real, dimension(-1:1)                  :: sigma

    type(particle),            intent(in)  :: targetNuc    ! real particles
    real,                      intent(in)  :: phi, theta
    real, dimension(0:3),      intent(in)  :: q
    real, dimension(-1:1,0:3), intent(out) :: k, pf


    type(particle)             :: targetNuc_free
    real, dimension(-1:1,0:3)  :: k_free, pf_free
    integer                    :: pionCharge
    logical, dimension(-1:1)   :: success
    real,dimension (-1:1)      :: sigmares
!     real, parameter               :: epsilon=10E-8
    logical :: success_res
    sigma=0.

    ! Check that target particle is nucleon:
    if (targetNuc%ID.ne.nucleon)   return


    ! Construct particle which is free and has the 3-momentum of the
    ! considered nucleon
    targetNuc_free=targetNuc
    targetNuc_free%mass=mN              ! m=m_0
    targetNuc_free%momentum(0)=freeEnergy(targetNuc_free) ! E=sqrt(p(1:3)^2+m_0^2)
    targetNuc_free%position(:)=1000.                      ! Far away the nucleus

    success=.false.
    pionChargeLoop: do pionCharge=targetNuc%charge-1,targetNuc%charge
       ! Full pion production cross section in the vacuum:
       sigma(pionCharge)=dSigmadOmega_k_med(targetNuc_free,pionCharge,phi,theta,q,k_free(pionCharge,:),pf_free(pionCharge,:)&
            & ,success(pionCharge))
       if (.not.success(pionCharge)) then
          ! Kinematics could not be established in the vacuum
          sigma(pionCharge)=0.
          return
       end if
       ! Evaluate kinematics in the medium:
       call getKinematics_photon(pionCharge,targetNuc,phi,theta,q,k(pionCharge,:),pf(pionCharge,:),success(pionCharge))
       !pf=pf_free
       !k=k_free
       if (.not.success(pionCharge)) then
          ! Kinematics could not be established in the medium
          sigma(pionCharge)=0.
       end if
    end do pionChargeLoop

    ! Resonance contribution:
    sigmaRes=dSdOmega_k_med_res(targetNuc_free,q,k_free,pf_free,success_res)
    if (.not.success_res) then
       write(*,*) 'sigmas =',sigma
       write(*,*) 'success=',success
       do pionCharge=targetNuc%charge-1,targetNuc%charge
          write(*,'(A,I3,2(A,4E10.2))') 'pionCharge', pionCharge,'k=', k_free(pionCharge,:), 'pf=', pf_free(pionCharge,:)
       end do
       write(*,*)
       write(*,*)
       write(*,*) line
    end if
    ! Subtract resonance contributions:
    sigma=sigma-sigmares

    ! Set those back to 0 zero where we couldn't establish the kinematics:
    do pionCharge=targetNuc%charge-1,targetNuc%charge
       if (.not.success(pionCharge)) sigma(pionCharge)=0.
    end do
  end function sigma_1pi_bg


  !****************************************************************************
  !****s* initLowPhoton/vecmesProduction
  ! NAME
  ! subroutine vecmesProduction(nuc, gamma, sigVM, sigTot, srts, srtFree,
  ! betaCM, mediumAtPosition, fState, flag)
  !
  ! PURPOSE
  ! photon-induced vector meson production on a nucleon (only neutral vector
  ! mesons)
  !****************************************************************************
  subroutine vecmesProduction(nuc, gamma, sigVM, sigTot, srts, srtFree, betaCM, mediumAtPosition, fState, flag)

    use random, only: rn
    use IdTable, only: nucleon,Delta,rho,omegaMeson,phi,JPsi
    use lorentzTrafo, only: lorentzCalcBeta
    use mediumDefinition, only: medium
    use dichteDefinition, only: dichte
    use densityModule, only: densityAt
    use master_2body, only: setKinematics
    use mediumModule, only: getMediumCutOff
    use RMF, only: getRMF_flag

    type(particle),intent(in):: nuc,gamma ! incoming nucleon, photon
    real,intent(in):: sigVM(0:8)          ! cross section for different channels
    real,intent(in):: sigTot              ! total cross section
    real,intent(in):: srts                ! sqrt(s)
    real,intent(in):: srtFree             ! sqrt(s) in vacuum
    real,intent(in):: betaCM(3)           ! beta for boost from CM to Lab frame
    type(medium),intent(in) :: mediumAtPosition
    type(particle),intent(out):: fState(2) ! final state particles
    logical,intent(out):: flag
    !!!!!
    type(particle):: iState(2) ! initial state part. (gamma + N)
    real:: x,betaLRF(3)
    type(dichte) :: density

    flag=.false.
    ! set up initial state
    iState=(/gamma,nuc/)
    ! set up final state
    fState%perturbative = .not. realrun
    fState%perweight=sigTot
    fState(1)%position=nuc%position
    fState(2)%position=nuc%position
    fState(1)%charge=nuc%charge
    fState(2)%charge=0
    !!! choose channel
    x=rn()*sigVM(0)
    if (x<sigVM(1)) then
      fState(1:2)%id = (/ nucleon, rho /)
    else if (x<sum(sigVM(1:2))) then
      fState(1:2)%id = (/ nucleon, omegaMeson /)
    else if (x<sum(sigVM(1:3))) then
      fState(1:2)%id = (/ nucleon, phi /)
    else if (x<sum(sigVM(1:4))) then
      fState(1:2)%id = (/ nucleon, JPsi /)
    else if (x<sum(sigVM(1:5))) then
      fState(1:2)%id = (/ Delta, rho /)
    else if (x<sum(sigVM(1:6))) then
      fState(1:2)%id = (/ Delta, omegaMeson /)
    else if (x<sum(sigVM(1:7))) then
      fState(1:2)%id = (/ Delta, phi /)
    else
      fState(1:2)%id = (/ Delta, JPsi /)
    end if
    !!! assign the masses
    density=densityAt(nuc%position)
    if (density%baryon(0).gt.getMediumCutOff()/100. .and. .not.getRMF_flag() ) then
      betaLRF = lorentzCalcBeta (density%baryon, 'vecmesProduction')
    else
      betaLRF=0.
    end if
    call setKinematics(srts,srtFree,betaLRF,-betaCM,mediumAtPosition,iState,fState,flag)
    if (gammaN_VN_formation) then
      fState(2)%in_Formation=.true.
      fState(2)%formationTime=-999 ! use old formation time concept by default
    end if
  end subroutine vecmesProduction



  !****************************************************************************
  !****s* initLowPhoton/transitionEvent
  ! NAME
  ! subroutine transitionEvent(nuc, flagOk, outPart, W, XS, pcm, beta, nout)
  !
  ! PURPOSE
  ! This routine creates inclusive FRITIOF events above E_gamma ~ 1.4 GeV,
  ! excluding the exclusive channels N+V and Delta+V, which are treated
  ! separately.
  !
  ! It is mainly used for inclusive production of vector mesons,
  ! such as V+N+pi, V+N+2pi, V+Delta+pi, V+Delta+2pi, etc.
  !
  ! Cf. Muehlich, Falter, Mosel: Inclusive omega photoproduction off nuclei
  ! Eur. Phys. J. A 20, 499-508 (2004)
  !****************************************************************************
  subroutine transitionEvent(nuc, flagOk, outPart, W, XS, pcm, beta, nout)

    use output, only: DoPR
    use random, only: rn
    use particleDefinition, only: particle,ResetNumberGuess
    use Coll_gammaN, only: DoColl_gammaN_Fr
    use IDTable, only: nucleon, pion, eta, rho, omegaMeson, phi, JPsi !, kaon, kaonBar

    type(particle),             intent(in)   :: nuc   ! incoming nucleon
    logical, intent(out)                     :: flagOk
    type(particle),dimension(:),intent(inout):: outPart  ! outgoing particles
    real,                       intent(in)   :: W
    real, dimension(0:3),       intent(in)   :: pcm
    real, dimension(1:3),       intent(in)   :: beta
    real,dimension(0:4),        intent(in)   :: XS  ! zero component is total XS
    integer, intent(out) :: nOut

    real :: pvmd
    integer :: iTyp,i,iTry,ID1,ID2,ch2
    logical :: eventOK

    flagOk = .false.

    pVMD = rn() * XS(0)

    call ResetNumberGuess()

    ! select incoming vector meson:
    if (pvmd <= XS(1)) then
       iTyp = 1
    else if (pvmd <= sum(XS(1:2))) then
       iTyp = 2
    else if (pvmd <= sum(XS(1:3))) then
       iTyp = 3
    else if (pvmd <= sum(XS(1:4))) then
       if (W.gt.4.04) then
          iTyp = 4
       else
          iTyp = 1
       end if
    else
       write(*,*) pVMD, XS
       write(*,*) 'Error in transitionEvent'
       stop
    end if

    !=== Do FRITIOF Event ===
    iTry = 0
    do
      iTry = iTry + 1
      eventOK = .false.

      if (iTry >= 100) then
        ! transitionEvent failed
        if (DoPR(3)) write(*,*) 'transitionEvent: iTry=100!  iTyp,W = ',iTyp,W,nuc%number
        !call PYLIST(2)
        !call writeFritiofCommons(9001)
        return
      end if

      call DoColl_gammaN_Fr(nuc,outPart,eventOK,W,0.,0.,pcm,beta,iTyp)
      if (.not.eventOK) cycle

      nOut = 0
      do i=1,size(outPart)
        if (outPart(i)%ID > 0) nOut = nOut+1
      end do

      ! the following exclusive events are treated separately:
      select case (nOut)

      case (2) ! two-particle final state
        ID1 = min(outPart(1)%ID,outPart(2)%ID) ! baryon
        ID2 = SUM(outPart(1:2)%ID) - ID1       ! meson
        select case (ID1)
        case (-1:0) ! invalid
          write(*,*) 'Ooops in transitionEvent (2)'
          stop
        case (1:2) ! nucleon or Delta in final state
          if (outPart(1)%ID==ID2) then
            ch2 = outPart(1)%charge
          else
            ch2 = outPart(2)%charge
          end if
          if ((singlePi .or. resonances) .and. ID1==nucleon .and. ID2==pion .and. W<2.0) then
            cycle  ! N+pi
          else if ((ID2==rho .or. ID2==omegaMeson .or. ID2==phi .or. ID2==JPsi) .and. ch2==0) then
            cycle  ! V+N,V+Delta
          end if
        end select

      case (3) ! three-particle final state
        select case (min(outPart(1)%ID,outPart(2)%ID,outPart(3)%ID))
        case (0)
          write(*,*) 'Ooops in transitionEvent (3)'
          stop
        case (1)
          ID2 = SUM(outPart(1:3)%ID) - 1
          if ((twoPi .or. resonances) .and. ID2==2*pion .and. W<2.1) cycle
          if (pi0eta .and. iD2==pion+eta .and. W<2.5)                cycle
!           if (ID2==kaon+kaonBar .and. (outPart(1)%ID==kaon .or. outPart(2)%ID==kaon .or. outPart(3)%ID==kaon)) &
!             cycle ! N+K+Kbar
        end select

      end select

      ! leave the loop:
      exit
    end do

    call setToDefault(outpart(nout+1:uBound(outpart,1)))

    ! Success: we have a real FRITIOF event
    flagOk = .true.
    if (debugFlag) write(*,*) "transitionEvent: ",nOut, outPart(1:nOut)%ID

  end subroutine transitionEvent



end module initLowPhoton
