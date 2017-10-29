!******************************************************************************
!****m* /neutrinoAnalysis
! NAME
! module neutrinoAnalysis
!
! PURPOSE
! * This module does the analysis of the output of neutrino induced processes.
!******************************************************************************
module neutrinoAnalysis

  use AnaEvent
  use histf90
  use hist2Df90
  use initNeutrino, only: max_Hist, includeHist, K2Hist, &
      get_init_namelist,OscLength,Osc,process_ID
  use AnaEventDefinition
  use CALLSTACK

  implicit none
  private

  Public:: neutrino_Analyze, cleanUp


  !****************************************************************************
  !****g* neutrinoAnalysis/detailed_diff_output
  ! SOURCE
  logical, save :: detailed_diff_output = .false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forPion
  ! SOURCE
  logical, save :: forPion =.true.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forEta
  ! SOURCE
  logical, save :: forEta =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forKaon
  ! SOURCE
  logical, save :: forKaon =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forKaonBar
  ! SOURCE
  logical, save :: forKaonBar =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

 !*****************************************************************************
  !****g* neutrinoAnalysis/forDmeson
  ! SOURCE
  logical, save :: forDmeson =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forDbar
  ! SOURCE
  logical, save :: forDbar =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forDs_plus
  ! SOURCE
  logical, save :: forDs_plus =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forDs_minus
  ! SOURCE
  logical, save :: forDs_minus =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forNucleon
  ! SOURCE
  logical, save :: forNucleon =.true.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forLambda
  ! SOURCE
  logical, save :: forLambda =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forSigmaResonance
  ! SOURCE
  logical, save :: forSigmaResonance =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forXi
  ! SOURCE
  logical, save :: forXi =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forOmegaResonance
  ! SOURCE
  logical, save :: forOmegaResonance =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/include_W_dist
  ! SOURCE
  !
  logical, save :: include_W_dist=.false.
  ! PURPOSE
  ! If .true. then the invariant mass distributions for events with 1 pion and
  ! 1 nucleon in the final state are produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/dW_Npi
  ! SOURCE
  !
  real, save  :: dW_Npi=0.02
  ! PURPOSE
  ! for dsigma/d(InvariantMass);
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Wmin_Npi
  ! SOURCE
  !
  real, save  :: Wmin_Npi=1.08
  ! PURPOSE
  ! for dsigma/d(InvariantMass);
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Wmax_Npi
  ! SOURCE
  !
  real, save  :: Wmax_Npi=1.6
  ! PURPOSE
  ! for dsigma/d(InvariantMass);
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/dW_mupi
  ! SOURCE
  !
  real, save  :: dW_mupi=0.04
  ! PURPOSE
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Wmin_mupi
  ! SOURCE
  !
  real, save  :: Wmin_mupi=0.24
  ! PURPOSE
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Wmax_mupi
  ! SOURCE
  !
  real, save  :: Wmax_mupi=1.2
  ! PURPOSE
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/dW_muN
  ! SOURCE
  !
  real, save  :: dW_muN=0.04
  ! PURPOSE
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Wmin_muN
  ! SOURCE
  !
  real, save  :: Wmin_muN=1.04
  ! PURPOSE
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Wmax_muN
  ! SOURCE
  !
  real, save  :: Wmax_muN=2.12
  ! PURPOSE
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/numin
  ! SOURCE
  !
  real, save  :: numin=0.
  ! PURPOSE
  ! for calorimetric analysis: values for transferred energy;
  ! only work if calorimetric_analysis is .true.
  ! set the min, max and bins for nu distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/numax
  ! SOURCE
  !
  real, save  :: numax=10.0
  ! PURPOSE
  ! for calorimetric analysis: values for transferred energy;
  ! only work if calorimetric_analysis is .true.
  ! set the min, max and bins for nu distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/nubin
  ! SOURCE
  !
  real, save  :: nubin=0.1
  ! PURPOSE
  ! for calorimetric analysis: values for transferred energy;
  ! only work if calorimetric_analysis is .true.
  ! set the min, max and bins for nu distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Enumin
  ! SOURCE
  !
  real, save  :: Enumin=0.
  ! PURPOSE
  ! for calorimetric analysis: values for neutrino energy;
  ! only work if calorimetric_analysis is .true.
  ! set the min, max and bins for nu distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Enumax
  ! SOURCE
  !
  real, save  :: Enumax=10.0
  ! PURPOSE
  ! for calorimetric analysis: values for neutrino energy;
  ! only work if calorimetric_analysis is .true.
  ! set the min, max and bins for nu distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Enubin
  ! SOURCE
  !
  real, save  :: Enubin=0.1
  ! PURPOSE
  ! for calorimetric analysis: values for neutrino energy;
  ! only work if calorimetric_analysis is .true.
  ! set the min, max and bins for nu distributions
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/EkinMin
  ! SOURCE
  real, save :: EkinMin=0.
  ! PURPOSE
  ! if detailed_diff_output is TRUE:
  ! Maximal kinetic energy for dsigma/dEkin for hadrons
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/EkinMax
  ! SOURCE
  real, save :: EkinMax=2.
  ! PURPOSE
  ! if detailed_diff_output is TRUE:
  ! Maximal kinetic energy for dsigma/dEkin for hadrons
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/dEkin
  ! SOURCE
  real, save :: dEkin=0.01
  ! PURPOSE
  ! if detailed_diff_output is TRUE:
  ! Delta(eKin) for dsigma/dEKin  for hadrons
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/EkinMin_lepton
  ! SOURCE
  real, save :: EkinMin_lepton=0.
  ! PURPOSE
  ! if detailed_diff_output is TRUE:
  ! Maximal kinetic energy for dsigma/dEkin for leptons
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/EkinMax_lepton
  ! SOURCE
  real, save :: EkinMax_lepton=2.
  ! PURPOSE
  !if detailed_diff_output is TRUE:
  ! Maximal kinetic energy for dsigma/dEkin for leptons
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/dEkin_lepton
  ! SOURCE
  real, save :: dEkin_lepton=0.01
  ! PURPOSE
  ! if detailed_diff_output is TRUE:
  ! Delta(eKin) for dsigma/dEKin  for leptons
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/cost_min
  ! SOURCE
  real, save :: cost_min= -1.
  ! PURPOSE
  !if detailed_diff_output is TRUE:
  ! Minimal cos(theta) of outgoing leptons
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/cost_max
  ! SOURCE
  real, save :: cost_max= +1
  ! PURPOSE
  !if detailed_diff_output is TRUE:
  ! Maximal cos(theta) of outgoing leptons
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/delta_cost
  ! SOURCE
  real, save :: delta_cost = 0.02
  ! PURPOSE
  !if detailed_diff_output is TRUE:
  ! stepsize of cos(theta) of outgoing leptons
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/kineticEnergyDetectionThreshold_nucleon
  ! SOURCE
  !
  real, save :: kineticEnergyDetectionThreshold_nucleon=0.0
  ! PURPOSE
  ! kineticEnergyDetectionThreshold
  ! lower detection threshold for nucleon kinetic energies
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/kineticEnergyDetectionThreshold_chargedpion
  ! SOURCE
  !
  real, save :: kineticEnergyDetectionThreshold_chargedpion=0.0
  ! PURPOSE
  ! kineticEnergyDetectionThreshold
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/kineticEnergyDetectionThreshold_neutralpion
  ! SOURCE
  !
  real, save :: kineticEnergyDetectionThreshold_neutralpion=0.0
  ! PURPOSE
  ! kineticEnergyDetectionThreshold
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/kineticEnergyDetectionThreshold_lepton
  ! SOURCE
  !
  real, save :: kineticEnergyDetectionThreshold_lepton=0.0
  ! PURPOSE
  ! kineticEnergyDetectionThreshold
  ! only lepton kinetic energies above this threshold can be detected
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/AngleUpperDetectionThresholdDegrees_nucleon
  ! SOURCE
  !
  real, save :: AngleUpperDetectionThresholdDegrees_nucleon=180.0
  ! PURPOSE
  ! nucleon angles up to this value can be detected
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/AngleUpperDetectionThresholdDegrees_chargedpion
  ! SOURCE
  !
  real, save :: AngleUpperDetectionThresholdDegrees_chargedpion=180.0
  ! PURPOSE
  ! charged pion angles up to this value can be detected
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/AngleUpperDetectionThresholdDegrees_neutralpion
  ! SOURCE
  !
  real, save :: AngleUpperDetectionThresholdDegrees_neutralpion=180.0
  ! PURPOSE
  ! neutral pion angles angles up to this value can be detected
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/AngleUpperDetectionThresholdDegrees_lepton
  ! SOURCE
  !
  real, save :: AngleUpperDetectionThresholdDegrees_lepton=180.0
  ! PURPOSE
  ! lepton angles up to this value can be detected
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/radialScale
  ! SOURCE
  !
  real, save :: radialScale=0.0
  ! PURPOSE
  ! If radial position of nucleon < radialScale, then the nucleon is assumed
  ! to be bound
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/inclusiveAnalysis
  ! SOURCE
  !
  logical, save ::  inclusiveAnalysis=.false.
  ! PURPOSE
  ! If .true. then we don't care whether particle has left the nucleus or not
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Fissum_analysis
  ! SOURCE
  logical, save ::  Fissum_analysis=.false.
  ! PURPOSE
  ! do analysis with cuts as needed for Fig 25 in
  ! Fissum et al, PRC 70, 034606 (2004)
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/ZeroPion_analysis
  ! SOURCE
  logical, save ::  ZeroPion_analysis=.false.
  ! PURPOSE
  ! produce output of xsec for various final states with 0 pions and 2 pions
  ! see file see sigma_0pions.dat  for the list of the final states
  !
  ! see files neutrino_0pions.dat,  neutrino_0pions_QE.dat,
  ! neutrino_0pions_Delta.dat, ... for output
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/calorimetric_analysis
  ! SOURCE
  logical, save ::  calorimetric_analysis=.false.
  ! PURPOSE
  ! do calorimetric energy-transfer and neutrino-energy reconstruction
  ! (for each QE, Delta, ...)  as in the MINOS experiment
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/reconstruct_neutrino_energy
  ! SOURCE
  logical, save ::  reconstruct_neutrino_energy=.false.
  ! PURPOSE
  ! reconstruct neutrino energy for final state in "specificEvent_analysis"
  ! NOTES
  ! .true. must be combined with specificEvent_analysis=.true. and
  ! at least one specific event .true.
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/specificEvent_analysis
  ! SOURCE
  logical, save ::  specificEvent_analysis=.false.
  ! PURPOSE
  ! do analysis for specific final states
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/no_pi
  ! SOURCE
  logical, save ::  no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=1,
  ! no_pi (for example, for QE-like MiniBooNE)
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/p_Xn_no_pi
  ! SOURCE
  logical, save ::  p_Xn_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=2
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/piplus
  ! SOURCE
  logical, save ::  piplus=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=3, 1 pi+ X nucleons
  ! mesons of other flavor
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/piplus_MULTI
  ! SOURCE
  logical, save ::  piplus_MULTI=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=4
  ! >=1 pi+  X other pions (incl pi+) X nucleons
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pi0
  ! SOURCE
  logical, save ::  pi0=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=5,
  ! 1 pi0 X nucleons, plus mesons of other flavor
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pi0_MULTI
  ! SOURCE
  logical, save ::  pi0_MULTI=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=6,
  ! >=1 pi0  X other pions X nucleons, (pi0 K2K)
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/piminus
  ! SOURCE
  logical, save ::  piminus=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=7
  ! 1 pi-  X other pions X nucleons
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

   !***************************************************************************
  !****g* neutrinoAnalysis/piminus_MULTI
  ! SOURCE
  logical, save ::  piminus_MULTI=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=8
  ! >=1 pi-  X other pions X nucleons
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pp_no_pi
  ! SOURCE
  logical, save ::  pp_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states:  specificEvent=9
  ! 2 protons, X neutrons, 0 pions
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pn_no_pi
  ! SOURCE
  logical, save ::  pn_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=10
  ! 1 neutron, 1 proton, 0 pions
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/nn_no_pi
  ! SOURCE
  logical, save ::  nn_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=11
  ! 2 neutrons, X protons, 0 pions
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pp_Xn_no_pi
  ! SOURCE
  logical, save ::  pp_Xn_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=12
  ! 2 protons, X neutrons, 0 pions
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/nn_Xp_no_pi
  ! SOURCE
  logical, save ::  nn_Xp_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=13
  ! 2 neutrons, X protons, 0 pions
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/ppp_Xn_no_pi
  ! SOURCE
  logical, save ::  ppp_Xn_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=14
  ! 3 protons, X neutrons, 0 pions
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pppp_Xn_no_pi
  ! SOURCE
  logical, save ::  pppp_Xn_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=15
  ! 4 protons, X neutrons, 0 pions
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/p_no_pi
  ! SOURCE
  logical, save ::  p_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=16
  ! 1 proton, 0 neutron, 0 pion
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/n_no_pi
  ! SOURCE
  logical, save ::  n_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=17
  ! 1 neutron, 0 proton, 0 pion
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Xn_no_pi
  ! SOURCE
  logical, save ::  Xn_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=18,
  ! 0 proton, X neutrons, 0 pions
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/excl_hadron
  ! SOURCE
  logical, save :: excl_hadron=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=19,20,21
  ! exclusive 1 pion, no other pions or other mesons of different flavor
  ! There could be still other mesons which are heavier than the D,
  ! Such events (very rare at DUNE energies) could be counted as exclusive
  ! single-meson cross section.
  ! This could be cured by extending the list of stable mesons
  !
  ! value can be changed in the namelist nl_specificEvent
  !****************************************************************************

  logical, save :: excl_pi0=.false.
  logical, save :: excl_piplus=.false.
  logical, save :: excl_piminus=.false.


  !****************************************************************************
  !****g* neutrinoAnalysis/QEp
  ! SOURCE
  logical ::  QEp=.false.
  ! PURPOSE
  ! if .true,
  ! do analysis for specific analysis for QE-like event with 1 mu, 0 pi, X p
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pcut
  ! SOURCE
  real, save ::  pcut = 0.0
  ! PURPOSE
  ! determines lower acceptance cut for outgoing protons
  !
  ! values can be changed in the namelist neutrinoAnalysis via the variable
  ! kineticEnergyDetectionThreshold_nucleon
  !****************************************************************************

  integer, parameter ::  max_SpeEvent=21
  ! maximum number of special Events for which detailed analysis is being done
  ! special events in file  'includeSpeEvent'

  !****************************************************************************
  !****g* neutrinoAnalysis/binsizeQ2
  ! SOURCE
  real, save    ::  binsizeQ2=0.01
  ! PURPOSE
  ! do analysis for specific final states:
  ! binning for reconstruction of Q2 and Enu
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/binsizeEnu
  ! SOURCE
  real, save    ::  binsizeEnu=0.02
  ! PURPOSE
  ! do analysis for specific final states:
  ! binning for reconstruction of Q2 and Enu
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/maxQ2
  ! SOURCE
  real, save    ::  maxQ2=5.0
  ! PURPOSE
  ! do analysis for specific final states:
  ! max values for reconstruction of Q2 and Enu
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/maxEnu
  ! SOURCE
  real, save    ::  maxEnu=5.0
  ! PURPOSE
  ! do analysis for specific final states:
  ! max values for reconstruction of Q2 and Enu
  !
  ! values can be changed in the namelist nl_specificEvent
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/outputEvents
  ! SOURCE
  !
  logical, save  :: outputEvents = .false.
  ! PURPOSE
  ! If .true. then all events are printed to the file 'FinalEvents.dat'.
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Xsection_analysis
  ! SOURCE
  !
  logical, save  :: Xsection_analysis = .false.
  ! PURPOSE
  ! If .true. then files "..._total_Xsection_..." are printed.
  !****************************************************************************



  logical, save :: initflag=.true.


  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists,dE_hists_QE,dE_hists_Delta, &
       & dE_hists_highRes,dE_hists_gen0,dE_hists_1piBG,dE_hists_2piBG,dE_hists_DIS,dE_hists_2p2h
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists_Multi,dE_hists_QE_Multi, &
       & dE_hists_Delta_Multi,dE_hists_highRes_Multi,dE_hists_gen0_Multi,dE_hists_1piBG_Multi,&
       & dE_hists_2piBG_Multi,dE_hists_DIS_Multi,dE_hists_2p2h_Multi
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists_1X,dE_hists_QE_1X, &
       & dE_hists_Delta_1X,dE_hists_highRes_1X,dE_hists_gen0_1X,dE_hists_1piBG_1X, &
       & dE_hists_2piBG_1X,dE_hists_DIS_1X, dE_hists_2p2h_1X
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists_2X,dE_hists_QE_2X, &
       & dE_hists_Delta_2X,dE_hists_highRes_2X,dE_hists_gen0_2X,dE_hists_1piBG_2X, &
       & dE_hists_2piBG_2X,dE_hists_DIS_2X,dE_hists_2p2h_2X

  ! used for ZeroPion_analysis=.true.  for kinetic energy distributions of nucleons
  ! in events with 0 pions
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists_0pions, dE_hists_QE_0pions, &
       & dE_hists_Delta_0pions,dE_hists_highRes_0pions, dE_hists_2p2h_0pions, &
       & dE_hists_1piBG_0pions, dE_hists_2piBG_0pions, dE_hists_DIS_0pions
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists_Multi_0pions, &
       & dE_hists_QE_Multi_0pions,dE_hists_Delta_Multi_0pions, dE_hists_highRes_Multi_0pions, &
       & dE_hists_2p2h_Multi_0pions,dE_hists_1piBG_Multi_0pions, dE_hists_2piBG_Multi_0pions, &
       & dE_hists_DIS_Multi_0pions
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists_1X_0pions,&
       & dE_hists_QE_1X_0pions, dE_hists_Delta_1X_0pions, &
       & dE_hists_highRes_1X_0pions, dE_hists_2p2h_1X_0pions, dE_hists_1piBG_1x_0pions, &
       & dE_hists_2piBG_1x_0pions,dE_hists_DIS_1x_0pions
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists_2X_0pions, &
       & dE_hists_QE_2X_0pions, dE_hists_Delta_2X_0pions, &
       & dE_hists_highRes_2X_0pions, dE_hists_2p2h_2X_0pions, dE_hists_1piBG_2x_0pions, &
       & dE_hists_2piBG_2x_0pions,dE_hists_DIS_2x_0pions


  ! used for include_W_dist=.true.
  type(histogram),save,dimension(-1:1) :: dW_Npi_hists, dW_mupi_hists, dW_muN_hists

  ! used for calorimetric_analysis=.true.
  type(histogram2D), dimension(0:max_Hist), save :: Ehad_versusNu_hist, Enurestored_versusEnu
  type(histogram), dimension(0:max_Hist), save :: dSigdNu_hist, dSigdEhad_hist, dSigdEnu,&
      & dSigdEnurestored


  ! used for specificEvent_analysis=.true.
  type(histogram),save, dimension(1:max_SpeEvent,0:max_Hist) :: dEnu_hist, dElepton_hist, &
      & dcoslepton_hist, dQ2lepton_hist, dQ2plepton_hist

  type(histogram2D), save, dimension(1:max_SpeEvent,0:max_Hist) :: dSigmaMC_EprimeCost_0pi


  type(histogram),save,dimension(1:numStableParts,-2:2) :: dTheta_hists,  &
       & dPhi_hists, dTheta_hists_QE,dPhi_hists_QE, &
       & dTheta_hists_Delta,dPhi_hists_Delta,dTheta_hists_highRes,  &
       & dPhi_hists_highRes, dTheta_hists_gen0,dPhi_hists_gen0, &
       & dTheta_hists_1piBG,dPhi_hists_1piBG,dTheta_hists_2piBG,dPhi_hists_2piBG

  type(histogram2D),save,dimension(1:numStableParts,-2:2) :: dOmega_hists, &
       & dOmega_hists_QE,dOmega_hists_Delta, &
       & dOmega_hists_highRes,dOmega_hists_gen0,dOmega_hists_1piBG,  &
       & dOmega_hists_2piBG

  type(histogram),save,dimension(1:numStableParts,-2:2) :: dEcostheta_hists, &
       & dEcostheta_hists_QE,dEcostheta_hists_Delta, &
       & dEcostheta_hists_highRes, dEcostheta_hists_gen0,dEcostheta_hists_1piBG, &
       & dEcostheta_hists_2piBG,dEcostheta_hists_DIS
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dEcostheta_hists_MULTI, &
       & dEcostheta_hists_QE_MULTI, &
       & dEcostheta_hists_Delta_MULTI, dEcostheta_hists_highRes_MULTI, &
       & dEcostheta_hists_gen0_MULTI, &
       & dEcostheta_hists_1piBG_MULTI,dEcostheta_hists_2piBG_MULTI,  &
       & dEcostheta_hists_DIS_MULTI
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dcostheta_hists, &
       & dcostheta_hists_QE, dcostheta_hists_Delta, &
       & dcostheta_hists_highRes, dcostheta_hists_gen0,dcostheta_hists_1piBG, &
       & dcostheta_hists_2piBG,dcostheta_hists_DIS
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dcostheta_hists_MULTI, &
       & dcostheta_hists_QE_MULTI, &
       & dcostheta_hists_Delta_MULTI, dcostheta_hists_highRes_MULTI, &
       & dcostheta_hists_gen0_MULTI,&
       & dcostheta_hists_1piBG_MULTI, dcostheta_hists_2piBG_MULTI, &
       & dcostheta_hists_DIS_MULTI

  type(histogram),save :: hNucleonVacuumMass

  logical :: initHists
  logical, dimension(1:max_speEvent), save :: includeSpeEvent




contains





  !****************************************************************************
  !****s* neutrinoAnalysis/readinput
  ! NAME
  ! subroutine readinput
  ! INPUTS
  ! NONE
  ! OUTPUT
  ! NONE
  ! PURPOSE
  ! This subroutine reads the namelist "neutrinoAnalysis".
  ! Only called once to initialize the module.
  !****************************************************************************
  subroutine readinput
    use output

    integer :: IOS

    !**************************************************************************
    !****n* neutrinoAnalysis/NeutrinoAnalysis
    ! NAME
    ! NAMELIST /NeutrinoAnalysis/
    ! PURPOSE
    ! This namelist includes:
    ! * detailed_diff_output
    ! * include_W_dist
    ! * kineticEnergyDetectionThreshold_lepton
    ! * AngleUpperDetectionThresholdDegrees_lepton
    ! * kineticEnergyDetectionThreshold_nucleon
    ! * AngleUpperDetectionThresholdDegrees_nucleon
    ! * kineticEnergyDetectionThreshold_chargedpion
    ! * AngleUpperDetectionThresholdDegrees_chargedpion
    ! * kineticEnergyDetectionThreshold_neutralpion
    ! * AngleUpperDetectionThresholdDegrees_neutralpion
    ! * inclusiveAnalysis
    ! * Fissum_analysis
    ! * ZeroPion_analysis
    ! * calorimetric_analysis
    ! * radialScale
    ! * reconstruct_neutrino_energy
    ! * outputEvents
    ! * specificEvent_analysis
    ! * Xsection_analysis
    !**************************************************************************
    NAMELIST /neutrinoAnalysis/ &
         detailed_diff_output, include_W_dist, &
         kineticEnergyDetectionThreshold_lepton, &
         AngleUpperDetectionThresholdDegrees_lepton, &
         kineticEnergyDetectionThreshold_nucleon, &
         AngleUpperDetectionThresholdDegrees_nucleon, &
         kineticEnergyDetectionThreshold_chargedpion, &
         AngleUpperDetectionThresholdDegrees_chargedpion,&
         kineticEnergyDetectionThreshold_neutralpion, &
         AngleUpperDetectionThresholdDegrees_neutralpion,&
         inclusiveAnalysis, Fissum_analysis,            &
         ZeroPion_analysis, calorimetric_analysis, &
         radialScale, reconstruct_neutrino_energy, &
         outputEvents, specificEvent_analysis, &
         Xsection_analysis

    !**************************************************************************
    !****n* neutrinoAnalysis/W_distributions
    ! NAME
    ! NAMELIST /W_distributions/
    ! PURPOSE
    ! This Namelist includes:
    ! * dW_Npi
    ! * Wmin_Npi
    ! * Wmax_Npi
    ! * dW_mupi
    ! * Wmin_mupi
    ! * Wmax_mupi
    ! * dW_muN
    ! * Wmin_muN
    ! * Wmax_muN
    !**************************************************************************
    NAMELIST /W_distributions/  dW_Npi,  Wmin_Npi,  Wmax_Npi, &
             dW_mupi, Wmin_mupi, Wmax_mupi, &
             dW_muN,  Wmin_muN,  Wmax_muN

    !**************************************************************************
    !****n* neutrinoAnalysis/nl_calorimetric_analysis
    ! NAME
    ! NAMELIST /nl_calorimetric_analysis/
    ! PURPOSE
    ! This Namelist includes:
    ! * numin
    ! * numax
    ! * nubin
    ! * Enumin
    ! * Enumax
    ! * Enubin
    !**************************************************************************
    NAMELIST /nl_calorimetric_analysis/ numin,numax,nubin,Enumin,Enumax,Enubin

    !**************************************************************************
    !****n* neutrinoAnalysis/nl_specificEvent
    ! NAME
    ! NAMELIST /nl_specificEvent/
    ! PURPOSE
    ! This namelist includes:
    ! * no_pi
    ! * p_Xn_no_pi
    ! * piplus
    ! * piplus_MULTI
    ! * pi0
    ! * pi0_MULTI
    ! * piminus
    ! * piminus_MULTI
    ! * pp_no_pi
    ! * pn_no_pi
    ! * nn_no_pi
    ! * pp_Xn_no_pi
    ! * nn_Xp_no_pi
    ! * ppp_Xn_no_pi
    ! * pppp_Xn_no_pi
    ! * p_no_pi
    ! * n_no_pi
    ! * Xn_no_pi
    ! * binsizeQ2
    ! * binsizeEnu
    ! * maxQ2
    ! * maxEnu
    ! * excl_hadron
    ! * QEp
    !**************************************************************************
    NAMELIST /nl_specificEvent/ no_pi, p_Xn_no_pi, piplus, piplus_MULTI, &
         pi0, pi0_MULTI, piminus, piminus_MULTI, &
         pp_no_pi, pn_no_pi, nn_no_pi, pp_Xn_no_pi, nn_Xp_no_pi, &
         ppp_Xn_no_pi, pppp_Xn_no_pi, p_no_pi, n_no_pi, Xn_no_pi, &
         binsizeQ2, binsizeEnu, maxQ2, maxEnu,excl_hadron, QEp

    !**************************************************************************
    !****n* neutrinoAnalysis/detailed_diff
    ! NAME
    ! NAMELIST /detailed_diff/
    ! PURPOSE
    ! This namelist includes:
    ! * EkinMin
    ! * EkinMax
    ! * dEkin
    ! * EkinMin_lepton
    ! * EkinMax_lepton
    ! * dEkin_lepton
    ! * forPion
    ! * forEta
    ! * forKaon
    ! * forKaonBar
    ! * forDmeson
    ! * forDbar
    ! * forDs_plus
    ! * forDs_minus
    ! * forNucleon
    ! * forLambda
    ! * forSigmaResonance
    ! * forXi
    ! * forOmegaResonance
    !**************************************************************************
    NAMELIST /detailed_diff/ ekinMin,ekinMax, dEkin,                         &
         &  ekinMin_lepton, ekinMax_lepton, dEkin_lepton,                    &
         &  cost_min,cost_max,delta_cost,                                    &
         &  forpion, foreta, forkaon, forkaonBar, forDmeson, forDbar,        &
         &  forDs_plus, forDs_minus,                                         &
         &  fornucleon, forLambda, forSigmaResonance,forXi,forOmegaResonance

    call Write_ReadingInput('neutrinoAnalysis',0)
    rewind(5)
    read(5,nml=neutrinoAnalysis,IOSTAT=IOS)
    call Write_ReadingInput('neutrinoAnalysis',0,IOS)
    write(*,*) '  detailed_diff_output  =', detailed_diff_output
    write(*,*) '  include_W_dist        =', include_W_dist
    write(*,*) '  calorimetric_analysis =', calorimetric_analysis
    write(*,*) '  specificEvent_analysis=', specificEvent_analysis
    write(*,*) '  reconstruct_neutrino_energy=', reconstruct_neutrino_energy
    write(*,nml=neutrinoAnalysis)
    call Write_ReadingInput('neutrinoAnalysis',1)

    if (detailed_diff_output) then
       call Write_ReadingInput('detailed_diff',0)
       rewind(5)
       read(5,nml=detailed_diff,IOSTAT=IOS)
       call Write_ReadingInput('detailed_diff',0,IOS)
       write(*,'(3(A,g14.5))') '  ekinMin= ',ekinMin, &
            &'  ekinMax= ',ekinMax, &
            &'  dEkin= ',dEkin, &
            &'  ekinmin_lepton=  ',ekinmin_lepton,&
            &'  ekinmax_lepton=  ',ekinmax_lepton,&
            &'  dEkin_lepton=  ', dEkin_lepton,&
            &'  cost_min=  ', cost_min, &
            &'  cost_max=  ', cost_max, &
            &'  delta_cost=  ', delta_cost
       write(*,'(10(A,L3,/))') &
            &'    pion=    ',forpion, &
            &'    Nucleon= ',fornucleon, &
            &'    Lambda=  ',forLambda, &
            &'    Sigma=   ',forSigmaResonance, &
            &'    Xi=      ',forXi, &
            &'    Omega    ',forOmegaResonance, &
            &'    eta=     ',foreta,&
            &'    kaon=    ',forkaon,&
            &'    kaonBar= ',forkaonBar,&
            &'    Dmeson=  ',forDmeson,&
            &'    Dbar=    ',forDbar,&
            &'    Ds+=     ',forDs_plus,&
            &'    Ds-=     ',forDs_minus

       call Write_ReadingInput('detailed_diff',1)
    end if

    if (include_W_dist) then
       call Write_ReadingInput('W_distributions',0)
       rewind(5)
       read(5,nml=W_distributions,IOSTAT=IOS)
       call Write_ReadingInput('W_distributions',0,IOS)
       write(*,'(3(A,g14.5,/))') &
            & '  dW_Npi = ',dW_Npi,&
            & '  dW_mupi= ',dW_mupi, &
            & '  dW_muN = ',dW_muN
       call Write_ReadingInput('W_distributions',1)
    end if

    if (calorimetric_analysis) then
       call Write_ReadingInput('nl_calorimetric_analysis',0)
       rewind(5)
       read(5,nml=nl_calorimetric_analysis,IOSTAT=IOS)
       call Write_ReadingInput('nl_calorimetric_analysis',0,IOS)
       write(*,'(3(A,g14.5,/))') &
            & '   numin= ',numin, &
            & '   numax= ',numax, &
            & '   nubin= ',nubin
       write(*,'(3(A,g14.5,/))') &
            & '   Enumin= ',Enumin, &
            & '   Enumax= ',Enumax, &
            & '   Enubin= ',Enubin
       call Write_ReadingInput('nl_calorimetric_analysis',1)
    end if

    if (specificEvent_analysis) then
       call Write_ReadingInput('nl_specificEvent',0)
       rewind(5)
       read(5,nml=nl_specificEvent,IOSTAT=IOS)
       call Write_ReadingInput('nl_specificEvent',0,IOS)
       write(*,'(A,/,20(A,L3,/))') '  specificEvent_analysis = .true. :', &
            & '    no_pi=        ', no_pi,&
            & '    p_Xn_no_pi=   ', p_Xn_no_pi, &
            & '    piplus=       ', piplus,&
            & '    piplus_MULTI= ', piplus_MULTI, &
            & '    pi0=          ', pi0, &
            & '    pi0_MULTI=    ', pi0_MULTI, &
            & '    piminus=      ', piminus, &
            & '    piminus_MULTI=', piminus_MULTI, &
            & '    pp_no_pi=     ', pp_no_pi, &
            & '    pn_no_pi=     ', pn_no_pi, &
            & '    nn_no_pi=     ', nn_no_pi, &
            & '    pp_Xn_no_pi=  ', pp_Xn_no_pi, &
            & '    nn_Xp_no_pi=  ', nn_Xp_no_pi, &
            & '    ppp_Xn_no_pi= ', ppp_Xn_no_pi, &
            & '    pppp_Xn_no_pi=', pppp_Xn_no_pi, &
            & '    p_no_pi=      ', p_no_pi, &
            & '    n_no_pi=      ', n_no_pi, &
            & '    Xn_no_pi=     ', Xn_no_pi, &
            & '    excl_hadron   ', excl_hadron, &
            & '    QEp           ', QEp
        if(QEp) then
          pcut = sqrt(kineticEnergyDetectionThreshold_nucleon**2  &
               &  + 2*0.938272*kineticEnergyDetectionThreshold_nucleon)
          write (*,*) 'pcut = ', pcut
        end if

 ! The switch excl_hadron selects truly exclusive 1-meson events
 ! ie. exactly 1 meson with fixed charged and no other mesons of any kind
 ! switch QEp selects QE-like events with 1 mu, 0 pi, X p


       call Write_ReadingInput('nl_specificEvent',1)
    end if

    call set_Exclusive(excl_hadron)
    if (excl_hadron .eqv. .true.) then
       excl_pi0 = .true.
       excl_piminus = .true.
       excl_piplus = .true.
       piplus = .false.
       pi0 = .false.
       piminus = .false.
    end if

    write(*,*) 'neutrinoAnalysis :',QEp,pcut
    call set_QElike(QEp,pcut)

    includeSpeEvent = (/no_pi,p_Xn_no_pi,piplus,piplus_MULTI,pi0,pi0_MULTI, &
         & piminus, piminus_MULTI,pp_no_pi, pn_no_pi, nn_no_pi, pp_Xn_no_pi,&
         & nn_Xp_no_pi, ppp_Xn_no_pi, pppp_Xn_no_pi,p_no_pi, n_no_pi, Xn_no_pi,&
         & excl_pi0,excl_piplus,excl_piminus/)

  end subroutine readinput





  subroutine cleanUp
    integer :: i,j
    do i=1,numStableParts
       do j=-2,2
          call RemoveHist(dE_hists(i,j))
          call RemoveHist(dE_hists_QE(i,j))
          call RemoveHist(dE_hists_Delta(i,j))
          call RemoveHist(dE_hists_highRes(i,j))
          call RemoveHist(dE_hists_1piBG(i,j))
          call RemoveHist(dE_hists_2piBG(i,j))
          call RemoveHist(dE_hists_DIS(i,j))
          call RemoveHist(dE_hists_gen0(i,j))
          call RemoveHist(dE_hists_Multi(i,j))
          call RemoveHist(dE_hists_QE_Multi(i,j))
          call RemoveHist(dE_hists_Delta_Multi(i,j))
          call RemoveHist(dE_hists_highRes_Multi(i,j))
          call RemoveHist(dE_hists_1piBG_Multi(i,j))
          call RemoveHist(dE_hists_DIS_Multi(i,j))
          call RemoveHist(dE_hists_gen0_Multi(i,j))
          !
          call RemoveHist(dE_hists_1X(i,j))
          call RemoveHist(dE_hists_QE_1X(i,j))
          call RemoveHist(dE_hists_Delta_1X(i,j))
          call RemoveHist(dE_hists_highRes_1X(i,j))
          call RemoveHist(dE_hists_1piBG_1X(i,j))
          call RemoveHist(dE_hists_2piBG_1X(i,j))
          call RemoveHist(dE_hists_DIS_1X(i,j))
          call RemoveHist(dE_hists_gen0_1X(i,j))
          call RemoveHist(dE_hists_2X(i,j))
          call RemoveHist(dE_hists_QE_2X(i,j))
          call RemoveHist(dE_hists_Delta_2X(i,j))
          call RemoveHist(dE_hists_highRes_2X(i,j))
          call RemoveHist(dE_hists_1piBG_2X(i,j))
          call RemoveHist(dE_hists_DIS_2X(i,j))
          call RemoveHist(dE_hists_gen0_2X(i,j))
          !
          call RemoveHist(dTheta_hists(i,j))
          call RemoveHist(dPhi_hists(i,j))
          call RemoveHist(dTheta_hists_QE(i,j))
          call RemoveHist(dPhi_hists_QE(i,j))
          call RemoveHist(dTheta_hists_Delta(i,j))
          call RemoveHist(dPhi_hists_Delta(i,j))
          call RemoveHist(dTheta_hists_highRes(i,j))
          call RemoveHist(dPhi_hists_highRes(i,j))
          call RemoveHist(dTheta_hists_1piBG(i,j))
          call RemoveHist(dPhi_hists_1piBG(i,j))
          call RemoveHist(dTheta_hists_2piBG(i,j))
          call RemoveHist(dPhi_hists_2piBG(i,j))
          call RemoveHist(dTheta_hists_gen0(i,j))
          call RemoveHist(dPhi_hists_gen0(i,j))
          call RemoveHist2d(dOmega_hists(i,j))
          call RemoveHist2d(dOmega_hists_QE(i,j))
          call RemoveHist2d(dOmega_hists_Delta(i,j))
          call RemoveHist2d(dOmega_hists_highRes(i,j))
          call RemoveHist2d(dOmega_hists_1piBG(i,j))
          call RemoveHist2d(dOmega_hists_gen0(i,j))
          !
          call RemoveHist(dEcostheta_hists(i,j))
          call RemoveHist(dEcostheta_hists_QE(i,j))
          call RemoveHist(dEcostheta_hists_Delta(i,j))
          call RemoveHist(dEcostheta_hists_highRes(i,j))
          call RemoveHist(dEcostheta_hists_1piBG(i,j))
          call RemoveHist(dEcostheta_hists_2piBG(i,j))
          call RemoveHist(dEcostheta_hists_DIS(i,j))
          call RemoveHist(dEcostheta_hists_gen0(i,j))
          call RemoveHist(dEcostheta_hists_MULTI(i,j))
          call RemoveHist(dEcostheta_hists_QE_MULTI(i,j))
          call RemoveHist(dEcostheta_hists_Delta_MULTI(i,j))
          call RemoveHist(dEcostheta_hists_highRes_MULTI(i,j))
          call RemoveHist(dEcostheta_hists_1piBG_MULTI(i,j))
          call RemoveHist(dEcostheta_hists_DIS_MULTI(i,j))
          call RemoveHist(dEcostheta_hists_gen0_MULTI(i,j))
          call RemoveHist(dcostheta_hists(i,j))
          call RemoveHist(dcostheta_hists_QE(i,j))
          call RemoveHist(dcostheta_hists_Delta(i,j))
          call RemoveHist(dcostheta_hists_highRes(i,j))
          call RemoveHist(dcostheta_hists_1piBG(i,j))
          call RemoveHist(dcostheta_hists_2piBG(i,j))
          call RemoveHist(dcostheta_hists_DIS(i,j))
          call RemoveHist(dcostheta_hists_gen0(i,j))
          call RemoveHist(dcostheta_hists_MULTI(i,j))
          call RemoveHist(dcostheta_hists_QE_MULTI(i,j))
          call RemoveHist(dcostheta_hists_Delta_MULTI(i,j))
          call RemoveHist(dcostheta_hists_highRes_MULTI(i,j))
          call RemoveHist(dcostheta_hists_1piBG_MULTI(i,j))
          call RemoveHist(dcostheta_hists_2piBG_MULTI(i,j))
          call RemoveHist(dcostheta_hists_DIS_MULTI(i,j))
          call RemoveHist(dcostheta_hists_gen0_MULTI(i,j))
          !
          call RemoveHist(dE_hists_0pions(i,j))
          call RemoveHist(dE_hists_QE_0pions(i,j))
          call RemoveHist(dE_hists_Delta_0pions(i,j))
          call RemoveHist(dE_hists_highRes_0pions(i,j))
          call RemoveHist(dE_hists_Multi_0pions(i,j))
          call RemoveHist(dE_hists_QE_Multi_0pions(i,j))
          call RemoveHist(dE_hists_Delta_Multi_0pions(i,j))
          call RemoveHist(dE_hists_highRes_Multi_0pions(i,j))
          call RemoveHist(dE_hists_1X_0pions(i,j))
          call RemoveHist(dE_hists_QE_1X_0pions(i,j))
          call RemoveHist(dE_hists_Delta_1X_0pions(i,j))
          call RemoveHist(dE_hists_highRes_1X_0pions(i,j))
          call RemoveHist(dE_hists_2X_0pions(i,j))
          call RemoveHist(dE_hists_QE_2X_0pions(i,j))
          call RemoveHist(dE_hists_Delta_2X_0pions(i,j))
          call RemoveHist(dE_hists_highRes_2X_0pions(i,j))
       end do
    end do

    do j=-1,1
       call RemoveHist(dW_Npi_hists(j))
       call RemoveHist(dW_mupi_hists(j))
       call RemoveHist(dW_muN_hists(j))
    end do

    do j=1,max_speEvent
       do i=1,max_Hist
          call RemoveHist(dEnu_hist(j,i))
          call RemoveHist(dElepton_hist(j,i))
          call RemoveHist(dcoslepton_hist(j,i))
          call RemoveHist(dQ2lepton_hist(j,i))
          call RemoveHist(dQ2plepton_hist(j,i))
          call RemoveHist2D(dSigmaMC_EprimeCost_0pi(j,i))
       end do
    end do

  end subroutine cleanUp






  !****************************************************************************
  !****s* neutrinoAnalysis/neutrino_Analyze
  ! NAME
  ! subroutine neutrino_Analyze(Particles,finalFlag,num_runs_sameEnergy)
  ! INPUTS
  ! * type(particle), intent(in),dimension(:,:)  :: Particles
  !   -- Particles which shall be analyzed
  ! * logical, intent(in) :: finalFlag -- if .true. than the final output
  !   for a series of calls will be done
  ! * integer             :: num_runs_sameEnergy
  ! NOTES
  ! * This subroutine produces output for neutrino-nucleus scattering.
  !****************************************************************************
  subroutine neutrino_Analyze(Particles,finalFlag,num_runs_sameEnergy)
    use initNeutrino, only: getFirstEventRange, getNeutrinoInfo, nuEXP
    use particleDefinition
    use AnaEventDefinition
    use IDTable, only: nucleon,pion
    use rotation, only: get_Phi_Theta
    use degRad_conversion, only: radian
    use output, only: intToChar
    use history, only: history_getGeneration
    use neutrino_IDTable
    use neutrinoInfoStorage, only:  neutrinoInfoStorage_clear
    use neutrinoProdInfo, only: neutrinoProdInfo_Get!, neutrinoProdInfo_clear
    use expNeutrinofluxes, only: CCQE_recQs, CCQE_recQs_Delta, CCQE_recEnergy, &
         & CCQE_recEnergy_Delta,K2K_recEnergy, K2K_recQs
    use initNeutrino, only: includeQE, includeDELTA, includeRES, include1pi, &
         & includeDIS, include2p2hQE, include2p2hDelta, include2pi
    use minkowski, only: abs4Sq
    use vector, only: absVec
    use constants, only: pi
    use MultiplicityAnalysis, only: Multiplicity_Reset, Multiplicity_AddEvent, &
        & Multiplicity_Write
    use ZeroPionAnalysis, only: event_sigma_0pions, event_dsigma_de_0pions

    type(particle), intent(in),dimension(:,:) ,target :: Particles
    logical, intent(in) :: finalFlag
    integer, intent(in) :: num_runs_sameEnergy

    ! Local variables:
    integer, dimension (1:2) :: firstEvents
    type(tAnaEvent), Allocatable, dimension(:) :: events,events_QE,    &
         & events_Delta,events_highRES, &
         & events_1piBG,events_DIS,events_2p2h, events_2piBG
    type(tAnaEvent), Allocatable, dimension(:) :: events_gen0,    &
         & events_gen1,events_gen2, &
         & events_gen3ormore ! A list of all events
    type(particle), POINTER :: particlePointer
    integer :: i,j,first

    type(particle), Allocatable, dimension(:),target :: lepton, struckNuc

    integer, parameter :: max_generation=3
    integer :: generation, prod_id


    real  :: dPhi   ! Delta(phi) for dsigma/dOmega
    real  :: dTheta ! Delta(theta) for dsigma/dOmega

    real :: theta, phi

    real :: ekin_lepton

    real :: raiseFlagVariable

    integer,save :: numberOfCalls=0
    integer,save :: numberOfFinals=0

    logical, dimension(1:numStableParts) :: printflags

    ! In these hists we save the information
    ! which the "AnaEvent" subroutines are returning
    real, dimension(1:dimSigma,1:2),save :: sigma, sigma_QE,sigma_Delta,  &
       & sigma_highRES,sigma_1piBG, &
       & sigma_2piBG,sigma_DIS,sigma_2p2h
    real, dimension(1:dimSigma,1:2),save :: sigma_gen0,sigma_gen1,sigma_gen2, &
       & sigma_gen3ormore

    ! In these hists we save the information
    ! which the "ZeroPionAnalysis" subroutines are returning:
    !  in particular channels  with "k" protons and "n" neutrons
    real, dimension(1:dimSigma,1:2),save :: sigma_0pions, sigma_QE_0pions, &
       & sigma_Delta_0pions, &
       & sigma_highRES_0pions, sigma_1piBG_0pions, &
       & sigma_2piBG_0pions, sigma_DIS_0pions, sigma_2p2h_0pions

    real :: tkin
    integer :: ntk
    real, parameter :: EcosthetaMin=0.
    real, parameter :: EcosthetaMax=2.
    real, parameter :: dEcostheta=0.01

    real,dimension(-1:1) :: tsigmapion=0.
    real,dimension(-1:1,0:200) :: tksigmapion=0.
    real,dimension(-1:1) :: tsigmanucleon=0.          ! in DIS antiprotons are also possible
                                                      ! in final states
    real,dimension(-1:1,0:200) :: tksigmanucleon=0.   ! in DIS antiproton are also possible

    real,dimension(-1:1),save :: sum_tsigmapion=0.
    real,dimension(-1:1,0:200),save :: sum_tksigmapion=0.
    real,dimension(-1:1),save :: sum_tsigmanucleon=0. ! in DIS antiproton are also possible
                                                      ! in final states
    real,dimension(-1:1,0:200),save :: sum_tksigmanucleon=0.  ! in DIS antiproton are also
                                                              ! possible in final states

    real,dimension(0:200) :: Emiss_Fissum=0.
    real,dimension(0:200),save ::  sum_Emiss_Fissum=0.

    integer :: numberofneutrons_tot=0
    integer :: numberofprotons_tot=0
    integer :: numberofANTIprotons_tot=0
    integer :: numberofneutrons_out=0
    integer :: numberofprotons_out=0
    integer :: numberofANTIprotons_out=0
    integer,dimension(-1:1,0:max_generation) :: numberofnucleons_out  ! in DIS antiproton are
                                                                      !also possible in final states
    integer,save :: sum_numberofneutrons_tot=0
    integer,save :: sum_numberofprotons_tot=0
    integer,save :: sum_numberofANTIprotons_tot=0
    integer,save :: sum_numberofneutrons_out=0
    integer,save :: sum_numberofprotons_out=0
    integer,save :: sum_numberofANTIprotons_out=0
    integer,save,dimension(-1:1,0:max_generation) :: sum_numberofnucleons_out

    ! reconstruction of kinematics : as in MiniBooNN, K2K:
    type(histogram),save   :: H_Q2_real(1:max_speEvent,0:max_Hist), &
       & H_Q2_rec(1:max_speEvent,0:max_Hist)
    type(histogram2D),save :: H_Q2_rec_versus_real(1:max_speEvent,0:max_Hist)
    type(histogram),save   :: H_enu_real(1:max_speEvent,0:max_Hist), &
       & H_enu_rec(1:max_speEvent,0:max_Hist)
    type(histogram2D),save :: H_enu_rec_versus_real(1:max_speEvent,0:max_Hist)

    real :: Enureal, Enurec, Q2real, Q2rec, dummy, perweight
    real, dimension(0:3) :: lepIn_mom, lep_mom, boson_mom, nuc_mom
    ! momenta of the ingoing and outgoing lepton and intermediate boson

    integer :: m, iHist, Chrg_Nuc

    character*(10) :: Prefix_MultAna
    character(100) :: filename
    character(13) :: filename1

    real :: L, Posc_mumu,Posc_mue,Posc_mue_max,Posc_mue_antimax
    ! used for oscillation analysis

    type(histogram),save :: Oscmumu_enu_real(1:max_speEvent,0:max_Hist), &
       & Oscmumu_enu_rec(1:max_speEvent,0:max_Hist)
    type(histogram),save :: Oscmuemax_enu_real(1:max_speEvent,0:max_Hist), &
       & Oscmuemax_enu_rec(1:max_speEvent,0:max_Hist)
    type(histogram),save :: Oscmue_enu_real(1:max_speEvent,0:max_Hist), &
       & Oscmue_enu_rec(1:max_speEvent,0:max_Hist)
    type(histogram),save :: Oscmueantimax_enu_real(1:max_speEvent,0:max_Hist), &
       & Oscmueantimax_enu_rec(1:max_speEvent,0:max_Hist)

    dPhi=radian(10.)   ! Delta(phi) for dsigma/dOmega
    dTheta=radian(10.) ! Delta(theta) for dsigma/dOmega

    if (initflag) then
       call readinput
       initflag=.false.
       numberOfCalls=0
       numberOfFinals=0

       sum_numberofneutrons_tot=0
       sum_numberofprotons_tot=0
       sum_numberofANTIprotons_tot=0
       sum_numberofneutrons_out=0
       sum_numberofprotons_out=0
       sum_numberofANTIprotons_out=0
       sum_numberofnucleons_out=0
       sum_tsigmapion=0.
       sum_tksigmapion=0.
       sum_tsigmanucleon=0.
       sum_tksigmanucleon=0.
       sum_Emiss_Fissum=0.

       printFlags=.false.

!  Now switches for printout of cross sections for various hadrons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  indices in printflags(i) refer to ordering in field particleIds
!  in module AnaEvent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if (forpion) printFlags(1)=.true. ! pion
       if (foreta) printFlags(2)=.true. ! eta
       if (forkaon) printFlags(3)=.true. ! kaon
       if (forkaonBar) printFlags(4)=.true. ! kaonBar
       if (forDmeson) printFlags(5)=.true. ! Dmeson
       if (forDbar) printFlags(6)=.true. ! Dbar meson
       if (forDs_plus) printFlags(7)=.true. ! Ds+ meson
       if (forDs_minus) printFlags(8)=.true. ! Ds- meson
       if (fornucleon) printFlags(9)=.true. ! nucleon
       if (forLambda) printFlags(10)=.true. ! Lambda
       if (forSigmaResonance) printFlags(11)=.true. ! Sigma
       if (forXi) printFlags(12)=.true. ! Xi baryon
       if (forOmegaResonance) printFlags(13)=.true. ! Omega baryon
       call set_particleIDs_flag(printFlags)


       if (Xsection_analysis) then

       !***********************************************************************
       !****o* neutrinoAnalysis/neutrino_total_Xsection_multiplicities.dat
       ! NAME
       ! file neutrino_total_Xsection_multiplicities.dat
       ! PURPOSE
       ! The file is produced in the runs with eventtype=5=neutrino .
       !
       ! The file shows the cross sections for  preselected final states
       !
       ! Units:
       ! * For process_ID=CC and NC the units 10^{-38} cm^2 for integrated xsec
       !   (10^{-38)cm^2/GeV for dsigma/dElepton,
       !   10^{-38)cm^2/GeV^2 for dsigma/dQ^2, and so on)
       ! * For process_ID=EM the units are nanobarns=10^{-33}cm^2
       !
       ! Columns:
       ! * #1: variable which was raised
       !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode,
       !   Elepton for nuXsectionMode=2=dSigmadQsdElepton  and so on)
       ! * #2-#119: see description in AnaEvent.f90, subroutine event_sigma
       !   OR in the output file sigma.dat
       !   The description of some columns is given below
       !   In each channel the outgoing lepton is presupposed
       !   (unless explicitely stated otherwise)
       !
       ! Some columns in detail:
       ! * #2: 1 pi-, no other pions of any charge and anything else
       ! * #3: 1 pi0, no other pions of any charge and anything else
       ! * #4: 1 pi+, no other pions of any charge and anything else
       ! * #5: 1 eta  and anything else
       ! * #14: 1 neutron, no other nucleons and anything else
       ! * #15: 1 proton,  no other nucleons and anything else
       ! * #68: 1 nucleon and 1 pion, no other pions or nucleons, anything else
       ! * #68: 1 proton (no other nucleons) and 0 pions and anything else
       !   (QE-like in Cherenkov detector)
       ! * #69: 0 pions and  anything else
       ! * #70: at least  1 pi-, any number of pi0 and/or pi+ , anything else,
       !   each event is counted once
       ! * #71: at least  1 pi0, any number of pi- and/or pi+ , anything else,
       !   each event is counted once
       ! * #72: at least  1 pi+, any number of pi- and/or pi0 , anything else,
       !   each event is counted once
       ! * #91: at least  1 pi0, any number of pi- and/or pi+ , anything else,
       !   each pi0 is counted
       ! * #96: no nucleons, anything else
       ! * #97: 5 or more nucleons, anything else
       ! * #120-#239:  errors to Columns #2-#119
       !***********************************************************************

       open(10,File='neutrino_total_Xsection_multiplicities.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)

       !***********************************************************************
       !****o* neutrinoAnalysis/neutrino_total_Xsection_multiplicities_QE.dat
       ! NAME
       ! file neutrino_total_Xsection_multiplicities_QE.dat
       ! PURPOSE
       ! The same as neutrino_total_Xsection_multiplicities.dat
       ! but for QE events (=the first interaction was quasielastic or
       ! elastic scattering)
       !***********************************************************************
       if (includeQE) then
       open(10,File='neutrino_total_Xsection_multiplicities_QE.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)
       end if

       !***********************************************************************
       !****o* neutrinoAnalysis/neutrino_total_Xsection_multiplicities_Delta.dat
       ! NAME
       ! file neutrino_total_Xsection_multiplicities_Delta.dat
       ! PURPOSE
       ! The same as neutrino_total_Xsection_multiplicities.dat
       ! but for Delta events (=the first interaction was production of the
       ! Delta resonance)
       !***********************************************************************
       if (includeDELTA) then
       open(10,File='neutrino_total_Xsection_multiplicities_Delta.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)
       end if

       !***********************************************************************
       !****o* neutrinoAnalysis/neutrino_total_Xsection_multiplicities_highRES.dat
       ! NAME
       ! file neutrino_total_Xsection_multiplicities_highRES.dat
       ! PURPOSE
       ! The same as neutrino_total_Xsection_multiplicities.dat
       ! but for highRES events (=the first interaction was production any
       ! resonance beyond Delta)
       !***********************************************************************
       if (includeRES) then
       open(10,File='neutrino_total_Xsection_multiplicities_highRES.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)
       end if

       !***********************************************************************
       !****o* neutrinoAnalysis/neutrino_total_Xsection_multiplicities_1piBG.dat
       ! NAME
       ! file neutrino_total_Xsection_multiplicities_1piBG.dat
       ! PURPOSE
       ! The same as neutrino_total_Xsection_multiplicities.dat
       ! but for 1piBG events (=the first interaction was background production of
       ! 1-pion final state)
       !***********************************************************************
       if (include1pi) then
       open(10,File='neutrino_total_Xsection_multiplicities_1piBG.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)
       end if

       !***********************************************************************
       !****o* neutrinoAnalysis/neutrino_total_Xsection_multiplicities_DIS.dat
       ! NAME
       ! file neutrino_total_Xsection_multiplicities_DIS.dat
       ! PURPOSE
       ! The same as neutrino_total_Xsection_multiplicities.dat
       ! but for DIS events (=the first interaction was DIS)
       !***********************************************************************
       if (includeDIS) then
       open(10,File='neutrino_total_Xsection_multiplicities_DIS.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)
       end if

       !***********************************************************************
       !****o* neutrinoAnalysis/neutrino_total_Xsection_multiplicities_2p2h.dat
       ! NAME
       ! file neutrino_total_Xsection_multiplicities_2p2h.dat
       ! PURPOSE
       ! The same as neutrino_total_Xsection_multiplicities.dat
       ! but for 2particle-2hole events
       !***********************************************************************
       if (include2p2hQE .or. include2p2hDelta) then
       open(10,File='neutrino_total_Xsection_multiplicities_2p2h.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)
       end if


       !***********************************************************************
       !****o* neutrinoAnalysis/neutrino_total_Xsection_multiplicities_gen0.dat
       ! NAME
       ! file neutrino_total_Xsection_multiplicities_gen0.dat
       ! PURPOSE
       ! The same as neutrino_total_Xsection_multiplicities.dat for particles
       ! of the 0th generation
       ! The definition of generation is given in history.f90, description
       ! of the module "history"
       !***********************************************************************
       open(10,File='neutrino_total_Xsection_multiplicities_gen0.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)

       !***********************************************************************
       !****o* neutrinoAnalysis/neutrino_total_Xsection_multiplicities_gen1.dat
       ! NAME
       ! file neutrino_total_Xsection_multiplicities_gen1.dat
       ! PURPOSE
       ! The same as neutrino_total_Xsection_multiplicities.dat for particles
       ! of the 1st generation
       ! The definition of generation is given in history.f90, description
       ! of the module "history"
       !***********************************************************************
       open(10,File='neutrino_total_Xsection_multiplicities_gen1.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)

       !***********************************************************************
       !****o* neutrinoAnalysis/neutrino_total_Xsection_multiplicities_gen2.dat
       ! NAME
       ! file neutrino_total_Xsection_multiplicities_gen2.dat
       ! PURPOSE
       ! The same as neutrino_total_Xsection_multiplicities.dat for particles
       ! of the 2nd generation
       ! The definition of generation is given in history.f90, description
       ! of the module "history"
       !***********************************************************************
       open(10,File='neutrino_total_Xsection_multiplicities_gen2.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)

       !***********************************************************************
       !****o* neutrinoAnalysis/neutrino_total_Xsection_multiplicities_gen3ormore.dat
       ! NAME
       ! file neutrino_total_Xsection_multiplicities_gen3ormore.dat
       ! PURPOSE
       ! The same as neutrino_total_Xsection_multiplicities.dat for particles
       ! of the 3rd-or-more generation
       ! The definition of generation is given in history.f90, description
       ! of the module "history"
       !***********************************************************************
       open(10,File='neutrino_total_Xsection_multiplicities_gen3ormore.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)

       open(10,File='neutrino_total_Xsection.dat')
       write(10,*)'# Total xsec for outgoing pions and nucleons:&
          & each pion/nucleon contribute to the xsec'
       write(10,*)'# 1: raiseVariable  2: sigma(pi-)    3:pi0    4:pi+  &
          & 5:antiproton 6:neutron, 7:proton'
       close(10)

       !***********************************************************************
       !****o* neutrinoAnalysis/neutrino_kinetic_energy_Xsection.dat
       ! NAME
       ! file neutrino_kinetic_energy_Xsection.dat
       ! PURPOSE
       ! The file is produced in the runs with eventtype=5=neutrino .
       !
       ! The file shows the kinetic energy differential cross sections after
       ! final state interactions for pions, protons and neutrons
       !
       ! Units:
       ! * For process_ID=CC and NC: 10^{-38} cm^2/GeV for integrated
       !   xsec and so on ..
       ! * For process_ID=EM: nanobarns=10^{-33}cm^2/ GeV and so on.
       !
       ! Columns:
       ! * #1: variable which was raised
       !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode,
       !   Elepton for nuXsectionMode=2=dSigmadQsdElepton and so on)
       ! * #2: kinetic energy in GeV
       ! * #3: xsec for events with at least one pi-  in the final state
       !   (coincides with column 2 in
       !   diff_XXX_dSigma_dEkin_pi_charge_-1_MULTI.dat)
       ! * #4: xsec for events with at least one pi0
       !   (coincides with column 2 in
       !   diff_XXX_dSigma_dEkin_pi_charge_+0_MULTI.dat)
       ! * #5: xsec for events with at least one pi+
       !   (coincides with column 2 in
       !   diff_XXX_dSigma_dEkin_pi_charge_+1_MULTI.dat)
       ! * #6: xsec for events with at least one neutron
       !   (coincides with column 2 in
       !   diff_XXX_dSigma_dEkin_N_charge_+0_MULTI.dat)
       ! * #7: xsec for events with at least one proton
       !   (coincides with column 2 in
       !   diff_XXX_dSigma_dEkin_N_charge_+1_MULTI.dat)
       !
       ! HERE  XXX=000,001, ... is the count of the raise variable
       !***********************************************************************
       open(10,File='neutrino_kinetic_energy_Xsection.dat')
       write(10,*)'# order of entries: 1:raise variable  2:ekin, 3:sig pi-, &
          & 4:sig pi0, 5:sig pi+, 6:sig neut, 7:sig prot'
       close(10)

       !***********************************************************************
       !****o* neutrinoAnalysis/neutrino_Xsection_numbers.dat
       ! NAME
       ! file neutrino_Xsection_numbers.dat
       ! PURPOSE
       ! The file is produced in the runs with eventtype=5=neutrino .
       !
       ! The file shows number of nucleons produced in one run before and
       ! after final state interactions
       ! (should be devided by numEnsembles to obtain the average multiplicity
       ! per target nucleus and divided further by target_A to obtain the
       ! average multiplicity per target nucleon)
       !
       ! Columns:
       ! * #1: variable which was raised
       !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode,
       !   Elepton for nuXsectionMode=2=dSigmadQsdElepton  and so on)
       ! * #2:  all the protons produced  before final state interactions (FSI)
       ! * #3: protons that made it out of the nucleus after FSI
       ! * #4: protons produced in generation 0 and made it out of the nucleus
       !   after FSI
       ! * #5: protons produced in generation 1 and made it out of the nucleus
       !   after FSI
       ! * #6: protons produced in generation 2 and made it out of the nucleus
       !   after FSI
       ! * #7: protons produced in generation 3-or-more and made it out of the
       !   nucleus after FSI
       ! * #8: all the neutrons produced before FSI
       ! * #9:  neutrons that made it out of the nucleus after FSI
       ! * #10: neutrons produced in generation 0 and made it out of the
       !   nucleus after FSI
       ! * #11: neutrons produced in generation 1 and made it out of the
       !   nucleus after FSI
       ! * #12: neutrons produced in generation 2 and made it out of the
       !   nucleus after FSI
       ! * #13: neutrons produced in generation 3-or-more and made it out of
       !   the nucleus after FSI
       ! * #14: all the antiprotons produced before FSI
       ! * #15: antiprotons that made it out of the nucleus after FSI
       ! * #16: antiprotons produced in generation 0 and made it out of the
       !   nucleus after FSI
       ! * #17: antiprotons produced in generation 1 and made it out of the
       !   nucleus after FSI
       ! * #18: antiprotons produced in generation 2 and made it out of the
       !   nucleus after FSI
       ! * #19: antiprotons produced in generation 3-or-more and made it out
       !   of the nucleus after FSI
       !***********************************************************************
       open(10,File='neutrino_Xsection_numbers.dat')
       write(10,*) '#1:raiseFlagVariable, #2:p_tot, #3:p_out, #4:p_out_gen0, &
            & #5:p_out_gen1, &
            & #6:p_out_gen2, #7:p_out_gen3ormore,&
            & #8:n_tot, #9:n_out, #10:n_out_gen0, #11:n_out_gen1 ,#12:n_out_gen2,&
            & #13:n_out_gen3ormore, &
            & #14:barp_tot, #15:barp_out, #16:barp_out_gen0, #17:barp_out_gen1, &
            & #18:barp_out_gen2,#19:bar p_out_gen3ormore'
       close(10)

       end if ! Xsection_analysis

       if (Fissum_analysis) then
          open(10,File='neutrino_Emiss_spectrum.dat')
          write(10,*)'# order of entries: ekin, sig prot'
          close(10)
       end if

       if (ZeroPion_analysis) then
          call WriteHeader10(.true., 'neutrino_0pions.dat')
          call WriteHeader10(includeQE,'neutrino_0pions_QE.dat')
          call WriteHeader10(includeDELTA,'neutrino_0pions_Delta.dat')
          call WriteHeader10(includeRES,'neutrino_0pions_highRES.dat')
          call WriteHeader10(include1pi,'neutrino_0pions_1piBG.dat')
          call WriteHeader10(include2pi,'neutrino_0pions_2piBG.dat')
          call WriteHeader10(includeDIS,'neutrino_0pions_DIS.dat')
          call WriteHeader10(include2p2hQE .or. include2p2hDelta,&
          & 'neutrino_0pions_2p2h.dat')
       end if

       if (reconstruct_neutrino_energy .and. specificEvent_Analysis) then
          if (binsizeQ2.lt.0.001) binsizeQ2=0.02
          if (binsizeEnu.lt.0.001) binsizeEnu=0.01

          do m=1, max_speEvent
             if (.not.includeSpeEvent(m)) cycle
             do iHist=0, max_Hist
                if (.not.includeHist(iHist)) cycle

                ! Q2 reconstruction
                call CreateHist(H_Q2_real(m,iHist),  &
                   &'true Q2 for a specific event',0.,maxQ2, binsizeQ2)
                call CreateHist(H_Q2_rec(m,iHist),  'reconstructed Q2  &
                   & for a specific event',0., maxQ2,binsizeQ2)
                call CreateHist2D(H_Q2_rec_versus_real(m,iHist), &
                     'reconstructed Q2 versus real Q2 for a specific event', &
                     (/0.,0./),(/maxQ2,maxQ2/),  (/binsizeQ2,binsizeQ2/))

                ! neutrino energy reconstruction
                if (nuExp>0) call CreateHist(H_enu_real(m,iHist), &
                   & 'true Enu for a specific event',0.,maxEnu,binsizeEnu)
                call CreateHist(H_enu_rec(m,iHist),  'reconstructed Enu&
                   & for a specific event', 0.,maxEnu,binsizeEnu)
                if (nuExp>0) call CreateHist2D(H_enu_rec_versus_real(m,iHist), &
                     'reconstructed Enu versus real Enu for a specific event', &
                     (/0.,0./),(/maxEnu,maxEnu/),  (/binsizeEnu,binsizeEnu/))

                ! oscillations:  nu_mu survival, nu_e appearence
                if (OSC(nuEXP)) then
                   call CreateHist(Oscmumu_enu_real(m,iHist), &
                        'nu_mu survival versus true energy',0.,maxEnu,binsizeEnu)
                   call CreateHist(Oscmumu_enu_rec(m,iHist) , &
                        'nu_mu survival versus reconstructed energy',0.,maxEnu, &
                        & binsizeEnu)
                   call CreateHist(Oscmuemax_enu_real(m,iHist), &
                        & 'nu_e appearence versus true energy for&
                        & delta_CP=pi/2',0., maxEnu,binsizeEnu)
                   call CreateHist(Oscmuemax_enu_rec(m,iHist),  &
                        &  'nu_e appearence versus reconstructed energy&
                        & for delta_CP=pi/2', 0.,maxEnu,binsizeEnu)
                   call CreateHist(Oscmue_enu_real(m,iHist), &
                        'nu_e appearence versus true energy',0.,maxEnu,binsizeEnu)
                   call CreateHist(Oscmue_enu_rec(m,iHist), &
                        'nu_e appearence versus reconstructed energy',0., &
                        & maxEnu,binsizeEnu)
                   call CreateHist(Oscmueantimax_enu_real(m,iHist), &
                        &  'nu_e appearence versus true energy for&
                        & delta_CP=-pi/2',0., maxEnu,binsizeEnu)
                   call CreateHist(Oscmueantimax_enu_rec(m,iHist),  &
                        &  'nu_e appearence versus reconstructed energy&
                        & for delta_CP=-pi/2', 0.,maxEnu,binsizeEnu)
                end if

             end do !iHist
          end do !m
       end if


       call CreateHist(hNucleonVacuumMass, 'mass of nucleons', 0.65,1.2,0.001)

       call Multiplicity_Reset

    end if

    call getNeutrinoInfo(raiseFlagVariable)

    numberOfCalls=numberOfCalls+1

    initHists = (numberOfCalls.eq.1) ! whether we initialize the histograms


    write(*,*) '################### NEUTRINO ANALYSIS STARTS ##################'
    write(*,*) ' number of calls: ', numberofCalls,' number of finals: ',&
              & numberofFinals
    if (inclusiveAnalysis) write(*,*) 'NOTE: we do inclusiveAnalysis!!!!!!!!!!!'
    write(*,'(a,F12.4)') 'kineticEnergyDetectionThreshold_nucleon ', &
       & kineticEnergyDetectionThreshold_nucleon
    write(*,'(a,F12.4)')'radialScale',radialScale




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! total (tsigma) and kinetic energy differential (tksigma) cross section
! for pions and nucleons

    !nullify everything
    tsigmapion=0.
    tksigmapion=0.
    tsigmanucleon=0.
    tksigmanucleon=0.
    numberofneutrons_tot=0
    numberofprotons_tot=0
    numberofANTIprotons_tot=0
    numberofneutrons_out=0
    numberofprotons_out=0
    numberofANTIprotons_out=0
    numberofnucleons_out=0
    Emiss_Fissum=0.

    Q2real=0.
    Enureal=0.
    Q2rec=0.
    Enurec=0.



! Now loops over all ensembles (dim=1) and all particles (dim=2)

    do i=lbound(Particles,dim=1),ubound(Particles,dim=1)
       do j=lbound(Particles,dim=2),ubound(Particles,dim=2)
          if (Particles(i,j)%ID.le.0) cycle

   ! The following 'if' selects only nucleons
          if (Particles(i,j)%ID.eq.1) then
             call AddHist(hNucleonVacuumMass,Particles(i,j)%mass, &
                         & Particles(i,j)%perweight)
             select case (Particles(i,j)%charge)
             case (-1)
                numberofANTIprotons_tot=numberofANTIprotons_tot+1
             case (0)
                numberofneutrons_tot=numberofneutrons_tot+1
             case (1)
                numberofprotons_tot=numberofprotons_tot+1
             end select
          end if

          if (IsBound(particles(i,j))) cycle
          if (IsBelowThreshold(particles(i,j))) cycle

   ! The following call to function neutrinoProdInfo_Get retrieves information
   ! on the production process for particle j in ensemble i

          if (.not.neutrinoProdInfo_Get(Particles(i,j)%firstEvent,prod_id,&
            &  dummy,lepIn_mom,lep_mom,boson_mom, nuc_mom, Chrg_Nuc)) then
             call TRACEBACK('error in getting production info')
          end if

 ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Now acceptance cuts in kinetic energy and angle for outgoing lepton
 ! Events are accepted only if
 ! ekin_lepton > threshold energy and lepton-angle < threshold angle
 !
          ekin_lepton=lep_mom(0)-sqrt( max(0.,abs4Sq(lep_mom)) )
          if (ekin_lepton.lt.kineticEnergyDetectionThreshold_lepton) cycle

          if ( lep_mom(3)/absVec(lep_mom(1:3)) &
          & .lt. cos(radian(AngleUpperDetectionThresholdDegrees_lepton)) ) cycle
 ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          tkin=Particles(i,j)%momentum(0)-Particles(i,j)%mass
          if (tkin.lt.0.) tkin = 0.  !only relevant if inclusiveAnalysis=.true.

          ntk=min(int(tkin/dEkin),200)

 !apply history information
          generation=history_getGeneration(Particles(i,j)%history)
          if (generation.gt.max_generation) generation=max_generation

          select case (Particles(i,j)%ID)
          case (pion)

             !total cross section
             tsigmapion(Particles(i,j)%charge)=tsigmapion(Particles(i,j)%charge) &
             & +Particles(i,j)%perweight

             !kinetic energy differential xsection
             tksigmapion(Particles(i,j)%charge,ntk)=&
                  & Particles(i,j)%perweight/dEkin   &
                  & +tksigmapion(Particles(i,j)%charge,ntk)

          case (nucleon)

             select case (Particles(i,j)%charge)
             case (-1)
                numberofANTIprotons_out=numberofANTIprotons_out+1
             case (0)
                numberofneutrons_out=numberofneutrons_out+1
             case (1)
                numberofprotons_out=numberofprotons_out+1
             end select

             numberofnucleons_out(Particles(i,j)%charge,generation)  &
             & =numberofnucleons_out(Particles(i,j)%charge,generation)+1

             !total cross section

             tsigmanucleon(Particles(i,j)%charge)  &
             & =tsigmanucleon(Particles(i,j)%charge)+Particles(i,j)%perweight

             !kinetic energy differential xsection
             tksigmanucleon(Particles(i,j)%charge,ntk)=&
                  &Particles(i,j)%perweight/dEkin  &
                  & +tksigmanucleon(Particles(i,j)%charge,ntk)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fissum analysis for O16(e,e'p) data, Phys.Rev.C70:034606,2004

             if (Fissum_analysis) then
                if (Particles(i,j)%charge.ne.1) cycle
                call get_phi_Theta(particles(i,j)%momentum(1:3),theta,phi)

                if (.not.(abs(cos(theta)-cos(pi/180.*38.45)).lt.0.02/2.)) cycle

                Emiss_Fissum(ntk)=Particles(i,j)%perweight/dEkin/0.02 &
                                 & +Emiss_Fissum(ntk)
             end if !Fissum analysis
! End Fissum analysis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          end select

       end do ! loop over particle vector
    end do ! loop over particle vector: ensemble loop





    ! summing over several runs at the same energy
    sum_numberofneutrons_tot=numberofneutrons_tot+sum_numberofneutrons_tot
    sum_numberofprotons_tot=numberofprotons_tot+sum_numberofprotons_tot
    sum_numberofANTIprotons_tot=numberofprotons_tot+sum_numberofANTIprotons_tot
    sum_numberofneutrons_out=numberofneutrons_out+sum_numberofneutrons_out
    sum_numberofprotons_out=numberofprotons_out+sum_numberofprotons_out
    sum_numberofANTIprotons_out=numberofprotons_out+sum_numberofANTIprotons_out
    sum_numberofnucleons_out=numberofnucleons_out+sum_numberofnucleons_out

    sum_tsigmapion=sum_tsigmapion+tsigmapion
    sum_tksigmapion=sum_tksigmapion+tksigmapion
    sum_tsigmanucleon=sum_tsigmanucleon+tsigmanucleon
    sum_tksigmanucleon=sum_tksigmanucleon+tksigmanucleon

    sum_Emiss_Fissum=sum_Emiss_Fissum+Emiss_Fissum


    if (Xsection_analysis) then

       open(10,File='neutrino_total_Xsection.dat',position='append')
       if (numberofcalls.ne.1) backspace(10)
       write(10,'(10g13.5)') raiseFlagVariable, &
            sum_tsigmapion/real(numberofcalls), &
            sum_tsigmanucleon/real(numberofcalls)
       close(10)

       open(11,File='neutrino_kinetic_energy_Xsection.dat',position='append')
       if (numberofcalls.ne.1) then
          do i=0,200
             backspace(11)
          end do
       end if
       do i=0,200
          write(11,'(10g13.5)') raiseFlagVariable, &
               (float(i)+0.5)*dEkin,sum_tksigmapion(-1,i)/real(numberofcalls), &
               sum_tksigmapion(0,i)/real(numberofcalls), &
               sum_tksigmapion(1,i)/real(numberofcalls), &
               sum_tksigmanucleon(0,i)/real(numberofcalls), &
               sum_tksigmanucleon(1,i)/real(numberofcalls)
       end do
       close(11)

       open(12,File='neutrino_Xsection_numbers.dat',position='append')
       if (numberofcalls.ne.1) backspace(12)
       write(12,'(30g13.5)') raiseFlagVariable, &
            sum_numberofprotons_tot/numberofcalls, &
            sum_numberofprotons_out/numberofcalls, &
            sum_numberofnucleons_out(1,:)/numberofcalls, &
            sum_numberofneutrons_tot/numberofcalls, &
            sum_numberofneutrons_out/numberofcalls, &
            sum_numberofnucleons_out(0,:)/numberofcalls, &
            sum_numberofANTIprotons_tot/numberofcalls, &
            sum_numberofANTIprotons_out/numberofcalls, &
            sum_numberofnucleons_out(-1,:)/numberofcalls
       close(12)

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fissum analysis for O16(e,e'p) data, Phys.Rev.C70:034606,2004

    if (Fissum_analysis) then
       open(10,File='neutrino_Emiss_spectrum.dat',position='append')
       if (numberofcalls.ne.1) then
          do i=0,200
             backspace(10)
          end do
       end if
       do i=0,200
          write(10,'(10g13.5)')raiseFlagVariable, &
             & (float(i)+0.5)*dEkin,sum_Emiss_Fissum(i)/real(numberofcalls)
       end do
       close(10)
    end if
! End Fissum analysis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call WriteHist(hNucleonVacuumMass,file="NucleonVacuumMass.dat",  &
                  & mul = 1.0/numberofcalls)

    if (finalFlag) then
       sum_numberofneutrons_tot=0
       sum_numberofprotons_tot=0
       sum_numberofANTIprotons_tot=0
       sum_numberofneutrons_out=0
       sum_numberofprotons_out=0
       sum_numberofANTIprotons_out=0
       sum_numberofnucleons_out=0
       sum_tsigmapion=0.
       sum_tksigmapion=0.
       sum_tsigmanucleon=0.
       sum_tksigmanucleon=0.
       sum_Emiss_Fissum=0.
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  !****************************************************************************
  ! now event analyses with more detailed information on multiplicities,
  ! generations and production process
  !
  ! (1) Setting up the particles into the events
  !
  ! This is done with the help of %firstEvent:
  ! Particles stemming from the same event get in the init the same %firstEvent
  ! entry.
  ! During the run %firstEvent stays constant and is inherited during
  ! collisions.

    firstEvents=getFirstEventRange()

    allocate(lepton(firstEvents(1):firstEvents(2)))
    allocate(struckNuc(firstEvents(1):firstEvents(2)))

    call setToDefault(lepton)
    call setToDefault(struckNuc)

    allocate(Events(firstEvents(1):firstEvents(2)))
    allocate(Events_QE(firstEvents(1):firstEvents(2)))
    allocate(Events_Delta(firstEvents(1):firstEvents(2)))
    allocate(Events_highRES(firstEvents(1):firstEvents(2)))
    allocate(Events_1piBG(firstEvents(1):firstEvents(2)))
    allocate(Events_2piBG(firstEvents(1):firstEvents(2)))
    allocate(Events_DIS(firstEvents(1):firstEvents(2)))
    allocate(Events_2p2h(firstEvents(1):firstEvents(2)))
    allocate(Events_gen0(firstEvents(1):firstEvents(2)))
    allocate(Events_gen1(firstEvents(1):firstEvents(2)))
    allocate(Events_gen2(firstEvents(1):firstEvents(2)))
    allocate(Events_gen3ormore(firstEvents(1):firstEvents(2)))

    do i=firstEvents(1),firstEvents(2)
       call event_init(events(i))
       call event_init(events_QE(i))
       call event_init(events_Delta(i))
       call event_init(events_highRES(i))
       call event_init(events_1piBG(i))
       call event_init(events_2piBG(i))
       call event_init(events_DIS(i))
       call event_init(events_2p2h(i))
       call event_init(events_gen0(i))
       call event_init(events_gen1(i))
       call event_init(events_gen2(i))
       call event_init(events_gen3ormore(i))
    end do


 !
 ! First get 'in' and 'out' lepton momenta and momentum transfer for each event;
 ! the loop runs over all events
 !

    do i=firstEvents(1),firstEvents(2)
       lepton(i)%firstEvent=i
       call get_init_namelist(outLepton_ID=lepton(i)%ID, &
                             & outLepton_charge=lepton(i)%charge)
       if (.not.neutrinoProdInfo_Get(i, prod_id,lepton(i)%perweight,lepIn_mom, &
          lepton(i)%momentum,boson_mom, struckNuc(i)%momentum, &
          & struckNuc(i)%charge)) then
          call TRACEBACK('error in getting perweight')
       end if


 ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Now acceptance cuts in kinetic energy and angle for outgoing lepton
 ! Events are accepted only if
 ! ekin_lepton > threshold energy and lepton-angle < threshold angle
 !

       ekin_lepton=lepton(i)%momentum(0)-sqrt( max(0.,abs4Sq(lepton(i)%momentum)) )
       if (ekin_lepton.lt.kineticEnergyDetectionThreshold_lepton) cycle
       if ( lepton(i)%momentum(3)/absVec(lepton(i)%momentum(1:3)) &
          & .lt. cos(radian(AngleUpperDetectionThresholdDegrees_lepton)) ) cycle

 ! lepton is added to event only if its energy and angle have been accepted
 ! cuts affect FinalEvents.dat
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 ! Analyze events separately for first interaction: QE, Delta, DIS, .....

 !  Now put leptons into specific event types (QE, Delta, ..)
       particlePointer=>lepton(i)

       call event_add(events(i),particlePointer)

 !  Now put struck nucleon into specific event types
       particlePointer=>struckNuc(i)
       struckNuc(i)%ID=1

       call event_add(events(i),particlePointer)

       select case (prod_id)
       case (1)
          call event_add(events_QE(i),particlePointer)
       case (2)
          call event_add(events_Delta(i),particlePointer)
       case (3:31)
          call event_add(events_highRES(i),particlePointer)
       case (32:33)
          call event_add(events_1piBG(i),particlePointer)
       case (34)
          call event_add(events_DIS(i),particlePointer)
       case (35:36)
          call event_add(events_2p2h(i),particlePointer)
       case (37)
          call event_add(events_2piBG(i),particlePointer)
       case default
          write(*,*) 'prod_id =', prod_id
          call TRACEBACK('strange prod_id')
       end select

       call event_add(events_gen0(i),particlePointer)
       call event_add(events_gen1(i),particlePointer)
       call event_add(events_gen2(i),particlePointer)
       call event_add(events_gen3ormore(i),particlePointer)

    end do



 !
 ! Now put all other particles into specific event types
 !

 ! Now loops over all ensembles (dim=1) and all particles (dim=2)

    do i=lbound(Particles,dim=1),ubound(Particles,dim=1)
       do j=lbound(Particles,dim=2),ubound(Particles,dim=2)
          if (Particles(i,j)%ID.le.0) cycle

          first=Particles(i,j)%firstEvent
 !  first is the number of the first event of particle j in ensemble i

          if (IsBound(particles(i,j))) cycle
          if (IsBelowThreshold(particles(i,j))) cycle

 ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Now acceptance cuts in kinetic energy and angle for outgoing lepton
 ! Events are accepted only if
 ! ekin_lepton > threshold energy and lepton-angle < threshold angle

          ekin_lepton=lepton(first)%momentum(0)  &
          &  - sqrt( max(0.,abs4Sq(lepton(first)%momentum)) )
          if (ekin_lepton.lt.kineticEnergyDetectionThreshold_lepton) cycle
          if ( lepton(first)%momentum(3)/absVec(lepton(first)%momentum(1:3)) &
          & .lt. cos(radian(AngleUpperDetectionThresholdDegrees_lepton)) ) cycle

 ! lepton is added to event only if its energy and angle have been accepted
 ! cuts affect FinalEvents.dat
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          particlePointer=>Particles(i,j)

    ! Add particle to the event with its firstEvent index.

          call event_add(events(first),particlePointer)

          !apply history information
          generation=history_getGeneration(Particles(i,j)%history)
          if (.not.neutrinoProdInfo_Get(Particles(i,j)%firstEvent,prod_id,dummy,lepIn_mom, &
                                       lep_mom,boson_mom,nuc_mom,Chrg_Nuc)) then
             call TRACEBACK('error in getting production info')
          end if

          if (generation.ge.max_generation) generation=max_generation

    ! to simplify output individual resonance contributions are all lumped together in
    ! the highRes component, this is controlled by array K2Hist in initneutrino.f90

    ! now generation is the real particle generation and prod_id contains
    ! the information, in which process the particle was produced (QE=1, Delta=2, highRES=3:31,
    ! BG=32:33)


          select case (prod_id)
          case (1)
             call event_add(events_QE(first),particlePointer)
          case (2)
             call event_add(events_Delta(first),particlePointer)
          case (3:31)
             call event_add(events_highRES(first),particlePointer)
          case (32:33)
             call event_add(events_1piBG(first),particlePointer)
          case (34)
             call event_add(events_DIS(first),particlePointer)
          case (35:36)
             call event_add(events_2p2h(first),particlePointer)
          case (37)
             call event_add(events_2piBG(first),particlePointer)
          case default
             write(*,*) 'prod_id =', prod_id
             call TRACEBACK('strange prod_id')
          end select

          select case (generation)
          case (0)
             call event_add(events_gen0(first),particlePointer)
          case (1)
             call event_add(events_gen1(first),particlePointer)
          case (2)
             call event_add(events_gen2(first),particlePointer)
          case (3)
             call event_add(events_gen3ormore(first),particlePointer)
          case default
             write(*,*) 'generation =', generation
             call TRACEBACK('strange generation')
          end select

       end do
    end do




! (2) Use the list "events" to evaluate total cross sections

    call event_sigma(events,sigma,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_QE,sigma_QE,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_Delta,sigma_Delta,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_highRES,sigma_highRES,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_1piBG,sigma_1piBG,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_2piBG,sigma_2piBG,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_DIS,sigma_DIS,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_2p2h,sigma_2p2h,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_gen0,sigma_gen0,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_gen1,sigma_gen1,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_gen2,sigma_gen2,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_gen3ormore,sigma_gen3ormore,initHists, &
         numberOfCalls, &
         identifier=raiseFlagVariable)

    if (Xsection_analysis) then
       call PrintVals10(.true., raiseFlagVariable, numberOfCalls, &
            'neutrino_total_Xsection_multiplicities.dat', sigma)
       call PrintVals10(includeQE, raiseFlagVariable, numberOfCalls, &
            'neutrino_total_Xsection_multiplicities_QE.dat', sigma_QE)
       call PrintVals10(includeDelta, raiseFlagVariable, numberOfCalls, &
            'neutrino_total_Xsection_multiplicities_Delta.dat', sigma_Delta)
       call PrintVals10(includeRes, raiseFlagVariable, numberOfCalls, &
            'neutrino_total_Xsection_multiplicities_highRES.dat', sigma_highRES)
       call PrintVals10(include1pi, raiseFlagVariable, numberOfCalls, &
            'neutrino_total_Xsection_multiplicities_1piBG.dat', sigma_1piBG)
       call PrintVals10(includeDIS, raiseFlagVariable, numberOfCalls, &
            'neutrino_total_Xsection_multiplicities_DIS.dat', sigma_DIS)
       call PrintVals10((include2p2hQE.or.include2p2hDelta), &
            raiseFlagVariable, numberOfCalls, &
            'neutrino_total_Xsection_multiplicities_2p2h.dat', sigma_2p2h)

       call PrintVals10(.true., raiseFlagVariable, numberOfCalls, &
            'neutrino_total_Xsection_multiplicities_gen0.dat', sigma_gen0)
       call PrintVals10(.true., raiseFlagVariable, numberOfCalls, &
            'neutrino_total_Xsection_multiplicities_gen1.dat', sigma_gen1)
       call PrintVals10(.true., raiseFlagVariable, numberOfCalls, &
            'neutrino_total_Xsection_multiplicities_gen2.dat', sigma_gen2)
       call PrintVals10(.true., raiseFlagVariable, numberOfCalls, &
            'neutrino_total_Xsection_multiplicities_gen3ormore.dat', &
            sigma_gen3ormore)
    end if



 ! extra channels for ZeroPion analysis

    if (ZeroPion_analysis) then
       call event_sigma_0pions(events,sigma_0pions,initHists,numberOfCalls, &
          & identifier=raiseFlagVariable)
       if (includeQE) call event_sigma_0pions(events_QE,sigma_QE_0pions, &
            initHists,numberOfCalls,identifier=raiseFlagVariable)
       if (includeDELTA) call event_sigma_0pions(events_Delta,sigma_Delta_0pions, &
            initHists,numberOfCalls,identifier=raiseFlagVariable)
       if (includeRES) call event_sigma_0pions(events_highRES,sigma_highRES_0pions, &
            initHists,numberOfCalls,identifier=raiseFlagVariable)
       if (includeDIS) call event_sigma_0pions(events_DIS,sigma_DIS_0pions, &
            initHists,numberOfCalls,identifier=raiseFlagVariable)
       if (include1pi) call event_sigma_0pions(events_1piBG,sigma_1piBG_0pions, &
            initHists,numberOfCalls,identifier=raiseFlagVariable)
       if (include2pi) call event_sigma_0pions(events_2piBG,sigma_2piBG_0pions, &
            initHists,numberOfCalls,identifier=raiseFlagVariable)
       if (include2p2hQE .or. include2p2hDelta) call event_sigma_0pions(events_2p2h, &
          & sigma_2p2h_0pions,initHists,numberOfCalls,identifier=raiseFlagVariable)

       call PrintVals10(.true., raiseFlagVariable, numberOfCalls, &
            'neutrino_0pions.dat', sigma_0pions)
       call PrintVals10(includeQE, raiseFlagVariable, numberOfCalls, &
            'neutrino_0pions_QE.dat', sigma_QE_0pions)
       call PrintVals10(includeDelta, raiseFlagVariable, numberOfCalls, &
            'neutrino_0pions_Delta.dat', sigma_Delta_0pions)
       call PrintVals10(includeRes, raiseFlagVariable, numberOfCalls, &
            'neutrino_0pions_highRES.dat', sigma_highRES_0pions)
       call PrintVals10(include1pi, raiseFlagVariable, numberOfCalls, &
            'neutrino_0pions_1piBG.dat', sigma_1piBG_0pions)
       call PrintVals10(include2pi, raiseFlagVariable, numberOfCalls, &
            'neutrino_0pions_2piBG.dat', sigma_2piBG_0pions)
       call PrintVals10(includeDIS, raiseFlagVariable, numberOfCalls, &
            'neutrino_0pions_DIS.dat', sigma_DIS_0pions)
       call PrintVals10((include2p2hQE .or. include2p2hDelta), raiseFlagVariable, numberOfCalls,&
            'neutrino_0pions_2p2h.dat', sigma_2p2h_0pions)
    end if


    if (calorimetric_analysis) then
       call event_hadronicEnergy(events,numin,numax,nubin,Ehad_versusNu_hist, &
            & dSigdNu_hist, dSigdEhad_hist, &
            Enumin, Enumax, Enubin, Enurestored_versusEnu, dSigdEnu, dSigdEnurestored)
    end if



    ! Now differential cross sections



    !**************************************************************************
    !****o* neutrinoAnalysis/diff_ZZZ_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat
    ! NAME
    ! file diff_ZZZ_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat
    !
    ! with:
    ! * ZZZ denotes the origin (the first interaction vertex) of the event
    !   (see description in neutrino.EprimeCostplaneSX.ZZZ.dat).
    !
    !   If ZZZ is missing the total cross section
    !   (sum over all primary events) is given, otherwise ZZZ=Delta, DIS, ....
    ! * XXX=000, 001, 002 --- the first, second, third and so on values of
    !   the "raised variable"
    !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode,
    !   Elepton for nuXsectionMode=2=dSigmadQsdElepton  and so on)
    ! * HADRON = pion, N, K, K~, Lambda, SigmaResonance, eta
    !   only those are shown, for which the switches in the namelist
    !   "detailed_diff" are set to .true.
    !   possible hadrons are contained in the field particleIDs defined in AnaEvent.f90
    ! * CHARGE = charge of the outgoing hadron
    !
    ! PURPOSE
    ! The file is produced in runs with:
    ! * eventtype=5=neutrino if switch "detailed_diff_output" in
    !   namelist "neutrinoAnalysis" the is set to .true.
    ! * eventtype=3=LowPhoto if switch "dE_switch" in
    !   namelist "LowElePhoto_Analysis" is set to .true.
    !
    ! The file shows the cross sections for ___1-HADRON___  final state
    ! versus kinetic energy of the outgoing hadron
    ! (one HADRON of a given CHARGE and no HADRONs with same flavor, but different charges;
    !  there could be additional hadrons with different flavor)
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV
    ! * For eventtype=5 and process_ID=EM: nanobarns=10^{-33}cm^2/GeV
    ! * For eventtype=3: microbarns=10^{-30}cm^2/GeV
    !
    ! Columns:
    ! * #1: kinetic energy of the outgoing HADRON of given CHARGE [GeV]
    ! * #2: dsi/dEkin  xsec
    ! * #3: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/diff_ZZZ_XXX_dSigma_dEkin_HADRON_charge_CHARGE_1X.dat
    ! NAME
    ! file diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE_1X.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat but
    ! for ___1-HADRON-X___  final state  (one HADRON of a given CHARGE and
    ! any number of  HADRONs of different charges)
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/diff_ZZZ_XXX_dSigma_dEkin_HADRON_charge_CHARGE_2X.dat
    ! NAME
    ! file diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE_2X.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat but
    ! for ___2-HADRON-X___  final state  (two HADRONs of a given CHARGE and
    ! any number of  HADRONs of different charges)
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/diff_ZZZ_XXX_dSigma_dEkin_HADRON_charge_CHARGE_MULTI.dat
    ! NAME
    ! file diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE_MULTI.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat but
    ! for ___MULTI-HADRON___  final state  (at least one HADRON of a given
    ! CHARGE and  any number of  HADRONs of different charges)
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/diff_ZZZ_XXX_dSigma_dEcostheta_HADRON_charge_CHARGE.dat
    ! NAME
    ! file diff_XXX_dSigma_dEcostheta_HADRON_charge_CHARGE.dat
    !
    ! file diff_XXX_dSigma_dEcostheta_HADRON_charge_CHARGE_1X.dat
    !
    ! file diff_XXX_dSigma_dEcostheta_HADRON_charge_CHARGE_2X.dat
    !
    ! file diff_XXX_dSigma_dEcostheta_HADRON_charge_CHARGE_MULTI.dat)
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat but
    ! dsigma/d(E(1-cosTheta))
    !
    ! The file is produced in runs with:
    ! * eventtype=5=neutrino if switch "detailed_diff_output"
    !   in namelist "neutrinoAnalysis" the is set to .true.
    !
    ! The file shows the cross sections dsigma/d(E(1-cosTheta)) versus
    ! E*(1-costheta), where E is the energy (full energy, not kinetic) of the
    ! outgoing hadron, costheta its polar scattering angle (recall here that
    ! in neutrino runs neutrinos are moving along z-direction)
    !
    ! Columns:
    ! * #1: E*(1-cosTheta) [GeV]
    ! * #2: dsi/d(E(1-costheta))  xsec
    ! * #3: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/diff_ZZZ_XXX_dSigma_dTheta_HADRON_charge_CHARGE.dat
    ! NAME
    ! file diff_XXX_dSigma_dTheta_HADRON_charge_CHARGE.dat
    !
    ! file diff_XXX_dSigma_dTheta_HADRON_charge_CHARGE_1X.dat
    !
    ! file diff_XXX_dSigma_dTheta_HADRON_charge_CHARGE_2X.dat
    !
    ! file diff_XXX_dSigma_dTheta_HADRON_charge_CHARGE_MULTI.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat but
    ! dsigma/dTheta
    !
    ! The file is produced in runs with:
    ! * eventtype=5=neutrino if switch "detailed_diff_output"
    !   in namelist "neutrinoAnalysis" the is set to .true.
    ! * eventtype=3=LowPhoto if switch "dTheta_switch" in
    !   namelist "LowElePhoto_Analysis" is set to .true.
    !
    ! The file shows the cross sections dsigma/dTheta versus Theta,
    ! where Theta is the polar scattering angle (in radians) of the outgoing
    ! hadron (recall here that in neutrino runs neutrinos are moving along
    ! z-direction)
    !
    ! Columns:
    ! * #1: Theta [radians]
    ! * #2: dsi/dTheta  xsec
    ! * #3: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************



    if (detailed_diff_output) then

       ! Here we fill the histograms with the appropriate event types

       call event_dSigma_dE(events,EkinMin,EkinMax,dEkin,&
            & 'diff_'//trim(intToChar(numberofFinals)),numberOfCalls, &
            & dE_hists,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_Multi,hists1X=dE_hists_1X,hists2X=dE_hists_2X)
       call event_dSigma_dOmega(events,dTheta,dPhi,&
            & 'diff_'//trim(intToChar(numberofFinals)),numberOfCalls,  &
            & dTheta_hists, dPhi_hists, dOmega_hists,initHists,sameFileNameIn=.true.)
       call event_dSigma_dEcostheta(events,EcosthetaMin,EcosthetaMax,dEcostheta,'diff_'// &
            & trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists,dcostheta_hists,&
            & dEcostheta_hists_MULTI,dcostheta_hists_MULTI, initHists,sameFileNameIn=.true.)

       if (includeQE) then
          call event_dSigma_dE(events_QE,EkinMin,EkinMax,dEkin,&
               & 'diff_QE_'//trim(intToChar(numberofFinals)),numberOfCalls, &
               & dE_hists_QE,initHists,sameFileNameIn=.true., histsMulti=dE_hists_QE_Multi,&
               & hists1X=dE_hists_QE_1X,hists2X=dE_hists_QE_2X)
          call event_dSigma_dOmega(events_QE,dTheta,dPhi,&
               & 'diff_QE_'//trim(intToChar(numberofFinals)),numberOfCalls,  &
               & dTheta_hists_QE, dPhi_hists_QE, dOmega_hists_QE,initHists,sameFileNameIn=.true.)
          call event_dSigma_dEcostheta(events_QE,EcosthetaMin,EcosthetaMax,dEcostheta, &
               &'diff_QE_'//trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists_QE, &
               & dcostheta_hists_QE,dEcostheta_hists_QE_MULTI,dcostheta_hists_QE_MULTI, &
               & initHists,sameFileNameIn=.true.)
       end if

       if (includeDELTA) then
          call event_dSigma_dE(events_Delta,EkinMin,EkinMax,dEkin,&
               & 'diff_Delta_'//trim(intToChar(numberofFinals)),numberOfCalls, &
               & dE_hists_Delta,initHists,sameFileNameIn=.true., histsMulti=dE_hists_Delta_Multi,&
               & hists1X=dE_hists_Delta_1X,hists2X=dE_hists_Delta_2X)
          call event_dSigma_dOmega(events_Delta,dTheta,dPhi,&
               & 'diff_Delta_'//trim(intToChar(numberofFinals)),numberOfCalls,  &
               & dTheta_hists_Delta, dPhi_hists_Delta, dOmega_hists_Delta,initHists, &
               & sameFileNameIn=.true.)
          call event_dSigma_dEcostheta(events_Delta,EcosthetaMin,EcosthetaMax,dEcostheta, &
               & 'diff_Delta_'// &
               & trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists_Delta, &
               & dcostheta_hists_Delta,&
               & dEcostheta_hists_Delta_MULTI,dcostheta_hists_Delta_MULTI,initHists, &
               & sameFileNameIn=.true.)
       end if

       if (includeRES) then
          call event_dSigma_dE(events_highRES,EkinMin,EkinMax,dEkin,&
               & 'diff_highRES_'//trim(intToChar(numberofFinals)),numberOfCalls,&
               & dE_hists_highRES,initHists,sameFileNameIn=.true., &
               & histsMulti=dE_hists_highRES_Multi,&
               & hists1X=dE_hists_highRES_1X,hists2X=dE_hists_highRES_2X)
          call event_dSigma_dOmega(events_highRES,dTheta,dPhi,&
               & 'diff_highRES_'//trim(intToChar(numberofFinals)),numberOfCalls,&
               & dTheta_hists_highRES, dPhi_hists_highRES, dOmega_hists_highRES,initHists, &
               & sameFileNameIn=.true.)
          call event_dSigma_dEcostheta(events_highRES,EcosthetaMin,EcosthetaMax,dEcostheta,&
               & 'diff_highRES_'// &
               & trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists_highRES, &
               & dcostheta_hists_highRES,&
               & dEcostheta_hists_highRES_MULTI,dcostheta_hists_highRES_MULTI,initHists, &
               & sameFileNameIn=.true.)
       end if

       if (include1pi) then
          call event_dSigma_dE(events_1piBG,EkinMin,EkinMax,dEkin,&
               & 'diff_1piBG_'//trim(intToChar(numberofFinals)),&
               & numberOfCalls,dE_hists_1piBG,initHists,sameFileNameIn=.true., &
               & histsMulti=dE_hists_1piBG_Multi,&
               & hists1X=dE_hists_1piBG_1X,hists2X=dE_hists_1piBG_2X)
          call event_dSigma_dOmega(events_1piBG,dTheta,dPhi,'diff_1piBG_'// &
               & trim(intToChar(numberofFinals)), &
               & numberOfCalls,dTheta_hists_1piBG, dPhi_hists_1piBG, dOmega_hists_1piBG, &
               & initHists,sameFileNameIn=.true.)
          call event_dSigma_dEcostheta(events_1piBG,EcosthetaMin,EcosthetaMax,dEcostheta, &
               & 'diff_1piBG_'// &
               & trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists_1piBG, &
               & dcostheta_hists_1piBG,&
               & dEcostheta_hists_1piBG_MULTI,dcostheta_hists_1piBG_MULTI,initHists, &
               & sameFileNameIn=.true.)
       end if

       if (include2pi) then
          call event_dSigma_dE(events_2piBG,EkinMin,EkinMax,dEkin,&
               & 'diff_2piBG_'//trim(intToChar(numberofFinals)),&
               & numberOfCalls,dE_hists_2piBG,initHists,sameFileNameIn=.true., &
               &  histsMulti=dE_hists_2piBG_Multi,&
               & hists1X=dE_hists_2piBG_1X,hists2X=dE_hists_2piBG_2X)
          call event_dSigma_dOmega(events_2piBG,dTheta,dPhi,'diff_2piBG_'//&
               &trim(intToChar(numberofFinals)), &
               & numberOfCalls,dTheta_hists_2piBG, dPhi_hists_2piBG, dOmega_hists_2piBG,&
               & initHists,sameFileNameIn=.true.)
          call event_dSigma_dEcostheta(events_2piBG,EcosthetaMin,EcosthetaMax,dEcostheta,&
               & 'diff_2piBG_'// &
               & trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists_2piBG, &
               & dcostheta_hists_2piBG,&
               & dEcostheta_hists_2piBG_MULTI,dcostheta_hists_2piBG_MULTI,initHists,&
               & sameFileNameIn=.true.)
       end if

       if (includeDIS) then
          call event_dSigma_dE(events_DIS,EkinMin,EkinMax,dEkin,&
               & 'diff_DIS_'//trim(intToChar(numberofFinals)), &
               & numberOfCalls,dE_hists_DIS,initHists,sameFileNameIn=.true., &
               & histsMulti=dE_hists_DIS_Multi,&
               & hists1X=dE_hists_DIS_1X,hists2X=dE_hists_DIS_2X)
          call event_dSigma_dEcostheta(events_DIS,EcosthetaMin,EcosthetaMax,dEcostheta,&
               & 'diff_DIS_'// &
               & trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists_DIS, &
               & dcostheta_hists_DIS,&
               & dEcostheta_hists_DIS_MULTI,dcostheta_hists_DIS_MULTI,initHists,&
               & sameFileNameIn=.true.)
       end if


       if (include2p2hQE .or. include2p2hDelta) then
          call event_dSigma_dE(events_2p2h,EkinMin,EkinMax,dEkin,&
               & 'diff_2p2h_'//trim(intToChar(numberofFinals)), &
               & numberOfCalls,dE_hists_2p2h,initHists,sameFileNameIn=.true., &
               & histsMulti=dE_hists_2p2h_Multi,&
               & hists1X=dE_hists_2p2h_1X,hists2X=dE_hists_2p2h_2X)
       end if



       call event_dSigma_dE(events_gen0,EkinMin,EkinMax,dEkin,&
            & 'diff_gen0_'//trim(intToChar(numberofFinals)), &
            & numberOfCalls,dE_hists_gen0,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_gen0_Multi,&
            & hists1X=dE_hists_gen0_1X,hists2X=dE_hists_gen0_2X)
       call event_dSigma_dOmega(events_gen0,dTheta,dPhi,'diff_gen0_'//&
            & trim(intToChar(numberofFinals)), &
            & numberOfCalls,dTheta_hists_gen0, dPhi_hists_gen0, dOmega_hists_gen0,initHists, &
            & sameFileNameIn=.true.)
       call event_dSigma_dEcostheta(events_gen0,EcosthetaMin,EcosthetaMax,dEcostheta, &
            & 'diff_gen0_'// &
            & trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists_gen0, &
            & dcostheta_hists_gen0,&
            & dEcostheta_hists_gen0_MULTI,dcostheta_hists_gen0_MULTI,initHists, &
            & sameFileNameIn=.true.)
    end if




    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_dSigma_dEkin_lepton_PPP.dat
    ! NAME
    ! file diff_XXX_dSigma_dEkin_lepton_PPP.dat
    !
    ! with:
    ! * XXX=000, 001, 002 --- the first, second, third and so on values of
    !   the "raised variable"
    !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode,
    !   Elepton for nuXsectionMode=2=dSigmadQsdElepton  and so on)
    ! * PPP = no_pi, p_Xn_no_pi, piplus, pi0, ....
    !   (see namelist nl_specificEvent for the full list)
    !   standing for the specific final state under consideration
    !
    ! PURPOSE
    ! The file is produced in runs with:
    ! * eventtype=5=neutrino if switch "specificEvent_Analysis" in namelist
    !   "neutrinoAnalysis" the is set to .true.
    !   and specific final states in the namelist "nl_specificEvent" are
    !   set to .true.
    !
    ! The file shows the cross sections for a specific final state versus kinetic energy of the
    ! outgoing lepton
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV
    ! * For eventtype=5 and process_ID=EM: nanobarns=10^{-33}cm^2/GeV
    !
    ! Columns:
    ! * #1: kinetic energy of the outgoing lepton  [GeV]
    ! * #2: dsi/dEkin  xsec
    ! * #3: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_dSigma_dQ2_lepton_PPP.dat
    ! NAME
    ! file diff_XXX_dSigma_dQ2_lepton_PPP.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_lepton_PPP.dat but for dsigma/dQ2
    !
    ! The file shows the cross sections for a specific final state versus Q2
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV^2
    ! * For eventtype=5 and process_ID=EM: nanobarns=10^{-33}cm^2/GeV^2
    !
    ! Columns:
    ! * #1: Q2 [GeV^2]  squared transfer momentum (Q2 = -q_mu \cdot q^\mu)
    ! * #2: dsi/dQ2  xsec
    ! * #3: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_dSigma_dQ2p_lepton_PPP.dat
    ! NAME
    ! file diff_XXX_dSigma_dQ2p_lepton_PPP.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_lepton_PPP.dat but for dsigma/dQ2
    !
    ! The file shows the cross sections for a final state with 1 mu, 0 pi ,
    ! and (at least) 1 p versus Q2
    ! here Q2 is calculated from the kinematics of the outgoing leading
    ! proton: Q2prot = (mProt - epsB)**2 - mProt**2  &
    !            & + 2*(mNeut - epsB)*(Tp + mProt - mNeut + epsB)
    ! This distribution can be converted into a kinetic energy distribution
    ! of 0pion events
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV^2
    ! * For eventtype=5 and process_ID=EM: nanobarns=10^{-33}cm^2/GeV^2
    !
    ! Columns:
    ! * #1: Q2 [GeV^2]  squared transfer momentum, from proton kinematics
    ! * #2: dsi/dQ2  xsec
    ! * #3: number of events contributed  (only for internal use,
    !   you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_dSigma_dcos_lepton_PPP.dat
    ! NAME
    ! file diff_XXX_dSigma_dcos_lepton_PPP.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_lepton_PPP.dat but for
    ! dsigma/dcos(theta_l)
    !
    ! The file shows the cross sections for a specific final state
    ! versus cos of the scattering angle of the outgoin lepton
    ! (with respect to neutrino direction)
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2
    ! * For eventtype=5 and process_ID=EM: nanobarns=10^{-33}cm^2
    !
    ! Columns:
    ! * #1: cos(theta_l)  cos of the scattering angle of the outgoing lepton
    ! * #2: dsi/dcos(theta_l)  xsec
    ! * #3: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_d2Sigma_dEdcost_lepton_no_pi.dat
    ! NAME
    ! file diff_XXX_dSigma_dEdcost_lepton_no_pi.dat
    !
    ! PURPOSE
    ! Double differential cross section d2sigma/(dEkin dcostheta) for outgoing
    ! lepton
    !
    ! The file shows the cross section for the outgoing lepton
    ! versus cos of the scattering angle of the outgoin lepton
    ! (with respect to neutrino direction) and its kinetic energy
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2
    ! * For eventtype=5 and process_ID=EM: nanobarns=10^{-33}cm^2
    !
    ! Columns:
    ! * #1: Ekin of the outgoing lepton (in GeV)
    ! * #2: cos(theta_l)  cos of the scattering angle of the outgoing lepton
    ! * #3: dsi/dcos(theta_l)  xsec
    ! * #4: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #5: xsec-statistical-error
    !**************************************************************************

 !  Calculate differential cross sections dsigma/dx for lepton with specific events
 !  Now loop over specific events, such as 0 pion, 1p + Xn, ...



    if (specificEvent_Analysis) then
        do m=1, max_speEvent           ! loop over specific events
          if (.not.includeSpeEvent(m)) cycle
!
          do iHist=0, 0         ! loop over first event types, iHist=0: total
!
            if (.not.includeHist(iHist)) cycle

            call event_dSigma_dLeptonVariables(events,EkinMin_lepton, &
                  & EkinMax_lepton,dEkin_lepton, &
                  & cost_min,cost_max,delta_cost, &
                  & 'diff_'//trim(intToChar(iHist)),numberOfCalls, &
                  & dEnu_hist(m,iHist), dElepton_hist(m,iHist),   &
                  & dcoslepton_hist(m,iHist), &
                  & dQ2lepton_hist(m,iHist),dQ2plepton_hist(m,iHist), &
                  & dSigmaMC_EprimeCost_0pi(m,iHist),&
                  & initHists, sameFileNameIn=.true., specificEvent=m)

          end do   !iHist
        end do  !m
    end if





    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_dSigma_dW_nucleon_pion_charge_CHARGE.dat
    ! NAME
    ! file diff_XXX_dSigma_dW_nucleon_pion_charge_CHARGE.dat
    !
    ! PURPOSE
    ! Similar to diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat
    ! (see notations XXX and CHARGE there) but
    ! gives pion-nucleon invariant mass (W) distribution
    ! (true invariant mass with the sum of the 4-momenta of the outgoing
    ! particles, W^2=(p_1out+p_2out)
    ! as opposed to widely used W2=mN2+2*mN*nu-Q2 from lepton kinematics )
    ! for events with ___1pion___  AND  ___1nucleon___ in the final state
    !
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_dSigma_dW_muon_nucleon_charge_CHARGE.dat
    ! NAME
    ! file diff_XXX_dSigma_dW_muon_nucleon_charge_CHARGE.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dW_nucleon_pion_charge_CHARGE.dat but
    ! for the  muon-nucleon invariant mass distribution
    !
    ! NOTES
    ! Note that CHARGE is still a pion charge
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_dSigma_dW_muon_pion_charge_CHARGE.dat
    ! NAME
    ! file diff_XXX_dSigma_dW_muon_pion_charge_CHARGE.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dW_nucleon_pion_charge_CHARGE.dat but
    ! for the  muon-pion invariant mass distribution
    !
    ! NOTES
    ! Note  that CHARGE is still a pion charge
    !**************************************************************************



    if (include_W_dist) then
       call event_dSigma_dInvMass(events,'diff_'//trim(intToChar(numberofFinals)),numberOfCalls, &
            & dW_Npi,Wmin_Npi,Wmax_Npi,dW_Npi_hists, dW_mupi,Wmin_mupi,Wmax_mupi,dW_mupi_hists, &
            & dW_muN,Wmin_muN,Wmax_muN,dW_muN_hists, initHists ,sameFileNameIn=.true.)
    end if


    ! kinetic energy distributions of nucleons in events with 0 pions
    if (detailed_diff_output .and. ZeroPion_analysis) then

       ! Here we initialize the histograms
       call event_dSigma_dE_0pions(events,EkinMin,EkinMax,dEkin,&
            & 'diff_'//trim(intToChar(numberofFinals)),numberOfCalls,&
            & dE_hists_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_Multi_0pions, &
            & hists1X=dE_hists_1X_0pions, hists2X=dE_hists_2X_0pions)

       if (includeQE) &
            & call event_dSigma_dE_0pions(events_QE,EkinMin,EkinMax,dEkin,&
            & 'diff_QE_'//trim(intToChar(numberofFinals)),&
            & numberOfCalls, dE_hists_QE_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_QE_Multi_0pions, &
            & hists1X=dE_hists_QE_1X_0pions, hists2X=dE_hists_QE_2X_0pions)

       if (includeDELTA) &
            & call event_dSigma_dE_0pions(events_Delta,EkinMin,EkinMax,dEkin,&
            & 'diff_Delta_'//trim(intToChar(numberofFinals)),&
            & numberOfCalls ,dE_hists_Delta_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_Delta_Multi_0pions, &
            & hists1X=dE_hists_Delta_1X_0pions, hists2X=dE_hists_Delta_2X_0pions)

       if (includeRES) &
            & call event_dSigma_dE_0pions(events_highRES,EkinMin,EkinMax,dEkin,&
            & 'diff_highRES_'//trim(intToChar(numberofFinals)),&
            & numberOfCalls, dE_hists_highRES_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_highRES_Multi_0pions, &
            & hists1X=dE_hists_highRES_1X_0pions, hists2X=dE_hists_highRES_2X_0pions)

       if (includeDIS) &
            & call event_dSigma_dE_0pions(events_DIS,EkinMin,EkinMax,dEkin,&
            & 'diff_DIS_'//trim(intToChar(numberofFinals)),&
            & numberOfCalls, dE_hists_DIS_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_DIS_Multi_0pions, &
            & hists1X=dE_hists_DIS_1X_0pions, hists2X=dE_hists_DIS_2X_0pions)

       if (include1pi) &
            & call event_dSigma_dE_0pions(events_1piBG,EkinMin,EkinMax,dEkin,&
            & 'diff_1piBG_'//trim(intToChar(numberofFinals)),&
            & numberOfCalls, dE_hists_1piBG_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_1piBG_Multi_0pions, &
            & hists1X=dE_hists_1piBG_1X_0pions, hists2X=dE_hists_1piBG_2X_0pions)

       if (include2pi) &
            & call event_dSigma_dE_0pions(events_2piBG,EkinMin,EkinMax,dEkin,&
            & 'diff_2piBG_'//trim(intToChar(numberofFinals)),&
            & numberOfCalls, dE_hists_2piBG_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_2piBG_Multi_0pions, &
            & hists1X=dE_hists_2piBG_1X_0pions, hists2X=dE_hists_2piBG_2X_0pions)

       if (include2p2hQE .or. include2p2hDelta) &
            & call event_dSigma_dE_0pions(events_2p2h,EkinMin,EkinMax,dEkin,&
            & 'diff_2p2h_'//trim(intToChar(numberofFinals)),&
            & numberOfCalls, dE_hists_2p2h_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_2p2h_Multi_0pions, &
            & hists1X=dE_hists_2p2h_1X_0pions, hists2X=dE_hists_2p2h_2X_0pions)

    end if

    !multiplicity output
    write(Prefix_MultAna,'(f8.4)') raiseFlagVariable
    do i=firstEvents(1),firstEvents(2)
       call Multiplicity_AddEvent(events(i))
    end do


    call Multiplicity_Write(Prefix_MultAna)

    if (outputEvents) then
       write(*,*) 'Writing events to file'

       !***********************************************************************
       !****o* neutrinoAnalysis/FinalEvents.dat
       ! NAME
       ! file FinalEvents.dat
       ! PURPOSE
       ! The file contains the positions and four-vectors of all outgoing particles from an event
       ! For a cross-section construction all of the events have to be weighted with the
       ! 'perweight' values in col. 5
       ! The output of this file is switched on by the switch 'outputEvents'.
       !
       ! Columns:
       ! * #1: run number (from 1 to num_runs_SameEnergy)
       ! * #2: event number (from 1 to ... less then target_A*numEnsembles)
       ! * #3: ID of the outgoing particle
       ! * #4: charge of the outgoing particle
       ! * #5: perweight of the event
       ! * #6-#8: position of the outgoing particle (set to 0 for outgoing lepton)
       ! * #9-#12: 4-momentum of the outgoing particle in GeV
       ! * #13: history of the outgoing particle (see history.f90 for the definition)
       ! * #14: production_ID (type of the first event: 1=QE, 2=Delta, 34=DIS)
       ! * #15: incoming neutrino energy
       !
       ! NOTES
       ! There is always in each event a particle with weight 0. This is the
       ! nucleon on which the initial interaction happened.
       !
       ! In the case of an initial 2p2h process the second initial-state nucleon
       ! is not written out. It is chosen to be at the same place as the first
       ! initial nucleon (the one with weight=0), with a randomly chosen
       ! momentum in the Fermi-sea.
       !
       ! For large target_A, numEnsembles and num_runs_SameEnergy, this file
       ! may become very large, e.g.:
       ! * target_A=12, 1000ens x 20runs, QE and Delta results in 41 Mb
       ! * target_A=56, 2000ens x 5runs, QE,Delta,highRES,1piBG,DIS results in 500 Mb
       !***********************************************************************
       filename='FinalEvents.dat'
       call event_dump(numberOfCalls,events,filename,writeNeutrinoProdID=.true.)

       ! for Sigma_MC run huge files, and they are not necessary
       ! because the same info can be obtained from 'FinalEvents.dat' using column 14

       !filename='FinalEvents_QE.dat'
       !call event_dump(numberOfCalls,events_QE,filename,writeNeutrinoProdID=.true.)

       !filename='FinalEvents_Delta.dat'
       !call event_dump(numberOfCalls,events_Delta,filename,writeNeutrinoProdID=.true.)

       !filename='FinalEvents_highRES.dat'
       !call event_dump(numberOfCalls,events_highRES,filename,writeNeutrinoProdID=.true.)

       !filename='FinalEvents_1piBG.dat'
       !call event_dump(numberOfCalls,events_1piBG,filename,writeNeutrinoProdID=.true.)

       !filename='FinalEvents_DIS.dat'
       !call event_dump(numberOfCalls,events_DIS,filename,writeNeutrinoProdID=.true.)
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








    ! RECONSTRUCTION

    if (reconstruct_neutrino_energy .and. specificEvent_Analysis) then

       do j=lbound(events,dim=1),ubound(events,dim=1)
          if (.not.neutrinoProdInfo_Get(j,prod_id,perweight,lepIn_mom,lep_mom,  &
            & boson_mom,nuc_mom,Chrg_Nuc)) then
             write(*,*) j,prod_id,perweight
             call TRACEBACK('error in getting perweight, stop')
          end if


          do m=1, max_speEvent
             if (.not.includeSpeEvent(m)) cycle
             if (.not.IfPass_SpecificEvent(m,events(j))) cycle
             Q2real=-abs4Sq(boson_mom)
             Enureal=lep_mom(0)+boson_mom(0)

             select case (m)
             case (1) ! events with 0 pions,
                Q2rec=CCQE_recQs(lep_mom)           !CCQE QE-like reconstruction
                Enurec=CCQE_recEnergy(lep_mom)

             case (2,7:16) ! events with 0 pions and 1 proton and X neutrons,
                Q2rec=K2K_recQs(lep_mom)                 ! K2K QE-like reconstruction
                Enurec=K2K_recEnergy(lep_mom)

             case (3,4) ! events with 1 pi+,  the same for 1pi0
                Q2rec=CCQE_recQs_Delta(lep_mom)         ! CCQE energy reconstruction assuming
                Enurec=CCQE_recEnergy_Delta(lep_mom)    ! Delta mass

             case (5,6) ! events with at least 1 pi0 (as K2K),   the same for at least 1pi+
                Q2rec=K2K_recQs(lep_mom,1.483)               ! K2K  energy reconstruction
                Enurec=K2K_recEnergy(lep_mom,1.483)          ! assuming DIS with W=1.483

             case default
             end select


             call AddHist(H_Q2_real(m,0),H_Q2_real(m,K2Hist(prod_id)), Q2real,perweight/ &
                & float(num_runs_sameEnergy))
             call AddHist(H_Q2_rec(m,0), H_Q2_rec(m,K2Hist(prod_id)),  Q2rec, perweight/ &
                & float(num_runs_sameEnergy))

             if (nuEXP>0) call AddHist(H_enu_real(m,0),H_enu_real(m,K2Hist(prod_id)), &
                  & Enureal,perweight/float(num_runs_sameEnergy))
             call AddHist(H_enu_rec(m,0), H_enu_rec(m,K2Hist(prod_id)),  Enurec, perweight/ &
                & float(num_runs_sameEnergy))

             call AddHist2D(H_Q2_rec_versus_real(m,0),H_Q2_rec_versus_real(m,K2Hist(prod_id)), &
                  & (/Q2real,Q2rec/),  perweight/float(num_runs_sameEnergy))
             if (nuEXP>0) call AddHist2D(H_enu_rec_versus_real(m,0), &
                  & H_enu_rec_versus_real(m,K2Hist(prod_id)), &
                  & (/Enureal,Enurec/),perweight/float(num_runs_sameEnergy))

             !! add oscillated muon disappearance  H_Enu_rec_muonDisappear
             !! add oscillated electron appearance for various CP violation phases
             if (OSC(nuEXP)) then
                L=OSCLENGTH(nuEXP)
                call oscillationProbability(Enureal,L,0.,Posc_mumu,Posc_mue,Posc_mue_max, &
                   & Posc_mue_antimax)

                call AddHist(Oscmumu_enu_real(m,0),Oscmumu_enu_real(m,K2Hist(prod_id)), &
                     &  Enureal,Posc_mumu*perweight/float(num_runs_sameEnergy))
                call AddHist(Oscmumu_enu_rec(m,0), Oscmumu_enu_rec(m,K2Hist(prod_id)),  &
                     & Enurec, Posc_mumu*perweight/float(num_runs_sameEnergy))

                call AddHist(Oscmue_enu_real(m,0),Oscmue_enu_real(m,K2Hist(prod_id)), &
                     & Enureal,Posc_mue*perweight/float(num_runs_sameEnergy))
                call AddHist(Oscmue_enu_rec(m,0), Oscmue_enu_rec(m,K2Hist(prod_id)),  &
                     & Enurec, Posc_mue*perweight/float(num_runs_sameEnergy))

                call AddHist(Oscmuemax_enu_real(m,0),Oscmuemax_enu_real(m,K2Hist(prod_id)), &
                     & Enureal,Posc_mue_max*perweight/float(num_runs_sameEnergy))
                call AddHist(Oscmuemax_enu_rec(m,0), Oscmuemax_enu_rec(m,K2Hist(prod_id)),  &
                     & Enurec, Posc_mue_max*perweight/float(num_runs_sameEnergy))

                call AddHist(Oscmueantimax_enu_real(m,0), &
                     & Oscmueantimax_enu_real(m,K2Hist(prod_id)),&
                     & Enureal,Posc_mue_antimax*perweight/float(num_runs_sameEnergy))
                call AddHist(Oscmueantimax_enu_rec(m,0), &
                     & Oscmueantimax_enu_rec(m,K2Hist(prod_id)),&
                     & Enurec, Posc_mue_antimax*perweight/float(num_runs_sameEnergy))
             end if

          end do ! m

       end do ! j


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! energy reconstruction output in files

       do m=1, max_speEvent
          if (.not.includeSpeEvent(m)) cycle
          call SpecificEvent_Name(m,filename1)
          do iHist=0, max_Hist
             if (.not.includeHist(iHist)) cycle

             call writeHist(H_Q2_real(m,iHist), file='reconstruction_Q2real_'//trim(filename1)//&
                & "."//trim(intToChar(iHist))//'.dat')
             call writeHist(H_Q2_rec(m,iHist),  file='reconstruction_Q2rec_'//trim(filename1)//&
                & "."//trim(intToChar(iHist))//'.dat')
             if (nuEXP>0) call writeHist(H_enu_real(m,iHist), &
                  & file='reconstruction_Enureal_'//trim(filename1)//"."//trim(intToChar(iHist))//&
                  & '.dat')
             call writeHist(H_enu_rec(m,iHist), file='reconstruction_Enurec_'//trim(filename1)// &
                  & "."//trim(intToChar(iHist))//'.dat')
             call writeHist2D_Gnuplot(H_Q2_rec_versus_real(m,iHist), &
                  & file='reconstruction_Q2_rec_versus_real_'//trim(filename1)//"."// &
                  & trim(intToChar(iHist))//'.dat')
             if (nuEXP>0)  call writeHist2D_Gnuplot(H_enu_rec_versus_real(m,iHist), &
                  & file='reconstruction_Enu_rec_versus_real_'//trim(filename1)//"."// &
                  & trim(intToChar(iHist))//'.dat')
          end do ! iHist
       end do ! m


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! oscillation versus true and reconstructed energy  output in files

       do m=1, max_speEvent
          if (.not.includeSpeEvent(m) .or. m>6) cycle
          call SpecificEvent_Name(m,filename1)
          do iHist=0, max_Hist
             if (.not.includeHist(iHist)) cycle


             if (OSC(nuEXP) ) then

                call writeHist(Oscmumu_enu_real(m,iHist), &
                     & file='oscillations_mumu_real_'//trim(filename1)//"."// &
                     & trim(intToChar(iHist))//'.dat')
                call writeHist(Oscmumu_enu_rec(m,iHist), &
                     & file='oscillations_mumu_rec_'//trim(filename1)//"."// &
                     & trim(intToChar(iHist))//'.dat')
                call writeHist(Oscmuemax_enu_real(m,iHist), &
                     & file='oscillations_mue_max_real_'//trim(filename1)//"."// &
                     & trim(intToChar(iHist))//'.dat')
                call writeHist(Oscmuemax_enu_rec(m,iHist), &
                     & file='oscillations_mue_max_rec_'//trim(filename1)//"."// &
                     & trim(intToChar(iHist))//'.dat')
                call writeHist(Oscmue_enu_real(m,iHist), &
                     & file='oscillations_mue_real_'//trim(filename1)//"."// &
                     & trim(intToChar(iHist))//'.dat')
                call writeHist(Oscmue_enu_rec(m,iHist), &
                     & file='oscillations_mue_rec_'//trim(filename1)//"."// &
                     & trim(intToChar(iHist))//'.dat')
                call writeHist(Oscmueantimax_enu_real(m,iHist), &
                     & file='oscillations_mue_antimax_real_'//trim(filename1)//"."// &
                     & trim(intToChar(iHist))//'.dat')
                call writeHist(Oscmueantimax_enu_rec(m,iHist), &
                     & file='oscillations_mue_antimax_rec_'//trim(filename1)// &
                     & "."//trim(intToChar(iHist))//'.dat')

             end if
          end do ! iHist
       end do ! m

    end if ! reconstruct_neutrino_energy .and. specificEvent_Analysis






    !**************************************************************************
    !****o* neutrinoAnalysis/reconstruction_Enureal_PPP.ZZZ.dat
    ! NAME
    ! file reconstruction_Enureal_PPP.ZZZ.dat
    !
    ! PURPOSE
    ! Similar to reconstruction_Enurec_PPP.ZZZ.dat, but for true neutrino energy
    !
    ! The file shows the event distribution (flux times xsec) versus true
    ! neutrino energy for a specific final state
    ! The file is produced only for the runs with a neutrino flux
    ! (in the namelist "neutrino_induced" nuXsectionMode > 10, nuExp>0)
    !
    ! Columns:
    ! * #1: true neutrino energy  [GeV]
    ! * #2: event distribution: normalized flux times xsec
    ! * #3: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #4: statistical-error of #2
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/reconstruction_Enurec_PPP.ZZZ.dat
    ! NAME
    ! file reconstruction_Enurec_PPP.ZZZ.dat
    !
    ! with:
    ! * PPP = no_pi, p_Xn_no_pi, piplus, pi0, ....
    !   (see namelist nl_specificEvent for the full list)
    !   standing for the specific final state under consideration
    ! * ZZZ=000 - 008 is the origin (the first interaction vertex) of the event
    !   (see description in  neutrino.EprimeCostplaneXS.ZZZ.dat)
    !
    ! PURPOSE
    ! The file shows the event distribution (normalized flux times xsec)
    ! versus reconstructed neutrino energy (see arXiv:1208.3678 [nucl-th])
    ! for a specific final state
    !
    ! Reconstruction method depends on a specific process;
    ! * generally for no_pion events it is based on QE-like kinematics
    ! * and for pion events it is based on on-shell-Delta-creation assumption
    !
    ! The file is produced in runs with:
    ! * eventtype=5=neutrino
    !   if switch "reconstruct_neutrino_energy" is set to .true.
    !   and switch "specificEvent_Analysis" in namelist "neutrinoAnalysis"
    !   is set to .true.
    !   and specific final states in the namelist "nl_specificEvent" are set
    !   to .true.
    !
    ! Units:
    ! * For eventtype=5  and process_ID=CC and NC: 10^{-38} cm^2/GeV
    ! * For eventtype=5  and process_ID=EM: nanobarns=10^{-33}cm^2/GeV
    ! * All x-sec per particle (1/A)
    !
    ! Columns:
    ! * #1: reconstructed neutrino energy  [GeV]
    ! * #2: event distribution: flux times xsec)
    ! * #3: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #4: statistical-error of #2
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/reconstruction_Enu_rec_versus_real_PPP.ZZZ.dat
    ! NAME
    ! file reconstruction_Enu_rec_versus_real_PPP.ZZZ.dat
    !
    ! PURPOSE
    ! Similar to reconstruction_Enurec_PPP.ZZZ.dat, but double-differential
    !
    ! The file shows the 2-D density of the flux times xsec
    ! (see arXiv:1208.3678 [nucl-th])
    ! versus true and reconstructed neutrino energies for a specific final state
    ! The file is produced only for the runs with a neutrino flux
    ! (in the namelist "neutrino_induced" nuXsectionMode > 10, nuExp>0)
    !
    ! Units:
    ! * For event_type=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV^2
    ! * All xsec per particle (1/A)
    !
    ! Columns:
    ! * #1: true neutrino energy  [GeV]
    ! * #2: reconstructed neutrino energy  [GeV]
    ! * #3: flux-folded xsec
    ! * #4: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #5: xsec-statistical-error
    !
    ! NOTES
    ! up to normalization this is "migration matrix" between true and
    ! reconstructed energies
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/reconstruction_Q2real_PPP.ZZZ.dat
    ! NAME
    ! file reconstruction_Q2real_PPP.ZZZ.dat
    !
    ! PURPOSE
    ! Similar to reconstruction_Enurec_PPP.ZZZ.dat, but for the true Q2
    !
    ! The file contains the flux averaged cross section dsigma/dQ2,
    ! i.e. integral  \int \Phi(E) dsigma/dQ2 (E) dE,
    ! versus true Q2 for a specific final state
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC:  10^{-38} cm^2/GeV^2
    ! * For eventtype=5 and process_ID=EM (one can run it, but makes no
    !   physical sense): nanobarns=10^{-33}cm^2/GeV^2
    ! * All xsec per particle (1/A)
    !
    ! Columns:
    ! * #1: true Q2  [GeV^2] (squared momentum transfer Q2 = - q_mu \cdot q^mu)
    ! * #2: flux-averaged dsigma/dQ2true
    ! * #3: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/reconstruction_Q2rec_PPP.ZZZ.dat
    ! NAME
    ! file reconstruction_Q2rec_PPP.ZZZ.dat
    !
    ! PURPOSE
    ! Similar to reconstruction_Enurec_PPP.ZZZ.dat, but for the reconstructed Q2
    ! (see arXiv:1208.3678 [nucl-th])
    !
    ! The file contains the flux averaged cross section dsigma/dQ2,
    ! i.e. integral  \int \Phi(E) dsigma/dQ2 (E) dE,
    ! versus reconstructed Q2 for a specific final state
    !
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV^2
    ! * For eventtype=5 and process_ID=EM (one can run it, but makes no
    !   physical sense): nanobarns=10^{-33}cm^2/GeV^2
    ! * All xsec per particle (1/A)
    !
    ! Columns:
    ! * #1: reconstructed Q2  [GeV^2] (squared momentum transfer Q2 = - q_mu \cdot q^mu)
    ! * #2: flux-averaged dsigma/dQ2rec
    ! * #3: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/reconstruction_Q2_rec_versus_real_PPP.ZZZ.dat
    ! NAME
    ! file reconstruction_Q2_rec_versus_real_PPP.ZZZ.dat
    !
    ! PURPOSE
    ! Similar to reconstruction_Enurec_PPP.ZZZ.dat, but double-differential
    ! in Q2
    !
    ! The file shows the dsigma/dQ2true dQ2rec distribution
    ! versus true and reconstructed Q^2 for a specific final state
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV^4
    ! * For eventtype=5 and process_ID=EM (one can run it, but makes no
    !   physical sense): nanobarns=10^{-33}cm^2/GeV^4
    !
    ! Columns:
    ! * #1: true Q2  [GeV^2]
    ! * #2: reconstructed Q2  [GeV^2]
    ! * #3: dsigma/dQ2true dQ2rec
    ! * #4: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #5: xsec-statistical-error
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/oscillations_UUU_real_PPP.ZZZ.dat
    ! NAME
    ! file oscillations_UUU_real_PPP.ZZZ.dat
    !
    ! with:
    ! * UUU =  "mumu", "mue", "mue_max", "mue_antimax"
    !   for oscillation of a given experimental flux
    !   for muon_neutrino survival ("mumu"),
    !   mu_e appearance with delta_CP=0 ("mue"),
    !   delta_CP=pi ("mue_max"),
    !   delta_CP=-pi ("mue_antimax")
    ! * PPP = no_pi, p_Xn_no_pi, piplus, pi0, ....
    !   (see namelist nl_specificEvent for the full list)
    !   standing for the specific final state under consideration
    ! * ZZZ=000 - 008 is the origin (the first interaction vertex) of the event
    !   (see description in neutrino.EprimeCostplaneXS.ZZZ.dat)
    !
    ! PURPOSE
    ! Similar to reconstruction_Enureal_PPP.ZZZ.dat, but for oscillated
    ! neutrino flux
    !
    ! Columns:
    ! * #1: true neutrino energy  [GeV]
    ! * #2: flux-folded xsec
    ! * #3: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/oscillations_UUU_rec_PPP.ZZZ.dat
    ! NAME
    ! file oscillations_UUU_rec_PPP.ZZZ.dat
    !
    ! with:
    ! * UUU =  "mumu", "mue", "mue_max", "mue_antimax"
    !   for oscillation of a given experimental flux
    !   for muon_neutrino survival ("mumu"),
    !   mu_e appearance with delta_CP=0 ("mue"),
    !   delta_CP=pi ("mue_max"),
    !   delta_CP=-pi ("mue_antimax")
    ! * PPP = no_pi, p_Xn_no_pi, piplus, pi0, ....
    !   (see namelist nl_specificEvent for the full list)
    !   standing for the specific final state under consideration
    ! * ZZZ=000 - 008 is the origin (the first interaction vertex) of the event
    !   (see description in neutrino.EprimeCostplaneXS.ZZZ.dat)
    !
    ! PURPOSE
    ! Similar to reconstruction_Enurec_PPP.ZZZ.dat, but for oscillated
    ! neutrino flux
    !
    ! Columns:
    ! * #1: reconstructed neutrino energy  [GeV]
    ! * #2: event rate: flux times xsec
    ! * #3: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************

    !##########################################################################
    ! energy-reconstruction/oscillation  ANALYSIS END
    !##########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! (3) Clear the list "events" to clear the memory
    do i=firstEvents(1),firstEvents(2)
       call event_clear(events(i))
       call event_clear(events_QE(i))
       call event_clear(events_Delta(i))
       call event_clear(events_highRES(i))
       call event_clear(events_1piBG(i))
       call event_clear(events_2piBG(i))
       call event_clear(events_DIS(i))
       call event_clear(events_2p2h(i))
       call event_clear(events_gen0(i))
       call event_clear(events_gen1(i))
       call event_clear(events_gen2(i))
       call event_clear(events_gen3ormore(i))
    end do
    deallocate(events)
    deallocate(Events_QE)
    deallocate(Events_Delta)
    deallocate(Events_highRES)
    deallocate(Events_1piBG)
    deallocate(Events_2piBG)
    deallocate(Events_DIS)
    deallocate(Events_2p2h)
    deallocate(Events_gen0)
    deallocate(Events_gen1)
    deallocate(Events_gen2)
    deallocate(Events_gen3ormore)


    if (finalflag) call Multiplicity_Reset


    !    call neutrinoProdInfo_clear
    call neutrinoInfoStorage_clear

    if (finalFlag) then
       numberOfCalls=0
       numberOfFinals=numberOfFinals+1
    end if

    write(*,*) '################### NEUTRINO ANALYSIS FINISHED #######################'

  end subroutine neutrino_Analyze






  !****************************************************************************
  !****f* neutrinoAnalysis/IsBound
  ! NAME
  ! logical function IsBound(part)
  ! PURPOSE
  ! return whether particle is bound or not
  !****************************************************************************
  logical function IsBound(part)
    use particleDefinition
    use potentialModule, only: potential_LRF
    use nucleusDefinition
    use nucleus, only: getTarget
    use vector, only: absVec

    type(particle), intent(in) :: part
    type(tnucleus),pointer :: TargetNuc
    real :: NucRadius

    IsBound=.false.

    targetNuc => getTarget()
    NucRadius=targetNuc%radius

    if (kineticEnergy(part)+potential_LRF(part).lt.0) IsBound=.true.

    if (absVec(part%position).lt.radialScale*NucRadius) IsBound=.true.

    !if flag "inclusiveAnalysis" is set to true, we take all particles as "unbound"
    if (inclusiveAnalysis) IsBound=.false.

  end function IsBound


  !****************************************************************************
  !****f* neutrinoAnalysis/IsBelowThreshold
  ! NAME
  ! logical function IsBelowThreshold(part)
  ! PURPOSE
  ! Returns .true. when a particle is below a given detection threshold
  ! This routine can be used to remove all events with kinetic energies
  ! and angles below specified detection thresholds, if an event is below
  ! threshold it is not added to the full event in subroutine neutrino_Analyze
  !****************************************************************************
  logical function IsBelowThreshold(part)
    use particleDefinition
    use vector, only: absVec
    use degRad_conversion, only: radian

    type(particle), intent(in) :: part
    real :: cos_theta

    IsBelowThreshold=.false.

    cos_theta = part%momentum(3)/absVec(part%momentum(1:3))

    !kinetic energy and angle threshold for nucleons
    if (part%id.eq.1 .and. ((part%momentum(0)-part%mass) .lt.                   &
         & kineticEnergyDetectionThreshold_nucleon .or.                        &
         & cos_theta .lt.                                                      &
         & cos(radian(AngleUpperDetectionThresholdDegrees_nucleon))))          &
         & IsBelowThreshold=.true.

    !kinetic energy and angle threshold for charged pions
    if (part%id.eq.101 .and. part%charge.ne.0.and.((part%momentum(0)-part%mass) &
         & .lt. kineticEnergyDetectionThreshold_chargedpion .or.               &
         & cos_theta .lt.                                                      &
         & cos(radian(AngleUpperDetectionThresholdDegrees_chargedpion))))      &
         & IsBelowThreshold=.true.

    !kinetic energy and angle threshold for neutral pions
    if (part%id.eq.101.and.part%charge.eq.0.and.((part%momentum(0)-part%mass).lt.&
         & kineticEnergyDetectionThreshold_neutralpion .or.                     &
         & cos_theta .lt.                                                       &
         & cos(radian(AngleUpperDetectionThresholdDegrees_neutralpion))))       &
         & IsBelowThreshold=.true.

  end function IsBelowThreshold


  !****************************************************************************
  !****s* neutrinoAnalysis/oscillationProbability
  ! NAME
  ! subroutine oscillationProbability(Enu,L,deltaCP,Posc_mumu,Posc_mue,Posc_mue_max, &
  !  & Posc_mue_antimax)
  ! PURPOSE
  ! Calculate Oscillation Probability
  !****************************************************************************
  subroutine oscillationProbability(Enu,L,deltaCP,Posc_mumu,Posc_mue,Posc_mue_max, &
     & Posc_mue_antimax)

    implicit none

    real, intent(in) :: enu ! neutrino energy  (GeV)
    real, intent(in)  :: L ! distance (in km) from the near to the far detector !
                           ! defined in initNeutrino
    real, optional, intent(in) :: deltaCP ! CP violating phase
    real, intent(out) :: Posc_mumu, Posc_mue, Posc_mue_max,Posc_mue_antimax

    !! parameters take from arXiv:1203.4090[hep-ex] Bishai, Diwan, ...
    !! "Neutrino Oscillations in the Precision Era"
    !! formulas are from the same reference, but adapted for V=0 (correspondingly hat{A}=0)
    !! sin(x) ~ x for small x has been used
    !! numerical factor 1.267 in oscillation expressions comes from conversion of
    !! units in deltaM2 * L/(4E)


    real, parameter :: deltaM2_32=2.5e-3,  sin_2theta23=1.0
    !! sin^2(2theta13)=0.1   theta13=9.2grad   sin_theta23=1./sqrt(2.)
    real, parameter :: sin_theta23=0.70710678, cos_theta23=0.70710678, sin_2theta13 = 0.316, &
       & cos_theta13 = 0.987
    real, parameter :: deltaM2_21 = 7.6e-5, sin_2theta12 = 0.927    !! theta12=34grad
    real :: deltaM2_31, sin_deltaCP, Posc_mutau, pid


    if (present(deltaCP)) then
       sin_deltaCP=sin(deltaCP)
    else
       sin_deltaCP=0.0
    end if

    Posc_mutau = sin_2theta23**2 * (sin(1.267*deltaM2_32*L/Enu))**2
    Posc_mumu  = 1.- Posc_mutau


    deltaM2_31 = deltaM2_32 + deltaM2_21

    if (present(deltaCP)) then
       sin_deltaCP=sin(deltaCP)
    else
       sin_deltaCP=1./sqrt(2.)
    end if

    !! In the following oscillation probabilities are calculated for
    !! deltaCP = pi/4,+pi/2 and -pi/2

    !! first switch for neutrino-antineutrino

    if (process_ID < -1) then
       pid = -1
    else if (process_ID > +1) then
       pid = +1
    else if (process_ID == +1 .or. process_ID == -1) then
       call TRACEBACK('oscillation for electrons makes no sense')
    end if

    Posc_mue = sin_theta23**2 * sin_2theta13**2 * (sin(1.267*deltaM2_31*L/Enu))**2 &
         &   + cos_theta23**2 * sin_2theta12**2 * (sin(1.267*deltaM2_21*L/Enu))**2 &
         &   + cos_theta13 * sin_2theta12 * sin_2theta13 * sin_2theta23  &
         &        * sin(1.267*deltaM2_21*L/Enu) * sin(1.267*deltaM2_31*L/Enu)&
         &        * (pid*sin_deltaCP * sin(1.267*deltaM2_31*L/Enu) + sqrt(1.-sin_deltaCP**2) &
         &   * cos(1.267*deltaM2_31*L/Enu) )


    sin_deltaCP=1.0
    Posc_mue_max = sin_theta23**2 * sin_2theta13**2 * (sin(1.267*deltaM2_31*L/Enu))**2 &
         &   + cos_theta23**2 * sin_2theta12**2 * (sin(1.267*deltaM2_21*L/Enu))**2 &
         &   + cos_theta13 * sin_2theta12 * sin_2theta13 * sin_2theta23  &
         &        * sin(1.267*deltaM2_21*L/Enu) * sin(1.267*deltaM2_31*L/Enu)&
         &        * (pid*sin_deltaCP * sin(1.267*deltaM2_31*L/Enu) + sqrt(1.-sin_deltaCP**2) &
         &   * cos(1.267*deltaM2_31*L/Enu) )


    sin_deltaCP=-1.0
    Posc_mue_antimax = sin_theta23**2 * sin_2theta13**2 * (sin(1.267*deltaM2_31*L/Enu))**2 &
         &   + cos_theta23**2 * sin_2theta12**2 * (sin(1.267*deltaM2_21*L/Enu))**2 &
         &   + cos_theta13 * sin_2theta12 * sin_2theta13 * sin_2theta23  &
         &        * sin(1.267*deltaM2_21*L/Enu) * sin(1.267*deltaM2_31*L/Enu)&
         &        * (pid*sin_deltaCP * sin(1.267*deltaM2_31*L/Enu) + sqrt(1.-sin_deltaCP**2) &
         &   * cos(1.267*deltaM2_31*L/Enu) )

  end subroutine oscillationProbability


  !****************************************************************************
  !****s* neutrinoAnalysis/PrintVals10
  ! NAME
  ! PURPOSE
  ! abbreviation used for many outputfiles: printing the vals stored in sigma,
  ! including its error
  !****************************************************************************
  subroutine PrintVals10(flag, raiseFlagVariable, numberOfCalls, fName, sigma)

    logical, intent(in) :: flag
    real, intent(in) :: raiseFlagVariable
    integer, intent(in) :: numberOfCalls
    character*(*), intent(in) :: fName
    real, dimension(1:dimSigma,1:2), intent(in) :: sigma

    if (.not.flag) return

    open(10,File=fName,position='append')
    if (numberOfCalls.ne.1) backspace(10)
    if (numberOfCalls.gt.1) then
       write(10,'(300E13.5)') raiseFlagVariable,sigma(2:,1)/float(numberOfCalls), &
            sqrt(max(0.,((sigma(2:,2)-sigma(2:,1)**2/float(numberOfCalls)) &
            / (float(numberOfCalls-1)*float(numberOfCalls)))))
    else
       write(10,'(300E13.5)') raiseFlagVariable,sigma(2:,1)
    end if
    close(10)

  end subroutine PrintVals10


  !****************************************************************************
  !****s* neutrinoAnalysis/PrintHeader10
  ! NAME
  ! PURPOSE
  ! abbreviation used for many outputfiles: printing a header line
  !****************************************************************************
  subroutine WriteHeader10(flag, fName)

    logical, intent(in) :: flag
    character*(*), intent(in) :: fName

    character*(*), parameter :: header = &
       "(' # total Xsections for defined finalstates (see sigma_0pions.dat)',/, &
       & ' # next 120 columns: error',/, &
       & ' # order of Xsections as in sigma_0pions.dat')"

    if (.not.flag) return

    open(10,file=fName)
    write(10,header)
    close(10)

  end subroutine WriteHeader10


end module neutrinoAnalysis
