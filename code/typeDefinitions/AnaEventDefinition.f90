!******************************************************************************
!****m* /AnaEventDefinition
! NAME
! module AnaEventDefinition
! PURPOSE
! Here type(tAnaEvent) is defined.
! Routines to work with this type are defined elsewhere, as e.g. in
! module AnaEvent
!******************************************************************************
module AnaEventDefinition

  use IDTABLE, only: pion, eta,kaon,kaonBar,DMeson, DBar, Ds_plus,Ds_minus, &
       & nucleon, lambda, sigmaResonance,Xi, OmegaResonance
! Hadrons listed here can be included in detailed analyses, see  namelist &detailed_diff

  use particlePointerListDefinition, only: tParticleList

  implicit none




  !****************************************************************************
  !****ig* AnaEvent/particleIDs
  ! PURPOSE
  ! The ID's of the particles we consider "stable". For these detailed analyses
  ! can be performed
  ! as determined by namelist &detailed_diff
  ! public is necessary because it is used in LArAnalysis
  !
  ! numStableParts = dimension of array particleIDs = number of stable particles
  ! for which final state analyses will be done.
  ! numStableMesons is number of long-lived mesons (width < 1E-4 GeV)
  !
  ! NOTE: all the stable particles also have to be listed in subroutine event_add
  ! which is contained in module AnaEvent
  !
  ! For higher energies more, heavier hadrons may appear. For these the
  ! array has to be extended together with numStableParts and numStableMesons
  ! SOURCE
  !
  integer, parameter, public :: numStableParts = 13
  integer, parameter, public :: numStableMesons = 8

  integer, dimension(1:numStableParts),parameter, public :: particleIDs=(/ &
       & pion, eta, kaon, kaonBar, DMeson, dBar,ds_plus,ds_minus,&
       & nucleon, lambda, sigmaResonance, Xi, OmegaResonance/)
  !****************************************************************************


  !****************************************************************************
  !****t* AnaEventDefinition/tAnaEvent
  ! NAME
  ! Type tAnaEvent
  ! PURPOSE
  ! Type definition for events.
  ! The first index in numberParticles runs over all stable hadrons
  ! used in the analysis + 1, (numStableParts + 1), 
  ! the second one over the possible charge states (-2 -> +2)
  !
  ! SOURCE
  !
  Type tAnaEvent
     sequence
     type(tParticleList) :: particleList                 ! particles in the event
     integer,dimension(1:14,-2:2) :: numberParticles =0  ! Counters for stable particles
                                                         ! (under strong decays)
  End Type tAnaEvent
  ! NOTES
  ! The field numberParticles includes the multiplicities of particles.
  ! Antiparticles are not counted!
  !
  ! 1st Index (see array particleIDs):
  ! * 1=pion
  ! * 2=eta
  ! * 3=kaon
  ! * 4=kaonBar
  ! * 5=dMeson
  ! * 6=dBar
  ! * 7=ds_plus
  ! * 8=ds_minus
  ! * 9=Nucleon
  ! * 10=Lambda
  ! * 11=Sigma
  ! * 12=Xi
  ! * 13=OmegaResonance
  ! * 14=Any other
  !
  ! 2nd Index:
  ! * Particle Charge (-2:2)
  !****************************************************************************



end module AnaEventDefinition
