!******************************************************************************
!****m* /eventtypes
! NAME
! module eventtypes
! PURPOSE
! Declare constants which define the different event classes:
! * 000 -- Elementary collisions (hadron+hadron, real particles)
! * 001 -- Heavy-Ion collisions (nucleus+nucleus, real particles)
! * 002 -- Pion-induced reactions (low energy, perturbative particles)
! * 003 -- Real-photon-induced reactions
! * 004 -- Lepton (i.e. virtual photon) induced reactions (low energy)
! * 005 -- Neutrino-induced reactions
! * 012 -- Pion- and nucleon-induced reactios (high energy, perturbative particles)
! * 014 -- Lepton (i.e. virtual photon) induced reactions (high energy)
! * 022 -- Transport of an external hadronic source from a data file
! * 031 -- Nucleons in a [box of nucleons] (continuous boundary conditions)
! * 032 -- Pions in a [box of nucleons] (continuous boundary conditions)
! * 033 -- Deltas in a [box of nucleons] (continuous boundary conditions)
! * 041 -- Box of particles
! * 100 -- Groundstate calculation
! * 200 -- Simple transport of a given particle
! * 300 -- Hadron-induced reactions (hadron+nucleus, real particles)
!******************************************************************************
module eventtypes

  implicit none

  integer, parameter :: elementary             = 0
  integer, parameter :: HeavyIon               = 1

  integer, parameter :: LoPion                 = 2
  integer, parameter :: RealPhoton             = 3
  integer, parameter :: LoLepton               = 4

  integer, parameter :: Neutrino               = 5

  integer, parameter :: HiPion                 = 12
  integer, parameter :: HiLepton               = 14

  integer, parameter :: ExternalSource         = 22

  integer, parameter :: InABox                 = 31
  integer, parameter :: InABox_pion            = 32
  integer, parameter :: InABox_delta           = 33

  integer, parameter :: Box                    = 41

  integer, parameter :: groundState            = 100
  integer, parameter :: transportGivenParticle = 200
  integer, parameter :: hadron                 = 300

end module eventtypes
