
module neutrino_IDTable
  implicit none

  public

  !values for flavor_ID
  integer, parameter :: electron  = 1
  integer, parameter :: muon      = 2
  integer, parameter :: taulepton = 3

  !values for XsectionMode
  integer, parameter :: integratedsigma=0
  integer, parameter :: dsigmadcosthetadelepton=1
  integer, parameter :: dsigmadQsdelepton=2
  integer, parameter :: dsigmadQs=3
  integer, parameter :: dsigmadcostheta=4
  integer, parameter :: dsigmadelepton=5
  integer, parameter :: dsigmaMC=6
  integer, parameter :: dsigmadW=7

  integer, parameter :: EXP_dSigmadEnu=10
  integer, parameter :: EXP_dsigmadcosthetadelepton=11
  integer, parameter :: EXP_dsigmadQsdelepton=12
  integer, parameter :: EXP_dsigmadQs=13
  integer, parameter :: EXP_dsigmadcostheta=14
  integer, parameter :: EXP_dsigmadelepton=15
  integer, parameter :: EXP_dsigmaMC=16
  integer, parameter :: EXP_dsigmadW=17

  character*(*), dimension(0:7), parameter :: sXsectionMode = (/&
       & "sigma            ", &
       & "dsigma/ dcost dE'", &
       & "dsigma/ dQ2 dE'  ", &
       & "dsigma/ dQ2      ", &
       & "dsigma/ dcost    ", &
       & "dsigma/ dE'      ", &
       & "sigma (MC)       ", &
       & "dsigma/ dW       "/)

  !values for nuExp
  integer, parameter :: MiniBooNE = 1
  integer, parameter :: K2K = 2
  integer, parameter :: Minos =3

  !values for onePion
  !n for outgoing channels with neutron
  !p for outgoing channels with proton
  integer,parameter :: onePionCH_n=32
  integer,parameter :: onePionCH_p=33

  integer,parameter :: DIS_CH=34

  integer,parameter :: QE2p2h = 35
  integer,parameter :: Delta2p2h = 36

  integer,parameter :: twoPion = 37




end module neutrino_IDTable
