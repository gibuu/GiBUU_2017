!******************************************************************************
!****m* /IdTable
! NAME
! Module IdTable
! PURPOSE
! Here the IDs of all particles are stored.
! SOURCE
!
Module IdTable
  implicit none
  public

  integer, parameter :: nres=30      !number of nucleon resonances
  integer, parameter :: nsres=21     !number of baryon resonances with s=-1
  integer, parameter,private :: ns2res=2     !number of baryon resonances with s=-2
  integer, parameter,private :: ns3res=1     !number of baryon resonances with s=-3
  integer, parameter,private :: ncres=6      !number of baryon resonances with charm c=+1
  integer, parameter :: nmes=22      !number of mesons
  integer, parameter :: nbar=1+nres+nsres+ns2res+ns3res+ncres !number of baryons
  integer, parameter,private :: idbmax=100   !maximal number of baryon resonances

  ! special ID codes
  integer, parameter :: EOV=-1     ! "end of vector", signals end of particle vector
  integer, parameter :: NOP= 0     ! "no particle"

  ! Nucleon Resonances
  integer, parameter :: nucleon=1
  integer, parameter :: delta=2
  integer, parameter :: P11_1440=3
  integer, parameter :: S11_1535=4
  integer, parameter :: S11_1650=5
  integer, parameter :: S11_2090=6
  integer, parameter :: D13_1520=7
  integer, parameter :: D13_1700=8
  integer, parameter :: D13_2080=9
  integer, parameter :: D15_1675=10
  integer, parameter :: G17_2190=11
  integer, parameter :: P11_1710=12
  integer, parameter :: P11_2100=13
  integer, parameter :: P13_1720=14
  integer, parameter :: P13_1900=15
  integer, parameter :: F15_1680=16
  integer, parameter :: F15_2000=17
  integer, parameter :: F17_1990=18
  integer, parameter :: S31_1620=19
  integer, parameter :: S31_1900=20
  integer, parameter :: D33_1700=21
  integer, parameter :: D33_1940=22
  integer, parameter :: D35_1930=23
  integer, parameter :: D35_2350=24
  integer, parameter :: P31_1750=25
  integer, parameter :: P31_1910=26
  integer, parameter :: P33_1600=27
  integer, parameter :: P33_1920=28
  integer, parameter :: F35_1750=29
  integer, parameter :: F35_1905=30
  integer, parameter :: F37_1950=31

  ! Baryons with strangeness s=-1
  integer, parameter :: Lambda=1+nres+1
  integer, parameter :: SigmaResonance=1+nres+2
  integer, parameter :: Sigma_1385=1+nres+3
  integer, parameter :: Lambda_1405=1+nres+4
  integer, parameter :: Lambda_1520=1+nres+5
  integer, parameter :: Lambda_1600=1+nres+6
  integer, parameter :: Lambda_1670=1+nres+7
  integer, parameter :: Lambda_1690=1+nres+8
  integer, parameter :: Lambda_1810=1+nres+9
  integer, parameter :: Lambda_1820=1+nres+10
  integer, parameter :: Lambda_1830=1+nres+11
  integer, parameter :: Sigma_1670=1+nres+12
  integer, parameter :: Sigma_1775=1+nres+13
  integer, parameter :: Sigma_2030=1+nres+14
  integer, parameter :: Lambda_1800=1+nres+15
  integer, parameter :: Lambda_1890=1+nres+16
  integer, parameter :: Lambda_2100=1+nres+17
  integer, parameter :: Lambda_2110=1+nres+18
  integer, parameter :: Sigma_1660=1+nres+19
  integer, parameter :: Sigma_1750=1+nres+20
  integer, parameter :: Sigma_1915=1+nres+21

  ! Baryons with strangeness s=-2
  integer, parameter :: Xi=1+nres+nsres+1
  integer, parameter :: XiStar=1+nres+nsres+2

  ! Baryons with strangeness s=-3
  integer, parameter :: OmegaResonance=1+nres+nsres+ns2res+1

  ! Baryons with charm c=+1
  integer, parameter :: Lambda_cPlus=1+nres+nsres+ns2res+ns3res+1
  integer, parameter :: Sigma_c=1+nres+nsres+ns2res+ns3res+2
  integer, parameter :: Sigma_cStar=1+nres+nsres+ns2res+ns3res+3
  integer, parameter :: Xi_c=1+nres+nsres+ns2res+ns3res+4
  integer, parameter :: Xi_cStar=1+nres+nsres+ns2res+ns3res+5
  integer, parameter :: Omega_c=1+nres+nsres+ns2res+ns3res+6


  ! Mesons
  integer, parameter :: pion=idbmax+1
  integer, parameter :: eta=idbmax+2
  integer, parameter :: rho=idbmax+3
  integer, parameter :: sigmaMeson=idbmax+4
  integer, parameter :: omegaMeson=idbmax+5
  integer, parameter :: etaPrime=idbmax+6
  integer, parameter :: phi=idbmax+7
  integer, parameter :: etaC=idbmax+8
  integer, parameter :: JPsi=idbmax+9
  integer, parameter :: kaon=idbmax+10     !(k+,k0)
  integer, parameter :: kaonBar=idbmax+11  !(k-,k0bar)
  integer, parameter :: kaonStar=idbmax+12
  integer, parameter :: kaonStarBar=idbmax+13
  integer, parameter :: dMeson=idbmax+14
  integer, parameter :: dBar=idbmax+15
  integer, parameter :: dStar=idbmax+16
  integer, parameter :: dStarBar=idbmax+17
  integer, parameter :: dS_plus=idbmax+18
  integer, parameter :: dS_minus=idbmax+19
  integer, parameter :: dSStar_plus=idbmax+20
  integer, parameter :: dSStar_minus=idbmax+21
  integer, parameter :: f2_1270=idbmax+22


  ! Gauge Bosons
  integer, parameter :: Wboson=997
  integer, parameter :: Zboson=998
  integer, parameter :: photon=999


  ! Leptons
  integer, parameter :: electron=901
  integer, parameter :: electronNeutrino=911 ! (-911 for antineutrino)

  integer, parameter :: muon=902
  integer, parameter :: muonNeutrino=912     ! (-912 for antineutrino)

  integer, parameter :: tau=903
  integer, parameter :: tauNeutrino=913      ! (-913 for antineutrino)

!******************************************************************************

contains


  !****************************************************************************
  !****f* IdTable/invalidID
  ! NAME
  ! logical function invalidID(id)
  ! INPUTS
  ! * integer :: ID
  ! RESULT
  ! * .true. if ID is invalid
  ! NOTES
  ! This functions is supposed to detect ID values which should *never* occur.
  !****************************************************************************
  pure logical function invalidID (id)
    integer, intent(in) :: id
    invalidID = (id<EOV .or. id>photon)
  end function

  !****************************************************************************
  !****f* IdTable/isPhoton
  ! NAME
  ! logical function isPhoton(ID)
  ! INPUTS
  ! * integer :: ID
  ! RESULT
  ! * .true. if ID corresponds to a photon
  !****************************************************************************
  pure logical function isPhoton(ID)
    integer, intent(in) :: ID
    isPhoton = (ID==photon)
  end function isPhoton


  !****************************************************************************
  !****f* IdTable/isLepton
  ! NAME
  ! logical function isLepton(ID)
  ! INPUTS
  ! * integer :: ID
  ! RESULT
  ! * .true. if ID corresponds to a lepton
  !****************************************************************************
  pure logical function isLepton (ID)
    integer, intent(in) :: ID
    isLepton = (ID>=electron .and. ID<=tauNeutrino)
  end function isLepton


  !****************************************************************************
  !****f* IdTable/isHadron
  ! NAME
  ! logical function isHadron(ID)
  ! INPUTS
  ! * integer :: ID
  ! RESULT
  ! * .true. if ID corresponds to a hadron
  !****************************************************************************
  pure logical function isHadron (ID)
    integer, intent(in) :: ID
    isHadron = (ID>=nucleon .and. ID<pion+nMes)
  end function isHadron


  !****************************************************************************
  !****f* IdTable/isBaryon
  ! NAME
  ! logical function isBaryon(ID)
  ! INPUTS
  ! * integer :: ID
  ! RESULT
  ! * .true. if ID corresponds to a baryon
  !****************************************************************************
  pure logical function isBaryon (ID)
    integer, intent(in) :: ID
    isBaryon = (ID>=nucleon .and. ID<=nbar)
  end function isBaryon


  !****************************************************************************
  !****s* IdTable/isBaryonResonance
  ! NAME
  ! logical function isBaryonResonance(ID)
  !
  ! PURPOSE
  ! Returns .true. if the ID corresponds to a baryonic resonance,
  ! i.e. any baryon except the nucleon.
  !
  ! INPUTS
  ! integer, intent(in) :: ID
  !****************************************************************************
  pure logical function isBaryonResonance (ID)
    integer,intent(in) :: ID
    isBaryonResonance = (ID>nucleon .and. ID<=nbar)
  end function isBaryonResonance


  !****************************************************************************
  !****s* IdTable/isHyperon
  ! NAME
  ! logical function isHyperon(ID)
  !
  ! PURPOSE
  ! Returns .true. if the ID corresponds to a hyperon,
  ! i.e. a baryon containing strange quarks but no charm.
  !
  ! INPUTS
  ! integer, intent(in) :: ID
  !****************************************************************************
  pure logical function isHyperon (ID)
    integer,intent(in) :: ID
    isHyperon = (ID>=Lambda .and. ID<=OmegaResonance)
  end function isHyperon


  !****************************************************************************
  !****f* IdTable/isMeson
  ! NAME
  ! logical function isMeson(ID)
  ! INPUTS
  ! * integer :: ID
  ! RESULT
  ! * .true. if ID corresponds to a meson
  !****************************************************************************
  pure logical function isMeson (ID)
    integer, intent(in) :: ID
    isMeson = (ID>=pion .and. ID<pion+nMes)
  end function isMeson


  !****************************************************************************
  !****s* IdTable/getAntiMeson
  ! NAME
  ! subroutine getAntiMeson(id, charge, antiID, antiCharge)
  ! PURPOSE
  ! This routine returnes for a given meson the ID and the charge
  ! of its antiparticle.
  ! INPUTS
  ! * integer :: id, charge         ! of particle
  ! OUTPUT
  ! * integer :: antiID, antiCharge ! Id and charge of particle
  !****************************************************************************
  subroutine getAntiMeson(id, charge, antiID, antiCharge)
    integer, intent(in) :: id, charge
    integer, intent(out) :: antiId, antiCharge

    if (.not.isMeson(id)) then
       write(*,*) 'Problem in IdTable/getAntiMeson'
       write(*,*) 'ID is no Meson :', id
       stop
    end if

    antiCharge=-charge

    select case (id)
    case (pion:JPsi,f2_1270)
       antiID=id
    case (kaon)
       antiID=kaonBar
    case (kaonBar)
       antiID=kaon
    case (kaonStar)
       antiID=kaonStarBar
    case (kaonStarBar)
       antiID=kaonStar
    case (dMeson)
       antiID=dBar
    case (dBar)
       antiID=dMeson
    case (dStar)
       antiID=dStarBar
    case (dStarBar)
       antiID=dStar
    case (ds_plus)
       antiID=ds_minus
    case (ds_minus)
       antiID=ds_plus
    case (dsStar_plus)
       antiID=dsStar_minus
    case (dsStar_minus)
       antiID=dsStar_plus
    case default
       write(*,*) 'Problem in particleProperties/getAntiMeson'
       write(*,*) 'Antiparticle for this particle is not yet implemented :', id
    end select

  end subroutine getAntiMeson


End Module IdTable
