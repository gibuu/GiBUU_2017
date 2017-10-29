!******************************************************************************
!****m* /particleProperties
! NAME
! module particleProperties
! PURPOSE
! Defines a database structure that contains the particle properties of all
! hadrons (baryons and mesons).
!******************************************************************************
Module ParticleProperties
  use IdTable, only: nMes, pion
  implicit none
  private

  integer, parameter, public :: nDecays = 9


  !****************************************************************************
  !****t* particleProperties/hadronProperties
  ! NAME
  ! type hadronProperties
  ! SOURCE
  !
  type, public :: hadronProperties
     integer, dimension(nDecays) :: decaysID         = 0         ! the channel id
     real, dimension(nDecays)    :: decays           = 0.        ! the branching ratio
     character(len=12)           :: name             = 'Invalid'
     character(len=14)           :: nameTeX          = 'Invalid'
     real                        :: mass             = 0.
     real                        :: width            = 0.
     real                        :: spin             = 0. ! J
     integer                     :: AngularMomentum  = 0  ! Angular Momentum of Nucleon-Pion or Lambda-Kaon final state.
     integer                     :: isoSpinTimes2    = 0  ! I
     integer                     :: strangeness      = 0
     integer                     :: charm            = 0
     integer                     :: rating           = 0
     logical                     :: propagated       = .false.
     logical                     :: usedForXsections = .false.
     integer                     :: stability        = 0
     real                        :: minmass          = 0.
  end type hadronProperties
  !
  ! NOTES
  ! Width of different decay channels for 2 Body Decays in the vacuum at
  ! the pole mass:
  ! The width are supposed to be the same within an isospin multiplet.
  !
  ! The additional flag 'stability' encodes on a bitwise level, how the
  ! particle may decay (cf. also master_1Body):
  ! * 1: particle may decay during run, if gammaMed > gammaCutOff
  ! * 2: particle may decay at the end of the run, if gammaMed > 1e-10
  ! * 4: particle may decay at the end via Jetset, if there the parameters
  !   allow for a decay.
  ! The default values are 0 or 3. You can change this in your jobcard via
  ! the flag StabilityFlag in the namelist "InitStability".
  !
  ! The decay channel IDs are positive for 2-body decays and negative for
  ! 3-body decays. Whether it is a meson decay or a baryon decay has to
  ! be decided by the property of the parent particle.
  !****************************************************************************



  !****************************************************************************
  !****g* particleProperties/hadron
  ! SOURCE
  type(hadronProperties), dimension(1:pion+nMes-1), save, public :: hadron
  ! PURPOSE
  ! Database of hadron properties.
  ! USAGE
  ! * To get e.g. the mass of the nucleon type "hadron(nucleon)%mass".
  !   The ID's are chosen according to IdTable.
  ! NOTES
  ! * Since the baryons cover the range 1:61 and the mesons 101:122, there
  !   is a 'hole' for 62:100.
  !****************************************************************************


  public :: InitParticleProperties
  public :: validCharge, validCharge_ID
  public :: PartName, TexName
  public :: isStrange, isCharmed, isNonExotic
  public :: get_rho_dilep
  public :: getAngularMomentum_meson, getAngularMomentum_baryon


  !****************************************************************************
  !****g* particleProperties/propagationSwitch
  ! SOURCE
  integer, save :: propagationSwitch = 3
  ! PURPOSE
  ! * 0 = propagate resonances with more than 1 star in their rating
  !   (irrespar=0 in old code)
  ! * 1 = propagate just the Delta (irrespar=2 in old code)
  ! * 2 = propagate no resonance (irrespar=3 in old code)
  ! * 3 = propagate all resonances (default)
  !****************************************************************************


  !****************************************************************************
  !****g* particleProperties/usageForXsectionSwitch
  ! SOURCE
  !
  integer, save :: usageForXsectionSwitch = 2
  ! PURPOSE
  ! * 0 = use resonances with more than 1 star in their rating for cross sections
  ! * 1 = use all resonances for cross sections
  ! * 2 = use all resonances besides the 1* star I=1/2 resonances (default)
  ! * 3 = use only the Delta
  !****************************************************************************


  !****************************************************************************
  !****g* particleProperties/rho_dilep
  ! SOURCE
  !
  logical, save :: rho_dilep = .false.
  ! PURPOSE
  ! If .false. (default), the rho meson width will be exclusively given by the
  ! 2pi decay and its minmass will be 2m_pi.
  ! If .true., the dilepton width will be included in the width and spectral
  ! function of the rho, and the minmass will be 2m_e.
  ! This is important for dilepton spectra, in order to get contributions from
  ! the rho below the 2pi threshold.
  !****************************************************************************


  !****************************************************************************
  !****g* particleProperties/FileNameDecayChannels
  ! PURPOSE
  ! The absolute filename of the file containing decay channel infos.
  !
  ! possible values:
  ! * if not set, default is '[path_To_Input]/DecayChannels.dat'
  ! * if given, but does not contain '/':
  !   default is '[path_To_Input]/[FileNameDecayChannels]'
  ! * otherwise: filename is absolute, including path
  !
  ! NOTE
  ! if you want to use the file 'XXX.dat' in the actual directory,
  ! give it as './XXX.dat'
  !
  ! SOURCE
  !
  character(300), save :: FileNameDecayChannels = ''
  !****************************************************************************


  !****************************************************************************
  !****f* particleProperties/PartName
  ! NAME
  ! function PartName(ID,IQ,isAnti)
  ! function PartName(ID)
  ! function PartName(Part)
  ! PURPOSE
  ! return a String containing name, charge etc. of the particle
  ! INPUTS
  ! * integer, intent(in) :: ID,IQ  -- ID and charge of particle
  ! * logical, intent(in) :: isAnti -- true if antiparticle, false if particle
  ! or :
  ! * type(particle) :: part -- particle
  ! USAGE
  ! * character*() = PartName(ID,IQ,isAnti)
  ! * character*() = PartName(Part)
  ! NOTES
  ! if only ID is given, the retruned name only contains the particle name,
  ! i.e. it is in a short version
  !****************************************************************************
  Interface PartName
     Module Procedure PartName1, PartName2, PartName3
  end Interface


  !****************************************************************************
  !****f* particleProperties/TeXName
  ! NAME
  ! function TeXName(ID,IQ,isAnti)
  ! function TeXName(ID)
  ! function TeXName(Part)
  ! PURPOSE
  ! return a String containing name, charge etc. of the particle.
  ! Fomramtting according TeX rules
  ! INPUTS
  ! * integer, intent(in) :: ID,IQ  -- ID and charge of particle
  ! * logical, intent(in) :: isAnti -- true if antiparticle, false if particle
  ! or :
  ! * type(particle) :: part -- particle
  ! NOTES
  ! if only ID is given, the retruned name only contains the particle name,
  ! i.e. it is in a short version
  !****************************************************************************
  Interface TexName
     Module Procedure TeXName1, TeXName2, TeXName3
  end Interface


contains


  logical function get_rho_dilep ()
    get_rho_dilep = rho_dilep
  end function


  !****************************************************************************
  !****s* particleProperties/InitParticleProperties
  ! NAME
  ! subroutine InitParticleProperties
  ! PURPOSE
  ! Define the variables baryon and meson with the Parameter Sets out of
  ! Manley and PDG
  !****************************************************************************
  subroutine InitParticleProperties
    use idtable
    use decayChannels, only: initDecayChannels, Print_DecayChannels

    call Init ! read jobcard

    ! Set Masses, Name, Index, ...
    call setParameters
    call decideOnPropagationUsage
    ! Set width of the decay channels
    call initDecayChannels
    call setThresholds
    call setDecays
    call updateMinmassesAndThresholds

    call ReadInitModify ! read additional jobcard

    ! Print out the properties and the decay channels of all particles
    call PrintParticleProperties

    call Print_DecayChannels

  contains


    !**************************************************************************
    !****s* InitParticleProperties/setParameters
    ! NAME
    ! subroutine setParameters
    ! PURPOSE
    ! Define the variables baryon and meson.
    ! Here all the parameters besides the decay channels
    ! are set.
    !**************************************************************************
    subroutine setParameters
      use constants, only: mN,mPi,mK,mElec

      ! --- Non-strange, non-charmed baryons

      hadron(nucleon) = initBaryon('N','N',mN,0.,0.5,4,1,0,0,1)
      hadron(nucleon)%minmass = 0.7

      hadron(Delta)    = initBaryon('Delta'     ,'\Delta'      ,1.232,0.118,1.5,4,3,0,0,1)
      hadron(P11_1440) = initBaryon('P_11(1440)','P_{11}(1440)',1.462,0.391,0.5,4,1,0,0,1)
      hadron(S11_1535) = initBaryon('S_11(1535)','S_{11}(1535)',1.534,0.151,0.5,3,1,0,0,0)
      hadron(S11_1650) = initBaryon('S_11(1650)','S_{11}(1650)',1.659,0.173,0.5,4,1,0,0,0)
      hadron(S11_2090) = initBaryon('S_11(2090)','S_{11}(2090)',1.928,0.414,0.5,1,1,0,0,0)
      hadron(D13_1520) = initBaryon('D_13(1520)','D_{13}(1520)',1.524,0.124,1.5,4,1,0,0,2)
      hadron(D13_1700) = initBaryon('D_13(1700)','D_{13}(1700)',1.737,0.249,1.5,1,1,0,0,2)
      hadron(D13_2080) = initBaryon('D_13(2080)','D_{13}(2080)',1.804,0.447,1.5,1,1,0,0,2)
      hadron(D15_1675) = initBaryon('D_15(1675)','D_{15}(1675)',1.676,0.159,2.5,4,1,0,0,2)
      hadron(G17_2190) = initBaryon('G_17(2190)','G_{17}(2190)',2.127,0.547,3.5,4,1,0,0,4)
      hadron(P11_1710) = initBaryon('P_11(1710)','P_{11}(1710)',1.717,0.478,0.5,1,1,0,0,1)
      hadron(P11_2100) = initBaryon('P_11(2100)','P_{11}(2100)',1.885,0.113,0.5,1,1,0,0,1)
      hadron(P13_1720) = initBaryon('P_13(1720)','P_{13}(1720)',1.717,0.383,1.5,1,1,0,0,1)
      hadron(P13_1900) = initBaryon('P_13(1900)','P_{13}(1900)',1.879,0.498,1.5,3,1,0,0,1)
      hadron(F15_1680) = initBaryon('F_15(1680)','F_{15}(1680)',1.684,0.139,2.5,4,1,0,0,3)
      hadron(F15_2000) = initBaryon('F_15(2000)','F_{15}(2000)',1.903,0.494,2.5,1,1,0,0,3)
      hadron(F17_1990) = initBaryon('F_17(1990)','F_{17}(1990)',2.086,0.535,3.5,2,1,0,0,3)
      hadron(S31_1620) = initBaryon('S_31(1620)','S_{31}(1620)',1.672,0.154,0.5,2,3,0,0,0)
      hadron(S31_1900) = initBaryon('S_31(1900)','S_{31}(1900)',1.920,0.263,0.5,3,3,0,0,0)
      hadron(D33_1700) = initBaryon('D_33(1700)','D_{33}(1700)',1.762,0.599,1.5,1,3,0,0,2)
      hadron(D33_1940) = initBaryon('D_33(1940)','D_{33}(1940)',2.057,0.460,1.5,1,3,0,0,2)
      hadron(D35_1930) = initBaryon('D_35(1930)','D_{35}(1930)',1.956,0.526,2.5,2,3,0,0,2)
      hadron(D35_2350) = initBaryon('D_35(2350)','D_{35}(2350)',2.171,0.264,2.5,2,3,0,0,2)
      hadron(P31_1750) = initBaryon('P_31(1750)','P_{31}(1750)',1.744,0.299,0.5,1,3,0,0,1)
      hadron(P31_1910) = initBaryon('P_31(1910)','P_{31}(1910)',1.882,0.239,0.5,4,3,0,0,1)
      hadron(P33_1600) = initBaryon('P_33(1600)','P_{33}(1600)',1.706,0.430,1.5,3,3,0,0,1)
      hadron(P33_1920) = initBaryon('P_33(1920)','P_{33}(1920)',2.014,0.152,1.5,1,3,0,0,1)
      hadron(F35_1750) = initBaryon('F_35(1750)','F_{35}(1750)',1.752,0.251,2.5,1,3,0,0,3)
      hadron(F35_1905) = initBaryon('F_35(1905)','F_{35}(1905)',1.881,0.327,2.5,3,3,0,0,3)
      hadron(F37_1950) = initBaryon('F_37(1950)','F_{37}(1950)',1.945,0.300,3.5,4,3,0,0,3)

      hadron(Delta:F37_1950)%minmass = mN + mPi

      ! --- Strange Baryons s=-1

      hadron(lambda)         = initBaryon('Lambda'      ,'\Lambda'      ,1.116,0.   ,0.5,4,0,-1,0,0)
      hadron(SigmaResonance) = initBaryon('Sigma'       ,'\Sigma'       ,1.189,0.   ,0.5,4,2,-1,0,0)
      hadron(Sigma_1385)     = initBaryon('Sigma(1385)' ,'\Sigma(1385)' ,1.385,0.036,1.5,4,2,-1,0,1)
      hadron(lambda_1405)    = initBaryon('Lambda(1405)','\Lambda(1405)',1.405,0.05 ,0.5,4,0,-1,0,0)
      hadron(lambda_1520)    = initBaryon('Lambda(1520)','\Lambda(1520)',1.520,0.016,1.5,4,0,-1,0,2)
      hadron(lambda_1600)    = initBaryon('Lambda(1600)','\Lambda(1600)',1.600,0.150,0.5,3,0,-1,0,1)
      hadron(lambda_1670)    = initBaryon('Lambda(1670)','\Lambda(1670)',1.670,0.035,0.5,4,0,-1,0,0)
      hadron(lambda_1690)    = initBaryon('Lambda(1690)','\Lambda(1690)',1.690,0.06 ,1.5,4,0,-1,0,2)
      hadron(lambda_1810)    = initBaryon('Lambda(1810)','\Lambda(1810)',1.810,0.15 ,0.5,3,0,-1,0,1)
      hadron(lambda_1820)    = initBaryon('Lambda(1820)','\Lambda(1820)',1.820,0.08 ,2.5,4,0,-1,0,3)
      hadron(lambda_1830)    = initBaryon('Lambda(1830)','\Lambda(1830)',1.830,0.095,2.5,4,0,-1,0,2)
      hadron(sigma_1670)     = initBaryon('Sigma(1670)' ,'\Sigma(1670)' ,1.670,0.06 ,1.5,4,2,-1,0,2)
      hadron(sigma_1775)     = initBaryon('Sigma(1775)' ,'\Sigma(1775)' ,1.775,0.12 ,2.5,4,2,-1,0,2)
      hadron(sigma_2030)     = initBaryon('Sigma(2030)' ,'\Sigma(2030)' ,2.030,0.18 ,3.5,4,2,-1,0,3)
      hadron(lambda_1800)    = initBaryon('Lambda(1800)','\Lambda(1800)',1.800,0.3  ,0.5,3,0,-1,0,0)
      hadron(lambda_1890)    = initBaryon('Lambda(1890)','\Lambda(1890)',1.890,0.1  ,1.5,4,0,-1,0,1)
      hadron(lambda_2100)    = initBaryon('Lambda(2100)','\Lambda(2100)',2.1  ,0.2  ,3.5,4,0,-1,0,4)
      hadron(lambda_2110)    = initBaryon('Lambda(2110)','\Lambda(2110)',2.110,0.2  ,2.5,3,0,-1,0,3)
      hadron(sigma_1660)     = initBaryon('Sigma(1660)' ,'\Sigma(1660)' ,1.660,0.1  ,0.5,3,2,-1,0,1)
      hadron(sigma_1750)     = initBaryon('Sigma(1750)' ,'\Sigma(1750)' ,1.750,0.09 ,0.5,3,2,-1,0,0)
      hadron(sigma_1915)     = initBaryon('Sigma(1915)' ,'\Sigma(1915)' ,1.915,0.12 ,2.5,4,2,-1,0,3)

      !The Angular Momenta the strange resonances in the partial wave analysis of the kaonBar-nucleon channel
      hadron(lambda:sigmaResonance)%minMass = mN + mPi

      ! --- Baryons with strangeness s=-2
      hadron(Xi)             = initBaryon('Xi'   ,'\Xi'   ,1.315,0.0   ,0.5,4,1,-2,0,0)
      hadron(XiStar)         = initBaryon('Xi*'  ,'\Xi^*' ,1.530,0.0095,1.5,4,1,-2,0,0)
      ! --- Baryons with strangeness s=-3
      hadron(OmegaResonance) = initBaryon('Omega','\Omega',1.672,0.0   ,1.5,4,0,-3,0,0)

      hadron(sigma_1385:OmegaResonance)%minMass = hadron(lambda)%mass + mPi

      ! --- Baryons with charm c=+1
      hadron(lambda_cPlus)   = initBaryon('Lambda_c','\Lambda_c' ,2.285,0.0  ,0.5,4,0, 0,1,0)
      hadron(sigma_c)        = initBaryon('Sigma_c' ,'\Sigma_c'  ,2.452,0.0  ,0.5,4,2, 0,1,0)
      hadron(Sigma_cStar)    = initBaryon('Sigma_c*','\Sigma_c^*',2.520,0.015,1.5,4,2, 0,1,0)
      hadron(Xi_c)           = initBaryon('Xi_c'    ,'\Xi_c'     ,2.466,0.0  ,0.5,3,1,-1,1,0)
      hadron(Xi_cStar)       = initBaryon('Xi_c*'   ,'\Xi_c^*'   ,2.645,0.004,1.5,3,1,-1,1,0)
      hadron(omega_c)        = initBaryon('Omega_c' ,'\Omega_c'  ,2.6975,0.0 ,0.5,3,0,-2,1,0)

      hadron(lambda_cPlus:sigma_c)%minMass = mN + mPi
      hadron(Sigma_cStar:omega_c)%minMass  = hadron(lambda_cPlus)%mass + mPi


      ! --- Mesons

      hadron(pion)         = initMeson('pi'   ,'\pi'           ,mPi    ,0.      ,0.,2, 0, 0, 0.)
      hadron(eta)          = initMeson('eta'  ,'\eta'          ,0.54785,1.30E-06,0.,0, 0, 0, 0.)

      if (rho_dilep) then
        hadron(rho)        = initMeson('rho'  ,'\rho'          ,0.7755 ,0.1491  ,1.,2, 0, 0, 2*mElec)
      else
        hadron(rho)        = initMeson('rho'  ,'\rho'          ,0.7755 ,0.1491  ,1.,2, 0, 0, 2*mPi)
      end if

      hadron(sigmaMeson)   = initMeson('sigma','\sigma'        ,0.800  ,0.500   ,0.,0, 0, 0, 2*mPi)
      hadron(omegaMeson)   = initMeson('omega','\omega'        ,0.7826 ,8.49E-3 ,1.,0, 0, 0, mPi)
      hadron(etaprime)     = initMeson('etaP' ,"\eta'"         ,0.95778,0.194E-3,0.,0, 0, 0, 0.)
      hadron(phi)          = initMeson('phi'  ,'\phi'          ,1.01945,4.26E-3 ,1.,0, 0, 0, 3*mPi)
      hadron(etaC)         = initMeson('eta_c','\eta_c'        ,2.980  ,28.E-3  ,0.,0, 0, 0, 0.)
      hadron(JPsi)         = initMeson('Jpsi' ,'J/\psi'        ,3.096916,0.      ,1.,0, 0, 0, 0.)

      hadron(kaon)         = initMeson('K'    ,'K'             ,mK     ,0.      ,0.,1, 1, 0, mK)
      hadron(kaonbar)      = initMeson('K~'   ,'\overline{K}'  ,mK     ,0.      ,0.,1,-1, 0, mK)
      hadron(kaonstar)     = initMeson('K*'   ,'K^*'           ,0.892  ,0.05    ,1.,1, 1, 0, mK+mPi)
      hadron(kaonstarbar)  = initMeson('K*~'  ,'\overline{K}^*',0.892  ,0.05    ,1.,1,-1, 0, mK+mPi)

      hadron(dMeson)       = initMeson('D'    ,'D'             ,1.867  ,0.      ,0.,1, 0, 1, 1.5)
      hadron(dbar)         = initMeson('D~'   ,'\overline{D}'  ,1.867  ,0.      ,0.,1, 0,-1, 1.5)
      hadron(dstar)        = initMeson('D*'   ,'D^*'           ,2.007  ,0.002   ,1.,1, 0, 1, 1.5)
      hadron(dStarbar)     = initMeson('D*~'  ,'\overline{D}^*',2.007  ,0.002   ,1.,1, 0,-1, 1.5)

      hadron(dS_Plus)      = initMeson('D_s+' ,'D_s^+'         ,1.969  ,0.      ,0.,0, 1, 1, 1.5)
      hadron(dS_Minus)     = initMeson('D_s-' ,'D_s^-'         ,1.969  ,0.      ,0.,0,-1,-1, 1.5)
      hadron(dsStar_Plus)  = initMeson('D_s*+','D_s^{*+}'      ,2.112  ,0.001   ,1.,0, 1, 1, 1.5)
      hadron(dsStar_Minus) = initMeson('D_s*-','D_s^{*-}'      ,2.112  ,0.001   ,1.,0,-1,-1, 1.5)
      hadron(f2_1270)      = initMeson('f_2(1270)','f_2(1270)' ,1.2754 ,0.1852  ,2.,0, 0, 0, 2*mPi)

    end subroutine setParameters


    !**************************************************************************
    !****s* InitParticleProperties/decideOnPropagationUsage
    ! NAME
    ! subroutine decideOnPropagationUsage
    ! PURPOSE
    ! Decides which baryons shall be explicitly propagated
    ! and used for the cross sections.
    !**************************************************************************
    subroutine decideOnPropagationUsage

      ! Decide who is used for crossSections:
      select case (usageForXsectionSwitch)
      case (0)
         ! All resonances with more than one star are used for crossSections
         hadron(1:nbar)%usedForXsections = (hadron(1:nbar)%rating>1)
      case (1)
         ! All resonances are used
         hadron(1:nbar)%usedForXsections = .true.
      case (2)
         ! All resonances besides the one-star I=1/2 resonances (default!)
         hadron(1:nbar)%usedForXsections = (hadron(1:nbar)%rating>1 .or. hadron(1:nbar)%isoSpinTimes2/=1)
      case (3)
         ! only Nucleon and Delta are used
         hadron(1:2)%usedForXsections    = .true.
         hadron(3:nbar)%usedForXsections = .false.
      case default
         write(*,*) "Bad value of usageForXsectionSwitch: ", usageForXsectionSwitch
         stop
      end select

      ! Decide who is propagated:
      select case (propagationSwitch)
      case (0)
         ! Everything with more than one star is propagated
         hadron(1:nbar)%propagated = (hadron(1:nbar)%rating>1)
      case (1)
         ! Only delta and nucleon are propagated
         hadron(1:2)%propagated    = .true.
         hadron(3:nbar)%propagated = .false.
      case (2)
         ! no resonance is propagated (only nucleon)
         hadron(1)%propagated      = .true.
         hadron(2:nbar)%propagated = .false.
      case (3)
         ! everybody is propagated (default!)
         hadron(1:nbar)%propagated = .true.
      case default
         write(*,*) "Bad value of propagationSwitch: ", propagationSwitch
         stop
      end select

    end subroutine decideOnPropagationUsage


    !**************************************************************************
    !****s* InitParticleProperties/setThresholds
    ! NAME
    ! subroutine setThresholds
    ! PURPOSE
    ! Set the thresholds of the decay channels, i.e. the lightest mass at which a resonance
    ! can decay into this channel (= sum of the daughter masses).
    !**************************************************************************
    subroutine setThresholds
      use decayChannels
      integer :: i,j

      Decay2bodyBaryon(:)%threshold = 0.
      Decay2bodyMeson(:)%threshold = 0.
      Decay3bodyMeson(:)%threshold = 0.

      do i=1,nDecay2bodyBaryon
        do j=1,2
          if (Decay2bodyBaryon(i)%stable(j)) then
            Decay2bodyBaryon(i)%threshold = Decay2bodyBaryon(i)%threshold + hadron(Decay2bodyBaryon(i)%ID(j))%mass
          else
            Decay2bodyBaryon(i)%threshold = Decay2bodyBaryon(i)%threshold + hadron(Decay2bodyBaryon(i)%ID(j))%minmass
          end if
        end do
      end do

      do i=1,nDecay2bodyMeson
        do j=1,2
          if (Decay2BodyMeson(i)%ID(j)==photon) then
            ! mass = 0.
          else if (Decay2bodyMeson(i)%stable(j)) then
            Decay2bodyMeson(i)%threshold = Decay2bodyMeson(i)%threshold + hadron(Decay2bodyMeson(i)%ID(j))%mass
          else
            Decay2bodyMeson(i)%threshold = Decay2bodyMeson(i)%threshold + hadron(Decay2bodyMeson(i)%ID(j))%minmass
          end if
        end do
      end do

      do i=1,nDecay3bodyMeson
        do j=1,3
          Decay3bodyMeson(i)%threshold = Decay3bodyMeson(i)%threshold + hadron(Decay3bodyMeson(i)%ID(j))%mass
        end do
      end do

    end subroutine setThresholds


    !**************************************************************************
    !****s* InitParticleProperties/setDecays
    ! NAME
    ! subroutine setDecays
    ! PURPOSE
    ! The decay ratios out of Manley & PDG are initialized.
    ! The numbering of the decay channels is according to
    ! the module decayChannels.
    !
    ! Here only the file "buuinput/DecayChannels.dat" is read in (or, if
    ! given explicitely via FileNameDecayChannels, its variant)
    !**************************************************************************
    subroutine setDecays
      use output, only: Write_ReadingInput
      use CALLSTACK, only: TRACEBACK

      integer :: ios, ID, iDummy

      call Write_ReadingInput(trim(FileNameDecayChannels),0)

      open(11,file=FileNameDecayChannels,status='UNKNOWN',iostat=ios)
      if (ios.ne.0) then
         call TRACEBACK("file '"//trim(FileNameDecayChannels)//"' not found.")
      end if
      do ID=1,pion+nMes-1
         read(11,*,iostat=ios) iDummy,hadron(ID)%decaysID, hadron(ID)%decays
         if (ios.ne.0) then
            call TRACEBACK("error while reading '"//trim(FileNameDecayChannels)//"'.")
         end if
      end do

      close(11)
      call Write_ReadingInput(trim(FileNameDecayChannels),1)

    end subroutine setDecays


    !**************************************************************************
    !****s* InitParticleProperties/updateMinmassesAndThresholds
    ! NAME
    ! subroutine updateMinmassesAndThresholds
    ! PURPOSE
    ! This routine updates the minmasses of all particles according to their
    ! decay channels (important e.g. when using rho_dilep). Then it updates
    ! the decay thresholds according to the updated minmasses.
    !**************************************************************************
    subroutine updateMinmassesAndThresholds
      use decayChannels
      integer :: ID,j,dID
      real :: minmass

      ! (1) update minmasses of mesons
      do ID=pion,pion+nMes-1
        minmass = hadron(ID)%minmass
        do j=1,nDecays
           dID = hadron(ID)%decaysID(j)
           select case (dID)
           case (0)
              exit
           case (1:) ! 2-Body-Decays
              minmass = min(minmass,Decay2BodyMeson(dID)%threshold)
           case (:-1) ! 3Body-Decays
              minmass = min(minmass,Decay3BodyMeson(-dID)%threshold)
           end select
        end do
        hadron(ID)%minmass = minmass
      end do

      ! (2) baryons
      do ID=nucleon,nbar
        minmass = hadron(ID)%minmass
        do j=1,nDecays
           dID = hadron(ID)%decaysID(j)
           select case (dID)
           case (:0)
              exit
           case (1:) ! 2-Body-Decays
              minmass = min(minmass,Decay2BodyBaryon(dID)%threshold)
           end select
        end do
        hadron(ID)%minmass = minmass
      end do

    end subroutine


    !**************************************************************************
    !****f* InitParticleProperties/initBaryon
    ! NAME
    ! type(baryonProperties) function initBaryon(name,nameTeX,mass,width,spin,rating,isospinTimes2,strangeness,charm,angMom)
    ! PURPOSE
    ! Constructor routine to simplify the initialization of a baryon.
    ! Sets also the flag %stability to a default value.
    !**************************************************************************
    function initBaryon(name,nameTeX,mass,width,spin,rating,isospinTimes2,strangeness,charm,angMom)

      type(hadronProperties) :: initBaryon ! return value

      character*(*),intent(in) :: name,nameTeX
      real,intent(in)          :: mass,width,spin
      integer,intent(in)       :: rating,isoSpinTimes2,strangeness,charm,angMom

      initBaryon%name=name
      initBaryon%nameTeX=nameTeX
      initBaryon%mass=mass
      initBaryon%width=width
      initBaryon%spin=spin
      initBaryon%rating=rating
      initBaryon%isoSpinTimes2=isoSpinTimes2
      initBaryon%strangeness=strangeness
      initBaryon%charm=charm
      initBaryon%angularMomentum=angMom
      initBaryon%propagated=.true.
      initBaryon%usedForXsections=.true.
      if (width.gt. 1e-10) then
         initBaryon%stability = 3
      else
         initBaryon%stability = 0
      end if

    end function initBaryon


    !**************************************************************************
    !****f* InitParticleProperties/initMeson
    ! NAME
    ! type(mesonProperties) function initMeson(name, nameTeX,
    ! mass, width, spin, isospinTimes2, strangeness, charm, minmass)
    ! PURPOSE
    ! Constructor routine to simplify the initialization of a meson.
    ! Sets also the flag %stability to a default value.
    !**************************************************************************
    function initMeson(name,nameTeX,mass,width,spin,isospinTimes2,strangeness,charm,minmass)

      type(hadronProperties) :: initMeson
      character*(*)      :: name
      character*(*)      :: nameTeX
      real               :: mass
      real               :: width
      real               :: spin
      integer            :: isoSpinTimes2
      integer            :: strangeness
      integer            :: charm
      real               :: minmass

      initMeson%name=name
      initMeson%nameTeX=nameTeX
      initMeson%mass=mass
      initMeson%width=width
      initMeson%spin=spin
      initMeson%isoSpinTimes2=isoSpinTimes2
      initMeson%strangeness=strangeness
      initMeson%charm=charm
      if (width.gt. 1e-10) then
         initMeson%stability = 3
      else
         initMeson%stability = 0
      end if

      initMeson%minmass = max(0.,minmass)

    end function initMeson

  end subroutine InitParticleProperties


  !****************************************************************************
  !****f* particleProperties/getAngularMomentum_meson
  ! NAME
  ! integer function getAngularMomentum_meson(DecayChannel, partID)
  ! PURPOSE
  ! Returns the angular momentum of the outgoing particles of a specific
  ! 2-body decay channel of a meson.
  ! INPUTS
  ! * integer :: DecayChannel -- The number of the 2-body decay channel of the
  !                              meson according to module decayChannels.
  ! * integer :: partID       -- The ID of the mother resonance
  ! NOTES
  ! The final state parity is given by (p_1)*(p_2)*(-1)^L
  !****************************************************************************
  function getAngularMomentum_meson(DecayChannel, partID) result (L)
    use CALLSTACK, only: TRACEBACK

    integer,intent(in) :: DecayChannel, partID
    integer :: L

    select case (DecayChannel)
    case (1,3,4,12)  ! pion Pion, kaon kaonBar , pion kaon,pion kaonBar
       ! All mesons are Spin=0 and parity -.
       ! Therefore angular momentum must equal to the spin of the resonance.
       L = int(hadron(partID)%spin)
    case (2,6)  ! rho+pi, pi+gamma
       ! L must be odd to conserve parity
       L = 1
    case (17,18)  ! pi+sigma, eta+sigma
       ! pi,eta: parity -1, sigma: parity +1 -> (-1)*(+1) = (-1)
       ! eg. eta, eta': parity -1
       L = 0
    case (5,7,8,9,10,11,13,14,15,16)
       ! those are assumed to have a constant width, therefore the angular
       ! momentum is not used and we just set it to zero
       L = 0
    case default
       write(*,*) DecayChannel,partID
       call TRACEBACK('Not yet implemented in getAngularMomentum_meson')
    end select

  end function getAngularMomentum_meson


  !****************************************************************************
  !****f* baryonWidthVacuum/getAngularMomentum_baryon
  ! NAME
  ! integer function getAngularMomentum_baryon(DecayChannel, partID)
  ! PURPOSE
  ! Returns the angular momentum of the outgoing particles of a specific
  ! 2-body decay channel of a baryon.
  ! INPUTS
  ! * integer :: DecayChannel -- The number of the 2-body decay channel of the
  !                              baryon according to module decayChannels.
  ! * integer :: partID       -- The ID of the mother resonance
  !****************************************************************************
  function getAngularMomentum_baryon(DecayChannel, partID) result (L)
    use DecayChannels, only: Decay2BodyBaryon
    use IdTable

    integer,intent(in) :: DecayChannel, partID
    integer :: L

    logical, parameter :: debug = .false.

    select case (DecayChannel)
    case (1:2,4,14)  ! Nucleon Pion, Nucleon eta , Lambda-Kaon, P11_1440 pion
       ! All baryons are Spin=1/2 and parity +
       ! All mesons are Spin=0 and parity -
       ! Therefore angular momentum should always look like in nucleon-pion.
       L = hadron(partID)%AngularMomentum
    case (13) !N sigma                  !Nucleon
       !sigma has parity + , therefore if L(pion Nucleon) =j pm 1/2
       !                             then L(sigma Nucleon)=j mp 1/2
       ! since they must differ by one to conserve parity!
       L = nint(2.*hadron(partID)%spin-hadron(partID)%AngularMomentum)
    case (3)
       ! Nucleon omega    spin 1/2, parity +  and spin = 1, parity=-
       ! J=S_omegaN + L
       ! To be improved: here we assume for simplicity that S_omegaN=1/2, then L will be the same as for pion Nucleon
       L = hadron(partID)%AngularMomentum
    case (5:12,15:18,22:26,30:34,36:40) ! Angular Momentum defined in the decay2Body field
       ! pion-Delta, rho-Nucleon, rho-Delta, pion Sigma^*, Kbar^* N, pion Lambda(1520)
       L = Decay2BodyBaryon(DecayChannel)%angularMomentum
    case (19:21,27,28) ! S=-1 decays to the lowest octet particles (pion Lambda, Kbar N, pion Sigma, eta Lambda, eta Sigma)
       L = hadron(partID)%AngularMomentum
    case (29)
       L = 3    ! Sigma(2030) --> Kbar Delta
    case (35)
       L = 1    ! Sigma(1670) --> pion Lambda(1405)
    case (41)
       L = 1    ! Sigma(2030) --> pion Lambda(1820)
    case (42)
       if (partID==lambda_2100) then
          L = 2   ! Lambda(2100) --> omega Lambda
       else if (partID==lambda_2110) then
          L = 1   ! Lambda(2110) --> omega Lambda
       else
          write(*,*) ' wrong decay channel of resonance: ', partID
          stop
       end if
    case (43:45) !s=-2 or charmed decay channels
       ! NEEDS TO BE IMPROVED
       L = 0
    case default
       write(*,*) 'Not yet implemented'
       stop
    end select

    if (debug) then
       if (L>4) then
          write(*,*) "L>4", L, "in getAngularMomentum_baryon"
          write(*,*) decayChannel, partID
          stop
       end if
    end if

  end function getAngularMomentum_baryon


  !****************************************************************************
  !****s* particleProperties/PrintParticleProperties
  ! NAME
  ! subroutine PrintParticleProperties
  ! PURPOSE
  ! Prints all particle properties and decay channels as latex tables.
  ! OUTPUT
  ! written to "GiBUU_database.tex"
  !****************************************************************************
  subroutine PrintParticleProperties
    !**************************************************************************
    !****o* particleProperties/GiBUU_database.tex
    ! NAME
    ! GiBUU_database.tex
    ! PURPOSE
    ! Lists the properties of the implemented baryons.
    !**************************************************************************
    open(10, file='GiBUU_database.tex')

    write(10,'(A)') "\documentclass[a4paper,10pt]{article}"
    write(10,'(A)') "\usepackage[left=1.5cm,right=1.5cm,top=1cm,bottom=1cm]{geometry}"
    write(10,'(A)') "\begin{document}"
    write(10,'(A)') "\pagestyle{empty}"

    call PrintMesonProperties
    call PrintBaryonProperties

    write(10,'(A)') "\end{document}"
    close(10)

  end subroutine PrintParticleProperties


  !****************************************************************************
  !****s* particleProperties/PrintBaryonProperties
  ! NAME
  ! subroutine PrintBaryonProperties
  ! PURPOSE
  ! Prints all baryon properties and decay channels as latex tables.
  !****************************************************************************
  subroutine PrintBaryonProperties
    use IdTable, only: nbar
    use DecayChannels, only: Decay2BodyBaryon

    integer :: ID, j, k, dID
    character(5), parameter :: stars='*****'
    character(14) :: name, sID

    write(*,*) 'Printing baryon properties to file'

    write(10,'(A)') "\section{Baryons}"
    write(10,'(A)') "\subsection{Baryon Properties}"

    write(10,'(A)') "\begin{tabular}{|lr|cccrccccc|} "
    write(10,'(A)') "\hline"
    write(10,'(A)') "\textbf{Name}&\textbf{ID}&\textbf{Mass}&\textbf{Width}&\textbf{Spin}&\textbf{Rating}&",  &
                    "\textbf{Isospin}&\textbf{Strange}&\textbf{Charm}&\textbf{Stability}&\textbf{min.Mass} \\"
    write(10,'(A)') "\hline"

    k=0
    do ID=1,nbar
      write(10,'("$",A14,"$ & ",I3,2("&  ",F5.3)," & ",F4.1," & ",A5,"&  ",F7.1,3("&  ",I6),"& ",F5.3,"\\")') &
            hadron(ID)%nameTeX, ID, hadron(ID)%mass, hadron(ID)%width, hadron(ID)%spin, stars(1:hadron(ID)%rating), &
            hadron(ID)%isoSpinTimes2/2., hadron(ID)%strangeness, hadron(ID)%charm, hadron(ID)%stability, hadron(ID)%minMass
      if (ID==31 .or. ID==55) write(10,'(A)') "\hline"
    end do

    write(10,'(A)') "\hline"
    write(10,'(A)') "\end{tabular}"
    write(10,'(A)') ""

    write(10,'(A)') "\subsection{Non-Strange Baryon Decays}"

    write(10,'(A)') "\begin{tabular}[t]{|lr|l|ll|l|} "
    write(10,'(A)') "\hline"
    write(10,'(A)') 'Name & ID & $\Gamma/\Gamma_{\rm tot}$ & $P_1$ & $P_2$ & L \\ '
    write(10,'(A)') "\hline"

    do ID=1,nbar
      name = hadron(ID)%nameTex
      write(sID,"(i3)") ID
      k = 0
      do j=1,nDecays
         dID = hadron(ID)%decaysID(j)
         select case (dID)
         case (0)
            exit

         case (1:) ! 2-Body-Decays
            write(10,'("$",A20,"$&",A3,"&",F8.5,"&$ ",A20,"$ &$ ",A20,"$& ",I1,"\\")') &
                 name, sID, hadron(ID)%decays(j), &
                 TeXName(Decay2BodyBaryon(dID)%id(1)), &
                 TeXName(Decay2BodyBaryon(dID)%id(2)), &
                 getAngularMomentum_baryon(dId,ID)
            name = " " ! reset for further lines
            sID = " "
            k = k + 1

         case (:-1) ! 3Body-Decays
            ! ---- not yet implemented here

         end select
      end do
      if (k>0) write(10,'(A)') "\hline"
      if (ID==16 .or. ID==32 .or. ID==46) then
         write(10,'(A)') "\end{tabular}"
         if (ID==32) then
           write(10,'(A)') " "
           write(10,'(A)') "\subsection{Strange/Charmed Baryon Decays}"
         end if
         write(10,'(A)') "\begin{tabular}[t]{|lr|l|ll|l|} "
         write(10,'(A)') "\hline"
      end if
    end do

    write(10,'(A)') "\end{tabular}"

  end subroutine PrintBaryonProperties


  !****************************************************************************
  !****s* particleProperties/PrintMesonProperties
  ! NAME
  ! subroutine PrintMesonProperties
  ! PURPOSE
  ! Prints all meson properties and decay channels as latex tables.
  !****************************************************************************
  subroutine PrintMesonProperties
    use DecayChannels, only: Decay2BodyMeson, Decay3BodyMeson

    integer       :: j, ID, dId, k
    character(14) :: name, sID

    write(*,*) 'Printing meson properties to file'

    write(10,'(A)') "\section{Mesons}"
    write(10,'(A)') "\subsection{Meson Properties}"

    write(10,'(A)') "\begin{tabular}{|lr|cccccccc|} "
    write(10,'(A)') "\hline"
    write(10,'(A)') "\textbf{Name}&\textbf{ID}&\textbf{Mass}&\textbf{Width}&\textbf{Spin}&",  &
                    "\textbf{Isospin}&\textbf{Strange}&\textbf{Charm}&\textbf{Stability}&\textbf{min.Mass} \\"
    write(10,'(A)') "\hline"

    do ID=pion,pion+nMes-1
      write(10,'("$",A14,"$  & ",I3,2(" & ",F6.4),2(" & ",F4.1),3(" & ",I2),"& ",F5.3,"\\")') &
            hadron(ID)%nameTeX, ID,  hadron(ID)%mass, hadron(ID)%width, hadron(ID)%spin, float(hadron(ID)%isoSpinTimes2)/2., &
            hadron(ID)%strangeness,  hadron(ID)%charm, hadron(ID)%stability, hadron(ID)%minMass
      if (ID==109 .or. ID==113 .or. ID==121) write(10,'(A)') "\hline"
    end do
    write(10,'(A)') "\hline"
    write(10,'(A)') "\end{tabular}"
    write(10,'(A)') ""

    write(10,'(A)') "\subsection{Meson Decays}"

    write(10,'(A)') "\begin{tabular} {|lr|l|lll|l|} "
    write(10,'(A)') "\hline"
    write(10,'(A)') 'Name & ID & $\Gamma/\Gamma_{\rm tot}$ & $P_1$ & $P_2$ & $P_3$ & L \\ '
    write(10,'(A)') "\hline"

    do ID=pion,pion+nMes-1
       name = hadron(ID)%nameTex
       write(sID,"(i3)") ID
       k = 0
       do j=1,nDecays
          dID = hadron(ID)%decaysID(j)
          select case (dID)
          case (0)
             exit

          case (1:) ! 2-Body-Decays
             write(10,'("$",A20,"$&",A3,"&",F8.4,"&$ ",A20,"$&$ ",A20,"$&$ ",A20,"$& ",I1,"\\")') &
                  name, sID, hadron(ID)%decays(j), &
                  TeXName(Decay2BodyMeson(dID)%id(1)), &
                  TeXName(Decay2BodyMeson(dID)%id(2)), &
                  " ", getAngularMomentum_meson(dId,ID)
             name = " "
             sID = " "
             k = k + 1
          case (:-1) ! 3Body-Decays
             write(10,'("$",A20,"$&",A3,"&",F8.4,"&$ ",A20,"$&$ ",A20,"$&$ ",A20,"$& ",A,"\\")') &
                  name, sID, hadron(ID)%decays(j), &
                  TeXName(Decay3BodyMeson(-dID)%id(1),Decay3BodyMeson(-dID)%charge(1),.false.), &
                  TeXName(Decay3BodyMeson(-dID)%id(2),Decay3BodyMeson(-dID)%charge(2),.false.), &
                  TeXName(Decay3BodyMeson(-dID)%id(3),Decay3BodyMeson(-dID)%charge(3),.false.), &
                  "---"
             name = " "
             sID = " "
             k = k + 1
          end select
       end do
       if (k>0) write(10,'(A)') "\hline"
    end do

    write(10,'(A)') "\end{tabular}"

  end subroutine PrintMesonProperties


  !****************************************************************************
  !****s* particleProperties/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Read out job card to initialize module parameters.
  !****************************************************************************
  subroutine init
    use inputGeneral, only: path_To_Input
    use output, only: Write_ReadingInput
    integer :: ios

    !**************************************************************************
    !****n* particleProperties/initDatabase
    ! NAME
    ! NAMELIST initDatabase
    ! PURPOSE
    ! Includes the switches:
    ! * propagationSwitch
    ! * usageForXsectionSwitch
    ! * rho_dilep
    ! * FileNameDecayChannels
    !**************************************************************************
    NAMELIST /initDatabase/ propagationSwitch, usageForXsectionSwitch, &
                            rho_dilep, FileNameDecayChannels

    call Write_ReadingInput('initDatabase',0)
    rewind(5)
    read(5,nml=initDatabase,IOSTAT=ios)
    call Write_ReadingInput('initDatabase',0,ios)

    write(*,*) 'propagationSwitch      : ', propagationSwitch
    write(*,*) 'usageForXsectionSwitch : ', usageForXsectionSwitch
    write(*,*) 'rho_dilep              : ', rho_dilep

    if (len_trim(FileNameDecayChannels)>0) then
       if (index(FileNameDecayChannels,"/")>0) then
          FileNameDecayChannels = trim(FileNameDecayChannels)
       else
          FileNameDecayChannels = trim(path_to_Input)//'/'//trim(FileNameDecayChannels)
       end if
       write(*,*) 'DecayChannels.dat      : ', trim(FileNameDecayChannels)
    else
       FileNameDecayChannels = trim(path_to_Input)//'/DecayChannels.dat'
       write(*,*) 'DecayChannels.dat      : -- hardwired --'
    end if

    call Write_ReadingInput('initDatabase',1)

  end subroutine init


  !****************************************************************************
  !****s* particleProperties/ReadInitModify
  ! NAME
  ! subroutine ReadInitModify
  ! PURPOSE
  ! Read in namelist "ModifyParticles" to initialize module parameters.
  !****************************************************************************
  subroutine ReadInitModify
    use output, only: Write_ReadingInput
    use IdTable, only: pion, nMes

    character(60) :: formI = '("  ",A," (ID=",i3,") changed: ",A," = ",i3," --> ",i3)'
    character(60) :: formR = '("  ",A," (ID=",i3,") changed: ",A," = ",f5.3," --> ",f5.3)'

    !**************************************************************************
    !****g* ReadInitModify/mass
    ! SOURCE
    real, dimension(1:pion+nMes-1) :: mass = -1.0
    ! PURPOSE
    ! Input array for modifications on the particle mass
    ! NOTES
    ! This array is intended to "input" values for the mass of the particles,
    ! which are different from the default. Therefore only entries, which are
    ! positive after reading the file are stored in the internal database.
    !**************************************************************************

    !**************************************************************************
    !****g* ReadInitModify/width
    ! SOURCE
    real, dimension(1:pion+nMes-1) :: width = -1.0
    ! PURPOSE
    ! Input array for modifications on the particle width
    ! NOTES
    ! This array is intended to "input" values for the width of the particles,
    ! which are different from the default. Therefore only entries, which are
    ! positive after reading the file are stored in the internal database.
    !**************************************************************************

    !**************************************************************************
    !****g* ReadInitModify/stabilityFlag
    ! SOURCE
    integer, dimension(1:pion+nMes-1) :: stabilityFlag = -1
    ! PURPOSE
    ! Input array for modifications on the particle stability
    ! NOTES
    ! This array is intended to "input" values for the stability of the
    ! particles, which are different from the default. Therefore only entries,
    ! which are >-1 after reading the file are stored in the internal database.
    !
    ! The index of the array is the particle ID. The value encodes on a bitwise
    ! level, how the particle may decay (cf. also master_1Body):
    ! * 1: particle may decay during run, if Gamma > gammaCutOff
    ! * 2: particle may decay at the end of the run, if Gamma > 0.
    ! * 4: particle may decay at the end via Jetset, if there the parameters allow for a decay.
    ! The default values are one of the following:
    ! * 0: particle may not decay at all (i.e. it is stable)
    ! * 3: particle may decay both during run and at the end (combination of 1 and 2)
    !**************************************************************************

    !**************************************************************************
    !****n* particleProperties/ModifyParticles
    ! NAME
    ! NAMELIST ModifyParticles
    ! PURPOSE
    ! Includes the switches:
    ! * mass
    ! * width
    ! * stabilityFlag
    !**************************************************************************
    NAMELIST /ModifyParticles/ mass, width, stabilityFlag

    integer :: ios,i,iOld
    real :: rOld

    call Write_ReadingInput('ModifyParticles',0)
    rewind(5)
    read(5,nml=ModifyParticles,IOSTAT=ios)
    call Write_ReadingInput('ModifyParticles',0,ios)

    do i=1,pion+nMes-1

       if (mass(i)>=0.) then
          rOld = hadron(i)%mass
          hadron(i)%mass = mass(i)
          if (rOld/=mass(i)) write(*,formR) TRIM(hadron(i)%name),i,'mass ',rOld,mass(i)
       end if

       if (width(i)>=0.) then
          rOld = hadron(i)%width
          hadron(i)%width = width(i)
          if (rOld/=width(i)) write(*,formR) TRIM(hadron(i)%name),i,'width',rOld,width(i)
       end if

       if (stabilityFlag(i)>-1) then
          iOld = hadron(i)%stability
          hadron(i)%stability = stabilityFlag(i)
          if (iOld/=stabilityFlag(i)) write(*,formI) TRIM(hadron(i)%name),i,'stability',iOld,stabilityFlag(i)
       end if

    end do

!     call TestOldVersion()

    call Write_ReadingInput('ModifyParticles',1)

!   contains
!
!     subroutine TestOldVersion
!       integer, dimension(1:pion+nMes-1) :: stabilityFlag = -1
!       integer :: ios
!       NAMELIST /initStability/ stabilityFlag
!       rewind(5)
!       read(5,nml=initStability,IOSTAT=IOS)
!
!       if (ios<0) return ! okay, namelist not in jobcard
!
!       write(*,*) 'Namelist "initStability" found in jobCard!'
!       write(*,*) 'Functionality has moved to "ModifyParticles".'
!       write(*,*) 'Stop.'
!       stop
!
!     end subroutine TestOldVersion


  end subroutine ReadInitModify


  !****************************************************************************
  !****f* particleProperties/isStrange
  ! NAME
  ! logical function isStrange(ID)
  ! INPUTS
  ! * integer :: ID
  ! RESULT
  ! * .true. if ID corresponds to a strange particle
  !****************************************************************************
  logical function isStrange(ID)
    use IdTable, only: isHadron
    integer, intent(in) :: ID
    if (isHadron(ID)) then
      isStrange = (hadron(ID)%strangeness/=0)
    else
      isStrange = .false.
    end if
  end function isStrange


  !****************************************************************************
  !****f* particleProperties/isCharmed
  ! NAME
  ! logical function isCharmed(ID)
  ! INPUTS
  ! * integer :: ID
  ! RESULT
  ! * .true. if ID corresponds to a charmed particle
  !****************************************************************************
  logical function isCharmed(ID)
    integer, intent(in) :: ID
    isCharmed = (hadron(ID)%charm/=0)
  end function isCharmed


  !****************************************************************************
  !****f* particleProperties/isNonExotic
  ! NAME
  ! logical function isNonExotic(ID)
  ! INPUTS
  ! * integer :: ID
  ! RESULT
  ! * .true. if ID corresponds to a non-exotic particle
  !****************************************************************************
  logical function isNonExotic(ID)
    integer, intent(in) :: ID
    isNonExotic=(.not.(isStrange(ID).or.isCharmed(ID)))
  end function isNonExotic


  !****************************************************************************
  !****f* particleProperties/validCharge
  ! NAME
  ! logical function validCharge (teilchen)
  ! PURPOSE
  ! This routine returns .true. if the charge is valid for the particle,
  ! and .false if not.
  ! This is done by checking the z-component of its isospin.
  !   Q=I_Z+(B+S+C+B+T)/2.
  !   => I_Z=Q-(B+S+C+B+T)/2.
  !   => I_Z=Q-(B+S+C)/2. if no bottom and top quarks are considered
  ! INPUTS
  ! type(particle), intent(in) :: teilchen   ! Regarded particle
  !****************************************************************************
  logical function validCharge (teilchen)
    use particleDefinition
    use IdTable, only: photon, isMeson, isBaryon
    use CALLSTACK, only: TRACEBACK

    type(particle), intent(in) :: teilchen

    real :: Iz, max_Iz
    integer :: Q
    real , parameter :: eps=0.001

    if (.not.teilchen%antiParticle) then
       Q=teilchen%charge
    else ! Antiparticle
       Q=-teilchen%charge
    end if

    if (isMeson(teilchen%id)) then
       Iz=Q-(hadron(teilchen%id)%strangeness+hadron(teilchen%id)%charm)/2.
    else if (isBaryon(teilchen%id)) then
       Iz=Q-(1+hadron(teilchen%id)%strangeness+hadron(teilchen%id)%charm)/2.
    else
       if (teilchen%Id==0 .or. teilchen%Id==photon) then
         validCharge=.true.
         return
       end if
       write(*,*) 'Error in validCharge: Particle is no meson and no baryon:', &
                  teilchen%id, teilchen%antiparticle, teilchen%charge
       call Traceback('stop')
    end if

    max_Iz = hadron(teilchen%id)%isoSpinTimes2/2.

    validCharge = (abs(Iz)<(max_Iz+eps))

  end function validCharge


  !****************************************************************************
  !****f* particleProperties/validCharge_ID
  ! NAME
  ! logical function validCharge_ID (ID, Q)
  ! PURPOSE
  ! This routine returns .true. if the charge is valid for the particle,
  ! and .false if not.
  ! This is done by checking the z-component of its isospin.
  !   Q=I_Z+(B+S+C+B+T)/2.
  !   => I_Z=Q-(B+S+C+B+T)/2.
  !   => I_Z=Q-(B+S+C)/2. if no bottom and top quarks are considered
  ! NOTES
  ! We assume that Q is the charge of the particle, not the antiparticle.
  ! INPUTS
  ! * integer, intent(in) :: ID  ! ID of Regarded particle
  ! * integer, intent(in) :: Q   ! Charge of Regarded particle
  !****************************************************************************
  logical function validCharge_ID (ID, Q)
    use IdTable, only: photon, isMeson, isBaryon
    use CALLSTACK, only: TRACEBACK

    integer, intent(in) :: ID, Q

    real :: Iz, max_Iz
    real , parameter :: eps=0.001

    if (isMeson(id)) then
       Iz=Q-(hadron(id)%strangeness+hadron(id)%charm)/2.
    else if (isBaryon(id)) then
       Iz=Q-(1+hadron(id)%strangeness+hadron(id)%charm)/2.
    else
       if (Id==0 .or. Id==photon) then
         validCharge_ID=.true.
         return
       end if
       write(*,*) 'Error in validCharge: Particle is no meson and no baryon:', id, Q
       call Traceback('stop')
    end if

    max_Iz = hadron(id)%isoSpinTimes2/2.

    validCharge_ID = (abs(Iz)<(max_Iz+eps))

  end function validCharge_ID


  !****************************************************************************
  ! cf. interface "PartName" :
  !****************************************************************************
  character(15) function PartName1 (ID, IQ, isAnti)
    use IdTable, only: isBaryon, isMeson, getAntiMeson, electron, muon, tau, Wboson
    integer, intent(in) :: ID,IQ
    logical, intent(in) :: isAnti

    integer :: iID, iIQ
    character*(*), dimension(-2:2),parameter :: NN = (/"--","- ","0 ","+ ", "++"/)

    if (isBaryon(ID)) then
      if (isAnti) then
         PartName1 = TRIM(hadron(ID)%name)//'~'//TRIM(NN(IQ))
      else
         PartName1 = TRIM(hadron(ID)%name)//TRIM(NN(IQ))
      end if
    else if (isMeson(ID)) then
      if (isAnti) then
         call getAntiMeson(ID,IQ, iID,iIQ)
      else
         iID = ID
         iIQ = IQ
      end if
      PartName1 = TRIM(hadron(iID)%name)//TRIM(NN(iIQ))
    else
      select case (ID)
      case (electron,muon,tau,Wboson)
         PartName1 = trim(PartName3(ID))//TRIM(NN(IQ))
      case default
         PartName1 = trim(PartName3(ID))
      end select
      if (isAnti) PartName1 = trim(PartName1) // "~"
    end if

  end function PartName1
  !-------------------------------------------------------------------------
  character(15) function PartName2 (Part)
    use particleDefinition
    type(particle), intent(in) :: Part
    PartName2 = PartName1 (Part%ID, Part%charge, Part%antiparticle)
  end function PartName2
  !-------------------------------------------------------------------------
  character(15) function PartName3 (ID)
    use IdTable
    integer, intent(in) :: ID
    if (isHadron(ID)) then
      PartName3 = TRIM(hadron(ID)%name)
    else
      select case (ID)
      case (electron)
         PartName3 = "e"
      case (muon)
         PartName3 = "mu"
      case (tau)
         PartName3 = "tau"
      case (electronNeutrino)
         PartName3 = "nu_e"
      case (muonNeutrino)
         PartName3 = "nu_mu"
      case (tauNeutrino)
         PartName3 = "nu_tau"
      case (Wboson)
         PartName3 = "W"
      case (Zboson)
         PartName3 = "Z0"
      case (photon)
         PartName3 = "gamma"
      case default
         PartName3 = "?????"
      end select
    end if
  end function PartName3


  !****************************************************************************
  ! cf. interface "TeXName" :
  !****************************************************************************
  character(25) function TeXName1 (ID, IQ, isAnti)
    use IdTable
    integer, intent(in) :: ID,IQ
    logical, intent(in) :: isAnti

    integer :: iID, iIQ
    character*(*), dimension(-2:2),parameter :: NN = (/"^{--}","^-   ","^0   ","^+   ", "^{++}"/)

    TeXName1 = 'XXX'
    if (ID<=0) return

    if (isBaryon(ID)) then

       if (isAnti) then
          TeXName1 = '\overline{'//TRIM(hadron(ID)%nameTeX)//'}'//TRIM(NN(IQ))
       else
          TeXName1 = TRIM(hadron(ID)%nameTeX)//TRIM(NN(IQ))
       end if

    else if (isMeson(ID)) then

       if (isAnti) then
          call getAntiMeson(ID,IQ, iID,iIQ)
       else
          iID = ID
          iIQ = IQ
       end if

       TeXName1 = TRIM(hadron(iID)%nameTeX)//TRIM(NN(iIQ))

    else
       select case (ID)
       case (electron)
          TeXName1 = "e"//TRIM(NN(IQ))
       case (muon)
          TeXName1 = "\mu"//TRIM(NN(IQ))
       case (tau)
          TeXName1 = "\tau"//TRIM(NN(IQ))

       case (electronNeutrino)
          if (isAnti) then
             TeXName1 = "\overline{\nu}_e"
          else
             TeXName1 = "\nu_e"
          end if
       case (muonNeutrino)
          if (isAnti) then
             TeXName1 = "\overline{\nu}_\mu"
          else
             TeXName1 = "\nu_\mu"
          end if
       case (tauNeutrino)
          if (isAnti) then
             TeXName1 = "\overline{\nu}_\tau"
          else
             TeXName1 = "\nu_\tau"
          end if

       case (Wboson)
          TeXName1 = "W"//TRIM(NN(IQ))
       case (Zboson)
          TeXName1 = "Z0"
       case (photon)
          TeXName1 = '\gamma'

       case default
          TeXName1 = '?????'
       end select
    end if

  end function TeXName1
  !-------------------------------------------------------------------------
  character(25) function TeXName2 (Part)
    use particleDefinition
    type(particle), intent(in) :: Part
    TeXName2 = TeXName1(Part%ID,Part%charge,Part%antiparticle)
  end function TeXName2
  !-------------------------------------------------------------------------
  character(25) function TeXName3 (ID)
    use IdTable, only: isHadron, photon
    integer, intent(in) :: ID

    TeXName3 = 'XXX'
    if (ID<=0) return

    if (isHadron(ID)) then
       TeXName3 = TRIM(hadron(ID)%nameTeX)
    else
       select case (ID)
       case (photon)
          TeXName3 = '\gamma'
       end select
    end if

  end function TeXName3


End Module ParticleProperties
