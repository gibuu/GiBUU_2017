!******************************************************************************
!****m* /Dilepton_Analysis
! NAME
! module Dilepton_Analysis
! PURPOSE
! This analysis module produces various dilepton spectra.
! Currently only di-electrons are supported.
!******************************************************************************
module Dilepton_Analysis

  use histMC, only: histogramMC
  use nucleusDefinition, only: tNucleus

  implicit none
  private

  !****************************************************************************
  !****g* Dilepton_Analysis/Enable
  ! SOURCE
  logical, save :: Enable = .false.
  ! PURPOSE
  ! If .true. the dilepton analysis will be performed, otherwise not.
  !****************************************************************************

  !****************************************************************************
  !****g* Dilepton_Analysis/Extra
  ! SOURCE
  logical, save :: Extra = .false.
  ! PURPOSE
  ! If .true. an extended analysis will be performed, writing out many extra
  ! histograms (beyond the basic ones: mass, pT and rapidity).
  !****************************************************************************

  !****************************************************************************
  !****g* Dilepton_Analysis/DeltaDalitz
  ! SOURCE
  integer, save :: DeltaDalitz = 2
  ! PURPOSE
  ! Choose between different parametrizations of the Delta Dalitz decay width
  ! (Delta -> N e+e-):
  ! * 1 = Wolf, http://inspirehep.net/record/306273
  ! * 2 = Krivoruchenko (default), http://inspirehep.net/record/555421
  ! * 3 = HadronTensor
  ! * 4 = Ernst, http://inspirehep.net/record/452782
  !****************************************************************************

  !****************************************************************************
  !****g* Dilepton_Analysis/DeltaDalitzFF
  ! SOURCE
  integer, save :: DeltaDalitzFF = 1
  ! PURPOSE
  ! Choose a parametrization of the electromagnetic N-Delta transition form
  ! factor for the Delta Dalitz decay (only used for DeltaDalitz = 2):
  ! * 1 = constant (default)
  ! * 2 = Dipole
  ! * 3 = MAID 2005
  ! * 4 = simple VMD
  ! * 5 = Wan/Iachello, Int. J. Mod. Phys. A 20 (2005) 1846, http://inspirehep.net/record/689265
  ! * 6 = Ramalho/Pena, Phys.Rev. D85 (2012) 113014, http://inspirehep.net/record/1114321
  !****************************************************************************

  !****************************************************************************
  !****g* Dilepton_Analysis/omegaDalitzFF
  ! SOURCE
  integer, save :: omegaDalitzFF = 1
  ! PURPOSE
  ! Choose between different parametrizations of the omega Dalitz decay
  ! (omega -> pi^0 e+e-) form factor:
  ! * 0 = constant
  ! * 1 = Effenberger/Bratkovskaya (default)
  ! * 2 = standard VMD
  ! * 3 = Terschluesen/Leupold
  !****************************************************************************

  !****************************************************************************
  !****g* Dilepton_Analysis/b_pi
  ! SOURCE
  real, save :: b_pi = 5.5
  ! PURPOSE
  ! This constant represents the b parameter in the form factor
  ! of the pi0 Dalitz decay (in GeV^-2), cf. Effenberger Diss. eq. (2.141).
  ! Originally taken from L.G. Landsberg, Phys. Rep. 128, 301 (1985).
  !****************************************************************************

  !****************************************************************************
  !****g* Dilepton_Analysis/lambda_eta
  ! SOURCE
  real, save :: lambda_eta = 0.716
  ! PURPOSE
  ! This constant represents the Lambda parameter in the form factor of the
  ! eta Dalitz decay in GeV.
  ! Values:
  ! * L.G. Landsberg, Phys. Rep. 128, 301 (1985): Lambda = 720 MeV
  ! * HADES pp@2.2, B. Spruck, Diss. (2008):      Lambda = 676 MeV
  ! * NA60, Arnaldi et al, PLB 677 (2009):        Lambda = 716 MeV (default)
  ! * CB/TAPS, Berghäuser et al, PLB 701 (2011):  Lambda = 722 MeV
  !****************************************************************************

  !****************************************************************************
  !****g* Dilepton_Analysis/etaPrimeDalitzFF
  ! SOURCE
  integer, save :: etaPrimeDalitzFF = 0
  ! PURPOSE
  ! Choose between different parametrizations of the eta' Dalitz decay
  ! (eta' -> e+e- gamma) form factor:
  ! * 0 = constant (default)
  ! * 1 = eta FF (cf. lambda_eta)
  ! * 2 = generic VMD
  ! * 3 = Genesis / Lepton-G
  ! * 4 = standard VMD (Terschluesen)
  !****************************************************************************

  !****************************************************************************
  !****g* Dilepton_Analysis/angDist
  ! SOURCE
  integer, save :: angDist = 1
  ! PURPOSE
  ! This switch determines the angular distribution of the pseudoscalar
  ! Dalitz decays P -> e+ e- gamma (with P=pi0,eta,etaPrime):
  ! * 0 = isotropic decay
  ! * 1 = anisotropic decay according to 1 + B*cos**2(theta) with B=1
  ! * 2 = the Dalitz decays of pi0 and eta will be done via Pythia.
  !****************************************************************************

  !****************************************************************************
  !****g* Dilepton_Analysis/brems
  ! SOURCE
  integer, save :: brems = 1
  ! PURPOSE
  ! This switch determines how the bremsstrahlung contribution is obtained:
  ! * 0 = none
  ! * 1 = soft-photon approximation (SPA)
  ! * 2 = according to the one-boson-exchange (OBE) model by R. Shyam
  !       (for NN bremstrahlung only, no em. form factors)
  ! * 3 = as 2, but with pion em. form factor (for pn)
  !****************************************************************************


  !****************************************************************************
  !****g* Dilepton_Analysis/nEvent
  ! SOURCE
  integer, save :: nEvent = 10
  ! PURPOSE
  ! Number of events to generate for each dilepton decay
  ! (to enhance statistics).
  !****************************************************************************

  !****************************************************************************
  !****g* Dilepton_Analysis/nEvent_BH
  ! SOURCE
  integer, save :: nEvent_BH = 1000
  ! PURPOSE
  ! Number of events for Bethe-Heitler simulation.
  ! BH typically needs a lot more statistics than the other dilepton channels.
  ! Therefore nEvent_BH should be much bigger than nEvent.
  !****************************************************************************


  !****************************************************************************
  !****g* Dilepton_Analysis/kp_cut
  ! SOURCE
  logical, save :: kp_cut = .false.
  ! PURPOSE
  ! Perform a cut on (k*p) in the dilepton analysis,
  ! where k is the photon 4-momentum, and p is the electron or positron
  ! 4-momentum.
  ! This is useful for suppressing the BH contribution. Cf. "kp_min".
  !****************************************************************************


  !****************************************************************************
  !****g* Dilepton_Analysis/kp_min
  ! SOURCE
  real, save :: kp_min = 0.01
  ! PURPOSE
  ! If kp_cut=.true. a cut on (k*p) is performed.
  ! kp_min determines the position of this cut.
  ! Only events with (k*p)>kp_min are taken into account.
  !****************************************************************************


  !****************************************************************************
  !****g* Dilepton_Analysis/filter
  ! SOURCE
  integer, save :: filter = 0
  ! PURPOSE
  ! If filter is nonzero, a filtering algorithm will be applied to the
  ! dilepton pairs, otherwise they will be written to the histograms
  ! unfiltered. For details on the filtering parameters see routine 'CS'.
  ! Choices:
  ! * 0 = no filter
  ! * 1 = DLS
  ! * 2 = HADES (simple cuts on polar angle, absolute momentum and opening angle)
  ! * 3 = HADES (full acceptance filter, using pair acceptance)
  ! * 4 = HADES (full acceptance filter, using single-particle acceptance)
  ! * 5 = g7/CLAS @ JLab
  ! * 6 = KEK E325 (cuts on rapidity, transverse momentum and opening angle)
  ! * 7 = JPARC E16
  !
  ! NOTES
  ! For filtering modes 3 and 4, the file containing the acceptance
  ! matrices must be specified (cf. hadesFilterFile).
  !****************************************************************************


  !****************************************************************************
  ! Detector resolution in GeV for the detector setups defined above
  ! (cf. "filter"). This is used to smear the mass spectrum with a Gaussian of
  ! the given width. Currently only used for filter #6 (KEK, 11 MeV)
  ! and #7 (JPARC, 5 MeV).
  real, parameter :: detector_res(6:7) = (/ 0.011, 0.005 /)
  !****************************************************************************


  !****************************************************************************
  !****g* Dilepton_Analysis/hadesFilterFile
  ! SOURCE
  character(len=200), save :: hadesFilterFile = ""
  ! PURPOSE
  ! This character string determines the location of the file containing
  ! the HADES acceptance matrices (filename with absolute or relative path).
  ! It has to be set for filtering modes 3 and 4.
  !****************************************************************************


  !****************************************************************************
  !****g* Dilepton_Analysis/p_lep_min
  ! SOURCE
  real, save :: p_lep_min = 0.
  ! PURPOSE
  ! This switch sets a lower bound on the lepton momentum.
  ! Only leptons with momenta larger than this threshold will pass the filter.
  ! This switch is only used for filtering mode 5 (JLab).
  !****************************************************************************


  !****************************************************************************
  !****g* Dilepton_Analysis/beta_gamma_cut
  ! SOURCE
  real, save :: beta_gamma_cut = 1.25
  ! PURPOSE
  ! This is an upper bound on the beta*gamma value of the lepton pair.
  ! Since beta*gamma = p/m, it cuts on slow pairs.
  !****************************************************************************


  !****************************************************************************
  !****g* Dilepton_Analysis/binsz
  ! SOURCE
  real, save :: binsz = 0.01
  ! PURPOSE
  ! This determines the bin size of the dilepton mass spectrum in GeV.
  ! Default is 10 MeV.
  !****************************************************************************


  !****************************************************************************
  !****g* Dilepton_Analysis/WriteEvents
  ! SOURCE
  integer, save :: WriteEvents = 0
  ! PURPOSE
  ! This switch decides whether we write out the simulated events.
  ! Possible values:
  ! * 0: Don't write events (default).
  ! * 1: We write out only the lepton pair information (including charge,
  !   four-momentum, perturbative weight, source channel and filter result).
  !   All of this will be written to a file called 'Dilepton_Events.dat'.
  ! * 2: As 1, but only writing exclusive events (NN->NNe+e-).
  ! * 3: We write out all produced particles in the event
  !   (including the lepton pair) to a file called 'Dilepton_FullEvents.dat'.
  ! * 4: As 3, but only writing exclusive events (NN->NNe+e-).
  ! * 5: As 4, but only writing out R->Ne+e- events (with R=N*,Delta*).
  !****************************************************************************


  !****************************************************************************
  !****g* Dilepton_Analysis/massBinning
  ! SOURCE
  real, dimension(1:4), save :: massBinning = (/ 0.150, 0.550, 9.999, 9.999 /)
  ! PURPOSE
  ! We produce several histograms (e.g. p,pT,mT,y,theta_cm) not only
  ! for the full mass range, but also for (up to 5) different mass bins.
  ! The borders of these bins are given by this array.
  !****************************************************************************


  !****************************************************************************
  !****g* Dilepton_Analysis/particle_source
  ! SOURCE
  logical, save :: particle_source = .true.
  ! PURPOSE
  ! This switch determines whether the mass spectrum will contain
  ! separate channels for different sources of particles,
  ! such as decays (R->rho N) or collisions (pi pi -> rho, K K -> phi).
  ! Currently this is only done for the rho and phi mesons.
  ! Note: If using this switch, the "sum" channel in the mass histogram should
  ! not be used, since the rho and phi contributions will enter twice.
  !****************************************************************************

  !logical, save :: CoulombCorrection=.false. ! Make Coulomb correction for e+e- pair

  type(histogramMC),save :: msigma,msigma_slow,msigma_fast,msigma_bgcut
  type(histogramMC),save :: msigma_pri,msigma_sec,ptsigma_pri,ptsigma_sec
  type(histogramMC),save :: esigma,bgsigma,oasigma,dsigma,partnum,partnum_noform,vm_mass
  ! spectra for different mass bins: 0=total; 1:5=different mass bins
  type(histogramMC),dimension(0:5),save :: psigma,ptsigma,mtsigma,ysigma,tcmsigma,helsigma
  type(histogramMC), save :: deltaMass                            ! mass distribution of the Delta
  type(histogramMC), save :: rhoMass                            ! mass distribution of the rho

  logical, save :: init = .true.
  real, save :: ProjectileEnergy = 0.
  real, save :: ProjectileEnergy_max = 0.
  integer, save :: NumberProjectiles = 1
  real, save :: ProjectileMomentum(0:3) = 0.
  real, dimension(1:3), save :: betaCM2Lab = 0.           ! boost vector from CM to Lab frame
  real, dimension(1:3), save :: betaTarg2Proj = 0.        ! boost vector from target to projectile
  real, dimension(1:8), save :: inclXS = 0.               ! inclusive particle production cross sections
  logical, save :: boostToLab = .false.                   ! do we have to boost to the lab frame?
  type(tNucleus), pointer, save :: targNuc
  real, dimension(1:2), save :: plep_cut

  ! array for compressing baryon channels, keeping only those that have a rho decay (others get a zero)
  integer, parameter :: bar_ch(1:31) = (/  0,  1,  0,  2,  3,  4,  5,  6,  7,  0, &
                                           8,  9, 10, 11, 12, 13, 14,  0, 15, 16, &
                                          17, 18,  0,  0,  0,  0,  0,  0, 19, 20, &
                                           0 /)

  public :: Dilep_Init,Dilep_UpdateProjectile,Dilep_Decays,Dilep_Brems,Dilep_BH,Dilep_write_CS

  public :: dGamma_dM_DeltaDalitz, DeltaWidth_gammaN
  public :: Gamma0_DeltaDalitz_Wolf, Gamma0_DeltaDalitz_Krivo, Gamma0_DeltaDalitz_Ernst
  public :: FF_omega_effe, FF_VMD, FF_omega_terschluesen
  public :: FF_eta, FF_etaprime_genesis, FF_etaprime_terschluesen
  public :: Delta_FF_MAID, Delta_FF_VMD, Delta_FF_Iachello, Delta_FF_Dipole

contains


  !****************************************************************************
  !****s* Dilepton_Analysis/Dilep_Init
  ! NAME
  ! subroutine Dilep_Init (Ekin_max, N, cms, beta)
  ! PURPOSE
  ! This callback routine passes some eventtype-specific information to the Dilepton_Analysis module.
  ! INPUTS
  ! * real, intent(in) :: Ekin_max                       - maximum energy to be expected during the simulation
  ! * integer, intent(in), optional :: N                 - number of projectile particles
  ! * logical, intent(in), optional :: cms               - is the collision performed in the center-of-mass frame?
  ! * real, dimension(1:3), intent(in), optional :: beta - boost vector to the laboratory frame
  !****************************************************************************
  subroutine Dilep_Init (Ekin_max, N, cms, beta)
    use inputGeneral, only: path_To_Input
    use HAFT_single, only: setFileName
    use HAFT_pair, only: setPairFileName

    real, intent(in) :: Ekin_max
    integer, intent(in), optional :: N
    logical, intent(in), optional :: cms
    real, dimension(1:3), intent(in), optional :: beta

    logical :: ex

    ProjectileEnergy_max = Ekin_max
    if (present(N)) NumberProjectiles = N
    if (present(cms)) boostToLab = cms
    if (present(beta)) betaCM2Lab = beta

    if (init) call readInput

    if (.not. Enable) return

    write(*,*) "Initializing dilepton analysis:"
    write(*,'(A,G12.5)')  "Ekin = ", Ekin_max
    write(*,'(A,I5)')     "number of projectiles = ", NumberProjectiles
    write(*,'(A,3G12.5)') "betaCM2Lab = ", betaCM2Lab
    write(*,'(A,L1)')      "boostToLab = ", boostToLab

    if (filter==4 .and. hadesFilterFile=="") then
      ! auto-choose filter file
      select case (targNuc%charge)
      case (1)  ! hydrogen
        select case (targNuc%mass)
        case (1)  ! proton
          if (abs(Ekin_max-3.5)<1E-3) then
            hadesFilterFile = trim(path_To_Input) // "/hades/HadesSingleAcc-p35p-APR07-effgt5-RK-v2.acc"
          else if (abs(Ekin_max-2.2)<1E-3) then
            hadesFilterFile = trim(path_To_Input) // "/hades/HadesSingleAcc-p22p-JAN04-effgt5-RK-v2.acc"
          else if (abs(Ekin_max-1.25)<1E-3) then
            hadesFilterFile = trim(path_To_Input) // "/hades/HadesSingleAcc-p125p-APR06-effgt5-RK-v2.acc"
          end if
        case (2)  ! deuteron
          if (abs(Ekin_max-1.25)<1E-3) then
            hadesFilterFile = trim(path_To_Input) // "/hades/HadesSingleAcc-d125p-MAY07-effgt5-RK-v2.acc"
          end if
        end select
      case (6)  ! carbon
        if (abs(Ekin_max-1.0)<1E-3) then
          hadesFilterFile = trim(path_To_Input) // "/hades/HadesSingleAcc-c10c-AUG04-effgt5-KP-v2.acc"
        else if (abs(Ekin_max-2.0)<1E-3) then
          hadesFilterFile = trim(path_To_Input) // "/hades/HadesSingleAcc-c20c-NOV02-effgt0-KP-v1.acc"
        end if
      case (17,18,19)  ! Cl,Ar,K
        if (abs(Ekin_max-1.756)<5E-3) then
          hadesFilterFile = trim(path_To_Input) // "/hades/HadesSingleAcc-ar17kcl-SEP05-effge5-HC-exp-flat-v1.acc"
        end if
      case (41)  ! niobium
        if (abs(Ekin_max-3.5)<1E-3) then
          hadesFilterFile = trim(path_To_Input) // "/hades/HadesSingleAcc-p35Nb-SEP08-effgt5-MVA-v1.acc"
        end if
      end select
    end if

    if (filter==3 .or. filter==4) then
      ! HADES filter: check if acceptance file exists
      Inquire(file=hadesFilterFile,exist=ex)
      if (.not. ex) then
        write(*,*) "Error: HADES filter file does not exist! ",hadesFilterFile
        stop
      end if

      select case (filter)
      case (3)
        write(*,*) "reading HADES pair acceptance from ",hadesFilterFile
        call setPairFileName(hadesFilterFile)
      case (4)
        write(*,*) "reading HADES single-particle acceptance from ",hadesFilterFile
        call setFileName(hadesFilterFile)
      end select
    end if

    if (filter == 4) then
      ! set up lepton momentum cuts for HADES
      if (abs(Ekin_max-3.5)<1E-3) then
         if (targNuc%charge == 41) then
            ! pNb@3.5
            plep_cut = (/ 0.1, 2.0 /)
         else
            ! pp@3.5
            plep_cut = (/ 0.08, 2.0 /)
         end if
      else if (abs(Ekin_max-2.2)<1E-3) then
        ! pp@2.2
        plep_cut = (/ 0.1, 2.0 /)
      else if (abs(Ekin_max-1.756)<5E-3 .or. targNuc%charge == 79 .or. index(hadesFilterFile,"ar17kcl")>0) then
        ! ArKCl@1.76, AuAu@1.23
        plep_cut = (/ 0.1, 1.1 /)
      else
        ! pp@1.25, np@1.25, CC@1, CC@2
        plep_cut = (/ 0.05, 1.8 /)
      end if
      write(*,'(A,2G12.5)') "lepton momentum cut = ", plep_cut
    end if

    write(*,*)
  end subroutine Dilep_Init


  subroutine Dilep_UpdateProjectile (Ekin, mom)
    real, intent(in) :: Ekin, mom(0:3)
    ProjectileEnergy = Ekin
    ProjectileMomentum = mom
    betaTarg2Proj = ProjectileMomentum(1:3)/ProjectileMomentum(0)
  end subroutine


  !****************************************************************************
  !****s* Dilepton_Analysis/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "DileptonAnalysis"
  ! and initializes all the histograms.
  !****************************************************************************
  subroutine readInput
    use IdTable, only: nucleon, F37_1950
    use output, only: Write_ReadingInput
    use histMC, only: CreateHistMC, CopyDesc
    use particleProperties, only: hadron
    use nucleus, only: getTarget

    integer, parameter :: nbins = 200
    integer :: nCh, ios, i
    real :: d
    character(40) :: mbstr(5)
    character(80) :: title

    !**************************************************************************
    !****n* Dilepton_Analysis/DileptonAnalysis
    ! NAME
    ! NAMELIST /DileptonAnalysis/
    ! PURPOSE
    ! Includes the switches:
    ! * Enable
    ! * Extra
    ! * DeltaDalitz
    ! * DeltaDalitzFF
    ! * omegaDalitzFF
    ! * b_pi
    ! * lambda_eta
    ! * etaPrimeDalitzFF
    ! * angDist
    ! * brems
    ! * nEvent
    ! * nEvent_BH
    ! * kp_cut
    ! * kp_min
    ! * binsz
    ! * filter
    ! * hadesFilterFile
    ! * WriteEvents
    ! * p_lep_min
    ! * beta_gamma_cut
    ! * massBinning
    ! * particle_source
    !**************************************************************************
    NAMELIST /DileptonAnalysis/ Enable, Extra, DeltaDalitz, DeltaDalitzFF, omegaDalitzFF, b_pi, lambda_eta, etaPrimeDalitzFF, &
                                angDist, brems, nEvent, nEvent_BH, kp_cut, kp_min, binsz, filter, hadesFilterFile, WriteEvents, &
                                p_lep_min, beta_gamma_cut, massBinning, particle_source

    call Write_ReadingInput('DileptonAnalysis',0)
    rewind(5)
    read(5,nml=DileptonAnalysis,IOSTAT=IOS)
    call Write_ReadingInput('DileptonAnalysis',0,IOS)

    write(*,*) 'Enable dilepton analysis ? ', Enable
    if (Enable) then
      write(*,*) 'Extra histograms         ? ', Extra
      write(*,*) 'DeltaDalitz              ? ', DeltaDalitz
      write(*,*) 'DeltaDalitzFF            ? ', DeltaDalitzFF
      write(*,*) 'omegaDalitzFF            ? ', omegaDalitzFF
      write(*,*) 'b_pi                     ? ', b_pi
      write(*,*) 'lambda_eta               ? ', lambda_eta
      write(*,*) 'etaPrimeDalitzFF         ? ', etaPrimeDalitzFF
      write(*,*) 'angDist                  ? ', angDist
      write(*,*) 'bremsstrahlung           ? ', brems
      write(*,*) 'nEvent                   ? ', nEvent
      write(*,*) 'nEvent_BH                ? ', nEvent_BH
      write(*,*) 'kp_cut                   ? ', kp_cut
      if (kp_cut) then
        write(*,*) 'kp_min                   ? ', kp_min
      end if
      write(*,*) 'filter                   ? ', filter
      select case (filter)
      case (3,4)  ! HADES
        write(*,*) 'HADES filter file        ? ', trim(hadesFilterFile)
      case (5)    ! JLab
        write(*,*) 'min. lepton momentum     ? ', p_lep_min
      end select
      write(*,*) 'beta*gamma cut           ? ', beta_gamma_cut
      write(*,*) 'binsz                    ? ', binsz
      write(*,*) 'write events             ? ', WriteEvents
      write(*,'(A,4f6.3)') ' mass binning             ? ', massBinning(1:4)
      write(*,*) 'particle_source          ? ', particle_source
    end if

    call Write_ReadingInput('DileptonAnalysis',1)

    if (.not. Enable) return


    targNuc => getTarget()

    if (particle_source) then
      nCh = 65
    else
      nCh = 17
    end if

    ! adjust bin size according to available energy
    d=ProjectileEnergy_max/nbins
    ! create histograms
    call CreateHistMC (msigma,     'Dilepton Cross Section dSigma/dM [microbarn/GeV]',                0.,3.2,binsz,nCh)
    call CreateHistMC (ptsigma,    'Dilepton Cross Section dSigma/dP_t [microbarn/GeV]',              0.,d*nbins,d,nCh)
    call CreateHistMC (ysigma,     'Dilepton Cross Section dSigma/dy [microbarn]',                    -10.,10.,0.1,nCh)
    if (Extra) then
      call CreateHistMC (msigma_slow,'Dilepton Cross Section dSigma/dM [microbarn/GeV] (low momentum)', 0.,3.2,binsz,nCh)
      call CreateHistMC (msigma_fast,'Dilepton Cross Section dSigma/dM [microbarn/GeV] (high momentum)',0.,3.2,binsz,nCh)
      call CreateHistMC (msigma_pri, 'Dilepton Cross Section dSigma/dM [microbarn/GeV] (primary)',      0.,3.2,binsz,nCh)
      call CreateHistMC (msigma_sec, 'Dilepton Cross Section dSigma/dM [microbarn/GeV] (secondary)',    0.,3.2,binsz,nCh)
      write(title,'(A,f4.2,A)')  "Dilepton Cross Section dSigma/dM [microbarn/GeV] (beta*gamma < ",beta_gamma_cut," )"
      call CreateHistMC (msigma_bgcut, title, 0.,3.2,binsz,nCh)
      call CreateHistMC (psigma,     'Dilepton Cross Section dSigma/dP_lab [microbarn/GeV]',            0.,d*nbins,d,nCh)
      call CreateHistMC (ptsigma_pri,'Dilepton Cross Section dSigma/dP_t [microbarn/GeV] (primary)',    0.,d*nbins,d,nCh)
      call CreateHistMC (ptsigma_sec,'Dilepton Cross Section dSigma/dP_t [microbarn/GeV] (secondary)',  0.,d*nbins,d,nCh)
      call CreateHistMC (mtsigma,    'Dilepton Cross Section dSigma/dM_t [microbarn/GeV]',              0.,d*nbins,d,nCh)
      call CreateHistMC (tcmsigma,   'Dilepton Cross Section dSigma/dTheta_cm [microbarn]',             0.,180.,1.,nCh)
      call CreateHistMC (helsigma,   'Dilepton Cross Section dSigma/dalpha [microbarn]',                0.,180.,1.,nCh)
      call CreateHistMC (esigma, 'Dilepton Cross Section dSigma/dE_kin [microbarn/GeV]',   0.,d*nbins,d,nCh)
      call CreateHistMC (bgsigma,'Dilepton Cross Section dSigma/d(beta*gamma) [microbarn]',0.,50.,0.1,nCh)
      call CreateHistMC (oasigma,'Dilepton Cross Section dSigma/dTheta_open [microbarn]',  0.,180.,1.,nCh)
      call CreateHistMC (dsigma, 'Dilepton Cross Section dSigma/dRho [microbarn]',         0.,0.2,0.002,nCh)
    end if
    call CreateHistMC (deltaMass,'mass distribution of the Delta',0.9,3.,binsz,2)
    write(mbstr(1),'(A,f5.3,A)')        " (mass bin #1: 0 - ",massBinning(1)," GeV)"
    write(mbstr(2),'(A,f5.3,A,f5.3,A)') " (mass bin #2: ",massBinning(1)," - ",massBinning(2)," GeV)"
    write(mbstr(3),'(A,f5.3,A,f5.3,A)') " (mass bin #3: ",massBinning(2)," - ",massBinning(3)," GeV)"
    write(mbstr(4),'(A,f5.3,A,f5.3,A)') " (mass bin #4: ",massBinning(3)," - ",massBinning(4)," GeV)"
    write(mbstr(5),'(A,f5.3,A)')        " (mass bin #5: ",massBinning(4)," - inf GeV)"
    do i=1,5
      ptsigma(i)%name  = trim(ptsigma(i)%name)  // mbstr(i)
      ysigma(i)%name   = trim(ysigma(i)%name)   // mbstr(i)
      if (Extra) then
        psigma(i)%name   = trim(psigma(i)%name)   // mbstr(i)
        mtsigma(i)%name  = trim(mtsigma(i)%name)  // mbstr(i)
        tcmsigma(i)%name = trim(tcmsigma(i)%name) // mbstr(i)
        helsigma(i)%name = trim(helsigma(i)%name) // mbstr(i)
      end if
    end do
    ! set channel descriptions
    msigma%xDesc='Dilepton Mass [GeV]'
    msigma%yDesc(1:17) = (/ "rho -> e+e-        ", "omega -> e+e-      ", "phi -> e+e-        ", "omega -> pi0 e+e-  " , &
                            "pi0 -> e+e- gamma  ", "eta -> e+e- gamma  ", "Delta -> N e+e-    ", "eta -> e+e-        " , &
                            "eta' -> e+e- gamma ", "N*(1520) -> N e+e- ", "pn -> pn e+e-      ", "pp -> pp e+e-      " , &
                            "pi- n -> pi- n e+e-", "pi- p -> pi- p e+e-", "pi+ n -> pi+ n e+e-", "pi+ p -> pi+ p e+e-" , &
                            "Bethe-Heitler      " /)
    if (particle_source) then
      msigma%yDesc(38:41) = (/ "rho (from BB coll.)  ", "rho (from mm coll.)  ", "rho (from mB coll.)  ", "rho (via meson decay)" /)
      do i=Nucleon,F37_1950
        if (bar_ch(i)==0) cycle
        msigma%yDesc(17+bar_ch(i)) = trim(hadron(i)%name) // " Dalitz (VMD)"
        msigma%yDesc(41+bar_ch(i)) = "rho (via " // trim(hadron(i)%name) // ")"
      end do
      msigma%yDesc(62:65) = (/ "phi (from pi-rho coll.)", "phi (from K-Kbar coll.)", &
                               "phi (from BB coll.)    ", "phi (from mB coll.)    " /)
      call CreateHistMC (rhoMass,'mass distribution of the rho meson',0.,2.,binsz,24)
      rhoMass%xDesc='rho mass [GeV]'
      rhoMass%yDesc(1:24) = msigma%yDesc(38:61)
    end if
    call CopyDesc (ptsigma(:), msigma)
    call CopyDesc (ysigma(:), msigma)
    ptsigma(:)%xDesc='Transverse Momentum [GeV]'
    ysigma(:)%xDesc='Rapidity y'
    if (Extra) then
      call CopyDesc (msigma_slow, msigma)
      call CopyDesc (msigma_fast, msigma)
      call CopyDesc (msigma_bgcut, msigma)
      call CopyDesc (msigma_pri, msigma)
      call CopyDesc (msigma_sec, msigma)
      call CopyDesc (ptsigma_pri, msigma)
      call CopyDesc (ptsigma_sec, msigma)
      call CopyDesc (psigma(:), msigma)
      call CopyDesc (mtsigma(:), msigma)
      call CopyDesc (tcmsigma(:), msigma)
      call CopyDesc (helsigma(:), msigma)
      call CopyDesc (esigma, msigma)
      call CopyDesc (bgsigma, msigma)
      call CopyDesc (oasigma, msigma)
      call CopyDesc (dsigma, msigma)
      psigma(:)%xDesc='Momentum [GeV]'
      ptsigma_pri%xDesc='Transverse Momentum [GeV]'
      ptsigma_sec%xDesc='Transverse Momentum [GeV]'
      mtsigma(:)%xDesc='Transverse Mass [GeV]'
      tcmsigma(:)%xDesc='Theta_cm [degrees]'
      helsigma(:)%xDesc='(reconstructed) helicity alpha [degrees]'
      esigma%xDesc='Kinetic Energy [GeV]'
      bgsigma%xDesc='beta*gamma = p/m'
      oasigma%xDesc='Opening Angle [degrees]'
      dsigma%xDesc='Density [fm^-3]'
    end if
    deltaMass%xDesc='Delta Mass [GeV]'
    deltaMass%yDesc(1:2) = (/ 'vacuum   ', 'in-medium' /)

    init = .false.
  end subroutine readInput


  subroutine boost2Lab (mom)
    ! boost from calc frame to lab (if necessary):
    ! for some systems, we don't directly work in the lab frame, therefore we have to transform the 4-momenta
    use inputGeneral, only: eventtype
    use eventtypes, only: HiPion, HeavyIon
    use lorentzTrafo, only: lorentz
    real, dimension(0:3), intent(inout) :: mom
    if (eventtype==HeavyIon .and. boostToLab) then
      ! for heavy ion runs, we might have to boost to the Lab frame first
      call lorentz (betaCM2Lab, mom, "boost2Lab(1)")
    else if (eventtype==HiPion .and. targNuc%mass==2 .and. abs(ProjectileEnergy-1.25)<1E-3) then
      ! for p+d@1.25 (HADES) we have to do some transformations, since we actually want d+p (deuteron projectile shot on proton target!)
      call lorentz (betaTarg2Proj, mom, "boost2Lab(2)")  ! boost to the projectile frame
      mom(3) = -mom(3)                                   ! invert the z-axis
    end if
  end subroutine


  !****************************************************************************
  !****s* Dilepton_Analysis/Dilep_Decays
  ! NAME
  ! subroutine Dilep_Decays (t, part, last)
  ! PURPOSE
  ! Subroutine for calculation of dilepton spectra
  ! (all cross sections are calculated in microbarn/GeV).
  ! INPUTS
  ! * real, intent(in) :: t        - time [fm]
  ! * type(particle) :: part(:,:)  - particle vector (real or perturbative)
  ! * integer,intent(in) :: last   - are we in the last timestep yet?
  !   0=during time evolution, before last timestep,
  !   1=last timestep, before/during forced decays,
  !   2=last timestep, after forced decays
  ! OUTPUT
  ! * the dilepton spectra are stored in the histograms msigma,psigma,esigma etc
  ! * at then end of the run, these histograms are written to files
  !   'Dilepton*.dat' etc (e.g. 'DileptonMass.dat'), cf. Dilep_Write_CS
  !****************************************************************************
  subroutine Dilep_Decays (t, part, last)
    use particleDefinition, only: particle
    use inputGeneral, only: delta_T,numEnsembles,num_Runs_SameEnergy,num_energies,eventtype
    use constants, only: pi,melec,alphaQED,hbarc,mN,mPi
    use IdTable, only: Delta, D13_1520, F37_1950, pion, rho, phi, omegaMeson, eta, etaPrime, Kaon, Kaonbar, photon, &
                       EOV, NOP, invalidID, isBaryon, isMeson
    use eventtypes, only: LoPion,RealPhoton,HiPion,HiLepton,HeavyIon,HadronInduced=>Hadron
    use particleProperties, only: hadron
    use random, only: rnFlat
    use mesonWidthVacuum, only: dileptonWidth
    use mesonWidthMedium, only: WidthMesonMedium, get_MediumSwitchMesons
    use baryonWidthMedium, only: WidthBaryonMedium
    use dichteDefinition, only: dichte
    use densityModule, only: densityAt
    use mediumDefinition, only: medium
    use mediumModule, only: mediumAt
    use histMC, only: AddHistMC
    use minkowski, only: abs4
    use lorentzTrafo, only: lorentz
    use history, only: history_getGeneration, history_getParents
    use master_1Body, only: decayParticle_rhoN
    ! input
    integer, intent(in) :: last
    real, intent(in) :: t
    type(particle), intent(in) :: part(:,:)
    ! constants
    real, parameter :: gtotpi0 = 7.836E-9          ! total pi0 width
    real, parameter :: gampi0 = gtotpi0 * 0.98823  ! partial width of pi0 -> gamma gamma (for pi0 Dalitz decay)
    real, parameter :: gameta = 5.11E-7            ! partial width of eta -> gamma gamma (for eta Dalitz decay)
    real, parameter :: etadirect_BR = 5.6E-6       ! eta->e+e- direct decay branching ratio (upper limit by HADES, PDG 2012)
    ! local variables
    integer :: i, j, k, id, ch, gen, parents(1:3), channel
    real :: gamma, mass, perw, mass2, gam, dgamdm, gtot, massma, massmi, weight, rmass, massvac
    real :: mom(0:3), momLRF(0:3), pout(0:3,1:2), betaLRF(1:3), rmom(0:3)
    type(dichte) :: density
    type(medium) :: mediumAtPosition
    real, save :: pw
    logical, save :: first=.true.
    type(particle) :: fs(1:2)

    if (.not. Enable) return

    ! initialize perturbative weight factors and histograms
    if (first) then
      pw = 1. / float(num_Runs_SameEnergy*num_energies*nevent)
      select case (eventtype)
      case (RealPhoton)       ! RealPhoton events are initialized in microbarn/A
        pw = pw * targNuc%mass
      case (HiPion,HiLepton)  ! these are initialized in millibarn
        pw = pw * 1000.
      case (HeavyIon,HadronInduced)
        pw = pw / float(numEnsembles)
      case (LoPion)
        ! this part is not guaranteed to be correct (has never really been used)
        pw = pw * 10./float(numEnsembles*(targNuc%mass))
      case default
        write(*,*) "Error: Dilepton analysis not available for this event type!", eventtype
        stop
      end select
      first=.false.
    end if

    call CountParticles(t,part)

    !!! loop for calculating dilepton cross sections
    do k=1,nevent
    do i=lbound(part,1),ubound(part,1) ! loop over ensembles
    do j=lbound(part,2),ubound(part,2) ! loop over all perturbative particles
      id = part(i,j)%ID
      ch = part(i,j)%charge
      if (id==EOV) exit ! end of current ensemble reached, jump to next one
      if (invalidID(id)) then
        write(*,*) "found bad ID: ",id    ! error: should not occur
        stop
      end if
      if (id==NOP .or. id==photon) cycle
      if (last==0 .and. part(i,j)%in_Formation) cycle
      if (id<=F37_1950) then
        if (part(i,j)%antiParticle) then
          if (ch/=0 .and. ch/=-1) cycle
        else
          if (ch/=0 .and. ch/=1) cycle
        end if
        if (bar_ch(id)==0) cycle
      else
        if (ch/=0) cycle
      end if

      mom = part(i,j)%momentum
      gen = history_getGeneration (part(i,j)%history)

      ! boost from calc frame to lab
      call boost2Lab (mom)

      mass  = abs4(mom)
      gamma = mom(0)/mass
      massvac = part(i,j)%mass
      if (part(i,j)%perturbative) then
        perw = part(i,j)%perWeight * pw
      else
        perw = pw
      end if
      density = densityAt (part(i,j)%position)
      mediumAtPosition = mediumAt (density, part(i,j)%position)
      if (mass<2.*melec) then
        write(*,*) "Warning: Bad mass in Dilep!",id,mass,part(i,j)%mass
        cycle
      end if
      if (last>0 .and. (id==rho .or. id==omegaMeson .or. id==phi) .and. mass<hadron(id)%minmass) then
        if (k==1) then
          write(*,*) "WARNING in Dilep: subthreshold-Particle in last timestep!"
          gtot=WidthMesonMedium(id,mass,mom,mediumAtPosition)
          write(*,'(4G12.5)') id,mass,sqrt(dot_product(mom(1:3),mom(1:3))),part(i,j)%offShellParameter
          write(*,'(3G12.5)') mediumAtPosition%density,mediumAtPosition%useMedium,gtot
        end if
        if (.not. (get_MediumSwitchMesons() .and. mediumAtPosition%useMedium)) cycle
      end if
      ! select particle decay
      select case (id)
      !========================================================================
      case (rho)
        ! rho0 -> e+e-
        gam = dileptonWidth (ID, mass)
        if (last<1) then
          weight = perw*gam/gamma*delta_T/hbarc
        else
          gtot = WidthMesonMedium(id,mass,mom,mediumAtPosition)
          weight = perw*gam/gtot
        end if
        pout = DilepSim2 (mom, mass)
        call CS(pout,1,weight,mediumAtPosition%density,gen,part,i,j)
        if (particle_source) then
          ! determine where it came from (e.g. N*/Delta* decay or pi-pi collision)
          parents = history_getParents (part(i,j)%history)
          if (parents(2) > 0) then
             ! 2-body collisions
             if (isBaryon(parents(1)) .and. isBaryon(parents(2))) then
               channel = -3  ! channel 17: BB
             else if (isMeson(parents(1)) .and. isMeson(parents(2))) then
               channel = -2  ! channel 18: mm
             else
               channel = -1  ! channel 19: mB
             end if
             !print *,"dil, rho-parents (2-body):", parents(1:3)
          else if (isMeson(parents(1))) then
             channel = 0   ! channel 20: meson decay
             !print *,"dil, rho-parents (meson dec.):", parents(1:3)
          else if (isBaryon(parents(1))) then
             channel = bar_ch(parents(1))   ! channel 21-40
             if (channel == 0) then
               write(*,*) "wrong baryon in dilepton analysis: ", parents(1)
               stop
             end if
          else
             write(*,*) "error in dilepton analysis: wrong rho parents!",parents(1:2)
             stop
          end if
          call CS(pout,41+channel,weight,mediumAtPosition%density,gen,part,i,j)
          call AddHistMC (rhoMass, mass, 4+channel, perw)
        end if
        inclXS(1) = inclXS(1) + perw
      !========================================================================
      case (omegaMeson)
        ! omega -> e+e-
        gam = dileptonWidth (ID, mass)
        if (last<1) then
          weight = perw*gam/gamma*delta_T/hbarc
        else
          gtot = WidthMesonMedium(id,mass,mom,mediumAtPosition)+gam
          weight = perw*gam/gtot
        end if
        pout = DilepSim2 (mom, mass)
        call CS(pout,2,weight,mediumAtPosition%density,gen,part,i,j)
        ! omega -> pi0 e+e-
        massma = mass - mPi
        massmi = 2.*melec
        mass2  = rnFlat(massmi,massma)
        dgamdm = dGamma_dM_omegaDalitz(mass,mass2) * (massma-massmi)
        if (last<1) then
          weight = perw * dgamdm/gamma * delta_T/hbarc
        else
          gtot = WidthMesonMedium(id,mass,mom,mediumAtPosition)+gam
          weight = perw * dgamdm/gtot
        end if
        pout = DilepSim3 (mom, mass2, mPi, 0)
        call CS(pout,4,weight,mediumAtPosition%density,gen,part,i,j)
        inclXS(2) = inclXS(2) + perw
      !========================================================================
      case (phi)
        ! phi -> e+e-
        gam = dileptonWidth (ID, mass)
        if (last<1) then
          weight = perw*gam/gamma*delta_T/hbarc
        else
          gtot = WidthMesonMedium(id,mass,mom,mediumAtPosition)+gam
          weight = perw*gam/gtot
        end if
        pout = DilepSim2 (mom, mass)
        call CS(pout,3,weight,mediumAtPosition%density,gen,part,i,j)
        if (particle_source) then
          ! determine where it came from (BB, mm or mB collisions)
          parents = history_getParents (part(i,j)%history)
          ! print *,"dil, phi-parents (2-body):", parents(1:3)
          if (parents(2) > 0) then
            ! 2-body collisions
            if (isMeson(parents(1)) .and. isMeson(parents(2))) then
              if (parents(1)==pion .and. parents(2)==rho) then
                channel = 0  ! channel 41: pi-rho
              else if (parents(1)==Kaon .and. parents(2)==Kaonbar) then
                channel = 1  ! channel 42: KKbar
              else
                write(*,*) "error in dilepton analysis: wrong phi parents!",parents(1:2)
                stop
              end if
            else if (isBaryon(parents(1)) .and. isBaryon(parents(2))) then
              channel = 2  ! channel 43: BB
            else
              channel = 3  ! channel 44: mB
            end if
          else
            write(*,*) "error in dilepton analysis: wrong phi parents!",parents(1:2)
            stop
          end if
          call CS(pout,62+channel,weight,mediumAtPosition%density,gen,part,i,j)
        end if
        inclXS(3) = inclXS(3) + perw
      !========================================================================
      case (pion)
        ! pi0 -> e+e- gamma
        if (last==1) cycle
        if (angDist==2) then
          pout = PseudoscalarDalitzPythia (id, mom)
          dgamdm = gtotpi0 * 1.174E-2               ! integrated width of pi0 -> e+e- gamma [GeV] = Gamma_tot * BR
        else
          massma = mPi
          massmi = 2.*melec
          mass2  = rnFlat(massmi,massma)
          pout = DilepSim3 (mom, mass2, 0., angDist)
          dgamdm = 4*alphaQED/3./pi*gampi0/mass2*(1.-mass2**2/massma**2)**3*(1.+b_pi*mass2**2)**2 * (massma-massmi)
        end if
        if (last == 0) then
          weight = perw*dgamdm/gamma*delta_T/hbarc
        else if (last == 2) then
          weight = perw * dgamdm/gtotpi0
        end if
        call CS(pout,5,weight,mediumAtPosition%density,gen,part,i,j)
        if (last>1) inclXS(4) = inclXS(4) + perw
      !========================================================================
      case (eta)
        ! eta -> e+e-
        gam = hadron(id)%width * etadirect_BR
        if (last < 1) then
          weight = perw*gam/gamma*delta_T/hbarc
        else
          weight = perw * etadirect_BR
        end if
        pout = DilepSim2 (mom, mass)
        call CS(pout,8,weight,mediumAtPosition%density,gen,part,i,j)
        ! eta -> e+e- gamma
        if (angDist==2) then
          pout = PseudoscalarDalitzPythia (id, mom)
          dgamdm = hadron(eta)%width * 6.9E-3         ! integrated width of eta -> e+e- gamma [GeV] = Gamma_tot * BR
        else
          massma = hadron(id)%mass
          massmi = 2.*melec
          mass2  = rnFlat(massmi,massma)
          pout = DilepSim3 (mom, mass2, 0., angDist)
          dgamdm = 4*alphaQED/3./pi*gameta/mass2*(1.-mass2**2/massma**2)**3 * FF_eta(mass2) * (massma-massmi)
        end if
        if (last < 1) then
          weight = perw*dgamdm/gamma*delta_T/hbarc
        else
          weight = perw*dgamdm/hadron(id)%width
        end if
        call CS(pout,6,weight,mediumAtPosition%density,gen,part,i,j)
        inclXS(5) = inclXS(5) + perw
      !========================================================================
      case (etaPrime)
        ! eta' -> e+e- gamma
        massma = hadron(id)%mass
        massmi = 2.*melec
        mass2  = rnFlat(massmi,massma)
        dgamdm = dGamma_dM_etaPrimeDalitz (massma, mass2) * (massma-massmi)
        if (last < 1) then
          weight = perw*dgamdm/gamma*delta_T/hbarc
        else
          weight = perw*dgamdm/hadron(id)%width
        end if
        pout = DilepSim3 (mom, mass2, 0., angDist)
        call CS(pout,9,weight,mediumAtPosition%density,gen,part,i,j)
        inclXS(6) = inclXS(6) + perw
      !========================================================================
      case (Delta:F37_1950)  ! baryon resonances
        if (id==Delta) then
          ! Delta -> e+e- N (radiative)
          massma = massvac - mN              ! maximum dilepton mass (vacuum mass difference of Delta and Nucleon)
          massmi = 2.*melec                  ! minimum dilepton mass
          mass2  = rnFlat(massmi,massma)     ! actual dilepton mass
          dgamdm = dGamma_dM_DeltaDalitz(massvac,mass2,ch) * (massma-massmi)
          if (last < 1) then
            weight = perw * dgamdm/gamma * delta_T/hbarc
          else
            momLRF = mom
            if (density%baryon(0)>1E-8) then
              betaLRF(1:3) = density%baryon(1:3)/density%baryon(0)
              call lorentz (betaLRF, momLRF, 'Dilep_Decays')
            end if
            gtot = max (WidthBaryonMedium (id, massvac, momLRF, mediumAtPosition), 0.001)
            weight = perw * dgamdm/gtot
          end if
          call AddHistMC (deltaMass, massvac, 1, weight)
          call AddHistMC (deltaMass,    mass, 2, weight)
          pout = DilepSim3 (mom, mass2, mN, 0)
          call CS(pout,7,weight,mediumAtPosition%density,gen,part,i,j)
          inclXS(7+ch) = inclXS(7+ch) + perw

          if (hadron(Delta)%decaysID(2)==0) cycle  ! no rho-N decay
        else if (id == D13_1520) then
          ! N*(1520) -> e+e- N
          massma = massvac - mN              ! maximum dilepton mass (vacuum mass difference of N* and Nucleon)
          massmi = 2.*melec                  ! minimum dilepton mass
          mass2  = rnFlat(massmi,massma)     ! actual dilepton mass
          dgamdm = dGamma_dM_N1520Dalitz(massvac,mass2) * (massma-massmi)
          if (last < 1) then
            weight = perw * dgamdm/gamma * delta_T/hbarc
          else
            momLRF = mom
            if (density%baryon(0)>1E-8) then
              betaLRF(1:3) = density%baryon(1:3)/density%baryon(0)
              call lorentz (betaLRF, momLRF, 'Dilep_Decays')
            end if
            gtot = max (WidthBaryonMedium (id, massvac, momLRF, mediumAtPosition), 0.001)
            weight = perw * dgamdm/gtot
          end if
          pout = DilepSim3 (mom, mass2, mN, 0)
          call CS(pout,10,weight,mediumAtPosition%density,gen,part,i,j)
        end if
        !======================================================================
        if (particle_source) then
          ! R -> N e+e- (Dalitz decays via rho)
          ! (1) do R -> N rho0
          if (massvac < mN + hadron(rho)%minmass) cycle
          if (id==Delta .and. k>1) cycle
          gam = decayParticle_rhoN (part(i,j), fs(1:2))
          if (all(fs(1:2)%id==0)) then
            write(*,*) "warning: R->rhoN->e+e-N failed!"
            cycle
          else if (fs(1)%id/=rho .or. fs(1)%charge/=0) then
            write(*,*) "error in Dilep_Decays: no rho! ", id, fs(1:2)%ID, ch, fs(1:2)%charge
            stop
          end if
          ! (2) do rho0 -> e+e-
          rmom = fs(1)%momentum
          rmass = abs4(rmom)
          gam = gam * dileptonWidth (fs(1)%ID,rmass) / WidthMesonMedium(fs(1)%ID,rmass,rmom,mediumAtPosition)
          pout = DilepSim2 (rmom, rmass)
          ! (3) put it together
          if (last<1) then
            weight = perw * gam/gamma * delta_T/hbarc
          else
            momLRF = mom
            if (density%baryon(0)>1E-8) then
              betaLRF(1:3) = density%baryon(1:3)/density%baryon(0)
              call lorentz (betaLRF, momLRF, 'Dilep_Decays')
            end if
            gtot = max (WidthBaryonMedium (id, massvac, momLRF, mediumAtPosition), 0.001)
            weight = perw * gam/gtot
          end if
          if (id==Delta) weight = weight * nEvent
          call CS(pout,17+bar_ch(id),weight,mediumAtPosition%density,gen,part,i,j)
        end if
      end select
    end do
    end do
    end do ! k=1,nevent

  end subroutine Dilep_Decays


  !****************************************************************************
  !****f* Dilepton_Analysis/dGamma_dM_omegaDalitz
  ! NAME
  ! real function dGamma_dM_omegaDalitz(m_D,M)
  ! PURPOSE
  ! This function calculates the mass differential decay width dGamma/dM of
  ! omega -> pi0 e+e-, neglecting the electron mass.
  ! See:
  ! * Effenberger, PhD thesis, eq. (2.143) and (2.144), http://inspirehep.net/record/1375881
  ! * Terschluesen/Leupold, Phys. Lett. B691 (2010) 191-201, http://inspirehep.net/record/847729
  ! INPUTS
  ! * real,intent(in)    :: m_D     ! mass of the omega
  ! * real,intent(in)    :: M       ! invariant mass of the gamma*
  !****************************************************************************
  real function dGamma_dM_omegaDalitz(m_w,M)
    use constants, only: pi,alphaQED,mPi
    real,intent(in) :: m_w            ! mass of the omega
    real,intent(in) :: M              ! invariant mass of the gamma*
    real,parameter  :: gdom=7.03E-04  ! decay width "omega -> pi0 gamma"
    dGamma_dM_omegaDalitz = 2*alphaQED*gdom/(3.*pi*M) * (max((1+M**2/(m_w**2-mPi**2))**2 - 4.*m_w**2*M**2 /  &
                            (m_w**2-mPi**2)**2,0.))**(3./2.)
    select case (omegaDalitzFF)
    case (0)
      return  ! constant form factor, do nothing!
    case (1)
      dGamma_dM_omegaDalitz = dGamma_dM_omegaDalitz * FF_omega_effe(M)
    case (2)
      dGamma_dM_omegaDalitz = dGamma_dM_omegaDalitz * FF_VMD(M)
    case (3)
      dGamma_dM_omegaDalitz = dGamma_dM_omegaDalitz * FF_omega_terschluesen(m_w,M)
    case default
      write(*,*) 'Error in DileptonAnalysis: omegaDalitzFF =',omegaDalitzFF
      stop
    end select
  end function dGamma_dM_omegaDalitz

  ! generic VMD form factor
  real function FF_VMD (M)
    real,intent(in) :: M
    real,parameter :: mV = 0.776, GV = 0.150
    FF_VMD = mV**4/((mV**2-M**2)**2+mV**2*GV**2)
  end function FF_VMD

  ! eta Dalitz form factor
  real function FF_eta (M)
    real,intent(in) :: M
    FF_eta = 1. / ((1.-M**2/lambda_eta**2))**2
  end function FF_eta

  ! eta-prime Dalitz form factor as used by Genesis
  ! (VMD-like, in agreement with Lepton-G data)
  real function FF_etaprime_genesis (M)
    real,intent(in) :: M
    real,parameter :: mV = 0.764, GV = 0.102
    FF_etaprime_genesis = mV**4/((mV**2-M**2)**2+mV**2*GV**2)
  end function FF_etaprime_genesis

  ! omega Dalitz form factor used by Effenberger and Bratkovskaya
  real function FF_omega_effe (M)
    real,intent(in) :: M
    real,parameter :: lamo=0.650, gamo=0.075
    FF_omega_effe = lamo**4/((lamo**2-M**2)**2+lamo**2*gamo**2)
  end function FF_omega_effe

  ! omega Dalitz form factor by Terschluesen/Leupold, parameter set P1, including rho width.
  ! cf. Terschluesen/Leupold, Phys. Lett. B691 (2010) 191-201, http://inspirehep.net/record/847729 , eq (4),(16).
  real function FF_omega_terschluesen (m_w, M)
    use constants, only: ii, mPi
    use particleProperties, only: hadron
    use IDTable, only: rho
    use mesonWidth, only: FullWidthMeson
    real,intent(in) :: m_w            ! mass of the omega
    real,intent(in) :: M              ! invariant mass of the gamma*
    real,parameter  :: hA = 2.32, bA = 0.19
    complex :: f
    real :: f0  ! for nomalization
    f = (8.*bA*mPi**2/m_w - hA*m_w*(1+M**2/m_w**2)) / (M**2-hadron(rho)%mass**2+ii*M*FullWidthMeson(rho,M))  ! f(q^2)
    f0 = (8.*bA*mPi**2/m_w - hA*m_w) / hadron(rho)%mass**2                                                   ! f(0)
    FF_omega_terschluesen = abs((f/f0)**2)
  end function FF_omega_terschluesen


  ! eta-prime Dalitz form factor (standard VMD), as given by C. Terschluesen, Diplomarbeit, eq. (4.75)
  real function FF_etaprime_terschluesen (M)
    use constants, only: pi, ii
    use IDTable, only: rho, omegaMeson, phi
    use particleProperties, only: hadron
    use mesonWidth, only: FullWidthMeson
    real,intent(in) :: M                     ! invariant mass of the gamma*
    real, parameter :: th = -19.5 * pi/180.  ! eta/eta-prime mixing angle
    complex :: f
    real :: f0  ! for nomalization
    f = 9*(sin(th)+sqrt(2.)*cos(th)) * hadron(rho)%mass**2/(M**2-hadron(rho)%mass**2+ii*M*FullWidthMeson(rho,M)) &
      + (sin(th)+sqrt(2.)*cos(th)) * &
      & hadron(omegaMeson)%mass**2/(M**2-hadron(omegaMeson)%mass**2+ii*M*FullWidthMeson(omegaMeson,M)) &
      - 2*(2*sin(th)-sqrt(2.)*cos(th)) * hadron(phi)%mass**2/(M**2-hadron(phi)%mass**2+ii*M*FullWidthMeson(phi,M))
    f0 = 10*(sin(th)+sqrt(2.)*cos(th)) - 2*(2*sin(th)-sqrt(2.)*cos(th))
    FF_etaprime_terschluesen = abs((f/f0)**2)
  end function FF_etaprime_terschluesen


  !****************************************************************************
  !****f* Dilepton_Analysis/dGamma_dM_etaPrimeDalitz
  ! NAME
  ! real function dGamma_dM_etaPrimeDalitz (m_etap, M)
  ! PURPOSE
  ! This function calculates the mass differential decay width dGamma/dM of
  ! eta' -> e+e- gamma, neglecting the electron mass.
  ! See:
  ! * Effenberger, PhD thesis, eq. (2.140).
  ! INPUTS
  ! * real,intent(in)    :: m_etap  ! mass of the eta prime
  ! * real,intent(in)    :: M       ! invariant mass of the gamma*
  !****************************************************************************
  real function dGamma_dM_etaPrimeDalitz (m_etap, M)
    use constants, only: pi,alphaQED
    real,intent(in) :: m_etap            ! mass of the eta prima
    real,intent(in) :: M                 ! invariant mass of the gamma*
    real,parameter :: gametap = 4.3e-06  ! eta' dalitz decay (partial width for eta' -> gamma gamma)
    dGamma_dM_etaPrimeDalitz = 4*alphaQED*gametap/(3*pi*M) * (1.-M**2/m_etap**2)**3
    select case (etaPrimeDalitzFF)
    case (0)
      return  ! constant form factor, do nothing!
    case (1)
      dGamma_dM_etaPrimeDalitz = dGamma_dM_etaPrimeDalitz * FF_eta(M)
    case (2)
      dGamma_dM_etaPrimeDalitz = dGamma_dM_etaPrimeDalitz * FF_VMD(M)
    case (3)
      dGamma_dM_etaPrimeDalitz = dGamma_dM_etaPrimeDalitz * FF_etaprime_genesis (M)
    case (4)
      dGamma_dM_etaPrimeDalitz = dGamma_dM_etaPrimeDalitz * FF_etaprime_terschluesen (M)
    case default
      write(*,*) 'Error in DileptonAnalysis: etaPrimeDalitzFF =',etaPrimeDalitzFF
      stop
    end select
  end function dGamma_dM_etaPrimeDalitz


  !****************************************************************************
  !****f* Dilepton_Analysis/dGamma_dM_DeltaDalitz
  ! NAME
  ! real function dGamma_dM_DeltaDalitz (W, q, charge)
  ! PURPOSE
  ! This function calculates the mass differential decay width dGamma/dM of
  ! Delta -> N e+e-. The charge of the Delta must be +1 or 0.
  ! See e.g. Krivoruchenko/Faessler, Phys. Rev. D65 (2001), 017502,
  ! equation (3) and (4).
  ! INPUTS
  ! * real,intent(in)    :: W       ! mass of the Delta
  ! * real,intent(in)    :: q       ! invariant mass of the gamma*
  ! * integer,intent(in) :: charge  ! charge of the Delta
  !****************************************************************************
  real function dGamma_dM_DeltaDalitz (W, q, charge)
    use constants, only: pi, alphaQED, melec

    real, intent(in)    :: W, q
    integer, intent(in) :: charge

    dGamma_dM_DeltaDalitz = DeltaWidth_gammaN (W, q, charge) &
                            * 2.*alphaQED/(3.*pi*q) * (1.+2.*melec**2/q**2) * sqrt(1.-4.*melec**2/q**2)

  end function dGamma_dM_DeltaDalitz


  !****************************************************************************
  !****f* Dilepton_Analysis/DeltaWidth_gammaN
  ! NAME
  ! real function DeltaWidth_gammaN (W, q, charge)
  ! PURPOSE
  ! This function calculates the decay width of Delta -> N gamma*.
  ! The charge of the Delta must be +1 or 0.
  ! See e.g. Krivoruchenko/Faessler, Phys. Rev. D65 (2001), 017502, equation (2).
  ! INPUTS
  ! * real,intent(in)    :: W       ! mass of the Delta
  ! * real,intent(in)    :: q       ! invariant mass of the gamma*
  ! * integer,intent(in) :: charge  ! charge of the Delta
  !****************************************************************************
  real function DeltaWidth_gammaN (W, q, charge)
    use IDTable, only: Delta
    use baryonWidth, only: baryonWidth_gammaN

    real, intent(in) :: W, q
    integer, intent(in) :: charge

    select case (DeltaDalitz)
    case (1)
      DeltaWidth_gammaN = Gamma0_DeltaDalitz_Wolf (W, q)
    case (2)
      DeltaWidth_gammaN = Gamma0_DeltaDalitz_Krivo (W, q)
    case (3)
      DeltaWidth_gammaN = baryonWidth_gammaN (Delta, W, q, charge)
    case (4)
      DeltaWidth_gammaN = Gamma0_DeltaDalitz_Ernst (W, q)
    case default
      write(*,*) 'Error in DileptonAnalysis: DeltaDalitz = ', DeltaDalitz
      stop
    end select

  end function DeltaWidth_gammaN


  !****************************************************************************
  !****f* Dilepton_Analysis/Gamma0_DeltaDalitz_Wolf
  ! NAME
  ! real function Gamma0_DeltaDalitz_Wolf (W, q)
  ! PURPOSE
  ! This function calculates the total decay rate of a Delta resonance going
  ! into a nucleon and a gamma*.
  ! See Wolf et al., Nucl. Phys. A517 (1990) 615-638, http://inspirehep.net/record/306273 ,
  ! equation (4.9).
  ! Effenberger gives a slightly modified formula in his thesis, equation (2.139).
  ! REMARKS
  ! The coupling constant 'g' in Wolf's paper is wrong. Effenberger has the
  ! correct value of g = 5.44. With a real photon and an on-shell Delta one
  ! should get the photonic decay width of Gamma_0 (0) = 0.72 MeV.
  ! INPUTS
  ! * real, intent(in) :: W   ! mass of the Delta
  ! * real, intent(in) :: q   ! invariant mass of the gamma*
  ! OUTPUT
  ! Gamma_0 in GeV.
  !****************************************************************************
  real function Gamma0_DeltaDalitz_Wolf (W, q)
    use constants, only: pi, alphaQED, mN
    use twoBodyTools, only: pCM_sqr
    real, intent(in) :: W, q
    real, parameter  :: g2 = 5.44**2  ! coupling constant
    real :: pfinal2,mt,ml,lambda,f2,q0
    pfinal2 = pCM_sqr (W**2,mN**2,q**2)  ! dilepton momentum in Delta rest frame
    if (pfinal2<0.) then
      write(*,*) 'Warning: Problems in Gamma0_DeltaDalitz! ',pfinal2,W,q
      pfinal2=0.
    end if
    q0 = sqrt(q**2+pfinal2)   ! dilepton energy in Delta rest frame
    f2 = (1.5*(W+mN)/mN/((mN+W)**2-q**2))**2
    ml = alphaQED*4.*pi*f2*g2*W**2/(9.*mN) * 4.*q**2*(W-mN-q0)
    mt = alphaQED*4.*pi*f2*g2*W**2/(9.*mN) * (q0**2*(5.*W-3.*(q0+mN))-q**2*(W+mN+q0))
    lambda = max(0.,W**4 + q**4 + mN**4 - 2.*(W**2*q**2+W**2*mN**2+q**2*mN**2))
    Gamma0_DeltaDalitz_Wolf = sqrt(lambda) / (16.*pi*W**2) * mN * (2*mt+ml)
  end function Gamma0_DeltaDalitz_Wolf


  !****************************************************************************
  !****f* Dilepton_Analysis/Gamma0_DeltaDalitz_Ernst
  ! NAME
  ! real function Gamma0_DeltaDalitz_Ernst (W, q)
  ! PURPOSE
  ! This function calculates the total decay rate of a Delta resonance going
  ! into a nucleon and a gamma*, assuming a constant form factor.
  ! See Ernst et al., Phys. Rev. C 58 (1998) 447-456, http://inspirehep.net/record/452782 ,
  ! equations (3),(12),(13).
  ! INPUTS
  ! * real,intent(in) :: W        ! mass of the Delta
  ! * real,intent(in) :: q        ! invariant mass of the gamma* (dilepton)
  ! OUTPUT
  ! Gamma_0 in GeV.
  !****************************************************************************
  real function Gamma0_DeltaDalitz_Ernst (W, q)
    use constants, only: pi, mN, electronChargeSQ
    use twoBodyTools, only: pCM
    use NDeltaFF_Ramalho, only: NDeltaTL
    real, intent(in) :: W, q
    real :: FF

    select case (DeltaDalitzFF)
    case (1)  ! constant
      FF = 3.029**2
    case (2)  ! Dipole
      FF = Delta_FF_Dipole (q**2)
    case (3)  ! MAID
      FF = Delta_FF_MAID (W, q**2)
    case (4)  ! VMD
      FF = Delta_FF_VMD (q**2)
    case (5)  ! Wan/Iachello (two-component quark model)
      FF = Delta_FF_Iachello (W, q**2)
    case (6)  ! Ramalho/Pena
      FF = NDeltaTL (dble(q**2), dble(W))
    case default
      write(*,*) 'Error DileptonAnalysis: DeltaDalitzFF = ', DeltaDalitzFF
      stop
    end select

    Gamma0_DeltaDalitz_Ernst = pCM(W,mN,q) / (8.*pi*W**2) &
                               * electronChargeSQ * FF * (W+mN)**2 / (4.*mN**2*((W+mN)**2-q**2)**2) &
                               * ((W-mN)**2 - q**2) * (7.*W**4 + 14.*W**2*q**2 + 3.*q**4 + 8.*W**3*mN &
                                                       + 2.*W**2*mN**2 + 6.*q**2*mN**2 + 3.*mN**4)
  end function


  !****************************************************************************
  !****f* Dilepton_Analysis/Gamma0_DeltaDalitz_Krivo
  ! NAME
  ! real function Gamma0_DeltaDalitz_Krivo (W, q)
  ! PURPOSE
  ! This function calculates the total decay rate of a Delta resonance going
  ! into a nucleon and a gamma*.
  ! See Krivoruchenko/Faessler, Phys. Rev. D65 (2001) 017502, http://inspirehep.net/record/555421 ,
  ! equation (2).
  ! REMARKS
  ! With a real photon and an on-shell Delta one should get the PDG value of
  ! Gamma_0 (0) = 0.66 MeV.
  ! INPUTS
  ! * real, intent(in) :: W        ! mass of the Delta
  ! * real, intent(in) :: q        ! invariant mass of the gamma* (dilepton)
  ! OUTPUT
  ! Gamma_0 in GeV.
  !****************************************************************************
  real function Gamma0_DeltaDalitz_Krivo (W, q)
    use constants, only: alphaQED, mN
    use NDeltaFF_Ramalho, only: NDeltaTL
    real, intent(in) :: W, q
    real :: FF

    select case (DeltaDalitzFF)
    case (1)  ! constant
      FF = 3.029**2
    case (2)  ! Dipole
      FF = Delta_FF_Dipole (q**2)
    case (3)  ! MAID
      FF = Delta_FF_MAID (W, q**2)
    case (4)  ! VMD
      FF = Delta_FF_VMD (q**2)
    case (5)  ! Wan/Iachello (two-component quark model)
      FF = Delta_FF_Iachello (W, q**2)
    case (6)  ! Ramalho/Pena
      FF = NDeltaTL (dble(q**2), dble(W))
    case default
      write(*,*) 'Error DileptonAnalysis: DeltaDalitzFF = ', DeltaDalitzFF
      stop
    end select

    Gamma0_DeltaDalitz_Krivo = alphaQED/16. * (W+mN)**2/(W**3*mN**2) * ((W+mN)**2-q**2)**0.5 * ((W-mN)**2-q**2)**1.5 * FF

  end function Gamma0_DeltaDalitz_Krivo


  !****************************************************************************
  !****f* Dilepton_Analysis/Delta_FF_MAID
  ! NAME
  ! real function Delta_FF_MAID (W,q2)
  ! PURPOSE
  ! This function uses the electromagnetic N-Delta transition form factors from
  ! MAID2005, as used in GiBUU electroproduction but with off-shell W and
  ! continued to negative Q^2. References:
  ! * D. Drechsel, O. Hanstein, S.S. Kamalov, L. Tiator; Nucl. Phys. A 645 (1999) 145, eq. (21) - (23).
  ! * M. Warns, H. Schroeder, W. Pfeil, H. Rollnik; Z. Phys. C 45 (1990) 627, eq. (2.16) - (2.17)
  ! INPUTS
  ! * real,intent(in) :: W          ! mass of the Delta
  ! * real,intent(in) :: q2         ! q^2 = m_ee^2 = - Q^2
  ! OUTPUT
  ! Squared form factor.
  !****************************************************************************
  real function Delta_FF_MAID (W,q2)
    real,intent(in) :: W           ! off-shell mass of the Delta
    real,intent(in) :: q2             ! q^2 = m_ee^2

    Delta_FF_MAID = MAID_FF (W,-q2) / MAID_FF(W,0.) * 3.029**2

  contains

    real function MAID_FF (W,Qs)
      use constants, only: alphaQED,pi,mN
      use helicityAmplitudes, only: HP33_MAID05
      real, intent(in) :: W, Qs
      real :: A1, A3, S1, G_M, G_E, G_C, K_C, K_W, Z, kcm0, kcm, DE, DM, DS

      ! get MAID form factors
      kcm0 = (W**2-mN**2)/(2.*W)
      kcm = sqrt( ((W**2-Qs-mN**2)/(2.*W))**2 + Qs)
      CALL HP33_MAID05 (dble(Qs),dble(kcm),dble(kcm0),dble(DE),dble(DM),dble(DS),dble(A1),dble(A3),dble(S1),0)
      ! convert to GeV^(-1/2)
      A1 =   A1*1E-3
      A3 =   A3*1E-3
      S1 = - S1*1E-3   ! we use different sign in S_1/2 amplitude
      ! convert to Sachs form factors, cf. Warns et al.
      K_W = kcm0
      K_C = kcm
      Z   = sqrt(mN*K_W/(4.*pi*alphaQED*W)*(1.+Qs/(mN+W)**2)) * mN/K_C
      G_M = -Z * (A1+A3*sqrt(3.))
      G_E =  Z * (A1-A3/sqrt(3.))
      G_C =  Z * 2.*W/K_C*sqrt(2.)*S1

      MAID_FF = (G_M**2 + 3.*G_E**2 + mN**2/(2*W**2)*G_C**2)
    end function

  end function Delta_FF_MAID


  ! simple Dipole form factor with standard parameter
  real function Delta_FF_Dipole (q2)
    real, intent(in) :: q2
    real, parameter :: Lambda = 0.71
    Delta_FF_Dipole = (1.-q2/Lambda)**(-4) * 3.029**2
  end function


  ! simple VMD form factor as used e.g. by Titov,Kämpfer,Bratkovskaya; Phys. Rev. C 51 (1995) 227
  real function Delta_FF_VMD (q2)
    real, intent(in) :: q2           ! invariant mass of gamma* (squared)
    real, parameter :: Mrho = 0.761
    real, parameter :: Grho = 0.118
    real :: f, f0
    f = (Mrho**4 / ((q2-Mrho**2)**2 + (Mrho*Grho)**2))**2
    f0 = (Mrho**4 / (Mrho**4 + (Mrho*Grho)**2))**2
    Delta_FF_VMD = 3.029**2 * f/f0
  end function


  !****************************************************************************
  !****f* Dilepton_Analysis/Delta_FF_Iachello
  ! NAME
  ! real function Delta_FF_Iachello (W, q2)
  ! PURPOSE
  ! This function calculates the electromagnetic N-Delta transition form factor
  ! according to the two-component quark model by Wan/Iachello. References:
  ! * Q. Wan, F. Iachello: A unified description of baryon electromagnetic form factors.
  !   Int. J. Mod. Phys. A 20 (2005) 1846, http://inspirehep.net/record/689265
  ! * F. Iachello, Q. Wan: Structure of the nucleon from electromagnetic timelike form factors.
  !   Phys. Rev. C 69 (2004) 055204, http://inspirehep.net/record/651033
  ! * I. Froehlich, F. Dohrmann, T. Galatyuk et al.: A versatile method for simulating pp->ppe+e- and dp->pne+e-p_spec reactions.
  !   Eur. Phys. J. A45 (2010) 401-411, http://inspirehep.net/record/832525
  ! INPUTS
  ! * real, intent(in) :: W      ! off-shell mass of the Delta [GeV]
  ! * real, intent(in) :: q2     ! q^2 = m_ee^2 [GeV^2]
  ! OUTPUT
  ! Squared form factor.
  !****************************************************************************
  real function Delta_FF_Iachello (W, q2)
    use constants, only: pi, ii, mN, mPi
    real, intent(in) :: W, q2
    real, parameter :: mu_p = 2.793
    real, parameter :: a = 0.29
    real, parameter :: theta = 53.*pi/180.
    real, parameter :: b = 1.2147
    real, parameter :: bp = 0.004
    complex :: FF

    FF = mu_p * 4./(3.*sqrt(2.)) * sqrt(2.*mN*W/(W**2+mN**2)) / (1.-a*q2*exp(ii*theta))**2 * (bp+b*F_rho())

    Delta_FF_Iachello = abs(FF**2)

  contains

    complex function F_rho ()
      real, parameter :: m_rho = 0.765
      real, parameter :: G_rho = 0.112
      real :: x

      x = q2/(4.*mPi**2)
      F_rho = ( m_rho**2 + 8.*G_rho*mPi/pi ) / ( m_rho**2 - q2 + 4.*mPi*(1.-x)*G_rho*(alpha(x)-ii*gamma_(x)) )

    end function

    real function alpha (x)
      real, intent(in) :: x
      if (x<0.) then ! space-like part
        alpha = 2./pi * sqrt((x-1.)/x) * log(sqrt(1.-x)+sqrt(-x))
      else if (x>1.) then
        alpha = 2./pi * sqrt((x-1.)/x) * log(sqrt(x-1.)+sqrt(x))
      else
        alpha = sqrt((1.-x)/x) * (1.-2./pi*acot(sqrt(x/(1-x))))
      end if
    end function

    real function gamma_ (x)
      real, intent(in) :: x
      if (x>1.) then
        gamma_ = sqrt((x-1.)/x)
      else
        gamma_ = 0.
      end if
    end function

    real function acot (x)
      real, intent(in) :: x
      acot = pi/2. - atan(x)
    end function

  end function Delta_FF_Iachello


  !****************************************************************************
  !****f* Dilepton_Analysis/dGamma_dM_N1520Dalitz
  ! NAME
  ! real function dGamma_dM_N1520Dalitz (W, q)
  ! PURPOSE
  ! This function calculates the mass differential decay width dGamma/dM of
  ! N*(1520) -> N e+e-. See Krivoruchenko/Martemyanov/Faessler/Fuchs,
  ! Ann. Phys. 296 (2002), 299-346, equation (III.22).
  ! INPUTS
  ! * real,intent(in)    :: W       ! mass of the N*
  ! * real,intent(in)    :: q       ! invariant mass of the gamma*
  !****************************************************************************
  real function dGamma_dM_N1520Dalitz (W, q)
    use constants, only: pi, alphaQED, melec, mN

    real, intent(in)    :: W, q

    real :: width_gammaN
    ! squared form factor, normalized to photon point (which is just a constant in QED approximation)
    real, parameter :: FF = 1.182

    ! N* -> N gamma*
    width_gammaN = alphaQED/16. * (W-mN)**2/(W**3*mN**2) * ((W-mN)**2-q**2)**0.5 * ((W+mN)**2-q**2)**1.5 * FF

    ! Dalitz width with lepton phase-space factors
    dGamma_dM_N1520Dalitz = width_gammaN &
                            * 2.*alphaQED/(3.*pi*q) * (1.+2.*melec**2/q**2) * sqrt(1.-4.*melec**2/q**2)

  end function dGamma_dM_N1520Dalitz


  !****************************************************************************
  !****s* Dilepton_Analysis/CountParticles
  ! NAME
  ! subroutine CountParticles(nt,pertPart)
  ! PURPOSE
  ! This subroutine simply counts the number of particles at each timestep.
  ! Currently only rho^0, omega, phi, pi^0, eta and Delta^(0,+) are counted.
  ! All these particles are counted separately, and in addition the particles
  ! which are not in formation are counted (for each species).
  ! The numbers are per ensemble, and averaged over all the runs.
  ! INPUTS
  ! * real t - time [fm]
  ! * type(particle) pertPart(:,:) - list of perturbative particles
  ! OUTPUT
  ! The data is stored in two multi-channel histograms:
  ! "partnum" and "partnum_noform".
  !****************************************************************************
  subroutine CountParticles(t,pertPart)
    use inputGeneral, only: delta_T,numTimeSteps,num_Runs_SameEnergy,numEnsembles
    use particleDefinition, only: particle
    use IdTable, only: Delta,pion,eta,rho,phi,omegaMeson,P11_1440,F37_1950
    use histMC, only: CreateHistMC, AddHistMC, CopyDesc
    use minkowski, only: abs4, abs4sq
    use output, only: writeParticle
    ! inputs
    real,intent(in):: t
    type(particle),intent(in):: pertPart(:,:)
    !!!
    integer::i,j
    real:: mass
    logical,save::first=.true.
    real,save::dn1=0,dn2=0
    if (first) then
      dn1=1./float(num_Runs_SameEnergy*numEnsembles)
      dn2=delta_T/float(num_Runs_SameEnergy*numEnsembles)
      call CreateHistMC (partnum, 'Particle Numbers (Total)', -0.5*delta_T, (numTimeSteps+1.5)*delta_T, delta_T, 12)
      partnum%xDesc = "time [fm]"
      partnum%yDesc(1:12) = (/ "pi-    ", "pi0    ", "pi+    ", "Delta- ", "Delta0 ", "Delta+ ", "Delta++", &
                               "eta    ", "rho0   ", "omega  ", "phi    ", "Res    " /)
      call CreateHistMC(partnum_noform,'Particle Numbers (Not In Formation)',-0.5*delta_T,(numTimeSteps+1.5)*delta_T,delta_T,12)
      call CopyDesc(partnum_noform,partnum)
      call CreateHistMC(vm_mass,'VM Mass Spectra',0.,2.0,0.01,3)
      vm_mass%xDesc='Mass [GeV]'
      vm_mass%yDesc(1:3) = (/ "rho0 ", "omega", "phi  " /)
      first=.false.
    end if
    do i=lbound(pertPart,1),ubound(pertPart,1) ! loop over ensembles
    do j=lbound(pertPart,2),ubound(pertPart,2) ! loop over all perturbative particles
      mass=abs4(pertPart(i,j)%momentum)
      if (pertPart(i,j)%ID>=rho .and. pertPart(i,j)%ID<=phi .and. .not.(mass>0.)) then
        write(*,*) "Warning: Bad mass in CountParticles!",pertPart(i,j)%ID,pertPart(i,j)%charge,pertPart(i,j)%mass,mass, &
                   pertPart(i,j)%momentum,abs4sq(pertPart(i,j)%momentum)
        call writeParticle(6,i,j,pertPart(i,j))
        stop
      end if

      select case (pertPart(i,j)%ID)
      case (pion)
        call AddHistMC(partnum,t,2+pertPart(i,j)%charge,dn1)
      case (Delta)
        call AddHistMC(partnum,t,5+pertPart(i,j)%charge,dn1)
      case (eta)
        call AddHistMC(partnum,t,8,dn1)
      case (rho)
        if (pertPart(i,j)%charge==0) then
          call AddHistMC(partnum,t,9,dn1)
          call AddHistMC(vm_mass,mass,1,dn2*pertPart(i,j)%perweight)
        end if
      case (omegaMeson)
        call AddHistMC(partnum,t,10,dn1)
        call AddHistMC(vm_mass,mass,2,dn2*pertPart(i,j)%perweight)
      case (phi)
        call AddHistMC(partnum,t,11,dn1)
        call AddHistMC(vm_mass,mass,3,dn2*pertPart(i,j)%perweight)
      case (P11_1440:F37_1950)
        call AddHistMC(partnum,t,12,dn1)
      end select

      if (pertPart(i,j)%in_Formation) cycle

      select case (pertPart(i,j)%ID)
      case (pion)
        call AddHistMC(partnum_noform,t,2+pertPart(i,j)%charge,dn1)
      case (Delta)
        call AddHistMC(partnum_noform,t,5+pertPart(i,j)%charge,dn1)
      case (eta)
        call AddHistMC(partnum_noform,t,8,dn1)
      case (rho)
        if (pertPart(i,j)%charge==0) call AddHistMC(partnum_noform,t,9,dn1)
      case (omegaMeson)
        call AddHistMC(partnum_noform,t,10,dn1)
      case (phi)
        call AddHistMC(partnum_noform,t,11,dn1)
      case (P11_1440:F37_1950)
        call AddHistMC(partnum_noform,t,12,dn1)
      end select

    end do
    end do
  end subroutine CountParticles


  !****************************************************************************
  !****s* Dilepton_Analysis/writePartNum
  ! NAME
  ! subroutine writePartNum()
  ! PURPOSE
  ! Writes the data collected by "CountParticles" to disk, creating three files:
  ! * PartNum.dat: Total particle numbers per timestep, including those which are still in formation
  ! * PartNum_Integrated.dat: Time-integrated particle numbers.
  ! * PartNum_NoForm: Numbers of particles which are not in formation.
  !****************************************************************************
  subroutine writePartNum
    use histMC, only: WriteHistMC,WriteHistMC_Integrated
    ! counted numbers of particles
    call WriteHistMC(partnum,'PartNum.dat',.false.)
    ! integrated numbers of particles
    call WriteHistMC_Integrated(partnum,'PartNumInt.dat',.false.)
    ! counted numbers of particles which are not in formation time
    call WriteHistMC(partnum_noform,'PartNum_NoForm.dat',.false.)
    ! vector meson mass spectrum
    call WriteHistMC(vm_mass,'VMmass.dat',.false.)
  end subroutine writePartNum


  logical function isExclusive (parts, pnr, iso)
    use particleDefinition
    use IdTable, only: NOP,EOV,nucleon
    type(particle), intent(in) :: parts(:) ! particle vector
    integer, intent(in) :: pnr,iso
    integer :: i, id
    isExclusive = .false.
    select case (iso)
    case (4,5,6,9)  ! meson Dalitz decays are not exclusive
      return
    case (8)        ! we also suppress the eta direct here
      return
    end select
    do i=lbound(parts,1),ubound(parts,1) ! loop over all particles
      id = parts(i)%ID
      if (id==NOP) cycle
      if (id==EOV) exit
      if (i/=pnr) then
        if (id/=nucleon) return
      end if
    end do
    isExclusive = .true.
  end function isExclusive


  !****************************************************************************
  !****s* Dilepton_Analysis/WriteFullEvent
  ! NAME
  ! subroutine WriteFullEvent (parts, dil, pnr, pw, iso)
  ! PURPOSE
  ! Print out a full dilepton event, including all hadrons and the lepton pair.
  ! The event will be written to the file "Dilepton_FullEvents.dat" in a
  ! format similar to the LesHouches format.
  !****************************************************************************
  subroutine WriteFullEvent (parts, dil, pnr, pw, iso, f)
    use particleDefinition
    use IdTable, only: electron, photon, pion, nucleon, EOV
    use inputGeneral, only: eventType
    use eventtypes, only: hiLepton
    use EventInfo_HiLep, only: EventInfo_HiLep_Get
    use eN_event, only: eNev_SetProcess, eNev_init_enQ!, eNev_Set_PhiLepton
    use eN_eventDefinition

    type(particle), intent(in) :: parts(:) ! particle vector
    real, intent(in) :: dil(0:3,1:2)       ! e+(e-) 4-momentum
    integer, intent(in) :: pnr
    real, intent(in) :: pw
    integer, intent(in) :: iso, f

    real, dimension(0:3) :: ptot, pmiss
    integer :: i,id,iFE,evtType
    real :: w, nu, Q2, eps, phiL, thetaL
    integer, parameter :: iFile = 24
    character(len=16), parameter :: f3 = '(2I4,4ES13.5,I4)'
    type(electronNucleon_event) :: evt
    logical :: flagOK

    if (particle_source .and. (iso==1 .or. iso==3)) return                ! prevent double counting of rhos and phis
    if (writeEvents>=4 .and. .not. isExclusive (parts, pnr, iso)) return  ! print only exclusive events
    if (writeEvents==5 .and. iso<=20) return                              ! print only R->Ne+e-

    !**************************************************************************
    !****o* Dilepton_Analysis/Dilepton_FullEvents.dat
    ! NAME
    ! file Dilepton_FullEvents.dat
    ! PURPOSE
    ! This file contains dilepton events (including all hadrons). It will only
    ! be written if enabled by the WriteEvents switch.
    ! Attention: This file can get very large for long runs!
    ! The format is XML-like and similar to the LesHouches format: Each event is
    ! contained in an <event> .... </event> block, with a header in the <event>
    ! line that contains the perturbative weight, source type (channel) and
    ! filter result (1=accept, 0=reject). Each line in the event describes one
    ! particle and has the following entries:
    ! * Column #1: Particle ID.
    ! * Column #2: Charge.
    ! * Column #3-6: 4-momentum (E,p_x,p_y,p_z) in GeV.
    ! * Column #7: Parent particle ID (if the particle came from a decay, 0 otherwise).
    !**************************************************************************
    open(iFile, file="Dilepton_FullEvents.dat", position="append")
    write(iFile, '(A,ES13.5,2I3)') '<event>', pw, iso, f

    ptot = 0.
    do i=lbound(parts,1),ubound(parts,1) ! loop over all particles
      if (parts(i)%ID==EOV) exit

      if (i /= pnr) then
         ! print hadronic remainder of the event
         call printStableParts (parts(i))
      else
         ! print dilepton pair
         ptot = ptot + dil(:,1) + dil(:,2)
         write(iFile, f3) -electron, 1, dil(0:3,1), parts(i)%ID
         write(iFile, f3) electron, -1, dil(0:3,2), parts(i)%ID
         iFE = parts(i)%firstEvent
         if ((iso>=4 .and. iso<=7) .or. iso==9 .or. iso==10) then
            ! Dalitz decays: print missing particle
            pmiss = parts(i)%momentum - dil(:,1) - dil(:,2)
            ptot = ptot + pmiss
            select case (iso)
            case (4)      ! omega Dalitz
               id = pion
            case (5,6,9)  ! pi0, eta, eta' Dalitz
               id = photon
            case (7,10)    ! Delta, N*(1520) Dalitz
               if (parts(i)%antiParticle) then
                  id = -nucleon
               else
                  id = nucleon
               end if
            end select
            write(iFile, f3) id, parts(i)%charge, pmiss(0:3), parts(i)%ID
         end if
      end if
    end do

    if (eventType == hiLepton) then
      if (EventInfo_HiLep_Get (0, iFE, w, nu, Q2, eps, evtType, phi_Lepton=phiL)) then
        call eNev_SetProcess(evt, 1,1)  ! set to EM and electron
        call eNev_init_enQ(evt, eps, nu, Q2, flagOK)
        thetaL =  atan2(evt%lepton_in%momentum(1),evt%lepton_in%momentum(3))
!         call eNev_Set_PhiLepton(evt, phiL)  ! possibly rotate by angle phi
        if (flagOK) write(iFile, f3) electron, -1, evt%lepton_out%momentum(0:3), 0
        write(iFile,'(A,1P,5e13.4,0P,I8)') '# 14 ', nu, Q2, eps, phiL, thetaL, evtType
      end if
    end if

    write(iFile, '(A)') '</event>'
    close(iFile)

  contains

    recursive subroutine printStableParts (p)
      use IdTable, only: isHadron
      use particleProperties, only: hadron
      use master_1Body, only: decayParticle
      use history, only: history_getParents

      type(particle),intent(in) :: p

      type(particle),dimension(1:3) :: fs
      logical :: cflag, pflag, stable
      integer :: j, parents(3), par, sgn

      if (isHadron(p%id)) then
         stable = hadron(p%id)%stability==0
      else
         stable = .true.
      end if

      if (stable) then
         ! print out
         ptot = ptot + p%momentum
         parents = history_getParents (p%history)
         if (parents(2) == 0) then
            par = parents(1)
         else
            par = 0
         end if
         if (p%antiparticle) then
           sgn = -1
         else
           sgn = 1
         end if
         write(iFile,f3) sgn*p%id, p%charge, p%momentum(0:3), par
      else
         ! force decay of unstable particles
         call decayParticle (p, fs, cflag, pflag, .true., 0.)
         do j=1,3
            if (fs(j)%ID>0) call printStableParts (fs(j))
         end do
      end if
    end subroutine

  end subroutine WriteFullEvent


  subroutine WriteEvent (pout, pw, iso, f)
    real, dimension(0:3,1:2), intent(in)::pout ! e+(e-) 4-momentum
    real, intent(in) :: pw                     ! perturbative weight
    integer, intent(in) :: iso                 ! source type
    integer, intent(in) :: f
    if (particle_source .and. (iso==1 .or. iso==3)) return  ! prevent double counting of rhos and phis
    !**************************************************************************
    !****o* Dilepton_Analysis/Dilepton_Events.dat
    ! NAME
    ! file Dilepton_Events.dat
    ! PURPOSE
    ! This file contains dilepton pair events (leptons only, no hadrons).
    ! It will only be written if enabled by the WriteEvents switch.
    ! Attention: This file can get very large for long runs!
    ! Entries:
    ! * Column #1: Projectile energy in GeV.
    ! * Column #2: Lepton charge (+ or -).
    ! * Column #3-5: Lepton 3-momentum (p_x,p_y,p_z) in GeV.
    ! * Column #6: Perturbative weight.
    ! * Column #7: Source type (channel).
    ! * Column #8: Filter result (1=accept, 0=reject).
    !**************************************************************************
    open(999,file="Dilepton_Events.dat",position="append")
    write(999,'(ES13.5,A2,4ES13.5,2I3)') ProjectileEnergy,"+",pout(1:3,1),pw,iso,f ! positron
    write(999,'(ES13.5,A2,4ES13.5,2I3)') ProjectileEnergy,"-",pout(1:3,2),pw,iso,f ! electron
    close(999)
  end subroutine


  real function get_helicity (pout)
    use lorentzTrafo, only: lorentz
    use minkowski, only: op_ang

    real, dimension(0:3,1:2), intent(in) :: pout          ! e+(e-) 4-momentum

    real, dimension(0:3) :: pel, ppos, pdil   ! four-momenta of electron, positron and dilepton
    real, dimension(1:3) :: beta

    ppos = pout(:,1)
    pel  = pout(:,2)
    pdil = pel+ppos

    beta = pdil(1:3)/pdil(0)

    pdil(1:3) = -pdil(1:3)  ! take backward dil. direction

    ! boost to dil. cm frame
    call lorentz (beta, ppos, "get_helicity(1)")
    call lorentz (beta, pdil, "get_helicity(2)")

    get_helicity = op_ang (ppos, pdil)

  end function get_helicity


  !****************************************************************************
  !****s* Dilepton_Analysis/CS
  ! NAME
  ! subroutine CS (pout_, iso, pw, dens, gen, part, enr, pnr)
  ! PURPOSE
  ! Calculates differential cross sections (in microbarn) and puts them into
  ! various histograms.
  ! INPUTS
  ! * real, dimension(0:3,1:2), intent(in) :: pout --- 4-momenta of the lepton pair (1=positron, 2=electron)
  ! * integer, intent(in) :: iso --- source type (channel)
  ! * real, intent(in) :: pw --- perturbative weight
  ! * real, intent(in) :: dens --- density at decay vertex
  ! * integer, intent(in) :: gen --- generation
  ! * type(particle), intent(in), optional :: part(:,:) --- particle vector
  ! * integer, intent(in), optional :: enr --- ensemble nr
  ! * integer, intent(in), optional :: pnr --- particle nr
  ! OUTPUT
  !  * cross section data is stored in histograms msigma, psigma, esigma, ...
  !****************************************************************************
  subroutine CS (pout_, iso, pw, dens, gen, part, enr, pnr)
    use particleDefinition
    use constants, only: pi
    use histMC, only: AddHistMC
    use minkowski, only: SP, abs3, op_ang
    use random, only: rnGauss
    use rotation, only: get_phi_Theta
    use lorentzTrafo, only: lorentz
    ! input variables
    real, dimension(0:3,1:2), intent(in)::pout_
    integer, intent(in) :: iso
    real, intent(in) :: pw
    real, intent(in) :: dens
    integer, intent(in) :: gen
    type(particle), intent(in), optional :: part(:,:)
    integer, intent(in), optional :: enr, pnr
    ! local variables
    real :: mass,plab,ekin,y,pt,mt,beta_gamma,angle,theta_cm,phi,hel
    real, dimension(0:3) :: ptot
    real, dimension(0:3,1:2) :: pout
    integer :: f   ! filtering result: 1=accept, 0=reject
    integer :: mb  ! mass bin

    pout = pout_   ! we may need to modify this locally (filtering!)
    ! k*p cut
    if (kp_cut .and. (SP(ProjectileMomentum,pout(:,1))<kp_min .or. SP(ProjectileMomentum,pout(:,2))<kp_min)) return
    !!! === #1: apply filter ===
    f = applyFilter (pout(:,1),pout(:,2))
    !!! === #2: write out (all) events ===
    select case (WriteEvents)
    case (1)
      call WriteEvent (pout, pw, iso, f)
    case (2)
      if (present(part)) then
        if (isExclusive(part(enr,:),pnr,iso)) call WriteEvent (pout, pw, iso, f)
      end if
    case (3,4,5)
      if (present(part)) call WriteFullEvent(part(enr,:), pout, pnr, pw, iso, f)
    end select
    ! dump rejected events
    if (f == 0) return
    !!! === #3: set up kinetic variables ===
    ptot = pout(:,1) + pout(:,2)                      ! total 4-momentum
    plab = abs3(ptot)                                 ! absolute 3-momentum
    mass = sqrt(ptot(0)**2-plab**2)                   ! invariant mass
    ekin = ptot(0)-mass                               ! kinetic energy
    y = 0.5*log((ptot(0)+ptot(3))/(ptot(0)-ptot(3)))  ! rapidity
    pt = sqrt(ptot(1)**2+ptot(2)**2)                  ! transverse momentum
    mt = sqrt(mass**2+pt**2)                          ! transverse mass
    beta_gamma = plab/mass                            ! beta*gamma=p/m
    angle = op_ang (pout(:,1),pout(:,2))              ! opening angle
    hel = get_helicity (pout)                         ! helicity
    ! boost to CM frame
    call lorentz (-betaCM2Lab, ptot, "CS")
    call get_phi_Theta (ptot(1:3), theta_cm, phi)
    theta_cm = theta_cm * 180./pi
    ! additional smearing
    if (filter==1) &
      mass = rnGauss (mass*0.1, mass)  ! DLS: 10% mass resolution
    !!! === #4: put cross sections into histograms (averaged over beam energies) ===
    call AddHistMC(msigma,mass,iso,pw)
    call AddHistMC (ptsigma(0), pt, iso, pw)
    call AddHistMC (ysigma(0), y, iso, pw)
    ! determine mass bin
    mb = 1
    do while (mb < 5)
      if (mass < massBinning(mb)) exit
      mb = mb + 1
    end do
    call AddHistMC (ptsigma(mb), pt, iso, pw)
    call AddHistMC (ysigma(mb), y, iso, pw)

    if (Extra) then
      if (plab < 0.8) then
        call AddHistMC(msigma_slow,mass,iso,pw)
      else
        call AddHistMC(msigma_fast,mass,iso,pw)
      end if
      if (beta_gamma < beta_gamma_cut) &
        call AddHistMC(msigma_bgcut,mass,iso,pw)
      if (gen>1) then
        call AddHistMC(msigma_sec,mass,iso,pw)
        call AddHistMC(ptsigma_sec,pt,iso,pw)
      else
        call AddHistMC(msigma_pri,mass,iso,pw)
        call AddHistMC(ptsigma_pri,pt,iso,pw)
      end if
      call AddHistMC (psigma(0), plab, iso, pw)
      call AddHistMC (mtsigma(0), mt, iso, pw)
      call AddHistMC (tcmsigma(0), theta_cm, iso, pw)
      call AddHistMC (helsigma(0), hel, iso, pw)
      call AddHistMC (helsigma(0), 180.-hel, iso, pw)
      call AddHistMC (esigma, ekin, iso, pw)
      call AddHistMC (bgsigma, beta_gamma, iso, pw)
      call AddHistMC (oasigma, angle, iso, pw)
      call AddHistMC (dsigma, dens, iso, pw)
      call AddHistMC (psigma(mb), plab, iso, pw)
      call AddHistMC (mtsigma(mb), mt, iso, pw)
      call AddHistMC (tcmsigma(mb), theta_cm, iso, pw)
      call AddHistMC (helsigma(mb), hel, iso, pw)
      call AddHistMC (helsigma(mb), 180.-hel, iso, pw)
    end if

  end subroutine CS


  !****************************************************************************
  !****f* Dilepton_Analysis/applyFilter
  ! NAME
  ! integer function applyFilter (p1,p2)
  ! PURPOSE
  ! Performs acceptance filtering and detector resolution smearing based on
  ! the value of the 'filter' switch.
  ! INPUTS
  !  * real,dimension(0:3)::p1,p2        ---  4-momenta of positron and electron
  ! OUTPUT
  !  * returns "1" if the events is accepted, "0" if rejected
  !  * p1 and p2 will be modified if resolution smearing is applied
  !****************************************************************************
  integer function applyFilter (p1,p2) result(f)
    use constants, only: pi
    use minkowski, only: abs4, abs3, op_ang
    use random, only: rn
    use rotation, only: get_phi_Theta
    use HAFT_single, only: getHadesAcceptance, smearHadesMomentum
    use HAFT_pair, only: getHadesPairAcceptance

    real,dimension(0:3),intent(inout) :: p1,p2     ! e+/e- 4-momentum
    real,dimension(0:3) :: ptot
    real :: mass,y,pt,p1abs,p2abs,angle,theta1,theta2,phi1,phi2,acc
    integer :: s1,s2,res

    interface
      ! interface for DLS acceptance filter
      logical function beta_filter (dataset_name, xm, pt, y, weight_cut)
        character*3 :: dataset_name
        real :: xm, pt, y, weight_cut
      end function
    end interface

    f = 1
    if (filter == 0) return

    ptot = p1 + p2                                    ! total 4-momentum
    mass = abs4(ptot)                                 ! invariant mass
    y = 0.5*log((ptot(0)+ptot(3))/(ptot(0)-ptot(3)))  ! rapidity
    pt = sqrt(ptot(1)**2+ptot(2)**2)                  ! transverse momentum

    p1abs = abs3(p1)       ! abs. mom. of positron
    p2abs = abs3(p2)       ! abs. mom. of electron
    angle = op_ang(p1,p2)  ! opening angle in degrees

    select case (filter)
    case (1)   ! DLS acceptance filter version 4.1
      if (.not. beta_filter ('95A', mass, pt, y, 1000.)) f=0
    case (2)  ! simple filter for HADES
      ! cut on polar angle (in degrees)
      theta1 = acos(p1(3)/p1abs) * 180./pi
      theta2 = acos(p2(3)/p2abs) * 180./pi
      if (theta1<18. .or. theta1>85.) f = 0
      if (theta2<18. .or. theta2>85.) f = 0
      ! cut on absolute momentum
      if (p1abs<0.1 .or. p2abs<0.1) f = 0
      ! cut on opening angle (in degrees)
      if (angle<9.) f = 0
    case (3)  ! full HADES acceptance filter (using pair acceptance) + mass resolution
      acc = getHadesPairAcceptance(mass,pt,y,-2)
      if (rn()>acc) f = 0
      ! detector resolution smearing
      call smearHadesMomentum(p1,3,2)                  ! positron
      call smearHadesMomentum(p2,3,3)                  ! electron
    case (4)  ! full HADES acceptance filter (using single-particle acceptance) + opening angle cut
      call get_phi_Theta (p1(1:3),theta1,phi1)
      call get_phi_Theta (p2(1:3),theta2,phi2)
      acc = getHadesAcceptance (2,p1abs,theta1*180./pi,phi1*180./pi,-2) * &   ! positron
            getHadesAcceptance (3,p2abs,theta2*180./pi,phi2*180./pi,-2)       ! electron
      if (rn()>acc) f = 0
      ! detector resolution smearing
      if (targNuc%charge == 6) then
        res = 1  ! C+C runs have lower resolution
      else
        res = 3  ! others have higher resolution
      end if
      call smearHadesMomentum (p1, res, 2)                ! positron
      call smearHadesMomentum (p2, res, 3)                ! electron
      ! cut on opening angle (in degrees)
      angle = op_ang(p1,p2)  ! opening angle in degrees
      if (angle<9.) f = 0
      ! cut on single lepton momenta
      p1abs = abs3(p1)       ! abs. mom. of positron
      p2abs = abs3(p2)       ! abs. mom. of electron
      if (p1abs<plep_cut(1) .or. p1abs>plep_cut(2)) f = 0
      if (p2abs<plep_cut(1) .or. p2abs>plep_cut(2)) f = 0
    case (5)  ! simple filter for g7/CLAS
      call get_phi_Theta (p1(1:3),theta1,phi1)
      call get_phi_Theta (p2(1:3),theta2,phi2)
      ! cut on polar angle (in degrees)
      theta1 = theta1 * 180./pi
      theta2 = theta2 * 180./pi
      if (theta1<5. .or. theta1>35.) f = 0
      if (theta2<5. .or. theta2>35.) f = 0
      ! different-sector cut
      s1 = int(phi1 * 3./pi)  ! each sector is 60 degrees wide (=pi/3)
      s2 = int(phi2 * 3./pi)
      if (s1==s2) f = 0
      ! momentum cut
      if (p1abs<p_lep_min .or. p2abs<p_lep_min) f = 0
    case (6)  ! KEK E325 acceptance filter, using same cuts as in Ozawa et al, PRL 86 (2001) 5019
      if (y<0.6 .or. y>2.2) f = 0               ! rapidity cut
      if (pt>1.5) f = 0                         ! transverse momentum cut
      if (angle<50. .or. angle>150.) f = 0      ! opening angle cut
      ! TODO: same-segment cuts
    case (7)   ! JPARC E16
      ! nothing here for now (only resolution smearing)
    end select

  end function applyFilter


  !****************************************************************************
  !****o* Dilepton_Analysis/DileptonMass.dat
  ! NAME
  ! file DileptonMass.dat
  ! PURPOSE
  ! Dilepton mass spectrum.
  !****************************************************************************
  !****o* Dilepton_Analysis/DileptonPlab.dat
  ! NAME
  ! file DileptonPlab.dat
  ! PURPOSE
  ! Dilepton momentum spectrum.
  !****************************************************************************
  !****o* Dilepton_Analysis/DileptonPt.dat
  ! NAME
  ! file DileptonPt.dat
  ! PURPOSE
  ! Dilepton transverse momentum spectrum.
  !****************************************************************************
  !****o* Dilepton_Analysis/DileptonMt.dat
  ! NAME
  ! file DileptonMt.dat
  ! PURPOSE
  ! Dilepton transverse mass spectrum.
  !****************************************************************************
  !****o* Dilepton_Analysis/DileptonY.dat
  ! NAME
  ! file DileptonY.dat
  ! PURPOSE
  ! Dilepton rapidity spectrum.
  !****************************************************************************
  !****o* Dilepton_Analysis/DileptonEkin.dat
  ! NAME
  ! file DileptonEkin.dat
  ! PURPOSE
  ! Dilepton kinetic energy spectrum.
  !****************************************************************************


  !****************************************************************************
  !****s* Dilepton_Analysis/Dilep_write_CS
  ! NAME
  ! subroutine Dilep_write_CS
  ! PURPOSE
  ! Writes out the cross section spectra (called at the end each run).
  ! INPUTS
  ! OUTPUT
  ! The spectra are written to files called
  ! 'DileptonMass.dat', 'DileptonPlab.dat', etc
  !****************************************************************************
  subroutine Dilep_write_CS
    use histMC, only: WriteHistMC, WriteHistMC_Gauss
    use inputGeneral, only: num_Runs_SameEnergy, num_energies

    integer :: i
    character(9) :: mbstr
    integer, save :: n = 0
    real :: m

    if (.not.enable) return

    n = n + 1                                      ! count runs
    m = num_Runs_SameEnergy*num_energies/float(n)  ! multiplication factor for writing histograms

    open(998,file="inclXS.dat")
    write(998,'(A)')       "inclusive cross sections [microbarn]: rho0, omega, phi, pi0, eta, eta', Delta0, Delta+"
    write(998,'(8ES13.5)') inclXS*m
    close(998)

    ! mass differential cross section
    call WriteHistMC (msigma,     'DileptonMass.dat',      mul=m)
    do i=0,5
      mbstr = ""
      if (i>0) write(mbstr,'(A,i1)') "_massBin",i
      call WriteHistMC (ptsigma(i), 'DileptonPt'   //trim(mbstr)//'.dat', mul=m)  ! transverse momentum
      call WriteHistMC (ysigma(i),  'DileptonY'    //trim(mbstr)//'.dat', mul=m)  ! rapidity
      if (Extra) then
        call WriteHistMC (psigma(i),  'DileptonPlab' //trim(mbstr)//'.dat', mul=m)  ! momentum
        call WriteHistMC (mtsigma(i), 'DileptonMt'   //trim(mbstr)//'.dat', mul=m)  ! transverse mass
        call WriteHistMC (tcmsigma(i),'DileptonTheta'//trim(mbstr)//'.dat', mul=m)  ! theta_cm
        call WriteHistMC (helsigma(i),'DileptonHelic'//trim(mbstr)//'.dat', mul=m)  ! helicity
      end if
    end do
    if (Extra) then
      call WriteHistMC (msigma_bgcut, 'DileptonMass_bgcut.dat', mul=m)
      if (filter>5) then
        call WriteHistMC_Gauss (msigma, 'DileptonMass_smeared.dat', detector_res(filter), mul=m)
        call WriteHistMC_Gauss (msigma_bgcut, 'DileptonMass_bgcut_smeared.dat', detector_res(filter), mul=m)
      end if
      call WriteHistMC (msigma_slow,'DileptonMass_slow.dat', mul=m)
      call WriteHistMC (msigma_fast,'DileptonMass_fast.dat', mul=m)
      call WriteHistMC (msigma_pri, 'DileptonMass_pri.dat',  mul=m)
      call WriteHistMC (msigma_sec, 'DileptonMass_sec.dat',  mul=m)
      call WriteHistMC (ptsigma_pri,'DileptonPt_pri.dat',    mul=m)
      call WriteHistMC (ptsigma_sec,'DileptonPt_sec.dat',    mul=m)
      call WriteHistMC (esigma, 'DileptonEkin.dat',      mul=m)  ! kinetic energy differential cross section
      call WriteHistMC (bgsigma,'DileptonBetaGamma.dat', mul=m)  ! beta*gamma
      call WriteHistMC (oasigma,'DileptonAngle.dat',     mul=m)  ! opening angle
      call WriteHistMC (dsigma, 'DileptonDensity.dat',   mul=m)  ! density
    end if
    call WriteHistMC (deltaMass, 'DeltaMass.dat', mul=m)                   ! Delta mass distribution
    if (particle_source) call WriteHistMC (rhoMass, 'rhoMass.dat', mul=m)  ! rho mass distribution
    ! particle numbers
    call writePartNum()
  end subroutine Dilep_write_CS


  !****************************************************************************
  !****f* Dilepton_Analysis/DilepSim2
  ! NAME
  ! function DilepSim2 (ptot, massdil) result(plep)
  ! PURPOSE
  ! Simulates a dilepton decay with 2-particle final state, using isotropic
  ! angular distribution.
  ! INPUTS
  !  * real, intent(in) :: ptot(0:3)    - four momentum of decaying particle
  !  * real, intent(in) :: massdil      - invariant mass of the dilepton pair
  ! OUTPUT
  !  * real, dimension(0:3,1:2) :: plep   - 4-vectors of the generated lepton pair
  !****************************************************************************
  function DilepSim2 (ptot, massdil) result(plep)
    use lorentzTrafo, only: lorentz
    use random, only: rnOmega
    use constants, only: melec
    use twoBodyTools, only: pCM
    ! input:
    real, intent(in) :: ptot(0:3) ! total momentum
    real, intent(in) :: massdil   ! dilepton mass at creation
    ! output:
    real, dimension(0:3,1:2) :: plep
    ! local variables
    real:: pout(0:3),pfinal
    pfinal = pCM (massdil, melec, melec)
    pout(1:3) = rnOmega() * pfinal
    pout(0) = sqrt(melec**2+pfinal**2)
    call lorentz (-ptot(1:3)/ptot(0), pout, "DilepSim2")
    plep(:,1) = pout
    plep(:,2) = ptot - pout
  end function DilepSim2


  !****************************************************************************
  !****f* Dilepton_Analysis/DilepSim3
  ! NAME
  ! function DilepSim3 (ptot, massdil, mass1, ang) result (plep)
  ! PURPOSE
  ! Simulates a dilepton decay with 3-particle final state (Dalitz decay).
  ! INPUTS
  !  * real, intent(in) :: ptot(0:3)    - four momentum of decaying particle
  !  * real, intent(in) :: massdil      - invariant mass of the dilepton pair
  !  * real, intent(in) :: mass1        - mass of 2nd particle
  !  * integer, intent(in) :: ang       - angular distribution of gamma* decay: 0=isotropic, 1,2=anisotropic
  ! OUTPUT
  !  * real, dimension(0:3,1:2) :: plep   - 4-vectors of the generated lepton pair
  !****************************************************************************
  function DilepSim3 (ptot, massdil, mass1, ang) result (plep)
    use lorentzTrafo, only: lorentz
    use random, only: rnOmega, rnOmega_anis
    use constants, only: melec
    use minkowski, only: abs4
    use rotation, only: rotateTo
    use twoBodyTools, only: pCM
    ! input:
    real, intent(in)    :: ptot(0:3)    ! total momentum
    real, intent(in)    :: massdil      ! dilepton invariant mass
    real, intent(in)    :: mass1        ! mass of 2nd particle
    integer, intent(in) :: ang          ! angular distr: 0=isotropic, 1,2=anisotropic
    ! output:
    real, dimension(0:3,1:2) :: plep
    ! local variables
    real, dimension(0:3) :: pdil,pout
    real, dimension(1:3) :: betacm
    real :: pfinal

    ! select decay angle of gamma* (isotropically)
    pfinal = pCM (abs4(ptot), massdil, mass1)
    pdil(1:3) = rnOmega() * pfinal
    pdil(0) = sqrt(massdil**2+pfinal**2)   ! pdil is now 4-momentum of gamma* in parent rest frame

    ! now do decay of gamma* into e+/e-
    pfinal = pCM (massdil, melec, melec)
    if (ang>0) then
      pout(1:3) = rnOmega_anis(1.) * pfinal                     ! anisotropic angular distribution
      pout(1:3) = rotateTo (pdil(1:3), pout(1:3))               ! rotate z-Axis on the momentum direction of the decaying meson
    else
      pout(1:3) = rnOmega() * pfinal   ! isotropic
    end if
    pout(0) = sqrt(melec**2+pfinal**2)   ! pout is now 4-momentum of one lepton in cm frame of dilepton

    ! boost to parent rest frame
    call lorentz (-pdil(1:3)/pdil(0), pout, "DilepSim3(1)")
    plep(:,1) = pout
    plep(:,2) = pdil - pout

    ! boost to lab frame
    betacm = ptot(1:3)/ptot(0)
    call lorentz (-betacm, plep(:,1), "DilepSim3(2)")
    call lorentz (-betacm, plep(:,2), "DilepSim3(3)")

  end function DilepSim3


  !****************************************************************************
  !****f* Dilepton_Analysis/PseudoscalarDalitzPythia
  ! NAME
  ! function PseudoscalarDalitzPythia (ID, ptot) result (pdil)
  ! PURPOSE
  ! Performs the Dalitz decay of a pseudoscalar meson (pi0 or eta) via Pythia.
  ! Cf. Pythia 6.4 manual, chapter 13.3.1 ("Strong and electromagnetic decays"),
  ! eq. (275) and (276).
  ! INPUTS
  !  * integer, intent(in) :: ID         --- ID of decaying particle (pi0 or eta)
  !  * real, dimension(0:3), intent(in)  --- four momentum of decaying particle
  ! OUTPUT
  !  * real,dimension(0:3,1:2) :: pdil     --- four vectors of the generated lepton pair (1=positron,2=electron)
  !****************************************************************************
  function PseudoscalarDalitzPythia (ID, ptot) result (pdil)
    use hadronFormation, only: useJetSetVec
    use IdTable, only: pion, eta
    use minkowski, only: abs4

    integer, intent(in) :: ID
    real, dimension(0:3), intent(in) :: ptot
    real, dimension(0:3,1:2) :: pdil

    logical, save :: first = .true.
    integer, external :: PYCOMP
    integer :: KF, KC, i, i_el, i_pos, decay_pi, decay_eta
    real :: BR_pi(2), BR_eta (8)

    COMMON /PYDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    integer MSTU,MSTJ
    double precision PARU,PARJ
    SAVE /PYDAT1/

    COMMON /PYDAT3/ MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
    integer MDCY,MDME,KFDP
    double precision BRAT
    SAVE /PYDAT3/

    COMMON /PYJETS/ N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

    if (first) then
      !call PYLIST(12)   ! print all Pythia decay channels
      !stop

      ! save Pythia parameters before modifying them
      KC = PYCOMP (111)  ! pi0
      decay_pi = MDCY(KC,1)
      BR_pi = BRAT(MDCY(KC,2):MDCY(KC,2)+MDCY(KC,3)-1)

      KC = PYCOMP (221)  ! eta
      decay_eta = MDCY(KC,1)
      BR_eta = BRAT(MDCY(KC,2):MDCY(KC,2)+MDCY(KC,3)-1)

      first = .false.
    end if

    ! reset PYTHIA arrays
    N = 1
    K(1,1:5) = 0
    P(1,1:5) = 0.0
    V(1,1:5) = 0.0

    ! get KC code
    select case (ID)
    case (pion)
      KF = 111
    case (eta)
      KF = 221
    case default
      write(*,*) "wrong ID in PseudoscalarDalitzPythia: ", ID
      stop
    end select
    KC = PYCOMP (KF)

    ! set particle into PYTHIA array
    K(1,1:2) = (/1,KF/)
    P(1,1:5) = (/ptot(1:3),ptot(0),abs4(ptot)/)

    ! set PYTHIA parameters for decays
    MSTJ(21) = 2       ! particle decay on/off
    MDCY(KC,1) = 1     ! particle allowed to decay

    ! set branching ratios
    BRAT(MDCY(KC,2):MDCY(KC,2)+MDCY(KC,3)-1) = 0.0    ! everything to zero
    select case (ID)
    case (pion)
      BRAT(555) = 1.0    ! dilepton dalitz decay: pi0 -> gamma e+ e-
    case (eta)
      BRAT(594) = 1.0    ! dilepton dalitz decay: eta -> pi0 e+ e-
    end select

    ! perform PYTHIA decay
    if (useJetSetVec) call GetJetsetVecINIT
    call PYEXEC

    !call PYLIST(2)  ! print Pythia particle vector

    ! clean up PYTHIA output
    call PYEDIT(1)
    K(1:N,4) = 0 ! normally nLead
    K(1:N,5) = 0 ! normally number of string

    ! find electron and positron
    i_el = 0
    i_pos = 0
    do i=1,N
      select case (K(i,2))
      case (11)  ! electron
        i_el = i
      case (-11) ! positron
        i_pos = i
      end select
    end do

    ! get lepton momenta
    if (i_pos>0) then
      pdil(:,1) = (/P(i_pos,4),P(i_pos,1:3)/)
    else
      write(*,*) "error in PseudoscalarDalitzPythia: No positron found!"
      call PYLIST(2)
      stop
    end if
    if (i_el>0) then
      pdil(:,2) = (/P(i_el,4),P(i_el,1:3)/)
    else
      write(*,*) "error in PseudoscalarDalitzPythia: No electron found!"
      call PYLIST(2)
      stop
    end if

    ! restore PYTHIA parameters
    MSTJ(21) = 0              ! particle decay on/off

    KC = PYCOMP (111)  ! pi0
    MDCY(KC,1) = decay_pi
    BRAT(MDCY(KC,2):MDCY(KC,2)+MDCY(KC,3)-1) = BR_pi

    KC = PYCOMP (221)  ! eta
    MDCY(KC,1) = decay_eta
    BRAT(MDCY(KC,2):MDCY(KC,2)+MDCY(KC,3)-1) = BR_eta

  end function PseudoscalarDalitzPythia


  !****************************************************************************
  !****s* Dilepton_Analysis/Dilep_Brems
  ! NAME
  ! subroutine Dilep_Brems (part, srtfree, XS_elast, XS_tot)
  ! PURPOSE
  ! Calculate Dilepton Production via Nucleon-Nucleon (or pion-Nucleon)
  ! Bremsstrahlung in Soft-Photon-Approximation (SPA).
  ! See: Wolf et al., Nucl. Phys. A517 (1990) 615-638, chapter 4.2
  ! INPUTS
  ! * type(particle), intent(in) :: part(1:2) --- colliding particles
  ! * real, intent(in) :: srtfree             --- CoM energy
  ! * real, intent(in) :: XS_elast --- elastic cross section in mb (NN or piN)
  ! * real, intent(in) :: XS_tot   --- total cross section in mb (NN or piN)
  ! TODO
  ! * Pauli blocking ?
  !****************************************************************************
  subroutine Dilep_Brems (part, srtfree, XS_elast, XS_tot)
    use particleDefinition
    use mediumDefinition
    use IDTable, only: pion
    use constants, only: melec, alphaQED, pi
    use random, only: rnFlat, rnOmega, rnPower
    use inputGeneral, only: num_energies, numEnsembles, num_Runs_SameEnergy, eventtype, path_To_Input
    use eventtypes, only: HiPion, HeavyIon, Hadron
    use mediumModule, only: mediumAt
    use lorentzTrafo, only: lorentz,lorentzCalcBeta
    use tabulation, only: table_2dim, readTable, getValue

    type(particle), intent(in), target :: part(1:2)
    real, intent(in) :: srtfree, XS_elast, XS_tot

    real, save :: pw
    logical, save :: first = .true.
    real :: m1,m2,s,XS_bar,M,M_max,M_min,q0,q0max,q,s2,R2_1,R2_2,weight,x(2),b(3)
    real, dimension(0:3) :: qvec
    real, dimension(0:3,1:2) :: pdil
    type(particle),pointer :: charged, nuc
    integer :: ch,k,oob
    type(medium) :: mediumAtPosition
    type(table_2dim), save :: OBE_pp, OBE_pn
    logical :: success

    if (.not. Enable .or. brems==0) return

    ! initialize perturbative weight factors
    if (first) then
      select case (eventtype)
      case (HiPion)
        pw = 1./float(num_Runs_SameEnergy*num_energies*nevent)
      case default  ! hadron, heavy-ion
        pw = 1./float(numEnsembles*num_Runs_SameEnergy*num_energies*nevent)
      end select
      if (brems==2) then
        success = readTable (OBE_pp, trim(path_To_Input)//"/brems-pp.bz2")
        success = readTable (OBE_pn, trim(path_To_Input)//"/brems-pn.bz2")
      else if (brems==3) then
        success = readTable (OBE_pp, trim(path_To_Input)//"/brems-pp.bz2")
        success = readTable (OBE_pn, trim(path_To_Input)//"/brems-pn-piFF.bz2")
      end if
      first = .false.
    end if

    ch = 0
    charged => null ()
    nuc     => null ()
    select case (part(1)%ID+part(2)%ID)
    case (2)
      ! Nucleon-Nucleon
      if (part(1)%antiParticle .or. part(2)%antiParticle) return
      select case (part(1)%charge+part(2)%charge)
      case (0)
        return ! nn
      case (1)
        ch = 11 ! pn
      case (2)
        ch = 12 ! pp
      end select
      if (part(1)%charge==1) then
        charged => part(1)
        nuc => part(2)
      else if (part(2)%charge==1) then
        charged => part(2)
        nuc => part(1)
      else
        return
      end if
    case (102)
      ! pion-Nucleon
      select case (part(1)%charge+part(2)%charge)
      case (-1)
        ch = 13   ! pi(-)n
      case (0)
        ch = 14   ! pi(-)p
      case (1)
        ch = 15   ! pi(+)n
      case (2)
        ch = 16   ! pi(+)p
      end select
      if (part(1)%ID==pion) then
        charged => part(1)
        nuc => part(2)
      else if (part(2)%ID==pion) then
        charged => part(2)
        nuc => part(1)
      else
        return
      end if
      if (nuc%antiParticle) return
      if (charged%charge==0) return   ! pi(0)
    case default
      return
    end select

    if (ch==0) then
      write(*,*) "problem in Dilep_Brems: ", part%ID, part%charge, part%antiParticle
      stop
    end if

    m1 = charged%mass ! mass of charged particle (pion or nucleon)
    m2 = nuc%mass     ! mass of 2nd particle (nucleon)
    s = srtfree**2    ! center of mass energy

    M_max = srtfree - m1 - m2  ! maximum dilepton mass
    M_min = 2*melec            ! minimum dilepton mass
    if (M_max < M_min) return

    mediumAtPosition = mediumAt ((part(1)%position+part(2)%position)/2.)

    do k=1,nevent ! loop to enhance statistics

      if ((brems==2 .or. brems==3) .and. (ch==11 .or. ch==12)) then

        ! use OBE results (table lookup & interpolation)

        M  = rnFlat (M_min, M_max)                 ! dilepton mass
        q0max = (s+M**2-(m1+m2)**2) / (2*srtfree)  ! max. energy
        q0 = rnPower (-1., M, q0max)               ! dilepton energy
        q = sqrt(q0**2-M**2)                       ! dilepton momentum

        if (ch==11) then
          weight = getValue ( (/REAL(M,4), Ebeam(s)/), OBE_pn, oob) * (M_max-M_min)
        else
          weight = getValue ( (/REAL(M,4), Ebeam(s)/), OBE_pp, oob) * (M_max-M_min)
        end if

      else

       ! soft-photon approximation (SPA)

        do  ! choose mass and energy independently
          x(1:2) = (/ rnFlat (M_min, M_max), rnFlat (M_min, M_max) /)
          M  = min (x(1), x(2))     ! dilepton mass
          q0 = max (x(1), x(2))     ! dilepton energy
          s2 = s+M**2-2*srtfree*q0  ! energy of two-nucleon system after radiation, must be > (m1+m2)**2
          if (s2>(m1+m2)**2) exit
        end do

        q = sqrt(q0**2-M**2)           ! dilepton momentum

        s2 = s+M**2-2*srtfree*q0       ! energy of two-nucleon system after radiation, must be > (m1+m2)**2
        R2_1 = sqrt(1-(m1+m2)**2/s)
        R2_2 = sqrt(1-(m1+m2)**2/s2)

        XS_bar = (s-(m1+m2)**2)/(2*m1**2) * XS_elast * 1E3  ! phase-space corrected cross section (in microbarn)

        weight = alphaQED**2/(6*pi**3) * XS_bar*q/(M*q0**2) * R2_2/R2_1

      end if

      if (eventtype==HeavyIon .or. eventtype==Hadron) then
        weight = weight * pw / (XS_tot*1E3)                                        ! real particles (divide out total XS in microbarn)
      else
        weight = weight * pw * max(part(1)%perweight,part(2)%perweight) / XS_tot   ! pert. particles (take care of perweight)
      end if

      ! construct dilepton 4-vector in CMS (isotropically)
      qvec(0) = q0
      qvec(1:3) = q * rnOmega()
      ! boost from NN cms frame to calc frame
      b = lorentzCalcBeta (part(1)%momentum+part(2)%momentum)
      call lorentz (-b, qvec, "Dilep_Brems")

      ! boost from calc frame to lab
      call boost2Lab (qvec)

      pdil = DilepSim2 (qvec, M)
      call CS (pdil, ch, weight, mediumAtPosition%density, 0)
    end do

  contains

    real(4) function Ebeam (s)
      use constants, only: mN
      real, intent(in) :: s
      Ebeam = (s-2*mN**2)/(2*mN) - mN
    end function

  end subroutine Dilep_Brems



  !****************************************************************************
  !****s* Dilepton_Analysis/Dilep_BH
  ! NAME
  ! subroutine Dilep_BH(k_in,nuc)
  ! PURPOSE
  ! Calculate Dilepton Production via the Bethe-Heitler Process.
  ! See Diploma Thesis J.Weil.
  ! INPUTS
  ! * real,intent(in) :: k_in(0:3) --- photon 4-momentum
  ! * type(particle),intent(in) :: nuc --- nucleon (momentum, position and charge must be set)
  !****************************************************************************
  subroutine Dilep_BH(k_in,nuc)
    use particleDefinition, only: particle
    use random, only: rn, rnOmega
    use constants, only: alphaQED, melec, pi, GeVSquared_times_mb, mN
    use minkowski, only: SP,abs4,abs4sq
    use lorentzTrafo, only: lorentz,lorentzCalcBeta
    use pauliBlockingModule, only: pauliBlocking
    use inputGeneral, only: numEnsembles,num_Runs_SameEnergy
    use mediumDefinition, only: medium
    use mediumModule, only: mediumAt
    real,intent(in) :: k_in(0:3)
    type(particle),intent(in) :: nuc
    real,dimension(0:3) :: k,p_i,p_f,p,p_p,p_f2
    real,dimension(1:3) :: b1,b2,b3
    real :: srts,p_f_abs,p_abs,q2,lw1,lw2,M,M_min,M_max,W1,W2,weight
    integer :: i
    type(medium)::mediumAtPosition
    real,dimension(0:3,1:2) :: pdil

    if (.not. Enable) return

    srts = abs4(k_in+nuc%momentum)
    M_max = srts-mN
    M_min = 2*melec

    mediumAtPosition=mediumAt(nuc%position)

    do i=1,nevent_BH ! loop to enhance statistics
      k = k_in
      p_i = nuc%momentum
      M = rn()*(M_max-M_min)+M_min

      ! (1) boost to COM frame
      b1 = lorentzCalcBeta (p_i+k)
      call lorentz(b1,k,"BH 1.1")
      call lorentz(b1,p_i,"BH 1.2")
      ! generate random p_f
      p_f_abs=sqrt((srts**2-mN**2-M**2)**2-4*mN**2*M**2)/(2*srts)
      p_f(1:3)=rnOmega()*p_f_abs
      p_f(0)=sqrt(sum(p_f(1:3)**2)+mN**2)

      ! (2) boost p_f back to nucleus frame, check Pauli blocking
      p_f2 = p_f
      call lorentz(-b1,p_f2,"BH 2.1")
      if (pauliBlocking(p_f2,nuc%position,nuc%charge)) cycle

      ! (3) boost to dilepton rest frame
      b2 = lorentzCalcBeta (-p_f(1:3), M)
      call lorentz(b2,k,"BH 3.1")
      call lorentz(b2,p_i,"BH 3.2")
      call lorentz(b2,p_f,"BH 3.3")
      ! generate random p (electron) and p_p (positron)
      p_abs=sqrt(M**2/4-melec**2)
      p(1:3)=rnOmega()*p_abs
      p(0)=sqrt(sum(p(1:3)**2)+melec**2)
      p_p=p
      p_p(1:3)=-p_p(1:3)

      ! (4) boost to nucleon rest frame
      b3 = lorentzCalcBeta (p_i)
      call lorentz(b3,k,"BH 4.1")
      call lorentz(b3,p_i,"BH 4.2")
      call lorentz(b3,p_f,"BH 4.3")
      call lorentz(b3,p,"BH 4.4")
      call lorentz(b3,p_p,"BH 4.5")
      ! calculate tensor contraction
      q2 = abs4sq(k-p-p_p)
      lw1 = melec**2*(q2+2*melec**2)*(1/SP(p,k)**2+1/SP(p_p,k)**2) + (4*melec**4-q2**2)/SP(p,k)/SP(p_p,k) &
            -2*(q2+2*melec**2+SP(p_p,k))/SP(p,k) -2*(q2+2*melec**2+SP(p,k))/SP(p_p,k)
      lw2 = -melec**2*(q2/2*(1-2*p(0)/mN)+2*p(0)**2)/SP(p_p,k)**2 &
            -melec**2*(2*(k(0)-p(0)+q2/2/mN)*(k(0)-p(0))+q2/2)/SP(p,k)**2 &
            -2*((melec**2-q2/2)*(2*p(0)*(p(0)-k(0))+q2/2*((k(0)-2*p(0))/mN+1))-q2*k(0)**2/2)/SP(p,k)/SP(p_p,k) &
            + (q2/mN*(mN+p(0)-k(0)-q2/2/mN)+SP(p,k))/SP(p_p,k) + (q2*(1-p(0)/mN)+SP(p_p,k))/SP(p,k)

      call struct_func(nuc%charge,-q2,W1,W2)
      weight = -(lw1*W1+lw2*W2)/q2**2 * alphaQED**3/(pi**2*4*k(0)*srts)*p_f_abs*p_abs*(4*pi)**2 * (1000./GeVSquared_times_mb) &
               / float(numEnsembles*num_Runs_SameEnergy*nevent_BH)

      ! (5) boost back to nucleus frame
      call lorentz(-b3,p,"BH 5.1")
      call lorentz(-b3,p_p,"BH 5.2")
      call lorentz(-b2,p,"BH 5.3")
      call lorentz(-b2,p_p,"BH 5.4")
      call lorentz(-b1,p,"BH 5.5")
      call lorentz(-b1,p_p,"BH 5.6")
      pdil(:,1)=p
      pdil(:,2)=p_p
      call CS(pdil,17,weight,mediumAtPosition%density,0)
    end do

  contains

    ! structure functions
    subroutine struct_func (ch, t, W1, W2)
      use constants, only: mN
      use FF_QE_nucleonScattering, only: formfactors_QE
      use leptonicID, only: EM
      integer,intent(in) :: ch    ! nucleon charge
      real,intent(in) :: t        ! t in GeV
      real,intent(out) :: W1,W2   ! structure funtions W_1, W_2
      real :: tau,G_E,G_M,F1,F2

      call formfactors_QE(t,EM,ch,F1,F2,GE=G_E,GM=G_M)

      tau = t/(4*mN**2)

      W1 = 2*mN * G_M**2 * tau
      W2 = 2*mN * (G_E**2+tau*G_M**2)/(1+tau)

    end subroutine

  end subroutine


end module Dilepton_Analysis
