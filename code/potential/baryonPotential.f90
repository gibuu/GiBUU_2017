!******************************************************************************
!****m* /baryonPotentialModule
! NAME
! module baryonPotentialModule
! PURPOSE
! Includes all information about the baryonic potentials.
! Note: Only the non-relativistic case of Skyrme-like mean fields is treated
! here. For relativistic mean fields, see RMF.f90.
!******************************************************************************
module baryonPotentialModule

  implicit none
  private


  !****************************************************************************
  !****g* baryonPotentialModule/EQS_Type
  ! PURPOSE
  ! Switch for equation of state for nucleon resonances with spin 1/2.
  ! Parameters for nucleon potentials:
  ! *  0 = nucleon potential is set to zero
  ! *  1 = soft,   momentum dependent, lambda = 2.130 (Teis PhD, K = 215 MeV)
  ! *  2 = hard,   momentum dependent, lambda = 2.126 (Teis PhD, K = 380 MeV)
  ! *  3 = soft,   momentum independent               (Teis PhD, K = 215 MeV)
  ! *  4 = hard,   momentum independent               (Teis PhD, K = 380 MeV)
  ! *  5 = medium, momentum dependent, lambda = 2.130 (Teis PhD, K = 290 MeV)
  ! *  6 = LDA potential (Birger Steinmueller)
  ! *  7 = Deuterium potential Argonne V18 (not usable for eventtypes 'heavyIon' and 'hadron')
  ! *  8 = LDA Potential Welke
  ! *  9 = Buss PhD, Set#1 (K = 220 MeV, momentum dependent)
  ! * 10 = Buss PhD, Set#2 (K = 220 MeV, momentum dependent)
  ! * 11 = Buss PhD, Set#3 (K = 220 MeV, momentum dependent)
  ! * 12 = Shanghai meeting 2014 (soft, momentum independent, K = 240 MeV)
  ! * 98 = use pre-stored values
  ! * 99 = variable Skyrme : E_bind, p_0, U_0, rho_0 must be defined!
  ! NOTES
  ! Can be set in namelist baryonPotential. References:
  ! * for 1-5,  see the PhD thesis of S. Teis, chapter 3.3.2 / table 3.1
  ! * for 9-11, see the PhD thesis of O. Buss, chapter 7.2.3 / table 7.1
  ! SOURCE
  !
  integer, save :: EQS_Type = 5
  !****************************************************************************


  !****************************************************************************
  ! Parameters of potential:
  ! Units: alpha [MeV], beta [MeV], c [MeV], lambda [1/fm], eta [fm^4]
  ! 'c' and 'lambda' determine the momemtum dependence.
  ! 'eta' concerns the surface term (cf. SurfacePotFlag)
  real, parameter, dimension(1:12) :: &
        alpha  = (/ -108.619,  -9.972, -286.99993 , -124.2611 , -29.253, 0., 0., 0.,  -22.9 ,    5.72 ,   56.7 , -209.2  /), &
        beta   = (/  136.779,  38.032,  233.6517  ,   71.03007,  57.248, 0., 0., 0.,   33.9 ,   13.6  ,    4.75,  156.4  /), &
        tau    = (/    1.259,   2.404,    1.225707,    2.00108,   1.76 , 0., 0., 0.,    1.82,    2.61 ,    3.69,    1.35 /), &
        lambda = (/    2.13 ,   2.126,    0.      ,    0.     ,   2.13 , 0., 0., 0.,    0.96,    0.807,    1.18,    0.   /), &
        c      = (/  -63.601, -63.601,    0.      ,    0.     , -63.516, 0., 0., 0., -102.  , -144.   , -145.  ,    0.   /), &
        eta    = (/    0.4  ,   0.3  ,    0.4     ,    0.27   ,   0.36 , 0., 0., 0.,    0.  ,    0.   ,    0.  ,    0.   /)
  !****************************************************************************


  !****************************************************************************
  !****g* baryonPotentialModule/symmetryPotFlag
  ! PURPOSE
  ! Switch for the asymmetry term in the nucleon potential.
  ! SOURCE
  !
  integer, save :: symmetryPotFlag = 0
  ! NOTES
  ! Can be set to different value in jobcard, namelist 'baryonPotential'.
  ! Possible values:
  ! * 0 = none (default)
  ! * 1 = linear (strength given by 'dsymm')
  ! * 2 = stiffer, Esym=Esym_rho_0*U^gamma=31.*U^gamma, gamma=2
  ! * 3 = stiff, linear increasing Esym=Esym_rho_0*U=31.*U
  ! * 4 = soft, U_c=3, can give negative Esym=Esym_rho_0*U*(U_c-U)/(U_c-1)
  !****************************************************************************


  !****************************************************************************
  !****g* baryonPotentialModule/symmetryPotFlag_Delta
  ! PURPOSE
  ! Switch for the asymmetry term in the Delta potential.
  ! SOURCE
  !
  logical, save :: symmetryPotFlag_Delta = .false.
  ! NOTES
  ! Can be set to different value in jobcard, namelist 'baryonPotential'.
  ! If .true., a symmetry potential will be used also for the Delta (but only
  ! if symmetryPotFlag>0). It is closely related to the symmetry potential of
  ! the nucleon.
  !****************************************************************************


  !****************************************************************************
  !****g* baryonPotentialModule/dsymm
  ! SOURCE
  !
  real, save :: dsymm = 0.03
  ! PURPOSE
  ! Parameter for symmetry potential in GeV.
  ! NOTES
  ! Can be set in namelist baryonPotential. Value is only used for
  ! symmetryPotFlag = 1.
  !****************************************************************************


  !****************************************************************************
  !****g* baryonPotentialModule/SurfacePotFlag
  ! PURPOSE
  ! Switch for the surface term in the nucleon potential.
  ! SOURCE
  !
  logical, save :: SurfacePotFlag = .false.
  ! NOTES
  ! * Do not use it together with yukawa!
  ! * Can be set to different value in jobcard, namelist 'baryonPotential'.
  !****************************************************************************


  !****************************************************************************
  !****g* baryonPotentialModule/noPerturbativePotential
  ! PURPOSE
  ! Switch for potential of perturbative particles.
  ! If .true. then perturbative baryons feel no potential.
  ! SOURCE
  !
  logical,save :: noPerturbativePotential=.false.
  ! NOTES
  ! Can be set in namelist baryonPotential.
  !****************************************************************************


  !****************************************************************************
  !****g* baryonPotentialModule/DeltaPot
  ! PURPOSE
  ! Switch for potential of spin=3/2 resonances:
  ! * 0 = no potential
  ! * 1 = nucleon (spin=1/2) potential times 2/3 [according to Ericson/Weise book]
  ! * 2 = 100 MeV * rho/rhoNull
  ! * 3 = nucleon (spin=1/2) potential
  ! SOURCE
  !
  integer, save :: DeltaPot = 1
  ! NOTES
  ! Can be set to different value in jobcard, namelist 'baryonPotential'.
  !****************************************************************************


  !****************************************************************************
  !****g* baryonPotentialModule/HypPot
  ! PURPOSE
  ! Switch for potential of hyperons:
  ! * 0 = no potential
  ! * 1 = nucleon (spin=1/2) potential times (3+S)/3 (i.e. according to the
  !   share of the light quarks)
  ! * 2 = nucleon (spin=1/2) potential
  ! SOURCE
  !
  integer, save :: HypPot = 1
  ! NOTES
  ! Can be set to different value in jobcard, namelist 'baryonPotential'
  !****************************************************************************


  !****************************************************************************
  !****g* baryonPotentialModule/rho_0
  ! PURPOSE
  !  Nuclear matter density for EQS_Type=99
  ! SOURCE
  !
  real, save :: rho_0=0.16
  ! NOTES
  ! * Units : fm^{-3}
  !****************************************************************************

  !****************************************************************************
  !****g* baryonPotentialModule/p_0
  ! PURPOSE
  ! momentum for which U(p_0,rho=rho_0)=0 for EQS_Type=99
  ! SOURCE
  !
  real, save :: p_0 =0.8
  ! NOTES
  ! * Units : GeV
  !****************************************************************************

  !****************************************************************************
  !****g* baryonPotentialModule/U_0
  ! PURPOSE
  ! U(p=0,rho=rho_0) for EQS_Type=99
  ! SOURCE
  !
  real, save :: U_0 =0.075
  ! NOTES
  ! * Units : GeV
  !****************************************************************************

  !****************************************************************************
  !****g* baryonPotentialModule/bindingEnergy
  ! PURPOSE
  ! Nuclear matter binding energy for EQS_Type=99
  ! SOURCE
  !
  real, save :: bindingEnergy=0.016
  ! NOTES
  ! * Units : GeV
  !****************************************************************************

  !****************************************************************************
  !****g* baryonPotentialModule/compressibility
  ! PURPOSE
  ! Nuclear matter compressibility for EQS_Type=99
  ! SOURCE
  !
  real, save :: compressibility=0.290
  ! NOTES
  ! * Units : GeV
  !****************************************************************************


  !****************************************************************************
  !****g* baryonPotentialModule/nLoopReAdjust
  ! PURPOSE
  ! number of iterations, if density is readjusted
  ! (cf. type(nucleus)%ReAdjustForConstBinding)
  ! SOURCE
  !
  integer, save :: nLoopReAdjust = 10
  ! NOTES
  ! It is necessary to reiterate (at least for momentum dependent potentials),
  ! since we calculate the potential for a given pF and then calculate for the
  ! radjusting a new pF.
  !****************************************************************************


  real, dimension(:), allocatable, save :: StorePotP, StorePotN, StoreRhoB
  real, dimension(:), allocatable, save :: StoreUbP, StoreUbN
  real, save :: StorePotDX

  logical, private :: initFlag=.true.

  integer, private, save :: EQS98MomDep ! EQS=98 is momentum dependent as EQS=...


  public :: baryonPotential
  public :: rearrangementPotential
  public :: getNoPertPot_baryon
  public :: getsymmetryPotFlag_baryon
  public :: getPotentialEQSType
  public :: rhoLaplace
  public :: HandPotentialToDensityStatic


contains


  !****************************************************************************
  subroutine init
    use output, only: Write_ReadingInput
    use yukawa, only: getYukawaFlag
    use CallStack, only: TRACEBACK
    use inputGeneral, only: eventtype
    use eventtypes, only: heavyIon, hadron

    integer :: ios

    !**************************************************************************
    !****n* baryonPotentialModule/baryonPotential
    ! NAME
    ! NAMELIST /baryonPotential/
    ! PURPOSE
    ! Includes the following switches:
    ! * EQS_Type
    ! * DeltaPot
    ! * HypPot
    ! * symmetryPotFlag
    ! * symmetryPotFlag_Delta
    ! * noPerturbativePotential
    ! * rho_0
    ! * p_0
    ! * U_0
    ! * bindingEnergy
    ! * compressibility
    ! * SurfacePotFlag
    ! * nLoopReAdjust
    ! * dsymm
    !**************************************************************************
    NAMELIST /baryonPotential/ EQS_Type, DeltaPot, HypPot, &
                               symmetryPotFlag, symmetryPotFlag_Delta, noPerturbativePotential, &
                               rho_0, p_0, U_0, bindingEnergy, compressibility, SurfacePotFlag, &
                               nLoopReAdjust, dsymm

    call Write_ReadingInput('baryonPotential',0)
    rewind(5)
    read(5,nml=baryonPotential,iostat=ios)
    call Write_ReadingInput('baryonPotential',0,ios)

    write(*,*) 'Equation of state       = ', EQS_Type
    write(*,*) 'DeltaPot                = ', DeltaPot
    write(*,*) 'HypPot                  = ', HypPot
    write(*,*) 'SurfacePotFlag          = ', SurfacePotFlag
    write(*,*) 'SymmetryPotFlag         = ', symmetryPotFlag
    write(*,*) 'SymmetryPotFlag_Delta   = ', symmetryPotFlag_Delta
    write(*,*) 'dsymm [GeV]             = ', dsymm
    write(*,*) 'noPerturbativePotential = ', noPerturbativePotential

    if (EQS_Type==99) then
       write(*,*) 'Nuclear matter density [1/fm^3]      = ', rho_0
       write(*,*) 'Nuclear matter compressibility [GeV] = ', compressibility
       write(*,*) 'Nuclear matter binding energy  [GeV] = ', bindingEnergy
       write(*,*) 'U(p=0,rho=rho_0)                     = ', U_0
       write(*,*) 'p_0 for which U(p_0,rho=rho_0)=0     : ', p_0
    end if

    call Write_ReadingInput('baryonPotential',1)

    if (getYukawaFlag() .and. SurfacePotFlag) &
       call TRACEBACK('yukawa and gradient contributions can not be done together !!!!')

    if (EQS_Type==7 .and. (eventtype==heavyIon .or. eventtype==hadron)) &
      call TRACEBACK ('Deuterium potential can not be used in real-particle mode!')

    initFlag = .false.

  end subroutine init


  !****************************************************************************
  !****f* baryonPotentialModule/symEn_nuc
  ! NAME
  ! real function symEn_nuc(Q, med)
  ! PURPOSE
  ! Returns the symmetry energy for a nucleon of charge Q at the given density.
  ! Different parametrizations can be chosen via the switch 'symmetryPotFlag'.
  ! INPUTS
  ! * integer :: Q          --- charge of the nucleon
  ! * type(medium) :: med   --- medium information
  ! OUTPUT
  ! * return value: symmetry energy in GeV
  !****************************************************************************
  real function symEn_nuc(Q, med) result (esymp)
    use mediumDefinition, only: medium
    use constants, only: rhoNull
    use CallStack, only: traceback
    integer, intent(in) :: Q
    type(medium), intent(in) :: med

    real :: rhoSum, rhoDiff, U, delta, iz
    real, parameter :: eps = 0.00001

    if (symmetryPotFlag==0) then
      esymp = 0.
      return
    end if

    rhoSum  = med%density
    rhoDiff = med%densityProton-med%densityNeutron

    U     = rhoSum/rhoNull
    delta = rhoDiff/(rhoSum+eps)

    iz = 2.*float(Q) - 1.

    select case (symmetryPotFlag)
    case (1)
      esymp = dsymm * rhoDiff/rhoNull * iz
    case (2)
      ! B.A. Li, Phys.Rev. C64 (2001) 054604, or arXiv:nucl-th/0108047
      ! U=(s_0*(gamma-1)U^gamma+4.2*U^(2/3.))delta**2+-(s_0*U^gamma-12.7*U^(2./3.))*delta, gamma=2.
      esymp = (31.*U**2+4.2*U**(2./3.)) * delta**2 + iz*(31.*U**2-12.7*U**(2./3.))*delta
      esymp = esymp / 1000.
    case (3)
      ! B.A. Li, Nucl.Phys. A708 (2002) 365-390, or arXiv:nucl-th/0206053
      ! U=rho/rho_0, Uasy=+-2.*(esym(rho_0)-12.7)*U*delta for stiff
      esymp = iz*2.*(31.-12.7)*U*delta
      esymp = esymp / 1000.
    case (4)
      ! Uasy=+-2.*(15.*U_c*U-15*U^2-12.7U^(2/3))*delta+(1./3.*12.7*U^(2/3)-15*U^2)*delta**2 for soft
      esymp = iz*2.*(15.5*3.*U-15.5*U**2-12.7*U**(2./3.))*delta + (1./3.*12.7* U**(2./3.)-15.5*U**2)*delta**2
      esymp = esymp / 1000.
    case default
      write(*,*) 'symmetryPotFlag: invalid value!', symmetryPotFlag
      call traceback()
    end select
  end function symEn_nuc


  !****************************************************************************
  !****f* baryonPotentialModule/symEn_Delta
  ! NAME
  ! real function symEn_Delta(Q, med)
  ! PURPOSE
  ! Returns the symmetry energy for a Delta of charge Q (related to symEn_nuc).
  ! INPUTS
  ! * integer :: Q          --- charge of the Delta
  ! * type(medium) :: med   --- medium information
  ! OUTPUT
  ! * return value: symmetry energy in GeV
  !****************************************************************************
  real function symEn_Delta(Q, med) result (esymp)
    use mediumDefinition, only: medium
    use CallStack, only: traceback
    integer, intent(in) :: Q
    type(medium), intent(in) :: med

    if (.not. symmetryPotFlag_Delta) then
      esymp = 0.
      return
    end if

    ! B.A. Li, Nucl.Phys. A708 (2002) 365-390, or arXiv:nucl-th/0206053
    select case (Q)
    case (2) ! Delta++
      esymp = symEn_nuc(1, med)
    case (1) ! Delta+
      esymp = 2./3.*symEn_nuc (1, med) + 1./3.*symEn_nuc(0, med)
    case (0) ! Delta0
      esymp = 2./3.*symEn_nuc (0, med) + 1./3.*symEn_nuc(1, med)
    case (-1) ! Delta-
      esymp = symEn_nuc(0, med)
    case default
      write(*,*) "symEn_Delta: wrong Delta charge!", Q
      call traceback()
    end select
  end function symEn_Delta


  !****************************************************************************
  !****f* baryonPotentialModule/BaryonPotential
  ! NAME
  ! function BaryonPotential(teilchen, med, positionNotSet, EQS_in)
  ! INPUTS
  ! * type(particle) :: teilchen    -- boosted to LRF
  ! * type(medium)   :: med         -- medium information
  ! * logical        :: positionNotSet --
  !   .true. :  %position of particle is not well defined
  ! * integer, OPTIONAL :: EQS_in --
  !   If present, then we use EQS_in as EQS type, if not present
  !   then EQS is chosen according to EQS_Type.
  !
  ! Some routines like Yukawa might need the position of the particle,
  ! and not only the densities.
  ! Therefore, the positionNotSet-flag is used to check whether the position
  ! is actually set. If e.g. Yukawa is used and the position is not set,
  ! then the code stops.
  !
  ! NOTES
  ! Baryon potential is defined as 0th component of a vector potential in the
  ! LRF.
  !****************************************************************************
  real function BaryonPotential(teilchen, med, positionNotSet, EQS_in)

    use constants, only: hbarc, rhoNull
    use particleDefinition
    use ParticleProperties, only: hadron
    use yukawa, only: getYukawaFlag, getYukawaAlpha, yukpot
    use IdTable, only: nucleon, Delta
    use deuterium_PL, only: get_DeuteriumPotential
    use CallStack, only: TRACEBACK
    use mediumDefinition
    use minkowski, only: abs3

    type(particle), intent(in)          :: teilchen
    type(medium)  , intent(in)          :: med
    logical       , intent(in)          :: positionNotSet
    integer       , intent(in),optional :: EQS_in

    real :: alphaRenorm ! alpha of Skyrme potential
    real :: skyrme      ! Skyrme part of potential
    real :: out, pAbs, sqrtR
    real, dimension(1:3) :: place
    integer :: EQS,i

    if (initFlag) call init

    if (present(EQS_in)) then
      EQS = EQS_in
    else
      EQS = EQS_Type
    end if

    BaryonPotential=0. ! default fallback value

    ! No potential for perturbative particles if wished
    if (noPerturbativePotential .and. teilchen%perturbative) return

    ! Check for quark content: no potential for charmed resonances
    if (hadron(teilchen%id)%charm/=0) return

    !**************************************************************************
    ! Define nucleon potential
    !**************************************************************************
    select case (EQS)
    case (0)
       return  !=> nucleon potential= potential of all baryons

    case (6)
       if (positionNotSet) &
            call TRACEBACK('the position of the particle must be known')
       baryonPotential=LDAPotential(teilchen)
       return  !=> nucleon potential= potential of all baryons

    case (7)
       baryonPotential=get_DeuteriumPotential(teilchen)
       return  !=> nucleon potential= potential of all baryons

    case (8)
       if (positionNotSet) &
            call TRACEBACK('the position of the particle must be known')
       baryonPotential=LDAPotentialWelke(teilchen)
       return  !=> nucleon potential= potential of all baryons

    case (99)       ! Variable Skyrme
       pAbs=abs3(Teilchen%momentum)
       out = variableSkyrme(med%density*hbarc**3, pabs)

    case (1:5,9:12)
       ! Skyrme potential is renormalized due to Yukawa contribution
       alphaRenorm=alpha(EQS)+getYukawaAlpha()*1000.

       skyrme = (alphaRenorm*med%density/rhoNull &
            + beta(EQS)*(med%density/rhoNull)**tau(EQS)) / 1000.

       if (getYukawaFlag()) then
          if (positionNotSet) &
               call TRACEBACK('the position of the particle must be known')
          place=teilchen%position(1:3)
          skyrme=skyrme+yukpot(place)  !add yukawa-Term
       end if

       if (SurfacePotFlag) then
          if (positionNotSet) &
               call TRACEBACK('the position of the particle must be known')
          skyrme = skyrme + SurfacePart(teilchen,eta(EQS))
       end if

       select case (EQS)
       case (3:4,12)   ! momentum independent
          out = skyrme
       case (1:2,5,9:11) ! momentum dependent
          pAbs=abs3(Teilchen%momentum)
          out = skyrme + momentumDependentPart(pAbs,c(EQS),lambda(EQS),med%density)
       end select

    case (98) ! use the tabulated values
       if (positionNotSet) &
            call TRACEBACK('the position of the particle must be known')

       place=teilchen%position(1:3)
       sqrtR = sqrt(Dot_Product(place,place))
       ! the following is necessary, due to possible overflow from nint:
       if (sqrtR.ge.float(ubound(StorePotP,1))*StorePotDX) then
          i = ubound(StorePotP,1)
       else
          i = min(nint(sqrtR/StorePotDX),ubound(StorePotP,1))
       end if
       if (teilchen%charge>0) then
          out = StorePotP(i)
       else
          out = StorePotN(i)
       end if

       if (EQS98MomDep>0) then
          pAbs=abs3(Teilchen%momentum)
!          out = out + momentumDependentPart(pAbs,c(EQS98MomDep),lambda(EQS98MomDep),med%density)
          out = out + momentumDependentPart(pAbs,c(EQS98MomDep),lambda(EQS98MomDep),StoreRhoB(i))
       end if

    case default
       write(*,*) 'ERROR: Wrong EQS_Type in baryonPotential'
       write(*,*) 'EQS=', EQS
       call TRACEBACK()

    end select

    !**************************************************************************
    ! Take care of different spins and isospins of baryons:
    !**************************************************************************
    if (hadron(teilchen%id)%strangeness==0) then
       ! non-strange baryons
       if (abs(hadron(teilchen%ID)%Spin-0.5)<0.001) then  ! nucleon resonance
          BaryonPotential = out
          if (teilchen%ID==nucleon) &
            BaryonPotential = BaryonPotential + symEn_nuc(teilchen%charge, med)
       else if (abs(hadron(teilchen%ID)%Spin-1.5)<0.001) then  ! Delta resonance
          select case (DeltaPot)
          case (0)
             BaryonPotential=0.
          case (1)
             BaryonPotential=out*2./3.
          case (2)
             BaryonPotential=0.1*med%density/rhoNull
          case (3)
             BaryonPotential=out
          case default
             write(*,*) "bad value for deltaPot: ", DeltaPot
             stop
          end select
          if (teilchen%ID==Delta) &
            BaryonPotential = BaryonPotential + symEn_Delta (teilchen%charge, med)
       else  ! Assume same potential as for nucleon for higher spin, spin > 3/2
          BaryonPotential=out
       end if
    else
       ! hyperons (S=-1,-2,-3)
       select case (HypPot)
       case (0)
         BaryonPotential = 0.
       case (1)
         ! scale according to light quark content
         BaryonPotential = out * (3.+hadron(teilchen%id)%strangeness)/3.
       case (2)
         BaryonPotential = out
       case default
         write(*,*) "bad value for HypPot: ", HypPot
         stop
       end select
    end if

  end function BaryonPotential


  !****************************************************************************
  !****if* baryonPotentialModule/variableSkyrme
  ! NAME
  ! real function variableSkyrme(rho,p)
  ! PURPOSE
  ! * This function evaluates a variable Skyrme mean field potential
  ! * It's variable in the sense, that the potential parameters are no longer
  !   fixed by the above tables "alpha", "beta",... They are evaluated based
  !   on the input values of rhoNull, p_0, u_0, bindingEnergy and
  !   compressibility!
  !
  ! INPUTS
  ! * real :: rho -- Density in GeV**3
  ! * real :: p   -- Momentumin GeV
  !
  ! OUTPUT
  ! * Single particle potential in GeV
  !
  ! NOTES
  ! * This function is initializing the potential parameters when its called
  !   for the first time.
  ! * See Oliver's Phd thesis appendix A.4
  ! * Stops the code if there is no solution for the parameters.
  !****************************************************************************
  real function variableSkyrme(rho,p)
    use skyrme, only: evaluate_skyrmeParameters, U
    use constants, only: hbarc

    real, intent(in) :: rho ! Density in GeV**3
    real, intent(in) :: p   ! Momentumin GeV
    real, save :: a,b,c, tau, lambda,rhoNull_GeV
    logical, save :: notInitialized=.true.
    logical :: success

    if (notInitialized) then
       write(*,*) 'Initializing skyrme parameters:'
       notInitialized=.false.
       rhoNull_GeV=rho_0*hbarc**3
       call evaluate_skyrmeParameters(rhoNull_GeV, p_0, u_0, bindingEnergy, &
            compressibility,A,B,C,tau,lambda,success)
       if (.not.success) then
          write(*,*) 'Could not establish Skyrme parameters!!! STOP'
          stop
       end if
    end if
    variableSkyrme= U(rho, p, rhoNull_GeV, A,B,C, tau, lambda)
  end function variableSkyrme


  !****************************************************************************
  !****if* baryonPotentialModule/momentumDependentPart
  ! NAME
  ! real function momentumDependentPart(pin,c,lambda,rho,pF_in)
  ! PURPOSE
  ! This function provides the analytical momdep. potential.
  ! INPUTS
  ! * real :: pIn       -- absolute momentum in GeV in LRF
  ! * real :: c         -- parameter of potential
  ! * real :: lambda    -- parameter of potential
  ! * real :: rho       -- baryon density in LRF
  ! * real, OPTIONAL :: pF_in -- value of fermi mom to use (in GeV)
  ! NOTES
  ! see effenberger, dr.-thesis, pages 14-16
  !****************************************************************************
  real function momentumDependentPart(pin,c,lambda,rho,pF_in)
    use constants, only: pi, rhoNull, hbarc

    !input
    real, intent(in) ::  c
    real, intent(in) ::  lambda
    real, intent(in) ::  pIn
    real, intent(in) ::  rho
    real, intent(in), optional :: pF_in

    integer,parameter ::  isum=5

    real      pfermi, p
    real      pot, t1, t2, t3, t4, t5, t6, t7
    real      xtest, temp, zwi
    integer   isu


    p = pin/hbarc !  convert p in units 1/fm
    if (p.lt.1.0e-10) then
       p = 1.0e-10
    end if

    if (present(pF_in)) then
       pfermi = pF_in/hbarc !  convert p in units 1/fm
    else
       pfermi = (3./2.*pi*pi*rho )**(1./3.)
    end if

    !further on everything is calculated in 1 / fm
    !the constant c,lambda converts then the potential in gev

    !determine whether the small momentum expansion has to be used
    xtest = 2.0*pfermi*p/(pfermi**2+lambda**2)

    if (xtest .gt. 1.0e-06) then

       t1 = pi*lambda**3*4.0/(2.0*pi)**3
       t2 = (pfermi**2+lambda**2-p**2)/(2.0*p*lambda)
       t3 = (p+pfermi)**2+lambda**2
       t4 = (p-pfermi)**2+lambda**2
       t5 = 2.0*pfermi/lambda
       t6 = (p+pfermi)/lambda
       t7 = (p-pfermi)/lambda
       pot = t2*log(t3/t4) + t5 - 2.0*(atan(t6)-atan(t7))
       pot = pot*t1
       pot = pot*2.0*c/rhoNull

    else

       t1   = pi*lambda**3*4.0/(2.0*pi)**3
       t2   = (pfermi**2 + lambda**2 - p**2)/lambda
       t3   = 0.0
       temp = 2.0*pfermi/(pfermi**2+lambda**2)
       do isu = 0, isum
          zwi = p**float(2*isu) * temp**(float(2*isu+1))
          t3 = t3 + zwi/float(2*isu+1)
       end do
       t5 = 2.0*pfermi/lambda
       t6 = (p+pfermi)/lambda
       t7 = (p-pfermi)/lambda

       pot = t2*t3+t5-2.0*(atan(t6)-atan(t7))
       pot = pot*t1
       pot = pot*2.0*c/rhoNull

    end if
    momentumDependentPart = pot/1000. !convert to GeV
    return
  end function momentumDependentPart


  !****************************************************************************
  !****f* baryonPotentialModule/rhoLaplace
  ! NAME
  ! function rhoLaplace(rvec,a)
  ! PURPOSE
  ! Calculates div(grad(rho))
  !****************************************************************************
  real function rhoLaplace(rvec,a)
    use dichteDefinition
    use densitymodule, only: densityAt
    use minkowski, only: abs4

    real::rho,rho1p,rho1m
    real, dimension(1:3), intent(in) :: a, rvec
    real, dimension(1:3) :: rvec1p,rvec1m
    integer :: j
    type(dichte) :: pos,pos1p,pos1m

    rhoLaplace=0.
    pos=densityAt(rvec)
    rho=abs4(pos%baryon)
    do j=1,3,1
       rvec1p=rvec
       rvec1m=rvec
       rvec1p(j)=rvec1p(j)+a(j)
       rvec1m(j)=rvec1m(j)-a(j)
       pos1p=densityAt(rvec1p)
       pos1m=densityAt(rvec1m)
       rho1p=abs4(pos1p%baryon)
       rho1m=abs4(pos1m%baryon)
       rhoLaplace=rhoLaplace+(rho1p-2*rho+rho1m)/(a(j)**2)
    end do
  end function rhoLaplace


  !****************************************************************************
  !****f* baryonPotentialModule/LDApotential(teilchen)
  ! NAME
  ! function LDApotential(teilchen)
  ! PURPOSE
  ! Calculates the Potential for a nucleus initialised with the LDA approach
  !****************************************************************************
  real function LDApotential(teilchen)
    use particleDefinition
    use dichteDefinition
    use densitymodule, only: densityAt,gridSpacing
    use nucDLDA, only: startcond
    use nucleusDefinition
    use nucleus, only: getTarget
    use constants, only: hbarc

    type(tnucleus), pointer, save :: targetNucleus
    type(particle), intent(in) :: teilchen

    type(dichte)::pos
    real::rho!,r
    real,save :: b1,b2,ck,a,rho0,b3,eta,E0  ! ,eof
    real,dimension(1:3),save::stepsize
    logical,save :: firsttime=.true.

    if (firsttime) then
       targetNucleus => getTarget()
!        e0f=targetNucleus%chemPot
       call startcond(rho0,E0,ck,b3,b1,b2,a,eta)
       stepsize=gridSpacing

       firsttime=.false.
    end if

    pos=densityAt(teilchen%position)
    rho=SQRT(pos%baryon(0)**2-Dot_Product(pos%baryon(1:3),pos%baryon(1:3)))
    !r=SQRT((teilchen%position(1))**2+(teilchen%position(2))**2+(teilchen%position(3))**2)
    !no spherical symmetry
!    LDApotential=(2.*b1*rho+7./3.*b2*rho**(4./3.)+8./3.*b3*rho**(5./3.)-2*a*rhoLaplace(teilchen%position,stepsize)-e0f*E0)*hbarc/1000.
    !kleiner Test
    LDApotential=(2.*b1*rho+7./3.*b2*rho**(4./3.)+8./3.*b3*rho**(5./3.)-2*a*rhoLaplace(teilchen%position,stepsize))*hbarc

  end function LDApotential


  !****************************************************************************
  !****f* baryonPotentialModule/LDApotentialWelke(teilchen)
  ! NAME
  ! function LDApotentialWelke(teilchen)
  ! PURPOSE
  ! Calculates the Potential for a nucleus initialised with the LDA approach
  ! and additionally a momentum dependent part
  !****************************************************************************
  real function LDApotentialWelke(teilchen)
    use particleDefinition
    use dichteDefinition
    use densitymodule, only: densityAt,gridSpacing
    use nucDLDA, only: startcondWelke
    use constants, only: hbarc

    type(particle), intent(in) :: teilchen

    type(dichte)::pos
    real::rho,pAbs!,r
    real,save :: b1,b2,ck,a,rho0,b3,eta,E0,lambdaLDA,alphaLDA,CLDA!,e0f
    real,dimension(1:3),save::stepsize
    logical,save :: firsttime=.true.

    if (firsttime) then
       call startcondWelke(rho0,E0,ck,b3,b1,b2,a,eta,lambdaLDA,alphaLDA,CLDA)
       stepsize=gridSpacing

       firsttime=.false.
    end if

    pos=densityAt(teilchen%position)
    rho=SQRT(pos%baryon(0)**2-Dot_Product(pos%baryon(1:3),pos%baryon(1:3)))
    pAbs=Sqrt(Dot_Product(Teilchen%momentum(1:3),Teilchen%momentum(1:3)))
    !r=SQRT((teilchen%position(1))**2+(teilchen%position(2))**2+(teilchen%position(3))**2)
    !no spherical symmetry
!    LDApotential=(2.*b1*rho+7./3.*b2*rho**(4./3.)+8./3.*b3*rho**(5./3.)-2*a*rhoLaplace(teilchen%position,stepsize)-e0f*E0)*hbarc/1000.
    !kleiner Test
    LDApotentialWelke=(2.*b1*rho+7./3.*b2*rho**(4./3.)+8./3.*b3*rho**(5./3.)-2*a*rhoLaplace(teilchen%position,stepsize)) &
                       *hbarc + momentumDependentPart(pAbs,CLDA,lambdaLDA,rho)


  end function LDApotentialWelke


  !****************************************************************************
  !****f* baryonPotentialModule/SurfacePart
  ! NAME
  ! function SurfacePart(teilchen,spar)
  ! PURPOSE
  ! Determines the surface contribution to the total baryon potential
  ! INPUTS
  ! * type(particle) :: teilchen -- particle, in a position of which
  !   grad(\nabla\rho) has to be calculated.
  ! * real           :: spar  -- Parameter of surface part of baryon potential
  ! OUTPUT
  ! real :: SurfacePart
  !****************************************************************************
  real function SurfacePart(teilchen,spar)
    use particleDefinition
    use dichteDefinition
    use constants, only: hbarc
    use densitymodule, only: gridSpacing

    type(particle), intent(in) :: teilchen
    real,           intent(in) :: spar

    SurfacePart = -2.* spar * rhoLaplace(teilchen%position,gridSpacing)
    SurfacePart = SurfacePart * hbarc
  end function SurfacePart


  !****************************************************************************
  !****f* baryonPotentialModule/rearrangementPotential
  ! NAME
  ! real function rearrangementPotential(teilchen, med)
  ! PURPOSE
  ! Returns the value of the rearrangement potential.
  ! INPUTS
  ! * type(particle) :: teilchen -- particle whose rearrangement potential
  !   should be calculated. It should be boosted to LRF.
  !   The position of the particle must be set!!!
  ! * type(medium) :: med -- density information
  ! OUTPUT
  ! function value
  ! NOTES
  ! Notation according to Teis Dr.-thesis. Pages 76-78.
  !
  ! This is *not* the rearrangement potential, but the expression
  !   -(1/2 U_b + U_r)
  ! which is needed to calculate the binding energy correctly.
  !****************************************************************************
  real function rearrangementPotential(teilchen, med)
    use particleDefinition
    use mediumDefinition
    use ParticleProperties, only: hadron
    use yukawa, only: getYukawaFlag, getYukawaAlpha, yukpot
    use idTable, only: nucleon
    use constants, only: rhoNull
    use CallStack, only: traceback

    type(particle), intent(in) :: teilchen
    type(medium), intent(in) :: med

    real :: rhoDiff,rhoSum      ! rho_p-rho_n, rho_p+rho_n
    real :: alphaRenorm         ! alpha of Skyrme potential
    real :: skyrme_b, skyrme_r  ! Skyrme part of rearrangement potential
    real :: pAbs, sqrtR
    integer :: EQS, i

    if (initFlag) call init

    rearrangementPotential = 0. ! default fallback value

    ! no potential for perturbative particles if wished
    if (noPerturbativePotential.and.teilchen%perturbative) return

    !Check for quark content:
    !No potential for s<-1 resonances and charmed resonances
    if (hadron(teilchen%id)%strangeness<-1) return
    if (hadron(teilchen%id)%charm/=0) return

    select case (EQS_Type)
    case (0,7,99)
       return
    case (1:5,9:12)
       ! Skyrme potential is renormalized due to Yukawa contribution
       alphaRenorm=alpha(EQS_Type)+getYukawaAlpha()*1000.

       skyrme_b = (alphaRenorm*(med%density/rhoNull) &
                + 2./(tau(EQS_Type)+1.)*beta(EQS_Type) &
                  *(med%density/rhoNull)**tau(EQS_Type))/1000.

       if (getYukawaFlag()) &
            skyrme_b = skyrme_b + yukpot(teilchen%position(1:3))

       if (SurfacePotFlag) &
            skyrme_b=skyrme_b + SurfacePart(teilchen,eta(EQS_Type))

       skyrme_r = ((tau(EQS_Type)-1.)/(tau(EQS_Type)+1.)*beta(EQS_Type) &
                  *(med%density/rhoNull)**tau(EQS_Type))/1000.

       rearrangementPotential = -skyrme_b/2-skyrme_r

       if (symmetryPotFlag>0 .and. teilchen%id==nucleon .and. .not.teilchen%antiparticle) then
          ! add symmetry potential
          if (symmetryPotFlag == 1) then
            rhoDiff = med%densityProton - med%densityNeutron
            rhoSum  = med%densityProton + med%densityNeutron
            rearrangementPotential = rearrangementPotential &
                                  + 0.5*dsymm*rhoDiff**2/(rhoNull*rhoSum) &
                                  - dsymm*rhoDiff/rhoNull*(float(teilchen%charge)-0.5)*2.
          else
            ! TODO: handle symmetryPotFlag>1
            call traceback("error: rearrangementPotential not implemented for symmetryPotFlag>1!")
          end if
       end if

       EQS = EQS_Type

    case (98)
       if (teilchen%id==nucleon .and. .not.teilchen%antiparticle) then
          sqrtR = sqrt(sum(teilchen%position(1:3)**2))
          ! the following is necessary, due to possible overflow from nint:
          if (sqrtR>=float(ubound(StorePotP,1))*StorePotDX) then
             i = ubound(StorePotP,1)
          else
             i = min(nint(sqrtR/StorePotDX),ubound(StorePotP,1))
          end if
          if (teilchen%charge==1) then
             rearrangementPotential = -(StorePotP(i)-StoreUbP(i)/2)
          else
             rearrangementPotential = -(StorePotN(i)-StoreUbN(i)/2)
          end if
       else
          ! do not know yet what to do here...
       end if
       EQS = EQS98MomDep

    end select

    select case (EQS)
    case (1:2,5,9:11) ! momentum dependent
       pAbs=Sqrt(Dot_Product(Teilchen%momentum(1:3),Teilchen%momentum(1:3)))
       rearrangementPotential = rearrangementPotential &
                              - momentumDependentPart(pAbs,c(EQS),lambda(EQS),med%density)/2
    end select

    if (hadron(teilchen%id)%strangeness==0) then
      if (abs(hadron(teilchen%ID)%Spin-1.5)<0.001) then  ! Delta resonance
         select case (deltaPot)
         case (1)
            rearrangementPotential=rearrangementPotential*2./3.
         case (2)
            rearrangementPotential=0.
         end select
      end if
    else
!     Multiply potential for s=-1 baryons by 2./3. according to the
!     share of the light quarks:
      rearrangementPotential=rearrangementPotential*2./3.
    end if

  end function rearrangementPotential


  !****************************************************************************
  !****f* baryonPotentialModule/HandPotentialToDensityStatic
  ! NAME
  ! subroutine HandPotentialToDensityStatic(nuc)
  ! PURPOSE
  ! This routine tabulates the proton and neutron density as a function of
  ! the radius. It also tabulates the Coulomb potential and the baryon
  ! potentials for protons and neutrons *for the corresponding Fermi momentum*.
  ! With this, the routine which actually readjusts the density according
  ! constant Fermi energy is called.
  ! This is repeated several times, until some convergence is believed to
  ! happen. (This is necessary, since the local thomas fermi momentum depends
  ! on the density.) (Also Coulomb depends on the density, but this is of
  ! minor importance here.)
  ! If the potential is not momentum dependent, no iteration would be
  ! necessary.
  !
  ! Some complications are due to the calculation of the baryon potential in
  ! baryonPotentialModule/BaryonPotential and the used tabulations.
  !
  ! INPUTS
  ! * type(tNucleus) :: nuc
  ! OUTPUT
  ! * The potential at input is frozen and stored in a r dependent grid
  ! * EQS_Type is set to 98
  ! * The static density is adjusted.
  !****************************************************************************
  subroutine HandPotentialToDensityStatic(nuc)
    use nucleusDefinition
    use densityStatic, only: ReAdjust
    use constants, only: hbarc, mN, pi
    use particleDefinition
    use mediumDefinition
    use coulomb, only: getCoulombFlag
    use output, only: Write_InitStatus

    type(tnucleus), pointer :: nuc

    integer :: i,iLoop
    real :: x, rP,rN, pF, y
    type(particle) :: Part
    real, dimension(:), allocatable :: potP, potN, potC
    type(medium)         :: med
    real :: IntCoulomb1, IntCoulomb2
    logical :: useCoulomb

    logical, parameter :: verbose = .false.

    call Write_InitStatus("Readjust density distribution",0)

    if (initFlag) call init

    useCoulomb = getCoulombFlag()

    call setToDefault(Part)

    allocate(StorePotP(0:nuc%MaxIndex))
    allocate(StorePotN(0:nuc%MaxIndex))
    allocate(StoreRhoB(0:nuc%MaxIndex))
    StorePotDX = nuc%dx

    allocate(PotP(0:nuc%MaxIndex))
    allocate(PotN(0:nuc%MaxIndex))
    allocate(PotC(0:nuc%MaxIndex))

    PotC = 0.0

    Part%ID=1 ! Nucleon
    Part%mass=mN
    Part%position = 0

    select case (EQS_Type)
    case (1:2,5,9:11)
      EQS98MomDep = EQS_Type
    case default
      EQS98MomDep = 0
    end select

    med%temperature    =0.
    med%useMedium      =.true.

    do iLoop=1,nLoopReAdjust
       write(*,*) 'Readjusting..., Iteration=',iLoop

       IntCoulomb1 = 0.0
       IntCoulomb2 = 0.0

       do i=0,nuc%MaxIndex
          x = i*nuc%dx
          Part%position(1) = x

          rP = nuc%densTab(i,1)
          rN = nuc%densTab(i,2)

          med%density        = rP+rN
          med%densityProton  = rP
          med%densityNeutron = rN

          ! Coulomb potential is calcuated by integrating the charge over the
          ! included volume:
          if (useCoulomb) then
             IntCoulomb1=IntCoulomb1 + x**2*rP
             if (i.gt.0) then
                IntCoulomb2=IntCoulomb2 + IntCoulomb1/x**2
             else
                IntCoulomb2=IntCoulomb2 + rP
             end if
             PotC(i) = IntCoulomb2 * (-4.0*pi*1./137.*hbarc) * nuc%dx**2 ! in GeV
          end if

          if (EQS_Type==98) then ! iLoop=2,3,...
             ! now we undo the rescaling, because we need the density only for
             ! calculating the fermi momentum
             rP = rP/nuc%facP
             rN = rN/nuc%facN
          end if


          pF = (3*pi**2*(rP+rN)/2)**(1./3.)*hbarc
          Part%momentum(1) = pF
          Part%momentum(0) = sqrt(Part%mass**2+pF**2)

          Part%charge = 0
!!$       pF = (3*pi**2*(rN))**(1./3.)*hbarc
!!$       Part%momentum(1) = pF
!!$       Part%momentum(0) = sqrt(Part%mass**2+pF**2)

          x = BaryonPotential(Part,med,.false.)
          PotN(i) = x
          if (iLoop.eq.1) then
             if (EQS98MomDep>0) then
                x = x - momentumDependentPart(pF,c(EQS98MomDep),lambda(EQS98MomDep),rP+rN)
             end if
             StorePotN(i) = x
          end if

          Part%charge = 1
!!$       pF = (3*pi**2*(rP))**(1./3.)*hbarc
!!$       Part%momentum(1) = pF
!!$       Part%momentum(0) = sqrt(Part%mass**2+pF**2)
          x = BaryonPotential(Part,med,.false.)
          PotP(i) = x
          if (iLoop.eq.1) then
             if (EQS98MomDep>0) then
                x = x - momentumDependentPart(pF,c(EQS98MomDep),lambda(EQS98MomDep),rP+rN)
             end if
             StorePotP(i) = x
             StoreRhoB(i) = rP+rN
          end if

!          write(*,'(6f13.5)') x,rN,rP,pF,PotN(i),PotP(i)

       end do

       if (useCoulomb) then
          x = nuc%MaxIndex*nuc%dx
          y = nuc%charge*(1./137.*hbarc)/x - PotC(nuc%MaxIndex)
          PotC = PotC + y
       end if

       ! print at least the output at the last step
       if ((verbose).or.(iLoop == nLoopReAdjust)) then
          call PlotPot
       end if

       call ReAdjust(nuc, PotP, PotN, PotC)

       EQS_Type = 98 ! we have to know this already in the next loop

    end do

    call CalculateRearrangePot

    deallocate(PotP)
    deallocate(PotN)
    deallocate(PotC)

!    write(*,*) 'stop'
!    stop

    call Write_InitStatus("Readjust density distribution",1)

  contains

    subroutine CalculateRearrangePot
      real :: VbigP, VbigN
      integer :: i
      real :: x, drB, rB, rB0

      allocate(StoreUbP(0:nuc%MaxIndex))
      allocate(StoreUbN(0:nuc%MaxIndex))

      if (verbose) then
         open(213,file='ReAdjust.Rearrange.dat', status='unknown')
         rewind(213)
      end if

      ! we calculate V(r) = /int dr drho/dr U(r)
      ! and set V(r) = rho^2/2 v(r)
      ! since U = dV/drho, we set U_b = rho*v and U_r = U-U_b

      VbigP = 0
      VbigN = 0
      rB0 = 0
      do i=nuc%MaxIndex-1,0,-1
         x = i*nuc%dx
         rB = nuc%densTab(i,1)+nuc%densTab(i,2)
         if (rB.eq.0) cycle
         drB = rB - rB0

         VbigP = VbigP+drB*StorePotP(i)
         StoreUbP(i) = rB*2*VbigP/rB**2
         VbigN = VbigN+drB*StorePotN(i)
         StoreUbN(i) = rB*2*VbigN/rB**2

         if (verbose) then
            write(213,'(i5,10f12.5)') i,x,rB, &
                 & StoreUbP(i)*1000,StoreUbN(i)*1000, &
                 & StorePotP(i)*1000,StorePotN(i)*1000
         end if
         rB0 = rB
      end do
      if (verbose) then
         close(213)
      end if
    end subroutine CalculateRearrangePot

    subroutine PlotPot
      use output, only: IntToChar
      use constants, only: rhoNull

      integer :: iMom, EQS
      integer, parameter :: nMom = 100
      real, parameter :: DeltaMom=0.020
      real :: xN, xP, xMom, xMomDep, rea

      open(113,file='ReAdjust.PlotPot.'//IntToChar(iLoop)//'.dat', status='unknown')
      rewind(113)

      write(113,'("# ",A)') "density and potential as function of position and momentum"
      write(113,'("# ",A)') "after the "//IntToChar(iLoop)//"th iteration step."
      write(113,'("#",A12,10A13)')
      write(113,'("#",A12,10A13)') 'x (pos)','p (mom)','potP', 'potN',&
           'momdep','densP','densN','rearrange'
      write(113,'("#",A12,10A13)')

      if (EQS_Type==98) then
        EQS = EQS98MomDep
      else
        EQS = EQS_Type
      end if

      xMomDep = 0
      rea = 0
      do i=0,nuc%MaxIndex,5
          x = i*nuc%dx
          Part%position(1) = x

          rP = nuc%densTab(i,1)
          rN = nuc%densTab(i,2)

          med%density        = rP+rN
          med%densityProton  = rP
          med%densityNeutron = rN

          do iMom=0,nMom
             xMom = iMom*DeltaMom
             Part%momentum(1) = xMom
             Part%momentum(0) = sqrt(Part%mass**2+xMom**2)

             Part%charge = 1
             xP = BaryonPotential(Part,med,.false.)

             if (EQS98MomDep>0) then
                xMomDep = momentumDependentPart(xMom,c(EQS98MomDep),lambda(EQS98MomDep),rP+rN)
             end if

             Part%charge = 0
             xN = BaryonPotential(Part,med,.false.)

             if (EQS>0) then
                rea = (alpha(EQS)*((rP+rN)/rhoNull) &
                     & + 2./(tau(EQS)+1.)*beta(EQS) &
                     & *((rP+rN)/rhoNull)**tau(EQS))/1000.
             end if

             write(113,'(10f13.5)') x,xMom,xP,xN,xMomDep,rP,rN,rea
          end do
          write(113,*)

       end do

      close(113)
    end subroutine PlotPot

  end subroutine HandPotentialToDensityStatic


  !****************************************************************************
  !****f* baryonPotentialModule/getNoPertPot_baryon
  ! NAME
  ! logical function getNoPertPot_baryon
  ! PURPOSE
  ! Returns the flag noPerturbativePotential
  !****************************************************************************
  logical function getNoPertPot_baryon()
    if (initFlag) call init
    getNoPertPot_baryon=noPerturbativePotential
  end function getNoPertPot_baryon


  !****************************************************************************
  !****f* baryonPotentialModule/getsymmetryPotFlag_baryon
  ! NAME
  ! logical function getsymmetryPotFlag_baryon
  ! PURPOSE
  ! Returns the flag symmetryPotFlag
  !****************************************************************************
  logical function getsymmetryPotFlag_baryon()
    if (initFlag) call init
    getsymmetryPotFlag_baryon=(symmetryPotFlag>0)
  end function getsymmetryPotFlag_baryon


  !****************************************************************************
  !****f* baryonPotentialModule/getPotentialEQSType
  ! NAME
  ! integer function getPotentialEQSType
  ! PURPOSE
  ! Returns the EQS type of the potential (0 = no potential)
  !****************************************************************************
  integer function getPotentialEQSType()
    if (initFlag) call init
    getPotentialEQSType = EQS_Type
  end function getPotentialEQSType


end module baryonPotentialModule
