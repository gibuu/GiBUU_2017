!******************************************************************************
!****m* /mesonWidthMedium
! NAME
! module mesonWidthMedium
! PURPOSE
! Basically this module calls the subroutines of the module "mesonWidth"
! and adds medium modifications to it.
! USES
! module mesonWidth
!******************************************************************************
module mesonWidthMedium
  use IDtable, only: pion,nMes
  use MassAssInfoDefinition

  implicit none
  private


  !****************************************************************************
  !****g* mesonWidthMedium/mediumSwitch
  ! SOURCE
  !
  integer, save :: mediumSwitch = 0
  ! PURPOSE
  ! Treatment of In-Medium Widths for mesons:
  ! * 0: Only vacuum widths are used.
  ! * 1: The collisional width is assumed to be constant
  !   (only density-dependent).
  ! * 2: The full tabulated in-medium width is used, as
  !   calculated via the collision term.
  !****************************************************************************

  !****************************************************************************
  !****g* mesonWidthMedium/Gamma_coll_rho
  ! SOURCE
  !
  real, save :: Gamma_coll_rho   = 0.150
  ! PURPOSE
  ! Collisional width for the rho meson in GeV.
  ! Only used if mediumSwitch = 1.
  !****************************************************************************

  !****************************************************************************
  !****g* mesonWidthMedium/Gamma_coll_omega
  ! SOURCE
  !
  real, save :: Gamma_coll_omega = 0.150
  ! PURPOSE
  ! Collisional width for the omega meson in GeV.
  ! Only used if mediumSwitch = 1.
  !****************************************************************************

  !****************************************************************************
  !****g* mesonWidthMedium/Gamma_coll_phi
  ! SOURCE
  !
  real, save :: Gamma_coll_phi   = 0.030
  ! PURPOSE
  ! Collisional width for the phi meson in GeV.
  ! Only used if mediumSwitch = 1.
  !****************************************************************************

  !****************************************************************************
  !****g* mesonWidthMedium/verboseInit
  ! SOURCE
  !
  logical, save :: verboseInit = .false.
  ! PURPOSE
  ! switch on/off informational messages during initialization
  !****************************************************************************

  !****************************************************************************
  !****g* mesonWidthMedium/allowMix
  ! SOURCE
  !
  logical, save :: allowMix = .false.
  ! PURPOSE
  ! switch on/off linear interpolation between bins in density while
  ! returning the tabulated values for MassAssInfo.
  !****************************************************************************

  logical, save:: initFlag=.true.

  logical,parameter :: debug=.false.

  integer, parameter :: nG1 =  50             ! number of density bins per rho0
  integer, parameter :: nG2 = 200             ! total number of density bins
  real, parameter :: rho_max = real(nG2)/nG1  ! maximum density which can be handled (in units of rho0)

  type(tMassAssInfoArray), dimension(pion:pion+nMes-1), save :: MassAssInfo_meson



!!$  !************************************************************************
!!$  !****f* mesonWidthMedium/partialWidthMesonMedium
!!$  ! NAME
!!$  ! real function partialWidthMesonMedium(ID, mass, IDout1,IDout2,IDout3, charge1,charge2,charge3, momLRF, mediumAtPos)
!!$  !
!!$  ! real function partialWidthMesonMedium(ID, mass, IDout1,IDout2, momLRF, mediumAtPos)
!!$  !
!!$  ! PURPOSE
!!$  ! Calculate the (partial) width energy dependent of all meson resonances
!!$  ! for a 2-body or a 3-body channel.
!!$  !
!!$  ! INPUTS
!!$  ! * integer :: ID -- id of resonance
!!$  ! * real    :: mass       -- p_mu p^mu = mass of the resonance (offshell)
!!$  ! * real, dimension(0:3) :: momLRF -- boost vector
!!$  ! * type(medium)         :: mediumAtPos -- medium information
!!$  !
!!$  ! with:
!!$  ! * integer :: IDout1, IDout2,IDout3 --
!!$  !   ID's of the decay products which one is interested in
!!$  ! * integer :: charge1, charge2, charge3  --
!!$  !   charges of the decay products which one is interested in
!!$  ! or:
!!$  ! * integer :: IDout1, IDout2 --
!!$  !   ID's of the decay products which one is interested in
!!$  !
!!$  ! USAGE
!!$  ! * 1.) first line if you are interested in a three-body process.
!!$  ! * 2.) second line if you are interested in a two-body process.
!!$  ! USES
!!$  ! Uses as an interface the module procedures partialWidthMesonMedium2body
!!$  ! and partialWidthMesonMedium3body
!!$  !************************************************************************
!!$  Interface partialWidthMesonMedium
!!$     Module Procedure partialWidthMesonMedium2body, partialWidthMesonMedium3body
!!$  End Interface

  public :: WidthMesonMedium
!!$  PUBLIC :: partialWidthMesonMedium
  public :: decayWidthMesonMedium
  public :: get_MediumSwitchMesons
  public :: GetMassAssInfo_Meson


contains


  !****************************************************************************
  !****s* mesonWidthMedium/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "width_Meson".
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput
    use CALLSTACK, only: TRACEBACK

    integer :: ios

    !**************************************************************************
    !****n* mesonWidthMedium/width_Meson
    ! NAME
    ! NAMELIST width_Meson
    ! PURPOSE
    ! Includes the switches:
    ! * mediumSwitch
    ! * Gamma_coll_rho
    ! * Gamma_coll_omega
    ! * Gamma_coll_phi
    ! * verboseInit
    ! * allowMix
    !**************************************************************************
    NAMELIST /width_Meson/ mediumSwitch, Gamma_coll_rho, Gamma_coll_omega, Gamma_coll_phi, &
                           verboseInit, allowMix

    call Write_ReadingInput('width_Meson',0)
    rewind(5)
    read(5,nml=width_Meson,IOSTAT=IOS)
    call Write_ReadingInput('width_Meson',0,IOS)

    write(*,*) 'Use in medium width for the mesons '  ,mediumSwitch

    select case (mediumSwitch)
    case (0)
       write(*,*) '    i.e. use vacuum widths'

    case (1)
       write(*,*) '    i.e. use constant collisional widths'
       write(*,*)
       write(*,*) 'Collisional width for RHO meson  : ',Gamma_coll_rho
       write(*,*) 'Collisional width for OMEGA meson: ',Gamma_coll_omega
       write(*,*) 'Collisional width for PHI meson  : ',Gamma_coll_phi
       write(*,*)
       write(*,*) 'interpolate between density bins : ',allowMix

    case (2)
       write(*,*) '    i.e. use full tabulated widths'

    case default
       call TRACEBACK()

    end select

    call InitMassAssInfo_Meson()

    call Write_ReadingInput('width_Meson',1)

    initFlag=.false.

  end subroutine readInput

  !****************************************************************************
  !****f* mesonWidthMedium/get_MediumSwitchMesons
  ! NAME
  ! logical function get_MediumSwitchMesons
  ! PURPOSE
  ! return true, if mediumSwitch>0
  !****************************************************************************
  logical function get_MediumSwitchMesons()
    if (initFlag) call readInput
    get_MediumSwitchMesons=(mediumSwitch>0)
  end function


!!$  !************************************************************************
!!$  ! c.f. interface 'partialWidthMesonMedium'
!!$  !************************************************************************
!!$  real function partialWidthMesonMedium3Body(ID,mass, &
!!$       & IDout1,IDout2,IDout3,charge1,charge2,charge3,&
!!$       & momLRF,mediumAtPos)
!!$
!!$    use particleProperties
!!$
!!$    real,intent(in) :: mass
!!$    integer,intent(in) :: ID
!!$    integer, intent(in) :: IDout1, IDout2,IDout3
!!$    integer, intent(in) :: charge1, charge2, charge3
!!$
!!$    ! Only Needed for InMediumModifications:
!!$    real,intent(in),dimension(0:3)  :: momLRF
!!$    type(medium),intent(in) :: mediumAtPos
!!$
!!$
!!$    If(initFlag) call readInput
!!$
!!$    ! (1) Check Input
!!$    if (.not.IsMeson(ID)) call TRACEBACK()
!!$
!!$    ! Since no medium modifications yet implemented :
!!$    partialwidthMesonMedium3Body=partialwidthMeson(ID,mass,IDout1,IDout2,IDout3,charge1,charge2,charge3)
!!$
!!$  end function partialwidthMesonMedium3Body
!!$  !-------------------------------------------------------------------------
!!$  real function partialWidthMesonMedium2Body(ID,mass,&
!!$       & IDout1,IDout2,&
!!$       & momLRF,mediumAtPos)
!!$
!!$    use particleProperties
!!$
!!$    real,intent(in) :: mass
!!$    integer,intent(in) :: ID
!!$    integer, intent(in) :: IDout1, IDout2
!!$
!!$    ! Only Needed for InMediumModifications:
!!$    real,intent(in),dimension(0:3)  :: momLRF
!!$    type(medium),intent(in) :: mediumAtPos
!!$
!!$    If(initFlag) call readInput
!!$
!!$    ! (1) Check Input
!!$    if (.not.IsMeson(ID)) call TRACEBACK()
!!$
!!$    ! Since no medium modifications yet implemented :
!!$    partialwidthMesonMedium2Body=partialwidthMeson(ID,mass,IDout1,IDout2)
!!$
!!$  end function partialwidthMesonMedium2Body



  !****************************************************************************
  !****s* mesonWidthMedium/decayWidthMesonMedium
  ! NAME
  ! real function decayWidthMesonMedium (ID, mass, ch, pauliFlag) result (decayWidth)
  ! PURPOSE
  ! This function returns the partial out widths with energy dependence.
  ! Also medium modifications of these out-width are included here. If the out-width has already
  ! Pauli-Blocking effects included, then the pauliFlag is set to .true. .
  ! INPUTS
  ! * integer              :: ID    -- ID of resonance
  ! * real                 :: mass  -- baremass of the resonance (offshell)
  ! * integer              :: ch    -- charge of resonance
  ! OUTPUT
  ! * real, dimension(:)   :: decayWidth  --
  !   array of ecay widths for all decay channels
  ! * logical :: pauliFlag    -- .true. if pauli-Blocking is already
  !   considered in the decay width description
  !****************************************************************************
  function decayWidthMesonMedium (ID, mass, ch, pauliFlag) result (decayWidth)
    use mediumDefinition
    use particleProperties, only: nDecays
    use mesonWidth, only: decayWidthMeson

    integer, intent(in)  :: ID
    real,    intent(in)  :: mass
    integer, intent(in)  :: ch
    logical, intent(out) :: pauliFlag
    real, dimension(1:nDecays) :: decayWidth

    ! Only Needed for InMediumModifications:
!     type(medium),intent(in) ::  mediumAtPos
!     real,intent(in),dimension(0:3) :: momLRF

    if (initFlag) call readInput

    ! No medium modification
    pauliFlag = .false.
    decayWidth = decayWidthMeson (ID, mass, ch)

  end function decayWidthMesonMedium


  !****************************************************************************
  !****s* mesonWidthMedium/WidthMesonMedium
  ! NAME
  ! real function WidthMesonMedium(ID,mass,momLRF,mediumAtPos)
  ! PURPOSE
  ! This function calculates the (partial) width energy dependent of
  ! all meson resonances with medium modifications.
  ! INPUTS
  ! * integer :: ID   -- id of resonance
  ! * real    :: mass -- p_mu p^mu = mass of the resonance (offshell)
  ! * real,dimension(0:3)  :: momLRF --
  ! * type(medium)         :: mediumAtPos --
  !****************************************************************************
  real function WidthMesonMedium(ID,mass,momLRF,mediumAtPos)
    use mediumDefinition
    use mesonWidth, only: FullWidthMeson
    real,intent(in) ::  mass
    integer,intent(in) :: ID
    real,intent(in),dimension(0:3)  :: momLRF
    type(medium),intent(in) :: mediumAtPos

    if (initFlag) call readInput

    WidthMesonMedium = FullWidthMeson(ID,mass) ! = Gamma_vac

    if (.not. mediumAtPos%useMedium) return

    WidthMesonMedium = WidthMesonMedium  & ! Gamma_tot = Gamma_vac + Gamma_coll
                       + WidthMesonMedium_GammaColl(ID,mass,momLRF,mediumAtPos)

  end function WidthMesonMedium


  real function GammaColl (ID)
    use IDTable, only: rho, omegaMeson, phi
    integer, intent(in) :: ID

    if (mediumSwitch == 0) then
      GammaColl = 0.
      return
    end if

    select case (ID)
    case (rho)
      GammaColl = Gamma_coll_rho
    case (omegaMeson)
      GammaColl = Gamma_coll_omega
    case (phi)
      GammaColl = Gamma_coll_phi
    case default
      GammaColl = 0.
    end select

  end function GammaColl


  !****************************************************************************
  !****s* mesonWidthMedium/WidthMesonMedium_GammaColl
  ! NAME
  ! real function WidthMesonMedium_GammaColl(ID,mass,momLRF,mediumAtPos)
  ! PURPOSE
  ! This function calculates the collisional width energy dependent of
  ! all meson resonances, density dependend.
  ! INPUTS
  ! * integer :: ID   -- id of resonance
  ! * real    :: mass -- p_mu p^mu = mass of the resonance (offshell)
  ! * real,dimension(0:3)  :: momLRF --
  ! * type(medium)         :: mediumAtPos --
  !****************************************************************************
  real function WidthMesonMedium_GammaColl (ID, mass, momLRF, mediumAtPos)
    use mediumDefinition
    use mesonWidthMedium_tables, only: get_inMediumWidth
    use constants, only: rhoNull

    integer,      intent(in) :: ID
    real,         intent(in) :: mass
    real,         intent(in) :: momLRF(0:3)
    type(medium), intent(in) :: mediumAtPos

    real :: dens, absP

    select case (mediumSwitch)
    case (0)
      WidthMesonMedium_GammaColl = 0.
    case (1) ! constant collisional width
       ! Collisional Width: Gamma = Gamma_0 * rho/rho_0
       dens = mediumAtPos%density/rhoNull
       if (Get_UseMassAssInfo()) then
         dens = min (dens, rho_max)   ! avoid problems in Heavy Ion Collisions
         if (mediumAtPos%density/rhoNull>rho_max) &
           write(*,*) "Warning in WidthMesonMedium_GammaColl: density too large!", mediumAtPos%density/rhoNull
       end if

       WidthMesonMedium_GammaColl = GammaColl(ID) * dens

    case (2) ! full in-medium width
       absP=sqrt(Dot_Product(momLRF(1:3),momLRF(1:3)))
       WidthMesonMedium_GammaColl = get_inMediumWidth(ID,absP,mass,mediumAtPos)

    end select

  end function WidthMesonMedium_GammaColl

  !****************************************************************************
  !****s* mesonWidthMedium/GetMassAssInfo_Meson
  ! NAME
  ! subroutine GetMassAssInfo_Meson(MassAssInfo, ID,momLRF,mediumAtPos)
  ! PURPOSE
  ! Set the information necessary for improved MassAss
  ! INPUTS
  ! * integer :: ID -- ID of particle
  ! * real, dimension(0:3) :: momLRF -- momentum of resonance in LRF
  ! * type(medium) :: mediumAtPos -- medium at position
  ! OUTPUT
  ! * type(tMassAssInfo) :: MassAssInfo -- the info to get
  ! NOTES
  ! * The type 'MassAssInfo' is given as 'intent(inout)' since we have
  !   to bother with the allocation of its components
  !****************************************************************************
  subroutine GetMassAssInfo_Meson(MassAssInfo, ID,momLRF,mediumAtPos,forceMix)
    use IdTable, only: isMeson
    use constants, only: rhoNull
    use mediumDefinition
    use MassAssInfoDefinition
    use ParticleProperties, only: hadron
    use CALLSTACK, only: TRACEBACK

    type(tMassAssInfo), intent(inout) :: MassAssInfo
    integer,intent(in) :: ID
    real,intent(in),dimension(0:3)  :: momLRF
    type(medium),intent(in) :: mediumAtPos
    real,intent(in),OPTIONAL :: forceMix

    real :: dens
    integer :: iG
    real :: mixG

    if (.not.(isMeson(ID))) then
       write(*,*) 'no meson: ',ID
       call TRACEBACK()
    end if

    call ResetMassAssInfo(MassAssInfo)

    MassAssInfo%Mass0 = hadron(ID)%mass
    MassAssInfo%Gamma0 = WidthMesonMedium(ID,hadron(ID)%mass,momLRF,mediumAtPos)

    if (MassAssInfo%Gamma0 .lt. 1e-3) then
       call SetMassAssInfo_Stable(MassAssInfo)
       return
    end if

    ! Now we have to set all the bin information:

    MassAssInfo%IsStable = .false.

    iG = 0
    mixG = 0.

    select case (mediumSwitch)
    case (0,1)
       ! vacuum width or constant collisional width

       if (GammaColl(ID)>0.) then
          dens = min (mediumAtPos%density/rhoNull, rho_max)
          iG = int(dens*nG1)
          if (allowMix) then
             mixG = dens - real(iG)/nG1
             if (present(forceMix)) mixG = forceMix
          end if
       end if

       call SetMassAssInfo_FromArray(MassAssInfo,MassAssInfo_meson(ID),iG,mixG)

    case (2) ! full in-medium width
       call TRACEBACK('Not yet implemented')
    end select


  end subroutine GetMassAssInfo_Meson

  !****************************************************************************
  !****s* mesonWidthMedium/InitMassAssInfo_Meson
  ! NAME
  ! subroutine InitMassAssInfo_Meson
  ! This routine initializes the array of information used by the
  ! improved massass routines.
  ! NOTES
  ! * TO BE DONE: at the moment only implemented for vacuum width and
  !   plus constant collisional width. Latter not yet in full glory.
  !****************************************************************************
  subroutine InitMassAssInfo_Meson
    use IDtable, only: pion,nMes
    use CALLSTACK, only: TRACEBACK

    integer :: ID

    select case (mediumSwitch)
    case (0,1)
       ! vacuum width or constant collisional width
       do ID=pion,pion+nMes-1
         call Init_ConstGammaColl(ID)
       end do

    case (2)
       ! full in-medium width
       call TRACEBACK('Not yet implemented')

    end select

  contains

    subroutine Init_ConstGammaColl (ID)
      use mesonWidth, only: GetMinMaxMass, GetMaxQ
      use ParticleProperties, only: hadron
      use mesonWidth, only: FullWidthMeson

      integer, intent(in) :: ID

      integer, parameter :: nB = 10
      real :: MinMass, MaxMass
      real, dimension(nB) :: BinM,BinW,BinQ,BinY
      integer :: iB
      integer :: nB1,nB2

      integer :: nG,iG
      real :: mass0,gamma0,Gamma,dGamma,Gcoll

      mass0 = hadron(ID)%mass
      gamma0 = FullWidthMeson(ID,mass0)
      if (gamma0<1e-3) return
      Gcoll = GammaColl(ID)

      if (verboseInit) write(*,*) '-----------------------'
      if (verboseInit) write(*,*) ID,mass0,gamma0,Gcoll

      nG = 0
      if (Gcoll>0.) nG=nG2

      call GetMinMaxMass(ID,MinMass,MaxMass, (mediumSwitch.ne.0) )

      ! We allocate memory in the global array:
      allocate(MassAssInfo_meson(ID)%W(0:nG,nB-1))
      allocate(MassAssInfo_meson(ID)%Q(0:nG,nB-1))
      allocate(MassAssInfo_meson(ID)%M(0:nG,nB))
      allocate(MassAssInfo_meson(ID)%Y(0:nG,nB))
      allocate(MassAssInfo_meson(ID)%n1(0:nG))
      allocate(MassAssInfo_meson(ID)%n2(0:nG))

      do iG=0,nG
         if (verboseInit) write(*,*)
         if (verboseInit) write(*,*) 'iG=',iG

         dGamma = Gcoll * real(iG)/nG1
         Gamma = Gamma0 + dGamma

         BinM = (/MinMass,mass0-3*gamma,mass0-2.5*gamma,mass0-2*gamma,&
              & mass0-gamma,mass0,&
           & mass0+gamma, mass0+2*gamma, mass0+20*gamma, MaxMass/)

         do iB=1,nB
            if (BinM(iB).lt.MinMass) BinM(iB) = MinMass
            if (BinM(iB).gt.MaxMass) BinM(iB) = MaxMass
         end do

         if (verboseInit) write(*,'(A10,10f13.4)') 'BinM',BinM

         BinY = 2*atan2(2*(BinM-mass0),gamma)


         do iB=1,nB-1
            BinW(iB) = BinY(iB+1)-BinY(iB)
         end do

         ! find minimal/maximal bin used:
         do iB=1,nB-1
            nB1 = iB
            if (BinW(iB) > 0.0) exit
         end do
         do iB=nB-1,1,-1
            nB2 = iB
            if (BinW(iB) > 0.0) exit
         end do


         if (verboseInit) write(*,'(A10,7(" "),10f13.4)') &
              & 'BinW',BinW(1:nb-1)/sum(BinW(1:nb-1))

         call GetMaxQ(ID,mass0,gamma0,dGamma, BinM,BinQ)

         if (verboseInit) write(*,'(A10,7(" "),10f13.4)') 'BinQ',BinQ(1:nb-1)
         if (verboseInit) write(*,*) 'nBin: ',nB1,nB2

         ! Now we store the information in the global array:

         MassAssInfo_meson(ID)%W(iG,1:nB-1) = BinW(1:nB-1)
         MassAssInfo_meson(ID)%Q(iG,1:nB-1) = BinQ(1:nB-1)
         MassAssInfo_meson(ID)%M(iG,1:nB)   = BinM(1:nB)
         MassAssInfo_meson(ID)%Y(iG,1:nB)   = BinY(1:nB)
         MassAssInfo_meson(ID)%n1(iG)       = nB1
         MassAssInfo_meson(ID)%n2(iG)       = nB2

      end do


    end subroutine Init_ConstGammaColl

  end subroutine InitMassAssInfo_Meson


end module mesonWidthMedium
