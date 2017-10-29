!******************************************************************************
!****m* /baryonWidthMedium
! NAME
! module baryonWidthMedium
! PURPOSE
! Basically this module calls the routines of baryonWidth and adds medium
! modifications to it.
! USES
! module baryonWidth
!******************************************************************************
module baryonWidthMedium
  use MassAssInfoDefinition, only: tMassAssInfoArray
  use IDtable, only: nbar

  implicit none
  private

  !****************************************************************************
  !****g* baryonWidthMedium/mediumSwitch
  ! SOURCE
  !
  logical, save :: mediumSwitch=.false.
  ! PURPOSE
  ! Switch on and off the in-medium width of all baryons at once.
  ! If .false., the vacuum width are used.
  !****************************************************************************

  !****************************************************************************
  !****g* baryonWidthMedium/mediumSwitch_Delta
  ! SOURCE
  !
  logical, save :: mediumSwitch_Delta=.false.
  ! PURPOSE
  ! Only meaningful if mediumSwitch=.true.:
  ! Switch on and off the in-medium width of the Delta. (.false.=off)
  !
  ! Note that in that case the Delta is treated specially: what is used
  ! for the in-medium width is determined by the flag in deltaWidth.
  ! This switch is not consistent with mediumSwitch_coll!
  !****************************************************************************

  !****************************************************************************
  !****g* baryonWidthMedium/mediumSwitch_coll
  ! SOURCE
  !
  logical, save :: mediumSwitch_coll=.false.
  ! PURPOSE
  ! Only meaningful if mediumSwitch=.true.:
  ! Use in-medium width according to collision term.
  ! NOTES
  ! if set to TRUE, then also UseOffShellPotentialBaryons (see module offShellPotential) must be .true.
  !****************************************************************************

  !****************************************************************************
  !****g* baryonWidthMedium/mediumSwitch_proton_neutron
  ! SOURCE
  !
  logical, save :: mediumSwitch_proton_neutron=.false.
  ! PURPOSE
  ! Only meaningful if mediumSwitch=.true.:
  ! Switch on and off the in-medium width of the proton and the neutron.
  ! (.false.=off)
  !
  ! Note that in that case the nucleons are treated specially.
  ! This switch is not consistent with mediumSwitch_coll!
  !****************************************************************************

  !****************************************************************************
  !****g* baryonWidthMedium/verboseInit
  ! SOURCE
  !
  logical, save :: verboseInit = .false.
  ! PURPOSE
  ! switch on/off informational messages during initialization
  !****************************************************************************

  !****************************************************************************
  !****g* baryonWidthMedium/verboseInitStop
  ! SOURCE
  !
  logical, save :: verboseInitStop = .false.
  ! PURPOSE
  ! Stop after informational messages during initialization or not.
  !****************************************************************************

  logical, save :: initFlag=.true.

  type(tMassAssInfoArray), dimension(1:nbar), save :: MassAssInfo_baryon

  public :: WidthBaryonMedium
  public :: partialWidthBaryonMedium
  public :: decayWidthBaryonMedium
  public :: get_MediumSwitch
  public :: get_MediumSwitch_coll
  public :: get_MediumSwitch_Delta
  public :: get_MediumSwitch_proton_neutron
  public :: GetMassAssInfo_Baryon

contains

  !****************************************************************************
  !****f* baryonWidthMedium/get_MediumSwitch
  ! NAME
  ! logical function get_MediumSwitch()
  ! PURPOSE
  ! Returns the value of mediumSwitch
  !****************************************************************************
  logical function get_MediumSwitch()
    if (initFlag) call readInput
    get_MediumSwitch=mediumSwitch
  end function get_MediumSwitch


  !****************************************************************************
  !****f* baryonWidthMedium/get_MediumSwitch_coll
  ! NAME
  ! logical function get_MediumSwitch_coll()
  ! PURPOSE
  ! Returns the value of mediumSwitch_coll
  !****************************************************************************
  logical function get_MediumSwitch_coll()
    if (initFlag) call readInput
    get_MediumSwitch_coll=mediumSwitch_coll
  end function get_MediumSwitch_coll


  !****************************************************************************
  !****f* baryonWidthMedium/get_MediumSwitch_Delta
  ! NAME
  ! logical function get_MediumSwitch_Delta()
  ! PURPOSE
  ! Returns the value of mediumSwitch_Delta
  !****************************************************************************
  logical function get_MediumSwitch_Delta()
    if (initFlag) call readInput
    get_MediumSwitch_Delta=mediumSwitch_Delta
  end function get_MediumSwitch_Delta


  !****************************************************************************
  !****f* baryonWidthMedium/get_MediumSwitch_proton_neutron
  ! NAME
  ! logical function get_MediumSwitch_proton_neutron()
  ! PURPOSE
  ! Returns the value of mediumSwitch_proton_neutron
  !****************************************************************************
  logical function get_MediumSwitch_proton_neutron()
    if (initFlag) call readInput
    get_MediumSwitch_proton_neutron=mediumSwitch_proton_neutron
  end function get_MediumSwitch_proton_neutron


  !****************************************************************************
  !****s* baryonWidthMedium/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "width_baryon".
  !****************************************************************************
  subroutine readInput
    use baryonWidthMedium_tables, only: get_deltaOset
    use output, only: Write_ReadingInput
    use MassAssInfoDefinition, only: Set_UseMassAssInfo

    integer :: ios ! checks file behavior

    !**************************************************************************
    !****n* baryonWidthMedium/width_Baryon
    ! NAME
    ! NAMELIST width_Baryon
    ! PURPOSE
    ! Includes the switches:
    ! * mediumSwitch
    ! * mediumSwitch_Delta
    ! * mediumSwitch_proton_neutron
    ! * mediumSwitch_coll
    ! * verboseInit
    ! * verboseInitStop
    !**************************************************************************
    NAMELIST /width_Baryon/ mediumSwitch, mediumSwitch_Delta, &
                            mediumSwitch_proton_neutron, mediumSwitch_coll, &
                            verboseInit, verboseInitStop

    call Write_ReadingInput('width_Baryon',0)
    rewind(5)
    read(5,nml=width_Baryon,IOSTAT=IOS)
    call Write_ReadingInput('width_Baryon',0,IOS)

    write(*,'(A,L8)') ' Use in medium width for the baryons         : ',mediumSwitch
    if (.not.mediumSwitch) then
       mediumSwitch_coll=.false.
       mediumSwitch_Delta=.false.
       mediumSwitch_proton_neutron=.false.
       write(*,'(A)') '   -> No in medium collisional width for the baryons'
       write(*,'(A)') '   -> No special in medium width for the Delta'
       write(*,'(A)') '   -> No special in medium width for the nucleons'
    end if

    write(*,'(A,L8)') ' Use medium width according to collision term: ',&
         & mediumSwitch_coll
    if (mediumSwitch_coll) then
       mediumSwitch_delta=get_deltaOset()
       mediumSwitch_proton_neutron=.false.
       write(*,'(A,L8)') ' Use special Delta width                     : ',&
            & mediumSwitch_Delta
       write(*,'(A,L8)') ' Use in medium width for the nucleons        : ',&
            & mediumSwitch_proton_neutron
    else
       write(*,'(A,L8)') ' Use in medium width for the Delta           : ',&
            & mediumSwitch_Delta
       write(*,'(A,L8)') ' Use in medium width for the nucleons        : ',&
            & mediumSwitch_proton_neutron
    end if

    call Write_ReadingInput('width_Baryon',1)

    ! Ensure we use the new MassAss only if implemented:
    if (mediumswitch) then
       if (mediumSwitch_coll) then
          call Set_UseMassAssInfo(.false.)
       else
          if (mediumSwitch_proton_neutron) call Set_UseMassAssInfo(.false.)
          if (mediumSwitch_Delta) then
!             write(*,*) 'we are working on it... You should set it to false.'
!             call Set_UseMassAssInfo(.false.)
          end if
       end if
    end if
    call InitMassAssInfo_Baryon()

    initFlag=.false.

  end subroutine readInput


  !****************************************************************************
  !****f* baryonWidthMedium/partialWidthBaryonMedium
  ! NAME
  ! real function partialWidthBaryonMedium(particleID,mass,inWidth,mesonID,baryonID,momentumLRF,mediumATposition,baryonMass,mesonMass)
  !
  ! PURPOSE
  ! This function calculates the partial width with energy dependence
  ! and medium modifications.
  ! If the medium modification is not yet implemented for a specific case,
  ! partialWidthBaryon (vacuum) is called!
  !
  ! INPUTS
  ! * integer   :: particleID    -- id of resonance
  ! * real      :: mass          -- bare mass of the resonance (offshell),
  !   not including the potentials!
  ! * logical   :: inWidth       -- .true.=> in-width (only important for channels
  !   with unstable particles); .false.=> out-width
  ! * integer  :: mesonID, baryonID -- ID's of the decay products which one is
  !   interested in
  ! * type(medium) :: mediumATposition  -- Medium information
  ! * real, dimension(0:3) :: momentumLRF -- Momentum im LRF of the resonance
  ! * real, OPTIONAL     :: baryonMass,mesonMass --
  !   Possibility to define the masses of the incoming baryon and meson,
  !   needed in the case of the In-Width if one of them is off-shell.
  !   Otherwise not relevant.
  !****************************************************************************
  real function partialWidthBaryonMedium(particleID,mass,inWidth,mesonID,baryonID,momLRF,mediumATposition,baryonMass,mesonMass)

    use mediumDefinition
    use deltaWidth, only: delta_nucleonPion
    use idTable, only: delta, nucleon, pion, isBaryon, isMeson
    use CALLSTACK, only: TRACEBACK
    use baryonWidth, only: partialWidthBaryon

    real,intent(in) :: mass
    integer,intent(in) :: particleID
    logical,intent(in) ::  inWidth
    integer, intent(in) :: mesonID, baryonID

    ! ONLY Needed for InMediumModifications:
    type(medium),intent(in) ::  mediumATposition
    real,intent(in),dimension(0:3) :: momLRF

    real, intent(in), optional :: baryonMass,mesonMass

    real :: dens     ! density neutron+proton

    !real :: rhop,rhoN,absP

    if (initFlag) call readInput

    ! (1) Check Input
    if (.not.isBaryon(baryonID)) then
       write(*,*) 'Wrong baryonID:', baryonID
       call TRACEBACK()
    end if

    if (.not.isMeson(mesonID)) then
       write(*,*) 'Wrong mesonID:', mesonID
       call TRACEBACK()
    end if

    if (.not.isBaryon(particleID)) then
       write(*,*) 'Wrong particleID:', particleID
       call TRACEBACK()
    end if

    if (mediumAtPosition%useMedium.and.mediumSwitch) then
       if (mediumSwitch_coll) then
          ! No special medium modifications!!
          if (present(baryonMass).and.Present(mesonMass)) then
             partialwidthBaryonMedium=partialwidthBaryon(particleID,mass,inWidth,mesonID,baryonID,mesonMass,baryonMass)
          else
             partialwidthBaryonMedium=partialwidthBaryon(particleID,mass,inWidth,mesonID,baryonID)
          end if
       else
          select case (particleID)

          case (Delta)
             ! modified P33(1232)
             if (mediumSwitch_Delta) then
                ! assuming onshell pions and nucleons
                if ( (mesonID.eq.pion).and.(baryonId.eq.nucleon) ) then
                   dens=mediumAtPosition%densityProton+mediumAtPosition%densityNeutron
                   partialWidthBaryonMedium=delta_nucleonPion(mass,momLRF(1:3),dens)
                else
                   partialWidthBaryonMedium=0.
                end if
             else
                if (present(baryonMass).and.Present(mesonMass)) then
                   partialwidthBaryonMedium=partialwidthBaryon(particleID,mass,inWidth,mesonID,baryonID,mesonMass,baryonMass)
                else
                   partialwidthBaryonMedium=partialwidthBaryon(particleID,mass,inWidth,mesonID,baryonID)
                end if
             end if

          case default
             ! Since no medium modifications yet implemented :
             if (present(baryonMass).and.Present(mesonMass)) then
                partialwidthBaryonMedium=partialwidthBaryon(particleID,mass,inWidth,mesonID,baryonID,mesonMass,baryonMass)
             else
                partialwidthBaryonMedium=partialwidthBaryon(particleID,mass,inWidth,mesonID,baryonID)
             end if
          end select
       end if
    else

       ! No medium modifications
       if (present(baryonMass).and.Present(mesonMass)) then
          partialwidthBaryonMedium=partialwidthBaryon(particleID,mass,inWidth,mesonID,baryonID,mesonMass,baryonMass)
       else
          partialwidthBaryonMedium=partialwidthBaryon(particleID,mass,inWidth,mesonID,baryonID)
       end if
    end if

  end function partialwidthBaryonMedium


  !****************************************************************************
  !****s* baryonWidthMedium/decayWidthBaryonMedium
  ! NAME
  ! subroutine decayWidthBaryonMedium(particleID,mass,momentumLRF,mediumATposition, decayWidth, pauliFlag)
  !
  ! PURPOSE
  ! This function returns the partial out-widths with energy dependence and
  ! medium modifications.
  ! If the out-width has already Pauli blocking effects included, then
  ! pauliFlag is set to .true. .
  !
  ! INPUTS
  ! * integer            :: particleID        -- id of resonance
  ! * real               :: mass              --
  !   bare mass of the resonance (offshell), not including potentials!
  ! * type(medium)       :: mediumATposition  -- Medium information
  ! * real,dimension(0:3):: momentumLRF       --
  !   Momentum im LRF of the resonance
  !
  ! OUTPUT
  ! * real, dimension(:)             :: decayWidth --
  !   array of decay widths for all decay channels
  ! * logical, intent(out)           :: pauliFlag  --
  !   If .true. then pauli-Blocking is already considered
  !   in the decay width description
  !****************************************************************************
  subroutine decayWidthBaryonMedium(particleID,mass,momLRF,mediumATposition, decayWidth, pauliFlag)

    use mediumDefinition
    use particleProperties, only: nDecays
    use baryonWidth, only: decayWidthBaryon

    real,intent(in) :: mass
    integer,intent(in) :: particleID
    real, intent(out), dimension(1:nDecays)  :: decayWidth
    logical, intent(out) :: pauliFlag

    ! ONLY Needed for InMediumModifications:
    type(medium),intent(in) ::  mediumATposition
    real,intent(in),dimension(0:3) :: momLRF

    if (initFlag) call readInput

    ! No medium modification
    pauliFlag=.false.
    decayWidth = decayWidthBaryon (particleID, mass)

  end subroutine decayWidthBaryonMedium


  !****************************************************************************
  !****f* baryonWidthMedium/WidthBaryonMedium
  ! NAME
  ! real function WidthBaryonMedium(particleID,mass,momentumLRF,mediumATposition,outOfBounds)
  !
  ! PURPOSE
  ! This function calculates the full width with energy dependence and
  ! medium modifications.
  ! If for a specific resonance, medium modifications are not yet implemented,
  ! then WidthBaryon is called and the vacuum width is used.
  !
  ! INPUTS
  ! * integer            :: particleID        -- id of resonance
  ! * real               :: mass              -- bare mass of the resonance
  !   (offshell), not including potentials
  ! * type(medium)       :: mediumATposition  -- Medium at decay vertex
  ! * real,dimension(0:3):: momentumLRF       -- momentum of resonance in LRF
  !
  ! OUTPUT
  ! * logical, OPTIONAL :: outOfBounds
  !
  ! NOTES
  ! For the Delta, the in-medium width 'delta_fullWidth' does not
  ! converge to the vacuum values 'FullWidthBaryon(2)' for large M
  ! in the limit rho->0 ! This is due to upper bound of the tabulation
  ! in both cases (for the inMedium cases, the width is kept constant
  ! for M>2GeV, while inVacuum the upper mass is larger.)
  !
  ! In order to get a smooth behaviour, we do not respect "medium%useMedium"
  ! anymore, if "mediumSwitch" is set: Also very tiny densities are treated
  ! as 'in-medium' and not as vacuum.
  !****************************************************************************
  real function WidthBaryonMedium(particleID,mass,momLRF,mediumATposition,outOfBounds)

    use mediumDefinition
    use idTable, only: delta,nucleon,nres,isBaryon
    use deltaWidth, only: delta_fullWidth
    use pn_medium_width, only: proton_width_medium
    use baryonWidth, only: FullWidthBaryon
    use baryonWidthMedium_tables, only: get_inMediumWidth
    use CallStack, only: TRACEBACK

    real,intent(in) ::  mass
    integer,intent(in) :: particleID
    type(medium),intent(in) ::  mediumATposition
    real,intent(in),dimension(0:3) :: momLRF
    logical, optional,intent(out) :: outOfBounds

    real :: rhoP,rhoN
    logical :: lBounds

    if (present(outOfBounds)) outOfBounds=.false.

    if (initFlag) call readInput

    ! Check Input
    if (.not.isBaryon(particleID)) then
       write(*,*) 'Wrong particleID:', particleID
       call TRACEBACK()
    end if

    ! Set the default value (as for vacuum):
    WidthBaryonMedium= FullWidthBaryon(particleID,mass)

    ! Now do the medium modifications:

    if (mediumSwitch) then

       if (mediumAtPosition%useMedium) then
          rhoP = mediumAtPosition%densityProton
          rhoN = mediumAtPosition%densityNeutron

          if (mediumSwitch_coll.and.particleID.le.(nres+1)) then
             !only implemented for non-strange resonances
             WidthBaryonMedium=get_inMediumWidth(particleID,momLRF,mass,&
                  & rhoN, rhoP, 3,outOfBounds=lBounds)
             if (present(outOfBounds)) outOfBounds=lBounds
          else
             select case (particleID)
             case (Delta)
                if (mediumSwitch_Delta) then
                   WidthBaryonMedium= delta_fullWidth(mass,momLRF(1:3),rhoP+rhoN)
                end if

             case (nucleon)
                if (mediumSwitch_proton_neutron) then
                   WidthBaryonMedium=proton_width_medium(momLRF,rhoP,rhoN)
                end if
             end select
          end if

       else

          if (mediumSwitch_coll.and.particleID.le.(nres+1)) then
             ! do nothing
          else
             select case (particleID)
             case (Delta)
                if (mediumSwitch_Delta) then
                   WidthBaryonMedium= delta_fullWidth(mass,momLRF(1:3),0.)
                end if
             end select
          end if
       end if

    end if

  end function WidthBaryonMedium


  !****************************************************************************
  !****s* baryonWidthMedium/GetMassAssInfo_Baryon
  ! NAME
  ! subroutine GetMassAssInfo_Baryon(MassAssInfo, ID,momLRF,mediumAtPos)
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
  subroutine GetMassAssInfo_Baryon(MassAssInfo, ID,momLRF,mediumAtPos)
    use MassAssInfoDefinition
    use mediumDefinition
    use IdTable, only: isBaryon
    use ParticleProperties, only: hadron
    use CALLSTACK, only: TRACEBACK
    use deltaWidth, only: GetRhoBin,delta_fullWidth
    use baryonWidth, only: FullWidthBaryon

    type(tMassAssInfo), intent(inout) :: MassAssInfo
    integer,intent(in) :: ID
    real,intent(in),dimension(0:3)  :: momLRF
    type(medium),intent(in) :: mediumAtPos

    integer :: iG
    real :: mixG

    if (.not.(isBaryon(ID))) then
       write(*,*) 'no baryon: ',ID
       call TRACEBACK()
    end if

    if (initFlag) call readInput

    call ResetMassAssInfo(MassAssInfo)
    MassAssInfo%Mass0 = hadron(ID)%mass
    MassAssInfo%Gamma0 = FullWidthBaryon(ID,hadron(ID)%mass)

    if (MassAssInfo%Gamma0 .lt. 1e-3) then
       call SetMassAssInfo_Stable(MassAssInfo)
       return
    end if

    ! Now we have to set all the bin information:

    iG = 0
    mixG = 0.0

    if (.not.mediumSwitch) then ! vacuum width

       call SetMassAssInfo_FromArray(MassAssInfo,MassAssInfo_baryon(ID),iG,mixG)

    else

       if (mediumSwitch_coll) then
          call TRACEBACK('Not yet implemented')
       else
          select case (ID)
          case (1)
             if (mediumSwitch_proton_neutron) then
                call TRACEBACK('Not yet implemented')
             end if

          case (2)
             if (mediumSwitch_Delta) then
                call GetRhoBin(mediumAtPos%density,iG)
                MassAssInfo%IsMomDep = .TRUE.
                MassAssInfo%Gamma0 = delta_fullWidth(MassAssInfo%Mass0,&
                     & (/0.0,0.0,0.0/),mediumAtPos%density)
             end if
          end select
       end if

       call SetMassAssInfo_FromArray(MassAssInfo,MassAssInfo_baryon(ID),iG,mixG)

    end if

  end subroutine GetMassAssInfo_Baryon


  !****************************************************************************
  !****s* baryonWidthMedium/InitMassAssInfo_Baryon
  ! NAME
  ! subroutine InitMassAssInfo_Baryon
  ! PURPOSE
  ! This routine initializes the array of information used by the
  ! improved massass routines.
  !****************************************************************************
  subroutine InitMassAssInfo_Baryon
    use ParticleProperties, only: hadron
    use CALLSTACK, only: TRACEBACK
    use deltaWidth, only: delta_fullWidth, GetRhoValue, GetMaxQ_Delta
    use baryonWidth, only: FullWidthBaryon, GetMaxQ
    use MassAssInfoDefinition, only: Get_UseMassAssInfo

    integer :: ID
    real :: mass0, gamma0

    if (.not.Get_UseMassAssInfo()) return

    if (.not.mediumSwitch) then
       do ID = 1,nbar
          mass0 = hadron(ID)%mass
          gamma0 = FullWidthBaryon(ID,mass0)
          if (gamma0.ge. 1e-3) call Init_Vacuum
       end do

    else

       if (mediumSwitch_coll) then
          call TRACEBACK('Not yet implemented')
       else
          ! === nucleon
          ID = 1
          if (mediumSwitch_proton_neutron) then
             call TRACEBACK('Not yet implemented')
          else
             mass0 = hadron(ID)%mass
             gamma0 = FullWidthBaryon(ID,mass0)
             if (gamma0.ge. 1e-3) call Init_Vacuum
          end if

          ! === Delta
          ID = 2
          if (mediumSwitch_Delta) then
             ! depends on density and momentum !!!
             call Init_Delta()
          else
             mass0 = hadron(ID)%mass
             gamma0 = FullWidthBaryon(ID,mass0)
             if (gamma0.ge. 1e-3) call Init_Vacuum
          end if

          ! === all the others
          do ID = 3,nbar
             mass0 = hadron(ID)%mass
             gamma0 = FullWidthBaryon(ID,mass0)
             if (gamma0.ge. 1e-3) call Init_Vacuum
          end do

       end if

    end if

    if (verboseInitStop) call TRACEBACK('VerboseInit: STOP NOW')


  contains

    subroutine Init_Vacuum
      integer, parameter :: nB = 11
      real :: MinMass, MaxMass
      real, dimension(nB) :: BinM,BinW,BinQ,BinY
      integer :: iB
      integer :: nB1,nB2

      if (verboseInit) write(*,*) '-----------------------'
      if (verboseInit) write(*,*) ID,mass0,gamma0

      MinMass = hadron(ID)%minmass
      MaxMass = 3.0

      ! We allocate memory in the global array:
      allocate(MassAssInfo_baryon(ID)%W(0:0,nB-1))
      allocate(MassAssInfo_baryon(ID)%Q(0:0,nB-1))
      allocate(MassAssInfo_baryon(ID)%M(0:0,nB))
      allocate(MassAssInfo_baryon(ID)%Y(0:0,nB))
      allocate(MassAssInfo_baryon(ID)%n1(0:0))
      allocate(MassAssInfo_baryon(ID)%n2(0:0))

      ! please note: Resonances 34,36,38 may get a very large
      ! MaxBWD in the order 40..80 for a maxmass~3GeV

      BinM = (/MinMass,mass0-2*gamma0,mass0-gamma0,mass0-0.5*gamma0,mass0,&
           & mass0+gamma0, mass0+2*gamma0, mass0+5*gamma0, mass0+10*gamma0, &
           & mass0+20*gamma0,MaxMass/)

      do iB=1,nB
         if (BinM(iB).lt.MinMass) BinM(iB) = MinMass
         if (BinM(iB).gt.MaxMass) BinM(iB) = MaxMass
      end do

      if (verboseInit) write(*,'(A10,20f13.4)') 'BinM',BinM

      BinY = 2*atan2(2*(BinM-mass0),gamma0)


      do iB=1,nB-1
         BinW(iB) = BinY(iB+1)-BinY(iB)
      end do
      do iB=1,nB-1
         nB1 = iB
         if (BinW(iB) > 0.0) exit
      end do
      do iB=nB-1,1,-1
         nB2 = iB
         if (BinW(iB) > 0.0) exit
      end do


      if (verboseInit) write(*,'(A10,7(" "),20f13.4)') 'BinW',BinW(1:nb-1)/sum(BinW(1:nb-1))

      call GetMaxQ(ID,mass0,gamma0,0.0, BinM,BinQ)

      if (verboseInit) write(*,'(A10,7(" "),20f13.4)') 'BinQ',BinQ(1:nb-1)
      if (verboseInit) write(*,*) 'nBin: ',nB1,nB2

      ! Now we store the information in the global array:

      MassAssInfo_baryon(ID)%W(0,1:nB-1) = BinW(1:nB-1)
      MassAssInfo_baryon(ID)%Q(0,1:nB-1) = BinQ(1:nB-1)
      MassAssInfo_baryon(ID)%M(0,1:nB)   = BinM(1:nB)
      MassAssInfo_baryon(ID)%Y(0,1:nB)   = BinY(1:nB)
      MassAssInfo_baryon(ID)%n1(0)       = nB1
      MassAssInfo_baryon(ID)%n2(0)       = nB2


    end subroutine Init_Vacuum

    subroutine Init_Delta
      use constants, only: rhoNull

      integer, parameter :: nB = 11
      real :: MinMass, MaxMass
      real, dimension(nB) :: BinM,BinW,BinQ,BinY
      integer :: iB
      integer :: nB1,nB2

      real    :: rhomin,drho,rho
      integer :: irho,nrho

      call GetRhoValue(rhomin,drho,nrho)

      mass0 = hadron(2)%mass
      gamma0= delta_fullWidth(mass0,(/0.0,0.0,0.0/),0.0)

      if (verboseInit) write(*,*) '-----------------------'
      if (verboseInit) write(*,*) ID,mass0,gamma0

      MinMass = hadron(ID)%minmass
      MaxMass = 3.0

      ! We allocate memory in the global array:
      allocate(MassAssInfo_baryon(ID)%W(0:nrho,nB-1))
      allocate(MassAssInfo_baryon(ID)%Q(0:nrho,nB-1))
      allocate(MassAssInfo_baryon(ID)%M(0:nrho,nB))
      allocate(MassAssInfo_baryon(ID)%Y(0:nrho,nB))
      allocate(MassAssInfo_baryon(ID)%n1(0:nrho))
      allocate(MassAssInfo_baryon(ID)%n2(0:nrho))

      do irho=0,nrho
         rho = rhoMin+drho*rhoNull*iRho
         gamma0= delta_fullWidth(mass0,(/0.0,0.0,0.0/),rho)

         if (verboseInit) write(*,*) '- - - - - - - - - - - -'
         if (verboseInit) write(*,*) 'rho',irho,rho,gamma0

         BinM = (/MinMass,mass0-2*gamma0,mass0-gamma0,mass0-0.5*gamma0,mass0,&
              & mass0+gamma0, mass0+2*gamma0, mass0+3*gamma0, mass0+5*gamma0, &
              & mass0+10*gamma0,MaxMass/)

         do iB=1,nB
            if (BinM(iB).lt.MinMass) BinM(iB) = MinMass
            if (BinM(iB).gt.MaxMass) BinM(iB) = MaxMass
         end do

         if (verboseInit) write(*,'(A10,20f13.4)') 'BinM',BinM


         BinY = 2*atan2(2*(BinM-mass0),gamma0)


         do iB=1,nB-1
            BinW(iB) = BinY(iB+1)-BinY(iB)
         end do
         do iB=1,nB-1
            nB1 = iB
            if (BinW(iB) > 0.0) exit
         end do
         do iB=nB-1,1,-1
            nB2 = iB
            if (BinW(iB) > 0.0) exit
         end do


         if (verboseInit) write(*,'(A10,7(" "),20f13.4)') 'BinW',BinW(1:nb-1)/sum(BinW(1:nb-1))

         call GetMaxQ_Delta(mass0,gamma0,rho, BinM,BinQ)

         if (verboseInit) write(*,'(A10,7(" "),20f13.4)') 'BinQ',BinQ(1:nb-1)
         if (verboseInit) write(*,*) 'nBin: ',nB1,nB2

         ! Now we store the information in the global array:

         MassAssInfo_baryon(ID)%W(irho,1:nB-1) = BinW(1:nB-1)
         MassAssInfo_baryon(ID)%Q(irho,1:nB-1) = BinQ(1:nB-1)
         MassAssInfo_baryon(ID)%M(irho,1:nB)   = BinM(1:nB)
         MassAssInfo_baryon(ID)%Y(irho,1:nB)   = BinY(1:nB)
         MassAssInfo_baryon(ID)%n1(irho)       = nB1
         MassAssInfo_baryon(ID)%n2(irho)       = nB2

      end do

    end subroutine Init_Delta

  end subroutine InitMassAssInfo_Baryon

end module baryonWidthMedium
