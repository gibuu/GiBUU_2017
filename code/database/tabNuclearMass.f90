!******************************************************************************
!****m* /tabNuclearMass
! NAME
! module tabNuclearMass
! PURPOSE
! * Provide database structure for nuclides
!   (mass, binding energy, ground state spin, name)
! * intitialize nuclear decay properties (unused/untested in this context)
! * provide access to single data points without loading the entire 200kB data file to RAM
! * Read table "path_To_Input/nndc.dat", which contains the nuclide properties:
!   "Ame2003 atomic mass evaluation (II)"
! * G.Audi, A.H. Wapsta, C. Thibautl, NPA729(2003)337-676,
!   original file: mass.mas03
! * Based on Theo Gaitanos 'testRun/auswerteTools/ClusterCode/code/InitGEM.f90'
! AUTHORS
!   Theo Gaitanos, Ivan Lappo-Danilevski
!******************************************************************************
Module tabNuclearMass

  implicit none
  private

  !****************************************************************************
  !****t* tabNuclearMass/Nuclide
  ! NAME
  ! type Nuclide
  ! SOURCE
  !
  type Nuclide
     integer                            :: ZahlN        !neutron number
     integer                            :: ZahlZ        !proton number
     real                               :: MassEx       !mass excess
     real                               :: Mass         !mass in GeV
     real                               :: Bind         !binding energy
     real                               :: spin         !groundstate spin
     character(7)                       :: Name         !name of element
     integer, allocatable, dimension(:) :: Ejectile     !index-position of ejectile
     integer, allocatable, dimension(:) :: Daughter     !index-position of daughter
     logical, allocatable, dimension(:) :: ChannelStatus !status of decay channel
  End Type Nuclide
  !
  ! NOTES
  ! Constructor for nuclide properties
  !****************************************************************************

  !****************************************************************************
  !****g* tabNuclearMass/NucMax
  ! NAME
  ! integer, parameter :: NucMax
  ! PURPOSE
  ! Maximum number of nuclides considered in evaporation model,
  ! number of entries in 'nndc.dat'
  ! NOTES
  ! DO NOT CHANGE THIS NUMBER!!!
  !****************************************************************************
  integer, parameter :: NucMax=2467

  !****************************************************************************
  !****g* tabNuclearMass/Elements
  ! NAME
  ! type(Nuclide), Allocatable,dimension (:),save :: Elements
  ! PURPOSE
  ! Database for elements properties
  ! USAGE
  ! * To get e.g. the mass excess of the nuclide IN type "Elements(IN)%MassEx".
  !   The ID's are choosen according to NNDC-Table.
  ! * The array indizes of "Elements" run from NucMax=1n up to
  !   NucMax=208Pb (IN=1,...,NucMax).
  !****************************************************************************
  type(Nuclide), Allocatable, dimension(:), save :: Elements

  !****************************************************************************
  !****g* tabNuclearMass/NucPosition
  ! NAME
  ! Character, Allocatable,dimension(:,:),save :: NucPosition
  ! PURPOSE
  ! (N,Z)-Matrix providing with the position (in the nndc.dat-file) of a
  ! stable (according N/Z-ratio) nuclide. Negative Value of a matrix
  ! element indicates not existing nulcides.
  !
  !****************************************************************************
  integer, Allocatable, dimension(:,:), save :: NucPosition

  public :: Nuclide  !Nuclide type: NuclideProperties constructor
  public :: Elements !Nuclides with their properties
  public :: NucPosition !Position index (N,Z)-array of each element
  public :: Init_Mass_Table,SetDecays, read_element, read_mass !subroutines
!   public :: End_Mass_Table
  public :: NucMax

contains

  !****************************************************************************
  subroutine IOControl(loc,status,routine,variable)
  !****************************************************************************
    integer,      intent(in) :: status,loc
    character(*), intent(in) :: routine,variable
    if (status /= 0) then
       write(*,*) 'from subroutine ',routine,' : '
       if (loc==1) then !allocations
          write(*,*) 'allocation of variable: ',variable ,'NOT successfull'
       end if
       if (loc==2) then!deallocations
          write(*,*) 'deallocation of variable: ',variable ,'NOT successfull'
       end if
       write(*,*) '!!! TERMINATION OF PROGRAM NOW !!!'
       STOP
    end if
  !****************************************************************************
  end subroutine IOControl
  !****************************************************************************

    !**************************************************************************
    !****s* tabNuclearMass/read_element
    ! NAME
    ! subroutine read_element
    !
    ! PURPOSE
    ! * Reading of data-tables from NNDC (mass excess, spin)
    ! * Initialize the "Elements(:)%" vector
    ! * Initialize the "NucPosition(:,:)" array
    ! INPUTS
    ! * integer, intent(in)                  :: Ni -- Neutron number of requested element
    ! * integer, intent(in)                  :: Zi -- Proton number
    ! OUTPUT
    ! * type(Nuclide),intent(out)            :: element -- information about the nucleus
    ! NOTES
    ! returns an element with %name='missing', when element is not found in table
    !**************************************************************************
    subroutine read_element(Ni,Zi,element)
    !**************************************************************************
      use inputGeneral, only: path_to_input
      use constants, only: massu
      integer, intent(in)                  :: Ni
      integer, intent(in)                  :: Zi
      type(Nuclide),intent(out)            :: element

      character(len=6) :: Name
      integer          :: ios,ioe,i,N,Z!,j,MaxN,MaxZ,status
      real             :: MassEx,Bind,Spin,Bind_exp,Bind_woCoul
      !---------------------------------------------------------------------
      open(Unit=10,File=trim(path_to_input)//'/nndc.dat', &
           & Status='old', Action='read',Iostat=ios)
      if (ios /= 0) then
         write(*,*) 'from tabNuclearMass: NNDC-File ',trim(path_to_input)//' /nndc.dat','Open failed: ios = ',ios
         write(*,*) 'Are you sure you are using the latest version of "buuinput"?',ios
         write(*,*) 'STOP!'
         STOP
      end if
      !---------------------------------------------------------------------
      ! Search element in file
      !---------------------------------------------------------------------
      element%Name = 'missing'
      rewind(10)
      do i=1,45
         read(10,*) !skip comments
      end do
      do i=1,NucMax !-->Maximum 2467 elements!!!
         read(10,5) N,Z,name,MassEx,Bind_exp,Spin,Bind,Bind_woCoul
         Element%ZahlN  = N
         Element%ZahlZ  = Z
         if ( (N==Ni) .and. (Z==Zi) ) then
           Element%MassEx = MassEx*0.000001 !keV > GeV
           Element%Mass   = MassEx*0.000001 + massu*(Ni+Zi) !keV > GeV
           Element%Bind   = Bind*0.000001   !keV > GeV
           Element%Spin   = Spin
           Element%Name   = name
           exit
         end if
      end do
5     format(i4,i5,1x,a6,f14.5,2x,f12.5,6x,f7.4,f14.5,f14.5)
      close(unit=10,status='keep',iostat=ioe)
      if (ioe /= 0) then
         write(*,*) 'from tabNuclearMass: NNDC-File ',trim(path_to_input)//'/nndc.dat',' Close failed: ios = ',ioe
         write(*,*) 'Are you sure you are using the latest version of "buuinput"?',ioe
         write(*,*) 'STOP!'
      end if
    !**************************************************************************
    end subroutine read_element !*******************************************
    !**************************************************************************

  !****************************************************************************
  !****f* tabNuclearMass/read_mass
  ! NAME
  ! function read_mass
  ! INPUTS
  ! * integer, intent(in)              :: Z -- target charge
  ! * integer, intent(in)              :: N -- target neutron number
  ! OUTPUT
  ! * real                             :: read_mass
  ! PURPOSE
  ! reads the tabulated mass (in GeV) for a given N,Z
  !****************************************************************************

  function read_mass(N,Z)

    integer, intent(in)      :: N
    integer, intent(in)      :: Z
    real :: read_mass
    type(Nuclide) :: element

    call read_element(N,Z,element)
    if (element%name == 'missing') then
         write(*,*) 'from tabNuclearMass: no tabulated nuclide with N,Z =',N,Z
         write(*,*) 'from tabNuclearMass: !!! Termination of program NOW !!!'
         STOP
    end if
    read_mass = element%mass
  end function read_mass

    !**************************************************************************
    !****s* tabNuclearMass/init_Mass_Table
    ! NAME
    ! subroutine init_Mass_Table
    !
    ! PURPOSE
    ! * Reading of data-tables from NNDC (mass excess, spin)
    ! * Initialize the "Elements(:)%" vector
    ! * Initialize the "NucPosition(:,:)" array
    ! NOTES
    ! * Initialization of (1:MaxN,1:MaxZ)-Matrix "NucPosition" is done here.
    ! * Variables "MaxN" and "MaxZ" give the max. possible value of
    !   neutron- and proton-content of a cluster. Its status
    !   (stable/unstable) according N/Z-ratio is initialized here.
    ! * NucPosition is needed in checking the correct N/Z-ratio of
    !   produced clusters, before starting the evaporation procedure.
    !**************************************************************************
    subroutine Init_Mass_Table
    !**************************************************************************
      use inputGeneral, only: path_to_input
      use constants, only: massu
      integer          :: ios,ioe,status,i,j,N,Z,MaxN,MaxZ
      character(len=6) :: Element
      real             :: MassEx,Bind,Spin,Bind_exp,Bind_woCoul
      !---------------------------------------------------------------------
      open(Unit=10,File=trim(path_to_input)//'/nndc.dat', &
           & Status='old', Action='read',Iostat=ios)
      if (ios /= 0) then
         write(*,*) 'from tabNuclearMass: NNDC-File ',trim(path_to_input)//'/nndc.dat',' Open failed: ios = ',ios
         write(*,*) 'Are you sure you are using the latest version of "buuinput"?',ios
         write(*,*) 'STOP!'
      end if
      !---------------------------------------------------------------------
      allocate(Elements(1:NucMax),STAT=status)
      call IOControl(1,status,'Init_GEM','Elements')
      !---------------------------------------------------------------------
      ! Construct type(Nuclide) Elements
      !---------------------------------------------------------------------
      rewind(10)
      do i=1,45
         read(10,*) !skip comments
      end do
      do i=1,NucMax !-->Maximum 2467 elements!!!
         read(10,5) N,Z,Element,MassEx,Bind_exp,Spin,Bind,Bind_woCoul
         Elements(i)%ZahlN  = N
         Elements(i)%ZahlZ  = Z
         Elements(i)%MassEx = MassEx*0.000001 !keV > GeV
         Elements(i)%Mass   = MassEx*0.000001 + massu*(N+Z) !keV > GeV
         Elements(i)%Bind   = Bind*0.000001   !keV > GeV
         Elements(i)%Spin   = Spin
         Elements(i)%Name   = Element
      end do
5     format(i4,i5,1x,a6,f14.5,2x,f12.5,6x,f7.4,f14.5,f14.5)
      close(unit=10,status='keep',iostat=ioe)
      if (ioe /= 0) then
         write(*,*) 'from tabNuclearMass: cannot close NNDC-File, IOSTAT = ',ioe
         write(*,*) 'STOP!'
         STOP
      end if
      !---------------------------------------------------------------------
      MaxN = 130  !maxval(Elements%ZahlN)
      MaxZ = 90   !maxval(Elements%ZahlZ)
      allocate(NucPosition(0:MaxN,0:MaxZ),STAT=status)
      call IOControl(1,status,'Init_GEM','NucPosition')
      !---------------------------------------------------------------------
      !- Initialize (N,Z)-matrix NucPosition
      !- Define first the unstable elements (negative value)
      !- Then define the stable elements in the sequence of nndc.dat-table
      !---------------------------------------------------------------------
      do i=0,MaxN
         do j=0,MaxZ
            NucPosition(i,j) = -1 !unstable elements
         end do
      end do
      do i=1,NucMax
         NucPosition(Elements(i)%ZahlN,Elements(i)%ZahlZ)= i !stable elements
      end do
    !**************************************************************************
    end subroutine Init_Mass_table !****************************************
    !**************************************************************************


    !**************************************************************************
    !****s* tabNuclearMass/SetDecays
    ! NAME
    ! subroutine SetDecays
    !
    ! PURPOSE
    ! * Defines valid decay channels for each nuclide
    ! * Important for Monte-Carlo decision of decaying cluster
    ! * The list of elements is taken from the mass.mas03-table (cf. init_Mass_Table)
    ! NOTES
    ! Tested only with these lines of code in GiBUU
    !     call Init_Mass_Table
    !     do i = 1,NucMax
    !       write (*,*) Elements(i)%name
    !     end do
    !     call SetDecays
    !     call End_Mass_Table
    !     write (*,*) 'allocation/deallocation successfull'
    !     STOP
    ! Thus allocation and deallocation works. Fully tested in other context by Theo Gaitanos
    !**************************************************************************
    subroutine SetDecays
    !**************************************************************************
      integer              :: i,j,k,Ndiff,Zdiff,channel,Ai,Aj
      integer              :: tempN,tempZ,MaxVal,MaxDaughter
      integer              :: EjectilePosition,DaughterPosition
      logical              :: flag
      !MaxEmission gives the max. number of decays. In original GEM-model
      !only up to 4He-dacays have been considered. Here we consider many
      !more...The choisc of MaxEmission just depends on the CPU-time consumed.
      !Its also possible to reduce it.
      integer, parameter   :: MaxEmission=100
      !---------------------------------------------------------------------
      ! neutron does not decay...
      allocate(Elements(1)%Ejectile(1:1))
      allocate(Elements(1)%Daughter(1:1))
      allocate(Elements(1)%ChannelStatus(1:1))
      Elements(1)%Ejectile(1)      = 1
      Elements(1)%Daughter(1)      = 1
      Elements(1)%ChannelStatus(1) = .false.
      ! proton does not decay...
      allocate(Elements(2)%Ejectile(1:1))
      allocate(Elements(2)%Daughter(1:1))
      allocate(Elements(2)%ChannelStatus(1:1))
      Elements(2)%Ejectile(1)      = 2
      Elements(2)%Daughter(1)      = 2
      Elements(2)%ChannelStatus(1) = .false.

      !---------------------------------------------------------------------
      ! * set the decay channels for each Element
      ! * Emissions up to NucMax=MaxEmission are considered, otherwise
      !   extremely time consuming...
      !---------------------------------------------------------------------
      NuclideLoop : do i=1,NucMax
         Ai = Elements(i)%ZahlN + Elements(i)%ZahlZ
         if (i.le.2) cycle !skeep n,p-->see also sequence in nndc.dat file
         channel = 1
         MaxDaughter  = i-1
         MaxVal       = i-1
         if (MaxVal > MaxEmission) MaxVal = MaxEmission
         allocate(Elements(i)%Ejectile(1:MaxVal))
         Elements(i)%Ejectile = 0 !--> set it undefined at init
         allocate(Elements(i)%Daughter(1:MaxVal))
         Elements(i)%Daughter = 0 !--> set it undefined at init
         allocate(Elements(i)%ChannelStatus(1:MaxVal))
         Elements(i)%ChannelStatus = .false. !--> set it undefined at init
         DecayChannels : do j=1,MaxVal
            Aj = Elements(j)%ZahlN + Elements(j)%ZahlZ
            if (Aj == Ai) cycle
            EjectilePosition = NucPosition(Elements(j)%ZahlN,Elements(j)%ZahlZ)
            Ndiff            = Elements(i)%ZahlN-Elements(j)%ZahlN
            Zdiff            = Elements(i)%ZahlZ-Elements(j)%ZahlZ
            flag = .false.
            Find1 : do k=1,MaxDaughter
               if (Elements(k)%ZahlN==Ndiff .and. Elements(k)%ZahlZ==Zdiff) then
                  flag  = .true.
                  tempN = Elements(k)%ZahlN
                  tempZ = Elements(k)%ZahlZ
                  exit
               end if
            end do Find1
            if (.not.flag) cycle
            DaughterPosition = NucPosition(tempN,tempZ)
            Elements(i)%Ejectile(channel)      = EjectilePosition
            Elements(i)%Daughter(channel)      = DaughterPosition
            Elements(i)%ChannelStatus(channel) = .true.
            channel = channel + 1
         end do DecayChannels
      end do NuclideLoop
    !**************************************************************************
    end subroutine SetDecays
    !**************************************************************************


!     !***********************************************************************
!     subroutine End_Mass_Table
!     !***********************************************************************
!       integer :: status,i
!
!       do i=1,NucMax
!          Deallocate(Elements(i)%Ejectile,STAT=status)
!          if(status /= 0) then
!             write(*,*) 'deallocation of variable "%Ejectile" NOT successfull'
!             write(*,*) '!!! TERMINATION NOW...!!!'
!             STOP
!          endif
!          Deallocate(Elements(i)%Daughter,STAT=status)
!          if(status /= 0) then
!             write(*,*) 'deallocation of variable "%Daughter" NOT successfull'
!             write(*,*) '!!! TERMINATION NOW...!!!'
!             STOP
!          endif
!       end do
!
!       Deallocate(Elements,STAT=status)
!       if(status /= 0) then
!          write(*,*) 'deallocation of variable "Elements" NOT successfull'
!          write(*,*) '!!! TERMINATION NOW...!!!'
!          STOP
!       endif
!
!       Deallocate(NucPosition,STAT=status)
!       if(status /= 0) then
!          write(*,*) 'deallocation of variable "NucPosition" NOT successfull'
!          write(*,*) '!!! TERMINATION NOW...!!!'
!          STOP
!       endif
!
!     !***********************************************************************
!     end subroutine End_Mass_Table !*****************************************
!     !***********************************************************************


  end Module tabNuclearMass
