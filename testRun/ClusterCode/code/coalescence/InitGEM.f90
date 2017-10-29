!***************************************************************************
!****m* /initGEM
! NAME
! module initGEM
! FUNCTION
! * Initialize database structure for nuclides and their decay properties
! * Read table "PathToClusterInput/nndc.dat", which contains the nuclide properties:
!   "Ame2003 atomic mass evaluation (II)" 
! * G.Audi, A.H. Wapsta, C. Thibautl, NPA729(2003)337-676, 
!   original file: mass.mas03
!***************************************************************************
Module InitGEM

  use WriteStatus, only : IOControl

  PRIVATE

  !*************************************************************************
  !****t* initGEM/Nuclide
  ! NAME
  ! type Nuclide
  ! SOURCE
  !
  type Nuclide
     integer                            :: ZahlN        !neutron number
     integer                            :: ZahlZ        !proton number
     real                               :: MassEx       !mass excess
     real                               :: Bind         !binding energy
     real                               :: spin         !groundstate spin
     character(6)                       :: Name         !name of element
     integer, allocatable, dimension(:) :: Ejectile     !index-position of ejectile
     integer, allocatable, dimension(:) :: Daughter     !index-position of daughter
     logical, allocatable, dimension(:) :: ChannelStatus !status of decay channel
  End Type Nuclide
  ! 
  ! NOTES
  ! Constructor for nuclide properties
  !*************************************************************************

  !*************************************************************************
  !****g* initGEM/NucMax
  ! NAME
  ! integer, parameter :: NucMax
  ! PURPOSE
  ! Maximum number of nuclides considered in evaporation model
  ! NOTES
  ! DO NOT CHANGE THIS NUMBER!!!
  !*************************************************************************
  integer, parameter :: NucMax=2467

  !*************************************************************************
  !****g* initGEM/Elements
  ! NAME
  ! type(Nuclide), Allocatable,dimension (:),save :: Elements
  ! PURPOSE
  ! Database for elements properties
  ! USAGE
  ! * To get e.g. the mass excess of the nuclide IN type "Elements(IN)%MassEx". 
  !   The ID's are choosen according to NNDC-Table.
  ! * The array indizes of "Elements" run from NucMax=1n up to 
  !   NucMax=208Pb (IN=1,...,NucMax).
  !*************************************************************************
  type(Nuclide), Allocatable, dimension(:), save :: Elements

  !*************************************************************************
  !****g* initGEM/NucPosition
  ! NAME
  ! Character, Allocatable,dimension(:,:),save :: NucPosition
  ! PURPOSE
  ! (N,Z)-Matrix providing with the position (in the nndc.dat-file) of a 
  ! stable (according N/Z-ratio) nuclide. Negative Value of a matrix 
  ! element indicates not existing nulcides.
  !
  !*************************************************************************
  integer, Allocatable, dimension(:,:), save :: NucPosition


  PUBLIC :: Nuclide  !Nuclide type: NuclideProperties constructor
  PUBLIC :: Elements !Nuclides with their properties
  PUBLIC :: NucPosition !Position index (N,Z)-array of each element
  PUBLIC :: Init_GEM,SetDecays,End_GEM !subroutines

  contains

    !***********************************************************************
    !****s* initGEM/init_GEM
    ! NAME
    ! subroutine init_GEM
    ! 
    ! PURPOSE
    ! * Reading of data-tables from NNDC (mass excess, spin)
    ! * Initialize the "Elements(:)%" vector
    ! * Initialize the "NucPosition(:,:)" array
    ! * Initialize Decaying channels for different (A,Z)-combinations
    ! NOTES
    ! * Initialization of (1:MaxN,1:MaxZ)-Matrix "NucPosition" is done here. 
    ! * Variables "MaxN" and "MaxZ" give the max. possible value of 
    !   neutron- and proton-content of a cluster. Its status 
    !   (stable/unstable) according N/Z-ratio is initialized here. 
    ! * NucPosition is needed in checking the correct N/Z-ratio of 
    !   produced clusters, before starting the evaporation procedure.
    !***********************************************************************
    subroutine Init_GEM
    !***********************************************************************
      use InputCoalescence, only : PathToClusterInput
      implicit none
      integer          :: ios,ioe,status,i,j,N,Z,MaxN,MaxZ
      character(len=6) :: Element
      real             :: MassEx,Bind,Spin,Bind_exp,Bind_woCoul
      !---------------------------------------------------------------------
      open(Unit=10,File=trim(PathToClusterInput)// 'nndc.dat', & 
           & Status='old', Action='read',Iostat=ios)
      if (ios /= 0) then
         write(*,*) 'from InputGEM: NNDC-File Open failed: ios = ',ios
         write(*,*) 'from InputGEM: !!! Termination of program NOW !!!'
         STOP
      endif
      !---------------------------------------------------------------------
      Allocate(Elements(1:NucMax),STAT=status)
      call IOControl(1,status,'Init_GEM','Elements')
      !---------------------------------------------------------------------
      ! Construct type(Nuclide) Elements
      !---------------------------------------------------------------------
      rewind(10)
      do i=1,45
         read(10,*) !skeep comments
      end do
      do i=1,NucMax !-->Maximum 2467 elements!!!
         read(10,5) N,Z,Element,MassEx,Bind_exp,Spin,Bind,Bind_woCoul
         Elements(i)%ZahlN  = N
         Elements(i)%ZahlZ  = Z
         Elements(i)%MassEx = MassEx*0.001 !MeV
         Elements(i)%Bind   = Bind*0.001   !MeV
         Elements(i)%Spin   = Spin
         Elements(i)%Name   = Element
      end do
5     format(i4,i5,1x,a6,f14.5,2x,f12.5,6x,f7.4,f14.5,f14.5)
      close(unit=10,status='keep',iostat=ioe)
      if(ioe /= 0) then
         write(*,*) 'from InputGEM: cannot close NNDC-File, IOSTAT = ',ioe
         write(*,*) 'from InputGEM: !!! Termination of program NOW !!!'
         STOP
      endif
      !---------------------------------------------------------------------
      MaxN = 130  !maxval(Elements%ZahlN)
      MaxZ = 90   !maxval(Elements%ZahlZ)
      Allocate(NucPosition(0:MaxN,0:MaxZ),STAT=status)
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
    !***********************************************************************
    end subroutine Init_GEM !***********************************************
    !***********************************************************************

    !***********************************************************************
    !****s* GEM_Evaporation/SetDecays
    ! NAME
    ! subroutine SetDecays
    ! 
    ! PURPOSE
    ! * Defines valid decay channels for each nuclide
    ! * Important for Monte-Carlo decision of decaying cluster
    ! * The list of elements is taken from the mass.mas03-table (cf. init_GEM)
    !***********************************************************************
    subroutine SetDecays
    !***********************************************************************
      implicit none
      integer              :: i,j,k,Ndiff,Zdiff,channel,Ai,Aj
      integer              :: tempN,tempZ,MaxVal,MaxDaughter
      integer              :: EjectilePosition,DaughterPosition
      logical              :: flag
      !MaxEmission gives the max. number of decays. In original GEM-model 
      !only up to 4He-dacays have been consideren. Here we consider many 
      !more...The choise of MaxEmission just depends on the CPU-time consumed. 
      !Its also possible to reduced it.
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
               endif
            end do Find1
            if (.not.flag) cycle
            DaughterPosition = NucPosition(tempN,tempZ)
            Elements(i)%Ejectile(channel)      = EjectilePosition
            Elements(i)%Daughter(channel)      = DaughterPosition
            Elements(i)%ChannelStatus(channel) = .true.
            channel = channel + 1
         end do DecayChannels
      end do NuclideLoop
    !***********************************************************************
    end subroutine SetDecays
    !***********************************************************************


    !***********************************************************************
    subroutine End_GEM
    !***********************************************************************
      implicit none
      integer :: status,i

      do i=1,NucMax
         Deallocate(Elements(i)%Ejectile,STAT=status)
         if(status /= 0) then
            write(*,*) 'deallocation of variable "%Ejectile" NOT successfull'
            write(*,*) '!!! TERMINATION NOW...!!!'
            STOP
         endif
         Deallocate(Elements(i)%Daughter,STAT=status)
         if(status /= 0) then
            write(*,*) 'deallocation of variable "%Daughter" NOT successfull'
            write(*,*) '!!! TERMINATION NOW...!!!'
            STOP
         endif
      end do

      Deallocate(Elements,STAT=status)
      if(status /= 0) then
         write(*,*) 'deallocation of variable "Elements" NOT successfull'
         write(*,*) '!!! TERMINATION NOW...!!!'
         STOP
      endif

      Deallocate(NucPosition,STAT=status)
      if(status /= 0) then
         write(*,*) 'deallocation of variable "NucPosition" NOT successfull'
         write(*,*) '!!! TERMINATION NOW...!!!'
         STOP
      endif


    !***********************************************************************
    end subroutine End_GEM !************************************************
    !***********************************************************************


  end Module initGEM



