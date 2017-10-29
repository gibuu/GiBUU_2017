!******************************************************************************
!****m* /VolumeElements
! NAME
! module VolumeElements
!
! PURPOSE
! This module defines all stuff necessary to seperate the whole
! interaction volume in different "volume elements" (also called "VE-cells")
!
! This is necessary for the implementation of "local ensemble" runs.
!
! INPUTS
! (none)
!******************************************************************************
module VolumeElements

  use particleDefinition
  use particlePointerListDefinition
  use particlePointerList

  implicit none

  private

  !****************************************************************************
  !****t* VolumeElements/tVolumeElements
  ! PURPOSE
  ! Discretize the whole possible space volume and hold in a 3D array
  ! Lists of particles, which are at the moment in a given coordinate cell
  !
  ! NOTES
  ! Die Arrays "zCoordFilled_xxx" sollen eine Abkürzung für die Loops
  ! über die z-Koordinate darstellen: Wenn für eine z-Koordinate
  ! für keine x- oder y-Zelle Einträge vorhanden sind, dann kann man
  ! diese z-Koord auch ganz schnell überspringen!
  !
  ! Denkbar wäre auch eine Verbesserung durch die Einführung von
  !   integer, dimension(3,2) :: iRange_Used
  ! wobei alle Schleifen statt über
  !   iRange(i,1)..iRange(i,2)
  ! über
  !   iRange_Used(i,1)..iRange_Used(i,2)
  ! laufen würden.
  ! Hierbei können aber nur zusammenhängende Bereiche benutzt werden,
  ! wodurch das benutzte Modell wiederum starken Auftrieb bekommt!!!
  !
  ! SOURCE
  !
  type tVolumeElements
     type(tParticleList), DIMENSION(:,:,:), ALLOCATABLE :: VE_real ! ~ 100 MB
     type(tParticleList), DIMENSION(:,:,:), ALLOCATABLE :: VE_pert ! ~ 100 MB
     real, dimension(3)      :: Delta
     real, dimension(3,2)    :: Range
     integer, dimension(3,2) :: iRange
     logical, dimension(:), ALLOCATABLE :: zCoordFilled_real
     logical, dimension(:), ALLOCATABLE :: zCoordFilled_pert
  end type tVolumeElements
  !****************************************************************************

  !****************************************************************************
  !****ig* VolumeElements/tVE
  ! PURPOSE
  ! The one and only instance of a tVolumeElements-object
  !
  ! SOURCE
  !
  type(tVolumeElements),save :: tVE
  !****************************************************************************

  !****************************************************************************
  !*****ig* VolumeElements/boxshift
  ! PURPOSE
  ! constant shift of the grid
  !
  ! This is needed e.g. for box calculaions, where one needs alignment of
  ! the grid to the box boundaries
  !
  ! SOURCE
  !
  real,save :: boxshift = 0.0
  !
  ! NOTES
  ! This value is automatically adjusted in the INIT routine by looking
  ! at inputGeneral/eventType
  !****************************************************************************

  integer, dimension(3) :: iPart
  type(tParticleListNode), POINTER :: pPert, pReal

  integer,save :: nEnsemble=0

  ! needed for pointer arithmetic:
  type(particle), POINTER, save :: pPert11,pPert12,pPert21
  type(particle), POINTER, save :: pReal11,pReal12,pReal21


  PUBLIC:: VolumeElements_INIT, VolumeElements_CLEAR, cleanUp
  PUBLIC:: VolumeElements_CLEAR_Pert, VolumeElements_CLEAR_Real
  PUBLIC:: VolumeElements_SETUP_Pert, VolumeElements_SETUP_Real
  PUBLIC:: VolumeElements_Statistics
  PUBLIC:: VolumeElements_InitGetPart_RealPert
  PUBLIC:: VolumeElements_InitGetPart_RealReal
  PUBLIC:: VolumeElements_GetPart_RealPert
  PUBLIC:: VolumeElements_GetPart_RealReal
  PUBLIC:: VolumeElements_boxSize
  PUBLIC:: VolumeElements_NukSearch
  PUBLIC:: VolumeElements_3Body

  logical, save :: initFlag=.true.


contains
  !****************************************************************************
  !****f* VolumeElements/VolumeElements_boxSize
  ! PURPOSE
  ! Returns size of box used for the local ensemble method, unit of fm^3
  ! OUTPUT
  ! (function value)
  !****************************************************************************
  real function VolumeElements_boxSize()

    if (initFlag) call VolumeElements_INIT

    VolumeElements_boxSize=tVE%Delta(1)*tVE%Delta(2)*tVE%Delta(3)
  end function VolumeElements_boxSize


  subroutine cleanUp
    call VolumeElements_CLEAR()
    if (allocated(tVE%VE_real)) deallocate(tVE%VE_real)
    if (allocated(tVE%VE_pert)) deallocate(tVE%VE_pert)
    if (allocated(tVE%zCoordFilled_real)) deallocate(tVE%zCoordFilled_real)
    if (allocated(tVE%zCoordFilled_pert)) deallocate(tVE%zCoordFilled_pert)
  end subroutine


  !****************************************************************************
  !****s* VolumeElements/VolumeElements_INIT
  ! NAME
  ! subroutine VolumeElements_INIT()
  !
  ! PURPOSE
  ! This routine initializes the tVE-instance.
  ! Initial sizes are set and all memory allocation is done.
  !
  ! NOTES
  ! The maximum volume size and also the volume elements size is still
  ! hard wired. maybe some more sophisticated init should be realized.
  !****************************************************************************
  subroutine VolumeElements_INIT()
    use densityModule, only: gridsize
    use inputGeneral, only: eventType
    use EventTypes

    integer :: i,j,k

    initFlag=.false. ! to remember that routine has been called

    select case (eventType)
    case (InABox,InABox_pion,InABox_delta,Box)
       boxshift = 0.5

       tVE%Range(:,1) = -gridsize(:) ! lower bound
       tVE%Range(:,2) =  gridsize(:) ! upper bound

    case default
       tVE%Range(1,1:2) = (/-40., 40./) ! x-Range
       tVE%Range(2,1:2) = (/-40., 40./) ! y-Range
       tVE%Range(3,1:2) = (/-40., 40./) ! z-Range

    end select


!    tVE%Delta = (/0.25,0.25,0.25/) ! x-, y-, z-Binning
    tVE%Delta = (/0.50,0.50,0.50/) ! x-, y-, z-Binning
!    tVE%Delta = (/1.00,1.00,1.00/) ! x-, y-, z-Binning

    do i=1,3
       tVE%iRange(i,1:2) = nint(tVE%Range(i,1:2)/tVE%Delta(i)+boxshift)
       tVE%iRange(i,1) = tVE%iRange(i,1)-1 ! for security
       tVE%iRange(i,2) = tVE%iRange(i,2)+1 ! for security
    end do

    write(*,'(79("#"))')
    write(*,'(A,3f9.4)') '  VE: Delta  = ',tVE%Delta
    write(*,'(A,6f9.4)') '  VE: Range  = ',tVE%Range
    write(*,'(A,6i9)')   '  VE: iRange = ',tVE%iRange
    write(*,*) ' VE: Size   = ',&
         (tVE%iRange(1,2)-tVE%iRange(1,1))* &
         (tVE%iRange(2,2)-tVE%iRange(2,1))* &
         (tVE%iRange(3,2)-tVE%iRange(3,1)), &
         ' Entries'
    write(*,'(79("#"))')

    allocate(tVE%VE_real(tVE%iRange(1,1):tVE%iRange(1,2),&
         tVE%iRange(2,1):tVE%iRange(2,2),&
         tVE%iRange(3,1):tVE%iRange(3,2)))

    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             call ParticleList_INIT(tVE%VE_real(i,j,k))
          end do
       end do
    end do

    allocate(tVE%zCoordFilled_real(tVE%iRange(3,1):tVE%iRange(3,2)))

    tVE%zCoordFilled_real = .FALSE.

    allocate(tVE%VE_pert(tVE%iRange(1,1):tVE%iRange(1,2),&
         & tVE%iRange(2,1):tVE%iRange(2,2),&
         & tVE%iRange(3,1):tVE%iRange(3,2)))

    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             call ParticleList_INIT(tVE%VE_pert(i,j,k))
          end do
       end do
    end do

    allocate(tVE%zCoordFilled_pert(tVE%iRange(3,1):tVE%iRange(3,2)))

    tVE%zCoordFilled_pert = .FALSE.

  end subroutine VolumeElements_INIT

  !****************************************************************************
  !****s* VolumeElements/VolumeElements_CLEAR_Pert
  ! NAME
  ! subroutine VolumeElements_CLEAR_Pert()
  !
  ! PURPOSE
  ! This routine resets the tVE-instance for perturbative particles
  !****************************************************************************
  subroutine VolumeElements_CLEAR_Pert()

    integer :: i,j,k

    if (initFlag) call VolumeElements_INIT

    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       if (.not.tVE%zCoordFilled_pert(k)) CYCLE
       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             call ParticleList_CLEAR(tVE%VE_pert(i,j,k))
          end do
       end do
    end do
    tVE%zCoordFilled_pert = .FALSE.

  end subroutine VolumeElements_CLEAR_Pert

  !****************************************************************************
  !****s* VolumeElements/VolumeElements_CLEAR_Real
  ! NAME
  ! subroutine VolumeElements_CLEAR_Real()
  !
  ! PURPOSE
  ! This routine resets the tVE-instance for real particles.
  !****************************************************************************
  subroutine VolumeElements_CLEAR_Real()

    integer :: i,j,k

    if (initFlag) call VolumeElements_INIT

    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       if (.not.tVE%zCoordFilled_real(k)) CYCLE
       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             call ParticleList_CLEAR(tVE%VE_real(i,j,k))
          end do
       end do
    end do
    tVE%zCoordFilled_real = .FALSE.

  end subroutine VolumeElements_CLEAR_Real

  !****************************************************************************
  !****s* VolumeElements/VolumeElements_CLEAR
  ! NAME
  ! subroutine VolumeElements_CLEAR()
  !
  ! PURPOSE
  ! This routine resets the tVE-instance, both for real and perturbative
  ! particles
  !****************************************************************************
  subroutine VolumeElements_CLEAR()
    call VolumeElements_CLEAR_Real
    call VolumeElements_CLEAR_Pert
  end subroutine VolumeElements_CLEAR


  !****************************************************************************
  !****s* VolumeElements/VolumeElements_SETUP_Pert
  ! NAME
  ! subroutine VolumeElements_SETUP_Pert(PartVec)
  !
  ! PURPOSE
  ! Build up the tVE structure of perturbative particles.
  !
  ! In order to introduce some randomness, particles are prepended/appended
  ! to the list according the outcome of a random generator.
  !****************************************************************************
  subroutine VolumeElements_SETUP_Pert(PartVec)
    use output
    use random, only: rn_truefalse

    type(particle), TARGET:: PartVec(:,:)

    integer :: i,j,k, ii(3)
    type(particle), POINTER :: pPart
    integer :: nPart

!    write(*,*) 'VolumeElements_SETUP_Pert: Start'

    if (initFlag) then
       call VolumeElements_Init
       call VolumeElements_Clear
    end if

    nPart = 0

    nEnsemble = Size(PartVec,dim=1)

    ensemble_loop : do i=1,Size(PartVec,dim=1)
       index_loop : do j=1,Size(PartVec,dim=2)
          if (PartVec(i,j)%ID < 0) exit index_loop
          if (PartVec(i,j)%ID == 0) cycle index_loop

          do k=1,3
             ii(k) = nint(PartVec(i,j)%position(k)/tVE%Delta(k)+boxshift)

             if ((ii(k)<=tVE%iRange(k,1)).or.(ii(k)>=tVE%iRange(k,2))) then
!                write(*,*) '#### Particle not in volume:',&
!                     & k,PartVec(i,j)%position(k)
                cycle index_loop
             end if
          end do
          pPart => PartVec(i,j)

          ! introduce some randomness:
          if (rn_truefalse()) then
             call ParticleList_APPEND(tVE%VE_Pert(ii(1),ii(2),ii(3)), pPart)
          else
             call ParticleList_PREPEND(tVE%VE_Pert(ii(1),ii(2),ii(3)), pPart)
          end if
          nPart = nPart + 1

          tVE%zCoordFilled_pert(ii(3)) = .TRUE.

       end do index_loop
    end do ensemble_loop

!    write(*,*) '...particles inserted: ',nPart

    ! pointers needed for pointer arithmetic:
    pPert11 => PartVec(1,1)
    if (Size(PartVec,dim=1)>1) then
      pPert21 => PartVec(2,1)
    else
      ! Only 1 ensemble
      pPert21 => PartVec(1,1) ! dummy value
    end if
    if (Size(PartVec,dim=2)>1) then
      pPert12 => PartVec(1,2)
    else
      ! Only 1 particle per ensemble. Routines probably all not working!!!
      pPert12 => PartVec(1,1) ! dummy value
    end if

    if (DoPR(1)) write(*,*) 'VolumeElements_SETUP_Pert: particles inserted: ',nPart

  end subroutine VolumeElements_SETUP_Pert


  !****************************************************************************
  !****s* VolumeElements/VolumeElements_SETUP_Real
  ! NAME
  ! subroutine VolumeElements_SETUP_Real(PartVec)
  !
  ! PURPOSE
  ! Build up the tVE structure of real particles.
  !
  ! In order to introduce some randomness, particles are prepended/appended
  ! to the list according the outcome of a random generator.
  !****************************************************************************
  subroutine VolumeElements_SETUP_Real(PartVec)
    use output
    use random, only: rn_truefalse

    type(particle), TARGET:: PartVec(:,:)

    integer :: i,j,k, ii(3)
    type(particle), POINTER :: pPart
    integer :: nPart

!    write(*,*) 'VolumeElements_SETUP_Real: Start'

    if (initFlag) then
       call VolumeElements_Init
       call VolumeElements_Clear
    end if


    nPart = 0

    ensemble_loop : do i=1,Size(PartVec,dim=1)
       index_loop : do j=1,Size(PartVec,dim=2)
          if (PartVec(i,j)%ID < 0) exit index_loop
          if (PartVec(i,j)%ID == 0) cycle index_loop

          do k=1,3
             ii(k) = nint(PartVec(i,j)%position(k)/tVE%Delta(k)+boxshift)

             if ((ii(k)<=tVE%iRange(k,1)).or.(ii(k)>=tVE%iRange(k,2))) then
!                write(*,*) '#### Particle not in volume:',&
!                     & k,PartVec(i,j)%position(k)
                cycle index_loop
             end if
          end do
          pPart => PartVec(i,j)

          ! introduce some randomness:
          if (rn_truefalse()) then
             call ParticleList_APPEND(tVE%VE_Real(ii(1),ii(2),ii(3)), pPart)
          else
             call ParticleList_PREPEND(tVE%VE_Real(ii(1),ii(2),ii(3)), pPart)
          end if
          nPart = nPart + 1

          tVE%zCoordFilled_real(ii(3)) = .TRUE.

       end do index_loop
    end do ensemble_loop

!    write(*,*) '...particles inserted: ',nPart

    ! pointers needed for pointer arithmetic:
    pReal11 => PartVec(1,1)
    if (Size(PartVec,dim=1)>1) then
      pReal21 => PartVec(2,1)
    else
      ! Only 1 ensemble.
      pReal21 => PartVec(1,1) ! dummy value
    end if
    if (Size(PartVec,dim=2)>1) then
      pReal12 => PartVec(1,2)
    else
      ! Only 1 particle per ensemble. Routines probably all not working!!!
      pReal12 => PartVec(1,1) ! dummy value
    end if

    if (DoPR(1)) write(*,*) 'VolumeElements_SETUP_Real: particles inserted: ',nPart

  end subroutine VolumeElements_SETUP_Real


  !****************************************************************************
  !****s* VolumeElements/VolumeElements_Statistics
  ! NAME
  ! subroutine VolumeElements_Statistics
  !
  ! PURPOSE
  ! This is a routine to produce some statistical informations about
  ! the elements in the tVE instance.
  !
  ! This routine is only for trial/documentational purposes.
  !****************************************************************************
  subroutine VolumeElements_Statistics

    integer :: i,j,k, ii, hh

    integer, dimension(-1:201) :: histPartR, histPartP

    integer, dimension(-1:201) :: histCollRR, histCollPR

    real :: NN(0:2)

    NN = 0.0

    hh = (tVE%iRange(1,2)-tVE%iRange(1,1)) * (tVE%iRange(2,2)-tVE%iRange(2,1))

    histPartR = 0
    histPartP = 0
    histCollRR = 0
    histCollPR = 0


    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       if (.not.tVE%zCoordFilled_real(k)) then
          histPartR(-1) = histPartR(-1) + hh
          CYCLE
       end if
       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             ii = min(tVE%VE_real(i,j,k)%nEntries, 201)
             histPartR(ii) = histPartR(ii)+1
             if (ii > 0) then
                NN = NN + (/ 1.0, ii*1.0, ii*(ii-1)*1.0 /)
             end if
          end do
       end do
    end do


    write(*,*) NN
    if (NN(0) > 0) then
       NN(1)=NN(1)/NN(0)
       NN(2)=NN(2)/NN(0)
    end if
    write(*,*) 'real-real: <N(N-1)>/(<N>(<N>-1)) :',&
         & NN(2)/( NN(1)*(NN(1)-1) )

    !    write(311,'(210(i9," "))') histPartR

!!$    do i=-1,100
!!$       write(411,*) i,histPartR(i)
!!$    end do
!!$    write(411,*)
!!$    write(411,*)



    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       if (.not.tVE%zCoordFilled_pert(k)) then
          histPartP(-1) = histPartP(-1) + hh
          CYCLE
       end if
       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             ii = min(tVE%VE_pert(i,j,k)%nEntries, 201)
             histPartP(ii) = histPartP(ii)+1
          end do
       end do
    end do

!    write(312,'(210(i9," "))') histPartP

    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       if (.not.tVE%zCoordFilled_real(k)) CYCLE

       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             ii = (tVE%VE_real(i,j,k)%nEntries * (tVE%VE_real(i,j,k)%nEntries - 1) )/2
             ii = max(min(ii,201),0)
             histCollRR(ii) = histCollRR(ii)+1
          end do
       end do
    end do

!    write(313,'(210(i9," "))') histCollRR


    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       if (.not.tVE%zCoordFilled_pert(k)) CYCLE
       if (.not.tVE%zCoordFilled_real(k)) CYCLE

       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             ii = min(tVE%VE_pert(i,j,k)%nEntries * tVE%VE_real(i,j,k)%nEntries, 201)
             histCollPR(ii) = histCollPR(ii)+1
          end do
       end do
    end do

!    write(314,'(210(i9," "))') histCollPR

!    write(*,*) 'VolumeElements_Statistics: End'


  end subroutine VolumeElements_Statistics


  !****************************************************************************
  !****if* VolumeElements/FindNextVE_RealPert
  ! NAME
  ! logical function FindNextVE_RealPert
  !
  ! PURPOSE
  ! This routine searches for the next volume element with non-vanishing
  ! number of real- and perturbative-test-particles.
  ! It remembers the values of the last call and finds really the "next"
  ! cell.
  !
  ! It proceeds via first increasing the x- and y-coordinates and then
  ! stepping to the next z-coordinate.
  !
  ! INPUTS
  ! * the static stored array iPart(1:3)
  !
  ! OUTPUT
  ! * the static stored array iPart(1:3) (changed!)
  ! * function value: .false. -> no more tVE-cells possible
  !****************************************************************************
  logical function FindNextVE_RealPert()

    FindNextVE_RealPert = .false. ! default return value
    do
       iPart(1)=iPart(1)+1
       if (iPart(1) > tVE%iRange(1,2)) then
          iPart(1) = tVE%iRange(1,1)
          iPart(2) = iPart(2)+1
          if (iPart(2) > tVE%iRange(2,2)) then
             iPart(2) = tVE%iRange(2,1)
             do
                iPart(3) = iPart(3)+1
                if (iPart(3) > tVE%iRange(3,2)) then
                   iPart(3) = tVE%iRange(3,1)
                   return ! --> false
                end if
                if (.not.tVE%zCoordFilled_pert(iPart(3))) cycle
                if (.not.tVE%zCoordFilled_real(iPart(3))) cycle
                exit
             end do

          end if
       end if

       if (tVE%VE_pert(iPart(1),iPart(2),iPart(3))%nEntries == 0) cycle
       if (tVE%VE_real(iPart(1),iPart(2),iPart(3))%nEntries == 0) cycle

       pReal => tVE%VE_real(iPart(1),iPart(2),iPart(3))%first
       pPert => tVE%VE_pert(iPart(1),iPart(2),iPart(3))%first

       FindNextVE_RealPert = .true.
       exit
    end do

  end function FindNextVE_RealPert

  !****************************************************************************
  !****if* VolumeElements/FindNextVE_RealReal
  ! NAME
  ! logical function FindNextVE_RealReal
  !
  ! PURPOSE
  ! This routine searches for the next volume element with non-vanishing
  ! number of real-test-particles.
  ! It remembers the values of the last call and finds really the "next"
  ! cell.
  !
  ! It proceeds via first increasing the x- and y-coordinates and then
  ! stepping to the next z-coordinate.
  !
  ! INPUTS
  ! * the static stored array iPart(1:3)
  !
  ! OUTPUT
  ! * the static stored array iPart(1:3) (changed!)
  ! * function value: .false. -> no more tVE-cells possible
  !****************************************************************************
  logical function FindNextVE_RealReal()

    FindNextVE_RealReal = .false. ! default return value
    do
       iPart(1)=iPart(1)+1
       if (iPart(1) > tVE%iRange(1,2)) then
          iPart(1) = tVE%iRange(1,1)
          iPart(2) = iPart(2)+1
          if (iPart(2) > tVE%iRange(2,2)) then
             iPart(2) = tVE%iRange(2,1)
             do
                iPart(3) = iPart(3)+1
                if (iPart(3) > tVE%iRange(3,2)) then
                   iPart(3) = tVE%iRange(3,1)
                   return ! --> false
                end if
                if (.not.tVE%zCoordFilled_real(iPart(3))) cycle
                exit
             end do

          end if
       end if

       if (tVE%VE_real(iPart(1),iPart(2),iPart(3))%nEntries < 2) cycle

       pReal => tVE%VE_real(iPart(1),iPart(2),iPart(3))%first

       FindNextVE_RealReal = .true.
       exit
    end do

  end function FindNextVE_RealReal


  !****************************************************************************
  !****s* VolumeElements/VolumeElements_InitGetPart_RealPert
  ! NAME
  ! subroutine VolumeElements_InitGetPart_RealPert
  !
  ! PURPOSE
  ! initialize the "GetPart_RealPert"-routines
  !
  ! NOTES
  ! this routine sets the x-,y-,z-indizes in such a way, that a call to
  ! "FindNextVE_RealPert" will start at the very first cell.
  !****************************************************************************
  subroutine VolumeElements_InitGetPart_RealPert
!    use callstack, only: traceBack

    NULLIFY(pPert,pReal)
    iPart(1:2) = tVE%iRange(1:2,2) ! really start with max value
    iPart(3) = tVE%iRange(3,1)-1
    if (.not.FindNextVE_RealPert()) then
       write(*,*) 'VolumeElements_InitGetPart_RealPert: warning'
!       call traceback
    end if

  end subroutine VolumeElements_InitGetPart_RealPert

  !****************************************************************************
  !****s* VolumeElements/VolumeElements_InitGetPart_RealReal
  ! NAME
  ! subroutine VolumeElements_InitGetPart_RealReal
  !
  ! PURPOSE
  ! initialize the "GetPart_RealPert"-routines
  !
  ! NOTES
  ! this routine sets the x-,y-,z-indizes in such a way, that a call to
  ! "FindNextVE_RealReal" will start at the very first cell.
  !****************************************************************************
  subroutine VolumeElements_InitGetPart_RealReal

    NULLIFY(pReal)
    iPart(1:2) = tVE%iRange(1:2,2) ! really start with max value
    iPart(3) = tVE%iRange(3,1)-1
    if (.not.FindNextVE_RealReal()) then
       write(*,*) 'VolumeElements_InitGetPart_RealReal: warning'
!       stop
    end if

  end subroutine VolumeElements_InitGetPart_RealReal



#ifdef f95
#define LOCcmd pointer
#elif defined nagfor

  integer function myLOC(x)
    use, intrinsic :: iso_c_binding
    implicit none
    type(particle), pointer :: x
    type(c_ptr) :: cptr
    integer(c_intptr_t) :: ptr
    cptr = C_LOC(x)
    ptr = transfer(cptr,ptr)
    myLOC = ptr
  end function

#define LOCcmd myLOC

#else
#define LOCcmd LOC
#endif


  !****************************************************************************
  !****f* VolumeElements/VolumeElements_GetPart_RealPert
  ! NAME
  ! logical function VolumeElements_GetPart_RealPert(Part1, Part2, nRealPart, iEns,iInd)
  !
  ! INPUTS
  !
  ! OUTPUT
  ! * type(particle), POINTER :: Part1, Part2 -- real and perturbative particle
  ! * integer :: nRealPart -- number of real particles in VE-cell
  ! * integer :: iEns,iInd -- coordinates of perturbative particle
  !
  ! PURPOSE
  ! This routine finds the next (possibly colliding?) pair
  ! of one perturbative particle and one real particle
  ! in the actual VE-cell.
  !
  ! If one stepped over all pert. particles in the given VE-cell, the next
  ! cell is choosen (cf. FindNextVE) and everything goes on.
  !
  ! If no "next VE-cell" is possible any more, the routine returns
  ! .false. as failure-indicator. (Otherwise always .true. is returned.)
  !
  ! NOTES
  ! * actually for a choosen VE-cell it returns the particle pairs
  !     (P_1,R_1), (P_2,R_2), ... (P_nPert, R_xxx)
  !   If there are more real particles than pert particles,
  !   "R_xxx" stands for "R_nPert". Otherwise, ie. if we have more
  !   pert particles than real particles, the loop restarts
  !   for the real particles
  !     (P_nReal,R_nReal), (P_nReal+1,R_1), ...
  !****************************************************************************
  logical function VolumeElements_GetPart_RealPert(Part1, Part2, nRealPart, iEns, iInd)
    type(particle), pointer :: Part1, Part2
    integer, intent(out)    :: nRealPart,iEns,iInd

    integer :: j2

    if (.not.associated(pPert)) then
       if (.not.FindNextVE_RealPert()) then
         VolumeElements_GetPart_RealPert = .false.
         return
      end if
    end if

    if (.not.associated(pReal)) then
       pReal => tVE%VE_real(iPart(1),iPart(2),iPart(3))%first
    end if

    Part1 => pReal%V
    Part2 => pPert%V
    nRealPart = tVE%VE_real(iPart(1),iPart(2),iPart(3))%nEntries

    j2 = int( (LOCcmd(Part2)-LOCcmd(pPert11)) / (LOCcmd(pPert12)-LOCcmd(pPert11)) )
    if (LOCcmd(pPert21) == LOCcmd(pPert11)) then
      iEns = 1
    else
      iEns = int((LOCcmd(Part2)-LOCcmd(pPert11)-j2*(LOCcmd(pPert12)-LOCcmd(pPert11)))/(LOCcmd(pPert21)-LOCcmd(pPert11))) + 1
    end if
    iInd=j2+1

    pPert => pPert%next
    pReal => pReal%next

    VolumeElements_GetPart_RealPert = .true.

  end function VolumeElements_GetPart_RealPert


  !****************************************************************************
  !****f* VolumeElements/VolumeElements_GetPart_RealReal
  ! NAME
  ! logical function VolumeElements_GetPart_RealReal(Part1, Part2, nRealPart, iEns1,iInd1, iEns2,iInd2)
  !
  ! INPUTS
  !
  ! OUTPUT
  ! * type(particle), POINTER :: Part1, Part2 -- real particles
  ! * integer :: nRealPart -- number of real particles in VE-cell
  ! * integer :: iEns1,iInd1 -- coordinates of particle1
  ! * integer :: iEns2,iInd2 -- coordinates of particle2
  !
  ! PURPOSE
  ! This routine finds the next (possibly colliding?) pair
  ! of two real particles  in the actual VE-cell.
  !
  ! If one stepped over all real particles in thi given VE-cell, the next
  ! cell s choosen (cf. FindNextVE_REalReal) and everything goes on.
  !
  ! If no "next VE-cell" is possible any more, the routine returns
  ! .false. as failure-indicator. (Otherwise always .true. is returned.)
  !
  ! NOTES
  ! * actually for a choosen VE-cell it returns the particle pairs
  !     (R_1,R_2), (R_3,R_4), ...
  !****************************************************************************
  logical function VolumeElements_GetPart_RealReal(Part1, Part2, nRealPart, iEns1,iInd1, iEns2,iInd2)

    type(particle), pointer :: Part1, Part2
    integer, intent(out)    :: nRealPart, iEns1,iInd1, iEns2,iInd2

    do
       if (.not.associated(pReal)) then
          if (.not.FindNextVE_RealReal()) then
             VolumeElements_GetPart_RealReal = .false.
             return
          end if
       end if

       Part1 => pReal%V
       pReal => pReal%next

       if (associated(pReal)) exit ! exit the loop
    end do

    Part2 => pReal%V
    pReal => pReal%next

    call GetEnsInd( Part1, iEns1, iInd1 )
    call GetEnsInd( Part2, iEns2, iInd2 )

    nRealPart = tVE%VE_real(iPart(1),iPart(2),iPart(3))%nEntries
    VolumeElements_GetPart_RealReal = .true.

  contains

    subroutine GetEnsInd(pPart, iEns,iInd)
      type(particle), POINTER, intent(in) :: pPart
      integer, intent(OUT) :: iEns, iInd

      integer :: j

      j = int( (LOCcmd(pPart)-LOCcmd(pReal11)) / (LOCcmd(pReal12)-LOCcmd(pReal11)) )
      if (LOCcmd(pReal21) == LOCcmd(pReal11)) then
        iEns = 1
      else
        iEns = int((LOCcmd(pPart)-LOCcmd(pReal11)-j*(LOCcmd(pReal12)-LOCcmd(pReal11)))/(LOCcmd(pReal21)-LOCcmd(pReal11))) + 1
      end if
      iInd=j+1
    end subroutine GetEnsInd

  end function VolumeElements_GetPart_RealReal


  !****************************************************************************
  !****s* VolumeElements/VolumeElements_NukSearch
  ! NAME
  ! subroutine VolumeElements_NukSearch(partIn,RadiusNukSearch,proton1,proton2,neutron1,neutron2,FlagOK)
  !
  ! PURPOSE
  ! This routine searches for two protons and two neutrons given in the
  ! "volume elements particle vector array" "VE_real" in the vicinty
  ! of the particle given by "partIn".
  !
  ! NOTES
  ! we are looking for 2 protons and 2 neutrons, i.e. for 4 nucleons, while
  ! only 2 nucleons are necessary: possible combinations are p+p, p+n, n+n.
  !
  ! INPUTS
  ! * type(particle)          :: partIn
  ! * real                    :: RadiusNukSearch
  !
  ! OUTPUT
  ! * type(particle), pointer :: proton1,proton2   -- Closest protons
  ! * type(particle), pointer :: neutron1,neutron2 -- Closest neutrons
  ! * logical                 :: FlagOK
  !
  !****************************************************************************
  subroutine VolumeElements_NukSearch(partIn,RadiusNukSearch,proton1,proton2,neutron1,neutron2,FlagOK)
    use GridOrdering
    use IDTable, only: nucleon
    use collisionNumbering, only: check_justCollided
    use constants, only: rhoNull
    use dichteDefinition
    use densityModule

    type(particle), intent(in) :: partIn
    real , intent(in) :: RadiusNukSearch
    type(particle), pointer :: proton1, proton2
    type(particle), pointer :: neutron1, neutron2
    logical, intent(out) :: FlagOK

    type(dichte) :: dens
    real,    save :: radiusMax
    logical, save :: DoInit_GridOrdering = .true.

    integer,dimension(0:1) :: nFound
    integer :: iRadius, iiRadius, nRadiusMax
    integer, save :: nRadius

    integer :: ii0(3), ii(3), k

    real, dimension(1:3) :: position,DeltaPos

    type(particle), pointer :: pPart
    type(tParticleListNode), pointer :: pNode

    FlagOK = .false.
    nFound = 0


    if (DoInit_GridOrdering) then
       call GridOrdering_INIT
       nRadius = ubound(nDistance,dim=1)-1

       DoInit_GridOrdering = .false.
    end if

    position = partIn%position
    dens = densityAt(position)

    if (dens%baryon(0).lt.5e-03) then  !If density too small, then no absorption
!       write(*,*) '#### VolumeElements_NukSearch: density too small:',&
!            & dens%baryon,'#####',position
       return
    end if

    radiusMax = min(RadiusNukSearch*(rhoNull/dens%baryon(0))**(1./3.),5.0)
    radiusMax = radiusMax/(float(nEnsemble))**(1./3.)

    nRadiusMax = (nint(radiusMax/minval(tVE%Delta(1:3))))**2 + 1

!    write(*,*) 'radiusMax...:',radiusMax,minval(tVE%Delta(1:3)),nRadiusMax,nRadius

    do k=1,3
       ii0(k) = nint(position(k)/tVE%Delta(k))

       if ((ii0(k)<=tVE%iRange(k,1)).or.(ii0(k)>=tVE%iRange(k,2))) then
          write(*,*) '#### VolumeElements_NukSearch: Particle not in volume:', &
               & k,partIn%position(:)
          return
       end if
    end do

    do iRadius=0,min(nRadius,nRadiusMax) ! = 0,1,2,3 ,...
       if (nDistance(iRadius,2) < 0) cycle ! for this radius no cells...

!       write(*,*) '#### iRadius = ',iRadius
!       if (iRadius>0) write(*,*) '#### VolumeElements_NukSearch: iRadius>0',iRadius

       call GridOrdering_RandomizeRadius(iRadius)

       ! loop over all possibilities to get iRadius=const:
       do iiRadius=nDistance(iRadius,1),nDistance(iRadius,2)

          ii(1:3) = ii0(1:3) + DeltaV(iDeltaV(iiRadius),1:3)
          pNode => tVE%VE_Real(ii(1),ii(2),ii(3))%first

          do
             if (.not.associated(pNode)) exit
             pPart => pNode%V
             pNode => pNode%next

             if (pPart%ID.ne.nucleon) cycle
             if (pPart%antiparticle) cycle
             if (check_justCollided(PartIn,pPart)) cycle
             if ((pPart%Charge < 0).or.(pPart%Charge>1)) cycle
             if (nFound(pPart%Charge) >= 2) cycle

             DeltaPos = position-pPart%position
             if (DOT_PRODUCT(DeltaPos,DeltaPos) > radiusMax**2) cycle

             select case (pPart%Charge)
             case (0)
                select case (nFound(0))
                case (0)
                   neutron1 => pPart
                case (1)
                   neutron2 => pPart
                case default
                   cycle
                end select
                nFound(0) = nFound(0)+1

             case (1)
                select case (nFound(1))
                case (0)
                   proton1 => pPart
                case (1)
                   proton2 => pPart
                case default
                   cycle
                end select
                nFound(1) = nFound(1)+1

             end select


             if (nFound(0)==2 .and. nFound(1)==2) then
                FlagOK = .TRUE.
                return
             end if

          end do
       end do
    end do

  end subroutine VolumeElements_NukSearch

  !****************************************************************************
  !****f* VolumeElements/VolumeElements_3Body
  ! NAME
  ! logical function VolumeElements_3Body(Decay3Body, isSameBool, doInit, iEns,iInd, Part1,Part2,Part3, scaleFak, mode)
  ! PURPOSE
  ! Find a triple of particles corresponding the given needs.
  !
  ! This routine starts with one cell, loops over all possibilities, then
  ! switchs to the next cell, and so on, until no triple can be found any more.
  ! Then the routine returns .false., otherwise .true.
  !
  ! INPUTS
  ! * type(tDecay3Body) :: Decay3Body -- the decay to study
  ! * integer :: isSameBool -- tricky number to hold the info, which particles
  !   are identical
  ! * integer :: mode -- select how nDo is choosen
  ! OUTPUT
  ! * integer, dimension(1:3) :: iEns -- coordinate 1 of the 3 particles
  ! * integer, dimension(1:3) :: iInd -- coordinate 2 of the 3 particles
  ! * type(particle), POINTER :: Part1, Part2, Part3 -- the particles
  ! * real :: scaleFak -- the factor nPossible/nDo
  !
  ! NOTES
  ! possible values for mode:
  ! * 1: nDo = nMax
  ! * 2: nDo = nMin
  !****************************************************************************
  logical function VolumeElements_3Body(Decay3Body, isSameBool, doInit, iEns,iInd, Part1,Part2,Part3, scaleFak, mode)
    use DecayChannels
    use particlePointerList
    use particlePointerListDefinition
    use callstack, only: traceBack

    type(tDecay3Body), intent(in) :: Decay3Body
    integer, intent(in) :: isSameBool
    logical, intent(inOut) :: doInit
    integer, intent(out), dimension(1:3) :: iEns,iInd
    type(particle), POINTER, intent(OUT) :: Part1, Part2, Part3
    real, intent(out) :: scaleFak ! see also FindFirst3
    integer, intent(in) :: mode

    type(tParticleList), save :: L1, L2, L3

    type(tParticleListNode), POINTER, save :: pNode1, pNode2, pNode3
    integer, save :: nDo
    real, save :: scaleFak_

    ! default return values:
    scaleFak = 0.0

    ! if it is the first call for all cells, then initialize
    if (doInit) then ! initialize
       NULLIFY(pReal)
       iPart(1:2) = tVE%iRange(1:2,2) ! really start with max value
       iPart(3) = tVE%iRange(3,1)-1
       doInit = .false.
    end if

    loop1: do
!       write(*,*) 'loop1'
       if (.not.associated(pReal)) then
          loop2: do
!             write(*,*) 'loop2'
             if (.not.FindNextVE_RealReal()) then ! failure, no more cells
                VolumeElements_3Body = .false.
                return
             end if
             if (FindFirst3()) exit loop2 ! particles found
          end do loop2
          exit loop1 ! exit the loop
       end if

       if (FindSecond3()) exit loop1 ! exit the loop
       NULLIFY(pReal)

    end do loop1

    call GetEnsInd( Part1, iEns(1), iInd(1) )
    call GetEnsInd( Part2, iEns(2), iInd(2) )
    call GetEnsInd( Part3, iEns(3), iInd(3) )

    scaleFak = scaleFak_

    VolumeElements_3Body = .true.
    return


  contains
    !**************************************************************************
    !****if* VolumeElements_3Body/FindFirst3
    ! NAME
    ! logical function FindFirst3()
    ! PURPOSE
    ! This internal routine has two purposes:
    ! * build up the particle lists L1, L2, L3 from the particles in the cell
    ! * return the first possible particle triple
    ! OUTPUT
    ! The return vaule is .true., if everything went smoothly.
    ! In addition, this routine sets the following internal variables:
    ! * Part1, Part2, Part3
    ! * pNode1, pNode2, pNode3 show to the next possible triple.
    ! * nMax indicates the number of triples to extract
    ! * scaleFak = nPossible/nMax
    ! This routine ensures, that Part1, Part2, Part3 do not accidentally
    ! point to a same particle.
    ! NOTES
    ! * a combinatorical factor is included in nPossible for indistinguishable
    !   particles
    ! * Particles are added to the lists randomly at the beginning or the end.
    !**************************************************************************
    logical function FindFirst3()

      use random, only: rn_truefalse

      type(tParticleListNode), POINTER :: pNode
      type(particle), POINTER :: pP
      integer :: nPossible, nMax, nMin

!      write(*,*) 'FindFirst'

      FindFirst3 = .false.

      scaleFak_ = 0.0 ! set default return value

      call ParticleList_CLEAR(L1)
      call ParticleList_CLEAR(L2)
      call ParticleList_CLEAR(L3)

      ! iterate over all particles in the cell, build up L1,L2,L3:

      pNode => pReal
      do
         if (.not.associated(pNode)) exit
         pP => pNode%V
         if ((pP%ID == Decay3Body%ID(1))&
              .and.(pP%Charge == Decay3Body%Charge(1))&
              .and.(pP%antiparticle.eqv.Decay3Body%isAnti(1))) then
            if (rn_truefalse()) then
               call ParticleList_APPEND(L1,pP)
            else
               call ParticleList_PREPEND(L1,pP)
            end if
         end if
         if ((pP%ID == Decay3Body%ID(2))&
              .and.(pP%Charge == Decay3Body%Charge(2))&
              .and.(pP%antiparticle.eqv.Decay3Body%isAnti(2))) then
            if (rn_truefalse()) then
               call ParticleList_APPEND(L2,pP)
            else
               call ParticleList_PREPEND(L2,pP)
            end if
         end if
         if ((pP%ID == Decay3Body%ID(3))&
              .and.(pP%Charge == Decay3Body%Charge(3))&
              .and.(pP%antiparticle.eqv.Decay3Body%isAnti(3))) then
            if (rn_truefalse()) then
               call ParticleList_APPEND(L3,pP)
            else
               call ParticleList_PREPEND(L3,pP)
            end if
         end if
         pNode => pNode%next
      end do

      if ((L1%nEntries==0).or.(L2%nEntries==0).or.(L3%nEntries==0)) then
         return  ! --> failure
      end if

      ! calculate combinatorical factors:

      select case (isSameBool)
      case (0) ! -- all particles different
         nPossible = L1%nEntries*L2%nEntries*L3%nEntries
         nMax = max(L1%nEntries,L2%nEntries,L3%nEntries)
         nMin = min(L1%nEntries,L2%nEntries,L3%nEntries)
      case (3) ! -- particle 1 and 2 the same
         nPossible = (L1%nEntries*(L1%nEntries-1)*L3%nEntries)/2
         nMax = max(L1%nEntries,L3%nEntries)
         nMin = min(L1%nEntries-1,L3%nEntries)
      case (6) ! -- particle 2 and 3 the same
         nPossible = (L2%nEntries*(L2%nEntries-1)*L1%nEntries)/2
         nMax = max(L2%nEntries,L1%nEntries)
         nMin = min(L2%nEntries-1,L1%nEntries)
      case (7) ! -- all particle the same
         nPossible = (L1%nEntries*(L1%nEntries-1)*(L1%nEntries-2))/6
         nMax = L1%nEntries
         nMin = max(0,L1%nEntries-2) ! could be < 0!
      case default
         call traceBack('not yet implemented')
      end select

      if (nPossible==0) return  ! --> failure

      nMax = min(nMax,nPossible)
      nMin = min(nMin,nPossible)

      select case (mode)
      case (1)
         nDo = nMax
      case (2)
         nDo = nMin
      case default
         call Traceback('wrong mode')
      end select

      if (nDo <= 0) return  ! --> failure

      scaleFak_ = real(nPossible) / nDo

!      write(*,'(4(A,i4,"  "),A,f9.3)') 'nPossible=',nPossible,'nDo=',nDo,'nMin=',nMin,'nMax=',nMax,'scaleFak=',scaleFak_

      ! set pointer to particle 1:

      pNode1 => L1%first
      Part1 => pNode1%V

      ! set pointer to particle 2:
      ! check, whether Part1 and Part2 point to the same particle

      pNode2 => L2%first
      do
         if (.not.associated(pNode2)) return ! --> failure
         Part2 => pNode2%V
         if (.not.associated(Part1,Part2)) exit ! okay
         pNode2 => pNode2%next
      end do

      ! set pointer to particle 3:
      ! check, whether Part1/2 and Part3 point to the same particle

      pNode3 => L3%first
      do
         if (.not.associated(pNode3)) return ! --> failure
         Part3 => pNode3%V
         if ((.not.associated(Part1,Part3)).and.(.not.ASSOCIATED(Part2,Part3))) exit ! okay
         pNode3 => pNode3%next
      end do

      FindFirst3 = .true.
      return

    end function FindFirst3

    !**************************************************************************
    !****if* VolumeElements_3Body/FindSecond3
    ! NAME
    ! logical function FindSecond3()
    ! PURPOSE
    ! This routine gives the 'next' triple.
    ! It loops over L1, L2 and L3, until it has genarated 'nDo' triples.
    ! This routine ensures, that Part1, Part2, Part3 do not accidentally
    ! point to a same particle.
    !**************************************************************************
    logical function FindSecond3()

      !      write(*,*) 'FindSecond'

      FindSecond3 = .false.

      do while (.not.FindSecond3)

         if (nDo==0) return ! --> failure
         nDo = nDo-1

         pNode1 => pNode1%next
         if (.not.associated(pNode1)) pNode1 => L1%first
         Part1 => pNode1%V

         pNode2 => pNode2%next
         do
            if (.not.associated(pNode2)) pNode2 => L2%first
            Part2 => pNode2%V
            if (.not.associated(Part1,Part2)) exit ! okay
            pNode2 => pNode2%next
         end do

         pNode3 => pNode3%next
         do
            if (.not.associated(pNode3)) pNode3 => L3%first
            Part3 => pNode3%V
            if ((.not.associated(Part1,Part3)).and.(.not.ASSOCIATED(Part2,Part3))) exit ! okay
            pNode3 => pNode3%next
         end do

         FindSecond3 = CheckFound3()
      end do
      return
    end function FindSecond3

    !**************************************************************************
    !****if* VolumeElements_3Body/CheckFound3()
    ! NAME
    ! logical function CheckFound3()
    ! PURPOSE
    ! This routine checks, whether the found triple is really according the
    ! input wishes. This is necessary, since a successfully performed
    ! 3->1 collision should delete the incoming mesons.
    !
    !**************************************************************************
    logical function CheckFound3()

      CheckFound3 = .false. ! default return value: failure

      if (Part1%ID /= Decay3Body%ID(1)) return
      if (Part2%ID /= Decay3Body%ID(2)) return
      if (Part3%ID /= Decay3Body%ID(3)) return

      if (Part1%Charge /= Decay3Body%Charge(1)) return
      if (Part2%Charge /= Decay3Body%Charge(2)) return
      if (Part3%Charge /= Decay3Body%Charge(3)) return

      if (Part1%antiParticle .neqv. Decay3Body%isAnti(1)) return
      if (Part2%antiParticle .neqv. Decay3Body%isAnti(2)) return
      if (Part3%antiParticle .neqv. Decay3Body%isAnti(3)) return

      CheckFound3 = .true.
    end function CheckFound3

    !**************************************************************************
    !****is* VolumeElements_3Body/GetEnsInd
    ! NAME
    ! subroutine GetEnsInd(pPart, iEns,iInd)
    ! PURPOSE
    ! needed for pointer arithmetics
    !**************************************************************************
    subroutine GetEnsInd(pPart, iEns,iInd)
      type(particle), POINTER, intent(in) :: pPart
      integer, intent(OUT) :: iEns, iInd

      integer :: j

      j  =int( (LOCcmd(pPart)-LOCcmd(pReal11))/(LOCcmd(pReal12)-LOCcmd(pReal11)))
      if (LOCcmd(pReal21) == LOCcmd(pReal11)) then
         iEns = 1
      else
         iEns=int( (LOCcmd(pPart)-LOCcmd(pReal11)-j*(LOCcmd(pReal12)-LOCcmd(pReal11)))/(LOCcmd(pReal21)-LOCcmd(pReal11)) )+1
      end if
      iInd=j+1
    end subroutine GetEnsInd

  end function VolumeElements_3Body

end module VolumeElements
