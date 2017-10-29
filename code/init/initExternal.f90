!******************************************************************************
!****m* /initExternal
! NAME
! module initExternal
! PURPOSE
! Initializes a hadronic system according to an external data file.
!******************************************************************************
module initExternal

  implicit none
  private

  !****************************************************************************
  !****g* externalSystem/inputFile
  ! SOURCE
  !
  character*100, save :: inputFile='source.inp'
  ! PURPOSE
  ! the name of the input file with hadrons to be propagated.
  !****************************************************************************

  !****************************************************************************
  !****g* externalSystem/DoPerturbative
  ! SOURCE
  logical, save :: DoPerturbative = .false.
  ! PURPOSE
  ! if true, the particles will be inserted into the perturbative particle
  ! vector, the real particles have to be initialized via some nucleus
  ! definition
  !****************************************************************************

  !****************************************************************************
  !****g* externalSystem/NumberingScheme
  ! SOURCE
  integer, save :: NumberingScheme = 1
  ! PURPOSE
  ! The way, how particles%event will be numbered:
  ! * 1: event = iPart, i.e. the particle number in the ensemble
  !   (historical, but does not work for fullensemble)
  ! * 2: event = -999 (should work for perturbative init)
  !****************************************************************************

  public :: initializeExternal, ExternalIsPerturbative

  logical, save :: initFlag=.true.

contains

  !****************************************************************************
  !****f* initExternal/ExternalIsPerturbative
  ! NAME
  ! logical function ExternalIsPerturbative()
  ! PURPOSE
  ! Returns the value of DoPerturbative
  !****************************************************************************
  logical function ExternalIsPerturbative()
    ExternalIsPerturbative = DoPerturbative
  end function ExternalIsPerturbative

  !****************************************************************************
  !****s* initExternal/initInput
  ! NAME
  ! subroutine initInput
  ! PURPOSE
  ! Reads input out of jobcard. Namelist 'externalSystem'.
  !****************************************************************************
  subroutine initInput
    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n* initExternal/externalSystem
    ! NAME
    ! NAMELIST externalSystem
    ! PURPOSE
    ! Includes the switches:
    ! * inputFile
    ! * DoPerturbative
    ! * NumberingScheme
    !**************************************************************************
    NAMELIST /externalSystem/ inputFile, DoPerturbative, NumberingScheme

    call Write_ReadingInput('externalSystem',0)
    rewind(5)
    read(5,nml=externalSystem)
    write(*,*) ' inputFile:  ', trim(inputFile)
    write(*,*) ' DoPerturbative : ', DoPerturbative
    write(*,*) ' NumberingScheme: ', NumberingScheme
    call Write_ReadingInput('externalSystem',1)

    initFlag = .false.

  end subroutine initInput

  !****************************************************************************
  !****s* initExternal/initializeExternal
  ! NAME
  ! subroutine initializeExternal(PartsReal,PartsPert)
  ! PURPOSE
  ! Read the particles from the file
  !****************************************************************************
  subroutine initializeExternal(PartsReal,PartsPert)
    use particleDefinition
    use insertion, only: GarbageCollection
    use output, only: Write_ReadingInput, Write_InitStatus
    use checks, only: ChecksSwitchRealRun

    type(particle), dimension(:,:), intent(inOut), target :: PartsReal
    type(particle), dimension(:,:), intent(inOut), target :: PartsPert
    type(particle), dimension(:,:), pointer :: pParts

    integer :: i,index,k,iens,IOS,count

    write(*,*)
    call Write_InitStatus('External',0)
    if (initFlag) call initInput

    call Write_ReadingInput(trim(inputFile),0)

    open(1,file=trim(inputFile),status='old')
!    open(2,file='source.chk',status='new')

    count = 0

    call ChecksSwitchRealRun(.not.DoPerturbative)
    if (DoPerturbative) then
       pParts => PartsPert
    else
       pParts => PartsReal
    end if

    ensemble_loop : do i=1,size(pParts,dim=1)

       index=0
       backspace(1)

       particle_loop : do

          index=index+1

          if (index.gt.size(pParts,dim=2)) then
            write(*,*) 'Particle vector too small. Stop in initializeExternal.'
            stop
          end if

          call setToDefault(pParts(i,index)) ! set particle to its default values
          call setNumber(pParts(i,index)) ! give each particle a unique number

          read(1,*,IOSTAT=IOS) pParts(i,index)%id,pParts(i,index)%charge,&
                            &pParts(i,index)%mass,&
                            &(pParts(i,index)%position(k), k=1,3),&
                            &(pParts(i,index)%momentum(k), k=1,3),iens

! 55        format(i4,1x,i2,1x,f5.3,3(1x,f8.3),3(1x,f8.3),1x,i5)


          if (IOS.lt.0) then ! E.o.f. is reached
            call setToDefault(pParts(i,index))
            exit ensemble_loop
          end if

          if (iens.ne.i) then
            call setToDefault(pParts(i,index))
            exit particle_loop
          end if

          ! Check input:
!!$          write(2,FMT=55) pParts(i,index)%id,pParts(i,index)%charge,&
!!$               &pParts(i,index)%mass,&
!!$               &(pParts(i,index)%position(k), k=1,3),&
!!$               &(pParts(i,index)%momentum(k), k=1,3),iens

          if (pParts(i,index)%id.lt.0) then
             pParts(i,index)%id=abs(pParts(i,index)%id)
             pParts(i,index)%antiparticle=.true.
          end if

          select case (NumberingScheme)
          case (1)
             pParts(i,index)%event=index ! This has to be changed in the full-ensemble mode
          case (2)
             pParts(i,index)%event=-999
          case (3)
             pParts(i,index)%event=count
          end select

          pParts(i,index)%perturbative = DoPerturbative

          count = count + 1

       end do particle_loop

    end do ensemble_loop

    call Write_ReadingInput(trim(inputFile),1)
    write(*,*) 'Particles read: ', count

    call GarbageCollection(pParts)

    call Write_InitStatus('External',1)

  end subroutine initializeExternal


end module initExternal
