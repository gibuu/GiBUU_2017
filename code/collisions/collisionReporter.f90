!******************************************************************************
!****m* /collisionReporter
! NAME
! module collisionReporter
!
! PURPOSE
! This module is for statistical informations about 2-body collisions.
! It provides some 2D grids in order to store the number of collisions
! as a function of time, sqrt(s) and 2-body event code.
!
! Event codes like in  collisionNumbering::ReportEventNumber:
! *  -3: BaB    (Baryon-Antibaryon-Annihilation)
! *  -2: Manni  (Meson-Baryon-Annihilation)
! *  -1: Elastic
! *   0: == Low Energy ==
! *   1: FRITIOF
! *   2: PYTHIA
!
! INPUTS
! The JobCard "collReporter"
!
!******************************************************************************

!******************************************************************************
!****n* collisionReporter/collReporter
! NAME
! NAMELIST /collReporter/
! PURPOSE
! Namelist for collisionReporter includes:
! * UseCollReporter
! * cR_sizeT
! * cR_sizeE
! * cR_DeltaT
! * cR_DeltaE
!******************************************************************************

module collisionReporter
  implicit none
  private

  !****************************************************************************
  !****g* collisionReporter/UseCollReporter
  ! PURPOSE
  ! Enable or disable the collision reporter.
  ! SOURCE
  logical,save :: UseCollReporter = .FALSE.
  !****************************************************************************

  !****************************************************************************
  !****g* collisionReporter/cR_sizeT
  ! PURPOSE
  ! Number of timestep bins.
  ! SOURCE
  integer,save :: cR_sizeT = 200
  !****************************************************************************

  !****************************************************************************
  !****g* collisionReporter/cR_sizeE
  ! PURPOSE
  ! Number of sqrt(s) bins.
  ! SOURCE
  integer,save :: cR_sizeE = 100
  !****************************************************************************

  !****************************************************************************
  !****g* collisionReporter/cR_DeltaT
  ! PURPOSE
  ! Size of timestep bins.
  ! SOURCE
  real,save    :: cR_DeltaT = 0.1
  !****************************************************************************

  !****************************************************************************
  !****g* collisionReporter/cR_DeltaE
  ! PURPOSE
  ! Size of sqrt(s) bins.
  ! SOURCE
  real,save    :: cR_DeltaE = 0.1
  !****************************************************************************

  logical,save :: initFlag=.true.

  integer,allocatable :: ArrayN(:,:,:)
  real,allocatable    :: ArrayW(:,:,:)


  Public:: cR_Add, cR_Write


contains

  !****************************************************************************
  !****s* collisionReporter/initInput
  ! NAME
  ! subroutine initInput
  !
  ! PURPOSE
  ! Initialize the "Collision Reporter"
  !
  ! INPUTS
  ! * NAMELIST collReporter
  !
  ! OUTPUT
  ! (none)
  !****************************************************************************
  subroutine initInput

    use output

    implicit none
    integer :: ios

    NAMELIST /collReporter/ UseCollReporter,cR_sizeT,cR_sizeE,cR_DeltaT,cR_DeltaE
    call Write_ReadingInput('collReporter',0)
    rewind(5)
    read(5,nml=collReporter,ioStat=ios)
    call Write_ReadingInput('collReporter',0,ios)

    write(*,*) '  Use Collision Reporter: ',UseCollReporter
    if (UseCollReporter) then
       write(*,*) '  number timesteps   =',cR_sizeT
       write(*,*) '  number sqrt(s)-bins=',cR_sizeE
       write(*,*) '  delta  timesteps   =',cR_deltaT,' fm'
       write(*,*) '  delta  sqrt(s)-bins=',cR_deltaE,' GeV'

       write(*,*)
       write(*,*) '  ===> MAX: ',cR_sizeT*cR_deltaT,' fm *',cR_sizeE*cR_deltaE,' GeV'

       allocate(ArrayN(-3:2,1:cR_sizeT,1:cR_sizeE))
       allocate(ArrayW(-3:2,1:cR_sizeT,1:cR_sizeE))

       ArrayN = 0
       ArrayW = 0.
    end if

    call Write_ReadingInput('collReporter',1)

  end subroutine initInput

  !****************************************************************************
  !****s* collisionReporter/cR_Add
  ! NAME
  ! subroutine cR_Add(sqrtS,time,code2,weight)
  !
  ! PURPOSE
  ! add a collision that occured to the history.
  !
  ! INPUTS
  ! * real    :: sqrtS  -- sqrt(S) of the collision
  ! * real    :: time   -- time, when the collisions happened
  ! * integer :: code2  --
  ! * real    :: weight --
  !
  ! OUTPUT
  ! (none)
  !****************************************************************************
  subroutine cR_Add(sqrtS,time,code2,weight)
    implicit none
    real,    intent(in) :: sqrtS,time,weight
    integer, intent(in) :: code2

    integer iX,iY

    if (initFlag) then
       call initInput
       initFlag=.false.
    end if
    if (.not.UseCollReporter) return

    iX = nint(time/cR_deltaT)+1
    iY = nint(sqrtS/cR_deltaE)+1
    if ((iX.lt.cR_sizeT).and.(iY.lt.cR_sizeE)) then
       ArrayN(code2,iX,iY) = ArrayN(code2,iX,iY)+1
       ArrayW(code2,iX,iY) = ArrayW(code2,iX,iY)+weight
    end if
  end subroutine cR_Add


  !****************************************************************************
  !****s* collisionReporter/cR_Write
  ! NAME
  ! subroutine cR_Write(nRun)
  !
  ! PURPOSE
  ! write out the information stored about collisions
  !
  ! INPUTS
  ! * integer :: nRun -- number of Runs
  !
  ! OUTPUT
  ! some files are (re-)written.
  !****************************************************************************
  subroutine cR_Write(nRun)
    use output

    implicit none

    integer, intent(in) :: nRun

    integer :: code2,iX,iY

    real :: x,y,z,zN,f

     if (initFlag) then
       call initInput
       initFlag=.false.
    end if
    if (.not.UseCollReporter) return

    do code2=-3,2

       open(141,file='CollRep.'//trim(intToChar(10+code2))//'.dat', status='unknown')
       open(142,file='CollRep.'//trim(intToChar(20+code2))//'.dat', status='unknown')

       rewind(141)
       rewind(142)

       do iX=1,cR_sizeT
          x = (real(iX-1))*cR_deltaT

          do iY=1,cR_sizeE

             f = 1
             if (iX==1) f = f*2
             if (iY==1) f = f*2

             y = (real(iY-1))*cR_deltaE

             z = f*ArrayW(code2,ix,iY)/(cR_DeltaT*cR_DeltaE*nRun)+1e-20
             zN= f*ArrayN(code2,ix,iY)/(cR_DeltaT*cR_DeltaE*nRun)+1e-20


             write(141,'(2f9.2,1P,1e12.5)')  x,y,z
             write(142,'(2f9.2,1P,1e12.5)')  x,y,zN

          end do
          write(141,*)
          write(142,*)
       end do

       close(141)
       close(142)

    end do

  end subroutine cR_Write

  !****************************************************************************

end module collisionReporter
