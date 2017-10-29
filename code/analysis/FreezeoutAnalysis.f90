!******************************************************************************
!****m* /FreezeoutAnalysis
! NAME
! module FreezeouAnalysis
!
! PURPOSE
! This module provides routines to do a 'freeze out' analysis, i.e. it yields
! access to the position of particles at their last interaction. Particles
! being subject to some potentials have some additional freeze out condition,
! e.g. when the baryon density drops to the value 0.2/fm^3.
!******************************************************************************
module FreezeoutAnalysis

  implicit none
  private

  public :: getFreezeoutAnalysis_Pert
  public :: getFreezeoutAnalysis_Real
  public :: DoFreezeoutAnalysisPerTime

  !****************************************************************************
  !****g* FreezeoutAnalysis/FreezeoutAnalysis_Pert
  ! SOURCE
  !
  logical, save :: FreezeoutAnalysis_Pert = .false.
  !
  ! PURPOSE
  ! Flag to do freeze out analysis for perturbative particles
  !****************************************************************************

  !****************************************************************************
  !****g* FreezeoutAnalysis/FreezeoutAnalysis_Real
  ! SOURCE
  !
  logical, save :: FreezeoutAnalysis_Real = .false.
  !
  ! PURPOSE
  ! Flag to do freeze out analysis for real particles
  !****************************************************************************

  !****************************************************************************
  !****g* FreezeoutAnalysis/potThreshold
  ! SOURCE
  !
  real, save :: potThreshold = 0.005
  !
  ! PURPOSE
  ! threshold value in GeV. If the absolute value of the potential is below
  ! this value, the particle is considered to be 'free', e.g. it 'escaped'
  !****************************************************************************

  logical, save :: init = .true.

contains

  !****************************************************************************
  !****f* FreezeoutAnalysis/getFreezeoutAnalysis_Pert
  ! NAME
  ! logical function getFreezeoutAnalysis_Pert()
  ! PURPOSE
  ! return the value of FreezeoutAnalysis_Pert
  !****************************************************************************
  logical function getFreezeoutAnalysis_Pert()
    if (init) call initInput
    getFreezeoutAnalysis_Pert = FreezeoutAnalysis_Pert
  end function getFreezeoutAnalysis_Pert

  !****************************************************************************
  !****f* FreezeoutAnalysis/getFreezeoutAnalysis_Real
  ! NAME
  ! logical function getFreezeoutAnalysis_Real()
  ! PURPOSE
  ! return the value of FreezeoutAnalysis_Real
  !****************************************************************************
  logical function getFreezeoutAnalysis_Real()
    if (init) call initInput
    getFreezeoutAnalysis_Real = FreezeoutAnalysis_Real
  end function getFreezeoutAnalysis_Real

  !****************************************************************************
  !****s* FreezeoutAnalysis/initInput
  ! NAME
  ! subroutine initInput
  ! PURPOSE
  ! Read namelist 'Freezeout' from jobcard.
  !****************************************************************************
  subroutine initInput
    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n* FreezeoutAnalysis/Freezeout
    ! NAME
    ! namelist /Freezeout/
    ! PURPOSE
    ! Namelist for FreezeoutAnalysis includes:
    ! * FreezeoutAnalysis_Pert
    ! * FreezeoutAnalysis_Real
    ! * potThreshold
    !**************************************************************************
    namelist /Freezeout/ FreezeoutAnalysis_Pert, FreezeoutAnalysis_Real, &
         potThreshold

    integer :: ios

    call Write_ReadingInput('Freezeout',0)
    rewind(5)
    read(5,nml=Freezeout,iostat=ios)

    call Write_ReadingInput('Freezeout',0,ios)

    write(*,*) 'Freezeout analysis (real,pert): ', &
               FreezeoutAnalysis_Real, FreezeoutAnalysis_Pert
    write(*,*) 'threshold potential: ',potThreshold,' GeV'
    call Write_ReadingInput('Freezeout',1)

    init = .false.

  end subroutine initInput

  !****************************************************************************
  !****s* FreezeoutAnalysis/DoFreezeoutAnalysisPerTime
  ! NAME
  ! subroutine DoFreezeoutAnalysisPerTime(iTime, Time, pertPart, realPart)
  ! PURPOSE
  ! Do the analysis after each time step
  !****************************************************************************
  subroutine DoFreezeoutAnalysisPerTime(iTime, Time, pertPart, realPart)
    use particleDefinition
    use particleProperties
    use PIL_freezeout, only: PIL_freezeout_PUT, PIL_freezeout_GET
    use potentialModule, only: potential_LRF

    integer, intent(in)       :: iTime
    real,    intent(in)       :: Time
    type(particle), intent(in), dimension(:,:),target :: realPart, pertPart

    integer :: i,j
    type(particle), POINTER :: pPart

    real, dimension(0:3) :: pos, posOld
    integer :: history, historyOld
    logical :: escaped, escapedOld
    logical :: changed
    real :: pot

    if (init) call initInput

    if (FreezeoutAnalysis_Real) then
       write(*,*) 'FreezeoutAnalysis_Real not yet implemented.'
       stop
    end if

    if (FreezeoutAnalysis_Pert) then

       do j=1,Size(pertPart,dim=2)
          do i=1,Size(pertPart,dim=1)
             pPart => pertPart(i,j)
             if (pPart%Id <  0) exit
             if (pPart%Id <= 0) cycle

             if (hadron(pPart%ID)%stability .ne. 0) cycle

             pos(0) = Time
             pos(1:3) = pPart%position
             history = pPart%history

             pot = potential_LRF(pPart)
!             escaped = .true. ! dummy !!!
             escaped = (abs(pot)<potThreshold)

             changed = .not.PIL_freezeout_GET(pPart%number,posOld,historyOld,escapedOld)

             if (.not.changed) then ! i.e. it is already in the list
!                write(*,*) 'got: ',pPart%number,posOld,historyOld
                if (history .ne. historyOld) changed = .true.
                if (escaped .neqv. escapedOld) changed = .true.
                if (.not.escaped) changed = .true.
             end if

             if (changed) then
!                write(*,*) 'set: ',pPart%number,pos,history
                call PIL_freezeout_PUT(pPart%number,pos,history,escaped)
             end if

          end do
       end do
    end if

  end subroutine DoFreezeoutAnalysisPerTime

end module FreezeoutAnalysis
