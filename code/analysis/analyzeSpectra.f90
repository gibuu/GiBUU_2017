!******************************************************************************
!****m* /analyzeSpectra
! NAME
! module analyzeSpectra
! PURPOSE
! generate spectra of different particle species, independent of event type
!******************************************************************************
module analyzeSpectra

  use histf90

  implicit none

  private

  public :: doAnalyzeSpectra

  !****************************************************************************
  !****g* analyzeSpectra/realID
  ! SOURCE
  logical, dimension(1:122), save :: realID = .false.
  ! PURPOSE
  ! Switch on/off the output for specific particle IDs from the real particles
  ! vector
  !****************************************************************************

  !****************************************************************************
  !****g* analyzeSpectra/pertID
  ! SOURCE
  logical, dimension(1:122), save :: pertID = .false.
  ! PURPOSE
  ! Switch on/off the output for specific particle IDs from the pert particles
  ! vector
  !****************************************************************************

  logical, save :: anyReal = .false.
  logical, save :: anyPert = .false.
  logical, save :: initFlag= .true.

contains

  !****************************************************************************
  !****s* analyzeSpectra/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input out of jobcard
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput
    use callstack, only: traceBack

    !**************************************************************************
    !****n* analyzeSpectra/AnalyzeSpectra
    ! NAME
    ! NAMELIST AnalyzeSpectra
    ! PURPOSE
    ! Includes the input parameters:
    ! * realID
    ! * pertID
    !**************************************************************************
    NAMELIST /AnalyzeSpectra/ realID, pertID

    integer :: ios, ID
    character(3),dimension(1:2) :: S

    call Write_ReadingInput('AnalyzeSpectra',0)
    rewind(5)
    read(5,nml=AnalyzeSpectra,iostat=ios)

    anyReal = any(realID)
    anyPert = any(pertID)

    if (anyReal.or.anyPert) then
       do ID=1,122
          if (realID(ID) .or. pertID(ID)) then
             if (realID(ID)) then
                S(1) = 'yes'
             else
                S(1) = '---'
             end if
             if (pertID(ID)) then
                S(2) = 'yes'
             else
                S(2) = '---'
             end if

             write(*,'("  ",i3,"  ",A3,"   ",A3)') ID, S(1), S(2)
          end if
       end do
    end if

    call Write_ReadingInput('AnalyzeSpectra',0,ios)
    call Write_ReadingInput('AnalyzeSpectra',1)

    initFlag = .false.

  end subroutine readInput


  subroutine doAnalyzeSpectra(realPart, pertPart, timestep, time, finalFlag)
    use particleDefinition
    use histMC
    use output, only: intToChar4

    type(particle), target, dimension(:,:), intent(in) :: realPart, pertPart
    integer, intent(in) :: timestep
    real, intent(in) :: time
    logical, intent(in) :: finalFlag

    type(histogramMC) :: hists
    integer :: i,j, nEns, nPart, ID
    type(particle), POINTER :: pPart
    real :: mom

    if (initFlag) call readInput
    if (.not.(anyReal.or.anyPert)) return

    call CreateHistMC(hists, "momentum", 0.0, 1.0, 0.01, 122)

    nEns  = size(realPart,dim=1)
    nPart = size(realPart,dim=2)

    if (anyReal) then
       do i=1,size(realPart,dim=1)
          do j=1,size(realPart,dim=2)
             pPart => realPart(i,j)
             if (pPart%Id <  0) exit
             if (pPart%Id <= 0) cycle

             if (realID(pPart%Id)) then
                mom = absMom(pPart)
                call AddHistMC(hists, mom, pPart%Id, 1.0)
             end if


          end do
       end do

       if (finalFlag) then
          call WriteHistMC(hists,"realSpectra_final.dat", mul = 1.0/nEns)
       else
          call WriteHistMC(hists,"realSpectra_"//intTochar4(timestep)//".dat", mul = 1.0/nEns)
       end if

!!$       do ID=1,122
!!$          if (realID(ID)) then
!!$          end if
!!$       end do

    end if


    if (anyPert) then
       if (anyReal) call ClearHistMC(hists)

       do i=1,size(pertPart,dim=1)
          do j=1,size(pertPart,dim=2)
             pPart => pertPart(i,j)
             if (pPart%Id <  0) exit
             if (pPart%Id <= 0) cycle

             if (pertID(pPart%Id)) then
!                mom = absMom(pPart)
                mom = pPart%momentum(0)
                call AddHistMC(hists, mom, pPart%Id, 1.0)
             end if

          end do
       end do

       if (finalFlag) then
          call WriteHistMC(hists,"pertSpectra_final.dat", mul = 1.0/nEns)
       else
          call WriteHistMC(hists,"pertSpectra_"//intTochar4(timestep)//".dat", mul = 1.0/nEns)
       end if

    end if

  end subroutine doAnalyzeSpectra

end module analyzeSpectra
