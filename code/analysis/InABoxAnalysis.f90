!******************************************************************************
!****m* /InABoxAnalysis
! NAME
! module InABoxAnalysis
! PURPOSE
! This module contains analysis routines for box simulations.
!******************************************************************************
module InABoxAnalysis

  implicit none

  private

  !****************************************************************************
  !****g* InABoxAnalysis/Enable
  ! SOURCE
  !
  logical, save :: Enable = .true.
  !
  ! PURPOSE
  ! Flag to enable or disable the box analysis alltogether.
  !****************************************************************************

  !****************************************************************************
  !****g* InABoxAnalysis/Interval
  ! SOURCE
  !
  integer, save :: Interval = 20
  !
  ! PURPOSE
  ! Interval for output, i.e. number of timesteps after which output is written.
  !****************************************************************************

  public :: DoInABoxAnalysisTime

contains

  !****************************************************************************
  !****s* InABoxAnalysis/DoInABoxAnalysisTime
  ! NAME
  ! subroutine DoInABoxAnalysisTime(realPart,timestep)
  ! PURPOSE
  ! Do the box analysis at a particular timestep.
  !****************************************************************************
  subroutine DoInABoxAnalysisTime(realPart,timestep)
    use particleDefinition
    use output, only: Write_ReadingInput, intTochar4
    use histMPf90
    use histf90
    use densityModule, only: gridsize

    type(particle), dimension(:,:), intent(in), target :: realPart
    integer, intent(in) :: timestep

    logical, save :: startFlag = .true.
    integer, save :: nHist
    type(histogramMP), save :: hist
    type(histogram), save :: zDens
    real, dimension(:,:), allocatable, save :: arrQ2_all, arrQ2_formed

    type(particle), pointer :: pPart
    integer :: nEns,nPart, i,j, iID, ios
    real :: wForm, Q2, mulfak

    !**************************************************************************
    !****n* InABoxAnalysis/InABoxAnalysis
    ! NAME
    ! NAMELIST InABoxAnalysis
    ! PURPOSE
    ! Includes the input parameters:
    ! * Enable
    ! * Interval
    !**************************************************************************
    NAMELIST /InABoxAnalysis/ Enable, Interval

    if (startFlag) then

       startFlag=.false.

       call Write_ReadingInput('InABoxAnalysis',0)
       rewind(5)
       read(5,nml=InABoxAnalysis,iostat=ios)
       write(*,*) ' Enable =', Enable
       write(*,*) ' Interval =', Interval
       call Write_ReadingInput('InABoxAnalysis',1,ios)

       if (.not. Enable) return

       call CreateHistMP(hist, "dN/dp", 0., 2., 0.01, 2)

       nHist = Map2HistMP_getN(2)
       allocate(arrQ2_all(0:nHist,2))
       allocate(arrQ2_formed(0:nHist,2))

       open(123,file="Q2_all_formed.dat", status="unknown")
       call WriteHistMP_Names(hist,123)
       close(123)
       open(123,file="N_all_formed.dat", status="unknown")
       call WriteHistMP_Names(hist,123)
       close(123)

       call CreateHist(zDens, "density along z direction", -gridsize(3), gridsize(3), 1.)
    end if

    if (.not. Enable) return

    call ClearHistMP(hist)
    arrQ2_all = 0.0
    arrQ2_formed = 0.0

    call ClearHist(zDens)

    nEns = size(realPart,dim=1)
    nPart = size(realPart,dim=2)

    mulfak = 1.0/nEns

    do i=1,nEns
       do j=1,nPart
          pPart => realPart(i,j)
          if (pPart%Id <  0) exit
          if (pPart%Id <= 0) cycle

          if (pPart%in_Formation) then
             wForm = 0.0
          else
             wForm = 1.0
          end if

          call AddHistMP(hist, pPart, absMom(pPart), 1.0, wForm)

          call AddHist(zDens, pPart%position(3), 1.)

          Q2 = 2*pPart%momentum(3)**2 - pPart%momentum(1)**2 - pPart%momentum(2)**2

          arrQ2_all(0,1) = arrQ2_all(0,1) + 1.0
          arrQ2_all(0,2) = arrQ2_all(0,2) + Q2

          arrQ2_formed(0,1) = arrQ2_formed(0,1) + wForm
          arrQ2_formed(0,2) = arrQ2_formed(0,2) + Q2*wForm

          iID = Map2HistMP(pPart, 2)
          if (iID>0) then
             arrQ2_all(iID,1) = arrQ2_all(iID,1) + 1.0
             arrQ2_all(iID,2) = arrQ2_all(iID,2) + Q2

             arrQ2_formed(iID,1) = arrQ2_formed(iID,1) + wForm
             arrQ2_formed(iID,2) = arrQ2_formed(iID,2) + Q2*wForm
          end if

       end do
    end do

    if (mod(timestep,Interval)==0) then
      call WriteHistMP(hist,  file='p_all_'   //intTochar4(timestep)//'.dat', mul=mulfak, iColumn=1)
      call WriteHistMP(hist,  file='p_formed_'//intTochar4(timestep)//'.dat', mul=mulfak, iColumn=3)
      call WriteHist  (zDens, file='zDensity.'//intTochar4(timestep)//'.dat', mul=mulfak)
    end if


    open(123,file="N_all_formed.dat", status="old",position='append')
    write(123,'(i11,1P,100E12.4,0P)') timestep, &
                                      arrQ2_all(1:nHist,1)*mulfak, arrQ2_all(0,1)*mulfak, &
                                      arrQ2_formed(1:nHist,1)*mulfak, arrQ2_formed(0,1)*mulfak
    close(123)

    do i=0,nHist
       if (arrQ2_all(i,1)>0.0)    arrQ2_all(i,2)    = arrQ2_all(i,2)    / arrQ2_all(i,1)
       if (arrQ2_formed(i,1)>0.0) arrQ2_formed(i,2) = arrQ2_formed(i,2) / arrQ2_formed(i,1)
    end do

    open(123,file="Q2_all_formed.dat", status="old",position='append')
    write(123,'(i11,1P,100E12.4,0P)') timestep, &
                                      arrQ2_all(1:nHist,2), arrQ2_all(0,2),&
                                      arrQ2_formed(1:nHist,2), arrQ2_formed(0,2)
    close(123)

  end subroutine DoInABoxAnalysisTime

end module InABoxAnalysis
