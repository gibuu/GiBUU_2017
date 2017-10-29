!******************************************************************************
!****m* /hist_multipleRuns
! NAME
! module hist_multipleRuns
! PURPOSE
! Encapsulate routines and datas for 1D Histogramms, in which you want to store the results of several sequentiell runs.
!
! Features of Histograms provided by this module:
! - store paramaters of the x-binning
! - enable two y-values (y and y2)
! - track keeping of under-/over-score the given extreme values of x.
! - provide simple-to-understand output routines (cf. WriteHist)
! - provide simple histogram arithmetic (not yet implemented)
! - many others...
!
! Every Histogram prints his own multicolumn output.
! A multicolumn output of many different histograms for the same x-value
! is not implemented. This is done by the module "histMPf90".
!
! INPUTS
! ...(This module needs no input)
!******************************************************************************
module hist_multipleRuns
  use histf90, only: histogram


  integer, parameter:: NameLength = 40



  !****************************************************************************
  !****t* hist_multipleRuns/histogram_mr
  ! NAME
  ! type histogram_mr
  ! PURPOSE
  ! Type definition to store all information for a multiple-run 1D Histogram.
  ! SOURCE
  !
  type histogram_mr
     type(histogram)  :: total                          ! total result of all runs
     type(histogram)  :: thisRun                        ! result of the present run
     integer          :: numRuns                        ! number of runs
     real, dimension (:) , Allocatable :: sumOfSquares  ! Sum (Squares of the result of a single run)
  end type histogram_mr
  !****************************************************************************




contains

  !****************************************************************************
  !****s* hist_multipleRuns/CreateHist_mr
  ! NAME
  ! subroutine CreateHist_mr(H, HName,x1,x2,bin)
  ! PURPOSE
  ! This is the Constructor of a 1D-Histogram!
  ! Allocate Memory for the entries and put additional variables to their
  ! default.
  ! INPUTS
  ! * type(histogram_mr) :: H         -- Histogramm to be created
  ! * character*(*)      :: HName     -- Name of Histogram
  ! * real               :: x1,x2,bin -- Minimal/maximal value for x-coordinate
  !                                      to be considered, bin-width
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine CreateHist_mr(H, HName,x1,x2,bin)
    use histf90, only: createHist
    implicit none
    type(histogram_mr),intent(inout) :: H
    character*(*),intent(in) :: HName
    real,intent(in) :: x1
    real,intent(in) :: x2
    real,intent(in) :: bin

    call createHist(H%total, HName,x1,x2,bin)
    call createHist(H%thisRun, HName,x1,x2,bin)
    H%numRuns=0
    if (allocated(H%sumofSquares)) deallocate(H%sumOFSquares)
    allocate(H%sumOFSquares(lbound(H%total%yval,dim=1):ubound(H%total%yval,dim=1)))
    H%sumOFSquares=0.
  end subroutine CreateHist_mr

  !****************************************************************************
  !****s* hist_multipleRuns/RemoveHist_mr
  ! NAME
  ! subroutine RemoveHist_mr(H)
  ! PURPOSE
  ! Free the allocated memory
  ! INPUTS
  ! * type(histogram_mr) :: H  -- Histogramm to be used
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine RemoveHist_mr(H)
    use histf90, only: removeHist
    implicit none
    type(histogram_mr),intent(inout) :: H

    call removeHist(H%total)
    call removeHist(H%thisRun)
    H%numRuns=0
    deallocate(H%sumOfSquares)
  end subroutine RemoveHist_mr

  !****************************************************************************
  !****s* hist_multipleRuns/AddHist_mr
  ! NAME
  ! subroutine AddHist_mr(H, x,y,y2)
  ! PURPOSE
  ! Add to the given histogram at the given x-value the weight y (y2).
  ! If y2 is given, also the second column is filled.
  ! INPUTS
  ! * type(histogram_mr) :: H  -- Histogramm to be used
  ! * real               :: x  -- x-value
  ! * real               :: y  -- weight to be added
  ! * real               :: y2 -- second weight to be added [OPTIONAL]
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine AddHist_mr(H, x,y,y2)
    use histf90, only: addhist
    implicit none
    type(histogram_mr),intent(inout) :: H
    real,intent(in)          :: x
    real,intent(in)          :: y
    real,intent(in),optional :: y2

    if (present(y2)) then
       call AddHist(H%thisRun,x,y,y2)
       call AddHist(H%total,x,y,y2)
    else
       call AddHist(H%thisRun,x,y)
       call AddHist(H%total,x,y)
    end if
  end subroutine AddHist_mr


  !****************************************************************************
  !****s* hist_multipleRuns/startRunHist_mr
  ! NAME
  ! subroutine startRunHist_mr(H)
  ! PURPOSE
  ! Tells the histogram that all next data belong to a new run
  ! INPUTS
  ! * type(histogram_mr) :: H  -- Histogramm to be used
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine startRunHist_mr(H)
    use histf90, only: clearHist
    implicit none
    type(histogram_mr),intent(inout) :: H

    call ClearHist(H%thisRun)
    H%numRuns=H%numRuns+1
  end subroutine startRunHist_mr

  !****************************************************************************
  !****s* hist_multipleRuns/endRunHist_mr
  ! NAME
  ! subroutine endRunHist_mr(H)
  ! PURPOSE
  ! Tells the histogram that the analysis of a run is finished
  ! INPUTS
  ! * type(histogram_mr) :: H  -- Histogramm to be used
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine endRunHist_mr(H)
    implicit none
    type(histogram_mr),intent(inout) :: H
    H%sumOfSquares=H%sumOfSquares+H%thisRun%yval(:,1)**2
  end subroutine endRunHist_mr

  !****************************************************************************
  !****s* hist_multipleRuns/WriteHist_mr
  ! NAME
  ! subroutine WriteHist_mr(H,iFile,file)
  ! PURPOSE
  ! Write out the histogram including error analysis.
  ! INPUTS
  ! * type(histogram),intent(in)          :: H
  ! * integer,        intent(in)          :: iFile
  ! * character(*),   intent(in),optional :: file
  !
  ! NOTES
  ! The Histogram Data is not affected!!!
  !
  !****************************************************************************
  subroutine WriteHist_mr(H,iFile,file)
    implicit none
    type(histogram_mr),intent(in)          :: H
    integer,        intent(in)          :: iFile
    character(*),   intent(in),optional :: file

    real,dimension(1:3) :: S
    integer             :: iBin
    real :: error

    if (present(file)) then
       open(iFile, file=file, status='unknown')
       rewind(iFile)
    end if


1000 format (1X,'#####',/,1X,'##### Histogram ',A40,/,1X,'#####')
1010 format (1P,1X,'##### Underflow: ',2E11.4,/,&
          & 1X,'##### Entries  : ',    2E11.4,/,&
          & 1X,'##### Overflow : ',    2E11.4,/,&
          & 1X,'#####     Extrema: [',E11.4,' ... ',E11.4,' ]'/,&
          & 1X,'#####     counted: [',E11.4,' ... ',E11.4,' ]'/,&
          & 1X,'#####',&
          & 0P)
2000 format (1X,1P,E14.6,3E12.4,0P)
!2001 format (1X,1P,E14.6,4E12.4,0P)

    if (.not.allocated(H%total%yVal)) return
    if (H%numRuns.lt.1) return

    write(iFile,1000) H%total%Name


    S = SUM(h%total%yval,dim=1) - h%total%yVal(-1,1:3) - h%total%yVal(0,1:3)


    write(iFile,1010) &
         & h%total%yVal(-1,1), h%total%yVal(-1,2), &
         & S(1)        , S(2)        , &
         & h%total%yVal( 0,1), h%total%yVal( 0,2), &
         & H%total%xExtreme , &
         & H%total%xMin      , H%total%xMax

    write(iFile,'(1X,"#####",/,1X,"#####",A,I4)') ' Number of runs=',H%numruns
    write(iFile,'(1X,"#####",/,1X,"#####",A)') ' x, y(1), y(2), error of y(1)'

    ibinLoop: do iBin=1,ubound(h%total%yval,dim=1)
       if (H%numRuns.gt.1) then
          error= 1/float(H%numRuns)/float(H%numRuns-1)*(h%sumOfSquares(ibin)-&
               & h%total%yVal(iBin,1)**2/float(H%numRuns))
          if (error.gt.0) then
             error=sqrt(error)
          else if (error.lt.-0.1) then
             write(*,*) 'error^2 lt -0.1 in writeHist_mr!',h%sumOfSquares(ibin),h%total%yVal(iBin,1)**2/float(H%numRuns), ibin
             stop 'hist_multipleRuns.f90'
          else
             error=0.
          end if
       else if (H%numRuns.eq.1) then
          error= h%total%yVal(iBin,1)
       else
          stop 'Error in writehist_mr'
       end if
       write(iFile,2000) &
            & H%total%xMin+(real(iBin)-0.5)*H%total%xBin,&
            & h%total%yVal(iBin,1)/H%total%xBin/float(H%numRuns),&
            & h%total%yVal(iBin,2), error/H%total%xBin
    end do ibinLoop

    if (present(file)) then
       close(iFile)
    end if

  end subroutine WriteHist_mr

end module hist_multipleRuns
