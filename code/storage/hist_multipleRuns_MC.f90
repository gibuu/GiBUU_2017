!******************************************************************************
!****m* /hist_multipleRuns_MC
! NAME
! module hist_multipleRuns_MC
!
! PURPOSE
! * Encapsulate routines and datas for 1D Histogramms, in which you want to store the results of several sequentiell runs.
! * Superior to hist_multipleRuns since it offers the possibility to have multiple columns in the histograms
!
! Features of Histograms provided by this module:
! - store paramaters of the x-binning
! - enable  multiple y-values
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
module hist_multipleRuns_MC
  use histMC, only: histogramMC
  implicit none

  integer, parameter:: NameLength = 40


  !****************************************************************************
  !****t* hist_multipleRuns_MC/histogram_mr_mc
  ! NAME
  ! type histogram_mr_mc
  ! PURPOSE
  ! Type definition to store all information for a multiple-run 1D Histogram.
  ! SOURCE
  !
  type histogram_mr_mc
     type(histogramMC)  :: total                          ! total result of all runs
     type(histogramMC)  :: thisRun                        ! result of the present run
     integer          :: numRuns                        ! number of runs
     real, dimension (:,:) , Allocatable :: sumOfSquares  ! Sum (Squares of the result of a single run)
  end type histogram_mr_mc
  !****************************************************************************


contains

  !****************************************************************************
  !****s* hist_multipleRuns_MC/CreateHist_mr_mc
  ! NAME
  ! subroutine CreateHist_mr_mc(H, HName,x1,x2,bin,nCH)
  ! PURPOSE
  ! This is the Constructor of a 1D-Histogram!
  ! Allocate Memory for the entries and put additional variables to their
  ! default.
  ! INPUTS
  ! * type(histogram_mr_mc) :: H         -- Histogramm to be created
  ! * character*(*)      :: HName     -- Name of Histogram
  ! * real               :: x1,x2,bin -- Minimal/maximal value for x-coordinate
  !                                      to be considered, bin-width
  ! * integer, intent(in) :: nCH      -- number of channels
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine CreateHist_mr_mc(H, HName,x1,x2,bin,nCh)
    use histMC, only: createHistMC

    type(histogram_mr_mc),intent(inout) :: H
    character*(*),intent(in) :: HName
    real,intent(in) :: x1
    real,intent(in) :: x2
    real,intent(in) :: bin
    integer, intent(in) :: nCH

    call createHistMC(H%total  , HName,x1,x2,bin,nCH)
    call createHistMC(H%thisRun, HName,x1,x2,bin,nCH)
    H%numRuns=0
    if (allocated(H%sumofSquares)) deallocate(H%sumOFSquares)
    allocate(H%sumOFSquares(lbound(H%total%yval,dim=1):ubound(H%total%yval,dim=1),1:nCH))
    H%sumOFSquares=0.
  end subroutine CreateHist_mr_mc

  !****************************************************************************
  !****s* hist_multipleRuns_MC/RemoveHist_mr_mc
  ! NAME
  ! subroutine RemoveHist_mr_mc(H)
  ! PURPOSE
  ! Free the allocated memory
  ! INPUTS
  ! * type(histogram_mr_mc) :: H  -- Histogramm to be used
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine RemoveHist_mr_MC(H)
    use histMC, only: removeHistMC
    type(histogram_mr_mc),intent(inout) :: H

    call removeHistMC(H%total)
    call removeHistMC(H%thisRun)
    H%numRuns=0
    deallocate(H%sumOfSquares)
  end subroutine RemoveHist_mr_MC

  !****************************************************************************
  !****s* hist_multipleRuns_MC/AddHist_mr_mc
  ! NAME
  ! subroutine AddHist_mr_mc(H, x,nCH,y)
  ! PURPOSE
  ! Add to the given histogram at the given x-value the weight y (y2).
  ! If y2 is given, also the second column is filled.
  ! INPUTS
  ! * type(histogram_mr_mc) :: H  -- Histogramm to be used
  ! * real               :: x  -- x-value
  ! * real               :: y  -- weight to be added
  ! * integer, intent(in) :: nCH  -- channel
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine AddHist_mr_mc(H, x,nCh,y)
    use histMC, only: addhistMC

    type(histogram_mr_mc),intent(inout) :: H
    real,intent(in)          :: x
    real,intent(in)          :: y
    integer, intent(in)      :: nCH

    call AddHistMC(H%thisRun,x,nCH,y)
    call AddHistMC(H%total  ,x,nCH,y)
  end subroutine AddHist_mr_mc


  !****************************************************************************
  !****s* hist_multipleRuns_MC/startRunHist_mr_mc
  ! NAME
  ! subroutine startRunHist_mr_mc(H)
  ! PURPOSE
  ! Tells the histogram that all next data belong to a new run
  ! INPUTS
  ! * type(histogram_mr_mc) :: H  -- Histogramm to be used
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine startRunHist_mr_mc(H)
    use histMC, only: clearHistMC
    type(histogram_mr_mc),intent(inout) :: H

    call ClearHistMC(H%thisRun)
    H%numRuns=H%numRuns+1
  end subroutine startRunHist_mr_mc

  !****************************************************************************
  !****s* hist_multipleRuns_MC/endRunHist_mr_mc
  ! NAME
  ! subroutine endRunHist_mr_mc(H)
  ! PURPOSE
  ! Tells the histogram that the analysis of a run is finished
  ! INPUTS
  ! * type(histogram_mr_mc) :: H  -- Histogramm to be used
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine endRunHist_mr_mc(H)
    type(histogram_mr_mc),intent(inout) :: H
    integer :: i
    do i=lbound(H%thisRun%yval,dim=1),ubound(H%thisRun%yval,dim=1)
       H%sumOfSquares(i,:)=H%sumOfSquares(i,:)+H%thisRun%yval(i,:)**2
    end do
  end subroutine endRunHist_mr_mc

  !****************************************************************************
  !****s* hist_multipleRuns_MC/WriteHist_mr_mc
  ! NAME
  ! subroutine WriteHist_mr_mc(H,iFile,file)
  ! PURPOSE
  ! Write out the histogram including error analysis.
  ! INPUTS
  ! * type(histogram_mr_mc),intent(in)          :: H
  ! * integer,              intent(in)          :: iFile
  ! * character(*),         intent(in),optional :: file
  !
  ! NOTES
  ! The Histogram Data is not affected!!!
  !
  !****************************************************************************
  subroutine WriteHist_mr_mc(H,iFile,file)
    use histMC

    type(histogram_mr_mc),intent(in)          :: H
    integer,              intent(in)          :: iFile
    character(*),         intent(in),optional :: file

    !real,dimension(1:3) :: S
    integer             :: iBin,nCH,j
    real ::  fak
    real , dimension (:), allocatable :: error
    character(30) :: f
    logical :: nonempty


    if (present(file)) then
       open(iFile, file=file, status='unknown')
       rewind(iFile)
    end if

    if (.not.allocated(H%total%yVal)) return
    if (H%numRuns.lt.1) return

    fak = 1./H%total%xBin

    open(iFile, file=file, status='unknown')
    rewind(iFile)

    nCh = size(H%total%yVal,dim=2)

    nonempty = WriteHeaderMC(H%total,iFile,fak,.false.)

    do ibin=1,nCh
       write(iFile,'(a,i2,a,i2,a,a)') ' ##### Column ',nch+2+iBin,&
            & ' (Error of Y-Channel ',iBin,'): Error of ',H%total%yDesc(iBin)
    end do
       write(iFile,'(A)') ' #####'

    ! write histogram data
    write(f,'(a,i2,a)') '(',2*nCh+2,'E12.4)'

    allocate(error(1:ubound(h%total%yval,dim=2)))
    ibinLoop: do iBin=1,ubound(h%total%yval,dim=1)
       if (H%numRuns.gt.1) then
          error= 1/float(H%numRuns)/float(H%numRuns-1)*(h%sumOfSquares(ibin,:)-&
               & h%total%yVal(iBin,:)**2/float(H%numRuns))
          do j=lbound(error,dim=1),ubound(error,dim=1)
             if (error(j).gt.0) then
                error(j)=sqrt(error(j))
             else if (error(j).lt.-0.1) then
                write(*,*) 'error^2 lt -0.1 in writeHist_mr!',h%sumOfSquares(ibin,:),h%total%yVal(iBin,:)**2/float(H%numRuns), ibin
                stop 'hist_multipleRuns.f90'
             else
                error(j)=0.
             end if
          end do
       else if (H%numRuns.eq.1) then
          error= h%total%yVal(iBin,:)
       else
          stop 'Error in writehist_mr_mc'
       end if
       write(iFile,fmt=f) H%total%xMin+(real(iBin)-0.5)*H%total%xBin, &
            sum(H%total%yVal(iBin,:))*fak/float(H%numRuns), &
            H%total%yVal(iBin,:)*fak/float(H%numRuns), error*fak
    end do ibinLoop

    if (present(file)) then
       close(iFile)
    end if

  end subroutine WriteHist_mr_mc

end module hist_multipleRuns_MC
