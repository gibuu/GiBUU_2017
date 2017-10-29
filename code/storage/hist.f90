!******************************************************************************
!****m* /histf90
! NAME
! module histf90
! PURPOSE
! Encapsulate all routines and datas for 1D Histograms.
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
module histf90

  implicit none
  private

  integer, parameter, public :: NameLength = 80


  !****************************************************************************
  !****t* histf90/histogram
  ! NAME
  ! type histogram
  ! PURPOSE
  ! Type definition to store all information for a 1D Histogram.
  ! SOURCE
  !
  type, public :: histogram
     real             :: xMin        ! smallest x-value
     real             :: xMax        ! largest x-value
     real             :: xBin        ! width of x-Bin
     real             :: xExtreme(2) ! extremes of used x-values
                                     ! (min,max)
     character*(NameLength) :: Name  ! name to be written
     real,allocatable :: yVal(:,:)   ! histogram values:
                                     ! (x), (yy1,yy2,yy3...)
                                     ! bin 0,-1 : extreme values !!!
     logical :: initialized = .false.! flag to indicate whether the histogram has been initialized
  end type histogram
  !****************************************************************************

  public :: CreateHist, RemoveHist, ClearHist, Addhist
  public :: WriteHist, ReadHist, FetchHist
  public :: WriteHist_Spline, WriteHist_BSpline
  public :: WriteHist_Gauss, WriteHist_Novo
  public :: WriteHist_Integrated
  public :: sumHist


  ! some format specifiers for writing out histograms
  character(60), parameter :: fmt1020 = "('### ',A,', #SubPoints=',i5,/,'###')"
  character(60), parameter :: fmt2000 = "(ES14.6,3ES12.4)"
  character(60), parameter :: fmt2001 = "(ES14.6,4ES12.4)"

  !****************************************************************************
  !****s* histf90/AddHist
  ! NAME
  ! subroutine AddHist(H, x,y,y2)
  !
  ! subroutine AddHist(H1,H2, x,y,y2)
  ! PURPOSE
  ! Add to the given histogram(s) at the given x-value the weight y (y2).
  ! If y2 is given, also the second column is filled.
  ! INPUTS
  ! * type(histogram) :: H  -- Histogram to be used
  ! or:
  ! * type(histogram) :: H1,H2  --  Histograms to be used
  ! * real            :: x  -- x-value
  ! * real            :: y  -- weight to be added
  ! * real            :: y2 -- second weight to be added [OPTIONAL]
  ! OUTPUT
  ! H is changed
  ! NOTES
  ! This routine is overloaded: Giving two histograms at the same time is
  ! like calling the routine twice.
  !
  ! This overloading is useful, if you want to fill some histogram in an
  ! array, but also keep the sum of all histograms of this array in the
  ! histogram at position 0: call AddHist( (/hist(0), hist(iH)/), ...)
  !****************************************************************************
  interface AddHist
     module procedure AddHist1,AddHist2
  end interface

contains

  !****************************************************************************
  !****s* histf90/ClearHist
  ! NAME
  ! subroutine ClearHist(H)
  ! PURPOSE
  ! Sets the histogram to zero again
  ! INPUTS
  ! * type(histogram) :: H         -- Histogram to be cleared
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine ClearHist(H)

    type(histogram),intent(inout) :: H

    H%xExtreme(1:2) = (/ 99e9, -99e9 /)
    H%yVal = 0.

  end subroutine ClearHist


  !****************************************************************************
  !****s* histf90/CreateHist
  ! NAME
  ! subroutine CreateHist (H, HName, xmin, xmax, bin)
  ! PURPOSE
  ! This is the constructor of a 1D-Histogram! It allocates memory for the
  ! entries and puts additional variables to their default.
  ! INPUTS
  ! * type(histogram) :: H          -- Histogram to be created
  ! * character*(*)   :: HName      -- Name of Histogram
  ! * real            :: xmin,xmax  -- Minimal/maximal value for x-coordinate to be considered
  ! * real            :: bin        -- bin width
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine CreateHist (H, HName, xmin, xmax, bin)
    use CALLSTACK, only: TRACEBACK

    type(histogram),intent(inout) :: H
    character*(*),intent(in) :: HName
    real,intent(in) :: xmin, xmax, bin

    integer :: L

    if (xmax<xmin) then
      write(*,*) 'Error: xmax < xmin ! '
      write(*,*) xmax,xmin,HName
      call TRACEBACK ()
    end if
    if (bin<=0.) then
      write(*,*) 'Error: bin width <= 0 ! '
      write(*,*) bin,HName
      call TRACEBACK ()
    end if

    H%xMin = xmin
    H%xMax = xmax
    H%xBin = bin
    H%xExtreme(1:2) = (/ 99e9, -99e9 /)

    if (len(HName) > NameLength) then
       H%Name = HName(1:NameLength)
    else
       H%Name = HName
    end if

    if (allocated(H%yVal)) deallocate(H%yVal)

    L = ceiling((xmax-xmin)/bin)
    allocate(H%yVal(-1:L,3))

    H%yVal = 0.
    H%initialized = .true.

  end subroutine CreateHist


  !****************************************************************************
  !****s* histf90/RemoveHist
  ! NAME
  ! subroutine RemoveHist(H)
  ! PURPOSE
  ! Free the allocated memory
  ! INPUTS
  ! * type(histogram) :: H  -- Histogram to be used
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine RemoveHist(H)

    type(histogram),intent(inout) :: H

    if (allocated(H%yVal)) deallocate(H%yVal)
    H%initialized = .false.

  end subroutine RemoveHist


  !****************************************************************************
  ! cf. interface AddHist
  !****************************************************************************
  subroutine AddHist1 (H, x, y, y2)

    type(histogram), intent(inout) :: H
    real, intent(in)               :: x, y
    real, intent(in), optional     :: y2

    integer :: iBin
    real :: yy

    yy = 0.
    if (present(y2)) yy = y2

    if (x<H%xExtreme(1)) H%xExtreme(1)=x
    if (x>H%xExtreme(2)) H%xExtreme(2)=x

    if (x>=H%xMin .and. x<H%xMax) then
       iBin = int( (x-H%xMin)/H%xBin )+1
    else if (x < H%xMin) then
       iBin = -1
    else
       iBin = 0
    end if

    H%yVal(iBin,1) = H%yVal(iBin,1)+y
    H%yVal(iBin,2) = H%yVal(iBin,2)+1.
    H%yVal(iBin,3) = H%yVal(iBin,3)+yy

  end subroutine AddHist1
  !---------------------------------------------------------------------------
  subroutine AddHist2 (H1, H2, x ,y , y2)

    type(histogram), intent(inout) :: H1,H2
    real, intent(in)               :: x, y
    real, intent(in), optional     :: y2

    if (present(y2)) then
       call AddHist1(H1, x,y,y2)
       call AddHist1(H2, x,y,y2)
    else
       call AddHist1(H1, x,y)
       call AddHist1(H2, x,y)
    end if

  end subroutine AddHist2


  !****************************************************************************
  !****s* histf90/WriteHeader
  ! NAME
  ! subroutine WriteHeader(H,iFile,mul)
  !
  ! PURPOSE
  ! A header is written:
  ! * Underflow: sum of all calls with a x-value less than x-min
  ! * Entries: sum of all entries listed in section "Data" below
  ! * Overflow: sum of all calls with a x-value above x-max
  ! * Average: average value of all bins
  ! * Extrema: the smallest/biggest x-values which had been added
  ! * counted: the range of x-values, which is considered in "Entries"
  !
  ! Summing "Underflow"+"Entries"+"Overflow" gives the number of ALL calls,
  ! i.e. the integral from -infty upto +infty.
  !
  ! Output of "Underflow","Entries","Overflow" is split into 3 columns,
  ! corresponding to the columns 2)..4) of the Data-Block, see below.
  !
  ! INPUTS
  ! * type(histogramMC) :: H -- Histogram
  ! * integer :: iFile -- file number to write to
  ! * real, optional :: mul
  ! OUTPUT
  ! write to file
  !
  !****************************************************************************
  subroutine WriteHeader (H, iFile, mul)
    type(histogram), intent(in) :: H
    integer, intent(in) :: iFile
    real, intent(in), optional :: mul

    real :: mulFak,x
    real,dimension(1:3) :: S,avg
    integer :: i,iBin,iBinMax

    character(60), parameter :: fmt1000 = "('###',/,'### Histogram: ',A,/,'###')"
    character(300),parameter :: fmt1010 = "('### Underflow: ',3ES11.4,/, &
                                           &'### Entries  : ',3ES11.4,/, &
                                           &'### Overflow : ',3ES11.4,/, &
                                           &'### Average  : ',3ES11.4,/, &
                                           &'###',/, &
                                           &'###     Extrema: [',ES11.4,' ... ',ES11.4,' ]'/, &
                                           &'###     counted: [',ES11.4,' ... ',ES11.4,' ]'/, &
                                           &'###')"

    mulFak = 1.0
    if (present(mul)) mulFak = mul

    ! write histogram name
    write(iFile,fmt1000) H%Name

    ! calculate total integrated entries
    iBinMax = ubound(H%yVal,dim=1)
    S = SUM(H%yVal(1:iBinMax,:),dim=1)

    ! calculate average
    avg = 0.
    do iBin=1,iBinMax
       x = H%xMin+(real(iBin)-0.5)*H%xBin
       avg = avg + H%yVal(iBin,:)*x
    end do
    do i=1,3
       if (S(i)>0.) then
          avg(i) = avg(i)/S(i)
       else
          avg(i) = 0.
       end if
    end do

    ! write header
    write(iFile,fmt1010) H%yVal(-1,1)*mulFak, H%yVal(-1,2), H%yVal(-1,3)*mulFak, &
                         S(1)*mulFak,         S(2),         S(3)*mulFak,         &
                         H%yVal( 0,1)*mulFak, H%yVal( 0,2), H%yVal( 0,3)*mulFak, &
                         avg(1:3), &
                         H%xExtreme, &
                         H%xMin, H%xMax

  end subroutine


  !****************************************************************************
  !****s* histf90/WriteHist
  ! NAME
  ! subroutine WriteHist (H, iFile, add, mul, DoAve, maxVal, H2, file, dump)
  ! PURPOSE
  ! Write out the histogram.
  !
  ! The entries are multiplied by 'mul' and 'add' is added.
  ! INPUTS
  ! * type(histogram) :: H     -- Histogram to be used
  ! * integer         :: iFile -- File number output to redirect [OPTIONAL]
  ! * real            :: add   -- factor to add      [OPTIONAL]
  ! * real            :: mul   -- factor to multiply [OPTIONAL]
  ! * logical         :: DoAve -- write also "average" (cf Notes) [OPTIONAL]
  ! * real            :: maxVal -- value used instead "divide by zero" [OPTIONAL]
  ! * type(histogram) :: H2    -- Histogram to divide by [OPTIONAL]
  ! * character*(*)   :: file  -- name of file to open and close [OPTIONAL]
  ! * logical         :: dump  -- if true, also dump it binary to file [OPTIONAL]
  ! OUTPUT
  ! Write to file number 'iFile'.
  !
  ! 4 columns are written in the data section:
  ! * 1) x-value (i.e. middle of bin)
  ! * 2) y-value
  ! * 3) number of entries that lead to y-value
  ! * 4) y2-value (if used, otherwise 0)
  !
  ! Columns 2)..4) are divided by the bin-width.
  !
  ! If the (optional) parameter "add" was given, this value is added
  ! to the written values of columns 2)...4).
  ! (E.g. this is used, if one wants to create logplots with xmgr/grace:
  ! the value "0.0" is not allowed as input and destroys everything.
  ! Therefore calling this routine with the argument "add=1e-20" prohibits
  ! writing of "0.0"-values and everything is fine.)
  !
  ! If the (optional) parameter "mul" was given, all written values
  ! of columns 2)...4) are multiplied by this value.
  ! (E.g. this is used, if one wants to divide the output by a
  ! "number of runs" value: "mul=1./NumberOfRuns".)
  !
  ! Column 4) "y2" provides you with a simple way, to not calculate only
  ! one histogram "y", but simultanousely a very closed connected one.
  ! (E.g. if you want to calculate ratios "y2"/"y")
  !
  ! If DoAve is given:
  ! Write out the histogram, but now also divide column 3 by column 1.
  ! If column 2 is zero, i.e. no entries, then parameter maxVal is
  ! written instead.
  ! This column output is NOT divided by bin widths.
  !
  ! Plotting programs normally neglect the header lines and
  ! columns 3) and 4) -- leaving us with the use of these routines
  ! as expected.
  !
  ! If the (optional) parameter H2 is given, the first all entries of the
  ! histogram are divided by the entries of histogram H2, column 1.
  !
  ! If the (optional) parameter 'file' is given, this routine first opens
  ! the file with this name (using stream number iFile), rewinds it. The
  ! file is closed at the end of the routine.
  !
  ! NOTES
  ! The histogram data is not affected!!!
  !****************************************************************************
  subroutine WriteHist (H, iFile_in, add, mul, DoAve, maxVal, H2, file, dump)

    type(histogram), intent(in)           :: H
    integer,         intent(in), optional :: iFile_in
    real,            intent(in), optional :: add, mul
    logical,         intent(in), optional :: DoAve
    real,            intent(in), optional :: maxVal
    type(histogram), intent(in), optional :: H2
    character*(*),   intent(in), optional :: file
    logical,         intent(in), optional :: dump

    integer             :: iBin, iBinMax, iFile
    real :: addFak, mulFak, Z, maxZ
    logical :: writeZ
    real, allocatable :: yVal(:,:)

    addFak = 0.
    mulFak = 1.
    maxZ = 99.0
    writeZ = .false.
    iFile = 65

    if (present(add)) addFak = add
    if (present(mul)) mulFak = mul
    if (present(maxVal)) maxZ = maxVal
    if (present(DoAve)) writeZ = DoAve
    if (present(iFile_in)) iFile = iFile_in

    if (.not.allocated(H%yVal)) return

    if (present(file)) then
       open(iFile, file=file, status='unknown')
       rewind(iFile)
    end if

    iBinMax = ubound(H%yVal,dim=1)

    allocate(yVal(-1:iBinMax,3))
    yVal = H%yVal

    if (present(H2)) then
       if (ubound(H2%yVal,dim=1).ne.iBinMax) then
          write(*,*) 'WriteHist: ERROR, dim not equal! STOP'
          stop
       end if

       do iBin=-1,iBinMax
          if (H2%yVal(iBin,1) .eq. 0.0) then
             if (yVal(iBin,1) .ne. 0.0) yVal(iBin,1) = maxZ
             if (yVal(iBin,3) .ne. 0.0) yVal(iBin,3) = maxZ
          else
             yVal(iBin,1) = yVal(iBin,1)/H2%yVal(iBin,1)
             yVal(iBin,3) = yVal(iBin,3)/H2%yVal(iBin,1) ! we divide by column 1!
          end if

       end do
    end if

    call WriteHeader(H,iFile,mulFak)

    if (writeZ) then
       do iBin=1,iBinMax
          if (yVal(iBin,2) > 0) then
            if (yVal(iBin,1) .eq. 0.0) then
               z = maxz
            else
               z = yVal(iBin,3)/yVal(iBin,1)
            end if
         else
            z = maxz
         end if

          write(iFile,fmt2001) &
            & H%xMin+(real(iBin)-0.5)*H%xBin,&
            & z, &
            & yVal(iBin,1)/H%xBin*mulFak+addFak,&
            & yVal(iBin,2),&
            & yVal(iBin,3)/H%xBin*mulFak+addFak
       end do
    else
       write(iFile,fmt2000) (&
            & H%xMin+(real(iBin)-0.5)*H%xBin,&
            & yVal(iBin,1)/H%xBin*mulFak+addFak,&
            & yVal(iBin,2),&
            & yVal(iBin,3)/H%xBin*mulFak+addFak , iBin=1,iBinMax)
    end if

    if (present(file)) close(iFile)

    if (present(dump)) then
       if (dump) then
          if (present(file)) then
             call DumpHist(H,file//".bin",iFile,addFak,mulFak)
          else
             call DumpHist(H,"DumpHist.bin",iFile,addFak,mulFak)
          end if
       end if
    end if

  end subroutine WriteHist


  !****************************************************************************
  !****s* histf90/WriteHist_Integrated
  ! NAME
  ! subroutine WriteHist_Integrated (H, file, backward, normalize)
  ! PURPOSE
  ! Write out the histogram, integrating the data over x.
  ! INPUTS
  ! * type(histogram) :: H         -- Histogram to be used
  ! * character*(*)   :: file      -- File name to write to
  ! * logical         :: backward  -- integrate from front or back?
  ! * logical         :: normalize -- normalize by total integral
  ! OUTPUT
  ! write to file
  ! NOTES
  ! The Histogram Data is not affected!!!
  !****************************************************************************
  subroutine WriteHist_Integrated (H, file, backward, normalize)

    type(histogram),intent(in) :: H
    character*(*),  intent(in) :: file
    logical,        intent(in) :: backward, normalize

    integer            :: iBin, iBinMax
    integer, parameter :: iFile = 65
    real :: total(1:3)

    if (.not.allocated(H%yVal)) return

    open(iFile, file=file, status='unknown')
    rewind(iFile)

    iBinMax = ubound(H%yVal,dim=1)

    call WriteHeader(H,iFile)

    if (normalize) then
      total = sum(H%yVal(1:iBinMax,1:3),dim=1)
    else
      total = 1.
    end if

    if (.not. backward) then
      do iBin=1,iBinMax
        write(iFile,fmt2000) H%xMin+(real(iBin)-0.5)*H%xBin, &
                             sum(H%yVal(1:iBin,1:3),dim=1)/total
      end do
    else
      do iBin=1,iBinMax
        write(iFile,fmt2000) H%xMin+(real(iBin)-0.5)*H%xBin, &
                             sum(H%yVal(iBin:iBinMax,1:3),dim=1)/total
      end do
    end if

    close(iFile)

  end subroutine WriteHist_Integrated


  !****************************************************************************
  !****s* histf90/WriteHist_Spline
  ! NAME
  ! subroutine WriteHist_Spline(H,file,add,mul,nPoints)
  ! PURPOSE
  ! Write out the histogram.
  !
  ! The entries are multiplied by 'mul' and 'add' is added.
  !
  ! The binning is subdivided into "nPoints" points per bin, where
  ! every point is calculated via cubic spline interpolation.
  ! (Default: nPoints=5)
  ! INPUTS
  ! * type(histogram) :: H     -- Histogram to be used
  ! * character*(*)   :: file  -- File name to write to
  ! * real            :: add   -- factor to add      [OPTIONAL]
  ! * real            :: mul   -- factor to multiply [OPTIONAL]
  ! * integer         :: nPoints -- number of sub-points [OPTIONAL]
  ! OUTPUT
  ! write to file number
  !
  ! cf. WriteHist
  ! NOTES
  ! The Histogram Data is not affected!!!
  !****************************************************************************
  subroutine WriteHist_Spline(H,file,add,mul,nPoints)
    use cl_splines

    type(histogram),intent(in)          :: H
    integer, parameter :: iFile = 66
    character*(*),  intent(in)          :: file
    real,           intent(in),optional :: add
    real,           intent(in),optional :: mul
    integer,        intent(in),optional :: nPoints

!    real,dimension(1:3) :: S
    integer             :: iBin, iBinMax

    real :: addFak
    real :: mulFak
    integer :: nPP

    real, allocatable :: XX(:),YY1(:),YY2(:),YY3(:)
    real :: x!, h1,h2,h3
    integer :: error
    logical :: success

    type(tspline),dimension(1:3) :: spline

    addFak = 0.
    mulFak = 1.
    nPP = 5

    if (present(add)) addFak = add
    if (present(mul)) mulFak = mul
    if (present(nPoints)) nPP = nPoints

    if (.not.allocated(H%yVal)) return

    open(iFile, file=file, status='unknown')
    rewind(iFile)

    call WriteHeader(H,iFile,mulFak)

    write(iFile,fmt1020) "SPLINE INTERPOLATION",nPP

    iBinMax = ubound(H%yVal,dim=1)

    allocate(XX(iBinMax))
    allocate(YY1(iBinMax))
    allocate(YY2(iBinMax))
    allocate(YY3(iBinMax))

    do iBin=1,iBinMax
       xx(iBin)  = H%xMin+(real(iBin)-0.5)*H%xBin
       YY1(iBin) = H%yVal(iBin,1)/H%xBin*mulFak+addFak
       YY2(iBin) = H%yVal(iBin,2)
       YY3(iBin) = H%yVal(iBin,3)/H%xBin*mulFak+addFak
    end do

    spline(1)=cl_initSpline(xx,YY1)
    spline(2)=cl_initSpline(xx,YY2)
    spline(3)=cl_initSpline(xx,YY3)

    do iBin=1,iBinMax*nPP
       x = H%xMin+(real(iBin)-0.5)*H%xBin/nPP
       write(iFile,fmt2000) x, cl_spline(spline(1),x,success,error),&
            & cl_spline(spline(2),x,success,error),cl_spline(spline(3),x,success,error)
    end do

    deallocate(XX, YY1,YY2,YY3)

    close(iFile)

  end subroutine WriteHist_Spline


  !****************************************************************************
  !****s* histf90/WriteHist_BSpline
  ! NAME
  ! subroutine WriteHist_BSpline(H,file,add,mul,nPoints)
  ! PURPOSE
  ! Write out the histogram.
  !
  ! The entries are multiplied by 'mul' and 'add' is added.
  !
  ! The binning is subdivided into "nPoints" points per bin, where
  ! every point is calculated via cubic B-spline interpolation.
  ! (Default: nPoints=5)
  ! INPUTS
  ! * type(histogram) :: H     -- Histogram to be used
  ! * integer         :: file  -- File name to write to
  ! * real            :: add   -- factor to add      [OPTIONAL]
  ! * real            :: mul   -- factor to multiply [OPTIONAL]
  ! * integer         :: nPoints -- number of sub-points [OPTIONAL]
  ! OUTPUT
  ! write to file number
  !
  ! cf. WriteHist
  ! NOTES
  ! The Histogram Data is not affected!!!
  !****************************************************************************
  subroutine WriteHist_BSpline(H,file,add,mul,nPoints)
    use spline, only: Bsplint2

    type(histogram),intent(in)          :: H
    integer, parameter :: iFile = 66
    character*(*),  intent(in)          :: file
    real,           intent(in),optional :: add
    real,           intent(in),optional :: mul
    integer,        intent(in),optional :: nPoints

    integer             :: iBin, iBinMax

    real :: addFak
    real :: mulFak
    integer :: nPP

    real, allocatable :: XX(:),YY1(:),YY2(:),YY3(:)
    real :: x, h1,h2,h3

    integer :: Ind(-1:2)
    real    :: Pref(-1:2)

    addFak = 0.
    mulFak = 1.
    nPP = 5

    if (present(add)) addFak = add
    if (present(mul)) mulFak = mul
    if (present(nPoints)) nPP = nPoints

    if (.not.allocated(H%yVal)) return

    open(iFile, file=file, status='unknown')
    rewind(iFile)

    call WriteHeader(H,iFile,mulFak)

    write(iFile,fmt1020) "B-SPLINE INTERPOLATION",nPP

    iBinMax = ubound(H%yVal,dim=1)

    allocate(XX(iBinMax))
    allocate(YY1(iBinMax))
    allocate(YY2(iBinMax))
    allocate(YY3(iBinMax))

    do iBin=1,iBinMax
       xx(iBin)  = H%xMin+(real(iBin)-0.5)*H%xBin
       YY1(iBin) = H%yVal(iBin,1)/H%xBin*mulFak+addFak
       YY2(iBin) = H%yVal(iBin,2)
       YY3(iBin) = H%yVal(iBin,3)/H%xBin*mulFak+addFak
    end do

    do iBin=1,iBinMax*nPP
       x = H%xMin+(real(iBin)-0.5)*H%xBin/nPP

       h1 = Bsplint2(xx,YY1, x, Ind,Pref, .false.)
       h2 = Bsplint2(xx,YY2, x, Ind,Pref, .true.)
       h3 = Bsplint2(xx,YY3, x, Ind,Pref, .true.)

       write(iFile,fmt2000) x, h1,h2,h3
    end do

    deallocate(XX, YY1,YY2,YY3)

    close(iFile)

  end subroutine WriteHist_BSpline


  !****************************************************************************
  !****s* histf90/WriteHist_Gauss
  ! NAME
  ! subroutine WriteHist_Gauss(H,file,width_in)
  ! PURPOSE
  ! Write out the histogram smeared with a gaussian distribution
  !
  ! INPUTS
  ! * type(histogram) :: H     -- Histogram to be used
  ! * character*(*)   :: file  -- name of file to write to
  ! * real, OPTIONAL  :: width_in -- width of gaussian
  ! OUTPUT
  ! write to file number
  !
  ! NOTES
  ! The Histogram Data is not affected!!!
  !****************************************************************************
  subroutine WriteHist_Gauss(H,file,width_in)
    use distributions, only: gauss

    type(histogram),intent(in) :: H
    integer, parameter :: iFile = 66
    character*(*),  intent(in) :: file
    real, optional :: width_in

    integer, parameter :: N = 1
    real,dimension(1:3) :: S
    integer             :: i,j,iBinMax
    real :: xi,xj,width

    ! set width of gaussian
    if (present(width_in)) then
      width = width_in
    else
      width = H%xBin ! default: bin width
    end if

    if (.not.allocated(H%yVal)) return

    open(iFile, file=file, status='unknown')
    rewind(iFile)

    call WriteHeader(H,iFile)

    write(iFile,"('### Gauss folding, sigma=',ES11.4,/,'###')") width

    iBinMax = ubound(H%yVal,dim=1)

    do i=1,iBinMax*N
      xi = H%xMin + (i+N/2.-N/2-1)*H%xBin/N
      S = 0.
      do j=max(1,nint((i-10*N*width/H%xBin)/N)),min(iBinMax,nint((i+10*N*width/H%xBin)/N))
        xj = H%xMin + (real(j)-0.5)*H%xBin
        S = S + H%yVal(j,:) * gauss(xi,xj,width)
      end do
      write(iFile,fmt2000) xi, S(1:3)
    end do

    close(iFile)

  end subroutine WriteHist_Gauss



  !****************************************************************************
  !****s* histf90/WriteHist_Novo
  ! NAME
  ! subroutine WriteHist_Novo(H,file,w,d)
  ! PURPOSE
  ! Write out the histogram smeared with a Novosibirsk distribution
  !
  ! INPUTS
  ! * type(histogram) :: H     -- Histogram to be used
  ! * character*(*)   :: file  -- name of file to write to
  ! * real            :: w,d   -- width and skewness of Novosibirsk function
  ! OUTPUT
  ! write to file number
  !
  ! NOTES
  ! The Histogram Data is not affected!!!
  !****************************************************************************
  subroutine WriteHist_Novo(H,file,w,d)
    use distributions, only: novo

    type(histogram),intent(in) :: H
    integer, parameter :: iFile = 67
    character*(*),  intent(in) :: file
    real, intent(in) :: w,d

    integer, parameter :: N = 1
    real,dimension(1:3) :: S
    integer             :: i,j,iBinMax
    real :: xi,xj

    if (.not.allocated(H%yVal)) return

    open(iFile, file=file, status='unknown')
    rewind(iFile)

    call WriteHeader(H,iFile)

    write(iFile,"('### Novosibirsk folding, sigma=',ES11.4,', tau=',ES11.4,/,'###')") w,d

    iBinMax = ubound(H%yVal,dim=1)

    do i=1,iBinMax*N
      xi = H%xMin + (i+N/2.-N/2-1)*H%xBin/N
      S = 0.
      do j=max(1,nint((i-10*N*w/H%xBin)/N)),min(iBinMax,nint((i+10*N*w/H%xBin)/N))
        xj = H%xMin + (real(j)-0.5)*H%xBin
        S = S + H%yVal(j,:) * novo(xi,xj,w,d)
      end do
      write(iFile,fmt2000) xi, S(1:3)
    end do

    close(iFile)

  end subroutine WriteHist_Novo



  !****************************************************************************
  !****f* histf90/ReadHist
  ! NAME
  ! function ReadHist (H, file) result (success)
  ! PURPOSE
  ! Read in a histogram from a file.
  ! INPUTS
  ! * type(histogram) :: H     -- Histogram to be used
  ! * character*(*)   :: file  -- name of file to read
  ! OUTPUT
  ! * H is changed
  ! * return value indicates success or failure
  ! NOTES
  ! Only data points are read (no under-/overflow etc).
  !****************************************************************************
  function ReadHist (H, file) result (success)

    type(histogram) :: H
    character*(*),  intent(in) :: file
    logical :: success

    integer, parameter :: iFile = 68
    character(len=100) :: line
    integer :: i,ios, head, numbins
    real :: x

    success = .false.

    if (H%initialized) call RemoveHist(H)

    open(iFile, file=file, status='old', iostat=ios)
    if (ios /= 0) return
    head = 0
    numbins = 0

    ! determine number of bins and header lines
    do
      read (iFile,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (line(2:2) == '#') then
        head = head + 1
      else
        numbins = numbins + 1
      end if
    end do

    allocate(H%yVal(-1:numbins,3))
    H%yVal = 0.

    rewind(iFile)
    ! discard header, get histogram name
    do i=1,head
      read (iFile,'(A)',iostat=ios) line
      if (ios /= 0) then
        close(iFile)
        call RemoveHist(H)
        return
      end if
      if (i==2) H%name = line(19:)
    end do

    ! read values
    do i=1,numbins
      read (iFile,fmt2000,iostat=ios) x,H%yval(i,1:3)
      if (ios /= 0) then
        close(iFile)
        call RemoveHist(H)
        return
      end if
      if (i==1) then
        H%xmin = x
      else if (i==numbins) then
        H%xmax = x
      end if
    end do

    ! set xbin, xmin, xmax
    H%xbin = (H%xmax - H%xmin)/(numbins-1)
    H%xmin = H%xmin - H%xbin/2.
    H%xmax = H%xmax - H%xbin/2.
    ! correct for bin width
    H%yval(:,1) = H%yval(:,1) * H%xbin

    close(iFile)
    H%initialized = .true.

    success = .true.

  end function ReadHist


  !****************************************************************************
  !****s* histf90/sumHist
  ! NAME
  ! subroutine sumHist (A, B, w)
  ! PURPOSE
  ! Perform a summation of two histograms: A = A + B.
  ! If a weight 'w' is given, do a weighted summation: A = A + w*B
  ! INPUTS
  ! * type(histogramMC) :: A, B         -- Histograms to be used
  ! * real, optional, intent(in) :: w   -- optional weighting factor
  ! OUTPUT
  ! A is changed.
  !****************************************************************************
  subroutine sumHist (A, B, w)
    type(histogram) :: A
    type(histogram), intent(in) :: B
    real, optional, intent(in) :: w

    if (present(w)) then
      A%yval(:,:) = A%yval(:,:) + w * B%yval(:,:)
    else
      A%yval(:,:) = A%yval(:,:) + B%yval(:,:)
    end if

  end subroutine sumHist


  !****************************************************************************
  !****s* histf90/DumpHist
  ! NAME
  ! subroutine DumpHist(H,file,iFile, add,mul)
  ! PURPOSE
  ! Write all the histogram information unformatted (i.e. binary) to a file
  !
  ! INPUTS
  ! * type(histogram) :: H     -- Histogram to be used
  ! * character*(*)   :: file  -- name of file to open and close
  ! * integer,OPTIONAL:: iFile -- File number output to redirect [OPTIONAL]
  ! * real            :: add   -- factor to add      [OPTIONAL]
  ! * real            :: mul   -- factor to multiply [OPTIONAL]
  ! OUTPUT
  ! H is written UNFORMATTED to the given file
  !
  !****************************************************************************
  subroutine DumpHist(H,file,iFile, add,mul)

    type(histogram),intent(in)          :: H
    character*(*),  intent(in)          :: file
    integer,        intent(in),optional :: iFile
    real,           intent(in),optional :: add,mul

    integer :: iF
    real :: addFak,mulFak
    logical :: WriteFaks

    iF=121
    if (present(iFile)) iF = iFile

    open(iF,file=file,status='UNKNOWN',form='UNFORMATTED')
    rewind(iF)

    write(iF) H%xMin,H%xMax,H%xBin,H%xExtreme
    write(iF) H%Name
    write(iF) H%yVal

    WriteFaks = .false.
    addFak = 0.0
    mulFak = 1.0
    if (present(mul)) then
       mulFak = mul
       WriteFaks = .true.
    end if
    if (present(add)) then
       addFak = add
       WriteFaks = .true.
    end if

    if (WriteFaks) write(iF) addFak,mulFak


    close(iF)

  end subroutine DumpHist

  !****************************************************************************
  !****s* histf90/FetchHist
  ! NAME
  ! subroutine FetchHist(H,file,iFile, add,mul,flagOK)
  ! PURPOSE
  ! Read in all the histogram information previously dumped unformatted
  ! (i.e. binary) to a file
  !
  ! INPUTS
  ! * character*(*)   :: file  -- name of file to open and close
  ! * integer,OPTIONAL:: iFile -- File number input to redirect [OPTIONAL]
  ! OUTPUT
  ! * type(histogram) :: H     -- Histogram to be used
  ! * real            :: add   -- factor to add      [OPTIONAL]
  ! * real            :: mul   -- factor to multiply [OPTIONAL]
  ! * logical         :: flagOK -- flag, if reading was okay [OPTIONAL]
  !
  ! H is read UNFORMATTED from the given file. Sizes are calculated as in
  ! CreateHist, also memory is allocated.
  !
  ! NOTES
  ! No checks about input are performed!
  !****************************************************************************
  subroutine FetchHist(H,file,iFile, add,mul, flagOK)

    type(histogram),intent(inout)       :: H
    character*(*),  intent(in)          :: file
    integer,        intent(in),optional :: iFile
    real,           intent(out),optional:: add,mul
    logical,        intent(out),optional:: flagOK

    integer :: iF
    integer :: L
    real :: addFak,mulFak
    integer :: ios

    iF=121
    if (present(iFile)) iF = iFile

    open(iF,file=file,status='UNKNOWN',form='UNFORMATTED',iostat=ios)
    if (ios.ne.0) then
       close(iF)
       if (present(flagOK)) flagOK=.false.
       return
    end if
    rewind(iF)

    read(iF,iostat=ios) H%xMin,H%xMax,H%xBin,H%xExtreme
    if (ios.ne.0) then
       close(iF)
       if (present(flagOK)) flagOK=.false.
       return
    end if
    read(iF) H%Name

    if (allocated(H%yVal)) deallocate(H%yVal)

    L = nint( (H%xMax-H%xMin)/H%xBin )+1
    allocate(H%yVal(-1:L,3))

    read(iF) H%yVal

    if (present(add).or.present(mul)) then
       addFak = 0.0
       mulFak = 1.0
       read(iF,iostat=ios) addFak,mulFak
       if (ios.ne.0) then
          write(*,*) 'FetchHist: old file version, no add/mul info!'
          addFak = 0.0
          mulFak = 1.0
       end if
       if (present(add)) add=addFak
       if (present(mul)) mul=mulFak
    end if

    close(iF)
    if (present(flagOK)) flagOK=.true.

  end subroutine FetchHist

  !****************************************************************************

end module histf90
