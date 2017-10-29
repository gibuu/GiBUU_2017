!******************************************************************************
!****m* /hist2Df90
! NAME
! module hist2Df90
! PURPOSE
! Encapsulate all routines and datas for 2D Histograms.
!
! Features of Histograms provided by this module:
! - store paramaters of the (x1,x2)-binning
! - enable two y-values (y and y2)
! - track keeping of under-/over-score the given extreme values of x.
! - provide simple-to-understand output routines
! - provide simple histogram arithmetic
! - many others...
!
! NOTES
! Programming is fairly similar to "module histf90"
!
! INPUTS
! ...(This module needs no input)
!******************************************************************************
module hist2Df90

  use histf90

  implicit none
  private

  !****************************************************************************
  !****t* hist2Df90/histogram2D
  ! NAME
  ! type histogram2D
  ! PURPOSE
  ! Type definition to store all information for a 2D Histogram.
  ! SOURCE
  !
  type histogram2D
     real             :: xMin(2)        ! smallest x-value: (x1,x2)
     real             :: xMax(2)        ! largest x-value:  (x1,x2)
     real             :: xBin(2)        ! width of x-Bin:   (x1,x2)
     real             :: xExtreme(2,2)  ! extremes of used x-values
                                        !  (x1,x2), (min,max)
     character*(NameLength) :: Name     ! name to be written
     real,allocatable :: yVal(:,:,:)    ! histogram values:
                                        !  (x1),(x2),(yy1,yy2,yy3)
  end type histogram2D
  !****************************************************************************

  public :: histogram2D
  public :: CreateHist2D
  public :: RemoveHist2D
  public :: CopyHist2D
  public :: AddHist2D
  public :: WriteHist2D_SYSTEM
  public :: WriteHist2D_Gnuplot
  public :: WriteHist2D_Gnuplot_Bspline
!   public :: DivideHist2D
  public :: ReadHist2D_Gnuplot
  public :: FetchHist2D
  public :: IntegrateHist2D
  public :: AverageHist2D

  !****************************************************************************
  !****s* hist2Df90/AddHist2D
  ! NAME
  ! subroutine AddHist2D(H, x,y,y2)
  !
  ! subroutine AddHist2D(H1,H2, x,y,y2)
  ! PURPOSE
  ! Add to the given histogram(s) at the given (x1,x2)-value the weight y (y2).
  ! If y2 is given, also the second column is filled.
  ! INPUTS
  ! * type(histogram2D) :: H  -- Histogram to be used
  ! or:
  ! * type(histogram2D) :: H1,H2  -- Histograms to be used
  ! * real,dimension(2) :: x  -- (x1,x2)-value
  ! * real              :: y  -- weight to be added
  ! * real              :: y2 -- second weight to be added [OPTIONAL]
  ! OUTPUT
  ! H is changed
  ! NOTES
  ! This routine is overloaded: Giving two histograms at the same time is
  ! like calling the routine twice.
  !****************************************************************************
  interface AddHist2D
     module procedure AddHist2D1,AddHist2D2
  end interface

contains

  !****************************************************************************
  !****s* hist2Df90/CreateHist2D
  ! NAME
  ! subroutine CreateHist2D(H, HName,x1,x2,bin,verbose)
  ! PURPOSE
  ! This is the Constructor of a 2D-Histogram!
  ! Allocate Memory for the entries and put additional variables to
  ! their default.

  ! INPUTS
  ! * type(histogram2D) :: H         -- Histogram to be created
  ! * character*(*)     :: HName     -- Name of Histogram
  ! * real,dimension(2) :: x1,x2,bin -- Minimal/maximal values for
  !                                   x-coordinates to be considered,
  !                                   bin-width
  ! * logical            :: verbose  -- switch for verbosity [OPTIONAL]
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine CreateHist2D(H, HName,x1,x2,Bin,verbose)

    type(histogram2D),intent(inout) :: H
    character*(*),intent(in) :: HName
    real,dimension(2),intent(in) :: x1
    real,dimension(2),intent(in) :: x2
    real,dimension(2),intent(in) :: bin
    logical, intent(in),optional :: verbose

    integer,dimension(2) :: L

    H%xMin = x1
    H%xMax = x2
    H%xBin = bin
    H%xExtreme(1:2,1) =  99e9
    H%xExtreme(1:2,2) = -99e9


    if (len(HName).gt.NameLength) then
       H%Name = HName(1:NameLength)
    else
       H%Name = HName
    end if

    if (allocated(H%yVal)) deallocate(H%yVal)

    L = nint( (x2-x1)/bin )+1
    allocate(H%yVal(-1:L(1),-1:L(2),3))

    H%yVal = 0.

    if (present(verbose)) then
       if (verbose) then
          write(*,*) '***** CreateHist2D:',hName,L
          write(*,*) '***** --x1        :',x1
          write(*,*) '***** --x2        :',x2
          write(*,*) '***** --bin       :',bin
       end if
    end if

  end subroutine CreateHist2D


  !****************************************************************************
  !****s* hist2Df90/WriteHist2D_SYSTEM
  ! NAME
  ! subroutine WriteHist2D_SYSTEM(H, iFile)
  ! PURPOSE
  ! Write out System information of the histogram
  ! INPUTS
  ! * type(histogram2D) :: H     -- Histogram to be used
  ! * integer      :: iFile -- File number output to redirect
  ! OUTPUT
  ! written to file iFile
  !****************************************************************************
  subroutine WriteHist2D_SYSTEM(H, iFile)

    type(histogram2D),intent(in)          :: H
    integer,          intent(in)          :: iFile

    write(iFile,'(A)') H%Name
    write(iFile,*) H%xMin
    write(iFile,*) H%xMax
    write(iFile,*) H%xBin

  end subroutine WriteHist2D_SYSTEM


  !****************************************************************************
  !****s* hist2Df90/RemoveHist2D
  ! NAME
  ! subroutine RemoveHist2D(H)
  ! PURPOSE
  ! Free the allocated memory
  ! INPUTS
  ! * type(histogram2D) :: H  -- Histogram to be used
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine RemoveHist2D(H)
    type(histogram2D),intent(inout) :: H

    if (allocated(H%yVal)) deallocate(H%yVal)

  end subroutine RemoveHist2D


  !****************************************************************************
  !****s* hist2Df90/CopyHist2D
  ! NAME
  ! subroutine CopyHist2D(H1,H2)
  ! PURPOSE
  ! Copies Histogram H1 into H2 (H2 is overwritten)
  ! INPUTS
  ! * type(histogram2D) :: H1,H2  -- Histograms to be used
  ! OUTPUT
  ! H2 is changed
  !****************************************************************************
  subroutine CopyHist2D(H1,H2)

    type(histogram2D),intent(inout) :: H1,H2

    integer,dimension(2) :: L

    call RemoveHist2D(H2)

    H2%xMin = H1%xMin
    H2%xMax = H1%xMax
    H2%xBin = H1%xBin
    H2%Name = H1%Name
    H2%xExtreme = H1%xExtreme

    L = nint( (H2%xMax-H2%xMin)/H2%xBin )+1
    allocate(H2%yVal(-1:L(1),-1:L(2),3))

    H2%yVal = H1%yVal

  end subroutine CopyHist2D


  !****************************************************************************
  ! cf. interface AddHist2D
  !****************************************************************************
  subroutine AddHist2D1(H, x,y,y2)

    type(histogram2D), intent(inout)       :: H
    real,dimension(2), intent(in)          :: x
    real,              intent(in)          :: y
    real,              intent(in),optional :: y2

    integer,dimension(2) :: iBin
    real :: yy
    integer :: i

    yy = 0.
    if (present(y2)) yy = y2

!   write(*,*) 'AddHist: ',x,y,yy


    do i=1,2

       !...extremes (counted):

       if (x(i)<H%xExtreme(i,1)) H%xExtreme(i,1)=x(i)
       if (x(i)>H%xExtreme(i,2)) H%xExtreme(i,2)=x(i)

       !...extremes (measured):

       if (x(i) >= H%xMin(i) .and. x(i) <= H%xMax(i)) then
          iBin(i) = int( (x(i)-H%xMin(i))/H%xBin(i) )+1
       else if (x(i) < H%xMin(i)) then
          iBin(i) = -1
       else
          iBin(i) =  0
       end if

    end do

    H%yVal(iBin(1),iBin(2),1) = H%yVal(iBin(1),iBin(2),1)+y
    H%yVal(iBin(1),iBin(2),2) = H%yVal(iBin(1),iBin(2),2)+1.
    H%yVal(iBin(1),iBin(2),3) = H%yVal(iBin(1),iBin(2),3)+yy

  end subroutine AddHist2D1


  subroutine AddHist2D2(H1,H2, x,y,y2)

    type(histogram2D), intent(inout)       :: H1,H2
    real,dimension(2), intent(in)          :: x
    real,              intent(in)          :: y
    real,              intent(in),optional :: y2

   if (present(y2)) then
       call AddHist2D1(H1, x,y,y2)
       call AddHist2D1(H2, x,y,y2)
    else
       call AddHist2D1(H1, x,y)
       call AddHist2D1(H2, x,y)
    end if
  end subroutine AddHist2D2


  !****************************************************************************
  !****s* hist2Df90/WriteHist2D_Gnuplot
  ! NAME
  ! subroutine WriteHist2D_Gnuplot(H,iFile,add,mul, iColumn,DoAve,MaxVal,H2,file,dump,SwapXY)
  ! PURPOSE
  ! Write out the 2D-histogram. Format is suitable as input to
  ! gnuplots "splot" command
  !
  ! The entries are multiplied by 'mul' and 'add' is added.
  !
  ! If iColumn is not given, all 3 columns are writen.
  !
  ! INPUTS
  ! * type(histogram2D) :: H     -- Histogram to be used
  ! * integer      :: iFile -- File number output to redirect [OPTIONAL]
  ! * real         :: add   -- factor to add      [OPTIONAL]
  ! * real         :: mul   -- factor to multiply [OPTIONAL]
  ! * integer      :: iColumn -- Column to write  [OPTIONAL]
  ! * logical      :: DoAve -- write also "average" (cf Notes) [OPTIONAL]
  ! * real         :: maxVal -- value used instead "divide by zero" [OPTIONAL]
  ! * type(histogram) :: H2    -- Histogram to divide by [OPTIONAL]
  ! * character*(*)   :: file  -- name of file to open and close [OPTIONAL]
  ! * logical         :: dump  -- if true, also dump it binary to file [OPTIONAL]
  ! * logical         :: SwapXY -- write data sorted according y:x:z instead of x:y:z [OPTIONAL]
  ! OUTPUT
  ! write to file number
  ! NOTES
  ! The Histogram Data is not affected!!!
  !
  ! if DoAve is given:
  ! Write out the 2D-histogram, but now also divide column 3 by column 1.
  ! If column 2 is zero, i.e. no entries, then parameter maxVal is
  ! written instead.
  ! This column output is NOT divided by bin widths.
  !
  !****************************************************************************
  subroutine WriteHist2D_Gnuplot(H,iFile_in,add,mul, iColumn,DoAve,MaxVal,H2,file,dump,SwapXY)

    type(histogram2D),intent(in)          :: H
    integer,          intent(in),optional :: iFile_in
    real,             intent(in),optional :: add
    real,             intent(in),optional :: mul
    integer,          intent(in),optional :: iColumn
    logical,          intent(in),optional :: DoAve
    real,             intent(in),optional :: maxVal
    type(histogram),  intent(in),optional :: H2
    character*(*),    intent(in),optional :: file
    logical,          intent(in),optional :: dump
    logical,          intent(in),optional :: SwapXY

    integer, dimension(2) :: iBinMax
    integer :: iBin1,iBin2, iFile

    real :: addFak
    real :: mulFak
    real :: Z, maxZ
    logical :: writeZ, SwapXY_
    real, dimension(2) :: xMin, xBin

    real,allocatable :: yVal(:,:,:)

    addFak = 0.
    mulFak = 1.
    maxZ = 99.0
    writeZ = .false.
    iFile = 65
     SwapXY_ = .false.

    if (present(add)) addFak = add
    if (present(mul)) mulFak = mul
    if (present(maxVal)) maxZ = maxVal
    if (present(DoAve)) writeZ = DoAve
    if (present(iFile_in)) iFile = iFile_in
    if (present(SwapXY)) SwapXY_ = SwapXY

    if (.not.allocated(H%yVal)) return

    if (present(file)) then
       open(iFile, file=file, status='unknown')
       rewind(iFile)
    end if

    mulFak = mulFak / (H%xBin(1)*H%xBin(2))

    iBinMax(1) = ubound(H%yVal,dim=1)
    iBinMax(2) = ubound(H%yVal,dim=2)

    allocate(yVal(-1:iBinMax(1),-1:iBinMax(2),3))
    yVal = H%yVal

    if (SwapXY_) then
       xMin(2:1:-1) = H%xMin
       xBin(2:1:-1) = H%xBin
    else
       xMin = H%xMin
       xBin = H%xBin
    end if

    if (present(H2)) then
       if (ubound(H2%yVal,dim=1).ne.iBinMax(1)) then
          write(*,*) 'WriteHist2D_Gnuplot: ERROR, dim not equal! STOP'
          write(*,*) '  histogram=',H%Name
          write(*,*) ubound(H2%yVal,dim=1),iBinMax(1)
          stop
       end if

       do iBin1=-1,iBinMax(1)
          if (H2%yVal(iBin1,1) .eq. 0.0) then
             do iBin2=-1,iBinMax(2)
                if (yVal(iBin1,iBin2,1) .ne. 0.0) yVal(iBin1,iBin2,1) = maxZ
                if (yVal(iBin1,iBin2,3) .ne. 0.0) yVal(iBin1,iBin2,3) = maxZ
             end do
          else
             do iBin2=-1,iBinMax(2)
                yVal(iBin1,iBin2,1) = yVal(iBin1,iBin2,1)/H2%yVal(iBin1,1)
                yVal(iBin1,iBin2,3) = yVal(iBin1,iBin2,3)/H2%yVal(iBin1,1) ! we divide by column 1!
             end do
          end if
       end do
    end if

    if (SwapXY_) then

       do iBin2=1,iBinMax(2)
          do iBin1=1,iBinMax(1)
             if (writeZ) then
                if (yVal(iBin1,iBin2,2) > 0) then
                   if (yVal(iBin1,iBin2,1) .eq. 0.0) then
                      Z = maxZ
                   else
                      Z = yVal(iBin1,iBin2,3)/yVal(iBin1,iBin2,1)
                   end if
                else
                   Z = maxZ
                end if
             end if

             if (present(iColumn)) then
                if (writeZ) then
                   write(iFile,1000) xMin+((/iBin2,iBin1/)*1.0-0.5)*xBin, Z, &
                        yVal(iBin1,iBin2,iColumn)*MulFak+AddFak
                else
                   write(iFile,1000) xMin+((/iBin2,iBin1/)*1.0-0.5)*xBin, &
                        yVal(iBin1,iBin2,iColumn)*MulFak+AddFak
                end if
             else
                if (writeZ) then
                   write(iFile,1000) xMin+((/iBin2,iBin1/)*1.0-0.5)*xBin, Z, &
                        yVal(iBin1,iBin2,1:3)*MulFak+AddFak
                else
                   write(iFile,1000) xMin+((/iBin2,iBin1/)*1.0-0.5)*xBin, &
                        yVal(iBin1,iBin2,1:3)*MulFak+AddFak
                end if
             end if
          end do
          write(iFile,*)
       end do

    else

       do iBin1=1,iBinMax(1)
          do iBin2=1,iBinMax(2)
             if (writeZ) then
                if (yVal(iBin1,iBin2,2) > 0) then
                   if (yVal(iBin1,iBin2,1) .eq. 0.0) then
                      Z = maxZ
                   else
                      Z = yVal(iBin1,iBin2,3)/yVal(iBin1,iBin2,1)
                   end if
                else
                   Z = maxZ
                end if
             end if

             if (present(iColumn)) then
                if (writeZ) then
                   write(iFile,1000) xMin+((/iBin1,iBin2/)*1.0-0.5)*xBin, Z, &
                        yVal(iBin1,iBin2,iColumn)*MulFak+AddFak
                else
                   write(iFile,1000) xMin+((/iBin1,iBin2/)*1.0-0.5)*xBin, &
                        yVal(iBin1,iBin2,iColumn)*MulFak+AddFak
                end if
             else
                if (writeZ) then
                   write(iFile,1000) xMin+((/iBin1,iBin2/)*1.0-0.5)*xBin, Z, &
                        yVal(iBin1,iBin2,1:3)*MulFak+AddFak
                else
                   write(iFile,1000) xMin+((/iBin1,iBin2/)*1.0-0.5)*xBin, &
                        yVal(iBin1,iBin2,1:3)*MulFak+AddFak
                end if
             end if
          end do
          write(iFile,*)
       end do

    end if

    if (present(file)) then
       close(iFile)
    end if

    if (present(dump)) then
       if (dump) then
          if (present(file)) then
             call DumpHist2D(H,file//".bin",iFile,addFak,mulFak)
          else
             call DumpHist2D(H,"DumpHist2D.bin",iFile,addFak,mulFak)
          end if
       end if
    end if

1000 FORMAT(1X,' ',1P,2E12.4,4E12.4)

  end subroutine WriteHist2D_Gnuplot

  !****************************************************************************
  !****s* hist2Df90/WriteHist2D_Gnuplot_Bspline
  ! NAME
  ! subroutine WriteHist2D_Gnuplot_Bspline(H,iFile,add,mul, iColumn,nPoints,DoAve,MaxVal)
  ! PURPOSE
  ! Write out the 2D-histogram. Format is suitable as input to
  ! gnuplots "splot" command
  !
  ! The entries are multiplied by 'mul' and 'add' is added.
  !
  ! If iColumn is not given, all 3 columns are writen.
  !
  ! INPUTS
  ! * type(histogram2D) :: H     -- Histogram to be used
  ! * integer           :: iFile -- File number output to redirect
  ! * real              :: add   -- factor to add      [OPTIONAL]
  ! * real              :: mul   -- factor to multiply [OPTIONAL]
  ! * integer           :: iColumn -- Column to write  [OPTIONAL]
  ! * integer           :: nPoints -- number of sub-points [OPTIONAL]
  ! * logical      :: DoAve -- write also "average" (cf Notes) [OPTIONAL]
  ! * real         :: maxVal -- value used instead "divide by zero" [OPTIONAL]
  ! OUTPUT
  ! write to file number
  ! NOTES
  ! The Histogram Data is not affected!!!
  !
  ! cf. WriteHist2D_Gnuplot
  !****************************************************************************
  subroutine WriteHist2D_Gnuplot_Bspline(H,iFile,add,mul, iColumn,nPoints,DoAve,MaxVal)
    use spline, only: Bsplint2

    type(histogram2D),intent(in)          :: H
    integer,          intent(in)          :: iFile
    real,             intent(in),optional :: add
    real,             intent(in),optional :: mul
    integer,          intent(in),optional :: iColumn
    integer,          intent(in),optional :: nPoints
    logical,          intent(in),optional :: DoAve
    real,             intent(in),optional :: maxVal

    integer, dimension(2) :: iBinMax
    integer :: iBin1,iBin2, iC,i
    real :: Z, maxZ
    logical :: writeZ

    real :: addFak
    real :: mulFak
    integer :: nPP

    real, allocatable :: XX(:),YY(:),ZZx(:),ZZy(:),ratio(:,:)
    integer :: Indx(-1:2),Indy(-1:2)
    real    :: Prefx(-1:2),Prefy(-1:2)
    real :: xVal,yVal,zVal(0:3)
    logical :: flag

    addFak = 0.
    mulFak = 1.
    nPP = 5
    maxZ = 99.0
    writeZ = .false.

    if (present(add)) addFak = add
    if (present(mul)) mulFak = mul
    if (present(nPoints)) nPP = nPoints
    if (present(maxVal)) maxZ = maxVal
    if (present(DoAve)) writeZ = DoAve

    mulFak = mulFak / (H%xBin(1)*H%xBin(2))

    iBinMax(1) = ubound(H%yVal,dim=1)
    iBinMax(2) = ubound(H%yVal,dim=2)

    allocate(XX(iBinMax(1)))
    allocate(ZZx(iBinMax(1)))
    allocate(YY(iBinMax(2)))
    allocate(ZZy(iBinMax(2)))
    allocate(ratio(iBinMax(1),iBinMax(2)))

    do iBin1=1,iBinMax(1)
       XX(iBin1) = H%xMin(1)+(iBin1*1.0-0.5)*H%xBin(1)
    end do
    do iBin2=1,iBinMax(2)
       YY(iBin2) = H%xMin(2)+(iBin2*1.0-0.5)*H%xBin(2)
    end do

    if (writeZ) then
       do iBin1=1,iBinMax(1)
          do iBin2=1,iBinMax(2)
             if (H%yVal(iBin1,iBin2,2) > 0) then
                if (H%yVal(iBin1,iBin2,1) .eq. 0.0) then
                   Z = maxZ
                else
                   Z = H%yVal(iBin1,iBin2,3)/H%yVal(iBin1,iBin2,1)
                end if
             else
                Z = maxZ
             end if
             ratio(iBin1,iBin2) = Z
          end do
       end do
    end if

    do iBin1=1,iBinMax(1)*nPP
       xVal = H%xMin(1)+(iBin1*1.0-0.5)*H%xBin(1)/nPP
       ZZx = 0.
       zVal(0) = Bsplint2(xx,ZZx, xVal, Indx,Prefx, .false.)

       do iBin2=1,iBinMax(2)*nPP
          yVal = H%xMin(2)+(iBin2*1.0-0.5)*H%xBin(2)/nPP

          if (writeZ) then
             flag = .false.
             do i=-1,2
                ZZy = ratio(Indx(i),:)
                ZZx(Indx(i)) = Bsplint2(yy,ZZy, yVal, Indy,Prefy, flag)
                flag = .true.
             end do
             zVal(0) = Bsplint2(xx,ZZx, xVal, Indx,Prefx, .true.)
          end if

          if (present(iColumn)) then
             iC = iColumn

             flag = .false.
             do i=-1,2
                ZZy = H%yVal(Indx(i),1:iBinMax(2),iC)*MulFak+AddFak
                ZZx(Indx(i)) = Bsplint2(yy,ZZy, yVal, Indy,Prefy, flag)
                flag = .true.
             end do
             zVal(iC) = Bsplint2(xx,ZZx, xVal, Indx,Prefx, .true.)
             if (writeZ) then
                write(iFile,1000) xVal,yVal,zVal(0),zVal(iC)
             else
                write(iFile,1000) xVal,yVal,zVal(iC)
             end if
          else
             do iC=1,3
                flag = .false.
                do i=-1,2
                   ZZy = H%yVal(Indx(i),1:iBinMax(2),iC)*MulFak+AddFak
                   ZZx(Indx(i)) = Bsplint2(yy,ZZy, yVal, Indy,Prefy, flag)
                   flag = .true.
                end do
                zVal(iC) = Bsplint2(xx,ZZx, xVal, Indx,Prefx, .true.)
             end do
             if (writeZ) then
                write(iFile,1000) xVal,yVal,zVal(0:3)
             else
                write(iFile,1000) xVal,yVal,zVal(1:3)
             end if
          end if
       end do
       write(iFile,*)
    end do

    deallocate(xx,yy,ZZx,ZZy,ratio)

1000 FORMAT(1X,' ',1P,2E12.4,4E12.4)

  end subroutine WriteHist2D_Gnuplot_Bspline


  !****************************************************************************
  !****s* hist2Df90/DivideHist2D
  ! NAME
  ! subroutine DivideHist2D(H1,H2, maxVal, MulWithBin)
  ! PURPOSE
  ! calculate H1 = H1/H2 on a bin-by-bin basis
  !
  ! INPUTS
  ! * type(histogram2D) :: H1,H2  -- Histograms to be used
  ! * real     :: maxVal      -- value used instead "divide by zero" [OPTIONAL]
  ! * logical  :: MulWithBin  -- flag: multiply with bin-width       [OPTIONAL]
  !
  ! Without explicit "MulWithBin" parameter, the default behaviour is as
  ! with "MulWithBin=.TRUE."
  !
  ! OUTPUT
  ! H1 is changed
  !
  ! NOTES
  ! This Dividing-Routine ignores all EXTREME values etc. !
  !
  ! The "MulWithBin=.TRUE." feature is provided, because during output the
  ! bin entry values are divided by the bin-width (as usual). Therefore if
  ! you divide a histogram by itself and writ it out, you would not get the
  ! bin-Values "1" but "1/bin-width". This flag here corrects for this.
  !****************************************************************************
!   subroutine DivideHist2D(H1,H2, maxVal, MulWithBin)
!
!     type(histogram2D),intent(inout)       :: H1
!     type(histogram2D),intent(in)          :: H2
!     real,             intent(in),optional :: maxVal
!     logical,          intent(in),optional :: MulWithBin
!
!     real :: maxZ
!     real :: BinFak
!
!     integer, dimension(2) :: iBinMax
!     integer :: iBin1,iBin2,i
!
!     maxZ = 99.0
!     BinFak = 1.0
!
!
!     if (present(maxVal)) maxZ = maxVal
!
!     if (present(MulWithBin)) then
!        if (MulWithBin) BinFak = (H1%xBin(1)*H1%xBin(2))
!     else
!        BinFak = (H1%xBin(1)*H1%xBin(2))
!     endif
!
!     iBinMax(1) = uBound(H1%yVal,dim=1)
!     iBinMax(2) = uBound(H1%yVal,dim=2)
!
!     if ((iBinMax(1).ne.uBound(H2%yVal,dim=1)).or.(iBinMax(2).ne.uBound(H2%yVal,dim=2))) then
!        write(*,*) '!!! DivideHist2D: bounds are different!'
!        write(*,*) iBinMax(1),iBinMax(2),uBound(H2%yVal,dim=1),uBound(H2%yVal,dim=1)
!        stop
!     endif
!
! !    H%yVal(iBin1,iBin2,1) = H%yVal(iBin1,iBin2,1)+y
! !    H%yVal(iBin1,iBin2,2) = H%yVal(iBin1,iBin2,2)+1.
! !    H%yVal(iBin1,iBin2,3) = H%yVal(iBin1,iBin2,3)+yy
!
!     do iBin1=1,iBinMax(1)
!        do iBin2=1,iBinMax(2)
!           do i=1,3
!              if (H2%yVal(iBin1,iBin2,i).eq.0.0) then
!                 if (H1%yVal(iBin1,iBin2,i).ne.0.0) H1%yVal(iBin1,iBin2,i) = maxZ
!              else
!                 H1%yVal(iBin1,iBin2,i)=H1%yVal(iBin1,iBin2,i)/H2%yVal(iBin1,iBin2,i)
!              endif
!           end do
!        end do
!     end do
!
!     H1%yVal = H1%yVal * BinFak
!
!   end subroutine DivideHist2D


  !****************************************************************************
  !****s* hist2Df90/ReadHist2D_Gnuplot
  ! NAME
  ! subroutine ReadHist2D_Gnuplot(H,iFile,add,mul, iColumn, DoAve)
  ! PURPOSE
  ! Read in a 2D-histogram from a file, which has been written in a gnuplot
  ! format.
  !
  !****************************************************************************
  subroutine ReadHist2D_Gnuplot(H,iFile,add,mul, iColumn, DoAve)

    type(histogram2D),intent(inout)       :: H
    integer,          intent(in)          :: iFile
    real,             intent(in),optional :: add
    real,             intent(in),optional :: mul
    integer,          intent(in),optional :: iColumn
    logical,          intent(in),optional :: DoAve

    integer, dimension(2) :: iBinMax
    integer :: iBin1,iBin2

    real :: addFak
    real :: mulFak

    real :: xDum, yDum,zDum(0:3)
    character*100 :: Buf
    logical :: writeZ

    addFak = 0.
    mulFak = 1.
    writeZ = .false.

    if (present(add)) addFak = add
    if (present(mul)) mulFak = mul
    if (present(DoAve)) writeZ = DoAve

    mulFak = mulFak / (H%xBin(1)*H%xBin(2))

    iBinMax(1) = ubound(H%yVal,dim=1)
    iBinMax(2) = ubound(H%yVal,dim=2)

    do iBin1=1,iBinMax(1)
       do iBin2=1,iBinMax(2)
          if (present(iColumn)) then
             if (writeZ) then
                read(iFile,*) xDum,yDum, zDum(0),zDum(iColumn)
             else
                read(iFile,*) xDum,yDum, zDum(iColumn)
             end if
             H%yVal(iBin1,iBin2,iColumn) = (zDum(iColumn)-AddFak)/MulFak
          else
             if (writeZ) then
                read(iFile,*) xDum,yDum, zDum(0:3)
             else
                read(iFile,*) xDum,yDum, zDum(1:3)
             end if
             H%yVal(iBin1,iBin2,1:3) = (zDum(1:3)-AddFak)/MulFak
          end if
       end do
       read(iFile,'(A)') Buf
    end do


  end subroutine ReadHist2D_Gnuplot


  !****************************************************************************
  !****s* hist2Df90/DumpHist2D
  ! NAME
  ! subroutine DumpHist2D(H,file,iFile)
  ! PURPOSE
  ! Write all the histogram information unformatted (i.e. binary) to a file
  !
  ! INPUTS
  ! * type(histogram2D):: H     -- Histogram to be used
  ! * character*(*)    :: file  -- name of file to open and close
  ! * integer,OPTIONAL :: iFile -- File number output to redirect [OPTIONAL]
  ! * real            :: add   -- factor to add      [OPTIONAL]
  ! * real            :: mul   -- factor to multiply [OPTIONAL]
  ! OUTPUT
  ! H is written UNFORMATTED to the given file
  !
  !****************************************************************************
  subroutine DumpHist2D(H,file,iFile, add,mul)

    type(histogram2D),intent(in)          :: H
    character*(*),    intent(in)          :: file
    integer,          intent(in),optional :: iFile
    real,             intent(in),optional :: add,mul

    real :: addFak,mulFak
    logical :: WriteFaks

    integer :: iF

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

  end subroutine DumpHist2D


  !****************************************************************************
  !****s* hist2Df90/FetchHist2D
  ! NAME
  ! subroutine FetchHist2D(H,file,iFile, add,mul, flagOK)
  ! PURPOSE
  ! Read in all the histogram information previously dumped unformatted
  ! (i.e. binary) to a file
  !
  ! INPUTS
  ! * character*(*)   :: file  -- name of file to open and close
  ! * integer,OPTIONAL:: iFile -- File number input to redirect [OPTIONAL]
  ! OUTPUT
  ! * type(histogram2D) :: H     -- Histogram to be used
  ! * logical           :: flagOK -- flag, if reading was okay [OPTIONAL]
  !
  ! H is read UNFORMATTED from the given file. Sizes are calculated as in
  ! CreateHist, also memory is allocated.
  !
  ! NOTES
  ! No checks about input are performed!
  !****************************************************************************
  subroutine FetchHist2D(H,file,iFile, add,mul, flagOK)

    type(histogram2D),intent(inout)       :: H
    character*(*),    intent(in)          :: file
    integer,          intent(in),optional :: iFile
    real,             intent(out),optional:: add,mul
    logical,          intent(out),optional:: flagOK

    integer :: iF
    integer,dimension(2) :: L
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

    if (allocated(H%yVal)) deallocate(H%yVal)

    read(iF,iostat=ios) H%xMin,H%xMax,H%xBin,H%xExtreme
    if (ios.ne.0) then
       close(iF)
       if (present(flagOK)) flagOK=.false.
       return
    end if
    read(iF) H%Name

    L = nint( (H%xMax-H%xMin)/H%xBin)+1
    allocate(H%yVal(-1:L(1),-1:L(2),3))

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

  end subroutine FetchHist2D


  !****************************************************************************
  !****s* hist2Df90/IntegrateHist2D
  ! NAME
  ! subroutine IntegrateHist2D(H2D,H,axis)
  ! PURPOSE
  ! Create a 1D histogram, which represents the integration along one axis
  ! of the 2D histogram
  ! INPUTS
  ! * type(histogram2D) :: H2D  -- 2D histogram to integrate
  ! * integer           :: axis -- axis to preserve (integration is done
  !   for the other axis)
  ! OUTPUT
  ! * type(histogram) :: H      -- 1D histogram
  !****************************************************************************
  subroutine IntegrateHist2D(H2D,H,axis)
    use CALLSTACK

    type(histogram2D),intent(in) :: H2D
    type(histogram),intent(inout) :: H
    integer,intent(in) :: axis

    if ((axis.lt.1).or.(axis.gt.2)) then
       write(*,*) 'axis=',axis
       call TRACEBACK()
    end if

    call CreateHist(H,&
         & "INT("//achar(51-axis)//"):"//H2D%Name, &
         & H2D%xmin(axis),H2D%xmax(axis),H2D%xBin(axis))

    H%xExtreme(:) = H2D%xExtreme(axis,:)

    H%yVal = SUM(H2D%yVal,dim=(3-axis))*H2D%xBin(3-axis)

  end subroutine IntegrateHist2D


  !****************************************************************************
  !****s* hist2Df90/AverageHist2D
  ! NAME
  ! subroutine AverageHist2D(H2D,H,axis)
  ! PURPOSE
  ! Create a 1D histogram, which represents the average along one axis
  ! of the 2D histogram
  ! INPUTS
  ! * type(histogram2D) :: H2D  -- 2D histogram to integrate
  ! * integer           :: axis -- axis to preserve (integration is done
  !   for the other axis)
  ! OUTPUT
  ! * type(histogram) :: H      -- 1D histogram
  !****************************************************************************
  subroutine AverageHist2D(H2D,H,axis)
    use CALLSTACK

    type(histogram2D),intent(in) :: H2D
    type(histogram),intent(inout) :: H
    integer,intent(in) :: axis

    integer :: L,i,j
    real,allocatable :: yVal2(:,:)
    real :: xval

    if ((axis.lt.1).or.(axis.gt.2)) then
       write(*,*) 'axis=',axis
       call TRACEBACK()
    end if

    call CreateHist(H,&
         & "AVE("//achar(51-axis)//"):"//H2D%Name, &
         & H2D%xmin(axis),H2D%xmax(axis),H2D%xBin(axis))

    H%xExtreme(:) = H2D%xExtreme(axis,:)

    ! NOTE, we can not use the SUM statement as in IntegrateHist2D,
    ! since this would also count the 'extreme' entries

    L = nint( (H%xMax-H%xMin)/H%xBin )+1
    allocate(yVal2(-1:L,3))
    yVal2 = 0.0

    do i=1,L
       do j=1,ubound(H2D%yval,dim=(3-axis))
          select case (axis)
          case (1)
             xval = H2D%xMin(2)+(j*1.0-0.5)*H2D%xBin(2)
             H%yval(i,:)=H%yval(i,:)+H2D%yval(i,j,:)*xval
             yval2(i,:) =yval2(i,:) +H2D%yval(i,j,:)
          case (2)
             xval = H2D%xMin(1)+(j*1.0-0.5)*H2D%xBin(1)
             H%yval(i,:)=H%yval(i,:)+H2D%yval(j,i,:)*xval
             yval2(i,:) =yval2(i,:) +H2D%yval(j,i,:)
          end select
       end do
    end do

    WHERE (yval2 /= 0.0) H%yval = H%yval/yval2

    deallocate(yval2)

    ! we multiply column 1,3 by bin-width, since this is divided
    ! out in the output routine:
    H%yval(:,1) = H%yval(:,1)*H%xBin
    H%yval(:,3) = H%yval(:,3)*H%xBin

  end subroutine AverageHist2D


end module hist2Df90
