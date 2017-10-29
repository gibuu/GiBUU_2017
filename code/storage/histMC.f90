!******************************************************************************
!****m* /histMC
! NAME
! module histMC
! PURPOSE
! Encapsulate all routines and data for 1D Histograms with multiple channels.
!
! Features of Histograms provided by this module:
! - store parameters of the x-binning
! - enable multiple y-values
! - track keeping of under-/over-score the given extreme values of x.
! - provide simple-to-understand output routines (cf. WriteHist)
! - provide simple histogram arithmetic (not yet implemented)
! - many others...
!
! Every Histogram prints its own multicolumn output.
! A multicolumn output of many different histograms for the same x-value
! is not implemented. This is done by the module "histMPf90".
!
! INPUTS
! ...(This module needs no input)
!******************************************************************************
module histMC

  implicit none
  private

  integer, parameter:: NameLength = 90


  !****************************************************************************
  !****t* histMC/histogramMC
  ! NAME
  ! type histogramMC
  ! PURPOSE
  ! Type definition to store all information for a multi-channel 1D Histogram.
  ! SOURCE
  !
  type, public :: histogramMC
     real             :: xMin        ! smallest x-value
     real             :: xMax        ! largest x-value
     real             :: xBin        ! width of x-Bin
     real             :: xExtreme(2) ! extremes of used x-values
                                     ! (min,max)
     character*(NameLength) :: Name  ! name to be written
     character*(NameLength) :: xDesc ! description of x-value
     character*(NameLength),allocatable :: yDesc(:) ! description of y-channels
     real,allocatable :: yVal(:,:)   ! histogramm values:
                                     ! (x), (yy1,yy2,yy3...)
                                     ! bin 0,-1 : underflow/overflow !!!
     logical :: initialized = .false.! flag to indicate whether the histogram has been initialized
  end type histogramMC
  !****************************************************************************

  public :: CreateHistMC, RemoveHistMC, clearHistMC, CopyDesc, AddHistMC
  public :: WriteHeaderMC, WriteHistMC, WriteHistMC_Integrated, &
       WriteHistMC_Gauss
  public :: ReadHistMC, FetchHistMC, sumHistMC,  rebinHistMC

contains


  !****************************************************************************
  !****s* histMC/ClearHistMC
  ! NAME
  ! subroutine ClearHistMC(H)
  ! PURPOSE
  ! Sets the histogram to zero again
  ! INPUTS
  ! * type(histogramMC) :: H         -- Histogramm to be cleared
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine ClearHistMC(H)
    type(histogramMC), intent(inout) :: H

    H%xExtreme(1:2) = (/ 99e9, -99e9 /)
    H%yVal(:,:)     = 0.

  end subroutine ClearHistMC


  !****************************************************************************
  !****s* histMC/CreateHistMC
  ! NAME
  ! subroutine CreateHistMC(H, HName, x1, x2, bin, nCh)
  ! PURPOSE
  ! This is the Constructor of a multi-channel 1D-Histogram!
  ! Allocate Memory for the entries and put additional variables to their
  ! default.
  ! INPUTS
  ! * type(histogramMC):: H         -- Histogramm to be created
  ! * character*(*)    :: HName     -- Name of Histogram
  ! * real             :: x1,x2,bin -- Minimal/maximal value for x-coordinate
  !   to be considered, bin-width
  ! * integer          :: nCh       -- Number of Channels
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  elemental subroutine CreateHistMC(H, HName, x1, x2, bin, nCh)
    type(histogramMC), intent(inout) :: H
    character*(*), intent(in) :: HName
    real, intent(in) :: x1, x2, bin
    integer, intent(in) :: nCh

    integer :: L

    H%xMin = x1
    H%xMax = x2
    H%xBin = bin
    H%xExtreme(1:2) = (/ 99e9, -99e9 /)

    if (len(HName) > NameLength) then
       H%Name = HName(1:NameLength)
    else
       H%Name = HName
    end if

    if (allocated(H%yVal)) deallocate(H%yVal)
    if (allocated(H%yDesc)) deallocate(H%yDesc)

    L = nint( (x2-x1)/bin )+1
    allocate(H%yVal(-1:L,nCh))
    allocate(H%yDesc(nCh))

    H%xDesc = ""
    H%yVal = 0.
    H%yDesc = ""
    H%initialized = .true.

  end subroutine CreateHistMC


  !****************************************************************************
  !****s* histMC/RemoveHistMC
  ! NAME
  ! subroutine RemoveHistMC(H)
  ! PURPOSE
  ! Free the allocated memory
  ! INPUTS
  ! * type(histogramMC) :: H  -- Histogramm to be used
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine RemoveHistMC(H)
    type(histogramMC),intent(inout) :: H

    if (allocated(H%yVal)) deallocate(H%yVal)
    if (allocated(H%yDesc)) deallocate(H%yDesc)
    H%initialized = .false.

  end subroutine RemoveHistMC


  !****************************************************************************
  !****s* histMC/AddHistMC
  ! NAME
  ! subroutine AddHistMC(H, x, nCh, y)
  ! PURPOSE
  ! Add to the given histogram at the given x-value the weight y (y2).
  ! If y2 is given, also the second column is filled.
  ! INPUTS
  ! * type(histogram) :: H  -- Histogramm to be used
  ! * real            :: x  -- x-value
  ! * integer         :: nCh-- number of Channel
  ! * real            :: y  -- weight to be added
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine AddHistMC(H, x, nCh, y)
    type(histogramMC)   :: H
    real, intent(in)    :: x
    integer, intent(in) :: nCh
    real, intent(in)    :: y

    integer :: iBin

    if (x >= H%xMin .and. x < H%xMax) then
       iBin = int( (x-H%xMin)/H%xBin )+1
    else if (x < H%xMin) then
       iBin = -1
    else
       iBin = 0
    end if

    if (x < H%xExtreme(1)) H%xExtreme(1)=x
    if (x > H%xExtreme(2)) H%xExtreme(2)=x

    H%yVal(iBin,nCh) = H%yVal(iBin,nCh)+y

  end subroutine AddHistMC


  !****************************************************************************
  !****s* histMC/RebinHistMC
  ! NAME
  ! subroutine RebinHistMC(H, N)
  ! PURPOSE
  ! Perform a rebinning of a present histogram, i.e. increase the bin size by an
  ! integer factor.
  ! INPUTS
  ! * type(histogramMC),intent(inout) :: H  --  histogram to be rebinned
  ! * integer, intent(in) :: N              --  rebinning factor
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine RebinHistMC(H, N)
    type(histogramMC),intent(inout) :: H
    integer, intent(in) :: N

    real, allocatable :: YY(:,:)
    integer :: i,iBinMax,nCh

    ! adjust bin width
    H%xBin = N*H%xBin

    ! do the rebinning
    iBinMax = ubound(H%yVal,dim=1)
    nCh = size(H%yVal,dim=2)

    allocate(YY(-1:iBinMax,nCh))

    YY = H%yVal

    deallocate(H%yVal)
    allocate(H%yVal(-1:iBinMax/N+1,nCh))

    H%yVal(-1,:) = YY(-1,:)  ! underflow
    H%yVal(0,:) = YY(0,:)    ! overflow  (unchanged)

    do i=1,iBinMax/N
      H%yVal(i,:) = sum( YY(N*(i-1)+1:N*i,:) , dim=1)
    end do

  end subroutine RebinHistMC


  !****************************************************************************
  !****s* histMC/CopyDesc
  ! NAME
  ! subroutine CopyDesc(H_dest,H_src)
  ! PURPOSE
  ! Copy the channel descriptions from one histogram to another.
  ! INPUTS
  ! * type(histogram) :: H_dest -- destination Histogramm
  ! * type(histogram) :: H_src  -- source Histogram
  ! OUTPUT
  ! H_dest is changed
  !****************************************************************************
  elemental subroutine CopyDesc(H_dest,H_src)
    type(histogramMC),intent(inout) :: H_dest
    type(histogramMC),intent(in) :: H_src

    H_dest%xDesc    = H_src%xDesc
    H_dest%yDesc(:) = H_src%yDesc(:)

  end subroutine CopyDesc


  !****************************************************************************
  !****s* histMC/WriteHeaderMC
  ! NAME
  ! logical function WriteHeaderMC(H, iFile, mul, FinalLineBreak)
  !
  ! PURPOSE
  ! A header is written to the file:
  ! * Underflow: sum of all calls with a x-value less than x-min
  ! * Entries: sum of all entries listed in section "Data" below
  ! * Overflow: sum of all calls with a x-value above x-max
  ! * Extrema: the smallest/biggest x-values which had been added
  ! * counted: the range of x-values, which is considered in "Entries"
  ! * descriptions of x-value and all y-channels
  !
  ! Summing "Underflow"+"Entries"+"Overflow" gives the number of ALL calls,
  ! i.e. the integral from -infty upto +infty.
  !
  ! Output of "Underflow","Entries","Overflow" shows only the sum of all channels.
  !
  ! INPUTS
  ! * type(histogramMC) :: H -- Histogramm
  ! * integer :: iFile -- file number to write to
  ! * real, optional :: mul
  ! * logical,optional :: FinalLineBreak   --
  !   .false.= No line break after the header
  ! OUTPUT
  ! * returns .true. if histogram is non-empty
  ! * header is written to file
  !****************************************************************************
  logical function WriteHeaderMC(H, iFile, mul, FinalLineBreak)
    type(histogramMC),intent(in) :: H
    integer :: iFile
    real,intent(in),optional :: mul
    logical,optional :: FinalLineBreak

    character(30) :: f
    integer :: iBin,iBinMax, nCh, iCh, maxlen
    real :: mulFak, x
    logical :: doFinalLineBreak
    real, dimension(:), allocatable :: S0,S1

    mulFak = 1.
    if (present(mul)) mulFak = mul

    if (H%xExtreme(2)<H%xExtreme(1)) then
       iBinMax = 0
    else
       iBinMax = min(ubound(H%yVal,dim=1), int((H%xExtreme(2)-H%xMin)/H%xBin)+1)
    end if
    nCh = size(H%yVal,dim=2)                     ! number of channels

    allocate (S0(0:nCh), S1(0:nCh))

    S0(1:nCh) = sum(H%yVal(1:iBinMax,:),dim=1)   ! sum of each channel
    S0(0) = sum(S0(1:nCh))                       ! total sum

    S1 = 0.
    do iBin=1,iBinMax
      x = H%xMin+(real(iBin)-0.5)*H%xBin
      S1(0) = S1(0) + sum(H%yVal(iBin,:))*x
      S1(1:nCh) = S1(1:nCh) + H%yVal(iBin,:)*x
    end do
    do iCh=0,nCh
      if (S0(iCh)>0.) then
        S1(iCh) = S1(iCh)/S0(iCh)
      else
        S1(iCh) = 0.
      end if
    end do

    write(iFile,'(a,/,a,a,/,a)')         '###','### Histogram: ', trim(H%Name), '###'
    write(iFile,'(a,ES12.4)')            '### Underflow: ', sum(H%yVal(-1,:))*mulFak
    write(iFile,'(a,ES12.4)')            '### Entries  : ', S0(0)*mulFak
    write(iFile,'(a,ES12.4)')            '### Overflow : ', sum(H%yVal(0,:))*mulFak
    write(iFile,'(a,ES11.4,a,ES11.4,a)') '###   Extrema: [',H%xExtreme(1),' ...',H%xExtreme(2),' ]'
    write(iFile,'(a,ES11.4,a,ES11.4,a)') '###   counted: [',H%xMin,' ...',H%xMax,' ]'

    if (allocated(H%yDesc)) then
      maxlen = max(24,len_trim(H%xDesc)-5,maxval(len_trim(H%yDesc)))
    else
      maxlen = 24
    end if

    write(f,'(a,i2,a)') '(a,a27,a',maxlen,',a12,a12)'
    write(iFile,fmt=f)  '###', '', '', 'integrated', ' averaged '
    write(f,'(a,i2,a)') '(a,a',5+maxlen,',2a12)'
    write(iFile,fmt=f)  '### Column  1 (X-Value): ', H%xDesc(1:maxlen+5), '----------', '----------'
    if (maxlen>24) then
      write(f,'(a,i2,a)') '(a,',maxlen-24,'X,2ES12.4)'
    else
      write(f,'(a)') '(a,2ES12.4)'
    end if
    write(iFile,fmt=f)  '### Column  2 (Y-Total): sum of all following channels', S0(0)*mulFak, S1(0)

    if (allocated(H%yDesc)) then
      write(f,'(a,i2,a)') '(a,i2,a,i2,a,a',maxlen,',2ES12.4)'
      do iCh=1,nCh
        write(iFile,fmt=f) '### Column ',iCh+2,' (Y-Channel ',iCh,'): ', H%yDesc(iCh)(1:maxlen), S0(iCh)*mulFak, S1(iCh)
      end do
    end if

    if (present(FinalLineBreak)) then
      doFinalLineBreak=FinalLineBreak
    else
      doFinalLineBreak=.true.
    end if
    if (doFinalLineBreak) then
      write(iFile,'(a,/)') '###'
    else
      write(iFile,'(a)') '###'
    end if

    WriteHeaderMC = (S0(0)>0.)    ! empty or not?

  end function WriteHeaderMC


  !****************************************************************************
  !****s* histMC/WriteHistMC
  ! NAME
  ! subroutine WriteHistMC(H,file,div,add,mul,dump)
  ! PURPOSE
  ! Write out the histogram.
  !
  ! The entries are multiplied by 'mul'.
  ! INPUTS
  ! * type(histogram) :: H     -- Histogramm to be used
  ! * character*(*)   :: file  -- name of file to open and close
  ! * real         :: div   -- if true, divide data by the bin width [OPTIONAL]
  ! * real         :: add   -- factor to add      [OPTIONAL]
  ! * real         :: mul   -- factor to multiply [OPTIONAL]
  ! * logical      :: dump  -- if true, also dump it binary to file [OPTIONAL]
  ! OUTPUT
  ! write to file
  !
  ! Several columns are written in the data section:
  ! * 1) x-value (i.e. middle of bin)
  ! * 2) total sum of all y-values
  ! * 3...) all single y-values
  !
  ! Columns are divided by the bin-width, if "div" is true (or omitted).
  ! (Default: div = .true.)
  !
  ! In addition a header with additional info is written to the file
  ! (see "WriteHeaderMC")
  !
  ! NOTES
  ! The Histogram Data is not affected!!!
  !****************************************************************************
  subroutine WriteHistMC(H,file,div,add,mul,dump)
    type(histogramMC),intent(in) :: H
    character*(*),    intent(in) :: file
    logical,          intent(in),optional :: div
    real,             intent(in),optional :: add
    real,             intent(in),optional :: mul
    logical,          intent(in),optional :: dump

    character(20) :: f
    integer       :: iBin, iBinMax, nCh
    real :: addFak, mulFak, fak
    integer,parameter:: iFile = 62

    fak = 1./H%xBin
    if (present(div) .and. .not. div) fak = 1.

    addFak = 0.
    mulFak = 1.
    if (present(add)) addFak = add
    if (present(mul)) mulFak = mul
    mulFak = mulFak*fak

    if (.not.allocated(H%yVal)) return

    open(iFile, file=file, status='unknown')
    rewind(iFile)

    write(*,*) H%xExtreme(1),H%xExtreme(2)

    if (H%xExtreme(2)<H%xExtreme(1)) then
       iBinMax = 0
    else
       iBinMax = min(ubound(H%yVal,dim=1), int((H%xExtreme(2)-H%xMin)/H%xBin)+1)
    end if
    nCh = size(H%yVal,dim=2)

    if (WriteHeaderMC(H,iFile,mulfak/fak)) then
      ! write histogram data
      write(f,'(a,i3,a)') '(',nCh+2,'ES12.4)'
      do iBin=1,iBinMax
        write(iFile,fmt=f) H%xMin+(real(iBin)-0.5)*H%xBin, &
                           sum(H%yVal(iBin,:))*mulFak+addFak, &
                           H%yVal(iBin,:)*mulFak+addFak
      end do
    else
      write(iFile,*) "(empty)"
    end if

    close(iFile)

    if (present(dump)) then
       if (dump) then
          call DumpHistMC(H,file//".bin",iFile,addFak,mulFak)
       end if
    end if

  end subroutine WriteHistMC


  !****************************************************************************
  !****s* histMC/WriteHistMC_Integrated
  ! NAME
  ! subroutine WriteHistMC_Integrated(H,file,div)
  ! PURPOSE
  ! Write out the histogram, integrating the data over x.
  ! INPUTS
  ! * type(histogram) :: H     -- Histogramm to be used
  ! * character*(*)   :: file  -- name of file to open and close
  ! * logical         :: div   -- if true, divide data by the bin width [OPTIONAL]
  ! OUTPUT
  ! write to file
  !****************************************************************************
  subroutine WriteHistMC_Integrated(H,file,div)
    type(histogramMC),intent(inout) :: H
    character*(*),intent(in) :: file
    logical,intent(in),optional :: div

    character(20) :: f
    integer       :: iBin, iBinMax, nCh
    real :: fak
    integer,parameter:: iFile = 63
    character(NameLength):: title

    fak = 1./H%xBin
    if (present(div) .and. .not. div) fak = 1.

    if (.not.allocated(H%yVal)) return

    open(iFile, file=file, status='unknown')
    rewind(iFile)

    iBinMax = min (ubound(H%yVal,dim=1), int((H%xExtreme(2)-H%xMin)/H%xBin)+1)
    nCh = size(H%yVal,dim=2)

    title = H%Name
    H%Name = trim(title) // " - Integrated"
    if (WriteHeaderMC(H,iFile)) then
      ! write histogram data
      write(f,'(a,i2,a)') '(',nCh+2,'ES12.4)'
      do iBin=1,iBinMax
        write(iFile,fmt=f) H%xMin+(real(iBin)-0.5)*H%xBin, &
                           sum(H%yVal(1:iBin,:))*fak,      &
                           sum(H%yVal(1:iBin,:),dim=1)*fak
      end do
    else
      write(iFile,*) "(empty)"
    end if

    H%Name = title
    close(iFile)

  end subroutine WriteHistMC_Integrated


  !****************************************************************************
  !****s* histMC/WriteHistMC_Spline
  ! NAME
  ! subroutine WriteHistMC_Spline(H,file,nPoints)
  ! PURPOSE
  ! Write out the histogram.
  !
  ! The binning is subdivided into "nPoints" points per bin, where
  ! every point is calculated via cubic spline interpolation.
  ! (Default: nPoints=5)
  ! INPUTS
  ! * type(histogram) :: H     -- Histogramm to be used
  ! * integer         :: iFile -- File number output to redirect
  ! * integer         :: nPoints -- number of sub-points [OPTIONAL]
  ! OUTPUT
  ! write to file number
  !
  ! cf. WriteHist
  ! NOTES
  ! The Histogram Data is not affected!!!
  !****************************************************************************
!   subroutine WriteHistMC_Spline(H,file,nPoints)
!     use cl_splines
!     type(histogramMC) :: H
!     character*(*),intent(in) :: file
!     integer,        intent(in),optional :: nPoints
!
!     integer,parameter:: iFile = 64
!     integer             :: iBin, iBinMax, nCh, ch, nPP=5,error
!     real :: x
!     real, allocatable :: XX(:),YY(:,:), y(:)
!     logical :: success
!     character(20) :: f
!     character(NameLength):: title
!
!     type(tspline),allocatable :: spline(:)
!
!     nPP = 5
!
!     if (present(nPoints)) nPP = nPoints
!
!      if (.not.allocated(H%yVal)) return
!
!     open(iFile, file=file, status='unknown')
!     rewind(iFile)
!
!     iBinMax = ubound(H%yVal,dim=1)
!     nCh = size(H%yVal,dim=2)
!
!     ! calculate splines
!     allocate(XX(iBinMax))
!     allocate(YY(iBinMax,0:nCh))
!     allocate(spline(0:nCh))
!
!     do iBin=1,iBinMax
!        xx(iBin)  = H%xMin+(real(iBin)-0.5)*H%xBin
!        YY(iBin,0) = sum(H%yVal(iBin,:))/H%xBin
!        do ch=1,nCh
!          YY(iBin,ch) = H%yVal(iBin,ch)/H%xBin
!        end do
!     end do
!
!     do ch=0,nCh
!       spline(ch)=cl_initSpline(xx,YY(:,nch))
!     end do
!
!     title = H%Name
!     H%Name = trim(title) // " - Splined"
!     call WriteHeaderMC(H,iFile)
!     H%Name = title
!
!     ! write histogram data
!     allocate(y(0:nCh))
!     write(f,'(a,i2,a)') '(',nCh+2,'ES12.4)'
!     do iBin=1,iBinMax*nPP
!        x = H%xMin+(real(iBin)-0.5)*H%xBin/nPP
!        do ch=0,nCh
!          y(ch)=cl_spline(spline(ch),x,success,error)
!        end do
!        write(iFile,fmt=f) x,y
!     end do
!     close(iFile)
!     deallocate(XX,YY,y)
!
!   end subroutine WriteHistMC_Spline


  !****************************************************************************
  !****s* histMC/WriteHistMC_Gauss
  ! NAME
  ! subroutine WriteHistMC_Gauss(H, file, width_in, mul)
  ! PURPOSE
  ! Write out the histogram, convoluted with a gaussian distribution.
  !
  ! INPUTS
  ! * type(histogram) :: H         -- Histogramm to be used
  ! * character*(*)   :: file      -- name of file to open and close
  ! * real, optional  :: width_in  -- width of gaussian
  ! * real, optional  :: mul       -- factor to multiply
  ! OUTPUT
  ! write to file
  !
  ! NOTES
  ! The histogram data is not affected!!!
  !****************************************************************************
  subroutine WriteHistMC_Gauss(H, file, width_in, mul)
    use distributions, only: gauss
    type(histogramMC) :: H
    character*(*), intent(in) :: file
    real, intent(in), optional :: width_in, mul

    character(20) :: f
    integer       :: i, iBinMax, nCh,j
    integer,parameter:: iFile = 65
    real,allocatable:: S(:)
    integer, parameter :: N = 1
    real :: xi,xj,width,mulFak
    character(NameLength):: title

    ! set width of gaussian
    if (present(width_in)) then
      width = width_in
    else
      width = H%xBin ! default: bin width
    end if

    if (present(mul)) then
      mulFak = mul
    else
      mulFak = 1.
    end if

    if (.not.allocated(H%yVal)) return

    open(iFile, file=file, status='unknown')
    rewind(iFile)

    iBinMax = min (ubound(H%yVal,dim=1), int((H%xExtreme(2)-H%xMin)/H%xBin)+1)
    nCh = size(H%yVal,dim=2)

    allocate(S(nCh))

    title = H%Name
    H%Name = trim(title) // " - Gauss-Convoluted"
    if (WriteHeaderMC(H,iFile,mulFak)) then
      ! write histogram data
      write(f,'(a,i2,a)') '(',nCh+2,'ES12.4)'
      do i=1,iBinMax*N
         xi = H%xMin + (i+N/2.-N/2-1)*H%xBin/N
         S = 0.
         do j=max(1,nint((i-10*N*width/H%xBin)/N)),min(iBinMax,nint((i+10*N*width/H%xBin)/N))
            xj = H%xMin + (real(j)-0.5)*H%xBin
            S = S + H%yVal(j,:) * mulFak * gauss(xi,xj,width)
         end do
         write(iFile,fmt=f) xi, sum(S), S
      end do
    else
      write(iFile,*) "(empty)"
    end if

    H%Name = title
    close(iFile)

  end subroutine WriteHistMC_Gauss


  !****************************************************************************
  !****f* histMC/ReadHistMC
  ! NAME
  ! function ReadHistMC(H, file, w) result (success)
  ! PURPOSE
  ! Read in a histogram from a file.
  ! INPUTS
  ! * type(histogramMC) :: H            -- Histogram to be used
  ! * character*(*)     :: file         -- name of file to read
  ! * real, optional    :: w   -- optional weigting factor
  ! OUTPUT
  ! * H is changed
  ! * return value indicates success or failure
  ! NOTES
  ! Only data points are read (no under-/overflow etc).
  !****************************************************************************
  function ReadHistMC(H, file, w) result (success)
    type(histogramMC) :: H
    character*(*),  intent(in) :: file
    real, optional, intent(in) :: w
    logical :: success

    integer, parameter :: iFile = 68
    character(len=300) :: line
    integer :: i,ios, head, numbins, nCh, linelen=-1
    real :: x,s
    character(20) :: f

    success = .false.

    if (H%initialized) call RemoveHistMC(H)

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

    ! correct for 1 blank line
    head = head + 1
    numbins = numbins - 1

    ! determine number of channels
    nCh = head - 13

    allocate(H%yVal(-1:numbins,nCh))
    allocate(H%yDesc(nCh))
    H%yVal = 0.

    rewind(iFile)
    ! discard header, read histogram name, read column descriptions
    do i=1,head
      read (iFile,'(A)',iostat=ios) line
      if (ios /= 0) then
        close(iFile)
        call RemoveHistMC(H)
        return
      end if
      if (i==2) then
        H%name = line(index(line,":")+2:)
      else if (i==10) then
        linelen = index(line,"----------")-1
        H%xDesc = line(index(line,":")+2:linelen)
      else if (i>11 .and. i<=11+nCh) then
        H%yDesc(i-11)=line(index(line,":")+2:linelen)
      end if
    end do

    ! read values
    write(f,'(a,i2,a)') '(',nCh+2,'ES12.4)'
    do i=1,numbins
      read (iFile,f,iostat=ios) x,s,H%yval(i,:)
      if (ios /= 0) then
        close(iFile)
        call RemoveHistMC(H)
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
    H%xExtreme(1:2) = (/ H%xMin, H%xMax /)  ! TODO: read extreme values from file header?
    ! correct for bin width
    H%yval(:,:) = H%yval(:,:) * H%xbin
    ! scale with weighting factor
    if (present(w)) H%yval(:,:) = H%yval(:,:) * w

    close(iFile)
    H%initialized = .true.

    success = .true.

  end function ReadHistMC


  !****************************************************************************
  !****s* histMC/sumHistMC
  ! NAME
  ! subroutine sumHistMC(A, B, w)
  ! PURPOSE
  ! Perform a summation of two histograms: A = A + B.
  ! If a weight 'w' is given, do a weighted summation: A = A + w*B
  ! INPUTS
  ! * type(histogramMC) :: A, B         -- Histograms to be used
  ! * real, optional    :: w   -- optional weighting factor
  ! OUTPUT
  ! A is changed.
  !****************************************************************************
  subroutine sumHistMC(A, B, w)
    type(histogramMC) :: A
    type(histogramMC), intent(in) :: B
    real, optional, intent(in) :: w

    real :: wgt
    integer :: nA, nB, nMin, nMax
    real, dimension(:,:), allocatable :: tmp

    if (present(w)) then
      wgt = w
    else
      wgt = 1.
    end if

    ! determine number of bins for both histograms
    nA = ubound(A%yVal,dim=1)
    nB = ubound(B%yVal,dim=1)
    nMin = min(nA, nB)
    nMax = max(nA, nB)

    if (nA>=nB) then
      ! A is large enough
      A%yVal(:nMin,:) = A%yVal(:nMin,:) + wgt * B%yVal(:nMin,:)
    else
      ! A is too small, must be resized
      call move_alloc(A%yVal, tmp)
      allocate (A%yVal(-1:nMax,lbound(tmp,dim=2):ubound(tmp,dim=2)))
      A%yVal(:nMin,:) = tmp(:nMin,:) + wgt * B%yVal(:nMin,:)
      A%yVal(nMin+1:nMax,:) = wgt * B%yVal(nMin+1:nMax,:)
      deallocate (tmp)
      A%xMax = B%xMax
    end if

    A%xExtreme(1) = min(A%xExtreme(1), B%xExtreme(1))
    A%xExtreme(2) = max(A%xExtreme(2), B%xExtreme(2))

  end subroutine sumHistMC


  !****************************************************************************
  !****s* histMC/DumpHistMC
  ! NAME
  ! subroutine DumpHistMC(H,file,iFile, add,mul)
  ! PURPOSE
  ! Write all the histogram information unformatted (i.e. binary) to a file
  !
  ! INPUTS
  ! * type(histogramMC) :: H     -- Histogramm to be used
  ! * character*(*)     :: file  -- name of file to open and close
  ! * integer,OPTIONAL  :: iFile -- File number output to redirect [OPTIONAL]
  ! * real            :: add   -- factor to add      [OPTIONAL]
  ! * real            :: mul   -- factor to multiply [OPTIONAL]
  ! OUTPUT
  ! H is written UNFORMATTED to the given file
  !
  !****************************************************************************
  subroutine DumpHistMC(H,file,iFile, add,mul)
    type(histogramMC),intent(in)          :: H
    character*(*),    intent(in)          :: file
    integer,          intent(in),optional :: iFile
    real,             intent(in),optional :: add,mul

    integer :: iF
    real :: addFak,mulFak
    logical :: WriteFaks

    iF=121
    if (present(iFile)) iF = iFile

    open(iF,file=file,status='UNKNOWN',form='UNFORMATTED')
    rewind(iF)

    write(iF) H%xMin,H%xMax,H%xBin,H%xExtreme
    write(iF) H%Name
    write(iF) H%xDesc
    write(iF) size(H%yDesc)
    write(iF) H%yDesc
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

  end subroutine DumpHistMC

  !****************************************************************************
  !****s* histMC/FetchHistMC
  ! NAME
  ! subroutine FetchHistMC(H,file,iFile, add,mul,flagOK)
  ! PURPOSE
  ! Read in all the histogram information previously dumped unformatted
  ! (i.e. binary) to a file
  !
  ! INPUTS
  ! * character*(*)   :: file  -- name of file to open and close
  ! * integer,OPTIONAL:: iFile -- File number input to redirect [OPTIONAL]
  ! OUTPUT
  ! * type(histogramMC) :: H     -- Histogramm to be used
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
  subroutine FetchHistMC(H,file,iFile, add,mul, flagOK)
    type(histogramMC),intent(inout)       :: H
    character*(*),    intent(in)          :: file
    integer,          intent(in),optional :: iFile
    real,             intent(out),optional:: add,mul
    logical,          intent(out),optional:: flagOK

    integer :: iF, nCh, L
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
    read(iF) H%xDesc
    read(iF) nCh

    if (allocated(H%yVal)) deallocate(H%yVal)
    if (allocated(H%yDesc)) deallocate(H%yDesc)

    L  = nint( (H%xMax-H%xMin)/H%xBin )+1
    allocate(H%yVal(-1:L,nCh))
    allocate(H%yDesc(nCh))

    read(iF) H%yDesc
    read(iF) H%yVal

    H%initialized = .true.

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

  end subroutine FetchHistMC


  !****************************************************************************
end module histMC
