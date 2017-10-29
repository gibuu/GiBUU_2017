!******************************************************************************
!****m* /histMC_avg
! NAME
! module histMC_avg
! PURPOSE
! This module implements a type 'histogramMC_avg' for 1D Histograms with
! multiple channels. It is similar to the type 'histogramMC' but in addition
! stores event counts for each bin and channel. In the end it outputs the
! average of the stored quantity over all entries in a particular bin.
! This is useful e.g. for the flow analysis in the module 'HeavyIonAnalysis'.
!******************************************************************************
module histMC_avg

  implicit none
  private

  integer, parameter:: NameLength = 90


  !****************************************************************************
  !****t* histMC_avg/histogramMC_avg
  ! NAME
  ! type histogramMC_avg
  ! PURPOSE
  ! Type definition to store all information for a multi-channel 1D Histogram
  ! with averaging capabilities.
  ! SOURCE
  !
  type, public :: histogramMC_avg
    real :: xMin, xMax, xBin                         ! smallest/largest x-value & bin width
    real :: xExtreme(2)                              ! extremes of used x-values (min,max)
    character*(NameLength) :: Name, xDesc            ! histogram name and description of x-value
    character*(NameLength), allocatable :: yDesc(:)  ! description of y-channels
    real, allocatable :: yVal(:,:)      ! histogramm values: (x), (y1,y2,y3,...)
                                        ! bin 0,-1 : underflow/overflow !!!
    integer, allocatable :: counts(:,:) ! counts: (x), (count1,count2,count3,...)
                                        ! bin 0,-1 : underflow/overflow counts !!!
    logical :: initialized = .false.    ! flag to indicate whether the histogram has been initialized
  end type histogramMC_avg
  !****************************************************************************

  public :: CreateHistMC_avg, CopyDesc_avg, AddHistMC_avg, WriteHistMC_avg
  public :: ReadHistMC_avg, sumHistMC_avg

contains


  !****************************************************************************
  !****s* histMC_avg/CreateHistMC_avg
  ! NAME
  ! subroutine CreateHistMC_avg (H, HName, x1, x2, bin, nCh)
  ! PURPOSE
  ! This is the constructor of a multi-channel 1D-Histogram.
  ! Allocate memory for the entries and put additional variables to their default.
  ! INPUTS
  ! * type(histogramMC_avg):: H     -- Histogramm to be created
  ! * character*(*)    :: HName     -- Name of histogram
  ! * real             :: x1,x2,bin -- Minimal/maximal value for x-coordinate to be considered, bin-width
  ! * integer          :: nCh       -- Number of channels
  ! OUTPUT
  ! H is changed.
  !****************************************************************************
  elemental subroutine CreateHistMC_avg (H, HName, x1, x2, bin, nCh)
    type(histogramMC_avg), intent(inout) :: H
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
    allocate(H%counts(-1:L,nCh))
    allocate(H%yDesc(nCh))

    H%yVal = 0.
    H%counts = 0
    H%xDesc = ""
    H%yDesc = ""
    H%initialized = .true.

  end subroutine CreateHistMC_avg


  !****************************************************************************
  !****s* histMC_avg/histMC_avg
  ! NAME
  ! subroutine RemoveHistMC_avg(H)
  ! PURPOSE
  ! Free the allocated memory.
  ! INPUTS
  ! * type(histogramMC_avg) :: H  -- Histogram to be used
  ! OUTPUT
  ! H is changed.
  !****************************************************************************
  subroutine RemoveHistMC_avg(H)
    type(histogramMC_avg),intent(inout) :: H

    if (allocated(H%yVal)) deallocate(H%yVal)
    if (allocated(H%yDesc)) deallocate(H%yDesc)
    if (allocated(H%counts)) deallocate(H%counts)
    H%initialized = .false.

  end subroutine RemoveHistMC_avg


  !****************************************************************************
  !****s* histMC_avg/AddHistMC_avg
  ! NAME
  ! subroutine AddHistMC_avg (H, x, nCh, y)
  ! PURPOSE
  ! Add the weight y to the given histogram at the given x-value and channel number,
  ! and increment the count for that bin.
  ! INPUTS
  ! * type(histogramMC_avg) :: H   -- Histogramm to be used
  ! * real :: x                    -- x-value
  ! * integer :: nCh               -- number of Channel
  ! * real :: y                    -- weight to be added
  ! OUTPUT
  ! H is changed.
  !****************************************************************************
  subroutine AddHistMC_avg (H, x, nCh, y)
    type(histogramMC_avg) :: H
    real,    intent(in)   :: x
    integer, intent(in)   :: nCh
    real,    intent(in)   :: y

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
    H%counts(iBin,nCh) = H%counts(iBin,nCh)+1

  end subroutine AddHistMC_avg


  !****************************************************************************
  !****s* histMC_avg/CopyDesc_avg
  ! NAME
  ! subroutine CopyDesc_avg (H_dest, H_src)
  ! PURPOSE
  ! Copy the channel descriptions from one histogram to another.
  ! INPUTS
  ! * type(histogramMC_avg) :: H_dest -- destination Histogramm
  ! * type(histogramMC_avg) :: H_src  -- source Histogram
  ! OUTPUT
  ! H_dest is changed
  !****************************************************************************
  elemental subroutine CopyDesc_avg (H_dest, H_src)
    type(histogramMC_avg),intent(inout) :: H_dest
    type(histogramMC_avg),intent(in) :: H_src

    H_dest%xDesc    = H_src%xDesc
    H_dest%yDesc(:) = H_src%yDesc(:)

  end subroutine CopyDesc_avg


  !****************************************************************************
  !****s* histMC_avg/WriteHeader
  ! NAME
  ! logical function WriteHeader (H, iFile, mul, FinalLineBreak)
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
  ! * type(histogramMC_avg) :: H -- Histogramm
  ! * integer :: iFile -- file number to write to
  ! * real, optional :: mul
  ! * logical,optional :: FinalLineBreak   -- .false.= No line break after the header
  ! OUTPUT
  ! * returns .true. if histogram is non-empty
  ! * header is written to file
  !****************************************************************************
  logical function WriteHeader (H, iFile, mul, FinalLineBreak)
    type(histogramMC_avg),intent(in) :: H
    integer :: iFile
    real,intent(in),optional :: mul
    logical,optional :: FinalLineBreak

    character(30) :: f
    integer :: iBin,iBinMax, nCh, iCh, maxlen
    real :: mulFak, x
    logical :: doFinalLineBreak
    real, dimension(:), allocatable :: S0,S1
    integer, dimension(:), allocatable :: cnt

    mulFak = 1.
    if (present(mul)) mulFak = mul

    iBinMax = min (ubound(H%yVal,dim=1), int((H%xExtreme(2)-H%xMin)/H%xBin)+1)
    nCh = size(H%yVal,dim=2)                     ! number of channels

    allocate (S0(0:nCh), S1(0:nCh), cnt(0:nCh))

    S0(1:nCh) = sum(H%yVal(1:iBinMax,:),dim=1)     ! sum of each channel
    S0(0) = sum(S0(1:nCh))                         ! total sum
    cnt(1:nCh) = sum(H%counts(1:iBinMax,:),dim=1)  ! counts per channel
    cnt(0) = sum(cnt(1:nCh))                       ! total counts

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

    WriteHeader = (cnt(0)>0)    ! empty or not?

  end function WriteHeader


  !****************************************************************************
  !****s* histMC_avg/WriteHistMC_avg
  ! NAME
  ! subroutine WriteHistMC_avg(H,file,mul)
  ! PURPOSE
  ! Write out the histogram. The entries are multiplied by 'mul'.
  ! INPUTS
  ! * type(histogramMC_avg) :: H   -- Histogramm to be used
  ! * character*(*) :: file        -- name of file to open and close
  ! * real :: mul                  -- factor to multiply [OPTIONAL]
  ! OUTPUT
  ! write to file
  !
  ! Several columns are written in the data section:
  ! * 1) x-value (i.e. middle of bin)
  ! * 2) total sum of all y-values
  ! * 3...) all single y-values (followed by their counts)
  !
  ! In addition a header with additional info is written to the file (see "WriteHeader").
  !
  ! NOTES
  ! The Histogram data is not affected.
  !****************************************************************************
  subroutine WriteHistMC_avg(H,file,mul)
    type(histogramMC_avg), intent(in) :: H
    character*(*), intent(in) :: file
    real, intent(in), optional :: mul

    character(20) :: f
    integer :: iBin, iBinMax, nCh, cnt
    real :: mulFak
    integer, parameter:: iFile = 62
    real, dimension(:), allocatable :: avg

    mulFak = 1.
    if (present(mul)) mulFak = mul

    if (.not.allocated(H%yVal)) return

    open(iFile, file=file, status='unknown')
    rewind(iFile)

    iBinMax = min (ubound(H%yVal,dim=1), int((H%xExtreme(2)-H%xMin)/H%xBin))
    nCh = size(H%yVal,dim=2)

    allocate (avg(0:nCh))

    if (WriteHeader(H,iFile,mulfak)) then
      ! write histogram data
      write(f,'(a,i2,a,i2,a)') '(',nCh+2,'ES12.4,',nCh,'I7)'
      do iBin=1,iBinMax
        cnt = sum(H%counts(iBin,:))
        if (cnt>0) then
          avg(0) = sum(H%yVal(iBin,:))*mulFak/cnt
          avg(1:nCh) = H%yVal(iBin,:)*mulFak/H%counts(iBin,:)
        else
          avg = 0.
        end if
        write(iFile,fmt=f) H%xMin+(real(iBin)-0.5)*H%xBin, &
                           avg(0:nCh), H%counts(iBin,:)
      end do
    else
      write(iFile,*) "(empty)"
    end if

    close(iFile)

  end subroutine WriteHistMC_avg


  !****************************************************************************
  !****f* histMC_avg/ReadHistMC_avg
  ! NAME
  ! function ReadHistMC_avg (H, file, w) result (success)
  ! PURPOSE
  ! Read in a histogram from a file.
  ! INPUTS
  ! * type(histogramMC_avg) :: H        -- Histogram to be used
  ! * character*(*)     :: file         -- name of file to read
  ! * real, optional, intent(in) :: w   -- optional weigting factor
  ! OUTPUT
  ! * H is changed
  ! * return value indicates success or failure
  ! NOTES
  ! Only data points are read (no under-/overflow etc).
  !****************************************************************************
  function ReadHistMC_avg (H, file, w) result (success)
    type(histogramMC_avg) :: H
    character*(*),  intent(in) :: file
    real, optional, intent(in) :: w
    logical :: success

    integer, parameter :: iFile = 68
    character(len=300) :: line
    integer :: i,ios, head, numbins, nCh, linelen=-1
    real :: x,s
    character(20) :: f

    success = .false.

    if (H%initialized) call RemoveHistMC_avg(H)

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
    allocate(H%counts(-1:numbins,nCh))
    allocate(H%yDesc(nCh))
    H%yVal = 0.

    rewind(iFile)
    ! discard header, read histogram name, read column descriptions
    do i=1,head
      read (iFile,'(A)',iostat=ios) line
      if (ios /= 0) then
        close(iFile)
        call RemoveHistMC_avg(H)
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
    write(f,'(a,i2,a,i2,a)') '(',nCh+2,'ES12.4,',nCh,'I7)'
    do i=1,numbins
      read (iFile,f,iostat=ios) x,s,H%yval(i,:),H%counts(i,:)
      if (ios /= 0) then
        close(iFile)
        call RemoveHistMC_avg(H)
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
    ! correct for counts
    H%yval(:,:) = H%yval(:,:) * H%counts(:,:)
    ! scale with weighting factor
    if (present(w)) H%yval(:,:) = H%yval(:,:) * w

    close(iFile)
    H%initialized = .true.

    success = .true.

  end function ReadHistMC_avg


  !****************************************************************************
  !****s* histMC_avg/sumHistMC_avg
  ! NAME
  ! subroutine sumHistMC_avg (A, B, w)
  ! PURPOSE
  ! Perform a summation of two histograms: A = A + B.
  ! If a weight 'w' is given, do a weighted summation: A = A + w*B
  ! INPUTS
  ! * type(histogramMC_avg) :: A, B     -- Histograms to be used
  ! * real, optional, intent(in) :: w   -- optional weighting factor
  ! OUTPUT
  ! A is changed.
  !****************************************************************************
  subroutine sumHistMC_avg (A, B, w)
    type(histogramMC_avg) :: A
    type(histogramMC_avg), intent(in) :: B
    real, optional, intent(in) :: w

    real :: wgt
    integer :: nA, nB, nMin, nMax
    real, dimension(:,:), allocatable :: tmp_val
    integer, dimension(:,:), allocatable :: tmp_cnt

    if (present(w)) then
      wgt = w
    else
      wgt = 1.
    end if

    ! determine number of bins for both histograms
    nA = ubound(A%yVal,dim=1)
    nB = ubound(B%yVal,dim=1)
    nMin = min (nA, nB)
    nMax = max (nA, nB)

    if (nA>=nB) then
      ! A is large enough
      A%yVal(:nMin,:) = A%yVal(:nMin,:) + wgt * B%yVal(:nMin,:)
      A%counts(:nMin,:) = A%counts(:nMin,:) + int(wgt * B%counts(:nMin,:))
    else
      ! A is too small, must be resized
      ! (1) y-values
      call move_alloc (A%yVal, tmp_val)
      allocate (A%yVal(-1:nMax,lbound(tmp_val,dim=2):ubound(tmp_val,dim=2)))
      A%yVal(:nMin,:) = tmp_val(:nMin,:) + wgt * B%yVal(:nMin,:)
      A%yVal(nMin+1:nMax,:) = wgt * B%yVal(nMin+1:nMax,:)
      deallocate (tmp_val)
      ! (2) counts
      call move_alloc (A%counts, tmp_cnt)
      allocate (A%counts(-1:nMax,lbound(tmp_cnt,dim=2):ubound(tmp_cnt,dim=2)))
      A%counts(:nMin,:) = tmp_cnt(:nMin,:) + int(wgt * B%counts(:nMin,:))
      A%counts(nMin+1:nMax,:) = int(wgt * B%counts(nMin+1:nMax,:))
      deallocate (tmp_cnt)
      A%xMax = B%xMax
    end if

    A%xExtreme(1) = min (A%xExtreme(1), B%xExtreme(1))
    A%xExtreme(2) = max (A%xExtreme(2), B%xExtreme(2))

  end subroutine sumHistMC_avg


end module histMC_avg
