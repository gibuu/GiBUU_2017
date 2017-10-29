!******************************************************************************
!****m* /hist2DMPf90
! NAME
! module hist2DMPf90
! PURPOSE
! Provide an array of 2D histogramms (module hist2Df90) according the
! multi particle scheme (module histMPf90)
! NOTES
! This is histogram2D with one additional dimension (=particle, encoded)
!******************************************************************************
module hist2DMPf90
  use histf90, only: NameLength
  implicit none

  public :: histogram2DMP

  !****************************************************************************
  !****t* hist2DMPf90/histogram2DMP
  ! NAME
  ! type histogram2DMP
  ! PURPOSE
  ! Type definition to store all information for 2D Histograms for multiple
  ! particles.
  ! SOURCE
  !
  type histogram2DMP
     real             :: xMin(2)        ! smallest x-value: (x1,x2)
     real             :: xMax(2)        ! largest x-value:  (x1,x2)
     real             :: xBin(2)        ! width of x-Bin:   (x1,x2)
     real,allocatable :: xExtreme(:,:,:)! extremes of used x-values
                                        !  (id), (x1,x2), (min,max)
     character*(NameLength) :: Name     ! name to be written
     real,allocatable :: yVal(:,:,:,:)  ! histogramm values:
                                        !  (id), (x1),(x2),(yy1,yy2,yy3)
     integer          :: iSet           ! particle set
  end type histogram2DMP
  !****************************************************************************

  private

  logical, save :: initFlag=.true.

  public :: CreateHist2DMP, AddHist2DMP
  public :: WriteHist2DMP_Gnuplot

contains

  !****************************************************************************
  !****s* hist2DMPf90/CreateHist2DMP
  ! NAME
  ! subroutine CreateHist2DMP(H, HName,x1,x2,bin, iPartSet,verbose)
  ! PURPOSE
  ! This is the Constructor of a multi-particle 2D-Histogram!
  ! Allocate Memory for the entries and put additional variables to their
  ! default.
  ! The parameter iPartSet specifies, which particles should be included.
  !
  ! INPUTS
  ! * type(histogram2DMP) :: H         -- 2D-Histogramm-Array to be created
  ! * character*(*)       :: HName     -- Name of Histogram
  ! * real,dimension(2)   :: x1,x2,bin -- Minimal/maximal value for x-coordinate
  !                                     to be considered, bin-width
  ! * integer             :: iPartSet  -- particle set to consider
  ! * logical             :: verbose   -- switch for verbosity [OPTIONAL]
  !
  ! possible values for "iPartSet" are:
  ! * 1: pions, kaons, nucleons (i.e. "stable" particles)
  ! * 2: all mesons (except charm), nucleons, Deltas
  ! * ... to be continued
  !
  !
  ! OUTPUT
  ! ---
  !****************************************************************************
  subroutine CreateHist2DMP(H, HName,x1,x2,bin, iPartSet,verbose)
    use histMPf90

    type(histogram2DMP),intent(inout) :: H
    character*(*),intent(in) :: HName
    real,dimension(2),intent(in) :: x1
    real,dimension(2),intent(in) :: x2
    real,dimension(2),intent(in) :: bin
    integer, intent(in) :: iPartSet
    logical, intent(in),optional :: verbose

    integer,dimension(2) :: L
    integer :: nHist


    if (initFlag) then
       call histMPf90_Init
       initFlag = .false.
    end if

    nHist = Map2HistMP_getN(iPartSet)
    if (nHist==0) then
       write(*,*) 'CreateHist2DMP: iPartSet invalid:',iPartSet
       stop
    end if
    H%iSet = iPartSet
    H%xMin = x1
    H%xMax = x2
    H%xBin = bin

    if (allocated(H%xExtreme)) deallocate(H%xExtreme)
    if (allocated(H%yVal)) deallocate(H%yVal)

    allocate(H%xExtreme(nHist,2,2))
    H%xExtreme(:,1:2,1) =  99e9
    H%xExtreme(:,1:2,2) = -99e9

    if (len(HName) > NameLength) then
       H%Name = HName(1:NameLength)
    else
       H%Name = HName
    end if

    L = nint( (x2-x1)/bin )+1
    allocate(H%yVal(nHist,-1:L(1),-1:L(2),3))

    H%yVal = 0.

    if (present(verbose)) then
       if (verbose) then
          write(*,*) '***** CreateHist2DMP:',hName,L
          write(*,*) '***** --iSet,nHist  :',iPartSet,nHist
          write(*,*) '***** --x1          :',x1
          write(*,*) '***** --x2          :',x2
          write(*,*) '***** --bin         :',bin
       end if
    end if

  end subroutine CreateHist2DMP

  !****************************************************************************
  !****s* hist2DMPf90/AddHist2DMP
  ! NAME
  ! subroutine AddHist2DMP(H, Part, x,y,y2)
  ! PURPOSE
  ! Add to the given histogram at the given x-value the weight y (y2).
  ! If y2 is given, also the second column is filled.
  ! ID and IC specify the particle.
  ! x is a 2D value.
  !
  ! INPUTS
  ! * type(histogram2DMP) :: H  -- Histogramm to be used
  ! * type(particle)    :: Part -- Particle
  ! * real,dimension(2) :: x  -- (x1,x2)-value
  ! * real              :: y  -- weight to be added
  ! * real              :: y2 -- second weight to be added [OPTIONAL]
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine AddHist2DMP(H, Part, x,y,y2)
    use histMPf90
    use particleDefinition

    type(histogram2DMP),intent(inout) :: H
    type(particle),intent(in) :: Part
    real,dimension(2), intent(in)   :: x
    real,intent(in)                 :: y
    real,intent(in),optional        :: y2

    integer :: iH
    real :: yy

    yy = 0.
    if (present(y2)) yy = y2

    iH = Map2HistMP(Part, H%iSet)
    if (iH < 1) return

    call AddHist2DMP_iH(H,iH,x,y,yy)

  end subroutine AddHist2DMP

  !****************************************************************************
  !****s* hist2DMPf90/AddHist2DMP_IH
  ! NAME
  ! subroutine AddHist2DMP_IH(H, iH, x,y,y2)
  ! PURPOSE
  ! Add to the given histogram at the given x-value the weight y (y2).
  ! If y2 is given, also the second column is filled.
  ! ID and IC specify the particle.
  ! x is a 2D value.
  !
  ! INPUTS
  ! * type(histogram2DMP) :: H  -- Histogramm to be used
  ! * integer             :: iH -- nr of Histogram
  ! * real,dimension(2)   :: x  -- (x1,x2)-value
  ! * real                :: y  -- weight to be added
  ! * real                :: y2 -- second weight to be added [OPTIONAL]
  ! OUTPUT
  ! H is changed
  ! NOTES
  ! this is an internal routine
  !****************************************************************************
  subroutine AddHist2DMP_IH(H, iH, x,y,y2)

    type(histogram2DMP),intent(inout) :: H
    integer, intent(in)               :: iH
    real,dimension(2), intent(in)     :: x
    real,intent(in)                   :: y
    real,intent(in),optional          :: y2

    integer,dimension(2) :: iBin
    real :: yy
    integer :: i

    if (iH < 1) return

    yy = 0.
    if (present(y2)) yy = y2


    do i=1,2

       !...extremes (counted):

       if (x(i).lt.H%xExtreme(iH,i,1)) H%xExtreme(iH,i,1)=x(i)
       if (x(i).gt.H%xExtreme(iH,i,2)) H%xExtreme(iH,i,2)=x(i)

       !...extremes (measured):

       if (x(i) < H%xMin(i)) then
          iBin(i) = -1
       else if (x(i) > H%xMax(i)) then
          iBin(i) =  0
       else
          iBin(i) = int( (x(i)-H%xMin(i))/H%xBin(i) )+1
       end if

    end do

    H%yVal(iH,iBin(1),iBin(2),1) = H%yVal(iH,iBin(1),iBin(2),1)+y
    H%yVal(iH,iBin(1),iBin(2),2) = H%yVal(iH,iBin(1),iBin(2),2)+1.
    H%yVal(iH,iBin(1),iBin(2),3) = H%yVal(iH,iBin(1),iBin(2),3)+yy

  end subroutine AddHist2DMP_IH

  !****************************************************************************
  !****s* hist2DMPf90/WriteHist2DMP_Gnuplot
  ! NAME
  ! subroutine WriteHist2DMP_Gnuplot(H,iFile,add,mul, iColumn,DoAve,MaxVal,file,dump)
  ! PURPOSE
  ! Write out the 2D-histogram. Format is suitable as input to
  ! gnuplots "splot" command
  !
  ! The entries are multiplied by 'mul' and 'add' is added.
  !
  ! If iColumn is not given, all 3 columns are writen.
  !
  ! INPUTS
  ! * type(histogram2D) :: H     -- Histogramm to be used
  ! * integer      :: iFile -- File number output to redirect [OPTIONAL]
  ! * real         :: add   -- factor to add      [OPTIONAL]
  ! * real         :: mul   -- factor to multiply [OPTIONAL]
  ! * integer      :: iColumn -- Column to write  [OPTIONAL]
  ! * logical      :: DoAve -- write also "average" (cf Notes) [OPTIONAL]
  ! * real         :: maxVal -- value used instead "divide by zero" [OPTIONAL]
  ! * character*(*)   :: file  -- name of file to open and close [OPTIONAL]
  ! * logical         :: dump  -- if true, also dump it binary to file [OPTIONAL]
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
  ! This routine does the same as WriteHist2D_Gnuplot for every single
  ! 2D-Histogram, but now looping over all stored particles.
  !
  !
  !****************************************************************************
  subroutine WriteHist2DMP_Gnuplot(H,iFile_in,add,mul, iColumn,DoAve,MaxVal,H2,file,dump)
    use histf90
    use histMPf90

    type(histogram2DMP),intent(in)        :: H
    integer,          intent(in),optional :: iFile_in
    real,             intent(in),optional :: add
    real,             intent(in),optional :: mul
    integer,          intent(in),optional :: iColumn
    logical,          intent(in),optional :: DoAve
    real,             intent(in),optional :: maxVal
    type(histogram),  intent(in),optional :: H2
    character*(*),    intent(in),optional :: file
    logical,          intent(in),optional :: dump

    integer, dimension(2) :: iBinMax
    integer :: iBin1,iBin2, iH, iFile

    real :: addFak
    real :: mulFak
    real :: Z, maxZ
    logical :: writeZ

    real,allocatable :: yVal(:,:,:)

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

    mulFak = mulFak / (H%xBin(1)*H%xBin(2))

    iBinMax(1) = ubound(H%yVal,dim=2)
    iBinMax(2) = ubound(H%yVal,dim=3)

    allocate(yVal(-1:iBinMax(1),-1:iBinMax(2),3))

    do iH=1,Map2HistMP_getN(H%iSet)

       yVal = H%yVal(iH,:,:,:)

       if (present(H2)) then
          if (ubound(H2%yVal,dim=1).ne.iBinMax(1)) then
             write(*,*) 'WriteHist2DMP_Gnuplot: ERROR, dim not equal! STOP'
             write(*,*) '  histogram=',H%Name
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
                   write(iFile,1000) H%xMin+((/iBin1,iBin2/)*1.0-0.5)*H%xBin, Z, &
                        yVal(iBin1,iBin2,iColumn)*MulFak+AddFak
                else
                   write(iFile,1000) H%xMin+((/iBin1,iBin2/)*1.0-0.5)*H%xBin, &
                        yVal(iBin1,iBin2,iColumn)*MulFak+AddFak
                end if
             else
                if (writeZ) then
                   write(iFile,1000) H%xMin+((/iBin1,iBin2/)*1.0-0.5)*H%xBin, Z, &
                        yVal(iBin1,iBin2,1:3)*MulFak+AddFak
                else
                   write(iFile,1000) H%xMin+((/iBin1,iBin2/)*1.0-0.5)*H%xBin, &
                        yVal(iBin1,iBin2,1:3)*MulFak+AddFak
                end if
             end if
          end do
          write(iFile,*)
       end do

       write(iFile,*)
       write(iFile,*)
    end do

    if (present(file)) then
       close(iFile)
    end if

    if (present(dump)) then
       if (dump) then
          if (present(file)) then
             call DumpHist2DMP(H,file//".bin",iFile)
          else
             call DumpHist2DMP(H,"DumpHist2DMP.bin",iFile)
          end if
       end if
    end if

1000 FORMAT(1X,' ',1P,2E12.4,4E12.4)

  end subroutine WriteHist2DMP_Gnuplot

  !****************************************************************************
  !****s* hist2DMPf90/DumpHist2DMP
  ! NAME
  ! subroutine DumpHist2DMP(H,file,iFile)
  ! PURPOSE
  ! Write all the histogram information unformatted (i.e. binary) to a file
  !
  ! INPUTS
  ! * type(histogram2DMP) :: H     -- Histogramm to be used
  ! * character*(*)       :: file  -- name of file to open and close
  ! * integer,OPTIONAL    :: iFile -- File number output to redirect [OPTIONAL]
  ! OUTPUT
  ! H is written UNFORMATTED to the given file
  !
  !****************************************************************************
  subroutine DumpHist2DMP(H,file,iFile)

    type(histogram2DMP),intent(in)          :: H
    character*(*),      intent(in)          :: file
    integer,            intent(in),optional :: iFile

    integer :: iF

    iF=121
    if (present(iFile)) iF = iFile

    open(iF,file=file,status='UNKNOWN',form='UNFORMATTED')
    rewind(iF)

    write(iF) H%iSet
    write(iF) H%xMin,H%xMax,H%xBin,H%xExtreme
    write(iF) H%Name
    write(iF) H%yVal

    close(iF)

  end subroutine DumpHist2DMP


  !****************************************************************************
  !****s* hist2DMPf90/FetchHist2DMP
  ! NAME
  ! subroutine FetchHist2DMP(H,file,iFile)
  ! PURPOSE
  ! Read in all the histogram information previously dumped unformatted
  ! (i.e. binary) to a file
  !
  ! INPUTS
  ! * character*(*)   :: file  -- name of file to open and close
  ! * integer,OPTIONAL:: iFile -- File number input to redirect [OPTIONAL]
  ! OUTPUT
  ! * type(histogram2DMP) :: H     -- Histogramm to be used
  !
  ! H is read UNFORMATTED from the given file. Sizes are calculated as in
  ! CreateHist, also memory is allocated.
  !
  ! NOTES
  ! No checks about input are performed!
  !****************************************************************************
!   subroutine FetchHist2DMP(H,file,iFile)
!     use histMPf90
!
!     type(histogram2DMP),intent(inout)       :: H
!     character*(*),  intent(in)          :: file
!     integer,        intent(in),optional :: iFile
!
!     integer :: iF
!     integer,dimension(2) :: L
!     integer :: nHist
!
!     if (initFlag) then
!        call histMPf90_Init
!        initFlag = .false.
!     end if
!
!     iF=121
!     if (present(iFile)) iF = iFile
!
!     open(iF,file=file,status='UNKNOWN',form='UNFORMATTED')
!     rewind(iF)
!
!     read(iF) H%iSet
!
!     nHist = Map2HistMP_getN(H%iSet)
!     if (nHist==0) then
!        write(*,*) 'FetchHist2DMP: H%iSet invalid:',H%iSet
!        stop
!     endif
!
!
!     if(ALLOCATED(H%xExtreme)) deallocate(H%xExtreme)
!     if(ALLOCATED(H%yVal)) deallocate(H%yVal)
!
!     allocate(H%xExtreme(nHist,2,2))
!
!     read(iF) H%xMin,H%xMax,H%xBin,H%xExtreme
!     read(iF) H%Name
!
!     L = nint( (H%xMax-H%xMin)/H%xBin )+1
!     allocate(H%yVal(nHist,-1:L(1),-1:L(2),3))
!
!     read(iF) H%yVal
!
!     close(iF)
!
!   end subroutine FetchHist2DMP



end module hist2DMPf90
