! This is a test case for histogram rebinning.
!
! It reads a present histogram from a file and produces a rebinned
! copy of it.  The parameters can be adjusted via the namelist
! "rebinHist", cf. jobcard "jobRebin".

program rebin_HistMC

  use output
  use histMC
  implicit none

  type(histogramMC) :: hist

  character(len=100) :: filename  ! histogram file to be read
  integer :: nPoints = 5          ! rebinning factor
  integer :: ios

  namelist /rebinHist/ filename, nPoints

  ! read namelist input
  call Write_ReadingInput('rebinHist',0)
  rewind(5)
  read(5,nml=rebinHist,iostat=ios)
  call Write_ReadingInput('rebinHist',0,ios)

  print *,"histogram file: ", trim(filename)
  print *,"rebinning = ", nPoints
  call Write_ReadingInput('rebinHist',1)

  ! read histogram file
  if (.not. ReadHistMC (hist, filename)) then
    print *,"error reading file: ", filename
    stop
  end if

  ! reproduce original histogram
  call WriteHistMC (hist, file=trim(filename)//".orig")

  ! do the rebinning
  call rebinHistMC (hist, nPoints)

  ! write rebinned histogram
  call WriteHistMC (hist, file=trim(filename)//".rebinned")

end program 
