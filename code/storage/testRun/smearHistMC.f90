! This is a test case for histogram smearing.
!
! It reads one or more  histograms from file and produces resolution smeared
! copies of them. The smearing is done by convoluting with a Gaussian function.
! The parameters can be adjusted via the namelist "smearHist",
! cf. jobcard "jobSmear".

program smear_HistMC

  use output
  use histMC
  implicit none

  type(histogramMC) :: hist

  integer, parameter :: NumFiles = 20
  integer, parameter :: strLen = 100
  character(len=strLen), dimension(1:NumFiles) :: fileNames = ""  ! histogram files to be read
  real :: sigma                                                   ! width parameter of the Gauss function
  integer :: ios,j

  namelist /smearHist/ fileNames, sigma

  ! read namelist input
  call Write_ReadingInput('smearHist',0)
  rewind(5)
  read(5,nml=smearHist,iostat=ios)
  call Write_ReadingInput('smearHist',0,ios)
  print *,"sigma     = ", sigma
  call Write_ReadingInput('smearHist',1)

  
  do j=1,NumFiles

    if (trim(fileNames(j)) == "") cycle
    print *,"file name: ", trim(fileNames(j))
    
    ! read histogram file
    if (ReadHistMC (hist, fileNames(j))) then

      ! reproduce original histogram
      call WriteHistMC (hist, file=trim(fileNames(j))//".plain")

      ! write out smeared histograms
      call WriteHistMC_Gauss (hist, trim(fileNames(j))//".gauss", sigma)
      
    else
    
      print *,"error when reading file: ", fileNames(j)
    
    end if
    
  end do

end program
