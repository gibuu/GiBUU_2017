! This is a test case for histogram integration.
!
! It reads a present histogram from a file and produces integrated
! versions of it. The parameters can be adjusted via the namelist
! "integrateHist", cf. jobcard "jobIntegrate".

program integrate_Hist

  use output
  use histf90
  implicit none

  type(histogram) :: hist

  integer, parameter :: NumFiles = 20
  integer, parameter :: strLen = 100
  character(len=strLen), dimension(1:NumFiles) :: fileNames = ""  ! histogram files to be read
  integer :: ios,j
  logical :: backward, normalize

  namelist /integrateHist/ fileNames, backward, normalize

  ! read namelist input
  call Write_ReadingInput('integrateHist',0)
  rewind(5)
  read(5,nml=integrateHist,iostat=ios)
  call Write_ReadingInput('integrateHist',0,ios)
  print *,"backward     = ", backward
  print *,"normalize    = ", normalize
  call Write_ReadingInput('integrateHist',1)

  do j=1,NumFiles
  
    if (trim(fileNames(j)) == "") cycle
    print *,"file name: ", trim(fileNames(j))

    ! read histogram file
    if (ReadHist (hist, fileNames(j))) then

      ! reproduce original histogram
      call WriteHist (hist, file=trim(fileNames(j))//".plain")

      ! write out integrated histograms
      call WriteHist_Integrated (hist, trim(fileNames(j))//".integrated", backward, normalize)
      
    else
    
      print *,"error when reading file: ", fileNames(j)
    
    end if

  end do

end program
 
