! This is a test case for histogram smearing.
!
! It reads a present histogram from a file and produces resolution smeared
! copies of it. The smearing is done by convoluting with a Gaussian or
! Novosibirsk function. The parameters can be adjusted via the namelist
! "smearHist", cf. jobcard "jobSmear".

program smear_Hist

  use output
  use histf90
  implicit none

  type(histogram) :: hist

  character(len=100) :: filename  ! histogram file to be read
  real :: sigma                   ! width parameter of the Gauss and Novosibirsk function
  real :: tau                     ! skewness parameter of the Novosibirsk function
  integer :: nPoints = 5          ! number of points for splining
  integer :: ios

  namelist /smearHist/ filename, sigma, tau, nPoints

  ! read namelist input
  call Write_ReadingInput('smearHist',0)
  rewind(5)
  read(5,nml=smearHist,iostat=ios)
  call Write_ReadingInput('smearHist',0,ios)
  print *,"histogram file: ", trim(filename)
  print *,"sigma     = ", sigma
  print *,"tau       = ", tau
  print *,"numPoints = ", nPoints
  call Write_ReadingInput('smearHist',1)

  ! read histogram file
  if (ReadHist (hist, filename)) then

    ! reproduce original histogram
    call WriteHist (hist, file=trim(filename)//".plain")

    ! write out smeared histograms
    call WriteHist_Gauss (hist, trim(filename)//".gauss", sigma)
    call WriteHist_Novo  (hist, trim(filename)//".novo",  sigma, tau)

    ! write out splined histograms
    call WriteHist_Spline  (hist, trim(filename)//".spline",  nPoints=nPoints)
    call WriteHist_BSpline (hist, trim(filename)//".bspline", nPoints=nPoints)
    
  else
  
    print *,"error when reading file: ", filename
  
  end if

end program
