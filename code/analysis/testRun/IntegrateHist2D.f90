!***************************************************
!* read in the 2D-Histogram form File 'H2D.dat.bin'
!* and writes the 2 1-dimensional histograms
!* 'H2D.x.dat', and 'H2D.y.dat', which represents
!* the integrals over the y or x component.
!****************************************************

program IntegrateHist2Dxy

  use histf90
  use hist2Df90

  IMPLICIT NONE

  type(histogram)   :: HHH
  type(histogram2D) :: H2D
  logical :: flagOK
  real :: add, mul

  call FetchHist2D(H2D,"H2D.dat.bin", add=add,mul=mul,flagOk=flagOK)

  if (.not.flagOK) then
     write(*,*) 'Error while reading "H2D.dat.bin". stop!'
     stop
  end if

  call IntegrateHist2D(h2D,HHH,1)
  call WriteHist(HHH,file="H2D.x.dat",add=add,mul=mul,dump=.true.)

  call IntegrateHist2D(h2D,HHH,2)
  call WriteHist(HHH,file="H2D.y.dat",add=add,mul=mul,dump=.true.)

  

  

end program IntegrateHist2Dxy
