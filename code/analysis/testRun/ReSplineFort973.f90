program main

  use hist2df90

  implicit none

  real, dimension(2) :: xmin, xmax, xbin
  type(histogram2D) :: HH(0:7)
  integer :: i,j
  character*200 :: BUF
  character*200, dimension(0:7)::  Name

  read(*,*) xmin
  read(*,*) xmax
  read(*,*) xbin
  do i=0,7
     read(*,'(A)') BUF
     name(i) = trim(BUF)

     write(*,*) 'Name: >',trim(Name(i)),'<'
  enddo

  do i=0,7
     call CreateHist2D(HH(i), Name(i)(1:lName(i)), xmin,xmax,xbin, .true.)
     open(20+i,file=trim(Name(i)),status="unknown")
     call ReadHist2D_Gnuplot(HH(i), 20+i, DoAve=.true.)
     close(20+i)

     Name(i) = trim(Name(i))//"_B"
     open(20+i,file=trim(Name(i)),status="unknown")
     call WriteHist2D_Gnuplot_Bspline(HH(i), 20+i, DoAve=.true.,maxval=0.0)
     close(20+i)
  enddo
     

  


end program main
