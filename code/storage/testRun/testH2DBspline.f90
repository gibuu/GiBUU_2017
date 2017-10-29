program main
  use hist2df90
  use random

  implicit none

  type(histogram2D),save :: HH,HH1
  integer :: i,ix,iy
  real :: x ,y, z

  call CreateHist2D(HH, 'HH', (/0.0,0.0/), (/10.0,20.0/) , (/1.0,1.0/), .true.)
  call CreateHist2D(HH1, 'HH1', (/0.0,0.0/), (/10.0,20.0/) , (/1.0,1.0/), .true.)
  
  do i=1,1000
     x = rn() * 10.
     y = rn() * 20.
     z = rn()
     
!!$     z = 0
!!$     if (x>3.) z=z+1
!!$     if (y>5.) z=z+1

     call AddHist2D(HH, (/x,y/) , z, 2.*z-1.)
  enddo

!!$  do ix = 0,20
!!$     x = ix*0.5
!!$     do iy=0,20
!!$        y = iy*0.5
!!$        z = 0
!!$        if (x>3.) z=z+1
!!$        if (y>5.) z=z+1
!!$        
!!$        call AddHist2D(HH, (/x,y/) , z, 2.*z-1.)
!!$     enddo
!!$  enddo

  call WriteHist2D_Gnuplot(HH, 21,DoAve=.true.,maxval=0.0)
  call WriteHist2D_Gnuplot_Bspline(HH, 22,DoAve=.true.,maxval=0.0)

  rewind(21)

  call ReadHist2D_Gnuplot(HH1, 21, DoAve=.true.)
  call WriteHist2D_Gnuplot(HH1, 23)


end program main
