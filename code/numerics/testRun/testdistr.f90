program main
  
  use random
  use histf90
  implicit none

  integer :: iMC,nMC=1000000, ix0, nx0=50
  real :: x, x2, x3, x0, x0max=0.25
  real, dimension(3) :: vec
  real :: S, S2, S3

  type(histogram) :: H1,H2,H3


  call CreateHist(H1, "x (=pT2)",0.0,4.0,0.02)
  call CreateHist(H2, "x*fermi",0.0,4.0,0.02)
  call CreateHist(H3, "fermi",0.0,4.0,0.02)

  S= 0.0
  S2 = 0.0
  S3 = 0.0

  do ix0=1,nx0
     S= 0.0
     S2 = 0.0
     S3 = 0.0

     do iMC=1,nMC
        x = rnExp(-3.0)
        S = S + x
        vec = rnOmega()
        x0 = ix0*0.010 *rn()**(1./3.)
        x2 = (sqrt(x)+x0*vec(1))**2 + (x0*vec(2))**2
        S2 = S2+x2
        x3 = (x0*vec(1))**2 + (x0*vec(2))**2
        S3 = S3+x3
     end do

     write(*,'(1P,5e14.5)') ix0*0.010,S/nMC, S2/nMC,S3/nMC
  end do



  do iMC=1,nMC
     x = rnExp(-3.0)

     call AddHist(H1,x,1.0)
     S = S+x

     x0 = x0max

     vec = rnOmega()

     x0 = x0max
     x0 = x0max *rn()**(1./3.)

!!$     do
!!$        vec(1) = 2*rn()-1
!!$        vec(2) = 2*rn()-1
!!$        vec(3) = 2*rn()-1
!!$        
!!$        if (vec(1)**2+vec(2)**2+vec(3)**2 .lt. 1.0) exit
!!$
!!$     end do
!!$     x0=x0max


     x2 = (sqrt(x)+x0*vec(1))**2 + (x0*vec(2))**2
     call AddHist(H2,x2,1.0)
     S2 = S2+x2


     x3 = (x0*vec(1))**2 + (x0*vec(2))**2
     call AddHist(H3,x3,1.0)

  enddo

  call WriteHist(H1,101,1e-20,1.0/nMC)
  call WriteHist(H2,102,1e-20,1.0/nMC)
  call WriteHist(H3,103,1e-20,1.0/nMC)

  write(*,*) S/nMC, S2/nMC

end program main
