program TryTransform
  use random
  use histf90

  integer :: iMC,nMC
  type(histogram) :: h1, h2

  real :: nu,M,Q2
  real :: W, gam,betgam
  real :: x1,x2, w2

  nMC = 100000

  call CreateHist(h1,"gen1",-1.0,1.0,0.1)
  call CreateHist(h2,"gen2",-1.0,1.0,0.1)
  
  M =.938
  Q2 = 1.5
  W = 1.7

  nu = (W**2-M**2+Q2)/(2*M)

  gam = (nu+M)/W
  betgam = sqrt(nu**2+Q2)/W

  do iMC = 1,nMC

     x1=2*rn()-1 !rnCos()

     call AddHist(h1,x1,1.0)

     x2 = (betgam+gam*x1)/(gam+betgam*x1)
     w2 = 1/(gam+betgam*x1)**2

     call AddHist(h2,x2,w2)

  end do

  call WriteHist(h1,131,add=0.0,mul=1.0/nMC,file="H1.txt")
  call WriteHist(h2,131,add=0.0,mul=1.0/nMC,file="H2.txt")

end program TryTransform
